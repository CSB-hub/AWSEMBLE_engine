"""
프라이머 설계 및 검증 (Basic 모듈)

이 모듈은 **기본적인 PCR 프라이머 설계 및 검증**만 제공합니다.

핵심 원칙:
    - Basic 모듈: 순수 프라이머 증폭 기능만
    - Assembly 모듈: 각 방법별 특수 로직 (overlap, 제한효소 등)
    - 코드 재사용: Assembly 클래스에서 이 모듈 호출

핵심 기능:
    - 기본 프라이머 설계 (특정 영역 증폭)
    - 프라이머 dimer 체크 (primer3-py)
    - GC clamp 검증
    - 특이성 체크 (off-target binding)
    - 종합적인 프라이머 품질 검증

라이브러리:
    - primer3-py: 프라이머 설계 및 dimer 체크 (업계 표준)
    - BioPython: Tm, GC% 계산

작성일: 2025-10-30
수정일: 2025-10-30 (Basic/Assembly 분리, primer3-py 도입)
"""

from typing import Dict, Tuple, Optional, List, Any

# 기존 모듈 import
from .sequence_analysis import calculate_tm, calculate_gc_content, get_typeIIS_params
from .secondary_structure import (
    analyze_primer_structures,
    analyze_target_accessibility,
    STRUCT_TEMP_C,
    DIMER_DG_CUTOFF_KCAL,
)

import primer3

#: primer3에서 받아올 후보 쌍 수. 5' 시작이 고정이면 길이만 변하므로
#  길이창 폭(보통 8)보다 넉넉하면 충분하다. 우리 축으로 고르려면 후보가 여럿이어야 한다.
PRIMER3_NUM_RETURN = 30

#: 길이 힌트 탐색(`_length_hint`) 전용 후보 수. 여기서 필요한 답은 "한 쌍이라도
#  창에 드는가" 뿐이라 후보를 많이 받을 이유가 없다. 설계 본체(위 30)와 달리
#  **고를 필요가 없으므로** 줄여서 실패 경로의 비용을 깎는다.
_HINT_NUM_RETURN = 8

#: 길이 힌트 탐색의 굵은 격자 간격(nt). 1nt씩 전수로 훑으면 길이창 끝까지
#  최대 ~20회 primer3를 돌려 체감이 멈춘다. 굵게 훑어 구간만 잡고 그 구간만
#  되짚으면 같은 답을 2~6회로 얻는다.
_HINT_COARSE_STEP = 4


# ─────────────────────────────────────────────────────────────────────────────
# 설계 범위를 못 맞췄을 때 **왜 그런지 / 무엇을 바꾸면 되는지** 를 만든다.
#
# 막연히 "범위를 넓혀 보세요"가 아니라 primer3가 직접 말해주는 사실과 실제 탐색으로
# 구한 숫자만 쓴다. 2026-11-05까지는 요청 범위를 못 맞추면 그냥 기본값으로 갈아타고
# "기본 범위로 설계했습니다"라고만 알려서, 유저는 **무엇을 바꿔야 할지 알 수 없었다.**
# ─────────────────────────────────────────────────────────────────────────────

def _explain_shortfall_axis(tm_fwd: float, tm_rev: float,
                            lo: float, hi: float) -> str:
    """어느 쪽이 왜 창을 못 맞췄는지 — **우리가 잰 값으로** 만든다.

    2026-11-05 오전까지는 primer3의 `PRIMER_LEFT_EXPLAIN`을 읽어 사유를 만들었다.
    지금은 primer3에 Tm 창을 걸지 않으므로(그쪽은 후보 생성만 한다) EXPLAIN에 Tm
    실패가 찍히지 않는다. 대신 **선택된 쌍의 실제 Tm**으로 판단한다 — 추론이 아니라
    보고하는 숫자 그대로라 시트와 어긋날 수가 없다.
    """
    f_out = tm_fwd < lo or tm_fwd > hi
    r_out = tm_rev < lo or tm_rev > hi
    if f_out and r_out:
        side = "양쪽 모두"
    elif r_out:
        side = "reverse 쪽이"
    elif f_out:
        side = "forward 쪽이"
    else:
        return ""
    outs = [t for t, o in ((tm_fwd, f_out), (tm_rev, r_out)) if o]
    if all(t < lo for t in outs):
        why = "Tm이 요청 범위보다 낮습니다(결합부가 짧아서)"
    elif all(t > hi for t in outs):
        why = "Tm이 요청 범위보다 높습니다(결합부가 길어서)"
    else:
        why = "한쪽은 낮고 한쪽은 높아 한 창에 같이 넣을 수 없습니다"
    return f"{side} 범위를 못 맞췄습니다 — {why}."


def _length_hint(settings: Dict, prm, lo: float, hi: float) -> str:
    """요청 Tm을 **맞출 수 있는 결합부 길이**를 실제 탐색으로 찾아 알려준다.

    길이창(`primer_len_min`~`primer_len_max`, GUI의 "Annealing length")이 Tm 도달을
    막는 숨은 제약인 경우가 많다. 상한을 올려 되는지 먼저 보고, 안 되면 하한을 내려 본다
    (낮은 Tm을 원할 땐 **줄여야** 한다 — 방향이 반대다).

    ⚠️ 판정은 **우리 축**(`calculate_tm`)으로 한다. primer3에 Tm 창을 걸어 판정하면
    축이 갈려서, 제안대로 바꿔도 경고가 안 사라지는 일이 생긴다.
    primer3 호출이 몇 번 더 들어가지만 **창을 못 맞춘 실패 경로에서만** 돈다.
    """
    from params import PRIMER3_MAX_SIZE

    def feasible(lmin, lmax):
        st = dict(settings)
        st.update({"PRIMER_MIN_SIZE": lmin, "PRIMER_MAX_SIZE": lmax,
                   "PRIMER_OPT_SIZE": min(max(prm.primer_len_opt, lmin), lmax),
                   "PRIMER_NUM_RETURN": _HINT_NUM_RETURN})
        try:
            r = primer3.bindings.design_primers(st, {})
        except Exception:
            return False
        for i in range(r.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            tf = calculate_tm(r[f"PRIMER_LEFT_{i}_SEQUENCE"], **prm.tm_kwargs())
            tr = calculate_tm(r[f"PRIMER_RIGHT_{i}_SEQUENCE"], **prm.tm_kwargs())
            if lo <= tf <= hi and lo <= tr <= hi:
                return True
        return False

    def first_feasible(cands, make_window):
        """`cands`를 앞에서부터 보며 **처음 되는 값**을 찾는다 — 굵게 훑고 되짚는 방식.

        1nt씩 전수로 보면 같은 답에 최대 ~20회 primer3가 들어간다. 길이창을 넓히면
        후보 집합은 커지기만 하므로 "되던 게 더 넓혀서 안 될" 일은 없다(단조). 그래서
        `_HINT_COARSE_STEP` 간격으로 먼저 되는 지점을 잡고, **직전 실패 지점과 그
        사이만** 1칸씩 되짚으면 전수와 같은 최소값을 2~6회로 얻는다.
        """
        hit = None
        last_no = -1
        # 격자를 **0번부터** 놓는다 — 실제로 필요한 변경은 보통 1~2nt라
        # 첫 칸이 곧 답인 경우가 가장 흔하고, 그 경우 호출 1회로 끝난다.
        for i in range(0, len(cands), _HINT_COARSE_STEP):
            if feasible(*make_window(cands[i])):
                hit = i
                break
            last_no = i
        if hit is None:
            # 굵은 격자가 끝까지 못 짚었을 수 있다(격자에 안 걸린 꼬리 구간) → 끝만 확인
            if not cands or not feasible(*make_window(cands[-1])):
                return None
            hit = len(cands) - 1
        for j in range(last_no + 1, hit):
            if feasible(*make_window(cands[j])):
                return cands[j]
        return cands[hit]

    lo0, hi0 = prm.primer_len_min, prm.primer_len_max
    _FLOOR = 10          # 이보다 짧으면 특이성이 깨져 PCR 자체가 성립하지 않는다

    lmax = first_feasible(list(range(hi0 + 1, PRIMER3_MAX_SIZE + 1)),
                          lambda v: (lo0, v))
    if lmax is not None:
        return (f"결합부 길이(Annealing length) **상한을 {hi0} → {lmax}** 로 올리면 "
                f"요청 범위를 맞출 수 있습니다.")
    lmin = first_feasible(list(range(lo0 - 1, _FLOOR - 1, -1)),
                          lambda v: (v, hi0))
    if lmin is not None:
        return (f"결합부 길이(Annealing length) **하한을 {lo0} → {lmin}** 로 내리면 "
                f"요청 범위를 맞출 수 있습니다(낮은 Tm은 짧게 가야 합니다).")
    return (f"결합부 길이를 {_FLOOR}~{PRIMER3_MAX_SIZE}nt 어디로 바꿔도 이 서열에서는 "
            f"요청 Tm에 닿지 않습니다 — Tm 범위를 조정해 주세요.")


def check_pair_tm_diff(fwd_tm, rev_tm, params=None):
    """fwd·rev **두 프라이머 사이**의 Tm 차이를 선언창(`pair_tm_diff_max`)과 비교한다.

    `tm_diff_max`와 혼동하지 않도록 — 그쪽은 **같은 프라이머**의 anneal Tm과 전체
    프라이머(5' tail 포함) Tm의 차이이고, GG 프라이머 길이를 정하는 knob이다.
    이쪽은 쌍 사이의 차이이고, **annealing 온도 하나로 두 프라이머를 돌릴 수 있는가**를
    본다: `ta_recommended()`가 Ta를 **낮은 쪽 Tm 기준**으로 정하므로 차이가 벌어지면
    뜨거운 쪽이 자기 Tm보다 훨씬 낮은 온도에서 돌아 비특이 결합이 늘어난다.

    설계기 두 곳(`design_basic_primers` · `design_goldengate_primers`)이 같은 자를
    쓰도록 여기 한 벌만 둔다.

    Args:
        fwd_tm: forward 프라이머 Tm (℃)
        rev_tm: reverse 프라이머 Tm (℃)
        params: DesignParams (None이면 DEFAULT_PARAMS)

    Returns:
        선언창을 벗어나면 경고 문자열, 아니면 None
    """
    from params import DEFAULT_PARAMS
    _p = params or DEFAULT_PARAMS
    if not isinstance(fwd_tm, (int, float)) or not isinstance(rev_tm, (int, float)):
        return None
    diff = abs(fwd_tm - rev_tm)
    if diff <= _p.pair_tm_diff_max:
        return None
    _ta = _p.ta_recommended(fwd_tm, rev_tm)
    _hot = "forward" if fwd_tm > rev_tm else "reverse"
    _msg = (f"fwd·rev Tm 차이가 설정한 상한보다 큽니다 "
            f"({diff:.1f}℃ > {_p.pair_tm_diff_max:g}℃ — "
            f"fwd {fwd_tm:.1f}℃ / rev {rev_tm:.1f}℃)")
    if _ta is not None:
        _msg += (f" — annealing 온도는 낮은 쪽 기준 {_ta:.1f}℃로 정해지므로 "
                 f"{_hot} 프라이머가 자기 Tm보다 낮은 온도에서 돌아 "
                 f"비특이 결합이 늘 수 있습니다")
    return _msg


class _DimerScan(list):
    """`check_cross_dimers`의 결과 — 걸린 쌍 리스트에 **검사한 쌍 수**를 달아 보낸다.

    빈 리스트만 돌려주면 "검사했는데 깨끗하다"와 "검사를 안 했다"가 구분이 안 된다.
    분모(`pairs_checked`)가 있어야 리포트가 "12쌍 중 0쌍"이라고 쓸 수 있다.
    (리스트를 상속하므로 기존 호출부의 `len()`·반복·인덱싱은 그대로 동작한다)
    """

    pairs_checked = 0


def check_cross_dimers(anneals, params=None):
    """한 반응에 같이 들어가는 프라이머들의 **결합부끼리** dimer를 전수 검사한다.

    ## 왜 결합부만 보는가

    one-pot Golden Gate는 부품 n개의 프라이머 2n개가 **같은 튜브**에 들어간다.
    4부품이면 8개 = 28쌍인데, 지금까지는 자기 쌍(fwd↔rev) 4쌍만 봤다.
    그러나 **전체 프라이머로 전수 검사하면 소음이 된다** — GG 프라이머는 전부 같은
    5' tail(`[패딩][인식서열][spacer]`)을 공유하고, 그 tail만으로 이미 기준을 넘는다
    (이 프로젝트 실측: tail 12bp 단독 self ΔG −7.44 / 역상보와 −10.82, 기준 −5.0).
    2026-10-20에 ΔG 경고를 화면에서 뺀 이유가 바로 이것이었다(51개 중 50개 = 98%에
    경고가 붙어 정보가 되지 못했다).

    그래서 **tail을 뺀 결합부끼리만** 비교한다. 결합부는 overhang·spacer 선정이나
    결합부 길이 재조정으로 **실제로 바꿀 수 있는** 값이라, 걸리면 손쓸 수 있다.

    ## 등급을 둘로 나눈다 — 실측으로 재조정했다

    처음엔 `qc_dimer_dg`(−5.0) 하나로 전부 차단했는데, 실측하니 **39행 중 24행(62%)이
    red**가 됐다(쌍 35건, ΔG 중앙 −5.85). 2026-10-20에 ΔG 경고를 화면에서 뺀 상황이
    그대로 재현된 것이다 — 기본값에서 과반이 뜨는 신호는 정의상 정보가 아니다.

    그래서 이 프로젝트가 **이미 가진 ΔG 어휘**를 그대로 쓴다(새 숫자를 만들지 않는다):

    | ΔG | 등급 | 처리 |
    |---|---|---|
    | ≤ `qc_dg_absolute` (−9.0) | 강한 dimer | **차단**(violation) — 실측 35건 중 3건 |
    | ≤ `qc_dimer_dg` (−5.0) | 중간 dimer | 기록만 — `cross_checks.csv`에 남기고 행은 안 건드린다 |

    차단을 `qc_dg_absolute`에 두는 근거: 그 필드가 원래 "위치·형태 불문 최상위 경고
    tier"로 정의돼 있고, `check_primer_dimer` 독스트링의 등급표에서도 −9가 "강한
    dimer" 경계다. tail 유래 dimer는 고칠 수 없어 경고로 두는데(`warn_only_dimer`),
    결합부끼리 붙는 것은 고칠 수 있고 증폭 효율을 직접 깎으므로(프라이머가 서로 붙어
    있으면 주형에 못 간다) 강한 등급은 차단으로 올린다.

    Args:
        anneals: `[(라벨, 결합부 서열), ...]` — 한 반응에 같이 들어가는 것만.
            라벨은 리포트에 그대로 나가므로 `"BBa_C0012_3_designed fwd"` 처럼
            사람이 읽을 수 있어야 한다.
        params: DesignParams (None이면 DEFAULT_PARAMS).
            `qc_dimer_dg`(기록 경계)와 `qc_dg_absolute`(차단 경계)를 쓴다.

    Returns:
        `_DimerScan` — 걸린 쌍 리스트 `[{'a','b','delta_g','blocking','detail'}, ...]`에
        `pairs_checked`(실제로 비교한 쌍 수)가 달려 있다. `blocking`이 True면 차단
        등급이다. ΔG가 더 음수인 순서로 정렬한다.

    Note:
        자기 쌍(fwd↔rev)은 **제외한다** — `validate_primer(check_dimer_with=...)`가
        이미 전체 프라이머 기준으로 본다. 여기서 또 보면 같은 쌍이 두 번 보고된다.
    """
    from params import DEFAULT_PARAMS
    from itertools import combinations
    _p = params or DEFAULT_PARAMS
    hits = []
    _n_pairs = 0        # 실제로 비교한 쌍 수 — "0건"의 분모가 된다
    items = [(lab, (seq or "").upper()) for lab, seq in anneals if (seq or "").strip()]
    for (la, sa), (lb, sb) in combinations(items, 2):
        # 같은 부품의 fwd↔rev는 이미 검사됐다 — 라벨의 부품 부분이 같으면 건너뛴다.
        if la.rsplit(" ", 1)[0] == lb.rsplit(" ", 1)[0]:
            continue
        _n_pairs += 1
        r = ValidatePrimer.check_primer_dimer(sa, sb, threshold_dg=_p.qc_dimer_dg)
        if r.get("has_dimer") and r.get("warning"):
            _dg = r.get("delta_g")
            hits.append({"a": la, "b": lb,
                         "delta_g": _dg,
                         "blocking": isinstance(_dg, (int, float))
                                     and _dg <= _p.qc_dg_absolute,
                         "detail": r["warning"]})
    hits.sort(key=lambda h: (h["delta_g"] if isinstance(h["delta_g"], (int, float))
                             else 0.0))
    scan = _DimerScan(hits)
    scan.pairs_checked = _n_pairs
    return scan


class PrimerDesign:
    @staticmethod
    def design_basic_primers(
        template: str,
        target_start: int = 0,
        target_end: Optional[int] = None,
        target_tm: float = None,
        primer_size_range: Tuple[int, int] = None,
        check_specificity: bool = True,
        anchor_ends: bool = False,
        params=None,
    ) -> Dict:
        """
        기본 PCR 프라이머 설계 (증폭만)

        ⭐ Basic 모듈의 핵심 함수!
        ⭐ Assembly 클래스에서 이 함수를 호출하여 사용!

        primer3-py를 사용하여 특정 영역을 증폭하는 프라이머를 설계합니다.
        - Tm, GC%, dimer, hairpin 모두 자동 검증
        - 특이성 체크 (off-target binding)

        Args:
            template: Template DNA 서열
            target_start: 증폭 시작 위치 (0-based)
            target_end: 증폭 종료 위치 (None이면 전체)
            target_tm: 목표 Tm (°C). 미지정 시 `params`의 Tm 범위 중앙값(`tm_opt`).
                명시하면 **폭은 유지한 채 범위를 그 값 중심으로 옮긴다**
                (조건이 바뀌면 Tm 분포가 폭은 그대로 이동한다는 실측에 따른 규약).
            primer_size_range: 프라이머 길이 범위 (min, max)
            check_specificity: 특이성 체크 여부 (기본값: True)
            anchor_ends: 프라이머 5' 시작을 target 양끝에 고정할지 (기본값: False).
                True면 amplicon = target 구간 정확히(경계 고정), 길이/Tm만 조건 기반 가변.
                False면 SEQUENCE_TARGET 방식(target을 포함하는 최적 amplicon — flanking 가능).
                Gibson general-prep처럼 insert fragment 말단을 target 경계에 맞춰야 할 때 True.

        Returns:
            {
                'forward': str,                 # Forward 프라이머 서열
                'reverse': str,                 # Reverse 프라이머 서열
                'fwd_tm': float,                # Forward Tm
                'rev_tm': float,                # Reverse Tm
                'fwd_gc': float,                # Forward GC%
                'rev_gc': float,                # Reverse GC%
                'fwd_specificity': Dict,        # Forward 특이성 결과
                'rev_specificity': Dict,        # Reverse 특이성 결과
                'amplicon_size': int,           # 증폭 산물 크기
                'warnings': List[str]           # 경고 메시지들
            }

        Example:
            >>> template = "ATCGATCGATCG" * 10
            >>> result = design_basic_primers(template, 0, 100, target_tm=60)
            >>> result['forward'] is not None
            True
            >>> result['amplicon_size'] > 0
            True

        Notes:
            - Assembly 모듈에서 이 함수 호출 후 overlap/제한효소 추가
            - 순수한 증폭 프라이머만 설계 (추가 로직 없음)
        """
    
        if not template:
            raise ValueError("Template 서열이 비어있습니다")

        if target_end is None:
            target_end = len(template)

        warnings = []

        # 설계 파라미터 (연구실별 값). 미지정 시 기존 하드코딩과 동일한 기본값.
        from params import DEFAULT_PARAMS
        prm = params if params is not None else DEFAULT_PARAMS
        # Tm 범위 — 설계 목표이자 판정 기준(단일 창). 자세한 배경은 params.py의 tm_min/tm_max.
        _tm_lo, _tm_hi = prm.tm_min, prm.tm_max
        if target_tm is None:
            target_tm = prm.tm_opt
        else:
            # 호출자가 목표를 직접 준 경우: 폭을 유지한 채 창을 그 값 중심으로 이동.
            _half = (prm.tm_max - prm.tm_min) / 2.0
            _tm_lo, _tm_hi = target_tm - _half, target_tm + _half
        if primer_size_range is None:
            primer_size_range = (prm.primer_len_min, prm.primer_len_max)

        # primer3 설정 — 전부 params에서 온다.
        #
        # ⚠️ **primer3에는 Tm 창을 걸지 않는다**(2026-11-05). primer3는 자기 열역학
        #    모델로 Tm을 재는데 우리 축은 NEB 재현(`calculate_tm`)이라 서로 다르다.
        #    여기서 `PRIMER_MIN_TM/MAX_TM`을 주면 **primer3 축으로 거르고 우리 축으로
        #    보고하게 되어** 유저가 선언한 창과 실제 선정 기준이 갈린다 — 바로 그 "자가
        #    두 개" 버그를 오늘 고쳤으므로 되살리지 않는다.
        #    → primer3는 **후보 생성**만 하고, 창 판정과 선택은 아래에서 우리 축으로 한다.
        #
        # 염 농도는 primer3의 내부 랭킹·2차구조 계산에만 쓰이므로 **1가 환산값**을 준다
        # (Na + K + Tris/2 — Biopython MeltingTemp의 합산 규칙과 같다).
        primer3_settings = {
            'SEQUENCE_ID': 'template',
            'SEQUENCE_TEMPLATE': template,
            'SEQUENCE_TARGET': [target_start, target_end - target_start],
            'PRIMER_OPT_SIZE': prm.primer_len_opt,
            'PRIMER_MIN_SIZE': primer_size_range[0],
            'PRIMER_MAX_SIZE': primer_size_range[1],
            'PRIMER_OPT_TM': target_tm,     # 랭킹 힌트로만 남긴다(판정은 우리 축)
            'PRIMER_MIN_TM': 0.0,           # 창은 열어 둔다 — 아래에서 우리가 고른다
            'PRIMER_MAX_TM': 100.0,
            'PRIMER_NUM_RETURN': PRIMER3_NUM_RETURN,
            'PRIMER_MIN_GC': prm.gc_min,
            'PRIMER_MAX_GC': prm.gc_max,
            'PRIMER_MAX_POLY_X': prm.max_poly_x,
            'PRIMER_SALT_MONOVALENT': prm.na_mM + prm.k_mM + prm.tris_mM / 2.0,
            'PRIMER_DNA_CONC': prm.primer_conc_nM,
            'PRIMER_PRODUCT_SIZE_RANGE': [[50, 1000]],  # Amplicon 크기 범위 확대
            'PRIMER_MAX_SELF_ANY': 12.0,  # Self-complementarity 완화
            'PRIMER_MAX_SELF_END': 8.0,   # 3' self-complementarity 완화
            'PRIMER_PAIR_MAX_COMPL_ANY': 12.0,  # Pair complementarity 완화
            'PRIMER_PAIR_MAX_COMPL_END': 8.0,
        }

        # anchor_ends: 프라이머 5' 시작을 target 양끝에 고정 → amplicon=target 정확.
        #   fwd 5'=target_start, rev 5'=target_end-1 (primer3 좌표는 0-based, target_end exclusive).
        #   길이/Tm은 안쪽으로만 가변 → insert 말단 외 서열 오염 없음(Gibson general-prep용).
        #   SEQUENCE_TARGET(포함만 요구 → flanking 허용)과 충돌하므로 제거.
        if anchor_ends:
            primer3_settings.pop('SEQUENCE_TARGET', None)
            primer3_settings['SEQUENCE_FORCE_LEFT_START'] = target_start
            primer3_settings['SEQUENCE_FORCE_RIGHT_START'] = target_end - 1
            # 양끝 고정 시 amplicon = target 전체이므로 product size range를 그 길이에 맞춘다.
            # (기본 [50,1000]은 module(>1000bp) 전체 증폭에서 primer3가 후보를 못 만들어 실패 — P7 구멍1)
            tlen = target_end - target_start
            primer3_settings['PRIMER_PRODUCT_SIZE_RANGE'] = [[tlen, tlen]]

        # primer3 실행
        try:
            result = primer3.bindings.design_primers(
                primer3_settings,
                {}  # 빈 global args (기본 파라미터 사용)
            )

            # 결과 확인
            if result['PRIMER_PAIR_NUM_RETURNED'] == 0:
                error_msg = "⚠️  적절한 프라이머를 찾을 수 없습니다."
                if 'PRIMER_ERROR' in result:
                    error_msg += f" (Primer3 오류: {result['PRIMER_ERROR']})"
                if 'PRIMER_LEFT_EXPLAIN' in result:
                    error_msg += f"\n  Forward 설명: {result['PRIMER_LEFT_EXPLAIN']}"
                if 'PRIMER_RIGHT_EXPLAIN' in result:
                    error_msg += f"\n  Reverse 설명: {result['PRIMER_RIGHT_EXPLAIN']}"
                warnings.append(error_msg)
                return {
                    'forward': None,
                    'reverse': None,
                    'warnings': warnings
                }

            # ── 프라이머 쌍 선택 — **우리 축(NEB 재현)으로 고른다** ──────────────
            # primer3는 후보만 만들었다(위에서 Tm 창을 열어 뒀다). 유저가 선언한 창
            # [tm_min, tm_max]은 여기서 적용한다. 그래야 "선언한 창 = 선정 기준 = 시트
            # 숫자"가 한 축으로 유지된다.
            _n = result['PRIMER_PAIR_NUM_RETURNED']
            _cands = []
            for _i in range(_n):
                _f = result[f'PRIMER_LEFT_{_i}_SEQUENCE']
                _r = result[f'PRIMER_RIGHT_{_i}_SEQUENCE']
                _tf = calculate_tm(_f, **prm.tm_kwargs())
                _tr = calculate_tm(_r, **prm.tm_kwargs())
                # 창 밖으로 벗어난 거리(양쪽 프라이머 합). 창 안이면 0.
                _dist = (max(0.0, _tm_lo - min(_tf, _tr))
                         + max(0.0, max(_tf, _tr) - _tm_hi))
                _cands.append((_i, _f, _r, _tf, _tr, _dist))

            _in_win = [c for c in _cands if c[5] <= 0.0]
            if _in_win:
                # 창 안에서는 **두 프라이머 Tm이 가장 비슷한 쌍**을 고른다. 한쪽만 뜨거우면
                # 낮은 쪽이 안 붙거나 높은 쪽이 비특이 결합한다(Ta는 낮은 쪽이 정한다).
                _pick = min(_in_win, key=lambda c: (abs(c[3] - c[4]), len(c[1]) + len(c[2])))
            else:
                # 창 안에 하나도 없으면 **가장 가까운 쌍**을 쓰고 알린다(조립은 진행).
                # 기본값으로 갈아타지 않는다 — 그러면 무엇을 요청하든 같은 프라이머가
                # 나와서 "값을 바꿔도 결과가 안 변한다"가 된다(2026-11-05 Y9).
                _pick = min(_cands, key=lambda c: (c[5], abs(c[3] - c[4])))
                _msg = (f"요청하신 Tm {_tm_lo:g}~{_tm_hi:g}℃ 안에서 프라이머 **쌍**을 "
                        f"만들지 못해 가장 가까운 쌍으로 설계했습니다 "
                        f"(fwd {_pick[3]:.1f}℃/{len(_pick[1])}nt · "
                        f"rev {_pick[4]:.1f}℃/{len(_pick[2])}nt).")
                _why = _explain_shortfall_axis(_pick[3], _pick[4], _tm_lo, _tm_hi)
                if _why:
                    _msg += f" {_why}"
                _hint = _length_hint(primer3_settings, prm, _tm_lo, _tm_hi)
                if _hint:
                    _msg += f" {_hint}"
                warnings.append(_msg)
                print(f"     ⚠️ {_msg}")

            _k = _pick[0]
            fwd_seq = result[f'PRIMER_LEFT_{_k}_SEQUENCE']
            rev_seq = result[f'PRIMER_RIGHT_{_k}_SEQUENCE']
            fwd_start = result[f'PRIMER_LEFT_{_k}'][0]
            rev_start = result[f'PRIMER_RIGHT_{_k}'][0]

            # 프라이머 특성 계산
            fwd_tm = calculate_tm(fwd_seq, **prm.tm_kwargs())
            rev_tm = calculate_tm(rev_seq, **prm.tm_kwargs())
            fwd_gc = calculate_gc_content(fwd_seq)
            rev_gc = calculate_gc_content(rev_seq)

            # 쌍 ΔTm 판정 — Ta는 낮은 쪽 Tm 기준이라 차이가 벌어지면 뜨거운 쪽이
            # 자기 Tm보다 낮은 온도에서 돈다(§check_pair_tm_diff).
            _pair_w = check_pair_tm_diff(fwd_tm, rev_tm, prm)
            if _pair_w:
                warnings.append(_pair_w)

            # Amplicon 크기 계산
            amplicon_size = rev_start - fwd_start + len(fwd_seq)

            # 특이성 체크
            fwd_specificity = None
            rev_specificity = None

            if check_specificity:
                fwd_specificity = ValidatePrimer.check_primer_specificity(
                    template, fwd_seq, target_start, target_end, is_reverse=False
                )
                rev_specificity = ValidatePrimer.check_primer_specificity(
                    template, rev_seq, target_start, target_end, is_reverse=True
                )

                # 특이성 경고 추가
                if fwd_specificity['warning']:
                    warnings.append(f"Forward: {fwd_specificity['warning']}")
                if rev_specificity['warning']:
                    warnings.append(f"Reverse: {rev_specificity['warning']}")

            # template target-site accessibility (별개 축, 경고만). primer가 붙을 template 자리가
            #   접혀 접근 저해되는지 — template은 못 바꾸는 축이라 fail 아님(첨가제/온도 권장).
            #   fwd 결합부위 = [fwd_start, fwd_start+len(fwd)] / rev = [rev_start-len+1, rev_start+1].
            fwd_access = analyze_target_accessibility(template, fwd_start, fwd_start + len(fwd_seq))
            rev_access = analyze_target_accessibility(
                template, max(0, rev_start - len(rev_seq) + 1), rev_start + 1)
            if fwd_access['warning']:
                warnings.append(f"Forward 결합부위: {fwd_access['warning']}")
            if rev_access['warning']:
                warnings.append(f"Reverse 결합부위: {rev_access['warning']}")

            return {
                'forward': fwd_seq,
                'reverse': rev_seq,
                'fwd_tm': fwd_tm,
                'rev_tm': rev_tm,
                'fwd_gc': fwd_gc,
                'rev_gc': rev_gc,
                'fwd_specificity': fwd_specificity,
                'rev_specificity': rev_specificity,
                'fwd_accessibility': fwd_access,
                'rev_accessibility': rev_access,
                'amplicon_size': amplicon_size,
                'tier': 'yellow' if warnings else 'green',
                'warnings': warnings
            }

        except Exception as e:
            warnings.append(f"❌ primer3 실행 오류: {e}")
            return {
                'forward': None,
                'reverse': None,
                'warnings': warnings
            }

    @classmethod
    def design_goldengate_primers(
        cls, # TODO: cls가 뭔지 파악
        insert: Any,
        left_overhang: str,
        right_overhang: str,
        restriction_enzyme: str = "BsaI",
        target_start: Optional[int] = None,
        target_end: Optional[int] = None,
        params=None,
    ) -> Dict:
        """Golden Gate Assembly용 프라이머 설계 및 완성된 서열 생성

        Args:
            insert: 증폭할 DNA (Dseqrecord)
            left_overhang: 왼쪽 overhang (4bp)
            right_overhang: 오른쪽 overhang (4bp)
            restriction_enzyme: 사용할 제한효소 (기본값: "BsaI", "BbsI" 지원)
            target_start/target_end: 증폭 타겟 좌표(호출자 주입). 미지정 시 부품 전체.
            ...
        """
        from Bio.Seq import Seq

        # 설계 파라미터(연구실 조건). prefix/suffix 조립부터 쓰이므로 함수 맨 앞에서 확정한다.
        from params import DEFAULT_PARAMS, THERMO_MAX_LEN
        prm = params if params is not None else DEFAULT_PARAMS

        # 1. 타겟 좌표 (C-1′: mode-blind) — 좌표는 호출자(prep)가 주입. 미지정 시 부품 전체(biopart 기본).
        #    feature 기반 타겟 결정은 general-prep(Gibson)의 책임으로 분리(설계함수는 좌표만 받음).
        if target_start is None:
            target_start = 0
        if target_end is None:
            target_end = len(insert.seq)

        # 2. 제한효소 인식부위 및 Spacer 설정
        # TODO: sequence_analysis의 master_enzyme과 연계해서 처리할 수는 없나?
        # B-1: 하드코딩 제거 → Biopython 단일 소스에서 동적 도출
        # 지원 효소는 params.SUPPORTED_ENZYMES 한 곳에서 관리한다(값은 Biopython이 준다).
        from params import SUPPORTED_ENZYMES
        enzyme_sites = {name: get_typeIIS_params(name) for name in SUPPORTED_ENZYMES}

        # ⚠️ 구버전은 모르는 효소가 오면 **말없이 BsaI로 바꿔** 설계했다.
        #    유저가 BsmBI를 고르고도 BsaI 프라이머를 받는, 알아채기 어려운 오류라 에러로 바꾼다.
        if restriction_enzyme not in enzyme_sites:
            raise ValueError(
                f"지원하지 않는 제한효소: {restriction_enzyme!r} "
                f"(지원: {list(SUPPORTED_ENZYMES)})"
            )

        site_info = enzyme_sites[restriction_enzyme]
        site_fwd = site_info["fwd"]
        site_rev = site_info["rev"]

        # Spacer 설정 (B-2): 효소별 고정 더미. get_typeIIS_params의 def_spacer가
        # BsaI='A'(1bp)/BbsI='AA'(2bp)로 spacer_len에 맞춰 산출됨. spacer는 절단 후 소멸하므로
        # vector에서 추출할 필요 없음 → _extract_spacer_from_vector 제거.
        spacer = site_info["def_spacer"]

        # 3. Amplicon 및 서열 조립 요소 준비
        target_seq = str(insert.seq[target_start:target_end])

        # 오버행 RC 처리 (절단 후 상보적 결합용)
        left_oh_seq = left_overhang.upper() if left_overhang else ""
        right_oh_rc = str(Seq(right_overhang).reverse_complement()).upper() if right_overhang else ""

        # 4. 5' 부착부(prefix/suffix) = 효소 site + spacer + overhang. 양끝 blunt(overhang 없음)은 부착부 없음.
        #   버그2 수정: suffix는 top strand에 쌓이므로 reverse 프라이머(=RC)에서 spacer가 RC(A)=T로
        #   보였다(표현 비대칭). suffix에 RC(spacer)를 쌓으면 RC 후 spacer가 A로 복원되어 양쪽
        #   프라이머 모두 site-A-overhang으로 통일된다. spacer nt는 절단 후 소멸하는 더미라
        #   조립 산물(SHA)은 불변 — reverse 프라이머/full_sequence 출력 표현만 정정된다.
        spacer_rc = str(Seq(spacer).reverse_complement())

        # 5' 여분 염기(pad): 제한효소는 인식서열이 단편 **말단**에 있으면 잘 못 자른다
        # (NEB "Cleavage Close to the End of DNA Fragments"). 구버전은 프라이머가
        # 인식서열로 바로 시작해 여분이 0bp였다 → PCR 산물 절단 효율 저하.
        #   forward: [pad][site][spacer][left_oh] ...
        #   reverse: 역상보이므로 suffix **끝에** RC(pad)를 붙여야 rev 프라이머가 [pad][site]로 시작한다.
        # pad는 절단 시 site·spacer와 함께 떨어져 나가므로 **조립 산물 서열은 불변**이다.
        pad = prm.resolve_site_pad(restriction_enzyme)
        pad_rc = str(Seq(pad).reverse_complement()) if pad else ""

        prefix = (pad + site_fwd + spacer + left_oh_seq) if left_oh_seq else ""
        suffix = (right_oh_rc + spacer_rc + site_rev + pad_rc) if right_oh_rc else ""

        # 5. C-3: amplicon 길이 자동구제. 두 Tm[① amplicon-only(1차 cycle 결합부) ② 전체 프라이머(2차~ 5'tail 포함)]
        #    차이를 TM_DIFF_MAX 이하로 수렴시키도록 amplicon을 늘린다.
        #    full_sequence는 target 전체 기반이라 amplicon 조정은 조립 서열에 영향 없음(프라이머 품질만 개선).
        #    [상한] 프라이머(=tail+amplicon)는 primer3 열역학 상한(PRIMER_MAX_LEN=60bp)을 넘지 않는다
        #      (60bp 초과 데이터를 primer3에 넣지 않음 — calc_hairpin/dimer 계산 불가 구간).
        #    [미수렴 정책] 60bp 안에서 TM_DIFF_MAX를 못 맞춰도 실패로 죽이지 않고,
        #      diff가 가장 작은(컷오프에 제일 근접한) 최적 amplicon을 채택한다(살짝 못 미친 것도 사용 가능).
        TM_DIFF_MAX = prm.tm_diff_max
        AMP_MIN = prm.amp_min
        # 터미널(tail 없는 끝단) 프라이머의 절대 Tm 목표 — 연구실 파라미터
        TERMINAL_TM_MIN = prm.terminal_tm_min
        # primer3 calc_hairpin/dimer 입력 상한(초과 시 계산 불가) — 실측 확인된 하드 리밋
        PRIMER_MAX_LEN = THERMO_MAX_LEN

        def _tm(seq: str) -> float:
            """이 설계 한 벌에 일관된 조건으로 Tm 계산(연구실 염/농도)."""
            return calculate_tm(seq, **prm.tm_kwargs())

        def _fit_amplicon(pool: str, tail: str, from_end: bool) -> dict:
            in_window = []   # 선언 범위 안에 든 후보들 (tail 있는 경우)
            # 프라이머 길이 = len(tail) + L 이므로 L ≤ 60-len(tail)로 amplicon 상한을 건다.
            amp_cap = max(1, PRIMER_MAX_LEN - len(tail))
            amp_max = min(len(pool), amp_cap)
            start = min(AMP_MIN, amp_max)
            best = None   # 목표 미달 시 채택할 최적 근사 후보
            for L in range(max(1, start), amp_max + 1):
                amp = pool[-L:] if from_end else pool[:L]
                full = (amp + tail) if from_end else (tail + amp)
                tm_amp = _tm(amp)                       # anneal(주형 결합부) Tm
                tm_full = _tm(full) if tail else tm_amp  # 전체 프라이머(5' tail 포함) Tm
                cand = {"amp": amp, "tm_amp": tm_amp, "tm_full": tm_full,
                        "diff": abs(tm_full - tm_amp), "len": L, "adjusted": L > start}
                if tail:
                    # ── 선정 규칙 (2026-09-23 재설계) ────────────────────────────
                    # ① **anneal Tm이 선언한 범위 [tm_min, tm_max] 안에 들 것** — 이것이 1순위.
                    #    1차 cycle은 anneal 구간만으로 결합하므로 이 값이 실제 Ta를 좌우한다.
                    # ② anneal Tm ↔ 전체 프라이머 Tm의 차(diff)는 **선호도**다(하드 게이트 아님).
                    #
                    # ⚠️ 왜 diff를 게이트에서 내렸나 (실측 근거):
                    #    GG 프라이머는 [고정 tail][anneal] 구조이고 tail(사이트+overhang)은
                    #    길이를 못 줄인다. diff를 5℃ 안으로 **먼저** 좁히려면 anneal을 길게
                    #    가야 하고(꼬리 기여를 희석해야 하므로), 그러면 Tm이 범위를 넘어간다.
                    #    그 결과 spacer 115개가 **전부** 판정에서 떨어져 GG 삽입 경로가
                    #    통째로 막혀 있었다(2026-09-23 실측 115/115).
                    #    반대로 창을 1순위로 두면 **어떤 범위를 선언해도 115/115가 들어간다.**
                    #    diff는 그때 결과로 따라 나오는 값이라 보고만 한다(리포트 노출).
                    #    낮은 범위를 선언할수록 diff가 커지는 것은 tail이 있는 한 물리적 귀결이다.
                    cand["tm_ok"] = prm.tm_min <= cand["tm_amp"] <= prm.tm_max
                    cand["diff_ok"] = cand["diff"] <= TM_DIFF_MAX
                    if cand["tm_ok"]:
                        in_window.append(cand)
                    # 창에 하나도 못 들 때를 대비한 최적 근사: 창 중앙에 가까운 것
                    _miss = abs(cand["tm_amp"] - prm.tm_opt)
                    if best is None or _miss < abs(best["tm_amp"] - prm.tm_opt):
                        best = cand
                else:
                    # 터미널 프라이머(tail 없음, 맞출 overhang 없는 끝단): amplicon==full이라
                    #   diff는 항상 0이라 무의미 → 절대 Tm(TERMINAL_TM_MIN)에 도달할 때까지 안쪽으로 연장.
                    #   (안 그러면 spacer 끝 15bp를 그대로 복제 → Tm이 낮은 짧은 프라이머가 됨)
                    if tm_amp >= TERMINAL_TM_MIN:
                        cand["ok"] = True
                        return cand   # 목표 도달: 그 길이 채택
                    if best is None or cand["tm_amp"] > best["tm_amp"]:
                        best = cand   # 미달: 가장 높은 Tm 후보(=가장 긴 것) 갱신
            if in_window:
                # 창 안에서 고른다: diff 컷을 만족하는 것 중 **가장 짧은**(합성 저렴) 것,
                # 없으면 diff가 가장 작은 것. (창 안이면 어느 쪽이든 Tm은 목표 범위다)
                good = [c for c in in_window if c["diff_ok"]]
                pick = good[0] if good else min(in_window, key=lambda c: c["diff"])
                pick["ok"] = True
                return pick
            if best is None:  # 빈 target 방어
                best = {"amp": "", "tm_amp": 0.0, "tm_full": 0.0, "diff": 0.0,
                        "len": 0, "adjusted": False}
            best["ok"] = False   # 목표 미달 → 최적 근사 프라이머(비차단 경고)
            return best

        fwd_fit = _fit_amplicon(target_seq, prefix, from_end=False)
        rev_fit = _fit_amplicon(target_seq, suffix, from_end=True)
        amplicon_fwd = fwd_fit["amp"]
        amplicon_end = rev_fit["amp"]

        # 6. 구조 설계 (Top strand 기준). overhang 없는 linear 양끝은 amplicon만(blunt) — C-2.
        start_structure = (prefix + amplicon_fwd) if left_oh_seq else amplicon_fwd
        end_structure = (amplicon_end + suffix) if right_oh_rc else amplicon_end

        # 7. 프라이머 구성 (Reverse 프라이머는 전체의 RC)
        fwd_primer = start_structure
        rev_primer = str(Seq(end_structure).reverse_complement())

        # 8. 완성된 PCR 부품 서열 (제한효소 cut 직전 전체) — target 전체 기반, amplicon 무관
        full_sequence = prefix + target_seq + suffix

        # 9. 프라이머 품질 검증(dimer/hairpin/GC) + Tm 수렴 결과 종합 (C-3, 비차단)
        validator = ValidatePrimer()
        # 5' 부착부(prefix=site+spacer+overhang / suffix RC)는 tail → length/Tm/GC는 amplicon만 검증.
        fwd_val = validator.validate_primer(fwd_primer, check_dimer_with=rev_primer,
                                            tail_len=len(prefix), params=prm)
        rev_val = validator.validate_primer(rev_primer, check_dimer_with=fwd_primer,
                                            tail_len=len(suffix), params=prm)

        all_violations = []
        if not fwd_val['passed']:
            all_violations.extend([f"Fwd: {v}" for v in fwd_val['violations']])
        if not rev_val['passed']:
            all_violations.extend([f"Rev: {v}" for v in rev_val['violations']])
        # C-3 Tm 수렴: 조정 로그 + 60bp 내 미수렴 시 최적 근사 프라이머를 경고로 안내(비차단).
        #   (하드 실패 아님 — 60bp 안 최적 근사를 채택. 살짝 못 미친 것도 사용 가능)
        tm_warnings = []
        for side, fit in (("Fwd", fwd_fit), ("Rev", rev_fit)):
            if fit["adjusted"]:
                print(f"     🔧 [{side}] amplicon {AMP_MIN}→{fit['len']}bp 조정 "
                      f"(Tm amplicon {fit['tm_amp']:.1f} / 전체 {fit['tm_full']:.1f}, diff {fit['diff']:.1f}°C)")
            if not fit["ok"]:
                unmet = []
                if not fit.get("tm_ok", True):
                    unmet.append(f"Tm 범위 [{prm.tm_min:.1f}, {prm.tm_max:.1f}]℃ 밖"
                                 f"(anneal {fit['tm_amp']:.1f}℃)")
                if not fit.get("diff_ok", True):
                    unmet.append(f"anneal↔전체 Tm 차 {fit['diff']:.1f}℃ > {TM_DIFF_MAX}℃")
                detail = " / ".join(unmet) if unmet else f"diff {fit['diff']:.1f}°C"
                tm_warnings.append(
                    f"{side}: 60bp 안에서 목표를 못 맞춰 최적 근사 채택 — {detail} "
                    f"(anneal {fit['tm_amp']:.1f}℃ / 전체 {fit['tm_full']:.1f}℃)")
                print(f"     ⚠️ [{side}] 60bp 내 최적 근사 프라이머 사용 — {detail}")

        # 10. template accessibility (별개 축, 경고만): 말단 고정 binding site(amplicon 양 말단)가
        #   insert 서열 내에서 접혀 접근 저해되는지. 위치는 못 바꾸는 축(말단 고정) → fail 아님.
        #   fwd binding = insert[target_start : +fwd_len] / rev binding = insert[target_end-rev_len : ].
        all_warnings = []
        tmpl_seq = str(insert.seq)
        fwd_access = analyze_target_accessibility(
            tmpl_seq, target_start, target_start + fwd_fit['len'])
        rev_access = analyze_target_accessibility(
            tmpl_seq, max(target_start, target_end - rev_fit['len']), target_end)
        if fwd_access['warning']:
            all_warnings.append(f"Fwd 결합부위: {fwd_access['warning']}")
        if rev_access['warning']:
            all_warnings.append(f"Rev 결합부위: {rev_access['warning']}")
        # 프라이머 자체 검증의 비차단 경고(tail/노출ss hairpin·dimer)도 합류
        all_warnings.extend([f"Fwd: {w}" for w in fwd_val.get('warnings', [])])
        all_warnings.extend([f"Rev: {w}" for w in rev_val.get('warnings', [])])
        all_warnings.extend(tm_warnings)   # C-3 Tm 미수렴(60bp 내 최적 근사) 경고 합류
        # 쌍 ΔTm — GG는 길이가 overhang·spacer에 묶여 Tm을 맞출 자유도가 적어
        # 실측 최악 케이스(11.68℃)가 이 경로에 몰렸다(§check_pair_tm_diff).
        _pair_w = check_pair_tm_diff(fwd_val['tm'], rev_val['tm'], prm)
        if _pair_w:
            all_warnings.append(_pair_w)

        # 구조 수치(dimer/hairpin ΔG) — tier에 영향 없이 리포트로만 나간다.
        all_structural = []
        all_structural.extend([f"Fwd: {n}" for n in fwd_val.get('structural_notes', [])])
        all_structural.extend([f"Rev: {n}" for n in rev_val.get('structural_notes', [])])

        # 전체 tier(🔴 red=차단 violation / 🟡 yellow=경고만 / 🟢 green): 3-tier 요약.
        tier = 'red' if all_violations else ('yellow' if all_warnings else 'green')

        return {
            'forward': fwd_primer,
            'reverse': rev_primer,
            'start_structure': start_structure,
            'end_structure': end_structure,
            'full_sequence': full_sequence,
            'target_seq': target_seq,
            'restriction_enzyme': restriction_enzyme,
            'fwd_tm': fwd_val['tm'],
            'rev_tm': rev_val['tm'],
            'fwd_gc': fwd_val['gc_content'],   # anneal(amplicon) 기준 GC% (0~1)
            'rev_gc': rev_val['gc_content'],
            'fwd_hairpin_dg': fwd_val['hairpin_dg_full'],  # 프라이머 self-hairpin ΔG(kcal/mol)
            'rev_hairpin_dg': rev_val['hairpin_dg_full'],
            'fwd_tm_amplicon': fwd_fit['tm_amp'],
            'rev_tm_amplicon': rev_fit['tm_amp'],
            # 조립용 프라이머는 [5' tail][anneal] 구조라 두 Tm이 크게 다르다.
            # 어느 쪽을 말하는지 혼동되면 연구실 기준과 어긋나므로 **둘 다** 내보낸다.
            #   *_tm / *_tm_amplicon : anneal(주형 결합부) 기준
            #   *_tm_full            : 주문할 프라이머 전체 기준
            'fwd_tm_full': fwd_fit['tm_full'],
            'rev_tm_full': rev_fit['tm_full'],
            # 이 설계에 쓰인 Tm 계산 조건(연구실 조건 추적용)
            'tm_na_mM': prm.na_mM,
            'tm_primer_nM': prm.primer_conc_nM,
            'site_5prime_pad': pad,
            'fwd_tm_diff': fwd_fit['diff'],
            'rev_tm_diff': rev_fit['diff'],
            'amplicon_fwd_len': fwd_fit['len'],
            'amplicon_rev_len': rev_fit['len'],
            # 2026-11-05(Y12) — **주형 결합부의 프라이머 서열**. 교차 dimer 검사가 쓴다.
            #   fwd_primer = prefix + amplicon_fwd      → 결합부 = amplicon_fwd
            #   rev_primer = RC(amplicon_end + suffix)  → 결합부 = RC(amplicon_end)
            #   ⚠️ 호출자가 프라이머를 슬라이스해 되찾게 두지 않는다 — Y10-b에서 같은
            #      실수를 했다(`_validate_primer_length`가 길이를 바꿔 슬라이스가 어긋났다).
            'fwd_anneal': amplicon_fwd,
            'rev_anneal': str(Seq(amplicon_end).reverse_complement()),
            'fwd_accessibility': fwd_access,
            'rev_accessibility': rev_access,
            'fwd_passed': fwd_val['passed'],   # fwd(=L 발판) 단독 통과 — spacer 선정 L 후보 판정용
            'rev_passed': rev_val['passed'],   # rev(=R 발판) 단독 통과 — spacer 선정 R 후보 판정용
            'passed': len(all_violations) == 0,
            'tier': tier,
            'violations': all_violations,
            'warnings': all_warnings,
            'structural_notes': all_structural,
        }



class ValidatePrimer:
    @staticmethod
    def check_primer_dimer(
        primer1: str,
        primer2: str,
        threshold_dg: float = DIMER_DG_CUTOFF_KCAL
    ) -> Dict:
        """
        프라이머 dimer 체크 (primer3-py 기반)

        primer3의 thermodynamic 계산을 사용하여 dimer를 정확히 체크합니다.
        - Homodimer (self-dimer): 같은 프라이머끼리
        - Heterodimer: 다른 프라이머끼리

        Args:
            primer1: 첫 번째 프라이머
            primer2: 두 번째 프라이머
            threshold_dg: Dimer 판단 임계값 (ΔG, kcal/mol, 기본값: -5.0)
                - ΔG < -5: 강한 dimer (문제)
                - -5 ≤ ΔG < -3: 중간 dimer (주의)
                - ΔG ≥ -3: 약한 dimer (괜찮음)

        Returns:
            {
                'has_dimer': bool,              # Dimer 존재 여부
                'delta_g': float,               # ΔG (kcal/mol)
                'tm': float,                    # Dimer Tm (°C)
                'warning': str                  # 경고 메시지
            }

        Example:
            >>> p1 = "ATCGATCGATCGATCG"
            >>> p2 = "CGATCGATCGATCGAT"
            >>> result = check_primer_dimer(p1, p2)
            >>> result['has_dimer']
            True
            >>> result['delta_g'] < -5
            True
        """
        if not primer1 or not primer2:
            raise ValueError("프라이머가 비어있습니다")

        primer1 = primer1.upper()
        primer2 = primer2.upper()

        # 하이브리드 경로(H1): primer3는 한쪽이라도 60bp를 넘으면 거부한다.
        # 조립용 프라이머는 [5' tail][anneal] 구조라 arm이 길어지면 이 구간에 쉽게 들어가고,
        # 그러면 dimer 축이 통째로 눈을 감는다(기존에는 에러 문구만 찍고 ΔG 0을 돌려줬다).
        # ⚠️ ViennaRNA 스케일은 primer3와 달라 **경고 전용**이다(has_dimer=False).
        from params import THERMO_MAX_LEN

        if max(len(primer1), len(primer2)) > THERMO_MAX_LEN:
            from basic.secondary_structure import analyze_long_duplex

            long_res = analyze_long_duplex(primer1, primer2)
            if long_res is not None:
                return long_res
            # ViennaRNA도 없으면 아래 primer3 경로로 내려가 기존 폴백(에러 → ΔG 0)을 탄다.

        try:
            # Heterodimer 계산
            result = primer3.bindings.calc_heterodimer(
                primer1,
                primer2,
                mv_conc=50.0,  # Na+ 농도
                dv_conc=0.0,   # Mg2+ 농도
                dntp_conc=0.0,
                dna_conc=250.0
            )
            # TODO: 위의 기본값 세팅의 기준 evidence

            delta_g = result.dg / 1000.0  # cal/mol → kcal/mol
            tm = result.tm

            # Dimer 판단
            has_dimer = delta_g < threshold_dg

            # 경고 메시지
            warning = ""
            if delta_g < -9:
                warning = f"⚠️⚠️  매우 강한 dimer (ΔG = {delta_g:.2f} kcal/mol, Tm = {tm:.1f}°C)"
            elif delta_g < -7:
                warning = f"⚠️  강한 dimer (ΔG = {delta_g:.2f} kcal/mol, Tm = {tm:.1f}°C)"
            elif delta_g < -5:
                warning = f"⚠️  중간 dimer (ΔG = {delta_g:.2f} kcal/mol, Tm = {tm:.1f}°C)"
            elif delta_g < -3:
                warning = f"약한 dimer (ΔG = {delta_g:.2f} kcal/mol, Tm = {tm:.1f}°C)"

            return {
                'has_dimer': has_dimer,
                'delta_g': delta_g,
                'tm': tm,
                'warning': warning
            }

        except Exception as e:
            print(f"❌ primer3 dimer 계산 오류: {e}")
            return {
                'has_dimer': False,
                'delta_g': 0.0,
                'tm': 0.0,
                'warning': f'Error: {e}'
            }

    @staticmethod
    def check_gc_clamp(
        primer: str,
        min_gc: int = 1,
        max_gc: int = 3
    ) -> Dict:
        """
        3' 말단 GC clamp 확인

        프라이머 3' 말단에 적절한 GC가 있는지 확인합니다.
        GC clamp는 프라이머 결합 안정성을 높입니다.

        Args:
            primer: 프라이머 서열
            min_gc: 최소 GC 개수 (기본값: 1)
            max_gc: 최대 GC 개수 (기본값: 3)

        Returns:
            {
                'has_gc_clamp': bool,       # GC clamp 존재 여부
                'gc_count': int,            # 3' 말단 5bp 내 GC 개수
                'last_3_bases': str,        # 마지막 3bp
                'warning': str              # 경고 메시지
            }

        Example:
            >>> check_gc_clamp("ATCGATCGGC")['has_gc_clamp']
            True
            >>> check_gc_clamp("ATCGATCGAA")['has_gc_clamp']
            False
        """
        if not primer:
            raise ValueError("프라이머가 비어있습니다")

        primer = primer.upper()

        # 3' 말단 5bp 내 GC 개수 계산
        terminal_5bp = primer[-5:]
        gc_count = terminal_5bp.count('G') + terminal_5bp.count('C')

        # 마지막 3bp
        last_3 = primer[-3:]

        # GC clamp 존재 여부
        has_gc_clamp = min_gc <= gc_count <= max_gc

        # 경고 메시지
        warning = ""
        if gc_count == 0:
            warning = "⚠️  3' 말단에 GC가 없습니다 (불안정한 결합)"
        elif gc_count > max_gc:
            warning = f"⚠️  3' 말단 GC가 너무 많습니다 ({gc_count}개, 권장: {min_gc}-{max_gc}개)"
        elif gc_count < min_gc:
            warning = f"⚠️  3' 말단 GC가 부족합니다 ({gc_count}개, 권장: {min_gc}-{max_gc}개)"

        return {
            'has_gc_clamp': has_gc_clamp,
            'gc_count': gc_count,
            'last_3_bases': last_3,
            'warning': warning
        }

    @staticmethod
    def check_primer_specificity(
        template: str,
        primer: str,
        target_start: int,
        target_end: int,
        is_reverse: bool = False
    ) -> Dict:
        """
        프라이머 특이성 체크 (Off-target binding)

        프라이머가 타겟 영역 외에 다른 곳에 결합하는지 확인합니다.

        Args:
            template: 전체 Template DNA 서열
            primer: 프라이머 서열
            target_start: 타겟 시작 위치
            target_end: 타겟 종료 위치
            is_reverse: Reverse 프라이머 여부

        Returns:
            {
                'is_specific': bool,            # 특이성 있음
                'target_match': bool,           # 타겟 위치에 결합
                'off_target_count': int,        # Off-target 결합 수
                'off_target_positions': List[int],  # Off-target 위치들
                'warning': str                  # 경고 메시지
            }

        Example:
            >>> template = "ATCGATCG" * 20
            >>> primer = "ATCGATCG"
            >>> result = check_primer_specificity(template, primer, 0, 50, False)
            >>> result['target_match']
            True
        """
        if not template or not primer:
            raise ValueError("서열이 비어있습니다")

        from Bio.Seq import Seq

        template = template.upper()
        primer = primer.upper()

        # Reverse 프라이머는 reverse complement로 변환
        if is_reverse:
            primer_seq = str(Seq(primer).reverse_complement())
        else:
            primer_seq = primer

        # 타겟 영역 확인
        target_region = template[target_start:target_end]

        # 프라이머가 template에서 결합하는 모든 위치 찾기
        positions = []
        search_start = 0

        while True:
            pos = template.find(primer_seq, search_start)
            if pos == -1:
                break
            positions.append(pos)
            search_start = pos + 1

        # 타겟 영역과 가장 가까운 위치 찾기
        if not positions:
            # 프라이머가 template에 없음
            target_match = False
            off_target_positions = []
        else:
            # Forward: 타겟 시작점에 가장 가까운 위치
            # Reverse: 타겟 끝점에 가장 가까운 위치
            if is_reverse:
                # Reverse complement가 타겟 끝 근처에 있어야 함
                closest_pos = min(positions, key=lambda p: abs(p - target_end))
                closest_distance = abs(closest_pos - target_end)
            else:
                # Forward가 타겟 시작 근처에 있어야 함
                closest_pos = min(positions, key=lambda p: abs(p - target_start))
                closest_distance = abs(closest_pos - target_start)

            # 가장 가까운 위치가 타겟 매치, 나머지는 off-target
            target_match = True  # 일단 프라이머가 template에 존재하면 True
            off_target_positions = [p for p in positions if p != closest_pos]

        off_target_count = len(off_target_positions)
        is_specific = target_match and off_target_count == 0

        # 경고 메시지
        warning = ""
        if not target_match:
            warning = f"⚠️⚠️  프라이머가 template에서 찾을 수 없습니다!"
        elif off_target_count > 0:
            warning = f"⚠️  Off-target 결합 {off_target_count}개 발견 (위치: {off_target_positions[:3]}{'...' if off_target_count > 3 else ''})"

        return {
            'is_specific': is_specific,
            'target_match': target_match,
            'off_target_count': off_target_count,
            'off_target_positions': off_target_positions,
            'warning': warning
        }

    @staticmethod
    def validate_primer(
        primer: str,
        min_tm: Optional[float] = None,
        max_tm: Optional[float] = None,
        min_gc: Optional[float] = None,
        max_gc: Optional[float] = None,
        check_dimer_with: Optional[str] = None,
        tail_len: int = 0,
        exposed_ss: Optional[str] = None,
        params=None,
    ) -> Dict:
        """
        프라이머 종합 검증

        프라이머의 품질을 종합적으로 검증합니다.

        검증 항목:
        - Tm 범위 (`tm_min`~`tm_max` — **유저 선언창**, 기본 55~65℃)
        - GC% 범위 (`gc_min`~`gc_max` — 유저 선언창, 기본 20~80%)
        - 결합부 길이 (`primer_len_min`~`primer_len_max` — 유저 선언창)
        - Self-dimer
        - Hairpin (2차 구조)
        - GC clamp

        Args:
            primer: 프라이머 서열
            min_tm: 최소 Tm (°C). None이면 `params.tm_min`(**유저 선언창** — 2026-11-05
                Y11에서 판정창을 없애고 선언창 하나로 합쳤다)
            max_tm: 최대 Tm (°C). None이면 `params.tm_max`
            min_gc: GC 하한 **분수**(0~1). None이면 `params.gc_min`(%)을
                100으로 나눠 쓴다. 차단창은 설계창(`gc_min`/`gc_max`)을 포함한다 —
                그 안쪽의 이상 구간 이탈은 차단이 아니라 경고로 나간다(아래)
            max_gc: GC 상한 **분수**(0~1). None이면 `params.gc_max`(%)/100
            check_dimer_with: 다른 프라이머 (dimer 체크용)
            tail_len: 5' 비-anneal tail 길이 (기본값: 0 = 일반 PCR).
                Gibson(vector homology)/GG(site+spacer+overhang)처럼 5'에 tail이 붙어
                프라이머 전체가 길어지는 경우:
                - length/Tm/GC는 anneal 부분(primer[tail_len:], 1차 cycle 결합부)으로만 검증
                  (GG design_goldengate_primers의 amplicon-only vs full Tm 관점과 일관)
                - self/cross dimer는 tail 유래(vector 의존, spacer로 해결 불가)가 지배적이라
                  **경고로만**(violations 제외, has_* 플래그는 유지).
                - hairpin은 출처를 나눠 처리한다(버그3 층2 재설계):
                    · anneal(spacer/증폭 영역) hairpin → 다른 spacer로 교체 가능 = **차단(violation)**
                    · exposed_ss(아래) hairpin → vector/method 의존 = **경고**('warnings')
                    · 위 둘로 설명 안 되는 전체-프라이머 구조 → **경고**(PCR 단계 신호)
                tail_len=0(순수 PCR)이면 anneal==full → hairpin/dimer 모두 기존대로 차단(하위호환).
            exposed_ss: ② 조립 시 외부로 노출되어 partner와 anneal해야 하는 말단 ss 구간 (기본 None).
                caller가 조립 method를 알고 잘라 넘긴다:
                - Gibson: T5가 노출시킨 homology 구간(보통 5' tail=vector homology) → primer[:tail_len]
                - Golden Gate: 절단 후 4nt sticky end라 접힘 무의미 → None(계산 생략)
                제공 시 primer3로 **직접** hairpin ΔG를 재서(추론 X) 강하면 경고에 명시한다.

        Returns:
            {
                'passed': bool,                 # 전체 검증 통과 여부
                'length': int,                  # 프라이머 길이
                'tm': float,                    # Tm (°C)
                'gc_content': float,            # GC% (0-1)
                'has_gc_clamp': bool,           # GC clamp 존재
                'has_self_dimer': bool,         # Self-dimer 존재
                'has_hairpin': bool,            # Hairpin 존재(전체 프라이머 기준)
                'has_hairpin_3prime': bool,     # 전체 프라이머 hairpin이 3'말단 관여
                'hairpin_dg_full': float,       # ① 전체 프라이머 self-hairpin ΔG(kcal/mol)
                'hairpin_dg_anneal': float,     #   anneal(교체가능 결합부) hairpin ΔG
                'hairpin_dg_exposed': Optional[float],  # ② 노출 말단 ss hairpin ΔG(미제공 시 None)
                'violations': List[str],        # 차단 위반 항목 리스트
                'warnings': List[str]           # 비차단 경고(노출 ss/tail hairpin, dimer 등)
            }

        Example:
            >>> primer = "TGCACCACCAACTGCTTAGC"  # Good primer
            >>> result = validate_primer(primer)
            >>> result['passed']
            True
            >>> len(result['violations']) == 0
            True
        """
        if not primer:
            raise ValueError("프라이머가 비어있습니다")

        primer = primer.upper()
        # anneal 부분 = 5' tail을 뺀 3' 영역(1차 cycle에 template 결합). tail_len=0이면 전체(일반 PCR).
        anneal = primer[tail_len:] if tail_len else primer
        violations = []
        warnings = []   # 비차단 경고(예: tail 유래 hairpin — spacer로 해결 불가한 vector/method 신호)

        # tail 프라이머(GG/Gibson, tail_len>0)는 anneal 영역이 설계상 짧아 length/Tm/GC가 낮게 나오는 게
        # 정상이다(1차 cycle만 3' 결합, 이후 full primer가 증폭 주도). 따라서 anneal 기준 length/Tm/GC는
        # 하드 차단(violation=red)이 아니라 경고(yellow)로 둔다 — dimer/hairpin의 tail 처리와 일관.
        # 순수 PCR(tail_len=0)은 anneal==full이라 기존대로 violation(red) 유지.
        anneal_gate = warnings if tail_len else violations
        _t = " (tail — 경고)" if tail_len else ""

        # ── 판정 기준은 **유저가 선언한 값** 하나다 (2026-11-05 Y11) ─────────────
        # 예전에는 세 축 모두 선언창 바깥에 판정창을 하나 더 뒀다(`qc_tm_*` 45~85℃,
        # `qc_gc_*` 20~80% + `qc_gc_ideal_*` 40~60%). 길이 축은 더 나빴다 — 판정값이
        # `15`/`30` **리터럴**이라 유저가 Q5(상한 26nt)를 골라도 경고는 30 기준으로
        # 떴다. 유저 설정과 경고가 따로 노는 상태였고, 그 결과 78행에 경고 154건
        # (행당 2.0건) 중 **대부분이 유저가 설정한 적 없는 기준** 위반이었다.
        #
        # → 판정도 선언창으로 한다. "내가 설정한 범위를 벗어나면 알려준다"가 되고,
        #   창이 하나뿐이라 어긋남이 구조적으로 불가능하다(§params.TM_INPUT_MIN).
        #   물리 경계는 입력 제한이 막으므로 여기서 또 막을 필요가 없다.
        #   호출자가 값을 직접 주면 그것이 우선한다(None일 때만 params에서 가져온다).
        from params import DEFAULT_PARAMS
        _p = params or DEFAULT_PARAMS
        min_tm = _p.tm_min if min_tm is None else min_tm
        max_tm = _p.tm_max if max_tm is None else max_tm
        # 선언창 gc_*는 **%**, 이 함수 내부 비교는 **분수**라 여기서 한 번만 환산한다.
        min_gc = (_p.gc_min / 100.0) if min_gc is None else min_gc
        max_gc = (_p.gc_max / 100.0) if max_gc is None else max_gc

        # ── 1. 길이 체크 — **tail 유무로 자가 갈린다** (2026-11-05 Y11) ───────────
        # 이 분기는 새 규칙이 아니다. 바로 아래 `anneal_gate`가 이미 같은 `tail_len`으로
        # 적/황을 가르고 있었고, 길이 축만 그 규칙을 안 따르고 있었다.
        #
        #   · tail_len == 0 (일반 PCR·끝단 프라이머)
        #       결합부 길이 = 프라이머 전체 길이이고, primer3가 Tm 창 안에서 길이를
        #       **자유롭게** 고를 수 있다. 그러니 유저 선언창(`primer_len_min/max`,
        #       GUI의 "Annealing length")을 지킬 수 있고, 못 지켰으면 알릴 값이 있다.
        #
        #   · tail_len > 0 (GG / Gibson 삽입 프라이머)
        #       tail이 길이·서열 **둘 다 고정**이다(GG는 [패딩][인식서열][spacer][overhang],
        #       Gibson은 벡터 homology 30bp). 그래서 결합부 길이는 유저가 고르는 값이
        #       아니라 `tm_diff_max` 수렴의 **종속 변수**다 — 실측하면 같은 설정에서도
        #       부품·프리셋에 따라 13~43nt로 흩어진다. 이 길이를 선언창으로 재면
        #       **유저가 숫자를 어떻게 바꿔도 못 없애는 경고**가 된다(실측 36건).
        #       → 여기서는 **실제로 막히는 제약**만 본다: 합성 주문 상한
        #         (`PRIMER_ORDER_MAX_NT`)과 결합부 특이성 바닥(`amp_min`).
        from params import PRIMER_ORDER_MAX_NT
        length = len(anneal)
        if tail_len:
            # 주문 단위는 **전체 프라이머**다(tail 포함). 결합부만 재면 못 잡는다.
            if len(primer) > PRIMER_ORDER_MAX_NT:
                anneal_gate.append(
                    f"주문할 프라이머가 합성 상한보다 깁니다 "
                    f"({len(primer)} nt > {PRIMER_ORDER_MAX_NT} nt){_t}")
            if length < _p.amp_min:
                anneal_gate.append(
                    f"주형 결합부가 특이성 하한보다 짧습니다 "
                    f"({length} bp < {_p.amp_min} bp — 짧은 부품이라 결합부가 tail에 "
                    f"흡수된 경우입니다){_t}")
        else:
            if length < _p.primer_len_min:
                anneal_gate.append(
                    f"프라이머가 설정한 길이 하한보다 짧습니다 "
                    f"({length} bp < {_p.primer_len_min} bp){_t}")
            elif length > _p.primer_len_max:
                anneal_gate.append(
                    f"프라이머가 설정한 길이 상한보다 깁니다 "
                    f"({length} bp > {_p.primer_len_max} bp){_t}")

        # 2. Tm 체크 — anneal 부분 기준 (1차 cycle 결합부)
        tm = calculate_tm(anneal, **_p.tm_kwargs())
        if tm < min_tm:
            anneal_gate.append(
                f"Tm이 설정한 하한보다 낮습니다 ({tm:.1f}°C < {min_tm:g}°C){_t}")
        elif tm > max_tm:
            anneal_gate.append(
                f"Tm이 설정한 상한보다 높습니다 ({tm:.1f}°C > {max_tm:g}°C){_t}")

        # 3. GC% 체크 — anneal 부분 기준, **선언창 1층**(2026-11-05 Y11)
        #    Y7에서 둔 경고층(`qc_gc_ideal_*` 40~60%)을 제거했다. 그 층은 유저가
        #    설정한 값이 아니면서 78행 중 39행(50%)에 경고를 달았고 — GC는 길이를
        #    늘려도 안 바뀌는 **템플릿의 성질**이라 경고를 봐도 실험자가 할 수 있는
        #    일이 없었다. 판정은 선언창(`gc_min/gc_max`) 하나로 한다.
        gc_content = calculate_gc_content(anneal)
        if gc_content < min_gc:
            anneal_gate.append(
                f"GC%가 설정한 하한보다 낮습니다 ({gc_content:.1%} < {min_gc:.1%}){_t}")
        elif gc_content > max_gc:
            anneal_gate.append(
                f"GC%가 설정한 상한보다 높습니다 ({gc_content:.1%} > {max_gc:.1%}){_t}")

        # 4. GC clamp 체크 (경고만, violations에 추가 안 함)
        gc_clamp_result = ValidatePrimer.check_gc_clamp(primer)
        has_gc_clamp = gc_clamp_result['has_gc_clamp']
        # GC clamp가 없어도 사용 가능하므로 경고만 출력

        # 5~7. dimer/hairpin (버그3 층2 재설계):
        #   [dimer] tail 프라이머(tail_len>0)의 self/cross dimer는 5' tail(vector homology /
        #     site+spacer+overhang) 유래가 지배적이고(P7 실증: fwd self-dimer 115/115 = vector tail
        #     공통 원인) spacer 선정으로 못 고치므로 → 경고로만 두고 차단 안 함(기존 정책 유지).
        #   [hairpin] 구버전은 dimer와 묶어 tail_len>0이면 통째로 껐지만, hairpin은 출처를 나눠야 한다:
        #     · anneal(spacer/증폭 영역) 자체 hairpin → 다른 spacer로 교체하면 풀림 = 선정으로 걸러야 함
        #       → **차단(violation)**. tail_len=0(순수 PCR)이면 anneal==full이라 동일하게 차단(하위호환).
        #     · tail 유래 hairpin(full엔 있으나 anneal엔 없음) → spacer로 해결 불가 = vector/method
        #       조합 신호 → **경고**(조용히 넘기지 않도록 warnings에 명시, 차단은 안 함).
        warn_only_dimer = tail_len > 0

        # ── ΔG는 화면 경고가 아니라 **내부 기록**으로 (2026-10-20) ──────────────
        # 실험자는 ΔG 숫자가 아니라 2차구조를 보고 "결합부가 뚫려 있나"로 판단한다
        # (연구실 실제 작업 방식). 그리고 ΔG 경고는 정보가 되지 못했다 — GG 프라이머는
        # 전부 같은 5' tail(`패딩+인식서열`)을 달고 있어서 tail만으로 이미 기준을 넘는다
        # (실측: tail 12bp 단독 self ΔG −7.44 / 역상보와 −10.82, 기준 −5.0).
        # 그 결과 51개 중 50개(98%)에 경고가 붙어 소음이 됐다.
        #
        # → 값은 버리지 않고 structural_notes에 남겨 리포트로 내보내고, tier 판정과
        #   화면 경고에서만 뺀다. 판정은 **결합부위 접근성**이 맡는다.
        structural_notes = []

        # 5. Self-dimer 체크
        self_dimer_result = ValidatePrimer.check_primer_dimer(
            primer, primer, threshold_dg=_p.qc_dimer_dg)
        has_self_dimer = self_dimer_result['has_dimer']
        if has_self_dimer and self_dimer_result['warning']:
            if warn_only_dimer:
                structural_notes.append(f"Self-dimer: {self_dimer_result['warning']}")
            else:
                violations.append(f"Self-dimer: {self_dimer_result['warning']}")

        # 6. Hairpin 체크 (2차 구조) — 세 관점(full/anneal/exposed)을 primer3로 직접 계산.
        hp = analyze_primer_structures(primer, tail_len=tail_len,
                                       exposed_ss=exposed_ss, temperature=STRUCT_TEMP_C,
                                       dg_3prime=_p.qc_hairpin_dg_3prime,
                                       dg_internal=_p.qc_hairpin_dg_internal,
                                       dg_absolute=_p.qc_dg_absolute)
        hp_full, hp_anneal, hp_exposed = hp['full'], hp['anneal'], hp['exposed']
        has_hairpin = hp_full['has_strong_structure']

        # 심각도 2단 tier 표기(🔴 strong ≤−9 / 🟡 caution ≤타입컷). raw ΔG도 함께 노출.
        _tier = {'strong': "🔴 strong", 'caution': "🟡 caution"}

        # anneal(spacer/증폭 결합부) 자체 hairpin. strong(ΔG ≤ −9)만 차단(red) — tier 체계와 일치.
        #   caution(≤ 타입컷 −2/−3, 그러나 > −9)은 경고(yellow)로 둔다.
        #   (구 로직은 has_strong_structure(≤−2/−3)로 violation을 걸어 caution까지 red로 찍혔음.)
        if hp_anneal['has_strong_structure']:
            kind = "3'-end" if hp_anneal['involves_3prime'] else "internal"
            _hp_msg = (f"Hairpin(anneal, {kind}, {_tier.get(hp_anneal['severity'], '')}) "
                       f"존재 (ΔG = {hp_anneal['delta_g']:.2f} kcal/mol)")
            if hp_anneal['severity'] == 'strong':
                violations.append(_hp_msg)
            else:
                warnings.append(_hp_msg)

        # (경고) ② 조립 시 노출되는 말단 ss(caller가 잘라 넘긴 exposed_ss)의 hairpin — 직접 측정.
        #   Gibson junction homology 등은 vector/method 의존이라 spacer로 못 고침 → 차단이 아니라
        #   vector/method 조합 재검토 신호로 명시(조용히 넘기지 않는다). ≤−9(strong)면 최상위 위험경고.
        if hp_exposed and hp_exposed['has_strong_structure']:
            kind = "3'-end" if hp_exposed['involves_3prime'] else "internal"
            structural_notes.append(
                f"Hairpin(노출 말단 ss, {kind}, {_tier.get(hp_exposed['severity'], '')}, "
                f"ΔG = {hp_exposed['delta_g']:.2f} kcal/mol) — 조립 시 junction anneal 저해. "
                f"spacer로 해결 불가, vector/method 조합 재검토 필요"
            )

        # (경고) anneal/노출ss 단독으로 설명 안 되는 전체-프라이머 구조(경계 걸침 등) = ① PCR 단계 신호.
        exposed_strong = bool(hp_exposed and hp_exposed['has_strong_structure'])
        if has_hairpin and not hp_anneal['has_strong_structure'] and not exposed_strong:
            structural_notes.append(
                f"Hairpin(전체 프라이머, {_tier.get(hp_full['severity'], '')}, "
                f"ΔG = {hp_full['delta_g']:.2f} kcal/mol) — "
                f"PCR 증폭 단계 구조(anneal/노출ss 단독으로는 설명 안 됨, 경계 걸침 등)"
            )

        # 7. 다른 프라이머와 dimer 체크
        if check_dimer_with:
            dimer_result = ValidatePrimer.check_primer_dimer(
                primer, check_dimer_with, threshold_dg=_p.qc_dimer_dg)
            if dimer_result['has_dimer'] and dimer_result['warning']:
                if warn_only_dimer:
                    structural_notes.append(f"프라이머 간 dimer: {dimer_result['warning']}")
                else:
                    violations.append(f"프라이머 간 dimer: {dimer_result['warning']}")

        # 전체 검증 통과 여부 + 3-tier 요약(🔴 red=차단 / 🟡 yellow=경고만 / 🟢 green)
        passed = len(violations) == 0
        tier = 'red' if violations else ('yellow' if warnings else 'green')

        return {
            'passed': passed,
            'tier': tier,
            # tier·화면 경고에는 안 쓰지만 버리지도 않는 구조 수치(dimer/hairpin ΔG)
            'structural_notes': structural_notes,
            'length': length,
            'tm': tm,
            'gc_content': gc_content,
            'has_gc_clamp': has_gc_clamp,
            'has_self_dimer': has_self_dimer,
            'has_hairpin': has_hairpin,
            'has_hairpin_3prime': hp_full['involves_3prime'],
            # 세 관점 hairpin ΔG(kcal/mol) 직접 측정값 — 나란히 조회용. exposed는 미제공 시 None.
            'hairpin_dg_full': hp_full['delta_g'],
            'hairpin_dg_anneal': hp_anneal['delta_g'],
            'hairpin_dg_exposed': hp_exposed['delta_g'] if hp_exposed else None,
            'violations': violations,
            'warnings': warnings
        }

# 함수 노출 (테스트 코드 호환성용)
design_basic_primers = PrimerDesign.design_basic_primers
check_primer_dimer = ValidatePrimer.check_primer_dimer
check_gc_clamp = ValidatePrimer.check_gc_clamp
check_primer_specificity = ValidatePrimer.check_primer_specificity
validate_primer = ValidatePrimer.validate_primer

# 모듈 메타데이터
__all__ = [
    'design_basic_primers',
    'check_primer_dimer',
    'check_gc_clamp',
    'check_primer_specificity',
    'validate_primer',
    'PrimerDesign',
    'ValidatePrimer'
]
