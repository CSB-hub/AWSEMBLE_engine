"""
DNA 2차 구조 분석

이 모듈은 DNA 서열의 2차 구조(hairpin)를 분석합니다.

핵심 개선사항:
    ❌ 기존 방식: 전체 DNA → RNA 변환 → 2차 구조 계산 (생물학적으로 부적절)
    ✅ 단일가닥 영역만 hairpin을 계산 (PCR degradation 영역 / 프라이머 등)

버그3 재설계(2026-07-02): seqfold → primer3 통일.
    구버전은 `seqfold.fold()` 반환("구조 요소 리스트")을 "대안 구조 리스트"로 오독해
    `structures[0].e`(첫 base-pair stack 하나)만 읽어 hairpin 탐지가 사실상 죽어 있었다.
    primer3-py는 이미 dimer 계산에 쓰는 의존성이고 ΔG·Tm·염/온도 조건을 모두 제공하므로
    hairpin도 `primer3.calc_hairpin`으로 통일한다(단일 소스, 새 의존성 없음).

생물학적 근거:
    - PCR 과정: 95°C 변성 → 단일가닥 DNA 생성
    - 50-60°C Annealing: 프라이머 결합 시 2차 구조 형성 가능
    - 이중가닥 DNA는 안정적이므로 2차 구조 형성 안 함

라이브러리:
    - primer3-py: hairpin/dimer 열역학 계산 (Nearest-Neighbor model, 업계 표준)

작성일: 2025-10-30
수정일: 2026-07-02 (버그3: seqfold 오독 경로 제거 → primer3 ΔG 게이트로 재작성)
"""

from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq

import primer3


# ─────────────────────────────────────────────────────────────────────────────
# 2차구조 ΔG 게이트 설정 (⚠️ wet-lab 확정 전 커뮤니티 관례 default)
#   - primer3 `.dg`는 cal/mol 반환 → 여기서는 /1000 후 kcal/mol로 통일해 비교한다
#     (기존 check_primer_dimer가 이미 /1000 kcal 방식이라 단일 방식으로 맞춤).
#   - 컷오프보다 "더 음수"(= 더 안정한 구조)면 문제로 판정한다.
#   - 정식 검증값은 primer3 `PRIMER_MAX_HAIRPIN_TH`(Tm 기반)뿐이고, 아래 ΔG 숫자는
#     커뮤니티 관례(premierbiosoft/oligopool 등)이지 단일 권위 표준이 아니다.
#     → 실제 PCR·조립 성공률에 직결되므로 wet-lab 합의로 확정 대상(정책 판단 영역).
# ─────────────────────────────────────────────────────────────────────────────
HAIRPIN_DG_3PRIME_KCAL = -2.0      # 3'말단 관여 hairpin: 가장 엄격 (3'이 sequester되면 mispriming/증폭 실패)
HAIRPIN_DG_INTERNAL_KCAL = -3.0    # 내부 hairpin: 상대적 관용
DIMER_DG_CUTOFF_KCAL = -5.0        # self/cross dimer (농도 의존)
DG_ABSOLUTE_STRICT_KCAL = -9.0     # 위치·형태 불문 강력 구조 배제(절대 컷)
HAIRPIN_3PRIME_WINDOW = 3          # 3'말단에서 몇 nt까지 "3'-end 관여"로 볼지 (구조 주석 끝 window)

# calc_* 온도/염 조건.
#   ⚠️ 관례 ΔG 컷오프는 대개 37°C 기준이라 분석 온도도 37°C로 맞춘다(컷오프와 정합).
#   실제 반응 온도(Ta)로 재평가할지는 정책 결정 대상(위 wet-lab 확정과 함께).
STRUCT_TEMP_C = 37.0
STRUCT_MV_CONC = 50.0    # monovalent(Na+) mM
STRUCT_DV_CONC = 0.0     # divalent(Mg2+) mM (기존 check_primer_dimer와 일치)
STRUCT_DNTP_CONC = 0.0   # dNTP mM
STRUCT_DNA_CONC = 250.0  # DNA nM (기존 코드와 일치)

# ─────────────────────────────────────────────────────────────────────────────
# template target-site accessibility (v1, warn-only) — ViennaRNA 미페어확률 기반
#   primer ΔG(self-hairpin/dimer)와 **별개 축**: primer가 붙을 template 자리가 이미 접혀
#   접근 불가한지 본다. 입력이 template(못 바꾸는 축)이라 **fail 아님 → 경고 + 조건 추천만**.
#   지표 = 결합부위 각 base의 미페어확률(unpaired probability, 길이 독립 0~1) 평균. 낮을수록 접힘.
#   ⚠️ 아래 임계·온도는 관례 default → wet-lab 확정 대상(hairpin 컷과 동일하게 정책 영역).
#   ※ 바이오파츠 어셈블은 결합 위치가 **말단 고정**(부품을 양 끝에서 증폭)이라 "위치 선택" 여지가 없다.
#     조정 가능한 건 GC%/Tm 기반 결합부위 **길이**뿐(_fit_amplicon/primer3가 담당). 따라서 accessibility는
#     "그 고정된 말단 binding site가 template 내에서 접힐 가능성"만 본다(위치 스코어링 아님).
# ─────────────────────────────────────────────────────────────────────────────
ANNEAL_TEMP_C = 55.0                    # 접근성 평가 온도(Ta 근사). hairpin의 STRUCT_TEMP_C(37)와 별개.
ACCESSIBILITY_FLANK = 100               # 결합부위 안쪽 문맥 nt(RNAplfold류 국소창 — 원거리 pairing 포착).
ACCESSIBILITY_MAX_WINDOW = 400          # 접힘 창 최대 길이(pf 비용 상한; 초과 시 결합부위 중심으로 clamp).
# 결합부위 평균 미페어확률 기준 (2026-10-20 실측 재조정).
#
# 구버전 값(0.50 / 0.25)의 문제: 문헌 파츠 102개 결합부위의 **중앙값이 0.53**이라
# 0.50은 중앙값 바로 옆이었다. "비정상"이 아니라 "보통"에 경고를 달아 47%가 울렸다.
#
#   분포: 최소 0.17 │ 25% 0.43 │ 중앙 0.53 │ 75% 0.70 │ 최대 0.76
#   기준별 경고율: 0.50→47% / 0.40→15% / 0.35→12% / 0.25→12% / 0.20→6%
#
# 0.25~0.40 구간이 비어 있어(0.25와 0.35의 개수가 같다) 분포에 자연스러운 틈이 있다.
# 그 틈 위를 주의선으로 잡으면 진짜 막힌 자리만 걸린다.
# 2차구조 그림을 함께 띄우므로 세부 판단은 그림이 맡고, 경고는 신호만 보낸다.
ACCESSIBILITY_MIN_UNPAIRED = 0.40       # 미만이면 접근성 낮음(🟡 주의) — 15%
ACCESSIBILITY_MIN_UNPAIRED_STRONG = 0.20  # 미만이면 강한 경고(🔴 위험) — 6%


def _structure_annotation(ascii_lines: Optional[List[str]], seq_len: int) -> str:
    """primer3 `ascii_structure_lines`에서 서열-정렬 구조 주석을 추출.

    calc_hairpin(output_structure=True)는 두 줄을 준다:
        'SEQ\\t//////////----\\\\...'  (구조 주석: '/'=5'페어, '\\'=3'페어, '-'=미페어)
        'STR\\tGCGC...'                (서열)
    주석은 서열과 위치가 1:1로 정렬된다. 구조 심볼('/','\\','-')로만 이뤄지고
    길이가 서열과 같은 토큰을 주석으로 채택한다(없으면 빈 문자열).
    """
    if not ascii_lines:
        return ''
    struct_chars = set('/\\- ')
    for line in ascii_lines:
        token = line.split('\t')[-1]
        if token and set(token) <= struct_chars and len(token) == seq_len:
            return token
    return ''


def _involves_3prime_end(annotation: str, window: int) -> bool:
    """구조 주석의 3'말단(끝) window nt 안에 페어('/' 또는 '\\')가 있으면 True."""
    if not annotation:
        return False
    tail = annotation[-window:] if window > 0 else annotation
    return any(c in '/\\' for c in tail)


def get_pcr_degradation_region(
    template: str,
    primer_length: int,
    extended_length: int = 50
) -> Tuple[str, str]:
    """
    PCR에서 단일가닥이 되는 영역 추출

    PCR 과정에서:
    1. 95°C 변성 → 전체 DNA가 단일가닥
    2. 50-60°C Annealing → 프라이머가 결합하는 영역은 이중가닥
    3. 프라이머 결합 영역 주변은 여전히 단일가닥

    이 함수는 프라이머 결합 영역 + 주변 영역을 추출합니다.

    Args:
        template: Template DNA 서열
        primer_length: 프라이머 길이
        extended_length: 프라이머 결합 영역 이외 추가로 포함할 길이 (bp)

    Returns:
        (forward_region, reverse_region)
        - forward_region: Forward 프라이머 결합 영역 + extended
        - reverse_region: Reverse 프라이머 결합 영역 + extended

    Example:
        >>> template = "ATCG" * 100  # 400 bp
        >>> fwd, rev = get_pcr_degradation_region(template, primer_length=20)
        >>> len(fwd) == 20 + 50  # primer + extended
        True
        >>> len(rev) == 20 + 50
        True
    """
    if not template:
        raise ValueError("Template 서열이 비어있습니다")

    if primer_length + extended_length > len(template):
        raise ValueError(
            f"Template 길이({len(template)} bp)가 "
            f"primer_length + extended_length({primer_length + extended_length} bp)보다 작습니다"
        )

    # Forward 영역: 5' 끝에서 primer_length + extended_length
    forward_region = template[:primer_length + extended_length]

    # Reverse 영역: 3' 끝에서 primer_length + extended_length
    reverse_region = template[-(primer_length + extended_length):]

    return forward_region, reverse_region


def analyze_single_strand_structure(
    sequence: str,
    temperature: float = STRUCT_TEMP_C,
    three_prime_window: int = HAIRPIN_3PRIME_WINDOW,
    dg_3prime: Optional[float] = None,
    dg_internal: Optional[float] = None,
    dg_absolute: Optional[float] = None,
) -> Dict:
    """
    단일가닥 DNA의 hairpin 2차 구조 분석 (primer3 `calc_hairpin` 기반)

    버그3 재작성: 구버전 seqfold 오독 경로를 제거하고 primer3 열역학으로 통일한다.
    `output_structure=True`로 얻은 구조 주석을 파싱해 hairpin이 3'말단에 관여하는지
    구분하고, 3'-end 관여 여부에 따라 다른 ΔG 컷오프를 적용한다.

    Args:
        sequence: 단일가닥 DNA 서열
        temperature: 분석 온도 (°C, 기본값 STRUCT_TEMP_C=37 — 관례 컷오프와 정합)
        three_prime_window: 3'말단에서 몇 nt까지 "3'-end 관여"로 볼지 (기본 HAIRPIN_3PRIME_WINDOW)

    Returns:
        {
            'delta_g': float,              # hairpin ΔG (kcal/mol, 음수 = 안정)
            'tm': float,                   # hairpin Tm (°C)
            'structure_found': bool,       # primer3가 구조를 찾았는지
            'involves_3prime': bool,       # hairpin이 3'말단(마지막 window nt)에 관여하는지
            'structure': str,              # primer3 구조 주석('/','\\','-'), 없으면 ''
            'severity': str,               # 심각도 tier: 'strong'(≤−9) / 'caution'(≤타입컷) / 'ok'
            'hairpins': List[Dict],        # (호환) 검출 hairpin 요약
            'has_strong_structure': bool   # ΔG ≤ 타입컷(3'-end −2 / internal −3)
        }

    Raises:
        ValueError: 빈 서열

    Example:
        >>> seq = "GCGCGAGCTCTTTTGAGCTCGCGC"  # 강한 3'말단 관여 hairpin
        >>> result = analyze_single_strand_structure(seq)
        >>> result['has_strong_structure']
        True
        >>> result['delta_g'] < -2  # 안정한 구조
        True
    """
    if not sequence:
        raise ValueError("서열이 비어있습니다")

    sequence = sequence.upper()

    # primer3는 'N'을 처리하지 못하므로 'A'로 치환 (에너지 계산용)
    # TODO: N을 A로 치환하는게 맞는지 확인필요. TGC와 A의 차이가 많이 나기 때문
    clean_sequence = sequence.replace('N', 'A')

    # ── 하이브리드 경로 (H1) ──────────────────────────────────────────────
    # primer3 `calc_hairpin`은 60bp를 넘으면 RuntimeError로 거부한다. 그 길이는
    # **합성 상한이 아니라 계산 상한**이라, 63nt 프라이머는 주문·실험이 되는데 우리만
    # 품질을 못 보는 상태가 된다. 그 구간은 ViennaRNA(DNA 파라미터)로 계산한다.
    #
    # ⚠️ **경고 전용**이다. ViennaRNA MFE와 primer3 hairpin ΔG는 스케일이 달라서
    #    기존 컷오프(−2/−3/−5/−9)를 그대로 적용할 수 없다. 그 컷오프들은 지금도
    #    wet-lab 확정 대기 상태라 재캘리브레이션 기준점이 없다. 그래서 이 경로는
    #    `has_strong_structure=False` / `severity='ok'`로 두어 **어떤 게이트도 움직이지
    #    않게** 하고, 값과 구조만 리포트로 넘긴다.
    from params import THERMO_MAX_LEN

    if len(clean_sequence) > THERMO_MAX_LEN:
        return _analyze_long_strand(sequence, clean_sequence, temperature,
                                    three_prime_window)

    result = primer3.calc_hairpin(
        clean_sequence,
        mv_conc=STRUCT_MV_CONC,
        dv_conc=STRUCT_DV_CONC,
        dntp_conc=STRUCT_DNTP_CONC,
        dna_conc=STRUCT_DNA_CONC,
        temp_c=temperature,
        output_structure=True,
    )

    if not result.structure_found:
        # 구조 없음 → 기본값
        return {
            'delta_g': 0.0,
            'tm': 0.0,
            'structure_found': False,
            'involves_3prime': False,
            'structure': '',
            'severity': 'ok',
            'hairpins': [],
            'has_strong_structure': False,
            'engine': 'primer3',
            'warn_only': False,
        }

    delta_g = result.dg / 1000.0  # cal/mol → kcal/mol (단위 통일)
    annotation = _structure_annotation(result.ascii_structure_lines, len(clean_sequence))
    involves_3prime = _involves_3prime_end(annotation, three_prime_window)

    # 3'말단 관여면 엄격한 컷(−2), 내부는 관용(−3).
    # ⚠️ 컷오프는 호출자가 주입할 수 있다(None이면 모듈 기본값). 주입 경로는
    #    `validate_primer` → `params.qc_hairpin_dg_*` 이며, **유저가 GUI에서 만지는 값이 아니다**
    #    (판정 축은 설계 목표와 분리 — params.py의 qc_* 설명 참조).
    _dg3 = HAIRPIN_DG_3PRIME_KCAL if dg_3prime is None else dg_3prime
    _dgi = HAIRPIN_DG_INTERNAL_KCAL if dg_internal is None else dg_internal
    _dga = DG_ABSOLUTE_STRICT_KCAL if dg_absolute is None else dg_absolute
    cutoff = _dg3 if involves_3prime else _dgi
    has_strong_structure = delta_g <= cutoff

    # 심각도 2단 tier (경고 등급): strong(≤ 절대컷 −9) / caution(≤ 타입컷) / ok(구조 약함).
    #   ⚠️ 구 로직은 has_strong에 '≤ −9'를 OR로 넣었으나 −9 < −2·−3이라 항상 타입컷에
    #   먼저 걸려 dead였다. 이제 −9는 kill이 아니라 '경고 최상위 tier'로만 쓴다(2축 정책).
    if delta_g <= _dga:
        severity = 'strong'
    elif has_strong_structure:
        severity = 'caution'
    else:
        severity = 'ok'

    # (호환) hairpin 요약 1건 — 위치 세부는 구조 주석(structure)에 담긴다.
    hairpins = [{
        'delta_g': delta_g,
        'tm': result.tm,
        'involves_3prime': involves_3prime,
        'sequence': sequence,
    }]

    return {
        'delta_g': delta_g,
        'tm': result.tm,
        'structure_found': True,
        'involves_3prime': involves_3prime,
        'structure': annotation,
        'severity': severity,
        'hairpins': hairpins,
        'has_strong_structure': has_strong_structure,
        'engine': 'primer3',
        'warn_only': False,
    }


def _vienna_dna_md(temperature: float):
    """ViennaRNA 모듈 + DNA 파라미터가 적용된 model detail을 준다.

    `params_load_DNA_Mathews2004()`는 **전역 설정**이라 폴딩 직전에 매번 세팅한다
    (`analyze_target_accessibility`도 같은 방식이다). 그 설정 코드를 한 곳에 모아두는 건
    "같은 일을 하는 코드를 두 벌 만들지 않는다"는 원칙 때문이다 — hairpin(단일가닥)과
    dimer(이중가닥) 두 경로가 같은 준비를 필요로 한다.

    Args:
        temperature: 폴딩 온도(°C)

    Returns:
        `(RNA 모듈, model_details)` 또는 ViennaRNA 미설치 시 None.
    """
    try:
        import RNA
    except ImportError:
        return None
    RNA.params_load_DNA_Mathews2004()   # DNA NN 파라미터(SantaLucia류) — 전역 설정
    md = RNA.md()
    md.temperature = temperature
    return RNA, md


def analyze_long_duplex(
    primer1: str,
    primer2: str,
    temperature: float = STRUCT_TEMP_C,
) -> Optional[Dict]:
    """primer3 상한을 넘는 두 서열의 duplex(dimer) — ViennaRNA cofold (H1, **경고 전용**).

    `primer3.calc_heterodimer`는 한쪽이라도 60bp를 넘으면 거부한다. 조립용 프라이머는
    `[5' tail][anneal]` 구조라 arm이 길어지면 이 구간에 쉽게 들어가고, 그러면 dimer 축이
    통째로 눈을 감는다. 그 구간을 ViennaRNA cofold로 메운다.

    Args:
        primer1: 서열 1
        primer2: 서열 2
        temperature: 분석 온도(°C)

    Returns:
        `{'has_dimer', 'delta_g', 'tm', 'warning', 'engine', 'warn_only'}`.
        ViennaRNA가 없으면 None(호출자가 기존 폴백을 쓰도록).

    Note:
        ΔG 스케일이 primer3와 다르다 — `has_dimer`는 **항상 False**로 두어 게이트를
        움직이지 않는다. 값은 경고 문구로만 나간다.
    """
    vienna = _vienna_dna_md(temperature)
    if vienna is None:
        return None
    RNA, md = vienna
    a = (primer1 or "").upper().replace("N", "A")
    b = (primer2 or "").upper().replace("N", "A")
    if not a or not b:
        return None

    # 두 가닥은 '&'로 이어 넘긴다. 리스트로 넘기면 첫 서열만 접혀 **단일가닥 결과**가
    # 조용히 나온다(실측: 상보 70mer 쌍이 ΔG 0으로 보고됨).
    fc = RNA.fold_compound(f"{a}&{b}", md)
    _structure, mfe = fc.mfe_dimer()
    delta_g = float(mfe)

    return {
        'has_dimer': False,          # 게이트를 움직이지 않는다(스케일 상이)
        'delta_g': delta_g,
        'tm': 0.0,                   # cofold MFE는 Tm을 주지 않는다
        'warning': (
            f"dimer 참고값(ViennaRNA cofold, {len(a)}/{len(b)}nt): "
            f"ΔG {delta_g:.2f} kcal/mol — primer3 계산 상한을 넘는 길이라 "
            f"**스케일이 달라 기존 컷오프와 비교하면 안 됩니다(경고 전용)**."
        ),
        'engine': 'viennarna',
        'warn_only': True,
    }


def _analyze_long_strand(
    sequence: str,
    clean_sequence: str,
    temperature: float,
    three_prime_window: int,
) -> Dict:
    """primer3 상한을 넘는 단일가닥의 2차구조 — ViennaRNA MFE (H1, **경고 전용**).

    `analyze_single_strand_structure`와 **같은 키**를 돌려주되 게이트를 움직이지 않는다.
    호출자는 `engine`/`warn_only`로 어느 엔진이 계산했는지 구분할 수 있다.

    DNA 파라미터(`params_load_DNA_Mathews2004`)를 쓰는 것은 `analyze_target_accessibility`와
    같은 방식이다 — 전역 설정이라 매 호출 직전에 세팅한다.

    Args:
        sequence: 원본 서열(리포트 표기용)
        clean_sequence: 'N'을 'A'로 치환한 계산용 서열
        temperature: 분석 온도(°C)
        three_prime_window: 3'말단 관여 판정 창(nt)

    Returns:
        `analyze_single_strand_structure`와 동일한 스키마.
        `severity='na'`는 ViennaRNA가 없어 계산을 건너뛴 경우다.

    Note:
        ΔG 스케일이 primer3와 다르다 — **기존 컷오프와 직접 비교하면 안 된다.**
    """
    from params import THERMO_MAX_LEN

    vienna = _vienna_dna_md(temperature)
    if vienna is None:
        # ViennaRNA 미설치 → 건너뛴다(파이프라인 무중단). accessibility와 같은 정책.
        return {
            'delta_g': 0.0,
            'tm': 0.0,
            'structure_found': False,
            'involves_3prime': False,
            'structure': '',
            'severity': 'na',
            'hairpins': [],
            'has_strong_structure': False,
            'engine': 'none',
            'warn_only': True,
            'note': (f"{len(clean_sequence)}bp는 primer3 계산 상한을 넘고 ViennaRNA도 "
                     f"없어 2차구조를 계산하지 못했습니다."),
        }

    RNA, md = vienna
    fc = RNA.fold_compound(clean_sequence, md)
    dot_bracket, mfe = fc.mfe()
    delta_g = float(mfe)

    # ViennaRNA는 dot-bracket을 주므로 3'말단 관여를 직접 읽을 수 있다
    # (primer3의 '/','\\','-' 주석과 형식이 달라 `_involves_3prime_end`를 쓸 수 없다).
    tail = dot_bracket[-three_prime_window:] if three_prime_window else ''
    involves_3prime = any(c in '()' for c in tail)
    paired = sum(1 for c in dot_bracket if c in '()')

    return {
        'delta_g': delta_g,
        'tm': 0.0,                      # ViennaRNA MFE는 Tm을 주지 않는다
        'structure_found': paired > 0,
        'involves_3prime': involves_3prime,
        'structure': dot_bracket,
        # ⚠️ 게이트를 움직이지 않는다 — 스케일이 달라 기존 컷오프를 적용할 수 없다.
        'severity': 'ok',
        'hairpins': [{
            'delta_g': delta_g,
            'tm': 0.0,
            'involves_3prime': involves_3prime,
            'sequence': sequence,
        }] if paired else [],
        'has_strong_structure': False,
        'engine': 'viennarna',
        'warn_only': True,
        'note': (f"{len(clean_sequence)}bp — primer3 계산 상한({THERMO_MAX_LEN}bp)을 "
                 f"넘어 ViennaRNA(DNA 파라미터)로 계산했습니다. MFE ΔG {delta_g:.2f} kcal/mol, "
                 f"페어 염기 {paired}개"
                 f"{', 3′말단 관여' if involves_3prime else ''}. "
                 f"**스케일이 달라 기존 ΔG 컷오프와 직접 비교하면 안 됩니다 — 경고 전용입니다.**"),
    }


def analyze_primer_structures(
    primer: str,
    tail_len: int = 0,
    exposed_ss: Optional[str] = None,
    temperature: float = STRUCT_TEMP_C,
    dg_3prime: Optional[float] = None,
    dg_internal: Optional[float] = None,
    dg_absolute: Optional[float] = None,
) -> Dict:
    """프라이머/조각의 hairpin을 물리적 단계가 다른 **세 관점**으로 각각 직접 계산.

    primer3만 사용해 각 구간의 hairpin ΔG를 나란히 뽑는다(추론이 아니라 측정).

    - full   : ① PCR 증폭 단계 — tail+anneal 전체가 자유 ssDNA일 때의 self-hairpin.
    - anneal : tail을 뺀 3' 결합부(spacer/증폭 영역). 다른 spacer로 교체 가능한 구간.
               tail_len=0이면 full과 동일.
    - exposed: ② 조립 단계 — T5 등이 노출시켜 partner와 anneal해야 하는 말단 ss 구간.
               caller가 method를 알고 잘라 넘긴다(Gibson=homology 구간 / GG=4nt sticky→None).
               (서열과 그 역상보는 self-hairpin ΔG가 동일하므로 노출부 서열을 그대로 넣으면 된다.)
               None이면 계산하지 않는다.

    Args:
        primer: 프라이머(또는 조각) 서열
        tail_len: 5' 비-anneal tail 길이 (anneal = primer[tail_len:])
        exposed_ss: 조립 시 외부 노출되는 말단 ss 서열(없으면 None → 계산 생략)
        temperature: 분석 온도 (기본 STRUCT_TEMP_C)

    Returns:
        {
            'full':    analyze_single_strand_structure(primer),
            'anneal':  analyze_single_strand_structure(anneal),
            'exposed': analyze_single_strand_structure(exposed_ss) 또는 None,
        }
    """
    _cut = dict(dg_3prime=dg_3prime, dg_internal=dg_internal, dg_absolute=dg_absolute)
    full = analyze_single_strand_structure(primer, temperature=temperature, **_cut)
    anneal = (
        analyze_single_strand_structure(primer[tail_len:], temperature=temperature, **_cut)
        if tail_len else full
    )
    exposed = (
        analyze_single_strand_structure(exposed_ss, temperature=temperature, **_cut)
        if exposed_ss else None
    )
    return {'full': full, 'anneal': anneal, 'exposed': exposed}


def analyze_target_accessibility(
    template: str,
    start: int,
    end: int,
    temp: float = ANNEAL_TEMP_C,
    flank: int = ACCESSIBILITY_FLANK,
) -> Dict:
    """primer가 붙을 template 자리([start,end))의 접근성 평가 (v1, 경고용).

    ViennaRNA로 결합부위+flank 문맥 창을 접어(DNA 파라미터, Ta) 결합부위 각 base의
    **미페어확률(unpaired probability)** 평균을 구한다. 낮을수록 자리가 접혀 primer 접근 저해.
    template은 못 바꾸는 축 → **fail 아님**(경고 + 첨가제/온도 조건 추천만).

    ViennaRNA 미설치(이미지 재빌드 전)면 조용히 skip(accessible=True, severity='na').

    Args:
        template: 전체 template DNA 서열
        start/end: 결합부위 좌표 (0-based, end exclusive)
        temp: 평가 온도 (기본 ANNEAL_TEMP_C, Ta 근사)
        flank: 결합부위 양쪽 문맥 nt

    Returns:
        {
            'mean_unpaired': float,   # 결합부위 평균 미페어확률(0~1). 높을수록 접근가능
            'accessible': bool,       # 임계 이상이면 True
            'severity': str,          # 'ok'/'caution'/'strong'/'na'(미설치)
            'temp': float,
            'warning': str            # 비차단 경고(접근성 낮을 때만)
        }
    """
    ok = {'mean_unpaired': 1.0, 'accessible': True, 'severity': 'ok', 'temp': temp, 'warning': ''}
    if not template or start >= end or start < 0:
        return ok
    try:
        import RNA
    except ImportError:
        # ViennaRNA 미설치(이미지 재빌드 전) → accessibility skip(파이프라인 무중단).
        return {**ok, 'severity': 'na'}

    # 결합부위 ± flank 창(국소 folding). 단 pf 비용 상한(MAX_WINDOW) 초과 시 결합부위 중심으로 clamp.
    eff_flank = flank
    if (end - start) + 2 * flank > ACCESSIBILITY_MAX_WINDOW:
        eff_flank = max(0, (ACCESSIBILITY_MAX_WINDOW - (end - start)) // 2)
    ws = max(0, start - eff_flank)
    we = min(len(template), end + eff_flank)
    window = template[ws:we].upper().replace('N', 'A')
    bind_lo, bind_hi = start - ws, end - ws   # 창 내 결합부위(0-based)
    if len(window) < 4 or bind_hi <= bind_lo:
        return ok

    RNA.params_load_DNA_Mathews2004()   # DNA NN 파라미터(SantaLucia류)
    md = RNA.md()
    md.temperature = temp
    fc = RNA.fold_compound(window, md)
    fc.pf()                             # partition function
    bpp = fc.bpp()                      # (n+1)x(n+1), 1-indexed 상삼각
    n = len(window)
    vals = []
    for i in range(bind_lo, bind_hi):
        pi = i + 1
        paired = sum(bpp[min(pi, j)][max(pi, j)] for j in range(1, n + 1) if j != pi)
        vals.append(max(0.0, 1.0 - paired))
    mean_unpaired = sum(vals) / len(vals) if vals else 1.0

    if mean_unpaired < ACCESSIBILITY_MIN_UNPAIRED_STRONG:
        severity, tier = 'strong', "🔴 strong"
    elif mean_unpaired < ACCESSIBILITY_MIN_UNPAIRED:
        severity, tier = 'caution', "🟡 caution"
    else:
        severity, tier = 'ok', ''
    accessible = severity == 'ok'
    warning = '' if accessible else (
        f"template 결합부위 접근성 낮음({tier}, 평균 미페어확률 {mean_unpaired:.2f} @ {temp}°C) — "
        f"자리가 접혀 프라이머 접근 저해 가능. template 고정이라 재설계 불가: "
        f"첨가제(DMSO/betaine)·어닐링 온도 조정 권장"
    )
    return {
        'mean_unpaired': mean_unpaired,
        'accessible': accessible,
        'severity': severity,
        'temp': temp,
        'warning': warning,
    }


def analyze_pcr_structure(
    template: str,
    forward_primer: str,
    reverse_primer: str,
    annealing_temp: float = 55.0
) -> Dict:
    """
    PCR degradation 영역의 2차 구조 분석

    ⭐ 핵심 개선!

    기존 방식 (잘못됨):
        전체 DNA → RNA 변환 → 2차 구조 계산

    새로운 방식 (올바름):
        1. Forward primer 결합 영역 추출
        2. Reverse primer 결합 영역 추출
        3. 각 영역의 2차 구조를 개별 계산
        4. 프라이머 자체의 2차 구조도 계산

    생물학적 근거:
        PCR 과정에서 95°C 변성 → 단일가닥
        50-60°C 프라이머 결합 시 2차 구조 가능
        이 단계에서만 2차 구조 체크가 의미 있음!

    Args:
        template: Template DNA
        forward_primer: Forward 프라이머
        reverse_primer: Reverse 프라이머
        annealing_temp: Annealing 온도 (°C)

    Returns:
        {
            'fwd_binding_region': str,      # Forward 결합 영역 서열
            'rev_binding_region': str,      # Reverse 결합 영역 서열
            'fwd_binding_dg': float,        # Forward 결합 영역 ΔG
            'rev_binding_dg': float,        # Reverse 결합 영역 ΔG
            'fwd_primer_dg': float,         # Forward 프라이머 ΔG
            'rev_primer_dg': float,         # Reverse 프라이머 ΔG
            'has_strong_structure': bool,   # 강한 구조 존재 여부
            'warnings': List[str]           # 경고 메시지
        }

    Raises:
        ValueError: 빈 서열 또는 프라이머

    Example:
        >>> template = "ATCG" * 100
        >>> fwd = "ATCGATCG"
        >>> rev = "CGATATAT"
        >>> result = analyze_pcr_structure(template, fwd, rev)
        >>> 'fwd_binding_dg' in result
        True
        >>> 'warnings' in result
        True
    """
    if not template:
        raise ValueError("Template 서열이 비어있습니다")
    if not forward_primer:
        raise ValueError("Forward 프라이머가 비어있습니다")
    if not reverse_primer:
        raise ValueError("Reverse 프라이머가 비어있습니다")

    warnings = []

    # PCR degradation 영역 추출
    fwd_region, rev_region = get_pcr_degradation_region(
        template,
        primer_length=max(len(forward_primer), len(reverse_primer)),
        extended_length=50
    )

    # Forward 결합 영역 2차 구조 분석
    fwd_binding_analysis = analyze_single_strand_structure(
        fwd_region,
        temperature=annealing_temp
    )

    # Reverse 결합 영역 2차 구조 분석
    rev_binding_analysis = analyze_single_strand_structure(
        rev_region,
        temperature=annealing_temp
    )

    # Forward 프라이머 자체의 2차 구조 분석
    fwd_primer_analysis = analyze_single_strand_structure(
        forward_primer,
        temperature=annealing_temp
    )

    # Reverse 프라이머 자체의 2차 구조 분석
    rev_primer_analysis = analyze_single_strand_structure(
        reverse_primer,
        temperature=annealing_temp
    )

    # 경고 메시지 생성
    if fwd_binding_analysis['has_strong_structure']:
        warnings.append(
            f"Forward 결합 영역에 강한 2차 구조 존재 (ΔG = {fwd_binding_analysis['delta_g']:.2f} kcal/mol)"
        )

    if rev_binding_analysis['has_strong_structure']:
        warnings.append(
            f"Reverse 결합 영역에 강한 2차 구조 존재 (ΔG = {rev_binding_analysis['delta_g']:.2f} kcal/mol)"
        )

    if fwd_primer_analysis['has_strong_structure']:
        warnings.append(
            f"Forward 프라이머에 강한 hairpin 존재 (ΔG = {fwd_primer_analysis['delta_g']:.2f} kcal/mol)"
        )

    if rev_primer_analysis['has_strong_structure']:
        warnings.append(
            f"Reverse 프라이머에 강한 hairpin 존재 (ΔG = {rev_primer_analysis['delta_g']:.2f} kcal/mol)"
        )

    # 강한 구조 존재 여부
    has_strong_structure = (
        fwd_binding_analysis['has_strong_structure'] or
        rev_binding_analysis['has_strong_structure'] or
        fwd_primer_analysis['has_strong_structure'] or
        rev_primer_analysis['has_strong_structure']
    )

    return {
        'fwd_binding_region': fwd_region,
        'rev_binding_region': rev_region,
        'fwd_binding_dg': fwd_binding_analysis['delta_g'],
        'rev_binding_dg': rev_binding_analysis['delta_g'],
        'fwd_primer_dg': fwd_primer_analysis['delta_g'],
        'rev_primer_dg': rev_primer_analysis['delta_g'],
        'has_strong_structure': has_strong_structure,
        'warnings': warnings
    }


# 모듈 메타데이터
__all__ = [
    'get_pcr_degradation_region',
    'analyze_single_strand_structure',
    'analyze_primer_structures',
    'analyze_target_accessibility',
    'analyze_pcr_structure',
    # 설정 상수(외부에서 컷오프/온도/염 조회·조정용)
    'HAIRPIN_DG_3PRIME_KCAL',
    'HAIRPIN_DG_INTERNAL_KCAL',
    'DIMER_DG_CUTOFF_KCAL',
    'DG_ABSOLUTE_STRICT_KCAL',
    'HAIRPIN_3PRIME_WINDOW',
    'STRUCT_TEMP_C',
]
