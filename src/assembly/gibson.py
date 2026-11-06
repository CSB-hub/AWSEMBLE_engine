"""
Gibson Assembly 구현

Homology-based overlap을 사용한 DNA 조립 방법입니다.

핵심 기능:
    - Vector linearization (제한효소 자동 선택)
    - 프라이머 설계 (Basic 모듈 + overlap 추가)
    - Overlap 길이 최적화 (GC%, Tm 기반)
    - 프라이머 특이성 검증
    - Circular/Linear assembly

워크플로우:
    1. 입력 검증 (insert/vector circular/linear)
    2. Vector linearization (SmaI → EcoRV → PvuII 등)
    3. 프라이머 설계 (Basic 모듈)
    4. Overlap 추가 (vector 3'/5' 끝 서열)
    5. 프라이머 검증 (Tm, GC%, 특이성)
    6. Assembly 시뮬레이션 (PyDNA)

참고:
    - KRIBB Line 396-539
    - archive/docs/development/PLAN_v2.md #1 GibsonCloning

작성일: 2025-11-06
"""

from typing import Tuple, List, Dict, Optional
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import Restriction
from .base import AssemblyBase

#: 3-tier 품질 라벨의 나쁜 순서. fragment 전체 tier는 fwd/rev 중 **나쁜 쪽**을 따른다.
_TIER_ORDER = ("red", "unknown", "yellow", "green")


def _worst_tier(*tiers: str) -> str:
    """여러 tier 중 가장 나쁜 것. 모르는 라벨은 `unknown` 취급."""
    known = [t for t in tiers if t in _TIER_ORDER]
    if not known:
        return "unknown"
    return min(known, key=_TIER_ORDER.index)


class GibsonCloning(AssemblyBase):
    """Gibson Assembly (homology-based overlap) 구현

    특징:
        - 15-25bp homology overlap 사용
        - Linear 및 Circular DNA 조립 지원
        - Vector 자동 linearization
        - 프라이머 자동 설계 및 검증

    Attributes:
        overlap_length: Overlap 길이 (기본값: 20bp)

    작성일: 2025-11-06
    """

# TODO: 전반적으로 DBTL 버전과 general 버전에서의 사용에 대해서 생각해야 할 것 같음.
# 처음에는 일반적인 gibson 코드를짜고 DBTL 버전에 응용하는 식으로만 했는데 실제  workflow 코드를 보고 해당 코드의 응용인지
# 따로 개발을 해야하는지 판단.

    #: homology arm 길이 **범위**를 담은 `DesignParams` 필드 이름 (min, max).
    #: 하위 클래스가 다른 축을 쓰면 이 이름만 바꾼다
    #: (`YeastHRCloning` → `yeast_homology_min`/`_max`).
    #:
    #: 길이를 하드코딩하거나 `__init__`을 다시 쓰지 않는 이유는, arm이
    #: `__init__`(pydna overlap limit) · `prepare_parts`(접합부 설계) 두 곳에서 쓰이기 때문이다
    #: — 한 곳만 바꾸면 "접합부는 45bp로 만들었는데 조립은 20bp로 받는" 어긋남이 생긴다.
    #: 둘 다 범위다 — padding을 없앤 뒤로 이 범위가 유일한 탈출구이기 때문이다.
    ARM_PARAMS = ("gibson_overlap", "gibson_overlap_max")

    #: arm이 엉뚱한 자리에 붙어 **유저 범위 안에서는 못 피할 때** 마지막으로 넘어갈 수
    #: 있는 길이 범위. 유저 범위(soft)와 달리 이건 방법 자체의 한계다 —
    #: 15bp 미만은 pydna `Assembly`가 겹침으로 인정하지 않고, 60bp를 넘으면 arm이
    #: 프라이머 5' tail이라 주문 한계(`PRIMER_ORDER_MAX_NT`)에 걸린다.
    #: 범위를 벗어나 설계하면 **그 사실을 note로 남긴다**(`_resolve_arm_crosstalk`).
    ARM_HARD_RANGE = (15, 60)

    #: arm을 무엇과 비교할지. in vitro Gibson은 exonuclease가 fragment **말단만**
    #: 깎으므로 arm이 다른 fragment의 몸통 한복판에 어닐링할 수가 없다 → arm끼리만 본다.
    #: in vivo(`YeastHRCloning`)는 상동 구간이면 어디든 재조합하므로 몸통까지 본다.
    CROSSTALK_SCOPE = "arm"

    @property
    def arm_range(self) -> Tuple[int, int]:
        """(최소, 최대) 접합부 길이. `ARM_PARAMS`가 가리키는 params 필드에서 읽는다."""
        lo_name, hi_name = self.ARM_PARAMS
        return int(getattr(self.params, lo_name)), int(getattr(self.params, hi_name))

    def __init__(self, overlap_length: int = None, params=None):
        """Gibson Assembly 초기화

        Args:
            overlap_length: homology arm 길이. 미지정 시 params의 `ARM_PARAM` 필드에서 가져온다.
                (연구실마다 다른 값 — 요청 문서 D-1에서 수집 중)
            params: 설계 파라미터 한 벌(`params.DesignParams`). 미지정 시 기본값.

        Note:
            명시적 `overlap_length`가 params보다 우선한다(호출부가 특정 길이를 강제할 때).
        """
        from params import DEFAULT_PARAMS
        self.params = params if params is not None else DEFAULT_PARAMS
        # pydna `Assembly(limit=...)`는 "이 길이 이상 겹치면 인정"이라 **범위의 최소값**을
        # 써야 한다. 최대값을 쓰면 더 짧게 고른 접합부를 조립이 인정하지 않는다.
        self.overlap_length = (
            overlap_length if overlap_length is not None
            else int(getattr(self.params, self.ARM_PARAMS[0]))
        )
        # vector를 열 효소. params에 있으면 그걸 쓰고, 없으면 자동 선택.
        self.linearization_enzyme = getattr(self.params, "linearization_enzyme", None)
        #: 마지막 `prepare_fragments` 설계 결과(anneal 서열·Tm·경고). 리포트가 읽는다.
        #  `prepare_fragments`를 아직 안 돌렸으면 None이다.
        self.last_primer_design = None

    def prepare_fragments(
        self,
        insert: Dseqrecord,
        vector: Dseqrecord,
        target_start: int = None,
        target_end: int = None,
        feature_type: str = None,
        feature_label: str = None
    ) -> Tuple[str, str]:
        """Gibson용 프라이머 설계 (Basic 모듈 + overlap 추가)

        워크플로우:
            1. 입력 검증
            2. Vector linearization
            # TODO: 3번도 general 인 경우에서만 필요하지 않나? 아니면 바이오파츠의 타입을 spacer인걸 감지해서 해당 부분의 nick을 체크해야하나?
            3. Target feature 자동 감지 (또는 명시적 지정)
            4. Basic 프라이머 설계
            5. Overlap 추가
            6. 프라이머 검증

        Args:
            insert: Insert DNA
            vector: Vector DNA
            target_start: Insert 증폭 시작 위치
                기본값: None (자동 감지)
            target_end: Insert 증폭 끝 위치
                기본값: None (자동 감지)
            feature_type: 자동 감지할 feature 타입 ("CDS", "gene" 등)
                기본값: None (CDS 우선 검색)
            feature_label: 자동 감지할 feature label (예: "EGFP")
                기본값: None (label 검색 안 함)

        Returns:
            (forward_primer, reverse_primer) 튜플

        Raises:
            ValueError: 입력이 유효하지 않을 때
            RuntimeError: Vector linearization 실패

        Notes:
            1. target_start/end가 None이면 feature 자동 감지
            2. pEGFP.gb: EGFP CDS(769-1485)만 증폭 (ori, promoter 보존)
            3. Basic 모듈로 순수 증폭 프라이머 설계
            4. Vector와의 homology overlap 추가
            5. 짧은 insert (<50bp) 프라이머 단축

            archive/docs/development/PLAN_v2.md "반드시 지켜야 할 목록" 준수:
            - ori, promoter 같은 부분 건드리지 않음
            - MCS 부위 (target feature)만 증폭

        Example:
            >>> # 자동 감지 (CDS 우선)
            >>> gibson = GibsonCloning(overlap_length=15)
            >>> fwd, rev = gibson.prepare_fragments(my_insert, my_vector)
            >>>
            >>> # 명시적 지정
            >>> fwd, rev = gibson.prepare_fragments(
            ...     my_insert, my_vector,
            ...     target_start=769, target_end=1485
            ... )
            >>>
            >>> # Label로 지정
            >>> fwd, rev = gibson.prepare_fragments(
            ...     my_insert, my_vector,
            ...     feature_label="EGFP"
            ... )
        """
        # 1. 입력 검증
        self._validate_inputs(insert, vector)

        # 2. Vector linearization
        linear_vector = self._linearize_vector(vector)

        # 3. Target feature 자동 감지 (None이면)
        if target_start is None or target_end is None:
            target_start, target_end = self._find_target_feature(
                insert, feature_type, feature_label
            )

        # 4. Basic 모듈로 기본 프라이머 설계 (D-1: design_basic_primers general 엔진)
        #    anchor_ends=True → insert anneal을 target(=증폭할 insert fragment) 경계에 고정.
        #    길이/Tm은 조건 기반 가변, amplicon=target 정확(말단 외 서열 오염 없음).
        #    vector MCS homology는 아래 5단계에서 이 anneal '앞에' 별도 부착.
        from basic.primer_design import PrimerDesign

        basic = PrimerDesign.design_basic_primers(
            template=str(insert.seq),
            target_start=target_start,
            target_end=target_end,
            anchor_ends=True,
            check_specificity=False,
            params=self.params,   # 연구실 조건(목표 Tm·길이·GC)을 primer3까지 전달
        )
        if not basic.get("forward") or not basic.get("reverse"):
            raise RuntimeError(
                f"general-prep primer design failed (insert={insert.name}, "
                f"target={target_start}-{target_end}): {basic.get('warnings')}"
            )
        basic_fwd = basic["forward"]      # insert 시작 anneal (target_start부터)
        basic_rev = basic["reverse"]      # insert 끝 anneal (primer3 right=이미 RC)

        # 5. Overlap 추가
        # Forward: 벡터의 3' 끝 + 기본 프라이머
        vector_3prime = str(linear_vector.seq[-self.overlap_length:])
        gibson_fwd = vector_3prime + basic_fwd

        # Reverse: 벡터의 5' 끝 역상보 + 기본 프라이머
        vector_5prime = str(linear_vector.seq[:self.overlap_length].reverse_complement())
        gibson_rev = vector_5prime + basic_rev

        # 6. 짧은 insert 처리
        gibson_fwd, gibson_rev = self._validate_primer_length(
            insert, gibson_fwd, gibson_rev
        )

        # 2026-11-05(Y10-b) — 설계기가 **이미 계산한** anneal Tm·경고를 밖으로 흘린다.
        #   `vector_primers.csv`가 Tm/Ta를 적으려면 anneal 구간의 Tm이 필요한데,
        #   완성된 프라이머에서 tail을 슬라이스로 되짚으면 `_validate_primer_length`가
        #   길이를 건드린 경우 **틀린 구간을 재게 된다.**
        #   반환 계약(`AssemblyBase`는 튜플)은 그대로 두고 인스턴스에 실어 보낸다 —
        #   같은 질문(이 프라이머의 anneal Tm은?)에 두 개의 답을 만들지 않기 위함이다.
        self.last_primer_design = {
            "fwd_anneal": basic_fwd, "rev_anneal": basic_rev,
            "fwd_tm": basic.get("fwd_tm"), "rev_tm": basic.get("rev_tm"),
            "warnings": list(basic.get("warnings") or []),
        }

        return gibson_fwd, gibson_rev

    # ====================================================================== #
    # part → module 준비 (G3)
    #
    # `prepare_fragments`는 "1 insert + 1 vector"용이라 part N개를 서로 잇는 데 쓸 수 없다.
    # 접합부 설계(누가 homology를 운반하나 / padding을 넣나)는 `workflow.junction_plan`이
    # 맡고, 여기서는 그 명세를 **실제 주문 산출물**(프라이머 / 어닐링 올리고)로 바꾼다.
    #
    # `YeastHRCloning`이 이 클래스를 상속하므로 yeast 경로도 이 메서드를 그대로 받는다.
    # ====================================================================== #

    #: 맨몸 PCR 설계가 실패했을 때 anneal 시작점을 안쪽으로 몇 nt까지 밀어볼지.
    #:
    #: `anchor_ends=True`는 프라이머 5' 시작을 부품 경계에 고정하므로 자유도가 길이뿐이고,
    #: 그 결과 달성 가능한 Tm이 띄엄띄엄해져 선언 범위를 **건너뛰는** 구간이 생긴다
    #: (실측 실패 사유가 한 부품에서 "low tm"과 "high tm"이 동시에 뜬다).
    #:
    #: 시작점을 k nt 안으로 밀고 **잃은 k nt를 프라이머 5' tail에 되살리면** amplicon은
    #: 부품 전체 그대로이면서 primer3에 자유도가 생긴다.
    #:
    #: 상한 16은 실측 곡선에서 고른 값이다(긴 파츠 150개, arm 20bp 양쪽 운반 최악 조건):
    #:
    #: | shift≤ | 설계 성공률 | 프라이머 중앙 | 최대 | >60nt |
    #: |---|---|---|---|---|
    #: | 0 | 64% | 21nt | 45nt | 0% |
    #: | 8 | 89% | 22nt | 52nt | 0% |
    #: | **16** | **97%** | 23nt | **58nt** | **0%** |
    #: | 24 | 99% | 23nt | 68nt | 1% |
    #:
    #: 16에서 멈추는 이유는 **품질 계산 가능성**이다. 24로 올리면 +2%p를 얻지만 일부
    #: 프라이머가 primer3 계산 상한(60bp)을 넘어 hairpin/dimer 판정이 경고 전용(H1)으로
    #: 떨어진다 — 성공률을 위해 품질 눈을 감는 거래라서 하지 않았다.
    #:
    #: `primer_len_max`를 25→35로 올리는 것은 shift가 있으면 +1%p뿐이라 기본값을 바꾸지 않았다.
    MAX_ANNEAL_SHIFT = 16

    def prepare_parts(self, parts: List[Dseqrecord], arm: int = None,
                      backbone: Optional[Dseqrecord] = None) -> Tuple[List[Dict], object]:
        """part 리스트 → module 조립용 fragment 설계 (Gibson/HR 공용).

        Args:
            parts: 조립 순서대로의 부품 Dseqrecord 리스트. 양 끝 spacer가 필요하면
                호출자가 미리 붙여서 넘긴다.
            arm: homology arm 길이(bp). None이면 params의 `ARM_PARAM` 필드
                (Gibson은 `gibson_overlap`, yeast HR은 `yeast_homology`).
            backbone: 같은 반응에 들어가는 벡터. `CROSSTALK_SCOPE == "arm+body"`일 때만
                쓴다 — arm이 백본의 반복 구간에 내려앉는지 보는 용도다. 프라이머가
                벡터에 붙는 off-target은 **보지 않는다**(부품에 프라이머를 붙이는
                설계라 벡터는 템플릿이 아니다).

        Returns:
            (design_results, plan) — `design_results`는 fragment별 설계 dict 리스트로
            `GoldenGateModule.prepare_fragments`의 반환과 같은 자리에서 쓰이도록
            `part_name`/`part_id`/`part_type`/`full_sequence` 키를 맞췄다.
            `plan`은 `workflow.junction_plan.JunctionPlan` (기대 서열·접합부 내역).

        Raises:
            RuntimeError: 어떤 fragment의 프라이머 설계가 shift 폴백까지 실패했을 때,
                또는 `_check_plan`이 방법별 금지 조건을 발견했을 때.
        """
        from workflow.junction_plan import JunctionPlanner

        planner = JunctionPlanner(params=self.params)
        lo, hi = self.arm_range
        if arm is not None:
            lo = hi = int(arm)
        extra = []
        if self.CROSSTALK_SCOPE == "arm+body" and backbone is not None:
            extra.append((f"{backbone.name} backbone", str(backbone.seq).upper()))
        plan = planner.plan(
            parts, arm_min=lo, arm_max=hi,
            hard_range=self.ARM_HARD_RANGE,
            crosstalk_scope=self.CROSSTALK_SCOPE,
            extra_seqs=extra,
        )
        # 방법별 추가 판정(훅).
        self._check_plan(plan)

        results: List[Dict] = []
        for frag in plan.fragments:
            if frag.production == "oligo":
                results.append(self._design_oligo_fragment(frag))
            else:
                results.append(self._design_pcr_fragment(frag))
        return results, plan

    def _check_plan(self, plan) -> None:
        """접합부 설계 결과에 대한 **방법별 추가 판정** (훅).

        지금은 기본 구현이 아무것도 하지 않는다. 예전에는 여기서 "비인접 fragment가
        서열을 공유하면" Gibson은 경고, yeast는 차단으로 갈렸는데, 2026-11-06에 판정
        자체를 `JunctionPlanner`로 옮기면서 그 차이가 사라졌다:

          - **arm이 엉뚱한 자리에 붙는 것**은 `_resolve_arm_crosstalk`가 길이를 다시
            골라 해결하고, 끝까지 안 되면 양쪽 트랙 모두 설계불가로 끊는다.
          - **몸통끼리의 공유**는 arm 길이로 고칠 수 없어서(부품을 바꿔야 한다)
            양쪽 모두 기록만 남긴다.

        비교 범위의 차이(`CROSSTALK_SCOPE`)가 in vitro/in vivo 구분을 가져간다.
        훅은 방법별 추가 판정이 또 생길 자리로 남겨 둔다.

        Args:
            plan: `workflow.junction_plan.JunctionPlan`
        """
        return None

    def _design_oligo_fragment(self, frag) -> Dict:
        """짧은 fragment → 어닐링 올리고 한 쌍.

        `basic.oligo_design.design_annealed_oligos`를 그대로 쓴다. Golden Gate와 달리
        4nt overhang이 없으므로 `left_oh`/`right_oh`는 비우고 **fragment 전체**를
        이중가닥으로 주문한다 — homology arm도 서열의 일부로 함께 합성된다.
        짧은 부품이 양쪽 arm을 드는 경우가 여기다(그 덕에 양옆 PCR은 맨몸이 된다).
        """
        from basic.oligo_design import design_annealed_oligos

        o = design_annealed_oligos(
            part_seq=frag.seq, left_oh="", right_oh="", name=frag.name,
            params=self.params,
        )
        return {
            "part_name": frag.name,
            "part_id": frag.part_ids[0] if frag.part_ids else frag.name,
            "part_type": frag.part_types[0] if frag.part_types else "",
            "merged_part_ids": ",".join(frag.part_ids),
            "production": "oligo",
            "full_sequence": frag.seq,
            "target_seq": frag.seq,
            "forward_primer": "",
            "reverse_primer": "",
            "oligo_top": o["top"],
            "oligo_bottom": o["bottom"],
            "oligo_len": o["top_len"],
            "duplex_grade": o["grade"],
            "note": o["note"],
            "left_arm": frag.left_arm,
            "right_arm": frag.right_arm,
            "anneal_shift": (0, 0),
            "warnings": list(frag.notes),
        }

    def _design_pcr_fragment(self, frag) -> Dict:
        """긴 fragment → PCR 프라이머 한 쌍 (필요하면 anneal shift 폴백).

        **템플릿은 `frag.template` — 실제로 존재하는 분자 하나다.** 템플릿에 없는 것은
        전부 프라이머 5' tail로 들어간다:

            fwd = [left_arm][prefix][template[:k]][anneal]
            rev = [rc(right_arm)][rc(suffix)][rc(template[-k:])][anneal(이미 RC)]

        `prefix`/`suffix`는 이웃에서 흡수한 짧은 부품이고, `template[:k]`는 anneal
        시작점을 안쪽으로 민 만큼 되살리는 염기다. 둘 다 합성으로 들어가므로
        **amplicon은 `frag.seq` 그대로** 유지된다 — shift도 흡수도 산물을 바꾸지 않는다.
        """
        from basic.primer_design import PrimerDesign
        from Bio.Seq import Seq

        template = frag.template
        attempts: List[str] = []
        basic = None
        shift = None
        for k in range(self.MAX_ANNEAL_SHIFT + 1):
            # 양쪽을 같은 k만큼 민다 — 실측 평균 필요 shift가 1nt 미만이라 대칭으로 충분하다.
            if len(template) - 2 * k < self.params.amp_min:
                attempts.append(
                    f"k={k}: the remaining amplicon ({len(template) - 2 * k}bp) is below amp_min"
                    f"({self.params.amp_min}bp) -> cannot shift further"
                )
                break
            try:
                cand = PrimerDesign.design_basic_primers(
                    template=template,
                    target_start=k,
                    target_end=len(template) - k,
                    anchor_ends=True,
                    check_specificity=False,
                    params=self.params,
                )
            except Exception as exc:                      # primer3 내부 거부
                attempts.append(f"k={k}: {type(exc).__name__} {str(exc)[:60]}")
                continue
            if cand.get("forward") and cand.get("reverse"):
                basic, shift = cand, k
                break
            attempts.append(f"k={k}: {'; '.join(cand.get('warnings') or ['no reason given'])[:80]}")

        from workflow import design_events as _ev

        if basic is None:
            _ev.record("primer", _ev.KIND_BLOCKED, frag.name,
                       "primer design failed at every anneal shift",
                       f"shift 0-{self.MAX_ANNEAL_SHIFT}nt, Tm window "
                       f"[{self.params.tm_min}, {self.params.tm_max}]")
            raise RuntimeError(
                f"Primer design failed for fragment '{frag.name}' (template {len(template)}bp) - "
                f"no candidate even after shifting the anneal start by 0-{self.MAX_ANNEAL_SHIFT}nt.\n"
                f"  Attempts: {' / '.join(attempts[:4])}\n"
                f"  -> Widen the declared Tm window ([{self.params.tm_min}, {self.params.tm_max}]) "
                f"or merge this part with a neighbour and order it as an oligo."
            )
        if shift:
            # 산물은 그대로다 — 밀린 만큼을 프라이머 5' tail로 되살리기 때문이다.
            # 그래도 "선언한 조건으로는 설계가 안 됐다"는 사실이므로 남긴다.
            _ev.record("primer", _ev.KIND_DETOUR, frag.name,
                       f"anneal start shifted inward by {shift}nt on both ends",
                       f"no primer passed at shift 0 ({attempts[0] if attempts else 'see log'})")

        recovered_left = template[:shift]
        recovered_right = template[len(template) - shift:] if shift else ""
        fwd = frag.left_tail + recovered_left + basic["forward"]
        rev = (str(Seq(recovered_right + frag.right_tail).reverse_complement())
               + basic["reverse"])

        # 주문할 **전체 프라이머**를 품질 검증한다(tail 포함). 길면 H1이 경고 전용으로 본다.
        fwd_qc = self._validate_assembly_primer(fwd, len(frag.left_tail) + shift, rev)
        rev_qc = self._validate_assembly_primer(rev, len(frag.right_tail) + shift, fwd)

        notes = list(frag.notes)
        if shift:
            notes.append(
                f"The anneal start was shifted inward by {shift}nt on both sides "
                f"(the lost bases are restored in the primer tail - the product is identical)"
            )
        # ⚠️ **설계기가 한 말을 버리지 않는다.** `design_basic_primers`는 요청한
        #    Tm/GC/길이 범위로 primer3가 후보를 하나도 못 내면 **기본 범위로 폴백**하고
        #    그 사실을 `warnings`로 돌려준다. 그런데 여기서 안 꺼내가는 바람에
        #    Build Log·주문 시트 어디에도 안 떴다 — 유저 눈엔 "파라미터를 바꿨는데
        #    결과가 그대로"로만 보였다(2026-11-05 실측: 운영 h101/h102가 Tm 창만
        #    55~65 → 60~65로 다른데 완전히 동일한 프라이머를 냈다).
        notes += list(basic.get("warnings") or [])

        return {
            "part_name": frag.name,
            "part_id": frag.part_ids[0] if frag.part_ids else frag.name,
            "part_type": frag.part_types[0] if frag.part_types else "",
            "merged_part_ids": ",".join(frag.part_ids),
            "production": "pcr",
            "full_sequence": frag.seq,
            "target_seq": frag.seq,
            "template_seq": template,
            "forward_primer": fwd,
            "reverse_primer": rev,
            "fwd_tm": basic.get("fwd_tm"),
            "rev_tm": basic.get("rev_tm"),
            "fwd_gc": basic.get("fwd_gc"),
            "rev_gc": basic.get("rev_gc"),
            "left_arm": frag.left_arm,
            "right_arm": frag.right_arm,
            "absorbed_prefix": frag.prefix,
            "absorbed_suffix": frag.suffix,
            "anneal_shift": (shift, shift),
            "fwd_tier": fwd_qc["tier"],
            "rev_tier": rev_qc["tier"],
            "fwd_qc_engine": fwd_qc["engine"],
            "rev_qc_engine": rev_qc["engine"],
            "tier": _worst_tier(fwd_qc["tier"], rev_qc["tier"]),
            "warnings": notes
                      + [f"[fwd] {w}" for w in fwd_qc["messages"]]
                      + [f"[rev] {w}" for w in rev_qc["messages"]],
        }

    def _validate_assembly_primer(self, primer: str, tail_len: int,
                                  other: str) -> Dict:
        """주문할 전체 프라이머의 품질 — tier + 메시지 + 어느 엔진이 봤는지.

        ⚠️ **이 판정은 보고용이고 설계를 막지 않는다.** `prepare_parts`는 tier가 `red`여도
        산출물을 내보낸다 — tier는 기존 3-tier 체계(`ValidatePrimer`)를 그대로 받아오고,
        그 체계의 dimer 컷오프(−5.0 kcal/mol)가 아직 **wet-lab 확정 대기**라서 현재 기준으로는
        red가 흔하다.

        실측(실제 iGEM 파츠, PCR fragment 273개): **red 41% / yellow 45% / green 14%**.
        red의 주된 사유는 tail이 없는 맨몸 프라이머의 self-dimer다(tail이 붙은 프라이머는
        기존 정책대로 경고로만 처리된다 — tail 유래 dimer는 spacer 선정으로 못 고치므로).
        즉 red를 "이 설계는 쓸 수 없다"로 읽으면 안 되고, 컷오프가 확정되기 전까지는
        **상대 비교용 라벨**로 봐야 한다.

        Args:
            primer: 검증할 프라이머(tail 포함 전체)
            tail_len: 5' tail 길이(arm + 복원염기). dimer 정책 판정에 쓰인다
            other: 짝 프라이머(cross-dimer 검사용)

        Returns:
            ``{'tier', 'messages', 'engine'}``.
            `engine`은 `'viennarna'`면 **경고 전용 구간**(60bp 초과)이라는 뜻이다.
        """
        from basic.primer_design import ValidatePrimer
        from params import THERMO_MAX_LEN

        engine = "viennarna" if len(primer) > THERMO_MAX_LEN else "primer3"
        try:
            v = ValidatePrimer.validate_primer(
                primer, check_dimer_with=other, tail_len=tail_len, params=self.params,
            )
        except Exception as exc:
            # 품질 검증 실패가 설계를 막지는 않는다 — 산출물은 내보내고 사실을 남긴다.
            return {
                "tier": "unknown",
                "messages": [f"Quality check failed ({type(exc).__name__}: {str(exc)[:60]}) "
                             f"- the primers were designed but could not be judged."],
                "engine": engine,
            }

        msgs = list(v.get("violations") or []) + list(v.get("warnings") or [])
        if engine == "viennarna":
            msgs.append(
                f"The primer is {len(primer)}nt, above the primer3 calculation limit ({THERMO_MAX_LEN}nt) "
                f"- hairpin/dimer are judged **warning-only**."
            )
        return {"tier": v.get("tier", "unknown"), "messages": msgs, "engine": engine}

    def amplify_insert(
        self,
        insert: Dseqrecord,
        forward_primer: str,
        reverse_primer: str
    ) -> Dseqrecord:
        """PCR 시뮬레이션으로 overlap 추가된 fragment 생성

        PyDNA pcr() 함수를 사용하여 프라이머로 insert를 증폭합니다.
        결과물은 vector와의 overlap이 포함된 PCR product입니다.

        Args:
            insert: Template DNA (insert)
            forward_primer: Forward primer (overlap 포함)
            reverse_primer: Reverse primer (overlap 포함)

        Returns:
            Overlap이 추가된 PCR product (Dseqrecord)

        Raises:
            RuntimeError: PCR 실패 시

        Notes:
            - PyDNA의 pcr() 함수 사용
            - 프라이머는 prepare_fragments()로 설계된 것 사용
            - 결과물은 assemble()에서 vector와 함께 조립

        Example:
            >>> # 1. 프라이머 설계
            >>> fwd, rev = gibson.prepare_fragments(insert, vector)
            >>>
            >>> # 2. PCR 수행
            >>> pcr_product = gibson.amplify_insert(insert, fwd, rev)
            >>>
            >>> # 3. Assembly
            >>> results = gibson.assemble([linear_vector, pcr_product])
        """
        from pydna.amplify import pcr

        try:
            pcr_product = pcr(forward_primer, reverse_primer, insert)
        except Exception as e:
            raise RuntimeError(f"PCR simulation failed: {e}")

        return pcr_product

    def assemble(
        self,
        fragments: List[Dseqrecord],
        assembly_type: str = "circular"
    ) -> List[Dseqrecord]:
        """Gibson Assembly 수행

        Args:
            fragments: 조립할 fragment 리스트 (vector + inserts)
            assembly_type: "circular" 또는 "linear"

        Returns:
            조립된 결과물 리스트

        Raises:
            ValueError: assembly_type이 유효하지 않을 때
            RuntimeError: Assembly 실패

        Notes:
            PyDNA의 Assembly 클래스 사용
            KRIBB Line 519 참고

        Example:
            >>> results = gibson.assemble([vector, insert], assembly_type="circular")
            >>> results[0].circular
            True
        """
        from pydna.assembly import Assembly
        from pydna.common_sub_strings import terminal_overlap

        # Assembly_type 검증
        if assembly_type not in ["circular", "linear"]:
            raise ValueError(
                f"assembly_type must be 'circular' or 'linear'. "
                f"Got: {assembly_type}"
            )

        # Assembly 수행 (limit=overlap_length)
        asm = Assembly(
            fragments,
            algorithm=terminal_overlap,
            limit=self.overlap_length
        )

        # Circular 또는 Linear 조립
        if assembly_type == "circular":
            results = asm.assemble_circular()
        else:
            results = asm.assemble_linear()

        if not results:
            raise RuntimeError("Assembly failed: no product was produced.")

        return results

    def validate(self, result: Dseqrecord) -> Dict:
        """Gibson Assembly 결과 검증

        Args:
            result: 조립 결과물

        Returns:
            검증 결과 Dict
            {
                'valid': bool,
                'warnings': List[str],
                'errors': List[str],
                'details': Dict
            }

        Example:
            >>> validation = gibson.validate(result)
            >>> if validation['valid']:
            ...     print("검증 통과!")
        """
        warnings = []
        errors = []

        # 1. 최소 길이 체크
        if len(result.seq) < 100:
            warnings.append(f"The assembly product is short: {len(result.seq)} bp")

        # 2. Feature 개수 체크
        if len(result.features) == 0:
            warnings.append("There are no features")

        return {
            'valid': len(errors) == 0,
            'warnings': warnings,
            'errors': errors,
            'details': {
                'length': len(result.seq),
                'circular': result.circular,
                'features': len(result.features)
            }
        }


    # I-4: 프라이머 품질 검증은 prep(ValidatePrimer, C-3)이 담당. validate()는 조립 결과물 고유 검증
    # (길이/feature/circular)만 수행 → 책임 분리 완료(중복 없음).

    # ============================================================================
    # Helper 메서드
    # ============================================================================

    def _validate_inputs(self, insert: Dseqrecord, vector: Dseqrecord):
        """입력 검증

        Args:
            insert: Insert DNA
            vector: Vector DNA

        Raises:
            ValueError: 입력이 유효하지 않을 때

        Checks:
            - Vector는 circular여야 함
            - Insert는 linear여야 함 (권장)
        """
        if not isinstance(insert, Dseqrecord):
            raise ValueError("insert must be a Dseqrecord")

        if not isinstance(vector, Dseqrecord):
            raise ValueError("vector must be a Dseqrecord")

        # Vector는 circular여야 함 (linearization 전)
        if not vector.circular:
            print("Warning: the vector is not circular. It may already be linearized.")

        # Insert는 linear 권장 (하지만 필수는 아님)
        if insert.circular:
            print("Warning: the insert is circular. Linear DNA is recommended.")

    @staticmethod
    def _usable_here(vector: Dseqrecord) -> list:
        """이 벡터를 정확히 한 곳에서 자르는 효소 목록(실패 메시지용)."""
        from basic.sequence_analysis import LINEARIZATION_ENZYMES
        return [n for n in LINEARIZATION_ENZYMES
                if len(getattr(Restriction, n).search(vector.seq)) == 1]

    def _linearize_vector(self, vector: Dseqrecord) -> Dseqrecord:
        """Vector linearization

        제한효소를 선택해 vector를 linearize합니다.

        후보 목록은 `basic.sequence_analysis.LINEARIZATION_ENZYMES` **한 곳**에 있다.
        (2026-10-20 이전엔 여기 5종을 따로 하드코딩했고 그중 4종은 master에 없었다 —
         그 결과 DB 벡터의 79%만 선형화됐다. 지금 목록으로는 99.5%.)

        Args:
            vector: Circular vector

        Returns:
            Linearized vector

        Raises:
            RuntimeError: Linearization 실패

        Example:
            >>> linear = self._linearize_vector(circular_vector)
            >>> linear.circular
            False
        """
        # 이미 linear면 그대로 반환
        if not vector.circular:
            return vector

        from basic.sequence_analysis import LINEARIZATION_ENZYMES

        # 유저가 효소를 지정했으면 그것만 쓴다(지정 ≠ 추천 — 틀리면 이유를 말해준다).
        if getattr(self, "linearization_enzyme", None):
            name = self.linearization_enzyme
            enzyme = getattr(Restriction, name, None)
            if enzyme is None:
                raise RuntimeError(f"Vector linearization failed: unknown enzyme '{name}'")
            n = len(enzyme.search(vector.seq))
            if n != 1:
                usable = self._usable_here(vector)
                hint = (f"Enzymes usable on this vector: {', '.join(usable)}" if usable
                        else "No candidate cuts this vector exactly once - Gibson insertion is not possible")
                raise RuntimeError(
                    f"Vector linearization failed: the enzyme {name} you specified cuts this vector at "
                    f"{n} sites (it must cut exactly once). {hint}"
                )
            print(f"Vector linearization OK: {name} (user-specified)")
            return vector.cut(enzyme)[0]

        for enzyme_name in LINEARIZATION_ENZYMES:
            enzyme = getattr(Restriction, enzyme_name)
            if len(enzyme.search(vector.seq)) == 1:
                print(f"Vector linearization OK: {enzyme_name}")
                return vector.cut(enzyme)[0]

        # 전부 실패 = 이 벡터는 Gibson 삽입 불가. 몇 번씩 자르는지까지 보여준다.
        detail = ", ".join(
            f"{n}({len(getattr(Restriction, n).search(vector.seq))}x)"
            for n in LINEARIZATION_ENZYMES
        )
        raise RuntimeError(
            f"Vector linearization failed: none of the {len(LINEARIZATION_ENZYMES)} candidates cuts "
            f"this vector exactly once -> **this vector cannot take a Gibson insertion**. "
            f"Use a vector with an MCS, or consider Golden Gate insertion. "
            f"[cut counts] {detail}"
        )

    def _validate_primer_length(
        self,
        insert: Dseqrecord,
        fwd: str,
        rev: str
    ) -> Tuple[str, str]:
        """프라이머가 insert보다 짧은지 검증 및 조정

        짧은 insert (<50bp)의 경우 프라이머를 단축시킵니다.

        Args:
            insert: Insert DNA
            fwd: Forward 프라이머
            rev: Reverse 프라이머

        Returns:
            조정된 (fwd, rev) 프라이머

        Notes:
            KRIBB Line 467-473 참고
        """
        while len(insert.seq) < len(fwd) + len(rev):
            print(
                f"Warning: insert {insert.name} ({len(insert.seq)} bp) is shorter than the "
                f"primer. Shortening the primer: "
                f"{len(fwd)} bp → {len(fwd)-1} bp"
            )
            fwd = fwd[:-1]
            rev = rev[:-1]

        return fwd, rev

    def _optimize_overhang_length(self, sequence: str) -> int:
        # ⚠️ **DEAD CODE (2026-09-23 확인).** 이 메서드를 부르는 곳이 전체 코드베이스에 없다
        #    (유일한 언급은 아래 docstring 예제). 그래서 내부 상수(MIN/MAX_LENGTH 등)는
        #    파라미터화 대상에서 제외했다 — 파라미터로 빼면 "동작하는 설정"처럼 보인다.
        #    살릴지 지울지는 Gibson homology 길이 정책(회신 대상)과 함께 결정한다.
        """GC%와 Tm을 기반으로 overlap 길이 최적화

        15-25bp 범위에서 가장 적절한 overlap 길이를 선택합니다.

        Args:
            sequence: Overlap 후보 서열 (벡터 끝 서열)

        Returns:
            최적화된 overlap 길이 (15-25bp)

        Notes:
            최적 조건:
                - GC% 40-60% 선호
                - Tm 50-60°C 선호
                - 15bp 미만은 비특이적 결합 위험
                - 25bp 초과는 불필요한 비용 증가

        Example:
            >>> optimized_length = self._optimize_overhang_length(vector_end_seq)
            >>> overlap = vector_end_seq[:optimized_length]
        """
        MIN_LENGTH = 15
        MAX_LENGTH = 25
        TARGET_GC_MIN = 40.0
        TARGET_GC_MAX = 60.0
        TARGET_TM_MIN = 50.0
        TARGET_TM_MAX = 60.0

        # 서열이 MAX_LENGTH보다 짧으면 그 길이 사용
        max_possible_length = min(len(sequence), MAX_LENGTH)

        best_length = self.overlap_length  # 기본값으로 초기화
        best_score = float('-inf')

        # 15-25bp 범위에서 최적 길이 탐색
        for length in range(MIN_LENGTH, max_possible_length + 1):
            subseq = sequence[:length]

            # GC% 계산
            gc_count = subseq.count('G') + subseq.count('C')
            gc_percent = (gc_count / length) * 100

            # Tm 계산 (Wallace rule: Tm = 4(G+C) + 2(A+T))
            at_count = subseq.count('A') + subseq.count('T')
            tm = 4 * gc_count + 2 * at_count

            # 점수 계산 (GC%와 Tm이 목표 범위에 가까울수록 높은 점수)
            gc_score = 0
            if TARGET_GC_MIN <= gc_percent <= TARGET_GC_MAX:
                gc_score = 100  # 최적 범위
            else:
                # 범위 벗어나면 거리에 비례하여 감점
                gc_distance = min(
                    abs(gc_percent - TARGET_GC_MIN),
                    abs(gc_percent - TARGET_GC_MAX)
                )
                gc_score = max(0, 100 - gc_distance * 2)

            tm_score = 0
            if TARGET_TM_MIN <= tm <= TARGET_TM_MAX:
                tm_score = 100  # 최적 범위
            else:
                # 범위 벗어나면 거리에 비례하여 감점
                tm_distance = min(
                    abs(tm - TARGET_TM_MIN),
                    abs(tm - TARGET_TM_MAX)
                )
                tm_score = max(0, 100 - tm_distance)

            # 총점 (GC%와 Tm 동등하게 중요)
            total_score = gc_score + tm_score

            if total_score > best_score:
                best_score = total_score
                best_length = length

        return best_length

    def _check_primer_specificity(
        self,
        primer: str,
        template: Dseqrecord
    ) -> bool:
        """프라이머가 template에서 단일 위치에만 결합하는지 확인

        Off-target 결합을 방지하기 위한 특이성 검증입니다.

        Args:
            primer: 검증할 프라이머 서열
            template: 검증 대상 template DNA

        Returns:
            True if 프라이머가 단일 위치에만 결합
            False if 다수 위치에 결합 (off-target 위험)

        Notes:
            간단한 구현:
                - 프라이머 전체 서열이 template에서 몇 번 등장하는지 확인
                - Forward와 Reverse complement 모두 체크
                - 1번만 등장하면 특이적

            향후 개선:
                - Partial match도 체크 (예: 70% 이상 일치)
                - Secondary structure 고려
                - BLAST 기반 검증

        Example:
            >>> is_specific = self._check_primer_specificity(fwd_primer, insert)
            >>> if not is_specific:
            ...     print("경고: 프라이머가 여러 위치에 결합합니다")
        """
        from Bio.Seq import Seq

        # Template 서열을 문자열로 변환
        template_str = str(template.seq).upper()
        primer_upper = primer.upper()

        # Forward 방향 검색
        fwd_count = template_str.count(primer_upper)

        # Reverse complement 방향 검색
        primer_rc = str(Seq(primer_upper).reverse_complement())
        rev_count = template_str.count(primer_rc)

        # 총 결합 횟수
        total_count = fwd_count + rev_count

        # 정확히 한 번만 결합하면 특이적
        if total_count == 1:
            return True
        elif total_count == 0:
            print(f"Warning: the primer does not bind the template")
            return False
        else:
            print(
                f"Warning: the primer binds the template {total_count} times "
                f"(off-target risk)"
            )
            return False

    def _find_target_feature(
        self,
        insert: Dseqrecord,
        feature_type: str = None,
        feature_label: str = None
    ) -> Tuple[int, int]:
        """Insert에서 증폭할 feature 자동 감지

        GenBank feature에서 CDS, gene 등을 찾아서 target 범위를 결정합니다.
        pEGFP.gb 같은 파일에서 EGFP CDS만 증폭하도록 자동으로 범위를 찾습니다.

        Args:
            insert: Insert DNA (features 포함)
            feature_type: 찾을 feature 타입 ("CDS", "gene" 등)
                기본값: None이면 "CDS" 우선 검색
            feature_label: 찾을 feature label (예: "EGFP", "GFP")
                기본값: None이면 label 검색 안 함

        Returns:
            (start, end) 튜플
            - 찾으면: feature의 실제 위치
            - 못 찾으면: (0, 전체 길이) fallback

        Notes:
            검색 우선순위:
            1. feature_label이 지정되면 → 해당 label 찾기
            2. feature_type이 지정되면 → 해당 타입 찾기
            3. 둘 다 None이면 → "CDS" 타입 찾기
            4. 모두 실패하면 → 전체 서열 사용 (0, len)

            archive/docs/development/PLAN_v2.md "반드시 지켜야 할 목록" 준수:
            - ori, promoter 같은 부분 건드리지 않음
            - 실제 insert fragment (CDS)만 증폭

        Example:
            >>> # pEGFP.gb: EGFP CDS는 769-1485
            >>> start, end = self._find_target_feature(insert_gfp)
            >>> print(f"Target: {start}-{end}")
            Target: 769-1485

            >>> # label로 명시적 지정
            >>> start, end = self._find_target_feature(
            ...     insert_gfp,
            ...     feature_label="EGFP"
            ... )
        """
        # biopart module/multimodule은 통째로 vector에 넣는 insert다(annotations["insert_whole"]).
        # 이 표시를 무시하면 아래 CDS 자동감지가 module 안의 CDS **부품**을 삽입 대상으로
        # 오인해 그 구간만 증폭한다 → promoter/RBS/terminator/spacer가 construct에서 유실.
        # (module part feature에 실제 GenBank type이 붙기 전엔 전부 misc_feature라
        #  감지가 늘 실패해서 우연히 전체 서열 fallback으로 맞았던 자리다.)
        if insert.annotations.get("insert_whole"):
            print(f"Using the whole insert (biopart module/multimodule): 0-{len(insert.seq)}")
            return 0, len(insert.seq)

        # Features가 없으면 전체 사용
        if not insert.features:
            print("Warning: the insert has no features. Using the whole sequence.")
            return 0, len(insert.seq)

        # 1. feature_label로 검색 (가장 높은 우선순위)
        if feature_label:
            for feature in insert.features:
                label = feature.qualifiers.get('label', [''])[0]
                if label and feature_label.lower() in label.lower():
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    print(f"Feature found (label={label}): {start}-{end} ({end-start} bp)")
                    return start, end

        # 2. feature_type으로 검색
        if feature_type:
            for feature in insert.features:
                if feature.type.lower() == feature_type.lower():
                    start = int(feature.location.start)
                    end = int(feature.location.end)
                    label = feature.qualifiers.get('label', [''])[0]
                    print(f"Feature found (type={feature_type}, label={label}): {start}-{end} ({end-start} bp)")
                    return start, end

        # 3. 기본값: "CDS" 타입 검색
        for feature in insert.features:
            if feature.type.lower() == "cds":
                start = int(feature.location.start)
                end = int(feature.location.end)
                label = feature.qualifiers.get('label', [''])[0]
                print(f"Feature found (CDS, label={label}): {start}-{end} ({end-start} bp)")
                return start, end

        # 4. Fallback: 전체 서열 사용
        print("Warning: no suitable feature found. Using the whole sequence.")
        return 0, len(insert.seq)


# 모듈 메타데이터
__all__ = ['GibsonCloning']
