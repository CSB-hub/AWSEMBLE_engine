"""part 리스트 → Gibson/HR fragment 명세 (접합부 설계기).

## 무엇을 푸는가

Golden Gate는 각 부품을 **독립적으로** 준비한다 — 부품마다 양 끝에 효소 인식서열과 4bp
overhang을 붙이면 끝이고, 이웃이 무엇인지 알 필요가 없다. homology 조립(Gibson / yeast HR)은
접합부가 **인접한 두 부품이 공유하는 서열 구간**이라 부품 하나만 보고는 설계할 수 없다.
이 모듈이 그 "이웃 관계"를 푼다.

## 불변식 — PCR 템플릿은 실제로 존재하는 분자 하나여야 한다

이게 설계 전체를 지배하는 규칙이다. 템플릿에 없는 서열(homology arm, 짧은 부품, 어떤
추가 서열이든)은 **합성으로만** 들어갈 수 있다. 들어갈 수 있는 자리는 둘뿐이다:

1. **프라이머 5' tail** — 주문 길이 상한(`PRIMER_ORDER_MAX_NT`)까지
2. **어닐링 올리고** — 합성 상한(`OLIGO_SYNTH_MAX_BP`)까지

> 2026-10-25 이전 구현은 짧은 부품을 **이웃의 body에 병합**해서 `lp + sp`를 한 fragment로
> 묶고 그걸 PCR 템플릿으로 삼았다. 그런 분자는 존재하지 않으므로 물리적으로 불가능한
> 설계였다(실측상 모듈의 63~88%가 이 경로를 탔다). 그 규칙을 폐기하고 위 2갈래로 바꿨다.

## 짧은 부품(sp) 처리 — tail 흡수 / 독립 올리고

arm보다 짧은 부품은 자기 몸에서 접합부 구간을 떼어줄 수 없다. 두 길이 있다:

| 방법 | 조건 | 결과 |
|---|---|---|
| **tail 흡수** | `arm + len(sp) + anneal <= PRIMER_ORDER_MAX_NT` | 이웃 PCR의 프라이머 tail에 실린다. fragment가 하나 줄어든다 |
| **독립 올리고** | `arm + len(sp) + arm <= OLIGO_SYNTH_MAX_BP` | sp가 **양쪽 arm을 모두 들고** 자기 올리고가 된다 |

독립 올리고 쪽이 특히 깔끔한데, sp가 양쪽을 다 들면 **양옆의 긴 부품은 tail이 전혀 없는
맨몸 PCR**이 되기 때문이다. 그리고 양옆 fragment는 서로 공유하는 구간이 없어 교차도 없다:

```
lp1 | sp | lp2   (arm 45)
  A = lp1                          맨몸 PCR
  B = lp1[-45:] + sp + lp2[:45]    올리고 106nt
  C = lp2                          맨몸 PCR
  A↔B 공유 = lp1[-45:] · B↔C 공유 = lp2[:45] · A↔C 공유 = 없음
```

## 접합부 길이는 범위에서 고른다 — padding은 없다

예전에는 접합부 품질이 미달이면 **padding**을 끼웠다. 그걸 없앴다. 이유가 둘이다:

1. **생물학적으로 위험하다.** padding은 품질 미달인 접합부면 어디든 들어갔는데, 전형
   모듈에서 그 자리는 **RBS–CDS 경계**다. 거기 30~45bp를 끼우면 번역 개시가 깨진다 —
   조립은 성공하는데 발현이 죽는, 가장 나쁜 실패다.
2. **필요가 없다.** 길이를 범위로 주면 품질을 통과하는 창을 찾을 수 있다(실측 §아래).

대신 `arm_min ~ arm_max`에서 **품질을 통과하는 가장 짧은 길이**를 접합부마다 고른다.
짧을수록 프라이머·올리고가 싸고, 결정적이라 재현된다.

**길이를 바꿔도 최종 서열은 바뀌지 않는다** — 접합부 구간은 공유일 뿐 삽입이 아니다.
(회귀 테스트에서 `gibson_overlap` 20→30으로 바꿔도 construct SHA가 동일했다.)
그래서 탐색이 산출물을 흔들지 않는다.

## 품질 기준 — 프라이머 축을 빌려 쓰지 않는다

접합부 구간은 프라이머가 아니다. 상보 가닥과 **전체 길이로** 어닐링하므로 40~60bp면
GC와 무관하게 Tm이 충분하다. 그런데 예전 구현은 프라이머 anneal 기준(GC 판정창
35~60%, hairpin ΔG > −3.0)을 그대로 적용했고, 그 컷오프들은 ~20nt 기준이라 **길다는
이유만으로** 탈락시켰다.

실측(실제 파츠 말단 창):

| 기준 | 고정 45bp 통과 | range 40~60 탈출 |
|---|---|---|
| 구 기준 (GC 35~60 + ΔG>−3) | 68% | 78% |
| **GC 20~80만** | **99.3%** | **100%** |

hairpin ΔG ≤ −3 비율도 창 길이에 끌려간다 — 20bp에서 1%, 60bp에서 **20%**.
그래서 지금은 **GC(`gc_min`/`gc_max`, 설계 목표 축)만 차단 기준**으로 쓰고,
2차구조는 경고로만 남긴다. 진짜 차단 기준은 **접합부 유일성**이다.

> 2026-11-05: **프라이머 판정축도 같은 결론으로 따라왔다.** 여기서 위 실측으로
> 35~60%를 버린 뒤에도 `ValidatePrimer.validate_primer`의 GC 판정창은
> 35~60%로 남아 있어서, 효모 AT-rich 끝단 프라이머(GC 28%)가 red로 떨어졌다.
> 지금은 차단창이 설계창을 포함(20~80%)하고, 관례적 이상 구간(40~60%,
> 2026-11-05(Y11)에는 그 2층마저 없애고 **유저 선언창 `gc_min/gc_max`** 하나로
> 합쳤다 — 기본값이 20~80%라 위 표의 "GC 20~80만" 행과 같은 기준이다.
"""
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence, Tuple

from pydna.dseqrecord import Dseqrecord

#: 비인접 fragment가 공유해도 되는 최대 연속 서열 길이(bp).
#: 이보다 길게 공유하면 의도하지 않은 자리끼리 재조합할 수 있다.
#: ⚠️ **wet-lab 확정 대상**. Gibson 트랙은 경고, yeast 트랙은 하드 차단으로 쓴다.
JUNCTION_UNIQUE_MAX_BP = 15

#: 접합부 2차구조 경고 임계(kcal/mol). **차단하지 않는다** — 창이 길수록 ΔG가 자동으로
#: 내려가므로(20bp 1% → 60bp 20%) 차단 기준으로 쓰면 길이를 벌할 뿐이다.
JUNCTION_HAIRPIN_WARN_DG = -8.0


@dataclass
class FragmentSpec:
    """조립에 들어갈 fragment 하나의 설계 명세.

    Attributes:
        name: fragment 이름(흡수된 부품이 있으면 `+`로 이어 붙인다)
        part_ids: 이 fragment가 최종 산물에 기여하는 부품 id들(흡수분 포함).
            **순서는 body 안에서의 배치 순서**다 — prefix 부품들, 템플릿 부품, suffix 부품들.
        part_types: `part_ids`와 같은 순서의 part_type
        part_lens: `part_ids`와 같은 순서의 부품 길이(bp). module feature 좌표를
            되찾는 데 쓴다 — 서열을 다시 탐색하지 않고 누적합으로 구간이 나온다.
        template: **PCR 템플릿이 되는 실제 분자**의 서열. 올리고 fragment는 빈 문자열
        prefix: 템플릿 앞에 프라이머 tail로 덧붙는 흡수 부품 서열
        suffix: 템플릿 뒤에 덧붙는 흡수 부품 서열
        left_arm: 왼쪽 접합부에서 **이 fragment가 운반하는** 구간(운반자가 아니면 빈 문자열)
        right_arm: 오른쪽 접합부에서 이 fragment가 운반하는 구간
        production: `"oligo"`(어닐링 합성) 또는 `"pcr"`(증폭)
        notes: 설계 중 남긴 메모(tail 흡수, 제작방식 등)
    """

    name: str
    part_ids: List[str]
    part_types: List[str]
    part_lens: List[int] = field(default_factory=list)
    template: str = ""
    prefix: str = ""
    suffix: str = ""
    left_arm: str = ""
    right_arm: str = ""
    production: str = "pcr"
    notes: List[str] = field(default_factory=list)

    @property
    def body(self) -> str:
        """이 fragment가 **최종 module 서열에 기여하는** 구간."""
        return self.prefix + self.template + self.suffix

    @property
    def seq(self) -> str:
        """주문/증폭할 fragment 전체 서열 (운반하는 arm 포함)."""
        return self.left_arm + self.body + self.right_arm

    @property
    def part_spans(self):
        """`body` 안에서 부품마다 차지하는 구간 — `(part_id, part_type, start, end)`.

        `part_ids` 순서가 곧 body 배치 순서이므로 길이 누적합이면 충분하다.
        길이 정보가 없는(옛 경로) fragment는 빈 리스트를 돌려준다 — 호출자는
        feature 기입을 건너뛰면 된다.
        """
        if len(self.part_lens) != len(self.part_ids):
            return []
        out, pos = [], 0
        for pid, ptype, n in zip(self.part_ids, self.part_types, self.part_lens):
            out.append((pid, ptype, pos, pos + n))
            pos += n
        return out

    @property
    def order_len(self) -> int:
        """주문 길이(nt). 올리고면 합성 길이, PCR이면 amplicon 길이."""
        return len(self.seq)

    @property
    def left_tail(self) -> str:
        """forward 프라이머의 5' tail (PCR일 때). 템플릿에 없는 앞쪽 전부."""
        return self.left_arm + self.prefix if self.production == "pcr" else ""

    @property
    def right_tail(self) -> str:
        """reverse 프라이머의 5' tail (정방향 표기). 템플릿에 없는 뒤쪽 전부."""
        return self.suffix + self.right_arm if self.production == "pcr" else ""

    def to_record(self) -> Dseqrecord:
        """`module_assembly` 입력으로 쓸 Dseqrecord로 변환.

        `part_id`/`part_type`은 module feature 기입(`annotate_module_parts`)이 읽으므로
        반드시 실어 보낸다 — 끊기면 feature가 전부 `misc_feature`로 떨어진다.
        """
        rec = Dseqrecord(self.seq, name=self.name)
        rec.annotations["molecule_type"] = "DNA"
        rec.annotations["part_id"] = self.part_ids[0] if self.part_ids else self.name
        rec.annotations["part_type"] = self.part_types[0] if self.part_types else ""
        rec.annotations["merged_part_ids"] = ",".join(self.part_ids)
        rec.annotations["merged_part_types"] = ",".join(self.part_types)
        rec.annotations["production"] = self.production
        # module feature 기입(`annotate_module_parts_from_layout`)이 읽는 배치 정보.
        # 서열을 다시 탐색하지 않고 좌표를 계산하기 위한 것 — arm은 인접 fragment와
        # **겹치는** 구간이라 module 안에서 한 번만 세어야 하므로 body 기준으로 적는다.
        rec.annotations["part_spans"] = self.part_spans
        # arm은 **한쪽 fragment만** 운반한다(운반자 규칙). 접합부 하나를 그리려면
        # 양옆 둘 중 비어 있지 않은 쪽을 봐야 하므로 좌·우를 모두 싣는다.
        rec.annotations["left_arm_len"] = len(self.left_arm)
        rec.annotations["right_arm_len"] = len(self.right_arm)
        rec.annotations["body_len"] = len(self.body)
        return rec


@dataclass
class JunctionSpec:
    """접합부 하나의 설계 결과.

    Attributes:
        index: 몇 번째 접합부인지(0-based, fragment i와 i+1 사이)
        left_name / right_name: 양쪽 fragment 이름
        window: 두 fragment가 공유하는 서열 구간
        length: 고른 구간 길이(bp). 범위에서 품질을 통과하는 **최소** 길이
        carrier: 구간을 운반하는 쪽 — `"left"` / `"right"`
        gc: `window`의 GC(%)
        hairpin_dg: `window`의 2차구조 ΔG. **경고용**이고 차단하지 않는다
        tried: 탐색에서 떨어진 길이와 사유(진단용)
    """

    index: int
    left_name: str
    right_name: str
    window: str
    length: int
    carrier: str
    gc: float = 0.0
    hairpin_dg: Optional[float] = None
    tried: List[str] = field(default_factory=list)


@dataclass
class JunctionPlan:
    """`JunctionPlanner.plan()`의 산출물.

    Attributes:
        fragments: 조립 순서대로의 fragment 명세
        junctions: 내부 접합부 명세(fragment가 N개면 N-1개)
        module_seq: 조립이 성공했을 때 기대되는 module 서열(= 서열 오라클)
        arm_min / arm_max: 사용한 접합부 길이 범위
        warnings: 차단은 아니지만 알려야 하는 사항
        max_crosstalk: 비인접 fragment **몸통** 쌍의 최대 공통서열 길이(bp)
        crosstalk_checked: 몸통 교차를 비교한 비인접 쌍의 수 — 0건의 분모다
        crosstalk_pairs: 임계를 넘은 (왼쪽 이름, 오른쪽 이름, 길이, 방향) 목록.
            방향은 `"forward"` 또는 `"reverse-complement"`.
            ⚠️ **차단하지 않는다**(2026-11-06). 몸통끼리의 공유는 arm이 끼지 않으면
            Gibson에서는 어닐링할 수가 없고(exonuclease가 말단만 깎는다), yeast에서도
            arm 기준 판정이 따로 돌기 때문이다. 중복 부품 경고가 소리 없이 사라지지
            않도록 **기록만** 남긴다.
        arm_checked: arm 기준으로 비교한 쌍의 수 — 0건의 분모다
        arm_hits: arm이 엉뚱한 자리에 붙을 뻔한 건과 그 처리 결과
            (`{a, b, bp, orient, verdict, detail}`). `verdict`는
            `redesigned` / `redesigned (outside the declared range)` / `BLOCKED`.
    """

    fragments: List[FragmentSpec]
    junctions: List[JunctionSpec]
    module_seq: str
    arm_min: int
    arm_max: int
    warnings: List[str] = field(default_factory=list)
    max_crosstalk: int = 0
    crosstalk_checked: int = 0
    crosstalk_pairs: List[Tuple[str, str, int, str]] = field(default_factory=list)
    arm_checked: int = 0
    arm_hits: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def arm(self) -> int:
        """대표 접합부 길이 — 실제로 쓰인 길이들의 최댓값(리포트 표기용)."""
        return max((j.length for j in self.junctions), default=self.arm_min)

    def records(self) -> List[Dseqrecord]:
        """`module_assembly(method="gibson")`에 그대로 넣을 수 있는 레코드 리스트."""
        return [f.to_record() for f in self.fragments]

    def summary(self) -> Dict[str, Any]:
        """리포트용 요약 집계."""
        return {
            "fragments": len(self.fragments),
            "oligo": sum(1 for f in self.fragments if f.production == "oligo"),
            "pcr": sum(1 for f in self.fragments if f.production == "pcr"),
            "absorbed": sum(1 for f in self.fragments if f.prefix or f.suffix),
            "junctions": len(self.junctions),
            "arm_lengths": sorted({j.length for j in self.junctions}),
            "max_order_len": max((f.order_len for f in self.fragments), default=0),
            "max_crosstalk": self.max_crosstalk,
            "module_bp": len(self.module_seq),
        }


class JunctionPlanner:
    """순서가 있는 부품 리스트를 homology 조립용 fragment 명세로 바꾼다.

    서열만 받아 명세를 내는 **순수 설계 함수**다 — DB에 접근하지 않는다.

    Example:
        >>> plan = JunctionPlanner(params=prm).plan(parts, arm_min=40, arm_max=60)
        >>> records = plan.records()          # module_assembly 입력
        >>> plan.module_seq                   # 기대 서열(오라클)
    """

    def __init__(self, params=None):
        from params import DEFAULT_PARAMS

        self.params = params if params is not None else DEFAULT_PARAMS

    # ------------------------------------------------------------------ #
    # 공개 API
    # ------------------------------------------------------------------ #

    def plan(
        self,
        parts: Sequence[Dseqrecord],
        arm_min: Optional[int] = None,
        arm_max: Optional[int] = None,
        *,
        hard_range: Optional[Tuple[int, int]] = None,
        crosstalk_scope: str = "arm",
        extra_seqs: Optional[Sequence[Tuple[str, str]]] = None,
    ) -> JunctionPlan:
        """부품 리스트 → fragment 명세.

        Args:
            parts: 조립 순서대로의 부품 Dseqrecord 리스트. 양 끝 spacer가 필요하면
                **호출자가 미리 붙여서** 넘긴다.
            arm_min / arm_max: 접합부 길이 범위(bp). None이면 `params.gibson_overlap`
                하나로 고정한다. 같은 값이면 탐색 없이 그 길이를 쓴다.
            hard_range: arm이 엉뚱한 자리에 붙어 **유저 범위 안에서는 못 피할 때**
                마지막으로 넘어갈 수 있는 길이 범위. None이면 범위 밖 재설계를 하지
                않는다(= 유저 범위에서 실패하면 바로 설계불가).
            crosstalk_scope: `"arm"`(arm끼리만) 또는 `"arm+body"`(arm ↔ 다른
                fragment 몸통까지). in vitro/in vivo 차이다 — `GibsonCloning`/
                `YeastHRCloning`가 각자 선언한다.
            extra_seqs: arm이 붙으면 안 되는 **추가 서열**. `(이름, 서열)` 목록으로,
                yeast 트랙이 백본을 넘긴다.

        Returns:
            `JunctionPlan`

        Raises:
            ValueError: 부품이 2개 미만이거나 빈 서열이 있을 때, 범위가 뒤집혔을 때.
            RuntimeError: 어떤 접합부도 범위 안에서 품질을 통과하지 못했을 때,
                짧은 부품을 tail·올리고 어느 쪽으로도 실을 수 없을 때,
                또는 arm이 엉뚱한 자리에 붙는 것을 끝까지 못 피했을 때.
        """
        from basic.sequence_analysis import validate_dna_sequence

        if len(parts) < 2:
            raise ValueError(
                f"Junction design needs at least two parts (got {len(parts)})."
            )

        amin = int(arm_min if arm_min is not None else self.params.gibson_overlap)
        amax = int(arm_max if arm_max is not None else amin)
        if amin < 1 or amax < amin:
            raise ValueError(f"Invalid junction length range: {amin}-{amax}")

        items: List[Tuple[str, str, str]] = []
        for p in parts:
            seq = str(p.seq).upper()
            if not validate_dna_sequence(seq):
                raise ValueError(
                    f"Part '{p.name}' contains non-ATGC bases (N/IUPAC) or an empty sequence "
                    f"-> junction design is not possible."
                )
            pid = p.annotations.get("part_id") or p.name
            items.append((str(pid), p.annotations.get("part_type", ""), seq))

        warnings: List[str] = []
        fragments = self._lay_out(items, amin, amax, warnings)
        junctions = self._design_junctions(fragments, amin, amax, warnings)
        self._check_order_limits(fragments, warnings)

        plan = JunctionPlan(
            fragments=fragments,
            junctions=junctions,
            module_seq="".join(f.body for f in fragments),
            arm_min=amin,
            arm_max=amax,
            warnings=warnings,
        )
        # 순서가 중요하다 — arm 재설계가 arm 서열을 바꾸므로 **먼저** 돌리고,
        # 몸통 기록은 그 결과와 무관하므로(몸통은 안 바뀐다) 뒤에 둔다.
        self._resolve_arm_crosstalk(plan, hard_range, crosstalk_scope, extra_seqs)
        self._check_body_crosstalk(plan)
        return plan

    # ------------------------------------------------------------------ #
    # 1단계: fragment 배치 — 짧은 부품을 tail/올리고로 보낸다
    # ------------------------------------------------------------------ #

    def _lay_out(
        self,
        items: List[Tuple[str, str, str]],
        amin: int,
        amax: int,
        warnings: List[str],
    ) -> List[FragmentSpec]:
        """부품 리스트 → fragment 명세 리스트.

        `amin`보다 짧은 부품은 자기 몸에서 접합부 구간을 떼어줄 수 없다. 연속한 짧은
        부품을 하나로 묶은 뒤, 이웃 PCR의 tail에 실을 수 있으면 싣고(fragment가 줄어든다)
        아니면 독립 올리고로 뺀다.

        ## 체인 맨 끝은 바깥쪽 arm이 없다 (2026-11-04)

        접합부는 fragment **사이**에만 생긴다(`_design_junctions`). 그래서 체인의 양
        바깥면은 어떤 arm도 운반하지 않는다. 두 곳에서 이 사실을 쓴다:

        1. **맨 끝 spacer는 길이 무관하게 tail로 보낸다.** `vector_homology` 모드의 양 끝
           spacer는 길이가 `params.gibson_overlap`에서 나오고 `amin`도 같은 값에서
           나온다 — 즉 `len >= amin`이 **항상 참**이라 "긴 부품"으로 분류돼 독립
           fragment가 되고, 그 다음 `should_use_oligo`에 걸려 올리고로 빠져 나갔다.
           그러면 벡터 homology가 프라이머에 안 실려 **프라이머 시트만으로는 주문이
           불가능**해진다(2026-11-04 실측: Gibson only가 fragment 3개 + 올리고 2개로
           나오고 프라이머 시트에는 insert 맨몸만 남았다).
        2. **맨 끝 run을 흡수할 때 `amax`를 더하지 않는다.** 더하면 주문 상한 계산이
           실제보다 arm 하나만큼 커져(30+15=45nt인데 45+30+15=90nt로 계산) 들어갈 수
           있는 tail을 못 들어간다고 판정한다.

        Args:
            items: (part_id, part_type, seq) 리스트
            amin / amax: 접합부 길이 범위
            warnings: 처리 내역을 적을 대상

        Returns:
            fragment 명세 리스트(아직 arm은 비어 있다)
        """
        from params import OLIGO_SYNTH_MAX_BP, PRIMER_ORDER_MAX_NT

        def _edge_spacer(i: int) -> bool:
            """체인 맨 끝의 spacer이고 **프라이머 tail에 실릴 길이**인가.

            길이 조건을 빼면 안 된다 — `part_type == "spacer"`는 "접합부 재료"라는
            뜻만 담고 길이는 말해주지 않아서, 3000bp를 spacer로 넘기는 호출자가 있으면
            tail에도 올리고에도 못 들어가 설계가 죽는다. 맨 끝은 바깥 arm이 없으므로
            tail 비용은 spacer 길이 + anneal 최소뿐이다.
            """
            if items[i][1] != "spacer" or i not in (0, len(items) - 1):
                return False
            return len(items[i][2]) + self.params.primer_len_min <= PRIMER_ORDER_MAX_NT

        # 연속한 짧은 부품을 run으로 묶는다.
        runs: List[Tuple[str, List[Tuple[str, str, str]]]] = []
        for i, it in enumerate(items):
            kind = "sp" if (len(it[2]) < amin or _edge_spacer(i)) else "lp"
            if kind == "sp" and runs and runs[-1][0] == "sp":
                runs[-1][1].append(it)
            else:
                runs.append((kind, [it]))

        if not any(k == "lp" for k, _ in runs):
            total = sum(len(i[2]) for _k, g in runs for i in g)
            raise ValueError(
                f"Every part is shorter than the junction length ({amin}bp) (total {total}bp) - "
                f"there is no junction to build. "
                + (f"Ordering the whole thing as a single oligo pair removes the need to assemble."
                   if total <= OLIGO_SYNTH_MAX_BP else
                   f"Shorten the junction length or mix in longer parts.")
            )

        frags: List[FragmentSpec] = []
        pending_sp: Optional[List[Tuple[str, str, str]]] = None

        def mk_lp(group):
            ids = [g[0] for g in group]
            return FragmentSpec(
                name="+".join(ids), part_ids=ids, part_types=[g[1] for g in group],
                part_lens=[len(g[2]) for g in group],
                template=group[0][2], production="pcr",
            )

        pending_at_head = False   # 대기 중인 짧은 run이 체인 맨 앞인가(바깥 arm 없음)

        for idx, (kind, group) in enumerate(runs):
            if kind == "sp":
                pending_sp = group
                pending_at_head = idx == 0
                continue

            frag = mk_lp(group)
            sp_len = sum(len(g[2]) for g in pending_sp) if pending_sp else 0
            # tail 여유: arm(최대) + sp + anneal 최소가 주문 상한 안에 드는가.
            # 체인 맨 앞이면 그 바깥면에는 접합부가 없으므로 arm 몫을 빼고 센다(§docstring).
            arm_cost = 0 if pending_at_head else amax
            tail_fits = (
                pending_sp is not None
                and arm_cost + sp_len + self.params.primer_len_min <= PRIMER_ORDER_MAX_NT
            )
            if pending_sp and tail_fits:
                # 이 PCR fragment의 앞쪽 tail에 싣는다 — fragment가 하나 줄어든다.
                frag.prefix = "".join(g[2] for g in pending_sp)
                frag.part_ids = [g[0] for g in pending_sp] + frag.part_ids
                frag.part_types = [g[1] for g in pending_sp] + frag.part_types
                frag.part_lens = [len(g[2]) for g in pending_sp] + frag.part_lens
                frag.name = "+".join(frag.part_ids)
                frag.notes.append(
                    f"Short part(s) {', '.join(g[0] for g in pending_sp)} ({sp_len}bp) were carried "
                    f"on the forward primer tail"
                )
                warnings.append(
                    f"tail absorption: {', '.join(g[0] for g in pending_sp)} -> {group[0][0]} primer"
                )
                pending_sp = None
            elif pending_sp:
                # 독립 올리고 — 양쪽 arm을 모두 들게 되므로 양옆 PCR은 맨몸이 된다.
                sp_frag = self._sp_oligo(pending_sp, amax, warnings)
                frags.append(sp_frag)
                pending_sp = None
            frags.append(frag)

        if pending_sp:
            # 마지막이 짧은 부품으로 끝났다 — 앞 PCR의 뒤쪽 tail에 싣거나 올리고로.
            # 이 run은 정의상 체인 맨 뒤이므로 바깥면에 접합부가 없다 → arm 몫 0(§docstring).
            sp_len = sum(len(g[2]) for g in pending_sp)
            prev = frags[-1]
            if (prev.production == "pcr"
                    and sp_len + self.params.primer_len_min <= PRIMER_ORDER_MAX_NT):
                prev.suffix = "".join(g[2] for g in pending_sp)
                prev.part_ids = prev.part_ids + [g[0] for g in pending_sp]
                prev.part_types = prev.part_types + [g[1] for g in pending_sp]
                prev.part_lens = prev.part_lens + [len(g[2]) for g in pending_sp]
                prev.name = "+".join(prev.part_ids)
                prev.notes.append(
                    f"Short part(s) {', '.join(g[0] for g in pending_sp)} ({sp_len}bp) were carried "
                    f"on the reverse primer tail"
                )
                warnings.append(
                    f"tail absorption: {', '.join(g[0] for g in pending_sp)} -> {prev.name} primer"
                )
            else:
                frags.append(self._sp_oligo(pending_sp, amax, warnings))

        # 제작방식 확정 — 짧은 PCR fragment는 올리고가 더 싸고 확실하다.
        from basic.oligo_design import should_use_oligo

        for f in frags:
            if f.production == "pcr" and should_use_oligo(len(f.body), self.params):
                f.production = "oligo"
                f.notes.append(f"{len(f.body)}bp - oligo synthesis is cheaper and more reliable")

        # fragment가 1개인 것은 **정상일 수 있다.** 양 끝 spacer가 벡터 homology면
        # module 자체가 단일 PCR 산물이고, 조립 상대는 선형화된 벡터다(= 조각 2개,
        # 접합부 2개인 교과서적 Gibson). 그 경우 접합부는 삽입 단계에서 생기므로
        # 여기서 0개인 게 맞다. 부품이 전부 흡수돼 **아무 것도 안 남은** 경우만 막는다.
        if not frags:
            raise ValueError(
                f"No fragment is left. "
                f"Check whether the parts are too short for the junction length ({amin}bp)."
            )
        return frags

    def _sp_oligo(self, group, amax: int, warnings: List[str]) -> FragmentSpec:
        """짧은 부품 run → **양쪽 arm을 드는** 독립 올리고 fragment."""
        from params import OLIGO_SYNTH_MAX_BP

        ids = [g[0] for g in group]
        body = "".join(g[2] for g in group)
        total = 2 * amax + len(body)
        if total > OLIGO_SYNTH_MAX_BP:
            raise RuntimeError(
                f"There is nowhere to carry short part(s) {', '.join(ids)} ({len(body)}bp) - "
                f"too long for a primer tail, and with arms on both sides the oligo would be "
                f"{total}nt, above the synthesis limit ({OLIGO_SYNTH_MAX_BP}nt). "
                f"Shorten the junction length, or merge this part with a longer one."
            )
        warnings.append(
            f"Standalone oligo: {', '.join(ids)} ({len(body)}bp) - up to {total}nt including both arms"
        )
        return FragmentSpec(
            name="+".join(ids), part_ids=ids, part_types=[g[1] for g in group],
            part_lens=[len(g[2]) for g in group],
            template="", prefix=body, production="oligo",
            notes=[f"Short part {len(body)}bp - both arms are synthesised with it"],
        )

    # ------------------------------------------------------------------ #
    # 2단계: 접합부 — 운반자 결정 + 길이 탐색
    # ------------------------------------------------------------------ #

    def _design_junctions(
        self,
        frags: List[FragmentSpec],
        amin: int,
        amax: int,
        warnings: List[str],
    ) -> List[JunctionSpec]:
        """인접 fragment 쌍마다 공유 구간의 **출처·길이**를 정하고 arm을 부착한다."""
        junctions: List[JunctionSpec] = []
        for i in range(len(frags) - 1):
            left, right = frags[i], frags[i + 1]
            carrier = self._pick_carrier(left, right)
            # 구간은 **운반자의 반대편 body**에서 떼어낸다.
            source = right.body if carrier == "left" else left.body

            chosen: Optional[Tuple[int, str, float, Optional[float]]] = None
            tried: List[str] = []
            for L in range(amin, amax + 1):
                if len(source) < L:
                    tried.append(f"{L}bp: the partner fragment is only {len(source)}bp")
                    continue
                window = source[:L] if carrier == "left" else source[-L:]
                ok, reason, gc, dg = self._eval_window(window)
                if ok:
                    chosen = (L, window, gc, dg)
                    break
                tried.append(f"{L}bp: {reason}")

            if chosen is None:
                raise RuntimeError(
                    f"At junction {i} ({left.name} <-> {right.name}) no length in {amin}-{amax}bp "
                    f"passed quality.\n"
                    f"  Tried: {'; '.join(tried[:6])}\n"
                    f"  -> Widen the junction length range, or change the part at this position."
                )

            L, window, gc, dg = chosen
            if dg is not None and dg <= JUNCTION_HAIRPIN_WARN_DG:
                warnings.append(
                    f"Junction {i} ({left.name}<->{right.name}) has secondary structure "
                    f"(dG {dg:.1f}, {L}bp). This does not block - if assembly is slow, "
                    f"check the temperature and additive conditions."
                )
            if carrier == "left":
                left.right_arm = window
            else:
                right.left_arm = window
            junctions.append(
                JunctionSpec(
                    index=i, left_name=left.name, right_name=right.name,
                    window=window, length=L, carrier=carrier, gc=gc, hairpin_dg=dg,
                    tried=tried,
                )
            )
        return junctions

    @staticmethod
    def _pick_carrier(left: FragmentSpec, right: FragmentSpec) -> str:
        """접합부 구간을 어느 쪽이 운반할지.

        올리고는 주문 서열에 넣기만 하면 되므로(프라이머가 아니라 서열) 비용이 없고,
        PCR은 프라이머 tail이 되어 주문 길이 상한의 압박을 받는다. 그래서
        **올리고 쪽을 우선 운반자로** 고르고, 둘 다 PCR이면 amplicon이 짧은 쪽에 맡긴다.
        """
        l_oligo = left.production == "oligo"
        r_oligo = right.production == "oligo"
        if l_oligo and not r_oligo:
            return "left"
        if r_oligo and not l_oligo:
            return "right"
        if l_oligo and r_oligo:
            return "right"   # 둘 다 자유 → 한쪽 고정으로 재현성 확보
        return "left" if len(left.body) <= len(right.body) else "right"

    def _eval_window(self, window: str) -> Tuple[bool, str, float, Optional[float]]:
        """접합부 구간이 **차단 기준**을 통과하는가.

        차단은 GC 하나뿐이다(`gc_min`/`gc_max` — 설계 목표 축). 2차구조는 창이 길수록
        ΔG가 자동으로 내려가므로 차단에 쓰지 않고 값만 돌려준다(호출자가 경고로 쓴다).
        유일성은 fragment가 다 만들어진 뒤에 전체를 보고 판정한다
        (`_resolve_arm_crosstalk` — 걸리면 그쪽이 길이를 다시 고른다).

        Returns:
            (통과 여부, 사유, GC%, 2차구조 ΔG)
        """
        gc, dg = self._window_metrics(window)
        if not (self.params.gc_min <= gc <= self.params.gc_max):
            return (False,
                    f"GC {gc:.0f}% is outside the window ({self.params.gc_min:.0f}-{self.params.gc_max:.0f}%)",
                    gc, dg)
        return True, "", gc, dg

    @staticmethod
    def _window_metrics(window: str) -> Tuple[float, Optional[float]]:
        """구간의 (GC%, 2차구조 ΔG). 계산 불가 시 ΔG는 None.

        둘 다 `basic/`의 기존 프리미티브를 쓴다 — 여기서 다시 구현하면 사본이 하나 늘어난다.
        """
        from basic.secondary_structure import analyze_single_strand_structure
        from basic.sequence_analysis import calculate_gc_content

        # ⚠️ `calculate_gc_content`는 BioPython `gc_fraction()` 그대로라 **0~1 비율**을 준다.
        gc = float(calculate_gc_content(window)) * 100.0
        try:
            dg = float(analyze_single_strand_structure(window)["delta_g"])
        except Exception:
            dg = None
        return gc, dg

    # ------------------------------------------------------------------ #
    # 3단계: 주문 한계·교차 검사
    # ------------------------------------------------------------------ #

    def _check_order_limits(self, frags: List[FragmentSpec], warnings: List[str]) -> None:
        """주문 길이 상한을 넘는 fragment가 있는지 본다.

        tail 흡수 판정은 `amax` 기준으로 미리 했지만, 실제로 고른 arm 길이가 짧으면
        여유가 생기고 길면 빠듯해진다. 확정된 arm으로 한 번 더 확인한다.
        """
        from params import OLIGO_SYNTH_MAX_BP, PRIMER_ORDER_MAX_NT

        for f in frags:
            if f.production == "oligo":
                if f.order_len > OLIGO_SYNTH_MAX_BP:
                    warnings.append(
                        f"{f.name}: the oligo order length {f.order_len}nt exceeds the synthesis limit"
                        f" ({OLIGO_SYNTH_MAX_BP}nt)."
                    )
                continue
            for side, tail in (("forward", f.left_tail), ("reverse", f.right_tail)):
                total = len(tail) + self.params.primer_len_min
                if total > PRIMER_ORDER_MAX_NT:
                    warnings.append(
                        f"{f.name}: the {side} primer is at least {total}nt, above the order limit"
                        f" ({PRIMER_ORDER_MAX_NT}nt) "
                        f"(tail {len(tail)}nt)."
                    )

    # ------------------------------------------------------------------ #
    # arm 기준 교차 판정 + 재설계 래더
    # ------------------------------------------------------------------ #

    def _arm_source_index(self, plan: JunctionPlan, i: int) -> int:
        """접합부 i의 구간을 **떼어 온** fragment의 인덱스.

        그 fragment의 몸통에는 구간이 100% 들어 있다 — 설계 의도이므로 비교에서 뺀다.
        """
        return i + 1 if plan.junctions[i].carrier == "left" else i

    @staticmethod
    def _shared_run(a: str, b: str) -> Tuple[int, str]:
        """두 서열의 최장 공통 연속구간 길이와 걸린 방향.

        ⚠️ **정방향과 역상보를 둘 다 본다** (2026-11-05 Y12). 정방향만 보면 교차의
        **절반을 놓친다** — homology 조립은 이중가닥끼리 상동성을 찾는 반응이라
        한쪽이 뒤집혀 들어가도 재조합이 성립한다(pydna가 원형 산물을 역상보 회전으로
        2개 돌려주는 이유 — `YeastHRCloning.assemble` 참고).
        """
        from basic.sequence_analysis import longest_common_substring_len
        from Bio.Seq import Seq as _Seq

        if not a or not b:
            return 0, "forward"
        fwd = longest_common_substring_len(a, b)
        rev = longest_common_substring_len(a, str(_Seq(b).reverse_complement()))
        return (max(fwd, rev), "reverse-complement" if rev > fwd else "forward")

    def _arm_targets(self, plan: JunctionPlan, i: int, scope: str,
                     extra_seqs) -> List[Tuple[str, str]]:
        """접합부 i의 arm이 **붙으면 안 되는** 서열 목록 — `(이름, 서열)`.

        - 다른 접합부의 arm: 거기 붙으면 fragment가 엉뚱한 짝과 이어진다.
        - (`arm+body`) 다른 fragment의 몸통·백본: 세포 안에서는 상동 구간이면 어디든
          재조합하므로 arm이 몸통 한복판에 내려앉을 수 있다. in vitro Gibson은
          exonuclease가 **말단만** 깎아서 몸통에는 붙을 수가 없으므로 보지 않는다.

        같은 fragment의 반대쪽 arm도 비교 대상이다 — 몸통이 짧아 양쪽 창이 겹치면
        그 fragment가 뒤집힌 짝으로 이어질 수 있다(길이를 줄이면 풀린다).
        """
        out: List[Tuple[str, str]] = []
        for j, jn in enumerate(plan.junctions):
            if j != i:
                out.append((f"junction {j} arm ({jn.left_name}<->{jn.right_name})", jn.window))
        if scope == "arm+body":
            src = self._arm_source_index(plan, i)
            for k, f in enumerate(plan.fragments):
                if k != src:
                    out.append((f"{f.name} body", f.body))
            for nm, sq in (extra_seqs or []):
                out.append((nm, sq))
        return out

    def _worst_arm_hit(self, window: str, targets) -> Tuple[int, str, str]:
        """창 하나를 대상 전부와 비교해 가장 나쁜 (bp, 대상 이름, 방향)."""
        worst, who, orient = 0, "", "forward"
        for nm, sq in targets:
            n, o = self._shared_run(window, sq)
            if n > worst:
                worst, who, orient = n, nm, o
        return worst, who, orient

    def _resolve_arm_crosstalk(self, plan: JunctionPlan,
                               hard_range: Optional[Tuple[int, int]],
                               scope: str, extra_seqs) -> None:
        """arm이 엉뚱한 자리에 붙는지 보고, 걸리면 **접합부 단위로 다시 설계한다**.

        래더(요청 2026-11-06):
          ① 유저가 선언한 범위 안에서 다른 길이로 다시 고른다.
          ② 안 되면 `hard_range` 안에서 범위 **밖** 길이까지 간다 — 산출물에
             "왜 범위를 벗어났는지" note로 남긴다.
          ③ 그래도 안 되면 **설계불가**. 어느 arm이 어디에 몇 bp로 붙는지 적는다.

        접합부 하나 때문에 전체 arm을 밀지 않는다 — 걸린 접합부만 바꾼다.
        그리고 한 접합부를 바꾸면 다른 접합부의 비교 대상도 바뀌므로, 바꾼 결과를
        반영하며 **순차로** 푼다(그리디). 임계는 `JUNCTION_UNIQUE_MAX_BP`.
        """
        from workflow import design_events as events

        juncs = plan.junctions
        if not juncs:
            return

        amin, amax = plan.arm_min, plan.arm_max
        hlo, hhi = hard_range if hard_range else (amin, amax)
        hlo, hhi = min(hlo, amin), max(hhi, amax)

        # 분모: 각 arm × 비교 대상 수. "0건"이 결과로 읽히려면 몇을 봤는지가 필요하다.
        plan.arm_checked = sum(
            len(self._arm_targets(plan, i, scope, extra_seqs)) for i in range(len(juncs)))

        for i, jn in enumerate(juncs):
            targets = self._arm_targets(plan, i, scope, extra_seqs)
            worst, who, orient = self._worst_arm_hit(jn.window, targets)
            if worst < JUNCTION_UNIQUE_MAX_BP:
                continue

            left, right = plan.fragments[i], plan.fragments[i + 1]
            source = right.body if jn.carrier == "left" else left.body
            label = f"junction {i} ({jn.left_name}<->{jn.right_name})"

            # ① 유저 범위 → ② 범위 밖(hard). 각 구간 안에서는 짧은 쪽부터 — 짧은 창이
            #    공유 구간을 벗어날 가능성이 높고 주문 길이에도 유리하다.
            in_range = [L for L in range(amin, amax + 1) if L != jn.length]
            out_range = [L for L in range(hlo, hhi + 1) if not (amin <= L <= amax)]
            picked = None
            for tier, lengths in (("in", in_range), ("out", out_range)):
                for L in lengths:
                    if len(source) < L:
                        continue
                    cand = source[:L] if jn.carrier == "left" else source[-L:]
                    ok, _reason, gc, dg = self._eval_window(cand)
                    if not ok:
                        continue
                    n, _w, _o = self._worst_arm_hit(cand, targets)
                    if n < JUNCTION_UNIQUE_MAX_BP:
                        picked = (tier, L, cand, gc, dg)
                        break
                if picked:
                    break

            if picked is None:
                plan.arm_hits.append({
                    "a": label, "b": who, "bp": worst, "orient": orient,
                    "verdict": "BLOCKED",
                    "detail": (f"No junction length in {hlo}-{hhi}bp kept this arm from "
                               f"binding {who}."),
                })
                events.record("junction", events.KIND_BLOCKED, label,
                              f"arm binds {who} ({worst}bp, {orient})",
                              f"no length in {hlo}-{hhi}bp avoided it")
                raise RuntimeError(
                    f"Design not possible - at {label} the homology arm also matches "
                    f"{who} over {worst}bp ({orient}, threshold {JUNCTION_UNIQUE_MAX_BP}bp), "
                    f"and no junction length between {hlo} and {hhi}bp avoided it.\n"
                    f"  The arm would anneal at that wrong position as readily as at the "
                    f"intended junction, so the assembly can give a rearranged product.\n"
                    f"  -> Replace the part at this position with a different one, or reorder "
                    f"the modules so the two matching positions are not in the same reaction."
                )

            tier, L, window, gc, dg = picked
            # 창을 바꿨으니 arm을 다시 붙인다 — 이전 창은 지운다(carrier는 그대로).
            if jn.carrier == "left":
                left.right_arm = window
            else:
                right.left_arm = window
            # ⚠️ 바꾸기 **전에** 원래 길이를 잡아 둔다 — 아래 경고 문구가 "원래
            #   몇 bp였는데 몇 bp로 바꿨다"를 말해야 하는데, 먼저 덮어쓰면 양쪽이
            #   같은 숫자로 찍힌다.
            old_len = jn.length
            jn.tried.append(f"{old_len}bp: the arm matched {who} over {worst}bp ({orient})")
            jn.window, jn.length, jn.gc, jn.hairpin_dg = window, L, gc, dg

            if tier == "in":
                verdict = "redesigned"
                why = f"redesigned to {L}bp within the declared range {amin}-{amax}bp"
                plan.warnings.append(
                    f"{label}: the {old_len}bp arm also matched {who} over {worst}bp "
                    f"({orient}) - the junction was redesigned to {L}bp inside the declared "
                    f"range ({amin}-{amax}bp). No action needed."
                )
                events.record("junction", events.KIND_DETOUR, label,
                              f"junction length -> {L}bp", f"the arm matched {who} ({worst}bp)")
            else:
                verdict = "redesigned (outside the declared range)"
                why = f"redesigned to {L}bp, outside the declared range {amin}-{amax}bp"
                plan.warnings.append(
                    f"{label}: no length in the declared range ({amin}-{amax}bp) kept the arm "
                    f"from matching {who} over {worst}bp ({orient}), so it was designed at "
                    f"{L}bp - **outside the range you set** (allowed limits {hlo}-{hhi}bp). "
                    f"The assembly is expected to work, but the junction length is not the one "
                    f"you declared; widen the range if you want this inside it."
                )
                events.record("junction", events.KIND_LIMIT, label,
                              f"junction length -> {L}bp (outside {amin}-{amax}bp)",
                              f"the arm matched {who} ({worst}bp) at every in-range length")
            plan.arm_hits.append({
                "a": label, "b": who, "bp": worst, "orient": orient,
                "verdict": verdict, "detail": why,
            })

    @staticmethod
    def _check_body_crosstalk(plan: JunctionPlan) -> None:
        """비인접 fragment의 **몸통**이 서열을 얼마나 공유하는지 기록한다.

        인접 쌍은 공유하는 것이 설계 의도이므로 제외한다. 몸통끼리 길게 공유하면
        (실측: 같은 부품을 두 자리에 쓰면 81bp) 세포 안에서는 의도하지 않은 자리끼리
        재조합할 수 있다.

        ⚠️ **차단하지 않는다**(2026-11-06). 2개를 바꿨다:
          - 비교 대상이 `fragment.seq`(arm 포함) → **`fragment.body`**. arm이 끼는
            위험은 `_resolve_arm_crosstalk`가 재설계로 처리하므로, 여기서 또 세면
            같은 사건이 두 번 보고된다.
          - yeast 하드 차단 → 기록. 몸통끼리의 공유는 **arm 길이로 고칠 수 없다**
            (부품을 바꿔야 한다). 고칠 수 없는 것으로 설계를 막는 대신, 원인과 조치를
            적어 남긴다. Gibson에서는 exonuclease가 말단만 깎으므로 애초에 어닐링
            자체가 안 된다.
        """
        frags = plan.fragments
        worst = 0
        checked = 0
        for i in range(len(frags)):
            for j in range(i + 2, len(frags)):
                checked += 1
                shared, orient = JunctionPlanner._shared_run(frags[i].body, frags[j].body)
                worst = max(worst, shared)
                if shared >= JUNCTION_UNIQUE_MAX_BP:
                    plan.crosstalk_pairs.append(
                        (frags[i].name, frags[j].name, shared, orient))
        plan.max_crosstalk = worst
        plan.crosstalk_checked = checked
        for a, b, n, orient in plan.crosstalk_pairs:
            plan.warnings.append(
                f"Non-adjacent fragments {a} <-> {b} share {n}bp of body sequence "
                f"({orient}, threshold {JUNCTION_UNIQUE_MAX_BP}bp). The homology arms were "
                f"kept clear of each other, so the assembly is not blocked - but in vivo "
                f"these positions can still recombine with each other. Check whether the "
                f"same part is used in two places."
            )


__all__ = [
    "FragmentSpec",
    "JunctionSpec",
    "JunctionPlan",
    "JunctionPlanner",
    "JUNCTION_UNIQUE_MAX_BP",
    "JUNCTION_HAIRPIN_WARN_DG",
]
