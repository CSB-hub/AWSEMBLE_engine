"""Yeast Homologous Recombination (HR) 클로닝.

## 2026-10-24: 스텁 → 구현

2026-09-23에 클래스 골격만 만들어 두고(거동이 Gibson과 동일) 확정 필요 4가지를 적어 뒀다.
그 4가지가 정해지고, 선행 병목이던 **part→module homology prep**이 G1~G3으로 열리면서
이제 yeast 전용 거동을 갖는다.

| 확정 항목 | 결정 |
|---|---|
| homology arm 길이 | **40~60bp 범위** (하드 30~60). `params.yeast_homology_min/_max` — Gibson과 **다른 축** |
| 벡터 선형화 방식 | **제한효소 절단**. DB의 yeast 셔틀 벡터 후보 31개가 **전부** 단일절단 효소를 4~13개씩 보유해 PCR 선형화를 만들 필요가 없었다 |
| 붙는 단계 | `part → module` (G3 prep 상속) + `module → vector` 삽입. 둘 다 같은 arm을 쓴다 |
| 검증 기준 | 서열 오라클(흉터가 없어 기대값을 분석적으로 계산 가능) + 산물 유일성(회전·역상보 정규화) + **접합부 유일성 하드 검사** + 벡터 적격성 |

## 왜 Gibson을 기반으로 하나

두 방법 모두 **서열 상동성으로 조각을 잇는다**. 제한효소 인식서열이 필요 없고
(sequence-independent), 조각 말단에 상동 구간을 붙여 설계한다는 점이 같다.
그래서 fragment 준비·프라이머 설계·산물 검증의 뼈대를 그대로 쓴다.
`ARM_PARAMS`만 바꾸면 접합부 설계기와 pydna overlap limit이 함께 yeast 길이로 움직인다.

## Gibson과 실제로 다른 점

| 축 | Gibson | Yeast HR | 코드 |
|---|---|---|---|
| 반응 장소 | 시험관 | **효모 세포 안(in vivo)** | `validate()` 문구 |
| homology 길이 | 15~25bp 고정 | **40~60bp 범위** (품질 통과하는 최소를 접합부마다) | `ARM_PARAMS` |
| 조립 확인 시점 | 조립 직후 | **형질전환·회수 후** | `validate()` 문구 |
| 접합부 교차 | 경고 | **차단** | `_check_plan()` |
| 벡터 요건 | 제약 없음 | **효모 복제원점 + 선택마커** | `check_yeast_vector()` |

### 접합부 교차를 왜 차단까지 하나

Gibson은 시험관 산물을 전기영동·시퀀싱으로 확인한 뒤 다음 단계로 간다. 엉뚱한 조립이
섞여도 **그 자리에서 보인다**. yeast HR은 형질전환 자체가 조립 단계라서, 잘못된 재조합은
콜로니를 키우고 플라스미드를 회수한 **뒤에야** 드러난다. 확인 시점이 뒤에 있으므로
설계 단계에서 막는 쪽이 맞다.

실측 근거: 같은 부품을 두 자리에 쓰면 비인접 fragment가 **81bp를 공유**한다
(`develop_log/ljw/261024_gibson_part_module_and_yeast_plan.md` §1-2).

### 산물은 "기대 산물"이다

in vivo라서 이 코드가 내놓는 조립 결과는 **시뮬레이션상 기대 산물**이다. 실제로는 효모에서
회수해 확인해야 한다. 리포트 문구가 그 전제를 유지한다(`IN_VIVO_NOTE`).
"""
from typing import Dict, List, Optional

from pydna.dseqrecord import Dseqrecord

from .gibson import GibsonCloning

#: 효모에서 플라스미드가 복제되게 하는 요소. feature label·설명에서 이 이름들을 찾는다.
#:   ARS/CEN  — 염색체 복제원점 + 중심체 (저복제수, 안정)
#:   2micron  — 2μ plasmid 복제원점 (고복제수)
YEAST_ORI_KEYWORDS = (
    "ARS", "ARSH4", "ARS1", "ARS4", "CEN", "CEN6", "CEN4", "CEN/ARS",
    "2micron", "2-micron", "2 micron", "2u ori", "2μ",
)

#: 효모 선택마커. 영양요구성(auxotrophic) + 항생제 저항성 둘 다 포함한다.
YEAST_MARKER_KEYWORDS = (
    "URA3", "HIS3", "LEU2", "TRP1", "ADE2", "MET15", "LYS2",
    "kanMX", "natMX", "hphMX", "KanMX", "G418", "nourseothricin", "hygromycin",
)

#: 벡터가 yeast 셔틀임을 시사하는 이름·설명 키워드(보조 증거).
YEAST_VECTOR_HINTS = (
    "pRS", "YEp", "YCp", "YIp", "pAUR", "cerevisiae", "yeast", "shuttle",
)


def check_yeast_vector(vector: Dseqrecord) -> Dict:
    """이 벡터가 효모에서 복제·선별될 수 있는가.

    yeast HR의 최종 산물은 효모 안에서 원형으로 만들어지고 **그 세포에서 복제·선별**되어야
    회수할 수 있다. 효모 복제원점이나 선택마커가 없는 벡터(예: 일반 E. coli 클로닝 벡터)를
    쓰면 조립 자체는 시뮬레이션상 성립하는데 **실험에서 콜로니가 나오지 않는다.**

    판정 재료는 **feature label + 레코드 이름/설명**이다. 서열 모티프는 쓰지 않는다 —
    ARS consensus(`WTTTAYRTTTW`)를 DB 벡터 373개로 실측했더니 yeast 후보군의 84%에서
    히트가 나오지만(중앙 1개) 비-yeast의 10%에서도 나오고, **밀도 상위 목록에 비-yeast가
    섞여** 개체 수준 판별력이 없었다. 11bp 축퇴 모티프는 수 kb 서열에서 우연히 나온다.

    Args:
        vector: 검사할 벡터 Dseqrecord

    Returns:
        ``{'verdict', 'ori', 'marker', 'hints', 'evidence', 'warning'}``

        - ``verdict``: ``'ok'``(복제원점·마커 둘 다 확인) / ``'partial'``(하나만) /
          ``'unknown'``(메타데이터가 없어 판단 불가)
        - ``warning``: 사용자에게 보여줄 문구. ``'ok'``면 빈 문자열

    Note:
        ``'unknown'``은 **"없다"가 아니라 "확인할 수 없다"** 다. DB에서 서열만 가져온
        벡터는 feature가 없어서 대개 여기에 해당한다 — 그래서 이 판정은 **차단이 아니라
        경고**다. GenBank 파일로 등록한 벡터는 feature가 있어 판정이 된다.
    """
    haystack: List[str] = []
    for feat in getattr(vector, "features", []) or []:
        for key in ("label", "gene", "note", "product", "standard_name"):
            for val in feat.qualifiers.get(key, []) or []:
                haystack.append(str(val))
        haystack.append(str(getattr(feat, "type", "")))

    ann = getattr(vector, "annotations", {}) or {}
    meta = [str(vector.name or ""), str(getattr(vector, "description", "") or "")]
    for key in ("comment", "keywords", "source", "organism"):
        v = ann.get(key)
        if isinstance(v, (list, tuple)):
            meta.extend(str(x) for x in v)
        elif v:
            meta.append(str(v))

    def _find(words, pool):
        low = [p.lower() for p in pool]
        return sorted({w for w in words if any(w.lower() in p for p in low)})

    searchable = haystack + meta
    ori = _find(YEAST_ORI_KEYWORDS, searchable)
    marker = _find(YEAST_MARKER_KEYWORDS, searchable)
    hints = _find(YEAST_VECTOR_HINTS, meta)

    evidence = []
    if ori:
        evidence.append(f"복제원점 후보: {', '.join(ori)}")
    if marker:
        evidence.append(f"선택마커 후보: {', '.join(marker)}")
    if hints:
        evidence.append(f"이름·설명 단서: {', '.join(hints)}")

    if ori and marker:
        verdict, warning = "ok", ""
    elif ori or marker or hints:
        missing = "선택마커" if ori else ("복제원점" if marker else "복제원점과 선택마커")
        verdict = "partial"
        warning = (
            f"⚠️ 벡터 '{vector.name}'에서 효모 {missing}를 확인하지 못했습니다 "
            f"({'; '.join(evidence)}). yeast HR 산물은 효모에서 복제·선별되어야 회수할 수 "
            f"있습니다 — 이 벡터가 효모 셔틀 벡터인지 확인하세요."
        )
    else:
        verdict = "unknown"
        warning = (
            f"⚠️ 벡터 '{vector.name}'에 feature·설명이 없어 효모 복제원점(ARS/CEN 또는 2μ)과 "
            f"선택마커를 **확인할 수 없습니다**(없다는 뜻은 아닙니다). yeast HR 산물은 효모에서 "
            f"복제·선별되어야 회수됩니다 — 셔틀 벡터인지 직접 확인하거나, feature가 있는 "
            f"GenBank 파일로 등록하면 자동 판정됩니다."
        )

    return {
        "verdict": verdict,
        "ori": ori,
        "marker": marker,
        "hints": hints,
        "evidence": evidence,
        "warning": warning,
    }


class YeastHRCloning(GibsonCloning):
    """Yeast HR 클로닝 — 세포 안에서 상동성으로 잇는다.

    `GibsonCloning`을 상속해 fragment 준비·프라이머 설계·조립 시뮬레이션을 공유하고,
    yeast 고유 거동 3가지만 재정의한다: **arm 길이 축 / 접합부 교차 차단 / in vivo 전제**.

    Attributes:
        ARM_PARAMS: `("yeast_homology_min", "yeast_homology_max")`.
            이 두 이름으로 `__init__`(pydna overlap limit)과 `prepare_parts`(접합부 설계)가
            함께 yeast 길이로 움직인다.

    Example:
        >>> y = AssemblyFactory.create_assembler("yeast_hr", params=prm)
        >>> y.overlap_length          # params.yeast_homology_min (기본 40)
        40
        >>> design, plan = y.prepare_parts(parts)     # 45bp 접합부로 설계
    """

    ARM_PARAMS = ("yeast_homology_min", "yeast_homology_max")

    #: 산출물·리포트에 유지해야 하는 in vivo 전제.
    IN_VIVO_NOTE = (
        "yeast HR은 형질전환 자체가 조립 단계입니다 — 이 결과는 **기대 산물**이며, "
        "실제로는 효모에서 플라스미드를 회수해 확인해야 합니다."
    )

    def _check_plan(self, plan) -> None:
        """접합부 교차를 **차단**한다 (Gibson은 경고까지만).

        in vivo라서 잘못된 재조합이 콜로니를 키우고 플라스미드를 회수한 **뒤에야**
        드러난다. 확인 시점이 뒤에 있으므로 설계 단계에서 막는다.

        Args:
            plan: `workflow.junction_plan.JunctionPlan`

        Raises:
            RuntimeError: 비인접 fragment가 임계 이상으로 서열을 공유할 때.
        """
        if not plan.crosstalk_pairs:
            return
        from workflow.junction_plan import JUNCTION_UNIQUE_MAX_BP

        detail = "\n".join(
            f"    · {a} ↔ {b}: {n}bp 공유 ({orient})"
            for a, b, n, orient in plan.crosstalk_pairs[:6]
        )
        raise RuntimeError(
            f"yeast HR 설계 중단 — 인접하지 않은 fragment가 서열을 공유합니다"
            f"(임계 {JUNCTION_UNIQUE_MAX_BP}bp, 최대 {plan.max_crosstalk}bp):\n{detail}\n"
            f"  효모 세포 안에서는 의도한 접합부 말고 **이 자리끼리도** 재조합할 수 있고, "
            f"형질전환이 곧 조립이라 잘못된 산물이 실험 뒤에야 드러납니다.\n"
            f"  → 같은 부품을 두 자리에 쓰고 있다면 **하나를 다른 부품으로 바꾸세요.** "
            f"접합부에 중립 서열을 끼우는 길은 없습니다 — RBS·CDS 경계처럼 간격이 기능을 "
            f"정하는 자리에 삽입되면 조립은 성공해도 발현이 죽기 때문입니다."
        )

    def assemble(self, fragments: List[Dseqrecord],
                 assembly_type: str = "circular") -> List[Dseqrecord]:
        """조립 시뮬레이션. 상동성 기반이라 Gibson과 같은 연산을 쓴다.

        Args:
            fragments: 조립할 fragment 리스트
            assembly_type: `"circular"`(기본 — 효모에서 회수할 플라스미드) 또는 `"linear"`

        Returns:
            조립 결과 리스트. pydna는 같은 원형 분자를 **역상보 회전으로 2개** 돌려주므로,
            유일성을 확인할 때는 회전·역상보로 정규화한 뒤 비교해야 한다.
        """
        return super().assemble(fragments, assembly_type=assembly_type)

    def validate(self, result: Dseqrecord, vector: Optional[Dseqrecord] = None) -> Dict:
        """결과 검증 — Gibson 검증 + in vivo 전제 + 벡터 적격성.

        Args:
            result: 조립 결과물
            vector: 삽입에 쓴 벡터. 주면 효모 복제원점·선택마커를 함께 점검한다
                (`check_yeast_vector`). 생략하면 그 점검은 건너뛴다.

        Returns:
            Gibson `validate()` 결과에 `warnings`가 추가되고, 벡터를 받았으면
            `details['yeast_vector']`에 판정 내역이 담긴다.
        """
        v = super().validate(result)
        v.setdefault("warnings", []).append(self.IN_VIVO_NOTE)

        if vector is not None:
            check = check_yeast_vector(vector)
            v.setdefault("details", {})["yeast_vector"] = check
            if check["warning"]:
                v["warnings"].append(check["warning"])
        return v


__all__ = [
    "YeastHRCloning",
    "check_yeast_vector",
    "YEAST_ORI_KEYWORDS",
    "YEAST_MARKER_KEYWORDS",
    "YEAST_VECTOR_HINTS",
]
