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
| arm 교차 비교 범위 | arm끼리만 | **arm ↔ 다른 몸통·백본까지** | `CROSSTALK_SCOPE` |
| 벡터 요건 | 제약 없음 | **효모 복제원점 + 선택마커** | `check_yeast_vector()` |

### arm 교차를 왜 몸통까지 보나

Gibson은 exonuclease가 fragment **말단만** 깎아 단일가닥을 노출시킨다. 그래서 arm이
다른 fragment의 몸통 한복판에 어닐링할 방법이 없고, 비교는 arm끼리로 충분하다.

효모 안에서는 그 제약이 없다 — 상동 구간이면 어디서든 재조합한다. arm이 다른
fragment의 몸통이나 백본의 반복 구간에 내려앉으면 재배열된 산물이 나오는데, 그게
형질전환·플라스미드 회수 **뒤에야** 드러난다. 그래서 비교 범위를 몸통·백본까지
넓힌다(`CROSSTALK_SCOPE = "arm+body"`).

걸렸을 때의 처리는 Gibson과 같다 — `JunctionPlanner._resolve_arm_crosstalk`가 접합부
길이를 다시 고르고(유저 범위 → 하드 범위), 끝까지 안 되면 설계불가로 끊는다.

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
        evidence.append(f"origin candidates: {', '.join(ori)}")
    if marker:
        evidence.append(f"selection marker candidates: {', '.join(marker)}")
    if hints:
        evidence.append(f"name/description hints: {', '.join(hints)}")

    if ori and marker:
        verdict, warning = "ok", ""
    elif ori or marker or hints:
        missing = "a selection marker" if ori else ("an origin of replication" if marker else "an origin of replication and a selection marker")
        verdict = "partial"
        warning = (
            f"Could not confirm {missing} for yeast in vector '{vector.name}' "
            f"({'; '.join(evidence)}). A yeast HR product can only be recovered if it "
            f"replicates and is selected in yeast - check that this is a yeast shuttle vector."
        )
    else:
        verdict = "unknown"
        warning = (
            f"Vector '{vector.name}' has no features or description, so a yeast origin "
            f"(ARS/CEN or 2-micron) and a selection marker **cannot be confirmed** (this does "
            f"not mean they are absent). A yeast HR product is only recovered if it replicates "
            f"and is selected in yeast - check that it is a shuttle vector yourself, or "
            f"register it as a GenBank file with features so this can be judged automatically."
        )

    return {
        "verdict": verdict,
        "ori": ori,
        "marker": marker,
        "hints": hints,
        "evidence": evidence,
        "warning": warning,
    }


from params import _YEAST_HOMOLOGY_RANGE


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
        "In yeast HR the transformation itself is the assembly step - this result is the "
        "**expected product**, and in practice you must recover the plasmid from yeast and verify it."
    )

    #: yeast HR arm의 하드 한계. 문헌이 보고한 범위이자 `params.validate()`가 거는
    #: 하드 리밋과 같은 값이다 — arm이 엉뚱한 자리에 붙어 유저 범위 안에서 못 피할 때
    #: 여기까지만 넘어간다(`_resolve_arm_crosstalk`가 note로 남긴다).
    ARM_HARD_RANGE = _YEAST_HOMOLOGY_RANGE

    #: in vivo라 **상동 구간이면 어디든** 재조합한다 → arm을 다른 arm뿐 아니라 다른
    #: fragment의 몸통·백본과도 비교한다. in vitro Gibson은 말단만 깎여서 몸통에는
    #: 붙을 수가 없으므로 `"arm"`만 본다.
    CROSSTALK_SCOPE = "arm+body"

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
