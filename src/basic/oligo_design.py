"""짧은 부품의 **올리고 어닐링** 제작안 — PCR 대신 올리고 2개를 주문해 붙인다.

## 왜 PCR이 아닌가

프로모터·RBS·터미네이터는 서열이 짧다(실측: Cello 파츠 기준 RBS 9~35bp,
터미네이터 34~60bp, 프로모터 63~80bp). PCR은 fwd/rev 프라이머가 겹치지 않고
앉을 자리가 있어야 하는데 이 길이에선 그게 안 나온다. 억지로 설계하면
**쓸 수 없는 프라이머**가 나온다 — 9bp 부품은 Tm 33.5~36.2℃로 어떤 프리셋 창에도
못 드는데 `passed=True`로 통과하고 있었다(2026-10-12 실측).

표준 올리고 합성이 100 base까지 되므로, 그 아래는 올리고 두 가닥을 주문해
어닐링하는 쪽이 비용·확실성 모두 낫다. CDS처럼 긴 부품만 PCR로 간다.

## 올리고 모양 (sticky end 미리 형성)

Golden Gate 사이트를 넣어 같이 자르는 대신, **overhang을 미리 만들어** 주문한다.
그러면 그 부품은 절단 단계를 건너뛰고 바로 ligation에 참여한다.

    Top:     5'-[좌 overhang][----- 부품 -----]          -3'
    Bottom:  3'-              [----- 부품' ----][우 overhang']-5'

즉 두 올리고 모두 길이가 `len(부품) + 4`이고, 어닐링하면 양 끝에 4nt 단일가닥이
남아 이웃 조각과 결합한다.

## 안정성 등급

어닐링한 두 가닥을 붙잡는 건 **부품 본체 길이뿐**이다 — overhang은 단일가닥이고
(그게 sticky end의 정의다), 패딩·인식서열은 애초에 붙이지 않는다.
one-pot GG는 37℃까지 올라가므로 듀플렉스가 그 온도에서 버텨야 한다.
"""
from typing import Dict, List, Optional

from Bio.Seq import Seq

#: 등급 이름 — 숫자가 아니라 이 문자열이 시트·UI에 그대로 나간다.
GRADE_SAFE = "안전"
GRADE_CAUTION = "주의"
GRADE_RISKY = "위험"


def _rc(seq: str) -> str:
    return str(Seq(seq.upper()).reverse_complement())


def grade_duplex(part_len: int, params=None) -> Dict[str, str]:
    """부품 길이로 어닐링 듀플렉스 안정성 등급을 매긴다.

    차단하지 않는다 — 등급과 이유만 돌려주고, 쓸지 말지는 실험자가 정한다.
    """
    from params import DEFAULT_PARAMS
    p = params if params is not None else DEFAULT_PARAMS
    if part_len >= p.duplex_safe_bp:
        return {"grade": GRADE_SAFE, "reason": ""}
    if part_len >= p.duplex_risky_bp:
        return {
            "grade": GRADE_CAUTION,
            "reason": (f"{part_len}bp — 듀플렉스 Tm이 대략 41~50℃로, one-pot GG의 "
                       f"37℃ 절단 단계 대비 여유가 적습니다."),
        }
    return {
        "grade": GRADE_RISKY,
        "reason": (f"{part_len}bp — 듀플렉스 Tm이 37℃ 안팎이거나 그 아래라 "
                   f"절단 단계에서 두 가닥이 풀릴 수 있습니다. 이 요소는 단독 부품으로 두지 말고 "
                   f"이웃 부품에 합쳐 넣는 것을 권합니다."),
    }


def design_annealed_oligos(
    part_seq: str,
    left_oh: str = "",
    right_oh: str = "",
    name: str = "",
    params=None,
) -> Dict[str, object]:
    """부품 하나를 어닐링 올리고 한 쌍으로 만든다.

    Args:
        part_seq: 부품 본체 서열(overhang 제외)
        left_oh / right_oh: 이 부품 양 옆의 4nt overhang. 모듈 바깥 끝은 빈 문자열.
        name: 부품 이름(시트 표기용)

    Returns:
        top / bottom (주문할 두 서열), 길이, 등급, 주의사항
    """
    part = part_seq.upper()
    lo = (left_oh or "").upper()
    ro = (right_oh or "").upper()

    top = lo + part
    bottom = _rc(ro) + _rc(part)

    g = grade_duplex(len(part), params)
    return {
        "name": name,
        "part_len": len(part),
        "top": top,
        "bottom": bottom,
        "top_len": len(top),
        "bottom_len": len(bottom),
        "left_oh": lo,
        "right_oh": ro,
        "grade": g["grade"],
        "note": g["reason"],
    }


def should_use_oligo(part_len: int, params=None) -> bool:
    """이 길이면 PCR 대신 올리고 어닐링을 권하는가."""
    from params import DEFAULT_PARAMS
    p = params if params is not None else DEFAULT_PARAMS
    return part_len < p.oligo_max_bp


def production_note(part_len: int, params=None) -> str:
    """PCR 시트에 남길 "왜 올리고를 권하는가" 한 줄.

    길이에 따라 **이유가 다르다** — 뭉뚱그려 "증폭이 어렵다"고 하면 거짓이 된다.
    80bp 부품은 PCR로 잘 되고(실측 Tm이 창 안에 든다), 올리고를 권하는 이유는
    비용·확실성이지 불가능해서가 아니다.
    """
    from params import DEFAULT_PARAMS
    p = params if params is not None else DEFAULT_PARAMS
    if part_len < p.duplex_safe_bp:
        return (f"{part_len}bp — PCR 프라이머 2개가 앉을 자리가 없어 증폭이 사실상 불가능합니다. "
                f"다만 이 길이는 어닐링 듀플렉스도 불안정하니 등급을 확인하세요.")
    if part_len < 2 * p.primer_len_opt:
        return (f"{part_len}bp — fwd/rev 프라이머가 겹쳐 품질이 떨어집니다. 올리고 어닐링을 권합니다.")
    return (f"{part_len}bp — PCR도 되지만 이 길이는 올리고 어닐링이 더 싸고 확실합니다.")


def build_oligo_rows(designs: List[dict], params=None) -> List[Dict[str, object]]:
    """prepare_fragments 결과에서 **올리고로 만들 부품**만 골라 주문 행을 만든다.

    designs: `GoldenGateModule.prepare_fragments`가 돌려준 설계 목록
             (target_seq / left_oh / right_oh / part_name 사용)
    """
    rows: List[Dict[str, object]] = []
    for d in designs:
        part = str(d.get("target_seq") or "")
        if not part or not should_use_oligo(len(part), params):
            continue
        o = design_annealed_oligos(
            part,
            left_oh=d.get("left_oh") or "",
            right_oh=d.get("right_oh") or "",
            name=str(d.get("part_name") or d.get("part_id") or ""),
            params=params,
        )
        rows.append({
            "부품": o["name"],
            "부품 길이(bp)": o["part_len"],
            "등급": o["grade"],
            "주의사항": o["note"],
            "Top 올리고 (5'→3')": o["top"],
            "Top 길이(nt)": o["top_len"],
            "Bottom 올리고 (5'→3')": o["bottom"],
            "Bottom 길이(nt)": o["bottom_len"],
            "좌 overhang": o["left_oh"],
            "우 overhang": o["right_oh"],
        })
    return rows
