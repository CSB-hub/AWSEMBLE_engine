"""중간 클로닝용 module GG 벡터 생성 (circular 조립 자산).

## 왜 필요한가

circular(중간 클로닝) 모드에서는 **spacer가 부품이 아니라 벡터에 박혀 있다.**
그래서 module 위치마다 (S_i, S_i+1) 쌍을 가진 전용 벡터가 필요한데, 그런 벡터가
아직 DB에 없다. 이 스크립트가 그 벡터를 설계·생성한다.

## 두 가지 backbone (2026-09-23)

1. **실사용 벡터에 삽입** (권장, `--backbone <벡터ID>`): DB의 실제 벡터를 받아
   그 **유일 절단 부위(= MCS로 추정)** 에 카세트를 끼운다. ori·마커가 실물이라
   습식에서 바로 쓸 수 있고, 산출물이 "이 벡터를 이 효소로 잘라 이걸 넣어라"라는
   실행 가능한 설계가 된다.
2. **합성 최소 backbone** (기본, backbone 미지정): 아래 `_BACKBONE_5/_3`.
   ori·마커가 없어 **습식에서 그대로 못 쓴다.** 구조 검토·테스트용이다.

## 구조

    [backbone]─[excision site]─[L spacer]─[asm site→][stuffer][←asm site]─[R spacer]─[excision site]─[backbone]

- **조립 효소(asm)** 사이트는 stuffer 쪽을 향한다 → 절단 시 인식서열이 stuffer와 함께 빠지고
  spacer 안쪽 끝에 overhang이 노출된다. 여기에 부품들이 들어간다.
- **excision 효소** 사이트는 **spacer 바깥**에 둔다. 안쪽에 두면 module을 꺼낼 때 spacer가
  벡터에 남아, 다음 단계(module끼리 Gibson 연결)의 homology가 사라진다.
- 조립 효소와 excision 효소는 반드시 달라야 한다. Type IIS 조립은 인식서열을 없애므로
  같은 효소로는 되꺼낼 수 없다.

## 표준 spacer 세트

경계를 순서대로 고정하면(S1→S2→S3…) 필요한 벡터는 **module 위치 개수만큼**이고
설계가 뭐든 재사용된다. palindromic overhang 없이 "경우의 수" 문제가 풀린다.
후보 선정 근거는 develop_log/ljw/260903_design_parameterization_plan.md 참고.

사용:
    python tools_dev/make_intermediate_vector.py --positions 3 --out-dir ./vectors
"""
import os
import sys
from typing import List, Optional

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from pydna.dseqrecord import Dseqrecord

from basic.sequence_analysis import get_typeIIS_params

# module 6개(경계 7개)까지 커버하는 표준 spacer 세트.
# 양 끝 20bp Tm이 58~64℃로 뭉치고(폭 6.2℃), 상호 최장공통 8bp(오조립 기준 15bp 미만).
STANDARD_SPACERS = ["S0092", "S0023", "S0063", "S0086", "S0032", "S0113", "S0029"]

# 최소 기능 backbone (ori/마커는 실제 설계 시 습식 쪽에서 확정)
_BACKBONE_5 = "TTTTTGCATCTAGACCCGGGATCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCC"
_BACKBONE_3 = "CCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC"
_STUFFER = "GGGGTACCCCTTTTGGGGAAAACCCCTTTTGGGGTACCCCAAAATTTTGGGGCCCCAAAA"


#: MCS 추정에 쓰는 범용 클로닝 효소. 벡터에서 **정확히 1회** 자르는 효소는
#: 관례상 MCS에 있다(= 그 자리에 넣으라고 만들어 둔 자리).
COMMON_CLONING_ENZYMES = [
    "EcoRI", "BamHI", "HindIII", "XbaI", "SalI", "PstI", "SphI",
    "KpnI", "NotI", "SpeI", "NheI", "SmaI", "ClaI", "XhoI",
]


def _rc(s: str) -> str:
    return str(Seq(s).reverse_complement()) if s else ""


def _site_positions(seq: str, enzyme: str) -> List[int]:
    """원형 서열에서 효소 절단 위치(원형 좌표, 중복 제거).

    ⚠️ 선형 검색만 하면 **origin을 걸치는 사이트를 놓친다**(V0001 실측 전례).
    2배 서열에서 찾고 원형 좌표로 접는다.
    """
    from Bio import Restriction
    s = str(seq).upper()
    if not s:
        return []
    hits = getattr(Restriction, enzyme).search(Seq(s + s))
    return sorted({h % len(s) for h in hits if h <= len(s)})


def check_backbone_free(backbone: Dseqrecord, enzymes: List[str]) -> None:
    """backbone에 조립·excision 효소 사이트가 **없는지** 확인. 있으면 에러.

    있으면 module을 꺼낼 때(또는 조립할 때) **backbone이 같이 잘린다.**
    DB 실측(2026-09-23): 벡터 387개 중 BsaI 없는 것 157개(41%), BsmBI 없는 것 131개(34%),
    둘 다 없는 것 50개.
    """
    bad = {e: _site_positions(backbone.seq, e) for e in enzymes}
    bad = {e: v for e, v in bad.items() if v}
    if bad:
        raise ValueError(
            f"backbone '{backbone.name}'에 사용할 효소 사이트가 있습니다: "
            + ", ".join(f"{e} {len(v)}곳{v}" for e, v in bad.items())
            + ". 이 벡터로는 module을 꺼낼 때 backbone도 잘립니다 — 다른 벡터를 고르거나 "
              "다른 효소 조합을 쓰세요(후보는 scan_db_candidates로 확인)."
        )


def find_unique_cloning_sites(backbone: Dseqrecord,
                              enzymes: Optional[List[str]] = None) -> List[tuple]:
    """벡터에서 **정확히 1회** 자르는 범용 효소 목록 → [(효소명, 위치), ...].

    유일 사이트는 관례상 MCS다. 여러 개면 호출자가 고른다(앞쪽이 더 흔한 선택).
    """
    out = []
    for e in (enzymes or COMMON_CLONING_ENZYMES):
        pos = _site_positions(backbone.seq, e)
        if len(pos) == 1:
            out.append((e, pos[0]))
    return out


def scan_db_candidates(db, asm_enzyme: str = "BsaI",
                       excision_enzyme: str = "BsmBI") -> List[dict]:
    """DB `vector`에서 중간 벡터로 쓸 수 있는 후보를 찾는다 (조회만).

    조건: 조립 효소·excision 효소 사이트가 **둘 다 없고**, 유일 절단 부위가 있을 것.
    """
    if not db._connect("Bioparts"):
        raise RuntimeError("Bioparts DB 접속 실패")
    try:
        cur = db.conn.cursor()
        cur.execute("SELECT ID, Name, Sequence FROM `vector` ORDER BY ID")
        rows = cur.fetchall()
        cur.close()
    finally:
        db._disconnect()

    out = []
    for vid, vname, vseq in rows:
        if not vseq:
            continue
        if _site_positions(vseq, asm_enzyme) or _site_positions(vseq, excision_enzyme):
            continue
        rec = Dseqrecord(str(vseq).upper(), name=vid, circular=True)
        uniq = find_unique_cloning_sites(rec)
        if uniq:
            out.append({"id": vid, "name": vname or vid, "len": len(vseq),
                        "sites": [e for e, _ in uniq]})
    return out


def build_intermediate_vector(
    l_spacer_seq: str,
    r_spacer_seq: str,
    asm_enzyme: str = "BsaI",
    excision_enzyme: str = "BsmBI",
    name: str = "pMOD",
    stuffer: str = _STUFFER,
    backbone: Optional[Dseqrecord] = None,
    insert_enzyme: Optional[str] = None,
) -> Dseqrecord:
    """중간 클로닝용 벡터 한 개 생성.

    Args:
        l_spacer_seq / r_spacer_seq: 벡터에 박을 좌/우 spacer 서열
        asm_enzyme: 부품 조립에 쓸 Type IIS 효소
        excision_enzyme: module을 꺼낼 효소 (조립 효소와 달라야 함)
        name: 벡터 이름
        stuffer: 조립 전 자리를 채우는 더미 (절단 시 인식서열과 함께 빠짐)
        backbone: 실사용 벡터(Dseqrecord). 주면 그 벡터의 **유일 절단 부위**에 카세트를
            끼운다. 생략하면 합성 최소 backbone을 쓴다(습식 사용 불가).
        insert_enzyme: 삽입 자리를 정할 효소. 생략 시 유일 절단 부위 중 첫 번째.

    Returns:
        circular Dseqrecord (feature 주석 포함). `annotations["cloning_strategy"]`에
        실사용 backbone일 때의 클로닝 방법이 문자열로 담긴다.

    Raises:
        ValueError: backbone에 조립/excision 효소 사이트가 있거나, 유일 절단 부위가 없거나,
            삽입 효소의 사이트가 카세트 안에도 있을 때.
    """
    if asm_enzyme == excision_enzyme:
        raise ValueError(
            f"조립 효소와 excision 효소가 같습니다({asm_enzyme}). "
            f"Type IIS 조립은 인식서열을 없애므로 같은 효소로는 module을 꺼낼 수 없습니다."
        )
    a = get_typeIIS_params(asm_enzyme)
    e = get_typeIIS_params(excision_enzyme)
    a_site, a_sp = a["fwd"], a["def_spacer"]
    e_site, e_sp = e["fwd"], e["def_spacer"]

    # 조립 효소는 stuffer를 향하게 배치한다(절단 시 인식서열이 stuffer와 함께 제거됨).
    #   ... [L spacer] [RC(spacer더미)][RC(asm site)] [stuffer] [asm site][spacer더미] [R spacer] ...
    left_arm = l_spacer_seq + _rc(a_sp) + _rc(a_site)
    right_arm = a_site + a_sp + r_spacer_seq

    # excision 효소는 spacer **바깥**에서 안쪽(spacer)을 향하게 배치한다.
    #   → 절단하면 [L spacer][부품들][R spacer]가 통째로 빠져 나온다.
    exc_left = e_site + e_sp
    exc_right = _rc(e_sp) + _rc(e_site)

    # 카세트 = 실제로 벡터에 새로 들어가는 조각
    cassette = exc_left + left_arm + stuffer + right_arm + exc_right

    strategy = None
    if backbone is None:
        # 합성 최소 backbone (ori/마커 없음 — 구조 검토용)
        seq = _BACKBONE_5 + cassette + _BACKBONE_3
        offset = len(_BACKBONE_5)
    else:
        # ── 실사용 벡터에 삽입 ──
        check_backbone_free(backbone, [asm_enzyme, excision_enzyme])
        uniq = find_unique_cloning_sites(backbone)
        if not uniq:
            raise ValueError(
                f"backbone '{backbone.name}'에 유일 절단 부위가 없어 삽입 위치를 정할 수 없습니다. "
                f"(검사한 효소: {COMMON_CLONING_ENZYMES})"
            )
        table = dict(uniq)
        if insert_enzyme is None:
            insert_enzyme, pos = uniq[0]
        elif insert_enzyme in table:
            pos = table[insert_enzyme]
        else:
            raise ValueError(
                f"{insert_enzyme}는 '{backbone.name}'에서 유일하게 자르지 않습니다 "
                f"(유일 절단 효소: {[e for e, _ in uniq]})."
            )
        # 삽입 효소가 카세트 안에도 있으면 삽입물이 잘린다.
        if _site_positions(cassette, insert_enzyme):
            raise ValueError(
                f"{insert_enzyme} 사이트가 삽입 카세트 안에도 있습니다 — 잘라 넣으면 "
                f"카세트가 함께 잘립니다. 다른 삽입 효소를 고르세요."
            )
        bb = str(backbone.seq).upper()
        seq = bb[:pos] + cassette + bb[pos:]
        offset = pos
        strategy = (f"{backbone.name}을(를) {insert_enzyme}로 절단(유일 사이트, {pos}번 위치) 후 "
                    f"{len(cassette)}bp 카세트를 삽입")

    vec = Dseqrecord(seq, name=name, circular=True)
    if backbone is not None:
        # backbone에 feature가 있으면 삽입 지점 뒤쪽을 카세트 길이만큼 민다.
        for f in backbone.features:
            try:
                st, en = int(f.location.start), int(f.location.end)
            except Exception:
                continue
            shift = len(cassette) if st >= offset else 0
            vec.features.append(
                SeqFeature(FeatureLocation(st + shift, en + shift, strand=f.location.strand),
                           type=f.type, qualifiers=dict(f.qualifiers))
            )
        vec.annotations["cloning_strategy"] = strategy

    # feature 주석 (설계 검토용)
    def add(start, end, ftype, label):
        vec.features.append(
            SeqFeature(FeatureLocation(start, end, strand=1), type=ftype,
                       qualifiers={"label": [label]})
        )
    i = offset
    add(i, i + len(exc_left), "misc_feature", f"{excision_enzyme} site (excision)")
    i += len(exc_left)
    add(i, i + len(l_spacer_seq), "misc_feature", "L spacer")
    i += len(l_spacer_seq)
    add(i, i + len(_rc(a_sp) + _rc(a_site)), "misc_feature", f"{asm_enzyme} site (assembly)")
    i += len(_rc(a_sp) + _rc(a_site))
    add(i, i + len(stuffer), "misc_feature", "stuffer (조립 시 제거)")
    i += len(stuffer)
    add(i, i + len(a_site + a_sp), "misc_feature", f"{asm_enzyme} site (assembly)")
    i += len(a_site + a_sp)
    add(i, i + len(r_spacer_seq), "misc_feature", "R spacer")
    i += len(r_spacer_seq)
    add(i, i + len(exc_right), "misc_feature", f"{excision_enzyme} site (excision)")
    return vec


def build_vector_set(
    spacer_seqs: List[str],
    asm_enzyme: str = "BsaI",
    excision_enzyme: str = "BsmBI",
    prefix: str = "pMOD",
    backbone: Optional[Dseqrecord] = None,
    insert_enzyme: Optional[str] = None,
) -> List[Dseqrecord]:
    """경계 spacer 목록 → 위치별 벡터 세트.

    spacer가 n+1개면 module 위치 n개용 벡터 n개를 만든다:
    (S1,S2) (S2,S3) ... 순으로 인접 module이 spacer를 공유한다.
    """
    if len(spacer_seqs) < 2:
        raise ValueError("경계 spacer가 2개 이상 필요합니다(module 1개 = 경계 2개).")
    return [
        build_intermediate_vector(
            spacer_seqs[i], spacer_seqs[i + 1], asm_enzyme, excision_enzyme,
            name=f"{prefix}{i + 1}",
            backbone=backbone, insert_enzyme=insert_enzyme,
        )
        for i in range(len(spacer_seqs) - 1)
    ]


def fetch_standard_spacers(db, ids: Optional[List[str]] = None) -> List[str]:
    """DB에서 표준 spacer 서열을 순서대로 가져온다."""
    ids = ids or STANDARD_SPACERS
    out = []
    for sid in ids:
        rec = db.get_part_by_id("Bioparts", "spacer", sid)
        if rec is None:
            raise ValueError(f"spacer를 찾을 수 없습니다: {sid}")
        out.append(str(rec.seq).upper())
    return out


def main():
    import argparse
    from utils.io import DatabaseImport, DNAFileWriter

    ap = argparse.ArgumentParser(description="중간 클로닝용 module GG 벡터 생성")
    ap.add_argument("--positions", type=int, default=3, help="module 위치 개수(벡터 개수)")
    ap.add_argument("--asm-enzyme", default="BsaI")
    ap.add_argument("--excision-enzyme", default="BsmBI")
    ap.add_argument("--out-dir", default="./intermediate_vectors")
    ap.add_argument("--backbone", default=None,
                    help="실사용 벡터 ID(DB `vector`). 생략하면 합성 최소 backbone(습식 사용 불가)")
    ap.add_argument("--insert-enzyme", default=None,
                    help="삽입 자리를 정할 효소. 생략 시 유일 절단 부위 중 첫 번째")
    ap.add_argument("--list-candidates", action="store_true",
                    help="중간 벡터로 쓸 수 있는 DB 벡터 후보만 출력하고 종료")
    args = ap.parse_args()

    db = DatabaseImport()

    if args.list_candidates:
        cands = scan_db_candidates(db, args.asm_enzyme, args.excision_enzyme)
        print(f"조립 {args.asm_enzyme} / excision {args.excision_enzyme} 기준 후보 {len(cands)}개")
        print(f"{'ID':<16}{'Name':<14}{'bp':>8}  유일 절단 부위")
        for c in cands:
            print(f"{c['id']:<16}{c['name'][:13]:<14}{c['len']:>8}  {', '.join(c['sites'][:6])}")
        return

    backbone = None
    if args.backbone:
        backbone = db.get_part_by_id("Bioparts", "vector", args.backbone)
        if backbone is None:
            raise SystemExit(f"벡터를 찾을 수 없습니다: {args.backbone}")
        backbone.name = args.backbone
        backbone = Dseqrecord(str(backbone.seq).upper(), name=args.backbone, circular=True)
        print(f"backbone: {args.backbone} ({len(backbone.seq)}bp)")
    else:
        print("⚠️ backbone 미지정 — 합성 최소 backbone을 씁니다(ori·마커 없음, 습식 사용 불가).")
        print("   실사용 벡터를 쓰려면 --backbone <ID> (후보는 --list-candidates)")
    need = args.positions + 1
    if need > len(STANDARD_SPACERS):
        raise SystemExit(
            f"표준 spacer 세트는 {len(STANDARD_SPACERS)}개(module {len(STANDARD_SPACERS)-1}개)까지입니다."
        )
    seqs = fetch_standard_spacers(db, STANDARD_SPACERS[:need])
    vectors = build_vector_set(seqs, args.asm_enzyme, args.excision_enzyme,
                               backbone=backbone, insert_enzyme=args.insert_enzyme)

    os.makedirs(args.out_dir, exist_ok=True)
    for v in vectors:
        path = os.path.join(args.out_dir, f"{v.name}.gb")
        DNAFileWriter.write_genbank(v, path)
        strat = v.annotations.get("cloning_strategy")
        print(f"  {v.name}: {len(v.seq)}bp → {path}")
        if strat:
            print(f"      ↳ {strat}")
    print(f"\n{len(vectors)}개 생성 (조립 {args.asm_enzyme} / excision {args.excision_enzyme})")


if __name__ == "__main__":
    main()
