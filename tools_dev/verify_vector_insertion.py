"""벡터 삽입 + multi-module 전수 시뮬레이션 + self-consistency 검증 (item4/5).

part→module(verify_selfconsistency.py)을 넘어, **module→vector / multimodule→vector**의
원형 construct 산출물(item5)을 실제로 생성하고 자기모순을 검증한다. 프로덕션 경로
(WorkflowBuilder.build_plasmid / build_multimodule_plasmid)를 그대로 구동한다.

시나리오:
  S1  igem 전체: 각 test parts → 1 module → vector(V0001, V0002) 삽입
  S2  KRIBB 7모듈(경로순) → 1 multimodule → 1 vector
  S3  KRIBB 7모듈 4+3 분할 → 서로 다른 2 vector
  S4  igem 랜덤 2~5모듈 multimodule → 1 vector 또는 2 vector 분할 (다양성)

construct 체크(gold standard 불요):
  circular(원형 닫힘) / parts_present(각 부품 서열 포함) / vector_present(벡터 서열 포함) /
  enzyme_clean(조립효소 BsaI·BbsI 인식서열 잔존 0)

산출물: {out}/{scenario}/{case}/ (build_* 가 constructs/·reports/primers.csv 등 기록) + {out}/_index.csv

실행:
    PYTHONPATH=src /home/csb/awsemble_env/bin/python tools_dev/verify_vector_insertion.py --out <dir> \
        [--scenarios S1 S2 S3 S4] [--s1-sample N] [--s4-sets N] [--seed N]
"""

import argparse
import contextlib
import csv
import io
import random
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

from Bio.Restriction import Restriction
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder

BIOPARTS = "Bioparts"
VECTORS = ["V0001", "V0002"]   # plB2(5550) / pUC19(2686), 둘 다 SmaI 유일 → gibson linearize OK

# 검증(S5)용 blunt-end 벡터 선형화 후보. gibson은 blunt 선형화가 정석 —
# sticky 커터는 homology 재조립 시 인식서열이 서열의존적으로 재생성돼 검증 불변식이 흐려지므로 제외.
# production _linearize_vector도 이미 blunt만 씀(SmaI/EcoRV/PvuII/StuI/DraI) → 이 리스트는 검증 다양성용.
BLUNT_ENZYMES = ["SmaI", "EcoRV", "PvuII", "StuI", "DraI", "ScaI", "HpaI",
                 "NaeI", "SspI", "ZraI", "FspI", "HincII", "EcoICRI"]
# production gibson._linearize_vector의 고정 우선순위(관찰용 미러 — src 미수정)
PROD_PRIORITY = ["SmaI", "EcoRV", "PvuII", "StuI", "DraI"]


def valid_gibson_enzymes(vector) -> List[str]:
    """벡터를 정확히 1번 자르는 blunt-end 후보 효소명 리스트(검증 강제 linearize용)."""
    out = []
    for e in BLUNT_ENZYMES:
        enz = getattr(Restriction, e)
        if enz.is_blunt() and len(enz.search(vector.seq)) == 1:
            out.append(e)
    return out


def production_enzyme(vector) -> Optional[str]:
    """production 우선순위-5가 이 벡터에 고를 효소(관찰). 미러일 뿐 src 호출 아님."""
    return next((e for e in PROD_PRIORITY
                 if len(getattr(Restriction, e).search(vector.seq)) == 1), None)


def linearize_with(vector, enzyme_name):
    """벡터를 지정 효소로 미리 잘라 linear Dseqrecord로 반환(하네스 강제 linearize).
    production _linearize_vector는 linear 입력을 그대로 반환하므로 효소를 강제할 수 있다."""
    lin = vector.cut(getattr(Restriction, enzyme_name))[0]
    lin.name = f"{vector.name}_{enzyme_name}"
    return lin

# KRIBB 메발론산 경로 모듈 순서 (공유 spacer 체인으로 확정: S0001→…→S0009)
KRIBB_ORDER = [
    "VA16-MvaK1-VA3", "VA3-MvaD-VA4", "VA4-MvaK2-VA5", "VA5-Idi-VA6",
    "VA6-MvaE-VA7", "VA7-MvaS-VA8", "VA8-IspA-VA2",
]
PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"


@dataclass
class Check:
    check: str
    status: str
    detail: str = ""


@dataclass
class CaseResult:
    case: str
    scenario: str
    status: str = "ok"          # ok / skip:<reason>
    error: Optional[str] = None
    enzyme: str = ""            # linearize에 사용된 효소(production 관찰 or 강제)
    checks: List[Check] = field(default_factory=list)

    @property
    def failed(self):
        return [c for c in self.checks if c.status == FAIL]


# ------------------------------------------------------------------ #
# DB 헬퍼
# ------------------------------------------------------------------ #
def get_part(db, pid, ptype, is_kribb):
    lid = ("K" + pid) if is_kribb else pid
    return db.get_part_by_id(BIOPARTS, ptype, lid) or \
        (db.get_part_by_id(BIOPARTS, ptype, pid) if is_kribb else None)


def get_vector(db, vid) -> Dseqrecord:
    v = db.get_part_by_id(BIOPARTS, "vector", vid)
    return Dseqrecord(str(v.seq), circular=True, name=vid)


def core_parts(db, test_id, is_kribb) -> Optional[List[Dseqrecord]]:
    """test_id의 core 부품(spacer 제외) Dseqrecord 리스트 (순서 보존). 누락 시 None."""
    ordered = db.get_ordered_parts_by_test_id(test_id)
    parts = []
    for pid, ptype in ordered:
        if ptype == "spacer":
            continue
        p = get_part(db, pid, ptype, is_kribb)
        if p is None:
            return None
        parts.append(p)
    return parts or None


# ------------------------------------------------------------------ #
# construct self-consistency 체크
# ------------------------------------------------------------------ #
def check_construct(construct: Dseqrecord, part_seqs: List[str], vector: Dseqrecord,
                    lin_enzyme: Optional[str], insert: Optional[Dseqrecord]) -> List[Check]:
    checks = []
    cs = str(construct.seq).upper()
    doubled = cs + cs   # circular 회전 대비

    # 1. 원형 닫힘
    checks.append(Check("circular", PASS if construct.circular else FAIL,
                        "" if construct.circular else "construct가 circular 아님"))

    # 2. 각 부품 서열이 construct에 존재(회전 허용)
    miss = [i for i, s in enumerate(part_seqs) if s and s.upper() not in doubled]
    checks.append(Check("parts_present", PASS if not miss else FAIL,
                        "" if not miss else f"부품 미검출 index={miss}/{len(part_seqs)}"))

    # 3. 벡터 서열 존재 (절단으로 회전될 수 있어 서명 슬라이스로 확인)
    vsig = str(vector.seq).upper()
    sig = vsig[:200] if len(vsig) >= 200 else vsig
    v_ok = sig in doubled or str(Seq(sig).reverse_complement()) in doubled
    checks.append(Check("vector_present", PASS if v_ok else FAIL,
                        "" if v_ok else "벡터 서명 서열 미검출"))

    # 4. insert가 linearization site에 실제로 삽입됐는지 — 벡터의 절단부를 '가로지르는' 연속
    #    window(cut 앞뒤 W bp)가 construct에 그대로 남아있지 않아야 한다(insert가 그 자리를 갈랐어야 함).
    #    ⚠️ 인식서열 count 방식은 junction 우연 재생성/degenerate 효소로 흐려지므로, 위치 기반이 견고.
    if lin_enzyme:
        enz = getattr(Restriction, lin_enzyme)
        sites = enz.search(vector.seq)   # 1-based 절단 위치
        if len(sites) == 1:
            vseq = str(vector.seq).upper()
            c = sites[0] - 1             # 0-based
            W = 15
            window = (vseq + vseq)[c + len(vseq) - W: c + len(vseq) + W]   # 원형 wrap 안전
            split = window not in doubled
            checks.append(Check("insertion_at_cut", PASS if split else FAIL,
                                "" if split else
                                f"{lin_enzyme} 절단부 window({len(window)}bp)가 construct에 온존 → 삽입 위치 이상"))
    return checks


def _insert_record(res):
    """build_* 결과에서 vector 삽입 직전 insert(module/multimodule) record 추출.
    run_plan results의 마지막 앞 항목(=construct 직전)이 insert."""
    rs = res.get("results") or []
    return rs[-2]["record"] if len(rs) >= 2 else (rs[0]["record"] if rs else None)


def _classify(msg: str) -> str:
    if "제한효소 충돌" in msg or "제3의 Type IIS" in msg:
        return "skip:enzyme_conflict"
    if "적합" in msg and "spacer" in msg:
        return "skip:spacer_select"
    if "linearization" in msg or "linearize" in msg:
        return "skip:vector_linearize"
    if "ATGC 외" in msg or "빈 서열" in msg:
        return "skip:invalid_base"
    return "skip:error"


# ------------------------------------------------------------------ #
# 실행 래퍼 (build_* 구동 + construct 체크)
# ------------------------------------------------------------------ #
def _prep_vector(vector, lin_enzyme):
    """lin_enzyme 지정 시 그 효소로 미리 잘라 linear 벡터 반환(강제), 없으면 circular 그대로
    (production이 우선순위-5로 선형화). 반환: (build에 넘길 벡터, 실제 사용 효소명, 원본 circular 벡터)."""
    if lin_enzyme:
        return linearize_with(vector, lin_enzyme), lin_enzyme, vector
    return vector, (production_enzyme(vector) or ""), vector


def run_single(db, case, scenario, parts, vector, out_dir, lin_enzyme=None) -> CaseResult:
    r = CaseResult(case, scenario)
    vec_in, r.enzyme, vcirc = _prep_vector(vector, lin_enzyme)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            res = WorkflowBuilder.build_plasmid(parts, vec_in, db, str(out_dir),
                                                module_method="golden_gate", insert_method="gibson")
        con = res.get("construct")
        if con is None:
            r.status = "skip:no_construct"
            return r
        r.checks = check_construct(con, [str(p.seq) for p in parts], vcirc, r.enzyme, _insert_record(res))
    except Exception as e:
        r.status = _classify(str(e))
        r.error = f"{type(e).__name__}: {str(e)[:150]}"
    return r


def run_multi(db, case, scenario, groups, vector, out_dir, lin_enzyme=None) -> CaseResult:
    r = CaseResult(case, scenario)
    vec_in, r.enzyme, vcirc = _prep_vector(vector, lin_enzyme)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            res = WorkflowBuilder.build_multimodule_plasmid(groups, vec_in, db, str(out_dir),
                                                            insert_method="gibson")
        con = res.get("construct")
        if con is None:
            r.status = "skip:no_construct"
            return r
        all_seqs = [str(p.seq) for g in groups for p in g]
        r.checks = check_construct(con, all_seqs, vcirc, r.enzyme, _insert_record(res))
    except Exception as e:
        r.status = _classify(str(e))
        r.error = f"{type(e).__name__}: {str(e)[:150]}"
    return r


# ------------------------------------------------------------------ #
# 시나리오
# ------------------------------------------------------------------ #
def scenario_S1(db, out, ids, results):
    """igem 각 test → 1 module → V0001, V0002 삽입."""
    root = out / "S1_single_module_vector"
    for k, tid in enumerate(ids, 1):
        parts = core_parts(db, tid, is_kribb=False)
        for vid in VECTORS:
            case = f"{tid}__{vid}"
            if parts is None:
                results.append(CaseResult(case, "S1", status="skip:db_missing"))
                continue
            results.append(run_single(db, case, "S1", parts, get_vector(db, vid), root / case))
        if k % 200 == 0:
            print(f"  S1 [{k}/{len(ids)}]")


def scenario_S2(db, out, results):
    """KRIBB 7모듈(경로순) → 1 multimodule → 각 vector."""
    root = out / "S2_kribb_7mod_1vector"
    groups = [core_parts(db, t, is_kribb=True) for t in KRIBB_ORDER]
    if any(g is None for g in groups):
        results.append(CaseResult("KRIBB_7mod", "S2", status="skip:db_missing"))
        return
    for vid in VECTORS:
        case = f"KRIBB_7mod__{vid}"
        results.append(run_multi(db, case, "S2", groups, get_vector(db, vid), root / case))


def scenario_S5(db, out, rep_tid, results):
    """효소 로바스트니스: 각 벡터 × 유효 Gibson 효소로 강제 linearize → 대표 module 삽입 검증.
    gibson 삽입이 SmaI 외 효소로도 동일하게 원형 construct를 만드는지(효소 무관성) 확인."""
    root = out / "S5_enzyme_robustness"
    parts = core_parts(db, rep_tid, is_kribb=False)
    if parts is None:
        results.append(CaseResult(f"{rep_tid}", "S5", status="skip:db_missing"))
        return
    for vid in VECTORS:
        vcirc = get_vector(db, vid)
        for e in valid_gibson_enzymes(vcirc):
            case = f"{rep_tid}__{vid}__{e}"
            results.append(run_single(db, case, "S5", parts, vcirc, root / case, lin_enzyme=e))


def scenario_S3(db, out, results):
    """KRIBB 7모듈 4+3 분할 → 서로 다른 2 vector."""
    root = out / "S3_kribb_split_2vector"
    groups = [core_parts(db, t, is_kribb=True) for t in KRIBB_ORDER]
    if any(g is None for g in groups):
        results.append(CaseResult("KRIBB_split", "S3", status="skip:db_missing"))
        return
    g4, g3 = groups[:4], groups[4:]   # 앞 4개 → V0001, 뒤 3개 → V0002
    results.append(run_multi(db, "KRIBB_first4__V0001", "S3", g4, get_vector(db, "V0001"), root / "first4_V0001"))
    results.append(run_multi(db, "KRIBB_last3__V0002", "S3", g3, get_vector(db, "V0002"), root / "last3_V0002"))


def scenario_S4(db, out, igem_ids, n_sets, seed, results):
    """igem 랜덤 2~5모듈 multimodule → 1 vector 또는 2 vector 분할 (다양성)."""
    root = out / "S4_igem_random_multimodule"
    rng = random.Random(seed)
    # core 부품이 확보되는 igem test만 후보로
    usable = []
    for tid in igem_ids:
        p = core_parts(db, tid, is_kribb=False)
        if p:
            usable.append((tid, p))
        if len(usable) >= 400:   # 후보 풀 충분(랜덤 조합 다양성 확보)
            break

    made = 0
    for i in range(n_sets):
        n_mod = rng.choice([2, 3, 4, 5])
        picks = rng.sample(usable, n_mod)
        groups = [p for _, p in picks]
        tags = "_".join(t.replace("BBa_", "") for t, _ in picks)
        split = (n_mod >= 2 and rng.random() < 0.5)   # 절반은 2벡터 분할
        if not split:
            vid = rng.choice(VECTORS)
            case = f"set{i:03d}_{n_mod}mod_1vec_{vid}_{tags}"[:120]
            results.append(run_multi(db, case, "S4", groups, get_vector(db, vid), root / f"set{i:03d}_1vec"))
        else:
            cut = rng.randint(1, n_mod - 1)
            gA, gB = groups[:cut], groups[cut:]
            vA, vB = rng.sample(VECTORS, 2)   # 서로 다른 2벡터
            # 그룹 크기 1이면 single-module 경로(build_plasmid), 아니면 multimodule
            for tag, g, v in ((f"A{len(gA)}", gA, vA), (f"B{len(gB)}", gB, vB)):
                case = f"set{i:03d}_split_{tag}_{v}"
                d = root / f"set{i:03d}_split" / tag
                if len(g) == 1:
                    results.append(run_single(db, case, "S4", g[0], get_vector(db, v), d))
                else:
                    results.append(run_multi(db, case, "S4", g, get_vector(db, v), d))
        made += 1
        if (i + 1) % 20 == 0:
            print(f"  S4 [{i+1}/{n_sets}]")


# ------------------------------------------------------------------ #
def main():
    ap = argparse.ArgumentParser(description="벡터 삽입 + multimodule 시뮬레이션·검증")
    ap.add_argument("--out", required=True)
    ap.add_argument("--scenarios", nargs="+", default=["S1", "S2", "S3", "S4", "S5"])
    ap.add_argument("--s1-sample", type=int, help="S1 igem 샘플 수(미지정=전수)")
    ap.add_argument("--s4-sets", type=int, default=80)
    ap.add_argument("--s5-rep", default="BBa_K256019", help="S5 대표 module test_id")
    ap.add_argument("--s1-ids-file", help="S1을 이 파일의 test_id(줄단위)로만 실행(실패분 재실행용)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    db = DatabaseImport()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    igem_ids = db.get_test_ids_by_institution("igem")
    if args.s1_ids_file:
        s1_ids = [x.strip() for x in Path(args.s1_ids_file).read_text().splitlines() if x.strip()]
    elif args.s1_sample:
        random.Random(args.seed).shuffle(igem_ids)
        s1_ids = igem_ids[:args.s1_sample]
    else:
        s1_ids = igem_ids

    results: List[CaseResult] = []
    if "S1" in args.scenarios:
        print(f"▶ S1 (igem {len(s1_ids)} × {len(VECTORS)}벡터)")
        scenario_S1(db, out, s1_ids, results)
    if "S2" in args.scenarios:
        print("▶ S2 (KRIBB 7모듈 → 1벡터)")
        scenario_S2(db, out, results)
    if "S3" in args.scenarios:
        print("▶ S3 (KRIBB 4+3 → 2벡터)")
        scenario_S3(db, out, results)
    if "S4" in args.scenarios:
        print(f"▶ S4 (igem 랜덤 multimodule {args.s4_sets}세트)")
        scenario_S4(db, out, db.get_test_ids_by_institution("igem"), args.s4_sets, args.seed, results)
    if "S5" in args.scenarios:
        print("▶ S5 (효소 로바스트니스: 벡터×유효 Gibson 효소 강제 linearize)")
        scenario_S5(db, out, args.s5_rep, results)

    # 집계 + 인덱스
    from collections import Counter
    per_check = {}
    rows = []
    for r in results:
        for c in r.checks:
            per_check.setdefault(c.check, Counter())[c.status] += 1
        rows.append({"scenario": r.scenario, "case": r.case, "status": r.status,
                     "enzyme": r.enzyme, "n_fail": len(r.failed), "error": r.error or "",
                     "fails": "; ".join(f"{c.check}:{c.detail}" for c in r.failed)})

    with open(out / "_index.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["scenario", "case", "status", "enzyme", "n_fail", "error", "fails"])
        w.writeheader(); w.writerows(rows)

    built = [r for r in results if r.status == "ok"]
    clean = [r for r in built if not r.failed]
    print("\n" + "=" * 60)
    print(f"총 {len(results)}케이스 | 조립성공 {len(built)} | 완전통과 {len(clean)} | "
          f"FAIL보유 {len(built)-len(clean)} | skip {len(results)-len(built)}")
    print("\n시나리오별:")
    for s in ["S1", "S2", "S3", "S4", "S5"]:
        rs = [r for r in results if r.scenario == s]
        if rs:
            b = sum(1 for r in rs if r.status == "ok")
            cl = sum(1 for r in rs if r.status == "ok" and not r.failed)
            print(f"  {s}: {len(rs)}케이스 조립성공 {b} 완전통과 {cl}")
    print("\nconstruct 체크별 (PASS/FAIL):")
    for name, cnt in sorted(per_check.items()):
        flag = "  ⬅️ FAIL" if cnt[FAIL] else ""
        print(f"  {name:16s} P{cnt[PASS]:<5d} F{cnt[FAIL]:<4d}{flag}")
    # 사용된 linearize 효소 분포(production 관찰 S1~S4 + 강제 S5)
    enz_dist = Counter(r.enzyme for r in results if r.status == "ok" and r.enzyme)
    print("\nlinearize 효소 분포:", dict(enz_dist))
    print("skip 사유:", dict(Counter(r.status for r in results if r.status != "ok")))
    print(f"\n📄 인덱스: {out/'_index.csv'}")
    sys.exit(1 if any(r.failed for r in built) else 0)


if __name__ == "__main__":
    main()
