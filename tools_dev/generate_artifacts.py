"""전수 조립 산출물 생성 — 실험자 산출물(gbk + primers.csv)을 코퍼스 전체에 대해 디스크 저장.

self-consistency 하네스는 TemporaryDirectory에서 검증만 하고 산출물을 버린다. 이 스크립트는
같은 프로덕션 경로(GoldenGateModule.prepare_fragments → _assemble_fragments →
write_primer_csv)를 돌려 **실제 gbk/primers.csv를 test_id별 디렉토리로 영구 저장**한다.

산출물 레이아웃 (계획 §1 item2/3/6):
    {ROOT}/{test_id}/
        parts/gg_prepared/{part_name}.gb   ← item2 prepared part (full_sequence)
        modules/assembled/{test_id}.gb     ← item3 module
        reports/primers.csv                ← item6 프라이머+조건
    {ROOT}/_index.csv                      ← test_id별 상태(ok/skip 사유·부품수·module 길이·효소)

조립 불가(제3효소 충돌/데이터누락/크래시)는 skip하고 _index.csv에 사유 기록.

실행:
    PYTHONPATH=src /home/csb/awsemble_env/bin/python tools_dev/generate_artifacts.py --all --out <dir>
    옵션: --inst igem|KRIBB  --sample N  --seed N
"""

import argparse
import contextlib
import csv
import io
import sys
import traceback
from pathlib import Path
from typing import List, Optional

from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileWriter
from workflow.biopart_assemble import BiopartWorkflow

BIOPARTS_SCHEMA = "Bioparts"


def _lookup_part(db, ptype, pid, is_kribb):
    """KRIBB는 Bioparts에 K 접두사로 적재됨(KP0019 등) → K+pid 우선, raw 폴백."""
    if is_kribb:
        return db.get_part_by_id(BIOPARTS_SCHEMA, ptype, "K" + pid) \
            or db.get_part_by_id(BIOPARTS_SCHEMA, ptype, pid)
    return db.get_part_by_id(BIOPARTS_SCHEMA, ptype, pid)


def _classify(msg: str) -> str:
    if "제한효소 충돌" in msg or "제3의 Type IIS" in msg:
        return "skip:enzyme_conflict"
    if "ATGC 외 염기" in msg or "빈 서열" in msg:
        return "skip:invalid_base"
    return "skip:error"


def generate_one(tid: str, db, wf: BiopartWorkflow, out_root: Path, is_kribb: bool) -> dict:
    """test_id 1건의 산출물 생성 → 상태 dict."""
    from assembly.golden_gate import GoldenGateModule

    ordered = db.get_ordered_parts_by_test_id(tid)
    if not ordered:
        return {"test_id": tid, "status": "skip:db_missing", "n_parts": 0, "module_len": "", "enzyme": ""}
    parts = []
    for pid, ptype in ordered:
        p = _lookup_part(db, ptype, pid, is_kribb)
        if p is None:
            return {"test_id": tid, "status": "skip:db_missing", "n_parts": len(ordered),
                    "module_len": "", "enzyme": ""}
        parts.append(p)

    try:
        gg = GoldenGateModule()
        with contextlib.redirect_stdout(io.StringIO()):
            design = gg.prepare_fragments(parts, assembly_type="linear")
    except Exception as e:
        return {"test_id": tid, "status": _classify(str(e)), "n_parts": len(parts),
                "module_len": "", "enzyme": ""}

    records = []
    for d in design:
        r = Dseqrecord(d["full_sequence"], name=f"{d['part_name']}_designed")
        r.annotations["target_seq"] = d.get("target_seq", "")
        r.annotations["restriction_enzyme"] = d.get("restriction_enzyme", "BsaI")
        r.annotations["part_id"] = d.get("part_id", "")
        r.annotations["part_type"] = d.get("part_type", "")
        records.append(r)

    with contextlib.redirect_stdout(io.StringIO()):
        module = wf._assemble_fragments(records, method="golden_gate", assembly_type="linear")
    if module is None:
        return {"test_id": tid, "status": "skip:assemble_fail", "n_parts": len(parts),
                "module_len": "", "enzyme": design[0].get("restriction_enzyme", "")}

    case = out_root / tid
    pp = case / "parts" / "gg_prepared"
    md = case / "modules" / "assembled"
    rp = case / "reports"
    for d_ in (pp, md, rp):
        d_.mkdir(parents=True, exist_ok=True)

    with contextlib.redirect_stdout(io.StringIO()):
        # item2: prepared part gbk
        for r in records:
            DNAFileWriter.write_genbank(r, str(pp / f"{r.name}.gb"))
        # item3: module gbk
        DNAFileWriter.write_genbank(module, str(md / f"{tid}.gb"))
        # item6: primers.csv (프로덕션 write_primer_csv 재사용)
        wf.project_dir = case
        wf.write_primer_csv(design)

    return {"test_id": tid, "status": "ok", "n_parts": len(parts),
            "module_len": len(module.seq), "enzyme": design[0].get("restriction_enzyme", "")}


def main():
    ap = argparse.ArgumentParser(description="전수 조립 산출물(gbk+primers.csv) 생성·저장")
    ap.add_argument("--out", required=True, help="산출물 루트 디렉토리")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--inst", choices=["igem", "KRIBB"])
    ap.add_argument("--sample", type=int)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    db = DatabaseImport()
    if args.inst:
        ids = db.get_test_ids_by_institution(args.inst)
    else:
        ids = db.get_all_test_ids_from_construct()
    if args.sample:
        import random
        random.seed(args.seed)
        ids = random.sample(ids, min(args.sample, len(ids)))

    out_root = Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)
    kribb_set = set(db.get_test_ids_by_institution("KRIBB"))

    # scratch 워크플로우 1개(재사용) — project_dir는 case별로 repoint(표준레이아웃 중복생성 회피)
    wf = BiopartWorkflow(project_dir=str(out_root / "_scratch"), parts=[])

    print(f"🏭 산출물 생성: {len(ids)}개 → {out_root}")
    from collections import Counter
    stat = Counter()
    rows = []
    for k, tid in enumerate(ids, 1):
        try:
            r = generate_one(tid, db, wf, out_root, tid in kribb_set)
        except Exception:
            r = {"test_id": tid, "status": "skip:crash", "n_parts": "", "module_len": "", "enzyme": ""}
            print(f"  ⚠️ {tid} 크래시:\n{traceback.format_exc(limit=2)}")
        rows.append(r)
        stat[r["status"].split(":")[0]] += 1
        if k % 500 == 0:
            print(f"  [{k}/{len(ids)}] ok={stat['ok']} skip={sum(v for s,v in stat.items() if s!='ok')}")

    with open(out_root / "_index.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["test_id", "status", "n_parts", "module_len", "enzyme"])
        w.writeheader()
        w.writerows(rows)

    print(f"\n✅ 완료: ok {stat['ok']} / skip {sum(v for s,v in stat.items() if s!='ok')}")
    print("상태 분해:", dict(Counter(r["status"] for r in rows)))
    print(f"📄 인덱스: {out_root/'_index.csv'}")


if __name__ == "__main__":
    main()
