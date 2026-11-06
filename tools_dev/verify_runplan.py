"""run_plan(build-plan) 경로 동치 스모크 — J-3/K-3 검증.

같은 prepared parts에 대해
  (A) shim 경로: run_step3_module_assembly(module_df)
  (B) build-plan 경로: run_full_workflow(plan=[{op:module_assembly, inputs:[...]}])
두 경로의 module 서열 SHA가 동일한지 확인한다.

build-plan 경로는 명시 module_id를 주지 않으므로 파일명은 결정적 해시(mod_<hash8>)지만,
**서열 자체**는 shim과 동일해야 한다(핵심 조립 로직 불변).

실행:
  docker run --rm --network host -v "$PWD/awsemble/.env:/app/.env:ro" \
    -v "$PWD/awsemble/src:/app/src" -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    awsemble:dev python tools_dev/verify_runplan.py
"""
import sys, hashlib, io, contextlib
from pathlib import Path

import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow as Workflow

SCHEMA = "Bioparts"
IDS = ["BBa_B1101", "BBa_C0630", "BBa_E0240", "BBa_C0480"]   # 2/3/4/7 parts


def _sha(seq: str) -> str:
    return hashlib.sha256(seq.lower().encode()).hexdigest()[:16]


def _design_records(test_id, db):
    ordered = db.get_ordered_parts_by_test_id(test_id)
    parts = [db.get_part_by_id(SCHEMA, ptype, pid) for pid, ptype in ordered]
    gg = GoldenGateModule()
    design = gg.prepare_fragments(parts, assembly_type="linear")
    records = []
    for res in design:
        r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
        r.annotations["target_seq"] = res.get("target_seq", "")
        r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
        records.append(r)
    return records


def run_one(test_id, db):
    records = _design_records(test_id, db)
    names = [r.name for r in records]

    # (A) shim 경로
    shim_dir = Path(f"/tmp/runplan/shim/{test_id}")
    wf_a = Workflow(project_dir=shim_dir, parts=records)
    mdf = pd.DataFrame({"id": [test_id] * len(names), "name": names})
    wf_a.run_step3_module_assembly(module_df=mdf, method="golden_gate")
    a_recs = DNAFileReader.read_genbank(shim_dir / "modules" / "assembled" / f"{test_id}.gb")
    sha_a = _sha(str(a_recs[0].seq)) if a_recs else "-"

    # (B) build-plan 경로
    plan_dir = Path(f"/tmp/runplan/plan/{test_id}")
    wf_b = Workflow(project_dir=plan_dir, parts=records)
    plan = [{"op": "module_assembly", "inputs": names,
             "params": {"method": "golden_gate", "assembly_type": "linear"}}]
    out = wf_b.run_full_workflow(plan=plan)
    sha_b = "-"
    if out:
        b_recs = DNAFileReader.read_genbank(out[0]["file"])
        sha_b = _sha(str(b_recs[0].seq)) if b_recs else "-"

    status = "match" if (sha_a == sha_b and sha_a != "-") else "MISMATCH"
    return (test_id, status, sha_a, sha_b, out[0]["id"] if out else "-")


def main():
    db = DatabaseImport()
    ids = sys.argv[1:] or IDS
    ok = True
    for tid in ids:
        with contextlib.redirect_stdout(io.StringIO()):
            r = run_one(tid, db)
        if r[1] != "match":
            ok = False
        print(f"{r[0]}\t{r[1]}\tshim={r[2]}\tplan={r[3]}\tplan_id={r[4]}", flush=True)
    print("RESULT:", "ALL_MATCH" if ok else "HAS_MISMATCH", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
