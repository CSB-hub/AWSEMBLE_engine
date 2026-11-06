"""QVS(Quick Verify Set) 검증 하네스 — 리팩토링 동작보존 확인용.

- production 경로(GoldenGateModule.prepare_fragments → workflow assemble → DNAValidator)를
  대표 iGEM Test_ID 15건에 대해 돌려 `id<TAB>status<TAB>seqSHA` 컴팩트 출력.
- 풀런과 충돌 안 나게 PROJECT_DIR=/tmp/qvs 로 격리, 공유 tests 디렉토리에 안 씀.
- behavior-preserving 단계(P0/P1)는 이 출력이 전후 100% 동일해야 통과.

실행(컨테이너):
  docker run --rm --network host -v "$PWD/awsemble/.env:/app/.env:ro" \
    -v "$PWD/awsemble/src:/app/src" -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    awsemble:dev python tools_dev/verify_qvs.py
"""
import sys, hashlib, io, contextlib, os
from pathlib import Path

import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAValidator
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow as Workflow

PROJECT_DIR = Path("/tmp/qvs/test_results")
SCHEMA = "Bioparts"

QVS_IDS = [
    "BBa_B1101",                      # 2 parts
    "BBa_C0630", "BBa_E0130", "BBa_E0241",   # 3
    "BBa_C0490", "BBa_E0240",         # 4
    "BBa_A340620", "BBa_B1201", "BBa_B1202", # 5
    "BBa_C0460", "BBa_C0500", "BBa_E0666",   # 6
    "BBa_C0480", "BBa_I0401", "BBa_I0402",   # 7
]


def run_one(test_id, db):
    """단일 Test_ID 조립+검증. (id, status, seq_sha, detail) 반환."""
    try:
        ordered = db.get_ordered_parts_by_test_id(test_id)
        if not ordered:
            return (test_id, "no_def", "-", "no ordered parts")
        cache = {}
        for pid, ptype in ordered:
            d = db.get_part_by_id(SCHEMA, ptype, pid)
            if d is None:
                return (test_id, "fail", "-", f"part missing: {pid}")
            cache[pid] = d
        expected = db.get_construct_sequence_by_test_id(test_id)

        parts = [cache[pid] for pid, _ in ordered]
        gg = GoldenGateModule()
        design = gg.prepare_fragments(parts, assembly_type="linear")

        records, names = [], []
        for res in design:
            r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
            r.annotations["target_seq"] = res.get("target_seq", "")
            r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
            r.annotations["part_id"] = res.get("part_id", "")
            r.annotations["part_type"] = res.get("part_type", "")
            records.append(r); names.append(r.name)

        wf = Workflow(project_dir=PROJECT_DIR, parts=records)
        module_df = pd.DataFrame({"id": [test_id] * len(names), "name": names})
        wf.run_step3_module_assembly(module_df=module_df, method="golden_gate")

        out = PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb"
        if not out.exists():
            return (test_id, "no_output", "-", "output not generated")

        recs = DNAFileReader.read_genbank(out)
        seq = str(recs[0].seq) if recs else ""
        sha = hashlib.sha256(seq.lower().encode()).hexdigest()[:16]

        # 성공 기준 = 부품 서열 + 입력 순서 동일 (check_parts_correct).
        # 전체 서열 exact-match는 overhang/method 차이로 달라질 수 있어 기준 아님.
        valid = DNAValidator(db).check_parts_correct(out, test_id)
        status = "pass" if valid else "fail"
        return (test_id, status, sha, "")
    except Exception as e:
        return (test_id, "error", "-", str(e).split("\n")[0])


def main():
    for s in ["parts/insert", "modules/assembled"]:
        (PROJECT_DIR / s).mkdir(parents=True, exist_ok=True)
    db = DatabaseImport()
    ids = sys.argv[1:] or QVS_IDS
    for tid in ids:
        # 조립 과정의 verbose stdout 억제 → 결과 한 줄만 출력
        with contextlib.redirect_stdout(io.StringIO()):
            r = run_one(tid, db)
        print(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}", flush=True)


if __name__ == "__main__":
    main()
