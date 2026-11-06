"""주어진 Test_ID의 조립 결과 '실제 서열 + seqSHA'를 출력 — regression results.tsv의
seqSHA(앞 16자)를 실제 산출 서열과 매칭/확인하는 용도.

seqSHA 정의: sha256(assembled_seq.lower()).hexdigest()[:16]  (verify_qvs/verify_full과 동일).

사용(컨테이너):
  docker run --rm --network host --user $(id -u):$(id -g) -e HOME=/tmp \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src:ro" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev:ro" \
    awsemble:dev python tools_dev/inspect_assembly.py BBa_A340620 [BBa_B1101 ...]
"""
import sys, io, hashlib, contextlib
from pathlib import Path

import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAValidator
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow as Workflow

PROJECT_DIR = Path("/tmp/inspect/test_results")
SCHEMA = "Bioparts"


def main():
    ids = [a for a in sys.argv[1:] if not a.startswith("--")]
    if not ids:
        print("Test_ID를 인자로 주세요 (예: BBa_A340620)", file=sys.stderr)
        return
    for s in ["parts/insert", "modules/assembled"]:
        (PROJECT_DIR / s).mkdir(parents=True, exist_ok=True)
    db = DatabaseImport()

    for tid in ids:
        ordered = db.get_ordered_parts_by_test_id(tid)
        seq, valid = "", False
        with contextlib.redirect_stdout(io.StringIO()):
            cache = {pid: db.get_part_by_id(SCHEMA, ptype, pid) for pid, ptype in ordered}
            parts = [cache[pid] for pid, _ in ordered if cache.get(pid) is not None]
            gg = GoldenGateModule()
            design = gg.prepare_fragments(parts, assembly_type="linear")
            recs, names = [], []
            for res in design:
                r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
                r.annotations["target_seq"] = res.get("target_seq", "")
                r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
                recs.append(r); names.append(r.name)
            wf = Workflow(project_dir=PROJECT_DIR, parts=recs)
            wf.run_step3_module_assembly(
                module_df=pd.DataFrame({"id": [tid] * len(names), "name": names}),
                method="golden_gate",
            )
            out = PROJECT_DIR / "modules" / "assembled" / f"{tid}.gb"
            if out.exists():
                recs2 = DNAFileReader.read_genbank(out)
                seq = str(recs2[0].seq) if recs2 else ""
                valid = DNAValidator(db).check_parts_correct(out, tid)

        sha = hashlib.sha256(seq.lower().encode()).hexdigest()[:16]
        print(f"# {tid}  status={'pass' if valid else 'fail'}  len={len(seq)}bp  seqSHA={sha}")
        print(f"  parts(순서): {[pid for pid, _ in ordered]}")
        print(seq)
        print()


if __name__ == "__main__":
    main()
