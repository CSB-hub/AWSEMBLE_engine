"""스모크: 양끝 spacer 부착 module → vector 삽입 (P7 구멍1 검증).

기존 스모크(smoke_module_vector.py)는 spacer 없는 module이 vector 삽입 시
프라이머 anneal 실패함을 보였다. 이번엔 module 양끝에 spacer(S0001/S0002)를 부착해
같은 vector_insertion(기존 코드)이 작동하는지 확인한다.

실행(컨테이너):
  docker run --rm --network host \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/smoke_spacer_insertion.py [TEST_ID ...]
"""
import sys
from pathlib import Path

import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow
from workflow.classical_assemble import ClassicalWorkflow

PROJECT_DIR = Path("/tmp/smoke_sp/biopart")
SCHEMA = "Bioparts"
SPACER_L, SPACER_R = "S0002", "S0005"          # anneal 적합 spacer (test_spacer_anneal 확인)
VECTORS = {
    "pUC19": ("/app/gfp_test/pUC19.gb", "gibson"),
    "pIB2": ("/app/gfp_test/pIB2.gbk", "gibson"),
    "pGGAselect": ("/app/gfp_test/snapgene_simulation/golden_gate/pGGAselect.gb", "golden_gate"),
}


def build_module(test_id, db):
    ordered = db.get_ordered_parts_by_test_id(test_id)
    parts = [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in ordered]
    gg = GoldenGateModule()
    design = gg.prepare_fragments(parts, assembly_type="linear")
    records, names = [], []
    for res in design:
        r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
        r.annotations["target_seq"] = res.get("target_seq", "")
        r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
        records.append(r); names.append(r.name)
    wf = BiopartWorkflow(project_dir=PROJECT_DIR, parts=records)
    wf.run_step3_module_assembly(module_df=pd.DataFrame({"id": [test_id]*len(names), "name": names}),
                                 method="golden_gate")
    return DNAFileReader.read_genbank(PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb")[0]


def attach_spacers(module, sL, sR):
    """module 양끝에 spacer 서열 부착 (검증용 단순 결합 — anneal 발판 효과 확인)."""
    seq = str(sL.seq) + str(module.seq) + str(sR.seq)
    m = Dseqrecord(seq, name="module_with_spacers")
    return m


def smoke_one(test_id, db):
    print(f"\n{'='*60}\n[{test_id}] 양끝 spacer 부착 → vector 삽입\n{'='*60}")
    module = build_module(test_id, db)
    sL = db.get_part_by_id(SCHEMA, "spacer", SPACER_L)
    sR = db.get_part_by_id(SCHEMA, "spacer", SPACER_R)
    mod_sp = attach_spacers(module, sL, sR)
    ms = str(mod_sp.seq).upper()
    bsai = ms.count("GGTCTC") + ms.count("GAGACC")   # BsaI(+RC) — GG 삽입 시 insert에 있으면 안 됨
    print(f"  module {len(module.seq)}bp + spacer({SPACER_L}/{SPACER_R}) → {len(mod_sp.seq)}bp "
          f"| BsaI site={bsai}")

    for vname, (vpath, method) in VECTORS.items():
        vector = DNAFileReader.read_genbank(vpath)[0]
        wf = ClassicalWorkflow(output_dir=f"/tmp/smoke_sp/classical_{vname}")
        try:
            results = wf.assemble(mod_sp, vector, method=method, validate=False)
            c = results[0]
            cs = str(c.seq).upper()
            contained = (cs + cs).find(str(mod_sp.seq).upper()) >= 0
            print(f"  ✅ [{vname}/{method}] {len(c.seq)}bp circular={c.circular} | insert 보존={contained}")
        except Exception as e:
            print(f"  ❌ [{vname}/{method}] {type(e).__name__}: {str(e)[:90]}")


def main():
    ids = sys.argv[1:] or ["BBa_C0630"]
    db = DatabaseImport()
    for tid in ids:
        try:
            smoke_one(tid, db)
        except Exception as e:
            import traceback; print(f"  ❌ [{tid}] 중단: {e}"); traceback.print_exc()


if __name__ == "__main__":
    main()
