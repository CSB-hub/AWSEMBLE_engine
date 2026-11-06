"""스모크: 정석 build_plasmid 파이프라인 (P7) — 메인코드 진입함수 작동 입증.

WorkflowBuilder.build_plasmid 한 호출로:
  parts(유저) + vector
    → 선정·부착(SpacerSelector) → prepare_fragments
    → run_plan([module_assembly(GG), vector_insertion(method)])
    → construct + reports/run_report.md

실행(컨테이너):
  docker run --rm --network host \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/smoke_spacer_pipeline.py [TEST_ID]
"""
import sys
import contextlib
import io

from utils.io import DatabaseImport, DNAFileReader
from workflow.builder import WorkflowBuilder

SCHEMA = "Bioparts"
VECTORS = {
    "pUC19": ("/app/gfp_test/pUC19.gb", "gibson"),
    "pIB2": ("/app/gfp_test/pIB2.gbk", "gibson"),
    "pGGAselect": ("/app/gfp_test/snapgene_simulation/golden_gate/pGGAselect.gb", "golden_gate"),
}


def main():
    test_id = sys.argv[1] if len(sys.argv) > 1 else "BBa_C0630"
    db = DatabaseImport()
    ordered = db.get_ordered_parts_by_test_id(test_id)
    parts = [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in ordered]
    print(f"[{test_id}] 부품 {len(parts)}개")

    for vname, (vpath, method) in VECTORS.items():
        vector = DNAFileReader.read_genbank(vpath)[0]
        print(f"\n{'='*56}\n[{vname} / {method}]")
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                out = WorkflowBuilder.build_plasmid(
                    parts, vector, db,
                    project_dir=f"/tmp/smoke_pipe/{vname}",
                    insert_method=method,
                )
            con = out["construct"]
            L, R = out["l_spacer"], out["r_spacer"]
            cs = str(con.seq).upper()
            sL = str(db.get_part_by_id(SCHEMA, "spacer", L).seq).upper()
            sR = str(db.get_part_by_id(SCHEMA, "spacer", R).seq).upper()
            ok = (cs + cs).find(sL) >= 0 and (cs + cs).find(sR) >= 0
            print(f"  ✅ L={L} R={R} | construct {len(con.seq)}bp circular={con.circular} "
                  f"| spacer보존={ok} | stage {len(out['results'])}개")
        except Exception as e:
            print(f"  ❌ {type(e).__name__}: {str(e)[:100]}")


if __name__ == "__main__":
    main()
