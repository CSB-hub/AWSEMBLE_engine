"""스모크: multi-module 정식 진입 build_multimodule_plasmid (P7).

KRIBB 패턴(공유 spacer로 module 연결)을 우리 Awsemble spacer로 재현 — igem/KRIBB core 부품 둘 다.
정답 대조가 아니라 메커니즘+구조(공유 spacer 연결 + circular construct) 입증.
부품은 spacer 제외 core(promoter/rbs/cds/terminator)부터 시작.

실행(컨테이너):
  docker run --rm --network host \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/smoke_multimodule.py [TEST_ID1 TEST_ID2 ...]
"""
import sys
import contextlib
import io

from utils.io import DatabaseImport, DNAFileReader
from workflow.builder import WorkflowBuilder

SCHEMA = "Bioparts"
VECTOR = "/app/gfp_test/pUC19.gb"


def parts_of(test_id, db):
    """test_id의 core 부품 — spacer(KRIBB KS 등)는 제외(우리 Awsemble spacer를 부착하므로)."""
    ordered = db.get_ordered_parts_by_test_id(test_id)
    core = [(pid, pt) for pid, pt in ordered if pt != "spacer"]
    return [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in core]


def main():
    ids = sys.argv[1:] if len(sys.argv) > 2 else ["BBa_C0630", "BBa_E0240"]
    db = DatabaseImport()
    vector = DNAFileReader.read_genbank(VECTOR)[0]
    groups = [parts_of(i, db) for i in ids]
    print(f"modules={ids} ({[len(g) for g in groups]}부품)")

    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = WorkflowBuilder.build_multimodule_plasmid(
                groups, vector, db, project_dir="/tmp/smoke_mm", insert_method="gibson")
        con = out["construct"]
        res = out["results"]
        mm = res[-2]["record"]   # multimodule (vector_insertion 직전)
        print(f"  boundaries(L,connectors,R)={out['boundaries']}")
        print(f"  ✅ multimodule {len(mm.seq)}bp | construct {len(con.seq)}bp "
              f"circular={con.circular} | stage {len(res)}개")
    except Exception as e:
        import traceback
        print(f"  ❌ {type(e).__name__}: {str(e)[:100]}")
        traceback.print_exc()


if __name__ == "__main__":
    main()
