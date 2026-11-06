"""스모크: biopart module → vector(pIB2) Gibson 삽입 갭 확인 (P7 첫 작업).

목적(데이터로 갭 확인):
  - 대표 module(기본 BBa_C0630)을 biopart 워크플로우로 조립해 modules/assembled/{id}.gb 생성.
  - 그 module을 insert로, pIB2를 vector로 ClassicalWorkflow.assemble(method="gibson") 시도.
  - module엔 CDS feature가 없으므로 _find_target_feature가 무엇을 잡는지, vector
    linearize/조립이 되는지를 그대로 노출한다(수리 전 현황 실측).

실행(컨테이너):
  docker run --rm --network host \
    -v "$PWD/awsemble/.env:/app/.env:ro" \
    -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/awsemble/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/smoke_module_vector.py [TEST_ID ...]
"""
import sys
from pathlib import Path

import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader
from assembly.golden_gate import GoldenGateModule
from assembly.factory import AssemblyFactory
from workflow.biopart_assemble import BiopartWorkflow
from workflow.classical_assemble import ClassicalWorkflow

PROJECT_DIR = Path("/tmp/smoke/biopart")
VECTOR_FILE = "/app/gfp_test/pIB2.gbk"
SCHEMA = "Bioparts"


def build_module(test_id, db):
    """biopart 워크플로우로 module 조립 → modules/assembled/{id}.gb 경로 반환."""
    ordered = db.get_ordered_parts_by_test_id(test_id)
    if not ordered:
        raise RuntimeError(f"no ordered parts: {test_id}")
    parts = []
    for pid, ptype in ordered:
        d = db.get_part_by_id(SCHEMA, ptype, pid)
        if d is None:
            raise RuntimeError(f"part missing: {pid}")
        parts.append(d)

    gg = GoldenGateModule()
    design = gg.prepare_fragments(parts, assembly_type="linear")
    records, names = [], []
    for res in design:
        r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
        r.annotations["target_seq"] = res.get("target_seq", "")
        r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
        records.append(r); names.append(r.name)

    wf = BiopartWorkflow(project_dir=PROJECT_DIR, parts=records)
    module_df = pd.DataFrame({"id": [test_id] * len(names), "name": names})
    wf.run_step3_module_assembly(module_df=module_df, method="golden_gate")
    out = PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb"
    if not out.exists():
        raise RuntimeError(f"module not generated: {test_id}")
    return out


def smoke_one(test_id, db):
    print(f"\n{'='*60}\n[{test_id}] module → pIB2 Gibson 삽입 스모크\n{'='*60}")

    # 1. module 조립
    module_path = build_module(test_id, db)
    module = DNAFileReader.read_genbank(module_path)[0]
    print(f"  module: {len(module.seq)} bp, circular={module.circular}, "
          f"features={len(module.features)}")
    ftypes = {}
    for f in module.features:
        ftypes[f.type] = ftypes.get(f.type, 0) + 1
    print(f"  module feature types: {ftypes}")
    has_cds = any(f.type.lower() == "cds" for f in module.features)
    print(f"  >> module에 CDS feature 있음? {has_cds}")

    # 2. _find_target_feature가 module에서 무엇을 잡는지 (갭 A의 핵심)
    gib = AssemblyFactory.create_assembler("gibson")
    ts, te = gib._find_target_feature(module, None, None)
    print(f"  >> _find_target_feature → ({ts}, {te}) / 길이 {te-ts} "
          f"(module 전체 {len(module.seq)} bp 대비)")

    # 3. vector(pIB2) 로드 + linearize 가능성
    vector = DNAFileReader.read_genbank(VECTOR_FILE)[0]
    print(f"  vector pIB2: {len(vector.seq)} bp, circular={vector.circular}")

    # 4. 실제 Gibson 삽입 시도 (현 코드 그대로 — 갭 노출)
    wf = ClassicalWorkflow(output_dir="/tmp/smoke/classical")
    try:
        results = wf.assemble(module, vector, method="gibson", validate=False)
        c = results[0]
        print(f"  ✅ 조립 성공: {len(c.seq)} bp, circular={c.circular}")
        # 구조 검증: module 전체 서열이 construct 안에 등장하는가?
        cs = str(c.seq).upper()
        ms = str(module.seq).upper()
        contained = (cs + cs).find(ms) >= 0
        print(f"  >> construct에 module 전체 서열 포함? {contained}")
    except Exception as e:
        import traceback
        print(f"  ❌ 조립 실패: {type(e).__name__}: {e}")
        traceback.print_exc()


def main():
    ids = sys.argv[1:] or ["BBa_C0630"]
    db = DatabaseImport()
    for tid in ids:
        try:
            smoke_one(tid, db)
        except Exception as e:
            import traceback
            print(f"  ❌ [{tid}] 스모크 중단: {type(e).__name__}: {e}")
            traceback.print_exc()


if __name__ == "__main__":
    main()
