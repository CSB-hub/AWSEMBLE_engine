"""165 spacer를 양끝 anneal 적합성으로 전수 테스트 → 통과 spacer로 vector_insertion 검증.

프라이머 설계 실패(high tm/hairpin)는 검증 시스템 정상 작동. 제약 완화 대신
anneal 발판으로 적합한 spacer(시작+끝 모두 프라이머 통과)를 165 풀에서 선별한다.

실행:
  docker run --rm --network host -v "$PWD/awsemble/.env:/app/.env:ro" \
    -v "$PWD/awsemble/src:/app/src" -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/test_spacer_anneal.py [TEST_ID]
"""
import sys, contextlib, io
from pathlib import Path
import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow
from workflow.classical_assemble import ClassicalWorkflow
from basic.primer_design import PrimerDesign

PROJECT_DIR = Path("/tmp/ts/biopart")
SCHEMA = "Bioparts"


def build_module(test_id, db):
    ordered = db.get_ordered_parts_by_test_id(test_id)
    parts = [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in ordered]
    gg = GoldenGateModule()
    with contextlib.redirect_stdout(io.StringIO()):
        design = gg.prepare_fragments(parts, assembly_type="linear")
        recs, names = [], []
        for r in design:
            d = Dseqrecord(r["full_sequence"], name=r["part_name"])
            d.annotations["target_seq"] = r.get("target_seq", "")
            d.annotations["restriction_enzyme"] = r.get("restriction_enzyme", "BsaI")
            recs.append(d); names.append(d.name)
        wf = BiopartWorkflow(project_dir=PROJECT_DIR, parts=recs)
        wf.run_step3_module_assembly(module_df=pd.DataFrame({"id": [test_id]*len(names), "name": names}),
                                     method="golden_gate")
    return DNAFileReader.read_genbank(PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb")[0]


def anneal_ok(seq):
    """[spacer+module+spacer] 양끝 anchor 프라이머가 설계되면 True (시작+끝 모두 적합)."""
    with contextlib.redirect_stdout(io.StringIO()):
        res = PrimerDesign.design_basic_primers(
            template=seq, target_start=0, target_end=len(seq),
            anchor_ends=True, check_specificity=False)
    return bool(res.get("forward") and res.get("reverse")), res


def main():
    test_id = sys.argv[1] if len(sys.argv) > 1 else "BBa_E0240"
    db = DatabaseImport()
    module = build_module(test_id, db)
    print(f"module [{test_id}] {len(module.seq)}bp")

    spacers = db._execute("Bioparts", "SELECT ID,Sequence FROM spacer WHERE Institution='Awsemble' ORDER BY ID")
    ok = []
    for sid, sseq in spacers:
        seq = sseq + str(module.seq) + sseq
        good, _ = anneal_ok(seq)
        if good:
            ok.append(sid)
    print(f"\n=== anneal 적합 spacer: {len(ok)}/{len(spacers)} ===")
    print("  적합 샘플:", ok[:25])

    if len(ok) >= 2:
        L, R = ok[0], ok[1]
        sL = db.get_part_by_id(SCHEMA, "spacer", L)
        sR = db.get_part_by_id(SCHEMA, "spacer", R)
        mod_sp = Dseqrecord(str(sL.seq) + str(module.seq) + str(sR.seq))
        print(f"\n=== vector_insertion 검증 (L={L}, R={R}) ===")
        for vname, vpath in [("pUC19", "/app/gfp_test/pUC19.gb"), ("pIB2", "/app/gfp_test/pIB2.gbk")]:
            vector = DNAFileReader.read_genbank(vpath)[0]
            wf = ClassicalWorkflow(output_dir=f"/tmp/ts/cl_{vname}")
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    results = wf.assemble(mod_sp, vector, method="gibson", validate=False)
                c = results[0]
                contained = (str(c.seq).upper()*2).find(str(mod_sp.seq).upper()) >= 0
                print(f"  ✅ [{vname}] {len(c.seq)}bp circular={c.circular} insert보존={contained}")
            except Exception as e:
                print(f"  ❌ [{vname}] {type(e).__name__}: {str(e)[:80]}")
    else:
        print("\n⚠️ 적합 spacer 부족 → 제작 규칙에 anneal 적합성(양끝 Tm/hairpin) 추가 필요")


if __name__ == "__main__":
    main()
