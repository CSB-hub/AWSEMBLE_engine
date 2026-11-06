"""P7 구멍1 스모크 종합 리포트 — 결과 영구 저장 (내일 코드 적용계획 검토용).

저장물(/out = develop_log/ljw/p7_smoke_results):
  - spacer_anneal.tsv     : 165 spacer 양끝 anneal 적합성 전수 (적합/부적합 + 실패이유 primer3 explain)
  - insertion_results.tsv : module×vector×method 삽입 결과 (성공/실패 + BsaI site + 이유)
  - constructs/*.gb       : 성공 construct (insert/spacer feature 기입 — P6 annotate)

실행:
  docker run --rm --network host --user $(id -u):$(id -g) -e HOME=/tmp \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    -v "$PWD/develop_log/ljw/p7_smoke_results:/out" \
    awsemble:dev python tools_dev/p7_smoke_report.py
"""
import contextlib, io
from pathlib import Path
import pandas as pd, primer3
from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAFileWriter
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow
from workflow.classical_assemble import ClassicalWorkflow
from assembly.common import annotate_inserted_region, remove_duplicate_features
from basic.primer_design import PrimerDesign

OUT = Path("/out")
SCHEMA = "Bioparts"
SL, SR = "S0002", "S0005"   # anneal 적합 (L/R)
MODULES = ["BBa_B1101", "BBa_C0630", "BBa_E0130", "BBa_E0241", "BBa_C0490",
           "BBa_E0240", "BBa_C0460", "BBa_C0480", "BBa_I0402"]
CASES = [("pUC19", "/app/gfp_test/pUC19.gb", "gibson"),
         ("pIB2", "/app/gfp_test/pIB2.gbk", "gibson"),
         ("pGGAselect", "/app/gfp_test/snapgene_simulation/golden_gate/pGGAselect.gb", "golden_gate")]


def build(tid, db):
    ordered = db.get_ordered_parts_by_test_id(tid)
    parts = [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in ordered]
    gg = GoldenGateModule()
    with contextlib.redirect_stdout(io.StringIO()):
        design = gg.prepare_fragments(parts, assembly_type="linear")
        recs = []
        for r in design:
            d = Dseqrecord(r["full_sequence"], name=r["part_name"])
            d.annotations["target_seq"] = r.get("target_seq", "")
            d.annotations["restriction_enzyme"] = r.get("restriction_enzyme", "BsaI")
            d.annotations["part_id"] = r.get("part_id", "")
            d.annotations["part_type"] = r.get("part_type", "")
            recs.append(d)
        wf = BiopartWorkflow(project_dir="/tmp/rep", parts=recs)
        wf.run_step3_module_assembly(
            module_df=pd.DataFrame({"id": [tid]*len(recs), "name": [d.name for d in recs]}),
            method="golden_gate")
    return DNAFileReader.read_genbank(f"/tmp/rep/modules/assembled/{tid}.gb")[0]


def anneal_detail(seq):
    with contextlib.redirect_stdout(io.StringIO()):
        res = PrimerDesign.design_basic_primers(template=seq, target_start=0, target_end=len(seq),
                                                anchor_ends=True, check_specificity=False)
    ok = bool(res.get("forward") and res.get("reverse"))
    if ok:
        return ok, res.get("fwd_tm"), res.get("rev_tm"), ""
    s = {"PRIMER_EXPLAIN_FLAG": 1, "SEQUENCE_FORCE_LEFT_START": 0,
         "SEQUENCE_FORCE_RIGHT_START": len(seq)-1,
         "PRIMER_PRODUCT_SIZE_RANGE": [[len(seq), len(seq)]],
         "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 25}
    out = primer3.bindings.design_primers({"SEQUENCE_ID": "x", "SEQUENCE_TEMPLATE": seq}, s)
    return ok, None, None, f"L[{out.get('PRIMER_LEFT_EXPLAIN')}] R[{out.get('PRIMER_RIGHT_EXPLAIN')}]"


def make_mod_sp(module, sL, sR):
    """spacer 부착 + module/spacer feature 좌표 이동 보존."""
    seq = str(sL.seq) + str(module.seq) + str(sR.seq)
    m = Dseqrecord(seq, name="module_with_spacers")
    off = len(sL.seq)
    for f in module.features:
        m.features.append(SeqFeature(
            FeatureLocation(int(f.location.start)+off, int(f.location.end)+off, f.location.strand),
            type=f.type, qualifiers=dict(f.qualifiers)))
    m.features.append(SeqFeature(FeatureLocation(0, len(sL.seq), 1),
                      type="misc_feature", qualifiers={"label": [f"spacer_L_{SL}"]}))
    m.features.append(SeqFeature(FeatureLocation(len(seq)-len(sR.seq), len(seq), 1),
                      type="misc_feature", qualifiers={"label": [f"spacer_R_{SR}"]}))
    return m


def main():
    db = DatabaseImport()
    (OUT / "constructs").mkdir(parents=True, exist_ok=True)

    # 1. spacer anneal 적합성 전수
    ref = build("BBa_C0630", db)
    spacers = db._execute("Bioparts", "SELECT ID,Sequence FROM spacer WHERE Institution='Awsemble' ORDER BY ID")
    nok = 0
    with open(OUT / "spacer_anneal.tsv", "w") as f:
        f.write("ID\tlen\tok\tfwd_tm\trev_tm\tfail_reason\n")
        for sid, sseq in spacers:
            ok, ft, rt, reason = anneal_detail(sseq + str(ref.seq) + sseq)
            nok += ok
            f.write(f"{sid}\t{len(sseq)}\t{ok}\t{ft}\t{rt}\t{reason}\n")
    print(f"[1] spacer anneal: {nok}/{len(spacers)} 적합 → spacer_anneal.tsv")

    # 2. module × vector × method 삽입 (feature 기입 construct 저장)
    sL = db.get_part_by_id(SCHEMA, "spacer", SL)
    sR = db.get_part_by_id(SCHEMA, "spacer", SR)
    wf = ClassicalWorkflow(output_dir="/tmp/cl")
    rows = []
    for mid in MODULES:
        module = build(mid, db)
        mod_sp_feat = make_mod_sp(module, sL, sR)
        mod_sp_seq = Dseqrecord(str(mod_sp_feat.seq), name="module_with_spacers")  # assemble용(feature 없음→전체 target)
        ms = str(mod_sp_feat.seq).upper()
        bsai = ms.count("GGTCTC") + ms.count("GAGACC")
        for vname, vpath, method in CASES:
            vector = DNAFileReader.read_genbank(vpath)[0]
            try:
                with contextlib.redirect_stdout(io.StringIO()):
                    results = wf.assemble(mod_sp_seq, vector, method=method, validate=False)
                    construct = results[0]
                    # feature 기입 (P6 annotate, target=insert 전체)
                    construct = annotate_inserted_region(construct, mod_sp_feat, 0, len(mod_sp_feat.seq))
                    construct = remove_duplicate_features(construct)
                    DNAFileWriter.write_genbank(construct, str(OUT / "constructs" / f"{mid}__{vname}__{method}.gb"))
                rows.append((mid, bsai, vname, method, "OK", len(construct.seq), len(construct.features), ""))
            except Exception as e:
                rows.append((mid, bsai, vname, method, "FAIL", 0, 0, str(e).split("\n")[0][:110]))
    with open(OUT / "insertion_results.tsv", "w") as f:
        f.write("module\tBsaI_site\tvector\tmethod\tstatus\tconstruct_len\tn_features\treason\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    ok = sum(1 for r in rows if r[4] == "OK")
    print(f"[2] insertion: {ok}/{len(rows)} 성공 → insertion_results.tsv + constructs/*.gb")
    for r in rows:
        print(f"    {r[0]:12} BsaI={r[1]} {r[2]:11} {r[3]:11} {r[4]:4} len={r[5]} feat={r[6]} {r[7]}")


if __name__ == "__main__":
    main()
