import sys, io, contextlib
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
from pydna.dseqrecord import Dseqrecord
from assembly.golden_gate import GoldenGateModule
from utils.io import DatabaseImport
from params import SUPPORTED_ENZYMES

db = DatabaseImport()
SPEC = [("promoter","BBa_J23100"),("rbs","BBa_B0032"),("cds","BBa_C0012"),("terminator","BBa_B0010")]
parts = [db.get_part_by_id("Bioparts", t, i) for t, i in SPEC]
print(f"부품 {len(parts)}개 로드\n")
ok = True
for enz in SUPPORTED_ENZYMES:
    gg = GoldenGateModule(enz)
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            design = gg.prepare_fragments(parts, assembly_type="linear")
            recs = [Dseqrecord(d["full_sequence"], name=d["part_name"]) for d in design]
            frags = [gg.extract_insert_fragment(r) for r in recs]
            if any(f is None for f in frags):
                raise RuntimeError("fragment 추출 실패")
            res = gg.assemble(frags, assembly_type="linear")
        used = design[0]["restriction_enzyme"]
        mod = res[0]
        # 인접 fragment 접합부 overhang 길이가 효소 파생값과 맞는지
        oh = GoldenGateModule(used).overhang_len
        junc_ok = all(str(frags[i].seq)[-oh:].upper() == str(frags[i+1].seq)[:oh].upper()
                      for i in range(len(frags)-1))
        switched = "" if used == enz else f"  (내부 사이트로 {used}로 전환됨)"
        print(f"  [{'OK ' if junc_ok else 'FAIL'}] {enz:<6} overhang={oh}  module={len(mod.seq)}bp  "
              f"접합부 {len(frags)-1}개 일치={junc_ok}{switched}")
        ok = ok and junc_ok
    except Exception as e:
        print(f"  [FAIL] {enz:<6} {type(e).__name__}: {str(e)[:120]}")
        ok = False

print("\n미지원 효소 차단 확인")
try:
    GoldenGateModule("EcoRI"); print("  [FAIL] EcoRI가 통과됨"); ok = False
except Exception as e:
    print(f"  [OK ] EcoRI 거부 — {type(e).__name__}")
from basic.primer_design import PrimerDesign
try:
    PrimerDesign.design_goldengate_primers(insert=parts[0], left_overhang="AAAA",
                                           right_overhang="TTTT", restriction_enzyme="EcoRI")
    print("  [FAIL] 모르는 효소가 조용히 BsaI로 대체됨"); ok = False
except ValueError as e:
    print(f"  [OK ] 조용한 BsaI 대체 제거 — {str(e)[:70]}")

print("\n=== CP2:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
