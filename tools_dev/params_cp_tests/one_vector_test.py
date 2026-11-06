"""V0001 하나로 module 2개를 만들면 말단이 서로 구분되는가 (multi-module 가능 여부)."""
import sys, io, contextlib, shutil, glob
from pathlib import Path
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
from utils.io import DatabaseImport, DNAFileReader
from workflow.builder import WorkflowBuilder
from params import DesignParams as _DP
# ⚠️ spacer_mode="library"로 고정한다. 2026-10-20에 기본값이 "vector_homology"로
#    바뀌었고 그 모드는 **최종 서열이 달라진다**(중립 spacer 흉터가 사라지므로).
#    이 테스트들은 기존 방식의 **하위호환 보증선**이므로 모드를 명시해 고정한다.
#    새 모드의 기준은 vector_homology_test.py가 따로 든다.
def DesignParams(**kw):
    kw.setdefault('spacer_mode', 'library')
    return _DP(**kw)

S=Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
ok=True
def chk(l,c,e=""):
    global ok; ok=ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))
db=DatabaseImport()
imv=db.get_part_by_id("Bioparts","vector","V0001").looped(); imv.name="V0001"
fin=db.get_part_by_id("Bioparts","vector","V0002").looped(); fin.name="V0002"

def run(tag, spec):
    d=S/f"onev_{tag}"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
    ps=[db.get_part_by_id("Bioparts",t,i) for t,i in spec]
    with contextlib.redirect_stdout(io.StringIO()) as b:
        res=WorkflowBuilder.build_plasmid(ps, fin, db, str(d), insert_method="gibson",
            params=DesignParams(module_assembly_type="circular"), intermediate_vector=imv)
    ex=DNAFileReader.read_genbank(glob.glob(str(d/"modules"/"excised"/"*.gb"))[0])[0]
    return res, ex, b.getvalue()

print("module A / module B — 같은 V0001 하나로")
resA, exA, logA = run("A", [("promoter","BBa_J23100"),("rbs","BBa_B0032")])
resB, exB, logB = run("B", [("cds","BBa_C0012"),("terminator","BBa_B0010")])
a,bq=str(exA.seq).upper(), str(exB.seq).upper()
print(f"  A: spacer {resA['l_spacer']}-{resA['r_spacer']}, excised {len(a)}bp")
print(f"  B: spacer {resB['l_spacer']}-{resB['r_spacer']}, excised {len(bq)}bp")
chk("A construct 산출", resA.get("construct") is not None)
chk("B construct 산출", resB.get("construct") is not None)
chk("spacer가 부품으로 선정됨(벡터 없이)", resA['l_spacer'] is not None)

print("\n말단이 spacer인가 (= module별 구분되는 homology)")
for tag,ex,res in (("A",a,resA),("B",bq,resB)):
    ls=str(db.get_part_by_id("Bioparts","spacer",res['l_spacer']).seq).upper()
    rs=str(db.get_part_by_id("Bioparts","spacer",res['r_spacer']).seq).upper()
    chk(f"{tag} 앞이 L spacer로 시작", ex.startswith(ls[:25]), ex[:25])
    chk(f"{tag} 뒤가 R spacer로 끝", ex.endswith(rs[-25:]), ex[-25:])

print("\n벡터 유래 말단이 남지 않았는가")
vs=str(imv.seq).upper()
chk("A 말단이 벡터 서열 아님", a[:25] not in vs)
chk("B 말단이 벡터 서열 아님", bq[:25] not in vs)

print("\nPCR 전환 경고")
chk("digest→PCR 안내 출력", "PCR(spacer 앵커)로 진행" in logA)

print("\n=== V0001 하나로 multi-module:", "가능" if ok else "불가", "===")
sys.exit(0 if ok else 1)
