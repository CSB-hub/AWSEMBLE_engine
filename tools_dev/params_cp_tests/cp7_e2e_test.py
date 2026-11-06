"""CP7 Stage C: build_plasmid circular 전 구간 E2E."""
import sys, io, contextlib, shutil, json
from pathlib import Path
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/tools_dev")
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from params import DesignParams
from make_intermediate_vector import build_vector_set, fetch_standard_spacers, STANDARD_SPACERS
S=Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
ok=True
def chk(l,c,e=""):
    global ok; ok=ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

db=DatabaseImport()
seqs=fetch_standard_spacers(db, STANDARD_SPACERS[:2])
imv=build_vector_set(seqs,"BsaI","BsmBI")[0]
SPEC=[("promoter","BBa_J23100"),("rbs","BBa_B0032"),("cds","BBa_C0012"),("terminator","BBa_B0010")]
parts=[db.get_part_by_id("Bioparts",t,i) for t,i in SPEC]
final=db.get_part_by_id("Bioparts","vector","pSB1A2")

print("circular 전 구간 실행")
d=S/"cp7_e2e"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
buf=io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        res=WorkflowBuilder.build_plasmid(parts, final, db, str(d),
            insert_method="gibson",
            params=DesignParams(module_assembly_type="circular"),
            intermediate_vector=imv)
except Exception as e:
    print(f"  [FAIL] 실행 중 예외: {type(e).__name__}: {e}")
    print("  --- 로그 마지막 ---"); print("\n".join(buf.getvalue().splitlines()[-8:]))
    sys.exit(1)

chk("construct 산출", res.get("construct") is not None)
con=res["construct"]
if con is not None:
    chk("최종 construct가 원형", getattr(con.seq,"circular",False), f"{len(con.seq)}bp")
    cs=str(con.seq).upper()*2
    for (t,i),p in zip(SPEC,parts):
        chk(f"{i} 포함", str(p.seq).upper() in cs)
    chk("최종 벡터 포함", str(final.seq).upper()[:40] in cs)
    # 중간 벡터 backbone은 남으면 안 된다 (excision으로 빠져야 함)
    from Bio.Seq import Seq
    imv_only = str(imv.seq).upper()[:30]
    chk("중간 벡터 backbone 잔류 없음", imv_only not in cs)

print("\n산출물")
for f in ["modules/assembled","modules/excised","constructs","reports"]:
    fp=d/f
    n=len(list(fp.glob('*'))) if fp.exists() else 0
    chk(f"{f} 생성", n>0, f"{n}개")

mf=d/"reports"/"manifest.json"
if mf.exists():
    m=json.loads(mf.read_text())
    ops=[e.get("op") for e in (m if isinstance(m,list) else m.get("artifacts",[]))]
    chk("계보에 module_excision 기록", "module_excision" in ops, str(ops))

print("\nlinear 경로 회귀 (같은 부품)")
d2=S/"cp7_e2e_lin"; shutil.rmtree(d2, ignore_errors=True); d2.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    res2=WorkflowBuilder.build_plasmid(parts, final, db, str(d2), insert_method="gibson",
                                       params=DesignParams())
chk("linear도 정상", res2.get("construct") is not None,
    f"{len(res2['construct'].seq)}bp / spacer {res2['l_spacer']}-{res2['r_spacer']}")

print("\n=== CP7 Stage C:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
