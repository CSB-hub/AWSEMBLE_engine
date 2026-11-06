"""실제 DB 벡터 회귀: V0001(중간 클로닝) + V0002(최종)로 circular 전 구간이 도는가."""
import sys, io, contextlib, shutil, json
from pathlib import Path
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from Bio import Restriction
from Bio.Seq import Seq
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from params import DesignParams
S=Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
db=DatabaseImport()
imv=db.get_part_by_id("Bioparts","vector","V0001")
fin=db.get_part_by_id("Bioparts","vector","V0002")
imv=imv.looped() if not imv.circular else imv
fin=fin.looped() if not fin.circular else fin
imv.name, fin.name = "V0001", "V0002"
SPEC=[("promoter","BBa_J23100"),("rbs","BBa_B0032"),("cds","BBa_C0012"),("terminator","BBa_B0010")]
parts=[db.get_part_by_id("Bioparts",t,i) for t,i in SPEC]

d=S/"v0001_run"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
buf=io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        res=WorkflowBuilder.build_plasmid(parts, fin, db, str(d), insert_method="gibson",
            params=DesignParams(module_assembly_type="circular"), intermediate_vector=imv)
except Exception as e:
    print(f"❌ 실패: {type(e).__name__}: {e}")
    print("--- 로그 마지막 12줄 ---"); print("\n".join(buf.getvalue().splitlines()[-12:]))
    sys.exit(1)

log=buf.getvalue()
con=res["construct"]
print(f"✅ construct {len(con.seq)}bp (circular={getattr(con.seq,'circular',False)})")
for ln in log.splitlines():
    if any(k in ln for k in ("excision","circular assembly","🧪","✂️","⚠️")):
        print("  "+ln.strip())

print("\n검증")
cs=str(con.seq).upper()*2
for (t,i),p in zip(SPEC,parts):
    print(f"  {i:<12} 포함: {str(p.seq).upper() in cs}")
print(f"  최종벡터 V0002 포함: {str(fin.seq).upper()[:50] in cs}")
n_bsa=len(Restriction.BsaI.search(Seq(str(con.seq))))
print(f"  BsaI 잔류: {n_bsa}개")

# 중간 산물 크기
import glob, os
for b in ("modules/assembled","modules/excised"):
    for f in glob.glob(str(d/b/"*.gb")):
        from utils.io import DNAFileReader
        r=DNAFileReader.read_genbank(f)[0]
        print(f"  {b}: {len(r.seq)}bp")

ok = (con is not None and getattr(con.seq,"circular",False)
      and all(str(p.seq).upper() in cs for p in parts))
print("\n=== 실제 벡터(V0001/V0002) circular:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
