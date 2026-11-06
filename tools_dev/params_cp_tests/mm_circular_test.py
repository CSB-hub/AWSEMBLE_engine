"""multi-module circular: 벡터 하나로 module N개, 인접 module이 spacer를 공유하는가."""
import sys, io, contextlib, shutil, glob, json
from pathlib import Path
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from utils.io import DatabaseImport, DNAFileReader
from workflow.builder import WorkflowBuilder
from params import DesignParams
S=Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
ok=True
def chk(l,c,e=""):
    global ok; ok=ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))
db=DatabaseImport()
imv=db.get_part_by_id("Bioparts","vector","V0001").looped(); imv.name="V0001"
fin=db.get_part_by_id("Bioparts","vector","V0002").looped(); fin.name="V0002"
gA=[db.get_part_by_id("Bioparts",t,i) for t,i in [("promoter","BBa_J23100"),("rbs","BBa_B0032")]]
gB=[db.get_part_by_id("Bioparts",t,i) for t,i in [("terminator","BBa_B0010")]]

d=S/"mm_circ"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
buf=io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        res=WorkflowBuilder.build_multimodule_plasmid([gA,gB], fin, db, str(d),
            insert_method="gibson",
            params=DesignParams(module_assembly_type="circular"),
            intermediate_vector=imv)
except Exception as e:
    print(f"  [FAIL] 예외: {type(e).__name__}: {e}")
    print("\n".join(buf.getvalue().splitlines()[-8:])); sys.exit(1)

b=res["boundaries"]
print(f"경계 spacer: {b}  (module 2개 = 경계 3개)")
chk("construct 산출", res.get("construct") is not None,
    f"{len(res['construct'].seq)}bp" if res.get("construct") else "")
chk("벡터 하나만 사용", res.get("intermediate_vector")=="V0001")
chk("경계가 3개이고 인접이 공유", len(b)==3 and b[1] is not None, str(b))

print("\nmodule별로 꺼낸 조각의 말단이 배정된 spacer인가")
exs=sorted(glob.glob(str(d/"modules"/"excised"/"*.gb")))
chk("excised 2개", len(exs)==2, f"{len(exs)}개")
seqs={}
for f in exs:
    r=DNAFileReader.read_genbank(f)[0]; seqs[Path(f).stem]=str(r.seq).upper()
# 경계 서열은 **저장된 parts/raw**에서 읽는다 — DB만 보면 안 된다.
# spacer_mode 기본값(vector_homology)에서 양 끝 경계는 벡터 절단부에서 그 자리에 만들어져
# DB에 ID가 없다(V0002_homL/homR). 풀에서 온 연결 spacer만 DB에 있다.
sp={}
for f in sorted(glob.glob(str(d/"parts"/"raw"/"*.gb"))):
    r=DNAFileReader.read_genbank(f)[0]
    sp[r.name]=str(r.seq).upper()
for x in b:
    if x not in sp:
        rec=db.get_part_by_id("Bioparts","spacer",x)
        if rec is None:
            print(f"  [FAIL] 경계 '{x}' 서열을 parts/raw에도 DB에도 없어 확인 불가")
            ok=False; sp[x]=""
        else:
            sp[x]=str(rec.seq).upper()
for k,v in seqs.items():
    start=[x for x in b if v.startswith(sp[x][:25])]
    end=[x for x in b if v.endswith(sp[x][-25:])]
    print(f"    {k[:14]:<16} {len(v):>5}bp  시작={start} 끝={end}")
    chk(f"{k[:10]} 말단이 경계 spacer", bool(start) and bool(end))

print("\n인접 module이 같은 spacer를 공유하는가 (Gibson 연결 조건)")
vals=list(seqs.values())
if len(vals)==2:
    shared=[x for x in b if (vals[0].endswith(sp[x][-25:]) and vals[1].startswith(sp[x][:25]))
            or (vals[1].endswith(sp[x][-25:]) and vals[0].startswith(sp[x][:25]))]
    chk("공유 spacer 존재", bool(shared), str(shared))

print("\n계보")
mf=d/"reports"/"manifest.json"
if mf.exists():
    m=json.loads(mf.read_text()); es=m if isinstance(m,list) else m.get("artifacts",[])
    ops=[e.get("op") for e in es]
    chk("module_excision 2회 + multimodule_assembly", ops.count("module_excision")==2 and "multimodule_assembly" in ops, str(ops))

print("\n=== multi-module circular:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
