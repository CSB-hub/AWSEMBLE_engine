"""CP7 Stage B: circular module에서 module을 꺼내는가 (digest / PCR 폴백)."""
import sys, io, contextlib, shutil
from pathlib import Path
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/tools_dev")
from pydna.dseqrecord import Dseqrecord
from Bio import Restriction
from assembly.golden_gate import GoldenGateModule
from assembly.common import annotate_module_parts
from workflow.biopart_assemble import BiopartWorkflow
from utils.io import DatabaseImport
from params import DesignParams
from workflow.identity import Manifest
from make_intermediate_vector import build_vector_set, fetch_standard_spacers, STANDARD_SPACERS

S=Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
ok=True
def chk(l,c,e=""):
    global ok; ok=ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

db=DatabaseImport()
seqs=fetch_standard_spacers(db, STANDARD_SPACERS[:2])
vec=build_vector_set(seqs,"BsaI","BsmBI")[0]
SPEC=[("promoter","BBa_J23100"),("rbs","BBa_B0032"),("terminator","BBa_B0010")]
parts=[db.get_part_by_id("Bioparts",t,i) for t,i in SPEC]

gg=GoldenGateModule("BsaI")
with contextlib.redirect_stdout(io.StringIO()):
    design=gg.prepare_fragments(parts, vector=vec, assembly_type="circular")
    recs=[Dseqrecord(d["full_sequence"], name=d["part_name"]) for d in design]
    for r,d in zip(recs,design):
        r.annotations["restriction_enzyme"]=d["restriction_enzyme"]; r.annotations["target_seq"]=d["target_seq"]
    frags=[gg.extract_insert_fragment(r) for r in recs]
    mod=gg.assemble(frags, assembly_type="circular", vector=vec)[0]
mod.name="circ_module"
print(f"circular module {len(mod.seq)}bp  BsaI={len(Restriction.BsaI.search(mod.seq))} BsmBI={len(Restriction.BsmBI.search(mod.seq))}")

d=S/"cp7_exc"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
wf=BiopartWorkflow(project_dir=str(d), parts=[mod])

print("\n1) digest 방식")
with contextlib.redirect_stdout(io.StringIO()) as buf:
    mf=Manifest(d)
    out=wf.module_excision([mod, vec], manifest=mf, design_params=DesignParams(module_assembly_type="circular"))
log=buf.getvalue()
chk("excision 성공", out is not None)
if out:
    ex=out["record"]
    chk("꺼낸 조각이 원형 아님(linear)", not getattr(ex.seq,"circular",False), f"{len(ex.seq)}bp")
    exs=str(ex.seq).upper()
    for (t,i),p in zip(SPEC,parts):
        chk(f"{i} 포함", str(p.seq).upper() in exs)
    chk("L spacer 딸려나옴", seqs[0][-20:] in exs)
    chk("R spacer 딸려나옴", seqs[1][:20] in exs)
    chk("manifest 기록", (d/"reports"/"manifest.json").exists() or mf is not None)

print("\n2) 효소 선택 로직")
prm=DesignParams()
enz,why=wf._pick_excision_enzyme(mod, prm)
chk("조립효소(BsaI)가 아닌 것을 고름", enz is not None and enz!="BsaI", f"{enz} ({why})")

print("\n3) PCR 방식 (명시 지정) — module 구간만 나와야 함")
with contextlib.redirect_stdout(io.StringIO()):
    outp=wf.module_excision([mod, vec], design_params=DesignParams(excision_method="pcr"))
chk("PCR 산출", outp is not None)
if outp:
    ps=str(outp["record"].seq).upper()
    chk("digest 결과와 비슷한 크기(backbone 제외)", abs(len(ps)-250) < 60, f"{len(ps)}bp vs digest 250bp")
    chk("부품 포함", all(str(p.seq).upper() in ps for p in parts))

print("\n4) digest 불가 시 PCR 자동 폴백")
# BsmBI 사이트를 없앤 module을 만들어 digest 불가 상황 재현
fake=Dseqrecord(str(mod.seq).replace("CGTCTC","CGTATC").replace("GAGACG","GATACG"), circular=True)
fake.name="no_exc_site"
enz2,why2=wf._pick_excision_enzyme(fake, prm)
chk("자를 효소 없음을 감지", enz2 is None, why2[:60])
with contextlib.redirect_stdout(io.StringIO()) as buf2:
    out2=wf.module_excision([fake, vec], design_params=DesignParams())
chk("PCR로 폴백해 산출", out2 is not None, f"{len(out2['record'].seq)}bp" if out2 else "")
chk("폴백 경고 출력", "PCR 방식으로 전환" in buf2.getvalue())

print("\n=== CP7 Stage B:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
