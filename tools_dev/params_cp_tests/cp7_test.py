"""CP7 Stage A: circular module 조립이 벡터 backbone을 포함해 성립하는가."""
import sys, io, contextlib, hashlib
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from pydna.dseqrecord import Dseqrecord
from Bio.Seq import Seq
from Bio.Restriction import BsaI
from assembly.golden_gate import GoldenGateModule, _rc
from basic.primer_design import PrimerDesign
from utils.io import DatabaseImport
ok=True
def chk(l,c,e=""):
    global ok; ok = ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

# 합성 중간 벡터 (BsaI 사이트가 stuffer를 향함)
OH_L,OH_R="AATG","GCTT"; site="GGTCTC"; sp="A"
bbA="TTTTTGCATCTAGACCCGGGATCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCC"
bbB="CCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAAC"
stuffer="GGGGTACCCCTTTTGGGGAAAACCCCTTTTGGGG"
vec=Dseqrecord(bbA+OH_L+_rc(sp)+_rc(site)+stuffer+site+sp+OH_R+bbB, name="destA", circular=True)

db=DatabaseImport()
parts=[db.get_part_by_id("Bioparts","promoter","BBa_J23100"),
       db.get_part_by_id("Bioparts","rbs","BBa_B0032"),
       db.get_part_by_id("Bioparts","terminator","BBa_B0010")]
gg=GoldenGateModule("BsaI")

def build(assembly_type, vector=None):
    with contextlib.redirect_stdout(io.StringIO()):
        design=gg.prepare_fragments(parts, vector=vector, assembly_type=assembly_type)
        recs=[Dseqrecord(d["full_sequence"], name=d["part_name"]) for d in design]
        for r,d in zip(recs,design):
            r.annotations["restriction_enzyme"]=d["restriction_enzyme"]
            r.annotations["target_seq"]=d["target_seq"]
        frags=[gg.extract_insert_fragment(r) for r in recs]
        return gg.assemble(frags, assembly_type=assembly_type, vector=vector)[0], frags

print("1) circular 조립 — 벡터 backbone 포함")
mod_c, frags_c = build("circular", vec)
chk("circular 산출", mod_c is not None)
chk("결과가 원형", getattr(mod_c.seq,"circular",False), f"{len(mod_c.seq)}bp")
bb=gg._get_backbone(vec)
exp=len(bb.seq)+sum(len(f.seq) for f in frags_c)-4*(len(frags_c)+1)
chk("길이가 예상과 일치", len(mod_c.seq)==exp, f"{len(mod_c.seq)}bp (예상 {exp})")

print("\n2) 부품 서열이 온전히 보존되는가 (원형이라 origin 넘김 고려)")
cs=(str(mod_c.seq).upper())*2
for p in parts:
    ps=str(p.seq).upper()
    chk(f"{p.name} 보존", ps in cs or _rc(ps) in cs)

print("\n3) 벡터 backbone도 포함되는가")
chk("backbone 서열 포함", str(bb.seq).upper()[10:40] in cs or _rc(str(bb.seq).upper())[10:40] in cs)

print("\n4) 조립 효소 사이트가 결과에 남지 않는가")
n=len(BsaI.search(Seq(str(mod_c.seq))))
chk("BsaI 잔류 0", n==0, f"{n}개")

print("\n5) linear 경로 회귀 (기존 동작 불변)")
mod_l,_=build("linear", None)
chk("linear 정상 산출", mod_l is not None, f"{len(mod_l.seq)}bp")
chk("linear는 원형 아님", not getattr(mod_l.seq,"circular",False))

print("\n6) vector 없이 circular 요청 → 명확한 에러")
try:
    build("circular", None); chk("에러 발생", False)
except ValueError as e:
    chk("에러 발생", True, str(e)[:60])

print("\n=== CP7 Stage A:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
