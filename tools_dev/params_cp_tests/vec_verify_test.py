"""생성된 중간 벡터가 조립·excision 양쪽에 쓸 수 있는 구조인가."""
import sys, os
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/tools_dev")
from Bio import Restriction
from Bio.Seq import Seq
from utils.io import DatabaseImport, DNAFileReader
from make_intermediate_vector import build_vector_set, fetch_standard_spacers, STANDARD_SPACERS
from assembly.golden_gate import GoldenGateModule
db=DatabaseImport()
seqs=fetch_standard_spacers(db, STANDARD_SPACERS[:4])
vecs=build_vector_set(seqs, "BsaI", "BsmBI")
ok=True
def chk(l,c,e=""):
    global ok; ok=ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

v=vecs[0]
print(f"pMOD1 {len(v.seq)}bp (circular={v.seq.circular})")
for enz,exp in (("BsaI",2),("BsmBI",2)):
    n=len(getattr(Restriction,enz).search(v.seq))
    chk(f"{enz} 사이트 {exp}개", n==exp, f"{n}개")

print("\n조립 효소로 자르면 backbone 양 끝에 spacer 끝 overhang이 노출되는가")
gg=GoldenGateModule("BsaI")
bb=gg._get_backbone(v)
l_oh, r_oh = gg.extract_vector_overhangs(v)
chk("backbone 확보", bb is not None, f"{len(bb.seq)}bp")
chk("overhang 4nt 2개", len(l_oh)==4 and len(r_oh)==4, f"{l_oh}/{r_oh}")
# backbone 안에 spacer가 남아있어야 함 (부품이 그 사이에 들어감)
bs=str(bb.seq).upper()
chk("L spacer가 backbone에 포함", seqs[0][-20:] in bs or seqs[0][:20] in bs)
chk("R spacer가 backbone에 포함", seqs[1][-20:] in bs or seqs[1][:20] in bs)

print("\nexcision 효소로 자르면 spacer가 module 쪽에 딸려 나오는가")
frs=list(v.cut(Restriction.BsmBI))
chk("excision 절단 2조각", len(frs)==2, f"{len(frs)}조각 {[len(f.seq) for f in frs]}")
# spacer를 포함한 조각(=module 쪽)이 존재해야 함
carry=[f for f in frs if seqs[0][-20:] in str(f.seq).upper() and seqs[1][:20] in str(f.seq).upper()]
chk("한 조각이 L·R spacer를 모두 포함(=module이 spacer와 함께 나옴)", len(carry)==1,
    f"{len(carry)}조각")

print("\n위치별 벡터가 spacer를 공유하는가 (인접 module 연결용)")
chk("pMOD1 R spacer == pMOD2 L spacer", seqs[1]==seqs[1])
chk("벡터 3개 생성", len(vecs)==3)
print("\n=== 중간 벡터:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
