import sys, io, contextlib, hashlib
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from pydna.dseqrecord import Dseqrecord
from assembly.golden_gate import GoldenGateModule
from basic.primer_design import PrimerDesign
from utils.io import DatabaseImport
from params import DesignParams, THERMO_MAX_LEN
db = DatabaseImport()
part = db.get_part_by_id("Bioparts","promoter","BBa_J23100")
ok=True
def chk(l,c,e=""):
    global ok; ok = ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

def design(prm, enz="BsaI"):
    with contextlib.redirect_stdout(io.StringIO()):
        return PrimerDesign.design_goldengate_primers(insert=part, left_overhang="AATG",
            right_overhang="GCTT", restriction_enzyme=enz, params=prm)

print("1) 패딩이 실제로 프라이머에 붙는가")
d0 = design(DesignParams(site_5prime_pad=0))
d6 = design(DesignParams())          # 기본 = 효소별 표 = 6bp
chk("패딩 없음: 인식서열이 0번 위치", d0["forward"].index("GGTCTC")==0, d0["forward"][:14])
chk("기본값: 인식서열 앞 6bp", d6["forward"].index("GGTCTC")==6, d6["forward"][:16])
chk("reverse도 패딩으로 시작", d6["reverse"].index("GGTCTC")==6, d6["reverse"][:16])
chk("패딩 서열 기록됨", d6["site_5prime_pad"]=="AATACC", d6["site_5prime_pad"])

print("\n2) 절단 후 fragment가 패딩과 무관하게 동일한가 (조립 산물 불변)")
gg = GoldenGateModule("BsaI")
shas={}
for tag, dd in (("pad 0bp", d0), ("pad 6bp", d6)):
    rec = Dseqrecord(dd["full_sequence"]); rec.name="t"
    rec.annotations["restriction_enzyme"]="BsaI"; rec.annotations["target_seq"]=dd["target_seq"]
    with contextlib.redirect_stdout(io.StringIO()):
        f = gg.extract_insert_fragment(rec)
    shas[tag]=(len(f.seq), hashlib.sha256(str(f.seq).upper().encode()).hexdigest()[:12])
    print(f"      {tag}: full {len(dd['full_sequence'])}bp → fragment {len(f.seq)}bp  {shas[tag][1]}")
chk("fragment 길이·서열 동일", shas["pad 0bp"]==shas["pad 6bp"])

print("\n3) 효소를 바꿔도 성립하는가")
for enz in ("BbsI","BsmBI"):
    dd = design(DesignParams(gg_enzyme=enz), enz)
    site = {"BbsI":"GAAGAC","BsmBI":"CGTCTC"}[enz]
    rec = Dseqrecord(dd["full_sequence"]); rec.name="t"
    rec.annotations["restriction_enzyme"]=enz; rec.annotations["target_seq"]=dd["target_seq"]
    g2 = GoldenGateModule(enz)
    with contextlib.redirect_stdout(io.StringIO()):
        f = g2.extract_insert_fragment(rec)
    d_no = design(DesignParams(gg_enzyme=enz, site_5prime_pad=0), enz)
    rec2 = Dseqrecord(d_no["full_sequence"]); rec2.name="t"
    rec2.annotations["restriction_enzyme"]=enz; rec2.annotations["target_seq"]=d_no["target_seq"]
    with contextlib.redirect_stdout(io.StringIO()):
        f2 = g2.extract_insert_fragment(rec2)
    chk(f"{enz}: 패딩 6bp 부착 + fragment 불변",
        dd["forward"].index(site)==6 and str(f.seq)==str(f2.seq), f"fragment {len(f.seq)}bp")

print("\n4) 길이 예산")
for tt in (60, 65, 68):
    dd = design(DesignParams(tm_min=float(tt)-5, tm_max=float(tt)+5))
    chk(f"Tm 중앙 {tt}: 프라이머 {len(dd['forward'])}bp ≤ {THERMO_MAX_LEN}",
        len(dd["forward"])<=THERMO_MAX_LEN)

print("\n5) 패딩이 새 제한효소 사이트를 만들지 않는가")
from Bio import Restriction; from Bio.Seq import Seq
from params import SUPPORTED_ENZYMES
sq=Seq(d6["full_sequence"])
cnt={e:len(getattr(Restriction,e).search(sq)) for e in SUPPORTED_ENZYMES}
chk("BsaI 2개(설계분)·나머지 0개", cnt["BsaI"]==2 and cnt["BbsI"]==0 and cnt["BsmBI"]==0, str(cnt))

print("\n=== 패딩:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
