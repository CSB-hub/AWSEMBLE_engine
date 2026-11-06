import sys, io, contextlib
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from basic.sequence_analysis import calculate_tm
from basic.primer_design import PrimerDesign
from utils.io import DatabaseImport
from params import DesignParams, DEFAULT_PARAMS
ok=True
def chk(l,c,e=""):
    global ok; ok = ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

seq = "ATGCGTAAAGGAGAAGAACTTTTCACTGGAGTTGTC"
print("1) Tm 계산 조건이 값에 실제로 반영")
t1 = calculate_tm(seq, Na=50.0, dnac1=250.0)
t2 = calculate_tm(seq, Na=100.0, dnac1=500.0)
chk("염/농도 바꾸면 Tm 변함", abs(t1-t2) > 0.3, f"{t1:.2f}℃ → {t2:.2f}℃")

print("\n2) 설계 함수가 params 조건을 사용")
db = DatabaseImport(); part = db.get_part_by_id("Bioparts","promoter","BBa_J23100")
def design(prm):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        return PrimerDesign.design_goldengate_primers(
            insert=part, left_overhang="AAAA", right_overhang="TTTT",
            restriction_enzyme="BsaI", params=prm)
d_def = design(None)
d_salt = design(DesignParams(na_mM=200.0, primer_conc_nM=1000.0))
# 2026-11-05: 기본값이 NEB Q5 조건으로 바뀌었다(Na 0 / 500nM).
chk("params 미지정 = 기본 동작(NEB Q5)",
    d_def["tm_na_mM"] == 0.0 and d_def["tm_primer_nM"] == 500.0)
chk("params 지정 시 조건 기록됨", d_salt["tm_na_mM"] == 200.0, f"Na={d_salt['tm_na_mM']}")
chk("조건 바꾸면 Tm 값도 달라짐", abs(d_def["fwd_tm"] - d_salt["fwd_tm"]) > 0.3,
    f"{d_def['fwd_tm']:.2f} vs {d_salt['fwd_tm']:.2f}")

print("\n3) full / anneal Tm 분리 노출")
chk("fwd_tm_full 존재", "fwd_tm_full" in d_def)
chk("anneal(꼬리 제외) < full(꼬리 포함)", d_def["fwd_tm"] < d_def["fwd_tm_full"],
    f"anneal {d_def['fwd_tm']:.1f} < full {d_def['fwd_tm_full']:.1f}")

print("\n4) 빡빡한 범위를 줘도 설계가 죽지 않는다(하드게이트 아님)")
tight = DesignParams(tm_min=67.0, tm_max=69.0, gc_min=45.0, gc_max=55.0)
d_tight = design(tight)
chk("빡빡한 params로도 설계가 실패하지 않음", d_tight.get("forward"))
# ⚠️ Tm 판정은 2026-09-23부터 유저가 선언한 범위를 그대로 쓴다(설계=판정 단일 창).
#    GC는 여전히 분리 축이다(gc_* 설계 / qc_gc_* 판정).
chk("tier는 여전히 산출됨", "tier" in d_tight, f"tier={d_tight.get('tier')}")

print("\n=== CP3:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
