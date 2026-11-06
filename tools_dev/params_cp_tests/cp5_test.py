"""CP5: build_plasmid에 다른 params를 주면 산출 프라이머가 실제로 달라지는가."""
import sys, shutil, hashlib, io, contextlib, csv
from pathlib import Path
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from params import DesignParams as _DP
# ⚠️ spacer_mode="library"로 고정한다. 2026-10-20에 기본값이 "vector_homology"로
#    바뀌었고 그 모드는 **최종 서열이 달라진다**(중립 spacer 흉터가 사라지므로).
#    이 테스트들은 기존 방식의 **하위호환 보증선**이므로 모드를 명시해 고정한다.
#    새 모드의 기준은 vector_homology_test.py가 따로 든다.
def DesignParams(**kw):
    kw.setdefault('spacer_mode', 'library')
    return _DP(**kw)


S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
db = DatabaseImport()
SPEC=[("promoter","BBa_J23100"),("rbs","BBa_B0032"),("cds","BBa_C0012"),("terminator","BBa_B0010")]
parts=[db.get_part_by_id("Bioparts",t,i) for t,i in SPEC]
vector=db.get_part_by_id("Bioparts","vector","pSB1C3")

def _mk():
    return DesignParams()

def run(tag, prm, parts_override=None):
    d = S/f"cp5_{tag}"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = WorkflowBuilder.build_plasmid(parts_override or parts, vector, db, str(d/"proj"),
                                            module_method="golden_gate", insert_method="gibson",
                                            params=prm)
    rows = list(csv.DictReader(open(d/"proj"/"reports"/"primers.csv", encoding="utf-8")))
    con = res["construct"]
    return {
        "primers": [(r["forward_primer"], r["reverse_primer"]) for r in rows],
        "tm": [(r["fwd_tm"], r["fwd_tm_full"]) for r in rows],
        "enzyme": sorted({r["restriction_enzyme"] for r in rows}),
        "sha": hashlib.sha256(str(con.seq).upper().encode()).hexdigest()[:16],
        "cond": (d/"proj"/"reports"/"design_conditions.md").exists(),
        "vp": int(list(csv.DictReader(open(d/"proj"/"reports"/"vector_primers.csv", encoding="utf-8")))[0]["fwd_mer"]),
        "na": rows[0].get("fwd_tm"),
        "tmfull_ok": all(r.get("fwd_tm_full") for r in rows),
    }

ok=True
def chk(l,c,e=""):
    global ok; ok = ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

print("1) 기본값 = 기존 동작")
base = run("base", DesignParams())
chk("construct SHA가 baseline과 동일", base["sha"] == "d69f339665bad36a", base["sha"])
chk("design_conditions.md 생성", base["cond"])
chk("효소 BsaI", base["enzyme"] == ["BsaI"], str(base["enzyme"]))

print("\n2) Tm 조건(염/농도)을 바꾸면 설계가 달라진다")
salt = run("salt", DesignParams(na_mM=200.0, primer_conc_nM=1000.0))
chk("프라이머 Tm 값이 달라짐", base["tm"] != salt["tm"], f"기본 {base['tm'][0]} vs 변경 {salt['tm'][0]}")
chk("프라이머 서열 자체가 달라짐", base["primers"] != salt["primers"])
chk("fwd_tm_full 채워짐(anneal/full 분리)", base["tmfull_ok"], f"예: {base['tm'][0]}")

print("\n3) 효소를 바꾸면 프라이머가 달라진다")
bsmbi = run("bsmbi", DesignParams(gg_enzyme="BsmBI"))
# BBa_C0012에 BsmBI 사이트가 있어 BsaI로 자동 전환되는 것이 정상 동작
chk("부품 내부 사이트로 안전한 효소 전환", bsmbi["enzyme"] == ["BsaI"], str(bsmbi["enzyme"]))
# 내부 사이트가 없는 부품만으로는 고른 효소가 그대로 쓰여야 한다
clean = run("clean", DesignParams(gg_enzyme="BsmBI"), parts_override=parts[:2]+[parts[3]])
chk("깨끗한 부품이면 고른 효소 사용", clean["enzyme"] == ["BsmBI"], str(clean["enzyme"]))

print("\n4) Gibson homology 길이를 바꾸면 삽입 프라이머가 달라진다")
# ⚠️ 양쪽 arm을 **명시**한다. 기본값에 기대면 기본값이 바뀌는 날(2026-10-24: 20→30)
#    두 쪽이 같아져 "길이를 바꿔도 프라이머가 안 변한다"는 거짓 실패가 난다.
ov20 = run("ov20", DesignParams(gibson_overlap=20))
ov30 = run("ov30", DesignParams(gibson_overlap=30))
# homology는 벡터에서 복제된 서열이라 최종 원형 산물(SHA)은 같은 게 정상 → 프라이머 길이로 판정
chk("vector 프라이머 길이가 +10 (20→30bp)", ov30["vp"] - ov20["vp"] == 10,
    f"{ov20['vp']} → {ov30['vp']}")
chk("construct SHA는 동일(정상)", ov30["sha"] == ov20["sha"])

print("\n=== CP5:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
