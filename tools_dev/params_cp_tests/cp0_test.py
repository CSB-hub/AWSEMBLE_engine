import sys
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from pydna.dseqrecord import Dseqrecord
from assembly.golden_gate import GoldenGateModule
from utils.io import DatabaseImport

gg = GoldenGateModule()
ok = True
def case(label, fn, expect_error):
    global ok
    try:
        r = fn(); got = ("정상", r)
    except ValueError as e:
        got = ("ValueError", str(e))
    except Exception as e:
        got = (type(e).__name__, str(e))
    good = (got[0] == "ValueError") if expect_error else (got[0] == "정상")
    ok = ok and good
    print(f"  [{'OK ' if good else 'FAIL'}] {label}")
    print(f"        → {got[0]}: {str(got[1])[:200]}")

print("1) BsaI 사이트 없는 vector (회귀 대상 케이스)")
no_site = Dseqrecord("ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC" * 5, name="no_bsai_vector", circular=True)
case("사이트 없음 → 명확한 ValueError", lambda: gg.extract_vector_overhangs(no_site), True)

print("\n2) 좌우 overhang이 같은 벡터 (방향성 클로닝 불가)")
# ⚠️ 2026-09-23: 벡터를 **원형으로 정규화**하면서 이 케이스의 성격이 바뀌었다.
#    구버전은 선형 절단이라 backbone이 4bp 미만이 되어 길이 가드에 걸렸다.
#    원형 절단에서는 길이가 충분해지는 대신 **좌우 overhang이 동일**해진다 —
#    사이트가 1개뿐인 벡터의 전형적 증상이고, 이게 진짜 막아야 할 조건이다.
tiny = Dseqrecord("GGTCTCAG", name="tiny_vector")
case("좌우 overhang 동일 → ValueError", lambda: gg.extract_vector_overhangs(tiny), True)

db = DatabaseImport()
v1 = db.get_part_by_id("Bioparts", "vector", "pSB1A2")   # BsaI 1개
case("pSB1A2(BsaI 1개) → 좌우 동일이므로 차단", lambda: gg.extract_vector_overhangs(v1), True)

print("\n3) 실제 vector — 동작 불변 확인")
v = db.get_part_by_id("Bioparts", "vector", "V0001")      # BsaI 2개 = GG 삽입 가능 조건
case("V0001(BsaI 2개) 정상 추출", lambda: gg.extract_vector_overhangs(v), False)

print("\n=== CP0:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
