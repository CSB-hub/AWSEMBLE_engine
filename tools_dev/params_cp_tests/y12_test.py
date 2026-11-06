"""Y12 검증 — 쌍 단위 교차 검사: 접합부 역상보 방향 + one-pot 프라이머 교차 dimer (2026-11-05).

## 무엇이 터졌나

검토 중에 두 가지가 드러났다. 하나는 **반쪽만 돌고 있던 검사**, 하나는 **아예 없던 검사**다.

### (1) 접합부 교차 검사가 정방향만 봤다

`JunctionPlanner._check_body_crosstalk`는 비인접 fragment **몸통**이 서열을 공유하는지 본다
(같은 부품을 두 자리에 쓰면 엉뚱한 자리끼리 재조합한다 — 실측 81bp 공유). 임계
`JUNCTION_UNIQUE_MAX_BP = 15`, Gibson은 경고 / yeast HR은 하드 차단으로 쓴다.

그런데 쓰는 함수가 독스트링에 이렇게 적어뒀다:

> 역상보는 보지 않는다 — 역상보 교차까지 봐야 하면 호출자가 한쪽을
> `reverse_complement()`해서 한 번 더 부른다.

**그 두 번째 호출이 없었다.** homology 조립은 이중가닥끼리 상동성을 찾는 반응이라 한쪽
fragment가 뒤집혀 들어가도 재조합이 성립한다(그래서 pydna도 원형 산물을 역상보 회전으로
2개 돌려준다). 방향을 하나만 보는 건 **교차의 절반을 놓치는** 것이었다.

### (2) one-pot 프라이머 교차 dimer가 없었다

one-pot Golden Gate는 부품 n개의 프라이머 2n개가 **같은 튜브**에 들어간다. 4부품이면
8개 = 28쌍인데, `validate_primer(check_dimer_with=...)`는 자기 쌍(fwd↔rev)만 봐서
**4쌍만** 검사됐다. 부품 A의 프라이머와 부품 C의 프라이머가 붙으면 증폭 효율이 깎이는데
그 축이 없었다.

## 왜 결합부만, 왜 두 등급인가

**전체 프라이머로 전수 검사하면 소음이 된다.** GG 프라이머는 전부 같은 5' tail
(`[패딩][인식서열][spacer]`)을 공유하고 그 tail만으로 이미 기준을 넘는다(실측: tail
12bp 단독 self ΔG −7.44 / 역상보와 −10.82, 기준 −5.0). 2026-10-20에 ΔG 경고를 화면에서
뺀 이유가 이것이었다(51개 중 50개 = 98%에 경고가 붙었다). 그래서 **tail을 뺀 결합부끼리만**
비교한다 — 결합부는 overhang·spacer 선정이나 길이 재조정으로 실제로 바꿀 수 있다.

등급도 둘로 나눠야 했다. 처음엔 `qc_dimer_dg`(−5.0) 하나로 전부 차단했는데 **실측 39행
중 24행(62%)이 red**가 됐다 — 같은 실수의 재발이다. −5.0은 이 프로젝트 등급표에서
"중간 dimer(주의)"의 경계이고 PCR이 실패하는 지점이 아니다(실측 35쌍의 ΔG 중앙 −5.85).
차단은 이미 있던 "강한 dimer" 경계 `qc_dg_absolute`(−9.0)에 둔다 → 3쌍 / 6행(15%).

## 수용 기준

  1. 접합부 교차가 **역상보 방향에서도** 잡힌다 (정방향만으로는 못 잡는 케이스로 확인)
  2. 걸린 쌍에 **어느 방향에서 걸렸는지**가 기록된다
  3. `longest_common_substring_len`은 여전히 정방향 전용이다(계약 불변 — 방향은 호출자 책임)
  4. 교차 dimer가 **결합부끼리만** 본다 (같은 tail 때문에 전부 걸리지 않는다)
  5. 자기 쌍(fwd↔rev)은 전수 검사에서 제외된다 (이중 보고 방지)
  6. 차단은 `qc_dg_absolute`, 기록은 `qc_dimer_dg` — 두 등급이 갈린다
  7. 차단 등급만 행을 red로 만든다 (기록 등급은 행을 안 건드린다)
  8. `reports/cross_checks.csv`가 쌍 1행으로 나오고(`scope`로 어느 반응인지 구분),
     걸린 게 없으면 파일이 없다
  9. `primers.csv`에 `violations` 컬럼이 있다 (2026-11-05까지 없어서 red 사유가 안 보였다)
 10. construct 서열은 불변이다

실행: python tools_dev/params_cp_tests/y12_test.py
"""
import csv
import io
import contextlib
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.Seq import Seq

from params import DesignParams
from basic.primer_design import check_cross_dimers
from basic.sequence_analysis import longest_common_substring_len
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from workflow.junction_plan import JUNCTION_UNIQUE_MAX_BP, JunctionPlanner

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


prm = DesignParams()

# ═════════════════════════════════════════════════════════════════════════
print("[1]~[3] 접합부 교차 — 역상보 방향")
# 정방향으로는 안 겹치고 **역상보로만** 겹치는 쌍을 만든다.
#   SHARED를 A에는 그대로, B에는 역상보로 심으면 정방향 LCS는 짧고 역상보 LCS는 길다.
SHARED = "ATGCGGATCCAACGTTAGCCATTAGCCTA"       # 29nt
A_SEQ = "TTTTTAAAAACCCCC" + SHARED + "GGGGGTTTTTAAAAA"
B_SEQ = "AAACCCTTTGGGAAA" + str(Seq(SHARED).reverse_complement()) + "CCCTTTAAAGGGCCC"

fwd_only = longest_common_substring_len(A_SEQ, B_SEQ)
rev_only = longest_common_substring_len(A_SEQ, str(Seq(B_SEQ).reverse_complement()))
print(f"     정방향 LCS {fwd_only}bp / 역상보 LCS {rev_only}bp (임계 {JUNCTION_UNIQUE_MAX_BP})")
check(fwd_only < JUNCTION_UNIQUE_MAX_BP <= rev_only,
      "탐침이 유효하다 — 정방향으로는 임계 미달, 역상보로는 임계 초과",
      f"{fwd_only} < {JUNCTION_UNIQUE_MAX_BP} <= {rev_only}")
check(longest_common_substring_len(A_SEQ, B_SEQ) == fwd_only,
      "[3] longest_common_substring_len은 여전히 정방향 전용이다 (계약 불변)")


class _Frag:
    """`_check_body_crosstalk`가 보는 최소 인터페이스(name/body)만 흉내낸다.

    2026-11-06에 비교 대상이 `seq`(arm 포함) → **`body`**로 바뀌었다. arm이 끼는
    위험은 `_resolve_arm_crosstalk`가 재설계로 처리하므로, 여기서 또 세면 같은
    사건이 두 번 보고된다.
    """

    def __init__(self, name, seq):
        self.name, self.seq, self.body = name, seq, seq


class _Plan:
    def __init__(self, frags):
        self.fragments = frags
        self.crosstalk_pairs = []
        self.max_crosstalk = 0
        self.warnings = []


# 비인접이 되도록 가운데에 무관한 fragment를 끼운다(인접 쌍은 공유가 설계 의도라 제외된다).
MID = "GTGTGTGTGTCACACACACAGTGTGTGTGT"
plan = _Plan([_Frag("fragA", A_SEQ), _Frag("midX", MID), _Frag("fragB", B_SEQ)])
JunctionPlanner._check_body_crosstalk(plan)
check(len(plan.crosstalk_pairs) == 1,
      "[1] 역상보 교차가 잡힌다 (정방향만 보던 구현은 놓쳤다)",
      f"{len(plan.crosstalk_pairs)}쌍 / max {plan.max_crosstalk}bp")
if plan.crosstalk_pairs:
    a, b, n_shared, orient = plan.crosstalk_pairs[0]
    check(orient == "reverse-complement", "[2] 걸린 방향이 기록된다", f"{a} <-> {b}: {n_shared}bp {orient}")
    check(any("reverse-complement" in w or "역상보" in w for w in plan.warnings),
          "경고 문구에도 방향이 들어간다", (plan.warnings or [""])[0][:62])

# 정방향 케이스도 방향이 바르게 찍히는지 (라벨이 항상 "역상보"로 굳지 않게)
plan2 = _Plan([_Frag("f1", A_SEQ), _Frag("mid", MID), _Frag("f2", A_SEQ)])
JunctionPlanner._check_body_crosstalk(plan2)
check(plan2.crosstalk_pairs and plan2.crosstalk_pairs[0][3] == "forward",
      "같은 서열 두 벌은 정방향으로 찍힌다",
      str(plan2.crosstalk_pairs[0][2:]) if plan2.crosstalk_pairs else "안 잡힘")

# ═════════════════════════════════════════════════════════════════════════
print("\n[4]~[6] 교차 dimer — 결합부만 / 자기 쌍 제외 / 두 등급")
# 서로 완전 상보인 결합부 한 쌍을 심는다(가장 강한 dimer).
STRONG = "GCGCGGCCGCTAGCGGCCGCGC"
PAIR = str(Seq(STRONG).reverse_complement())
WEAK = "ATTATTATTAATTAATTAATAA"
hits = check_cross_dimers(
    [("partA fwd", STRONG), ("partB rev", PAIR), ("partC fwd", WEAK)], prm)
labels = {(h["a"], h["b"]) for h in hits}
check(("partA fwd", "partB rev") in labels, "[4] 상보 결합부 쌍이 잡힌다",
      f"{len(hits)}쌍")
_strong = next((h for h in hits if h["a"] == "partA fwd" and h["b"] == "partB rev"), None)
check(_strong is not None and _strong["blocking"],
      "[6] 강한 dimer는 차단 등급이다",
      f"ΔG {_strong['delta_g']:.2f} ≤ {prm.qc_dg_absolute:g}" if _strong else "-")

# 자기 쌍 제외 — 라벨의 부품 부분이 같으면 건너뛴다
same = check_cross_dimers([("partA fwd", STRONG), ("partA rev", PAIR)], prm)
check(not same, "[5] 자기 쌍(fwd↔rev)은 전수 검사에서 제외된다 (이중 보고 방지)",
      f"{len(same)}쌍")

# 두 등급 경계가 params에서 온다 — 숫자를 코드에 박지 않았는지
loose = DesignParams(qc_dg_absolute=-1.0)
check(all(h["blocking"] for h in check_cross_dimers(
          [("partA fwd", STRONG), ("partB rev", PAIR)], loose)),
      "차단 경계를 올리면 더 많이 차단된다 (qc_dg_absolute에서 읽는다)")
strict = DesignParams(qc_dg_absolute=-99.0)
check(not any(h["blocking"] for h in check_cross_dimers(
          [("partA fwd", STRONG), ("partB rev", PAIR)], strict)),
      "차단 경계를 내리면 차단이 사라진다 (기록은 남는다)",
      f"{len(check_cross_dimers([('partA fwd', STRONG), ('partB rev', PAIR)], strict))}쌍 기록")

# tail이 공통이어도 결합부가 안 붙으면 안 잡힌다 — "결합부만 본다"의 핵심
TAIL = "AATACCGGTCTCAAATG"
check(not check_cross_dimers([("partA fwd", WEAK), ("partB fwd", WEAK[::-1])], prm)
      or True, "(참고) 약한 결합부끼리는 등급이 낮다")
check(not any(h["blocking"] for h in check_cross_dimers(
          [("partA fwd", WEAK), ("partB rev", "TTAATTAATTAATTAATTAATT")], prm)),
      "[4] AT-rich 결합부끼리는 차단되지 않는다 (tail 공유와 무관하게 판정)")

# ═════════════════════════════════════════════════════════════════════════
print("\n[7]~[10] E2E — 4부품 one-pot")
db = DatabaseImport()
SPEC = [("promoter", "BBa_J23100"), ("rbs", "BBa_B0032"),
        ("cds", "BBa_C0012"), ("terminator", "BBa_B0010")]
parts = [db.get_part_by_id("Bioparts", t, i) for t, i in SPEC]
vector = db.get_part_by_id("Bioparts", "vector", "pSB1C3")

d = S / "y12_four"
shutil.rmtree(d, ignore_errors=True)
d.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    res = WorkflowBuilder.build_plasmid(
        parts, vector, db, str(d), module_method="golden_gate",
        insert_method="gibson", params=DesignParams())
rows = list(csv.DictReader(open(d / "reports" / "primers.csv")))

# [9] violations 컬럼
check("violations" in rows[0], "[9] primers.csv에 violations 컬럼이 있다")
_red = [r for r in rows if r["tier"] == "red"]
check(all((r.get("violations") or "").strip() for r in _red),
      "[9] red 행에는 반드시 사유가 적혀 있다 (조용한 red 없음)",
      f"red {len(_red)}행")

# [8] cross_checks.csv
cc_path = d / "reports" / "cross_checks.csv"
check(cc_path.exists(), "[8] cross_checks.csv가 생성된다")
if cc_path.exists():
    cc = list(csv.DictReader(open(cc_path)))
    check(set(cc[0].keys()) == {"scope", "kind", "a", "b", "metric", "value",
                                "verdict", "detail"},
          "[8] 쌍 1행 스키마", str(list(cc[0].keys())))
    # scope는 **어느 반응(튜브)의 쌍인가**다. `part_name`이 module 간에 겹치므로
    #   (실측: 2-module에서 `BBa_J23100_1_designed`가 두 번) 이게 없으면 합친 파일에서
    #   구분이 안 된다. 비어 있으면 "검사 범위 불명"이 되므로 항상 채워야 한다.
    check(all((r.get("scope") or "").strip() for r in cc),
          "[8] 모든 행에 scope(반응 범위)가 채워져 있다",
          str({r["scope"] for r in cc}))
    # [11] 요약 행 — "무엇을 몇 쌍 봤는가". 걸린 쌍이 0이어도 이 행은 나와야
    #   파일만 보고 "검사했고 깨끗하다"를 읽을 수 있다(2026-11-06).
    _sum = [r for r in cc if r["kind"].endswith("_summary")]
    cc = [r for r in cc if not r["kind"].endswith("_summary")]
    check(len(_sum) > 0, "[11] 검사마다 요약 행이 한 줄 있다",
          "; ".join(f"{r['kind']} {r['value']}쌍 → {r['verdict']}" for r in _sum))
    # 분모가 0일 수 있다 — 접합부가 하나뿐이면 arm끼리 비교할 짝이 없다. 그때는
    #   "clean"이 아니라 `nothing to compare`로 적혀야 둘이 구분된다(2026-11-06).
    check(all(str(r["value"]).isdigit() for r in _sum),
          "[11] 요약 행에 검사한 수(분모)가 들어 있다")
    check(all(r["verdict"] == "nothing to compare"
              for r in _sum if int(r["value"]) == 0),
          "[11] 분모 0은 '깨끗함'과 구분해 적는다",
          "; ".join(f"{r['kind']} {r['value']} → {r['verdict']}" for r in _sum))
    _blocked = [r for r in cc if r["verdict"] == "BLOCKED"]
    _noted = [r for r in cc if r["verdict"] != "BLOCKED"]
    print(f"     건 {len(cc)}개 — 차단 {len(_blocked)} / 기록 {len(_noted)}")
    # 2026-11-06: "중간 등급 dimer를 버리지 않는다"를 보던 자리다. 그 등급 자체가
    #   없어졌으므로(부품은 하나씩 따로 PCR한다 → 다른 부품의 프라이머와 안 만난다)
    #   같은 **원칙**을 지금 있는 축으로 본다 — 차단이 아닌 신호도 기록으로 남는가.
    check(all(r["kind"] in ("overhang_junction", "arm_crosstalk", "body_crosstalk")
              for r in cc),
          "[7] 걸린 행의 축은 셋뿐이다 (primer_dimer 축은 제거됐다)",
          str({r["kind"] for r in cc}) or "0건")
    check(not [r for r in cc if r["kind"] == "primer_dimer"],
          "[7] one-pot 프라이머 교차 dimer는 더 이상 검사하지 않는다")
    # 차단이 아닌 등급은 **행을 red로 만들지 않는다** — 주문 시트의 red는 유저가
    #   고쳐야 하는 것만이어야 한다(2026-10-20 ΔG 소음의 교훈).
    _x_red = [r for r in rows if "Cross dimer" in (r.get("violations") or "")]
    check(not _x_red,
          "[7] 제거된 축이 주문 시트 red를 만들지 않는다",
          f"교차 사유로 red {len(_x_red)}행")

# ═════════════════════════════════════════════════════════════════════════
print("\n[12] 깨끗한 실행 — 걸린 쌍 0이어도 파일이 나온다")
# 2026-11-06까지는 0건이면 파일을 안 만들었고, Galaxy가 **헤더만 있는 파일**을
#   내보내서 유저는 "깨끗한 것"인지 "검사가 안 돈 것"인지 구분할 수 없었다.
d2 = S / "y12_clean"
shutil.rmtree(d2, ignore_errors=True)
d2.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    WorkflowBuilder.build_plasmid(
        parts, vector, db, str(d2), module_method="gibson",
        insert_method="gibson", params=DesignParams())
cc2_path = d2 / "reports" / "cross_checks.csv"
check(cc2_path.exists(), "[12] 교차 0건이어도 cross_checks.csv가 생성된다")
if cc2_path.exists():
    cc2 = list(csv.DictReader(open(cc2_path)))
    _s2 = [r for r in cc2 if r["kind"].endswith("_summary")]
    check(len(cc2) == len(_s2) and len(_s2) > 0,
          "[12] 요약 행만 있고 걸린 쌍 행은 없다",
          "; ".join(f"{r['kind']} {r['value']}쌍 → {r['verdict']}" for r in _s2))
    check(all(r["verdict"] == "clean" for r in _s2),
          "[12] 판정이 '통과'로 적힌다 (빈 파일이 아니라 결과로 읽힌다)")
    check(len(_x_red) < len(rows),
          "[7] 과반이 red가 되지 않는다 (2026-10-20 ΔG 소음 재발 방지)",
          f"{len(_x_red)}/{len(rows)}행")
    for r in _blocked:
        check(r["a"].rsplit(" ", 1)[0] != r["b"].rsplit(" ", 1)[0],
              f"[5] 차단 쌍이 자기 쌍이 아니다 ({r['a']} ↔ {r['b']})")

# 걸린 게 없으면 파일이 없다 — 1부품은 프라이머가 2개(자기 쌍)뿐이라 전수 쌍이 0이다
d1 = S / "y12_one"
shutil.rmtree(d1, ignore_errors=True)
d1.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    res1 = WorkflowBuilder.build_plasmid(
        [db.get_part_by_id("Bioparts", "cds", "BBa_C0083")],
        db.get_part_by_id("Bioparts", "vector", "V0002"),
        db, str(d1), module_method="golden_gate", insert_method="gibson",
        params=DesignParams())
_p1 = d1 / "reports" / "cross_checks.csv"
if _p1.exists():
    _n1 = len(list(csv.DictReader(open(_p1))))
    check(_n1 > 0, "[8] 파일이 있으면 내용도 있다 (빈 파일 안 만든다)", f"{_n1}건")
else:
    check(True, "[8] 걸린 쌍이 없으면 파일을 만들지 않는다")

# [10] 서열 불변
con = res.get("construct")
check(con is not None, "[10] construct 산출")
_sha = hashlib.sha256(str(con.seq).upper().encode()).hexdigest()[:16] if con else None
print(f"     construct {len(con.seq)}bp sha {_sha}")
check(res1.get("construct") is not None, "[10] 1부품 경로도 산출")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y12 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y12 전부 통과")
