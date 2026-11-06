"""Y4 검증 — 조립 효소를 **한 곳에서** 정하는가 + homology 삽입이 벡터를 선형화하는가
(2026-10-25에 실제로 터진 버그 두 건).

## 무엇이 터졌나

**(A) yeast HR 삽입이 100% 실패했다.** `WorkflowBase.assemble`의 분기가
`method == "gibson"`에만 걸려 있어서, yeast_hr은 그 블록을 통째로 건너뛰고 맨 아래
`else`로 떨어졌다. 거기서는 **선형화하지 않은 원형 벡터**를 그대로 pydna에 넘긴다.
원형 벡터는 말단이 없어 homology를 물릴 자리가 없으므로, 예외 없이 빈 결과가 나오고
"Assembly 실패: 결과물이 생성되지 않았습니다"로만 보였다.

**(B) Golden Gate 효소가 세 군데서 따로 정해졌다.** `params.gg_enzyme`(유저 선택),
`prepare_fragments`의 자동 전환(**부품만** 봄), 벡터 삽입(**BsaI 하드코딩**). one-pot이라
전부 같은 효소여야 하는데 서로 모르는 채 각자 정해서, DB 실측 기준 될 수 있었던 설계의
32.7%p를 놓치고 있었다(GG 가능 벡터 202개 × 무작위 4-part module 3000회:
공통 효소 존재 68.7% vs 실제 성공 36.0%).

## 수용 기준

  1. homology 기반 assembler는 `ARM_PARAMS`를 갖고, Golden Gate는 안 갖는다
     (`assemble`의 분기가 이 훅으로 갈린다 — 방법 이름으로 가르면 (A)가 재발한다)
  2. `select_gg_enzyme`이 **벡터 조건(2곳 이상)과 부품 조건(0개)의 교집합**을 쓴다
  3. 유저가 고른 효소가 교집합에 있으면 **바꾸지 않는다**
  4. 교집합에 없으면 전환하고 `switched=True` + 사유를 남긴다
  5. 교집합이 비면 `enzyme=None` + **부품 쪽을 가리키는** 안내를 준다
  6. `vector=None`(part→module만)이면 벡터 조건을 걸지 않는다
  7. 벡터 삽입이 BsaI를 하드코딩하지 않는다 — 결정된 효소를 읽는다

실행: python tools_dev/params_cp_tests/y4_test.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from pydna.dseqrecord import Dseqrecord

from assembly.factory import AssemblyFactory
from assembly.golden_gate import VECTOR_MIN_SITES, select_gg_enzyme
from params import DesignParams, SUPPORTED_ENZYMES

FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


def spaced(site, n, pad=600):
    """인식서열 `site`를 `n`번 가진 합성 서열. 사이는 GC가 치우치지 않게 채운다."""
    filler = "ACGTTGCAATCGGATCCATGCAGTTACCGGATTACGCATGACCTTAAGCTGACGTACGTT"
    out = []
    for i in range(n):
        out.append(filler[: pad // max(n, 1)])
        out.append(site)
    out.append(filler[: pad // max(n, 1)])
    return "".join(out)


print("=" * 72)
print("Y4 — 효소 단일 결정 + homology 삽입의 벡터 선형화")
print("=" * 72)

prm = DesignParams()

print("\n[1] 분기 훅 — homology 기반만 ARM_PARAMS를 갖는다")
for method, want in (("gibson", True), ("yeast_hr", True), ("golden_gate", False)):
    a = AssemblyFactory.create_assembler(method, params=prm)
    has = bool(getattr(a, "ARM_PARAMS", None))
    print(f"  {method:<12} ARM_PARAMS={getattr(a, 'ARM_PARAMS', None)}")
    check(has == want, f"{method}: ARM_PARAMS 보유 == {want}")

print("\n[2] 벡터 선형화 분기가 '방법 이름'이 아니라 훅으로 갈린다")
src = (ROOT / "src" / "workflow" / "base.py").read_text(encoding="utf-8")
check('if getattr(assembler, "ARM_PARAMS", None):' in src,
      "assemble()이 ARM_PARAMS 훅으로 분기한다")
check(not re.search(r'if method\.lower\(\) == "gibson":', src),
      "gibson 이름으로만 거는 분기가 남아 있지 않다")

print("\n[3] 벡터 삽입이 BsaI를 하드코딩하지 않는다")
code = re.sub(r"#.*", "", src)          # 주석의 BsaI 언급은 제외
check("BsaI" not in code, "base.py 실행 코드에 BsaI 리터럴 없음",
      ", ".join(l.strip()[:40] for l in code.splitlines() if "BsaI" in l))

print(f"\n[4] 교집합 — 벡터는 {VECTOR_MIN_SITES}곳 이상 / 부품은 0개")
# BsaI 2곳 + BbsI 2곳을 가진 벡터, BsaI를 내부에 가진 부품 → BbsI만 남아야 한다
vector = Dseqrecord(spaced("GGTCTC", 2) + spaced("GAAGAC", 2), name="v_test", circular=True)
part_clean = Dseqrecord(spaced("TTTTTT", 0), name="clean")
part_bsai = Dseqrecord(spaced("GGTCTC", 1), name="has_BsaI")

r = select_gg_enzyme([part_clean, part_bsai], vector, preferred="BsaI")
print(f"  벡터사이트={r['vector_sites']}  후보={r['candidates']}  선택={r['enzyme']}")
check(r["enzyme"] == "BbsI", "BsaI가 부품에 있으면 BbsI로 전환", str(r["enzyme"]))
check(r["switched"], "전환을 switched=True로 알린다")
check("has_BsaI" in r["reason"], "사유가 막는 부품 이름을 가리킨다", r["reason"][:60])
check(r["blocked_by"]["BsaI"] == ["has_BsaI"], "blocked_by가 부품을 기록한다")

print("\n[5] 유저 선택이 교집합에 있으면 바꾸지 않는다")
r = select_gg_enzyme([part_clean], vector, preferred="BbsI")
check(r["enzyme"] == "BbsI" and not r["switched"], "BbsI 유지", str(r["enzyme"]))
r = select_gg_enzyme([part_clean], vector, preferred="BsaI")
check(r["enzyme"] == "BsaI" and not r["switched"], "BsaI 유지(부품에 없으므로)", str(r["enzyme"]))

print("\n[6] 벡터가 못 여는 효소는 후보에서 빠진다")
# BsaI만 2곳인 벡터 — 부품이 깨끗해도 BbsI/BsmBI는 벡터를 못 연다
v_bsai = Dseqrecord(spaced("GGTCTC", 2), name="v_bsai_only", circular=True)
r = select_gg_enzyme([part_clean], v_bsai, preferred="BbsI")
print(f"  벡터사이트={r['vector_sites']}  후보={r['candidates']}  선택={r['enzyme']}")
check(r["candidates"] == ["BsaI"], "BsaI만 후보", str(r["candidates"]))
check(r["enzyme"] == "BsaI" and r["switched"], "BbsI → BsaI 전환")
check("0곳만 잘림" in r["reason"] or "곳만 잘림" in r["reason"],
      "사유가 벡터 절단 횟수를 말한다", r["reason"][:70])

print("\n[7] 1곳만 자르는 벡터는 Golden Gate 불가 (좌우 overhang이 같아짐)")
v_one = Dseqrecord(spaced("GGTCTC", 1), name="v_one_site", circular=True)
r = select_gg_enzyme([part_clean], v_one, preferred="BsaI")
check(r["enzyme"] is None, "효소 없음", str(r["enzyme"]))
check("불가능" in r["reason"], "불가능하다고 명시", r["reason"][:60])

print("\n[8] 교집합이 비면 고칠 수 있는 쪽(부품)을 가리킨다")
r = select_gg_enzyme([Dseqrecord(spaced("GGTCTC", 1), name="p_bad")], v_bsai, preferred="BsaI")
check(r["enzyme"] is None, "효소 없음")
check("p_bad" in r["reason"], "막는 부품 이름이 안내에 들어간다", r["reason"][:80])
check("도메스티케이션" in r["reason"], "도메스티케이션을 권고한다")

print("\n[9] vector=None이면 부품 조건만 본다 (part→module 경로)")
r = select_gg_enzyme([part_bsai], None, preferred="BsaI")
check(r["enzyme"] in [e for e in SUPPORTED_ENZYMES if e != "BsaI"],
      "벡터 없이도 부품 충돌을 피해 고른다", str(r["enzyme"]))
check(r["vector_sites"] == {}, "벡터 조건을 걸지 않았다")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y4 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y4 전부 통과")
