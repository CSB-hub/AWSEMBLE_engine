"""Y3 검증 — arm 축이 **삽입 방법**을 따라가는가 (2026-10-24에 실제로 터진 버그).

## 무엇이 터졌나

`spacer_mode="vector_homology"`는 module 양 끝에 벡터 절단부 서열을 미리 박아 넣는다.
그 서열이 **나중에 삽입 조립의 overlap**이 된다. 그런데 `VectorHomologySpacer`가
`params.gibson_overlap`(30bp)을 직접 읽고 있었고, yeast 삽입은 `yeast_homology`(45bp)를
요구하므로 pydna가 45bp 겹침을 못 찾아 **"Assembly 실패: 결과물이 생성되지 않았습니다"**로
떨어졌다. 조용한 실패여서 원인이 바로 안 보였다.

## 수용 기준

  1. `_insert_arm()`이 삽입 방법별로 맞는 arm을 돌려준다
  2. Golden Gate 삽입은 homology를 안 쓰므로 `gibson_overlap`으로 폴백한다
  3. `VectorHomologySpacer(arm=…)`가 그 길이로 spacer를 만든다
  4. arm 미지정 시 기존 동작(`gibson_overlap`)이 유지된다 — 하위호환
  5. module 양 끝이 **삽입 방법이 요구하는 길이만큼** 벡터 homology를 들고 있다
실행: python tools_dev/params_cp_tests/y3_test.py
"""
import sys

sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")

from utils.io import DatabaseImport
from params import DesignParams
from workflow.builder import WorkflowBuilder
from workflow.vector_homology_spacer import VectorHomologySpacer, module_carries_homology

FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


print("=" * 72)
print("Y3 — arm 축이 삽입 방법을 따라가는가")
print("=" * 72)

prm = DesignParams()
prm.validate()

print(f"\n[1] _insert_arm() — gibson_overlap={prm.gibson_overlap} / "
      f"yeast_homology={prm.yeast_homology_min}~{prm.yeast_homology_max}")
for method, expect in (("gibson", prm.gibson_overlap),
                       ("yeast_hr", prm.yeast_homology_min),
                       ("golden_gate", prm.gibson_overlap)):
    got = WorkflowBuilder._insert_arm(method, prm)
    print(f"  {method:<12} → {got}bp (기대 {expect})")
    check(got == expect, f"{method} arm = {expect}")

print("\n[2] 알 수 없는 방법은 gibson_overlap으로 폴백")
check(WorkflowBuilder._insert_arm("nonexistent_method", prm) == prm.gibson_overlap,
      "미지원 방법 → gibson_overlap 폴백")

print("\n[3] VectorHomologySpacer(arm=…)가 그 길이로 만든다")
db = DatabaseImport()
vector = db.get_part_by_id("Bioparts", "vector", "pSB1C3")
for arm in (20, 30, 45):
    l, r = VectorHomologySpacer(params=prm, arm=arm).make(vector)
    print(f"  arm={arm:>3} → L {len(l.seq)}bp / R {len(r.seq)}bp")
    check(len(l.seq) == arm and len(r.seq) == arm, f"arm={arm} spacer 길이 일치")

print("\n[4] arm 미지정 = 기존 동작 (하위호환)")
l0, r0 = VectorHomologySpacer(params=prm).make(vector)
check(len(l0.seq) == prm.gibson_overlap,
      f"미지정 → gibson_overlap({prm.gibson_overlap}bp)", f"{len(l0.seq)}bp")

print("\n[5] 삽입 방법별로 module 양 끝이 요구 길이를 채우는가")
SPEC = [("promoter", "BBa_J23100"), ("rbs", "BBa_B0032"),
        ("cds", "BBa_C0012"), ("terminator", "BBa_B0010")]
parts = [db.get_part_by_id("Bioparts", t, i) for t, i in SPEC]

for insert_method in ("gibson", "yeast_hr"):
    arm = WorkflowBuilder._insert_arm(insert_method, prm)
    parts_sp = VectorHomologySpacer(params=prm, arm=arm).attach(parts, vector)
    lin = VectorHomologySpacer(params=prm, arm=arm).linearize(vector)
    # attach 결과를 이어붙인 것이 module 본체가 된다(조립 전 서열 기준 점검).
    joined = "".join(str(p.seq) for p in parts_sp)
    from pydna.dseqrecord import Dseqrecord
    ok = module_carries_homology(Dseqrecord(joined), lin, arm)
    print(f"  {insert_method:<10} arm {arm}bp → module 양 끝 homology {ok}")
    check(ok, f"{insert_method}: module 양 끝이 {arm}bp homology를 들고 있음")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y3 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y3 전부 통과")
