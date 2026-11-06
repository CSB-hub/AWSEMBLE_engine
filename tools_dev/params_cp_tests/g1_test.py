"""G1 검증 — 접합부 설계기 (병합 폐기 → tail/올리고 2갈래, 길이 탐색, padding 없음).

핵심 불변식
    **PCR 템플릿은 실제로 존재하는 분자 하나여야 한다.**
    템플릿에 없는 서열(arm, 흡수된 짧은 부품)은 프라이머 tail 또는 올리고로만 들어간다.

수용 기준
  1. 전형 모듈이 fragment 명세를 받는다
  2. arm보다 짧은 부품은 **이웃 body에 병합되지 않는다** — tail 흡수 또는 독립 올리고
  3. 모든 PCR fragment의 `template`이 **입력 부품 중 하나와 정확히 일치**한다 (핵심)
  4. 주문 한계: 프라이머 ≤ PRIMER_ORDER_MAX_NT, 올리고 ≤ OLIGO_SYNTH_MAX_BP
  5. 공유 구간이 양쪽에 같은 서열로 존재하고, 길이가 범위 안이다
  6. 서열 오라클: module_seq == 부품 순서 이어붙임 (**padding 없음 = 흉터 0**)
  7. 길이 탐색이 **품질 통과하는 최소**를 고른다
  8. 재현성
실행: python tools_dev/params_cp_tests/g1_test.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from pydna.dseqrecord import Dseqrecord

from params import OLIGO_SYNTH_MAX_BP, PRIMER_ORDER_MAX_NT, DesignParams
from workflow.junction_plan import JunctionPlanner

FAIL = []


def check(cond, msg):
    print(("  ✅ " if cond else "  ❌ ") + msg)
    if not cond:
        FAIL.append(msg)


def part(pid, ptype, seq):
    r = Dseqrecord(seq, name=pid)
    r.annotations["part_id"] = pid
    r.annotations["part_type"] = ptype
    r.annotations["molecule_type"] = "DNA"
    return r


def seeded(n, seed, gc=0.5):
    rng = random.Random(seed)
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT")
                   for _ in range(n))


def typical(tag=""):
    """전형 모듈 — rbs(16bp)가 arm보다 짧아 tail/올리고 대상이 된다."""
    return [
        part("P" + tag, "promoter", seeded(200, 11)),
        part("R" + tag, "rbs", seeded(16, 12)),
        part("C" + tag, "cds", seeded(900, 13)),
        part("T" + tag, "terminator", seeded(140, 14)),
    ]


print("=" * 72)
print("G1 — 접합부 설계기")
print("=" * 72)

prm = DesignParams()
prm.validate()
plan = JunctionPlanner(params=prm).plan(typical())

print(f"\n[1] 전형 모듈 명세 (arm {plan.arm_min}~{plan.arm_max}bp)")
for f in plan.fragments:
    print(f"  {f.name:<16} {f.production:<6} template {len(f.template):>4}bp  "
          f"body {len(f.body):>4}bp  주문 {f.order_len:>4}nt  "
          f"L_arm {len(f.left_arm):>2} R_arm {len(f.right_arm):>2}"
          + (f"  [흡수 {len(f.prefix) + len(f.suffix)}bp]" if f.prefix or f.suffix else ""))
check(len(plan.fragments) >= 2, f"fragment 생성됨 ({len(plan.fragments)}개)")

print("\n[2] 짧은 부품이 이웃 body에 병합되지 않는다")
inputs = {str(p.annotations["part_id"]): str(p.seq).upper() for p in typical()}
merged_body = [f.name for f in plan.fragments
               if f.production == "pcr" and f.template not in inputs.values()]
check(not merged_body, f"PCR fragment의 template이 모두 입력 부품 그대로 (위반: {merged_body})")
absorbed = [f.name for f in plan.fragments if f.prefix or f.suffix]
oligo_sp = [f.name for f in plan.fragments if f.production == "oligo" and not f.template]
print(f"  tail 흡수 {absorbed} · 독립 올리고 {oligo_sp}")
check(bool(absorbed or oligo_sp), "짧은 부품이 tail 또는 올리고로 처리됨")

print("\n[3] ★ PCR 템플릿이 실제 존재하는 분자인가 (원래 버그)")
bad = []
for f in plan.fragments:
    if f.production != "pcr":
        continue
    if f.template not in inputs.values():
        bad.append((f.name, len(f.template)))
check(not bad, f"모든 PCR 템플릿이 입력 부품과 일치 (위반: {bad})")

print("\n[4] 주문 한계")
over_p, over_o = [], []
for f in plan.fragments:
    if f.production == "oligo":
        if f.order_len > OLIGO_SYNTH_MAX_BP:
            over_o.append((f.name, f.order_len))
    else:
        for side, tail in (("fwd", f.left_tail), ("rev", f.right_tail)):
            if len(tail) + prm.primer_len_min > PRIMER_ORDER_MAX_NT:
                over_p.append((f.name, side, len(tail)))
check(not over_p, f"프라이머 ≤ {PRIMER_ORDER_MAX_NT}nt (초과: {over_p})")
check(not over_o, f"올리고 ≤ {OLIGO_SYNTH_MAX_BP}nt (초과: {over_o})")

print("\n[5] 공유 구간")
mismatch = []
for i in range(len(plan.fragments) - 1):
    left, right = plan.fragments[i], plan.fragments[i + 1]
    j = plan.junctions[i]
    if not (left.seq.endswith(j.window) and right.seq.startswith(j.window)):
        mismatch.append(i)
    print(f"  접합부{i}: {j.left_name[:14]} ↔ {j.right_name[:14]}  "
          f"{j.length}bp  운반 {j.carrier}  GC {j.gc:.0f}%")
check(not mismatch, f"좌측 끝 == 우측 시작 == window (불일치: {mismatch})")
check(all(plan.arm_min <= j.length <= plan.arm_max for j in plan.junctions),
      "모든 길이가 범위 안")

print("\n[6] 서열 오라클 — padding 없음(흉터 0)")
raw = "".join(str(p.seq).upper() for p in typical())
check(plan.module_seq == raw,
      f"module_seq == 부품 순서 이어붙임 ({len(plan.module_seq)}bp == {len(raw)}bp)")

print("\n[7] 길이 탐색이 '품질 통과하는 최소'를 고른다")
# 창은 **운반자의 반대편** body에서 떼어낸다. 둘 다 PCR이면 body가 짧은 쪽이 운반하므로,
# LOWGC를 뒤쪽 이웃보다 길게 만들어야 LOWGC의 끝이 창이 된다.
#   끝 40bp = poly-A(GC 0%) → 30~49bp에서 전부 탈락
#   50bp = G 10개 + A 40개 → GC 20% → 통과
crafted = [
    part("GOOD", "cds", seeded(400, 31)),
    part("LOWGC", "cds", seeded(300, 32) + "G" * 12 + "A" * 40),   # 352bp
    part("TAIL", "cds", seeded(200, 33)),                          # 200bp (더 짧게)
]
pl_w = JunctionPlanner(params=DesignParams()).plan(crafted, arm_min=30, arm_max=60)
for j in pl_w.junctions:
    print(f"  접합부{j.index}: 선택 {j.length}bp (GC {j.gc:.0f}%), 떨어진 길이 {len(j.tried)}개")
    if j.tried:
        print(f"      첫 실패: {j.tried[0]}")
    check(len(j.tried) == j.length - 30,
          f"접합부{j.index}: 선택 길이 이전을 모두 시도했고 그중 통과가 없었음")
grew = [j for j in pl_w.junctions if j.length > 30]
check(bool(grew), f"최소에서 떨어지면 길이를 늘려 통과시킴 (늘어난 접합부 {len(grew)}개)")

print("\n[8] 재현성")
a = JunctionPlanner(params=DesignParams()).plan(typical())
b = JunctionPlanner(params=DesignParams()).plan(typical())
check([f.seq for f in a.fragments] == [f.seq for f in b.fragments], "같은 입력 → 같은 fragment")
check([j.length for j in a.junctions] == [j.length for j in b.junctions], "같은 입력 → 같은 길이")

print("\n[9] sp가 lp 사이에 낀 경우 — 양옆은 맨몸 PCR")
mid = [
    part("L1", "cds", seeded(600, 21)),
    part("S1", "rbs", seeded(16, 22)),
    part("L2", "cds", seeded(700, 23)),
]
pm = JunctionPlanner(params=DesignParams(yeast_homology_min=40, yeast_homology_max=60)).plan(
    mid, arm_min=40, arm_max=60)
for f in pm.fragments:
    print(f"  {f.name:<10} {f.production:<6} fwd_tail {len(f.left_tail):>3}nt "
          f"rev_tail {len(f.right_tail):>3}nt  주문 {f.order_len}nt")
naked = [f.name for f in pm.fragments
         if f.production == "pcr" and not f.left_tail and not f.right_tail]
check(bool(naked), f"tail 없는 맨몸 PCR이 존재 ({naked})")

print("\n" + "=" * 72)
if FAIL:
    print(f"G1 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("G1 전부 통과")
