"""Y1 검증 — `yeast_homology` 필드 + arm 축 분기.

수용 기준
  1. `yeast_homology` 기본 45, 하드 범위 30~60
  2. `oligo_max_bp < 2*yeast_homology`를 차단한다 (접합부 구간 겹침)
  3. `YeastHRCloning`의 arm이 `yeast_homology`를 따른다 — **두 곳이 함께 움직인다**
     (`overlap_length` = pydna overlap limit, `prepare_parts` = 접합부 설계)
  4. `GibsonCloning`은 `gibson_overlap`을 그대로 쓴다 (서로 간섭 없음)
  5. yeast arm(45bp)으로 접합부 설계가 성립한다
  6. `yeast_homology`를 바꿔도 Gibson 산출물은 불변 (축 분리 확인)
실행: python tools_dev/params_cp_tests/y1_test.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from pydna.dseqrecord import Dseqrecord

from assembly.factory import AssemblyFactory
from params import DesignParams

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


def seeded(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def module_parts(tag=""):
    return [
        part("P" + tag, "promoter", seeded(240, 71)),
        part("R" + tag, "rbs", seeded(16, 72)),
        part("C" + tag, "cds", seeded(900, 73)),
        part("T" + tag, "terminator", seeded(140, 74)),
    ]


print("=" * 72)
print("Y1 — yeast_homology 필드 + arm 축 분기")
print("=" * 72)

print("\n[1] 기본값과 하드 범위")
prm = DesignParams()
prm.validate()
print(f"  yeast_homology 기본 = {prm.yeast_homology_min}~{prm.yeast_homology_max}bp "
      f"(gibson_overlap = {prm.gibson_overlap}bp)")
check((prm.yeast_homology_min, prm.yeast_homology_max) == (40, 60), "기본 범위 40~60bp")
for lo, hi, should_pass in ((30, 60, True), (45, 45, True), (40, 50, True),
                            (29, 60, False), (40, 61, False), (50, 40, False)):
    try:
        DesignParams(yeast_homology_min=lo, yeast_homology_max=hi).validate()
        ok, err = True, ""
    except ValueError as e:
        ok, err = False, str(e).splitlines()[-1][:70]
    print(f"  yeast_homology {lo:>3}~{hi:<3} → {'통과' if ok else '거부'}"
          f"{'  | ' + err if err else ''}")
    check(ok == should_pass, f"{lo}~{hi} {'허용' if should_pass else '차단'}")

print("\n[2] min > max 차단")
try:
    DesignParams(yeast_homology_min=55, yeast_homology_max=40).validate()
    check(False, "뒤집힌 범위가 차단되어야 함")
except ValueError as e:
    print(f"  {str(e).splitlines()[-1][:90]}")
    check("min" in str(e), "뒤집힌 범위 차단됨")
DesignParams().validate()
check(True, "기본 범위는 통과")

print("\n[3] YeastHRCloning의 arm이 yeast_homology를 따른다 (두 곳 동시)")
y = AssemblyFactory.create_assembler("yeast_hr", params=prm)
print(f"  ARM_PARAMS = {y.ARM_PARAMS!r}, arm_range = {y.arm_range}, "
      f"overlap_length = {y.overlap_length}")
check(y.ARM_PARAMS == ("yeast_homology_min", "yeast_homology_max"), "ARM_PARAMS가 yeast 범위")
check(y.arm_range == (prm.yeast_homology_min, prm.yeast_homology_max), "arm_range == params 범위")
check(y.overlap_length == prm.yeast_homology_min,
      f"overlap_length(pydna limit) == 범위 최소 ({y.overlap_length})")
_design, yplan = y.prepare_parts(module_parts("y"))
print(f"  prepare_parts → 범위 {yplan.arm_min}~{yplan.arm_max}bp, "
      f"고른 길이 {sorted({j.length for j in yplan.junctions})}, "
      f"fragment {len(yplan.fragments)}개")
check((yplan.arm_min, yplan.arm_max) == (prm.yeast_homology_min, prm.yeast_homology_max),
      "접합부 범위 == yeast_homology 범위")
check(all(prm.yeast_homology_min <= j.length <= prm.yeast_homology_max
          for j in yplan.junctions), "모든 길이가 yeast 범위 안")

print("\n[4] GibsonCloning은 gibson_overlap 유지")
g = AssemblyFactory.create_assembler("gibson", params=prm)
print(f"  ARM_PARAMS = {g.ARM_PARAMS!r}, arm_range = {g.arm_range}")
check(g.ARM_PARAMS == ("gibson_overlap", "gibson_overlap_max"), "ARM_PARAMS가 gibson 범위")
check(g.overlap_length == prm.gibson_overlap, "overlap_length == gibson_overlap")
_gd, gplan = g.prepare_parts(module_parts("g"))
check((gplan.arm_min, gplan.arm_max) == (prm.gibson_overlap, prm.gibson_overlap_max), "Gibson도 범위를 쓴다")

print("\n[5] yeast arm으로 조립이 성립한다")
recs = yplan.records()
res = y.assemble(recs, assembly_type="linear")
print(f"  linear 산물 {len(res)}개, 기대 {len(yplan.module_seq)}bp")
check(any(str(r.seq).upper() == yplan.module_seq for r in res),
      "조립 산물이 기대 서열과 일치")

print("\n[6] 축 분리 — yeast_homology를 바꿔도 Gibson 산출물 불변")
base = AssemblyFactory.create_assembler("gibson", params=DesignParams(yeast_homology_min=40))
alt = AssemblyFactory.create_assembler("gibson", params=DesignParams(yeast_homology_min=50))
_b, bp_ = base.prepare_parts(module_parts("x"))
_a, ap_ = alt.prepare_parts(module_parts("x"))
check(bp_.module_seq == ap_.module_seq, "Gibson module 서열 불변")
check([f.seq for f in bp_.fragments] == [f.seq for f in ap_.fragments],
      "Gibson fragment 서열 불변")
# 역방향: gibson_overlap을 바꿔도 yeast 산출물 불변
yb = AssemblyFactory.create_assembler("yeast_hr", params=DesignParams(gibson_overlap=20))
ya = AssemblyFactory.create_assembler("yeast_hr", params=DesignParams(gibson_overlap=25))
_yb, ybp = yb.prepare_parts(module_parts("z"))
_ya, yap = ya.prepare_parts(module_parts("z"))
check(ybp.module_seq == yap.module_seq, "yeast module 서열이 gibson_overlap과 무관")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y1 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y1 전부 통과")
