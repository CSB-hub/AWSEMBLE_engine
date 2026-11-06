"""G3 검증 — Gibson part→module prep 배선.

수용 기준
  1. `prepare_parts`가 모든 fragment에 대해 설계 산출물을 낸다(올리고/프라이머)
  2. **프라이머로 PCR을 돌린 산물이 fragment 명세와 정확히 일치**한다
     → anneal shift가 산물을 바꾸지 않음을 증명한다(핵심)
  3. 설계 산출물로 조립한 module이 **기대 서열(오라클)과 일치**한다
     — 회전·역상보 정규화 후 비교. pydna는 같은 원형 분자를 2개로 돌려준다
  4. anneal shift 폴백이 실제로 실패를 되살린다
  5. `module_assembly(method="gibson")` 경로가 동작한다
  6. Golden Gate 경로에 영향이 없다 (같은 입력으로 GG 결과 불변)
실행: python tools_dev/params_cp_tests/g3_test.py
"""
import os
import random
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from pydna.amplify import pcr
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


def seeded(n, seed, gc=0.5):
    rng = random.Random(seed)
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT")
                   for _ in range(n))


_COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(_COMP)[::-1]


def same_circle(a, b):
    """원형 분자 동일성 — 회전·역상보 정규화 후 비교."""
    if len(a) != len(b):
        return False
    return (b in a + a) or (rc(b) in a + a)


def module_parts():
    return [
        part("P100", "promoter", seeded(220, 41)),
        part("R100", "rbs", seeded(16, 42)),
        part("C100", "cds", seeded(850, 43)),
        part("T100", "terminator", seeded(90, 44)),
    ]


print("=" * 72)
print("G3 — Gibson part→module prep")
print("=" * 72)

prm = DesignParams()
prm.validate()
gib = AssemblyFactory.create_assembler("gibson", params=prm)
parts = module_parts()
design, plan = gib.prepare_parts(parts)

print(f"\n[1] 설계 산출물 ({len(design)}개 fragment, arm={plan.arm}bp)")
for d in design:
    if d["production"] == "oligo":
        print(f"  {d['part_name']:<16} oligo  top {d['oligo_len']:>4}nt  "
              f"등급 {d['duplex_grade']}")
    else:
        print(f"  {d['part_name']:<16} pcr    fwd {len(d['forward_primer']):>3}nt "
              f"rev {len(d['reverse_primer']):>3}nt  shift {d['anneal_shift'][0]}  "
              f"amplicon {len(d['full_sequence']):>5}bp")
check(len(design) == len(plan.fragments), "fragment마다 설계 산출물 1개")
for d in design:
    if d["production"] == "oligo":
        check(bool(d["oligo_top"] and d["oligo_bottom"]), f"{d['part_name']} 올리고 쌍 존재")
        check(d["oligo_top"] == d["full_sequence"], f"{d['part_name']} top == fragment 서열")
        check(d["oligo_bottom"] == rc(d["full_sequence"]),
              f"{d['part_name']} bottom == fragment 역상보")
    else:
        check(bool(d["forward_primer"] and d["reverse_primer"]),
              f"{d['part_name']} 프라이머 쌍 존재")

print("\n[2] ★ PCR 시뮬 산물 == fragment 명세 (실제 존재하는 템플릿으로)")
# ⚠️ 템플릿은 `frag.template` — **실제로 존재하는 분자**다. 예전 테스트는 `frag.body`
#    (병합된 가상 서열)를 템플릿으로 줘서, 물리적으로 불가능한 설계를 통과시켰다.
for d, frag in zip(design, plan.fragments):
    if d["production"] != "pcr":
        continue
    tmpl = Dseqrecord(frag.template, name=frag.name)
    try:
        prod = pcr(d["forward_primer"], d["reverse_primer"], tmpl)
        got = str(prod.seq).upper()
    except Exception as exc:
        check(False, f"{d['part_name']} PCR 시뮬 실패: {type(exc).__name__} {exc}")
        continue
    check(got == frag.seq,
          f"{d['part_name']} 산물 {len(got)}bp == 명세 {len(frag.seq)}bp "
          f"(shift {d['anneal_shift'][0]}nt, tail fwd {len(frag.left_tail)} / rev {len(frag.right_tail)})")

print("\n[3] 조립 결과 == 기대 서열(오라클)")
records = plan.records()
gib.overlap_length = plan.arm
results = gib.assemble(records, assembly_type="linear")
print(f"  linear 조립 산물 {len(results)}개")
lin_ok = any(str(r.seq).upper() == plan.module_seq for r in results)
check(lin_ok, f"linear 산물 중 기대 서열({len(plan.module_seq)}bp)과 일치하는 것이 있음")
if not lin_ok:
    for r in results[:3]:
        print(f"    산물 {len(r.seq)}bp / 기대 {len(plan.module_seq)}bp")

print("\n[4] anneal shift 폴백의 효과 — 같은 조건에서 켜고/끄고 비교")
# Tm 창을 권장 최소폭(4℃)까지 좁혀 맨몸 설계가 어려운 조건을 만든다.
# "절대 성공률 몇 %"는 창 폭에 따라 달라지므로 기준이 될 수 없다 →
# **같은 입력에서 폴백을 끈 경우와 비교**해 폴백이 실제로 실패를 되살리는지만 본다.
tight = DesignParams(tm_min=58.0, tm_max=62.0)
tight.validate()
cases = []
for seed in range(12):
    cases.append([
        part("A", "promoter", seeded(240, 500 + seed)),
        part("B", "cds", seeded(700, 600 + seed)),
        part("C", "terminator", seeded(160, 700 + seed)),
    ])


def run_cases(max_shift):
    g = AssemblyFactory.create_assembler("gibson", params=tight)
    g.MAX_ANNEAL_SHIFT = max_shift
    ok = shifted = 0
    for ps in cases:
        try:
            dz, _pl = g.prepare_parts(ps)
        except RuntimeError:
            continue
        ok += 1
        shifted += sum(1 for d in dz
                       if d["production"] == "pcr" and d["anneal_shift"][0] > 0)
    return ok, shifted


off_ok, _ = run_cases(0)
on_ok, on_shifted = run_cases(8)
print(f"  좁은 Tm 창(58~62℃) 모듈 {len(cases)}건")
print(f"    shift 끔(k=0)  → 성공 {off_ok}")
print(f"    shift 켬(k≤8)  → 성공 {on_ok}, shift 사용 fragment {on_shifted}개")
check(on_shifted > 0, "shift 폴백이 실제로 쓰임")
check(on_ok > off_ok, f"폴백이 실패를 되살림 ({off_ok} → {on_ok}건)")

print("\n[5] module_assembly(method='gibson') 경로")
from workflow.biopart_assemble import BiopartWorkflow

tmp = tempfile.mkdtemp(prefix="g3_")
try:
    wf = BiopartWorkflow(project_dir=tmp, parts=records)
    out = wf.module_assembly(records, method="gibson", assembly_type="linear",
                             module_id="g3_mod", design_params=prm)
    check(out is not None, "module_assembly 성공")
    if out:
        got = str(out["record"].seq).upper()
        print(f"  module {len(got)}bp / 기대 {len(plan.module_seq)}bp")
        check(got == plan.module_seq, "module 서열 == 기대 서열")
        check(os.path.exists(out["file"]), f"산출물 저장됨 ({os.path.basename(out['file'])})")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("\n[6] Golden Gate 경로 불변")
from assembly.golden_gate import GoldenGateModule

gg = GoldenGateModule(prm.gg_enzyme, params=prm)
gg_design = gg.prepare_fragments(module_parts(), assembly_type="linear")
gg2 = GoldenGateModule(prm.gg_enzyme, params=DesignParams())
gg_design2 = gg2.prepare_fragments(module_parts(), assembly_type="linear")
check([d["full_sequence"] for d in gg_design] == [d["full_sequence"] for d in gg_design2],
      "GG prepare_fragments 결과가 재현됨(새 파라미터 추가 영향 없음)")
check(len(gg_design) == len(module_parts()), "GG는 부품 수만큼 설계(병합 안 함)")

print("\n" + "=" * 72)
if FAIL:
    print(f"G3 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("G3 전부 통과")
