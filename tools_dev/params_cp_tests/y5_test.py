"""Y5 검증 — homology 경로의 **주문 시트와 module feature**
(2026-10-25에 실제로 터진 버그).

## 무엇이 터졌나

조립 서열은 정확했는데 **그걸 실제로 만들 방법이 주문서에 없었다.**

1. `write_primer_csv`가 올리고 컬럼을 안 적었다. `_design_oligo_fragment`가 Top/Bottom
   가닥까지 계산해 `design`에 실어 보내는데 CSV에서 통째로 버려졌다.
2. 올리고 주문 시트는 `design`이 아니라 **module의 feature**에서 부품을 뽑았다.
3. 그런데 feature 기입이 `if method == "golden_gate"`로 막혀 있어 homology 경로의
   module은 feature가 **0개**였다. → 올리고 시트가 빈 파일.

결과: vector homology spacer 2개가 프라이머 시트에도, 올리고 시트에도 안 나왔다.

## 수용 기준

  1. homology module에 부품 feature가 기입되고 **좌표가 서열과 일치**한다
  2. part_type이 GenBank 표준 feature type으로 반영된다(promoter/CDS/terminator)
  3. 접합부 feature가 설계기(`JunctionPlan.junctions`)의 window와 **정확히 일치**한다
     — 운반자가 left든 right든(경계의 좌/우가 갈린다)
  4. `FragmentSpec.part_spans`가 흡수된 부품까지 순서대로 담는다
  5. `write_primer_csv`에 올리고 컬럼이 있고 올리고 fragment에 가닥이 채워진다
  6. Golden Gate 경로의 feature 기입은 예전 그대로다(회귀 방지)
  7. 길이 0 조각이 섞인 location을 정규화해, 작은 feature가 지도를 덮거나
     중복 판정을 빠져나가지 않는다

실행: python tools_dev/params_cp_tests/y5_test.py
"""
import csv
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from pydna.dseqrecord import Dseqrecord

from assembly.factory import AssemblyFactory
from params import DesignParams
from workflow.biopart_assemble import BiopartWorkflow
from workflow.builder import WorkflowBuilder

FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


def mk(name, ptype, seq):
    r = Dseqrecord(seq, name=name)
    r.annotations["part_type"] = ptype
    r.annotations["part_id"] = name
    return r


def synth(n, seed):
    """GC ~50%의 **서로 다른** 합성 부품.

    ⚠️ 같은 단위를 반복하면 부품끼리 서열을 공유하게 되고, yeast HR의 교차 차단
       (`_check_plan`)이 설계를 중단시킨다 — 그건 이 테스트가 보려는 게 아니다.
       seed마다 다른 난수열을 쓴다.
    """
    import random
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))

print("=" * 72)
print("Y5 — homology 경로의 주문 시트와 module feature")
print("=" * 72)

prm = DesignParams()
parts = [
    mk("P1", "promoter", synth(300, 11)),
    mk("C1", "cds", synth(600, 22)),
    mk("T1", "terminator", synth(250, 33)),
]

print("\n[1] part_spans — 배치 정보가 실린다")
a = AssemblyFactory.create_assembler("yeast_hr", params=prm)
design, plan = a.prepare_parts(parts)
spans_ok = all(len(f.part_lens) == len(f.part_ids) for f in plan.fragments)
check(spans_ok, "모든 fragment가 part_ids와 같은 길이의 part_lens를 갖는다")
for f in plan.fragments:
    sp = f.part_spans
    check(len(sp) == len(f.part_ids), f"{f.name}: span {len(sp)}개 == 부품 {len(f.part_ids)}개")
    check(all(f.body[s:e] for _i, _t, s, e in sp), f"{f.name}: span이 body 안에 든다")

print("\n[2] 올리고 fragment에 Top/Bottom이 실린다")
oli = [d for d in design if d.get("production") == "oligo"]
pcr = [d for d in design if d.get("production") != "oligo"]
print(f"  design {len(design)}건 — 올리고 {len(oli)} / PCR {len(pcr)}")
check(all(d.get("oligo_top") and d.get("oligo_bottom") for d in oli),
      "올리고 fragment는 두 가닥을 모두 갖는다")
check(all(d.get("forward_primer") and d.get("reverse_primer") for d in pcr),
      "PCR fragment는 프라이머 두 개를 모두 갖는다")

print("\n[3] primers.csv에 올리고 컬럼이 있다")
with tempfile.TemporaryDirectory() as td:
    wf = BiopartWorkflow(project_dir=td, parts=[])
    out = wf.write_primer_csv(design, "primers.csv")
    rows = list(csv.DictReader(open(out, encoding="utf-8")))
    cols = set(rows[0].keys())
    for c in ("production", "oligo_top", "oligo_bottom", "oligo_len", "duplex_grade"):
        check(c in cols, f"컬럼 '{c}' 존재")
    csv_oli = [r for r in rows if r["production"] == "oligo"]
    check(len(csv_oli) == len(oli), f"CSV 올리고 행 {len(csv_oli)} == design {len(oli)}")
    check(all(r["oligo_top"] and r["oligo_bottom"] for r in csv_oli),
          "CSV 올리고 행에 가닥이 비어 있지 않다")

print("\n[4] homology module에 feature가 기입된다 (좌표·type 검증)")
for method in ("yeast_hr", "gibson"):
    a = AssemblyFactory.create_assembler(method, params=prm)
    _d, plan = a.prepare_parts(parts)
    recs = plan.records()
    with tempfile.TemporaryDirectory() as td:
        wf = BiopartWorkflow(project_dir=td, parts=[])
        mod = wf.module_assembly(recs, method=method, assembly_type="linear", design_params=prm)
    r = mod["record"]
    s = str(r.seq).upper()
    labels = {(f.qualifiers.get("label") or [""])[0]: f for f in r.features}
    print(f"  {method}: module {len(s)}bp feature {len(r.features)}개")
    check(len(r.features) > 0, f"{method}: feature가 0개가 아니다")
    for p in parts:
        f = labels.get(p.name)
        check(f is not None, f"{method}: {p.name} feature 존재")
        if f is None:
            continue
        got = s[int(f.location.start):int(f.location.end)]
        check(got == str(p.seq).upper(), f"{method}: {p.name} 좌표가 서열과 일치")
    check(labels.get("P1") is not None and labels["P1"].type == "promoter",
          f"{method}: part_type이 GenBank type으로 (promoter)",
          labels.get("P1").type if labels.get("P1") else "")
    check(labels.get("C1") is not None and labels["C1"].type == "CDS",
          f"{method}: part_type이 GenBank type으로 (CDS)",
          labels.get("C1").type if labels.get("C1") else "")

    # 접합부 feature ↔ 설계기 window
    hom = [f for f in r.features if (f.qualifiers.get("label") or [""])[0] == "homology"]
    check(len(hom) == len(plan.junctions),
          f"{method}: 접합부 feature {len(hom)} == 설계기 {len(plan.junctions)}")
    for i, f in enumerate(hom):
        got = s[int(f.location.start):int(f.location.end)]
        check(got == plan.junctions[i].window.upper(),
              f"{method}: 접합부 #{i} 서열 일치 (운반={plan.junctions[i].carrier})")

print("\n[5] Golden Gate 경로는 예전 그대로 (회귀 방지)")
gg_recs, _gd = WorkflowBuilder._prep_module_records(parts, params=prm, module_method="golden_gate")
with tempfile.TemporaryDirectory() as td:
    wf = BiopartWorkflow(project_dir=td, parts=[])
    gg = wf.module_assembly(gg_recs, method="golden_gate", assembly_type="linear", design_params=prm)
check(gg is not None and len(gg["record"].features) > 0,
      "golden_gate module에도 여전히 feature가 기입된다",
      str(len(gg["record"].features)) if gg else "조립 실패")

print("\n[6] 길이 0 조각이 섞인 location 정규화")
# pydna가 환형 조립에서 join{[n:n],[0:k]} 같은 퇴화된 조각을 내보낸다.
# 그대로 두면 GenBank에 join(n^1,1..k)로 적히고, 숫자만 긁어 범위를 잡는 쪽(클라이언트
# 벡터맵)이 1..n으로 읽어 **작은 feature가 지도를 통째로 덮는다**(2026-10-25 실측).
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from assembly.common import normalize_location, remove_duplicate_features

degenerate = CompoundLocation([FeatureLocation(5381, 5381, strand=1),
                               FeatureLocation(0, 30, strand=1)])
norm = normalize_location(degenerate)
check(int(norm.start) == 0 and int(norm.end) == 30,
      "길이 0 조각을 떼고 [0,30)으로 좁힌다", f"{int(norm.start)}-{int(norm.end)}")
check(not isinstance(norm, CompoundLocation), "조각이 하나면 CompoundLocation을 푼다")

real_wrap = CompoundLocation([FeatureLocation(5352, 5381, strand=1),
                              FeatureLocation(0, 30, strand=1)])
check(normalize_location(real_wrap) is real_wrap,
      "진짜 원점 wrap(두 조각 다 길이>0)은 건드리지 않는다")
check(normalize_location(FeatureLocation(10, 20, strand=1)).start == 10,
      "단일 구간은 그대로")

# 정규화 덕에 같은 구간의 중복이 실제로 합쳐지는가
rec = Dseqrecord("A" * 5381, circular=True)
quals = {"label": ["spacer"]}
rec.features = [SeqFeature(degenerate, type="misc_feature", qualifiers=dict(quals)),
                SeqFeature(FeatureLocation(0, 30, strand=1), type="misc_feature",
                           qualifiers=dict(quals))]
remove_duplicate_features(rec)
check(len(rec.features) == 1,
      "퇴화 join과 같은 구간의 단일 location이 중복으로 합쳐진다", f"{len(rec.features)}개")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y5 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y5 전부 통과")
