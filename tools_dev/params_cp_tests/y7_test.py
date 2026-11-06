"""Y7 검증 — GC 판정창을 설계창에서 분리하고 2층으로 나눴다 (2026-11-05).

## 무엇이 터졌나

Y6(2026-11-04)에서 **Tm** 판정창을 설계창에서 분리했는데, **똑같은 모순이 GC 축에
그대로 남아 있었다.** 설계창은 [20, 80]%(primer3 입력)인데 판정창은 [35, 60]%라,
primer3가 **설계 목표대로** 집어온 프라이머가 판정에서 떨어졌다.

게다가 GC 지적은 `anneal_gate`를 타서 tail 유무로 등급이 갈렸다 — mCherry rev의
실제 anneal 구간(`TTATTTATACAGTTCATCCATACCA`, GC 28.0%)으로 실측:

| | 끝단(tail 없음) | tail 있음 |
|---|---|---|
| before (차단 35~60) | **🔴 red** | 🟡 yellow |
| after (차단 20~80 + 경고 40~60) | 🟡 yellow | 🟡 yellow |

**같은 서열, 같은 GC인데 arm 운반자가 누구냐에 따라 적/황이 갈렸다.**

## 왜 Tm처럼 그냥 넓히지 않았나

GC는 Tm과 사정이 하나 다르다 — **설계자가 고칠 수 없는 값**이다. 끝단 프라이머는
`terminal_tm_min`까지 길이를 늘려 Tm은 맞출 수 있지만, 늘려도 그 구간의 GC%는 거의
그대로다(템플릿의 성질). 그래서 차단창만 설계창까지 넓히면 GC 신호가 통째로 사라진다.

→ **층을 둘로 나눴다**(`duplex_safe_bp`/`duplex_risky_bp`와 같은 구조):
  · `qc_gc_min/max`       — 차단창(red). 설계창을 포함한다.
  · `qc_gc_ideal_min/max` — 경고창(yellow). 차단 안 하되 초록으로 숨기지도 않는다.

35%가 틀린 컷오프라는 건 추측이 아니라 **이 프로젝트의 자체 실측**이다 —
`workflow/junction_plan` docstring "품질 기준" 절: 접합부 축에서 GC 35~60%는 실제
파츠 말단의 32%를 떨어뜨렸고(통과 68%), 20~80%는 99.3%가 통과했다. 접합부 축은
그때 옮겼는데 프라이머 판정축만 35%에 남아 있었다.

## 수용 기준

  1. 기본값에서 차단창이 설계창을 포함한다 + 소프트 경고가 하나도 안 뜬다
  2. GC 28% anneal이 tail 유무와 **무관하게** 같은 등급이다(적/황 갈림 제거)
  3. 설계창 안(20~80%)의 GC는 차단되지 않는다
  4. 설계창 밖(>80%, <20%)은 여전히 차단된다
  5. 이상 구간(40~60%) 안쪽은 GC 지적이 없다(경고도 소음이 되지 않는다)
  6. 경고창이 차단창 밖으로 나가면 `validate()`가 **하드 에러**
     (경고 없이 차단되는 구간 방지 — validate()는 워크플로우 진입점마다 호출된다)
  7. 차단창이 설계창을 못 덮으면 `validate()`가 경고한다
  8. 호출자가 min_gc/max_gc를 직접 주면 경고창이 그 안쪽으로 클램프된다
  9. GC %/분수 단위 실수가 경고창 필드에서도 막힌다
 10. 실제 케이스(mCherry + pRS425, yeast HR) E2E에서 GC 차단이 없다

실행: python tools_dev/params_cp_tests/y7_test.py
"""
import io
import contextlib
import shutil
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from pydna.dseqrecord import Dseqrecord

from basic.primer_design import ValidatePrimer
from params import DesignParams
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
SANDBOX = (ROOT.parents[2] / "galaxy_dev" / "BioCAD_galaxy" / "tools" / "awsemble"
           / "user_parts_sandbox_2.db")
FAIL = []

# mCherry rev의 실제 anneal 구간 — GC 28.0%, 설계창 [20,80] 안쪽이다.
AT_ANNEAL = "TTATTTATACAGTTCATCCATACCA"
# pRS425_homR (40bp) — tail 유무 비교용. 이 자체가 GC 판정에 들어가면 안 된다.
TAIL40 = "AGCTCCACCGCGGTGGCGGCCGCTCTAGAACTAGTGGATC"


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


def gc_tier(primer, tail_len=0, params=None):
    """GC **축만** 보고 등급을 돌려준다(다른 축의 지적은 세지 않는다)."""
    r = ValidatePrimer().validate_primer(primer, tail_len=tail_len, params=params)
    v = [m for m in r["violations"] if "GC%" in m]
    w = [m for m in r["warnings"] if "GC%" in m]
    return ("red" if v else ("yellow" if w else "green")), r["gc_content"], v + w


def synth(n, gc_pct):
    """목표 GC%에 맞춘 결정적 서열 — 반복을 피해 다른 축 지적을 줄인다."""
    n_gc = round(n * gc_pct / 100.0)
    pat_gc, pat_at = "GCCG", "ATTA"
    out = []
    for i in range(n):
        src = pat_gc if i < n_gc else pat_at
        out.append(src[i % 4])
    return "".join(out)


prm = DesignParams()

print(f"[1] 기본값: 차단창이 설계창을 포함 + 소프트 경고 없음")
check(prm.qc_gc_min <= prm.gc_min and prm.qc_gc_max >= prm.gc_max,
      "차단창 ⊇ 설계창",
      f"차단 [{prm.qc_gc_min}, {prm.qc_gc_max}] ⊇ 설계 [{prm.gc_min}, {prm.gc_max}]")
check(prm.qc_gc_min <= prm.qc_gc_ideal_min <= prm.qc_gc_ideal_max <= prm.qc_gc_max,
      "경고창 ⊆ 차단창",
      f"경고 [{prm.qc_gc_ideal_min}, {prm.qc_gc_ideal_max}]")
warns = prm.validate()
check(not warns, "기본값에서 소프트 경고가 하나도 없다 (소음 제거)", str(warns))

print(f"\n[2] GC 28% anneal — tail 유무로 등급이 갈리지 않는다")
t_bare, gc_bare, m_bare = gc_tier(AT_ANNEAL, 0)
t_tail, gc_tail, m_tail = gc_tier(TAIL40 + AT_ANNEAL, len(TAIL40))
print(f"  끝단(tail 없음): {t_bare:<7} GC {gc_bare:.1%}")
print(f"  tail 있음      : {t_tail:<7} GC {gc_tail:.1%}")
check(t_bare == t_tail, "운반자와 무관하게 같은 등급", f"{t_bare} == {t_tail}")
check(t_bare == "yellow", "둘 다 yellow (차단 아님, 숨김도 아님)", t_bare)
check(abs(gc_bare - gc_tail) < 1e-9, "tail이 GC 계산에 섞이지 않았다",
      f"{gc_bare:.3f} / {gc_tail:.3f}")
# before 재현: 차단창만 되돌리면 끝단이 red가 된다 = 이 테스트가 지키는 회귀다
old = DesignParams(qc_gc_min=35.0, qc_gc_max=60.0)
t_old, _, _ = gc_tier(AT_ANNEAL, 0, params=old)
check(t_old == "red", "차단창을 35%로 되돌리면 끝단이 red로 돌아간다(회귀 감시)", t_old)

print(f"\n[3]~[5] GC 밴드별 등급")
BANDS = [(15.0, "red", "설계창 밖(너무 낮음)"),
         (28.0, "yellow", "설계창 안 / 이상 구간 밖"),
         (50.0, "green", "이상 구간 안"),
         (70.0, "yellow", "설계창 안 / 이상 구간 밖"),
         (90.0, "red", "설계창 밖(너무 높음)")]
for pct, want, why in BANDS:
    seq = synth(40, pct)
    got, gc, msgs = gc_tier(seq, 0)
    check(got == want, f"GC {gc:.0%} → {want} ({why})",
          f"got={got}" + (f" / {msgs[0][:48]}" if msgs else ""))

print(f"\n[6] 경고창이 차단창 밖이면 하드 에러")
# 하드 검증은 `validate()`가 한다 — 워크플로우 진입점마다 호출된다
# (`builder` build_plasmid/build_multimodule_plasmid/insert, `base.run_plan`,
#  galaxy `vector_insertion_helper`). 생성자는 검증하지 않는다.
for kw in ({"qc_gc_ideal_min": 10.0}, {"qc_gc_ideal_max": 95.0}):
    try:
        DesignParams(**kw).validate()
        check(False, f"{kw} → 에러가 나야 함")
    except ValueError as e:
        check("경고창" in str(e) and "차단창" in str(e), f"{kw} 차단",
              str(e).splitlines()[-1][:70])

print(f"\n[7] 차단창이 설계창을 못 덮으면 경고")
w = DesignParams(qc_gc_min=35.0, qc_gc_max=60.0).validate()
hit = next((x for x in w if "GC 차단창" in x), None)
check(hit is not None, "포함 관계 위반을 경고한다", (hit or "(경고 없음)")[:70])
check(not any("GC 차단창" in x for x in DesignParams().validate()),
      "기본값(20~80 ⊇ 20~80)은 경고 없음")

print(f"\n[8] 호출자가 차단창을 직접 주면 경고창이 클램프된다")
# min_gc=0.45로 좁히면 경고창 하한 40%는 45%로 당겨져야 한다
# → GC 50%는 green, GC 47%는 차단창 안이지만 경고창 밖이라 yellow가 아니라 green이어야
#   하는 게 아니라, 45% 하한이 적용돼 47%는 green, 43%는 red(차단창 밖)다.
t43, gc43, _ = gc_tier(synth(40, 43.0), 0, params=prm)
r_clamped = ValidatePrimer().validate_primer(synth(40, 43.0), min_gc=0.45, params=prm)
check(any("너무 낮습니다" in m for m in r_clamped["violations"]),
      "호출자 하한(45%) 밖은 차단된다", f"GC {r_clamped['gc_content']:.1%}")
check(not any("이상 구간" in m for m in r_clamped["warnings"]),
      "차단과 경고가 겹쳐 두 번 나오지 않는다")

print(f"\n[9] GC 단위 실수(0~1 분수) 방어가 경고창에도 걸린다")
for kw in ({"qc_gc_ideal_min": 0.4}, {"qc_gc_ideal_max": 0.6}):
    try:
        DesignParams(**kw).validate()
        check(False, f"{kw} → 에러가 나야 함")
    except ValueError as e:
        check("단위가" in str(e), f"{kw} 차단", str(e).splitlines()[-1][:70])

print(f"\n[10] 실제 케이스 E2E — mCherry + pRS425 (yeast HR)")
if not SANDBOX.exists():
    check(False, f"샌드박스 DB 없음: {SANDBOX}")
else:
    con = sqlite3.connect(str(SANDBOX))

    def grab(table, pid, circular=False):
        seq = con.execute(f"select sequence from {table} where ID=?", (pid,)).fetchone()[0]
        r = Dseqrecord(seq.strip().upper(), circular=circular, name=pid)
        r.annotations.update({"part_id": pid, "part_type": table, "molecule_type": "DNA"})
        return r

    cds, vec = grab("cds", "mCherry"), grab("vector", "pRS425", circular=True)
    d = S / "y7_mcherry"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        res = WorkflowBuilder.build_plasmid(
            [cds], vec, DatabaseImport(), str(d),
            module_method="gibson", insert_method="yeast_hr", params=DesignParams())
    construct = res.get("construct")
    check(construct is not None, "construct 산출")
    import csv as _csv
    rows = list(_csv.DictReader(open(d / "reports" / "primers.csv")))
    check(len(rows) == 1, "fragment 1종", f"{len(rows)}종")
    row = rows[0]
    viol = row.get("violations") or ""
    check("GC%" not in viol, "GC로 차단된 프라이머 없음", f"violations={viol or '없음'}")
    check(float(row["rev_gc"]) < prm.qc_gc_ideal_min,
          "rev는 여전히 이상 구간 밖이다(신호가 사라지지 않았다)",
          f"rev_gc={row['rev_gc']}% < {prm.qc_gc_ideal_min}%")
    check("이상 구간보다 낮습니다" in (row.get("warnings") or ""),
          "그 사실이 경고로 남아 있다")
    check(row.get("tier") == "yellow", "tier=yellow", f"tier={row.get('tier')}")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y7 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y7 전부 통과")
