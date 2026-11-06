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

⚠️ **2026-11-05(Y11)에 처방이 뒤집혔다.** 위 진단(35%가 틀린 컷오프였다)은 그대로
유효하지만, 처방인 "층을 둘로 나눈다"는 **폐기**했다 — 같은 2층 구조가 그 뒤 길이
축에서 세 번째 같은 버그를 냈고(판정값이 `15`/`30` 리터럴이라 유저가 Q5 상한 26nt를
골라도 경고는 30 기준으로 떴다), 무엇보다 **유저가 설정한 적 없는 기준으로 경고가
떴다**(GC 이상구간 경고가 78행 중 39행, 50%). GC는 템플릿의 성질이라 그 경고를 봐도
실험자가 할 수 있는 일이 없었다.

→ Y11의 처방: **판정 기준은 유저 선언창 하나**(`gc_min/gc_max` 기본 20~80%).
  경고는 "내가 설정한 범위를 벗어났을 때"만 뜬다. `qc_gc_*`·`qc_gc_ideal_*` 제거.

  1. 판정 기준이 선언창과 **같은 숫자**다 (두 창이 어긋날 여지가 없다)
  2. 기본값에서 소프트 경고가 하나도 안 뜬다
  3. GC 28% anneal(원래 Y7이 yellow로 만든 값)이 이제 **조용하다** — 선언창 안이므로
  4. 선언창 밖(>80%, <20%)은 여전히 지적된다
  5. tail이 GC 계산에 섞이지 않는다 (운반자가 값을 바꾸지 않는다)
  6. 선언창을 좁히면 지적이 **따라 움직인다** (기준이 유저 입력에 묶여 있다)
  7. 호출자가 min_gc/max_gc를 직접 주면 그것이 우선한다
  8. GC %/분수 단위 실수가 선언창 필드에서 막힌다
  9. 실제 케이스(mCherry + pRS425, yeast HR) E2E에서 GC 차단이 없다

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
# 픽스처는 이 저장소 밖(Galaxy 트리)에 있다. 재구성 전에는 형제 디렉토리였지만
# 지금은 깊이가 달라 맞지 않는다 — 알려진 자리를 순서대로 찾는다.
_SB_NAME = "user_parts_sandbox_2.db"
SANDBOX = next(
    (p for p in (
        ROOT.parents[2] / "galaxy_dev" / "BioCAD_galaxy" / "tools" / "awsemble" / _SB_NAME,
        Path("/home/csb/Desktop/BioCAD/Awsemble/galaxy_dev/BioCAD_galaxy/tools/awsemble") / _SB_NAME,
    ) if p.exists()),
    ROOT.parents[2] / "galaxy_dev" / "BioCAD_galaxy" / "tools" / "awsemble" / _SB_NAME)
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

print(f"[1] 판정 기준이 선언창과 같은 숫자다 + 기본값에 소프트 경고 없음")
check(not hasattr(prm, "qc_gc_min") and not hasattr(prm, "qc_gc_ideal_min"),
      "GC 2층(차단창·경고창)이 제거됐다")
check((prm.gc_min, prm.gc_max) == (20.0, 80.0),
      "선언창 하나만 남았다", f"[{prm.gc_min:g}, {prm.gc_max:g}]%")
# 판정이 실제로 그 숫자를 쓰는지 — 경계 바로 안/밖으로 확인한다(필드만 보면 못 잡는다).
_in = gc_tier(synth(40, 22.0), 0)
_out = gc_tier(synth(40, 15.0), 0)
check(_in[0] == "green" and _out[0] == "red",
      "선언창 경계가 실제 판정 경계다", f"22%→{_in[0]} / 15%→{_out[0]}")
warns = prm.validate()
check(not warns, "기본값에서 소프트 경고가 하나도 없다 (소음 제거)", str(warns))

print(f"\n[2] GC 28% anneal — 이제 조용하다 (Y7에선 yellow였다)")
t_bare, gc_bare, m_bare = gc_tier(AT_ANNEAL, 0)
t_tail, gc_tail, m_tail = gc_tier(TAIL40 + AT_ANNEAL, len(TAIL40))
print(f"  끝단(tail 없음): {t_bare:<7} GC {gc_bare:.1%}")
print(f"  tail 있음      : {t_tail:<7} GC {gc_tail:.1%}")
check(t_bare == t_tail, "운반자와 무관하게 같은 등급", f"{t_bare} == {t_tail}")
check(not any("GC" in x for x in m_bare + m_tail),
      "GC 지적이 없다 — 선언창 20~80% 안이므로",
      str([x[:40] for x in m_bare + m_tail] or "지적 없음"))
check(abs(gc_bare - gc_tail) < 1e-9, "tail이 GC 계산에 섞이지 않았다",
      f"{gc_bare:.3f} / {gc_tail:.3f}")
# before 재현: 선언창을 Y7 이전 값으로 좁히면 끝단이 red로 돌아간다 = 회귀 감시.
#   Y7에선 `qc_gc_min`을 되돌려 재현했는데, 그 필드가 사라졌으므로 **선언창**으로 한다.
#   기준이 유저 입력에 묶여 있다는 증거이기도 하다([6]과 같은 성질).
old_win = DesignParams(gc_min=35.0, gc_max=60.0)
t_old, _, _ = gc_tier(AT_ANNEAL, 0, params=old_win)
check(t_old == "red", "선언창을 35~60%로 좁히면 끝단이 red가 된다(회귀 감시)", t_old)

print(f"\n[3]~[5] GC 밴드별 등급 — 선언창(20~80%) 하나로 갈린다")
# ⚠️ Y7의 밴드 표에서 28%·70%가 yellow → **green**으로 바뀌었다. 이상 구간(40~60%)
#   경고층을 없앴기 때문이고, 그게 이 작업의 목적이다(유저가 설정한 적 없는 기준).
BANDS = [(15.0, "red", "선언창 밖(너무 낮음)"),
         (28.0, "green", "선언창 안 — Y7에선 yellow였다"),
         (50.0, "green", "선언창 안"),
         (70.0, "green", "선언창 안 — Y7에선 yellow였다"),
         (90.0, "red", "선언창 밖(너무 높음)")]
for pct, want, why in BANDS:
    seq = synth(40, pct)
    got, gc, msgs = gc_tier(seq, 0)
    check(got == want, f"GC {gc:.0%} → {want} ({why})",
          f"got={got}" + (f" / {msgs[0][:48]}" if msgs else ""))

print(f"\n[6] 선언창을 좁히면 지적이 따라 움직인다 (기준이 유저 입력에 묶여 있다)")
# 같은 서열(GC 70%)을 창만 바꿔 두 번 본다 — 2층 시절엔 유저가 창을 바꿔도
# 판정창이 그대로라 등급이 안 움직였다. 그게 Y11이 고친 것이다.
seq70 = synth(40, 70.0)
check(gc_tier(seq70, 0)[0] == "green", "기본창(20~80)에선 green")
narrowed = DesignParams(gc_min=20.0, gc_max=60.0)
t_n, _, m_n = gc_tier(seq70, 0, params=narrowed)
check(t_n == "red", "상한을 60%로 내리면 같은 서열이 red", t_n)
check(any("you set" in x for x in m_n),
      "메시지가 '설정한 기준'임을 밝힌다", (m_n or [""])[0][:58])

print(f"\n[7] 호출자가 창을 직접 주면 그것이 우선한다")
r_called = ValidatePrimer().validate_primer(synth(40, 43.0), min_gc=0.45, params=prm)
check(any("below the floor you set" in m for m in r_called["violations"]),
      "호출자 하한(45%) 밖은 차단된다", f"GC {r_called['gc_content']:.1%}")
check(len([m for m in r_called["violations"] + r_called["warnings"] if "GC" in m]) == 1,
      "GC 지적이 한 번만 나온다 (2층 시절 차단+경고 중복 제거)")

print(f"\n[8] GC 단위 실수(0~1 분수) 방어가 선언창에도 걸린다")
for kw in ({"gc_min": 0.4}, {"gc_max": 0.6}):
    try:
        DesignParams(**kw).validate()
        check(False, f"{kw} → 에러가 나야 함")
    except ValueError as e:
        check("the unit is" in str(e), f"{kw} 차단", str(e).splitlines()[-1][:70])

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
    # ⚠️ Y11 — rev의 GC는 여전히 낮지만(28%대) 그건 **선언창 20~80% 안**이므로
    #   더 이상 경고가 아니다. 숫자는 시트에 그대로 남아 실험자가 볼 수 있다.
    _rev_gc = float(row["rev_gc"])
    check(prm.gc_min <= _rev_gc <= prm.gc_max,
          "rev GC가 선언창 안이다", f"rev_gc={_rev_gc}% ⊆ {prm.gc_min:g}~{prm.gc_max:g}%")
    check("GC%" not in (row.get("warnings") or ""),
          "선언창 안이므로 GC 경고가 없다 (유저가 설정한 적 없는 기준으로 경고하지 않는다)",
          (row.get("warnings") or "경고 없음")[:60])
    check(_rev_gc < 40.0,
          "숫자 자체는 시트에 남아 있다 (신호를 지운 게 아니라 경고에서 뺀 것)",
          f"rev_gc={_rev_gc}%")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y7 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y7 전부 통과")
