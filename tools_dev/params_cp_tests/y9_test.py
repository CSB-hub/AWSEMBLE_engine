"""Y9 검증 — 설계 범위를 못 맞추면 **가장 가까운 값**으로 설계하고 이유를 알린다 (2026-11-05).

## 무엇이 터졌나

유저가 Tm 창을 55~65 → 60~65로 바꿔도 **완전히 동일한 프라이머**가 나왔다.
파라미터는 엔진까지 제대로 도달했는데(Build Log에 찍힘), 원인은 폴백 정책이었다:

요청 범위로 primer3가 후보를 못 내면 `DEFAULT_PARAMS`로 **갈아타서** 재설계했다.
그건 "근접"이 아니라 그냥 기본값이라, **무엇을 요청하든 같은 프라이머**가 나왔다.
실측(BBa_C0083 + V0002):

| 요청 Tm | before | 요청과의 거리 |
|---|---|---|
| 60~65 | 59.4℃ | 0.6 (우연히 가까움) |
| 70~75 | 59.4℃ | 10.6 |
| 45~50 | 59.4℃ | 9.4 — **더 뜨거운** 걸 준다 |

게다가 그 폴백 경고조차 `gibson.py`가 안 꺼내가서 유저에게 안 보였다(Y8에서 해결).

## 해법

**Tm 창만 열고 목표(`OPT_TM`)는 유저 값 그대로** 둔다. primer3의 목적함수가
|Tm − OPT_TM|에 벌점을 주므로 0번 후보가 곧 "요청에 가장 가까운 쌍"이다 —
새 탐색 로직이 필요 없다.

그리고 후보가 0이 되는 흔한 원인은 **쌍 제약**이다: fwd는 창을 만족하는데 rev가
못 미치면 쌍이 깨지고 되는 쪽까지 버려진다. primer3의 EXPLAIN이
`LEFT ok 1 / RIGHT ok 0`으로 정확히 말해주므로 추론하지 않고 그대로 옮긴다.
필요한 결합부 길이도 **실제 탐색으로 구해** 숫자로 알려준다.

## 수용 기준

  1. 요청 창이 다르면 **결과 프라이머가 다르다** (핵심 증상 해소)
  2. 창 안에서 찾을 수 있으면 경고 없이 그대로 (소음 방지)
  3. 못 맞추면 경고에 **요청 창 / 실제 값 / 어느 쪽이 걸렸는지**가 들어간다
  4. 길이 제안이 **실제로 동작하는 숫자**다 — 그 값으로 다시 돌리면 경고가 사라진다
  5. 길이로 불가능하면 그렇다고 말한다(엉뚱한 숫자를 제시하지 않는다)
  6. 근접값은 기본값보다 요청에 **더 가깝다**
  7. construct 서열은 불변이다

실행: python tools_dev/params_cp_tests/y9_test.py
"""
import io
import contextlib
import csv
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from params import DesignParams
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
FAIL = []
db = DatabaseImport()
CDS = db.get_part_by_id("Bioparts", "cds", "BBa_C0083")
VEC = db.get_part_by_id("Bioparts", "vector", "V0002")


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


def run(tag, **kw):
    d = S / f"y9_{tag}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        res = WorkflowBuilder.build_plasmid(
            [CDS], VEC, db, str(d), module_method="gibson", insert_method="gibson",
            params=DesignParams(**kw))
    row = list(csv.DictReader(open(d / "reports" / "primers.csv")))[0]
    con = res.get("construct")
    return {
        "fwd": row["forward_primer"], "rev": row["reverse_primer"],
        "fwd_tm": float(row["fwd_tm"]), "rev_tm": float(row["rev_tm"]),
        "warn": row.get("warnings") or "",
        "sha": hashlib.sha256(str(con.seq).upper().encode()).hexdigest()[:16] if con else None,
    }


# ⚠️ 2026-11-05 Y10 — 탐침 창을 60~65 → **64~68** 로 갈았다.
#   Y9를 쓸 당시 Tm 축은 primer3였고 이 서열의 fwd가 59.4℃라 60~65는 "못 맞추는 창"
#   이었다. Y10에서 축을 NEB Tm Calculator 재현으로 바꾸자 같은 프라이머가 62.3℃가
#   되어 60~65를 **정당하게 만족**한다 — 그러면 55~65와 결과가 같은 게 맞고, 이 창으론
#   더 이상 폴백 경로를 건드리지 못한다. 그래서 지금 축에서 길이창 상한(25nt)으로는
#   못 닿는 64~68로 옮겼다(상한 26nt면 닿는다 → [4]의 길이 제안이 살아있다).
#   테스트가 보는 성질(요청 창이 다르면 결과가 다르다)은 그대로다.
print("[1] 요청 창이 다르면 결과가 다르다 (핵심 증상)")
base = run("base", tm_min=55.0, tm_max=65.0)
hot = run("hot", tm_min=64.0, tm_max=68.0)
cold = run("cold", tm_min=45.0, tm_max=50.0)
print(f"  55~65 → fwd {base['fwd_tm']:.2f}℃   "
      f"64~68 → fwd {hot['fwd_tm']:.2f}℃   45~50 → fwd {cold['fwd_tm']:.2f}℃")
check(hot["fwd"] != base["fwd"], "64~68 요청이 기본창과 다른 프라이머를 낸다",
      f"{base['fwd_tm']:.2f} → {hot['fwd_tm']:.2f}")
check(cold["fwd"] != base["fwd"], "45~50 요청이 기본창과 다른 프라이머를 낸다",
      f"{base['fwd_tm']:.2f} → {cold['fwd_tm']:.2f}")
check(cold["fwd_tm"] < base["fwd_tm"], "낮은 창을 요청하면 **더 낮은** Tm이 나온다",
      f"{cold['fwd_tm']:.2f} < {base['fwd_tm']:.2f}")
check(hot["fwd_tm"] > base["fwd_tm"], "높은 창을 요청하면 **더 높은** Tm이 나온다",
      f"{hot['fwd_tm']:.2f} > {base['fwd_tm']:.2f}")

print("\n[2] 창 안에서 되면 경고 없음 (소음 방지)")
check("가장 가까운 쌍" not in base["warn"], "기본창 실행엔 근접 경고가 없다")
check(55.0 <= base["fwd_tm"] <= 65.0, "실제로 창 안이다", f"{base['fwd_tm']:.2f}℃")

print("\n[3] 못 맞추면 이유가 적힌다")
for tag, r, lo, hi in (("hot", hot, 64, 68), ("cold", cold, 45, 50)):
    w = next((x for x in r["warn"].split(" | ") if "가장 가까운 쌍" in x), "")
    check(bool(w), f"{lo}~{hi}: 근접 경고가 있다")
    check(f"{lo}~{hi}" in w, f"{lo}~{hi}: 요청 창이 적혀 있다")
    check("reverse" in w or "forward" in w or "양쪽" in w,
          f"{lo}~{hi}: 어느 쪽이 걸렸는지 적혀 있다",
          next((s for s in ("reverse", "forward", "양쪽") if s in w), "-"))

print("\n[4] 길이 제안이 실제로 동작하는 숫자다")
w = next(x for x in hot["warn"].split(" | ") if "가장 가까운 쌍" in x)
import re
m = re.search(r"상한을 (\d+) → (\d+)", w)
check(m is not None, "상한 제안이 들어 있다", w[-70:] if not m else f"{m.group(1)} → {m.group(2)}")
if m:
    fixed = run("fixed", tm_min=64.0, tm_max=68.0, primer_len_max=int(m.group(2)))
    check("가장 가까운 쌍" not in fixed["warn"],
          f"제안대로 상한을 {m.group(2)}로 올리면 경고가 사라진다",
          f"fwd {fixed['fwd_tm']:.2f}℃ / rev {fixed['rev_tm']:.2f}℃")
    check(64.0 <= fixed["fwd_tm"] <= 68.0 and 64.0 <= fixed["rev_tm"] <= 68.0,
          "fwd·rev 둘 다 요청 창 안으로 들어온다")

m2 = re.search(r"하한을 (\d+) → (\d+)", next(
    x for x in cold["warn"].split(" | ") if "가장 가까운 쌍" in x))
check(m2 is not None, "낮은 창에는 **하한** 제안이 나온다(방향이 반대다)",
      f"{m2.group(1)} → {m2.group(2)}" if m2 else "-")

print("\n[5] 길이로 불가능하면 그렇다고 말한다")
far = run("far", tm_min=70.0, tm_max=75.0)
w3 = next(x for x in far["warn"].split(" | ") if "가장 가까운 쌍" in x)
check("어디로 바꿔도" in w3, "도달 불가를 명시한다", w3[-60:])
check("상한을" not in w3, "엉뚱한 길이 숫자를 제시하지 않는다")

print("\n[6] 근접값이 기본값보다 요청에 가깝다")
# 기본창 결과(= 옛 폴백이 주던 값)와 비교한다.
for tag, r, lo, hi in (("hot", hot, 64, 68), ("cold", cold, 45, 50)):
    def dist(tm):
        return 0.0 if lo <= tm <= hi else (lo - tm if tm < lo else tm - hi)
    d_new, d_old = dist(r["fwd_tm"]), dist(base["fwd_tm"])
    check(d_new <= d_old, f"{lo}~{hi}: 근접값이 더 가깝다",
          f"{d_new:.2f}℃ <= {d_old:.2f}℃ (옛 폴백)")

print("\n[7] construct 서열은 불변")
shas = {base["sha"], hot["sha"], cold["sha"], far["sha"]}
check(len(shas) == 1, "Tm 창이 달라도 construct는 같다 (프라이머만 바뀐다)", str(shas))

print("\n" + "=" * 72)
if FAIL:
    print(f"Y9 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y9 전부 통과")
