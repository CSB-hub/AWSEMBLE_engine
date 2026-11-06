"""Y11 검증 — 판정 기준을 **유저 선언창 하나**로 합치고, ΔTm(fwd−rev)를 신설했다 (2026-11-05).

## 무엇이 터졌나

Y6(Tm)·Y7(GC)에서 "설계창 바깥에 판정창을 하나 더 둔다"는 처방을 썼다. 그 구조가
2026-11-04부터 **같은 버그를 세 번** 냈고, 세 번째(길이 축)에서는 판정값이 아예
리터럴이었다:

```python
if length < 15:   ...  # ← 유저 설정과 무관한 매직넘버
elif length > 30: ...
```

그래서 유저가 Q5(결합부 상한 26nt)를 골라도 경고는 **30 기준**으로 떴다. 설정과
경고가 따로 노는 상태였고, 실측하면 소음이 지배적이었다 — 프라이머 78행에 경고
154건(행당 2.0건) 중 대부분이 **유저가 설정한 적 없는 기준** 위반이었다:

| 분류 | 건수 | 실체 |
|---|---|---|
| 결합부 길이 15~30bp 밖 | 36 | 하드코딩 리터럴 |
| GC 이상구간 40~60% 밖 | 39 | `qc_gc_ideal_*` (템플릿 성질이라 고칠 수 없다) |
| 결합부 Tm < 45℃ | 22 | `qc_tm_min` (tail에 흡수된 13bp 결합부) |
| tail 흡수 / 올리고 권장 | 33 | 경고가 아니라 설계 설명 |
| template 접근성 | 18 | 주형 고정이라 재설계 불가 |
| 요청 범위 미달(근접값) | 4 | ← 실제로 유저가 고칠 수 있는 유일한 축 |
| hairpin / dimer | 2 | ΔG 하나로 PCR이 안 걸리는 신호 |

그리고 **ΔTm(fwd−rev)가 아예 없었다.** `tm_diff_max`가 있지만 그건 *같은 프라이머*의
결합부↔전체 Tm 차(= GG 프라이머 길이를 정하는 knob)로, 쌍 사이의 차가 아니다.
실측(조성 3 × 경로 2 × 프리셋 3 = 69쌍): 중앙 0.89℃인데 **최대 11.68℃**, 5℃ 초과 17%.
`ta`는 **낮은 쪽 Tm 기준**으로 정해지므로 그 Ta에서 뜨거운 쪽 프라이머가 자기 Tm보다
훨씬 낮은 온도에서 돌아 비특이 결합이 늘어난다.

## 해법

**창을 하나로 합친다.** 판정도 유저가 선언한 값(`tm_*` `gc_*` `primer_len_*`)으로 하고,
판정창이 맡던 "물리적으로 PCR이 성립하지 않는 구간 차단"은 **입력 단계**로 옮긴다
(`TM_INPUT_MIN/MAX` = 35~80℃, GUI·XML·`validate()` 세 곳에 같은 숫자). 창이 하나뿐이면
어긋남이 구조적으로 불가능하다.

길이 축만 예외가 필요했다 — **`tail_len`으로 자를 가른다**:

- `tail_len == 0` (일반 PCR·끝단): 결합부 = 프라이머 전체이고 primer3가 Tm 창 안에서
  길이를 자유롭게 고를 수 있다 → 유저 선언창으로 판정한다.
- `tail_len > 0` (GG·Gibson·yeast): tail이 길이·서열 **둘 다 고정**이라 결합부 길이는
  `tm_diff_max` 수렴의 **종속 변수**다(같은 설정에서도 부품·프리셋에 따라 13~43nt로
  흩어진다). 선언창으로 재면 **유저가 숫자를 어떻게 바꿔도 못 없애는 경고**가 된다.
  → 실제로 막히는 제약만 본다: 합성 주문 상한(`PRIMER_ORDER_MAX_NT`, 전체 프라이머
  기준)과 결합부 특이성 바닥(`amp_min`).

이 분기는 새 규칙이 아니다 — `validate_primer`의 `anneal_gate`가 이미 같은 `tail_len`으로
적/황을 가르고 있었고, 길이 축만 그 규칙에서 빠져 있었다.

## 수용 기준

  1. 판정창 필드(`qc_tm_*` `qc_gc_*` `qc_gc_ideal_*`)가 제거됐다
  2. 세 축(Tm·GC·길이)의 판정이 **유저 선언값에 묶여 있다** — 창을 바꾸면 지적이 따라 온다
  3. Tm 입력 제한이 `validate()`에서 막힌다 (GUI를 거치지 않는 호출의 방어선)
  4. 길이 판정이 `tail_len`으로 갈린다 — tail 없으면 선언창, 있으면 주문 상한·`amp_min`
  5. ΔTm(fwd−rev) 선언창(`pair_tm_diff_max`)이 있고, 초과하면 경고에 Ta 근거가 적힌다
  6. `tm_diff` 컬럼이 주문 시트에 있고 `ta`와 함께 읽을 수 있다
  7. 경고 총량이 실제로 줄었다 — 프리셋 3종 전부에서 "유저가 설정한 적 없는 기준"
     경고(GC 이상구간·결합부 Tm 45℃·길이 리터럴)가 0건
  8. construct 서열은 불변이다 (판정 축만 건드렸다)

실행: python tools_dev/params_cp_tests/y11_test.py
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

from params import (DesignParams, PRIMER_ORDER_MAX_NT, TM_INPUT_MIN,
                    TM_INPUT_MAX)
from basic.primer_design import ValidatePrimer, check_pair_tm_diff
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


prm = DesignParams()

# ═════════════════════════════════════════════════════════════════════════
print("[1] 판정창 필드가 제거됐다 (창은 선언창 하나뿐이다)")
for f in ("qc_tm_min", "qc_tm_max", "qc_gc_min", "qc_gc_max",
          "qc_gc_ideal_min", "qc_gc_ideal_max"):
    check(not hasattr(prm, f), f"{f} 제거됨")
# ΔG 컷오프는 **범위가 아니라 물리량**이라 남는다 — 유저가 선언하는 값이 아니고,
# hairpin/dimer는 선언창과 무관하게 알려야 하는 신호다.
for f in ("qc_dimer_dg", "qc_hairpin_dg_3prime", "qc_hairpin_dg_internal"):
    check(hasattr(prm, f), f"{f}는 남아 있다 (범위가 아니라 ΔG 물리량)")

# ═════════════════════════════════════════════════════════════════════════
print("\n[2] 세 축의 판정이 유저 선언값에 묶여 있다")
# 같은 프라이머를 창만 바꿔 두 번 본다. 2층 시절엔 유저가 창을 바꿔도 판정창이
# 그대로라 등급이 안 움직였다 — 그게 Y11이 고친 것이다.
# 탐침 둘. **기본값에서 깨끗한** 것이 있어야 "창을 바꾸면 지적이 생긴다"를
# 깨끗한 상태와 비교할 수 있다(실측으로 고른 값 — Tm 62.2℃ / GC 40.9% / 22nt).
PR = "ATGACCATGATTACGGATTCAC"              # 22nt, tail 없음, 기본창에서 무경고
PR_LONG = "ATGGCGGATAAAGAACTTAAATTGGCATGC"  # 30nt — 길이창(25)만 넘는 탐침


def msgs(p, axis, primer=PR, **kw):
    v = ValidatePrimer().validate_primer(primer, params=p, **kw)
    return [m for m in v["violations"] + v["warnings"] if axis in m]


_tm = ValidatePrimer().validate_primer(PR, params=prm)["tm"]
_gc = ValidatePrimer().validate_primer(PR, params=prm)["gc_content"] * 100

# Tm 축 — 창을 프라이머 Tm 위로 올리면 지적이 생긴다
hot = DesignParams(tm_min=min(TM_INPUT_MAX, _tm + 4), tm_max=TM_INPUT_MAX)
check(not msgs(prm, "Tm이"), "Tm: 기본창 안이면 조용하다", f"Tm {_tm:.1f}℃")
check(bool(msgs(hot, "Tm is")), "Tm: 창을 올리면 지적된다",
      f"창 {hot.tm_min:g}~{hot.tm_max:g}")

# GC 축 — 선언창을 프라이머 GC 위로 좁히면 지적이 생긴다
gc_narrow = DesignParams(gc_min=20.0, gc_max=max(20.1, _gc - 5))
check(not msgs(prm, "GC%"), "GC: 기본창 안이면 조용하다", f"GC {_gc:.1f}%")
check(bool(msgs(gc_narrow, "GC%")), "GC: 창을 좁히면 지적된다",
      f"창 20~{gc_narrow.gc_max:g}%")

# 길이 축 — 선언창을 좁히면 지적이 생긴다 (리터럴 30이 아니라 설정값을 쓴다)
len_narrow = DesignParams(primer_len_min=18, primer_len_opt=20, primer_len_max=25)
_lm = msgs(len_narrow, "length ceiling", primer=PR_LONG)
check(bool(_lm), f"길이: 상한 25로 좁히면 {len(PR_LONG)}nt가 지적된다",
      (_lm or [""])[0][:58])
check("25 bp" in (_lm or [""])[0],
      "메시지 숫자가 **설정값**이다 (리터럴 30이 아니다)", (_lm or [""])[0][-28:])
check(not msgs(DesignParams(primer_len_max=30), "length ceiling", primer=PR_LONG),
      "상한 30이면 같은 프라이머가 조용하다 (기준이 설정을 따라 움직인다)")

# ═════════════════════════════════════════════════════════════════════════
print("\n[3] Tm 입력 제한 — 판정이 아니라 입력 단계에서 막힌다")
check((TM_INPUT_MIN, TM_INPUT_MAX) == (35.0, 80.0),
      "입력 허용 범위 35~80℃", f"{TM_INPUT_MIN:g}~{TM_INPUT_MAX:g}")
for bad, why in ((dict(tm_min=34.9, tm_max=60.0), "하한이 35℃ 아래"),
                 (dict(tm_min=60.0, tm_max=80.1), "상한이 80℃ 위"),
                 (dict(tm_min=20.0, tm_max=90.0), "양쪽 모두 밖")):
    try:
        DesignParams(**bad).validate()
        check(False, f"{why} → 에러가 나야 함")
    except ValueError as e:
        check("allowed" in str(e), f"{why}: 막힌다",
              str(e).splitlines()[-1][:58])
# 경계값은 통과해야 한다 (off-by-one으로 정상 입력을 막으면 안 된다)
try:
    DesignParams(tm_min=TM_INPUT_MIN, tm_max=TM_INPUT_MAX).validate()
    check(True, "경계값(35~80)은 통과한다")
except ValueError as e:
    check(False, "경계값(35~80)은 통과한다", str(e)[:60])

# ═════════════════════════════════════════════════════════════════════════
print("\n[4] 길이 판정이 tail_len으로 갈린다")
# tail 있는 쪽: 선언창(상한 25)을 한참 넘는 40nt 결합부여도 선언창으로는 안 잡는다.
#   GG 결합부는 `tm_diff_max` 수렴의 종속 변수라 유저가 못 움직이는 축이기 때문이다.
TAIL = "GGTCTCAAATG"          # 11nt (인식서열+spacer+overhang 모사)
LONG_ANNEAL = "ATGGCGGATAAAGAACTTAAATTGGCATGCGGATCCAACA"   # 40nt
tight_len = DesignParams(primer_len_min=18, primer_len_opt=20, primer_len_max=25)
with_tail = msgs(tight_len, "length", primer=TAIL + LONG_ANNEAL, tail_len=len(TAIL))
no_tail = msgs(tight_len, "length", primer=LONG_ANNEAL, tail_len=0)
check(not with_tail, "tail 있으면 선언창 상한으로 지적하지 않는다",
      str([m[:40] for m in with_tail] or "지적 없음"))
check(bool(no_tail), "tail 없으면 같은 40nt가 선언창 상한으로 지적된다",
      (no_tail or [""])[0][:52])

# tail 있는 쪽의 자: 주문 상한(전체 프라이머 기준) + amp_min(결합부 바닥)
check(PRIMER_ORDER_MAX_NT == 80, "주문 상한 80nt", str(PRIMER_ORDER_MAX_NT))
over = ValidatePrimer().validate_primer(
    TAIL + "A" * (PRIMER_ORDER_MAX_NT + 5 - len(TAIL)),
    params=prm, tail_len=len(TAIL))
check(any("synthesis limit" in m for m in over["violations"] + over["warnings"]),
      f"전체 {PRIMER_ORDER_MAX_NT + 5}nt는 주문 상한으로 지적된다")
short = ValidatePrimer().validate_primer(TAIL + "ATGCGGATCCAAC",  # 결합부 13nt
                                         params=prm, tail_len=len(TAIL))
_sm = [m for m in short["warnings"] if "specificity floor" in m]
check(bool(_sm), f"결합부 13nt < amp_min {prm.amp_min}은 지적된다", (_sm or [""])[0][:56])
check("absorbed into the tail" in (_sm or [""])[0],
      "사유가 적혀 있다 (짧은 부품이라 결합부가 tail에 흡수됨)")
ok_tail = ValidatePrimer().validate_primer(TAIL + LONG_ANNEAL,
                                           params=prm, tail_len=len(TAIL))
check(not any("synthesis limit" in m or "specificity floor" in m
              for m in ok_tail["violations"] + ok_tail["warnings"]),
      "제약 안쪽이면 길이 지적이 없다 (상시 경고 제거)")

# ═════════════════════════════════════════════════════════════════════════
print("\n[5] ΔTm(fwd−rev) 선언창")
check(prm.pair_tm_diff_max == 5.0, "pair_tm_diff_max 기본 5℃",
      str(prm.pair_tm_diff_max))
check(check_pair_tm_diff(62.0, 60.0, prm) is None, "창 안(2.0℃)이면 None")
w = check_pair_tm_diff(71.1, 59.4, DesignParams.from_preset("phusion"))
check(w is not None, "창 밖(11.7℃)이면 경고")
check("11.7" in w and "5" in w, "실측값과 설정값이 둘 다 적혀 있다", (w or "")[:56])
check("annealing temperature" in w and "lower one" in w,
      "Ta 근거가 적혀 있다 (왜 문제인지 알 수 있다)")
check("forward" in w, "어느 쪽이 뜨거운지 적혀 있다")
# `tm_diff_max`와 **다른 축**이라는 걸 고정한다 — 이름이 비슷해 합쳐지기 쉽다.
check(prm.tm_diff_max == 5.0 and prm.pair_tm_diff_max == 5.0,
      "두 값은 기본값이 같지만 서로 다른 필드다")
loose = DesignParams(pair_tm_diff_max=15.0)
check(check_pair_tm_diff(71.1, 59.4, loose) is None,
      "pair_tm_diff_max만 올리면 ΔTm 경고가 사라진다 (tm_diff_max와 독립)")
tight_pair = DesignParams(tm_diff_max=0.1)
check(check_pair_tm_diff(62.0, 60.0, tight_pair) is None,
      "tm_diff_max를 조여도 ΔTm 판정은 안 바뀐다 (역방향 독립)")
try:
    DesignParams(pair_tm_diff_max=-1.0).validate()
    check(False, "음수는 에러가 나야 함")
except ValueError as e:
    check("pair_tm_diff_max" in str(e), "음수를 막는다",
          str(e).splitlines()[-1][:52])

# ═════════════════════════════════════════════════════════════════════════
print("\n[6]~[8] E2E — 주문 시트 컬럼 / 경고 총량 / 서열 불변")
db = DatabaseImport()
SPEC = [("promoter", "BBa_J23100"), ("rbs", "BBa_B0032"),
        ("cds", "BBa_C0012"), ("terminator", "BBa_B0010")]
parts = [db.get_part_by_id("Bioparts", t, i) for t, i in SPEC]
vector = db.get_part_by_id("Bioparts", "vector", "pSB1C3")

# 유저가 설정한 적 없는 기준으로 떴던 **옛 문구들** — 이제 하나도 나오면 안 된다.
# ⚠️ 숫자로 잡으면 안 된다 — 새로 생긴 정당한 경고가 같은 숫자를 쓸 수 있다.
#    실제로 `amp_min`(=15) 경고가 "13 bp < 15 bp"라고 찍혀서 처음엔 오검출됐다.
#    옛 판정층은 문구가 달랐으므로(아래) 문구로 잡는 것이 정확하다.
GONE = ("이상 구간보다",                      # qc_gc_ideal_* (39건)
        "Tm이 너무 낮습니다",                  # qc_tm_min 45℃ (22건)
        "Tm이 너무 높습니다",                  # qc_tm_max 85℃
        "프라이머(anneal)가 너무 깁니다",       # 길이 리터럴 30 (36건)
        "프라이머(anneal)가 너무 짧습니다")     # 길이 리터럴 15
shas = {}
for preset in ("taq", "phusion", "q5"):
    d = S / f"y11_{preset}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        res = WorkflowBuilder.build_plasmid(
            parts, vector, db, str(d), module_method="golden_gate",
            insert_method="gibson", params=DesignParams.from_preset(preset))
    rows = list(csv.DictReader(open(d / "reports" / "primers.csv")))
    hdr = rows[0].keys()

    # [6] 컬럼
    check("tm_diff" in hdr, f"{preset}: primers.csv에 tm_diff 컬럼이 있다")
    check(list(hdr).index("tm_diff") == list(hdr).index("ta") + 1,
          f"{preset}: tm_diff가 ta 바로 옆이다 (같이 읽어야 뜻이 통한다)")
    _pairs = [(r["ta"], r["tm_diff"]) for r in rows
              if (r.get("forward_primer") or "").strip()]
    check(all(a and t for a, t in _pairs),
          f"{preset}: 프라이머가 있는 행은 ta·tm_diff가 둘 다 채워진다",
          f"{len(_pairs)}행")
    _oligo = [r for r in rows if not (r.get("forward_primer") or "").strip()]
    check(all(not (r["tm_diff"] or "").strip() for r in _oligo),
          f"{preset}: 올리고 행은 tm_diff가 비어 있다", f"{len(_oligo)}행")
    # 숫자가 실제로 |fwd−rev|인가 (컬럼만 있고 값이 엉뚱하면 더 나쁘다)
    for r in rows:
        if not (r["tm_diff"] or "").strip():
            continue
        want = abs(float(r["fwd_tm"]) - float(r["rev_tm"]))
        check(abs(float(r["tm_diff"]) - want) < 0.051,
              f"{preset}/{r['part_name'][:16]}: tm_diff = |fwd−rev|",
              f"{r['tm_diff']} ≈ {want:.1f}")

    # [7] 경고 총량 — 사라져야 할 문구가 하나도 없다
    allw = " | ".join(r.get("warnings") or "" for r in rows)
    for g in GONE:
        check(g not in allw, f"{preset}: '{g}' 경고가 없다",
              next((x[:44] for x in allw.split(" | ") if g in x), ""))
    # 반대로 **남아야 하는** 경고도 확인한다 — 소음을 지우면서 신호까지 지우면
    # 더 나쁘다. BBa_B0032는 짧아서 결합부가 tail에 흡수돼 13nt가 되고, 그건
    # `amp_min`(유저 설정값) 미달이라 알려야 한다.
    check("specificity floor" in allw,
          f"{preset}: amp_min 미달 경고는 남아 있다 (신호를 지운 게 아니다)",
          next((x[:52] for x in allw.split(" | ") if "specificity floor" in x), "없음"))
    n_warn = sum(len([x for x in (r.get("warnings") or "").split(" | ") if x.strip()])
                 for r in rows)
    print(f"     {preset}: 프라이머 {len(rows)}행 · 경고 {n_warn}건 "
          f"(행당 {n_warn / max(len(rows), 1):.2f})")

    # [8] 서열
    con = res.get("construct")
    shas[preset] = hashlib.sha256(str(con.seq).upper().encode()).hexdigest()[:16] if con else None

check(len(set(shas.values())) == 1,
      "프리셋이 달라도 construct는 같다 (판정·Tm 축만 건드렸다)", str(shas))
check(all(v for v in shas.values()), "세 프리셋 전부 construct 산출")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y11 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y11 전부 통과")
