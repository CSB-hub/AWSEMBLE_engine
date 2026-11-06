"""H1 검증 — 하이브리드 열역학 (≤60bp primer3 / >60bp ViennaRNA-DNA, 경고 전용).

수용 기준
  1. ≤60bp 경로가 **완전히 불변**이다 (엔진 primer3, 모든 키·값 동일)
  2. >60bp가 예외 없이 값을 돌려준다 (기존에는 primer3가 RuntimeError)
  3. >60bp는 **게이트를 움직이지 않는다** — has_strong_structure=False, severity='ok'
  4. engine/warn_only 키로 어느 엔진이 계산했는지 구분된다
  5. 60bp 경계가 정확하다 (60=primer3, 61=viennarna)
  6. 강한 구조(완전 상보 헤어핀)를 실제로 검출한다 (ΔG가 충분히 음수)
  7. 같은 입력 → 같은 결과 (재현성)
  8. `validate_primer`가 >60bp 프라이머에서 더 이상 터지지 않는다
실행: python tools_dev/params_cp_tests/h1_test.py
"""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

import primer3

from basic.secondary_structure import (
    STRUCT_DNA_CONC,
    STRUCT_DNTP_CONC,
    STRUCT_DV_CONC,
    STRUCT_MV_CONC,
    STRUCT_TEMP_C,
    analyze_long_duplex,
    analyze_single_strand_structure,
)
from params import THERMO_MAX_LEN

FAIL = []


def check(cond, msg):
    print(("  ✅ " if cond else "  ❌ ") + msg)
    if not cond:
        FAIL.append(msg)


_COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(_COMP)[::-1]


def seeded(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


print("=" * 72)
print(f"H1 — 하이브리드 열역학 (경계 {THERMO_MAX_LEN}bp)")
print("=" * 72)

print("\n[1] ≤60bp 경로 불변 (primer3)")
short_cases = [seeded(n, 900 + n) for n in (18, 25, 40, 55, 60)]
for s in short_cases:
    r = analyze_single_strand_structure(s)
    # primer3를 직접 불러 같은 값이 나오는지 대조한다 — 하이브리드 분기가 끼어들어
    # 짧은 경로의 숫자를 바꾸지 않았음을 보인다.
    # ⚠️ 염/농도/온도를 모듈과 **같은 조건**으로 줘야 한다. primer3 기본값으로 부르면
    #    ΔG가 달라져서 "코드가 바뀐 것처럼" 보인다(실제로는 조건이 다른 것).
    direct = primer3.calc_hairpin(
        s, mv_conc=STRUCT_MV_CONC, dv_conc=STRUCT_DV_CONC,
        dntp_conc=STRUCT_DNTP_CONC, dna_conc=STRUCT_DNA_CONC,
        temp_c=STRUCT_TEMP_C, output_structure=True,
    )
    expect_dg = (direct.dg / 1000.0) if direct.structure_found else 0.0
    ok = (r["engine"] == "primer3" and abs(r["delta_g"] - expect_dg) < 1e-9
          and r["warn_only"] is False)
    print(f"  {len(s):>3}bp  engine={r['engine']:<10} ΔG {r['delta_g']:>7.2f} "
          f"(primer3 직접 {expect_dg:>7.2f})  severity={r['severity']}")
    check(ok, f"{len(s)}bp → primer3 경로, ΔG 일치")

print("\n[2] >60bp가 예외 없이 값을 돌려준다")
long_cases = [seeded(n, 700 + n) for n in (61, 80, 120, 200)]
for s in long_cases:
    # 기존 거동 확인: primer3는 이 길이를 거부한다
    raised = False
    try:
        primer3.calc_hairpin(s, temp_c=37.0)
    except Exception:
        raised = True
    r = analyze_single_strand_structure(s)
    print(f"  {len(s):>3}bp  primer3 거부={raised}  engine={r['engine']:<10} "
          f"ΔG {r['delta_g']:>8.2f}  구조 페어 {r['structure'].count('(')}")
    check(raised, f"{len(s)}bp는 primer3가 거부함(하이브리드가 필요한 구간)")
    check(r["engine"] == "viennarna", f"{len(s)}bp → ViennaRNA 경로")
    check(isinstance(r["delta_g"], float), f"{len(s)}bp ΔG 값 반환")

print("\n[3] >60bp는 게이트를 움직이지 않는다 (경고 전용)")
for s in long_cases:
    r = analyze_single_strand_structure(s)
    check(r["has_strong_structure"] is False,
          f"{len(s)}bp has_strong_structure=False")
    check(r["severity"] == "ok", f"{len(s)}bp severity='ok'")
    check(r["warn_only"] is True, f"{len(s)}bp warn_only=True")
    check("경고 전용" in r.get("note", ""), f"{len(s)}bp note에 경고 전용 명시")

print("\n[4] 키 스키마 호환")
keys_short = set(analyze_single_strand_structure(short_cases[0]))
keys_long = set(analyze_single_strand_structure(long_cases[0]))
common = {"delta_g", "tm", "structure_found", "involves_3prime", "structure",
          "severity", "hairpins", "has_strong_structure", "engine", "warn_only"}
check(common <= keys_short, f"짧은 경로가 공통 키를 모두 가짐 (누락 {common - keys_short})")
check(common <= keys_long, f"긴 경로가 공통 키를 모두 가짐 (누락 {common - keys_long})")

print("\n[5] 경계 정확성")
r60 = analyze_single_strand_structure(seeded(60, 11))
r61 = analyze_single_strand_structure(seeded(61, 11))
check(r60["engine"] == "primer3", f"{THERMO_MAX_LEN}bp → primer3")
check(r61["engine"] == "viennarna", f"{THERMO_MAX_LEN + 1}bp → viennarna")

print("\n[6] 강한 구조 검출")
# 완전 상보 스템 + 루프 → 매우 안정한 헤어핀이어야 한다.
stem = seeded(35, 77)
hairpin = stem + "AAAA" + rc(stem)          # 74bp
flat = "A" * 74                              # 구조 없음
rh = analyze_single_strand_structure(hairpin)
rf = analyze_single_strand_structure(flat)
print(f"  완전상보 헤어핀 74bp → ΔG {rh['delta_g']:.2f}, 페어 {rh['structure'].count('(')}")
print(f"  polyA 74bp        → ΔG {rf['delta_g']:.2f}, 페어 {rf['structure'].count('(')}")
check(rh["delta_g"] < rf["delta_g"] - 10, "헤어핀이 무구조보다 훨씬 안정하게 계산됨")
check(rh["structure_found"] is True, "헤어핀에서 구조 검출")
check(rh["structure"].count("(") >= 25, "스템 염기쌍이 다수 검출됨")

print("\n[7] 재현성")
a = analyze_single_strand_structure(long_cases[1])
b = analyze_single_strand_structure(long_cases[1])
check(a == b, "같은 입력 → 같은 결과 dict")

print("\n[9] dimer 축도 >60bp를 계산한다 (경고 전용)")
from basic.primer_design import ValidatePrimer as _VP

a = seeded(70, 55)
comp = rc(a)                                  # 완전 상보 → 강한 duplex
non = "A" * 70                                # 비상보 대조군
d_strong = _VP.check_primer_dimer(a, comp)
d_weak = _VP.check_primer_dimer(non, non)
print(f"  상보 70/70nt  → engine={d_strong.get('engine')} ΔG {d_strong['delta_g']:.2f}")
print(f"  비상보 70/70nt → engine={d_weak.get('engine')} ΔG {d_weak['delta_g']:.2f}")
check(d_strong.get("engine") == "viennarna", ">60bp dimer가 ViennaRNA로 계산됨")
check(d_strong["delta_g"] < d_weak["delta_g"] - 10,
      "상보 쌍이 비상보보다 훨씬 안정하게 계산됨")
check(d_strong["has_dimer"] is False, ">60bp dimer는 게이트를 움직이지 않음")
check("경고 전용" in d_strong.get("warning", ""), "경고 문구에 전용 표시")
# ≤60bp는 primer3 경로 그대로
d_short = _VP.check_primer_dimer(seeded(25, 1), seeded(25, 2))
check("engine" not in d_short or d_short.get("engine") != "viennarna",
      "≤60bp dimer는 primer3 경로 유지")
check(analyze_long_duplex("ACGT" * 20, "ACGT" * 20) is not None,
      "analyze_long_duplex 직접 호출 가능")

print("\n[8] validate_primer가 >60bp에서 터지지 않는다")
from basic.primer_design import ValidatePrimer

long_primer = seeded(72, 31)
try:
    v = ValidatePrimer.validate_primer(long_primer, tail_len=45)
    print(f"  72nt 프라이머 검증 → passed={v.get('passed')} "
          f"violations={len(v.get('violations') or [])} "
          f"warnings={len(v.get('warnings') or [])}")
    check(True, "예외 없이 검증 완료")
    check(not any("hairpin" in str(x).lower() for x in (v.get("violations") or [])),
          "긴 프라이머 hairpin이 차단(violation)으로 올라가지 않음")
except Exception as exc:
    check(False, f"예외 발생: {type(exc).__name__}: {exc}")

print("\n" + "=" * 72)
if FAIL:
    print(f"H1 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("H1 전부 통과")
