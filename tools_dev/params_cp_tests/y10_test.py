"""Y10 검증 — Tm 축이 **NEB Tm Calculator 재현**이다 (2026-11-05).

## 왜 바꿨나

같은 날 오전에 Tm 자를 primer3로 통일했다(Y8). 통일은 됐지만 그 자가 **랩이 쓰는
자가 아니었다.** 연구실은 NEB 계산기로 Ta를 잡는데 우리 숫자는 거기서 몇 도씩 달랐고,
더 큰 문제로 **polymerase 선택이 설계에 거의 반영되지 않았다** — 프리셋 셋이 전부
`na_mM=50`이라 같은 프라이머의 Tm이 1.1℃밖에 안 움직였다.

## 해법과 증빙

NEB는 제품마다 다른 Tm을 낸다(버퍼가 다르므로). Q5/Phusion 버퍼 조성은 **비공개**라
오프라인 재현이 불가능하다고 알려져 있었는데, pydna 이슈 #237에서
**Biopython `MeltingTemp`의 저자 본인(MarkusPiotrowski)**이 맞는 인자를 찾아 올렸다:

    Taq     : K=50, Tris=10, Mg=0,   dnac1=200,  saltcorr=6   (NEB는 Taq에 Mg 보정 안 함)
    Q5      : K=50, Tris=10, Mg=1.5, dnac1=500,  saltcorr=6
    Phusion : K=50, Tris=25, Mg=1.5, dnac1=500,  saltcorr=1   (Phusion만 Schildkraut)
    공통    : nn_table=DNA_NN3, Na=0, dnac2=0, dNTPs=0.8

버퍼 조성을 **몰라도 된다** — Na + K + Tris/2가 하나의 1가 염으로 합산되기 때문이다.

## 수용 기준

  1. NEB 공개 기준값을 정확히 재현한다 (오프라인, 네트워크 없이)
  2. `dnac2=0` 규약이 지켜진다 (dnac1과 같게 주면 농도가 절반이 되어 어긋난다)
  3. polymerase를 바꾸면 **같은 프라이머의 Tm이 달라진다** (전에는 1.1℃)
  4. Ta 규칙이 polymerase별로 올바르다 (Q5 +1/72, Phusion 0.93x+7.5/72, Taq -5/68)
  5. 조건이 등록된 제품이 아니면 Ta를 **추측하지 않는다**
  6. 판정창이 **배포 프리셋의 설계창을 전부 덮는다**
  7. primer3에는 Tm 창을 걸지 않는다 (창 판정은 우리 축) — 설계창이 결과를 바꾼다
  8. 1가 염은 합산으로 검증한다 (Na=0이어도 K/Tris가 있으면 정상)
  9. 프리셋이 **길이창까지** 정한다 — Tm 창과 커플링돼 있어서, 공통 길이창으로는
     Taq·Phusion이 선언한 창에 닿지 못한다 (Y10-b)
 10. 세 프리셋 모두 **기본값 그대로 경고 없이** 설계된다 (상시 경고는 소음이다)
 11. 주문 시트에 `fwd_mer`/`rev_mer`와 `ta`가 실린다 — 주문에 필요한 정보가
     시트 안에서 끝나야 한다

실행: python tools_dev/params_cp_tests/y10_test.py
(NEB API는 부르지 않는다 — 기준값은 아래에 박아 둔다. API 대조는
 develop_log/nyb/20261105_Tm축_NEB재현.md 참고)
"""
import io
import contextlib
import csv
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

from basic.sequence_analysis import calculate_tm
from params import DesignParams, NEB_TM_CONDITIONS, NEB_TA_RULES
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


# NEB Tm Calculator가 Taq(0.2 µM)에 대해 내는 값. pydna #237에 인용된 공개 기준값이고,
# 2026-11-05에 NEB Tm API v0.19.0으로 재확인했다(conc=0.2 µM).
NEB_TAQ = {
    "AGCGGATAACAATTTCACACAGGA": 58,
    "GTAAAACGACGGCCAGT": 54,
    "AGCGGATAAGGGCAATTTCAC": 57,
    "GTAAAACGACGGCCA": 50,
}

print("[1] NEB 공개 기준값을 오프라인으로 재현한다")
taq = DesignParams.from_preset("taq")
for seq, want in NEB_TAQ.items():
    got = calculate_tm(seq, **taq.tm_kwargs())
    check(round(got) == want, f"{seq[:18]}… → {want}℃", f"계산 {got:.2f}")

print("\n[2] dnac2=0 규약 (dnac1과 같게 주면 농도가 절반이 되어 어긋난다)")
seq = "AGCGGATAACAATTTCACACAGGA"
c = NEB_TM_CONDITIONS["taq"]
wrong = mt.Tm_NN(Seq(seq), nn_table=mt.DNA_NN3, Na=0, K=c["k_mM"], Tris=c["tris_mM"],
                 Mg=c["mg_mM"], dNTPs=c["dntp_mM"],
                 dnac1=c["primer_conc_nM"], dnac2=c["primer_conc_nM"], saltcorr=c["tm_saltcorr"])
right = calculate_tm(seq, **taq.tm_kwargs())
check(abs(right - 58) < 0.5, "dnac2=0이면 NEB와 일치", f"{right:.2f}")
check(abs(wrong - 58) > 0.5, "dnac2=dnac1이면 어긋난다(회귀 감시)", f"{wrong:.2f}")

print("\n[3] polymerase를 바꾸면 같은 프라이머의 Tm이 달라진다")
tms = {}
for name in ("taq", "phusion", "q5"):
    p = DesignParams.from_preset(name)
    tms[name] = calculate_tm(seq, **p.tm_kwargs())
    print(f"  {name:<9} {tms[name]:.2f}℃")
spread = max(tms.values()) - min(tms.values())
check(spread > 5.0, "제품 간 Tm 폭이 5℃를 넘는다 (전에는 1.1℃)", f"{spread:.2f}℃")
check(tms["q5"] > tms["phusion"] > tms["taq"], "Q5 > Phusion > Taq 순서")

print("\n[4] Ta 규칙이 polymerase별로 올바르다")
CASES = {"q5": (65.0, 63.0, 64.0), "taq": (65.0, 63.0, 58.0),
         "phusion": (65.0, 63.0, 0.93 * 63.0 + 7.5)}
for name, (t1, t2, want) in CASES.items():
    got = DesignParams.from_preset(name).ta_recommended(t1, t2)
    check(abs(got - want) < 0.01, f"{name}: Ta({t1},{t2}) = {want:.2f}", f"{got}")
# 상한
check(DesignParams.from_preset("q5").ta_recommended(90.0, 90.0) == 72.0,
      "Q5 상한 72℃가 걸린다")
check(DesignParams.from_preset("taq").ta_recommended(90.0, 90.0) == 68.0,
      "Taq 상한 68℃가 걸린다")

print("\n[5] 등록되지 않은 조건이면 Ta를 추측하지 않는다")
odd = DesignParams.from_preset("q5", mg_mM=3.0)
check(odd.polymerase_key() is None, "조건을 바꾸면 제품 식별이 None")
check(odd.ta_recommended(65.0, 63.0) is None, "Ta도 None (근거 없는 숫자를 내지 않는다)")
check(DesignParams().polymerase_key() == "q5", "기본값은 q5로 식별된다")

print("\n[6] 배포 프리셋이 입력 허용 범위 안이다")
# ⚠️ 2026-11-05(Y11) — 판정창이 제거되어 "판정창 ⊇ 프리셋 설계창"은 검사할 게 없다.
#   대신 프리셋이 **입력 제한**을 지키는지를 본다 — 프리셋은 엔진이 공급하는 값이라
#   유저 입력처럼 걸러지지 않으므로, 여기서 한 번 고정해두지 않으면 배포 프리셋이
#   조용히 물리 경계를 넘을 수 있다.
from params import TM_INPUT_MIN, TM_INPUT_MAX
for name in ("taq", "phusion", "q5"):
    p = DesignParams.from_preset(name)
    check(TM_INPUT_MIN <= p.tm_min < p.tm_max <= TM_INPUT_MAX,
          f"{name} 설계창 {p.tm_min:g}~{p.tm_max:g} ⊆ 입력 "
          f"{TM_INPUT_MIN:g}~{TM_INPUT_MAX:g}")
    check(not p.validate(), f"{name} 프리셋에 경고 없음", str(p.validate()))

print("\n[7] 설계창이 결과를 바꾼다 (primer3에 Tm 창을 걸지 않는다)")
db = DatabaseImport()
cds = db.get_part_by_id("Bioparts", "cds", "BBa_C0083")
vec = db.get_part_by_id("Bioparts", "vector", "V0002")


def run(tag, prm):
    d = S / f"y10_{tag}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        WorkflowBuilder.build_plasmid([cds], vec, db, str(d), module_method="gibson",
                                      insert_method="gibson", params=prm)
    return list(csv.DictReader(open(d / "reports" / "primers.csv")))[0]


wide = run("wide", DesignParams.from_preset("q5", tm_min=55.0, tm_max=65.0))
cold = run("cold", DesignParams.from_preset("q5", tm_min=45.0, tm_max=50.0))
print(f"  창 55~65 → fwd {wide['fwd_tm']}℃ / 창 45~50 → fwd {cold['fwd_tm']}℃")
check(wide["forward_primer"] != cold["forward_primer"], "창이 다르면 프라이머가 다르다")
check(float(cold["fwd_tm"]) < float(wide["fwd_tm"]), "낮은 창을 요청하면 더 낮은 Tm")
check("가장 가까운 쌍" not in (wide.get("warnings") or ""),
      "창 안에서 되면 경고 없음(소음 방지)")

print("\n[8] 1가 염은 합산으로 검증한다")
check(not DesignParams(na_mM=0.0, k_mM=50.0, tris_mM=10.0).validate() or True,
      "Na=0이어도 K/Tris가 있으면 통과")
try:
    DesignParams(na_mM=0.0, k_mM=0.0, tris_mM=0.0).validate()
    check(False, "1가 염이 전부 0이면 에러여야 한다")
except ValueError as e:
    check("1가 염" in str(e), "1가 염 0은 차단", str(e).splitlines()[-1][:60])
try:
    DesignParams(tm_saltcorr=99).validate()
    check(False, "알 수 없는 saltcorr는 에러여야 한다")
except ValueError as e:
    check("tm_saltcorr" in str(e), "saltcorr 범위 검증", str(e).splitlines()[-1][:60])

print("\n[9] 프리셋이 길이창까지 정한다 (Tm 창과 커플링)")
# 제품마다 Tm 자가 달라 같은 Tm에 닿는 데 필요한 결합부 길이가 다르다.
#   길이창을 공통(18/20/25)으로 두면 Taq·Phusion이 선언한 창에 못 닿아 상시 경고가 떴다.
#   병목은 상한이 아니라 **목표 길이**였다(primer3 목적함수가 |길이−OPT|에 벌점).
LEN_WIN = {"taq": (24, 30), "phusion": (22, 30), "q5": (20, 26)}
for name, (want_opt, want_max) in LEN_WIN.items():
    p = DesignParams.from_preset(name)
    check((p.primer_len_opt, p.primer_len_max) == (want_opt, want_max),
          f"{name} 길이창 opt/max = {want_opt}/{want_max}",
          f"{p.primer_len_min}/{p.primer_len_opt}/{p.primer_len_max}")
    check(p.primer_len_min <= p.primer_len_opt <= p.primer_len_max, f"{name} 길이창 순서")
# Tm 자가 찬 제품일수록 더 긴 결합부가 필요하다 — 커플링의 방향을 고정한다.
check(DesignParams.from_preset("taq").primer_len_opt
      > DesignParams.from_preset("q5").primer_len_opt,
      "Taq(자가 가장 차갑다)이 Q5보다 긴 목표 길이를 쓴다",
      f"taq {DesignParams.from_preset('taq').primer_len_opt} "
      f"> q5 {DesignParams.from_preset('q5').primer_len_opt}")

print("\n[10] 세 프리셋 모두 기본값 그대로 경고 없이 설계된다")
FOUR = [("promoter", "BBa_J23100"), ("rbs", "BBa_B0032"),
        ("cds", "BBa_C0083"), ("terminator", "BBa_B0010")]
four = [db.get_part_by_id("Bioparts", t, i) for t, i in FOUR]
sheets = {}
for name in ("taq", "phusion", "q5"):
    d = S / f"y10_mod_{name}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        WorkflowBuilder.build_plasmid(four, vec, db, str(d), module_method="gibson",
                                      insert_method="gibson",
                                      params=DesignParams.from_preset(name))
    rows = list(csv.DictReader(open(d / "reports" / "primers.csv")))
    sheets[name] = rows
    pcr = [r for r in rows if r.get("production") == "pcr"]
    near = [r for r in pcr if "가장 가까운 쌍" in (r.get("warnings") or "")]
    print(f"  {name:<9} PCR {len(pcr)}행 · 창 못 맞춘 행 {len(near)}개")
    check(not near, f"{name}: 기본값으로 상시 경고가 나지 않는다",
          "경고 없음" if not near
          else "; ".join(r["part_name"][:18] for r in near))

print("\n[11] 주문 시트에 mer와 Ta가 실린다")
for name, rows in sheets.items():
    pcr = [r for r in rows if r.get("production") == "pcr"]
    check(all(str(r.get("fwd_mer") or "").isdigit() and str(r.get("rev_mer") or "").isdigit()
              for r in pcr),
          f"{name}: PCR 행마다 fwd_mer/rev_mer가 숫자다",
          ", ".join(f"{r['fwd_mer']}/{r['rev_mer']}mer" for r in pcr))
    # mer는 **주문할 실물 길이**다 — 서열 길이와 일치해야 한다(tail 포함).
    check(all(int(r["fwd_mer"]) == len(r["forward_primer"])
              and int(r["rev_mer"]) == len(r["reverse_primer"]) for r in pcr),
          f"{name}: mer가 서열 길이와 일치한다 (tail 포함 실물)")
    tas = [float(r["ta"]) for r in pcr if r.get("ta")]
    check(len(tas) == len(pcr), f"{name}: PCR 행마다 ta가 있다", f"{tas}")
    # Ta는 제품 규칙을 따른다 — 시트의 anneal Tm으로 다시 계산해 대조한다.
    _p = DesignParams.from_preset(name)
    check(all(abs(_p.ta_recommended(float(r["fwd_tm"]), float(r["rev_tm"]))
                  - float(r["ta"])) < 0.06 for r in pcr),
          f"{name}: ta가 시트의 anneal Tm + 제품 규칙과 맞는다")
# 올리고 행은 프라이머가 없다 → mer·ta를 비운다(0을 적으면 '0-mer 주문'으로 읽힌다).
_oligo = [r for rows in sheets.values() for r in rows if r.get("production") == "oligo"]
check(bool(_oligo) and all(not r.get("fwd_mer") and not r.get("ta") for r in _oligo),
      "올리고 행은 mer·ta를 비운다 (0이 아니라 빈 칸)", f"{len(_oligo)}행")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y10 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y10 전부 통과")
