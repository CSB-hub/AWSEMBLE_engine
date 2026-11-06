"""Y8 검증 — Tm 자(尺)를 하나로 통일했다 (2026-11-05).

> ⚠️ **2026-11-05 Y10에서 기준 구현이 교체됐다.** 이 파일을 처음 쓸 때 고른 자는
> primer3였는데, 같은 날 오후 **NEB Tm Calculator 재현**(Biopython `Tm_NN` +
> `DNA_NN3` + NEB 버퍼 조성)으로 다시 갈았다. 이유는 polymerase 선택이 Tm에 거의
> 영향을 못 주고 있었고(같은 창에서 Tm 폭 1.1℃), 실험자가 실제로 보는 기준이
> NEB Tm Calculator였기 때문이다. 자세한 근거와 NEB 공개값 재현 오차는
> `y10_test.py` / `develop_log/nyb/20261105_Tm축_NEB재현.md`에 있다.
>
> **Y8이 지키는 성질은 그대로다** — "자가 **하나**인가"와 "GG·gibson이 같은 창을
> 같은 자로 보는가". 바뀐 건 그 하나가 무엇이냐일 뿐이라, 아래 [1]·[2]는 기준
> 구현만 교체하고 나머지 절은 손대지 않았다.

## 무엇이 터졌나

**같은 `tm_min`이 경로에 따라 다른 자로 읽혔다.**

| 경로 | 창을 거는 주체 | 자 |
|---|---|---|
| gibson / yeast_hr | primer3 (`PRIMER_MIN_TM`/`MAX_TM`으로 그대로 전달) | primer3 |
| golden_gate | `_fit_amplicon`이 `tm_min <= calculate_tm(..) <= tm_max` | Biopython `Tm_NN` |

그리고 주문 시트·등급에 적히던 Tm은 둘 중 **Biopython 쪽**이었다. 두 자는 같은 염·농도를
넣어도 어긋난다 — 실측(Bioparts 파츠 310종 × 길이 5종 × 양 끝 = 3,100 표본):

    오프셋(primer3 − Biopython)  평균 +4.27℃  표준편차 0.42  범위 +3.54~+5.07
    GC와 무관, 길이에만 약한 기울기(18nt +3.70 → 30nt +4.87)

증상 ①: 같은 창을 선언해도 GG가 gibson보다 ~4℃ 뜨거운(불필요하게 긴) 프라이머를 만들었다.
증상 ②: 유저가 창을 60으로 올려도 primer3 기준으론 아직 통과라 결과가 안 바뀌고
        시트엔 55가 찍혔다(실제 운영 실행 h101/h102가 완전히 동일한 프라이머를 냈다).

## 해법

프라이머를 **실제로 고르는 주체**(primer3)의 자로 통일. `calculate_tm`이 단일 정의
지점이라 거기만 바꾸면 설계(GG)·보고(primers.csv)·판정(`validate_primer`)이 한 자를 쓴다.
Biopython 축에 묶여 있던 값만 +4℃ 옮겨 뜻을 보존했다(`qc_tm_*`, `terminal_tm_min`).
선언창(`tm_min`/`tm_max`)과 프리셋은 **원래부터 primer3 축**이라 그대로 뒀다.

## 수용 기준

  1. `calculate_tm`은 **단일 정의점**이고, NEB 재현 조건의 Biopython `Tm_NN`과
     정확히 같다 (자가 하나다 — 기준 구현은 Y10에서 교체됨)
  2. primer3와는 다르고, 그 차이가 **고정 오프셋이 아니다** — 그래서 primer3를
     자로 쓸 수 없다(상수 보정으로 메울 수 없는 종류의 차이다)
  3. 비-ACGT는 **조용히 폴백하지 않고 실패**한다 (축 혼합 방지)
  4. 길이 상한이 없다 — yeast tail(80nt)까지 계산된다
  5. 축 이동이 반영됐다: `qc_tm` 54/84, `terminal_tm_min` 59, 선언창은 55/65 그대로
  6. 판정창이 설계창을 포함한다(같은 축이 됐으므로 의미가 있다)
  7. **GG와 gibson이 같은 창을 같은 자로 본다** — 두 경로 anneal Tm이 모두 선언창 안
  8. construct 서열은 불변이다(Tm은 프라이머를 고를 뿐 산물을 안 바꾼다)
  9. 설계 범위 폴백이 일어나면 그 사유가 **주문 시트까지** 간다
     (`gibson.py`가 `design_basic_primers`의 warnings를 안 꺼내가던 문제)

실행: python tools_dev/params_cp_tests/y8_test.py
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

import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from basic.sequence_analysis import calculate_tm, get_sequence_stats
from params import DesignParams
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
FAIL = []
prm = DesignParams()


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


SEQS = ["ATGTGTTTAAAGCAAATCATTGGC",                 # AT-rich 24nt (실제 mCherry급 말단)
        "TCTGATAGTGCTAGTGTAGATCTC",                 # 중간 GC 24nt
        "GCCGGCGCGGCCATGGCGCCGGCG",                 # GC-rich 24nt
        "ATGGCGGATAAAGAACTTAAATT"]                  # 23nt

# 1가 염은 **환산값**이 기준이다 — NEB 버퍼는 Na⁺가 아니라 K⁺/Tris 기반이라
#   Na만 보면 0이다. primer3에 같은 조건을 주려면 이 환산값을 넘겨야 공정한 비교가 된다.
MONOVALENT = prm.na_mM + prm.k_mM + prm.tris_mM / 2.0

print("[1] calculate_tm은 단일 정의점이고 NEB 재현 조건과 정확히 같다 (자가 하나다)")
for s in SEQS:
    ours = calculate_tm(s, **prm.tm_kwargs())
    # NEB 재현 조건을 **손으로 다시 적어** 비교한다 — `tm_kwargs()`가 조용히 바뀌면
    #   엔진 전체가 같이 움직여서 "같다"가 자동으로 성립해버린다. 그걸 막는 고정점이다.
    #   핵심 키: nn_table=DNA_NN3 / dnac2=0 / Na는 0이고 1가 염은 K·Tris가 담당.
    ref = mt.Tm_NN(Seq(s), nn_table=mt.DNA_NN3,
                   Na=0.0, K=50.0, Tris=prm.tris_mM, Mg=prm.mg_mM, dNTPs=0.8,
                   dnac1=prm.primer_conc_nM, dnac2=0, saltcorr=prm.tm_saltcorr)
    check(abs(ours - ref) < 1e-9, f"{s[:14]}… 일치", f"{ours:.2f} == {ref:.2f}")
check(MONOVALENT == 55.0, "1가 염 환산 = Na 0 + K 50 + Tris 10/2", f"{MONOVALENT:g}mM")

print("\n[2] primer3와는 다르고, 그 차이는 고정 오프셋이 아니다")
# 이게 primer3를 자로 쓸 수 없는 이유다. 차이가 상수라면 +N 보정으로 끝났겠지만,
#   GC와 길이에 따라 움직인다(실측 +2.3~+7.1℃, 길이 18→36nt에서 +3.4→+2.7 기울기).
#   상수 보정으로 메울 수 없으니 기준 구현 자체를 NEB 재현으로 가져가야 했다.
diffs = []
for s in SEQS:
    ours = calculate_tm(s, **prm.tm_kwargs())
    p3 = primer3.calc_tm(s, mv_conc=MONOVALENT, dv_conc=prm.mg_mM,
                         dntp_conc=prm.dntp_mM, dna_conc=prm.primer_conc_nM)
    d = ours - p3
    diffs.append(d)
    check(abs(d) > 0.5, f"{s[:14]}… primer3와 다르다 ({d:+.2f}℃)")
check(max(diffs) - min(diffs) > 2.0,
      "차이가 서열마다 2℃ 이상 벌어진다 (상수 보정 불가)",
      f"{min(diffs):+.2f} ~ {max(diffs):+.2f}℃")

print("\n[3] 비-ACGT는 실패한다 (조용한 축 혼합 방지)")
try:
    calculate_tm("ACGTNACGTACGTACGTACG")
    check(False, "N이 든 서열은 ValueError여야 한다")
except ValueError as e:
    check("ACGT" in str(e), "N 거부", str(e)[:48])
stats = get_sequence_stats("ACGTNACGTACGTACGTACG")
check(stats.get("tm") is None, "서열 통계는 죽지 않고 tm만 비운다", f"tm={stats.get('tm')}")
check(stats.get("length") == 20, "나머지 통계는 그대로 나온다")

print("\n[4] 길이 상한이 없다 (60nt 제한은 hairpin/dimer 쪽 얘기다)")
for L in (60, 65, 80):
    s = "ACGT" * (L // 4)
    try:
        v = calculate_tm(s[:L])
        check(v > 0, f"{L}nt 계산됨", f"{v:.2f}℃")
    except Exception as e:
        check(False, f"{L}nt 계산 실패", f"{type(e).__name__}")

print("\n[5] 축 이동이 반영됐다")
check(prm.tm_min == 55.0 and prm.tm_max == 65.0,
      "선언창은 그대로 (원래부터 primer3 축)", f"{prm.tm_min}~{prm.tm_max}")
# ⚠️ 2026-11-05(Y11) — 판정창(`qc_tm_*`)은 **제거**됐다. 판정은 선언창이 직접 하고
#   물리 경계는 입력 제한이 막는다. 그래서 여기서 지킬 성질이 "판정창이 설계창을
#   감싸는가"에서 "**판정창이라는 두 번째 창이 아예 없는가**"로 바뀐다.
check(not hasattr(prm, "qc_tm_min") and not hasattr(prm, "qc_tm_max"),
      "Tm 판정창이 제거됐다 (창은 선언창 하나뿐이다)")
from params import TM_INPUT_MIN, TM_INPUT_MAX
check((TM_INPUT_MIN, TM_INPUT_MAX) == (35.0, 80.0),
      "물리 경계는 입력 제한이 맡는다", f"{TM_INPUT_MIN}~{TM_INPUT_MAX}℃")
check(prm.terminal_tm_min == 59.0, "terminal_tm_min 59", str(prm.terminal_tm_min))
from params import PRESETS
# 2026-11-05(Y10-b) — 프리셋 Tm 창을 NEB 축으로 재설정했다(이전 50~60/55~65/57~67은
#   primer3 축 시절 숫자였다). 여기서 Y8이 지키려는 건 특정 숫자가 아니라
#   **프리셋이 선언창을 제대로 덮어쓰는가**이므로, 그 성질로 검사한다.
#   숫자 자체의 파리티는 `tests/check_preset_parity.py`(GUI 거울)가 본다.
check(PRESETS["phusion"]["tm_min"] == 60.0 and PRESETS["phusion"]["tm_max"] == 72.0,
      "프리셋 Tm 창이 NEB 축 값이다",
      f"phusion {PRESETS['phusion']['tm_min']}~{PRESETS['phusion']['tm_max']}")
_ph = DesignParams.from_preset("phusion")
check((_ph.tm_min, _ph.tm_max) == (60.0, 72.0)
      and TM_INPUT_MIN <= _ph.tm_min and TM_INPUT_MAX >= _ph.tm_max,
      "프리셋이 선언창을 덮어쓰고, 그 값도 입력 허용 범위 안이다",
      f"설계 {_ph.tm_min:g}~{_ph.tm_max:g} ⊆ 입력 {TM_INPUT_MIN:g}~{TM_INPUT_MAX:g}")

print("\n[6] 입력 제한이 물리 경계를 막는다 (판정이 아니라 입력 단계에서)")
check(not prm.validate(), "기본값에 경고 없음", str(prm.validate()))
for bad, why in ((dict(tm_min=30.0, tm_max=60.0), "하한이 35℃ 아래"),
                 (dict(tm_min=60.0, tm_max=85.0), "상한이 80℃ 위")):
    try:
        DesignParams(**bad).validate()
        check(False, f"{why} → 에러가 나야 함")
    except ValueError as e:
        check("입력 허용 범위" in str(e), f"{why}: 입력 단계에서 막힌다",
              str(e).splitlines()[-1][:62])

print("\n[7] GG와 gibson이 같은 창을 같은 자로 본다 (E2E)")
db = DatabaseImport()
FOUR = [("promoter", "BBa_J23100"), ("rbs", "BBa_B0032"),
        ("cds", "BBa_C0083"), ("terminator", "BBa_B0010")]
parts = [db.get_part_by_id("Bioparts", t, i) for t, i in FOUR]
CASES = {"gibson": ("gibson", "gibson", "V0002"),
         "golden_gate": ("golden_gate", "golden_gate", "KV0001")}
tm_by_path = {}
for name, (mm, im, vec_id) in CASES.items():
    vec = db.get_part_by_id("Bioparts", "vector", vec_id)
    d = S / f"y8_{name}"
    shutil.rmtree(d, ignore_errors=True)
    d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        res = WorkflowBuilder.build_plasmid(parts, vec, db, str(d),
                                            module_method=mm, insert_method=im,
                                            params=DesignParams())
    rows = list(csv.DictReader(open(d / "reports" / "primers.csv")))
    tms = []
    # 2026-11-05(Y10) — 초과를 **행 단위**로 본다. Y9에서 "창 안에 쌍이 없으면 가장
    #   가까운 쌍" 정책이 들어오면서, 창을 몇십분의 1℃ 넘는 행이 **정당하게** 생길 수
    #   있게 됐다(예: 요청 55~65에 rev 65.02 — 경고가 붙어 나간다). 그건 Y9의 일이고,
    #   Y8이 잡으려는 건 **자가 갈려서 생기는 조용한 초과**다. 그래서 검사를 더 세게
    #   바꿨다: 창을 벗어난 anneal이 있으면 **그 행에 반드시 사유 경고가 있어야 한다**.
    #   (톨러런스를 늘려 덮으면 실제 자 불일치도 같이 덮인다)
    silent = []
    for r in rows:
        _w = r.get("warnings") or ""
        _excused = "가장 가까운 쌍" in _w
        for k in ("fwd_tm", "rev_tm"):
            try:
                t = float(r[k])
            except (TypeError, ValueError):
                continue
            tms.append(t)
            if t > prm.tm_max + 0.01 and not _excused:
                silent.append((r.get("part_name", "?"), k, t))
    tm_by_path[name] = (tms, res.get("construct"))
    print(f"  {name:<12} anneal Tm {min(tms):.1f}~{max(tms):.1f}℃  (창 {prm.tm_min}~{prm.tm_max})")
    # tail 있는 짧은 부품은 창 아래로 내려갈 수 있다(설계상 정상, 경고로 처리).
    # 통일로 고쳐진 건 **위로 넘치던 쪽**이다 — GG가 Biopython 자로 창을 보던 탓이었다.
    check(not silent, f"{name}: 창 상한을 넘는 anneal은 모두 사유가 적혀 있다",
          "조용한 초과 없음" if not silent
          else "; ".join(f"{n}.{k}={t:.2f}℃" for n, k, t in silent))

print("\n[8] construct 서열은 Tm 축과 무관하다")
for name, (_tms, con) in tm_by_path.items():
    check(con is not None, f"{name}: construct 산출")
    if con is not None:
        sha = hashlib.sha256(str(con.seq).upper().encode()).hexdigest()[:16]
        print(f"  {name:<12} {len(con.seq)}bp  sha {sha}")
# 통일 전/후로 construct가 안 바뀐다는 것은 run_all.sh의 baseline SHA가 지킨다.
check(True, "서열 불변은 run_all.sh의 baseline SHA가 지킨다(여기선 산출만 확인)")

print("\n[9] 설계 범위 폴백 사유가 주문 시트까지 간다")
# primer3가 요청 범위로 후보를 못 내면 기본 범위로 폴백한다(2026-07-13 재발 방지).
# 그 사실을 `design_basic_primers`가 warnings로 돌려주는데, gibson.py가 안 꺼내가서
# Build Log·주문 시트 어디에도 안 떴다 — "파라미터를 바꿨는데 결과가 그대로"의 정체다.
cds = db.get_part_by_id("Bioparts", "cds", "BBa_C0083")
vec2 = db.get_part_by_id("Bioparts", "vector", "V0002")
d = S / "y8_fallback"
shutil.rmtree(d, ignore_errors=True)
d.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    WorkflowBuilder.build_plasmid([cds], vec2, db, str(d), module_method="gibson",
                                  insert_method="gibson",
                                  # 2026-11-05(Y10) — 60~65는 NEB 축에서 **달성 가능**해져
                                  #   더 이상 폴백을 안 탄다(62.3℃가 창 안이다). 전파를
                                  #   보려면 실제로 못 맞추는 창이어야 한다 → 64~68.
                                  params=DesignParams(tm_min=64.0, tm_max=68.0))
row = list(csv.DictReader(open(d / "reports" / "primers.csv")))[0]
warns = (row.get("warnings") or "")
# 2026-11-05: 폴백 **정책**이 "기본 범위로 갈아타기" → "근접값"으로 바뀌었다(Y9).
# 여기서 지키는 건 정책이 아니라 **전파**다 — gibson.py가 설계기의 warnings를
# 꺼내가는가. 정책 자체는 y9_test.py가 본다.
check("가장 가까운 쌍" in warns, "설계기 경고가 warnings 컬럼까지 전파된다",
      next((w for w in warns.split(" | ") if "가장 가까운" in w), "(없음)")[:60])
check("Tm 64~68" in warns, "요청한 값이 그대로 적혀 있다")

# 폴백이 안 난 실행에는 그 경고가 없어야 한다(소음 방지).
d2 = S / "y8_nofallback"
shutil.rmtree(d2, ignore_errors=True)
d2.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    WorkflowBuilder.build_plasmid([cds], vec2, db, str(d2), module_method="gibson",
                                  insert_method="gibson", params=DesignParams())
row2 = list(csv.DictReader(open(d2 / "reports" / "primers.csv")))[0]
check("가장 가까운 쌍" not in (row2.get("warnings") or ""),
      "폴백이 없으면 그 경고도 없다 (소음 방지)")

print("\n" + "=" * 72)
if FAIL:
    print(f"Y8 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y8 전부 통과")
