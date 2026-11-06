"""Y6 검증 — 양 끝 spacer가 주문서에 제대로 나오는가 (2026-11-04에 실제로 터진 3건).

## 무엇이 터졌나

**(A) vector homology spacer가 "긴 부품"으로 분류돼 올리고로 빠졌다.**
`VectorHomologySpacer.arm`과 `JunctionPlanner`의 `amin`이 **같은 변수**
(`params.gibson_overlap`)에서 나온다. 그래서 `_lay_out`의 `len >= amin` 비교가
**항상 참**이 되어 spacer가 독립 fragment가 되고, 그 다음 `should_use_oligo`에 걸려
`production="oligo"`로 바뀌었다. 결과(실측, Gibson only 1 CDS + V0001):

  · fragment 3개(올리고 2 + PCR 1) — 정석은 PCR 1개
  · 프라이머 시트에 **insert 맨몸만** 남고 벡터 homology가 어디에도 없음
  · 올리고 시트가 불필요하게 생성됨

tail 흡수로 가는 길도 막혀 있었다 — `tail_fits`가 `amax`를 더해서(45+30+15=90nt)
주문 상한 80nt를 넘는다고 판정했다. 그런데 체인 맨 끝은 **바깥쪽 arm이 없다**
(접합부는 fragment 사이에만 생긴다) → 실제 비용은 30+15=45nt다.

**(B) 멀티모듈이 `spacer_mode`를 아예 보지 않았다.** GUI·XML 기본값은
`vector_homology`("벡터 맞춤(권장)")인데 `build_multimodule_plasmid`는 양 끝 경계를
`SpacerSelector`로만 골라, 벡터와 무관한 라이브러리 spacer(S0002…)가 양 끝에 붙었다.
→ 흉터 114bp(57×2)가 남고 삽입 PCR이 따로 필요했다.

**(C) Tm 판정창이 설계창을 겸해서** tail 없는 끝단 프라이머만 Tm 미달이 red로 올라갔다
(tail 있는 쪽은 같은 미달이 yellow). 판정창을 `qc_tm_*`로 다시 분리했다 —
단 **판정창이 설계창을 포함**하게 두어 2026-07-13 '전부 red'의 재발 경로를 막는다.

## 수용 기준

  1. 길이가 `amin`과 같은 체인 맨 끝 spacer는 tail로 흡수된다(독립 fragment가 아니다)
  2. tail 상한을 넘는 긴 spacer는 흡수하지 않는다(길이 무시하면 설계가 죽는다 — y2 회귀)
  3. 부품 1개 + 양 끝 spacer → fragment **1개**가 허용된다(상대는 선형 벡터다)
  4. Gibson only E2E: 올리고 0개 + 프라이머가 벡터 homology를 들고 있다
  5. 멀티모듈 `vector_homology`: 양 끝은 벡터 유래, **연결 spacer는 라이브러리**
  6. 멀티모듈 `library`: 양 끝도 풀에서 (기존 동작 회귀)
  7. `validate_primer`의 Tm 판정은 `qc_tm_*`를 쓴다(설계창 `tm_*`가 아니다)
  8. 판정창이 설계창을 못 덮으면 `validate()`가 경고한다

실행: python tools_dev/params_cp_tests/y6_test.py
"""
import io
import contextlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from pydna.dseqrecord import Dseqrecord

from basic.primer_design import ValidatePrimer
from params import DesignParams
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from workflow.junction_plan import JunctionPlanner

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/"
         "2f76d81c-815a-4c70-bfd6-00e46c3dcbd8/scratchpad")
FAIL = []


def check(cond, msg, extra=""):
    print(("  ✅ " if cond else "  ❌ ") + msg + (f"  ({extra})" if extra else ""))
    if not cond:
        FAIL.append(msg)


def seeded(n, seed):
    """재현 가능한 의사난수 서열 — 반복이 없어야 접합부 교차 검사를 통과한다."""
    import random
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def part(pid, ptype, seq):
    r = Dseqrecord(seq, name=pid)
    r.annotations["part_id"] = pid
    r.annotations["part_type"] = ptype
    r.annotations["molecule_type"] = "DNA"
    return r


prm = DesignParams()
AMIN = prm.gibson_overlap

print(f"[1] 길이 == amin({AMIN}bp)인 맨 끝 spacer는 tail로 흡수된다")
parts = [
    part("L_hom", "spacer", seeded(AMIN, 11)),
    part("CDS", "cds", seeded(900, 12)),
    part("R_hom", "spacer", seeded(AMIN, 13)),
]
plan = JunctionPlanner(params=prm).plan(parts)
names = [f.name for f in plan.fragments]
check(len(plan.fragments) == 1, "fragment 1개로 합쳐짐", f"{len(plan.fragments)}개: {names}")
check(all(f.production == "pcr" for f in plan.fragments),
      "올리고로 빠진 fragment가 없음",
      str([f.production for f in plan.fragments]))
f0 = plan.fragments[0]
check(f0.prefix == str(parts[0].seq).upper(), "L spacer가 forward tail(prefix)로 들어감")
check(f0.suffix == str(parts[2].seq).upper(), "R spacer가 reverse tail(suffix)로 들어감")
check(plan.module_seq == "".join(str(p.seq).upper() for p in parts),
      "module 서열은 그대로(흡수가 산물을 바꾸지 않는다)")

print("\n[2] tail 상한을 넘는 긴 spacer는 흡수하지 않는다 (y2 회귀)")
# part_type이 'spacer'여도 3000bp면 tail에도 올리고에도 못 들어간다 → 긴 부품으로 봐야 한다.
big = [
    part("VEC", "spacer", seeded(3000, 21)),
    part("P", "promoter", seeded(250, 22)),
    part("C", "cds", seeded(800, 23)),
]
try:
    bplan = JunctionPlanner(params=prm).plan(big)
    got = [(f.name, f.production) for f in bplan.fragments]
    check(any(n == "VEC" for n, _ in got), "3000bp spacer는 독립 fragment로 남음", str(got))
except Exception as e:
    check(False, "3000bp spacer로 설계가 죽지 않음", f"{type(e).__name__}: {e}")

print("\n[3] 맨 끝 흡수는 amax 몫을 세지 않는다")
# amax + sp + anneal최소 = 45+30+18 = 93 > 80(주문 상한)이지만, 맨 끝은 arm이 없으므로
# 실제 비용은 30+18 = 48 ≤ 80이다. [1]이 흡수에 성공한 것이 곧 이 증거다.
check(prm.gibson_overlap_max + AMIN + prm.primer_len_min > 80,
      "옛 식이라면 상한을 넘었을 조건임을 확인",
      f"{prm.gibson_overlap_max}+{AMIN}+{prm.primer_len_min}"
      f"={prm.gibson_overlap_max + AMIN + prm.primer_len_min}nt > 80nt")
check(len(plan.fragments) == 1, "그래도 [1]에서 흡수됐다 = arm 몫을 안 셌다")

print("\n[4] Gibson only E2E — 올리고 0개 + 프라이머에 벡터 homology")
db = DatabaseImport()
cds = db.get_part_by_id("Bioparts", "cds", "BBa_C0020")
vec = db.get_part_by_id("Bioparts", "vector", "V0001")
d = S / "y6_gib"
shutil.rmtree(d, ignore_errors=True)
d.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    res = WorkflowBuilder.build_plasmid(
        [cds], vec, db, str(d), module_method="gibson", insert_method="gibson",
        params=DesignParams(),
    )
con = res.get("construct")
check(con is not None, "construct 산출")
if con is not None:
    check(len(con.seq) == len(vec.seq) + len(cds.seq),
          "흉터 없이 벡터+insert 길이",
          f"{len(con.seq)}bp = {len(vec.seq)}+{len(cds.seq)}")
import csv as _csv
rows = list(_csv.DictReader(open(d / "reports" / "primers.csv")))
check(len(rows) == 1, "주문할 fragment 1종", f"{len(rows)}종")
check(all((r.get("production") or "") != "oligo" for r in rows),
      "올리고로 빠진 fragment 없음")
if rows:
    fwd = (rows[0].get("forward_primer") or "").upper()
    rev = (rows[0].get("reverse_primer") or "").upper()
    vs = str(vec.seq).upper()
    # fwd 5' tail = 선형 벡터의 3' 끝, rev 5' tail = 벡터 5' 끝의 역상보.
    # 선형화 자리를 모르므로 "앞 arm만큼이 벡터(원형) 안에 있는가"로 본다.
    arm = prm.gibson_overlap
    check(fwd[:arm] in vs + vs[:arm], "fwd 5' tail이 벡터 서열", fwd[:arm])
    from Bio.Seq import Seq
    rtail = str(Seq(rev[:arm]).reverse_complement())
    check(rtail in vs + vs[:arm], "rev 5' tail(RC)이 벡터 서열", rtail)
    check(len(fwd) <= 80 and len(rev) <= 80,
          "둘 다 주문 상한 80nt 이내", f"fwd {len(fwd)}nt / rev {len(rev)}nt")

print("\n[5] 멀티모듈 vector_homology — 양 끝은 벡터, 연결은 라이브러리")
gA = [db.get_part_by_id("Bioparts", t, i) for t, i in
      [("promoter", "BBa_J23100"), ("cds", "BBa_C0012")]]
gB = [db.get_part_by_id("Bioparts", t, i) for t, i in
      [("promoter", "BBa_J23101"), ("terminator", "BBa_B0010")]]
fin = db.get_part_by_id("Bioparts", "vector", "V0002")
d2 = S / "y6_mm_vh"
shutil.rmtree(d2, ignore_errors=True)
d2.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    r2 = WorkflowBuilder.build_multimodule_plasmid(
        [gA, gB], fin, db, str(d2), insert_method="gibson",
        params=DesignParams(spacer_mode="vector_homology"),
    )
b2 = r2["boundaries"]
check(len(b2) == 3, "경계 3개(module 2개)", str(b2))
check(b2[0].startswith(fin.name[:11]) and b2[-1].startswith(fin.name[:11]),
      "양 끝 경계가 벡터에서 만들어졌다", f"{b2[0]} / {b2[-1]}")
check(not b2[1].startswith(fin.name[:11]),
      "연결 경계는 라이브러리 spacer다(module↔module은 벡터 접점이 아니다)", b2[1])
check(r2.get("construct") is not None, "construct 산출",
      f"{len(r2['construct'].seq)}bp" if r2.get("construct") else "")

print("\n[6] 멀티모듈 library — 양 끝도 풀에서 (기존 동작)")
d3 = S / "y6_mm_lib"
shutil.rmtree(d3, ignore_errors=True)
d3.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    r3 = WorkflowBuilder.build_multimodule_plasmid(
        [gA, gB], fin, db, str(d3), insert_method="gibson",
        params=DesignParams(spacer_mode="library"),
    )
b3 = r3["boundaries"]
check(all(not x.startswith(fin.name[:11]) for x in b3),
      "경계 전부가 라이브러리 spacer", str(b3))
check(len(set(b3)) == len(b3), "경계끼리 중복 없음", str(b3))

print("\n[7] Tm 판정은 qc_tm_* 를 쓴다 (설계창 tm_* 가 아니다)")
# 설계창을 [70,75]로 좁혀도, 판정창 [50,80] 안에 드는 프라이머는 Tm으로 떨어지지 않아야 한다.
pr = "ATGGCGGATAAAGAACTTAAATT"
narrow = DesignParams(tm_min=70.0, tm_max=75.0, qc_tm_min=50.0, qc_tm_max=80.0)
v = ValidatePrimer().validate_primer(pr, params=narrow)
# 필터는 **"Tm" 부분문자열이 아니라 실제 지적 문구**로 잡는다 — 다른 축의 메시지가
# 설명으로 "Tm"을 언급할 수 있다(GC 경고문의 "(Tm은 별도 판정)"에 2026-11-05 실제로 걸렸다).
tm_msgs = [m for m in v["violations"] + v["warnings"] if "Tm이" in m]
check(not tm_msgs, "설계창 밖이어도 판정창 안이면 Tm 지적 없음",
      f"Tm {v['tm']:.1f}℃ / 설계창 70~75 / 판정창 50~80" + (f" → {tm_msgs}" if tm_msgs else ""))
tight = DesignParams(tm_min=55.0, tm_max=65.0, qc_tm_min=60.0, qc_tm_max=80.0)
v2 = ValidatePrimer().validate_primer(pr, params=tight)
check(any("Tm이" in m for m in v2["violations"] + v2["warnings"]),
      "판정창을 올리면 같은 프라이머가 Tm으로 지적된다", f"Tm {v2['tm']:.1f}℃ / 판정창 60~80")

print("\n[8] 판정창이 설계창을 못 덮으면 경고")
warns = DesignParams(tm_min=55.0, tm_max=65.0, qc_tm_min=60.0, qc_tm_max=80.0).validate()
check(any("판정창" in w for w in warns), "포함 관계 위반을 경고한다",
      next((w for w in warns if "판정창" in w), "(경고 없음)"))
clean = DesignParams().validate()
check(not any("판정창" in w for w in clean), "기본값(50~80 ⊃ 55~65)은 경고 없음")

print("=" * 72)
print(f"Y6 {'통과' if not FAIL else f'실패 {len(FAIL)}건'}")
sys.exit(0 if not FAIL else 1)
