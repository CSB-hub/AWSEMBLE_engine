"""Y2 검증 — yeast HR 고유 3종 (교차 차단 / 벡터 적격성 / in vivo 전제).

수용 기준
  1. 반복 부품(비인접 fragment 교차)을 **차단**한다 — Gibson은 경고만
  2. 차단 메시지가 원인과 해결 방법을 알려준다
  3. padding으로 접합부를 다르게 만들면 통과한다 (해결 경로가 실제로 동작)
  4. `check_yeast_vector`가 feature 기반으로 ok/partial/unknown을 가린다
  5. 메타데이터 없는 벡터는 **unknown(차단 아님)** — "없다"가 아니라 "확인 불가"
  6. `validate()`가 in vivo 전제와 벡터 경고를 warnings에 싣는다
  7. 정상 설계는 서열 오라클을 만족하고 원형 산물이 **정규화 후 유일**하다
  8. `insert_method="yeast_hr"`이 vector_insertion 경로에서 동작한다
실행: python tools_dev/params_cp_tests/y2_test.py
"""
import os
import random
import shutil
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "src"))

from Bio.SeqFeature import FeatureLocation, SeqFeature
from pydna.dseqrecord import Dseqrecord

from assembly.factory import AssemblyFactory
from assembly.yeast_hr import YeastHRCloning, check_yeast_vector
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


def seeded(n, seed):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


_COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(_COMP)[::-1]


def same_circle(a, b):
    """원형 분자 동일성 — 회전·역상보 정규화."""
    return len(a) == len(b) and ((b in a + a) or (rc(b) in a + a))


print("=" * 72)
print("Y2 — yeast HR 고유 3종")
print("=" * 72)

prm = DesignParams()
prm.validate()

# 같은 부품(REP)을 두 자리에 쓴 설계 — 비인접 fragment가 서열을 공유하게 된다.
REP = seeded(320, 81)
repeat_parts = [
    part("A", "promoter", seeded(260, 82)),
    part("REP1", "cds", REP),
    part("B", "cds", seeded(260, 83)),
    part("REP2", "cds", REP),
]

print("\n[1] 반복 부품 — arm이 같아져서 양쪽 다 차단")
# 2026-11-06 정책 변경: 판정이 `_check_plan`(방법별 등급)에서 `JunctionPlanner`의
#   **arm 기준 + 재설계 래더**로 옮겨졌다. 같은 부품을 두 자리에 쓰면 두 접합부의
#   arm이 **동일한 서열**이 되므로 어떤 길이로도 못 피한다 → Gibson도 차단된다.
#   (Gibson에서 arm은 실제로 어닐링하는 구간이다 — 경고로 둘 수 있는 축이 아니다.)
#   방법별 차이는 "차단/경고"가 아니라 **비교 범위**로 남았다(`CROSSTALK_SCOPE`).
gib = AssemblyFactory.create_assembler("gibson", params=prm)
yeast = AssemblyFactory.create_assembler("yeast_hr", params=prm)
blocks = {}
for method, asm in (("gibson", gib), ("yeast_hr", yeast)):
    try:
        asm.prepare_parts(repeat_parts)
        blocks[method] = ""
    except RuntimeError as e:
        blocks[method] = str(e)
    print(f"  {method}: {'차단됨' if blocks[method] else '차단 안 됨(문제)'}")
check(blocks["yeast_hr"], "yeast_hr은 반복 부품 설계를 차단한다")
check(blocks["gibson"],
      "Gibson도 차단한다 — arm이 동일하면 어떤 길이로도 못 피한다")
check(gib.CROSSTALK_SCOPE == "arm" and yeast.CROSSTALK_SCOPE == "arm+body",
      "방법별 차이는 등급이 아니라 비교 범위로 남았다")

print("\n[2] 차단 메시지 품질")
msg = blocks["yeast_hr"]
if msg:
    print("  " + "\n  ".join(msg.splitlines()[:3]))
    check("homology arm also matches" in msg, "무엇이 문제인지(arm이 다른 자리에 붙음) 설명")
    check("no junction length between" in msg, "왜 차단인지(래더를 다 돌았다) 설명")
    check("Replace the part at this position" in msg, "해결 방법(부품 교체)을 알려줌")

print("\n[3] 반복 부품을 교체하면 통과 (padding 경로는 삭제됨)")
y_pad = AssemblyFactory.create_assembler("yeast_hr", params=DesignParams())
# 교차의 원인은 부품 본체가 같은 것이므로, 해결은 부품 교체뿐이다.
fixed_parts = [
    part("A", "promoter", seeded(260, 82)),
    part("REP1", "cds", REP),
    part("B", "cds", seeded(260, 83)),
    part("REP3", "cds", seeded(320, 84)),   # 다른 부품으로 교체
]
try:
    f_design, f_plan = y_pad.prepare_parts(fixed_parts)
    ok = True
    print(f"  부품 교체 후: 설계 성공, 교차 최대 {f_plan.max_crosstalk}bp, "
          f"arm {f_plan.arm}bp")
except RuntimeError as e:
    ok = False
    print(f"  실패: {str(e).splitlines()[0]}")
check(ok, "반복을 없애면 yeast 설계가 통과한다")

print("\n[4] check_yeast_vector — feature 기반 판정")


def vec_with_features(name, labels):
    v = Dseqrecord(seeded(4000, 91), name=name, circular=True)
    v.annotations["molecule_type"] = "DNA"
    for i, lab in enumerate(labels):
        v.features.append(SeqFeature(
            FeatureLocation(100 * (i + 1), 100 * (i + 1) + 60),
            type="misc_feature", qualifiers={"label": [lab]},
        ))
    return v


cases = [
    ("ori+marker", vec_with_features("pRS316_like", ["CEN6/ARSH4", "URA3"]), "ok"),
    ("ori만", vec_with_features("v_ori", ["ARSH4"]), "partial"),
    ("marker만", vec_with_features("v_mk", ["LEU2"]), "partial"),
    ("2micron+marker", vec_with_features("YEp_like", ["2micron ori", "HIS3"]), "ok"),
    ("무관 feature", vec_with_features("pUC19_like", ["lacZ", "AmpR"]), "unknown"),
    ("feature 없음", Dseqrecord(seeded(4000, 92), name="V0001", circular=True), "unknown"),
]
for label, vec, expect in cases:
    r = check_yeast_vector(vec)
    print(f"  {label:<16} → {r['verdict']:<8} ori={r['ori']} marker={r['marker']}")
    check(r["verdict"] == expect, f"{label} → {expect}")

print("\n[5] 메타데이터 없음은 unknown (차단 아님)")
bare = Dseqrecord(seeded(4000, 92), name="V0001", circular=True)
r = check_yeast_vector(bare)
check(r["verdict"] == "unknown", "verdict=unknown")
check("cannot be confirmed" in r["warning"], "'확인 불가'로 표현(없다고 단정 안 함)")
check("does not\n            mean they are absent" in r["warning"].replace("  ", " ") or "not mean they are absent" in r["warning"], "오해를 막는 문구 포함")
# 이름 단서만 있는 경우 → partial (보조 증거)
hinted = Dseqrecord(seeded(4000, 93), name="pRS315", circular=True)
rh = check_yeast_vector(hinted)
print(f"  이름 단서만('pRS315') → {rh['verdict']}, hints={rh['hints']}")
check(rh["verdict"] == "partial", "이름 단서는 partial로 올려준다")

print("\n[6] validate()가 in vivo 전제와 벡터 경고를 싣는다")
construct = Dseqrecord(seeded(5000, 94), name="con_test", circular=True)
construct.annotations["molecule_type"] = "DNA"
v_novec = yeast.validate(construct)
check(any(YeastHRCloning.IN_VIVO_NOTE in w for w in v_novec["warnings"]),
      "in vivo 전제가 warnings에 있음")
v_withvec = yeast.validate(construct, vector=bare)
check("yeast_vector" in v_withvec["details"], "details에 벡터 판정 내역")
check(any("cannot be confirmed" in w for w in v_withvec["warnings"]),
      "벡터 경고가 warnings에 추가됨")
good_vec = vec_with_features("pRS316_like", ["CEN6/ARSH4", "URA3"])
v_good = yeast.validate(construct, vector=good_vec)
check(v_good["details"]["yeast_vector"]["verdict"] == "ok", "적격 벡터는 ok")
check(sum(1 for w in v_good["warnings"] if "vector" in w.lower()) == 0,
      "적격 벡터면 벡터 경고 없음")

print("\n[7] 서열 오라클 + 원형 산물 정규화 후 유일")
# 난수 서열은 프라이머 설계가 실패할 수 있다(실측 모듈 성공률 83%). 이 절은 **오라클과
# 유일성**을 보는 것이므로, 설계가 성립하는 서열 조합을 먼저 찾아서 그 위에서 검사한다.
# (설계 성공률 자체는 tools_dev/junction_scan.py가 측정한다)


def first_working(make_parts, tries=30):
    """설계가 성립하는 첫 입력을 찾아 (parts, design, plan)을 돌려준다."""
    for s in range(tries):
        ps = make_parts(s)
        try:
            d, p = yeast.prepare_parts(ps)
        except RuntimeError:
            continue
        return ps, d, p, s
    return None, None, None, None


clean_parts, design, plan, used_seed = first_working(lambda s: [
    part("P", "promoter", seeded(250, 1100 + s)),
    part("R", "rbs", seeded(16, 1200 + s)),
    part("C", "cds", seeded(800, 1300 + s)),
    part("T", "terminator", seeded(150, 1400 + s)),
])
check(plan is not None, "설계가 성립하는 모듈을 찾음")
if plan is None:
    print("  설계 가능한 조합을 찾지 못해 [7]~[8] 생략")
    print("=" * 72)
    print(f"Y2 실패 {len(FAIL)}건")
    sys.exit(1)
print(f"  seed #{used_seed}: fragment {len(plan.fragments)}개, arm {plan.arm}bp, "
      f"arm {plan.arm}bp, "
      f"module {len(plan.module_seq)}bp")
lin = yeast.assemble(plan.records(), assembly_type="linear")
check(any(str(r.seq).upper() == plan.module_seq for r in lin), "linear 오라클 일치")

# 원형: module 양 끝이 이어지도록 벡터 역할 fragment를 추가한다.
circ_parts, cdesign, cplan, cseed = first_working(lambda s: [
    part("VEC", "spacer", seeded(3000, 1500 + s))] + clean_parts)
check(cplan is not None, "벡터 포함 원형 설계가 성립하는 조합을 찾음")
if cplan is None:
    print("=" * 72)
    print(f"Y2 실패 {len(FAIL)}건")
    sys.exit(1)
# 마지막 fragment의 오른쪽 끝과 첫 fragment의 왼쪽 끝을 잇는 arm을 수동 부착
first, last = cplan.fragments[0], cplan.fragments[-1]
closing = first.body[:cplan.arm]
last.right_arm = closing
circ = yeast.assemble(cplan.records(), assembly_type="circular")
seqs = [str(r.seq).upper() for r in circ]
uniq = []
for s in seqs:
    if not any(same_circle(s, u) for u in uniq):
        uniq.append(s)
print(f"  circular 산물 {len(seqs)}개 → 정규화 후 {len(uniq)}개 "
      f"(길이 {sorted({len(s) for s in seqs})})")
check(len(uniq) == 1, "정규화 후 산물 유일")
# ⚠️ 기대 서열은 `module_seq`**그 자체**다. 닫는 arm(`closing`)은 첫 fragment의 머리
#    서열이라 원형에서 **공유되어 한 번만** 나타난다 — 더해서 세면 arm만큼 중복 계산된다.
#    이게 homology 조립이 흉터를 남기지 않는다는 성질의 직접적인 확인이다.
expect_circ = cplan.module_seq
check(same_circle(uniq[0], expect_circ),
      f"원형 산물이 기대 서열과 일치 ({len(expect_circ)}bp, 흉터 0)")
check(len(uniq[0]) == len(cplan.module_seq),
      f"닫는 arm {cplan.arm}bp가 중복되지 않음")

print("\n[8] insert_method='yeast_hr' vector_insertion 경로")
from workflow.biopart_assemble import BiopartWorkflow

module_rec = Dseqrecord(plan.module_seq, name="mod_y2")
module_rec.annotations["molecule_type"] = "DNA"
module_rec.annotations["insert_whole"] = "1"
# 벡터: module 양 끝 arm을 들고 있는 선형 벡터(= vector_homology 모드와 같은 형태)
vec_lin = Dseqrecord(
    plan.module_seq[-plan.arm:] + seeded(3200, 106) + plan.module_seq[:plan.arm],
    name="vec_y2",
)
vec_lin.annotations["molecule_type"] = "DNA"
tmp = tempfile.mkdtemp(prefix="y2_")
try:
    wf = BiopartWorkflow(project_dir=tmp, parts=[module_rec])
    out = wf.vector_insertion([module_rec, vec_lin], method="yeast_hr",
                              assembly_type="circular", design_params=prm)
    check(out is not None, "vector_insertion(method='yeast_hr') 성공")
    if out:
        got = str(out["record"].seq).upper()
        print(f"  construct {len(got)}bp (module {len(plan.module_seq)}bp + "
              f"vector backbone)")
        check(out["record"].circular, "원형 construct")
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print("\n" + "=" * 72)
if FAIL:
    print(f"Y2 실패 {len(FAIL)}건")
    for m in FAIL:
        print(f"  - {m}")
    sys.exit(1)
print("Y2 전부 통과")
