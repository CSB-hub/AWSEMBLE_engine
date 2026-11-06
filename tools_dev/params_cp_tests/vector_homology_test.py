"""spacer_mode="vector_homology" 기준선 — 기존 방식과 **무엇이 달라야 하는지**를 고정한다.

기존 테스트들(cp5/one_vector/regression_sha)은 library 모드로 고정해 하위호환을 지킨다.
이 테스트는 새 모드가 **의도한 차이**를 실제로 만들어내는지 본다:

  1. module 양 끝이 벡터 선형화 자리의 서열과 정확히 일치한다
  2. 그래서 삽입 PCR이 생략된다 (삽입 프라이머가 아예 없다)
  3. 중립 spacer 흉터가 최종 construct에서 사라진다 (library 대비 짧아진다)
  4. 부품은 전부 보존된다
"""
import contextlib
import io as _io
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "src"))

from pydna.dseqrecord import Dseqrecord          # noqa: E402
from params import DesignParams                   # noqa: E402
from utils.io import DatabaseImport               # noqa: E402
from workflow.builder import WorkflowBuilder      # noqa: E402
from workflow.vector_homology_spacer import (     # noqa: E402
    VectorHomologySpacer, module_carries_homology,
)

SPEC = [("promoter", "BBa_J23100"), ("rbs", "BBa_B0032"),
        ("cds", "BBa_C0012"), ("terminator", "BBa_B0010")]
VECTOR = "pSB1C3"

fails = []


def ok(cond, label, detail=""):
    print(f"  [{'OK ' if cond else 'FAIL'}] {label}" + (f" — {detail}" if detail else ""))
    if not cond:
        fails.append(label)


def run(mode, work):
    db = DatabaseImport()
    parts = [db.get_part_by_id("Bioparts", t, i) for t, i in SPEC]
    vector = db.get_part_by_id("Bioparts", "vector", VECTOR)
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    buf = _io.StringIO()
    with contextlib.redirect_stdout(buf):
        res = WorkflowBuilder.build_plasmid(
            parts, vector, db, str(work / "proj"),
            module_method="golden_gate", insert_method="gibson",
            params=DesignParams(spacer_mode=mode))
    return res, buf.getvalue(), vector, db


base = Path("/tmp/vh_test")
print("1) vector_homology 모드")
res_v, log_v, vector, db = run("vector_homology", base / "vh")
mod_v = res_v["results"][0]["record"]
con_v = res_v["construct"]

lin = VectorHomologySpacer(params=DesignParams()).linearize(vector)
# arm 길이는 params에서 읽는다 — 숫자를 박으면 기본값이 바뀌는 날 거짓 실패가 난다.
ok(module_carries_homology(mod_v, lin, DesignParams().gibson_overlap),
   "module 양 끝 = 벡터 선형화 자리 서열",
   f"{str(mod_v.seq)[:12]}… / …{str(mod_v.seq)[-12:]}")
ok("skipping the insertion PCR" in log_v, "삽입 PCR이 생략됨")

vp = work_csv = base / "vh/proj/reports/vector_primers.csv"
if vp.exists():
    body = vp.read_text(encoding="utf-8")
    ok("No insertion primers needed" in body, "리포트가 '프라이머 없음'을 명시")
else:
    ok(False, "vector_primers.csv 생성")

print("\n2) library 모드와 비교")
res_l, log_l, _, _ = run("library", base / "lib")
con_l = res_l["construct"]
ok("삽입 PCR 생략" not in log_l, "library는 기존대로 PCR 수행")
diff = len(con_l.seq) - len(con_v.seq)
ok(diff > 0, "흉터가 제거되어 construct가 짧아짐",
   f"library {len(con_l.seq)}bp → vector_homology {len(con_v.seq)}bp (−{diff}bp)")

print("\n3) 부품 보존")
seqs = {i: str(db.get_part_by_id("Bioparts", t, i).seq).upper() for t, i in SPEC}
cv = str(con_v.seq).upper()
for pid, sq in seqs.items():
    ok(sq in cv, f"{pid} 보존")

print()
if fails:
    print(f"=== vector_homology: 실패 {len(fails)}건 — {', '.join(fails)}")
    sys.exit(1)
print("=== vector_homology: 통과")
