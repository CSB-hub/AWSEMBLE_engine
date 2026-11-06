"""CP6: 파라미터가 다르면 artifact id가 갈리고 서로 덮어쓰지 않는가."""
import sys, shutil, json, io, contextlib
from pathlib import Path
sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parents[2] / "src"))
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from params import DesignParams

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
db = DatabaseImport()
SPEC=[("promoter","BBa_J23100"),("rbs","BBa_B0032"),("cds","BBa_C0012"),("terminator","BBa_B0010")]
parts=[db.get_part_by_id("Bioparts",t,i) for t,i in SPEC]
vector=db.get_part_by_id("Bioparts","vector","pSB1C3")
ok=True
def chk(l,c,e=""):
    global ok; ok = ok and c
    print(f"  [{'OK ' if c else 'FAIL'}] {l}" + (f" — {e}" if e else ""))

def run(proj, prm):
    d = S/"cp6"/proj; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
    with contextlib.redirect_stdout(io.StringIO()):
        res = WorkflowBuilder.build_plasmid(parts, vector, db, str(d),
                                            module_method="golden_gate",
                                            insert_method="gibson", params=prm)
    mf = json.loads((d/"reports"/"manifest.json").read_text(encoding="utf-8"))
    entries = mf if isinstance(mf, list) else mf.get("artifacts", list(mf.values())[0] if mf else [])
    return d, res, entries

print("1) 같은 입력 2회 → 같은 id (결정성)")
d1, r1, e1 = run("a1", DesignParams())
d2, r2, e2 = run("a2", DesignParams())
ids1 = [x["id"] for x in e1]; ids2 = [x["id"] for x in e2]
chk("id 동일", ids1 == ids2, f"{ids1[:2]}")

print("\n2) 파라미터만 다르면 → 다른 id")
d3, r3, e3 = run("b1", DesignParams(tm_min=67.0, tm_max=77.0))
ids3 = [x["id"] for x in e3]
chk("module id가 갈림", ids1[0] != ids3[0], f"{ids1[0]} vs {ids3[0]}")
chk("construct id가 갈림", ids1[-1] != ids3[-1], f"{ids1[-1]} vs {ids3[-1]}")

print("\n3) 같은 프로젝트에 두 파라미터 → 덮어쓰지 않고 둘 다 남는가")
d = S/"cp6"/"shared"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
with contextlib.redirect_stdout(io.StringIO()):
    WorkflowBuilder.build_plasmid(parts, vector, db, str(d), params=DesignParams())
    WorkflowBuilder.build_plasmid(parts, vector, db, str(d), params=DesignParams(tm_min=67.0, tm_max=77.0))
mods = sorted((d/"modules"/"assembled").glob("*.gb"))
cons = sorted((d/"constructs").glob("*.gb"))
chk("module 파일 2개 공존", len(mods) == 2, f"{[m.name for m in mods]}")
chk("construct 파일 2개 공존", len(cons) == 2, f"{len(cons)}개")

print("\n4) manifest에 설계 파라미터 기록")
has = [x for x in e1 if isinstance(x.get("params"), dict) and "design_params" in x["params"]]
chk("design_params 기록됨", len(has) > 0, f"{len(has)}/{len(e1)} 엔트리")
if has:
    dp = has[0]["params"]["design_params"]
    # 2026-11-05 Y10 — Tm 축을 NEB Tm Calculator 재현으로 바꾸면서 완충액 기본값이 바뀌었다.
    #   예전: na_mM=50 (Biopython 관용 기본값)
    #   지금: NEB Q5 조성 → na_mM=0, k_mM=50, tris_mM=10 (1가 염은 K/Tris가 담당)
    # 여기서 보는 건 "완충액 축이 manifest로 복원되는가"이므로 기준값만 갈아끼운다.
    chk("Tm 조건 복원 가능",
        dp.get("na_mM") == 0.0 and dp.get("k_mM") == 50.0
        and (dp.get("tm_min"), dp.get("tm_max")) == (55.0, 65.0),
        f"tm=[{dp.get('tm_min')},{dp.get('tm_max')}], na={dp.get('na_mM')}, "
        f"k={dp.get('k_mM')}, enzyme={dp.get('gg_enzyme')}")

print("\n=== CP6:", "통과" if ok else "실패", "===")
sys.exit(0 if ok else 1)
