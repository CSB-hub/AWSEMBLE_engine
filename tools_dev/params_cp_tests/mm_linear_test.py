"""CP7의 record.name 변경이 multimodule(linear) 경로를 깼는지 확인."""
import sys, io, contextlib, shutil
from pathlib import Path
sys.path.insert(0, "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src")
from utils.io import DatabaseImport
from workflow.builder import WorkflowBuilder
from params import DesignParams
S=Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
db=DatabaseImport()
fin=db.get_part_by_id("Bioparts","vector","pSB1A2")
gA=[db.get_part_by_id("Bioparts",t,i) for t,i in [("promoter","BBa_J23100"),("rbs","BBa_B0032")]]
gB=[db.get_part_by_id("Bioparts",t,i) for t,i in [("terminator","BBa_B0010")]]
d=S/"mm_lin"; shutil.rmtree(d, ignore_errors=True); d.mkdir(parents=True)
buf=io.StringIO()
try:
    with contextlib.redirect_stdout(buf):
        res=WorkflowBuilder.build_multimodule_plasmid([gA,gB], fin, db, str(d),
             insert_method="gibson", params=DesignParams())
    print(f"✅ multimodule(linear) 정상 — construct {len(res['construct'].seq)}bp")
    print(f"   boundaries={res['boundaries']}")
except Exception as e:
    print(f"❌ 깨짐: {type(e).__name__}: {e}")
    print("--- 로그 마지막 6줄 ---")
    print("\n".join(buf.getvalue().splitlines()[-6:]))
    sys.exit(1)
sys.exit(0)
