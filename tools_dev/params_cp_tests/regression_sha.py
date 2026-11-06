"""오늘 수정 전(baseline) / 후(fixed) 코드로 같은 파이프라인을 돌려 서열 SHA 비교."""
import sys, shutil, hashlib, json
from pathlib import Path

WHICH = sys.argv[1]
SRC = {"baseline": "/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad/baseline/src",
       "fixed": "/home/csb/Desktop/BioCAD/Awsemble/awsemble_dev/BioCAD_private/awsemble/src"}[WHICH]
sys.path.insert(0, SRC)
from utils.io import DatabaseImport, DNAFileReader, DNAFileWriter
from workflow.builder import WorkflowBuilder

S = Path("/tmp/claude-1000/-home-csb-Desktop-BioCAD/171b8caf-fb9b-4b74-bfce-a9ecd5beada0/scratchpad")
work = S/f"cmp_{WHICH}"; shutil.rmtree(work, ignore_errors=True); work.mkdir(parents=True)
db = DatabaseImport()
SPEC=[("promoter","BBa_J23100"),("rbs","BBa_B0032"),("cds","BBa_C0012"),("terminator","BBa_B0010")]
files=[]
for i,(t,pid) in enumerate(SPEC,1):
    rec=db.get_part_by_id("Bioparts",t,pid); fp=work/f"{i}_{pid}.genbank"
    DNAFileWriter.write_genbank(rec,fp); files.append(fp)
parts=[DNAFileReader.read_genbank(str(f))[0] for f in files]
vector=db.get_part_by_id("Bioparts","vector","pSB1C3")
from params import DesignParams
# ⚠️ spacer_mode="library"로 고정한다. 2026-10-20에 기본값이 "vector_homology"로
#    바뀌었고 그 모드는 **최종 서열이 달라진다**(중립 spacer 흉터가 사라지므로).
#    이 테스트들은 기존 방식의 **하위호환 보증선**이므로 모드를 명시해 고정한다.
#    새 모드의 기준은 vector_homology_test.py가 따로 든다.
# ⚠️ gibson_overlap=20도 같은 이유로 고정한다. 2026-10-24에 기본값이 30으로 바뀌었고
#    이 값은 벡터 삽입 homology로 **최종 서열에 들어간다**. 고정하지 않으면 baseline SHA가
#    기본값을 바꿀 때마다 이동해 "코드가 서열을 바꿨나"를 더 이상 묻지 못하게 된다.
#    새 기본값(30)의 커버리지는 g1~g3·y1~y2가 DesignParams() 기본으로 든다.
res=WorkflowBuilder.build_plasmid(parts,vector,db,str(work/"proj"),
                                  module_method="golden_gate",insert_method="gibson",
                                  params=DesignParams(spacer_mode="library",
                                                      gibson_overlap=20))
def sha(r): return hashlib.sha256(str(r.seq).upper().encode()).hexdigest()[:16]
mod = res["results"][0]["record"]; con = res["construct"]
out = {"which": WHICH,
       "l_spacer": res["l_spacer"], "r_spacer": res["r_spacer"],
       "module_len": len(mod.seq), "module_sha": sha(mod),
       "construct_len": len(con.seq), "construct_sha": sha(con),
       "module_feature_types": sorted({f.type for f in mod.features}),
       "construct_feature_types": sorted({f.type for f in con.features})}
print("RESULT " + json.dumps(out, ensure_ascii=False))
