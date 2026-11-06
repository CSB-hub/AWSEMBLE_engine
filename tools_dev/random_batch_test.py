
import sys, hashlib, io, contextlib, os
import random
from pathlib import Path
import pandas as pd
from datetime import datetime

# Add src to path
sys.path.append(os.path.abspath("Awsemble/awsemble_dev/BioCAD_private/awsemble/src"))

from utils.io import DatabaseImport, DNAFileReader, DNAValidator
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow as Workflow
from pydna.dseqrecord import Dseqrecord

PROJECT_DIR = Path("/tmp/random_test/test_results")
SCHEMA = "Bioparts"

def get_random_test_ids(db, count=10):
    query = """
        SELECT DISTINCT Test_ID 
        FROM test_part 
        WHERE Mapping = 'o' 
        AND Test_ID != 'BBa_K1679013'
        AND Test_ID NOT IN (
            SELECT Test_ID FROM test_part WHERE Mapping != 'o' OR Mapping IS NULL
        )
    """
    rows = db._execute("Test", query)
    all_ids = [r[0] for r in rows]
    return random.sample(all_ids, min(count, len(all_ids)))

def run_one(test_id, db):
    try:
        ordered = db.get_ordered_parts_by_test_id(test_id)
        if not ordered:
            return None, "No ordered parts"
        cache = {}
        for pid, ptype in ordered:
            d = db.get_part_by_id(SCHEMA, ptype, pid)
            if d is None:
                return None, f"Part missing: {pid}"
            cache[pid] = d

        parts = [cache[pid] for pid, _ in ordered]
        gg = GoldenGateModule()
        # prepare_fragments already prints details
        design = gg.prepare_fragments(parts, assembly_type="linear")

        records = []
        for res in design:
            r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
            r.annotations["target_seq"] = res.get("target_seq", "")
            r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
            r.annotations["part_id"] = res.get("part_id", "")
            r.annotations["part_type"] = res.get("part_type", "")
            records.append(r)

        wf = Workflow(project_dir=PROJECT_DIR, parts=records)
        names = [r.name for r in records]
        module_df = pd.DataFrame({"id": [test_id] * len(names), "name": names})
        wf.run_step3_module_assembly(module_df=module_df, method="golden_gate")

        out = PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb"
        if not out.exists():
            return None, "Output file not found"
            
        return out, "Success"
    except Exception as e:
        return None, str(e)

def main():
    if PROJECT_DIR.exists():
        import shutil
        shutil.rmtree(PROJECT_DIR)
    PROJECT_DIR.mkdir(parents=True, exist_ok=True)
    
    db = DatabaseImport()
    test_ids = get_random_test_ids(db, 10)
    
    log_content = []
    log_content.append(f"# Random Batch Test Report - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log_content.append(f"\n- **Selected Test IDs**: {', '.join(test_ids)}")
    log_content.append("\n---")

    for tid in test_ids:
        log_content.append(f"\n## Test ID: {tid}")
        
        # Capture stdout to log details of prep and assembly
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            out_path, status = run_one(tid, db)
        
        execution_log = f.getvalue()
        if execution_log:
            log_content.append("\n### Execution Log")
            log_content.append("```text")
            log_content.append(execution_log.strip())
            log_content.append("```")

        if out_path:
            recs = DNAFileReader.read_genbank(out_path)
            if recs:
                rec = recs[0]
                log_content.append(f"\n- **Status**: ✅ Success")
                log_content.append(f"- **Result File**: `{out_path.name}`")
                log_content.append(f"- **Total Length**: {len(rec.seq)} bp")
                log_content.append("\n#### Features Found")
                log_content.append("| Type | Label | Location | Note |")
                log_content.append("| :--- | :--- | :--- | :--- |")
                
                for feat in rec.features:
                    label = feat.qualifiers.get("label", ["-"])[0]
                    note = feat.qualifiers.get("note", ["-"])[0]
                    log_content.append(f"| {feat.type} | {label} | {feat.location} | {note} |")
            else:
                log_content.append(f"\n- **Status**: ❌ Failed to read GenBank")
        else:
            log_content.append(f"\n- **Status**: ❌ Failed ({status})")
        
        log_content.append("\n---")

    print("\n".join(log_content))

if __name__ == "__main__":
    main()
