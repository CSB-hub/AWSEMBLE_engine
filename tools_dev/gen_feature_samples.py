"""P6 feature 기입 결과물 샘플 생성 (팀 공유용).

- module feature (annotate_module_parts): 대표 module을 biopart 조립 → 부품 본체 + 접합 overhang feature.
- construct feature (annotate_inserted_region): pEGFP의 GFP CDS를 vector에 삽입 → insert 부품 feature.
산출물은 /out(호스트 마운트)에 저장돼 그대로 남는다.

실행(컨테이너):
  docker run --rm --network host --user $(id -u):$(id -g) -e HOME=/tmp \
    -v "$PWD/awsemble/.env:/app/.env:ro" \
    -v "$PWD/awsemble/src:/app/src:ro" -v "$PWD/awsemble/tools_dev:/app/tools_dev:ro" \
    -v "$PWD/archive/test/gfp_test:/app/data/gfp_test:ro" \
    -v "$PWD/develop_log/ljw/p6_feature_samples:/out" \
    awsemble:dev python tools_dev/gen_feature_samples.py
"""
from pathlib import Path
import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow
from workflow.classical_assemble import ClassicalWorkflow

OUT = Path("/out")
DATA = Path("/app/data/gfp_test")
SCHEMA = "Bioparts"

# 대표 module (부품 수 3~7 다양하게 — baseline QVS 셋에서)
MODULE_IDS = ["BBa_C0630", "BBa_C0460", "BBa_E0240", "BBa_C0480", "BBa_I0402"]

# construct (insert=pEGFP의 GFP CDS, P5 골든셋과 동일 입력)
CONSTRUCT_CASES = [
    ("gibson", "snapgene_simulation/gibson/pEGFP.gb", "pIB2.gbk"),
    ("golden_gate", "snapgene_simulation/golden_gate/pEGFP.gb",
     "snapgene_simulation/golden_gate/pGGAselect.gb"),
]


def feat_summary(gb_path):
    """gbk feature 요약 한 줄."""
    from utils.io import DNAFileReader
    rec = DNAFileReader.read_genbank(gb_path)[0]
    from collections import Counter
    types = Counter(f.type for f in rec.features)
    labels = [f.qualifiers.get("label", [""])[0] for f in rec.features if f.qualifiers.get("label")]
    return f"{len(rec.seq)}bp / features={len(rec.features)} {dict(types)} / labels={labels[:8]}"


def build_modules(db):
    print("=== module feature 샘플 (annotate_module_parts) ===")
    wf_dir = OUT / "modules_work"
    for tid in MODULE_IDS:
        ordered = db.get_ordered_parts_by_test_id(tid)
        if not ordered:
            print(f"  [{tid}] no def — skip"); continue
        parts = [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in ordered]
        if any(p is None for p in parts):
            print(f"  [{tid}] part missing — skip"); continue
        gg = GoldenGateModule()
        design = gg.prepare_fragments(parts, assembly_type="linear")
        records, names = [], []
        for res in design:
            r = Dseqrecord(res["full_sequence"], name=f"{res['part_name']}_designed")
            r.annotations["target_seq"] = res.get("target_seq", "")
            r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
            r.annotations["part_id"] = res.get("part_id", "")
            r.annotations["part_type"] = res.get("part_type", "")
            records.append(r); names.append(r.name)
        wf = BiopartWorkflow(project_dir=wf_dir, parts=records)
        mdf = pd.DataFrame({"id": [tid] * len(names), "name": names})
        wf.run_step3_module_assembly(module_df=mdf, method="golden_gate")
        src = wf_dir / "modules" / "assembled" / f"{tid}.gb"
        dst = OUT / "modules" / f"{tid}.gb"
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.exists():
            dst.write_bytes(src.read_bytes())
            print(f"  [{tid}] {feat_summary(dst)}")
        else:
            print(f"  [{tid}] 생성 실패")


def build_constructs():
    print("\n=== construct feature 샘플 (annotate_inserted_region) ===")
    for method, insert_f, vector_f in CONSTRUCT_CASES:
        wf = ClassicalWorkflow(output_dir=str(OUT / f"constructs_work_{method}"))
        try:
            results = wf.run(str(DATA / insert_f), str(DATA / vector_f),
                             method=method, validate=False)
        except Exception as e:
            print(f"  [{method}] 실패: {e}"); continue
        if not results:
            print(f"  [{method}] no output"); continue
        # constructs/ 에 생성된 gb를 OUT/constructs로 복사
        work = OUT / f"constructs_work_{method}" / "constructs"
        gbs = sorted(work.glob("*.gb")) if work.exists() else []
        dstdir = OUT / "constructs"; dstdir.mkdir(parents=True, exist_ok=True)
        for g in gbs:
            dst = dstdir / f"{method}_{g.name}"
            dst.write_bytes(g.read_bytes())
            print(f"  [{method}] {dst.name}: {feat_summary(dst)}")


def main():
    db = DatabaseImport()
    build_modules(db)
    build_constructs()
    print("\n저장 위치: develop_log/ljw/p6_feature_samples/{modules,constructs}/")


if __name__ == "__main__":
    main()
