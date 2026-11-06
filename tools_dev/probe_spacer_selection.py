"""[1] 선정 프로토콜 설계 실증 — anneal 적합성이 (a) 조합 의존인가 (b) proxy concat으로 충분한가.

검증 2건:
  (1) 조합 의존성: module A,B의 165-spacer 적합집합이 다른가 → 다르면 조합별 캐싱 정당.
  (2) proxy 타당성: 가조립 module 서열 vs 부품 원본 concat 의 적합집합이 같은가
      → 같으면 선정 dry-run을 조립 전 concat(proxy)로 해도 무방(닭-달걀 회피).

실행:
  docker run --rm --network host -v "$PWD/awsemble/.env:/app/.env:ro" \
    -v "$PWD/awsemble/src:/app/src" -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    awsemble:dev python tools_dev/probe_spacer_selection.py BBa_E0240 BBa_C0630
"""
import sys, contextlib, io
from pathlib import Path
import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader
from assembly.golden_gate import GoldenGateModule
from workflow.biopart_assemble import BiopartWorkflow
from basic.primer_design import PrimerDesign

PROJECT_DIR = Path("/tmp/probe/biopart")
SCHEMA = "Bioparts"


def assembled_module_seq(test_id, db):
    """실제 GG 가조립 module 서열 (test_spacer_anneal.build_module 동등)."""
    ordered = db.get_ordered_parts_by_test_id(test_id)
    parts = [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in ordered]
    gg = GoldenGateModule()
    with contextlib.redirect_stdout(io.StringIO()):
        design = gg.prepare_fragments(parts, assembly_type="linear")
        recs, names = [], []
        for r in design:
            d = Dseqrecord(r["full_sequence"], name=r["part_name"])
            d.annotations["target_seq"] = r.get("target_seq", "")
            d.annotations["restriction_enzyme"] = r.get("restriction_enzyme", "BsaI")
            recs.append(d); names.append(d.name)
        wf = BiopartWorkflow(project_dir=PROJECT_DIR, parts=recs)
        wf.run_step3_module_assembly(
            module_df=pd.DataFrame({"id": [test_id]*len(names), "name": names}),
            method="golden_gate")
    m = DNAFileReader.read_genbank(PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb")[0]
    return str(m.seq)


def proxy_concat_seq(test_id, db):
    """부품 원본 서열 단순 concat (조립 전 = 선정 시점에 알 수 있는 정보)."""
    ordered = db.get_ordered_parts_by_test_id(test_id)
    return "".join(str(db.get_part_by_id(SCHEMA, pt, pid).seq) for pid, pt in ordered)


def fit_set(module_seq, spacers):
    """[spacer+module+spacer] 양끝 anchor 프라이머가 설계되는 spacer ID 집합."""
    ok = set()
    for sid, sseq in spacers:
        seq = sseq + module_seq + sseq
        with contextlib.redirect_stdout(io.StringIO()):
            res = PrimerDesign.design_basic_primers(
                template=seq, target_start=0, target_end=len(seq),
                anchor_ends=True, check_specificity=False)
        if res.get("forward") and res.get("reverse"):
            ok.add(sid)
    return ok


def main():
    test_ids = sys.argv[1:] or ["BBa_E0240", "BBa_C0630"]
    db = DatabaseImport()
    spacers = db._execute(
        "Bioparts", "SELECT ID,Sequence FROM spacer WHERE Institution='Awsemble' ORDER BY ID")
    print(f"spacer pool: {len(spacers)}")

    sets_assembled = {}
    for tid in test_ids:
        a = assembled_module_seq(tid, db)
        p = proxy_concat_seq(tid, db)
        fa = fit_set(a, spacers)
        fp = fit_set(p, spacers)
        sets_assembled[tid] = fa
        print(f"\n[{tid}] assembled={len(a)}bp proxy={len(p)}bp")
        print(f"  적합(assembled)={len(fa)}  적합(proxy)={len(fp)}")
        print(f"  (2)proxy 타당성: assembled△proxy = {len(fa ^ fp)} "
              f"{'✅동일' if fa == fp else '⚠️상이 '+str(sorted(fa ^ fp)[:10])}")

    if len(test_ids) >= 2:
        ids = list(sets_assembled)
        a0, a1 = sets_assembled[ids[0]], sets_assembled[ids[1]]
        print(f"\n(1)조합 의존성: {ids[0]} vs {ids[1]}")
        print(f"  |A|={len(a0)} |B|={len(a1)} 공통={len(a0&a1)} "
              f"A\\B={len(a0-a1)} B\\A={len(a1-a0)}")
        print(f"  → {'⚠️ 조합 의존 있음(조합별 캐싱 정당)' if a0 != a1 else '✅ 조합 독립(spacer 고유속성)'}")
        # 첫 적합 spacer가 조합마다 같은지 (선정 결과 동일성)
        first = lambda s: sorted(s)[:2]
        print(f"  첫 적합 L,R: {ids[0]}={first(a0)}  {ids[1]}={first(a1)}")


if __name__ == "__main__":
    main()
