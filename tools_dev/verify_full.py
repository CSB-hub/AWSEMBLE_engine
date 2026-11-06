"""풀 회귀 (실제 성공기준 = check_parts_correct: 부품 서열+입력순서) — iGEM 전수.

test_all.py는 성공기준을 "전체 서열 exact-match(compare_sequences)"로 두어
overhang/method 차이 때문에 전건 fail(0%)이 나온다. 본 하네스는 ljw가 확정한
실제 기준(부품+순서)으로 평가해 역사적 baseline(성공 ~9222)과 비교 가능하게 한다.

- 옛/새 코드 모두에서 동작 (workflow import 폴백).
- per-id TSV(id<TAB>status<TAB>seqSHA<TAB>error)를 stdout으로, 총계를 stderr로.
- PROJECT_DIR=/tmp/fullreg 로 격리. LIMIT 환경변수로 소규모 검증 가능.

실행:
  docker run --rm --network host --user $(id -u):$(id -g) -e HOME=/tmp \
    -v .env:/app/.env:ro -v tools_dev:/app/tools_dev:ro [-v src:/app/src:ro] \
    awsemble:dev python tools_dev/verify_full.py > out.tsv 2> out.summary
"""
import sys, os, hashlib, io, contextlib
from pathlib import Path
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAValidator
from assembly.golden_gate import GoldenGateModule
try:
    from workflow.biopart_assemble import BiopartWorkflow as Workflow
except ImportError:  # 옛 코드(pre-P0) 호환
    from workflow.dbtl import DBTLWorkflow as Workflow

PROJECT_DIR = Path("/tmp/fullreg/test_results")
SCHEMA = "Bioparts"
MAX_WORKERS = 20


def run_one(test_id, ordered, parts_cache):
    """단일 Test_ID 조립 + check_parts_correct. (id, status, sha, error)."""
    try:
        # 주의: verbose stdout 억제는 main()에서 전역 1회 처리(스레드 안전).
        #       redirect_stdout는 스레드별로 안전하지 않아 여기서 쓰지 않는다.
        parts = []
        for pid, ptype in ordered:
            if pid not in parts_cache:
                return (test_id, "fail", "-", f"part missing: {pid}")
            parts.append(parts_cache[pid])
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
        wf = Workflow(project_dir=PROJECT_DIR, parts=records)
        mdf = pd.DataFrame({"id": [test_id] * len(names), "name": names})
        wf.run_step3_module_assembly(module_df=mdf, method="golden_gate")
        out = PROJECT_DIR / "modules" / "assembled" / f"{test_id}.gb"
        if not out.exists():
            return (test_id, "no_output", "-", "output not generated")
        recs = DNAFileReader.read_genbank(out)
        seq = str(recs[0].seq) if recs else ""
        sha = hashlib.sha256(seq.lower().encode()).hexdigest()[:16]
        valid = DNAValidator(DatabaseImport()).check_parts_correct(out, test_id)
        return (test_id, "pass" if valid else "fail", sha, "")
    except Exception as e:
        msg = str(e).split("\n")[0]
        # "제한효소 충돌"은 에러가 아니라 의도된 예외처리(BsmBI 등 다른 방법 제안) → 별도 분류
        if "제한효소 충돌" in str(e):
            return (test_id, "alt_method", "-", msg)
        # E-1 입력검증 거부(비-ATGC/빈 서열) → 의도된 불량입력 거부 → 별도 분류
        if "E-1" in str(e) or "ATGC 외" in str(e):
            return (test_id, "invalid_input", "-", msg)
        return (test_id, "error", "-", msg)


def main():
    for s in ["parts/insert", "modules/assembled"]:
        (PROJECT_DIR / s).mkdir(parents=True, exist_ok=True)

    # TSV는 원래 stdout으로 출력하고, 조립/검증의 verbose stdout은 전역 1회 devnull로 억제
    # (스레드 시작 전에 한 번만 교체 → 스레드 안전). 진행/총계는 stderr로.
    real_out = sys.stdout
    sys.stdout = open(os.devnull, "w")
    try:
        db = DatabaseImport()
        ids = sorted(db.get_test_ids_by_institution("igem"))
        limit = os.getenv("LIMIT")
        if limit:
            ids = ids[:int(limit)]

        # 사전 로딩 (정의 + 부품)
        defs, need = {}, set()
        for tid in ids:
            op = db.get_ordered_parts_by_test_id(tid)
            if op:
                defs[tid] = op
                need.update(op)
        cache = {}
        for pid, ptype in need:
            d = db.get_part_by_id(SCHEMA, ptype, pid)
            if d is not None:
                cache[pid] = d
        print(f"[preload] ids={len(ids)} defs={len(defs)} parts={len(cache)}", file=sys.stderr, flush=True)

        results = []
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
            futs = {ex.submit(run_one, tid, defs[tid], cache): tid for tid in defs}
            for i, f in enumerate(futs, 1):
                results.append(f.result())
                if i % 500 == 0:
                    print(f"[progress] {i}/{len(defs)}", file=sys.stderr, flush=True)

        results.sort()
        for r in results:
            print(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}", file=real_out, flush=True)
        c = Counter(r[1] for r in results)
        print(f"[TOTAL] {len(results)} {dict(c)}", file=sys.stderr, flush=True)
    finally:
        sys.stdout.close()
        sys.stdout = real_out


if __name__ == "__main__":
    main()
