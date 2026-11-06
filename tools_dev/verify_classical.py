"""Classical(general) 경로 검증 하네스 — P5 gfp_test 골든셋.

ClassicalWorkflow(insert+vector→circular construct)를 SnapGene 정답과 대조한다.
insert template = pEGFP(eGFP 함유 plasmid) → prepare_fragments가 그 안의 GFP CDS만 증폭.
  - Gibson:      pEGFP.gb + pIB2.gbk      → vs pIB2-GFP-gibson.gb       (5550+717=6267bp)
  - Golden Gate: pEGFP.gb + pGGAselect.gb → vs pGGAselect-GFP-golden.gb (2854bp, BsaI dropout)

circular 비교는 pydna `seguid()`(circular Dseqrecord면 cdseguid= → 회전·strand 불변) 사용.
출력: `method<TAB>status<TAB>len(res/ans)<TAB>detail`  status=match|MISMATCH|error|no_output

실행(컨테이너, DB 불필요·파일기반):
  docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$PWD/awsemble/src:/app/src:ro" -v "$PWD/awsemble/tools_dev:/app/tools_dev:ro" \
    -v "$PWD/archive/test/gfp_test:/app/data/gfp_test:ro" \
    awsemble:dev python tools_dev/verify_classical.py
"""
import sys, io, contextlib
from pathlib import Path

from utils.io import DNAFileReader
from workflow.classical_assemble import ClassicalWorkflow

DATA = Path("/app/data/gfp_test")

# (method, insert template, vector, 정답 construct)
# insert = pEGFP(eGFP 함유) → prepare_fragments가 GFP CDS만 증폭해 vector에 삽입.
CASES = [
    ("gibson", "snapgene_simulation/gibson/pEGFP.gb", "pIB2.gbk",
     "snapgene_simulation/gibson/pIB2-GFP-gibson.gb"),
    ("golden_gate", "snapgene_simulation/golden_gate/pEGFP.gb",
     "snapgene_simulation/golden_gate/pGGAselect.gb",
     "snapgene_simulation/golden_gate/pGGAselect-GFP-golden.gb"),
]


def _read_one(path: Path):
    """fasta/genbank 단일 레코드 읽기."""
    p = Path(path)
    if p.suffix.lower() in (".fa", ".fasta"):
        recs = DNAFileReader.read_fasta(p)
    else:
        recs = DNAFileReader.read_genbank(p)
    return recs[0] if recs else None


def run_case(method, insert_f, vector_f, answer_f):
    # build-plan 경로: run() → run_plan([vector_insertion]) 전체를 검증
    wf = ClassicalWorkflow(output_dir=f"/tmp/classical/{method}")
    answer = _read_one(DATA / answer_f)
    ans_len = len(answer.seq)
    ans_guid = answer.seguid()

    try:
        results = wf.run(str(DATA / insert_f), str(DATA / vector_f), method=method, validate=False)
    except Exception as e:
        return (method, "error", f"-/{ans_len}", str(e).split("\n")[0][:70])

    if not results:
        return (method, "no_output", f"-/{ans_len}", "no result")

    res = results[0]
    status = "match" if res.seguid() == ans_guid else "MISMATCH"
    return (method, status, f"{len(res.seq)}/{ans_len}", "")


def main():
    ok = True
    for case in CASES:
        with contextlib.redirect_stdout(io.StringIO()):
            r = run_case(*case)
        if r[1] != "match":
            ok = False
        print(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}", flush=True)
    print("RESULT:", "ALL_MATCH" if ok else "HAS_FAIL", flush=True)


if __name__ == "__main__":
    main()
