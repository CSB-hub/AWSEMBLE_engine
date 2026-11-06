"""util_function_add 테스트: KRIBB+igem 섞은 module — 서열입력 vs gbk입력 동일성.

같은 부품(KRIBB promoter/rbs/terminator + igem cds)을 두 경로로 입력해 동일 construct가
나오는지 검증한다. "DB에 없다고 가정" → gbk 추출 시 db=None(매칭 생략, 유저 부품 그대로).

결과: awsemble/tests/user_input_test/mixed/

실행(컨테이너):
  docker run --rm --network host --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/awsemble/tests:/app/tests" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/test_part_intake_mixed.py
"""
import contextlib
import io
from pathlib import Path

from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAFileWriter
from utils.part_intake import parts_from_sequences, parts_from_genbank
from workflow.builder import WorkflowBuilder

OUT = Path("/app/tests/user_input_test/mixed")
SCHEMA = "Bioparts"
VECTOR = "/app/gfp_test/pUC19.gb"
GBK_TYPE = {"promoter": "promoter", "rbs": "RBS", "cds": "CDS", "terminator": "terminator"}
# KRIBB + igem 섞은 single module (KRIBB promoter/rbs/terminator + igem cds)
MIX = [("KP0019", "promoter"), ("KR0021", "rbs"), ("BBa_C0620", "cds"), ("KT0022", "terminator")]


def make_gbk(mix_meta, db, out_path):
    full, feats, pos = "", [], 0
    for pid, pt in mix_meta:
        seq = str(db.get_part_by_id(SCHEMA, pt, pid).seq)
        full += seq
        feats.append((pos, pos + len(seq), GBK_TYPE.get(pt, pt), pid))
        pos += len(seq)
    rec = Dseqrecord(full, name="mixgbk")
    for s, e, t, lbl in feats:
        rec.features.append(SeqFeature(FeatureLocation(s, e), type=t, qualifiers={"label": [lbl]}))
    DNAFileWriter.write_genbank(rec, str(out_path))


def main():
    db = DatabaseImport()
    OUT.mkdir(parents=True, exist_ok=True)
    vector = DNAFileReader.read_genbank(VECTOR)[0]
    seqs = [str(db.get_part_by_id(SCHEMA, pt, pid).seq) for pid, pt in MIX]
    types = [pt for _, pt in MIX]
    names = [pid for pid, _ in MIX]
    L = ["# KRIBB+igem 섞은 module — 서열입력 vs gbk입력 동일성\n",
         f"- 부품(섞음): {[(pid, pt) for pid, pt in MIX]}\n"]

    # 1) 서열 직접 입력 → 조립
    p1 = parts_from_sequences(seqs, types, names)
    with contextlib.redirect_stdout(io.StringIO()):
        o1 = WorkflowBuilder.build_plasmid(p1, vector, db, project_dir=str(OUT / "seq"), insert_method="gibson")
    c1 = o1["construct"]
    L.append(f"\n## 1) 서열 직접 입력\n- construct {len(c1.seq)}bp circular={c1.circular} "
             f"(spacer {o1['l_spacer']}/{o1['r_spacer']})\n")

    # 2) gbk 입력 (DB 없다 가정 = db=None 매칭 생략) → 조립
    gbk = OUT / "mix.gb"
    make_gbk(MIX, db, gbk)
    ext = parts_from_genbank(gbk, db=None)
    p2 = [e["part"] for e in ext]
    with contextlib.redirect_stdout(io.StringIO()):
        o2 = WorkflowBuilder.build_plasmid(p2, vector, db, project_dir=str(OUT / "gbk"), insert_method="gibson")
    c2 = o2["construct"]
    L.append(f"\n## 2) gbk 입력 (DB 없다 가정, db=None)\n- 추출 {len(ext)}부품 "
             f"status={[e['status'] for e in ext]} name={[e['part'].name for e in ext]}\n"
             f"- construct {len(c2.seq)}bp circular={c2.circular} (spacer {o2['l_spacer']}/{o2['r_spacer']})\n")

    # 동일성
    same = str(c1.seq).upper() == str(c2.seq).upper()
    L.append(f"\n## 동일성 검증\n- 서열입력 vs gbk입력 construct: **{'동일 ✅' if same else '상이 ❌'}** "
             f"(len {len(c1.seq)} vs {len(c2.seq)})\n")

    # 참고) gbk + db 매칭 (DB에 있을 때 어떻게 잡히는지)
    ext_db = parts_from_genbank(gbk, db)
    L.append("\n## 참고) gbk + db 매칭 결과\n")
    L.append("  | label | part_type | status | db_id |\n  |---|---|---|---|\n")
    for e in ext_db:
        L.append(f"  | {e['label']} | {e['part_type']} | {e['status']} | {e['db_id']} |\n")

    (OUT / "report.md").write_text("".join(L), encoding="utf-8")
    print("".join(L))


if __name__ == "__main__":
    main()
