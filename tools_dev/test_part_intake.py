"""util_function_add 테스트: part_intake (유저 입력 → 부품) 검증.

  1) 유저 서열 직접 입력 → Dseqrecord 부품 → build_plasmid 조립
  2) DB 부품으로 합성한 gbk → feature 추출 + DB 매칭(matched) → 조립
  3) 실제 plasmid(pEGFP.gb) → feature 추출 + DB 매칭(미등록=new)

결과: awsemble/tests/user_input_test/  (report.md + 합성 gbk + construct)

실행(컨테이너):
  docker run --rm --network host --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/awsemble/tests:/app/tests" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/test_part_intake.py
"""
import contextlib
import io
from pathlib import Path

from Bio.SeqFeature import SeqFeature, FeatureLocation
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAFileWriter
from utils.part_intake import parts_from_sequences, parts_from_genbank
from workflow.builder import WorkflowBuilder

OUT = Path("/app/tests/user_input_test")
SCHEMA = "Bioparts"
VECTOR = "/app/gfp_test/pUC19.gb"
GBK_TYPE = {"promoter": "promoter", "rbs": "RBS", "cds": "CDS", "terminator": "terminator"}


def make_test_gbk(test_id, db, out_path):
    """DB 부품(spacer 제외)을 이어붙이고 part_type을 GenBank feature로 단 plasmid gbk 합성."""
    core = [(pid, pt) for pid, pt in db.get_ordered_parts_by_test_id(test_id) if pt != "spacer"]
    full, feats, pos = "", [], 0
    for pid, pt in core:
        seq = str(db.get_part_by_id(SCHEMA, pt, pid).seq)
        full += seq
        feats.append((pos, pos + len(seq), GBK_TYPE.get(pt, pt), pid))
        pos += len(seq)
    rec = Dseqrecord(full, name="testgbk")
    for s, e, t, lbl in feats:
        rec.features.append(SeqFeature(FeatureLocation(s, e), type=t, qualifiers={"label": [lbl]}))
    DNAFileWriter.write_genbank(rec, str(out_path))
    return core


def _table(rows_dicts):
    L = ["  | label | feature→part_type | status | db_id |\n  |---|---|---|---|\n"]
    for e in rows_dicts:
        L.append(f"  | {str(e['label'])[:22]} | {e['part_type']} | {e['status']} | {e['db_id']} |\n")
    return "".join(L)


def main():
    db = DatabaseImport()
    OUT.mkdir(parents=True, exist_ok=True)
    vector = DNAFileReader.read_genbank(VECTOR)[0]
    L = ["# util_function_add — part_intake 테스트 결과\n"]

    # === 1) 서열 직접 입력 → 부품 → 조립 ===
    tid = "BBa_C0630"
    core = [(pid, pt) for pid, pt in db.get_ordered_parts_by_test_id(tid) if pt != "spacer"]
    seqs = [str(db.get_part_by_id(SCHEMA, pt, pid).seq) for pid, pt in core]
    types = [pt for _, pt in core]
    parts = parts_from_sequences(seqs, types=types, names=[pid for pid, _ in core])
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = WorkflowBuilder.build_plasmid(parts, vector, db,
                                                project_dir=str(OUT / "seq_input"), insert_method="gibson")
        con = out["construct"]
        L.append(f"\n## 1) 서열 직접 입력 → 조립\n- 입력 부품 {len(parts)}개 {types}\n"
                 f"- ✅ construct {len(con.seq)}bp circular={con.circular} (L/R spacer {out['l_spacer']}/{out['r_spacer']})\n")
    except Exception as e:
        L.append(f"\n## 1) 서열 직접 입력 → 조립\n- ❌ {type(e).__name__}: {str(e)[:90]}\n")

    # === 2) 합성 gbk(DB 부품) → feature 추출 + 매칭 → 조립 ===
    gbk = OUT / "synth_from_db.gb"
    make_test_gbk(tid, db, gbk)
    extracted = parts_from_genbank(gbk, db)
    L.append(f"\n## 2) 합성 gbk feature 추출 + DB 매칭\n- gbk: `{gbk.name}` (feature {len(extracted)}개)\n")
    L.append(_table(extracted))
    matched = [e["part"] for e in extracted if e["status"] == "matched"]
    L.append(f"- matched {len(matched)}/{len(extracted)}\n")
    if len(matched) == len(extracted) and matched:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                out2 = WorkflowBuilder.build_plasmid(matched, vector, db,
                                                     project_dir=str(OUT / "gbk_input"), insert_method="gibson")
            c2 = out2["construct"]
            L.append(f"- ✅ 전부 matched → 조립: construct {len(c2.seq)}bp circular={c2.circular}\n")
        except Exception as e:
            L.append(f"- ❌ 조립 실패: {type(e).__name__}: {str(e)[:90]}\n")

    # === 3) 실제 plasmid → feature 추출/매칭 (미등록 케이스) ===
    real_gbk = "/app/gfp_test/addgene-plasmid-54762-sequence-212808.gbk"
    ex3 = parts_from_genbank(real_gbk, db)
    by_status = {}
    for e in ex3:
        by_status[e["status"]] = by_status.get(e["status"], 0) + 1
    L.append(f"\n## 3) 실제 plasmid (addgene-54762) feature 추출\n- feature {len(ex3)}개, 상태 {by_status}\n")
    L.append(_table([e for e in ex3 if e["status"] != "unknown_type"][:15]))

    (OUT / "report.md").write_text("".join(L), encoding="utf-8")
    print("".join(L))


if __name__ == "__main__":
    main()
