"""golden→gibson 최종 풀 테스트 하네스 (P7).

part→module(Golden Gate) / module 이후(Gibson) 루트를 igem 전수 + multi-module + KRIBB로
돌려 construct·report를 만든다. 결과는 seqSHA가 아니라 **전체 정보**(결합 구성/프라이머/전체 서열)로 기록.
메인코드 변경 없음 — 기존 WorkflowBuilder.build_plasmid / build_multimodule_plasmid를 그대로 호출.

모드:
  --mode sample : 케이스 1개(igem BBa_C0630 × pUC19 single)만 → report 형식 확인용
  --mode full   : igem 9999×벡터3 single(병렬) + igem 랜덤 multi(2·3) + KRIBB(pathway7, 4+3)

실행(컨테이너):
  docker run --rm --network host \
    -v "$PWD/awsemble/.env:/app/.env:ro" -v "$PWD/awsemble/src:/app/src" \
    -v "$PWD/awsemble/tools_dev:/app/tools_dev" \
    -v "$PWD/awsemble/tests:/app/tests" \
    -v "$PWD/archive/test/gfp_test:/app/gfp_test:ro" \
    awsemble:dev python tools_dev/fulltest_golden_to_gibson.py --mode sample
"""
import os
import sys
import csv
import random
import argparse
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from utils.io import DatabaseImport, DNAFileReader, DNAFileWriter
from workflow.builder import WorkflowBuilder

# 워커 스레드별 DB connection (mysql connection은 스레드 안전이 아님 → 공유 시 segfault)
_tls = threading.local()


def _db():
    if not hasattr(_tls, "conn"):
        _tls.conn = DatabaseImport()
    return _tls.conn

SCHEMA = "Bioparts"
BASE = Path("/app/tests/golden_to_gibson")
VECTORS = {
    "pUC19": "/app/gfp_test/pUC19.gb",
    "pIB2": "/app/gfp_test/pIB2.gbk",
    "addgene54762": "/app/gfp_test/addgene-plasmid-54762-sequence-212808.gbk",
}
KRIBB_IDS = ["VA16-MvaK1-VA3", "VA3-MvaD-VA4", "VA4-MvaK2-VA5", "VA5-Idi-VA6",
             "VA6-MvaE-VA7", "VA7-MvaS-VA8", "VA8-IspA-VA2"]
MAX_WORKERS = 20
MULTI_N = 10        # 랜덤 multi 세트 수 (2-module, 3-module 각)
SEED = 20260616

INDEX_FIELDS = ["case_id", "group", "inputs", "vector", "method", "L", "R",
                "n_modules", "construct_len", "circular", "status", "reason"]


def _p(*a):
    """진행 로그를 stderr로 (build의 대량 stdout과 분리)."""
    print(*a, file=sys.stderr, flush=True)


# index.csv 증분 기록 (완료마다 flush → 중간에 죽어도 진행분 보존)
_idx_f = None
_idx_w = None


def _index_open():
    global _idx_f, _idx_w
    _idx_f = open(BASE / "index.csv", "w", newline="", encoding="utf-8")
    _idx_w = csv.DictWriter(_idx_f, fieldnames=INDEX_FIELDS)
    _idx_w.writeheader()
    _idx_f.flush()


def _index_append(row):
    _idx_w.writerow(row)
    _idx_f.flush()


def parts_of(test_id, db):
    """test_id의 core 부품 메타[(pid,type)]와 Dseqrecord 리스트 (spacer 제외)."""
    ordered = db.get_ordered_parts_by_test_id(test_id)
    meta = [(pid, pt) for pid, pt in ordered if pt != "spacer"]
    return meta, [db.get_part_by_id(SCHEMA, pt, pid) for pid, pt in meta]


def _sseq(db, sid):
    return str(db.get_part_by_id(SCHEMA, "spacer", sid).seq)


def _csv_table(path, title):
    """csv → markdown 표 라인 리스트 (파일 없으면 빈 리스트)."""
    if not path.exists():
        return []
    rows = list(csv.reader(open(path, encoding="utf-8")))
    if not rows:
        return []
    out = [f"\n## {title}\n\n", "| " + " | ".join(rows[0]) + " |\n",
           "|" + "---|" * len(rows[0]) + "\n"]
    for r in rows[1:]:
        out.append("| " + " | ".join(r) + " |\n")
    return out


def write_report(case_dir, info, construct):
    """report.md — 결합정보 + 조립결과(전체 서열) + 프라이머(primers.csv 표). SHA 없이 전체."""
    L = [f"# {info['case_id']}\n\n## 결합 정보\n"]
    L.append(f"- 입력 module: {', '.join(info['modules'])}\n")
    L.append("- 부품 구성:\n\n  | module | part_id | type |\n  |---|---|---|\n")
    for mod, pid, pt in info["parts"]:
        L.append(f"  | {mod} | {pid} | {pt} |\n")
    L.append("- spacer:\n\n  | 위치 | id | 서열 |\n  |---|---|---|\n")
    for label, sid, seq in info["spacers"]:
        L.append(f"  | {label} | {sid} | {seq} |\n")
    L.append(f"- vector: {info['vector']} ({info['vector_len']}bp)\n")
    L.append(f"- method: {info['method']}\n")
    L.append(f"\n## 조립 결과\n- construct: {len(construct.seq)}bp, circular={construct.circular}\n")
    L.append(f"\n### 전체 서열\n```\n{str(construct.seq)}\n```\n")

    L += _csv_table(case_dir / "reports" / "primers.csv", "프라이머 (module 조립, Golden Gate)")
    L += _csv_table(case_dir / "reports" / "vector_primers.csv", "프라이머 (vector 삽입, Gibson)")

    DNAFileWriter.write_genbank(construct, str(case_dir / "construct.gb"))
    (case_dir / "report.md").write_text("".join(L), encoding="utf-8")


def _row(case_id, group, inputs, vector, method, L="", R="", n=1,
         clen="", circular="", status="ok", reason=""):
    return dict(case_id=case_id, group=group, inputs=inputs, vector=vector, method=method,
                L=L, R=R, n_modules=n, construct_len=clen, circular=circular,
                status=status, reason=reason)


def run_single(test_id, vname, vpath, sub="igem/single"):
    db = _db()
    case_id = f"{test_id}__{vname}"
    case_dir = BASE / sub / case_id
    try:
        meta, parts = parts_of(test_id, db)
        vector = DNAFileReader.read_genbank(vpath)[0]
        out = WorkflowBuilder.build_plasmid(
            parts, vector, db, project_dir=str(case_dir), insert_method="gibson")
        con = out["construct"]
        Lsp, Rsp = out["l_spacer"], out["r_spacer"]
        info = dict(case_id=case_id, modules=[test_id],
                    parts=[(test_id, pid, pt) for pid, pt in meta],
                    spacers=[("L", Lsp, _sseq(db, Lsp)), ("R", Rsp, _sseq(db, Rsp))],
                    vector=vname, vector_len=len(vector.seq), method="gibson")
        write_report(case_dir, info, con)
        return _row(case_id, "igem-single", test_id, vname, "gibson", Lsp, Rsp, 1,
                    len(con.seq), con.circular, "ok")
    except Exception as e:
        return _row(case_id, "igem-single", test_id, vname, "gibson",
                    status="fail", reason=f"{type(e).__name__}: {str(e)[:90]}")


def run_multi(set_id, test_ids, vname, vpath, group, sub):
    db = _db()
    case_id = f"{set_id}__{vname}"
    case_dir = BASE / sub / case_id
    try:
        groups, parts_meta = [], []
        for tid in test_ids:
            meta, parts = parts_of(tid, db)
            groups.append(parts)
            parts_meta += [(tid, pid, pt) for pid, pt in meta]
        vector = DNAFileReader.read_genbank(vpath)[0]
        out = WorkflowBuilder.build_multimodule_plasmid(
            groups, vector, db, project_dir=str(case_dir), insert_method="gibson")
        con = out["construct"]
        bnd = out["boundaries"]   # [L, *connectors, R]
        labels = ["L"] + [f"연결{i+1}" for i in range(len(bnd) - 2)] + ["R"]
        info = dict(case_id=case_id, modules=test_ids, parts=parts_meta,
                    spacers=[(lab, sid, _sseq(db, sid)) for lab, sid in zip(labels, bnd)],
                    vector=vname, vector_len=len(vector.seq), method="gibson")
        write_report(case_dir, info, con)
        return _row(case_id, group, "+".join(test_ids), vname, "gibson",
                    bnd[0], bnd[-1], len(test_ids), len(con.seq), con.circular, "ok")
    except Exception as e:
        return _row(case_id, group, "+".join(test_ids), vname, "gibson", n=len(test_ids),
                    status="fail", reason=f"{type(e).__name__}: {str(e)[:90]}")


def write_index(rows):
    """증분 기록된 index.csv를 닫고 최종 집계 출력."""
    if _idx_f:
        _idx_f.close()
    ok = sum(1 for r in rows if r["status"] == "ok")
    _p(f"\n=== index.csv: {len(rows)}건 (ok {ok} / fail {len(rows)-ok}) ===")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["sample", "pre", "full"], default="sample")
    args = ap.parse_args()
    # build_plasmid의 대량 stdout을 프로세스 전역으로 억제(멀티스레드 redirect race 방지).
    # 진행 로그는 _p()로 stderr에 따로 출력한다.
    sys.stdout = open(os.devnull, "w")
    db = _db()   # main 스레드용 (igem 목록 조회). 워커는 각자 _db()

    if args.mode == "sample":
        r = run_single("BBa_C0630", "pUC19", VECTORS["pUC19"])
        _p(f"sample: {r['case_id']} status={r['status']} "
           f"construct={r['construct_len']}bp circular={r['circular']}")
        _p(f"  → {BASE}/igem/single/{r['case_id']}/report.md")
        return

    if args.mode == "pre":
        # 풀테스트 전 3가지 multi 케이스만 먼저 검증
        rows = []
        igem = db.get_test_ids_by_institution("igem")
        rng = random.Random(SEED)
        # 1) igem 랜덤 2 module → multimodule → pUC19
        tids = rng.sample(igem, 2)
        _p(f"[pre] 1) igem 랜덤 2-module {tids} → pUC19")
        rows.append(run_multi("pre_igem2", tids, "pUC19", VECTORS["pUC19"], "igem-multi", "igem/multi"))
        # 2) KRIBB 7-module pathway → pUC19 (하나의 벡터)
        _p("[pre] 2) KRIBB 7-module pathway → pUC19")
        rows.append(run_multi("pathway7", KRIBB_IDS, "pUC19", VECTORS["pUC19"], "kribb-pathway7", "kribb"))
        # 3) KRIBB 4+3 쪼개서 2개의 벡터 (4→pUC19, 3→pIB2)
        _p("[pre] 3) KRIBB 4→pUC19, 3→pIB2")
        rows.append(run_multi("kribb_4", KRIBB_IDS[:4], "pUC19", VECTORS["pUC19"], "kribb-multi43", "kribb/multi_4_3"))
        rows.append(run_multi("kribb_3", KRIBB_IDS[4:], "pIB2", VECTORS["pIB2"], "kribb-multi43", "kribb/multi_4_3"))
        _p("\n=== pre 결과 ===")
        for r in rows:
            _p(f"  {r['case_id']}: {r['status']} {r['construct_len']}bp "
               f"n_modules={r['n_modules']} {r['reason']}")
        return

    # ---- full ----
    rows = []
    igem = db.get_test_ids_by_institution("igem")
    _index_open()
    _p(f"[full] igem single: {len(igem)} × 벡터 {len(VECTORS)} (병렬 {MAX_WORKERS})")
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        futs = [ex.submit(run_single, tid, vn, vp)
                for tid in igem for vn, vp in VECTORS.items()]
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            rows.append(r)
            _index_append(r)
            if i % 500 == 0:
                _p(f"  igem single {i}/{len(futs)}")

    _p("[full] igem multi (랜덤 2·3 module)")
    _tls.conn = DatabaseImport()   # single(장시간) 동안 idle된 main db 재생성(2013 방지)
    rng = random.Random(SEED)
    for k in range(MULTI_N):
        for size in (2, 3):
            tids = rng.sample(igem, size)
            for vn, vp in VECTORS.items():
                r = run_multi(f"set{size}_{k:02d}", tids, vn, vp, "igem-multi", "igem/multi")
                rows.append(r)
                _index_append(r)

    _p("[full] KRIBB pathway7 + multi_4_3")
    for vn, vp in VECTORS.items():
        r = run_multi("pathway7", KRIBB_IDS, vn, vp, "kribb-pathway7", "kribb")
        rows.append(r)
        _index_append(r)
    r = run_multi("kribb_4", KRIBB_IDS[:4], "pUC19", VECTORS["pUC19"], "kribb-multi43", "kribb/multi_4_3")
    rows.append(r)
    _index_append(r)
    r = run_multi("kribb_3", KRIBB_IDS[4:], "pUC19", VECTORS["pUC19"], "kribb-multi43", "kribb/multi_4_3")
    rows.append(r)
    _index_append(r)

    write_index(rows)


if __name__ == "__main__":
    main()
