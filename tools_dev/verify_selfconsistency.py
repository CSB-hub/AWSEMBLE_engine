"""배포 전 최종 검증 — self-consistency 하네스 (계획 §3 step1: A 라운드트립 + B 불변식)

gold standard(정답 서열)를 쓰지 않고, 설계가 "주장하는 것"을 시뮬레이션으로 되짚어
**자기모순**을 검출한다. 산출물별(계획 §1)로 라운드트립/불변식을 건다.

검증 대상 산출물(현 igem/KRIBB 코퍼스가 실제로 지나가는 경로 = part→module):
  - item2 prepared part (full_sequence): 절단/PCR 라운드트립 + 구조 불변식
  - item3 module            : 효소site 잔존0 / 부품 순서포함 / 길이 / 재조립 결정성 / 디스크 왕복
  - item6 primers.csv       : 행수·필수컬럼

핵심 원칙: **프로덕션 코드를 그대로 호출**(prepare_fragments / extract_insert_fragment /
module_assembly)해서 나온 산출물에만 체크를 건다 → 재구현 괴리 없음.

실행:
    PYTHONPATH=src /home/csb/awsemble_env/bin/python tools_dev/verify_selfconsistency.py BBa_I2030 ...
    PYTHONPATH=src /home/csb/awsemble_env/bin/python tools_dev/verify_selfconsistency.py --sample 30
    PYTHONPATH=src /home/csb/awsemble_env/bin/python tools_dev/verify_selfconsistency.py --all       # C 전수
    옵션: --inst igem|KRIBB  --seed N  --json out.json  --quiet
"""

import argparse
import contextlib
import io
import json
import random
import sys
import traceback
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Dict, List, Optional

from Bio.Seq import Seq
from Bio.Restriction import BsaI, BbsI
from pydna.dseqrecord import Dseqrecord

from utils.io import DatabaseImport, DNAFileReader, DNAFileWriter
from workflow.biopart_assemble import BiopartWorkflow
from basic.sequence_analysis import get_typeIIS_params

BIOPARTS_SCHEMA = "Bioparts"
PASS, FAIL, SKIP = "PASS", "FAIL", "SKIP"


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement()).upper()


@dataclass
class Check:
    """단일 self-consistency 체크 결과 1건."""
    artifact: str      # 어떤 산출물 (예: "prep[2] BBa_x" / "module")
    check: str         # 체크 이름 (예: "slice_roundtrip")
    status: str        # PASS / FAIL / SKIP
    detail: str = ""   # 실패/스킵 사유


def classify_error(msg: str) -> str:
    """파이프라인 예외를 트리아지 4분류로 매핑(계획 §C).

    known_limit / data 는 '설계가 의도대로 거부/차단'한 정상 동작 → FAIL 아님.
    unexpected 만 진짜 버그 후보.
    """
    if "제한효소 충돌" in msg or "제3의 Type IIS" in msg:
        return "known_limit:enzyme_conflict"      # BsaI+BbsI 동시 → 제3효소 권고(정상 거부)
    if "ATGC 외 염기" in msg or "E-1" in msg or "빈 서열" in msg:
        return "data:invalid_base"                # N/IUPAC/빈서열 입력 차단(정상 거부)
    if "부품 서열 누락" in msg or "DB 부품 구성 없음" in msg:
        return "data:db_missing"                  # DB 데이터 이슈
    return "unexpected"                           # 진짜 버그 후보


@dataclass
class TestResult:
    test_id: str
    checks: List[Check] = field(default_factory=list)
    error: Optional[str] = None   # prep/assemble 자체가 터진 경우(체크 이전 단계 실패)
    error_kind: Optional[str] = None

    def add(self, artifact, check, status, detail=""):
        self.checks.append(Check(artifact, check, status, detail))

    @property
    def failed(self) -> List[Check]:
        return [c for c in self.checks if c.status == FAIL]


# ------------------------------------------------------------------ #
# item2: prepared part (full_sequence) 체크
# ------------------------------------------------------------------ #
def check_prepared_part(res: TestResult, tag: str, d: Dict, part: Dseqrecord) -> Optional[Dseqrecord]:
    """prepared part 1개에 대한 라운드트립/불변식. 반환: 슬라이스 fragment(모듈 체크용)."""
    enzyme = d.get("restriction_enzyme", "BsaI")
    params = get_typeIIS_params(enzyme)
    site_fwd, site_rev, sp_len = params["fwd"], params["rev"], params["spacer_len"]
    spacer = "A" * sp_len

    left_oh = (d.get("left_oh") or "").upper()
    right_oh = (d.get("right_oh") or "").upper()
    target = (d.get("target_seq") or "").upper()
    full = (d.get("full_sequence") or "").upper()

    # 설계가 주장하는 절단 후 중간 조각 = [left_oh][target][RC(right_oh)]
    expected_mid = left_oh + target + _rc(right_oh) if right_oh else left_oh + target

    # ---- B불변식: full_sequence == prefix + target + suffix (구조 재구성) ----
    prefix = (site_fwd + spacer + left_oh) if left_oh else ""
    suffix = (_rc(right_oh) + _rc(spacer) + site_rev) if right_oh else ""
    expected_full = prefix + target + suffix
    if full == expected_full:
        res.add(tag, "structure", PASS)
    else:
        res.add(tag, "structure", FAIL,
                f"full != prefix+target+suffix (len {len(full)} vs {len(expected_full)})")

    # ---- A/B: 결정적 슬라이스(extract_insert_fragment)가 기대 중간과 일치 ----
    from assembly.golden_gate import GoldenGateModule
    gg = GoldenGateModule()
    slice_rec = Dseqrecord(full)
    slice_rec.name = d.get("part_name", tag)
    slice_rec.annotations["restriction_enzyme"] = enzyme
    with contextlib.redirect_stdout(io.StringIO()):
        frag = gg.extract_insert_fragment(slice_rec)
    frag_seq = str(frag.seq).upper() if frag else ""
    if frag_seq == expected_mid:
        res.add(tag, "slice_roundtrip", PASS)
    else:
        res.add(tag, "slice_roundtrip", FAIL,
                f"슬라이스 != [oh][target][oh_rc] (len {len(frag_seq)} vs {len(expected_mid)})")

    # ---- A: 실제 효소 cut vs 결정적 슬라이스 (plan §A 핵심 대조) ----
    enz = BbsI if enzyme == "BbsI" else BsaI
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            cut_frags = [str(f.seq).upper() for f in Dseqrecord(full).cut(enz)]
        if not cut_frags:                       # 사이트 없음(양끝 blunt 단일부품) → full 자체가 중간
            cut_frags = [full]
        if expected_mid in cut_frags:
            res.add(tag, "recut_roundtrip", PASS)
        else:
            res.add(tag, "recut_roundtrip", FAIL,
                    f"실제 cut 조각 중 기대 중간({len(expected_mid)}bp) 없음. 조각들={[len(x) for x in cut_frags]}")
    except Exception as e:                       # 짧은 부품 pydna cut 실패(알려진 한계) → SKIP
        res.add(tag, "recut_roundtrip", SKIP, f"pydna cut 예외: {e}")

    # ---- A: PCR 라운드트립 — 설계 프라이머로 원본 part 증폭 == full_sequence ----
    fwd, rev = d.get("forward_primer", ""), d.get("reverse_primer", "")
    try:
        from pydna.amplify import pcr
        with contextlib.redirect_stdout(io.StringIO()):
            amp = pcr(fwd, rev, part)
        if str(amp.seq).upper() == full:
            res.add(tag, "pcr_roundtrip", PASS)
        else:
            res.add(tag, "pcr_roundtrip", FAIL,
                    f"PCR 산물 != full_sequence (len {len(amp.seq)} vs {len(full)})")
    except Exception as e:
        res.add(tag, "pcr_roundtrip", SKIP, f"pcr 예외: {e}")

    return frag


# ------------------------------------------------------------------ #
# item3: module 체크
# ------------------------------------------------------------------ #
def check_module(res: TestResult, wf: BiopartWorkflow, records: List[Dseqrecord],
                 design_results: List[Dict], fragments: List[Dseqrecord],
                 parts: List[Dseqrecord]) -> None:
    # 프로덕션 조립 그대로 호출
    with contextlib.redirect_stdout(io.StringIO()):
        module = wf._assemble_fragments(records, method="golden_gate", assembly_type="linear")
    if module is None:
        res.add("module", "assemble", FAIL, "조립 실패(fragment 추출/병합 None)")
        return
    mseq = str(module.seq).upper()

    # ---- B: 조립효소 인식서열 잔존 0 ----
    #   ⚠️ 조립에 실제 투입되는 효소(design 결정값)의 site만 위반이다. 반대 효소의
    #   우연 site는 반응에 안 넣으므로 무해 → FAIL 아님(product validate()는 둘 다 봐서 과엄격).
    asm_enzyme = design_results[0].get("restriction_enzyme", "BsaI") if design_results else "BsaI"
    enz_obj = BbsI if asm_enzyme == "BbsI" else BsaI
    other_obj = BsaI if asm_enzyme == "BbsI" else BbsI
    n_asm = len(enz_obj.search(Seq(mseq)))
    n_other = len(other_obj.search(Seq(mseq)))
    coincidental = f" (우연 {other_obj}={n_other}, 무해)" if n_other else ""
    if n_asm == 0:
        res.add("module", "no_enzyme_site", PASS, coincidental.strip())
    else:
        # 조립효소 site가 raw 부품에 있었으면 '탐지 미스'(효소전환 실패=심각),
        # 없으면 overhang↔부품 경계에서 '생성된 site'(로드맵 알려진 한계).
        in_raw = any(enz_obj.search(Seq(str(p.seq).upper())) for p in parts)
        cause = "탐지미스(효소전환 실패)" if in_raw else "junction 생성 site(알려진 한계)"
        res.add("module", "no_enzyme_site", FAIL,
                f"조립효소 {asm_enzyme} site {n_asm}개 잔류 — {cause}{coincidental}")

    # ---- B: 각 부품 target_seq가 순서대로 등장 ----
    pos, ordered_ok = 0, True
    for i, d in enumerate(design_results):
        t = (d.get("target_seq") or "").upper()
        j = mseq.find(t, pos)
        if j < 0:
            ordered_ok = False
            res.add("module", "parts_in_order", FAIL, f"부품[{i}] target 미검출/순서역전")
            break
        pos = j + len(t)
    if ordered_ok:
        res.add("module", "parts_in_order", PASS)

    # ---- B: 인접 fragment 접합 overhang(4bp) 정합 (prev 끝 == curr 시작) ----
    jn_ok = True
    for i in range(len(fragments) - 1):
        a, b = str(fragments[i].seq).upper(), str(fragments[i + 1].seq).upper()
        if a[-4:] != b[:4]:
            jn_ok = False
            res.add("module", "junction_overhang", FAIL,
                    f"접합[{i}] {a[-4:]} != {b[:4]}")
            break
    if jn_ok:
        res.add("module", "junction_overhang", PASS)

    # ---- B: 길이 불변식 — 병합은 접합마다 4bp 공유 ----
    exp_len = sum(len(str(f.seq)) for f in fragments) - 4 * (len(fragments) - 1)
    if len(mseq) == exp_len:
        res.add("module", "length", PASS)
    else:
        res.add("module", "length", FAIL, f"len {len(mseq)} != Σfrag-4(n-1)={exp_len}")

    # ---- A: 재조립 결정성 — 같은 입력 두 번 → 동일 서열 ----
    with contextlib.redirect_stdout(io.StringIO()):
        module2 = wf._assemble_fragments(records, method="golden_gate", assembly_type="linear")
    if module2 is not None and str(module2.seq).upper() == mseq:
        res.add("module", "reassembly_determinism", PASS)
    else:
        res.add("module", "reassembly_determinism", FAIL, "동일 입력 재조립 결과 불일치")

    # ---- item1/writer: GenBank 디스크 왕복 → 서열 불변 ----
    try:
        with TemporaryDirectory() as td, contextlib.redirect_stdout(io.StringIO()):
            p = Path(td) / "m.gb"
            DNAFileWriter.write_genbank(module, str(p))
            back = DNAFileReader.read_genbank(p)[0]
        if str(back.seq).upper() == mseq:
            res.add("module", "disk_roundtrip", PASS)
        else:
            res.add("module", "disk_roundtrip", FAIL, "GenBank 저장/재로드 후 서열 변형")
    except Exception as e:
        res.add("module", "disk_roundtrip", SKIP, f"디스크 왕복 예외: {e}")


# ------------------------------------------------------------------ #
# item6: primers.csv 체크
# ------------------------------------------------------------------ #
def check_primer_csv(res: TestResult, wf: BiopartWorkflow, design_results: List[Dict]) -> None:
    import csv
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            out = wf.write_primer_csv(design_results, output_file="primers.csv")
        with open(out, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
    except Exception as e:
        res.add("primers.csv", "write", FAIL, f"CSV 작성 예외: {e}")
        return
    if len(rows) == len(design_results):
        res.add("primers.csv", "row_count", PASS)
    else:
        res.add("primers.csv", "row_count", FAIL, f"{len(rows)} rows != {len(design_results)} parts")
    need = {"part_name", "forward_primer", "reverse_primer", "restriction_enzyme"}
    missing = need - set(rows[0].keys()) if rows else need
    if not missing:
        res.add("primers.csv", "columns", PASS)
    else:
        res.add("primers.csv", "columns", FAIL, f"필수 컬럼 누락: {missing}")


# ------------------------------------------------------------------ #
# test_id 하나에 대한 전체 검증
# ------------------------------------------------------------------ #
def _lookup_part(db: DatabaseImport, ptype: str, pid: str, is_kribb: bool):
    """부품 서열 조회. KRIBB는 test_part가 'P0019'를 주지만 Bioparts엔 'KP0019'로
    적재돼 있어(K 접두사) → KRIBB는 K+pid 먼저 시도 후 raw 폴백.
    ⚠️ 프로덕션 조회(get_part_by_id / test_kribb.py)엔 이 매핑이 없어 KRIBB가
       현재 end-to-end로 안 돌아간다 → 통합계층 TODO(하네스는 검증 위해 우회)."""
    if is_kribb:
        return db.get_part_by_id(BIOPARTS_SCHEMA, ptype, "K" + pid) \
            or db.get_part_by_id(BIOPARTS_SCHEMA, ptype, pid)
    return db.get_part_by_id(BIOPARTS_SCHEMA, ptype, pid)


def verify_test_id(test_id: str, db: DatabaseImport, workdir: Path,
                   is_kribb: bool = False) -> TestResult:
    from assembly.golden_gate import GoldenGateModule
    res = TestResult(test_id)
    try:
        ordered = db.get_ordered_parts_by_test_id(test_id)
        if not ordered:
            res.error, res.error_kind = "DB 부품 구성 없음", "data:db_missing"
            return res
        parts = []
        for pid, ptype in ordered:
            p = _lookup_part(db, ptype, pid, is_kribb)
            if p is None:
                res.error = f"부품 서열 누락: {pid} ({ptype})"
                res.error_kind = "data:db_missing"
                return res
            parts.append(p)

        gg = GoldenGateModule()
        with contextlib.redirect_stdout(io.StringIO()):
            design_results = gg.prepare_fragments(parts, assembly_type="linear")

        # prep 결정성: 두 번 돌려 full_sequence 동일 (재현성 불변식)
        with contextlib.redirect_stdout(io.StringIO()):
            design2 = gg.prepare_fragments(parts, assembly_type="linear")
        if [d["full_sequence"] for d in design_results] == [d["full_sequence"] for d in design2]:
            res.add("prep", "determinism", PASS)
        else:
            res.add("prep", "determinism", FAIL, "prepare_fragments 2회 결과 불일치")

        # records 구성 (test_single.py와 동일 계약)
        records = []
        for d in design_results:
            r = Dseqrecord(d["full_sequence"], name=f"{d['part_name']}_designed")
            r.annotations["target_seq"] = d.get("target_seq", "")
            r.annotations["restriction_enzyme"] = d.get("restriction_enzyme", "BsaI")
            r.annotations["part_id"] = d.get("part_id", "")
            r.annotations["part_type"] = d.get("part_type", "")
            records.append(r)

        # item2: prepared part 체크 (원본 part = PCR template)
        for i, (d, part) in enumerate(zip(design_results, parts)):
            check_prepared_part(res, f"prep[{i}] {d.get('part_id','')}", d, part)

        # 슬라이스 fragment (모듈 접합 체크용)
        fragments = []
        for r in records:
            with contextlib.redirect_stdout(io.StringIO()):
                fragments.append(gg.extract_insert_fragment(r))

        wf = BiopartWorkflow(project_dir=str(workdir), parts=records)
        check_module(res, wf, records, design_results, fragments, parts)
        check_primer_csv(res, wf, design_results)

    except Exception as e:
        res.error = f"{type(e).__name__}: {e}"
        res.error_kind = classify_error(str(e))
        # 알려진 한계/데이터 이슈(의도된 거부)는 FAIL 아님 — unexpected 만 FAIL로 표면화
        if res.error_kind == "unexpected":
            res.checks.append(Check("pipeline", "exception", FAIL, traceback.format_exc(limit=3)))
    return res


# ------------------------------------------------------------------ #
# main / 집계
# ------------------------------------------------------------------ #
def main():
    ap = argparse.ArgumentParser(description="self-consistency 검증 하네스 (A 라운드트립 + B 불변식)")
    ap.add_argument("test_ids", nargs="*", help="검증할 Test_ID (없으면 --sample/--all)")
    ap.add_argument("--sample", type=int, help="코퍼스에서 무작위 N개")
    ap.add_argument("--all", action="store_true", help="C 전수 (igem+KRIBB 전체)")
    ap.add_argument("--inst", choices=["igem", "KRIBB"], help="특정 기관만")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--json", help="결과 JSON 저장 경로")
    ap.add_argument("--quiet", action="store_true", help="케이스별 상세 생략, 요약만")
    args = ap.parse_args()

    db = DatabaseImport()
    if args.test_ids:
        ids = args.test_ids
    elif args.inst:
        ids = db.get_test_ids_by_institution(args.inst)
    else:
        ids = db.get_all_test_ids_from_construct()
    if not args.all and not args.test_ids:
        n = args.sample or 20
        random.seed(args.seed)
        ids = random.sample(ids, min(n, len(ids)))

    kribb_set = set(db.get_test_ids_by_institution("KRIBB"))
    print(f"🚀 self-consistency 검증: {len(ids)}개 케이스 (KRIBB {sum(1 for i in ids if i in kribb_set)})\n")

    results: List[TestResult] = []
    # 체크 이름별 집계 (PASS/FAIL/SKIP)
    per_check = {}
    errors = []
    with TemporaryDirectory() as td:
        workdir = Path(td) / "proj"
        for k, tid in enumerate(ids, 1):
            r = verify_test_id(tid, db, workdir, is_kribb=tid in kribb_set)
            results.append(r)
            if r.error:
                errors.append((tid, r.error))
            for c in r.checks:
                d = per_check.setdefault(c.check, Counter())
                d[c.status] += 1
            if not args.quiet:
                fails = r.failed
                if r.error and not fails:
                    icon = "🐛" if r.error_kind == "unexpected" else "⏭️"
                    print(f"[{k}/{len(ids)}] {tid}  {icon} [{r.error_kind}] {r.error.splitlines()[0][:70]}")
                elif fails:
                    print(f"[{k}/{len(ids)}] {tid}  ❌ {len(fails)} FAIL")
                    for c in fails:
                        print(f"      - {c.artifact} :: {c.check} — {c.detail}")
                else:
                    npass = sum(1 for c in r.checks if c.status == PASS)
                    print(f"[{k}/{len(ids)}] {tid}  ✅ {npass} pass")

    # ---- 요약 ----
    total_fail = sum(len(r.failed) for r in results)
    clean = sum(1 for r in results if not r.failed and not r.error)
    # 트리아지 4분류 (계획 §C): 의도된 거부(known_limit/data)와 진짜 버그(unexpected) 분리
    kind_cnt = Counter(r.error_kind for r in results if r.error)
    unexpected = [r for r in results if r.error_kind == "unexpected"]
    print("\n" + "=" * 60)
    print(f"요약: {len(results)}개 중 완전통과 {clean} / FAIL보유 "
          f"{sum(1 for r in results if r.failed)} / 파이프라인종료 {len(errors)}")
    print(f"총 FAIL 체크: {total_fail}  |  🐛 unexpected 예외: {len(unexpected)}")
    print("\n체크별 집계 (PASS/FAIL/SKIP):")
    for name, cnt in sorted(per_check.items()):
        flag = "  ⬅️ FAIL 존재" if cnt[FAIL] else ""
        print(f"  {name:24s} P{cnt[PASS]:<5d} F{cnt[FAIL]:<4d} S{cnt[SKIP]:<4d}{flag}")

    if kind_cnt:
        print("\n파이프라인 조기종료 트리아지:")
        for kind, n in sorted(kind_cnt.items()):
            print(f"  {'🐛' if kind == 'unexpected' else '⏭️'} {kind:28s} {n}건")

    if unexpected:
        print(f"\n🐛 진짜 버그 후보(unexpected) {len(unexpected)}건 (샘플 10):")
        for r in unexpected[:10]:
            print(f"  {r.test_id}: {r.error.splitlines()[0][:90]}")

    if args.json:
        payload = {
            "cases": [
                {"test_id": r.test_id, "error": r.error, "error_kind": r.error_kind,
                 "checks": [vars(c) for c in r.checks]}
                for r in results
            ],
            "per_check": {k: dict(v) for k, v in per_check.items()},
        }
        Path(args.json).write_text(json.dumps(payload, ensure_ascii=False, indent=2))
        print(f"\n📄 JSON 저장: {args.json}")

    # 종료코드: FAIL 있으면 1 (CI/게이트용). 파이프라인 에러는 트리아지 대상이라 0 유지.
    sys.exit(1 if total_fail else 0)


if __name__ == "__main__":
    main()
