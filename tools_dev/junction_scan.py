"""접합부 설계기 전수조사 — 실제 파츠로 명세 생성·품질·주문길이·교차를 측정한다.

## 왜 필요한가

`params_cp_tests/g1_test.py`·`g2_test.py`는 **불변식**을 본다(병합 규칙, 운반자 규칙,
서열 오라클). 그런데 "실제 파츠에서 padding이 몇 %나 발동하는가", "올리고 주문 길이가
합성 상한을 넘는 조합이 있는가"는 불변식이 아니라 **분포**라서 전수로 봐야 한다.

## 데이터 출처

`archive/docs/analysis/Bioparts_DB/bioparts_data/*.csv` (iGEM 파츠). 운영 MySQL이 아니라
리포 안의 스냅샷을 쓴다 — 설계기는 서열만 받는 순수 함수라 DB 없이 돌 수 있어야 하고,
그래야 전수조사가 싸다(파츠 말단 전수 0.2초).

## 측정 비용 (실측)

| 축 | 단위 |
|---|---|
| 접합부 품질(GC+hairpin) | 9µs/창 → 파츠 11,120개 양끝 전수 0.2초 |
| 명세 생성 + 비인접 교차 | 9ms/모듈 |
| 프라이머 설계(지배 비용) | 22ms/벌 |

조합 전체(≈4.2조)는 불가능하므로 **파츠 말단은 전수, 모듈은 표본**으로 간다.

사용:
    python tools_dev/junction_scan.py              # 기본 표본
    python tools_dev/junction_scan.py --modules 3000 --primer-sample 300
"""
import argparse
import csv
import os
import random
import statistics
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from pydna.dseqrecord import Dseqrecord

from params import OLIGO_SYNTH_MAX_BP, DesignParams
from workflow.junction_plan import JUNCTION_UNIQUE_MAX_BP, JunctionPlanner

DATA_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "archive", "docs", "analysis", "Bioparts_DB", "bioparts_data",
)
FILES = {
    "promoter": "promoter_igem.csv",
    "rbs": "rbs_igem.csv",
    "cds": "cds_igem.csv",
    "terminator": "terminator_igem.csv",
}
MODULE_ORDER = ("promoter", "rbs", "cds", "terminator")


def load_parts():
    """타입별 (part_id, 서열) 목록. ATGC 외 염기는 걸러낸다."""
    out = {}
    for ptype, fname in FILES.items():
        rows = []
        path = os.path.join(DATA_DIR, fname)
        with open(path, encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh):
                raw = (r.get("Sequence") or "").upper()
                seq = "".join(c for c in raw if c in "ACGT")
                # 모호코드가 섞인 파츠는 조립 입력 자체가 거부되므로(E-1) 제외한다.
                if seq and len(seq) == len(raw.strip()):
                    rows.append((r["ID"], seq))
        out[ptype] = rows
    return out


def as_record(pid, ptype, seq):
    rec = Dseqrecord(seq, name=pid)
    rec.annotations["part_id"] = pid
    rec.annotations["part_type"] = ptype
    rec.annotations["molecule_type"] = "DNA"
    return rec


def scan_junction_quality(parts, prm):
    """파츠 말단 arm창 **전수** — padding이 얼마나 발동할지의 상한."""
    from workflow.junction_plan import JunctionPlanner as JP

    arm = prm.gibson_overlap
    planner = JP(params=prm)
    print(f"\n[A] 파츠 말단 창 전수 (최소 길이 {arm}bp 기준 — 범위 탐색 전 상한)")
    print(f"  {'타입':<12}{'검사창':>8}{'GC밖':>12}{'hairpin':>12}{'둘중하나':>12}")
    grand = bad_grand = 0
    for ptype in MODULE_ORDER:
        n = bad = bad_gc = bad_hp = 0
        for _pid, seq in parts[ptype]:
            if len(seq) < arm:
                continue
            for win in (seq[:arm], seq[-arm:]):
                n += 1
                ok, reason, _gc, _dg = planner._eval_window(win)
                if not ok:
                    bad += 1
                    if "GC" in reason:
                        bad_gc += 1
                    if "hairpin" in reason:
                        bad_hp += 1
        if not n:
            continue
        grand += n
        bad_grand += bad
        print(f"  {ptype:<12}{n:>8}{bad_gc:>7} ({100*bad_gc/n:>3.0f}%)"
              f"{bad_hp:>7} ({100*bad_hp/n:>3.0f}%){bad:>7} ({100*bad/n:>3.0f}%)")
    print(f"  {'합계':<12}{grand:>8}{'':>12}{'':>12}{bad_grand:>7} "
          f"({100*bad_grand/grand:>3.0f}%)")
    return bad_grand / grand if grand else 0.0


def scan_modules(parts, prm, n_modules, seed=20261024):
    """모듈 표본 — 명세 생성 성공률, 병합·padding·주문길이·교차 분포."""
    rng = random.Random(seed)
    planner = JunctionPlanner(params=prm)
    stats = {
        "ok": 0, "fail": 0, "frag": [], "oligo": [], "pcr": [], "merged": [],
"junctions": [], "order_len": [], "crosstalk": [],
        "over_synth": 0, "over_crosstalk": 0,
    }
    errors = Counter()
    t0 = time.perf_counter()
    for _ in range(n_modules):
        recs = []
        for ptype in MODULE_ORDER:
            pid, seq = rng.choice(parts[ptype])
            recs.append(as_record(pid, ptype, seq))
        try:
            # 조립기와 같은 **범위**를 넘긴다 — 기본값은 단일 길이라
            # 범위 탐색(= padding을 대신하는 유일한 탈출구)이 꺼진 채로 측정된다.
            plan = planner.plan(recs, arm_min=prm.gibson_overlap,
                                arm_max=prm.gibson_overlap_max)
        except Exception as exc:
            stats["fail"] += 1
            errors[f"{type(exc).__name__}: {str(exc)[:70]}"] += 1
            continue
        s = plan.summary()
        stats["ok"] += 1
        stats["frag"].append(s["fragments"])
        stats["oligo"].append(s["oligo"])
        stats["pcr"].append(s["pcr"])
        stats["merged"].append(s["absorbed"])
        stats["junctions"].append(s["junctions"])
        stats["crosstalk"].append(s["max_crosstalk"])
        for f in plan.fragments:
            stats["order_len"].append(f.order_len)
            if f.production == "oligo" and f.order_len > OLIGO_SYNTH_MAX_BP:
                stats["over_synth"] += 1
        if s["max_crosstalk"] >= JUNCTION_UNIQUE_MAX_BP:
            stats["over_crosstalk"] += 1
    dt = time.perf_counter() - t0

    print(f"\n[B] 모듈 표본 {n_modules}건 ({dt:.1f}초, {1000*dt/n_modules:.0f}ms/모듈)")
    print(f"  명세 생성 성공 {stats['ok']} / 실패 {stats['fail']}")
    for msg, cnt in errors.most_common(5):
        print(f"    실패 {cnt}건: {msg}")
    if stats["ok"]:
        def med(key):
            return statistics.median(stats[key])
        total_j = sum(stats["junctions"])
        print(f"  fragment 중앙 {med('frag'):.0f}개 "
              f"(올리고 {med('oligo'):.0f} / PCR {med('pcr'):.0f})")
        print(f"  tail 흡수가 일어난 모듈: "
              f"{100*sum(1 for m in stats['merged'] if m)/stats['ok']:.0f}%")
        print(f"  주문 길이: 중앙 {med('order_len'):.0f}nt / "
              f"최대 {max(stats['order_len'])}nt / 합성상한({OLIGO_SYNTH_MAX_BP}nt) 초과 "
              f"{stats['over_synth']}건")
        print(f"  비인접 교차: 중앙 {med('crosstalk'):.0f}bp / "
              f"최대 {max(stats['crosstalk'])}bp / 임계({JUNCTION_UNIQUE_MAX_BP}bp) 이상 "
              f"{stats['over_crosstalk']}모듈 "
              f"({100*stats['over_crosstalk']/stats['ok']:.0f}%)")
    return stats


def scan_bare_pcr(parts, prm, n_sample, seed=7):
    """긴 파츠 **맨몸 PCR** 설계 성공률 — `primer_len_max`를 바꿔가며 측정.

    네 알고리즘의 전제가 "긴 파츠는 tail 없이 증폭"이라, 이 성공률이 곧 경로의 가용성이다.
    `anchor_ends=True`는 프라이머 5' 시작을 파츠 경계에 고정하므로 자유도가 길이뿐이고,
    `primer_len` 범위가 좁으면 후보가 말라 Tm 하한에서 전멸한다.
    """
    from basic.primer_design import PrimerDesign

    rng = random.Random(seed)
    longs = [(pid, s) for ptype in ("promoter", "cds", "terminator")
             for pid, s in parts[ptype] if len(s) >= prm.oligo_max_bp]
    sample = rng.sample(longs, min(n_sample, len(longs)))
    print(f"\n[C] 맨몸 PCR 설계 성공률 (긴 파츠 {len(sample)}개 표본, "
          f"중앙 {statistics.median([len(s) for _i, s in sample]):.0f}bp)")
    print(f"  {'primer_len_max':>15}{'성공':>10}{'실패':>10}{'성공률':>9}{'소요':>9}")
    results = {}
    for lmax in (25, 30, 32, 35):
        p = DesignParams(primer_len_max=lmax,
                         primer_len_opt=min(prm.primer_len_opt, lmax))
        p.validate()
        ok = fail = 0
        t0 = time.perf_counter()
        for _pid, seq in sample:
            try:
                r = PrimerDesign.design_basic_primers(
                    template=seq, target_start=0, target_end=len(seq),
                    anchor_ends=True, check_specificity=False, params=p)
                if r.get("forward") and r.get("reverse"):
                    ok += 1
                else:
                    fail += 1
            except Exception:
                fail += 1
        dt = time.perf_counter() - t0
        results[lmax] = ok / (ok + fail) if (ok + fail) else 0.0
        print(f"  {lmax:>15}{ok:>10}{fail:>10}{100*results[lmax]:>8.0f}%{dt:>8.1f}s")
    return results


def scan_prep_e2e(parts, prm, n_modules, seed=913):
    """`prepare_parts` end-to-end 성공률 — shift 폴백까지 포함한 실제 가용성.

    [C]는 `design_basic_primers`를 날로 부른 값이고, 이쪽은 **실제 경로**
    (접합부 설계 → 병합 → 운반자 → padding → 프라이머/올리고)를 통과한 비율이다.
    병합이 짧은 부품을 흡수하고 shift가 Tm 사각지대를 메우므로 [C]보다 높아야 한다.
    """
    from assembly.factory import AssemblyFactory

    rng = random.Random(seed)
    gib = AssemblyFactory.create_assembler("gibson", params=prm)
    ok = fail = 0
    shifted_frags = pcr_frags = total_frags = 0
    shift_hist = Counter()
    errors = Counter()
    t0 = time.perf_counter()
    for _ in range(n_modules):
        recs = []
        for ptype in MODULE_ORDER:
            pid, seq = rng.choice(parts[ptype])
            recs.append(as_record(pid, ptype, seq))
        try:
            design, _plan = gib.prepare_parts(recs)
        except Exception as exc:
            fail += 1
            errors[f"{type(exc).__name__}: {str(exc)[:80]}"] += 1
            continue
        ok += 1
        for d in design:
            total_frags += 1
            if d["production"] == "pcr":
                pcr_frags += 1
                k = d["anneal_shift"][0]
                shift_hist[k] += 1
                if k:
                    shifted_frags += 1
    dt = time.perf_counter() - t0

    print(f"\n[D] prepare_parts end-to-end {n_modules}건 "
          f"({dt:.0f}초, {1000*dt/n_modules:.0f}ms/모듈)")
    print(f"  모듈 성공 {ok} / 실패 {fail}  → 성공률 {100*ok/(ok+fail):.0f}%")
    # ⚠️ 모듈 성공률 < fragment 성공률이다 — 모듈은 **모든** fragment가 성공해야 하므로
    #    fragment 성공률의 거듭제곱으로 떨어진다(PCR fragment 2개면 p²).
    for msg, cnt in errors.most_common(4):
        print(f"    실패 {cnt}건: {msg.splitlines()[0]}")
    if pcr_frags:
        print(f"  PCR fragment {pcr_frags}개(전체 {total_frags}개) 중 "
              f"shift 사용 {shifted_frags}개 ({100*shifted_frags/pcr_frags:.0f}%)")
        print("  shift 분포: " + ", ".join(
            f"k={k}:{v}" for k, v in sorted(shift_hist.items())))
    return ok / (ok + fail) if (ok + fail) else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modules", type=int, default=1000,
                    help="모듈 표본 수 (기본 1000)")
    ap.add_argument("--primer-sample", type=int, default=100,
                    help="맨몸 PCR 측정 표본 수 (기본 100). 0이면 생략")
    ap.add_argument("--arm", type=int, default=None, help="arm 길이 (기본 params 값)")
    ap.add_argument("--prep-modules", type=int, default=200,
                    help="prepare_parts end-to-end 표본 수 (기본 200). 0이면 생략")
    args = ap.parse_args()

    prm = DesignParams() if args.arm is None else DesignParams(gibson_overlap=args.arm)
    prm.validate()

    parts = load_parts()
    total = sum(len(v) for v in parts.values())
    print("=" * 72)
    print(f"접합부 설계기 전수조사 — 파츠 {total}개 "
          f"({' / '.join(f'{k} {len(v)}' for k, v in parts.items())})")
    print(f"arm={prm.gibson_overlap}bp  oligo_max_bp={prm.oligo_max_bp}  "
          f"합성상한={OLIGO_SYNTH_MAX_BP}nt  교차임계={JUNCTION_UNIQUE_MAX_BP}bp")
    print("=" * 72)

    scan_junction_quality(parts, prm)
    stats = scan_modules(parts, prm, args.modules)
    if args.primer_sample:
        scan_bare_pcr(parts, prm, args.primer_sample)
    if args.prep_modules:
        scan_prep_e2e(parts, prm, args.prep_modules)

    print("\n" + "=" * 72)
    # padding을 없앤 뒤로 **접합부 품질을 범위 안에서 못 맞추는 조합**은 실패가 정상이다
    #   (예전엔 padding으로 가렸고, 그 padding이 RBS–CDS 경계를 깨뜨렸다).
    #   그래서 "0건"이 아니라 **비율**로 본다. 합성상한 초과는 여전히 0건이어야 한다.
    FAIL_RATE_MAX = 0.05
    verdict = []
    total = stats["ok"] + stats["fail"]
    rate = stats["fail"] / total if total else 0.0
    if rate > FAIL_RATE_MAX:
        verdict.append(f"명세 생성 실패율 {100*rate:.1f}% > 허용 {100*FAIL_RATE_MAX:.0f}%")
    if stats["over_synth"]:
        verdict.append(f"합성상한 초과 {stats['over_synth']}건")
    print(f"명세 생성 실패율: {100*rate:.1f}% ({stats['fail']}/{total}) "
          f"— 범위로도 품질을 못 맞추는 조합")
    print("수용 기준 위반: " + (", ".join(verdict) if verdict else "없음"))
    return 1 if verdict else 0


if __name__ == "__main__":
    sys.exit(main())
