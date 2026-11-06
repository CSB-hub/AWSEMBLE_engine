"""Spacer 선정·부착 (P7 구멍1) — module을 vector에 삽입하기 위한 anneal 발판.

biopart module은 양 끝이 임의 부품 서열이라 vector 삽입 프라이머가 anneal할 발판이 없다.
양 끝에 spacer(vector 비종속 독립 서열)를 부착하면 프라이머가 spacer에 anneal해 삽입이 된다.

책임 분리(260614_P7_spacer_design.md §6.5):
  - 1차(spacer 고유 적합성: anneal 발판 + 조립효소 site 없음)는 **DB 자산에서 완료**
    (gen_spacers 제작규칙 + 165→115 정제). 이 모듈은 그 적합풀을 조회만 한다.
  - 2차(vector 조합 적합성)는 여기서 — vector가 정해질 때 후보 spacer를 두 관문으로 검증:
      ① 프라이머 품질  : 실제 삽입 프라이머를 validate_primer(tail_len)로 검증(hairpin/dimer 등)
      ② 실제 조립 dry-run: vector_insertion을 proxy로 돌려 삽입이 성사되는지
    둘 다 통과한 spacer만 채택 → "주문 가능한 프라이머" + "조립 성공" 모두 보장.

부착은 part_prep 정석: 선정한 spacer를 부품 리스트 양 끝에 끼워(prepare_fragments가 overhang
자동 할당) module = [L-spacer]-부품들-[R-spacer]로 조립한다.
"""
import contextlib
import io
from collections import Counter
from typing import List, Optional, Tuple

from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord

SCHEMA = "Bioparts"
SPACER_INSTITUTION = "Awsemble"
GG_ENZYME = "BsaI"   # vector insertion GG 기본 효소(내부 site 시 BbsI 전환은 prepare_fragments 책임)
DEGRADED_TOPN = 20   # graceful degradation 시 조립 dry-run할 후보 상위 N개(폭주 방지)

# 조립 dry-run 총 시도 상한. 체계적 실패(벡터 불호환 등)는 **모든 쌍이 같은 이유로** 실패하는데
# 그걸 구분할 방법이 없어 clean 후보가 많으면 115×114회까지 돈다(실측 14분+).
# 상한에 걸리면 **그 사실을 에러 메시지에 명시**한다 — 조용히 포기하지 않는다.
MAX_INSERTION_DRYRUNS = 400

# ✅ Golden Gate vector insertion: 2026-09-23 해결. 그 전까지는 GG fwd amplicon
#    (= L spacer 시작부)의 Tm이 판정 상한을 넘어 **L 후보 0개**로 선정이 늘 실패했다.
#    원인은 고GC spacer가 아니라 `tm_diff_max`였다 — tail이 고정 길이라 diff를 먼저
#    좁히면 amplicon이 길어지고 Tm이 상한을 넘었다. `_fit_amplicon`이 Tm을 선언 범위
#    [tm_min, tm_max]에 맞추고 diff를 결과로 보고하도록 바꿔 **115/115 통과**한다.
#    경위: develop_log/ljw/260923_engine_followup_plan.md (구 분석: 260616_P7_미해결_...)


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


class SpacerSelector:
    """vector 조합으로 적합 L/R spacer를 선정(2차)하고 part_prep용으로 부착한다.

    Args:
        db: DatabaseImport 인스턴스 (적합풀 조회 + spacer 서열 조회)
    """

    def __init__(self, db):
        self.db = db
        self._pool_ids: Optional[List[str]] = None
        # 직전 select()의 비차단 경고(graceful degradation 등). 호출자가 조회 가능.
        self.last_selection_warnings: List[str] = []

    def pool_ids(self) -> List[str]:
        """DB의 적합 spacer 풀 ID(오름차순). 1차 적합성은 자산에서 이미 보장됨."""
        if self._pool_ids is None:
            rows = self.db._execute(
                SCHEMA,
                "SELECT ID FROM spacer WHERE Institution=%s ORDER BY ID",
                (SPACER_INSTITUTION,),
            )
            self._pool_ids = [r[0] for r in rows]
        return self._pool_ids

    # ------------------------------------------------------------------ #
    # 2차 관문 ① 프라이머 품질 — side('fwd'=L 발판 / 'rev'=R 발판)별 단독 판정.
    #   fwd는 L spacer에, rev는 R spacer에만 의존하므로 L/R 후보를 독립 평가할 수 있다.
    # ------------------------------------------------------------------ #
    def _primer_side_eval(
        self, insert: Dseqrecord, vector: Dseqrecord, method: str, side: str, params=None
    ) -> Tuple[bool, float]:
        """실제 삽입 프라이머 한쪽(fwd|rev)의 (품질통과여부, 접힘점수)를 반환.

        접힘점수 = anneal 영역 hairpin ΔG(kcal/mol). 0에 가까울수록 덜 접힘(좋음).
        graceful degradation에서 clean 후보가 없을 때 '가장 덜 접히는' 후보를 고르는 데 쓴다.
        프라이머 설계 자체가 실패하면 (False, -inf)로 최하위 정렬.
        """
        from basic.primer_design import PrimerDesign, ValidatePrimer

        try:
            with contextlib.redirect_stdout(io.StringIO()):
                if method == "gibson":
                    from assembly.gibson import GibsonCloning
                    gib = GibsonCloning(params=params)
                    fwd, rev = gib.prepare_fragments(insert, vector)
                    primer = fwd if side == "fwd" else rev
                    # ② 조립 시 노출되는 말단 ss = vector homology(=5' tail, overlap_length bp).
                    #   직접 계산해 hairpin ΔG를 재도록 넘긴다(경고 판정용, passed엔 영향 없음).
                    # ⚠️ params를 넘겨야 **선정 기준과 실제 조립 기준이 같은 조건**이 된다.
                    #    안 넘기면 염/농도가 기본값으로 고정돼, 랩이 조건을 바꿔도 선정만 옛 조건으로 돈다.
                    v = ValidatePrimer.validate_primer(
                        primer,
                        tail_len=gib.overlap_length,
                        exposed_ss=primer[:gib.overlap_length],
                        params=params,
                    )
                    return v["passed"], v["hairpin_dg_anneal"]
                elif method == "golden_gate":
                    from assembly.golden_gate import GoldenGateModule
                    gg = GoldenGateModule(
                        params.gg_enzyme if params is not None else "BsaI", params=params)
                    v_left, v_right = gg.extract_vector_overhangs(vector)
                    # ⚠️ 효소도 params를 따라야 한다. 여기만 GG_ENZYME(BsaI) 고정이면
                    #    랩이 BsmBI를 골라도 **선정은 BsaI 프라이머로 판정**하고 조립은 BsmBI로 돈다.
                    design = PrimerDesign.design_goldengate_primers(
                        insert,
                        left_overhang=v_right,
                        right_overhang=_rc(v_left),
                        restriction_enzyme=(params.gg_enzyme if params is not None else GG_ENZYME),
                        params=params,
                    )
                    passed = design["fwd_passed"] if side == "fwd" else design["rev_passed"]
                    # GG 설계 결과에는 anneal ΔG가 없어 랭킹은 중립(0.0)으로 둔다.
                    # 통과 여부(passed)는 정상 판정된다 — 선정 자체는 동작한다.
                    return passed, 0.0
        except Exception:
            return False, float("-inf")
        return False, float("-inf")

    # ------------------------------------------------------------------ #
    # 2차 관문 ② 실제 조립 dry-run (선정 통과 == 실제 조립 성공)
    # ------------------------------------------------------------------ #
    def _insertion_ok(self, insert: Dseqrecord, vector: Dseqrecord, method: str, params=None) -> bool:
        """proxy insert를 실제 vector_insertion 경로로 조립해보고 성공하면 True.

        프라이머 품질(①)을 통과해도 sticky-end 불호환 등으로 조립이 안 될 수 있으므로,
        실제 조립까지 돌려 "삽입 성사"를 직접 확인한다(같은 코드 경로 = 선정==조립 일치).

        ⚠️ **실패 사유를 버리지 않는다** (2026-09-23): 구버전은 `except Exception: return False`로
        모든 예외를 삼켜서
          · 최종 에러가 "어떤 쌍도 실패"라는 일반 문구뿐이라 **원인이 표면에 안 나왔고**
            (GG 조립 버그 진단이 오래 걸린 직접 원인),
          · 오타·시그니처 불일치 같은 **진짜 버그가 "부적합 spacer"로 위장**됐다.
        이제 사유를 `self._insertion_failures`에 집계해 `select()`의 에러에 싣는다.
        """
        from workflow.classical_assemble import ClassicalWorkflow

        cw = ClassicalWorkflow(output_dir="/tmp/spacer_select")
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                results = cw.assemble(insert, vector, method=method, validate=False,
                                      design_params=params)
            if results:
                return True
            self._record_insertion_failure("the assembly produced nothing (no exception)")
            return False
        except Exception as e:
            self._record_insertion_failure(f"{type(e).__name__}: {e}")
            return False

    def _record_insertion_failure(self, reason: str) -> None:
        """조립 실패 사유를 집계. 메시지가 길면 앞부분만 서명으로 쓴다."""
        if not hasattr(self, "_insertion_failures"):
            self._insertion_failures = Counter()
        self._insertion_failures[reason[:160]] += 1

    def _failure_summary(self, limit: int = 3) -> str:
        """가장 흔한 실패 사유 몇 개를 사람이 읽을 형태로."""
        fails = getattr(self, "_insertion_failures", None)
        if not fails:
            return "(no failure reason recorded)"
        total = sum(fails.values())
        lines = [f"    - {cnt}/{total} times: {reason}" for reason, cnt in fails.most_common(limit)]
        if len(fails) > limit:
            lines.append(f"    - (and {len(fails) - limit} other kinds)")
        return "\n".join(lines)

    def select(
        self,
        core_seq: str,
        vector: Dseqrecord,
        method: str,
        params=None,
    ) -> Tuple[str, str]:
        """vector dry-run으로 L/R spacer 선정 (적합풀 순회, 두 관문 통과 첫 둘 = L,R).

        proxy module = [spacer + core_seq + spacer]로 평가한다. anneal 적합성은 spacer 끝
        ~20bp에만 의존하므로(P7 §6.5 실증) 조립 전 부품 concat(core_seq)으로 충분하다.

        Args:
            core_seq: 부품들을 이은 module 코어 서열(proxy). spacer 제외.
            vector: 삽입 대상 vector Dseqrecord
            method: 조립 방법 ("gibson" | "golden_gate")

        Returns:
            (L_spacer_id, R_spacer_id) — 서로 다른 두 적합 spacer.

        Raises:
            RuntimeError: 적합풀 자체가 2개 미만이거나, (degradation까지 해도) 어떤 쌍도 실제
                조립에 성사되지 않을 때(module/vector 조합 자체 문제 — 예: GG insert에 조립효소
                site). "clean 후보 부족"만으로는 raise하지 않는다(graceful degradation).

        Notes:
            proxy 양 끝은 **서로 다른** spacer여야 한다(같으면 프라이머가 양쪽에 anneal → PCR 비특이).
            fwd는 L spacer, rev는 R spacer에만 의존하므로 L/R 후보를 독립으로 평가한다.

            graceful degradation (2축 정책): hairpin은 "고칠 수 있는(spacer 교체 가능)" 축이라
            원칙은 clean 후보 채택이지만, **clean 후보가 0개여도 하드 실패하지 않고** 가장 덜
            접히는(anneal ΔG가 0에 가까운) 후보를 경고와 함께 사용한다("재설계해도 컷 못 넘기면
            경고 후 사용"). 실제 조립(②)이 성사되는 한 진행하며, 조립 자체가 전부 실패할 때만 raise.
        """
        # ── 불변 조건 먼저 (2026-10-20) ─────────────────────────────────────
        # spacer를 바꿔도 달라지지 않는 조건은 **여기서 끝낸다.** 예전엔 이 검사가 없어
        # 벡터를 열 수 없는 경우에도 spacer 조합 380회를 전부 돌고 나서야 같은 이유로
        # 실패했다(실측 6초, 380/380 동일 사유). spacer 탓이 아닌 실패를 spacer 탐색으로
        # 찾으려 한 셈이다.
        if method == "gibson":
            from basic.sequence_analysis import usable_linearization_enzymes
            want = getattr(params, "linearization_enzyme", None) if params else None
            usable = usable_linearization_enzymes(str(vector.seq),
                                                  circular=bool(getattr(vector, "circular", False)))
            if not usable:
                raise RuntimeError(
                    f"No enzyme cuts vector '{vector.name}' exactly once "
                    f"-> **this vector cannot take a Gibson insertion**. Changing the spacer "
                    f"will not help, so the search stops here. Use a vector with an MCS, or "
                    f"consider Golden Gate."
                )
            if want and want not in usable:
                raise RuntimeError(
                    f"The linearization enzyme you specified, {want}, does not cut vector "
                    f"'{vector.name}' exactly once. Enzymes usable on this vector: "
                    f"{', '.join(usable)}. This has nothing to do with the spacer, so the "
                    f"search stops here."
                )

        pool = self.pool_ids()
        if len(pool) < 2:
            raise RuntimeError(f"The eligible pool is too small ({len(pool)}). Check the spacer load in the DB.")

        seqs = {sid: str(self.db.get_part_by_id(SCHEMA, "spacer", sid).seq) for sid in pool}
        ref = seqs[pool[0]]   # 반대쪽 자리 임시 spacer (fwd/rev는 자기 쪽 spacer에만 의존)

        self.last_selection_warnings = []
        self._insertion_failures = Counter()   # 조립 실패 사유 집계(에러 메시지에 실음)

        # 각 후보를 side별로 평가: (품질통과, 접힘ΔG). sid가 L이면 fwd, R이면 rev.
        l_eval, r_eval = {}, {}
        for sid in pool:
            s = seqs[sid]
            l_eval[sid] = self._primer_side_eval(
                Dseqrecord(s + core_seq + ref, name=f"L_{sid}"), vector, method, "fwd", params=params)
            r_eval[sid] = self._primer_side_eval(
                Dseqrecord(ref + core_seq + s, name=f"R_{sid}"), vector, method, "rev", params=params)

        def ranked(evald):
            """clean(passed) 후보는 **pool 순서 그대로**(기존 선정 동작 보존). clean이 0개면
            (degradation) 전체를 덜 접히는(anneal ΔG 0에 가까운) 순 상위 DEGRADED_TOPN개만
            (조립 dry-run 폭주 방지)."""
            clean = [s for s in pool if evald[s][0]]
            if clean:
                return clean, False
            order = sorted(pool, key=lambda s: evald[s][1], reverse=True)  # ΔG 0에 가까운순
            return order[:DEGRADED_TOPN], True

        l_list, l_deg = ranked(l_eval)
        r_list, r_deg = ranked(r_eval)
        n_l_clean = sum(1 for s in pool if l_eval[s][0])
        n_r_clean = sum(1 for s in pool if r_eval[s][0])

        # ② 서로 다른 L,R 쌍을 실제 조립으로 확정 (프라이머 품질 통과 == 조립 성공 보장은 아니므로).
        #    순위 상위(덜 접히는) 쌍부터 시도해 성사되는 첫 쌍 채택.
        tried = 0
        capped = False
        for l_id in l_list:
            if capped:
                break
            for r_id in r_list:
                if r_id == l_id:
                    continue
                if tried >= MAX_INSERTION_DRYRUNS:
                    capped = True
                    break
                tried += 1
                proxy = Dseqrecord(seqs[l_id] + core_seq + seqs[r_id], name="proxy")
                if self._insertion_ok(proxy, vector, method, params=params):
                    if l_deg or r_deg:
                        w = (
                            f"Spacer selection fell back gracefully - not enough clean candidates"
                            f"(L clean={n_l_clean}, R clean={n_r_clean}/{len(pool)}). "
                            f"using the least-folded candidates: L={l_id}(dG {l_eval[l_id][1]:.2f}) "
                            f"R={r_id}(dG {r_eval[r_id][1]:.2f}). Assembly will work, but consider additives or a review."
                        )
                        self.last_selection_warnings.append(w)
                        print(w)
                    return l_id, r_id

        cap_note = (
            f" **Stopped at the attempt limit of {MAX_INSERTION_DRYRUNS}** (candidates remain)."
            if capped else ""
        )
        raise RuntimeError(
            f"No L/R pair assembled successfully (including the degraded fallback) - "
            f"vector={vector.name}, method={method}, {tried} attempts.{cap_note}\n"
            f"  Most common failure reasons:\n{self._failure_summary()}\n"
            f"  -> If every pair fails for the same reason, the problem is the "
            f"**module/vector combination itself**, not the spacer (for Golden Gate, check "
            f"the number of assembly-enzyme sites in the vector and sites inside the insert, "
            f"or switch to Gibson)."
        )

    def attach(
        self,
        parts: List[Dseqrecord],
        l_spacer_id,
        r_spacer_id,
    ) -> List[Dseqrecord]:
        """선정한 spacer를 부품 리스트 양 끝에 끼운다(part_prep 정석).

        반환 리스트를 prepare_fragments에 넣으면 spacer에도 overhang이 자동 할당되어
        module = [L-spacer]-부품들-[R-spacer]로 조립된다.

        Args:
            parts: 유저가 고른 부품 Dseqrecord 리스트(순서 의미 있음)
            l_spacer_id / r_spacer_id: spacer **ID**(str, 적합풀에서 고른 것) 또는
                이미 만들어진 **Dseqrecord**. `vector_homology` 모드의 양 끝 spacer는
                벡터 절단부에서 그 자리에서 만들어지므로 DB에 ID가 없어 record로 들어온다.

        Returns:
            [L_spacer, *parts, R_spacer] Dseqrecord 리스트
        """
        def _one(x) -> Dseqrecord:
            if not isinstance(x, Dseqrecord):
                return self.db.get_part_by_id(SCHEMA, "spacer", x)
            # **module마다 새 객체를 준다.** 멀티모듈에서 같은 경계 record를 두 module에
            # 그대로 넘기면 prepare_fragments가 annotations(overhang 등)를 덮어써서
            # 뒤 module의 설계가 앞 module 값을 물고 간다.
            rec = Dseqrecord(str(x.seq).upper(), name=x.name)
            rec.annotations.update(dict(x.annotations))
            return rec

        return [_one(l_spacer_id), *parts, _one(r_spacer_id)]


__all__ = ["SpacerSelector"]
