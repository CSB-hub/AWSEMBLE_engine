"""
Workflow Builder (통합 인터페이스)

워크플로우 생성을 위한 통합 인터페이스를 제공합니다.

사용자는 이 클래스를 통해 Biopart 또는 Classical 워크플로우를 생성합니다.

사용 예제:
    >>> # Classical 워크플로우
    >>> workflow = WorkflowBuilder.create_classical(output_dir="results")
    >>> workflow.run("gfp.fasta", "pUC19.gb", method="gibson")
    >>>
    >>> # Biopart 워크플로우 (parts + build plan 기반)
    >>> workflow = WorkflowBuilder.create_biopart(
    ...     project_dir="projects/my_project",
    ...     parts=prepared_records
    ... )
    >>> workflow.run_full_workflow(plan=build_plan)

작성일: 2025-10-30
"""

from typing import List, Optional

from pydna.dseqrecord import Dseqrecord

from .classical_assemble import ClassicalWorkflow


class _CrossAware(list):
    """설계 결과 리스트 + **쌍 단위** 교차 검사 결과를 같이 운반한다.

    `primers.csv`는 행이 부품 하나라 "부품 A ↔ 부품 C" 같은 쌍 결과를 담을 자리가
    없다. 그렇다고 `_prep_module_records`의 반환 서명을 바꾸면 호출부 3곳이 전부
    움직이고, 쌍 결과는 **run 단위**라 design 행에 넣는 것도 맞지 않다.
    → 리스트 자체에 붙인다. 호출부는 지금처럼 `records, design = ...`로 받고,
      리포트를 쓸 때만 `getattr(design, "cross_dimers", [])`로 꺼낸다.

    Attributes:
        cross_dimers: one-pot 프라이머 결합부끼리 dimer (Golden Gate 경로)
        crosstalk_pairs: 비인접 fragment 공유 서열 (Gibson / yeast HR 경로)
        scans: `[(scope, kind, 검사한 쌍 수), ...]` — **무엇을 몇 쌍 봤는지**.
            걸린 쌍이 0일 때 "검사했고 깨끗하다"를 말하려면 분모가 있어야 한다.
    """

    cross_dimers = ()
    crosstalk_pairs = ()
    scans = ()


class WorkflowBuilder:
    """
    워크플로우 생성 통합 인터페이스

    사용자는 이 클래스를 통해 워크플로우를 생성합니다.
    Factory 패턴을 사용하여 적절한 워크플로우 인스턴스를 반환합니다.

    Example:
        >>> # Classical 워크플로우 생성
        >>> workflow = WorkflowBuilder.create_classical(output_dir="results")
        >>>
        >>> # Biopart 워크플로우 생성 (parts + build plan 기반)
        >>> workflow = WorkflowBuilder.create_biopart(
        ...     project_dir="projects/my_project",
        ...     parts=prepared_records
        ... )

    작성일: 2025-10-30
    """

    @staticmethod
    def create_classical(output_dir: str = "output") -> ClassicalWorkflow:
        """
        Classical 워크플로우 생성

        간단한 클로닝 워크플로우를 생성합니다.
        Insert + Vector → Assembled Plasmid

        Args:
            output_dir: 출력 디렉토리 (기본값: "output")

        Returns:
            ClassicalWorkflow 인스턴스

        Example:
            >>> workflow = WorkflowBuilder.create_classical(output_dir="results")
            >>> workflow.run(
            ...     insert_file="gfp.fasta",
            ...     vector_file="plasmid.gb",
            ...     method="gibson"
            ... )

        Notes:
            - GFP 삽입 같은 간단한 클로닝에 적합
            - GUI 없이 script로 사용
            - 빠른 테스트 및 검증용
        """
        return ClassicalWorkflow(output_dir=output_dir)

    @staticmethod
    def create_biopart(
        project_dir: str,
        parts: List[Dseqrecord],
        vector_file: Optional[str] = None
    ):
        """Biopart 워크플로우 생성 (K-2: parts + build plan 기반).

        구 Excel 설계 파일(design_file) 시그니처를 제거하고, 이미 준비된
        prepared part(Dseqrecord) 리스트를 직접 주입한다. 조립 순서/방법 같은
        조합 정의는 고정 step이 아니라 run_full_workflow(plan)의 **build plan**
        으로 표현한다(입력원 무관 — CLI/galaxy가 build plan 주입).

        Args:
            project_dir: 프로젝트 디렉토리
            parts: prepared part Dseqrecord 리스트 (DB 조회 + prepare_fragments 산출)
            vector_file: 벡터 GenBank 파일명 (선택, vector_insertion stage는 P7)

        Returns:
            BiopartWorkflow 인스턴스

        Example:
            >>> workflow = WorkflowBuilder.create_biopart(
            ...     project_dir="projects/my_project",
            ...     parts=prepared_records,
            ... )
            >>> workflow.run_full_workflow(plan=[
            ...     {"op": "module_assembly", "inputs": [r.name for r in prepared_records],
            ...      "params": {"method": "golden_gate"}},
            ... ])
        """
        from .biopart_assemble import BiopartWorkflow
        return BiopartWorkflow(project_dir=project_dir, parts=parts, vector_file=vector_file)

    @staticmethod
    def _save_parts(project_dir: str, parts_sp: List[Dseqrecord], records: List[Dseqrecord]) -> None:
        """parts/raw(spacer 부착 원본 부품) + parts/gg_prepared(overhang 붙은 prepared) 저장."""
        from pathlib import Path
        from utils.io import DNAFileWriter

        base = Path(project_dir)
        for sub, recs in (("parts/raw", parts_sp), ("parts/gg_prepared", records)):
            d = base / sub
            d.mkdir(parents=True, exist_ok=True)
            for r in recs:
                DNAFileWriter.write_genbank(r, str(d / f"{r.name}.gb"))

    @staticmethod
    def _save_gibson_prepared(project_dir: str, insert: Dseqrecord, vector: Dseqrecord, level: str,
                              params=None) -> None:
        """vector 삽입(gibson) 직전 insert PCR amplicon → {level}/gibson_prepared,
        삽입 프라이머(fwd/rev) → reports/vector_primers.csv. (저장 전용 1회 계산)

        Args:
            level: "modules"(single module 삽입) | "multimodules"(multimodule 삽입)
        """
        import csv as _csv
        from pathlib import Path
        from utils.io import DNAFileWriter
        from assembly.gibson import GibsonCloning

        from .vector_homology_spacer import module_carries_homology

        gib = GibsonCloning(params=params)

        # module이 벡터 homology를 이미 들고 있으면 **삽입 프라이머가 존재하지 않는다.**
        # 리포트에도 그 사실을 적는다 — 빈 칸을 남기면 "설계가 실패했나"로 읽힌다.
        _lin = gib._linearize_vector(vector)
        if module_carries_homology(insert, _lin, gib.overlap_length):
            base = Path(project_dir)
            gp = base / level / "gibson_prepared"
            gp.mkdir(parents=True, exist_ok=True)
            # PCR을 안 하므로 module 자체가 그대로 삽입 조각이다.
            DNAFileWriter.write_genbank(insert, str(gp / f"{insert.name}.gb"))
            rp = base / "reports"
            rp.mkdir(parents=True, exist_ok=True)
            with open(rp / "vector_primers.csv", "w", newline="", encoding="utf-8") as f:
                w = _csv.writer(f)
                w.writerow(["insert", "method", "forward", "reverse",
                            "fwd_mer", "rev_mer", "fwd_tm", "rev_tm", "ta",
                            "tm_diff", "amplicon_len", "note"])
                w.writerow([insert.name, "gibson", "", "", 0, 0, "", "", "", "",
                            len(insert.seq),
                            "No insertion primers needed - both module ends already carry vector homology "
                            "(spacer_mode=vector_homology). There are no insertion primers to order."])
            return

        fwd, rev = gib.prepare_fragments(insert, vector)
        pcr = gib.amplify_insert(insert, fwd, rev)
        pcr.name = insert.name

        base = Path(project_dir)
        gp = base / level / "gibson_prepared"
        gp.mkdir(parents=True, exist_ok=True)
        DNAFileWriter.write_genbank(pcr, str(gp / f"{insert.name}.gb"))

        rp = base / "reports"
        rp.mkdir(parents=True, exist_ok=True)
        # 2026-11-05(Y10-b) — 주문 시트에 필요한 건 서열만이 아니다.
        #   `fwd_mer`/`rev_mer`: 주문 단위(몇-mer). 서열을 눈으로 세게 하면 오발주가 난다.
        #   `fwd_tm`/`rev_tm`  : **anneal 기준**(tail 제외) — 설계기가 계산한 값을 그대로 쓴다.
        #   `ta`               : 권장 annealing 온도. 제품 규칙(`NEB_TA_RULES`)이 정한다.
        #                        프리셋 조성과 다르면 추측하지 않고 비운다.
        #   `tm_diff`          : fwd·rev Tm 차이(Y11). `ta`가 **낮은 쪽 기준**이라
        #                        이 값이 크면 뜨거운 쪽이 자기 Tm보다 낮은 온도에서 돈다.
        _dsg = getattr(gib, "last_primer_design", None) or {}
        _tf, _tr = _dsg.get("fwd_tm"), _dsg.get("rev_tm")
        _prm = gib.params
        _ta = None
        if isinstance(_tf, (int, float)) and isinstance(_tr, (int, float)):
            _ta = _prm.ta_recommended(_tf, _tr)
        _td = abs(_tf - _tr) if isinstance(_tf, (int, float)) and isinstance(_tr, (int, float)) else None

        def _r(x, nd=2):
            return round(x, nd) if isinstance(x, (int, float)) else ""

        with open(rp / "vector_primers.csv", "w", newline="", encoding="utf-8") as f:
            w = _csv.writer(f)
            w.writerow(["insert", "method", "forward", "reverse",
                        "fwd_mer", "rev_mer", "fwd_tm", "rev_tm", "ta",
                        "tm_diff", "amplicon_len"])
            w.writerow([insert.name, "gibson", fwd, rev, len(fwd), len(rev),
                        _r(_tf), _r(_tr), _r(_ta, 1), _r(_td, 1), len(pcr.seq)])

    @staticmethod
    def _insert_arm(insert_method: str, params) -> int:
        """삽입 조립에 쓸 homology arm 길이 — **삽입 방법의 축**을 따른다.

        조립 방법마다 arm 필드가 다르다(`GibsonCloning.ARM_PARAMS`,
        `YeastHRCloning.ARM_PARAMS`). module 양 끝에 미리 박아 넣는
        homology(`VectorHomologySpacer`)는 **나중에 삽입 조립의 overlap이 되므로**,
        그 길이를 삽입 방법에서 읽어와야 한다.

        Args:
            insert_method: 삽입 조립 방법 ("gibson" | "golden_gate" | "yeast_hr")
            params: 설계 파라미터 한 벌

        Returns:
            arm 길이(bp). 방법이 homology를 안 쓰거나(Golden Gate) 알 수 없으면
            `params.gibson_overlap`로 떨어진다.
        """
        from assembly.factory import AssemblyFactory

        try:
            assembler = AssemblyFactory.create_assembler(insert_method, params=params)
        except Exception:
            return int(params.gibson_overlap)
        return int(getattr(assembler, "overlap_length", params.gibson_overlap))

    @staticmethod
    def _prep_module_records(parts_sp: List[Dseqrecord], params=None, vector=None,
                             assembly_type: str = "linear",
                             module_method: str = "golden_gate"):
        """spacer 부착된 부품 리스트 → prep → (prepared records, design).

        build_plasmid/build_multimodule_plasmid 공통. `module_method`에 따라 prep이 갈린다:

        - `"golden_gate"` — `GoldenGateModule.prepare_fragments`. 부품마다 양 끝에
          인식서열+spacer+overhang을 붙인다. **부품 개수 == fragment 개수**.
        - `"gibson"` / `"yeast_hr"` — `GibsonCloning.prepare_parts`(G3). 접합부가
          인접 두 부품이 **공유하는 구간**이라 부품을 하나씩 볼 수 없다. arm보다 짧은
          부품은 이웃에 병합되므로 **부품 개수 != fragment 개수**가 될 수 있다.

        design 원본도 함께 반환 → primers.csv(write_primer_csv)에 사용.

        Args:
            parts_sp: 양 끝에 spacer가 부착된 부품 Dseqrecord 리스트
            params: 설계 파라미터
            vector: GG circular(중간 클로닝)에서 backbone으로 참여할 벡터
            assembly_type: "linear" | "circular"
            module_method: module 조립 방법

        Returns:
            (records, design) — records: module_assembly 입력 / design: prep 원본
        """
        from pydna.dseqrecord import Dseqrecord as _Dseqrecord
        from assembly.golden_gate import GoldenGateModule

        from params import DEFAULT_PARAMS
        prm = params if params is not None else DEFAULT_PARAMS
        # 하드 리밋 위반은 여기서 즉시 막는다. 안 막으면 primer3 안쪽에서 원인 불명 메시지로
        # 터지거나(길이/열역학 상한), module_assembly_type 오타처럼 **조용히 다른 경로로** 샌다.
        # 소프트 경고는 차단하지 않고 reports/design_conditions.md에 남긴다(base.py).
        prm.validate()
        if module_method.lower() in ("gibson", "yeast_hr"):
            # homology 경로(G3): 접합부 설계기가 병합·운반자·padding을 정한다.
            # GG와 달리 레코드를 여기서 재조립하지 않는다 — `JunctionPlan.records()`가
            # part_id/part_type/병합 내역까지 실어 주므로 그대로 쓴다.
            from assembly.factory import AssemblyFactory

            assembler = AssemblyFactory.create_assembler(module_method, params=prm)
            design, plan = assembler.prepare_parts(parts_sp)
            for w in plan.warnings:
                print(f"  ⚠️  {w}")
            # 기대 서열(오라클)을 design에 실어 보낸다 — 조립 결과 대조에 쓸 수 있다.
            for d in design:
                d["module_seq_expected"] = plan.module_seq
            # 비인접 fragment 공유 서열(쌍 단위) — reports/cross_checks.csv로 나간다.
            #   yeast는 `YeastHRCloning._check_plan`이 이미 하드 차단했으므로 여기까지
            #   오면 Gibson 트랙이다. 쌍 결과는 run 단위라 design 행에 담을 자리가
            #   없어서 **리스트 객체에 속성으로** 붙인다(호출부 3곳의 서명을 안 바꾼다).
            design = _CrossAware(design)
            design.crosstalk_pairs = list(plan.crosstalk_pairs)
            design.scans = [("", "junction_crosstalk",
                             getattr(plan, "crosstalk_checked", 0))]
            return plan.records(), design

        # 효소와 설계 조건을 그대로 넘긴다(안 넘기면 기본값으로 돌아 유저 입력이 무시된다).
        gg = GoldenGateModule(prm.gg_enzyme, params=prm)
        design = gg.prepare_fragments(parts_sp, vector=vector, assembly_type=assembly_type)
        records = []
        for res in design:
            r = _Dseqrecord(res["full_sequence"], name=res["part_name"])
            r.annotations["target_seq"] = res.get("target_seq", "")
            r.annotations["restriction_enzyme"] = res.get("restriction_enzyme", "BsaI")
            # prepare_fragments가 보존해둔 원본 부품 id/type을 prepared record로 넘긴다.
            # (여기서 새 Dseqrecord를 만들면서 안 넘기면 원본 annotations가 통째로 끊겨,
            #  module feature 기입(_assemble_fragments→annotate_module_parts)에서
            #  part_type이 빈 값이 돼 전부 misc_feature로 떨어진다.)
            r.annotations["part_id"] = res.get("part_id", "")
            r.annotations["part_type"] = res.get("part_type", "")
            records.append(r)
        design = _CrossAware(design)
        design.cross_dimers = list(getattr(gg, "last_cross_dimers", []) or [])
        design.scans = [("", "primer_dimer",
                         getattr(gg, "last_cross_dimer_pairs", 0))]
        return records, design

    @staticmethod
    def _write_cross_checks(wf, design, params=None) -> int:
        """`reports/cross_checks.csv`를 쓴다 — **검사가 돌았으면** 결과가 0이어도 쓴다.

        쌍 단위 결과는 `_CrossAware` 리스트에 속성으로 실려 온다(그 클래스 독스트링 참고).

        ⚠️ 2026-11-06까지는 걸린 쌍이 0이면 파일을 안 만들었다. 의도는 "빈 파일은
        깨끗함과 미검사가 구분 안 된다"였는데, Galaxy 출력은 파일이 있어야 해서
        **헤더만 있는 파일**이 유저에게 나갔고 결국 같은 모호함이 생겼다. 지금은
        반대로 간다 — 검사가 돌았으면 `..._summary` 행에 "N쌍 중 0쌍"을 적어 파일
        자체가 답하게 한다. 검사가 아예 없었을 때만 파일이 없다.

        Args:
            wf: WorkflowBase 인스턴스(출력 경로를 안다)
            design: `_prep_module_records`가 돌려준 설계 결과 리스트
            params: 요약 행에 적을 기준값(ΔG 경계 등)

        Returns:
            쓴 행 수 (0이면 파일 없음)
        """
        xd = list(getattr(design, "cross_dimers", []) or [])
        xt = list(getattr(design, "crosstalk_pairs", []) or [])
        sc = list(getattr(design, "scans", []) or [])
        # 단일 module 경로는 반응이 하나다 — scope를 비워두면 "범위 불명"으로 읽힌다.
        xd = [h if h.get("scope") else {**h, "scope": "module 1/1"} for h in xd]
        xt = [t if len(t) > 4 else tuple(t) + ("module 1/1",) for t in xt]
        sc = [s if s[0] else ("module 1/1",) + tuple(s[1:]) for s in sc]
        if not xd and not xt and not sc:
            return 0
        out = wf.project_dir / "reports" / "cross_checks.csv"
        cnt = wf.write_cross_checks_csv(out, cross_dimers=xd, crosstalk_pairs=xt,
                                        scans=sc, params=params)
        if cnt:
            print(f"  Cross-checks: {cnt} rows -> reports/cross_checks.csv")
        return cnt

    @staticmethod
    def build_plasmid(
        parts: List[Dseqrecord],
        vector: Dseqrecord,
        db,
        project_dir: str,
        module_method: str = "golden_gate",
        insert_method: str = "gibson",
        params=None,
        intermediate_vector=None,
    ) -> dict:
        """P7 정식 플라스미드 설계 파이프라인 (parts + vector → construct + report).

        유저가 고른 부품(parts)과 vector를 받아 한 흐름으로 처리한다:
          1. spacer 선정·부착(part_prep): vector 조합 2차 선정(SpacerSelector) → 양 끝 spacer 부착
          2. prepare_fragments: spacer 포함 부품에 GG overhang 자동 할당
          3. run_plan: module_assembly(GG) → vector_insertion(insert_method) 정식 경로
          4. save_report: reports/run_report.md (manifest 계보)

        Args:
            parts: 유저가 고른 부품 Dseqrecord 리스트 (DB 조회 산물, 순서 의미)
            vector: 삽입 대상 vector Dseqrecord
            db: DatabaseImport (spacer 적합풀/서열 조회용)
            project_dir: 프로젝트 디렉토리
            module_method: module 조립 방법 ("golden_gate" | "gibson" | "yeast_hr")
            insert_method: vector 삽입 방법 ("gibson" | "golden_gate")

        Returns:
            {"workflow", "results"(run_plan 반환), "l_spacer", "r_spacer", "construct"}

        Raises:
            RuntimeError: 적합 spacer 부족(SpacerSelector.select) 등.

        Notes:
            module 조립은 golden_gate가 기본이고, **gibson/yeast_hr 경로가 2026-10-24에
            열렸다**(G1~G3). 전용 prep은 `GibsonCloning.prepare_parts` +
            `workflow.junction_plan`이며, GG와 달리 arm보다 짧은 부품이 이웃에 병합되므로
            **부품 개수와 fragment 개수가 다를 수 있다.**
            계획·실측: develop_log/ljw/261024_gibson_part_module_and_yeast_plan.md
            GG vector 삽입은 2026-09-23에 해결됐다(구 §B의 amplicon Tm 이슈).
            원인은 `tm_diff_max`였고, `_fit_amplicon`이 Tm을 선언 범위에 맞추도록 바꿔
            spacer 115개 전부 통과한다.
        """
        from .spacer_select import SpacerSelector
        from .identity import derived_id
        from .biopart_assemble import BiopartWorkflow

        from params import DEFAULT_PARAMS
        prm = params if params is not None else DEFAULT_PARAMS
        # 하드 리밋 위반은 여기서 즉시 막는다. 안 막으면 primer3 안쪽에서 원인 불명 메시지로
        # 터지거나(길이/열역학 상한), module_assembly_type 오타처럼 **조용히 다른 경로로** 샌다.
        # 소프트 경고는 차단하지 않고 reports/design_conditions.md에 남긴다(base.py).
        prm.validate()

        # ── 조립 효소를 **한 번만** 결정한다 ────────────────────────────────
        # Golden Gate는 one-pot이라 벡터와 모든 부품이 같은 효소로 잘려야 하는데,
        # 두 쪽 조건이 정반대다(벡터는 사이트가 있어야 하고, 부품은 없어야 한다).
        # 예전엔 이 결정이 세 군데에 흩어져 서로 모르는 채 각자 정했다 — 자세한 배경은
        # `assembly.golden_gate.select_gg_enzyme` 독스트링.
        #
        # ⚠️ **벡터 조건은 삽입이 Golden Gate일 때만 건다.** Module based(삽입=Gibson)는
        #    벡터를 효소로 열지 않으므로 벡터에 Type IIS 사이트가 없어도 된다. 여기에
        #    벡터를 넘기면 멀쩡한 메인 경로가 효소를 바꾸게 되고 artifact id까지 흔들린다.
        if "golden_gate" in (module_method, insert_method):
            from assembly.golden_gate import select_gg_enzyme
            pick = select_gg_enzyme(
                parts,
                vector=vector if insert_method == "golden_gate" else None,
                preferred=prm.gg_enzyme,
            )
            if pick['enzyme'] is None:
                raise RuntimeError(pick['reason'])
            if pick['switched']:
                print(f"     ⚠️ {pick['reason']}")
                prm = prm.with_(gg_enzyme=pick['enzyme'])

        if prm.module_assembly_type == "circular":
            # 중간 클로닝: module을 중간 벡터 안에서 circular로 조립했다가 꺼내 최종 벡터에 넣는다.
            # spacer는 부품이 아니라 **중간 벡터에 박혀 있으므로** SpacerSelector를 쓰지 않는다.
            if intermediate_vector is None:
                raise ValueError(
                    "The circular (intermediate cloning) mode needs an intermediate_vector. "
                    "Specify a module Golden Gate vector that carries the spacers "
                    "(tools_dev/make_intermediate_vector.py can create one)."
                )
            return WorkflowBuilder._build_plasmid_circular(
                parts, vector, intermediate_vector, db, project_dir, insert_method, prm
            )

        # 1) spacer 선정·부착 (part_prep 정석)
        selector = SpacerSelector(db)
        core_seq = "".join(str(p.seq) for p in parts)
        # ⚠️ spacer 선정은 프라이머 품질로 후보를 거르므로 같은 파라미터를 넘겨야
        #    "선정 기준"과 "실제 조립 기준"이 어긋나지 않는다.
        if getattr(prm, "spacer_mode", "library") == "vector_homology":
            # 벡터 절단부 서열을 그대로 양 끝에 붙인다 → 조합 탐색도, 삽입 PCR도 없다.
            # ⚠️ arm은 **삽입 방법의 축**을 따라야 한다 — 이 spacer가 곧 삽입 조립의
            #    overlap이 되므로, 길이가 어긋나면 조립이 조용히 실패한다(§_insert_arm).
            from .vector_homology_spacer import VectorHomologySpacer
            parts_sp = VectorHomologySpacer(
                params=prm, arm=WorkflowBuilder._insert_arm(insert_method, prm),
            ).attach(parts, vector)
            l_spacer, r_spacer = parts_sp[0].name, parts_sp[-1].name
        else:
            l_spacer, r_spacer = selector.select(core_seq, vector, insert_method, params=prm)
            parts_sp = selector.attach(parts, l_spacer, r_spacer)

        # 2) prepare_fragments → module_assembly용 prepared records (공통 헬퍼)
        records, design = WorkflowBuilder._prep_module_records(
            parts_sp, params=prm, module_method=module_method)
        WorkflowBuilder._save_parts(project_dir, parts_sp, records)   # parts/raw + parts/gg_prepared

        # 3) run_plan 정식 경로: module_assembly → vector_insertion (공유 stage)
        wf = BiopartWorkflow(project_dir=project_dir, parts=records + [vector])
        # run_plan의 module_assembly stage가 만드는 id와 **같은 재료**로 계산해야 plan의
        # inputs 참조가 맞는다(어긋나면 KeyError: build plan 입력 artifact를 찾을 수 없음).
        _pk = prm.canonical()
        module_id = derived_id("module", [r.name for r in records] + [module_method, "linear", _pk])
        plan = [
            {"op": "module_assembly", "inputs": [r.name for r in records],
             "params": {"method": module_method, "assembly_type": "linear"}},
            {"op": "vector_insertion", "inputs": [module_id, vector.name],
             "params": {"method": insert_method, "assembly_type": "circular"}},
        ]
        results = wf.run_plan(plan, design_params=prm)
        if insert_method == "gibson" and results:
            # modules/gibson_prepared(vector 삽입 직전 amplicon) + reports/vector_primers.csv
            WorkflowBuilder._save_gibson_prepared(project_dir, results[0]["record"], vector, "modules",
                                                  params=prm)
        wf.write_primer_csv(design)   # module 조립 프라이머 → reports/primers.csv (주문용)
        WorkflowBuilder._write_cross_checks(wf, design, params=prm)
        wf.save_report()

        return {
            "workflow": wf,
            "results": results,
            "l_spacer": l_spacer,
            "r_spacer": r_spacer,
            "construct": results[-1]["record"] if results else None,
        }

    @staticmethod
    def _build_plasmid_circular(parts, vector, intermediate_vector, db, project_dir,
                                insert_method, prm) -> dict:
        """중간 클로닝 경로: 부품 → (중간 벡터 안 circular module) → 꺼냄 → 최종 벡터 삽입.

        linear 경로와의 차이:
          - spacer를 부품 리스트에 끼우지 않는다(중간 벡터가 제공).
          - module_assembly가 circular로 돌고 중간 벡터가 조립에 참여한다.
          - module_excision stage가 하나 더 붙는다(digest 기본, 불가 시 PCR 폴백).
        """
        from .biopart_assemble import BiopartWorkflow
        from .identity import derived_id
        from .spacer_select import SpacerSelector

        # 0) spacer 처리 — 벡터가 spacer를 갖고 있지 않으면 **부품 양 끝에 붙인다**.
        #    이렇게 하면 벡터 하나로 모든 module 위치를 커버할 수 있다(위치별 전용 벡터 불필요).
        #    꺼낼 때 spacer에 프라이머를 앉히면 양 끝이 spacer가 되어 인접 module과 구분되는
        #    Gibson homology가 된다. (digest로 꺼내면 벡터 서열이 딸려 나와 모든 module이 같아진다)
        spacer_recs = []
        if prm.intermediate_vector_has_spacers:
            parts_sp = list(parts)
            l_spacer = r_spacer = None
        else:
            selector = SpacerSelector(db)
            core_seq = "".join(str(p.seq) for p in parts)
            if getattr(prm, "spacer_mode", "library") == "vector_homology":
                from .vector_homology_spacer import VectorHomologySpacer
                # arm은 삽입 방법의 축을 따른다 — 위 build_plasmid와 같은 이유.
                parts_sp = VectorHomologySpacer(
                    params=prm, arm=WorkflowBuilder._insert_arm(insert_method, prm),
                ).attach(parts, vector)
                l_spacer, r_spacer = parts_sp[0].name, parts_sp[-1].name
            else:
                l_spacer, r_spacer = selector.select(core_seq, vector, insert_method, params=prm)
                parts_sp = selector.attach(parts, l_spacer, r_spacer)
            spacer_recs = [parts_sp[0], parts_sp[-1]]   # excision PCR 앵커

        # 1) prepare_fragments(circular) — 중간 벡터의 overhang이 체인 양 끝을 담당
        records, design = WorkflowBuilder._prep_module_records(
            parts_sp, params=prm, vector=intermediate_vector, assembly_type="circular"
        )
        WorkflowBuilder._save_parts(project_dir, parts_sp, records)

        # 2) run_plan: module_assembly(circular) → module_excision → vector_insertion
        wf = BiopartWorkflow(
            project_dir=project_dir,
            parts=records + [intermediate_vector, vector] + spacer_recs,
        )
        pk = prm.canonical()
        module_id = derived_id(
            "module",
            [r.name for r in records] + [intermediate_vector.name, "golden_gate", "circular", pk],
        )
        excised_id = derived_id("module", [module_id, "excision", prm.excision_method, pk])
        plan = [
            {"op": "module_assembly",
             "inputs": [r.name for r in records] + [intermediate_vector.name],
             "params": {"method": "golden_gate", "assembly_type": "circular"}},
            {"op": "module_excision",
             # spacer를 부품으로 붙였으면 그 spacer가 PCR 앵커가 된다(양 끝이 spacer가 되도록).
             "inputs": ([module_id, intermediate_vector.name]
                        + [r.name for r in spacer_recs]),
             "params": {"method": prm.excision_method, "assembly_type": "linear"}},
            {"op": "vector_insertion",
             "inputs": [excised_id, vector.name],
             "params": {"method": insert_method, "assembly_type": "circular"}},
        ]
        results = wf.run_plan(plan, design_params=prm)
        if insert_method == "gibson" and len(results) >= 2:
            WorkflowBuilder._save_gibson_prepared(
                project_dir, results[-2]["record"], vector, "modules", params=prm
            )
        wf.write_primer_csv(design)
        WorkflowBuilder._write_cross_checks(wf, design, params=prm)
        wf.save_report()

        return {
            "workflow": wf,
            "results": results,
            "intermediate_vector": intermediate_vector.name,
            "l_spacer": l_spacer,
            "r_spacer": r_spacer,
            "construct": results[-1]["record"] if results else None,
        }

    @staticmethod
    def build_multimodule_plasmid(
        module_groups: List[List[Dseqrecord]],
        vector: Dseqrecord,
        db,
        project_dir: str,
        insert_method: str = "gibson",
        params=None,
        intermediate_vector=None,
        module_method: str = "golden_gate",
    ) -> dict:
        """P7 multi-module 플라스미드 설계 (module 그룹들 → 공유 spacer 연결 → vector 삽입).

        KRIBB 패턴: 인접 module이 공유 spacer로 이어진다(module_i의 R == module_{i+1}의 L).
        그 공유 spacer(57bp)가 곧 multimodule_assembly(gibson)의 homology가 되어 연결된다.
          - 양 끝 spacer: vector 접점 → SpacerSelector.select
          - 연결 spacer: 적합풀에서 인접마다 다른 것 배정(15-mer 직교로 mis-assembly 방지)

        Args:
            module_groups: 각 module의 부품 Dseqrecord 리스트의 리스트 (2개 이상).
                기존 spacer가 있는 부품은 제외하고 core 부품(promoter/rbs/cds/terminator)부터 준다.
            vector: 삽입 대상 vector Dseqrecord
            db: DatabaseImport (spacer 적합풀/서열 조회)
            project_dir: 프로젝트 디렉토리
            insert_method: vector 삽입 + module 간 연결 방법 ("gibson" | "golden_gate")

        Returns:
            {"workflow","results","boundaries"([L,*connectors,R]),"construct"}

        Raises:
            ValueError: module_groups < 2 (단일은 build_plasmid).
            RuntimeError: 적합 spacer/연결 spacer 부족.
        """
        from .spacer_select import SpacerSelector
        from .identity import derived_id
        from .biopart_assemble import BiopartWorkflow

        n = len(module_groups)
        if n < 2:
            raise ValueError("module_groups needs two or more entries (use build_plasmid for a single module).")

        # 1) 양 끝 spacer(vector 접점) + 연결 spacer(인접 공유) 배정
        from params import DEFAULT_PARAMS
        prm = params if params is not None else DEFAULT_PARAMS
        # 하드 리밋 위반은 여기서 즉시 막는다. 안 막으면 primer3 안쪽에서 원인 불명 메시지로
        # 터지거나(길이/열역학 상한), module_assembly_type 오타처럼 **조용히 다른 경로로** 샌다.
        # 소프트 경고는 차단하지 않고 reports/design_conditions.md에 남긴다(base.py).
        prm.validate()

        selector = SpacerSelector(db)
        core_seq = "".join(str(p.seq) for g in module_groups for p in g)
        # ── 양 끝(벡터 접점)은 spacer_mode를 따른다 ────────────────────────────
        # ⚠️ 2026-11-04까지 이 경로는 `prm.spacer_mode`를 **아예 보지 않았다.** GUI는
        #    "벡터 맞춤(권장)"을 기본값으로 내놓는데 멀티모듈만 조용히 라이브러리
        #    spacer(S0002…)를 붙여, 양 끝이 벡터와 무관한 서열로 나왔다.
        #    **연결 spacer(module↔module)는 라이브러리가 맞다** — 벡터 접점이 아니고,
        #    15-mer 직교성으로 인접 module 간 mis-assembly를 막는 것이 그 목적이다.
        if getattr(prm, "spacer_mode", "library") == "vector_homology":
            from .vector_homology_spacer import VectorHomologySpacer
            l_outer, r_outer = VectorHomologySpacer(
                params=prm, arm=WorkflowBuilder._insert_arm(insert_method, prm),
            ).make(vector)
            # 양 끝이 풀에서 나오지 않으므로 연결 spacer는 풀 전체에서 고른다.
            connectors = list(selector.pool_ids())[: n - 1]
        else:
            l_outer, r_outer = selector.select(core_seq, vector, insert_method, params=prm)
            connectors = [s for s in selector.pool_ids()
                          if s not in (l_outer, r_outer)][: n - 1]
        if len(connectors) < n - 1:
            raise RuntimeError(f"Not enough connector spacers ({len(connectors)}/{n - 1}).")
        # attach에는 record/ID를 그대로 넘기고, 리포트·테스트가 읽는 boundaries는 **이름**으로.
        bound_items = [l_outer] + connectors + [r_outer]   # 경계 n+1개 (record 또는 ID)
        boundaries = [b.name if hasattr(b, "name") else b for b in bound_items]

        # 2) 각 module: 양 끝 경계 spacer 부착 + prepare_fragments (공통 헬퍼)
        circular = prm.module_assembly_type == "circular"
        if circular and intermediate_vector is None:
            raise ValueError(
                "A circular (intermediate cloning) multi-module run needs an intermediate_vector. "
                "The spacers attach as parts, so one vector covers every module position."
            )
        asm_type = "circular" if circular else "linear"

        records_list, all_designs, spacer_pairs = [], [], []
        _xd, _xt, _sc = [], [], []   # 쌍 단위 교차 검사 결과·분모(module별로 모은다)
        for i, g in enumerate(module_groups):
            parts_sp = selector.attach(g, bound_items[i], bound_items[i + 1])
            recs, design = WorkflowBuilder._prep_module_records(
                parts_sp, params=prm,
                vector=intermediate_vector if circular else None,
                assembly_type=asm_type,
                module_method=module_method,
            )
            records_list.append(recs)
            spacer_pairs.append([parts_sp[0], parts_sp[-1]])   # excision PCR 앵커
            all_designs.extend(design)
            # 쌍 결과는 리스트 속성에 실려 오므로 extend로는 안 따라온다 — 따로 합친다.
            #   ⚠️ module마다 **다른 튜브**다(교차 dimer는 one-pot 범위 안에서만 의미가
            #     있다). module 간 쌍을 새로 만들지 않고, **어느 반응인지 scope로 태깅**해
            #     모은다 — `part_name`이 module 간에 겹쳐서(실측: 2-module에서
            #     `BBa_J23100_1_designed`가 두 번) 태그 없이는 합친 파일에서 구분이 안 된다.
            _scope = f"module {i + 1}/{len(module_groups)}"
            for _h in (getattr(design, "cross_dimers", []) or []):
                _xd.append({**_h, "scope": _scope})
            for _t in (getattr(design, "crosstalk_pairs", []) or []):
                _xt.append(tuple(_t) + (_scope,))
            for _s in (getattr(design, "scans", []) or []):
                _sc.append((_scope,) + tuple(_s[1:]))
            WorkflowBuilder._save_parts(project_dir, parts_sp, recs)   # parts/raw + parts/gg_prepared

        # 3) run_plan: module_assembly×N (→ module_excision×N) → multimodule_assembly → vector_insertion
        all_records = [r for recs in records_list for r in recs]
        extra = [intermediate_vector] if circular else []
        extra += [r for pair in spacer_pairs for r in pair] if circular else []
        wf = BiopartWorkflow(project_dir=project_dir, parts=all_records + [vector] + extra)

        pk = prm.canonical()
        module_ids = [
            derived_id(
                "module",
                [r.name for r in recs]
                + ([intermediate_vector.name] if circular else [])
                + [module_method, asm_type, pk],
            )
            for recs in records_list
        ]

        plan = []
        for recs, mid in zip(records_list, module_ids):
            plan.append({
                "op": "module_assembly",
                "inputs": ([r.name for r in recs]
                           + ([intermediate_vector.name] if circular else [])),
                "params": {"method": module_method, "assembly_type": asm_type},
            })

        if circular:
            # 각 module을 중간 벡터에서 꺼낸다. spacer 앵커로 PCR하면 말단이 spacer가 되어
            # 인접 module과 구분되는 Gibson homology가 된다.
            excised_ids = [
                derived_id("module", [mid, "excision", prm.excision_method, pk])
                for mid in module_ids
            ]
            for mid, pair in zip(module_ids, spacer_pairs):
                plan.append({
                    "op": "module_excision",
                    "inputs": [mid, intermediate_vector.name] + [r.name for r in pair],
                    "params": {"method": prm.excision_method, "assembly_type": "linear"},
                })
            join_inputs = excised_ids
        else:
            join_inputs = module_ids

        # ⚠️ multimodule_assembly의 id는 **입력 record의 name**으로 계산된다.
        #    _save_artifact가 artifact id를 record.name에 심으므로 그 id를 그대로 쓴다.
        #    (구버전은 name이 "assembled_module" 고정이라 그 값을 썼는데, CP7에서 name 규칙이
        #     바뀌면서 예측이 어긋나 "입력 artifact를 찾을 수 없습니다"로 깨졌다)
        mm_id = derived_id("multimodule", list(join_inputs) + ["gibson", "linear", pk])

        plan.append({"op": "multimodule_assembly", "inputs": list(join_inputs),
                     "params": {"method": "gibson", "assembly_type": "linear"}})
        plan.append({"op": "vector_insertion", "inputs": [mm_id, vector.name],
                     "params": {"method": insert_method, "assembly_type": "circular"}})

        results = wf.run_plan(plan, design_params=prm)
        if insert_method == "gibson" and len(results) >= 2:
            # multimodules/gibson_prepared(vector 삽입 직전 amplicon) + reports/vector_primers.csv
            # results[-2] = multimodule_assembly 출력(vector_insertion 직전)
            WorkflowBuilder._save_gibson_prepared(project_dir, results[-2]["record"], vector, "multimodules",
                                                  params=prm)
        wf.write_primer_csv(all_designs)   # 모든 module 조립 프라이머 → reports/primers.csv (주문용)
        _merged = _CrossAware(all_designs)
        _merged.cross_dimers = _xd
        _merged.crosstalk_pairs = _xt
        _merged.scans = _sc
        WorkflowBuilder._write_cross_checks(wf, _merged, params=prm)
        wf.save_report()

        return {
            "workflow": wf,
            "results": results,
            "boundaries": boundaries,
            "intermediate_vector": intermediate_vector.name if intermediate_vector else None,
            "construct": results[-1]["record"] if results else None,
        }

    @staticmethod
    def get_available_workflows() -> list:
        """
        사용 가능한 워크플로우 목록 반환

        Returns:
            워크플로우 이름 리스트

        Example:
            >>> workflows = WorkflowBuilder.get_available_workflows()
            >>> print(workflows)
            ['classical', 'biopart']
        """
        return ['classical', 'biopart']


__all__ = ['WorkflowBuilder']
