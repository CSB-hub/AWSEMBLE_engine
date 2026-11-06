"""
Workflow 공통 기능

Biopart와 Classical 워크플로우가 공유하는 기본 기능을 제공합니다.

핵심 기능:
    - 디렉토리 생성 및 관리
    - 결과 파일 저장
    - 워크플로우 리포트 생성

작성일: 2025-10-30
"""

import csv
import json
from pathlib import Path
from typing import Dict, List, Any, Optional
from pydna.dseqrecord import Dseqrecord

from workflow.identity import Manifest
from assembly.factory import AssemblyFactory


class WorkflowBase:
    """
    Workflow의 공통 기능
    # TODO: 저장 디렉토리명을 고정해야 할 것 같음.
    # parts/gg_prepared, modules/gibson_prepared 이런식으로.
    # 핵심 목적은 프로젝트명은 <user_input_proj_name>_<날짜_시간>은 고정으로 박고 하위에 생성되는 디렉토리들의 종류들을 고정으로 해서 workflow 설정에 따라서 항상
    # 동일한 형태의 디렉토리와 파일명들이 생기게 하는거야. 

    Biopart와 Classical 워크플로우가 공유하는 기능을 제공합니다.

    Args:
        project_dir: 프로젝트 디렉토리 경로

    Attributes:
        project_dir: 프로젝트 디렉토리 Path 객체
        results: 워크플로우 실행 결과 저장

    Example:
        >>> workflow = WorkflowBase(project_dir="my_project")
        >>> workflow.create_directories(["parts", "modules"])

    작성일: 2025-10-30
    """

    def __init__(self, project_dir: str = None):
        """
        WorkflowBase 초기화

        Args:
            project_dir: 프로젝트 디렉토리 경로 (None이면 현재 디렉토리)
        """
        self.project_dir = Path(project_dir) if project_dir else Path.cwd()
        self.results: Dict[str, Any] = {}

    def create_directories(self, subdirs: List[str]) -> None:
        """
        필요한 디렉토리 생성

        프로젝트 디렉토리 아래에 하위 디렉토리들을 생성합니다.
        이미 존재하는 디렉토리는 건너뜁니다.

        Args:
            subdirs: 생성할 하위 디렉토리 리스트

        Example:
            >>> workflow.create_directories([
            ...     "parts/insert",
            ...     "parts/withvector",
            ...     "modules/insert"
            ... ])

        Notes:
            exist_ok=True로 이미 존재하는 디렉토리도 에러 없이 처리
        """
        for subdir in subdirs:
            dir_path = self.project_dir / subdir
            dir_path.mkdir(parents=True, exist_ok=True)

    # J-2: 고정 표준 레이아웃 (모든 워크플로우 공통 버킷). config에 따라 분기만 생성하되
    #      어휘는 고정 → 같은 설정=같은 트리. 계보는 reports/manifest.json이 보유.
    STANDARD_LAYOUT = [
        "parts/raw", "parts/gg_prepared", "parts/gibson_prepared",
        "modules/gg_prepared", "modules/gibson_prepared", "modules/assembled",
        "multimodules/gg_prepared", "multimodules/gibson_prepared", "multimodules/assembled",
        "constructs", "reports", "logs",
    ]

    def create_standard_layout(self) -> None:
        """J-2 고정 버킷 레이아웃 생성 (parts/modules/multimodules/constructs/reports/logs)."""
        self.create_directories(self.STANDARD_LAYOUT)

    # J-3/K-3: build plan 실행 엔진 (진입 일원화).
    #   build plan = stage 연산의 리스트/DAG. 각 항목:
    #     {"op": <stage 메서드명>, "inputs": [artifact_id...], "params": {"method", "assembly_type", ...}}
    #   평평한 버킷엔 계보가 없으므로, artifact registry(id→record)로 DAG 입력을 resolve하고
    #   계보는 reports/manifest.json(Manifest)에 기록한다.
    def run_plan(self, plan: List[Dict[str, Any]], design_params=None) -> List[Dict[str, Any]]:
        """build plan을 실행 — 각 stage op를 메서드로 dispatch.

        op 이름이 곧 stage 메서드명이다(module_assembly / multimodule_assembly ...).
        stage는 `(records, method=, assembly_type=, manifest=)` 시그니처를 따르고
        `{"id","file","record"}` 또는 None을 반환한다. 반환된 record는 registry에
        등록되어 후속 stage의 입력(inputs id)으로 참조된다.

        Args:
            plan: build plan (stage 연산 리스트/DAG)
            design_params: 설계 파라미터 한 벌(`params.DesignParams`). 각 stage로 전달된다.
                미지정 시 기본값(= 기존 하드코딩 동작).

        Returns:
            성공한 stage별 결과 dict 리스트

        Raises:
            ValueError: 알 수 없는 stage op
            KeyError: build plan 입력 artifact를 registry에서 찾지 못함
        """
        # 소프트 경고를 리포트에 남기기 위해 보관한다(화면이 아니라 산출물에 남긴다 —
        # 기본값도 경고를 내므로 매 실행 stdout에 찍으면 아무도 안 읽게 된다).
        self._design_params = design_params
        manifest = Manifest(self.project_dir)

        # artifact registry(id→record): 초기 part records 등록 (있으면)
        registry: Dict[str, Dseqrecord] = {}
        for p in getattr(self, "parts", None) or []:
            registry[p.name] = p

        results: List[Dict[str, Any]] = []
        for step in plan:
            op = step["op"]
            input_ids = step.get("inputs", [])
            params = dict(step.get("params", {}))
            method = params.pop("method", "golden_gate")
            assembly_type = params.get("assembly_type", "linear")

            stage = getattr(self, op, None)
            if stage is None or not callable(stage):
                raise ValueError(f"Unknown stage op: {op}")

            records = [registry.get(i) for i in input_ids]
            if None in records:
                missing = input_ids[records.index(None)]
                raise KeyError(f"Build plan input artifact not found: {missing}")

            out = stage(records, method=method, assembly_type=assembly_type,
                        manifest=manifest, design_params=design_params)
            if out:
                registry[out["id"]] = out["record"]
                results.append(out)

        manifest.write()
        return results

    def save_results(
        self,
        sequences: List[Dseqrecord],
        output_dir: str,
        prefix: str = "seq",
        format: str = "genbank"
    ) -> List[Path]:
        """
        결과 파일 저장

        DNA 서열들을 GenBank 또는 FASTA 형식으로 저장합니다.

        Args:
            sequences: 저장할 DNA 서열 리스트 (Dseqrecord)
            output_dir: 출력 디렉토리 (project_dir 기준 상대 경로)
            prefix: 파일명 prefix (기본값: "seq")
            format: 파일 형식 ("genbank" 또는 "fasta")

        Returns:
            저장된 파일 경로 리스트

        Example:
            >>> sequences = [dseq1, dseq2, dseq3]
            >>> paths = workflow.save_results(
            ...     sequences,
            ...     output_dir="parts/insert",
            ...     prefix="part",
            ...     format="genbank"
            ... )
            >>> print(paths)
            [Path('parts/insert/part_0001.gb'), ...]

        Notes:
            utils/io.py의 write_genbank, write_fasta 사용
        """
        from awsemble.utils.io import DNAFileWriter

        output_path = self.project_dir / output_dir
        output_path.mkdir(parents=True, exist_ok=True)

        saved_paths = []

        for i, seq in enumerate(sequences):
            filename = f"{prefix}_{i+1:04d}"

            if format == "genbank":
                file_path = output_path / f"{filename}.gb"
                DNAFileWriter.write_genbank(seq, file_path)
            elif format == "fasta":
                file_path = output_path / f"{filename}.fasta"
                DNAFileWriter.write_fasta([seq], file_path)
            else:
                raise ValueError(f"Unsupported format: {format}. Only 'genbank' or 'fasta' is allowed.")

            saved_paths.append(file_path)

        return saved_paths

    # K-4: 실행 결과 리포트 — reports/manifest.json(계보의 진실의 원천)을 소비.
    def generate_report(self) -> str:
        """K-4: manifest 소비 실행 결과 리포트 생성.

        reports/manifest.json을 읽어 stage(op)별 집계 + artifact 계보 표를
        Markdown으로 만든다. (실험 프로토콜 생성기는 후순위 — 별개 산출물.)

        Returns:
            Markdown 형식의 리포트 문자열
        """
        manifest_path = self.project_dir / "reports" / Manifest.FILENAME
        lines = ["# Workflow run report\n", f"\n- Project: `{self.project_dir}`\n"]

        if not manifest_path.exists():
            lines.append("\n(No manifest.json - no stage was run or recorded.)\n")
            return "".join(lines)

        with open(manifest_path, encoding="utf-8") as f:
            artifacts = json.load(f).get("artifacts", [])

        # stage(op)별 집계
        by_op: Dict[str, int] = {}
        for a in artifacts:
            by_op[a["op"]] = by_op.get(a["op"], 0) + 1

        lines.append("\n## Summary\n")
        lines.append(f"- Artifacts: {len(artifacts)}\n")
        for op, n in by_op.items():
            lines.append(f"- {op}: {n}\n")

        lines.append("\n## Artifact lineage\n")
        lines.append("| id | level | method | inputs | file |\n")
        lines.append("|---|---|---|---|---|\n")
        for a in artifacts:
            inputs = ", ".join(a.get("inputs", []))
            lines.append(
                f"| {a['id']} | {a['level']} | {a.get('method', '')} | {inputs} | {a['file']} |\n"
            )

        return "".join(lines)

    def write_primer_csv(self, design_results: List[Dict[str, Any]], output_file: str = "primers.csv") -> Path:
        """prep(prepare_fragments)의 primer 파라미터를 reports/primers.csv로 저장 (P4 필수).

        Golden Gate prepare_fragments가 반환하는 design 리스트의 프라이머 정보
        (fwd/rev 서열, Tm, amplicon 길이, 효소, 검증 통과 여부)를 CSV로 기록한다.

        Args:
            design_results: `GoldenGateModule.prepare_fragments` 또는
                `GibsonCloning.prepare_parts` 반환 리스트. 각 항목 키: part_name,
                forward_primer, reverse_primer, fwd_tm, rev_tm,
                full_sequence(=amplicon), restriction_enzyme, passed.
                homology 경로는 `production`/`oligo_top`/`oligo_bottom`도 싣는다 —
                올리고로 만드는 fragment는 프라이머가 없으므로 이게 주문 정보다.
            output_file: 출력 파일명 (reports/ 아래)

        Returns:
            저장된 CSV 경로
        """
        out = self.project_dir / "reports" / output_file
        out.parent.mkdir(parents=True, exist_ok=True)

        # 실험자 산출물(item 6): 프라이머 서열·Tm·GC% + 2차구조(self-hairpin ΔG) +
        #   template 결합부위 접근성(accessibility) + 전체 tier + 경고까지 담는다.
        fields = [
            "part_name", "restriction_enzyme",
            "left_oh_id", "left_oh", "right_oh_id", "right_oh",
            "forward_primer", "reverse_primer",
            # 2026-11-05(Y10-b) — **주문 단위**가 몇-mer다. 서열 칸을 눈으로 세게 하면
            #   오발주로 이어진다(특히 tail이 붙어 50nt를 넘는 조립용 프라이머).
            #   ⚠️ 주문할 **실물 전체 길이**다 — tail 포함. anneal 길이는 `fwd_tm`/`rev_tm`이
            #      가리키는 구간이고, 둘이 다르다는 건 `design_conditions.md`가 설명한다.
            "fwd_mer", "rev_mer",
            # ⚠️ Tm은 두 가지를 **둘 다** 낸다. 조립용 프라이머는 [5' tail][anneal] 구조라
            #    "Tm 65℃"가 어느 쪽을 말하는지에 따라 완전히 다른 프라이머를 가리킨다.
            #      fwd_tm / rev_tm           : anneal(주형 결합부) 기준 — 값이 낮게 나오는 게 정상
            #      fwd_tm_full / rev_tm_full : 주문할 프라이머 전체 기준
            "fwd_tm", "rev_tm", "fwd_tm_full", "rev_tm_full",
            # 2026-11-05(Y10) — **권장 annealing 온도**. NEB 축으로 간 이유가 이것이다:
            #   Tm을 NEB와 같게 재현해도 그 값을 PCR 조건으로 옮기는 규칙(Ta)이 없으면
            #   유저는 여전히 "그래서 몇 도로 돌리냐"를 손으로 계산해야 한다.
            #   효소별 규칙은 `params.NEB_TA_RULES`(Q5 min+1 / Phusion 0.93·min+7.5 / Taq min−5).
            #   ⚠️ 기준은 **anneal Tm**(꼬리 제외)이다 — 주형에 실제로 붙는 부분이 결합을 정한다.
            "ta",
            # 2026-11-05(Y11) — **fwd·rev Tm 차이**. 위 `ta`와 **붙여서 읽어야** 뜻이
            #   통한다: Ta는 낮은 쪽 Tm 기준으로 정해지므로, 이 값이 크면 그 Ta에서
            #   뜨거운 쪽 프라이머가 자기 Tm보다 훨씬 낮은 온도에서 돌아 비특이 결합이
            #   늘어난다. 선언창은 `pair_tm_diff_max`(기본 5℃)이고 초과하면 경고가 붙는다.
            #   ⚠️ `tm_diff_max`(같은 프라이머의 anneal↔전체)와 **다른 축**이다.
            "tm_diff",
            "fwd_gc", "rev_gc",
            "fwd_hairpin_dg", "rev_hairpin_dg",
            "fwd_accessibility", "rev_accessibility",
            # 올리고 어닐링으로 만드는 fragment는 프라이머가 없다. 그 fragment의
            # **주문할 실물**은 Top/Bottom 두 가닥이므로 여기 같이 적는다.
            #   ⚠️ 빼두면 설계는 끝났는데 주문서에 안 실려 "만들 방법이 없는 fragment"가 된다
            #      (2026-10-25: vector homology spacer 2개가 이렇게 통째로 누락돼 있었다).
            "production", "oligo_top", "oligo_bottom", "oligo_len", "duplex_grade",
            "amplicon_len", "tier", "passed",
            # 2026-11-05(Y12) — **차단 사유**. 2026-11-05까지 이 컬럼이 아예 없어서
            #   tier=red / passed=False가 떠도 **시트에는 이유가 없었다**(엔진은 계속
            #   `violations`를 채우고 있었는데 리포트가 안 꺼내갔다). 차단은 경고보다
            #   중한 신호인데 경고만 보이고 차단은 안 보이던 상태였다.
            "violations", "warnings",
        ]

        def _pct(x):
            return round(x * 100, 1) if isinstance(x, (int, float)) else ""

        def _num(x, nd=2):
            return round(x, nd) if isinstance(x, (int, float)) else ""

        def _mer(seq):
            """주문할 프라이머 길이(nt) — 올리고 행처럼 프라이머가 없으면 비운다."""
            return len(seq) if isinstance(seq, str) and seq else ""

        def _ta(tm_f, tm_r):
            """권장 Ta — 두 anneal Tm으로 효소 규칙을 적용한다.

            파라미터가 NEB 프리셋 조성과 정확히 일치하지 않으면(랩에서 염을 바꿨을 때)
            `ta_recommended`가 None을 돌려준다. 그때는 **추측해서 숫자를 적지 않고 비운다**
            — 틀린 Ta는 없는 Ta보다 나쁘다.
            """
            _p = getattr(self, "_design_params", None)
            if _p is None or not isinstance(tm_f, (int, float)) or not isinstance(tm_r, (int, float)):
                return ""
            v = _p.ta_recommended(tm_f, tm_r)
            return round(v, 1) if v is not None else ""

        def _tm_diff(tm_f, tm_r):
            """fwd·rev Tm 차이(℃) — 올리고 행처럼 한쪽이 없으면 비운다."""
            if not isinstance(tm_f, (int, float)) or not isinstance(tm_r, (int, float)):
                return ""
            return round(abs(tm_f - tm_r), 1)

        def _access(a):
            # accessibility dict → "severity(미페어확률)" 예: "ok(0.78)" / "strong(0.18)" / "na"
            if not isinstance(a, dict):
                return ""
            sev = a.get("severity", "")
            mu = a.get("mean_unpaired")
            return f"{sev}({mu:.2f})" if isinstance(mu, (int, float)) and sev != "na" else sev

        with open(out, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            for d in design_results:
                writer.writerow({
                    "part_name": d.get("part_name", ""),
                    "restriction_enzyme": d.get("restriction_enzyme", ""),
                    "left_oh_id": d.get("left_oh_id", ""),
                    "left_oh": d.get("left_oh", ""),
                    "right_oh_id": d.get("right_oh_id", ""),
                    "right_oh": d.get("right_oh", ""),
                    "forward_primer": d.get("forward_primer", ""),
                    "reverse_primer": d.get("reverse_primer", ""),
                    "fwd_mer": _mer(d.get("forward_primer")),
                    "rev_mer": _mer(d.get("reverse_primer")),
                    "fwd_tm": _num(d.get("fwd_tm")),
                    "rev_tm": _num(d.get("rev_tm")),
                    "fwd_tm_full": _num(d.get("fwd_tm_full")),
                    "rev_tm_full": _num(d.get("rev_tm_full")),
                    "ta": _ta(d.get("fwd_tm"), d.get("rev_tm")),
                    "tm_diff": _tm_diff(d.get("fwd_tm"), d.get("rev_tm")),
                    "fwd_gc": _pct(d.get("fwd_gc")),
                    "rev_gc": _pct(d.get("rev_gc")),
                    "fwd_hairpin_dg": _num(d.get("fwd_hairpin_dg")),
                    "rev_hairpin_dg": _num(d.get("rev_hairpin_dg")),
                    "fwd_accessibility": _access(d.get("fwd_accessibility")),
                    "rev_accessibility": _access(d.get("rev_accessibility")),
                    # production이 없는 설계(Golden Gate)는 전부 PCR이다.
                    "production": d.get("production", "pcr"),
                    "oligo_top": d.get("oligo_top", ""),
                    "oligo_bottom": d.get("oligo_bottom", ""),
                    "oligo_len": d.get("oligo_len", ""),
                    "duplex_grade": d.get("duplex_grade", ""),
                    "amplicon_len": len(d.get("full_sequence", "")),
                    "tier": d.get("tier", ""),
                    "passed": d.get("passed", ""),
                    "violations": " | ".join(d.get("violations", []) or []),
                    "warnings": " | ".join(d.get("warnings", []) or []),
                })

        self._write_design_conditions(design_results)
        return out

    def write_cross_checks_csv(self, out, overhang_checks=None, arm_hits=None,
                               crosstalk_pairs=None, scans=None, params=None):
        """접합부 단위 검사 결과를 `reports/cross_checks.csv`로 쓴다.

        **왜 별 파일인가**: `primers.csv`는 행이 부품 하나다. 그런데 이 검사들은 전부
        **접합부/쌍 단위**라 담을 자리가 없다 — "부품 A와 C 사이"는 A의 속성도 C의
        속성도 아니다. 억지로 한쪽 행에 넣으면 다른 쪽을 보는 사람이 놓친다.
        그래서 접합부를 1행으로 두는 파일을 따로 두고, `primers.csv`의 양쪽 행에는
        **이 파일을 가리키는 포인터**를 남긴다.

        Args:
            out: 출력 경로
            overhang_checks: Golden Gate 접합부 overhang 판정
                (`GoldenGateModule.last_overhang_checks` — 제한효소 사이트/codon 교체·차단)
            arm_hits: homology arm이 엉뚱한 자리에 붙을 뻔한 건과 처리 결과
                (`JunctionPlan.arm_hits` — 재설계 또는 설계불가)
            crosstalk_pairs: `JunctionPlan.crosstalk_pairs`
                (비인접 fragment **몸통** 공유 — 기록 전용)
            scans: `[(scope, kind, 검사한 수), ...]` — **무엇을 몇 건 봤는지**.
                검사마다 요약 행을 한 줄 먼저 쓴다(아래 설명).
            params: DesignParams — 요약 행에 적을 기준값. None이면 기본값.

        Returns:
            쓴 행 수. 0이면 파일을 만들지 않는다.

        Note:
            **요약 행을 왜 넣는가** (2026-11-06). 처음엔 걸린 쌍이 0이면 파일을 아예
            안 만들었다. 그런데 Galaxy 출력은 파일이 있어야 해서 헤더만 있는 파일이
            나갔고, 유저가 받아 보면 "검사했는데 깨끗한 것"인지 "검사가 안 돈 것"인지
            구분할 수 없었다. 그래서 `scans`가 있으면 **걸린 게 0이어도** 파일을
            쓰고, 검사마다 `..._summary` 행으로 "N건 중 0건"을 적는다. 분모가 있어야
            0이 결과로 읽힌다.

            **`primer_dimer` 행은 2026-11-06에 없어졌다.** 부품을 하나씩 따로 PCR하므로
            다른 부품의 프라이머와 같은 튜브에 있을 일이 없다 — 일어나지 않는 일을
            검사해 차단하고 있었다. 한 쌍의 fwd↔rev dimer는 `primers.csv`의 violations에
            그대로 남는다.
        """
        import csv as _csv
        rows = []
        for h in (overhang_checks or []):
            rows.append({
                "scope": h.get("scope", ""),
                "kind": "overhang_junction",
                "a": h.get("junction", ""), "b": "",
                "metric": h.get("metric", ""),
                "value": h.get("value", ""),
                "verdict": h.get("verdict", ""),
                "detail": h.get("detail", ""),
            })
        for h in (arm_hits or []):
            rows.append({
                "scope": h.get("scope", ""),
                "kind": "arm_crosstalk",
                "a": h.get("a", ""), "b": h.get("b", ""),
                "metric": f"shared run (bp, {h.get('orient', '')})",
                "value": h.get("bp", ""),
                "verdict": h.get("verdict", ""),
                "detail": h.get("detail", ""),
            })
        for item in (crosstalk_pairs or []):
            a, b, shared, orient = item[:4]
            rows.append({
                "scope": item[4] if len(item) > 4 else "",
                "kind": "body_crosstalk",
                "a": a, "b": b,
                "metric": f"shared run (bp, {orient})",
                "value": shared,
                # 차단하지 않는다 — arm 길이로 고칠 수 없는 축이다(부품을 바꿔야 한다).
                "verdict": "noted",
                "detail": "Non-adjacent fragment bodies share sequence. The homology arms "
                          "were kept clear of each other, so the assembly is not blocked - "
                          "but in vivo these positions can still recombine. Check whether "
                          "the same part is used in two places.",
            })
        # ── 요약 행 — 검사마다 한 줄, 걸린 게 0이어도 쓴다 ────────────────────
        from params import DEFAULT_PARAMS
        from workflow.junction_plan import JUNCTION_UNIQUE_MAX_BP
        _p = params or DEFAULT_PARAMS
        _why = {
            "overhang_junction": (
                "Every part-to-part junction was checked for an overhang that completes an "
                "enzyme recognition sequence with the flanking part ends, and - where the "
                "junction is translated - for a stop or start codon created by the overhang "
                "bases. A junction that failed was given a different overhang from the same "
                "pool in a single pass."),
            "arm_crosstalk": (
                f"Each homology arm was compared against the other arms (and, for in vivo "
                f"assembly, against the other fragment bodies and the backbone) in both "
                f"forward and reverse-complement orientation, threshold "
                f"{JUNCTION_UNIQUE_MAX_BP}bp. A junction that matched elsewhere was "
                f"redesigned at a different length; if no length avoided it, the design "
                f"stopped."),
            "body_crosstalk": (
                f"The longest run shared between **non-adjacent fragment bodies** was "
                f"compared in both orientations (threshold {JUNCTION_UNIQUE_MAX_BP}bp). "
                f"Adjacent pairs are excluded because sharing is the design intent there. "
                f"This is recorded, not blocked - arm length cannot fix it."),
        }
        summary = []
        for scope, kind, n_pairs in (scans or []):
            n_hit = sum(1 for r in rows
                        if r["scope"] == scope and r["kind"] == kind)
            summary.append({
                "scope": scope,
                "kind": f"{kind}_summary",
                "a": "(all junctions)" if kind == "overhang_junction" else "(all pairs)",
                "b": "",
                "metric": "junctions checked" if kind == "overhang_junction" else "pairs checked",
                "value": n_pairs,
                # 분모가 0이면 "깨끗하다"가 아니라 **비교할 것이 없었다**는 뜻이다
                # (접합부가 하나뿐이면 arm끼리 비교할 짝이 없다). 둘을 같은 말로
                # 적으면 분모를 넣은 이유가 사라진다.
                "verdict": ("nothing to compare" if not n_pairs
                            else "clean" if n_hit == 0 else f"{n_hit} need review"),
                "detail": _why.get(kind, ""),
            })
        rows = summary + rows

        if not rows:
            return 0
        fields = ["scope", "kind", "a", "b", "metric", "value", "verdict", "detail"]
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = _csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
        return len(rows)

    def write_design_events_csv(self, out=None) -> int:
        """`reports/design_events.csv` — 설계 중의 **우회·컷오프 기록**.

        `cross_checks.csv`가 "접합부를 검사한 결과"라면 이 파일은 "엔진이 무엇을
        우회했는지"다. 효소 전환, overhang 교체, anneal shift, spacer graceful
        degradation처럼 **결과물만 보면 알 수 없는 판단**이 여기 남는다
        (`workflow.design_events` 독스트링 참고 — 그 전에는 stdout으로만 나갔다).

        Args:
            out: 출력 경로. None이면 `reports/design_events.csv`.

        Returns:
            기록된 사건 수(0이어도 파일은 쓴다).
        """
        from . import design_events

        path = out or (self.project_dir / "reports" / design_events.FILENAME)
        return design_events.write_csv(path)

    def _write_design_conditions(self, design_results: List[Dict[str, Any]]) -> Optional[Path]:
        """reports/design_conditions.md — 이 설계에 쓰인 조건을 명시한다.

        **왜 필요한가**: 같은 프라이머라도 Tm 계산 도구·염 농도가 다르면 값이 몇 도씩 다르다.
        조건을 안 적어두면 실험자가 "우리 계산기는 72℃인데 왜 68℃냐"를 확인할 방법이 없다.
        (연구실마다 값이 다르다는 게 이 파라미터화의 출발점이었다)
        """
        if not design_results:
            return None
        d0 = design_results[0]
        na = d0.get("tm_na_mM")
        conc = d0.get("tm_primer_nM")
        enzymes = sorted({d.get("restriction_enzyme", "") for d in design_results} - {""})
        prm = getattr(self, "_design_params", None)

        # 2026-11-05(Y10) — Tm 축이 **NEB Tm Calculator 재현**으로 바뀌었다.
        #   그 전엔 "Biopython Tm_NN + Na 50mM"이라고만 적어서, 유저가 NEB에서 본 값과
        #   다를 때 "조건 차이일 수 있다"는 말밖에 못 했다. 지금은 같은 값이 나오므로
        #   **무엇을 재현했고 오차가 얼마인지**를 숫자로 적는다 — 그게 신뢰의 근거다.
        _poly = prm.polymerase_key() if prm is not None else None
        _mono = (prm.na_mM + prm.k_mM + prm.tris_mM / 2.0) if prm is not None else None
        #: NEB 공개값과의 실측 오차 — `tools_dev/params_cp_tests/y10_test.py`가 재검증한다.
        _NEB_RMSE = {"taq": "0.29", "phusion": "0.62", "q5": "0.35"}

        lines = [
            "# Conditions used for this design",
            "",
            "## Tm calculation - reproducing the NEB Tm Calculator",
            "",
            "| Item | Value |",
            "|---|---|",
            "| Method | Nearest-Neighbor (SantaLucia 1998 / Biopython `Tm_NN`, `DNA_NN3`) |",
        ]
        if prm is not None:
            lines += [
                f"| Reproducing | {_poly or '(custom buffer - not an NEB preset)'} |",
                f"| Monovalent salt (converted) | {_mono:g} mM = Na+ {prm.na_mM:g} + K+ {prm.k_mM:g}"
                f" + Tris {prm.tris_mM:g}/2 |",
                f"| Mg²⁺ / dNTP | {prm.mg_mM:g} mM / {prm.dntp_mM:g} mM |",
                f"| Primer concentration | {prm.primer_conc_nM:g} nM |",
                f"| Salt correction | saltcorr={prm.tm_saltcorr}"
                f" ({'Schildkraut - NEB applies this to Phusion only' if prm.tm_saltcorr == 1 else 'Owczarzy 2008'}) |",
            ]
            if _poly in _NEB_RMSE:
                lines.append(f"| Error vs the published NEB values | RMSE {_NEB_RMSE[_poly]}C"
                             f" (checked on 18-36nt primers) |")
        else:
            lines += [
                f"| Na+ concentration | {na} mM |" if na is not None else "| Na+ concentration | (not recorded) |",
                f"| Primer concentration | {conc} nM |" if conc is not None
                else "| Primer concentration | (not recorded) |",
            ]
        lines += [
            "",
            "> These conditions are tuned to give **the same values as the NEB Tm Calculator**.",
            "> The monovalent salt in NEB buffers is K+/Tris rather than Na+, so Na+ is 0"
            " - there is salt; the converted value in the table above is what applies.",
            "",
        ]

        # 권장 Ta — NEB 축으로 간 **목적**이다. Tm만 맞춰도 "몇 도로 돌리냐"가 없으면
        #   유저가 손으로 환산해야 하고, 그러면 축을 맞춘 의미가 절반 사라진다.
        if prm is not None and _poly is not None:
            # 지연 import — `params`는 `basic/`을 쓰고 `workflow/`는 양쪽을 쓰므로
            #   모듈 최상단에서 끌면 import 순환이 생긴다(기존 관례와 동일).
            from params import NEB_TA_RULES
            _kind, _off, _cap = NEB_TA_RULES[_poly]
            _rule = ("0.93 × min(Tm) + 7.5" if _kind == "phusion"
                     else f"min(Tm) {'+' if _off >= 0 else '−'} {abs(_off):g}")
            lines += [
                "## Recommended annealing temperature (the `ta` column in `primers.csv`)",
                "",
                "| Enzyme | Rule | Cap |",
                "|---|---|---|",
                f"| {_poly} | {_rule} | {_cap:g}℃ |",
                "",
                "> `min(Tm)` is the **lower of the fwd/rev annealing Tm values**"
                " - whichever lets go first decides the PCR.",
                "> These are NEB's published Ta rules verbatim, so the NEB site gives the same answer.",
                "",
            ]

        lines += [
            "## The two kinds of Tm in primers.csv",
            "",
            "An assembly primer is `[5' tail][template-binding region]`, so the two values differ a lot.",
            "",
            "| Column | Meaning |",
            "|---|---|",
            "| `fwd_tm` / `rev_tm` | the **annealing region** only - tail excluded. A low value is normal |",
            "| `fwd_tm_full` / `rev_tm_full` | the **whole primer you order** - tail included |",
            "",
        ]

        # 2026-11-05(Y11) — 판정 기준이 **전부 유저 선언창 하나**로 바뀌었다.
        # 그러면 "경고가 왜 떴나"의 답이 항상 "내가 설정한 범위를 벗어났다"이므로,
        # 설명할 것은 등급 2층의 뜻이 아니라 **어떤 범위가 쓰였는가**다.
        if prm is not None:
            lines += [
                "## Warning criteria - the single window you set",
                "",
                "A primer warning appears **only when a value falls outside the window below**."
                " Previously there was a second internal judgment window outside your settings"
                " (Tm 45-85C, GC ideal band 40-60%, annealing length fixed at 15-30nt), so"
                " warnings could come from criteria you never set.",
                "",
                "| Axis | Your setting | Compared against |",
                "|---|---|---|",
                f"| Tm | {prm.tm_min:g} - {prm.tm_max:g} C | annealing Tm |",
                f"| Annealing length | {prm.primer_len_min} - {prm.primer_len_max} nt | "
                f"annealing length (tail excluded) |",
                f"| GC | {prm.gc_min:g} - {prm.gc_max:g} % | annealing-region GC% |",
                f"| fwd/rev Tm difference | {prm.pair_tm_diff_max:g} C or less | "
                f"the `tm_diff` column |",
                "",
                "> The annealing-region GC% is **a property of the template**, so extending the"
                " primer barely changes it (it does raise the Tm). When a GC warning appears,"
                " judge it together with whether the Tm is sufficient rather than reworking"
                " the primer.",
                "",
                "> Two things are reported **regardless** of the window above: hairpin and dimer."
                " They are not a matter of range - dG alone says the PCR will not run - so"
                " they are reported even when every value sits inside your settings.",
                "",
                "## Read `tm_diff` and `ta` in primers.csv together",
                "",
                f"The recommended annealing temperature (`ta`) is set from the **lower Tm**."
                f" So when `tm_diff` is large, at that Ta the **hotter primer runs well below"
                f" its own Tm** and non-specific binding increases."
                f" The limit you set is {prm.pair_tm_diff_max:g}C.",
                "",
            ]

        lines += [
            "## Assembly",
            "",
            f"| Restriction enzymes used | {', '.join(enzymes) if enzymes else '(none)'} |",
            "|---|---|",
            "",
        ]

        # 설계 파라미터 소프트 경고 — 차단은 아니지만 "왜 이런 값이 나왔나"의 단서다.
        # 하드 위반은 여기 오기 전에 build_plasmid 진입부에서 이미 에러로 막힌다.
        if prm is not None:
            try:
                warns = prm.validate()
            except ValueError as e:      # 방어적 — 진입부에서 이미 걸렀어야 한다
                warns = [f"(validation failed) {e}"]
            if warns:
                lines += [
                    "## Design parameter warnings",
                    "",
                    "> These are **notes, not blocks** - the design used values outside the recommended range.",
                    "",
                ] + [f"- {w}" for w in warns] + [""]
        out = self.project_dir / "reports" / "design_conditions.md"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text("\n".join(lines), encoding="utf-8")
        return out

    def get_result(self, step_name: str) -> Any:
        """
        특정 단계의 결과 반환

        Args:
            step_name: 단계 이름 (예: "step1", "step2")

        Returns:
            해당 단계의 결과 (저장되어 있으면)

        Raises:
            KeyError: 해당 단계 결과가 없는 경우

        Example:
            >>> result = workflow.get_result("step1")
        """
        if step_name not in self.results:
            raise KeyError(f"No result for step '{step_name}'")
        return self.results[step_name]

    def save_report(self, output_file: str = "run_report.md") -> Path:
        """K-4: 실행 결과 리포트를 reports/ 아래 파일로 저장.

        Args:
            output_file: 출력 파일명 (reports/ 아래, 기본값 "run_report.md")

        Returns:
            저장된 파일 경로
        """
        report = self.generate_report()
        report_path = self.project_dir / "reports" / output_file
        report_path.parent.mkdir(parents=True, exist_ok=True)

        with open(report_path, 'w', encoding='utf-8') as f:
            f.write(report)

        return report_path

    # ------------------------------------------------------------------ #
    # P7: vector insertion (insert/module/multimodule + vector → circular construct).
    #   Classical/Biopart 공유 stage → 공통 부모에 둔다(형제 모두 run_plan으로 dispatch).
    #   assemble: method별 조립 핵심 / _ligate_circular: GG sticky-end / vector_insertion: stage 계약.
    # ------------------------------------------------------------------ #
    def assemble(
        self,
        insert: Dseqrecord,
        vector: Dseqrecord,
        method: str = "gibson",
        validate: bool = True,
        design_params=None,
        **kwargs
    ) -> List[Dseqrecord]:
        """insert + vector를 지정 방법으로 circular 조립.

        Args:
            insert: Insert DNA (Dseqrecord) — 일반 insert / biopart module / multimodule
            vector: Vector 서열 (Dseqrecord)
            method: "gibson" 또는 "golden_gate"
            validate: 검증 수행 여부
            **kwargs: gibson=overlap_length / golden_gate=enzyme

        Returns:
            조립된 plasmid 서열 리스트 (Dseqrecord)

        Notes:
            검증 실패 시 경고만 출력하고 조립은 진행.
        """
        print(f"=== Vector Insertion ({method}) ===\n")

        # Assembler 생성 — 설계 파라미터(연구실 조건)를 그대로 넘긴다.
        # 넘기지 않으면 assembler가 기본값으로 돌아 "입력했는데 반영이 안 되는" 상태가 된다.
        if design_params is not None:
            kwargs.setdefault("params", design_params)
            # golden_gate는 **효소까지** 넘겨야 한다. 안 넘기면 assembler가 기본 BsaI로
            # 돌아, 앞에서 `select_gg_enzyme`이 벡터·부품을 보고 정한 결과가 버려진다
            # (BsaI 사이트가 없는 벡터는 그대로 "자를 수 없다"로 죽었다).
            if method.lower() == "golden_gate":
                kwargs.setdefault("enzyme", design_params.gg_enzyme)
        assembler = AssemblyFactory.create_assembler(method, **kwargs)

        # 검증
        if validate:
            print("Validating...")
            validation = assembler.validate(insert)
            if not validation['valid']:
                print("Validation warnings:")
                for warning in validation.get('warnings', []):
                    print(f"  - {warning}")
                for error in validation.get('errors', []):
                    print(f"  ❌ {error}")
                print()

        # 프라이머 준비 및 조립 — **homology 기반 방법 공통**(Gibson · yeast HR).
        #
        # ⚠️ 여기를 `method == "gibson"`으로 좁히면 yeast HR이 이 블록을 통째로 건너뛰어
        #    **벡터가 선형화되지 않은 채** 맨 아래 else로 떨어진다. 원형 벡터는 말단이
        #    없어서 pydna가 homology를 물릴 자리를 못 찾고, 예외 없이 빈 결과가 나온다
        #    (2026-10-25 실측: yeast HR 삽입이 100% "결과물이 생성되지 않았습니다"로 죽음).
        #
        # 그래서 **방법 이름이 아니라 "homology arm을 쓰는가"로 가른다** — `ARM_PARAMS` 훅.
        # 새 homology 기반 방법을 추가해도 이 분기를 다시 고칠 필요가 없다.
        if getattr(assembler, "ARM_PARAMS", None):
            label = assembler.__class__.__name__
            # module 양 끝이 이미 벡터 homology면 **삽입 PCR이 필요 없다.**
            #   플래그가 아니라 서열을 직접 대조한다 — 실제로 맞을 때만 건너뛰므로
            #   모드 설정이 꼬여도 잘못된 산물이 나오지 않는다.
            from workflow.vector_homology_spacer import module_carries_homology
            _lin = assembler._linearize_vector(vector)
            if module_carries_homology(insert, _lin, assembler.overlap_length):
                print("The module already carries vector homology - skipping the insertion PCR")
                print(f"  Linearized vector: {len(_lin.seq)} bp / insert: {len(insert.seq)} bp\n")
                print(f"Assembling {label}...")
                return assembler.assemble([_lin, insert], assembly_type="circular")

            print("Designing primers...")
            fwd, rev = assembler.prepare_fragments(insert, vector)
            print(f"  Forward: {fwd[:30]}..." if len(fwd) > 30 else f"  Forward: {fwd}")
            print(f"  Reverse: {rev[:30]}..." if len(rev) > 30 else f"  Reverse: {rev}")
            print()

            # PCR 시뮬레이션
            print("Simulating PCR...")
            pcr_product = assembler.amplify_insert(insert, fwd, rev)
            print(f"  PCR product: {len(pcr_product.seq)} bp")
            print()

            # Vector 선형화
            print("Linearizing the vector...")
            linear_vector = assembler._linearize_vector(vector)
            print(f"  Linearized vector: {len(linear_vector.seq)} bp")
            print()

            # homology 조립 (Gibson / yeast HR)
            print(f"Assembling {label}...")
            results = assembler.assemble([linear_vector, pcr_product], assembly_type="circular")
        elif method.lower() == "golden_gate":
            # Golden Gate vector insertion (P5): vector를 조립 효소로 잘라 backbone을 얻고,
            # insert를 backbone overhang에 맞춰 같은 효소 프라이머로 증폭·절단한 뒤
            # sticky-end ligation으로 circular 조립. (구 select= 호출 = 구멍2 교체)
            #
            # ⚠️ biopart module 조립(prepare_fragments/GoldenGateModule.assemble)은
            #    불변 — 여기서는 vector_insertion 전용 레시피를 직접 구성한다.
            #
            # 효소는 **여기서 고르지 않는다.** one-pot이라 벡터·부품·module이 전부 같은
            # 효소여야 하므로 `select_gg_enzyme`이 앞에서 한 번 정하고, 이 블록은
            # assembler가 들고 있는 그 결정을 읽기만 한다.
            from Bio.Seq import Seq
            from basic.primer_design import PrimerDesign

            gg_enz = assembler.enzyme            # Biopython 효소 객체
            gg_name = assembler.enzyme_name      # 프라이머 설계용 이름
            print(f"Golden Gate vector insertion... (enzyme: {gg_name})")

            # 1. insert에서 target CDS 좌표 (Gibson과 동일한 feature 자동감지 재사용)
            gibson_helper = AssemblyFactory.create_assembler("gibson")
            target_start, target_end = gibson_helper._find_target_feature(insert, None, None)
            target = insert[target_start:target_end]
            target.name = insert.name

            # 2. vector를 조립 효소로 잘라 backbone overhang 추출
            #    insert overhang = left:v_right / right:RC(v_left)
            #    (sticky-end 가닥 정합 — P5 GG 결정. RC로 backbone과 맞물림)
            # 벡터를 원형으로 정규화 — DB가 주는 vector는 circular=False다.
            # 선형으로 자르면 backbone 한쪽이 blunt가 되어 sticky-end ligation이 불가능하다
            # (2026-09-23 실측: V0001 선형 3조각/blunt vs 원형 2조각/양끝 sticky).
            vector = assembler._as_circular(vector)
            v_left, v_right = assembler.extract_vector_overhangs(vector)
            rc = lambda s: str(Seq(s).reverse_complement())

            # 3. insert를 backbone overhang에 맞춰 같은 효소 프라이머로 증폭 → 절단
            design = PrimerDesign.design_goldengate_primers(
                target,
                left_overhang=v_right,
                right_overhang=rc(v_left),
                restriction_enzyme=gg_name,
            )
            amplicon = Dseqrecord(design["full_sequence"])
            insert_fragment = max(amplicon.cut(gg_enz), key=lambda f: len(f.seq))
            backbone = max(vector.cut(gg_enz), key=lambda f: len(f.seq))
            print(f"  Insert fragment: {len(insert_fragment.seq)} bp | Backbone: {len(backbone.seq)} bp")

            # 4. sticky-end ligation → circular construct
            results = [self._ligate_circular(backbone, insert_fragment)]
        else:
            # 기타 조립 방법 (향후 추가)
            print("Simulating assembly...")
            results = assembler.assemble([vector, insert], assembly_type="circular")

        print(f"{len(results)} constructs produced\n")
        return results

    def _ligate_circular(self, backbone: Dseqrecord, insert_fragment: Dseqrecord) -> Dseqrecord:
        """sticky-end ligation으로 circular construct 생성 (Golden Gate vector insertion).

        backbone과 insert를 pydna `+`(sticky-end 호환 ligation)로 이어 circular화한다.
        문자열 병합과 달리 pydna는 4bp overhang의 **가닥(complementary)까지** 따져
        호환될 때만 결합하므로, 호환되는 방향/순서 조합을 탐색해 첫 성공을 반환한다.

        Args:
            backbone: 조립 효소로 절단된 vector backbone fragment (sticky-end 보존)
            insert_fragment: 같은 효소로 절단된 insert fragment (sticky-end 보존)

        Returns:
            circular Dseqrecord

        Raises:
            RuntimeError: 어떤 방향 조합으로도 sticky-end가 호환되지 않을 때
        """
        for a, b in [
            (backbone, insert_fragment),
            (insert_fragment, backbone),
            (backbone, insert_fragment.reverse_complement()),
            (insert_fragment.reverse_complement(), backbone),
        ]:
            try:
                return (a + b).looped()
            except Exception:
                continue
        raise RuntimeError(
            "Golden Gate circular ligation failed: the backbone and insert sticky ends are not compatible."
        )

    def vector_insertion(
        self,
        records: List[Dseqrecord],
        method: str = "gibson",
        assembly_type: str = "circular",
        vi_id: Optional[str] = None,
        manifest=None,
        design_params=None,
    ) -> Optional[Dict[str, Any]]:
        """stage: 1 insert + 1 vector → 1 circular construct (run_plan dispatch 호환).

        records=[insert, vector]. assemble로 circular 조립 후 constructs/에 저장하고
        manifest에 계보(level="construct")를 기록한다. module_assembly와 같은
        stage 계약 `(records, method=, assembly_type=, manifest=) → {id,file,record}|None`.
        (P7 Step4 = biopart module/multimodule을 vector에 삽입하는 이 stage의 일반화.)

        Args:
            records: [insert, vector] (run_plan이 inputs id를 resolve해 전달)
            method: 조립 방법 ("gibson" | "golden_gate")
            assembly_type: 조립 형태 (vector insertion은 "circular")
            vi_id: 명시 id. None이면 결정적 해시 con_<hash8>.
            manifest: 계보 기록 대상

        Returns:
            {"id","file","record"} 또는 실패 시 None.
        """
        from utils.io import DNAFileWriter
        from workflow.identity import derived_id

        insert, vector = records[0], records[1]
        results = self.assemble(insert, vector, method=method, validate=False,
                                design_params=design_params)
        if not results:
            return None
        construct = results[0]

        # P6: 소실된 insert 부품 feature를 construct에 기입 후 중복 제거 + source 정정
        #     (조립 불변 — 후처리. annotate로 복원 → remove_duplicate로 vector와 겹친 중복 정리)
        try:
            from assembly.common import annotate_inserted_region, remove_duplicate_features
            from Bio.SeqFeature import FeatureLocation, SeqFeature
            helper = AssemblyFactory.create_assembler("gibson")
            ts, te = helper._find_target_feature(insert, None, None)
            construct = annotate_inserted_region(construct, insert, ts, te)
            construct = remove_duplicate_features(construct)
            # source feature를 construct 전체로 보정(있으면 좌표 정정, 없으면 추가).
            #   gibson: pydna가 vector source(옛 vector 길이)를 물려받아 insert 구간을 안 덮음 → 정정
            #   golden_gate: sticky-end ligation은 source를 안 만듦 → [0,len] 추가
            n = len(construct.seq)
            sources = [f for f in construct.features if f.type == "source"]
            if sources:
                for f in sources:
                    f.location = FeatureLocation(0, n, strand=f.location.strand)
            else:
                construct.features.insert(0, SeqFeature(FeatureLocation(0, n, strand=1), type="source"))
        except Exception as e:
            print(f"  Skipped writing insert features: {e}")

        if vi_id is None:
            # 설계 파라미터가 해시 재료에 없으면 다른 조건으로 만든 construct가 같은 id를 갖는다.
            _pk = design_params.canonical() if design_params is not None else ""
            vi_id = derived_id(
                "construct", [insert.name, vector.name, method, assembly_type, _pk]
            )

        out_dir = self.project_dir / "constructs"
        out_dir.mkdir(parents=True, exist_ok=True)
        output_path = out_dir / f"{vi_id}.gb"
        DNAFileWriter.write_genbank(construct, str(output_path))

        if manifest is not None:
            manifest.add(
                id=vi_id,
                level="construct",
                file=f"constructs/{vi_id}.gb",
                op="vector_insertion",
                method=method,
                inputs=[insert.name, vector.name],
                params={
                    "assembly_type": assembly_type,
                    **({"design_params": design_params.to_dict()} if design_params is not None else {}),
                },
            )
        return {"id": vi_id, "file": output_path, "record": construct}


__all__ = ['WorkflowBase']
