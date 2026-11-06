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
                raise ValueError(f"알 수 없는 stage op: {op}")

            records = [registry.get(i) for i in input_ids]
            if None in records:
                missing = input_ids[records.index(None)]
                raise KeyError(f"build plan 입력 artifact를 찾을 수 없습니다: {missing}")

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
                raise ValueError(f"지원하지 않는 형식: {format}. 'genbank' 또는 'fasta'만 가능합니다.")

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
        lines = ["# Workflow 실행 결과 리포트\n", f"\n- 프로젝트: `{self.project_dir}`\n"]

        if not manifest_path.exists():
            lines.append("\n(manifest.json 없음 — 실행/기록된 stage가 없습니다.)\n")
            return "".join(lines)

        with open(manifest_path, encoding="utf-8") as f:
            artifacts = json.load(f).get("artifacts", [])

        # stage(op)별 집계
        by_op: Dict[str, int] = {}
        for a in artifacts:
            by_op[a["op"]] = by_op.get(a["op"], 0) + 1

        lines.append("\n## 요약\n")
        lines.append(f"- 총 artifact: {len(artifacts)}개\n")
        for op, n in by_op.items():
            lines.append(f"- {op}: {n}개\n")

        lines.append("\n## Artifact 계보\n")
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
            # ⚠️ Tm은 두 가지를 **둘 다** 낸다. 조립용 프라이머는 [5' tail][anneal] 구조라
            #    "Tm 65℃"가 어느 쪽을 말하는지에 따라 완전히 다른 프라이머를 가리킨다.
            #      fwd_tm / rev_tm           : anneal(주형 결합부) 기준 — 값이 낮게 나오는 게 정상
            #      fwd_tm_full / rev_tm_full : 주문할 프라이머 전체 기준
            "fwd_tm", "rev_tm", "fwd_tm_full", "rev_tm_full", "fwd_gc", "rev_gc",
            "fwd_hairpin_dg", "rev_hairpin_dg",
            "fwd_accessibility", "rev_accessibility",
            # 올리고 어닐링으로 만드는 fragment는 프라이머가 없다. 그 fragment의
            # **주문할 실물**은 Top/Bottom 두 가닥이므로 여기 같이 적는다.
            #   ⚠️ 빼두면 설계는 끝났는데 주문서에 안 실려 "만들 방법이 없는 fragment"가 된다
            #      (2026-10-25: vector homology spacer 2개가 이렇게 통째로 누락돼 있었다).
            "production", "oligo_top", "oligo_bottom", "oligo_len", "duplex_grade",
            "amplicon_len", "tier", "passed", "warnings",
        ]

        def _pct(x):
            return round(x * 100, 1) if isinstance(x, (int, float)) else ""

        def _num(x, nd=2):
            return round(x, nd) if isinstance(x, (int, float)) else ""

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
                    "fwd_tm": _num(d.get("fwd_tm")),
                    "rev_tm": _num(d.get("rev_tm")),
                    "fwd_tm_full": _num(d.get("fwd_tm_full")),
                    "rev_tm_full": _num(d.get("rev_tm_full")),
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
                    "warnings": " | ".join(d.get("warnings", []) or []),
                })

        self._write_design_conditions(design_results)
        return out

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

        lines = [
            "# 이 설계에 쓰인 조건",
            "",
            "## Tm 계산",
            "",
            "| 항목 | 값 |",
            "|---|---|",
            "| 방식 | Nearest-Neighbor (SantaLucia 1998, Biopython `Tm_NN`) |",
            f"| Na⁺ 농도 | {na} mM |" if na is not None else "| Na⁺ 농도 | (미기록) |",
            f"| 프라이머 농도 | {conc} nM |" if conc is not None else "| 프라이머 농도 | (미기록) |",
            "",
            "> 계산 도구가 다르면 같은 서열도 Tm이 몇 도씩 다르게 나옵니다.",
            "> NEB/SnapGene 등에서 보신 값과 차이가 있다면 위 조건 차이일 가능성이 큽니다.",
            "",
            "## primers.csv의 Tm 두 종류",
            "",
            "조립용 프라이머는 `[5' 꼬리][주형 결합부]` 구조라 두 값이 크게 다릅니다.",
            "",
            "| 컬럼 | 의미 |",
            "|---|---|",
            "| `fwd_tm` / `rev_tm` | **결합부(anneal)** 기준 — 꼬리 제외. 낮게 나오는 것이 정상 |",
            "| `fwd_tm_full` / `rev_tm_full` | **주문할 프라이머 전체** 기준 — 꼬리 포함 |",
            "",
        ]

        # GC 등급은 **2층**이다(차단/경고). 설명이 없으면 yellow가 뜬 이유를 알 수 없다 —
        # "등급을 숨기지 않는다"의 나머지 절반이다. 표의 숫자는 실제 적용값에서 읽는다.
        if prm is not None:
            lines += [
                "## primers.csv의 GC 등급",
                "",
                "`fwd_gc` / `rev_gc`도 **결합부(anneal) 기준**입니다. 등급은 2층입니다.",
                "",
                "| 결합부 GC% | 등급 | 뜻 |",
                "|---|---|---|",
                f"| {prm.qc_gc_min:g} ~ {prm.qc_gc_ideal_min:g} 또는 "
                f"{prm.qc_gc_ideal_max:g} ~ {prm.qc_gc_max:g} | 🟡 주의 | "
                f"관례적 이상 구간 밖입니다. **AT-rich/GC-rich 주형이면 정상**일 수 있습니다 |",
                f"| {prm.qc_gc_ideal_min:g} ~ {prm.qc_gc_ideal_max:g} | 🟢 양호 | 이상 구간 |",
                f"| {prm.qc_gc_min:g} 미만 또는 {prm.qc_gc_max:g} 초과 | 🔴 차단 | "
                f"설계 목표 범위 밖입니다 |",
                "",
                "> 결합부 GC%는 **주형의 성질**이라 프라이머를 늘려도 거의 바뀌지 않습니다"
                " (늘리면 Tm은 올라갑니다). 그래서 이상 구간을 벗어난 것 자체는 차단하지 않고,",
                "> Tm이 위 조건에서 충분한지를 함께 보고 판단하시면 됩니다.",
                "",
            ]

        lines += [
            "## 조립",
            "",
            f"| 사용 제한효소 | {', '.join(enzymes) if enzymes else '(없음)'} |",
            "|---|---|",
            "",
        ]

        # 설계 파라미터 소프트 경고 — 차단은 아니지만 "왜 이런 값이 나왔나"의 단서다.
        # 하드 위반은 여기 오기 전에 build_plasmid 진입부에서 이미 에러로 막힌다.
        if prm is not None:
            try:
                warns = prm.validate()
            except ValueError as e:      # 방어적 — 진입부에서 이미 걸렀어야 한다
                warns = [f"(검증 실패) {e}"]
            if warns:
                lines += [
                    "## 설계 파라미터 경고",
                    "",
                    "> 아래는 **차단이 아니라 참고**입니다. 권장 범위를 벗어난 값으로 설계했다는 뜻입니다.",
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
            raise KeyError(f"단계 '{step_name}'의 결과가 없습니다")
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
            print("검증 중...")
            validation = assembler.validate(insert)
            if not validation['valid']:
                print("⚠️ 검증 경고:")
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
                print("Vector homology가 module에 이미 있습니다 — 삽입 PCR 생략")
                print(f"  선형화된 Vector: {len(_lin.seq)} bp / insert: {len(insert.seq)} bp\n")
                print(f"{label} 조립 수행 중...")
                return assembler.assemble([_lin, insert], assembly_type="circular")

            print("프라이머 설계 중...")
            fwd, rev = assembler.prepare_fragments(insert, vector)
            print(f"  Forward: {fwd[:30]}..." if len(fwd) > 30 else f"  Forward: {fwd}")
            print(f"  Reverse: {rev[:30]}..." if len(rev) > 30 else f"  Reverse: {rev}")
            print()

            # PCR 시뮬레이션
            print("PCR 시뮬레이션 중...")
            pcr_product = assembler.amplify_insert(insert, fwd, rev)
            print(f"  PCR 산물: {len(pcr_product.seq)} bp")
            print()

            # Vector 선형화
            print("Vector 선형화 중...")
            linear_vector = assembler._linearize_vector(vector)
            print(f"  선형화된 Vector: {len(linear_vector.seq)} bp")
            print()

            # homology 조립 (Gibson / yeast HR)
            print(f"{label} 조립 수행 중...")
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
            print(f"Golden Gate vector insertion... (효소: {gg_name})")

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
            print("조립 시뮬레이션 중...")
            results = assembler.assemble([vector, insert], assembly_type="circular")

        print(f"✅ {len(results)}개 construct 생성\n")
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
            "Golden Gate circular ligation 실패: backbone-insert sticky-end가 호환되지 않습니다."
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
            print(f"  ⚠️  insert feature 기입 생략: {e}")

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
