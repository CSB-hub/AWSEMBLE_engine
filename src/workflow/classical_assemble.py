"""
Classical 조립 워크플로우 (일반 클로닝, 구 Simple)

간단한 클로닝 워크플로우를 제공합니다.

핵심 기능:
    - Insert + Vector → Assembled Plasmid
    - Gibson 또는 Golden Gate 선택
    - 빠른 테스트 및 검증
    - GUI 없이 script로 사용

사용 예제:
    >>> workflow = ClassicalWorkflow(output_dir="results")
    >>> workflow.run(
    ...     insert_file="gfp.fasta",
    ...     vector_file="plasmid.gb",
    ...     method="gibson"
    ... )

작성일: 2025-10-30
"""

from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
from pydna.dseqrecord import Dseqrecord

from workflow.base import WorkflowBase


# TODO: 이건 galaxy에 얹을때 입력방식이 완전히 달라질거 같은데... galaxy에 classicial mode 얹을 때 어떻게 얹을지 생각
# 260606_todo.md에서 작성했던 dead code, 로직 아키텍쳐 확인 등을 체크
class ClassicalWorkflow(WorkflowBase):
    """
    간단한 클로닝 워크플로우

    Input: Insert + Plasmid (fasta/gbk)
    → Gibson or Golden Gate
    → Output: Assembled plasmid

    특징:
        - 빠른 테스트용
        - GUI 없이 script로 사용
        - GFP 삽입 같은 간단한 클로닝에 적합

    Args:
        output_dir: 출력 디렉토리 (기본값: "output")

    Example:
        >>> workflow = ClassicalWorkflow(output_dir="results")
        >>> results = workflow.run(
        ...     insert_file="gfp.fasta",
        ...     vector_file="pUC19.gb",
        ...     method="gibson"
        ... )
        >>> print(f"{len(results)}개 construct 생성됨")

    작성일: 2025-10-30
    """

    def __init__(self, output_dir: str = "output"):
        """
        ClassicalWorkflow 초기화

        Args:
            output_dir: 출력(=프로젝트) 디렉토리 경로. P5에서 project_dir과 통일하고
                고정 표준 레이아웃(constructs/reports/...)을 생성한다.
        """
        super().__init__(project_dir=output_dir)
        self.output_dir = Path(output_dir)   # 하위호환 별칭
        self.create_standard_layout()         # J-2 고정 버킷 (constructs/reports 포함)

    def load_sequences(
        self,
        insert_file: str,
        vector_file: str
    ) -> Tuple[Dseqrecord, Dseqrecord]:
        """
        fasta/gbk 파일 읽기

        FASTA 또는 GenBank 형식의 파일을 읽어 Dseqrecord로 반환합니다.

        Args:
            insert_file: Insert DNA 파일 (fasta or gbk)
            vector_file: Vector 파일 (fasta or gbk)

        Returns:
            (insert_dseqrecord, vector_dseqrecord)

        Raises:
            FileNotFoundError: 파일이 존재하지 않는 경우
            ValueError: 지원하지 않는 파일 형식

        Example:
            >>> insert, vector = workflow.load_sequences(
            ...     "gfp.fasta",
            ...     "pUC19.gb"
            ... )
        """
        from utils.io import DNAFileReader

        # Insert 파일 읽기
        insert_path = Path(insert_file)
        if not insert_path.exists():
            raise FileNotFoundError(f"Insert file not found: {insert_file}")

        if insert_path.suffix.lower() in ['.gb', '.gbk', '.genbank']:
            insert_list = DNAFileReader.read_genbank(insert_file)
            insert = insert_list[0] if insert_list else None
        elif insert_path.suffix.lower() in ['.fa', '.fasta', '.fna']:
            insert_list = DNAFileReader.read_fasta(insert_file)
            insert = insert_list[0] if insert_list else None
        else:
            raise ValueError(f"Unsupported insert file format: {insert_path.suffix}")

        if insert is None:
            raise ValueError(f"The insert file has no sequence: {insert_file}")

        # Vector 파일 읽기
        vector_path = Path(vector_file)
        if not vector_path.exists():
            raise FileNotFoundError(f"Vector file not found: {vector_file}")

        if vector_path.suffix.lower() in ['.gb', '.gbk', '.genbank']:
            vector_list = DNAFileReader.read_genbank(vector_file)
            vector = vector_list[0] if vector_list else None
        elif vector_path.suffix.lower() in ['.fa', '.fasta', '.fna']:
            vector_list = DNAFileReader.read_fasta(vector_file)
            vector = vector_list[0] if vector_list else None
        else:
            raise ValueError(f"Unsupported vector file format: {vector_path.suffix}")

        if vector is None:
            raise ValueError(f"The vector file has no sequence: {vector_file}")

        return insert, vector

    def export(
        self,
        sequences: List[Dseqrecord],
        prefix: str = "construct",
        format: str = "both"
    ) -> List[Path]:
        """
        결과 저장

        조립된 서열들을 파일로 저장합니다.

        Args:
            sequences: 저장할 서열 리스트 (Dseqrecord)
            prefix: 파일명 prefix (기본값: "construct")
            format: "gbk", "fasta", or "both" (기본값: "both")

        Returns:
            저장된 파일 경로 리스트

        Example:
            >>> paths = workflow.export(
            ...     results,
            ...     prefix="gfp_plasmid",
            ...     format="both"
            ... )

        Notes:
            파일명 형식: {prefix}_0001.gb, {prefix}_0001.fasta, ...
        """
        from utils.io import DNAFileWriter

        # format 검증
        valid_formats = ["gbk", "fasta", "both"]
        if format not in valid_formats:
            raise ValueError(
                f"Invalid format: {format}. "
                f"Must be one of {valid_formats}"
            )

        print(f"Saving results to: {self.output_dir}/")

        saved_paths = []

        for i, seq in enumerate(sequences):
            filename = f"{prefix}_{i+1:04d}"

            if format in ["gbk", "both"]:
                gbk_path = self.output_dir / f"{filename}.gb"
                DNAFileWriter.write_genbank(seq, gbk_path)
                saved_paths.append(gbk_path)
                print(f"  ✅ {gbk_path.name}")

            if format in ["fasta", "both"]:
                fasta_path = self.output_dir / f"{filename}.fasta"
                DNAFileWriter.write_fasta([seq], fasta_path)
                saved_paths.append(fasta_path)
                print(f"  ✅ {fasta_path.name}")

        print(f"\nSaved {len(saved_paths)} files")
        return saved_paths

    def run(
        self,
        insert_file: str,
        vector_file: str,
        method: str = "gibson",
        output_prefix: str = "construct",
        output_format: str = "both",
        validate: bool = True,
        **kwargs
    ) -> List[Dseqrecord]:
        """
        전체 워크플로우 실행 (한 번에)

        파일 읽기 → 조립 → 저장을 한 번에 수행합니다.

        Args:
            insert_file: Insert 파일 경로 (fasta or gbk)
            vector_file: Vector 파일 경로 (fasta or gbk)
            method: "gibson" 또는 "golden_gate"
            output_prefix: 출력 파일명 prefix
            output_format: "gbk", "fasta", or "both"
            validate: 검증 수행 여부
            **kwargs: 조립 방법별 추가 파라미터

        Returns:
            조립된 plasmid 서열 리스트

        Example:
            >>> workflow = ClassicalWorkflow(output_dir="results")
            >>> results = workflow.run(
            ...     insert_file="gfp.fasta",
            ...     vector_file="plasmid.gb",
            ...     method="gibson",
            ...     output_prefix="gfp_plasmid"
            ... )

        Notes:
            [P5 shim] 파일기반 진입을 build-plan(vector_insertion stage)으로 위임한다.
            load → run_plan([vector_insertion]) → construct. 출력은 constructs/{id}.gb,
            계보는 reports/manifest.json. (구 export(output_dir 평면 저장)는 stage가 대체.)
        """
        print("=" * 60)
        print(f"  Classical Cloning Workflow - {method.upper()}")
        print("=" * 60)
        print()

        # 1. 파일 읽기 → run_plan artifact registry seed
        print("Loading files...")
        insert, vector = self.load_sequences(insert_file, vector_file)
        print(f"  Insert: {len(insert.seq)} bp | Vector: {len(vector.seq)} bp\n")
        self.parts = [insert, vector]

        # 2. build-plan 실행 (vector_insertion stage)
        plan = [{
            "op": "vector_insertion",
            "inputs": [insert.name, vector.name],
            "params": {"method": method, "assembly_type": "circular"},
        }]
        outs = self.run_plan(plan)

        print()
        print("=" * 60)
        print("  Workflow complete")
        print("=" * 60)

        return [o["record"] for o in outs]


__all__ = ['ClassicalWorkflow']
