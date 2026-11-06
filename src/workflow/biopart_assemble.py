"""
Biopart 조립 워크플로우 (구 DBTL — Design-Build-Test-Learn 방법론 기반)

KRIBB 6단계 조합 조립 워크플로우를 제공합니다.

핵심 기능:
    - Excel 설계 파일 → GenBank 변환 (Step 1)
    - Part-Vector 클로닝 (Step 2)
    - Module 조립 (Step 3)
    - Module-Vector 클로닝 (Step 4)
    - Multi-module 조립 (Step 5)
    - Pathway 조립 (Step 6)
    - 각 단계별 조립 방법 선택 가능 (Gibson/Golden Gate)

사용 예제:
    >>> workflow = BiopartWorkflow(
    ...     project_dir="haseong_240704",
    ...     design_file="assembly_design.xlsx",
    ...     vector_file="addgene-plasmid-50005.gb"
    ... )
    >>> workflow.run_full_workflow(
    ...     step2_method="gibson",
    ...     step3_method="golden_gate",
    ...     step4_method="gibson",
    ...     step5_method="gibson",
    ...     step6_method="gibson"
    ... )

작성일: 2025-01-13
"""

from pathlib import Path 
from typing import List, Dict, Any, Optional
import copy
import json
import pandas as pd
from pydna.dseqrecord import Dseqrecord

from workflow.base import WorkflowBase
from workflow.identity import Manifest, derived_id
from assembly.factory import AssemblyFactory
from assembly.common import remove_duplicate_features
from utils.io import (
    DNAFileReader, 
    DNAFileWriter, 
    DatabaseImport
)

# TODO: builder와 base했던 작성했던 내용들 중복

def _params_key(design_params) -> str:
    """artifact id 해시 재료로 쓸 설계 파라미터 문자열.

    None이면 빈 문자열 — **기존(파라미터 도입 이전) id를 그대로 유지**하기 위해서다.
    (파라미터를 안 쓰는 호출은 예전과 같은 id/파일명을 갖는다 = 회귀 없음)
    """
    return design_params.canonical() if design_params is not None else ""


def _params_record(design_params) -> dict:
    """manifest에 남길 설계 파라미터(사람이 읽고 재현하는 용도)."""
    return {"design_params": design_params.to_dict()} if design_params is not None else {}


class BiopartWorkflow(WorkflowBase):
    """
    Biopart 조립 워크플로우 (구 DBTL — Design-Build-Test-Learn 방법론)

    KRIBB 프로젝트의 6단계 조합 조립 워크플로우를 제공합니다.
    각 단계별로 독립 실행 가능하며, 조립 방법을 자유롭게 선택할 수 있습니다.

    워크플로우 구조:
        Step 1: Excel → GenBank (Part 설계)
        Step 2: Part + Vector → Part-Vector (Gibson/Golden Gate)
        Step 3: Parts → Module (Gibson/Golden Gate)
        Step 4: Module + Vector → Module-Vector (Gibson/Golden Gate)
        Step 5: Modules → Multi-module (Gibson/Golden Gate)
        Step 6: Multi-modules → Pathway (Gibson/Golden Gate)

    Args:
        project_dir: 프로젝트 디렉토리 경로
        vector_file: Vector GenBank 파일명

    Attributes:
        project_dir: 프로젝트 디렉토리 Path 객체
        vector_file: 벡터 파일 경로
        common: 공통 유틸리티 모듈
        part_prep: Part 준비 모듈

    Example:
        >>> workflow = BiopartWorkflow(
        ...     project_dir="my_project",
        ...     vector_file="vector.gb"
        ... )
        >>> workflow.run_step1_part_design()
        >>> workflow.run_step2_part_cloning(method="gibson")

    작성일: 2025-01-13
    """

    def __init__(
        self,
        project_dir: str,
        parts: List[Dseqrecord],
        vector_file: Optional[str] = None
    ):
        """
        BiopartWorkflow 초기화

        Args:
            project_dir: 프로젝트 디렉토리 경로
            parts: Dseqrecord 객체들의 리스트
            vector_file: Vector GenBank 파일명 (선택 사항)
        """
        super().__init__(project_dir)
        self.parts = parts
        self.vector_file = None
        if vector_file:
            self.vector_file = self.project_dir / vector_file
            # 필수 파일 존재 확인
            if not self.vector_file.exists():
                raise FileNotFoundError(f"벡터 파일을 찾을 수 없습니다: {self.vector_file}")

        # 디렉토리 구조 생성
        self._setup_directories()

    def _setup_directories(self) -> None:
        """J-2 고정 표준 레이아웃 생성 (WorkflowBase.create_standard_layout 위임).

        parts/modules/multimodules/constructs/reports/logs + {method}_prepared/assembled.
        (구 parts/insert·modules/inserts·pathways 구조 → 표준 버킷으로 교체)
        """
        self.create_standard_layout()

    def run_step1_part_design(self) -> None:
        """[DEPRECATED no-op shim — P4-J3] Step 1: Part Design.

        Part 준비는 워크플로우 앞단(DB 조회 + prepare_fragments)에서 끝나므로 내부 구현 없음.
        하네스(test_single/test_all/test_kribb) 호환용으로 남겨둔 no-op.
        하네스를 build-plan API로 이전하면 제거 예정.
        """
        pass

    # ------------------------------------------------------------------ #
    # J-3: 조합형 stage (run_plan 엔진이 dispatch). run_step1~6 폐기 → stage.
    #   - _assemble_fragments: method별 조립 핵심 (불변 로직)
    #   - _save_artifact: 저장 + manifest 기록 공통화
    #   - module_assembly / multimodule_assembly: 레벨만 다른 조합 연산
    # ------------------------------------------------------------------ #

    def _assemble_fragments(
        self,
        parts: List[Dseqrecord],
        method: str,
        assembly_type: str,
        design_params=None,
        vector: Optional[Dseqrecord] = None,
    ) -> Optional[Dseqrecord]:
        """조립 핵심 로직 (method별 fragment 처리 + assemble).

        핵심 조립 로직은 **불변**(구 run_step3 본문 그대로). golden_gate는
        각 부품에서 유사도 기반으로 insert fragment를 추출한 뒤 결합하고,
        그 외(gibson 등)는 부품을 바로 결합한다.

        Args:
            parts: 입력 Dseqrecord 리스트 (part 또는 하위 module/multimodule)
            method: 조립 방법 ("golden_gate" | "gibson" ...)
            assembly_type: 조립 형태 ("linear" | "circular")

        Returns:
            조립된 Dseqrecord. 추출/조립 실패 시 None.
        """
        # 설계 파라미터(연구실 조건)를 assembler에 전달.
        # golden_gate는 params.gg_enzyme으로 효소까지 결정된다(부품 내부 사이트 회피는 그 안에서).
        akw = {}
        if design_params is not None:
            akw["params"] = design_params
            if method.lower() == "golden_gate":
                akw["enzyme"] = design_params.gg_enzyme
        assembler = AssemblyFactory.create_assembler(method, **akw)

        if method.lower() == "golden_gate":
            # 각 부품에서 assemble에 사용할 fragment를 추출합니다.
            # (GoldenGateModule 내부의 유사도 기반 추출 로직 사용)
            fragments = []
            for dseq in parts:
                best_frag = assembler.extract_insert_fragment(dseq)
                if best_frag:
                    fragments.append(best_frag)

            if len(fragments) != len(parts):
                print(f"  ⚠️  Fragment 추출 실패. {len(parts)}개 중 {len(fragments)}개만 확보됨. 건너뜀.\n")
                return None
            target = fragments
            if assembly_type == "circular":
                # circular은 벡터 backbone이 조립에 참여한다(중간 클로닝).
                akw2 = {"vector": vector}
            else:
                akw2 = {}
        else:  # gibson 등 다른 조립 방법
            target = parts

        # fragment가 하나면 **이을 것이 없다** — 그 조각이 곧 module이다.
        # Gibson only(부품 1개)에서 양 끝 vector homology가 프라이머 tail로 흡수되면
        # module 전체가 단일 PCR 산물이 되고, 조립 상대는 선형화된 벡터다(삽입 단계).
        # pydna에 조각 하나를 넘기면 겹칠 상대가 없어 빈 결과가 나오므로 조립만 건너뛴다.
        #
        # ⚠️ **여기서 함수를 빠져나가면 안 된다.** 2026-11-04에 `return target[0]`으로
        #    끊었더니 아래 **P6 feature 기입까지 통째로 건너뛰어** module/construct에
        #    파츠 feature가 하나도 안 남았다(실측: 부품 4개는 8개, 부품 1개는 0개).
        #    그 결과 벡터맵에서 파츠를 눌러도 자리를 못 잡았다("not in the case you are
        #    viewing"). 조립은 건너뛰되 **후처리는 똑같이 태운다.**
        if len(target) == 1 and assembly_type == "linear":
            print(f"  · fragment 1개({target[0].name}) — 이을 접합부가 없어 그대로 module로 씁니다")
            # 조립 경로의 `results[0]`는 **새 객체**다. 여기서 입력을 그대로 쓰면 아래
            # feature 기입이 prepared part 레코드를 직접 건드리게 되므로 사본을 쓴다.
            record = copy.deepcopy(target[0])
        else:
            try:
                results = assembler.assemble(target, assembly_type=assembly_type,
                                             **(akw2 if method.lower() == "golden_gate" else {}))
            except Exception as e:
                print(f"  ⚠️  조립 실패: {e}\n")
                return None

            if not results:
                print("  ⚠️  조립 실패\n")
                return None

            record = results[0]

        # 조립 산물에 효소 사이트가 남았는지 알린다(차단 X).
        # 접합부에서 인식서열이 새로 생기는 경우가 실제로 있다(vector overhang + 인접 부품).
        if method.lower() == "golden_gate":
            try:
                assembler.check_residual_sites(record, "module")
            except Exception:
                pass

        # P6: module에 part 본체 + 접합부 feature 기입 (조립 불변 — 후처리, 서열 불변)
        #
        # ⚠️ 경로마다 좌표를 얻는 방법이 다르다. 예전에는 이 블록 전체가
        #    `method == "golden_gate"` 하나로 막혀 있어서 homology 경로(gibson·yeast_hr)는
        #    feature가 **0개**로 나왔다. 그 결과 올리고 주문 시트(module feature에서
        #    부품을 뽑는다)가 통째로 비고, 워크플로우 창에서 파츠를 눌러도 벡터맵이
        #    자리를 못 잡았다(2026-10-25).
        if method.lower() in ("gibson", "yeast_hr"):
            # homology 경로: 접합부 길이가 자리마다 다르고 한 fragment가 부품 여러 개를
            # 담을 수 있어 exact match로는 못 되찾는다. 대신 JunctionPlan이 실어 보낸
            # 배치(part_spans/body_len)로 **좌표를 계산**한다.
            try:
                from assembly.common import annotate_module_parts_from_layout
                annotate_module_parts_from_layout(record, parts)
            except Exception as e:
                print(f"  ⚠️  part feature 기입 생략: {e}")
        elif method.lower() == "golden_gate":
            try:
                import re
                from assembly.common import annotate_module_parts
                # 부품 id/type은 annotations로 전파됨(io.get_part_by_id→prep res→records).
                # 없으면(옛 경로) prep 네이밍("<part>_<i>_designed"+하네스 재부착)에서 정리.
                names, types = [], []
                for p in parts:
                    pid = p.annotations.get("part_id") or re.sub(r"(_\d+)?(_designed)+$", "", p.name)
                    names.append(pid)
                    types.append(p.annotations.get("part_type", ""))
                # overhang 길이는 조립에 쓴 효소가 정한다(BsaI 4 / SapI 3) — 4 고정 금지.
                annotate_module_parts(
                    record, fragments, names, types,
                    overhang_len=getattr(assembler, "overhang_len", 4),
                )
            except Exception as e:
                print(f"  ⚠️  part feature 기입 생략: {e}")

        return record

    def _save_artifact(
        self,
        record: Dseqrecord,
        *,
        level: str,
        bucket: str,
        artifact_id: str,
        op: str,
        method: str,
        inputs: List[Dseqrecord],
        assembly_type: str,
        manifest: Optional[Manifest] = None,
        design_params=None,
    ) -> Dict[str, Any]:
        """조립 결과를 버킷에 저장하고 manifest에 계보를 기록.

        Args:
            record: 저장할 조립 결과 Dseqrecord
            level: artifact 레벨 ("module" | "multimodule")
            bucket: project_dir 기준 출력 버킷 (예: "modules/assembled")
            artifact_id: artifact id (파일명 stem이 됨 → id == 파일 stem)
            op: 생성 연산 (manifest op)
            method: 조립 방법
            inputs: 입력 artifact Dseqrecord 리스트 (manifest inputs=이름)
            assembly_type: 조립 형태
            manifest: 계보 기록 대상 (None이면 기록 생략)

        Returns:
            {"id", "file", "record"} dict
        """
        # 후속 stage가 이 record를 입력으로 받을 때 id를 알 수 있도록 이름에 심는다.
        # (구버전은 name이 "assembled_module" 고정이라 build plan이 예측한 id와 어긋나
        #  "입력 artifact를 찾을 수 없습니다"로 실패했다)
        record.name = str(artifact_id)

        out_dir = self.project_dir / bucket
        out_dir.mkdir(parents=True, exist_ok=True)
        output_path = out_dir / f"{artifact_id}.gb"
        DNAFileWriter.write_genbank(record, str(output_path))

        if manifest is not None:
            manifest.add(
                id=artifact_id,
                level=level,
                file=f"{bucket}/{artifact_id}.gb",
                op=op,
                method=method,
                inputs=[p.name for p in inputs],
                # 해시는 충돌을 막고, manifest는 "무슨 값으로 만들었는지"를 읽게 한다 — 둘 다 필요.
                params={"assembly_type": assembly_type, **_params_record(design_params)},
            )

        return {"id": artifact_id, "file": output_path, "record": record}

    def module_assembly(
        self,
        parts: List[Dseqrecord],
        method: str = "golden_gate",
        assembly_type: str = "linear",
        module_id: Optional[str] = None,
        manifest: Optional[Manifest] = None,
        design_params=None,
    ) -> Optional[Dict[str, Any]]:
        """stage: N parts → 1 module.

        현재 구현의 핵심 조합 연산(구 run_step3). 출력 modules/assembled/{id}.gb.

        Args:
            parts: 조립할 prepared part Dseqrecord 리스트 (순서 의미 있음)
            method: 조립 방법
            assembly_type: 조립 형태
            module_id: 명시 id (shim/하네스=test_id). None이면 결정적 해시 mod_<hash8>.
            manifest: 계보 기록 대상

        Returns:
            {"id", "file", "record"} 또는 실패 시 None.
        """
        # circular(중간 클로닝)은 마지막 record를 중간 벡터로 받는다.
        # (vector_insertion의 records=[insert, vector] 관례와 동일)
        vector = None
        if assembly_type == "circular":
            if len(parts) < 2:
                print("  ⚠️  circular 조립에는 부품 + 중간 벡터가 필요합니다.")
                return None
            parts, vector = list(parts[:-1]), parts[-1]

        record = self._assemble_fragments(parts, method, assembly_type, design_params, vector=vector)
        if record is None:
            return None
        # biopart module은 vector에 **통째로** 들어간다. 표시가 없으면 vector_insertion의
        # CDS 자동감지(_find_target_feature)가 module 안의 CDS 부품을 삽입 대상으로 오인해
        # 그 구간만 증폭한다(= promoter/RBS/terminator/spacer 유실).
        record.annotations["insert_whole"] = "1"

        if module_id is None:
            # build-plan 경로: 안정입력의 결정적 해시 (part 이름 + method + 형태)
            # ⚠️ 설계 파라미터를 해시 재료에 넣지 않으면 Tm 65로 만든 산출물과 Tm 72로 만든
            #    산출물이 **같은 id**를 갖고 파일({id}.gb)을 서로 덮어쓴다.
            #    (일부 파라미터만 골라 넣으면 누락 시 조용한 덮어쓰기가 되므로 전부 넣는다.
            #     무관한 파라미터가 바뀌어도 id가 갈리는 건 파일이 늘 뿐이라 안전한 쪽이다.)
            module_id = derived_id(
                "module",
                [p.name for p in parts]
                + ([vector.name] if vector is not None else [])
                + [method, assembly_type, _params_key(design_params)],
            )

        return self._save_artifact(
            record,
            level="module",
            bucket="modules/assembled",
            artifact_id=str(module_id),
            op="module_assembly",
            method=method,
            inputs=parts,
            assembly_type=assembly_type,
            manifest=manifest,
            design_params=design_params,
        )

    def module_excision(
        self,
        records: List[Dseqrecord],
        method: str = "digest",
        assembly_type: str = "linear",
        exc_id: Optional[str] = None,
        manifest: Optional[Manifest] = None,
        design_params=None,
    ) -> Optional[Dict[str, Any]]:
        """stage: circular module(중간 벡터에 든 것) → 꺼낸 linear module.

        **왜 필요한가**: 중간 클로닝에서 module은 벡터 안에 원형으로 들어가 있다. 다음 단계
        (module끼리 Gibson 연결 / 최종 벡터 삽입)에 쓰려면 벡터에서 꺼내야 한다.
        Type IIS 조립은 인식서열을 없애므로 **조립 효소로는 되꺼낼 수 없고**, 중간 벡터가
        제공하는 제2 효소로 자르거나(digest) spacer에 프라이머를 앉혀 증폭한다(pcr).

        꺼낸 조각에는 양 끝 spacer가 딸려 나와야 한다 — 그 spacer가 다음 단계의 homology다.
        (중간 벡터에서 excision 사이트를 spacer **바깥**에 두는 이유)

        Args:
            records: [circular module]
            method: 미사용(stage 계약 호환). 실제 방식은 design_params.excision_method
            assembly_type: 미사용(stage 계약 호환)
            exc_id: 명시 id. None이면 결정적 해시
            manifest: 계보 기록 대상
            design_params: 설계 파라미터

        Returns:
            {"id","file","record"} 또는 실패 시 None
        """
        from params import DEFAULT_PARAMS

        prm = design_params if design_params is not None else DEFAULT_PARAMS
        # records=[circular module, 중간 vector]. vector는 **어느 조각이 backbone인지** 가리는 데 쓴다.
        # (길이로 추측하면 backbone이 작은 벡터에서 뒤집힌다 — 실제로 그렇게 틀렸었다)
        module = records[0]
        vector = records[1] if len(records) > 1 else None
        # spacer를 부품으로 붙인 경우 그 spacer record들이 뒤에 온다 → PCR 앵커로 쓴다.
        spacers = list(records[2:]) if len(records) > 2 else []
        want = prm.excision_method

        # spacer를 부품으로 붙였다면 digest는 벡터 서열을 달고 나온다(모든 module이 같은 말단이 되어
        # 인접 module 순서를 구분할 수 없다). 이 경우 PCR이 맞다.
        if spacers and want == "digest":
            print("     ⚠️ spacer를 부품으로 붙인 구성입니다 — digest로 꺼내면 벡터 서열이 "
                  "말단에 남아 module 간 구분이 사라집니다. PCR(spacer 앵커)로 진행합니다.")
            want = "pcr"

        enzyme, why = self._pick_excision_enzyme(module, prm, vector)
        if want == "digest" and enzyme is None:
            # 실측상 10% 정도는 쓸 효소가 없다 → 조립을 죽이지 않고 pcr로 내린다.
            print(f"     ⚠️ digest용 효소를 찾지 못해 PCR 방식으로 전환합니다 ({why})")
            want = "pcr"

        if want == "digest":
            record = self._excise_by_digest(module, enzyme, vector)
            used = f"digest({enzyme})"
        else:
            record = self._excise_by_pcr(module, vector, prm, spacers)
            used = "pcr"
        if record is None:
            return None
        record.name = "excised_module"

        if exc_id is None:
            # id는 **요청한 방식**(prm.excision_method)으로 만든다. 실제 선택된 효소는
            # 자동 폴백으로 달라질 수 있어 build plan이 미리 알 수 없다.
            exc_id = derived_id(
                "module",
                [module.name, "excision", prm.excision_method, _params_key(design_params)],
            )
        return self._save_artifact(
            record,
            level="module",
            bucket="modules/excised",
            artifact_id=str(exc_id),
            op="module_excision",
            method=used,
            inputs=[module],
            assembly_type="linear",
            manifest=manifest,
            design_params=design_params,
        )

    def _pick_excision_enzyme(self, module: Dseqrecord, prm,
                              vector: Optional[Dseqrecord] = None):
        """module을 꺼낼 효소 선택.

        ⚠️ "module 안에 사이트가 2개"만 보면 안 된다 — 부품 내부에 우연히 2개가 있는 효소를
        벡터의 excision 효소로 오인해 엉뚱한 위치를 자른다(실제로 BBa_C0012의 BbsI 2개를
        그렇게 오인했다). **중간 벡터가 실제로 제공하는 효소**인지 먼저 확인한다.

        조건:
          (1) 조립 효소가 아니다 (Type IIS 조립은 인식서열을 없애므로 되꺼낼 수 없다)
          (2) **중간 벡터**에 사이트가 정확히 2개 (= 그 벡터가 제공하는 excision 효소)
          (3) **module**에도 정확히 2개 (부품 내부에 추가 사이트가 없다)

        Returns:
            (효소명 또는 None, 사유 문자열)
        """
        from Bio import Restriction
        from params import SUPPORTED_ENZYMES

        cands = [prm.excision_enzyme] if prm.excision_enzyme else [
            e for e in SUPPORTED_ENZYMES if e != prm.gg_enzyme
        ]
        counts = {}
        for name in cands:
            if not name:
                continue
            in_mod = len(getattr(Restriction, name).search(module.seq))
            in_vec = (len(getattr(Restriction, name).search(vector.seq))
                      if vector is not None else in_mod)
            counts[name] = f"vector {in_vec} / module {in_mod}"
            if in_vec == 2 and in_mod == 2:
                return name, f"{name} (벡터·module 모두 사이트 2개)"
        return None, f"쓸 수 있는 효소 없음 {counts}"

    @staticmethod
    def _backbone_marker(vector: Optional[Dseqrecord], asm_enzyme: str = "BsaI") -> str:
        """circular module 안에서 **벡터 backbone**을 식별할 표지 서열.

        ⚠️ 벡터 서열의 중앙을 쓰면 안 된다 — 거기는 stuffer라서 조립 때 제거되어
        module 어디에도 없다(그렇게 짰다가 표지를 못 찾아 module 전체를 반환했다).
        조립 효소로 자른 뒤 **남는 backbone**의 중앙을 써야 한다.

        Returns:
            표지 서열(30bp). 구할 수 없으면 빈 문자열.
        """
        if vector is None:
            return ""
        try:
            from assembly.golden_gate import GoldenGateModule
            backbone = GoldenGateModule(asm_enzyme)._get_backbone(vector)
        except Exception:
            return ""
        bb = str(backbone.seq).upper()
        if len(bb) < 60:
            return ""
        mid = len(bb) // 2
        return bb[mid - 15:mid + 15]

    def _excise_by_digest(self, module: Dseqrecord, enzyme_name: str,
                          vector: Optional[Dseqrecord] = None) -> Optional[Dseqrecord]:
        """excision 효소로 잘라 module 조각(벡터 backbone이 **아닌** 쪽)을 얻는다."""
        from Bio import Restriction

        try:
            frags = list(module.cut(getattr(Restriction, enzyme_name)))
        except Exception as e:
            print(f"  ⚠️  excision 절단 실패({enzyme_name}): {e}")
            return None
        if len(frags) < 2:
            print(f"  ⚠️  excision 절단 결과가 {len(frags)}조각 — 꺼낼 수 없습니다.")
            return None

        marker = self._backbone_marker(vector, self._asm_enzyme_for(vector))
        if marker:
            # backbone 표지를 **포함하지 않는** 조각이 module 쪽이다.
            cands = [f for f in frags if marker not in str(f.seq).upper()]
            if len(cands) == 1:
                module_frag = cands[0]
            else:
                module_frag = max(frags, key=lambda f: len(f.seq))
                print(f"  ⚠️  backbone 판별 모호({len(cands)}조각) — 가장 긴 조각을 module로 봅니다.")
        else:
            # 벡터를 못 받은 경우: ori/마커가 있는 backbone이 보통 더 길다는 가정.
            module_frag = min(frags, key=lambda f: len(f.seq))
            print("  ⚠️  vector 미지정 — 길이로 module 조각을 추정합니다(부정확할 수 있음).")

        print(f"  ✂️  excision({enzyme_name}): {len(frags)}조각 → module {len(module_frag.seq)}bp")
        return module_frag

    def _excise_by_pcr(self, module: Dseqrecord, vector: Optional[Dseqrecord],
                       prm, spacers: Optional[List[Dseqrecord]] = None) -> Optional[Dseqrecord]:
        """module 구간만 증폭해 꺼낸다 (부품 내부 제한효소 사이트 제약이 없는 방식).

        구간을 길이 산술로 구하면 안 된다 — spacer가 벡터 쪽에 속하는지 module 쪽에 속하는지가
        섞여 어긋난다(실제로 그렇게 틀렸다). **실제 프라이머가 앉을 자리**를 앵커로 쓴다:
        중간 벡터에서 excision 사이트 안쪽 구간의 양 끝(= L/R spacer의 시작·끝)을 뽑아,
        조립된 module에서 그 앵커 사이를 증폭 구간으로 잡는다. digest 결과와 같은 구간이 된다.
        """
        try:
            anchors = self._spacer_anchors(spacers) or self._module_anchors(vector, prm)
            seq = str(module.seq).upper()
            n = len(seq)
            if not anchors:
                print("  ⚠️  PCR excision: 앵커를 찾지 못해 module 전체를 증폭합니다.")
                region = seq
            else:
                fwd_anchor, rev_anchor = anchors
                doubled = seq * 2
                i = doubled.find(fwd_anchor)
                j = doubled.find(rev_anchor, i + len(fwd_anchor)) if i >= 0 else -1
                if i < 0 or j < 0 or (j + len(rev_anchor) - i) > n:
                    print("  ⚠️  PCR excision: 앵커 위치가 모호해 module 전체를 증폭합니다.")
                    region = seq
                else:
                    region = doubled[i: j + len(rev_anchor)]
            linear = Dseqrecord(region, name=module.name)
            print(f"  🧪 PCR excision: {len(region)}bp 증폭")
            return linear
        except Exception as e:
            print(f"  ⚠️  PCR 방식 excision 실패: {e}")
            return None

    @staticmethod
    def _spacer_anchors(spacers: Optional[List[Dseqrecord]], k: int = 25):
        """부품으로 붙인 spacer를 PCR 앵커로 쓴다.

        양 끝이 spacer가 되어야 인접 module과 **구분되는** Gibson homology가 된다.
        벡터 쪽 앵커를 쓰면 모든 module의 말단이 같아져 순서를 못 정한다.
        """
        if not spacers or len(spacers) < 2:
            return None
        l = str(spacers[0].seq).upper()
        r = str(spacers[-1].seq).upper()
        if len(l) < k or len(r) < k:
            return None
        return l[:k], r[-k:]

    def _module_anchors(self, vector: Optional[Dseqrecord], prm, k: int = 25):
        """중간 벡터에서 module 구간의 양 끝 앵커(프라이머가 앉을 자리)를 뽑는다.

        excision 효소로 벡터를 자르면 [module 자리(spacer·stuffer 포함)]와 [backbone]이 나온다.
        전자의 앞/뒤 k bp가 곧 L/R spacer의 시작·끝이며, 조립 후에도 module에 그대로 남는다.

        Returns:
            (fwd_anchor, rev_anchor) 또는 None
        """
        from Bio import Restriction

        if vector is None:
            return None
        enz, _ = self._pick_excision_enzyme(vector, prm)
        if enz is None:
            return None
        try:
            frags = list(vector.cut(getattr(Restriction, enz)))
        except Exception:
            return None
        if len(frags) < 2:
            return None
        marker = self._backbone_marker(vector, prm.gg_enzyme)
        cands = [f for f in frags if marker and marker not in str(f.seq).upper()]
        mod_side = cands[0] if len(cands) == 1 else max(frags, key=lambda f: len(f.seq))
        ms = str(mod_side.seq).upper()
        if len(ms) < 2 * k:
            return None
        return ms[:k], ms[-k:]

    @staticmethod
    def _vector_backbone_len(vector: Optional[Dseqrecord], asm_enzyme: str = "BsaI") -> int:
        """조립 효소로 자른 뒤 남는 backbone 길이."""
        if vector is None:
            return 0
        try:
            from assembly.golden_gate import GoldenGateModule
            return len(GoldenGateModule(asm_enzyme)._get_backbone(vector).seq)
        except Exception:
            return 0

    def _asm_enzyme_for(self, vector: Optional[Dseqrecord]) -> str:
        """digest 경로에서 backbone을 뽑을 때 쓸 조립 효소."""
        return getattr(self, "_last_asm_enzyme", "BsaI")

    def multimodule_assembly(
        self,
        modules: List[Dseqrecord],
        method: str = "golden_gate",
        assembly_type: str = "linear",
        mm_id: Optional[str] = None,
        manifest: Optional[Manifest] = None,
        design_params=None,
    ) -> Optional[Dict[str, Any]]:
        """stage: N (modules|multimodules) → 1 multimodule.

        module_assembly의 **재귀 일반화** — 입력 레벨만 다르고 조립 메커니즘은 동일
        (_assemble_fragments 공유). 재귀 입력(multimodule) 허용. 저위험 신규.
        출력 multimodules/assembled/{id}.gb.

        Args:
            modules: 조립할 module/multimodule Dseqrecord 리스트 (순서 의미 있음)
            method: 조립 방법
            assembly_type: 조립 형태
            mm_id: 명시 id. None이면 결정적 해시 mm_<hash8>.
            manifest: 계보 기록 대상

        Returns:
            {"id", "file", "record"} 또는 실패 시 None.
        """
        record = self._assemble_fragments(modules, method, assembly_type, design_params)
        if record is None:
            return None
        record.annotations["insert_whole"] = "1"   # module_assembly와 동일 이유

        if mm_id is None:
            mm_id = derived_id(
                "multimodule",
                [m.name for m in modules] + [method, assembly_type, _params_key(design_params)],
            )

        return self._save_artifact(
            record,
            level="multimodule",
            bucket="multimodules/assembled",
            artifact_id=str(mm_id),
            op="multimodule_assembly",
            method=method,
            inputs=modules,
            assembly_type=assembly_type,
            manifest=manifest,
            design_params=design_params,
        )

    def run_step3_module_assembly(self, module_df: pd.DataFrame, method: str = "golden_gate") -> None:
        """[DEPRECATED shim — P4-J3] 새 module_assembly stage 호출 래퍼.

        ⚠️ 테스트 하네스(test_all/test_single/verify_full) 호환용 얇은 래퍼다.
        조립/저장/manifest는 module_assembly stage가 담당(로직 0줄 위임) → 동작 불변.
        하네스를 build-plan API(run_plan)로 이전하면 **제거** 예정
        (260612_P4_design.md "추적: shim 제거" 참조).

        Args:
            module_df: 'id'(모듈 id) / 'name'(part 이름) 컬럼을 가진 DataFrame
            method: 조립 방법
        """
        try:
            if module_df.empty:
                print("✅ 모듈 정의가 제공되지 않아 모듈 조립을 건너뜁니다.\n")
                self.results['step3'] = {'method': method, 'module_count': 0, 'files': []}
                return

            if 'id' not in module_df.columns:
                raise ValueError("module_df DataFrame에 'id' 컬럼이 있어야 합니다.")

            module_df['id'] = module_df['id'].ffill()

            # 메모리에 저장된 self.parts를 이름으로 쉽게 찾을 수 있도록 딕셔너리로 변환
            parts_dict = {part.name: part for part in self.parts}
            manifest = Manifest(self.project_dir)  # J-2: 계보 기록(부가 산출물)
            saved_files = []

            for module_id in module_df['id'].unique():
                module_rows = module_df[module_df['id'] == module_id]
                part_names = [row['name'] for _, row in module_rows.iterrows() if pd.notna(row['name'])]

                if not part_names:
                    print(f"  ⚠️  {module_id}: Part 이름을 찾을 수 없습니다. 건너뜀.\n")
                    continue

                # part_names 순서에 따라 메모리에서 Dseqrecord 객체를 가져옴
                part_dseqs = [parts_dict.get(name) for name in part_names]
                if None in part_dseqs:
                    missing_part_index = part_dseqs.index(None)
                    missing_part_name = part_names[missing_part_index]
                    print(f"  ❌ Part '{missing_part_name}'을(를) 메모리에서 찾을 수 없습니다. 워크플로우를 중단합니다.\n")
                    continue

                # shim 경로: 호출자가 준 module_id(하네스=test_id)를 그대로 사용 → 파일/SHA 불변
                out = self.module_assembly(
                    part_dseqs, method=method, assembly_type="linear",
                    module_id=module_id, manifest=manifest,
                )
                if out:
                    saved_files.append(out["file"])

            manifest.write()
            out_dir = self.project_dir / "modules" / "assembled"
            self.results['step3'] = {'method': method, 'module_count': len(saved_files), 'output_dir': str(out_dir), 'files': [f.name for f in saved_files]}

        except Exception as e:
            print(f"\n❌ Step 3 실패: {e}")
            import traceback
            traceback.print_exc()
            raise

    # 구 run_step2/4/5/6(빈 골격)는 J-3에서 폐기 → 조합형 stage로 대체.
    #   step2(part_cloning)/step4(module_cloning) → vector_insertion (P7/Step4)
    #   step5(multimodule)/step6(pathway)         → multimodule_assembly
    # vector_insertion stage는 P7(Step4)에서 추가 예정.

    def run_full_workflow(self, plan: Optional[List[Dict[str, Any]]] = None, **kwargs) -> List[Dict[str, Any]]:
        """K-3: build plan 실행으로 전체 워크플로우 수행 (진입 일원화).

        구 run_step1~6 순차 호출 대신, build plan(stage 연산 리스트)을
        run_plan 엔진에 위임한다. 입력원(CLI/galaxy)은 build plan만 주입하면 된다.

        Args:
            plan: build plan. 각 항목 {"op", "inputs":[artifact_id...], "params":{...}}.
                  None이면 실행할 stage가 없으므로 빈 결과 반환.
            **kwargs: (구 stepN_method 호환 자리 — 현재 미사용)

        Returns:
            stage별 실행 결과 dict 리스트 (run_plan 반환)
        """
        if not plan:
            print("⚠️  build plan이 제공되지 않았습니다. 실행할 stage가 없습니다.")
            return []
        return self.run_plan(plan)

__all__ = ['BiopartWorkflow']