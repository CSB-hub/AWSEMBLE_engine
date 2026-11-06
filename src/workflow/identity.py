"""Workflow artifact id 체계 + manifest (J-2)

평평한 고정 버킷(parts/modules/multimodules/...)에는 계보 정보가 없으므로,
reports/manifest.json이 artifact 간 DAG(부모 id)를 보유하는 **진실의 원천**이다.

id 체계 (모드별):
    - Biopart: 파생 artifact(prepared fragment / module / multimodule) id =
      안정입력(part DB id + method + overhang)의 **결정적 짧은 해시**.
      같은 입력 → 같은 id → 재현성/계보 추적 가능.
    - Classical: 시퀀셜 temp id(`part_0001` ...). DB anchor가 없으므로 해시 불가.

핵심 조립 로직(prepare_fragments/assemble)은 이 모듈에 의존하지 않는다.
id/manifest는 오케스트레이션 계층이 조립 결과에 **덧붙이는** 부가 산출물이다.

작성일: 2026-06-12 (P4-J2)
"""

import hashlib
import json
from pathlib import Path
from typing import Dict, List, Any, Optional, Sequence

# 레벨별 id 접두사 (승인 스펙 3, + P5 construct)
LEVEL_PREFIX = {
    "prepared_fragment": "pf",
    "module": "mod",
    "multimodule": "mm",
    "construct": "con",   # vector_insertion 결과 (insert+vector→circular)
}

_HASH_LEN = 8       # sha256[:8]
_DELIM = "|"        # 정규화 직렬화 구분자


def normalize_inputs(stable_inputs: Sequence[Any]) -> str:
    """안정입력을 결정적 문자열로 직렬화.

    part 순서는 조립상 의미가 있으므로 **보존**한다(정렬하지 않음).
    각 항목을 str로 변환 후 고정 구분자로 연결해, 같은 논리적 입력이
    항상 동일한 문자열이 되도록 한다.

    Args:
        stable_inputs: 정규화할 안정입력 리스트
            (예: [part_db_id..., method, overhang...] — 호출자가 필드 순서 고정)

    Returns:
        직렬화된 정규화 문자열
    """
    return _DELIM.join(str(x) for x in stable_inputs)


def derived_id(level: str, stable_inputs: Sequence[Any]) -> str:
    """결정적 짧은 해시 id 생성 (Biopart 파생 artifact용).

    Args:
        level: artifact 레벨 ("prepared_fragment" | "module" | "multimodule")
        stable_inputs: 안정입력 리스트(part DB id + method + overhang 등).
            part 순서는 보존됨(normalize_inputs 참조).

    Returns:
        `<prefix>_<hash8>` 형식 id (예: "mod_a1b2c3d4")

    Raises:
        KeyError: 알 수 없는 level
    """
    prefix = LEVEL_PREFIX[level]
    norm = normalize_inputs(stable_inputs)
    digest = hashlib.sha256(norm.encode("utf-8")).hexdigest()[:_HASH_LEN]
    return f"{prefix}_{digest}"


class TempIdAllocator:
    """Classical 모드용 시퀀셜 temp id 발급기 (`part_0001` ...).

    DB id 같은 안정 anchor가 없는 Classical 경로에서, 접두사별로
    1부터 증가하는 4자리 zero-pad id를 발급한다.

    Example:
        >>> alloc = TempIdAllocator()
        >>> alloc.next("part")
        'part_0001'
        >>> alloc.next("part")
        'part_0002'
        >>> alloc.next("mod")
        'mod_0001'
    """

    def __init__(self) -> None:
        self._counters: Dict[str, int] = {}

    def next(self, prefix: str) -> str:
        """prefix별 다음 시퀀셜 id 반환."""
        n = self._counters.get(prefix, 0) + 1
        self._counters[prefix] = n
        return f"{prefix}_{n:04d}"


class Manifest:
    """reports/manifest.json — artifact 계보(DAG)의 진실의 원천.

    artifact별 1엔트리. 필드(승인 스펙 2): id, level, file, op, method, inputs, params.
    평평한 버킷엔 계보가 없으므로 inputs(부모 id)가 DAG를 형성한다.

    Args:
        project_dir: 프로젝트 디렉토리 (manifest는 그 아래 reports/manifest.json)

    Example:
        >>> m = Manifest(project_dir)
        >>> m.add(id="mod_a1b2c3d4", level="module",
        ...       file="modules/assembled/mod_a1b2c3d4.gb",
        ...       op="module_assembly", method="golden_gate",
        ...       inputs=["BBa_B0034", "BBa_C0062"],
        ...       params={"assembly_type": "linear"})
        >>> m.write()
    """

    FILENAME = "manifest.json"

    def __init__(self, project_dir) -> None:
        self.project_dir = Path(project_dir)
        self.path = self.project_dir / "reports" / self.FILENAME
        self.entries: List[Dict[str, Any]] = []

    def add(
        self,
        *,
        id: str,
        level: str,
        file: str,
        op: str,
        method: str,
        inputs: Sequence[str],
        params: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """artifact 엔트리 1건 추가.

        Args:
            id: artifact id (해시 또는 temp)
            level: artifact 레벨 (part/module/multimodule/construct ...)
            file: project_dir 기준 상대 경로
            op: 생성 연산 (예: "module_assembly")
            method: 조립 방법 (golden_gate/gibson). 해당 없으면 "" 가능
            inputs: 부모 artifact id 리스트 (계보 DAG의 간선)
            params: 연산 파라미터 (assembly_type 등)

        Returns:
            추가된 엔트리 dict
        """
        entry = {
            "id": id,
            "level": level,
            "file": file,
            "op": op,
            "method": method,
            "inputs": list(inputs),
            "params": params or {},
        }
        self.entries.append(entry)
        return entry

    def write(self) -> Path:
        """reports/manifest.json으로 직렬화 저장.

        Returns:
            저장된 manifest.json 경로
        """
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.path, "w", encoding="utf-8") as f:
            json.dump({"artifacts": self.entries}, f, indent=2, ensure_ascii=False)
        return self.path


__all__ = [
    "LEVEL_PREFIX",
    "normalize_inputs",
    "derived_id",
    "TempIdAllocator",
    "Manifest",
]
