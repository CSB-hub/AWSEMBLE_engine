"""
awemble.assembly - DNA Assembly 메서드 패키지

DNA 조립을 위한 다양한 방법(Gibson, Golden Gate 등)을 제공합니다.

핵심 컴포넌트:
    - AssemblyBase: 모든 Assembly 방법의 추상 베이스 클래스
    - GibsonCloning: Homology-based overlap 조립 (15-25bp)
    - GoldenGateModule: Type IIS 제한효소 기반 조립 (BsaI, 4bp overhang)
    - AssemblyFactory: Assembly 인스턴스 생성 팩토리
    - 공통 유틸리티: remove_duplicate_features, adjust_feature_positions

사용 예제:
    >>> from awsemble.assembly import AssemblyFactory
    >>>
    >>> # Gibson Assembly
    >>> gibson = AssemblyFactory.create_assembler("gibson", overlap_length=20)
    >>> fwd, rev = gibson.prepare_fragments(insert, vector)
    >>> results = gibson.assemble([vector, insert], assembly_type="circular")
    >>>
    >>> # Golden Gate Assembly
    >>> gg = AssemblyFactory.create_assembler("golden_gate")
    >>> fragments = gg.prepare_fragments([part1, part2, part3])
    >>> results = gg.assemble(fragments, assembly_type="linear")

변경 사항 (v2.0.0):
    - prepare_fragments() 반환 타입 통일 제거 (Gibson: Tuple[str, str], Golden Gate: List[Dseqrecord])
    - calculate_overlap_region() 제거 → gibson.py 내부로 통합
    - create_directories() 제거 → Workflow 모듈로 이동
    - Vector linearization 자동화 (Gibson)
    - BsaI cutting 자동화 (Golden Gate)
    - 프라이머 특이성 검증 추가

작성일: 2025-11-06 (v2.0.0 리팩토링)
"""

# 베이스 클래스
from .base import AssemblyBase

# 조립 방법 구현
from .gibson import GibsonCloning
from .golden_gate import GoldenGateModule

# Factory 패턴
from .factory import AssemblyFactory

# 공통 유틸리티
from .common import (
    remove_duplicate_features,
    adjust_feature_positions,
)

# Public API
__all__ = [
    # 베이스 클래스
    'AssemblyBase',

    # Assembly 방법
    'GibsonCloning',
    'GoldenGateModule',

    # 팩토리
    'AssemblyFactory',

    # 유틸리티
    'remove_duplicate_features',
    'adjust_feature_positions',
]

__version__ = '2.0.0'
__author__ = 'awemble team'
