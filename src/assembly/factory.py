"""
Assembly Factory 패턴

Assembly 방법 선택을 문자열로 추상화하여 인스턴스를 생성합니다.

핵심 개념:
    - Factory 패턴으로 "gibson" 또는 "golden_gate" 문자열로 선택
    - 새로운 Assembly 방법 추가 시 한 곳에서만 수정 (개방-폐쇄 원칙)
    - 사용자는 구체 클래스 import 불필요

작성일: 2025-11-06
"""

from typing import Dict, Type, List
from .base import AssemblyBase

# Lazy import로 순환 참조 방지
# gibson, golden_gate는 나중에 import


class AssemblyFactory:
    """Assembly 인스턴스 생성 팩토리

    Factory 패턴으로 Assembly 방법 선택을 추상화합니다.

    Attributes:
        _assemblers: 등록된 Assembly 클래스 딕셔너리

    작성일: 2025-11-06
    """
    # TODO: 어떤 목적과 어떤 의미? line:31-33
    # Lazy loading으로 순환 참조 방지
    _assemblers: Dict[str, Type[AssemblyBase]] = {}

    @classmethod
    def _register_assemblers(cls):
        """Assembly 클래스 등록 (lazy loading)"""
        if not cls._assemblers:
            from .gibson import GibsonCloning
            from .golden_gate import GoldenGateModule
            from .yeast_hr import YeastHRCloning

            cls._assemblers = {
                'gibson': GibsonCloning,
                'golden_gate': GoldenGateModule,
                # 2026-10-24 구현 완료(구 스텁). Gibson 상속 + yeast 고유 3가지:
                #   arm 축(params.yeast_homology_min/_max 40~60bp) / 접합부 교차 차단 / in vivo 전제.
                #   GUI(툴 XML) 노출은 별도 트랙. 자세한 내용은 yeast_hr.py 모듈 docstring.
                'yeast_hr': YeastHRCloning,
            }

    @staticmethod
    def create_assembler(method_type: str, **kwargs) -> AssemblyBase:
        # TODO: primer_design.py와 연계가 잘 되어 있나?
        """Assembly 인스턴스 생성

        Args:
            method_type: "gibson" 또는 "golden_gate"
            **kwargs: Assembler 초기화 파라미터
                - Gibson: overlap_length (기본값: 20)
                - Golden Gate: enzyme (기본값: "BsaI")

        Returns:
            AssemblyBase 인스턴스

        Raises:
            ValueError: 지원하지 않는 method_type

        Example:
            >>> # Gibson Assembly
            >>> gibson = AssemblyFactory.create_assembler("gibson")
            >>> gibson = AssemblyFactory.create_assembler("gibson", overlap_length=15)
            >>>
            >>> # Golden Gate Assembly
            >>> gg = AssemblyFactory.create_assembler("golden_gate")
            >>> gg = AssemblyFactory.create_assembler("golden_gate", enzyme="BbsI")
        """
        # Lazy loading
        AssemblyFactory._register_assemblers()

        method_type = method_type.lower()

        if method_type not in AssemblyFactory._assemblers:
            available = ', '.join(AssemblyFactory._assemblers.keys())
            raise ValueError(
                f"Unknown assembly method: '{method_type}'. "
                f"Available methods: {available}"
            )

        assembler_class = AssemblyFactory._assemblers[method_type]
        return assembler_class(**kwargs)

    @staticmethod
    def get_available_methods() -> List[str]:
        """사용 가능한 Assembly 방법 목록 반환

        Returns:
            Assembly 방법 문자열 리스트

        Example:
            >>> methods = AssemblyFactory.get_available_methods()
            >>> print(methods)
            ['gibson', 'golden_gate']
        """
        # Lazy loading
        AssemblyFactory._register_assemblers()

        return list(AssemblyFactory._assemblers.keys())


# 모듈 메타데이터
__all__ = ['AssemblyFactory']