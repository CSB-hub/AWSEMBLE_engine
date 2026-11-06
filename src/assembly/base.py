"""
Assembly 추상 베이스 클래스

모든 DNA Assembly 방법(Gibson, Golden Gate 등)의 공통 인터페이스를 정의합니다.

핵심 개념:
    - prepare_fragments(): 각 Assembly 방법에 따라 다른 준비 과정
      * Gibson: 프라이머 설계 → Tuple[str, str] 반환
      * Golden Gate: Restriction cutting → List[Dseqrecord] 반환
    - assemble(): DNA 조립 수행
    - validate(): 결과 검증

작성일: 2025-11-06
"""

from abc import ABC, abstractmethod
from typing import List, Dict, Any
from pydna.dseqrecord import Dseqrecord


class AssemblyBase(ABC):
    """모든 Assembly 방법의 추상 베이스 클래스

    모든 Assembly 클래스는 이 인터페이스를 구현해야 함.

    Attributes:
        준비 과정, 조립, 검증 메서드를 추상 메서드로 정의

    작성일: 2025-11-06
    """

    @abstractmethod
    def prepare_fragments(self, **kwargs) -> Any:
        """Fragment 준비 (프라이머 설계, overlap 추가, cutting 등)

        각 Assembly 방법에 따라 다른 처리를 수행합니다.

        Args:
            **kwargs: 각 Assembly 방법별로 다른 파라미터
            # TODO: general 버전에 dependent 한 거 같음. 초반에 설계한 코드라서 그런듯.
            # biopart용과 호환되도록 개편해야하나?
                - Gibson: insert, vector, target_start, target_end
                - Golden Gate: dseqs (cutting할 DNA 리스트)

        Returns:
            각 Assembly 방법에 따라 다름:
            # TODO: prepare_fragments 함수의 args와 동일한 고민
                - Gibson: Tuple[str, str] (forward_primer, reverse_primer)
                - Golden Gate: List[Dseqrecord] (cut fragments)

        Notes:
            - Gibson은 프라이머 설계 후 사용자가 PCR 수행
            - Golden Gate는 Restriction cutting으로 Fragment 준비
            - 반환 타입이 다른 이유: 근본적으로 다른 워크플로우

        Example:
            >>> # Gibson
            >>> gibson = GibsonCloning()
            >>> fwd, rev = gibson.prepare_fragments(insert=my_insert, vector=my_vector)
            >>>
            >>> # Golden Gate
            >>> gg = GoldenGateModule()
            >>> fragments = gg.prepare_fragments(dseqs=[vector, insert])
        """
        pass

    @abstractmethod
    def assemble(self, fragments: List[Dseqrecord], **kwargs) -> List[Dseqrecord]:
        """Assembly 수행

        Args:
            fragments: 조립할 fragment 리스트
            **kwargs: 조립 옵션
                # TODO: golden은?
                - assembly_type: "circular" 또는 "linear" (Gibson)
                - 기타 방법별 옵션

        Returns:
            조립된 결과물 리스트 (여러 후보 가능)

        Raises:
            ValueError: 조립 실패 시

        Example:
            >>> # Gibson circular assembly
            >>> results = gibson.assemble([vector, insert], assembly_type="circular")
            >>> results[0].circular
            True
            >>>
            >>> # Golden Gate linear assembly
            >>> results = gg.assemble(fragments)
            >>> len(results) > 0
            True
        """
        pass

    @abstractmethod
    def validate(self, result: Dseqrecord) -> Dict:
        """Assembly 결과 검증

        Args:
            result: 조립 결과물

        Returns:
            검증 결과 Dict:
            {
                'valid': bool,           # 전체 검증 통과 여부
                'warnings': List[str],   # 경고 메시지
                'errors': List[str],     # 에러 메시지
                'details': Dict          # 상세 정보 (길이, features 등)
            }

        Example:
            >>> validation = gibson.validate(result)
            >>> if not validation['valid']:
            ...     print("에러:", validation['errors'])
            >>> else:
            ...     if validation['warnings']:
            ...         print("경고:", validation['warnings'])
        """
        pass


# 모듈 메타데이터
__all__ = ['AssemblyBase']
