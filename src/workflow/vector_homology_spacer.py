"""module 양 끝을 **벡터 선형화 자리의 서열**로 만든다 — 삽입 PCR이 필요 없어진다.

## 기존 방식과 무엇이 다른가

기존(`SpacerSelector`)은 양 끝에 **vector 비종속 중립 서열**을 붙였다. 그 spacer는
삽입 프라이머가 anneal할 발판일 뿐이고, 벡터와의 homology는 프라이머의 5' tail이
따로 들고 왔다. 그래서:

  · spacer 115개 중 L/R 조합을 탐색해야 했다(실측 최대 380회 dry-run)
  · 중립 spacer 57bp × 2 = **114bp가 최종 플라스미드에 흉터로 남았다**
  · 삽입용 프라이머 2개를 따로 주문해야 했다

여기서는 spacer 자리에 **벡터 절단부 양옆 서열을 그대로** 넣는다. 그러면 module이
homology를 이미 들고 있으므로:

  · 조합 탐색이 사라진다 — (vector, 효소)가 정해지면 spacer도 정해진다
  · Gibson에서 homology는 벡터 서열과 **병합**되므로 흉터가 남지 않는다
  · 삽입 PCR이 불필요하다

## 어느 서열을 어느 쪽에 붙이나

`GibsonCloning.prepare_fragments`의 규약을 그대로 따른다. 거기서 증폭된 insert는

    vector[-overlap:] + insert + vector[:overlap]

가 되므로(fwd tail = 벡터 3' 끝, rev tail = 벡터 5' 끝의 역상보), module의
**왼쪽 끝 = 선형 벡터의 3' 끝**, **오른쪽 끝 = 선형 벡터의 5' 끝**이면 같은 산물이 된다.
"""
from typing import Optional, Tuple

from pydna.dseqrecord import Dseqrecord


class VectorHomologySpacer:
    """선형화된 벡터의 양 끝에서 homology arm을 떼어 module용 spacer로 만든다."""

    def __init__(self, params=None, arm: Optional[int] = None):
        """
        Args:
            params: 설계 파라미터 한 벌. None이면 `DEFAULT_PARAMS`.
            arm: 떼어낼 homology 길이(bp). None이면 `params.gibson_overlap`.

                ⚠️ **삽입에 쓸 조립 방법의 arm을 넘겨야 한다.** 이 spacer는 module 서열에
                그대로 박혀서 나중에 삽입 조립의 overlap이 된다. 삽입 방법이 다른 축을
                쓰면(yeast HR은 `yeast_homology_min` 40bp) 길이가 어긋나 **조립이 조용히
                실패**한다 — pydna가 요구 길이만큼 겹치는 구간을 못 찾는다.
                실제로 2026-10-24에 이 경로로 터졌다: module 양 끝 30bp,
                yeast 삽입 요구 45bp → "Assembly 실패: 결과물이 생성되지 않았습니다".
        """
        from params import DEFAULT_PARAMS
        self.params = params if params is not None else DEFAULT_PARAMS
        self.arm = int(arm) if arm is not None else int(self.params.gibson_overlap)

    def linearize(self, vector: Dseqrecord) -> Dseqrecord:
        """유저가 고른(또는 자동 선택된) 효소로 벡터를 연다.

        `GibsonCloning._linearize_vector`를 그대로 쓴다 — 삽입 단계와 **같은 자리**에서
        잘라야 homology가 맞는다. 여기서 따로 자르면 미묘하게 어긋난다.
        """
        from assembly.factory import AssemblyFactory
        g = AssemblyFactory.create_assembler("gibson", params=self.params)
        return g._linearize_vector(vector)

    def make(self, vector: Dseqrecord) -> Tuple[Dseqrecord, Dseqrecord]:
        """(L spacer, R spacer) — 각각 벡터의 3' 끝 / 5' 끝 서열.

        Raises:
            RuntimeError: 벡터가 homology arm 두 개를 뗄 만큼 길지 않을 때.
        """
        ov = self.arm
        lin = self.linearize(vector)
        seq = str(lin.seq).upper()
        if len(seq) < 2 * ov:
            raise RuntimeError(
                f"vector가 너무 짧습니다({len(seq)}bp) — homology arm {ov}bp 두 개를 뗄 수 없습니다."
            )

        l_seq = seq[-ov:]   # 선형 벡터의 3' 끝 → module 왼쪽 끝
        r_seq = seq[:ov]    # 선형 벡터의 5' 끝 → module 오른쪽 끝

        # GenBank LOCUS 이름은 16자 제한이라 잘린다. 접미사(_homL/_homR)를 먼저 확보하고
        # **앞쪽(벡터명)을 줄인다** — 뒤를 자르면 L과 R이 같은 이름이 돼 리포트에서 구분이 안 된다.
        base = (vector.name or "vector")[:11]
        l_rec = Dseqrecord(l_seq, name=f"{base}_homL")
        r_rec = Dseqrecord(r_seq, name=f"{base}_homR")
        for rec in (l_rec, r_rec):
            rec.annotations["part_type"] = "spacer"
            rec.annotations["molecule_type"] = "DNA"
        return l_rec, r_rec

    def attach(self, parts, vector: Dseqrecord):
        """[L_hom, *parts, R_hom] — SpacerSelector.attach와 같은 계약."""
        l_rec, r_rec = self.make(vector)
        return [l_rec, *parts, r_rec]


def module_carries_homology(insert: Dseqrecord, linear_vector: Dseqrecord,
                            overlap: int) -> bool:
    """module 양 끝이 이미 벡터 homology인가 — 맞으면 삽입 PCR을 건너뛸 수 있다.

    플래그가 아니라 **서열을 직접 대조**한다. 어떤 경로로 만들어졌든 실제로 맞을 때만
    건너뛰므로, 모드 설정이 꼬여도 잘못된 산물이 나오지 않는다.
    """
    ins = str(insert.seq).upper()
    vec = str(linear_vector.seq).upper()
    if len(ins) < 2 * overlap or len(vec) < 2 * overlap:
        return False
    return ins.startswith(vec[-overlap:]) and ins.endswith(vec[:overlap])


__all__ = ["VectorHomologySpacer", "module_carries_homology"]
