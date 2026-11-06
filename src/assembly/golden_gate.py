"""
Golden Gate Assembly 구현

Type IIS 제한효소 (BsaI)를 사용한 DNA 조립 방법입니다.

핵심 기능:
    - BsaI/BbsI cutting (1bp spacer + 4bp overhang)
    - 부품 내부 BsaI 사이트 존재 시 조립 전체 BbsI 자동 전환
    - Overhang 호환성 체크
    - Circular/Linear assembly

워크플로우:
    1. 입력 검증 (DNA sequences)
    2. BsaI/BbsI cutting 수행
    3. Overhang 호환성 확인
    4. Assembly 시뮬레이션 (PyDNA)

참고:
    - KRIBB Line 230-396
    - archive/docs/development/PLAN_v2.md #2 GoldenGateModule

작성일: 2025-11-06
"""

from typing import List, Dict, Tuple, Optional
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from Bio.Restriction import BsaI, BbsI
from Bio.Seq import Seq
from .base import AssemblyBase


def _rc(seq: str) -> str:
    """역상보. overhang 방향을 다룰 때 반복적으로 쓰인다."""
    return str(Seq(seq).reverse_complement()) if seq else ""

# Basic 및 Utils 모듈 임포트
from basic.primer_design import PrimerDesign
from utils.io import DatabaseImport
from basic.sequence_analysis import check_restriction_site_in_part_sequence, get_typeIIS_params, validate_dna_sequence


#: Golden Gate로 벡터를 열려면 인식 서열이 **몇 곳 이상** 있어야 하는가.
#  1곳만 자르면 좌우 overhang이 서로 역상보가 되어 insert 없이 도로 붙고,
#  방향성 클로닝이 성립하지 않는다. (클라이언트 `linearizationEnzymes.ts`와 같은 규칙)
VECTOR_MIN_SITES = 2


def select_gg_enzyme(parts: List[Dseqrecord], vector: Optional[Dseqrecord] = None,
                     preferred: Optional[str] = None) -> Dict:
    """조립 전체가 쓸 Type IIS 효소를 **한 번에** 결정한다.

    ## 왜 한 곳에서 정해야 하는가

    Golden Gate는 one-pot 반응이라 벡터와 모든 부품이 **같은 효소**로 잘려야 한다.
    그런데 두 쪽의 요구 조건이 정반대다:

    - **벡터**: 그 효소 사이트가 `VECTOR_MIN_SITES`곳 이상 (열어서 좌우 다른 overhang)
    - **부품**: 그 효소 사이트가 **0개** (있으면 의도치 않은 자리에서 잘린다)

    예전에는 이 결정이 세 군데에 흩어져 있었다 — `params.gg_enzyme`(유저 선택),
    `GoldenGateModule.prepare_fragments`의 자동 전환(**부품만** 보고 고름),
    `WorkflowBase.assemble`의 벡터 삽입(**BsaI 하드코딩**). 서로 모르는 채 각자 정해서,
    DB 실측 기준 **될 수 있었던 설계의 32.7%p를 놓치고 있었다**
    (GG 가능 벡터 202개 × 무작위 4-part module 3000회: 공통 효소 존재 68.7% vs 성공 36.0%).

    ## 고르는 순서

    1. 두 조건의 **교집합**을 구한다. 순차로 "벡터 먼저 → 부품"을 하면 벡터 기준으로
       고른 효소가 어떤 부품 안에 들어 있어 되돌아와야 하는 경우가 생긴다.
    2. 교집합에 `preferred`가 있으면 그걸 쓴다 — 유저가 고른 값을 말없이 바꾸지 않는다.
    3. 없으면 `SUPPORTED_ENZYMES` 순서로 첫 번째를 쓰고 `switched=True`로 알린다.
    4. 교집합이 비면 `enzyme=None`. 호출자가 안내 문구(`reason`)를 띄운다.

    교집합이 빌 때 안내는 **부품 쪽을 가리킨다** — 벡터 제약은 고칠 수 없지만
    (사이트를 새로 넣으려면 다시 클로닝해야 한다) 부품은 도메스티케이션으로 내부
    사이트를 지울 수 있기 때문이다.

    Args:
        parts: 조립에 들어가는 부품 Dseqrecord 리스트
        vector: 삽입 대상 벡터. None이면 **부품 조건만** 본다 — module만 만들고
            벡터에 넣지 않는 경로(part→module)에서는 벡터 제약이 없다.
        preferred: 유저가 고른 효소(`params.gg_enzyme`). 교집합에 있으면 그대로 쓴다.

    Returns:
        ``{'enzyme', 'candidates', 'switched', 'reason', 'vector_sites', 'blocked_by'}``

        - ``enzyme``: 결정된 효소 이름. 불가능하면 ``None``
        - ``candidates``: 두 조건을 모두 만족하는 효소 목록
        - ``switched``: ``preferred``와 다른 효소를 골랐는가
        - ``reason``: 사람이 읽을 설명(불가능하거나 전환했을 때만 채운다)
        - ``vector_sites``: 효소별 벡터 절단 횟수 (vector가 None이면 빈 dict)
        - ``blocked_by``: 효소별로 내부에 그 사이트를 가진 부품 이름들
    """
    from params import SUPPORTED_ENZYMES
    from basic.sequence_analysis import find_restriction_sites

    enzymes = list(SUPPORTED_ENZYMES)

    # 1) 부품 조건 — 내부에 인식 서열이 있으면 그 효소는 못 쓴다
    blocked_by: Dict[str, List[str]] = {e: [] for e in enzymes}
    for p in parts:
        res = check_restriction_site_in_part_sequence(str(p.seq), enzymes=enzymes)
        for e in enzymes:
            if res['found_sites'].get(e):
                blocked_by[e].append(p.name)
    part_ok = [e for e in enzymes if not blocked_by[e]]

    # 2) 벡터 조건 — 2곳 이상 잘라야 열린다
    vector_sites: Dict[str, int] = {}
    if vector is not None:
        found = find_restriction_sites(str(vector.seq).upper())
        vector_sites = {e: len(found.get(e, [])) for e in enzymes}
        vec_ok = [e for e in enzymes if vector_sites[e] >= VECTOR_MIN_SITES]
    else:
        vec_ok = enzymes

    candidates = [e for e in enzymes if e in part_ok and e in vec_ok]

    if not candidates:
        bits = []
        if vector is not None:
            usable = [f"{e}({vector_sites[e]}곳)" for e in enzymes if vector_sites[e] >= VECTOR_MIN_SITES]
            vname = getattr(vector, "name", "vector")
            bits.append(
                f"벡터 '{vname}'를 {VECTOR_MIN_SITES}곳 이상 자르는 효소: "
                + (", ".join(usable) if usable else "없음")
            )
        for e in enzymes:
            if blocked_by[e]:
                bits.append(f"{e}는 부품 내부에 사이트가 있음: {', '.join(blocked_by[e])}")
        return {
            'enzyme': None, 'candidates': [], 'switched': False,
            'vector_sites': vector_sites, 'blocked_by': blocked_by,
            'reason': (
                "Golden Gate 조립이 불가능합니다 — 벡터와 모든 부품을 동시에 만족하는 "
                f"Type IIS 효소가 없습니다.\n   - " + "\n   - ".join(bits)
                + "\n   → 부품 내부 사이트를 제거(도메스티케이션)하거나, 다른 벡터를 "
                  "고르거나, Gibson·Module based 방식을 쓰세요."
            ),
        }

    if preferred in candidates:
        return {'enzyme': preferred, 'candidates': candidates, 'switched': False,
                'reason': "", 'vector_sites': vector_sites, 'blocked_by': blocked_by}

    chosen = candidates[0]
    why = []
    if preferred and blocked_by.get(preferred):
        why.append(f"부품 내부에 {preferred} 사이트가 있음({', '.join(blocked_by[preferred])})")
    if preferred and vector is not None and vector_sites.get(preferred, 0) < VECTOR_MIN_SITES:
        why.append(f"벡터가 {preferred}로 {vector_sites.get(preferred, 0)}곳만 잘림")
    return {
        'enzyme': chosen, 'candidates': candidates, 'switched': bool(preferred),
        'vector_sites': vector_sites, 'blocked_by': blocked_by,
        'reason': (
            f"조립 효소를 {preferred} → {chosen}로 전환했습니다"
            + (f" ({'; '.join(why)})" if why else "")
            + f". 사용 가능 후보: {candidates}"
        ) if preferred else f"조립 효소로 {chosen}를 선택했습니다. 후보: {candidates}",
    }


class GoldenGateModule(AssemblyBase):
    """Golden Gate Assembly (Type IIS 제한효소) 구현

    특징:
        - BsaI/BbsI 제한효소 지원
        - 내부 BsaI 사이트 감지 시 조립 전체 BbsI 자동 전환
        - 1bp spacer + 4bp overhang 자동 추출
        - Overhang 기반 순서 조립

    Attributes:
        enzyme: 기본 제한효소 (기본값: BsaI)

    작성일: 2025-11-06
    """

    def __init__(self, enzyme=BsaI, params=None):
        """Golden Gate Assembly 초기화

        Args:
            enzyme: Type IIS 제한효소. 이름(str, 예: "BsmBI") 또는 Biopython 효소 객체.
                기본값 BsaI. 지원 목록은 `params.SUPPORTED_ENZYMES`.

        Note:
            **overhang 길이와 spacer 길이는 파라미터가 아니라 효소에서 파생된다.**
            (BsaI/BbsI/BsmBI는 4nt, SapI는 3nt) — 코드 곳곳에 4를 박아두면 SapI에서 깨진다.
        """
        from Bio import Restriction
        from basic.sequence_analysis import get_typeIIS_params

        # 이름으로 받으면 Biopython 객체로 정규화 (galaxy/params는 문자열로 넘긴다)
        from params import SUPPORTED_ENZYMES
        name = enzyme if isinstance(enzyme, str) else str(enzyme)
        # Biopython에 존재하는지가 아니라 **지원 목록에 있는지**로 검증한다.
        # (EcoRI처럼 Type IIS가 아닌 효소도 Biopython엔 있으므로 hasattr 검사만으론 못 거른다)
        if name not in SUPPORTED_ENZYMES:
            raise ValueError(
                f"지원하지 않는 제한효소: {name!r} (지원: {list(SUPPORTED_ENZYMES)}). "
                f"Type IIS 효소만 Golden Gate에 쓸 수 있습니다."
            )
        enzyme = getattr(Restriction, name)

        self.enzyme = enzyme
        self.enzyme_name = name
        # 설계 파라미터 한 벌(연구실 조건). 미지정 시 기본값 = 기존 하드코딩 동작.
        # 효소는 이 객체가 이미 확정했으므로 params의 gg_enzyme과 다를 수 있다(내부 사이트 회피 등).
        from params import DEFAULT_PARAMS
        self.params = params if params is not None else DEFAULT_PARAMS
        _ep = get_typeIIS_params(self.enzyme_name)
        self.overhang_len = int(_ep["overhang_len"])   # 효소 파생 (BsaI 4 / SapI 3)
        self.spacer_len = int(_ep["spacer_len"])       # 효소 파생 (BsaI 1 / BbsI 2)
        # one-pot 교차 dimer 결과(2026-11-05 Y12). `prepare_fragments`가 채우고
        #   리포트(`reports/cross_checks.csv`)가 읽는다. 쌍 단위 결과라 부품 행에
        #   담을 자리가 없어서 **설계기가 들고 있다가 넘기는** 방식을 쓴다
        #   (`last_primer_design`과 같은 패턴 — Y10-b).
        self.last_cross_dimers = []
        #: 교차 dimer로 **비교한 쌍 수**. 걸린 쌍이 0일 때 "검사했고 깨끗하다"를
        #  말하려면 분모가 필요하다(`_DimerScan.pairs_checked`).
        self.last_cross_dimer_pairs = 0

    # TODO: 현재 기능상으로는 항상 랜덤으로 overhang을 설계해서 나옴.
    # 관계자분 피드백으로는 물리적으로 파츠가 존재할때 overhang 재활용 관점이 굉장히 중요하다고하심.
    #: 접합부 문맥을 볼 때 양옆에서 가져오는 염기 수.
    #  인식서열이 최대 8bp(NotI)라, overhang 4bp와 걸쳐도 8nt면 충분히 덮는다.
    _JUNCTION_FLANK = 8

    def _junction_creates_site(self, left_seq: str, oh: str, right_seq: str,
                               enzymes) -> bool:
        """이 overhang을 끼우면 접합부에 **없던 제한효소 사이트가 생기는가**.

        overhang 자체에는 사이트가 없어도, 양옆 부품의 말단과 만나면서 인식서열이
        완성될 수 있다. 실제로 PhlF|L3S2P11 접합부에서 overhang + 터미네이터 첫 염기가
        `GGTCTC`(BsaI)를 만들어, **원팟 반응 내내 그 자리가 다시 잘리는** 산물이 나왔다
        (2026-10-12 실측). one-pot Golden Gate는 비가역성이 핵심이라 이게 깨지면
        수율이 떨어지고 조각이 엉뚱하게 재결합할 수 있다.

        규약 혼동을 피하려고 overhang과 그 역상보를 **둘 다** 끼워 본다.
        """
        from basic.sequence_analysis import find_restriction_sites
        l = (left_seq or "")[-self._JUNCTION_FLANK:].upper()
        r = (right_seq or "")[:self._JUNCTION_FLANK].upper()
        base = {e: len(find_restriction_sites(l + r, e)[e]) for e in enzymes}
        for cand in {oh.upper(), _rc(oh.upper())}:
            joined = l + cand + r
            if any(len(find_restriction_sites(joined, e)[e]) > base[e] for e in enzymes):
                return True
        return False

    def _avoid_site_creating_overhangs(self, internal_ohs: List[str],
                                       parts, enzymes) -> List[str]:
        """접합부에 사이트를 만드는 overhang을 쓰지 않는 후보로 바꾼다.

        교체 비용은 거의 없다 — 실측상 최악의 접합부에서도 Awsemble overhang 165개 중
        159개가 깨끗했고, 나머지 접합부는 전부 165/165였다. 그래도 전부 걸리면
        **조용히 진행하지 않고** 에러로 알린다(원팟에서는 "덜 나쁜 것"도 결국 잘린다).
        """
        if not internal_ohs:
            return internal_ohs
        pool = [o.upper() for o in getattr(self, "_oh_id_map", {}).keys()]
        used = {o.upper() for o in internal_ohs}
        out = list(internal_ohs)

        for j, oh in enumerate(internal_ohs):
            left = str(parts[j].seq)
            right = str(parts[j + 1].seq)
            if not self._junction_creates_site(left, oh, right, enzymes):
                continue
            alt = next(
                (c for c in pool
                 if c not in used and _rc(c) not in used
                 and not self._junction_creates_site(left, c, right, enzymes)),
                None,
            )
            if alt is None:
                raise RuntimeError(
                    f"접합부 '{parts[j].name}|{parts[j+1].name}'에서 어떤 overhang 후보도 "
                    f"제한효소 사이트 생성을 피하지 못했습니다(후보 {len(pool)}개). "
                    f"원팟 Golden Gate에서는 그 자리가 반응 내내 다시 잘리므로 진행하지 않습니다. "
                    f"부품 말단을 조금 옮기거나 다른 효소를 쓰세요."
                )
            print(f"     🔄 접합부 {parts[j].name}|{parts[j+1].name}: overhang {oh} → {alt} "
                  f"(원래 것은 제한효소 사이트를 만듦)")
            used.discard(oh.upper())
            used.add(alt)
            out[j] = alt
        return out

    def _get_high_fidelity_overhang_set(self, count: int) -> List[str]:
        """DB에서 Awsemble 오버행 후보군을 가져와 최적의 서열 세트를 반환합니다.

        **효소의 overhang 길이와 같은 후보만 쓴다.** 길이가 다르면 조립 단계에서
        "Overhang 불일치"라는 원인 불명 에러로 터지므로(예: SapI 3nt인데 후보는 4nt),
        여기서 미리 걸러 실제 원인을 알린다.
        """
        db = DatabaseImport()

        def get_rc(seq: str):
            return "".join([{'A':'T','T':'A','G':'C','C':'G'}[b] for b in seq[::-1]])

        def get_hd(s1, s2):
            return sum(c1 != c2 for c1, c2 in zip(s1, s2))

        pool_list = []
        if db._connect("Bioparts"):
            try:
                cursor = db.conn.cursor()
                cursor.execute("SELECT ID, Sequence FROM `overhang` WHERE `institution` = 'Awsemble' ORDER BY Sequence")
                rows = cursor.fetchall()
                pool_list = [r[1] for r in rows]
                self._oh_id_map = {str(r[1]).upper(): r[0] for r in rows}  # 서열→id (report용)
                cursor.close()
            finally:
                db._disconnect()

        if not pool_list:
            pool_list = sorted(["ATGC", "GCTA", "CGAT", "TAGC", "CCAG", "GGTC", "AATG", "TTAC"])

        # 효소가 요구하는 길이의 후보만 남긴다.
        oh_len = self.overhang_len
        sized = [q for q in pool_list if len(q) == oh_len]
        if not sized:
            have = sorted({len(q) for q in pool_list})
            raise ValueError(
                f"{self.enzyme_name}는 {oh_len}nt overhang이 필요한데, 사용 가능한 overhang 후보가 "
                f"{have}nt뿐입니다(후보 {len(pool_list)}개). "
                f"{oh_len}nt overhang 세트를 DB(`overhang` 테이블, institution='Awsemble')에 "
                f"등록해야 이 효소를 쓸 수 있습니다."
            )
        pool_list = sized

        selected_seqs = []
        for stage in range(1, 4):
            for seq in pool_list:
                if len(selected_seqs) >= count: break
                seq_u = seq.upper()
                if seq_u in selected_seqs: continue
                rc_u = get_rc(seq_u)
                if all((get_hd(seq_u, s) >= (2 if stage < 3 else 1) and 
                        (stage > 1 or get_hd(rc_u, s) >= 2)) for s in selected_seqs):
                    selected_seqs.append(seq_u)
        return selected_seqs

    def prepare_fragments(
        self,
        parts: List[Dseqrecord],
        vector: Optional[Dseqrecord] = None,
        assembly_type: str = "linear"
    ) -> List[Dict]:
        """조립 준비 (오버행 할당 및 프라이머 설계)"""
        if not parts:
            raise ValueError("parts 리스트가 비어있습니다")

        num_parts = len(parts)

        # 1. [핵심 로직] 부품 내부 제한효소 사이트 스캔 → 사용할 효소 결정
        #    조립에 쓸 효소의 인식 서열이 부품 **내부**에 있으면 조립 중 의도치 않게 잘린다.
        #    구버전은 BsaI/BbsI 두 개만 보고 BsaI가 있으면 말없이 BbsI로 바꿨다.
        #    효소를 유저가 고르는 구조가 되면서 지원 효소 전체를 스캔하도록 일반화한다.
        from params import SUPPORTED_ENZYMES

        found: Dict[str, List[str]] = {name: [] for name in SUPPORTED_ENZYMES}
        for p in parts:
            # E-1: 입력 부품 검증 — ATGC만 허용(N/IUPAC·빈 서열 거부). 모호코드는 조립 불가로 차단.
            if not validate_dna_sequence(str(p.seq)):
                raise ValueError(
                    f"입력 부품 '{p.name}'에 ATGC 외 염기(N/IUPAC) 또는 빈 서열 포함 → 조립 불가 (E-1 입력검증)"
                )
            res_check = check_restriction_site_in_part_sequence(
                str(p.seq), enzymes=list(SUPPORTED_ENZYMES)
            )
            for name in SUPPORTED_ENZYMES:
                if res_check['found_sites'].get(name):
                    found[name].append(p.name)

        chosen = self.enzyme_name
        if not found.get(chosen):
            target_enzyme = chosen
        else:
            # 고른 효소가 부품 내부에 있음 → 그 효소로는 조립 불가. 내부 사이트가 없는 효소로 대체.
            # ⚠️ 이 자동 대체는 결과를 유저 모르게 바꾸므로 **크게 알린다**.
            #    (요청 문서 C-1에서 "고른 효소만 쓰고 충돌은 에러로" 방향의 동의를 받는 중)
            alts = [n for n in SUPPORTED_ENZYMES if n != chosen and not found.get(n)]
            if not alts:
                detail = "\n".join(
                    f"   - {n} 발견: {', '.join(v)}" for n, v in found.items() if v
                )
                error_msg = (
                    f"❌ 제한효소 충돌: 지원 효소 {list(SUPPORTED_ENZYMES)} 전부가 부품 내부에 "
                    f"인식 서열을 갖고 있어 Golden Gate 조립이 불가능합니다.\n{detail}\n"
                    f"   - 해결책: 해당 부품 서열에서 인식 서열을 제거(도메스티케이션)하거나 다른 조립 방식을 쓰세요."
                )
                print(error_msg)
                raise RuntimeError(error_msg)
            target_enzyme = alts[0]
            print(
                f"     ⚠️ 선택한 {chosen} 인식 서열이 부품 내부에 있습니다"
                f"({', '.join(found[chosen])}) → 조립 전체를 {target_enzyme}로 전환합니다."
            )

        # TODO: 위와 동일하게 overhang 재사용 고려
        # 또한 현재 part끼리만 결합하여 module을 제작할 때는 linear만 고려되는데 circular가 어케 나온거지?

        # 2. 오버행 할당 명단(all_ohs) 생성
        if assembly_type == "circular":
            if not vector:
                raise ValueError("Circular 조립을 위해서는 vector가 필요합니다.")
            v_left, v_right = self.extract_vector_overhangs(vector)
            internal_ohs = self._get_high_fidelity_overhang_set(num_parts - 1)
            internal_ohs = self._avoid_site_creating_overhangs(
                internal_ohs, parts, list(SUPPORTED_ENZYMES))
            # ⚠️ 벡터 쪽 fusion site는 **RC**해서 넣는다.
            #   extract_vector_overhangs는 backbone의 top-strand 말단 서열을 그대로 준다.
            #   부품 fragment는 규약상 left_oh=RC(fusion)을 노출하므로, 그대로 쓰면 부품의
            #   말단 overhang이 backbone의 노출 overhang과 **상보가 아니라 동일**해져 결합이 안 된다
            #   (pydna: "sticky ends not compatible"). RC를 한 번 더 적용해야 맞물린다.
            all_ohs = [_rc(v_right)] + internal_ohs + [_rc(v_left)]
        else:
            internal_ohs = self._get_high_fidelity_overhang_set(num_parts - 1) if num_parts > 1 else []
            internal_ohs = self._avoid_site_creating_overhangs(
                internal_ohs, parts, list(SUPPORTED_ENZYMES))
            all_ohs = [""] + internal_ohs + [""]

        # 3. 통합 설계 수행
        from Bio.Seq import Seq as _Seq
        oh_id_map = getattr(self, '_oh_id_map', {})   # 서열→overhang id (report용)
        design_results = []
        print(f"\n✨ Golden Gate 부품 설계 상세 (모드: {assembly_type}, Enzyme: {target_enzyme})")
        for i, part in enumerate(parts):
            # [FIX] 인접 부품의 노출 overhang이 상보(RC)가 되도록 left_oh만 RC로 할당.
            #   part i는 bottom 5'overhang = right_oh(공유 fusion 그대로),
            #   part i+1은 top 5'overhang = left_oh → 실제 ligation하려면 둘이 RC여야 한다.
            #   따라서 left_oh = RC(공유 fusion all_ohs[i]). (report id는 RC 전 fusion 기준 = fusion_left)
            fusion_left = all_ohs[i]
            left_oh = str(_Seq(fusion_left).reverse_complement()).upper() if fusion_left else ""
            right_oh = all_ohs[i+1]

            design = PrimerDesign.design_goldengate_primers(
                insert=part,
                left_overhang=left_oh,
                right_overhang=right_oh,
                restriction_enzyme=target_enzyme,
                params=self.params,
            )

            unique_name = f"{part.name}_{i}_designed"
            
            # [수정] BbsI(12bp 이상)를 고려하여 로그 출력 길이를 효소에 맞춰 조절
            # 부착부 = [인식서열][spacer][overhang] → 길이를 효소에서 파생 (구: BbsI 12 / 그 외 11 고정)
            _tp = get_typeIIS_params(target_enzyme)
            show_len = len(_tp["fwd"]) + _tp["spacer_len"] + _tp["overhang_len"]
            actual_left = design['start_structure'][:show_len] if left_oh else "(None)"
            actual_right = design['end_structure'][-show_len:] if right_oh else "(None)"

            print(f"  📦 [{i}] {unique_name}:")
            print(f"     L-End: {actual_left} | R-End: {actual_right}")

            design_results.append({
                'part_name': unique_name,
                # 원본 part id/type 보존 — module feature(P6)에서 label=part_id, GenBank type=part_type.
                # (part_name은 prep 추적용 유니크명 "<id>_<i>_designed"라 feature label로 부적합)
                'part_id': part.name,
                'part_type': part.annotations.get('part_type', ''),
                'forward_primer': design['forward'],
                'reverse_primer': design['reverse'],
                'start_structure': design['start_structure'],
                'end_structure': design['end_structure'],
                'full_sequence': design['full_sequence'],
                'target_seq': design['target_seq'],
                'restriction_enzyme': target_enzyme, # 결정된 효소 정보 저장
                # overhang 추적(report용): 인접 부품은 같은 fusion site(id) 공유 → 정합 눈으로 확인.
                #   left_oh 서열은 RC(부착 실제값), id는 RC 전 fusion(fusion_left) 기준이라 인접 right와 id 일치.
                'left_oh': left_oh,
                'left_oh_id': oh_id_map.get(fusion_left.upper(), "") if fusion_left else "",
                'right_oh': right_oh,
                'right_oh_id': oh_id_map.get(right_oh.upper(), "") if right_oh else "",
                'fwd_tm': design['fwd_tm'],
                'rev_tm': design['rev_tm'],
                # 2026-11-05(Y12) — 주형 결합부(tail 제외). 교차 dimer 검사가 쓴다.
                'fwd_anneal': design.get('fwd_anneal'),
                'rev_anneal': design.get('rev_anneal'),
                # 조립용 프라이머는 [5' tail][anneal] 구조라 두 Tm이 크게 다르다.
                # 위 *_tm은 anneal 기준, 아래 *_tm_full은 주문할 프라이머 전체 기준.
                'fwd_tm_full': design.get('fwd_tm_full'),
                'rev_tm_full': design.get('rev_tm_full'),
                'tm_na_mM': design.get('tm_na_mM'),
                'tm_primer_nM': design.get('tm_primer_nM'),
                # 리포트(primers.csv)용 품질 지표 — GC%/self-hairpin/template accessibility/tier/경고
                'fwd_gc': design.get('fwd_gc'),
                'rev_gc': design.get('rev_gc'),
                'fwd_hairpin_dg': design.get('fwd_hairpin_dg'),
                'rev_hairpin_dg': design.get('rev_hairpin_dg'),
                'fwd_accessibility': design.get('fwd_accessibility'),
                'rev_accessibility': design.get('rev_accessibility'),
                'tier': design.get('tier'),
                'warnings': design.get('warnings', []),
                'passed': design['passed'],
                'violations': design['violations'],
                'structural_notes': design.get('structural_notes', [])
            })

        # ── one-pot 교차 dimer (2026-11-05 Y12) ────────────────────────────
        # 여기가 올바른 자리인 이유: `prepare_fragments` 한 번이 곧 **한 튜브**다.
        #   부품 n개의 프라이머 2n개가 같이 들어가므로 전수 쌍의 범위가 정확히 이것이다.
        #   (부품 하나씩 설계하는 `design_goldengate_primers`는 자기 쌍밖에 모른다)
        self.last_cross_dimers = self._check_cross_dimers(design_results)
        return design_results

    def _check_cross_dimers(self, design_results: List[Dict]) -> List[Dict]:
        """같은 튜브에 들어가는 프라이머들의 **결합부끼리** dimer를 전수 검사한다.

        걸린 쌍은 ① 연루된 **양쪽 행**의 violations에 넣고(시트만 보고 주문하는 사람이
        놓치면 안 된다) ② 전체 목록을 돌려준다(`reports/cross_checks.csv`로 나간다).

        왜 결합부만 보는지, 왜 경고가 아니라 차단인지는
        `basic.primer_design.check_cross_dimers` 독스트링에 있다.

        Args:
            design_results: `prepare_fragments`가 만든 설계 결과 리스트(제자리 수정됨)

        Returns:
            `check_cross_dimers`의 결과 리스트
        """
        from basic.primer_design import check_cross_dimers

        anneals = []
        for d in design_results:
            for side in ("fwd", "rev"):
                seq = d.get(f"{side}_anneal")
                if seq:
                    anneals.append((f"{d['part_name']} {side}", seq))
        hits = check_cross_dimers(anneals, self.params)
        self.last_cross_dimer_pairs = getattr(hits, "pairs_checked", 0)
        if not hits:
            return []

        # 라벨("<part_name> fwd") → 그 부품의 행
        #   ⚠️ **차단 등급(blocking)만 행에 적는다.** 중간 등급까지 적으면 기본값에서
        #     과반이 red가 된다(실측 39행 중 24행 = 62%). 중간 등급은
        #     `reports/cross_checks.csv`에 남아 필요한 사람이 볼 수 있다.
        by_name = {d["part_name"]: d for d in design_results}
        blocking = [h for h in hits if h.get("blocking")]
        for h in blocking:
            for lab, other in ((h["a"], h["b"]), (h["b"], h["a"])):
                row = by_name.get(lab.rsplit(" ", 1)[0])
                if row is None:
                    continue
                row.setdefault("violations", []).append(
                    f"Cross dimer: the {lab.rsplit(' ', 1)[1]} annealing region binds {other} "
                    f"strongly (dG {h['delta_g']:.2f} kcal/mol, blocking cutoff "
                    f"{self.params.qc_dg_absolute:g}) - they share a tube, so these "
                    f"primers bind each other instead of the template. See "
                    f"reports/cross_checks.csv for the full pair list."
                )
                row["passed"] = False
                row["tier"] = "red"
        if blocking:
            print(f"\n  [BLOCKED] one-pot cross dimer: {len(blocking)} pairs blocked "
                  f"({len(hits)} recorded) - reports/cross_checks.csv")
        elif hits:
            print(f"\n  [NOTE] one-pot cross dimer: {len(hits)} pairs recorded "
                  f"(none at blocking grade) - reports/cross_checks.csv")
        for h in hits[:4]:
            print(f"     {'blocked' if h.get('blocking') else 'noted'} "
                  f"{h['a']} ↔ {h['b']}: ΔG {h['delta_g']:.2f}")
        return hits

    # TODO: 이 코드가 내가 알기로는 제한효소를 자른뒤에 실제 삽입시켜야하는 서열을 찾는거로 알고 있는데
    # 당연히 overhang 사이에 있는걸 고르도록 해야하지않나? 이게 왜 필요하지?

    def extract_insert_fragment(self, dseq: Dseqrecord) -> Optional[Dseqrecord]:
        """prepared fragment 추출 — prep이 만든 full_sequence를 알려진 오프셋으로 슬라이스 (H-4).

        full_sequence = [site][spacer][left_oh] + target + [right_oh_rc][spacer_rc][site_rev]
        → 절단 후 중간 조각 = [left_oh][target][right_oh_rc]. prep이 구조를 이미 알므로
        효소 재절단(pydna cut) 대신 결정적 슬라이스로 산출한다.
        효과: 짧은 부품의 pydna cut 실패(no_output) 회피 + 비결정 유사도매칭 제거(재현성↑).

        [site][spacer] 길이는 get_typeIIS_params로 효소별 동적 도출(BsaI=7, BbsI=8).
        linear 양 끝(blunt)은 site가 없으므로 startswith/endswith로 자동 판별→오프셋 0 (C-2).
        (구 재절단+유사도 매칭은 _extract_insert_fragment_by_recut에 백업 보존 — junction-site로
         조각 수 가변 문제가 재발하면 fallback으로 사용.)
        """
        seq = str(dseq.seq)
        seq_u = seq.upper()
        enzyme_name = dseq.annotations.get('restriction_enzyme', "BsaI")
        params = get_typeIIS_params(enzyme_name)
        site_fwd = params['fwd']
        site_rev = params['rev']
        trim = len(site_fwd) + params['spacer_len']  # [site][spacer] 길이

        # 인식서열 앞에 5' 여분(pad)이 붙을 수 있으므로 startswith/endswith로는 못 찾는다.
        # (pad가 붙었는데 startswith가 실패하면 트림이 0이 되어 **자르지 않은 통짜**가 나가고,
        #  조립 단계에서 overhang 불일치로 뒤늦게 터진다 — 조용히 틀리는 경로였다)
        # → 양 끝 근처에서 인식서열 위치를 찾아 pad 길이만큼 함께 트림한다.
        #   부품은 조립 효소 사이트가 없도록 이미 걸러졌으므로(prepare_fragments),
        #   말단 근처에서 찾은 사이트는 설계가 붙인 것이 확실하다.
        MAX_PAD = 20   # 여분이 이보다 길 일은 없다(길이 예산상 불가) — 오탐 방지 상한

        i = seq_u.find(site_fwd)
        front = (i + trim) if 0 <= i <= MAX_PAD else 0

        j = seq_u.rfind(site_rev)
        tail_pad = len(seq) - (j + len(site_rev)) if j >= 0 else -1
        back = (len(seq) - j + params['spacer_len']) if 0 <= tail_pad <= MAX_PAD else 0

        frag_seq = seq[front: len(seq) - back]

        frag = Dseqrecord(frag_seq)
        frag.name = dseq.name
        print(f"  🔪 [{dseq.name}] prepared fragment 슬라이스 (Enzyme: {enzyme_name}, "
              f"front={front} back={back} → {len(frag_seq)}bp)")
        return frag

    # 백업(ljw 지시): junction-site로 조각 수가 가변되는 문제가 추후 재발하면
    # 아래 재절단+유사도 매칭 방식으로 폴백. 현재 핵심 경로는 슬라이스(extract_insert_fragment)가 담당, 본 메서드는 미사용.
    def _extract_insert_fragment_by_recut(self, dseq: Dseqrecord) -> Optional[Dseqrecord]:
        """[백업] 제한효소 재절단 후 원본과 가장 유사한 단편 선택 (구 extract_insert_fragment)."""
        from difflib import SequenceMatcher

        orig_seq = dseq.annotations.get('target_seq', "").upper()
        if not orig_seq:
            orig_seq = str(dseq.seq).upper()

        target_enzyme_name = dseq.annotations.get('restriction_enzyme', "BsaI")
        target_enzyme = BbsI if target_enzyme_name == "BbsI" else BsaI

        part_frags = self._cut_with_enzyme(dseq, target_enzyme)
        if not part_frags:
            return None

        best_frag = None
        max_relative_score = 0
        params = get_typeIIS_params(target_enzyme_name)
        site_fwd = params["fwd"]
        site_rev = params["rev"]

        for frag in part_frags:
            frag_seq_str = str(frag.seq).upper()
            analysis_seq = frag_seq_str
            t_len = len(site_fwd) + params["spacer_len"]
            if frag_seq_str.startswith(site_fwd):
                analysis_seq = frag_seq_str[:t_len]
            elif frag_seq_str.endswith(site_rev):
                analysis_seq = frag_seq_str[-t_len:]

            display_len = len(analysis_seq)
            actual_ratio = SequenceMatcher(None, analysis_seq, orig_seq, autojunk=False).ratio()
            perfect_ratio = (2.0 * len(orig_seq)) / (display_len + len(orig_seq))
            relative_score = actual_ratio / perfect_ratio if perfect_ratio > 0 else 0
            if relative_score > max_relative_score:
                max_relative_score = relative_score
                best_frag = frag

        if best_frag and max_relative_score >= 0.8:
            return best_frag
        return max(part_frags, key=lambda x: len(x.seq))
        
    def amplify_insert(
        self,
        insert: Dseqrecord,
        forward_primer: str,
        reverse_primer: str
    ) -> Dseqrecord:
        """Insert DNA를 PCR로 증폭 (프라이머 부착)"""
        from pydna.amplify import pcr
        try:
            pcr_product = pcr(forward_primer, reverse_primer, insert)
            pcr_product.name = insert.name
            # [중요] 원본의 annotations(특히 target_seq, restriction_enzyme)을 PCR 결과물로 복사
            pcr_product.annotations.update(insert.annotations)
            return pcr_product
        except Exception as e:
            raise RuntimeError(f"PCR 시뮬레이션 실패: {e}")

    def assemble(
            self,
            fragments: List[Dseqrecord],
            assembly_type: str = "linear",
            vector: Optional[Dseqrecord] = None,
        ) -> List[Dseqrecord]:
            """Golden Gate Assembly 수행.

            - linear: 부품 fragment를 문자열 병합(검증된 기존 경로, 서열 불변).
            - circular: **벡터 backbone이 조립에 참여**한다. 구버전은 부품만으로 원을 닫으려 해
              양 끝 일치 검사가 사실상 통과 불가였고(한 번도 실행된 적 없음), 그래서
              실제 절단 sticky end로 pydna 결합을 쓴다.

            Args:
                fragments: 부품 prepared fragment (슬라이스 표현)
                assembly_type: "linear" | "circular"
                vector: circular일 때 필수 — 조립 효소로 잘라 backbone을 얻는다
            """
            from Bio.Seq import Seq
            
            if assembly_type not in ["circular", "linear"]:
                raise ValueError(f"assembly_type 오류: {assembly_type}")

            if not fragments:
                raise ValueError("조립할 Fragment가 비어 있습니다.")

            try:
                # 첫 번째 조각으로 조립 시작
                combined = fragments[0]

                for i in range(1, len(fragments)):
                    prev_part = combined
                    curr_part = fragments[i]
                    
                    # 물리적 끝단 추출 (Top strand 기준)
                    oh = self.overhang_len   # 효소 파생 — 4 고정 아님(SapI는 3)
                    prev_right = str(prev_part.seq)[-oh:].upper()
                    curr_left = str(curr_part.seq)[:oh].upper()

                    # prep 슬라이스 표현: fragment = [left_oh][body][RC(right_oh)].
                    # 인접 부품은 같은 fusion site를 공유하되 left_oh를 RC로 할당하므로(prepare_fragments),
                    # prev 끝(= RC(right_oh))과 curr 시작(= left_oh)이 **같은 값**이어야 정상이다.
                    # (생물학적으로는 part i가 bottom 5'overhang=right_oh, part i+1이 top 5'overhang=RC을
                    #  노출해 상보결합 — top strand 표현에서는 동일 4bp로 보인다.)
                    if prev_right != curr_left:
                        error_msg = (
                            f"Overhang 불일치: {prev_part.name} 끝({prev_right})과 "
                            f"{curr_part.name} 시작({curr_left})이 일치하지 않습니다."
                        )
                        raise ValueError(error_msg)

                    # 슬라이스 fragment는 sticky-end가 없으므로 string-merge로 병합한다.
                    # 공유 fusion 4bp는 prev 끝(= RC(right_oh) = junction top strand)에서 1회 유지하고
                    # curr 선두 4bp(left_oh = RC(right_oh), 같은 값)를 제거. prep 할당이 인접 overhang을
                    # RC로 잡아주므로(left_oh=RC), prev_right == curr_left 가 성립한다.
                    new_seq = str(prev_part.seq) + str(curr_part.seq)[oh:]
                    combined = Dseqrecord(Dseq(new_seq))
                    combined.name = f"joined_{i}"

                if assembly_type == "circular":
                    result = self._close_with_vector(fragments, vector)
                else:
                    result = combined

                result.name = "assembled_module"
                result.description = f"Golden Gate {assembly_type} assembly"
                return [result]

            except Exception as e:
                if "Overhang 불일치" in str(e):
                    raise e
                raise RuntimeError(f"Assembly 시뮬레이션 실패: {e}")

    def validate(self, result: Dseqrecord) -> Dict:
        """결과물에 조립용 효소 사이트가 남아있지 않은지 검증.

        구버전은 BsaI/BbsI만 검사해서 BsmBI/SapI로 조립하면 아무것도 안 봤고,
        `warnings` 조립식의 연산자 우선순위 버그로 BsaI가 잡히면 BbsI 경고가 사라졌다.
        지원 효소 전체를 검사하도록 정정한다(조립에 쓴 효소 잔류가 가장 중요).
        """
        from Bio import Restriction
        from Bio.Seq import Seq
        from params import SUPPORTED_ENZYMES

        # ⚠️ 원형 산물은 **origin을 걸치는 사이트**를 선형 검색으로 놓친다.
        #    (V0001 실측: 접합부에서 새로 생긴 BsaI가 origin에 걸려 0개로 보고됐다)
        #    2배 서열에서 찾고 원형 좌표로 접어 중복을 제거한다.
        seq_str = str(result.seq).upper()
        n = len(seq_str)
        is_circular = bool(getattr(result.seq, "circular", False))
        target = Seq(seq_str + seq_str) if is_circular else Seq(seq_str)

        def _find(enz_name):
            hits = getattr(Restriction, enz_name).search(target)
            if not is_circular:
                return list(hits)
            return sorted({(h - 1) % n for h in hits})

        warnings = []
        residual = {}
        for name in SUPPORTED_ENZYMES:
            sites = _find(name)
            if sites:
                residual[name] = len(sites)
                mark = " ← 조립에 사용한 효소" if name == self.enzyme_name else ""
                warnings.append(f"{name} 사이트 {len(sites)}개 잔류{mark}")

        return {
            # 조립에 쓴 효소가 결과물에 남아있으면 재절단 위험 → invalid
            'valid': self.enzyme_name not in residual,
            'warnings': warnings,
            'residual_sites': residual,
        }

    def check_residual_sites(self, result: Dseqrecord, label: str = "조립 결과") -> Dict:
        """조립 산물에 효소 사이트가 남았는지 확인하고 **경고를 출력**한다(차단하지 않음).

        Golden Gate 산물은 인식서열이 제거되어 재절단되지 않는 것이 정상이다. 그런데
        **접합부에서 인식서열이 새로 생길 수 있다** — vector overhang과 인접 부품 서열이
        우연히 이어져 만들어지는 경우다(V0001 + 특정 부품 조합에서 실측됨).
        조립 자체는 성립하므로 막지 않되, 다음 단계에서 같은 효소를 쓰면 문제가 되므로 알린다.
        """
        v = self.validate(result)
        if v["residual_sites"]:
            print(f"     ⚠️ {label}에 제한효소 사이트가 남아 있습니다: {v['residual_sites']}")
            if not v["valid"]:
                print(f"        조립에 쓴 {self.enzyme_name}가 재생성되었습니다 — "
                      f"이 산물로 같은 효소를 다시 쓰면 잘립니다(접합부 우연 생성 가능).")
        return v

    # ============================================================================
    # Helper 메서드
    # ============================================================================

    # TODO: 아래 코드의 목적성
    @staticmethod
    def _as_circular(vector: Dseqrecord) -> Dseqrecord:
        """벡터를 **원형으로 정규화**한다. 이미 원형이면 그대로 돌려준다.

        ⚠️ **왜 필요한가 (2026-09-23 실측)**: `DatabaseImport.get_part_by_id`는 vector를
        `circular=False`로 준다(플라스미드인데 선형). 선형으로 자르면 조각 수와 말단이 달라진다.

            V0001 선형 절단 -> 3조각. 최대 조각(2334bp)의 3' 말단이 **blunt**
            V0001 원형 절단 -> 2조각. backbone(4305bp) 양끝이 모두 sticky (5'TTTC / 3'CGGT)

        blunt 말단은 sticky end와 결합할 수 없으므로 GG vector insertion이 **항상 실패**했다.
        덤으로 `extract_vector_overhangs`의 반환값도 선형일 때 엉뚱한 값이 나온다
        (v_right: 선형 CGTC(무의미) vs 원형 ACCG(정답)).

        CP7 circular 경로는 호출부가 이미 `looped()`로 넘기고 있어 영향이 없다(정규화는 무해).
        """
        return vector if getattr(vector, "circular", False) else vector.looped()

    def _cut_with_enzyme(self, dseq: Dseqrecord, enzyme) -> List[Dseqrecord]:
        """제한효소로 DNA를 절단합니다 (pydna 예외 처리 포함)."""
        try:
            # pydna의 cut()은 인식 사이트가 겹치거나 말단에 너무 가까우면 예외를 발생시킬 수 있음
            frags = list(dseq.cut(enzyme))
            for f in frags:
                f.name = dseq.name
            return frags
        except Exception as e:
            # 에러 발생 시 빈 리스트 반환 (상위 extract_insert_fragment에서 처리됨)
            print(f"     ⚠️  [{dseq.name}] 절단 시뮬레이션 예외 발생: {e}")
            return []

    def _cut_with_bsai(self, dseq: Dseqrecord) -> List[Dseqrecord]:
        """하위 호환성 유지용"""
        return self._cut_with_enzyme(dseq, BsaI)

    def _close_with_vector(self, fragments: List[Dseqrecord], vector: Optional[Dseqrecord]) -> Dseqrecord:
        """부품 fragment들을 **벡터 backbone과 함께** 원형으로 결합 (circular 조립).

        문자열 병합이 아니라 pydna sticky-end ligation을 쓴다. 문자열로는 overhang의
        가닥(상보성)을 따질 수 없어 backbone 접점에서 방향이 어긋나기 때문이다.

        Raises:
            ValueError: vector 미지정
            RuntimeError: sticky end가 호환되지 않아 결합 실패
        """
        if vector is None:
            raise ValueError(
                "circular 조립에는 vector가 필요합니다 — module이 들어갈 중간 벡터를 지정하세요."
            )
        backbone = self._get_backbone(vector)
        sticky = [self.to_sticky(f) for f in fragments]

        acc = backbone
        for i, frag in enumerate(sticky):
            try:
                acc = acc + frag
            except Exception as e:
                raise RuntimeError(
                    f"circular 조립 실패: {i + 1}번째 부품({frag.name})의 sticky end가 "
                    f"직전 조각과 맞지 않습니다. overhang 배정을 확인하세요. ({e})"
                )
        try:
            result = acc.looped()
        except Exception as e:
            raise RuntimeError(
                f"circular 조립 실패: 마지막 부품과 vector backbone이 닫히지 않습니다. "
                f"벡터 overhang과 체인 양 끝 배정을 확인하세요. ({e})"
            )
        result.name = "assembled_module"
        result.description = f"Golden Gate circular assembly (vector={vector.name})"
        return result

    def _get_backbone(self, vector: Dseqrecord) -> Dseqrecord:
        """벡터를 조립 효소로 잘라 backbone(가장 긴 조각)을 얻는다.

        ⚠️ 원형으로 정규화한 뒤 자른다 — 선형으로 자르면 조각이 하나 더 생기고
        최대 조각의 한쪽 말단이 blunt가 되어 sticky-end ligation이 불가능해진다
        (`_as_circular` 참조).
        """
        vector = self._as_circular(vector)
        frags = self._cut_with_enzyme(vector, self.enzyme)
        if not frags:
            raise ValueError(
                f"vector '{vector.name}'를 {self.enzyme_name}로 자를 수 없어 circular 조립이 "
                f"불가능합니다. {self.enzyme_name} 사이트를 가진 중간 벡터가 필요합니다."
            )
        return max(frags, key=len)

    def to_sticky(self, frag: Dseqrecord) -> Dseqrecord:
        """양 끝 overhang을 한 줄에 편 문자열 표현 → 실제 5' 돌출을 가진 Dseqrecord.

        `extract_insert_fragment`는 pydna 절단 실패(짧은 부품)를 피하려고 결정적 슬라이스를
        쓴다. 그 결과는 `[left_oh][core][RC(right_oh)]`를 **한 가닥에 편** 표현이라 blunt이고,
        sticky-end ligation(circular 조립)에는 그대로 쓸 수 없다.

        실제 절단 산물은 각 가닥이 (전체 길이 − overhang)이고 양 끝에 5' 돌출이 있다.
        아래 변환이 pydna 절단 결과와 **완전히 일치**함을 확인했다(watson/crick/ovhg 모두).

        Returns:
            5' 돌출을 가진 Dseqrecord. overhang보다 짧으면 원본을 그대로 돌려준다.
        """
        oh = self.overhang_len
        seq = str(frag.seq).upper()
        if len(seq) <= 2 * oh:
            return frag
        sticky = Dseq(watson=seq[:-oh], crick=_rc(seq[oh:]), ovhg=-oh)
        out = Dseqrecord(sticky)
        out.name = frag.name
        out.annotations = dict(getattr(frag, "annotations", {}) or {})
        return out

    def extract_vector_overhangs(self, vector: Dseqrecord) -> Tuple[str, str]:
        """Vector를 **선택된 효소**로 자른 후 노출되는 오버행을 추출합니다.

        (구버전은 `self.enzyme`과 무관하게 항상 BsaI로 잘랐다 — BsmBI 등을 골라도 BsaI 결과가
         나오는 잠복 버그였고, 효소를 유저가 고르게 되면서 실제 문제가 되므로 정정했다.)

        Raises:
            ValueError: vector에 인식 사이트가 없어(= 이 vector로는 Golden Gate 조립 불가)
                절단 결과가 비었거나, backbone이 overhang보다 짧은 경우.
        """
        # 플라스미드는 원형이다. 선형으로 자르면 조각 수·말단이 달라져 조용히 틀린
        # overhang이 나가고, backbone 한쪽이 blunt가 되어 ligation이 불가능해진다.
        vector = self._as_circular(vector)
        frags = self._cut_with_enzyme(vector, self.enzyme)

        # frags가 비면 아래 max()가 "ValueError: max() iterable argument is empty"로 죽는다.
        # 원인이 전혀 드러나지 않는 에러라, 실제 원인을 그대로 알린다.
        # _cut_with_enzyme는 사이트가 없을 때뿐 아니라 **pydna 절단 예외도 삼키고 []를 반환**하므로
        # (사이트가 겹치거나 말단에 너무 가까운 경우) 두 경우를 함께 안내한다.
        if not frags:
            raise ValueError(
                f"vector '{vector.name}'를 {self.enzyme_name}로 자를 수 없어 Golden Gate 조립이 "
                f"불가능합니다. {self.enzyme_name} 인식 서열이 없거나, 사이트가 겹치거나 말단에 "
                f"너무 가까워 절단 시뮬레이션이 실패한 경우입니다"
                f"(후자면 바로 위 '절단 시뮬레이션 예외' 로그를 확인하세요). "
                f"→ {self.enzyme_name} 사이트를 가진 vector를 고르거나 다른 효소/조립 방식을 쓰세요."
            )

        backbone = max(frags, key=len)
        full_seq = str(backbone.seq).upper()

        # backbone이 overhang **2배**보다 짧으면 앞/뒤 슬라이스가 겹치거나 같은 조각을
        # 돌려주며 에러 없이 **조용히 틀린** 오버행이 나간다 → 명시적으로 막는다.
        # ⚠️ 기준이 `< oh`였는데 그건 부족하다: backbone이 정확히 oh면 [:oh]와 [-oh:]가
        #    **같은 서열**이 되어 좌우 overhang이 동일하게 나간다(2026-09-23, 원형 정규화로 드러남).
        #    좌우가 구분되려면 최소 2*oh가 필요하다.
        oh = self.overhang_len
        if len(full_seq) < 2 * oh:
            raise ValueError(
                f"vector '{vector.name}'의 backbone이 {len(full_seq)}bp로 너무 짧아 "
                f"{self.enzyme_name} overhang {oh}bp를 좌우로 구분해 추출할 수 없습니다 "
                f"(최소 {2 * oh}bp 필요)."
            )

        left, right = full_seq[:oh], full_seq[-oh:]

        # ── 좌우 overhang이 같으면 방향성 클로닝이 성립하지 않는다 ──────────────
        # 인식 사이트가 **1개뿐인** 벡터가 전형적이다. 원형을 한 번 자르면 양 끝이 같은
        # 사이트에서 나오므로 좌우 overhang이 동일해진다. 그러면
        #   · insert가 양방향으로 들어가고(방향 제어 불가)
        #   · 벡터가 insert 없이 자기들끼리 다시 닫힌다.
        # 조용히 진행하면 spacer 선정이 수천 번 조립 dry-run을 돌다 실패하므로
        # (2026-09-23 실측 8분+) 여기서 원인을 밝혀 즉시 막는다.
        #
        # DB 실측(387개, BsaI 기준): 사이트 0개 157 / **1개 151(전부 좌우 동일)** /
        # 2개 50 / 3개 이상 29 → 실제로 GG 삽입에 쓸 수 있는 벡터는 79개(20%)뿐이다.
        if left == right:
            raise ValueError(
                f"vector '{vector.name}'는 {self.enzyme_name} 절단 시 좌우 overhang이 "
                f"모두 '{left}'로 같아 Golden Gate 삽입에 쓸 수 없습니다 "
                f"(방향성 클로닝 불가 · 벡터 self-ligation). "
                f"{self.enzyme_name} 사이트가 **2곳 이상**인 vector를 고르거나 Gibson을 쓰세요."
            )
        return left, right

__all__ = ['GoldenGateModule']
