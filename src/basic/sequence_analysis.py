"""
서열 분석 기초 기능

이 모듈은 프로젝트 전체에서 재사용되는 DNA 서열 분석 기능을 제공합니다.

핵심 원칙:
    - DRY (Don't Repeat Yourself): 각 기능은 이 파일에서만 구현
    - 라이브러리 활용: BioPython, pytrf 등 검증된 라이브러리 사용
    - 모든 GC%, Tm 계산은 여기서만!
    - 다른 모듈에서는 이 모듈을 import해서 사용

작성일: 2025-10-30
수정일: 2025-10-30 (라이브러리 기반으로 재작성)
"""

import re
from typing import Union, List, Dict, Optional

# BioPython imports
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Restriction import Restriction

# pytrf import
import pytrf


def get_typeIIS_params(enzyme_name: str) -> Dict[str, object]:
    """Type IIS 제한효소의 설계 파라미터를 Biopython에서 동적 도출 (B-1).

    하드코딩 enzyme_sites 딕셔너리 대체. elucidate()(예: 'GGTCTCN^NNNN_N')의
    site/`^`(절단점)/`_`(상보가닥 절단점) 위치를 파싱해 설계 파라미터를 산출한다.
    단일 소스(Biopython)로 BsaI/BbsI/BsmBI/SapI 등 Type IIS를 일관 처리.

    Args:
        enzyme_name: Type IIS 효소명 (예: "BsaI", "BbsI")

    Returns:
        {'fwd': 인식서열, 'rev': 역상보, 'spacer_len': 인식~절단 더미 nt 수,
         'def_spacer': 'A'*spacer_len, 'overhang_len': 점착말단 길이}
    """
    enz = getattr(Restriction, enzyme_name)
    site = str(enz.site)
    el = enz.elucidate()
    spacer_len = el.index('^') - len(site)
    overhang_len = el.index('_') - el.index('^') - 1
    return {
        'fwd': site,
        'rev': str(Seq(site).reverse_complement()),
        'spacer_len': spacer_len,
        'def_spacer': 'A' * spacer_len,
        'overhang_len': overhang_len,
    }


def validate_dna_sequence(sequence: str) -> bool:
    """
    DNA 서열 유효성 검증

    검증 항목 (E-1):
        - 빈 서열 체크
        - ATGC만 허용 (N·IUPAC 모호코드 R/Y/S/W 등 모두 거부)
          → 무오류 우선: 모호코드를 받으면 절단/상보성/Tm 등에서 예외 분기가 폭증하므로
            1차 배포에서는 ATGC-only로 엄격히 거부. (길이 기반 거부는 추후 별도.)

    Args:
        sequence: DNA 서열

    Returns:
        유효하면 True, 아니면 False

    Example:
        >>> validate_dna_sequence("ATCG")
        True
        >>> validate_dna_sequence("ATCGN")
        False
        >>> validate_dna_sequence("")
        False
    """
    if not sequence:
        return False

    # 대소문자 무관, ATGC만 허용 (N 제거 — E-1)
    valid_bases = set('ATCGatcg')
    return all(base in valid_bases for base in sequence)


def validate_cds(sequence: str) -> dict:
    """
    CDS 서열 유효성 검증

    검증 항목:
        1. ATG start codon 존재 여부
        2. 인프레임 종결코돈(TAA, TAG, TGA) 존재 여부
        3. ORF 길이가 3의 배수인지 (ATG ~ 종결코돈 끝)
        4. ORF 길이가 전체 서열의 70% 이상인지 
        # TODO: 4번째 체크항목의 evidence 확인

    Args:
        sequence: CDS DNA 서열

    Returns:
        {
            'has_start_codon': bool,   # ATG 존재 여부
            'start_pos': int,          # ATG 위치 (0-based), 없으면 -1
            'has_stop_codon': bool,    # 인프레임 종결코돈 존재 여부
            'stop_codon': str,         # 발견된 종결코돈 ('TAA'/'TAG'/'TGA'), 없으면 None
            'stop_pos': int,           # 종결코돈 시작 위치 (0-based), 없으면 -1
            'is_valid_length': bool,   # ORF 길이 % 3 == 0
            'is_dominant_orf': bool,   # ORF 길이 >= 전체 서열의 70%
            'orf_ratio': float,        # ORF / 전체 서열 비율 (디버깅용), 없으면 0.0
            'is_valid_cds': bool       # 위 4가지 모두 통과 여부
        }

    Example:
        >>> validate_cds("ATGAAATAA")
        {'has_start_codon': True, 'start_pos': 0, 'has_stop_codon': True,
         'stop_codon': 'TAA', 'stop_pos': 6, 'is_valid_length': True,
         'is_dominant_orf': True, 'orf_ratio': 1.0, 'is_valid_cds': True}
        >>> validate_cds("AAATAA")
        {'has_start_codon': False, 'start_pos': -1, 'has_stop_codon': False,
         'stop_codon': None, 'stop_pos': -1, 'is_valid_length': False,
         'is_dominant_orf': False, 'orf_ratio': 0.0, 'is_valid_cds': False}
    """
    STOP_CODONS = {"TAA", "TAG", "TGA"}

    seq = sequence.upper()

    # 1. ATG start codon 확인
    start_pos = seq.find("ATG")
    has_start = start_pos != -1

    # 2. 인프레임 종결코돈 확인 (ATG 기준 reading frame 유지, 3 step씩 탐색)
    stop_codon = None
    stop_pos = -1

    if has_start:
        sub_seq = seq[start_pos:]
        for i in range(0, len(sub_seq) - 2, 3):
            codon = sub_seq[i:i + 3]
            if codon in STOP_CODONS:
                stop_codon = codon
                stop_pos = start_pos + i
                break

    has_stop = stop_codon is not None

    # 3. ORF 길이 검증 (ATG 시작 ~ 종결코돈 끝, 3의 배수)
    is_valid_length = False
    orf_ratio = 0.0

    if has_start and has_stop:
        orf_length = (stop_pos - start_pos) + 3  # 종결코돈 포함
        is_valid_length = (orf_length % 3 == 0)

        # 4. ORF 비율 검증 (전체 서열의 70% 이상)
        orf_ratio = orf_length / len(seq) if len(seq) > 0 else 0.0

    is_dominant_orf = orf_ratio >= 0.7

    return {
        'has_start_codon': has_start,
        'start_pos': start_pos,
        'has_stop_codon': has_stop,
        'stop_codon': stop_codon,
        'stop_pos': stop_pos,
        'is_valid_length': is_valid_length,
        'is_dominant_orf': is_dominant_orf,
        'orf_ratio': round(orf_ratio, 3),
        'is_valid_cds': has_start and has_stop and is_valid_length and is_dominant_orf
    }


def calculate_gc_content(
    sequence: str,
    window: Optional[int] = None
) -> Union[float, List[float]]:
    """
    GC 함량 계산 (BioPython 기반)

    ⭐ BioPython의 gc_fraction() 함수를 사용합니다.
    ⭐ 프로젝트 전체에서 재사용!

    Args:
        sequence: DNA 서열
        window: 윈도우 크기
            - None: 전체 서열의 GC%
            - int: 슬라이딩 윈도우로 계산

    Returns:
        - window=None: 전체 GC% (0.0 ~ 1.0)
        - window=int: 각 윈도우의 GC% 리스트

    Raises:
        ValueError: 빈 서열이거나 윈도우 크기가 서열보다 클 때

    Example:
        >>> calculate_gc_content("ATCG")
        0.5
        >>> calculate_gc_content("GGCC")
        1.0
        >>> calculate_gc_content("AATT")
        0.0
        >>> calculate_gc_content("ATCGATCG", window=4)
        [0.5, 0.5, 0.5, 0.5, 0.5]

    Notes:
        - BioPython의 gc_fraction()은 이미 0-1 범위를 반환
        - Legacy 코드와 호환성 유지
    """
    if not sequence:
        raise ValueError("The sequence is empty")

    sequence = sequence.upper()

    # 전체 GC% 계산 (BioPython 사용)
    if window is None:
        # BioPython gc_fraction은 이미 0-1 범위를 반환
        return gc_fraction(sequence)

    # 윈도우 기반 GC% 계산
    if window > len(sequence):
        raise ValueError(
            f"The window size ({window}) is larger than the sequence length ({len(sequence)})"
        )

    gc_percentages = []
    for i in range(len(sequence) - window + 1):
        window_seq = sequence[i:i + window]
        # BioPython gc_fraction을 각 윈도우에 적용
        gc_percentages.append(gc_fraction(window_seq))

    return gc_percentages


def calculate_tm(
    sequence: str,
    Na: float = 0.0,
    dnac1: float = 500.0,
    K: float = 50.0,
    Tris: float = 10.0,
    Mg: float = 1.5,
    dNTPs: float = 0.8,
    saltcorr: int = 6,
) -> float:
    """
    Tm (melting temperature) 계산 (BioPython 기반)

    ⭐ BioPython의 Tm_NN() 함수를 사용합니다 (Nearest-Neighbor method).
    ⭐ 프로젝트 전체에서 재사용!

    ⚠️ **NEB Tm Calculator를 재현한다** (2026-11-05). 아래 "왜 NEB 축인가" 참고.

    Args:
        sequence: DNA 서열
        Na: Na+ 농도 (mM). NEB 버퍼는 K+/Tris 기반이라 **0이 기본**이다
        dnac1: 프라이머 농도 (nM, 기본 500 = Q5 프로토콜 0.5 µM)
        K: K+ 농도 (mM)
        Tris: Tris 농도 (mM) — 1가 염 환산에 Tris/2가 더해진다
        Mg: Mg2+ 농도 (mM). **Taq 경로는 0** — NEB가 Taq에 Mg 보정을 하지 않는다
        dNTPs: dNTP 총 농도 (mM) — Mg2+를 킬레이트해 유효 Mg를 낮춘다
        saltcorr: 염보정 방식 (6=Owczarzy, 1=Schildkraut). **Phusion만 1**
    Returns:
        Tm (°C)

    Raises:
        ValueError: 빈 서열

    Example:
        >>> tm = calculate_tm("ATCGATCGATCGATCG")
        >>> 50 < tm < 60  # 대략적인 범위
        True
        >>> tm = calculate_tm("GGCCGGCCGGCC")
        >>> 60 < tm < 75  # GC-rich는 높은 Tm
        True
        >>> tm = calculate_tm("ATCGATCGATCGATCG", Na=100, dnac1=500)
        >>> 52 < tm < 62  # 높은 salt 농도는 Tm 증가
        True

    Notes:
        ## 왜 NEB 축인가 — **자가 하나여야 하고, 그 자는 랩이 쓰는 자여야 한다**

        이 함수는 하루에 두 번 바뀌었다. 기록을 남긴다.

        ### 1단계 (2026-11-05 오전): Biopython → primer3

        같은 `tm_min`이 **경로마다 다른 자로 읽히고 있었다**:
          · gibson / yeast_hr → primer3에 `PRIMER_MIN_TM`으로 그대로 전달 → primer3 축
          · golden_gate       → `_fit_amplicon`이 `calculate_tm`과 직접 비교 → Biopython 축
        두 자는 같은 염·농도에서도 +4.27℃ 어긋났다(실측 3,100 표본). 그래서 **프라이머를
        고르는 주체**인 primer3의 자로 통일했다.

        ### 2단계 (2026-11-05 오후): primer3 → **NEB Tm Calculator 재현**

        통일은 됐지만 그 자가 **랩이 쓰는 자가 아니었다.** 연구실은 NEB 계산기로 Ta를
        잡는데, 우리 숫자는 거기서 몇 도씩 달랐다. 게다가 polymerase 선택이 설계에
        거의 반영되지 않았다 — 프리셋 셋이 전부 `na_mM=50`이라 같은 프라이머의 Tm이
        1.1℃밖에 안 움직였다.

        NEB는 **제품마다 다른 Tm**을 낸다(버퍼가 다르므로). 그런데 Q5/Phusion 버퍼
        조성은 **영업비밀로 비공개**라 오랫동안 오프라인 재현이 불가능하다고 알려져 있었다.

        ### 해법과 증빙

        pydna 이슈 #237(https://github.com/pydna-group/pydna/issues/237)에서
        **Biopython `MeltingTemp`의 저자 본인(MarkusPiotrowski)**이 맞는 인자를 찾아 올렸다.
        핵심 통찰:
          · 버퍼 조성을 **몰라도 된다** — Na + K + Tris/2 가 하나의 1가 염 농도로 합산되므로
            비공개인 Q5/Phusion도 K/Tris 조합으로 같은 유효 농도를 만들면 값이 맞는다.
          · **NEB는 Taq에 Mg2+ 보정을 하지 않는다** → Taq만 `Mg=0`
            (Mg=0이면 Biopython은 dNTPs도 무시한다).
          · **Phusion만 Schildkraut 염보정**(`saltcorr=1`), 나머지는 Owczarzy(`6`).
            NEB 자신의 API 문서와 일치한다 — "For Phusion DNA Polymerases, the salt
            correction of Schildkraut is used."

        조건값은 `params.NEB_TM_CONDITIONS`에 한 벌만 둔다(프리셋도 거기서 끌어온다).

        ### 자체 검증 (2026-11-05, NEB Tm API v0.19.0)

        실제 Bioparts 파츠 말단 160개(길이 18~36nt)를 NEB API와 대조:

            제품       RMSE     최대오차   ±1℃ 내   Ta 일치
            Taq        0.29℃    0.50      100%     100%
            Phusion    0.62℃    1.45       88%      95%
            Q5         0.35℃    0.72      100%     100%

        ⚠️ **재검증할 때의 함정**: NEB API의 `conc` 파라미터는 **nM이 아니라 µM**이다.
           nM로 알고 500을 넣으면 1000배 과대평가라 Tm이 ~8℃ 높게 나오고, 그걸
           "NEB 고유 모델"로 오해하게 된다(2026-11-05에 실제로 겪었다).
           Q5/Phusion은 `conc=0.5`, Taq은 `conc=0.2`로 보낼 것.

        ⚠️ 재현이 안 되는 구현을 쓰지 말 것: `polymerase-tm`(PyPI)은 "0℃ 오차 재현"을
           표방하지만 Owczarzy에 **길이 항이 없어** 길이가 길수록 벌어진다(잔차 길이상관
           −0.95). 반대로 primer3/SantaLucia는 길이 항이 과해 +0.95로 뒤집힌다.
           Biopython의 위 설정만이 양쪽 다 상관이 사라진다(Taq −0.08, Q5 +0.02).

        ## 비-ACGT는 **조용히 넘기지 않는다**

        모호 염기가 든 프라이머는 어차피 주문할 수 없고, 여기서 다른 방식으로 폴백하면
        그 서열만 다른 자로 재게 되어 "자가 두 개"가 되살아난다. 그래서 명시적으로
        실패시킨다. Bioparts 11,738종은 전부 ACGT라 실제 경로엔 영향이 없고,
        서열 통계처럼 설계 판단이 아닌 자리는 호출부가 감싼다.

    """
    if not sequence:
        raise ValueError("The sequence is empty")

    sequence = sequence.upper()
    bad = set(sequence) - set("ACGT")
    if bad:
        raise ValueError(
            f"Tm is computed on ACGT only - disallowed bases {sorted(bad)} are present. "
            f"(A primer containing ambiguous bases cannot be ordered.)"
        )

    # ⚠️ `dnac2=0`이 **중요하다.** Biopython은 비-자기상보 서열에서 유효 농도를
    #    `dnac1 - dnac2/2`로 쓴다. dnac2를 dnac1과 같게 주면 농도가 절반이 되어
    #    NEB와 어긋난다(2026-11-05에 이것 때문에 한참 헤맸다).
    return mt.Tm_NN(
        Seq(sequence),
        nn_table=mt.DNA_NN3,      # SantaLucia 1998 unified (Allawi & SantaLucia 1997)
        Na=Na, K=K, Tris=Tris, Mg=Mg, dNTPs=dNTPs,
        dnac1=dnac1, dnac2=0,
        saltcorr=saltcorr,
    )


# ── 효소 목록 (단일 정의) ────────────────────────────────────────────────────
# ⚠️ 선형화 효소를 여기 말고 다른 곳에 또 적지 않는다.
#    2026-10-20까지 `assembly/gibson._linearize_vector`가 자체 목록 5종
#    (SmaI/EcoRV/PvuII/StuI/DraI)을 하드코딩하고 있었고, 그중 4종은 이 파일에
#    **있지도 않았다.** 목록이 두 벌로 갈라진 결과 DB 벡터 387개 중 79%만
#    선형화가 됐다(아래 목록을 쓰면 99.5%). 동료가 "사용 불가"로 막힌 pACYC
#    벡터도 실제로는 쓸 수 있는 효소가 7개 있었다.

# Golden Gate용 Type IIS (인식서열 바깥을 자름)
GOLDEN_GATE_ENZYMES = [
    'BsaI', 'BbsI', 'BsmBI',
    'Esp3I',   # BsmBI 동의효소
    'PaqCI',   # AarI 동의효소, 7bp 인식 → unintended cut 확률 낮음
    'SapI',    # 7bp 인식, 3bp overhang, tRNA scaffold 클로닝
]

# Gibson / Yeast HR 벡터 선형화용. **정확히 한 곳만 자르는 것**을 골라 쓴다.
# blunt/sticky를 섞어 둔 이유: Gibson은 말단을 깎아내므로 양쪽 다 쓸 수 있고,
# 후보가 넓을수록 쓸 수 있는 벡터가 늘어난다(실측 5종 79% → 19종 99.5%).
LINEARIZATION_ENZYMES = [
    # MCS에 흔한 6bp 회문 (sticky)
    'EcoRI', 'BamHI', 'HindIII', 'XhoI', 'SalI', 'PstI', 'XbaI',
    'SpeI',    # BioBrick 표준
    'NheI',    # SpeI compatible overhang
    'AvrII',   # NheI/SpeI compatible
    'KpnI', 'ClaI', 'SphI',
    'NotI',    # 8bp 인식, pET/pUC 계열에 범용
    # blunt cutter
    'SmaI', 'EcoRV', 'PvuII', 'StuI', 'DraI',
]


def usable_linearization_enzymes(vector_seq: str, circular: bool = True) -> List[str]:
    """이 벡터를 **정확히 한 곳에서** 자르는 선형화 효소 목록.

    Gibson/Yeast HR 삽입은 벡터를 한 번만 잘라 열어야 한다. 두 번 이상 자르면
    조각이 나고, 0번이면 열리지 않는다. 유저에게 선택지를 보여줄 때 쓴다.
    """
    seq = Seq(vector_seq.upper())
    out = []
    for name in LINEARIZATION_ENZYMES:
        try:
            if len(getattr(Restriction, name).search(seq, linear=not circular)) == 1:
                out.append(name)
        except Exception:
            continue
    return out


# IUPAC 모호코드 → 정규식. 현재 쓰는 효소들은 모두 비모호 서열이지만,
# 사용자가 Biopython의 임의 효소명을 넘길 수 있으므로 방어적으로 처리한다.
_IUPAC_TO_REGEX = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
    'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}


def _recognition_positions(sequence: str, site: str) -> List[int]:
    """서열에서 **인식서열이 나타나는 0-based 시작 위치**를 양쪽 가닥에서 찾는다.

    ⚠️ Biopython `Restriction.search()`를 쓰면 안 된다 — 그건 **절단 위치**를 돌려준다.
       Type IIS는 인식서열 **바깥**을 자르므로, 서열 끝 근처의 인식서열은 절단점이
       서열 밖으로 나가 **탐지되지 않는다**. 단독 조각만 보면 그 판정이 맞지만,
       부품은 옆 부품과 이어진 뒤에 잘리므로 우리가 묻는 건 "지금 자를 수 있나"가
       아니라 **"인식서열이 있나"**다.
       (2026-10-10 실측: BsaI/BsmBI는 양 끝 11nt, BbsI는 12nt가 맹점이었고
        운영 Bioparts 11227개 중 27개가 "깨끗함"으로 오판되고 있었다.)

    역가닥은 인식서열의 역상보를 정방향 서열에서 찾는 것과 같다. 회문(EcoRI 등)은
    두 패턴이 같으므로 중복을 제거한다. 겹치는 반복도 놓치지 않도록 lookahead를 쓴다.
    """
    seq = sequence.upper()
    site = site.upper()
    patterns = {site, str(Seq(site).reverse_complement()).upper()}

    found = set()
    for pat in patterns:
        regex = "".join(_IUPAC_TO_REGEX.get(ch, re.escape(ch)) for ch in pat)
        # (?=...) — 겹쳐 나오는 경우도 전부 잡는다.
        for m in re.finditer(f"(?=({regex}))", seq):
            found.add(m.start())
    return sorted(found)


def find_restriction_sites(
    sequence: str,
    enzyme: Optional[str] = None
) -> Dict[str, List[int]]:
    """
    제한효소 인식 서열 탐색

    BioPython의 Restriction 모듈을 사용하여 제한효소 사이트를 찾습니다.

    Args:
        sequence: DNA 서열
        enzyme: 특정 효소명 (None이면 주요 효소들 검색)
            예: "BsaI", "BbsI", "EcoRI", "BamHI" 등

    Returns:
        {enzyme_name: [position1, position2, ...]}
        위치는 **인식서열이 시작하는 0-based index**(절단 위치가 아니다).
        정방향·역방향 가닥을 모두 찾으며, 역가닥 히트는 역상보 서열이 나타나는 위치로 보고한다.

    Raises:
        ValueError: 빈 서열이거나 알 수 없는 효소

    Example:
        >>> seq = "ATCGGGTCTCATCG"
        >>> find_restriction_sites(seq, "BsaI")
        {"BsaI": [4]}

        >>> find_restriction_sites(seq)  # 여러 효소
        {"BsaI": [4], "EcoRI": [], ...}
    """
    if not sequence:
        raise ValueError("The sequence is empty")

    sequence = sequence.upper()
    bio_seq = Seq(sequence)

    main_enzymes = {name: getattr(Restriction, name)
                    for name in (GOLDEN_GATE_ENZYMES + LINEARIZATION_ENZYMES)}

    # 특수: 특정 효소만 검색
    if enzyme:
        lookup = {name.upper(): (name, obj) for name, obj in main_enzymes.items()}
        enzyme_upper = enzyme.upper()

        if enzyme_upper in lookup:
            _, enzyme_obj = lookup[enzyme_upper]
        else:
            try:
                enzyme_obj = getattr(Restriction, enzyme)
            except AttributeError:
                raise ValueError(f"Unknown restriction enzyme: {enzyme}")

        return {enzyme: _recognition_positions(sequence, enzyme_obj.site)}

    # 기본: 전체 검색
    return {name: _recognition_positions(sequence, obj.site) for name, obj in main_enzymes.items()}


def check_restriction_site_in_part_sequence(
    sequence: str,
    enzymes: List[str] = None
) -> dict:
    """
    Part 서열 내부의 제한효소 인식 서열 검증

    Golden Gate 조립에 사용하는 제한효소(BsaI, BbsI 등)의 인식 서열이
    part 내부에 존재하면 조립 시 의도치 않은 절단이 발생하므로 에러로 처리.

    find_restriction_sites()를 내부 엔진으로 사용하며,
    "발견됨 = 나쁜 것"으로 해석하는 part 검증 전용 래퍼 함수.

    Args:
        sequence: Part DNA 서열
        enzymes: 체크할 제한효소 목록 (None이면 기본값 사용)
            기본값: ["BsaI", "BbsI"] (Golden Gate에서 사용하는 효소들)

    Returns:
        {
            'is_clean': bool,          # 모든 효소에서 site 없음 = True (good)
            'found_sites': {           # 발견된 효소별 위치 목록
                'BsaI': [pos1, ...],
                'BbsI': []
            }
        }

    Example:
        >>> check_restriction_site_in_part_sequence("ATCGGGTCTCATCG")
        {'is_clean': False, 'found_sites': {'BsaI': [4], 'BbsI': []}}
        >>> check_restriction_site_in_part_sequence("ATCGATCGATCG")
        {'is_clean': True, 'found_sites': {'BsaI': [], 'BbsI': []}}
    """
    if enzymes is None:
        enzymes = ["BsaI", "BbsI"]

    found_sites = {}
    for enzyme in enzymes:
        sites = find_restriction_sites(sequence, enzyme)
        found_sites[enzyme] = sites.get(enzyme, [])

    is_clean = all(len(positions) == 0 for positions in found_sites.values())

    return {
        'is_clean': is_clean,
        'found_sites': found_sites
    }


# 임시로 남겨놓음 -> 잘못된 파라미터 호출
'''
def find_repeats(
    sequence: str,
    min_length: int = 20,
    max_length: int = 100,
    min_copy_number: int = 2
) -> List[Dict]:
    """
    Tandem repeat 검출

    pytrf 라이브러리를 사용하여 tandem repeat를 검출합니다.

    Args:
        sequence: DNA 서열
        min_length: 최소 repeat 길이 (bp)
        max_length: 최대 repeat 길이 (bp)
        min_copy_number: 최소 반복 횟수

    Returns:
        [
            {
                'motif': str,          # Repeat 모티프 서열
                'length': int,         # Repeat 총 길이 (bp)
                'copy_number': int,    # 반복 횟수
                'start': int,          # 시작 위치 (0-based)
                'end': int             # 끝 위치 (0-based)
            },
            ...
        ]

    Raises:
        ValueError: 빈 서열

    Example:
        >>> seq = "ATCGATCGATCGATCG"
        >>> find_repeats(seq, min_length=8)
        [{'motif': 'ATCG', 'length': 16, 'copy_number': 4, ...}]
    """
    if not sequence:
        raise ValueError("The sequence is empty")

    sequence = sequence.upper()

    # pytrf로 repeat 검출
    try:
        repeats = pytrf.STRFinder(
            sequence,
            min_motif_length=1,
            max_motif_length=max_length
        )
    except Exception:
        # pytrf 에러 발생 시 빈 리스트 반환
        return []

    # 결과 필터링 및 포매팅
    results = []
    for repeat in repeats:
        # 최소 길이 체크
        if repeat.length < min_length:
            continue

        # 최소 반복 횟수 체크
        if repeat.copy_number < min_copy_number:
            continue

        results.append({
            'motif': repeat.motif,
            'length': repeat.length,
            'copy_number': repeat.copy_number,
            'start': repeat.start,
            'end': repeat.end
        })

    return results
'''

def find_repeats(
    sequence: str,
    min_length: int = 20,
    max_length: int = 100,
    min_copy_number: int = 2
) -> List[Dict]:
    """
    Tandem repeat 검출

    pytrf 라이브러리를 사용하여 tandem repeat를 검출합니다.

    Args:
        sequence: DNA 서열
        min_length: 최소 repeat 길이 (bp)
        max_length: 최대 repeat 길이 (bp)
        min_copy_number: 최소 반복 횟수
        # TODO: min=20, max=100, copy_nubmer=2로 지정한 기본값의 evidence 확인

    Returns:
        [
            {
                'motif': str,          # Repeat 모티프 서열
                'length': int,         # Repeat 총 길이 (bp)
                'copy_number': int,    # 반복 횟수
                'start': int,          # 시작 위치 (0-based)
                'end': int             # 끝 위치 (0-based)
            },
            ...
        ]

    Raises:
        ValueError: 빈 서열

    Example:
        >>> seq = "ATCGATCGATCGATCG"
        >>> find_repeats(seq, min_length=8)
        [{'motif': 'ATCG', 'length': 16, 'copy_number': 4, ...}]
    """
    if not sequence:
        raise ValueError("The sequence is empty")

    sequence = sequence.upper()

    # pytrf로 repeat 검출
    try:
        # pytrf.STRFinder는 (chrom, seq, mono, di, tri, tetra, penta, hexa) 파라미터를 사용
        # 각 motif 길이별 최소 반복 횟수 설정
        repeats = pytrf.STRFinder(
            "seq",  # 서열 이름 (임의)
            sequence,  # DNA 서열
            mono=min_copy_number,   # 1bp 반복
            di=min_copy_number,     # 2bp 반복
            tri=min_copy_number,    # 3bp 반복
            tetra=min_copy_number,  # 4bp 반복
            penta=min_copy_number,  # 5bp 반복
            hexa=min_copy_number    # 6bp 반복
        )
    except Exception:
        # pytrf 에러 발생 시 빈 리스트 반환
        return []

    # 결과 필터링 및 포매팅
    results = []
    for rep in repeats:
        # 최소 길이 체크
        if rep.length < min_length:
            continue

        # 최소 반복 횟수 체크
        if rep.repeat < min_copy_number:
            continue

        # motif 길이가 max_length를 초과하면 스킵
        if rep.type > max_length:
            continue

        results.append({
            'motif': rep.motif,
            'length': rep.length,
            'copy_number': rep.repeat,  # pytrf에서는 'repeat' 속성
            'start': rep.start - 1,  # pytrf는 1-based, 0-based로 변환
            'end': rep.end - 1  # pytrf는 1-based, 0-based로 변환
        })

    return results

def get_sequence_stats(sequence: str) -> Dict:
    """
    서열 통계 정보 반환

    여러 분석을 한 번에 수행하여 통계 정보를 제공합니다.

    Args:
        sequence: DNA 서열

    Returns:
        {
            'length': int,
            'gc_content': float,
            'tm': float,
            'base_composition': {'A': int, 'T': int, 'C': int, 'G': int},
            'restriction_sites': Dict[str, List[int]],
            'has_repeats': bool
        }

    Example:
        >>> stats = get_sequence_stats("ATCGATCG")
        >>> stats['length']
        8
        >>> stats['gc_content']
        0.5
    """
    if not sequence:
        raise ValueError("The sequence is empty")

    sequence = sequence.upper()

    # 염기 조성
    base_composition = {
        'A': sequence.count('A'),
        'T': sequence.count('T'),
        'C': sequence.count('C'),
        'G': sequence.count('G'),
        'N': sequence.count('N')
    }

    return {
        'length': len(sequence),
        'gc_content': calculate_gc_content(sequence),
        # 서열 통계는 **설계 판단이 아니다.** 모호 염기(N)가 있으면 Tm만 비우고
        # 나머지 통계는 그대로 준다(calculate_tm은 그 경우 명시적으로 실패한다).
        'tm': (calculate_tm(sequence) if not (set(sequence.upper()) - set("ACGT")) else None),
        'base_composition': base_composition,
        'restriction_sites': find_restriction_sites(sequence),
        'has_repeats': len(find_repeats(sequence, min_length=20)) > 0
    }


def longest_common_substring_len(seq_a: str, seq_b: str) -> int:
    """두 서열의 **가장 긴 공통 연속 부분서열** 길이.

    homology 기반 조립(Gibson / yeast HR)에서 "의도한 접합부 말고도 서로 붙을 수 있는
    구간이 있는가"를 묻는 데 쓴다. 의도한 접합부는 인접 fragment끼리만 공유해야 하고,
    **비인접 fragment가 긴 구간을 공유하면 엉뚱한 자리끼리 재조합**한다.
    (실측: 같은 부품을 두 자리에 쓰면 비인접 fragment가 81bp를 공유했다)

    길이에 대한 이분탐색 + k-mer 집합으로 찾는다. 공통 부분서열의 존재는 길이에 대해
    단조(길이 k가 있으면 k-1도 있다)라서 이분탐색이 성립한다.

    Args:
        seq_a: 서열 A
        seq_b: 서열 B

    Returns:
        가장 긴 공통 연속 부분서열의 길이. 공통 구간이 없으면 0.

    Note:
        대소문자를 구분하지 않는다. 역상보는 보지 않는다 — 역상보 교차까지 봐야 하면
        호출자가 한쪽을 `reverse_complement()`해서 한 번 더 부른다.

    Example:
        >>> longest_common_substring_len("AAATTTGGG", "CCCTTTGAA")
        4
    """
    a = (seq_a or "").upper()
    b = (seq_b or "").upper()
    if not a or not b:
        return 0

    lo, hi, best = 1, min(len(a), len(b)), 0
    while lo <= hi:
        k = (lo + hi) // 2
        kmers_a = {a[i:i + k] for i in range(len(a) - k + 1)}
        if any(b[i:i + k] in kmers_a for i in range(len(b) - k + 1)):
            best, lo = k, k + 1
        else:
            hi = k - 1
    return best


# 모듈 메타데이터
__all__ = [
    'validate_dna_sequence',
    'validate_cds',
    'calculate_gc_content',
    'calculate_tm',
    'find_restriction_sites',
    'check_restriction_site_in_part_sequence',
    'find_repeats',
    'longest_common_substring_len',
    'get_sequence_stats'
]
