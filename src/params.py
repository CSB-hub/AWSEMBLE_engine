"""설계 파라미터 — 유저(연구실)가 입력하는 값이 모이는 **단일 지점**.

## 왜 이 모듈이 있나

Tm·GC·길이·homology·효소는 연구실마다 다르다(Q5를 쓰는 랩과 Taq을 쓰는 랩의 Tm 기준이 다르고,
Gibson homology 길이도 랩마다 다르다). 그런데 이 값들은 지금까지 엔진 5계층에 흩어져
하드코딩돼 있었고, 각 계층에 인자는 있는데 **아무도 값을 넘기지 않아** 전부 기본값으로만 돌았다.
즉 랩별로 다른 값을 넣을 방법 자체가 없었다.

`DesignParams`는 그 값들이 모이는 한 곳이다. **연구실 회신이 오면 이 파일의 기본값만 고친다.**
다른 파일을 뒤져야 한다면 이 설계는 실패한 것이다.

## 하드 리밋 vs 소프트 권장

- **하드**: 벗어나면 코드가 깨지거나 물리적으로 불가능한 값 → `validate()`가 **에러**를 낸다.
  실측으로 확인한 값만 넣는다(교과서적 이상치를 하드로 두지 않는다).
- **소프트**: 벗어나도 동작하는 값 → **경고만** 반환한다.
  "툴은 max·min만 제한한다"는 방침의 실체가 이 구분이다.

⚠️ 소프트 경고의 경계값은 **연구실 회신 전까지의 잠정치**다(참고용 일반 범위).
회신이 오면 이 파일의 `_SOFT_*` 상수만 갱신한다.

## 기본값

**현재 하드코딩돼 있던 값 그대로**다. 아무것도 넘기지 않으면 기존과 100% 동일하게 동작한다.
"""
from dataclasses import dataclass, fields, replace
from typing import Any, Dict, List, Optional

# ─────────────────────────────────────────────────────────────────────────────
# 하드 리밋 — 실측으로 확인한 값 (2026-09-03, primer3 2.2.0)
# ─────────────────────────────────────────────────────────────────────────────

# primer3 내장 최대치. 초과 시 OSError: "PRIMER_MAX_SIZE exceeds built-in maximum of 36"
PRIMER3_MAX_SIZE = 36

# primer3 열역학 계산(calc_hairpin/calc_heterodimer) 입력 상한.
# 초과 시 RuntimeError: "At least one sequence must be equal to or shorter than 60bp ..."
# 조립용 프라이머는 [5' tail][anneal] 구조라 tail+anneal이 이 값을 넘으면 품질 계산이 불가능하다.
# ⚠️ 이건 **primer3의 계산 상한**이고 합성 상한이 아니다 — 63nt 프라이머는 주문·실험이 되지만
#    우리가 hairpin/dimer를 측정하지 못한다. 이 구간은 `basic.secondary_structure`의
#    하이브리드 경로(>60bp → ViennaRNA-DNA, 경고 전용)가 담당한다.
THERMO_MAX_LEN = 60

# 올리고 **합성** 상한 (nt). `oligo_max_bp`(어느 길이에서 올리고를 권하는가)와 다른 축이다.
# 일반 올리고는 60nt급이고 그 위는 ultramer 등급이며, 이 값을 넘으면 유전자 합성으로 넘어간다.
# junction 설계에서 부품에 homology arm을 붙이면 주문 길이가 늘어나므로(실측: arm 20bp에서
# rbs 92%·terminator 77%가 100nt 이하, promoter 짧은 것 37%가 100~200nt), 병합이 이 상한을
# 넘기면 병합하지 않고 분리한다.
# ⚠️ `DesignParams` 필드가 아니라 모듈 상수다 — 연구실이 run마다 바꾸는 값이 아니고,
#    필드로 올리면 `canonical()`에 들어가 기존 산출물의 artifact id가 전부 이동한다.
OLIGO_SYNTH_MAX_BP = 200

# 주문할 **전체 프라이머**(5' tail + anneal)의 길이 상한(nt). 연구실 주문 기준(2026-10-25).
# ⚠️ `THERMO_MAX_LEN`(60)과 다른 축이다 — 저건 primer3가 **계산**할 수 있는 한계고,
#    이건 **주문**할 수 있는 한계다. 60~80nt 구간은 주문은 되지만 hairpin/dimer가
#    경고 전용으로만 판정된다(`basic.secondary_structure`의 하이브리드 경로).
# 짧은 부품을 이웃 프라이머 tail에 실을지(tail 흡수) 올리고로 뺄지를 이 값이 가른다.
PRIMER_ORDER_MAX_NT = 80

# ─────────────────────────────────────────────────────────────────────────────
# 소프트 권장 경계 — 회신 전 잠정치 (CP9에서 실제 값으로 교체)
# ─────────────────────────────────────────────────────────────────────────────
_SOFT_TM_RANGE = (45.0, 75.0)     # 선언 범위가 이 밖이면 경고
_SOFT_TM_WIDTH = (4.0, 15.0)      # 선언 범위의 폭(max-min) 권장 구간
_SOFT_PRIMER_LEN_MIN = 15
# GC 설계창의 소프트 경계. **기본값(20~80%)에서 경고가 뜨면 안 된다** — 2026-11-05까지
# (30, 70)이어서 아무 설정도 바꾸지 않은 실행에서도 "범위가 넓습니다" 경고가 떴다.
# 기본값에서 뜨는 경고는 정의상 정보가 아니라 소음이다. 설계창이 넓은 건 의도된 것이고
# (primer3가 AT-rich 영역에서 후보를 못 찾는 걸 막는다), 품질은 Tm·판정창이 맡는다.
_SOFT_GC_RANGE = (20.0, 80.0)
_SOFT_GIBSON_OVERLAP_RANGE = (15, 40)

# ─────────────────────────────────────────────────────────────────────────────
# Tm 선언값의 **하드 입력 제한** (2026-11-05 Y11)
# ─────────────────────────────────────────────────────────────────────────────
# 예전 구조는 유저 설계창(`tm_min/tm_max`) **바깥에 판정창을 하나 더** 뒀다
# (`qc_tm_min/max` 45~85℃, `qc_gc_ideal_*` 40~60%). 그 구조가 2026-11-04부터 같은
# 버그를 세 번 냈다 — Y6(Tm) · Y7(GC) · Y11(길이). 두 창이 어긋나면 ① 설계 목표대로
# 만든 프라이머를 판정이 떨어뜨리거나 ② 유저가 좁힌 창을 판정이 못 보고 넘겼다.
# 길이 축은 아예 판정 쪽이 `15`/`30` **리터럴**이어서 유저가 Q5(상한 26nt)를 골라도
# 경고는 30 기준으로 떴다 — 유저 설정과 경고가 따로 노는 상태였다.
#
# → 창을 **하나로** 합친다. 판정도 유저가 선언한 값(`tm_*` `gc_*` `primer_len_*`)으로
#   한다. 판정창이 맡던 "물리적으로 PCR이 성립하지 않는 구간 차단"은 **입력 단계**로
#   옮긴다 — 애초에 못 넣게 하면 판정에서 걸러낼 필요가 없고, 창이 하나뿐이면
#   어긋남이 구조적으로 불가능하다.
#
# 경계 근거: 35℃ 아래는 annealing 단계에서 프라이머가 주형에 붙어 있지 못하고,
# 80℃ 위는 extension 온도(68~72℃)를 넘어 annealing 단계 자체가 뜻을 잃는다.
# ⚠️ 두 값 모두 wet-lab 확정 대상이다.
TM_INPUT_MIN = 35.0
TM_INPUT_MAX = 80.0

# yeast HR homology arm의 **하드** 범위. Gibson(소프트 권장)과 달리 하드로 둔 이유는
# 범위를 벗어나면 그 방법이 성립하지 않기 때문이다 — 너무 짧으면 세포의 재조합 기구가
# 상동성을 인식하지 못하고, 너무 길면 fragment가 이웃 서열로 뒤덮여 조립 단위가 흐려진다.
# ⚠️ 경계값은 문헌 기반 잠정치다 — wet-lab 확정 대상.
_YEAST_HOMOLOGY_RANGE = (30, 60)

ASSEMBLY_TYPES = ("linear", "circular")

# circular(중간 클로닝)에서 module을 중간 벡터로부터 꺼내는 방법.
#   digest : 다른 Type IIS 효소로 잘라낸다. Type IIS 조립은 인식서열을 없애므로
#            **조립 효소로는 되꺼낼 수 없고** 중간 벡터가 제2 효소 사이트를 제공해야 한다.
#            실측: 실제 module 120개 중 90%에서 가능(부품 내부에 없는 효소가 존재).
#   pcr    : spacer에 프라이머를 앉혀 증폭한다. 부품 내부 사이트 제약이 없다.
# 기본은 digest이고, 부품 내부에 쓸 수 있는 효소가 없으면 pcr로 자동 폴백한다.
EXCISION_METHODS = ("digest", "pcr")

# ── 부품 타입 (단일 정의) ────────────────────────────────────────────────────
# 조립 슬롯에 들어갈 수 있는 부품 타입. **튜플 순서가 곧 기본 배치 순서**다.
#
# `custom_part`: 유저가 **미리 융합해 둔 조합 파츠**(promoter-RBS 세트, terminator 2연 등).
#   조합의 경우의 수가 너무 많아 타입을 일일이 만들 수 없으므로 범용 타입 하나로 받는다.
#   쓰는 이유는 **내부 경계가 Golden Gate로 바뀌지 않는 것**이다 — 엔진은 부품 1개를
#   fragment 1개로 보고 경계에만 overhang을 만들므로(`prepare_fragments`), 융합 서열
#   내부는 건드리지 않는다. 즉 타입만 추가하면 요구가 충족된다.
ASSEMBLY_PART_TYPES = ("promoter", "rbs", "cds", "terminator", "custom_part")

# 공용 Bioparts(MySQL)에 **테이블이 있는** 부품. 서열로 부품을 역검색할 때의 범위다.
# `custom_part`가 여기 없는 건 의도다 — 유저가 각자 만든 조합이라 유저별 sandbox에만
# 두고 공용 DB에는 테이블을 만들지 않는다(개인 계정에 귀속).
# `spacer`는 조립 슬롯이 아니지만 서열 역검색 대상이라 포함된다.
MYSQL_PART_TABLES = ("promoter", "rbs", "cds", "terminator", "spacer")

# 공용 DB 없이 유저 sandbox에만 존재하는 타입.
SANDBOX_ONLY_PART_TYPES = tuple(t for t in ASSEMBLY_PART_TYPES if t not in MYSQL_PART_TABLES)

# 지원하는 Type IIS 제한효소.
# 각 효소의 **overhang 길이와 spacer 길이는 여기 적지 않는다** — Biopython의 절단부위 정의
# (elucidate)에서 파생시킨다(`basic.sequence_analysis.get_typeIIS_params`).
#   BsaI  GGTCTCN^NNNN_N   → spacer 1, overhang 4
#   BbsI  GAAGACNN^NNNN_N  → spacer 2, overhang 4
#   BsmBI CGTCTCN^NNNN_N   → spacer 1, overhang 4
#   SapI  GCTCTTCN^NNN_N   → spacer 1, overhang **3**  ← 4 고정 가정이 깨지는 케이스
# 효소를 늘리려면 이 튜플에 이름만 추가하면 된다(값은 Biopython이 준다).
SUPPORTED_ENZYMES = ("BsaI", "BbsI", "BsmBI")

# ── 아직 못 쓰는 효소: 코드가 아니라 **데이터**가 없어서 막혀 있다 ──────────────
# SapI는 3nt overhang이 필요한데 DB `overhang` 테이블(institution='Awsemble')은
# 185개가 전부 4bp다. 4bp 세트를 3bp로 자르면 "서로 오결합하지 않도록 검증된 조합"이라는
# 성질이 깨지므로 잘라 쓰면 안 된다(2026-09-04 판단).
#
# 추가 방법: 해당 길이의 fidelity 검증 세트를 DB에 등록한 뒤
#            SUPPORTED_ENZYMES에 이름만 넣으면 된다. 코드 변경은 필요 없다
#            (overhang·spacer 길이는 Biopython에서 파생되고,
#             `_get_high_fidelity_overhang_set`이 길이가 맞는 후보만 고른다).
ENZYMES_PENDING_OVERHANG_SET = {"SapI": 3}   # {효소명: 필요한 overhang 길이(nt)}

# ── 인식서열 5' 앞에 붙이는 여분 염기 ────────────────────────────────────────
# 제한효소는 인식서열이 DNA 단편 **말단에 붙어 있으면 절단 효율이 크게 떨어진다**
# (NEB "Cleavage Close to the End of DNA Fragments"). 구버전 프라이머는 인식서열이
# 5' 맨 끝에서 시작해 여분 염기가 0bp였다 → PCR 산물을 자를 때 효율 저하.
#
# ⚠️ 아래 값은 **잠정치(6bp)** 다. NEB 표의 효소별 수치를 받으면 이 표만 교체하면 된다
#    (2026-09-04, 공동연구자 피드백 반영. 실제 수치 확인 전).
SITE_PAD_DEFAULT = {"BsaI": 6, "BbsI": 6, "BsmBI": 6}
SITE_PAD_FALLBACK = 6

# 패딩에 쓸 기본 서열. 아무 염기나 쓰면 안 된다 —
#   · 인접 인식서열과 합쳐져 **새 제한효소 사이트**가 생길 수 있고
#   · self-complementary면 hairpin이 되며
#   · spacer 더미가 'A'*n이라 패딩도 A면 poly-A 구간이 생긴다.
# 6-mer 4096개를 전수 검사해(GC 33~50% · 3연속 없음 · self-comp 아님 ·
# 지원효소 사이트 신규 생성 없음) 통과한 1921개 중 A/T 시작 후보를 골랐다.
SITE_PAD_SEQ_DEFAULT = "AATACC"


# ── polymerase 프리셋 ────────────────────────────────────────────────────────
# 프리셋이 정하는 것은 주로 **Tm의 '단위'**(na_mM·primer_conc_nM)다. 같은 서열도
# 염·농도가 다르면 Tm이 몇 도씩 달라지므로, polymerase를 바꾸면 이 두 값이 먼저 바뀐다.
#
# ⚠️ **담는 값의 기준**: 벤더 프로토콜로 방어 가능한 값만 담는다. 근거가 약한 값을
#    프리셋에 박으면 유저가 "랩 조건"으로 오인한다. 그래서 프리셋에 없는 필드는
#    전부 기본값을 그대로 상속한다(프리셋 = 부분 덮어쓰기).
#
# ⚠️ **Tm 범위는 잠정치다.** 벤더가 문서화하는 것은 Ta(annealing temperature)이고
#    우리가 받는 값은 설계 Tm 범위라 서로 같지 않다(Ta는 polymerase마다 Tm에
#    가감하는 규칙이 다르다). 우리는 Ta→Tm 변환을 하지 않기로 했으므로(2026-09-23 결정)
#    여기 값은 회신이 오면 교체할 자리표시다.
#
# 폭(max−min)은 **10℃로 통일**했다. 실측상 조건(염·농도)이 바뀌면 Tm 분포는
# **폭은 그대로인 채 통째로 이동**한다(2026-09-23: 농도 200→500nM에서 전 서열이 +0.9℃,
# 염 50→200mM에서 +7.2℃ 이동, 분포 폭 14.5℃ 불변). 그래서 polymerase가 바꾸는 것은
# 범위의 **위치**이지 폭이 아니다.
#
# 회신이 오면 **이 딕셔너리만 고친다** — params.py 단일 지점 원칙.
# ─────────────────────────────────────────────────────────────────────────────
# NEB Tm Calculator 재현 조건
#
# ## 왜 이 값들인가 — 증빙
#
# NEB Tm Calculator(https://tmcalculator.neb.com)는 **제품마다 다른 Tm**을 낸다.
# 버퍼 조성이 다르기 때문인데, Q5/Phusion 버퍼 조성은 **NEB 영업비밀로 비공개**다.
# 그래서 오랫동안 오프라인 재현이 안 되는 것으로 알려져 있었다.
#
# 해법은 pydna 이슈 #237(https://github.com/pydna-group/pydna/issues/237)에 있다.
# **Biopython `MeltingTemp`의 저자 본인(MarkusPiotrowski)**이 직접 맞는 인자를 찾아 올렸다.
# 핵심 통찰 두 가지:
#   · 버퍼 조성을 **몰라도 된다.** Na + K + Tris/2 가 하나의 1가 염 농도로 합산되므로,
#     비공개인 Q5/Phusion도 K/Tris 조합으로 같은 유효 농도를 만들면 값이 맞는다.
#   · **NEB는 Taq에 Mg2+ 보정을 하지 않는다.** 그래서 Taq만 Mg=0이다.
#     (Mg=0이면 Biopython은 dNTPs도 무시한다)
#     Phusion만 염보정이 Schildkraut(saltcorr=1)이고 나머지는 Owczarzy(6)인데,
#     이는 NEB 자신의 API 문서와도 일치한다 —
#     "For Phusion DNA Polymerases, the salt correction of Schildkraut is used."
#
# ## 자체 검증 (2026-11-05, NEB Tm API v0.19.0)
#
# 실제 Bioparts 파츠 말단 160개(길이 18~36nt)를 NEB API와 대조:
#
#     제품       RMSE     최대오차   ±1℃ 내   Ta 일치
#     Taq        0.29℃    0.50      100%     100%
#     Phusion    0.62℃    1.45       88%      95%
#     Q5         0.35℃    0.72      100%     100%
#
# ⚠️ **함정**: NEB API의 `conc` 파라미터는 **nM이 아니라 µM**이다. nM로 알고 500을
#    넣으면 1000배 과대평가라 Tm이 ~8℃ 높게 나온다(2026-11-05에 실제로 겪었다).
#    재검증할 때 반드시 µM으로 보낼 것 — Q5/Phusion은 0.5, Taq은 0.2.
# ─────────────────────────────────────────────────────────────────────────────

#: polymerase별 Tm 계산 조건. `DesignParams` 필드와 같은 이름을 쓴다.
NEB_TM_CONDITIONS: Dict[str, Dict[str, Any]] = {
    # Taq: 버퍼 조성이 **공개**돼 있다(10mM Tris + 50mM KCl + 1.5mM MgCl2).
    #      그런데 NEB 계산기는 Taq에 Mg 보정을 하지 않으므로 mg_mM=0으로 둔다.
    "taq":     {"na_mM": 0.0, "k_mM": 50.0, "tris_mM": 10.0, "mg_mM": 0.0,
                "dntp_mM": 0.8, "tm_saltcorr": 6, "primer_conc_nM": 200.0},
    # Phusion: HF/GC 버퍼 조성 비공개. Schildkraut 염보정(NEB 문서가 명시).
    "phusion": {"na_mM": 0.0, "k_mM": 50.0, "tris_mM": 25.0, "mg_mM": 1.5,
                "dntp_mM": 0.8, "tm_saltcorr": 1, "primer_conc_nM": 500.0},
    # Q5: 버퍼 조성 비공개. Owczarzy 염보정 + Mg 보정 포함.
    "q5":      {"na_mM": 0.0, "k_mM": 50.0, "tris_mM": 10.0, "mg_mM": 1.5,
                "dntp_mM": 0.8, "tm_saltcorr": 6, "primer_conc_nM": 500.0},
}

#: polymerase별 권장 어닐링 온도(Ta) 규칙. (종류, 오프셋, 상한℃)
#  NEB 계산기가 쓰는 규칙이며 위 검증에서 Ta 일치율 95~100%로 확인했다.
#  Ta는 **두 프라이머 중 낮은 Tm**을 기준으로 한다.
NEB_TA_RULES: Dict[str, tuple] = {
    "taq":     ("offset", -5.0, 68.0),   # min(Tm1,Tm2) - 5, 최대 68℃
    "phusion": ("phusion", 0.0, 72.0),   # 0.93 * min(Tm1,Tm2) + 7.5, 최대 72℃
    "q5":      ("offset", +1.0, 72.0),   # min(Tm1,Tm2) + 1, 최대 72℃
}


PRESETS: Dict[str, Dict[str, Any]] = {
    # ⚠️ 버퍼 조건은 **`NEB_TM_CONDITIONS`에서 끌어온다** — 같은 값을 두 곳에 적지 않는다.
    #    (2026-11-05 이전에는 셋 다 na_mM=50 '잠정'이라 효소를 바꿔도 Tm이 거의 안 움직였다)
    #    여기 남는 건 **Tm 선언창**뿐이고, 그건 여전히 wet-lab 확정 대상이다.
    # 2026-11-05 Y10-b — Tm 창을 NEB 축에 맞게 재설정했다.
    #   이전 창(50~60 / 55~65 / 57~67)은 **primer3 축 시절 숫자**였다. 축을 NEB로 옮기면
    #   같은 프라이머의 Tm이 제품별로 8.86℃까지 벌어지는데, 창은 그대로 둬서 "Taq은
    #   창이 너무 뜨겁고 Q5는 너무 차가운" 상태였다.
    #   NEB 축에서 제품별로 실제 도달 가능한 구간 + NEB가 권고하는 Ta 대역을 기준으로 잡았다.
    #     · Taq      55~65  (Ta = min−5 → 50~60℃, Taq의 일반적 Ta 대역)
    #     · Phusion  60~72  (Ta = 0.93·min+7.5 → 63.3~74.5, 상한 72로 캡)
    #     · Q5       60~72  (Ta = min+1 → 61~73, 상한 72로 캡)
    #   Phusion·Q5가 같은 창인 건 우연이 아니다 — 둘 다 high-fidelity 고온 효소라
    #   NEB 권고 Ta 대역이 겹친다. 갈리는 건 Tm 계산 쪽이다(saltcorr·Tris 조성).
    #
    # **길이창도 프리셋이 정한다** — Tm 창과 커플링돼 있기 때문이다.
    #   제품마다 Tm 자가 다르므로(같은 프라이머가 Taq 54.5 / Q5 62.3℃) 같은 Tm에 닿는 데
    #   필요한 결합부 길이가 다르다. 길이창을 전 제품 공통(18/20/25)으로 두면 Tm 창만
    #   제품별로 갈려서, **Taq·Phusion은 선언한 창에 닿지 못해 상시 경고**가 떴다.
    #
    #   특히 병목은 `primer_len_max`가 아니라 **`primer_len_opt`**였다. primer3의 목적함수가
    #   |길이 − OPT_SIZE|에 벌점을 주므로, 상한만 30으로 올려도 반환되는 후보 30쌍이
    #   20nt 근처에 몰려 26~28nt 후보가 상위에 안 올라온다. 실측(4부품 module, AT-rich
    #   접합부 포함):
    #     Taq  opt20/max25 → 경고 1/2 · opt20/max30 → 1/2 · **opt24/max30 → 0/2**
    #   즉 상한을 36까지 올려도 opt를 안 올리면 해결되지 않았다.
    #
    #   값은 "경고 0이 되는 최소"로 잡았다 — 길이는 합성 단가이므로 넉넉하게 둘 이유가 없다.
    "taq":     {**NEB_TM_CONDITIONS["taq"],
                "tm_min": 55.0, "tm_max": 65.0,       # ⚠️ 창은 잠정 (wet-lab 확정 대상)
                "primer_len_opt": 24, "primer_len_max": 30},
    "phusion": {**NEB_TM_CONDITIONS["phusion"],
                "tm_min": 60.0, "tm_max": 72.0,       # ⚠️ 창은 잠정
                "primer_len_opt": 22, "primer_len_max": 30},
    "q5":      {**NEB_TM_CONDITIONS["q5"],
                "tm_min": 60.0, "tm_max": 72.0,       # ⚠️ 창은 잠정
                "primer_len_opt": 20, "primer_len_max": 26},
}


@dataclass(frozen=True)
class DesignParams:
    """설계 파라미터 한 벌. **불변(frozen)** — 변형은 `with_(...)`로 새 인스턴스를 만든다.

    불변으로 둔 이유: 이 객체가 `derived_id` 해시 재료로 쓰이므로, 파이프라인 중간에
    누가 값을 바꾸면 id와 실제 산출물이 어긋난다.
    """

    # ── PCR 조건 (Tm 계산의 전제) ────────────────────────────────────────────
    # ⚠️ 이 값들이 Tm의 '단위'를 정한다. 같은 서열도 염/농도가 다르면 Tm이 몇 도씩 달라진다.
    #    "Tm 70℃"라는 말이 랩마다 다른 프라이머를 가리키는 원인.
    #
    # 2026-11-05: **NEB Tm Calculator 재현 조건으로 교체**(`NEB_TM_CONDITIONS` 참고).
    #    polymerase를 고르면 이 값들이 따라오고, 그 결과 Tm이 실제로 달라진다.
    #    (그 전에는 셋 다 na_mM=50으로 눌려 있어 효소 선택이 Tm에 거의 영향이 없었다)
    na_mM: float = 0.0             # Na+ 농도 (mM) — NEB 버퍼는 K+/Tris 기반이라 0이다
    k_mM: float = 50.0             # K+ 농도 (mM)
    tris_mM: float = 10.0          # Tris 농도 (mM) — 1가 염 환산에 Tris/2가 더해진다
    mg_mM: float = 1.5             # Mg2+ 농도 (mM) — Taq 경로는 0(아래 설명)
    dntp_mM: float = 0.8           # dNTP 총 농도 (mM) — Mg2+를 킬레이트해 유효 Mg를 낮춘다
    tm_saltcorr: int = 6           # 염보정 방식 (Biopython MeltingTemp) — 6=Owczarzy, 1=Schildkraut
    # 기본값은 **Q5 조건**이다(위 버퍼 필드와 한 벌). 조립용 증폭은 proofreading
    # 효소를 쓰는 것이 표준이고(Taq 오류율의 ~1/280), 실사용도 Q5/Phusion에 몰린다.
    primer_conc_nM: float = 500.0  # 프라이머(올리고) 농도 (nM) — Q5 프로토콜 0.5 µM

    # ── 프라이머 설계 목표 ───────────────────────────────────────────────────
    # ── 프라이머 Tm 범위 — **설계 목표이자 판정 기준** (2026-09-23 통합) ──────
    # 구버전은 `target_tm ∓ tm_tolerance`(설계창)와 `qc_tm_min/max`(판정창)를 **따로** 뒀다.
    # 그 둘이 어긋나 있어서(기본값 [45,75] vs [50,65]) **설계 목표대로 만든 프라이머가
    # 판정에서 떨어지는** 구간이 있었고, 유저가 목표 Tm을 올릴수록 그 구간이 커졌다
    # (target 70℃면 설계는 70을 겨냥하는데 판정은 65에서 잘랐다).
    #
    # 이제 **범위 하나**가 둘 다를 맡는다. 유저는 자기 조건에서 쓰는 Tm 범위를 선언하고,
    # 설계는 그 안을 겨냥하며 판정도 같은 범위로 한다 → 두 창이 어긋날 수가 없다.
    #
    # ⚠️ 유저 입력이 판정을 겸해도 되는 이유: 2026-07-13의 "전부 red" 실패는
    #    **설계 목표와 판정 기준이 다른 값이었기 때문**에 생겼다. 같은 범위를 쓰면
    #    "내가 겨냥한 곳에 만들어 놓고 내가 떨어뜨리는" 모순 자체가 없다.
    #    2026-11-05(Y11): GC·길이 축도 같은 처방으로 합쳤다. 판정창(`qc_tm_*`,
    #    `qc_gc_*`, `qc_gc_ideal_*`)은 제거하고 물리 경계는 입력 제한
    #    (`TM_INPUT_MIN/MAX`)으로 옮겼다 — 위 상수 블록 참고.
    tm_min: float = 55.0
    tm_max: float = 65.0
    primer_len_min: int = 18
    primer_len_opt: int = 20
    primer_len_max: int = 25
    gc_min: float = 20.0
    gc_max: float = 80.0
    max_poly_x: int = 8            # 같은 염기 연속 허용 개수
    # ⚠️ **primer3 축** (2026-11-05 통일). 이 값은 `_fit_amplicon`이 `calculate_tm`
    #    결과와 직접 비교하는 절대 Tm이라 Biopython 축에 묶여 있었다. 자를 primer3로
    #    통일하면서 뜻을 보존하려고 +4℃ 옮겼다(실측 오프셋 +4.27, 자세한 근거는
    #    `basic/sequence_analysis.calculate_tm` 독스트링).
    terminal_tm_min: float = 59.0  # module 끝단(tail 없는) 프라이머를 이 Tm까지 연장

    # ── Golden Gate (부품 → 모듈) ────────────────────────────────────────────
    # overhang 길이와 spacer 길이는 여기 없다 — **효소가 정하는 값**이라 자유 파라미터가 아니다.
    # `enzyme_params()`로 Biopython에서 파생시킨다. 효소만 고르면 나머지는 따라온다.
    gg_enzyme: str = "BsaI"
    # 인식서열 5' 앞 여분 염기. None이면 효소별 기본값(SITE_PAD_DEFAULT), 숫자면 그 값을 강제.
    site_5prime_pad: Optional[int] = None
    site_5prime_pad_seq: str = SITE_PAD_SEQ_DEFAULT

    # ── Gibson homology arm ──────────────────────────────────────────────────
    # module → vector 삽입과 **part → module 접합부가 이 값을 공유한다.**
    # 따로 두지 않은 이유: 한 플라스미드 안에서 접합부 길이를 섞을 이유가 없고,
    # 필드를 늘리면 `canonical()`에 들어가 기존 산출물의 artifact id가 전부 이동한다.
    # 2026-10-24: 기본값 20 → 30 (연구실 지정). 소프트 권장범위(15~40) 안이고
    # `arm + primer_len_min = 48bp`로 primer3 계산 상한(60) 안쪽이라 품질 게이트가 유지된다.
    # ⚠️ 이 값은 **최종 산물 서열에 들어간다**(벡터 삽입 homology) → baseline SHA가 이동했다.
    #    구 baseline(arm 20): module 601809584fd54212 / construct d69f339665bad36a
    gibson_overlap: int = 30       # homology arm 길이 — 범위의 **최소**이자 기본
    # 접합부 품질이 최소 길이에서 떨어지면 여기까지 늘려본다(통과하는 가장 짧은 길이 채택).
    # ⚠️ padding을 없앤 뒤로 **이 범위가 유일한 탈출구**다. 단일 길이로 고정하면
    #    실측상 모듈의 2.3%가 어떤 수로도 통과하지 못해 설계가 막힌다.
    #    길이는 공유 구간일 뿐 삽입이 아니라서 **최종 서열을 바꾸지 않는다**(주문서만 바뀜).
    gibson_overlap_max: int = 45

    # ── Yeast HR homology arm ────────────────────────────────────────────────
    # yeast에서 세포 안(in vivo) 재조합으로 잇는 homology 길이. **Gibson과 다른 축이다** —
    # Gibson은 시험관에서 exonuclease가 깎아낸 단일가닥이 어닐링하면 되므로 15~25bp로 충분하고,
    # yeast HR은 세포의 재조합 기구가 상동성을 인식해야 해서 문헌상 수십 bp 규모를 쓴다.
    # 그래서 `gibson_overlap`을 재사용하지 않고 별도 필드로 둔다.
    #
    # ⚠️ 상한은 **`oligo_max_bp`의 절반**에 묶인다 — PCR fragment가 자기 몸에서 좌우 접합부
    #    구간을 떼어주므로 `body >= 2*arm`이어야 두 구간이 겹치지 않는다(§junction_plan).
    #    기본 `oligo_max_bp=100`이면 실사용 상한이 50bp다. 더 긴 arm을 쓰려면 `oligo_max_bp`도
    #    함께 올려야 한다.
    #  min == max면 그 길이로 고정, 다르면 **품질을 통과하는 가장 짧은 길이**를 접합부마다
    #  따로 고른다. 길이는 공유 구간일 뿐 삽입이 아니라서 **최종 서열을 바꾸지 않는다**
    #  (바뀌는 건 주문서뿐) — 그래서 탐색이 안전하다.
    yeast_homology_min: int = 40
    yeast_homology_max: int = 60


    # ── module 양 끝 spacer를 어떻게 만들 것인가 ──────────────────────────────
    #  "vector_homology" (기본) — 선형화 자리의 양옆 서열을 그대로 양 끝에 붙인다.
    #      · module이 homology를 이미 들고 있어 **삽입 PCR이 필요 없다**
    #      · Gibson에서 homology는 벡터 서열과 병합되므로 **흉터가 남지 않는다**
    #        (library 방식은 중립 spacer 57bp × 2 = 114bp가 최종 플라스미드에 그대로 남는다)
    #      · spacer가 (vector, 효소)로 결정되므로 조합 탐색이 사라진다(구 최대 380회)
    #  "library" — 기존 방식. DB의 중립 spacer 풀에서 L/R을 골라 붙인다.
    #      랩에 이미 만들어 둔 spacer 실물이 있을 때 쓴다.
    spacer_mode: str = "vector_homology"

    # Gibson/Yeast HR에서 vector를 열 제한효소. None이면 LINEARIZATION_ENZYMES 순서로 자동.
    # **정확히 한 곳만 자르는 효소**여야 한다 — 유저가 고른 값이 그렇지 않으면 몇 곳에서
    # 자르는지 알려주고 중단한다(조용히 다른 효소로 바꾸지 않는다).
    linearization_enzyme: Optional[str] = None

    # ── 프라이머 길이 수렴 (설계 목표) ─────────────────────────────────────────
    # ⚠️ tm_diff_max는 **GG 프라이머 길이를 실제로 정하는 knob**이다.
    #    anneal(주형 결합부) Tm과 전체 프라이머(5' tail 포함) Tm의 차이를 이 값 이하로
    #    수렴시키려고 amplicon을 늘린다. 값을 키우면 프라이머가 짧아지고 작아지면 길어진다.
    tm_diff_max: float = 5.0
    # ⚠️ 위 `tm_diff_max`와 **다른 축**이다. 이쪽은 fwd와 rev **두 프라이머 사이**의
    #    Tm 차이다(같은 프라이머의 anneal↔full이 아니다).
    #    왜 필요한가: annealing 온도는 `ta_recommended()`가 **낮은 쪽 Tm 기준**으로
    #    정한다. 두 Tm이 벌어지면 그 Ta에서 **뜨거운 쪽 프라이머가 자기 Tm보다 훨씬
    #    낮은 온도에서 돌아** 비특이 결합이 늘어난다. 실측(2026-11-05, 조성 3 × 경로 2
    #    × 프리셋 3 = 69쌍): 중앙 0.89℃인데 최대 **11.68℃**였고 5℃ 초과가 17%였다.
    #    최악 케이스는 GG 벡터 homology 프라이머에 몰렸다(길이가 overhang·spacer에
    #    묶여 Tm을 맞출 자유도가 적다).
    #    ⚠️ 경계값은 통상 권장치(≤5℃) 기반 잠정치다 — wet-lab 확정 대상.
    pair_tm_diff_max: float = 5.0
    # anneal 최소 길이. 이보다 짧으면 비특이 결합 위험이 커진다.
    # ⚠️ 고GC 구간에서 Tm을 낮출 수 있는 **유일한 레버의 바닥**이기도 하다(GG 삽입 경로).
    amp_min: int = 15

    # ── 부품 제작 방식 경계 (올리고 어닐링 vs PCR) ────────────────────────────
    # 짧은 부품은 PCR로 증폭할 게 아니라 **올리고 2개를 어닐링**해서 만드는 것이 표준이다.
    # 표준 올리고 합성이 100 base까지 되므로 그 아래는 올리고가 비용·확실성 모두 낫다.
    # 실측(2026-10-12, Cello 문헌 파츠)에서 PCR 경로의 문제가 전부 이 구간에 몰려 있었다:
    #   9bp → Tm 33.5~36.2℃, 25bp → 49.7~50.8℃ (전 프리셋 창 밖인데 passed=True)
    oligo_max_bp: int = 100        # 이 미만이면 올리고 주문시트를 기본으로 제시

    # ── 어닐링 듀플렉스 안정성 등급 ──────────────────────────────────────────
    # one-pot Golden Gate는 16~37℃를 오간다. 어닐링한 두 가닥을 붙잡는 건 **부품 본체**
    # 길이뿐이다(overhang은 단일가닥이고, 패딩·인식서열은 조립 전에 잘려나간다).
    # 실측 듀플렉스 Tm (GC 40~60%):
    #     10bp → 22~32℃   37℃에서 풀림
    #     15bp → 41~50℃   여유 4~13℃
    #     20bp → 52~60℃   여유 15℃+
    # 차단하지는 않는다 — 시뮬레이션은 되고, 판단은 실험자 몫이다. 대신 **등급을 숨기지 않는다.**
    duplex_safe_bp: int = 20       # 이상이면 안전
    duplex_risky_bp: int = 15      # 이상~safe 미만이면 주의 / 미만이면 위험

    # ── 품질 판정 축 (tier) ───────────────────────────────────────────────────
    # ⚠️ **위의 '설계 목표'와 의도적으로 분리된 값이다.**
    #    설계 목표(gc_min/max 등)는 "어떤 프라이머를 만들까"를 정하고,
    #    여기 qc_* 값은 "만들어진 프라이머가 쓸 만한가"를 판정한다.
    #    둘을 합치면 유저가 판정 기준을 낮춰 통과시키는 길이 열리므로 합치지 않는다.
    #    → 그래서 **GUI에 노출하지 않는다.** 코드 내부에서 경로별로 다른 값을 주기 위한 필드다.
    #
    # ⚠️ 아래 값은 전부 **wet-lab 확정 전 잠정치**다(7/30 요청 문서 대상).
    #    ΔG 컷오프는 커뮤니티 관례(37℃ 기준)이지 단일 권위 표준이 아니다.
    # ── Tm 판정창 — 설계창(tm_min/tm_max)과 **분리**, 단 반드시 설계창을 포함한다 ──
    # 2026-09-23에 둘을 하나로 합친 건 "설계 목표와 판정 기준이 다른 값"이라 설계대로
    # 만든 프라이머가 판정에서 떨어지던 모순(2026-07-13 '전부 red') 때문이었다.
    # 그런데 합친 결과 **설계창이 판정창까지 겸하게 되어** 부작용이 생겼다:
    #   · 설계창을 넓히면 선정 규칙이 "창 안에서 가장 짧은 것"을 집으므로(primer_design
    #     `_fit_amplicon`) 58℃ 20-mer를 만들 수 있는 자리에서 50℃ 15-mer가 선택된다
    #     → 판정을 느슨하게 하려고 창을 넓히면 **프라이머가 실제로 나빠진다.**
    #   · 반대로 좁게 두면, tail 없는 module 끝단 프라이머(anneal==full)의 Tm 미달이
    #     violation(red)으로 올라간다. tail 있는 쪽은 같은 미달이 warning(yellow)이라
    #     (`primer_design` anneal_gate) **같은 접합부가 arm 운반자에 따라 적/황으로 갈렸다.**
    #
    # 그래서 다시 분리하되, **판정창이 설계창을 완전히 포함**하게 둔다. 그러면 "겨냥한
    # 곳에 만들어 놓고 내가 떨어뜨리는" 모순은 구조적으로 불가능하다 — 2026-07-13의
    # 실패는 판정창이 설계창보다 **좁거나 어긋나서** 생긴 것이고, 포함 관계면 그 방향이
    # 막힌다. `validate()`가 포함 관계를 검사해 어긋나면 경고한다.
    # ⚠️ Tm·GC 판정창(`qc_tm_*` / `qc_gc_*` / `qc_gc_ideal_*`)은 **2026-11-05(Y11)에
    #    제거**했다. 판정은 유저 선언창(`tm_min/tm_max`, `gc_min/gc_max`)이 직접 맡고,
    #    물리 경계는 입력 제한(`TM_INPUT_MIN/MAX`)이 맡는다 — 파일 상단 상수 블록 참고.
    #    아래 ΔG 컷오프는 **범위가 아니라 물리량**이라 그대로 남는다: 유저가 선언하는
    #    값이 아니고(설정 화면에 없다), hairpin/dimer는 ΔG 하나로 PCR이 안 걸리는
    #    신호라 "선언창을 벗어났나"와 무관하게 알려야 한다.
    qc_hairpin_dg_3prime: float = -2.0    # 3'말단 관여 hairpin (kcal/mol, 더 음수면 문제)
    qc_hairpin_dg_internal: float = -3.0  # 내부 hairpin
    qc_dimer_dg: float = -5.0             # self/cross dimer
    qc_dg_absolute: float = -9.0          # 위치·형태 불문 최상위 경고 tier

    # ── 조립 형태 ────────────────────────────────────────────────────────────
    # part → module 조립을 linear로 할지, 벡터를 끼워 circular(중간 클로닝)로 할지.
    module_assembly_type: str = "linear"

    # ── 중간 클로닝 (circular일 때만) ────────────────────────────────────────
    # 중간 벡터가 spacer를 이미 갖고 있는가.
    #   False(기본): 벡터에 spacer가 없고 **부품 양 끝에 spacer를 붙여** 조립한다.
    #                → 벡터 하나로 모든 module 위치를 커버할 수 있다(V0001 같은 기성 GG 벡터 사용 가능).
    #                  꺼낼 때는 spacer에 프라이머를 앉혀야 양 끝이 spacer가 된다(= module별 구분되는 homology).
    #   True: 벡터에 spacer가 박혀 있다(전용 벡터). 위치마다 다른 벡터가 필요하다.
    intermediate_vector_has_spacers: bool = False

    excision_method: str = "digest"
    # module을 꺼낼 효소. None이면 부품 내부에 없는 효소를 자동 선택하고,
    # 그런 효소가 없으면 pcr로 폴백한다.
    excision_enzyme: Optional[str] = None

    # ------------------------------------------------------------------ #
    # 검증
    # ------------------------------------------------------------------ #
    def validate(self) -> List[str]:
        """하드 리밋 위반은 에러, 소프트 이탈은 경고 리스트로 반환.

        Returns:
            경고 문자열 리스트 (비었으면 권장 범위 안).

        Raises:
            ValueError: 하드 리밋 위반 — 그대로 두면 코드가 깨지거나 물리적으로 불가능한 값.
        """
        errs: List[str] = []

        # ── 하드: 프라이머 길이 ──
        if self.primer_len_max > PRIMER3_MAX_SIZE:
            errs.append(
                f"primer_len_max={self.primer_len_max} exceeds the primer3 built-in maximum {PRIMER3_MAX_SIZE} "
                f"(primer3 fails with 'PRIMER_MAX_SIZE exceeds built-in maximum')."
            )
        if self.primer_len_min < 1:
            errs.append(f"primer_len_min={self.primer_len_min}: must be 1 or more.")
        if self.primer_len_min > self.primer_len_max:
            errs.append(
                f"primer_len_min({self.primer_len_min}) > primer_len_max({self.primer_len_max})."
            )
        if not (self.primer_len_min <= self.primer_len_opt <= self.primer_len_max):
            errs.append(
                f"primer_len_opt({self.primer_len_opt}) is outside min-max"
                f"({self.primer_len_min}-{self.primer_len_max})."
            )

        # ── 하드: Tm / 농도 ──
        if self.tm_min <= 0 or self.tm_max <= 0:
            errs.append(f"tm_min/tm_max({self.tm_min}/{self.tm_max}): must be greater than 0.")
        elif self.tm_min >= self.tm_max:
            errs.append(
                f"tm_min({self.tm_min}) >= tm_max({self.tm_max}): the Tm window must be "
                f"[low, high] and wider than 0."
            )
        # 선언창이 곧 판정창이므로(Y11) 물리 경계는 **여기서** 막는다.
        # 판정에서 걸러내는 게 아니라 애초에 못 넣게 한다 — GUI·Galaxy XML의 입력
        # min/max와 같은 숫자이고, 이 검사가 GUI를 거치지 않는 호출의 방어선이다.
        if self.tm_min < TM_INPUT_MIN or self.tm_max > TM_INPUT_MAX:
            errs.append(
                f"The Tm window [{self.tm_min:g}, {self.tm_max:g}]C falls outside the allowed "
                f"input range [{TM_INPUT_MIN:g}, {TM_INPUT_MAX:g}]C - below {TM_INPUT_MIN:g}C "
                f"a primer cannot stay on the template during annealing, and above "
                f"{TM_INPUT_MAX:g}C it passes the extension temperature, so the annealing "
                f"step loses its meaning."
            )
        # ⚠️ 1가 염은 **합산**으로 본다(2026-11-05). NEB 버퍼는 Na+가 아니라 K+/Tris
        #    기반이라 `na_mM`만 보면 정상 설정(Na=0, K=50, Tris=10)을 오류로 막는다.
        #    Biopython MeltingTemp의 합산 규칙과 같다: Na + K + Tris/2.
        _monovalent = self.na_mM + self.k_mM + self.tris_mM / 2.0
        for _n, _v in (("na_mM", self.na_mM), ("k_mM", self.k_mM),
                       ("tris_mM", self.tris_mM), ("mg_mM", self.mg_mM),
                       ("dntp_mM", self.dntp_mM)):
            if _v < 0:
                errs.append(f"{_n}={_v}: cannot be negative.")
        if _monovalent <= 0:
            errs.append(
                f"Monovalent salt is 0 (Na {self.na_mM} + K {self.k_mM} + "
                f"Tris/2 {self.tris_mM / 2.0} = {_monovalent}) - Tm cannot be computed."
            )
        if self.tm_saltcorr not in (1, 2, 3, 4, 5, 6, 7):
            errs.append(
                f"tm_saltcorr={self.tm_saltcorr}: must be a Biopython MeltingTemp salt-correction "
                f"method number (1-7). NEB reproduction uses 6 (Owczarzy) or 1 (Schildkraut)."
            )
        if self.primer_conc_nM <= 0:
            errs.append(f"primer_conc_nM={self.primer_conc_nM}: must be greater than 0 (Tm cannot be computed).")

        # ── 하드: GC ──
        if not (0.0 <= self.gc_min <= 100.0) or not (0.0 <= self.gc_max <= 100.0):
            errs.append(f"gc_min/gc_max({self.gc_min}/{self.gc_max}) must be between 0 and 100.")
        elif self.gc_min > self.gc_max:
            errs.append(f"gc_min({self.gc_min}) > gc_max({self.gc_max}).")

        # ── 하드: Gibson homology ──
        if self.gibson_overlap < 1:
            errs.append(f"gibson_overlap={self.gibson_overlap}: must be 1 or more.")

        # part→module 접합부: PCR fragment는 자기 몸에서 좌우 접합부 구간을 떼어주므로
        # body >= 2*arm이어야 두 구간이 겹치지 않는다. PCR 그룹은 정의상 oligo_max_bp 이상이라
        # 그 전제를 여기서 본다 — 깨지면 양옆 fragment가 서로 붙는 경로가 열린다.
        if self.gibson_overlap_max < self.gibson_overlap:
            errs.append(
                f"gibson_overlap_max({self.gibson_overlap_max}) < "
                f"gibson_overlap({self.gibson_overlap}) - the range is inverted."
            )

        # ── 하드: Yeast HR homology ──
        # 생물학적 범위(문헌)와 구조적 상한(oligo_max_bp/2)을 따로 본다 — 둘이 서로
        # 다른 이유로 묶이므로 메시지도 분리해야 사용자가 무엇을 고쳐야 할지 안다.
        lo, hi = _YEAST_HOMOLOGY_RANGE
        if not (lo <= self.yeast_homology_min <= hi) or not (lo <= self.yeast_homology_max <= hi):
            errs.append(
                f"yeast_homology {self.yeast_homology_min}~{self.yeast_homology_max}: "
                f"must be within {lo}-{hi}bp (the yeast HR range reported in the literature)."
            )
        elif self.yeast_homology_min > self.yeast_homology_max:
            errs.append(
                f"yeast_homology_min({self.yeast_homology_min}) > "
                f"max({self.yeast_homology_max})."
            )

        # ── 하드: 기타 ──
        if self.gg_enzyme not in SUPPORTED_ENZYMES:
            errs.append(
                f"gg_enzyme={self.gg_enzyme!r}: unsupported enzyme "
                f"(supported: {list(SUPPORTED_ENZYMES)}). To add one, put its name in params.SUPPORTED_ENZYMES."
            )
        if self.site_5prime_pad is not None and self.site_5prime_pad < 0:
            errs.append(f"site_5prime_pad={self.site_5prime_pad}: cannot be negative.")
        if self.site_5prime_pad_seq and set(self.site_5prime_pad_seq.upper()) - set("ATGC"):
            errs.append(
                f"site_5prime_pad_seq={self.site_5prime_pad_seq!r}: only ATGC is allowed."
            )
        if self.max_poly_x < 1:
            errs.append(f"max_poly_x={self.max_poly_x}: must be 1 or more.")
        if self.excision_method not in EXCISION_METHODS:
            errs.append(
                f"excision_method={self.excision_method!r}: must be one of {EXCISION_METHODS}."
            )
        if self.excision_enzyme is not None and self.excision_enzyme not in SUPPORTED_ENZYMES:
            errs.append(
                f"excision_enzyme={self.excision_enzyme!r}: unsupported enzyme "
                f"(supported: {list(SUPPORTED_ENZYMES)})."
            )
        if self.excision_enzyme is not None and self.excision_enzyme == self.gg_enzyme:
            errs.append(
                f"excision_enzyme is the same as the assembly enzyme ({self.gg_enzyme}). Type IIS "
                f"assembly removes the recognition site, so the same enzyme cannot excise the module."
            )
        if self.module_assembly_type not in ASSEMBLY_TYPES:
            errs.append(
                f"module_assembly_type={self.module_assembly_type!r}: must be one of {ASSEMBLY_TYPES}."
            )

        # ── 하드: 프라이머 길이 수렴 / 판정 축 ──
        if self.tm_diff_max < 0:
            errs.append(f"tm_diff_max={self.tm_diff_max}: cannot be negative.")
        if self.pair_tm_diff_max < 0:
            errs.append(f"pair_tm_diff_max={self.pair_tm_diff_max}: cannot be negative.")
        if self.spacer_mode not in ("vector_homology", "library"):
            errs.append(
                f"spacer_mode={self.spacer_mode!r}: only 'vector_homology' or 'library' is allowed."
            )
        if self.linearization_enzyme:
            from basic.sequence_analysis import LINEARIZATION_ENZYMES
            if self.linearization_enzyme not in LINEARIZATION_ENZYMES:
                errs.append(
                    f"linearization_enzyme={self.linearization_enzyme!r}: not in the supported list "
                    f"({', '.join(LINEARIZATION_ENZYMES)})"
                )
        if self.amp_min < 1:
            errs.append(f"amp_min={self.amp_min}: must be 1 or more.")
        elif self.amp_min > THERMO_MAX_LEN:
            errs.append(
                f"amp_min={self.amp_min}bp exceeds the thermodynamic calculation limit {THERMO_MAX_LEN}bp "
                f"- primer quality cannot be computed at this length."
            )
        for _n, _v in (("qc_hairpin_dg_3prime", self.qc_hairpin_dg_3prime),
                       ("qc_hairpin_dg_internal", self.qc_hairpin_dg_internal),
                       ("qc_dimer_dg", self.qc_dimer_dg),
                       ("qc_dg_absolute", self.qc_dg_absolute)):
            if _v > 0:
                errs.append(
                    f"{_n}={_v}: a dG cutoff must be 0 or below "
                    f"(more negative means a more stable structure; anything below the cutoff is flagged)."
                )

        # ── 하드: GC 단위 실수 방어 (% vs 0~1 분수) ──
        # ⚠️ 판정 컷오프는 원래 코드에서 **분수**(0.35/0.60)였고 필드로 빼면서 %로 통일했다.
        #    옛 값을 그대로 넣으면 0.35%가 되어 **게이트가 사실상 꺼진 채로 조용히 통과**한다.
        #    0 초과 1 미만은 GC 백분율로 의미가 없으므로 단위 실수로 보고 차단한다.
        for _n, _v in (("gc_min", self.gc_min), ("gc_max", self.gc_max)):
            if 0.0 < _v < 1.0:
                errs.append(
                    f"{_n}={_v}: the unit is **%** (0-100). This looks like a 0-1 fraction "
                    f"- if you meant {_v * 100:.0f}, enter that number."
                )

        # ── 하드: 패딩 서열이 **새 제한효소 사이트를 만들지 않는가** ──
        # 패딩은 [pad][인식서열][spacer][overhang] 순으로 인식서열 **바로 5'** 에 붙는다.
        # 경계에서 새 사이트가 생기면 PCR 산물이 의도치 않은 곳에서 잘린다.
        # 기본값 AATACC는 6-mer 4096개를 전수 검사해 고른 값이지만, 유저가 바꾸면 검증이 없었다.
        if self.site_5prime_pad_seq and not (set(self.site_5prime_pad_seq.upper()) - set("ATGC")):
            errs.extend(self._check_pad_creates_site())

        if errs:
            raise ValueError("Design parameter error:\n  - " + "\n  - ".join(errs))

        # ── 소프트: 경고만 (차단하지 않음) ──
        warns: List[str] = []
        # ⚠️ 여기 있던 "판정창 ⊇ 설계창" 포함 검사 2건(Tm·GC)은 **제거**했다 —
        #    Y11에서 창을 하나로 합쳐 어긋남이 구조적으로 불가능해졌다. 포함 관계를
        #    검사해야 하는 상태 자체가 사라진 것이 이 작업의 결과다.
        lo, hi = _SOFT_TM_RANGE
        if self.tm_min < lo or self.tm_max > hi:
            warns.append(
                f"The Tm window [{self.tm_min}, {self.tm_max}]C is outside the usual range ({lo}-{hi}C)."
            )
        wlo, whi = _SOFT_TM_WIDTH
        width = self.tm_max - self.tm_min
        if width < wlo:
            warns.append(
                f"The Tm window is only {width:.1f}C wide (recommended {wlo:g}C or more) - "
                f"there may be no primer on this sequence that fits."
            )
        elif width > whi:
            warns.append(
                f"The Tm window is {width:.1f}C wide (recommended {whi:g}C or less) - "
                f"primers at the two ends would run under different conditions."
            )
        if self.primer_len_min < _SOFT_PRIMER_LEN_MIN:
            warns.append(
                f"primer_len_min={self.primer_len_min}bp is short (recommended {_SOFT_PRIMER_LEN_MIN}bp or more)."
            )
        glo, ghi = _SOFT_GC_RANGE
        if self.gc_min < glo or self.gc_max > ghi:
            warns.append(
                f"The GC window ({self.gc_min}-{self.gc_max}%) is wider than the usual range ({glo}-{ghi}%)."
            )
        olo, ohi = _SOFT_GIBSON_OVERLAP_RANGE
        if not (olo <= self.gibson_overlap <= ohi):
            warns.append(
                f"gibson_overlap={self.gibson_overlap}bp is outside the usual range ({olo}-{ohi}bp)."
            )

        # homology는 프라이머의 5' tail이 되므로 tail + anneal이 THERMO_MAX_LEN을 넘을 수 있다.
        # ⚠️ 2026-10-24(H1)로 **하드 에러 → 경고**로 내렸다. 그 전에는 "품질 계산이 아예
        #    불가능"해서 막았지만, 지금은 >60bp를 ViennaRNA로 계산한다(`secondary_structure`).
        #    다만 스케일이 달라 **경고 전용**이라 차단 게이트가 사라지는 건 사실이므로 알린다.
        #
        # yeast_homology는 여기서 보지 않는다 — 기본값(45)이 이미 이 구간이라 yeast를 쓰지
        # 않는 실행에도 경고가 떠서 소음이 된다. 긴 homology는 yeast HR 방법의 전제이고,
        # 실제로 상한을 넘는 프라이머만 `GibsonCloning._validate_assembly_primer`가
        # 길이와 함께 정확히 지적한다.
        if self.gibson_overlap + self.primer_len_min > THERMO_MAX_LEN:
            warns.append(
                f"gibson_overlap({self.gibson_overlap}) + primer_len_min({self.primer_len_min}) = "
                f"{self.gibson_overlap + self.primer_len_min}bp가 primer3 계산 상한 "
                f"exceeds {THERMO_MAX_LEN}bp - for this combination hairpin/dimer are judged "
                f"**warning-only** (ViennaRNA). The blocking gate does not run."
            )
        _type_cut = min(self.qc_hairpin_dg_3prime, self.qc_hairpin_dg_internal)
        if self.qc_dg_absolute > _type_cut:
            warns.append(
                f"qc_dg_absolute({self.qc_dg_absolute}) is less negative than the hairpin type cut "
                f"({_type_cut}) - weak structures would be reported at the top tier (strong)."
            )
        # ── 소프트: 패딩 서열 자체의 형질 (차단하진 않지만 프라이머 품질에 영향) ──
        warns.extend(self._check_pad_quality())
        return warns

    # ------------------------------------------------------------------ #
    # 패딩 서열 검증 (CP-C)
    # ------------------------------------------------------------------ #
    def _pad_variants(self) -> List[tuple]:
        """실제로 쓰일 (효소, 붙는 패딩, 인식서열) 조합.

        `resolve_site_pad`가 길이에 맞춰 자르고 반복하므로, **입력 서열 그대로가 아니라
        실제로 붙는 서열**을 검사해야 한다. 조립 효소는 부품 내부 사이트 회피로 바뀔 수
        있으므로 지원 효소 전체를 본다.
        """
        from basic.sequence_analysis import get_typeIIS_params
        out = []
        for enz in SUPPORTED_ENZYMES:
            try:
                site = str(get_typeIIS_params(enz)["fwd"])
            except Exception:
                continue
            out.append((enz, self.resolve_site_pad(enz).upper(), site.upper()))
        return out

    def _check_pad_creates_site(self) -> List[str]:
        """패딩이 인식서열과 만나 **새 사이트**를 만드는지. 만들면 하드 에러.

        탐색은 `basic.sequence_analysis.find_restriction_sites`가 **단일 구현**이다.
        (예전엔 여기서 따로 문자열 탐색을 했는데, 같은 일을 하는 코드가 두 벌이 되면서
         한쪽만 고쳐지는 일이 생겼다 — 2026-10-10에 sequence_analysis로 합쳤다.)
        인식서열이 서열 끝에 있을 때 Biopython `search()`가 놓치는 함정은 그쪽에 문서화돼 있다.
        """
        from basic.sequence_analysis import find_restriction_sites, get_typeIIS_params

        errs: List[str] = []
        seen = set()
        for enz, pad, site in self._pad_variants():
            if not pad:
                continue
            joined = pad + site
            for probe in SUPPORTED_ENZYMES:
                # 인식서열 단독에 이미 있는 개수보다 늘어나면 패딩이 만든 것이다.
                before = len(find_restriction_sites(site, probe).get(probe, []))
                after = len(find_restriction_sites(joined, probe).get(probe, []))
                if after <= before:
                    continue

                # 아래는 **메시지 표기용**이다. 생겼는지 여부의 판정은 위에서 끝났다.
                try:
                    pf = str(get_typeIIS_params(probe)["fwd"]).upper()
                    pr = str(get_typeIIS_params(probe)["rev"]).upper()
                except Exception:
                    pf = pr = ""
                if pf and joined.count(pf) > site.count(pf):
                    strand, pat = "forward", pf
                else:
                    strand, pat = "reverse", pr

                key = (probe, strand)
                if key in seen:
                    continue
                seen.add(key)
                errs.append(
                    f"site_5prime_pad_seq={self.site_5prime_pad_seq!r}: placed before the recognition "
                    f"sequence it **creates a new {probe} site ({strand} {pat})** "
                    f"(e.g. {pad}|{site} for {enz}) - the PCR product would be cut where you did not intend."
                )
        return errs

    def _check_pad_quality(self) -> List[str]:
        """패딩 서열 자체의 형질 — 경고만(동작은 한다)."""
        from Bio.Seq import Seq
        warns: List[str] = []
        seq = (self.site_5prime_pad_seq or "").upper()
        if not seq or (set(seq) - set("ATGC")):
            return warns
        # ① self-complementary면 프라이머 5' 끝에서 hairpin이 된다.
        if seq == str(Seq(seq).reverse_complement()):
            warns.append(
                f"site_5prime_pad_seq={seq!r} is self-complementary - it forms a hairpin at the "
                f"primer 5' end."
            )
        # ② 동종 염기 3연속 — 기본값 선정 때 배제했던 조건. spacer 더미('A'*n)와 이어지면
        #    poly-A 구간이 길어진다.
        for i in range(len(seq) - 2):
            if seq[i] == seq[i + 1] == seq[i + 2]:
                warns.append(
                    f"site_5prime_pad_seq={seq!r} has three identical bases in a row ({seq[i] * 3}) - "
                    f"synthesis and amplification quality may suffer."
                )
                break
        # ③ GC 극단 — 기본값 선정 기준은 33~50%였다.
        gc = 100.0 * sum(c in "GC" for c in seq) / len(seq)
        if not (20.0 <= gc <= 60.0):
            warns.append(
                f"site_5prime_pad_seq={seq!r} has {gc:.0f}% GC "
                f"(the default was chosen within 33-50%)."
            )
        return warns

    # ------------------------------------------------------------------ #
    # 효소에서 파생되는 값 (자유 파라미터가 아님)
    # ------------------------------------------------------------------ #
    def enzyme_params(self) -> Dict[str, Any]:
        """선택한 효소의 설계 파라미터를 Biopython에서 파생.

        Returns:
            {'fwd', 'rev', 'spacer_len', 'def_spacer', 'overhang_len'}
        """
        from basic.sequence_analysis import get_typeIIS_params  # 지연 import (순환 방지)
        return get_typeIIS_params(self.gg_enzyme)

    def tm_kwargs(self) -> Dict[str, Any]:
        """`basic.sequence_analysis.calculate_tm`에 넘길 버퍼 조건 한 벌.

        인자가 7개라 호출부마다 풀어 쓰면 한 곳만 빠뜨려도 **그 자리만 다른 자로**
        재게 된다 — 2026-11-05에 바로 그 "자가 두 개" 문제를 고쳤으므로 되살리지 않는다.
        """
        return {
            "Na": self.na_mM, "dnac1": self.primer_conc_nM,
            "K": self.k_mM, "Tris": self.tris_mM, "Mg": self.mg_mM,
            "dNTPs": self.dntp_mM, "saltcorr": self.tm_saltcorr,
        }

    def ta_recommended(self, tm1: float, tm2: float) -> Optional[float]:
        """NEB 권장 어닐링 온도. 조건이 등록된 polymerase가 아니면 None.

        **두 프라이머 중 낮은 Tm**을 기준으로 한다(NEB 규칙). 규칙과 검증 근거는
        `NEB_TA_RULES` 주석 참고 — 자체 검증에서 일치율 95~100%였다.
        """
        name = self.polymerase_key()
        if name is None:
            return None
        kind, off, cap = NEB_TA_RULES[name]
        lo = min(tm1, tm2)
        val = 0.93 * lo + 7.5 if kind == "phusion" else lo + off
        return min(val, cap)

    def polymerase_key(self) -> Optional[str]:
        """현재 버퍼 조건이 어느 NEB polymerase와 같은지. 없으면 None.

        유저가 advanced에서 값을 직접 바꾸면 어느 제품도 아니게 되는데, 그때
        Ta를 특정 제품 규칙으로 계산하면 **근거 없는 숫자**가 된다 → None을 돌려
        호출부가 "권장 Ta 없음"으로 처리하게 한다.
        """
        for name, cond in NEB_TM_CONDITIONS.items():
            if all(abs(getattr(self, k) - v) < 1e-9 for k, v in cond.items()):
                return name
        return None

    @property
    def tm_opt(self) -> float:
        """primer3가 요구하는 '목표 Tm' — 선언 범위의 중앙값.

        구버전 `target_tm`의 자리다. 범위를 선언하면 목표는 따라 나오므로 필드가 아니다.
        """
        return (self.tm_min + self.tm_max) / 2.0

    @property
    def overhang_len(self) -> int:
        """점착말단 길이 — 효소가 정한다 (BsaI/BbsI/BsmBI 4, SapI 3)."""
        return int(self.enzyme_params()["overhang_len"])

    def resolve_site_pad(self, enzyme: Optional[str] = None) -> str:
        """인식서열 앞에 실제로 붙일 여분 서열.

        Args:
            enzyme: 효소명. 생략 시 `gg_enzyme`.
                (부품 내부 사이트 회피로 조립 효소가 바뀔 수 있어 호출자가 실제 효소를 넘긴다)

        Returns:
            붙일 서열. 길이 0이면 빈 문자열.
        """
        enz = enzyme or self.gg_enzyme
        n = self.site_5prime_pad
        if n is None:
            n = SITE_PAD_DEFAULT.get(enz, SITE_PAD_FALLBACK)
        if n <= 0:
            return ""
        seq = self.site_5prime_pad_seq
        if not seq:
            return ""
        # 인식서열에 가까운 쪽(서열의 3' 끝)을 남긴다 — 사이트 신규생성 검사를 그 기준으로 했으므로.
        rep_seq = seq * ((n // len(seq)) + 1)
        return rep_seq[-n:]

    @property
    def site_spacer_len(self) -> int:
        """인식서열과 절단점 사이 더미 길이 — 효소가 정한다 (BsaI 1, BbsI 2)."""
        return int(self.enzyme_params()["spacer_len"])

    # ------------------------------------------------------------------ #
    # 직렬화
    # ------------------------------------------------------------------ #
    def canonical(self) -> str:
        """해시·비교용 정규화 문자열.

        `derived_id`가 artifact id를 만들 때 재료로 쓴다. 파라미터가 이 재료에 안 들어가면
        **다른 파라미터로 만든 산출물이 같은 id를 갖고 서로 덮어쓴다.**

        같은 값이면 항상 같은 문자열이어야 하므로 필드명 정렬 + 부동소수 자리수 고정.
        """
        parts = []
        for f in sorted(fields(self), key=lambda x: x.name):
            v = getattr(self, f.name)
            parts.append(f"{f.name}={v:.6f}" if isinstance(v, float) else f"{f.name}={v}")
        return ";".join(parts)

    def to_dict(self) -> Dict[str, Any]:
        """galaxy가 DB에 저장할 때 쓰는 표현 (엔진↔galaxy 직렬화 계약)."""
        return {f.name: getattr(self, f.name) for f in fields(self)}

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "DesignParams":
        """`to_dict()` 결과로부터 복원. 빠진 키는 기본값을 쓴다.

        모르는 키는 **조용히 무시하지 않고 에러**를 낸다 — 오타나 버전 불일치를 삼키면
        "입력했는데 반영이 안 되는" 가장 찾기 어려운 버그가 된다.
        """
        known = {f.name for f in fields(cls)}
        unknown = set(data) - known
        if unknown:
            raise ValueError(
                f"Unknown parameter: {sorted(unknown)} (allowed: {sorted(known)})"
            )
        return cls(**data)

    @classmethod
    def from_preset(cls, name: str, **overrides: Any) -> "DesignParams":
        """polymerase 프리셋으로 한 벌 만들기. 프리셋에 없는 필드는 기본값을 상속한다.

        Args:
            name: `PRESETS`의 키 (대소문자·앞뒤공백 무시). 예: "q5", "Phusion", "taq"
            **overrides: 프리셋 위에 덮어쓸 값. 오타는 `from_dict`와 같은 규칙으로 에러.

        Returns:
            DesignParams

        Raises:
            ValueError: 모르는 프리셋 이름이거나, `overrides`에 모르는 키가 있을 때.

        Examples:
            >>> DesignParams.from_preset("q5")
            >>> DesignParams.from_preset("q5", gibson_overlap=25)   # 프리셋 + 랩 사정
        """
        key = (name or "").strip().lower()
        if key not in PRESETS:
            raise ValueError(
                f"Unknown preset: {name!r} (available: {sorted(PRESETS)}). "
                f"To add one, put it in params.PRESETS."
            )
        # 프리셋 → overrides 순으로 덮어쓰고, 나머지는 from_dict가 기본값으로 채운다.
        # (from_dict를 거치므로 모르는 키 검사가 한 곳에서만 일어난다)
        return cls.from_dict({**PRESETS[key], **overrides})

    def with_(self, **changes: Any) -> "DesignParams":
        """일부 값만 바꾼 새 인스턴스 (frozen이라 직접 대입 불가)."""
        return replace(self, **changes)


# 아무것도 지정하지 않았을 때 쓰는 한 벌 (= 기존 하드코딩 동작)
DEFAULT_PARAMS = DesignParams()
