# Bioparts DB 서열 방향(strand) 정합성 — 참조 문서

> 작성 2026-07-09. **자주 헷갈리는 주제**라 확정 근거를 남긴다.
> 핵심 질문: "DB 부품 서열은 전부 정방향(+)인가? part_intake의 서열 매칭에 역상보(RC) 처리가 필요한가?"
> **결론: DB는 전부 +(top strand 5'→3'). RC 매칭 불필요(part_intake fix5 드롭).**

---

## 0. 30초 요약
- DNA 서열은 **항상 5'→3'(top strand)로 적는다.** "3'→5' 저장" 같은 건 존재하지 않는다.
- 크롤러가 SynBioHub의 `<sbol:elements>`(= 그 부품의 top strand 서열)를 **verbatim** 저장 → DB `Sequence` 컬럼 = **그 엔트리의 +가닥**.
- 따라서 **모든 DB 파트/모듈 = 자기 자신의 정방향(+)**. 어떤 엔트리도 "자기 역상보"로 저장돼 있지 않다.
- RC(역상보) 변환은 **저장이 아니라 "읽을 때"** 일어난다 — gbk의 feature가 `strand=-1`이면 Biopython `feat.extract`가 자동으로 RC해 그 feature의 5'→3'로 돌려준다. → 매칭이 방향 무관하게 성립.
- 그래서 part_intake의 **RC 매칭(fix5)은 불필요 → 드롭**. (오히려 RC 자동매칭은 다른 부품 오인식 위험)

---

## 1. 크롤러가 서열을 얻는 방식 (근거)
출처: SynBioHub(iGEM/Cello/Bacillus/Cryptococcus 등) + BioPartsDB. 형식 SBOL/FASTA.
크롤 코드: `archive/docs/analysis/Bioparts_DB/code/`, `.../data_analyze/analyze.py`.

| 스크립트 | 서열 출처 | 방향 처리 |
|---|---|---|
| `fetch_collection_fasta.py` | SynBioHub `/fasta` 엔드포인트 | `r.text` **그대로** (정규화 없음) |
| `analyze.py` (SynBioHubFinalScraper) | SBOL XML `<sbol:elements>` | `seq_tag.text` **그대로** |
| `fetch_*_metadata.py` | (서열 아님 — Version/날짜만) | — |

→ **어디에도 RC/방향 정규화 코드가 없다.** DB 서열 = SBOL이 정의한 부품 canonical 서열 그대로.
→ 실측: DB의 `BBa_K1480002.gb` ORIGIN이 SBOL `BBa_K1480002_sequence`의 `<sbol:elements>`와 **글자까지 일치**(verbatim 저장 확인).

---

## 2. "reverse"라고 적혀 있어도 서열은 정방향인 이유 (핵심 뉘앙스)

SBOL에는 방향 필드가 **두 개**이고 서로 다를 수 있다:
- `sbol:orientation` = `inline`(정방향, 같은 가닥) | `reverseComplement`(역가닥) — **SBOL 표준의 진짜 방향**.
- `igem:direction` = `forward` | `reverse` — iGEM 옛 위키 레거시 라벨(부정확할 수 있음).

### 실측 사례: BBa_K1480002 안의 서브파트 BBa_E1010 (mRFP1)
- XML 라벨: `igem:direction = **reverse**`  ← 얼핏 "역방향"으로 보임
- 그러나: `sbol:orientation = **inline**`(정방향), 그리고 **실제 composite 서열에 E1010이 정방향으로 박혀 있음**.
- 검증(코드): `E1010 정방향 ∈ composite → True` / `E1010 역상보 ∈ composite → False`.
- .gb 표기도 `CDS 84..758`(complement 없음 = forward).

**결론**: `reverse`(igem:direction)는 **iGEM→SBOL 변환의 부정확한 잔재 라벨**일 뿐, **저장 서열의 방향을 뒤집지 않는다.** 심판은 라벨이 아니라 **실제 서열**이고, 서열은 정방향이다. (mRFP1은 프로모터+RBS 뒤에서 발현돼야 하므로 정방향이 생물학적으로도 맞다.)

### "reverse"의 진짜 의미
어떤 feature가 **진짜** 역가닥(`reverseComplement`)이면 → 부모 top strand 서열(여전히 5'→3')엔 그 feature의 **역상보(RC)**가 나타난다. 서열이 3'→5'가 되는 게 아니라 top strand에 RC 형태로 박힌다.

---

## 3. RC(자동변환)는 "읽을 때"만 일어난다
- **저장**: 항상 그 엔트리의 top strand 5'→3'. (예: E1010 엔트리 `Sequence` = `atggcttcctcc...` 정방향)
- **읽기(gbk 추출)**: gbk의 feature가 `strand=-1`(진짜 역가닥)이면 `feat.extract(rec)`가 그 구간 top strand를 **RC해서 feature 자신의 5'→3'로 되돌려 준다.** → part_intake는 이 정규화된 5'→3'로 DB(전부 +)와 매칭 → **방향 무관하게 정확히 matched.**
- 실측(part_intake): 프로모터 BBa_I0500을 gbk에 **마이너스 가닥(strand=-1)**으로 넣어도 `parts_from_genbank` → `matched, BBa_I0500`.

---

## 4. 그래서 part_intake fix5(RC 매칭)를 드롭한 이유
- DB가 전부 +, 읽기 쪽 `feat.extract`가 strand를 정규화 → **RC 매칭을 추가할 필요가 없다.**
- 오히려 서열을 무조건 RC까지 조회하면 **다른 부품을 오인식**할 위험(예: `RC(BBa_I0500)`이 우연히 별개 부품 `BBa_K607040`과 일치). 방향이 기능적으로 의미 있는 파트(promoter/rbs/cds/terminator)에서 blind RC 매칭은 위험.

---

## 5. 남은 이론적 리스크 (탐지 한계)
- DB 파트 테이블엔 **strand/orientation 컬럼이 없다**(`ID/Name/Sequence/Description/Institution/Version/Date`). 방향 메타데이터 자체가 저장 안 됨.
- 유일한 이론적 불일치: 서로 다른 소스가 **같은 생물학적 요소를 다른 canonical 방향으로 정의**한 경우(교차 소스) — DB 내부 불일치가 아니라 소스 간 명명 차이, 드묾.
- 실측 감사가 필요하면: DB에서 **"부품 서열 & 그 역상보가 둘 다 존재"하는 쌍**을 스캔(방향 중복 신호). (예: I0500↔K607040)

---

## 6. 한 줄 결론
**Bioparts DB의 모든 파트/모듈 서열 = 그 엔트리의 top strand 5'→3'(=+). 저장은 항상 정방향, RC는 gbk 읽을 때만. → 매칭은 안전하고 RC 매칭은 불필요.**
