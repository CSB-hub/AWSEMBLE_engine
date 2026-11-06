# CLAUDE.md

이 파일은 Claude Code(claude.ai/code)가 **현행 `awsemble/src/` 코드**를 다룰 때 참고하는 가이드입니다.

> 📌 이 문서는 **DB(Bioparts/Test 스키마) 기반의 현재 코드**를 기준으로 합니다.
> 과거 Excel 설계 파일 기반(v1/v2, `assemble_pattern_modules.py`,
> `PartPreparationModule` 등)의 레거시 설명은 `archive/CLAUDE.md`를 참고하세요.

## User 요청사항
- 파일들을 생성할 때는 한글로 작성해줘. (문서화 파일과 코드 주석 모두 해당)
- 존댓말이 아니고 반말로 대답해줘.
- 대답할 때 "내가"라는 단어를 사용해줘. 예: "내가 도와줄게"

## 문서화 규칙 (중요)
- **`docs/AWSEMBLE_GUIDE.md` = 프로젝트 현황 단일 진실원천(SSOT).** `develop_log/`에 새 작업·결과를 기록하거나 갱신할 때는 **같은 커밋에서 `docs/AWSEMBLE_GUIDE.md`의 §6 진행상황·로드맵도 함께 갱신**한다. (가이드만 보면 "지금 어디까지 / 무엇이 남았나"가 항상 파악되는 상태를 유지)

## 코드 재활용 원칙 (중요)

**같은 일을 하는 코드를 두 벌 만들지 않는다.** 문서에 SSOT가 있듯 코드에도 있어야 한다.

- 기능을 추가하기 전에 **먼저 찾는다**: `basic/`(서열 프리미티브), `assembly/common.py`(조립 공통), `utils/io.py`(DB·파일 입출력)에 이미 있는지.
- 비슷한 게 있는데 조금 다르면, **인자로 차이를 흡수**할 수 있는지 먼저 본다. 새 함수는 마지막 수단이다.
- **계층을 지킨다.** 서열 자체를 다루는 일은 `basic/`이 집이고 `params.py`·`assembly/`는 **소비자**다. 소비자 쪽에 프리미티브를 구현하면 그게 곧 두 번째 사본이 된다.
- 버그를 고칠 때는 **같은 일을 하는 다른 코드가 있는지 먼저 확인**한다. 한쪽만 고치면 안 고친 쪽이 어디서 쓰이는지 모른 채 남는다.

**왜 중요한가 — 2026-10-10에 실제로 터진 사고:**

제한효소 인식서열 탐색이 `params.py`와 `basic/sequence_analysis.py`에 **두 벌** 있었다.
2026-09-23에 "Biopython `Restriction.search()`는 절단 위치를 돌려주므로 서열 끝쪽 인식서열을
놓친다"는 걸 발견하고 **`params.py` 쪽만** 고쳤다. `sequence_analysis.py`는 그대로 뒀는데,
하필 **실제 조립 경로(`assembly/golden_gate.py`의 효소 자동 전환)가 쓰는 게 안 고친 쪽**이었다.

결과: 운영 Bioparts 11227개 중 35개가 효소 판정을 잘못 받고 있었다. 내부에 인식서열이 있는데
"깨끗하다"고 판정 → 그 효소로 조립 진행 → **에러 없이 조용히 잘린 결과**가 나올 수 있었다.
사본이 하나였다면 2026-09-23의 수정 한 번으로 끝났을 일이다.

---

## 프로젝트 개요

**awsemble**는 조합 대사 경로 설계를 위한 DNA 조립 자동화 프레임워크다. Golden Gate Assembly와 Gibson Assembly를 사용해 DNA 구조물을 자동 설계·조립하는 객체지향 도구를 제공한다.

가장 큰 특징은 **부품(Part) 정보를 Excel이 아니라 MySQL 데이터베이스에서 직접 가져온다**는 점이다. 표준화된 부품 레지스트리(Bioparts 스키마)와 조립 정의(Test 스키마)에서 서열과 조립 순서를 끌어와, 오버행·프라이머를 자동 설계하고 조립한 뒤, DB의 정답 서열과 대조 검증하는 한 사이클을 자동화한다.

생물학적 맥락에서 각 계층의 의미:
- **Part**: 프로모터, RBS, CDS, 터미네이터 같은 최소 유전 요소
- **Module**: Part들을 이어 만든 발현 카세트 (예: 프로모터+RBS+CDS+터미네이터)
- **Multi-module / Pathway**: Module들을 이어 만든 대사 경로 전체 구조물

---

## 현행 디렉토리 구조

```
awsemble/
├── .env                      # MySQL 접속 정보 (HOST/PORT/USER/PASSWORD)
├── requirements.txt
├── src/
│   ├── basic/                # 기본 생물학 분석 (DNA 검증, Tm/GC, 프라이머)
│   │   ├── sequence_analysis.py
│   │   ├── secondary_structure.py
│   │   └── primer_design.py
│   ├── assembly/             # 조립 방법 (Strategy + Factory 패턴)
│   │   ├── base.py           # AssemblyBase (추상 클래스)
│   │   ├── factory.py        # AssemblyFactory
│   │   ├── common.py         # feature 처리 유틸
│   │   ├── golden_gate.py    # GoldenGateModule
│   │   └── gibson.py         # GibsonCloning
│   ├── workflow/             # 워크플로우 오케스트레이션
│   │   ├── base.py           # WorkflowBase
│   │   ├── biopart_assemble.py           # BiopartWorkflow (핵심)
│   │   ├── classical_assemble.py         # ClassicalWorkflow
│   │   └── builder.py        # WorkflowBuilder
│   └── utils/
│       └── io.py             # DB 접속, 파일 I/O, 검증
└── tests/
    ├── igem/test_single.py   # 단일 Test_ID 실행/검증
    ├── igem/test_all.py      # iGEM 전수(9999) 병렬 테스트
    └── kribb/test_kribb.py   # KRIBB 기관 Test_ID 실행/검증
```

---

## 핵심 아키텍처

### 클래스 계층 구조 (현행)

```
basic/
  ├── PrimerDesign / ValidatePrimer        # 프라이머 설계·검증
  └── (함수) calculate_tm, calculate_gc_content,
            find_restriction_sites,
            check_restriction_site_in_part_sequence,
            analyze_single_strand_structure ...

assembly/
  AssemblyBase (ABC: prepare_fragments / assemble / validate)
    ├── GoldenGateModule   # Type IIS 제한효소 (BsaI/BbsI)
    └── GibsonCloning      # 상동성 overlap (15-25bp)
  AssemblyFactory          # "gibson" | "golden_gate" 문자열로 인스턴스 생성

workflow/
  WorkflowBase
    ├── BiopartWorkflow       # 6단계 (현재 Step3 Module 조립이 핵심 구현)
    └── ClassicalWorkflow     # insert + vector 단순 클로닝
  WorkflowBuilder          # create_classical / create_biopart 진입점

utils/
  DatabaseImport           # MySQL(Bioparts/Test) 조회
  DNAFileReader / DNAFileWriter   # FASTA/GenBank I/O
  DNAValidator             # 조립 결과 ↔ DB 정답 대조 검증
```

### 디자인 패턴
- **Factory**: `AssemblyFactory.create_assembler("golden_gate" | "gibson")`로 조립 방법 선택을 문자열로 추상화. 새 방법 추가 시 `factory.py` 한 곳만 수정 (개방-폐쇄 원칙).
- **Strategy / Dependency Injection**: 워크플로우는 구체 클래스가 아니라 `AssemblyBase` 인터페이스에 의존. 런타임에 조립 방법 교체 가능.
- **Lazy import**: `factory.py`에서 순환 참조를 피하기 위해 assembler를 지연 로딩.

---

## 데이터베이스 구조 (핵심)

부품·조립 정의·정답 서열이 모두 MySQL에 있다. 접속 정보는 `awsemble/.env`:

```bash
HOST=...        # DB 호스트
PORT=3306
USER=...
PASSWORD=...
```

### 스키마 구성

**`Bioparts` 스키마** — 부품 서열 저장소
- `part_type`별 테이블(예: promoter, rbs, cds, terminator ...): 컬럼 `ID`, `sequence`
- `overhang` 테이블: Golden Gate용 high-fidelity 오버행 후보군 (`institution = 'Awsemble'`)

**`Test` 스키마** — 조립 정의 + 정답
- `test_construct`: `Test_ID`, `Institution`(igem/KRIBB), `sequence`(정답 서열)
- `test_part`: `Test_ID`, `part_id`, `part_type`, `part_order`(조립 순서)

### DB 접근 진입점 (`utils/io.py`의 `DatabaseImport`)
| 메서드 | 역할 |
|---|---|
| `get_ordered_parts_by_test_id(test_id)` | `test_part`에서 순서대로 `(part_id, part_type)` 목록 |
| `get_part_by_id(schema, table, part_id)` | `Bioparts`에서 실제 서열 → `Dseqrecord` |
| `get_construct_sequence_by_test_id(test_id)` | 정답 서열 조회 (검증용) |
| `get_all_test_ids_from_construct()` | igem/KRIBB Test_ID 전체 |
| `get_test_ids_by_institution(inst)` | 특정 기관 Test_ID |
| `get_all_overhang_ids()` | Awsemble 오버행 ID 목록 |

---

## Biopart 워크플로우 (현행)

`workflow/biopart_assemble.py`의 `BiopartWorkflow`는 Design-Build-Test-Learn 6단계를 메서드로 구현한다. **핵심 차이점: 부품을 Excel에서 읽지 않고, 이미 만들어진 `Dseqrecord` 리스트를 생성자에 주입받는다.**

```python
workflow = BiopartWorkflow(project_dir=..., parts=final_dseq_records)
```

### 단계별 현황
- `run_step1_part_design()` — **`pass` (빈 스텝)**. Part 준비는 워크플로우 *앞단*(DB 조회 + `GoldenGateModule.prepare_fragments` 자동 설계)에서 끝나므로 내부 구현 없음.
- `run_step2_part_cloning()` — 현재 스킵.
- **`run_step3_module_assembly(module_df, method)` — 핵심 구현부.** `module_df`(컬럼: `id`, `name`)를 받아, 메모리의 `self.parts`(이름→Dseqrecord 딕셔너리)에서 부품을 순서대로 꺼내 조립한다.
  - `method="golden_gate"`: 각 부품을 `extract_insert_fragment()`로 절단·추출 → `assemble()`
  - `method="gibson"`: `part_dseqs`를 바로 `assemble()`
  - 결과는 `modules/inserts/{module_id}.gb`로 저장
- `run_step4`~`run_step6` — 골격만 있고 본문 미구현(`# ...`).

### 표준 실행 흐름 (tests/igem/test_single.py 기준)
```python
db = DatabaseImport()

# 1. DB에서 부품 목록(순서 보존) + 서열 조회
ordered = db.get_ordered_parts_by_test_id(test_id)        # [(part_id, part_type), ...]
parts = [db.get_part_by_id("Bioparts", ptype, pid) for pid, ptype in ordered]

# 2. 오버행·프라이머 자동 설계 → PCR 산물(full_sequence) 생성
gg = GoldenGateModule()
design = gg.prepare_fragments(parts, assembly_type="linear")

# 3. 설계 결과를 Dseqrecord로 변환 (원본 서열/효소 정보를 annotations에 보관)
records = []
for res in design:
    r = Dseqrecord(res['full_sequence'], name=f"{res['part_name']}_designed")
    r.annotations['target_seq'] = res['target_seq']           # 절단 후 조각 선택용 대조군
    r.annotations['restriction_enzyme'] = res['restriction_enzyme']
    records.append(r)

# 4. 워크플로우 Step3 실행
wf = BiopartWorkflow(project_dir=PROJECT_DIR, parts=records)
module_df = pd.DataFrame({"id": [test_id]*len(records), "name": [r.name for r in records]})
wf.run_step3_module_assembly(module_df=module_df, method="golden_gate")

# 5. DB 정답과 대조 검증
DNAValidator(db).check_parts_correct(output_gb_path, test_id)
```

---

## Golden Gate Assembly 구현 상세 (`assembly/golden_gate.py`)

생물학적 원리: BsaI(`GGTCTC`)/BbsI(`GAAGAC`) 같은 **Type IIS 제한효소**는 인식 서열 *바깥*을 잘라, 인식 서열은 떨어져 나가고 4bp **점착성 말단(overhang)**만 남긴다. 이 4bp가 상보적인 조각끼리만 ligation되므로 여러 부품을 한 번에 순서대로 조립할 수 있다.

### 주요 메서드 ↔ 생물학적 의미
| 메서드 | 역할 |
|---|---|
| `prepare_fragments(parts, vector, assembly_type)` | 부품 전수 스캔 → 효소 결정 → 오버행 할당 → 프라이머 설계 |
| `_get_high_fidelity_overhang_set(count)` | DB overhang 풀에서 서로 해밍거리가 먼 세트 선택 (mis-ligation 방지) |
| `extract_insert_fragment(dseq)` | 효소 절단 후 여러 조각 중 원본과 가장 유사한 insert를 선택 |
| `assemble(fragments, assembly_type)` | 인접 조각의 4bp 말단이 역상보(RC)인지 확인 후 결합 |
| `validate(result)` | 결과물에 BsaI/BbsI 인식 서열이 남았는지 확인 |
| `extract_vector_overhangs(vector)` | 벡터를 BsaI로 잘라 양 끝 오버행 추출 (circular용) |

### 효소 자동 전환 로직 (중요)
부품 *내부*에 조립용 효소의 인식 서열이 우연히 들어있으면 의도치 않은 위치에서 잘려 조립이 깨진다. `prepare_fragments()`는 조립에 들어가는 **모든 부품을 전수 스캔**(`check_restriction_site_in_part_sequence`)하여:
- 내부에 BsaI 사이트 발견 → **조립 전체를 BbsI로 전환**
- BsaI·BbsI가 **둘 다** 발견 → `RuntimeError` (BsmBI 등 제3 효소 권고)

효소별 파라미터는 `primer_design.py`의 `enzyme_sites`에 정의 (BsaI: spacer 1bp, BbsI: spacer 2bp — 절단 칸 수가 다르기 때문).

---

## Gibson Assembly 구현 상세 (`assembly/gibson.py`)

생물학적 원리: 제한효소 없이, exonuclease가 단일가닥을 노출시키고 인접 fragment의 **상동 overlap(15-40bp)**끼리 어닐링한 뒤 polymerase+ligase가 봉합한다.

| 메서드 | 역할 |
|---|---|
| `prepare_fragments(insert, vector, ...)` | insert 증폭 프라이머에 벡터 말단 서열을 덧붙여 overlap 생성 |
| `amplify_insert(insert, fwd, rev)` | pydna `pcr()`로 PCR 시뮬레이션 |
| `assemble(fragments, assembly_type)` | pydna `Assembly`(terminal_overlap)로 overlap 자동 탐색·조립 |
| `_validate_primer_length(...)` | 짧은 insert(<50bp)에서 프라이머가 템플릿보다 길어지는 문제 보정 |

---

## 프라이머 설계·검증 (`basic/primer_design.py`)

`PrimerDesign`/`ValidatePrimer`는 "이 프라이머로 PCR이 실제로 될까?"를 코드화한 생물학적 품질 게이트다. **primer3-py 라이브러리(`import primer3`)에 의존**한다.

- `PrimerDesign.design_basic_primers(...)` — primer3로 순수 증폭 프라이머 설계
- `PrimerDesign.design_goldengate_primers(...)` — 양 끝에 `[효소 인식서열]+[spacer]+[overhang]`을 붙여 절단 후 정확한 overhang이 노출되도록 설계. `start_structure`/`end_structure`/`full_sequence` 반환
- `ValidatePrimer.validate_primer(...)` — Tm/GC%/self-dimer/hairpin/길이 종합 검증 → `passed`/`violations`
- 보조: `check_primer_dimer`(ΔG 기반 dimer), `check_gc_clamp`(3' 안정성), `check_primer_specificity`(off-target)

---

## 서열 분석 (`basic/sequence_analysis.py`)

- `validate_dna_sequence`, `validate_cds` — 입력 서열/ORF 유효성 검증
- `calculate_gc_content`, `calculate_tm` — GC%, 융해 온도
- `find_restriction_sites(seq, enzyme)` — Type IIS(BsaI/BbsI/BsmBI/Esp3I/PaqCI/SapI) + 범용 효소(EcoRI/BamHI/HindIII/XbaI/SpeI/NheI...) 마스터 목록 검색
- `check_restriction_site_in_part_sequence(seq)` — "부품 내부에 BsaI/BbsI가 있으면 나쁜 것"으로 해석하는 part 검증 래퍼 (Golden Gate 효소 전환의 입력)
- `find_repeats`, `get_sequence_stats` — 반복서열, 통계

---

## 검증 시스템 (`utils/io.py`의 `DNAValidator`)

`check_parts_correct(gbk_path, test_id)`:
1. 조립 결과 GenBank를 읽음
2. `get_ordered_parts_by_test_id`로 DB의 부품 순서 조회
3. 각 부품의 DB 서열이 결과물 안에 **순서대로** 등장하는지 위치 추적 검증

테스트 스크립트는 추가로 `get_construct_sequence_by_test_id`로 가져온 정답 서열과 `pairwise2` 정렬로 100% 일치 여부까지 비교한다.

---

## 환경 설정 및 실행

### 의존성 (`requirements.txt` 핵심)
- `pydna==5.5.3` — DNA 조립 시뮬레이션
- `biopython==1.85` — 서열/제한효소 조작
- `mysql-connector-python==9.5.0` — DB 접속
- `python-dotenv==1.1.1` — `.env` 로드
- `primers==0.5.10`, primer3 — 프라이머 설계 (※ `import primer3` 사용)
- 검증/분석: `dnachisel`, `pytrf`, `pydivsufsort` (※ 2차구조는 `primer3`로 통일 — 구 `seqfold`는 API 오독으로 hairpin 탐지가 죽어 있어 제거, 아래 secondary_structure 항목 참고)
- (레거시 흔적) `pandas`, `openpyxl` — 현재는 DataFrame 구성·일부 I/O에만 사용

### 실행
```bash
# 단일 Test_ID
cd awsemble && python tests/igem/test_single.py <Test_ID>

# iGEM 전수 병렬 테스트 (ThreadPoolExecutor, MAX_WORKERS=20)
python tests/igem/test_all.py

# KRIBB 기관
python tests/kribb/test_kribb.py
```
> `src/`가 import 경로에 잡혀야 한다(`from utils.io import ...`, `from workflow.biopart_assemble import ...`). PYTHONPATH에 `awsemble/src`를 추가하거나 해당 위치에서 실행.

---

## 중요한 구현 세부사항

### 파일 순서와 재현성 (필수)
`os.listdir()`/`glob`은 비결정적 순서를 반환할 수 있으므로 항상 `sorted()`로 고정한다. (`DNAFileReader.read_genbank_dir`은 `sorted(directory.glob("*.gb"))` 사용.) 조합 조립 결과의 재현성을 위한 프로젝트 핵심 원칙.

### annotations로 원본 정보 운반
Golden Gate에서 절단 후 올바른 insert 조각을 고르려면 원본 서열이 필요하다. `prepare_fragments` → Dseqrecord 변환 시 `annotations['target_seq']`, `annotations['restriction_enzyme']`에 보관하고, `amplify_insert`/`extract_insert_fragment`가 이를 읽어 유사도 비교(`SequenceMatcher(autojunk=False)`)에 사용한다.

### Feature 중복 제거 / 위치 조정
- `assembly/common.py`의 `remove_duplicate_features(dseq)` — (start,end,type,label) 기준 중복 feature 제거
- `adjust_feature_positions(dseq, shift)` — 프라이머로 길이가 변하면 feature 위치를 shift

### 알려진 이슈 / TODO (개발 로그 기준)
`develop_log/nyb/debugging/` 참고. 최신 전수 테스트(260319): 성공 9222 / 추출 실패 155 / 제3효소 추천 621.
- [ ] 제한효소 site + overhang이 겹쳐 의도치 않은 site가 생기는 케이스 (`BBa_K561005`, `BBa_M31790`)
- [ ] 짧은 부품 절단 시 pydna `Sequence is linear, origin can only be shifted` 에러 (`BBa_B0034`) — 안정적 해법 분석 필요
- [ ] overhang 불일치 재현성 (`BBa_J22203`) 원인 분석
- [ ] Circular 모드 전 과정 통합 테스트 / 대규모 조립 overhang 충돌 방지 검증

---

## 코딩 표준 및 개발 가이드

### 네이밍 (PEP 8)
```python
gene_sequence = "ATCG"          # 변수/함수: snake_case
class GoldenGateModule: ...      # 클래스: PascalCase
MAX_SEQUENCE_LENGTH = 10000      # 상수: UPPER_SNAKE_CASE
self._connection = None          # protected: _leading_underscore
```

### 매직 넘버 제거
```python
# 나쁜 예
if relative_score >= 0.8: ...
# 좋은 예
SIMILARITY_THRESHOLD = 0.8
if relative_score >= SIMILARITY_THRESHOLD: ...
```

### Docstring (Google 스타일)
```python
def design_goldengate_primers(insert, left_overhang, right_overhang, ...):
    """Golden Gate Assembly용 프라이머 설계 및 완성된 서열 생성

    Args:
        insert: 증폭할 DNA (Dseqrecord)
        left_overhang: 왼쪽 overhang (4bp)
        ...

    Returns:
        설계 결과 Dict ('forward','reverse','full_sequence', ...)
    """
```

### 타입 힌트 / 에러 처리
```python
from typing import List, Dict, Optional, Tuple

def get_part_by_id(schema_name: str, table_name: str, part_id: str) -> Optional[Dseqrecord]: ...

# 구체적인 예외 + 의미 있는 메시지
raise ValueError(
    f"Overhang 불일치: {prev.name} 끝({prev_right})과 "
    f"{curr.name} 시작({curr_left})이 상보적이지 않습니다."
)
```

---

## 객체지향 설계 원칙 (SOLID) — 현행 코드 예시

### Single Responsibility
조립(`GoldenGateModule`), 부품 조회(`DatabaseImport`), 파일 I/O(`DNAFileReader/Writer`), 검증(`DNAValidator`)이 각각 분리되어 있다.

### Open/Closed
```python
# assembly/factory.py — 새 조립 방법 추가 시 이 dict만 확장
cls._assemblers = {
    'gibson': GibsonCloning,
    'golden_gate': GoldenGateModule,
    # 'slic': SLICAssembly,   # 기존 코드 수정 없이 확장
}
```

### Dependency Inversion
```python
# 워크플로우는 구체 클래스가 아닌 AssemblyBase 인터페이스에 의존
assembler = AssemblyFactory.create_assembler(method)   # 런타임 교체 가능
results = assembler.assemble(fragments, assembly_type="linear")
```

### 새 Assembly 메서드 추가 절차
1. `AssemblyBase` 상속, `prepare_fragments()`/`assemble()`/`validate()` 구현
2. `assembly/factory.py`의 `_assemblers`에 등록
3. 필요한 효소/파라미터는 `basic/`에 공통화

---

## 생명과학 소프트웨어 특수 요구사항

### 서열 검증
모든 외부 서열 입력은 `validate_dna_sequence()`로 검증한다 (유효 염기 `ATCGN`, 빈 값/길이 체크).

### 표준 포맷
GenBank 저장 시 `molecule_type` annotation 필요 (pydna `Dseqrecord.write()`가 처리). FASTA/GenBank I/O는 `utils/io.py`로 일원화.

### 재현성
- 파일 순서는 항상 `sorted()`
- DB 조회 쿼리도 `ORDER BY`로 순서 고정 (`part_order`, `ID`)
- 분석 파라미터는 결과와 함께 기록

---

## 테스트 & 품질 관리

현재 정식 pytest suite는 없고, `tests/`의 실행 스크립트가 DB 기반 통합 테스트 역할을 한다.
- `tests/igem/test_single.py` — 단일 케이스 디버깅용
- `tests/igem/test_all.py` — 전수 회귀 테스트 (성공/추출실패/조립실패/검증실패 집계, 로그 저장)
- `tests/kribb/test_kribb.py` — KRIBB 케이스

### 커밋 전 체크리스트
- [ ] 디버그 `print`/임시 코드 제거
- [ ] 사용하지 않는 import 제거
- [ ] Docstring 작성 (한글)
- [ ] 관련 Test_ID로 회귀 테스트 통과 확인

---

## Git Commit Guidelines

모든 commit 메시지는 다음 3단계 형식을 따른다.

### 1단계: 변경 유형별 요약 (맨 위)
```
ADD: 추가된 기능/파일 요약
MOD: 수정된 코드 요약
FIX: 수정된 버그/에러 요약
```
- 여러 유형은 줄바꿈으로 구분, 각 요약은 한 줄, 파일명 포함 권장

### 2단계: 상세 내용 (# ## 계층)
```markdown
# 주요 변경사항
## 신규 추가 (ADD)
- 파일명: 설명
## 수정 사항 (MOD)
- 함수명: 개선 사항
## 버그 수정 (FIX)
- 에러 설명: 해결 방법
## 검증 결과
- 테스트 항목 / 실행 결과
```

### 3단계: 파일별 라인 변경 분류
```markdown
### 파일명.py
**ADD:**
- Line 100-120: 새로운 함수 추가
**MOD:**
- Line 50-80: 기존 함수 리팩토링
**FIX:**
- Line 150: import 경로 수정
```

### 커밋 전략
- 하나의 커밋은 하나의 논리적 변경 (ADD/MOD/FIX 구분되게 분리 권장)
- 첫 줄 요약만 읽어도 변경 내용을 파악할 수 있게 작성

### ⚠️ 금지 사항
- 커밋 메시지에 다음 footer를 포함하지 않는다:
  ```
  🤖 Generated with [Claude Code](https://claude.com/claude-code)
  Co-Authored-By: Claude <noreply@anthropic.com>
  ```

---

## 참고
- 레거시(Excel 기반 v1/v2) 설명: `archive/CLAUDE.md`
- 개발/디버깅 기록: `develop_log/nyb/`
- 외부 라이브러리: pydna(https://github.com/BjornFJohansson/pydna), Biopython, primer3-py
