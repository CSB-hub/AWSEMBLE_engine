# AWSEMBLE_GUIDE — awsemble 프로젝트 단일 가이드

> **이 문서 = awsemble의 단일 진실원천(SSOT).** 코드↔생물학 매핑 + 작업 진행상황 + 로드맵을
> 한 곳에서 본다. `develop_log`에 새 작업이 기록되면 §6 진행상황·로드맵도 같은 커밋에서 갱신한다.

- **버전 / 최종 갱신**: 2026-06-25 (P7 + golden→gibson 전수 + 저장소 정리 반영)
- **현재 상태**: part→module(Golden Gate) / module→vector(Gibson) 루트 동작. igem 전수 조립 ok 28,178. 미해결은 §6.3 로드맵.
- **2026-11-06부터 `Awsemble_csbtools`가 엔진의 메인 저장소**다 — `awsemble_dev`·Galaxy의 사본에서 동기화하지 않는다(그쪽에 없는 수정이 여기 있다).
- **기준**: 현행 `awsemble/src/` 코드 (DB 기반). 본 문서가 구가이드(archive)와 충돌하면 본 문서를 따른다.
- **구버전 스냅샷**: `archive/docs/`(awemble-완벽-가이드 v1, DBTL_WORKFLOW_GUIDE)

---

## 1. 개요

**awsemble**는 조합 대사 경로 설계를 위한 DNA 조립 자동화 프레임워크다. 표준 부품(Part)을 골라 발현 카세트(Module), 나아가 대사 경로(Pathway)로 조립하는 과정을 자동화한다.

v1과의 가장 큰 차이는 **부품 입력이 Excel → 데이터베이스(MySQL)로 전환**된 점이다. 부품 서열은 `Bioparts` 스키마에서, 오버행 후보군도 DB에서 직접 가져온다.

### 생물학적 계층

| 계층 | 생물학적 의미 | 예시 |
|---|---|---|
| **Part** | 최소 유전 요소 | 프로모터, RBS, CDS, 터미네이터 |
| **Module** | Part를 이어 만든 발현 카세트 | promoter + RBS + CDS + terminator |
| **Multi-module / Pathway** | Module을 이어 만든 대사 경로 전체 | 여러 효소 발현 카세트의 집합 |

---

## 2. 두 가지 흐름 — 유저 흐름 vs 개발 검증 흐름

> ⚠️ **반드시 구분할 것.** "DB 정답 검증"은 awsemble를 *개발·검증*하기 위한 회귀테스트 단계다.
> 실제 유저가 쓸 때는 정답이라는 개념이 없고, **조립 결과(output)가 나오면 거기서 끝**이다.

### 2-A. 기본(유저) 흐름 — 제품의 본질

```
[부품 목록 + 조립 순서]            ← 유저 입력 (무엇을 어떤 순서로 조립할지)
        │
        ▼
[부품 서열 조회]                   Bioparts 스키마에서 서열 로드
   DatabaseImport.get_part_by_id()
        │
        ▼
[오버행 · 프라이머 자동 설계]      효소 결정 + overhang 할당 + 프라이머 구조 생성
   GoldenGateModule.prepare_fragments()
   └─ PrimerDesign.design_goldengate_primers()
        │
        ▼
[조립 시뮬레이션]                  절단 → 조각 선택 → ligation
   BiopartWorkflow.run_step3_module_assembly()
   └─ GoldenGateModule.extract_insert_fragment() / assemble()
        │
        ▼
[output: GenBank 저장]  ★유저 흐름은 여기서 종료★
   DNAFileWriter.write_genbank()  →  modules/inserts/<id>.gb
```

### 2-B. 개발 검증 흐름 — 기본 흐름 + α (회귀테스트 전용)

기본 흐름 위에 **"기댓값과의 대조"** 단계가 얹힌다. 이 부분은 `tests/`에서만 쓰이고, 유저 배포 경로에는 포함되지 않는다.

```
(추가 입력) Test 스키마 활용
   ├─ get_ordered_parts_by_test_id()      # 테스트 케이스가 정의한 부품 목록·순서
   └─ get_construct_sequence_by_test_id() # ★정답 서열(test_construct.sequence)★

(추가 출력) 조립 결과 ↔ 정답 대조
   ├─ DNAValidator.check_parts_correct()  # 부품이 순서대로 들어갔는지 위치 추적
   └─ compare_sequences() (pairwise2)     # 정답과 100% 일치 여부
```

### 두 흐름에서 Test 스키마의 위치

| DB 요소 | 유저 흐름에서 | 개발 검증 흐름에서 |
|---|---|---|
| `Bioparts.<part_type>.sequence` | ✅ 부품 서열 소스 (실사용에도 유효) | 동일 |
| `Bioparts.overhang` | ✅ overhang 후보 풀 | 동일 |
| `Test.test_part` (부품 목록·순서) | ❌ 유저 입력으로 대체됨 | ✅ 테스트 케이스 정의 |
| `Test.test_construct.sequence` (정답) | ❌ 존재하지 않는 개념 | ✅ 기댓값 |

> 정리: **Bioparts = 부품 라이브러리(제품·테스트 공용)**, **Test = 테스트 케이스 + 정답(개발 전용)**.

---

## 3. 전체 아키텍처

```
basic/                  ← 생물학 기본 연산 (PCR 성패를 좌우하는 계산)
  ├─ sequence_analysis.py   DNA 검증, Tm/GC, 제한효소 탐색
  ├─ secondary_structure.py 2차 구조(hairpin)
  └─ primer_design.py       PrimerDesign / ValidatePrimer

assembly/               ← 조립 전략 (Strategy + Factory)
  ├─ base.py    AssemblyBase (ABC)
  ├─ factory.py AssemblyFactory  ("gibson" | "golden_gate")
  ├─ common.py  feature 중복 제거 / 위치 조정
  ├─ golden_gate.py GoldenGateModule (Type IIS)
  └─ gibson.py      GibsonCloning (상동 overlap)

workflow/               ← 단계 오케스트레이션
  ├─ base.py    WorkflowBase (디렉토리/저장/리포트)
  ├─ biopart_assemble.py    BiopartWorkflow (Module 조립 핵심)
  ├─ classical_assemble.py  ClassicalWorkflow (insert+vector 클로닝)
  ├─ spacer_select.py       SpacerSelector (P7 spacer 2단계 선정)
  └─ builder.py WorkflowBuilder (진입점)

utils/
  └─ io.py  DatabaseImport / DNAFileReader / DNAFileWriter / DNAValidator
```

의존 방향: `workflow → assembly → basic`, 그리고 `workflow·assembly → utils`. 상위 계층은 하위 계층의 **추상 인터페이스**(`AssemblyBase`)에만 의존한다(의존성 역전).

---

## 4. 계층별 코드 ↔ 생물학 매핑

### 4.1 `basic/` — PCR 성패를 좌우하는 기본 연산

#### sequence_analysis.py
| 함수 | 생물학적 의미 |
|---|---|
| `validate_dna_sequence`, `validate_cds` | 입력 서열/ORF가 유효한지 (ATCGN, 시작·종결코돈) |
| `calculate_tm` | 프라이머-템플릿이 결합·분리되는 융해 온도 → PCR annealing 온도 결정 |
| `calculate_gc_content` | GC 함량(40–60% 권장) → 결합 안정성 |
| `find_restriction_sites` | Type IIS(BsaI/BbsI/BsmBI/Esp3I/PaqCI/SapI) + 범용 효소 마스터 검색. 근거: [restriction_enzyme_selection.md](reference/restriction_enzyme_selection.md) |
| `check_restriction_site_in_part_sequence` | "부품 내부에 조립용 효소 site가 있으면 나쁜 것"으로 해석하는 검증 래퍼 → Golden Gate 효소 자동 전환의 입력 |
| `find_repeats`, `get_sequence_stats` | 반복서열·통계 |

#### secondary_structure.py
- `analyze_single_strand_structure` — 프라이머가 스스로 접혀(hairpin) 결합력을 떨어뜨리는지 ΔG로 평가.

#### primer_design.py — `PrimerDesign` / `ValidatePrimer`
프라이머 설계는 "이 프라이머로 PCR이 실제로 될까?"를 코드화한 품질 게이트다. **primer3-py(`import primer3`) 의존** (requirements: `primer3-py==2.2.0`).

| 메서드 | 역할 / 생물학 |
|---|---|
| `design_basic_primers` | primer3로 순수 증폭 프라이머 설계 (Tm/GC/dimer/hairpin 자동 고려) |
| `design_goldengate_primers` | 양 끝에 `[효소 인식서열]+[spacer]+[overhang]`을 붙여, 절단 시 정확한 4bp overhang이 노출되도록 설계. `start_structure`/`end_structure`/`full_sequence` 반환 |
| `_extract_spacer_from_vector` | 벡터의 효소 site 뒤 spacer 길이 추출 (BsaI 1bp / BbsI 2bp) |
| `ValidatePrimer.validate_primer` | Tm·GC·self-dimer·hairpin·길이 종합 검증 → `passed`/`violations` |
| `check_primer_dimer` | 프라이머 간 결합(ΔG) → 반응 방해 여부 |
| `check_gc_clamp` | 3' 말단 GC → 결합 안정성 |
| `check_primer_specificity` | off-target binding → 비특이 증폭 방지 |

> **효소별 파라미터 차이의 생물학적 이유**: BsaI는 (1/5), BbsI는 (2/6) 절단 특성 → spacer가 각각 1bp/2bp. `enzyme_sites` 딕셔너리에 정의.

### 4.2 `assembly/` — 두 가지 조립 메커니즘

`AssemblyBase`(ABC)가 `prepare_fragments` / `assemble` / `validate` 세 메서드를 강제하고, `AssemblyFactory.create_assembler("gibson"|"golden_gate")`가 문자열로 구현체를 선택한다(개방-폐쇄 원칙).

#### Golden Gate (golden_gate.py) — "제한효소가 잘라 붙이는 방식"
**원리**: BsaI(`GGTCTC`)·BbsI(`GAAGAC`) 같은 Type IIS 효소는 인식 서열 *바깥*을 잘라, 인식 서열은 떨어져 나가고 4bp **점착성 말단(overhang)**만 남긴다. 상보적 overhang끼리만 ligation → 여러 부품을 한 번에 순서대로 조립.

| 메서드 | 생물학적 의미 |
|---|---|
| `prepare_fragments` | 전 부품 스캔 → 효소 결정 → overhang 할당 → 프라이머 설계 |
| `_get_high_fidelity_overhang_set` | DB overhang 풀에서 서로 해밍거리 먼 세트 선택 → mis-ligation 방지 ([overhang_guide.md](reference/overhang_guide.md)) |
| `extract_insert_fragment` | 절단 후 여러 조각 중 원본과 가장 유사한 insert 선택 (`SequenceMatcher(autojunk=False)`) |
| `assemble` | 인접 조각의 4bp 말단이 역상보(RC)인지 확인 후 결합. circular은 양 끝 overhang도 상보적이어야 닫힘 |
| `validate` | 결과물에 BsaI/BbsI site 잔류 여부 (남으면 다음 반응에서 재절단 위험) |
| `extract_vector_overhangs` | 벡터를 잘라 양 끝 overhang 추출 (circular용) |

**효소 자동 전환 (핵심 안전장치)**: 부품 *내부*에 조립용 효소 site가 우연히 있으면 의도치 않은 위치에서 잘린다. `prepare_fragments`가 전 부품을 스캔하여 — 내부 BsaI 발견 시 **조립 전체를 BbsI로 전환**, BsaI·BbsI가 둘 다 있으면 `RuntimeError`(BsmBI 등 제3효소 권고).

#### Gibson (gibson.py) — "상동서열이 들러붙는 방식"
**원리**: 제한효소 없이, exonuclease가 단일가닥을 노출시키고 인접 fragment의 **상동 overlap(15–40bp)**끼리 어닐링 → polymerase+ligase가 봉합.

| 메서드 | 생물학적 의미 |
|---|---|
| `prepare_fragments` | insert 증폭 프라이머에 벡터 말단 서열을 덧붙여 overlap 생성 |
| `amplify_insert` | pydna `pcr()`로 PCR 시뮬레이션 |
| `assemble` | pydna `Assembly(terminal_overlap)`로 overlap 자동 탐색·조립 |
| `_validate_primer_length` | 짧은 insert(<50bp)에서 프라이머가 템플릿보다 길어지는 문제 보정 |

#### common.py
- `remove_duplicate_features` — vector·insert 양쪽의 동일 feature(예: promoter) 중복 제거
- `adjust_feature_positions` — 프라이머로 길이 변동 시 feature 위치 shift

### 4.3 `workflow/` — 단계 오케스트레이션

`BiopartWorkflow`는 Design-Build-Test-Learn 6단계를 메서드로 구현한다. **부품을 Excel에서 읽지 않고, 이미 만들어진 `Dseqrecord` 리스트를 생성자에 주입받는다.**

| 단계 | 메서드 | 현황 |
|---|---|---|
| Step 1 | `run_step1_part_design` | **`pass`** — Part 준비는 워크플로우 앞단(DB조회 + `prepare_fragments`)에서 끝남 |
| Step 2 | `run_step2_part_cloning` | 현재 스킵 |
| **Step 3** | `run_step3_module_assembly` | **핵심 구현부.** `module_df`(id, name) 기반으로 메모리의 `self.parts`에서 부품을 꺼내 조립 → `modules/inserts/<id>.gb` |
| Step 4–6 | `run_step4`~`run_step6` | 골격만 존재, 본문 미구현 |

- `WorkflowBase` — `create_directories` / `save_results` / `generate_report` 등 공통 기능
- `ClassicalWorkflow` — insert+vector 단순 클로닝 (1회 조립)
- `WorkflowBuilder` — `create_classical` / `create_biopart` 진입점

### 4.4 `utils/io.py` — DB·파일·검증

| 클래스 | 역할 |
|---|---|
| `DatabaseImport` | MySQL(Bioparts/Test) 조회. `.env`의 `HOST/PORT/USER/PASSWORD` 사용 |
| `DNAFileReader` | FASTA/GenBank 읽기 (`read_genbank_dir`은 `sorted()`로 순서 고정) |
| `DNAFileWriter` | GenBank/FASTA 쓰기 (유저 흐름의 최종 output) |
| `DNAValidator` | **개발 검증 전용** — 조립 결과 ↔ DB 정답 대조 |

**DatabaseImport 주요 메서드**
| 메서드 | 흐름 | 역할 |
|---|---|---|
| `get_part_by_id(schema, table, id)` | 유저+개발 | Bioparts에서 서열 → Dseqrecord |
| `get_all_overhang_ids` / overhang 조회 | 유저+개발 | overhang 후보 풀 |
| `get_ordered_parts_by_test_id` | 개발 | 테스트 케이스의 부품 목록·순서 |
| `get_construct_sequence_by_test_id` | 개발 | 정답 서열 |
| `get_all_test_ids_from_construct`, `get_test_ids_by_institution` | 개발 | 테스트 대상 목록 |

---

## 5. 표준 실행 흐름 (코드 워크스루)

기준: [tests/igem/test_single.py](../tests/igem/test_single.py). **★ 표시가 개발 검증 흐름(유저 배포 경로엔 없음)**.

```python
db = DatabaseImport()

# [개발] 테스트 케이스에서 부품 목록·순서 로드 (유저 흐름이면 유저 입력으로 대체)
ordered = db.get_ordered_parts_by_test_id(test_id)        # ★입력 소스
parts   = [db.get_part_by_id("Bioparts", ptype, pid) for pid, ptype in ordered]

# [유저+개발] 오버행·프라이머 자동 설계 → PCR 산물(full_sequence)
gg = GoldenGateModule()
design = gg.prepare_fragments(parts, assembly_type="linear")

# [유저+개발] 설계 결과를 Dseqrecord로 (원본/효소 정보를 annotations에 운반)
records = []
for res in design:
    r = Dseqrecord(res['full_sequence'], name=f"{res['part_name']}_designed")
    r.annotations['target_seq'] = res['target_seq']           # 절단 조각 선택용 대조군
    r.annotations['restriction_enzyme'] = res['restriction_enzyme']
    records.append(r)

# [유저+개발] Step3 조립 → output GenBank  ★유저 흐름은 여기서 종료★
wf = BiopartWorkflow(project_dir=PROJECT_DIR, parts=records)
module_df = pd.DataFrame({"id": [test_id]*len(records), "name": [r.name for r in records]})
wf.run_step3_module_assembly(module_df=module_df, method="golden_gate")
#   → modules/inserts/<test_id>.gb

# ★[개발 검증 전용] 정답 대조 — 유저 흐름엔 없음
DNAValidator(db).check_parts_correct(output_gb_path, test_id)   # 부품 순서·위치
expected = db.get_construct_sequence_by_test_id(test_id)         # 정답 서열
mismatches = compare_sequences(assembled_seq, expected)          # pairwise2 100% 일치 확인
```

**annotations로 원본 정보 운반**: Golden Gate에서 절단 후 올바른 insert 조각을 고르려면 원본 서열이 필요하다. `prepare_fragments` → Dseqrecord 변환 시 `target_seq`/`restriction_enzyme`를 annotations에 넣고, `extract_insert_fragment`가 이를 읽어 유사도 비교에 사용한다.

---

## 6. 진행상황 · 로드맵 (단일 진실원천)

> `develop_log`에 새 작업이 기록될 때마다 이 절을 갱신한다. 상세는 각 출처 링크 참조.

### 6.1 구현 완료 (계층)
- `basic/` 서열 분석·프라이머 설계/검증
- `assembly/` Golden Gate(효소 자동 전환, high-fidelity overhang, 절단 조각 선택) + Gibson
- `workflow/` BiopartWorkflow(Module 조립), ClassicalWorkflow, **WorkflowBuilder.build_plasmid / build_multimodule_plasmid**, **SpacerSelector**(P7 2단계 선정)
- `utils/io.py` DB 연동·I/O·검증, **part_intake.py**(유저 입력 레이어)
- **`params.py` 설계 파라미터 단일 지점**(2026-09-04): Tm 조건·프라이머 목표·효소·homology·
  5' 패딩이 모두 여기 모인다. 연구실 회신이 오면 **이 파일만 고친다**.
  overhang·spacer는 필드가 아니라 **효소에서 파생**(Biopython 절단부위 정의)
- **Tm 창 통합 + 프리셋 + 판정 축 파라미터화**(2026-09-23): 설계창과 판정창을
  `[tm_min, tm_max]` **하나로** 합쳤다(구 `target_tm ∓ tm_tolerance` / `qc_tm_min/max`).
  polymerase 프리셋(`PRESETS` + `from_preset()`), 판정 축 `qc_*`(GUI 비노출),
  패딩 서열 자동 검증, `validate()`의 파이프라인 배선
- **`YeastHRCloning`**: 2026-09-23 스텁 → **2026-10-24 구현 완료**. `GibsonCloning`을
  상속해 fragment 준비·프라이머 설계를 공유하고 yeast 고유 3가지만 재정의한다 —
  **arm 축**(`ARM_PARAM = "yeast_homology"`, 기본 45bp / 하드 30~60),
  **접합부 교차 하드 차단**(Gibson은 경고. in vivo라 잘못된 재조합이 회수 후에야 드러난다),
  **벡터 적격성**(`check_yeast_vector` — 효모 복제원점 ARS/CEN·2μ + 선택마커).
  산물은 **기대 산물**이라는 in vivo 전제를 `validate()`가 유지한다.
  벡터 선형화는 **효소 경로로 충분**하다 — DB yeast 셔틀 후보 31개가 전부 단일절단 효소를
  4~13개씩 보유한다. GUI(툴 XML) 노출은 별도 트랙
- **인식서열 탐색 단일화**(2026-10-10): `find_restriction_sites`가 Biopython
  `Restriction.search()`(= **절단 위치**)를 쓰다가 서열 끝쪽 인식서열을 놓치고 있었다.
  인식서열을 양쪽 가닥에서 문자열로 찾도록 고치고, `params.py`의 중복 구현을 제거해
  **`basic/sequence_analysis.py` 한 벌로 합쳤다.** 조립 경로의 효소 자동 전환이 이 함수를
  쓰므로, 운영 Bioparts 11227개 중 35개가 받던 잘못된 효소 판정이 함께 정정됐다
- **`part → module` Gibson 경로 개통**(2026-10-24): `workflow/junction_plan.py` 신규 +
  `GibsonCloning.prepare_parts()`. GG는 부품을 하나씩 준비할 수 있지만 homology 조립은
  접합부가 **인접 두 부품이 공유하는 구간**이라 이웃 관계를 풀어야 한다.
  긴 부품은 맨몸 PCR, **짧은 부품이 homology arm을 운반**하고(올리고 주문 서열에 넣으면
  프라이머가 아니라 서열이라 열역학 상한을 안 받는다), arm보다 짧은 부품은 이웃에 병합한다.
  품질 미달 접합부(실측 26%)에는 **padding을 자동 삽입**한다.
  arm은 `gibson_overlap`을 vector 삽입과 **공유**한다(기본 20bp, 설정 가능)
- **anneal shift 폴백**(2026-10-24): 맨몸 PCR 설계 실패의 원인이 `anchor_ends`로 5' 시작이
  고정돼 **달성 가능한 Tm이 선언 범위를 건너뛰는** 것이었다(한 부품에서 "low tm"과
  "high tm"이 동시에 발생). 시작점을 k nt 안으로 밀고 잃은 염기를 프라이머 tail로 되살려
  **산물을 바꾸지 않으면서** 설계 성공률 65% → **97%**(`MAX_ANNEAL_SHIFT = 16`)
- **하이브리드 열역학**(2026-10-24): primer3의 60bp는 **합성 상한이 아니라 계산 상한**이라,
  긴 조립 프라이머의 hairpin/dimer를 우리만 못 보고 있었다. >60bp는
  ViennaRNA(DNA Mathews 파라미터)로 계산한다 — **스케일이 달라 경고 전용**이고 어떤 게이트도
  움직이지 않는다. hairpin(`_analyze_long_strand`)과 dimer(`analyze_long_duplex`, cofold) 양쪽

### 6.2 작업 타임라인 (develop_log/ljw)
| 단계 | 내용 | 결과 | 출처 |
|---|---|---|---|
| **리팩토링 P0~P6** | 모드 rename(Biopart/Classical), workflow 재정비(고정버킷+manifest+stage+run_plan), 검증 레이어 | 풀회귀 **9377 / 621(제3효소) / 1(데이터)** 유지(회귀 0) | [260611_refactor_before_after.md](../../develop_log/ljw/plasmid_assemble_draft/260611_refactor_before_after.md), [260611_full_regression_results.md](../../develop_log/ljw/plasmid_assemble_draft/260611_full_regression_results.md) |
| **P7 spacer/vector** | spacer 2단계 선정(anneal + 조합 dry-run) + vector 삽입(Gibson) | Awsemble spacer **115** 채택 | [260614_P7_spacer_design.md](../../develop_log/ljw/plasmid_assemble_draft/260614_P7_spacer_design.md), [260614_P7_design.md](../../develop_log/ljw/plasmid_assemble_draft/260614_P7_design.md) |
| **golden→gibson 전수** | part→module(GG) + module→vector(Gibson) 전수 | 30,062건 **ok 28,178** / alt 1,875 / fail 9 | [260624_golden_to_gibson_results.md](../../develop_log/ljw/plasmid_assemble_draft/260624_golden_to_gibson_results.md) |
| **part_intake** | 서열/genbank → 부품 추출 + DB 매칭 | 서열·gbk 입력 동일 construct 검증 | [part_intake_result.md](../../develop_log/ljw/util_function_add/part_intake_result.md) |
| **리팩 잔버그1** | 버그1 gbk feature circular wrap(origin 넘는 부품 feature 소실) 수정 | `join()` 보존 (회귀 0, annotation) | [260625_bug_fix-add_function.md](../../develop_log/ljw/plasmid_assemble_draft/260625_bug_fix-add_function.md) |
| **gbk 산출물 수정** | part_type 유실(파일 왕복 + prep 두 곳) → COMMENT 마커 왕복 + prep 전달. 파생 회귀(Gibson CDS 자동감지가 module 내부 CDS만 증폭)를 `insert_whole`로 해결 | 서열 불변(SHA 동일), feature type만 정정 | `galaxy_dev` repo: `develop_log/ljw/20260903_part_type_roundtrip_and_insert_whole_regression.md` (repo가 달라 링크 대신 경로 표기) |
| **설계 파라미터화 CP0~CP6** | 하드코딩된 Tm·GC·길이·homology·효소를 `params.py`로 모으고 엔진 5계층에 배선. 효소 객체화(BsmBI 추가)·overhang 파생·derived_id/manifest 반영 | 파라미터 미지정 시 **baseline SHA 불변**. 부수 버그 5건 수정 | [260903_design_parameterization_plan.md](../../develop_log/ljw/260903_design_parameterization_plan.md) |
| **5' 패딩** | 인식서열이 프라이머 말단에 붙어 절단 효율이 떨어지던 문제(NEB) — 여분 염기 6bp 부착. `extract_insert_fragment` 슬라이스도 함께 수정 | 절단 fragment 불변(SHA 동일) | 동상 |
| **CP7 circular** | part→module 중간 클로닝(단일+multi). 벡터 backbone이 조립에 참여하도록 정정(overhang RC 어긋남), module_excision stage 신설, spacer를 부품으로 붙여 벡터 1개로 커버 | 실제 벡터 V0001/V0002로 전 구간 통과. linear SHA 불변. 검증 16종 | [260904_cp0_cp7_total_summary.md](../../develop_log/ljw/260904_cp0_cp7_total_summary.md) |
| **엔진 후속 CP-A~F** | 프리셋·상수 파라미터화·패딩 검증·yeast 스텁·**GG 삽입 경로 해결**·중간 벡터를 실사용 벡터 기반으로 | GG 삽입 spacer 통과 **0/115 → 115/115**. baseline SHA 불변 | [260923_engine_followup_plan.md](../../develop_log/ljw/260923_engine_followup_plan.md) |
| **리팩 잔버그2·3** | 버그2 GG spacer nt A 고정 + 버그3 hairpin 게이트 seqfold→primer3 재작성(2축 정책: −9 tier화·graceful degradation) | QVS 회귀 0. 컷오프 wet-lab 확정 대상 | [260708_hairpin_cutoff_gate.md](../../develop_log/ljw/plasmid_assemble_draft/260708_hairpin_cutoff_gate.md) |
| **인식서열 탐색 단일화** | 두 벌이던 인식서열 탐색을 `sequence_analysis.py` 한 벌로. 끝쪽 사이트 누락(BsaI/BsmBI 양 끝 11nt, BbsI 12nt) 수정 — 조립 경로의 효소 자동 전환이 영향받고 있었다 | DB 11227개 중 **35개 효소 판정 정정**, 성공→실패 전환 **0건**. CP 16개 통과, SHA 불변 | [261010_restriction_site_single_source.md](../../develop_log/ljw/261010_restriction_site_single_source.md) |
| **G1~G3 + H1: part→module Gibson** | 접합부 설계기(병합·운반자·padding) 신설 + `prepare_parts` + anneal shift 폴백 + 하이브리드 열역학(>60bp ViennaRNA, 경고 전용) | 파츠 11,120개 전수: 명세 누락 **0**, 합성상한 초과 **0**, padding 발동 27%, 모듈 e2e **83%**. 검증 20종 + baseline SHA 불변 | [261024_gibson_part_module_and_yeast_plan.md](../../develop_log/ljw/261024_gibson_part_module_and_yeast_plan.md) |
| **Y1~Y2: yeast HR 구현** | `ARM_PARAM` 축 분기(`yeast_homology` 45bp) + 접합부 교차 **하드 차단** + 벡터 적격성(`check_yeast_vector`) + in vivo 전제. 부수로 열역학 정책 통일·조립 프라이머 품질 검증 보강 | arm 45 전수: 병합 88%, padding 29%, e2e **80%**, 명세 실패 1.7%. 검증 **22종 + baseline SHA 불변** | 동상 §7-6·7-7 |
| **Y4: 효소 단일 결정 + yeast 삽입 선형화** | `WorkflowBase.assemble`의 분기를 방법 이름 → `ARM_PARAMS` 훅으로 바꿔 yeast HR도 벡터를 선형화하게 하고, Golden Gate 효소를 `select_gg_enzyme`(벡터 2곳 이상 ∩ 부품 0개)이 **한 번만** 결정하도록 통합. 벡터 삽입의 BsaI 하드코딩 제거 | yeast 삽입 **0% → 동작**(8223bp). GG는 BsaI 없는 벡터에서 BbsI 자동 전환. 실측상 놓치던 설계 **32.7%p** 회수. 검증 **23종 + baseline SHA 불변** | `galaxy_dev` repo: `develop_log/nyb/20261025_조립효소_단일결정과_yeast_선형화.md` |
| **Y5: homology 주문 시트·module feature** | `write_primer_csv`에 올리고 컬럼을 실어 설계된 Top/Bottom이 버려지지 않게 하고, 올리고 시트를 module feature 대신 `design`에서 읽게 통일. homology module feature 기입(`annotate_module_parts_from_layout`)을 `JunctionPlan` 배치 기반 **좌표 계산**으로 신설 | 주문 시트가 fragment를 빠짐없이 분담(올리고 0 bytes → 2건). module feature 0 → 9개, 접합부가 설계기 window와 완전 일치. 부수로 sandbox `part_type` 유실·이름 과다 절단 2건 수정. 검증 **24종 + baseline SHA 불변** | `galaxy_dev` repo: `develop_log/nyb/20261025_homology_주문시트와_module_feature.md` |

| **Y6: 양 끝 spacer 주문서·Tm 판정창** | ① vector homology spacer 길이와 접합부 `amin`이 **같은 변수**(`gibson_overlap`)라 `len >= amin`이 항상 참 → 독립 fragment → 올리고로 빠지던 것을, 체인 맨 끝은 바깥 arm이 없다는 사실로 tail 흡수로 돌림 ② `build_multimodule_plasmid`가 `spacer_mode`를 아예 안 보던 것(GUI 기본값은 `vector_homology`인데 멀티모듈만 라이브러리 spacer를 붙였다) 수정 — 양 끝만 벡터 유래, 연결 spacer는 라이브러리 유지 ③ Tm 판정창을 `qc_tm_*`(50~80℃)로 재분리, **설계창을 포함**하는 조건을 `validate()`가 검사 | Gibson only: fragment 3개(올리고 2)+insert 맨몸 프라이머 → **fragment 1개, 올리고 0개**, 프라이머가 벡터 homology 30bp를 들고 주문 가능(51/52nt). 흉터 0(5943bp = 5550+393). 멀티모듈 5개: 양 끝이 벡터 유래로 바뀌어 **흉터 114bp 제거**(12019 → 11905bp) + 삽입 PCR 불필요. tier red → yellow. 검증 **25종 + baseline SHA 불변**. **사후 회귀(2026-11-05 수정)**: 단축 경로가 `return`으로 빠져나가 그 아래 P6 파츠 feature 기입까지 건너뛰어 부품 1개 module의 construct feature가 `source` 하나뿐이었다(부품 4개는 9개) → 벡터맵에서 파츠를 눌러도 자리를 못 잡음. 조립만 건너뛰도록 고치고 `y6_test.py` [4-1]로 고정 | `galaxy_dev` repo: `develop_log/nyb/20261104_양끝_spacer_주문시트_수정.md` |
| **Y7: GC 판정창 분리** | Y6에서 Tm을 분리한 뒤 **같은 모순이 GC 축에 남아 있던 것**을 해결. 설계창 `[20,80]`%(primer3 입력) vs 판정창 `[35,60]`%라 설계 목표대로 만든 프라이머가 떨어졌고, GC 지적이 `anneal_gate`를 타서 **tail 없는 끝단만 red**가 됐다. GC는 길이로 보상되지 않아(템플릿의 성질) 차단창만 넓히면 신호가 사라지므로 **2층**으로 나눴다 — 차단창 `qc_gc_*`(20~80%, 설계창 ⊇ 강제)과 경고창 `qc_gc_ideal_*`(40~60%, `anneal_gate` 미사용). 곁에 있던 소음도 제거(`_SOFT_GC_RANGE` (30,70)→(20,80) — 기본값에서 소프트 경고가 떴다) | 같은 anneal(mCherry rev, GC 28.0%)이 끝단 **red → yellow**, tail 있는 쪽은 yellow 유지 → **운반자에 따른 적/황 갈림 소멸**. GC 밴드: 15%→red / 28%→yellow / 50%→green / 70%→yellow / 90%→red. `DesignParams().validate()`가 빈 리스트. yellow가 늘어나므로 `design_conditions.md`에 **GC 등급 설명 절**을 신설(숫자는 실제 적용값에서 읽음). 35%가 틀린 컷오프라는 건 자체 실측 근거가 있다(접합부 축: GC 35~60 통과 68% vs 20~80 통과 99.3%). **서열은 불변** — 판정 축만 건드렸으므로 기존 MY AWSEMBLE 기록과 그대로 일치. 검증 **26종 + baseline SHA 불변** | `galaxy_dev` repo: `develop_log/nyb/20261105_GC_판정창_분리.md` |
| **Y8: Tm 자 primer3 통일** | 운영에서 Tm을 55~65→60~65로 바꿔도 **완전히 동일한 프라이머**가 나온 건을 추적. 파라미터는 엔진까지 제대로 갔고, 원인은 **같은 `tm_min`이 경로마다 다른 자로 읽힌 것**이었다 — gibson·yeast는 `PRIMER_MIN_TM`으로 primer3에 넘기고(primer3 축), GG `_fit_amplicon`은 `calculate_tm`과 직접 비교(Biopython 축). 시트에 적히던 건 Biopython 쪽. 실측 3,100표본에서 오프셋 **+4.27℃**(σ0.42, GC 무관). `calculate_tm`을 primer3 `calc_tm`으로 통일하고 Biopython 축에 묶였던 `qc_tm_*`(50/80→54/84)·`terminal_tm_min`(55→59)만 +4 이동. 선언창·프리셋은 원래부터 primer3 축이라 무변경 → GUI·클라이언트·파리티 그대로 | **gibson·yeast 프라이머 바이트 동일**(보고 Tm만 +4.2), **GG만 교정**되어 anneal 3~9nt 단축(전에는 primer3 기준 ~69℃로 창 상한 65를 넘는 과하게 긴 프라이머를 만들고 있었다). construct 서열 **전부 불변**. 비-ACGT는 조용한 축 혼합을 막으려 명시적 실패(Bioparts 11,738종 전수 ACGT 확인). 검증 **27종 + baseline SHA 불변**. **후속 2건 동시 수정**: ① `design_basic_primers`의 폴백 경고를 `gibson.py`가 안 꺼내가 Build Log·주문 시트 어디에도 안 뜨던 것(이번 건을 유저가 못 알아챈 직접적 이유) → 주의사항 칸 + Build Log 노출, ② 올리고 폴백이 `if not _oligos:` 하나라 "엔진이 0개로 정함"과 "CSV 못 읽음"을 구분 못 해 **이미 프라이머 tail에 실린 30bp spacer를 올리고로 이중 주문**하던 것 → 조립 방법으로 분기. 회귀 `y8[9]` + `tests/check_oligo_sheet_gating.py` | `galaxy_dev` repo: `develop_log/nyb/20261105_Tm축_primer3_통일.md` |
| **Y9: 설계 범위 근접값 + 이유·길이 제안** | "값을 바꿔도 결과가 안 변한다"의 **직접 원인 수정**. 요청 범위로 primer3 후보가 0이면 `DEFAULT_PARAMS`로 갈아타 재설계했는데, 그건 근접이 아니라 기본값이라 **무엇을 요청하든 같은 프라이머**가 나왔다(60~65/70~75/45~50 전부 59.4℃ — 45~50은 더 *뜨거운* 값). 또 후보 0의 진짜 원인은 **쌍 제약**이었다 — fwd는 되는데 rev가 못 미쳐 쌍이 깨지고 되는 쪽까지 버려졌다(EXPLAIN `LEFT ok 1 / RIGHT ok 0`). → **Tm 창만 열고 `OPT_TM`은 유저 값 유지**(primer3 목적함수가 |Tm−OPT|에 벌점 → 0번이 곧 근접값, 새 탐색 로직 불필요) + EXPLAIN으로 사유 생성 + **실제 탐색으로 필요한 결합부 길이 계산** | 55~65→59.42 / 60~65→**60.86** / 45~50→**48.50**로 **요청마다 달라짐**. 60~65·45~50은 근접값이 오히려 요청 창 **안**으로 들어옴(버려졌던 쪽 회수). 길이 제안은 빈말이 아님 — 제안대로 상한 25→26 적용 시 fwd·rev 둘 다 60.86℃로 창 안, 경고 소멸. 70~75는 36nt로도 불가라고 명시. 창 안에서 되는 실행엔 경고 없음(소음 방지). construct 서열 불변. 검증 **28종 + baseline SHA 불변** | `galaxy_dev` repo: `develop_log/nyb/20261105_설계범위_근접값_제안.md` |

| **Y10: Tm 자 NEB Tm Calculator 재현** | Y8(primer3 통일) 직후, **자가 하나가 됐어도 그 하나가 실험자가 보는 자가 아니라는** 더 깊은 문제가 드러났다. GUI에서 Taq/Phusion/Q5를 고를 수 있는데 프리셋이 건드리는 값이 `primer_conc_nM`뿐이라 **효소를 바꿔도 Tm이 사실상 안 움직였다**(같은 창 55~65에서 Tm 폭 **1.1℃**). 그리고 Tm을 맞추는 **목적**인 annealing 온도(Ta)가 산출물에 아예 없었다. → 기준 구현을 **NEB Tm Calculator 재현**으로 교체: Biopython `Tm_NN` + `DNA_NN3` + NEB 버퍼 조성(`NEB_TM_CONDITIONS`). 핵심 키 — Taq `K50/Tris10/Mg0/dnac1 200/saltcorr6`(NEB는 Taq에 Mg 보정 안 함) · Q5 `K50/Tris10/Mg1.5/500/6` · Phusion `K50/Tris25/Mg1.5/500/saltcorr1`(Phusion만 Schildkraut) · 공통 `Na=0, dnac2=0, dNTPs=0.8`. `dnac2=0`이 결정적(Biopython은 유효농도를 `dnac1−dnac2/2`로 계산해 기본값이면 농도가 절반으로 깎인다). primer3는 **후보 생성으로 역할 축소**(Tm 창 0~100, `NUM_RETURN=30`) 하고 **선택은 우리 자로** 한다 | **polymerase가 설계에 실제로 반영됨** — 같은 창 55~65에서 Tm 폭 1.1℃ → **8.86℃**(Taq 54.5 / Phusion 58.4 / Q5 62.3), Ta 49.5 / 61.8 / 63.3. NEB 공개값 오프라인 재현 **RMSE Taq 0.29 · Phusion 0.62 · Q5 0.35℃**(18~36nt). `primers.csv`에 **`ta` 컬럼 신설**(`NEB_TA_RULES` — Q5 `min+1`/72 · Phusion `0.93·min+7.5`/72 · Taq `min−5`/68; 프리셋 조성과 다르면 추측하지 않고 비운다) + `design_conditions.md`에 재현 조건·오차·Ta 규칙 절 신설. `qc_tm` 54/84→**45/85**, `na_mM` 50→**0**(1가 염은 `Na+K+Tris/2`=55mM 환산이 기준), `primer_conc_nM` 250→**500**. construct 서열 **전부 불변**(baseline SHA 일치). **primer3를 자로 쓸 수 없는 이유를 고정점으로 박았다** — NEB와의 차이가 고정 오프셋이 아니라 서열별 +2.3~+7.1℃(GC 의존)·길이 기울기(18nt +3.4→36nt +2.7)라 `+N` 보정으로 메울 수 없다(`y8_test.py`[2]). **성능 회귀 동시 해결**: `_length_hint`가 길이창을 1nt씩 전수로 훑어 실패 경로마다 primer3 호출 ~20회(cp5 3분+) → 굵은 격자(4nt)+구간 되짚기 & 후보 30→8로 **2~6회**(실패 경로 0.1초). 검증 **29종 + baseline SHA 불변** + galaxy 파리티 4종 | `galaxy_dev` repo: `develop_log/nyb/20261105_Tm축_NEB재현.md` |

| **Y10-b: 프리셋을 NEB 키트 기준으로 확정** | Y10으로 Tm 자는 NEB가 됐지만 **프리셋 Tm 창이 primer3 축 시절 숫자(50~60/55~65/57~67)로 남아** 있었다. NEB 축으로 재설정: **Taq 55~65 · Phusion 60~72 · Q5 60~72**. 그 과정에서 **Tm 창과 길이창이 커플링돼 있다**는 걸 찾았다 — 제품마다 Tm 자가 달라(같은 프라이머가 Taq 54.5 / Q5 62.3℃) 같은 Tm에 닿는 데 필요한 결합부 길이가 다른데, 길이창은 전 제품 공통(18/20/25)이라 **Taq·Phusion이 선언한 창에 못 닿아 상시 경고**가 떴다. 실측으로 병목이 `primer_len_max`가 아니라 **`primer_len_opt`**임을 확인(primer3 목적함수가 |길이−OPT|에 벌점 → 상한을 36까지 올려도 후보 30쌍이 20nt 근처에 몰린다: `opt20/max36 → 경고 1/2` vs `opt24/max30 → 0/2`). → **길이창(opt/max)도 프리셋이 정하게** 했다(Taq 24/30 · Phusion 22/30 · Q5 20/26 — '경고 0이 되는 최소'. 길이는 곧 합성 단가라 넉넉히 둘 이유가 없다) | 세 제품 **모두 경고 0** (4부품 module, AT-rich 접합부 포함). 주문 시트 보강: `primers.csv`에 **`fwd_mer`/`rev_mer`**(주문 단위 — 서열을 눈으로 세게 하면 오발주가 난다. tail 포함 실물 길이), `vector_primers.csv`에 **mer·anneal Tm·`ta`** 신설(`gibson.prepare_fragments`가 설계기의 anneal Tm을 `last_primer_design`으로 흘려줌 — 완성 프라이머에서 tail을 슬라이스로 되짚으면 `_validate_primer_length`가 길이를 건드린 경우 틀린 구간을 재게 된다). GUI: **"Default"(제품 미지정) 옵션 제거 + Taq 기본 선택**(Tm 축이 제품 조성을 따라가므로 미지정은 어느 계산기와도 대조 불가), `Primer parameter **(NEB kit standard)**` 라벨, 옵션에 결합부 길이·Ta 규칙 표기. `cp5_test.py`의 `fwd_len` → `fwd_mer` 동반 수정. 검증 **29종 + baseline SHA 불변** + galaxy 파리티 4종 | `galaxy_dev` repo: `develop_log/nyb/20261105_Tm축_NEB재현.md` §10 |
| **Y11: 판정창 제거 — 경고 기준을 유저 선언창 하나로 + ΔTm(fwd−rev) 신설** | Y10-b 직후 주문 시트 경고가 **행당 2.0건**이었는데, 실측 분류해보니 유저가 실제로 고칠 수 있는 건 78행 154건 중 **6건(4%)**뿐이었다. 나머지는 **유저가 설정한 적 없는 기준** 위반이었다 — GC 이상구간 `qc_gc_ideal_*` 39건(GC는 주형의 성질이라 봐도 할 수 있는 일이 없다) · 결합부 Tm `qc_tm_min` 45℃ 22건 · 그리고 길이 36건은 판정값이 아예 **리터럴 `15`/`30`**이라 Q5(상한 26nt)를 골라도 경고는 30 기준으로 떴다. **2026-11-04부터 같은 버그 세 번째**(Y6 Tm · Y7 GC · Y11 길이)여서 처방 자체를 뒤집었다 — 판정창(`qc_tm_*`/`qc_gc_*`/`qc_gc_ideal_*`)을 **제거**하고 판정을 선언창(`tm_*`/`gc_*`/`primer_len_*`)이 직접 맡게 했다. 판정창이 맡던 물리 경계는 **입력 단계**로 이동 — `TM_INPUT_MIN/MAX` 35~80℃를 **세 곳**(엔진 `validate()` · Galaxy XML `min`/`max` · Vue `tmInputError`)에 같이 둔다. 길이 축만 `tail_len`으로 분기한다 — tail 없는 프라이머는 선언창으로, GG·Gibson·yeast 삽입 프라이머는 tail이 길이·서열 고정이라 결합부가 `tm_diff_max` 수렴의 **종속 변수**(실측 13~43nt로 흩어진다)여서 선언창으로 재면 **유저가 못 없애는 경고**가 된다 → 실제 제약(주문 상한 80nt + `amp_min`)만 본다. 동시에 **ΔTm(fwd−rev)가 축 자체로 없던 것**을 신설(`pair_tm_diff_max` 5℃) — `tm_diff_max`는 *같은* 프라이머의 결합부↔전체 차라 다른 축이다 | 경고 **154건 → 64건**(행당 2.0 → **1.12**): GC 이상구간 39→**0**, 길이 36→**6**(`amp_min` 미달만 — 진짜 짧은 것), 결합부 Tm 22→**0**, ΔTm **0→12건 신설**. ΔTm 실측 69쌍은 중앙 0.89℃인데 **최대 11.68℃**·5℃ 초과 17%였고 최악 케이스가 **GG 벡터 homology 프라이머에 몰렸다**(길이가 overhang·spacer에 묶여 Tm 자유도가 적다). `ta`가 **낮은 쪽 Tm 기준**이라 Δ가 크면 뜨거운 쪽이 자기 Tm보다 ~12℃ 낮은 온도에서 돈다 — 경고에 그 근거를 적는다. `primers.csv`/`vector_primers.csv`에 **`tm_diff` 컬럼**을 `ta` 바로 옆에 신설. **신호를 지운 게 아니라 경고에서 뺀 것** — GC 숫자는 컬럼에 그대로 남고, `amp_min` 미달은 유저 설정값 위반이라 **남겼다**(`y11_test.py`[7]이 고정). 기각한 안 2개도 기록: GG가 `primer_len_max`를 지키게 하면 `tm_diff_max` 초과가 62→**76%**로 늘고, 하한을 41nt로 박으면 41이 상수가 아니어서 경고가 **부호만 바뀌어** 돌아온다. construct 서열 **전부 불변**. 검증 **30종 + baseline SHA 불변** + galaxy 파리티 4종. ⚠️ **클라이언트 빌드 필요**(입력 검증·help 문구가 화면 변경) | `galaxy_dev` repo: `develop_log/nyb/20261105_판정창_선언창_통합.md` |
| **Y12: 쌍 단위 교차 검사 — 접합부 역상보 방향 + one-pot 프라이머 교차 dimer** | "우리가 만든 프라이머끼리 dimer가 안 생기게, homology arm끼리 안 붙게 돼 있는가"를 확인하다 **반쪽만 돌던 검사 하나**와 **아예 없던 검사 하나**를 찾았다. ① 접합부 교차 검사(`JUNCTION_UNIQUE_MAX_BP=15`bp, Gibson 경고 / yeast 하드 차단)는 **이미 있었는데 정방향만** 봤다 — `longest_common_substring_len`이 독스트링에 "역상보는 호출자가 한 번 더 부른다"고 적어뒀는데 그 호출이 없어서 **교차의 절반을 놓쳤다**(homology 조립은 이중가닥끼리 상동성을 찾아 한쪽이 뒤집혀도 재조합이 성립한다 — pydna가 원형 산물을 역상보 회전으로 2개 돌려주는 이유와 같다). ② one-pot GG는 4부품이면 프라이머 8개 = **28쌍**이 같은 튜브에 들어가는데 자기 쌍 4쌍만 검사됐다 → `check_cross_dimers()` 신설. 전체 프라이머로 전수 검사하면 **같은 5' tail 때문에 소음이 되므로**(실측 tail 12bp 단독 ΔG −7.44 / 역상보 −10.82, 2026-10-20에 98% 소음으로 ΔG를 화면에서 뺀 전례) **tail을 뺀 결합부끼리만** 본다 — 결합부는 overhang·spacer 선정으로 실제로 바꿀 수 있다 | 역상보 검출 확인(탐침: 정방향 6bp / 역상보 29bp — 정방향만으론 임계 15 미달이라 **옛 구현이 놓치던 케이스**), `crosstalk_pairs`에 걸린 **방향**까지 기록. 교차 dimer는 **임계 재조정이 핵심이었다** — 처음 `qc_dimer_dg`(−5.0)로 전부 차단하니 **39행 중 24행(62%)이 red**가 돼서 2026-10-20 소음이 재발했다(쌍 35건 ΔG 중앙 −5.85 — −5.0은 "중간 dimer(주의)" 경계이고 PCR이 실패하는 지점이 아니다). **새 숫자를 만들지 않고** 이미 있던 "강한 dimer" 경계 `qc_dg_absolute`(−9.0)로 차단을 옮겨 **62% → 15%**(차단 3쌍/6행), 중간 32쌍은 버리지 않고 기록. 쌍 결과를 담을 면이 없어서 **`reports/cross_checks.csv` 신설**(쌍 1행 — `primers.csv`는 행이 부품 하나라 "A ↔ C"를 담을 자리가 없다) + 차단 쌍은 양쪽 행에 포인터 + **Build Log 승격**(`_log_cross_checks` — 둘 다 에러 없이 조용히 실패하는 종류라서, Y8의 폴백 승격과 같은 처방). **별 건 발견**: `primers.csv`에 **`violations` 컬럼이 없어서** 엔진이 채운 차단 사유가 한 번도 시트에 안 나왔다(tier=red만 뜨고 이유가 없었다 — y7의 `violations` 검사가 공허하게 통과한 이유). 같은 커밋에서 신설. construct 서열 **불변**. 검증 **31종 + baseline SHA 불변** + 파리티 4종 | `galaxy_dev` repo: `develop_log/nyb/20261105_쌍단위_교차검사.md` |
| **Y13: 교차 검사 CSV가 "깨끗함"을 말하게 + 결과창 경고 문장 줄바꿈** | 유저가 실제 실행에서 둘을 지적했다. ① 다운로드한 `cross_checks.csv`에 **헤더만** 있었다 — Y12에서 "빈 파일은 깨끗함과 미검사가 구분이 안 되므로 0건이면 파일을 만들지 않는다"고 정했는데, Galaxy `<data>` 출력은 **파일이 반드시 있어야 해서** helper가 헤더만 있는 파일을 썼다 → **피하려던 모호함이 그대로 유저 손에 도착**했다(그 실행은 homology module 경로였고 결과는 정상 통과였다). 판단을 뒤집어 **검사가 돌았으면 0건이어도 파일을 쓰고**, 검사마다 `<kind>_summary` 행에 **분모(검사한 쌍 수)**를 적게 했다 — "0쌍 걸림"만으로는 검사 범위를 알 수 없고, 분모가 있어야 0이 **결과로** 읽힌다. ② 결과창 경고가 세 문장(무엇이/왜/어떻게 고치나)이 한 줄로 붙어 **실제 조치인 세 번째 문장이 묻혔다** — 엔진에 `\n`을 넣으면 `runLogSummary.ts`가 로그를 **줄 단위로** 읽어 경고를 모으므로 두 번째 줄부터 떨어져 나간다 → **표시 층에서만** 끊는다 | 분모를 담을 자리 2개 신설: `_DimerScan.pairs_checked`(list 상속이라 기존 호출부 무영향) · `JunctionPlan.crosstalk_checked`. 실측 3경로 — GG 3부품 **40쌍/5건** · GG 4부품 **60쌍/6건** · homology 4부품 **1쌍/0건(통과)**. Build Log도 `_summary` 행을 분리하고(안 그러면 깨끗한 실행에서도 ⚠️가 뜬다) 대신 **✅ 줄로 무엇을 몇 쌍 봤는지** 남긴다 — 아무 줄도 안 나오면 빈 CSV와 같은 문제다. 파일이 **없는** 경우는 이제 쌍 단위 검사가 하나도 안 돈 경로뿐이다. 경고 줄바꿈은 `warnLines()` + `white-space: pre-line`, 마침표 **앞 글자가 숫자가 아닐 때만** 끊어 `62.6℃`·`ΔG -9.0`을 보존(실제 경고 2종으로 확인). construct 서열 **불변**. 검증 **31종 + baseline SHA 불변**, y12에 기준 [11]·[12] 추가. ⚠️ **클라이언트 빌드 필요** | `galaxy_dev` repo: `develop_log/nyb/20261106_교차검사_요약행_경고_줄바꿈.md` |
| **Y14: 검사 범위를 반응 단위로 재정렬 — 교차 dimer 제거 · overhang codon 신설 · arm 재설계 래더 · 우회 기록** | 네 건이 **한 뿌리**였다 — "같은 반응에 실제로 무엇이 들어가는가"를 코드가 틀리게 잡고 있었다. ① **교차 dimer 검사 제거**: Y12에서 "`prepare_fragments` 한 번이 곧 한 튜브"라고 보고 부품 n개의 프라이머 2n개를 전수 쌍으로 비교해 차단했는데, 실제로는 **부품을 하나씩 따로 PCR**한다 — 다른 부품의 프라이머와 만날 일이 없다(올리고도 쌍별로 따로 어닐링). 일어나지 않는 일로 설계를 막고 있었다. `_check_cross_dimers`와 `primer_dimer` 행을 걷어내고, 프리미티브 `check_cross_dimers`는 "조립 경로는 부르지 않는다"를 독스트링에 박아 남겼다(여러 쌍을 한 튜브에 모으는 프로토콜이 생기면 그때 맞는 자리). ② **overhang frame/codon 신설**: GG overhang은 산물에 그대로 남는 **삽입 서열**(산물 접합부 = 왼쪽 부품 + overhang + 오른쪽 부품)인데 그게 만드는 codon을 아무도 안 봤다(`validate_cds`는 정의만 있고 **호출부 0곳**). 접합부 양쪽 `part_type`으로 번역 여부를 가려 **CDS↔CDS는 in-frame stop**, **RBS/UTR→CDS는 새 ATG**만 본다 — 문맥 전체에서 TAA/TGA/TAG를 찾으면 임의 서열에서도 거의 항상 걸리므로 **overhang 염기를 쓰는 codon만** 센다. CDS↔CDS에 3의 배수가 아닌 overhang(BsaI 4nt)이 끼면 **서열이 뭐든 프레임시프트**라 교체로 못 고치는 길이 문제 → 설계불가 + SapI(3nt) 안내. 교체는 기존 사이트 회피와 **한 패스로 합쳤다**(따로 돌리면 뒤 패스가 앞 조건을 깬다) — `_avoid_site_creating_overhangs` → `_refine_overhangs`. ③ **arm 기준 교차 + 재설계 래더**: 비교 대상이 `fragment.seq`(arm 포함) 전체라 주로 걸리는 건 "같은 부품을 두 자리에 썼다"(실측 81bp)였는데 **그건 arm 길이로 고칠 수 없다** — 재설계 래더를 붙이면 헛돌다 늘 설계불가로 떨어진다. 그래서 두 신호를 분리: **arm vs (다른 arm [+ in vivo면 다른 몸통·백본])**은 재설계로 고쳐지는 축 → 래더(① 유저 범위 → ② `ARM_HARD_RANGE` 안에서 범위 밖 + note → ③ 설계불가), **몸통 vs 몸통**은 `body`로 재서 **기록만**(Gibson은 exonuclease가 말단만 깎아 몸통끼리는 어닐링 자체가 안 되고, yeast는 arm 축이 따로 돈다). in vitro/in vivo 차이는 "경고 vs 차단" 등급이 아니라 **비교 범위**(`CROSSTALK_SCOPE`)로 옮겼다 — `YeastHRCloning._check_plan` 하드 차단은 삭제. ④ **우회·컷오프 기록 신설**: `build_log.txt`는 **백엔드가** 쓰고 엔진은 `print()`만 하는데 stdout 리다이렉트가 없어서, 효소 자동 전환·overhang 교체·anneal shift·spacer graceful degradation 같은 **"우회했다"는 사실이 케이스 산출물에 한 줄도 안 남았다**(journal로만 나갔다). `workflow/design_events.py`(모듈 전역 기록 — 사건 생성 지점이 5단 깊이로 흩어져 인자로 꿰면 경로 전체의 서명이 바뀐다) + `reports/design_events.csv` + 백엔드가 케이스로 복사하고 build log에 이어 붙인다 | **교차 dimer 차단 소멸** — "일어나지 않는 반응"으로 막히는 설계가 없어졌다(Y12의 차단 3쌍/6행이 애초에 범위 오류였다). `cross_checks.csv`의 축이 `overhang_junction` / `arm_crosstalk` / `body_crosstalk`로 재편되고 **분모 행은 유지**(0건과 미검사를 구분). **Gibson도 차단으로 올라간다** — 같은 부품을 두 자리에 쓰면 두 접합부 arm이 **동일 서열**이 되어 어떤 길이로도 못 피하므로(Gibson에서 arm은 실제로 어닐링하는 구간이다) 경고로 둘 수 있는 축이 아니었다: `y2_test.py`[1]을 그 계약으로 교체. **codon 검사는 지금까지 완전 공백**이었다(frame/codon 검사가 `workflow/`·`assembly/` 어디에도 없었다). 래더가 범위를 벗어나면 **왜 벗어났는지 note에 적는다** — 조용히 다른 길이로 돌지 않는다. `design_events.csv`는 **0건이어도 쓴다**(`clean` 행 — 빈 파일 모호함 재발 방지, Y13과 같은 처방). 산출물에 `output_events` 추가 + Metadata zip 포함. **construct 서열 전부 불변** — arm은 산물에 들어가지 않으므로(`module_seq`는 body의 연결) 길이를 다시 골라도 서열이 안 바뀐다: **baseline SHA 일치** 확인. 검증 **32종**(y2·y12 기대값 갱신). 부수로 **실행 기록 영구 보관**(`awsemble_log/` — 워크플로우 양식 + 파츠 리스트 + 파라미터 + build log, 사람별 디렉토리 + 평면 CSV 2벌; `data/runs`는 게스트 보존기간 뒤 지워지므로 추적이 끊겼다). ⚠️ **클라이언트 빌드 필요**(`runLogSummary.ts`가 새 머리표를 읽는다) | `Awsemble_csbtools` 작업(2026-11-06) — 이 저장소가 메인이다 |
| **Y15: module 층이 없는 조립법에 multi-module이 흘러들던 경로 차단** | `build_multimodule_plasmid`는 인접 module을 **라이브러리 spacer**(S0002…)로 잇는다 — "module을 따로 만들어 one-pot으로 합친다"는 Golden Gate 패턴의 부품이다. homology 경로(`gibson`/`yeast_hr`)는 부품이 한 통에서 접합부끼리 직접 붙으므로 중간 spacer가 들어갈 자리가 없다. 그런데 `module_method`를 **아무도 검사하지 않아서**, yeast_hr로 5-module을 넘기면 경고 한 줄 없이 카세트 사이에 spacer가 박힌 construct가 나왔다(실측 β-carotene 15조각: `['YEplac195_homL','S0002','S0005','S0006','S0009','YEplac195_homR']`, 17211bp). 입구가 셋이었다 — ① 엔진이 안 막았고, ② 백엔드가 config를 그대로 넘겼고, ③ **화면은 Module 탭을 숨기기만 했다**(`modules` state는 남아 전송된다 → Module based로 5개 만들어 두고 Yeast로 바꿔 Run하면 그대로 들어간다. 게다가 Module 2 이후 부품에는 손도 닿지 않는다). 세 층 모두에 가드: 엔진 `ValueError`, 백엔드 `RunInputError`(입력 단계), 프런트는 **순서대로 한 Module로 합쳐 전송**(`runModules`) + 합쳐진다는 경고문 + Module이 여럿이면 탭을 남기고 추가·복제만 감춘다 | 같은 부품을 한 Module(15슬롯)로 돌리면 **16983bp** — 정확히 228bp(4×57) 차이로, 그만큼이 유저가 요청한 적 없는 서열이었다. 더 컸던 부수 효과: module으로 쪼개면 arm 교차 검사가 module 안에서만 돌아 **module당 6쌍 × 5 = 30쌍**만 봤다 — 실제 반응은 한 통이므로 비교해야 할 쌍은 **378쌍**이고, 348쌍이 조용히 빠져 있었다(재실행에서 378쌍 clean + 몸통 교차 `PADH1↔BTS1` 16bp note 1건). construct에 `S000*` 0건 확인. 체크포인트 32종 + 회귀 SHA 불변 | `Awsemble_csbtools` 작업(2026-11-06) |

### 6.3 로드맵 / 미해결

> **남은 작업 전수 목록은 [260904_remaining_work_inventory.md](../../develop_log/ljw/260904_remaining_work_inventory.md)** 참조.
> 아래는 그중 주요 항목만.

- ~~CP7 part→module circular~~ **완료(2026-09-04)**: 중간 클로닝 전 구간 동작(단일 + multi-module).
  `assemble(circular)`을 pydna sticky-end ligation으로 바꾸고 벡터 backbone을 조립에 참여시켰다.
  `module_excision` stage 신설(digest 기본 / PCR 폴백).
  **기본 경로는 spacer를 부품으로 붙이는 방식**이라 기성 GG 벡터(V0001) 하나로 모든 module 위치를
  커버한다 — 위치별 전용 벡터가 필요 없다. 전용 벡터를 쓰려면
  `intermediate_vector_has_spacers=True` + `tools_dev/make_intermediate_vector.py`.
  → **총정리: [260904_cp0_cp7_total_summary.md](../../develop_log/ljw/260904_cp0_cp7_total_summary.md)**
- **galaxy 배선 G1~G5**: 툴 XML 파라미터 섹션·helper argparse·프리셋 저장소·GUI 패널.
  **이게 없으면 랩마다 다른 값을 각자 넣을 수 없다**(엔진은 준비 완료, 입력 통로가 없는 상태).
  ⚠️ 2026-09-23 스키마 변경 반영 필요: `target_tm`/`tm_tolerance` → `tm_min`/`tm_max`.
- **회신 대기**: NEB 효소별 패딩 수치(현재 잠정 6bp), 연구실 설계 파라미터, 7/30 tier 경계.
- **SapI**: overhang 3nt 세트가 DB에 없어 보류(4bp를 자르면 fidelity 보장이 깨짐).
  세트 등록 시 코드 변경 없이 `SUPPORTED_ENZYMES`에 이름 한 줄.
  **실험 검증 이후 기능 업데이트에서 다루기로 결정(2026-09-23)** — `Biopart_database` 쪽 작업이다.
- ~~**GG amplicon Tm**(module→vector 삽입 Golden Gate)~~ **해결(2026-09-23)**:
  원인은 6/16 진단(절대 Tm 미고려)도 `amp_min` 바닥도 아니었고 **`tm_diff_max`가 유일한
  구속 조건**이었다. tail(인식서열+overhang)이 고정 길이라 diff를 먼저 좁히면 Tm이 범위를
  넘어가 spacer 115개가 **전부** 탈락했다. `_fit_amplicon`의 선정 순서를 뒤집어
  **Tm을 선언 범위에 맞추고 diff는 결과로 보고**하게 바꿔 115/115 통과.
  유저가 선언한 범위를 GG 프라이머가 실제로 따라간다(구버전은 유저 입력이 무효였다).
  **조립 단계의 2차 버그도 같이 해결**: DB가 vector를 `circular=False`로 주는데 그대로 잘라
  backbone 한쪽이 blunt가 돼 sticky-end ligation이 전패했다 → `_as_circular` 정규화.
  실측 `build_plasmid(V0001, golden_gate)` 실패(487초) → **성공 5720bp(10초)**.
  ⚠️ **GG 삽입 가능 벡터는 387개 중 79개(20%)뿐이다** — 사이트 0개 157 / 1개 151(좌우
  overhang 동일 → 방향성 클로닝 불가) / 2개 이상 79. 사이트 1개 벡터는 가드로 즉시 차단한다.
  → 원 분석: [260616_P7_미해결_GG_amplicon_Tm.md](../../develop_log/ljw/plasmid_assemble_draft/260616_P7_미해결_GG_amplicon_Tm.md)
- **part→module gibson 통로**: 여전히 미검증(보류). part가 GG overhang 전제로 준비돼
  homology overlap이 없어 전용 prep 설계가 선행한다.
- **중간 벡터**: `make_intermediate_vector.py`가 **실사용 벡터**를 받아 유일 절단 부위(MCS 추정)에
  카세트를 삽입한다(2026-09-23). DB 후보 **45개**(BsaI·BsmBI 둘 다 없고 유일 사이트 보유).
  `--list-candidates`로 목록 확인. 합성 최소 backbone은 구조 검토용으로만 남겼다.
- ~~**part→module Gibson 경로 미검증**~~ **완료(2026-10-24)**: 전용 prep
  (`workflow/junction_plan.py` + `GibsonCloning.prepare_parts`)을 만들어 개통했다.
  실제 파츠 모듈 end-to-end 83%, 검증 20종 + baseline SHA 불변.
  남은 실패는 Tm 선언 범위(연구실 회신 대기)에 묶여 있다.
- ~~**yeast HR 미구현**~~ **완료(2026-10-24)**: 확정 4가지(arm 45bp / 효소 선형화 /
  part→module + module→vector 두 단계 / 검증 기준)를 정하고 구현했다.
  arm 45bp 전수조사: 병합 88% · padding 29% · prep e2e 80% ·
  명세 생성 실패 1.7%(부품 전체가 arm보다 짧아 합쳐지는 경우).
  남은 것은 **galaxy 배선**뿐: [261024_gibson_part_module_and_yeast_plan.md](../../develop_log/ljw/261024_gibson_part_module_and_yeast_plan.md) §7-8
- **열역학 하드 차단 → 경고**(2026-10-24): `arm + primer_len_min > 60`이 하드 에러였는데,
  H1이 >60bp를 계산 가능하게 만들어 **경고로 내렸다**(차단 게이트가 사라지는 건 사실이라
  반드시 알린다). yeast arm 45bp가 이 구간에 필연적으로 들어가므로 정책을 통일한 것이다.
- **조립 프라이머 품질 검증 보강**(2026-10-24): `design_basic_primers`는 anneal만 보는데
  tail이 붙은 **전체 프라이머**를 아무도 검증하지 않고 있었다. Gibson/yeast prep이
  tier·경고·판정 엔진을 산출물에 싣는다(GG 경로와 동등).
- **접합부 설계 재설계**(2026-10-25): 전제 두 개가 틀려서 고쳤다.
  ① **PCR 템플릿이 존재하지 않는 분자였다** — 짧은 부품을 이웃 body에 병합해 `lp+sp`를
  템플릿으로 삼았다(모듈의 63~88%). 폐기하고 **프라이머 tail 흡수 / 양쪽 arm을 드는
  독립 올리고** 2갈래로. 불변식: *PCR 템플릿은 실제로 존재하는 분자 하나*.
  ② **padding이 RBS–CDS 경계를 깨뜨렸다** — 품질 미달 접합부면 어디든 들어갔는데 전형
  모듈에서 그 자리가 RBS–CDS다. 조립은 되고 발현이 죽는 실패라 **기능 자체를 삭제**.
  대신 접합부 길이를 **범위**로 주고(Gibson 30~45 / yeast 40~60) 품질을 통과하는 가장
  짧은 길이를 자리마다 고른다. 길이는 공유 구간이라 **최종 서열을 바꾸지 않는다**.
  ③ **품질 게이트가 프라이머 축을 빌려 썼다** — hairpin ΔG≤−3 비율이 창 20bp에서 1%,
  60bp에서 20%로 길이에 끌려간다. 차단은 GC(`gc_min/gc_max` 20~80) 하나로 내리고
  2차구조는 경고로. 실측 통과율 68% → **범위 탐색 시 100%**.
- **2차구조/accessibility 게이트** (코드 완료, 정리: [260708_hairpin_cutoff_gate.md](../../develop_log/ljw/plasmid_assemble_draft/260708_hairpin_cutoff_gate.md)): hairpin/dimer 2축 게이트 + template accessibility(ViennaRNA 미페어확률, GG/Gibson 말단 binding 검사, warn-only) + 3-tier `tier` 필드까지 구현·QVS 회귀 0. **blocked(외부 의존)**: 컷오프/임계·온도 **wet-lab 확정**(−2/−3/−5/−9, accessibility 0.5/0.25, 온도) / **캘리브레이션 훅**(리퀴드핸들러 연동 후).
- **배포 전 최종 내부로직 버그 스윕**(다음 세션): basic/ 파라미터 전반(Tm/GC/ΔG/accessibility) + 조립 로직을 실험자 인계 전 마지막으로 점검. 이후 입력데이터 처리 기능(assemble와 별개) 추가.
- **KRIBB id 불일치**: `test_part`(무prefix `P0019`) vs `Bioparts`(K-prefix `KP0019`) → 추후 DB 정리(P→KP/R→KR/C→KC/T→KT/S→KS)
- **igem-single 진짜 fail 4건**(특정 벡터에서만 = vector 삽입 단계) — 후속 분석
- **part_intake 후속**: `new` 부품 임시 DB 등록 + galaxy 연동
- 과거 알려진 케이스: `BBa_K561005`·`BBa_M31790`(RE+overhang 겹침), `BBa_B0034`(짧은 부품 cut), `BBa_J22203`(overhang 재현성), `BBa_K122016`(부품 `BBa_K122010`=1bp 데이터 불량 → error 1건)

### 6.4 파츠 생성 자산 (별도 repo)
DB 자산(spacer/overhang) 생성 코드는 **`Biopart_database/Part_generation/`** 에 있다(awsemble repo 밖).
- `awsemble_codes.py` — awsemble 엔진 어댑터(facade). docker `/app/src` 마운트로 src 자동 연동.
- `gen_spacers.py`(spacer 생성·검증), `spacer_anneal_filter.py`(적합 필터, DELETE 대상 산출), `overhang_automation.py`(overhang fidelity 세트 생성).
- 실행: `docker run --network host -v $AWSEMBLE/awsemble/src:/app/src:ro -v $BIOPART_DB/Part_generation:/app/part_gen awsemble:dev python /app/part_gen/gen_spacers.py`

---

## 7. 문서 맵 — 어떤 내용은 어디를 볼까

| 주제 | 위치 |
|---|---|
| 코드↔생물학 매핑 + 진행상황 (본 문서, **SSOT**) | `docs/AWSEMBLE_GUIDE.md` |
| overhang 설계 원칙·fidelity 알고리즘 | `docs/reference/overhang_guide.md` (+ `standard_overhangs.csv`) |
| **조립 로직 전체 (Workflow logic 4종·접합부 설계·파라미터 모델)** | **`docs/reference/assembly_logic.md`** |
| 제한효소 목록 선택 근거 | `docs/reference/restriction_enzyme_selection.md` |
| 파츠(spacer/overhang) 생성 코드 | `Biopart_database/Part_generation/` (§6.4) |
| 코딩 표준·DB 스키마·커밋 규칙 | `awsemble/CLAUDE.md` |
| 작업 상세 로그 | `develop_log/ljw/`, `develop_log/nyb/` |
| 구버전 가이드 스냅샷 (v1, DBTL) | `archive/docs/` |
| legacy 설계 문서 / Excel v1·v2 | `archive/docs/development/`, `archive/CLAUDE.md` |

### `src/` 내부 plan 문서 정리 (완료, 2026-06-24)
작성 시점이 오래되어 현행과 어긋난 legacy 설계 문서들. **전부 `archive/docs/development/`로 이동 완료** (src 내부 .md 0개, docstring 참조도 archive 경로로 갱신).

| 파일 (현 위치) | 비고 |
|---|---|
| `archive/docs/development/PLAN.md` | 초기 계획(BsaI만 가정), 성공기준 대부분 달성 |
| `archive/docs/development/PLAN_v2.md` | 확정 설계, 거의 구현 완료 |
| `archive/docs/development/WORKFLOW_ANALYSIS.md` | `dbtl.py 미구현` 시점 분석 → 현재 Step3 구현됨 |
| `archive/docs/development/dbtl_plan.md` | Excel 기반 Step3–6 계획 (현재 DB 기반) |
| `archive/docs/development/plan.md` | Excel/`excel_to_gbk` 가정 (레거시) |
