# AWSEMBLE_engine

플라스미드 조립 설계 엔진. Golden Gate / Gibson / Yeast HR 조립을 설계하고
프라이머·올리고를 뽑는 순수 파이썬 계층이다. **웹도 GUI도 없다** — 두 제품이
이걸 `sys.path`에 올려 그대로 import한다.

```
AWSEMBLE_engine  ──┬─→ AWSEMBLE (csbtools.dev/awsemble/)   ← submodule
                   └─→ Galaxy awsemble 툴 (Spark 8080)      ← tools/awsemble/src로 복사
```

## 구조

| 경로 | 내용 |
|---|---|
| `src/basic/` | 서열 프리미티브 — Tm·GC, 프라이머 설계, 2차구조, 구조 리포트 |
| `src/assembly/` | 조립 전략 (Strategy + Factory) — `golden_gate` · `gibson` · `yeast_hr` |
| `src/workflow/` | 오케스트레이션 — 접합부 계획, spacer 선정, 벡터 상동성, 빌더 |
| `src/utils/` | DB 접속(Bioparts MySQL), GenBank/FASTA I/O, 검증 |
| `src/params.py` | `DesignParams` 한 벌 + `PRESETS`(NEB 키트 기준) |
| `docs/` | `AWSEMBLE_GUIDE.md`가 프로젝트 현황 단일 진실원천(SSOT) |
| `tools_dev/` | 체크포인트 테스트 스크립트 (`params_cp_tests/run_all.sh`) |

진입점은 `src/`다. `PYTHONPATH`에 `src`를 넣고 `from workflow.builder import ...` 식으로 쓴다.

## 이 저장소에 **없는** 것

- **`.env`** — Bioparts MySQL 접속 정보(HOST/PORT/USER/PASSWORD). 커밋하지 않는다.
  쓰는 쪽에서 `awsemble_engine/.env`로 직접 놓거나 환경변수로 넘긴다.
- **`tests/`** — 11G짜리 실행 산출물. 원본은 `BioCAD/Awsemble/awsemble_dev`에 있고
  `Backup/backup.sh`의 `ARTIFACTS`로 보관된다. `CLAUDE.md`가 참조하는
  `tests/igem/test_all.py` 등은 거기서 돌린다.
- **`develop_log/`** — 작업 기록 3.8G. 같은 위치에 있다.

## 상류 — 단일 진실원천은 awsemble_dev다

엔진 코드를 **여기서 처음 쓰지 않는다.** `BioCAD/Awsemble/awsemble_dev/BioCAD_private`
(`BioCAD_private` 저장소의 `dev` 브랜치)에서 개발하고, 그 `awsemble/` 트리를 여기로
동기화한다. 그래야 Galaxy 쪽 사본과 축이 어긋나지 않는다.

동기화할 때는 **커밋 메시지에 상류 커밋 해시를 적는다** — 이 저장소는 히스토리를
새로 시작했으므로 git 조상 관계로는 상류와 대조할 수 없다. 해시가 유일한 단서다.

```
MOD: 엔진 동기화 — synced from awsemble_dev <sha>
```

## 주의 — 소비자와 맞춰야 하는 값

`params.py`의 `PRESETS`(Tm 창·길이창·버퍼 조건)는 **GUI가 같은 값을 거울처럼 들고 있다**
(AWSEMBLE의 `awsemble_frontend/src/lib/workflowLogic.ts`, Galaxy의 `WorkflowPanel.vue`).
한쪽만 바꾸면 화면이 거짓말을 한다 — 고른 조건과 실제로 돌아간 조건이 달라진다.
`PRESETS`를 손대면 소비자 쪽도 같이 고친다.

## 개발 규칙

`CLAUDE.md` 참고. 요점 둘:

- **같은 일을 하는 코드를 두 벌 만들지 않는다.** 2026-10-10에 제한효소 탐색이 두 벌이라
  Bioparts 11227개 중 35개가 효소 판정을 잘못 받은 사고가 있었다.
- `develop_log`에 기록을 남기면 **같은 커밋에서** `docs/AWSEMBLE_GUIDE.md` §6도 갱신한다.
