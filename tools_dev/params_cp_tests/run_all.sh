#!/bin/bash
# 설계 파라미터화(CP0~CP6) 검증 일괄 실행.
#   각 CP의 수용 기준 + "파라미터 미지정 시 서열 불변" 회귀를 함께 확인한다.
# 사용: bash tools_dev/params_cp_tests/run_all.sh
PY=${PY:-/home/csb/awsemble_env/bin/python}
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE/../.." || exit 1     # awsemble 프로젝트 루트
fail=0
#   g1·g3: part→module homology 경로(접합부 설계기·prep), h1: 하이브리드 열역학,
#   y1~y2: yeast HR(arm 축 분기 / 교차 차단·벡터 적격성·in vivo 전제),
#   y3: arm 축이 삽입 방법을 따라가는가(vector_homology spacer 길이 어긋남 회귀),
#   y4: 조립 효소 단일 결정(벡터∩부품) + homology 삽입의 벡터 선형화 회귀,
#   y5: homology 경로의 주문 시트(올리고 컬럼)와 module feature 기입 회귀.
#   (2026-10-24 추가 — develop_log/ljw/261024_gibson_part_module_and_yeast_plan.md)
#   y6: 양 끝 spacer의 tail 흡수 + 멀티모듈 spacer_mode + Tm 판정창 분리 회귀.
#   (2026-11-04 추가 — develop_log/nyb/20261104_양끝_spacer_주문시트_수정.md)
#   y7: GC 판정창 분리(차단창 ⊇ 설계창) + 이상 구간 경고층 + 적/황 비대칭 제거 회귀.
#   (2026-11-05 추가 — develop_log/nyb/20261105_GC_판정창_분리.md)
#   y8: Tm 자를 primer3로 통일(GG는 Biopython, gibson은 primer3로 갈려 있었다) 회귀.
#   (2026-11-05 추가 — develop_log/nyb/20261105_Tm축_primer3_통일.md)
#   y9: 설계 범위를 못 맞추면 기본값이 아니라 **근접값** + 이유·길이 제안 회귀.
#   (2026-11-05 추가 — develop_log/nyb/20261105_설계범위_근접값_제안.md)
#   y10: Tm 축을 NEB Tm Calculator 재현으로 통일 — polymerase 특성이 설계에 실제로
#        반영되는가 + NEB 공개값 오프라인 재현 회귀.
#   (2026-11-05 추가 — develop_log/nyb/20261105_Tm축_NEB재현.md)
#   y11: 판정 기준을 **유저 선언창 하나**로 합침(판정창 qc_tm_*/qc_gc_* 제거, 물리 경계는
#        입력 제한 TM_INPUT_MIN/MAX로) + 길이 판정을 tail_len으로 분기 + ΔTm(fwd−rev) 신설.
#   (2026-11-05 추가 — develop_log/nyb/20261105_판정창_선언창_통합.md)
#   y12: 쌍 단위 교차 검사 — 접합부 교차에 **역상보 방향** 추가(정방향만 봐서 절반을
#        놓치고 있었다) + one-pot 프라이머 **결합부끼리** 교차 dimer 신설(차단은
#        qc_dg_absolute, 기록은 qc_dimer_dg 2등급) + reports/cross_checks.csv.
#   (2026-11-05 추가 — develop_log/nyb/20261105_쌍단위_교차검사.md)
for t in cp0 cp1 cp2 cp3 cp5 cp6 pad cp7 cp7_excise cp7_e2e vec_verify real_vector one_vector mm_linear mm_circular vector_homology g1 g3 h1 y1 y2 y3 y4 y5 y6 y7 y8 y9 y10 y11 y12; do
  printf "%-6s " "$t"
  if timeout 1800 "$PY" "$HERE/${t}_test.py" >/tmp/${t}_out.txt 2>&1; then echo "통과"
  else echo "실패 (로그: /tmp/${t}_out.txt)"; fail=1; fi
done
# 회귀 SHA는 **library 모드 기준**이다 — 기본값(vector_homology)은 흉터가 빠져 서열이 다르다.
printf "%-6s " "회귀"
out=$(timeout 900 "$PY" "$HERE/regression_sha.py" fixed 2>&1 | grep '^RESULT')
if echo "$out" | grep -q '601809584fd54212' && echo "$out" | grep -q 'd69f339665bad36a'; then
  echo "통과 (baseline SHA 일치)"
else
  echo "실패 — 서열이 바뀌었습니다"; echo "$out"; fail=1
fi
exit $fail
