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
for t in cp0 cp1 cp2 cp3 cp5 cp6 pad cp7 cp7_excise cp7_e2e vec_verify real_vector one_vector mm_linear mm_circular vector_homology g1 g3 h1 y1 y2 y3 y4 y5 y6 y7; do
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
