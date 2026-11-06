#!/usr/bin/env bash
# 단계별 풀회귀 실행 → develop_log/ljw/regression/<stage>/<YYYYMMDD_HHMMSS>/ 에 저장.
#   results.tsv  : per-id (id<TAB>status<TAB>seqSHA<TAB>detail)  status=pass|alt_method|error|no_output|fail
#   summary.txt  : preload/progress/TOTAL 집계
#
# 사용:
#   awsemble/tools_dev/run_regression.sh <stage> [--old]
#     <stage> : 단계 라벨 (예: old, b1, h4, c1, p3 ...)
#     --old   : src 마운트 없이 이미지 baked(pre-refactor) 코드로 실행 (옛 baseline용)
# 예:
#   awsemble/tools_dev/run_regression.sh p3            # 현재 host src로
#   awsemble/tools_dev/run_regression.sh old --old     # 옛 코드 baseline
set -euo pipefail
STAGE="${1:?stage 라벨 필요 (예: h4, c1, p3)}"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"   # BioCAD_private
TS="$(date +%Y%m%d_%H%M%S)"
OUT="$ROOT/develop_log/ljw/regression/$STAGE/$TS"
mkdir -p "$OUT"

SRC_MOUNT=(-v "$ROOT/awsemble/src:/app/src:ro")
[ "${2:-}" = "--old" ] && SRC_MOUNT=()      # baked 옛 코드 사용

docker run --rm --network host --user "$(id -u):$(id -g)" -e HOME=/tmp -e PYTHONDONTWRITEBYTECODE=1 \
  -v "$ROOT/awsemble/.env:/app/.env:ro" \
  -v "$ROOT/awsemble/tools_dev:/app/tools_dev:ro" \
  "${SRC_MOUNT[@]}" \
  awsemble:dev python tools_dev/verify_full.py > "$OUT/results.tsv" 2> "$OUT/summary.txt"

echo "saved → $OUT"
grep TOTAL "$OUT/summary.txt" || true
