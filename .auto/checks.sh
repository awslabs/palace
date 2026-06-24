#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git diff --check
LAST_MEASURE="${PALACE_AR_LAST_MEASURE:-/home/ubuntu/workspace/palace_ar_paraview_stream_last_measure.env}"
if [[ ! -f "$LAST_MEASURE" ]]; then
  echo "Missing last measure status: $LAST_MEASURE" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$LAST_MEASURE"
if [[ "${SCALAR_OK:-0}" != "1" ]]; then
  echo "Scalar output check failed: max_rel=${SCALAR_MAX_REL:-unknown}, max_abs=${SCALAR_MAX_ABS:-unknown}, files=${SCALAR_COMPARED_FILES:-unknown}" >&2
  echo "Log: ${LOG:-unknown}" >&2
  exit 1
fi
