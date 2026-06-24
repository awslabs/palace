#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git diff --check
LAST_MEASURE="${PALACE_AR_LAST_MEASURE:-/home/ubuntu/workspace/palace_ar_paraview_proxy_last_measure.env}"
if [[ ! -f "$LAST_MEASURE" ]]; then
  echo "Missing last measure status: $LAST_MEASURE" >&2
  exit 1
fi
# shellcheck source=/dev/null
source "$LAST_MEASURE"
if [[ "${PROXY_OK:-0}" != "1" ]]; then
  echo "Proxy writer benchmark failed. Log: ${LOG:-unknown}" >&2
  exit 1
fi
