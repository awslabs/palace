#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
RUN_ID="$(git rev-parse --short HEAD)-$(date +%Y%m%d%H%M%S)"
CASE_SET="${PALACE_AR_CASES:-cpw}"
exec .auto/remote_measure_payload.sh "$RUN_ID" "$CASE_SET"
