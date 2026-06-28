#!/usr/bin/env bash
set -euo pipefail
# Diagnostic/no-fallback objective for the NC boundary trace-plan experiment.
# This intentionally opts back into the libCEED nonconforming ParaView point-field path so
# autoresearch can improve that implementation instead of winning via the legacy fallback.
export PALACE_CEED_NONCONFORMING_PARAVIEW=1
export PALACE_SURFACE_PROFILE=1
exec "$(dirname "${BASH_SOURCE[0]}")/measure.sh" "$@"
