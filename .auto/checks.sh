#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git diff --check
# Keep this fast. Heavy correctness/performance checks live in .auto/measure.sh because
# they need the remote GPU logs and generated scalar outputs.
