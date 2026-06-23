#!/usr/bin/env bash
set -euo pipefail
cd "$(git rev-parse --show-toplevel)"
git diff --check

# Processor-boundary gate: every kept experiment must also complete the CPW
# postprocessing workload on 16 CPU MPI ranks. The primary metric remains the
# one-rank GPU CPW postprocessing time, but multi-rank correctness is a hard
# pass/fail condition because ghost/processor-boundary surfaces must stay on
# the safe path.
if [[ -x .auto/verify_mpi_cpu.sh ]]; then
  PALACE_AR_MPI_CPU_RANKS=16 .auto/verify_mpi_cpu.sh
fi
