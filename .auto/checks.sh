#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
git diff --check

if [[ "${PALACE_AR_SKIP_RUNTIME_CHECKS:-0}" == "1" ]]; then
  echo "Skipping runtime validation because PALACE_AR_SKIP_RUNTIME_CHECKS=1"
  exit 0
fi

SPACK_SETUP="${SPACK_SETUP:-/home/ubuntu/spack-development/spack/share/spack/setup-env.sh}"
BUILD_JOBS="${PALACE_AR_BUILD_JOBS:-32}"
GPU_NP="${PALACE_AR_GPU_NP:-$(nvidia-smi -L | wc -l)}"
CPU_NP="${PALACE_AR_CPU_NP:-$(nproc)}"
VALIDATION_CASES="${PALACE_AR_VALIDATION_CASES:-spheres rings cpw/lumped_uniform cpw/wave_uniform antenna/antenna_short_dipole}"

source "$SPACK_SETUP"
spack -e "$ROOT" install -j"$BUILD_JOBS" >/tmp/palace_ar_checks_build.log 2>&1
PREFIX="$(spack -e "$ROOT" location -i local.palace)"
WRAPPER="$PREFIX/bin/palace"
if [[ ! -x "$WRAPPER" ]]; then
  echo "ERROR missing Palace wrapper: $WRAPPER" >&2
  tail -80 /tmp/palace_ar_checks_build.log >&2 || true
  exit 2
fi

export PATH="$HOME/.local/bin:$PATH"
export OMP_NUM_THREADS=1
# Avoid default OpenMPI CPU binding starvation when the concretized MPI is OpenMPI.
export OMPI_MCA_hwloc_base_binding_policy=none
export PRTE_MCA_hwloc_default_binding_policy=none

julia --project=test/examples -e 'using Pkg; Pkg.instantiate()' >/tmp/palace_ar_julia_instantiate.log 2>&1

echo "VALIDATION gpu_regression np=$GPU_NP cases=$VALIDATION_CASES"
julia --project=test/examples --color=yes test/examples/runtests.jl \
  --palace-test "$WRAPPER" \
  --num-proc-test "$GPU_NP" \
  --palace-device GPU \
  --test-cases $VALIDATION_CASES \
  2>&1 | tee /tmp/palace_ar_gpu_regression.log

echo "VALIDATION cpu_regression np=$CPU_NP cases=$VALIDATION_CASES"
julia --project=test/examples --color=yes test/examples/runtests.jl \
  --palace-test "$WRAPPER" \
  --num-proc-test "$CPU_NP" \
  --palace-device CPU \
  --test-cases $VALIDATION_CASES \
  2>&1 | tee /tmp/palace_ar_cpu_regression.log
