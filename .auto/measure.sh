#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

SPACK_ROOT="${SPACK_ROOT:-/home/ubuntu/spack-development/spack}"
# shellcheck source=/dev/null
source "$SPACK_ROOT/share/spack/setup-env.sh"

export PALACE_AR_BUILD_JOBS="${PALACE_AR_BUILD_JOBS:-32}"
export PALACE_PROXY_MPI_NP="${PALACE_PROXY_MPI_NP:-2}"
RUN_ROOT="${PALACE_AR_RUN_ROOT:-/home/ubuntu/workspace/palace_ar_paraview_proxy}"
POST_ROOT="${PALACE_AR_POST_ROOT:-/home/ubuntu/workspace/palace_ar_paraview_proxy_latest}"
LAST_MEASURE="${PALACE_AR_LAST_MEASURE:-/home/ubuntu/workspace/palace_ar_paraview_proxy_last_measure.env}"
SHA="$(git rev-parse --short HEAD)"
RUN_ID="${SHA}-$(date +%Y%m%dT%H%M%SZ)"
OUT_DIR="$RUN_ROOT/runs/$RUN_ID"
LOG="$OUT_DIR/proxy.log"
BUILD_LOG="$OUT_DIR/build.log"
CONCRETIZE_LOG="$OUT_DIR/concretize.log"
PROXY_OUTPUT_DIR="$POST_ROOT/proxy-output"

mkdir -p "$RUN_ROOT/runs"
# Keep only the current proxy artifact tree. This benchmark is much smaller than
# the full solver workload, but still large enough that accumulation obscures
# disk-use failures during long autonomous runs.
find "$RUN_ROOT/runs" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
rm -rf "$POST_ROOT"
mkdir -p "$OUT_DIR" "$POST_ROOT"

if [[ ! -f spack.lock ]] || ! grep -q '"name":"catch2"' spack.lock; then
  spack -e "$ROOT" concretize -f --test root >"$CONCRETIZE_LOG" 2>&1 || {
    tail -160 "$CONCRETIZE_LOG" >&2 || true
    exit 1
  }
fi

: >"$BUILD_LOG"
if spack -e "$ROOT" find -lv catch2 >/dev/null 2>&1; then
  spack -e "$ROOT" install --test=root --only package --overwrite --keep-stage -y -j"$PALACE_AR_BUILD_JOBS" \
    >>"$BUILD_LOG" 2>&1 || {
      tail -180 "$BUILD_LOG" >&2 || true
      exit 1
    }
else
  spack -e "$ROOT" install --test=root --overwrite --keep-stage -y -j"$PALACE_AR_BUILD_JOBS" \
    >>"$BUILD_LOG" 2>&1 || {
      tail -180 "$BUILD_LOG" >&2 || true
      exit 1
    }
fi

PALACE_PREFIX="$(spack -e "$ROOT" location -i palace)"
UNIT_TESTS="$PALACE_PREFIX/bin/palace-unit-tests"
if [[ ! -x "$UNIT_TESTS" ]]; then
  echo "Missing palace-unit-tests at $UNIT_TESTS" >&2
  tail -120 "$BUILD_LOG" >&2 || true
  exit 1
fi

rm -rf "$PROXY_OUTPUT_DIR"
mkdir -p "$PROXY_OUTPUT_DIR"

set +e
PALACE_PROXY_OUTPUT_DIR="$PROXY_OUTPUT_DIR" \
  PALACE_PROXY_ORDER="${PALACE_PROXY_ORDER:-3}" \
  PALACE_PROXY_REFINE_STEPS="${PALACE_PROXY_REFINE_STEPS:-2}" \
  PALACE_PROXY_REFINE_PROB="${PALACE_PROXY_REFINE_PROB:-0.3}" \
  PALACE_PROXY_SEED="${PALACE_PROXY_SEED:-20260624}" \
  PALACE_PROXY_COMPRESSION="${PALACE_PROXY_COMPRESSION:-1}" \
  spack -e "$ROOT" build-env palace -- mpirun -n "$PALACE_PROXY_MPI_NP" \
    "$UNIT_TESTS" --device "${PALACE_PROXY_DEVICE:-cuda}" \
    --backend "${PALACE_PROXY_CEED_BACKEND:-/gpu/cuda/magma}" \
    "[paraview-writer-proxy]" --skip-benchmarks \
    >"$LOG" 2>&1
run_status=$?
set -e
if [[ $run_status -ne 0 ]]; then
  tail -180 "$LOG" >&2 || true
  exit "$run_status"
fi

metric() {
  local key="$1"
  awk -F= -v key="$key" '$1 == "METRIC " key {v=$2} END {if (v == "") exit 1; print v}' "$LOG"
}

proxy_paraview_seconds="$(metric proxy_paraview_seconds)"
proxy_auto2_like_seconds="$(metric proxy_auto2_like_seconds)"
proxy_legacy_coefficient_seconds="$(metric proxy_legacy_coefficient_seconds)"
proxy_materialize_seconds="$(metric proxy_materialize_seconds)"
proxy_materialized_save_seconds="$(metric proxy_materialized_save_seconds)"
proxy_direct_buffer_seconds="$(metric proxy_direct_buffer_seconds)"
proxy_direct_buffer_eval_seconds="$(metric proxy_direct_buffer_eval_seconds)"
proxy_direct_buffer_save_seconds="$(metric proxy_direct_buffer_save_seconds)"
proxy_materialized_vs_legacy="$(metric proxy_materialized_vs_legacy)"
proxy_direct_vs_materialized="$(metric proxy_direct_vs_materialized)"
proxy_legacy_bytes="$(metric proxy_legacy_bytes)"
proxy_materialized_bytes="$(metric proxy_materialized_bytes)"
proxy_direct_buffer_bytes="$(metric proxy_direct_buffer_bytes)"
proxy_mpi_ranks="$(metric proxy_mpi_ranks)"
proxy_elements="$(metric proxy_elements)"
proxy_e_true_dofs="$(metric proxy_e_true_dofs)"
proxy_b_true_dofs="$(metric proxy_b_true_dofs)"
proxy_scalar_true_dofs="$(metric proxy_scalar_true_dofs)"
proxy_vector_true_dofs="$(metric proxy_vector_true_dofs)"

cat >"$LAST_MEASURE" <<EOF2
PROXY_OK=1
RUN_ID=$RUN_ID
LOG=$LOG
BUILD_LOG=$BUILD_LOG
PROXY_OUTPUT_DIR=$PROXY_OUTPUT_DIR
PROXY_PARAVIEW_SECONDS=$proxy_paraview_seconds
PROXY_AUTO2_LIKE_SECONDS=$proxy_auto2_like_seconds
PROXY_LEGACY_COEFFICIENT_SECONDS=$proxy_legacy_coefficient_seconds
PROXY_MATERIALIZE_SECONDS=$proxy_materialize_seconds
PROXY_MATERIALIZED_SAVE_SECONDS=$proxy_materialized_save_seconds
PROXY_DIRECT_BUFFER_SECONDS=$proxy_direct_buffer_seconds
PROXY_DIRECT_BUFFER_EVAL_SECONDS=$proxy_direct_buffer_eval_seconds
PROXY_DIRECT_BUFFER_SAVE_SECONDS=$proxy_direct_buffer_save_seconds
PROXY_MATERIALIZED_VS_LEGACY=$proxy_materialized_vs_legacy
PROXY_DIRECT_VS_MATERIALIZED=$proxy_direct_vs_materialized
PROXY_LEGACY_BYTES=$proxy_legacy_bytes
PROXY_MATERIALIZED_BYTES=$proxy_materialized_bytes
PROXY_DIRECT_BUFFER_BYTES=$proxy_direct_buffer_bytes
PROXY_MPI_RANKS=$proxy_mpi_ranks
PROXY_ELEMENTS=$proxy_elements
PROXY_E_TRUE_DOFS=$proxy_e_true_dofs
PROXY_B_TRUE_DOFS=$proxy_b_true_dofs
PROXY_SCALAR_TRUE_DOFS=$proxy_scalar_true_dofs
PROXY_VECTOR_TRUE_DOFS=$proxy_vector_true_dofs
EOF2

printf 'METRIC proxy_paraview_seconds=%s\n' "$proxy_paraview_seconds"
printf 'METRIC proxy_auto2_like_seconds=%s\n' "$proxy_auto2_like_seconds"
printf 'METRIC proxy_legacy_coefficient_seconds=%s\n' "$proxy_legacy_coefficient_seconds"
printf 'METRIC proxy_materialize_seconds=%s\n' "$proxy_materialize_seconds"
printf 'METRIC proxy_materialized_save_seconds=%s\n' "$proxy_materialized_save_seconds"
printf 'METRIC proxy_direct_buffer_seconds=%s\n' "$proxy_direct_buffer_seconds"
printf 'METRIC proxy_direct_buffer_eval_seconds=%s\n' "$proxy_direct_buffer_eval_seconds"
printf 'METRIC proxy_direct_buffer_save_seconds=%s\n' "$proxy_direct_buffer_save_seconds"
printf 'METRIC proxy_materialized_vs_legacy=%s\n' "$proxy_materialized_vs_legacy"
printf 'METRIC proxy_direct_vs_materialized=%s\n' "$proxy_direct_vs_materialized"
printf 'METRIC proxy_legacy_bytes=%s\n' "$proxy_legacy_bytes"
printf 'METRIC proxy_materialized_bytes=%s\n' "$proxy_materialized_bytes"
printf 'METRIC proxy_direct_buffer_bytes=%s\n' "$proxy_direct_buffer_bytes"
printf 'METRIC proxy_mpi_ranks=%s\n' "$proxy_mpi_ranks"
printf 'METRIC proxy_elements=%s\n' "$proxy_elements"
printf 'METRIC proxy_e_true_dofs=%s\n' "$proxy_e_true_dofs"
printf 'METRIC proxy_b_true_dofs=%s\n' "$proxy_b_true_dofs"
printf 'METRIC proxy_scalar_true_dofs=%s\n' "$proxy_scalar_true_dofs"
printf 'METRIC proxy_vector_true_dofs=%s\n' "$proxy_vector_true_dofs"
printf 'METRIC proxy_ok=1\n'
printf 'ARTIFACT log=%s\n' "$LOG"
printf 'ARTIFACT build_log=%s\n' "$BUILD_LOG"
printf 'ARTIFACT proxy_output_dir=%s\n' "$PROXY_OUTPUT_DIR"
