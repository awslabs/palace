#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

SPACK_ROOT="${SPACK_ROOT:-/home/ubuntu/spack-development/spack}"
# shellcheck source=/dev/null
source "$SPACK_ROOT/share/spack/setup-env.sh"

export PALACE_AR_BUILD_JOBS="${PALACE_AR_BUILD_JOBS:-32}"
export PALACE_AR_GPU_NP="${PALACE_AR_GPU_NP:-8}"
RUN_ROOT="${PALACE_AR_RUN_ROOT:-/home/ubuntu/workspace/palace_ar_paraview_stream}"
BASELINE_DIR="${PALACE_AR_BASELINE_DIR:-$RUN_ROOT/baseline_csv}"
POST_ROOT="${PALACE_AR_POST_ROOT:-/home/ubuntu/workspace/palace_ar_paraview_stream_latest_postpro}"
REFERENCE_POST_ROOT="${PALACE_AR_REFERENCE_POST_ROOT:-/home/ubuntu/workspace/palace_ar_paraview_stream_origin_main_reference_postpro}"
REFERENCE_POST_DIR="$REFERENCE_POST_ROOT/transmon_amr"
LAST_MEASURE="${PALACE_AR_LAST_MEASURE:-/home/ubuntu/workspace/palace_ar_paraview_stream_last_measure.env}"
SHA="$(git rev-parse --short HEAD)"
RUN_ID="${SHA}-$(date +%Y%m%dT%H%M%SZ)"
OUT_DIR="$RUN_ROOT/runs/$RUN_ID"
CONFIG="$OUT_DIR/transmon_amr_gpu_paraview.json"
LOG="$OUT_DIR/palace.log"
BUILD_LOG="$OUT_DIR/build.log"
POST_DIR="$POST_ROOT/transmon_amr"
mkdir -p "$RUN_ROOT/runs"
# Keep scalar CSV references plus, when seeded separately, an origin/main
# ParaView/GridFunction reference tree at REFERENCE_POST_ROOT. Large output from
# each candidate run is written to POST_ROOT on the instance NVMe and overwritten
# each run; do not accumulate one multi-GB postpro tree per trial.
find "$RUN_ROOT/runs" -mindepth 1 -maxdepth 1 -type d -exec rm -rf {} +
rm -rf "$POST_ROOT"
mkdir -p "$OUT_DIR" "$POST_ROOT"

if [[ ! -f spack.lock ]]; then
  spack -e "$ROOT" concretize -f >"$OUT_DIR/concretize.log" 2>&1 || {
    tail -120 "$OUT_DIR/concretize.log" >&2
    exit 1
  }
fi

spack -e "$ROOT" install --only package --overwrite -y -j"$PALACE_AR_BUILD_JOBS" \
  >"$BUILD_LOG" 2>&1 || {
    tail -160 "$BUILD_LOG" >&2
    exit 1
  }

python3 - "$CONFIG" "$POST_DIR" <<'PY'
import json
import pathlib
import sys
cfg_path = pathlib.Path(sys.argv[1])
out_dir = pathlib.Path(sys.argv[2])
root = pathlib.Path.cwd()
with open(root / "examples/transmon/transmon_amr.json") as f:
    cfg = json.load(f)
cfg["Problem"]["Output"] = str(out_dir)
cfg.setdefault("Problem", {}).setdefault("OutputFormats", {})["Paraview"] = True
cfg["Problem"]["OutputFormats"]["GridFunction"] = True
cfg.setdefault("Solver", {})["Device"] = "GPU"
cfg_path.write_text(json.dumps(cfg, indent=2) + "\n")
PY

PALACE_PREFIX="$(spack -e "$ROOT" location -i palace)"
PALACE_BIN="$PALACE_PREFIX/bin/palace"
rm -rf "$POST_DIR"

set +e
spack -e "$ROOT" build-env palace -- bash -lc 'cd "$1/examples/transmon" && "$2" -np "$3" "$4"' _ "$ROOT" "$PALACE_BIN" "$PALACE_AR_GPU_NP" "$CONFIG" \
  >"$LOG" 2>&1
run_status=$?
set -e
if [[ $run_status -ne 0 ]]; then
  tail -180 "$LOG" >&2
  exit "$run_status"
fi

parse_time() {
  local kind="$1"
  awk -v kind="$kind" '
    /^Elapsed Time Report/ {in_time=1; next}
    /^Peak Memory/ {in_time=0}
    !in_time {next}
    kind == "Paraview" && $1 == "Paraview" {v=$4}
    kind == "Postprocessing" && $1 == "Postprocessing" {v=$4}
    kind == "Grid function" && $1 == "Grid" && $2 == "function" {v=$5}
    kind == "Operator Construction" && $1 == "Operator" && $2 == "Construction" {v=$5}
    kind == "Total" && $1 == "Total" {v=$4}
    END { if (v == "") exit 1; print v }
  ' "$LOG"
}

paraview="$(parse_time "Paraview")"
postprocessing="$(parse_time "Postprocessing")"
gridfunction="$(parse_time "Grid function")"
operator_construction="$(parse_time "Operator Construction")"
total="$(parse_time "Total")"

compare_output="$(python3 - "$POST_DIR" "$BASELINE_DIR" <<'PY'
import csv
import pathlib
import re
import shutil
import sys
post = pathlib.Path(sys.argv[1])
base = pathlib.Path(sys.argv[2])
float_re = re.compile(r'^[+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?$')
files = sorted(p for p in post.rglob('*.csv') if p.is_file())
if not files:
    print('scalar_ok=0 scalar_max_rel=inf scalar_max_abs=inf scalar_compared_files=0')
    sys.exit(0)
if not base.exists():
    base.mkdir(parents=True, exist_ok=True)
    for p in files:
        dest = base / p.relative_to(post)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dest)
    print(f'scalar_ok=1 scalar_max_rel=0 scalar_max_abs=0 scalar_compared_files={len(files)}')
    sys.exit(0)

def nums(path):
    vals = []
    with open(path, newline='') as f:
        for row in csv.reader(f):
            for cell in row:
                s = cell.strip()
                if float_re.match(s):
                    vals.append(float(s))
    return vals
max_rel = 0.0
max_abs = 0.0
compared = 0
ok = 1
for p in files:
    rel = p.relative_to(post)
    q = base / rel
    if not q.exists():
        ok = 0
        continue
    a = nums(q)
    b = nums(p)
    if len(a) != len(b):
        ok = 0
        continue
    compared += 1
    for x, y in zip(a, b):
        da = abs(y - x)
        dr = da / max(abs(x), 1.0e-30)
        max_abs = max(max_abs, da)
        max_rel = max(max_rel, dr)
# Output-only changes should be bitwise or roundoff-close in scalar CSVs. Use a loose
# gate so harmless eigensolver residual jitter does not block output-writer experiments.
if max_abs > 1.0e-7 and max_rel > 1.0e-4:
    ok = 0
print(f'scalar_ok={ok} scalar_max_rel={max_rel:.17g} scalar_max_abs={max_abs:.17g} scalar_compared_files={compared}')
PY
)"

eval "$compare_output"
cat >"$LAST_MEASURE" <<EOF2
SCALAR_OK=${scalar_ok:-0}
SCALAR_MAX_REL=${scalar_max_rel:-inf}
SCALAR_MAX_ABS=${scalar_max_abs:-inf}
SCALAR_COMPARED_FILES=${scalar_compared_files:-0}
RUN_ID=$RUN_ID
LOG=$LOG
POST_DIR=$POST_DIR
REFERENCE_POST_DIR=$REFERENCE_POST_DIR
EOF2

printf 'METRIC transmon_amr_paraview_seconds=%s\n' "$paraview"
printf 'METRIC transmon_amr_total_seconds=%s\n' "$total"
printf 'METRIC transmon_amr_postprocessing_seconds=%s\n' "$postprocessing"
printf 'METRIC transmon_amr_gridfunction_seconds=%s\n' "$gridfunction"
printf 'METRIC transmon_amr_operator_construction_seconds=%s\n' "$operator_construction"
printf 'METRIC scalar_max_rel=%s\n' "${scalar_max_rel:-inf}"
printf 'METRIC scalar_max_abs=%s\n' "${scalar_max_abs:-inf}"
printf 'METRIC scalar_compared_files=%s\n' "${scalar_compared_files:-0}"
printf 'METRIC scalar_ok=%s\n' "${scalar_ok:-0}"
printf 'ARTIFACT log=%s\n' "$LOG"
printf 'ARTIFACT post_dir=%s\n' "$POST_DIR"
if [[ -d "$REFERENCE_POST_DIR" ]]; then
  printf 'ARTIFACT reference_post_dir=%s\n' "$REFERENCE_POST_DIR"
fi
