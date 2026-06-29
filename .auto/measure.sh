#!/usr/bin/env bash
set -euo pipefail

ROOT=/home/ubuntu/palace
SPACK_SETUP=/home/ubuntu/spack/share/spack/setup-env.sh
API=/tmp/libcuda_api_counter.so
TAG=${PALACE_MEASURE_TAG:-$(git -C "$ROOT" rev-parse --short HEAD)-$(date +%Y%m%d%H%M%S)}
CFG=/tmp/palace_auto_transmon_ncfalse_${TAG}.json
LOG=${PALACE_MEASURE_LOG:-/tmp/palace_auto_transmon_ncfalse_${TAG}.log}
EXIT=${PALACE_MEASURE_EXIT:-/tmp/palace_auto_transmon_ncfalse_${TAG}.exit}
SRC=${PALACE_MEASURE_SRC:-/tmp/transmon_amr_branch_fixed.json}
OUT=/tmp/palace_auto_transmon_ncfalse_${TAG}_postpro

parse_metrics() {
  local log=$1
  local exit_code=${2:-999}
  python3 - "$log" "$exit_code" <<'PY'
import re, sys
log_path=sys.argv[1]
exit_code=int(sys.argv[2])
try:
    text=open(log_path, errors='replace').read()
except FileNotFoundError:
    text=''
oom=1 if ('CUDA_ERROR_OUT_OF_MEMORY' in text or 'OUT_OF_MEMORY' in text or 'Out of memory' in text) else 0
completed=1 if re.search(r'Completed\s+\d+\s+iterations of adaptive mesh refinement', text) else 0
# Palace prints multiple timing tables. Use the last Total for successful runs.
totals=[float(x) for x in re.findall(r'^Total\s+([0-9]+(?:\.[0-9]+)?)\s+', text, re.M)]
metric=totals[-1] if completed and exit_code==0 and totals else 999999.0
initial_elements=0
refined_elements=0
m=re.search(r'^\s*elements\s+(\d+)\s+', text, re.M)
if m:
    initial_elements=int(m.group(1))
m=re.search(r'Conforming mesh refinement added\s+\d+\s+elements \(initial = \d+, final = (\d+)\)', text)
if m:
    refined_elements=int(m.group(1))
bdr_elems=[int(x) for x in re.findall(r'SurfaceFunctional profile kind=BDR_[^\n]* elems=(\d+)', text)]
bdr_max=max(bdr_elems) if bdr_elems else 0
# Initial timing table before AMR iteration 1, if present.
pre=text.split('Adaptive mesh refinement (AMR) iteration 1:')[0]
def last_timer(label, block):
    vals=[float(x) for x in re.findall(r'^\s*' + re.escape(label) + r'\s+([0-9]+(?:\.[0-9]+)?)\s+', block, re.M)]
    return vals[-1] if vals else 0.0
nvrtc=0
modload=0
m=re.search(r'\[cuda_api_counter\] nvrtcCompileProgram=(\d+)', text)
if m: nvrtc=int(m.group(1))
m=re.search(r'\[cuda_api_counter\] cuModuleLoadData=(\d+)', text)
if m: modload=int(m.group(1))
print(f'METRIC transmon_ncfalse_penalty_s={metric}')
print(f'METRIC transmon_ncfalse_oom={oom}')
print(f'METRIC transmon_ncfalse_exit={exit_code}')
print(f'METRIC transmon_ncfalse_completed={completed}')
print(f'METRIC transmon_ncfalse_initial_elements={initial_elements}')
print(f'METRIC transmon_ncfalse_refined_elements={refined_elements}')
print(f'METRIC transmon_ncfalse_bdr_elems_max={bdr_max}')
print(f'METRIC transmon_ncfalse_initial_paraview_s={last_timer("Paraview", pre)}')
print(f'METRIC transmon_ncfalse_initial_gridfunction_s={last_timer("Grid function", pre)}')
print(f'METRIC transmon_ncfalse_nvrtc={nvrtc}')
print(f'METRIC transmon_ncfalse_cuModuleLoadData={modload}')
print(f'LOG {log_path}')
PY
}

if [[ -n "${PALACE_MEASURE_EXISTING_LOG:-}" ]]; then
  parse_metrics "$PALACE_MEASURE_EXISTING_LOG" "${PALACE_MEASURE_EXISTING_EXIT:-1}"
  exit 0
fi

cd "$ROOT"
source "$SPACK_SETUP"
# Spack only; keep parallelism at -j4.
yes | spack -e "$ROOT" install --overwrite -j4 >/tmp/palace_auto_spack_${TAG}.log 2>&1
PREFIX=$(spack -e "$ROOT" location -i local.palace)
BIN=$PREFIX/bin/palace-x86_64.bin

python3 - "$SRC" "$CFG" "$OUT" <<'PY'
import json, sys
src, dst, out = sys.argv[1:4]
with open(src) as f:
    d=json.load(f)
d.setdefault('Problem', {})['Output']=out
d.setdefault('Problem', {}).setdefault('OutputFormats', {})['Paraview']=True
d['Problem']['OutputFormats']['GridFunction']=True
d.setdefault('Solver', {})['Device']='GPU'
ref=d.setdefault('Model', {}).setdefault('Refinement', {})
ref['Nonconformal']=False
ref['MaxIts']=1
# Make mesh absolute if the source was relative.
mesh=d.setdefault('Model', {}).get('Mesh', '')
if mesh and not mesh.startswith('/'):
    d['Model']['Mesh']='/home/ubuntu/palace/examples/transmon/' + mesh
with open(dst,'w') as f:
    json.dump(d,f,indent=2)
PY

rm -f "$LOG" "$EXIT"
set +e
if [[ -f "$API" ]]; then
  env PALACE_SURFACE_PROFILE=1 LD_PRELOAD="$API" "$BIN" "$CFG" >"$LOG" 2>&1
else
  env PALACE_SURFACE_PROFILE=1 "$BIN" "$CFG" >"$LOG" 2>&1
fi
rc=$?
set -e
echo "$rc" >"$EXIT"
parse_metrics "$LOG" "$rc"
# Keep the script exit code zero so crash/OOM experiments are logged with penalty metrics.
exit 0
