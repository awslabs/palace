#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
RUN_ROOT=${PALACE_AR_RUN_ROOT:-/home/ubuntu/scratch/palace_nsight_autoresearch}
LATEST=${PALACE_AR_LATEST:-$RUN_ROOT/latest}
NP=${PALACE_AR_GPU_NP:-8}
JOBS=${PALACE_AR_BUILD_JOBS:-96}
mkdir -p "$LATEST" "$LATEST/profile"
rm -rf "$LATEST/output" "$LATEST/profile"/*
mkdir -p "$LATEST/output" "$LATEST/profile"

source /home/ubuntu/spack-development/spack/share/spack/setup-env.sh
cd "$ROOT"

if [[ "${PALACE_AR_SKIP_BUILD:-0}" != "1" ]]; then
  spack -e "$ROOT" install --only package --overwrite --keep-stage -y -j "$JOBS" \
    > "$LATEST/build.log" 2>&1
fi

PREFIX=$(spack -e "$ROOT" location -i palace)
CUDA_PREFIX=$(spack -e "$ROOT" location -i cuda)
NSYS_BIN="$CUDA_PREFIX/bin/nsys"
PALACE_BIN="$PREFIX/bin/palace-x86_64.bin"
CONFIG="$LATEST/transmon_amr1_pv_only.json"
PALACE_LOG="$LATEST/palace.log"
NSYS_STATS="$LATEST/nsys_stats.txt"

python3 - "$ROOT/examples/transmon/transmon_amr.json" "$CONFIG" "$LATEST/output" <<'PY'
import json, sys
src, dst, out = sys.argv[1:]
cfg = json.load(open(src))
cfg.setdefault('Problem', {})['Output'] = out
cfg['Problem'].setdefault('OutputFormats', {})['Paraview'] = True
cfg['Problem']['OutputFormats']['GridFunction'] = False
cfg.setdefault('Solver', {})['Device'] = 'GPU'
cfg.setdefault('Model', {}).setdefault('Refinement', {})['MaxIts'] = 1
json.dump(cfg, open(dst, 'w'), indent=2)
PY

cat > "$LATEST/rank_wrapper.sh" <<'EOS'
#!/usr/bin/env bash
set -euo pipefail
rank=${PMI_RANK:-0}
if [[ "$rank" == "0" ]]; then
  exec "$NSYS_BIN" profile \
    --trace=cuda,nvtx,osrt,mpi \
    --mpi-impl=mpich \
    --sample=none \
    --osrt-file-access=true \
    --force-overwrite=true \
    --stats=false \
    --show-output=true \
    --output "$PROFROOT/rank0" \
    "$PALACE_BIN" "$PALACE_CFG"
else
  exec "$PALACE_BIN" "$PALACE_CFG"
fi
EOS
chmod +x "$LATEST/rank_wrapper.sh"

export NSYS_BIN PROFROOT="$LATEST/profile" PALACE_BIN PALACE_CFG="$CONFIG"
set +e
(
  cd "$ROOT/examples/transmon"
  spack -e "$ROOT" build-env palace -- mpirun -n "$NP" "$LATEST/rank_wrapper.sh"
) > "$PALACE_LOG" 2>&1
run_rc=$?
set -e

report=$(ls "$LATEST"/profile/rank0*.nsys-rep 2>/dev/null | head -1 || true)
if [[ -n "$report" ]]; then
  "$NSYS_BIN" stats --force-export=true \
    --report cuda_api_sum,cuda_gpu_sum,cuda_gpu_mem_time_sum,cuda_gpu_mem_size_sum,osrt_sum,syscall_sum,mpi_event_sum \
    --format column "$report" > "$NSYS_STATS" 2>&1 || true
else
  : > "$NSYS_STATS"
fi

python3 - "$PALACE_LOG" "$LATEST/output" "$report" "$NSYS_STATS" "$run_rc" <<'PY'
import os, pathlib, re, subprocess, sys
log_path, out_dir, report, stats_path, rc_s = sys.argv[1:]
rc = int(rc_s)
text = pathlib.Path(log_path).read_text(errors='replace') if pathlib.Path(log_path).exists() else ''
stats = pathlib.Path(stats_path).read_text(errors='replace') if pathlib.Path(stats_path).exists() else ''

def last_timer_row(name):
    val = None
    pat = re.compile(r"^\s*" + re.escape(name) + r"\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$")
    for line in text.splitlines():
        m = pat.match(line)
        if m:
            val = tuple(float(x) for x in m.groups())
    return val

def last_memory_row(name):
    val = None
    in_mem = False
    pat = re.compile(r"^\s*" + re.escape(name) + r"\s+([0-9.]+)([KMG]?)\s+([0-9.]+)([KMG]?)\s+([0-9.]+)([KMG]?)\s*$")
    for line in text.splitlines():
        if line.startswith('Peak Memory'):
            in_mem = True
            continue
        if in_mem:
            m = pat.match(line)
            if m:
                val = m.groups()
    return val

def to_gb(groups, idx):
    if not groups:
        return 0.0
    n = float(groups[idx]); unit = groups[idx+1]
    if unit == 'K': return n / (1024*1024)
    if unit == 'M': return n / 1024
    if unit == 'G': return n
    return n / (1024**3)

def du_bytes(path):
    try:
        return int(subprocess.check_output(['du','-sb',path], text=True).split()[0])
    except Exception:
        return 0

pv = last_timer_row('Paraview')
total = last_timer_row('Total')
post = last_timer_row('Postprocessing')
op = last_timer_row('Operator Construction')
mem_pv = last_memory_row('Paraview')
report_bytes = os.path.getsize(report) if report and os.path.exists(report) else 0
print(f"METRIC paraview_seconds={pv[2] if pv else 999999.0}")
print(f"METRIC total_seconds={total[2] if total else 999999.0}")
print(f"METRIC postprocessing_seconds={post[2] if post else 0.0}")
print(f"METRIC operator_construction_seconds={op[2] if op else 0.0}")
print(f"METRIC paraview_memory_hwm_gb={to_gb(mem_pv,4) if mem_pv else 0.0}")
print(f"METRIC output_bytes={du_bytes(out_dir)}")
print(f"METRIC nsys_report_bytes={report_bytes}")
print(f"METRIC nsys_stats_ready={1 if stats.strip() else 0}")
print(f"METRIC cuda_log={1 if 'Device configuration: cuda,cpu' in text else 0}")
print(f"METRIC ceed_gpu={1 if 'libCEED backend: /gpu/cuda' in text else 0}")
print(f"METRIC run_rc={rc}")
print(f"PROFILE palace_log={log_path}")
print(f"PROFILE nsys_stats={stats_path}")
print(f"PROFILE nsys_report={report or 'NONE'}")
if rc != 0:
    print('RUN_FAILED')
    sys.exit(1)
if 'Device configuration: cuda,cpu' not in text or 'libCEED backend: /gpu/cuda' not in text:
    print('GPU_CHECK_FAILED')
    sys.exit(2)
PY
