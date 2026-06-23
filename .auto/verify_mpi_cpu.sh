#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
RANKS="${PALACE_AR_MPI_CPU_RANKS:-16}"
RUN_ID="${PALACE_AR_CHECK_RUN_ID:-$(git rev-parse --short HEAD)-mpi${RANKS}cpu-$(date +%Y%m%d%H%M%S)}"
SPACK_SETUP="/home/ubuntu/spack/share/spack/setup-env.sh"

if [[ ! -f "$SPACK_SETUP" ]]; then
  echo "SKIP mpi_cpu_check: missing $SPACK_SETUP"
  exit 0
fi

source "$SPACK_SETUP"
cd "$ROOT"

PREFIX="$(spack -e "$ROOT" location -i local.palace 2>/dev/null || true)"
BIN="${PREFIX:+$PREFIX/bin/palace-x86_64.bin}"
if [[ -z "$PREFIX" || ! -x "$BIN" ]]; then
  echo "mpi_cpu_check: installed Palace binary missing; building first" >&2
  spack -e "$ROOT" install -j4 >/tmp/palace_ar_mpi${RANKS}_cpu_build_${RUN_ID}.log 2>&1
  PREFIX="$(spack -e "$ROOT" location -i local.palace)"
  BIN="$PREFIX/bin/palace-x86_64.bin"
fi
if [[ ! -x "$BIN" ]]; then
  echo "ERROR mpi_cpu_check: missing Palace binary: $BIN" >&2
  exit 2
fi
MPIEXEC="$(command -v mpirun || command -v mpiexec || true)"
if [[ -z "$MPIEXEC" ]]; then
  echo "ERROR mpi_cpu_check: mpirun/mpiexec not found" >&2
  exit 2
fi

python3 - "$RUN_ID" "$RANKS" <<'PY'
import json, re, sys
from pathlib import Path
run_id=sys.argv[1]
ranks=sys.argv[2]
root=Path('/home/ubuntu/palace')

def load_jsonc(path):
    s=Path(path).read_text()
    out=[]; ins=False; esc=False; i=0
    while i < len(s):
        c=s[i]
        if ins:
            out.append(c)
            if esc: esc=False
            elif c=='\\': esc=True
            elif c=='"': ins=False
            i+=1; continue
        if c=='"': ins=True; out.append(c); i+=1; continue
        if c=='/' and i+1 < len(s) and s[i+1]=='/':
            while i < len(s) and s[i] not in '\r\n': i+=1
            continue
        out.append(c); i+=1
    return json.loads(re.sub(r',\s*([}\]])', r'\1', ''.join(out)))

j=load_jsonc(root/'examples/cpw/cpw_lumped_uniform.json')
j.setdefault('Problem', {})['Output']=f'postpro/autoresearch_cpw_mpi{ranks}_cpu_{run_id}'
j['Problem'].setdefault('OutputFormats', {})['Paraview']=True
j['Problem']['OutputFormats']['GridFunction']=True
j.setdefault('Solver', {})['Device']='CPU'
j['Solver'].setdefault('Driven', {})['Samples']=[{'Type':'Point','Freq':[17.0],'SaveStep':1}]
j['Solver']['Driven']['Save']=[17.0]
out=Path(f'/tmp/palace_ar_cpw_mpi{ranks}_cpu_{run_id}.json')
out.write_text(json.dumps(j, indent=2))
print(out)
PY

CFG="/tmp/palace_ar_cpw_mpi${RANKS}_cpu_${RUN_ID}.json"
LOG="/tmp/palace_ar_cpw_mpi${RANKS}_cpu_${RUN_ID}.log"
rm -rf "$ROOT/examples/cpw/postpro/autoresearch_cpw_mpi${RANKS}_cpu_${RUN_ID}"
set +e
(cd "$ROOT/examples/cpw" && env PALACE_SURFACE_PROFILE=1 "$MPIEXEC" -np "$RANKS" "$BIN" "$CFG") >"$LOG" 2>&1
rc=$?
set -e
echo "CHECK_LOG cpw_mpi${RANKS}_cpu $LOG rc=$rc"
if [[ $rc -ne 0 ]]; then
  tail -120 "$LOG" >&2 || true
  exit "$rc"
fi

python3 - "$LOG" "$RANKS" <<'PY'
import re, sys
from pathlib import Path
path=Path(sys.argv[1]); ranks=sys.argv[2]
lines=path.read_text(errors='replace').splitlines()
errors=[]
for needle in ['MFEM abort', 'Segmentation fault', 'MPI_ABORT', 'CUDA error', 'CUDA_ERROR']:
    if any(needle in line for line in lines):
        errors.append(f'found {needle}')

def last_float(pattern):
    val=None; rx=re.compile(pattern)
    for line in lines:
        m=rx.search(line)
        if m: val=float(m.group(1))
    return val

def group_sum(kind):
    total=0; rx=re.compile(rf'SurfaceFunctional profile kind={kind} groups=(\d+)')
    for line in lines:
        m=rx.search(line)
        if m: total += int(m.group(1))
    return total
metrics={
  f'cpw_mpi{ranks}_cpu_total': last_float(r'^Total\s+([0-9.]+)\s'),
  f'cpw_mpi{ranks}_cpu_postprocessing': last_float(r'^Postprocessing\s+([0-9.]+)\s'),
  f'cpw_mpi{ranks}_cpu_paraview': last_float(r'^\s+Paraview\s+([0-9.]+)\s') or 0.0,
  f'cpw_mpi{ranks}_cpu_gridfunction': last_float(r'^\s+Grid function\s+([0-9.]+)\s') or 0.0,
  f'cpw_mpi{ranks}_cpu_farfields': last_float(r'^\s+Far Fields\s+([0-9.]+)\s') or 0.0,
  f'cpw_mpi{ranks}_cpu_surface_flux_groups': group_sum('SURFACE_FLUX'),
  f'cpw_mpi{ranks}_cpu_interface_epr_groups': group_sum('INTERFACE_EPR'),
  f'cpw_mpi{ranks}_cpu_farfield_groups': group_sum('FARFIELD'),
}
if metrics[f'cpw_mpi{ranks}_cpu_postprocessing'] is None:
    errors.append('missing CPU MPI postprocessing timer')
for k,v in metrics.items():
    if v is None:
        continue
    if isinstance(v, int):
        print(f'METRIC {k}={v}')
    else:
        print(f'METRIC {k}={float(v):.6f}')
if errors:
    for e in errors:
        print(f'ERROR mpi_cpu_check: {e}', file=sys.stderr)
    sys.exit(3)
PY
