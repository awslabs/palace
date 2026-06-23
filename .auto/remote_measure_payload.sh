#!/usr/bin/env bash
set -euo pipefail

RUN_ID="${1:?run id required}"
CASE_SET="${2:-cpw}"
REMOTE_ROOT="/home/ubuntu/palace"
SPACK_SETUP="/home/ubuntu/spack/share/spack/setup-env.sh"
API="/tmp/libcuda_api_counter.so"

source "$SPACK_SETUP"
cd "$REMOTE_ROOT"

spack -e "$REMOTE_ROOT" install -j4 >/tmp/palace_ar_build_${RUN_ID}.log 2>&1
PREFIX="$(spack -e "$REMOTE_ROOT" location -i local.palace)"
BIN="$PREFIX/bin/palace-x86_64.bin"
if [[ ! -x "$BIN" ]]; then
  echo "ERROR missing Palace binary: $BIN" >&2
  tail -80 /tmp/palace_ar_build_${RUN_ID}.log >&2 || true
  exit 2
fi
if [[ ! -f "$API" ]]; then
  echo "ERROR missing CUDA API counter: $API" >&2
  exit 2
fi

python3 - "$RUN_ID" "$CASE_SET" <<'PY'
import json, re, sys
from pathlib import Path
run_id=sys.argv[1]
case_set=sys.argv[2]
root=Path('/home/ubuntu/palace')

def labels_for(case_set):
    if case_set == 'cpw':
        return [('cpw', 'examples/cpw/cpw_lumped_uniform.json', 'cpw')]
    if case_set == 'all':
        return [
            ('spheres', 'examples/spheres/spheres.json', 'static'),
            ('rings', 'examples/rings/rings.json', 'static'),
            ('cpw', 'examples/cpw/cpw_lumped_uniform.json', 'cpw'),
        ]
    raise SystemExit(f'unknown PALACE_AR_CASES={case_set!r}; expected cpw or all')

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

def write_case(label, rel, kind):
    j=load_jsonc(root/rel)
    j.setdefault('Problem', {})['Output']=f'postpro/autoresearch_{label}_{run_id}'
    j['Problem'].setdefault('OutputFormats', {})['Paraview']=True
    j['Problem']['OutputFormats']['GridFunction']=True
    j.setdefault('Solver', {})['Device']='GPU'
    if kind == 'cpw':
        j['Solver'].setdefault('Driven', {})['Samples']=[{'Type':'Point','Freq':[17.0],'SaveStep':1}]
        j['Solver']['Driven']['Save']=[17.0]
    out=Path(f'/tmp/palace_ar_{label}_{run_id}.json')
    out.write_text(json.dumps(j, indent=2))

for label, rel, kind in labels_for(case_set):
    write_case(label, rel, kind)
PY

run_case() {
  local label="$1"
  local dir="$2"
  local cfg="/tmp/palace_ar_${label}_${RUN_ID}.json"
  local log="/tmp/palace_ar_${label}_${RUN_ID}.log"
  rm -rf "$dir/postpro/autoresearch_${label}_${RUN_ID}"
  set +e
  (cd "$dir" && env PALACE_SURFACE_PROFILE=1 LD_PRELOAD="$API" "$BIN" "$cfg") >"$log" 2>&1
  local rc=$?
  set -e
  echo "CASE_LOG $label $log rc=$rc"
  if [[ $rc -ne 0 ]]; then
    tail -120 "$log" >&2 || true
    return $rc
  fi
}

rc=0
case "$CASE_SET" in
  cpw)
    run_case cpw "$REMOTE_ROOT/examples/cpw" || rc=$?
    ;;
  all)
    run_case spheres "$REMOTE_ROOT/examples/spheres" || rc=$?
    run_case rings "$REMOTE_ROOT/examples/rings" || rc=$?
    run_case cpw "$REMOTE_ROOT/examples/cpw" || rc=$?
    ;;
  *)
    echo "ERROR unknown PALACE_AR_CASES=$CASE_SET" >&2
    exit 2
    ;;
esac

python3 - "$RUN_ID" "$rc" "$CASE_SET" <<'PY'
import re, sys
from pathlib import Path
run_id=sys.argv[1]
rc=int(sys.argv[2])
case_set=sys.argv[3]
labels=['cpw'] if case_set == 'cpw' else ['spheres','rings','cpw']
metrics={}
errors=[]

def last_float(pattern, lines):
    val=None
    rx=re.compile(pattern)
    for line in lines:
        m=rx.search(line)
        if m: val=float(m.group(1))
    return val

def last_int(pattern, lines):
    val=None
    rx=re.compile(pattern)
    for line in lines:
        m=rx.search(line)
        if m: val=int(m.group(1))
    return val

def group_sum(kind, lines):
    total=0
    rx=re.compile(rf'SurfaceFunctional profile kind={kind} groups=(\d+)')
    for line in lines:
        m=rx.search(line)
        if m: total += int(m.group(1))
    return total

for label in labels:
    path=Path(f'/tmp/palace_ar_{label}_{run_id}.log')
    if not path.exists():
        errors.append(f'missing log {path}')
        continue
    lines=path.read_text(errors='replace').splitlines()
    for needle in ['MFEM abort', 'CUDA error', 'CUDA_ERROR']:
        if any(needle in line for line in lines):
            errors.append(f'{label}: found {needle}')
    metrics[f'{label}_total']=last_float(r'^Total\s+([0-9.]+)\s', lines)
    metrics[f'{label}_postprocessing']=last_float(r'^Postprocessing\s+([0-9.]+)\s', lines)
    metrics[f'{label}_paraview']=last_float(r'^\s+Paraview\s+([0-9.]+)\s', lines) or 0.0
    metrics[f'{label}_gridfunction']=last_float(r'^\s+Grid function\s+([0-9.]+)\s', lines) or 0.0
    metrics[f'{label}_diskio']=last_float(r'^Disk IO\s+([0-9.]+)\s', lines) or 0.0
    metrics[f'{label}_operator_construction']=last_float(r'^Operator Construction\s+([0-9.]+)\s', lines) or 0.0
    metrics[f'{label}_farfields']=last_float(r'^\s+Far Fields\s+([0-9.]+)\s', lines) or 0.0
    metrics[f'{label}_nvrtc']=last_int(r'nvrtcCompileProgram=(\d+)', lines) or 0
    metrics[f'{label}_cuModuleLoadData']=last_int(r'cuModuleLoadData=(\d+)', lines) or 0
    metrics[f'{label}_surface_flux_groups']=group_sum('SURFACE_FLUX', lines)
    metrics[f'{label}_interface_epr_groups']=group_sum('INTERFACE_EPR', lines)
    metrics[f'{label}_farfield_groups']=group_sum('FARFIELD', lines)

needed=[f'{label}_postprocessing' for label in labels]
for k in needed:
    if metrics.get(k) is None:
        errors.append(f'missing metric {k}')
if not errors:
    metrics['surface_score_seconds']=sum(float(metrics[k]) for k in needed)

if 'surface_score_seconds' in metrics:
    print(f"METRIC surface_score_seconds={metrics['surface_score_seconds']:.6f}")
for k in sorted(metrics):
    if k == 'surface_score_seconds':
        continue
    v=metrics[k]
    if v is None:
        continue
    if isinstance(v, int):
        print(f'METRIC {k}={v}')
    else:
        print(f'METRIC {k}={float(v):.6f}')

if errors:
    for e in errors:
        print(f'ERROR {e}', file=sys.stderr)
    sys.exit(3)
if rc != 0:
    print(f'ERROR one or more cases failed rc={rc}', file=sys.stderr)
    sys.exit(rc)
PY
