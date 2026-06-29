#!/usr/bin/env bash
set -euo pipefail
LATEST=${PALACE_AR_LATEST:-/home/ubuntu/scratch/palace_nsight_autoresearch/latest}
python3 - "$LATEST/output" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
required_main = ['E_real','E_imag','B_real','B_imag','U_e','U_m','S']
for cycle in ['Cycle000001', 'Cycle000002']:
    pvtu = root / 'paraview' / 'eigenmode' / cycle / 'data.pvtu'
    if not pvtu.exists():
        raise SystemExit(f'missing {pvtu}')
    txt = pvtu.read_text(errors='replace')
    missing = [name for name in required_main if f'Name="{name}"' not in txt]
    if missing:
        raise SystemExit(f'{pvtu} missing fields: {missing}')
final = root / 'paraview' / 'eigenmode' / 'Cycle000003' / 'data.pvtu'
if final.exists():
    txt = final.read_text(errors='replace')
    for name in ['Rank', 'Indicator']:
        if f'Name="{name}"' not in txt:
            raise SystemExit(f'{final} missing {name}')
size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
if size < 650_000_000:
    raise SystemExit(f'output unexpectedly small: {size}')
PY
