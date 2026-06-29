#!/usr/bin/env bash
set -euo pipefail
LATEST=${PALACE_AR_LATEST:-/home/ubuntu/scratch/palace_nsight_autoresearch/latest}
if [[ -n "${PALACE_AR_PROFILE_EXAMPLE:-}" || -f "$LATEST/profile_example.txt" ]]; then
python3 - "$LATEST/output" "$LATEST/profile_example.txt" <<'PY'
from pathlib import Path
import sys
root = Path(sys.argv[1])
marker = Path(sys.argv[2])
if marker.exists():
    print(f'checking profiled example output for {marker.read_text().strip()}')
if not root.exists():
    raise SystemExit(f'missing output directory: {root}')
size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
if size < 50_000:
    raise SystemExit(f'profiled example output unexpectedly small: {size}')
if not any(root.rglob('*.csv')) and not any(root.rglob('*.pvd')) and not any(root.rglob('*.pvtu')) and not any(root.rglob('*.vtu')):
    raise SystemExit(f'profiled example produced no CSV or ParaView files under {root}')
PY
exit 0
fi
python3 - "$LATEST/output" <<'PY'
from pathlib import Path
import re
import sys
root = Path(sys.argv[1])

def tag_block(txt, tag):
    m = re.search(rf'<{tag}\b[^>]*>(.*?)</{tag}>', txt, flags=re.S)
    return m.group(1) if m else None

required_main = ['E_real','E_imag','B_real','B_imag','U_e','U_m','S']
for cycle in ['Cycle000001', 'Cycle000002']:
    pvtu = root / 'paraview' / 'eigenmode' / cycle / 'data.pvtu'
    if not pvtu.exists():
        raise SystemExit(f'missing {pvtu}')
    txt = pvtu.read_text(errors='replace')
    point_data = tag_block(txt, 'PPointData')
    if point_data is None:
        raise SystemExit(f'{pvtu} missing PPointData block')
    missing = [name for name in required_main if f'Name="{name}"' not in point_data]
    if missing:
        raise SystemExit(f'{pvtu} PPointData missing fields: {missing}')

final = root / 'paraview' / 'eigenmode' / 'Cycle000003' / 'data.pvtu'
if final.exists():
    txt = final.read_text(errors='replace')
    point_data = tag_block(txt, 'PPointData')
    cell_data = tag_block(txt, 'PCellData')
    if point_data is None or cell_data is None:
        raise SystemExit(f'{final} missing PPointData/PCellData blocks')
    for name in ['Rank', 'Indicator']:
        if f'Name="{name}"' in point_data:
            raise SystemExit(f'{final} incorrectly writes {name} as PPointData')
        if f'Name="{name}"' not in cell_data:
            raise SystemExit(f'{final} PCellData missing {name}')
    pieces = sorted(final.parent.glob('proc*.vtu'))
    if not pieces:
        raise SystemExit(f'{final.parent} has no VTU pieces')
    for piece in pieces:
        ptxt = piece.read_text(errors='replace')
        point_data = tag_block(ptxt, 'PointData')
        cell_data = tag_block(ptxt, 'CellData')
        if point_data is None or cell_data is None:
            raise SystemExit(f'{piece} missing PointData/CellData blocks')
        for name in ['Rank', 'Indicator']:
            if f'Name="{name}"' in point_data:
                raise SystemExit(f'{piece} incorrectly writes {name} as PointData')
            if f'Name="{name}"' not in cell_data:
                raise SystemExit(f'{piece} CellData missing {name}')

size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
if size < 650_000_000:
    raise SystemExit(f'output unexpectedly small: {size}')
PY
