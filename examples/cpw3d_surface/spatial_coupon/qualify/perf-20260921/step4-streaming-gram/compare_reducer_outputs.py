#!/usr/bin/env python3
"""Compare two reducer output directories' surface-/domain-response-matrix.csv per key."""
import csv, json, sys
from pathlib import Path
def load(path):
    rows = [[c.strip() for c in r] for r in csv.reader(open(path, newline=""))]
    header, body = rows[0], rows[1:]
    keys = [i for i, h in enumerate(header) if h in ("interface", "edge", "R (m)", "basis_i", "basis_j")]
    vals = [i for i in range(len(header)) if i not in keys]
    return {tuple(r[i] for i in keys): [r[i] for i in vals] for r in body}, [header[i] for i in vals], [tuple(r[i] for i in keys) for r in body]
out = {}
for name in ("domain-response-matrix.csv", "surface-response-matrix.csv"):
    a, va, oa = load(Path(sys.argv[1]) / name); b, vb, ob = load(Path(sys.argv[2]) / name)
    assert set(a) == set(b) and va == vb, name
    same = total = 0; worst = (0.0, None); largest = 0.0
    for k in a:
        for col, x, y in zip(va, a[k], b[k]):
            total += 1; fx, fy = float(x), float(y); largest = max(largest, abs(fx), abs(fy))
            if x == y: same += 1
            else:
                rel = abs(fx - fy) / max(abs(fx), abs(fy))
                if rel > worst[0]: worst = (rel, k + (col, x, y))
    # relative to the largest entry of the same column family as well
    out[name] = {"Rows": len(a), "Entries": total, "BitIdenticalEntries": same, "MaxRelativeDifferencePerEntry": worst[0],
                 "MaxRelativeDifferenceAt": worst[1], "RowOrderIdentical": oa == ob}
    print(name, json.dumps(out[name]))
if len(sys.argv) > 3:
    Path(sys.argv[3]).write_text(json.dumps(out, indent=2) + "\n")
