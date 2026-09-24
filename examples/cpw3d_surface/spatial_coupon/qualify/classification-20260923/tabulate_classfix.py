#!/usr/bin/env python3
"""Tabulate the classification-fix AMR runs: per cycle SA/MS/MA raw + self-consistent corrected and C, vs the r5nm p5 reference."""
import csv, json, os, sys
REF = {"SA": 9.913492750223e-05, "MS": 0.0002127369678677, "MA": 4.000217518061e-06, "C": 8.704808384636999e-14}
def read_csv(p):
    with open(p) as f:
        r = csv.reader(f); hdr = [h.strip() for h in next(r)]; row = [x.strip() for x in next(r)]
    return dict(zip(hdr, row))
def cycle(d):
    q = read_csv(os.path.join(d, "surface-Q-corrected.csv"))
    c = read_csv(os.path.join(d, "terminal-C.csv"))
    out = {}
    for i, k in ((1, "SA"), (2, "MS"), (3, "MA")):
        out[k] = {"raw": float(q[f"p_surf raw[{i}]"]), "sc": float(q[f"p_surf corrected[{i}]"])}
    e_raw = float(q["E_elec raw (J)"]); e_cor = float(q["E_elec corrected (J)"])
    out["C_raw"] = float(c["C[i][1] (F)"]); out["C_corr"] = out["C_raw"] * e_cor / e_raw
    out["spread"] = float(q["max trace closure spread"]); out["conf"] = int(float(q["confidence pass"]))
    return out
def run(root):
    cycles = {}
    for n in range(1, 11):
        d = os.path.join(root, "postpro", f"iteration{n:02d}")
        if os.path.isfile(os.path.join(d, "surface-Q-corrected.csv")): cycles[f"c{n}"] = cycle(d)
    log = open(os.path.join(root, "palace.log")).read()
    if "Completed 10 iterations" in log: cycles["c11"] = cycle(os.path.join(root, "postpro"))
    return cycles
def pct(v, k): return 100.0 * (v / REF[k] - 1.0)
if __name__ == "__main__":
    runs = {r: run(r) for r in sys.argv[1:]}
    for r, cs in runs.items():
        print(f"== {r}")
        for c, v in cs.items():
            print(f" {c:>3}: SA {pct(v['SA']['raw'],'SA'):+6.2f}/{pct(v['SA']['sc'],'SA'):+6.2f}  MS {pct(v['MS']['raw'],'MS'):+6.2f}/{pct(v['MS']['sc'],'MS'):+6.2f}  MA {pct(v['MA']['raw'],'MA'):+6.2f}/{pct(v['MA']['sc'],'MA'):+6.2f}  C {pct(v['C_raw'],'C'):+5.2f}/{pct(v['C_corr'],'C'):+5.2f}  conf {v['conf']} spread {v['spread']:.4f}")
    if len(runs) >= 2:
        a, b = list(runs.values())[:2]
        worst = 0.0
        for c in a:
            if c not in b: continue
            for k in ("SA", "MS", "MA"):
                for f in ("raw", "sc"):
                    x, y = a[c][k][f], b[c][k][f]; worst = max(worst, abs(x - y) / abs(y))
            for k in ("C_raw", "C_corr"):
                worst = max(worst, abs(a[c][k] - b[c][k]) / abs(b[c][k]))
        print(f"max relative difference over common cycles ({', '.join(c for c in a if c in b)}): {worst:.3e}")
    json.dump(runs, open("tabulation.json", "w"), indent=1)
