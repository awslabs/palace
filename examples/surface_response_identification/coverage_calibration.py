#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Calibration of the coverage report's predicted error bands (readiness STAGE D; DESIGN.md section 5).

Reads a RECORDS file (JSON) naming, per solver family (ES / Maxwell), the graded window cells of record
(the cards' headline discrepancies per Type: SA / MS / MA_sharp / C) with each window's run postpro and
preflight manifest, extracts the predictors with coverage_report.py (the RAW share, the class shares by
LENGTH, the correction ratio g, the last AMR step), fits the three pre-registered bias models per Type

    M0  b = beta0                                   (the median discrepancy)
    M1  b = beta0 + beta1 x (last ft step, pt)      (the depth lever)
    M2  b = beta0 + beta1 x SEC share + beta2 x corner share   (the class-mix lever; LENGTH shares)

by least squares (closed-form normal equations, pure Python, bitwise reproducible), takes H_t = the 0.90
empirical quantile of the in-sample residuals, validates LEAVE-ONE-WINDOW-OUT (USER 1002 (d): coverage
>= 90 %), SELECTS per Type the model with the smallest LOO 90 % half-width among those reaching 90 %
(else M0 WIDENED by the smallest factor k that reaches it), and writes the band file coverage_report.py
reads, with the LOO coverage, its binomial (Wilson) 90 % CI and n per Type, the class census (n_c), the
median g per Type, every cell with its predictors and in/out-of-band verdict, and the sha256 of every input
(decision 1049: the calibration set is the DEVELOPMENT set; the LOO coverage is an in-sample-family estimate).

Records file schema:
{"families": {"ES": {"calibrate": [CELL, ...], "test": [CELL, ...]}, "Maxwell": {...}}}
CELL = {"window": str, "postpro": dir, "manifest": file, "cells": {"SA": pct, "MS": pct, "MA_sharp": pct, "C": pct},
        "ladder_run": dir (optional), "reference_order": int (optional), "ma_reading": "sharp" (optional),
        "flags": [str] (information), "note": str}
"""

import argparse
import datetime
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_report as cr  # noqa: E402

TARGET_COVERAGE = 0.90
RESIDUAL_QUANTILE = 0.90
Z90 = 1.6448536269514722  # the two-sided 90 % normal quantile (the Wilson interval)
CARD_TYPES = ("SA", "MS", "MA_sharp", "C")
BAND_TYPE = {"SA": "SA", "MS": "MS", "MA_sharp": "MA", "C": "C"}
MODELS = ("M0", "M1", "M2")


def empirical_quantile(values, q):
    """The q-quantile by sorted order: the ceil(q n)-th smallest (the 0.90 of 10 values = the 9th)."""
    if not values:
        return 0.0
    s = sorted(values)
    k = max(1, min(len(s), int(math.ceil(q * len(s) - 1e-12))))
    return s[k - 1]


def median(values):
    if not values:
        return None
    s = sorted(values)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def wilson_ci(successes, n, z=Z90):
    """The Wilson score interval for a binomial proportion (pure arithmetic; a 90 % CI with z = 1.645)."""
    if n == 0:
        return [0.0, 1.0]
    p = successes / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    return [max(0.0, centre - half), min(1.0, centre + half)]


def solve_normal_equations(rows, y):
    """Least squares beta for rows (list of predictor lists incl. the 1 column) by Gaussian elimination
    with partial pivoting on X^T X beta = X^T y; a singular system falls back to the mean (M0)."""
    p = len(rows[0])
    A = [[0.0] * p for _ in range(p)]
    b = [0.0] * p
    for x, yi in zip(rows, y):
        for i in range(p):
            b[i] += x[i] * yi
            for j in range(p):
                A[i][j] += x[i] * x[j]
    for col in range(p):
        pivot = max(range(col, p), key=lambda r: abs(A[r][col]))
        if abs(A[pivot][col]) < 1e-14 * (1.0 + abs(A[col][col])):
            return None
        A[col], A[pivot] = A[pivot], A[col]
        b[col], b[pivot] = b[pivot], b[col]
        for r in range(p):
            if r == col:
                continue
            f = A[r][col] / A[col][col]
            for c in range(col, p):
                A[r][c] -= f * A[col][c]
            b[r] -= f * b[col]
    return [b[i] / A[i][i] for i in range(p)]


def design_row(model, pred):
    if model == "M0":
        return [1.0]
    if model == "M1":
        step = pred.get("last_step_pt")
        return [1.0, step if step is not None else 0.0]
    return [1.0, pred.get("sec_share", 0.0), pred.get("corner_share", 0.0)]


def fit_model(model, cells, type_key):
    """-> {"coefficients", "residuals", "H_pt"} on the cells (each with "pred" and "d" per type_key)."""
    rows = [design_row(model, c["pred"][type_key]) for c in cells]
    y = [c["d"][type_key] for c in cells]
    beta = None
    if model == "M0":
        beta = [median(y)]
    else:
        beta = solve_normal_equations(rows, y)
        if beta is None:
            beta = [median(y)] + [0.0] * (len(rows[0]) - 1)
    predicted = [sum(bi * xi for bi, xi in zip(beta, x)) for x in rows]
    residuals = [yi - pi for yi, pi in zip(y, predicted)]
    names = ("beta0", "beta1", "beta2")
    return {"coefficients": {names[i]: beta[i] for i in range(len(beta))}, "residuals": residuals,
            "H_pt": empirical_quantile([abs(r) for r in residuals], RESIDUAL_QUANTILE),
            "median_residual_pt": median([abs(r) for r in residuals])}


def loo_evaluate(model, cells, type_key):
    """Leave-one-window-out: for each cell the fit on the others, the prediction, H_(-w); covered iff
    |d_w - b_(-w)(x_w)| <= H_(-w). -> coverage, the LOO 90 % half-width, the per-cell records, the widening k."""
    records = []
    for i, cell in enumerate(cells):
        others = cells[:i] + cells[i + 1:]
        fit = fit_model(model, others, type_key)
        beta = [fit["coefficients"].get(k, 0.0) for k in ("beta0", "beta1", "beta2")]
        x = design_row(model, cell["pred"][type_key])
        predicted = sum(bi * xi for bi, xi in zip(beta, x))
        resid = cell["d"][type_key] - predicted
        H = fit["H_pt"]
        records.append({"window": cell["window"], "d_pt": cell["d"][type_key], "predicted_pt": predicted,
                        "loo_residual_pt": resid, "H_loo_pt": H, "covered": abs(resid) <= H,
                        "ratio": (abs(resid) / H) if H > 0 else (0.0 if resid == 0.0 else float("inf"))})
    n = len(records)
    covered = sum(1 for r in records if r["covered"])
    half_width = empirical_quantile([abs(r["loo_residual_pt"]) for r in records], RESIDUAL_QUANTILE)
    # the smallest factor k with LOO coverage >= 90 %: the ceil(0.9 n)-th smallest residual / H_(-w) ratio
    k = max(1.0, empirical_quantile([r["ratio"] for r in records], TARGET_COVERAGE)) if n else 1.0
    return {"coverage": covered / n if n else None, "covered": covered, "n": n,
            "ci90": wilson_ci(covered, n) if n else None, "loo_half_width_pt": half_width, "widen_k": k,
            "cells": records}


def select_model(cells, type_key):
    """DESIGN 5.2 SELECTION: the smallest LOO 90 % half-width among the models with LOO coverage >= 90 %;
    none -> M0 widened by k. -> the band entry for this Type plus every candidate's report."""
    candidates = {}
    for model in MODELS:
        if model != "M0" and len(cells) < len(design_row(model, cells[0]["pred"][type_key])) + 2:
            continue
        fit = fit_model(model, cells, type_key)
        loo = loo_evaluate(model, cells, type_key)
        candidates[model] = {"fit": fit, "loo": loo}
    reaching = [m for m, c in candidates.items() if (c["loo"]["coverage"] or 0.0) >= TARGET_COVERAGE]
    if reaching:
        chosen = min(reaching, key=lambda m: (candidates[m]["loo"]["loo_half_width_pt"], MODELS.index(m)))
        widened = 1.0
    else:
        chosen = "M0"
        widened = candidates["M0"]["loo"]["widen_k"]
    c = candidates[chosen]
    loo = c["loo"]
    if widened > 1.0:
        covered = sum(1 for r in loo["cells"] if abs(r["loo_residual_pt"]) <= widened * r["H_loo_pt"])
        loo_reported = dict(loo, coverage=covered / loo["n"], covered=covered, ci90=wilson_ci(covered, loo["n"]),
                            note=f"WIDENED x {widened:.3f}")
    else:
        loo_reported = loo
    entry = {"model": chosen, "coefficients": c["fit"]["coefficients"], "H_pt": c["fit"]["H_pt"],
             "median_residual_pt": c["fit"]["median_residual_pt"], "widened": widened, "n": len(cells),
             "loo_coverage": {k: v for k, v in loo_reported.items() if k != "cells"},
             "candidates": {m: {"H_pt": cc["fit"]["H_pt"], "coefficients": cc["fit"]["coefficients"],
                                "loo_coverage": cc["loo"]["coverage"], "loo_half_width_pt": cc["loo"]["loo_half_width_pt"],
                                "ci90": cc["loo"]["ci90"]} for m, cc in candidates.items()},
             "loo_cells": loo_reported["cells"]}
    return entry


def cell_predictors(cell, bands=None):
    """Run coverage_report.build_report on the cell's run; -> the report and the per-Type predictors."""
    report = cr.build_report(cell["postpro"], cell["manifest"], bands, cell.get("discovery"), cell.get("ladder_run"),
                             cell.get("reference_order"), cell.get("ma_reading", "sharp"), cell.get("window"))
    pred = {}
    for card_type in CARD_TYPES:
        t = BAND_TYPE[card_type]
        e = report["block_A"]["energy"].get(t, {})
        L = report["block_A"]["length"].get(t, {})
        lad = report["block_D"]["ladder"]["types"].get(t, {})
        p = cr.class_predictors(L)
        p.update({"last_step_pt": lad.get("last_step_pt"), "raw_share_energy": e.get("raw_share_total"),
                  "raw_share_length": L.get("raw_share_total"), "g": e.get("g"), "exact_share_energy": e.get("exact_share"),
                  "unclaimed_share": e.get("unclaimed_share"), "classes_share_length": L.get("classes_share", {}),
                  "classes_share_energy": e.get("classes_share", {}), "energy_mode": report["energy_mode"]})
        if card_type == "C":
            # C has no interface Type: the SA predictors stand in for the depth lever only
            sa = report["block_D"]["ladder"]["types"].get("SA", {})
            p = {"last_step_pt": sa.get("last_step_pt"), "sec_share": 0.0, "corner_share": 0.0}
        pred[card_type] = p
    return report, pred


def class_census(cells):
    """n_c = the number of calibration windows with the class share > 5 % (LENGTH) on any Type."""
    census = {}
    for cell in cells:
        seen = set()
        for card_type in ("SA", "MS", "MA_sharp"):
            for cls, share in cell["pred"][card_type].get("classes_share_length", {}).items():
                if share > cr.CLASS_MATERIAL_SHARE:
                    seen.add(cls)
        for cls in seen:
            census[cls] = census.get(cls, 0) + 1
    return census


DEV_SET_CAVEAT = ("the calibration set is the DEVELOPMENT set (gap G8): the LOO coverage is an in-sample-family "
                  "estimate; the out-of-family test is stage E's blind windows (decision 1049)")


def calibrate_family(family, spec, log):
    """Fit one family (ES / Maxwell) from its records: -> (band entry per Type, the loaded cells)."""
    cells = []
    for rec in spec.get("calibrate", []):
        report, pred = cell_predictors(rec)
        d = {ct: rec["cells"][ct] for ct in CARD_TYPES if rec["cells"].get(ct) is not None}
        cells.append({"window": rec["window"], "pred": pred, "d": d, "report": report, "record": rec})
        log(f"  {family} {rec['window']}: {report['energy_mode']} raw SA {cr.fmt_pct(pred['SA']['raw_share_energy'])} "
            f"g SA {pred['SA']['g']} last step SA {cr.fmt_pt(pred['SA']['last_step_pt'])}")
    types = {}
    for card_type in CARD_TYPES:
        usable = [c for c in cells if card_type in c["d"]]
        if len(usable) < 3:
            continue
        entry = select_model(usable, card_type)
        gs = [c["pred"][card_type].get("g") for c in usable if c["pred"][card_type].get("g") is not None
              and c["pred"][card_type].get("energy_mode") == "ENERGY"]
        entry["median_g"] = median(gs) if gs else 0.0
        entry["n_energy_mode"] = len(gs)
        types[BAND_TYPE[card_type]] = entry
        log(f"  {family} {card_type}: {entry['model']} H {entry['H_pt']:.3f} pt widened {entry['widened']:.3f} "
            f"LOO {entry['loo_coverage']['coverage']:.3f} CI {entry['loo_coverage']['ci90']} n {entry['n']}")
    return {"types": types, "class_census": class_census(cells), "n": len(cells)}, cells


def evaluate_with_bands(family, cells, test_specs, bands, log):
    """Second pass: every calibration + test cell re-reported WITH the band file; in / out of band per Type with
    the fired named terms (the reason), the energy-vs-length share agreement."""
    out = []
    for group, items in (("calibrate", cells), ("test", test_specs)):
        for item in items:
            rec = item["record"] if group == "calibrate" else item
            report, pred = cell_predictors(rec, bands)
            entry = {"window": rec["window"], "group": group, "label": report["label"]["label"], "energy_mode": report["energy_mode"],
                     "types": {}, "fired": [r["term"] for r in report["block_B"] if r["fired"]],
                     "postpro": rec["postpro"], "manifest": rec["manifest"], "note": rec.get("note")}
            for card_type in CARD_TYPES:
                d = (rec.get("cells") or {}).get(card_type)
                band = (report.get("block_C") or {}).get(BAND_TYPE[card_type])
                if d is None or not band or not band.get("calibrated"):
                    continue
                inside = band["low_pt"] <= d <= band["high_pt"]
                entry["types"][card_type] = {"d_pt": d, "band": [band["low_pt"], band["high_pt"]], "in_band": inside,
                                             "bias_pt": band.get("bias_pt"), "gap_pt": band.get("gap_pt"),
                                             "extensions": band.get("extensions"), "predictors": pred[card_type]}
                if not inside:
                    log(f"  OUT {family} {rec['window']} {card_type}: {d:+.2f} vs [{band['low_pt']:+.2f}, {band['high_pt']:+.2f}] "
                        f"fired {entry['fired']}")
            out.append(entry)
    return out


def share_agreement(cells):
    """ENERGY-mode windows: per Type the (energy share - length share) of the Exact row and of the SEC classes."""
    rows = []
    for c in cells:
        for card_type in ("SA", "MS", "MA_sharp"):
            p = c["pred"][card_type]
            if p.get("energy_mode") != "ENERGY":
                continue
            e_sec = sum(v for k, v in p.get("classes_share_energy", {}).items() if k.startswith("SpatialEdgeCluster") or k == "JunctionLead")
            rows.append({"window": c["window"], "type": BAND_TYPE[card_type], "raw_energy": p.get("raw_share_energy"),
                         "raw_length": p.get("raw_share_length"), "sec_energy": e_sec, "sec_length": p.get("sec_share")})
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--records", required=True, help="the records JSON (the cells of record with their runs)")
    parser.add_argument("--output", required=True, help="the band file to write (coverage-bands-<tag>.json)")
    parser.add_argument("--tag", default="v1")
    parser.add_argument("--fit-date", help="ISO date stamped in the file (default: today, UTC)")
    parser.add_argument("--lean-output", help="also write the band file WITHOUT the per-cell evaluation (the copy "
                        "committed as coverage_bands_v1.json: the fitted numbers + the provenance shas)")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    log = (lambda s: None) if args.quiet else (lambda s: print(s, flush=True))
    records = cr.read_json(args.records)
    inputs = [{"path": os.path.abspath(args.records), "sha256": cr.sha256_file(args.records)}]
    for fam in records.get("families", {}).values():
        for rec in fam.get("calibrate", []) + fam.get("test", []):
            for key in ("manifest",):
                inputs.append({"path": os.path.abspath(rec[key]), "sha256": cr.sha256_file(rec[key])})
            for name in ("surface-Q-corrected.csv", "surface-Q-edge.csv", "surface-response-model-energy.csv",
                         "surface-response-model-raw-energy.csv"):
                p = os.path.join(rec["postpro"], name)
                if os.path.isfile(p):
                    inputs.append({"path": os.path.abspath(p), "sha256": cr.sha256_file(p)})
    for src in records.get("cards", []):
        inputs.append({"path": src["path"], "sha256": cr.sha256_file(src["path"]) if os.path.isfile(src["path"]) else None,
                       "decision": src.get("decision")})
    bands = {"version": args.tag, "tool": "coverage_calibration.py", "fit_date": args.fit_date or datetime.date.today().isoformat(),
             "target_coverage": TARGET_COVERAGE, "residual_quantile": RESIDUAL_QUANTILE, "dev_set_caveat": DEV_SET_CAVEAT,
             "families": {}, "provenance": {"inputs": inputs, "n_per_family": {}}}
    loaded = {}
    for family, spec in records.get("families", {}).items():
        log(f"family {family}: {len(spec.get('calibrate', []))} calibration + {len(spec.get('test', []))} test cells")
        fam_bands, cells = calibrate_family(family, spec, log)
        bands["families"][family] = fam_bands
        bands["provenance"]["n_per_family"][family] = fam_bands["n"]
        loaded[family] = cells
    bands["_path"] = args.output
    bands["_sha256"] = ""
    evaluation = {}
    for family, spec in records.get("families", {}).items():
        evaluation[family] = {"cells": evaluate_with_bands(family, loaded[family], spec.get("test", []), bands, log),
                              "share_agreement": share_agreement(loaded[family])}
        for t, entry in bands["families"][family]["types"].items():
            in_band = [c["types"][ct]["in_band"] for c in evaluation[family]["cells"] if c["group"] == "calibrate"
                       for ct in c["types"] if BAND_TYPE[ct] == t]
            entry["full_band_coverage_calibrate"] = (sum(in_band) / len(in_band)) if in_band else None
    bands.pop("_path")
    bands.pop("_sha256")
    bands["evaluation"] = evaluation
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(bands, f, indent=1, sort_keys=True)
        f.write("\n")
    log(f"wrote {args.output}")
    if args.lean_output:
        lean = json.loads(json.dumps(bands))
        lean.pop("evaluation")
        for fam in lean["families"].values():
            for entry in fam["types"].values():
                entry.pop("loo_cells", None)
        lean["full_record"] = {"path": os.path.abspath(args.output), "sha256": cr.sha256_file(args.output)}
        with open(args.lean_output, "w") as f:
            json.dump(lean, f, indent=1, sort_keys=True)
            f.write("\n")
        log(f"wrote {args.lean_output} (lean)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
