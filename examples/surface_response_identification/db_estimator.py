# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The DomainBoundary local-density estimator g_db (boundary-cut DESIGN 2.1; decisions 442 /
454): the pre-solve bracket of the within-R energy of the DomainBoundary-excluded cells and the
re-reading of old records. After F-DB-a the runtime MEASURES that energy
(surface-response-domain-boundary-energy.csv), so the estimator serves as the expectation a
record is checked against, never as the verdict term.

    python3 db_estimator.py --run RUN_DIR --radius R_METRES [--length-unit L0_METRES]
                            [--output out.json]

RUN_DIR holds postpro/palace.json (the operator record: SurfaceResponse.Diagnostics.
DomainBoundaryExclusions.RawPortions = the excluded cells' raw claims, one per GEOMETRIC cell:
the co-located first-order split patches of one cell are already one portion),
postpro/surface-response-patches.csv and postpro/surface-response-patch-energy.csv (PatchEnergy:
true; the applied cells' ft energies) and the surface-Q-corrected.csv (the window ft per
interface).

Rule (DESIGN 2.1): g_db = sum over geometric cells of cell length x the ft density of the
nearest applied cells of the SAME feature within 3 R along the edge (the three nearest), times
the configuration factor [1, 2.5] (perpendicular / obtuse 1-1.5, a convex acute spike up to
2.5: the measured 1.45-2.2x), divided by the window's ft. No applied cell of the same feature
within 10 R -> UNBRACKETED (never a foreign-feature density, never the window mean). The
measured share (the DB CSV) is printed beside the bracket when present.
"""

import argparse
import csv
import json
import math
import os
import sys

NEAREST_CELLS = 3
CONFIGURATION_FACTOR = (1.0, 2.5)
SEARCH_OVER_R = 3.0
UNBRACKETED_OVER_R = 10.0


def _rows(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle)
        header = [h.strip() for h in next(reader)]
        rows = []
        for row in reader:
            if not row or not "".join(row).strip():
                continue
            rows.append([x.strip() for x in row])
    return header, rows


def _patch_energies(run):
    """The ft (evaluation 0, first source) energies per applied patch, by patch index."""
    header, rows = _rows(os.path.join(run, "postpro", "surface-response-patch-energy.csv"))
    column = {h: i for i, h in enumerate(header)}
    energy_columns = {h: i for i, h in enumerate(header) if h.startswith("fabricated surface energy[")}
    sources = sorted({int(float(r[column["source"]])) for r in rows})
    out = {}
    for r in rows:
        if int(float(r[column["source"]])) != sources[0] or int(float(r[column["evaluation"]])) != 0:
            continue
        patch = int(round(float(r[column["patch"]]))) - 1
        out[patch] = {
            "feature": int(round(float(r[column["feature"]]))),
            "origin": [float(r[column[f"origin {a} (m)"]]) for a in "xyz"],
            "cell": (float(r[column["cell begin (m)"]]), float(r[column["cell end (m)"]])),
            "energy": {h.split("[")[1].rstrip("] (J)"): float(r[i]) for h, i in energy_columns.items()},
        }
    return out, sources[0]


def _window_ft(run, source):
    header, rows = _rows(os.path.join(run, "postpro", "surface-Q-corrected.csv"))
    column = {h: i for i, h in enumerate(header)}
    ft = {}
    for h, i in column.items():
        if h.startswith("E_surf postprocessed fixed-trace["):
            interface = h.split("[")[1].split("]")[0]
            for r in rows:
                if int(round(float(r[0]))) == source:
                    ft[interface] = float(r[i])
    return ft


def _measured_share(run, source):
    path = os.path.join(run, "postpro", "surface-response-domain-boundary-energy.csv")
    if not os.path.isfile(path):
        return None
    header, rows = _rows(path)
    column = {h: i for i, h in enumerate(header)}
    shares = {}
    for r in rows:
        if int(r[column["source"]]) == source and int(r[column["evaluation"]]) == 0 and r[column["type"]] == "Total":
            for h, i in column.items():
                if h.startswith("domainboundary share["):
                    shares[h.split("[")[1].rstrip("]")] = float(r[i])
    return shares


def estimate(run, radius, length_unit):
    """radius: the matching radius R in metres; length_unit: metres per record length unit
    (the configuration's Model.L0; the record's RawPortions are in mesh length units, the
    patch-energy CSV in metres)."""
    with open(os.path.join(run, "postpro", "palace.json")) as handle:
        record = json.load(handle)
    diagnostics = record["SurfaceResponse"]["Diagnostics"]
    exclusions = diagnostics["DomainBoundaryExclusions"]
    portions = exclusions.get("RawPortions", {}).get("Portions", [])
    scale = 1.0  # the patch-energy CSV is in metres already
    energies, source = _patch_energies(run)
    window_ft = _window_ft(run, source)
    interfaces = sorted(window_ft)
    # The applied cells' densities per feature: energy / cell length at the cell centre.
    by_feature = {}
    for patch, data in energies.items():
        length = (data["cell"][1] - data["cell"][0]) * scale
        if length <= 0.0:
            continue
        by_feature.setdefault(data["feature"], []).append(
            {"centre": [c * scale for c in data["origin"]], "density": {t: e / length for t, e in data["energy"].items()}}
        )
    cells = []
    bracket = {t: [0.0, 0.0] for t in interfaces}
    unbracketed = []
    for portion in portions:
        length = float(portion["Length"]) * length_unit
        centre = [0.5 * (a + b) * length_unit for a, b in zip(portion["P0"], portion["P1"])]
        feature = portion["Feature"]
        pool = by_feature.get(feature, [])
        distances = sorted((math.dist(centre, c["centre"]), c) for c in pool)
        near = [c for d, c in distances if d <= SEARCH_OVER_R * radius][:NEAREST_CELLS]
        if not near:
            near = [c for d, c in distances if d <= UNBRACKETED_OVER_R * radius][:NEAREST_CELLS]
        entry = {"Feature": feature, "Type": portion["Type"], "Length": length, "Patches": portion["Patches"]}
        if not near:
            entry["Status"] = "UNBRACKETED"
            unbracketed.append(entry)
            cells.append(entry)
            continue
        density = {}
        for t in interfaces:
            values = [c["density"].get(t, 0.0) for c in near]
            density[t] = sum(values) / len(values)
        entry["Status"] = "bracketed"
        entry["NearestAppliedDistance"] = distances[0][0] if distances else None
        entry["LocalDensity"] = density
        entry["Energy"] = {t: [density[t] * length * f for f in CONFIGURATION_FACTOR] for t in interfaces}
        for t in interfaces:
            bracket[t][0] += entry["Energy"][t][0]
            bracket[t][1] += entry["Energy"][t][1]
        cells.append(entry)
    share = {t: [bracket[t][k] / window_ft[t] if window_ft[t] else float("nan") for k in (0, 1)] for t in interfaces}
    return {
        "Run": run,
        "Source": source,
        "MatchingRadius": radius,
        "GeometricCells": len(portions),
        "Bracketed": len(cells) - len(unbracketed),
        "Unbracketed": len(unbracketed),
        "Status": "UNBRACKETED" if unbracketed else ("bracketed" if cells else "no DomainBoundary cells"),
        "ConfigurationFactor": list(CONFIGURATION_FACTOR),
        "Rule": "boundary-cut DESIGN 2.1: cell length x the ft density of the nearest applied cells of the SAME feature within 3 R, x [1, 2.5]; no own-feature cell within 10 R -> UNBRACKETED",
        "EnergyBracket": bracket,
        "ShareBracket": share,
        "MeasuredShare": _measured_share(run, source),
        "Cells": cells,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", required=True, help="a run directory with postpro/")
    parser.add_argument("--radius", type=float, required=True, help="the matching radius R (m)")
    parser.add_argument("--length-unit", type=float, default=None, help="metres per record length unit (default: Model.L0 of postpro/config_resolved.json, else 1)")
    parser.add_argument("--output", help="write the JSON record here (else stdout)")
    args = parser.parse_args(argv)
    length_unit = args.length_unit
    if length_unit is None:
        resolved = os.path.join(args.run, "postpro", "config_resolved.json")
        length_unit = 1.0
        if os.path.isfile(resolved):
            with open(resolved) as handle:
                length_unit = float(json.load(handle).get("Model", {}).get("L0", 1.0))
    result = estimate(args.run, args.radius, length_unit)
    text = json.dumps(result, indent=2)
    if args.output:
        with open(args.output, "w") as handle:
            handle.write(text + "\n")
    else:
        sys.stdout.write(text + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
