#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Compare two validation manifests of one geometry (validate_window_mesh.jl output), e.g. the
recorded transmon reference mesh and its regeneration from the exported polygon set:

    python3 compare_mesh_validations.py RECORDED.validation.json REGENERATED.validation.json \\
        [--merge 6+9=6] [--area-tolerance 1e-6] [--output compare.json] [--markdown compare.md]

Reports node / tetrahedron / surface-triangle counts, per-attribute element counts, areas and
volumes with their relative differences, adjacency errors and quality. ``--merge A+B=C`` sums
attributes A and B of the SECOND manifest into attribute C before comparing (the recorded
transmon mesh carries the substrate backside inside attribute 6; the mesher writes it as 9).
Exit status 0 when every compared area and volume agrees within ``--area-tolerance``
(relative) and both meshes have zero adjacency errors and zero nonpositive tetrahedra.
"""

import argparse
import json
import sys


def merged(values, merges):
    out = dict(values)
    for sources, target in merges:
        total = 0.0
        for source in sources:
            total += out.pop(source, 0.0)
        out[target] = out.get(target, 0.0) + total
    return out


def relative(a, b):
    scale = max(abs(a), abs(b), 1e-300)
    return (b - a) / scale


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("recorded")
    parser.add_argument("regenerated")
    parser.add_argument("--merge", action="append", default=[], help="A+B=C on the second manifest")
    parser.add_argument("--area-tolerance", type=float, default=1e-6)
    parser.add_argument("--output")
    parser.add_argument("--markdown")
    args = parser.parse_args()
    a = json.load(open(args.recorded))
    b = json.load(open(args.regenerated))
    merges = []
    for rule in args.merge:
        sources, target = rule.split("=")
        merges.append((sources.split("+"), target))

    rows = []
    worst = 0.0
    for key, label in (("surface_area_um2", "area"), ("volume_um3", "volume"), ("surface_attribute_counts", "faces"),
                       ("volume_attribute_counts", "tets")):
        va = a.get(key, {})
        vb = merged(b.get(key, {}), merges) if key.startswith("surface") else b.get(key, {})
        for attribute in sorted(set(va) | set(vb), key=lambda s: int(s)):
            x, y = va.get(attribute), vb.get(attribute)
            rel = relative(x, y) if x is not None and y is not None else None
            rows.append({"Quantity": label, "Attribute": attribute, "Recorded": x, "Regenerated": y, "Relative": rel})
            if label in ("area", "volume") and rel is not None:
                worst = max(worst, abs(rel))
    totals = [{"Quantity": name, "Recorded": a.get(name), "Regenerated": b.get(name),
               "Relative": relative(a.get(name, 0), b.get(name, 0))}
              for name in ("nodes", "tetrahedra", "surface_triangles")]
    adjacency_ok = all(v == 0 for v in a["surface_adjacency_errors"].values()) and \
        all(v == 0 for v in b["surface_adjacency_errors"].values())
    positive_ok = a["mesh_quality"]["nonpositive"] == 0 and b["mesh_quality"]["nonpositive"] == 0
    verdict = worst <= args.area_tolerance and adjacency_ok and positive_ok
    result = {"Recorded": a["mesh"], "Regenerated": b["mesh"], "Merges": args.merge, "Totals": totals, "Rows": rows,
              "WorstRelativeAreaOrVolume": worst, "AreaTolerance": args.area_tolerance,
              "AdjacencyErrorsZero": adjacency_ok, "NonpositiveZero": positive_ok,
              "Quality": {"Recorded": a["mesh_quality"], "Regenerated": b["mesh_quality"]}, "Pass": verdict}
    if args.output:
        json.dump(result, open(args.output, "w"), indent=2)
    lines = ["| quantity | attribute | recorded | regenerated | relative |", "|---|---|---:|---:|---:|"]
    for t in totals:
        lines.append(f"| {t['Quantity']} | | {t['Recorded']} | {t['Regenerated']} | {t['Relative']:+.3e} |")
    for r in rows:
        rel = "" if r["Relative"] is None else f"{r['Relative']:+.3e}"
        lines.append(f"| {r['Quantity']} | {r['Attribute']} | {r['Recorded']} | {r['Regenerated']} | {rel} |")
    lines.append("")
    lines.append(f"Worst relative area / volume difference: {worst:.3e} (tolerance {args.area_tolerance:g}); "
                 f"adjacency errors zero: {adjacency_ok}; nonpositive zero: {positive_ok}; PASS: {verdict}")
    text = "\n".join(lines)
    print(text)
    if args.markdown:
        open(args.markdown, "w").write(text + "\n")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
