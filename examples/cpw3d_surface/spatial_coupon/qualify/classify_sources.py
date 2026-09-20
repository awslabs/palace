#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Geometric classes of a coupon's sources and per-class statistics of comparisons
(gallery-physics-06b classify_sources.py, case-independent: the classes come from
locate_sources.py rows, the zero-trace knots from the basis contract, the conductor
terminals from the reference config's TerminalAttributes entries - never a hand list).

Classes, in priority order:
  ZeroTrace (apex on the metal cross-section at the cut; contract ZeroTraceIndices) -> raw view only
  junction rings (apex where a metal edge meets the cut, not on the bottom / top face)
  junction columns on the bottom / top faces
  narrow hats w < NARROW_WIDTH next to a junction (junction closer than JUNCTION_REACH)
  near-junction hats (w >= NARROW_WIDTH, junction closer than JUNCTION_REACH)
  narrow hats w < NARROW_WIDTH away from every junction (isolated slivers)
  box 3D corners (bottom / top face corners)
  wide hats, bottom / top faces
  wide hats, metal-top ring
  wide hats, substrate / trench rings
  conductor terminal (a TerminalAttributes source: no trace apex)
Hat width w = distance from the apex to the nearest apex on the same z level.  The
control sources of a qualification are chosen by class (`choose_controls`): a
deterministic geometric choice - one source per class in the priority order above,
cycling until the requested count, the lowest index of each class first.

usage: classify_sources.py --locations source-locations.csv --out-csv PATH --out-md PATH
       [--zero-trace I ...] [--terminal I ...] [label=comparison.json ...]
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics

NARROW_WIDTH = 0.1
JUNCTION_REACH = 0.6
CLASS_ZERO_TRACE = "ZeroTrace (apex on PEC at the cut)"
CLASS_JUNCTION_RINGS = "junction rings (metal edge meets cut)"
CLASS_JUNCTION_COLUMNS = "junction columns on the bottom/top faces"
CLASS_NARROW_JUNCTION = f"narrow hats w < {NARROW_WIDTH:g} um next to a junction"
CLASS_NEAR_JUNCTION = f"near-junction hats ({NARROW_WIDTH:g} <= w, junction < {JUNCTION_REACH:g} um)"
CLASS_NARROW_ISOLATED = f"narrow hats w < {NARROW_WIDTH:g} um away from every junction (isolated slivers)"
CLASS_BOX_CORNERS = "box 3D corners (bottom/top faces)"
CLASS_WIDE_FACES = "wide hats, bottom/top faces"
CLASS_WIDE_METAL_TOP = "wide hats, metal-top ring"
CLASS_WIDE_SUBSTRATE = "wide hats, substrate/trench rings"
CLASS_TERMINAL = "conductor terminal (TerminalAttributes source)"
CLASS_ORDER = [CLASS_WIDE_FACES, CLASS_BOX_CORNERS, CLASS_WIDE_METAL_TOP, CLASS_WIDE_SUBSTRATE, CLASS_NEAR_JUNCTION,
               CLASS_NARROW_JUNCTION, CLASS_NARROW_ISOLATED, CLASS_JUNCTION_COLUMNS, CLASS_JUNCTION_RINGS,
               CLASS_TERMINAL, CLASS_ZERO_TRACE]
WIDE_CLASSES = (CLASS_WIDE_FACES, CLASS_BOX_CORNERS, CLASS_WIDE_METAL_TOP, CLASS_WIDE_SUBSTRATE)
OFFSET_KEYS = ("E_rel", "p_MA_rel", "p_MS_rel", "p_SA_rel")


def read_locations(path):
    with open(path, newline="") as stream:
        return {int(row["index"]): row for row in csv.DictReader(stream)}


def hat_width(locations, i):
    ri = locations[i]
    xi, yi, zi = float(ri["x"]), float(ri["y"]), float(ri["z"])
    best = math.inf
    for j, rj in locations.items():
        if j == i or abs(float(rj["z"]) - zi) > 1e-9:
            continue
        best = min(best, math.hypot(float(rj["x"]) - xi, float(rj["y"]) - yi))
    return best


def classify_all(locations, zero_trace, terminals=()):
    """Index -> (class, hat width) for every located source and every terminal."""
    zero_trace = set(zero_trace)
    classes = {}
    for i in sorted(set(locations) | set(terminals)):
        if i in terminals and i not in locations:
            classes[i] = (CLASS_TERMINAL, math.nan)
            continue
        row = locations[i]
        width = hat_width(locations, i)
        on_face = row["z_role"] in ("bottom", "top")
        distance = float(row["junction_distance_um"]) if row["junction_distance_um"] else math.inf
        if i in zero_trace:
            name = CLASS_ZERO_TRACE
        elif row["metal_edge_junction_here"] == "1" and not on_face:
            name = CLASS_JUNCTION_RINGS
        elif row["metal_edge_junction_here"] == "1":
            name = CLASS_JUNCTION_COLUMNS
        elif distance < JUNCTION_REACH and width < NARROW_WIDTH:
            name = CLASS_NARROW_JUNCTION
        elif distance < JUNCTION_REACH:
            name = CLASS_NEAR_JUNCTION
        elif width < NARROW_WIDTH:
            name = CLASS_NARROW_ISOLATED
        elif row["box_corner_3d"] == "1":
            name = CLASS_BOX_CORNERS
        elif on_face:
            name = CLASS_WIDE_FACES
        elif row["z_role"] == "metal-top":
            name = CLASS_WIDE_METAL_TOP
        else:
            name = CLASS_WIDE_SUBSTRATE
        classes[i] = (name, width)
    return classes


def choose_controls(classes, count, free):
    """`count` control sources by class: cycle over the classes in CLASS_ORDER (free
    sources only, the zero-trace class excluded), taking the lowest-index member of each
    class not yet chosen, until `count` sources are chosen (every class represented
    before any class contributes a second source); sorted by index."""
    members = {name: [i for i in sorted(free) if classes.get(i, (None,))[0] == name]
               for name in CLASS_ORDER if name != CLASS_ZERO_TRACE}
    chosen = []
    while len(chosen) < count and any(members.values()):
        for name in CLASS_ORDER:
            if name in members and members[name] and len(chosen) < count:
                chosen.append(members[name].pop(0))
    if len(chosen) < count:
        raise ValueError(f"only {len(chosen)} free sources for {count} controls")
    return sorted(chosen)


def class_statistics(classes, comparison):
    """Per class: n, members and the |offset| median / max, signed median, worst source
    and within-1/2/5% counts of every offset key of a compare_matrices summary."""
    per_source = comparison["PerSource"]
    summary = {}
    for name in CLASS_ORDER:
        members = [i for i, (class_name, _) in classes.items() if class_name == name and str(i) in per_source]
        if not members:
            continue
        record = {"n": len(members), "members": members}
        for key in OFFSET_KEYS:
            values = {i: per_source[str(i)][key] for i in members
                      if key in per_source[str(i)] and math.isfinite(per_source[str(i)][key])}
            if not values:
                continue
            worst = max(values, key=lambda i: abs(values[i]))
            record[key] = {"n": len(values), "abs_median": statistics.median(abs(v) for v in values.values()),
                           "abs_max": abs(values[worst]), "signed_median": statistics.median(values.values()),
                           "worst_source": worst, "worst_offset": values[worst],
                           "within_1pct": sum(abs(v) < 0.01 for v in values.values()),
                           "within_2pct": sum(abs(v) < 0.02 for v in values.values()),
                           "within_5pct": sum(abs(v) < 0.05 for v in values.values())}
        summary[name] = record
    return summary


def fmt(x):
    return "n/a" if x is None or not math.isfinite(x) else f"{x:.2e}"


def write_class_report(classes, locations, comparisons, out_csv, out_md):
    """source-classes.csv (index, class, width, apex, offsets per comparison) and the
    per-class Markdown / JSON statistics; returns {label: class statistics}."""
    out_csv, out_md = Path(out_csv), Path(out_md)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["index", "class", "hat_width_um", "x", "y", "z"]
                        + [f"{label}:{key}" for label, _ in comparisons for key in OFFSET_KEYS])
        for i, (name, width) in classes.items():
            location = locations.get(i, {"x": "", "y": "", "z": ""})
            row = [i, name, f"{width:.4f}" if math.isfinite(width) else "", location["x"], location["y"], location["z"]]
            for _, comparison in comparisons:
                offsets = comparison["PerSource"].get(str(i), {})
                row += [offsets.get(key, "") for key in OFFSET_KEYS]
            writer.writerow(row)
    lines = ["# Per-class statistics (|rel diff| median / max; signed median in brackets; n within 1% / 2% / 5% for p_MA)", ""]
    summary = {}
    for label, comparison in comparisons:
        statistics_by_class = class_statistics(classes, comparison)
        summary[label] = statistics_by_class
        lines += [f"## {label}", "",
                  "| class | n | E_ii median / max [signed median] | p_SA median / max [signed median] | p_MS median / max [signed median] "
                  "| p_MA median / max [signed median] | p_MA <1% / <2% / <5% | sources (p_MA rel) |",
                  "|---|---:|---|---|---|---|---|---|"]
        for name, record in statistics_by_class.items():
            cells = [name, str(record["n"])]
            for key in ("E_rel", "p_SA_rel", "p_MS_rel", "p_MA_rel"):
                stats = record.get(key)
                cells.append("n/a" if stats is None else
                             f"{fmt(stats['abs_median'])} / {fmt(stats['abs_max'])} [{stats['signed_median']:+.2e}]")
            ma = record.get("p_MA_rel")
            cells.append(f"{ma['within_1pct']} / {ma['within_2pct']} / {ma['within_5pct']}" if ma else "n/a")
            ranked = sorted(record["members"], key=lambda i: -abs(comparison["PerSource"][str(i)].get("p_MA_rel", 0.0)))
            cells.append(", ".join(f"{i} ({comparison['PerSource'][str(i)].get('p_MA_rel', math.nan):+.1e})" for i in ranked[:6]))
            lines.append("| " + " | ".join(cells) + " |")
        lines.append("")
    out_md.write_text("\n".join(lines) + "\n")
    out_md.with_suffix(".json").write_text(json.dumps(summary, indent=1) + "\n")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--locations", required=True)
    parser.add_argument("--out-csv", required=True)
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--zero-trace", type=int, action="append", default=[])
    parser.add_argument("--terminal", type=int, action="append", default=[])
    parser.add_argument("comparisons", nargs="*", metavar="label=comparison.json")
    args = parser.parse_args(argv)
    locations = read_locations(args.locations)
    classes = classify_all(locations, args.zero_trace, args.terminal)
    comparisons = [(item.split("=", 1)[0], json.loads(Path(item.split("=", 1)[1]).read_text())) for item in args.comparisons]
    write_class_report(classes, locations, comparisons, args.out_csv, args.out_md)
    print(Path(args.out_md).read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
