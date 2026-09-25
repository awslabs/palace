#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Scaling series for the geometry identification: one chip-like cell tiled N x N.

Every cell (64 x 64 um) carries one row segment of a CPW (2 um trace, 3 um gaps, 6 um ground
flanks; the trace follows an S-bend of radius 12 um = 6 R, the curved pair class) plus a
rounded island (fillets), a cross (sub-2R concave corners: clusters), a bent bar (90 deg
corners), a straight pair at 1.9 um (translational DifferentConductorGap), a 3 um strip
(SameConductorStrip) and a tight arc bar (CurvedEdge). The CPW of a row is continuous across
the row's N cells and is cut by the truncation box at both ends (port cuts), so the trace and
flank chains grow with N while the number of rows grows with N too: the per-stage wall times
of the identification must grow ~linearly with the number of cells (n log n acceptable); a
chain-length dependence shows up as a super-linear series.

    python3 -m surface_response_identification.scaling_layouts --output DIR --tiles 1 2 4 8 \\
        [--palace BIN] [--ranks 1] [--no-generate] [--rerun]

Writes DIR/meshes/tile-N.msh2 (Gmsh, via synthetic_layouts.jl), DIR/tile-N/np1/{config.json,
palace.log, postpro/}, and DIR/scaling.json with the parsed stage lines
(`Identification <stage>: ... (<seconds> s)`) and the fitted exponent of every stage.
"""

import argparse
import json
import math
import os
import re
import sys
import time

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import synthetic_layouts as S  # noqa: E402
from .preflight_config import preflight_config  # noqa: E402

PITCH = 64.0
DEPTH = 8.0
LC_FINE = 0.75
LC_FAR = 8.0


def translate(lp, dx, dy):
    return S.loop([(x + dx, y + dy) for x, y in lp["Points"]], {k: (cx + dx, cy + dy) for k, (cx, cy) in lp["Arcs"].items()})


def offset_polyline(centreline, distance):
    """Exact offset of a polyline by `distance` (positive = left) with vertex bisectors."""
    result = []
    n = len(centreline)
    for i, p in enumerate(centreline):
        if i == 0:
            d = centreline[1] - p
        elif i == n - 1:
            d = p - centreline[i - 1]
        else:
            d0 = p - centreline[i - 1]
            d1 = centreline[i + 1] - p
            n0 = np.array([-d0[1], d0[0]]) / np.linalg.norm(d0)
            n1 = np.array([-d1[1], d1[0]]) / np.linalg.norm(d1)
            b = n0 + n1
            b /= np.linalg.norm(b)
            result.append(p + (distance / (b @ n0)) * b)
            continue
        normal = np.array([-d[1], d[0]]) / np.linalg.norm(d)
        result.append(p + distance * normal)
    return result


def bar_along(centreline, offset, width):
    right = offset_polyline(centreline, offset - 0.5 * width)
    left = offset_polyline(centreline, offset + 0.5 * width)
    return S.loop([tuple(p) for p in right + left[::-1]])


def cell_centreline(x0, y0, radius=12.0, sweep=30.0, step=5.0):
    """The CPW centreline of one cell: straight, S-bend up (two arcs of +-sweep), straight,
    S-bend down, straight; starts at (x0, y0) heading +x and ends PITCH further at the same y
    heading +x, so that consecutive cells chain into one continuous route."""
    steps = max(1, int(round(sweep / step)))
    turn = math.radians(sweep) / steps
    chord = 2.0 * radius * math.sin(0.5 * turn)  # chord of one polyline step of the arc

    def s_bend(sign):
        """Displacements of an S-bend: an arc turning by sign * sweep then one turning back."""
        moves = []
        heading = 0.0
        for direction in (sign, -sign):
            for _ in range(steps):
                mid = heading + 0.5 * direction * turn
                moves.append(chord * np.array([math.cos(mid), math.sin(mid)]))
                heading += direction * turn
        return moves

    up, down = s_bend(+1.0), s_bend(-1.0)
    straight = (PITCH - sum(m[0] for m in up + down)) / 3.0
    moves = [np.array([straight, 0.0]), *up, np.array([straight, 0.0]), *down, np.array([straight, 0.0])]
    points = [np.array([x0, y0])]
    for move in moves:
        points.append(points[-1] + move)
    # The two S-bends cancel exactly in y and the straights fill the pitch in x: snap the
    # roundoff (1e-14) of the end point so that consecutive cells share one vertex.
    points[-1] = np.array([x0 + PITCH, y0])
    return points


def cell_features(cx, cy):
    """Static features of one cell centred at (cx, cy), outside the CPW band |y - cy| <= 13.3
    (the flanks reach 10 um from the centreline, the S-bends rise 3.2 um); every feature is
    more than 2R = 4 um from the CPW, from its neighbours and from the next cell's features."""
    sheets = []
    # Top band (y in [18, 30]): rounded island, cross, bent bar.
    sheets.append(S.sheet(S.GROUND, translate(S.rounded_rectangle(5.0, 4.0, 1.5), cx - 20.0, cy + 24.0)))
    sheets.append(S.sheet(S.GROUND, translate(S.cross_shape(5.0, 3.0), cx - 2.0, cy + 24.0)))
    sheets.append(S.sheet(S.GROUND, translate(S.bent_bar(4.0, 8.0, 90.0), cx + 17.0, cy + 20.0)))
    # Bottom band (y in [-30, -18]): straight pair at 1.9 um, strip of width 3, tight arc bar.
    sheets.append(S.sheet(S.GROUND, S.rectangle(cx - 26.0, cy - 25.0, cx - 16.0, cy - 21.0)))
    sheets.append(S.sheet(S.GROUND, S.rectangle(cx - 14.1, cy - 25.0, cx - 4.1, cy - 21.0)))
    sheets.append(S.sheet(S.GROUND, S.rectangle(cx + 10.0, cy - 24.5, cx + 24.0, cy - 21.5)))
    sheets.append(S.sheet(S.GROUND, translate(S.arc_bar(3.0, 5.0, 90.0, 10.0, lead=4.0), cx + 0.5, cy - 28.5)))
    return sheets


def tiled_layout(columns, rows=None, name=None, step=5.0):
    """rows x columns cells (rows = columns by default); `step` is the S-bend chord angle."""
    rows = columns if rows is None else rows
    half_x, half_y = 0.5 * columns * PITCH, 0.5 * rows * PITCH
    sheets = []
    for j in range(rows):
        cy = -half_y + (j + 0.5) * PITCH
        centreline = [np.array([-half_x, cy])]
        for i in range(columns):
            centreline.extend(cell_centreline(-half_x + i * PITCH, cy, step=step)[1:])
        # Ground flanks at offsets +-(1 + 3 + 3): 3 um gaps to the 2 um trace.
        sheets.append(S.sheet(S.GROUND, bar_along(centreline, 0.0, 2.0)))
        sheets.append(S.sheet(S.GROUND, bar_along(centreline, 7.0, 6.0)))
        sheets.append(S.sheet(S.GROUND, bar_along(centreline, -7.0, 6.0)))
        for i in range(columns):
            cx = -half_x + (i + 0.5) * PITCH
            sheets.extend(cell_features(cx, cy))
    return S.layout(name or f"tile-{columns}", sheets, half_x=half_x, half_y=half_y, depth=DEPTH, height=DEPTH, lc_fine=LC_FINE, lc_far=LC_FAR, notes=f"{rows} x {columns} chip-like cells of {PITCH:g} um: {rows} CPW rows of {columns} cells (S-bends of radius 12 um, {step:g} deg chords) plus {rows * columns} feature groups")


def row_layout(n, step=2.0):
    """One row of n cells: the three CPW chains grow with n (chain-length scaling) while the
    number of chains stays fixed; finer S-bend chords (2 deg) make the chains longer."""
    return tiled_layout(n, rows=1, name=f"row-{n}", step=step)


STAGE_LINE = re.compile(r"^\s*Identification (?P<stage>[^:]+): (?P<counts>.*) \((?P<seconds>[0-9.]+) s\)\s*$")
TOTAL_LINE = re.compile(r"^\s*Identification total: (?P<seconds>[0-9.]+) s\s*$")


def parse_stages(log_path):
    stages = {}
    counts = {}
    total = None
    elements = None
    peak = None
    with open(log_path) as log:
        for line in log:
            if (m := STAGE_LINE.match(line)):
                stages[m.group("stage")] = float(m.group("seconds"))
                counts[m.group("stage")] = m.group("counts")
            elif (m := TOTAL_LINE.match(line)):
                total = float(m.group("seconds"))
            elif line.startswith(" elements"):
                elements = int(line.split()[-1])
            elif "Estimated peak per-node memory usage" in line and "Total" in line:
                peak = line.split("Total")[-1].strip().split()[0]
    return {"Stages": stages, "Counts": counts, "Total": total, "Elements": elements, "PeakMemory": peak}


def fit_exponent(sizes, seconds, floor=0.02):
    """Least-squares slope of log(seconds) vs log(cells) over the points above `floor` s."""
    pairs = [(n, t) for n, t in zip(sizes, seconds) if t is not None and t > floor]
    if len(pairs) < 2:
        return None
    x = np.log([p[0] for p in pairs])
    y = np.log([p[1] for p in pairs])
    return float(np.polyfit(x, y, 1)[0])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", required=True)
    parser.add_argument("--tiles", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--rows", action="store_true", help="one-row layouts row-N of N cells instead of N x N tiles")
    parser.add_argument("--palace", default=S.DEFAULT_PALACE)
    parser.add_argument("--library", default=S.DEFAULT_LIBRARY)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1])
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--no-generate", action="store_true")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--timeout", type=float, default=3600.0)
    args = parser.parse_args(argv)

    os.makedirs(args.output, exist_ok=True)
    mesh_directory = os.path.join(args.output, "meshes")
    layouts = [row_layout(n) if args.rows else tiled_layout(n) for n in args.tiles]
    missing = [lay for lay in layouts if not os.path.exists(os.path.join(mesh_directory, lay["Name"] + ".msh2"))]
    if missing and not args.no_generate:
        code, seconds, log = S.generate_meshes(missing, mesh_directory, args.julia)
        print(f"generated {len(missing)} meshes: exit {code} in {seconds:.1f} s ({log})", flush=True)
    results = []
    for n, lay in zip(args.tiles, layouts):
        mesh_path = os.path.join(mesh_directory, lay["Name"] + ".msh2")
        if not os.path.exists(mesh_path):
            print(json.dumps({"Tiles": n, "Error": "mesh not generated"}), flush=True)
            continue
        for ranks in args.ranks:
            directory = os.path.join(args.output, lay["Name"], f"np{ranks}")
            os.makedirs(directory, exist_ok=True)
            config = preflight_config(mesh_path, [S.GROUND], [], [S.SUBSTRATE_AIR], args.library, os.path.join(directory, "postpro"))
            config_path = os.path.join(directory, "config.json")
            with open(config_path, "w") as target:
                json.dump(config, target, indent=2)
            log_path = os.path.join(directory, "palace.log")
            manifest_path = os.path.join(directory, "postpro", "surface-response-requirements.json")
            record = {"Tiles": n, "Cells": n if args.rows else n * n, "Ranks": ranks, "Mesh": mesh_path, "Directory": directory}
            if args.rerun or not os.path.exists(manifest_path):
                started = time.time()
                record["ExitCode"] = S.run_preflight(args.palace, config_path, ranks, log_path, args.timeout)
                record["WallSeconds"] = time.time() - started
            record.update(parse_stages(log_path))
            if os.path.exists(manifest_path):
                with open(manifest_path) as source:
                    manifest = json.load(source)
                identification = manifest.get("Identification", {})
                record["GeometryDigest"] = identification.get("GeometryDigest")
                record["Features"] = len(identification.get("Features", []))
                record["Segments"] = len(identification.get("Segments", []))
            results.append(record)
            print(json.dumps({k: record[k] for k in ("Tiles", "Ranks", "Elements", "Segments", "Features", "Total", "Stages") if k in record}), flush=True)
    # Exponents per stage over the np = ranks[0] series.
    series = [r for r in results if r["Ranks"] == args.ranks[0] and r.get("Total") is not None]
    exponents = {}
    if len(series) >= 2:
        sizes = [r["Cells"] for r in series]
        stage_names = sorted({s for r in series for s in r["Stages"]})
        for stage in stage_names:
            exponents[stage] = fit_exponent(sizes, [r["Stages"].get(stage) for r in series])
        exponents["total"] = fit_exponent(sizes, [r["Total"] for r in series])
    with open(os.path.join(args.output, "scaling.json"), "w") as target:
        json.dump({"Results": results, "Exponents": exponents, "Pitch": PITCH, "LcFine": LC_FINE, "LcFar": LC_FAR}, target, indent=1)
    print(json.dumps({"Exponents": exponents}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
