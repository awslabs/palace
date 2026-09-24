#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Synthetic stress layouts with exact oracles for the geometry identification.

A layout is a set of zero-thickness metal sheets (plan-view polygons with optional circular
arcs and holes) on the process plane of a substrate / vacuum box, optionally a facing sheet on
a second plane and vertical walls. This module (1) defines the stress suite of the
identification plan (parallel edges through R and 2R, corners, junction-like polygons, an
edge ending at the truncation boundary, sub-2R islands, apertures, polyline arcs, a taper, a
facing layer and a wall, the rounded 8 x 6 um island and its aperture), (2) writes the
specification synthetic_layouts.jl meshes with Gmsh, (3) computes the layout oracle from the
polygons (perimeter length, corner list with angles and convexity, parallel edge pairs with
separations and the expected translational class, non-parallel interactions, arcs), and (4)
runs the preflight matrix (ranks x uniform levels) with the audit on every mesh and compares
the audit perimeter and the manifest with the oracle.

    python3 -m surface_response_identification.synthetic_layouts --output DIR \\
        [--names a b ...] [--ranks 1 2] [--uniform-levels 0 1] [--library seed.json]

The classifier is never modified here: every disagreement is recorded as the observed rule.
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time
from collections import Counter

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import audit, manifest as M, perimeter as P  # noqa: E402
from .msh2 import read_msh2  # noqa: E402
from .preflight_config import preflight_config  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
DEFAULT_PALACE = os.path.join(REPO, "build", "bin", "palace")
DEFAULT_LIBRARY = os.path.join(REPO, "examples", "transmon", "benchmark", "transmon_surface_process_seed.json")
# A committed Version-3 library with an isolated-edge model (and a 2 um same-conductor
# strip): with it the classifier records every isolated segment instead of one representative.
ISOLATED_LIBRARY = os.path.join(REPO, "examples", "transmon", "benchmark", "corrected-p1", "library", "process-library.json")
GENERATOR = os.path.join(HERE, "synthetic_layouts.jl")

RADIUS = 2.0
GROUND = 4
SECOND_CONDUCTOR = 5
SUBSTRATE_AIR = 8
CORNER_TURN_TOLERANCE_DEGREES = P.CORNER_ANGLE_TOLERANCE_DEGREES


# ----------------------------------------------------------------------------- geometry ----


def loop(points, arcs=None):
    """A closed loop: points counter-clockwise for outer boundaries; arcs maps the index of an
    end vertex to the arc centre (the arc runs from the previous vertex to that vertex)."""
    return {"Points": [(float(x), float(y)) for x, y in points], "Arcs": {int(k): (float(v[0]), float(v[1])) for k, v in (arcs or {}).items()}}


def sheet(attribute, outer, holes=(), z=0.0):
    return {"Attribute": attribute, "Z": z, "Loops": [outer, *holes]}


def layout(name, sheets, walls=(), half_x=30.0, half_y=30.0, depth=15.0, height=15.0, lc_fine=1.0, lc_far=6.0, notes=""):
    return {"Name": name, "HalfX": half_x, "HalfY": half_y, "Depth": depth, "Height": height, "LcFine": lc_fine, "LcFar": lc_far, "Sheets": list(sheets), "Walls": list(walls), "Notes": notes}


def rectangle(x0, y0, x1, y1):
    return loop([(x0, y0), (x1, y0), (x1, y1), (x0, y1)])


def rounded_rectangle(half_x, half_y, radius):
    """Counter-clockwise rounded rectangle: 4 lines joined by 4 quarter arcs."""
    points = [
        (half_x - radius, -half_y), (half_x, -half_y + radius),  # arc 1 -> 2
        (half_x, half_y - radius), (half_x - radius, half_y),  # arc 3 -> 4
        (-half_x + radius, half_y), (-half_x, half_y - radius),  # arc 5 -> 6
        (-half_x, -half_y + radius), (-half_x + radius, -half_y),  # arc 7 -> 8
    ]
    centres = {2: (half_x - radius, -half_y + radius), 4: (half_x - radius, half_y - radius), 6: (-half_x + radius, half_y - radius), 8: (-half_x + radius, -half_y + radius)}
    return loop(points, centres)


def bent_bar(width, arm, interior_angle_degrees):
    """A bar of the given width with two arms of the given length meeting at a bend whose
    outer (convex) corner has the given interior angle; the inner corner is concave with
    the same interior angle on the metal side measured as 360 - angle."""
    theta = math.radians(interior_angle_degrees)
    # Centreline: from A along -x to the bend at the origin, then along direction rotated by (180 - theta).
    turn = math.pi - theta
    d1 = np.array([1.0, 0.0])
    d2 = np.array([math.cos(turn), math.sin(turn)])
    n1 = np.array([0.0, 1.0])  # left normal of d1
    n2 = np.array([-d2[1], d2[0]])
    h = 0.5 * width
    start = -arm * d1
    end = arm * d2
    # Offset intersections at the bend (left = inner side of a left turn).
    bisector = n1 + n2
    scale = h / (n1 @ bisector / np.linalg.norm(bisector)) / np.linalg.norm(bisector)
    inner = bisector * scale
    outer = -inner
    points = [start - h * n1, outer, end - h * n2, end + h * n2, inner, start + h * n1]
    return loop([tuple(p) for p in points])


def arc_bar(width, radius, sweep_degrees, step_degrees, lead=6.0):
    """A bar following a circular arc discretised as a polyline with the given turning angle
    per vertex, with straight leads at both ends."""
    steps = max(1, int(round(sweep_degrees / step_degrees)))
    step = math.radians(sweep_degrees) / steps
    h = 0.5 * width
    centre = np.array([0.0, radius])
    centreline = [centre + radius * np.array([math.sin(k * step), -math.cos(k * step)]) for k in range(steps + 1)]
    d_start = np.array([1.0, 0.0])
    d_end = np.array([math.cos(steps * step), math.sin(steps * step)])
    centreline = [centreline[0] - lead * d_start, *centreline, centreline[-1] + lead * d_end]
    # Offset polyline by +-h using the vertex bisector normals (exact offset of the polyline).
    def offset(sign):
        result = []
        n = len(centreline)
        for i, p in enumerate(centreline):
            if i == 0:
                d = centreline[1] - p
                normal = np.array([-d[1], d[0]]) / np.linalg.norm(d)
                result.append(p + sign * h * normal)
            elif i == n - 1:
                d = p - centreline[i - 1]
                normal = np.array([-d[1], d[0]]) / np.linalg.norm(d)
                result.append(p + sign * h * normal)
            else:
                d0 = p - centreline[i - 1]
                d1 = centreline[i + 1] - p
                n0 = np.array([-d0[1], d0[0]]) / np.linalg.norm(d0)
                n1 = np.array([-d1[1], d1[0]]) / np.linalg.norm(d1)
                b = n0 + n1
                b /= np.linalg.norm(b)
                result.append(p + sign * (h / (b @ n0)) * b)
        return result

    right = offset(-1.0)
    left = offset(1.0)
    return loop([tuple(p) for p in right + left[::-1]])


def trapezoid(bottom_width, top_width, height):
    return loop([(-0.5 * bottom_width, 0.0), (0.5 * bottom_width, 0.0), (0.5 * top_width, height), (-0.5 * top_width, height)])


def t_shape(bar_length, bar_width, stem_length, stem_width):
    """Bar along x at the top, stem down from its middle; counter-clockwise."""
    hb, hs = 0.5 * bar_length, 0.5 * stem_width
    return loop([(-hs, -stem_length), (hs, -stem_length), (hs, 0.0), (hb, 0.0), (hb, bar_width), (-hb, bar_width), (-hb, 0.0), (-hs, 0.0)])


def cross_shape(arm_length, arm_width):
    a, h = arm_length, 0.5 * arm_width
    return loop([(-h, -a), (h, -a), (h, -h), (a, -h), (a, h), (h, h), (h, a), (-h, a), (-h, h), (-a, h), (-a, -h), (-h, -h)])


# ------------------------------------------------------------------------------- suite ----


def stress_suite():
    layouts = []
    separations = [1.8, 1.9, 1.95, 2.0, 2.05, 2.1, 3.9, 4.0, 4.1]
    for s in separations:
        tag = f"{s:g}".replace(".", "p")
        # Gap between two ground rectangles (same conductor).
        layouts.append(layout(f"gap-same-{tag}", [sheet(GROUND, rectangle(-14.0 - s / 2, -6.0, -s / 2, 6.0)), sheet(GROUND, rectangle(s / 2, -6.0, 14.0 + s / 2, 6.0))], notes=f"SameConductorGap at separation {s}"))
        # Gap between ground and a second conductor.
        layouts.append(layout(f"gap-different-{tag}", [sheet(GROUND, rectangle(-14.0 - s / 2, -6.0, -s / 2, 6.0)), sheet(SECOND_CONDUCTOR, rectangle(s / 2, -6.0, 14.0 + s / 2, 6.0))], notes=f"DifferentConductorGap at separation {s}"))
        # Strip of width s.
        layouts.append(layout(f"strip-{tag}", [sheet(GROUND, rectangle(-12.0, -s / 2, 12.0, s / 2))], lc_fine=min(1.0, s / 2), notes=f"SameConductorStrip at separation {s}"))
    for angle in [30, 45, 60, 90, 120, 135, 150, 170, 175, 178]:
        layouts.append(layout(f"corner-{angle}", [sheet(GROUND, bent_bar(6.0, 16.0, angle))], notes=f"bent bar: one convex and one concave corner of interior angle {angle} deg, four 90 deg arm ends"))
    layouts.append(layout("tee-stem3", [sheet(GROUND, t_shape(24.0, 6.0, 14.0, 3.0))], notes="T: two concave 90 deg corners 3 um apart (< 2R), stem edges at separation 3"))
    layouts.append(layout("tee-stem6", [sheet(GROUND, t_shape(24.0, 6.0, 14.0, 6.0))], notes="T: two concave 90 deg corners 6 um apart (> 2R)"))
    layouts.append(layout("cross-arm3", [sheet(GROUND, cross_shape(12.0, 3.0))], notes="X: four concave corners at 3 and 4.24 um, arm edges at separation 3"))
    layouts.append(layout("edge-to-boundary", [sheet(GROUND, rectangle(-30.0, -3.0, 6.0, 3.0))], notes="strip of width 6 entering from the x = -30 truncation boundary: two chain endpoints at the truncation"))
    layouts.append(layout("island-3x3", [sheet(GROUND, rectangle(-1.5, -1.5, 1.5, 1.5))], lc_fine=0.5, notes="isolated island smaller than 2R: every edge within 2R of every other"))
    for hole in [1.0, 3.0, 10.0]:
        tag = f"{hole:g}"
        layouts.append(layout(f"hole-{tag}", [sheet(GROUND, rectangle(-30.0, -30.0, 30.0, 30.0), holes=[rectangle(-hole / 2, -hole / 2, hole / 2, hole / 2)])], lc_fine=min(1.0, hole / 3), notes=f"square aperture of side {hole} in a ground plane reaching the truncation boundary"))
    for radius, sweep in [(5.0, 90.0), (20.0, 90.0), (50.0, 45.0), (250.0, 15.0)]:
        for step in [1.0, 5.0, 20.0]:
            if step > sweep:
                continue
            extent = radius * math.sin(math.radians(sweep)) + 12.0
            half = max(30.0, math.ceil(extent + 4.0))
            lc = 1.0 if radius < 100 else 2.0
            layouts.append(
                layout(f"arc-r{radius:g}-step{step:g}", [sheet(GROUND, arc_bar(3.0, radius, sweep, step))], half_x=half, half_y=half, lc_fine=lc, lc_far=max(6.0, half / 5.0), notes=f"bar of width 3 on a polyline arc of radius {radius}, sweep {sweep} deg, {step} deg per vertex")
            )
    layouts.append(layout("arc-r20-step5-fine", [sheet(GROUND, arc_bar(3.0, 20.0, 90.0, 5.0))], half_x=36.0, half_y=36.0, lc_fine=0.5, lc_far=6.0, notes="arc-r20-step5 at half the mesh size (mesh independence)"))
    layouts.append(layout("taper-10", [sheet(GROUND, trapezoid(24.0, 24.0 - 2.0 * 20.0 * math.tan(math.radians(10.0)), 20.0))], notes="trapezoid with two 10 deg taper edges: corners 80 and 100 deg"))
    layouts.append(layout("facing-layers", [sheet(GROUND, rectangle(-10.0, -6.0, 10.0, 6.0)), sheet(GROUND, rectangle(-10.0, -6.0, 10.0, 6.0), z=3.0)], notes="two facing metal sheets 3 um apart in one PEC attribute: cross-layer class, excluded by decision 73(3)"))
    layouts.append(layout("vertical-wall", [sheet(GROUND, rectangle(-10.0, -6.0, 10.0, 6.0))], walls=[(GROUND, 0.0, -6.0, 0.0, 6.0, 4.0)], notes="a vertical metal wall standing across the sheet: non-planar and non-manifold classes"))
    layouts.append(layout("island-rounded-8x6", [sheet(GROUND, rounded_rectangle(4.0, 3.0, 0.5))], half_x=12.0, half_y=12.0, depth=8.0, height=8.0, lc_fine=0.25, lc_far=2.0, notes="survey geometry 2: 8 x 6 um island with 0.5 um fillets (perimeter 24 + pi)"))
    layouts.append(layout("aperture-rounded-8x6", [sheet(GROUND, rectangle(-12.0, -12.0, 12.0, 12.0), holes=[rounded_rectangle(4.0, 3.0, 0.5)])], half_x=12.0, half_y=12.0, depth=8.0, height=8.0, lc_fine=0.25, lc_far=2.0, notes="survey geometry 2, aperture variant: the same shape as a hole in a ground plane"))
    return layouts


# ----------------------------------------------------------------------- specification ----


def write_specification(layouts, path):
    lines = []
    for lay in layouts:
        lines.append(f"layout {lay['Name']}")
        lines.append(f"box {lay['HalfX']!r} {lay['HalfY']!r} {lay['Depth']!r} {lay['Height']!r}")
        lines.append(f"size {lay['LcFine']!r} {lay['LcFar']!r}")
        for sh in lay["Sheets"]:
            lines.append(f"polygon {sh['Attribute']}" if sh["Z"] == 0.0 else f"sheet {sh['Attribute']} {sh['Z']!r}")
            for index, lp in enumerate(sh["Loops"]):
                if index > 0:
                    lines.append("hole")
                for k, (x, y) in enumerate(lp["Points"], start=1):
                    if k in lp["Arcs"]:
                        cx, cy = lp["Arcs"][k]
                        lines.append(f"arc {x!r} {y!r} {cx!r} {cy!r}")
                    else:
                        lines.append(f"v {x!r} {y!r}")
                lines.append("close")
        for attribute, x0, y0, x1, y1, h in lay["Walls"]:
            lines.append(f"wall {attribute} {x0!r} {y0!r} {x1!r} {y1!r} {h!r}")
        lines.append("end")
    with open(path, "w") as target:
        target.write("\n".join(lines) + "\n")
    return path


# -------------------------------------------------------------------------------- oracle ----


def cross2(a, b):
    return float(a[0] * b[1] - a[1] * b[0])


def point_segment_distance(point, start, end):
    d = end - start
    t = float(np.clip(((point - start) @ d) / (d @ d), 0.0, 1.0))
    return float(np.linalg.norm(point - (start + t * d)))


def within_interaction(distance, radius):
    """The classifier's quantized decision: within 2R unless within half a length quantum
    of 2R (kDecisionLengthQuantumRelativeToMatchingRadius = 1e-8)."""
    return distance < 2.0 * radius - 0.5 * 1.0e-8 * radius


def design_event_cores(all_edges, corner_points, radius, samples=400):
    """Event cores of SURFACE-RESPONSE-IDENTIFICATION.md (b) 3 from the polygon edges: sampled
    points on an edge that are within 2R of a point on a non-parallel edge of another chain,
    excluding through-vertex pairs (both points within 2R of a vertex the two chains share).
    Edges joined by a sub-threshold vertex (not a classifier corner) are one chain."""
    n = len(all_edges)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def is_corner(point):
        return any(np.linalg.norm(point - q) < 1.0e-9 for q in corner_points)

    shared_vertices = {}
    for i in range(n):
        for j in range(i + 1, n):
            for p in (all_edges[i]["Start"], all_edges[i]["End"]):
                for q in (all_edges[j]["Start"], all_edges[j]["End"]):
                    if np.allclose(p, q):
                        shared_vertices.setdefault((i, j), []).append(p)
                        if not is_corner(p):
                            parent[find(i)] = find(j)
    interaction = 2.0 * radius - 0.5 * 1.0e-8 * radius
    cores = []
    ts = np.linspace(0.0, 1.0, samples + 1)
    for i in range(n):
        for j in range(n):
            if i == j or find(i) == find(j):
                continue
            a, b = all_edges[i], all_edges[j]
            if abs(float(a["Tangent"] @ b["Tangent"])) >= 1.0 - 1.0e-8:
                continue
            if P.segment_distance(np.append(a["Start"], 0.0), np.append(a["End"], 0.0), np.append(b["Start"], 0.0), np.append(b["End"], 0.0)) >= interaction:
                continue
            pa = a["Start"][None, :] + ts[:, None] * (a["End"] - a["Start"])[None, :]
            pb = b["Start"][None, :] + ts[:, None] * (b["End"] - b["Start"])[None, :]
            close = np.linalg.norm(pa[:, None, :] - pb[None, :, :], axis=2) < interaction
            for v in shared_vertices.get((min(i, j), max(i, j)), []):
                za = np.linalg.norm(pa - v, axis=1) < interaction
                zb = np.linalg.norm(pb - v, axis=1) < interaction
                close &= ~(za[:, None] & zb[None, :])
            hits = pa[close.any(axis=1)]
            if len(hits):
                cores.append(hits)
    # Two corners closer than 2R (overlapping windows) are an event of their own.
    for i, p in enumerate(corner_points):
        for q in corner_points[i + 1 :]:
            if np.linalg.norm(p - q) < interaction:
                cores.append(np.array([p, q]))
    return np.concatenate(cores) if cores else np.zeros((0, 2))


def _on_box(point, lay, tolerance=1.0e-9):
    x, y = point
    return abs(abs(x) - lay["HalfX"]) < tolerance or abs(abs(y) - lay["HalfY"]) < tolerance


def loop_edges(lp, hole):
    """Straight edges and arcs of a loop with the inward (metal-side) normal. Outer loops are
    counter-clockwise so the metal is on the left of every edge; holes are the reverse."""
    points = lp["Points"]
    n = len(points)
    edges = []
    arcs = []
    for i in range(n):
        a = np.array(points[i])
        b = np.array(points[(i + 1) % n])
        end_index = (i + 1) % n + 1
        if end_index in lp["Arcs"]:
            centre = np.array(lp["Arcs"][end_index])
            r = float(np.linalg.norm(a - centre))
            a0 = math.atan2(*(a - centre)[::-1])
            a1 = math.atan2(*(b - centre)[::-1])
            sweep = a1 - a0
            if sweep > math.pi:
                sweep -= 2.0 * math.pi
            if sweep < -math.pi:
                sweep += 2.0 * math.pi
            arcs.append({"Start": a, "End": b, "Centre": centre, "Radius": r, "SweepDegrees": math.degrees(sweep), "Length": r * abs(sweep)})
        else:
            d = b - a
            length = float(np.linalg.norm(d))
            t = d / length
            left = np.array([-t[1], t[0]])
            edges.append({"Start": a, "End": b, "Tangent": t, "Inward": -left if hole else left, "Length": length})
    return edges, arcs


def _tangent_at_vertex(lp, index, incoming):
    """Unit tangent of the loop at vertex index (0-based) on the incoming or outgoing side."""
    points = lp["Points"]
    n = len(points)
    p = np.array(points[index])
    if incoming:
        arc_key = index + 1
        if arc_key in lp["Arcs"]:
            centre = np.array(lp["Arcs"][arc_key])
            radial = p - centre
            previous = np.array(points[(index - 1) % n])
            # Tangent direction consistent with the traversal sense of the arc.
            sense = np.sign(cross2(previous - centre, p - centre))
            t = sense * np.array([-radial[1], radial[0]])
        else:
            t = p - np.array(points[(index - 1) % n])
    else:
        next_index = (index + 1) % n
        arc_key = next_index + 1
        if arc_key in lp["Arcs"]:
            centre = np.array(lp["Arcs"][arc_key])
            radial = p - centre
            following = np.array(points[next_index])
            sense = np.sign(cross2(p - centre, following - centre))
            t = sense * np.array([-radial[1], radial[0]])
        else:
            t = np.array(points[next_index]) - p
    return t / np.linalg.norm(t)


def oracle(lay, radius=RADIUS, corner_turn_tolerance=CORNER_TURN_TOLERANCE_DEGREES):
    """Exact layout features of the process-plane sheets: perimeter length (physical, i.e.
    not on the truncation box), corner list, parallel edge pairs within 2R with separations
    and the expected translational class, non-parallel interactions, arcs, and the excluded
    off-plane metal."""
    all_edges = []
    corners = []
    arcs = []
    physical_length = 0.0
    truncation_length = 0.0
    excluded = {"CrossLayerSheets": 0, "Walls": len(lay["Walls"])}
    for sh in lay["Sheets"]:
        if sh["Z"] != 0.0:
            excluded["CrossLayerSheets"] += 1
            continue
        for hole_index, lp in enumerate(sh["Loops"]):
            hole = hole_index > 0
            edges, loop_arcs = loop_edges(lp, hole)
            for e in edges:
                on_box = _on_box(e["Start"], lay) and _on_box(e["End"], lay) and (abs(abs(e["Start"][0]) - lay["HalfX"]) < 1e-9 and abs(abs(e["End"][0]) - lay["HalfX"]) < 1e-9 or abs(abs(e["Start"][1]) - lay["HalfY"]) < 1e-9 and abs(abs(e["End"][1]) - lay["HalfY"]) < 1e-9)
                e["Truncation"] = on_box
                e["Attribute"] = sh["Attribute"]
                e["Conductor"] = 0 if sh["Attribute"] == GROUND else 1
                if on_box:
                    truncation_length += e["Length"]
                else:
                    physical_length += e["Length"]
                    all_edges.append(e)
            for a in loop_arcs:
                physical_length += a["Length"]
                sweep = abs(a["SweepDegrees"])
                convex = (a["SweepDegrees"] > 0) != hole
                tangent_distance = a["Radius"] * math.tan(math.radians(0.5 * sweep))
                rounded = a["Radius"] < radius and tangent_distance < radius
                arcs.append(
                    {
                        "Radius": a["Radius"],
                        "SweepDegrees": round(a["SweepDegrees"], 9),
                        "Length": a["Length"],
                        "BendRadiusOverR": a["Radius"] / radius,
                        "TangentDistance": tangent_distance,
                        "Convex": bool(convex),
                        "Expected": (("ConvexCorner" if convex else "ConcaveCorner") + f"@{round(180.0 - sweep, 6)} radius {a['Radius']}") if rounded else "straight chain (bend radius or tangent distance >= R)",
                        "ExpectedRoundedCorner": rounded,
                    }
                )
            points = lp["Points"]
            for index, point in enumerate(points):
                if _on_box(point, lay):
                    continue
                t_in = _tangent_at_vertex(lp, index, True)
                t_out = _tangent_at_vertex(lp, index, False)
                cosine = float(np.clip(t_in @ t_out, -1.0, 1.0))
                turn = math.degrees(math.acos(cosine))
                if turn < 1.0e-9:
                    continue
                cross = cross2(t_in, t_out)
                # Metal on the left for outer loops: a left turn (cross > 0) is convex.
                convex = (cross > 0) != hole
                interior = 180.0 - turn if convex else 180.0 + turn
                corners.append(
                    {
                        "Point": [round(float(point[0]), 9), round(float(point[1]), 9)],
                        "TurnDegrees": round(turn, 9),
                        "InteriorAngleDegrees": round(interior, 9),
                        "Convex": bool(convex),
                        "ClassifierCorner": turn > corner_turn_tolerance + 1.0e-9,
                        "Expected": ("ConvexCorner" if convex else "ConcaveCorner") if turn > corner_turn_tolerance + 1.0e-9 else "straight (turn <= 30 deg tolerance)",
                    }
                )
    # Parallel and non-parallel pairs within 2R (non-adjacent straight physical edges).
    pairs = []
    nonparallel = []
    for i in range(len(all_edges)):
        for j in range(i + 1, len(all_edges)):
            a, b = all_edges[i], all_edges[j]
            shared = any(np.allclose(p, q) for p in (a["Start"], a["End"]) for q in (b["Start"], b["End"]))
            distance = P.segment_distance(np.append(a["Start"], 0.0), np.append(a["End"], 0.0), np.append(b["Start"], 0.0), np.append(b["End"], 0.0))
            if distance > 2.0 * radius * (1.0 + 1.0e-9):
                continue
            cosine = abs(float(a["Tangent"] @ b["Tangent"]))
            if cosine < 1.0 - 1.0e-8 and not within_interaction(distance, radius):
                continue
            if cosine < 1.0 - 1.0e-8:
                if not shared:
                    nonparallel.append({"Distance": round(distance, 9), "AngleDegrees": round(math.degrees(math.acos(min(1.0, cosine))), 6)})
                continue
            if shared:
                continue
            separation = abs(float((b["Start"] - a["Start"]) @ a["Inward"]))
            # Overlap of the projections along the common tangent.
            t = a["Tangent"]
            ia = sorted([float(a["Start"] @ t), float(a["End"] @ t)])
            ib = sorted([float(b["Start"] @ t), float(b["End"] @ t)])
            overlap = min(ia[1], ib[1]) - max(ia[0], ib[0])
            if overlap <= 1.0e-9:
                continue
            metal_between = float((b["Start"] - a["Start"]) @ a["Inward"]) > 0 and float((a["Start"] - b["Start"]) @ b["Inward"]) > 0
            if metal_between:
                expected = "SameConductorStrip"
            elif a["Conductor"] == b["Conductor"]:
                expected = "SameConductorGap"
            else:
                expected = "DifferentConductorGap"
            quantum = 1.0e-8 * radius
            pairs.append(
                {
                    "EdgeA": a,
                    "OverlapStart": max(ia[0], ib[0]),
                    "OverlapEnd": min(ia[1], ib[1]),
                    "Separation": round(separation, 9),
                    "OverlapLength": round(overlap, 9),
                    "Expected": expected,
                    "Within2R": within_interaction(separation, radius),
                    "AtExactly2R": abs(separation - 2.0 * radius) <= 0.5 * quantum,
                    "AtExactlyR": abs(separation - radius) <= 0.5 * quantum,
                }
            )
    corner_points = [np.array(c["Point"]) for c in corners if c["ClassifierCorner"]]
    corner_pairs_within_2r = sum(1 for i in range(len(corner_points)) for j in range(i + 1, len(corner_points)) if np.linalg.norm(corner_points[i] - corner_points[j]) <= 2.0 * radius * (1 + 1e-9))
    # Design rule (SURFACE-RESPONSE-IDENTIFICATION.md (b) 3): a corner joins a cluster when
    # its 2R through-vertex zone reaches a cluster region, i.e. an event core lies within 3R.
    cores = design_event_cores(all_edges, corner_points, radius)
    for c in corners:
        if not c["ClassifierCorner"]:
            c["Standalone"] = None
            continue
        point = np.array(c["Point"])
        core_distance = float(np.linalg.norm(cores - point, axis=1).min()) if len(cores) else math.inf
        c["CoreDistance"] = None if math.isinf(core_distance) else round(core_distance, 6)
        c["Standalone"] = not core_distance < 3.0 * radius
        if not c["Standalone"]:
            c["Expected"] += " (an interaction event core within 3R: member of a spatial cluster)"
    # A parallel pair is a feature when part of its overlap survives the cluster regions
    # (distance to a core >= R) and the corner windows (R along the edge from a corner).
    for pair in pairs:
        if not pair["Within2R"]:
            pair["Survives"] = False
            continue
        a = pair["EdgeA"]
        ts = np.linspace(pair["OverlapStart"], pair["OverlapEnd"], 201)
        points = a["Start"][None, :] + (ts[:, None] - float(a["Start"] @ a["Tangent"])) * a["Tangent"][None, :]
        survives = np.ones(len(ts), dtype=bool)
        if len(cores):
            survives &= np.linalg.norm(points[:, None, :] - cores[None, :, :], axis=2).min(axis=1) >= radius
        for end in (a["Start"], a["End"]):
            if any(np.linalg.norm(end - q) < 1.0e-9 for q in corner_points):
                survives &= np.linalg.norm(points - end, axis=1) >= radius
        pair["Survives"] = bool(survives.any())
        del pair["EdgeA"]
    return {
        "Layout": lay["Name"],
        "Notes": lay["Notes"],
        "MatchingRadius": radius,
        "PhysicalPerimeterLength": physical_length,
        "TruncationLength": truncation_length,
        "StraightEdges": len(all_edges),
        "Corners": corners,
        "ClassifierCornerCount": sum(1 for c in corners if c["ClassifierCorner"]),
        "StandaloneCornerCount": sum(1 for c in corners if c.get("Standalone")),
        "ExpectedRoundedCorners": sum(1 for a in arcs if a["ExpectedRoundedCorner"]),
        "SubThresholdTurns": sum(1 for c in corners if not c["ClassifierCorner"]),
        "CornerPairsWithin2R": corner_pairs_within_2r,
        "ParallelPairs": pairs,
        "NonparallelInteractions": nonparallel,
        "Arcs": arcs,
        "Excluded": excluded,
    }


# -------------------------------------------------------------------------------- runner ----


def generate_meshes(layouts, directory, julia="julia"):
    os.makedirs(directory, exist_ok=True)
    spec = write_specification(layouts, os.path.join(directory, "layouts.spec"))
    log = os.path.join(directory, "generate.log")
    started = time.time()
    with open(log, "w") as target:
        completed = subprocess.run([julia, GENERATOR, spec, directory], stdout=target, stderr=subprocess.STDOUT)
    return completed.returncode, time.time() - started, log


def run_preflight(palace, config_path, ranks, log_path, timeout=600.0):
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    for name in list(environment):
        if name.startswith("PALACE_RESPONSE_"):
            del environment[name]
    with open(log_path, "w") as log:
        completed = subprocess.run([palace, "-np", str(ranks), "--surface-response-preflight", config_path], stdout=log, stderr=subprocess.STDOUT, env=environment, timeout=timeout)
    return completed.returncode


def compare_with_oracle(orc, audit_result, manifest):
    """Oracle vs the audit's mesh perimeter and vs the manifest: per-invariant pass / fail."""
    census = audit_result["Perimeter"]
    checks = {}
    mesh_length = census["LengthByClass"].get("PHYSICAL", 0.0)
    has_arcs = bool(orc["Arcs"])
    length_tolerance = 1.0e-6 * max(1.0, orc["PhysicalPerimeterLength"]) if not has_arcs else 0.02 * sum(a["Length"] for a in orc["Arcs"])
    checks["A6-perimeter-length"] = {"Oracle": orc["PhysicalPerimeterLength"], "Mesh": mesh_length, "Tolerance": length_tolerance, "Pass": abs(mesh_length - orc["PhysicalPerimeterLength"]) <= length_tolerance, "Meaning": "polyline perimeter exact; arcs shorter by the chord defect"}
    oracle_corners = [c for c in orc["Corners"] if c["ClassifierCorner"]]
    mesh_corners = [c for c in census["Corners"] if c["Kind"] == "CORNER"]
    matched = 0
    angle_mismatch = []
    unmatched_mesh = []
    for mc in mesh_corners:
        point = np.array(mc["Point"][:2])
        hit = next((oc for oc in oracle_corners if np.linalg.norm(np.array(oc["Point"]) - point) < 1.0e-6), None)
        if hit is None:
            unmatched_mesh.append(mc)
        elif abs(hit["InteriorAngleDegrees"] - mc["InteriorAngleDegrees"]) > 1.0e-6 and abs(hit["TurnDegrees"] - mc["TurnDegrees"]) > 1.0e-6:
            angle_mismatch.append({"Oracle": hit, "Mesh": mc})
        else:
            matched += 1
    checks["A6-corner-list-mesh"] = {
        "OracleCorners": len(oracle_corners),
        "MeshCorners": len(mesh_corners),
        "Matched": matched,
        "AngleMismatch": angle_mismatch,
        "UnmatchedMesh": unmatched_mesh[:20],
        "Pass": matched == len(oracle_corners) == len(mesh_corners) and not angle_mismatch,
        "Meaning": "mesh corners beyond the oracle's are arc discretisation vertices above the 30 deg tolerance",
    }
    summary = M.summarize(manifest)
    manifest_sharp = Counter()
    manifest_rounded = Counter()
    manifest_topologies = Counter()
    for r in manifest["Requirements"]:
        manifest_topologies[r["Topology"]] += int(r["Count"])
        if r["Topology"] in ("ConvexCorner", "ConcaveCorner"):
            angle = round(float(r["Geometry"].get("AngleDegrees", float("nan"))), 6)
            corner_radius = float(r["Geometry"].get("CornerRadius", 0.0))
            if corner_radius > 0.0:
                manifest_rounded[f"{r['Topology']}@{angle} radius {corner_radius}"] += int(r["Count"])
            else:
                manifest_sharp[(r["Topology"], angle)] += int(r["Count"])
    standalone = [c for c in oracle_corners if c["Standalone"]]
    oracle_sharp = Counter((c["Expected"], round(c["InteriorAngleDegrees"] if c["Convex"] else 360.0 - c["InteriorAngleDegrees"], 6)) for c in standalone)
    oracle_rounded = Counter(a["Expected"] for a in orc["Arcs"] if a["ExpectedRoundedCorner"])
    checks["A6-corner-list-manifest"] = {
        "OracleClassifierCorners": len(oracle_corners),
        "OracleStandaloneCorners": len(standalone),
        "OracleAbsorbedByClusters": len(oracle_corners) - len(standalone),
        "OracleSharpByAngle": {f"{k[0]}@{k[1]}": v for k, v in sorted(oracle_sharp.items())},
        "ManifestSharpByAngle": {f"{k[0]}@{k[1]}": v for k, v in sorted(manifest_sharp.items())},
        "OracleRounded": dict(oracle_rounded),
        "ManifestRounded": dict(manifest_rounded),
        "ManifestClusters": manifest_topologies.get("SpatialEdgeCluster", 0),
        "Pass": manifest_sharp == oracle_sharp and manifest_rounded == oracle_rounded,
        "Meaning": "standalone corners (no interaction event core within 3R) must appear as sharp corner records with the wedge angle; fillet arcs with radius < R as rounded corners; corners with a core within 3R are members of spatial clusters",
    }
    # Quantized decisions (kDecisionLengthQuantum): a separation within half a quantum of 2R is
    # AT the threshold and not within it, so pairs at exactly 2R are expected to be isolated.
    oracle_pairs = Counter((p["Expected"], round(p["Separation"], 6)) for p in orc["ParallelPairs"] if p["Within2R"] and p.get("Survives", True))
    mesh_separations = census["Interactions"]["ParallelSeparations"]
    manifest_pairs = Counter()
    for r in manifest["Requirements"]:
        if r["Topology"] in M.TRANSLATIONAL and r["Topology"] != "IsolatedEdge":
            manifest_pairs[(r["Topology"], round(float(r["Geometry"].get("Separation", r["Geometry"].get("Width", float("nan")))), 6))] += 1
    expected_pairs_set = {(k[0], k[1]) for k in oracle_pairs}
    manifest_pairs_set = set(manifest_pairs)
    checks["A6-parallel-pairs"] = {
        "OraclePairs": {f"{k[0]}@{k[1]}": v for k, v in sorted(oracle_pairs.items())},
        "OracleAtExactly2R": sum(1 for p in orc["ParallelPairs"] if p["AtExactly2R"]),
        "OracleAtExactlyR": sum(1 for p in orc["ParallelPairs"] if p["AtExactlyR"]),
        "MeshParallelSeparations": mesh_separations,
        "ManifestPairClasses": {f"{k[0]}@{k[1]}": v for k, v in sorted(manifest_pairs.items())},
        "ManifestTopologies": dict(manifest_topologies),
        "Pass": expected_pairs_set == manifest_pairs_set,
        "Meaning": "expected class and separation of every parallel pair within 2R vs the manifest's translational records; pairs at exactly 2R are the knife edge",
    }
    checks["A6-nonparallel-interactions"] = {"Oracle": len(orc["NonparallelInteractions"]), "Mesh": census["Interactions"]["NonparallelPairs"], "Pass": None, "Meaning": "record only: non-parallel pairs are omitted by the classifier"}
    checks["A6-excluded-classes"] = {
        "Oracle": orc["Excluded"],
        "MeshLengthByClass": {k: census["LengthByClass"].get(k, 0.0) for k in ("NONPLANAR", "CROSS_LAYER", "NONMANIFOLD")},
        "Pass": (orc["Excluded"]["CrossLayerSheets"] > 0) == (census["LengthByClass"].get("CROSS_LAYER", 0.0) > 0) and (orc["Excluded"]["Walls"] > 0) == (census["LengthByClass"].get("NONPLANAR", 0.0) > 0),
    }
    return checks


def run_suite(args):
    layouts = stress_suite()
    if args.names:
        layouts = [lay for lay in layouts if lay["Name"] in set(args.names)]
    os.makedirs(args.output, exist_ok=True)
    mesh_directory = os.path.join(args.output, "meshes")
    missing = [lay for lay in layouts if not os.path.exists(os.path.join(mesh_directory, lay["Name"] + ".msh2"))]
    if missing and not args.no_generate:
        code, seconds, log = generate_meshes(missing, mesh_directory, args.julia)
        print(f"generated {len(missing)} meshes: exit {code} in {seconds:.1f} s ({log})", flush=True)
    results = []
    for lay in layouts:
        mesh_path = os.path.join(mesh_directory, lay["Name"] + ".msh2")
        orc = oracle(lay)
        record = {"Layout": lay["Name"], "Notes": lay["Notes"], "Oracle": orc, "Cells": [], "Mesh": mesh_path}
        with open(os.path.join(args.output, lay["Name"] + ".oracle.json"), "w") as target:
            json.dump(orc, target, indent=1, default=str)
        if not os.path.exists(mesh_path):
            record["Error"] = "mesh not generated"
            results.append(record)
            print(json.dumps({"Layout": lay["Name"], "Error": record["Error"]}), flush=True)
            continue
        mesh = read_msh2(mesh_path)
        record["MeshNodes"] = int(len(mesh.coordinates))
        attributes = sorted({int(t) for t in mesh.physical_tags(2)})
        ground = [GROUND] if GROUND in attributes else []
        terminals = [[SECOND_CONDUCTOR]] if SECOND_CONDUCTOR in attributes else []
        sa = [SUBSTRATE_AIR] if SUBSTRATE_AIR in attributes else []
        digests = {}
        for label, library in args.libraries.items():
            for levels in args.uniform_levels:
                for ranks in args.ranks:
                    name = f"{label}-u{levels}-np{ranks}"
                    directory = os.path.join(args.output, lay["Name"], name)
                    os.makedirs(directory, exist_ok=True)
                    config = preflight_config(mesh_path, ground, terminals, sa, library, os.path.join(directory, "postpro"), uniform_levels=levels)
                    config_path = os.path.join(directory, "config.json")
                    with open(config_path, "w") as target:
                        json.dump(config, target, indent=2)
                    log_path = os.path.join(directory, "palace.log")
                    manifest_path = os.path.join(directory, "postpro", "surface-response-requirements.json")
                    cell = {"Name": name, "Library": label, "UniformLevels": levels, "Ranks": ranks}
                    if not (os.path.exists(manifest_path) and not args.rerun):
                        try:
                            cell["ExitCode"] = run_preflight(args.palace, config_path, ranks, log_path, args.timeout)
                        except subprocess.TimeoutExpired:
                            cell["ExitCode"] = "timeout"
                    else:
                        cell["ExitCode"] = 0
                    if os.path.exists(manifest_path):
                        manifest = M.load_manifest(manifest_path)
                        cell["DigestFull"] = M.canonical_digest(manifest)
                        cell["DigestGeometryOnly"] = M.canonical_digest(manifest, geometry_only_counts=True)
                        cell["GeometryDigest"] = manifest.get("Identification", {}).get("GeometryDigest")
                        cell["Log"] = M.parse_palace_log(log_path)
                        # Version 2: the identification's GeometryDigest is the A3 / A4 / A5 identity.
                        digests.setdefault((label, levels), {})[ranks] = cell["GeometryDigest"] or cell["DigestFull"]
                        if levels == 0:
                            audit_args = argparse.Namespace(mesh=mesh_path, config=config_path, manifest=manifest_path, log=log_path, compare=None, radius=None, corner_tolerance=CORNER_TURN_TOLERANCE_DEGREES, output_prefix=None)
                            result = audit.run_audit(audit_args)
                            with open(os.path.join(directory, "audit.json"), "w") as target:
                                json.dump(result, target, indent=1, default=str)
                            with open(os.path.join(directory, "audit.md"), "w") as target:
                                target.write(audit.render_markdown(result))
                            cell["Gates"] = {g["Gate"]: g["Status"] for g in result["Gates"]}
                            cell["OracleChecks"] = compare_with_oracle(orc, result, manifest)
                            cell["ManifestSummary"] = {k: {"Count": v["Count"], "Length": round(v["TotalEdgeLength"], 6), "Missing": v["Missing"]} for k, v in M.summarize(manifest)["ByTopology"].items()}
                            cell["Clusters"] = M.summarize(manifest)["Clusters"]
                    else:
                        cell["Error"] = "no manifest"
                        if os.path.exists(log_path):
                            with open(log_path, errors="replace") as source:
                                text = source.read()
                            marker = text.find("Verification failed")
                            cell["Abort"] = text[marker : marker + 400] if marker >= 0 else text[-600:]
                    record["Cells"].append(cell)

        def manifest_of(label, levels, ranks):
            return M.load_manifest(os.path.join(args.output, lay["Name"], f"{label}-u{levels}-np{ranks}", "postpro", "surface-response-requirements.json"))

        # Invariants across cells: A4 ranks (per library and level), A5 refinement (per
        # library, geometry-only digest), A3 libraries (geometry-only, at level 0).
        record["A4-rank-determinism"] = all(len(set(v.values())) == 1 for v in digests.values()) if digests else None
        record["A5-refinement-invariance"] = {}
        record["A3-library-independence"] = {}
        first_ranks = args.ranks[0]
        labels = list(args.libraries)
        def feature_diff(a, b):
            ia, ib = a.get("Identification"), b.get("Identification")
            if not (ia and ib):
                return None
            ca = Counter((f["Type"], f["Hash"][:12]) for f in ia["Features"])
            cb = Counter((f["Type"], f["Hash"][:12]) for f in ib["Features"])
            return {"GeometryDigestIdentical": ia["GeometryDigest"] == ib["GeometryDigest"], "OnlyA": sorted(f"{t}:{h}x{n}" for (t, h), n in (ca - cb).items()), "OnlyB": sorted(f"{t}:{h}x{n}" for (t, h), n in (cb - ca).items())}

        for label in labels:
            if (label, 0) in digests and (label, 1) in digests:
                a, b = manifest_of(label, 0, first_ranks), manifest_of(label, 1, first_ranks)
                diff = M.diff_manifests(a, b, geometry_only_counts=True)
                record["A5-refinement-invariance"][label] = {"GeometryOnlyIdentical": diff["Identical"], "Added": [e["Topology"] for e in diff["Added"]], "Removed": [e["Topology"] for e in diff["Removed"]], "Changed": [(e["A"]["Topology"], round(e["A"]["TotalEdgeLength"], 6), round(e["B"]["TotalEdgeLength"], 6)) for e in diff["Changed"]], "Identification": feature_diff(a, b)}
        for label in labels[1:]:
            if (labels[0], 0) in digests and (label, 0) in digests:
                a, b = manifest_of(labels[0], 0, first_ranks), manifest_of(label, 0, first_ranks)
                diff = M.diff_manifests(a, b, geometry_only_counts=True)
                record["A3-library-independence"][f"{labels[0]} vs {label}"] = {"GeometryOnlyIdentical": diff["Identical"], "Added": [e["Topology"] for e in diff["Added"]], "Removed": [e["Topology"] for e in diff["Removed"]], "Changed": [(e["A"]["Topology"], round(e["A"]["TotalEdgeLength"], 6), round(e["B"]["TotalEdgeLength"], 6)) for e in diff["Changed"]], "Identification": feature_diff(a, b)}
        results.append(record)
        with open(os.path.join(args.output, lay["Name"] + ".result.json"), "w") as target:
            json.dump(record, target, indent=1, default=str)
        print(json.dumps(summary_row(record), default=str), flush=True)
    write_summary(results, args.output)
    return results


def summary_row(record):
    cells = record.get("Cells", [])
    audited = [c for c in cells if c.get("OracleChecks")]
    by_library = {}
    for c in audited:
        by_library.setdefault(c["Library"], c)
    return {
        "Layout": record["Layout"],
        "Nodes": record.get("MeshNodes"),
        "Exit": sorted({str(c.get("ExitCode")) for c in cells}),
        "A4": record.get("A4-rank-determinism"),
        "A5": {k: (v["Identification"]["GeometryDigestIdentical"] if v.get("Identification") else v["GeometryOnlyIdentical"]) for k, v in (record.get("A5-refinement-invariance") or {}).items()},
        "A3": {k: (v["Identification"]["GeometryDigestIdentical"] if v.get("Identification") else v["GeometryOnlyIdentical"]) for k, v in (record.get("A3-library-independence") or {}).items()},
        "Oracle": {label: {k.replace("A6-", ""): v["Pass"] for k, v in c["OracleChecks"].items() if v["Pass"] is not None} for label, c in by_library.items()},
        "GatesFailing": {label: [k for k, v in c["Gates"].items() if v != "PASS"] for label, c in by_library.items()},
        "Topologies": {label: {k: v["Count"] for k, v in c["ManifestSummary"].items()} for label, c in by_library.items()},
        "Error": record.get("Error"),
    }


def write_summary(results, output):
    lines = ["# Synthetic stress layouts: identification vs oracle", ""]
    lines.append("| Layout | nodes | exit | A4 ranks | A5 refine | A3 libraries | A6 oracle (per library) | audit gates failing | manifest topologies |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for r in results:
        row = summary_row(r)
        oracle_text = "; ".join(f"{label}: " + ", ".join(f"{k}={'P' if v else 'F'}" for k, v in checks.items()) for label, checks in row["Oracle"].items())
        gates_text = "; ".join(f"{label}: {', '.join(g.replace('A1-', '').replace('A2-', '') for g in gates)}" for label, gates in row["GatesFailing"].items())
        topology_text = "; ".join(f"{label}: " + ", ".join(f"{k}:{v}" for k, v in t.items()) for label, t in row["Topologies"].items())
        lines.append(f"| {row['Layout']} | {row['Nodes']} | {','.join(row['Exit'])} | {row['A4']} | {row['A5']} | {row['A3']} | {oracle_text} | {gates_text} | {topology_text} |")
    with open(os.path.join(output, "summary.md"), "w") as target:
        target.write("\n".join(lines) + "\n")
    with open(os.path.join(output, "summary.json"), "w") as target:
        json.dump([summary_row(r) for r in results], target, indent=1, default=str)
    print("\n".join(lines))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", help="output directory (meshes/, per-layout cells, summary)")
    parser.add_argument("--names", nargs="*", help="subset of layouts")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--uniform-levels", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--library", action="append", help="label=path (repeatable; default seed=transmon seed, isolated=corrected-p1 library)")
    parser.add_argument("--palace", default=DEFAULT_PALACE)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--no-generate", action="store_true")
    parser.add_argument("--rerun", action="store_true")
    parser.add_argument("--list", action="store_true", help="print the layout names and exit")
    parser.add_argument("--write-spec", help="write the specification file and exit")
    args = parser.parse_args(argv)
    if args.list:
        for lay in stress_suite():
            print(lay["Name"], "-", lay["Notes"])
        return 0
    if args.write_spec:
        write_specification(stress_suite() if not args.names else [lay for lay in stress_suite() if lay["Name"] in set(args.names)], args.write_spec)
        return 0
    if not args.output:
        parser.error("--output is required to run the suite")
    args.libraries = dict(entry.split("=", 1) for entry in args.library) if args.library else {"seed": DEFAULT_LIBRARY, "isolated": ISOLATED_LIBRARY}
    for ranks in args.ranks:
        if ranks > 6 or ranks > (os.cpu_count() or 1):
            raise SystemExit(f"{ranks} ranks exceeds the local limit")
    run_suite(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
