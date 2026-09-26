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
from collections import Counter, defaultdict

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import audit, manifest as M, perimeter as P  # noqa: E402
from .msh2 import read_msh2  # noqa: E402
from .preflight_config import preflight_config  # noqa: E402
from .signature_library import build_signature_library  # noqa: E402

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
# Curved-edge chain rule (SURFACE-RESPONSE-IDENTIFICATION.md (b) 7, Identification.Conventions):
# a bend whose radius is below STRAIGHT_BEND_RADIUS_OVER_R x R is a curved class; two chains pair
# along a bend when their closest-point separation is constant within PAIR_SEPARATION_TOLERANCE.
STRAIGHT_BEND_RADIUS_OVER_R = 10.0
PAIR_SEPARATION_TOLERANCE = 0.05


# ----------------------------------------------------------------------------- geometry ----


def loop(points, arcs=None):
    """A closed loop: points counter-clockwise for outer boundaries; arcs maps the index of an
    end vertex to the arc centre (the arc runs from the previous vertex to that vertex)."""
    return {"Points": [(float(x), float(y)) for x, y in points], "Arcs": {int(k): (float(v[0]), float(v[1])) for k, v in (arcs or {}).items()}}


def sheet(attribute, outer, holes=(), z=0.0):
    return {"Attribute": attribute, "Z": z, "Loops": [outer, *holes]}


def layout(name, sheets, walls=(), half_x=30.0, half_y=30.0, depth=15.0, height=15.0, lc_fine=1.0, lc_far=6.0, notes="", bend=None, ports=(), expected=None, frame_normal=None):
    """bend = {"Radius": centreline radius, "Width": bar width} records the design intent of a
    polyline arc bar so that the oracle can state the expected curvature class; ports = non-metal
    port faces (attribute, loop) on z = 0 written as LumpedPort boundaries of the preflight
    configuration; expected = the stated expectation of a decision-82 oracle case
    (check_expected): feature counts by type, exclusion lengths by class, vertex types."""
    return {"Name": name, "HalfX": half_x, "HalfY": half_y, "Depth": depth, "Height": height, "LcFine": lc_fine, "LcFar": lc_far, "Sheets": list(sheets), "Walls": list(walls), "Ports": [{"Attribute": a, "Loops": [lp]} for a, lp in ports], "Notes": notes, "Bend": bend, "Expected": expected, "FrameNormal": frame_normal}


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


def arc_bar(width, radius, sweep_degrees, step_degrees, lead=6.0, offset=0.0, lead_end=None, tee=None):
    """A bar following a circular arc discretised as a polyline with the given turning angle
    per vertex, with straight leads at both ends. offset shifts the bar across the centreline:
    it occupies the offsets [offset - width / 2, offset + width / 2] of the centreline polyline
    (positive = left = towards the arc centre); two bars of opposite offsets about one
    centreline face each other across a gap whose chords are exactly parallel at the design
    separation, like an offset path. lead_end is the length of the end lead (default lead);
    tee = (bar_length, bar_width) ends the bar in a perpendicular cross-bar (a T junction)."""
    steps = max(1, int(round(sweep_degrees / step_degrees)))
    step = math.radians(sweep_degrees) / steps
    h = 0.5 * width
    centre = np.array([0.0, radius])
    lead_end = lead if lead_end is None else lead_end
    centreline = [centre + radius * np.array([math.sin(k * step), -math.cos(k * step)]) for k in range(steps + 1)]
    d_start = np.array([1.0, 0.0])
    d_end = np.array([math.cos(steps * step), math.sin(steps * step)])
    centreline = [centreline[0] - lead * d_start, *centreline, centreline[-1] + lead_end * d_end]
    # Offset polyline by +-h using the vertex bisector normals (exact offset of the polyline).
    def offset_polyline(sign):
        distance = offset + sign * h
        result = []
        n = len(centreline)
        for i, p in enumerate(centreline):
            if i == 0:
                d = centreline[1] - p
                normal = np.array([-d[1], d[0]]) / np.linalg.norm(d)
                result.append(p + distance * normal)
            elif i == n - 1:
                d = p - centreline[i - 1]
                normal = np.array([-d[1], d[0]]) / np.linalg.norm(d)
                result.append(p + distance * normal)
            else:
                d0 = p - centreline[i - 1]
                d1 = centreline[i + 1] - p
                n0 = np.array([-d0[1], d0[0]]) / np.linalg.norm(d0)
                n1 = np.array([-d1[1], d1[0]]) / np.linalg.norm(d1)
                b = n0 + n1
                b /= np.linalg.norm(b)
                result.append(p + (distance / (b @ n0)) * b)
        return result

    right = offset_polyline(-1.0)
    left = offset_polyline(1.0)
    if tee is None:
        return loop([tuple(p) for p in right + left[::-1]])
    # T junction: a cross-bar of the given length and width across the end of the bar,
    # perpendicular to the end lead (the bar's right side runs into the cross-bar's right arm).
    bar_length, bar_width = tee
    end = centreline[-1] + offset * np.array([-d_end[1], d_end[0]])
    n_end = np.array([-d_end[1], d_end[0]])  # left normal of the end lead
    cross = [end - 0.5 * bar_length * n_end, end - 0.5 * bar_length * n_end + bar_width * d_end, end + 0.5 * bar_length * n_end + bar_width * d_end, end + 0.5 * bar_length * n_end]
    return loop([tuple(p) for p in right + cross + left[::-1]])


def trapezoid(bottom_width, top_width, height):
    return loop([(-0.5 * bottom_width, 0.0), (0.5 * bottom_width, 0.0), (0.5 * top_width, height), (-0.5 * top_width, height)])


def slot_shape(arm_length, half_height, width, depth):
    """U-shaped sheet: two arms of the given length either side of a slot of the given width
    and depth (open at the top, y = +half_height), joined by a bar below the slot; counter-
    clockwise."""
    hw = 0.5 * width
    x_out = arm_length + hw
    y_top = half_height
    y_slot = y_top - depth
    y_bottom = y_slot - half_height
    return loop([(-x_out, y_bottom), (x_out, y_bottom), (x_out, y_top), (hw, y_top), (hw, y_slot), (-hw, y_slot), (-hw, y_top), (-x_out, y_top)])


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
        # Gap between the two arms of one connected sheet (a slot of width s and depth 12 in a
        # U-shaped ground): conductor identity is metal connectivity (phase 3), so two disjoint
        # rectangles would be different conductors even under one attribute.
        layouts.append(layout(f"gap-same-{tag}", [sheet(GROUND, slot_shape(14.0, 6.0, s, 12.0))], notes=f"SameConductorGap at separation {s}: slot in one connected sheet (the slot ends are corner clusters below 2R)"))
        # Gap between two disjoint sheets (different conductors by connectivity; the second
        # carries its own attribute so the label-based tagging agrees).
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
                layout(f"arc-r{radius:g}-step{step:g}", [sheet(GROUND, arc_bar(3.0, radius, sweep, step))], half_x=half, half_y=half, lc_fine=lc, lc_far=max(6.0, half / 5.0), notes=f"bar of width 3 on a polyline arc of radius {radius}, sweep {sweep} deg, {step} deg per vertex", bend={"Radius": radius, "Width": 3.0})
            )
    layouts.append(layout("arc-r20-step5-fine", [sheet(GROUND, arc_bar(3.0, 20.0, 90.0, 5.0))], half_x=36.0, half_y=36.0, lc_fine=0.5, lc_far=6.0, notes="arc-r20-step5 at half the mesh size (mesh independence)", bend={"Radius": 20.0, "Width": 3.0}))
    # A gap of exactly 2R (and 2R +/- 1e-3 R) between two concentric 8 um bars (4 R: the far
    # corners of a bar end are beyond the 3R vertex-join reach of the near corners) along
    # bends of 50 and 250 um at three discretisations: the interaction decision uses the curve
    # separation, so every discretisation gives the straight-pair answer (exactly 2R and
    # 2R + 1e-3 R: isolated edges; 2R - 1e-3 R: a DifferentConductorGap). Mid-chord dips of
    # the polylines below 2R must not create events (DS-SCT-001: 4 um CPW gaps at R = 2 um).
    for radius, sweep in [(50.0, 45.0), (250.0, 15.0)]:
        for step in [1.0, 5.0, 15.0]:
            for gap_tag, gap in [("2R", 2.0 * RADIUS), ("2Rminus", 2.0 * RADIUS - 1.0e-3 * RADIUS), ("2Rplus", 2.0 * RADIUS + 1.0e-3 * RADIUS)]:
                width = 8.0
                extent = (radius + 0.5 * gap + width) * math.sin(math.radians(sweep)) + 12.0
                half = max(30.0, math.ceil(extent + 4.0))
                lc = 1.0 if radius < 100 else 2.0
                # Both bars are offsets of the gap's centreline (radius): the facing edges are
                # exactly parallel chords at the design gap along the bend (an offset path).
                layouts.append(
                    layout(
                        f"gap-bend-r{radius:g}-{gap_tag}-step{step:g}",
                        [sheet(GROUND, arc_bar(width, radius, sweep, step, offset=0.5 * gap + 0.5 * width)), sheet(GROUND, arc_bar(width, radius, sweep, step, offset=-0.5 * gap - 0.5 * width))],
                        half_x=half, half_y=half, lc_fine=lc, lc_far=max(6.0, half / 5.0),
                        notes=f"two concentric 8 um bars with a gap of {gap:g} ({gap_tag}) along a polyline arc of radius {radius}, sweep {sweep} deg, {step} deg per vertex: the straight-pair answer at every discretisation",
                        bend={"Radius": radius, "Width": width, "Gap": gap},
                    )
                )
    layouts.append(layout("taper-10", [sheet(GROUND, trapezoid(24.0, 24.0 - 2.0 * 20.0 * math.tan(math.radians(10.0)), 20.0))], notes="trapezoid with two 10 deg taper edges: corners 80 and 100 deg"))
    # Pairs along bends with divergent ends (the local constancy rule): an 8 um centre bar on a
    # 50 um bend with two 8 um flanks at a 3 um gap, ending in a T junction (the flanks stop
    # 3 um before the cross-bar: the gap pairs end there, the T is a cluster on each side) and
    # ending at a port cut (every bar reaches the truncation box). Tapers: a slow one (gap
    # 3.0 -> 3.15 um, 5 % over 20 R: locally constant, one pair) and a fast one (a bar at 80
    # deg converging to 1 um: events, a cluster, no pair).
    for step in [5.0, 15.0]:
        centre = arc_bar(8.0, 50.0, 90.0, step, lead=6.0, lead_end=9.0, tee=(40.0, 6.0))
        flanks = [arc_bar(8.0, 50.0, 90.0, step, lead=6.0, lead_end=6.0, offset=sign * (3.0 + 8.0)) for sign in (1.0, -1.0)]
        layouts.append(layout(f"cpw-tee-r50-step{step:g}", [sheet(GROUND, centre), sheet(GROUND, flanks[0]), sheet(GROUND, flanks[1])], half_x=84.0, half_y=90.0, lc_fine=1.0, lc_far=12.0, notes=f"8 um centre bar on a 50 um bend ({step} deg per vertex) with 8 um flanks at a 3 um gap ending in a T junction: gap pairs along the bend; each flank end has two corners 3 um below the cross-bar (its side edges meet the cross-bar's bottom edge at 90 deg): two clusters per side", bend={"Radius": 50.0, "Width": 8.0, "Gap": 3.0, "Clusters": 4}))
        half = 60.0
        port = [arc_bar(8.0, 50.0, 45.0, step, lead=half, lead_end=6.0, offset=o) for o in (0.0, 11.0, -11.0)]
        layouts.append(layout(f"cpw-port-r50-step{step:g}", [sheet(GROUND, p) for p in port], half_x=half, half_y=half, lc_fine=1.0, lc_far=12.0, notes=f"8 um centre bar on a 50 um bend ({step} deg per vertex) with 8 um flanks at a 3 um gap, every bar cut by the truncation box at x = -{half:g} (a port cut): gap pairs up to the cut", bend={"Radius": 50.0, "Width": 8.0, "Gap": 3.0}))
    layouts.append(layout("taper-slow", [sheet(GROUND, rectangle(-20.0, -8.0, 20.0, -2.0)), sheet(GROUND, loop([(-20.0, 1.0), (20.0, 1.15), (20.0, 7.15), (-20.0, 7.0)]))], notes="two 6 um bars whose gap tapers from 3.0 to 3.15 um over 40 um (5 % over 20 R): locally constant, one DifferentConductorGap; the bar ends face each other (two corner-pair clusters)", bend={"Gap": 3.075, "Clusters": 2}))
    d80 = np.array([math.cos(math.radians(80.0)), math.sin(math.radians(80.0))])
    n80 = np.array([-d80[1], d80[0]])
    p0 = np.array([10.0, -1.0])
    layouts.append(layout("taper-fast", [sheet(GROUND, rectangle(-20.0, -6.0, 20.0, -2.0)), sheet(GROUND, loop([tuple(p0), tuple(p0 + 20.0 * d80), tuple(p0 + 20.0 * d80 + 4.0 * n80), tuple(p0 + 4.0 * n80)]))], notes="a 4 um bar at 80 deg converging to 1 um above a straight bar: no constant separation, a cluster at the convergence (the bar's corner within 2R of the straight edge), no pair", bend={"Clusters": 1}))
    layouts.append(layout("facing-layers", [sheet(GROUND, rectangle(-10.0, -6.0, 10.0, 6.0)), sheet(GROUND, rectangle(-10.0, -6.0, 10.0, 6.0), z=3.0)], notes="two facing metal sheets 3 um apart in one PEC attribute: cross-layer class, excluded by decision 73(3)"))
    layouts.append(layout("vertical-wall", [sheet(GROUND, rectangle(-10.0, -6.0, 10.0, 6.0))], walls=[(GROUND, 0.0, -6.0, 0.0, 6.0, 4.0)], notes="a vertical metal wall standing across the sheet: non-planar and non-manifold classes"))
    layouts += decision_82_suite()
    layouts.append(layout("island-rounded-8x6", [sheet(GROUND, rounded_rectangle(4.0, 3.0, 0.5))], half_x=12.0, half_y=12.0, depth=8.0, height=8.0, lc_fine=0.25, lc_far=2.0, notes="survey geometry 2: 8 x 6 um island with 0.5 um fillets (perimeter 24 + pi)"))
    layouts.append(layout("aperture-rounded-8x6", [sheet(GROUND, rectangle(-12.0, -12.0, 12.0, 12.0), holes=[rounded_rectangle(4.0, 3.0, 0.5)])], half_x=12.0, half_y=12.0, depth=8.0, height=8.0, lc_fine=0.25, lc_far=2.0, notes="survey geometry 2, aperture variant: the same shape as a hole in a ground plane"))
    return layouts


PORT = 9


def fillet_points(corner, d_in, d_out, radius, chords, perturb=None):
    """Inscribed polyline of the fillet replacing the sharp corner at `corner` between the
    incoming direction d_in and the outgoing direction d_out (unit vectors): the tangent
    points and `chords` chords on the circle of the given radius tangent to both arms.
    `perturb` = (rng, relative) moves the interior chord vertices radially by up to the
    relative fraction of the radius (a seeded mesh perturbation; the tangent points stay)."""
    d_in = np.asarray(d_in, dtype=float)
    d_out = np.asarray(d_out, dtype=float)
    corner = np.asarray(corner, dtype=float)
    turn = math.atan2(d_in[0] * d_out[1] - d_in[1] * d_out[0], d_in @ d_out)  # signed
    t = radius * math.tan(abs(turn) / 2.0)  # tangent length
    ta = corner - t * d_in
    tb = corner + t * d_out
    n_in = np.array([-d_in[1], d_in[0]]) * (1.0 if turn > 0 else -1.0)  # toward the centre
    centre = ta + radius * n_in
    points = []
    for k in range(chords + 1):
        phi = abs(turn) * k / chords
        # rotate (ta - centre) about the centre by phi in the turn direction
        v = ta - centre
        c, s_ = math.cos(phi), math.sin(phi)
        sign = 1.0 if turn > 0 else -1.0
        r = np.array([c * v[0] - sign * s_ * v[1], sign * s_ * v[0] + c * v[1]])
        if perturb is not None and 0 < k < chords:
            rng, relative = perturb
            r = r * (1.0 + relative * (2.0 * rng.random() - 1.0))
        points.append(tuple(centre + r))
    assert np.allclose(points[-1], tb, atol=1.0e-9 * max(1.0, radius))
    return points


def filleted_polygon(vertices, radius, chords, perturb=None, fillet=None):
    """Replace every corner of the counter-clockwise polygon (or those where fillet(i) is
    true) by an inscribed fillet polyline of `chords` chords."""
    n = len(vertices)
    out = []
    for i in range(n):
        p = np.asarray(vertices[i], dtype=float)
        if fillet is not None and not fillet(i):
            out.append(tuple(p))
            continue
        a = np.asarray(vertices[(i - 1) % n], dtype=float)
        b = np.asarray(vertices[(i + 1) % n], dtype=float)
        d_in = (p - a) / np.linalg.norm(p - a)
        d_out = (b - p) / np.linalg.norm(b - p)
        out.extend(fillet_points(p, d_in, d_out, radius, chords, perturb))
    return loop(out)


FILLET_RATIOS = [0.1, 0.25, 0.5, 0.9, 1.1, 2.0, 5.0, 9.0, 11.0, 20.0]
FILLET_TURNS = [45, 90, 135, 180]
FILLET_CHORDS = [2, 4, 8, 16]


def fillet_chord_length(radius, turn_degrees, chords):
    return 2.0 * radius * math.sin(math.radians(turn_degrees) / (2.0 * chords))


def fillet_suite(ratios=FILLET_RATIOS, turns=FILLET_TURNS, chords=FILLET_CHORDS, perturbed_chords=8):
    """Decision 82(3) mesh-independence gate: a bar (width 4R) bent by the turn with both bend
    corners filleted at rho = ratio x R (45 / 90 / 135 deg), or a strip of width 2 rho ending in
    a semicircle (180 deg), each meshed with 2 / 4 / 8 / 16 chords per fillet plus a seeded
    radial perturbation (1 % of the chord) of the chord vertices at 8 chords; and the DS-SCT-002
    cross-shaped pattern (1 um fillets = 0.5 R with two 45 deg chords). Expectation per case:
    rho < R -> one rounded corner per fillet (total turn, rho / R), no CurvedEdge; R <= rho <
    10 R -> one CurvedEdge per fillet with RadiusOverR = rho / R exactly; rho >= 10 R ->
    straight-like isolated edges; the isolated / curved features are one per chain group so
    the counts do not depend on the chords. The gate holds for the chord counts whose chord is
    shorter than 2R (a coarser polyline is a different geometry at the scale of R: its kinks
    are real corners)."""
    R = RADIUS
    layouts = []
    for q in ratios:
        rho = q * R
        for turn in turns:
            # Perturbation: 1 % of the chord length, radially (a mesh-noise scale: the joint
            # turns keep their sign; a perturbation on the scale of rho is a different arc).
            variants = [(n, None) for n in chords] + [(perturbed_chords, 0.01 * fillet_chord_length(rho, turn, perturbed_chords) / rho)]
            for n, perturbation in variants:
                tag = f"fillet-q{q:g}-t{turn}-c{n}" + ("-perturbed" if perturbation else "")
                perturb = (np.random.default_rng(20260925), perturbation) if perturbation else None
                if turn < 180:
                    width = 4.0 * R
                    tangent = rho * math.tan(math.radians(turn) / 2.0)
                    arm = max(16.0, tangent + 4.0 * R + 6.0)
                    theta = 180.0 - turn
                    base = bent_bar(width, arm, theta)["Points"]
                    # corners 1 (outer, convex) and 4 (inner, concave) are the bend corners
                    poly = filleted_polygon(base, rho, n, perturb, fillet=lambda i: i in (1, 4))
                    xs = [p[0] for p in poly["Points"]]; ys = [p[1] for p in poly["Points"]]
                    half_x = math.ceil(max(abs(v) for v in xs) + 6.0); half_y = math.ceil(max(abs(v) for v in ys) + 6.0)
                    corner_rounded = q < 1.0
                    # Acute bend (interior 45 deg): the arms interact beyond the through-arc zone
                    # (points at a < R / sin(22.5 deg) = 2.61 R from the virtual corner are within 2R of
                    # each other; the zone covers a <= rho tan(67.5 deg) + 2R), i.e. for
                    # rho / R < (2.61 - 2) / 2.414 = 0.254: then the multiset is the sharp corner-45
                    # layout's (2 clusters absorb two bar-end corners), the two rounded corners intact.
                    acute_clusters = turn == 135 and 2.414 * q + 2.0 < 1.0 / math.sin(math.radians(22.5))
                    if corner_rounded:
                        features = {"ConvexCorner": 3, "ConcaveCorner": 1, "IsolatedEdge": 4, "SpatialEdgeCluster": 2} if acute_clusters else {"ConvexCorner": 5, "ConcaveCorner": 1, "IsolatedEdge": 6}
                    elif q < 10.0:
                        features = {"ConvexCorner": 4, "CurvedEdge": 2, "IsolatedEdge": 4}
                    else:
                        features = {"ConvexCorner": 4, "IsolatedEdge": 4}
                    expected = {
                        # a rounded corner separates its arms like a sharp one (two isolated edges per side); a bend continues the edge (one)
                        "Features": features,
                        "Exclusions": {},
                        "CornerSignatures": {f"ConvexCorner@{theta:g}@{q:g}": 1, f"ConcaveCorner@{theta:g}@{q:g}": 1} if corner_rounded else {},
                        "CurvedRadii": {f"{q:g}": 2} if (1.0 <= q < 10.0) else {},
                        "BendAnnotation": q if q >= 10.0 else None,
                    }
                    lay = layout(tag, [sheet(GROUND, poly)], half_x=half_x, half_y=half_y, lc_fine=min(1.0, max(0.05, rho)), lc_far=max(6.0, half_x / 5.0), notes=f"bar of width 4R bent by {turn} deg, both bend corners filleted with rho = {q:g} R as {n} chords" + (" (chord vertices perturbed radially by 1 % of the chord)" if perturbation else ""), expected=expected)
                else:
                    # U-turn: strip of width 2 rho from the box edge x = -half_x to x = 0, semicircular end.
                    half_x = math.ceil(max(30.0, 6.0 * R + 2.0 * rho)); half_y = math.ceil(max(30.0, rho + 8.0 * R))
                    # Semicircle about (0, 0) from (0, -rho) through (rho, 0) to (0, rho), n chords.
                    arc = []
                    for k in range(n + 1):
                        phi = -math.pi / 2.0 + math.pi * k / n
                        r = rho
                        if perturb is not None and 0 < k < n:
                            r *= 1.0 + perturb[1] * (2.0 * perturb[0].random() - 1.0)
                        arc.append((r * math.cos(phi), r * math.sin(phi)))
                    poly = loop([(-half_x, -rho)] + arc + [(-half_x, rho)])
                    corner_rounded = q < 1.0
                    expected = {
                        "Features": ({"ConvexCorner": 1, "SameConductorStrip": 1} if corner_rounded else ({"CurvedEdge": 1, "IsolatedEdge": 1} if q < 10.0 else {"IsolatedEdge": 1})),
                        "CornerSignatures": {f"ConvexCorner@0@{q:g}": 1} if corner_rounded else {},
                        "CurvedRadii": {f"{q:g}": 1} if (1.0 <= q < 10.0) else {},
                        "BendAnnotation": q if q >= 10.0 else None,
                        "Vertices": {"TruncationCut": 2},
                    }
                    lay = layout(tag, [sheet(GROUND, poly)], half_x=half_x, half_y=half_y, lc_fine=min(1.0, max(0.05, rho)), lc_far=max(6.0, half_x / 5.0), notes=f"strip of width 2 rho = {2 * q:g} R ending in a semicircle of {n} chords (a U-turn)" + (" (chord vertices perturbed radially by 1 % of the chord)" if perturbation else ""), expected=expected)
                lay["Fillet"] = {"RatioOverR": q, "TurnDegrees": turn, "Chords": n, "Perturbed": bool(perturbation), "ChordOverR": fillet_chord_length(rho, turn, n) / R}
                layouts.append(lay)
    # DS-SCT-002: a cross of 6 um arms with 1 um fillets meshed as two 45 deg chords.
    cross = cross_shape(12.0, 6.0)["Points"]
    for n in (2, 4):
        poly = filleted_polygon(cross, 1.0, n)
        layouts.append(layout(f"fillet-cross-sct002-c{n}", [sheet(GROUND, poly)], lc_fine=0.5, notes=f"DS-SCT-002 pattern: cross of 6 um arms, every corner a 1 um fillet (0.5 R) as {n} chords", expected={
            # every rounded corner separates its arms: 12 rounded corners (0.5 R) and 12 straight edges
            # between them, of which the 4 arm-end edges (6 - 2 x 1 um) lie inside the two adjacent
            # corner windows (2 x R): 8 isolated-edge features of 9 - 2 - 4 = 3 um
            "Features": {"ConvexCorner": 8, "ConcaveCorner": 4, "IsolatedEdge": 8},
            "CornerSignatures": {"ConvexCorner@90@0.5": 8, "ConcaveCorner@90@0.5": 4},
            "CurvedRadii": {},
            "Exclusions": {},
        }))
    return layouts


def decision_82_suite():
    """Oracle cases of the decision-82 rules with their expectation stated before running
    (`Expected`, checked by check_expected): (1) one interaction distance and no cross-plane
    features — two 20 x 12 um sheets (perimeter 64 um each) on planes 1.5R / 2.0R / 2.4R apart,
    overlapping and offset by 8 um in x, plus two slot patterns (corner clusters) 2.4R apart
    whose sites are within the former 3R vertex join but beyond 2R; (5) a lumped port bridging
    the 4 um gap between two 6 um leads."""
    layouts = []
    R = RADIUS
    box = rectangle(-10.0, -6.0, 10.0, 6.0)
    one_sheet = {"ConvexCorner": 4, "IsolatedEdge": 4}
    # EdgeFrameNormal +z on every interface so that the sheet in the vacuum (one material on
    # both sides) is identified like the flip chip's top metal instead of being
    # UndeterminedProcessSide.
    up = [0.0, 0.0, 1.0]
    for tag, z in [("1p5R", 1.5 * R), ("2R", 2.0 * R), ("2p4R", 2.4 * R)]:
        within = z < 2.0 * R
        # Overlapping: every edge point of either sheet is at distance z from the other sheet's
        # face -> CrossLayer over the whole perimeter of both sheets below 2R (strict), nothing
        # at 2R and beyond: two independent rectangles, no cross-plane pair or cluster.
        layouts.append(layout(f"planes-{tag}-overlap", [sheet(GROUND, box), sheet(GROUND, box, z=z)], frame_normal=up, notes=f"two identical 20 x 12 sheets on planes {z:g} um = {z / R:g} R apart: {'whole perimeter CrossLayer' if within else 'two independent rectangles (no cross-plane feature)'}", expected={
            "Features": {} if within else {k: 2 * v for k, v in one_sheet.items()},
            "Exclusions": {"CrossLayer": 128.0} if within else {},
            "ExcludedVertices": 8 if within else 0,
            "CrossPlaneFeatures": 0,
            "PlanesRule": "no feature may hold portions on both planes",
        }))
        # Offset by 8 um in x: the top sheet covers x in [-2, 18]. Below 2R every bottom-edge
        # point within sqrt((2R)^2 - z^2) in plan view of the top sheet's face is CrossLayer
        # (analytic zone boundary), the far ends stay planar features; at and beyond 2R nothing.
        reach = math.sqrt(max(0.0, (2.0 * R) ** 2 - z * z))
        zone = 2.0 * (10.0 + 2.0 + reach) + 12.0  # two long edges from x = -2 - reach to 10, the end edge
        layouts.append(layout(f"planes-{tag}-offset", [sheet(GROUND, box), sheet(GROUND, rectangle(-2.0, -6.0, 18.0, 6.0), z=z)], half_x=36.0, frame_normal=up, notes=f"two 20 x 12 sheets offset by 8 um on planes {z / R:g} R apart: {'CrossLayer zones where the sheets overlap in plan view (+ the in-plane reach), planar features on the far ends' if within else 'two independent rectangles'}", expected={
            # below 2R: the far-end corners of both sheets and, per sheet, the end edge + the
            # planar remainders of the two long chains (one IsolatedEdge feature per chain)
            "Features": {"ConvexCorner": 4, "IsolatedEdge": 6} if within else {k: 2 * v for k, v in one_sheet.items()},
            "Exclusions": {"CrossLayer": 2.0 * zone} if within else {},
            "ExcludedVertices": 4 if within else 0,
            "CrossPlaneFeatures": 0,
            "PlanesRule": "no feature may hold portions on both planes",
            "Tolerance": 0.05,
        }))
    # Two slot patterns (a 3 um slot in a U: its two 90 deg concave corners are 3 um < 2R apart
    # -> one corner cluster per slot end, as gap-same-3) on planes 2.4R apart, the top one shifted
    # by 2 um in y: the top cluster's sites are 5.2 um = 2.6R (< the former 3R join) from the
    # bottom cores in 3D, beyond 2R. Expected: the clusters of the two planes stay separate.
    z = 2.4 * R
    shifted = slot_shape(14.0, 6.0, 3.0, 12.0)
    shifted = loop([(x, y + 2.0) for x, y in shifted["Points"]])
    layouts.append(layout("planes-2p4R-clusters", [sheet(GROUND, slot_shape(14.0, 6.0, 3.0, 12.0)), sheet(GROUND, shifted, z=z)], frame_normal=up, notes="two 3 um slots (each: one SameConductorGap, two corner clusters at the slot end and mouth, 4 convex corners, 5 isolated edges, as gap-same-*) on planes 2.4R apart, the top one shifted by 2 um in y: its sites are 2.6R from the bottom cores (inside the former 3R join, beyond 2R) and its slot edges lie above the bottom slot's (a cross-plane translational pair before the plane rule)", expected={
        "Features": {"SpatialEdgeCluster": 4, "SameConductorGap": 2, "ConvexCorner": 8, "IsolatedEdge": 10},
        "CrossPlaneFeatures": 0,
        "Exclusions": {},
        "ExcludedVertices": 0,
    }))
    # A lumped port bridging the 4 um gap between two 6 um wide leads (x in [-20, -2] and
    # [2, 20]): the lead ends bordering the port (2 x 6 um) are the Port exclusion, their
    # corners are PortCut, the leads' long edges run to the cut; the far ends keep their corners.
    layouts.append(layout("port-bridge", [sheet(GROUND, rectangle(-20.0, -3.0, -2.0, 3.0)), sheet(GROUND, rectangle(2.0, -3.0, 20.0, 3.0))], ports=[(PORT, rectangle(-2.0, -3.0, 2.0, 3.0))], notes="a lumped port face bridging the 4 um gap between two 6 um leads: the lead ends along the port are the Port exclusion (12 um), their four vertices PortCut, no corner / endpoint / pair there", expected={
        "Features": {"ConvexCorner": 4, "IsolatedEdge": 6},
        "Exclusions": {"Port": 12.0},
        "Vertices": {"PortCut": 4, "ConvexCorner": 4},
    }))
    return layouts


def arc_band(r_in, r_out, sweep_degrees=90.0, lead=12.0, cap=True):
    """A band between the radii r_in < r_out along a circular path about the origin from the
    angle -90 deg (the point (0, -r), heading +x) counter-clockwise over sweep_degrees, with
    straight tangent leads of the given length at both ends; true circular arcs (Gmsh) split
    into pieces below 180 deg. Counter-clockwise outline. cap=False leaves the end edges to be
    cut by the truncation box (the caller sizes the box)."""
    sweep = math.radians(sweep_degrees)
    a0 = -0.5 * math.pi
    a1 = a0 + sweep
    t0 = np.array([math.cos(a0 + 0.5 * math.pi), math.sin(a0 + 0.5 * math.pi)])  # heading at the start
    t1 = np.array([math.cos(a1 + 0.5 * math.pi), math.sin(a1 + 0.5 * math.pi)])
    def at(r, a):
        return (r * math.cos(a), r * math.sin(a))
    pieces = max(1, int(math.ceil(sweep_degrees / 90.0 - 1.0e-9)))
    points, arcs = [], {}
    # Outer side forward: start lead, arc, end lead.
    p = np.array(at(r_out, a0))
    points.append(tuple(p - lead * t0))
    points.append(tuple(p))
    for k in range(1, pieces + 1):
        points.append(at(r_out, a0 + sweep * k / pieces))
        arcs[len(points)] = (0.0, 0.0)  # 1-based: the arc ends at this vertex
    q = np.array(at(r_out, a1))
    points.append(tuple(q + lead * t1))
    # Inner side backward.
    q = np.array(at(r_in, a1))
    points.append(tuple(q + lead * t1))
    points.append(tuple(q))
    for k in range(pieces - 1, -1, -1):
        points.append(at(r_in, a0 + sweep * k / pieces))
        arcs[len(points)] = (0.0, 0.0)
    p = np.array(at(r_in, a0))
    points.append(tuple(p - lead * t0))
    return loop(points, arcs)


def hairpin(rho, gap, half_y):
    """A strip of width w = 2 rho - gap folded through a semicircle of centreline radius rho
    about the origin (the fold below y = 0), its two legs running up to y = half_y (the
    truncation): the inner edges of the legs face each other across `gap` at x = +-gap / 2, the
    inner fold has radius gap / 2 and the outer fold 2 rho - gap / 2."""
    r_in = 0.5 * gap
    r_out = 2.0 * rho - 0.5 * gap
    points = [(-r_out, half_y), (-r_out, 0.0), (0.0, -r_out), (r_out, 0.0), (r_out, half_y),
              (r_in, half_y), (r_in, 0.0), (0.0, -r_in), (-r_in, 0.0), (-r_in, half_y)]
    arcs = {3: (0.0, 0.0), 4: (0.0, 0.0), 8: (0.0, 0.0), 9: (0.0, 0.0)}  # 1-based end vertices
    return loop(points, arcs)


def kinked_hairpin(gap, width, half_y, turns_degrees=(22.0, 27.0, 24.0, 29.0, 23.0, 28.0, 27.0), chords=(0.35, 0.65, 0.45, 0.8, 0.5, 0.7)):
    """A strip whose inner edge folds through 180 deg along a NON-circular polyline (unequal
    chords and turns, every turn below the corner threshold, the total 180 deg): no arc fits
    it (the arc rule's 5 % circle test fails), so the fold stays inside the chain and the two
    inner legs, `gap` apart, are one chain facing itself — the self-pairing case. The outer
    edge is the exact offset polyline at `width` (> 2R: no strip pair). The fold's chords are
    scaled so that the polyline closes on the leg separation."""
    heading = -0.5 * math.pi
    pts = [np.array([0.0, 0.0])]
    turns = [t * 180.0 / sum(turns_degrees) for t in turns_degrees]
    for c, t in zip(chords, turns[:-1]):
        heading += math.radians(t)
        pts.append(pts[-1] + c * np.array([math.cos(heading), math.sin(heading)]))
    scale = gap / (pts[-1][0] - pts[0][0])
    pts = [p * scale for p in pts]
    shift = np.array([-0.5 * gap, 0.0]) - pts[0]
    inner = [np.array([-0.5 * gap, half_y])] + [p + shift for p in pts] + [np.array([0.5 * gap, half_y])]
    # Exact offset polyline by width to the right of the travel direction (away from the slot).
    outer = []
    n = len(inner)
    for i, p in enumerate(inner):
        if i == 0 or i == n - 1:
            d = inner[1] - inner[0] if i == 0 else inner[-1] - inner[-2]
            d = d / np.linalg.norm(d)
            outer.append(p + width * np.array([d[1], -d[0]]))
        else:
            d0 = p - inner[i - 1]
            d1 = inner[i + 1] - p
            n0 = np.array([d0[1], -d0[0]]) / np.linalg.norm(d0)
            n1 = np.array([d1[1], -d1[0]]) / np.linalg.norm(d1)
            b = n0 + n1
            b = b / np.linalg.norm(b)
            outer.append(p + (width / (b @ n0)) * b)
    return loop([tuple(p) for p in outer] + [tuple(p) for p in inner[::-1]])


def stack_suite():
    """Decision 82(2) oracle cases (multi-edge cross-section stacks) with their expectation
    stated before running. Straight stacks: bars along x from -30 to 30 (end corners inside
    the box: one corner cluster per end) between ground bars that reach the truncation box
    (no ground corners); curved stacks: bands along a 90 deg arc with 12 um leads (arc_band);
    the taper, the U-ring, the wide ground, the hairpins and the kinked (self-facing) hairpin.
    R = 2 um: a 2 um gap / trace is 1R, 4 um exactly 2R (no interaction, strict rule)."""
    R = RADIUS
    layouts = []
    L, half = 30.0, 40.0

    def straight(name, traces, ground_below_gap, ground_above_gap, offsets, isolated, notes, strips=0, extra=None):
        """strips: recomposed SameConductorStrip features next to the end clusters (the trace
        alone where its ground partners are cluster material while the trace edges, exactly R
        from the cores, are not); extra: further recomposed features at the stack ends
        ({"ParallelEdgeCluster": [offsets, ...], "<PairType>": count})."""
        sheets = []
        y = 0.0
        for width, gap in traces:
            sheets.append(sheet(GROUND, rectangle(-L, y, L, y + width)))
            y += width + gap
        y_top = y - traces[-1][1]
        if ground_below_gap is not None:
            sheets.append(sheet(GROUND, rectangle(-half, -8.0 - ground_below_gap, half, -ground_below_gap)))
        if ground_above_gap is not None:
            sheets.append(sheet(GROUND, rectangle(-half, y_top + ground_above_gap, half, y_top + ground_above_gap + 8.0)))
        extra = dict(extra or {})
        extra_stacks = extra.pop("ParallelEdgeCluster", [])
        features = {"ParallelEdgeCluster": 1 + len(extra_stacks), "SpatialEdgeCluster": 2, "IsolatedEdge": isolated, **extra}
        if strips:
            features["SameConductorStrip"] = strips
        expected = {"Features": features, "Offsets": {"ParallelEdgeCluster": [offsets] + extra_stacks}, "FacingGates": True}
        layouts.append(layout(f"stack-{name}", sheets, half_x=half, half_y=half, lc_fine=1.0, lc_far=6.0, notes=notes, expected=expected))

    # k = 3: ground | 2 | trace 2: edges G-top, T-bottom, T-top at 0 / 1R / 2R (the trace's top
    # edge faces nothing but is the strip partner of its bottom edge). Isolated edges are one
    # feature per chain: the ground's far edge and its near edge beyond the end clusters (2 per
    # ground). Where the end cluster holds the ground edge (R around its cores) but not the
    # trace edges (exactly R away: the knife edge), the trace pair is recomposed as a
    # SameConductorStrip (the stack-end rule); with a ground on both sides the trace edges
    # between the two ground claims are recomposed the same way. Decision 85(2) (2026-09-26):
    # a member alone next to a cluster (its partners cluster material, itself within 2R of
    # the cluster's claimed perimeter across) JOINS the cluster instead of standing as an
    # isolated ClusterNeighbour edge — the isolated counts below are restated accordingly.
    straight("k3-2-2", [(2.0, 0.0)], 2.0, None, [0.0, 1.0, 2.0], 2, "ground | 2 um gap | 2 um trace: a 3-edge stack at 0 / 1 / 2 R over the straight run, one corner cluster per trace end, the ground's far edge and its near edge beyond the clusters isolated; the trace strip recomposed next to the clusters", strips=1)
    straight("k3-1p5-3", [(3.0, 0.0)], 1.5, None, [0.0, 0.75, 2.25], 2, "ground | 1.5 um gap | 3 um trace: 3-edge stack at 0 / 0.75 / 2.25 R (asymmetric); the trace's bottom edge is 0.75 R from the ground's cluster cores (inside the end cluster over a shorter reach than the ground itself), its top edge 1.5 R from them: the top edge alone next to the clusters joins them (decision 85(2); was an isolated ClusterNeighbour), the trace strip recomposed between the two claim ends", strips=1)
    straight("k3-3-1", [(1.0, 0.0)], 3.0, None, [0.0, 1.5, 2.0], 2, "ground | 3 um gap | 1 um trace: 3-edge stack at 0 / 1.5 / 2 R", strips=1)
    straight("k4-2-2-2", [(2.0, 0.0)], 2.0, 2.0, [0.0, 1.0, 2.0, 3.0], 4, "ground | 2 | trace 2 | 2 | ground (the DS-SCT-001 flux line): a 4-edge stack at 0 / 1 / 2 / 3 R; isolated: two far edges + the near edges beyond the end clusters", strips=1)
    straight("k4-1-1p5-3", [(1.5, 0.0)], 1.0, 3.0, [0.0, 0.5, 1.25, 2.75], 4, "ground | 1 | trace 1.5 | 3 | ground: asymmetric 4-edge stack at 0 / 0.5 / 1.25 / 2.75 R; stack ends (cluster claims of unequal reach: lower ground to 2 R + sqrt(15) um, trace bottom 0.5 R from its cores, upper ground to 2 R + sqrt(7) um): a 3-edge stack (trace | 3 | ground) then a DifferentConductorGap (trace top | ground); the trace top alone beyond joins the cluster (decision 85(2))", strips=0, extra={"ParallelEdgeCluster": [[0.0, 0.75, 2.25]], "DifferentConductorGap": 1})
    straight("k5-2-2-2-2", [(2.0, 2.0), (2.0, 0.0)], 2.0, None, [0.0, 1.0, 2.0, 3.0, 4.0], 2, "ground | 2 | trace 2 | 2 | trace 2: 5-edge stack at 0 .. 4 R; each trace's end edge makes events on the other trace's facing long edge (within 2R of it, not through-vertex), so the inner edges are cluster material to 2 R + sqrt(12) um while the outer edges (exactly R from the cores) join the clusters (decision 85(2); were isolated ClusterNeighbours): isolated = the ground's two edges", strips=0)
    straight("k6-2-2-2-2-2", [(2.0, 2.0), (2.0, 0.0)], 2.0, 2.0, [0.0, 1.0, 2.0, 3.0, 4.0, 5.0], 4, "ground | 2 | trace | 2 | trace | 2 | ground: 6-edge stack at 0 .. 5 R; the traces' inner edges and both grounds are cluster material at the ends, the traces' outer edges there join the clusters (decision 85(2)): isolated = the grounds' four edges", strips=0)
    # Around 2R: the second gap at 3.9 um joins the stack (k = 4), at exactly 4.0 = 2R and at
    # 4.1 um it does not (k = 3, the upper ground's near edge isolated; its end corners are
    # beyond 2R of the trace corners, so it takes no part in the end clusters).
    straight("k4-2-2-3p9", [(2.0, 0.0)], 2.0, 3.9, [0.0, 1.0, 2.0, 3.95], 4, "ground | 2 | trace 2 | 3.9 | ground: the 3.9 um gap (< 2R) makes a 4-edge stack at 0 / 1 / 2 / 3.95 R; the upper ground's cluster claim is shorter (its cores within 2R of the trace end edge span sqrt(16 - 3.9^2) um): a 3-edge stack trace | 3.9 | ground between the two claims, the trace strip alone before it", strips=1, extra={"ParallelEdgeCluster": [[0.0, 1.0, 2.95]]})
    straight("k3-2-2-4p0", [(2.0, 0.0)], 2.0, 4.0, [0.0, 1.0, 2.0], 4, "ground | 2 | trace 2 | 4.0 = 2R | ground: the knife edge — exactly 2R does not interact: a 3-edge stack, the upper ground's edges isolated", strips=1)
    straight("k3-2-2-4p1", [(2.0, 0.0)], 2.0, 4.1, [0.0, 1.0, 2.0], 4, "ground | 2 | trace 2 | 4.1 | ground: 3-edge stack, the upper ground isolated", strips=1)
    straight("k3-wide-ground", [(2.0, 0.0)], 2.0, 6.0, [0.0, 1.0, 2.0], 4, "a 3-edge stack (ground | 2 | trace 2) next to a wide ground edge 6 um = 3R above the trace: the wide edge is isolated over its whole length (no third-edge facing)", strips=1)

    # Curved stacks: bands about the origin from -90 deg over 90 deg with 12 um leads; the
    # innermost edge at r0 = rho (the tightest bend of the stack, RadiusOverR = rho / R for the
    # curved class when rho < 10 R). Traces 3 um wide (1.5 R), gaps 2 um; the ground band 8 um.
    for q in (3.0, 8.0, 30.0):
        rho = q * R
        curved = q < STRAIGHT_BEND_RADIUS_OVER_R
        r_far = rho + 3.0 + 2.0 + 8.0
        far_curved = r_far < STRAIGHT_BEND_RADIUS_OVER_R * R
        extent = r_far + 12.0 + 4.0
        hx = math.ceil(extent) + 2.0
        # k = 3: trace + ground: stacks straight on the two leads and curved along the arc (or
        # one straight stack when rho >= 10 R); the ground's far edge: isolated leads + curved
        # arc (or one isolated edge); its end edges beyond the corner cluster and the far
        # corner window: one isolated piece each; far corners convex.
        # The straight stack is ONE feature on both leads (one feature per signature, class and
        # member chains); the ground's far edge one isolated chain (+ a curved section) and
        # its two end edges isolated between the cluster and the far corner's window; the
        # trace's inner edge next to the end clusters (its partner, the trace's outer edge,
        # inside the cluster; itself exactly R from the cores) JOINS the clusters (decision
        # 85(2); was an isolated ClusterNeighbour): isolated = far edge + two end edges.
        k3 = [0.0, 1.5, 2.5]
        stacks = {"ParallelEdgeCluster": 1, "CurvedParallelEdgeCluster": 1} if curved else {"ParallelEdgeCluster": 1}
        far = ({"IsolatedEdge": 1 + 2, "CurvedEdge": 1} if far_curved else {"IsolatedEdge": 1 + 2})
        expected = {"Features": {**stacks, "SpatialEdgeCluster": 2, "ConvexCorner": 2, **far},
                    "Offsets": {t: [k3] * n for t, n in stacks.items()}, "FacingGates": True}
        if curved:
            expected["StackRadii"] = {"CurvedParallelEdgeCluster": [q]}
        layouts.append(layout(f"stack-curved-k3-rho{q:g}", [sheet(GROUND, arc_band(rho, rho + 3.0)), sheet(GROUND, arc_band(rho + 5.0, rho + 13.0))], half_x=hx, half_y=hx, lc_fine=1.0, lc_far=max(6.0, hx / 6.0), notes=f"3 um trace (inner radius {q:g} R) | 2 um gap | 8 um ground band along a 90 deg bend with 12 um leads: 3-edge stack at 0 / 1.5 / 2.5 R, {'curved along the bend (RadiusOverR ' + f'{q:g}' + ') and straight on the leads' if curved else 'straight-like throughout (one feature, bend annotation)'}", bend={"Radius": rho}, expected=expected))
        # k = 4: two 3 um traces 2 um apart, nothing else: 0 / 1.5 / 2.5 / 4 R; the outer
        # edges next to the end clusters join them (decision 85(2)): nothing isolated.
        k4 = [0.0, 1.5, 2.5, 4.0]
        expected = {"Features": {**stacks, "SpatialEdgeCluster": 2}, "Offsets": {t: [k4] * n for t, n in stacks.items()}, "FacingGates": True}
        if curved:
            expected["StackRadii"] = {"CurvedParallelEdgeCluster": [q]}
        hx4 = math.ceil(rho + 8.0 + 12.0 + 4.0) + 2.0
        layouts.append(layout(f"stack-curved-k4-rho{q:g}", [sheet(GROUND, arc_band(rho, rho + 3.0)), sheet(GROUND, arc_band(rho + 5.0, rho + 8.0))], half_x=hx4, half_y=hx4, lc_fine=1.0, lc_far=max(6.0, hx4 / 6.0), notes=f"two 3 um traces 2 um apart along a 90 deg bend (inner radius {q:g} R): 4-edge stack at 0 / 1.5 / 2.5 / 4 R", bend={"Radius": rho}, expected=expected))
        # k = 5: two traces + ground: 0 / 1.5 / 2.5 / 4 / 5 R.
        k5 = [0.0, 1.5, 2.5, 4.0, 5.0]
        r_far5 = rho + 8.0 + 2.0 + 8.0
        far_curved5 = r_far5 < STRAIGHT_BEND_RADIUS_OVER_R * R
        far = ({"IsolatedEdge": 1 + 2, "CurvedEdge": 1} if far_curved5 else {"IsolatedEdge": 1 + 2})
        expected = {"Features": {**stacks, "SpatialEdgeCluster": 2, "ConvexCorner": 2, **far}, "Offsets": {t: [k5] * n for t, n in stacks.items()}, "FacingGates": True}
        if curved:
            expected["StackRadii"] = {"CurvedParallelEdgeCluster": [q]}
        hx5 = math.ceil(r_far5 + 12.0 + 4.0) + 2.0
        layouts.append(layout(f"stack-curved-k5-rho{q:g}", [sheet(GROUND, arc_band(rho, rho + 3.0)), sheet(GROUND, arc_band(rho + 5.0, rho + 8.0)), sheet(GROUND, arc_band(rho + 10.0, rho + 18.0))], half_x=hx5, half_y=hx5, lc_fine=1.0, lc_far=max(6.0, hx5 / 6.0), notes=f"two 3 um traces and an 8 um ground band, 2 um gaps, along a 90 deg bend (inner radius {q:g} R): 5-edge stack at 0 / 1.5 / 2.5 / 4 / 5 R", bend={"Radius": rho}, expected=expected))
        # k = 6: three traces: 0 / 1.5 / 2.5 / 4 / 5 / 6.5 R; the outer edges next to the
        # end clusters join them (decision 85(2)): nothing isolated.
        k6 = [0.0, 1.5, 2.5, 4.0, 5.0, 6.5]
        expected = {"Features": {**stacks, "SpatialEdgeCluster": 2}, "Offsets": {t: [k6] * n for t, n in stacks.items()}, "FacingGates": True}
        if curved:
            expected["StackRadii"] = {"CurvedParallelEdgeCluster": [q]}
        hx6 = math.ceil(rho + 13.0 + 12.0 + 4.0) + 2.0
        layouts.append(layout(f"stack-curved-k6-rho{q:g}", [sheet(GROUND, arc_band(rho, rho + 3.0)), sheet(GROUND, arc_band(rho + 5.0, rho + 8.0)), sheet(GROUND, arc_band(rho + 10.0, rho + 13.0))], half_x=hx6, half_y=hx6, lc_fine=1.0, lc_far=max(6.0, hx6 / 6.0), notes=f"three 3 um traces 2 um apart along a 90 deg bend (inner radius {q:g} R): 6-edge stack at 0 / 1.5 / 2.5 / 4 / 5 / 6.5 R", bend={"Radius": rho}, expected=expected))

    # Taper: a 2 um trace between ground edges whose gap tapers from 6 um (3R: isolated ground
    # edges, the trace's own strip pair) to 2 um (a 4-edge stack) through 8 deg kinks over 14.2 um
    # per side (variation 14 % per R: not locally constant -> one cluster over the taper joining
    # both sides through the trace). All bars reach the truncation in x.
    dx = 4.0 / math.tan(math.radians(8.0))
    ground_lower = loop([(-half, -14.0), (half, -14.0), (half, -2.0), (dx / 2, -2.0), (-dx / 2, -6.0), (-half, -6.0)])
    ground_upper = loop([(-half, 8.0), (-dx / 2, 8.0), (dx / 2, 4.0), (half, 4.0), (half, 16.0), (-half, 16.0)])
    layouts.append(layout("stack-taper-wide-to-narrow", [sheet(GROUND, rectangle(-half, 0.0, half, 2.0)), sheet(GROUND, ground_lower), sheet(GROUND, ground_upper)], half_x=half, half_y=half, lc_fine=1.0, lc_far=6.0, notes="2 um trace between ground edges tapering from a 6 um gap (isolated ground edges + the trace strip pair at 1R) to 2 um (4-edge stack at 0 / 1 / 2 / 3 R) through 8 deg kinks: the taper is not locally constant and is one cluster; the wide-end kinks read curved (8 deg over 1 R: radius 7 R)", expected={
        "FeaturesSubset": {"ParallelEdgeCluster": 1, "SameConductorStrip": 1, "SpatialEdgeCluster": 1},
        "Offsets": {"ParallelEdgeCluster": [[0.0, 1.0, 2.0, 3.0]]}, "FacingGates": True}))

    # U-ring (the DS-SCT-001 / DS-SCT-002 flux-line end): a 2 um trace from the top truncation
    # into a rounded-rectangle ring (strip 2 um; hole 12 x 10 um with 1 um = 0.5 R corners, outer
    # corners 3 um), the ground 2 um outside (corners 5 um) and along the feed. Stacks: the feed
    # (4 edges at 0 / 1 / 2 / 3 R) and the ring's left / right / bottom sides (3 edges: ground |
    # ring outer | ring inner at 0 / 1 / 2 R) between the rounded corners' windows; the bends
    # are curved gaps (ground / ring outer, radii 5 and 3 um), the inner corners rounded
    # concave corners, the junction a cluster.
    ring_hole = rounded_rectangle(6.0, 5.0, 1.0)
    ring_outer_pts = rounded_rectangle(8.0, 7.0, 3.0)
    def with_feed(rr, x_feed, y_top):
        """Insert the feed (x = +-x_feed up to y_top) into the top edge of a rounded rectangle."""
        pts = rr["Points"]; arcs = rr["Arcs"]  # arcs keyed by 1-based end vertex
        # The top edge runs from the 4th point (right-top arc end) to the 5th (left-top arc start).
        out = []; new_arcs = {}
        for k, p in enumerate(pts, start=1):
            out.append(p)
            if k in arcs:
                new_arcs[len(out)] = arcs[k]
            if k == 4:
                out += [(x_feed, p[1]), (x_feed, y_top), (-x_feed, y_top), (-x_feed, p[1])]
        return loop(out, new_arcs)
    ring_metal = with_feed(ring_outer_pts, 1.0, half)
    # The ground: the box minus the feed slot and the ring's surroundings as ONE loop (a hole
    # touching the truncation boundary is not a valid sheet): clockwise around the ring.
    gr = rounded_rectangle(10.0, 9.0, 5.0)
    g_pts, g_arcs = gr["Points"], gr["Arcs"]  # 1-based arc end keys
    ground_pts = [(-half, -half), (half, -half), (half, half), (3.0, half), (3.0, 9.0)]
    ground_arcs = {}
    # Around the ring clockwise: from the 4th point (5, 9) back through 3, 2, 1, 8, 7, 6, 5.
    order = [4, 3, 2, 1, 8, 7, 6, 5]
    for i, k in enumerate(order):
        ground_pts.append(g_pts[k - 1])
        # Going backwards, the arc that ended at vertex k (key k) now runs from k to k - 1:
        # it ends at the NEXT appended vertex.
        if k in g_arcs and i + 1 < len(order):
            ground_arcs[len(ground_pts) + 1] = g_arcs[k]
    ground_pts += [(-3.0, 9.0), (-3.0, half), (-half, half)]
    ground_metal = loop(ground_pts, ground_arcs)
    layouts.append(layout("stack-u-ring", [sheet(GROUND, ring_metal, holes=[ring_hole]), sheet(GROUND, ground_metal)], half_x=half, half_y=half, lc_fine=0.5, lc_far=6.0, notes="flux-line end: 2 um trace feeding a 2 um ring (hole 12 x 10 um, 0.5 R inner corners, 3 um outer corners) inside a ground at 2 um (5 um corners): the feed is a 4-edge stack, the ring's three free sides 3-edge stacks between the rounded corners' windows, curved ground / ring gaps around the bends, four rounded concave corners, the junction a cluster; recomposed 2-edge gaps where the inner edge is in a corner window", expected={
        "FeaturesSubset": {"ParallelEdgeCluster": 4},
        "Offsets": {"ParallelEdgeCluster": [[0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0], [0.0, 1.0, 2.0]]}, "FacingGates": True}))

    # Hairpins (supervisor addition): a strip of width 2 rho - g folded through a semicircle of
    # centreline radius rho, legs to the truncation; inner fold radius g / 2 (a rounded concave
    # corner of 180 deg turn below R, a bend at and above R), outer fold 2 rho - g / 2.
    for q in (1.5, 3.0, 8.0):
        rho = q * R
        for gq in (1.0, 1.8, 2.2):
            g = gq * R
            w = 2.0 * rho - g
            r_in, r_out = 0.5 * g, 2.0 * rho - 0.5 * g
            hx = math.ceil(r_out) + (9.0 if q < 8 else 14.0)
            # The fold's outer arc reaches y = -r_out: the box must hold it (review m8: at rho =
            # 8 R the 31 um fold crossed the former 30 um half-height and Gmsh's boolean left a
            # mesh MFEM could not load).
            hy = max(30.0, hx)
            corner = r_in < R
            outer_curved = r_out < STRAIGHT_BEND_RADIUS_OVER_R * R
            strip = w < 2.0 * R
            gap_pair = g < 2.0 * R
            features = Counter()
            notes = []
            if corner:
                # Decision 85(2): the outer fold within 2R of the corner's arc (strip w = 2 rho - g
                # below 2R; at w = 2R exactly the inscribed outer chords dip below 2R: the
                # recorded knife edge, mesh dependent) joins the corner, which is then a cluster.
                outer_joins = w < 2.0 * R + 1.0e-9
                features["SpatialEdgeCluster" if outer_joins else "ConcaveCorner"] += 1
                notes.append(f"inner fold {r_in / R:g} R < R: one rounded concave corner (180 deg turn) claiming the arc and R along each leg" + ("; the outer fold within 2R of its arc joins it: a cluster (decision 85(2))" if outer_joins else ""))
                # The legs are rigid runs: the translational rule pairs them from the corner's
                # window on (no through-vertex zone for exactly parallel rigid runs; the corner
                # window is the only exclusion) as a gap (w >= 2R) or, with the strip pair, a
                # 4-edge stack outer | strip | inner | gap | inner | strip | outer.
                if strip:
                    features["ParallelEdgeCluster"] += 1
                    notes.append(f"4-edge stack at 0 / {w / R:g} / {(w + g) / R:g} / {(2 * w + g) / R:g} R from the corner window on; the outer legs opposite the window are isolated (one outer chain)")
                else:
                    features["SameConductorGap"] += 1
                    notes.append(f"inner legs {gq:g} R apart: one SameConductorGap from the corner window on")
                if not (outer_joins and strip):
                    features["IsolatedEdge"] += 1  # the outer chain (both legs)
                if outer_curved:
                    if not (outer_joins and w < 1.5 * R):
                        features["CurvedEdge"] += 1
                        notes.append(f"outer fold {r_out / R:g} R: a curved edge (concentric with the corner arc" + (", its chords within 2R of the arc absorbed by the cluster" if outer_joins else ", a vertex neighbour") + ")")
                    else:
                        notes.append(f"outer fold {r_out / R:g} R entirely within 2R of the corner arc: absorbed by the cluster; the outer legs are stack members")
                else:
                    notes.append(f"outer chain straight-like ({r_out / R:g} R >= 10 R): one isolated edge incl. the fold")
            else:
                if strip:
                    features["SameConductorStrip"] += 1
                    features["CurvedSameConductorStrip"] += 1
                    notes.append(f"strip {w / R:g} R pairs along the whole U (straight legs one feature, curved fold radius {r_in / R:g} R)")
                else:
                    features["IsolatedEdge"] += 2  # inner and outer chains (legs)
                    features["CurvedEdge"] += 1 + (1 if outer_curved else 0)
                    notes.append(f"legs {gq:g} R apart beyond 2R: inner chain isolated legs + curved fold ({r_in / R:g} R); outer chain isolated legs" + (f" + curved fold ({r_out / R:g} R)" if outer_curved else f" (fold {r_out / R:g} R straight-like)"))
                    notes.append("the inner fold faces itself within 2R across the U: SelfNeighbourhood (below pi R of arc length), no pair")
            expected = {"Features": dict(features), "FacingGates": True}
            if corner and strip:
                expected["Offsets"] = {"ParallelEdgeCluster": [[0.0, w / R, (w + g) / R, (2 * w + g) / R]]}
            layouts.append(layout(f"hairpin-rho{q:g}-g{gq:g}".replace(".", "p"), [sheet(GROUND, hairpin(rho, g, hy))], half_x=hx, half_y=hy, lc_fine=(0.5 if r_in < 2.5 else 1.0) if q < 8 else 0.7, lc_far=6.0 if q < 8 else 5.0, notes=f"hairpin strip of width {w:g} um (centreline radius {q:g} R, legs {gq:g} R apart): " + "; ".join(notes), bend={"Radius": rho}, expected=expected))

    # The self-facing chain: a kinked hairpin whose inner fold is a non-circular polyline
    # (no arc fits: the fold stays in the chain) with the legs 1.5 R apart: the inner chain
    # pairs with itself beyond the pi R neighbourhood (SameConductorGap), the fold a curved
    # edge; the outer edge (offset 2.5 R: no strip pair) isolated legs + a curved fold.
    layouts.append(layout("hairpin-kinked-g1p5", [sheet(GROUND, kinked_hairpin(1.5 * R, 2.5 * R, 30.0))], half_x=20.0, half_y=30.0, lc_fine=0.3, lc_far=6.0, notes="kinked hairpin: the inner edge folds through 180 deg along a non-circular polyline (turns 18-32 deg, unequal chords) with the legs 1.5 R apart: one chain facing itself -> SameConductorGap beyond pi R of arc length; the fold's end faces the far leg within 2R beyond pi R of arc at a non-constant separation: events, a cluster (decision 82(2) self events point-wise, 2026-09-26); the rest of the fold a CurvedEdge; outer edge 2.5 R away: two isolated legs + a curved fold", expected={
        "FeaturesSubset": {"SameConductorGap": 1, "SpatialEdgeCluster": 1}, "FacingGates": True}))
    return layouts


def check_expected(manifest, expected, radius=RADIUS):
    """Compare a manifest with a stated expectation: feature counts by type, exclusion lengths by
    class (relative tolerance `Tolerance`, default 1e-6 of the perimeter), vertex types, and the
    plane rule (no feature with portions on two planes)."""
    ident = manifest["Identification"]
    checks = {}
    features = Counter(f["Type"] for f in ident["Features"])
    if "Features" in expected:
        checks["Features"] = {"Expected": dict(expected["Features"]), "Manifest": dict(features), "Pass": features == Counter(expected["Features"])}
    if "Exclusions" in expected:
        recorded = defaultdict(float)
        for e in ident["Exclusions"]:
            if e["Class"] != "TruncationCut":
                recorded[e["Class"]] += float(e["Length"])
        tolerance = expected.get("Tolerance", 1.0e-6) * max(float(ident["Totals"]["PerimeterLength"]), radius)
        ok = set(recorded) == set(expected["Exclusions"]) and all(abs(recorded[k] - v) <= tolerance for k, v in expected["Exclusions"].items())
        checks["Exclusions"] = {"Expected": dict(expected["Exclusions"]), "Manifest": dict(recorded), "Pass": ok}
    if "Vertices" in expected:
        types = Counter(v["Type"] for v in ident["Vertices"])
        checks["Vertices"] = {"Expected": dict(expected["Vertices"]), "Manifest": dict(types), "Pass": all(types.get(k, 0) == v for k, v in expected["Vertices"].items())}
    if "ExcludedVertices" in expected:
        n = sum(1 for v in ident["Vertices"] if v["Type"] == "Excluded")
        checks["ExcludedVertices"] = {"Expected": expected["ExcludedVertices"], "Manifest": n, "Pass": n == expected["ExcludedVertices"]}
    if "FeaturesSubset" in expected:
        checks["FeaturesSubset"] = {"Expected": dict(expected["FeaturesSubset"]), "Manifest": {k: features.get(k, 0) for k in expected["FeaturesSubset"]}, "Pass": all(features.get(k, 0) == v for k, v in expected["FeaturesSubset"].items())}
    if "Offsets" in expected:
        # Offsets / R of every pair / stack feature of the type (sorted lists), compared as a
        # multiset with a tolerance of 0.02 R (the chord readings of concentric polylines).
        result = {}
        ok = True
        def same(a, b):
            return len(a) == len(b) and all(abs(x - y) <= 0.02 for x, y in zip(a, b))
        for ftype, wanted in expected["Offsets"].items():
            found = [sorted(e["OffsetOverR"] for e in f["Signature"]["Edges"]) for f in ident["Features"] if f["Type"] == ftype]
            # The signature is canonical: either orientation of the stated offsets may be the
            # one serialised (the lexicographically smaller); accept both, greedily matched.
            remaining = list(found)
            match = len(found) == len(wanted)
            for v in wanted:
                v = sorted(v)
                mirrored = sorted(v[-1] - x for x in v)
                hit = next((i for i, a in enumerate(remaining) if same(a, v) or same(a, mirrored)), None)
                if hit is None:
                    match = False
                else:
                    remaining.pop(hit)
            result[ftype] = {"Expected": [sorted(v) for v in wanted], "Manifest": found, "Pass": match}
            ok = ok and match
        checks["Offsets"] = {"ByType": result, "Pass": ok}
    if "StackRadii" in expected:
        result = {}
        ok = True
        for ftype, wanted in expected["StackRadii"].items():
            found = sorted(f["Signature"]["RadiusOverR"] for f in ident["Features"] if f["Type"] == ftype)
            match = len(found) == len(wanted) and all(abs(a - b) <= 0.05 * b for a, b in zip(found, sorted(wanted)))
            result[ftype] = {"Expected": sorted(wanted), "Manifest": found, "Pass": match}
            ok = ok and match
        checks["StackRadii"] = {"ByType": result, "Pass": ok}
    if "CornerSignatures" in expected:
        found = Counter(f"{f['Type']}@{f['Signature']['AngleDegrees']:g}@{f['Signature']['CornerRadiusOverR']:g}" for f in ident["Features"] if f["Type"] in ("ConvexCorner", "ConcaveCorner") and f["Signature"].get("CornerRadiusOverR", 0.0) > 0.0)
        checks["CornerSignatures"] = {"Expected": dict(expected["CornerSignatures"]), "Manifest": dict(found), "Pass": found == Counter(expected["CornerSignatures"])}
    if "CurvedRadii" in expected:
        found = Counter(f"{f['Signature']['RadiusOverR']:g}" for f in ident["Features"] if f["Type"] == "CurvedEdge")
        checks["CurvedRadii"] = {"Expected": dict(expected["CurvedRadii"]), "Manifest": dict(found), "Pass": found == Counter(expected["CurvedRadii"])}
    if expected.get("BendAnnotation") is not None:
        annotations = [f["BendRadiusOverR"] for f in ident["Features"] if f["Type"] == "IsolatedEdge" and f.get("BendRadiusOverR") is not None]
        target = expected["BendAnnotation"]
        checks["BendAnnotation"] = {"Expected": target, "Manifest": annotations, "Pass": bool(annotations) and all(abs(a - target) <= 0.02 * target for a in annotations)}
    if "CrossPlaneFeatures" in expected:
        normal = np.array(ident["ReferenceProcessNormal"], dtype=float)
        offsets = [round(float(np.array(s["Key"][0]) @ normal), 6) for s in ident["Segments"]]
        spanning = 0
        for f in ident["Features"]:
            planes = {offsets[int(p[0])] for p in f.get("Portions", [])}
            spanning += len(planes) > 1
        checks["CrossPlaneFeatures"] = {"Expected": 0, "Manifest": spanning, "Pass": spanning == 0}
    return checks


# ----------------------------------------------------------------------- specification ----


def write_specification(layouts, path):
    lines = []
    for lay in layouts:
        lines.append(f"layout {lay['Name']}")
        lines.append(f"box {lay['HalfX']!r} {lay['HalfY']!r} {lay['Depth']!r} {lay['Height']!r}")
        lines.append(f"size {lay['LcFine']!r} {lay['LcFar']!r}")
        if lay.get("Algorithm"):
            lines.append(f"algorithm {int(lay['Algorithm'])}")
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
        for port in lay.get("Ports", []):
            lines.append(f"port {port['Attribute']}")
            for x, y in port["Loops"][0]["Points"]:
                lines.append(f"v {x!r} {y!r}")
            lines.append("close")
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


def design_chains(all_edges, corner_points):
    """Chain root of every polygon edge (edges joined by a sub-threshold vertex, i.e. not a
    classifier corner, are one chain) and the vertices shared by edge pairs."""
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
    return [find(i) for i in range(n)], shared_vertices


def _polyline_distance(points, edges):
    """Closest-point distance from every point (n x 2) to a set of straight edges."""
    result = np.full(len(points), np.inf)
    for e in edges:
        d = e["End"] - e["Start"]
        t = np.clip(((points - e["Start"]) @ d) / float(d @ d), 0.0, 1.0)
        feet = e["Start"][None, :] + t[:, None] * d[None, :]
        result = np.minimum(result, np.linalg.norm(points - feet, axis=1))
    return result


def _chain_order(all_edges, members):
    """Edges of one chain ordered along the chain with their orientation, start position and
    whether the chain bends at either of the edge's joints: [(edge index, forward, x_start,
    bent)], walking from an end (or anywhere on a closed chain)."""
    key = lambda p: (round(float(p[0]), 9), round(float(p[1]), 9))
    incident = {}
    for i in members:
        for p in (all_edges[i]["Start"], all_edges[i]["End"]):
            incident.setdefault(key(p), []).append(i)
    ends = [k for k, v in incident.items() if len(v) == 1]
    start_vertex = ends[0] if ends else key(all_edges[members[0]]["Start"])
    order = []
    used = set()
    vertex = start_vertex
    x = 0.0
    while True:
        candidates = [i for i in incident.get(vertex, []) if i not in used]
        if not candidates:
            break
        i = candidates[0]
        used.add(i)
        e = all_edges[i]
        forward = key(e["Start"]) == vertex
        order.append([i, forward, x, False])
        x += float(e["Length"])
        vertex = key(e["End"] if forward else e["Start"])
        if vertex == start_vertex:
            break
    # A joint bends when the consecutive tangents are not collinear (the windowed curvature is
    # then nonzero on both adjacent edges); every edge records the larger turn at its joints.
    turns = [0.0] * len(order)
    for k in range(len(order) - (0 if order and vertex == start_vertex and len(order) > 1 else 1)):
        i, j = order[k][0], order[(k + 1) % len(order)][0]
        cosine = float(all_edges[i]["Tangent"] @ all_edges[j]["Tangent"])
        if abs(cosine) < 1.0 - 1.0e-9:
            turn = math.acos(max(-1.0, min(1.0, abs(cosine))))
            order[k][3] = True
            order[(k + 1) % len(order)][3] = True
            turns[k] = max(turns[k], turn)
            turns[(k + 1) % len(order)] = max(turns[(k + 1) % len(order)], turn)
    return [(*o, turns[k]) for k, o in enumerate(order)]


def _closest_on_chain(points, all_edges, order):
    """Closest-point distance from every point to an ordered chain, with the chain position of
    the foot and the length of the edge carrying it."""
    best = np.full(len(points), np.inf)
    foot_x = np.zeros(len(points))
    foot_chord = np.zeros(len(points))
    foot_turn = np.zeros(len(points))
    for i, forward, x0, bent, turn in order:
        e = all_edges[i]
        d = e["End"] - e["Start"]
        t = np.clip(((points - e["Start"]) @ d) / float(d @ d), 0.0, 1.0)
        feet = e["Start"][None, :] + t[:, None] * d[None, :]
        dist = np.linalg.norm(points - feet, axis=1)
        better = dist < best
        best[better] = dist[better]
        along = t if forward else 1.0 - t
        foot_x[better] = x0 + along[better] * float(e["Length"])
        foot_chord[better] = float(e["Length"]) if bent else 0.0  # window half-width source: chord where the chain bends, else R
        foot_turn[better] = turn
    return best, foot_x, foot_chord, foot_turn


def design_bent_pairs(all_edges, corner_points, radius, samples=None):
    """Pairs along bends (SURFACE-RESPONSE-IDENTIFICATION.md (b) 7, phase 3): two chains that
    are not both single straight edges pair where their closest-point separation is LOCALLY
    constant. Every chain is sampled (at most R / 2 apart, at least 17 per edge) over its
    candidate facing region (within the reach 2R (1 + tolerance) of the other chain, not
    beyond either of its ends, outside the 2R zones of shared vertices); a sample is constant
    when the distances within R of it along its own chain vary by at most the tolerance; its
    curve separation is the smaller of the two directional maxima over windows of half-width
    max(R, local chord) where the chain bends and R on straight edges, about the sample and
    about its foot (the chord reading C; the inscribed reading is C / cos(turn / 2) with the
    larger local joint turn); it interacts iff both readings are within 2R, i.e. C / cos(turn / 2)
    is within 2R (the straight-pair answer). Constant portions claim their chord-level
    interactions (no event cores); the rest (tees, port ends, fast tapers, acute arms) keep
    the event rule. Returns {(rootA, rootB): stats with the constant sample points}."""
    roots, shared_vertices = design_chains(all_edges, corner_points)
    chains = {}
    for i, r in enumerate(roots):
        chains.setdefault(r, []).append(i)
    interaction = 2.0 * radius - 0.5 * 1.0e-8 * radius
    reach = 2.0 * radius * (1.0 + PAIR_SEPARATION_TOLERANCE) - 0.5 * 1.0e-8 * radius

    def chain_ends(members):
        # Vertices used by exactly one edge of the chain, with the outward tangent.
        ends = []
        for i in members:
            e = all_edges[i]
            for point, outward in ((e["Start"], -e["Tangent"]), (e["End"], e["Tangent"])):
                others = [j for j in members if j != i and any(np.allclose(point, q) for q in (all_edges[j]["Start"], all_edges[j]["End"]))]
                if not others:
                    ends.append((point, outward))
        return ends

    def facing_samples(order_a, order_b, ends_b, shared):
        # Samples of chain A within the reach of chain B: (point, x, own chord, distance, foot
        # x, foot chord, larger joint turn of the own and the foot edge).
        points, xs, chords, turns = [], [], [], []
        for i, forward, x0, bent, turn in order_a:
            e = all_edges[i]
            n = max(16, int(math.ceil(2.0 * float(e["Length"]) / radius)))
            ts = np.linspace(0.0, 1.0, n + 1)
            pa = e["Start"][None, :] + ts[:, None] * (e["End"] - e["Start"])[None, :]
            along = ts if forward else 1.0 - ts
            points.append(pa)
            xs.append(x0 + along * float(e["Length"]))
            chords.append(np.full(len(ts), float(e["Length"]) if bent else 0.0))
            turns.append(np.full(len(ts), turn))
        points = np.concatenate(points)
        xs = np.concatenate(xs)
        chords = np.concatenate(chords)
        turns = np.concatenate(turns)
        d, foot_x, foot_chord, foot_turn = _closest_on_chain(points, all_edges, order_b)
        keep = d < reach
        for point, outward in ends_b:
            keep &= ((points - point) @ outward) <= 0.0
        for v in shared:
            keep &= np.linalg.norm(points - v, axis=1) >= interaction
        return points[keep], xs[keep], chords[keep], d[keep], foot_x[keep], foot_chord[keep], np.maximum(turns, foot_turn)[keep]

    def window_max(xs, d, centre, half):
        sel = (xs >= centre - half) & (xs <= centre + half)
        return (float(d[sel].max()), float(d[sel].min())) if sel.any() else (0.0, np.inf)

    result = {}
    rs = sorted(chains)
    for ia, ra in enumerate(rs):
        for rb in rs[ia + 1 :]:
            ma, mb = chains[ra], chains[rb]
            if len(ma) == 1 and len(mb) == 1 and abs(float(all_edges[ma[0]]["Tangent"] @ all_edges[mb[0]]["Tangent"])) >= 1.0 - 1.0e-8:
                continue  # two exactly parallel straight edges: the translational rule
            shared = [v for (i, j), vs in shared_vertices.items() if roots[i] != roots[j] and {roots[i], roots[j]} == {ra, rb} for v in vs]
            order_a, order_b = _chain_order(all_edges, ma), _chain_order(all_edges, mb)
            sa = facing_samples(order_a, order_b, chain_ends(mb), shared)
            sb = facing_samples(order_b, order_a, chain_ends(ma), shared)
            if len(sa[0]) == 0 or len(sb[0]) == 0:
                continue
            constant_points, interacting_points, separations = [], [], []
            for own, other in ((sa, sb), (sb, sa)):
                points, xs, chords, d, foot_x, foot_chord, turns = own
                for k in range(len(points)):
                    hi, lo = window_max(xs, d, xs[k], radius)
                    constant = hi - lo <= PAIR_SEPARATION_TOLERANCE * lo * (1.0 + 1.0e-9)
                    if not constant:
                        continue
                    w_own, _ = window_max(xs, d, xs[k], max(radius, chords[k]))
                    w_other, _ = window_max(other[1], other[3], foot_x[k], max(radius, foot_chord[k]))
                    # Chord reading C and inscribed reading C / cos(turn / 2): the pair interacts
                    # only when both are within 2R (the recorded discretisation ambiguity).
                    w = min(w_own, w_other) if w_other > 0.0 else w_own
                    constant_points.append(points[k])
                    if within_interaction(w / math.cos(0.5 * turns[k]), radius):
                        interacting_points.append(points[k])
                        separations.append(w)
            if not constant_points:
                continue
            result[(ra, rb)] = {
                "ConstantSamples": len(constant_points),
                "InteractingSamples": len(interacting_points),
                "Samples": len(sa[0]) + len(sb[0]),
                "Interacting": bool(interacting_points),
                "Separation": float(np.mean(separations)) if separations else None,
                "MinSeparation": float(min(separations)) if separations else None,
                "MaxSeparation": float(max(separations)) if separations else None,
                "EdgesA": len(ma),
                "EdgesB": len(mb),
                "ConstantPoints": np.array(constant_points),
            }
    return result


def design_event_cores(all_edges, corner_points, radius, samples=200, bent_pairs=None):
    """Event cores of SURFACE-RESPONSE-IDENTIFICATION.md (b) 3 from the polygon edges: sampled
    points on an edge that are within 2R of a point on a non-parallel edge of another chain,
    excluding through-vertex pairs (either point within 2R of a vertex the two chains share)
    and pairs of chains that pair along a bend (described by the pair feature, (b) 7).
    Edges joined by a sub-threshold vertex (not a classifier corner) are one chain."""
    n = len(all_edges)
    roots, shared_vertices = design_chains(all_edges, corner_points)

    def find(i):
        return roots[i]

    bent = bent_pairs or {}
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
            # Points inside a locally constant portion of the pair along a bend (within the
            # sample spacing R / 2 of a constant sample of either chain) are not events.
            record = bent.get((min(find(i), find(j)), max(find(i), find(j))))
            if record is not None and len(record["ConstantPoints"]):
                constant = record["ConstantPoints"]
                paired_a = distance_to_cores(pa, constant) <= 0.5 * radius
                paired_b = distance_to_cores(pb, constant) <= 0.5 * radius
                close &= ~(paired_a[:, None] | paired_b[None, :])
            for v in shared_vertices.get((min(i, j), max(i, j)), []):
                za = np.linalg.norm(pa - v, axis=1) < interaction
                zb = np.linalg.norm(pb - v, axis=1) < interaction
                close &= ~(za[:, None] | zb[None, :])  # either point inside the zone
            hits = pa[close.any(axis=1)]
            if len(hits):
                cores.append(hits)
    # Two corners closer than 2R (overlapping windows) are an event of their own.
    for i, p in enumerate(corner_points):
        for q in corner_points[i + 1 :]:
            if np.linalg.norm(p - q) < interaction:
                cores.append(np.array([p, q]))
    if not cores:
        return np.zeros((0, 2))
    # Decimate the sampled cores to a grid of R / 200 (polyline arcs produce O(10^5) hits).
    points = np.concatenate(cores)
    quantum = radius / 200.0
    return np.unique(np.round(points / quantum), axis=0) * quantum


def distance_to_cores(points, cores, chunk=20000):
    """Minimum distance from every point to the sampled cores (chunked: O(points x cores))."""
    result = np.full(len(points), np.inf)
    for start in range(0, len(cores), chunk):
        block = cores[start : start + chunk]
        result = np.minimum(result, np.linalg.norm(points[:, None, :] - block[None, :, :], axis=2).min(axis=1))
    return result


def off_plane_obstacles(lay):
    """Metal off the process plane: facing sheets (polygon at height Z) and walls (their foot
    segment on the plane; the wall stands on the plane so the nearest wall point to a plane
    point lies on the foot)."""
    obstacles = []
    for sh in lay["Sheets"]:
        if sh["Z"] != 0.0:
            for lp in sh["Loops"][:1]:
                obstacles.append(("sheet", float(sh["Z"]), [np.array(p, dtype=float) for p in lp["Points"]]))
    for wall in lay["Walls"]:
        _, x0, y0, x1, y1, _ = wall
        obstacles.append(("wall", 0.0, [np.array([x0, y0], dtype=float), np.array([x1, y1], dtype=float)]))
    return obstacles


def _point_polygon_distance_2d(point, polygon):
    inside = True
    n = len(polygon)
    for i in range(n):
        a, b, c = polygon[i], polygon[(i + 1) % n], polygon[(i + 2) % n]
        d = b - a
        normal = np.array([-d[1], d[0]])
        if float((c - b) @ normal) * float((point - a) @ normal) < 0.0:
            inside = False
            break
    if inside:
        return 0.0
    best = math.inf
    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        d = b - a
        t = min(1.0, max(0.0, float((point - a) @ d) / float(d @ d)))
        best = min(best, float(np.linalg.norm(point - (a + t * d))))
    return best


def obstacle_distance(point, obstacles):
    """Distance from a plane point to the nearest off-plane metal (facing sheets, walls)."""
    best = math.inf
    for kind, z, geometry in obstacles:
        if kind == "sheet":
            planar = _point_polygon_distance_2d(point, geometry)
            best = min(best, math.hypot(planar, z))
        else:
            a, b = geometry
            d = b - a
            t = min(1.0, max(0.0, float((point - a) @ d) / float(d @ d)))
            best = min(best, float(np.linalg.norm(point - (a + t * d))))
    return best


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
    for sheet_index, sh in enumerate(lay["Sheets"]):
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
                # Conductor identity is metal connectivity (phase 3): every sheet is its own
                # conductor (the layouts' sheets are disjoint polygons).
                e["Conductor"] = sheet_index
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
            edge_indices = (i, j)
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
                    "EdgeIndices": edge_indices,
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
    # Two corners are an event of their own when within 2R (the classifier's quantized strict
    # decision: exactly 2R is not within).
    corner_pairs_within_2r = sum(1 for i in range(len(corner_points)) for j in range(i + 1, len(corner_points)) if within_interaction(float(np.linalg.norm(corner_points[i] - corner_points[j])), radius))
    # Design rule (SURFACE-RESPONSE-IDENTIFICATION.md (b) 7): chains that pair along a bend
    # are pair features; their cross-chord interactions are not events. The expected
    # curvature class follows the layout's design bend (inner side radius vs 10 R, decision 75).
    bent_pairs = design_bent_pairs(all_edges, corner_points, radius)
    roots, _ = design_chains(all_edges, corner_points)
    bend = lay.get("Bend")
    bent_pair_records = []
    for (ra, rb), stats in sorted(bent_pairs.items()):
        record = {k: v for k, v in stats.items() if k != "ConstantPoints"}
        record["Chains"] = [int(ra), int(rb)]
        if not record["Interacting"]:
            record["ExpectedClasses"] = []
        else:
            # The pair class: two sheets -> different-conductor gap; one sheet -> the design's
            # bar (its two sides are a strip) or a same-conductor gap.
            a = all_edges[next(i for i in range(len(all_edges)) if roots[i] == ra)]
            b = all_edges[next(i for i in range(len(all_edges)) if roots[i] == rb)]
            if a["Conductor"] != b["Conductor"]:
                base = "DifferentConductorGap"
            else:
                base = "SameConductorStrip" if bend and "Gap" not in bend else "SameConductorGap"
            curved = False
            if bend and "Radius" in bend:
                # A bar: its two sides are a strip of the design width (inner side radius =
                # radius - width / 2); a gap between bars: the facing sides at radius -/+ gap / 2.
                inner = bend["Radius"] - 0.5 * bend["Width"] if "Gap" not in bend else bend["Radius"] - 0.5 * bend["Gap"]
                record["DesignInnerRadiusOverR"] = inner / radius
                curved = inner / radius < STRAIGHT_BEND_RADIUS_OVER_R
            record["Curved"] = curved
            record["ExpectedClasses"] = sorted({base} | ({"Curved" + base} if curved else set()))
            record["ExpectedSeparation"] = (bend or {}).get("Gap", (bend or {}).get("Width", record["Separation"]))
        bent_pair_records.append(record)
    # Design rule (SURFACE-RESPONSE-IDENTIFICATION.md (b) 3): a corner joins a cluster when
    # its 2R through-vertex zone reaches a cluster region, i.e. an event core lies within 3R.
    cores = design_event_cores(all_edges, corner_points, radius, bent_pairs=bent_pairs)
    for c in corners:
        if not c["ClassifierCorner"]:
            c["Standalone"] = None
            continue
        point = np.array(c["Point"])
        core_distance = float(distance_to_cores(point[None, :], cores)[0]) if len(cores) else math.inf
        c["CoreDistance"] = None if math.isinf(core_distance) else round(core_distance, 6)
        c["Standalone"] = not core_distance < 3.0 * radius
        if not c["Standalone"]:
            c["Expected"] += " (an interaction event core within 3R: member of a spatial cluster)"
    # Decision 73(3): metal off the process plane (a facing sheet, a wall) within 2R excludes
    # the corners and the edge portions it reaches (Identifier::ClassifyPlanes CrossLayer).
    obstacles = off_plane_obstacles(lay)
    excluded_corners = 0
    for c in corners:
        if not c["ClassifierCorner"]:
            continue
        distance = obstacle_distance(np.array(c["Point"]), obstacles)
        c["OffPlaneDistance"] = None if math.isinf(distance) else round(distance, 9)
        c["ExcludedCrossLayer"] = distance < 2.0 * radius * (1.0 - 1.0e-9)
        if c["ExcludedCrossLayer"]:
            c["Standalone"] = False
            c["Expected"] += " (metal off the plane within 2R: excluded vertex)"
            excluded_corners += 1
    cross_layer_length = 0.0
    if obstacles:
        for e in all_edges:
            length = float(e["Length"])
            interval = P._sublevel_interval(lambda t: obstacle_distance(e["Start"] + t * (e["End"] - e["Start"]), obstacles), 0.0, 1.0, 2.0 * radius)
            if interval is not None:
                cross_layer_length += (interval[1] - interval[0]) * length
    excluded["CrossLayerLength"] = cross_layer_length
    excluded["ExcludedCorners"] = excluded_corners
    # Corners of the metal off the process plane (facing sheets, wall tops): geometric corners
    # of the perimeter census (excluded from the identification, never manifest records).
    off_plane_corners = []
    for sh in lay["Sheets"]:
        if sh["Z"] == 0.0:
            continue
        for hole_index, lp in enumerate(sh["Loops"]):
            hole = hole_index > 0
            points = lp["Points"]
            for index, point in enumerate(points):
                if _on_box(point, lay):
                    continue
                t_in = _tangent_at_vertex(lp, index, True)
                t_out = _tangent_at_vertex(lp, index, False)
                turn = math.degrees(math.acos(float(np.clip(t_in @ t_out, -1.0, 1.0))))
                if turn <= corner_turn_tolerance + 1.0e-9:
                    continue
                convex = (cross2(t_in, t_out) > 0) != hole
                off_plane_corners.append({"Point": [round(float(point[0]), 9), round(float(point[1]), 9), float(sh["Z"])], "TurnDegrees": round(turn, 9), "InteriorAngleDegrees": round(180.0 - turn if convex else 180.0 + turn, 9), "Convex": bool(convex), "Expected": "Excluded (facing sheet)"})
    for _, x0, y0, x1, y1, height in lay["Walls"]:
        for x, y in ((x0, y0), (x1, y1)):
            off_plane_corners.append({"Point": [float(x), float(y), float(height)], "TurnDegrees": 90.0, "InteriorAngleDegrees": 90.0, "Convex": True, "Expected": "Excluded (wall top)"})
    # A parallel pair is a feature when part of its overlap survives the cluster regions
    # (distance to a core >= R) and the corner windows (R along the edge from a corner).
    for pair in pairs:
        i, j = pair.pop("EdgeIndices")
        pair["InBentPair"] = (min(roots[i], roots[j]), max(roots[i], roots[j])) in bent_pairs
        if not pair["Within2R"]:
            pair["Survives"] = False
            continue
        a = pair["EdgeA"]
        ts = np.linspace(pair["OverlapStart"], pair["OverlapEnd"], 201)
        points = a["Start"][None, :] + (ts[:, None] - float(a["Start"] @ a["Tangent"])) * a["Tangent"][None, :]
        survives = np.ones(len(ts), dtype=bool)
        if len(cores):
            survives &= distance_to_cores(points, cores) >= radius
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
        "OffPlaneCorners": off_plane_corners,
        "ClassifierCornerCount": sum(1 for c in corners if c["ClassifierCorner"]),
        "StandaloneCornerCount": sum(1 for c in corners if c.get("Standalone")),
        "ExpectedRoundedCorners": sum(1 for a in arcs if a["ExpectedRoundedCorner"]),
        "SubThresholdTurns": sum(1 for c in corners if not c["ClassifierCorner"]),
        "CornerPairsWithin2R": corner_pairs_within_2r,
        "ParallelPairs": pairs,
        "BentPairs": bent_pair_records,
        "NonparallelInteractions": nonparallel,
        "Arcs": arcs,
        "Excluded": excluded,
        "Bend": bend,
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
    # Mesh corners are geometric (turn > 30 deg) whether or not the vertex is excluded: the
    # in-plane corners (z = 0) and the corners of the metal off the plane, matched in 3D.
    census_corners = [dict(c, Point=[c["Point"][0], c["Point"][1], 0.0]) for c in oracle_corners] + list(orc.get("OffPlaneCorners", []))
    mesh_corners = [c for c in census["Corners"] if c["Kind"] == "CORNER"]
    matched = 0
    angle_mismatch = []
    unmatched_mesh = []
    for mc in mesh_corners:
        point = np.array(mc["Point"][:3], dtype=float)
        hit = next((oc for oc in census_corners if np.linalg.norm(np.array(oc["Point"], dtype=float) - point) < 1.0e-6), None)
        if hit is None:
            unmatched_mesh.append(mc)
        elif abs(hit["InteriorAngleDegrees"] - mc["InteriorAngleDegrees"]) > 1.0e-6 and abs(hit["TurnDegrees"] - mc["TurnDegrees"]) > 1.0e-6:
            angle_mismatch.append({"Oracle": hit, "Mesh": mc})
        else:
            matched += 1
    checks["A6-corner-list-mesh"] = {
        "OracleCorners": len(oracle_corners),
        "OracleOffPlaneCorners": len(census_corners) - len(oracle_corners),
        "MeshCorners": len(mesh_corners),
        "Matched": matched,
        "AngleMismatch": angle_mismatch,
        "UnmatchedMesh": unmatched_mesh[:20],
        "Pass": matched == len(census_corners) == len(mesh_corners) and not angle_mismatch,
        "Meaning": "mesh corners beyond the oracle's (in-plane + off-plane metal) are arc discretisation vertices above the 30 deg tolerance",
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
    oracle_pairs = Counter((p["Expected"], round(p["Separation"], 6)) for p in orc["ParallelPairs"] if p["Within2R"] and p.get("Survives", True) and not p.get("InBentPair"))
    mesh_separations = census["Interactions"]["ParallelSeparations"]
    manifest_pairs = Counter()
    for r in manifest["Requirements"]:
        if r["Topology"] in M.TRANSLATIONAL and r["Topology"] not in ("IsolatedEdge", "CurvedEdge"):
            manifest_pairs[(r["Topology"], round(float(r["Geometry"].get("Separation", r["Geometry"].get("Width", float("nan")))), 6))] += 1
    expected_pairs_set = {(k[0], k[1]) for k in oracle_pairs}
    manifest_pairs_set = set(manifest_pairs)
    # Records whose constant portion is significant (more than 2R of chain: > 4 samples at
    # <= R / 2 spacing) or interacting; a perpendicular edge end against a side gives a couple
    # of trivially constant samples at its piece end.
    bent_records = [r for r in orc.get("BentPairs", []) if r["Interacting"] or r["ConstantSamples"] > 4]
    if bent_records:
        # Pairs along bends: the manifest's pair classes must be exactly the expected classes
        # (straight-like strip, plus the curved strip when the design inner radius is below
        # STRAIGHT_BEND_RADIUS_OVER_R x R) with separations within the pair tolerance of the design width; the chord-wise
        # parallel pairs inside a bent pair are not separate features.
        expected_classes = set()
        for record in bent_records:
            expected_classes |= set(record.get("ExpectedClasses", []))
        manifest_classes = {k[0] for k in manifest_pairs_set}
        widths = [record["ExpectedSeparation"] for record in bent_records if "ExpectedSeparation" in record]
        separations_ok = all(any(abs(k[1] - w) <= PAIR_SEPARATION_TOLERANCE * w for w in widths) for k in manifest_pairs_set) if widths else True
        pass_pairs = (manifest_classes == (expected_classes | {k[0] for k in expected_pairs_set})) and separations_ok and all(k in manifest_pairs_set for k in expected_pairs_set)
    else:
        pass_pairs = expected_pairs_set == manifest_pairs_set
    checks["A6-parallel-pairs"] = {
        "OraclePairs": {f"{k[0]}@{k[1]}": v for k, v in sorted(oracle_pairs.items())},
        "OracleBentPairs": bent_records,
        "OracleAtExactly2R": sum(1 for p in orc["ParallelPairs"] if p["AtExactly2R"]),
        "OracleAtExactlyR": sum(1 for p in orc["ParallelPairs"] if p["AtExactlyR"]),
        "MeshParallelSeparations": mesh_separations,
        "ManifestPairClasses": {f"{k[0]}@{k[1]}": v for k, v in sorted(manifest_pairs.items())},
        "ManifestTopologies": dict(manifest_topologies),
        "Pass": pass_pairs,
        "Meaning": "expected class and separation of every parallel pair within 2R vs the manifest's translational records (pairs along bends: expected classes with separations within the pair tolerance); pairs at exactly 2R are the knife edge",
    }
    if bent_records:
        interacting = [r for r in bent_records if r["Interacting"]]
        bar_layout = bool(orc.get("Bend")) and "Width" in orc["Bend"] and "Gap" not in orc["Bend"]
        if interacting and bar_layout:
            # A pair along a bend leaves no isolated or curved-edge remainder on its chains and
            # the bar ends are the only clusters (one per group of corners within 2R).
            checks["A6-bent-pair-classes"] = {
                "ManifestTopologies": dict(manifest_topologies),
                "ExpectedClusters": orc["CornerPairsWithin2R"],
                "Pass": manifest_topologies.get("IsolatedEdge", 0) == 0 and manifest_topologies.get("CurvedEdge", 0) == 0 and manifest_topologies.get("SpatialEdgeCluster", 0) == orc["CornerPairsWithin2R"],
                "Meaning": "the two sides of a constant-width bend are fully paired; the strip ends (two corners within 2R) are the only spatial clusters",
            }
        else:
            # A gap between two bars along a bend, decided on the curve separation like a
            # straight pair: within 2R -> one gap class with the design separation and the
            # corner pairs across the gap as the only clusters; at or beyond 2R -> no pair
            # class, no cluster, the sides are isolated edges (no curved edge: inner radius >=
            # 10 R). The chord dips below 2R of the polylines must not change this.
            separations = sorted({round(k[1], 6) for k in manifest_pairs_set})
            expected_separation = sorted({round(r["ExpectedSeparation"], 6) for r in interacting if "ExpectedSeparation" in r})
            expected_clusters = (orc.get("Bend") or {}).get("Clusters", orc["CornerPairsWithin2R"])
            separations_match = len(separations) == len(expected_separation) and all(abs(m - e) <= PAIR_SEPARATION_TOLERANCE * e for m, e in zip(separations, expected_separation))
            checks["A6-bent-pair-classes"] = {
                "ManifestTopologies": dict(manifest_topologies),
                "ManifestPairSeparations": separations,
                "ExpectedSeparations": expected_separation,
                "ExpectedClusters": expected_clusters,
                "Interacting": bool(interacting),
                "Pass": manifest_topologies.get("SpatialEdgeCluster", 0) == expected_clusters
                and manifest_topologies.get("CurvedEdge", 0) == 0
                and (separations_match if interacting else (not separations and manifest_topologies.get("IsolatedEdge", 0) > 0)),
                "Meaning": "pairs with divergent ends / tapers / the 2R threshold: the straight-pair answer at every discretisation (within 2R: the pair classes at the design separations within the pair tolerance; at or beyond 2R: isolated edges, no pair); clusters only at the corner pairs / tees",
            }
    checks["A6-nonparallel-interactions"] = {"Oracle": len(orc["NonparallelInteractions"]), "Mesh": census["Interactions"]["NonparallelPairs"], "Pass": None, "Meaning": "record only: non-parallel pairs are omitted by the classifier"}
    # Decision 73(3) exclusions: the manifest's CrossLayer record must carry the analytic
    # length of the perimeter within 2R of the off-plane metal, and its excluded vertices the
    # corners within 2R; walls give NonPlanar / NonManifold records.
    recorded = Counter()
    identification = manifest.get("Identification") or {}
    for e in identification.get("Exclusions", []):
        recorded[e["Class"]] += float(e["Length"])
    manifest_excluded_vertices = sum(1 for v in identification.get("Vertices", []) if v["Type"] == "Excluded")
    expected_cross_layer = orc["Excluded"]["CrossLayerLength"]
    checks["A6-excluded-classes"] = {
        "Oracle": orc["Excluded"],
        "MeshLengthByClass": {k: census["LengthByClass"].get(k, 0.0) for k in ("NONPLANAR", "FOLD", "CROSS_LAYER", "NONMANIFOLD", "EMBEDDED")},
        "ManifestRecorded": dict(recorded),
        "ManifestExcludedVertices": manifest_excluded_vertices,
        "Pass": abs(recorded.get("CrossLayer", 0.0) - expected_cross_layer) <= 1.0e-6 * max(1.0, expected_cross_layer)
        and (orc["Excluded"]["Walls"] > 0) == (recorded.get("NonPlanar", 0.0) > 0)
        and (orc["Excluded"]["Walls"] > 0) == (recorded.get("NonManifold", 0.0) > 0)
        and (orc["Excluded"]["CrossLayerSheets"] > 0) == (recorded.get("UndeterminedProcessSide", 0.0) > 0)
        and manifest_excluded_vertices >= orc["Excluded"]["ExcludedCorners"],
        "Meaning": "CrossLayer length = analytic perimeter within 2R of facing sheets / walls; walls -> NonPlanar + NonManifold; a facing sheet's own edges -> CrossLayer or UndeterminedProcessSide",
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
        crack_digests = {}
        cracks = [True] + ([False] if args.crack_false else [])
        libraries = dict(args.libraries)
        if args.signature_library:
            libraries["signature"] = os.path.join(args.output, lay["Name"], "signature-library.json")
        for label, library in libraries.items():
            if label == "signature" and args.signature_library:
                first_label = next(iter(args.libraries))
                source_path = os.path.join(args.output, lay["Name"], f"{first_label}-u0-np{args.ranks[0]}", "postpro", "surface-response-requirements.json")
                if not os.path.exists(source_path):
                    record["Cells"].append({"Name": "signature", "Library": "signature", "Error": "no manifest to build the signature library from"})
                    continue
                with open(source_path) as source:
                    signature_library = build_signature_library(json.load(source), name="signature-only")
                with open(library, "w") as target:
                    json.dump(signature_library, target, indent=1)
            for levels in args.uniform_levels:
              for crack in cracks:
                for ranks in args.ranks:
                    if not crack and (levels != 0 or ranks != args.ranks[0]):
                        continue
                    name = f"{label}-u{levels}-np{ranks}" if crack else f"{label}-u{levels}-crackfalse-np{ranks}"
                    directory = os.path.join(args.output, lay["Name"], name)
                    os.makedirs(directory, exist_ok=True)
                    config = preflight_config(mesh_path, ground, terminals, sa, library, os.path.join(directory, "postpro"), uniform_levels=levels, crack=crack, ports=[port["Attribute"] for port in lay.get("Ports", [])], frame_normal=lay.get("FrameNormal"))
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
                        if crack:
                            digests.setdefault((label, levels), {})[ranks] = cell["GeometryDigest"] or cell["DigestFull"]
                        else:
                            crack_digests[label] = cell["GeometryDigest"] or cell["DigestFull"]
                        if levels == 0:
                            audit_args = argparse.Namespace(mesh=mesh_path, config=config_path, manifest=manifest_path, log=log_path, compare=None, radius=None, corner_tolerance=CORNER_TURN_TOLERANCE_DEGREES, output_prefix=None)
                            result = audit.run_audit(audit_args)
                            with open(os.path.join(directory, "audit.json"), "w") as target:
                                json.dump(result, target, indent=1, default=str)
                            with open(os.path.join(directory, "audit.md"), "w") as target:
                                target.write(audit.render_markdown(result))
                            cell["Gates"] = {g["Gate"]: g["Status"] for g in result["Gates"]}
                            cell["Patches"] = {k: v for k, v in (result.get("Patches") or {}).items() if k != "Models"}
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
        # Crack independence: CrackInternalBoundaryElements false reproduces the digest (and
        # the per-signature lengths through the audit gates of its own cell).
        record["CrackIndependence"] = {label: crack_digests[label] == digests.get((label, 0), {}).get(args.ranks[0]) for label in crack_digests} if crack_digests else None
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
            la, lb = {}, {}
            for f in ia["Features"]:
                la.setdefault(f["Hash"], []).append(float(f["Length"]))
            for f in ib["Features"]:
                lb.setdefault(f["Hash"], []).append(float(f["Length"]))
            radius = float(ia["MatchingRadius"])
            mismatch = [h[:12] for h in set(la) | set(lb) if len(la.get(h, [])) != len(lb.get(h, [])) or any(abs(x - y) > 1.0e-6 * max(abs(x), radius) for x, y in zip(sorted(la.get(h, [])), sorted(lb.get(h, []))))]
            return {"GeometryDigestIdentical": ia["GeometryDigest"] == ib["GeometryDigest"] and not mismatch, "OnlyA": sorted(f"{t}:{h}x{n}" for (t, h), n in (ca - cb).items()), "OnlyB": sorted(f"{t}:{h}x{n}" for (t, h), n in (cb - ca).items()), "LengthMismatch": mismatch}

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
        "Crack": record.get("CrackIndependence"),
        "A5": {k: (v["Identification"]["GeometryDigestIdentical"] if v.get("Identification") else v["GeometryOnlyIdentical"]) for k, v in (record.get("A5-refinement-invariance") or {}).items()},
        "A3": {k: (v["Identification"]["GeometryDigestIdentical"] if v.get("Identification") else v["GeometryOnlyIdentical"]) for k, v in (record.get("A3-library-independence") or {}).items()},
        "Oracle": {label: {k.replace("A6-", ""): v["Pass"] for k, v in c["OracleChecks"].items() if v["Pass"] is not None} for label, c in by_library.items()},
        "GatesFailing": {label: [k for k, v in c["Gates"].items() if v != "PASS"] for label, c in by_library.items()},
        "PatchCoverage": {label: round(c["Patches"]["CoveredFractionOfAssigned"], 9) for label, c in by_library.items() if c.get("Patches") and c["Patches"].get("CoveredFractionOfAssigned") is not None},
        "Topologies": {label: {k: v["Count"] for k, v in c["ManifestSummary"].items()} for label, c in by_library.items()},
        "Error": record.get("Error"),
    }


def write_summary(results, output):
    lines = ["# Synthetic stress layouts: identification vs oracle", ""]
    lines.append("| Layout | nodes | exit | A4 ranks | crack | A5 refine | A3 libraries | A6 oracle (per library) | audit gates failing | patch coverage | manifest topologies |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for r in results:
        row = summary_row(r)
        oracle_text = "; ".join(f"{label}: " + ", ".join(f"{k}={'P' if v else 'F'}" for k, v in checks.items()) for label, checks in row["Oracle"].items())
        gates_text = "; ".join(f"{label}: {', '.join(g.replace('A1-', '').replace('A2-', '') for g in gates)}" for label, gates in row["GatesFailing"].items())
        topology_text = "; ".join(f"{label}: " + ", ".join(f"{k}:{v}" for k, v in t.items()) for label, t in row["Topologies"].items())
        coverage_text = "; ".join(f"{label}: {v}" for label, v in row["PatchCoverage"].items())
        lines.append(f"| {row['Layout']} | {row['Nodes']} | {','.join(row['Exit'])} | {row['A4']} | {row['Crack']} | {row['A5']} | {row['A3']} | {oracle_text} | {gates_text} | {coverage_text} | {topology_text} |")
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
    parser.add_argument("--crack-false", action="store_true", help="add a CrackInternalBoundaryElements=false cell (level 0, first rank count) per library and compare its digest")
    parser.add_argument("--signature-library", action="store_true", help="per layout, add a 'signature' library built from the first library's level-0 manifest (signature-only models for every feature: the fully matched patch dry run must cover the whole assigned perimeter)")
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
