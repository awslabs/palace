#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Window inventory of a version-2 identification manifest (validation plan stage 0, E0).

    python3 -m surface_response_identification.window_query \\
        --manifest postpro/surface-response-requirements.json \\
        --window S1 500 700 -200 -30 [--window ...] [--margin 19] [--weights weights.json] \\
        --output out.json [--csv out.csv]

For every window (x0, x1, y0, y1 in the manifest's mesh units; every plane z, or ``--z zmin zmax``) the
tool clips every segment of the segment table to the box and reports, from the clipped portions only:

* the perimeter P inside the box, split into assigned / excluded (per exclusion class) and per plane (z);
* every feature with a portion inside: Type, signature Hash, plane, length inside, total length, whether the
  feature also has length outside the box (a cut-through feature) and whether any of its portions lies within
  ``--margin`` of a wall (the plan's exclusion margin, default 10 R);
* the class x key x plane table (count of features, length inside) and the proxy energy per key = length
  inside x per-length weight (``--weights``: {Type: weight}; defaults = the plan's transmon proxies:
  IsolatedEdge 1, CurvedEdge 1, SameConductorStrip / SameConductorGap 2.4, ParallelEdgeCluster 3,
  ConvexCorner 14.6, ConcaveCorner 0.01, SpatialEdgeCluster 8; recorded in the output);
* the metal bodies (conductors) whose perimeter enters the box, from the perimeter loops of the chip: every
  closed loop of PHYSICAL perimeter segments (assigned or excluded as TruncationCut / SimulationBoundary /
  Port / CrossLayer / Untargeted / UndeterminedProcessSide; NonManifold and NonPlanar segments are walls,
  not perimeter) on one plane has its metal side decided by the sharp corner features on it (on a loop with
  metal inside a polygon-convex vertex is a ConvexCorner of the metal, with metal outside a ConcaveCorner;
  majority vote), a loop without sharp corners takes the opposite side of its enclosing loop, and a
  parentless one is a hole of the plane's unbounded ground (a flip chip's L1 ground has no outer boundary in
  the perimeter table). A loop with metal inside is a body (an island, or a plane's outer boundary); a hole
  belongs to the nearest enclosing body, else to the plane's ground (body index < 0). Bodies are reported with the plan's terminal assignment:
  every body is ground (0 V) or a terminal; the proposed quoted terminal is the island entirely inside the
  box with the largest proxy energy inside, else none (the CPW centre trace open-terminated, chosen at cut
  time); a body whose perimeter crosses a wall is flagged ``CutByWall`` (bridged to ground at the wall
  unless it is the terminal); a perimeter component that does not close (interrupted by NonPlanar / NonManifold
  edges: airbridge spans, wirebond feet) is reported as ``OpenChain`` with no body assignment. Minority corner
  votes on a loop are reported as ``MetalSideDisagreements``;
* the class-(B) mean-separation prediction of plan (b)-4 for every pair / stack with a portion inside: the
  mean separation of the window-interior portions of side 0 to the other sides (sampled at 0.5 R) vs the
  chip's signature separation (``SeparationOverR`` / ``OffsetOverR``), and whether the window reading would
  change the signature key (difference beyond the 1e-3 R signature tolerance);
* the cut-arc exclusions of the E1 comparison region (``E1CutArcExclusions``, VALIDATION-PLAN (h)-9, decision 214 (iii)): a chip
  arc (``Arcs``, the segments' ``Arc`` id) is CUT when one of the three window-introduced lines crosses it: a window wall, the
  setback line (the open-terminated terminal's metal is clipped to the setback box) or the edge of a wall BRIDGE strip (the metal
  rectangles ``window_polygons.py`` lays over the gaps next to a bridged conductor's wall interval, 3 R into the box: the chip
  edges inside a strip vanish, so the arm of a cut-adjacent joint ends on the strip's edge, not on the wall). The bridges are
  taken from the window's polygon set (``--polygons``: ``WindowPolygons.Bridges``, the one source of that geometry); without it
  the clip region is the box alone, which is conservative (a bridge edge only SHORTENS a wall-cut arm, so more joints are kept).
  The identification's own concyclicity rule predicts the window's reading of a cut arc: at a cut-adjacent joint the arm is the
  cut chord, whose far end is off the circle, so the joint stays an end joint only when that arm meets the circle's tangent
  within the joint noise rule (CornerRule: ``(c / 2) tan(kink / 4) < JointNoiseSagittaOverR x R`` on the shorter of the arm and
  the first chord, kink = half the chord's central angle); a dropped cut-adjacent joint becomes a corner (``DroppedJoint``: 3 R
  along the chain on both sides of it excluded); the bend survives on the kept joints only when at least four remain (the
  least-squares clause), else every kept joint becomes a corner (``Collapse``: the arc's chords and 3 R along both arms
  excluded); an arc whose joints are all noise-level reads straight in both meshes (``None``). The 3 R walk along the chain
  stops short at a branching vertex, a chain end or after the first chord of another arc (less exclusion: fails closed); every
  such stop is recorded per arc in ``MarginTruncations``. Limit: a CLOSED arc cut by a wall whose inside run turns more than
  180 degrees cannot be one arc in the window (the identification's <= 180-degree clause) while the tool predicts from the kept
  joints alone (``None`` when they are all tangent-noise); such a window reading surfaces as a ``PredictionMismatch`` (none on
  record: S1b 10211 keeps 25 of 161 joints, 56 degrees). The pieces are consumed by ``segment_identity.py --exclude`` (class
  ``Excluded``, reported, never a defect) and the prediction is checked there against the observed classes. ``SetbackHints``
  lists, for the terminal's collapsing arcs, the setback that would cut the arc away whole (information for the stage-2 window
  placement; nothing is moved).

Coordinates are the manifest's (absolute chip coordinates); a segment clipped to the box contributes only
its inside length. The tool reads the manifest once for any number of windows.
"""

import argparse
import collections
import json
import math
import sys

DEFAULT_WEIGHTS = {"IsolatedEdge": 1.0, "CurvedEdge": 1.0, "SameConductorStrip": 2.4, "SameConductorGap": 2.4,
                   "DifferentConductorGap": 2.4, "CurvedSameConductorStrip": 2.4, "CurvedSameConductorGap": 2.4,
                   "CurvedDifferentConductorGap": 2.4, "ParallelEdgeCluster": 3.0, "ConvexCorner": 14.6,
                   "ConcaveCorner": 0.01, "SpatialEdgeCluster": 8.0, "Endpoint": 14.6, "Junction": 14.6}
PERIMETER_EXCLUSIONS = ("TruncationCut", "SimulationBoundary", "Port", "PortCut", "CrossLayer", "Untargeted",
                        "UndeterminedProcessSide")
SIGNATURE_TOLERANCE_OVER_R = 1e-3
POINT_DECIMALS = 6
LEAST_SQUARES_BEND_JOINTS = 4  # ArcRule: a least-squares bend needs four joints (three points are always concyclic)
JOINT_NOISE_SAGITTA_OVER_R_DEFAULT = 0.05  # CornerRule's JointNoiseSagittaOverR when the manifest carries no Conventions
CUT_ARC_MARGIN_OVER_R = 3.0  # (h)-1: 3 R around every window-introduced feature


def clip_interval(p0, p1, box):
    """[t0, t1] in [0, 1] of the parametric segment p0 + t (p1 - p0) inside the box, else None."""
    t0, t1 = 0.0, 1.0
    for axis in range(2):
        lo, hi = box[2 * axis], box[2 * axis + 1]
        d = p1[axis] - p0[axis]
        if abs(d) < 1e-15:
            if p0[axis] < lo or p0[axis] > hi:
                return None
            continue
        ta, tb = (lo - p0[axis]) / d, (hi - p0[axis]) / d
        if ta > tb:
            ta, tb = tb, ta
        t0, t1 = max(t0, ta), min(t1, tb)
        if t0 > t1:
            return None
    return t0, t1


def point_key(p):
    return tuple(round(float(c), POINT_DECIMALS) for c in p[:3])


def _segment_distance(px, py, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    l2 = dx * dx + dy * dy
    if l2 <= 0.0:
        return math.hypot(px - a[0], py - a[1])
    t = max(0.0, min(1.0, ((px - a[0]) * dx + (py - a[1]) * dy) / l2))
    return math.hypot(px - a[0] - t * dx, py - a[1] - t * dy)


def _point_in_polygon(x, y, polygon):
    inside = False
    n = len(polygon)
    j = n - 1
    for i in range(n):
        xi, yi = polygon[i]
        xj, yj = polygon[j]
        if (yi > y) != (yj > y):
            xcross = xj + (y - yj) * (xi - xj) / (yi - yj)
            if x < xcross:
                inside = not inside
        j = i
    return inside


class PerimeterLoops:
    """Closed loops of physical perimeter segments per plane, with containment depth and body assignment."""

    def __init__(self, segments, features, radius):
        self.segments = segments
        self.radius = radius
        self.segment_loop = {}
        self.loops = {}  # index -> dict: Plane, Segments, Polygon, Area, BBox, Depth, Parent, Body, Closed, MetalInside
        self.ground_bodies = {}
        self.metal_side_checks = 0
        self.metal_side_disagreements = []
        self.bodies = None  # window_extract body records (chip conductors) when built from an extract
        self.footprints = []
        if segments is not None:
            self._build(features)

    @classmethod
    def from_extract(cls, manifest):
        """The CHIP's loop / body model restricted to a window extract (``window_extract.py``): the loops keep their
        chip indices (sparse), ``segment_loop`` maps extract segment indices, ``bodies`` holds the extract's body
        records (Kind, Plane, chip Conductor, ConductorBodies, ConductorGround) and ``footprints`` the bump footprints."""
        we = manifest["WindowExtract"]
        model = cls(None, None, float(manifest["Identification"]["MatchingRadius"]))
        model.segments = manifest["Identification"]["Segments"]
        for rec in we["Loops"]:
            loop = dict(rec)
            idx = loop.pop("Index")
            loop["Segments"] = [j for j in rec["Segments"] if j >= 0]
            for j in loop["Segments"]:
                model.segment_loop[j] = idx
            model.loops[idx] = loop
        model.ground_bodies = {("Ground", float(p)): v for p, v in we.get("GroundBodies", {}).items()}
        model.bodies = {int(k): v for k, v in we["Bodies"].items()}
        model.footprints = we["Footprints"]
        return model

    def _build(self, features):
        by_plane = collections.defaultdict(dict)  # plane -> point -> [segment indices]
        physical = []
        for i, s in enumerate(self.segments):
            ex = s.get("Exclusion")
            if ex is not None and ex["Class"] not in PERIMETER_EXCLUSIONS:
                continue
            physical.append(i)
            z = round(float(s["Key"][0][2]), 3)
            for p in s["Key"]:
                by_plane[z].setdefault(point_key(p), []).append(i)
        seen = set()
        for i in physical:
            if i in seen:
                continue
            z = round(float(self.segments[i]["Key"][0][2]), 3)
            adjacency = by_plane[z]
            # walk the component (segments sharing points); a loop is a component where every point has degree 2
            stack = [i]
            component = []
            seen.add(i)
            while stack:
                j = stack.pop()
                component.append(j)
                for p in self.segments[j]["Key"]:
                    for k in adjacency[point_key(p)]:
                        if k not in seen:
                            seen.add(k)
                            stack.append(k)
            polygon, closed = self._order(component, adjacency)
            area = 0.0
            for a in range(len(polygon)):
                x0, y0 = polygon[a]
                x1, y1 = polygon[(a + 1) % len(polygon)]
                area += x0 * y1 - x1 * y0
            xs = [p[0] for p in polygon]
            ys = [p[1] for p in polygon]
            loop = {"Plane": z, "Segments": component, "Polygon": polygon, "Area": 0.5 * area, "Closed": closed,
                    "BBox": [min(xs), max(xs), min(ys), max(ys)], "Depth": 0, "Parent": None, "Body": None}
            for j in component:
                self.segment_loop[j] = len(self.loops)
            self.loops[len(self.loops)] = loop
        self._containment()
        self._metal_side(features)
        self._bodies()

    def _order(self, component, adjacency):
        """Polyline through the component's points (an ordered traversal; degree-2 loops close exactly)."""
        comp = set(component)
        start = component[0]
        k0, k1 = point_key(self.segments[start]["Key"][0]), point_key(self.segments[start]["Key"][1])
        polygon = [k0[:2], k1[:2]]
        used = {start}
        current, previous = k1, k0
        closed = False
        while True:
            nxt = None
            for j in adjacency[current]:
                if j in comp and j not in used:
                    nxt = j
                    break
            if nxt is None:
                break
            used.add(nxt)
            a, b = point_key(self.segments[nxt]["Key"][0]), point_key(self.segments[nxt]["Key"][1])
            other = b if a == current else a
            if other == k0 and len(used) == len(comp):
                closed = True
                break
            polygon.append(other[:2])
            previous, current = current, other
        if len(used) != len(comp):
            closed = False
        return polygon, closed

    def _containment(self):
        by_plane = collections.defaultdict(list)
        for idx, loop in self.loops.items():
            by_plane[loop["Plane"]].append(idx)
        for plane, indices in by_plane.items():
            indices.sort(key=lambda k: abs(self.loops[k]["Area"]))  # smallest first
            for pos, idx in enumerate(indices):
                loop = self.loops[idx]
                x, y = loop["Polygon"][0]
                # containment of the loop's first vertex in every larger CLOSED loop (vertices never coincide
                # across loops); an open component (interrupted by walls) is neither contained nor a container
                parent = None
                if not loop["Closed"]:
                    loop["Parent"] = None
                    continue
                for cand in indices[pos + 1:]:
                    c = self.loops[cand]
                    if not c["Closed"]:
                        continue
                    bb = c["BBox"]
                    if x < bb[0] or x > bb[1] or y < bb[2] or y > bb[3]:
                        continue
                    if _point_in_polygon(x, y, c["Polygon"]):
                        parent = cand
                        break
                loop["Parent"] = parent
            for idx in indices:
                depth, p = 0, self.loops[idx]["Parent"]
                while p is not None:
                    depth += 1
                    p = self.loops[p]["Parent"]
                self.loops[idx]["Depth"] = depth

    def _bodies(self):
        """Body of a loop: itself when the metal is inside (an island, or a plane's outer boundary), else the
        nearest enclosing loop with metal inside; a hole with no such ancestor belongs to the plane's unbounded
        ground (body key ("Ground", plane): the L1 ground of a flip chip whose outer boundary is not in the
        perimeter table). Open components are their own record."""
        self.ground_bodies = {}
        for idx, loop in self.loops.items():
            if not loop["Closed"] or loop["MetalInside"]:
                loop["Body"] = idx
                continue
            p = loop["Parent"]
            while p is not None and not self.loops[p]["MetalInside"]:
                p = self.loops[p]["Parent"]
            if p is not None:
                loop["Body"] = p
            else:
                key = ("Ground", loop["Plane"])
                if key not in self.ground_bodies:
                    self.ground_bodies[key] = -1 - len(self.ground_bodies)
                loop["Body"] = self.ground_bodies[key]

    def _metal_side(self, features):
        """Metal side of every closed loop: the sharp corner features on it vote (on a loop with metal inside a
        polygon-convex vertex is a ConvexCorner of the metal, with metal outside a ConcaveCorner); a loop
        without sharp corners takes the opposite side of its parent (parity), and a parentless one is a hole
        of the plane's unbounded ground. Minority votes are reported as MetalSideDisagreements."""
        vertex_index = {}
        for loop_index, loop in self.loops.items():
            loop["MetalInside"] = None
            loop["CornerVotes"] = [0, 0]  # [metal inside, metal outside]
            if not loop["Closed"] or len(loop["Polygon"]) < 3:
                continue
            for k, p in enumerate(loop["Polygon"]):
                vertex_index[(loop["Plane"], p)] = (loop_index, k)
        votes = []
        for f in features:
            if f["Type"] not in ("ConvexCorner", "ConcaveCorner") or not f.get("Frame"):
                continue
            if float(f.get("Signature", {}).get("CornerRadiusOverR", 0.0) or 0.0) > 0.0:
                continue
            origin = f["Frame"]["Origin"]
            hit = vertex_index.get((round(float(origin[2]), 3), point_key(origin)[:2]))
            if hit is None:
                continue
            loop_index, k = hit
            loop = self.loops[loop_index]
            poly = loop["Polygon"]
            prev_p, v, next_p = poly[k - 1], poly[k], poly[(k + 1) % len(poly)]
            cross = (v[0] - prev_p[0]) * (next_p[1] - v[1]) - (v[1] - prev_p[1]) * (next_p[0] - v[0])
            polygon_convex = (cross > 0.0) == (loop["Area"] > 0.0)
            metal_inside = polygon_convex == (f["Type"] == "ConvexCorner")
            loop["CornerVotes"][0 if metal_inside else 1] += 1
            votes.append((loop_index, metal_inside, f["Id"], origin, f["Type"]))
            self.metal_side_checks += 1
        for loop in self.loops.values():
            inside, outside = loop["CornerVotes"]
            if inside or outside:
                loop["MetalInside"] = inside >= outside
        # loops without votes: parity from the nearest decided ancestor (top down: parents before children)
        order = sorted(self.loops, key=lambda i: self.loops[i]["Depth"])
        for idx in order:
            loop = self.loops[idx]
            if not loop["Closed"] or loop["MetalInside"] is not None:
                continue
            p = loop["Parent"]
            loop["MetalInside"] = (not self.loops[p]["MetalInside"]) if p is not None else False
        for loop_index, metal_inside, fid, origin, ftype in votes:
            if metal_inside != self.loops[loop_index]["MetalInside"]:
                self.metal_side_disagreements.append({"Feature": fid, "Type": ftype, "Loop": loop_index,
                                                      "Depth": self.loops[loop_index]["Depth"], "Origin": origin})


def pair_sides(feature, segments):
    """Portions grouped by side for pairs (``Sides`` per portion) or stacks (``Sides`` per portion)."""
    sides = feature.get("Sides")
    portions = feature.get("Portions", [])
    if sides is None or len(sides) != len(portions):
        return None
    groups = collections.defaultdict(list)
    for side, portion in zip(sides, portions):
        groups[side].append(portion)
    return groups


def portion_points(portion, segments):
    seg = segments[portion[0]]
    p0, p1 = seg["Key"][0], seg["Key"][1]
    length = float(seg["Length"])
    if length <= 0.0:
        return None
    t0, t1 = portion[1] / length, portion[2] / length
    a = [p0[k] + t0 * (p1[k] - p0[k]) for k in range(3)]
    b = [p0[k] + t1 * (p1[k] - p0[k]) for k in range(3)]
    return a, b


def mean_separation_prediction(feature, inside_portions, segments, radius):
    """Window-interior mean separation of side 0 to the nearest other side vs the chip signature."""
    groups = pair_sides(feature, segments)
    if not groups or len(groups) < 2:
        return None
    signature = feature.get("Signature", {})
    chip = signature.get("SeparationOverR")
    if chip is None and signature.get("Edges"):
        offsets = sorted(float(e.get("OffsetOverR", 0.0)) for e in signature["Edges"])
        chip = offsets[1] - offsets[0] if len(offsets) > 1 else None
    if chip is None:
        return None
    inside = {tuple(p) for p in inside_portions}
    side_keys = sorted(groups)
    partner = []
    for key in side_keys[1:]:
        for portion in groups[key]:
            pts = portion_points(portion, segments)
            if pts:
                partner.append(pts)
    if not partner:
        return None
    samples, total = [], 0.0
    spacing = 0.5 * radius
    for portion in groups[side_keys[0]]:
        if tuple(portion) not in inside:
            continue
        pts = portion_points(portion, segments)
        if not pts:
            continue
        a, b = pts
        length = math.dist(a[:2], b[:2])
        n = max(1, int(math.ceil(length / spacing)))
        for k in range(n):
            t = (k + 0.5) / n
            x, y = a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])
            d = min(_segment_distance(x, y, q[0], q[1]) for q in partner)
            samples.append(d)
    if not samples:
        return None
    mean = sum(samples) / len(samples)
    chip_sep = chip * radius
    return {"ChipSeparation": chip_sep, "WindowInteriorMean": mean, "WindowInteriorMin": min(samples),
            "WindowInteriorMax": max(samples), "Samples": len(samples),
            "DifferenceOverR": (mean - chip_sep) / radius,
            "KeyMayChange": abs(mean - chip_sep) / radius > SIGNATURE_TOLERANCE_OVER_R}


def joint_turn_is_noise(prev, vertex, nxt, noise_sagitta):
    """CornerRule: a joint turning by t between two pieces is noise when (c / 2) tan(t / 4) of the SHORTER piece c is below
    JointNoiseSagittaOverR x R (``noise_sagitta``); such a joint is never a corner."""
    d0 = (vertex[0] - prev[0], vertex[1] - prev[1])
    d1 = (nxt[0] - vertex[0], nxt[1] - vertex[1])
    l0, l1 = math.hypot(*d0), math.hypot(*d1)
    if l0 <= 0.0 or l1 <= 0.0:
        return True
    cos_t = max(-1.0, min(1.0, (d0[0] * d1[0] + d0[1] * d1[1]) / (l0 * l1)))
    return (min(l0, l1) / 2.0) * math.tan(math.acos(cos_t) / 4.0) < noise_sagitta


class ClipRegion:
    """The clip region of one plane of a window: the box (the setback box for the open-terminated terminal) minus the wall bridge
    strips of that plane (``window_polygons.py`` ``WindowPolygons.Bridges``: rectangles [x0, x1, y0, y1] of metal over the gaps
    next to a bridged conductor's wall interval, BridgeWidth = 3 R into the box). A chip chord vertex strictly inside a strip is
    swallowed by the bridge metal, so the chord is cut on the strip's edge; a vertex on a box wall or on a strip edge is inside."""

    def __init__(self, box, bridges=(), tol=1e-9):
        self.box = tuple(float(v) for v in box)
        self.bridges = [tuple(float(v) for v in r[:4]) for r in bridges]
        self.tol = tol

    def in_box(self, p):
        b, tol = self.box, self.tol
        return b[0] - tol <= p[0] <= b[1] + tol and b[2] - tol <= p[1] <= b[3] + tol

    def bridge_index(self, p):
        """The bridge strip holding ``p`` strictly inside, else None."""
        tol = self.tol
        for k, r in enumerate(self.bridges):
            if r[0] + tol < p[0] < r[1] - tol and r[2] + tol < p[1] < r[3] - tol:
                return k
        return None

    def contains(self, p):
        return self.in_box(p) and self.bridge_index(p) is None

    def exit(self, inside_point, outside_point):
        """Where the chord inside -> outside leaves the region and on what: (point, ``"x = ..."`` / ``"y = ..."`` for a box wall,
        ``"bridge k"`` for the first bridge strip the chord enters)."""
        t_exit, source = 1.0, "chord end"
        for axis, value in ((0, self.box[0]), (0, self.box[1]), (1, self.box[2]), (1, self.box[3])):
            if (inside_point[axis] - value) * (outside_point[axis] - value) < 0.0:
                t = (value - inside_point[axis]) / (outside_point[axis] - inside_point[axis])
                if t < t_exit:
                    t_exit, source = t, f"{'xy'[axis]} = {value:g}"
        for k, r in enumerate(self.bridges):
            interval = clip_interval(inside_point, outside_point, r)
            if interval is not None and interval[1] > interval[0] + self.tol and interval[0] < t_exit:
                t_exit, source = interval[0], f"bridge {k}"
        return (inside_point[0] + t_exit * (outside_point[0] - inside_point[0]), inside_point[1] + t_exit * (outside_point[1] - inside_point[1])), source


def cut_arm_is_tangent_noise(vertex, outside_neighbour, first_inside_neighbour, arc_radius, clip, noise_sagitta):
    """ArcRule end-joint test at a cut-adjacent joint: its arm is the cut chord (vertex -> the ``clip`` region's exit), which meets
    the circle's tangent at the vertex at half the chord's central angle; noise on the shorter of the arm and the first chord
    keeps the joint. Returns the record of the test: {Joint, ExitOn, Arm, KinkDegrees, ImpliedSagitta, Kept}."""
    exit_point, source = clip.exit(vertex, outside_neighbour)
    arm = math.dist(vertex[:2], exit_point)
    first_chord = math.dist(vertex[:2], first_inside_neighbour[:2]) if first_inside_neighbour is not None else arm
    kink = math.asin(max(-1.0, min(1.0, math.dist(vertex[:2], outside_neighbour[:2]) / (2.0 * arc_radius))))
    sagitta = (min(arm, first_chord) / 2.0) * math.tan(kink / 4.0)
    return {"Joint": [vertex[0], vertex[1]], "ExitOn": source, "Arm": arm, "KinkDegrees": math.degrees(kink),
            "ImpliedSagitta": sagitta, "Kept": sagitta < noise_sagitta}


def _arc_vertices(segment_indices, segments):
    """The arc's chord vertices in chain order (a closed arc returns first == last)."""
    adjacency = collections.defaultdict(list)
    for i in segment_indices:
        a, b = tuple(segments[i]["Key"][0][:2]), tuple(segments[i]["Key"][1][:2])
        adjacency[a].append(b)
        adjacency[b].append(a)
    ends = [v for v, n in adjacency.items() if len(n) == 1]
    start = ends[0] if ends else min(adjacency)
    out, prev = [start], None
    while True:
        nxt = [v for v in adjacency[out[-1]] if v != prev]
        if not nxt:
            break
        prev = out[-1]
        out.append(nxt[0])
        if nxt[0] == start or len(out) > len(segment_indices) + 1:
            break
    return out


def _inside_runs(vertices, clip):
    """Maximal runs of consecutive chord vertices inside the ``clip`` region: [(run, cut-adjacent vertices, their outside
    neighbours)]."""
    closed = len(vertices) > 1 and vertices[0] == vertices[-1]
    ring = vertices[:-1] if closed else vertices
    n = len(ring)
    inside = [clip.contains(v) for v in ring]
    if all(inside):
        return [(ring, [], [])]
    start = 0
    if closed:
        while inside[start]:
            start += 1
    runs, current = [], []
    for k in range(n):
        idx = (start + k) % n
        if inside[idx]:
            current.append(idx)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    out = []
    for run in runs:
        cut_adjacent, outside = [], []
        if closed or run[0] != 0:
            cut_adjacent.append(ring[run[0]])
            outside.append(ring[(run[0] - 1) % n])
        if (closed or run[-1] != n - 1) and run[-1] != run[0]:
            cut_adjacent.append(ring[run[-1]])
            outside.append(ring[(run[-1] + 1) % n])
        out.append(([ring[k] for k in run], cut_adjacent, outside))
    return out


def cut_arc_exclusions(ident, conductor_of_segment, box, setback_box, terminal, bridges=(), margin_over_r=CUT_ARC_MARGIN_OVER_R):
    """E1 cut-arc exclusions of one window (VALIDATION-PLAN (h)-9; module docstring): every chip arc cut by a wall (the
    setback line for the open-terminated ``terminal``) or by the edge of a wall bridge strip (``bridges``: [(plane z,
    rectangle [x0, x1, y0, y1])] from the polygon set's ``WindowPolygons.Bridges``), its predicted window reading (``Collapse``
    / ``DroppedJoint`` / ``None``), the pieces to exclude (absolute coordinates, ``margin_over_r`` R along the chain) and the
    setback hints."""
    segments = ident["Segments"]
    arcs = ident.get("Arcs", [])
    radius = float(ident["MatchingRadius"])
    conventions = ident.get("Conventions") or {}
    noise_over_r = float(conventions.get("JointNoiseSagittaOverR", JOINT_NOISE_SAGITTA_OVER_R_DEFAULT))
    noise_sagitta = noise_over_r * radius
    margin = margin_over_r * radius
    at_point = collections.defaultdict(list)
    by_arc = collections.defaultdict(list)
    for i, s in enumerate(segments):
        for p in s["Key"]:
            at_point[(point_key(p[:2]), round(float(p[2]), 3))].append(i)
        if s.get("Arc") is not None:
            by_arc[s["Arc"]].append(i)
    records, hints = [], []
    total_excluded = 0.0

    def piece(i, s0, s1):
        s = segments[i]
        p0, p1 = s["Key"][0], s["Key"][1]
        length = float(s["Length"])
        t0, t1 = (s0 / length, s1 / length) if length > 0.0 else (0.0, 1.0)
        return [[p0[k] + t0 * (p1[k] - p0[k]) for k in range(3)], [p0[k] + t1 * (p1[k] - p0[k]) for k in range(3)]]

    def along_chain(vertex, plane, forbidden, budget, pieces, truncations, through_arc=None):
        """Follow the chain from a vertex through segments not in ``forbidden`` for ``budget`` um, collecting pieces. The walk
        stops short (less exclusion: fails closed) at a branching vertex, at a chain end, or after the first chord of an arc
        other than ``through_arc``; the stop is appended to ``truncations`` with the margin not taken."""
        v, seen = vertex, set(forbidden)
        while budget > 1e-9:
            candidates = [j for j in at_point[(point_key(v[:2]), plane)] if j not in seen]
            if len(candidates) != 1:
                truncations.append({"At": [v[0], v[1], v[2] if len(v) > 2 else plane], "MarginNotTaken": budget,
                                    "Reason": "chain end" if not candidates else f"branching vertex ({len(candidates)} segments)"})
                return
            j = candidates[0]
            seen.add(j)
            s = segments[j]
            length = float(s["Length"])
            take = min(length, budget)
            starts_here = point_key(s["Key"][0][:2]) == point_key(v[:2])
            pieces.append(piece(j, 0.0, take) if starts_here else piece(j, length - take, length))
            budget -= take
            v = s["Key"][1] if starts_here else s["Key"][0]
            if s.get("Arc") is not None and s["Arc"] != through_arc and budget > 1e-9:
                truncations.append({"At": [v[0], v[1], v[2]], "MarginNotTaken": budget, "Reason": f"next arc {s['Arc']}"})
                return

    for arc_id, arc_segments in sorted(by_arc.items()):
        conductors = {conductor_of_segment(i) for i in arc_segments}
        is_terminal = terminal is not None and conductors == {terminal}
        clip_kind = "setback" if (is_terminal and setback_box is not None) else "wall"
        plane = round(float(segments[arc_segments[0]]["Key"][0][2]), 3)
        plane_bridges = [r for z, r in bridges if abs(float(z) - plane) < 1e-6]
        clip = ClipRegion(setback_box if clip_kind == "setback" else box, plane_bridges)
        vertices = _arc_vertices(arc_segments, segments)
        runs = _inside_runs(vertices, clip)
        n_inside = sum(len(run) for run, _, _ in runs)
        if n_inside == 0 or (len(runs) == 1 and not runs[0][1] and n_inside == len(set(vertices))):
            continue  # entirely outside or entirely inside: not cut
        arc = arcs[arc_id]
        arc_radius = float(arc["Radius"])
        closed = len(vertices) > 1 and vertices[0] == vertices[-1]
        ring = vertices[:-1] if closed else vertices
        interior = range(len(ring)) if closed else range(1, len(ring) - 1)
        real_turns = any(not joint_turn_is_noise(ring[k - 1], ring[k], ring[(k + 1) % len(ring)], noise_sagitta) for k in interior)
        dropped, kept, cut_adjacent_records = [], [], []
        for run, cut_adjacent, outside in runs:
            d = []
            for v, u in zip(cut_adjacent, outside):
                others = [w for w in run if w != v]
                nearest = min(others, key=lambda q: math.dist(q, v)) if others else None
                test = cut_arm_is_tangent_noise(v, u, nearest, arc_radius, clip, noise_sagitta)
                cut_adjacent_records.append(test)
                if not test["Kept"]:
                    d.append(v)
            dropped.append(d)
            kept.append(len(run) - len(d))
        lines = []
        for axis, value in ((0, clip.box[0]), (0, clip.box[1]), (1, clip.box[2]), (1, clip.box[3])):
            values = [v[axis] for v in ring]
            if min(values) < value - 1e-9 < max(values):
                lines.append({"Line": f"{'xy'[axis]} = {value:g}", "OutwardToFree": value - min(values), "InwardToFree": max(values) - value})
        # the third window-introduced cut: a cut arm ending on a bridge strip's edge, or a chord vertex swallowed by a strip
        swallowed = {clip.bridge_index(v) for v in ring if clip.in_box(v)} - {None}
        bridge_cuts = sorted(swallowed | {int(t["ExitOn"].split()[1]) for t in cut_adjacent_records if t["ExitOn"].startswith("bridge")})
        if bridge_cuts:
            clip_kind += "+bridge"
        pieces, truncations = [], []
        if not real_turns:
            prediction, reason = "None", "every joint of the arc is noise under the CornerRule: straight edges in both readings"
        elif any(k < LEAST_SQUARES_BEND_JOINTS for k in kept):
            prediction = "Collapse"
            reason = f"fewer than {LEAST_SQUARES_BEND_JOINTS} joints kept inside the clip region: every kept joint reads as a corner"
            for i in arc_segments:
                pieces.append(piece(i, 0.0, float(segments[i]["Length"])))
            if not closed:
                for end in (vertices[0], vertices[-1]):
                    along_chain(end, plane, arc_segments, margin, pieces, truncations)
        elif any(dropped):
            prediction = "DroppedJoint"
            reason = "a cut-adjacent joint whose cut arm is not tangent within the noise rule leaves the bend and reads as a corner"
            for d in dropped:
                for v in d:
                    # the full margin along the chain on each side of the dropped joint (its own chords included), as the
                    # Collapse branch takes it along the arms
                    at_v = at_point[(point_key(v[:2]), plane)]
                    for j in at_v:
                        along_chain(v, plane, [k for k in at_v if k != j], margin, pieces, truncations, through_arc=arc_id)
        else:
            prediction, reason = "None", "every cut-adjacent arm is tangent within the noise rule and at least four joints are kept"
        excluded_length = sum(math.dist(p[0][:2], p[1][:2]) for p in pieces)
        total_excluded += excluded_length
        records.append({"Arc": arc_id, "Kind": arc["Kind"], "Radius": arc_radius, "Joints": arc["Joints"],
                        "TurnDegrees": arc["TurnDegrees"], "Plane": plane, "Conductor": sorted(conductors, key=str),
                        "Terminal": is_terminal, "Clip": clip_kind, "Lines": lines,
                        "BridgeCuts": [{"Bridge": k, "Rectangle": list(clip.bridges[k])} for k in bridge_cuts],
                        "CutAdjacent": cut_adjacent_records,
                        "Segments": arc_segments, "InsideJoints": [len(run) for run, _, _ in runs],
                        "DroppedCutAdjacent": [len(d) for d in dropped], "KeptJoints": kept, "RealTurns": real_turns,
                        "Prediction": prediction, "Reason": reason, "ExcludedLength": excluded_length, "Pieces": pieces,
                        "MarginTruncations": truncations})
        if is_terminal and prediction == "Collapse" and setback_box is not None:
            for line in lines:
                hints.append({"Arc": arc_id, "Line": line["Line"], "InwardToFree": line["InwardToFree"],
                              "SetbackToFree": (setback_box[0] - box[0]) + line["InwardToFree"],
                              "Note": "the setback that cuts the arc away whole (the terminal then ends on the arm beyond it); "
                                      "information for the window placement, nothing is moved"})
    return {"Rule": "VALIDATION-PLAN (h)-9 cut-arc exclusion: Collapse = the arc's chords + margin along both arms; DroppedJoint = "
                    "margin along the chain on both sides of the dropped cut-adjacent joint; None = nothing; the margin walk stops "
                    "short at a branching vertex / chain end / the next arc (MarginTruncations)",
            "ClipModel": "box (the setback box for the open-terminated terminal) minus the wall bridge strips of the plane",
            "Bridges": [{"Plane": z, "Rectangle": list(r[:4])} for z, r in bridges],
            "MarginOverR": margin_over_r, "JointNoiseSagittaOverR": noise_over_r, "LeastSquaresBendJoints": LEAST_SQUARES_BEND_JOINTS,
            "CutArcs": len(records), "Predictions": dict(collections.Counter(r["Prediction"] for r in records)),
            "ExcludedLength": total_excluded, "Arcs": records, "SetbackHints": hints}


def polygon_set_bridges(polygon_set):
    """(window name, [(plane z, bridge rectangle)]) of a ``window_polygons.py`` polygon set (``WindowPolygons.Bridges`` named by
    plane, resolved to the plane's ``SurfaceZ``)."""
    z_of_plane = {p["Name"]: float(p["SurfaceZ"]) for p in polygon_set["Planes"]}
    bridges = [(z_of_plane[b["Plane"]], [float(v) for v in b["Rectangle"]]) for b in polygon_set["WindowPolygons"]["Bridges"]]
    return polygon_set["Name"], bridges


def inventory(manifest, windows, margin, weights, z_range=None, bridges=None):
    """``bridges``: {window name: [(plane z, bridge rectangle)]} (``polygon_set_bridges``) for the cut-arc clip model; a window
    without an entry is clipped by the box alone (conservative: no bridge edge shortens a wall-cut arm)."""
    ident = manifest["Identification"]
    radius = float(ident["MatchingRadius"])
    segments = ident["Segments"]
    features = ident["Features"]
    feature_by_id = {f["Id"]: f for f in features}
    from_extract = "WindowExtract" in manifest
    loops = PerimeterLoops.from_extract(manifest) if from_extract else PerimeterLoops(segments, features, radius)
    if margin is None:
        margin = 10.0 * radius
    results = {"MatchingRadius": radius, "Margin": margin, "Weights": weights, "Windows": {},
               "BodyModel": "chip (window_extract: bump-joined conductors)" if from_extract else "manifest (per plane, no bumps)",
               "PerimeterLoops": {"Count": len(loops.loops), "Closed": sum(1 for l in loops.loops.values() if l["Closed"]),
                                  "MetalSideChecks": loops.metal_side_checks,
                                  "MetalSideDisagreements": loops.metal_side_disagreements[:20],
                                  "MetalSideDisagreementCount": len(loops.metal_side_disagreements)}}
    open_setback = 3.0 * radius + margin
    for name, box in windows.items():
        x0, x1, y0, y1 = box
        inner = (x0 + margin, x1 - margin, y0 + margin, y1 - margin)
        setback_box = (x0 + open_setback, x1 - open_setback, y0 + open_setback, y1 - open_setback)
        per_feature = collections.defaultdict(lambda: {"Inside": 0.0, "InMargin": 0.0, "Portions": []})
        excluded = collections.defaultdict(float)
        per_plane = collections.defaultdict(lambda: {"Assigned": 0.0, "Excluded": 0.0})
        loops_inside = collections.defaultdict(lambda: {"Inside": 0.0, "InsideSetback": 0.0, "Crossing": False, "PortSegments": []})
        for i, s in enumerate(segments):
            p0, p1 = s["Key"][0], s["Key"][1]
            z = round(float(p0[2]), 3)
            if z_range and not (z_range[0] <= z <= z_range[1]):
                continue
            clip = clip_interval(p0, p1, (x0, x1, y0, y1))
            if clip is None:
                continue
            length = float(s["Length"])
            t0, t1 = clip
            a, b = t0 * length, t1 * length
            if b - a <= 0.0 and not (t0 == 0.0 and t1 == 1.0):
                continue
            inner_clip = clip_interval(p0, p1, inner)
            inner_a, inner_b = (inner_clip[0] * length, inner_clip[1] * length) if inner_clip else (0.0, 0.0)
            loop_index = loops.segment_loop.get(i)
            if loop_index is not None:
                entry = loops_inside[loop_index]
                entry["Inside"] += b - a
                setback_clip = clip_interval(p0, p1, setback_box)
                if setback_clip is not None:
                    entry["InsideSetback"] += (setback_clip[1] - setback_clip[0]) * length
                if (b - a) < length * (1.0 - 1e-9):
                    entry["Crossing"] = True
                if s.get("Exclusion", {}).get("Class") == "Port":
                    entry["PortSegments"].append(i)
            if "Exclusion" in s:
                excluded[s["Exclusion"]["Class"]] += b - a
                per_plane[z]["Excluded"] += b - a
                continue
            for s0, s1, fid in s["Portions"]:
                lo, hi = max(a, s0), min(b, s1)
                if hi <= lo:
                    continue
                entry = per_feature[fid]
                entry["Inside"] += hi - lo
                entry["Portions"].append([i, lo, hi])
                in_lo, in_hi = max(inner_a, lo), min(inner_b, hi)
                entry["InMargin"] += (hi - lo) - max(0.0, in_hi - in_lo)
                per_plane[z]["Assigned"] += hi - lo
            for s0, s1, ex_index in s.get("ExcludedPortions", []):
                lo, hi = max(a, s0), min(b, s1)
                if hi > lo:
                    cls = ident["Exclusions"][ex_index]["Class"]
                    excluded[cls] += hi - lo
                    per_plane[z]["Excluded"] += hi - lo
        rows = []
        by_key = collections.defaultdict(lambda: {"Features": 0, "Inside": 0.0, "Proxy": 0.0, "Planes": set(), "Cut": 0})
        predictions = []
        for fid, entry in per_feature.items():
            f = feature_by_id[fid]
            total = float(f.get("Length", 0.0))
            z = round(float(f["Frame"]["Origin"][2]), 3) if f.get("Frame") else None
            weight = weights.get(f["Type"], 1.0)
            cut = entry["Inside"] < total * (1.0 - 1e-9)
            row = {"Feature": fid, "Type": f["Type"], "Hash": f["Hash"][:12], "Plane": z, "LengthInside": entry["Inside"],
                   "Length": total, "CutByWall": cut, "LengthInMargin": entry["InMargin"], "Proxy": entry["Inside"] * weight,
                   "Origin": f["Frame"]["Origin"] if f.get("Frame") else None}
            sig = f.get("Signature", {})
            for k in ("EdgeCount", "AngleDegrees", "CornerRadiusOverR", "SeparationOverR", "RadiusOverR"):
                if k in sig:
                    row[k] = sig[k]
            rows.append(row)
            key = (f["Type"], f["Hash"][:12])
            by_key[key]["Features"] += 1
            by_key[key]["Inside"] += entry["Inside"]
            by_key[key]["Proxy"] += entry["Inside"] * weight
            by_key[key]["Planes"].add(z)
            by_key[key]["Cut"] += int(cut)
            if f["Type"] in ("ParallelEdgeCluster", "SameConductorStrip", "SameConductorGap", "DifferentConductorGap",
                             "CurvedSameConductorStrip", "CurvedSameConductorGap", "CurvedDifferentConductorGap"):
                pred = mean_separation_prediction(f, entry["Portions"], segments, radius)
                if pred:
                    pred.update({"Feature": fid, "Type": f["Type"], "Hash": f["Hash"][:12], "CutByWall": cut})
                    predictions.append(pred)
        rows.sort(key=lambda r: -r["Proxy"])
        keys = [{"Type": k[0], "Hash": k[1], "Features": v["Features"], "LengthInside": v["Inside"], "Proxy": v["Proxy"],
                 "Planes": sorted(p for p in v["Planes"] if p is not None), "CutFeatures": v["Cut"]}
                for k, v in sorted(by_key.items(), key=lambda kv: -kv[1]["Proxy"])]
        # bodies
        bodies = collections.defaultdict(lambda: {"Loops": [], "PerimeterInside": 0.0, "PerimeterInsideSetback": 0.0,
                                                  "Crossing": False, "Plane": None, "Kind": None, "Proxy": 0.0})
        def body_kind(body_index):
            if body_index < 0:
                return "Ground"
            body = loops.loops[body_index]
            if not body["Closed"]:
                return "OpenChain"
            return "Ground" if body["Parent"] is None else "Island"

        for loop_index, entry in loops_inside.items():
            loop = loops.loops[loop_index]
            body_index = loop["Body"]
            b = bodies[body_index]
            b["Loops"].append(loop_index)
            b["PerimeterInside"] += entry["Inside"]
            b["PerimeterInsideSetback"] += entry["InsideSetback"]
            b["Crossing"] = b["Crossing"] or entry["Crossing"]
            b["Plane"] = loop["Plane"]
            b["Kind"] = body_kind(body_index)
            if body_index >= 0:
                body = loops.loops[body_index]
                b["Area"] = abs(body["Area"])
                b["BBox"] = body["BBox"]
                b["Closed"] = body["Closed"]
            else:
                b["Area"], b["BBox"], b["Closed"] = None, None, None
        # proxy energy per body from the features inside (feature -> first portion's segment -> loop -> body)
        for row in rows:
            f = feature_by_id[row["Feature"]]
            if f.get("Portions"):
                loop_index = loops.segment_loop.get(f["Portions"][0][0])
                if loop_index is not None and loops.loops[loop_index]["Body"] in bodies:
                    bodies[loops.loops[loop_index]["Body"]]["Proxy"] += row["Proxy"]
        body_rows = []
        for body_index, b in bodies.items():
            entirely_inside = (not b["Crossing"]) and b["Kind"] == "Island"
            body_rows.append({"Body": body_index, "Kind": b["Kind"], "Plane": b["Plane"], "Loops": len(b["Loops"]),
                              "PerimeterInside": b["PerimeterInside"], "PerimeterInsideSetback": b["PerimeterInsideSetback"],
                              "CutByWall": b["Crossing"], "EntirelyInside": entirely_inside, "Proxy": b["Proxy"], "Area": b.get("Area"),
                              "BBox": b.get("BBox"), "LoopClosed": b.get("Closed")})
        body_rows.sort(key=lambda r: (r["Kind"] != "Ground", r["Kind"] == "OpenChain", -r["Proxy"]))
        conductor_rows, footprint_rows = conductor_table(loops, body_rows, box)
        port_lines = port_closed_loops(loops, loops_inside, body_kind)
        terminal, assignment, excitation = terminal_assignment(body_rows, conductor_rows, radius, margin, box, port_lines)
        # E1 comparison region (supervisor reply to decision 184 / 190, 2026-10-01): 3 R around every WINDOW-INTRODUCED feature
        # (walls, wall bridges, open-terminated ends) is excluded -> for an open-terminated terminal the margin is OpenSetback + 3 R
        e1_margin = max(margin, open_setback + 3.0 * radius) if excitation["Kind"] == "OpenTerminatedTrace" else margin
        e1_region = (x0 + e1_margin, x1 - e1_margin, y0 + e1_margin, y1 - e1_margin)

        def conductor_of_segment(i):
            loop_index = loops.segment_loop.get(i)
            if loop_index is None:
                return None
            body = loops.loops[loop_index]["Body"]
            rec = loops.bodies.get(body) if loops.bodies else None
            return rec["Conductor"] if rec else body

        cut_arcs = cut_arc_exclusions(ident, conductor_of_segment, box,
                                      setback_box if excitation["Kind"] == "OpenTerminatedTrace" else None, terminal,
                                      bridges=(bridges or {}).get(name, ()))
        results["Windows"][name] = {
            "Box": list(box), "ComparisonRegion": list(inner),
            "E1ComparisonRegion": list(e1_region), "E1Margin": e1_margin, "E1CutArcExclusions": cut_arcs,
            "Perimeter": {"Total": sum(v["Assigned"] + v["Excluded"] for v in per_plane.values()),
                          "Assigned": sum(v["Assigned"] for v in per_plane.values()),
                          "Excluded": dict(sorted(excluded.items())),
                          "PerPlane": {str(z): v for z, v in sorted(per_plane.items())}},
            "FeatureCount": len(rows), "CutFeatureCount": sum(1 for r in rows if r["CutByWall"]),
            "ByType": {t: {"Features": sum(1 for r in rows if r["Type"] == t),
                           "LengthInside": sum(r["LengthInside"] for r in rows if r["Type"] == t),
                           "Proxy": sum(r["Proxy"] for r in rows if r["Type"] == t)}
                       for t in sorted({r["Type"] for r in rows})},
            "Keys": keys, "Features": rows,
            "Bodies": body_rows, "Conductors": conductor_rows, "ConductorCount": len(conductor_rows),
            "BumpFootprints": footprint_rows, "PortClosedLoops": port_lines,
            "ProposedTerminal": terminal, "Excitation": excitation,
            "TerminalAssignment": assignment,
            "MeanSeparationPredictions": predictions,
        }
    return results


def conductor_table(loops, body_rows, box):
    """Bodies of the window grouped into conductors: with a chip model (``PerimeterLoops.from_extract``) a
    conductor = the bodies joined by bump footprints anywhere on the chip (decision 184 M1), else every body is
    its own conductor. A conductor is Ground when one of its bodies (seen or unseen) is a ground, EntirelyInside
    when every body of it is an island entirely inside the window and the chip joins no body the window does not
    see. Also the bump footprints meeting the window (plane, height, bodies joined, whole / straddling); a closed
    NonManifold loop without a partner on the other plane is not a bump (the box wall outline) and is skipped."""
    groups = collections.OrderedDict()
    for r in body_rows:
        rec = loops.bodies.get(r["Body"]) if loops.bodies else None
        root = rec["Conductor"] if rec else r["Body"]
        g = groups.setdefault(root, {"Conductor": root, "Bodies": [], "Planes": set(), "Kind": None, "Proxy": 0.0,
                                      "PerimeterInside": 0.0, "PerimeterInsideSetback": 0.0, "CutByWall": False,
                                      "ChipBodies": rec["ConductorBodies"] if rec else 1,
                                      "ChipGround": rec["ConductorGround"] if rec else r["Kind"] == "Ground"})
        g["Bodies"].append(r["Body"])
        g["Planes"].add(r["Plane"])
        g["Proxy"] += r["Proxy"]
        g["PerimeterInside"] += r["PerimeterInside"]
        g["PerimeterInsideSetback"] += r.get("PerimeterInsideSetback", r["PerimeterInside"])
        g["CutByWall"] = g["CutByWall"] or r["CutByWall"]
        kinds = {kind for kind in (g["Kind"], r["Kind"]) if kind}
        g["Kind"] = "Ground" if "Ground" in kinds or g["ChipGround"] else ("OpenChain" if "OpenChain" in kinds else "Island")
    conductor_rows = []
    for g in groups.values():
        seen_all = len(g["Bodies"]) >= g["ChipBodies"]
        g["UnseenBodies"] = max(0, g["ChipBodies"] - len(g["Bodies"]))
        g["EntirelyInside"] = g["Kind"] == "Island" and not g["CutByWall"] and seen_all
        g["Planes"] = sorted(p for p in g["Planes"] if p is not None)
        conductor_rows.append(g)
    conductor_rows.sort(key=lambda g: (g["Kind"] != "Ground", g["Kind"] == "OpenChain", -g["Proxy"]))
    footprint_rows = []
    for fp in loops.footprints:
        bb = fp["BBox"]
        if fp.get("Partner") is None:  # a closed NonManifold loop without a partner is not a bump (the box wall outline)
            continue
        if bb[1] < box[0] or bb[0] > box[1] or bb[3] < box[2] or bb[2] > box[3]:
            continue
        partner_body = None
        if fp.get("Partner") is not None:
            partner = [q for q in loops.footprints if q.get("Index") == fp["Partner"]]
            partner_body = partner[0]["Body"] if partner else "unseen"
        footprint_rows.append({"Index": fp.get("Index"), "Plane": fp["Plane"], "Height": fp.get("Height"), "Body": fp["Body"],
                               "PartnerBody": partner_body, "Centroid": fp["Centroid"], "Area": fp["Area"],
                               "Inside": bool(fp.get("Inside")), "Straddles": bool(fp.get("Straddles"))})
    return conductor_rows, footprint_rows


def port_closed_loops(loops, loops_inside, body_kind):
    """Body-model convention (review of decision 195, m2): a lumped-port patch is not metal, so the metal edge along it is a
    ``Port`` segment; when it lies on a closed DIELECTRIC-inside loop, the metal beyond the port reads as the loop's
    exterior body. That metal is either the ground itself (the ground-side edges of a JJ port: SCT loop 326 / CTX loops 39,
    91) or a port-terminated LINE fused with the ground by the port closing its gap loop (C4's XY drive line above pad 43,
    loop 42; S5's stub below the ROA_IN port, loop 22; both read as the unbounded ground -2) -- the body model cannot tell
    the two apart. Returns one record per such loop meeting the window: the loop, the body (and kind) the metal beyond the
    port reads as, the Port segments; the terminal rule fails closed when such a loop is the only metal that could have
    been a candidate."""
    lines = []
    for loop_index, entry in sorted(loops_inside.items()):
        if not entry["PortSegments"]:
            continue
        loop = loops.loops[loop_index]
        if not loop["Closed"] or loop["MetalInside"] is not False:
            continue  # metal inside (a JJ island closed through its port edges) is a body of its own
        lines.append({"Loop": loop_index, "Plane": loop["Plane"], "ReadAsBody": loop["Body"], "ReadAsKind": body_kind(loop["Body"]),
                      "PortSegments": entry["PortSegments"]})
    return lines


def terminal_assignment(body_rows, conductor_rows, radius, margin, box, port_lines=None):
    """Decision-180 terminal rule per conductor: the island conductor entirely inside the box with the largest
    proxy energy is the quoted terminal (1 V); else the cut non-ground conductor (a CPW centre trace / cut island)
    with the largest proxy is the terminal, open-terminated inside the wall (its metal is clipped to the setback box, 3 R +
    margin inside every wall, so EVERY cut of the terminal is open-terminated — ``window_polygons.py`` bridges no cut of the
    terminal and leaves its neighbours unbridged on the terminal's side); every other conductor is ground, a cut one bridged to
    ground at the wall by an explicit metal polygon. An UNCUT island conductor with bump-joined bodies the window does
    not see (``EntirelyInside`` False without ``CutByWall``) is not a terminal candidate of either kind: its extent is
    unknown (``UncutWithUnseenBodies``). A window with no non-ground conductor has no excitation (S3 class: move /
    resize / drop); when a port-terminated line read as ground (``port_lines``, m2 convention) is then the only metal
    that could have been a candidate, the rule FAILS CLOSED with the diagnostic (``Flag``). Returns (terminal
    conductor root or None, per-body assignment, excitation)."""
    islands = [g for g in conductor_rows if g["EntirelyInside"]]
    cut_traces = [g for g in conductor_rows if g["Kind"] == "Island" and g["CutByWall"]]
    uncut_unseen = [g for g in conductor_rows if g["Kind"] == "Island" and not g["CutByWall"] and not g["EntirelyInside"]]
    ground_port_lines = [l for l in (port_lines or []) if l["ReadAsKind"] == "Ground"]  # m2: possibly a line fused with ground
    open_setback = 3.0 * radius + margin
    # a cut conductor whose metal edges all lie within the setback band vanishes when open-terminated (O4's wall-hugging
    # trace): not a realisable terminal
    realisable_traces = [g for g in cut_traces if g.get("PerimeterInsideSetback", g["PerimeterInside"]) > 0.0]
    vanishing = [g["Conductor"] for g in cut_traces if g not in realisable_traces]
    if islands:
        terminal = max(islands, key=lambda g: g["Proxy"])
        excitation = {"Kind": "Island", "Conductor": terminal["Conductor"], "Bodies": terminal["Bodies"], "Planes": terminal["Planes"],
                      "Realisable": True, "Rule": "whole island conductor (bump-joined bodies included) at 1 V"}
    elif realisable_traces:
        terminal = max(realisable_traces, key=lambda g: g["Proxy"])
        excitation = {"Kind": "OpenTerminatedTrace", "Conductor": terminal["Conductor"], "Bodies": terminal["Bodies"],
                      "Planes": terminal["Planes"], "Realisable": True, "OpenSetback": open_setback,
                      "Rule": f"cut non-ground conductor with the largest proxy, open-terminated {open_setback:g} um inside the wall"}
    else:
        terminal = None
        if cut_traces:
            rule = f"every cut non-ground conductor lies within {open_setback:g} um of a wall: move / resize / drop the window"
        elif uncut_unseen:
            rule = "the only non-ground conductors are uncut islands whose chip conductor has bump-joined bodies the window does " \
                   "not see (extent unknown): move / resize the window to see the whole conductor, or drop it"
        else:
            rule = "no non-ground conductor in the window: move / resize / drop the window"
        excitation = {"Kind": "None", "Conductor": None, "Bodies": [], "Planes": [], "Realisable": False, "Rule": rule}
        if ground_port_lines:
            excitation["Flag"] = "PortTerminatedLineReadAsGround"
            excitation["Rule"] = (f"FAIL CLOSED: the metal beyond the Port segments of {len(ground_port_lines)} gap loop(s) "
                                  f"{[l['Loop'] for l in ground_port_lines]} reads as a ground body and is the only metal that could "
                                  "have been a candidate; the body model cannot tell a port-terminated line (its gap loop closed by "
                                  "the port) from the ground: fuse the port in the chip model or move the window; " + rule)
    if vanishing:
        excitation["VanishingUnderSetback"] = vanishing
    if uncut_unseen:
        excitation["UncutWithUnseenBodies"] = [g["Conductor"] for g in uncut_unseen]
    if port_lines:
        excitation["PortClosedLoops"] = [{"Loop": l["Loop"], "ReadAsBody": l["ReadAsBody"], "ReadAsKind": l["ReadAsKind"]}
                                             for l in port_lines]
    terminal_bodies = set(terminal["Bodies"]) if terminal else set()
    assignment = []
    for r in body_rows:
        if r["Body"] in terminal_bodies:
            role = "Terminal (1 V, quoted)" if excitation["Kind"] == "Island" else \
                f"Terminal (1 V, quoted): trace open-terminated {open_setback:g} um inside the wall"
        elif r["Kind"] == "OpenChain":
            role = "unresolved: perimeter interrupted by walls / non-manifold edges (body from the mesh connectivity at cut time)"
        elif r["Kind"] == "Ground":
            role = "Ground (0 V)"
        elif r["CutByWall"]:
            role = "Ground (0 V): cut island / trace bridged to ground at the wall"
        else:
            conductor = [g for g in conductor_rows if r["Body"] in g["Bodies"]][0]
            if conductor["Kind"] == "Ground" or conductor["CutByWall"]:
                role = "Ground (0 V): island bump-joined to a ground / cut conductor"
            elif not conductor["EntirelyInside"]:
                role = "Ground (0 V): uncut island whose chip conductor has bump-joined bodies the window does not see (not a terminal candidate)"
            else:
                role = "Ground (0 V): non-excited island (capacitance-matrix convention)"
        assignment.append({"Body": r["Body"], "Kind": r["Kind"], "Plane": r["Plane"], "Role": role})
    return (terminal["Conductor"] if terminal else None), assignment, excitation


def write_csv(results, path):
    import csv
    with open(path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["Window", "Feature", "Type", "Hash", "Plane", "LengthInside", "Length", "CutByWall", "LengthInMargin",
                         "Proxy", "OriginX", "OriginY", "OriginZ", "EdgeCount", "AngleDegrees", "CornerRadiusOverR",
                         "SeparationOverR", "RadiusOverR"])
        for name, w in results["Windows"].items():
            for r in w["Features"]:
                o = r["Origin"] or [None, None, None]
                writer.writerow([name, r["Feature"], r["Type"], r["Hash"], r["Plane"], f"{r['LengthInside']:.6f}", f"{r['Length']:.6f}",
                                 int(r["CutByWall"]), f"{r['LengthInMargin']:.6f}", f"{r['Proxy']:.4f}", o[0], o[1], o[2],
                                 r.get("EdgeCount", ""), r.get("AngleDegrees", ""), r.get("CornerRadiusOverR", ""),
                                 r.get("SeparationOverR", ""), r.get("RadiusOverR", "")])


def markdown(results):
    lines = []
    for name, w in results["Windows"].items():
        p = w["Perimeter"]
        lines.append(f"### {name} box {w['Box']} (R = {results['MatchingRadius']}, margin {results['Margin']})")
        lines.append(f"P inside = {p['Total']:.1f} (assigned {p['Assigned']:.1f}; excluded {json.dumps({k: round(v, 1) for k, v in p['Excluded'].items()})}); "
                     f"per plane {json.dumps({z: {k: round(v, 1) for k, v in d.items()} for z, d in p['PerPlane'].items()})}")
        lines.append(f"features {w['FeatureCount']} ({w['CutFeatureCount']} cut by a wall); conductors {w['ConductorCount']}; proposed terminal body {w['ProposedTerminal']}")
        lines.append("| Type | features | length inside | proxy |")
        lines.append("|---|---|---|---|")
        for t, d in w["ByType"].items():
            lines.append(f"| {t} | {d['Features']} | {d['LengthInside']:.1f} | {d['Proxy']:.1f} |")
        lines.append("| body | kind | plane | perimeter inside | cut by wall | entirely inside | proxy | role |")
        lines.append("|---|---|---|---|---|---|---|---|")
        roles = {a["Body"]: a["Role"] for a in w["TerminalAssignment"]}
        for b in w["Bodies"]:
            lines.append(f"| {b['Body']} | {b['Kind']} | {b['Plane']} | {b['PerimeterInside']:.1f} | {b['CutByWall']} | {b['EntirelyInside']} | {b['Proxy']:.1f} | {roles[b['Body']]} |")
        lines.append("| key type | hash | features | length inside | proxy | planes | cut |")
        lines.append("|---|---|---|---|---|---|---|")
        for k in w["Keys"][:40]:
            lines.append(f"| {k['Type']} | {k['Hash']} | {k['Features']} | {k['LengthInside']:.1f} | {k['Proxy']:.1f} | {k['Planes']} | {k['CutFeatures']} |")
        if w["MeanSeparationPredictions"]:
            lines.append("| pair / stack | type | hash | chip separation | window mean | min | max | key may change | cut |")
            lines.append("|---|---|---|---|---|---|---|---|---|")
            for m in w["MeanSeparationPredictions"]:
                lines.append(f"| {m['Feature']} | {m['Type']} | {m['Hash']} | {m['ChipSeparation']:.4f} | {m['WindowInteriorMean']:.4f} | {m['WindowInteriorMin']:.4f} | {m['WindowInteriorMax']:.4f} | {m['KeyMayChange']} | {m['CutByWall']} |")
        lines.append("")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--window", nargs=5, action="append", metavar=("NAME", "X0", "X1", "Y0", "Y1"), required=True)
    parser.add_argument("--margin", type=float, default=None, help="exclusion margin from every wall (default 10 R)")
    parser.add_argument("--z", nargs=2, type=float, default=None, metavar=("ZMIN", "ZMAX"))
    parser.add_argument("--weights", default=None, help="JSON {Type: per-length proxy weight}")
    parser.add_argument("--polygons", action="append", default=[], metavar="POLYGON_SET",
                        help="window_polygons.py polygon set (its Name selects the window): the wall bridge strips of the cut-arc clip model")
    parser.add_argument("--output", required=True)
    parser.add_argument("--csv", default=None)
    parser.add_argument("--markdown", default=None)
    args = parser.parse_args(argv)
    weights = dict(DEFAULT_WEIGHTS)
    if args.weights:
        weights.update(json.load(open(args.weights)))
    windows = {w[0]: tuple(float(v) for v in w[1:]) for w in args.window}
    manifest = json.load(open(args.manifest))
    bridges = dict(polygon_set_bridges(json.load(open(path))) for path in args.polygons)
    unknown = sorted(set(bridges) - set(windows))
    if unknown:
        parser.error(f"--polygons for windows not requested: {unknown}")
    results = inventory(manifest, windows, args.margin, weights, args.z, bridges)
    results["Manifest"] = args.manifest
    results["PolygonSets"] = list(args.polygons)
    with open(args.output, "w") as out:
        json.dump(results, out, indent=1, default=list)
    if args.csv:
        write_csv(results, args.csv)
    text = markdown(results)
    if args.markdown:
        with open(args.markdown, "w") as out:
            out.write(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
