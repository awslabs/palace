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
  change the signature key (difference beyond the 1e-3 R signature tolerance).

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
        self.loops = []  # dicts: Plane, Segments, Polygon, Area, BBox, Depth, Parent, Body, Closed
        self._build(features)

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
            self.loops.append(loop)
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
        for idx, loop in enumerate(self.loops):
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
        for idx, loop in enumerate(self.loops):
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
        self.metal_side_checks = 0
        self.metal_side_disagreements = []
        vertex_index = {}
        for loop_index, loop in enumerate(self.loops):
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
        for loop in self.loops:
            inside, outside = loop["CornerVotes"]
            if inside or outside:
                loop["MetalInside"] = inside >= outside
        # loops without votes: parity from the nearest decided ancestor (top down: parents before children)
        order = sorted(range(len(self.loops)), key=lambda i: self.loops[i]["Depth"])
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


def inventory(manifest, windows, margin, weights, z_range=None):
    ident = manifest["Identification"]
    radius = float(ident["MatchingRadius"])
    segments = ident["Segments"]
    features = ident["Features"]
    feature_by_id = {f["Id"]: f for f in features}
    loops = PerimeterLoops(segments, features, radius)
    if margin is None:
        margin = 10.0 * radius
    results = {"MatchingRadius": radius, "Margin": margin, "Weights": weights, "Windows": {},
               "PerimeterLoops": {"Count": len(loops.loops), "Closed": sum(1 for l in loops.loops if l["Closed"]),
                                  "MetalSideChecks": loops.metal_side_checks,
                                  "MetalSideDisagreements": loops.metal_side_disagreements[:20],
                                  "MetalSideDisagreementCount": len(loops.metal_side_disagreements)}}
    for name, box in windows.items():
        x0, x1, y0, y1 = box
        inner = (x0 + margin, x1 - margin, y0 + margin, y1 - margin)
        per_feature = collections.defaultdict(lambda: {"Inside": 0.0, "InMargin": 0.0, "Portions": []})
        excluded = collections.defaultdict(float)
        per_plane = collections.defaultdict(lambda: {"Assigned": 0.0, "Excluded": 0.0})
        loops_inside = collections.defaultdict(lambda: {"Inside": 0.0, "Crossing": False})
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
                if (b - a) < length * (1.0 - 1e-9):
                    entry["Crossing"] = True
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
        bodies = collections.defaultdict(lambda: {"Loops": [], "PerimeterInside": 0.0, "Crossing": False, "Plane": None,
                                                  "Kind": None, "Proxy": 0.0})
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
                              "PerimeterInside": b["PerimeterInside"], "CutByWall": b["Crossing"],
                              "EntirelyInside": entirely_inside, "Proxy": b["Proxy"], "Area": b.get("Area"),
                              "BBox": b.get("BBox"), "LoopClosed": b.get("Closed")})
        body_rows.sort(key=lambda r: (r["Kind"] != "Ground", r["Kind"] == "OpenChain", -r["Proxy"]))
        islands_inside = [r for r in body_rows if r["EntirelyInside"]]
        terminal = max(islands_inside, key=lambda r: r["Proxy"]) if islands_inside else None
        assignment = []
        for r in body_rows:
            if terminal is not None and r["Body"] == terminal["Body"]:
                role = "Terminal (1 V, quoted)"
            elif r["Kind"] == "Ground":
                role = "Ground (0 V)"
            elif r["Kind"] == "OpenChain":
                role = "unresolved: perimeter interrupted by walls / non-manifold edges (body from the mesh connectivity at cut time)"
            elif r["CutByWall"]:
                role = "Ground (0 V): cut island / trace bridged to ground at the wall"
            else:
                role = "Ground (0 V): non-excited island (capacitance-matrix convention)"
            assignment.append({"Body": r["Body"], "Kind": r["Kind"], "Plane": r["Plane"], "Role": role})
        results["Windows"][name] = {
            "Box": list(box), "ComparisonRegion": list(inner),
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
            "Bodies": body_rows, "ConductorCount": len(body_rows),
            "ProposedTerminal": terminal["Body"] if terminal else None,
            "TerminalAssignment": assignment,
            "MeanSeparationPredictions": predictions,
        }
    return results


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
    parser.add_argument("--output", required=True)
    parser.add_argument("--csv", default=None)
    parser.add_argument("--markdown", default=None)
    args = parser.parse_args(argv)
    weights = dict(DEFAULT_WEIGHTS)
    if args.weights:
        weights.update(json.load(open(args.weights)))
    windows = {w[0]: tuple(float(v) for v in w[1:]) for w in args.window}
    manifest = json.load(open(args.manifest))
    results = inventory(manifest, windows, args.margin, weights, args.z)
    results["Manifest"] = args.manifest
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
