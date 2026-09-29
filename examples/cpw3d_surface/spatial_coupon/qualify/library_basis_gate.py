#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Library trace-basis gate (corner-family review 2026-09-29, supervisor decision 132): no
FREE trace basis function of any coupon may have support on the PEC part of the coupon's
matching-box contour. A free hat that straddles the metal edge imposes conflicting Dirichlet
data on the coupon solve (the trace on the box says V != 0 where the metal says 0) and puts a
near-singular field into the 2 nm MS / MA layers: the fabricated Q_MS diagonal of that knot
was 6-7 orders above a physical free knot on the 75 / 105 / 120 / 150 / 165-degree corner
coupons of the lane-2 layout, and the device's MS / MA over-corrected 4-6x while the domain
defect (fab - thin) cancelled the artefact. The gate is geometric (seconds, no solve) and is
evaluated from the model's basis points, its trace mesh where it has one, and the coupon
geometry the library records:

* corner coupons (ConvexCorner / ConcaveCorner with ContourGroups): on every box ring with
  PEC knots (ZeroTraceIndices), the two crossings of the metal arms (the first along +x, the
  second at Angle counterclockwise; a rounded corner's arms are straight where they cross the
  box) must be PEC knots and no free knot may lie on the closed metal footprint (the sector
  between the arms for a convex corner, its complement for a concave one). The same rule as
  the runtime's fail-closed library check (palace/models/cornertracebasis.cpp);
* spatial (3D box) coupons with a TraceMesh and a PlanViewBoundary (SpatialEdgeCluster,
  Endpoint, Junction): on every trace ring at a level that carries conductor vertices, every
  crossing of a conductor's plan-view polygon with the box perimeter must be a conductor
  vertex (basis 0, conductor > 0) or a PEC knot, and no free knot may lie inside a polygon;
* two-dimensional cross-section coupons (IsolatedEdge, CurvedEdge, the gap pairs): the metal
  band on the box side(s) (x = -R for an edge, both sides for a gap; 0 <= y <= MetalThickness)
  must contain no free knot, and the knots adjacent to the band must lie outside it (the
  generator's junction knots are constrained and absent from the basis); strips and
  parallel-edge clusters whose metal does not reach the box sides pass by construction and
  are recorded as such; a two-dimensional model whose metal pattern on the contour cannot be
  derived from its record is NotEvaluable, which FAILS the library verdict (fail closed).

The DIAGONAL sanity check of the review accompanies the geometry for every family (corner
family: the 90-degree node of the convexity, else the least-turn node, is the clean
reference; curvature families: the straight anchor): max over the free knots of the
fabricated within-R Q_MS / Q_MA diagonal must lie within DiagonalFactor (10) of the reference
node's; the factor is recorded per model.

    python3 library_basis_gate.py process-library.json [--json out.json]

Verdict Passed / Failed with every finding; exit code 1 on Failed.
"""
import argparse
import csv
import json
import math
import os
import sys

DIAGONAL_FACTOR = 10.0
FRACTION_TOLERANCE = 1.0e-6
PLAN_VIEW_SCALE = 1.0e-9  # PlanViewBoundary integers = coordinate / R x 1e9
CORNER_TOPOLOGIES = ("ConvexCorner", "ConcaveCorner")
SPATIAL_TOPOLOGIES = ("SpatialEdgeCluster", "Endpoint", "Junction")
EDGE_2D = ("IsolatedEdge", "CurvedEdge")
GAP_2D = ("SameConductorGap", "DifferentConductorGap", "CurvedSameConductorGap",
          "CurvedDifferentConductorGap")
STRIP_2D = ("SameConductorStrip", "CurvedSameConductorStrip")
CLUSTER_2D = ("ParallelEdgeCluster", "CurvedParallelEdgeCluster")


def resolve(library_path, library, value):
    if value is None:
        return None
    if os.path.isabs(value) and os.path.exists(value):
        return value
    for base in (os.path.dirname(os.path.abspath(library_path)), library.get("Root")):
        if base and os.path.exists(os.path.join(base, value)):
            return os.path.join(base, value)
    return None


def read_points(path):
    with open(path) as source:
        reader = csv.reader(source)
        next(reader)
        return [tuple(float(v) for v in row[:3]) for row in reader if row and row[0].strip()]


def read_trace_vertices(path):
    with open(path) as source:
        rows = list(csv.DictReader(source))
    vertices = []
    for row in rows:
        vertices.append({
            "point": (float(row["x"]), float(row["y"]), float(row["z"])),
            "basis": int(row["basis"]), "conductor": int(row["conductor"]),
            "parent_a": int(row.get("parent_a") or 0), "parent_b": int(row.get("parent_b") or 0),
        })
    return vertices


def read_surface_diagonal(path, size):
    """Per interface: {basis (1-based): within-R diagonal Q_ii} of a compact surface file
    (columns interface, edge, R (m), basis_i, basis_j, Q_ij (J)) or the whole-box Q_total for
    a file without the within-R column (translational coupons)."""
    with open(path) as source:
        rows = list(csv.DictReader(source))
    if not rows:
        return {}
    column = "Q_ij (J)" if "Q_ij (J)" in rows[0] else "Q_total_ij (J)"
    radii = sorted({float(r["R (m)"]) for r in rows if "R (m)" in r and r["R (m)"]})
    radius = radii[-1] if radii else None  # the matching radius = the largest recorded
    diagonal = {}
    for row in rows:
        if radius is not None and abs(float(row["R (m)"]) - radius) > 1.0e-9 * radius:
            continue
        i, j = int(row["basis_i"]), int(row["basis_j"])
        if i != j or i > size:
            continue
        key = (int(row["interface"]), i)
        diagonal[key] = diagonal.get(key, 0.0) + float(row[column])
    result = {}
    for (interface, i), value in diagonal.items():
        result.setdefault(interface, {})[i] = value
    return result


# ---------------------------------------------------------------------------- square rings

def square_fraction(half_width, point, tolerance):
    x, y = point[0], point[1]
    if abs(x + half_width) <= tolerance and y <= 0.0:
        c = -y
    elif abs(y + half_width) <= tolerance:
        c = half_width + (x + half_width)
    elif abs(x - half_width) <= tolerance:
        c = 3.0 * half_width + (y + half_width)
    elif abs(y - half_width) <= tolerance:
        c = 5.0 * half_width + (half_width - x)
    elif abs(x + half_width) <= tolerance:
        c = 7.0 * half_width + (half_width - y)
    else:
        return None
    return (c / (8.0 * half_width)) % 1.0


def square_point(half_width, z, fraction):
    c = (fraction % 1.0) * 8.0 * half_width
    if c < half_width:
        return (-half_width, -c, z)
    if c < 3.0 * half_width:
        return (-half_width + c - half_width, -half_width, z)
    if c < 5.0 * half_width:
        return (half_width, -half_width + c - 3.0 * half_width, z)
    if c < 7.0 * half_width:
        return (half_width - (c - 5.0 * half_width), half_width, z)
    return (-half_width, half_width - (c - 7.0 * half_width), z)


def arm_crossings(half_width, angle):
    """Points where the two arms (along +x and at `angle` counterclockwise) meet the square."""
    points = []
    for direction in ((1.0, 0.0), (math.cos(angle), math.sin(angle))):
        scale = half_width / max(abs(direction[0]), abs(direction[1]))
        p = [scale * direction[0], scale * direction[1]]
        for d in range(2):
            if abs(abs(p[d]) - half_width) <= 1.0e-12 * half_width:
                p[d] = math.copysign(half_width, p[d])
        points.append(tuple(p))
    return points


def on_corner_footprint(point, angle, convex, tolerance):
    first = point[1]
    second = point[0] * math.sin(angle) - point[1] * math.cos(angle)
    in_wedge = first >= -tolerance and second >= -tolerance
    if convex:
        return in_wedge
    on_arm = (abs(first) <= tolerance and second >= 0.0) or (abs(second) <= tolerance and first >= 0.0)
    return (not in_wedge) or on_arm


# ---------------------------------------------------------------------------- corner gate

def corner_geometry(model, points):
    angle = math.radians(float(model.get("AngleDegrees", model.get("Angle"))))
    convex = model["Topology"] == "ConvexCorner"
    groups = model.get("ContourGroups")
    zero = set(int(i) - 1 for i in model.get("ZeroTraceIndices", []))
    findings = []
    if not groups:
        return ["NotEvaluable: corner model without ContourGroups"], None
    offset = 0
    metal_rings = 0
    for size in groups:
        ring = points[offset:offset + size]
        half_width = max(max(abs(p[0]), abs(p[1])) for p in ring)
        tolerance = 1.0e-6 * half_width
        # A box ring of the coupon frame: every knot on the perimeter of the centred square.
        centred = all(abs(max(abs(p[0]), abs(p[1])) - half_width) <= tolerance for p in ring)
        metal = any((offset + i) in zero for i in range(size))
        if metal and centred:
            metal_rings += 1
            for crossing in arm_crossings(half_width, angle):
                hit = any((offset + i) in zero and math.hypot(ring[i][0] - crossing[0], ring[i][1] - crossing[1]) <= tolerance
                          for i in range(size))
                if not hit:
                    findings.append(
                        f"ring at z = {ring[0][2]:.4g}: metal arm crosses the box at "
                        f"({crossing[0]:.4g}, {crossing[1]:.4g}) with no PEC knot there "
                        "(a free hat has support on the PEC contour)")
            for i in range(size):
                if (offset + i) not in zero and on_corner_footprint(ring[i], angle, convex, tolerance):
                    findings.append(f"free knot {offset + i + 1} lies on the metal footprint")
        offset += size
    if metal_rings == 0:
        findings.append("NotEvaluable: no centred box ring with PEC knots")
    return findings, zero


def corner_diagonal(library_path, library, models, points_of):
    """Diagonal sanity check per convexity family: max free-knot fabricated Q diagonal per
    interface type vs the reference node (90 degrees, else the least turn)."""
    results = {}
    for convexity in CORNER_TOPOLOGIES:
        family = [m for m in models if m["Topology"] == convexity and float(m.get("CornerRadius", 0.0)) == 0.0
                  and m.get("FabricatedSurfaceMatrix")]
        if len(family) < 2:
            continue
        def angle_of(m):
            return float(m.get("AngleDegrees", m.get("Angle")))
        reference = min(family, key=lambda m: (abs(angle_of(m) - 90.0), angle_of(m)))
        maxima = {}
        for m in family:
            path = resolve(library_path, library, m["FabricatedSurfaceMatrix"])
            size = len(points_of[m["Name"]])
            zero = set(int(i) for i in m.get("ZeroTraceIndices", []))
            interfaces = {int(i["Coupon"]): i["Type"] for i in m.get("Interfaces", [])}
            diagonal = read_surface_diagonal(path, size)
            per_type = {}
            for interface, values in diagonal.items():
                free = [v for k, v in values.items() if k not in zero]
                if free:
                    per_type[interfaces.get(interface, str(interface))] = max(free)
            maxima[m["Name"]] = per_type
        for m in family:
            factors = {}
            for interface_type, value in maxima[m["Name"]].items():
                ref = maxima[reference["Name"]].get(interface_type)
                if ref and ref > 0.0:
                    factors[interface_type] = value / ref
            results[m["Name"]] = {"Reference": reference["Name"], "Factors": factors,
                                  "Passed": all(f <= DIAGONAL_FACTOR for t, f in factors.items() if t in ("MS", "MA"))}
    return results


# ---------------------------------------------------------------------------- spatial gate

def point_in_polygon(x, y, polygon):
    inside = False
    n = len(polygon)
    for i in range(n):
        (x0, y0), (x1, y1) = polygon[i], polygon[(i + 1) % n]
        if (y0 > y) != (y1 > y):
            xi = x0 + (y - y0) * (x1 - x0) / (y1 - y0)
            if x < xi:
                inside = not inside
    return inside


def polygons_from_segments(segments, scale):
    """Closed loops from an unordered segment soup."""
    adjacency = {}
    for (a, b) in segments:
        a = (round(a[0] * scale, 9), round(a[1] * scale, 9))
        b = (round(b[0] * scale, 9), round(b[1] * scale, 9))
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)
    seen = set()
    loops = []
    for start in adjacency:
        if start in seen:
            continue
        loop = [start]
        seen.add(start)
        previous, current = None, start
        while True:
            nxt = [n for n in adjacency[current] if n != previous]
            if not nxt:
                break
            nxt = nxt[0]
            if nxt == start:
                break
            if nxt in seen:
                break
            loop.append(nxt)
            seen.add(nxt)
            previous, current = current, nxt
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def segment_perimeter_crossings(polygon, xmin, xmax, ymin, ymax, tolerance):
    """Points where the polygon boundary crosses the rectangle perimeter."""
    crossings = []
    n = len(polygon)
    for i in range(n):
        (x0, y0), (x1, y1) = polygon[i], polygon[(i + 1) % n]
        dx, dy = x1 - x0, y1 - y0
        for xs in (xmin, xmax):
            if abs(dx) > tolerance:
                t = (xs - x0) / dx
                if -tolerance <= t <= 1.0 + tolerance:
                    y = y0 + t * dy
                    if ymin - tolerance <= y <= ymax + tolerance:
                        crossings.append((xs, min(max(y, ymin), ymax)))
        for ys in (ymin, ymax):
            if abs(dy) > tolerance:
                t = (ys - y0) / dy
                if -tolerance <= t <= 1.0 + tolerance:
                    x = x0 + t * dx
                    if xmin - tolerance <= x <= xmax + tolerance:
                        crossings.append((min(max(x, xmin), xmax), ys))
    return crossings


def spatial_geometry(library_path, library, model, points):
    findings = []
    trace = model.get("TraceMesh")
    boundary = model.get("PlanViewBoundary")
    if not trace or not boundary:
        return ["NotEvaluable: spatial model without TraceMesh / PlanViewBoundary"]
    radius = float(library["MatchingRadius"])
    vertices = read_trace_vertices(resolve(library_path, library, trace["Vertices"]))
    zero = set(int(i) for i in model.get("ZeroTraceIndices", []))
    polygons = []
    for conductor in boundary:
        for loop in polygons_from_segments([(s[0], s[1]) for s in conductor["Segments"]], PLAN_VIEW_SCALE * radius):
            polygons.append(loop)
    if not polygons:
        return ["NotEvaluable: PlanViewBoundary has no closed loop"]
    xs = [v["point"][0] for v in vertices]
    ys = [v["point"][1] for v in vertices]
    xmin, xmax, ymin, ymax = min(xs), max(xs), min(ys), max(ys)
    tolerance = 1.0e-6 * radius
    levels = sorted({round(v["point"][2], 9) for v in vertices})
    metal_levels = 0
    for level in levels:
        ring = [v for v in vertices if abs(v["point"][2] - level) <= tolerance and (
            abs(v["point"][0] - xmin) <= tolerance or abs(v["point"][0] - xmax) <= tolerance or
            abs(v["point"][1] - ymin) <= tolerance or abs(v["point"][1] - ymax) <= tolerance)]
        constrained = [v for v in ring if v["conductor"] > 0 or v["basis"] in zero]
        if not constrained:
            continue
        metal_levels += 1
        for polygon in polygons:
            for crossing in segment_perimeter_crossings(polygon, xmin, xmax, ymin, ymax, tolerance):
                if not any(math.hypot(v["point"][0] - crossing[0], v["point"][1] - crossing[1]) <= 10.0 * tolerance
                           for v in constrained):
                    findings.append(
                        f"level z = {level:.4g}: metal boundary crosses the box at "
                        f"({crossing[0]:.4g}, {crossing[1]:.4g}) with no conductor / PEC vertex")
        for v in ring:
            if v["basis"] > 0 and v["basis"] not in zero and v["conductor"] == 0:
                if any(point_in_polygon(v["point"][0], v["point"][1], polygon) for polygon in polygons):
                    findings.append(f"level z = {level:.4g}: free knot {v['basis']} lies inside the metal")
    if metal_levels == 0:
        findings.append("NotEvaluable: no trace ring level with conductor / PEC vertices")
    return findings


# ---------------------------------------------------------------------------- 2D gate

def translational_geometry(library, model, points):
    thickness = float(library.get("Fabrication", {}).get("MetalThickness", 0.0))
    if thickness <= 0.0:
        return ["NotEvaluable: library without Fabrication.MetalThickness"]
    radius = float(library["MatchingRadius"])
    tolerance = 1.0e-6 * radius
    topology = model["Topology"]
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    xmin, xmax = min(xs), max(xs)
    if topology in STRIP_2D:
        separation = float(model.get("Separation", 0.0))
        if separation < 2.0 * radius:
            return ["ByConstruction: strip metal (width Separation < 2R) does not reach the box sides"]
        return ["NotEvaluable: strip wider than the box"]
    if topology in CLUSTER_2D:
        return ["NotEvaluable: parallel-edge cluster metal pattern not derived from the record"]
    if topology in EDGE_2D:
        sides = [xmin]  # canonical frame: metal at x < 0
    elif topology in GAP_2D:
        sides = [xmin, xmax]
    else:
        return [f"NotEvaluable: two-dimensional topology {topology}"]
    findings = []
    for side in sides:
        band = [p for p in points if abs(p[0] - side) <= tolerance and -tolerance <= p[1] <= thickness + tolerance]
        if band:
            findings.append(f"free knot(s) on the metal band x = {side:.4g}, 0 <= y <= {thickness:g}: "
                            f"{[(round(p[0], 4), round(p[1], 4)) for p in band]}")
        column = sorted(p[1] for p in points if abs(p[0] - side) <= tolerance)
        below = [y for y in column if y < -tolerance]
        above = [y for y in column if y > thickness + tolerance]
        if not below or not above:
            findings.append(f"side x = {side:.4g}: no free knots on both sides of the metal band")
    return findings


# ---------------------------------------------------------------------------- driver

def evaluate(library_path):
    library = json.load(open(library_path))
    models = library.get("Models", [])
    points_of = {}
    record = {"Gate": "TraceBasis", "Library": os.path.abspath(library_path), "Name": library.get("Name"),
              "Rule": "no free trace basis function with support on the PEC part of the box contour "
                      "(corner-family review 2026-09-29); NotEvaluable fails the verdict",
              "DiagonalFactor": DIAGONAL_FACTOR, "Models": [], "Verdict": "Passed"}
    for model in models:
        entry = {"Name": model["Name"], "Topology": model["Topology"], "Findings": [], "Class": None}
        points_path = resolve(library_path, library, model.get("BasisPoints"))
        if not points_path:
            entry["Findings"].append("NotEvaluable: BasisPoints file missing")
        else:
            points = read_points(points_path)
            points_of[model["Name"]] = points
            topology = model["Topology"]
            if topology in CORNER_TOPOLOGIES and model.get("ContourGroups"):
                entry["Class"] = "corner"
                entry["Findings"], _ = corner_geometry(model, points)
            elif topology in SPATIAL_TOPOLOGIES:
                entry["Class"] = "spatial"
                entry["Findings"] = spatial_geometry(library_path, library, model, points)
            else:
                entry["Class"] = "translational"
                entry["Findings"] = translational_geometry(library, model, points)
        entry["Passed"] = not any(not f.startswith("ByConstruction") for f in entry["Findings"])
        record["Models"].append(entry)
    diagonal = corner_diagonal(library_path, library, models, points_of)
    for entry in record["Models"]:
        if entry["Name"] in diagonal:
            entry["Diagonal"] = diagonal[entry["Name"]]
            if not diagonal[entry["Name"]]["Passed"]:
                entry["Passed"] = False
                entry["Findings"].append(
                    "free-knot fabricated Q diagonal above DiagonalFactor x the reference node: "
                    f"{ {k: round(v, 3) for k, v in diagonal[entry['Name']]['Factors'].items()} }")
    if any(not e["Passed"] for e in record["Models"]):
        record["Verdict"] = "Failed"
    return record


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("library")
    parser.add_argument("--json")
    args = parser.parse_args()
    record = evaluate(args.library)
    for entry in record["Models"]:
        status = "PASS" if entry["Passed"] else "FAIL"
        extra = ""
        if "Diagonal" in entry:
            extra = " diagonal vs %s: %s" % (entry["Diagonal"]["Reference"],
                                             {k: "%.3g" % v for k, v in entry["Diagonal"]["Factors"].items()})
        print(f"{status} {entry['Name']} [{entry['Class']}]{extra}")
        for finding in entry["Findings"]:
            print(f"     - {finding}")
    print(f"VERDICT {record['Verdict']} ({record['Name']})")
    if args.json:
        with open(args.json, "w") as out:
            json.dump(record, out, indent=1)
    return 0 if record["Verdict"] == "Passed" else 1


if __name__ == "__main__":
    sys.exit(main())
