#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Pre-build element-count estimate of a Gmsh-only production case (supervisor decision
45(b) follow-up, P1 headroom gate): the recipe's size laws integrated over the coupon
before anything is built, so a case that cannot fit the element cap fails closed in
seconds instead of after a build.

The estimate is the size-field integral N = (1 / TetrahedraPerCubicSize) x integral of
dV / h(x)^3 over the coupon, with h the composed production size field, evaluated
component by component from the frozen inputs only (signature, boundary, process, mesh
recipe, trace basis, build options):

- far field: the coupon box volume at FarSize;
- prism tubes: every straight metal edge not on the box carries a top and a bottom tube,
  Layers = length / TangentialSize (bounded by FarSize), prisms per layer = Sectors +
  2 x Sectors x (Rings - 1), pyramids per layer = Sectors (the recipe's ring rule and
  sector angle); the tube band (NormalSize at the tube surface, FarGrowth beyond) is a
  cylindrical shell integral per unit length on the dielectric side;
- trace basis: every basis triangle with requested size s = TraceBasisSizeRatio x
  minimum altitude below FarSize contributes the half-space Steiner shell integral of
  the volume law (s + FarGrowth x distance) over its area and perimeter - the needle
  cost the gate exists for (a needle's altitude sets s, its length the perimeter);
- corner balls and tube caps: the graded shells from CornerSize to NormalSize inside
  CornerIsotropyRadius, then the corner exterior law to FarSize, one ball per semantic
  corner and per tube cap (a cap centre is a graded point of the same law);
- junction / band curves: the cut-surface junction lines of the etch trench are a
  producer outcome (the footprint is not derivable from the inputs), so their band-law
  integral is charged per unit length on a recorded proxy length: the metal-loop
  perimeter inside the box; the census JunctionCurves length is recorded next to the
  proxy after every build.

TetrahedraPerCubicSize is the one dimensionless model constant: the number of
tetrahedra Gmsh realizes per h^3 of prescribed volume (an equilateral tetrahedron of
edge h has volume h^3 / (6 sqrt 2) = 0.1179 h^3; Delaunay meshes of a graded field
realize a mean edge below the prescription, so the constant is calibrated).  It is
recorded in the manifest (ProductionRecipe.BuildCostEstimate) together with the
verified roots it was calibrated on and their estimate / actual ratios; it changes no
mesh and no gate: the gate is the unchanged MaximumElements, applied to the estimate
before the build.

usage: estimate_build_cost.py CASE_ID [--manifest PATH] [--actual BUILD_CENSUS]
"""
import argparse
import csv
import json
import math
from pathlib import Path
import tomllib

import numpy as np

HERE = Path(__file__).resolve().parent
TUBE_SECTOR_SPAN_DEGREES = 270.0   # the dielectric side of a metal edge (3 / 4 turn)
DIELECTRIC_FRACTION = 0.75         # of a ball / cylinder around a metal edge or corner
TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING = 0.5
BAND_RADIAL_GROWTH = 1.0
BAND_PROTECTED_DISTANCE_OVER_NORMAL = 2.0
QUADRATURE_POINTS = 4000


def read_edges(signature):
    with Path(signature).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    edges = []
    for row in rows:
        edges.append({"slot": int(row["Slot"]), "conductor": int(row["Conductor"]),
                      "point": np.array([float(row[k]) for k in ("Px", "Py", "Pz")]),
                      "gap": np.array([float(row[k]) for k in ("Gx", "Gy", "Gz")]),
                      "tangent": np.array([float(row[k]) for k in ("Tx", "Ty", "Tz")]),
                      "interval": (float(row["S0"]), float(row["S1"])),
                      "normal_sign": int(float(row["Nz"])),
                      "vertex_arm": bool(int(float(row.get("VertexArm", 0) or 0)))})
    return edges


def edge_chains(edges):
    """Decision 47: collinear touching rows of one metal edge (same slot / conductor are
    not in the estimator's rows, so same plane, tangent, gap and no vertex arm) form a
    chain; returns {row index: (union start, union end) in the row's own coordinate}."""
    tolerance = 1e-9 * max(1.0, max(float(np.abs(edge["point"]).max()) for edge in edges))
    ends = [[edge["point"] + s * edge["tangent"] for s in edge["interval"]] for edge in edges]

    def same_line(i, j):
        a, b = edges[i], edges[j]
        if a["vertex_arm"] or b["vertex_arm"] or a.get("slot") != b.get("slot") or a.get("conductor") != b.get("conductor"):
            return False
        if (np.abs(a["tangent"] - b["tangent"]).max() > tolerance or np.abs(a["gap"] - b["gap"]).max() > tolerance or
                abs(a["point"][2] - b["point"][2]) > tolerance):
            return False
        delta = b["point"] - a["point"]
        return np.linalg.norm(delta - np.dot(delta, a["tangent"]) * a["tangent"]) <= tolerance

    unions, used = {}, set()
    for i in range(len(edges)):
        if i in used or edges[i]["vertex_arm"]:
            continue
        chain, used, changed = [i], used | {i}, True
        while changed:
            changed = False
            for j in range(len(edges)):
                if j in used or not same_line(i, j):
                    continue
                if any(np.linalg.norm(a - b) <= tolerance for k in chain for a in ends[k] for b in ends[j]):
                    chain.append(j); used.add(j); changed = True
        if len(chain) < 2:
            continue
        origin, tangent = edges[chain[0]]["point"], edges[chain[0]]["tangent"]
        coordinate = lambda k, s: float(np.dot(edges[k]["point"] - origin, tangent)) + s
        u0 = min(coordinate(k, edges[k]["interval"][0]) for k in chain)
        u1 = max(coordinate(k, edges[k]["interval"][1]) for k in chain)
        for k in chain:
            shift = coordinate(k, 0.0)
            unions[k] = (u0 - shift, u1 - shift)
    return unions


def coupon_box(edges, radius, metal_thickness, overetch):
    """The mesher's coupon box (mesh_spatial_coupon.jl coupon_bounds / extended_interval,
    with the decision-47 chain rule for CAD-subdivided edges and the per-layer-sign
    vertical padding of decision 48)."""
    points = []
    unions = edge_chains(edges)
    for index, edge in enumerate(edges):
        first, second = edge["interval"]
        extension = 2.0 * radius
        tolerance = 1e-10 * radius
        if index in unions:
            u0, u1 = unions[index]
            if (u1 - u0) / 2.0 >= radius - tolerance:
                if abs(first - u0) <= tolerance:
                    first -= extension
                if abs(second - u1) <= tolerance:
                    second += extension
        elif edge["vertex_arm"]:
            if abs(first) <= tolerance:
                second += extension
            elif abs(second) <= tolerance:
                first -= extension
        else:
            if first <= -radius + tolerance:
                first -= extension
            if second >= radius - tolerance:
                second += extension
        for s in (first, second):
            boundary = edge["point"] + s * edge["tangent"]
            for side in (-1.0, 1.0):
                points.append(boundary + side * radius * edge["gap"])
    points = np.asarray(points)
    lower = points.min(axis=0) - radius
    upper = points.max(axis=0) + radius
    # Vertical padding per process layer sign (decision 48): Overetch on the substrate
    # side (-Nz), MetalThickness on the metal side (+Nz) of every row's plane.
    lower[2] = min(lower[2], min(edge["point"][2] - radius - (overetch if edge["normal_sign"] > 0 else metal_thickness)
                                 for edge in edges))
    upper[2] = max(upper[2], max(edge["point"][2] + radius + (metal_thickness if edge["normal_sign"] > 0 else overetch)
                                 for edge in edges))
    return lower, upper


def read_loops(boundary):
    with Path(boundary).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    loops = {}
    for row in rows:
        loops.setdefault(int(row["Loop"]), []).append(
            ((float(row["X"]), float(row["Y"])), row["Class"], float(row["Plane"])))
    return list(loops.values())


def on_box(point, lower, upper, tolerance):
    return any(abs(point[d] - lower[d]) <= tolerance or abs(point[d] - upper[d]) <= tolerance
               for d in range(2))


def metal_sides(loops, lower, upper, tolerance):
    """Straight metal loop sides that are not box faces (each carries two tubes) and the
    number of tube caps: two per side end at a semantic (Physical) corner."""
    sides, caps = [], 0
    for loop in loops:
        n = len(loop)
        for i in range(n):
            (p, class_p, _), (q, class_q, _) = loop[i], loop[(i + 1) % n]
            same_face = any((abs(p[d] - lower[d]) <= tolerance and abs(q[d] - lower[d]) <= tolerance) or
                            (abs(p[d] - upper[d]) <= tolerance and abs(q[d] - upper[d]) <= tolerance)
                            for d in range(2))
            if same_face:
                continue
            sides.append(math.dist(p, q))
            caps += 2 * sum(1 for point, cls in ((p, class_p), (q, class_q))
                            if cls == "Physical" and not on_box(point, lower, upper, tolerance))
    return sides, caps


def tube_rings(edge_size, ratio, overetch, metal_thickness, corner_radius):
    """Largest K with r_K + h_K <= min(Overetch, MetalThickness / 2, CornerIsotropyRadius)."""
    bound = min(overetch, metal_thickness / 2.0, corner_radius)
    sizes, radius = [], 0.0
    size = edge_size
    while radius + size + size <= bound or not sizes:
        if radius + size > bound:
            break
        sizes.append(size)
        radius += size
        size *= ratio
    return sizes, radius


def band_law(r, normal, far, growth):
    protected = BAND_PROTECTED_DISTANCE_OVER_NORMAL * normal
    return np.minimum(far, normal + BAND_RADIAL_GROWTH * np.minimum(r, protected) +
                      growth * np.maximum(r - protected, 0.0))


def shell_integral(r0, r1, h, weight):
    """Integral of weight(r) / h(r)^3 over [r0, r1] by composite midpoint quadrature."""
    if r1 <= r0:
        return 0.0
    r = r0 + (np.arange(QUADRATURE_POINTS) + 0.5) * (r1 - r0) / QUADRATURE_POINTS
    return float(np.sum(weight(r) / h(r) ** 3) * (r1 - r0) / QUADRATURE_POINTS)


def corner_ball_law(r, corner_size, ratio, normal, corner_radius, far, growth):
    """Graded shells from CornerSize (x ratio per shell) to NormalSize inside the ball,
    then the corner exterior law NormalSize + FarGrowth x (r - CornerIsotropyRadius)."""
    sizes, radii, size, radius = [], [], corner_size, 0.0
    while size < normal:
        radius += size
        sizes.append(size); radii.append(radius)
        size *= ratio
    h = np.full_like(r, normal, dtype=float)
    for shell_radius, shell_size in zip(reversed(radii), reversed(sizes)):
        h = np.where(r < shell_radius, shell_size, h)
    outside = r > corner_radius
    h = np.where(outside, np.minimum(far, normal + growth * (r - corner_radius)), h)
    return h


def trace_basis_integral(vertices_path, triangles_path, ratio, far, growth):
    vertices = {}
    with Path(vertices_path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            vertices[int(row["vertex"])] = np.array([float(row[k]) for k in ("x", "y", "z")])
    integral, narrow, minimum_altitude, areas_narrow = 0.0, 0, math.inf, 0.0
    with Path(triangles_path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        a, b, c = (vertices[int(row[k])] for k in ("vertex_i", "vertex_j", "vertex_k"))
        area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
        longest = max(np.linalg.norm(b - a), np.linalg.norm(c - b), np.linalg.norm(a - c))
        altitude = 2.0 * area / longest
        minimum_altitude = min(minimum_altitude, altitude)
        s = ratio * altitude
        if s < far:
            narrow += 1
            areas_narrow += area
            perimeter = np.linalg.norm(b - a) + np.linalg.norm(c - b) + np.linalg.norm(a - c)
            # Half-space Steiner volume at distance d from the triangle: A d + (pi / 4) P d^2,
            # so the shell weight is A + (pi / 2) P d under the law s + FarGrowth x d.
            integral += shell_integral(0.0, (far - s) / growth, lambda d: s + growth * d,
                                       lambda d: area + 0.5 * math.pi * perimeter * d)
    return {"Triangles": len(rows), "NarrowTriangles": narrow, "NarrowArea": areas_narrow,
            "MinimumBasisAltitude": minimum_altitude, "MinimumRequestedSize": ratio * minimum_altitude,
            "Integral": integral}


def estimate(paths, options, tetrahedra_per_cubic_size):
    """Component integrals and the element estimate of one case from its frozen inputs."""
    process = tomllib.loads(Path(paths["Process"]).read_text())
    recipe = json.loads(Path(paths["MeshRecipe"]).read_text())
    radius, thickness, overetch = process["Radius"], process["MetalThickness"], process["Overetch"]
    normal = recipe["NormalSizeOverThickness"] * thickness
    corner_radius = normal * recipe["TangentialSizeOverNormalSize"]
    far = recipe["FarSizeOverRadius"] * radius
    requested_tangent = float(options["--lc-tangent"])
    tangent = min(requested_tangent, far)               # the coupon-scale size bound
    edge_size, growth_ratio = float(options["--edge-size"]), float(options["--edge-growth-ratio"])
    corner_size, far_growth = float(options["--corner-size"]), float(options["--far-growth"])
    edges = read_edges(paths["Signature"])
    lower, upper = coupon_box(edges, radius, thickness, overetch)
    tolerance = 1e-8 * radius
    extent = upper - lower
    volume = float(np.prod(extent))
    loops = read_loops(paths["Boundary"])
    sides, caps = metal_sides(loops, lower, upper, tolerance)
    corners = sum(1 for loop in loops for _, cls, _ in loop if cls == "Physical")
    ring_sizes, tube_radius = tube_rings(edge_size, growth_ratio, overetch, thickness, corner_radius)
    sectors = int(round(TUBE_SECTOR_SPAN_DEGREES / 30.0))
    rings = len(ring_sizes)
    tube_length = 2.0 * sum(sides)
    layers = sum(2 * math.ceil(side / tangent) for side in sides)
    prisms = layers * (sectors + 2 * sectors * (rings - 1))
    pyramids = layers * sectors
    offset = tube_radius + TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING * ring_sizes[-1]
    reach = offset + (far - normal) / far_growth
    tube_band = tube_length * shell_integral(
        offset, reach, lambda r: np.minimum(far, normal + far_growth * (r - offset)),
        lambda r: DIELECTRIC_FRACTION * 2.0 * math.pi * r)
    ball_reach = corner_radius + (far - normal) / far_growth
    ball = shell_integral(0.0, ball_reach,
                          lambda r: corner_ball_law(r, corner_size, growth_ratio, normal, corner_radius, far, far_growth),
                          lambda r: DIELECTRIC_FRACTION * 4.0 * math.pi * r ** 2)
    corner_balls = (corners + caps) * ball
    perimeter = sum(math.dist(loop[i][0], loop[(i + 1) % len(loop)][0])
                    for loop in loops for i in range(len(loop)))
    band_reach = BAND_PROTECTED_DISTANCE_OVER_NORMAL * normal + (far - 3.0 * normal) / far_growth
    band_per_length = shell_integral(0.0, band_reach, lambda r: band_law(r, normal, far, far_growth),
                                     lambda r: 2.0 * math.pi * r)
    junction = perimeter * band_per_length
    far_field = volume / far ** 3
    trace = None
    if "TraceVertices" in paths:
        trace = trace_basis_integral(paths["TraceVertices"], paths["TraceTriangles"],
                                     float(options["--trace-basis-size-ratio"]), far, far_growth)
    integrals = {"FarField": far_field, "TubeBand": tube_band, "CornerBallsAndCaps": corner_balls,
                 "JunctionProxy": junction, "TraceBasis": trace["Integral"] if trace else 0.0}
    total_integral = sum(integrals.values())
    tetrahedra = total_integral / tetrahedra_per_cubic_size
    return {"Rule": ("N = FarField + TubeBand + CornerBallsAndCaps + JunctionProxy + TraceBasis size-field "
                     "integrals of dV / h^3 divided by TetrahedraPerCubicSize, plus the tube prisms and "
                     "pyramids; compared with MaximumElements before the build (fail closed)"),
            "Box": {"Lower": lower.tolist(), "Upper": upper.tolist(), "Volume": volume},
            "Sizes": {"NormalSize": normal, "TangentialSize": tangent, "RequestedTangentialSize": requested_tangent,
                      "FarSize": far, "CornerIsotropyRadius": corner_radius, "EdgeSize": edge_size,
                      "CornerSize": corner_size, "GrowthRatio": growth_ratio, "FarGrowth": far_growth},
            "Tubes": {"Sides": len(sides), "Length": tube_length, "Layers": layers, "Rings": rings,
                      "Sectors": sectors, "TubeRadius": tube_radius, "Prisms": prisms, "Pyramids": pyramids},
            "Corners": corners, "Caps": caps, "LoopPerimeter": perimeter,
            "TraceBasis": trace,
            "Integrals": integrals, "TotalIntegral": total_integral,
            "TetrahedraPerCubicSize": tetrahedra_per_cubic_size,
            "EstimatedTetrahedra": tetrahedra, "EstimatedPrisms": prisms, "EstimatedPyramids": pyramids,
            "EstimatedElements": tetrahedra + prisms + pyramids}


def case_paths(manifest, manifest_path, case):
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    directory = Path(case["Source"]["Directory"])
    directory = directory if directory.is_absolute() else repository / directory
    return {role: (repository / entry["RepositoryPath"] if entry.get("RepositoryPath") else directory / entry["Name"])
            for role, entry in case["Source"]["Files"].items()}


def gate(manifest, manifest_path, case):
    """The estimate of a manifest case against the manifest's element cap; a case whose
    estimate exceeds the cap is reported Passed false (the caller fails closed)."""
    recipe = manifest["ProductionRecipe"]
    options = dict(recipe["BuildCommandOptions"])
    model = recipe["BuildCostEstimate"]
    result = estimate(case_paths(manifest, manifest_path, case), options, model["TetrahedraPerCubicSize"])
    cap = manifest["Gates"]["MaximumElements"]
    result.update({"MaximumElements": cap, "EstimateOverCap": result["EstimatedElements"] / cap,
                   "Passed": result["EstimatedElements"] <= cap})
    return result


def actual_counts(census_path):
    census = json.loads(Path(census_path).read_text())
    quality = census["PrismTubes"]["Quality"]
    return {kind: quality[kind]["Count"] for kind in ("Tetrahedron", "Prism", "Pyramid")}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_id")
    parser.add_argument("--manifest", type=Path, default=HERE / "geometry-independence-suite.json")
    parser.add_argument("--actual", type=Path, help="a build census to compare the estimate with")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    case = next(item for item in manifest["Cases"] if item["Id"] == args.case_id)
    result = gate(manifest, manifest_path, case)
    if args.actual is not None:
        counts = actual_counts(args.actual)
        result["Actual"] = {**counts, "Elements": sum(counts.values()),
                            "EstimateOverActualTetrahedra": result["EstimatedTetrahedra"] / counts["Tetrahedron"],
                            "EstimateOverActualElements": result["EstimatedElements"] / sum(counts.values())}
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    shown = {k: result[k] for k in ("EstimatedTetrahedra", "EstimatedPrisms", "EstimatedPyramids",
                                    "EstimatedElements", "EstimateOverCap", "Passed", "Integrals")}
    if "Actual" in result:
        shown["Actual"] = result["Actual"]
    print(json.dumps(shown, indent=1))
    raise SystemExit(0 if result["Passed"] else 1)


if __name__ == "__main__":
    main()
