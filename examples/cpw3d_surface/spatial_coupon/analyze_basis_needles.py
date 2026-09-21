#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Needle cost of a coupon's trace basis (supervisor decision 54b, analysis only).

Decomposes the pre-build element estimate of a coupon source directory
(estimate_build_cost.estimate under the production manifest's recipe and options) into
its components - tube prisms and pyramids, tube band, corner balls and caps, junction
band proxy, far volume and the trace-basis volume law - and splits the trace-basis law
by basis triangle class: a NEEDLE has minimum altitude below
NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE x its shortest edge (trace_basis.py; the fan
triangles of the box's top and bottom faces), a RIGHT SLIVER is any other triangle
whose requested size TraceBasisSizeRatio x altitude lies below FarSize (the strip
triangles of the box side faces), the rest is WIDE (no trace cost).  The needle share is
the needle tetrahedra over the estimated elements.

Given a built root (--census, --mesh) it also attributes every tetrahedron of the actual
mesh to the size law governing at its centroid (the smallest of the tube band law, the
corner / cap ball law, the junction band law, the trace-basis law per narrow triangle
and FarSize), so the estimate's decomposition can be read against the built mesh.

Then it evaluates two basis constructions for DEVICE coupons (whose basis is the
library's own, generate_spatial_response --basis-only), never the gallery references:

(a) DELAUNAY: re-triangulate every box face whose basis vertices all lie on the face
    boundary (the ear-clipped caps) by Lawson edge flips to the Delaunay (max-min-angle)
    triangulation of the same vertices (generate_spatial_response.delaunay_flip_cap, the
    producer's --cap-triangulation delaunay) - the sources are unchanged, the needles
    become right slivers whose altitude is the local ring spacing, and the estimate is
    recomputed exactly under the unchanged recipe.  Face selection differs from the
    producer's: the producer flips the two caps only, this tool flips every box face
    whose vertices all lie on the face boundary - identical for the three-level boxes
    of every device and gallery coupon (the side faces carry mid-level vertices), not
    for a two-level box whose side faces would be flipped here and not by the producer
    (the recorded exact agreement, 7f03 3,135,224, is a three-level box);
(b) GRADED REFINEMENT (modelled, not constructed): a basis whose local edge h(x) on
    every face is the ring spacing graded away from the ring vertices at slope
    BasisGrowth (h = min_i spacing_i + BasisGrowth x |x - v_i|, capped at the face
    diagonal); the added sources are the interior vertex count
    2 / (sqrt 3 h^2) integrated over the faces, and the trace-basis cost is the
    half-space shell integral of the volume law with s = TraceBasisSizeRatio x
    (sqrt 3 / 2) h (the altitude of an equilateral triangle of edge h) integrated over
    the faces.  The physics cost index is Sources x Elements (per-source solve time
    scales with the mesh), reported relative to the current basis.

usage: analyze_basis_needles.py SOURCE_DIR [--manifest PATH] [--census BUILD_CENSUS --mesh MESH]
                                 [--basis-growth G ...] [--output JSON]
"""
import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import estimate_build_cost as cost
from generate_spatial_response import delaunay_flip_cap
from trace_basis import NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE

HERE = Path(__file__).resolve().parent
SOURCE_ROLES = {"Signature": "mesh-signature.csv", "Boundary": "plan-view-boundary.csv",
                "Mask": "plan-view-mask.csv", "Process": "process.toml",
                "TraceVertices": "trace-vertices.csv", "TraceTriangles": "trace-triangles.csv"}
EQUILATERAL_ALTITUDE_OVER_EDGE = math.sqrt(3.0) / 2.0
VERTICES_PER_EDGE_SQUARED = 2.0 / math.sqrt(3.0)     # vertex density of an equilateral mesh of edge h
FACE_QUADRATURE = 400                                  # points per face side for the graded-basis model


def read_basis(directory):
    vertices = {}
    with (Path(directory) / "trace-vertices.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            vertices[int(row["vertex"])] = np.array([float(row[k]) for k in ("x", "y", "z")])
    with (Path(directory) / "trace-triangles.csv").open(newline="") as stream:
        triangles = [tuple(int(row[k]) for k in ("vertex_i", "vertex_j", "vertex_k"))
                     for row in csv.DictReader(stream)]
    return vertices, triangles


def triangle_geometry(a, b, c):
    lengths = (np.linalg.norm(b - a), np.linalg.norm(c - b), np.linalg.norm(a - c))
    area = 0.5 * np.linalg.norm(np.cross(b - a, c - a))
    altitude = 2.0 * area / max(lengths) if max(lengths) > 0.0 else 0.0
    return area, lengths, altitude


def classify(vertices, triangles, ratio, far):
    """Per triangle: class, requested size, area, perimeter, altitude."""
    rows = []
    for triangle in triangles:
        a, b, c = (vertices[v] for v in triangle)
        area, lengths, altitude = triangle_geometry(a, b, c)
        s = ratio * altitude
        if s >= far:
            kind = "Wide"
        elif altitude < NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE * min(lengths):
            kind = "Needle"
        else:
            kind = "RightSliver"
        rows.append({"Triangle": triangle, "Class": kind, "RequestedSize": s, "Area": area,
                     "Perimeter": float(sum(lengths)), "Altitude": altitude,
                     "ShortestEdge": float(min(lengths)), "LongestEdge": float(max(lengths))})
    return rows


def trace_integral_by_class(rows, far, growth):
    integrals = {"Needle": 0.0, "RightSliver": 0.0, "Wide": 0.0}
    for row in rows:
        if row["Class"] == "Wide":
            continue
        s, area, perimeter = row["RequestedSize"], row["Area"], row["Perimeter"]
        integrals[row["Class"]] += cost.shell_integral(
            0.0, (far - s) / growth, lambda d: s + growth * d,
            lambda d: area + 0.5 * math.pi * perimeter * d)
    return integrals


def source_paths(directory):
    directory = Path(directory)
    paths = {role: directory / name for role, name in SOURCE_ROLES.items()}
    paths["MeshRecipe"] = HERE / "testdata" / "generality-mesh-recipe.json"
    return paths


def estimate_decomposition(paths, options, model, rows, far, growth):
    result = cost.estimate(paths, options, model)
    per_cubic = result["TetrahedraPerCubicSize"]
    by_class = trace_integral_by_class(rows, far, growth)
    components = {name: value / per_cubic for name, value in result["Integrals"].items()}
    components["TraceBasisNeedles"] = by_class["Needle"] / per_cubic
    components["TraceBasisRightSlivers"] = by_class["RightSliver"] / per_cubic
    components["TubePrisms"] = result["EstimatedPrisms"]
    components["TubePyramids"] = result["EstimatedPyramids"]
    elements = result["EstimatedElements"]
    return {"EstimatedElements": elements, "EstimatedTetrahedra": result["EstimatedTetrahedra"],
            "Components": components,
            "Shares": {name: value / elements for name, value in components.items()},
            "NeedleShare": components["TraceBasisNeedles"] / elements,
            "Classes": {kind: sum(1 for row in rows if row["Class"] == kind)
                        for kind in ("Needle", "RightSliver", "Wide")},
            "MinimumAltitude": {kind: min([row["Altitude"] for row in rows if row["Class"] == kind], default=None)
                                for kind in ("Needle", "RightSliver")},
            "Box": result["Box"], "Sizes": result["Sizes"]}


# --- (a) Delaunay re-triangulation of the boundary-only box faces -----------------------

def box_faces(vertices):
    points = np.array(list(vertices.values()))
    lower, upper = points.min(axis=0), points.max(axis=0)
    tolerance = 1e-9 * float(np.abs(points).max())
    faces = []
    for axis in range(3):
        for value in (lower[axis], upper[axis]):
            faces.append((axis, value, tolerance))
    return faces, lower, upper


def face_of(triangle, vertices, faces):
    for index, (axis, value, tolerance) in enumerate(faces):
        if all(abs(vertices[v][axis] - value) <= tolerance for v in triangle):
            return index
    return None


def delaunay_faces(vertices, triangles):
    """Re-triangulate the faces whose vertices all lie on the face boundary; returns the
    new triangle list and the list of re-triangulated face indices.  The producer flips
    the two caps only: on a three-level box the selections coincide (the side faces have
    mid-level vertices), on a two-level box this tool would also flip the side faces."""
    faces, lower, upper = box_faces(vertices)
    grouped = {}
    for triangle in triangles:
        grouped.setdefault(face_of(triangle, vertices, faces), []).append(triangle)
    result = []
    retriangulated = []
    for face, members in grouped.items():
        if face is None:
            result.extend(members)
            continue
        axis, value, tolerance = faces[face]
        others = [d for d in range(3) if d != axis]
        used = sorted({v for t in members for v in t})
        on_boundary = all(any(abs(vertices[v][d] - lower[d]) <= tolerance or abs(vertices[v][d] - upper[d]) <= tolerance
                              for d in others) for v in used)
        if not on_boundary:
            result.extend(members)
            continue
        projected = np.zeros((max(vertices) + 1, 3))
        for v in used:
            projected[v, :2] = (vertices[v][others[0]], vertices[v][others[1]])
        block = [tuple(t) for t in members]
        delaunay_flip_cap(block, projected, 0)
        result.extend(block)
        retriangulated.append(face)
    return result, retriangulated


# --- (b) graded-basis model -------------------------------------------------------------

def graded_basis_model(vertices, triangles, ratio, far, growth, basis_growth, per_cubic):
    """Sources and trace-basis tetrahedra of a basis whose face triangulations are graded
    from the ring spacing at slope basis_growth (see the module docstring)."""
    faces, lower, upper = box_faces(vertices)
    ring_vertices = sorted(vertices)
    # Spacing at a ring vertex: the shortest basis edge incident to it.
    incident = {v: math.inf for v in ring_vertices}
    for t in triangles:
        for k in range(3):
            a, b = t[k], t[(k + 1) % 3]
            length = float(np.linalg.norm(vertices[a] - vertices[b]))
            incident[a] = min(incident[a], length)
            incident[b] = min(incident[b], length)
    interior_sources = 0.0
    trace_integral = 0.0
    surface_shell = lambda s: (1.0 / (2.0 * growth)) * (1.0 / s ** 2 - 1.0 / far ** 2) if s < far else 0.0
    for axis, value, tolerance in faces:
        others = [d for d in range(3) if d != axis]
        face_vertices = [v for v in ring_vertices if abs(vertices[v][axis] - value) <= tolerance]
        if len(face_vertices) < 3:
            continue
        spacing = np.array([incident[v] for v in face_vertices])
        anchors = np.array([[vertices[v][d] for d in others] for v in face_vertices])
        cap = float(math.hypot(upper[others[0]] - lower[others[0]], upper[others[1]] - lower[others[1]]))
        u = np.linspace(lower[others[0]], upper[others[0]], FACE_QUADRATURE, endpoint=False)
        w = np.linspace(lower[others[1]], upper[others[1]], FACE_QUADRATURE, endpoint=False)
        du = (upper[others[0]] - lower[others[0]]) / FACE_QUADRATURE
        dw = (upper[others[1]] - lower[others[1]]) / FACE_QUADRATURE
        grid = np.stack(np.meshgrid(u + 0.5 * du, w + 0.5 * dw, indexing="ij"), axis=-1).reshape(-1, 2)
        h = np.full(len(grid), cap)
        for anchor, s0 in zip(anchors, spacing):
            h = np.minimum(h, s0 + basis_growth * np.linalg.norm(grid - anchor, axis=1))
        interior_sources += float(np.sum(VERTICES_PER_EDGE_SQUARED / h ** 2) * du * dw)
        s = ratio * EQUILATERAL_ALTITUDE_OVER_EDGE * h
        shell = np.where(s < far, (1.0 / (2.0 * growth)) * (1.0 / s ** 2 - 1.0 / far ** 2), 0.0)
        trace_integral += float(np.sum(shell) * du * dw)
    return {"BasisGrowth": basis_growth, "RingSources": len(ring_vertices),
            "ModelledSources": len(ring_vertices) + interior_sources,
            "TraceBasisTetrahedra": trace_integral / per_cubic}


# --- actual attribution of a built mesh -------------------------------------------------

def segment_distance(points, segments):
    """Distance from every point to the nearest of the segments (n x 6 array)."""
    best = np.full(len(points), np.inf)
    for segment in segments:
        a, b = np.asarray(segment[:3]), np.asarray(segment[3:])
        d = b - a
        t = np.clip(((points - a) @ d) / max(float(d @ d), 1e-300), 0.0, 1.0)
        best = np.minimum(best, np.linalg.norm(points - a - t[:, None] * d, axis=1))
    return best


def point_triangle_distance(points, a, b, c):
    ab, ac, ap = b - a, c - a, points - a
    n = np.cross(ab, ac)
    nn = float(n @ n)
    d00, d01, d11 = float(ab @ ab), float(ab @ ac), float(ac @ ac)
    d20, d21 = ap @ ab, ap @ ac
    denominator = d00 * d11 - d01 * d01
    v = (d11 * d20 - d01 * d21) / denominator
    w = (d00 * d21 - d01 * d20) / denominator
    inside = (v >= 0.0) & (w >= 0.0) & (v + w <= 1.0)
    plane = np.abs(ap @ n) / math.sqrt(nn)
    edge = np.full(len(points), np.inf)
    for first, second in ((a, b), (b, c), (c, a)):
        d = second - first
        t = np.clip(((points - first) @ d) / float(d @ d), 0.0, 1.0)
        edge = np.minimum(edge, np.linalg.norm(points - first - t[:, None] * d, axis=1))
    return np.where(inside, plane, edge)


def attribute_actual(census_path, mesh_path, ratio, chunk=50000):
    import meshio  # noqa: E402  (the repository's mesh readers depend on it)
    census = json.loads(Path(census_path).read_text())
    tubes = census["PrismTubes"]
    normal, far, growth = tubes["NormalSize"], tubes["FarSize"], tubes["FarGrowth"]
    corner_radius = census["CornerIsotropyRadius"]
    section = tubes["Section"]
    offset = section["Radius"] + section["PyramidHeight"]
    axes = [row["StartPoint"] + row["EndPoint"] for row in tubes["Tubes"]]
    balls = [c for c in census["SemanticCorners"]] + [r["Point"] for r in tubes["CapRegions"]["Regions"]]
    bands = census["JunctionCurves"]["Segments"] + tubes["BandCurves"]["Segments"]
    sizing = census["TraceBasisSizing"]
    narrow = []
    for triangle in sizing["MeshFrameTriangles"]:
        a, b, c = (np.array(p) for p in triangle)
        area, lengths, altitude = triangle_geometry(a, b, c)
        s = ratio * altitude
        if s < far:
            kind = "Needle" if altitude < NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE * min(lengths) else "RightSliver"
            narrow.append((a, b, c, s, kind))
    mesh = meshio.read(mesh_path)
    counts = {"Prism": 0, "Pyramid": 0}
    tets = None
    for block in mesh.cells:
        if block.type == "tetra":
            tets = block.data if tets is None else np.vstack([tets, block.data])
        elif block.type == "wedge":
            counts["Prism"] += len(block.data)
        elif block.type == "pyramid":
            counts["Pyramid"] += len(block.data)
    laws = ["TubeBand", "CornerBallsAndCaps", "JunctionBand", "TraceBasisNeedles", "TraceBasisRightSlivers", "Far"]
    attributed = {law: 0 for law in laws}
    corner_law = lambda r: cost.corner_ball_law(r, tubes["InnerSize"], tubes["GrowthRatio"], normal, corner_radius, far, growth)
    for start in range(0, len(tets), chunk):
        centroids = mesh.points[tets[start:start + chunk]].mean(axis=1)
        tube = np.minimum(far, normal + growth * np.maximum(segment_distance(centroids, axes) - offset, 0.0))
        ball = corner_law(np.min(np.stack([np.linalg.norm(centroids - np.asarray(b), axis=1) for b in balls]), axis=0))
        band = cost.band_law(segment_distance(centroids, bands), normal, far, growth)
        needle = np.full(len(centroids), far)
        sliver = np.full(len(centroids), far)
        for a, b, c, s, kind in narrow:
            value = np.minimum(far, s + growth * point_triangle_distance(centroids, a, b, c))
            if kind == "Needle":
                needle = np.minimum(needle, value)
            else:
                sliver = np.minimum(sliver, value)
        stacked = np.stack([tube, ball, band, needle, sliver, np.full(len(centroids), far)])
        governing = np.argmin(stacked, axis=0)
        # Where every law has saturated at FarSize the far field governs (ties go to it).
        governing[stacked.min(axis=0) >= far * (1.0 - 1e-9)] = len(laws) - 1
        for index, law in enumerate(laws):
            attributed[law] += int(np.sum(governing == index))
    total = len(tets) + counts["Prism"] + counts["Pyramid"]
    return {"Tetrahedra": len(tets), "TubePrisms": counts["Prism"], "TubePyramids": counts["Pyramid"],
            "Elements": total, "Attributed": attributed,
            "Shares": {law: value / total for law, value in attributed.items()},
            "NeedleShare": attributed["TraceBasisNeedles"] / total,
            "Rule": "each tetrahedron is attributed to the size law achieving the smallest size at its "
                    "centroid among the tube band law (NormalSize at the tube surface + FarGrowth x "
                    "distance), the corner / cap ball law, the junction band law on the junction and "
                    "band curves, the trace-basis law of every narrow basis triangle (by class) and "
                    "FarSize; prisms and pyramids are the tubes"}


def analyze(directory, manifest_path, basis_growths, census=None, mesh=None):
    manifest = json.loads(Path(manifest_path).read_text())
    recipe = manifest["ProductionRecipe"]
    options = dict(recipe["BuildCommandOptions"])
    model = recipe["BuildCostEstimate"]["TetrahedraPerCubicSize"]
    paths = source_paths(directory)
    mesh_recipe = json.loads(paths["MeshRecipe"].read_text())
    import tomllib
    process = tomllib.loads(paths["Process"].read_text())
    far = mesh_recipe["FarSizeOverRadius"] * process["Radius"]
    ratio, growth = float(options["--trace-basis-size-ratio"]), float(options["--far-growth"])
    vertices, triangles = read_basis(directory)
    rows = classify(vertices, triangles, ratio, far)
    current = estimate_decomposition(paths, options, model, rows, far, growth)
    # (a) Delaunay faces: the estimate with the re-triangulated basis (same vertices).
    flipped, faces = delaunay_faces(vertices, triangles)
    flipped_rows = classify(vertices, flipped, ratio, far)
    by_class = trace_integral_by_class(flipped_rows, far, growth)
    trace_current = current["Components"]["TraceBasis"]
    trace_flipped = (by_class["Needle"] + by_class["RightSliver"]) / model
    delaunay = {"RetriangulatedFaces": faces, "Sources": len(vertices),
                "Classes": {kind: sum(1 for row in flipped_rows if row["Class"] == kind)
                            for kind in ("Needle", "RightSliver", "Wide")},
                "MinimumAltitude": min(row["Altitude"] for row in flipped_rows),
                "TraceBasisTetrahedra": trace_flipped,
                "TraceBasisNeedles": by_class["Needle"] / model,
                "EstimatedElements": current["EstimatedElements"] - trace_current + trace_flipped}
    delaunay["ElementsOverCurrent"] = delaunay["EstimatedElements"] / current["EstimatedElements"]
    delaunay["PhysicsCostOverCurrent"] = delaunay["ElementsOverCurrent"]
    # (b) graded model.
    graded = []
    for basis_growth in basis_growths:
        record = graded_basis_model(vertices, triangles, ratio, far, growth, basis_growth, model)
        record["EstimatedElements"] = current["EstimatedElements"] - trace_current + record["TraceBasisTetrahedra"]
        record["ElementsOverCurrent"] = record["EstimatedElements"] / current["EstimatedElements"]
        record["SourcesOverCurrent"] = record["ModelledSources"] / len(vertices)
        record["PhysicsCostOverCurrent"] = record["ElementsOverCurrent"] * record["SourcesOverCurrent"]
        graded.append(record)
    result = {"Version": 1, "SourceDirectory": str(Path(directory).resolve()),
              "Manifest": str(Path(manifest_path).resolve()), "MaximumElements": manifest["Gates"]["MaximumElements"],
              "Sources": len(vertices), "Triangles": len(triangles),
              "NeedleRule": f"minimum altitude < {NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE} x shortest edge",
              "Current": current, "Delaunay": delaunay, "GradedModel": graded}
    if census is not None and mesh is not None:
        result["Actual"] = attribute_actual(census, mesh, ratio)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source_dir")
    parser.add_argument("--manifest", type=Path, default=HERE / "geometry-independence-suite.json")
    parser.add_argument("--census", type=Path)
    parser.add_argument("--mesh", type=Path)
    parser.add_argument("--basis-growth", type=float, nargs="+", default=[0.5, 1.0])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if (args.census is None) != (args.mesh is None):
        parser.error("--census and --mesh go together")
    result = analyze(args.source_dir, args.manifest, args.basis_growth, args.census, args.mesh)
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2, default=str) + "\n")
    current = result["Current"]
    print(f"estimate {current['EstimatedElements']:.0f} elements; needle share {current['NeedleShare']:.3f}; "
          f"classes {current['Classes']}")
    for name, share in sorted(current["Shares"].items(), key=lambda item: -item[1]):
        print(f"  {name:24s} {current['Components'][name]:12.0f}  {share:6.3f}")
    d = result["Delaunay"]
    print(f"delaunay faces {d['RetriangulatedFaces']}: classes {d['Classes']}, min altitude {d['MinimumAltitude']:.4g}, "
          f"estimate {d['EstimatedElements']:.0f} ({d['ElementsOverCurrent']:.3f} x current), sources unchanged")
    for g in result["GradedModel"]:
        print(f"graded basis growth {g['BasisGrowth']}: sources {g['ModelledSources']:.0f} ({g['SourcesOverCurrent']:.2f} x), "
              f"estimate {g['EstimatedElements']:.0f} ({g['ElementsOverCurrent']:.3f} x), physics cost {g['PhysicsCostOverCurrent']:.2f} x")
    if "Actual" in result:
        a = result["Actual"]
        print(f"actual {a['Elements']} elements; needle share {a['NeedleShare']:.3f}; shares "
              + ", ".join(f"{k} {v:.3f}" for k, v in a["Shares"].items()))


if __name__ == "__main__":
    main()
