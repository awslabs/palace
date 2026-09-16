#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Bound trace basis of a coupon and the cut-surface size rule it prescribes.

The trace basis (basis-contract.json with trace-vertices.csv and
trace-triangles.csv of the frozen campaign inputs) is the closed box
triangulation whose vertex hat functions are the response sources.  Its
vertices are stored in the process-library frame; the mesh frame follows from
the bound process library exactly as the campaign producer derived it
(generate_spatial_response.frame_from_geometry): local z is the fabrication
normal, local x the gap direction of the first edge.

Size rule (TraceBasisSizeRatio, the only parameter, dimensionless): on the cut
surface the local element size must not exceed TraceBasisSizeRatio times the
shortest edge of the basis triangle containing the point.  The hat of a basis
vertex varies linearly over every incident triangle, so its support is resolved
where the hat actually varies only when the whole triangle is discretized at
that scale; the per-edge alternative would resolve only the edges.  Away from
the surface the size grows with the recipe's own grading law up to the far
size, so far from narrow hats nothing changes.
"""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def process_frame(library, model_name):
    """Rows of the mesh (local) frame in process-library coordinates for one model.

    Local z is the process normal and local x the gap direction of the model's
    first edge (SpatialEdgeCluster) or arm; local y completes a right-handed
    frame.  A local point is `frame @ canonical`.
    """
    models = [model for model in library.get("Models", []) if model.get("Name") == model_name]
    if len(models) != 1:
        raise ValueError("Process library must contain the trace basis model exactly once")
    model = models[0]
    entries = model.get("Edges") if model.get("Topology") == "SpatialEdgeCluster" else model.get("Arms")
    if not isinstance(entries, list) or not entries:
        raise ValueError("Process library model has no complete edge geometry")
    normal = np.asarray(entries[0].get("ProcessNormal"), dtype=float).reshape(-1)
    gap = np.asarray(entries[0].get("GapDirection"), dtype=float).reshape(-1)
    if (normal.shape != (3,) or gap.shape != (3,) or not np.all(np.isfinite(normal)) or
            not np.all(np.isfinite(gap)) or np.linalg.norm(normal) <= 0 or
            np.linalg.norm(gap) <= 0):
        raise ValueError("Process library edge frame vectors are invalid")
    normal = normal / np.linalg.norm(normal)
    gap = gap / np.linalg.norm(gap)
    if abs(np.dot(normal, gap)) > 1e-9:
        raise ValueError("ProcessNormal and GapDirection must be orthogonal")
    axis_y = np.cross(normal, gap)
    axis_y /= np.linalg.norm(axis_y)
    return np.vstack((gap, axis_y, normal))


def _read_rows(path, columns):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or not set(columns) <= set(rows[0]):
        raise ValueError(f"{Path(path).name} lacks the columns {sorted(columns)}")
    return rows


def load_trace_basis(contract_path, vertices_path, triangles_path, process_library_path,
                     tolerance=1e-8):
    """Read, frame and validate the bound trace basis in the mesh (source-local) frame.

    Returns a dict with `Points` (n, 3), `Triangles` (m, 3) zero-based,
    `Lower`/`Upper` (the contract's mesh-frame box), `Frame`, `Model`, `Basis`
    (source index per vertex, 0 for none) and the input SHA-256 digests.  Fails
    closed when the counts differ from the contract geometry, when a vertex is
    not on the box surface, or when a triangle is degenerate or not in one box
    face plane (the contract asserts OffBoxTriangles 0).
    """
    contract = json.loads(Path(contract_path).read_text())
    geometry = contract.get("Geometry") if isinstance(contract, dict) else None
    if (not isinstance(contract, dict) or contract.get("Version") != 1 or
            not isinstance(contract.get("Model"), str) or not isinstance(geometry, dict) or
            geometry.get("OffBoxTriangles") != 0 or
            geometry.get("ClosedOrientedSurface") is not True):
        raise ValueError("Unsupported trace basis contract")
    lower = np.asarray(geometry.get("Lower"), dtype=float).reshape(-1)
    upper = np.asarray(geometry.get("Upper"), dtype=float).reshape(-1)
    if (lower.shape != (3,) or upper.shape != (3,) or not np.all(np.isfinite(lower)) or
            not np.all(np.isfinite(upper)) or np.any(upper <= lower)):
        raise ValueError("Trace basis contract box is invalid")
    library = json.loads(Path(process_library_path).read_text())
    frame = process_frame(library, contract["Model"])
    vertex_rows = _read_rows(vertices_path, ("vertex", "x", "y", "z", "basis", "conductor"))
    triangle_rows = _read_rows(triangles_path, ("triangle", "vertex_i", "vertex_j", "vertex_k"))
    if [int(row["vertex"]) for row in vertex_rows] != list(range(1, len(vertex_rows) + 1)):
        raise ValueError("Trace vertices must be numbered contiguously from one")
    if [int(row["triangle"]) for row in triangle_rows] != list(range(1, len(triangle_rows) + 1)):
        raise ValueError("Trace triangles must be numbered contiguously from one")
    canonical = np.array([[float(row[name]) for name in ("x", "y", "z")] for row in vertex_rows])
    basis = np.array([int(row["basis"]) for row in vertex_rows])
    triangles = np.array([[int(row[name]) - 1 for name in ("vertex_i", "vertex_j", "vertex_k")]
                          for row in triangle_rows], dtype=int)
    if (len(canonical) != geometry.get("Vertices") or len(triangles) != geometry.get("Triangles") or
            not np.all(np.isfinite(canonical))):
        raise ValueError("Trace basis files differ from the contract geometry counts")
    if (triangles.min() < 0 or triangles.max() >= len(canonical) or
            np.any(triangles[:, [0, 1, 2]] == triangles[:, [1, 2, 0]])):
        raise ValueError("Trace triangle connectivity is invalid")
    points = canonical @ frame.T
    scale = float(np.max(upper - lower))
    on_face = (np.abs(points - lower) <= tolerance * scale) | (np.abs(points - upper) <= tolerance * scale)
    inside = np.all(points >= lower - tolerance * scale, axis=1) & np.all(points <= upper + tolerance * scale, axis=1)
    if not np.all(inside & np.any(on_face, axis=1)):
        raise ValueError("Trace basis vertices are not on the contract box surface")
    xyz = points[triangles]
    if np.any(np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1) <= 0):
        raise ValueError("Degenerate trace basis triangle")
    if not np.all(np.any(np.all(on_face[triangles], axis=1), axis=1)):
        raise ValueError("Trace basis triangle is not in one box face plane")
    return {"Points": points, "Triangles": triangles, "Lower": lower, "Upper": upper,
            "Frame": frame, "Model": contract["Model"], "Basis": basis,
            "InputSHA256": {"BasisContract": sha256(contract_path),
                            "TraceVertices": sha256(vertices_path),
                            "TraceTriangles": sha256(triangles_path),
                            "ProcessLibrary": sha256(process_library_path)}}


def transform_trace_basis(basis, matrix):
    """The basis placed by a homogeneous rigid transform (row-major 4x4)."""
    matrix = np.asarray(matrix, dtype=float).reshape(4, 4)
    points = np.asarray(basis["Points"], dtype=float)
    return {**basis, "Points": points @ matrix[:3, :3].T + matrix[:3, 3]}


def triangle_edge_lengths(points, triangles):
    xyz = np.asarray(points, dtype=float)[np.asarray(triangles, dtype=int)]
    return np.stack([np.linalg.norm(xyz[:, i] - xyz[:, j], axis=1)
                     for i, j in ((0, 1), (1, 2), (2, 0))], axis=1)


def unique_edges(points, triangles):
    """The (k, 2, 3) endpoint coordinates of every distinct basis edge, in sorted
    vertex-index order."""
    triangles = np.asarray(triangles, dtype=int)
    edges = np.unique(np.sort(triangles[:, [(0, 1), (1, 2), (2, 0)]].reshape(-1, 2), axis=1), axis=0)
    return np.asarray(points, dtype=float)[edges]


def unique_edge_lengths(points, triangles):
    edges = unique_edges(points, triangles)
    return np.linalg.norm(edges[:, 0] - edges[:, 1], axis=1)


def requested_sizes(points, triangles, ratio):
    """TraceBasisSizeRatio times the shortest edge of every basis triangle."""
    if not np.isfinite(ratio) or ratio <= 0:
        raise ValueError("TraceBasisSizeRatio must be a positive finite dimensionless number")
    return ratio * triangle_edge_lengths(points, triangles).min(axis=1)


def point_triangle_distances(query, triangle):
    """Euclidean distance from every query point to one triangle (a, b, c)."""
    query = np.asarray(query, dtype=float).reshape(-1, 3)
    a, b, c = np.asarray(triangle, dtype=float).reshape(3, 3)
    ab, ac = b - a, c - a
    normal = np.cross(ab, ac)
    area2 = np.dot(normal, normal)
    if not area2 > 0:
        raise ValueError("Degenerate trace basis triangle")
    d = query - a
    # Barycentric coordinates of the plane projection.
    daa, dab, dbb = np.dot(ab, ab), np.dot(ab, ac), np.dot(ac, ac)
    dpa, dpb = d @ ab, d @ ac
    denominator = daa * dbb - dab * dab
    v = (dbb * dpa - dab * dpb) / denominator
    w = (daa * dpb - dab * dpa) / denominator
    inside = (v >= 0) & (w >= 0) & (v + w <= 1)
    result = np.abs(d @ normal) / np.sqrt(area2)
    if not np.all(inside):
        outside = ~inside
        best = np.full(int(outside.sum()), np.inf)
        for first, last in ((a, b), (b, c), (c, a)):
            vector = last - first
            delta = query[outside] - first
            parameter = np.clip((delta @ vector) / np.dot(vector, vector), 0.0, 1.0)
            best = np.minimum(best, np.linalg.norm(delta - parameter[:, None] * vector, axis=1))
        result[outside] = best
    return result


def trace_basis_sizes(query, basis, ratio, far_size, growth):
    """Isotropic size prescribed by the trace basis at every query point.

    size(x) = min(far_size, min over basis triangles T with ratio * minedge(T) <
    far_size of ratio * minedge(T) + growth * distance(x, T)): the surface rule
    inside each narrow triangle, the recipe's grading law away from it, and the
    far size everywhere else (`growth` is the metric stage's effective far
    growth, so the blend is the existing far/grading law).
    """
    query = np.asarray(query, dtype=float).reshape(-1, 3)
    if not np.all(np.isfinite([far_size, growth])) or far_size <= 0 or growth <= 0:
        raise ValueError("Invalid far size or growth")
    points, triangles = np.asarray(basis["Points"], dtype=float), np.asarray(basis["Triangles"], dtype=int)
    requested = requested_sizes(points, triangles, ratio)
    sizes = np.full(len(query), float(far_size))
    xyz = points[triangles]
    for index in np.flatnonzero(requested < far_size):
        # A triangle cannot lower the size where even its bounding box is too far.
        lower, upper = xyz[index].min(axis=0), xyz[index].max(axis=0)
        bound = np.linalg.norm(np.maximum(0.0, np.maximum(lower - query, query - upper)), axis=1)
        candidates = np.flatnonzero(requested[index] + growth * bound < sizes)
        if not len(candidates):
            continue
        distance = point_triangle_distances(query[candidates], xyz[index])
        sizes[candidates] = np.minimum(sizes[candidates], requested[index] + growth * distance)
    return sizes


def basis_statistics(basis, ratio, far_size):
    """Recorded, dimensionless-rule statistics of the bound basis against the far size."""
    points, triangles = basis["Points"], basis["Triangles"]
    edges = unique_edge_lengths(points, triangles)
    requested = requested_sizes(points, triangles, ratio)
    xyz = np.asarray(points)[np.asarray(triangles)]
    areas = 0.5 * np.linalg.norm(np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]), axis=1)
    narrow = requested < far_size
    return {"Vertices": int(len(points)), "Triangles": int(len(triangles)),
            "UniqueEdges": int(len(edges)),
            "MinimumBasisEdge": float(edges.min()), "MedianBasisEdge": float(np.median(edges)),
            "BasisEdgesBelowFarSize": int(np.sum(edges < far_size)),
            "TrianglesBelowFarSize": int(narrow.sum()),
            "AreaBelowFarSize": float(areas[narrow].sum()),
            "MinimumRequestedSize": float(requested.min()),
            "RequestedSizes": requested.tolist()}


def cut_surface_size_report(mesh_points, cut_triangles, basis, ratio, far_size, tolerance=1e-8):
    """Achieved cut-surface sizes of a mesh against the basis rule (reported, not gated).

    The size of a cut triangle is its longest edge; for every basis triangle
    with a requested size below the far size, the cut triangles whose centroid
    lies in it are measured along the direction of the basis triangle's
    shortest edge (the direction across which its hats vary fastest), and the
    largest extent over the requested size is the compliance ratio.  A basis
    triangle narrower than its requested size may centre no cut triangle; then
    only the count of cut triangles touching it is recorded.
    """
    mesh_points = np.asarray(mesh_points, dtype=float)
    cut_triangles = np.asarray(cut_triangles, dtype=int)
    xyz = mesh_points[cut_triangles]
    longest = triangle_edge_lengths(mesh_points, cut_triangles).max(axis=1)
    centroids = xyz.mean(axis=1)
    points, triangles = np.asarray(basis["Points"], dtype=float), np.asarray(basis["Triangles"], dtype=int)
    requested = requested_sizes(points, triangles, ratio)
    lengths = triangle_edge_lengths(points, triangles)
    scale = float(np.max(basis["Upper"] - basis["Lower"]))
    rows, worst = [], 0.0
    for index in np.flatnonzero(requested < far_size):
        corners = points[triangles[index]]
        inside = np.flatnonzero(point_triangle_distances(centroids, corners) <= tolerance * scale)
        touching = np.stack([point_triangle_distances(xyz[:, k], corners) <= tolerance * scale
                             for k in range(3)], axis=1).any(axis=1)
        shortest = int(np.argmin(lengths[index]))
        first, last = corners[shortest], corners[(shortest + 1) % 3]
        direction = (last - first) / np.linalg.norm(last - first)
        extent = np.ptp(xyz[inside] @ direction, axis=1) if len(inside) else np.zeros(0)
        achieved = float(extent.max()) if len(inside) else None
        compliance = None if achieved is None else achieved / float(requested[index])
        worst = max(worst, compliance or 0.0)
        rows.append({"BasisTriangle": int(index) + 1, "RequestedSize": float(requested[index]),
                     "CutTriangles": int(len(inside)), "CutTrianglesTouching": int(touching.sum()),
                     "MaximumExtentAlongShortestEdge": achieved,
                     "ExtentOverRequested": compliance})
    return {"CutTriangles": int(len(cut_triangles)),
            "SizeMeasure": "longest edge of each cut-surface triangle",
            "Minimum": float(longest.min()), "Median": float(np.median(longest)),
            "Maximum": float(longest.max()),
            "ComplianceMeasure": "largest extent of the cut triangles centred in a basis "
                                 "triangle along its shortest edge direction over the "
                                 "requested size (reported, not gated)",
            "MaximumExtentOverRequested": worst,
            "BasisTrianglesBelowFarSize": rows}
