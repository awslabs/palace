#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Read-only census of the seeded transverse edge layer in a seed, adapted or restored mesh.

Measures, against the recipe's recorded edge-layer spans, the tetrahedra by
distance from the spans (count, transverse edge sizes per shell, tangential
edge size, scaled Jacobian minimum and the cells below 0.01 / 0.02) and the
surface rows (transverse edge sizes per shell, tangential spacing), and, per
semantic corner ball, the edges and cells per shell of the recipe's corner law
(corner_ball_census).  Reported, not a gate: the physical gates are applied by
the audits.
"""
import argparse
import json
from pathlib import Path

import meshio
import numpy as np

from edge_volume_metric import segment_distances

EDGE_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
SHELLS_UM = ((0.0, 0.002), (0.002, 0.005), (0.005, 0.01), (0.01, 0.025), (0.025, 0.05))


def _span_distance(points, spans):
    distance = np.full(len(points), np.inf)
    tangent = np.zeros((len(points), 3))
    for span in spans:
        r, t = segment_distances(points, span)
        closer = r < distance
        distance[closer] = r[closer]
        tangent[closer] = t
    return distance, tangent


def _scaled_jacobian(xyz):
    jacobian = np.stack((xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0], xyz[:, 3] - xyz[:, 0]), axis=2)
    determinant = np.linalg.det(jacobian)
    return determinant / np.prod(np.linalg.norm(jacobian, axis=1), axis=1), determinant


def _edge_statistics(points, edges, distance, tangent):
    """Transverse and tangential edge sizes of a distinct edge set by the shell of
    the edge midpoint's distance from the spans."""
    vector = points[edges[:, 1]] - points[edges[:, 0]]
    length = np.linalg.norm(vector, axis=1)
    along = np.abs(np.einsum("ij,ij->i", vector, tangent[edges[:, 0]]))
    transverse = np.sqrt(np.maximum(length**2 - along**2, 0.0))
    middle = 0.5 * (distance[edges[:, 0]] + distance[edges[:, 1]])
    rows = {}
    for low, high in SHELLS_UM:
        shell = (middle > low) & (middle <= high)
        transverse_shell = shell & (transverse > along)
        tangential_shell = shell & (along >= transverse)
        if not shell.any():
            continue
        rows[f"{low * 1000:.0f}-{high * 1000:.0f}nm"] = {
            "Edges": int(shell.sum()),
            "TransverseP50Um": float(np.median(transverse[transverse_shell])) if transverse_shell.any() else None,
            "TransverseP90Um": float(np.percentile(transverse[transverse_shell], 90)) if transverse_shell.any() else None,
            "TangentialP50Um": float(np.median(along[tangential_shell])) if tangential_shell.any() else None}
    return rows


def edge_layer_census(mesh, recipe):
    layer = recipe.get("EdgeLayer")
    if not isinstance(layer, dict):
        raise ValueError("Recipe records no edge layer")
    spans = np.asarray(layer["Spans"], dtype=float).reshape(-1, 6)
    points = np.asarray(mesh.points, dtype=float)
    tetrahedra = np.concatenate([c.data for c in mesh.cells if c.type == "tetra"])
    triangles = np.concatenate([c.data for c in mesh.cells if c.type == "triangle"])
    distance, tangent = _span_distance(points, spans)
    xyz = points[tetrahedra]
    centroid_distance, _ = _span_distance(xyz.mean(axis=1), spans)
    scaled, determinant = _scaled_jacobian(xyz)
    if np.any(determinant <= 0):
        raise ValueError("Census received a nonpositive tetrahedron")
    lengths = np.stack([np.linalg.norm(xyz[:, i] - xyz[:, j], axis=1) for i, j in EDGE_PAIRS], axis=1)
    report = {"Tetrahedra": int(len(tetrahedra)), "Vertices": int(len(points)),
              "SurfaceTriangles": int(len(triangles)), "Spans": int(len(spans)),
              "TotalSpanLength": float(layer["TotalSpanLength"]), "EdgeSize": float(layer["EdgeSize"]),
              "GrowthRatio": float(layer["GrowthRatio"]), "Aspect": float(layer.get("Aspect", np.nan)),
              "RowOffsets": list(layer["RowOffsets"]), "LayerThickness": float(layer["LayerThickness"])}
    cells = {}
    for radius in (0.005, 0.01, 0.05, 0.1):
        inside = centroid_distance <= radius
        cells[f"Within{radius}"] = {
            "Count": int(inside.sum()),
            "MinimumScaledJacobian": float(scaled[inside].min()) if inside.any() else None,
            "CellsBelow0.01": int(np.sum(scaled[inside] < 0.01)),
            "CellsBelow0.02": int(np.sum(scaled[inside] < 0.02)),
            "MaximumEdgeAspect": float((lengths[inside].max(axis=1) / lengths[inside].min(axis=1)).max()) if inside.any() else None,
            "ShortestEdgeP50Um": float(np.median(lengths[inside].min(axis=1))) if inside.any() else None}
    report["LayerCells"] = cells
    band = tetrahedra[centroid_distance <= 0.05]
    edges = np.unique(np.sort(band[:, EDGE_PAIRS].reshape(-1, 2), axis=1), axis=0)
    report["TetrahedronEdgesByShell"] = _edge_statistics(points, edges, distance, tangent)
    surface = triangles[_span_distance(points[triangles].mean(axis=1), spans)[0] <= layer["LayerThickness"] + layer["EdgeSize"]]
    surface_edges = np.unique(np.sort(surface[:, [(0, 1), (1, 2), (0, 2)]].reshape(-1, 2), axis=1), axis=0)
    report["SurfaceTrianglesInLayer"] = int(len(surface))
    report["SurfaceEdgesByShell"] = _edge_statistics(points, surface_edges, distance, tangent)
    report["VerticesOnSpans"] = int(np.sum(distance <= 1e-9))
    return report


def corner_ball_census(mesh, recipe):
    """Per semantic corner and per shell of the recipe's corner law (CornerGrading
    ShellRadii, or the single NormalSize ball): the tetrahedron edges by midpoint
    distance (count, P50/P90 length against the law's size) and the cells by
    centroid, with the scaled-Jacobian minimum and the cells by decade.  Reported,
    not a gate."""
    corners = np.asarray(recipe["TruePhysicalCorners"], dtype=float).reshape(-1, 3)
    radius = float(recipe["CornerIsotropyRadius"])
    normal = float(recipe["NormalSize"])
    grading = recipe.get("CornerGrading")
    if isinstance(grading, dict):
        boundaries = [0.0] + [float(value) for value in grading["ShellRadii"]]
        sizes = [float(value) for value in grading["ShellSizes"]]
    else:
        boundaries, sizes = [0.0, radius], [normal]
    points = np.asarray(mesh.points, dtype=float)
    tetrahedra = np.concatenate([c.data for c in mesh.cells if c.type == "tetra"])
    xyz = points[tetrahedra]
    centroids = xyz.mean(axis=1)
    scaled, determinant = _scaled_jacobian(xyz)
    rows = []
    for corner in corners:
        ball = np.flatnonzero(np.linalg.norm(centroids - corner, axis=1) <= radius)
        edges = np.unique(np.sort(tetrahedra[ball][:, EDGE_PAIRS].reshape(-1, 2), axis=1), axis=0)
        middle = np.linalg.norm(0.5 * (points[edges[:, 0]] + points[edges[:, 1]]) - corner, axis=1)
        length = np.linalg.norm(points[edges[:, 1]] - points[edges[:, 0]], axis=1)
        cell_distance = np.linalg.norm(centroids[ball] - corner, axis=1)
        shells = []
        for inner, outer, size in zip(boundaries[:-1], boundaries[1:], sizes):
            shell = (middle > inner) & (middle <= outer)
            cells = (cell_distance > inner) & (cell_distance <= outer)
            shells.append({"InnerRadiusUm": inner, "OuterRadiusUm": outer, "TargetSizeUm": size,
                           "Edges": int(shell.sum()), "Cells": int(cells.sum()),
                           "EdgeP50Um": float(np.median(length[shell])) if shell.any() else None,
                           "EdgeP90Um": float(np.percentile(length[shell], 90)) if shell.any() else None,
                           "EdgesOverSqrt2TargetSize": int(np.sum(length[shell] > np.sqrt(2.0) * size)),
                           "MinimumScaledJacobian": float(scaled[ball][cells].min()) if cells.any() else None})
        rows.append({"Point": corner.tolist(), "Cells": int(len(ball)),
                     "MinimumScaledJacobian": float(scaled[ball].min()) if len(ball) else None,
                     "CellsBelow0.01": int(np.sum(scaled[ball] < 0.01)),
                     "CellsBelow0.02": int(np.sum(scaled[ball] < 0.02)),
                     "PositiveOrientation": bool(np.all(determinant[ball] > 0)), "Shells": shells})
    return {"Radius": radius, "CornerSize": (float(grading["CornerSize"]) if isinstance(grading, dict) else None),
            "ShellRadii": boundaries[1:], "ShellSizes": sizes, "Corners": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("recipe", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mesh, recipe = meshio.read(args.mesh), json.loads(args.recipe.read_text())
    report = edge_layer_census(mesh, recipe)
    report["CornerBalls"] = corner_ball_census(mesh, recipe)
    report["Mesh"] = str(args.mesh.resolve())
    text = json.dumps(report, indent=2)
    print(text)
    if args.output is not None:
        args.output.write_text(text + "\n")


if __name__ == "__main__":
    main()
