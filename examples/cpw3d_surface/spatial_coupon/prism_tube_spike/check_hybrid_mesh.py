#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Independent check of a written hybrid Gmsh 2.2 mesh (prism-tube SPIKE, decision 37).

Reads the file with meshio and reports, without gating: element counts by type, node
count, positive orientation of every volume element (tets, prisms, pyramids: every
tetrahedron of the standard decomposition has positive volume), conformity (every
interior face is shared by exactly two volume elements with matching node sets, quad
faces by a prism and a prism or pyramid, triangular faces by two elements), boundary
coverage (every physical surface element is a face of exactly one volume element for
exterior labels or of two for interior interfaces), physical surface areas, volumes
per material and the SHA256 of the file. Writes a JSON report.

usage: check_hybrid_mesh.py mesh.msh report.json [--expect-area LABEL=AREA ...]
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import meshio
import numpy as np

TET_FACES = [(0, 1, 2), (0, 3, 1), (1, 3, 2), (0, 2, 3)]
PRISM_TRI_FACES = [(0, 1, 2), (3, 4, 5)]
PRISM_QUAD_FACES = [(0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)]
PYRAMID_QUAD_FACES = [(0, 1, 2, 3)]
PYRAMID_TRI_FACES = [(0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)]
TYPE_NAMES = {"tetra": "Tetrahedron", "wedge": "Prism", "pyramid": "Pyramid"}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def tet_volumes(points, conn):
    a, b, c, d = (points[conn[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)) / 6.0


def decomposition(kind, conn):
    """Tetrahedra of the standard decomposition of each element (list of index tuples)."""
    if kind == "tetra":
        return [(0, 1, 2, 3)]
    if kind == "wedge":
        return [(0, 1, 2, 3), (1, 2, 3, 4), (2, 3, 4, 5)]
    if kind == "pyramid":
        return [(0, 1, 2, 4), (0, 2, 3, 4), (0, 1, 3, 4), (1, 2, 3, 4)]
    raise ValueError(kind)


def element_faces(kind, conn):
    tri, quad = [], []
    if kind == "tetra":
        tri = TET_FACES
    elif kind == "wedge":
        tri, quad = PRISM_TRI_FACES, PRISM_QUAD_FACES
    elif kind == "pyramid":
        tri, quad = PYRAMID_TRI_FACES, PYRAMID_QUAD_FACES
    return tri, quad


# Corner (vertex) frames of the linear elements: at each listed vertex the three
# edges to its neighbours span the local Jacobian. Scaled Jacobian = normalized
# triple product (production tetrahedron_scaled_jacobian at every vertex);
# Jacobian condition = largest / smallest singular value of that edge matrix
# (production tetrahedron_aspect).
CORNER_FRAMES = {
    "tetra": [(0, 1, 2, 3), (1, 2, 0, 3), (2, 0, 1, 3), (3, 0, 2, 1)],
    "wedge": [(0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 1, 5), (3, 5, 4, 0), (4, 3, 5, 1), (5, 4, 3, 2)],
    "pyramid": [(0, 1, 3, 4), (1, 2, 0, 4), (2, 3, 1, 4), (3, 0, 2, 4)],
}


def corner_quality(points, kind, conn, batch=200000):
    """Per element: minimum scaled Jacobian and maximum Jacobian condition over corners."""
    scaled = np.full(len(conn), np.inf)
    condition = np.zeros(len(conn))
    for start in range(0, len(conn), batch):
        block = conn[start:start + batch]
        for (v, a, b, c) in CORNER_FRAMES[kind]:
            p0 = points[block[:, v]]
            e = np.stack([points[block[:, a]] - p0, points[block[:, b]] - p0,
                          points[block[:, c]] - p0], axis=2)          # (n, 3, 3) columns = edges
            triple = np.einsum("ij,ij->i", e[:, :, 0], np.cross(e[:, :, 1], e[:, :, 2]))
            norms = np.linalg.norm(e, axis=1).prod(axis=1)
            sj = np.abs(triple) / norms
            sigma = np.linalg.svd(e, compute_uv=False)
            cond = sigma[:, 0] / np.maximum(sigma[:, -1], 1e-300)
            scaled[start:start + batch] = np.minimum(scaled[start:start + batch], sj)
            condition[start:start + batch] = np.maximum(condition[start:start + batch], cond)
    return scaled, condition


def percentiles(values):
    return {"Min": float(values.min()), "P1": float(np.percentile(values, 1)),
            "P50": float(np.percentile(values, 50)), "Max": float(values.max())}


def face_key(nodes):
    return tuple(sorted(int(n) for n in nodes))


def surface_areas(points, cells, cell_data):
    areas = defaultdict(float)
    counts = defaultdict(Counter)
    for block, tags in zip(cells, cell_data):
        if block.type == "triangle":
            a, b, c = (points[block.data[:, i]] for i in range(3))
            area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)
        elif block.type == "quad":
            a, b, c, d = (points[block.data[:, i]] for i in range(4))
            area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1) + \
                0.5 * np.linalg.norm(np.cross(c - a, d - a), axis=1)
        else:
            continue
        for tag in np.unique(tags):
            mask = tags == tag
            areas[int(tag)] += float(area[mask].sum())
            counts[int(tag)][block.type] += int(mask.sum())
    return areas, counts


def check(path, expected_areas):
    mesh = meshio.read(path)
    points = mesh.points
    physical = mesh.cell_data["gmsh:physical"]
    report = {"Mesh": str(path), "SHA256": sha256(path), "Nodes": int(len(points)),
              "VolumeElements": {}, "SurfaceElements": {}}
    faces = defaultdict(list)          # face key -> list of (kind, index)
    quad_owner_kinds = defaultdict(list)
    volumes = defaultdict(float)
    total = 0
    nonpositive = 0
    for block, tags in zip(mesh.cells, physical):
        if block.type not in TYPE_NAMES:
            continue
        conn = block.data
        total += len(conn)
        minimum = np.full(len(conn), np.inf)
        positive = np.ones(len(conn), dtype=bool)
        signed_total = np.zeros(len(conn))
        for tet in decomposition(block.type, conn):
            v = tet_volumes(points, conn[:, list(tet)])
            positive &= v > 0.0
            minimum = np.minimum(minimum, v)
            signed_total += v
        if block.type == "pyramid":
            # the four tetrahedra double-cover the pyramid
            signed_total *= 0.5
        nonpositive += int((~positive).sum())
        for tag in np.unique(tags):
            volumes[int(tag)] += float(signed_total[tags == tag].sum())
        scaled, condition = corner_quality(points, block.type, conn)
        report["VolumeElements"][TYPE_NAMES[block.type]] = {
            "Count": int(len(conn)), "NonPositive": int((~positive).sum()),
            "MinimumVolume": float(signed_total.min()), "MinimumSubTetVolume": float(minimum.min()),
            "ScaledJacobian": percentiles(scaled), "ScaledJacobianBelow0.01": int((scaled < 0.01).sum()),
            "JacobianCondition": percentiles(condition),
            "JacobianConditionAbove1000": int((condition > 1000.0).sum())}
        tri, quad = element_faces(block.type, conn)
        for index, element in enumerate(conn):
            for f in tri:
                faces[face_key(element[list(f)])].append((block.type, index))
            for f in quad:
                key = face_key(element[list(f)])
                faces[key].append((block.type, index))
                quad_owner_kinds[key].append(block.type)
    report["TotalVolumeElements"] = total
    report["NonPositiveVolumeElements"] = nonpositive
    report["MaterialVolumes"] = {str(k): v for k, v in sorted(volumes.items())}
    # conformity: every face is shared by 1 (boundary) or 2 elements
    sharing = Counter(len(v) for v in faces.values())
    over_shared = sum(count for n, count in sharing.items() if n > 2)
    boundary_faces = {k for k, v in faces.items() if len(v) == 1}
    report["FaceSharing"] = {str(k): v for k, v in sorted(sharing.items())}
    report["FacesSharedByMoreThanTwo"] = over_shared
    quad_pairs = Counter(tuple(sorted(TYPE_NAMES[k] for k in kinds)) for kinds in quad_owner_kinds.values())
    report["QuadFaceOwners"] = {"+".join(k): v for k, v in sorted(quad_pairs.items())}
    # hanging nodes: a node of a boundary face of one element that lies strictly
    # inside a face of another would break sharing; sharing counts already catch
    # mismatched faces (they would appear as unshared interior faces). Interior
    # unshared faces are those boundary faces not covered by a physical surface.
    areas, counts = surface_areas(points, mesh.cells, physical)
    report["SurfaceAreas"] = {str(k): v for k, v in sorted(areas.items())}
    report["SurfaceElements"] = {str(k): dict(v) for k, v in sorted(counts.items())}
    labeled = {}
    uncovered = 0
    coverage = defaultdict(Counter)
    for block, tags in zip(mesh.cells, physical):
        if block.type not in ("triangle", "quad"):
            continue
        for element, tag in zip(block.data, tags):
            key = face_key(element)
            owners = len(faces.get(key, []))
            coverage[int(tag)][owners] += 1
            if owners == 0:
                uncovered += 1
    report["SurfaceCoverage"] = {str(k): {str(o): n for o, n in sorted(v.items())}
                                 for k, v in sorted(coverage.items())}
    report["SurfaceElementsNotOnVolumeFaces"] = uncovered
    labeled_keys = set()
    for block, tags in zip(mesh.cells, physical):
        if block.type in ("triangle", "quad"):
            for element in block.data:
                labeled_keys.add(face_key(element))
    report["UnlabeledBoundaryFaces"] = int(len(boundary_faces - labeled_keys))
    report["ExpectedAreas"] = {}
    for label, expected in expected_areas.items():
        measured = areas.get(label, 0.0)
        report["ExpectedAreas"][str(label)] = {
            "Expected": expected, "Measured": measured,
            "RelativeError": abs(measured - expected) / expected if expected else None}
    report["Conforming"] = over_shared == 0 and uncovered == 0 and \
        report["UnlabeledBoundaryFaces"] == 0
    report["PositivelyOriented"] = nonpositive == 0
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mesh")
    parser.add_argument("report")
    parser.add_argument("--expect-area", action="append", default=[],
                        help="LABEL=AREA reference (um^2), reported, not gated")
    args = parser.parse_args()
    expected = {}
    for item in args.expect_area:
        label, value = item.split("=")
        expected[int(label)] = float(value)
    report = check(args.mesh, expected)
    Path(args.report).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("TotalVolumeElements", "VolumeElements",
                                             "NonPositiveVolumeElements", "FaceSharing",
                                             "QuadFaceOwners", "SurfaceAreas",
                                             "UnlabeledBoundaryFaces", "Conforming",
                                             "PositivelyOriented", "SHA256")}, indent=1))


if __name__ == "__main__":
    main()
