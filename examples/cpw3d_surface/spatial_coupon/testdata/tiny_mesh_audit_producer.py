#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Write a real two-material Gmsh mesh for evidence-chain regression tests."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import meshio
import numpy as np


def produce(output, transform, scale=1.0):
    matrix = np.asarray(transform, dtype=float).reshape(4, 4)
    points = scale * np.array([[0., 0., 0.], [.08, 0., 0.], [0., .001, 0.],
                               [0., 0., .001], [0., 0., -.001], [10., 0., 0.],
                               [10.001, 0., 0.], [10., .001, 0.], [10., 0., .001]])
    tetrahedra = [[0, 1, 2, 3], [0, 2, 1, 4], [5, 6, 7, 8]]
    triangles = [[0, 1, 3], [0, 2, 3], [1, 2, 3],
                 [0, 2, 4], [0, 1, 4], [2, 1, 4], [0, 1, 2],
                 [5, 6, 7], [5, 8, 6], [5, 7, 8], [6, 8, 7]]
    extra_points = []
    for i in range(1, 6):
        for x in (i, i + 1):
            offset = len(points) + len(extra_points)
            extra_points.extend(([x, i, 0], [x + .08, i, 0],
                                 [x, i + .001, 0], [x, i, .001]))
            tetrahedra.append([offset, offset + 1, offset + 2, offset + 3])
            triangles.extend(([offset, offset + 1, offset + 2],
                              [offset, offset + 3, offset + 1],
                              [offset, offset + 2, offset + 3],
                              [offset + 1, offset + 3, offset + 2]))
    points = np.vstack((points, scale * np.asarray(extra_points)))
    homogeneous = np.column_stack((points, np.ones(len(points))))
    points = (homogeneous @ matrix.T)[:, :3]
    tetrahedra = np.asarray(tetrahedra)
    triangles = np.asarray(triangles)
    mesh = meshio.Mesh(points, [("triangle", triangles), ("tetra", tetrahedra)],
                       cell_data={"gmsh:physical": [
                                      np.array([1, 1, 1, 2, 2, 2, 3] +
                                               [1] * (len(triangles) - 7)),
                                      np.array([1, 7] + [1] * (len(tetrahedra) - 2))],
                                  "gmsh:geometrical": [np.ones(len(triangles), dtype=int),
                                                       np.ones(len(tetrahedra), dtype=int)]},
                       field_data={"substrate": np.array([1, 3]),
                                   "vacuum": np.array([7, 3]),
                                   "surface_1": np.array([1, 2]),
                                   "surface_2": np.array([2, 2]),
                                   "surface_3": np.array([3, 2])})
    meshio.write(output, mesh, file_format="gmsh22", binary=False)


def junction_curves(scale):
    """The fixture's cut-surface/material-interface junction lines: the three edges of
    the interface triangle [0, 1, 2] shared with the cut triangles of the first tet."""
    a, b, c = [0.0, 0.0, 0.0], [.08 * scale, 0.0, 0.0], [0.0, .001 * scale, 0.0]
    segments = [[*a, *b], [*b, *c], [*c, *a]]
    return {"Count": len(segments),
            "TotalLength": sum(math.dist(segment[:3], segment[3:]) for segment in segments),
            "Segments": segments, "Rule": "fixture"}


def corner_census(output, contract_path, radius, isotropic_size, etch_boundary=None,
                  scale=1.0):
    """Recorded corner-ball census of the fixture seed (schema of the production seeder)."""
    contract = json.loads(contract_path.read_text())
    corners = contract["SemanticCorners"]
    output.write_text(json.dumps({
        "Version": 1, "Frame": "SourceLocal", "SemanticCorners": corners,
        "CornerIsotropyRadius": radius, "IsotropicSize": isotropic_size,
        "EtchBoundary": "producer-default" if etch_boundary is None else str(etch_boundary),
        "EtchBoundarySHA256": None if etch_boundary is None else
                              hashlib.sha256(etch_boundary.read_bytes()).hexdigest(),
        # Simplified etch footprint polygons (schema of the production seeder): the
        # fixture footprint is one square whose duplicated vertex was merged.
        "FootprintCollinearTolerance": 1e-6,
        "FootprintSimplification": {"Rule": "fixture", "Polygons": 1, "RemovedVertices": 1,
                                    "MaximumRelativeDeviation": 0.0},
        "FootprintPolygons": [{"Conductor": 1, "Plane": 0.0, "Hole": False,
                               "Points": [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]],
                               "Simplification": {"OriginalVertices": 5, "Vertices": 4,
                                                  "RemovedVertexCount": 1,
                                                  "RemovedVertexIndices": [3],
                                                  "MaximumDeviation": 0.0,
                                                  "MaximumDeviationLocalScale": 10.0,
                                                  "MaximumRelativeDeviation": 0.0,
                                                  "Tolerance": 1e-6}}],
        "JunctionCurves": junction_curves(scale),
        "InterfaceAreaUnits": "um^2",
        # One row per contract boundary label, as written in the fixture seed.
        "InterfaceAreas": [{"Attribute": 1, "Name": "surface_1", "Triangles": 7,
                            "Area": 2.0},
                           {"Attribute": 2, "Name": "surface_2", "Triangles": 3, "Area": 2.0},
                           {"Attribute": 3, "Name": "surface_3", "Triangles": 1, "Area": 1.0}],
        "Corners": [{"Corner": index, "Point": corner, "BallEdges": 0}
                    for index, corner in enumerate(corners)],
        "LongitudinalFaces": [{"Surface": 1, "Triangles": 2, "InteriorNodes": 0,
                               "InteriorNodesAwayFromCorners": 0, "FullHeightTriangles": 0,
                               "FullHeightTrianglesAwayFromCorners": 0}]}, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path); parser.add_argument("transform", type=Path)
    parser.add_argument("--scale", type=float, default=1.0)
    # Bound source inputs: the fixture only requires them to exist, mirroring the
    # production seeder's positional signature and --mask/--boundary consumption.
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--boundary", type=Path, required=True)
    # Corner isotropy contract of the production seeder: contract corners, ball
    # radius, isotropic size and the recorded census artifact.
    parser.add_argument("--semantic-contract", type=Path, required=True)
    parser.add_argument("--corner-isotropy-radius", type=float, required=True)
    parser.add_argument("--lc-fine", type=float, required=True)
    parser.add_argument("--corner-census", type=Path, required=True)
    # The device etch footprint is consumed only when the case binds one.
    parser.add_argument("--etch-boundary", type=Path)
    args = parser.parse_args()
    if args.etch_boundary is not None and not args.etch_boundary.is_file():
        parser.error("the bound retained etch footprint must exist")
    if any(not path.is_file() for path in (args.signature, args.mask, args.boundary,
                                           args.semantic_contract)):
        parser.error("signature, mask, boundary, and semantic contract inputs must exist")
    if args.corner_census.exists():
        parser.error("corner census output must be fresh")
    produce(args.output, json.loads(args.transform.read_text()), args.scale)
    corner_census(args.corner_census, args.semantic_contract, args.corner_isotropy_radius,
                  args.lc_fine, args.etch_boundary, args.scale)


if __name__ == "__main__":
    main()
