#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Write a real two-material Gmsh mesh for evidence-chain regression tests."""
import argparse
import json
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path); parser.add_argument("transform", type=Path)
    parser.add_argument("--scale", type=float, default=1.0)
    args = parser.parse_args()
    produce(args.output, json.loads(args.transform.read_text()), args.scale)


if __name__ == "__main__":
    main()
