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
                               [0., 0., .001], [0., 0., -.001]])
    homogeneous = np.column_stack((points, np.ones(len(points))))
    points = (homogeneous @ matrix.T)[:, :3]
    tetrahedra = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
    triangles = np.array([[0, 1, 3], [0, 2, 3], [1, 2, 3],
                          [0, 2, 4], [0, 1, 4], [2, 1, 4], [0, 1, 2]])
    mesh = meshio.Mesh(points, [("triangle", triangles), ("tetra", tetrahedra)],
                       cell_data={"gmsh:physical": [np.array([1, 1, 1, 2, 2, 2, 3]),
                                                    np.array([1, 7])],
                                  "gmsh:geometrical": [np.ones(7, dtype=int),
                                                       np.ones(2, dtype=int)]})
    meshio.write(output, mesh, file_format="gmsh22", binary=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path); parser.add_argument("transform", type=Path)
    parser.add_argument("--scale", type=float, default=1.0)
    args = parser.parse_args()
    produce(args.output, json.loads(args.transform.read_text()), args.scale)


if __name__ == "__main__":
    main()
