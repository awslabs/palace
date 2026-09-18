#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests of check_hybrid_mesh.py on a synthetic prism + pyramid + tetrahedron mesh."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import meshio
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import check_hybrid_mesh as chm  # noqa: E402


def synthetic_mesh(flip_pyramid=False, duplicate_tet=False):
    # Prism (0,1,2 | 3,4,5) with the quad face (1,2,5,4) closed by a pyramid whose
    # apex 6 lies outside; a tetrahedron on the pyramid face (1,2,6)... kept simple:
    # prism + pyramid sharing one quad, plus a tetrahedron on pyramid face (2,5,6).
    points = np.array([
        [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0],
        [1.0, 1.0, 0.5],                       # pyramid apex outside the quad (1,2,5,4)
        [1.0, 1.5, 1.5],                       # tetra apex on pyramid face (2,5,6)
    ])
    prism = np.array([[0, 1, 2, 3, 4, 5]])
    pyramid = np.array([[1, 2, 5, 4, 6]]) if not flip_pyramid else np.array([[1, 4, 5, 2, 6]])
    tets = [[2, 5, 6, 7]] + ([[2, 5, 6, 7]] if duplicate_tet else [])
    tetra = np.array(tets)
    tri = np.array([[0, 1, 2], [3, 4, 5]])     # labeled prism triangle faces
    quad = np.array([[0, 1, 4, 3], [2, 0, 3, 5]])   # labeled prism quad faces
    cells = [meshio.CellBlock("wedge", prism), meshio.CellBlock("pyramid", pyramid),
             meshio.CellBlock("tetra", tetra), meshio.CellBlock("triangle", tri),
             meshio.CellBlock("quad", quad)]
    physical = [np.array([2]), np.array([2]), np.full(len(tetra), 2), np.array([7, 6001]),
                np.array([6001, 6001])]
    return meshio.Mesh(points, cells, cell_data={"gmsh:physical": physical})


class CheckHybridMeshTest(unittest.TestCase):
    def write(self, mesh):
        directory = Path(tempfile.mkdtemp())
        path = directory / "synthetic.msh"
        meshio.write(path, mesh, file_format="gmsh22", binary=True)
        return path, directory / "report.json"

    def test_conforming_positive_mesh(self):
        path, report = self.write(synthetic_mesh())
        result = chm.check(path, {6001: 3.0})
        self.assertTrue(result["PositivelyOriented"])
        self.assertEqual(result["VolumeElements"]["Prism"]["Count"], 1)
        self.assertEqual(result["VolumeElements"]["Pyramid"]["Count"], 1)
        self.assertEqual(result["QuadFaceOwners"], {"Prism": 2, "Prism+Pyramid": 1})
        self.assertEqual(result["FacesSharedByMoreThanTwo"], 0)
        # unlabeled exterior faces exist (the synthetic mesh labels only three faces)
        self.assertGreater(result["UnlabeledBoundaryFaces"], 0)
        self.assertFalse(result["Conforming"])
        self.assertAlmostEqual(result["ExpectedAreas"]["6001"]["Measured"], 2.0 + 0.5)
        self.assertIn("ScaledJacobian", result["VolumeElements"]["Prism"])
        self.assertGreater(result["VolumeElements"]["Prism"]["JacobianCondition"]["Max"], 1.0)

    def test_flipped_pyramid_is_reported(self):
        path, report = self.write(synthetic_mesh(flip_pyramid=True))
        result = chm.check(path, {})
        self.assertFalse(result["PositivelyOriented"])
        self.assertEqual(result["VolumeElements"]["Pyramid"]["NonPositive"], 1)

    def test_duplicate_tetrahedron_breaks_sharing(self):
        path, report = self.write(synthetic_mesh(duplicate_tet=True))
        result = chm.check(path, {})
        self.assertGreater(result["FacesSharedByMoreThanTwo"], 0)
        self.assertFalse(result["Conforming"])

    def test_command_line_writes_report(self):
        path, report = self.write(synthetic_mesh())
        subprocess.run([sys.executable, str(ROOT / "check_hybrid_mesh.py"), str(path), str(report),
                        "--expect-area", "6001=2.5"], check=True, capture_output=True)
        data = json.loads(report.read_text())
        self.assertEqual(data["TotalVolumeElements"], 3)
        self.assertEqual(len(data["SHA256"]), 64)


if __name__ == "__main__":
    unittest.main()
