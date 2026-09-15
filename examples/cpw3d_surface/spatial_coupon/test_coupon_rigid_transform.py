#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Rigid source-production contract and mesh covariance tests."""
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import meshio
import numpy as np

from transform_coupon_source_contract import (
    transform_semantic_contract, transform_vector, transformed_supports,
    validate_rigid_transform,
)

REPO = HERE.parents[2]
MESHER = HERE / "mesh_spatial_coupon.jl"
SOURCE = HERE / "testdata" / "four-edge-9d2cb9bbb3fe"
ANGLE = 0.63
ROTATE_Z = [math.cos(ANGLE), -math.sin(ANGLE), 0.0, 0.0,
            math.sin(ANGLE), math.cos(ANGLE), 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0]
IDENTITY = [1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
TILT = 0.41
TILTED_TRANSLATED = [
    math.cos(TILT) * math.cos(ANGLE), -math.cos(TILT) * math.sin(ANGLE),
    math.sin(TILT), 1.2,
    math.sin(ANGLE), math.cos(ANGLE), 0.0, -0.7,
    -math.sin(TILT) * math.cos(ANGLE), math.sin(TILT) * math.sin(ANGLE),
    math.cos(TILT), 0.9,
    0.0, 0.0, 0.0, 1.0]
MULTISLOT = {
    "signature": HERE / "testdata/six-edge-cluster-signature.csv",
    "mask": HERE / "testdata/six-edge-cluster-mask.csv",
    "boundary": HERE / "testdata/six-edge-cluster-boundary.csv",
    "semantic": HERE / "testdata/six-edge-semantic.json",
}


class RigidContractTest(unittest.TestCase):
    def test_source_vectors_and_semantic_points_transform_without_label_changes(self):
        matrix = validate_rigid_transform(ROTATE_Z)
        contract = json.loads((SOURCE / "semantic-contract.json").read_text())
        transformed = transform_semantic_contract(contract, matrix)
        self.assertEqual(transformed["VolumeMaterials"], contract["VolumeMaterials"])
        self.assertEqual(transformed["BoundaryLabels"], contract["BoundaryLabels"])
        self.assertEqual(transformed["ProtectedSupports"], contract["ProtectedSupports"])
        expected = [math.cos(ANGLE) * 10.0, math.sin(ANGLE) * 10.0, 0.0]
        np.testing.assert_allclose(transformed["SemanticCorners"][0], expected,
                                   rtol=0.0, atol=1e-14)
        supports = transformed_supports(SOURCE, matrix)
        first = supports["Edges"][0]
        for name in ("GapDirection", "TangentDirection", "ProcessNormalDirection"):
            self.assertAlmostEqual(np.linalg.norm(first[name]), 1.0, places=14)
        original_normal = [0.0, 0.0, 1.0]
        np.testing.assert_allclose(first["ProcessNormalDirection"],
                                   transform_vector(matrix, original_normal),
                                   rtol=0.0, atol=1e-14)
        self.assertEqual(len(supports["Edges"]), 4)

    def test_nonrigid_singular_and_reflecting_transforms_are_rejected(self):
        for transform in (
            [2.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0],
            [1.0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0],
            [-1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0],
        ):
            with self.subTest(transform=transform):
                with self.assertRaises(ValueError):
                    validate_rigid_transform(transform)


class RigidProducerIntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.julia = shutil.which(os.environ.get("JULIA", "julia"))
        if cls.julia is None:
            raise unittest.SkipTest("Julia is not available")
        cls.project = REPO / "test" / "examples"
        probe = subprocess.run(
            [cls.julia, "--startup-file=no", f"--project={cls.project}", "-e",
             "import Gmsh: gmsh; gmsh.initialize(); gmsh.finalize()"],
            capture_output=True, text=True, check=False)
        if probe.returncode:
            raise unittest.SkipTest("The test Julia project has no Gmsh dependency")

    def produce(self, root, name, transform=None, source=None, ownership=False):
        output = root / f"{name}.msh"
        source = source or {
            "signature": SOURCE / "mesh-signature.csv",
            "mask": SOURCE / "plan-view-mask.csv",
            "boundary": SOURCE / "plan-view-boundary.csv"}
        command = [
            self.julia, "--startup-file=no", f"--project={self.project}", str(MESHER),
            str(source["signature"]), "fabricated", str(output),
            "--mask", str(source["mask"]),
            "--boundary", str(source["boundary"]),
            "--radius", "2", "--metal-thickness", "0.1", "--overetch", "0.05",
            "--sidewall-angle", "90", "--top-radius", "0", "--bottom-radius", "0",
            "--lc-fine", "0.2", "--lc-tangent", "0.4", "--lc-far", "0.6",
            "--mesh-order", "1", "--max-nodes", "1000000", "--max-elements", "1000000",
        ]
        if ownership:
            command += ["--interface-ownership-report", str(root / f"{name}-ownership.csv")]
        if transform is not None:
            command += ["--rigid-transform", ",".join(format(value, ".17g")
                                                        for value in transform)]
        result = subprocess.run(command, cwd=REPO, capture_output=True, text=True,
                                check=False, timeout=180)
        if result.returncode:
            self.fail(f"producer failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")
        return output

    def test_identity_equivalence_and_rotation_covariance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = self.produce(root, "baseline")
            identity = self.produce(root, "identity", IDENTITY)
            rotated = self.produce(root, "rotated", ROTATE_Z)
            self.assertEqual(hashlib.sha256(baseline.read_bytes()).hexdigest(),
                             hashlib.sha256(identity.read_bytes()).hexdigest())
            left, right = meshio.read(baseline), meshio.read(rotated)
            self.assertEqual([cell.type for cell in left.cells],
                             [cell.type for cell in right.cells])
            for first, second in zip(left.cells, right.cells):
                np.testing.assert_array_equal(first.data, second.data)
            for first, second in zip(left.cell_data["gmsh:physical"],
                                     right.cell_data["gmsh:physical"]):
                np.testing.assert_array_equal(first, second)
            matrix = np.asarray(ROTATE_Z).reshape(4, 4)
            expected = left.points @ matrix[:3, :3].T + matrix[:3, 3]
            np.testing.assert_allclose(right.points, expected, rtol=0.0, atol=2e-14)
            tetrahedra = np.concatenate([cell.data for cell in left.cells
                                         if cell.type == "tetra"])
            def determinants(points):
                xyz = points[tetrahedra]
                return np.linalg.det(np.stack((xyz[:, 1] - xyz[:, 0],
                                               xyz[:, 2] - xyz[:, 0],
                                               xyz[:, 3] - xyz[:, 0]), axis=2))
            before, after = determinants(left.points), determinants(right.points)
            self.assertTrue(np.all(before * after > 0.0))
            np.testing.assert_allclose(after, before, rtol=2e-12, atol=1e-14)

    def test_multislot_tilted_translation_preserves_exact_slot_conductor_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            baseline = self.produce(root, "multislot-local", source=MULTISLOT,
                                    ownership=True)
            transformed = self.produce(root, "multislot-transformed", TILTED_TRANSLATED,
                                       source=MULTISLOT, ownership=True)
            left, right = meshio.read(baseline), meshio.read(transformed)
            for first, second in zip(left.cells, right.cells):
                np.testing.assert_array_equal(first.data, second.data)
            for first, second in zip(left.cell_data["gmsh:physical"],
                                     right.cell_data["gmsh:physical"]):
                np.testing.assert_array_equal(first, second)
            matrix = np.asarray(TILTED_TRANSLATED).reshape(4, 4)
            np.testing.assert_allclose(
                right.points, left.points @ matrix[:3, :3].T + matrix[:3, 3],
                rtol=0.0, atol=3e-14)
            expected = {item["Attribute"] for item in
                        json.loads(MULTISLOT["semantic"].read_text())["BoundaryLabels"]}
            actual = set(np.concatenate([
                values for cell, values in zip(right.cells,
                    right.cell_data["gmsh:physical"]) if cell.type == "triangle"]))
            # The fixture has no exposed substrate-air patch for the optional
            # 3000/3001 roles; every geometrically present declared label is exact.
            self.assertEqual(actual, expected - {3000, 3001})
            self.assertTrue({5001, 5101, 5002, 5102,
                             6001, 6101, 6002, 6102} <= actual)

    def test_julia_option_rejects_nonrigid_and_singular_matrices(self):
        invalid = (
            [2.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0],
            [1.0, 0, 0, 0, 0, 0.0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 1.0],
        )
        for transform in invalid:
            value = ",".join(map(str, transform))
            code = f'include(raw"{MESHER}"); parse_rigid_transform("{value}")'
            result = subprocess.run(
                [self.julia, "--startup-file=no", f"--project={self.project}", "-e", code],
                cwd=REPO, capture_output=True, text=True, check=False, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Rigid transform", result.stderr)


if __name__ == "__main__":
    unittest.main()
