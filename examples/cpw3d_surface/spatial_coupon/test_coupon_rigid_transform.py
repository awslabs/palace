#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Rigid source-production contract and mesh covariance tests."""
import csv
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

from general_mesh_audit_producer import _ownership_report
from prepare_edge_metric_scout import cluster_planar_supports, validate_transformed_supports
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

    def test_metric_support_validation_rejects_tampered_transformed_edge(self):
        matrix = validate_rigid_transform(ROTATE_Z)
        contract = transform_semantic_contract(
            json.loads((SOURCE / "semantic-contract.json").read_text()), matrix)
        contract["CanonicalTransformSHA256"] = "1" * 64
        contract["SourceSemanticContractSHA256"] = "2" * 64
        supports = transformed_supports(
            SOURCE, matrix, transform_sha256="1" * 64, semantic_sha256="2" * 64)
        segments = []
        for edge in supports["Edges"]:
            point = np.asarray(edge["Point"])
            tangent = np.asarray(edge["TangentDirection"])
            segments.append([point + edge["Interval"][0] * tangent,
                             point + edge["Interval"][1] * tangent])
        validate_transformed_supports(supports, contract, segments)
        tampered_hash = json.loads(json.dumps(supports))
        tampered_hash["CanonicalTransformSHA256"] = "3" * 64
        with self.assertRaisesRegex(ValueError, "provenance differs"):
            validate_transformed_supports(tampered_hash, contract, segments)
        tampered = json.loads(json.dumps(supports))
        tampered["Edges"][0]["Point"][0] += 0.01
        with self.assertRaisesRegex(ValueError, "seed-derived features"):
            validate_transformed_supports(tampered, contract, segments)

    def test_numerically_noisy_coplanar_supports_are_clustered(self):
        attributes = np.array([1, 1, 1, 2])
        normals = np.array([[1.0, 0.0, 0.0],
                            [1.0, 6e-9, 0.0],
                            [1.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0]])
        points = np.array([[9.8333333333, 0.0, 0.0],
                           [9.8333333340, 1.0, 0.0],
                           [9.8333340, 0.0, 0.0],
                           [9.8333333333, 0.0, 0.0]])
        planes, patch = cluster_planar_supports(attributes, normals, points)
        self.assertEqual(len(planes), 3)
        np.testing.assert_array_equal(patch, [0, 0, 1, 2])

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
            self.assertEqual(actual, expected)
            self.assertTrue({5001, 5101, 5002, 5102,
                             6001, 6101, 6002, 6102} <= actual)

            def ownership(name):
                report_path = root / f"{name}-ownership.csv"
                with report_path.open(newline="") as stream:
                    report = list(csv.DictReader(stream))
                with (Path(str(report_path) + ".elements.csv")).open(newline="") as stream:
                    elements = list(csv.DictReader(stream))
                with (Path(str(report_path) + ".quadrature.csv")).open(newline="") as stream:
                    quadrature = list(csv.DictReader(stream))
                return report, elements, quadrature

            local_report, local_elements, local_quadrature = ownership("multislot-local")
            global_report, global_elements, global_quadrature = ownership(
                "multislot-transformed")
            self.assertEqual(local_report, global_report)
            self.assertEqual(local_elements, global_elements)
            self.assertEqual(local_quadrature, global_quadrature)
            self.assertEqual(sum(int(row["elements"]) for row in local_report),
                             len(local_elements))
            element_ids = [row["element"] for row in local_elements]
            self.assertEqual(len(element_ids), len(set(element_ids)))
            summary = local_report[0]
            self.assertEqual((summary["quadrature_rule"],
                              int(summary["quadrature_order"])), ("Gauss4", 4))
            self.assertEqual(int(summary["quadrature_positive_weights"]), 1)
            self.assertEqual(int(summary["quadrature_unmatched"]), 0)
            self.assertEqual(int(summary["quadrature_overlaps"]), 0)
            self.assertLessEqual(float(summary["quadrature_relative_closure"]),
                                 float(summary["quadrature_closure_tolerance"]))
            semantic = json.loads(MULTISLOT["semantic"].read_text())
            expected_owners = sorted(item["Attribute"] for item in semantic["BoundaryLabels"]
                                     if item["Role"] not in semantic["CutSurfaceRoles"])
            self.assertEqual([int(row["attribute"]) for row in local_quadrature],
                             expected_owners)
            ownership_audit = _ownership_report(
                root / "multislot-local-ownership.csv",
                root / "multislot-local-ownership.csv.quadrature.csv", semantic)
            self.assertEqual(ownership_audit["ResponseOwnership"]["OwnerAttributes"],
                             expected_owners)
            self.assertLessEqual(ownership_audit["ResponseOwnership"]
                                 ["OwnerPartitionRelativeClosure"], 1e-12)
            owned_measure = sum(float(row["measure"]) for row in local_quadrature)
            self.assertAlmostEqual(owned_measure,
                                   float(summary["quadrature_whole_measure"]), places=11)
            quadrature_path = root / "multislot-local-ownership.csv.quadrature.csv"
            original_quadrature = quadrature_path.read_text()
            bad_rows = list(local_quadrature)
            for name, rows in (
                    ("changed", [*bad_rows[:-1], {**bad_rows[-1], "attribute": "9999"}]),
                    ("missing", bad_rows[:-1]),
                    ("duplicate", [*bad_rows[:-1], {**bad_rows[-1],
                                                     "attribute": bad_rows[0]["attribute"]}])):
                with self.subTest(owner_mutation=name):
                    with quadrature_path.open("w", newline="") as stream:
                        writer = csv.DictWriter(stream, fieldnames=["attribute", "measure"])
                        writer.writeheader(); writer.writerows(rows)
                    with self.assertRaises(ValueError):
                        _ownership_report(root / "multislot-local-ownership.csv",
                                          quadrature_path, semantic)
            quadrature_path.write_text(original_quadrature)
            # Whole-element ambiguity is preserved as a non-authoritative diagnostic;
            # it is not response ownership and must not drop or duplicate a triangle.
            self.assertGreater(sum(int(row["unresolved_elements"])
                                   for row in local_report), 0)

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
