#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Rigid source-production contract and mesh covariance tests."""
import csv
import hashlib
import json
import math
import os
import re
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
        """The metric's PlanarSupports use the shared plane-equivalence rule: roundoff
        noise on one CAD plane is one support, a parallel plane offset by a
        resolvable fraction of the local size is another, and labels never merge."""
        attributes = np.array([1, 1, 1, 2])
        def triangle(x, y, z, tilt=0.0):
            corners = np.array([[x, y, z], [x, y + .1, z], [x, y, z + .1]])
            corners[:, 0] += [0.0, 0.0, tilt * .1]
            return corners
        xyz = np.array([triangle(9.8333333333, 0.0, 0.0),
                        triangle(9.8333333340, 1.0, 0.0, tilt=6e-9),
                        triangle(9.8333340, 0.0, 0.0),
                        triangle(9.8333333333, 0.0, 0.0)])
        normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
        normals /= np.linalg.norm(normals, axis=1)[:, None]
        planes, patch = cluster_planar_supports(attributes, normals, xyz)
        self.assertEqual(len(planes), 3)
        np.testing.assert_array_equal(patch, [0, 0, 1, 2])
        np.testing.assert_array_equal(planes[:, 0], [1, 1, 2])
        np.testing.assert_allclose(planes[:, 4], [9.8333333333, 9.8333340, 9.8333333333])
        # Same planes in a rotated frame: same clustering, rotated representatives.
        angle = 0.63
        rotation = np.array([[math.cos(angle), -math.sin(angle), 0.0],
                             [math.sin(angle), math.cos(angle), 0.0], [0.0, 0.0, 1.0]])
        rotated_planes, rotated_patch = cluster_planar_supports(
            attributes, normals @ rotation.T, xyz @ rotation.T)
        np.testing.assert_array_equal(rotated_patch, patch)
        np.testing.assert_allclose(rotated_planes[:, 1:4], planes[:, 1:4] @ rotation.T, atol=1e-15)
        np.testing.assert_allclose(rotated_planes[:, 4], planes[:, 4], rtol=1e-14)

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

    def produce(self, root, name, transform=None, source=None, ownership=False,
                corner_isotropy=None, expect_failure=None):
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
        if corner_isotropy is not None:
            command += corner_isotropy
        result = subprocess.run(command, cwd=REPO, capture_output=True, text=True,
                                check=False, timeout=180)
        if expect_failure is not None:
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(expect_failure, result.stderr)
            self.assertFalse(output.exists())
            return None
        if result.returncode:
            self.fail(f"producer failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}")
        return output

    def corner_isotropy(self, root, name, transform=None, radius="0.4", census=True):
        """Seed corner-isotropy options for a contract placed by `transform`."""
        contract = json.loads((SOURCE / "semantic-contract.json").read_text())
        if transform is not None:
            contract = transform_semantic_contract(contract, validate_rigid_transform(transform))
        path = root / f"{name}-semantic.json"; path.write_text(json.dumps(contract))
        options = ["--semantic-contract", str(path), "--corner-isotropy-radius", radius]
        if census:
            options += ["--corner-census", str(root / f"{name}-census.json")]
        return options

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

    def test_seed_corner_isotropy_is_rotation_covariant_and_recorded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plain = self.produce(root, "plain")
            identity = self.produce(root, "identity", IDENTITY,
                                    corner_isotropy=self.corner_isotropy(root, "identity"))
            rotated = self.produce(root, "rotated", ROTATE_Z,
                                   corner_isotropy=self.corner_isotropy(root, "rotated", ROTATE_Z))
            left, right = meshio.read(identity), meshio.read(rotated)
            for first, second in zip(left.cells, right.cells):
                np.testing.assert_array_equal(first.data, second.data)
            for first, second in zip(left.cell_data["gmsh:physical"],
                                     right.cell_data["gmsh:physical"]):
                np.testing.assert_array_equal(first, second)
            matrix = np.asarray(ROTATE_Z).reshape(4, 4)
            np.testing.assert_allclose(right.points,
                                       left.points @ matrix[:3, :3].T + matrix[:3, 3],
                                       rtol=0.0, atol=3e-14)
            # The corner ball changed the seed: the plain seed is a different mesh.
            self.assertNotEqual(hashlib.sha256(plain.read_bytes()).hexdigest(),
                                hashlib.sha256(identity.read_bytes()).hexdigest())
            census = json.loads((root / "identity-census.json").read_text())
            rotated_census = json.loads((root / "rotated-census.json").read_text())
            contract = json.loads((SOURCE / "semantic-contract.json").read_text())
            self.assertEqual(census["Version"], 1)
            self.assertEqual(census["Frame"], "SourceLocal")
            self.assertEqual(census["SemanticCorners"], contract["SemanticCorners"])
            self.assertEqual(census["CornerIsotropyRadius"], 0.4)
            self.assertEqual(census["IsotropicSize"], 0.2)
            self.assertEqual(len(census["Corners"]), len(contract["SemanticCorners"]))
            self.assertGreater(census["CornerIsotropicLongitudinalCurves"], 0)
            self.assertLessEqual(census["CornerIsotropicLongitudinalCurves"],
                                 census["LongitudinalCurves"])
            for row in census["Corners"]:
                self.assertGreater(row["BallEdges"], 0)
                self.assertGreaterEqual(row["FractionOverSqrt2IsotropicSize"], 0.0)
                self.assertLessEqual(row["EdgeMaximum"], 0.4)
                self.assertGreater(row["IncidentMaximumAspect"], 1.0)
            # Ridge-to-ridge faces are censused (interior row, full-height triangles).
            # With this test's coarse tangential size the 0.1 sidewalls are legitimately
            # full-height; the recipe-size property is covered by the Julia unit tests.
            self.assertGreater(census["CornerLawReach"], census["CornerIsotropyRadius"])
            self.assertGreater(len(census["LongitudinalFaces"]), 0)
            for row in census["LongitudinalFaces"]:
                self.assertGreaterEqual(row["LongitudinalCurves"], 2)
                self.assertLessEqual(row["FullHeightTrianglesAwayFromCorners"],
                                     row["FullHeightTriangles"])
                self.assertLessEqual(row["FullHeightTriangles"], row["Triangles"])
                self.assertEqual(len(row["InteriorNodeHistogramAlongEdge"]),
                                 census["LongitudinalFaceHistogramBins"])
                self.assertEqual(sum(row["InteriorNodeHistogramAlongEdge"]),
                                 row["InteriorNodesAwayFromCorners"])
            # The census is source-local: only the placement differs.
            self.assertEqual(rotated_census["RigidTransform"], ROTATE_Z)
            self.assertEqual(census["RigidTransform"], IDENTITY)
            for key in ("CornerIsotropicLongitudinalCurves", "LongitudinalCurves",
                        "Sqrt2IsotropicSize", "LongitudinalFaces"):
                self.assertEqual(census[key], rotated_census[key])
            # Pulling the rotated contract corners back leaves roundoff only.
            np.testing.assert_allclose(rotated_census["SemanticCorners"],
                                       census["SemanticCorners"], rtol=0.0, atol=1e-14)
            for row, rotated_row in zip(census["Corners"], rotated_census["Corners"]):
                self.assertEqual(set(row), set(rotated_row))
                for key, value in row.items():
                    if isinstance(value, int):
                        self.assertEqual(value, rotated_row[key], key)
                    else:
                        np.testing.assert_allclose(rotated_row[key], value, rtol=1e-9,
                                                   atol=1e-13, err_msg=key)
            metadata = json.loads((root / "identity.msh.metadata.json").read_text())
            self.assertEqual(metadata["CornerIsotropyRadius"], 0.4)
            self.assertEqual(metadata["SemanticCornerCount"], len(contract["SemanticCorners"]))

    def test_longitudinal_face_census_detects_misaligned_ridge_rows(self):
        """Julia unit tests: misaligned ridge rows leave full-height triangles the census
        counts; the grid-preserving corner law leaves none and keeps the interior row."""
        result = subprocess.run(
            [self.julia, "--startup-file=no", f"--project={self.project}",
             str(HERE / "test_longitudinal_face_census.jl")],
            cwd=REPO, capture_output=True, text=True, check=False, timeout=600)
        self.assertEqual(result.returncode, 0,
                         f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    def test_footprint_simplification_julia_unit_tests(self):
        """Julia unit tests: near-collinear duplicate footprint vertices are merged within
        the shared tolerance, a genuine 5-degree bend is kept, and the simplification is
        covariant under in-plane rigid motion."""
        result = subprocess.run(
            [self.julia, "--startup-file=no", f"--project={self.project}",
             str(HERE / "test_footprint_simplification.jl")],
            cwd=REPO, capture_output=True, text=True, check=False, timeout=600)
        self.assertEqual(result.returncode, 0,
                         f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")

    def test_footprint_collinear_tolerance_is_the_shared_coplanar_tolerance(self):
        from edge_volume_metric import COPLANAR_TOLERANCE
        source = MESHER.read_text()
        match = re.search(r"^const FOOTPRINT_COLLINEAR_TOLERANCE = ([0-9.eE+-]+)$", source,
                          re.MULTILINE)
        self.assertIsNotNone(match)
        self.assertEqual(float(match.group(1)), COPLANAR_TOLERANCE)

    def test_device_and_default_footprints_are_simplified_and_recorded(self):
        """The bound device footprint's near-collinear vertices (retained-etch loop 2,
        vertices 8-10) are merged before CAD face creation and recorded; the producer
        default collars are recorded with nothing to merge."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            etch = SOURCE / "retained-etch.csv"
            self.produce(root, "device", IDENTITY,
                         corner_isotropy=self.corner_isotropy(root, "device") +
                         ["--etch-boundary", str(etch)])
            census = json.loads((root / "device-census.json").read_text())
            from edge_volume_metric import COPLANAR_TOLERANCE
            self.assertEqual(census["FootprintCollinearTolerance"], COPLANAR_TOLERANCE)
            self.assertEqual(census["EtchBoundarySHA256"],
                             hashlib.sha256(etch.read_bytes()).hexdigest())
            polygons = census["FootprintPolygons"]
            self.assertEqual(len(polygons), 2)
            self.assertEqual(census["FootprintSimplification"]["Polygons"], 2)
            with etch.open(newline="") as stream:
                loops = {}
                for row in csv.DictReader(stream):
                    loops.setdefault(int(row["Loop"]), []).append([float(row["X"]), float(row["Y"])])
            by_size = sorted(polygons, key=lambda polygon: polygon["Simplification"]["OriginalVertices"])
            second = by_size[0]
            self.assertEqual(second["Simplification"]["OriginalVertices"], len(loops[2]))
            self.assertEqual(second["Simplification"]["RemovedVertexIndices"], [8, 9, 10])
            self.assertEqual(second["Points"], [point for index, point in enumerate(loops[2], 1)
                                                if index not in (8, 9, 10)])
            self.assertEqual(by_size[1]["Simplification"]["OriginalVertices"], len(loops[1]))
            for polygon in polygons:
                record = polygon["Simplification"]
                self.assertLessEqual(record["MaximumRelativeDeviation"], COPLANAR_TOLERANCE)
                self.assertLessEqual(record["MaximumDeviation"],
                                     COPLANAR_TOLERANCE * record["MaximumDeviationLocalScale"])
                self.assertEqual(record["Tolerance"], COPLANAR_TOLERANCE)
                self.assertFalse(polygon["Hole"])
            self.assertEqual(census["FootprintSimplification"]["RemovedVertices"],
                             sum(p["Simplification"]["RemovedVertexCount"] for p in polygons))
            # The producer default: one collar polygon per conductor loop, nothing merged.
            self.produce(root, "default", IDENTITY,
                         corner_isotropy=self.corner_isotropy(root, "default"))
            default = json.loads((root / "default-census.json").read_text())
            self.assertEqual(default["EtchBoundary"], "producer-default")
            self.assertGreater(len(default["FootprintPolygons"]), 0)
            self.assertEqual(default["FootprintSimplification"]["RemovedVertices"], 0)
            for polygon in default["FootprintPolygons"]:
                self.assertEqual(polygon["Simplification"]["RemovedVertexCount"], 0)
                self.assertGreaterEqual(len(polygon["Points"]), 3)

    def test_seed_corner_isotropy_fails_closed_on_placement_or_missing_options(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # Contract placed by the rotation but seed generated at identity.
            self.produce(root, "misplaced", IDENTITY,
                         corner_isotropy=self.corner_isotropy(root, "misplaced", ROTATE_Z),
                         expect_failure="differs from the seed rigid transform")
            self.produce(root, "no-census",
                         corner_isotropy=self.corner_isotropy(root, "no-census", census=False),
                         expect_failure="together")
            self.produce(root, "no-radius",
                         corner_isotropy=self.corner_isotropy(root, "no-radius", radius="0"),
                         expect_failure="together")
            absent = json.loads((SOURCE / "semantic-contract.json").read_text())
            absent["SemanticCorners"] = [[0.5, 0.5, 0.0]]
            path = root / "absent-semantic.json"; path.write_text(json.dumps(absent))
            self.produce(root, "absent", corner_isotropy=[
                "--semantic-contract", str(path), "--corner-isotropy-radius", "0.4",
                "--corner-census", str(root / "absent-census.json")],
                expect_failure="absent from the seed CAD")

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

    def test_multislot_census_interface_areas_are_the_written_slot_labels(self):
        """The census measures the labels the seed is written with: after the multi-slot
        postprocess these are the contract's slot/conductor labels (not the CAD faces),
        and their areas equal the ownership report's whole-element areas."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contract = json.loads(MULTISLOT["semantic"].read_text())
            path = root / "multislot-semantic.json"; path.write_text(json.dumps(contract))
            self.produce(root, "multislot-census", source=MULTISLOT, ownership=True,
                         corner_isotropy=["--semantic-contract", str(path),
                                          "--corner-isotropy-radius", "0.4",
                                          "--corner-census", str(root / "multislot-census.json")])
            census = json.loads((root / "multislot-census.json").read_text())
            rows = {row["Attribute"]: row for row in census["InterfaceAreas"]}
            self.assertEqual(len(rows), len(census["InterfaceAreas"]))
            self.assertEqual(set(rows), {item["Attribute"] for item in contract["BoundaryLabels"]})
            with (root / "multislot-census-ownership.csv").open(newline="") as stream:
                report = {int(row["attribute"]): row for row in csv.DictReader(stream)}
            for attribute, row in report.items():
                self.assertEqual(rows[attribute]["Triangles"], int(row["elements"]))
                self.assertAlmostEqual(rows[attribute]["Area"], float(row["area"]),
                                       delta=1e-9 * float(row["area"]))

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
