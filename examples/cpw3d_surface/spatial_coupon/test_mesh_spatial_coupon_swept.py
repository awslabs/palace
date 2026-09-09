#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Executable invariants for straight and masked spatial swept meshers."""

import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
MESHER = ROOT / "mesh_spatial_coupon_swept.jl"
SIGNATURE = ROOT / "testdata" / "straight-edge-signature.csv"
TWO_EDGE_SIGNATURE = ROOT / "testdata" / "two-edge-transition-signature.csv"
TWO_EDGE_MASK = ROOT / "testdata" / "two-edge-transition-mask.csv"
TWO_EDGE_BOUNDARY = ROOT / "testdata" / "two-edge-transition-boundary.csv"
MULTISLOT_SIGNATURE = ROOT / "testdata" / "two-edge-multislot-signature.csv"
MULTISLOT_MASK = ROOT / "testdata" / "two-edge-multislot-mask.csv"
MULTISLOT_BOUNDARY = ROOT / "testdata" / "two-edge-multislot-boundary.csv"
SIX_EDGE_SIGNATURE = ROOT / "testdata" / "six-edge-cluster-signature.csv"
SIX_EDGE_MASK = ROOT / "testdata" / "six-edge-cluster-mask.csv"
SIX_EDGE_BOUNDARY = ROOT / "testdata" / "six-edge-cluster-boundary.csv"
MASKED_ONE_EDGE_SIGNATURE = ROOT / "testdata" / "one-edge-masked-signature.csv"
MASKED_ONE_EDGE_MASK = ROOT / "testdata" / "one-edge-masked-mask.csv"
MASKED_ONE_EDGE_BOUNDARY = ROOT / "testdata" / "one-edge-masked-boundary.csv"


class SweptSpatialCouponMesherTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.julia = shutil.which(os.environ.get("JULIA", "julia"))
        if cls.julia is None:
            raise unittest.SkipTest("Julia is not available")
        requested = os.environ.get("PALACE_JULIA_PROJECT")
        candidates = [Path(requested)] if requested else [
            REPO / "examples",
            REPO / "test" / "examples",
        ]
        cls.project = None
        for candidate in candidates:
            result = subprocess.run(
                [
                    cls.julia,
                    f"--project={candidate}",
                    "-e",
                    "import Gmsh: gmsh; gmsh.initialize(); gmsh.finalize()",
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode == 0:
                cls.project = candidate
                break
        if cls.project is None:
            raise unittest.SkipTest(
                "No candidate Julia project has an instantiated Gmsh dependency"
            )

    def run_mesher(
        self,
        root,
        name,
        *,
        kind="thin",
        normal=0.2,
        tangent=1.0,
        max_nodes=500_000,
        max_elements=2_000_000,
        check=True,
    ):
        mesh = root / f"{name}.msh"
        command = [
            self.julia,
            f"--project={self.project}",
            str(MESHER),
            str(SIGNATURE),
            kind,
            str(mesh),
            "--radius",
            "1",
            "--metal-thickness",
            "0.2",
            "--overetch",
            "0.1",
            "--sidewall-angle",
            "80",
            "--top-radius",
            "0.02",
            "--bottom-radius",
            "0.02",
            "--lc-normal",
            str(normal),
            "--lc-tangent",
            str(tangent),
            "--lc-far",
            "0.4",
            "--mesh-order",
            "2",
            "--max-nodes",
            str(max_nodes),
            "--max-elements",
            str(max_elements),
        ]
        result = subprocess.run(
            command,
            cwd=REPO,
            capture_output=True,
            text=True,
            check=False,
        )
        if check and result.returncode:
            self.fail(
                f"Mesher failed with exit {result.returncode}:\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        metadata = None
        metadata_path = Path(str(mesh) + ".metadata.json")
        if metadata_path.is_file():
            metadata = json.loads(metadata_path.read_text())
        return result, mesh, metadata

    def run_edge_cluster_mesher(
        self,
        root,
        name,
        *,
        kind,
        normal,
        tangent=0.5,
        signature=TWO_EDGE_SIGNATURE,
        mask=TWO_EDGE_MASK,
        boundary=TWO_EDGE_BOUNDARY,
    ):
        mesh = root / f"{name}.msh"
        command = [
            self.julia,
            f"--project={self.project}",
            str(MESHER),
            str(signature),
            kind,
            str(mesh),
            "--mask",
            str(mask),
            "--boundary",
            str(boundary),
            "--radius",
            "2",
            "--metal-thickness",
            "0.1",
            "--overetch",
            "0.05",
            "--sidewall-angle",
            "90",
            "--top-radius",
            "0",
            "--bottom-radius",
            "0",
            "--lc-normal",
            str(normal),
            "--lc-tangent",
            str(tangent),
            "--lc-far",
            "0.5",
            "--process-core-width",
            "0.4",
            "--normal-growth-ratio",
            "1.4",
            "--mesh-order",
            "2",
            "--max-nodes",
            "2000000",
            "--max-elements",
            "4000000",
        ]
        result = subprocess.run(
            command,
            cwd=REPO,
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode:
            self.fail(
                f"Edge-cluster mesher failed with exit {result.returncode}:\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        metadata = json.loads(Path(str(mesh) + ".metadata.json").read_text())
        return mesh, metadata

    def assert_common_metadata(self, metadata, *, kind, normal, tangent):
        self.assertEqual(metadata["Version"], 2)
        self.assertEqual(metadata["MeshingMode"], "Swept")
        self.assertEqual(metadata["SweptTopology"], "StraightEdgePrisms")
        self.assertEqual(metadata["Fabricated"], kind == "fabricated")
        self.assertEqual(
            metadata["MetalSurfacePartition"],
            "InterfaceSlotAndConductor",
        )
        self.assertEqual(metadata["InputEdgeCount"], 1)
        self.assertEqual(metadata["FineSize"], normal)
        self.assertEqual(metadata["NormalSize"], normal)
        self.assertEqual(metadata["TangentialSize"], tangent)
        self.assertEqual(metadata["VolumeAttributes"], [1, 2])
        self.assertEqual(metadata["NonmanifoldFaceCount"], 0)
        self.assertEqual(metadata["TransitionVolumeElementCount"], 0)
        self.assertEqual(
            metadata["SweptVolumeElementCount"],
            metadata["VolumeElementCount"],
        )
        self.assertGreater(metadata["NodeCount"], 0)
        self.assertGreater(metadata["VolumeElementCount"], 0)
        self.assertGreater(metadata["BoundaryFaceCount"], 0)
        self.assertGreater(metadata["InteriorFaceCount"], 0)
        self.assertGreater(metadata["MinimumScaledJacobian"], 0.0)
        self.assertTrue(
            all(name.startswith("Prism") for name in metadata["VolumeElementTypes"])
        )
        self.assertLessEqual(
            metadata["MeasuredMaximumNearFeatureNormalSpacing"],
            1.05 * normal,
        )
        self.assertLessEqual(
            metadata["MeasuredMaximumLongitudinalSpacing"],
            1.05 * tangent,
        )

    def test_normal_and_tangent_coordinates_are_independent(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, _, base = self.run_mesher(root, "n20-t100")
            _, _, normal_refined = self.run_mesher(
                root, "n10-t100", normal=0.1
            )
            _, _, tangent_refined = self.run_mesher(
                root, "n20-t050", tangent=0.5
            )

            self.assert_common_metadata(
                base, kind="thin", normal=0.2, tangent=1.0
            )
            self.assert_common_metadata(
                normal_refined, kind="thin", normal=0.1, tangent=1.0
            )
            self.assert_common_metadata(
                tangent_refined, kind="thin", normal=0.2, tangent=0.5
            )
            self.assertEqual(
                base["LongitudinalCoordinates"],
                normal_refined["LongitudinalCoordinates"],
            )
            self.assertEqual(
                base["MeasuredNearFeatureNormalSpacings"],
                tangent_refined["MeasuredNearFeatureNormalSpacings"],
            )
            self.assertGreater(
                normal_refined["CrossSectionNodeCount"],
                base["CrossSectionNodeCount"],
            )
            self.assertGreater(
                tangent_refined["LongitudinalStationCount"],
                base["LongitudinalStationCount"],
            )

    def test_fabricated_attributes_rounding_and_quality(self):
        with tempfile.TemporaryDirectory() as directory:
            _, _, metadata = self.run_mesher(
                Path(directory),
                "fabricated",
                kind="fabricated",
                normal=0.1,
            )
            self.assert_common_metadata(
                metadata,
                kind="fabricated",
                normal=0.1,
                tangent=1.0,
            )
            self.assertEqual(
                metadata["SurfaceAttributes"],
                [1, 3100, 5001, 6001],
            )

    def test_masked_edge_cluster_transition_mesh_refines_subquadratically(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, thin = self.run_edge_cluster_mesher(
                root,
                "thin-n50",
                kind="thin",
                normal=0.05,
            )
            _, fabricated = self.run_edge_cluster_mesher(
                root,
                "fabricated-n50",
                kind="fabricated",
                normal=0.05,
            )
            _, refined = self.run_edge_cluster_mesher(
                root,
                "fabricated-n25",
                kind="fabricated",
                normal=0.025,
            )
            _, tangent_refined = self.run_edge_cluster_mesher(
                root,
                "fabricated-tangent-refined",
                kind="fabricated",
                normal=0.05,
                tangent=0.25,
            )
            _, multislot = self.run_edge_cluster_mesher(
                root,
                "fabricated-multislot",
                kind="fabricated",
                normal=0.05,
                signature=MULTISLOT_SIGNATURE,
                mask=MULTISLOT_MASK,
                boundary=MULTISLOT_BOUNDARY,
            )

            _, masked_one_edge = self.run_edge_cluster_mesher(
                root,
                "thin-masked-one-edge",
                kind="thin",
                normal=0.05,
                signature=MASKED_ONE_EDGE_SIGNATURE,
                mask=MASKED_ONE_EDGE_MASK,
                boundary=MASKED_ONE_EDGE_BOUNDARY,
            )
            _, six_edge = self.run_edge_cluster_mesher(
                root,
                "fabricated-six-edge",
                kind="fabricated",
                normal=0.05,
                signature=SIX_EDGE_SIGNATURE,
                mask=SIX_EDGE_MASK,
                boundary=SIX_EDGE_BOUNDARY,
            )

        for metadata in (thin, fabricated, refined, tangent_refined, multislot):
            self.assertEqual(metadata["Version"], 2)
            self.assertEqual(metadata["MeshingMode"], "ExplicitEdgeClusterTransition")
            self.assertEqual(metadata["InputEdgeCount"], 2)
            self.assertEqual(metadata["VolumeAttributes"], [1, 2])
            self.assertEqual(metadata["VolumeElementTypes"], ["Prism 18"])
            self.assertGreater(metadata["MinimumScaledJacobian"], 0.0)
            self.assertLessEqual(
                metadata["MeasuredMaximumLongitudinalSpacing"],
                0.5,
            )
        self.assertEqual(masked_one_edge["MeshingMode"], "ExplicitEdgeClusterTransition")
        self.assertEqual(masked_one_edge["InputEdgeCount"], 1)
        self.assertEqual(masked_one_edge["SurfaceAttributes"], [1, 3000, 4001])
        self.assertEqual(six_edge["MeshingMode"], "ExplicitEdgeClusterTransition")
        self.assertEqual(six_edge["InputEdgeCount"], 6)
        self.assertEqual(six_edge["UnmatchedInternalPlanEdgeCount"], 0)
        self.assertEqual(six_edge["NonmanifoldFaceCount"], 0)
        self.assertEqual(
            six_edge["SurfaceAttributes"],
            [1, 3100, 3101, 5001, 5002, 5101, 5102, 6001, 6002, 6101, 6102],
        )
        self.assertEqual(thin["SurfaceAttributes"], [1, 3000, 4001])
        self.assertEqual(
            fabricated["SurfaceAttributes"],
            [1, 3100, 5001, 6001],
        )
        self.assertEqual(
            multislot["SurfaceAttributes"],
            [1, 3001, 3100, 3101, 5001, 5101, 6001, 6101],
        )
        self.assertEqual(fabricated["NormalCoordinates"][1], 0.05)
        self.assertEqual(refined["NormalCoordinates"][1], 0.025)
        self.assertGreater(refined["NodeCount"], fabricated["NodeCount"])
        self.assertLess(refined["NodeCount"] / fabricated["NodeCount"], 2.0)
        self.assertEqual(
            fabricated["NormalCoordinates"],
            tangent_refined["NormalCoordinates"],
        )
        self.assertLessEqual(
            tangent_refined["MeasuredMaximumLongitudinalSpacing"],
            0.25,
        )
        self.assertGreater(tangent_refined["NodeCount"], fabricated["NodeCount"])

    def test_budget_failure_removes_partial_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            result, mesh, metadata = self.run_mesher(
                Path(directory),
                "over-budget",
                normal=0.2,
                tangent=1.0,
                max_elements=1,
                check=False,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exceeds element budget", result.stderr)
            self.assertFalse(mesh.exists())
            self.assertIsNone(metadata)


if __name__ == "__main__":
    unittest.main()
