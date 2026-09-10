#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import os
import shutil
import subprocess
import tempfile
import textwrap
import tomllib
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[2]
MESHER = ROOT / "mesh_graded_tet_experiment.jl"
RELABEL = ROOT / "relabel_frozen_interface_mesh.jl"
TESTDATA = ROOT / "testdata"


class FrozenInterfaceRelabelTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.julia = shutil.which(os.environ.get("JULIA", "julia"))
        if cls.julia is None:
            raise unittest.SkipTest("Julia is not available")
        requested = os.environ.get("PALACE_JULIA_PROJECT")
        candidates = [Path(requested)] if requested else [REPO / "examples", REPO / "test" / "examples"]
        cls.project = None
        for candidate in candidates:
            result = subprocess.run(
                [cls.julia, f"--project={candidate}", "-e", "import Gmsh: gmsh; gmsh.initialize(); gmsh.finalize()"],
                capture_output=True,
                text=True,
                check=False,
            )
            if result.returncode == 0:
                cls.project = candidate
                break
        if cls.project is None:
            raise unittest.SkipTest("No instantiated Julia Gmsh project")

    def test_serialized_multislot_relabel_preserves_frozen_mesh(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            inputs = {
                "two-edge-multislot-signature.csv": "mesh-signature.csv",
                "two-edge-multislot-mask.csv": "plan-view-mask.csv",
                "two-edge-multislot-boundary.csv": "plan-view-boundary.csv",
            }
            for source, target in inputs.items():
                shutil.copy2(TESTDATA / source, root / target)
            (root / "process.toml").write_text(textwrap.dedent("""
                Units = "um"
                Radius = 0.5
                MetalThickness = 0.1
                Overetch = 0.05
                SidewallAngle = 90.0
                TopRounding = 0.0
                TrenchRounding = 0.0
            """))
            study = root / "study.toml"
            study.write_text(textwrap.dedent("""
                MaxElements = 500000
                MaxNodes = 200000
                OptimizeVolume = true
                [[Variants]]
                Name = "candidate"
                MinimumSize = 0.05
                NearGrowth = 1.0
                FarGrowth = 2.0
                TransitionDistance = 0.04
                MaximumSize = 0.4
            """))
            geometry = root / "geometry.msh"
            common = [
                self.julia,
                f"--project={self.project}",
                str(MESHER),
                str(root),
                "thin",
            ]
            subprocess.run(
                [*common, str(geometry), "1", "0.05", "0.4", "--process", str(root / "process.toml"), "--geometry-only"],
                check=True,
                capture_output=True,
                text=True,
            )
            family = root / "family.msh"
            environment = {**os.environ, "TET_GEOMETRY_ORDER": "1", "JULIA_NUM_THREADS": "1"}
            subprocess.run(
                [
                    *common,
                    str(family),
                    "1",
                    "0.05",
                    "0.4",
                    "--process",
                    str(root / "process.toml"),
                    "--volume-study",
                    str(study),
                    "--reference-measures",
                    str(geometry) + ".cad-measures.csv",
                ],
                check=True,
                capture_output=True,
                text=True,
                env=environment,
            )
            expected = root / "expected.csv"
            expected.write_text(
                "dimension,attribute,measure\n"
                "2,1,0\n2,3000,0\n2,3001,0\n2,4001,0\n2,4101,0\n"
                "3,1,0\n3,2,0\n"
            )
            output = root / "relabeled.msh"
            subprocess.run(
                [
                    self.julia,
                    f"--project={self.project}",
                    str(RELABEL),
                    str(root),
                    "thin",
                    str(root / "family-candidate.msh"),
                    str(output),
                    "--expected-measures",
                    str(expected),
                ],
                check=True,
                capture_output=True,
                text=True,
                env=environment,
            )
            with output.open("rb") as stream:
                self.assertEqual(stream.readline().strip(), b"$MeshFormat")
                self.assertEqual(stream.readline().strip(), b"2.2 1 8")
            metadata = tomllib.loads((Path(str(output) + ".relabel.toml")).read_text())
            self.assertEqual(metadata["ActualInterfaceAttributes"], [3000, 3001, 4001, 4101])
            self.assertTrue(metadata["InMemoryNodeAndElementIDsPreserved"])
            self.assertTrue(metadata["SerializedGeometryAndConnectivityPreserved"])
            self.assertTrue(metadata["SerializedRoundTripVerified"])
            certificate = tomllib.loads((Path(str(output) + ".partition-certificate.toml")).read_text())
            self.assertEqual(certificate["MeshSHA256"], metadata["OutputSHA256"])


if __name__ == "__main__":
    unittest.main()
