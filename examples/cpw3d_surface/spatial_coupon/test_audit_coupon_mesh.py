#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Small analytic-prism tests for the standalone, read-only MFEM mesh audit.

Set PALACE_MESH_AUDIT_BINARY to the compiled audit_coupon_mesh executable.
"""
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


@unittest.skipUnless(os.environ.get("PALACE_MESH_AUDIT_BINARY"), "mesh audit binary not set")
class MeshAuditTest(unittest.TestCase):
    def test_perfect_and_stretched_prisms(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            boundary = root / "boundary.csv"
            boundary.write_text(
                "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
                "1,1,1,0,0,Physical,0,0\n"
                "1,2,1,0,0,Physical,1,0\n"
                "1,3,1,0,0,Physical,0.5,0.8660254037844386\n"
            )
            for height, expected in [(1.0, 1.0), (0.002, 500.0)]:
                with self.subTest(height=height):
                    points = [
                        (0, 0, 0), (1, 0, 0), (0.5, math.sqrt(3) / 2, 0),
                        (0, 0, height), (1, 0, height),
                        (0.5, math.sqrt(3) / 2, height),
                    ]
                    mesh = root / "test.msh"
                    mesh.write_text(
                        "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n6\n"
                        + "".join(
                            f"{i} {x} {y} {z}\n"
                            for i, (x, y, z) in enumerate(points, 1)
                        )
                        + "$EndNodes\n$Elements\n1\n"
                        "1 6 2 1 1 1 2 3 4 5 6\n$EndElements\n"
                    )
                    output = root / "report.json"
                    subprocess.run(
                        [os.environ["PALACE_MESH_AUDIT_BINARY"], str(mesh),
                         str(boundary), "0", str(output)],
                        check=True, timeout=15, capture_output=True, text=True,
                    )
                    report = json.loads(output.read_text())
                    self.assertEqual(report["Elements"], 1)
                    self.assertAlmostEqual(
                        report["KappaPercentiles"]["100.000000"] / expected, 1.0,
                        places=12,
                    )
                    self.assertAlmostEqual(
                        report["Plan"]["MinimumAnglePercentilesDegrees"]["0.000000"],
                        60.0, places=10,
                    )
                    self.assertEqual(report["Plan"]["UniqueTriangles"], 1)
                    self.assertEqual(report["NonpositiveCenterJacobians"], 0)
                    self.assertLess(report["MaximumVolumeDiscrepancy"], 1e-12)
                    self.assertEqual(
                        sum(r["Elements"] for r in report["DistanceJointHistogram"]), 1
                    )


if __name__ == "__main__":
    unittest.main()
