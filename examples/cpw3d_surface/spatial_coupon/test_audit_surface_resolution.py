# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


BINARY = os.environ.get("AUDIT_SURFACE_RESOLUTION_BIN")


@unittest.skipUnless(BINARY, "Set AUDIT_SURFACE_RESOLUTION_BIN to the compiled auditor")
class SurfaceResolutionTest(unittest.TestCase):
    def audit(self, height, attributes=(5001, 6001)):
        # Only the unit-length edge shared by faces 5 and 6 is physical; the
        # remaining edges touch matching attribute 1 and must be excluded.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mesh = root / "one.mesh"
            mesh.write_text(
                "MFEM mesh v1.0\n\ndimension\n3\nelements\n1\n1 4 0 1 2 3\n"
                f"boundary\n4\n{attributes[0]} 2 0 2 1\n{attributes[1]} 2 0 1 3\n"
                "1 2 0 3 2\n1 2 1 2 3\nvertices\n4\n3\n"
                f"0 0 0\n1 0 0\n0.5 {height} 0\n0.5 0 {height}\n"
            )
            out = root / "audit.json"
            subprocess.run([BINARY, str(mesh), str(out)], check=True, timeout=15)
            return json.loads(out.read_text())

    def test_known_normal_and_tangential_lengths(self):
        result = self.audit(0.002)
        self.assertEqual(result["PhysicalMeshEdges"], 1)
        for family in ("5", "6"):
            row = result["Families"][family]
            self.assertEqual(row["Samples"], 1)
            self.assertAlmostEqual(row["AltitudePercentiles"]["50.000000"], 0.002)
            self.assertAlmostEqual(row["TangentEdgeLengthPercentiles"]["50.000000"], 1)
            self.assertEqual(row["IncidentEdgeLengthFractionWithAltitudeAtMost2nm"], 1)

    def test_extreme_aspect_ratio_does_not_round_altitude_to_zero(self):
        result = self.audit(1e-12)
        for family in ("5", "6"):
            self.assertAlmostEqual(
                result["Families"][family]["AltitudePercentiles"]["50.000000"] / 1e-12,
                1,
            )

    def test_matching_boundary_exclusion(self):
        result = self.audit(0.002, (1, 6001))
        self.assertEqual(result["PhysicalMeshEdges"], 0)
        self.assertEqual(result["Families"], {})


if __name__ == "__main__":
    unittest.main()
