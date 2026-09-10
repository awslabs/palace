# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import unittest

BINARY = os.environ.get("AUDIT_MESH_MEASURES_BIN")


@unittest.skipUnless(BINARY, "Set AUDIT_MESH_MEASURES_BIN to the compiled auditor")
class MeshMeasuresTest(unittest.TestCase):
    def run_mesh(self, faces):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            mesh = root / "one.mesh"
            mesh.write_text(
                "MFEM mesh v1.0\ndimension\n3\nelements\n1\n1 4 0 1 2 3\n"
                f"boundary\n{len(faces)}\n"+"\n".join(faces)+
                "\nvertices\n4\n3\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n"
            )
            output = root / "result.json"
            run = subprocess.run([BINARY, str(mesh), str(output)], timeout=15,
                                 capture_output=True, text=True)
            return run.returncode, json.loads(output.read_text()) if output.exists() else None

    def test_affine_measure_and_closed_boundary(self):
        status, result = self.run_mesh(["1 2 0 2 1", "1 2 0 1 3", "1 2 0 3 2", "1 2 1 2 3"])
        self.assertEqual(status, 0)
        self.assertTrue(result["BoundaryCoverageChecked"])
        self.assertAlmostEqual(result["MaterialVolumes"]["1"], 1/6)
        self.assertAlmostEqual(result["BoundaryAreas"]["1"], 1.5+math.sqrt(3)/2)
        self.assertLess(result["MaximumRelativeQuadratureDifference"], 1e-12)

    def test_missing_exterior_face_is_rejected(self):
        status, _ = self.run_mesh(["1 2 0 2 1", "1 2 0 1 3", "1 2 0 3 2"])
        self.assertNotEqual(status, 0)


if __name__ == "__main__":
    unittest.main()
