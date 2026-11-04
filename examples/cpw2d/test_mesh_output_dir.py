# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""OI-CPW2D-MESHER-WRITES-REPO-TREE (decision 1247; fix round P2): the cpw2d cross-section
meshers (mesh_surface_calibration.jl, mesh_thin.jl, mesh_thick.jl) take an `output_dir`
keyword; the default (the mesher's own directory, i.e. the repository tree) is unchanged, so
a caller that needs the mesh elsewhere passes the directory instead of moving the file after
the fact (the fix-round P0 narrow-gap D3 incident, STEP 25 of its REPORT)."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
MESH_DIR = HERE / "mesh"


class MesherOutputDirTest(unittest.TestCase):
    def test_every_cpw2d_mesher_declares_the_output_dir_keyword_with_the_unchanged_default(self):
        for name in ("mesh_surface_calibration.jl", "mesh_thin.jl", "mesh_thick.jl"):
            text = (MESH_DIR / name).read_text()
            self.assertIn("output_dir::String", text, name)
            self.assertIn("gmsh.write(joinpath(output_dir, filename))", text, name)
            # The default is the file's own directory (the record behaviour).
            self.assertRegex(text, r"output_dir::String\s*=\s*@__DIR__", name)
            self.assertNotIn("gmsh.write(joinpath(@__DIR__, filename))", text, name)

    def test_the_calibration_mesher_writes_into_the_requested_directory(self):
        julia = shutil.which(os.environ.get("JULIA", "julia"))
        if julia is None:
            self.skipTest("Julia is not available")
        probe = subprocess.run([julia, "--startup-file=no", "-e", "import Gmsh"], capture_output=True,
                               text=True, timeout=600)
        if probe.returncode != 0:
            self.skipTest("the Gmsh Julia package is not available: " + probe.stderr.strip()[-200:])
        before = sorted(os.listdir(MESH_DIR))
        with tempfile.TemporaryDirectory() as temp:
            # Julia string literals are double-quoted (a single-quoted literal is a Char).
            mesher = json.dumps(str(MESH_DIR / "mesh_surface_calibration.jl"))
            script = (f'include({mesher})\n'
                      f'generate_cpw2d_mesh(; w_ground=60.0, h_substrate=60.0, h_vacuum=60.0, lc_fine=0.5, '
                      f'lc_far=30.0, mesh_order=1, filename="unit-test-output-dir.msh", output_dir={json.dumps(temp)})\n')
            env = dict(os.environ, JULIA_NUM_THREADS="1")
            result = subprocess.run([julia, "--startup-file=no", "-t", "1", "-e", script], capture_output=True,
                                    text=True, timeout=900, env=env)
            self.assertEqual(result.returncode, 0, result.stderr[-2000:])
            written = Path(temp) / "unit-test-output-dir.msh"
            self.assertTrue(written.is_file(), result.stdout[-2000:])
            self.assertGreater(written.stat().st_size, 0)
        # Nothing was written into the repository tree.
        self.assertEqual(sorted(os.listdir(MESH_DIR)), before)


if __name__ == "__main__":
    unittest.main()
