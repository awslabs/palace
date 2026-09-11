# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import tomllib
import unittest

ROOT = Path(__file__).resolve().parent


class EtchFootprintTest(unittest.TestCase):
    def test_retained_floor_recreates_the_same_cad_geometry(self):
        julia = shutil.which("julia")
        if not julia:
            self.skipTest("Julia not available")
        project = Path(os.environ.get("PALACE_JULIA_PROJECT", ROOT.parents[2] / "test/examples"))
        check = subprocess.run([julia, f"--project={project}", "-e", "using Gmsh"],
                               capture_output=True, text=True)
        if check.returncode:
            self.skipTest("No instantiated Julia Gmsh project")
        environment = dict(os.environ, TET_GEOMETRY_ORDER="1", JULIA_NUM_THREADS="1",
                           OPENBLAS_NUM_THREADS="1")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for source, target in (("signature", "mesh-signature"), ("mask", "plan-view-mask"),
                                   ("boundary", "plan-view-boundary")):
                shutil.copy2(ROOT / f"testdata/one-edge-masked-{source}.csv", root / f"{target}.csv")
            process = root / "process.toml"
            process.write_text('Units="um"\nRadius=0.5\nMetalThickness=0.1\nOveretch=0.05\n'
                               'SidewallAngle=90.0\nTopRounding=0.0\nTrenchRounding=0.0\n')
            base = [julia, "--startup-file=no", f"--project={project}"]
            mesher = [*base, str(ROOT / "mesh_graded_tet_experiment.jl"), str(root), "fabricated"]
            for name, extra in (("cad", ["--geometry-only"]),
                                ("full", ["--reference-measures", str(root / "cad.msh.cad-measures.csv")])):
                subprocess.run([*mesher, str(root / f"{name}.msh"), "1", ".05", ".4",
                                "--process", str(process), *extra],
                               env=environment, check=True, capture_output=True, text=True)
            footprint = root / "footprint.csv"
            subprocess.run([*base, str(ROOT / "export_etch_footprint.jl"), str(root / "full.msh"),
                            "-.05", "0", str(footprint)], env=environment,
                           check=True, capture_output=True, text=True)
            subprocess.run([*mesher, str(root / "reconstructed.msh"), "1", ".05", ".4",
                            "--process", str(process), "--etch-boundary", str(footprint),
                            "--geometry-only", "--reference-measures", str(root / "cad.msh.cad-measures.csv")],
                           env=environment, check=True, capture_output=True, text=True)
            report = tomllib.loads(Path(str(footprint) + ".provenance.toml").read_text())
            self.assertGreater(report["Vertices"], 2)
            self.assertAlmostEqual(report["FloorArea"], report["LoopArea"], places=10)
            self.assertTrue((root / "reconstructed.msh.geometry-gate.txt").is_file())


if __name__ == "__main__":
    unittest.main()
