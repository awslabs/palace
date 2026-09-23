#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("prepare_surface_response_electrostatic.py")
SOURCE = Path(__file__).with_name("transmon_surface_amr.json")


def prepare(*extra):
    with tempfile.TemporaryDirectory() as tmp:
        library = Path(tmp, "process-library.json")
        library.write_text("{}\n")
        output = Path(tmp, "config.json")
        command = [sys.executable, str(SCRIPT), "--mesh", str(Path(tmp, "device.msh2")), "--library", str(library),
                   "--output", str(output), "--postpro", str(Path(tmp, "postpro")), "--source", str(SOURCE), *extra]
        result = subprocess.run(command, capture_output=True, text=True)
        config = json.loads(output.read_text()) if result.returncode == 0 else None
        return result, config


class PrepareSurfaceResponseElectrostaticTest(unittest.TestCase):
    def test_default_is_a_single_solve_without_the_source_edge_refinement(self):
        result, config = prepare("--order", "3")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(config["Model"]["Refinement"], {"MaxIts": 0, "UniformLevels": 0, "SerialUniformLevels": 0})
        for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
            self.assertEqual(dielectric["EdgeDistances"], [2.0])
            self.assertNotIn("EdgeRefinement", dielectric)
        self.assertEqual(config["Solver"]["Order"], 3)
        correction = config["Solver"]["Electrostatic"]["ResponseCorrection"]
        self.assertEqual((correction["TraceCoupling"], correction["MortarOversampling"]), ("Collocated", 2))

    def test_amr_writes_the_recorded_prism_era_refinement_block(self):
        result, config = prepare("--order", "5", "--amr-max-its", "10", "--amr-max-size", "120000000", "--save-adapt-mesh",
                                 "--trace-coupling", "SurfaceMortar", "--mortar-oversampling", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(config["Model"]["Refinement"], {
            "MaxIts": 10, "Tol": 0.01, "UpdateFraction": 0.7, "MaximumImbalance": 1.15, "MaxNCLevels": 8, "Nonconformal": True,
            "SaveAdaptIterations": True, "SaveAdaptMesh": True, "MaxSize": 120000000, "UniformLevels": 0, "SerialUniformLevels": 0})
        correction = config["Solver"]["Electrostatic"]["ResponseCorrection"]
        self.assertEqual((correction["TraceCoupling"], correction["MortarOversampling"]), ("SurfaceMortar", 2))
        for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
            self.assertNotIn("EdgeRefinement", dielectric)

    def test_edge_refinement_is_an_explicit_opt_in_with_its_radius_in_the_distances(self):
        result, config = prepare("--amr-max-its", "2", "--amr-tol", "0.05", "--edge-refinement", "0.2", "3")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(config["Model"]["Refinement"]["Tol"], 0.05)
        self.assertNotIn("MaxSize", config["Model"]["Refinement"])
        for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
            self.assertEqual(dielectric["EdgeDistances"], [0.2, 2.0])
            self.assertEqual(dielectric["EdgeRefinement"],
                             {"Radius": 0.2, "ElementsPerRadius": 3, "OuterRadiusFactor": 2.0, "CoreIndicatorWeight": 0.0})
        result, config = prepare("--edge-refinement", "2.0", "1")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeDistances"], [2.0])

    def test_rejects_inconsistent_amr_and_tube_requests(self):
        for extra in (["--save-adapt-mesh"], ["--amr-max-size", "5"], ["--amr-max-its", "-1"], ["--amr-max-its", "2", "--amr-tol", "0"],
                      ["--edge-refinement", "2.5", "3"], ["--edge-refinement", "0.2", "2.5"], ["--edge-refinement", "0.2", "0"]):
            result, config = prepare(*extra)
            self.assertEqual(result.returncode, 2, extra)
            self.assertIsNone(config)


if __name__ == "__main__":
    unittest.main()
