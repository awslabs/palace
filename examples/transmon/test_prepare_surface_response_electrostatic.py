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
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "surface_response_identification"))
import developer_form  # noqa: E402
EXAMPLE = Path(__file__).with_name("transmon_surface_amr.json")
# The example has no EdgeRefinement block since its tube default changed (S4-a item (6)), so the
# test writes its own --source: the example with one EdgeRefinement block per dielectric, which
# the script must strip unless --edge-refinement is given (review r2157 m-2).
SOURCE_EDGE_REFINEMENT = {"Radius": 0.2, "ElementsPerRadius": 3, "OuterRadiusFactor": 2.0, "CoreIndicatorWeight": 0.0}


def write_source(directory):
    config = json.loads(EXAMPLE.read_text())
    for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        dielectric["EdgeRefinement"] = dict(SOURCE_EDGE_REFINEMENT)
    source = Path(directory, "source.json")
    source.write_text(json.dumps(config, indent=2) + "\n")
    return source


def prepare(*extra):
    with tempfile.TemporaryDirectory() as tmp:
        library = Path(tmp, "process-library.json")
        library.write_text("{}\n")
        output = Path(tmp, "config.json")
        source = write_source(tmp)
        command = [sys.executable, str(SCRIPT), "--mesh", str(Path(tmp, "device.msh2")), "--library", str(library),
                   "--output", str(output), "--postpro", str(Path(tmp, "postpro")), "--source", str(source), *extra]
        result = subprocess.run(command, capture_output=True, text=True)
        config = json.loads(output.read_text()) if result.returncode == 0 else None
        return result, config


class PrepareSurfaceResponseElectrostaticTest(unittest.TestCase):
    def test_the_source_carries_an_edge_refinement_block_per_dielectric(self):
        with tempfile.TemporaryDirectory() as tmp:
            dielectrics = json.loads(write_source(tmp).read_text())["Boundaries"]["Postprocessing"]["Dielectric"]
        self.assertGreater(len(dielectrics), 0)
        for dielectric in dielectrics:
            self.assertEqual(dielectric["EdgeRefinement"], SOURCE_EDGE_REFINEMENT)

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
        # The v1 form (the default; the v1 cut): post-processing only, no parked key.
        self.assertEqual(correction["CorrectionMode"], "PostprocessOnly")
        self.assertEqual(developer_form.parked_config_settings(config), [])
        self.assertNotIn("Developer", correction)
        self.assertIn("the v1 form", result.stderr)

    def test_developer_form_is_the_record(self):
        # The record: the prism-library device run transmon_final_library_device_20260830/
        # config-p5-exact-amr11-full-mortar-over2-4r8g.json (cluster /data/home/simlap; 2026-09-08), whose block is
        # (in its own key order) Library, TargetInterfaces [1, 2, 3], UnmatchedPolicy Error, CorrectionMode Both,
        # SolveTol 1e-06, TraceCoupling SurfaceMortar, MortarOversampling 2, TranslationalDomainCorrection FixedTrace.
        # This writer (written after that run) orders the keys its own way: the recorded developer form carries the
        # record's settings and "Developer": true last.
        recorded = {"TargetInterfaces": [1, 2, 3], "UnmatchedPolicy": "Error", "CorrectionMode": "Both", "SolveTol": 1e-06,
                    "TraceCoupling": "SurfaceMortar", "MortarOversampling": 2, "TranslationalDomainCorrection": "FixedTrace"}
        result, config = prepare("--developer", "--trace-coupling", "SurfaceMortar")
        self.assertEqual(result.returncode, 0, result.stderr)
        correction = config["Solver"]["Electrostatic"]["ResponseCorrection"]
        self.assertEqual({k: v for k, v in correction.items() if k not in ("Library", "Developer")}, recorded)
        self.assertEqual(list(correction)[-1], "Developer")
        self.assertIs(correction["Developer"], True)
        self.assertIn("the recorded developer form", result.stderr)

    def test_parked_flags_need_the_developer_switch(self):
        for mode in ("SelfConsistent", "Both"):
            result, _ = prepare("--correction-mode", mode)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("--developer", result.stderr)
            result, config = prepare("--correction-mode", mode, "--developer")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(config["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"], mode)

    def test_amr_writes_the_recorded_prism_era_refinement_block(self):
        result, config = prepare("--order", "5", "--amr-max-its", "10", "--amr-max-size", "120000000", "--save-adapt-mesh",
                                 "--trace-coupling", "SurfaceMortar", "--mortar-oversampling", "2")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(config["Model"]["Refinement"], {
            "MaxIts": 10, "Tol": 0.01, "UpdateFraction": 0.7, "MaximumImbalance": 1.15, "MaxNCLevels": 8, "Nonconformal": True,
            "SaveAdaptIterations": True, "SaveAdaptMesh": True, "MaxSize": 120000000, "MaxSizePredicted": True, "UniformLevels": 0,
            "SerialUniformLevels": 0})
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
