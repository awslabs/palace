#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE = Path(__file__).with_name("run_corrected_eigenmode.py")
SPEC = importlib.util.spec_from_file_location("corrected_benchmark", MODULE)
BENCHMARK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCHMARK)


class CorrectedEigenmodeBenchmarkTest(unittest.TestCase):
    def test_checked_in_baseline_and_library(self):
        baseline = json.loads(BENCHMARK.DEFAULT_BASELINE.read_text())
        self.assertEqual(baseline["SchemaVersion"], 1)
        self.assertEqual(baseline["Benchmark"], BENCHMARK.NAME)
        self.assertFalse(baseline["Qualified"])
        self.assertEqual(baseline["FEMOrder"], 1)
        self.assertEqual(sorted(baseline["Ranks"]), ["1", "2", "4"])
        for rank in baseline["Ranks"].values():
            self.assertEqual(rank["Workload"]["Patches"], 6264)
            self.assertEqual(rank["Workload"]["TraceCoefficients"], 69212)
            self.assertEqual(rank["Workload"]["PointQueries"], 378597)
        BENCHMARK.verify_fixtures(
            baseline["Fixtures"],
            BENCHMARK.TRANSMON / "transmon_surface_coarse.json",
            BENCHMARK.TRANSMON / "mesh/transmon_surface_p1.msh2",
            BENCHMARK.DEFAULT_LIBRARY,
        )
        library = json.loads(BENCHMARK.DEFAULT_LIBRARY.read_text())
        self.assertTrue(
            all(
                model["BoundaryLawQualification"]["Status"] == "Unqualified"
                for model in library["Models"]
            )
        )
        isolated_basis = (
            BENCHMARK.DEFAULT_LIBRARY.parent
            / "models/001-isolated-edge/basis-points.csv"
        )
        strip_basis = (
            BENCHMARK.DEFAULT_LIBRARY.parent
            / "models/002-same-conductor-strip-2um/basis-points.csv"
        )
        self.assertEqual(len(isolated_basis.read_text().splitlines()) - 1, 11)
        self.assertEqual(len(strip_basis.read_text().splitlines()) - 1, 12)

    def test_snapshot_comparison(self):
        expected = {"x.csv": [{"m": 1.0, "value": 2.0, "tiny": 1.0e-16}]}
        observed = {"x.csv": [{"m": 1.0, "value": 2.0001, "tiny": 1.0001e-16}]}
        self.assertEqual(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20), []
        )
        observed["x.csv"][0]["m"] = 2.0
        self.assertTrue(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20)
        )
        # An appended coverage column is not compared; any other new column is a change.
        observed["x.csv"][0]["m"] = 1.0
        observed["x.csv"][0]["coverage status"] = "OK"
        self.assertEqual(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20), []
        )
        observed["x.csv"][0]["another column"] = 1.0
        self.assertTrue(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20)
        )
        del observed["x.csv"][0]["another column"]
        del observed["x.csv"][0]["value"]
        self.assertTrue(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20)
        )

    def test_prepare_config_is_p1_and_fixed_mesh(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary = Path(temporary)
            output = temporary / "config.json"
            mesh = BENCHMARK.prepare_config(
                BENCHMARK.TRANSMON / "transmon_surface_coarse.json",
                BENCHMARK.DEFAULT_LIBRARY,
                output,
                temporary / "postpro",
            )
            config = json.loads(output.read_text())
            self.assertTrue(mesh.is_file())
            self.assertEqual(config["Solver"]["Order"], 1)
            self.assertEqual(config["Model"]["Refinement"]["MaxIts"], 0)
            self.assertEqual(config["Model"]["Refinement"]["UniformLevels"], 0)
            self.assertEqual(config["Solver"]["Eigenmode"]["N"], 1)
            self.assertFalse(config["Problem"]["OutputFormats"]["Paraview"])
            for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
                self.assertNotIn("EdgeRefinement", dielectric)
                self.assertEqual(dielectric["EdgeDistances"], [2.0])
                self.assertFalse(dielectric["LocalizeEdgeEnergy"])

    def prepare(self, temporary, developer):
        output = Path(temporary) / ("developer.json" if developer else "v1.json")
        BENCHMARK.prepare_config(
            BENCHMARK.TRANSMON / "transmon_surface_coarse.json",
            BENCHMARK.DEFAULT_LIBRARY,
            output,
            Path(temporary) / "postpro",
            developer,
        )
        return json.loads(output.read_text())

    def test_default_config_is_the_v1_form(self):
        # The v1 cut (developer_form.py): post-processing only, no parked key.
        with tempfile.TemporaryDirectory() as temporary:
            config = self.prepare(temporary, False)
        correction = config["Solver"]["SurfaceResponseCorrection"]
        self.assertEqual(list(correction), ["Library", "TargetInterfaces", "UnmatchedPolicy"])
        self.assertEqual(BENCHMARK.developer_form.parked_config_settings(config), [])

    def test_developer_form_is_the_record(self):
        # The record: the block the baseline's runs were made with (prepare_config of e4195bac6c, this file before
        # the cut: Library, TargetInterfaces [1, 2, 3], UnmatchedPolicy Warn, SolveTol 1e-6 under the default of
        # record CorrectionMode Both). The recorded developer form is that block with the default stated and
        # "Developer": true; every other key of the config is the v1 form's.
        with tempfile.TemporaryDirectory() as temporary:
            developer = self.prepare(temporary, True)
            v1 = self.prepare(temporary, False)
        correction = developer["Solver"]["SurfaceResponseCorrection"]
        self.assertEqual(
            list(correction.items()),
            [
                ("Library", str(BENCHMARK.DEFAULT_LIBRARY.resolve())),
                ("TargetInterfaces", [1, 2, 3]),
                ("UnmatchedPolicy", "Warn"),
                ("SolveTol", 1.0e-6),
                ("CorrectionMode", "Both"),
                ("Developer", True),
            ],
        )
        del developer["Solver"]["SurfaceResponseCorrection"]
        del v1["Solver"]["SurfaceResponseCorrection"]
        self.assertEqual(developer, v1)

    def test_v1_form_compares_the_kept_columns(self):
        expected = {"x.csv": [{"m": 1.0, "p_surf corrected": 2.0, "p_surf self-consistent": 3.0,
                               "E_surf corrected fixed-flux (J)": 4.0}]}
        observed = {"x.csv": [{"m": 1.0, "p_surf corrected": 2.0}]}
        self.assertEqual(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20, parked_columns=False), []
        )
        self.assertTrue(BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20))
        observed["x.csv"][0]["p_surf corrected"] = 2.1
        self.assertTrue(
            BENCHMARK.compare_snapshots(expected, observed, 1.0e-3, 1.0e-20, parked_columns=False)
        )


if __name__ == "__main__":
    unittest.main()
