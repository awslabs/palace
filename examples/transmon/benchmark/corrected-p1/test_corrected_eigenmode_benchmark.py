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


if __name__ == "__main__":
    unittest.main()
