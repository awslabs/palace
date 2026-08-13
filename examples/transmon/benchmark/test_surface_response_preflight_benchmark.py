#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("run_surface_response_preflight.py")
SPEC = importlib.util.spec_from_file_location("surface_response_benchmark", MODULE_PATH)
BENCHMARK = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCHMARK)


class SurfaceResponsePreflightBenchmarkTest(unittest.TestCase):
    def setUp(self):
        self.library = Path(__file__).with_name("transmon_surface_process_seed.json")
        self.manifest = {
            "Version": 1,
            "Complete": False,
            "Library": {
                "Path": str(self.library.resolve()),
                "Name": "process",
                "MatchingRadius": 2.0,
            },
            "LengthUnit": "mesh",
            "Summary": {
                "Counts": {"Exact": 0, "Interpolated": 0, "Missing": 3},
                "TotalEdgeLengths": {
                    "Exact": 0.0,
                    "Interpolated": 0.0,
                    "Missing": 6.0,
                },
            },
            "Statistics": {"PairChecks": 100},
            "Requirements": [
                {
                    "Topology": "IsolatedEdge",
                    "Status": "Missing",
                    "Count": 2,
                    "TotalEdgeLength": 4.0,
                },
                {
                    "Topology": "ConvexCorner",
                    "Status": "Missing",
                    "Count": 1,
                    "TotalEdgeLength": 2.0,
                },
            ],
        }

    def test_canonical_manifest_ignores_statistics(self):
        first = BENCHMARK.canonical_sha256(self.manifest, self.library)
        changed = copy.deepcopy(self.manifest)
        changed["Statistics"] = {"PairChecks": 1}
        self.assertEqual(
            first, BENCHMARK.canonical_sha256(changed, self.library)
        )

    def test_canonical_manifest_detects_requirement_changes(self):
        changed = copy.deepcopy(self.manifest)
        changed["Requirements"][0]["Count"] += 1
        self.assertNotEqual(
            BENCHMARK.canonical_sha256(self.manifest, self.library),
            BENCHMARK.canonical_sha256(changed, self.library),
        )

    def test_canonical_manifest_rejects_wrong_or_relative_library(self):
        changed = copy.deepcopy(self.manifest)
        changed["Library"]["Path"] = "relative.json"
        with self.assertRaisesRegex(ValueError, "must be absolute"):
            BENCHMARK.canonical_manifest(changed, self.library)
        changed["Library"]["Path"] = str(MODULE_PATH.resolve())
        with self.assertRaisesRegex(ValueError, "expected"):
            BENCHMARK.canonical_manifest(changed, self.library)

    def test_explicit_mpi_command(self):
        command = BENCHMARK.build_command(
            Path("/usr/bin/mpiexec"),
            ["--bind-to", "core"],
            "-n",
            4,
            Path("/repo/build/bin/palace"),
            Path("/tmp/config.json"),
        )
        self.assertEqual(
            command,
            [
                "/usr/bin/mpiexec",
                "--bind-to",
                "core",
                "-n",
                "4",
                "/repo/build/bin/palace",
                "--serial",
                "--surface-response-preflight",
                "/tmp/config.json",
            ],
        )

    def test_manifest_summary_groups_topologies(self):
        summary = BENCHMARK.manifest_summary(self.manifest)
        self.assertEqual(summary["RequirementRecords"], 2)
        self.assertEqual(summary["ByTopology"]["IsolatedEdge"]["Count"], 2)
        self.assertEqual(
            summary["ByTopology"]["ConvexCorner"]["TotalEdgeLength"], 2.0
        )

    def test_prepare_config_is_bounded(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "source.json"
            mesh = tmp / "mesh.msh"
            mesh.write_text("fixture")
            source.write_text(
                json.dumps(
                    {
                        "Problem": {"Output": "unused", "Type": "Eigenmode"},
                        "Model": {"Mesh": str(mesh), "Refinement": {"MaxIts": 4}},
                        "Boundaries": {
                            "Postprocessing": {
                                "Dielectric": [
                                    {
                                        "Index": 1,
                                        "EdgeDistances": [1.0, 2.0],
                                        "LocalizeEdgeEnergy": True,
                                    }
                                ]
                            }
                        },
                        "Solver": {"Eigenmode": {"Save": 2}},
                    }
                )
            )
            output = tmp / "config.json"
            config, observed_mesh = BENCHMARK.prepare_config(
                source, self.library, output, tmp / "postpro"
            )
            self.assertEqual(observed_mesh, mesh.resolve())
            self.assertEqual(config["Model"]["Refinement"]["MaxIts"], 0)
            self.assertEqual(config["Model"]["Refinement"]["UniformLevels"], 0)
            self.assertFalse(config["Problem"]["OutputFormats"]["Paraview"])
            dielectric = config["Boundaries"]["Postprocessing"]["Dielectric"][0]
            self.assertEqual(dielectric["EdgeDistances"], [2.0])
            self.assertFalse(dielectric["LocalizeEdgeEnergy"])
            self.assertFalse(dielectric["SaveLocalEdgeEnergy"])

    def test_prepare_config_rejects_edge_refinement(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            mesh = tmp / "mesh.msh"
            mesh.write_text("fixture")
            source = tmp / "source.json"
            source.write_text(
                json.dumps(
                    {
                        "Problem": {"Output": "unused", "Type": "Eigenmode"},
                        "Model": {"Mesh": str(mesh), "Refinement": {}},
                        "Boundaries": {
                            "Postprocessing": {
                                "Dielectric": [
                                    {
                                        "Index": 1,
                                        "EdgeRefinement": {
                                            "Radius": 2.0,
                                            "ElementsPerRadius": 1,
                                        },
                                    }
                                ]
                            }
                        },
                        "Solver": {"Eigenmode": {"Save": 0}},
                    }
                )
            )
            with self.assertRaisesRegex(ValueError, "forbids EdgeRefinement"):
                BENCHMARK.prepare_config(
                    source, self.library, tmp / "config.json", tmp / "postpro"
                )

    def test_checked_in_baseline_is_complete(self):
        baseline = json.loads(
            Path(__file__)
            .with_name("surface-response-preflight-baseline.json")
            .read_text()
        )
        self.assertEqual(baseline["SchemaVersion"], 1)
        self.assertEqual(baseline["Benchmark"], BENCHMARK.BENCHMARK_NAME)
        self.assertEqual(len(baseline["CanonicalManifestSha256"]), 64)
        self.assertEqual(baseline["Summary"]["Counts"]["Missing"], 215)
        self.assertEqual(len(baseline["CanonicalManifest"]["Requirements"]), 22)


if __name__ == "__main__":
    unittest.main()
