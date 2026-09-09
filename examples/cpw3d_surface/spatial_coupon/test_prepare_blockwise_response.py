#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import json
import subprocess
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("prepare_blockwise_response.py")


class PrepareBlockwiseResponseTest(unittest.TestCase):
    def test_split_preserves_source_indices_and_writes_reducer(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.json"
            config.write_text(
                json.dumps(
                    {
                        "Problem": {"Type": "Electrostatic", "Output": "old"},
                        "Boundaries": {
                            "PrescribedPotential": [
                                {
                                    "Index": index,
                                    "Attributes": [1],
                                    "DataFile": f"{index}.csv",
                                }
                                for index in (3, 7, 11, 19, 23)
                            ]
                        },
                        "Solver": {
                            "Electrostatic": {
                                "ResponseMatrix": True,
                                "AggregateResponseMatrix": True,
                            }
                        },
                    }
                )
                + "\n"
            )
            output = root / "generation"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    str(config),
                    "--output",
                    str(output),
                    "--execution-mode",
                    "blocks",
                    "--sources-per-block",
                    "2",
                    "--reduction-block-size",
                    "3",
                    "--ranks",
                    "8",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            manifest = json.loads((output / "blockwise-manifest.json").read_text())
            self.assertEqual(manifest["SourceCount"], 5)
            self.assertEqual(manifest["BlockCount"], 3)
            self.assertEqual(
                [entry["Sources"] for entry in manifest["Blocks"]],
                [[3, 7], [11, 19], [23]],
            )
            block = json.loads((output / "blocks/block-0001.json").read_text())
            self.assertEqual(
                [entry["Index"] for entry in block["Boundaries"]["PrescribedPotential"]],
                [11, 19],
            )
            self.assertTrue(block["Solver"]["Electrostatic"]["ResponseMatrix"])
            reducer = json.loads((output / "reducer.json").read_text())
            self.assertEqual(
                [entry["Index"] for entry in reducer["Boundaries"]["PrescribedPotential"]],
                [3, 7, 11, 19, 23],
            )
            self.assertTrue(reducer["Solver"]["Electrostatic"]["ResponseMatrix"])
            reducer_script = (output / "run-reducer.sh").read_text()
            self.assertIn("PALACE_RESPONSE_REDUCE_ONLY=1", reducer_script)
            self.assertIn("PALACE_RESPONSE_BLOCK_SIZE=3", reducer_script)

    def test_streaming_reuses_one_process_and_releases_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.json"
            config.write_text(
                json.dumps(
                    {
                        "Problem": {"Type": "Electrostatic"},
                        "Boundaries": {
                            "PrescribedPotential": [
                                {"Index": index, "Attributes": [1]}
                                for index in (2, 5, 9, 14)
                            ]
                        },
                        "Solver": {"Electrostatic": {}},
                    }
                )
                + "\n"
            )
            output = root / "streaming"
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    str(config),
                    "--output",
                    str(output),
                    "--ranks",
                    "6",
                    "--local",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            manifest = json.loads((output / "blockwise-manifest.json").read_text())
            self.assertEqual(manifest["ExecutionMode"], "streaming")
            self.assertTrue(manifest["Local"])
            self.assertEqual(manifest["BlockCount"], 1)
            self.assertEqual(manifest["Blocks"][0]["Sources"], [2, 5, 9, 14])
            script = (output / "run-blocks.sh").read_text()
            self.assertIn("PALACE_RESPONSE_ARCHIVE_ONLY=1", script)
            self.assertNotIn("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS=1", script)
            self.assertNotIn("PBS_NODEFILE", script)
            self.assertIn('"$MPIEXEC" -n 6', script)

    def test_rejects_duplicate_source_indices(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "config.json"
            config.write_text(
                json.dumps(
                    {
                        "Boundaries": {
                            "PrescribedPotential": [
                                {"Index": 1},
                                {"Index": 1},
                            ]
                        }
                    }
                )
                + "\n"
            )
            result = subprocess.run(
                ["python3", str(SCRIPT), str(config), "--output", str(root / "out")],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unique and positive", result.stderr)


if __name__ == "__main__":
    unittest.main()
