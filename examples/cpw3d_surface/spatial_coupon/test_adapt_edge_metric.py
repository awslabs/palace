# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Fail-closed publication tests for the optional native MMG adapter."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


EXECUTABLE = os.environ.get("EDGE_METRIC_ADAPTER_EXE")


@unittest.skipUnless(EXECUTABLE and Path(EXECUTABLE).is_file(),
                     "EDGE_METRIC_ADAPTER_EXE is not available")
class EdgeMetricAdapterPublicationTest(unittest.TestCase):
    def invoke(self, root, output):
        return subprocess.run(
            [EXECUTABLE, str(root / "missing.meshb"), str(root / "missing.f64"),
             str(root / "missing.txt"), str(output), ".01", ".1", "1.2"],
            capture_output=True, text=True, timeout=10)

    def test_existing_solution_file_is_not_mutated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "accepted.meshb"
            solution = Path(str(output) + ".sol")
            solution.write_bytes(b"immutable solution")
            result = self.invoke(root, output)
            self.assertEqual(result.returncode, 2)
            self.assertIn("Refusing to overwrite", result.stderr)
            self.assertEqual(solution.read_bytes(), b"immutable solution")
            self.assertFalse(output.exists())

    def test_unsavable_existing_solution_directory_is_not_mutated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "accepted.meshb"
            solution = Path(str(output) + ".sol")
            solution.mkdir()
            result = self.invoke(root, output)
            self.assertEqual(result.returncode, 2)
            self.assertIn("Refusing to overwrite", result.stderr)
            self.assertTrue(solution.is_dir())
            self.assertFalse(output.exists())

    def test_existing_mesh_and_rejected_attempt_are_not_mutated(self):
        for suffix in ("", ".rejected.meshb"):
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                output = root / "accepted.meshb"
                protected = Path(str(output) + suffix)
                protected.write_bytes(b"immutable attempt")
                result = self.invoke(root, output)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(protected.read_bytes(), b"immutable attempt")


if __name__ == "__main__":
    unittest.main()
