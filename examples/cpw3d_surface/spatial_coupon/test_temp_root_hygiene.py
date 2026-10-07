# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The near-key and spatial-qualification suites leave the temporary root EMPTY: every suite runs in a
subprocess whose TMPDIR is a fresh empty directory, and that directory holds nothing afterwards.
(test_nearkey_predictor's module-level mkdtemp leaked one nearkey-activation-* directory per run,
qualify/test_spatial_qualification's two setUp mkdtemp fixtures one spatial-qualification-* directory per
test: 55 + 235 directories under the host's temp root when this test was written.)"""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent

# (working directory, unittest module names) of every suite checked; each runs once per invocation.
SUITES = (
    (HERE, ("test_nearkey_predictor", "test_nearkey_window_report", "test_nearkey_reuse")),
    (HERE / "qualify", ("test_spatial_qualification",)),
)


class TempRootHygieneTest(unittest.TestCase):
    def test_near_key_and_qualification_suites_leave_the_temp_root_empty(self):
        for directory, modules in SUITES:
            with self.subTest(modules=modules):
                with tempfile.TemporaryDirectory(prefix="temp-root-hygiene-") as temp_root:
                    environment = {**os.environ, "TMPDIR": temp_root, "PYTHONDONTWRITEBYTECODE": "1"}
                    result = subprocess.run([sys.executable, "-m", "unittest", *modules], cwd=directory, env=environment,
                                            text=True, capture_output=True)
                    self.assertEqual(result.returncode, 0, result.stderr[-4000:])
                    self.assertIn("OK", result.stderr.splitlines()[-1])
                    left = sorted(os.listdir(temp_root))
                    self.assertEqual(left, [], f"{modules} left {left} under the temp root")


if __name__ == "__main__":
    unittest.main()
