#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import importlib.util
import math
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("replay_device_trace.py")
SPEC = importlib.util.spec_from_file_location("replay_device_trace", MODULE_PATH)
REPLAY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPLAY)


class ReplayDeviceTraceTest(unittest.TestCase):
    def test_relative_error_is_symmetric_and_scaled(self):
        self.assertAlmostEqual(REPLAY.relative_error(2.0, 2.2, 10.0), 0.2 / 2.2)
        self.assertAlmostEqual(REPLAY.relative_error(2.2, 2.0, 10.0), 0.2 / 2.2)
        self.assertEqual(REPLAY.relative_error(0.0, 0.0, 1.0), 0.0)
        self.assertAlmostEqual(
            REPLAY.relative_error(0.0, 1.0e-15, 1.0),
            0.1,
        )

    def test_relative_error_rejects_nonfinite_values(self):
        self.assertEqual(REPLAY.relative_error(math.nan, 1.0, 1.0), math.inf)
        self.assertEqual(REPLAY.relative_error(1.0, math.inf, 1.0), math.inf)
        self.assertEqual(REPLAY.relative_error(1.0, 1.0, math.nan), math.inf)


if __name__ == "__main__":
    unittest.main()
