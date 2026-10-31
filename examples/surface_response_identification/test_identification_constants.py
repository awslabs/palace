#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""C++ == Python: every matching-band constant of identification_constants.py equals its constexpr in
palace/models/surfaceresponseidentification.hpp (identification fix round DESIGN section 4 (f) / A.10 G-6)."""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import identification_constants as ic  # noqa: E402


class TestIdentificationConstants(unittest.TestCase):
    def test_header_exists(self):
        self.assertTrue(os.path.isfile(ic.CXX_HEADER), ic.CXX_HEADER)

    def test_every_constant_equals_the_cxx(self):
        for py_name, cxx_name in ic.CXX_NAMES.items():
            with self.subTest(py_name):
                self.assertEqual(getattr(ic, py_name), ic.cxx_constant(cxx_name), f"{py_name} != {cxx_name}")

    def test_cxx_constant_fails_closed_on_an_unknown_name(self):
        with self.assertRaises(ValueError):
            ic.cxx_constant("kNoSuchConstantOfTheIdentification")

    def test_pair_band_is_the_max(self):
        self.assertEqual(ic.pair_near_match_band(1.0e-3, 1.9), 5.0e-3 * 1.9)
        self.assertEqual(ic.pair_near_match_band(0.5, 1.9), 0.5)


if __name__ == "__main__":
    unittest.main()
