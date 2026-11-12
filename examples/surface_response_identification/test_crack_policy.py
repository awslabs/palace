#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of crack_policy (decision 1570): one shared default; legacy writes true; default omits."""

import argparse
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import crack_policy as C  # noqa: E402


class CrackPolicyTest(unittest.TestCase):
    def test_shared_default_is_legacy_in_this_integration(self):
        self.assertEqual(C.CRACK_POLICY_DEFAULT, "legacy")
        self.assertEqual(C.CRACK_POLICIES, ("legacy", "default"))

    def test_legacy_writes_true_and_default_omits(self):
        self.assertEqual(C.apply_crack_policy({"Mesh": "m"}, "legacy"), {"Mesh": "m", "CrackInternalBoundaryElements": True})
        self.assertEqual(C.apply_crack_policy({"Mesh": "m"}, "default"), {"Mesh": "m"})

    def test_refuses_unknown_policy_and_a_preset_key(self):
        with self.assertRaises(ValueError):
            C.apply_crack_policy({}, "true")
        with self.assertRaises(ValueError):
            C.apply_crack_policy({"CrackInternalBoundaryElements": False}, "legacy")

    def test_argument_defaults_to_the_shared_constant(self):
        parser = argparse.ArgumentParser()
        C.add_crack_policy_argument(parser)
        self.assertEqual(parser.parse_args([]).crack_policy, C.CRACK_POLICY_DEFAULT)
        self.assertEqual(parser.parse_args(["--crack-policy", "default"]).crack_policy, "default")
        with self.assertRaises(SystemExit):
            parser.parse_args(["--crack-policy", "none"])


if __name__ == "__main__":
    unittest.main()
