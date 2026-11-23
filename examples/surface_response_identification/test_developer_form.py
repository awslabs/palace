#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of developer_form (the v1 cut): the v1 form drops the parked options; the recorded developer form states
the defaults of record and appends the switch; the parked options mirror Palace's configfile.cpp."""

import os
import re
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import developer_form as F  # noqa: E402

CONFIGFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "palace", "utils", "configfile.cpp")


class DeveloperFormTest(unittest.TestCase):
    def test_v1_block_drops_the_self_consistent_solve(self):
        recorded = {"Library": "l", "CorrectionMode": "Both", "SolveTol": 1e-6, "TranslationalDomainCorrection": "FixedTrace",
                    "PatchConstruction": "Features"}
        self.assertEqual(F.v1_block(recorded), {"Library": "l", "CorrectionMode": "PostprocessOnly", "PatchConstruction": "Features"})
        self.assertEqual(F.v1_block({"Library": "l", "SolveTol": 1e-6}), {"Library": "l"})
        self.assertEqual(F.parked_block_settings(F.v1_block(recorded)), [])
        with self.assertRaises(ValueError):
            F.v1_block({"Library": "l", "Developer": True})

    def test_v1_block_drops_the_mortar_coupling(self):
        recorded = {"Library": "l", "TraceCoupling": "SurfaceMortar", "MortarOversampling": 2}
        self.assertEqual(F.v1_block(recorded), {"Library": "l", "TraceCoupling": "Collocated"})
        self.assertEqual(F.parked_block_settings(recorded), ["TraceCoupling SurfaceMortar", "MortarOversampling"])
        self.assertEqual(F.parked_block_settings({"Models": []}), ["the hand-placed form (no Library)"])

    def test_developer_block_states_the_defaults_of_record_and_the_switch_last(self):
        self.assertEqual(list(F.developer_block({"Library": "l", "SolveTol": 1e-6}).items()),
                         [("Library", "l"), ("SolveTol", 1e-6), ("CorrectionMode", "Both"), ("Developer", True)])
        self.assertEqual(list(F.developer_block({"Library": "l", "CorrectionMode": "SelfConsistent"}).items()),
                         [("Library", "l"), ("CorrectionMode", "SelfConsistent"), ("Developer", True)])
        self.assertEqual(list(F.with_switch({"Developer": False, "Library": "l"}).items()), [("Library", "l"), ("Developer", True)])

    def test_parked_settings_of_a_block_and_a_config(self):
        block = {"Library": "l", "CorrectionMode": "Both", "SolveTol": 1e-6}
        self.assertEqual(F.parked_block_settings(block), ["CorrectionMode Both", "SolveTol"])
        config = {"Solver": {"Electrostatic": {"ResponseCorrection": dict(block, Developer=True)}}}
        self.assertEqual(F.parked_config_settings(config), ["Solver.Electrostatic.ResponseCorrection: CorrectionMode Both",
                                                            "Solver.Electrostatic.ResponseCorrection: SolveTol",
                                                            "Solver.Electrostatic.ResponseCorrection: Developer"])
        self.assertEqual(F.parked_config_settings({"Solver": {"SurfaceResponseCorrection": {"Library": "l"}}}), [])
        self.assertEqual(F.form_name(False), "v1 form")
        self.assertEqual(F.form_name(True), "recorded developer form")

    def test_device_run_crack_policy_per_form(self):
        self.assertEqual(F.device_run_crack_policy(None, False), "default")
        self.assertEqual(F.device_run_crack_policy(None, True), "legacy")
        self.assertEqual(F.device_run_crack_policy("default", True), "default")
        with self.assertRaises(ValueError):
            F.device_run_crack_policy("legacy", False)
        cracked = {"Model": {"CrackInternalBoundaryElements": True}, "Solver": {"SurfaceResponseCorrection": {"Library": "l"}}}
        self.assertEqual(F.parked_config_settings(cracked), ["Model.CrackInternalBoundaryElements true"])
        del cracked["Solver"]  # no correction run: the key is not parked
        self.assertEqual(F.parked_config_settings(cracked), [])

    def test_parked_options_mirror_palace(self):
        # Un-parking an option deletes its PARKED_OPTIONS entry and its RequireDeveloperOption line in Palace together.
        with open(CONFIGFILE) as source:
            text = source.read()
        guards = set(re.findall(r'RequireDeveloperOption\(\s*data,\s*"((?:[^"\\]|\\.)*)"\)', text))
        for setting, values, keys, _ in F.PARKED_OPTIONS.values():
            for value in values:
                self.assertIn(f'\\"{setting}\\": \\"{value}\\"', guards)
            for key in keys:  # read with the option only, else named among the ignored keys
                self.assertIn(f'"{key}"', text)
        self.assertIn("ignored_keys.push_back", text)


if __name__ == "__main__":
    unittest.main()
