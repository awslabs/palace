#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of corner_validation/prepare_maxwell_validation.py's two forms (the v1 cut, developer_form.py)."""

import json
import unittest
from pathlib import Path

import prepare_maxwell_validation as PREPARE

developer_form = PREPARE.developer_form
# No run record of this writer is on the development host: the record is the config this writer wrote at e4195bac6c
# (before the cut, crack policy legacy) for the arguments below (testdata/, sha256 a30a36543a02e1a2...).
RECORD = Path(__file__).resolve().parent / "testdata" / "record-main-e4195bac6c-maxwell-thin-corrected.json"
ARGUMENTS = (Path("/records/postpro/thin-corrected"), Path("/records/island-thin.msh"), 2, 5.0, False,
             Path("/records/process-library.json"), 1.9)


class PrepareCornerMaxwellValidationForms(unittest.TestCase):
    def test_default_corrected_config_carries_no_parked_key(self):
        config = PREPARE.config(*ARGUMENTS)
        self.assertEqual(developer_form.parked_config_settings(config), [])
        self.assertNotIn("Developer", config["Solver"]["SurfaceResponseCorrection"])
        self.assertNotIn("CrackInternalBoundaryElements", config["Model"])

    def test_developer_form_is_the_record(self):
        # The recorded developer form = the record byte for byte with the default of record "CorrectionMode": "Both"
        # stated and "Developer": true.
        text = RECORD.read_text()
        record = json.loads(text)
        self.assertEqual(json.dumps(record, indent=2) + "\n", text)
        config = PREPARE.config(*ARGUMENTS, crack_policy=developer_form.device_run_crack_policy(None, True), developer=True)
        record["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both"
        record["Solver"]["SurfaceResponseCorrection"]["Developer"] = True
        self.assertEqual(json.dumps(config, indent=2) + "\n", json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
