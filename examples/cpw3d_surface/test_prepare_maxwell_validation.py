#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of prepare_maxwell_validation.py's two forms (the v1 cut, developer_form.py)."""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "prepare_maxwell_validation.py"
sys.path.insert(0, str(HERE.parent / "surface_response_identification"))
import developer_form  # noqa: E402

# The record: v1/lanes/scope/driven/prep/thin_corrected.json (sha256 673c69b2d8910e2a..., written by this script at
# e4195bac6c with its then default crack policy legacy; copied to testdata/).
RECORD = HERE / "testdata" / "record-scope-driven-prep-thin_corrected.json"


def prepare(directory, *extra):
    directory = Path(directory)
    meshes = directory / "mesh"
    meshes.mkdir()
    for name in ("cpw3d_surface_maxwell_thin.msh", "cpw3d_surface_maxwell_fabricated.msh"):
        (meshes / name).write_text("")
    library = directory / "process-library.json"
    library.write_text("{}\n")
    output = directory / "prep"
    command = [sys.executable, str(SCRIPT), "--output", str(output), "--library", str(library),
               "--mesh-directory", str(meshes), *extra]
    result = subprocess.run(command, capture_output=True, text=True)
    return result, output, meshes, library


class PrepareMaxwellValidationForms(unittest.TestCase):
    def test_default_corrected_config_carries_no_parked_key(self):
        with tempfile.TemporaryDirectory() as directory:
            result, output, _, _ = prepare(directory)
            self.assertEqual(result.returncode, 0, result.stderr)
            corrected = json.loads((output / "thin_corrected.json").read_text())
        self.assertEqual(developer_form.parked_config_settings(corrected), [])
        self.assertNotIn("Developer", corrected["Solver"]["SurfaceResponseCorrection"])
        self.assertIn("the v1 form", result.stderr)
        with tempfile.TemporaryDirectory() as directory:
            result, _, _, _ = prepare(directory, "--crack-policy", "legacy")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--developer", result.stderr)

    def test_developer_form_is_the_record(self):
        # The recorded developer form = the record byte for byte (its paths mapped to this test's) with the default of
        # record "CorrectionMode": "Both" stated and "Developer": true.
        text = RECORD.read_text()
        record = json.loads(text)
        self.assertEqual(json.dumps(record, indent=2) + "\n", text)
        with tempfile.TemporaryDirectory() as directory:
            result, output, meshes, library = prepare(directory, "--developer")
            self.assertEqual(result.returncode, 0, result.stderr)
            corrected = (output / "thin_corrected.json").read_text()
        record["Problem"]["Output"] = str(output / "thin_corrected")
        record["Model"]["Mesh"] = str(meshes / "cpw3d_surface_maxwell_thin.msh")
        record["Solver"]["SurfaceResponseCorrection"]["Library"] = str(library)
        record["Solver"]["SurfaceResponseCorrection"]["CorrectionMode"] = "Both"
        record["Solver"]["SurfaceResponseCorrection"]["Developer"] = True
        self.assertEqual(corrected, json.dumps(record, indent=2) + "\n")
        self.assertIn("the recorded developer form", result.stderr)


if __name__ == "__main__":
    unittest.main()
