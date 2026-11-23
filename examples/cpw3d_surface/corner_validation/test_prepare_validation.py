# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of prepare_validation.py (the corner validation / corner-family device-check configs)."""

import json
import tempfile
import unittest
from pathlib import Path

import prepare_validation as PREPARE


class PrepareValidationTest(unittest.TestCase):
    def test_fabricated_edge_lines_are_the_sa_perimeter(self):
        # The fabricated device is a 3D metal slab whose every edge is a fold between metal
        # faces (MS bottom, MA sidewalls, MA top): the automatic perimeter extraction retains no
        # one-sided metal edge on it ("No physical metal perimeter was found"), so the fabricated
        # reference names its edge lines as the SA perimeter minus the outer box, as the
        # fabricated corner coupon does; the thin corrected run keeps the automatic perimeter.
        for fabricated in (False, True):
            config = PREPARE.config(
                Path("postpro"), Path("mesh.msh"), 3, 3, fabricated, Path("lib.json"), 1.9, 1
            )
            dielectrics = config["Boundaries"]["Postprocessing"]["Dielectric"]
            self.assertEqual([entry["Type"] for entry in dielectrics], ["SA", "MS", "MA"])
            for entry in dielectrics:
                self.assertEqual(entry["EdgeDistances"], [1.9])
                self.assertEqual(entry["EdgeExcludeAttributes"], PREPARE.OUTER_ATTRIBUTES)
                if fabricated:
                    self.assertNotIn("AutomaticEdges", entry)
                    self.assertEqual(entry["EdgeAttributes"], [PREPARE.SA_ATTRIBUTE])
                    self.assertNotIn("EdgeRefinement", entry)
                else:
                    self.assertTrue(entry["AutomaticEdges"])
                    self.assertNotIn("EdgeAttributes", entry)
                    self.assertEqual(entry["EdgeRefinement"]["ElementsPerRadius"], 1)
            self.assertEqual(dielectrics[0]["Attributes"], [PREPARE.SA_ATTRIBUTE])
            self.assertEqual(
                config["Boundaries"]["Terminal"][0]["Attributes"], [2, 4] if fabricated else [2]
            )
            self.assertEqual(
                "ResponseCorrection" in config["Solver"]["Electrostatic"], not fabricated
            )

    def test_process_from_library_reads_the_fabrication_record(self):
        with tempfile.TemporaryDirectory() as directory:
            library = Path(directory) / "process-library.json"
            library.write_text(
                json.dumps(
                    {
                        "Fabrication": {
                            "SubstratePermittivity": 11.45,
                            "InterfaceLayers": {"MS": {"Permittivity": 11.45, "Thickness": 0.002}},
                        }
                    }
                )
            )
            interfaces, substrate = PREPARE.process_from_library(library)
        self.assertEqual(substrate, 11.45)
        self.assertEqual(interfaces["MS"], (11.45, 3.0e-4, 0.002))
        self.assertEqual(interfaces["SA"], PREPARE.INTERFACES["SA"])


# The record: corner-basis-refinement-20260930/device-check-convex/square-s10-p3/thin-corrected.json (sha256
# a8d967a8fc6e89b3..., 2026-09-30; copied to testdata/), written with the library's Fabrication layers below.
RECORD = Path(__file__).resolve().parent / "testdata" / "record-square-s10-p3-thin-corrected.json"
RECORD_PROCESS = ({"SA": (4.0, 0.002, 0.002), "MS": (11.45, 0.0003, 0.002), "MA": (10.0, 0.03, 0.002)}, 11.45)


def record_arguments(record):
    block = record["Solver"]["Electrostatic"]["ResponseCorrection"]
    dielectric = record["Boundaries"]["Postprocessing"]["Dielectric"][0]
    return (Path(record["Problem"]["Output"]), Path(record["Model"]["Mesh"]), record["Solver"]["Order"],
            record["Model"]["Refinement"]["MaxIts"], False, Path(block["Library"]), dielectric["EdgeDistances"][0],
            dielectric["EdgeRefinement"]["ElementsPerRadius"], RECORD_PROCESS)


class ResponseFormsTest(unittest.TestCase):
    """The v1 cut (developer_form.py): the default thin-corrected config is the v1 form; --developer the record."""

    def test_default_thin_corrected_config_carries_no_parked_key(self):
        config = PREPARE.config(*record_arguments(json.loads(RECORD.read_text())))
        self.assertEqual(PREPARE.developer_form.parked_config_settings(config), [])
        self.assertNotIn("Developer", config["Solver"]["Electrostatic"]["ResponseCorrection"])
        self.assertNotIn("CrackInternalBoundaryElements", config["Model"])
        # The protocol's mortar trace lift is parked: the v1 form collocates.
        self.assertEqual(config["Solver"]["Electrostatic"]["ResponseCorrection"],
                         {"Library": "/data/home/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/corner-basis-refinement-20260930/"
                                     "libraries/corner-basis-refinement-R1p9-convex/process-library.json",
                          "TargetInterfaces": [1, 2, 3], "UnmatchedPolicy": "Error", "TraceCoupling": "Collocated"})

    def test_developer_form_is_the_record(self):
        # The recorded developer form = the record byte for byte plus the keys that state its behaviour on today's
        # binary: Model.CrackInternalBoundaryElements true (the record's binary cracked an unset key; crack policy
        # legacy), the default of record "CorrectionMode": "Both" and "Developer": true.
        text = RECORD.read_text()
        record = json.loads(text)
        self.assertEqual(json.dumps(record, indent=2) + "\n", text)
        config = PREPARE.config(*record_arguments(record), crack_policy="legacy", developer=True)
        record["Model"]["CrackInternalBoundaryElements"] = True
        record["Solver"]["Electrostatic"]["ResponseCorrection"]["CorrectionMode"] = "Both"
        record["Solver"]["Electrostatic"]["ResponseCorrection"]["Developer"] = True
        self.assertEqual(json.dumps(config, indent=2) + "\n", json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main()
