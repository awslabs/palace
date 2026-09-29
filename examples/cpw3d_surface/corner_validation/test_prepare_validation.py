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


if __name__ == "__main__":
    unittest.main()
