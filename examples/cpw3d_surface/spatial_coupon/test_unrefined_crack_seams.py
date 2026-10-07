# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Mesher design round 2 SEAM, supervisor decision 368 (3): the solve-config writer sets
Model.RefineCrackElements false exactly when the mesh's build census records
ThinSheetSeams.UnrefinedCrackSeams (the pinched seams of the thin tips sharper than the tip
bisector's minimum opening), fail closed on a fabricated case, a non-Dirichlet sheet or an
inconsistent census."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import case_inputs  # noqa: E402

CONFIG = {"Model": {"Mesh": "canonical.msh", "L0": 1e-6, "Refinement": {"MaxIts": 0}},
          "Boundaries": {"Ground": {"Attributes": [1]},
                         "PrescribedPotential": [{"Index": 1, "Attributes": [4001], "DataFile": "t.csv"}]}}
CENSUS = {"InterfaceAreas": [{"Attribute": 1, "Area": 1.0}, {"Attribute": 3000, "Area": 1.0},
                             {"Attribute": 4001, "Area": 0.5}],
          "ThinSheetSeams": {"Count": 6, "Edges": [[[0, 0, 0], [1, 0, 0]]] * 6,
                             "UnrefinedCrackSeams": {"Rule": "... Model.RefineCrackElements false ...",
                                                     "Count": 6, "RefineCrackElements": False,
                                                     "Tips": [{"Point": [0, 0, 0], "OpeningDegrees": 22.5}]}}}


class UnrefinedCrackSeamsTest(unittest.TestCase):
    def write(self, census):
        directory = Path(tempfile.mkdtemp())
        mesh = directory / "canonical.msh"
        mesh.write_text("mesh")
        if census is not None:
            (directory / case_inputs.BUILD_CENSUS_NAME).write_text(json.dumps(census))
        return mesh

    def test_thin_case_with_recorded_seams_runs_with_refine_crack_elements_false(self):
        config = copy.deepcopy(CONFIG)
        record = case_inputs.apply_unrefined_crack_seams(config, self.write(CENSUS), False)
        self.assertIs(config["Model"]["RefineCrackElements"], False)
        self.assertEqual(record["Count"], 6)
        self.assertEqual(record["SheetAttributes"], [4001])
        self.assertIn("Dirichlet", record["Rule"])
        self.assertIn("decision 368", record["Rule"])

    def test_no_census_or_no_seams_leaves_the_default(self):
        config = copy.deepcopy(CONFIG)
        # A fabricated case without a census beside its mesh: nothing to apply (a fabricated
        # build never carries the record).
        self.assertIsNone(case_inputs.apply_unrefined_crack_seams(config, self.write(None), True))
        clean = copy.deepcopy(CENSUS)
        clean["ThinSheetSeams"] = {"Count": 0, "Edges": [], "UnrefinedCrackSeams": None}
        self.assertIsNone(case_inputs.apply_unrefined_crack_seams(config, self.write(clean), False))
        # A fabricated build records no thin sheet (ThinSheetSeams null).
        fabricated = copy.deepcopy(CENSUS)
        fabricated["ThinSheetSeams"] = None
        self.assertIsNone(case_inputs.apply_unrefined_crack_seams(config, self.write(fabricated), True))
        self.assertNotIn("RefineCrackElements", config["Model"])

    def test_thin_case_without_a_census_fails_closed(self):
        # Decision 392 MINOR-5 (decision 368 (3) "fail closed"): a thin case whose mesh has no
        # build-census.json beside it is refused - its UnrefinedCrackSeams record cannot be
        # read, so Model.RefineCrackElements is not silently left at the default.
        config = copy.deepcopy(CONFIG)
        with self.assertRaisesRegex(case_inputs.CaseInputError, "thin case has no build census"):
            case_inputs.apply_unrefined_crack_seams(config, self.write(None), False)
        self.assertNotIn("RefineCrackElements", config["Model"])
        self.assertIn("without a build census", case_inputs.UNREFINED_CRACK_SEAMS_RULE)

    def test_fail_closed(self):
        with self.assertRaisesRegex(case_inputs.CaseInputError, "fabricated case carries UnrefinedCrackSeams"):
            case_inputs.apply_unrefined_crack_seams(copy.deepcopy(CONFIG), self.write(CENSUS), True)
        seams_without_record = copy.deepcopy(CENSUS)
        seams_without_record["ThinSheetSeams"]["UnrefinedCrackSeams"] = None
        with self.assertRaisesRegex(case_inputs.CaseInputError, "without an UnrefinedCrackSeams record"):
            case_inputs.apply_unrefined_crack_seams(copy.deepcopy(CONFIG), self.write(seams_without_record), False)
        inconsistent = copy.deepcopy(CENSUS)
        inconsistent["ThinSheetSeams"]["UnrefinedCrackSeams"]["Count"] = 5
        with self.assertRaisesRegex(case_inputs.CaseInputError, "disagrees with ThinSheetSeams"):
            case_inputs.apply_unrefined_crack_seams(copy.deepcopy(CONFIG), self.write(inconsistent), False)
        # The sheet must be a Dirichlet boundary of the config (a PEC sheet, F.7).
        non_pec = copy.deepcopy(CONFIG)
        non_pec["Boundaries"]["PrescribedPotential"][0]["Attributes"] = [4002]
        with self.assertRaisesRegex(case_inputs.CaseInputError, "not all Dirichlet boundaries"):
            case_inputs.apply_unrefined_crack_seams(non_pec, self.write(CENSUS), False)


if __name__ == "__main__":
    unittest.main()
