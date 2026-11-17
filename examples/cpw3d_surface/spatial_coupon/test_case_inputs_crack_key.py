# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""case_inputs.config_differences and Model.CrackInternalBoundaryElements (decision 1582 (2); SL9 review
MAJOR-1 of the upstream integration): a config written now by the legacy crack policy (the key true) equals
a STORED reference made before the policy (the key absent = true on the pre-integration binary), and only
that pair is equal. The stored references are read, never written (their bytes are checked)."""
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "qualify"))
import case_inputs  # noqa: E402
from crack_policy import CRACK_KEY, apply_crack_policy  # noqa: E402

# In-tree stored producer configs (registered gallery inputs 06 / 10 / 05, decision 1405 era).
STORED = sorted((HERE / "testdata").glob("*/spatial_fabricated.json"))
ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", HERE.parents[3] / "coupon-accuracy-assessment-20260913"))
# A stored campaign config of record (the four-edge physics-13 reference ran it), when the assessment tree is present.
CAMPAIGN = ASSESSMENT / "four-edge-physics-13" / "reference" / "case-07-fabricated" / "worker.json"


def legacy_copy(config):
    """The stored config as a legacy-policy writer spells it now: the same Model with the key written."""
    out = copy.deepcopy(config)
    model = {k: v for k, v in out["Model"].items() if k != CRACK_KEY}
    out["Model"] = apply_crack_policy(model, "legacy")
    return out


class ConfigDifferencesCrackKeyTest(unittest.TestCase):
    def check_stored(self, path):
        before = hashlib.sha256(path.read_bytes()).hexdigest()
        reference = json.loads(path.read_text())
        self.assertNotIn(CRACK_KEY, reference["Model"], f"{path} is not a pre-policy stored reference")
        derived = legacy_copy(reference)
        self.assertIs(derived["Model"][CRACK_KEY], True)
        # true vs absent: equal, in both argument orders.
        self.assertEqual(case_inputs.config_differences(derived, reference), [])
        self.assertEqual(case_inputs.config_differences(reference, derived), [])
        # Only that pair: false vs absent, true vs false, a non-bool value vs absent still differ.
        off = copy.deepcopy(derived)
        off["Model"][CRACK_KEY] = False
        self.assertEqual(case_inputs.config_differences(off, reference),
                         [f"/Model/{CRACK_KEY}: only in derived"])
        self.assertEqual(case_inputs.config_differences(reference, off),
                         [f"/Model/{CRACK_KEY}: only in reference"])
        self.assertEqual(case_inputs.config_differences(derived, off), [f"/Model/{CRACK_KEY}: True vs False"])
        one = copy.deepcopy(derived)
        one["Model"][CRACK_KEY] = 1
        self.assertEqual(case_inputs.config_differences(one, reference), [f"/Model/{CRACK_KEY}: only in derived"])
        # Another difference beside the key is still reported.
        other = copy.deepcopy(derived)
        other["Model"]["L0"] = 2.0 * reference["Model"]["L0"]
        self.assertEqual(len(case_inputs.config_differences(other, reference)), 1)
        # The stored reference is unchanged (read only).
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), before)

    def test_in_tree_stored_references(self):
        self.assertGreaterEqual(len(STORED), 3)
        for path in STORED:
            with self.subTest(path=str(path.relative_to(HERE))):
                self.check_stored(path)

    @unittest.skipUnless(CAMPAIGN.is_file(), f"needs the stored campaign config {CAMPAIGN}")
    def test_stored_campaign_worker_config(self):
        self.check_stored(CAMPAIGN)


if __name__ == "__main__":
    unittest.main()
