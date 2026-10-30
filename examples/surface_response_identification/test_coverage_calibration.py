#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for coverage_calibration.py (DESIGN.md section 6): a planted linear law -> M1 / M2 selected, noise
only -> M0, the LOO coverage counting, the widening factor, the Wilson CI, bitwise reproducibility."""

import json
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_calibration as cc  # noqa: E402


def make_cells(n, law, noise=0.3, seed=7):
    """n windows with predictors (last step, SEC share, corner share) and d = law(pred) + a noise pattern: two
    outliers of +-noise (windows 0 / 1) and a +-0.001 pattern elsewhere, so a planted law's LOO coverage is
    exactly (n - 2) / n = 0.90 at n = 20 (the two outliers fall out of their own leave-one-out band)."""
    rng = random.Random(seed)
    cells = []
    for i in range(n):
        pred = {"last_step_pt": rng.uniform(0.0, 1.0), "sec_share": rng.uniform(0.0, 0.5), "corner_share": rng.uniform(0.0, 0.2)}
        eps = (noise if i == 0 else -noise) if i < 2 else 0.001 * ((i % 3) - 1)
        d = law(pred) + eps
        cells.append({"window": f"w{i}", "pred": {"SA": pred}, "d": {"SA": d}})
    return cells


class FitTest(unittest.TestCase):
    def test_quantile_and_median(self):
        self.assertEqual(cc.empirical_quantile([5, 1, 4, 2, 3, 6, 7, 8, 9, 10], 0.9), 9)
        self.assertEqual(cc.empirical_quantile([1.0], 0.9), 1.0)
        self.assertEqual(cc.empirical_quantile([], 0.9), 0.0)
        self.assertEqual(cc.median([3, 1, 2]), 2)
        self.assertEqual(cc.median([4, 1, 2, 3]), 2.5)

    def test_wilson_ci(self):
        lo, hi = cc.wilson_ci(18, 20)
        self.assertLess(lo, 0.9)
        self.assertGreater(hi, 0.9)
        self.assertAlmostEqual(lo, 0.7383369536731332, places=12)
        self.assertEqual(cc.wilson_ci(0, 0), [0.0, 1.0])

    def test_normal_equations_recover_a_plane(self):
        rows = [[1.0, x, y] for x, y in ((0, 0), (1, 0), (0, 1), (1, 1), (2, 3))]
        y = [2.0 + 3.0 * r[1] - 1.5 * r[2] for r in rows]
        beta = cc.solve_normal_equations(rows, y)
        self.assertAlmostEqual(beta[0], 2.0, places=10)
        self.assertAlmostEqual(beta[1], 3.0, places=10)
        self.assertAlmostEqual(beta[2], -1.5, places=10)
        self.assertIsNone(cc.solve_normal_equations([[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]], [1.0, 2.0, 3.0]))

    def test_planted_depth_lever_selects_m1(self):
        cells = make_cells(20, lambda p: -3.0 + 2.0 * p["last_step_pt"])
        entry = cc.select_model(cells, "SA")
        self.assertEqual(entry["model"], "M1")
        self.assertAlmostEqual(entry["coefficients"]["beta1"], 2.0, delta=0.2)
        self.assertEqual(entry["widened"], 1.0)
        self.assertAlmostEqual(entry["loo_coverage"]["coverage"], 0.9)
        self.assertEqual(entry["loo_coverage"]["n"], 20)
        self.assertLess(entry["loo_coverage"]["loo_half_width_pt"], entry["candidates"]["M0"]["loo_half_width_pt"])

    def test_planted_class_mix_selects_m2(self):
        cells = make_cells(20, lambda p: -2.0 + 4.0 * p["sec_share"] - 6.0 * p["corner_share"])
        entry = cc.select_model(cells, "SA")
        self.assertEqual(entry["model"], "M2")
        self.assertAlmostEqual(entry["coefficients"]["beta1"], 4.0, delta=0.3)

    def test_noise_only_selects_m0(self):
        cells = make_cells(20, lambda p: -3.0, seed=3)
        entry = cc.select_model(cells, "SA")
        self.assertEqual(entry["model"], "M0")
        self.assertAlmostEqual(entry["coefficients"]["beta0"], -3.0, delta=0.01)

    def test_loo_counting_and_widening(self):
        # 10 windows at -3 exactly plus two outliers: the LOO coverage of M0 counts them out; widening reaches 90 %
        cells = [{"window": f"w{i}", "pred": {"SA": {"last_step_pt": 0.0, "sec_share": 0.0, "corner_share": 0.0}},
                  "d": {"SA": -3.0 + 0.01 * (i % 3)}} for i in range(10)]
        cells.append({"window": "out1", "pred": cells[0]["pred"], "d": {"SA": -1.0}})
        cells.append({"window": "out2", "pred": cells[0]["pred"], "d": {"SA": -5.0}})
        loo = cc.loo_evaluate("M0", cells, "SA")
        self.assertEqual(loo["n"], 12)
        self.assertEqual(loo["covered"], 10)
        self.assertAlmostEqual(loo["coverage"], 10.0 / 12.0)
        self.assertGreater(loo["widen_k"], 1.0)
        entry = cc.select_model(cells, "SA")
        self.assertEqual(entry["model"], "M0")
        self.assertGreater(entry["widened"], 1.0)
        self.assertGreaterEqual(entry["loo_coverage"]["coverage"], 0.9)
        self.assertIn("WIDENED", entry["loo_coverage"]["note"])

    def test_bitwise_reproducible(self):
        cells = make_cells(15, lambda p: -3.0 + 1.0 * p["last_step_pt"], seed=11)
        a = json.dumps(cc.select_model(cells, "SA"), sort_keys=True)
        b = json.dumps(cc.select_model(list(cells), "SA"), sort_keys=True)
        self.assertEqual(a, b)

    def test_class_census(self):
        cells = [{"pred": {"SA": {"classes_share_length": {"IsolatedEdge": 0.5, "ConvexCorner": 0.06}},
                           "MS": {"classes_share_length": {"IsolatedEdge": 0.5}}, "MA_sharp": {"classes_share_length": {}}}},
                 {"pred": {"SA": {"classes_share_length": {"IsolatedEdge": 0.5, "ConvexCorner": 0.01}},
                           "MS": {"classes_share_length": {}}, "MA_sharp": {"classes_share_length": {}}}}]
        self.assertEqual(cc.class_census(cells), {"IsolatedEdge": 2, "ConvexCorner": 1})


if __name__ == "__main__":
    unittest.main()
