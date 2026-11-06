#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for coverage_calibration.py (DESIGN.md section 6): a planted linear law -> M1 / M2 selected, noise
only -> M0, the LOO coverage counting, the widening factor, the Wilson CI, bitwise reproducibility."""

import json
import os
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_calibration as cc  # noqa: E402
import coverage_report as cr  # noqa: E402


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
        # 12 windows at -3 exactly plus two outliers: the LOO coverage of M0 counts them out (12 / 14 = 0.857 >= 0.85); widening reaches 90 %
        cells = [{"window": f"w{i}", "pred": {"SA": {"last_step_pt": 0.0, "sec_share": 0.0, "corner_share": 0.0}},
                  "d": {"SA": -3.0 + 0.01 * (i % 3)}} for i in range(12)]
        cells.append({"window": "out1", "pred": cells[0]["pred"], "d": {"SA": -1.0}})
        cells.append({"window": "out2", "pred": cells[0]["pred"], "d": {"SA": -5.0}})
        loo = cc.loo_evaluate("M0", cells, "SA")
        self.assertEqual(loo["n"], 14)
        self.assertEqual(loo["covered"], 12)
        self.assertAlmostEqual(loo["coverage"], 12.0 / 14.0)
        self.assertGreater(loo["widen_k"], 1.0)
        entry = cc.select_model(cells, "SA")
        self.assertEqual(entry["model"], "M0")
        self.assertGreater(entry["widened"], 1.0)
        self.assertGreaterEqual(entry["loo_coverage"]["coverage"], 0.9)
        self.assertIn("WIDENED", entry["loo_coverage"]["note"])

    def test_uncalibrated_below_0_85(self):
        # bands v6 rule (5) (DESIGN-v3 s7; decision 1233 (6)): 10 windows + two outliers = 10 / 12 = 0.833 < 0.85 -> UNCALIBRATED, no widening
        cells = [{"window": f"w{i}", "pred": {"SA": {"last_step_pt": 0.0, "sec_share": 0.0, "corner_share": 0.0}},
                  "d": {"SA": -3.0 + 0.01 * (i % 3)}} for i in range(10)]
        cells.append({"window": "out1", "pred": cells[0]["pred"], "d": {"SA": -1.0}})
        cells.append({"window": "out2", "pred": cells[0]["pred"], "d": {"SA": -5.0}})
        entry = cc.select_model(cells, "SA")
        self.assertTrue(entry["uncalibrated"])
        self.assertEqual(entry["uncalibrated_below"], 0.85)
        self.assertNotIn("model", entry)
        self.assertIn("0.833", entry["reason"])
        self.assertEqual(sorted(entry["candidates"]), ["M0", "M1", "M2", "M3"])
        # the reader: an uncalibrated Type carries no band and labels LOW
        bands = {"families": {"ES": {"types": {"SA": entry}}}}
        row = cr.band_table(bands, "ES", {"types": {"SA": {"raw_share_total": 0.0, "classes_share": {}, "exact_share": 1.0, "raw_share": {}}}},
                            {"types": {"SA": {}}}, {}, [], False)
        self.assertFalse(row["SA"]["calibrated"])
        self.assertTrue(row["SA"]["uncalibrated"])
        label = cr.label_run({"types": {"SA": {"raw_share_total": 0.0, "classes_share": {}, "exact_share": 1.0, "raw_share": {}, "scp_owner_credit": None}}},
                             {"types": {"SA": {}}}, {"SA": {"status": "quotable"}}, [], {}, bands["families"]["ES"], False, "sharp")
        self.assertEqual(label["label"], "LOW")
        self.assertTrue(any("UNCALIBRATED Type" in r for r in label["reasons_low"]))

    def test_m3_corner_share_candidate(self):
        # bands v6 (decision 1233 (6) (b)): a planted corner-share law selects M3 (the smaller LOO half-width than M2's two levers)
        cells = make_cells(20, lambda p: -1.0 + 8.0 * p["corner_share"])
        entry = cc.select_model(cells, "SA")
        self.assertEqual(entry["model"], "M3")
        self.assertAlmostEqual(entry["coefficients"]["beta1"], 8.0, delta=0.5)   # the two +-0.3 outliers of make_cells pull the slope
        self.assertEqual(list(entry["coefficients"]), ["beta0", "beta1"])
        self.assertAlmostEqual(cr.predict_bias(entry, {"corner_share": 0.1}), entry["coefficients"]["beta0"] + 0.1 * entry["coefficients"]["beta1"])
        # the sensitivity row: a multi-source cell's aggregated last step swapped in -> an information re-fit
        cells[0]["pred"]["SA"]["aggregated"] = {"last_step_pt": 0.5, "raw_share_energy": 0.0, "n_sources": 4}
        cells[0]["pred"]["SA"]["source"] = 3
        sens = cc.aggregated_sensitivity(cells, "SA", entry)
        self.assertEqual(sens["multi_source_windows"], ["w0"])
        self.assertEqual(sens["model"], "M3")
        self.assertIsNone(cc.aggregated_sensitivity(cells, "SA", {"uncalibrated": True}))

    def test_c_is_fitted_with_m0_only(self):
        # DESIGN 5.5 / MAJOR-4: a planted depth lever on C must NOT select M1; the C band applied (beta0 +- H) is the one validated
        cells = make_cells(20, lambda p: 0.1 - 0.5 * p["last_step_pt"])
        for c in cells:
            c["pred"]["C"] = c["pred"].pop("SA")
            c["d"]["C"] = c["d"].pop("SA")
        entry = cc.select_model(cells, "C", models=("M0",))
        self.assertEqual(entry["model"], "M0")
        self.assertEqual(list(entry["candidates"]), ["M0"])
        self.assertEqual(list(entry["coefficients"]), ["beta0"])
        self.assertEqual(entry["coefficients"]["beta0"], cc.median([c["d"]["C"] for c in cells]))
        self.assertEqual(cc.select_model(cells, "C")["model"], "M1")  # the unrestricted selection would have taken the lever
        # the widening note discloses the unwidened LOO coverage (MINOR-2)
        cells = [{"window": f"w{i}", "pred": {"C": {"last_step_pt": 0.0, "sec_share": 0.0, "corner_share": 0.0}},
                  "d": {"C": 0.1 + 0.01 * (i % 3)}} for i in range(12)]
        cells += [{"window": "out1", "pred": cells[0]["pred"], "d": {"C": 1.0}}, {"window": "out2", "pred": cells[0]["pred"], "d": {"C": -1.0}}]
        entry = cc.select_model(cells, "C", models=("M0",))
        self.assertGreater(entry["widened"], 1.0)
        self.assertIn("unwidened LOO coverage is 12 / 14", entry["loo_coverage"]["note"])
        self.assertEqual(entry["loo_coverage"]["coverage_unwidened"], entry["candidates"]["M0"]["loo_coverage"])
        self.assertGreaterEqual(entry["loo_coverage"]["coverage"], 0.9)

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

    def test_v6_band_file_kept_beside(self):
        """data/coverage_bands_v6.json (fix round P1; DESIGN-v3 s7; decisions 1230 / 1233 (6) / 1236): the default band file until v6.1. The ES
        block (n = 28: the 20 v5 inputs + the 8 stage-E ES windows; W31 excluded, named) fits M0 on every Type (no candidate incl. M3
        reached 0.90 unwidened), its C band holds the measured +0.15..+0.33; the Maxwell block is UNCALIBRATED on every Type (the 12
        PROXY inputs: best unwidened LOO 10 / 12 = 0.833 < 0.85; 1236 (A): not carried from v5) until the ENERGY-form v6.1."""
        here = os.path.dirname(os.path.abspath(cr.__file__))
        v6 = cr.load_bands(os.path.join(here, "data", "coverage_bands_v6.json"))
        # the default moved to v6.1 at the P4 re-freeze (test_v6_1_band_file_is_the_default); v6 is kept beside
        self.assertEqual(os.path.basename(cr.DEFAULT_BANDS), "coverage_bands_v6.1.json")
        self.assertEqual((v6["version"], v6["models"], v6["uncalibrated_below"]), ("v6", ["M0", "M1", "M2", "M3"], 0.85))
        es = v6["families"]["ES"]
        self.assertEqual((es["n"], es["class_census"]["DifferentConductorGap"], es["class_census"]["ConvexCorner"]), (28, 5, 2))
        for t in ("SA", "MS", "MA", "C"):
            self.assertEqual(es["types"][t]["model"], "M0", t)
            self.assertGreaterEqual(es["types"][t]["loo_coverage"]["coverage"], 0.9, t)
            self.assertLess(es["types"][t]["loo_coverage"]["coverage_unwidened"], 0.9, t)
        c = es["types"]["C"]
        self.assertLess(c["coefficients"]["beta0"] - c["H_pt"] * c["widened"], 0.15)
        self.assertGreater(c["coefficients"]["beta0"] + c["H_pt"] * c["widened"], 0.32)
        self.assertEqual([e["window"] for e in v6["provenance"]["excluded"]], ["W31"])
        mx = v6["families"]["Maxwell"]
        self.assertTrue(mx["form"].startswith("PROXY"))
        for t in ("SA", "MS", "MA", "C"):
            self.assertTrue(mx["types"][t]["uncalibrated"], t)
            self.assertIn("0.833", mx["types"][t]["reason"])

    def test_v6_1_band_file_is_the_default(self):
        """data/coverage_bands_v6.1.json (fix round P4 re-freeze; REFREEZE-DESIGN s1 DATA-1b / s3 step 5; decisions 1236 (A) / 1364): the
        default band file. The ES block is v6's, bit-identical (the 28 cells unchanged; 7 pruned palace.json read through substitutes);
        the Maxwell block adds the three stage-E ENERGY-form cells (n 12 -> 15): SA / MS CALIBRATED (M0 widened to >= 0.90 LOO),
        MA_sharp / C stay UNCALIBRATED under rule (5) (0.733 / 0.786 < 0.85); the 1364 record lines are in the provenance."""
        here = os.path.dirname(os.path.abspath(cr.__file__))
        v6 = cr.load_bands(os.path.join(here, "data", "coverage_bands_v6.json"))
        v61 = cr.load_bands(os.path.join(here, "data", "coverage_bands_v6.1.json"))
        self.assertEqual(os.path.basename(cr.DEFAULT_BANDS), "coverage_bands_v6.1.json")
        self.assertEqual((v61["version"], v61["models"], v61["uncalibrated_below"]), ("v6.1", ["M0", "M1", "M2", "M3"], 0.85))
        self.assertEqual(v61["families"]["ES"], v6["families"]["ES"])
        mx = v61["families"]["Maxwell"]
        self.assertEqual(mx["n"], 15)
        self.assertTrue(mx["form"].startswith("ENERGY"))
        for t in ("SA", "MS"):
            self.assertEqual(mx["types"][t]["model"], "M0", t)
            self.assertFalse(mx["types"][t].get("uncalibrated"), t)
            self.assertGreaterEqual(mx["types"][t]["loo_coverage"]["coverage"], 0.9, t)
            self.assertLess(mx["types"][t]["loo_coverage"]["coverage_unwidened"], 0.9, t)
        for t in ("MA", "C"):
            self.assertTrue(mx["types"][t]["uncalibrated"], t)
        self.assertIn("0.733", mx["types"]["MA"]["reason"])
        self.assertIn("0.786", mx["types"]["C"]["reason"])
        self.assertEqual(len(v61["provenance"]["record_lines_1364"]), 5)
        self.assertEqual([e["window"] for e in v61["provenance"]["excluded"]], ["W31"])

    def test_v5_band_file_restates_v3s_fit(self):
        """data/coverage_bands_v5.json (covfix lane, decisions 1125 / 1126 (4) / 1128 / 1129; the default until v6): its fitted
        numbers and LOO coverage + CI per family x Type equal v3's (the predictors are untouched by the reader fixes); median_g
        moves (g reads the credited raw); the note says which numbers move and why."""
        here = os.path.dirname(os.path.abspath(cr.__file__))
        v5 = cr.load_bands(os.path.join(here, "data", "coverage_bands_v5.json"))
        v3 = cr.load_bands(os.path.join(here, "data", "coverage_bands_v3.json"))
        self.assertEqual(v5["version"], "v5")
        self.assertIn("IDENTICAL to v3", v5.get("note", ""))
        for family, fam in v3["families"].items():
            self.assertEqual(fam["class_census"], v5["families"][family]["class_census"])
            for t, e3 in fam["types"].items():
                e5 = v5["families"][family]["types"][t]
                for key in ("model", "coefficients", "H_pt", "widened", "n", "median_residual_pt", "n_energy_mode"):
                    self.assertEqual(e3[key], e5[key], f"{family} {t} {key}")
                self.assertEqual({k: v for k, v in e3["loo_coverage"].items() if k != "cells"},
                                 {k: v for k, v in e5["loo_coverage"].items() if k != "cells"}, f"{family} {t} loo")
                if t != "C":
                    self.assertNotEqual(e3["median_g"], e5["median_g"], f"{family} {t} median_g")

    def test_note_option_is_stamped(self):
        records = {"families": {}}
        tmp = tempfile.mkdtemp(prefix="covcal-")
        try:
            rec = os.path.join(tmp, "records.json")
            with open(rec, "w") as f:
                json.dump(records, f)
            out = os.path.join(tmp, "bands.json")
            cc.main(["--records", rec, "--output", out, "--tag", "t", "--note", "why the numbers move", "--quiet"])
            with open(out) as f:
                self.assertEqual(json.load(f)["note"], "why the numbers move")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
