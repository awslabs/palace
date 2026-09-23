# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""qualify/refit_cost_model.py: the committed cost model is the refit of the recorded
2026-09-22 device library qualification (decision 64a), conservative on every stage of
every coupon of that run, with the previous (physics-11, b28 at b = 6) model kept and
bound by digest; the parsers of Palace's reports."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
import estimate_stages  # noqa: E402
import refit_cost_model  # noqa: E402

RECORDS = HERE / "qualify" / "device-library-20260922"
PREVIOUS_SHA256 = "22ee226d91c081f9b9275b020e7f68220941f27b32829245d54344da8e8f5ee9"   # the model PBS 47080-47106 were planned with


class ReportParsersTest(unittest.TestCase):

    def test_elapsed_time_report_average_column_and_linear_solve(self):
        text = ("Elapsed Time Report (s)           Min.        Max.        Avg.\n"
                "==============================================================\n"
                "Initialization                   0.535       1.827       1.808\n"
                "  Mesh Preprocessing             2.288       3.572       2.299\n"
                "Linear Solve                     2.238       5.284       3.677\n"
                "  Setup                          2.275       2.300       2.286\n"
                "  Preconditioner                46.377      50.692      48.875\n"
                "  Coarse Solve                   0.352       1.649       0.758\n"
                "Estimation                       0.125       0.252       0.156\n"
                "  Solve                          2.206       8.169       8.135\n"
                "Postprocessing                  22.446      28.476      22.531\n"
                "--------------------------------------------------------------\n"
                "Total                          101.191     101.359     101.288\n")
        timers = refit_cost_model.parse_elapsed_time_report(text)
        self.assertEqual(timers["Preconditioner"], 48.875)
        self.assertEqual(timers["Solve"], 8.135)     # the estimator's solve is not the linear solve
        self.assertNotIn("Total", timers)
        self.assertAlmostEqual(refit_cost_model.linear_solve_seconds(timers), 3.677 + 2.286 + 48.875 + 0.758)
        with self.assertRaises(refit_cost_model.RefitError):
            refit_cost_model.parse_elapsed_time_report("nothing")

    def test_palace_memory_figures_and_rounding_up(self):
        self.assertEqual(refit_cost_model.memory_gb("48.2G"), 48.2)
        self.assertAlmostEqual(refit_cost_model.memory_gb("512M"), 0.5)
        self.assertAlmostEqual(refit_cost_model.memory_gb("1.5T"), 1536.0)
        self.assertEqual(refit_cost_model.round_up(1.23451, 3), 1.235)
        self.assertEqual(refit_cost_model.round_up(2.0, 2), 2.0)
        self.assertEqual(refit_cost_model.archive_gb({"SizesBeforeDeletion": "51G\t/a/p4/archive\n4.2G\t/a/p5/archive\n"}), 51.0)
        self.assertIsNone(refit_cost_model.archive_gb({"SizesBeforeDeletion": ""}))


@unittest.skipUnless((RECORDS / "library-qualification.json").is_file(), "the 2026-09-22 device library records are not present")
class RecordedRefitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.previous = HERE / "qualify" / estimate_stages.PREVIOUS_COST_MODEL.name
        cls.model, cls.record = refit_cost_model.refit(RECORDS / "library-qualification.json", previous_path=cls.previous,
                                                       previous_kept=cls.previous)
        cls.committed = json.loads(estimate_stages.COST_MODEL.read_text())

    def test_previous_model_is_kept_byte_for_byte(self):
        self.assertEqual(refit_cost_model.sha256(self.previous), PREVIOUS_SHA256)
        self.assertEqual(self.committed["Previous"], {"Path": self.previous.name, "SHA256": PREVIOUS_SHA256,
                                                      "Provenance": self.model["Previous"]["Provenance"]})
        previous = estimate_stages.load_cost_model(self.previous)
        self.assertEqual((previous["MeasuredBlockSize"], previous["Stages"]["p4"]["Sources"]), (6, 80))

    def test_committed_model_is_the_refit_of_the_recorded_run(self):
        self.assertEqual(self.committed, self.model)

    def test_refit_passes_the_closed_form_self_check_and_measured_block_size_48(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "cost-model.json"
            path.write_text(json.dumps(self.model))
            loaded = estimate_stages.load_cost_model(path)
        for name, check in loaded["ClosedFormCheck"].items():
            self.assertEqual(check["ClosedForm"], check["Measured"], name)
        self.assertEqual(sorted(loaded["Stages"]), ["p3", "p4", "p5"])
        self.assertEqual(loaded["MeasuredBlockSize"], 48)
        self.assertEqual(loaded["Provenance"]["FrozenExecutable"], "170439c4a9fc5d5ce329310812055be5fb83a4a7f288024b57b3b83551cbe70b")
        self.assertEqual(loaded["MeasuredMesh"]["SHA256"],
                         loaded["Provenance"]["Coupons"][loaded["Provenance"]["ReferenceCase"]]["MeshSHA256"])
        # Kept calibrations (not re-measurable at one block size), stated as kept.
        self.assertEqual((loaded["ReducerEvaluationFraction"], loaded["ReducerResidentFieldGBPerMillionH1"]), (0.97, 0.068))
        self.assertIn("KEPT", loaded["ReducerEvaluationFractionRule"])
        self.assertIn("KEPT", loaded["ReducerResidentFieldGBRule"])

    def test_refit_is_conservative_on_every_stage_and_job_of_every_coupon(self):
        check = self.record["SelfCheck"]
        self.assertTrue(check["Conservative"])
        self.assertGreaterEqual(check["MinimumStageEstimateOverActual"], 1.0 - 1e-9)
        coupons = self.model["Provenance"]["Coupons"]
        self.assertEqual(len(coupons), 6)
        for case in coupons:
            jobs = check[case]["Jobs"]
            self.assertEqual(set(jobs), set(coupons[case]["Jobs"]), case)
            for job, record in jobs.items():
                for name, ratio in record["PerStage"].items():
                    self.assertGreaterEqual(ratio, 1.0 - 1e-9, (case, job, name))
                self.assertGreaterEqual(record["Estimate1xOverStageWall"], 1.0 - 1e-9, (case, job))
                self.assertGreater(record["EstimateWorstWithMarginOverJob"], 1.0, (case, job))
            self.assertGreater(check[case]["EstimateWorstWithMarginOverActualNodeSeconds"], 1.0, case)
            self.assertEqual(self.model["SelfCheck"]["EstimateWorstWithMarginOverActual"][case],
                             check[case]["EstimateWorstWithMarginOverActualNodeSeconds"])
        # Every coupon ran split (policy speed): the single-job decision is recorded, not exercised.
        self.assertEqual(set(self.model["SelfCheck"]["FitsOneJob"]), set(coupons))
        # The previous model was not conservative on this run (the reason for the refit).
        self.assertLess(self.model["SelfCheck"]["PreviousModel"]["MinimumStageEstimateOverActual"], 1.0)
        # The reference coupon sets at least one rate at exactly its measurement (ratio 1).
        reference = self.model["Provenance"]["ReferenceCase"]
        set_by = {item["SetBy"] for stage in self.model["Provenance"]["SetBy"].values() for item in stage.values()}
        self.assertIn(reference, set_by)
        self.assertTrue(set_by <= set(coupons) | {"previous model (scaled)"})

    def test_every_rate_covers_the_largest_scaled_measurement(self):
        for order, stage in self.model["Stages"].items():
            for case, coupon in self.model["Provenance"]["Coupons"].items():
                measured = coupon["Stages"].get(order)
                if measured is None:
                    continue
                ratio = stage["H1"] / measured["H1"]
                self.assertGreaterEqual(stage["SecondsPerPCGIteration"], measured["SecondsPerPCGIteration"] * ratio - 1e-9, (order, case))
                self.assertGreaterEqual(stage["MeanPCGIterations"], measured["MeanPCGIterations"] - 1e-9, (order, case))
                self.assertGreaterEqual(stage["MaxPCGIterations"], measured["MaxPCGIterations"], (order, case))
                self.assertGreaterEqual(stage["WorkerNonSourceSeconds"], measured["WorkerNonSourceSeconds"] * ratio - 1e-9, (order, case))
                self.assertGreaterEqual(stage["WorkerPalacePeakGB"], measured["WorkerPalacePeakGB"] * ratio - 1e-9, (order, case))
                self.assertGreaterEqual(stage["ReducerPalacePeakGB"], measured["ReducerPalacePeakGB"] * ratio - 1e-9, (order, case))
        local = self.model["LocalEdge"]
        for case, coupon in self.model["Provenance"]["Coupons"].items():
            ratio = local["H1"] / coupon["LocalEdge"]["H1"]
            self.assertGreaterEqual(local["PalaceTotalSeconds"], coupon["LocalEdge"]["WallSeconds"] * ratio - 1e-9, case)


if __name__ == "__main__":
    unittest.main()
