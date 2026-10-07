# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""qualify/refit_cost_model.py: the kept device model (cost-model-device-20260922.json) is
the refit of the recorded 2026-09-22 device library qualification (decision 64a),
conservative on every stage of every coupon of that run, with the previous (physics-11,
b28 at b = 6) model kept and bound by digest; the committed cost-model.json is the
decision-457 (2) / 458 measured refit of the 30 stage-2 coupons (qualify/stage2-20261006:
mean rates with per-stage safety factors from the residuals, the reducer node-used line
refit, the measured policy factors, the old-vs-new cap table), with the device model kept
and bound by digest; the parsers of Palace's reports."""
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
STAGE2_RECORDS = HERE / "qualify" / "stage2-20261006"
PREVIOUS_SHA256 = "22ee226d91c081f9b9275b020e7f68220941f27b32829245d54344da8e8f5ee9"   # the model PBS 47080-47106 were planned with
DEVICE_SHA256 = "f695da1b0aad66322a2c081e3ca0e3e46d94649b6741b47cab2e1a10996a91a0"     # the model every stage-2 qualify run planned with


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
        # The ReducerPeakStreaming block (USER decision 2026-09-22 (B), calibrated on
        # library-device-thin-01's node-used peaks) is carried, not refit: the committed model
        # is the refit of the recorded run with that block carried.
        cls.model, cls.record = refit_cost_model.refit(RECORDS / "library-qualification.json", previous_path=cls.previous,
                                                       previous_kept=cls.previous, streaming_path=estimate_stages.DEVICE_COST_MODEL)
        cls.committed = json.loads(estimate_stages.DEVICE_COST_MODEL.read_text())

    def test_previous_model_is_kept_byte_for_byte(self):
        self.assertEqual(refit_cost_model.sha256(self.previous), PREVIOUS_SHA256)
        self.assertEqual(self.committed["Previous"], {"Path": self.previous.name, "SHA256": PREVIOUS_SHA256,
                                                      "Provenance": self.model["Previous"]["Provenance"]})
        previous = estimate_stages.load_cost_model(self.previous)
        self.assertEqual((previous["MeasuredBlockSize"], previous["Stages"]["p4"]["Sources"]), (6, 80))

    def test_committed_model_is_the_refit_of_the_recorded_run(self):
        self.assertEqual(refit_cost_model.sha256(estimate_stages.DEVICE_COST_MODEL), DEVICE_SHA256)
        self.assertEqual(self.committed, self.model)
        # Without a streaming source the physics-11 previous model yields a model without the
        # block (the pre-streaming reducer peak term), otherwise identical.
        without, _ = refit_cost_model.refit(RECORDS / "library-qualification.json", previous_path=self.previous,
                                            previous_kept=self.previous)
        self.assertNotIn("ReducerPeakStreaming", without)
        self.assertEqual({k: v for k, v in self.model.items() if k != "ReducerPeakStreaming"}, without)

    def test_streaming_block_is_carried_unchanged_and_bound_to_the_frozen_executable(self):
        streaming = self.committed["ReducerPeakStreaming"]
        self.assertEqual(streaming["Executable"], self.committed["Provenance"]["FrozenExecutable"])
        self.assertIn("carried", streaming["CarriedRule"])
        self.assertEqual((streaming["NodeBaselineGiB"], streaming["NodeUsedGiBPerMillionH1"],
                          streaming["NodeUsedGiBPerMillionH1PerSource"], streaming["SafetyFactor"]), (60.0, 2.2, 0.015, 1.5))
        # Calibrated on the thin run's nine reducers, the line covers every measured reducer of
        # the fabricated run this model is refit from by the same safety factor (18 stages).
        stages = [stage for coupon in self.committed["Provenance"]["Coupons"].values() for stage in coupon["Stages"].values()]
        self.assertEqual(len(stages), 18)
        for stage in stages:
            estimate_gib = (estimate_stages.streaming_reducer_peak_gb(self.committed, stage["H1"], stage["Sources"])
                            / self.committed["PalaceGBPerGiB"])
            self.assertGreaterEqual(estimate_gib, streaming["SafetyFactor"] * stage["ReducerNodeUsedGiB"], stage)
        with tempfile.TemporaryDirectory() as tmp:
            other = Path(tmp) / "other.json"
            other.write_text(json.dumps({"ReducerPeakStreaming": {**streaming, "Executable": "0" * 64}}))
            with self.assertRaises(refit_cost_model.RefitError):
                refit_cost_model.refit(RECORDS / "library-qualification.json", previous_path=self.previous,
                                       previous_kept=self.previous, streaming_path=other)
            other.write_text(json.dumps({}))
            with self.assertRaises(refit_cost_model.RefitError):
                refit_cost_model.refit(RECORDS / "library-qualification.json", previous_path=self.previous,
                                       previous_kept=self.previous, streaming_path=other)

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


@unittest.skipUnless(STAGE2_RECORDS.is_dir(), "the stage-2 qualify records are not present")
class MeasuredRefitTest(unittest.TestCase):
    """Decision 457 (2) / 458: the committed cost-model.json is refit_measured of the 30
    stage-2 coupons (coupons-a, coupons-bc, pair 5; 90 jobs), with the device model kept."""

    @classmethod
    def setUpClass(cls):
        cls.runs = sorted(STAGE2_RECORDS.glob("**/library-qualification.json"))
        cls.model, cls.record = refit_cost_model.refit_measured(
            cls.runs, previous_path=estimate_stages.DEVICE_COST_MODEL, previous_kept=estimate_stages.DEVICE_COST_MODEL,
            label="stage-2 2026-10-05/06 (coupons-a, coupons-bc, pair 5)")
        cls.committed = json.loads(estimate_stages.COST_MODEL.read_text())

    def test_committed_model_is_the_measured_refit_with_the_device_model_kept(self):
        self.assertEqual(len(self.runs), 30)
        self.assertEqual(self.committed, self.model)
        self.assertEqual(self.committed["Version"], 3)
        self.assertEqual(self.committed["Previous"]["Path"], estimate_stages.DEVICE_COST_MODEL.name)
        self.assertEqual(self.committed["Previous"]["SHA256"], DEVICE_SHA256)
        self.assertEqual(len(self.committed["Provenance"]["Coupons"]), 30)
        self.assertEqual(len(self.committed["Provenance"]["Runs"]), 30)
        self.assertEqual(sorted(self.committed["Provenance"]["FrozenExecutables"]),
                         ["8357b113c16646dd80b8262757a010b2b7671df62875289113a0d2963e3bc333",
                          "cc7c4091fa47bde8739fe57678d1b72a22a4aee51b0de8da18ccb9912a5aec88"])
        self.assertEqual(sorted(self.committed["Stages"]), ["p3", "p4", "p5"])
        estimate_stages.load_cost_model(estimate_stages.COST_MODEL)   # the closed-form check of the reference mesh

    def test_mean_rates_with_per_stage_safety_factors_cover_every_measured_stage(self):
        for order, stage in self.committed["Stages"].items():
            for kind in ("Worker", "Reducer"):
                self.assertGreaterEqual(stage["SafetyFactor"][kind], 1.0, (order, kind))
                self.assertLessEqual(stage["Residuals"][kind]["Max"], stage["SafetyFactor"][kind] + 1e-9, (order, kind))
                self.assertLess(stage["SafetyFactor"][kind], 1.6, (order, kind))   # the measured spread, not a constant
            # The mean rate sits inside the measured spread (not the largest, as the device refit).
            set_by = self.committed["Provenance"]["SetBy"][order]["SecondsPerPCGIteration"]
            self.assertEqual(set_by["Rule"], "mean")
            self.assertLess(set_by["Value"], set_by["Largest"][0])
            self.assertGreater(set_by["Value"], set_by["Statistics"]["Min"])
        local = self.committed["LocalEdge"]
        self.assertGreaterEqual(local["SafetyFactor"]["LocalEdge"], local["Residuals"]["LocalEdge"]["Max"] - 1e-9)
        self.assertTrue(self.committed["SelfCheck"]["ConservativeAtOwnPCG"])
        # estimate_stages applies the factors: a stage of the reference mesh at the reference
        # source count costs the mean rates x the safety factor.
        counts = self.committed["MeasuredMesh"]["EntityCounts"]
        p4 = self.committed["Stages"]["p4"]
        stage = estimate_stages.estimate_stage(self.committed, 4, p4["Sources"], counts, 48)
        self.assertAlmostEqual(stage["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"],
                               (p4["WorkerNonSourceSeconds"] + p4["Sources"] * p4["MeanPerSourceSeconds"]) * p4["SafetyFactor"]["Worker"], places=6)
        self.assertAlmostEqual(stage["ReducerSecondsEstimate"], p4["ReducerPalaceSeconds"] * p4["SafetyFactor"]["Reducer"], places=6)
        self.assertEqual(stage["SafetyFactor"], p4["SafetyFactor"])
        # Memory stays the largest scaled measurement (fail-closed).
        for name in ("WorkerPalacePeakGB", "ReducerPalacePeakGB"):
            self.assertEqual(self.committed["Provenance"]["SetBy"]["p4"][name]["Rule"], "largest")

    def test_reducer_node_used_line_is_refit_with_a_residual_safety_factor(self):
        streaming = self.committed["ReducerPeakStreaming"]
        self.assertEqual(len(streaming["Measured"]), 60)   # 30 coupons x (p3 + p4 + p5) reducers
        self.assertGreaterEqual(streaming["SafetyFactor"], streaming["Residuals"]["Max"] - 1e-9)
        self.assertLess(streaming["SafetyFactor"], 1.5)   # the previous constant
        for row in streaming["Measured"]:
            estimate_gib = estimate_stages.streaming_reducer_peak_gb(self.committed, row["H1"], row["Sources"]) / self.committed["PalaceGBPerGiB"]
            self.assertGreaterEqual(estimate_gib, row["NodeUsedGiB"] - 1e-6, row["Label"])
        # The per-source term is the resident archived potentials: ~8 bytes per H1 DOF per source
        # (0.0075 GiB per million H1 per source) within a factor of two.
        self.assertGreater(streaming["NodeUsedGiBPerMillionH1PerSource"], 0.005)
        self.assertLess(streaming["NodeUsedGiBPerMillionH1PerSource"], 0.015)
        self.assertEqual(streaming["Previous"]["SafetyFactor"], 1.5)

    def test_no_measured_stage_of_record_exceeds_its_new_cap_and_the_pessimism_is_stated(self):
        caps = self.record["CapCheck"]
        self.assertEqual(caps["StagesBelowNewCap"], [])
        self.assertGreaterEqual(caps["MinimumMarginNewCapOverWall"], 1.0)
        self.assertEqual(caps["MarginStatistics"]["Count"], len(caps["Rows"]))
        self.assertTrue(all(row["CapOfRecordSeconds"] for row in caps["Rows"]))   # every stored stage has its cap of record
        self.assertEqual(self.committed["SelfCheck"]["CapCheck"]["MinimumMarginNewCapOverWall"], caps["MinimumMarginNewCapOverWall"])
        pessimism = self.committed["SelfCheck"]["Pessimism"]
        self.assertLess(pessimism["ThisModel"]["Max"], pessimism["PreviousModel"]["Min"])
        self.assertGreater(pessimism["ThisModel"]["Min"], 1.0)
        # The measured policy: the largest PCG factor covers the largest measured coupon-mean
        # ratio with headroom; the margin covers the measured runner overhead.
        self.assertEqual(self.committed["PCGFactors"], [1.0, 1.3, 1.5])
        self.assertEqual((self.committed["PreflightAndMarginFactor"], self.committed["PreflightSeconds"]), (1.10, 60))
        self.assertLess(self.committed["PolicyRule"]["MeasuredOverheads"]["Max"], 1.10)
        self.assertIn("previous [1.0, 1.5, 2.0]", self.committed["PolicyRule"]["PCGFactors"])


if __name__ == "__main__":
    unittest.main()
