# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""qualify: the gate evaluation on the stored physics-11 (four-edge) and gallery-06b
(three-edge) response CSVs reproduces the RESULTS.md class counts and passes; the
estimator reproduces the recorded gallery-06b stage estimate; the source locations and
classes are the recorded ones; the plan caps, the control choice and the verdict rules."""
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
import build_plan  # noqa: E402
import classify_sources  # noqa: E402
import compare_matrices  # noqa: E402
import estimate_stages  # noqa: E402
import gates  # noqa: E402
import locate_sources  # noqa: E402
import ma_ms_offsets  # noqa: E402
import p_sequence  # noqa: E402
import reference_campaign  # noqa: E402
import summarize_cost  # noqa: E402

ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", HERE.parents[3] / "coupon-accuracy-assessment-20260913"))
# The two recorded campaigns: (directory, graded_v2 input key, stage prefix, controls, RESULTS.md counts).
CAMPAIGNS = {
    "four-edge-physics-11": {"Key": "07", "Prefix": "va", "Controls": [1, 7, 23, 26, 34, 35, 48, 80], "Free": 60,
                             "E1": 60, "SA": (35, 44, 58), "MS": (51, 60, 60), "MA": (32, 55, 60), "Strongest": (15, 20),
                             "NodeHours": 0.508},
    "gallery-physics-06b": {"Key": "06", "Prefix": "g06b", "Controls": [14, 25, 26, 33, 52, 99, 133, 134], "Free": 95,
                            "E1": 93, "SA": (47, 69, 95), "MS": (77, 94, 95), "MA": (68, 91, 95), "Strongest": (19, 20),
                            "NodeHours": 0.974},
}


def campaign_available(name):
    return (ASSESSMENT / name / "results" / "main" / "status.json").is_file()


class StoredCampaign:
    """The recorded campaign's inputs, reference and results as the qualify modules read them."""

    def __init__(self, name):
        self.name = name
        self.spec = CAMPAIGNS[name]
        self.root = ASSESSMENT / name
        self.inputs = self.root / "reference" / f"inputs-{self.spec['Key']}"
        self.reference = self.root / "reference" / f"case-{self.spec['Key']}-fabricated" / "reducer"
        self.zero_trace = json.loads((self.inputs / "basis-contract.json").read_text())["ZeroTraceIndices"]
        self.results = self.root / "results" / "main"

    def stage(self, suffix):
        return self.results / f"{self.spec['Prefix']}-{suffix}" / "reducer"

    def locations_and_classes(self, out_dir):
        etch = self.inputs / "retained-etch.csv"
        rows, _ = locate_sources.locate_directory(self.inputs / "traces", self.inputs / "plan-view-boundary.csv",
                                                  Path(out_dir) / "source-locations.csv",
                                                  retained_etch=str(etch) if etch.is_file() else None)
        locations = {row["index"]: {key: str(value) for key, value in row.items()} for row in rows}
        return locations, classify_sources.classify_all(locations, self.zero_trace)


@unittest.skipUnless(all(campaign_available(name) for name in CAMPAIGNS), "the assessment campaigns are not available")
class StoredCampaignGateTest(unittest.TestCase):
    """The verdicts on the stored CSVs reproduce the RESULTS.md class counts."""

    def evaluate(self, campaign):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(tmp, True))
        locations, classes = campaign.locations_and_classes(tmp)
        comparison = compare_matrices.compare(campaign.reference, campaign.stage("p4"),
                                              zero_trace_indices=campaign.zero_trace, locations=locations)
        summary = compare_matrices.summary_record(comparison)
        sequence = p_sequence.p_sequence({"low": campaign.stage("p3-control"), "main": campaign.stage("p4"),
                                          "high": campaign.stage("p5-control"), "ref": campaign.reference},
                                         campaign.spec["Controls"])
        table, digest = gates.load_gates()
        record = gates.evaluate(table, comparison=summary, classes={i: name for i, (name, _) in classes.items()},
                                ref_pma=ma_ms_offsets.reference_p_ma(campaign.reference), p_sequence_summary=sequence,
                                reference_order=4, gates_sha256=digest)
        return record, summary, classes

    def check_counts(self, name):
        campaign = StoredCampaign(name)
        spec = campaign.spec
        record, summary, classes = self.evaluate(campaign)
        self.assertEqual(record["Verdict"], gates.VERDICT_PASSED, record["Reason"])
        self.assertEqual(record["ReferenceAnchor"], "vs p4 anchor")
        self.assertEqual(len(record["GatesSHA256"]), 64)
        self.assertEqual(len(summary["Free"]), spec["Free"])
        energy = record["Gates"]["E"]
        self.assertEqual((energy["AllFree"]["within_1pct"], energy["AllFree"]["n"]), (spec["E1"], spec["Free"]))
        self.assertEqual(energy["WithinBound"], energy["WideSources"])
        self.assertEqual(energy["FailingSources"], [])
        sa = record["Gates"]["p_SA"]["Free"]
        self.assertEqual((sa["within_1pct"], sa["within_2pct"], sa["within_5pct"]), spec["SA"])
        ms = record["Gates"]["p_MS"]["Free"]
        self.assertEqual((ms["within_1pct"], ms["within_2pct"], ms["within_5pct"]), spec["MS"])
        ma = record["Gates"]["p_MA"]
        self.assertEqual((ma["Free"]["within_1pct"], ma["Free"]["within_2pct"], ma["Free"]["within_5pct"]), spec["MA"])
        self.assertEqual((ma["Strongest"]["within_1pct"], ma["Strongest"]["within_2pct"]), spec["Strongest"])
        self.assertEqual(len(ma["StrongestSources"]), 20)
        self.assertTrue(record["Gates"]["PSequenceControls"]["Passed"])
        self.assertEqual(sorted(int(i) for i in record["Gates"]["PSequenceControls"]["Controls"]), spec["Controls"])
        # The classes partition the sources; every wide-class source is in the free view.
        by_class = {}
        for i, (name_, _) in classes.items():
            by_class.setdefault(name_, []).append(i)
        self.assertEqual(sum(len(v) for v in by_class.values()), len(summary["Sources"]))
        self.assertEqual(len(by_class[classify_sources.CLASS_ZERO_TRACE]), len(campaign.zero_trace))
        for name_ in classify_sources.WIDE_CLASSES:
            for i in by_class.get(name_, []):
                self.assertIn(i, summary["Free"])
        return record, classes

    def test_four_edge_physics_11_reproduces_results_counts(self):
        self.check_counts("four-edge-physics-11")

    def test_gallery_06b_reproduces_results_counts_and_class_table(self):
        record, classes = self.check_counts("gallery-physics-06b")
        counts = {}
        for name, _ in classes.values():
            counts[name] = counts.get(name, 0) + 1
        # RESULTS.md per-class table: 26 wide bottom/top, 8 box corners, 7 metal-top ring, 24
        # substrate/trench rings, 12 junction columns, 6 junction rings, 12 narrow hats (the
        # 8 junction-adjacent slivers and the 4 isolated slivers 0.5 um from a junction).
        self.assertEqual(counts[classify_sources.CLASS_WIDE_FACES], 26)
        self.assertEqual(counts[classify_sources.CLASS_BOX_CORNERS], 8)
        self.assertEqual(counts[classify_sources.CLASS_WIDE_METAL_TOP], 7)
        self.assertEqual(counts[classify_sources.CLASS_WIDE_SUBSTRATE], 24)
        self.assertEqual(counts[classify_sources.CLASS_JUNCTION_COLUMNS], 12)
        self.assertEqual(counts[classify_sources.CLASS_JUNCTION_RINGS], 6)
        self.assertEqual(counts[classify_sources.CLASS_NARROW_JUNCTION], 12)
        self.assertEqual(counts[classify_sources.CLASS_ZERO_TRACE], 40)
        self.assertEqual(classes[25][0], classify_sources.CLASS_NARROW_JUNCTION)
        self.assertLess(classes[25][1], 0.1)

    def test_gallery_06_before_the_sliver_fix_fails_the_ma_gate(self):
        """gallery-physics-06 (the mesh before decision 44) is the recorded negative: the
        isolated slivers 25 / 133 at MA -6.9 / -7.3% break 'every free source within 5%'
        (26 / 134 at -3.4 / -4.5% stay within it)."""
        run = ASSESSMENT / "gallery-physics-06" / "results" / "main" / "g06-p4" / "reducer"
        if not (run / "domain-response-matrix.csv").is_file():
            self.skipTest("gallery-physics-06 results not available")
        campaign = StoredCampaign("gallery-physics-06b")
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(tmp, True))
        locations, classes = campaign.locations_and_classes(tmp)
        summary = compare_matrices.summary_record(compare_matrices.compare(
            campaign.reference, run, zero_trace_indices=campaign.zero_trace, locations=locations))
        table, digest = gates.load_gates()
        record = gates.evaluate(table, comparison=summary, classes={i: name for i, (name, _) in classes.items()},
                                ref_pma=ma_ms_offsets.reference_p_ma(campaign.reference), reference_order=4, gates_sha256=digest)
        self.assertEqual(record["Verdict"], gates.VERDICT_FAILED)
        self.assertFalse(record["GatesPassed"]["p_MA"])
        self.assertTrue(record["GatesPassed"]["E"])   # the failing E sources are the narrow class
        self.assertEqual(record["Gates"]["E"]["AllFree"]["within_1pct"], 91)
        self.assertEqual(set(record["Gates"]["p_MA"]["FreeFailing"]), {25, 133})

    def test_locations_reproduce_the_recorded_gallery_06b_table(self):
        campaign = StoredCampaign("gallery-physics-06b")
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(tmp, True))
        locations, _ = campaign.locations_and_classes(tmp)
        recorded = compare_matrices.read_locations(campaign.root / "comparison" / "source-locations.csv")
        self.assertEqual(set(locations), set(recorded))
        for i, row in recorded.items():
            for key in ("lateral", "box_corner_3d", "adjacent_surface_at_z0", "metal_meets_cut_here",
                        "metal_edge_junction_here", "junction_distance_um", "near_trench_cut_junction"):
                self.assertEqual(locations[i][key], row[key], (i, key))

    def test_estimate_reproduces_the_recorded_gallery_06b_stage_estimate(self):
        recorded = json.loads((ASSESSMENT / "gallery-physics-06b" / "preflight" / "stage-estimate.json").read_text())
        counts = {"Vertices": 428366, "Edges": 2489399, "TriangleFaces": 3630370, "QuadFaces": 276012,
                  "Tetrahedra": 1659373, "Prisms": 172692, "Pyramids": 13284}
        estimate = estimate_stages.estimate(counts, [("p4-135", 4, 135), ("p5-8", 5, 8), ("p3-8", 3, 8)],
                                            local_edge=("local-edge", 4, 8))
        self.assertTrue(estimate["FitsOneJob"])
        for name in ("p4-135", "p5-8", "p3-8", "local-edge"):
            for factor in ("1.0", "1.5", "2.0"):
                self.assertAlmostEqual(estimate["Stages"][name]["ByPCGFactor"][factor]["StageSecondsEstimate"],
                                       recorded["Stages"][name]["ByPCGFactor"][factor]["StageSecondsEstimate"], delta=1.0)
        for factor in ("1.0", "2.0"):
            self.assertAlmostEqual(estimate["JobSecondsEstimateByPCGFactor"][factor],
                                   recorded["JobSecondsEstimateByPCGFactor"][factor], delta=2.0)
        self.assertAlmostEqual(estimate["MaxPalacePeakGBEstimate"], recorded["MaxPalacePeakGBEstimate"], delta=0.5)
        self.assertEqual(estimate["Mesh"]["H1ByOrder"]["p4"], 24844050)

    def test_cost_summary_reproduces_the_node_hours_of_the_suite_table(self):
        for name, spec in CAMPAIGNS.items():
            campaign = StoredCampaign(name)
            status = json.loads((campaign.results / "status.json").read_text())
            sources = 80 if name == "four-edge-physics-11" else 135
            summary = summarize_cost.summarize(status, full_sources=sources, nodes=1)
            main = summary["Stages"][f"{spec['Prefix']}-p4"]
            self.assertAlmostEqual(main["NodeHours"], spec["NodeHours"], places=3)
            self.assertEqual(len(main["Sources"]), sources)
            self.assertGreater(summary["JobNodeHours"], main["NodeHours"])

    def test_reference_campaign_resolves_by_basis_contract_digest(self):
        campaign = StoredCampaign("four-edge-physics-11")
        digest = reference_campaign.sha256(campaign.inputs / "basis-contract.json")
        resolved = reference_campaign.resolve(campaign.root / "reference", digest)
        self.assertEqual(resolved["Key"], "07")
        self.assertEqual(resolved["ReferenceOrder"], 4)
        self.assertEqual(resolved["LinearTol"], 1e-10)
        self.assertEqual(len(resolved["Sources"]), 80)
        self.assertEqual(resolved["ZeroTraceIndices"], campaign.zero_trace)
        self.assertEqual(resolved["Results"], str(campaign.reference))
        self.assertIn("worker.json", resolved["Config"])
        self.assertIsNone(reference_campaign.resolve(campaign.root / "reference", "0" * 64))
        self.assertIsNone(reference_campaign.resolve(None, digest))


class GateRuleTest(unittest.TestCase):
    """The verdict rules on synthetic offsets (no campaign needed)."""

    def synthetic(self, e=0.002, ma=0.005, ms=0.004, sa=0.01, n=30):
        free = list(range(1, n + 1))
        per_source = {str(i): {"E_rel": e, "p_MA_rel": ma, "p_MS_rel": ms, "p_SA_rel": sa} for i in free}
        comparison = {"PerSource": per_source, "Free": free, "Sources": free, "ZeroTrace": []}
        classes = {i: classify_sources.CLASS_WIDE_FACES for i in free}
        ref_pma = {i: 1.0 / i for i in free}
        return comparison, classes, ref_pma

    def sequence(self, step=0.001):
        return {1: {name: {"seq": {"d_low": 3 * step, "d_high": step, "r": 1 / 3, "p_inf": None}}
                    for name in p_sequence.OBSERVABLES}}

    def test_passed_failed_and_pending(self):
        table, digest = gates.load_gates()
        comparison, classes, ref_pma = self.synthetic()
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma,
                                p_sequence_summary=self.sequence(), reference_order=4, gates_sha256=digest)
        self.assertEqual(record["Verdict"], gates.VERDICT_PASSED)
        self.assertEqual(record["ReferenceAnchor"], "vs p4 anchor")
        self.assertEqual(len(record["Gates"]["p_MA"]["StrongestSources"]), 20)
        # One wide source beyond 1% in E fails the E gate only.
        comparison["PerSource"]["7"]["E_rel"] = 0.011
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma,
                                p_sequence_summary=self.sequence(), reference_order=4, gates_sha256=digest)
        self.assertEqual(record["Verdict"], gates.VERDICT_FAILED)
        self.assertEqual(record["Gates"]["E"]["FailingSources"], [7])
        self.assertTrue(record["GatesPassed"]["p_MA"])
        # The same source in a narrow class does not count against the wide-class E gate.
        classes[7] = classify_sources.CLASS_NARROW_ISOLATED
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma,
                                p_sequence_summary=self.sequence(), reference_order=4, gates_sha256=digest)
        self.assertTrue(record["GatesPassed"]["E"])
        self.assertEqual(record["Gates"]["E"]["AllFree"]["within_1pct"], 29)
        # A strongest-MA source beyond 2% fails p_MA; a p-step beyond the bound fails the controls.
        comparison["PerSource"]["1"]["p_MA_rel"] = 0.021
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma,
                                p_sequence_summary=self.sequence(step=0.02), reference_order=4, gates_sha256=digest)
        self.assertFalse(record["GatesPassed"]["p_MA"])
        self.assertEqual(record["Gates"]["p_MA"]["StrongestFailing"], [1])
        self.assertFalse(record["GatesPassed"]["PSequenceControls"])
        self.assertEqual(record["Gates"]["PSequenceControls"]["Failing"][0]["Observable"], "E")
        # Without a reference only the controls are evaluated: PendingQualification, never Passed.
        record = gates.evaluate(table, comparison=None, p_sequence_summary=self.sequence(), gates_sha256=digest)
        self.assertEqual(record["Verdict"], gates.VERDICT_PENDING)
        self.assertEqual(set(record["Gates"]), {"PSequenceControls"})
        self.assertTrue(record["GatesPassed"]["PSequenceControls"])

    def test_strongest_set_is_every_free_source_when_fewer_than_twenty(self):
        table, digest = gates.load_gates()
        comparison, classes, ref_pma = self.synthetic(n=12)
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma, reference_order=5,
                                gates_sha256=digest)
        self.assertEqual(len(record["Gates"]["p_MA"]["StrongestSources"]), 12)
        self.assertEqual(record["ReferenceAnchor"], "vs p5 anchor")

    def test_frozen_gate_table_is_the_confirmed_one(self):
        table, digest = gates.load_gates()
        self.assertEqual(table["Gates"]["E"]["MaximumAbsoluteRelativeOffset"], 0.01)
        self.assertEqual(table["Gates"]["p_MA"]["StrongestCount"], 20)
        self.assertEqual(table["Gates"]["p_MA"]["MaximumAbsoluteRelativeOffsetStrongest"], 0.02)
        self.assertEqual(table["Gates"]["p_MA"]["MaximumAbsoluteRelativeOffsetFree"], 0.05)
        self.assertEqual(table["Gates"]["p_MS"]["MaximumAbsoluteRelativeOffsetFree"], 0.05)
        self.assertAlmostEqual(table["Gates"]["p_SA"]["MinimumFractionWithin2Percent"], 2 / 3)
        self.assertEqual(table["Gates"]["p_SA"]["MinimumFractionWithin5Percent"], 0.9)
        self.assertEqual(table["Gates"]["PSequenceControls"]["MaximumAbsoluteEnergyStep"], 0.01)
        self.assertEqual(table["WideClasses"], list(classify_sources.WIDE_CLASSES))
        self.assertEqual(digest, reference_campaign.sha256(gates.GATES_FILE))


class ControlsPlanAndEstimateRuleTest(unittest.TestCase):
    def test_controls_by_class_cycle_over_the_classes(self):
        classes = {1: (classify_sources.CLASS_WIDE_FACES, 1.0), 2: (classify_sources.CLASS_WIDE_FACES, 1.0),
                   3: (classify_sources.CLASS_BOX_CORNERS, 1.0), 4: (classify_sources.CLASS_JUNCTION_RINGS, 1.0),
                   5: (classify_sources.CLASS_ZERO_TRACE, 1.0), 6: (classify_sources.CLASS_NARROW_JUNCTION, 0.05)}
        free = [1, 2, 3, 4, 6]
        self.assertEqual(classify_sources.choose_controls(classes, 4, free), [1, 3, 4, 6])
        self.assertEqual(classify_sources.choose_controls(classes, 5, free), [1, 2, 3, 4, 6])
        with self.assertRaises(ValueError):
            classify_sources.choose_controls(classes, 6, free)

    def test_cap_rule_and_stage_order(self):
        stage = {"ByPCGFactor": {"1.0": {"WorkerSecondsEstimate": 1000.0, "StageSecondsEstimate": 1400.0},
                                 "2.0": {"WorkerSecondsEstimate": 2000.0, "StageSecondsEstimate": 2400.0}},
                 "ReducerSecondsEstimate": 400.0}
        self.assertEqual(build_plan.stage_caps(stage, "worker", ["1.0", "2.0"], 20000), (4200, 1200))
        self.assertEqual(build_plan.stage_caps(stage, "reducer", ["1.0", "2.0"], 20000), (900, 600))
        self.assertEqual(build_plan.stage_caps(stage, "local-edge", ["1.0", "2.0"], 20000), (4800, 1500))
        self.assertEqual(build_plan.stage_caps(stage, "worker", ["1.0", "2.0"], 3000), (3000, 1200))
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        estimate = {"Stages": {"p4-10": stage, "le": stage}, "JobSecondsEstimateByPCGFactor": {},
                    "JobSecondsEstimateWithPreflightAndMargin": {}, "MaxPalacePeakGBEstimate": 1.0, "FitsOneJob": True,
                    "Decision": "fits"}
        plan = build_plan.build_plan(case_id="c", remote_case_root="/r/run/c", mesh={"Remote": "/r/m.msh", "SHA256": "a" * 64, "Local": "/l"},
                                     stage_layout=[{"Prefix": "x-p4", "Kind": "response", "EstimateKey": "p4-10"},
                                                   {"Prefix": "x-p4-local-edge", "Kind": "local-edge", "EstimateKey": "le"}],
                                     estimate=estimate, config_digests={"x-p4": {"worker.json": "b" * 64, "reducer.json": "c" * 64},
                                                                        "x-p4-local-edge": {"config.json": "d" * 64}},
                                     trace_pins={"/r/run/c/inputs/traces/basis-0001.csv": "e" * 64}, profile=profile,
                                     binary="/r/bin", binary_sha256="f" * 64, mpiexec="/r/mpi", purpose="test", factors=["1.0", "2.0"])
        self.assertEqual([s["Name"] for s in plan["Stages"]], ["x-p4-worker", "x-p4-reducer", "x-p4-local-edge"])
        self.assertEqual(plan["Stages"][1]["Requires"], ["x-p4-worker"])
        self.assertEqual(plan["Stages"][0]["Environment"]["PALACE_RESPONSE_ARCHIVE_ONLY"], "1")
        self.assertEqual(plan["Stages"][1]["Environment"]["PALACE_RESPONSE_BLOCK_SIZE"], "6")
        self.assertEqual(plan["Stages"][2]["Environment"], {})
        self.assertEqual(len(plan["PinnedSHA256"]), 5)
        self.assertEqual(plan["BinarySHA256"], "f" * 64)
        script = build_plan.render_job_script(profile=profile, remote_root="/r", remote_case_root="/r/run/c",
                                              runner="/r/run/run_stages.py", job_name="j", walltime_seconds=21600)
        self.assertIn("#PBS -l walltime=06:00:00", script)
        self.assertIn('python3 "/r/run/run_stages.py" "$D/plan.json"', script)

    def test_cost_model_closed_form_is_consistent_and_a_too_large_coupon_fails_closed(self):
        model = estimate_stages.load_cost_model()
        for name, check in model["ClosedFormCheck"].items():
            self.assertEqual(check["ClosedForm"], check["Measured"], name)
        counts = {key: value * 6 for key, value in model["MeasuredMesh"]["EntityCounts"].items()}
        estimate = estimate_stages.estimate(counts, [("p4-400", 4, 400)])
        self.assertFalse(estimate["FitsOneJob"])
        self.assertIn("does NOT fit", estimate["Decision"])
        with self.assertRaises(ValueError):
            estimate_stages.estimate(counts, [("p6-1", 6, 1)])

    def test_stage_estimate_scales_with_sources_and_dofs(self):
        model = estimate_stages.load_cost_model()
        counts = model["MeasuredMesh"]["EntityCounts"]
        same = estimate_stages.estimate_stage(model, 4, 80, counts)
        self.assertAlmostEqual(same["DOFRatioVsMeasured"], 1.0)
        self.assertAlmostEqual(same["ReducerSecondsEstimate"], model["Stages"]["p4"]["ReducerPalaceSeconds"], places=6)
        self.assertAlmostEqual(same["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"],
                               model["Stages"]["p4"]["WorkerNonSourceSeconds"] + 80 * model["Stages"]["p4"]["MeanPerSourceSeconds"], places=6)
        more = estimate_stages.estimate_stage(model, 4, 160, counts)
        self.assertGreater(more["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"], same["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"])
        self.assertGreater(more["ReducerPairs"], same["ReducerPairs"])

    def test_compare_offsets_and_free_view(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(tmp, True))
        for name, scale in (("ref", 1.0), ("run", 1.01)):
            (tmp / name).mkdir()
            (tmp / name / "domain-response-matrix.csv").write_text(
                "basis_i,basis_j,Q_ij (J)\n1,1,2.0\n1,2,0.5\n2,2,%g\n1,3,0.1\n2,3,0.2\n3,3,0.0\n" % (4.0 * scale))
            rows = ["interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_ij normal (J),Q_ij tangential (J),Q_total_ij (J),"
                    "Q_total_ij normal (J),Q_total_ij tangential (J)"]
            for interface in (1, 2, 3):
                for i, j, q in ((1, 1, 0.2), (1, 2, 0.05), (2, 2, 0.4 * scale), (3, 3, 0.0), (1, 3, 0.0), (2, 3, 0.0)):
                    rows.append(f"{interface},1,2e-06,{i},{j},{q},{q},{q},{q},{q},{q}")
            (tmp / name / "surface-response-matrix.csv").write_text("\n".join(rows) + "\n")
        comparison = compare_matrices.compare(tmp / "ref", tmp / "run", zero_trace_indices=[])
        self.assertEqual(comparison["ZeroTrace"], [3])     # zero diagonal energy
        self.assertEqual(comparison["Free"], [1, 2])
        self.assertAlmostEqual(comparison["PerSource"][2]["E_rel"], 0.01)
        self.assertAlmostEqual(comparison["PerSource"][2]["p_MA_rel"], 0.0)   # both scale: participation unchanged
        self.assertAlmostEqual(comparison["PerSource"][1]["E_rel"], 0.0)
        summary = compare_matrices.write_comparison(tmp / "ref", tmp / "run", tmp / "out", zero_trace_indices=[1])
        self.assertEqual(summary["ZeroTrace"], [1, 3])
        self.assertTrue((tmp / "out.md").is_file() and (tmp / "out.csv").is_file())
        self.assertTrue(math.isfinite(summary["PerSource"]["2"]["E_rel"]))


if __name__ == "__main__":
    unittest.main()
