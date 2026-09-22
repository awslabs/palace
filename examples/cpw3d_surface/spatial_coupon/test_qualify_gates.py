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
import case_inputs  # noqa: E402
import classify_sources  # noqa: E402
import compare_matrices  # noqa: E402
import estimate_stages  # noqa: E402
import gates  # noqa: E402
import locate_sources  # noqa: E402
import ma_ms_offsets  # noqa: E402
import p_sequence  # noqa: E402
import reference_campaign  # noqa: E402
import remote  # noqa: E402
import summarize_cost  # noqa: E402

ASSESSMENT = Path(os.environ.get("COUPON_ASSESSMENT_ROOT", HERE.parents[3] / "coupon-accuracy-assessment-20260913"))
# The recorded references postprocess MA / MS / SA at the indices 1 / 2 / 3 (read from their configs).
THREE_INTERFACES = {1: "MA", 2: "MS", 3: "SA"}
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
        self.interface_types = case_inputs.interface_types(json.loads((self.inputs / "spatial_fabricated.json").read_text()))

    def stage(self, suffix):
        return self.results / f"{self.spec['Prefix']}-{suffix}" / "reducer"

    def locations_and_classes(self, out_dir):
        etch = self.inputs / "retained-etch.csv"
        rows, _ = locate_sources.locate_directory(self.inputs / "traces", self.inputs / "plan-view-boundary.csv",
                                                  Path(out_dir) / "source-locations.csv",
                                                  signature_path=self.inputs / "mesh-signature.csv",
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
        self.assertEqual(campaign.interface_types, THREE_INTERFACES)
        comparison = compare_matrices.compare(campaign.reference, campaign.stage("p4"),
                                              zero_trace_indices=campaign.zero_trace, locations=locations,
                                              interface_types=campaign.interface_types)
        summary = compare_matrices.summary_record(comparison)
        sequence = p_sequence.p_sequence({"low": campaign.stage("p3-control"), "main": campaign.stage("p4"),
                                          "high": campaign.stage("p5-control"), "ref": campaign.reference},
                                         campaign.spec["Controls"], campaign.interface_types)
        table, digest = gates.load_gates()
        record = gates.evaluate(table, comparison=summary, classes={i: name for i, (name, _) in classes.items()},
                                ref_pma=ma_ms_offsets.reference_p_ma(campaign.reference, campaign.interface_types),
                                p_sequence_summary=sequence, reference_order=4, gates_sha256=digest)
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
            campaign.reference, run, zero_trace_indices=campaign.zero_trace, locations=locations,
            interface_types=campaign.interface_types))
        table, digest = gates.load_gates()
        record = gates.evaluate(table, comparison=summary, classes={i: name for i, (name, _) in classes.items()},
                                ref_pma=ma_ms_offsets.reference_p_ma(campaign.reference, campaign.interface_types),
                                reference_order=4, gates_sha256=digest)
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
        # The recorded estimate ran at the physics-11 model's measured block size 6 with the
        # pair-scaled reducer formula (before decision 62(1)): reproduced as recorded with
        # that model (kept as PREVIOUS_COST_MODEL since the decision-64a refit).
        model = {**estimate_stages.load_cost_model(estimate_stages.PREVIOUS_COST_MODEL), "ReducerEvaluationFraction": 0.0}
        estimate = estimate_stages.estimate(counts, [("p4-135", 4, 135), ("p5-8", 5, 8), ("p3-8", 3, 8)],
                                            local_edge=("local-edge", 4, 8), block_size=6, model=model)
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

    def test_undeclared_interface_is_not_applicable(self):
        """A reference that postprocesses MA and MS only (gallery case 10): the p_SA gate and the
        p_SA p-sequence observable are NotApplicable, recorded, never a failure; the declared
        interfaces are still gated; without the interface list everything is gated."""
        table, digest = gates.load_gates()
        comparison, classes, ref_pma = self.synthetic()
        for record in comparison["PerSource"].values():
            del record["p_SA_rel"]
        sequence = self.sequence()
        sequence[1]["p_SA"]["seq"] = {"d_low": None, "d_high": None, "r": None, "p_inf": None}
        gated = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma, p_sequence_summary=sequence,
                               reference_order=5, gates_sha256=digest)
        self.assertEqual(gated["Verdict"], gates.VERDICT_FAILED)
        self.assertFalse(gated["GatesPassed"]["p_SA"])
        self.assertFalse(gated["GatesPassed"]["PSequenceControls"])
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma, p_sequence_summary=sequence,
                                reference_order=5, gates_sha256=digest, interfaces=["MA", "MS"], gated_order=5)
        self.assertEqual(record["Verdict"], gates.VERDICT_PASSED)
        self.assertEqual(record["Reason"], "every gate passed (p_SA not applicable)")
        self.assertEqual(record["NotApplicable"], ["p_SA"])
        self.assertTrue(record["Gates"]["p_SA"]["NotApplicable"])
        self.assertIn("no SA interface", record["Gates"]["p_SA"]["Reason"])
        self.assertNotIn("NotApplicable", record["Gates"]["p_MA"])
        self.assertEqual(record["Gates"]["PSequenceControls"]["NotApplicableObservables"], ["p_SA"])
        self.assertEqual(set(record["Gates"]["PSequenceControls"]["Controls"]["1"]), {"E", "p_MA", "p_MS"})
        self.assertEqual(record["GatedOrder"], 5)
        self.assertEqual(record["ReferenceAnchor"], "vs p5 anchor")
        # A declared interface that fails still fails.
        comparison["PerSource"]["1"]["p_MA_rel"] = 0.021
        record = gates.evaluate(table, comparison=comparison, classes=classes, ref_pma=ref_pma, p_sequence_summary=sequence,
                                reference_order=5, gates_sha256=digest, interfaces=["MA", "MS"])
        self.assertEqual(record["Reason"], "failing gates ['p_MA'] (p_SA not applicable)")
        self.assertEqual(reference_campaign.interfaces(
            {"Boundaries": {"Postprocessing": {"Dielectric": [{"Type": "MS"}, {"Type": "MA"}]}}}), ["MA", "MS"])
        self.assertEqual(reference_campaign.interfaces({"Boundaries": {}}), [])
        self.assertEqual(gates.applicable_observables(None), list(p_sequence.GATED_OBSERVABLES))

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


class RemoteCommandTest(unittest.TestCase):
    def test_upload_commands_are_portable_rsync(self):
        """The macOS openrsync client has no --mkpath: directories are created by ssh, the
        rsync commands carry only -a; a single file (the runner) is copied as a file."""
        self.assertEqual(remote.upload_command("h", "/l/case", "/r/run/case"), ["rsync", "-a", "/l/case/", "h:/r/run/case/"])
        self.assertEqual(remote.upload_file_command("h", "/l/run_stages.py", "/r/run/run_stages.py"),
                         ["rsync", "-a", "/l/run_stages.py", "h:/r/run/run_stages.py"])
        self.assertEqual(remote.fetch_command("h", "/r/run/case/main", "/l/results/main"),
                         ["rsync", "-a", "--exclude", "archive/", "--exclude", "tmp/", "h:/r/run/case/main/", "/l/results/main/"])


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
        self.assertEqual(plan["Stages"][1]["Environment"]["PALACE_RESPONSE_BLOCK_SIZE"], "48")
        self.assertEqual(plan["ReducerBlockSize"], 48)
        self.assertEqual(build_plan.DEFAULT_REDUCER_BLOCK_SIZE, 48)
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
        # On the physics-11 model (80 sources measured at block size 6): its own stage reproduces itself.
        model = estimate_stages.load_cost_model(estimate_stages.PREVIOUS_COST_MODEL)
        counts = model["MeasuredMesh"]["EntityCounts"]
        same = estimate_stages.estimate_stage(model, 4, 80, counts, model["MeasuredBlockSize"])
        self.assertAlmostEqual(same["DOFRatioVsMeasured"], 1.0)
        self.assertAlmostEqual(same["ReducerSecondsEstimate"], model["Stages"]["p4"]["ReducerPalaceSeconds"], places=6)
        self.assertEqual(same["ReducerBlockPairs"], 105)
        self.assertEqual(same["ReducerSourceEvaluations"], 80 * 14)
        self.assertAlmostEqual(same["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"],
                               model["Stages"]["p4"]["WorkerNonSourceSeconds"] + 80 * model["Stages"]["p4"]["MeanPerSourceSeconds"], places=6)
        more = estimate_stages.estimate_stage(model, 4, 160, counts, model["MeasuredBlockSize"])
        self.assertGreater(more["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"], same["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"])
        self.assertGreater(more["ReducerPairs"], same["ReducerPairs"])

    def test_run_stages_parses_the_streaming_reducer_progress(self):
        """Decision 62(4): the streaming reducer logs one line per source and the per-rank
        sample counts; run_stages.parse_log records them beside the block-pair progress."""
        import importlib.util
        spec = importlib.util.spec_from_file_location("run_stages", HERE / "qualify" / "run_stages.py")
        run_stages = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(run_stages)
        text = (" Archived response reduction: 8 sources, 3 interfaces, quadrature samples per rank min 64665, "
                "max 100476, total 487776 (streaming; block size 48 recorded)\n"
                " Archived response source 1/8\n Archived response source 8/8\n"
                " Archived response reduction complete: 8 sources streamed once\n"
                "PCG solver converged in 12 iterations\n")
        parsed = run_stages.parse_log(text)
        self.assertEqual(parsed["StreamedSourcesProgress"], [[1, 8], [8, 8]])
        self.assertEqual(parsed["ReductionSamples"], {"Sources": 8, "Interfaces": 3, "PerRankMin": 64665,
                                                      "PerRankMax": 100476, "Total": 487776})
        self.assertEqual(parsed["BlockPairsProgress"], [])
        self.assertEqual(parsed["PCG"], [12])

    def test_reducer_estimate_shrinks_with_the_block_size_evaluation_part_only(self):
        """Decision 62(1): at block size b the reducer evaluates every source ceil(N / b)
        times; the evaluation fraction of the measured pair seconds scales with it, the
        Gram part with the source pairs, the setup not at all (on the physics-11 model, the
        one measured at block size 6)."""
        model = estimate_stages.load_cost_model(estimate_stages.PREVIOUS_COST_MODEL)
        counts = model["MeasuredMesh"]["EntityCounts"]
        measured = model["Stages"]["p4"]
        at6 = estimate_stages.estimate_stage(model, 4, 80, counts, 6)
        at48 = estimate_stages.estimate_stage(model, 4, 80, counts, 48)
        default = estimate_stages.estimate_stage(model, 4, 80, counts)
        self.assertEqual(default["ReducerBlockSize"], build_plan.DEFAULT_REDUCER_BLOCK_SIZE)
        self.assertEqual((at48["ReducerBlockPairs"], at48["ReducerSourceEvaluations"]), (3, 160))
        parts6, parts48 = at6["ReducerSecondsEstimateParts"], at48["ReducerSecondsEstimateParts"]
        self.assertAlmostEqual(parts6["Setup"], parts48["Setup"])
        self.assertAlmostEqual(parts6["Gram"], parts48["Gram"])
        self.assertAlmostEqual(parts48["Evaluation"], parts6["Evaluation"] * 160 / 1120)
        self.assertAlmostEqual(parts6["Evaluation"], measured["ReducerPairSeconds"] * model["ReducerEvaluationFraction"])
        self.assertAlmostEqual(sum(parts48.values()), at48["ReducerSecondsEstimate"])
        self.assertLess(at48["ReducerSecondsEstimate"], at6["ReducerSecondsEstimate"])
        self.assertGreater(at48["ReducerSecondsEstimate"], parts6["Setup"] + parts6["Gram"])
        # The calibration (PBS 46685, two-edge p4 at b = 6 / 48): evaluation fraction 0.97,
        # 0.068 GB of reducer peak per resident field per million H1 DOFs (84 more fields at b = 48).
        self.assertEqual((model["ReducerEvaluationFraction"], model["ReducerResidentFieldGBPerMillionH1"]), (0.97, 0.068))
        self.assertAlmostEqual(at6["ReducerPalacePeakGBEstimate"], measured["ReducerPalacePeakGB"])
        self.assertAlmostEqual(at48["ReducerPalacePeakGBEstimate"] - at6["ReducerPalacePeakGBEstimate"],
                               84 * 0.068 * measured["H1"] / 1e6)
        self.assertAlmostEqual(at48["ReducerResidentFieldsGBEstimate"], 84 * 0.068 * measured["H1"] / 1e6)
        # Reproduces the measured two-edge p4 pair (78 sources, H1 7.97M): 86.9 -> 15.3 s, 45.6 -> 91.1 GB.
        two_edge = {"Vertices": 91_000, "Edges": 600_000, "TriangleFaces": 1_000_000, "QuadFaces": 40_000,
                    "Tetrahedra": 470_000, "Prisms": 44_000, "Pyramids": 3_000}
        e6, e48 = (estimate_stages.estimate_stage(model, 4, 78, two_edge, b) for b in (6, 48))
        self.assertAlmostEqual(e6["ReducerSecondsEstimateParts"]["Evaluation"] / e48["ReducerSecondsEstimateParts"]["Evaluation"],
                               1014 / 156)
        whole = estimate_stages.estimate(counts, [("p4-80", 4, 80)], block_size=48, model=model)
        self.assertEqual(whole["ReducerBlockSize"], 48)
        self.assertEqual(whole["CostModel"]["MeasuredBlockSize"], 6)
        with self.assertRaises(ValueError):
            estimate_stages.estimate_stage(model, 4, 80, counts, 0)
        with self.assertRaises(ValueError):
            build_plan.reducer_environment(0)
        self.assertEqual(build_plan.reducer_environment(7), {"PALACE_RESPONSE_REDUCE_ONLY": "1", "PALACE_RESPONSE_BLOCK_SIZE": "7"})
        # summarize_cost reads the block size the reducer stage ran at.
        if not campaign_available("four-edge-physics-11"):
            self.skipTest("four-edge-physics-11 campaign not available")
        status = json.loads((StoredCampaign("four-edge-physics-11").results / "status.json").read_text())
        summary = summarize_cost.summarize(status, full_sources=80, nodes=1)
        main = summary["Stages"]["va-p4"]
        self.assertEqual((main["ReducerBlockSize"], main["ReducerBlockPairs"], main["ReducerSourceEvaluations"]), (6, 105, 1120))
        for stage in status["Stages"]:
            stage["Environment"].pop("PALACE_RESPONSE_BLOCK_SIZE", None)
        main = summarize_cost.summarize(status, full_sources=80, nodes=1, block_size=48)["Stages"]["va-p4"]
        self.assertEqual((main["ReducerBlockSize"], main["ReducerBlockPairs"]), (48, 3))

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
        comparison = compare_matrices.compare(tmp / "ref", tmp / "run", zero_trace_indices=[], interface_types=THREE_INTERFACES)
        self.assertEqual(comparison["ZeroTrace"], [3])     # zero diagonal energy
        self.assertEqual(comparison["Free"], [1, 2])
        self.assertAlmostEqual(comparison["PerSource"][2]["E_rel"], 0.01)
        self.assertAlmostEqual(comparison["PerSource"][2]["p_MA_rel"], 0.0)   # both scale: participation unchanged
        self.assertAlmostEqual(comparison["PerSource"][1]["E_rel"], 0.0)
        summary = compare_matrices.write_comparison(tmp / "ref", tmp / "run", tmp / "out", zero_trace_indices=[1],
                                                    interface_types=THREE_INTERFACES)
        self.assertEqual(summary["ZeroTrace"], [1, 3])
        self.assertTrue((tmp / "out.md").is_file() and (tmp / "out.csv").is_file())
        self.assertTrue(math.isfinite(summary["PerSource"]["2"]["E_rel"]))
        # The interface index -> type map comes from the config, never from the indices: with
        # the reference's indices labeled MS / MA / SA (index 1 = MS) the participations follow
        # the labels; two interfaces of one type sum; an unlabeled index fails closed.
        relabeled = compare_matrices.compare(tmp / "ref", tmp / "run", zero_trace_indices=[],
                                             interface_types={1: "MS", 2: "MA", 3: "SA"})
        self.assertAlmostEqual(relabeled["PerSource"][1]["p_MS_ref"], comparison["PerSource"][1]["p_MA_ref"])
        two_ma = compare_matrices.compare(tmp / "ref", tmp / "run", zero_trace_indices=[],
                                          interface_types={1: "MA", 2: "MA", 3: "SA"})
        self.assertAlmostEqual(two_ma["PerSource"][1]["p_MA_ref"], 2 * comparison["PerSource"][1]["p_MA_ref"])
        self.assertNotIn("p_MS_ref", two_ma["PerSource"][1])
        with self.assertRaises(ValueError):
            compare_matrices.compare(tmp / "ref", tmp / "run", zero_trace_indices=[], interface_types={1: "MA", 2: "MS"})
        with self.assertRaises(ValueError):
            compare_matrices.compare(tmp / "ref", tmp / "run", zero_trace_indices=[])
        self.assertAlmostEqual(ma_ms_offsets.reference_p_ma(tmp / "ref", {1: "MA", 2: "MA", 3: "SA"})[1], 2 * 0.2 / 2.0)
        sequence = p_sequence.observables(tmp / "ref", {1: "MS", 2: "MA", 3: "SA"})
        self.assertAlmostEqual(sequence[1]["p_MS"], 0.2 / 2.0)
        self.assertAlmostEqual(p_sequence.observables(tmp / "ref", {1: "MA", 2: "MA", 3: "SA"})[1]["Q_MA"], 0.4)

    def test_downward_and_multiple_layers_stop_the_source_location(self):
        """A synthetic downward trace set (Nz = -1: the metal below the plane) and a two-layer
        signature are refused with a recorded ScopeGuard reason before any role is assigned;
        the same traces with an upward single-layer signature locate."""
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(tmp, True))
        traces = tmp / "traces"
        traces.mkdir()
        # Apexes on the x = -1 side of a box [-1, 1]^2 x [-2.05, 2.1]: levels bottom, trench
        # (metal below the plane), plane, metal-top (-0.1), top.
        levels = [-2.05, -0.1, 0.0, 0.05, 2.1]
        apexes = [(-1.0, 0.0, z) for z in levels] + [(1.0, 0.0, -2.05), (1.0, 0.0, 2.1)]
        for index, (x, y, z) in enumerate(apexes, start=1):
            lines = ["x,y,z,V,triangle"]
            for k, (px, py, pz) in enumerate(apexes):
                lines.append(f"{px},{py},{pz},{1.0 if k == index - 1 else 0.0},{k // 3 + 1}")
            (traces / f"basis-{index:04d}.csv").write_text("\n".join(lines) + "\n")
        (tmp / "plan-view-boundary.csv").write_text("Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
                                                    "1,1,1,0.0,0,Continuation,-1.0,-1.0\n1,2,1,0.0,0,Physical,-1.0,0.5\n"
                                                    "1,3,1,0.0,0,Physical,-0.5,0.5\n1,4,1,0.0,0,Continuation,-0.5,-1.0\n")
        header = "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm\n"
        (tmp / "downward.csv").write_text(header + "1,0,1,-0.5,0,0,1,0,0,0,1,0,-1,-1,1,0\n")
        (tmp / "two-layers.csv").write_text(header + "1,0,1,-0.5,0,0,1,0,0,0,1,0,1,-1,1,0\n2,1,2,-0.5,0,0.6,1,0,0,0,1,0,-1,-1,1,0\n")
        (tmp / "upward.csv").write_text(header + "1,0,1,-0.5,0,0,1,0,0,0,1,0,1,-1,1,0\n")
        for name, guard in (("downward.csv", "DownwardLayers"), ("two-layers.csv", "MultipleLayers")):
            with self.assertRaises(locate_sources.UnsupportedSourceGeometry) as stop:
                locate_sources.locate_directory(traces, tmp / "plan-view-boundary.csv", tmp / f"{name}.locations.csv",
                                                signature_path=tmp / name)
            self.assertEqual(stop.exception.guard, guard)
            self.assertIn(f"ScopeGuard[{guard}]", str(stop.exception))
            self.assertFalse((tmp / f"{name}.locations.csv").exists())
        rows, geometry = locate_sources.locate_directory(traces, tmp / "plan-view-boundary.csv", tmp / "up.csv",
                                                         signature_path=tmp / "upward.csv")
        self.assertEqual(len(rows), 7)
        self.assertEqual(geometry["layers"], [(0.0, 1)])
        self.assertEqual(locate_sources.signature_layers([{"Pz": "0", "Nz": ""}]), [(0.0, 1)])

    def test_pending_qualification_never_masks_a_failed_p_sequence(self):
        table, digest = gates.load_gates()
        passing = {"1": {name: {"seq": {"d_high": 0.001, "d_low": 0.002, "r": 0.5}} for name in p_sequence.GATED_OBSERVABLES}}
        record = gates.evaluate(table, comparison=None, p_sequence_summary=passing, gates_sha256=digest)
        self.assertEqual(record["Verdict"], gates.VERDICT_PENDING)
        failing = json.loads(json.dumps(passing))
        failing["1"]["E"]["seq"]["d_high"] = 0.02
        record = gates.evaluate(table, comparison=None, p_sequence_summary=failing, gates_sha256=digest)
        self.assertEqual(record["Verdict"], gates.VERDICT_FAILED)
        self.assertIn("PSequenceControls", record["Reason"])
        self.assertEqual(gates.evaluate(table, comparison=None, p_sequence_summary={}, gates_sha256=digest)["Verdict"],
                         gates.VERDICT_FAILED)

    def test_source_class_thresholds_are_the_frozen_tables(self):
        table, _ = gates.load_gates()
        self.assertEqual(classify_sources.thresholds_of_gates(table), (classify_sources.NARROW_WIDTH, classify_sources.JUNCTION_REACH))
        self.assertEqual(classify_sources.thresholds_of_gates(table), (0.1, 0.6))
        with self.assertRaises(ValueError):
            classify_sources.classify_all({}, [], thresholds=(0.2, 0.6))


if __name__ == "__main__":
    unittest.main()
