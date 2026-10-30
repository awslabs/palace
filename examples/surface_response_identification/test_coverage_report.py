#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for coverage_report.py (readiness STAGE D; DESIGN.md section 6): a synthetic ES postpro fixture
with KNOWN numbers (the ENERGY partition, the shared-label split, the RAW rows, the Unclaimed remainder, the
ladder / DEPTH rule, the label rules, the PROXY path, the preflight form, the missing-input exit code), plus a
record test on the sct002-S1 stretchfix-regrade run (skipped when the records are not present)."""

import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import coverage_report as cr  # noqa: E402

RECORD_S1 = ("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004/v3h/"
             "stretchfix-regrade/runs/j2/sct002-S1-thin-P4-sfr/postpro")
RECORD_S1_MANIFEST = ("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004/v3h/"
                      "stretchfix-regrade/runs/j2/sct002-S1-thin-preflight-sfr-j2/postpro/surface-response-requirements.json")

# The fixture: 3 interfaces (1 SA, 2 MS, 3 MA), 3 models (m1 Exact isolated edge, m2 Exact SEC-3 grown, m3
# Interpolated concave corner), the shared label m1+m2 split by weight length 3 : 1.
E_RAW = {1: 10.0e-16, 2: 20.0e-16, 3: 2.0e-16}
E_OUT = {1: 1.0e-16, 2: 2.0e-16, 3: 0.2e-16}
W = {i: E_RAW[i] - E_OUT[i] for i in E_RAW}
RAW_LABELS = {  # within-R raw energy per label per interface
    "m1+m2": {1: 6.2e-16, 2: 12.4e-16, 3: 1.24e-16},
    "m3": {1: 2.0e-16, 2: 4.0e-16, 3: 0.4e-16},
    "unclaimed": {1: 0.8e-16, 2: 1.6e-16, 3: 0.16e-16},
}
UNCOVERED = {1: 0.3e-16, 2: 0.6e-16, 3: 0.06e-16}  # a recorded RAW class inside the unclaimed integral
FAB = {1: {1: 6.5e-16, 2: 12.5e-16, 3: 1.5e-16}, 2: {1: 0.2e-16, 2: 0.3e-16, 3: 0.1e-16}, 3: {1: 1.9e-16, 2: 3.9e-16, 3: 0.5e-16}}


def _dump(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f)


def _wide_row(values):
    return ",".join(f"{v:+.12e}" if isinstance(v, float) else str(v) for v in values)


def write_surface_q(dirpath, scale=1.0, e_elec=1.0e-12):
    """surface-Q-corrected.csv (ES wide form) with ft = raw x (1 + 0.1 x scale) per interface."""
    header = ["i", "E_elec raw (J)", "E_elec postprocessed fixed-trace (J)", "E_elec postprocessed fixed-flux (J)",
              "E_elec corrected (J)", "max trace closure spread", "response-weighted local trace closure spread",
              "trace-closure response fraction above limit", "confidence pass"]
    values = [1.0, e_elec, e_elec, e_elec, float("nan"), 0.1, 0.1, 0.99, 0]
    for i in (1, 2, 3):
        header += [f"E_surf raw[{i}] (J)", f"p_surf raw[{i}]", f"Q_surf raw[{i}]",
                   f"E_surf postprocessed fixed-trace[{i}] (J)", f"p_surf postprocessed fixed-trace[{i}]",
                   f"Q_surf postprocessed fixed-trace[{i}]", f"E_surf postprocessed fixed-flux[{i}] (J)",
                   f"p_surf postprocessed fixed-flux[{i}]", f"Q_surf postprocessed fixed-flux[{i}]",
                   f"E_surf corrected[{i}] (J)", f"p_surf corrected[{i}]", f"Q_surf corrected[{i}]", f"trace closure spread[{i}]"]
        ft = E_RAW[i] * (1.0 + 0.1 * scale)
        values += [E_RAW[i], 1.0, 1.0, ft, 1.0, 1.0, ft * 1.01, 1.0, 1.0, float("nan"), float("nan"), float("nan"), 0.1]
    with open(os.path.join(dirpath, "surface-Q-corrected.csv"), "w") as f:
        f.write(",".join(header) + "\n" + _wide_row(values) + "\n")
    with open(os.path.join(dirpath, "surface-Q-edge.csv"), "w") as f:
        f.write("i,exc,interface,R (m),E_out (J),p_out,E_ann (J),p_ann\n")
        for i in (1, 2, 3):
            f.write(_wide_row([1.0, 0, i, 1.9e-6, E_OUT[i], 0.0, 0.0, 0.0]) + "\n")


def write_fixture(root, raw_partition=True, iterations=("quotable",)):
    """A synthetic ES postpro + manifest. iterations: 'quotable' (steps shrinking) or 'depth' (not shrinking)."""
    post = os.path.join(root, "postpro")
    os.makedirs(post, exist_ok=True)
    config = {"Problem": {"Type": "Electrostatic"},
              "Boundaries": {"Postprocessing": {"Dielectric": [
                  {"Index": 1, "Type": "SA", "EdgeDistances": [1.9], "LocalizeEdgeEnergy": False},
                  {"Index": 2, "Type": "MS", "EdgeDistances": [1.9], "LocalizeEdgeEnergy": False},
                  {"Index": 3, "Type": "MA", "EdgeDistances": [1.9], "LocalizeEdgeEnergy": False}]}},
              "Solver": {"Electrostatic": {"ResponseCorrection": {"TargetInterfaces": [1, 2, 3], "PatchEnergy": raw_partition}}}}
    _dump(config, os.path.join(post, "config_resolved.json"))
    palace = {"Problem": {"DegreesOfFreedom": 1000}, "SurfaceResponse": {"ModelCatalog": [
        {"Index": 1, "Name": "isolated-edge", "Topology": "isolated edge"},
        {"Index": 2, "Name": "spatialedgecluster_edgecount-3_a5e14f9da16f", "Topology": "spatial edge cluster"},
        {"Index": 3, "Name": "concave-corner-150deg@corner-angle157.5-cubic", "Topology": "concave corner"}],
        "Diagnostics": {}}}
    _dump(palace, os.path.join(post, "palace.json"))
    write_surface_q(post, scale=1.0)
    with open(os.path.join(post, "surface-response-model-energy.csv"), "w") as f:
        f.write("source,evaluation,model,patch count,patch weight,domain correction (J)," + ",".join(
            f"fabricated surface energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
        for m in (1, 2, 3):
            for ev in (0, 1):
                f.write(_wide_row([1.0, float(ev), float(m), 10.0, 1.0, 0.0] + [FAB[m][i] * (1.0 + 0.01 * ev) for i in (1, 2, 3)]) + "\n")
    if raw_partition:
        with open(os.path.join(post, "surface-response-model-raw-energy.csv"), "w") as f:
            f.write("source,evaluation,label," + ",".join(f"raw energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
            for label, by_i in RAW_LABELS.items():
                f.write(f"1,0,{label}," + ",".join(f"{by_i[i]:+.12e}" for i in (1, 2, 3)) + "\n")
        with open(os.path.join(post, "surface-response-model-raw-groups.csv"), "w") as f:
            f.write("label,model,patches,weight length (m)\nm1+m2,1,30,+3.0e-06\nm1+m2,2,3,+1.0e-06\nm3,3,4,+4.0e-07\n")
    with open(os.path.join(post, "surface-response-uncovered-energy.csv"), "w") as f:
        f.write("source,evaluation,type,portions,length (m),clipped by spatial support (m)," + ",".join(
            f"uncovered raw energy[{i}] (J),uncovered share[{i}]" for i in (1, 2, 3)) + "\n")
        f.write("1,0,CurvedEdge,2,+5.0e-06,+0.0," + ",".join(f"{UNCOVERED[i]:+.12e},0.0" for i in (1, 2, 3)) + "\n")
    with open(os.path.join(post, "error-indicators.csv"), "w") as f:
        f.write("Norm,Minimum,Maximum,Mean\n+1.0e-01,0,0,0\n")
    steps = {"quotable": (0.4, 0.9, 0.98), "depth": (0.4, 0.55, 0.7)}[iterations[0]]
    for n, sc in enumerate(steps, start=1):
        d = os.path.join(post, f"iteration{n:02d}")
        os.makedirs(d, exist_ok=True)
        write_surface_q(d, scale=sc)
    manifest = {"Requirements": [
        {"Hash": "f" * 64, "Topology": "IsolatedEdge", "Status": "Exact", "TotalEdgeLength": 300.0, "Signature": {},
         "Interfaces": [{"Type": "SA"}, {"Type": "MS"}, {"Type": "MA"}], "SelectedModels": [{"Name": "isolated-edge"}]},
        {"Hash": "a5e14f9da16f" + "0" * 52, "Topology": "SpatialEdgeCluster", "Status": "Exact", "TotalEdgeLength": 100.0,
         "Signature": {"EdgeCount": 3}, "Interfaces": [{"Type": "SA"}, {"Type": "MS"}, {"Type": "MA"}],
         "SelectedModels": [{"Name": "spatialedgecluster_edgecount-3_a5e14f9da16f"}]},
        {"Hash": "c" * 64, "Topology": "ConcaveCorner", "Status": "Interpolated", "TotalEdgeLength": 80.0, "Signature": {},
         "Interfaces": [{"Type": "SA"}, {"Type": "MS"}, {"Type": "MA"}],
         "SelectedModels": [{"Name": "concave-corner-150deg@corner-angle157.5-cubic"}]},
        {"Hash": "d" * 64, "Topology": "CurvedEdge", "Status": "Missing", "TotalEdgeLength": 20.0, "Signature": {},
         "Interfaces": [{"Type": "SA"}, {"Type": "MS"}, {"Type": "MA"}], "SelectedModels": []}],
        "Summary": {"Uncovered": {"ByType": {"CurvedEdge": {"Length": 20.0}}, "Length": 20.0},
                    "SpatialClusterPriority": {"OutsideLength": 0.0}, "DomainBoundary": {"RawPortionLength": 0.0},
                    "Counts": {"Exact": 2, "Interpolated": 1, "Missing": 1}},
        "Identification": {"Totals": {"AssignedLength": 500.0}, "Exclusions": [{"Class": "NonManifold", "Count": 4, "Length": 40.0}],
                           "Diagnostics": {"GeometryFallbacks": {"Entries": []}}}}
    mpath = os.path.join(root, "surface-response-requirements.json")
    _dump(manifest, mpath)
    return post, mpath


def write_bands(path, census=None):
    bands = {"version": "test", "fit_date": "2026-10-30", "dev_set_caveat": "test", "provenance": {},
             "families": {"ES": {"class_census": census if census is not None else
                                 {"IsolatedEdge": 5, "SpatialEdgeCluster-2-3": 5, "ConcaveCorner": 5},
                                 "types": {t: {"model": "M0", "coefficients": {"beta0": -3.0}, "H_pt": 1.0, "widened": 1.0,
                                               "median_g": -0.1, "n": 5, "loo_coverage": {"coverage": 0.9, "ci90": [0.7, 0.97]}}
                                           for t in ("SA", "MS", "MA", "C")}}}}
    _dump(bands, path)
    return path


class EnergyPartitionTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covrep-")
        self.post, self.manifest = write_fixture(self.tmp)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_within_r_and_partition_sum_to_W(self):
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        self.assertEqual(report["energy_mode"], "ENERGY")
        for t, i in (("SA", 1), ("MS", 2), ("MA", 3)):
            e = report["block_A"]["energy"][t]
            self.assertAlmostEqual(e["W"], W[i], delta=1e-28)
            self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)

    def test_shared_label_split_by_weight_length(self):
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        # m1+m2 = 6.2e-16 split 3 : 1 -> isolated edge 4.65e-16, SEC 1.55e-16 (both Exact); m3 Interpolated 2e-16
        self.assertAlmostEqual(e["models"]["Exact"], 6.2e-16, delta=1e-28)
        self.assertAlmostEqual(e["models"]["Interpolated"], 2.0e-16, delta=1e-28)
        self.assertAlmostEqual(e["classes"]["IsolatedEdge"], 4.65e-16, delta=1e-28)
        self.assertAlmostEqual(e["classes"]["SpatialEdgeCluster-2-3"], 1.55e-16, delta=1e-28)
        self.assertAlmostEqual(e["classes"]["ConcaveCorner"], 2.0e-16, delta=1e-28)
        self.assertAlmostEqual(e["grown_share"], 1.55e-16 / W[1], places=9)

    def test_raw_rows_and_unclaimed_remainder(self):
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["MS"]
        self.assertAlmostEqual(e["raw"]["Missing"], UNCOVERED[2], delta=1e-28)
        # Unclaimed = the unclaimed label minus the recorded RAW rows (never hidden): 1.6e-16 - 0.6e-16
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 1.6e-16 - 0.6e-16, delta=1e-28)
        self.assertAlmostEqual(e["raw_share_total"], 1.6e-16 / W[2], places=9)
        # g = (sum fab - sum raw of the models) / sum raw: (12.5 + 0.3 + 3.9 - 16.4) / 16.4
        self.assertAlmostEqual(e["g"], (16.7e-16 - 16.4e-16) / 16.4e-16, places=9)

    def test_length_partition(self):
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        L = report["block_A"]["length"]["SA"]
        self.assertAlmostEqual(L["total"], 500.0)
        self.assertAlmostEqual(L["models"]["Exact"], 400.0)
        self.assertAlmostEqual(L["models"]["Interpolated"], 80.0)
        self.assertAlmostEqual(L["raw"]["Missing"], 20.0)
        self.assertAlmostEqual(L["missing_sub"]["buildable"], 20.0)
        self.assertAlmostEqual(report["block_A"]["length_bookkeeping_miss_um"], 0.0)

    def test_ladder_quotable(self):
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        lad = report["block_D"]["ladder"]["types"]["SA"]
        self.assertEqual(lad["solves"], 4)
        # ft = raw x (1 + 0.1 s), s = 0.4, 0.9, 0.98, 1.0 -> the last step 0.002 / 1.1 = 0.18 pt <= 0.3: quotable
        self.assertAlmostEqual(lad["steps_pt"][-1], 100.0 * 0.002 / 1.1, places=6)
        self.assertEqual(lad["status"], "quotable")
        self.assertAlmostEqual(lad["step_ratio"], 0.25, places=9)
        self.assertEqual(lad["tail_pt"], 0.0)

    def test_depth_rule_unbounded_and_tail_sign(self):
        shutil.rmtree(self.tmp)
        self.post, self.manifest = write_fixture(self.tmp, iterations=("depth",))
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        lad = report["block_D"]["ladder"]["types"]["SA"]
        # steps 0.15 / 0.15 / 0.30 (x 10 / 1.1 pt): ratio 2.0 >= 0.7 -> DEPTH-UNBOUNDED, tail = 3 x the last step, positive
        self.assertEqual(lad["status"], "DEPTH-UNBOUNDED")
        self.assertAlmostEqual(lad["tail_pt"], 3.0 * lad["last_step_pt"], places=9)
        self.assertGreater(lad["tail_pt"], 0.0)
        depth = [r for r in report["block_B"] if r["term"] == "Maxwell DEPTH"][0]
        self.assertTrue(depth["fired"])
        bands = cr.load_bands(write_bands(os.path.join(self.tmp, "bands.json")))
        report = cr.build_report(self.post, self.manifest, bands, ma_reading="sharp")
        band = report["block_C"]["SA"]
        # the band is extended DOWNWARD by the tail when the ft still rises
        self.assertLess(band["low_pt"], -3.0 - 1.0 - lad["tail_pt"] + 1e-9)
        self.assertIn("LOW", report["label"]["label"])

    def test_depth_rule_converging(self):
        solves = [{"ft": {"SA": v}} for v in (0.90, 0.96, 0.995, 0.997)]
        rule = cr.depth_rule(solves)["SA"]
        self.assertEqual(rule["status"], "quotable")
        self.assertAlmostEqual(rule["last_step_pt"], 100.0 * 0.002 / 0.997, places=9)
        self.assertEqual(rule["tail_pt"], 0.0)
        solves = [{"ft": {"SA": v}} for v in (0.90, 0.95, 0.975, 0.9875)]
        rule = cr.depth_rule(solves)["SA"]
        self.assertEqual(rule["status"], "DEPTH")  # last step 1.27 pt > 0.3, ratio 0.5 < 0.7: the Richardson tail
        self.assertAlmostEqual(rule["tail_pt"], rule["last_step_pt"] / 0.5, places=9)


class LabelAndFormsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covrep-")
        self.post, self.manifest = write_fixture(self.tmp)
        self.bands = cr.load_bands(write_bands(os.path.join(self.tmp, "bands.json")))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_label_medium_by_raw_share_and_low_by_uncalibrated_class(self):
        report = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")
        # the Unclaimed row (5.6 % of W_t) > 1 % and Exact 68.9 % < 95 % break HIGH; every class calibrated and
        # the RAW total 8.9 % <= 10 %: MEDIUM
        self.assertEqual(report["label"]["label"], "MEDIUM", report["label"])
        self.assertTrue(any("Unclaimed" in r for r in report["label"]["reasons_not_high"]))
        self.assertTrue(any("Exact share" in r for r in report["label"]["reasons_not_high"]))
        # an uncalibrated class above 5 % -> LOW (MINOR-7)
        bands = cr.load_bands(write_bands(os.path.join(self.tmp, "bands2.json"), census={"IsolatedEdge": 5}))
        report = cr.build_report(self.post, self.manifest, bands, ma_reading="sharp")
        self.assertEqual(report["label"]["label"], "LOW")
        self.assertTrue(any("UNCALIBRATED class ConcaveCorner" in r for r in report["label"]["reasons_low"]))

    def test_label_low_when_ma_read_raw(self):
        report = cr.build_report(self.post, self.manifest, self.bands, ma_reading="raw")
        self.assertEqual(report["label"]["label"], "LOW")
        self.assertTrue(any("MA read raw" in r for r in report["label"]["reasons_low"]))

    def test_label_high(self):
        # remove the unclaimed remainder and the uncovered row: everything Exact / Interpolated... Exact must be >= 95 %
        post, manifest = write_fixture(os.path.join(self.tmp, "high"))
        with open(os.path.join(post, "surface-response-model-raw-energy.csv"), "w") as f:
            f.write("source,evaluation,label,raw energy[1] (J),raw energy[2] (J),raw energy[3] (J)\n")
            f.write("1,0,m1+m2,+8.8e-16,+17.6e-16,+1.76e-16\n1,0,m3,+0.2e-16,+0.4e-16,+0.04e-16\n")
        os.remove(os.path.join(post, "surface-response-uncovered-energy.csv"))
        for d in cr.iteration_dirs(post):
            write_surface_q(d, scale={"iteration01": 0.9, "iteration02": 0.99, "iteration03": 0.998}[os.path.basename(d)])
        report = cr.build_report(post, manifest, self.bands, ma_reading="sharp")
        self.assertEqual(report["label"]["label"], "HIGH", report["label"])
        self.assertAlmostEqual(report["block_A"]["energy"]["SA"]["exact_share"], 8.8 / 9.0, places=9)

    def test_proxy_path_caps_at_medium(self):
        post, manifest = write_fixture(os.path.join(self.tmp, "proxy"), raw_partition=False)
        report = cr.build_report(post, manifest, self.bands, ma_reading="sharp")
        self.assertEqual(report["energy_mode"], "PROXY")
        self.assertNotEqual(report["label"]["label"], "HIGH")
        self.assertTrue(any("PROXY" in r for r in report["label"]["reasons_not_high"]))
        e = report["block_A"]["energy"]["SA"]
        # the model rows = W_t x the LENGTH share; the recorded Missing row keeps its own energy; Unclaimed = the remainder
        self.assertAlmostEqual(e["models"]["Exact"], W[1] * 400.0 / 500.0, delta=1e-28)
        self.assertAlmostEqual(e["raw"]["Missing"], UNCOVERED[1], delta=1e-28)
        self.assertAlmostEqual(sum(e["models"].values()) + sum(e["raw"].values()), W[1], delta=1e-28)
        band = report["block_C"]["SA"]
        self.assertEqual(band["g_source"], "the calibration median (no raw partition)")
        md = cr.render_markdown(report)
        self.assertIn("ENERGY = PROXY", md)

    def test_gap_term_one_sided(self):
        report = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")
        band = report["block_C"]["MS"]
        g = report["block_A"]["energy"]["MS"]["g"]
        self.assertGreater(g, 0.0)
        self.assertAlmostEqual(band["gap_pt"], 100.0 * (1.6 / 18.0) * g, places=9)
        self.assertAlmostEqual(band["low_pt"], -3.0 - 1.0, places=9)  # g > 0: the gap extends the upper edge only
        self.assertAlmostEqual(band["high_pt"], -3.0 + 1.0 + band["gap_pt"] + 0.0, places=9)

    def test_named_residuals_grown_and_p3(self):
        report = cr.build_report(self.post, self.manifest, self.bands, reference_order=3, ma_reading="sharp")
        rows = {r["term"]: r for r in report["block_B"]}
        self.assertTrue(rows["grown-coupon (T5) excess"]["fired"])
        grown_ma = report["block_A"]["energy"]["MA"]["grown_share"]
        self.assertAlmostEqual(rows["grown-coupon (T5) excess"]["extensions"]["MA"]["upper_pt"], 4.30 * grown_ma, places=9)
        self.assertTrue(rows["p3-REFERENCE / u_p"]["fired"])
        band = report["block_C"]["SA"]
        self.assertAlmostEqual(band["low_pt"], -3.0 - 1.0 - cr.P3_UP_PCT, places=9)
        self.assertFalse(rows["Maxwell INDUCTIVE"]["fired"])
        self.assertFalse(rows["the cut convention (reference side)"]["fired"])

    def test_preflight_form(self):
        report = cr.build_preflight_report(self.manifest)
        self.assertEqual(report["form"], "preflight")
        self.assertEqual(report["label"]["label"], "READY")
        self.assertAlmostEqual(report["ready_criterion_1"]["SA"]["uncorrected_share_by_length"], 20.0 / 500.0)
        self.assertIn("not available on a preflight", report["block_C"])
        self.assertIn("READY", cr.render_markdown(report))

    def test_missing_input_exit_code_2(self):
        os.remove(os.path.join(self.post, "surface-Q-corrected.csv"))
        rc = cr.main([self.post, "--manifest", self.manifest, "--no-bands", "--output", os.path.join(self.tmp, "out")])
        self.assertEqual(rc, 2)
        rc = cr.main([self.post, "--manifest", os.path.join(self.tmp, "absent.json"), "--no-bands"])
        self.assertEqual(rc, 2)

    def test_cli_writes_json_and_md(self):
        out = os.path.join(self.tmp, "out")
        rc = cr.main([self.post, "--manifest", self.manifest, "--bands", self.bands["_path"], "--output", out, "--ma-reading", "sharp"])
        self.assertEqual(rc, 0)
        with open(os.path.join(out, "coverage-report.json")) as f:
            report = json.load(f)
        self.assertEqual(report["bands"]["sha256"], self.bands["_sha256"])
        self.assertTrue(os.path.isfile(os.path.join(out, "coverage-report.md")))
        self.assertTrue(cr.summary_line(report).startswith("COVERAGE MEDIUM: SA exact"))


@unittest.skipUnless(os.path.isdir(RECORD_S1) and os.path.isfile(RECORD_S1_MANIFEST), "the records of record are not present")
class RecordS1Test(unittest.TestCase):
    def test_sct002_s1_partition(self):
        report = cr.build_report(RECORD_S1, RECORD_S1_MANIFEST, None, ma_reading="sharp")
        self.assertEqual(report["energy_mode"], "ENERGY")
        e = report["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["W"], 7.3893e-16, delta=1e-20)   # E_surf raw[1] - E_out[1] of the final solve
        self.assertLess(abs(e["bookkeeping_miss_rel"]), 1e-9)
        self.assertAlmostEqual(e["unclaimed_share"], 0.0531, delta=0.0005)  # OI-COVERAGE-EXCLUSION-ENERGY: every raw class 0
        self.assertEqual(report["block_D"]["ladder"]["types"]["SA"]["solves"], 9)
        self.assertEqual(report["block_D"]["ladder"]["types"]["SA"]["status"], "quotable")


if __name__ == "__main__":
    unittest.main()
