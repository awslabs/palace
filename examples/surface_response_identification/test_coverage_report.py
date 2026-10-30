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


def _load(path):
    with open(path) as f:
        return json.load(f)


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
        depth = [r for r in report["block_B"] if r["term"] == "DEPTH (ES / Maxwell)"][0]
        self.assertTrue(depth["fired"])
        self.assertFalse(any(r["term"] == "Maxwell DEPTH" for r in report["block_B"]))  # MINOR-4: the row fires on ES ladders too
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

    def test_gap_term_one_sided_in_the_thin_bias_direction(self):
        # decision 1052 (the PLAN's reading): g > 0 -> the uncorrected share misses the library's raise -> thin + library
        # LOW -> the band extends DOWNWARD only; g < 0 -> upward only; no raw partition (PROXY) -> two-sided
        report = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")
        band = report["block_C"]["MS"]
        g = report["block_A"]["energy"]["MS"]["g"]
        self.assertGreater(g, 0.0)
        self.assertAlmostEqual(band["gap_pt"], 100.0 * (1.6 / 18.0) * g, places=9)
        self.assertAlmostEqual(band["low_pt"], -3.0 - 1.0 - band["gap_pt"], places=9)
        self.assertAlmostEqual(band["high_pt"], -3.0 + 1.0, places=9)
        # g < 0: scale the fabricated energies below the raw ones
        post, manifest = write_fixture(os.path.join(self.tmp, "neg"))
        with open(os.path.join(post, "surface-response-model-energy.csv"), "w") as f:
            f.write("source,evaluation,model,patch count,patch weight,domain correction (J)," + ",".join(
                f"fabricated surface energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
            for m in (1, 2, 3):
                f.write(_wide_row([1.0, 0.0, float(m), 10.0, 1.0, 0.0] + [0.5 * FAB[m][i] for i in (1, 2, 3)]) + "\n")
        report = cr.build_report(post, manifest, self.bands, ma_reading="sharp")
        band = report["block_C"]["MS"]
        self.assertLess(report["block_A"]["energy"]["MS"]["g"], 0.0)
        self.assertAlmostEqual(band["low_pt"], -3.0 - 1.0, places=9)
        self.assertAlmostEqual(band["high_pt"], -3.0 + 1.0 + band["gap_pt"], places=9)
        post, manifest = write_fixture(os.path.join(self.tmp, "proxy2"), raw_partition=False)
        band = cr.build_report(post, manifest, self.bands, ma_reading="sharp")["block_C"]["MS"]
        self.assertAlmostEqual(band["high_pt"] - band["low_pt"], 2.0 + 2.0 * band["gap_pt"], places=9)

    def test_named_residuals_grown_and_p3_on_the_graded_type_only(self):
        report = cr.build_report(self.post, self.manifest, self.bands, reference_order=3, ma_reading="sharp")
        rows = {r["term"]: r for r in report["block_B"]}
        self.assertTrue(rows["grown-coupon (T5) excess"]["fired"])
        grown_ma = report["block_A"]["energy"]["MA"]["grown_share"]
        self.assertAlmostEqual(rows["grown-coupon (T5) excess"]["extensions"]["MA"]["upper_pt"], 4.30 * grown_ma, places=9)
        self.assertTrue(rows["p3-REFERENCE / u_p"]["fired"])
        # DESIGN 1.2 / 5.1: u_p widens the p3-graded Type (MA) only; SA / MS are graded against the thick ES best
        self.assertEqual(set(rows["p3-REFERENCE / u_p"]["extensions"]), {"MA"})
        plain = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")["block_C"]
        self.assertAlmostEqual(report["block_C"]["SA"]["low_pt"], plain["SA"]["low_pt"], places=9)
        self.assertAlmostEqual(report["block_C"]["MS"]["high_pt"], plain["MS"]["high_pt"], places=9)
        self.assertAlmostEqual(report["block_C"]["MA"]["low_pt"], plain["MA"]["low_pt"] - cr.P3_UP_PCT, places=9)
        self.assertAlmostEqual(report["block_C"]["MA"]["high_pt"], plain["MA"]["high_pt"] + cr.P3_UP_PCT, places=9)
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


class AnchorRuleAndFallbacksTest(unittest.TestCase):
    """Decision 1052: a manifest-absent catalog model (the kappa0.1 first-order curvature node) is credited to the
    requirement of the manifest-known model it shares a raw label with; the geometry-fallbacks CSV's own-energy column."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covrep-")
        self.post, self.manifest = write_fixture(self.tmp)
        palace = _load(os.path.join(self.post, "palace.json"))
        palace["SurfaceResponse"]["ModelCatalog"] += [
            {"Index": 4, "Name": "curved-edge-concave-kappa0.1", "Topology": "curved edge"},
            {"Index": 5, "Name": "curved-edge-convex-kappa0.1", "Topology": "curved edge"}]
        _dump(palace, os.path.join(self.post, "palace.json"))
        # m1+m4: the isolated edge (manifest Exact) shares a label with the kappa node m4; m5 never shares a label
        with open(os.path.join(self.post, "surface-response-model-raw-energy.csv"), "w") as f:
            f.write("source,evaluation,label," + ",".join(f"raw energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
            f.write("1,0,m1+m4,+5.0e-16,+10.0e-16,+1.0e-16\n1,0,m2,+1.0e-16,+2.0e-16,+0.2e-16\n1,0,m3,+2.0e-16,+4.0e-16,+0.4e-16\n")
            f.write("1,0,m5,+0.2e-16,+0.4e-16,+0.04e-16\n1,0,unclaimed,+0.8e-16,+1.6e-16,+0.16e-16\n")
        with open(os.path.join(self.post, "surface-response-model-raw-groups.csv"), "w") as f:
            f.write("label,model,patches,weight length (m)\nm1+m4,1,30,+1.0e-06\nm1+m4,4,30,+1.0e-06\nm2,2,3,+1.0e-06\nm3,3,4,+4.0e-07\nm5,5,2,+1.0e-08\n")
        manifest = _load(self.manifest)
        manifest["Requirements"][3] = {"Hash": "d" * 64, "Topology": "CurvedEdge", "Status": "Interpolated", "TotalEdgeLength": 20.0,
                                       "Signature": {}, "Interfaces": [{"Type": "SA"}, {"Type": "MS"}, {"Type": "MA"}],
                                       "SelectedModels": [{"Name": "isolated-edge@convex-kappa0.38-cubic"}]}
        manifest["Summary"]["Uncovered"] = {"ByType": {}, "Length": 0.0}
        _dump(manifest, self.manifest)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_anchor_rule_credits_the_anchor_requirement(self):
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        # m1+m4 (5e-16) -> IsolatedEdge / Exact in full (the topology fallback would have said CurvedEdge / Interpolated);
        # m5 has no anchor -> the topology fallback (CurvedEdge Interpolated), recorded
        self.assertAlmostEqual(e["classes"]["IsolatedEdge"], 5.0e-16, delta=1e-28)
        self.assertAlmostEqual(e["models"]["Exact"], 6.0e-16, delta=1e-28)
        self.assertAlmostEqual(e["classes"]["CurvedEdge"], 0.2e-16, delta=1e-28)
        self.assertAlmostEqual(e["mapped_models_share"], (2.5e-16 + 0.2e-16) / W[1], places=9)
        mapped = {m["model"]: m for m in report["block_A"]["anchor_mapped_models"]}
        self.assertEqual(mapped[4]["via"], "anchor")
        self.assertEqual(mapped[4]["anchor_model"], 1)
        self.assertEqual((mapped[4]["status"], mapped[4]["class"]), ("Exact", "IsolatedEdge"))
        self.assertEqual(mapped[5]["via"], "topology")
        self.assertEqual((mapped[5]["status"], mapped[5]["class"]), ("Interpolated", "CurvedEdge"))
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        self.assertIn("ANCHOR rule", cr.render_markdown(report))

    def test_anchor_rule_is_consistent_without_a_curved_edge_requirement(self):
        manifest = _load(self.manifest)
        manifest["Requirements"] = manifest["Requirements"][:3]
        manifest["Identification"]["Totals"]["AssignedLength"] = 480.0
        _dump(manifest, self.manifest)
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["classes"]["IsolatedEdge"], 5.0e-16 + 0.2e-16, delta=1e-28)  # m4 by anchor; m5 by topology
        mapped = {m["model"]: m for m in report["block_A"]["anchor_mapped_models"]}
        self.assertEqual(mapped[4]["via"], "anchor")
        self.assertEqual(mapped[5]["mapped_to"], "IsolatedEdge")  # no anchor, no CurvedEdge requirement: the next topology

    def _write_raw_rows(self, rows):
        with open(os.path.join(self.post, "surface-response-model-raw-energy.csv"), "w") as f:
            f.write("source,evaluation,label," + ",".join(f"raw energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
            f.writelines(rows)
        with open(os.path.join(self.post, "surface-response-model-raw-groups.csv"), "w") as f:
            f.write("label,model,patches,weight length (m)\n")
            for row in rows:
                label = row.split(",")[2]
                for m in label.split("+"):
                    if m.startswith("m"):
                        f.write(f"{label},{m[1:]},3,+1.0e-06\n")

    def test_anchor_tie_break_is_independent_of_the_raw_csv_row_order(self):
        """Decision 1053: the kappa node m4 shares a ZERO-energy label with the Interpolated corner m3 (visited first) and
        a 5e-16 label with the Exact isolated edge m1: the straight-like member anchors it whatever the row order."""
        rows = ["1,0,m3+m4,+0.0e+00,+0.0e+00,+0.0e+00\n", "1,0,m1+m4,+5.0e-16,+10.0e-16,+1.0e-16\n",
                "1,0,m2,+1.0e-16,+2.0e-16,+0.2e-16\n", "1,0,m3,+2.0e-16,+4.0e-16,+0.4e-16\n",
                "1,0,m5,+0.2e-16,+0.4e-16,+0.04e-16\n", "1,0,unclaimed,+0.8e-16,+1.6e-16,+0.16e-16\n"]
        reports = []
        for order in (rows, list(reversed(rows)), rows[1:] + rows[:1]):
            self._write_raw_rows(order)
            reports.append(cr.build_report(self.post, self.manifest, None, ma_reading="sharp"))
        for report in reports:
            mapped = {m["model"]: m for m in report["block_A"]["anchor_mapped_models"]}
            self.assertEqual((mapped[4]["via"], mapped[4]["anchor_model"], mapped[4]["class"]), ("anchor", 1, "IsolatedEdge"))
            self.assertAlmostEqual(report["block_A"]["energy"]["SA"]["classes"]["IsolatedEdge"], 5.0e-16, delta=1e-28)
            # the attribution (model -> anchor / status / class) is identical; the record list keeps one entry per label
            # visited, in CSV order (the v2 report form of record), so it is compared as a set
            attribution = {(m["model"], m["via"], m["anchor_model"], m["status"], m["class"]) for m in report["block_A"]["anchor_mapped_models"]}
            self.assertEqual(attribution, {(4, "anchor", 1, "Exact", "IsolatedEdge"), (5, "topology", None, "Interpolated", "CurvedEdge")})
            for kind in ("classes", "models"):
                for T in ("SA", "MS", "MA"):
                    self.assertEqual(sorted(report["block_A"]["energy"][T][kind]), sorted(reports[0]["block_A"]["energy"][T][kind]))
                    for name, value in report["block_A"]["energy"][T][kind].items():
                        # the attribution is identical; only the summation-order roundoff (~1e-16 relative) may differ
                        self.assertAlmostEqual(value, reports[0]["block_A"]["energy"][T][kind][name], delta=1e-28)

    def test_anchor_tie_break_straight_like_then_shared_energy_then_index(self):
        catalog = {1: {"Name": "edge-a", "Topology": "isolated edge"}, 2: {"Name": "gap-b", "Topology": "same-conductor gap"},
                   3: {"Name": "corner-c", "Topology": "concave corner"}, 7: {"Name": "kappa", "Topology": "curved edge"},
                   8: {"Name": "kappa-2", "Topology": "curved edge"}, 9: {"Name": "kappa-3", "Topology": "curved edge"}}
        known = {"edge-a": {"Status": "Exact", "Class": "IsolatedEdge", "Topology": "IsolatedEdge"},
                 "gap-b": {"Status": "Exact", "Class": "SameConductorGap", "Topology": "SameConductorGap"},
                 "corner-c": {"Status": "Interpolated", "Class": "ConcaveCorner", "Topology": "ConcaveCorner"}}
        labels = {"m3+m7": {"SA": 9.0e-16}, "m1+m7": {"SA": 1.0e-16}, "m2+m7": {"SA": 2.0e-16},  # straight-like beats the corner's energy; gap-b's 2e-16 beats edge-a's 1e-16
                  "m1+m8": {"SA": 1.0e-16}, "m2+m8": {"SA": 1.0e-16},  # equal shared energy -> the lowest index (1)
                  "m3+m9": {"SA": 0.0}}  # no straight-like member -> the corner
        members = {label: [int(m[1:]) for m in label.split("+")] for label in labels}
        partition = {"labels": labels, "members": members, "weights": {}}
        self.assertEqual(cr.anchor_models(partition, catalog, known), {7: 2, 8: 1, 9: 3})
        shuffled = {"labels": dict(reversed(list(labels.items()))), "members": dict(reversed(list(members.items()))), "weights": {}}
        self.assertEqual(cr.anchor_models(shuffled, catalog, known), {7: 2, 8: 1, 9: 3})
        self.assertTrue(cr.is_straight_like("IsolatedEdge") and cr.is_straight_like("same-conductor strip") and cr.is_straight_like("DifferentConductorGap"))
        self.assertFalse(cr.is_straight_like("ConcaveCorner") or cr.is_straight_like("spatial edge cluster") or cr.is_straight_like(None))

    def test_geometry_fallbacks_own_energy_column(self):
        with open(os.path.join(self.post, "surface-response-geometry-fallbacks.csv"), "w") as f:
            f.write("source,class,type,count,length (m),raw energy (J),within-R energy (J),share\n")
            f.write("1,OverlapExcluded,MS,2,7.9e-06,3.0e-18,1.8e-15,1.67e-03\n")
        report = cr.build_report(self.post, self.manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["MS"]
        self.assertAlmostEqual(e["raw"]["GeometryFallbacks:OverlapExcluded"], 3.0e-18, delta=1e-30)
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 1.6e-16 - 0.6e-16 - 3.0e-18, delta=1e-28)
        self.assertEqual(report["block_A"]["energy"]["SA"]["raw"]["GeometryFallbacks:OverlapExcluded"], 0.0)  # an MS-only row


MX_RAW = {1: 4.0e-16, 2: 10.0e-16}   # Maxwell fixture: interfaces 1 SA, 2 MS; the first (f, exc) block only
MX_OUT = {1: 0.4e-16, 2: 1.0e-16}
MX_LOOP = {1: 0.1e-16, 2: 0.3e-16}


def write_maxwell_surface_q(dirpath, scale=1.0, e_norm=1.0e-12):
    """surface-Q-corrected.csv (Maxwell long form: one row per interface per (f, exc)) + surface-Q-edge.csv."""
    header = ["f (GHz)", "exc", "interface", "E_norm raw (J)", "E_norm corrected (J)", "E_norm corrected fixed-flux (J)",
              "E_norm self-consistent (J)", "E_surf raw (J)", "p_surf raw", "Q_surf raw", "E_surf corrected (J)", "p_surf corrected",
              "Q_surf corrected", "E_surf corrected fixed-flux (J)", "p_surf corrected fixed-flux", "Q_surf corrected fixed-flux",
              "E_surf self-consistent (J)", "p_surf self-consistent", "Q_surf self-consistent", "trace closure spread"]
    with open(os.path.join(dirpath, "surface-Q-corrected.csv"), "w") as f:
        f.write(",".join(header) + "\n")
        for exc in (1, 2):
            for i in (1, 2):
                ft = MX_RAW[i] * (1.0 + 0.1 * scale) * (1.0 + 0.5 * (exc - 1))
                f.write(_wide_row([1.0, exc, i, e_norm, e_norm, e_norm, float("nan"), MX_RAW[i], 1.0, 1.0, ft, 1.0, 1.0, ft * 1.01, 1.0, 1.0,
                                   float("nan"), float("nan"), float("nan"), 0.1]) + "\n")
    with open(os.path.join(dirpath, "surface-Q-edge.csv"), "w") as f:
        f.write("f (GHz),exc,interface,R (m),E_out (J),p_out,E_ann (J),p_ann\n")
        for exc in (1, 2):
            for i in (1, 2):
                f.write(_wide_row([1.0, exc, i, 1.9e-6, MX_OUT[i], 0.0, 0.0, 0.0]) + "\n")


def write_maxwell_fixture(root, e_mag_over_e_elec=11.5, ladder=False):
    """A synthetic Maxwell (Driven) follower: one solve (MaxIts 0) unless `ladder`, an R-50 excited port, domain-E.csv
    with E_mag / E_elec, the loop-remainder row, the confidence + loop-fit CSVs; a manifest with one Exact isolated edge."""
    post = os.path.join(root, "postpro")
    os.makedirs(post, exist_ok=True)
    config = {"Problem": {"Type": "Driven"},
              "Boundaries": {"Postprocessing": {"Dielectric": [{"Index": 1, "Type": "SA", "EdgeDistances": [1.9]},
                                                               {"Index": 2, "Type": "MS", "EdgeDistances": [1.9]}]},
                             "LumpedPort": [{"Index": 1, "R": 50.0, "L": 0.0, "C": 0.0, "Excitation": True}]},
              "Solver": {"SurfaceResponseCorrection": {"TargetInterfaces": [1, 2]}}}
    _dump(config, os.path.join(post, "config_resolved.json"))
    _dump({"Problem": {"DegreesOfFreedom": 2000}, "SurfaceResponse": {"ModelCatalog": [
        {"Index": 1, "Name": "isolated-edge", "Topology": "isolated edge"}]}}, os.path.join(post, "palace.json"))
    write_maxwell_surface_q(post, scale=1.0)
    with open(os.path.join(post, "surface-response-model-energy.csv"), "w") as f:
        f.write("source,evaluation,model,patch count,patch weight,domain correction (J),fabricated surface energy[1] (J),"
                "fabricated surface energy[2] (J),f (GHz),exc\n")
        f.write(_wide_row([1.0, 0.0, 1.0, 10.0, 1.0, 0.0, 3.3e-16, 8.5e-16, 1.0, 1]) + "\n")
        f.write(_wide_row([1.0, 0.0, 1.0, 10.0, 1.0, 0.0, 9.9e-16, 2.5e-15, 1.0, 2]) + "\n")
    with open(os.path.join(post, "surface-response-model-raw-energy.csv"), "w") as f:
        f.write("source,evaluation,label,raw energy[1] (J),raw energy[2] (J),f (GHz),exc\n")
        f.write(_wide_row([1.0, 0.0, "m1", 3.2e-16, 8.2e-16, 1.0, 1]) + "\n" + _wide_row([1.0, 0.0, "unclaimed", 0.4e-16, 0.8e-16, 1.0, 1]) + "\n")
        f.write(_wide_row([1.0, 0.0, "m1", 9.6e-16, 2.46e-15, 1.0, 2]) + "\n" + _wide_row([1.0, 0.0, "unclaimed", 1.2e-16, 2.4e-16, 1.0, 2]) + "\n")
    with open(os.path.join(post, "surface-response-model-raw-groups.csv"), "w") as f:
        f.write("label,model,patches,weight length (m)\nm1,1,30,+3.0e-06\n")
    with open(os.path.join(post, "surface-response-loop-remainder-energy.csv"), "w") as f:
        f.write("source,evaluation,model,patches,weighted claim length (m),remainder energy[1] (J),remainder energy[2] (J),f (GHz),exc\n")
        f.write(_wide_row([1.0, 0.0, "Total", 3.0, 2.0e-6, MX_LOOP[1], MX_LOOP[2], 1.0, 1]) + "\n")
        f.write(_wide_row([1.0, 0.0, "Total", 3.0, 2.0e-6, 3.0 * MX_LOOP[1], 3.0 * MX_LOOP[2], 1.0, 2]) + "\n")
    with open(os.path.join(post, "domain-E.csv"), "w") as f:
        f.write("f (GHz),E_elec (J),E_mag (J),E_cap (J),E_ind (J)\n" + _wide_row([1.0, 1.0e-12, e_mag_over_e_elec * 1.0e-12, 0.0, 0.0]) + "\n")
    with open(os.path.join(post, "surface-response-confidence.csv"), "w") as f:
        f.write("f (GHz),exc,loop residual,response-weighted loop residual,max trace closure spread,"
                "response-weighted local trace closure spread,confidence pass\n" + _wide_row([1.0, 1, 0.01, 0.001, 0.2, 0.07, 0]) + "\n")
    with open(os.path.join(post, "surface-response-loop-fit.csv"), "w") as f:
        f.write("source,evaluation,f (GHz),exc,max loop share\n" + _wide_row([1.0, 0.0, 1.0, 1, 0.01]) + "\n")
    with open(os.path.join(post, "error-indicators.csv"), "w") as f:
        f.write("Norm,Minimum,Maximum,Mean\n+1.0e-01,0,0,0\n")
    if ladder:
        for n, sc in enumerate((0.4, 0.9, 0.98), start=1):
            d = os.path.join(post, f"iteration{n:02d}")
            os.makedirs(d, exist_ok=True)
            write_maxwell_surface_q(d, scale=sc)
            with open(os.path.join(d, "error-indicators.csv"), "w") as f:
                f.write("Norm,Minimum,Maximum,Mean\n" + _wide_row([1.0 / n, 0, 0, 0]) + "\n")
    manifest = {"Requirements": [
        {"Hash": "f" * 64, "Topology": "IsolatedEdge", "Status": "Exact", "TotalEdgeLength": 300.0, "Signature": {},
         "Interfaces": [{"Type": "SA"}, {"Type": "MS"}], "SelectedModels": [{"Name": "isolated-edge"}]}],
        "Summary": {"Uncovered": {"ByType": {}, "Length": 0.0}, "SpatialClusterPriority": {"OutsideLength": 0.0},
                    "DomainBoundary": {"RawPortionLength": 0.0}, "Counts": {"Exact": 1}},
        "Identification": {"Totals": {"AssignedLength": 300.0}, "Exclusions": [], "Diagnostics": {"GeometryFallbacks": {"Entries": []}}}}
    mpath = os.path.join(root, "surface-response-requirements.json")
    _dump(manifest, mpath)
    return post, mpath


def write_maxwell_bands(path):
    bands = {"version": "test", "fit_date": "2026-10-30", "dev_set_caveat": "test", "provenance": {},
             "families": {"Maxwell": {"class_census": {"IsolatedEdge": 5},
                                      "types": {t: {"model": "M0", "coefficients": {"beta0": -3.0}, "H_pt": 1.0, "widened": 1.2,
                                                    "median_g": -0.1, "n": 12, "loo_coverage": {"coverage": 0.917, "ci90": [0.7, 0.98]},
                                                    "candidates": {"M0": {"loo_coverage": 0.833}}}
                                                for t in ("SA", "MS", "MA", "C")}}}}
    _dump(bands, path)
    return path


class MaxwellFixtureTest(unittest.TestCase):
    """DESIGN section 6: INDUCTIVE from a Maxwell fixture; the loop-remainder row; the C class; FORM; --ladder-run."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covrep-mx-")
        self.post, self.manifest = write_maxwell_fixture(os.path.join(self.tmp, "follower"))
        self.bands = cr.load_bands(write_maxwell_bands(os.path.join(self.tmp, "bands.json")))

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_inductive_loop_remainder_and_c_class(self):
        report = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")
        self.assertEqual(report["family"], "Maxwell")
        flags = report["block_D"]["flags"]
        self.assertTrue(flags["inductive"])
        self.assertAlmostEqual(flags["E_mag_over_E_elec"], 11.5, places=9)
        self.assertEqual(flags["c_class"], "R-50")
        self.assertFalse(flags["form"])
        rows = {r["term"]: r for r in report["block_B"]}
        self.assertTrue(rows["Maxwell INDUCTIVE"]["fired"])
        self.assertTrue(rows["the C form terms"]["fired"])
        self.assertFalse(rows["Maxwell FORM"]["fired"])
        self.assertTrue(any("INDUCTIVE" in r for r in report["label"]["reasons_not_high"]))
        # the first (f, exc) block only: W_SA = 4.0 - 0.4; the loop remainder row = 0.1e-16; Unclaimed = 0.4 - 0.1
        e = report["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["W"], MX_RAW[1] - MX_OUT[1], delta=1e-28)
        self.assertAlmostEqual(e["raw"]["LoopRemainder"], MX_LOOP[1], delta=1e-28)
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.4e-16 - MX_LOOP[1], delta=1e-28)
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        self.assertAlmostEqual(e["g"], (3.3e-16 - 3.2e-16) / 3.2e-16, places=9)
        # a follower (one solve) has no ladder of its own
        self.assertEqual(report["block_D"]["ladder"]["types"]["SA"]["status"], "ladder: not in this run")
        self.assertEqual(report["block_C"]["C"]["c_class"], "R-50")
        self.assertFalse(report["block_C"]["C"]["form"])
        # MINOR-2: the widening is disclosed with the unwidened LOO coverage beside it
        self.assertAlmostEqual(report["block_C"]["SA"]["loo_coverage_unwidened"], 0.833)
        md = cr.render_markdown(report)
        self.assertIn("x 1.2000 WIDENED", md)
        self.assertIn("unwidened 0.833", md)

    def test_form_flag_labels_the_c_column(self):
        report = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp", port_form="r2a-pembed")
        self.assertTrue(report["block_D"]["flags"]["form"])
        rows = {r["term"]: r for r in report["block_B"]}
        self.assertTrue(rows["Maxwell FORM"]["fired"])
        self.assertTrue(report["block_C"]["C"]["form"])
        self.assertIn("the C column is FORM", cr.render_markdown(report))
        # the Type bands do not move with FORM
        plain = cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")
        self.assertEqual(report["block_C"]["SA"]["low_pt"], plain["block_C"]["SA"]["low_pt"])
        rc = cr.main([self.post, "--manifest", self.manifest, "--bands", self.bands["_path"], "--output", os.path.join(self.tmp, "out"),
                      "--ma-reading", "sharp", "--port-form", "r2a-pembed"])
        self.assertEqual(rc, 0)
        # the config's Length override is the other trigger (DESIGN 1.2)
        config = _load(os.path.join(self.post, "config_resolved.json"))
        config["Boundaries"]["LumpedPort"][0]["Length"] = 2.0
        _dump(config, os.path.join(self.post, "config_resolved.json"))
        self.assertTrue(cr.build_report(self.post, self.manifest, self.bands, ma_reading="sharp")["block_D"]["flags"]["form"])

    def test_ladder_run_borrows_the_amr_ladder(self):
        amr, _ = write_maxwell_fixture(os.path.join(self.tmp, "amr"), ladder=True)
        report = cr.build_report(self.post, self.manifest, self.bands, ladder_dir=amr, ma_reading="sharp")
        lad = report["block_D"]["ladder"]
        self.assertTrue(lad["borrowed"])
        self.assertEqual(lad["types"]["SA"]["solves"], 4)
        self.assertEqual(lad["types"]["SA"]["status"], "quotable")
        self.assertAlmostEqual(lad["types"]["SA"]["steps_pt"][-1], 100.0 * 0.002 / 1.1, places=6)
        self.assertTrue(report["block_D"]["flags"]["indicator_falls"])
        self.assertIn("ladder borrowed from", cr.render_markdown(report))

    def test_non_inductive_when_e_mag_small(self):
        post, manifest = write_maxwell_fixture(os.path.join(self.tmp, "quiet"), e_mag_over_e_elec=0.01)
        report = cr.build_report(post, manifest, self.bands, ma_reading="sharp")
        self.assertFalse(report["block_D"]["flags"]["inductive"])
        self.assertFalse({r["term"]: r for r in report["block_B"]}["Maxwell INDUCTIVE"]["fired"])

    def test_c_band_must_be_m0(self):
        bands = _load(self.bands["_path"])
        bands["families"]["Maxwell"]["types"]["C"]["model"] = "M1"
        path = os.path.join(self.tmp, "bands-m1.json")
        _dump(bands, path)
        with self.assertRaises(ValueError):
            cr.build_report(self.post, self.manifest, cr.load_bands(path), ma_reading="sharp")


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
