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


def write_surface_q(dirpath, scale=1.0, e_elec=1.0e-12, sources=(1,)):
    """surface-Q-corrected.csv (ES wide form) with ft = raw x (1 + 0.1 x scale) per interface; one row per source (source s
    carries the energies x s: the v2 per-source fixture)."""
    header = ["i", "E_elec raw (J)", "E_elec postprocessed fixed-trace (J)", "E_elec postprocessed fixed-flux (J)",
              "E_elec corrected (J)", "max trace closure spread", "response-weighted local trace closure spread",
              "trace-closure response fraction above limit", "confidence pass"]
    for i in (1, 2, 3):
        header += [f"E_surf raw[{i}] (J)", f"p_surf raw[{i}]", f"Q_surf raw[{i}]",
                   f"E_surf postprocessed fixed-trace[{i}] (J)", f"p_surf postprocessed fixed-trace[{i}]",
                   f"Q_surf postprocessed fixed-trace[{i}]", f"E_surf postprocessed fixed-flux[{i}] (J)",
                   f"p_surf postprocessed fixed-flux[{i}]", f"Q_surf postprocessed fixed-flux[{i}]",
                   f"E_surf corrected[{i}] (J)", f"p_surf corrected[{i}]", f"Q_surf corrected[{i}]", f"trace closure spread[{i}]"]
    rows = []
    for src in sources:
        values = [float(src), src * e_elec, src * e_elec, src * e_elec, float("nan"), 0.1, 0.1, 0.99, 0]
        for i in (1, 2, 3):
            ft = src * E_RAW[i] * (1.0 + 0.1 * scale)
            values += [src * E_RAW[i], 1.0, 1.0, ft, 1.0, 1.0, ft * 1.01, 1.0, 1.0, float("nan"), float("nan"), float("nan"), 0.1]
        rows.append(_wide_row(values))
    with open(os.path.join(dirpath, "surface-Q-corrected.csv"), "w") as f:
        f.write(",".join(header) + "\n" + "\n".join(rows) + "\n")
    with open(os.path.join(dirpath, "surface-Q-edge.csv"), "w") as f:
        f.write("i,exc,interface,R (m),E_out (J),p_out,E_ann (J),p_ann\n")
        for src in sources:
            for i in (1, 2, 3):
                f.write(_wide_row([float(src), 0, i, 1.9e-6, src * E_OUT[i], 0.0, 0.0, 0.0]) + "\n")


def write_fixture(root, raw_partition=True, iterations=("quotable",), total_row=False, scp_dropped=None, export_rows=None,
                  raw_labels=None, uncovered=None, scp_in_box=None, unpatched_length=None, sources=(1,)):
    """A synthetic ES postpro + manifest. iterations: 'quotable' (steps shrinking) or 'depth' (not shrinking).
    total_row: the binary's `Total` sum row appended to the uncovered CSV (the records of record carry it; decision 1125).
    scp_dropped: Identification.Diagnostics.SpatialClusterPriority.Patches entries (the F2 record; decision 1126 (2)).
    export_rows: surface-response-exclusion-energy.csv lines after the header (class,type,...; decision 1129).
    raw_labels / uncovered: overrides of RAW_LABELS / UNCOVERED.
    scp_in_box: {interface: J} -> the SCP clipped-correction CSV (one dropped patch row; raw whole = in-box + outside 1e-18,
    raw outside = 1e-18): the F2 credit's CAP (fix round 2A (i)); None = no CSV.
    unpatched_length: an UnpatchedPairs record of that length (um) on the manifest (the no-energy RAW class of 2A (iii)).
    sources: the ES sources written to every energy CSV; source s > 1 carries the raw labels and fabricated energies x s (the
    v2 per-source reading, decision 1209)."""
    raw_labels = RAW_LABELS if raw_labels is None else raw_labels
    uncovered = UNCOVERED if uncovered is None else uncovered
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
    write_surface_q(post, scale=1.0, sources=sources)
    with open(os.path.join(post, "surface-response-model-energy.csv"), "w") as f:
        f.write("source,evaluation,model,patch count,patch weight,domain correction (J)," + ",".join(
            f"fabricated surface energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
        for src in sources:
            for m in (1, 2, 3):
                for ev in (0, 1):
                    f.write(_wide_row([float(src), float(ev), float(m), 10.0, 1.0, 0.0] + [src * FAB[m][i] * (1.0 + 0.01 * ev) for i in (1, 2, 3)]) + "\n")
    if raw_partition:
        with open(os.path.join(post, "surface-response-model-raw-energy.csv"), "w") as f:
            f.write("source,evaluation,label," + ",".join(f"raw energy[{i}] (J)" for i in (1, 2, 3)) + "\n")
            for src in sources:
                for label, by_i in raw_labels.items():
                    f.write(f"{src},0,{label}," + ",".join(f"{src * by_i[i]:+.12e}" for i in (1, 2, 3)) + "\n")
        with open(os.path.join(post, "surface-response-model-raw-groups.csv"), "w") as f:
            f.write("label,model,patches,weight length (m)\nm1+m2,1,30,+3.0e-06\nm1+m2,2,3,+1.0e-06\nm3,3,4,+4.0e-07\n")
    if scp_in_box is not None:
        with open(os.path.join(post, cr.SCP_CLIPPED_CSV), "w") as f:
            f.write("# DERIVED (test fixture)\n")
            f.write("source,evaluation,patch,model,feature,weight," + ",".join(
                f"fab energy[{i}] (J),raw whole energy[{i}] (J),raw outside energy[{i}] (J),ratio[{i}],clipped correction[{i}] (J)" for i in (1, 2, 3)) + "\n")
            for src in sources:
                for ev in (0, 1):
                    f.write(f"{src},{ev},6,4,7,+1.0," + ",".join(
                        f"+0.0,{src * scp_in_box[i] + 1.0e-18:+.12e},+1.0e-18,+1.0,+0.0" for i in (1, 2, 3)) + "\n")
    with open(os.path.join(post, "surface-response-uncovered-energy.csv"), "w") as f:
        f.write("source,evaluation,type,portions,length (m),clipped by spatial support (m)," + ",".join(
            f"uncovered raw energy[{i}] (J),uncovered share[{i}]" for i in (1, 2, 3)) + "\n")
        f.write("1,0,CurvedEdge,2,+5.0e-06,+0.0," + ",".join(f"{uncovered[i]:+.12e},0.0" for i in (1, 2, 3)) + "\n")
        if total_row:
            f.write("1,0,Total,2,+5.0e-06,+0.0," + ",".join(f"{uncovered[i]:+.12e},0.0" for i in (1, 2, 3)) + "\n")
    if export_rows is not None:
        with open(os.path.join(post, "surface-response-exclusion-energy.csv"), "w") as f:
            f.write("source,class,type,class portions,class length (m),raw energy (J),interface raw energy (J),share\n")
            for row in export_rows:
                f.write("1," + row + "\n")
    with open(os.path.join(post, "error-indicators.csv"), "w") as f:
        f.write("Norm,Minimum,Maximum,Mean\n+1.0e-01,0,0,0\n")
    steps = {"quotable": (0.4, 0.9, 0.98), "depth": (0.4, 0.55, 0.7)}[iterations[0]]
    for n, sc in enumerate(steps, start=1):
        d = os.path.join(post, f"iteration{n:02d}")
        os.makedirs(d, exist_ok=True)
        write_surface_q(d, scale=sc, sources=sources)
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
    if unpatched_length is not None:
        manifest["Summary"]["UnpatchedPairs"] = {"Count": 1, "Length": unpatched_length}
        manifest["Identification"]["Diagnostics"]["UnpatchedPairs"] = {
            "Count": 1, "Length": unpatched_length, "Features": [{"Id": 7, "Type": "SameConductorGap", "Length": unpatched_length}]}
        manifest["Identification"]["Features"] = [{"Id": 7, "Hash": "e" * 64}]
    if scp_dropped is not None:
        manifest["Identification"]["Diagnostics"]["SpatialClusterPriority"] = {
            "Count": len(scp_dropped), "OutsideLength": sum(p.get("OutsideLength", 0.0) for p in scp_dropped), "Patches": scp_dropped,
            "WindowLength": sum(p["WindowLength"] for p in scp_dropped)}
        manifest["Summary"]["SpatialClusterPriority"]["OutsideLength"] = sum(p.get("OutsideLength", 0.0) for p in scp_dropped)
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
        # the record form (postoperator.cpp:2662-2769): the m-rows, then the binary's `Total` sum row
        f.write("source,evaluation,model,patches,weighted claim length (m),remainder energy[1] (J),remainder energy[2] (J),f (GHz),exc\n")
        f.write(_wide_row([1.0, 0.0, "m1", 3.0, 2.0e-6, MX_LOOP[1], MX_LOOP[2], 1.0, 1]) + "\n")
        f.write(_wide_row([1.0, 0.0, "Total", 3.0, 2.0e-6, MX_LOOP[1], MX_LOOP[2], 1.0, 1]) + "\n")
        f.write(_wide_row([1.0, 0.0, "m1", 3.0, 2.0e-6, 3.0 * MX_LOOP[1], 3.0 * MX_LOOP[2], 1.0, 2]) + "\n")
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
        # OI-COVERAGE-UNCLAIMED-INCELL (decisions 1124 / 1125 / 1126 (2)): S1's 5.31 % unclaimed remainder (8 SCP-dropped vertex
        # coupons, 32.29 um in-box, 4 owner clusters, all Exact) is CREDITED to the owners; the printed Unclaimed is 0
        credit = e["scp_owner_credit"]
        self.assertAlmostEqual(credit["share"], 0.0531, delta=0.0005)
        self.assertEqual((credit["dropped"], len(credit["owners"])), (8, 4))
        self.assertAlmostEqual(credit["in_box_length_um"], 32.287, delta=0.001)
        self.assertTrue(all(o["status"] == "Exact" for o in credit["owners"]))
        self.assertAlmostEqual(e["unclaimed_share"], 0.0, delta=1e-12)
        self.assertAlmostEqual(e["raw_share_total"], 0.0, delta=1e-12)   # every recorded raw class 0 on S1; no export
        self.assertEqual(report["block_D"]["ladder"]["types"]["SA"]["solves"], 9)
        self.assertEqual(report["block_D"]["ladder"]["types"]["SA"]["status"], "quotable")


def _strip_json_comments(text):
    """// and /* */ comments outside strings and trailing commas removed (the C++ PreprocessFile of the schema test)."""
    out = []
    i = 0
    in_string = False
    while i < len(text):
        c = text[i]
        if in_string:
            out.append(c)
            if c == "\\" and i + 1 < len(text):
                out.append(text[i + 1])
                i += 1
            elif c == '"':
                in_string = False
        elif c == '"':
            in_string = True
            out.append(c)
        elif text.startswith("//", i):
            j = text.find("\n", i)
            i = len(text) if j < 0 else j
            continue
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            i = len(text) if j < 0 else j + 2
            continue
        else:
            out.append(c)
        i += 1
    # the trailing commas Palace's preprocessor tolerates
    import re
    return re.sub(r",(\s*[}\]])", r"\1", "".join(out))


class TestNearMatchedStatusRows(unittest.TestCase):
    """DESIGN sections 4 / 5 (decisions 1075 Q1 / 1079 / 1096 (2)): NearMatched and GeometricNearMatched requirements
    are CORRECTED rows of their own in the LENGTH partition, never Missing and never folded into Exact."""

    def test_length_partition_rows(self):
        reqs = [
            {"Hash": "1" * 64, "Topology": "IsolatedEdge", "Status": "Exact", "TotalEdgeLength": 100.0, "Signature": {},
             "Interfaces": [{"Type": "SA"}], "SelectedModels": [{"Name": "isolated-edge"}]},
            {"Hash": "2" * 64, "Topology": "SameConductorStrip", "Status": "NearMatched", "TotalEdgeLength": 10.0,
             "Signature": {}, "Interfaces": [{"Type": "SA"}], "SelectedModels": [{"Name": "same-conductor-strip-2.91785um"}]},
            {"Hash": "3" * 64, "Topology": "ParallelEdgeCluster", "Status": "GeometricNearMatched", "TotalEdgeLength": 30.0,
             "Signature": {}, "Interfaces": [{"Type": "SA"}], "SelectedModels": [{"Name": "parallel-edge-cluster-af972e5a4a5d"}]},
            {"Hash": "4" * 64, "Topology": "CurvedEdge", "Status": "Missing", "TotalEdgeLength": 60.0, "Signature": {},
             "Interfaces": [{"Type": "SA"}], "SelectedModels": []}]
        manifest = {"Requirements": reqs, "Summary": {"Counts": {"Exact": 1, "NearMatched": 1, "GeometricNearMatched": 1, "Missing": 1}},
                    "Identification": {"Totals": {"AssignedLength": 200.0}, "Exclusions": [], "Diagnostics": {}}}
        tmp = tempfile.mkdtemp(prefix="cov-nearmatched-")
        try:
            mpath = os.path.join(tmp, "surface-response-requirements.json")
            _dump(manifest, mpath)
            loaded = cr.load_manifest(mpath)
            L = cr.length_partition(loaded, ["SA"])["types"]["SA"]
            self.assertAlmostEqual(L["models_share"]["Exact"], 0.5)
            self.assertAlmostEqual(L["models_share"]["NearMatched"], 0.05)
            self.assertAlmostEqual(L["models_share"]["GeometricNearMatched"], 0.15)
            self.assertAlmostEqual(L["raw_share"]["Missing"], 0.3)
            self.assertAlmostEqual(L["raw_share_total"], 0.3)
            self.assertEqual(loaded["models"]["same-conductor-strip-2.91785um"]["Status"], "NearMatched")
            self.assertIn("NearMatched", cr.MODEL_STATUSES)
            self.assertIn("GeometricNearMatched", cr.MODEL_STATUSES)
        finally:
            shutil.rmtree(tmp)


class TestFixRoundRecordRows(unittest.TestCase):
    """The identification fix round's record classes in the report (decision 1096 (4)): the RAW rows UnpatchedPairs (the
    feature's length leaves its requirement's corrected row; 1077 Q3) and UnclusterableComponent (joins the total), the NAMED
    rows (TaperCrossings / TaperMeanKeys / PlacedLocalDeviation / CornerArmOnArc / RemainderBelowBar), the exclusion-energy
    export's Excluded:* energy rows, the ConductorConsistency summary and the GeometricNearPairs join."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="cov-fixround-")
        reqs = [
            {"Hash": "a" * 64, "Topology": "IsolatedEdge", "Status": "Exact", "TotalEdgeLength": 100.0, "Signature": {},
             "Interfaces": [{"Type": "SA"}], "SelectedModels": [{"Name": "isolated-edge"}]},
            {"Hash": "b" * 64, "Topology": "SameConductorGap", "Status": "Exact", "TotalEdgeLength": 20.0, "Signature": {},
             "Interfaces": [{"Type": "SA"}], "SelectedModels": [{"Name": "same-conductor-gap-2um"}]}]
        manifest = {"Requirements": reqs, "Summary": {"Counts": {"Exact": 2}, "UnpatchedPairs": {"Count": 1, "Length": 8.0},
                                                      "TaperCrossings": {"Count": 3, "Length": 1.5, "Treatment": "ApproximatedByIsolatedEdge"},
                                                      "TaperMeanKeys": {"Count": 1, "Length": 20.0, "Treatment": "KeyedAtTheMean"},
                                                      "PlacedLocalDeviation": {"Count": 1, "Length": 6.0}, "CornerArmOnArc": {"Count": 2, "Length": 4.0}},
                    "Library": {"GeometricNearPairs": [{"Class": "ParallelEdgeCluster", "ModelA": "pec-A", "ModelB": "pec-B", "MaxOffsetDeviationOverR": 0.003}]},
                    "Identification": {"Totals": {"AssignedLength": 120.0},
                                       "Exclusions": [{"Class": "TruncationCut", "Count": 2, "Length": 30.0}, {"Class": "UnclusterableComponent", "Count": 5, "Length": 80.0}],
                                       "Features": [{"Id": 7, "Hash": "b" * 64, "Type": "SameConductorGap", "Length": 8.0, "Match": {"Status": "Unpatched"}}],
                                       "Diagnostics": {"UnpatchedPairs": {"Count": 1, "Length": 8.0, "Features": [{"Id": 7, "Type": "SameConductorGap", "Length": 8.0}]},
                                                       "RemainderBelowBar": {"Count": 2, "Length": 0.01}}}}
        self.mpath = os.path.join(self.tmp, "surface-response-requirements.json")
        _dump(manifest, self.mpath)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_length_rows(self):
        L = cr.length_partition(cr.load_manifest(self.mpath), ["SA"])
        t = L["types"]["SA"]
        self.assertAlmostEqual(t["total"], 200.0)                 # 120 assigned + 80 UnclusterableComponent
        self.assertAlmostEqual(t["models"]["Exact"], 112.0)        # 120 - the 8 um unpatched
        self.assertAlmostEqual(t["raw"]["UnpatchedPairs"], 8.0)
        self.assertAlmostEqual(t["raw"]["UnclusterableComponent"], 80.0)
        self.assertAlmostEqual(t["raw_share_total"], 88.0 / 200.0)
        self.assertEqual(sorted(t["named"]), ["CornerArmOnArc", "PlacedLocalDeviation", "RemainderBelowBar", "TaperCrossings", "TaperMeanKeys"])
        self.assertAlmostEqual(t["named"]["TaperMeanKeys"]["Share"], 0.1)
        self.assertEqual(t["named"]["TaperCrossings"]["Treatment"], "ApproximatedByIsolatedEdge")
        self.assertEqual(L["unpatched_pairs"][0]["Status"], "Exact")
        self.assertEqual(L["identification_raw_exclusions"], {"UnclusterableComponent": 80.0})
        self.assertEqual(len(L["geometric_near_pairs"]), 1)

    def test_exclusion_energy_export_rows(self):
        post = os.path.join(self.tmp, "postpro")
        os.makedirs(post)
        with open(os.path.join(post, "surface-response-exclusion-energy.csv"), "w") as f:
            f.write("source,class,type,class portions,class length (m),raw energy (J),interface raw energy (J),share\n")
            f.write("1,Claimed,SA,0,0,9e-17,1e-16,0.9\n1,NonManifold,SA,4,1e-5,6e-18,1e-16,0.06\n")
            f.write("1,Unclaimed-other,SA,2,2e-6,4e-18,1e-16,0.04\n1,Total,SA,6,1.2e-5,1e-17,1e-16,0.1\n")
            f.write("2,NonManifold,SA,4,1e-5,7e-18,1e-16,0.07\n")   # a second source block: not read
        rows = cr.exclusion_energy_rows(post)
        self.assertTrue(rows["present"])
        self.assertEqual(sorted(rows["classes"]), ["NonManifold", "Unclaimed-other"])
        self.assertAlmostEqual(rows["classes"]["NonManifold"]["SA"], 6e-18)
        self.assertAlmostEqual(rows["claimed"]["SA"], 9e-17)
        raw_rows, locations = cr.raw_class_energies(post, {1: {"Type": "SA"}}, [1], "ES")
        self.assertAlmostEqual(raw_rows["Excluded:NonManifold"]["SA"], 6e-18)
        self.assertEqual(locations["Excluded:NonManifold"][0]["portions"], 4.0)
        self.assertFalse(cr.exclusion_energy_rows(self.tmp)["present"])

    def test_conductor_summary_and_near_pairs_join(self):
        palace = {"SurfaceResponse": {"Diagnostics": {"ConductorConsistency": {
            "Count": 1, "TestedPatches": 5, "UnprobedModels": ["m-x"],
            "Records": [{"FloorApplied": True, "SelfConsistentPreGate": True, "Maxwell": True}, {"FloorApplied": False}]}}}}
        cc = cr.conductor_consistency_summary(palace)
        self.assertEqual((cc["excluded"], cc["tested"], cc["floor_applied"], cc["pre_gate_records"], cc["unprobed_models"]), (1, 5, 1, 1, ["m-x"]))
        palace["SurfaceResponse"]["Diagnostics"]["ConductorConsistency"]["FloorAppliedPatches"] = 3
        self.assertEqual(cr.conductor_consistency_summary(palace)["floor_applied"], 3)   # the leaf wins over the records
        self.assertFalse(cr.conductor_consistency_summary({})["present"])
        pairs = cr.near_pairs_rows([{"ModelA": "pec-A", "ModelB": "pec-B", "Class": "ParallelEdgeCluster"}],
                                   {"Pairs": [{"ModelA": "pec-B", "ModelB": "pec-A", "CorrectionDifference": {"SA": {"Max": 0.03}}, "MaxOverTypes": 0.03, "Flagged": True}]})
        self.assertTrue(pairs[0]["Flagged"])
        self.assertAlmostEqual(pairs[0]["MaxOverTypes"], 0.03)
        rows = cr.named_residuals("ES", {"types": {"SA": {}}}, {"types": {}}, {}, {}, None, True, {}, cc, pairs)
        by_key = {r["term"]: r for r in rows}
        self.assertTrue(by_key["SelfConsistentPreGate (the eigen conductor-consistency gate before the energies)"]["fired"])
        self.assertTrue(by_key["FloorAppliedPatches (the conductor-consistency floor)"]["fired"])
        self.assertTrue(by_key["UnprobedModels / ConductorConsistency exclusions"]["fired"])
        near = by_key["GeometricNearPairs (the library's in-band model pairs)"]
        self.assertTrue(near["fired"])
        self.assertIn("FLAGGED 1", near["detail"])


class TestCovfixRound(unittest.TestCase):
    """The pre-freeze coverage fixes (decisions 1125 / 1126 / 1128 / 1129): (1) the raw-class CSVs' Total row is not a class;
    (2) F2: the Unclaimed remainder credited to the SCP-dropped coupons' OWNER clusters by in-box length; (3) UnprobedModels =
    applied, in no RAW row; (4) the exclusion export's classes = a second partition beside the physical one (a C4-type record
    with the label 0 and an O4-type record with the label above the export) and the per-class gap term."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covfix-")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_total_row_of_a_raw_class_csv_is_not_a_class(self):
        post, manifest = write_fixture(self.tmp, total_row=True)
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        for t, i in (("SA", 1), ("MS", 2), ("MA", 3)):
            e = report["block_A"]["energy"][t]
            self.assertAlmostEqual(e["raw"]["Missing"], UNCOVERED[i], delta=1e-28)          # main: 2 x UNCOVERED
            self.assertAlmostEqual(e["raw"]["Unclaimed"], RAW_LABELS["unclaimed"][i] - UNCOVERED[i], delta=1e-28)
        self.assertEqual(report["block_A"]["energy_locations"]["Missing"], [{"type": "CurvedEdge", "portions": 2.0, "length_m": 5e-06}])

    def test_total_row_of_the_loop_remainder_csv_is_not_a_class(self):
        # review d1126 MAJOR-1: the Maxwell loop-remainder CSV's `model` column carries the m-rows and then the binary's `Total`
        # sum row (postoperator.cpp:2769); the LoopRemainder RAW row = the m-rows once (ctx003-C1 Maxwell read 2 x the Total)
        post, manifest = write_maxwell_fixture(os.path.join(self.tmp, "follower"))
        with open(os.path.join(post, "surface-response-loop-remainder-energy.csv"), "w") as f:
            f.write("source,evaluation,model,patches,weighted claim length (m),remainder energy[1] (J),remainder energy[2] (J),f (GHz),exc\n")
            for exc, scale in ((1, 1.0), (2, 3.0)):
                f.write(_wide_row([1.0, 0.0, "m12", 2.0, 1.5e-6, 0.6 * scale * MX_LOOP[1], 0.6 * scale * MX_LOOP[2], 1.0, exc]) + "\n")
                f.write(_wide_row([1.0, 0.0, "m21", 1.0, 0.5e-6, 0.4 * scale * MX_LOOP[1], 0.4 * scale * MX_LOOP[2], 1.0, exc]) + "\n")
                f.write(_wide_row([1.0, 0.0, "Total", 3.0, 2.0e-6, scale * MX_LOOP[1], scale * MX_LOOP[2], 1.0, exc]) + "\n")
                f.write(_wide_row([1.0, 1.0, "Total", 3.0, 2.0e-6, scale * MX_LOOP[1], scale * MX_LOOP[2], 1.0, exc]) + "\n")
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        unclaimed_label = {1: 0.4e-16, 2: 0.8e-16}        # the fixture's `unclaimed` rows of surface-response-model-raw-energy.csv
        for t, i in (("SA", 1), ("MS", 2)):
            e = report["block_A"]["energy"][t]
            self.assertAlmostEqual(e["raw"]["LoopRemainder"], MX_LOOP[i], delta=1e-30)       # head d8db64814e: 2 x MX_LOOP
            self.assertAlmostEqual(e["raw"]["Unclaimed"], unclaimed_label[i] - MX_LOOP[i], delta=1e-30)
            self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        # a clean window writes the Total row alone, with 0 (postoperator.cpp:2664): no LoopRemainder energy
        with open(os.path.join(post, "surface-response-loop-remainder-energy.csv"), "w") as f:
            f.write("source,evaluation,model,patches,weighted claim length (m),remainder energy[1] (J),remainder energy[2] (J),f (GHz),exc\n")
            f.write(_wide_row([1.0, 0.0, "Total", 0.0, 0.0, 0.0, 0.0, 1.0, 1]) + "\n")
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        self.assertEqual(e["raw"]["LoopRemainder"], 0.0)
        self.assertAlmostEqual(e["raw"]["Unclaimed"], unclaimed_label[1], delta=1e-30)
    SCP = [{"Patch": 6, "Feature": 7, "Model": "convex-corner-45deg", "Topology": "convex corner", "WindowLength": 3.8, "OutsideLength": 0.0,
            "Owner": {"Model": "spatialedgecluster_edgecount-3_a5e14f9da16f", "Patch": 0, "Feature": 1}},
           {"Patch": 9, "Feature": 8, "Model": "convex-corner-90deg", "Topology": "convex corner", "WindowLength": 3.8, "OutsideLength": 2.8,
            "Owner": {"Model": "concave-corner-150deg@corner-angle157.5-cubic", "Patch": 1, "Feature": 2}}]

    def test_f2_unclaimed_remainder_credited_to_the_owner_clusters_by_in_box_length(self):
        # v2 (fix round 2A (i)): the credit needs the dropped coupons' in-box energy (the clipped-correction CSV) >= the remainder
        post, manifest = write_fixture(self.tmp, total_row=True, scp_dropped=self.SCP, scp_in_box={1: 1.0e-16, 2: 2.0e-16, 3: 0.2e-16})
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        remainder = RAW_LABELS["unclaimed"][1] - UNCOVERED[1]          # 0.5e-16: the label minus the in-label Missing row
        credit = e["scp_owner_credit"]
        self.assertFalse(credit["capped"])
        self.assertAlmostEqual(credit["cap_J"], 1.0e-16, delta=1e-28)
        self.assertAlmostEqual(credit["energy"], remainder, delta=1e-28)
        self.assertAlmostEqual(credit["in_box_length_um"], 3.8 + 1.0)   # the second coupon's 2.8 um outside stay the SCP RAW row
        owners = {o["model"]: o for o in credit["owners"]}
        self.assertAlmostEqual(owners["spatialedgecluster_edgecount-3_a5e14f9da16f"]["energy"], remainder * 3.8 / 4.8, delta=1e-28)
        self.assertEqual((owners["spatialedgecluster_edgecount-3_a5e14f9da16f"]["status"], owners["spatialedgecluster_edgecount-3_a5e14f9da16f"]["class"]),
                         ("Exact", "SpatialEdgeCluster-2-3"))
        self.assertEqual(owners["concave-corner-150deg@corner-angle157.5-cubic"]["status"], "Interpolated")
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.0, delta=1e-28)
        self.assertAlmostEqual(e["models"]["Exact"], 6.2e-16 + remainder * 3.8 / 4.8, delta=1e-28)
        self.assertAlmostEqual(e["classes"]["ConcaveCorner"], 2.0e-16 + remainder * 1.0 / 4.8, delta=1e-28)
        self.assertAlmostEqual(e["raw_share_total"], UNCOVERED[1] / W[1], places=9)
        # g over the credited raw: fab 8.6e-16 vs raw 8.2e-16 + 0.5e-16
        self.assertAlmostEqual(e["g"], (8.6e-16 - 8.7e-16) / 8.7e-16, places=9)
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        self.assertEqual(len(report["block_A"]["scp_dropped"]), 2)
        self.assertIn("F2 OWNER CREDIT", cr.render_markdown(report))
        self.assertIn("not separately measured", cr.render_markdown(report))

    def test_scp_owner_credit_above_1_pct_caps_the_label_at_medium(self):
        # decision 1133 (ii): the F2 credit rests on an unmeasured premise; a credit > 1 % of W_t on any Type limits the label
        # to MEDIUM. A run HIGH on every other rule (Exact 95.9 % with the credit, no RAW row > 1 %, the ladder quotable)...
        bands = cr.load_bands(write_bands(os.path.join(self.tmp, "bands.json")))
        zero = {1: 0.0, 2: 0.0, 3: 0.0}
        labels = {"m1+m2": {1: 8.0e-16, 2: 16.0e-16, 3: 1.6e-16}, "m3": {1: 0.2e-16, 2: 0.4e-16, 3: 0.04e-16},
                  "unclaimed": {1: 0.8e-16, 2: 1.6e-16, 3: 0.16e-16}}
        in_box = {1: 1.0e-16, 2: 2.0e-16, 3: 0.2e-16}
        post, manifest = write_fixture(self.tmp, raw_labels=labels, uncovered=zero, scp_dropped=self.SCP, scp_in_box=in_box)
        for d in cr.iteration_dirs(post):
            write_surface_q(d, scale={"iteration01": 0.9, "iteration02": 0.99, "iteration03": 0.998}[os.path.basename(d)])
        report = cr.build_report(post, manifest, bands, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["scp_owner_credit"]["share"], 0.8e-16 / W[1], places=9)        # 8.9 % of W_t credited
        self.assertGreaterEqual(e["exact_share"], cr.HIGH_EXACT_SHARE)
        self.assertEqual(report["label"]["reasons_low"], [])
        self.assertEqual(report["label"]["label"], "MEDIUM", report["label"])                     # head 299adda955: HIGH
        self.assertEqual([r for r in report["label"]["reasons_not_high"] if "SCP owner credit" not in r], [])
        self.assertTrue(any(r.startswith("SA: SCP owner credit") and "> 1 % of W_t" in r and "premise unmeasured" in r
                            for r in report["label"]["reasons_not_high"]), report["label"])
        self.assertIn("SCP owner credit", cr.render_markdown(report))
        # ...and the same run with the credit at 0.56 % of W_t (the label sum unchanged) is HIGH
        shutil.rmtree(self.tmp)
        small = {"m1+m2": {1: 8.75e-16, 2: 17.5e-16, 3: 1.75e-16}, "m3": labels["m3"], "unclaimed": {1: 0.05e-16, 2: 0.1e-16, 3: 0.01e-16}}
        post, manifest = write_fixture(self.tmp, raw_labels=small, uncovered=zero, scp_dropped=self.SCP, scp_in_box=in_box)
        for d in cr.iteration_dirs(post):
            write_surface_q(d, scale={"iteration01": 0.9, "iteration02": 0.99, "iteration03": 0.998}[os.path.basename(d)])
        report = cr.build_report(post, manifest, bands, ma_reading="sharp")
        self.assertAlmostEqual(report["block_A"]["energy"]["SA"]["scp_owner_credit"]["share"], 0.05e-16 / W[1], places=9)
        self.assertEqual(report["label"]["label"], "HIGH", report["label"])

    def test_f2_no_credit_without_a_drop_or_with_a_negative_remainder(self):
        post, manifest = write_fixture(self.tmp)
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        self.assertEqual(e["scp_owner_credit"]["energy"], 0.0)
        self.assertAlmostEqual(e["raw"]["Unclaimed"], RAW_LABELS["unclaimed"][1] - UNCOVERED[1], delta=1e-28)
        shutil.rmtree(self.tmp)
        big = {1: 0.9e-16, 2: 1.8e-16, 3: 0.18e-16}      # the recorded Missing row above the label: a negative remainder
        post, manifest = write_fixture(self.tmp, scp_dropped=self.SCP, uncovered=big)
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        self.assertEqual(e["scp_owner_credit"]["energy"], 0.0)
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.8e-16 - 0.9e-16, delta=1e-28)

    def test_unprobed_models_is_applied_energy_not_a_raw_row(self):
        rows = ["Claimed,SA,0,0,9e-17,1e-16,0.9", "UnprobedModels,SA,3,2e-6,5e-18,1e-16,0.05", "CornerArmTrim,SA,2,1e-6,3e-18,1e-16,0.03",
                "Total,SA,5,3e-6,1e-16,1e-16,1.0"]
        post, manifest = write_fixture(self.tmp, export_rows=rows)
        export = cr.exclusion_energy_rows(post)
        self.assertEqual(sorted(export["classes"]), ["CornerArmTrim"])
        self.assertAlmostEqual(export["applied_unprobed"]["SA"], 5e-18)
        raw_rows, _ = cr.raw_class_energies(post, {1: {"Type": "SA"}}, [1], "ES")
        self.assertNotIn("Excluded:UnprobedModels", raw_rows)                     # main: a RAW row in the gap
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        self.assertAlmostEqual(report["block_A"]["energy"]["SA"]["unprobed_applied"], 5e-18)
        self.assertAlmostEqual(report["block_A"]["exclusion_energy_export"]["applied_unprobed_J"]["SA"], 5e-18)
        self.assertIn("UnprobedModels (the export's class of APPLIED patches", cr.render_markdown(report))

    def test_export_classes_are_a_second_partition_beside_the_physical_one(self):
        # a C4-type record: the 'unclaimed' label 0 and no Missing energy, yet CornerArmTrim > 0 in the export
        zero = {1: 0.0, 2: 0.0, 3: 0.0}
        labels = dict(RAW_LABELS, unclaimed=zero, m3={1: 2.8e-16, 2: 5.6e-16, 3: 0.56e-16})   # the labels still sum to W_t
        rows = ["Claimed,SA,0,0,9e-17,1e-16,0.9", "CornerArmTrim,SA,2,1e-6,4e-18,1e-16,0.04", "NonManifold,SA,1,1e-6,1e-18,1e-16,0.01",
                "Total,SA,3,2e-6,9.5e-17,1e-16,0.95"]
        post, manifest = write_fixture(self.tmp, raw_labels=labels, uncovered=zero, export_rows=rows)
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.0, delta=1e-28)                   # main: -5e-18 (negative)
        self.assertAlmostEqual(e["raw"]["Excluded:CornerArmTrim"], 4e-18, delta=1e-28)
        self.assertAlmostEqual(e["raw_share_total"], 5e-18 / W[1], places=9)
        self.assertAlmostEqual(e["raw_share_physical"], 0.0, places=12)
        self.assertEqual(e["excluded_shares"], {"CornerArmTrim": 4e-18 / W[1], "NonManifold": 1e-18 / W[1]})
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)     # the physical partition sums to W_t
        shutil.rmtree(self.tmp)
        # an O4-type record: the label above the export; the remainder is the label minus the in-label rows only
        post, manifest = write_fixture(self.tmp, export_rows=rows, total_row=True)
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.8e-16 - 0.3e-16, delta=1e-28)       # main: 0.5e-16 - 5e-18
        self.assertAlmostEqual(e["raw_share_total"], (0.8e-16 + 5e-18) / W[1], places=9)
        self.assertAlmostEqual(e["raw_share_physical"], 0.8e-16 / W[1], places=9)

    def test_gap_term_per_class_reduces_to_total_raw_without_the_export(self):
        bands = cr.load_bands(write_bands(os.path.join(self.tmp, "bands.json")))
        post, manifest = write_fixture(self.tmp, total_row=True)
        row = cr.build_report(post, manifest, bands, ma_reading="sharp")["block_C"]["SA"]
        self.assertEqual(row["gap_parts"]["classes"], {})
        self.assertAlmostEqual(row["gap_pt"], 100.0 * row["raw_share"] * abs(row["g"]), places=9)   # the TOTAL raw x |g| form
        self.assertAlmostEqual(row["gap_parts"]["r_other"], row["raw_share"], places=12)
        shutil.rmtree(self.tmp)
        rows = ["Claimed,SA,0,0,9e-17,1e-16,0.9", "CornerArmTrim,SA,2,1e-6,4e-18,1e-16,0.04", "TruncationCut,SA,1,1e-6,1e-18,1e-16,0.01",
                "UnprobedModels,SA,3,2e-6,5e-18,1e-16,0.05"]
        post, manifest = write_fixture(self.tmp, total_row=True, export_rows=rows)
        report = cr.build_report(post, manifest, bands, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        row = report["block_C"]["SA"]
        g_corner = (1.9e-16 - 2.0e-16) / 2.0e-16          # m3: fab 1.9 vs raw 2.0 (the only corner-class model)
        self.assertAlmostEqual(e["g_corner"], g_corner, places=9)
        parts = row["gap_parts"]["classes"]
        self.assertEqual(sorted(parts), ["CornerArmTrim", "TruncationCut"])                 # UnprobedModels in no part
        self.assertAlmostEqual(parts["CornerArmTrim"]["b_c"], abs(g_corner), places=9)
        self.assertAlmostEqual(parts["TruncationCut"]["b_c"], abs(row["g"]), places=9)
        expected = 100.0 * (e["raw_share_physical"] * abs(row["g"]) + 4e-18 / W[1] * abs(g_corner) + 1e-18 / W[1] * abs(row["g"]))
        self.assertAlmostEqual(row["gap_pt"], expected, places=9)
        self.assertAlmostEqual(row["gap_parts"]["r_other"], 0.8e-16 / W[1], places=9)        # the label-based raw share
        self.assertGreater(row["gap_pt"], 100.0 * (row["raw_share"] - 5e-18 / W[1]) * abs(row["g"]) - 1e-12)

class TestBandFileLocation(unittest.TestCase):
    """OI-SCHEMA-TEST-COVERAGE-BANDS (decision 1065): test-schema's "Example Configs" case parses EVERY *.json
    directly under examples/<dir>/ as a Palace config, so the committed band files live in data/ (the default
    path of the tool) and no non-config JSON may sit at the package level."""

    def test_default_band_file_is_in_data_and_loads(self):
        here = os.path.dirname(os.path.abspath(cr.__file__))
        self.assertEqual(os.path.dirname(cr.DEFAULT_BANDS), os.path.join(here, "data"))
        self.assertTrue(os.path.isfile(cr.DEFAULT_BANDS), cr.DEFAULT_BANDS)
        bands = cr.load_bands(cr.DEFAULT_BANDS)
        self.assertIn("families", bands)
        self.assertTrue(bands.get("version"))

    def test_no_json_at_the_package_level(self):
        here = os.path.dirname(os.path.abspath(cr.__file__))
        stray = sorted(f for f in os.listdir(here) if f.endswith(".json"))
        self.assertEqual(stray, [], "a non-config JSON at examples/<dir>/ level fails test-schema 'Example Configs'")

    def test_examples_level_json_are_palace_configs(self):
        """The test-schema walk emulated (test/unit/test-schema.cpp: every *.json directly under examples/<dir>/,
        comments stripped, parsed, validated): the config schema requires exactly the five top-level sections."""
        examples = os.path.dirname(os.path.dirname(os.path.abspath(cr.__file__)))
        schema_path = os.path.join(os.path.dirname(examples), "scripts", "schema", "config-schema.json")
        with open(schema_path) as f:
            schema = json.load(f)
        required = set(schema["required"])
        allowed = set(schema["properties"])
        found = 0
        for entry in sorted(os.listdir(examples)):
            sub = os.path.join(examples, entry)
            if not os.path.isdir(sub):
                continue
            for name in sorted(os.listdir(sub)):
                if not name.endswith(".json"):
                    continue
                found += 1
                with open(os.path.join(sub, name)) as f:
                    text = f.read()
                config = json.loads(_strip_json_comments(text))
                keys = set(config)
                self.assertTrue(required <= keys <= allowed, f"{sub}/{name}: top-level keys {sorted(keys)}")
        self.assertGreater(found, 0)


if __name__ == "__main__":
    unittest.main()


class TestFixRoundV2(unittest.TestCase):
    """Fix round P1 (decisions 1209 / 1230 / 1233 (2); DESIGN-v3 s2 2A + s6): the F2 credit CAP, the Unclaimed (unattributed) row,
    the no-energy RAW classes' ESTIMATE, the measured `unpatched` label, and the PER-SOURCE reading of a multi-source ES run."""
    SCP = TestCovfixRound.SCP

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covrep-")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_f2_credit_capped_at_the_in_box_energy(self):
        # the remainder 0.5e-16 (SA) exceeds the dropped coupons' in-box energy 0.2e-16: the credit is CAPPED there, the rest is
        # the RAW row Unclaimed (unattributed) and enters raw_share_total
        post, manifest = write_fixture(self.tmp, total_row=True, scp_dropped=self.SCP, scp_in_box={1: 0.2e-16, 2: 0.4e-16, 3: 0.04e-16})
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        remainder = RAW_LABELS["unclaimed"][1] - UNCOVERED[1]
        credit = e["scp_owner_credit"]
        self.assertTrue(credit["capped"])
        self.assertAlmostEqual(credit["uncapped_energy"], remainder, delta=1e-28)
        self.assertAlmostEqual(credit["energy"], 0.2e-16, delta=1e-28)
        self.assertAlmostEqual(credit["cap_J"], 0.2e-16, delta=1e-28)
        self.assertIn(cr.SCP_CLIPPED_CSV, credit["cap_source"])
        owners = {o["model"]: o for o in credit["owners"]}
        self.assertAlmostEqual(owners["spatialedgecluster_edgecount-3_a5e14f9da16f"]["energy"], 0.2e-16 * 3.8 / 4.8, delta=1e-28)
        self.assertAlmostEqual(e["unclaimed_unattributed"], remainder - 0.2e-16, delta=1e-28)
        self.assertAlmostEqual(e["raw"]["Unclaimed"], remainder - 0.2e-16, delta=1e-28)
        self.assertAlmostEqual(e["raw_share_total"], (UNCOVERED[1] + remainder - 0.2e-16) / W[1], places=9)
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        self.assertEqual(report["tool"], "coverage_report.py v2")
        md = cr.render_markdown(report)
        self.assertIn("RAW: Unclaimed (unattributed)", md)

    def test_f2_credit_without_the_clipped_csv_is_capped_at_the_scp_raw_row(self):
        # no clipped-correction CSV and no SCP RAW energy -> the cap is 0: nothing is credited (fail-closed); the whole remainder is
        # Unclaimed (unattributed)
        post, manifest = write_fixture(self.tmp, total_row=True, scp_dropped=self.SCP)
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        remainder = RAW_LABELS["unclaimed"][1] - UNCOVERED[1]
        self.assertTrue(e["scp_owner_credit"]["capped"])
        self.assertEqual(e["scp_owner_credit"]["energy"], 0.0)
        self.assertIn("SCP RAW row", e["scp_owner_credit"]["cap_source"])
        self.assertAlmostEqual(e["raw"]["Unclaimed"], remainder, delta=1e-28)
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)

    def test_unpatched_pairs_estimate_from_the_unattributed_remainder(self):
        # 2A (iii): UnpatchedPairs (a no-energy RAW class of the manifest) receives the unattributed remainder as an ESTIMATE
        # (the only no-energy class -> the whole remainder, an upper bound: W31's case); Unclaimed reads 0
        post, manifest = write_fixture(self.tmp, total_row=True, unpatched_length=8.0)
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        e = report["block_A"]["energy"]["SA"]
        remainder = RAW_LABELS["unclaimed"][1] - UNCOVERED[1]
        self.assertAlmostEqual(e["raw"]["UnpatchedPairs"], remainder, delta=1e-28)
        self.assertIn("UnpatchedPairs", e["raw_estimated"])
        self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.0, delta=1e-28)
        self.assertAlmostEqual(e["raw_share_total"], (UNCOVERED[1] + remainder) / W[1], places=9)
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        self.assertAlmostEqual(report["block_A"]["length"]["SA"]["raw"]["UnpatchedPairs"], 8.0)
        md = cr.render_markdown(report)
        self.assertIn("(ESTIMATE)", md)
        self.assertNotIn("(in Unclaimed)", md.split("### MS")[0])

    def test_measured_unpatched_label_supersedes_the_estimate(self):
        # 2B: a binary exporting the `unpatched` label -> the MEASURED UnpatchedPairs RAW row; the unclaimed label stays the remainder
        labels = dict(RAW_LABELS, unpatched={1: 0.1e-16, 2: 0.2e-16, 3: 0.02e-16})
        labels["m1+m2"] = {1: 6.1e-16, 2: 12.2e-16, 3: 1.22e-16}
        post, manifest = write_fixture(self.tmp, total_row=True, raw_labels=labels, unpatched_length=8.0)
        e = cr.build_report(post, manifest, None, ma_reading="sharp")["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e["raw"]["UnpatchedPairs"], 0.1e-16, delta=1e-28)
        self.assertNotIn("UnpatchedPairs", e["raw_estimated"])
        self.assertAlmostEqual(e["raw"]["Unclaimed"], RAW_LABELS["unclaimed"][1] - UNCOVERED[1], delta=1e-28)
        self.assertAlmostEqual(e["models"]["Exact"], 6.1e-16, delta=1e-28)
        self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)

    def test_multi_source_run_read_per_source(self):
        # OI-COVERAGE-READER-MULTI-SOURCE-LAST-WINS (1209): source 2 carries every energy x 2; v1 let its rows overwrite source 1's
        # models (an Exact share > 100 % on W10). v2: the report reads ONE source (the first by default; --source), carries every
        # source's blocks and the union.
        post, manifest = write_fixture(self.tmp, total_row=True, sources=(1, 2))
        report = cr.build_report(post, manifest, None, ma_reading="sharp")
        self.assertEqual((report["source"], report["sources_in_run"]), (1, [1, 2]))
        e1 = report["block_A"]["energy"]["SA"]
        self.assertAlmostEqual(e1["W"], W[1], delta=1e-28)
        self.assertAlmostEqual(e1["models"]["Exact"], 6.2e-16, delta=1e-28)             # head 1928a4f2: source 2's 12.4e-16 over W_t of source 1
        self.assertLess(abs(e1["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
        self.assertEqual(sorted(report["sources"]), ["1", "2"])
        s2 = report["sources"]["2"]["energy"]["SA"]
        self.assertAlmostEqual(s2["exact_share"], e1["exact_share"], places=12)
        u = report["source_union"]["types"]["SA"]
        self.assertEqual(u["n_sources"], 2)
        self.assertAlmostEqual(u["exact_share_range"][0], e1["exact_share"], places=12)
        second = cr.build_report(post, manifest, None, ma_reading="sharp", source=2, all_sources=False)
        self.assertEqual(second["source"], 2)
        self.assertAlmostEqual(second["block_A"]["energy"]["SA"]["W"], 2.0 * W[1], delta=1e-28)
        self.assertAlmostEqual(second["block_A"]["energy"]["SA"]["models"]["Exact"], 12.4e-16, delta=1e-28)
        self.assertNotIn("sources", second)
        # the ladder follows the source too
        self.assertEqual(report["block_D"]["ladder"]["types"]["SA"]["solves"], second["block_D"]["ladder"]["types"]["SA"]["solves"])


STAGE_E = "/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/readiness-20261029/stageE"
STAGE_E_MANIFESTS = os.path.join(STAGE_E, "thin", "manifests", "es-loadcheck-v3-k-2d-E")
STAGE_E_NEAR_PAIRS = os.path.join(STAGE_E, "lib", "library", "s2-r1p9-v3-k-2d-E", "geometric-near-pairs.json")


def _value_differences(old, new, ignore, path=""):
    """The leaves of `old` whose value differs in `new` (added keys of `new` are not differences; NaN == NaN)."""
    out = []
    if isinstance(old, dict) and isinstance(new, dict):
        for k, v in old.items():
            if k not in new:
                out.append(("removed", path + "/" + k))
            else:
                out += _value_differences(v, new[k], ignore, path + "/" + k)
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            out.append(("length", path))
        else:
            for i, (a, b) in enumerate(zip(old, new)):
                out += _value_differences(a, b, ignore, f"{path}[{i}]")
    elif old != new and path not in ignore and not (isinstance(old, float) and isinstance(new, float) and old != old and new != new):
        out.append(("value", path, old, new))
    return out


@unittest.skipUnless(os.path.isdir(os.path.join(STAGE_E, "thin", "runs")) and os.path.isdir(STAGE_E_MANIFESTS), "the stage-E records are not present")
class TestStageERecordsV2(unittest.TestCase):
    """The v2 reader on the stage-E records of 1215: a single-source window value-identical to its v1 report of record (every
    pre-existing leaf; the tool string / paths / label aside), the multi-source windows' Exact share <= 100 % per source (v1 of
    record: W12 117 % / W13 114 %), W31's UnpatchedPairs ESTIMATE (~16 % of W_t) with the F2 credit capped at the clipped CSV."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="covrep-")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def report(self, w, **kwargs):
        manifest = os.path.join(STAGE_E_MANIFESTS, f"{w}-new-s2-r1p9-v3-k-2d-E", "postpro", "surface-response-requirements.json")
        bands = cr.load_bands(cr.DEFAULT_BANDS)
        return cr.build_report(os.path.join(STAGE_E, "thin", "runs", f"{w}-thin-corr", "postpro"), manifest, bands, **kwargs)

    def test_single_source_window_value_identical_to_the_report_of_record(self):
        record = os.path.join(STAGE_E, "thin", "predicted", "W06", "coverage-report.json")
        if not os.path.isfile(record) or not os.path.isfile(STAGE_E_NEAR_PAIRS):
            self.skipTest("the W06 predicted report of record is not present")
        old = _load(record)
        if "v5" not in str(old.get("bands", {}).get("version")) or "v5" not in cr.DEFAULT_BANDS:
            self.skipTest("the report of record and the default bands are not both v5")
        new = self.report("W06", near_pairs_path=STAGE_E_NEAR_PAIRS)
        self.assertEqual(new["tool"], "coverage_report.py v2")
        self.assertEqual((new["source"], new["sources_in_run"]), (1, [1]))
        self.assertNotIn("sources", new)
        ignore = {"/tool", "/bands/path", "/block_D/ladder/source", "/postpro", "/label_name"}
        self.assertEqual(_value_differences(old, new, ignore), [])

    def test_multi_source_windows_exact_share_per_source(self):
        for w, n in (("W12", 4), ("W13", 3)):
            new = self.report(w)
            self.assertEqual(new["sources_in_run"], list(range(1, n + 1)))
            for s, v in new["sources"].items():
                for t in new["types"]:
                    self.assertLessEqual(v["energy"][t]["exact_share"], 1.0 + 1e-9, f"{w} source {s} {t}")
                    self.assertGreater(v["energy"][t]["exact_share"], 0.99, f"{w} source {s} {t}")
            self.assertEqual(new["source_union"]["types"]["SA"]["n_sources"], n)
            record = os.path.join(STAGE_E, "thin", "predicted", w, "coverage-report.json")
            if os.path.isfile(record):
                self.assertGreater(_load(record)["block_A"]["energy"]["SA"]["exact_share"], 1.1)   # the v1 LAST-WINS over-read

    def test_w31_unpatched_pairs_estimate_and_the_capped_credit(self):
        new = self.report("W31")
        for t in new["types"]:
            e = new["block_A"]["energy"][t]
            self.assertTrue(e["scp_owner_credit"]["capped"])
            self.assertIn(cr.SCP_CLIPPED_CSV, e["scp_owner_credit"]["cap_source"])
            self.assertLess(e["scp_owner_credit"]["energy"] / e["W"], 1e-6)
            self.assertIn("UnpatchedPairs", e["raw_estimated"])
            self.assertGreater(e["raw_share"]["UnpatchedPairs"], 0.15)
            self.assertLess(e["raw_share"]["UnpatchedPairs"], 0.18)
            self.assertAlmostEqual(e["raw"]["Unclaimed"], 0.0, delta=1e-30)
            self.assertGreater(e["raw_share_total"], 0.22)
            self.assertLess(abs(e["bookkeeping_miss_rel"]), cr.ENERGY_BOOKKEEPING_REL_TOL)
