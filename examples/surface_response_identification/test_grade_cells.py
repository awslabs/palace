#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Tests of grade_cells.py (the scorer of record, fix round R-7): the ported letter rules on fixtures, the fail-closed cell rule
over sources / frequencies (decisions 1210 / 1213), the Q1 static-limit C form, a synthetic ES window and a synthetic Maxwell
window graded end to end against a values reference, and -- when the records are present -- the REPRODUCTION of the 47 stage-E
cells of decision 1215 (bit-exact per cell and per source / frequency), the Q1 letters on the two ported Maxwell windows and the
995 Maxwell card's 7 re-scored windows (letters identical, |dx| < 1e-6 pt; the 4 carried pembed-form windows read from the card
JSON as CARRIED: supervisor decision 1233 (1))."""

import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import grade_cells as gc  # noqa: E402

RECORDS = "/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913"
E = os.path.join(RECORDS, "readiness-20261029", "stageE")
V3H = os.path.join(RECORDS, "stage2-20261004", "v3h")
LIBRARY_E = os.path.join(E, "lib", "library", "s2-r1p9-v3-k-2d-E", "process-library.json")
SHELLS = os.path.join(E, "grade", "records", "configs")
ES_WINDOWS = ("W01p", "W03", "W05", "W06", "W10", "W11", "W12", "W13", "W31")
MX_REFERENCES = {"W06": {5.2: "W06-mx-r10-p4-a2", 5.7: "W06-mx-r10-p4-a2", 6.2: "W06-mx-r10-p4-a2"},
                 "W05": {5.2: "W05-mx-r10-p4-a3-f5p2", 5.7: "W05-mx-r10-p4-a3-f5p7", 6.2: "W05-mx-r10-p4-a3-f6p2"},
                 "FE6": {5.2: "FE6-mx-r10-p4-a2-f5p2", 5.7: "FE6-mx-r10-p4-a2-f5p7", 6.2: "FE6-mx-r10-p4-a2-f6p2"}}
# The 995 card: the 7 RE-SCORED windows (junction-driven followers; the E_cap renormalisation of the card tools) and the 4 CARRIED
# pembed-form windows (read from the card JSON, decision 1233 (1)).
CARD_995 = os.path.join(V3H, "mx-stretchfix-rescore", "results", "scorecard-11windows-mx-v3j-stretchfix.json")
RESCORED_995 = ("sct002-S1", "ctx003-C1", "ctx003-C2", "ctx003-C2p", "ctx003-C3", "sct002-S3p", "sct002-S4")
CARRIED_995 = ("ctx003-C4", "osc003-O2", "sct002-S5", "sct002-S6")
LIBRARY_V3J = os.path.join(V3H, "boxgrow", "campaign", "records", "library", "s2-r1p9-v3-j", "process-library.json")
REFRESH = os.path.join(RECORDS, "stage2-20261004", "reference-refresh-17", "results", "refresh.json")


def write_csv(path, header, rows):
    with open(path, "w") as handle:
        handle.write(",".join(header) + "\n")
        for r in rows:
            handle.write(",".join(str(v) for v in r) + "\n")


def synthetic_es_run(root, sources, p_ft, c_ft, c_off, e_elec=1.0e-12):
    """A one-cycle ES run: interfaces 1 / 2 / 3 (SA / MS / MA), one isolated-edge model carrying half of the MA energy."""
    post = os.path.join(root, "postpro")
    os.makedirs(post)
    with open(os.path.join(post, "config_resolved.json"), "w") as handle:
        json.dump({"Boundaries": {"Postprocessing": {"Dielectric": [{"Index": 1, "Type": "SA"}, {"Index": 2, "Type": "MS"}, {"Index": 3, "Type": "MA"}]}}}, handle)
    with open(os.path.join(post, "palace.json"), "w") as handle:
        json.dump({"SurfaceResponse": {"ModelCatalog": [{"Index": 0, "Name": "isolated-edge"}]}}, handle)
    hq = ["i", "E_elec raw (J)", "E_elec postprocessed fixed-trace (J)"]
    for k in (1, 2, 3):
        hq += [f"p_surf raw[{k}]", f"p_surf postprocessed fixed-trace[{k}]", f"E_surf raw[{k}] (J)", f"E_surf postprocessed fixed-trace[{k}] (J)"]
    rows = []
    for s in sources:
        r = [s, e_elec, e_elec]
        for k in (1, 2, 3):
            p = p_ft[s][k]
            r += [0.9 * p, p, 0.9 * p * e_elec, p * e_elec]
        rows.append(r)
    write_csv(os.path.join(post, "surface-Q-corrected.csv"), hq, rows)
    write_csv(os.path.join(post, "surface-Q-edge.csv"), ["i", "interface", "R (m)", "E_out (J)"],
              [[s, k, 1.9e-6, 0.0] for s in sources for k in (1, 2, 3)])
    hc = ["i"] + [f"C fixed-trace[i][{j}] (F)" for j in sources]
    # the off-diagonal entries equal the reference's (-0.1 x the reference diagonal c_off): the untargeted rows read 0
    write_csv(os.path.join(post, "terminal-C-corrected.csv"), hc, [[s] + [c_ft[s] if j == s else c_off[s] for j in sources] for s in sources])
    hr = ["i"] + [f"C[i][{j}] (F)" for j in sources]
    write_csv(os.path.join(post, "terminal-C.csv"), hr, [[s] + [0.95 * c_ft[s] if j == s else c_off[s] for j in sources] for s in sources])
    write_csv(os.path.join(post, "surface-response-model-energy.csv"),
              ["source", "evaluation", "model", "fabricated surface energy[1] (J)", "fabricated surface energy[2] (J)", "fabricated surface energy[3] (J)"],
              [[s, ev, 0, 0.0, 0.0, 0.5 * p_ft[s][3] * e_elec] for s in sources for ev in (0, 1)])
    return root


def es_values_reference(sources, p_ref, c_ref, ma_sharp_ref):
    return {"members": {str(s): {"SA_p": p_ref[s][1], "MS_p": p_ref[s][2], "MA_p": p_ref[s][3], "MA_sharp_p": ma_sharp_ref[s],
                                 "C_F": c_ref[s], "E_elec_J": 1.0e-12, "C_offdiag_F": {str(j): -0.1 * c_ref[s] for j in sources if j != s}}
                        for s in sources}}


class LetterRules(unittest.TestCase):
    def test_verdict(self):
        self.assertEqual(gc.verdict(4.9, 5.0), "PASS")
        self.assertEqual(gc.verdict(-5.1, 5.0), "FAIL")
        self.assertEqual(gc.verdict(1.0, 5.0, [-6.0, 2.0]), "INCONCLUSIVE (ft inside, bracket outside)")
        self.assertEqual(gc.verdict(1.0, 5.0, [-4.0, 2.0]), "PASS")
        self.assertEqual(gc.verdict(float("nan"), 5.0), "n/a")

    def test_model_class_and_tails(self):
        self.assertEqual(gc.model_class("convex-corner-90deg"), "corner")
        self.assertEqual(gc.model_class("parallel-edge-cluster-95b13e7fdc7c-pieces3"), "stack (ParallelEdgeCluster)")
        self.assertEqual(gc.model_class("spatialedgecluster_edgecount-10_13c7cc02fd9b"), "spatial cluster 10-edge")
        self.assertEqual(gc.model_class("different-conductor-gap-3um"), "gap pair")
        self.assertEqual(gc.tail_of("convex-corner-90deg", {})[0], gc.TAILS["corner"])
        self.assertEqual(gc.tail_of("isolated-edge", {})[0], gc.TAILS["isolated"])
        self.assertEqual(gc.tail_of("spatialedgecluster_edgecount-3_abc@1", {"spatialedgecluster_edgecount-3_abc": (0.02, 0.01, 0.03)})[0], (0.02, 0.01, 0.03))
        self.assertEqual(gc.tail_of("spatialedgecluster_edgecount-3_zzz", {})[0], gc.SPATIAL_DEFAULT_TAIL)

    def test_ma_sharp_convention_b(self):
        sharp = gc.ma_sharp_convention_b(2.0, {0: 1.0}, {0: "isolated-edge"}, {})
        self.assertAlmostEqual(sharp["device_sharp_E_J"], 1.0 + 1.0678)
        self.assertEqual(sharp["bracket_E_J"], [1.0 + 1.06, 1.0 + 1.0726])
        spatial = gc.ma_sharp_convention_b(2.0, {0: 1.0}, {0: "spatialedgecluster_edgecount-3_abc"}, {"spatialedgecluster_edgecount-3_abc": (0.03, 0.02, 0.04)})
        self.assertEqual(spatial["bracket_E_J"], [2.0, 1.0 + 1.06])  # the spatial bar [0, 2 x record]


class CellRules(unittest.TestCase):
    def grades(self, xs, letters):
        return {m: {"x_pct": {"SA": x}, "letters": {"SA": l}, "band": {"SA": {"in": x < 1.0}}} for m, (x, l) in zip(xs, zip(xs, letters))}

    def test_fail_closed_over_members(self):
        g = self.grades([0.5, -5.5, 2.0], ["PASS", "FAIL", "PASS"])
        cell = gc.aggregate_cell([0.5, -5.5, 2.0], g, "SA", "rule")
        self.assertEqual(cell["letter"], "FAIL")
        self.assertEqual(cell["x_pct"], -5.5)
        self.assertEqual(cell["x_range_pct"], [-5.5, 2.0])
        self.assertEqual((cell["n_pass"], cell["n_fail"], cell["n_members"]), (2, 1, 3))
        self.assertFalse(cell["in_band"])
        g = self.grades([0.5, 0.2], ["PASS", "INCONCLUSIVE (ft inside, bracket outside)"])
        self.assertTrue(gc.aggregate_cell([0.5, 0.2], g, "SA", "rule")["letter"].startswith("INCONCLUSIVE"))
        g = self.grades([0.5, 0.2], ["PASS", "PASS"])
        cell = gc.aggregate_cell([0.5, 0.2], g, "SA", "rule")
        self.assertEqual(cell["letter"], "PASS")
        self.assertTrue(cell["in_band"])

    def test_static_limit_fit(self):
        freqs = [5.2, 5.7, 6.2]
        c_s, a = 2.0e-13, 0.004
        fit = gc.static_limit_fit(freqs, [c_s * (1.0 + a * f * f) for f in freqs])
        self.assertAlmostEqual(fit["C_s"] / c_s, 1.0, places=12)
        self.assertAlmostEqual(fit["a_per_GHz2"], a, places=12)
        self.assertLess(fit["max_residual_pct"], 1e-9)
        with self.assertRaises(ValueError):
            gc.static_limit_fit([5.2], [1.0])


class SyntheticWindows(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_es_two_sources_fail_closed(self):
        sources = [1, 2]
        p_ref = {1: {1: 1.0e-4, 2: 2.0e-4, 3: 1.0e-5}, 2: {1: 1.0e-4, 2: 2.0e-4, 3: 1.0e-5}}
        # source 1 inside every bar; source 2's SA 6 % high -> the SA cell FAILS fail-closed, the others PASS
        p_ft = {1: {1: 1.01e-4, 2: 2.0e-4, 3: 1.0e-5}, 2: {1: 1.06e-4, 2: 2.0e-4, 3: 1.0e-5}}
        c_ref = {1: 1.0e-14, 2: 2.0e-14}
        c_ft = {1: 1.001e-14, 2: 2.004e-14}
        # the device MA_sharp (B) = (E_MA - F) + F (1 + 0.0678) with F = 0.5 E_MA -> p_MA (1 + 0.0339): the reference set so x = 0
        ma_sharp_ref = {s: p_ref[s][3] * (1.0 + 0.5 * 0.0678) for s in sources}
        run = synthetic_es_run(os.path.join(self.tmp, "run"), sources, p_ft, c_ft, {s: -0.1 * c_ref[s] for s in sources})
        values = os.path.join(self.tmp, "ref.json")
        with open(values, "w") as handle:
            json.dump(es_values_reference(sources, p_ref, c_ref, ma_sharp_ref), handle)
        out = os.path.join(self.tmp, "grade.json")
        rc = gc.main(["--form", "es", "--device", run, "--reference-values", values, "--out", out, "--label", "syn"])
        self.assertEqual(rc, 0)
        with open(out) as handle:
            g = json.load(handle)
        self.assertEqual(g["members"], [1, 2])
        self.assertEqual(g["cells"]["SA"]["letter"], "FAIL")
        self.assertEqual(g["cells"]["SA"]["worst_member"], 2)
        self.assertAlmostEqual(g["cells"]["SA"]["x_pct"], 6.0, places=9)
        self.assertEqual(g["cells"]["SA"]["letters_by_member"], {"1": "PASS", "2": "FAIL"})
        self.assertEqual(g["cells"]["MS"]["letter"], "PASS")
        self.assertEqual(g["cells"]["MA_sharp"]["letter"], "PASS")
        self.assertAlmostEqual(g["cells"]["MA_sharp"]["x_pct"], 0.0, places=9)
        self.assertEqual(g["cells"]["C"]["letter"], "PASS")
        self.assertAlmostEqual(g["cells"]["C"]["x_pct"], 0.2, places=9)
        # the off-diagonal C rows (910): identical on both sides -> 0 absolute deviation
        off = g["per_member"]["1"]["grade"]["C_offdiag"]
        self.assertEqual(list(off), ["2"])
        self.assertAlmostEqual(off["2"]["abs_over_Cii_pct"], 0.0)
        # a single source graded with --sources 1 -> the single-source rule
        rc = gc.main(["--form", "es", "--device", run, "--reference-values", values, "--out", out, "--sources", "1"])
        self.assertEqual(rc, 0)
        with open(out) as handle:
            g1 = json.load(handle)
        self.assertEqual(g1["cells"]["SA"]["rule"], "single source")
        self.assertEqual(g1["cells"]["SA"]["letter"], "PASS")

    def test_missing_input_exit_2(self):
        rc = gc.main(["--form", "es", "--device", os.path.join(self.tmp, "absent"), "--reference-values", os.path.join(self.tmp, "x.json"), "--out",
                      os.path.join(self.tmp, "g.json")])
        self.assertEqual(rc, 2)

    def synthetic_mx_run(self, root, freqs, p_ft, v2, e_norm=1.0e-12, e_cap=0.0):
        post = os.path.join(root, "postpro")
        os.makedirs(post)
        with open(os.path.join(post, "config_resolved.json"), "w") as handle:
            json.dump({"Boundaries": {"Postprocessing": {"Dielectric": [{"Index": 1, "Type": "SA"}, {"Index": 2, "Type": "MS"}, {"Index": 3, "Type": "MA"}]}}}, handle)
        with open(os.path.join(post, "palace.json"), "w") as handle:
            json.dump({"SurfaceResponse": {"ModelCatalog": [{"Index": 0, "Name": "isolated-edge"}]}}, handle)
        hq = ["f (GHz)", "exc", "interface", "E_norm raw (J)", "E_norm corrected (J)", "E_surf raw (J)", "p_surf raw", "E_surf corrected (J)", "p_surf corrected",
              "E_surf corrected fixed-flux (J)", "p_surf corrected fixed-flux"]
        rows = []
        for F in freqs:
            for k in (1, 2, 3):
                p = p_ft[F][k]
                rows.append([F, 1, k, e_norm, e_norm, 0.9 * p * e_norm, 0.9 * p, p * e_norm, p, p * e_norm, p])
        write_csv(os.path.join(post, "surface-Q-corrected.csv"), hq, rows)
        write_csv(os.path.join(post, "surface-Q-edge.csv"), ["f (GHz)", "interface", "R (m)", "E_out (J)"], [[F, k, 1.9e-6, 0.0] for F in freqs for k in (1, 2, 3)])
        write_csv(os.path.join(post, "domain-E.csv"), ["f (GHz)", "E_elec (J)", "E_mag (J)", "E_cap (J)"], [[F, e_norm, 0.5 * e_norm, e_cap] for F in freqs])
        write_csv(os.path.join(post, "port-V.csv"), ["f (GHz)", "Re{V[1]} (V)", "Im{V[1]} (V)"], [[F, v2[F] ** 0.5, 0.0] for F in freqs])
        write_csv(os.path.join(post, "surface-response-model-energy.csv"),
                  ["f (GHz)", "source", "evaluation", "model", "fabricated surface energy[1] (J)", "fabricated surface energy[2] (J)", "fabricated surface energy[3] (J)"],
                  [[F, 1, ev, 0, 0.0, 0.0, 0.5 * p_ft[F][3] * e_norm] for F in freqs for ev in (0, 1)])
        return root

    def test_maxwell_c_forms(self):
        freqs = [5.2, 5.7, 6.2]
        p = {1: 1.0e-4, 2: 2.0e-4, 3: 1.0e-5}
        p_ft = {F: dict(p) for F in freqs}
        c_s_thin, c_s_thick, a_thin, a_thick = 2.004e-13, 2.0e-13, 0.0060, 0.0055
        e_norm = 1.0e-12
        v2 = {F: 2.0 * e_norm / (c_s_thin * (1.0 + a_thin * F * F)) for F in freqs}
        run = self.synthetic_mx_run(os.path.join(self.tmp, "mx"), freqs, p_ft, v2)
        ref = {"members": {str(F): {"SA_p": p[1], "MS_p": p[2], "MA_p": p[3], "MA_sharp_p": p[3] * (1.0 + 0.5 * 0.0678), "E_elec_J": e_norm,
                                    "C_eff_F": c_s_thick * (1.0 + a_thick * F * F)} for F in freqs}}
        values = os.path.join(self.tmp, "ref.json")
        with open(values, "w") as handle:
            json.dump(ref, handle)
        out = os.path.join(self.tmp, "g.json")
        # direct: C_eff(f) grows past the 0.5 % bar at the upper frequencies -> FAIL fail-closed
        self.assertEqual(gc.main(["--form", "maxwell", "--device", run, "--reference-values", values, "--c-form", "direct", "--out", out]), 0)
        with open(out) as handle:
            direct = json.load(handle)
        self.assertEqual(direct["cells"]["C"]["letter"], "FAIL")
        self.assertEqual(direct["cells"]["C"]["worst_member"], 6.2)
        self.assertEqual(direct["cells"]["SA"]["letter"], "PASS")
        self.assertEqual(direct["cells"]["SA"]["rule"][:28], "every frequency must PASS (f")
        # static-limit (Q1): C_s thin / C_s thick - 1 = +0.2 % -> PASS; C_eff(f) printed beside
        self.assertEqual(gc.main(["--form", "maxwell", "--device", run, "--reference-values", values, "--c-form", "static-limit", "--out", out]), 0)
        with open(out) as handle:
            sl = json.load(handle)
        self.assertEqual(sl["cells"]["C"]["letter"], "PASS")
        self.assertAlmostEqual(sl["cells"]["C"]["x_pct"], 0.2, places=6)
        self.assertEqual(sorted(sl["cells"]["C"]["C_eff_by_f_pct"]), ["5.2", "5.7", "6.2"])
        self.assertEqual(sl["cells"]["C"]["direct_form_beside"]["letter"], "FAIL")
        # information (F-40): not counted
        self.assertEqual(gc.main(["--form", "maxwell", "--device", run, "--reference-values", values, "--c-form", "information", "--out", out]), 0)
        with open(out) as handle:
            info = json.load(handle)
        self.assertFalse(info["cells"]["C"]["counted"])
        self.assertTrue(info["cells"]["C"]["letter"].startswith("INFORMATION (F-40)"))
        # the E_cap renormalisation (the junction-driven form): p = E / (E_norm - E_cap)
        run2 = self.synthetic_mx_run(os.path.join(self.tmp, "mx2"), [1.0], {1.0: dict(p)}, {1.0: 1.0}, e_norm=e_norm, e_cap=0.5 * e_norm)
        ref2 = {"members": {"1.0": {"SA_p": 2.0 * p[1], "MS_p": 2.0 * p[2], "MA_p": 2.0 * p[3], "MA_sharp_p": 2.0 * p[3] * (1.0 + 0.5 * 0.0678), "E_elec_J": e_norm}}}
        with open(values, "w") as handle:
            json.dump(ref2, handle)
        self.assertEqual(gc.main(["--form", "maxwell", "--device", run2, "--reference-values", values, "--c-form", "information", "--frequencies", "1.0",
                                  "--p-denominator", "E_norm-minus-E_cap", "--out", out]), 0)
        with open(out) as handle:
            cap = json.load(handle)
        for T in ("SA", "MS", "MA_sharp"):
            self.assertAlmostEqual(cap["cells"][T]["x_pct"], 0.0, places=9)
        self.assertEqual(gc.main(["--form", "maxwell", "--device", run2, "--reference-values", values, "--c-form", "information", "--frequencies", "1.0", "--out", out]), 0)
        with open(out) as handle:
            nocap = json.load(handle)
        self.assertAlmostEqual(nocap["cells"]["SA"]["x_pct"], -50.0, places=9)


def load(path):
    with open(path) as handle:
        return json.load(handle)


def refresh_values(window, path):
    """The Richardson reference of record (refresh.json, decision 423) as a grade_cells values file at 1 GHz."""
    w = load(REFRESH)["Windows"][window]
    e_elec = w["Levels"]["r2p5p5"]["E_elec_J"]
    member = {"SA_p": w["Richardson"]["SA"]["E_inf_J"] / e_elec, "MS_p": w["Richardson"]["MS"]["E_J"] / e_elec,
              "MA_p": w["Richardson"]["MA"]["E_inf_J"] / e_elec, "MA_sharp_p": w["MA_sharp"]["r2p5p5"]["MA_sharp_J"] / e_elec,
              "C_F": w["Richardson"]["C"]["C_F"], "E_elec_J": e_elec}
    with open(path, "w") as handle:
        json.dump({"source": REFRESH, "window": window, "members": {"1.0": member}}, handle)
    return path


@unittest.skipUnless(os.path.isdir(os.path.join(E, "thin", "runs")) and os.path.isfile(LIBRARY_E), "the stage-E records are not present")
class StageERecords(unittest.TestCase):
    """The 47 cells of decision 1215 reproduced from the stored CSVs: per cell and per source / frequency x_pct, letter, in_band."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def es_args(self, w):
        args = ["--form", "es", "--device", os.path.join(E, "thin", "runs", f"{w}-thin-corr"), "--reference", os.path.join(E, "ref", "runs", "ref", f"{w}-ref-r10-p4"),
                "--shells-dir", SHELLS, "--library", LIBRARY_E, "--band", os.path.join(E, "thin", "predicted", w, "coverage-report.json"), "--label", w,
                "--out", os.path.join(self.tmp, f"{w}.json")]
        if not os.path.isfile(os.path.join(E, "thin", "runs", f"{w}-thin-corr", "postpro", "palace.json")):
            args += ["--model-catalog", os.path.join(E, "grade", "records", f"{w}-thin-corr-model-catalog.json")]
        return args

    def mx_args(self, w, c_form):
        args = ["--form", "maxwell", "--device", os.path.join(E, "mxthin", "runs", f"{w}-mxthin-corr"), "--shells-dir", SHELLS, "--library", LIBRARY_E,
                "--band", os.path.join(E, "mxthin", "predicted", w, "coverage-report.json"), "--c-form", c_form, "--label", f"{w}-mx",
                "--out", os.path.join(self.tmp, f"{w}-mx-{c_form}.json")]
        for F, run in MX_REFERENCES[w].items():
            args += ["--reference-at", f"{F}={os.path.join(E, 'ref', 'runs', 'mx', run)}"]
        return args

    def assert_cells_identical(self, got, record, members_key):
        for T in ("SA", "MS", "MA_sharp", "C"):
            g, r = got["cells"][T], record["cells"][T]
            self.assertEqual((g["x_pct"], g["letter"], g["in_band"]), (r["x_pct"], r["letter"], r["in_band"]), f"{record['window']} {T}")
        for m, v in record[members_key].items():
            g = got["per_member"][m]["grade"]
            self.assertEqual(g["x_pct"], v["grade"]["x_pct"], f"{record['window']} member {m}")
            self.assertEqual(g["letters"], v["grade"]["letters"], f"{record['window']} member {m}")

    def test_es_36_cells(self):
        for w in ES_WINDOWS:
            if w == "W10":
                continue  # the 46-source window is its own test (the longest read)
            self.assertEqual(gc.main(self.es_args(w)), 0)
            self.assert_cells_identical(load(os.path.join(self.tmp, f"{w}.json")), load(os.path.join(E, "grade", "results", f"{w}-grade.json")), "per_source")

    def test_es_w10_46_sources(self):
        self.assertEqual(gc.main(self.es_args("W10")), 0)
        got = load(os.path.join(self.tmp, "W10.json"))
        self.assertEqual(len(got["members"]), 46)
        self.assert_cells_identical(got, load(os.path.join(E, "grade", "results", "W10-grade.json")), "per_source")

    def test_maxwell_11_cells_direct_form(self):
        for w, c_form in (("W06", "direct"), ("W05", "direct"), ("FE6", "information")):
            self.assertEqual(gc.main(self.mx_args(w, c_form)), 0)
            self.assert_cells_identical(load(os.path.join(self.tmp, f"{w}-mx-{c_form}.json")), load(os.path.join(E, "grade", "results", f"{w}-mx-grade.json")), "per_frequency")

    def test_maxwell_q1_static_limit(self):
        """USER 1230 Q1: the static-limit C is the graded cell on the two ported lines; the pre-registered form change moves exactly these two letters."""
        expected = {"W06": 0.330861738517374, "W05": 0.3165595747347316}
        for w, x in expected.items():
            self.assertEqual(gc.main(self.mx_args(w, "static-limit")), 0)
            cell = load(os.path.join(self.tmp, f"{w}-mx-static-limit.json"))["cells"]["C"]
            self.assertEqual(cell["letter"], "PASS")
            self.assertAlmostEqual(cell["x_pct"], x, places=9)
            self.assertTrue(cell["in_band"])
            self.assertEqual(cell["direct_form_beside"]["letter"], "FAIL")
            self.assertEqual(sorted(cell["C_eff_by_f_pct"]), ["5.2", "5.7", "6.2"])


@unittest.skipUnless(os.path.isfile(CARD_995) and os.path.isfile(LIBRARY_V3J) and os.path.isfile(REFRESH), "the 995 card records are not present")
class MaxwellCard995(unittest.TestCase):
    """The 995 Maxwell card (43 P / 1 I / 0 F; G2): the 7 RE-SCORED windows re-scored by grade_cells.py from the stored junction-driven
    followers (p = E / (E_norm - E_cap), the A1.1.1 renormalisation of the card tools) against refresh.json's Richardson references
    -> the card's letters, |dx| < 1e-6 pt; the 4 CARRIED pembed-form windows are read from the card JSON as CARRIED (decision 1233 (1))."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.card = {row["window"]: row for row in load(CARD_995)["rows"]}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_seven_rescored_windows(self):
        runs = os.path.join(V3H, "mx-stretchfix-rescore", "runs")
        for w in RESCORED_995:
            run = os.path.join(runs, f"mxp-{w}-driven-p4-1GHz-jjdrive-follower-v3j-mxsf")
            if not os.path.isdir(run):
                self.skipTest(f"{run} is not present")
            values = refresh_values(w, os.path.join(self.tmp, f"ref-{w}.json"))
            out = os.path.join(self.tmp, f"{w}.json")
            rc = gc.main(["--form", "maxwell", "--device", run, "--reference-values", values, "--library", LIBRARY_V3J, "--frequencies", "1.0",
                          "--p-denominator", "E_norm-minus-E_cap", "--c-form", "information", "--label", w, "--out", out])
            self.assertEqual(rc, 0)
            got = load(out)["cells"]
            for T in ("SA", "MS", "MA_sharp"):
                card = self.card[w]["cells"][T]
                self.assertEqual(got[T]["letter"][0], card["letter"], f"{w} {T}")
                self.assertLess(abs(got[T]["x_pct"] - card["value"]), 1e-6, f"{w} {T}: {got[T]['x_pct']} vs the card {card['value']}")

    def test_four_carried_windows_read_from_the_card(self):
        for w in CARRIED_995:
            row = self.card[w]
            self.assertTrue(row["status"].startswith("CARRIED"), w)
            for T in ("SA", "MS", "MA_sharp"):
                self.assertEqual(row["cells"][T]["letter"], "P", f"{w} {T}")


if __name__ == "__main__":
    unittest.main()
