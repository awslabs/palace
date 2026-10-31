#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Unit tests for geometric_near_pairs.py (identification fix round DESIGN A.4 / A.5 (b); decisions 1079 Q1 (ii) / 1097
Q4): synthetic surface matrices with KNOWN correction differences (the Frobenius / trace / top-10 metrics, the flag
over the bar, the manifest pair list as the authority, the ProbeForm enumeration labelled) + the v3-j record test
(the probes' calibration of record reproduced; skipped when the library mirror is absent)."""

import json
import os
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import geometric_near_pairs as gnp  # noqa: E402
import identification_constants as ic  # noqa: E402

V3J = ("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004/v3h/stretchfix/local/library/"
       "s2-r1p9-v3-j-local/process-library-local.json")
HEADER = "interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_ij normal (J)\n"


def write_matrix(path, entries):
    """entries: {(interface, i, j): Q}; written with the production header (columns read by NAME)."""
    with open(path, "w") as f:
        f.write(HEADER)
        for (iface, i, j), q in sorted(entries.items()):
            f.write(f"{iface},0,1.9e-6,{i},{j},{q:.16e},0\n")


def make_model(root, name, topology, fab, thin, extra):
    os.makedirs(os.path.join(root, name), exist_ok=True)
    write_matrix(os.path.join(root, name, "fab.csv"), fab)
    write_matrix(os.path.join(root, name, "thin.csv"), thin)
    model = {"Name": name, "Topology": topology, "FabricatedSurfaceMatrix": f"{name}/fab.csv", "ThinSurfaceMatrix": f"{name}/thin.csv",
             "Interfaces": [{"Type": "SA", "Coupon": 1}, {"Type": "MS", "Coupon": 2}]}
    model.update(extra)
    return model


class TestGeometricNearPairs(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="gnp-")
        R = 1.9
        # D_A on SA = diag(1, 2, 3) (fab - thin); D_B = D_A x 1.01 (1 % off); D_C = D_A x 1.05 (5 % off); MS = D_A on all
        base = {(1, 1, 1): 1.0, (1, 2, 2): 2.0, (1, 3, 3): 3.0, (2, 1, 1): 0.5, (2, 2, 2): 0.5, (2, 3, 3): 0.5}
        thin = {k: 0.1 for k in base}

        def fab(scale_sa):
            return {k: (v * scale_sa if k[0] == 1 else v) + 0.1 for k, v in base.items()}
        edges = lambda d: [{"Offset": 0.0, "GapDirection": 1, "Conductor": 1}, {"Offset": 2.0 + d, "GapDirection": -1, "Conductor": 1},
                           {"Offset": 4.0, "GapDirection": 1, "Conductor": 1}]
        models = [make_model(self.tmp, "pec-A", "ParallelEdgeCluster", fab(1.0), thin, {"Edges": edges(0.0)}),
                  make_model(self.tmp, "pec-B", "ParallelEdgeCluster", fab(1.01), thin, {"Edges": edges(0.005 * R)}),   # in band
                  make_model(self.tmp, "pec-C", "ParallelEdgeCluster", fab(1.05), thin, {"Edges": edges(0.009 * R)}),   # in band, over the bar
                  make_model(self.tmp, "pec-D", "ParallelEdgeCluster", fab(1.0), thin, {"Edges": edges(0.02 * R)})]     # outside the band
        self.library = os.path.join(self.tmp, "process-library.json")
        with open(self.library, "w") as f:
            json.dump({"Version": 3, "Name": "synthetic", "MatchingRadius": R, "Models": models}, f)

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def test_probe_form_pairs_and_the_metrics(self):
        record = gnp.near_pairs_record(self.library)
        pairs = {(p["ModelA"], p["ModelB"]): p for p in record["Pairs"]}
        self.assertEqual(sorted(pairs), [("pec-A", "pec-B"), ("pec-A", "pec-C"), ("pec-B", "pec-C")])   # D is 0.02 R away from A, 0.015 from B, 0.011 from C
        ab = pairs[("pec-A", "pec-B")]
        self.assertEqual(ab["Form"], "ProbeForm")
        self.assertAlmostEqual(ab["MaxOffsetDeviationOverR"], 0.005)
        self.assertAlmostEqual(ab["CorrectionDifference"]["SA"]["Frobenius"], 0.01, places=9)
        self.assertAlmostEqual(ab["CorrectionDifference"]["SA"]["Trace"], 0.01, places=9)
        self.assertAlmostEqual(ab["CorrectionDifference"]["SA"]["Top10"], 0.01, places=9)
        self.assertAlmostEqual(ab["CorrectionDifference"]["MS"]["Max"], 0.0, places=12)
        self.assertFalse(ab["Flagged"])
        ac = pairs[("pec-A", "pec-C")]
        self.assertAlmostEqual(ac["CorrectionDifference"]["SA"]["Frobenius"], 0.05, places=9)
        self.assertTrue(ac["Flagged"])
        self.assertEqual(record["Flagged"], [["pec-A", "pec-C", ac["MaxOverTypes"]], ["pec-B", "pec-C", pairs[("pec-B", "pec-C")]["MaxOverTypes"]]])
        self.assertIn("NOT the C++", record["PairSource"])
        self.assertEqual(record["CorrectionBar"], ic.GEOMETRIC_NEAR_PAIR_CORRECTION_BAR)

    def test_manifest_pairs_are_the_authority(self):
        manifest = os.path.join(self.tmp, "manifest.json")
        with open(manifest, "w") as f:
            json.dump({"Library": {"GeometricNearPairs": [
                {"Class": "ParallelEdgeCluster", "ModelA": "pec-A", "ModelB": "pec-D", "MaxOffsetDeviationOverR": 0.02},
                {"Class": "ParallelEdgeCluster", "ModelA": "pec-A", "ModelB": "absent-model", "MaxOffsetDeviationOverR": 0.0}]}}, f)
        record = gnp.near_pairs_record(self.library, manifest)
        self.assertEqual([(p["ModelA"], p["ModelB"], p["Form"]) for p in record["Pairs"]],
                         [("pec-A", "pec-D", "Manifest"), ("pec-A", "absent-model", "Manifest")])
        self.assertAlmostEqual(record["Pairs"][0]["CorrectionDifference"]["SA"]["Frobenius"], 0.0, places=12)
        self.assertIsNone(record["Pairs"][1]["CorrectionDifference"])
        self.assertIn("absent-model", record["Pairs"][1]["Note"])
        self.assertIn("decision 1079", record["PairSource"])

    def test_different_basis_sizes_use_trace_and_spectrum(self):
        Da = {"SA": {(1, 1): 1.0, (2, 2): 2.0}}
        Db = {"SA": {(1, 1): 1.0, (2, 2): 2.0, (3, 3): 0.0}}
        diff = gnp.correction_difference(Da, 2, Db, 3)["SA"]
        self.assertFalse(diff["SameBasis"])
        self.assertNotIn("Frobenius", diff)
        self.assertAlmostEqual(diff["Trace"], 0.0)
        self.assertAlmostEqual(diff["Top10"], 0.0)

    def test_main_writes_the_record(self):
        out = os.path.join(self.tmp, "out", "pairs.json")
        rc = gnp.main(["--library", self.library, "--out", out])
        self.assertEqual(rc, 0)
        with open(out) as f:
            record = json.load(f)
        self.assertEqual(len(record["Pairs"]), 3)


class TestV3jRecord(unittest.TestCase):
    """The v3-j calibration of record (probes/v3j_fb2_response_calibration_v2.out line 18 / 35; v3j_pec_nearpair_calibration_v1.out
    line 20): reproduced by the tool's ProbeForm enumeration to 3 significant digits."""

    @unittest.skipUnless(os.path.isfile(V3J), "the v3-j library mirror is not present")
    def test_v3j_pairs_of_record(self):
        with open(V3J) as f:
            first = json.load(f)["Models"][0]
        matrix = gnp.resolve(first["FabricatedSurfaceMatrix"], V3J)
        if not os.path.isfile(matrix):
            self.skipTest("the v3-j surface matrices are not present")
        record = gnp.near_pairs_record(V3J)
        pairs = {(p["ModelA"][-12:], p["ModelB"][-12:]): p for p in record["Pairs"]}
        self.assertEqual(sorted(pairs), [("3255f4cca186", "bf97438bcd9d"), ("7d88636452b7", "5bda1bc70e63"),
                                         ("a2f53cd8630e", "7d88636452b7"), ("af972e5a4a5d", "186d141371db")])
        self.assertAlmostEqual(pairs[("a2f53cd8630e", "7d88636452b7")]["CorrectionDifference"]["SA"]["Frobenius"], 1.028e-3, delta=1e-6)
        self.assertAlmostEqual(pairs[("7d88636452b7", "5bda1bc70e63")]["CorrectionDifference"]["SA"]["Frobenius"], 8.485e-3, delta=1e-6)
        self.assertAlmostEqual(pairs[("af972e5a4a5d", "186d141371db")]["CorrectionDifference"]["SA"]["Frobenius"], 9.935e-3, delta=1e-6)
        self.assertAlmostEqual(pairs[("3255f4cca186", "bf97438bcd9d")]["CorrectionDifference"]["SA"]["Frobenius"], 2.367e-4, delta=1e-7)
        self.assertEqual(record["Flagged"], [])   # every in-band pair of v3-j is under the 2 % bar (A.4 (c) / A.5 (b))


if __name__ == "__main__":
    unittest.main()
