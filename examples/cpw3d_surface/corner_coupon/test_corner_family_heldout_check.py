#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The corner family's held-out interpolation gate (corner_family_heldout_check.py; gate
CornerFamilyInterpolation of spatial_coupon/qualify/qualification-gates.json, USER decision
149 (5): participation-referenced 0.5 % of the fabricated energy, the defect-referenced form
reported alongside) on synthetic coupon directories whose matrices are polynomials of the
angle: a cubic segment stencil reproduces a cubic exactly, a perturbed held-out coupon fails
only when the perturbation exceeds the gate in the participation-referenced form, and a
held-out angle the stencil rule refuses fails."""

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CHECK = load("corner_family_heldout_check")
SIZE = 4
GATES = json.load(open(CHECK.GATES_FILE))["Gates"][CHECK.GATE]


def matrices(angle, scale=1.0):
    """Cubic-in-angle symmetric positive matrices: thin and fabricated domain and per-interface
    (fabricated - thin small: a 5 % defect of the SA / MS / MA matrices)."""
    t = (angle - 100.0) / 50.0
    base = np.eye(SIZE) + 0.1 * np.ones((SIZE, SIZE))
    thin = base * (1.0 + 0.3 * t + 0.2 * t**2 + 0.1 * t**3)
    fab = (1.05 * thin + 0.02 * t**3 * base) * scale
    interfaces = {k: (thin * factor, fab * factor) for k, factor in ((1, 1e-3), (2, 2e-3), (3, 5e-5))}
    return thin, fab, interfaces


def write_coupon(root, topology, angle, connectivity, scale=1.0):
    d = root / f"{topology}-{angle:g}-{connectivity}"
    for kind in ("thin", "fabricated"):
        (d / "postpro" / kind).mkdir(parents=True)
    thin, fab, interfaces = matrices(angle, scale)
    for kind, domain in (("thin", thin), ("fabricated", fab)):
        rows = ["basis_i,basis_j,Q_ij (J)"]
        agg = ["interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_total_ij (J)"]
        for i in range(SIZE):
            for j in range(i, SIZE):
                rows.append(f"{i + 1},{j + 1},{domain[i, j]:.16e}")
                for k, (thin_k, fab_k) in interfaces.items():
                    q = (thin_k if kind == "thin" else fab_k)[i, j]
                    agg.append(f"{k},1,1.9e-06,{i + 1},{j + 1},{q:.16e},{q:.16e}")
        (d / "postpro" / kind / "domain-response-matrix.csv").write_text("\n".join(rows) + "\n")
        (d / "postpro" / kind / "surface-response-matrix-aggregate.csv").write_text("\n".join(agg) + "\n")
    spec = {"Topology": topology, "AngleDegrees": angle}
    model = {"Name": f"{topology}-{angle:g}", "TraceBasis": {}}
    if connectivity is not None:
        spec["ConnectivityAngleDegrees"] = connectivity
        model["TraceBasis"]["ConnectivityAngleDegrees"] = connectivity
    (d / "coupon-spec.json").write_text(json.dumps(spec))
    (d / "process-library.json").write_text(json.dumps({"Models": [model]}))
    np.savetxt(d / "heldout-coefficients.csv", np.array([0.3, -0.2, 0.5, 0.1]), delimiter=",",
               header="coefficient_V", comments="")
    return d


class CornerFamilyHeldoutCheckTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def segment(self, scale=1.0, heldout_angle=112.0, heldout_connectivity=None):
        # The convex [90, 135] segment (keys 112.5): four nodes, one held-out coupon.
        root = Path(tempfile.mkdtemp(dir=self.root))
        dirs = [write_coupon(root, "convex", a, 112.5) for a in (90.0, 105.0, 120.0, 135.0)]
        dirs.append(write_coupon(root, "convex", heldout_angle, heldout_connectivity, scale))
        return dirs

    def run_check(self, dirs, nodes=None):
        out = self.root / "out.csv"
        record = self.root / "out.json"
        command = [sys.executable, str(ROOT / "corner_family_heldout_check.py"), str(out), *map(str, dirs),
                   "--json", str(record)]
        if nodes is not None:
            command += ["--nodes", *nodes]
        result = subprocess.run(command, capture_output=True, text=True)
        return result.returncode, json.loads(record.read_text()), result.stdout

    def test_gate_definition(self):
        self.assertEqual(GATES["MaximumParticipationReferencedResidual"], 0.005)
        self.assertIn("participation-referenced", GATES["Statement"])
        self.assertIn("DefectReferenced", GATES["ReportedAlongside"])

    def test_cubic_segment_reproduces_a_cubic_held_out_coupon(self):
        code, record, _ = self.run_check(self.segment())
        self.assertEqual(code, 0, record)
        self.assertTrue(record["Passed"])
        family = record["Families"]["convex"]
        self.assertEqual(len(family["Nodes"]), 4)
        self.assertEqual(len(family["HeldOut"]), 1)
        entry = family["HeldOut"][0]
        self.assertEqual(entry["Rule"], "cubic")
        self.assertEqual(entry["Stencil"], [90.0, 105.0, 120.0, 135.0])
        for value in entry["ParticipationReferencedPercent"].values():
            self.assertAlmostEqual(value, 0.0, places=9)
        for value in entry["DefectReferencedPercent"].values():
            self.assertAlmostEqual(value, 0.0, places=7)
        self.assertEqual(record["GatesSHA256"], CHECK.hashlib.sha256(CHECK.GATES_FILE.read_bytes()).hexdigest())
        self.assertEqual(record["MaximumParticipationReferencedResidual"], 0.005)

    def test_participation_referenced_residual_decides_the_verdict(self):
        # A held-out coupon whose fabricated matrices are 0.4 % off: participation-referenced
        # -0.4 % passes the 0.5 % gate although the DEFECT-referenced SA residual (the
        # correction is 5 % of the matrix) is ~ -8 %: reported, not gated.
        code, record, _ = self.run_check(self.segment(scale=1.004))
        entry = record["Families"]["convex"]["HeldOut"][0]
        self.assertEqual(code, 0)
        self.assertAlmostEqual(entry["ParticipationReferencedPercent"]["SA"], -100.0 * 0.004 / 1.004, places=6)
        self.assertLess(entry["DefectReferencedPercent"]["SA"], -5.0)
        self.assertTrue(record["Passed"])
        # 0.6 % off fails.
        code, record, _ = self.run_check(self.segment(scale=1.006))
        self.assertEqual(code, 1)
        self.assertFalse(record["Passed"])
        self.assertFalse(record["Families"]["convex"]["HeldOut"][0]["Passed"])

    def test_refused_held_out_angle_fails(self):
        # 80 deg lies below the segment's node range: the stencil rule refuses (no
        # extrapolation) and the family cannot pass.
        code, record, _ = self.run_check(self.segment(heldout_angle=80.0))
        self.assertEqual(code, 1)
        family = record["Families"]["convex"]
        self.assertEqual(len(family["Refused"]), 1)
        self.assertIn("sharper", family["Refused"][0]["Reason"])
        self.assertFalse(family["Passed"])

    def test_nodes_option_stamps_legacy_coupons(self):
        # Recorded legacy coupons (no connectivity record) become the segment's nodes through
        # --nodes ANGLE:CONNECTIVITY; the coupon left out is held out.
        dirs = [write_coupon(self.root, "convex", a, None) for a in (90.0, 105.0, 120.0, 135.0)]
        dirs.append(write_coupon(self.root, "convex", 112.0, None))
        code, record, _ = self.run_check(dirs, nodes=["90:112.5", "105:112.5", "120:112.5", "135:112.5"])
        self.assertEqual(code, 0, record)
        self.assertEqual([n["ConnectivityAngleDegrees"] for n in record["Families"]["convex"]["Nodes"]], [112.5] * 4)
        self.assertEqual(record["Families"]["convex"]["HeldOut"][0]["AngleDegrees"], 112.0)
        # Without --nodes a legacy family has no nodes: every coupon is held out and refused.
        code, record, _ = self.run_check(dirs)
        self.assertEqual(code, 1)
        self.assertEqual(len(record["Families"]["convex"]["Refused"]), 5)


if __name__ == "__main__":
    unittest.main()
