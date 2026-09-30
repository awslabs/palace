#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The corner family's held-out interpolation gate (corner_family_heldout_check.py; gate
CornerFamilyInterpolation of spatial_coupon/qualify/qualification-gates.json, USER decision
149 (5): participation-referenced 0.5 % of the fabricated energy, the defect-referenced form
reported alongside) on synthetic coupon directories whose matrices are polynomials of the
angle: a cubic segment stencil reproduces a cubic exactly, a perturbed held-out coupon fails
only when the perturbation exceeds the gate in the participation-referenced form, and a
held-out angle the stencil rule refuses fails. The held-out trace is recorded (TraceForm):
the cache's coefficients are classified band | option-c against the generator's forms at
the basis points, --trace option-c recomputes the option-(c) coefficients from
basis-points.csv + the coupon spec, and the fail-closed guards (a coefficient vector of the
wrong length, a zero fabricated held-out energy) error out."""

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
GENERATOR = load("generate_corner_response")
SIZE = 4
RADIUS, THICKNESS = 1.9, 0.1
# Four basis points on the matching box: two on the metal rings (z = 0 and z = THICKNESS, on the
# free arc of a 90-degree convex corner) and two on the cap ring at z = -RADIUS.
BASIS_POINTS = np.array([[-RADIUS, 0.4 * RADIUS, 0.0], [0.5 * RADIUS, -RADIUS, THICKNESS],
                         [-RADIUS, 0.0, -RADIUS], [RADIUS, RADIUS, -RADIUS]])
GATES = json.load(open(CHECK.GATES_FILE))["Gates"][CHECK.GATE]


def matrices(angle, scale=1.0, ring_perturbation=0.0):
    """Cubic-in-angle symmetric positive matrices: thin and fabricated domain and per-interface
    (fabricated - thin small: a 5 % defect of the SA / MS / MA matrices). ring_perturbation
    adds that fraction of the base matrix on the block of the two metal-ring basis functions
    (BASIS_POINTS[:2]) of the fabricated matrices only: invisible to a trace that is zero there."""
    t = (angle - 100.0) / 50.0
    base = np.eye(SIZE) + 0.1 * np.ones((SIZE, SIZE))
    thin = base * (1.0 + 0.3 * t + 0.2 * t**2 + 0.1 * t**3)
    ring_block = np.zeros((SIZE, SIZE))
    ring_block[:2, :2] = base[:2, :2]
    fab = (1.05 * thin + 0.02 * t**3 * base) * scale + ring_perturbation * ring_block
    interfaces = {k: (thin * factor, fab * factor) for k, factor in ((1, 1e-3), (2, 2e-3), (3, 5e-5))}
    return thin, fab, interfaces


def write_coupon(root, topology, angle, connectivity, scale=1.0, trace=None, basis_points=None,
                 ring_perturbation=0.0, trace_basis=None):
    """trace = the recorded heldout-coefficients.csv (default a fixed vector); basis_points =
    basis-points.csv and the spec entries the generator's held-out traces need; trace_basis =
    the model's TraceBasis record (default an empty record, a connectivity angle when given)."""
    d = root / f"{topology}-{angle:g}-{connectivity}"
    for kind in ("thin", "fabricated"):
        (d / "postpro" / kind).mkdir(parents=True)
    thin, fab, interfaces = matrices(angle, scale, ring_perturbation)
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
    if basis_points is not None:
        spec["MatchingRadius"] = RADIUS
        spec["Fabrication"] = {"metal_thickness": THICKNESS}
        np.savetxt(d / "basis-points.csv", basis_points, delimiter=",", header="x,y,z", comments="", fmt="%.16e")
    model = {"Name": f"{topology}-{angle:g}", "TraceBasis": dict(trace_basis or {})}
    if connectivity is not None:
        spec["ConnectivityAngleDegrees"] = connectivity
        model["TraceBasis"]["ConnectivityAngleDegrees"] = connectivity
    (d / "coupon-spec.json").write_text(json.dumps(spec))
    (d / "process-library.json").write_text(json.dumps({"Models": [model]}))
    np.savetxt(d / "heldout-coefficients.csv", np.array([0.3, -0.2, 0.5, 0.1]) if trace is None else trace,
               delimiter=",", header="coefficient_V", comments="", fmt="%.16e")
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

    def run_check(self, dirs, nodes=None, trace=None):
        out = self.root / "out.csv"
        record = self.root / "out.json"
        if record.exists():
            record.unlink()
        command = [sys.executable, str(ROOT / "corner_family_heldout_check.py"), str(out), *map(str, dirs),
                   "--json", str(record)]
        if nodes is not None:
            command += ["--nodes", *nodes]
        if trace is not None:
            command += ["--trace", trace]
        result = subprocess.run(command, capture_output=True, text=True)
        if not record.exists():
            return result.returncode, None, result.stdout + result.stderr
        return result.returncode, json.loads(record.read_text()), result.stdout

    def family_with_basis_points(self, heldout_trace, ring_perturbation=0.0, heldout_angle=112.0):
        # The convex [90, 135] segment with basis points and the generator's spec entries on
        # every coupon; the held-out coupon's recorded coefficients are heldout_trace and its
        # fabricated matrices carry ring_perturbation on the metal-ring block.
        root = Path(tempfile.mkdtemp(dir=self.root))
        dirs = [write_coupon(root, "convex", a, 112.5, basis_points=BASIS_POINTS) for a in (90.0, 105.0, 120.0, 135.0)]
        dirs.append(write_coupon(root, "convex", heldout_angle, None, trace=heldout_trace, basis_points=BASIS_POINTS,
                                 ring_perturbation=ring_perturbation))
        return dirs

    def test_gate_definition(self):
        self.assertEqual(GATES["MaximumParticipationReferencedResidual"], 0.005)
        self.assertIn("participation-referenced", GATES["Statement"])
        self.assertIn("DefectReferenced", GATES["ReportedAlongside"])
        # USER decision 161 (2): the option-(c) traces gate the family.
        self.assertEqual(GATES["GatingTrace"], "option-c")

    def test_all_rings_family_is_one_segment_with_plain_node_specs(self):
        # An AllRingsFollowMetal family (no connectivity records, no events): the nodes are
        # the --nodes ANGLE list, the stencil the cubic sliding window; without --nodes the
        # nodes cannot be told from the held-out coupons (an error, not a silent verdict); a
        # run on the band trace is not the gating one (Gating false), the option-c run is.
        record_all = CHECK.generator.trace_basis_rule(16, None, CHECK.generator.REFINED_RULE)
        root = Path(tempfile.mkdtemp(dir=self.root))
        dirs = [write_coupon(root, "convex", a, None, trace_basis=record_all) for a in (75.0, 90.0, 105.0, 120.0, 135.0)]
        dirs.append(write_coupon(root, "convex", 82.5, None, trace_basis=record_all))
        code, record, _ = self.run_check(dirs, nodes=["75", "90", "105", "120", "135"])
        self.assertEqual(code, 0, record)
        entry = record["Families"]["convex"]["HeldOut"][0]
        self.assertEqual((entry["AngleDegrees"], entry["Rule"], entry["Stencil"]), (82.5, "cubic", [75.0, 90.0, 105.0, 120.0]))
        self.assertIsNone(entry["ConnectivityAngleDegrees"])
        self.assertEqual(record["GatingTrace"], "option-c")
        self.assertFalse(record["Gating"])  # the recorded (band) trace of these synthetic caches
        code, record, output = self.run_check(dirs)
        self.assertEqual(code, 1)
        self.assertIsNone(record)
        self.assertIn("pass --nodes ANGLE", output)

    def test_option_c_run_is_the_gating_one(self):
        heldout_trace = CHECK.generator.heldout_potential(BASIS_POINTS, RADIUS, THICKNESS, 112.0, "convex")
        code, record, _ = self.run_check(self.family_with_basis_points(heldout_trace), trace="option-c")
        self.assertEqual(code, 0, record)
        self.assertEqual(record["TraceForms"], ["option-c"])
        self.assertTrue(record["Gating"])
        code, record, _ = self.run_check(self.family_with_basis_points(heldout_trace), trace="recorded")
        self.assertEqual(code, 0, record)
        self.assertTrue(record["Gating"])  # the recorded coefficients ARE the option-(c) ones here
        band = CHECK.generator.metal_band_cutoff(BASIS_POINTS, RADIUS, THICKNESS) * CHECK.generator.heldout_polynomial(BASIS_POINTS, RADIUS)
        code, record, _ = self.run_check(self.family_with_basis_points(band), trace="recorded")
        self.assertEqual(record["TraceForms"], ["band"])
        self.assertFalse(record["Gating"])
        # The band trace recomputed on caches that record option-(c) coefficients (the verdict
        # reported alongside the gating one).
        code, record, _ = self.run_check(self.family_with_basis_points(heldout_trace), trace="band")
        self.assertEqual(code, 0, record)
        self.assertEqual((record["TraceSource"], record["TraceForms"], record["Gating"]), ("band", ["band"], False))

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

    def test_recorded_trace_form_is_classified_and_recorded(self):
        # A recorded coefficient file equal to the generator's band-cutoff trace at the basis
        # points is stamped "band", the option-(c) trace "option-c", anything else
        # "recorded-unclassified"; a coupon without basis points is unclassified too.
        band = GENERATOR.metal_band_cutoff(BASIS_POINTS, RADIUS, THICKNESS) * GENERATOR.heldout_polynomial(BASIS_POINTS, RADIUS)
        option_c = GENERATOR.heldout_potential(BASIS_POINTS, RADIUS, THICKNESS, 112.0, "convex")
        self.assertEqual(list(band[:2]), [0.0, 0.0])  # blind on both metal rings
        self.assertTrue(np.all(option_c[:2] != 0.0))  # the free metal-ring knots excited
        for trace, form in ((band, "band"), (option_c, "option-c"), (np.array([0.3, -0.2, 0.5, 0.1]), "recorded-unclassified")):
            code, record, stdout = self.run_check(self.family_with_basis_points(trace))
            self.assertEqual(code, 0, stdout)
            family = record["Families"]["convex"]
            self.assertEqual(record["TraceSource"], "recorded")
            self.assertEqual(record["TraceForms"], [form])
            self.assertEqual(family["TraceForms"], [form])
            self.assertEqual(family["HeldOut"][0]["TraceForm"], form)
            self.assertIn(f"held-out trace form {form}", stdout)
            self.assertIn("on the recorded held-out trace", stdout)
        code, record, _ = self.run_check(self.segment(scale=1.004))
        self.assertEqual(record["TraceForms"], ["recorded-unclassified"])
        with open(self.root / "out.csv") as handle:
            rows = list(CHECK.csv.reader(handle))
        self.assertEqual(rows[0][-1], "trace_form")
        self.assertEqual({r[-1] for r in rows[1:]}, {"recorded-unclassified"})

    def test_option_c_trace_is_recomputed_from_the_basis_points(self):
        # The held-out coupon's fabricated matrices are perturbed on the metal-ring block only
        # (the review's MAJOR-1 mechanism): on the recorded band trace, zero on both metal
        # rings, the family reproduces the coupon exactly and PASSES; --trace option-c ignores
        # the recorded coefficients, recomputes the option-(c) trace from the basis points and
        # sees the perturbation: FAIL, with the residuals of a cache whose recorded file IS the
        # option-(c) trace.
        band = GENERATOR.metal_band_cutoff(BASIS_POINTS, RADIUS, THICKNESS) * GENERATOR.heldout_polynomial(BASIS_POINTS, RADIUS)
        option_c = GENERATOR.heldout_potential(BASIS_POINTS, RADIUS, THICKNESS, 112.0, "convex")
        code_band, record_band, _ = self.run_check(self.family_with_basis_points(band, ring_perturbation=0.1))
        code_c, record_c, stdout = self.run_check(self.family_with_basis_points(band, ring_perturbation=0.1), trace="option-c")
        code_ref, record_ref, _ = self.run_check(self.family_with_basis_points(option_c, ring_perturbation=0.1))
        residuals = lambda r: r["Families"]["convex"]["HeldOut"][0]["ParticipationReferencedPercent"]  # noqa: E731
        self.assertEqual(code_band, 0)
        self.assertEqual(record_band["TraceForms"], ["band"])
        for value in residuals(record_band).values():
            self.assertAlmostEqual(value, 0.0, places=9)
        self.assertEqual((code_c, code_ref), (1, 1), stdout)
        self.assertEqual(record_c["TraceSource"], "option-c")
        self.assertEqual(record_c["TraceForms"], ["option-c"])
        self.assertIn("on the option-c held-out trace", stdout)
        self.assertFalse(record_c["Passed"])
        for quantity, value in residuals(record_c).items():
            self.assertAlmostEqual(value, residuals(record_ref)[quantity], places=9)
        self.assertLess(residuals(record_c)["SA"], -0.5)
        # Without basis points the option-(c) trace cannot be recomputed: fail closed.
        code, record, stdout = self.run_check(self.segment(), trace="option-c")
        self.assertNotEqual(code, 0)
        self.assertIsNone(record)
        self.assertIn("basis-points.csv", stdout)

    def test_fail_closed_guards(self):
        # A held-out coefficient vector of the wrong length is an error, not a truncation.
        code, record, stdout = self.run_check(self.segment() + [write_coupon(
            self.root, "convex", 118.0, None, trace=np.array([0.3, -0.2, 0.5, 0.1, 0.7]))])
        self.assertNotEqual(code, 0)
        self.assertIsNone(record)
        self.assertIn("5 held-out coefficients for 4 basis functions", stdout)
        # A zero fabricated held-out energy leaves the residual undefined: an error.
        dirs = self.segment()
        np.savetxt(dirs[-1] / "heldout-coefficients.csv", np.zeros(SIZE), delimiter=",", header="coefficient_V", comments="")
        code, record, stdout = self.run_check(dirs)
        self.assertNotEqual(code, 0)
        self.assertIsNone(record)
        self.assertIn("zero held-out energy", stdout)


if __name__ == "__main__":
    unittest.main()
