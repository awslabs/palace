# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""spatial_qualification: the (F) qualification upgrade (decisions 282 / 285 (5), restated by
decision 299) on synthetic energies, matrices, traces and records - the criteria (closure, the
Domain twin-consistency, the identity), the thin surface p-steps as information, the
Palace-output readers and the twin class map, the dense trace family, the resolution-aware
gate, the reference box hook and the status transitions."""
import json
import math
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import spatial_qualification as sq  # noqa: E402


class CriteriaTest(unittest.TestCase):
    def test_closure_and_twin_consistency_residuals(self):
        thin = {"SA": 1.0, "MS": 2.0, "Domain": 10.0}
        de = {"SA": 0.10, "MS": -0.20, "Domain": 0.50}
        fab = {"SA": 1.10, "MS": 1.80, "Domain": 10.6}
        closure = sq.dense_closure(thin, de, fab)
        self.assertTrue(closure["SA"]["Passed"])
        self.assertTrue(closure["MS"]["Passed"])
        self.assertAlmostEqual(closure["Domain"]["Residual"], 0.1 / 10.6)
        self.assertTrue(closure["Domain"]["Passed"])
        # Twin-consistency judges the Domain class only (decision 299 (1a)): a 3.8 % domain
        # defect step fails, a 30 % MS step of the thin twin is not judged.
        twin = sq.domain_twin_consistency({"SA": 0.10, "MS": -0.8, "Domain": 0.9}, de, fab)
        self.assertEqual(sorted(twin), ["Domain"])
        self.assertFalse(twin["Domain"]["Passed"])  # |0.9 - 0.5| / 10.6 = 3.8 %
        self.assertAlmostEqual(twin["Domain"]["Residual"], 0.4 / 10.6)
        self.assertTrue(sq.domain_twin_consistency({"Domain": 0.6}, de, fab)["Domain"]["Passed"])
        with self.assertRaises(sq.SpatialQualificationError):
            sq.domain_twin_consistency({"SA": 0.1}, {"SA": 0.1}, {"SA": 1.0})
        # A zero reference energy with a nonzero residual is infinite, never silently passed.
        self.assertFalse(sq.dense_closure({"SA": 1.0}, {"SA": 0.0}, {"SA": 0.0})["SA"]["Passed"])

    def test_thin_surface_p_steps_are_information(self):
        steps = sq.thin_surface_p_steps({"MS": 2.0, "SA": 1.0, "Domain": 10.0}, {"MS": 2.3, "SA": 0.98, "Domain": 10.1},
                                        {"MS": 1.0, "SA": 1.1, "Domain": 10.6})
        self.assertEqual(sorted(steps), ["MS", "SA"])  # the Domain is a criterion, not a step record
        self.assertAlmostEqual(steps["MS"]["Step"], 0.3)
        self.assertAlmostEqual(steps["MS"]["RelativeToThin"], 0.3 / 2.3)
        self.assertAlmostEqual(steps["MS"]["RelativeToFabricated"], 0.3)
        self.assertAlmostEqual(steps["SA"]["RelativeToFabricated"], -0.02 / 1.1)
        for value in steps.values():
            self.assertNotIn("Passed", value)

    def test_quadratic_form_dense_and_sparse_agree(self):
        q = np.array([[2.0, 0.5, 0.0], [0.5, 1.0, -0.25], [0.0, -0.25, 3.0]])
        t = [1.0, -2.0, 0.5]
        dense = sq.quadratic_form(q, t)
        sparse = sq.quadratic_form({(1, 1): 2.0, (1, 2): 0.5, (2, 2): 1.0, (2, 3): -0.25, (3, 3): 3.0}, t)
        self.assertAlmostEqual(dense, float(np.array(t) @ q @ np.array(t)))
        self.assertAlmostEqual(sparse, dense)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.quadratic_form({(1, 4): 1.0}, t)

    def test_matrix_identity(self):
        q = {"SA": {(1, 1): 2.0, (1, 2): 0.5, (2, 2): 1.0}}
        t = [1.0, 3.0]
        exact = 2.0 + 2 * 0.5 * 3.0 + 9.0
        result = sq.matrix_identity({"SA": exact}, q, t)
        self.assertTrue(result["SA"]["Passed"])
        self.assertEqual(result["SA"]["Residual"], 0.0)
        off = sq.matrix_identity({"SA": exact * (1 + 1e-5)}, q, t)
        self.assertFalse(off["SA"]["Passed"])
        with self.assertRaises(sq.SpatialQualificationError):
            sq.matrix_identity({"MS": 1.0}, q, t)

    def test_gate(self):
        diagnostics = {"Count": 0, "Tolerance": 0.02,
                       "Records": [{"Model": "a", "MaxRatio": 7.7e-8, "Excluded": False},
                                   {"Model": "a", "MaxRatio": 2.0e-7, "Excluded": False},
                                   {"Model": "b", "MaxRatio": 5.0e-2, "Excluded": True}]}
        self.assertTrue(sq.conductor_consistency_gate(diagnostics, "a")["Passed"])
        self.assertAlmostEqual(sq.conductor_consistency_gate(diagnostics, "a")["MaxRatio"], 2.0e-7)
        self.assertFalse(sq.conductor_consistency_gate(diagnostics, "b")["Passed"])
        missing = sq.conductor_consistency_gate(diagnostics, "c")
        self.assertFalse(missing["Passed"])
        self.assertIn("untestable", missing["Reason"])
        failing = {"Count": 1, "Tolerance": 0.02, "Records": [{"Model": "a", "MaxRatio": 3e-2, "Excluded": True}]}
        self.assertFalse(sq.conductor_consistency_gate(failing, "a")["Passed"])

    def test_resolution_aware_gate(self):
        """Decision 299 (2): the A6 box-3 readings (2.2e-7 at c0 / Order 3, 6.3e-7 at Order 5, 5.4e-6 at c7
        / Order 5) all pass the 1e-4 acceptance (the c0-only 1e-6 of decision 282 failed from c2 on);
        fictitious metal (D3-A 0.20-0.43) fails; every supplied solve must pass; the production
        orders not supplied are named."""
        def solve(ratio, order, excluded=False):
            return {"Diagnostics": {"Count": 1 if excluded else 0, "Tolerance": 0.02,
                                    "Records": [{"Model": "m", "MaxRatio": ratio, "Excluded": excluded}]},
                    "Order": order, "Source": f"p{order}"}
        gate = sq.resolution_aware_gate([solve(2.2e-7, 3), solve(4.8e-7, 4), solve(6.3e-7, 5)], "m")
        self.assertTrue(gate["Passed"])
        self.assertEqual(gate["Orders"], [3, 4, 5])
        self.assertEqual(gate["ProductionOrdersMissing"], [])
        self.assertAlmostEqual(gate["MaxRatio"], 6.3e-7)
        self.assertEqual(gate["Tolerance"], 1e-4)
        refined = sq.resolution_aware_gate([solve(5.4e-6, 5)], "m")
        self.assertTrue(refined["Passed"])
        self.assertEqual(refined["ProductionOrdersMissing"], [4])
        self.assertFalse(sq.resolution_aware_gate([solve(2.2e-7, 4), solve(0.2, 5, excluded=True)], "m")["Passed"])
        self.assertFalse(sq.resolution_aware_gate([solve(2.0e-4, 4)], "m")["Passed"])
        self.assertFalse(sq.resolution_aware_gate([solve(1e-7, 4)], "other")["Passed"])
        with self.assertRaises(sq.SpatialQualificationError):
            sq.resolution_aware_gate([], "m")

    def test_reference_box_closure_marker_and_bracket(self):
        result = sq.reference_box_closure(1.03, 0.98, 0.04, validated_class=True)
        self.assertAlmostEqual(result["Reference"], 1.0)
        self.assertAlmostEqual(result["Ratio"], 1.03)
        self.assertEqual(result["Marker"], sq.REFERENCE_MARKER_VALIDATED)
        self.assertTrue(result["Passed"])
        self.assertAlmostEqual(result["Bracket"][0], 1.03 / 1.02)
        self.assertAlmostEqual(result["Bracket"][1], 1.03 / 0.98)
        self.assertFalse(sq.reference_box_closure(1.08, 1.0, 0.0, True)["Passed"])
        self.assertTrue(sq.reference_box_closure(1.08, 1.0, 0.0, False)["Passed"])  # new class: 0.10

    def test_status_transitions(self):
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=True, identity_passed=True, gate_passed=True),
                         sq.STATUS_QUALIFIED)
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=True, identity_passed=True, gate_passed=True,
                                                 window_validated=True), sq.STATUS_WINDOW_VALIDATED)
        for dense, identity, gate in ((False, True, True), (True, False, True), (True, True, False)):
            self.assertEqual(sq.qualification_status(sq.STATUS_QUALIFIED, dense_passed=dense, identity_passed=identity,
                                                     gate_passed=gate), sq.STATUS_FAILED)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.qualification_status("Accepted", dense_passed=True, identity_passed=True, gate_passed=True)
        # Decisions 474 / 477 (1): an unjudged p-sequence Type refuses the Qualified / WindowValidated
        # lift (PendingQualification stays), while a failing (a) / identity / gate is still Failed.
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=True, identity_passed=True, gate_passed=True,
                                                 unjudged_types=["p_SA"]), sq.STATUS_PENDING)
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=True, identity_passed=True, gate_passed=True,
                                                 window_validated=True, unjudged_types=["p_SA"]), sq.STATUS_PENDING)
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=False, identity_passed=True, gate_passed=True,
                                                 unjudged_types=["p_SA"]), sq.STATUS_FAILED)


def write_csv(path, header, rows):
    lines = [",".join(header)] + [",".join(f"{v:+.12e}" if isinstance(v, float) else str(v) for v in row) for row in rows]
    Path(path).write_text("\n".join(lines) + "\n")


class PalaceOutputReadersTest(unittest.TestCase):
    """The within-R energies of a dense-trace run and the p4 matrices from Palace's CSVs."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="spatial-qualification-"))
        self.classes = {1: "SA", 2: "MS", 3: "MA"}
        postpro = self.tmp / "postpro"
        postpro.mkdir()
        write_csv(postpro / "domain-E.csv", ["i", "E_elec (J)", "E_mag (J)"], [[1, 10.0, 0.0], [2, 20.0, 0.0]])
        write_csv(postpro / "surface-Q.csv", ["i", "p_surf[1]", "Q_surf[1]", "p_surf[2]", "Q_surf[2]", "p_surf[3]", "Q_surf[3]"],
                  [[1, 0.1, 1.0, 0.2, 1.0, 0.05, 1.0], [2, 0.3, 1.0, 0.1, 1.0, 0.02, 1.0]])
        rows = []
        for i in (1, 2):
            for interface, e_out in ((1, 0.25), (2, 0.5), (3, 0.0)):
                rows.append([i, 0, interface, 2e-6, e_out * i, 0.0, 0.1, 0.0])
                rows.append([i, 0, interface, 1e-6, e_out * i + 0.1, 0.0, 0.1, 0.0])  # a smaller radius: ignored
        write_csv(postpro / "surface-Q-edge.csv", ["i", "exc", "interface", "R (m)", "E_out (J)", "p_out", "E_ann (J)", "p_ann"], rows)
        self.postpro = postpro

    def test_within_r_energies(self):
        energies = sq.within_r_energies(self.postpro, self.classes)
        self.assertEqual(sorted(energies), [1, 2])
        self.assertAlmostEqual(energies[1]["Domain"], 10.0)
        self.assertAlmostEqual(energies[1]["SA"], 0.1 * 10.0 - 0.25)   # p_surf E_elec - E_out at the largest R
        self.assertAlmostEqual(energies[1]["MS"], 0.2 * 10.0 - 0.5)
        self.assertAlmostEqual(energies[2]["MA"], 0.02 * 20.0 - 0.0)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.within_r_energies(self.postpro, {1: "SA", 2: "MS"})

    def test_response_matrices(self):
        write_csv(self.postpro / "domain-response-matrix.csv", ["basis_i", "basis_j", "Q_ij (J)"],
                  [[1, 1, 2.0], [1, 2, 0.5], [2, 2, 1.0]])
        rows = []
        for interface in (1, 2, 3):
            for (i, j, q) in ((1, 1, 1.0 * interface), (1, 2, 0.1), (2, 2, 0.5)):
                rows.append([interface, 0, 2e-6, i, j, q, 0.0, 0.0, q, 0.0, 0.0])
                rows.append([interface, 7, 2e-6, i, j, 99.0, 0.0, 0.0, 99.0, 0.0, 0.0])  # a per-edge group: ignored
        write_csv(self.postpro / "surface-response-matrix.csv",
                  ["interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)",
                   "Q_total_ij (J)", "Q_total_ij normal (J)", "Q_total_ij tangential (J)"], rows)
        matrices = sq.response_matrices(self.postpro, self.classes)
        self.assertEqual(matrices["Domain"][(1, 2)], 0.5)
        self.assertEqual(matrices["MS"][(1, 1)], 2.0)
        self.assertEqual(matrices["SA"][(2, 2)], 0.5)
        difference = sq.matrix_difference(matrices, {cls: {k: 0.5 * v for k, v in m.items()} for cls, m in matrices.items()})
        self.assertEqual(difference["MA"][(1, 1)], 1.5)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.matrix_difference(matrices, {"Domain": {}})


    def test_twin_interface_classes_union(self):
        """Decision 293 (2): the fabricated twin (2 MS, 3 SA, 4.. the MA shells) and the thin twin (1 MA,
        2 MS, 3 SA) carry different indices; the fabricated map alone cannot read the thin runs
        (interface 1 has no class), the union reads both; a conflicting index fails closed."""
        def config(entries):
            return {"Boundaries": {"Postprocessing": {"Dielectric": [{"Index": i, "Type": t} for i, t in entries]}}}
        fabricated = config([(2, "MS"), (3, "SA"), (4, "MA"), (5, "MA")])
        thin = config([(1, "MA"), (2, "MS"), (3, "SA")])
        with self.assertRaises(sq.SpatialQualificationError):
            sq.within_r_energies(self.postpro, sq.interface_classes(fabricated))  # interface 1 of the thin run: no class
        classes = sq.twin_interface_classes({"fabricated-p4": fabricated, "fabricated-p5": fabricated,
                                             "thin-p4": thin, "thin-p5": thin})
        self.assertEqual(classes, {1: "MA", 2: "MS", 3: "SA", 4: "MA", 5: "MA"})
        energies = sq.within_r_energies(self.postpro, classes)
        self.assertAlmostEqual(energies[1]["MA"], 0.1 * 10.0 - 0.25)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.twin_interface_classes({"fabricated-p4": fabricated, "thin-p4": config([(2, "SA")])})

class DenseTracesTest(unittest.TestCase):
    """The T2 family on a tiny synthetic basis, the T1 reader on the production trace format,
    the representable excitation and the trace file / config writers."""

    def setUp(self):
        # A box [-2, 2] x [-1, 1] x [-0.5, 0.5] (R = 1) with 4 basis knots on the top cap rim,
        # two conductor-labelled vertices (conductors 1 and 2) and a slave vertex.
        points = np.array([[-2.0, -1.0, 0.5], [2.0, -1.0, 0.5], [2.0, 1.0, 0.5], [-2.0, 1.0, 0.5],
                           [0.0, -1.0, 0.0], [0.0, 1.0, 0.0], [-2.0, 0.0, 0.0]])
        self.basis = {"Points": points, "Triangles": np.array([[0, 1, 4], [1, 2, 5], [2, 3, 6], [3, 0, 6]]),
                      "Basis": np.array([1, 2, 3, 4, 0, 0, 0]),
                      "Lower": np.array([-2.0, -1.0, -0.5]), "Upper": np.array([2.0, 1.0, 0.5]), "Frame": np.identity(3)}
        self.labels = np.array([0, 0, 0, 0, 1, 2, 0])

    def test_synthetic_family(self):
        traces = sq.synthetic_traces(self.basis, self.labels, 1.0, {1: np.zeros(3), 2: np.array([0.0, 1.0, -1.0])})
        self.assertEqual([t["Name"] for t in traces], ["state-2", "line-charge-x0", "line-charge-x1", "line-charge-y0",
                                                       "line-charge-y1"])
        self.assertEqual(traces[0]["Coefficients"], [0.0, 0.0, 0.0, 0.0, 1.0])
        for trace in traces[1:]:
            self.assertEqual(len(trace["Coefficients"]), 5)
            self.assertEqual(trace["Coefficients"][-1], 0.0)  # grounded conductors
            knots = trace["Coefficients"][:4]
            self.assertAlmostEqual(min(knots), 0.0)
            self.assertAlmostEqual(max(knots), 1.0)
        # The x0 line charge sits at x = -7: the knots at x = -2 are closer (higher potential).
        x0 = traces[1]["Coefficients"]
        self.assertGreater(x0[0], x0[1])
        self.assertGreater(x0[3], x0[2])

    def single_conductor_basis(self):
        """A single-conductor basis whose in-metal vertex IS a basis knot (the ZeroTrace
        construction of generate_spatial_response for conductor_count == 1): the 4 rim knots
        plus knot 5 at (-2, 0, 0) inside conductor 1."""
        basis = dict(self.basis)
        basis["Basis"] = np.array([1, 2, 3, 4, 0, 0, 5])
        return basis, np.array([0, 0, 0, 0, 1, 1, 1])

    def test_synthetic_family_vanishes_at_pec_knots(self):
        """Decision 360 (b*): the line-charge value is 0 at every knot inside a conductor; the
        free knots keep the values of the unit-range normalisation over ALL knots (so a
        two-conductor coupon, whose PEC knots have no column, is unchanged byte for byte)."""
        basis, labels = self.single_conductor_basis()
        traces = sq.synthetic_traces(basis, labels, 1.0, {1: np.zeros(3)})
        self.assertEqual([t["Name"] for t in traces], ["line-charge-x0", "line-charge-x1", "line-charge-y0",
                                                       "line-charge-y1"])
        knots = np.array([[-2.0, -1.0], [2.0, -1.0], [2.0, 1.0], [-2.0, 1.0], [-2.0, 0.0]])
        for trace, source in zip(traces, ([-7.0, 0.0], [7.0, 0.0], [0.0, -6.0], [0.0, 6.0])):
            self.assertEqual(len(trace["Coefficients"]), 5)
            self.assertEqual(trace["Coefficients"][4], 0.0)  # the PEC knot: exactly zero
            phi = -np.log(np.linalg.norm(knots - np.asarray(source), axis=1))
            expected = (phi - phi.min()) / np.ptp(phi)  # normalised over all 5 knots, PEC knot included
            np.testing.assert_allclose(trace["Coefficients"][:4], expected[:4], rtol=0, atol=1e-15)
        # The x0 line charge is nearest to the PEC knot (-2, 0): without the fix that knot would
        # carry the maximum 1.0; with it the free rim knots are unchanged.
        self.assertLess(max(traces[0]["Coefficients"][:4]), 1.0)

    def test_write_dense_trace_fails_closed_at_pec_knots(self):
        """Decision 360 (b*): a hat prescribing a non-zero potential at a conductor vertex is a
        Dirichlet jump along the PEC cross-section boundary and is refused; zero there (the
        synthetic family, the runtime-zeroed device trace) is written."""
        basis, labels = self.single_conductor_basis()
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(sq.SpatialQualificationError) as caught:
                sq.write_dense_trace(Path(tmp) / "jump.csv", basis, labels, [0.1, 0.2, 0.3, 0.4, 0.5])
            self.assertIn("PEC", str(caught.exception))
            self.assertIsInstance(caught.exception, sq.DenseTraceJumpError)
            self.assertFalse((Path(tmp) / "jump.csv").exists())
            path, terminal, scale = sq.write_dense_trace(Path(tmp) / "ok.csv", basis, labels, [0.1, 0.2, 0.3, 0.4, 0.0])
            self.assertEqual((terminal, scale), (None, 1.0))
            rows = [line.split(",") for line in Path(path).read_text().splitlines() if line and not line.startswith("x")]
            values = {tuple(round(float(v), 9) for v in row[:3]): float(row[3]) for row in rows}
            self.assertEqual(values[(-2.0, 0.0, 0.0)], 0.0)
            self.assertAlmostEqual(values[(-2.0, -1.0, 0.5)], 0.1)
            for trace in sq.synthetic_traces(basis, labels, 1.0, {1: np.zeros(3)}):
                sq.write_dense_trace(Path(tmp) / f"{trace['Name']}.csv", basis, labels, trace["Coefficients"])

    @staticmethod
    def device_traces_csv(path, patch, coefficients, states):
        """A surface-response-traces.csv of one placed patch (model 1, first excitation): the
        basis coefficients 1..N then the conductor states {conductor: value}."""
        rows = [[1, patch, 1, k, 0, value] for k, value in enumerate(coefficients, start=1)]
        rows += [[1, patch, 1, 0, conductor, value] for conductor, value in sorted(states.items())]
        path.write_text("i,patch,model,coefficient,conductor state,value (V)\n" +
                        "\n".join(",".join(f"{v:.6e}" for v in row) for row in rows) + "\n")

    def traces_arguments(self, tmp, device_traces):
        return types.SimpleNamespace(source=str(Path(tmp) / "source"), output=str(Path(tmp) / "out"),
                                     device_traces=str(device_traces), model_index=1, fabricated_mesh=None,
                                     thin_mesh=None, run_config=None, thin_run_config=None, orders=[4, 5])

    def test_command_traces_aborts_on_pec_jump(self):
        """Review MAJOR-1 of decision 369: a T1 device trace with a non-zero potential at a
        conductor knot ABORTS the traces command (DenseTraceJumpError is not the recorded
        Unsupported limitation): no dense-traces.json, so the (F) cannot stamp from T2 alone."""
        basis, labels = self.single_conductor_basis()
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.object(sq, "load_basis", return_value=(basis, labels, {"Name": "m"}, 1.0, {1: np.zeros(3)})):
            device = Path(tmp) / "surface-response-traces.csv"
            self.device_traces_csv(device, 401, [0.1, -0.1, 0.2, -0.7, -0.7025], {})
            with self.assertRaises(sq.DenseTraceJumpError) as caught:
                sq.command_traces(self.traces_arguments(tmp, device))
            self.assertIn("PEC", str(caught.exception))
            self.assertFalse((Path(tmp) / "out" / "dense-traces.json").exists())
            # The runtime-zeroed T1 (0 at the PEC knot) is written with the T2 family.
            self.device_traces_csv(device, 401, [0.1, -0.1, 0.2, -0.7, 0.0], {})
            self.assertEqual(sq.command_traces(self.traces_arguments(tmp, device)), 0)
            record = json.loads((Path(tmp) / "out" / "dense-traces.json").read_text())
            self.assertEqual([t["Name"] for t in record["Traces"]],
                             ["line-charge-x0", "line-charge-x1", "line-charge-y0", "line-charge-y1", "device-patch-401"])
            self.assertEqual(record["Unsupported"], [])

    def test_command_traces_records_unrepresentable_trace(self):
        """The pre-existing limitation stays a record, not an abort: a T1 trace with two non-zero
        conductor states is written to Unsupported and the command completes."""
        labels = np.array([0, 0, 0, 0, 1, 2, 3])  # conductors 2 and 3: two state lifts
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.object(sq, "load_basis", return_value=(self.basis, labels, {"Name": "m"}, 1.0,
                                                                  {1: np.zeros(3), 2: np.array([0.0, 1.0, -1.0]),
                                                                   3: np.array([-2.0, 0.0, 0.0])})):
            device = Path(tmp) / "surface-response-traces.csv"
            self.device_traces_csv(device, 401, [0.1, 0.2, 0.3, 0.4], {2: 1.0, 3: 2.0})
            self.assertEqual(sq.command_traces(self.traces_arguments(tmp, device)), 0)
            record = json.loads((Path(tmp) / "out" / "dense-traces.json").read_text())
            self.assertEqual([t["Name"] for t in record["Traces"]],
                             ["state-2", "state-3", "line-charge-x0", "line-charge-x1", "line-charge-y0", "line-charge-y1"])
            self.assertEqual([t["Name"] for t in record["Unsupported"]], ["device-patch-401"])
            self.assertIn("TerminalAttributes hold one volt", record["Unsupported"][0]["Unsupported"])
            self.assertFalse((Path(tmp) / "out" / "traces" / "device-patch-401.csv").exists())

    def test_representable_excitation(self):
        self.assertEqual(sq.representable_trace([0.5, 0.25, 0.0, 0.0], [2, 3]), ([0.5, 0.25, 0.0, 0.0], None, 1.0))
        scaled, terminal, scale = sq.representable_trace([1.0, 2.0, 4.0, 0.0], [2, 3])
        self.assertEqual((terminal, scale), (2, 16.0))
        self.assertEqual(scaled, [0.25, 0.5, 1.0, 0.0])
        with self.assertRaises(sq.SpatialQualificationError):
            sq.representable_trace([1.0, 2.0, 4.0, -1.0], [2, 3])

    def test_write_dense_trace_and_config(self):
        with tempfile.TemporaryDirectory() as tmp:
            path, terminal, scale = sq.write_dense_trace(Path(tmp) / "trace.csv", self.basis, self.labels,
                                                         [0.1, 0.2, 0.3, 0.4, 2.0])
            self.assertEqual((terminal, scale), (2, 4.0))
            rows = [line.split(",") for line in Path(path).read_text().splitlines() if line and not line.startswith("x")]
            values = {tuple(round(float(v), 9) for v in row[:3]): float(row[3]) for row in rows}  # per vertex (x, y, z)
            self.assertAlmostEqual(values[(-2.0, -1.0, 0.5)], 0.05)  # knot 1 scaled by 1 / 2
            self.assertAlmostEqual(values[(0.0, -1.0, 0.0)], 0.0)    # conductor 1 (the reference)
            self.assertAlmostEqual(values[(0.0, 1.0, 0.0)], 1.0)     # conductor 2: the terminal at one volt
            self.assertAlmostEqual(values[(-2.0, 0.0, 0.0)], 0.0)    # a slave vertex (no knot, no conductor)
            with self.assertRaises(sq.SpatialQualificationError):
                sq.write_dense_trace(Path(tmp) / "bad.csv", self.basis, self.labels, [0.1, 0.2])
            reference = {"Model": {"Mesh": "fab.msh"}, "Problem": {"Output": "x"},
                         "Solver": {"Order": 4, "Electrostatic": {"ResponseMatrix": True}},
                         "Boundaries": {"PrescribedPotential": [{"Index": 1, "Attributes": [1], "DataFile": "basis-0001.csv"},
                                                                {"Index": 2, "Attributes": [1], "TerminalAttributes": [5002, 6002],
                                                                 "DataFile": "conductor-2.csv"}]}}
            config = sq.dense_trace_config(reference, "thin.msh", "out/postpro",
                                           [{"Path": "a.csv", "Terminal": None}, {"Path": "b.csv", "Terminal": 2}], 5)
            self.assertEqual(config["Solver"]["Order"], 5)
            self.assertNotIn("ResponseMatrix", config["Solver"]["Electrostatic"])
            sources = config["Boundaries"]["PrescribedPotential"]
            self.assertEqual([s["Index"] for s in sources], [1, 2])
            self.assertNotIn("TerminalAttributes", sources[0])
            self.assertEqual(sources[1]["TerminalAttributes"], [5002, 6002])
            with self.assertRaises(sq.SpatialQualificationError):
                sq.dense_trace_config(reference, "thin.msh", "out", [{"Path": "c.csv", "Terminal": 3}], 5)

    def test_device_traces_reader(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "surface-response-traces.csv"
            rows = []
            for patch, model in ((409, 2), (410, 3)):
                for k, value in enumerate((0.5, 0.25, -0.1), start=1):
                    rows.append([1, patch, model, k, 0, value * patch])
                rows.append([1, patch, model, 4, 1, 7.0])
                rows.append([2, patch, model, 1, 0, 99.0])  # the second excitation: ignored
            path.write_text("i,patch,model,coefficient,conductor state,value (V)\n" +
                            "\n".join(",".join(f"{v:.6e}" for v in row) for row in rows) + "\n")
            traces = sq.device_traces(path, 2)
            self.assertEqual(sorted(traces), [409])
            np.testing.assert_allclose(traces[409], [0.5 * 409, 0.25 * 409, -0.1 * 409, 7.0])
            self.assertEqual(sorted(sq.device_traces(path, 3)), [410])


class EvaluationTest(unittest.TestCase):
    """evaluate_trace / evaluate on consistent synthetic energies: a trace whose energies obey
    the model exactly qualifies; a trace with a 3 % closure defect fails the coupon; the gate
    and the window reference drive the status."""

    def setUp(self):
        self.q_fab = {"SA": {(1, 1): 2.0, (1, 2): 0.5, (2, 2): 1.0}, "Domain": {(1, 1): 10.0, (1, 2): 0.0, (2, 2): 5.0}}
        self.q_thin = {"SA": {(1, 1): 1.8, (1, 2): 0.5, (2, 2): 0.9}, "Domain": {(1, 1): 9.5, (1, 2): 0.0, (2, 2): 4.9}}
        self.t = [1.0, 2.0]
        self.fab_p4 = {cls: sq.quadratic_form(self.q_fab[cls], self.t) for cls in self.q_fab}
        self.thin_p4 = {cls: sq.quadratic_form(self.q_thin[cls], self.t) for cls in self.q_thin}
        # p5 within 0.5 % of p4 on the fabricated twin and on the thin domain: closure and the Domain
        # twin-consistency hold; the thin SA steps 6 % (the decision-66 log-divergence): information.
        self.fab_p5 = {cls: v * 1.005 for cls, v in self.fab_p4.items()}
        self.thin_p5 = {"SA": self.thin_p4["SA"] * 1.06, "Domain": self.thin_p4["Domain"] * 1.004}
        self.gate = {"Passed": True, "Probed": 2, "Count": 0, "MaxRatio": 7.7e-8, "Orders": [4, 5],
                     "ProductionOrdersMissing": []}

    def test_consistent_trace_qualifies(self):
        trace = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                  thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertTrue(trace["Passed"])
        self.assertAlmostEqual(trace["Energies"]["ModelCorrection"]["SA"], self.fab_p4["SA"] - self.thin_p4["SA"])
        self.assertEqual(trace["MatrixIdentity"]["SA"]["Residual"], 0.0)
        self.assertEqual(sorted(trace["TwinConsistency"]), ["Domain"])
        self.assertAlmostEqual(trace["ThinSurfacePSteps"]["SA"]["RelativeToThin"], 0.06 / 1.06)
        record = sq.evaluate([trace], gate=self.gate)
        self.assertEqual(record["Status"], sq.STATUS_QUALIFIED)
        self.assertTrue(record["TwinConsistencyPassed"] and record["ClosurePassed"])
        self.assertEqual(record["Tolerances"]["GateMaxRatio"], 1e-4)
        self.assertEqual(record["Families"], ["T2"])
        window = [sq.reference_box_closure(1.02, 1.0, 0.0, True) | {"Class": "SA", "Window": "S1p"}]
        self.assertEqual(sq.evaluate([trace], gate=self.gate, reference_boxes=window)["Status"], sq.STATUS_WINDOW_VALIDATED)
        off_marker = [sq.reference_box_closure(1.2, 1.0, 0.0, True)]
        self.assertEqual(sq.evaluate([trace], gate=self.gate, reference_boxes=off_marker)["Status"], sq.STATUS_QUALIFIED)
        self.assertEqual(sq.evaluate([trace], gate={"Passed": False, "Probed": 0})["Status"], sq.STATUS_FAILED)
        # The model's Qualification.UnjudgedTypes (decisions 474 / 477 (1)) refuses the lift with the reason recorded.
        unjudged = sq.evaluate([trace], gate=self.gate, reference_boxes=window, unjudged_types=["p_SA"])
        self.assertEqual(unjudged["Status"], sq.STATUS_PENDING)
        self.assertEqual(unjudged["UnjudgedTypes"], ["p_SA"])
        self.assertEqual(unjudged["UnjudgedTypesRule"], sq.UNJUDGED_TYPES_RULE)
        self.assertTrue(unjudged["DensePassed"] and unjudged["GatePassed"])
        self.assertIsNone(record["UnjudgedTypesRule"])
        self.assertEqual(record["UnjudgedTypes"], [])

    def test_defective_trace_fails_the_coupon(self):
        defective_p5 = dict(self.fab_p5)
        defective_p5["SA"] = self.fab_p4["SA"] * 1.03  # a 3 % p-convergence defect of the correction on SA
        trace = sq.evaluate_trace("line-charge-x0", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4,
                                  fab_p5=defective_p5, thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertFalse(trace["Closure"]["SA"]["Passed"])
        self.assertTrue(trace["Closure"]["Domain"]["Passed"])
        self.assertFalse(trace["Passed"])
        good = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                 thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        record = sq.evaluate([good, trace], gate=self.gate)
        self.assertFalse(record["DensePassed"])
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        identity_defect = sq.evaluate_trace("state-2", "T2", self.t, fab_p4={cls: v * (1 + 1e-5) for cls, v in self.fab_p4.items()},
                                            thin_p4=self.thin_p4, fab_p5=self.fab_p5, thin_p5=self.thin_p5,
                                            q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertFalse(identity_defect["MatrixIdentity"]["SA"]["Passed"])
        self.assertEqual(sq.evaluate([identity_defect], gate=self.gate)["Status"], sq.STATUS_FAILED)
        # The thin twin's domain step of 3 % (a p-unstable domain defect) fails the twin-consistency.
        unstable_thin_p5 = {"SA": self.thin_p5["SA"], "Domain": self.thin_p4["Domain"] * 0.97}
        unstable = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                     thin_p5=unstable_thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertFalse(unstable["TwinConsistency"]["Domain"]["Passed"])
        self.assertTrue(all(v["Passed"] for v in unstable["Closure"].values()))
        record = sq.evaluate([unstable], gate=self.gate)
        self.assertFalse(record["TwinConsistencyPassed"])
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.evaluate([], gate=self.gate)

    def test_stamp_library_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            path.write_text(json.dumps({"Models": [{"Name": "m", "LibraryQualified": False}]}))
            trace = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                      thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
            record = sq.evaluate([trace], gate=self.gate)
            library = sq.stamp_library_status(path, "m", record)
            model = library["Models"][0]
            self.assertEqual(model["QualificationStatus"], sq.STATUS_QUALIFIED)
            self.assertTrue(model["LibraryQualified"])
            self.assertEqual(model["SpatialQualification"]["Families"], ["T2"])
            self.assertTrue(model["SpatialQualification"]["TwinConsistencyPassed"])
            self.assertEqual(model["SpatialQualification"]["Gate"]["Orders"], [4, 5])
            self.assertEqual(model["SpatialQualification"]["Traces"][0]["Name"], "state-2")
            with self.assertRaises(sq.SpatialQualificationError):
                sq.stamp_library_status(path, "other", record)
            # An unjudged Type stamps PendingQualification, LibraryQualified stays False (decision 477 (1)).
            unjudged = sq.evaluate([trace], gate=self.gate, unjudged_types=["p_MS"])
            model = sq.stamp_library_status(path, "m", unjudged)["Models"][0]
            self.assertEqual(model["QualificationStatus"], sq.STATUS_PENDING)
            self.assertFalse(model["LibraryQualified"])
            self.assertEqual(model["SpatialQualification"]["UnjudgedTypes"], ["p_MS"])


if __name__ == "__main__":
    unittest.main()
