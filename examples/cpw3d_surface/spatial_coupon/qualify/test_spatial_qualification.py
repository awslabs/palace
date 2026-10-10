# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""spatial_qualification: the (F) qualification upgrade (decisions 282 / 285 (5), restated by
decision 299) on synthetic energies, matrices, traces and records - the criteria (closure, the
Domain twin-consistency, the identity), the thin surface p-steps as information, the
Palace-output readers and the twin class map, the dense trace family, the resolution-aware
gate, the reference box hook and the status transitions."""
import json
import math
import shutil
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

    def test_exclusion_coverage_gate(self):
        """Decision 657 (b): the EXCLUSION-COVERAGE gate on the palace.json record of a coupon run. The
        32b0 shape (its 23 foreign segments, coverage ~0.74 under the 1.9-nm rule) FAILS by name; the
        fixed shape (>= 0.99) passes; a missing record, a missing interface and an over-exclusion
        (> 1.01) fail closed; a config without EdgeExcludeSegments is NotApplicable and passes.
        Decision 665: the 0.99 floor is also met by the end-allowance coverage (Excluded + Straddling)
        / Listed - a short, many-ended straight-foreign model with a CORRECT exclusion (raw 0.918, 14
        ends inside ~50-nm perimeter segments) passes, the census's 0.9929 (7d88) and 1.0000 (448693)
        pass, 32b0 on main's rule (raw 0.7764; decision 667 (1): its chain's 2 ends lie on the box, 0
        straddlers, 21 removed runs on 1 chain -> 0.7764; 665's run-neighbour form read ~0.84) still
        fails, the upper bound stays on the raw coverage, a record without StraddlingPerimeterLength (a
        binary before decision 665) is read on the raw coverage alone, and the 7 patterns of the 667
        pre-registration (REPORT 12.1) read as pre-registered on the record's lengths."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            segments = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]] * 23
            config = {"Boundaries": {"Postprocessing": {"Dielectric": [
                {"Index": 1, "Type": "SA", "EdgeExcludeSegments": segments},
                {"Index": 2, "Type": "MS", "EdgeExcludeSegments": segments},
                {"Index": 3, "Type": "MA"}]}}}

            def write(run, coverages, listed=32.852, straddling=None, runs=None, chains=None):
                postpro = tmp / run
                postpro.mkdir(exist_ok=True)
                entries = {str(index): {"Coverage": value, "ExcludeSegmentsLength": listed, "ExcludedPerimeterLength": listed * value,
                                        "ExcludeSegmentsCount": 23, "ExcludedPerimeterCount": 400, "Tolerance": 1.9e-3,
                                        "MaxSagitta": 4.24e-3} for index, value in coverages.items()}
                if straddling is not None:
                    for entry in entries.values():
                        entry["StraddlingPerimeterLength"] = straddling
                        entry["EndAllowanceCoverage"] = (entry["ExcludedPerimeterLength"] + straddling) / listed
                        if runs is not None:
                            entry["RemovedRuns"] = runs
                            entry["ListedChains"] = chains
                (postpro / "palace.json").write_text(json.dumps({"SurfaceResponse": {"Diagnostics": {"EdgeExclusion": {
                    "Interfaces": entries, "Rule": "r"}}}}))
                return postpro

            defective = sq.exclusion_coverage_run("fabricated-p4", config, write("fabricated-p4", {"1": 0.737, "2": 0.737}))
            self.assertFalse(defective["Passed"])
            self.assertEqual(defective["Required"], [1, 2])
            self.assertIn("EXCLUSION-COVERAGE FAIL", defective["Reason"])
            self.assertAlmostEqual(defective["MinCoverage"], 0.737)
            self.assertAlmostEqual(defective["MinEndAllowanceCoverage"], 0.737)  # no straddling field: the raw coverage
            fixed = sq.exclusion_coverage_run("thin-p4", config, write("thin-p4", {"1": 0.9995, "2": 1.0003}))
            self.assertTrue(fixed["Passed"])
            self.assertIsNone(fixed["Reason"])
            self.assertEqual(sorted(fixed["Interfaces"]), ["1", "2"])
            over = sq.exclusion_coverage_run("fabricated-p5", config, write("fabricated-p5", {"1": 1.0, "2": 1.02}))
            self.assertFalse(over["Passed"])
            incomplete = sq.exclusion_coverage_run("thin-p5", config, write("thin-p5", {"1": 1.0}))
            self.assertFalse(incomplete["Passed"])
            self.assertIn("incomplete", incomplete["Reason"])
            (tmp / "legacy").mkdir()
            (tmp / "legacy" / "palace.json").write_text(json.dumps({"Problem": {"MeshElements": 1}}))
            missing = sq.exclusion_coverage_run("fabricated-matrices", config, tmp / "legacy")
            self.assertFalse(missing["Passed"])
            self.assertIn("record missing", missing["Reason"])
            not_applicable = sq.exclusion_coverage_run("thin-matrices", {"Boundaries": {"Postprocessing": {"Dielectric": [
                {"Index": 1, "Type": "SA"}]}}}, tmp / "legacy")
            self.assertTrue(not_applicable["Passed"])
            self.assertEqual(not_applicable["Required"], [])
            gate = sq.exclusion_coverage_gate([defective, fixed, over, incomplete, missing, not_applicable])
            self.assertFalse(gate["Passed"])
            self.assertTrue(gate["Applicable"])
            self.assertEqual(sorted(gate["Failures"]), ["fabricated-matrices", "fabricated-p4", "fabricated-p5", "thin-p5"])
            self.assertAlmostEqual(gate["MinCoverage"], 0.737)
            self.assertAlmostEqual(gate["MaxCoverage"], 1.02)
            self.assertEqual(gate["Rule"], sq.EXCLUSION_COVERAGE_RULE)
            passing = sq.exclusion_coverage_gate([fixed, not_applicable])
            self.assertTrue(passing["Passed"] and passing["Applicable"])
            empty = sq.exclusion_coverage_gate([not_applicable])
            self.assertTrue(empty["Passed"])
            self.assertFalse(empty["Applicable"])
            # Decision 665, the end allowance. A short, many-ended straight-foreign model: 7 pieces of
            # 6.0 um in total, 14 ends each keeping a ~50-nm perimeter segment with ~35 nm of deficit:
            # raw 0.918 (< 0.99 with a CORRECT exclusion), straddling 0.7 um -> 1.035: PASS.
            short = sq.exclusion_coverage_run("short-p4", config, write("short-p4", {"1": 0.918, "2": 0.918}, listed=6.0,
                                                                        straddling=0.7))
            self.assertTrue(short["Passed"])
            self.assertAlmostEqual(short["MinCoverage"], 0.918)
            self.assertAlmostEqual(short["MinEndAllowanceCoverage"], (6.0 * 0.918 + 0.7) / 6.0)
            self.assertAlmostEqual(short["Interfaces"]["1"]["EndAllowanceCoverage"], (6.0 * 0.918 + 0.7) / 6.0)
            self.assertEqual(short["Interfaces"]["1"]["StraddlingPerimeterLength"], 0.7)
            # The census readings: 7d88 thin (0.9929, 14 ends) and 448693 (1.0000) pass.
            census_min = sq.exclusion_coverage_run("7d88-thin-p4", config, write("7d88-thin-p4", {"1": 0.9929, "2": 0.9929}, listed=12.376,
                                                                                straddling=0.7))
            self.assertTrue(census_min["Passed"])
            arc_foreign = sq.exclusion_coverage_run("448693-p4", config, write("448693-p4", {"1": 1.0, "2": 1.0}, listed=17.855,
                                                                              straddling=0.1))
            self.assertTrue(arc_foreign["Passed"])
            # 32b0 on main's rule (decision 667 (1)): raw 0.7764, the chain's 2 ends on the box -> 0 straddlers,
            # 21 removed runs on 1 listed chain (information, carried): 0.7764, FAIL.
            main_32b0 = sq.exclusion_coverage_run("32b0-main-p4", config, write("32b0-main-p4", {"1": 0.7764, "2": 0.7764},
                                                                                straddling=0.0, runs=21, chains=1))
            self.assertFalse(main_32b0["Passed"])
            self.assertAlmostEqual(main_32b0["MinEndAllowanceCoverage"], 0.7764, places=6)
            self.assertEqual(main_32b0["Interfaces"]["1"]["RemovedRuns"], 21)
            self.assertEqual(main_32b0["Interfaces"]["1"]["ListedChains"], 1)
            self.assertIn("end-allowance coverage min 0.7764", main_32b0["Reason"])
            self.assertIsNone(short["Interfaces"]["1"]["RemovedRuns"])  # a 665-era record: no runs / chains
            # 665's run-neighbour reading of the same defect (straddlers ~2.1 um -> ~0.84) also fails.
            main_32b0_665 = sq.exclusion_coverage_run("32b0-main-665-p4", config, write("32b0-main-665-p4", {"1": 0.7764, "2": 0.7764},
                                                                                        straddling=2.1))
            self.assertFalse(main_32b0_665["Passed"])
            self.assertAlmostEqual(main_32b0_665["MinEndAllowanceCoverage"], (32.852 * 0.7764 + 2.1) / 32.852, places=6)
            self.assertIn("end-allowance coverage min 0.84", main_32b0_665["Reason"])
            # The 7 patterns of decision 667 (1) (REPORT 12.1; the C++ / census-v3 tests read the same values from the
            # masks): the perimeter of 20 segments of 0.05, the piece listed over all of them (listed 1.0).
            patterns = {"correct": (1.0, 0.0, 1, True), "alternating": (0.5, 0.05, 10, False), "holes1": (0.8, 0.0, 5, False),
                        "holes2": (0.7, 0.0, 4, False), "holes3": (0.7, 0.0, 3, False), "hole4": (0.8, 0.0, 2, False),
                        "like32b0": (0.3, 0.0, 3, False)}
            for name, (raw, straddling, runs, expected) in patterns.items():
                reading = sq.exclusion_coverage_run(name, config, write(name, {"1": raw, "2": raw}, listed=1.0,
                                                                       straddling=straddling, runs=runs, chains=1))
                self.assertEqual(reading["Passed"], expected, name)
                self.assertAlmostEqual(reading["MinEndAllowanceCoverage"], raw + straddling, msg=name)
                self.assertEqual(reading["Interfaces"]["1"]["RemovedRuns"], runs)
            self.assertAlmostEqual(sq.exclusion_coverage_run("alternating", config, tmp / "alternating")["MinEndAllowanceCoverage"], 0.55)
            # The upper bound stays on the raw coverage: an over-exclusion is not excused by the allowance.
            over_allowed = sq.exclusion_coverage_run("over-p4", config, write("over-p4", {"1": 1.02, "2": 1.0}, straddling=0.1))
            self.assertFalse(over_allowed["Passed"])
            # A record without StraddlingPerimeterLength (a binary before decision 665): the raw coverage alone.
            legacy_short = sq.exclusion_coverage_run("legacy-short-p4", config, write("legacy-short-p4", {"1": 0.918, "2": 0.918}, listed=6.0))
            self.assertFalse(legacy_short["Passed"])
            self.assertAlmostEqual(legacy_short["MinEndAllowanceCoverage"], 0.918)
            allowance_gate = sq.exclusion_coverage_gate([short, census_min, arc_foreign, fixed])
            self.assertTrue(allowance_gate["Passed"])
            self.assertAlmostEqual(allowance_gate["MinCoverage"], 0.918)
            self.assertAlmostEqual(allowance_gate["MinEndAllowanceCoverage"], 0.9995)
            self.assertAlmostEqual(sq.end_allowance_coverage({"Coverage": 0.5, "ExcludeSegmentsLength": 2.0,
                                                              "ExcludedPerimeterLength": 1.0, "StraddlingPerimeterLength": 0.5}),
                                   0.75)
        # The status: a failing exclusion coverage is Failed whatever the other criteria read.
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=True, identity_passed=True, gate_passed=True,
                                                 exclusion_passed=False), sq.STATUS_FAILED)
        self.assertEqual(sq.qualification_status(sq.STATUS_PENDING, dense_passed=True, identity_passed=True, gate_passed=True,
                                                 exclusion_passed=True), sq.STATUS_QUALIFIED)

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
        self.addCleanup(shutil.rmtree, self.tmp, True)
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


class CommandEvaluateTest(unittest.TestCase):
    """`evaluate` end to end on a synthetic run tree whose energies obey the p4 matrices exactly
    (the Palace CSV formats of PalaceOutputReadersTest): the library's Qualification.Verdict
    seeds the status (decisions 485 / 487 (a)) and the identity partition is recorded
    (decisions 482 / 487 (d))."""

    Q_FAB = {"SA": {(1, 1): 2.0, (1, 2): 0.5, (2, 2): 1.0}, "Domain": {(1, 1): 10.0, (1, 2): 0.0, (2, 2): 5.0}}
    Q_THIN = {"SA": {(1, 1): 1.8, (1, 2): 0.5, (2, 2): 0.9}, "Domain": {(1, 1): 9.5, (1, 2): 0.0, (2, 2): 4.9}}
    TRACE = [1.0, 2.0]
    MODEL = "spatialedgecluster_edgecount-8_e7f44561bbf0"

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="spatial-qualification-evaluate-"))
        self.addCleanup(shutil.rmtree, self.tmp, True)
        fab_p4 = {cls: sq.quadratic_form(self.Q_FAB[cls], self.TRACE) for cls in self.Q_FAB}
        thin_p4 = {cls: sq.quadratic_form(self.Q_THIN[cls], self.TRACE) for cls in self.Q_THIN}
        energies = {"fabricated-p4": fab_p4, "thin-p4": thin_p4,
                    "fabricated-p5": {cls: v * 1.005 for cls, v in fab_p4.items()},
                    "thin-p5": {"SA": thin_p4["SA"] * 1.06, "Domain": thin_p4["Domain"] * 1.004}}
        configs = {}
        for name, per_class in energies.items():
            postpro = self.tmp / name / "postpro"
            postpro.mkdir(parents=True)
            write_csv(postpro / "domain-E.csv", ["i", "E_elec (J)", "E_mag (J)"], [[1, per_class["Domain"], 0.0]])
            write_csv(postpro / "surface-Q.csv", ["i", "p_surf[1]", "Q_surf[1]"], [[1, per_class["SA"] / per_class["Domain"], 1.0]])
            write_csv(postpro / "surface-Q-edge.csv", ["i", "exc", "interface", "R (m)", "E_out (J)", "p_out", "E_ann (J)", "p_ann"],
                      [[1, 0, 1, 2e-6, 0.0, 0.0, 0.0, 0.0]])
            config = {"Problem": {"Output": str(postpro)}, "Model": {"Mesh": str(self.tmp / "identity.msh")},
                      "Solver": {"Order": int(name[-1])},
                      "Boundaries": {"Postprocessing": {"Dielectric": [{"Index": 1, "Type": "SA"}]},
                                     "PrescribedPotential": [{"Index": 1, "DataFile": str(self.tmp / "trace.csv")}]}}
            (self.tmp / name / "config.json").write_text(json.dumps(config))
            configs[name] = str(self.tmp / name / "config.json")
        for name, matrices in (("fab-matrices", self.Q_FAB), ("thin-matrices", self.Q_THIN)):
            directory = self.tmp / name
            directory.mkdir()
            write_csv(directory / "domain-response-matrix.csv", ["basis_i", "basis_j", "Q_ij (J)"],
                      [[i, j, q] for (i, j), q in matrices["Domain"].items()])
            write_csv(directory / "surface-response-matrix.csv",
                      ["interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)",
                       "Q_total_ij (J)", "Q_total_ij normal (J)", "Q_total_ij tangential (J)"],
                      [[1, 0, 2e-6, i, j, q, 0.0, 0.0, q, 0.0, 0.0] for (i, j), q in matrices["SA"].items()])
        (self.tmp / "dense-traces.json").write_text(json.dumps(
            {"Version": 1, "Model": self.MODEL, "Configs": configs, "Orders": [4, 5], "Unsupported": [],
             "Traces": [{"Name": "state-2", "Family": "T2", "Coefficients": self.TRACE, "EnergyScale": 1.0}]}))
        gate_dir = self.tmp / "c0-p4"
        gate_dir.mkdir()
        (gate_dir / "palace.json").write_text(json.dumps({"SurfaceResponse": {"Diagnostics": {"ConductorConsistency": {
            "Count": 0, "Tolerance": 0.02, "Records": [{"Model": self.MODEL, "MaxRatio": 7.7e-8, "Excluded": False}]}}}}))
        (gate_dir / "config_resolved.json").write_text(json.dumps({"Solver": {"Order": 4}}))
        self.library = self.tmp / "process-library.json"
        self.record = self.tmp / "spatial-qualification.json"

    def write_library(self, model):
        self.library.write_text(json.dumps({"Version": 3, "Models": [{"Name": self.MODEL, "LibraryQualified": False, **model}]}))

    def run_evaluate(self, *extra):
        return sq.main(["evaluate", "--dense-traces", str(self.tmp / "dense-traces.json"), "--fabricated-matrices",
                        str(self.tmp / "fab-matrices"), "--thin-matrices", str(self.tmp / "thin-matrices"), "--gate",
                        str(self.tmp / "c0-p4" / "palace.json"), "--library", str(self.library), "--output", str(self.record), *extra])

    def test_control_verdict_seeds_the_status(self):
        # The reproduced b-batch1 interim case: Verdict Failed, no QualificationStatus -> the old
        # tool stamped Qualified; now Failed, rc 1, the rule recorded.
        self.write_library({"Qualification": {"Verdict": "Failed", "Record": "qualification.json", "ReferenceAnchor": None,
                                              "Order": 4, "UnjudgedTypes": []}})
        self.assertEqual(self.run_evaluate(), 1)
        record = json.loads(self.record.read_text())
        self.assertTrue(record["DensePassed"] and record["IdentityPassed"] and record["GatePassed"])
        self.assertEqual((record["PreviousStatus"], record["Status"], record["ControlVerdict"]),
                         (sq.STATUS_PENDING, sq.STATUS_FAILED, "Failed"))
        self.assertEqual(record["ControlVerdictRule"], sq.CONTROL_VERDICT_RULE)
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertEqual((model["QualificationStatus"], model["LibraryQualified"]), (sq.STATUS_FAILED, False))
        self.assertEqual(model["SpatialQualification"]["ControlVerdict"], "Failed")
        # The production verdict (--reference none, controls passed) lifts to Qualified as before.
        self.write_library({"Qualification": {"Verdict": "PendingQualification", "Record": "qualification.json",
                                              "ReferenceAnchor": None, "Order": 4, "UnjudgedTypes": []},
                            "QualificationStatus": sq.STATUS_PENDING})
        self.assertEqual(self.run_evaluate(), 0)
        record = json.loads(self.record.read_text())
        self.assertEqual((record["Status"], record["ControlVerdict"], record["ControlVerdictRule"]),
                         (sq.STATUS_QUALIFIED, "PendingQualification", None))
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertEqual((model["QualificationStatus"], model["LibraryQualified"]), (sq.STATUS_QUALIFIED, True))
        # A model without a driver verdict (no Qualification block) is left to the criteria.
        self.write_library({})
        self.assertEqual(self.run_evaluate("--dry-run"), 0)
        self.assertIsNone(json.loads(self.record.read_text())["ControlVerdict"])

    def test_exclusion_coverage_gate_reads_every_run_of_both_twins(self):
        """Decision 657 (b): with EdgeExcludeSegments in the configs the evaluate command reads the
        EdgeExclusion record of the four dense runs AND the two basis-matrix runs; the 32b0 reading
        (0.737 on both twins) is Failed by name, the fixed reading Qualified, a run without the record
        (a binary before the decision) Failed by name."""
        segments = [[0.0, 0.0, 0.0, 1.0, 0.0, 0.0]] * 23
        for name in ("fabricated-p4", "thin-p4", "fabricated-p5", "thin-p5"):
            path = self.tmp / name / "config.json"
            config = json.loads(path.read_text())
            config["Boundaries"]["Postprocessing"]["Dielectric"][0]["EdgeExcludeSegments"] = segments
            path.write_text(json.dumps(config))

        def write_records(coverage, skip=(), listed=32.852, straddling=None):
            for name in ("fabricated-p4", "thin-p4", "fabricated-p5", "thin-p5", "fab-matrices", "thin-matrices"):
                postpro = self.tmp / name / ("postpro" if name.endswith("p4") or name.endswith("p5") else "")
                target = postpro / "palace.json"
                if name in skip:
                    target.write_text(json.dumps({"Problem": {"MeshElements": 1}}))
                    continue
                entry = {"Coverage": coverage, "ExcludeSegmentsLength": listed, "ExcludedPerimeterLength": listed * coverage,
                         "ExcludeSegmentsCount": 23, "ExcludedPerimeterCount": 400, "Tolerance": 1.9e-3, "MaxSagitta": 4.24e-3}
                if straddling is not None:
                    entry["StraddlingPerimeterLength"] = straddling
                    entry["EndAllowanceCoverage"] = (entry["ExcludedPerimeterLength"] + straddling) / listed
                target.write_text(json.dumps({"SurfaceResponse": {"Diagnostics": {"EdgeExclusion": {"Rule": "r", "Interfaces": {
                    "1": entry}}}}}))
        qualification = {"Verdict": "PendingQualification", "Record": "q.json", "ReferenceAnchor": None, "Order": 4, "UnjudgedTypes": []}
        # The 32b0 reading on main's rule: Failed, the gate named with its runs.
        write_records(0.737)
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate(), 1)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        self.assertFalse(record["ExclusionCoveragePassed"])
        self.assertTrue(record["DensePassed"] and record["IdentityPassed"] and record["GatePassed"])
        self.assertEqual(sorted(record["ExclusionCoverage"]["Failures"]),
                         ["fabricated-matrices", "fabricated-p4", "fabricated-p5", "thin-matrices", "thin-p4", "thin-p5"])
        self.assertAlmostEqual(record["ExclusionCoverage"]["MinCoverage"], 0.737)
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertEqual((model["QualificationStatus"], model["LibraryQualified"]), (sq.STATUS_FAILED, False))
        self.assertEqual(model["SpatialQualification"]["ExclusionCoverage"]["Passed"], False)
        # The fixed reading: Qualified, six runs recorded.
        write_records(0.9998)
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate(), 0)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_QUALIFIED)
        self.assertEqual(len(record["ExclusionCoverage"]["Runs"]), 6)
        self.assertTrue(all(run["Passed"] and run["Required"] == [1] for run in record["ExclusionCoverage"]["Runs"]))
        # One run without the record (the matrices reduced by an older binary): Failed by name.
        write_records(0.9998, skip=("fab-matrices",))
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate(), 1)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        self.assertEqual(list(record["ExclusionCoverage"]["Failures"]), ["fabricated-matrices"])
        self.assertIn("record missing", record["ExclusionCoverage"]["Failures"]["fabricated-matrices"])
        # Decision 665: a short, many-ended straight-foreign model with a correct exclusion (raw 0.918 of
        # 6.0 um, straddlers 0.7 um -> 1.035) is Qualified on the end allowance; 32b0 on main's rule
        # (decision 667 (1): 0 straddlers, its chain's ends on the box -> 0.7764) stays Failed.
        write_records(0.918, listed=6.0, straddling=0.7)
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate(), 0)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_QUALIFIED)
        self.assertAlmostEqual(record["ExclusionCoverage"]["MinCoverage"], 0.918)
        self.assertAlmostEqual(record["ExclusionCoverage"]["MinEndAllowanceCoverage"], (6.0 * 0.918 + 0.7) / 6.0)
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertAlmostEqual(model["SpatialQualification"]["ExclusionCoverage"]["MinEndAllowanceCoverage"], (6.0 * 0.918 + 0.7) / 6.0)
        write_records(0.7764, straddling=0.0)
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate(), 1)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        self.assertAlmostEqual(record["ExclusionCoverage"]["MinEndAllowanceCoverage"], 0.7764)

    def test_identity_partition_is_recorded_and_mandatory_for_a_multi_node_reducer(self):
        """Decisions 482 / 487 (d): the loop end's 4-node reducer against a 1-node identity twin read
        MS 1.6-5.3e-6 (> 1e-6); the identity twin now takes the reducer's partition and the (F) record
        says so.  A multi-node model needs --identity-plan with the matching partition (absent or
        differing: fail closed); a one-node model records "not supplied" without it."""
        qualification = {"Verdict": "PendingQualification", "Record": "q.json", "ReferenceAnchor": None, "Order": 4, "UnjudgedTypes": []}
        identity_config = str(self.tmp / "fabricated-p4" / "config.json")

        def plan(nodes, ranks, reducer_nodes=None, run="fabricated-p4", config=identity_config):
            path = self.tmp / f"plan-{nodes}-{ranks}-{run}.json"
            path.write_text(json.dumps({"Case": self.MODEL, "Stages": [{"Name": f"{run}-dense", "Config": config}],
                                        "IdentityPartition": {"Run": run, "Nodes": nodes, "Ranks": ranks, "MatchesReducer": True,
                                                              "Reducer": {"Nodes": reducer_nodes or nodes, "Ranks": ranks, "Job": "reducer",
                                                                          "PBSJobID": "57910.h"}}}))
            return str(path)
        # One-node reducer, no plan: recorded as not supplied; the lift proceeds.
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate(), 0)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_QUALIFIED)
        self.assertEqual(record["IdentityPartition"]["Supplied"], False)
        self.assertEqual(record["IdentityPartition"]["Note"], "not supplied, one-node reducer")
        self.assertEqual(record["IdentityPartition"]["Rule"], sq.IDENTITY_PARTITION_RULE)
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertEqual(model["SpatialQualification"]["IdentityPartition"]["Supplied"], False)
        # One-node reducer with a one-node plan: recorded; a 4-node plan against it fails closed.
        self.write_library({"Qualification": qualification})
        self.assertEqual(self.run_evaluate("--identity-plan", plan(1, 192)), 0)
        record = json.loads(self.record.read_text())
        self.assertEqual((record["IdentityPartition"]["Supplied"], record["IdentityPartition"]["Nodes"], record["IdentityPartition"]["Ranks"]),
                         (True, 1, 192))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "differs from the reducer"):
            self.run_evaluate("--identity-plan", plan(4, 768))
        # The loop-end shape: a 4-node reducer (MultiNodeReduction) - the plan is mandatory ...
        multi = {"Qualification": qualification, "MultiNodeReduction": {"Nodes": 4, "Ranks": 768, "Note": "n"}}
        self.write_library(multi)
        self.record.unlink()
        with self.assertRaisesRegex(sq.SpatialQualificationError, "mandatory"):
            self.run_evaluate()
        self.assertFalse(self.record.exists())   # nothing evaluated, nothing stamped
        # ... a 1-node twin (the decision-482 reading) fails closed, as does a partition-matched twin
        # whose plan recorded another reducer, another run or another coupon's config ...
        with self.assertRaisesRegex(sq.SpatialQualificationError, "differs from the reducer"):
            self.run_evaluate("--identity-plan", plan(1, 192))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "differs from the reducer"):
            self.run_evaluate("--identity-plan", plan(4, 768, reducer_nodes=2))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "names the run"):
            self.run_evaluate("--identity-plan", plan(4, 768, run="fabricated-p5"))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "another coupon"):
            self.run_evaluate("--identity-plan", plan(4, 768, config="/elsewhere/fabricated-p4/config.json"))
        no_partition = self.tmp / "plan-none.json"
        no_partition.write_text(json.dumps({"Case": self.MODEL, "Stages": [{"Name": "fabricated-p4-dense", "Config": identity_config}]}))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "no IdentityPartition"):
            self.run_evaluate("--identity-plan", str(no_partition))
        # ... and the reducer's partition is recorded in the (F) record and the stamped library.
        self.assertEqual(self.run_evaluate("--identity-plan", plan(4, 768)), 0)
        record = json.loads(self.record.read_text())
        self.assertEqual(record["Status"], sq.STATUS_QUALIFIED)
        partition = record["IdentityPartition"]
        self.assertEqual({key: partition[key] for key in ("Supplied", "Run", "Nodes", "Ranks", "ReducerNodes", "ReducerRanks", "ReducerJob", "MatchesReducer")},
                         {"Supplied": True, "Run": "fabricated-p4", "Nodes": 4, "Ranks": 768, "ReducerNodes": 4, "ReducerRanks": 768,
                          "ReducerJob": "57910.h", "MatchesReducer": True})
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertEqual(model["SpatialQualification"]["IdentityPartition"]["Nodes"], 4)
        # A model recording Nodes only (a library written before the Ranks field) is checked on nodes.
        self.write_library({"Qualification": qualification, "MultiNodeReduction": {"Nodes": 4, "Note": "n"}})
        self.assertEqual(self.run_evaluate("--identity-plan", plan(4, 768)), 0)
        self.assertEqual(json.loads(self.record.read_text())["IdentityPartition"]["ReducerRanks"], 768)

    def test_identity_record_reads_the_executed_twin(self):
        """Decision 495 (4): when the identity twin job's status.json is present (the sibling of the plan, or
        --identity-status) its EXECUTED node count (the runner's Nodes host list) must be the planned one and
        the run's stage complete; the executed partition is recorded; a 1-node execution of a 4-node plan (the
        decision-482 reading) or an incomplete stage fails closed."""
        qualification = {"Verdict": "PendingQualification", "Record": "q.json", "ReferenceAnchor": None, "Order": 4, "UnjudgedTypes": []}
        self.write_library({"Qualification": qualification, "MultiNodeReduction": {"Nodes": 4, "Ranks": 768, "Note": "n"}})
        job_dir = self.tmp / "le-dense" / "fabricated-p4"
        job_dir.mkdir(parents=True)
        plan = job_dir / "plan.json"
        plan.write_text(json.dumps({"Case": self.MODEL, "Stages": [{"Name": "fabricated-p4-dense", "Config": str(self.tmp / "fabricated-p4" / "config.json")}],
                                    "IdentityPartition": {"Run": "fabricated-p4", "Nodes": 4, "Ranks": 768, "MatchesReducer": True,
                                                          "Reducer": {"Nodes": 4, "Ranks": 768, "Job": "reducer", "PBSJobID": "57910.h"}}}))
        hosts = ["ip-10-0-0-1", "ip-10-0-0-2", "ip-10-0-0-3", "ip-10-0-0-4"]
        status = {"Version": 2, "PBSJobID": "58043.h", "Host": hosts[0], "Nodes": hosts, "State": "complete",
                  "Stages": [{"Name": "fabricated-p4-dense", "State": "complete", "WallSeconds": 1234.5}]}
        # No status next to the plan: the plan alone (Executed None).
        self.assertEqual(self.run_evaluate("--identity-plan", str(plan)), 0)
        self.assertIsNone(json.loads(self.record.read_text())["IdentityPartition"]["Executed"])
        # The sibling status.json of a 4-node execution: recorded.
        (job_dir / "status.json").write_text(json.dumps(status))
        self.write_library({"Qualification": qualification, "MultiNodeReduction": {"Nodes": 4, "Ranks": 768, "Note": "n"}})
        self.assertEqual(self.run_evaluate("--identity-plan", str(plan)), 0)
        executed = json.loads(self.record.read_text())["IdentityPartition"]["Executed"]
        self.assertEqual((executed["Nodes"], executed["Hosts"], executed["PBSJobID"], executed["StageState"], executed["WallSeconds"]),
                         (4, hosts, "58043.h", "complete", 1234.5))
        model = json.loads(self.library.read_text())["Models"][0]
        self.assertEqual(model["SpatialQualification"]["IdentityPartition"]["Executed"]["Nodes"], 4)
        # A 1-node execution (the runner writes no Nodes list, Host only) of the 4-node plan: fail closed.
        one_node = self.tmp / "status-1node.json"
        one_node.write_text(json.dumps({**status, "Nodes": None, "PBSJobID": "58012.h"}))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "EXECUTED on 1 node"):
            self.run_evaluate("--identity-plan", str(plan), "--identity-status", str(one_node))
        # An incomplete / absent stage, and a named status that does not exist: fail closed.
        incomplete = self.tmp / "status-incomplete.json"
        incomplete.write_text(json.dumps({**status, "State": "incomplete", "Stages": [{"Name": "fabricated-p4-dense", "State": "timed-out"}]}))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "no complete stage"):
            self.run_evaluate("--identity-plan", str(plan), "--identity-status", str(incomplete))
        with self.assertRaisesRegex(sq.SpatialQualificationError, "does not exist"):
            self.run_evaluate("--identity-plan", str(plan), "--identity-status", str(self.tmp / "missing.json"))


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
        # A coupon without EdgeExcludeSegments: the EXCLUSION-COVERAGE gate is not applicable and passes.
        self.coverage = sq.exclusion_coverage_gate([])

    def test_consistent_trace_qualifies(self):
        trace = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                  thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertTrue(trace["Passed"])
        self.assertAlmostEqual(trace["Energies"]["ModelCorrection"]["SA"], self.fab_p4["SA"] - self.thin_p4["SA"])
        self.assertEqual(trace["MatrixIdentity"]["SA"]["Residual"], 0.0)
        self.assertEqual(sorted(trace["TwinConsistency"]), ["Domain"])
        self.assertAlmostEqual(trace["ThinSurfacePSteps"]["SA"]["RelativeToThin"], 0.06 / 1.06)
        record = sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate)
        self.assertEqual(record["Status"], sq.STATUS_QUALIFIED)
        self.assertTrue(record["TwinConsistencyPassed"] and record["ClosurePassed"])
        self.assertEqual(record["Tolerances"]["GateMaxRatio"], 1e-4)
        self.assertEqual(record["Families"], ["T2"])
        window = [sq.reference_box_closure(1.02, 1.0, 0.0, True) | {"Class": "SA", "Window": "S1p"}]
        self.assertEqual(sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, reference_boxes=window)["Status"], sq.STATUS_WINDOW_VALIDATED)
        off_marker = [sq.reference_box_closure(1.2, 1.0, 0.0, True)]
        self.assertEqual(sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, reference_boxes=off_marker)["Status"], sq.STATUS_QUALIFIED)
        self.assertEqual(sq.evaluate([trace], exclusion_coverage=self.coverage, gate={"Passed": False, "Probed": 0})["Status"], sq.STATUS_FAILED)
        # The model's Qualification.UnjudgedTypes (decisions 474 / 477 (1)) refuses the lift with the reason recorded.
        unjudged = sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, reference_boxes=window, unjudged_types=["p_SA"])
        self.assertEqual(unjudged["Status"], sq.STATUS_PENDING)
        self.assertEqual(unjudged["UnjudgedTypes"], ["p_SA"])
        self.assertEqual(unjudged["UnjudgedTypesRule"], sq.UNJUDGED_TYPES_RULE)
        self.assertTrue(unjudged["DensePassed"] and unjudged["GatePassed"])
        self.assertIsNone(record["UnjudgedTypesRule"])
        self.assertEqual(record["UnjudgedTypes"], [])
        # The EXCLUSION-COVERAGE gate (decision 657 (b)) is recorded and drives the status; it is required.
        self.assertTrue(record["ExclusionCoveragePassed"])
        self.assertFalse(record["ExclusionCoverage"]["Applicable"])
        self.assertEqual(record["Tolerances"]["ExclusionCoverageMin"], 0.99)
        failing = sq.exclusion_coverage_gate([{"Run": "fabricated-p4", "Required": [1], "Passed": False, "MinCoverage": 0.737,
                                               "MaxCoverage": 0.737, "Reason": "EXCLUSION-COVERAGE FAIL", "Interfaces": {}}])
        failed = sq.evaluate([trace], exclusion_coverage=failing, gate=self.gate)
        self.assertEqual(failed["Status"], sq.STATUS_FAILED)
        self.assertFalse(failed["ExclusionCoveragePassed"])
        self.assertTrue(failed["DensePassed"] and failed["GatePassed"])
        self.assertEqual(failed["ExclusionCoverage"]["Failures"], {"fabricated-p4": "EXCLUSION-COVERAGE FAIL"})
        with self.assertRaisesRegex(sq.SpatialQualificationError, "EXCLUSION-COVERAGE"):
            sq.evaluate([trace], exclusion_coverage=None, gate=self.gate)

    def test_failed_control_verdict_is_never_lifted(self):
        """Decisions 485 / 487 (a): the b-batch1 interim libraries of S3p / S2p read Qualification.Verdict
        Failed (control 31 p_SA) next to QualificationStatus Qualified, because the (F) path seeded the
        transition from QualificationStatus alone.  The verdict now refuses the lift: Failed whatever
        the dense traces, the identity and the gate read; PendingQualification / Passed proceed."""
        trace = sq.evaluate_trace("device-patch-9", "T1", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                  thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        window = [sq.reference_box_closure(1.02, 1.0, 0.0, True) | {"Class": "SA", "Window": "S3p"}]
        failed = sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, reference_boxes=window, control_verdict="Failed")
        self.assertEqual(failed["Status"], sq.STATUS_FAILED)
        self.assertTrue(failed["DensePassed"] and failed["IdentityPassed"] and failed["GatePassed"])
        self.assertEqual(failed["ControlVerdict"], "Failed")
        self.assertEqual(failed["ControlVerdictRule"], sq.CONTROL_VERDICT_RULE)
        for verdict, status in (("PendingQualification", sq.STATUS_WINDOW_VALIDATED), ("Passed", sq.STATUS_WINDOW_VALIDATED),
                                (None, sq.STATUS_WINDOW_VALIDATED)):
            record = sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, reference_boxes=window, control_verdict=verdict)
            self.assertEqual(record["Status"], status, verdict)
            self.assertEqual(record["ControlVerdict"], verdict)
            self.assertIsNone(record["ControlVerdictRule"])
        with self.assertRaisesRegex(sq.SpatialQualificationError, "unknown control verdict"):
            sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, control_verdict="Qualified")
        # A current status of Qualified (a stamped library) is no shield either.
        self.assertEqual(sq.qualification_status(sq.STATUS_QUALIFIED, dense_passed=True, identity_passed=True, gate_passed=True,
                                                 control_verdict="Failed"), sq.STATUS_FAILED)
        # The reproduced interim case: the S3p model with the driver's Failed verdict and no
        # QualificationStatus (PendingQualification by default) stamps Failed, LibraryQualified False.
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            model = {"Name": "spatialedgecluster_edgecount-8_e7f44561bbf0", "LibraryQualified": False,
                     "Qualification": {"Verdict": "Failed", "Record": "qualification.json", "ReferenceAnchor": None, "Order": 4,
                                       "UnjudgedTypes": []}}
            path.write_text(json.dumps({"Models": [model]}))
            stamped = sq.stamp_library_status(path, model["Name"], failed)["Models"][0]
            self.assertEqual(stamped["QualificationStatus"], sq.STATUS_FAILED)
            self.assertFalse(stamped["LibraryQualified"])
            self.assertEqual(stamped["SpatialQualification"]["ControlVerdict"], "Failed")
            self.assertEqual(stamped["SpatialQualification"]["Status"], sq.STATUS_FAILED)

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
        record = sq.evaluate([good, trace], exclusion_coverage=self.coverage, gate=self.gate)
        self.assertFalse(record["DensePassed"])
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        identity_defect = sq.evaluate_trace("state-2", "T2", self.t, fab_p4={cls: v * (1 + 1e-5) for cls, v in self.fab_p4.items()},
                                            thin_p4=self.thin_p4, fab_p5=self.fab_p5, thin_p5=self.thin_p5,
                                            q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertFalse(identity_defect["MatrixIdentity"]["SA"]["Passed"])
        self.assertEqual(sq.evaluate([identity_defect], exclusion_coverage=self.coverage, gate=self.gate)["Status"], sq.STATUS_FAILED)
        # The thin twin's domain step of 3 % (a p-unstable domain defect) fails the twin-consistency.
        unstable_thin_p5 = {"SA": self.thin_p5["SA"], "Domain": self.thin_p4["Domain"] * 0.97}
        unstable = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                     thin_p5=unstable_thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
        self.assertFalse(unstable["TwinConsistency"]["Domain"]["Passed"])
        self.assertTrue(all(v["Passed"] for v in unstable["Closure"].values()))
        record = sq.evaluate([unstable], exclusion_coverage=self.coverage, gate=self.gate)
        self.assertFalse(record["TwinConsistencyPassed"])
        self.assertEqual(record["Status"], sq.STATUS_FAILED)
        with self.assertRaises(sq.SpatialQualificationError):
            sq.evaluate([], exclusion_coverage=self.coverage, gate=self.gate)

    def test_stamp_library_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "process-library.json"
            path.write_text(json.dumps({"Models": [{"Name": "m", "LibraryQualified": False}]}))
            trace = sq.evaluate_trace("state-2", "T2", self.t, fab_p4=self.fab_p4, thin_p4=self.thin_p4, fab_p5=self.fab_p5,
                                      thin_p5=self.thin_p5, q_fab=self.q_fab, q_thin=self.q_thin)
            record = sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate)
            library = sq.stamp_library_status(path, "m", record)
            model = library["Models"][0]
            self.assertEqual(model["QualificationStatus"], sq.STATUS_QUALIFIED)
            self.assertTrue(model["LibraryQualified"])
            self.assertEqual(model["SpatialQualification"]["Families"], ["T2"])
            self.assertTrue(model["SpatialQualification"]["TwinConsistencyPassed"])
            self.assertEqual(model["SpatialQualification"]["Gate"]["Orders"], [4, 5])
            self.assertEqual(model["SpatialQualification"]["Traces"][0]["Name"], "state-2")
            self.assertEqual(model["SpatialQualification"]["ExclusionCoverage"]["Passed"], True)
            self.assertEqual(model["SpatialQualification"]["ExclusionCoverage"]["Applicable"], False)
            with self.assertRaises(sq.SpatialQualificationError):
                sq.stamp_library_status(path, "other", record)
            # An unjudged Type stamps PendingQualification, LibraryQualified stays False (decision 477 (1)).
            unjudged = sq.evaluate([trace], exclusion_coverage=self.coverage, gate=self.gate, unjudged_types=["p_MS"])
            model = sq.stamp_library_status(path, "m", unjudged)["Models"][0]
            self.assertEqual(model["QualificationStatus"], sq.STATUS_PENDING)
            self.assertFalse(model["LibraryQualified"])
            self.assertEqual(model["SpatialQualification"]["UnjudgedTypes"], ["p_MS"])


if __name__ == "__main__":
    unittest.main()
