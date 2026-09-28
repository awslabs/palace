#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Domain-defect convergence gate of qualify_surface_response.py (USER decision 117(3)):
the gate is decided by the change of the held-out correction energy c^T D c between
consecutive orders, the former Frobenius change of D is recorded next to it."""

import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

CPW2D = Path(__file__).resolve().parent
QUALIFIER = CPW2D / "qualify_surface_response.py"


def load(name):
    spec = importlib.util.spec_from_file_location(name, CPW2D / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


QUALIFY = load("qualify_surface_response")

SIZE = 4
COEFFICIENTS = np.array([0.4, 0.3, 0.2, 0.01])


def write_matrix(path, matrix):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write("  basis_i,  basis_j,                   Q_ij (J)\n")
        for i in range(SIZE):
            for j in range(i, SIZE):
                stream.write(f" {i + 1:.2e}, {j + 1:.2e}, {matrix[i, j]:+.12e}\n")


def write_surface_matrix(path, matrix):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(
            "interface,     edge,  basis_i,  basis_j,                   Q_ij (J),"
            "             Q_total_ij (J)\n"
        )
        for i in range(SIZE):
            for j in range(i, SIZE):
                stream.write(
                    f" 1.00e+00, 1.00e+00, {i + 1:.2e}, {j + 1:.2e}, "
                    f"{matrix[i, j]:+.12e}, {matrix[i, j]:+.12e}\n"
                )


def write_coupon(root, thin, fabricated):
    """A one-model coupon directory as generate_edge_response.py lays it out, with
    exact held-out solves (the direct energies equal c^T Q c) so only the convergence
    checks decide the verdict."""
    root.mkdir(parents=True, exist_ok=True)
    surface = 0.01 * np.eye(SIZE)
    for kind, matrix in (("thin", thin), ("fabricated", fabricated)):
        postpro = root / "postpro" / f"edge_{kind}"
        write_matrix(postpro / "domain-response-matrix.csv", matrix)
        write_surface_matrix(postpro / "surface-response-matrix.csv", surface)
        heldout = root / "postpro" / f"heldout_edge_{kind}"
        heldout.mkdir(parents=True, exist_ok=True)
        domain = float(COEFFICIENTS @ matrix @ COEFFICIENTS)
        interface = float(COEFFICIENTS @ surface @ COEFFICIENTS)
        (heldout / "domain-E.csv").write_text(
            "        i,                 E_elec (J)\n"
            f" 1.00e+00,        {domain:+.12e}\n"
        )
        (heldout / "surface-Q.csv").write_text(
            "        i,                  p_surf[1]\n"
            f" 1.00e+00,        {interface / domain:+.12e}\n"
        )
    (root / "heldout_coefficients.csv").write_text(
        "coefficient_V\n" + "".join(f"{float(c)!r}\n" for c in COEFFICIENTS)
    )
    (root / "process-library.json").write_text(
        json.dumps(
            {
                "Version": 1,
                "Models": [
                    {
                        "Name": "isolated-edge",
                        "Topology": "IsolatedEdge",
                        "FabricatedMatrix": "postpro/edge_fabricated/domain-response-matrix.csv",
                        "ThinMatrix": "postpro/edge_thin/domain-response-matrix.csv",
                        "FabricatedSurfaceMatrix": "postpro/edge_fabricated/surface-response-matrix.csv",
                        "ThinSurfaceMatrix": "postpro/edge_thin/surface-response-matrix.csv",
                        "Interfaces": [{"Type": "SA", "Coupon": 1}],
                    }
                ],
            }
        )
    )


def symmetric(diagonal, off=0.0):
    matrix = np.diag(np.asarray(diagonal, dtype=float))
    matrix[0, 1] = matrix[1, 0] = off
    return matrix


def qualify(current, previous):
    output = current / "qualification.json"
    completed = subprocess.run(
        [sys.executable, str(QUALIFIER), str(current), "--previous", str(previous), "--output", str(output)],
        capture_output=True,
        text=True,
    )
    if not output.is_file():
        raise AssertionError(completed.stderr)
    return completed.returncode, json.loads(output.read_text())


def defect_check(report):
    return next(c for c in report["ConvergenceChecks"] if c["Quantity"] == "domain-defect")


class DomainDefectGateTest(unittest.TestCase):
    def test_noisy_defect_matrix_with_converged_correction_energy_passes(self):
        """The decision-95 / 97 pathology: D changes by ~50 % in Frobenius norm through an
        edge hat the held-out trace barely weights while c^T D c changes by 0.3 %."""
        thin = symmetric([10.0, 10.0, 10.0, 10.0])
        previous_fab = thin + symmetric([0.05, 0.02, 0.0, 0.0])
        current_fab = thin + symmetric([0.05 * 1.003, 0.02, 0.0, 0.5])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_coupon(root / "p3", thin, previous_fab)
            write_coupon(root / "p4", thin, current_fab)
            code, report = qualify(root / "p4", root / "p3")
        check = defect_check(report)
        self.assertEqual(report["Version"], 2)
        self.assertEqual(code, 0)
        self.assertTrue(report["Passed"])
        self.assertTrue(check["Passed"])
        c = COEFFICIENTS
        expected = c @ (current_fab - thin) @ c
        previous = c @ (previous_fab - thin) @ c
        self.assertAlmostEqual(check["HeldoutCorrectionEnergy"], expected)
        self.assertAlmostEqual(check["PreviousHeldoutCorrectionEnergy"], previous)
        self.assertAlmostEqual(
            check["HeldoutCorrectionEnergyChangePercent"],
            100.0 * abs(expected - previous) / abs(expected),
        )
        self.assertLess(check["HeldoutCorrectionEnergyChangePercent"], 1.0)
        self.assertEqual(check["HeldoutCorrectionEnergyLimitPercent"], 5.0)
        former = check["PreviousDefinition"]
        self.assertGreater(former["MatrixChangePercent"], 50.0)
        self.assertFalse(former["Passed"])
        self.assertEqual(former["MatrixLimitPercent"], 5.0)
        self.assertIn("held-out correction energy", check["Definition"])
        self.assertIn("Frobenius", former["Definition"])

    def test_converged_matrix_with_moving_correction_energy_fails(self):
        """The converse: a defect matrix converged in norm whose held-out correction
        energy still moves by 8 % fails the gate at the recorded 5 %."""
        thin = symmetric([10.0, 10.0, 10.0, 10.0])
        previous_fab = thin + symmetric([1.0, 1.0, 1.0, 1.0])
        current_fab = thin + symmetric([1.08, 1.08, 1.08, 1.08])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_coupon(root / "p3", thin, previous_fab)
            write_coupon(root / "p4", thin, current_fab)
            code, report = qualify(root / "p4", root / "p3")
        check = defect_check(report)
        self.assertEqual(code, 1)
        self.assertFalse(report["Passed"])
        self.assertFalse(check["Passed"])
        self.assertAlmostEqual(check["HeldoutCorrectionEnergyChangePercent"], 100.0 * 0.08 / 1.08)
        self.assertAlmostEqual(check["PreviousDefinition"]["MatrixChangePercent"], 100.0 * 0.08 / 1.08)
        self.assertFalse(check["PreviousDefinition"]["Passed"])

    def test_limit_option_applies_to_the_energy_change(self):
        thin = symmetric([10.0, 10.0, 10.0, 10.0])
        previous_fab = thin + symmetric([1.0, 1.0, 1.0, 1.0])
        current_fab = thin + symmetric([1.08, 1.08, 1.08, 1.08])
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_coupon(root / "p3", thin, previous_fab)
            write_coupon(root / "p4", thin, current_fab)
            output = root / "p4" / "loose.json"
            completed = subprocess.run(
                [
                    sys.executable, str(QUALIFIER), str(root / "p4"), "--previous", str(root / "p3"),
                    "--output", str(output), "--max-domain-defect-change", "10",
                ],
                capture_output=True,
                text=True,
            )
            report = json.loads(output.read_text())
        self.assertEqual(completed.returncode, 0, completed.stderr)
        check = defect_check(report)
        self.assertTrue(check["Passed"])
        self.assertEqual(check["HeldoutCorrectionEnergyLimitPercent"], 10.0)
        self.assertEqual(check["PreviousDefinition"]["MatrixLimitPercent"], 10.0)

    def test_domain_defect_check_uses_the_active_trace_only(self):
        """A ZeroTraceIndices knot is removed from both the matrices and the coefficients."""
        thin = symmetric([10.0, 10.0, 10.0, 10.0])
        previous = {"responses": {"thin": {"domain": thin}, "fabricated": {"domain": thin + symmetric([1.0, 1.0, 1.0, 100.0])}}}
        current = {
            "responses": {"thin": {"domain": thin}, "fabricated": {"domain": thin + symmetric([1.0, 1.0, 1.0, 1.0])}},
            "coefficients": COEFFICIENTS,
        }
        active = np.array([0, 1, 2])
        check = QUALIFY.domain_defect_check(current, previous, active, 5.0, 5.0)
        self.assertTrue(check["Passed"])
        self.assertEqual(check["HeldoutCorrectionEnergyChangePercent"], 0.0)
        self.assertAlmostEqual(check["HeldoutCorrectionEnergy"], float(COEFFICIENTS[:3] @ np.eye(3) @ COEFFICIENTS[:3]))


if __name__ == "__main__":
    unittest.main()
