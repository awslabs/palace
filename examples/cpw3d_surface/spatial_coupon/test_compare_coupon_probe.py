#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "compare_coupon_probe", Path(__file__).with_name("compare_coupon_probe.py"))
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class ProbeComparisonTest(unittest.TestCase):
    def test_thin_does_not_read_or_gate_raw_spr(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/"a", Path(tmp)/"b"
            a.mkdir(); b.mkdir()
            for p in (a, b):
                (p/"domain-E.csv").write_text("i,E_elec (J)\n1,1e-15\n")
            (a/"surface-Q.csv").write_text("i,p_surf[1]\n1,1e-4\n")
            (b/"surface-Q.csv").write_text("i,p_surf[1]\n1,1e-2\n")
            report = module.compare(a, b, "thin")
            self.assertTrue(report["Passed"])
            self.assertFalse(report["ThinRawSPRGated"])
            self.assertEqual(len(report["Checks"]), 1)
            self.assertFalse(module.compare(a, b, "fabricated")["Passed"])

    def test_domain_error_is_not_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/"a", Path(tmp)/"b"
            a.mkdir(); b.mkdir()
            (a/"domain-E.csv").write_text("i,E_elec (J)\n1,1e-15\n")
            (b/"domain-E.csv").write_text("i,E_elec (J)\n1,1.01e-15\n")
            self.assertFalse(module.compare(a, b, "thin")["Passed"])

    def test_archived_matrix_probe_matches_ordinary_energy(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/"a", Path(tmp)/"b"
            a.mkdir(); b.mkdir()
            (a/"domain-E.csv").write_text("i,E_elec (J)\n1,1e-15\n")
            (a/"surface-Q.csv").write_text("i,p_surf[1]\n1,1e-4\n")
            (b/"domain-response-matrix.csv").write_text("basis_i,basis_j,Q_ij (J)\n1,1,1e-15\n")
            (b/"surface-response-matrix.csv").write_text(
                "basis_i,basis_j,interface,Q_ij (J),Q_total_ij (J)\n1,1,1,8e-20,1e-19\n")
            self.assertTrue(module.compare(a, b, "fabricated")["Passed"])
            self.assertEqual(len(module.compare(a, b, "thin")["Checks"]), 1)
            (b/"surface-response-matrix.csv").write_text(
                "basis_i,basis_j,interface,Q_total_ij (J)\n1,1,1,1e-19\n")
            with self.assertRaises(ValueError):
                module.compare(a, b, "fabricated")
            (b/"domain-response-matrix.csv").write_text(
                "basis_i,basis_j,Q_ij (J)\n1,1,1e-15\n1,1,1e-15\n")
            with self.assertRaises(ValueError):
                module.compare(a, b, "thin")

    def test_engineering_and_strict_are_distinct_diagnostics(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/"a", Path(tmp)/"b"
            a.mkdir(); b.mkdir()
            for p in (a, b):
                (p/"domain-E.csv").write_text("i,E_elec (J)\n1,1e-15\n")
            (a/"surface-Q.csv").write_text("i,p_surf[1]\n1,1e-4\n")
            (b/"surface-Q.csv").write_text("i,p_surf[1]\n1,1.003e-4\n")
            report = module.compare(a, b, "fabricated")
            self.assertTrue(report["Passed"])
            self.assertFalse(report["LibraryQualified"])
            self.assertFalse(module.compare(a, b, "fabricated", *module.PROFILES["strict"])["Passed"])
            with self.assertRaises(ValueError):
                module.compare(a, b, "fabricated", surface_tol=float("nan"))

    def test_empty_data_fails_closed(self):
        with tempfile.TemporaryDirectory() as tmp:
            a, b = Path(tmp)/"a", Path(tmp)/"b"
            a.mkdir(); b.mkdir()
            for p in (a, b):
                (p/"domain-E.csv").write_text("i,E_elec (J)\n")
            with self.assertRaises(ValueError):
                module.compare(a, b, "thin")


if __name__ == "__main__":
    unittest.main()
