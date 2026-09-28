#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The compact corner library file carries the within-R energy the library adds to a
device (`Q_ij (J)` summed over the physical edges) next to the whole-box `Q_total_ij (J)`
and the matching radius; the within-R matrices are what the finalizer returns."""

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FINALIZE = load("finalize_corner_response")


class AggregateLocalizedSurfaceMatrixTest(unittest.TestCase):
    def write_source(self, path, radii):
        header = [
            "interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)",
            "Q_ij normal (J)", "Q_ij tangential (J)", "Q_total_ij (J)",
            "Q_total_ij normal (J)", "Q_total_ij tangential (J)",
        ]
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(header)
            for radius in radii:
                for interface in (1, 2):
                    for edge in (1, 2):
                        for i in (1, 2):
                            for j in range(i, 3):
                                within = 1.0e-20 * (interface + 0.1 * edge + 0.01 * (i + j))
                                total = 3.0 * within
                                writer.writerow(
                                    [interface, edge, f"{radius:.6e}", i, j, within, 0.0,
                                     within, total, 0.0, total]
                                )

    def test_within_r_and_whole_box_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "surface-response-matrix.csv"
            destination = Path(directory) / "surface-response-matrix-aggregate.csv"
            self.write_source(source, [1.9e-6])
            matrices = FINALIZE.aggregate_localized_surface_matrix(source, destination, 2)
            with destination.open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(
                list(rows[0]),
                ["interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)",
                 "Q_total_ij (J)"],
            )
            self.assertEqual(len(rows), 2 * 3)
            for row in rows:
                self.assertEqual(int(row["edge"]), 1)
                self.assertAlmostEqual(float(row["R (m)"]), 1.9e-6, delta=1e-18)
                interface = int(row["interface"])
                i, j = int(row["basis_i"]), int(row["basis_j"])
                within = sum(1.0e-20 * (interface + 0.1 * edge + 0.01 * (i + j)) for edge in (1, 2))
                self.assertAlmostEqual(float(row["Q_ij (J)"]), within, delta=1e-32)
                self.assertAlmostEqual(float(row["Q_total_ij (J)"]), 3.0 * within, delta=1e-32)
                self.assertAlmostEqual(matrices[interface][i - 1, j - 1], within, delta=1e-32)
                self.assertAlmostEqual(matrices[interface][j - 1, i - 1], within, delta=1e-32)

    def test_refuses_several_radii(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "surface-response-matrix.csv"
            destination = Path(directory) / "aggregate.csv"
            self.write_source(source, [1.9e-6, 2.0e-6])
            with self.assertRaises(ValueError):
                FINALIZE.aggregate_localized_surface_matrix(source, destination, 2)


if __name__ == "__main__":
    unittest.main()
