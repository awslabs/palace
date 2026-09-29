#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""library_basis_gate: no free trace basis function with support on the PEC part of the box
contour (corner-family review 2026-09-29). A corner coupon of the lane-2 angle-independent
layout at 165 degrees fails (the second arm crosses the ring at no PEC knot), the trace basis
rule's layout passes; a two-dimensional edge coupon with a knot on the metal band fails, one
whose knots straddle the band passes; a strip passes by construction; a parallel-edge
cluster is NotEvaluable (fails closed)."""
import csv
import importlib.util
import json
import os
import sys
import tempfile
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import library_basis_gate as gate  # noqa: E402

CORNER_ROOT = os.path.join(HERE, "..", "..", "corner_coupon")


def load_generator():
    spec = importlib.util.spec_from_file_location(
        "generate_corner_response", os.path.join(CORNER_ROOT, "generate_corner_response.py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GENERATOR = load_generator()
R, T, OE = 1.9, 0.1, 0.05


def write_points(path, points):
    with open(path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["x", "y", "z"])
        writer.writerows(points)


def write_surface(path, size, diagonal):
    with open(path, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["interface", "edge", "R (m)", "basis_i", "basis_j", "Q_ij (J)", "Q_total_ij (J)"])
        for i in range(1, size + 1):
            writer.writerow([2, 1, 1.9e-6, i, i, diagonal[i - 1], diagonal[i - 1]])


class LibraryBasisGateTest(unittest.TestCase):
    def corner_model(self, root, name, angle, rule_layout):
        if rule_layout:
            surface = GENERATOR.build_surface(R, 8, T, OE, angle_degrees=angle, topology="convex")
            points = np.asarray(surface.knot_points)
            zero = surface.zero_trace_indices()
            groups = surface.contour_groups
        else:
            surface = GENERATOR.build_surface(R, 8, T, OE)
            points = np.asarray(surface.knot_points)
            zero = np.flatnonzero(GENERATOR.pec_trace_mask(points, R, angle, 0.0, T, 90.0, "convex")).tolist()
            groups = surface.contour_groups
        write_points(os.path.join(root, name + ".csv"), points)
        diagonal = np.full(len(points), 1.0e-19)
        write_surface(os.path.join(root, name + "-surface.csv"), len(points), diagonal)
        return {"Name": name, "Topology": "ConvexCorner", "Angle": angle, "AngleDegrees": angle,
                "CornerRadius": 0.0, "BasisPoints": name + ".csv",
                "FabricatedSurfaceMatrix": name + "-surface.csv",
                "ContourGroups": groups, "ZeroTraceIndices": [i + 1 for i in zero],
                "Interfaces": [{"Type": "MS", "Coupon": 2}]}

    def library(self, root, models):
        path = os.path.join(root, "process-library.json")
        json.dump({"Version": 3, "Name": "test", "MatchingRadius": R,
                   "Fabrication": {"MetalThickness": T}, "Models": models}, open(path, "w"))
        return path

    def test_corner_lane2_layout_fails_and_rule_layout_passes(self):
        with tempfile.TemporaryDirectory() as root:
            record = gate.evaluate(self.library(root, [
                self.corner_model(root, "c90", 90.0, False),
                self.corner_model(root, "c165-lane2", 165.0, False),
                self.corner_model(root, "c165-rule", 165.0, True),
                self.corner_model(root, "c120-rule", 120.0, True)]))
            by_name = {m["Name"]: m for m in record["Models"]}
            self.assertTrue(by_name["c90"]["Passed"])
            self.assertFalse(by_name["c165-lane2"]["Passed"])
            self.assertTrue(any("no PEC knot" in f for f in by_name["c165-lane2"]["Findings"]))
            self.assertTrue(by_name["c165-rule"]["Passed"], by_name["c165-rule"]["Findings"])
            self.assertTrue(by_name["c120-rule"]["Passed"], by_name["c120-rule"]["Findings"])
            self.assertEqual(record["Verdict"], "Failed")
            self.assertEqual(by_name["c165-rule"]["Diagonal"]["Reference"], "c90")

    def test_diagonal_factor(self):
        with tempfile.TemporaryDirectory() as root:
            models = [self.corner_model(root, "c90", 90.0, True), self.corner_model(root, "c165", 165.0, True)]
            surface = GENERATOR.build_surface(R, 8, T, OE, angle_degrees=165.0, topology="convex")
            diagonal = np.full(72, 1.0e-19)
            free = [i for i in range(72) if i not in surface.zero_trace_indices()]
            diagonal[free[0]] = 5.0e-18  # 50x the reference's free-knot maximum
            write_surface(os.path.join(root, "c165-surface.csv"), 72, diagonal)
            record = gate.evaluate(self.library(root, models))
            by_name = {m["Name"]: m for m in record["Models"]}
            self.assertFalse(by_name["c165"]["Passed"])
            self.assertAlmostEqual(by_name["c165"]["Diagonal"]["Factors"]["MS"], 50.0)

    def test_two_dimensional_models(self):
        with tempfile.TemporaryDirectory() as root:
            good = [(-R, -0.3, 0.0), (-R, 0.4, 0.0), (0.0, -R, 0.0), (R, 0.0, 0.0), (0.0, R, 0.0)]
            bad = good + [(-R, 0.05, 0.0)]
            write_points(os.path.join(root, "good.csv"), good)
            write_points(os.path.join(root, "bad.csv"), bad)
            write_points(os.path.join(root, "strip.csv"), good)
            write_points(os.path.join(root, "cluster.csv"), good)
            record = gate.evaluate(self.library(root, [
                {"Name": "edge", "Topology": "IsolatedEdge", "BasisPoints": "good.csv"},
                {"Name": "edge-bad", "Topology": "IsolatedEdge", "BasisPoints": "bad.csv"},
                {"Name": "strip", "Topology": "SameConductorStrip", "Separation": 2.0, "BasisPoints": "strip.csv"},
                {"Name": "cluster", "Topology": "ParallelEdgeCluster", "BasisPoints": "cluster.csv"}]))
            by_name = {m["Name"]: m for m in record["Models"]}
            self.assertTrue(by_name["edge"]["Passed"])
            self.assertFalse(by_name["edge-bad"]["Passed"])
            self.assertTrue(by_name["strip"]["Passed"])
            self.assertFalse(by_name["cluster"]["Passed"])
            self.assertTrue(by_name["cluster"]["Findings"][0].startswith("NotEvaluable"))
            self.assertEqual(record["Verdict"], "Failed")


if __name__ == "__main__":
    unittest.main()
