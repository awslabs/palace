#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""analyze_basis_needles.py on the repository testdata (supervisor decision 54b)."""
import unittest
from pathlib import Path

import numpy as np

import analyze_basis_needles as needles
import estimate_build_cost as cost

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "geometry-independence-suite.json"


class AnalyzeBasisNeedlesTest(unittest.TestCase):
    def test_ten_edge_decomposition_matches_the_estimator_and_the_delaunay_caps_remove_the_needles(self):
        result = needles.analyze(HERE / "testdata" / "ten-edge-6791f1c84123", MANIFEST, [2.0])
        current = result["Current"]
        # The component sum is the estimator's total; the trace-basis law is the sum of its classes.
        components = current["Components"]
        self.assertAlmostEqual(sum(components[k] for k in ("FarField", "TubeBand", "CornerBallsAndCaps", "JunctionProxy",
                                                            "TraceBasis", "TubePrisms", "TubePyramids")),
                               current["EstimatedElements"], places=6)
        self.assertAlmostEqual(components["TraceBasisNeedles"] + components["TraceBasisRightSlivers"],
                               components["TraceBasis"], places=6)
        self.assertEqual(current["Classes"]["Needle"], 34)
        self.assertGreater(current["NeedleShare"], 0.2)
        self.assertLess(current["MinimumAltitude"]["Needle"], 0.006)
        # Delaunay caps: same sources, the needle tets vanish to within 2% of the estimate, the
        # estimate falls; the graded model adds sources faster than it removes elements.
        delaunay = result["Delaunay"]
        self.assertEqual(delaunay["Sources"], result["Sources"])
        self.assertEqual(sorted(delaunay["RetriangulatedFaces"]), [2, 3])
        self.assertLess(delaunay["TraceBasisNeedles"], 0.02 * current["EstimatedElements"])
        self.assertLess(delaunay["ElementsOverCurrent"], 0.8)
        self.assertGreater(delaunay["MinimumAltitude"], 0.03)
        graded = result["GradedModel"][0]
        self.assertGreater(graded["SourcesOverCurrent"], 2.0)
        self.assertGreater(graded["PhysicsCostOverCurrent"], 1.0)

    def test_four_edge_has_no_needles_and_a_needle_free_basis_is_a_fixed_point(self):
        result = needles.analyze(HERE / "testdata" / "four-edge-9d2cb9bbb3fe", MANIFEST, [2.0])
        self.assertEqual(result["Current"]["Classes"]["Needle"], 0)
        self.assertEqual(result["Current"]["NeedleShare"], 0.0)
        self.assertEqual(result["Current"]["Components"]["TraceBasisNeedles"], 0.0)

    def test_classification_of_a_needle_a_right_sliver_and_a_wide_triangle(self):
        vertices = {1: np.array([0.0, 0.0, 0.0]), 2: np.array([0.05, 0.0, 0.0]), 3: np.array([10.0, 0.5, 0.0]),
                    4: np.array([0.0, 2.0, 0.0]), 5: np.array([5.0, 0.0, 0.0]), 6: np.array([0.0, 5.0, 0.0])}
        rows = needles.classify(vertices, [(1, 2, 3), (1, 2, 4), (1, 5, 6)], 0.5, 0.16)
        self.assertEqual([row["Class"] for row in rows], ["Needle", "RightSliver", "Wide"])
        # 2 x area / longest edge of the (0, 0), (0.05, 0), (0, 2) right sliver: 0.05 x 2 / hypot(0.05, 2).
        self.assertAlmostEqual(rows[1]["Altitude"], 0.1 / np.hypot(0.05, 2.0))
        integrals = needles.trace_integral_by_class(rows, 0.16, 0.5)
        self.assertGreater(integrals["Needle"], integrals["RightSliver"])
        self.assertEqual(integrals["Wide"], 0.0)
        # The class integral is the estimator's per-triangle shell integral.
        row = rows[1]
        expected = cost.shell_integral(0.0, (0.16 - row["RequestedSize"]) / 0.5, lambda d: row["RequestedSize"] + 0.5 * d,
                                       lambda d: row["Area"] + 0.5 * np.pi * row["Perimeter"] * d)
        self.assertAlmostEqual(integrals["RightSliver"], expected)


if __name__ == "__main__":
    unittest.main()
