# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest

import numpy as np

from restore_planar_metric_mesh import _quality_repair, _tetra_quality


class PlanarMetricQualityRepairTest(unittest.TestCase):
    def test_repair_improves_corner_and_sidewall_without_leaving_support(self):
        points = np.array([[0., 0., 0.], [0., 0., .05], [.05, 0., 0.],
                           [.01, .0002, .01]])
        original = points.copy()
        tetrahedra = np.array([[0, 1, 2, 3]])
        supports = {100: {"Normal": [0., 1., 0.], "Offset": 0.,
                          "Attribute": 6001}}
        node_supports = {0: {100}, 1: {100}, 2: {100}}
        recipe = {
            "SemanticContract": {
                "CutSurfaceRoles": ["matching-surface"],
                "BoundaryLabels": [
                    {"Role": "physical-sidewall", "Attribute": 6001}
                ]
            },
            "TruePhysicalCorners": [[0., 0., 0.]]
        }
        before_scaled, before_aspect, _ = _tetra_quality(points, tetrahedra)
        report = _quality_repair(points, tetrahedra, node_supports, supports,
                                 recipe, .01, 4., .01875)
        after_scaled, after_aspect, determinant = _tetra_quality(points, tetrahedra)
        self.assertLess(before_scaled[0], .02)
        self.assertGreater(before_aspect[0], 4.)
        self.assertGreaterEqual(after_scaled[0], .01)
        self.assertLessEqual(after_aspect[0], 4.)
        self.assertGreater(determinant[0], 0.)
        self.assertTrue(np.all(points[:3, 1] == 0.))
        self.assertTrue(np.array_equal(points[0], original[0]))
        self.assertGreater(points[3, 1], original[3, 1])
        self.assertEqual(report["MaximumSupportConstraintError"], 0.)
        self.assertEqual(report["QualityRepairVertices"], 3)


if __name__ == "__main__":
    unittest.main()
