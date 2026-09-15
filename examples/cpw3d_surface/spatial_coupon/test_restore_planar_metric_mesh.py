# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest

import numpy as np

from restore_planar_metric_mesh import (
    _bounded_offset,
    _movement_basis,
    _quality_repair,
    _tetra_quality,
)


class PlanarMetricQualityRepairTest(unittest.TestCase):
    def test_parameterization_bounds_free_surface_and_ridge_vertices(self):
        maximum = .01875
        supports = {
            100: {"Normal": [1., 0., 0.]},
            101: {"Normal": [0., 1., 0.]},
        }
        node_supports = {1: {100}, 2: {100, 101}}
        fixed_nodes = set()
        parameters = np.array([.75, .75, .75])

        free = _movement_basis(0, node_supports, supports, fixed_nodes)
        free_offset = _bounded_offset(free, parameters, maximum)
        surface = _movement_basis(1, node_supports, supports, fixed_nodes)
        surface_offset = _bounded_offset(surface, parameters[:2], maximum)
        ridge = _movement_basis(2, node_supports, supports, fixed_nodes)
        ridge_offset = _bounded_offset(ridge, parameters[:1], maximum)

        self.assertEqual(free.shape, (3, 3))
        self.assertLessEqual(np.linalg.norm(free_offset), maximum)
        self.assertGreater(np.linalg.norm(parameters), 1.)
        self.assertLessEqual(np.linalg.norm(surface_offset), maximum)
        self.assertAlmostEqual(float(np.dot(surface_offset, [1., 0., 0.])), 0.)
        self.assertLessEqual(np.linalg.norm(ridge_offset), maximum)
        self.assertAlmostEqual(float(np.dot(ridge_offset, [1., 0., 0.])), 0.)
        self.assertAlmostEqual(float(np.dot(ridge_offset, [0., 1., 0.])), 0.)
        self.assertEqual(_movement_basis(2, node_supports, supports, {2}).shape,
                         (3, 0))

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
        displacement = np.linalg.norm(points - original, axis=1)
        self.assertLessEqual(displacement.max(), .01875)
        self.assertEqual(report["MaximumFinalQualityDisplacementUm"],
                         displacement.max())
        self.assertLessEqual(report["MaximumQualityStepDisplacementUm"], .01875)
        self.assertEqual(report["QualityDisplacementBoundUm"], .01875)
        self.assertEqual(report["MaximumSupportConstraintError"], 0.)
        self.assertEqual(report["QualityRepairVertices"], 3)

    def test_cumulative_bound_when_vertices_are_repaired_in_multiple_passes(self):
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
        maximum = .01875

        report = _quality_repair(points, tetrahedra, node_supports, supports,
                                 recipe, .6, 4., maximum)
        displacement = np.linalg.norm(points - original, axis=1)

        # The corner pass is followed by a low-quality-component pass over the
        # same tetrahedron, so its movable vertices participate more than once.
        self.assertEqual(report["QualityRepairComponents"], 1)
        self.assertLessEqual(displacement.max(), maximum)
        self.assertEqual(report["MaximumFinalQualityDisplacementUm"],
                         displacement.max())
        self.assertLessEqual(report["MaximumQualityStepDisplacementUm"], maximum)
        self.assertTrue(np.all(points[:3, 1] == 0.))


if __name__ == "__main__":
    unittest.main()
