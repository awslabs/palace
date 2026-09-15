# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest
from types import SimpleNamespace
from unittest import mock

import meshio
import numpy as np

import restore_planar_metric_mesh
from restore_planar_metric_mesh import (
    DISPLACEMENT_ROUNDOFF_TOLERANCE,
    _bounded_offset,
    _movement_basis,
    _quality_repair,
    _tetra_quality,
    _transactional_quality_commit,
    _within_displacement_bound,
    restore_in_source_frame,
)


def _single_tetrahedron_repair_case(apex_height):
    points = np.array([[0., 0., 0.], [0., 0., .05], [.05, 0., 0.],
                       [.01, apex_height, .01]])
    tetrahedra = np.array([[0, 1, 2, 3]])
    supports = {100: {"Normal": [0., 1., 0.], "Offset": 0., "Attribute": 6001}}
    node_supports = {0: {100}, 1: {100}, 2: {100}}
    recipe = {
        "SemanticContract": {
            "CutSurfaceRoles": ["matching-surface"],
            "BoundaryLabels": [{"Role": "physical-sidewall", "Attribute": 6001}]},
        "TruePhysicalCorners": [[0., 0., 0.]]}
    return points, tetrahedra, supports, node_supports, recipe


class PlanarMetricQualityRepairTest(unittest.TestCase):
    def test_tilted_restoration_is_performed_in_source_local_frame(self):
        points = np.array([[0., .001, 0.], [0., -.001, .05], [.05, .0005, 0.],
                           [.01, .02, .01]])
        cells = [("triangle", np.array([[0, 1, 2]])),
                 ("tetra", np.array([[0, 1, 2, 3]]))]
        refs = [np.array([100]), np.array([1])]
        angle = .47
        rotation = np.array([[np.cos(angle), 0., np.sin(angle)], [0., 1., 0.],
                             [-np.sin(angle), 0., np.cos(angle)]])
        translation = np.array([1.2, -.4, .7])
        matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = translation
        transformed = meshio.Mesh(points @ rotation.T + translation, cells,
                                  cell_data={"medit:ref": refs})
        normal = rotation @ np.array([0., 1., 0.])
        recipe = {
            "PlanarSupports": {"100": {"Normal": normal.tolist(),
                "Offset": float(np.dot(normal, translation)), "Attribute": 6001}},
            "TruePhysicalCorners": [(np.array([0., 0., 0.]) @ rotation.T +
                                      translation).tolist()],
            "SemanticContract": {"RigidTransform": matrix.reshape(-1).tolist(),
                "VolumeMaterials": [{"Attribute": 1, "Material": "vacuum"}],
                "BoundaryLabels": [{"Attribute": 6001, "Role": "side",
                                    "AdjacentMaterials": [1]}],
                "CutSurfaceRoles": []}}
        local, published, report = restore_in_source_frame(
            transformed, recipe, .01)
        self.assertEqual(report["RestorationFrame"], "SourceLocal")
        self.assertLess(np.max(np.abs(local.points[:3, 1])), 1e-12)
        np.testing.assert_allclose(published.points,
                                   local.points @ rotation.T + translation,
                                   rtol=0., atol=2e-15)

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

    def test_rejected_overlapping_move_rolls_back_against_original_global_floors(self):
        points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.],
                           [0., 0., -1.]])
        tetrahedra = np.array([[0, 1, 2, 3], [0, 2, 1, 4]])
        original = points.copy()
        floor = np.minimum(_tetra_quality(original, tetrahedra)[0], .5)
        candidate = points.copy()
        candidate[0] = [.99, .99, .99] # degrades/inverts shared incident cells
        accepted, _, _ = _transactional_quality_commit(
            points, candidate, np.array([0]), tetrahedra, np.array([0, 1]), floor)
        self.assertFalse(accepted)
        np.testing.assert_array_equal(points, original)
        improving = points.copy(); improving[0] = [.01, .01, 0.]
        accepted, _, _ = _transactional_quality_commit(
            points, improving, np.array([0]), tetrahedra, np.array([0, 1]), floor)
        self.assertTrue(accepted)
        np.testing.assert_array_equal(points[1:], original[1:])

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

    def test_saturated_vertex_roundoff_is_within_bound_but_true_overshoot_is_not(self):
        bound = .01875
        saturated = np.nextafter(bound, np.inf)
        self.assertGreater(saturated, bound)
        self.assertTrue(_within_displacement_bound([bound, saturated, 0.], bound))
        self.assertFalse(_within_displacement_bound(
            [bound * (1. + 10. * DISPLACEMENT_ROUNDOFF_TOLERANCE)], bound))
        self.assertFalse(_within_displacement_bound([np.nan], bound))
        self.assertFalse(_within_displacement_bound([np.inf], bound))
        # A clamped offset lands on the ball; recomputing its magnitude must not raise.
        directions = np.eye(3)
        offset = _bounded_offset(directions, np.array([.75, .75, .75]), bound)
        self.assertLessEqual(np.linalg.norm(offset),
                             bound * (1. + DISPLACEMENT_ROUNDOFF_TOLERANCE))
        self.assertGreaterEqual(np.linalg.norm(offset), bound * (1. - 1e-12))

    def test_saturated_corner_repair_does_not_raise_on_the_displacement_ball(self):
        # A tiny bound forces every movable vertex to saturate on the ball.
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0002)
        original = points.copy()
        bound = 1e-4
        with self.assertRaisesRegex(ValueError, "Semantic-corner quality repair failed"):
            _quality_repair(points, tetrahedra, node_supports, supports, recipe,
                            .01, 4., bound)
        displacement = np.linalg.norm(points - original, axis=1)
        self.assertTrue(_within_displacement_bound(displacement, bound))

    def test_true_overshoot_is_a_rejected_component_not_an_exception(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0002)
        original = points.copy()
        bound = .01875

        def overshooting_offset(directions, parameters, maximum_displacement):
            return 2. * maximum_displacement * directions[:, 0]

        with mock.patch.object(restore_planar_metric_mesh, "_bounded_offset",
                               overshooting_offset):
            with self.assertRaisesRegex(ValueError,
                                        "Semantic-corner quality repair failed"):
                _quality_repair(points, tetrahedra, node_supports, supports, recipe,
                                .01, 4., bound)
        # Every overshooting candidate was rejected: nothing moved, and the
        # final corner gate (not the bound check) reported the failure.
        np.testing.assert_array_equal(points, original)

    def test_corner_result_above_target_is_rejected_not_committed(self):
        # Initial corner aspect 3.88 is above the 3.8 repair target but within the
        # 4.0 gate, so a rejected repair still yields a report.
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.014)
        original = points.copy()
        before = float(_tetra_quality(points, tetrahedra)[1].max())
        self.assertGreater(before, 3.8); self.assertLess(before, 4.)

        def stalled_least_squares(residual, initial, **_):
            return SimpleNamespace(x=np.asarray(initial, dtype=float))

        with mock.patch.object(restore_planar_metric_mesh, "least_squares",
                               stalled_least_squares):
            report = _quality_repair(points, tetrahedra, node_supports, supports,
                                     recipe, .01, 4., .01875)
        self.assertEqual(report["RejectedCornerRepairs"], 1)
        self.assertEqual(report["CornerAspectsBefore"], [before])
        self.assertEqual(report["CornerAspectsAfter"], [before])
        self.assertEqual(report["QualityRepairVertices"], 0)
        np.testing.assert_array_equal(points, original)
        unpatched = _quality_repair(points, tetrahedra, node_supports, supports,
                                    recipe, .01, 4., .01875)
        self.assertEqual(unpatched["RejectedCornerRepairs"], 0)
        self.assertLessEqual(unpatched["CornerAspectsAfter"][0], 3.8 * (1. + 1e-6))


if __name__ == "__main__":
    unittest.main()
