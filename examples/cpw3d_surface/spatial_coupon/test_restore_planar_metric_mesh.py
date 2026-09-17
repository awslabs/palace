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
    _collapse_corner_ball_vertices,
    _movement_basis,
    _pinned_vertices,
    _quality_repair,
    _tetra_quality,
    _transactional_quality_commit,
    _within_displacement_bound,
    frozen_edge_layer_vertices,
    local_size_bounds,
    required_tetrahedron_vertices,
    restore_in_source_frame,
)
from mesh_array_io import read_medit_binary, read_mesh


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
        "TruePhysicalCorners": [[0., 0., 0.]], "CornerIsotropyRadius": .1}
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
            "TruePhysicalCorners": [[0., 0., 0.]],
            "CornerIsotropyRadius": .1
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
            "TruePhysicalCorners": [[0., 0., 0.]],
            "CornerIsotropyRadius": .1
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

    def test_stalled_corner_solver_is_a_rejected_repair_with_unchanged_points(self):
        # Initial corner aspect 3.88 is above the 3.8 repair target but within the
        # 4.0 gate; a solver that returns its initial point moves nothing, so the
        # repair is recorded as rejected and the report is still produced.
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
        self.assertEqual(report["CornerRepairOutcomes"][0]["Outcome"], "rejected")
        self.assertEqual(report["CornerAspectsBefore"], [before])
        self.assertEqual(report["CornerAspectsAfter"], [before])
        self.assertEqual(report["QualityRepairVertices"], 0)
        np.testing.assert_array_equal(points, original)
        unpatched = _quality_repair(points, tetrahedra, node_supports, supports,
                                    recipe, .01, 4., .01875)
        self.assertEqual(unpatched["RejectedCornerRepairs"], 0)
        self.assertEqual(unpatched["CornerRepairOutcomes"][0]["Outcome"], "target-reached")
        self.assertLessEqual(unpatched["CornerAspectsAfter"][0], 3.8 * (1. + 1e-6))

    @staticmethod
    def _nudging_least_squares(parameter_index, step, calls=1):
        """Solver stand-in: nudge one parameter on the first `calls` passes, then stall."""
        remaining = [calls]
        def nudge(residual, initial, **_):
            value = np.asarray(initial, dtype=float).copy()
            if remaining[0] > 0:
                remaining[0] -= 1
                value[parameter_index] += step
            return SimpleNamespace(x=value)
        return nudge

    def test_corner_result_within_gate_but_above_target_is_committed_and_reported(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0134)
        before = float(_tetra_quality(points, tetrahedra)[1].max())
        self.assertGreater(before, 4.)
        # Active vertices are [1, 2, 3] with 2 + 2 + 3 parameters; index 5 is the
        # apex y offset.  One nudge lands in (3.8, 4.0]: gate satisfied, target missed.
        with mock.patch.object(restore_planar_metric_mesh, "least_squares",
                               self._nudging_least_squares(5, .0006 / .01875)):
            report = _quality_repair(points, tetrahedra, node_supports, supports,
                                     recipe, .01, 4., .01875)
        after = report["CornerAspectsAfter"][0]
        self.assertGreater(after, 3.8); self.assertLessEqual(after, 4.)
        self.assertEqual(report["CornerRepairOutcomes"][0]["Outcome"],
                         "target-missed-gate-satisfied")
        self.assertEqual(report["CornerRepairsTargetMissedGateSatisfied"], 1)
        self.assertEqual(report["RejectedCornerRepairs"], 0)

    def test_corner_result_above_gate_is_rejected_and_the_final_gate_reports(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.01)
        original = points.copy()
        with mock.patch.object(restore_planar_metric_mesh, "least_squares",
                               self._nudging_least_squares(5, .0005 / .01875)):
            with self.assertRaisesRegex(ValueError, "Semantic-corner quality repair failed"):
                _quality_repair(points, tetrahedra, node_supports, supports,
                                recipe, .01, 4., .01875)
        np.testing.assert_array_equal(points, original)

    def test_corner_passes_alternate_one_ring_and_two_ring_inside_the_corner_ball(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0132)
        points = np.vstack((points, [[.03, .03, .03]]))
        tetrahedra = np.vstack((tetrahedra, [[1, 2, 3, 4]]))
        original = points.copy()
        parameter_counts = []

        def nudge(residual, initial, **_):
            parameter_counts.append(len(initial))
            value = np.asarray(initial, dtype=float).copy()
            value[5] += .0002 / .01875  # apex y offset; same slot in both rings
            return SimpleNamespace(x=value)

        with mock.patch.object(restore_planar_metric_mesh, "least_squares", nudge):
            report = _quality_repair(points, tetrahedra, node_supports, supports,
                                     recipe, .01, 4., .01875)
        # One ring: vertices 1, 2 (2 dof each) and 3 (3 dof); two ring adds vertex 4.
        self.assertEqual(parameter_counts, [7, 10, 7, 10])
        outcome = report["CornerRepairOutcomes"][0]
        self.assertEqual(outcome["Passes"], 4)
        self.assertEqual(outcome["Outcome"], "target-missed-gate-satisfied")
        self.assertLessEqual(report["CornerAspectsAfter"][0], 4.)
        self.assertLess(report["CornerAspectsAfter"][0], report["CornerAspectsBefore"][0])
        np.testing.assert_allclose(points[3, 1], original[3, 1] + 4 * .0002)
        # Vertex 4 at |x| = 0.052 falls outside a 0.051 ball; vertices 1, 2 (0.05) stay.
        out_of_ball = recipe | {"CornerIsotropyRadius": .051}
        parameter_counts.clear(); points[:] = original
        with mock.patch.object(restore_planar_metric_mesh, "least_squares", nudge):
            _quality_repair(points, tetrahedra, node_supports, supports,
                            out_of_ball, .01, 4., .01875)
        self.assertEqual(parameter_counts, [7, 7, 7, 7])
        # A 0.03 ball excludes the one-ring vertices 1, 2 (0.05) from the one-ring
        # pass too; only the apex (0.0195, 3 dof) moves in every pass.
        small_ball = recipe | {"CornerIsotropyRadius": .03}
        parameter_counts.clear(); points[:] = original

        def nudge_apex_only(residual, initial, **_):
            parameter_counts.append(len(initial))
            value = np.asarray(initial, dtype=float).copy()
            value[1] += .0002 / .01875  # apex y offset is the only 3-dof block
            return SimpleNamespace(x=value)

        with mock.patch.object(restore_planar_metric_mesh, "least_squares", nudge_apex_only):
            _quality_repair(points, tetrahedra, node_supports, supports,
                            small_ball, .01, 4., .01875)
        self.assertEqual(parameter_counts, [3, 3, 3, 3])

    def test_best_chained_candidate_is_committed_when_a_later_pass_regresses(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0134)
        original = points.copy()
        steps = iter([.0006 / .01875, -.0003 / .01875])

        def improve_then_regress(residual, initial, **_):
            value = np.asarray(initial, dtype=float).copy()
            value[5] += next(steps, 0.)
            return SimpleNamespace(x=value)

        with mock.patch.object(restore_planar_metric_mesh, "least_squares",
                               improve_then_regress):
            report = _quality_repair(points, tetrahedra, node_supports, supports,
                                     recipe, .01, 4., .01875)
        outcome = report["CornerRepairOutcomes"][0]
        # Pass 2 regressed (3.96 > 3.88) and stalled the chain; the pass-1
        # candidate is committed and AchievedAspect is its aspect, not the last one.
        self.assertEqual(outcome["Passes"], 2)
        self.assertEqual(outcome["Outcome"], "target-missed-gate-satisfied")
        np.testing.assert_allclose(points[3, 1], original[3, 1] + .0006)
        self.assertEqual(outcome["AchievedAspect"], report["CornerAspectsAfter"][0])
        self.assertAlmostEqual(report["CornerAspectsAfter"][0], 3.87936898, places=6)

    def test_chained_pass_initial_point_is_clipped_to_the_optimizer_bounds(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0002)
        real_least_squares = restore_planar_metric_mesh.least_squares
        initials = []

        def saturate_then_solve(residual, initial, **options):
            initials.append(np.asarray(initial, dtype=float).copy())
            if len(initials) == 1:
                value = np.zeros_like(initial); value[5] = .75  # apex y on the bound
                return SimpleNamespace(x=value)
            return real_least_squares(residual, initial, **options)

        def roundoff_offset(directions, parameters, maximum_displacement):
            # A few ulps of roundoff on a saturated offset, as re-deriving the
            # parameters of a clamped vertex can produce.
            return _bounded_offset(directions, parameters, maximum_displacement) * (
                1. + 4. * np.finfo(float).eps)

        with mock.patch.object(restore_planar_metric_mesh, "least_squares",
                               saturate_then_solve), \
                mock.patch.object(restore_planar_metric_mesh, "_bounded_offset",
                                  roundoff_offset):
            report = _quality_repair(points, tetrahedra, node_supports, supports,
                                     recipe, .01, 4., .01875)
        self.assertGreaterEqual(len(initials), 2)
        self.assertGreater(.75 * (1. + 4. * np.finfo(float).eps), .75)
        self.assertTrue(all(np.all(np.abs(initial) <= .75) for initial in initials[1:]))
        self.assertLessEqual(report["CornerAspectsAfter"][0], 4.)

    def test_corner_and_component_without_movable_vertices_are_recorded_rejections(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.014)
        original = points.copy()
        before = float(_tetra_quality(points, tetrahedra)[1].max())
        self.assertGreater(before, 3.8); self.assertLess(before, 4.)
        # Minimum scaled 0.4 makes the 0.70 cell a low-quality component (target
        # 0.8) while the final 0.4 gate still holds; every vertex is pinned.
        report = _quality_repair(points, tetrahedra, node_supports, supports, recipe,
                                 .4, 4., .01875, pinned_nodes=frozenset({0, 1, 2, 3}))
        outcome = report["CornerRepairOutcomes"][0]
        self.assertEqual(outcome, {"Passes": 0, "AchievedAspect": before,
                                   "Outcome": "rejected"})
        self.assertEqual(report["RejectedCornerRepairs"], 1)
        self.assertEqual(report["QualityRepairComponents"], 1)
        self.assertEqual(report["RejectedQualityRepairComponents"], 1)
        self.assertEqual(report["CornerAspectsAfter"], [before])
        self.assertEqual(report["QualityRepairVertices"], 0)
        np.testing.assert_array_equal(points, original)

    def test_pinned_vertices_are_matched_natively_and_never_move(self):
        points, tetrahedra, supports, node_supports, recipe = \
            _single_tetrahedron_repair_case(.0002)
        pinned = recipe | {"PinnedVertices": [{"Point": points[3].tolist(),
                                               "Kind": "reference-turn"}]}
        self.assertEqual(_pinned_vertices(points, pinned), frozenset({3}))
        with self.assertRaisesRegex(ValueError, "Pinned vertex is absent"):
            _pinned_vertices(points, recipe | {"PinnedVertices": [
                {"Point": [.5, .5, .5], "Kind": "reference-turn"}]})
        original = points.copy()
        try:
            report = _quality_repair(points, tetrahedra, node_supports, supports, pinned,
                                     .01, 4., .01875, pinned_nodes=frozenset({3}))
            self.assertEqual(report["PinnedVerticesFixed"], 1)
        except ValueError as error:
            self.assertIn("Semantic-corner quality repair failed", str(error))
        np.testing.assert_array_equal(points[3], original[3])
        with self.assertRaisesRegex(ValueError, "corner isotropy radius"):
            _quality_repair(points, tetrahedra, node_supports, supports,
                            {k: v for k, v in recipe.items() if k != "CornerIsotropyRadius"},
                            .01, 4., .01875)


def _corner_ball_with_inserted_vertex(height, ring=8, radius=.025):
    """A corner vertex at the origin, a ring of supported vertices at the isotropic
    size, a supported apex, and one free interior vertex on the axis at `height`
    (MMG's sub-hmin insertion when height < NormalSize)."""
    angles = 2. * np.pi * np.arange(ring) / ring
    points = np.vstack([[0., 0., 0.],
                        np.column_stack([radius * np.cos(angles), radius * np.sin(angles),
                                         np.full(ring, radius)]),
                        [0., 0., 2. * radius], [0., 0., height]])
    corner, apex, free = 0, ring + 1, ring + 2
    tetrahedra = []
    for i in range(ring):
        a, b = 1 + i, 1 + (i + 1) % ring
        tetrahedra.append([corner, a, b, free])
        tetrahedra.append([free, a, b, apex])
    tetrahedra = np.array(tetrahedra)
    _, _, determinant = _tetra_quality(points, tetrahedra)
    assert np.all(determinant > 0)
    supports = {100: {"Normal": [0., 0., 1.], "Offset": 0., "Attribute": 6001}}
    node_supports = {node: {100} for node in range(ring + 2)}
    recipe = {"NormalSize": radius, "CornerIsotropyRadius": 4. * radius,
              "TruePhysicalCorners": [[0., 0., 0.]],
              "SemanticContract": {"CutSurfaceRoles": [], "BoundaryLabels": []}}
    return points, tetrahedra, supports, node_supports, recipe, free


class CornerBallCollapseTest(unittest.TestCase):
    def test_sub_hmin_free_vertex_is_collapsed_onto_the_corner_with_valid_cavity(self):
        points, tetrahedra, supports, node_supports, recipe, free = \
            _corner_ball_with_inserted_vertex(.012)
        refs = np.arange(len(tetrahedra))
        before = float(_tetra_quality(points, tetrahedra[np.any(tetrahedra == 0, axis=1)])[1].max())
        new_points, new_tetrahedra, new_refs, vertex_map, collapsed = \
            _collapse_corner_ball_vertices(points, tetrahedra, refs, node_supports,
                                           frozenset(), recipe, .02)
        self.assertEqual(collapsed[0]["CollapsedVertices"], 1)
        self.assertLess(collapsed[0]["AspectAfter"], before)
        self.assertEqual(vertex_map[free], -1)
        self.assertEqual(len(new_points), len(points) - 1)
        # The cells joining the free vertex to the corner vanish; the others remap.
        self.assertEqual(len(new_tetrahedra), len(tetrahedra) // 2)
        self.assertTrue(np.all(new_refs % 2 == 1))
        scaled, _, determinant = _tetra_quality(new_points, new_tetrahedra)
        self.assertTrue(np.all(determinant > 0))
        self.assertGreaterEqual(float(scaled.min()), .02)
        # Every supported (boundary) vertex survives at its coordinates, unmoved.
        for node in node_supports:
            self.assertGreaterEqual(vertex_map[node], 0)
            np.testing.assert_array_equal(new_points[vertex_map[node]], points[node])
        # The corner keeps exactly the remapped cells; their aspect is the report's.
        incident = new_tetrahedra[np.any(new_tetrahedra == vertex_map[0], axis=1)]
        self.assertEqual(len(incident), len(tetrahedra) // 2)
        self.assertAlmostEqual(float(_tetra_quality(new_points, incident)[1].max()),
                               collapsed[0]["AspectAfter"])

    def test_vertex_at_or_above_hmin_and_supported_or_pinned_vertices_are_kept(self):
        for height, pinned, supported in ((.03, frozenset(), False),
                                          (.012, frozenset(), True),
                                          (.012, None, False)):
            points, tetrahedra, supports, node_supports, recipe, free = \
                _corner_ball_with_inserted_vertex(height)
            if supported:
                node_supports[free] = {100}
            pinned_nodes = frozenset({free}) if pinned is None else pinned
            new_points, new_tetrahedra, new_refs, vertex_map, collapsed = \
                _collapse_corner_ball_vertices(points, tetrahedra, np.arange(len(tetrahedra)),
                                               node_supports, pinned_nodes, recipe, .02)
            self.assertEqual(collapsed[0]["CollapsedVertices"], 0)
            np.testing.assert_array_equal(new_tetrahedra, tetrahedra)
            np.testing.assert_array_equal(vertex_map, np.arange(len(points)))
            self.assertEqual(len(new_points), len(points))

    def test_corner_ball_collapse_is_rotation_covariant_and_fails_closed(self):
        points, tetrahedra, supports, node_supports, recipe, free = \
            _corner_ball_with_inserted_vertex(.012)
        angle = .61
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                             [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
        tilt = np.array([[1., 0., 0.], [0., np.cos(.3), -np.sin(.3)],
                         [0., np.sin(.3), np.cos(.3)]])
        rotation = tilt @ rotation
        shift = np.array([2., -1., .5])
        moved_recipe = dict(recipe, TruePhysicalCorners=[(rotation @ [0., 0., 0.] + shift).tolist()])
        plain = _collapse_corner_ball_vertices(points, tetrahedra, np.arange(len(tetrahedra)),
                                               node_supports, frozenset(), recipe, .02)
        moved = _collapse_corner_ball_vertices(points @ rotation.T + shift, tetrahedra,
                                               np.arange(len(tetrahedra)), node_supports,
                                               frozenset(), moved_recipe, .02)
        np.testing.assert_array_equal(moved[1], plain[1])
        np.testing.assert_array_equal(moved[3], plain[3])
        self.assertEqual(moved[4][0]["CollapsedVertices"], plain[4][0]["CollapsedVertices"])
        np.testing.assert_allclose(moved[0], plain[0] @ rotation.T + shift, rtol=0., atol=1e-14)
        for missing in ("NormalSize", "CornerIsotropyRadius"):
            with self.assertRaises(ValueError):
                _collapse_corner_ball_vertices(points, tetrahedra, np.arange(len(tetrahedra)),
                                               node_supports, frozenset(),
                                               {k: v for k, v in recipe.items() if k != missing},
                                               .02)
        with self.assertRaisesRegex(ValueError, "absent"):
            _collapse_corner_ball_vertices(points, tetrahedra, np.arange(len(tetrahedra)),
                                           node_supports, frozenset(),
                                           dict(recipe, TruePhysicalCorners=[[1., 1., 1.]]), .02)


class LocalSizeBoundTest(unittest.TestCase):
    """Restoration bounds are relative to the local prescribed size at each vertex."""

    @staticmethod
    def recipe():
        return {"PhysicalSegments": [[-1., 0., 0., 1., 0., 0.]],
                "JunctionSegments": {"Segments": []}, "NormalSize": .025, "FarSize": .16,
                "RadialGrowth": 1., "ProtectedDistance": .05, "FarGrowth": .5,
                "TruePhysicalCorners": [[-1., 0., 0.]], "CornerIsotropyRadius": .1,
                "EdgeLayer": {"Spans": [[0., 0., 0., 1., 0., 0.]], "EdgeSize": .001,
                              "GrowthRatio": 2., "LayerThickness": .031}}

    def test_bounds_follow_the_local_size_and_layer_surface_vertices_are_frozen(self):
        recipe = self.recipe()
        points = np.array([[.5, 0., .001], [-.5, 0., .001], [.5, 0., .5], [.5, 0., .03],
                           [.5, 0., .033], [-.5, 0., .02]])
        # EdgeSize law on the span (2 nm); everywhere else the NormalSize cap (the
        # band law and the far size never tighten a bound below production's).
        np.testing.assert_allclose(local_size_bounds(points, recipe, .25),
                                   [.0005, .00625, .00625, .00625, .00625, .00625])
        with self.assertRaises(ValueError):
            local_size_bounds(points, recipe, 0.)
        node_supports = {0: {100}, 1: {100}, 2: {100}, 3: {100}, 4: {100}}
        # Supported vertices within LayerThickness + EdgeSize (32 nm) of a span are
        # frozen; the plain-band vertex and the vertex beyond the footprint are not,
        # nor is an interior (unsupported) vertex.
        self.assertEqual(frozen_edge_layer_vertices(points, node_supports, recipe), {0, 3})
        self.assertEqual(frozen_edge_layer_vertices(points, {5: {100}}, recipe), frozenset())
        plain = {k: v for k, v in recipe.items() if k != "EdgeLayer"}
        self.assertEqual(frozen_edge_layer_vertices(points, node_supports, plain), frozenset())
        np.testing.assert_allclose(local_size_bounds(points, plain, .75), [.75 * .025] * 6)

    def test_quality_repair_honors_per_vertex_bounds_and_fixes_frozen_vertices(self):
        points, tetrahedra, supports, node_supports, recipe = _single_tetrahedron_repair_case(.0002)
        original = points.copy()
        # A per-vertex bound: the apex may move up to 0.01875, the surface vertices
        # only 1e-6, so the sidewall vertices stay within their tiny ball.
        bound = np.array([1e-6, 1e-6, 1e-6, .01875])
        report = _quality_repair(points, tetrahedra, node_supports, supports, recipe,
                                 .01, 4., bound)
        displacement = np.linalg.norm(points - original, axis=1)
        self.assertTrue(np.all(displacement <= bound * (1. + DISPLACEMENT_ROUNDOFF_TOLERANCE)))
        self.assertGreater(displacement[3], 1e-6)
        self.assertEqual(report["QualityDisplacementBoundUm"]["MinimumUm"], 1e-6)
        self.assertEqual(report["QualityDisplacementBoundUm"]["MaximumUm"], .01875)
        self.assertLessEqual(report["MaximumFinalQualityDisplacementOverLocalBound"],
                             1. + DISPLACEMENT_ROUNDOFF_TOLERANCE)
        # Frozen edge-layer surface vertices never move, even with a large bound.
        points, tetrahedra, supports, node_supports, recipe = _single_tetrahedron_repair_case(.0002)
        original = points.copy()
        report = _quality_repair(points, tetrahedra, node_supports, supports, recipe,
                                 .01, 4., .01875, frozen_nodes=frozenset({1, 2}))
        np.testing.assert_array_equal(points[[0, 1, 2]], original[[0, 1, 2]])
        self.assertEqual(report["FrozenEdgeLayerVertices"], 2)
        self.assertEqual(report["QualityDisplacementBoundUm"], .01875)
        with self.assertRaises(ValueError):
            _quality_repair(points, tetrahedra, node_supports, supports, recipe, .01, 4.,
                            np.array([.01, .01]))

    def test_corner_collapse_threshold_is_the_local_size(self):
        points, tetrahedra, supports, node_supports, recipe, free = \
            _corner_ball_with_inserted_vertex(.012)
        refs = np.arange(len(tetrahedra))
        # With a local size below the vertex's shortest edge nothing is collapsed;
        # with the recipe NormalSize (default) the sub-hmin vertex is.
        small = np.full(len(points), .005)
        kept = _collapse_corner_ball_vertices(points, tetrahedra, refs, node_supports,
                                              frozenset(), recipe, .02, small)
        self.assertEqual(kept[4][0]["CollapsedVertices"], 0)
        collapsed = _collapse_corner_ball_vertices(points, tetrahedra, refs, node_supports,
                                                   frozenset(), recipe, .02)
        self.assertEqual(collapsed[4][0]["CollapsedVertices"], 1)
        with self.assertRaises(ValueError):
            _collapse_corner_ball_vertices(points, tetrahedra, refs, node_supports, frozenset(),
                                           recipe, .02, small[:-1])


if __name__ == "__main__":
    unittest.main()

class RequiredTetrahedraTest(unittest.TestCase):
    """MMG's kept-verbatim seed cells are read from the native output and frozen."""

    @staticmethod
    def write_medit_binary(path, points, triangles, tetrahedra, required, extra_keyword=True):
        """A version-2 binary Medit mesh in MMG's keyword order with a keyword the
        tools do not consume (Corners) before the cells and RequiredTetrahedra after."""
        import struct
        chunks = [struct.pack("<ii", 1, 2), struct.pack("<i", 3)]
        position = 8 + 4 + 4 + 4
        chunks.append(struct.pack("<ii", position, 3))
        def field(code, payload):
            nonlocal position
            header_size = 4 + 4
            position += header_size + len(payload)
            chunks.append(struct.pack("<ii", code, position)); chunks.append(payload)
        field(4, struct.pack("<i", len(points)) + b"".join(
            struct.pack("<dddi", *point, 0) for point in points))
        if extra_keyword:
            field(13, struct.pack("<ii", 1, 1))
        field(8, struct.pack("<i", len(tetrahedra)) + b"".join(
            struct.pack("<iiiii", *(v + 1 for v in cell), 7) for cell in tetrahedra))
        if required is not None:
            field(12, struct.pack("<i", len(required)) + b"".join(
                struct.pack("<i", index + 1) for index in required))
        field(6, struct.pack("<i", len(triangles)) + b"".join(
            struct.pack("<iiii", *(v + 1 for v in tri), 100) for tri in triangles))
        chunks.append(struct.pack("<ii", 54, 0))
        path.write_bytes(b"".join(chunks))

    def test_native_reader_skips_unknown_keywords_and_flags_required_cells(self):
        import tempfile
        from pathlib import Path
        points, tetrahedra, supports, node_supports, recipe, free = \
            _corner_ball_with_inserted_vertex(.012)
        triangles = np.array([[0, 1, 2], [1, 2, 9]])
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "adapted.meshb"
            self.write_medit_binary(path, points, triangles, tetrahedra, [1, 3])
            mesh = read_mesh(path)
            self.assertEqual([block.type for block in mesh.cells], ["tetra", "triangle"])
            np.testing.assert_array_equal(mesh.points, points)
            np.testing.assert_array_equal(mesh.cells[0].data, tetrahedra)
            np.testing.assert_array_equal(mesh.cells[1].data, triangles)
            np.testing.assert_array_equal(mesh.cell_data["medit:ref"][1], [100, 100])
            flags = mesh.cell_data["medit:required"][0]
            self.assertEqual(flags.tolist(), [0, 1, 0, 1] + [0] * (len(tetrahedra) - 4))
            nodes, count = required_tetrahedron_vertices(mesh)
            self.assertEqual(count, 2)
            self.assertEqual(nodes, frozenset(int(v) for v in tetrahedra[[1, 3]].ravel()))
            # meshio's own reader mis-parses the RequiredTetrahedra keyword; ours is exact.
            self.write_medit_binary(path.with_name("plain.meshb"), points, triangles, tetrahedra, None)
            plain = read_medit_binary(path.with_name("plain.meshb"))
            self.assertEqual(int(plain.cell_data["medit:required"][0].sum()), 0)
            self.assertEqual(required_tetrahedron_vertices(plain), (frozenset(), 0))
            self.write_medit_binary(path.with_name("bad.meshb"), points, triangles, tetrahedra,
                                    [len(tetrahedra)])
            with self.assertRaisesRegex(ValueError, "out of range"):
                read_medit_binary(path.with_name("bad.meshb"))

    def test_required_vertices_are_neither_collapsed_nor_moved(self):
        points, tetrahedra, supports, node_supports, recipe, free = \
            _corner_ball_with_inserted_vertex(.012)
        refs = np.arange(len(tetrahedra))
        # The sub-hmin free vertex belongs to a required cell: no collapse.
        _, new_tetrahedra, _, vertex_map, collapsed = _collapse_corner_ball_vertices(
            points, tetrahedra, refs, node_supports, frozenset(), recipe, .02,
            frozen_nodes=frozenset({free}))
        self.assertEqual(collapsed[0]["CollapsedVertices"], 0)
        np.testing.assert_array_equal(new_tetrahedra, tetrahedra)
        # Through the restoration: the required cells' vertices are frozen in the
        # repair and reported; without the flag the same vertex is collapsed.
        triangles = np.array([[0, 1, 2]])
        normal = np.cross(points[1] - points[0], points[2] - points[0])
        normal /= np.linalg.norm(normal)
        plane = {"Normal": normal.tolist(), "Offset": float(normal @ points[0]), "Attribute": 6001}
        recipe = {**recipe, "PlanarSupports": {"100": plane}, "FarSize": .1,
                  "PhysicalSegments": [], "PinnedVertices": [],
                  "SemanticContract": {"CutSurfaceRoles": [], "BoundaryLabels": [],
                                       "VolumeMaterials": [{"Material": "vacuum", "Attribute": 7}]}}
        def mesh(required):
            flags = np.zeros(len(tetrahedra), dtype=np.int32); flags[list(required)] = 1
            return meshio.Mesh(points.copy(), [("triangle", triangles), ("tetra", tetrahedra)],
                               cell_data={"medit:ref": [np.array([100]), np.full(len(tetrahedra), 7)],
                                          "medit:required": [np.zeros(1, dtype=np.int32), flags]})
        _, _, report = restore_in_source_frame(mesh([0]), recipe, (.25, "local"), .01, 20., (.75, "local"))
        self.assertEqual(report["RequiredTetrahedra"], 1)
        self.assertEqual(report["RequiredVertices"], 4)
        self.assertEqual(report["CollapsedCornerVertices"], 0)
        self.assertLess(report["MaximumRequiredVertexCorrectionUm"], 1e-15)  # plane roundoff only
        _, _, report = restore_in_source_frame(mesh([]), recipe, (.25, "local"), .01, 20., (.75, "local"))
        self.assertEqual(report["RequiredTetrahedra"], 0)
        self.assertGreater(report["CollapsedCornerVertices"], 0)
