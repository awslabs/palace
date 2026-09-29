#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The corner coupon's trace basis rule (supervisor decision on the corner-family review,
2026-09-29): on the rings that meet the metal the knots are the two metal-arm crossings and
the metal-interior knot (PEC) plus five free knots at equal fractions of the free arc, so
that every node of the angle family has the same zero set and like-to-like free knots; the
box corners that are no knot are slave vertices of the trace triangulation; no free hat has
support on the PEC part of the box contour. The 90-degree convex node reproduces the lane-2
8-knot layout byte-identically."""

import importlib.util
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GENERATOR = load("generate_corner_response")
RADIUS, THICKNESS, OVERETCH = 1.9, 0.1, 0.05
ANGLES = (75.0, 82.5, 90.0, 105.0, 120.0, 135.0, 142.5, 150.0, 165.0, 172.5, 180.0)


def lane2_surface():
    """The lane-2 (angle-independent) layout: every ring the fixed knots, no slaves."""
    return GENERATOR.build_surface(RADIUS, 8, THICKNESS, OVERETCH)


def pec_mask(angle, topology):
    return lambda points: GENERATOR.pec_trace_mask(
        points, RADIUS, angle, 0.0, THICKNESS, 90.0, topology
    )


class TraceBasisRuleTest(unittest.TestCase):
    def test_knot_semantics_do_not_depend_on_the_angle(self):
        for topology in ("convex", "concave"):
            zero_sets = set()
            for angle in ANGLES:
                surface = GENERATOR.build_surface(
                    RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology
                )
                self.assertEqual(surface.basis_size, 72)
                self.assertEqual(surface.contour_groups, [8] * 9)
                zero_sets.add(tuple(surface.zero_trace_indices()))
                points = np.asarray(surface.knot_points)
                # The rule's zero set is exactly the set of knots on the PEC footprint.
                np.testing.assert_array_equal(
                    np.asarray(surface.knot_zero), pec_mask(angle, topology)(points)
                )
                # Rings that do not meet the metal keep the fixed layout.
                fixed = lane2_surface()
                for ring in (0, 1, 2, 5, 6, 7, 8):
                    np.testing.assert_array_equal(
                        points[8 * ring:8 * ring + 8],
                        np.asarray(fixed.knot_points)[8 * ring:8 * ring + 8],
                    )
                # The trace surface lies on the box: every vertex has |x| or |y| = R
                # (outer rings) or lies on an inner cap ring.
                vertices = surface.vertex_points()
                on_box = np.isclose(np.abs(vertices[:, :2]).max(axis=1), RADIUS, atol=1.0e-9)
                inner = np.isclose(np.abs(vertices[:, 2]), RADIUS, atol=1.0e-12)
                self.assertTrue(np.all(on_box | inner))
                # The hats form a partition of unity at every vertex (slaves included).
                hats = np.array([surface.hat_values(k) for k in range(surface.basis_size)])
                np.testing.assert_allclose(hats.sum(axis=0), 1.0)
                # Every triangle is nondegenerate and every vertex is used.
                triangles = surface.triangles
                areas = 0.5 * np.linalg.norm(
                    np.cross(
                        vertices[triangles[:, 1]] - vertices[triangles[:, 0]],
                        vertices[triangles[:, 2]] - vertices[triangles[:, 0]],
                    ),
                    axis=1,
                )
                self.assertGreater(areas.min(), 1.0e-6)
                self.assertEqual(len(np.unique(triangles)), len(vertices))
            self.assertEqual(len(zero_sets), 1, topology)
        self.assertEqual(
            GENERATOR.zero_slots("convex"), [4, 5, 6]
        )  # crossing 1, metal interior, crossing 2 (0 / 45 / 90 deg of the 90 node)
        self.assertEqual(GENERATOR.zero_slots("concave"), [0, 6, 7])

    def test_ninety_degree_convex_node_is_the_lane2_layout(self):
        surface = GENERATOR.build_surface(
            RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=90.0, topology="convex"
        )
        fixed = lane2_surface()
        np.testing.assert_array_equal(
            np.asarray(surface.knot_points), np.asarray(fixed.knot_points)
        )
        np.testing.assert_array_equal(surface.triangles, fixed.triangles)
        self.assertEqual(surface.slaves, [])
        self.assertEqual(
            [index + 1 for index in surface.zero_trace_indices()],
            [29, 30, 31, 37, 38, 39],
        )

    def test_free_hats_have_no_support_on_the_pec_contour(self):
        for topology in ("convex", "concave"):
            for angle in ANGLES:
                surface = GENERATOR.build_surface(
                    RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology
                )
                self.assertEqual(
                    GENERATOR.free_hat_pec_support(
                        np.asarray(surface.knot_points),
                        surface.contour_groups,
                        surface.zero_trace_indices(),
                        pec_mask(angle, topology),
                        surface.slaves,
                    ),
                    [],
                    (topology, angle),
                )

    def test_gate_reports_the_lane2_layout_defect(self):
        # The corner-family review's root cause: with the angle-independent 8-knot rings the
        # second arm of a 165-degree convex corner crosses the left side between the free
        # knot (-R, 0) (1-based 25) and the PEC knot (-R, R) (32); 73 % of that segment is
        # PEC. 90 / 135 / 180 degrees pass (the arms lie on knot rays).
        fixed = lane2_surface()
        points = np.asarray(fixed.knot_points)
        for angle, expected in ((165.0, [(24, 31, 0.73), (32, 39, 0.73)]),
                                (120.0, [(31, 30, 0.58), (39, 38, 0.58)])):
            mask = pec_mask(angle, "convex")
            offending = GENERATOR.free_hat_pec_support(
                points, fixed.contour_groups, np.flatnonzero(mask(points)), mask
            )
            self.assertEqual(
                [(free, fixed_knot, round(fraction, 2)) for free, fixed_knot, fraction in offending],
                expected,
            )
        for angle in (90.0, 135.0, 180.0):
            mask = pec_mask(angle, "convex")
            self.assertEqual(
                GENERATOR.free_hat_pec_support(
                    points, fixed.contour_groups, np.flatnonzero(mask(points)), mask
                ),
                [],
            )

    def test_connect_rings_by_fraction_reproduces_connect_rings(self):
        fractions = np.arange(8) / 8.0
        expected = []
        for index in range(8):  # the lane-2 connect_rings
            following = (index + 1) % 8
            expected.append((index, following, 8 + following))
            expected.append((index, 8 + following, 8 + index))
        actual = []
        GENERATOR.connect_rings_by_fraction(actual, 0, fractions, 8, fractions)
        self.assertEqual(actual, expected)
        # One extra vertex on the second ring adds one triangle to its quad.
        actual = []
        GENERATOR.connect_rings_by_fraction(
            actual, 0, fractions, 8, np.sort(np.append(fractions, 0.3))
        )
        self.assertEqual(len(actual), 17)
        self.assertEqual(
            sorted(set(v for t in actual for v in t)), list(range(17))
        )

    def test_arm_crossings(self):
        first, second = GENERATOR.arm_crossing_fractions(RADIUS, 90.0)
        self.assertAlmostEqual(first, 0.5)
        self.assertAlmostEqual(second, 0.75)
        first, second = GENERATOR.arm_crossing_fractions(RADIUS, 180.0)
        self.assertAlmostEqual(second, 1.0)
        first, second = GENERATOR.arm_crossing_fractions(RADIUS, 165.0)
        point = GENERATOR.square_perimeter_point(RADIUS, 0.0, second)
        self.assertAlmostEqual(point[0], -RADIUS)
        self.assertAlmostEqual(point[1], RADIUS * np.tan(np.deg2rad(15.0)))


if __name__ == "__main__":
    unittest.main()
