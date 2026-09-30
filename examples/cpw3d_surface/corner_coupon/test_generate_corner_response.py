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

    def test_ring_points_agree_with_square_perimeter_point(self):
        """Review m6: the generator places fixed-fraction knots with square_ring's
        coordinates, the C++ runtime (SquarePerimeterPoint, the square_perimeter_point
        formula) by the fraction alone; the two agree to floating-point rounding (~1e-15 R
        of either arithmetic; MatchCornerFamily compares at 1e-9 R) at the fixed fractions
        and use the same formula elsewhere. The eight fixed points are the table
        pinned by the C++ unit test SurfaceResponseOperatorCornerTraceBasis."""
        R = RADIUS
        table = [(-R, 0.0), (-R, -R), (0.0, -R), (R, -R), (R, 0.0), (R, R), (0.0, R), (-R, R)]
        fixed = GENERATOR.square_ring(R, THICKNESS, 8)
        for k, (x, y) in enumerate(table):
            np.testing.assert_allclose(fixed[k], (x, y, THICKNESS), rtol=0, atol=1.0e-14 * R)
            np.testing.assert_allclose(
                GENERATOR.square_perimeter_point(R, THICKNESS, k / 8),
                (x, y, THICKNESS),
                rtol=0,
                atol=1.0e-14 * R,
            )
        for topology in ("convex", "concave"):
            for angle in ANGLES + (153.435, 100.0, 170.0):
                layout = GENERATOR.metal_ring_layout(R, angle, topology, 8)
                points = GENERATOR.ring_points(R, 0.0, layout, 8)
                for (fraction, _, _), point in zip(layout, points):
                    expected = GENERATOR.square_perimeter_point(R, 0.0, fraction)
                    self.assertLessEqual(
                        max(abs(a - b) for a, b in zip(point, expected)),
                        GENERATOR.FIXED_FRACTION_IDENTITY * 8.0 * R,
                    )

    def test_knot_coincidence_snaps_a_free_knot_onto_a_box_corner(self):
        """Review m3: a free knot within KNOT_COINCIDENCE_FRACTION of a box corner's
        fraction takes it exactly (no slave there); outside the band the corner keeps its
        slave at least 8e-6 R from the knot. Convex free 2 sits at (-R, -R) (fraction 1/8)
        when the second arm crosses the left side at fraction 15/16 (153.435 degrees); the
        arm's crossing moves by 8 R per unit fraction, free 2 by two thirds of that. The
        C++ rule (kKnotCoincidenceFraction, cornertracebasis.cpp) uses the same value."""
        R = RADIUS
        self.assertEqual(GENERATOR.KNOT_COINCIDENCE_FRACTION, 1.0e-6)
        # The cross-language check is part of the contract (merge review m-b): the C++
        # source must be present, a missing file fails instead of skipping silently.
        source = (ROOT / "../../../palace/models/cornertracebasis.cpp").resolve()
        self.assertTrue(source.is_file(), f"C++ rule source missing: {source}")
        line = [
            text for text in source.read_text().splitlines()
            if text.startswith("constexpr double kKnotCoincidenceFraction")
        ]
        self.assertEqual(len(line), 1)
        self.assertEqual(float(line[0].split("=")[1].strip(" ;")), 1.0e-6)
        for offset, snapped in ((7.5e-7, True), (3.0e-6, False)):
            angle = np.degrees(np.arctan2(0.5 * R - 8.0 * R * offset, -R))
            layout = GENERATOR.metal_ring_layout(R, angle, "convex", 8)
            slaves = [fraction for fraction, kind, _ in layout if kind == "slave"]
            free2 = [fraction for fraction, kind, slot in layout if slot == 0]
            self.assertEqual(len(free2), 1)
            self.assertEqual(len(slaves), 3 if snapped else 4)
            points = GENERATOR.ring_points(R, 0.0, layout, 8)
            (point,) = [
                tuple(p) for (_, _, slot), p in zip(layout, points) if slot == 0
            ]
            if snapped:
                self.assertEqual(free2[0], 0.125)
                self.assertNotIn(0.125, slaves)
                self.assertEqual(point, (-R, -R, 0.0))
            else:
                self.assertIn(0.125, slaves)
                self.assertAlmostEqual(free2[0], 0.125 + offset * 2.0 / 3.0, places=12)
                self.assertEqual(point[1], -R)
                self.assertGreaterEqual(point[0] - (-R), 8.0e-6 * R)
            self.assertEqual(
                GENERATOR.free_hat_pec_support(
                    *self._surface_arguments(angle, "convex")
                ),
                [],
            )

    def test_connectivity_angle_fixes_the_band_triangulation_over_a_segment(self):
        """Corner-qualification block 2026-09-29: with a segment connectivity angle the
        bands next to the metal rings are merged in the order of the rule's layout at that
        angle, so the triangulation is the same at every angle of the segment (here across
        the convex free-1 passage of the side midpoint (-R, 0) at 141.34 deg, where the
        perimeter-ordered merge flips a quad diagonal), the knots stay at the angle's own
        positions, and the slave table is unchanged. Without one the surface is byte-identical
        to the perimeter-ordered merge (the recorded coupons). A connectivity angle across a
        knot-corner passage (135 deg, the crossing at (-R, R)) is refused: the band triangles
        would fold."""
        keyed = {}
        for angle in (137.0, 141.0, 142.0, 150.0, 153.0):
            surface = GENERATOR.build_surface(
                RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=angle, topology="convex",
                connectivity_angle_degrees=144.2,
            )
            plain = GENERATOR.build_surface(
                RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=angle, topology="convex"
            )
            np.testing.assert_array_equal(surface.knot_points, plain.knot_points)
            self.assertEqual(surface.slaves, plain.slaves)
            keyed[angle] = set(map(tuple, np.sort(surface.triangles, axis=1).tolist()))
            self.assertEqual(
                GENERATOR.free_hat_pec_support(
                    np.asarray(surface.knot_points), surface.contour_groups,
                    surface.zero_trace_indices(), pec_mask(angle, "convex"), surface.slaves
                ),
                [],
            )
            plain_set = set(map(tuple, np.sort(plain.triangles, axis=1).tolist()))
            # Perimeter order = the keyed order once the knot has passed the midpoint.
            self.assertEqual(plain_set == keyed[angle], angle > 141.3402)
        self.assertEqual(len({frozenset(t) for t in keyed.values()}), 1)
        with self.assertRaises(ValueError):
            GENERATOR.build_surface(
                RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=130.0, topology="convex",
                connectivity_angle_degrees=144.2,
            )
        rule = GENERATOR.trace_basis_rule(8, 144.2)
        self.assertEqual(rule["ConnectivityAngleDegrees"], 144.2)
        self.assertNotIn("ConnectivityAngleDegrees", GENERATOR.trace_basis_rule(8))

    def _surface_arguments(self, angle, topology):
        surface = GENERATOR.build_surface(
            RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology
        )
        return (
            np.asarray(surface.knot_points),
            surface.contour_groups,
            surface.zero_trace_indices(),
            pec_mask(angle, topology),
            surface.slaves,
        )


class HeldoutPotentialTest(unittest.TestCase):
    """USER decision 149 (6), option (c): the held-out self-check potential excites the free
    knots; its cutoff vanishes only on the PEC part of the metal rings."""

    def metal_ring_free_knots(self, surface):
        points = np.asarray(surface.knot_points)
        on_metal_ring = np.any(
            np.abs(points[:, 2][:, None] - np.asarray([0.0, THICKNESS])) <= 1.0e-12, axis=1
        ) & (np.max(np.abs(points[:, :2]), axis=1) >= RADIUS - 1.0e-12)
        return on_metal_ring & ~np.asarray(surface.knot_zero)

    def test_perimeter_arc_distance(self):
        arc = (0.5, 0.75)
        self.assertEqual(GENERATOR.perimeter_arc_distance(0.6, arc), 0.0)
        self.assertEqual(GENERATOR.perimeter_arc_distance(0.5, arc), 0.0)
        self.assertEqual(GENERATOR.perimeter_arc_distance(0.75 + 1.0e-13, arc), 0.0)
        self.assertEqual(GENERATOR.perimeter_arc_distance(0.5 - 1.0e-13, arc), 0.0)
        self.assertAlmostEqual(GENERATOR.perimeter_arc_distance(0.85, arc), 0.1)
        self.assertAlmostEqual(GENERATOR.perimeter_arc_distance(0.4, arc), 0.1)
        self.assertAlmostEqual(GENERATOR.perimeter_arc_distance(0.125, arc), 0.375)
        wrapped = (0.9, 1.2)  # the concave metal arc wraps through fraction 0
        self.assertEqual(GENERATOR.perimeter_arc_distance(0.1, wrapped), 0.0)
        self.assertAlmostEqual(GENERATOR.perimeter_arc_distance(0.3, wrapped), 0.1)

    def test_heldout_ring_layout_adds_the_crossings_as_pec_vertices(self):
        for topology in ("convex", "concave"):
            for angle in ANGLES:
                layout = GENERATOR.heldout_ring_layout(RADIUS, angle, topology, 32)
                fractions = [fraction for fraction, _, _ in layout]
                self.assertEqual(fractions, sorted(fractions))
                self.assertEqual([slot for _, _, slot in layout], list(range(len(layout))))
                first, second = GENERATOR.arm_crossing_fractions(RADIUS, angle)
                for crossing in (first, second % 1.0):
                    self.assertTrue(any(abs(f - crossing) <= 1.0e-12 for f in fractions))
                # 90 / 135 / 180 degrees: the arms lie on fixed rays; other angles add
                # the second crossing.
                self.assertEqual(len(layout), 32 if angle in (90.0, 135.0, 180.0) else 33)
                points = GENERATOR.ring_points(RADIUS, 0.0, layout, 32)
                np.testing.assert_array_equal(
                    np.asarray([kind == "zero" for _, kind, _ in layout]),
                    pec_mask(angle, topology)(points),
                )

    def test_heldout_trace_vanishes_on_the_pec_part_only(self):
        for topology in ("convex", "concave"):
            for angle in ANGLES:
                fine = GENERATOR.build_surface(
                    RADIUS, 32, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology,
                    cap_centers=True, crossing_vertices=True,
                )
                self.assertEqual(fine.slaves[-2:], [((0.0, 0.0, RADIUS), None, None, None),
                                                    ((0.0, 0.0, -RADIUS), None, None, None)])
                points = fine.vertex_points()
                pec = pec_mask(angle, topology)(points)
                self.assertGreater(np.count_nonzero(pec), 0)
                values = GENERATOR.heldout_potential(points, RADIUS, THICKNESS, angle, topology)
                cutoff = GENERATOR.heldout_cutoff(points, RADIUS, angle, topology, THICKNESS)
                np.testing.assert_array_equal(values[pec], 0.0)
                self.assertTrue(np.all(cutoff[~pec] > 0.0), (topology, angle))
                self.assertTrue(np.all(cutoff <= 1.0))
                # Far from the metal band (the outer rings at z = -/+ R and the caps) the
                # potential is the bare polynomial.
                far = np.abs(np.abs(points[:, 2]) - RADIUS) <= 1.0e-12
                np.testing.assert_allclose(cutoff[far], 1.0)
                # Every fine vertex of the metal rings is PEC or free by the same rule as the
                # coupon's PEC mask, so the piecewise-linear trace is zero on the whole PEC
                # part of the box (its boundary, the crossings, is a vertex) and nonzero on
                # every triangle touching a free vertex.
                on_metal_ring = np.any(
                    np.abs(points[:, 2][:, None] - np.asarray([0.0, THICKNESS])) <= 1.0e-12,
                    axis=1,
                )
                self.assertEqual(np.count_nonzero(on_metal_ring), 2 * (32 if angle in (90.0, 135.0, 180.0) else 33))

    def test_heldout_coefficients_excite_every_free_knot(self):
        for topology in ("convex", "concave"):
            for angle in ANGLES:
                surface = GENERATOR.build_surface(
                    RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology
                )
                points = np.asarray(surface.knot_points)
                zero = np.asarray(surface.knot_zero)
                coefficients = GENERATOR.heldout_potential(points, RADIUS, THICKNESS, angle, topology)
                cutoff = GENERATOR.heldout_cutoff(points, RADIUS, angle, topology, THICKNESS)
                np.testing.assert_array_equal(coefficients[zero], 0.0)
                self.assertTrue(np.all(cutoff[~zero] > 0.0), (topology, angle))
                free_on_metal_rings = self.metal_ring_free_knots(surface)
                self.assertEqual(np.count_nonzero(free_on_metal_rings), 2 * GENERATOR.FREE_KNOTS)
                # The blind spot of the former band cutoff: exactly zero on both metal rings.
                legacy = GENERATOR.metal_band_cutoff(points, RADIUS, THICKNESS)
                np.testing.assert_array_equal(legacy[free_on_metal_rings], 0.0)
                # Option (c): every free knot of the metal rings carries a nonzero
                # coefficient; the knot nearest a crossing (a sixth of the free arc away, at
                # least 0.21 R here) is at least a third excited.
                self.assertGreater(np.min(cutoff[free_on_metal_rings]), 1.0 / 3.0, (topology, angle))
                self.assertTrue(np.all(coefficients[free_on_metal_rings] != 0.0))

    def test_crossing_vertices_need_an_angle_and_no_connectivity(self):
        with self.assertRaises(ValueError):
            GENERATOR.build_surface(RADIUS, 32, THICKNESS, OVERETCH, crossing_vertices=True)
        with self.assertRaises(ValueError):
            GENERATOR.build_surface(
                RADIUS, 32, THICKNESS, OVERETCH, angle_degrees=120.0, topology="convex",
                connectivity_angle_degrees=112.5, crossing_vertices=True,
            )
        # The straight-arm form of the 2D generators: the cutoff at a point of the free arc
        # depends on its distance to the nearest crossing and to the band.
        points = np.asarray([[RADIUS, -RADIUS / 3.0, 0.0], [RADIUS, 0.0, RADIUS / 3.0]])
        cutoff = GENERATOR.heldout_cutoff(points, RADIUS, 90.0, "convex", THICKNESS)
        np.testing.assert_allclose(cutoff, [1.0, (1.0 - 3.0 * THICKNESS / RADIUS) ** 2 * (3.0 - 2.0 * (1.0 - 3.0 * THICKNESS / RADIUS))])


class RefinedRuleTest(unittest.TestCase):
    """The AllRingsFollowMetal rule (corner-basis refinement, USER decision 161 (1),
    2026-09-30): every ring follows the metal fractions (5 metal-interior + 9 graded free
    knots, the extra ring at t + d, cap centre slaves), PEC on the metal rings only, no events;
    the held-out reference surface decoupled from the basis levels. Pinned to the C++ rule
    (test/unit/test-cornerbasisrefinement.cpp) by the same fractions."""

    RULE = GENERATOR.REFINED_RULE

    def test_rule_record_and_parameters(self):
        rule = self.RULE
        self.assertEqual((rule.layout, rule.metal_interior_knots, rule.free_knots, rule.ring_size),
                         ("AllRingsFollowMetal", 5, 9, 16))
        self.assertEqual(rule.free_knot_grading, (1.0 / 3.0, 2.0 / 3.0))
        record = GENERATOR.trace_basis_rule(16, None, rule)
        self.assertEqual(record["RingLayout"], "AllRingsFollowMetal")
        self.assertEqual(record["RingSize"], 16)
        self.assertEqual(record["FreeKnotGrading"], [1.0 / 3.0, 2.0 / 3.0])
        self.assertEqual(GENERATOR.TraceBasisRule.from_record(record), rule)
        self.assertEqual(GENERATOR.TraceBasisRule.from_record(GENERATOR.trace_basis_rule(8)), GENERATOR.LEGACY_RULE)
        with self.assertRaises(ValueError):
            GENERATOR.trace_basis_rule(16, 112.5, rule)  # no events: no connectivity angle
        with self.assertRaises(ValueError):
            GENERATOR.trace_basis_rule(8, None, rule)
        with self.assertRaises(ValueError):
            GENERATOR.TraceBasisRule("MetalRingsOnly", 1, 5, (1.0 / 3.0,))
        with self.assertRaises(ValueError):
            GENERATOR.TraceBasisRule("AllRingsFollowMetal", 5, 3, (1.0 / 3.0, 2.0 / 3.0))
        self.assertEqual(
            [round(z, 6) for z in rule.levels(RADIUS, THICKNESS, OVERETCH)],
            [round(z, 6) for z in (-RADIUS, -RADIUS / 3.0, -OVERETCH, 0.0, THICKNESS, THICKNESS + OVERETCH, RADIUS / 3.0, RADIUS)],
        )

    def test_layout_pins_shared_with_the_cpp_rule(self):
        # The same numbers as CornerRefinedRuleLayout in test-cornerbasisrefinement.cpp.
        pins = {
            ("convex", 120.0): ([0.905502116982, 0.990696208596, 0.07589030021, 0.161084391824,
                                 0.246278483438, 0.331472575053, 0.416666666667, 0.458333333333,
                                 0.5, 0.553694797275, 0.60738959455, 0.661084391824,
                                 0.714779189099, 0.768473986374, 0.822168783649, 0.863835450315],
                                [8, 9, 10, 11, 12, 13, 14]),
            ("concave", 105.0): ([0.5, 0.541666666667, 0.583333333333, 0.602804497065,
                                  0.622275660796, 0.641746824527, 0.661217988258, 0.680689151989,
                                  0.700160315721, 0.741826982387, 0.783493649054, 0.902911374212,
                                  0.022329099369, 0.141746824527, 0.261164549685, 0.380582274842],
                                 [0, 10, 11, 12, 13, 14, 15]),
        }
        for (topology, angle), (fractions, zero_slots) in pins.items():
            layout = GENERATOR.rule_ring_layout(RADIUS, angle, topology, self.RULE, True)
            by_slot = {slot: (f, kind) for f, kind, slot in layout if kind != "slave"}
            self.assertEqual(len(by_slot), 16)
            self.assertEqual(sum(1 for _, kind, _ in layout if kind == "slave"), 4)
            for slot in range(16):
                self.assertAlmostEqual(by_slot[slot][0], fractions[slot], places=9, msg=(topology, angle, slot))
                self.assertEqual(by_slot[slot][1] == "zero", slot in zero_slots, (topology, angle, slot))
            self.assertEqual(GENERATOR.zero_slots(topology, self.RULE), zero_slots)
            off_metal = GENERATOR.rule_ring_layout(RADIUS, angle, topology, self.RULE, False)
            self.assertTrue(all(kind != "zero" for _, kind, _ in off_metal))

    def test_every_ring_follows_the_metal_fractions(self):
        for topology in ("convex", "concave"):
            for angle in ANGLES + (78.0, 176.0):
                surface = GENERATOR.build_surface(
                    RADIUS, 16, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology, rule=self.RULE
                )
                points = np.asarray(surface.knot_points)
                self.assertEqual(len(points), 160)
                self.assertEqual(surface.contour_groups, [16] * 10)
                self.assertEqual(sum(surface.knot_zero), 14)
                # PEC knots on the two metal rings only, at the rule's zero slots.
                zero = np.asarray(surface.knot_zero).reshape(10, 16)
                for ring in range(10):
                    z = points[16 * ring, 2]
                    metal = abs(z) < 1e-12 or abs(z - THICKNESS) < 1e-12
                    self.assertEqual(list(np.flatnonzero(zero[ring])),
                                     GENERATOR.zero_slots(topology, self.RULE) if metal else [])
                # Identical perimeter fractions on every ring (outer and cap): regular columns.
                fractions = []
                for ring in range(10):
                    half_width = np.max(np.abs(points[16 * ring: 16 * ring + 16, :2]))
                    fractions.append(sorted(
                        GENERATOR.square_perimeter_fraction(half_width, p) for p in points[16 * ring: 16 * ring + 16]
                    ))
                for ring in range(1, 10):
                    np.testing.assert_allclose(fractions[ring], fractions[0], atol=1e-12)
                # Slaves: the box corners of every ring that are no knot, plus the two cap
                # centres at the mean of the crossing knots; a partition of unity.
                centres = [s for s in surface.slaves if abs(s[0][0]) < 1e-12 and abs(s[0][1]) < 1e-12]
                self.assertEqual(len(centres), 2)
                first, second = GENERATOR.arm_crossing_fractions(RADIUS, angle)
                for point, parent_a, parent_b, weight in centres:
                    self.assertEqual(weight, 0.5)
                    for parent, crossing in ((parent_a, first), (parent_b, second % 1.0)):
                        self.assertAlmostEqual(
                            GENERATOR.square_perimeter_fraction(RADIUS / 3.0, points[parent]), crossing, places=12)
                for _, parent_a, parent_b, weight in surface.slaves:
                    self.assertIsNotNone(parent_a)
                    self.assertTrue(0.0 <= weight <= 1.0)
                # No degenerate triangle; no free hat on the PEC contour.
                vertices = surface.vertex_points()
                tri = vertices[surface.triangles]
                areas = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
                self.assertGreater(areas.min(), 1e-8, (topology, angle))
                self.assertEqual(
                    GENERATOR.free_hat_pec_support(points, surface.contour_groups, np.flatnonzero(surface.knot_zero),
                                                   pec_mask(angle, topology), surface.slaves), [])
                # The rule's zero set is the PEC footprint of both coupons.
                np.testing.assert_array_equal(pec_mask(angle, topology)(points), np.asarray(surface.knot_zero))

    def test_no_connectivity_angle_or_cap_centers_option(self):
        with self.assertRaises(ValueError):
            GENERATOR.build_surface(RADIUS, 16, THICKNESS, OVERETCH, angle_degrees=120.0, topology="convex",
                                    rule=self.RULE, connectivity_angle_degrees=112.5)
        with self.assertRaises(ValueError):
            GENERATOR.build_surface(RADIUS, 8, THICKNESS, OVERETCH, angle_degrees=120.0, topology="convex",
                                    rule=self.RULE)

    def test_heldout_reference_levels_are_decoupled_from_the_basis(self):
        levels = GENERATOR.heldout_reference_levels(RADIUS, THICKNESS, OVERETCH)
        # Every OveretchDepth across both R / 3 ramps, the t + R / 3 level, the standard levels.
        for z in (-RADIUS, -RADIUS / 3.0, -OVERETCH, 0.0, THICKNESS, THICKNESS + RADIUS / 3.0, RADIUS / 3.0, RADIUS):
            self.assertTrue(any(abs(z - level) < 1e-12 for level in levels), z)
        k = 1
        while THICKNESS + k * OVERETCH < THICKNESS + RADIUS / 3.0 - 1e-12:
            self.assertTrue(any(abs(THICKNESS + k * OVERETCH - level) < 1e-12 for level in levels), k)
            k += 1
        self.assertEqual(k, 13)  # 0.15 .. 0.70
        k = 2
        while -k * OVERETCH > -RADIUS / 3.0 + 1e-12:
            self.assertTrue(any(abs(-k * OVERETCH - level) < 1e-12 for level in levels), k)
            k += 1
        self.assertEqual(k, 13)  # -0.10 .. -0.60
        self.assertEqual(len(levels), 31)
        self.assertEqual(GENERATOR.HELDOUT_REFERENCE_RING_SIZE, 64)
        fine = GENERATOR.build_surface(RADIUS, 64, THICKNESS, OVERETCH, angle_degrees=120.0, topology="convex",
                                       cap_centers=True, crossing_vertices=True, levels=levels)
        heights = sorted({round(p[2], 9) for p in fine.vertex_points() if np.max(np.abs(p[:2])) > RADIUS - 1e-9})
        self.assertEqual(heights, [round(z, 9) for z in levels])
        values = GENERATOR.heldout_potential(fine.vertex_points(), RADIUS, THICKNESS, 120.0, "convex")
        pec = pec_mask(120.0, "convex")(fine.vertex_points())
        np.testing.assert_array_equal(values[pec], 0.0)
        self.assertTrue(np.all(values[~pec] != 0.0))


if __name__ == "__main__":
    unittest.main()
