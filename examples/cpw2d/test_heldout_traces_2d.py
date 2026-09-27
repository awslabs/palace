#!/usr/bin/env python3

"""Held-out traces of the 2D response generators: the polynomial blends to the potential of
the conductor at every cut (a live terminal's potential, zero at the ground), so the trace
is a continuous field the hat basis can represent; ground-only coupons are unchanged."""

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np

CPW2D = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, CPW2D / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CLUSTER = load("generate_edge_cluster_response")
PAIR = load("generate_edge_pair_response")

RADIUS = 1.9
METAL_THICKNESS = 0.1


def polynomial(points):
    x = points[:, 0] / RADIUS
    y = points[:, 1] / RADIUS
    return 0.35 + 0.20 * x - 0.15 * y + 0.08 * x * y + 0.06 * y * y


def smoothstep(distance):
    coordinate = np.clip(distance / (RADIUS / 3.0), 0.0, 1.0)
    return coordinate * coordinate * (3.0 - 2.0 * coordinate)


def cut_distance(points, boundary_x):
    vertical = np.maximum.reduce(
        (-points[:, 1], points[:, 1] - METAL_THICKNESS, np.zeros(len(points)))
    )
    return np.hypot(points[:, 0] - boundary_x, vertical)


def hat_interpolant(output, coefficients, trace_paths):
    values = None
    for coefficient, path in zip(coefficients, trace_paths):
        trace = np.loadtxt(path, delimiter=",", skiprows=1)
        values = coefficient * trace[:, 3] if values is None else values + coefficient * trace[:, 3]
    return values


class ClusterHeldoutTraceTest(unittest.TestCase):
    def build(self, offsets, directions, conductors, basis_size=96):
        root = Path(self.directory.name) / "-".join(str(d) for d in directions).replace("-", "m")
        root.mkdir()
        traces, conductor_traces, _, _, points, cuts = CLUSTER.write_bases(
            root,
            np.asarray(offsets, dtype=float),
            np.asarray(directions),
            np.asarray(conductors),
            RADIUS,
            METAL_THICKNESS,
            basis_size,
            13 * basis_size,
        )
        CLUSTER.write_heldout(
            root, traces, conductor_traces, points, cuts, RADIUS, METAL_THICKNESS
        )
        heldout = np.loadtxt(root / "heldout_trace.csv", delimiter=",", skiprows=1)
        coefficients = np.atleast_1d(
            np.loadtxt(root / "heldout_coefficients.csv", delimiter=",", skiprows=1)
        )
        return root, traces, conductor_traces, cuts, points, heldout, coefficients

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.directory.cleanup()

    @staticmethod
    def former_cutoff_to_zero_error(root, traces, conductor_traces, cuts, points):
        # The recorded construction: cutoff polynomial (to zero at every cut) plus the
        # conductor traces at 1 V; its hat interpolation error on the same knots.
        distance = np.min(
            np.column_stack([cut_distance(points, cut[3]) for cut in cuts]), axis=1
        )
        values = smoothstep(distance) * polynomial(points)
        basis_points = np.loadtxt(root / "basis_points.csv", delimiter=",", skiprows=1)
        knot_distance = np.min(
            np.column_stack([cut_distance(basis_points, cut[3]) for cut in cuts]), axis=1
        )
        coefficients = list(smoothstep(knot_distance) * polynomial(basis_points))
        for path in conductor_traces:
            coefficients.append(1.0)
            values = values + np.loadtxt(path, delimiter=",", skiprows=1)[:, 3]
        interpolant = hat_interpolant(root, coefficients, list(traces) + list(conductor_traces))
        return np.abs(values - interpolant).max()

    def test_live_cut_trace_is_continuous_and_carries_the_conductor_potential(self):
        # k4 2 / 2 / 2 um, three conductors: the right metal (conductor 3) is a terminal at
        # 1 V meeting the box at the right side (the recorded -12 % MS "reconstruction
        # error" coupon).
        root, traces, conductor_traces, cuts, points, heldout, coefficients = self.build(
            [0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1], [1, 2, 2, 3]
        )
        self.assertEqual(len(cuts), 2)
        self.assertEqual([cut[2] for cut in cuts], [1, 3])
        self.assertEqual(len(coefficients), len(traces) + len(conductor_traces))
        np.testing.assert_array_equal(coefficients[len(traces) :], [1.0, 1.0])
        values = heldout[:, 3]
        # Continuous along the contour: the largest step between consecutive samples is
        # the polynomial's own variation over one sample spacing (no 1 V step at a cut).
        steps = np.abs(np.diff(values))
        self.assertLess(steps.max(), 0.02)
        # On the live cut the trace is the conductor's potential; next to it the blend.
        right_x = cuts[1][3]
        on_right_cut = cut_distance(points, right_x) <= 1.0e-12
        self.assertGreater(on_right_cut.sum(), 0)
        np.testing.assert_allclose(values[on_right_cut], 1.0, rtol=0.0, atol=1.0e-12)
        distance = cut_distance(points, right_x)
        near = (distance > 1.0e-12) & (distance < RADIUS / 3.0)
        expected = smoothstep(distance[near]) * polynomial(points[near]) + (
            1.0 - smoothstep(distance[near])
        )
        np.testing.assert_allclose(values[near], expected, rtol=0.0, atol=1.0e-12)
        # On the ground cut and next to it the trace is the plain cutoff polynomial.
        left_x = cuts[0][3]
        distance = cut_distance(points, left_x)
        near_left = distance < RADIUS / 3.0
        np.testing.assert_allclose(
            values[near_left],
            smoothstep(distance[near_left]) * polynomial(points[near_left]),
            rtol=0.0,
            atol=1.0e-12,
        )
        # The hat interpolant of the coefficients deviates from the trace by the
        # interpolation error of a 0.2-0.4 V blend, several times below the error of the
        # former cutoff-to-zero construction (a 1 V step against the metal) on the same hats.
        interpolant = hat_interpolant(root, coefficients, list(traces) + list(conductor_traces))
        error = np.abs(values - interpolant).max()
        former = self.former_cutoff_to_zero_error(root, traces, conductor_traces, cuts, points)
        self.assertLess(error, 0.25 * former)
        self.assertLess(error, 0.02)

    def test_ground_only_stack_trace_is_unchanged(self):
        # The threshold-type stack: the live strip lies inside the box, the only cut is the
        # ground's; the trace is the recorded cutoff polynomial with a zero conductor trace.
        root, traces, conductor_traces, cuts, points, heldout, coefficients = self.build(
            [0.0, 3.7715, 7.543], [1, -1, 1], [1, 2, 2]
        )
        self.assertEqual([cut[2] for cut in cuts], [1])
        distance = cut_distance(points, cuts[0][3])
        np.testing.assert_array_equal(
            heldout[:, 3], smoothstep(distance) * polynomial(points)
        )
        basis_points = np.loadtxt(root / "basis_points.csv", delimiter=",", skiprows=1)
        np.testing.assert_array_equal(
            coefficients[: len(traces)],
            smoothstep(cut_distance(basis_points, cuts[0][3])) * polynomial(basis_points),
        )
        np.testing.assert_array_equal(coefficients[len(traces) :], [1.0])
        conductor = np.loadtxt(conductor_traces[0], delimiter=",", skiprows=1)
        np.testing.assert_array_equal(conductor[:, 3], 0.0)

    def test_probe_traces_stay_cut_off_to_ground(self):
        # The probe solves ground every conductor: the probes keep the cutoff to zero at
        # the live cut too.
        root, traces, conductor_traces, cuts, points, heldout, coefficients = self.build(
            [0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1], [1, 2, 2, 3]
        )
        probe = np.loadtxt(root / "probe_01.csv", delimiter=",", skiprows=1)
        distance = np.minimum(
            cut_distance(points, cuts[0][3]), cut_distance(points, cuts[1][3])
        )
        np.testing.assert_array_equal(probe[:, 3], smoothstep(distance) * np.ones(len(points)))


class PairHeldoutTraceTest(unittest.TestCase):
    def build(self, separation, different_conductors, strip, basis_size=96):
        root = Path(self.directory.name) / f"{different_conductors}-{strip}"
        root.mkdir()
        traces, conductor_trace, _ = PAIR.write_bases(
            root,
            separation,
            RADIUS,
            METAL_THICKNESS,
            basis_size,
            13 * basis_size,
            different_conductors,
            strip,
        )
        PAIR.write_heldout(
            root, traces, conductor_trace, separation, RADIUS, METAL_THICKNESS, strip
        )
        heldout = np.loadtxt(root / "heldout_trace.csv", delimiter=",", skiprows=1)
        coefficients = np.atleast_1d(
            np.loadtxt(root / "heldout_coefficients.csv", delimiter=",", skiprows=1)
        )
        return root, traces, conductor_trace, heldout, coefficients

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.directory.cleanup()

    def test_same_conductor_gap_and_strip_traces_are_unchanged(self):
        separation = 2.0
        half_width = 0.5 * separation + RADIUS
        root, traces, conductor_trace, heldout, coefficients = self.build(
            separation, False, False
        )
        self.assertIsNone(conductor_trace)
        points = heldout[:, :3]
        distance = np.minimum(
            cut_distance(points, -half_width), cut_distance(points, half_width)
        )
        np.testing.assert_array_equal(heldout[:, 3], smoothstep(distance) * polynomial(points))
        self.assertEqual(len(coefficients), len(traces))
        root, traces, conductor_trace, heldout, coefficients = self.build(
            separation, False, True
        )
        np.testing.assert_array_equal(heldout[:, 3], polynomial(heldout[:, :3]))

    def test_different_conductor_gap_blends_to_the_terminal_potential(self):
        separation = 2.0
        half_width = 0.5 * separation + RADIUS
        root, traces, conductor_trace, heldout, coefficients = self.build(
            separation, True, False
        )
        self.assertIsNotNone(conductor_trace)
        self.assertEqual(coefficients[-1], 0.17)
        points, values = heldout[:, :3], heldout[:, 3]
        right = cut_distance(points, half_width)
        on_cut = right <= 1.0e-12
        self.assertGreater(on_cut.sum(), 0)
        np.testing.assert_allclose(values[on_cut], 0.17, rtol=0.0, atol=1.0e-12)
        near = (right > 1.0e-12) & (right < RADIUS / 3.0)
        np.testing.assert_allclose(
            values[near],
            smoothstep(right[near]) * polynomial(points[near])
            + (1.0 - smoothstep(right[near])) * 0.17,
            rtol=0.0,
            atol=1.0e-12,
        )
        left = cut_distance(points, -half_width)
        near_left = left < RADIUS / 3.0
        np.testing.assert_allclose(
            values[near_left],
            smoothstep(left[near_left]) * polynomial(points[near_left]),
            rtol=0.0,
            atol=1.0e-12,
        )
        self.assertLess(np.abs(np.diff(values)).max(), 0.02)
        interpolant = hat_interpolant(root, coefficients, list(traces) + [conductor_trace])
        error = np.abs(values - interpolant).max()
        # The recorded construction (cutoff to zero, conductor trace at 0.17 V) on the same
        # knots: a 0.17 V step against the metal the hats cannot follow.
        distance = np.minimum(left, right)
        former_values = smoothstep(distance) * polynomial(points) + 0.17 * np.loadtxt(
            conductor_trace, delimiter=",", skiprows=1
        )[:, 3]
        basis_points = np.loadtxt(root / "basis_points.csv", delimiter=",", skiprows=1)
        knot_distance = np.minimum(
            cut_distance(basis_points, -half_width), cut_distance(basis_points, half_width)
        )
        former_coefficients = list(smoothstep(knot_distance) * polynomial(basis_points)) + [0.17]
        former = np.abs(
            former_values
            - hat_interpolant(root, former_coefficients, list(traces) + [conductor_trace])
        ).max()
        # A 0.17 V terminal: the blend amplitude (polynomial - 0.17) is 0.7 of the former
        # step, so the interpolation error shrinks by that ratio, not more.
        self.assertLess(error, former)
        self.assertLess(error, 0.8 * former)


if __name__ == "__main__":
    unittest.main()
