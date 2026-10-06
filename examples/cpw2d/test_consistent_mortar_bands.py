#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The constrained metal band of every straight coupon generator against the runtime's
consistent-mortar rule (decision 404 D1, SurfaceResponseOperator DeriveConsistentMortarBands).

Each generator's write_bases inserts the knots where the fabricated metal meets the matching
contour, (x_face, 0) and (x_face, MetalThickness), constrains the band between them and
publishes the FREE knots only. The runtime inserts the band vertices back from the topology
(IsolatedEdge / CurvedEdge: the face x = min; SameConductorGap: both faces; a ParallelEdge-
Cluster: the face of an outer edge whose gap points inward; strips: none) and
Fabrication.MetalThickness, between the two consecutive published knots on either side of the
band. This test states the generators' side of that contract: the constrained knots are
exactly the rule's band ends on the rule's faces, no published knot lies on a band, the band
sits between consecutive published knots of the closed contour, strips constrain nothing,
and the open-contour (different-conductor) models anchor their band knots on conductor
references (outside the rule, recorded as such by the runtime)."""

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


EDGE = load("generate_edge_response")
PAIR = load("generate_edge_pair_response")
CLUSTER = load("generate_edge_cluster_response")

RADIUS = 1.9
METAL_THICKNESS = 0.1
BASIS_SIZE = 96
TOLERANCE = 1.0e-9 * RADIUS


def runtime_rule_faces(topology, edges=None):
    """The faces the runtime rule derives (True = x min, False = x max)."""
    if topology in ("IsolatedEdge", "CurvedEdge"):
        return [True]
    if topology in ("SameConductorGap", "CurvedSameConductorGap"):
        return [True, False]
    if topology in ("SameConductorStrip", "CurvedSameConductorStrip"):
        return []
    if topology == "ParallelEdgeCluster":
        faces = []
        if edges[0][1] == 1:
            faces.append(True)
        if edges[-1][1] == -1:
            faces.append(False)
        return faces
    raise ValueError(topology)


def runtime_rule_insertion(points, faces, t=METAL_THICKNESS):
    """The runtime rule on the PUBLISHED knots (closed contour, listed order): per face the
    two consecutive knots on either side of the band [0, t]; raises like the runtime fails
    closed when a published knot lies on the band or the band is not between consecutive
    knots."""
    n = len(points)
    x_min, x_max = points[:, 0].min(), points[:, 0].max()
    insertions = []
    for left in faces:
        x_face = x_min if left else x_max
        face = [i for i in range(n) if abs(points[i, 0] - x_face) <= TOLERANCE]
        face.sort(key=lambda i: points[i, 1])
        on_band = [i for i in face if -TOLERANCE <= points[i, 1] <= t + TOLERANCE]
        if on_band:
            raise AssertionError(f"published knot(s) {on_band} on the band of x = {x_face}")
        below = [i for i in face if points[i, 1] < -TOLERANCE]
        above = [i for i in face if points[i, 1] > t + TOLERANCE]
        if not below or not above:
            raise AssertionError(f"no knots on both sides of the band of x = {x_face}")
        a, b = below[-1], above[0]
        if (a + 1) % n != b and (b + 1) % n != a:
            raise AssertionError(f"knots {a} and {b} are not consecutive on x = {x_face}")
        insertions.append((x_face, a, b))
    return insertions


class BandConstraintTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)

    def tearDown(self):
        self.directory.cleanup()

    def published(self, output):
        return np.atleast_2d(np.loadtxt(output / "basis_points.csv", delimiter=",", skiprows=1))

    @staticmethod
    def constrained_knots(all_knots, published):
        """The generator's constrained knots = the contour knots it did not publish (a
        uniform knot a rounding error away from a junction counts once)."""
        constrained = []
        for knot in all_knots:
            if np.any(np.all(np.abs(published - knot) <= TOLERANCE, axis=1)):
                continue
            if any(np.all(np.abs(other - knot) <= TOLERANCE) for other in constrained):
                continue
            constrained.append(knot)
        return np.asarray(constrained).reshape(-1, 3)

    def check_bands(self, published, constrained, faces):
        insertions = runtime_rule_insertion(published, faces)
        self.assertEqual(len(insertions), len(faces))
        expected = []
        for x_face, _, _ in insertions:
            expected += [(x_face, 0.0), (x_face, METAL_THICKNESS)]
        self.assertEqual(len(constrained), len(expected), (constrained, expected))
        for x, y in expected:
            hit = np.any(
                (np.abs(constrained[:, 0] - x) <= TOLERANCE)
                & (np.abs(constrained[:, 1] - y) <= TOLERANCE)
            )
            self.assertTrue(hit, f"constrained knot ({x}, {y}) missing: {constrained}")
        return insertions

    def test_isolated_edge(self):
        output = self.root / "isolated"
        output.mkdir()
        EDGE.write_bases(output, RADIUS, METAL_THICKNESS, BASIS_SIZE, 13 * BASIS_SIZE)
        published = self.published(output)
        # Every contour knot the generator considered: the uniform knots plus the junctions.
        perimeter = 8.0 * RADIUS
        spacing = perimeter / BASIS_SIZE
        lower = 7.0 * RADIUS
        distances = np.unique(
            np.concatenate(
                (
                    (lower % spacing + spacing * np.arange(BASIS_SIZE)) % perimeter,
                    [lower - METAL_THICKNESS, lower],
                )
            )
        )
        all_knots = np.asarray([EDGE.contour_point(d, RADIUS) for d in distances])
        constrained = self.constrained_knots(all_knots, published)
        insertions = self.check_bands(published, constrained, runtime_rule_faces("IsolatedEdge"))
        # The record's reading: 95 free knots of 97, the band between the knots at
        # (-R, +spacing) and (-R, -spacing) (the knot at (-R, 0) is itself constrained).
        self.assertEqual(len(published), 95)
        x_face, a, b = insertions[0]
        self.assertAlmostEqual(x_face, -RADIUS, places=12)
        self.assertAlmostEqual(abs(published[a, 1]), spacing, places=9)
        self.assertAlmostEqual(abs(published[b, 1]), spacing, places=9)

    def pair(self, name, different_conductors, strip, separation=2.0):
        output = self.root / name
        output.mkdir()
        PAIR.PLACEMENT = PAIR.PairPlacement(0.0, "convex", separation, strip)
        paths, conductor_trace, open_paths = PAIR.write_bases(
            output, separation, RADIUS, METAL_THICKNESS, BASIS_SIZE, 13 * BASIS_SIZE,
            different_conductors, strip,
        )
        half_width = 0.5 * separation + RADIUS
        perimeter = 4.0 * (half_width + RADIUS)
        spacing = perimeter / BASIS_SIZE
        if strip:
            distances = spacing * np.arange(BASIS_SIZE)
        else:
            right_lower = 2.0 * half_width + RADIUS
            left_lower = right_lower + 0.5 * perimeter
            distances = np.unique(
                np.concatenate(
                    (
                        (right_lower % spacing + spacing * np.arange(BASIS_SIZE)) % perimeter,
                        [right_lower, right_lower + METAL_THICKNESS,
                         left_lower - METAL_THICKNESS, left_lower],
                    )
                )
            )
        all_knots = np.asarray([PAIR.contour_point(d, half_width, RADIUS) for d in distances])
        published = self.published(output)
        return published, self.constrained_knots(all_knots, published), open_paths

    def test_same_conductor_gap_constrains_both_faces(self):
        published, constrained, open_paths = self.pair("gap", False, False)
        self.assertEqual(open_paths, [])
        insertions = self.check_bands(published, constrained, runtime_rule_faces("SameConductorGap"))
        self.assertEqual(sorted(x for x, _, _ in insertions), [-(1.0 + RADIUS), 1.0 + RADIUS])

    def test_same_conductor_strip_constrains_nothing(self):
        published, constrained, open_paths = self.pair("strip", False, True)
        self.assertEqual(open_paths, [])
        self.assertEqual(len(constrained), 0)
        self.assertEqual(len(published), BASIS_SIZE)
        self.assertEqual(runtime_rule_faces("SameConductorStrip"), [])

    def test_different_conductor_gap_anchors_the_bands_on_open_paths(self):
        # The band knots are conductor anchors (StartConductor / EndConductor of the two open
        # paths), not zero vertices: outside the consistent-mortar rule, recorded by the
        # runtime as such. The generator still constrains exactly the two bands.
        published, constrained, open_paths = self.pair("different", True, False)
        self.assertEqual(len(open_paths), 2)
        self.assertEqual({p["StartConductor"] for p in open_paths}, {1})
        self.assertEqual({p["EndConductor"] for p in open_paths}, {2})
        self.assertEqual(len(constrained), 4)
        for x in (-(1.0 + RADIUS), 1.0 + RADIUS):
            for y in (0.0, METAL_THICKNESS):
                self.assertTrue(np.any((np.abs(constrained[:, 0] - x) <= TOLERANCE)
                                       & (np.abs(constrained[:, 1] - y) <= TOLERANCE)))

    def cluster(self, name, offsets, directions, conductors):
        output = self.root / name
        output.mkdir()
        traces, conductor_traces, references, open_paths, points, cuts = CLUSTER.write_bases(
            output, np.asarray(offsets, dtype=float), np.asarray(directions),
            np.asarray(conductors), RADIUS, METAL_THICKNESS, BASIS_SIZE, 13 * BASIS_SIZE,
        )
        xmin, xmax = offsets[0] - RADIUS, offsets[-1] + RADIUS
        width = xmax - xmin
        perimeter = 2.0 * (width + 2.0 * RADIUS)
        spacing = perimeter / BASIS_SIZE
        endpoints = [value for cut in cuts for value in cut[:2]]
        offset = cuts[0][1] % spacing if cuts else 0.0
        distances = np.unique(
            np.concatenate(((offset + spacing * np.arange(BASIS_SIZE)) % perimeter, endpoints))
        )
        all_knots = np.asarray([CLUSTER.contour_point(d, xmin, xmax, RADIUS) for d in distances])
        published = self.published(output)
        edges = list(zip(offsets, directions, conductors))
        return published, self.constrained_knots(all_knots, published), open_paths, edges

    def test_single_conductor_stack_constrains_the_inward_gap_faces(self):
        # S1p's 4-edge stack (gap, metal, gap: the outer edges' gaps point inward): both faces.
        published, constrained, open_paths, edges = self.cluster(
            "stack4", [0.0, 2.0, 4.0, 6.0], [1, -1, 1, -1], [1, 1, 1, 1])
        self.assertEqual(open_paths, [])
        faces = runtime_rule_faces("ParallelEdgeCluster", edges)
        self.assertEqual(faces, [True, False])
        self.check_bands(published, constrained, faces)
        # The metal strip between the first two edges: a 3-edge stack whose first edge's gap
        # points outward meets the metal on the right face only.
        published, constrained, open_paths, edges = self.cluster(
            "stack3", [0.0, 2.0, 4.0], [-1, 1, -1], [1, 1, 1])
        self.assertEqual(open_paths, [])
        faces = runtime_rule_faces("ParallelEdgeCluster", edges)
        self.assertEqual(faces, [False])
        self.check_bands(published, constrained, faces)

    def test_two_conductor_stack_anchors_its_bands_on_open_paths(self):
        published, constrained, open_paths, edges = self.cluster(
            "stack-two-conductors", [0.0, 2.0, 4.0], [-1, 1, -1], [1, 1, 2])
        self.assertEqual(len(open_paths), 2)
        self.assertEqual(len(constrained), 2)  # the right face band of conductor 2

    def test_runtime_rule_refuses_a_knot_on_the_band(self):
        output = self.root / "isolated-refused"
        output.mkdir()
        EDGE.write_bases(output, RADIUS, METAL_THICKNESS, BASIS_SIZE, 13 * BASIS_SIZE)
        published = self.published(output)
        # A library that listed the band knot (-R, 0) as free disagrees with its constraint.
        tampered = np.vstack((published, [[-RADIUS, 0.0, 0.0]]))
        with self.assertRaises(AssertionError):
            runtime_rule_insertion(tampered, runtime_rule_faces("IsolatedEdge"))


if __name__ == "__main__":
    unittest.main()
