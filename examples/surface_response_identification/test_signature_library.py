# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Canonical chording of SpatialEdgeCluster arc portions (option A)."""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from surface_response_identification import signature_library as L  # noqa: E402


def arc_portion(center, radius, theta0, theta1, gap_radial, conductor=1):
    p0 = (center[0] + radius * math.cos(theta0), center[1] + radius * math.sin(theta0))
    p1 = (center[0] + radius * math.cos(theta1), center[1] + radius * math.sin(theta1))
    tm = 0.5 * (theta0 + theta1)
    m = (center[0] + radius * math.cos(tm), center[1] + radius * math.sin(tm))
    a, b = sorted([p0, p1])
    return {"P": [a[0], a[1], b[0], b[1]], "Arc": [center[0], center[1], m[0], m[1]], "GapRadial": gap_radial, "Conductor": conductor, "Interfaces": ["MetalAir"], "Law": "{}"}


class ClusterChordingTest(unittest.TestCase):
    R = 2.0

    def test_straight_portion_is_one_edge(self):
        signature = {"Portions": [{"P": [0.0, 0.0, 1.5, 0.0], "Gap": [0.0, 1.0], "Conductor": 1}]}
        edges = L.cluster_plan_view_edges(signature, self.R)
        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0]["P0"], (0.0, 0.0))
        self.assertEqual(edges[0]["P1"], (3.0, 0.0))
        self.assertEqual(edges[0]["Gap"], (0.0, 1.0))

    def test_quarter_arc_chord_count_and_geometry(self):
        # A 90 deg arc of radius 0.5 R: 90 / 5 = 18 chords by the angular step; the length
        # rule (pi / 4 R over 0.25 R = 4 chords) is weaker.
        signature = {"Portions": [arc_portion((0.0, 0.0), 0.5, 0.0, 0.5 * math.pi, +1)]}
        edges = L.cluster_plan_view_edges(signature, self.R)
        self.assertEqual(len(edges), 18)
        for k, edge in enumerate(edges):
            for point in (edge["P0"], edge["P1"]):
                self.assertAlmostEqual(math.hypot(*point), 0.5 * self.R, places=12)
            self.assertEqual(edge["Chord"], k)
            # Outward gap (metal inside the circle): the gap direction is the outward radial.
            mid = (0.5 * (edge["P0"][0] + edge["P1"][0]), 0.5 * (edge["P0"][1] + edge["P1"][1]))
            self.assertGreater(edge["Gap"][0] * mid[0] + edge["Gap"][1] * mid[1], 0.0)
        # The chording runs from the lexicographically smaller end to the larger one.
        ends = sorted([edges[0]["P0"], edges[-1]["P1"]])
        self.assertAlmostEqual(ends[0][1], 0.5 * self.R, places=12)
        self.assertAlmostEqual(ends[1][0], 0.5 * self.R, places=12)

    def test_length_rule_refines_a_large_radius(self):
        # A 20 deg arc of radius 10 R: 4 chords by the angle, 10 R x 0.349 / 0.25 R = 14 by length.
        signature = {"Portions": [arc_portion((3.0, -2.0), 10.0, 0.3, 0.3 + math.radians(20.0), -1)]}
        edges = L.cluster_plan_view_edges(signature, self.R)
        self.assertEqual(len(edges), 14)
        # Inward gap (metal outside the circle).
        edge = edges[5]
        mid = (0.5 * (edge["P0"][0] + edge["P1"][0]) - 3.0 * self.R, 0.5 * (edge["P0"][1] + edge["P1"][1]) + 2.0 * self.R)
        self.assertLess(edge["Gap"][0] * mid[0] + edge["Gap"][1] * mid[1], 0.0)

    def test_side_of_the_arc_from_the_midpoint(self):
        # The same two end points, the midpoint on the other side: the long way round.
        short = arc_portion((0.0, 0.0), 1.0, -0.25 * math.pi, 0.25 * math.pi, +1)
        long_way = dict(short, Arc=[0.0, 0.0, -1.0 * self.R / self.R, 0.0])
        n_short = len(L.cluster_plan_view_edges({"Portions": [short]}, self.R))
        n_long = len(L.cluster_plan_view_edges({"Portions": [long_way]}, self.R))
        self.assertEqual(n_short, 18)
        self.assertEqual(n_long, 54)

    def test_closed_circle(self):
        portion = arc_portion((0.0, 0.0), 0.3, 0.0, 2.0 * math.pi, +1)
        edges = L.cluster_plan_view_edges({"Portions": [portion]}, self.R)
        self.assertEqual(len(edges), 72)
        self.assertAlmostEqual(edges[0]["P0"][0], edges[-1]["P1"][0], places=12)
        self.assertAlmostEqual(edges[0]["P0"][1], edges[-1]["P1"][1], places=12)

    def test_chording_is_mesh_independent_by_construction(self):
        # Two serialisations of one arc (the builder's input) give identical edges.
        a = arc_portion((1.0, 1.0), 0.75, 0.1, 1.3, +1)
        b = dict(a)
        self.assertEqual(L.cluster_plan_view_edges({"Portions": [a]}, self.R), L.cluster_plan_view_edges({"Portions": [b]}, self.R))




class SpatialSupportContractTest(unittest.TestCase):
    """Contract v3 (decision 282): the claims-derived box and the Context entries."""

    # The 4-edge lead-end cluster of the C++ unit case SurfaceResponseIdentificationSpatialSupportContract
    # (pad edge claimed over 6.873 R, a 2-um lead ending 1 um above it; R = 2): the C++ record's ClaimsBox.
    SIGNATURE = {"Type": "SpatialEdgeCluster", "EdgeCount": 4, "Portions": [
        {"Conductor": 1, "Gap": [-1.0, 0.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.218559, -0.5, -0.218559, 0.5]},
        {"Conductor": 1, "Gap": [0.0, -1.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.218559, -0.5, 2.281441, -0.5]},
        {"Conductor": 1, "Gap": [0.0, 1.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.218559, 0.5, 2.281441, 0.5]},
        {"Conductor": 2, "Gap": [1.0, 0.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.718559, -3.436492, -0.718559, 3.436492]}],
        "Vertices": [{"P": [-0.218559, -0.5], "TurnDegrees": 90.0, "Type": "ConvexCorner"}, {"P": [-0.218559, 0.5], "TurnDegrees": 90.0, "Type": "ConvexCorner"}]}
    CPP_BOX = [-3.218559, -6.436492, 5.281441, 6.436492]

    def test_box_matches_the_cpp_record(self):
        box = L.cluster_support_box(self.SIGNATURE)
        for value, expected in zip(box, self.CPP_BOX):
            self.assertAlmostEqual(value, expected, places=9)

    def test_box_rule_on_one_free_portion(self):
        # A lone 2 R portion: both ends free (|half| = 1 >= R) -> continued by 2R, widened by R, padded by R.
        signature = {"Portions": [{"P": [-1.0, 0.0, 1.0, 0.0], "Gap": [0.0, 1.0], "Conductor": 1}]}
        self.assertEqual(L.cluster_support_box(signature), [-4.0, -2.0, 4.0, 2.0])
        # A connected short end (a vertex there) keeps its length: no continuation on that side.
        signature["Vertices"] = [{"P": [-0.5, 0.0], "Type": "ConvexCorner", "TurnDegrees": 90.0}]
        signature["Portions"][0]["P"] = [-0.5, 0.0, 0.5, 0.0]
        self.assertEqual(L.cluster_support_box(signature), [-1.5, -2.0, 4.0, 2.0])

    def test_context_entries_follow_the_claims(self):
        signature = dict(self.SIGNATURE)
        signature["Box"] = list(self.CPP_BOX)
        signature["Context"] = [{"Chain": True, "Conductor": 1, "Gap": [0.0, -1.0], "Interfaces": ["MS"], "Law": "{}", "P": [2.281441, -0.5, 5.281441, -0.5]},
                                {"Chain": False, "Conductor": 3, "Gap": [0.0, 1.0], "Interfaces": ["MS"], "Law": "{}", "P": [-3.0, 5.0, -2.0, 5.0]}]
        claims = L.cluster_plan_view_edges(signature, 2.0)
        self.assertEqual(len(claims), 4)
        self.assertTrue(all("Context" not in e for e in claims))
        edges = L.cluster_plan_view_edges(signature, 2.0, include_context=True)
        self.assertEqual(len(edges), 6)
        self.assertEqual([e.get("Context", False) for e in edges], [False] * 4 + [True] * 2)
        self.assertEqual([e["Chain"] for e in edges[4:]], [True, False])
        self.assertEqual([e["Portion"] for e in edges[4:]], [0, 1])
        self.assertEqual(edges[4]["P0"], (2.0 * 2.281441, -1.0))
        self.assertEqual(L.conductor_count({"Type": "SpatialEdgeCluster", "Signature": signature}), 3)
        self.assertEqual(L.conductor_count({"Type": "SpatialEdgeCluster", "Signature": self.SIGNATURE}), 2)


if __name__ == "__main__":
    unittest.main()
