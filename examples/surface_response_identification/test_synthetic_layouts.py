#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the synthetic stress layouts: oracle on shapes with known answers, the
specification writer, and the polygon builders."""

import math
import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import synthetic_layouts as S  # noqa: E402


class OracleTest(unittest.TestCase):
    def test_rectangle_pair_at_separation_2_is_a_same_conductor_gap_with_absorbed_corners(self):
        lay = S.layout("gap", [S.sheet(S.GROUND, S.rectangle(-15.0, -6.0, -1.0, 6.0)), S.sheet(S.GROUND, S.rectangle(1.0, -6.0, 15.0, 6.0))])
        orc = S.oracle(lay)
        self.assertAlmostEqual(orc["PhysicalPerimeterLength"], 2 * (2 * 14.0 + 2 * 12.0))
        self.assertEqual(orc["ClassifierCornerCount"], 8)
        # The four corners at the gap are within 2R of the facing edge: absorbed into a cluster.
        self.assertEqual(orc["StandaloneCornerCount"], 4)
        self.assertTrue(all(c["Convex"] and c["InteriorAngleDegrees"] == 90.0 for c in orc["Corners"]))
        pairs = orc["ParallelPairs"]
        self.assertEqual([(p["Expected"], p["Separation"], p["Within2R"], p["AtExactlyR"]) for p in pairs], [("SameConductorGap", 2.0, True, True)])
        self.assertEqual(pairs[0]["OverlapLength"], 12.0)

    def test_pair_at_exactly_2r_is_at_the_threshold_and_not_within(self):
        lay = S.layout("gap4", [S.sheet(S.GROUND, S.rectangle(-16.0, -6.0, -2.0, 6.0)), S.sheet(S.SECOND_CONDUCTOR, S.rectangle(2.0, -6.0, 16.0, 6.0))])
        orc = S.oracle(lay)
        pair = orc["ParallelPairs"][0]
        self.assertEqual(pair["Expected"], "DifferentConductorGap")
        self.assertTrue(pair["AtExactly2R"])
        self.assertFalse(pair["Within2R"])
        self.assertEqual(orc["StandaloneCornerCount"], 8)

    def test_strip_edges_have_metal_between_them(self):
        lay = S.layout("strip", [S.sheet(S.GROUND, S.rectangle(-12.0, -1.0, 12.0, 1.0))])
        orc = S.oracle(lay)
        self.assertEqual([p["Expected"] for p in orc["ParallelPairs"]], ["SameConductorStrip"])
        self.assertEqual(orc["ParallelPairs"][0]["Separation"], 2.0)

    def test_bent_bar_corner_angles_and_convexity(self):
        for angle in (30.0, 90.0, 135.0):
            lay = S.layout("bar", [S.sheet(S.GROUND, S.bent_bar(6.0, 16.0, angle))])
            orc = S.oracle(lay)
            # The arm ends are 16 um from the bend; the bend corners lie within 3 / sin(15 deg) = 11.6 um of the origin.
            bend = [c for c in orc["Corners"] if abs(c["TurnDegrees"] - (180.0 - angle)) < 1.0e-9 and np.linalg.norm(c["Point"]) < 14.0]
            self.assertEqual(len(bend), 2, angle)
            convex = [c for c in bend if c["Convex"]]
            concave = [c for c in bend if not c["Convex"]]
            self.assertEqual((len(convex), len(concave)), (1, 1))
            self.assertAlmostEqual(convex[0]["InteriorAngleDegrees"], angle)
            self.assertAlmostEqual(concave[0]["InteriorAngleDegrees"], 360.0 - angle)
            ends = [c for c in orc["Corners"] if abs(c["TurnDegrees"] - 90.0) < 1.0e-9 and np.linalg.norm(c["Point"]) > 14.0]
            self.assertEqual(len(ends), 4)
            self.assertEqual(len(orc["Corners"]), 6)

    def test_turn_of_exactly_30_degrees_is_straight_for_the_classifier(self):
        lay = S.layout("bar150", [S.sheet(S.GROUND, S.bent_bar(6.0, 16.0, 150.0))])
        orc = S.oracle(lay)
        self.assertEqual(orc["ClassifierCornerCount"], 4)
        self.assertEqual(orc["SubThresholdTurns"], 2)

    def test_truncation_edge_is_not_physical_perimeter(self):
        lay = S.layout("edge", [S.sheet(S.GROUND, S.rectangle(-30.0, -3.0, 6.0, 3.0))])
        orc = S.oracle(lay)
        self.assertAlmostEqual(orc["TruncationLength"], 6.0)
        self.assertAlmostEqual(orc["PhysicalPerimeterLength"], 36.0 * 2 + 6.0)
        self.assertEqual(orc["ClassifierCornerCount"], 2)

    def test_hole_edges_face_a_gap_and_the_corners_are_concave(self):
        lay = S.layout("hole", [S.sheet(S.GROUND, S.rectangle(-30.0, -30.0, 30.0, 30.0), holes=[S.rectangle(-5.0, -5.0, 5.0, 5.0)])])
        orc = S.oracle(lay)
        self.assertAlmostEqual(orc["PhysicalPerimeterLength"], 40.0)
        self.assertAlmostEqual(orc["TruncationLength"], 240.0)
        self.assertEqual(len(orc["Corners"]), 4)
        self.assertTrue(all(not c["Convex"] and c["InteriorAngleDegrees"] == 270.0 for c in orc["Corners"]))
        self.assertEqual(orc["ParallelPairs"], [])
        small = S.oracle(S.layout("hole1", [S.sheet(S.GROUND, S.rectangle(-30.0, -30.0, 30.0, 30.0), holes=[S.rectangle(-0.5, -0.5, 0.5, 0.5)])]))
        self.assertEqual([p["Expected"] for p in small["ParallelPairs"]], ["SameConductorGap", "SameConductorGap"])
        self.assertEqual(small["StandaloneCornerCount"], 0)

    def test_rounded_island_has_tangent_arcs_and_no_sharp_corners(self):
        lay = S.layout("island", [S.sheet(S.GROUND, S.rounded_rectangle(4.0, 3.0, 0.5))])
        orc = S.oracle(lay)
        self.assertAlmostEqual(orc["PhysicalPerimeterLength"], 24.0 + math.pi)
        self.assertEqual(orc["ClassifierCornerCount"], 0)
        self.assertEqual(len(orc["Arcs"]), 4)
        self.assertTrue(all(a["ExpectedRoundedCorner"] and a["Convex"] for a in orc["Arcs"]))
        self.assertEqual(orc["Arcs"][0]["Expected"], "ConvexCorner@90.0 radius 0.5")
        aperture = S.oracle(S.layout("aperture", [S.sheet(S.GROUND, S.rectangle(-12.0, -12.0, 12.0, 12.0), holes=[S.rounded_rectangle(4.0, 3.0, 0.5)])], half_x=12.0, half_y=12.0))
        self.assertTrue(all(not a["Convex"] for a in aperture["Arcs"]))
        self.assertAlmostEqual(aperture["PhysicalPerimeterLength"], 24.0 + math.pi)

    def test_arc_bar_polyline_turns_and_offset_width(self):
        lp = S.arc_bar(3.0, 20.0, 90.0, 5.0)
        points = np.array(lp["Points"])
        n = len(points) // 2
        right, left = points[:n], points[n:][::-1]
        # The offset polylines are 3 um apart at the lead ends (width) and every opposite
        # pair of vertices lies on a common normal through the centreline vertex.
        self.assertTrue(np.allclose(np.linalg.norm(right[[0, -1]] - left[[0, -1]], axis=1), 3.0))
        centre = np.array([0.0, 20.0])
        radii_right = np.linalg.norm(right[1:-1] - centre, axis=1)
        radii_left = np.linalg.norm(left[1:-1] - centre, axis=1)
        self.assertTrue(np.all(radii_right > radii_left))
        orc = S.oracle(S.layout("arc", [S.sheet(S.GROUND, lp)], half_x=40.0, half_y=40.0))
        # Chord polyline of a 90 deg arc in 5 deg steps: 2.5 deg where the chords meet the
        # straight leads, 5 deg between chords, 90 deg at the four lead ends.
        turns = sorted({round(c["TurnDegrees"], 6) for c in orc["Corners"]})
        self.assertEqual(turns, [2.5, 5.0, 90.0])
        self.assertEqual(orc["SubThresholdTurns"], 2 * (17 + 2))
        self.assertEqual(orc["ClassifierCornerCount"], 4)

    def test_arc_bar_sides_pair_along_the_bend_with_the_design_curvature_class(self):
        # Curved-edge chain rule: the two sides of a constant-width polyline bend are a pair
        # (separation constant within 5 %: 3 / cos(10 deg) at 20 deg per vertex), their
        # cross-chord interactions are not events (the two bar ends are the only cores), and
        # the expected class follows the inner-side design radius vs 20 R: r = 5 (1.75 R)
        # curved, r = 50 (24.25 R) straight-like.
        for radius, sweep, step, curved in [(5.0, 90.0, 20.0, True), (50.0, 45.0, 20.0, False), (250.0, 15.0, 5.0, False)]:
            lay = S.layout("arc", [S.sheet(S.GROUND, S.arc_bar(3.0, radius, sweep, step))], half_x=300.0, half_y=300.0, bend={"Radius": radius, "Width": 3.0})
            orc = S.oracle(lay)
            self.assertEqual(len(orc["BentPairs"]), 1, (radius, step))
            record = orc["BentPairs"][0]
            self.assertEqual(record["Curved"], curved)
            self.assertLessEqual(record["MaxSeparation"] - record["MinSeparation"], S.PAIR_SEPARATION_TOLERANCE * record["MinSeparation"])
            self.assertEqual(record["ExpectedClasses"], ["CurvedSameConductorStrip", "SameConductorStrip"] if curved else ["SameConductorStrip"])
            # Every chord-wise parallel pair belongs to the bent pair; the end corners are
            # cluster members (two corners 3 um apart), no other corner exists.
            self.assertTrue(all(p["InBentPair"] for p in orc["ParallelPairs"]))
            self.assertEqual(orc["CornerPairsWithin2R"], 2)
            self.assertEqual(orc["StandaloneCornerCount"], 0)

    def test_acute_corner_arms_are_not_a_pair_along_a_bend(self):
        # The arms of a 30 deg corner come within 2R of each other with a separation growing
        # along the arm: no constant separation, hence events (a spatial cluster), not a pair.
        lay = S.layout("corner", [S.sheet(S.GROUND, S.bent_bar(6.0, 16.0, 30.0))])
        orc = S.oracle(lay)
        self.assertEqual(orc["BentPairs"], [])

    def test_facing_layer_and_wall_are_excluded_classes(self):
        lay = S.layout("facing", [S.sheet(S.GROUND, S.rectangle(-10.0, -6.0, 10.0, 6.0)), S.sheet(S.GROUND, S.rectangle(-10.0, -6.0, 10.0, 6.0), z=3.0)], walls=[(S.GROUND, 0.0, -6.0, 0.0, 6.0, 4.0)])
        orc = S.oracle(lay)
        self.assertEqual(orc["Excluded"], {"CrossLayerSheets": 1, "Walls": 1})
        self.assertAlmostEqual(orc["PhysicalPerimeterLength"], 64.0)


class SpecificationTest(unittest.TestCase):
    def test_writer_emits_every_statement_kind(self):
        lay = S.layout("demo", [S.sheet(S.GROUND, S.rounded_rectangle(4.0, 3.0, 0.5)), S.sheet(S.GROUND, S.rectangle(-12.0, -12.0, 12.0, 12.0), holes=[S.rectangle(-1.0, -1.0, 1.0, 1.0)], z=3.0)], walls=[(S.GROUND, 0.0, -6.0, 0.0, 6.0, 4.0)], half_x=12.0, lc_fine=0.25)
        with tempfile.TemporaryDirectory() as directory:
            path = S.write_specification([lay], os.path.join(directory, "spec.txt"))
            with open(path) as source:
                text = source.read()
        self.assertIn("layout demo\nbox 12.0 30.0 15.0 15.0\nsize 0.25 6.0\npolygon 4\n", text)
        self.assertEqual(text.count("arc "), 4)
        self.assertIn("sheet 4 3.0\n", text)
        self.assertIn("hole\n", text)
        self.assertIn("wall 4 0.0 -6.0 0.0 6.0 4.0\n", text)
        self.assertTrue(text.endswith("end\n"))

    def test_suite_names_are_unique_and_cover_the_plan(self):
        layouts = S.stress_suite()
        names = [lay["Name"] for lay in layouts]
        self.assertEqual(len(names), len(set(names)))
        for prefix, count in (("gap-same-", 9), ("gap-different-", 9), ("strip-", 9), ("corner-", 10), ("hole-", 3), ("arc-", 12)):
            self.assertEqual(sum(1 for n in names if n.startswith(prefix)), count, prefix)
        for name in ("tee-stem3", "tee-stem6", "cross-arm3", "edge-to-boundary", "island-3x3", "taper-10", "facing-layers", "vertical-wall", "island-rounded-8x6", "aperture-rounded-8x6"):
            self.assertIn(name, names)


if __name__ == "__main__":
    unittest.main()
