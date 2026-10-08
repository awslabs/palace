# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""generate_spatial_response: the conductor labelling of the trace knots against the plan-view
mask (decision 360 (c*): a tolerance of 4 x the 1e-9 R plan quantum, inclusive, fail-closed to
PEC), the labelling validator (no FREE knot within the tolerance of a metal outline; the near-
outline inspection record) and the mask-frame check (design round 2 review MINOR-5)."""
from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cluster_signature_geometry as csg  # noqa: E402
import generate_spatial_response as gsr  # noqa: E402

RADIUS = 1.9
THICKNESS = 0.1
# Conductor 1 covers x in [-1, 0], y in [-1, 1] on the plane z = 0 (two mask facets).
FACETS = [{"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [0.0, -1.0], [0.0, 1.0]]},
          {"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [0.0, 1.0], [-1.0, 1.0]]}]
EDGES = [{"Conductor": 1, "Point": [-1.0, 0.0, 0.0], "ProcessNormal": [0.0, 0.0, 1.0], "Tangent": [0.0, 1.0, 0.0],
          "GapDirection": [-1.0, 0.0, 0.0], "Interval": [-1.0, 1.0], "VertexArm": False}]
QUANTUM = 1.0e-9 * RADIUS


def labels_of(points):
    return gsr.conductor_at_points(np.asarray(points, dtype=float), EDGES, RADIUS, THICKNESS, 90.0, FACETS)


class ConductorLabelToleranceTest(unittest.TestCase):
    def test_tolerance_is_four_plan_quanta_inclusive(self):
        self.assertEqual(gsr.CONDUCTOR_LABEL_TOLERANCE_OVER_R, 4.0e-9)
        inside = [-0.5, 0.0, 0.0]
        on_outline = [-1.0, 0.5, 0.0]
        rounding = [-1.0 - 0.3 * QUANTUM, -1.0 - 0.3 * QUANTUM, THICKNESS]  # a corner knot ~0.4 quanta outside (M2)
        three_quanta = [-1.0 - 3.0 * QUANTUM, 0.5, 0.0]
        four_quanta = [-1.0 - 4.0 * QUANTUM, 0.5, THICKNESS]                 # exactly the tolerance: inclusive
        five_quanta = [-1.0 - 5.0 * QUANTUM, 0.5, 0.0]
        above_metal = [-0.5, 0.0, 0.5]
        labels = labels_of([inside, on_outline, rounding, three_quanta, four_quanta, five_quanta, above_metal])
        self.assertEqual(labels.tolist(), [1, 1, 1, 1, 1, 0, 0])

    def test_edge_only_signature_without_facets_is_unchanged(self):
        # The edge-only signature puts the metal on the side opposite the gap direction (x >= -1 here).
        points = np.array([[-1.0 - 3.0 * QUANTUM, 0.5, 0.0], [-0.5, 0.0, 0.0], [-1.0 - 0.5 * 1.0e-10 * RADIUS, 0.0, 0.0]])
        labels = gsr.conductor_at_points(points, EDGES, RADIUS, THICKNESS, 90.0, [])
        self.assertEqual(labels.tolist(), [0, 1, 1])  # tolerance 1e-10 R, no plan-view quantum involved


class NearOutlineFreeKnotsTest(unittest.TestCase):
    def test_validator_binds_the_labels_and_records_the_inspection_band(self):
        points = np.array([[-1.0 - 3.0 * QUANTUM, 0.5, 0.0],        # within the tolerance: must be labelled
                           [-1.0 - 5.0 * QUANTUM, 0.5, THICKNESS],  # free, inside the 1e-3 R inspection band
                           [-1.0 - 2.0e-3 * RADIUS, 0.5, 0.0],      # free, beyond the band
                           [-1.0 - 3.0 * QUANTUM, 0.5, 0.5],        # above the metal: not judged
                           [-0.5, 0.0, 0.0]])                       # in the metal
        labels = labels_of(points)
        self.assertEqual(labels.tolist(), [1, 0, 0, 0, 1])
        violations, inspect = gsr.near_outline_free_knots(points, labels, EDGES, RADIUS, THICKNESS, FACETS)
        self.assertEqual(violations, [])
        self.assertEqual([k["Vertex"] for k in inspect], [2])
        self.assertAlmostEqual(inspect[0]["Distance"], 5.0 * QUANTUM, places=15)
        # A label set that leaves a knot within the tolerance free is a violation (fail closed).
        stale = labels.copy()
        stale[0] = 0
        violations, inspect = gsr.near_outline_free_knots(points, stale, EDGES, RADIUS, THICKNESS, FACETS)
        self.assertEqual([k["Vertex"] for k in violations], [1])
        self.assertEqual([k["Vertex"] for k in inspect], [2])
        self.assertEqual(gsr.near_outline_free_knots(points, stale, EDGES, RADIUS, THICKNESS, []), ([], []))


class MaskFrameTest(unittest.TestCase):
    def test_mask_inside_the_box_passes_and_a_foreign_frame_fails_closed(self):
        lower, upper = np.array([-1.0, -1.0, -2.0]), np.array([3.0, 1.0, 2.0])
        self.assertEqual(gsr.validate_mask_frame(FACETS, lower, upper, RADIUS), [-1.0, -1.0, 0.0, 1.0])
        self.assertIsNone(gsr.validate_mask_frame([], lower, upper, RADIUS))
        slack = gsr.CONDUCTOR_LABEL_TOLERANCE_OVER_R * RADIUS
        self.assertIsNotNone(gsr.validate_mask_frame(FACETS, lower + [0.5 * slack, 0.0, 0.0], upper, RADIUS))
        rotated = [{"Conductor": 1, "Plane": 0.0, "Points": [[2.0, 3.0], [12.0, 3.0], [12.0, 9.0]]}]
        with self.assertRaises(ValueError) as caught:
            gsr.validate_mask_frame(rotated, lower, upper, RADIUS)
        self.assertIn("not in the mesh frame", str(caught.exception))
        with self.assertRaises(ValueError):
            gsr.validate_mask_frame(FACETS, lower + [2.0 * slack, 0.0, 0.0], upper, RADIUS)

    def test_generator_stops_before_any_output_on_a_foreign_frame(self):
        """Decision 369 review MINOR-2: the frame check runs before the first output file, so a
        mis-framed mask leaves generation-failure.json (Stage PlanViewMaskFrame) and nothing
        else: no partial mesh-signature.csv next to it."""
        law = '{"Type":"PEC"}'
        portions = [{"Conductor": 1, "Gap": [0.0, 1.0], "Interfaces": ["MA", "MS", "SA"], "Law": law, "P": [-2.0, 0.0, 0.0, 0.0]},
                    {"Conductor": 1, "Gap": [-1.0, 0.0], "Interfaces": ["MA", "MS", "SA"], "Law": law, "P": [0.0, 0.0, 0.0, 2.0]}]
        signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 2, "Portions": portions,
                     "Vertices": [{"P": [0.0, 0.0], "TurnDegrees": 90.0, "Type": "ConcaveCorner"}]}
        cluster = {"Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": 2, "Signature": signature}, "Signature": signature,
                   "Interfaces": [{"Slot": 0, "Type": t, "Target": n} for n, t in ((3, "MA"), (2, "MS"), (1, "SA"))],
                   "BoundaryCondition": {"Type": "PEC"}}
        coupon, _ = csg.cluster_coupon(cluster, RADIUS, THICKNESS, 0.05)
        for facet in coupon["Geometry"]["PlanViewFacets"]:  # the mask of another frame: shifted by 1 um
            facet["Points"] = [[p[0] + 1.0, p[1]] + list(p[2:]) for p in facet["Points"]]
        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "coupon.json").write_text(json.dumps(coupon))
            result = subprocess.run([sys.executable, str(HERE / "generate_spatial_response.py"), str(Path(tmp) / "coupon.json"),
                                     "--output", str(Path(tmp) / "out"), "--radius", str(RADIUS), "--metal-thickness",
                                     str(THICKNESS), "--overetch-depth", "0.05", "--sidewall-angle", "90", "--top-rounding", "0",
                                     "--trench-rounding", "0", "--model-name", "m", "--basis-only"],
                                    capture_output=True, text=True, timeout=120)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("not in the mesh frame", result.stderr)
            failure = json.loads((Path(tmp) / "out" / "generation-failure.json").read_text())
            self.assertEqual(failure["Stage"], "PlanViewMaskFrame")
            self.assertEqual(sorted(p.name for p in (Path(tmp) / "out").iterdir()), ["generation-failure.json"])


CENSUS_B = HERE / "testdata" / "census-b-signatures" / "signatures.json"


def census_loops(prefix, radius=RADIUS):
    """main's path on a census-b signature: cluster_coupon -> normalize_geometry -> the boundary
    loops with the protected set of signature_arcs -> tag_arc_boundary_loops from the same arcs."""
    record = json.loads(CENSUS_B.read_text())[prefix]
    coupon, _ = csg.cluster_coupon(record, radius, THICKNESS, 0.05)
    frame, edges, facets = gsr.normalize_geometry(coupon, radius)
    lower, upper = gsr.coupon_bounds(edges, radius, THICKNESS, 0.05, coupon["Geometry"].get("SupportBox"))
    arcs = gsr.signature_arcs(coupon, radius)
    loops = gsr.plan_view_boundary_loops(facets, radius, lower, upper,
                                         gsr.classified_continuation_segments(coupon["Geometry"], frame, radius),
                                         protected=gsr.protected_arc_vertex_keys(arcs, frame, radius))
    moved = gsr.reconcile_mask_with_boundary(facets, loops, radius)
    count = gsr.tag_arc_boundary_loops(loops, coupon, frame, radius, arcs=arcs)
    return coupon, arcs, loops, moved, count


class ProtectedArcVerticesTest(unittest.TestCase):
    """Round 3 classes (2) / (3), fix 2b = R1 (DESIGN-part-M 1.3, part G G.3.3 (i); decision 510):
    the near-collinear merge of plan_view_boundary_loops never pops a chord vertex or an end of
    a kept rebuilt arc; an unprotected vertex within JUNCTION_TANGENT_ANGLE is popped as before
    (every straight loop bitwise)."""

    def facets_with_kink(self, kink):
        # Conductor 1 covers [-1, 0] x [-1, 1] with the right side bent at (kink, 0): a turn of
        # ~2 kink rad between the two halves of that side.
        return [{"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [0.0, -1.0], [kink, 0.0]]},
                {"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [kink, 0.0], [0.0, 1.0]]},
                {"Conductor": 1, "Plane": 0.0, "Points": [[-1.0, -1.0], [0.0, 1.0], [-1.0, 1.0]]}]

    def test_unprotected_small_turn_is_popped_and_a_protected_one_kept(self):
        lower, upper = np.array([-3.0, -3.0, -1.0]), np.array([3.0, 3.0, 1.0])
        facets = self.facets_with_kink(5.0e-6)             # turn 1e-5 rad <= JUNCTION_TANGENT_ANGLE
        loops = gsr.plan_view_boundary_loops(facets, RADIUS, lower, upper)
        self.assertEqual(len(loops), 1)
        self.assertEqual(len(loops[0]["Points"]), 4)        # popped: the straight loop of today
        self.assertEqual(loops, gsr.plan_view_boundary_loops(facets, RADIUS, lower, upper, protected=frozenset()))
        kink_key = gsr._arc_vertex_key((5.0e-6, 0.0), RADIUS)
        kept = gsr.plan_view_boundary_loops(facets, RADIUS, lower, upper, protected=frozenset([kink_key]))
        self.assertEqual(len(kept[0]["Points"]), 5)
        self.assertIn(kink_key, {gsr._arc_vertex_key(p, RADIUS) for p in kept[0]["Points"]})
        # A protected key elsewhere protects nothing; a turn above the angle is kept either way.
        other = gsr.plan_view_boundary_loops(facets, RADIUS, lower, upper,
                                             protected=frozenset([gsr._arc_vertex_key((0.5, 0.5), RADIUS)]))
        self.assertEqual(len(other[0]["Points"]), 4)
        corner = gsr.plan_view_boundary_loops(self.facets_with_kink(1.0e-3), RADIUS, lower, upper)
        self.assertEqual(len(corner[0]["Points"]), 5)

    def test_straight_coupons_have_no_protected_keys(self):
        self.assertEqual(gsr.signature_arcs({"Topology": "SpatialEdgeCluster", "Geometry": {}}, RADIUS), [])
        self.assertEqual(gsr.protected_arc_vertex_keys([], np.identity(3), RADIUS), frozenset())

    def test_efe678516aa0_same_circle_joint_vertex_is_kept_and_tagged(self):
        # Class (2) (part M 1.1): the joint vertex of claim arc 2 and context arc 6 (one circle of
        # 2240 R; chord-to-chord turn 3.26e-5 rad at a signature turn of 1.49e-8) was popped by the
        # merge and tag_arc_boundary_loops failed closed ("chord 8 is not exactly one side"); under
        # the protected set every chord of every arc is a loop side.
        coupon, arcs, loops, moved, count = census_loops("efe678516aa0")
        self.assertEqual(count, 6)
        tagged = sum(arc is not None for loop in loops for arc in loop["Arcs"])
        self.assertEqual(tagged, sum(len(arc["Vertices"]) - 1 for arc in arcs))
        # The merge still pops the device plan's straight-straight pseudo-corners (the 1e-6 R
        # grid's sub-1e-4 turns, reconciled into the mask as before); no arc vertex is among them.
        protected = gsr.protected_arc_vertex_keys(arcs, np.identity(3), RADIUS)
        self.assertGreater(len(moved), 0)
        self.assertFalse({gsr._arc_vertex_key(source, RADIUS) for source, _ in moved} & protected)

    def test_baf9dacceb51_demoted_context_arcs_are_straight_sides_of_the_boundary(self):
        # Class (3) (part G G.3.1): the C3 trio's context arcs 9 / 10 (R_arc 1.1e5 .. 1.4e5 R) are
        # DEMOTED by rebuilt_arcs (fix (3)(ii)): no tag, their chords straight boundary sides; the
        # eight device arcs are protected and tagged. The generator of record stops earlier at the
        # MINOR-5 margin (test_cluster_signature_geometry.ChainAndDemotionTest); the boundary is
        # read here with the bound lifted (a probe of the geometry).
        bound = csg.ARC_SMOOTH_JOINT_TURN_BOUND
        try:
            csg.ARC_SMOOTH_JOINT_TURN_BOUND = 1.0
            coupon, arcs, loops, moved, count = census_loops("baf9dacceb51")
        finally:
            csg.ARC_SMOOTH_JOINT_TURN_BOUND = bound
        self.assertEqual(count, 8)
        self.assertEqual([arc["ArcId"] for arc in arcs if arc["Demoted"]], [9, 10])
        tagged = sum(arc is not None for loop in loops for arc in loop["Arcs"])
        self.assertEqual(tagged, sum(len(arc["Vertices"]) - 1 for arc in arcs if not arc["Demoted"]))
        self.assertFalse({tag["ArcId"] for loop in loops for tag in loop["Arcs"] if tag is not None} & {9, 10})
        self.assertEqual({tag["ArcChain"] for loop in loops for tag in loop["Arcs"] if tag is not None}, {1, 2, 3, 4})
        # A demoted arc's chord is a straight side: its ends are loop vertices.
        keys = {gsr._arc_vertex_key(p, RADIUS) for loop in loops for p in loop["Points"]}
        for arc in arcs:
            if arc["Demoted"]:
                self.assertTrue(all(gsr._arc_vertex_key(v, RADIUS) in keys for v in arc["Vertices"]))

    def test_boundary_csv_carries_the_arc_chain_column_only_when_tagged(self):
        # R7 (510 MINOR-6): ArcChain joins the tagged columns; a straight boundary's CSV is unchanged.
        self.assertEqual(gsr.ARC_BOUNDARY_COLUMNS[-1], "ArcChain")
        coupon, arcs, loops, _, count = census_loops("efe678516aa0")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "plan-view-boundary.csv"
            gsr.write_plan_view_boundary(path, loops)
            lines = path.read_text().splitlines()
            self.assertEqual(lines[0], "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y,ArcId,ArcCx,ArcCy,ArcR,ArcSign,JointTurn,JointSmooth,ArcChain")
            chains = {line.split(",")[-1] for line in lines[1:] if line.split(",")[8] != ""}
            self.assertEqual(chains, {"0", "1", "2"})
            self.assertTrue(all(line.endswith(",") for line in lines[1:] if line.split(",")[8] == ""))
            straight = [{"Conductor": 1, "Plane": 0.0, "Hole": False, "Classes": ["Physical"] * 4,
                         "Points": [[-1.0, -1.0], [0.0, -1.0], [0.0, 1.0], [-1.0, 1.0]]}]
            gsr.write_plan_view_boundary(path, straight)
            self.assertEqual(path.read_text().splitlines()[0], "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y")


if __name__ == "__main__":
    unittest.main()
