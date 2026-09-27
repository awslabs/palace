# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""placement_check: the A10 placement gates on synthetic manifests whose patch frames are
known (a placed model passes; a mirrored / offset frame — the 2394fdb0c class — fails)."""
import math
import os
import sys
import unittest

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import placement_check as PC  # noqa: E402

R = 2.0
LAW = '{"Type":"PEC"}'


def frame(origin, angle_degrees, mirror=False):
    """Patch frame rows U, V, W: U at angle_degrees in the x-y plane, W = +z (V = W x U,
    negated for a mirrored frame)."""
    a = math.radians(angle_degrees)
    u = np.asarray([math.cos(a), math.sin(a), 0.0])
    w = np.asarray([0.0, 0.0, 1.0])
    v = np.cross(w, u) * (-1.0 if mirror else 1.0)
    return np.asarray(origin, dtype=float), np.asarray([u, v, w])


def patch(feature, model, origin, axes, segment=-1, index=0):
    return {"Patch": index, "Feature": feature, "Topology": "x", "Model": model, "ModelIndex": 0, "Weight": 1.0, "ModelWeight": 1.0,
            "QuadratureWeight": 1.0, "SideFactor": 1.0, "CouponDepth": 0.0, "Segment": segment, "S0": 0.0, "S1": 0.0, "Origin": origin, "Axes": axes}


class FingerEndCluster(unittest.TestCase):
    """A rounded finger end in the mesh (end edge at x = 10 facing +x, two 0.5 R fillets, R
    along both sides), placed at a 90 deg rotated, translated frame; the model is the
    canonical-frame description (end edge along y at x = 0) with chorded arcs. The lower
    side is 0.5 R longer, so a mirrored frame is a different placement."""

    def setUp(self):
        # Canonical geometry (units of R): end edge x = 0, |y| <= 0.25; fillet centres
        # (-0.5, +-0.25) radius 0.5; sides y = +-0.75 for x in [-1.5, -0.5] (top) and
        # [-2.0, -0.5] (bottom).
        c45 = 0.5 * math.cos(math.pi / 4)
        self.signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 5, "Portions": [
            {"P": [0.0, -0.25, 0.0, 0.25], "Gap": [1.0, 0.0], "Conductor": 1, "Interfaces": ["SA"], "Law": LAW},
            {"P": [-0.5, 0.75, 0.0, 0.25], "Arc": [-0.5, 0.25, -0.5 + c45, 0.25 + c45], "GapRadial": 1, "Conductor": 1, "Interfaces": ["SA"], "Law": LAW},
            {"P": [-0.5, -0.75, 0.0, -0.25], "Arc": [-0.5, -0.25, -0.5 + c45, -0.25 - c45], "GapRadial": 1, "Conductor": 1, "Interfaces": ["SA"], "Law": LAW},
            {"P": [-1.5, 0.75, -0.5, 0.75], "Gap": [0.0, 1.0], "Conductor": 1, "Interfaces": ["SA"], "Law": LAW},
            {"P": [-2.0, -0.75, -0.5, -0.75], "Gap": [0.0, -1.0], "Conductor": 1, "Interfaces": ["SA"], "Law": LAW}], "Vertices": []}
        # The feature in the mesh: the canonical frame rotated by 90 deg about +z and moved
        # to (10, 20): mesh = origin + x U + y V.
        self.origin, self.axes = frame((10.0, 20.0, 0.0), 90.0)
        segments, arcs, portions = [], [], []

        def mesh(p):
            return (self.origin + np.asarray([p[0] * R, p[1] * R, 0.0]) @ self.axes).tolist()

        def add_segment(a, b, arc=None):
            length = R * math.dist(a, b)
            entry = {"Key": [mesh(a), mesh(b)], "Length": length, "Chain": 0}
            if arc is not None:
                entry["Arc"] = arc
            segments.append(entry)
            portions.append([len(segments) - 1, 0.0, length])

        add_segment((0.0, -0.25), (0.0, 0.25))
        for sign, arc_index in ((1.0, 0), (-1.0, 1)):
            center = (-0.5, sign * 0.25)
            arcs.append({"Center": mesh(center), "Radius": 0.5 * R})
            # 6 mesh chords per quarter circle (15 deg), from the end edge to the side.
            for k in range(6):
                t0, t1 = sign * math.radians(15.0 * k), sign * math.radians(15.0 * (k + 1))
                add_segment((center[0] + 0.5 * math.cos(t0), center[1] + 0.5 * math.sin(t0)),
                            (center[0] + 0.5 * math.cos(t1), center[1] + 0.5 * math.sin(t1)), arc=arc_index)
        add_segment((-0.5, 0.75), (-1.5, 0.75))
        add_segment((-0.5, -0.75), (-2.0, -0.75))  # the longer side: the cluster is chiral
        self.identification = {"MatchingRadius": R, "Segments": segments, "Arcs": arcs, "Features": [
            {"Id": 3, "Type": "SpatialEdgeCluster", "Signature": self.signature, "Hash": "h", "Chirality": 1, "Length": 1.0, "Portions": portions,
             "Frame": {"Origin": self.origin.tolist(), "Axes": self.axes.tolist()}, "Match": {"Status": "Matched", "Model": "finger"}}]}
        self.library = {"MatchingRadius": R, "Models": [{"Name": "finger", "Topology": "SpatialEdgeCluster", "Signature": self.signature}]}

    def test_signature_only_model_in_the_feature_frame_passes(self):
        gates, summary = PC.placement_gates(self.identification, [patch(3, "finger", self.origin, self.axes)], self.library, R)
        by_name = {g["Gate"]: g for g in gates}
        self.assertEqual(by_name["A10-placement-clusters"]["Status"], "PASS")
        self.assertLess(by_name["A10-placement-clusters"]["Detail"]["WorstDeviationOverR"], 1.0e-9)
        self.assertEqual(by_name["A10-placement-clusters"]["Detail"]["Features"], 1)
        self.assertEqual(by_name["A10-placement-corners"]["Detail"]["Features"], 0)  # nothing to place: PASS

    def test_model_edges_scaled_from_library_units_pass(self):
        # Library in other length units (R = 1.9): Edges = the chorded portions x 1.9.
        library = {"MatchingRadius": 1.9, "Models": [{"Name": "finger", "Topology": "SpatialEdgeCluster", "Signature": self.signature, "Edges": []}]}
        from . import signature_library
        for edge in signature_library.cluster_plan_view_edges(self.signature, 1.9):
            gap = np.asarray(edge["Gap"]) / np.linalg.norm(edge["Gap"])
            tangent = np.asarray([gap[1], -gap[0]])
            p0, p1 = np.asarray(edge["P0"]), np.asarray(edge["P1"])
            forward = float(np.dot(p1 - p0, tangent)) > 0.0
            length = float(np.linalg.norm(p1 - p0))
            library["Models"][0]["Edges"].append({"Point": [p0[0], p0[1], 0.0], "GapDirection": [gap[0], gap[1], 0.0], "ProcessNormal": [0.0, 0.0, 1.0],
                                                  "Interval": [0.0, length] if forward else [-length, 0.0], "Conductor": 1})
        gates, _ = PC.placement_gates(self.identification, [patch(3, "finger", self.origin, self.axes)], library, R)
        cluster = next(g for g in gates if g["Gate"] == "A10-placement-clusters")
        self.assertEqual(cluster["Status"], "PASS")
        self.assertLess(cluster["Detail"]["WorstDeviationOverR"], 1.0e-9)

    def test_mirrored_or_offset_frame_fails(self):
        mirrored_origin, mirrored_axes = frame((10.0, 20.0, 0.0), 90.0, mirror=True)
        gates, _ = PC.placement_gates(self.identification, [patch(3, "finger", mirrored_origin, mirrored_axes)], self.library, R)
        cluster = next(g for g in gates if g["Gate"] == "A10-placement-clusters")
        self.assertEqual(cluster["Status"], "FAIL")
        self.assertGreater(cluster["Detail"]["WorstDeviationOverR"], 0.4)
        offset_origin = self.origin + np.asarray([0.0, 0.01 * R, 0.0])
        gates, _ = PC.placement_gates(self.identification, [patch(3, "finger", offset_origin, self.axes)], self.library, R)
        cluster = next(g for g in gates if g["Gate"] == "A10-placement-clusters")
        self.assertEqual(cluster["Status"], "FAIL")
        # A rotation by the chord step (5 deg) about the origin moves every point off the arcs
        # and the straight portions.
        _, rotated_axes = frame((10.0, 20.0, 0.0), 95.0)
        gates, _ = PC.placement_gates(self.identification, [patch(3, "finger", self.origin, rotated_axes)], self.library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-clusters")["Status"], "FAIL")


class CornerStackPair(unittest.TestCase):
    def test_sharp_corner(self):
        # A 90 deg convex corner at (5, 5): arms toward +x and +y (R long); the feature frame
        # x = first arm (+x), y = +y.
        origin, axes = frame((5.0, 5.0, 0.0), 0.0)
        segments = [{"Key": [[5.0, 5.0, 0.0], [5.0 + R, 5.0, 0.0]], "Length": R, "Chain": 0},
                    {"Key": [[5.0, 5.0, 0.0], [5.0, 5.0 + R, 0.0]], "Length": R, "Chain": 0}]
        identification = {"MatchingRadius": R, "Segments": segments, "Arcs": [], "Features": [
            {"Id": 1, "Type": "ConvexCorner", "Signature": {"AngleDegrees": 90.0}, "Hash": "c", "Chirality": 1, "Length": 2 * R,
             "Portions": [[0, 0.0, R], [1, 0.0, R]], "Frame": {"Origin": origin.tolist(), "Axes": axes.tolist()}, "Match": {"Status": "Matched", "Model": "corner"}}]}
        library = {"MatchingRadius": 1.9, "Models": [{"Name": "corner", "Topology": "ConvexCorner", "Angle": 90.0, "CornerRadius": 0.0}]}
        gates, _ = PC.placement_gates(identification, [patch(1, "corner", origin, axes)], library, R)
        corner = next(g for g in gates if g["Gate"] == "A10-placement-corners")
        self.assertEqual(corner["Status"], "PASS")
        # The frame's first arm along +y (the arms swapped): the second arm would be at -x.
        _, swapped = frame((5.0, 5.0, 0.0), 90.0)
        gates, _ = PC.placement_gates(identification, [patch(1, "corner", origin, swapped)], library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-corners")["Status"], "FAIL")
        # A 60 deg model on the 90 deg feature.
        library["Models"][0]["Angle"] = 60.0
        gates, _ = PC.placement_gates(identification, [patch(1, "corner", origin, axes)], library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-corners")["Status"], "FAIL")

    def test_rounded_corner_arc_on_the_fillet_circle(self):
        # A 90 deg corner rounded at r = 0.5 R: virtual corner at the origin, tangent points
        # at (r, 0) and (0, r), fillet centre (r, r); two mesh chords on the arc.
        r = 0.5 * R
        origin, axes = frame((0.0, 0.0, 0.0), 0.0)
        mid = (r + r * math.cos(math.radians(225.0)), r + r * math.sin(math.radians(225.0)))
        segments = [{"Key": [[r, 0.0, 0.0], [R, 0.0, 0.0]], "Length": R - r, "Chain": 0},
                    {"Key": [[0.0, r, 0.0], [0.0, R, 0.0]], "Length": R - r, "Chain": 0},
                    {"Key": [[r, 0.0, 0.0], [mid[0], mid[1], 0.0]], "Length": math.dist((r, 0.0), mid), "Chain": 0, "Arc": 0},
                    {"Key": [[mid[0], mid[1], 0.0], [0.0, r, 0.0]], "Length": math.dist((r, 0.0), mid), "Chain": 0, "Arc": 0}]
        identification = {"MatchingRadius": R, "Segments": segments, "Arcs": [{"Center": [r, r, 0.0], "Radius": r}], "Features": [
            {"Id": 1, "Type": "ConvexCorner", "Signature": {"AngleDegrees": 90.0, "CornerRadiusOverR": 0.5}, "Hash": "c", "Chirality": 1, "Length": 1.0,
             "Portions": [[s, 0.0, segments[s]["Length"]] for s in range(4)], "Frame": {"Origin": origin.tolist(), "Axes": axes.tolist()},
             "Match": {"Status": "Matched", "Model": "rounded"}}]}
        library = {"MatchingRadius": R, "Models": [{"Name": "rounded", "Topology": "ConvexCorner", "Angle": 90.0, "CornerRadius": r}]}
        gates, _ = PC.placement_gates(identification, [patch(1, "rounded", origin, axes)], library, R)
        corner = next(g for g in gates if g["Gate"] == "A10-placement-corners")
        self.assertEqual(corner["Status"], "PASS")
        library["Models"][0]["CornerRadius"] = 0.0
        gates, _ = PC.placement_gates(identification, [patch(1, "rounded", origin, axes)], library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-corners")["Status"], "FAIL")

    def test_stack_and_pair_sides(self):
        # Three parallel edges along +x at y = 0, 1, 2.5 (offsets 0, 0.5 R, 1.25 R), one
        # longitudinal patch per side with the origin on side 0 and U = +y; chirality -1
        # reverses the side order.
        segments = [{"Key": [[0.0, y, 0.0], [10.0, y, 0.0]], "Length": 10.0, "Chain": 0} for y in (0.0, 1.0, 2.5)]
        sides = [0, 1, 2]
        feature = {"Id": 7, "Type": "ParallelEdgeCluster", "Signature": {}, "Hash": "s", "Chirality": 1, "Length": 30.0,
                   "Portions": [[s, 0.0, 10.0] for s in range(3)], "Sides": sides, "Frame": {"Origin": [0.0, 0.0, 0.0], "Axes": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]},
                   "Match": {"Status": "Matched", "Model": "stack"}}
        identification = {"MatchingRadius": R, "Segments": segments, "Arcs": [], "Features": [feature]}
        library = {"MatchingRadius": 1.9, "Models": [{"Name": "stack", "Topology": "ParallelEdgeCluster",
                                                      "Edges": [{"Offset": 0.0}, {"Offset": 0.5 * 1.9}, {"Offset": 1.25 * 1.9}]}]}
        _, axes = frame((0.0, 0.0, 0.0), 90.0)   # U = +y (increasing offset), W = +z
        patches = [patch(7, "stack", np.asarray([x, 0.0, 0.0]), axes, segment=0, index=i) for i, x in enumerate((2.5, 5.0, 7.5))]
        gates, _ = PC.placement_gates(identification, patches, library, R)
        stack = next(g for g in gates if g["Gate"] == "A10-placement-stacks")
        self.assertEqual(stack["Status"], "PASS")
        self.assertEqual(stack["Detail"]["Checks"], 9)
        # Mirrored feature (chirality -1): the patch origin is on the top side, U = -y.
        feature["Chirality"] = -1
        _, down = frame((0.0, 0.0, 0.0), -90.0)
        library["Models"][0]["Edges"] = [{"Offset": 0.0}, {"Offset": 0.75 * 1.9}, {"Offset": 1.25 * 1.9}]
        patches = [patch(7, "stack", np.asarray([x, 2.5, 0.0]), down, segment=2, index=i) for i, x in enumerate((2.5, 5.0))]
        gates, _ = PC.placement_gates(identification, patches, library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-stacks")["Status"], "PASS")
        # An offset off by 6 % of itself (beyond the pair rule's 5 %) fails; 2e-3 R (0.3 %) does not.
        library["Models"][0]["Edges"][1]["Offset"] = (0.75 + 0.002) * 1.9
        gates, _ = PC.placement_gates(identification, patches, library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-stacks")["Status"], "PASS")
        library["Models"][0]["Edges"][1]["Offset"] = 0.75 * 1.06 * 1.9
        gates, _ = PC.placement_gates(identification, patches, library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-stacks")["Status"], "FAIL")
        # A pair: sides y = 0 and y = 1, separation 0.5 R, origin at the midline.
        pair = {"Id": 8, "Type": "SameConductorGap", "Signature": {"SeparationOverR": 0.5}, "Hash": "p", "Chirality": 1, "Length": 20.0,
                "Portions": [[0, 0.0, 10.0], [1, 0.0, 10.0]], "Sides": [0, 1], "Frame": feature["Frame"], "Match": {"Status": "Matched", "Model": "pair"}}
        identification["Features"] = [pair]
        library = {"MatchingRadius": R, "Models": [{"Name": "pair", "Topology": "SameConductorGap", "Separation": 1.0}]}
        patches = [patch(8, "pair", np.asarray([4.0, 0.5, 0.0]), axes, segment=0, index=0)]
        gates, _ = PC.placement_gates(identification, patches, library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-pairs")["Status"], "PASS")
        library["Models"][0]["Separation"] = 1.04  # within the pair rule's 5 %
        gates, _ = PC.placement_gates(identification, patches, library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-pairs")["Status"], "PASS")
        library["Models"][0]["Separation"] = 1.2
        gates, _ = PC.placement_gates(identification, patches, library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-pairs")["Status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
