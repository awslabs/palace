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
# The PRODUCTION placement record of the rerun-3b J1-e osc003-O3 thin preflight: feature 90 =
# the mirror-formed cluster d7c875318447 placed on its REAL half (decision 557 (4)), whose IMAGE
# edges / context ends beyond the mirror plane x = -830 failed the not-mirror-aware A10 by
# construction (J1-E10 / F-R3b-4; decision 616 MINOR-8; v3-f E-4).
MIRROR_FORMED_FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "testdata", "mirror-formed-placement-o3-d7c875318447.json")
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

    def test_arc_aware_device_perimeter_distance_and_the_context_gate(self):
        # Block (b) DESIGN section 2 (a): a context entry END cut at a box face lies on the
        # FITTED circle, up to the chord sagitta (0.5 R x (1 - cos 7.5 deg) = 4.3e-3 R here)
        # from the device polyline: the chord reading refuses it, the arc reading accepts it.
        c = (-0.5, 0.25)
        t = math.radians(7.5)  # the middle of the first 15-deg mesh chord
        on_circle = self.origin + np.asarray([(c[0] + 0.5 * math.cos(t)) * R, (c[1] + 0.5 * math.sin(t)) * R, 0.0]) @ self.axes
        self.assertLess(PC.device_perimeter_distance(self.identification, on_circle), 1.0e-9 * R)
        chords_only = {**self.identification, "Segments": [{k: v for k, v in s.items() if k != "Arc"} for s in self.identification["Segments"]]}
        sagitta = 0.5 * R * (1.0 - math.cos(t))
        self.assertAlmostEqual(PC.device_perimeter_distance(chords_only, on_circle), sagitta, delta=1.0e-6 * R)
        # A point radially 0.01 R off the circle inside the arc's angular range reads 0.01 R;
        # one outside the chord's angular range reads the distance to the chord end.
        off = self.origin + np.asarray([(c[0] + 0.51 * math.cos(t)) * R, (c[1] + 0.51 * math.sin(t)) * R, 0.0]) @ self.axes
        self.assertAlmostEqual(PC.device_perimeter_distance(self.identification, off), 0.01 * R, delta=1.0e-9 * R)
        a = np.asarray([0.0, 0.0, 0.0]); b = np.asarray([1.0, 0.0, 0.0]); centre = np.asarray([0.5, -0.5, 0.0])
        rho = math.hypot(0.5, 0.5)
        beyond = centre + rho * np.asarray([math.cos(math.radians(170.0)), math.sin(math.radians(170.0)), 0.0])  # a is at 135 deg, b at 45
        self.assertAlmostEqual(PC._arc_chord_distance(beyond, a, b, centre, rho), min(np.linalg.norm(beyond - a), np.linalg.norm(beyond - b)), places=12)
        self.assertAlmostEqual(PC._arc_chord_distance(np.asarray([0.5, 0.2, 0.3]), a, b, centre, rho), math.hypot(0.7 - rho, 0.3), places=12)
        self.assertAlmostEqual(PC._arc_chord_distance(np.asarray([0.5, 0.2, 0.0]), a, b, np.asarray([2.0, 0.0, 0.0]), 1.0), 0.2, places=12)  # degenerate: the chord
        # The context gate: a contract-3 model whose Context entry ends at the arc's middle
        # (on the circle) passes; the same end moved 2e-3 R radially fails.
        def library(end):
            signature = dict(self.signature)
            signature["Box"] = [-3.0, -2.0, 1.0, 2.0]
            signature["Context"] = [{"P": [c[0] + 0.5 * math.cos(math.radians(45.0)), c[1] + 0.5 * math.sin(math.radians(45.0)), end[0], end[1]],
                                     "Arc": [c[0], c[1], c[0] + 0.5 * math.cos(math.radians(26.25)), c[1] + 0.5 * math.sin(math.radians(26.25))],
                                     "GapRadial": 1, "Conductor": 1, "Interfaces": ["SA"], "Law": LAW, "Chain": True}]
            return {"MatchingRadius": R, "Models": [{"Name": "finger", "Topology": "SpatialEdgeCluster", "Signature": signature}]}
        exact_end = (c[0] + 0.5 * math.cos(t), c[1] + 0.5 * math.sin(t))
        gates, summary = PC.placement_gates(self.identification, [patch(3, "finger", self.origin, self.axes)], library(exact_end), R)
        context = {g["Gate"]: g for g in gates}["A10-context"]
        self.assertEqual(context["Status"], "PASS")
        self.assertEqual(context["Detail"]["Checks"], 2)
        self.assertEqual(context["Detail"]["ArcEnds"], 2)
        self.assertLess(context["Detail"]["WorstDeviationOverR"], 1.0e-9)
        self.assertEqual(summary["Context"]["Features"], 1)
        moved_end = (c[0] + 0.502 * math.cos(t), c[1] + 0.502 * math.sin(t))
        gates, _ = PC.placement_gates(self.identification, [patch(3, "finger", self.origin, self.axes)], library(moved_end), R)
        context = {g["Gate"]: g for g in gates}["A10-context"]
        self.assertEqual(context["Status"], "FAIL")
        self.assertAlmostEqual(context["Detail"]["WorstDeviationOverR"], 2.0e-3, places=6)
        # A claims-only model (no Box): nothing to check.
        gates, _ = PC.placement_gates(self.identification, [patch(3, "finger", self.origin, self.axes)], self.library, R)
        self.assertEqual({g["Gate"]: g for g in gates}["A10-context"]["Status"], "PASS")
        self.assertEqual({g["Gate"]: g for g in gates}["A10-context"]["Detail"]["Features"], 0)

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

    def test_curved_pair_convexity_side(self):
        # A concentric curved gap (inner arc radius 3 R, outer 3.5 R, separation 0.5 R) on
        # 1-chord-per-5-deg polylines, side 0 = the inner arc; a curvature-family runtime
        # model `<anchor>@convex-kappa...` resolves to its anchor and its e1 (origin - s/2 U)
        # must sit on the INNER circle (Convex gap edge = metal inside the bend). With the
        # geometrically right frame (U outward, e1 on side 0) a coupon of the WRONG convexity
        # (the concave blend, the concave first-order node) is the defect the side checks
        # cannot see; a runtime model of an unknown anchor is not evaluable.
        centre = np.asarray([0.0, 0.0, 0.0])
        segments, arcs, portions, sides = [], [], [], []
        for side, radius in enumerate((3.0 * R, 3.5 * R)):
            arcs.append({"Center": centre.tolist(), "Radius": radius, "RadiusOverR": radius / R, "TurnDegrees": 40.0, "Kind": "Bend", "Joints": 8, "Segments": 8})
            for k in range(8):
                a0, a1 = math.radians(5 * k), math.radians(5 * (k + 1))
                p0 = [radius * math.cos(a0), radius * math.sin(a0), 0.0]
                p1 = [radius * math.cos(a1), radius * math.sin(a1), 0.0]
                length = float(np.linalg.norm(np.asarray(p1) - p0))
                segments.append({"Key": [p0, p1], "Length": length, "Chain": side, "Arc": side})
                portions.append([len(segments) - 1, 0.0, length])
                sides.append(side)
        feature = {"Id": 3, "Type": "CurvedSameConductorGap", "Signature": {"SeparationOverR": 0.5, "RadiusOverR": 3.0, "Convexity": "Convex"},
                   "Hash": "c", "Chirality": 0, "Length": 40.0, "Portions": portions, "Sides": sides,
                   "Frame": {"Origin": [0.0, 0.0, 0.0], "Axes": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]},
                   "Match": {"Status": "Matched", "Model": "gap@convex-kappa0.333333333-cubic"}}
        identification = {"MatchingRadius": R, "Segments": segments, "Arcs": arcs, "Features": [feature]}
        library = {"MatchingRadius": R, "Models": [{"Name": "gap", "Topology": "SameConductorGap", "Separation": 0.5 * R},
                                                   {"Name": "gap-concave-0.1", "Topology": "CurvedSameConductorGap", "Separation": 0.5 * R, "Kappa": 0.1, "Convexity": "Concave"}]}
        # Patch at 20 deg: origin on the centreline radius 3.25 R, U radially outward (e1 inner).
        angle = math.radians(20.0)
        origin = np.asarray([3.25 * R * math.cos(angle), 3.25 * R * math.sin(angle), 0.0])
        _, outward = frame((0.0, 0.0, 0.0), 20.0)
        _, inward = frame((0.0, 0.0, 0.0), 200.0)
        gates, _ = PC.placement_gates(identification, [patch(3, "gap@convex-kappa0.333333333-cubic", origin, outward, segment=0)], library, R)
        pairs = next(g for g in gates if g["Gate"] == "A10-placement-pairs")
        self.assertEqual(pairs["Status"], "PASS")
        self.assertEqual(pairs["Detail"]["ConvexityChecks"], 1)
        self.assertEqual(pairs["Detail"]["ConvexityNotEvaluable"], 0)
        for wrong in ("gap@concave-kappa0.333333333-cubic", "gap-concave-0.1"):
            gates, _ = PC.placement_gates(identification, [patch(3, wrong, origin, outward, segment=0)], library, R)
            pairs = next(g for g in gates if g["Gate"] == "A10-placement-pairs")
            self.assertEqual(pairs["Status"], "FAIL")
            self.assertEqual(pairs["Detail"]["Examples"][0]["Defects"][0]["Point"], "patch 0 convexity")
        # An inverted frame (U inward: e1 on the outer arc) fails the side check itself.
        gates, _ = PC.placement_gates(identification, [patch(3, "gap@convex-kappa0.333333333-cubic", origin, inward, segment=0)], library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-pairs")["Status"], "FAIL")
        # A runtime model whose anchor is not in the library is listed, never passed.
        gates, _ = PC.placement_gates(identification, [patch(3, "other@convex-kappa0.3-cubic", origin, outward, segment=0)], library, R)
        self.assertEqual(next(g for g in gates if g["Gate"] == "A10-placement-evaluable")["Status"], "FAIL")
        # e1 on no fitted arc (a chorded bend without an arc: the straight-like first-order
        # patches of a polyline meander): the side checks pass, the convexity is counted not
        # evaluable and never checked (review accuracy-fixes m6).
        straight_segments = [{"Key": [[0.0, y, 0.0], [10.0, y, 0.0]], "Length": 10.0, "Chain": k} for k, y in enumerate((0.0, 0.5 * R))]
        straight_feature = {"Id": 4, "Type": "SameConductorGap", "Signature": {"SeparationOverR": 0.5}, "Hash": "n", "Chirality": 1, "Length": 20.0,
                            "Portions": [[0, 0.0, 10.0], [1, 0.0, 10.0]], "Sides": [0, 1], "Frame": feature["Frame"],
                            "Match": {"Status": "Matched", "Model": "gap@convex-kappa0.05-cubic"}}
        straight_identification = {"MatchingRadius": R, "Segments": straight_segments, "Arcs": [], "Features": [straight_feature]}
        _, up = frame((0.0, 0.0, 0.0), 90.0)
        gates, _ = PC.placement_gates(straight_identification, [patch(4, "gap@convex-kappa0.05-cubic", np.asarray([5.0, 0.25 * R, 0.0]), up, segment=0)], library, R)
        pairs = next(g for g in gates if g["Gate"] == "A10-placement-pairs")
        self.assertEqual(pairs["Status"], "PASS")
        self.assertEqual(pairs["Detail"]["ConvexityChecks"], 0)
        self.assertEqual(pairs["Detail"]["ConvexityNotEvaluable"], 1)


class MirrorFormedPlacement(unittest.TestCase):
    """v3-f E-4 (decision 616 MINOR-8): the A10 clusters / context gates read a mirror-formed
    placement's REAL half (the operator's MirrorFormedPlacement EdgeWeights = RealEdges, cross-
    checked against the ENTRY stamp) and label the image half by name; on the O3 record the
    reading of record (FAIL by construction on the image edge 2.53 R / the image context end
    2.30 R beyond the plane) becomes PASS with the image edges / ends listed. A record without
    the planes or the weights is listed, never silently passed; a misplaced REAL edge still
    fails."""

    def setUp(self):
        import json
        with open(MIRROR_FORMED_FIXTURE) as source:
            self.fixture = json.load(source)
        self.identification = self.fixture["Identification"]
        self.library = self.fixture["Library"]
        self.radius = float(self.fixture["MatchingRadius"])
        self.rows = []
        for row in self.fixture["PatchRows"]:
            entry = dict(row)
            for key in ("Patch", "Feature", "ModelIndex", "Segment"):
                entry[key] = int(row[key])
            for key in ("Weight", "ModelWeight", "QuadratureWeight", "SideFactor", "CouponDepth", "S0", "S1"):
                entry[key] = float(row[key])
            entry["Origin"] = np.asarray([float(row[k]) for k in ("OriginX", "OriginY", "OriginZ")])
            entry["Axes"] = np.asarray([[float(row[f"Axis{a}{d}"]) for d in "XYZ"] for a in "UVW"])
            self.rows.append(entry)

    def gates(self, identification=None, library=None):
        gates, summary = PC.placement_gates(identification or self.identification, self.rows, library or self.library, self.radius)
        return {g["Gate"]: g for g in gates}, summary

    def test_fixture_is_the_o3_record_and_its_reading_of_record(self):
        feature = self.identification["Features"][0]
        self.assertEqual((feature["Id"], feature["Type"], feature["Hash"][:12]), (90, "SpatialEdgeCluster", "d7c875318447"))
        self.assertEqual(feature["Match"]["Status"], "Matched")
        self.assertTrue(feature["Match"]["Note"].startswith("mirror-formed cluster configuration placed on its REAL half"))
        self.assertEqual(feature["Mirror"]["Contract"]["RealPortions"], [0, 3])
        self.assertEqual(feature["Mirror"]["Contract"]["Planes"], [0])
        model = self.library["Models"][0]
        self.assertEqual(model["Name"], feature["Match"]["Model"])
        self.assertEqual([edge["Weight"] for edge in model["Edges"]], [1.0] * 15 + [0.0, 0.0, 1.0])
        placement = self.identification["Diagnostics"]["MirrorFormedPlacement"]["Patches"][0]
        self.assertEqual((placement["Feature"], placement["Patch"], placement["RealEdges"], placement["ImageEdges"]), (90, 7347, 16, 2))
        self.assertEqual(placement["EdgeWeights"], [edge["Weight"] for edge in model["Edges"]])
        self.assertEqual(self.rows[0]["Patch"], 7347)
        plane = self.identification["Diagnostics"]["MirrorBand"]["Planes"][0]
        self.assertEqual(plane["Normal"], [-1.0, 0.0, 0.0])
        self.assertAlmostEqual(plane["Offset"], 830.0, places=5)
        # The not-mirror-aware reading of record: the image edges 15 / 16 and the image context
        # ends 1 / 2 / 3-begin, all at world x < -830 (beyond the plane).
        of_record = self.fixture["PlacementOfRecord"]
        self.assertEqual((of_record["Clusters"]["Status"], of_record["Context"]["Status"]), ("FAIL", "FAIL"))
        self.assertAlmostEqual(of_record["Clusters"]["WorstDeviationOverR"], 2.5281348249990288)
        self.assertAlmostEqual(of_record["Context"]["WorstDeviationOverR"], 2.3038284115964296)
        self.assertEqual([d["Point"] for d in of_record["Clusters"]["Defects"]], ["edge 15 begin", "edge 15 end", "edge 16 begin", "edge 16 end"])
        self.assertEqual([d["Point"] for d in of_record["Context"]["Defects"]],
                         ["context 1 begin", "context 1 end", "context 2 begin", "context 2 end", "context 3 begin"])
        inside = PC.mirror_plane_inside(self.identification, [0])
        for defect in of_record["Clusters"]["Defects"] + of_record["Context"]["Defects"]:
            self.assertLess(inside(np.asarray(defect["Mapped"])), -PC.SIGNATURE_TOLERANCE_OVER_R * self.radius)

    def test_real_half_checked_image_half_labelled(self):
        by_name, summary = self.gates()
        clusters = by_name["A10-placement-clusters"]
        self.assertEqual(clusters["Status"], "PASS")
        self.assertEqual(clusters["Detail"]["Features"], 1)
        self.assertEqual((clusters["Detail"]["MirrorFormedFeatures"], clusters["Detail"]["ImageEdges"]), (1, 2))
        self.assertLess(clusters["Detail"]["WorstDeviationOverR"], PC.SIGNATURE_TOLERANCE_OVER_R)
        # 16 real edges x 2 ends + the 4 claimed points (the straight portion's 2 ends + the arc's 2 ends)
        self.assertEqual(clusters["Detail"]["Checks"], 32 + 4)
        context = by_name["A10-context"]
        self.assertEqual(context["Status"], "PASS")
        self.assertEqual((context["Detail"]["MirrorFormedFeatures"], context["Detail"]["ImageEnds"], context["Detail"]["Checks"]), (1, 5, 5))
        self.assertLess(context["Detail"]["WorstDeviationOverR"], PC.SIGNATURE_TOLERANCE_OVER_R)
        self.assertEqual(summary["NotEvaluable"], [])
        self.assertNotIn("A10-placement-evaluable", by_name)
        # The split by name: the placement record's RealEdges (15 arc chords + edge 17) and the
        # image edges 15 / 16, whose mapped ends are the defect points of record.
        split, reason = PC.mirror_formed_placement(self.identification, self.identification["Features"][0], self.library["Models"][0], self.rows[0])
        self.assertIsNone(reason)
        self.assertEqual((split["Record"], split["RealEdges"], split["ImageEdges"], split["RealPortions"], split["Planes"]),
                         ("MirrorFormedPlacement", list(range(15)) + [17], [15, 16], [0, 3], [0]))
        pairs, source = PC.model_cluster_edges(self.library["Models"][0], float(self.library["MatchingRadius"]))
        self.assertEqual(source, "Edges")
        scale = self.radius / float(self.library["MatchingRadius"])
        origin, axes = self.rows[0]["Origin"], self.rows[0]["Axes"]
        of_record = {d["Point"]: d["Mapped"] for d in self.fixture["PlacementOfRecord"]["Clusters"]["Defects"]}
        for k in split["ImageEdges"]:
            np.testing.assert_allclose(origin + (scale * pairs[k][0]) @ axes, of_record[f"edge {k} begin"], atol=1e-9)
            np.testing.assert_allclose(origin + (scale * pairs[k][1]) @ axes, of_record[f"edge {k} end"], atol=1e-9)
        _, context_entries = PC.context_gate(self.identification, self.rows, self.library, self.radius)
        self.assertEqual(context_entries[0]["ImageEnds"], ["context 1 begin", "context 1 end", "context 2 begin", "context 2 end", "context 3 begin"])
        self.assertEqual(context_entries[0]["Defects"], [])
        self.assertEqual(context_entries[0]["Checks"], 5)

    def test_without_the_mirror_record_the_reading_of_record_returns(self):
        """The same inputs read as an ORDINARY placement (the feature's Mirror record removed):
        the not-mirror-aware FAILs of record, to the digit — the fix changes nothing else."""
        import copy
        plain = copy.deepcopy(self.identification)
        del plain["Features"][0]["Mirror"]
        by_name, _ = self.gates(identification=plain)
        clusters, context = by_name["A10-placement-clusters"], by_name["A10-context"]
        self.assertEqual((clusters["Status"], context["Status"]), ("FAIL", "FAIL"))
        self.assertAlmostEqual(clusters["Detail"]["WorstDeviationOverR"], self.fixture["PlacementOfRecord"]["Clusters"]["WorstDeviationOverR"], places=12)
        self.assertAlmostEqual(context["Detail"]["WorstDeviationOverR"], self.fixture["PlacementOfRecord"]["Context"]["WorstDeviationOverR"], places=12)
        self.assertEqual([d["Point"] for d in clusters["Detail"]["Examples"][0]["Defects"]],
                         [d["Point"] for d in self.fixture["PlacementOfRecord"]["Clusters"]["Defects"]])
        self.assertEqual([d["Point"] for d in context["Detail"]["Examples"][0]["Defects"]],
                         [d["Point"] for d in self.fixture["PlacementOfRecord"]["Context"]["Defects"]])
        self.assertEqual((clusters["Detail"]["MirrorFormedFeatures"], context["Detail"]["MirrorFormedFeatures"]), (0, 0))

    def test_misplaced_real_edge_still_fails_and_incomplete_records_are_listed(self):
        import copy
        # A real edge moved by 0.01 R along its gap direction (library units): the real half is checked.
        library = copy.deepcopy(self.library)
        edge = library["Models"][0]["Edges"][17]
        self.assertEqual(edge["Weight"], 1.0)
        edge["Point"] = [p + 0.01 * float(self.library["MatchingRadius"]) * g for p, g in zip(edge["Point"], edge["GapDirection"])]
        by_name, _ = self.gates(library=library)
        self.assertEqual(by_name["A10-placement-clusters"]["Status"], "FAIL")
        self.assertGreater(by_name["A10-placement-clusters"]["Detail"]["WorstDeviationOverR"], 0.009)
        # An image edge moved the same way: labelled, so no defect (the image half is not checked).
        library = copy.deepcopy(self.library)
        edge = library["Models"][0]["Edges"][15]
        self.assertEqual(edge["Weight"], 0.0)
        edge["Point"] = [p + 0.01 * float(self.library["MatchingRadius"]) * g for p, g in zip(edge["Point"], edge["GapDirection"])]
        by_name, _ = self.gates(library=library)
        self.assertEqual(by_name["A10-placement-clusters"]["Status"], "PASS")
        # The placement record's weights disagreeing with the ENTRY stamp: not evaluable, listed.
        identification = copy.deepcopy(self.identification)
        identification["Diagnostics"]["MirrorFormedPlacement"]["Patches"][0]["EdgeWeights"][15] = 1.0
        by_name, summary = self.gates(identification=identification)
        self.assertEqual(by_name["A10-placement-clusters"]["Detail"]["Features"], 0)
        self.assertEqual(by_name["A10-placement-evaluable"]["Status"], "FAIL")
        self.assertIn("differ from the model's ENTRY stamp", summary["NotEvaluable"][0]["Reason"])
        self.assertIn("differ from the model's ENTRY stamp", by_name["A10-context"]["Detail"]["Examples"][0]["Defects"][0]["Reason"])
        # No placement record (an older manifest): the ENTRY stamp alone carries the split.
        identification = copy.deepcopy(self.identification)
        del identification["Diagnostics"]["MirrorFormedPlacement"]
        by_name, summary = self.gates(identification=identification)
        self.assertEqual((by_name["A10-placement-clusters"]["Status"], by_name["A10-context"]["Status"]), ("PASS", "PASS"))
        split, reason = PC.mirror_formed_placement(identification, identification["Features"][0], self.library["Models"][0], self.rows[0])
        self.assertEqual((reason, split["Record"], split["ImageEdges"]), (None, "EntryStamp", [15, 16]))
        # Neither record nor stamp: listed, never silently passed.
        library = copy.deepcopy(self.library)
        for edge in library["Models"][0]["Edges"]:
            del edge["Weight"]
        by_name, summary = self.gates(identification=identification, library=library)
        self.assertEqual(by_name["A10-placement-evaluable"]["Status"], "FAIL")
        self.assertIn("without a MirrorFormedPlacement record or an ENTRY stamp", summary["NotEvaluable"][0]["Reason"])
        self.assertEqual(by_name["A10-context"]["Status"], "FAIL")
        # The planes missing from the record: the context gate fails by name (the image side cannot be told).
        identification = copy.deepcopy(self.identification)
        del identification["Diagnostics"]["MirrorBand"]
        by_name, _ = self.gates(identification=identification)
        self.assertEqual(by_name["A10-placement-clusters"]["Status"], "PASS")
        self.assertEqual(by_name["A10-context"]["Status"], "FAIL")
        self.assertEqual(by_name["A10-context"]["Detail"]["Examples"][0]["Defects"][0]["Point"], "mirror planes")
        self.assertIsNone(PC.mirror_plane_inside(identification, [0]))



if __name__ == "__main__":
    unittest.main()
