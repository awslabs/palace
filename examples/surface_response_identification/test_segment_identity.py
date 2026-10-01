# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""segment_identity (plan (b)-4, E1) on hand-built version-2 manifests of one layout under two segment
splits: identical readings pass whatever the split; a changed Type, a lost feature, an exclusion facing an
assignment and a changed cluster edge count are class D; a stack whose mean separation moved is class B; a
curved edge whose hash moved is class C; the comparison region excludes what lies outside it."""
import json
import math
import os
import sys
import tempfile
import unittest

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import segment_identity as SI  # noqa: E402
from . import window_query as WQ  # noqa: E402

R = 2.0
LAW = '{"Type":"PEC"}'


def manifest(pieces, features):
    """pieces: list of (p0, p1, feature id or ("Exclusion", class)); features: id -> record fields."""
    segments = []
    for p0, p1, owner in pieces:
        length = math.dist(p0, p1)
        seg = {"Key": [[p0[0], p0[1], 0.0], [p1[0], p1[1], 0.0]], "Length": length, "Chain": 0}
        if isinstance(owner, tuple):
            seg["Exclusion"] = {"Class": owner[1], "Reason": "test"}
        elif isinstance(owner, list):
            seg["Portions"] = [[s0, s1, fid] for s0, s1, fid in owner]
        else:
            seg["Portions"] = [[0.0, length, owner]]
        segments.append(seg)
    records = []
    for fid, spec in sorted(features.items()):
        record = {"Id": fid, "Type": spec["Type"], "Hash": spec.get("Hash", spec["Type"] + "-hash"), "Signature": dict(spec.get("Signature", {}), Type=spec["Type"]),
                  "ExactParameters": spec.get("Exact", True), "Length": 1.0, "Portions": [], "Vertices": [], "Chirality": 1,
                  "Frame": {"Origin": [0.0, 0.0, 0.0], "Axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}}
        records.append(record)
    return {"Version": 2, "Identification": {"Version": 2, "MatchingRadius": R, "Segments": segments, "Features": records,
                                             "Exclusions": [], "Vertices": [], "Totals": {}, "GeometryDigest": "x"}}


ISO = {"Type": "IsolatedEdge", "Hash": "iso"}
CORNER = {"Type": "ConvexCorner", "Hash": "cvx", "Signature": {"AngleDegrees": 90.0, "CornerRadiusOverR": 0.0}}
REGION = (-1.0, 101.0, -1.0, 101.0)


class SegmentIdentity(unittest.TestCase):
    def test_identical_layout_under_two_splits_passes(self):
        # an edge from (0,0) to (100,0): corner claims R at both ends, the isolated edge in between
        a = manifest([((0.0, 0.0), (100.0, 0.0), [(0.0, R, 1), (R, 100.0 - R, 0), (100.0 - R, 100.0, 2)])],
                     {0: ISO, 1: CORNER, 2: CORNER})
        b = manifest([((0.0, 0.0), (37.0, 0.0), [(0.0, R, 1), (R, 37.0, 0)]), ((37.0, 0.0), (100.0, 0.0), [(0.0, 63.0 - R, 0), (63.0 - R, 63.0, 2)])],
                     {0: ISO, 1: CORNER, 2: CORNER})
        r = SI.run(a, b, REGION)
        self.assertEqual(r["Verdict"], "PASS")
        self.assertEqual(set(r["Forward"]), {"Identical"})
        self.assertEqual(set(r["Reverse"]), {"Identical"})
        self.assertAlmostEqual(r["Forward"]["Identical"]["Length"], 100.0, places=9)

    def test_type_change_is_class_D(self):
        a = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: ISO})
        b = manifest([((0.0, 0.0), (100.0, 0.0), [(0.0, 60.0, 0), (60.0, 100.0, 1)])],
                     {0: ISO, 1: {"Type": "CurvedEdge", "Hash": "crv", "Signature": {"RadiusOverR": 3.0}}})
        r = SI.run(a, b, REGION)
        self.assertEqual(r["Verdict"], "FAIL")
        self.assertAlmostEqual(r["Forward"]["D"]["Length"], 40.0, places=9)
        self.assertEqual(r["Forward"]["D"]["Examples"][0]["Why"], "type changed")

    def test_straight_class_hash_change_is_class_A(self):
        a = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: ISO})
        b = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: {"Type": "IsolatedEdge", "Hash": "iso-other"}})
        r = SI.run(a, b, REGION)
        self.assertEqual(r["Verdict"], "FAIL")
        self.assertAlmostEqual(r["Forward"]["A"]["Length"], 100.0, places=9)

    def test_lost_perimeter_and_exclusion_are_class_D(self):
        a = manifest([((0.0, 50.0), (100.0, 50.0), 0), ((0.0, 0.0), (100.0, 0.0), 0)], {0: ISO})
        b = manifest([((0.0, 0.0), (100.0, 0.0), ("Exclusion", "Port"))], {0: ISO})
        r = SI.run(a, b, REGION, examples=1000)
        self.assertEqual(r["Verdict"], "FAIL")
        whys = {e["Why"] for e in r["Forward"]["D"]["Examples"]}
        self.assertIn("exclusion vs assignment or a different exclusion class", whys)
        self.assertIn("no perimeter within tolerance on the other side", whys)
        self.assertAlmostEqual(r["Forward"]["D"]["Length"], 200.0, places=9)

    def test_stack_parameter_change_is_class_B_with_deltas(self):
        stack = {"Type": "ParallelEdgeCluster", "Hash": "st1", "Signature": {"Edges": [{"Conductor": 1, "OffsetOverR": 0.0}, {"Conductor": 1, "OffsetOverR": 1.0}]}}
        moved = {"Type": "ParallelEdgeCluster", "Hash": "st2", "Signature": {"Edges": [{"Conductor": 1, "OffsetOverR": 0.0}, {"Conductor": 1, "OffsetOverR": 1.01}]}}
        a = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: stack})
        b = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: moved})
        r = SI.run(a, b, REGION)
        self.assertEqual(r["Verdict"], "PASS")
        self.assertAlmostEqual(r["Forward"]["B"]["Length"], 100.0, places=9)
        self.assertEqual(r["Forward"]["B"]["Examples"][0]["Deltas"], {"Edges.[1].OffsetOverR": [1.0, 1.01]})
        # a changed conductor topology of the stack is not a parameter change
        other = {"Type": "ParallelEdgeCluster", "Hash": "st3", "Signature": {"Edges": [{"Conductor": 1, "OffsetOverR": 0.0}, {"Conductor": 2, "OffsetOverR": 1.0}]}}
        r = SI.run(a, manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: other}), REGION)
        self.assertEqual(r["Verdict"], "FAIL")
        self.assertIn("D", r["Forward"])

    def test_curved_edge_hash_change_and_displaced_chord_are_class_C(self):
        curved = {"Type": "CurvedEdge", "Hash": "c1", "Signature": {"RadiusOverR": 3.0}}
        curved2 = {"Type": "CurvedEdge", "Hash": "c2", "Signature": {"RadiusOverR": 3.001}}
        a = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: curved})
        b = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: curved2})
        r = SI.run(a, b, REGION)
        self.assertEqual(r["Verdict"], "PASS")
        self.assertAlmostEqual(r["Forward"]["C"]["Length"], 100.0, places=9)
        # the same curve re-meshed 0.3 R off the chord: found within the curve tolerance -> C
        b_off = manifest([((0.0, 0.3 * R), (100.0, 0.3 * R), 0)], {0: curved})
        r = SI.run(a, b_off, REGION)
        self.assertEqual(r["Verdict"], "PASS")
        self.assertIn("chord displaced", r["Forward"]["C"]["Examples"][0]["Why"])
        # a straight edge 0.3 R off is lost perimeter (D)
        a_iso = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: ISO})
        b_iso = manifest([((0.0, 0.3 * R), (100.0, 0.3 * R), 0)], {0: ISO})
        self.assertEqual(SI.run(a_iso, b_iso, REGION)["Verdict"], "FAIL")

    def test_cluster_edge_count_change_is_class_D(self):
        c3 = {"Type": "SpatialEdgeCluster", "Hash": "k3", "Signature": {"EdgeCount": 3, "Portions": [{"P": [0, 0, 1, 1]}]}}
        c4 = {"Type": "SpatialEdgeCluster", "Hash": "k4", "Signature": {"EdgeCount": 4, "Portions": [{"P": [0, 0, 1, 1]}]}}
        a = manifest([((0.0, 0.0), (10.0, 0.0), 0)], {0: c3})
        b = manifest([((0.0, 0.0), (10.0, 0.0), 0)], {0: c4})
        r = SI.run(a, b, REGION)
        self.assertEqual(r["Verdict"], "FAIL")
        self.assertEqual(r["Forward"]["D"]["Examples"][0]["Why"], "cluster edge count changed")

    def test_region_and_offset(self):
        a = manifest([((0.0, 0.0), (100.0, 0.0), 0)], {0: ISO})
        b = manifest([((10.0, 5.0), (110.0, 5.0), 0)], {0: {"Type": "IsolatedEdge", "Hash": "iso-other"}})
        # without the offset nothing coincides; with offset (-10, -5) the pieces coincide and differ in hash (A);
        # a region covering x in [0, 50] sees half the length
        r = SI.run(a, b, (0.0, 50.0, -1.0, 1.0), offset=(-10.0, -5.0, 0.0))
        self.assertAlmostEqual(r["Forward"]["A"]["Length"], 50.0, places=6)
        self.assertNotIn("Identical", r["Forward"])
        r = SI.run(a, b, (0.0, 50.0, -1.0, 1.0))
        self.assertEqual(r["Forward"]["D"]["Examples"][0]["Why"], "no perimeter within tolerance on the other side")

    def test_topology_key_strips_parameters_only(self):
        sig = {"Type": "SameConductorStrip", "SeparationOverR": 1.05, "Edges": [{"Conductor": 1, "GapSide": -1, "OffsetOverR": 0.0}]}
        self.assertEqual(json.loads(SI.topology_key(sig)), {"Type": "SameConductorStrip", "Edges": [{"Conductor": 1, "GapSide": -1}]})
        self.assertEqual(SI.parameters(sig), {"SeparationOverR": 1.05, "Edges.[0].OffsetOverR": 0.0})

    def cut_arc_layout(self):
        """The chip (A): a chain lead -> 3-joint 90-deg bend (r 14, two 10.7-um chords) -> trail, one IsolatedEdge (the bend
        inside its chain); the window (B) cut 2.2 um beyond the start joint: the two joints that remain read as 135 / 157.5-deg
        corners (the collapse the arc rule predicts); region = x >= -8 (the wall at x = -11.8 plus the margin)."""
        j0, j1, j2 = (-14.0, 0.0), (-14.0 * math.cos(math.radians(45)), 14.0 * math.sin(math.radians(45))), (0.0, 14.0)
        lead, trail = (-14.0, -30.0), (30.0, 14.0)
        chord = math.dist(j0, j1)
        a = manifest([(lead, j0, 0), (j0, j1, 0), (j1, j2, 0), (j2, trail, 0)], {0: ISO})
        for i in (1, 2):
            a["Identification"]["Segments"][i]["Arc"] = 0
        a["Identification"]["Arcs"] = [{"Center": [0.0, 0.0, 0.0], "Joints": 3, "Kind": "Bend", "Radius": 14.0, "RadiusOverR": 7.0,
                                        "Segments": 2, "TurnDegrees": 90.0, "MaxChordSagittaOverR": 0.5}]
        a["Identification"]["Conventions"] = {"JointNoiseSagittaOverR": 0.05}
        corner1 = {"Type": "ConvexCorner", "Hash": "c135", "Signature": {"AngleDegrees": 135.0, "CornerRadiusOverR": 0.0}}
        corner2 = {"Type": "ConvexCorner", "Hash": "c157", "Signature": {"AngleDegrees": 157.5, "CornerRadiusOverR": 0.0}}
        b = manifest([(j0, j1, [(0.0, chord - R, 0), (chord - R, chord, 1)]),
                      (j1, j2, [(0.0, R, 1), (R, chord - R, 0), (chord - R, chord, 2)]),
                      (j2, trail, [(0.0, R, 2), (R, 30.0, 0)])], {0: ISO, 1: corner1, 2: corner2})
        box = (-11.8, 100.0, -100.0, 100.0)
        exclusions = WQ.cut_arc_exclusions(a["Identification"], lambda i: 5, box, None, None)
        return a, b, (-8.0, 100.0, -100.0, 100.0), exclusions

    def test_cut_arc_exclusion_is_reported_never_a_defect(self):
        # VALIDATION-PLAN (h)-9: without the exclusions the two corners are class D (2 x 2 R per direction); with them the
        # cut arc's chords + 3 R along the arms are class Excluded, the per-arc record keeps the would-be classes, PASS
        a, b, region, exclusions = self.cut_arc_layout()
        self.assertEqual(exclusions["Arcs"][0]["Prediction"], "Collapse")
        plain = SI.run(a, b, region)
        self.assertEqual(plain["Verdict"], "FAIL")
        self.assertAlmostEqual(plain["Forward"]["D"]["Length"], 2 * R, delta=0.1)  # the 157.5-deg corner (the 135-deg one is at x < -8)
        self.assertNotIn("CutArcs", plain)
        r = SI.run(a, b, region, exclusions=exclusions)
        self.assertEqual(r["Verdict"], "PASS")
        self.assertNotIn("D", r["Forward"])
        self.assertNotIn("D", r["Reverse"])
        j1x = a["Identification"]["Segments"][2]["Key"][0][0]
        chord2 = a["Identification"]["Segments"][2]["Length"]
        inside_region = chord2 * (0.0 - (-8.0)) / (0.0 - j1x)  # the first chord and the start of the second lie at x < -8
        self.assertAlmostEqual(r["Forward"]["Excluded"]["Length"], inside_region + 3.0 * R, delta=0.3)  # + 3 R along the trail arm
        self.assertAlmostEqual(r["Reverse"]["Excluded"]["Length"], r["Forward"]["Excluded"]["Length"], delta=0.3)
        self.assertEqual(r["Summary"]["A->B"]["Excluded"]["Samples"], r["Forward"]["Excluded"]["Samples"])
        arcs = r["CutArcs"]["Arcs"]
        self.assertEqual([(x["Arc"], x["Prediction"], x["MismatchLength"]) for x in arcs], [(0, "Collapse", 0.0)])
        self.assertAlmostEqual(arcs[0]["WouldBe"]["D"], 2 * 2 * R, delta=0.2)  # the corner, both directions
        self.assertEqual(r["CutArcs"]["PredictionMismatches"], [])
        self.assertAlmostEqual(r["CutArcs"]["ExcludedLength"]["A->B"], r["Forward"]["Excluded"]["Length"], places=12)
        self.assertIn("Cut-arc exclusions", SI.markdown(r))
        self.assertIn("| Excluded |", SI.markdown(r).replace("| A->B | Excluded |", "| Excluded |"))

    def test_prediction_mismatch_is_flagged_and_stays_a_defect(self):
        # a cut arc the rule predicts unchanged ("None", no pieces) that reads differently: the D samples on its segments are
        # reported as PredictionMismatch and the verdict stays FAIL
        a, b, region, exclusions = self.cut_arc_layout()
        arc = exclusions["Arcs"][0]
        arc.update({"Prediction": "None", "Pieces": [], "ExcludedLength": 0.0})
        r = SI.run(a, b, region, exclusions=exclusions)
        self.assertEqual(r["Verdict"], "FAIL")
        self.assertNotIn("Excluded", r["Forward"])
        self.assertEqual(r["CutArcs"]["PredictionMismatches"], [0])
        self.assertGreater(r["CutArcs"]["PredictionMismatchLength"], 0.0)
        self.assertTrue(any("PredictionMismatch" in e["Why"] for e in r["Forward"]["D"]["Examples"]))

    def test_cli_exclude_reads_the_window_record(self):
        a, b, region, exclusions = self.cut_arc_layout()
        with tempfile.TemporaryDirectory() as tmp:
            pa, pb, pe, out = (os.path.join(tmp, n) for n in ("a.json", "b.json", "e0.json", "out.json"))
            json.dump(a, open(pa, "w"))
            json.dump(b, open(pb, "w"))
            json.dump({"Windows": {"W": {"E1CutArcExclusions": exclusions}}}, open(pe, "w"))
            rc = SI.main(["--a", pa, "--b", pb, "--region", *[str(v) for v in region], "--exclude", pe, "W", "--output", out])
            self.assertEqual(rc, 0)
            result = json.load(open(out))
            self.assertEqual((result["Verdict"], result["Exclude"]), ("PASS", [pe, "W"]))
            self.assertEqual(result["CutArcs"]["PredictionMismatches"], [])


if __name__ == "__main__":
    unittest.main()
