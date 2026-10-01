# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""window_query (plan stage 0, E0) on a hand-built version-2 manifest: a ground with a square hole (the
CPW-gap analogue) holding a square island, plus a 2-um strip pair. Windows that hold the island whole /
cut it / miss it give the expected perimeter, features, bodies, terminal proposal and mean-separation
prediction."""
import json
import math
import os
import sys
import tempfile
import unittest

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import window_query as WQ  # noqa: E402

R = 2.0
LAW = '{"Type":"PEC"}'


def square(x0, y0, x1, y1, z=0.0):
    """Four segments of an axis-aligned square (keys lexicographically ordered as the manifest does)."""
    corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
    segments = []
    for i in range(4):
        a, b = corners[i], corners[(i + 1) % 4]
        if (b[0], b[1]) < (a[0], a[1]):
            a, b = b, a
        segments.append({"Key": [[a[0], a[1], z], [b[0], b[1], z]], "Length": math.dist(a, b), "Chain": i})
    return segments, corners


class Builder:
    def __init__(self):
        self.segments, self.features = [], []

    def feature(self, ftype, signature, origin, length, portions, **extra):
        fid = len(self.features)
        record = {"Id": fid, "Type": ftype, "Signature": dict(signature, Type=ftype), "Hash": f"{ftype}-{fid:04d}" + "0" * 40,
                  "Chirality": 1, "ExactParameters": True, "Length": length, "Portions": portions, "Vertices": [],
                  "Frame": {"Origin": [origin[0], origin[1], origin[2] if len(origin) > 2 else 0.0], "Axes": [[1, 0, 0], [0, 1, 0], [0, 0, 1]]}}
        record.update(extra)
        self.features.append(record)
        return fid

    def add_loop(self, segments, corners, corner_type, exclusion=None):
        """A square loop: whole segments as one isolated edge each (or excluded), corners as vertex features
        claiming R on both arms."""
        base = len(self.segments)
        for s in segments:
            s = dict(s)
            if exclusion:
                s["Exclusion"] = {"Class": exclusion, "Reason": "test"}
            else:
                s["Portions"] = []
            self.segments.append(s)
        if exclusion:
            return
        for i, c in enumerate(corners):
            # arms: the two segments touching corner c
            arms = [base + k for k in range(4) if any(tuple(p[:2]) == c for p in segments[k]["Key"])]
            portions = []
            for si in arms:
                key = self.segments[si]["Key"]
                if tuple(key[0][:2]) == c:
                    portions.append([si, 0.0, R])
                else:
                    portions.append([si, self.segments[si]["Length"] - R, self.segments[si]["Length"]])
            fid = self.feature(corner_type, {"AngleDegrees": 90.0, "CornerRadiusOverR": 0.0, "Interfaces": ["SA"], "Law": LAW},
                               (c[0], c[1], segments[0]["Key"][0][2]), 2 * R, portions)
            for si, s0, s1 in portions:
                self.segments[si]["Portions"].append([s0, s1, fid])
        for k in range(4):
            si = base + k
            length = self.segments[si]["Length"]
            fid = self.feature("IsolatedEdge", {"Interfaces": ["SA"], "Law": LAW}, self.segments[si]["Key"][0][:2], length - 2 * R,
                               [[si, R, length - R]])
            self.segments[si]["Portions"].append([R, length - R, fid])
            self.segments[si]["Portions"].sort()

    def manifest(self):
        return {"Version": 2, "Identification": {"Version": 2, "MatchingRadius": R, "Segments": self.segments, "Features": self.features,
                                                 "Vertices": [], "Exclusions": [{"Class": "SimulationBoundary", "Reason": "t", "Count": 4, "Length": 800.0}],
                                                 "Totals": {}, "GeometryDigest": "x", "Conventions": {}}}


def build():
    b = Builder()
    outer_segments, outer_corners = square(-100.0, -100.0, 100.0, 100.0)
    b.add_loop(outer_segments, outer_corners, None, exclusion="SimulationBoundary")  # ground outer boundary
    hole_segments, hole_corners = square(-20.0, -20.0, 20.0, 20.0)
    b.add_loop(hole_segments, hole_corners, "ConcaveCorner")  # hole in the ground: concave metal corners
    island_segments, island_corners = square(-10.0, -10.0, 10.0, 10.0)
    b.add_loop(island_segments, island_corners, "ConvexCorner")  # island: convex metal corners
    # a 2-um strip pair far away: two parallel isolated edges at y = 60 and 62, x in [40, 80]
    for y in (60.0, 62.0):
        b.segments.append({"Key": [[40.0, y, 0.0], [80.0, y, 0.0]], "Length": 40.0, "Chain": 99, "Portions": []})
    s0, s1 = len(b.segments) - 2, len(b.segments) - 1
    fid = b.feature("SameConductorStrip", {"SeparationOverR": 1.0, "Interfaces": ["SA"], "Law": LAW,
                                           "Edges": [{"OffsetOverR": 0.0}, {"OffsetOverR": 1.0}]}, (40.0, 60.0), 40.0,
                    [[s0, 0.0, 40.0], [s1, 0.0, 40.0]], Sides=[0, 1])
    b.segments[s0]["Portions"].append([0.0, 40.0, fid])
    b.segments[s1]["Portions"].append([0.0, 40.0, fid])
    return b.manifest()


class WindowInventory(unittest.TestCase):
    def setUp(self):
        self.manifest = build()

    def test_loops_bodies_and_metal_side(self):
        ident = self.manifest["Identification"]
        loops = WQ.PerimeterLoops(ident["Segments"], ident["Features"], R)
        closed = [l for l in loops.loops.values() if l["Closed"]]
        self.assertEqual(len(closed), 3)
        depths = sorted(l["Depth"] for l in closed)
        self.assertEqual(depths, [0, 1, 2])
        # the hole (concave metal corners) belongs to the plane's ground: the excluded outer loop carries no
        # corner feature, so the ground is the synthetic body -1; the island (convex corners) is its own body
        by_depth = {l["Depth"]: l for l in closed}
        self.assertEqual([l["MetalInside"] for l in (by_depth[0], by_depth[1], by_depth[2])], [False, False, True])
        self.assertEqual(by_depth[1]["Body"], -1)
        self.assertEqual(by_depth[2]["Body"], [k for k, l in loops.loops.items() if l is by_depth[2]][0])
        self.assertEqual(loops.metal_side_checks, 8)
        self.assertEqual(loops.metal_side_disagreements, [])

    def test_window_holding_the_island(self):
        r = WQ.inventory(self.manifest, {"W": (-30.0, 30.0, -30.0, 30.0)}, None, WQ.DEFAULT_WEIGHTS)
        w = r["Windows"]["W"]
        self.assertAlmostEqual(w["Perimeter"]["Total"], 160.0 + 80.0, places=9)  # hole 4 x 40 + island 4 x 20
        self.assertAlmostEqual(w["Perimeter"]["Assigned"], 240.0, places=9)
        self.assertEqual(w["Perimeter"]["Excluded"], {})
        self.assertEqual(w["ConductorCount"], 2)
        kinds = {b["Kind"] for b in w["Bodies"]}
        self.assertEqual(kinds, {"Ground", "Island"})
        island = [b for b in w["Bodies"] if b["Kind"] == "Island"][0]
        self.assertTrue(island["EntirelyInside"])
        self.assertFalse(island["CutByWall"])
        self.assertEqual(w["ProposedTerminal"], island["Body"])
        roles = {a["Kind"]: a["Role"] for a in w["TerminalAssignment"]}
        self.assertTrue(roles["Island"].startswith("Terminal"))
        self.assertTrue(roles["Ground"].startswith("Ground"))
        self.assertEqual(w["CutFeatureCount"], 0)
        self.assertEqual(w["ByType"]["ConvexCorner"]["Features"], 4)
        self.assertEqual(w["ByType"]["ConcaveCorner"]["Features"], 4)
        self.assertEqual(w["ByType"]["IsolatedEdge"]["Features"], 8)
        # proxy: convex corners 4 x 2R x 14.6, isolated 8 edges (16 + 36) x 1, concave 4 x 2R x 0.01
        expected = 4 * 2 * R * 14.6 + (4 * 16.0 + 4 * 36.0) + 4 * 2 * R * 0.01
        self.assertAlmostEqual(sum(t["Proxy"] for t in w["ByType"].values()), expected, places=9)

    def test_window_cutting_the_island(self):
        r = WQ.inventory(self.manifest, {"W": (0.0, 30.0, -30.0, 30.0)}, 2.0, WQ.DEFAULT_WEIGHTS)
        w = r["Windows"]["W"]
        island = [b for b in w["Bodies"] if b["Kind"] == "Island"][0]
        self.assertTrue(island["CutByWall"])
        self.assertFalse(island["EntirelyInside"])
        # decision-180 rule (M2): no whole island -> the cut non-ground conductor is the terminal, open-terminated
        # 3 R + margin inside the wall
        self.assertEqual(w["ProposedTerminal"], island["Body"])
        self.assertEqual(w["Excitation"]["Kind"], "OpenTerminatedTrace")
        self.assertAlmostEqual(w["Excitation"]["OpenSetback"], 3.0 * R + 2.0, places=12)
        # the E1 region excludes 3 R around the open end: margin = OpenSetback + 3 R = 14 (the E0 region keeps the margin 2)
        self.assertAlmostEqual(w["E1Margin"], 14.0, places=12)
        self.assertEqual(w["E1ComparisonRegion"], [14.0, 16.0, -16.0, 16.0])
        self.assertEqual(w["ComparisonRegion"], [2.0, 28.0, -28.0, 28.0])
        roles = {a["Kind"]: a["Role"] for a in w["TerminalAssignment"]}
        self.assertIn("open-terminated", roles["Island"])
        # island perimeter inside: the x = 10 side (20) + halves of the y = +-10 sides (10 each) = 40
        self.assertAlmostEqual(island["PerimeterInside"], 40.0, places=9)
        # the isolated edges along y = +-10 are cut by the wall at x = 0 -> CutByWall on those features
        cut = [f for f in w["Features"] if f["CutByWall"]]
        self.assertTrue(all(f["Type"] == "IsolatedEdge" for f in cut))
        self.assertEqual(len(cut), 4)  # two island sides + two hole sides
        # the margin (2 um from the walls) accounts every isolated portion within 2 um of x = 0 / y = +-30 etc.
        self.assertGreater(sum(f["LengthInMargin"] for f in w["Features"]), 0.0)

    def test_wall_hugging_cut_island_is_not_a_realisable_terminal(self):
        # the island's portion inside (5, 35) lies within the 8-um setback band of the x = 5 wall: open-terminating it
        # removes all its metal -> no excitation (move / resize / drop), the candidate recorded as vanishing
        r = WQ.inventory(self.manifest, {"W": (5.0, 35.0, -30.0, 30.0)}, 2.0, WQ.DEFAULT_WEIGHTS)
        w = r["Windows"]["W"]
        island = [b for b in w["Bodies"] if b["Kind"] == "Island"][0]
        self.assertTrue(island["CutByWall"])
        self.assertAlmostEqual(island["PerimeterInsideSetback"], 0.0, places=12)
        self.assertEqual(w["Excitation"]["Kind"], "None")
        self.assertEqual(w["Excitation"]["VanishingUnderSetback"], [island["Body"]])
        self.assertIsNone(w["ProposedTerminal"])

    def test_window_without_the_island(self):
        r = WQ.inventory(self.manifest, {"W": (30.0, 90.0, 50.0, 70.0)}, 2.0, WQ.DEFAULT_WEIGHTS)
        w = r["Windows"]["W"]
        self.assertEqual(w["ByType"], {"SameConductorStrip": {"Features": 1, "LengthInside": 80.0, "Proxy": 80.0 * 2.4}})
        # the strip's two chains are open (no closed loop): two OpenChain records, no body, no terminal
        self.assertEqual(w["ConductorCount"], 2)
        self.assertEqual({b["Kind"] for b in w["Bodies"]}, {"OpenChain"})
        self.assertTrue(all(a["Role"].startswith("unresolved") for a in w["TerminalAssignment"]))
        self.assertIsNone(w["ProposedTerminal"])
        pred = w["MeanSeparationPredictions"]
        self.assertEqual(len(pred), 1)
        self.assertAlmostEqual(pred[0]["ChipSeparation"], 2.0)
        self.assertAlmostEqual(pred[0]["WindowInteriorMean"], 2.0, places=9)
        self.assertFalse(pred[0]["KeyMayChange"])

    def test_mean_separation_prediction_flags_a_taper(self):
        # tilt the second strip edge into a taper: separation 2.0 -> 3.0 over the run
        m = json.loads(json.dumps(self.manifest))
        seg = [s for s in m["Identification"]["Segments"] if s["Key"][0][1] == 62.0][0]
        seg["Key"][1][1] = 63.0
        seg["Length"] = math.dist(seg["Key"][0][:2], seg["Key"][1][:2])
        r = WQ.inventory(m, {"W": (30.0, 90.0, 50.0, 70.0)}, 2.0, WQ.DEFAULT_WEIGHTS)
        pred = r["Windows"]["W"]["MeanSeparationPredictions"][0]
        self.assertTrue(pred["KeyMayChange"])
        self.assertGreater(pred["WindowInteriorMax"], pred["WindowInteriorMin"] + 0.5)

    def test_clip_interval(self):
        self.assertEqual(WQ.clip_interval((0.0, 0.0), (10.0, 0.0), (2.0, 5.0, -1.0, 1.0)), (0.2, 0.5))
        self.assertIsNone(WQ.clip_interval((0.0, 5.0), (10.0, 5.0), (2.0, 5.0, -1.0, 1.0)))
        self.assertEqual(WQ.clip_interval((3.0, -5.0), (3.0, 5.0), (2.0, 5.0, -1.0, 1.0)), (0.4, 0.6))

    def test_cli_writes_json_csv_markdown(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "m.json")
            json.dump(self.manifest, open(path, "w"))
            out = os.path.join(tmp, "out.json")
            csv_path = os.path.join(tmp, "out.csv")
            md = os.path.join(tmp, "out.md")
            rc = WQ.main(["--manifest", path, "--window", "W", "-30", "30", "-30", "30", "--output", out, "--csv", csv_path, "--markdown", md])
            self.assertEqual(rc, 0)
            result = json.load(open(out))
            self.assertIn("W", result["Windows"])
            self.assertEqual(len(open(csv_path).read().splitlines()), 1 + result["Windows"]["W"]["FeatureCount"])
            self.assertIn("proposed terminal body", open(md).read())


if __name__ == "__main__":
    unittest.main()
