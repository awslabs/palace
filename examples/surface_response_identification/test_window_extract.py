# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""window_extract (lane W, decision 184 M1): a two-plane hand-built manifest - the plane-0 ground / hole / island of
test_window_query plus a plane-10 ground (outer boundary in the perimeter table, convex corners) with a hole and an island, a bump footprint (closed
NonManifold loops on both planes, 4 um square) landing on the plane-0 island and on the plane-10 GROUND, and a second
bump joining the plane-10 island to nothing (partner missing). The extract must keep the chip's loops, pair the
footprints, join the plane-0 island into the ground conductor, and window_query on the extract must then refuse the
island as a terminal (M1) and fall back to the decision-180 rule (M2)."""
import json
import os
import sys
import tempfile
import unittest

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import window_extract as WE  # noqa: E402
from . import window_query as WQ  # noqa: E402
from .test_window_query import R, Builder, square  # noqa: E402

Z2 = 10.0


def bump(b, x0, y0, z, size=4.0):
    for s in square(x0, y0, x0 + size, y0 + size, z)[0]:
        s = dict(s)
        s["Exclusion"] = {"Class": "NonManifold", "Reason": "bump"}
        b.segments.append(s)


def build_two_planes():
    b = Builder()
    outer, outer_corners = square(-100.0, -100.0, 100.0, 100.0)
    b.add_loop(outer, outer_corners, None, exclusion="SimulationBoundary")
    b.add_loop(*square(-20.0, -20.0, 20.0, 20.0), "ConcaveCorner")
    b.add_loop(*square(-10.0, -10.0, 10.0, 10.0), "ConvexCorner")  # plane-0 island, body A
    b.add_loop(*square(-100.0, -100.0, 100.0, 100.0, Z2), "ConvexCorner")  # plane-10 ground outline (metal inside)
    b.add_loop(*square(40.0, 40.0, 80.0, 80.0, Z2), "ConcaveCorner")
    b.add_loop(*square(50.0, 50.0, 70.0, 70.0, Z2), "ConvexCorner")  # plane-10 island, body B
    bump(b, -2.0, -2.0, 0.0)  # on island A
    bump(b, -2.0, -2.0, Z2)  # partner: on the plane-10 ground -> A is ground
    bump(b, 58.0, 58.0, Z2)  # on island B, no partner -> B stays an island
    return b.manifest()


class WindowExtract(unittest.TestCase):
    def setUp(self):
        self.manifest = build_two_planes()
        self.ident = self.manifest["Identification"]
        self.loops = WQ.PerimeterLoops(self.ident["Segments"], self.ident["Features"], R)

    def test_footprints_and_conductors(self):
        fps = WE.chip_footprints(self.ident["Segments"], self.loops.loops, self.loops.ground_bodies)
        self.assertEqual(len(fps), 3)
        self.assertTrue(all(f["Closed"] for f in fps))
        paired = [f for f in fps if f["Partner"] is not None]
        self.assertEqual(len(paired), 2)
        self.assertEqual({f["Height"] for f in paired}, {Z2})
        by_plane = {(f["Plane"], round(f["Centroid"][0])): f for f in fps}
        island_a = by_plane[(0.0, 0)]
        self.assertEqual(island_a["BodyBy"], "Loop")
        self.assertEqual(WE.body_kind(self.loops.loops, island_a["Body"]), "Island")
        ground_2 = by_plane[(Z2, 0)]
        self.assertEqual(ground_2["BodyBy"], "Loop")
        self.assertEqual(WE.body_kind(self.loops.loops, ground_2["Body"]), "Ground")
        conductors, joins = WE.chip_conductors(self.loops.loops, fps)
        self.assertEqual(joins, 1)
        self.assertEqual(conductors.find(island_a["Body"]), conductors.find(ground_2["Body"]))
        facts = WE.conductor_facts(self.loops.loops, fps, conductors)
        self.assertTrue(facts[conductors.find(island_a["Body"])]["Ground"])
        self.assertEqual(facts[conductors.find(island_a["Body"])]["Bodies"], 2)
        island_b = by_plane[(Z2, 60)]
        self.assertFalse(facts[conductors.find(island_b["Body"])]["Ground"])

    def test_extract_then_inventory(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "m.json")
            json.dump(self.manifest, open(path, "w"))
            out = os.path.join(tmp, "x")
            summary = os.path.join(tmp, "s.json")
            rc = WE.main(["--manifest", path, "--chip", "t", "--window", "A", "-30", "30", "-30", "30",
                          "--window", "B", "45", "75", "45", "75", "--halo", "5", "--output-dir", out, "--summary", summary])
            self.assertEqual(rc, 0)
            s = json.load(open(summary))
            self.assertEqual((s["Footprints"], s["FootprintsPaired"], s["BumpJoins"]), (3, 2, 1))
            a = json.load(open(os.path.join(out, "t-A.json")))
            b = json.load(open(os.path.join(out, "t-B.json")))
        # window A: the island A is bump-joined to the plane-10 ground -> a ground conductor, no terminal (M1)
        self.assertEqual(a["WindowExtract"]["Window"], "A")
        self.assertTrue(all(WQ.clip_interval(seg["Key"][0], seg["Key"][1], (-35, 35, -35, 35)) is not None
                            for seg in a["Identification"]["Segments"]))
        self.assertEqual(len(a["WindowExtract"]["Footprints"]), 2)
        ra = WQ.inventory(a, {"A": (-30.0, 30.0, -30.0, 30.0)}, 2.0, WQ.DEFAULT_WEIGHTS)["Windows"]["A"]
        island = [r for r in ra["Bodies"] if r["Kind"] == "Island"][0]
        self.assertTrue(island["EntirelyInside"])  # geometrically whole ...
        conductor = [g for g in ra["Conductors"] if island["Body"] in g["Bodies"]][0]
        self.assertEqual(conductor["Kind"], "Ground")  # ... but a ground under the bump rule
        self.assertIsNone(ra["ProposedTerminal"])
        self.assertEqual(ra["Excitation"]["Kind"], "None")
        self.assertFalse(ra["Excitation"]["Realisable"])
        roles = {r["Body"]: r["Role"] for r in ra["TerminalAssignment"]}
        self.assertIn("bump-joined", roles[island["Body"]])
        self.assertEqual(len(ra["BumpFootprints"]), 2)
        self.assertEqual({f["Height"] for f in ra["BumpFootprints"]}, {Z2})
        # window B: island B has an unpaired footprint -> still an island, the quoted terminal
        rb = WQ.inventory(b, {"B": (45.0, 75.0, 45.0, 75.0)}, 2.0, WQ.DEFAULT_WEIGHTS)["Windows"]["B"]
        self.assertEqual(rb["Excitation"]["Kind"], "Island")
        island_b = [r for r in rb["Bodies"] if r["Kind"] == "Island"][0]
        self.assertEqual(rb["ProposedTerminal"], island_b["Body"])
        # the extract's loop model is the chip's: the plane-10 hole loop still has its chip parent / depth
        hole = [l for l in b["WindowExtract"]["Loops"] if l["Plane"] == Z2 and not l["MetalInside"] and l["Depth"] == 1]
        self.assertEqual(len(hole), 1)


if __name__ == "__main__":
    unittest.main()
