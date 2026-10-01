# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""window_polygons on a hand-built single-plane mesh: ground | gap | trace 1 | gap | trace 2 | gap | ground (rectangles
of two triangles each, nodes shared along the common edges), the matching extract (segments = the six metal edges,
loops / bodies: grounds -1, traces 1 and 2) and an E0 excitation. Trace 1 is the open-terminated terminal (metal stopped
``OpenSetback`` inside the walls), trace 2 is a cut non-ground conductor -> bridged to the right ground along both
walls (the gap 30 < x < 40 becomes metal over the bridge width), the gap next to the terminal is left open; the
verification against the extract passes (every boundary edge off the walls lies on a segment)."""
import json
import math
import os
import sys
import tempfile
import unittest

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import window_polygons as WP  # noqa: E402

METAL, GAP = 5, 8
COLUMNS = [(-100.0, -10.0, METAL), (-10.0, -5.0, GAP), (-5.0, 5.0, METAL), (5.0, 20.0, GAP), (20.0, 30.0, METAL),
           (30.0, 40.0, GAP), (40.0, 100.0, METAL)]
Y0, Y1 = -100.0, 100.0


JITTER = math.pi * 1e-5  # a coordinate offset with more than 6 decimals (3.14159e-5 um)
ISLAND = (8.0, 17.0, -20.0, 20.0)  # an uncut island inside the gap column 5 < x < 20 (island=True)


def build(tmp, jitter=0.0, island=False):
    """Column edges shifted by ``jitter`` (every metal edge then sits on coordinates with more than 6 decimals); with
    ``island`` the gap column 5..20 is re-triangulated around a 9 x 40 metal island (the gap rectangles below / above / left /
    right of it share the island's corner nodes)."""
    node_tag = {}

    def tag(x, y):
        return node_tag.setdefault((x, y), len(node_tag) + 1)

    nodes, xyz, attribute = [], [], []

    def rectangle(x0, x1, y0, y1, attr):
        corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
        for tri in ((0, 1, 2), (0, 2, 3)):
            nodes.append([tag(*corners[i]) for i in tri])
            xyz.append([[corners[i][0], corners[i][1], 0.0] for i in tri])
            attribute.append(attr)

    ix0, ix1, iy0, iy1 = (v + jitter for v in ISLAND)
    for x0, x1, attr in COLUMNS:
        x0, x1 = x0 + jitter, x1 + jitter
        if island and attr == GAP and x0 < ix0 < x1:
            rectangle(x0, x1, Y0, iy0, GAP)
            rectangle(x0, x1, iy1, Y1, GAP)
            rectangle(x0, ix0, iy0, iy1, GAP)
            rectangle(ix1, x1, iy0, iy1, GAP)
            rectangle(ix0, ix1, iy0, iy1, METAL)
            continue
        rectangle(x0, x1, Y0, Y1, attr)
    tris = os.path.join(tmp, "tris.npz")
    np.savez(tris, nodes=np.asarray(nodes, np.int64), xyz=np.asarray(xyz, float), attribute=np.asarray(attribute, np.int32),
             windows=np.zeros((1, 4)), names=np.asarray("{}"), mesh=np.asarray("synthetic"))
    edges = [x + jitter for x in (-10.0, -5.0, 5.0, 20.0, 30.0, 40.0)]
    segments = [{"Key": [[x, Y0, 0.0], [x, Y1, 0.0]], "Length": Y1 - Y0, "Chain": i, "Portions": []} for i, x in enumerate(edges)]
    loops = [{"Index": 0, "Plane": 0.0, "Closed": False, "Body": -1, "Segments": [0], "MetalInside": None, "Parent": None, "Depth": 0, "Area": 0.0, "BBox": [-10, -10, Y0, Y1], "Polygon": [], "CornerVotes": [0, 0]},
             {"Index": 1, "Plane": 0.0, "Closed": True, "Body": 1, "Segments": [1, 2], "MetalInside": True, "Parent": 0, "Depth": 1, "Area": 2000.0, "BBox": [-5, 5, Y0, Y1], "Polygon": [], "CornerVotes": [0, 0]},
             {"Index": 2, "Plane": 0.0, "Closed": True, "Body": 2, "Segments": [3, 4], "MetalInside": True, "Parent": 0, "Depth": 1, "Area": 2000.0, "BBox": [20, 30, Y0, Y1], "Polygon": [], "CornerVotes": [0, 0]},
             {"Index": 3, "Plane": 0.0, "Closed": False, "Body": -1, "Segments": [5], "MetalInside": None, "Parent": None, "Depth": 0, "Area": 0.0, "BBox": [40, 40, Y0, Y1], "Polygon": [], "CornerVotes": [0, 0]}]
    bodies = {"-1": {"Kind": "Ground", "Conductor": -1, "ConductorBodies": 1, "ConductorGround": True, "Plane": 0.0},
              "1": {"Kind": "Island", "Conductor": 1, "ConductorBodies": 1, "ConductorGround": False, "Plane": 0.0},
              "2": {"Kind": "Island", "Conductor": 2, "ConductorBodies": 1, "ConductorGround": False, "Plane": 0.0}}
    if island:
        corners = [(ix0, iy0), (ix1, iy0), (ix1, iy1), (ix0, iy1)]
        for k in range(4):
            a, b = corners[k], corners[(k + 1) % 4]
            segments.append({"Key": [[a[0], a[1], 0.0], [b[0], b[1], 0.0]], "Length": math.dist(a, b), "Chain": 6 + k, "Portions": []})
        loops.append({"Index": 4, "Plane": 0.0, "Closed": True, "Body": 3, "Segments": [6, 7, 8, 9], "MetalInside": True, "Parent": 0,
                      "Depth": 1, "Area": 360.0, "BBox": [ix0, ix1, iy0, iy1], "Polygon": [list(c) for c in corners], "CornerVotes": [4, 0]})
        bodies["3"] = {"Kind": "Island", "Conductor": 3, "ConductorBodies": 1, "ConductorGround": False, "Plane": 0.0}
    extract = {"Version": 2, "Identification": {"MatchingRadius": 2.0, "Segments": segments, "Features": []},
               "WindowExtract": {"Loops": loops, "Bodies": bodies, "Footprints": [], "GroundBodies": {"0.0": -1}}}
    extract_path = os.path.join(tmp, "extract.json")
    json.dump(extract, open(extract_path, "w"))
    e0_path = os.path.join(tmp, "e0.json")
    json.dump({"Windows": {"W": {"Excitation": {"Kind": "OpenTerminatedTrace", "Conductor": 1, "Bodies": [1], "Planes": [0.0],
                                                 "OpenSetback": 10.0, "Realisable": True}}}}, open(e0_path, "w"))
    chip_path = os.path.join(tmp, "chip.json")
    json.dump({"Name": "synthetic", "Planes": [{"Name": "L1", "SurfaceZ": 0.0, "Facing": "up", "Attributes": [METAL], "SubstrateThickness": 525.0}],
               "Bump": [], "Exterior": [], "Process": {"MetalThickness": 0.1, "Overetch": 0.05}, "Vacuum": {"Below": 0.0, "Above": 1000.0}},
              open(chip_path, "w"))
    return tris, chip_path, extract_path, e0_path


class WindowPolygons(unittest.TestCase):
    def test_clip_and_split(self):
        tri = [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0)]
        half = WP.clip_half_plane(tri, 0, 5.0, False)
        self.assertAlmostEqual(WP.signed_area(half), 50.0 - 12.5)
        parts = WP.split_by_line(tri, 1, 2.0)
        self.assertEqual(len(parts), 2)
        self.assertAlmostEqual(sum(WP.signed_area(p) for p in parts), 50.0)

    def test_termination_rule_and_verification(self):
        with tempfile.TemporaryDirectory() as tmp:
            tris, chip, extract, e0 = build(tmp)
            out, ver = os.path.join(tmp, "out.json"), os.path.join(tmp, "ver.json")
            rc = WP.main(["--triangles", tris, "--chip", chip, "--extract", extract, "--e0", e0, "--window", "W", "-50", "50", "-50", "50",
                          "--bridge-width", "6", "--halo", "10", "--output", out, "--verification", ver])
            self.assertEqual(rc, 0)
            result = json.load(open(out))
            verification = json.load(open(ver))
        self.assertTrue(verification["Passed"])
        self.assertEqual(len(verification["Warnings"]), 2)  # trace 2 next to the terminal along both walls: not bridged there
        self.assertTrue(all("adjacent to the terminal" in w for w in verification["Warnings"]))
        plane = verification["Planes"]["L1"]
        self.assertEqual(plane["UnmatchedLength"], 0.0)
        self.assertEqual(len(plane["Bridges"]), 2)  # gap 30..40 at both walls; the gap 5..20 next to the terminal stays open
        # the bridge runs 3 R = 6 um into trace 2's and into the right ground's wall intervals (runs of contact, not vertices)
        self.assertEqual(sorted((b[0], b[1]) for b in plane["Bridges"]), [(24.0, 46.0), (24.0, 46.0)])
        self.assertEqual(result["Terminals"], ["trace_1"])
        polygons = result["Planes"][0]["Polygons"]
        self.assertEqual(len(polygons), 3)
        by_label = {}
        for p in polygons:
            xs = [v[0] for v in p["Outer"]]
            ys = [v[1] for v in p["Outer"]]
            by_label.setdefault(p["Conductor"], []).append((min(xs), max(xs), min(ys), max(ys), abs(WP.signed_area([tuple(v) for v in p["Outer"]]))))
        self.assertEqual(by_label["trace_1"], [(-5.0, 5.0, -40.0, 40.0, 800.0)])  # open-terminated 10 um inside both walls
        grounds = sorted(by_label["ground"])
        self.assertEqual(grounds[0][:4], (-50.0, -10.0, -50.0, 50.0))
        # trace 2 + two bridges + the right ground fused into one ground polygon (outer 30 x 100) whose enclosed gap
        # (10 x 88) is a hole, emitted with the outer's (counter-clockwise) orientation
        self.assertEqual(grounds[1][:4], (20.0, 50.0, -50.0, 50.0))
        self.assertAlmostEqual(grounds[1][4], 3000.0)
        fused = [p for p in polygons if p["Conductor"] == "ground" and max(v[0] for v in p["Outer"]) == 50.0][0]
        self.assertEqual(len(fused["Holes"]), 1)
        self.assertAlmostEqual(WP.signed_area([tuple(v) for v in fused["Holes"][0]]), 10.0 * 88.0)
        self.assertGreater(WP.signed_area([tuple(v) for v in fused["Outer"]]), 0.0)
        self.assertEqual(sum(len(p["Holes"]) for p in polygons), 1)
        self.assertEqual(result["Version"], 1)
        self.assertEqual(result["Box"], {"X": [-50.0, 50.0], "Y": [-50.0, 50.0]})
        # decision 190 (b): the smallest gap between two polygons of the plane is the 5-um slot next to the terminal
        self.assertAlmostEqual(plane["IntraPlaneClearance"]["Min"], 5.0)
        self.assertEqual(plane["IntraPlaneClearance"]["Violations"], [])
        self.assertEqual(verification["CrossPlaneCoincidence"]["Runs"], 0)  # one plane

    def test_decision_190_checks(self):
        box = (0.0, 100.0, 0.0, 100.0)
        square = lambda x0, y0, x1, y1: {"Outer": [[x0, y0], [x1, y0], [x1, y1], [x0, y1]], "Holes": []}
        # (b) two polygons 0.01 um apart are a violation at 0.05 um, a 1-um slot is not
        best, violations = WP.intra_plane_clearance([square(10, 10, 20, 20), square(20.01, 10, 30, 20)], box, 0.05)
        self.assertAlmostEqual(best, 0.01)
        self.assertEqual(len(violations), 8)  # the two vertices of each facing side, each against the facing + the adjacent edge
        best, violations = WP.intra_plane_clearance([square(10, 10, 20, 20), square(21, 10, 30, 20)], box, 0.05)
        self.assertAlmostEqual(best, 1.0)
        self.assertEqual(violations, [])
        # (c) the L2 edge y = 50.03 from x 30 to 70 is nominally coincident with the L1 edge y = 50 from x 20 to 60: the run is
        # the overlap 30..60 (30 um) at offset 0.03; the far edges of both squares coincide with nothing
        l1 = [square(20, 10, 60, 50)]
        l2 = [square(30, 50.03, 70, 90)]
        c = WP.cross_plane_coincidence({"L1": l1, "L2": l2}, box, 0.5)
        self.assertEqual(c["Runs"], 1)
        self.assertAlmostEqual(c["Length"], 30.0)
        self.assertAlmostEqual(c["MaxOffset"], 0.03)
        self.assertEqual(WP.intra_plane_clearance(l1, box, 0.05), (None, []))  # one polygon: nothing within the search radius
        self.assertEqual(WP.clearance_label(None), "> 10 um")  # reported as "> search radius", not as a missing value
        self.assertEqual(WP.clearance_label(0.25), 0.25)

    def test_full_precision_vertices(self):
        """Review of decision 195, M1: a chip mesh node is emitted with its original double (bitwise), never rounded; only the
        clipping-introduced vertices lie off the nodes, and those lie exactly on a wall / setback / bridge line."""
        with tempfile.TemporaryDirectory() as tmp:
            tris, chip, extract, e0 = build(tmp, jitter=JITTER, island=True)
            out, ver = os.path.join(tmp, "out.json"), os.path.join(tmp, "ver.json")
            rc = WP.main(["--triangles", tris, "--chip", chip, "--extract", extract, "--e0", e0, "--window", "W", "-50", "50", "-50", "50",
                          "--bridge-width", "6", "--halo", "10", "--output", out, "--verification", ver])
            self.assertEqual(rc, 0)
            result = json.load(open(out))
            verification = json.load(open(ver))
            data = np.load(tris)
        self.assertTrue(verification["Passed"])
        chip_vertices = {(float(p[0]), float(p[1])) for tri in data["xyz"] for p in tri}
        self.assertTrue(all(round(x, 6) != x for x, _ in chip_vertices if abs(x) < 50.0))  # the edges carry > 6 decimals
        plane = verification["Planes"]["L1"]
        lines = {-50.0, 50.0, -40.0, 40.0} | {v for b in plane["Bridges"] for v in b[:4]}
        emitted, on_nodes = [], 0
        for poly in result["Planes"][0]["Polygons"]:
            for loop in [poly["Outer"]] + poly["Holes"]:
                for v in loop:
                    v = tuple(v)
                    emitted.append(v)
                    if v in chip_vertices:
                        on_nodes += 1
                    else:
                        self.assertTrue(v[0] in lines or v[1] in lines, f"vertex {v} is neither a chip node nor on a clip line")
        self.assertGreater(on_nodes, 0)
        self.assertTrue(any(round(v[0], 6) != v[0] for v in emitted))  # unrounded doubles are emitted
        # the uncut island (a non-excited conductor, label ground) is exactly the chip's four nodes, bitwise
        island = [p for p in result["Planes"][0]["Polygons"] if min(v[0] for v in p["Outer"]) > 5.0 and max(v[0] for v in p["Outer"]) < 20.0]
        self.assertEqual(len(island), 1)
        corners = {(ISLAND[0] + JITTER, ISLAND[2] + JITTER), (ISLAND[1] + JITTER, ISLAND[2] + JITTER),
                   (ISLAND[1] + JITTER, ISLAND[3] + JITTER), (ISLAND[0] + JITTER, ISLAND[3] + JITTER)}
        self.assertEqual({tuple(v) for v in island[0]["Outer"]}, corners)
        self.assertTrue(corners <= chip_vertices)
        self.assertEqual(island[0]["Conductor"], "ground")
        # verification scope (m1): the setback strips exclude the terminal's edges only (its end caps on y = +-40: 2 x 10 um);
        # the bridge strips exclude the bridges' inner edges (y = +-44, 2 x 10 um; the trace-2 / ground edges inside the
        # rectangles are interior to the union); trace 2's and the grounds' edges inside the setback band are matched against
        # the chip (x = -10 / 20 over the full 100 um, x = 30 / 40 over 88 um outside the bridges; +-1e-6 strip tolerance)
        self.assertEqual(plane["SetbackExclusionScope"], "terminal conductor edges only")
        self.assertAlmostEqual(plane["ExcludedLength (bridges / terminal setback)"], 2 * 10.0 + 2 * 10.0, places=4)
        self.assertAlmostEqual(plane["MatchedLength"], 100.0 + 80.0 + 80.0 + 100.0 + 2 * 88.0 + 2 * (9.0 + 40.0), places=4)
        self.assertEqual(plane["UncutSegments"], 4)  # the island's four edges, covered over their whole length
        self.assertEqual(plane["IntraPlaneClearance"]["Label"], plane["IntraPlaneClearance"]["Min"])
        self.assertAlmostEqual(plane["IntraPlaneClearance"]["Min"], 3.0)  # island to trace 1 / trace 2


if __name__ == "__main__":
    unittest.main()
