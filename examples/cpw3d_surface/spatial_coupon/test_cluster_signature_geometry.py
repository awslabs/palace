# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""cluster_signature_geometry: the spatial coupon geometry of a version-2 SpatialEdgeCluster
signature (the v2 cluster contract) on synthetic signatures with known masks."""
import json
import math
from pathlib import Path
import sys
import unittest

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cluster_signature_geometry as csg  # noqa: E402

INTERFACES = [{"Slot": 0, "Type": "MA", "Target": 3}, {"Slot": 0, "Type": "MS", "Target": 2}, {"Slot": 0, "Type": "SA", "Target": 1}]
LAW = '{"Type":"PEC"}'


def portion(p, gap, conductor=1):
    return {"Conductor": conductor, "Gap": list(gap), "Interfaces": ["MA", "MS", "SA"], "Law": LAW, "P": list(p)}


def record(portions, vertices=()):
    signature = {"Type": "SpatialEdgeCluster", "EdgeCount": len(portions), "Portions": portions, "Vertices": list(vertices)}
    return {"Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": len(portions), "Signature": signature},
            "Signature": signature, "Interfaces": INTERFACES, "BoundaryCondition": {"Type": "PEC"}}


def mask_area(coupon, conductor):
    total = 0.0
    for facet in coupon["Geometry"]["PlanViewFacets"]:
        if facet["Conductor"] != conductor:
            continue
        pts = np.asarray(facet["Points"])[:, :2]
        total += 0.5 * abs(sum(a[0] * b[1] - a[1] * b[0] for a, b in zip(pts, np.roll(pts, -1, axis=0))))
    return total


class ConcaveCornerClusterTest(unittest.TestCase):
    """Two portions meeting at a concave corner of the metal at the origin (arms toward -x
    and +y, the gap above the horizontal arm and left of the vertical one), free far ends:
    the metal is the box minus its quadrant x <= 0, y >= 0."""
    R = 1.9

    def setUp(self):
        self.record = record([portion((-2.0, 0.0, 0.0, 0.0), (0.0, 1.0)), portion((0.0, 0.0, 0.0, 2.0), (-1.0, 0.0))],
                             [{"P": [0.0, 0.0], "TurnDegrees": 90.0, "Type": "ConcaveCorner"}])
        self.coupon, self.edges = csg.cluster_coupon(self.record, self.R, 0.1, 0.05)

    def test_free_ends_are_lengthened_and_vertex_ends_kept(self):
        rows = self.coupon["Geometry"]["Edges"]
        self.assertEqual(len(rows), 2)
        for row in rows:
            begin, end = row["Interval"]
            # One end at the corner (half length 1.9 = R from the midpoint), the free end at
            # least R from the midpoint: the box rule then extends it by 2R.
            self.assertAlmostEqual(min(-begin, end), 0.5 * 2.0 * self.R, places=9)
            self.assertGreaterEqual(max(-begin, end), self.R - 1e-12)

    def test_model_edges_reproduce_the_portions_exactly(self):
        for edge, entry in zip(self.edges, self.record["Signature"]["Portions"]):
            tangent = np.cross(np.asarray(edge["GapDirection"]), np.asarray(edge["ProcessNormal"]))
            ends = sorted(tuple(np.round((np.asarray(edge["Point"]) + s * tangent)[:2], 9)) for s in edge["Interval"])
            expected = sorted([tuple(np.round(np.asarray(entry["P"][:2]) * self.R, 9)), tuple(np.round(np.asarray(entry["P"][2:]) * self.R, 9))])
            self.assertEqual(ends, expected)

    def test_mask_is_the_box_minus_the_gap_quadrant(self):
        geometry = self.coupon["Geometry"]
        self.assertTrue(geometry["PlanViewFacets"])
        self.assertEqual({f["Conductor"] for f in geometry["PlanViewFacets"]}, {1})
        # The mask area equals the generator's box minus the gap quadrant: recompute the box
        # in the generator's frame (a rotation about z of the canonical frame; corner at 0).
        import generate_spatial_response as spatial_generator
        frame, edges, _ = spatial_generator.normalize_geometry({**self.coupon, "Geometry": {**geometry, "PlanViewFacets": []}}, self.R)
        lower, upper = spatial_generator.coupon_bounds(edges, self.R, 0.1, 0.05)
        gaps = [frame[:2, :2] @ np.asarray(p["Gap"]) for p in self.record["Signature"]["Portions"]]
        widths = []
        for gap in gaps:
            axis = int(np.argmax(np.abs(gap)))
            widths.append((upper[axis] - 0.0) if gap[axis] > 0 else (0.0 - lower[axis]))
        box_area = (upper[0] - lower[0]) * (upper[1] - lower[1])
        self.assertAlmostEqual(mask_area(self.coupon, 1), box_area - widths[0] * widths[1], places=6)
        boundary = geometry["PlanViewBoundary"]
        self.assertEqual([b["Conductor"] for b in boundary], [1])
        self.assertEqual(len(boundary[0]["Segments"]), 6)   # an L: two physical, four box sides


class TwoConductorFingersTest(unittest.TestCase):
    """Two parallel strips of different conductors (a 1 um gap): each conductor's mask is
    its own strip; the signature's conductor labels are kept."""
    R = 1.9

    def test_two_conductors_two_strips(self):
        rec = record([portion((-1.0, -2.0, -1.0, 2.0), (1.0, 0.0), 1), portion((-2.0, -2.0, -2.0, 2.0), (-1.0, 0.0), 1),
                      portion((0.0, -2.0, 0.0, 2.0), (-1.0, 0.0), 2), portion((1.0, -2.0, 1.0, 2.0), (1.0, 0.0), 2)])
        coupon, edges = csg.cluster_coupon(rec, self.R, 0.1, 0.05)
        self.assertEqual({f["Conductor"] for f in coupon["Geometry"]["PlanViewFacets"]}, {1, 2})
        area_1, area_2 = mask_area(coupon, 1), mask_area(coupon, 2)
        self.assertAlmostEqual(area_1, area_2, places=6)
        # Each strip is R wide and spans the whole box height (free ends extended).
        import generate_spatial_response as spatial_generator
        frame, local_edges, _ = spatial_generator.normalize_geometry({**coupon, "Geometry": {**coupon["Geometry"], "PlanViewFacets": []}}, self.R)
        lower, upper = spatial_generator.coupon_bounds(local_edges, self.R, 0.1, 0.05)
        height = (upper - lower)[1] if abs(frame[0][0]) > 0.5 else (upper - lower)[0]
        self.assertAlmostEqual(area_1, self.R * height, places=6)
        self.assertEqual([e["Conductor"] for e in edges], [1, 1, 2, 2])


class FailClosedTest(unittest.TestCase):
    R = 2.0

    def test_crossing_chains_are_refused(self):
        rec = record([portion((-2.0, 0.0, 2.0, 0.0), (0.0, 1.0)), portion((0.0, -2.0, 0.0, 2.0), (1.0, 0.0))])
        with self.assertRaises(csg.SignatureGeometryError):
            csg.cluster_coupon(rec, self.R, 0.1, 0.05)

    def test_inconsistent_metal_sides_are_refused(self):
        # Two parallel edges whose gaps both point inward: the strip between them would be
        # non-metal for both and the outside metal for both - the outside faces then hold
        # metal of the same conductor on both sides of the box, consistent; flipping one gap
        # makes the middle face metal for one edge and gap for the other.
        rec = record([portion((-1.0, -2.0, -1.0, 2.0), (1.0, 0.0)), portion((1.0, -2.0, 1.0, 2.0), (1.0, 0.0))])
        with self.assertRaises(csg.SignatureGeometryError):
            csg.cluster_coupon(rec, self.R, 0.1, 0.05)

    def test_portion_interfaces_must_match_a_slot(self):
        rec = record([portion((-2.0, 0.0, 0.0, 0.0), (0.0, 1.0)), portion((0.0, 0.0, 0.0, 2.0), (1.0, 0.0))],
                     [{"P": [0.0, 0.0], "TurnDegrees": 90.0, "Type": "ConcaveCorner"}])
        rec["Signature"]["Portions"][0]["Interfaces"] = ["MS"]
        with self.assertRaises(csg.SignatureGeometryError):
            csg.cluster_coupon(rec, self.R, 0.1, 0.05)


class TriangulationTest(unittest.TestCase):
    def test_l_shape_ear_clipping_preserves_area(self):
        polygon = [(0, 0), (2, 0), (2, 1), (1, 1), (1, 2), (0, 2)]
        triangles = csg.triangulate(polygon)
        area = sum(0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) for a, b, c in triangles)
        self.assertEqual(len(triangles), 4)
        self.assertAlmostEqual(area, 3.0)


class TransmonSignatureTest(unittest.TestCase):
    """The recorded transmon 3-edge cluster signature at R 1.9 (a 2 um strip meeting the
    ground plane at a concave corner): one conductor, a six-segment boundary."""

    def test_recorded_three_edge_cluster(self):
        signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 3,
                     "Portions": [portion((-0.364831, -0.360984, -0.364831, 2.339595), (-1.0, 0.0)),
                                  portion((-2.364831, -0.360984, -0.364831, -0.360984), (0.0, 1.0)),
                                  portion((0.6878, -3.061562, 0.6878, 2.339595), (1.0, 0.0))],
                     "Vertices": [{"P": [-0.364831, -0.360984], "TurnDegrees": 90.0, "Type": "ConcaveCorner"}]}
        rec = {"Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": 3, "Signature": signature}, "Signature": signature,
               "Interfaces": INTERFACES, "BoundaryCondition": {"Type": "PEC"}}
        coupon, edges = csg.cluster_coupon(rec, 1.9, 0.1, 0.05)
        boundary = coupon["Geometry"]["PlanViewBoundary"]
        self.assertEqual(len(boundary), 1)
        self.assertEqual(len(boundary[0]["Segments"]), 6)
        states = csg.end_states(csg.portions_from_signature(signature, 1.9), csg.vertex_points(signature, 1.9), 1.9)
        self.assertEqual(states, [(False, True), (True, False), (True, True)])
        self.assertEqual(json.dumps(coupon["Geometry"]["Signature"], sort_keys=True), json.dumps(signature, sort_keys=True))


class ArcPortionTest(unittest.TestCase):
    """A rounded finger end (option A): the end edge, two 0.5 R fillet arcs (metal inside) and
    R along both sides. The arcs are chorded at 5 deg / 0.25 R (18 chords per quarter circle):
    every chord end lies on its circle, the chord's gap is radial, the model Edges are the
    chords and the mask covers the conductor."""
    R = 1.0

    def signature(self):
        c45 = 0.5 * math.cos(math.pi / 4)
        top = {"Conductor": 1, "Interfaces": ["MA", "MS", "SA"], "Law": LAW, "P": [-0.5, 0.75, 0.0, 0.25],
               "Arc": [-0.5, 0.25, -0.5 + c45, 0.25 + c45], "GapRadial": 1}
        bottom = {"Conductor": 1, "Interfaces": ["MA", "MS", "SA"], "Law": LAW, "P": [-0.5, -0.75, 0.0, -0.25],
                  "Arc": [-0.5, -0.25, -0.5 + c45, -0.25 - c45], "GapRadial": 1}
        return {"Type": "SpatialEdgeCluster", "EdgeCount": 5,
                "Portions": [portion((0.0, -0.25, 0.0, 0.25), (1.0, 0.0)), top, bottom,
                             portion((-1.5, 0.75, -0.5, 0.75), (0.0, 1.0)), portion((-1.5, -0.75, -0.5, -0.75), (0.0, -1.0))],
                "Vertices": []}

    def test_arcs_are_chorded_on_their_circles(self):
        portions = csg.portions_from_signature(self.signature(), self.R)
        chords = [p for p in portions if p["Portion"] in (1, 2)]
        self.assertEqual(len(portions), 3 + 2 * 18)
        self.assertEqual(len(chords), 36)
        for chord in chords:
            center = np.asarray([-0.5, 0.25 if chord["Portion"] == 1 else -0.25])
            for end in (chord["P0"], chord["P1"]):
                self.assertAlmostEqual(float(np.linalg.norm(end - center)), 0.5, places=12)
            middle = 0.5 * (chord["P0"] + chord["P1"])
            radial = (middle - center) / np.linalg.norm(middle - center)
            self.assertAlmostEqual(float(np.dot(radial, chord["Gap"])), 1.0, places=12)
        # Consecutive chords connect; the arc ends meet the end edge and the sides.
        states = csg.end_states(portions, [], self.R)
        self.assertEqual(sum(free for state in states for free in state), 2)  # the two far side ends

    def test_model_edges_are_the_chords_and_the_mask_covers_the_conductor(self):
        signature = self.signature()
        rec = {"Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": 5, "Signature": signature}, "Signature": signature,
               "Interfaces": INTERFACES, "BoundaryCondition": {"Type": "PEC"}}
        coupon, edges = csg.cluster_coupon(rec, self.R, 0.1, 0.05)
        self.assertEqual(len(edges), 3 + 2 * 18)
        self.assertEqual(len(csg.model_edges(rec, self.R)), len(edges))
        for edge in edges:
            gap = np.asarray(edge["GapDirection"][:2])
            tangent = np.asarray([gap[1], -gap[0]])
            point = np.asarray(edge["Point"][:2])
            for s in edge["Interval"]:
                q = point + s * tangent
                on_arc = min(abs(np.linalg.norm(q - np.asarray([-0.5, y])) - 0.5) for y in (0.25, -0.25))
                on_straight = min(abs(q[0]) if abs(q[1]) <= 0.25 + 1e-9 else np.inf, abs(abs(q[1]) - 0.75) if q[0] <= -0.5 + 1e-9 else np.inf)
                self.assertLess(min(on_arc, on_straight), 1.0e-9)
        self.assertEqual({f["Conductor"] for f in coupon["Geometry"]["PlanViewFacets"]}, {1})
        self.assertGreater(mask_area(coupon, 1), 0.0)
        self.assertEqual(json.dumps(coupon["Geometry"]["Signature"], sort_keys=True), json.dumps(signature, sort_keys=True))

    def test_straight_portion_without_gap_is_refused(self):
        signature = self.signature()
        del signature["Portions"][0]["Gap"]
        with self.assertRaises(csg.SignatureGeometryError):
            csg.portions_from_signature(signature, self.R)


if __name__ == "__main__":
    unittest.main()


def context_piece(p, gap, conductor, chain):
    return {**portion(p, gap, conductor), "Chain": chain}


class DevicePlanCouponTest(unittest.TestCase):
    """Spatial-support contract v3 (decisions 282 / 285): a signature carrying Box + Context
    (the S1p pattern: a lead B on a pad A whose corner lies inside the box, a foreign lead C)
    is built in its own frame inside its Box from the claims plus the context pieces: no
    straight continuation past the pad corner (no fictitious metal), the foreign lead present
    in the mask and listed as ForeignEdges (rule B5), the context rows flagged for the mesher."""
    R = 2.0

    def setUp(self):
        claims = [portion((-3.4, 0.0, 3.4, 0.0), (0.0, 1.0), 1),                 # pad A top edge
                  portion((-0.5, 0.5, 0.5, 0.5), (0.0, -1.0), 2),                # lead B end edge
                  portion((-0.5, 0.5, -0.5, 2.5), (-1.0, 0.0), 2),               # lead B left side
                  portion((0.5, 0.5, 0.5, 2.5), (1.0, 0.0), 2)]                  # lead B right side
        context = [context_piece((3.4, 0.0, 5.0, 0.0), (0.0, 1.0), 1, True),     # pad edge to the corner
                   context_piece((5.0, -2.5, 5.0, 0.0), (1.0, 0.0), 1, True),    # pad right edge (the chain turns)
                   context_piece((-6.4, 0.0, -3.4, 0.0), (0.0, 1.0), 1, True),   # pad edge left continuation
                   context_piece((-0.5, 2.5, -0.5, 6.0), (-1.0, 0.0), 2, True),  # lead B sides to the face
                   context_piece((0.5, 2.5, 0.5, 6.0), (1.0, 0.0), 2, True),
                   context_piece((5.5, 1.0, 6.1, 1.0), (0.0, -1.0), 3, False),   # foreign lead C end edge
                   context_piece((5.5, 1.0, 5.5, 6.0), (-1.0, 0.0), 3, False),   # its sides to the face
                   context_piece((6.1, 1.0, 6.1, 6.0), (1.0, 0.0), 3, False)]
        vertices = [{"P": [-0.5, 0.5], "TurnDegrees": 90.0, "Type": "ConvexCorner"},
                    {"P": [0.5, 0.5], "TurnDegrees": 90.0, "Type": "ConvexCorner"}]
        self.signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 4, "Portions": claims, "Vertices": vertices,
                          "Box": [-6.4, -2.5, 6.4, 6.0], "Context": context}
        self.record = {"Topology": "SpatialEdgeCluster", "Geometry": {"EdgeCount": 4, "Signature": self.signature},
                       "Signature": self.signature, "Interfaces": INTERFACES, "BoundaryCondition": {"Type": "PEC"}}
        self.coupon, self.edges = csg.cluster_coupon(self.record, self.R, 0.1, 0.05)

    def device_edges(self):
        edges = []
        for entry in self.signature["Portions"] + self.signature["Context"]:
            p = [v * self.R for v in entry["P"]]
            edges.append((np.asarray(p[:2]), np.asarray(p[2:])))
        return edges

    def on_device_edge_or_box(self, a, b):
        """Every sample of the segment lies on some device edge (the boundary merges collinear
        device edges into one segment) or on a box face."""
        x0, y0, x1, y1 = self.coupon["Geometry"]["SupportBox"]
        tol = 1e-7 * self.R
        edges = self.device_edges()
        for t in np.linspace(0.0, 1.0, 11):
            point = a + t * (b - a)
            on_box = any(abs(point[c] - v) <= tol for c, v in ((0, x0), (0, x1), (1, y0), (1, y1)))
            on_edge = False
            for (p, q) in edges:
                d = q - p
                L = np.linalg.norm(d)
                s = np.clip(np.dot(point - p, d) / L**2, 0.0, 1.0)
                if np.linalg.norm(point - (p + s * d)) <= tol:
                    on_edge = True
                    break
            if not (on_box or on_edge):
                return False
        return True

    def test_box_frame_and_rows(self):
        geometry = self.coupon["Geometry"]
        self.assertEqual(geometry["SupportBox"], [v * self.R for v in self.signature["Box"]])
        self.assertEqual(geometry["EdgeCount"], 4)
        rows = geometry["Edges"]
        self.assertEqual(len(rows), 4 + 8)
        self.assertEqual(geometry["ContextEdgeCount"], 8)
        # The claims rows are exact (no lengthening: the continuation is in the context).
        for row, claim in zip(rows[:4], self.signature["Portions"]):
            self.assertNotIn("Context", row)
            length = math.dist(claim["P"][:2], claim["P"][2:]) * self.R
            self.assertAlmostEqual(row["Interval"][1] - row["Interval"][0], length, places=9)
        for row, piece in zip(rows[4:], self.signature["Context"]):
            self.assertTrue(row["Context"])
            self.assertEqual(row["Chain"], piece["Chain"])
            self.assertEqual(row["Conductor"], piece["Conductor"])
        # The model's Edges stay the exact claims (A10).
        self.assertEqual(len(self.edges), 4)

    def test_mask_is_the_device_plan_without_fictitious_metal(self):
        facets = self.coupon["Geometry"]["PlanViewFacets"]
        self.assertEqual({f["Conductor"] for f in facets}, {1, 2, 3})
        # Pad A: the box below y = 0 left of the corner x = 5 R (no metal past the corner).
        self.assertAlmostEqual(mask_area(self.coupon, 1), (5.0 + 6.4) * 2.5 * self.R**2, places=6)
        self.assertAlmostEqual(mask_area(self.coupon, 2), 1.0 * (6.0 - 0.5) * self.R**2, places=6)
        self.assertAlmostEqual(mask_area(self.coupon, 3), 0.6 * (6.0 - 1.0) * self.R**2, places=6)
        # Every boundary segment of the mask lies on a device edge or on a box face.
        for component in self.coupon["Geometry"]["PlanViewBoundary"]:
            for first, second in component["Segments"]:
                a = np.asarray(first[:2], dtype=float) * 1e-9 * self.R
                b = np.asarray(second[:2], dtype=float) * 1e-9 * self.R
                self.assertTrue(self.on_device_edge_or_box(a, b), (first, second))

    def test_foreign_edges_are_the_chain_false_pieces(self):
        foreign = self.coupon["Geometry"]["ForeignEdges"]
        self.assertEqual(len(foreign), 3)
        expected = [[v * self.R for v in piece["P"]] for piece in self.signature["Context"] if not piece["Chain"]]
        for segment, piece in zip(foreign, expected):
            self.assertEqual(segment[2], 0.0)
            self.assertEqual(segment[5], 0.0)
            np.testing.assert_allclose([segment[0], segment[1], segment[3], segment[4]], piece, atol=1e-12)
        self.assertEqual(csg.foreign_edge_segments(csg.context_from_signature(self.signature, self.R)), foreign)

    def test_generator_consumes_the_device_plan(self):
        """generate_spatial_response.py --signature-only on the coupon: the identity frame, the
        signature's box, the context rows in mesh-signature.csv, the mask and the classified
        boundary (every face-crossing segment a Continuation)."""
        import subprocess
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "coupon.json").write_text(json.dumps(self.coupon, indent=1) + "\n")
            command = [sys.executable, str(HERE / "generate_spatial_response.py"), str(tmp / "coupon.json"), "--output", str(tmp),
                       "--radius", str(self.R), "--metal-thickness", "0.1", "--overetch-depth", "0.05", "--sidewall-angle", "90",
                       "--top-rounding", "0", "--trench-rounding", "0", "--ring-size", "16", "--cap-triangulation", "delaunay",
                       "--cap-interior-spacing", "0.5", "--order", "1", "--model-name", "device-plan", "--signature-only"]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            rows = (tmp / "mesh-signature.csv").read_text().splitlines()
            self.assertTrue(rows[0].endswith(",Context,Chain"))
            self.assertEqual(len(rows), 1 + 12)
            self.assertEqual([row.split(",")[-2:] for row in rows[1:5]], [["0", "0"]] * 4)
            self.assertEqual([row.split(",")[-2:] for row in rows[5:]], [["1", "1"]] * 5 + [["1", "0"]] * 3)
            boundary = (tmp / "plan-view-boundary.csv").read_text().splitlines()
            self.assertGreater(len(boundary), 1)
            import generate_spatial_response as generator
            frame, edges, facets = generator.normalize_geometry(self.coupon, self.R)
            np.testing.assert_allclose(frame, np.identity(3))
            lower, upper = generator.coupon_bounds(edges, self.R, 0.1, 0.05, self.coupon["Geometry"]["SupportBox"])
            np.testing.assert_allclose(lower[:2], [-6.4 * self.R, -2.5 * self.R])
            np.testing.assert_allclose(upper[:2], [6.4 * self.R, 6.0 * self.R])
            self.assertEqual(sum(1 for e in edges if e.get("Context")), 8)


class LegacyByteIdentityTest(unittest.TestCase):
    """A claims-only (contract-2) signature regenerates the pre-v3 builder's generator inputs
    byte for byte (testdata/legacy-byte-identity: coupon.json, mesh-signature.csv,
    plan-view-mask.csv and plan-view-boundary.csv written by the builder at a79b6af748)."""

    def test_legacy_inputs_are_byte_identical(self):
        import subprocess
        import tempfile
        fixture = HERE / "testdata" / "legacy-byte-identity"
        record = json.loads((fixture / "signature.json").read_text())
        coupon, _ = csg.cluster_coupon(record, 1.9, 0.1, 0.05)
        self.assertEqual(json.dumps(coupon, indent=1) + "\n", (fixture / "coupon.json").read_text())
        self.assertNotIn("SupportBox", coupon["Geometry"])
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "coupon.json").write_text(json.dumps(coupon, indent=1) + "\n")
            command = [sys.executable, str(HERE / "generate_spatial_response.py"), str(tmp / "coupon.json"), "--output", str(tmp),
                       "--radius", "1.9", "--metal-thickness", "0.1", "--overetch-depth", "0.05", "--sidewall-angle", "90",
                       "--top-rounding", "0", "--trench-rounding", "0", "--ring-size", "16", "--cap-triangulation", "delaunay",
                       "--cap-interior-spacing", "0.5", "--order", "1", "--model-name", "legacy-fixture", "--signature-only"]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            for name in ("mesh-signature.csv", "plan-view-mask.csv", "plan-view-boundary.csv"):
                self.assertEqual((tmp / name).read_bytes(), (fixture / name).read_bytes(), name)
