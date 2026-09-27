# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""cluster_signature_geometry: the spatial coupon geometry of a version-2 SpatialEdgeCluster
signature (the v2 cluster contract) on synthetic signatures with known masks."""
import json
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


if __name__ == "__main__":
    unittest.main()
