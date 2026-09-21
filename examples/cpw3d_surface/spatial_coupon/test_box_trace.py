# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import csv
from pathlib import Path
import unittest
import numpy as np
from box_trace import complete_box_trace, validate_box_trace
from generate_spatial_response import (NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE, build_matching_surface,
                                       cap_triangulation_report, delaunay_flip_cap)


def triangle_altitudes(points, triangles):
    """(minimum altitude, shortest edge) of every triangle."""
    rows = []
    for triangle in triangles:
        a, b, c = (points[v] for v in triangle)
        lengths = (np.linalg.norm(b - a), np.linalg.norm(c - b), np.linalg.norm(a - c))
        rows.append((np.linalg.norm(np.cross(b - a, c - a)) / max(lengths), min(lengths)))
    return rows


class BoxTraceTest(unittest.TestCase):
    def test_missing_corner_is_added_without_dropping_old_nodes(self):
        xy=[(-8.,-8.),(2.,-8.),(6.,-7.),(6.,8.),(-8.,8.)]
        points=np.asarray([(x,y,z) for z in (-2.05,-.05,0.,.1,2.1) for x,y in xy])
        new,triangles,groups,report=complete_box_trace(points)
        self.assertEqual(report['AddedVertices'],5)
        self.assertEqual(groups,[6]*5)
        self.assertTrue(np.array_equal(new[report['RetainedVertexMap']],points))
        self.assertEqual(report['OffBoxTriangles'],0)
        self.assertTrue(report['ClosedOrientedSurface'])
        self.assertLess(report['MaximumRelativeFaceAreaDifference'],1e-12)

    def test_non_square_generator_always_includes_all_corners(self):
        bounds=np.asarray([[-8.,-8.,-2.05],[6.,8.,2.1]])
        points,triangles,_=build_matching_surface(bounds,[-2.05,0.,2.1],8,[],2.,.1,90.,[])
        result=validate_box_trace(points,triangles,bounds)
        self.assertTrue(result['ClosedOrientedSurface'])
        for z in (-2.05,0.,2.1):
            self.assertIn((6.,-8.,z),set(map(tuple,points)))

    def test_corner_cut_is_rejected(self):
        points=np.asarray([[0.,0.,0.],[1.,.1,0.],[0.,0.,1.],[1.,1.,1.]])
        with self.assertRaisesRegex(ValueError,'matching-box face'):
            validate_box_trace(points,np.asarray([[0,1,2]]))

    def test_delaunay_caps_remove_the_needle_ears_of_the_ten_edge_basis(self):
        # The gallery ten-edge basis (testdata; process-library frame, y vertical): its two
        # ear-clipped caps carry 17 needles each (a 56 nm ring edge joined to a vertex
        # 19 um away along the side: altitude 5.5 nm). Lawson flips to the Delaunay caps
        # keep every vertex and the closed orientation; the altitudes below FarSize /
        # TraceBasisSizeRatio (0.32 um) become ring spacings (times the sine of the join
        # angle), a 5.5 nm -> 32 nm minimum.
        directory = Path(__file__).parent / 'testdata' / 'ten-edge-6791f1c84123'
        with (directory / 'trace-vertices.csv').open(newline='') as stream:
            rows = list(csv.DictReader(stream))
        points = np.asarray([[float(row['x']), float(row['z']), float(row['y'])] for row in rows])   # caps in xy
        with (directory / 'trace-triangles.csv').open(newline='') as stream:
            triangles = [tuple(int(row[k]) - 1 for k in ('vertex_i', 'vertex_j', 'vertex_k')) for row in csv.DictReader(stream)]
        bounds = np.asarray([points.min(axis=0), points.max(axis=0)])
        needles = lambda rows: sum(a < NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE * e for a, e in rows if a < 0.32)
        self.assertEqual(needles(triangle_altitudes(points, triangles)), 34)
        flipped = list(triangles)
        flips = 0
        for level in bounds[:, 2]:
            cap = [i for i, t in enumerate(triangles) if all(abs(points[v, 2] - level) < 1e-9 for v in t)]
            self.assertEqual(cap, list(range(cap[0], cap[-1] + 1)))            # one contiguous cap block
            rest = flipped[:cap[0]]
            block = flipped[cap[0]:cap[-1] + 1]
            flips += delaunay_flip_cap(block, points, 0)
            flipped = rest + block + flipped[cap[-1] + 1:]
        self.assertGreater(flips, 0)
        self.assertEqual(len(flipped), len(triangles))
        self.assertEqual(sorted(set(v for t in flipped for v in t)), list(range(len(points))))
        for connectivity in (triangles, flipped):
            self.assertTrue(validate_box_trace(points, np.asarray(connectivity), bounds)['ClosedOrientedSurface'])
        after = triangle_altitudes(points, flipped)
        # The residual needles by the 0.6 rule are the 0.107 um cluster edges at the box's
        # +x side joined obliquely to the nearest far vertex (altitude 0.107 x sin 35 deg =
        # 0.062, marginally below 0.6 x 0.107) and their +z mirror: 8 of 34, at 11x the
        # ear-clipped minimum altitude.
        self.assertLessEqual(needles(after), 8)
        self.assertLess(min(a for a, _ in triangle_altitudes(points, triangles)), 0.006)
        self.assertGreater(min(a for a, _ in after), 0.03)
        self.assertGreater(min(a for a, e in after if a < 0.32 and a < NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE * e), 0.06)
        # Side strips are untouched; the Delaunay caps are a fixed point of the flips.
        on_cap = lambda t: any(all(abs(points[v, 2] - level) < 1e-9 for v in t) for level in bounds[:, 2])
        sides = lambda connectivity: [t for t in connectivity if not on_cap(t)]
        self.assertEqual(sides(flipped), sides(triangles))
        for level in bounds[:, 2]:
            cap = [t for t in flipped if all(abs(points[v, 2] - level) < 1e-9 for v in t)]
            self.assertEqual(delaunay_flip_cap(cap, points, 0), 0)
        report = cap_triangulation_report(points, flipped)
        self.assertAlmostEqual(report['MinimumAltitude'], min(a for a, _ in after))

    def test_build_matching_surface_delaunay_keeps_vertices_and_closure(self):
        bounds = np.asarray([[-8., -8., -2.05], [6., 8., 2.1]])
        clipped = build_matching_surface(bounds, [-2.05, 0., 2.1], 8, [], 2., .1, 90., [])
        delaunay = build_matching_surface(bounds, [-2.05, 0., 2.1], 8, [], 2., .1, 90., [], cap_triangulation='delaunay')
        self.assertTrue(np.array_equal(clipped[0], delaunay[0]) and clipped[2] == delaunay[2])
        self.assertEqual(len(clipped[1]), len(delaunay[1]))
        self.assertTrue(validate_box_trace(delaunay[0], delaunay[1], bounds)['ClosedOrientedSurface'])
        with self.assertRaisesRegex(ValueError, 'cap triangulation'):
            build_matching_surface(bounds, [-2.05, 2.1], 8, [], 2., .1, 90., [], cap_triangulation='fan')


if __name__=='__main__':unittest.main()
