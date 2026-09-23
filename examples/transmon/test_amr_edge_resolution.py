#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import importlib.util
import tempfile
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("amr_edge_resolution.py")
SPEC = importlib.util.spec_from_file_location("amr_edge_resolution", MODULE_PATH)
TOOL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TOOL)

# Two root tetrahedra sharing the face 0-1-2 in the plane z = 0 (a "metal" sheet of attribute 5), the
# upper one refined once (ref_type 7: the four corner children are enough for the reader, which only
# needs leaves), and the refined vertices 6..11 as midpoints via vertex_parents. Node 4 is the apex
# below, 5 the apex above. The boundary section lists root faces, as in the meshes Palace writes: the
# non-metal face 1-2-5 (attribute 8) shares the hypotenuse with the metal sheet; the faces 0-1-4 and
# 0-2-4 (attribute 3) lie on the bounding box y = 0 and x = 0.
NC_MESH = """MFEM NC mesh v1.0

dimension
3

# rank attr geom ref_type nodes/children
elements
6
0 1 4 0 0 1 2 4
0 2 4 7 2 3 4 5
0 2 4 0 0 6 7 8
0 2 4 0 6 1 9 10
0 2 4 0 7 9 2 11
0 2 4 0 8 10 11 5

# attr geom nodes
boundary
4
5 2 0 1 2
8 2 1 2 5
3 2 0 1 4
3 2 0 2 4

# vert_id p1 p2
vertex_parents
6
6 0 1
7 0 2
8 0 5
9 1 2
10 1 5
11 2 5

# top-level node coordinates
coordinates
6
3
0 0 0
10 0 0
0 10 0
0 0 0
0 0 -10
0 0 10
"""


class AmrEdgeResolutionTest(unittest.TestCase):
    def write(self, text):
        self.tmp = tempfile.TemporaryDirectory()
        path = Path(self.tmp.name, "adapted.mesh")
        path.write_text(text)
        return path

    def tearDown(self):
        if hasattr(self, "tmp"):
            self.tmp.cleanup()

    def test_reads_leaves_boundary_and_midpoint_vertices(self):
        leaves, boundary, coordinates = TOOL.read_nc_mesh(self.write(NC_MESH))
        self.assertEqual(leaves.shape, (5, 4))
        self.assertEqual(boundary.shape, (4, 4))
        self.assertEqual(coordinates[6].tolist(), [5.0, 0.0, 0.0])
        self.assertEqual(coordinates[10].tolist(), [5.0, 0.0, 5.0])
        self.assertEqual(coordinates[11].tolist(), [0.0, 5.0, 5.0])

    def test_metal_edges_are_shared_with_non_metal_faces_off_the_bounding_box(self):
        leaves, boundary, coordinates = TOOL.read_nc_mesh(self.write(NC_MESH))
        segments, truncation = TOOL.metal_edge_segments(boundary, {5}, coordinates)
        # Edge 0-1 lies on the bounding-box face y = 0 (shared with the face 0-1-4), edge 0-2 on the
        # face x = 0 (shared with 0-2-4); only the hypotenuse 1-2 is a physical metal edge.
        self.assertEqual(segments.tolist(), [[1, 2]])
        self.assertEqual(sorted(truncation.tolist()), [[0, 1], [0, 2]])

    def test_band_statistics_follow_the_centroid_distance(self):
        result = TOOL.analyze(self.write(NC_MESH), {5}, [2.5, 0.5], 2.5)
        self.assertEqual(result["LeafElements"], 5)
        self.assertEqual(result["MetalEdgeSegments"], 1)
        self.assertAlmostEqual(result["MetalEdgeLength"], 10.0 * 2 ** 0.5)
        self.assertAlmostEqual(result["GlobalMinLongestEdge"], 5.0 * 2 ** 0.5)
        # The child tetrahedra 6-1-9-10 and 7-9-2-11 touch the hypotenuse; their centroids
        # (6.25, 1.25, 1.25) and (1.25, 6.25, 1.25) are sqrt(2.5^2 / 2 + 1.25^2) = 2.165 from the
        # edge 1-2: inside the 2.5 band only. The lower root tetrahedron's centroid is 4.33 away.
        by_radius = {band["Radius"]: band for band in result["Bands"]}
        self.assertEqual(by_radius[2.5]["Elements"], 2)
        self.assertEqual(by_radius[0.5]["Elements"], 0)
        self.assertAlmostEqual(by_radius[2.5]["MinLongestEdge"], 5.0 * 2 ** 0.5)
        self.assertEqual(result["TouchingEdge"]["Elements"], 2)

    def test_rejects_meshes_that_are_not_nonconforming_tetrahedra(self):
        with self.assertRaises(ValueError):
            TOOL.read_nc_mesh(self.write("MFEM mesh v1.0\n\ndimension\n3\n"))
        with self.assertRaises(ValueError):
            TOOL.read_nc_mesh(self.write(NC_MESH.replace("0 1 4 0 0 1 2 4", "0 1 5 0 0 1 2 4")))


if __name__ == "__main__":
    unittest.main()
