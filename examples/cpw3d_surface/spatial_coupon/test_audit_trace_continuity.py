# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest
import numpy as np
from collections import Counter
from audit_trace_continuity import audit
from generate_spatial_response import cap_ring, connect_rings


class TraceContinuityTest(unittest.TestCase):
    def test_conforming_linear_trace(self):
        triangles = {1: [((0., 0., 0.), 0.), ((1., 0., 0.), 1.), ((0., 1., 0.), 0.)],
                     2: [((0., 0., 0.), 0.), ((.5, 0., 0.), .5), ((0., 0., 1.), 0.)]}
        self.assertTrue(audit(triangles)["ContinuousAtAllVertices"])

    def test_caps_retain_collinear_boundary_nodes(self):
        xy = [(0., 0.), (.2, 0.), (1., 0.), (1., 1.), (.8, 1.),
              (.5, 1.), (0., 1.), (0., .5), (0., .25)]
        size = len(xy)
        points = np.asarray([(x, y, z) for z in (0., 1.) for x, y in xy])
        triangles = []
        connect_rings(triangles, 0, size, size)
        cap_ring(triangles, points, 0, size, True)
        cap_ring(triangles, points, size, size, False)
        self.assertEqual(len(triangles), 4 * size - 4)
        counts = Counter(tuple(sorted((tri[k], tri[(k + 1) % 3])))
                         for tri in triangles for k in range(3))
        self.assertTrue(all(count == 2 for count in counts.values()))
        for source in range(len(points)):
            trace = {i + 1: [(tuple(points[v]), float(v == source)) for v in tri]
                     for i, tri in enumerate(triangles)}
            self.assertTrue(audit(trace)["ContinuousAtAllVertices"], source)

    def test_side_cap_tjunction_jump(self):
        triangles = {1: [((0., 0., 0.), 0.), ((1., 0., 0.), 0.), ((0., 1., 0.), 0.)],
                     2: [((0., 0., 0.), 0.), ((.5, 0., 0.), 1.), ((0., 0., 1.), 0.)]}
        result = audit(triangles)
        self.assertFalse(result["ContinuousAtAllVertices"])
        self.assertAlmostEqual(result["Worst"]["AbsoluteJump"], 1.)
        self.assertEqual(result["Worst"]["Point"], [.5, 0., 0.])


if __name__ == "__main__":
    unittest.main()
