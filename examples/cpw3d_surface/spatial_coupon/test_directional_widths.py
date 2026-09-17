# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The band anisotropy statistics exclude a recorded edge layer's cells and report
them separately (supervisor decision 31)."""
import unittest

import meshio
import numpy as np

from audit_edge_metric_mesh import directional_widths


def band_mesh():
    """Two rows of tetrahedra along a segment on the x axis: a 'layer' row with
    centroids ~2 nm from the segment (1 nm x 10 nm cells) and a 'band' row with
    centroids ~30 nm from it (25 nm cubes split into tets)."""
    points, cells = [], []
    def add_box(x0, dx, y0, dy, z0, dz):
        base = len(points)
        for z in (z0, z0 + dz):
            for y in (y0, y0 + dy):
                for x in (x0, x0 + dx):
                    points.append([x, y, z])
        for tet in ([0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7], [0, 4, 5, 7], [0, 4, 6, 7]):
            cells.append([base + i for i in tet])
    for k in range(40):
        add_box(1.0 + .01 * k, .01, 0., .001, 0., .001)          # layer cells, 1 nm transverse
        add_box(1.0 + .025 * k, .025, .02, .025, .02, .025)      # band cells, isotropic 25 nm
    cells = np.asarray(cells)
    return meshio.Mesh(np.asarray(points, dtype=float), [("tetra", cells)],
                       cell_data={"gmsh:physical": [np.ones(len(cells), dtype=int)]})


class DirectionalWidthsTest(unittest.TestCase):
    def test_edge_layer_cells_are_excluded_and_reported(self):
        mesh = band_mesh()
        recipe = {"PhysicalSegments": [[0., 0., 0., 3., 0., 0.]], "TruePhysicalCorners": [[-9., 0., 0.]],
                  "NormalSize": .025, "TangentialSize": .1}
        plain = directional_widths(mesh, recipe)
        self.assertEqual(plain["3.0"]["Cells"], 480)
        self.assertEqual(plain["3.0"]["ExcludedEdgeLayerCells"], 0)
        self.assertNotIn("EdgeLayer", plain["3.0"])
        # The layer cells pull the plain tangential P50 below the band's 25 nm ...
        self.assertLess(plain["3.0"]["WidthsTangentialTransverse1Transverse2"][1][0], .02)
        layered = directional_widths(mesh, recipe, layer_spans=[[1., 0., 0., 1.4, 0., 0.]],
                                     layer_reach=.005)
        band = layered["3.0"]
        self.assertEqual(band["Cells"], 240)
        self.assertEqual(band["ExcludedEdgeLayerCells"], 240)
        # ... and the band statistics without them are the isotropic 25 nm cells.
        self.assertAlmostEqual(band["WidthsTangentialTransverse1Transverse2"][1][0], .025)
        self.assertLessEqual(band["WidthsTangentialTransverse1Transverse2"][2][1], .025 * np.sqrt(3) + 1e-12)
        layer = band["EdgeLayer"]
        self.assertEqual(layer["Cells"], 240); self.assertEqual(layer["Reach"], .005)
        self.assertAlmostEqual(layer["WidthsTangentialTransverse1Transverse2"][1][0], .01)
        self.assertLess(layer["WidthsTangentialTransverse1Transverse2"][2][2], .0015)
        with self.assertRaises(ValueError):
            directional_widths(mesh, recipe, layer_spans=[[1., 0., 0., 1.4, 0., 0.]], layer_reach=0.)


if __name__ == "__main__":
    unittest.main()
