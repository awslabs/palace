# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest

import meshio
import numpy as np

from edge_layer_census import edge_layer_census


class EdgeLayerCensusTest(unittest.TestCase):
    def test_census_counts_layer_cells_edges_and_surface_rows_against_the_recipe_spans(self):
        # A ridge along x at z = 0.1 (the span), one row 4 nm below it on the sidewall
        # (x-z plane) and two tets: one in the layer, one far away.
        points = np.array([[0., 0., .1], [.0125, 0., .1], [0., 0., .096], [.0125, 0., .096],
                           [.005, .003, .098], [1., 1., 1.], [1.1, 1., 1.], [1., 1.1, 1.], [1., 1., 1.1]])
        cells = [("triangle", np.array([[0, 1, 2], [1, 3, 2]])),
                 ("tetra", np.array([[0, 1, 2, 4], [5, 6, 7, 8]]))]
        recipe = {"EdgeLayer": {"Spans": [[0., 0., .1, .1, 0., .1]], "EdgeSize": .004, "GrowthRatio": 2.,
                                "Aspect": 4., "RowOffsets": [.004, .012, .028], "LayerThickness": .028,
                                "TotalSpanLength": .1}}
        report = edge_layer_census(meshio.Mesh(points, cells), recipe)
        self.assertEqual(report["Tetrahedra"], 2)
        self.assertEqual(report["LayerCells"]["Within0.01"]["Count"], 1)
        self.assertEqual(report["LayerCells"]["Within0.1"]["Count"], 1)
        self.assertEqual(report["LayerCells"]["Within0.01"]["CellsBelow0.01"], 0)
        self.assertEqual(report["SurfaceTrianglesInLayer"], 2)
        self.assertEqual(report["VerticesOnSpans"], 2)
        # The surface row is 4 nm below the ridge: transverse edges of 4 nm in the
        # 2-5 nm shell, tangential spacing 12.5 nm.
        shell = report["SurfaceEdgesByShell"]["2-5nm"]
        self.assertAlmostEqual(shell["TransverseP50Um"], .004)
        self.assertAlmostEqual(shell["TangentialP50Um"], .0125)
        with self.assertRaisesRegex(ValueError, "no edge layer"):
            edge_layer_census(meshio.Mesh(points, cells), {})


if __name__ == "__main__":
    unittest.main()
