# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest

import meshio
import numpy as np

from edge_layer_census import corner_ball_census, edge_layer_census


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

    def test_corner_ball_census_reports_edges_and_cells_per_shell_of_the_corner_law(self):
        # A corner at the origin: a 4 nm cell at the corner (first shell), a 16 nm
        # cell centred 20 nm out (third shell) and a far cell outside the ball.
        points = np.array([[0., 0., 0.], [.004, 0., 0.], [0., .004, 0.], [0., 0., .004],
                           [.02, 0., 0.], [.036, 0., 0.], [.02, .016, 0.], [.02, 0., .016],
                           [1., 1., 1.], [1.1, 1., 1.], [1., 1.1, 1.], [1., 1., 1.1]])
        mesh = meshio.Mesh(points, [("tetra", np.array([[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]))])
        recipe = {"TruePhysicalCorners": [[0., 0., 0.]], "CornerIsotropyRadius": .1, "NormalSize": .025,
                  "CornerGrading": {"CornerSize": .004, "GrowthRatio": 2., "ShellRadii": [.004, .012, .028, .1],
                                    "ShellSizes": [.004, .008, .016, .025]}}
        report = corner_ball_census(mesh, recipe)
        corner = report["Corners"][0]
        self.assertEqual(corner["Cells"], 2)
        self.assertTrue(corner["PositiveOrientation"])
        self.assertEqual([shell["Cells"] for shell in corner["Shells"]], [1, 0, 1, 0])
        self.assertEqual([shell["TargetSizeUm"] for shell in corner["Shells"]], [.004, .008, .016, .025])
        self.assertEqual(corner["Shells"][0]["Edges"], 6)
        self.assertAlmostEqual(corner["Shells"][0]["EdgeP50Um"], .004 * (1 + np.sqrt(2)) / 2)
        self.assertEqual(corner["Shells"][0]["EdgesOverSqrt2TargetSize"], 0)
        self.assertEqual([shell["Edges"] for shell in corner["Shells"]], [6, 0, 4, 2])
        # Without a grading the ball is one NormalSize shell.
        plain = corner_ball_census(mesh, {k: v for k, v in recipe.items() if k != "CornerGrading"})
        self.assertIsNone(plain["CornerSize"])
        self.assertEqual([shell["Cells"] for shell in plain["Corners"][0]["Shells"]], [2])


if __name__ == "__main__":
    unittest.main()
