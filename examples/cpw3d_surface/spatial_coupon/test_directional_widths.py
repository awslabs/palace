# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The band anisotropy statistics exclude a recorded edge layer's cells and report
them separately (supervisor decision 31); a layer covering the whole one-NormalSize
band makes the design gate not applicable by construction and records the
layer-adjacent band informationally (supervisor decision 35)."""
import unittest

import meshio
import numpy as np

from audit_edge_metric_mesh import (ANISOTROPY_GATE_APPLIED, ANISOTROPY_GATE_NOT_APPLICABLE,
                                    LAYER_ADJACENT_BAND_RULE, achieved_anisotropy,
                                    directional_widths)


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

    def test_achieved_anisotropy_record_and_layer_covered_band(self):
        mesh = band_mesh()
        recipe = {"PhysicalSegments": [[0., 0., 0., 3., 0., 0.]], "TruePhysicalCorners": [[-9., 0., 0.]],
                  "NormalSize": .025, "TangentialSize": .1}
        layer = {"EdgeSize": .001, "Aspect": 4.}
        # No layer: the first non-empty cutoff is the gated sample, as before.
        plain = achieved_anisotropy(directional_widths(mesh, recipe), .025)
        self.assertEqual((plain["Gate"], plain["Samples"], plain["ExcludedEdgeLayerCells"]),
                         (ANISOTROPY_GATE_APPLIED, 240, 0))       # the unrecorded layer rows
        self.assertEqual(plain["DistanceCutoff"], .025)
        self.assertLess(plain["TangentialP50"], .02)
        self.assertIsNone(plain["EdgeLayer"]); self.assertIsNone(plain["LayerAdjacentBand"])
        # With the layer recorded and band cells outside it within the gated cutoff, the
        # gate applies to those cells and the adjacent band is still recorded.
        spans = [[1., 0., 0., 1.4, 0., 0.]]
        widths = directional_widths(mesh, recipe, layer_spans=spans, layer_reach=.005)
        self.assertEqual(widths["1.0"]["Cells"], 0)                # band row centroids ~30 nm
        self.assertEqual(widths["1.0"]["ExcludedEdgeLayerCells"], 240)
        self.assertIsNone(widths["1.0"]["NearestSpanVertexDistance"])
        outer = widths["3.0"]["NearestSpanVertexDistance"]
        self.assertAlmostEqual(outer["Minimum"], .02 * np.sqrt(2))   # band row corner to the span
        self.assertAlmostEqual(outer["Maximum"], np.sqrt(.575**2 + 2 * .02**2))  # last box past the span end
        covered = achieved_anisotropy(widths, .025, layer)
        self.assertEqual(covered["Gate"], ANISOTROPY_GATE_NOT_APPLICABLE)
        self.assertEqual((covered["Samples"], covered["DistanceCutoff"], covered["ExcludedEdgeLayerCells"]),
                         (0, .025, 240))
        self.assertIsNone(covered["TangentialP50"])
        adjacent = covered["LayerAdjacentBand"]
        self.assertEqual((adjacent["Cells"], adjacent["DistanceCutoff"], adjacent["Rule"]),
                         (240, .025 * 3., LAYER_ADJACENT_BAND_RULE))
        self.assertAlmostEqual(adjacent["TangentialP50"], .025)
        self.assertEqual(adjacent["TransverseP90OverNormalSize"],
                         max(adjacent["Transverse1P90"], adjacent["Transverse2P90"]) / .025)
        self.assertEqual(adjacent["NearestSpanVertexDistance"], outer)
        self.assertEqual((covered["EdgeLayer"]["Cells"], covered["EdgeLayer"]["EdgeSize"],
                          covered["EdgeLayer"]["Aspect"]), (240, .001, 4.))
        self.assertAlmostEqual(covered["EdgeLayer"]["TangentialP50"], .01)
        # A layer that leaves band cells within one NormalSize keeps the gate applied to
        # them (here a wider band reaches the isotropic row), the layer still reported.
        partial = directional_widths(mesh, {**recipe, "NormalSize": .05}, layer_spans=spans,
                                     layer_reach=.005)
        applied = achieved_anisotropy(partial, .05, layer)
        self.assertEqual(applied["Gate"], ANISOTROPY_GATE_APPLIED)
        self.assertEqual((applied["Samples"], applied["DistanceCutoff"], applied["ExcludedEdgeLayerCells"]),
                         (partial["1.0"]["Cells"], .05, 240))
        self.assertGreater(applied["Samples"], 0)                  # the nearer tets of each box
        self.assertAlmostEqual(applied["TangentialP50"], .025)
        self.assertEqual(applied["LayerAdjacentBand"]["Cells"], 240)
        self.assertEqual(applied["EdgeLayer"]["Cells"], 240)
        with self.assertRaises(ValueError):
            achieved_anisotropy(directional_widths(
                mesh, {**recipe, "TruePhysicalCorners": [[1.5, 0., 0.]]}), .025)

    def test_mesh_quality_reports_layer_and_outside_regions(self):
        # Decision 32: the quality audit reports the recorded layer's cells (orientation,
        # edge aspect, scaled Jacobian as a diagnostic) and the cells outside it
        # separately; the whole-mesh statistics are unchanged.
        from edge_volume_metric import EDGE_LAYER_QUALITY_RULE
        from general_mesh_audit_producer import _tetra_quality
        mesh = band_mesh()
        # Orient every box tetrahedron positively (the fixture's split is unoriented).
        cells = mesh.cells[0].data
        xyz = mesh.points[cells]
        negative = np.linalg.det(np.stack((xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0],
                                           xyz[:, 3] - xyz[:, 0]), axis=2)) < 0
        cells[negative] = cells[negative][:, [0, 2, 1, 3]]
        plain = _tetra_quality(mesh)
        self.assertTrue(plain["PositiveOrientation"])
        self.assertIsNone(plain["EdgeLayer"]); self.assertNotIn("OutsideEdgeLayer", plain)
        quality = _tetra_quality(mesh, [[1., 0., 0., 1.4, 0., 0.]], .005)
        self.assertEqual({k: v for k, v in quality.items() if k not in ("EdgeLayer", "OutsideEdgeLayer")},
                         {k: v for k, v in plain.items() if k != "EdgeLayer"})
        layer, outside = quality["EdgeLayer"], quality["OutsideEdgeLayer"]
        self.assertEqual(layer["Cells"], 240); self.assertEqual(outside["Samples"], 240)
        self.assertEqual(layer["Cells"] + outside["Samples"], quality["Samples"])
        self.assertTrue(layer["PositiveOrientation"]); self.assertGreater(layer["MinimumDeterminant"], 0.)
        # 1 nm x 10 nm x 1 nm boxes split into tets: longest edge sqrt(102) nm over a
        # height between the box height (1 nm) and its diagonal split (1/sqrt(2) nm).
        diagonal = np.sqrt(.01**2 + 2 * .001**2)
        self.assertGreaterEqual(layer["MaximumEdgeAspect"], diagonal / .001 - 1e-9)
        self.assertLessEqual(layer["MaximumEdgeAspect"], diagonal / (.001 / np.sqrt(2)) + 1e-9)
        self.assertLess(layer["MinimumScaledJacobian"], .01)          # diagnostic
        self.assertEqual(sum(layer["CellsByScaledJacobianDecade"].values()), 240)
        self.assertEqual(layer["Reach"], .005); self.assertEqual(layer["Rule"], EDGE_LAYER_QUALITY_RULE)
        self.assertGreater(outside["MinimumScaledJacobian"], .1)       # isotropic 25 nm cubes


if __name__ == "__main__":
    unittest.main()
