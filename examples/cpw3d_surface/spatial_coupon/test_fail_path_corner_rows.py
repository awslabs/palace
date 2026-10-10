# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""F5T R-7 (supervisor decisions 674 (2) / 675): the audit producer's CornerNeighborhoods rows
of the fail-path legacy corners named by the build census - the order-invariant kappa_reg
RECOMPUTED over the corner-incident tetrahedra of the written mesh (never copied from the
census), Measure RegularCondition, FailPathMeasure true and the stored-order vertex-frame
reading; every other row byte-identical to the rows of before; a point both invariant and
fail-path fails closed."""
import math
import unittest

import meshio
import numpy as np

from general_mesh_audit_producer import CORNER_MEASURES, _point_aspects
from mesh_stage_contract import fail_path_legacy_corners


def corner_mesh():
    # Two tetrahedra sharing the corner point [0, 0, 0] (one right-angled, one skewed) and
    # a third one away from it (another corner at [10, 0, 0]).
    points = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
                         [0.3, 0.3, -0.05], [-0.4, 0.5, 0.2],
                         [10.0, 0.0, 0.0], [11.0, 0.0, 0.0], [10.0, 1.0, 0.0], [10.0, 0.0, 1.0]], dtype=float)
    cells = np.asarray([[0, 1, 2, 3], [4, 0, 5, 3], [6, 7, 8, 9]])
    return meshio.Mesh(points, [("tetra", cells)], cell_data={"gmsh:physical": [np.ones(len(cells), dtype=int)]})


W_INVERSE = np.linalg.inv(np.array([[1.0, 0.5, 0.5], [0.0, math.sqrt(3.0) / 2.0, math.sqrt(3.0) / 6.0],
                                    [0.0, 0.0, math.sqrt(2.0 / 3.0)]]))


def condition(matrix):
    singular = np.linalg.svd(matrix, compute_uv=False)
    return singular[0] / singular[-1]


def measures(points, cell):
    xyz = points[cell]
    jacobian = np.stack((xyz[1] - xyz[0], xyz[2] - xyz[0], xyz[3] - xyz[0]), axis=1)
    return condition(jacobian), condition(jacobian @ W_INVERSE)


class FailPathCornerRowsTest(unittest.TestCase):
    CORNERS = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]

    def test_plain_rows_are_unchanged_without_fail_path_corners(self):
        mesh = corner_mesh()
        rows = _point_aspects(mesh, self.CORNERS)
        self.assertEqual(rows, _point_aspects(mesh, self.CORNERS, None, None))
        self.assertEqual(rows, _point_aspects(mesh, self.CORNERS, [], []))
        frames = [measures(mesh.points, cell)[0] for cell in mesh.cells[0].data]
        self.assertEqual(rows, [{"Point": self.CORNERS[0], "MaximumAspect": max(frames[:2]), "Cells": 2},
                                {"Point": self.CORNERS[1], "MaximumAspect": frames[2], "Cells": 1}])
        for row in rows:
            self.assertNotIn("Measure", row)
            self.assertNotIn("FailPathMeasure", row)

    def test_fail_path_corner_row_recomputes_kappa_reg_and_keeps_the_frame_reading(self):
        mesh = corner_mesh()
        rows = _point_aspects(mesh, self.CORNERS, fail_path=[self.CORNERS[0]])
        frames = [measures(mesh.points, cell) for cell in mesh.cells[0].data]
        self.assertEqual(rows[0], {"Point": self.CORNERS[0],
                                   "MaximumAspect": max(frames[0][1], frames[1][1]),
                                   "Cells": 2, "Measure": CORNER_MEASURES["Invariant"],
                                   "MaximumVertexFrameAspect": max(frames[0][0], frames[1][0]),
                                   "FailPathMeasure": True})
        self.assertEqual(rows[0]["Measure"], "RegularCondition")
        # The other corner's row is the row of before, byte for byte.
        self.assertEqual(rows[1], _point_aspects(mesh, self.CORNERS)[1])
        # The reading is order-invariant: the skewed cell relabelled gives the same kappa_reg.
        relabelled = corner_mesh()
        relabelled.cells[0].data[1] = [0, 3, 4, 5]
        again = _point_aspects(relabelled, self.CORNERS, fail_path=[self.CORNERS[0]])
        self.assertAlmostEqual(again[0]["MaximumAspect"], rows[0]["MaximumAspect"], places=12)
        # ... while the stored-order frame reading changes with the relabelling (the lottery).
        self.assertNotAlmostEqual(again[0]["MaximumVertexFrameAspect"], rows[0]["MaximumVertexFrameAspect"], places=6)

    def test_invariant_rows_keep_their_form_and_a_corner_cannot_be_both(self):
        mesh = corner_mesh()
        rows = _point_aspects(mesh, self.CORNERS, invariant=[self.CORNERS[0]])
        self.assertEqual(rows[0]["Measure"], "RegularCondition")
        self.assertNotIn("FailPathMeasure", rows[0])
        self.assertEqual(rows, _point_aspects(mesh, self.CORNERS, [self.CORNERS[0]], []))
        with self.assertRaisesRegex(ValueError, "both invariant and a fail-path legacy corner"):
            _point_aspects(mesh, self.CORNERS, invariant=[self.CORNERS[0]], fail_path=[self.CORNERS[0]])

    def test_fail_path_corners_are_named_by_the_census_rows(self):
        optimization = {"CornerMeasures": [
            {"Point": [1.0, 2.0, 0.0], "Kind": "Legacy", "After": 3.9},
            {"Point": [3.0, 4.0, 0.0], "Kind": "Legacy", "After": 4.55,
             "FailPathMeasure": {"Measure": "RegularCondition", "After": 3.78, "Passed": True}},
            {"Point": [5.0, 6.0, 0.0], "Kind": "Invariant", "After": 3.2}]}
        self.assertEqual(fail_path_legacy_corners(optimization), [[3.0, 4.0, 0.0]])
        self.assertEqual(fail_path_legacy_corners({"CornerAspectsAfter": [3.0]}), [])
        self.assertEqual(fail_path_legacy_corners({}), [])
        self.assertEqual(fail_path_legacy_corners(None), [])


if __name__ == "__main__":
    unittest.main()
