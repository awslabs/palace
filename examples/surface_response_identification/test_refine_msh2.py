#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the external uniform refinement: element counts, volume and area
conservation, positive orientation, perimeter length and corners unchanged."""

import os
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import perimeter as P, refine_msh2 as R  # noqa: E402
from surface_response_identification.msh2 import read_msh2  # noqa: E402
from surface_response_identification.test_audit import CONFIG, island_mesh, write_msh2  # noqa: E402


def tet_volumes(points, tets):
    a, b, c, d = (points[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


class RefineTest(unittest.TestCase):
    def test_counts_volume_orientation_and_perimeter_are_conserved(self):
        nodes, elements, names = island_mesh()
        with tempfile.TemporaryDirectory() as directory:
            source = os.path.join(directory, "m.msh2")
            refined = os.path.join(directory, "r.msh2")
            write_msh2(source, nodes, elements, names, True)
            counts = R.refine_file(source, refined)
            coarse = read_msh2(source)
            fine = read_msh2(refined)
            self.assertEqual(counts["Elements"][4], 8 * len(coarse.physical_tags(4)))
            self.assertEqual(counts["Elements"][2], 4 * len(coarse.physical_tags(2)))
            self.assertEqual(sorted(fine.physical_names.items()), sorted(coarse.physical_names.items()))
            coarse_volume = tet_volumes(coarse.coordinates, coarse.corner_indices(4))
            fine_volume = tet_volumes(fine.coordinates, fine.corner_indices(4))
            self.assertTrue(np.all(fine_volume > 0.0))
            self.assertAlmostEqual(np.abs(coarse_volume).sum(), fine_volume.sum())
            # Physical tags of the children follow the parents.
            self.assertEqual(sorted(set(fine.physical_tags(2))), sorted(set(coarse.physical_tags(2))))
            before = P.extract_perimeter(coarse, CONFIG)
            after = P.extract_perimeter(fine, CONFIG)
            for kind in ("PHYSICAL", "TRUNCATION", "NONPLANAR", "CROSS_LAYER", "NONMANIFOLD"):
                self.assertAlmostEqual(before.length(kind), after.length(kind), msg=kind)
            self.assertEqual(len([e for e in after.edges if e.kind == "PHYSICAL"]), 2 * len([e for e in before.edges if e.kind == "PHYSICAL"]))
            self.assertEqual(sum(1 for v in after.vertices if v.physical_kind == "CORNER"), sum(1 for v in before.vertices if v.physical_kind == "CORNER"))
            self.assertEqual(after.chains, before.chains)


if __name__ == "__main__":
    unittest.main()
