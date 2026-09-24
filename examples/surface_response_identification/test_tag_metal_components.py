#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests of the N-component metal tagger on a synthetic binary MSH 2.2 file."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from surface_response_identification import tag_metal_components as T  # noqa: E402
from surface_response_identification.msh2 import read_msh2  # noqa: E402
from surface_response_identification.test_audit import write_msh2  # noqa: E402


class TagMetalComponentsTest(unittest.TestCase):
    def test_two_floating_components_become_islands_in_node_order(self):
        # Metal (6): a grounded strip touching the port (7), and two floating squares.
        nodes = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),  # 1-4 grounded strip
                 (2, 0, 0), (2, 1, 0),  # 5-6 port
                 (5, 0, 0), (6, 0, 0), (6, 1, 0), (5, 1, 0),  # 7-10 island A
                 (8, 0, 0), (9, 0, 0), (9, 1, 0), (8, 1, 0),  # 11-14 island B
                 (4, 4, -1)]  # 15 tet apex
        elements = [
            (2, 6, (1, 2, 3)), (2, 6, (1, 3, 4)),
            (2, 7, (2, 5, 6)), (2, 7, (2, 6, 3)),
            (2, 6, (11, 12, 13)), (2, 6, (11, 13, 14)),  # island B listed first
            (2, 6, (7, 8, 9)), (2, 6, (7, 9, 10)),
            (4, 1, (1, 2, 3, 15)),
        ]
        names = [(3, 1, "substrate"), (2, 6, "metal"), (2, 7, "wave_port_1")]
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "m.msh2")
            write_msh2(path, nodes, elements, names, binary=True)
            with open(path, "rb") as source:
                data = source.read()
            result, counts = T.tag_components(data, metal=6, ground_adjacent={7}, first_island=9)
            self.assertEqual(counts["Components"], 3)
            self.assertEqual(counts["GroundComponents"], 1)
            self.assertEqual(counts["IslandAttributes"], [9, 10])
            # Island order follows the smallest node id, not the element order.
            self.assertEqual([s["SmallestNode"] for s in counts["Islands"]], [7, 11])
            output = os.path.join(directory, "out.msh2")
            with open(output, "wb") as target:
                target.write(result)
            mesh = read_msh2(output)
            self.assertEqual(mesh.physical_names[(2, 9)], "island_1")
            self.assertEqual(mesh.physical_names[(2, 10)], "island_2")
            tags = list(mesh.physical_tags(2))
            self.assertEqual(tags.count(6), 2)
            self.assertEqual(tags.count(9), 2)
            self.assertEqual(tags.count(10), 2)
            corners = mesh.corner_indices(2)
            island_1_nodes = {int(n) for row, tag in zip(corners, tags) if tag == 9 for n in row}
            self.assertEqual(island_1_nodes, {6, 7, 8, 9})  # 0-based indices of nodes 7-10
            # Everything but the physical names and the re-tagged records is byte-identical.
            self.assertEqual(len(result), len(data) + len(b'2 9 "island_1"\n2 10 "island_2"\n'))

    def test_no_floating_component_leaves_the_file_unchanged_apart_from_the_count(self):
        nodes = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (2, 0, 0), (2, 1, 0), (4, 4, -1)]
        elements = [(2, 6, (1, 2, 3)), (2, 6, (1, 3, 4)), (2, 7, (2, 5, 6)), (2, 7, (2, 6, 3)), (4, 1, (1, 2, 3, 7))]
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "m.msh2")
            write_msh2(path, nodes, elements, [(3, 1, "substrate"), (2, 6, "metal"), (2, 7, "port")], binary=True)
            with open(path, "rb") as source:
                data = source.read()
            result, counts = T.tag_components(data, metal=6, ground_adjacent={7}, first_island=9)
            self.assertEqual(counts["IslandAttributes"], [])
            self.assertEqual(result, data)

    def test_metal_triangles_duplicating_a_port_face_are_dropped(self):
        nodes = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (2, 0, 0), (2, 1, 0), (4, 4, -1)]
        elements = [(2, 6, (1, 2, 3)), (2, 6, (1, 3, 4)), (2, 6, (2, 5, 6)), (2, 7, (2, 5, 6)), (2, 7, (2, 6, 3)), (4, 1, (1, 2, 3, 7))]
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "m.msh2")
            write_msh2(path, nodes, elements, [(3, 1, "substrate"), (2, 6, "metal"), (2, 7, "port")], binary=True)
            with open(path, "rb") as source:
                data = source.read()
            result, dropped = T.drop_metal_duplicates(data, metal=6)
            self.assertEqual(dropped, 1)
            output = os.path.join(directory, "out.msh2")
            with open(output, "wb") as target:
                target.write(result)
            mesh = read_msh2(output)
            self.assertEqual(sorted(mesh.physical_tags(2)), [6, 6, 7, 7])
            self.assertEqual(len(mesh.physical_tags(4)), 1)
            # Re-reading and tagging the repaired file works (consecutive element numbers).
            again, counts = T.tag_components(result, metal=6, ground_adjacent={7}, first_island=9)
            self.assertEqual(counts["Components"], 1)
            self.assertEqual(again, result)

    def test_drop_attributes_removes_the_substrate_air_faces(self):
        nodes = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (2, 0, 0), (2, 1, 0), (4, 4, -1)]
        elements = [(2, 6, (1, 2, 3)), (2, 6, (1, 3, 4)), (2, 8, (2, 5, 6)), (2, 8, (2, 6, 3)), (4, 1, (1, 2, 3, 7))]
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "m.msh2")
            write_msh2(path, nodes, elements, [(3, 1, "substrate"), (2, 6, "metal"), (2, 8, "substrate_air")], binary=True)
            with open(path, "rb") as source:
                data = source.read()
            result, dropped = T.drop_attributes(data, {8})
            self.assertEqual(dropped, 2)
            output = os.path.join(directory, "out.msh2")
            with open(output, "wb") as target:
                target.write(result)
            mesh = read_msh2(output)
            self.assertEqual(sorted(mesh.physical_tags(2)), [6, 6])

    def test_add_interface_group_covers_the_uncovered_material_interface_faces(self):
        # Two tets (substrate 1 below, vacuum 2 above) sharing the face 1-2-3 at z = 0 which
        # has no boundary element; the metal triangle 1-2-4 is a different face.
        nodes = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, -1), (0, 0, 1)]
        elements = [(2, 6, (1, 2, 4)), (4, 1, (1, 2, 3, 5)), (4, 2, (1, 2, 3, 6))]
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, "m.msh2")
            write_msh2(path, nodes, elements, [(3, 1, "substrate"), (3, 2, "vacuum"), (2, 6, "metal")], binary=True)
            with open(path, "rb") as source:
                data = source.read()
            result, added = T.add_material_interface_group(data, 8)
            self.assertEqual(added, 1)
            output = os.path.join(directory, "out.msh2")
            with open(output, "wb") as target:
                target.write(result)
            mesh = read_msh2(output)
            self.assertEqual(mesh.physical_names[(2, 8)], "substrate_air")
            tags = list(mesh.physical_tags(2))
            self.assertEqual(sorted(tags), [6, 8])
            corners = mesh.corner_indices(2)
            self.assertEqual(sorted(int(n) for n in corners[tags.index(8)]), [0, 1, 2])
            again, added_again = T.add_material_interface_group(result, 9)
            self.assertEqual(added_again, 0)


if __name__ == "__main__":
    unittest.main()
