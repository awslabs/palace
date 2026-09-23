#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import importlib.util
import struct
import unittest
from pathlib import Path

MODULE_PATH = Path(__file__).with_name("tag_island_conductor.py")
SPEC = importlib.util.spec_from_file_location("tag_island_conductor", MODULE_PATH)
TOOL = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TOOL)


def binary_msh2(nodes, triangles, names):
    """A binary MSH 2.2 file: nodes {tag: (x, y, z)}, triangles [(physical, elementary, (a, b, c))]."""
    out = bytearray(b"$MeshFormat\n2.2 1 8\n" + struct.pack("<i", 1) + b"\n$EndMeshFormat\n")
    out += f"$PhysicalNames\n{len(names)}\n".encode()
    for dimension, tag, name in names:
        out += f'{dimension} {tag} "{name}"\n'.encode()
    out += b"$EndPhysicalNames\n"
    out += f"$Nodes\n{len(nodes)}\n".encode()
    for tag, (x, y, z) in nodes.items():
        out += struct.pack("<iddd", tag, x, y, z)
    out += b"\n$EndNodes\n"
    out += f"$Elements\n{len(triangles)}\n".encode()
    out += struct.pack("<3i", TOOL.TRIANGLE, len(triangles), 2)
    for number, (physical, elementary, (a, b, c)) in enumerate(triangles, start=1):
        out += struct.pack("<6i", number, physical, elementary, a, b, c)
    out += b"\n$EndElements\n"
    return bytes(out)


class TagIslandConductorTest(unittest.TestCase):
    def setUp(self):
        # Ground: quad 1-2-3-4 (two triangles); port 6 triangle 3-4-5 touches the ground; feedline
        # quad 5-6-7-8 touches the port only; island: quad 9-10-11-12, isolated.
        self.nodes = {i: (float(i), 0.0, 0.0) for i in range(1, 13)}
        self.triangles = [
            (5, 1, (1, 2, 3)), (5, 1, (1, 3, 4)),
            (6, 2, (3, 4, 5)),
            (5, 3, (4, 5, 6)), (5, 3, (5, 7, 6)),
            (5, 4, (9, 10, 11)), (5, 4, (9, 11, 12)),
        ]
        self.names = [(2, 5, "metal"), (2, 6, "port_1"), (3, 1, "substrate")]

    def test_island_is_the_component_touching_no_port(self):
        data = binary_msh2(self.nodes, self.triangles, self.names)
        result, counts = TOOL.tag_island(data, metal=5, ground_adjacent={6}, island=9)
        self.assertEqual(counts, {"MetalTriangles": 6, "IslandTriangles": 2, "GroundTriangles": 4, "GroundComponents": 2})
        count, offset = TOOL.parse_sections(result)
        elements = TOOL.read_elements(result, count, offset)
        self.assertEqual([e[3] for e in elements], [5, 5, 6, 5, 5, 9, 9])
        self.assertEqual([e[4] for e in elements], [t[2] for t in self.triangles])
        self.assertIn(b'2 9 "island"\n$EndPhysicalNames', result)
        self.assertIn(b"$PhysicalNames\n4\n", result)
        self.assertEqual(result[result.index(b"$Nodes"):result.index(b"$Elements")],
                         data[data.index(b"$Nodes"):data.index(b"$Elements")])

    def test_two_isolated_components_stop(self):
        triangles = self.triangles + [(5, 5, (2, 8, 12))]
        data = binary_msh2(self.nodes, triangles, self.names)
        with self.assertRaises(SystemExit):
            TOOL.tag_island(data, metal=5, ground_adjacent={6}, island=9)

    def test_existing_island_attribute_stops(self):
        data = binary_msh2(self.nodes, self.triangles, self.names + [(2, 9, "taken")])
        with self.assertRaises(SystemExit):
            TOOL.tag_island(data, metal=5, ground_adjacent={6}, island=9)


if __name__ == "__main__":
    unittest.main()
