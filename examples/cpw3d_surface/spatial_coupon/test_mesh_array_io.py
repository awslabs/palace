# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""mesh_array_io.read_msh22_binary: the direct MSH 2.2 binary reader returns exactly the
meshio.Mesh meshio.read returns (decision 62 step 3, proposal 4), both for Gmsh's
one-header-per-element layout and meshio's one-header-per-block layout, and declines
(None -> meshio fallback) any other layout."""
from pathlib import Path
import struct
import sys
import tempfile
import unittest

import meshio
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mesh_array_io import read_mesh, read_msh22_binary  # noqa: E402


def mixed_mesh():
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.], [1., 1., 0.],
                       [0., 0., -1.], [1., 0., -1.], [0., 1., -1.], [2., 0., 0.], [2., 1., 0.], [1.5, .5, 1.]])
    cells = [("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)),
             ("quad", np.array([[1, 8, 9, 4]], dtype=np.int32)),
             ("tetra", np.array([[0, 1, 2, 3], [1, 4, 2, 3]], dtype=np.int32)),
             ("wedge", np.array([[5, 6, 7, 0, 1, 2]], dtype=np.int32)),
             ("pyramid", np.array([[1, 8, 9, 4, 10]], dtype=np.int32))]
    cell_data = {"gmsh:physical": [np.array([7, 7], dtype=np.int32), np.array([6001], dtype=np.int32),
                                   np.array([1, 2], dtype=np.int32), np.array([1], dtype=np.int32),
                                   np.array([2], dtype=np.int32)],
                 "gmsh:geometrical": [np.array([11, 12], dtype=np.int32), np.array([13], dtype=np.int32),
                                      np.array([21, 21], dtype=np.int32), np.array([22], dtype=np.int32),
                                      np.array([23], dtype=np.int32)]}
    field_data = {"matching_surface": np.array([7, 2]), "ma_shell_top_1_6001": np.array([6001, 2]),
                  "substrate": np.array([1, 3]), "vacuum": np.array([2, 3])}
    return meshio.Mesh(points, cells, cell_data=cell_data, field_data=field_data)


def per_element_headers(data):
    """Rewrite meshio's per-block $Elements headers (type, count, tags) into Gmsh's
    per-element headers (type, 1, tags)."""
    node_count = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 15: 1}
    start = data.index(b"$Elements\n")
    header_end = data.index(b"\n", start + len(b"$Elements\n"))
    count = int(data[start + len(b"$Elements\n"):header_end])
    position = header_end + 1
    out = bytearray(data[:position])
    read = 0
    while read < count:
        element_type, block, tags = struct.unpack_from("<3i", data, position)
        position += 12
        width = 4 * (1 + tags + node_count[element_type])
        for _ in range(block):
            out += struct.pack("<3i", element_type, 1, tags) + data[position:position + width]
            position += width
        read += block
    out += data[position:]
    return bytes(out)


def assert_same_mesh(test, left, right):
    test.assertTrue(np.array_equal(left.points, right.points))
    test.assertEqual(left.points.dtype, right.points.dtype)
    test.assertEqual([c.type for c in left.cells], [c.type for c in right.cells])
    for a, b in zip(left.cells, right.cells):
        test.assertTrue(np.array_equal(a.data, b.data))
        test.assertEqual(a.data.dtype, b.data.dtype)
    test.assertEqual(list(left.cell_data), list(right.cell_data))
    for key in left.cell_data:
        test.assertEqual(len(left.cell_data[key]), len(right.cell_data[key]))
        for a, b in zip(left.cell_data[key], right.cell_data[key]):
            test.assertTrue(np.array_equal(a, b))
            test.assertEqual(a.dtype, b.dtype)
    test.assertEqual(list(left.field_data), list(right.field_data))
    for key in left.field_data:
        test.assertTrue(np.array_equal(left.field_data[key], right.field_data[key]))
    test.assertEqual(dict(left.point_data), dict(right.point_data))


class DirectMsh22ReaderTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.block_layout = self.root / "block.msh"
        meshio.write(self.block_layout, mixed_mesh(), file_format="gmsh22", binary=True)
        self.element_layout = self.root / "element.msh"
        self.element_layout.write_bytes(per_element_headers(self.block_layout.read_bytes()))

    def tearDown(self):
        self.temporary.cleanup()

    def test_both_header_layouts_equal_meshio(self):
        for path in (self.block_layout, self.element_layout):
            direct = read_msh22_binary(path)
            self.assertIsNotNone(direct)
            assert_same_mesh(self, direct, meshio.read(path))
            assert_same_mesh(self, read_mesh(path), meshio.read(path))
        self.assertNotEqual(self.block_layout.read_bytes(), self.element_layout.read_bytes())

    def test_other_layouts_decline_and_fall_back(self):
        ascii_path = self.root / "ascii.msh"
        meshio.write(ascii_path, mixed_mesh(), file_format="gmsh22", binary=False)
        self.assertIsNone(read_msh22_binary(ascii_path))
        assert_same_mesh(self, read_mesh(ascii_path), meshio.read(ascii_path))
        data = self.element_layout.read_bytes()
        trailing = self.root / "trailing.msh"
        trailing.write_bytes(data + b"$NodeData\n0\n$EndNodeData\n")
        self.assertIsNone(read_msh22_binary(trailing))
        truncated = self.root / "truncated.msh"
        truncated.write_bytes(data[:-40])
        self.assertIsNone(read_msh22_binary(truncated))
        renumbered = bytearray(data)
        start = data.index(b"$Nodes\n")
        header_end = data.index(b"\n", start + len(b"$Nodes\n"))
        struct.pack_into("<i", renumbered, header_end + 1, 5)  # the first node tag is not 1
        (self.root / "renumbered.msh").write_bytes(bytes(renumbered))
        self.assertIsNone(read_msh22_binary(self.root / "renumbered.msh"))


if __name__ == "__main__":
    unittest.main()
