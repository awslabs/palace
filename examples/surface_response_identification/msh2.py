# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Minimal Gmsh MSH 2.2 reader (ASCII and binary) for the identification audit.

Only what the audit needs is read: node coordinates, element type / physical tag / node
ids, and the physical names. Elements of every dimension are kept so that boundary faces
can be tested for exteriority against the volume elements.
"""

import struct
from dataclasses import dataclass, field

import numpy as np

NODES_PER_TYPE = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6, 10: 9, 11: 10, 15: 1}
# Corner-node count of the linear element with the same shape (high-order elements are
# reduced to their corner nodes).
CORNER_NODES = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 2, 9: 3, 10: 4, 11: 4, 15: 1}
ELEMENT_DIMENSION = {1: 1, 2: 2, 3: 2, 4: 3, 5: 3, 6: 3, 7: 3, 8: 1, 9: 2, 10: 2, 11: 3, 15: 0}
LINE, TRIANGLE, QUAD, TET = 1, 2, 3, 4
# First- and second-order variants of the simplex types.
TRIANGLE_TYPES = (2, 9)
TETRAHEDRON_TYPES = (4, 11)


@dataclass
class Msh2:
    """Nodes as an (N, 3) array indexed by `node_index[tag]`; elements grouped by type."""

    coordinates: np.ndarray
    node_index: dict
    # type -> (physical tags (M,), corner node tags (M, k))
    elements: dict = field(default_factory=dict)
    physical_names: dict = field(default_factory=dict)  # (dim, tag) -> name
    # type -> all node tags (M, n) incl. the high-order nodes (Gmsh order)
    all_nodes: dict = field(default_factory=dict)

    def corner_indices(self, element_type):
        """(M, k) array of 0-based node indices for the given element type."""
        _, tags = self.elements[element_type]
        lookup = np.vectorize(self.node_index.__getitem__, otypes=[np.int64])
        return lookup(tags) if tags.size else tags.reshape(0, CORNER_NODES[element_type])

    def node_indices(self, element_type):
        """(M, n) array of 0-based node indices incl. the high-order nodes (Gmsh order: the
        corners, then the mid-edge nodes 0-1, 1-2, 2-0 of a type-9 triangle)."""
        tags = self.all_nodes[element_type]
        lookup = np.vectorize(self.node_index.__getitem__, otypes=[np.int64])
        return lookup(tags) if tags.size else tags.reshape(0, NODES_PER_TYPE[element_type])

    def physical_tags(self, element_type):
        return self.elements[element_type][0]

    def has(self, element_type):
        return element_type in self.elements and self.elements[element_type][0].size > 0


def _section(data, name):
    start = data.index(f"${name}\n".encode()) + len(name) + 2
    end = data.index(f"$End{name}".encode(), start)
    return start, end


ASCII_CHUNK_BYTES = 64 << 20


def ascii_element_lines(data, start, end):
    """Yield (token width, int64 array (m, width)) groups of equal-width lines of an ASCII
    element section, in chunks of ASCII_CHUNK_BYTES (bytes or mmap input)."""
    position = start
    while position < end:
        stop = min(end, position + ASCII_CHUNK_BYTES)
        if stop < end:
            stop = data.find(b"\n", stop) + 1
        chunk = data[position:stop]
        position = stop
        lines = np.array(chunk.split(b"\n"))
        lines = lines[np.char.str_len(lines) > 0]
        widths = np.char.count(lines, b" ") + 1
        for width in np.unique(widths):
            subset = lines[widths == width]
            yield int(width), np.fromstring(b" ".join(subset), dtype=np.int64, sep=" ").reshape(len(subset), width)


def read_msh2(path):
    with open(path, "rb") as source:
        data = source.read()
    fmt_start, fmt_end = _section(data, "MeshFormat")
    version, file_type, data_size = data[fmt_start:fmt_end].split()[:3]
    if version != b"2.2":
        raise ValueError(f"{path}: expected MSH 2.2, found {version.decode()}")
    binary = file_type == b"1"
    if binary and data_size != b"8":
        raise ValueError(f"{path}: binary MSH 2.2 must use 8-byte floats")

    physical_names = {}
    if b"$PhysicalNames\n" in data:
        start, end = _section(data, "PhysicalNames")
        lines = data[start:end].decode().splitlines()
        for line in lines[1 : 1 + int(lines[0])]:
            dim, tag, name = line.split(maxsplit=2)
            physical_names[(int(dim), int(tag))] = name.strip().strip('"')

    start, end = _section(data, "Nodes")
    count_end = data.index(b"\n", start)
    node_count = int(data[start:count_end])
    if binary:
        record = np.dtype([("tag", "<i4"), ("xyz", "<f8", (3,))])
        nodes = np.frombuffer(data, dtype=record, count=node_count, offset=count_end + 1)
        tags = nodes["tag"].astype(np.int64)
        coordinates = np.array(nodes["xyz"], dtype=np.float64)
    else:
        rows = np.array(data[count_end + 1 : end].split(), dtype=np.float64).reshape(node_count, 4)
        tags = rows[:, 0].astype(np.int64)
        coordinates = rows[:, 1:4].copy()
    node_index = {int(tag): index for index, tag in enumerate(tags)}

    start, end = _section(data, "Elements")
    count_end = data.index(b"\n", start)
    element_count = int(data[start:count_end])
    per_type = {}
    if binary:
        position = count_end + 1
        seen = 0
        while seen < element_count:
            element_type, following, tag_count = struct.unpack_from("<3i", data, position)
            position += 12
            nodes_per_element = NODES_PER_TYPE[element_type]
            width = 1 + tag_count + nodes_per_element
            block = np.frombuffer(data, dtype="<i4", count=following * width, offset=position)
            block = block.reshape(following, width)
            position += following * width * 4
            seen += following
            physical = block[:, 1] if tag_count else np.zeros(following, dtype=np.int32)
            corners = block[:, 1 + tag_count : 1 + tag_count + CORNER_NODES[element_type]]
            all_nodes = block[:, 1 + tag_count : 1 + tag_count + nodes_per_element]
            per_type.setdefault(element_type, []).append((physical.astype(np.int64), corners.astype(np.int64), all_nodes.astype(np.int64)))
    else:
        # Vectorized: consecutive lines of equal token width are parsed together (a
        # chip-scale ASCII mesh has 10^7 element lines).
        for width, rows in ascii_element_lines(data, count_end + 1, end):
            for element_type in np.unique(rows[:, 1]):
                element_type = int(element_type)
                block = rows[rows[:, 1] == element_type]
                tag_count = int(block[0, 2])
                if not (block[:, 2] == tag_count).all() or width != 3 + tag_count + NODES_PER_TYPE[element_type]:
                    raise ValueError(f"{path}: inconsistent ASCII element lines of type {element_type} and width {width}")
                physical = block[:, 3] if tag_count else np.zeros(len(block), dtype=np.int64)
                corners = block[:, 3 + tag_count : 3 + tag_count + CORNER_NODES[element_type]]
                all_nodes = block[:, 3 + tag_count : 3 + tag_count + NODES_PER_TYPE[element_type]]
                per_type.setdefault(element_type, []).append((physical.copy(), corners.copy(), all_nodes.copy()))
    elements = {}
    all_nodes = {}
    for element_type, blocks in per_type.items():
        elements[element_type] = (
            np.concatenate([b[0] for b in blocks]),
            np.concatenate([b[1] for b in blocks]),
        )
        all_nodes[element_type] = np.concatenate([b[2] for b in blocks])
    return Msh2(coordinates=coordinates, node_index=node_index, elements=elements, physical_names=physical_names, all_nodes=all_nodes)
