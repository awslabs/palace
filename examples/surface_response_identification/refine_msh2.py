#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Uniformly refine a first-order MSH 2.2 mesh outside Palace (same geometry, every element
split once): tetrahedra 1 -> 8 (corner tetrahedra plus the inner octahedron cut along one
diagonal), triangles 1 -> 4, lines 1 -> 2, with the edge midpoints as new nodes. The output
is a binary MSH 2.2 file with the same physical tags, so the identification of the refined
mesh can be audited against the refined mesh itself (invariant A5) instead of relying on
Palace's internal UniformLevels, which leaves no mesh file for the audit.

    python3 -m surface_response_identification.refine_msh2 --mesh M.msh2 --output M_r1.msh2 [--levels 1]
"""

import argparse
import struct
import sys
import os

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from .msh2 import LINE, TET, TRIANGLE, read_msh2  # noqa: E402

POINT = 15


def refine_once(coordinates, elements):
    """coordinates: (n, 3); elements: {type: (physical (m,), nodes (m, k) 0-based)} of
    first-order simplices -> refined (coordinates, elements) in the same format."""
    supported = {LINE, TRIANGLE, TET, POINT}
    unsupported = set(elements) - supported
    if unsupported:
        raise NotImplementedError(f"only first-order simplices are refined (found element types {sorted(unsupported)})")
    points = [row for row in coordinates]
    midpoints = {}

    def midpoint(a, b):
        key = (a, b) if a < b else (b, a)
        index = midpoints.get(key)
        if index is None:
            index = len(points)
            points.append(0.5 * (coordinates[a] + coordinates[b]))
            midpoints[key] = index
        return index

    result = {}
    if POINT in elements:
        result[POINT] = elements[POINT]
    if LINE in elements:
        physical, nodes = elements[LINE]
        out_physical, out_nodes = [], []
        for tag, (a, b) in zip(physical, nodes):
            m = midpoint(a, b)
            out_nodes += [(a, m), (m, b)]
            out_physical += [tag, tag]
        result[LINE] = (np.array(out_physical), np.array(out_nodes))
    if TRIANGLE in elements:
        physical, nodes = elements[TRIANGLE]
        out_physical, out_nodes = [], []
        for tag, (a, b, c) in zip(physical, nodes):
            ab, bc, ca = midpoint(a, b), midpoint(b, c), midpoint(c, a)
            out_nodes += [(a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca)]
            out_physical += [tag] * 4
        result[TRIANGLE] = (np.array(out_physical), np.array(out_nodes))
    if TET in elements:
        physical, nodes = elements[TET]
        out_physical, out_nodes = [], []
        for tag, (a, b, c, d) in zip(physical, nodes):
            ab, ac, ad = midpoint(a, b), midpoint(a, c), midpoint(a, d)
            bc, bd, cd = midpoint(b, c), midpoint(b, d), midpoint(c, d)
            children = [
                (a, ab, ac, ad), (ab, b, bc, bd), (ac, bc, c, cd), (ad, bd, cd, d),
                # Inner octahedron cut along the diagonal ab - cd.
                (ab, cd, ac, ad), (ab, cd, ad, bd), (ab, cd, bd, bc), (ab, cd, bc, ac),
            ]
            out_nodes += children
            out_physical += [tag] * 8
        out_nodes = np.array(out_nodes)
        # Positive orientation of every child (the parent's orientation is not assumed).
        p = np.array(points)
        v = np.einsum("ij,ij->i", np.cross(p[out_nodes[:, 1]] - p[out_nodes[:, 0]], p[out_nodes[:, 2]] - p[out_nodes[:, 0]]), p[out_nodes[:, 3]] - p[out_nodes[:, 0]])
        flip = v < 0
        out_nodes[flip, 2], out_nodes[flip, 3] = out_nodes[flip, 3].copy(), out_nodes[flip, 2].copy()
        result[TET] = (np.array(out_physical), out_nodes)
    return np.array(points), result


def write_binary_msh2(path, coordinates, elements, physical_names):
    out = bytearray(b"$MeshFormat\n2.2 1 8\n" + struct.pack("<i", 1) + b"\n$EndMeshFormat\n")
    if physical_names:
        out += f"$PhysicalNames\n{len(physical_names)}\n".encode()
        for (dimension, tag), name in sorted(physical_names.items()):
            out += f'{dimension} {tag} "{name}"\n'.encode()
        out += b"$EndPhysicalNames\n"
    out += f"$Nodes\n{len(coordinates)}\n".encode()
    tags = np.arange(1, len(coordinates) + 1, dtype=np.int32)
    node_block = np.zeros(len(coordinates), dtype=[("tag", "<i4"), ("xyz", "<f8", 3)])
    node_block["tag"] = tags
    node_block["xyz"] = coordinates
    out += node_block.tobytes() + b"\n$EndNodes\n"
    total = sum(len(v[0]) for v in elements.values())
    out += f"$Elements\n{total}\n".encode()
    number = 1
    for element_type, (physical, nodes) in sorted(elements.items()):
        nodes = np.asarray(nodes)
        count = len(physical)
        out += struct.pack("<3i", int(element_type), count, 2)
        block = np.empty((count, 3 + nodes.shape[1]), dtype="<i4")
        block[:, 0] = np.arange(number, number + count)
        block[:, 1] = physical
        block[:, 2] = physical
        block[:, 3:] = nodes + 1
        out += block.tobytes()
        number += count
    out += b"\n$EndElements\n"
    with open(path, "wb") as target:
        target.write(bytes(out))


def refine_file(mesh_path, output_path, levels=1):
    mesh = read_msh2(mesh_path)
    coordinates = mesh.coordinates
    elements = {t: (np.asarray(mesh.physical_tags(t)), mesh.corner_indices(t)) for t in mesh.elements}
    for _ in range(levels):
        coordinates, elements = refine_once(coordinates, elements)
    write_binary_msh2(output_path, coordinates, elements, mesh.physical_names)
    return {"Nodes": int(len(coordinates)), "Elements": {int(t): int(len(v[0])) for t, v in elements.items()}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--levels", type=int, default=1)
    args = parser.parse_args(argv)
    print(refine_file(args.mesh, args.output, args.levels))
    return 0


if __name__ == "__main__":
    sys.exit(main())
