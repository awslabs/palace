#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tag the transmon island as its own boundary attribute in a binary Gmsh MSH 2.2 mesh.

The checked-in thin transmon mesh (mesh/transmon_surface_p1.msh2) carries the whole metal
sheet as one physical surface (5 "metal"); the electrostatic correction config written by
prepare_surface_response_electrostatic.py needs the island as terminal attribute 9 and the
rest of the sheet, with the port patches 6 and 7, as ground. This tool splits the metal
physical surface into edge-connected components and re-tags the one component that touches
none of the ground-adjacent attributes (the port patches) as the island attribute; every
other component touches a port patch and stays ground (the ground plane, and the feedline
centre conductor that the port patches join to it). The lumped-element sheet (4) bridges the
island and the ground and is not part of the metal surface. Everything else in the file
(nodes, element numbering, elementary tags, other sections) is copied byte for byte.
"""

import argparse
import hashlib
import struct
from collections import defaultdict
from pathlib import Path

NODES_PER_TYPE = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 15: 1}
TRIANGLE = 2


def parse_sections(data):
    """Byte offsets of the $Elements payload and the $PhysicalNames block of a binary MSH 2.2 file."""
    header_end = data.index(b"\n", data.index(b"$MeshFormat\n") + len(b"$MeshFormat\n"))
    version, file_type, data_size = data[data.index(b"$MeshFormat\n") + len(b"$MeshFormat\n"):header_end].split()
    if version != b"2.2" or file_type != b"1" or data_size != b"8":
        raise SystemExit(f"expected a binary MSH 2.2 file with 8-byte floats, found {version} {file_type} {data_size}")
    elements_tag = data.index(b"$Elements\n")
    count_end = data.index(b"\n", elements_tag + len(b"$Elements\n"))
    element_count = int(data[elements_tag + len(b"$Elements\n"):count_end])
    return element_count, count_end + 1


def read_elements(data, element_count, offset):
    """(type, byte offset of the element record, physical tag, node ids) for every element."""
    elements = []
    position = offset
    while len(elements) < element_count:
        element_type, following, tag_count = struct.unpack_from("<3i", data, position)
        position += 12
        nodes_per_element = NODES_PER_TYPE[element_type]
        record = struct.Struct(f"<{1 + tag_count + nodes_per_element}i")
        for _ in range(following):
            values = record.unpack_from(data, position)
            elements.append((element_type, position, tag_count, values[1] if tag_count else 0,
                             values[1 + tag_count:]))
            position += record.size
    if data[position:position + len(b"\n$EndElements")] != b"\n$EndElements":
        raise SystemExit("the element block did not end at $EndElements")
    return elements


def edge_connected_components(triangles):
    """Components of triangles (index -> node triple) connected through shared edges."""
    by_edge = defaultdict(list)
    for index, nodes in triangles.items():
        for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
            by_edge[(min(a, b), max(a, b))].append(index)
    component = {}
    components = []
    for start in triangles:
        if start in component:
            continue
        members, stack = [], [start]
        component[start] = len(components)
        while stack:
            index = stack.pop()
            members.append(index)
            nodes = triangles[index]
            for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
                for other in by_edge[(min(a, b), max(a, b))]:
                    if other not in component:
                        component[other] = len(components)
                        stack.append(other)
        components.append(members)
    return components


def tag_island(data, *, metal, ground_adjacent, island):
    element_count, offset = parse_sections(data)
    elements = read_elements(data, element_count, offset)
    metal_triangles = {i: e[4] for i, e in enumerate(elements) if e[0] == TRIANGLE and e[3] == metal}
    if not metal_triangles:
        raise SystemExit(f"no triangles carry the metal physical tag {metal}")
    ground_nodes = {node for e in elements if e[0] == TRIANGLE and e[3] in ground_adjacent for node in e[4]}
    if not ground_nodes:
        raise SystemExit(f"no triangles carry the ground-adjacent tags {sorted(ground_adjacent)}")
    components = edge_connected_components(metal_triangles)
    touching = [any(node in ground_nodes for index in members for node in metal_triangles[index]) for members in components]
    if touching.count(False) != 1:
        raise SystemExit(f"metal surface {metal} has {len(components)} edge-connected components of which "
                         f"{touching.count(False)} touch none of the ground-adjacent attributes {sorted(ground_adjacent)}; "
                         "expected exactly one (the island)")
    island_members = components[touching.index(False)]
    mutable = bytearray(data)
    for index in island_members:
        _, position, tag_count, _, _ = elements[index]
        if tag_count < 1:
            raise SystemExit("island element without a physical tag")
        struct.pack_into("<i", mutable, position + 4, island)
    names_tag = data.index(b"$PhysicalNames\n")
    count_end = data.index(b"\n", names_tag + len(b"$PhysicalNames\n"))
    names_end = data.index(b"$EndPhysicalNames")
    count = int(data[names_tag + len(b"$PhysicalNames\n"):count_end])
    if f'2 {island} '.encode() in data[names_tag:names_end]:
        raise SystemExit(f"physical surface {island} already exists")
    new_names = (f"{count + 1}\n".encode() + data[count_end + 1:names_end] + f'2 {island} "island"\n'.encode())
    result = bytes(mutable[:names_tag + len(b"$PhysicalNames\n")]) + new_names + bytes(mutable[names_end:])
    return result, {"MetalTriangles": len(metal_triangles), "IslandTriangles": len(island_members),
                    "GroundTriangles": len(metal_triangles) - len(island_members),
                    "GroundComponents": len(components) - 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", type=Path, required=True, help="binary MSH 2.2 input")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metal", type=int, default=5, help="physical surface of the whole metal sheet")
    parser.add_argument("--ground-adjacent", type=int, nargs="+", default=[6, 7],
                        help="physical surfaces that touch the ground component only (the port patches)")
    parser.add_argument("--island", type=int, default=9, help="new physical surface of the island")
    args = parser.parse_args()
    data = args.mesh.read_bytes()
    result, counts = tag_island(data, metal=args.metal, ground_adjacent=set(args.ground_adjacent), island=args.island)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(result)
    print(f"{args.output} island {args.island}: {counts['IslandTriangles']} of {counts['MetalTriangles']} metal triangles "
          f"(ground keeps {counts['GroundTriangles']} in {counts['GroundComponents']} components joined by the port patches); "
          f"input sha256 {hashlib.sha256(data).hexdigest()}; "
          f"output sha256 {hashlib.sha256(result).hexdigest()}")


if __name__ == "__main__":
    main()
