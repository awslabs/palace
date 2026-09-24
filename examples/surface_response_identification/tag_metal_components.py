#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tag every floating metal component of a binary Gmsh MSH 2.2 mesh as its own boundary attribute.

The classifier gives every PEC attribute conductor 0 (surfaceresponseoperator.cpp GetConductor;
survey blocker 1), so a mesh whose whole metallisation is one physical surface cannot yield
DifferentConductorGap features. This tool generalises examples/transmon/tag_island_conductor.py
to N islands: the metal physical surface is split into edge-connected components (first- or
second-order triangles); every component which touches none of the ground-adjacent attributes
(the port patches, the exterior boundary, ...) is an island and is re-tagged as a new physical
surface FIRST, FIRST + 1, ... in the order of its smallest node id (deterministic), named
island_1, island_2, ...; everything else stays in the metal attribute. Nodes, element numbering
and elementary tags are copied byte for byte.

    python3 -m surface_response_identification.tag_metal_components --mesh M.msh2 \\
        --output M_islands.msh2 --metal 6 --ground-adjacent 3 7 8 [--first-island 9] \\
        [--drop-metal-duplicates]

--drop-metal-duplicates removes every metal triangle whose corner nodes coincide with a
triangle of another attribute (port sheets drawn on top of the metal: Palace rejects a face
with two boundary elements, geodata.cpp GetFaceToBdrElementMap); the element section is then
rewritten with consecutive element numbers.
"""

import argparse
import hashlib
import json
import struct
import sys
from collections import defaultdict
from pathlib import Path

NODES_PER_TYPE = {1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6, 10: 9, 11: 10, 15: 1}
TRIANGLE_TYPES = (2, 9)


def parse_sections(data):
    header_start = data.index(b"$MeshFormat\n") + len(b"$MeshFormat\n")
    header_end = data.index(b"\n", header_start)
    version, file_type, data_size = data[header_start:header_end].split()
    if version != b"2.2" or file_type != b"1" or data_size != b"8":
        raise SystemExit(f"expected a binary MSH 2.2 file with 8-byte floats, found {version} {file_type} {data_size}")
    elements_tag = data.index(b"$Elements\n")
    count_end = data.index(b"\n", elements_tag + len(b"$Elements\n"))
    element_count = int(data[elements_tag + len(b"$Elements\n"):count_end])
    return element_count, count_end + 1


def read_elements(data, element_count, offset):
    """(type, byte offset of the record, tag count, physical tag, node ids) per element."""
    elements = []
    position = offset
    while len(elements) < element_count:
        element_type, following, tag_count = struct.unpack_from("<3i", data, position)
        position += 12
        record = struct.Struct(f"<{1 + tag_count + NODES_PER_TYPE[element_type]}i")
        for _ in range(following):
            values = record.unpack_from(data, position)
            elements.append((element_type, position, tag_count, values[1] if tag_count else 0, values[1 + tag_count:]))
            position += record.size
    if data[position:position + len(b"\n$EndElements")] != b"\n$EndElements":
        raise SystemExit("the element block did not end at $EndElements")
    return elements


def triangle_edges(nodes):
    a, b, c = nodes[:3]
    return ((min(a, b), max(a, b)), (min(b, c), max(b, c)), (min(c, a), max(c, a)))


def edge_connected_components(triangles):
    """Components of triangles (index -> node tuple, corner nodes first) connected through
    shared corner edges, each sorted; the list is ordered by the smallest node id."""
    by_edge = defaultdict(list)
    for index, nodes in triangles.items():
        for edge in triangle_edges(nodes):
            by_edge[edge].append(index)
    component = {}
    components = []
    for start in sorted(triangles):
        if start in component:
            continue
        members, stack = [], [start]
        component[start] = len(components)
        while stack:
            index = stack.pop()
            members.append(index)
            for edge in triangle_edges(triangles[index]):
                for other in by_edge[edge]:
                    if other not in component:
                        component[other] = len(components)
                        stack.append(other)
        components.append(sorted(members))
    components.sort(key=lambda members: min(min(triangles[i][:3]) for i in members))
    return components


def drop_metal_duplicates(data, *, metal):
    """Rewrite $Elements without the metal triangles duplicating a face of another attribute."""
    element_count, offset = parse_sections(data)
    elements = read_elements(data, element_count, offset)
    faces = defaultdict(set)
    for element_type, _, _, physical, nodes in elements:
        if element_type in TRIANGLE_TYPES:
            faces[tuple(sorted(nodes[:3]))].add(physical)

    def duplicate(element):
        element_type, _, _, physical, nodes = element
        return element_type in TRIANGLE_TYPES and physical == metal and len(faces[tuple(sorted(nodes[:3]))]) > 1

    return rewrite_elements(data, elements, offset, duplicate)


def drop_attributes(data, attributes):
    """Rewrite $Elements without the surface elements of the given physical attributes (to
    mimic a production mesh without a substrate_air group)."""
    element_count, offset = parse_sections(data)
    elements = read_elements(data, element_count, offset)
    return rewrite_elements(data, elements, offset, lambda e: e[0] in TRIANGLE_TYPES and e[3] in attributes)


TETRAHEDRON_TYPES = (4, 11)
TETRAHEDRON_FACES = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def add_material_interface_group(data, attribute, name="substrate_air"):
    """Add first-order triangles with the given physical attribute on every face shared by
    two volume elements of different material attributes which carries no boundary element
    yet (the substrate / vacuum interface outside the metal: a production mesh without a
    substrate_air group). Second-order tetrahedra contribute their corner faces."""
    element_count, offset = parse_sections(data)
    elements = read_elements(data, element_count, offset)
    existing = {tuple(sorted(e[4][:3])) for e in elements if e[0] in TRIANGLE_TYPES}
    faces = defaultdict(set)
    for element_type, _, _, physical, nodes in elements:
        if element_type in TETRAHEDRON_TYPES:
            for face in TETRAHEDRON_FACES:
                faces[tuple(sorted(nodes[i] for i in face))].add(physical)
    new_faces = sorted(face for face, materials in faces.items() if len(materials) > 1 and face not in existing)
    if not new_faces:
        return data, 0
    names_tag = data.index(b"$PhysicalNames\n")
    count_end = data.index(b"\n", names_tag + len(b"$PhysicalNames\n"))
    names_end = data.index(b"$EndPhysicalNames")
    if f"2 {attribute} ".encode() in data[names_tag:names_end]:
        raise SystemExit(f"physical surface {attribute} already exists")
    count = int(data[names_tag + len(b"$PhysicalNames\n"):count_end])
    names_block = f"{count + 1}\n".encode() + data[count_end + 1:names_end] + f'2 {attribute} "{name}"\n'.encode()
    elements_tag = data.index(b"$Elements\n")
    end = data.index(b"\n$EndElements", offset)
    payload = bytearray(data[offset:end])
    record = struct.Struct("<6i")
    payload += struct.pack("<3i", 2, len(new_faces), 2)
    for k, face in enumerate(new_faces):
        payload += record.pack(element_count + k + 1, attribute, attribute, *face)
    result = data[:names_tag + len(b"$PhysicalNames\n")] + names_block + data[names_end:elements_tag] + f"$Elements\n{element_count + len(new_faces)}\n".encode() + bytes(payload) + data[end:]
    return result, len(new_faces)


def rewrite_elements(data, elements, offset, drop):
    dropped = 0
    kept = []
    for element in elements:
        if drop(element):
            dropped += 1
            continue
        kept.append(element)
    if dropped == 0:
        return data, 0
    blocks = []
    for element in kept:
        if blocks and blocks[-1][0] == (element[0], element[2]):
            blocks[-1][1].append(element)
        else:
            blocks.append(((element[0], element[2]), [element]))
    payload = bytearray()
    number = 1
    for (element_type, tag_count), members in blocks:
        payload += struct.pack("<3i", element_type, len(members), tag_count)
        record = struct.Struct(f"<{1 + tag_count + NODES_PER_TYPE[element_type]}i")
        for _, position, _, _, _ in members:
            values = list(record.unpack_from(data, position))
            values[0] = number
            number += 1
            payload += record.pack(*values)
    elements_tag = data.index(b"$Elements\n")
    end = data.index(b"\n$EndElements", offset)
    result = data[:elements_tag] + f"$Elements\n{len(kept)}\n".encode() + bytes(payload) + data[end:]
    return result, dropped


def tag_components(data, *, metal, ground_adjacent, first_island):
    element_count, offset = parse_sections(data)
    elements = read_elements(data, element_count, offset)
    metal_triangles = {i: e[4] for i, e in enumerate(elements) if e[0] in TRIANGLE_TYPES and e[3] == metal}
    if not metal_triangles:
        raise SystemExit(f"no triangles carry the metal physical tag {metal}")
    ground_nodes = {node for e in elements if e[0] in TRIANGLE_TYPES and e[3] in ground_adjacent for node in e[4]}
    if not ground_nodes:
        raise SystemExit(f"no triangles carry the ground-adjacent tags {sorted(ground_adjacent)}")
    components = edge_connected_components(metal_triangles)
    islands = [members for members in components if not any(node in ground_nodes for index in members for node in metal_triangles[index])]
    names_tag = data.index(b"$PhysicalNames\n")
    count_end = data.index(b"\n", names_tag + len(b"$PhysicalNames\n"))
    names_end = data.index(b"$EndPhysicalNames")
    count = int(data[names_tag + len(b"$PhysicalNames\n"):count_end])
    mutable = bytearray(data)
    new_names = b""
    summary = []
    for k, members in enumerate(islands):
        attribute = first_island + k
        if f"2 {attribute} ".encode() in data[names_tag:names_end]:
            raise SystemExit(f"physical surface {attribute} already exists")
        for index in members:
            _, position, tag_count, _, _ = elements[index]
            if tag_count < 1:
                raise SystemExit("metal element without a physical tag")
            struct.pack_into("<i", mutable, position + 4, attribute)
        new_names += f'2 {attribute} "island_{k + 1}"\n'.encode()
        summary.append({"Attribute": attribute, "Triangles": len(members), "SmallestNode": min(min(metal_triangles[i][:3]) for i in members)})
    names_block = f"{count + len(islands)}\n".encode() + data[count_end + 1:names_end] + new_names
    result = bytes(mutable[:names_tag + len(b"$PhysicalNames\n")]) + names_block + bytes(mutable[names_end:])
    return result, {
        "MetalTriangles": len(metal_triangles),
        "Components": len(components),
        "GroundComponents": len(components) - len(islands),
        "Islands": summary,
        "IslandAttributes": [s["Attribute"] for s in summary],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", type=Path, required=True, help="binary MSH 2.2 input")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metal", type=int, required=True, help="physical surface of the whole metal sheet")
    parser.add_argument("--ground-adjacent", type=int, nargs="+", required=True, help="physical surfaces that touch the grounded components only (ports, exterior boundary)")
    parser.add_argument("--first-island", type=int, default=9, help="attribute of the first island; the k-th island gets first + k - 1")
    parser.add_argument("--drop-metal-duplicates", action="store_true", help="remove metal triangles coinciding with a triangle of another attribute")
    parser.add_argument("--drop-attributes", type=int, nargs="*", default=[], help="remove the surface elements of these attributes (e.g. the substrate_air group)")
    parser.add_argument("--add-interface-group", type=int, help="add triangles with this attribute on every material-interface face without a boundary element (a substrate_air group)")
    args = parser.parse_args(argv)
    data = args.mesh.read_bytes()
    dropped = 0
    if args.drop_metal_duplicates:
        data, dropped = drop_metal_duplicates(data, metal=args.metal)
    dropped_attributes = 0
    if args.drop_attributes:
        data, dropped_attributes = drop_attributes(data, set(args.drop_attributes))
    added_interface = 0
    if args.add_interface_group:
        data, added_interface = add_material_interface_group(data, args.add_interface_group)
    result, counts = tag_components(data, metal=args.metal, ground_adjacent=set(args.ground_adjacent), first_island=args.first_island)
    counts["DroppedMetalDuplicates"] = dropped
    counts["DroppedAttributeElements"] = dropped_attributes
    counts["AddedInterfaceElements"] = added_interface
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(result)
    counts["InputSha256"] = hashlib.sha256(args.mesh.read_bytes()).hexdigest()
    counts["OutputSha256"] = hashlib.sha256(result).hexdigest()
    counts["Output"] = str(args.output)
    print(json.dumps(counts))
    return 0


if __name__ == "__main__":
    sys.exit(main())
