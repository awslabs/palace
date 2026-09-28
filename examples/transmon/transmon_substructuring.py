# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Prepare the transmon mesh for substructuring (see transmon_substructuring_*.json).

Substructuring needs the mesh split into a region (solved every run) and an environment
(condensed once), and a capacitance matrix needs one boundary attribute per conductor. The
DeviceLayout-generated mesh (mesh/transmon.msh2) has neither, so this script:

  1. splits the single "metal" boundary attribute (5) into its connected conductors:
       5 = ground plane (with the readout resonator, which is shorted to ground),
       8 = feedline center trace, 9 = transmon island;
  2. assigns the tetrahedra whose centroid lies in a box around the qubit to the region
     (substrate 1, vacuum 2) and the rest to the environment (substrate 10, vacuum 11).

Usage (from this directory):

    python3 transmon_substructuring.py [--box XMIN XMAX YMIN YMAX ZMIN ZMAX]

writes mesh/transmon_substructuring.msh2. Any conforming split works: the interface between
the region and the environment follows the element faces.
"""

import argparse
import collections
import os
import struct

# Nodes per Gmsh element type (the types that can appear in this mesh).
NODES_PER_ELEMENT = {1: 2, 2: 3, 4: 4, 8: 3, 9: 6, 11: 10, 15: 1}
TRIANGLE_TYPES = (2, 9)
TETRAHEDRON_TYPES = (4, 11)


def read_msh2(path):
    """Read a Gmsh 2.2 mesh (ASCII or binary) into (physical names, nodes, elements)."""
    data = open(path, "rb").read()
    pos = 0

    def line():
        nonlocal pos
        end = data.index(b"\n", pos)
        text = data[pos:end].decode()
        pos = end + 1
        return text

    names, nodes, elements, binary = [], {}, [], False
    while pos < len(data):
        section = line().strip()
        if section == "$MeshFormat":
            binary = int(line().split()[1]) == 1
            if binary:
                pos += 4  # endianness check integer
                line()
            line()
        elif section == "$PhysicalNames":
            names = [line() for _ in range(int(line()))]
            line()
        elif section == "$Nodes":
            count = int(line())
            for _ in range(count):
                if binary:
                    i, x, y, z = struct.unpack_from("<iddd", data, pos)
                    pos += 28
                else:
                    tok = line().split()
                    i, x, y, z = int(tok[0]), *map(float, tok[1:4])
                nodes[i] = (x, y, z)
            if binary:
                line()
            line()
        elif section == "$Elements":
            count, read = int(line()), 0
            while read < count:
                if binary:
                    etype, n, ntags = struct.unpack_from("<iii", data, pos)
                    pos += 12
                    k = NODES_PER_ELEMENT[etype]
                    for _ in range(n):
                        rec = struct.unpack_from("<%di" % (1 + ntags + k), data, pos)
                        pos += 4 * (1 + ntags + k)
                        elements.append([etype, list(rec[1 : 1 + ntags]), list(rec[1 + ntags :])])
                    read += n
                else:
                    tok = list(map(int, line().split()))
                    etype, ntags = tok[1], tok[2]
                    elements.append([etype, tok[3 : 3 + ntags], tok[3 + ntags :]])
                    read += 1
            if binary:
                line()
            line()
        elif section:
            while line().strip() != "$End" + section[1:]:
                pass
    return names, nodes, elements


def write_msh2(path, names, nodes, elements):
    """Write a binary Gmsh 2.2 mesh."""
    with open(path, "wb") as f:
        f.write(b"$MeshFormat\n2.2 1 8\n")
        f.write(struct.pack("<i", 1))
        f.write(b"\n$EndMeshFormat\n$PhysicalNames\n%d\n" % len(names))
        for name in names:
            f.write(name.encode() + b"\n")
        f.write(b"$EndPhysicalNames\n$Nodes\n%d\n" % len(nodes))
        for i, (x, y, z) in nodes.items():
            f.write(struct.pack("<iddd", i, x, y, z))
        f.write(b"\n$EndNodes\n$Elements\n%d\n" % len(elements))
        for num, (etype, tags, nd) in enumerate(elements, 1):
            f.write(struct.pack("<iii", etype, 1, len(tags)))
            f.write(struct.pack("<%di" % (1 + len(tags) + len(nd)), num, *tags, *nd))
        f.write(b"\n$EndElements\n")


def split_conductors(elements, metal_attribute):
    """Connected components (by shared corner nodes) of the triangles with metal_attribute,
    largest first."""
    parent = {}

    def find(a):
        while parent.setdefault(a, a) != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    metal = [
        e
        for e, (t, tags, _) in enumerate(elements)
        if t in TRIANGLE_TYPES and tags[0] == metal_attribute
    ]
    for e in metal:
        corners = elements[e][2][:3]
        for v in corners[1:]:
            parent[find(corners[0])] = find(v)
    components = collections.defaultdict(list)
    for e in metal:
        components[find(elements[e][2][0])].append(e)
    return sorted(components.values(), key=len, reverse=True)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--input", default=os.path.join(here, "mesh", "transmon.msh2"))
    parser.add_argument(
        "--output", default=os.path.join(here, "mesh", "transmon_substructuring.msh2")
    )
    parser.add_argument(
        "--box",
        type=float,
        nargs=6,
        default=[-300.0, 300.0, -200.0, 800.0, -300.0, 300.0],
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
        help="region box around the qubit, in mesh units (um)",
    )
    args = parser.parse_args()

    names, nodes, elements = read_msh2(args.input)

    ground, feedline, island = split_conductors(elements, metal_attribute=5)
    for e in feedline:
        elements[e][1][0] = 8
    for e in island:
        elements[e][1][0] = 9

    box = list(zip(args.box[0::2], args.box[1::2]))
    env_attribute = {1: 10, 2: 11}  # substrate, vacuum
    count = collections.Counter()
    for etype, tags, nd in elements:
        if etype in TETRAHEDRON_TYPES:
            centroid = [sum(nodes[v][d] for v in nd[:4]) / 4.0 for d in range(3)]
            if not all(lo <= c <= hi for c, (lo, hi) in zip(centroid, box)):
                tags[0] = env_attribute[tags[0]]
            count[tags[0]] += 1

    names = [name for name in names if not name.startswith("2 5 ")] + [
        '2 5 "ground"',
        '2 8 "feedline"',
        '2 9 "island"',
        '3 10 "substrate_environment"',
        '3 11 "vacuum_environment"',
    ]
    write_msh2(args.output, names, nodes, elements)
    print(
        "Wrote %s: region %d tets (substrate %d, vacuum %d), environment %d tets"
        % (
            args.output,
            count[1] + count[2],
            count[1],
            count[2],
            count[10] + count[11],
        )
    )


if __name__ == "__main__":
    main()
