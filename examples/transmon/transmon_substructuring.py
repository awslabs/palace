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
the region and the environment follows the element faces. To redesign the region after the
offline run, label the conductors of the new design's mesh with --input/--output and build the
redesigned mesh with ../substructuring/remesh_region.py (see transmon_substructuring_redesign.json
and the substructuring documentation).
"""

import argparse
import collections
import os
import sys

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "substructuring")
)
from remesh_region import TETRAHEDRA, TRIANGLES, face_components, read_msh2, write_msh2


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
        if t in TRIANGLES and tags[0] == metal_attribute
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
        if etype in TETRAHEDRA:
            centroid = [sum(nodes[v][d] for v in nd[:4]) / 4.0 for d in range(3)]
            if not all(lo <= c <= hi for c, (lo, hi) in zip(centroid, box)):
                tags[0] = env_attribute[tags[0]]
            count[tags[0]] += 1
    if face_components(elements, env_attribute.keys()) > 1:
        sys.exit("the region of the box is not face-connected: change --box")

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
