# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Prepare the qubit lattice mesh (qubit_lattice.jl) for substructuring (see
qubit_lattice_*.json).

The DeviceLayout mesh has a single "metal" boundary and no region/environment split. This
script

  1. splits the metal into its conductors: 4 = ground plane, 4 + k = pad k (k = 1, 2, ...,
     ordered by rows of increasing y, then by increasing x);
  2. assigns the tetrahedra whose centroid lies in a box around the middle qubit to the
     region (substrate 1, vacuum 2) and the rest to the environment (substrate 10, vacuum
     11).

Usage (from this directory):

    python3 qubit_lattice_substructuring.py [--input mesh/qubit_lattice.msh2]
        [--output mesh/qubit_lattice_substructuring.msh2]
        [--box XMIN XMAX YMIN YMAX ZMIN ZMAX]
"""

import argparse
import collections
import os

from remesh_region import TETRAHEDRA, TRIANGLES, read_msh2, write_msh2

EXTERIOR, GROUND, PAD = 3, 4, 4  # pad k: PAD + k
REGION = {"substrate": 1, "vacuum": 2}
ENVIRONMENT = {"substrate": 10, "vacuum": 11}


def conductors(elements, metal):
    """Connected components (by shared corner nodes) of the metal triangles."""
    parent = {}

    def find(a):
        while parent.setdefault(a, a) != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    tri = [
        e for e, (t, tags, _) in enumerate(elements) if t in TRIANGLES and tags[0] == metal
    ]
    for e in tri:
        corners = elements[e][2][:3]
        for v in corners[1:]:
            parent[find(corners[0])] = find(v)
    groups = collections.defaultdict(list)
    for e in tri:
        groups[find(elements[e][2][0])].append(e)
    return list(groups.values())


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--input", default=os.path.join(here, "mesh", "qubit_lattice.msh2"))
    parser.add_argument(
        "--output", default=os.path.join(here, "mesh", "qubit_lattice_substructuring.msh2")
    )
    parser.add_argument(
        "--box",
        type=float,
        nargs=6,
        default=[-500.0, 500.0, -500.0, 500.0, -400.0, 400.0],
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
        help="region box, in mesh units (um); for re-meshing with remesh_region.py it must "
        "not reach the bottom of the substrate",
    )
    args = parser.parse_args()

    names, nodes, elements = read_msh2(args.input)
    attr = {name.split('"')[1]: int(name.split()[1]) for name in names}

    # Conductors: the ground plane (the largest) and the pads.
    parts = sorted(conductors(elements, attr["metal"]), key=len, reverse=True)
    ground, pads = parts[0], parts[1:]

    def center(part):
        c = [nodes[v] for e in part for v in elements[e][2][:3]]
        return tuple(0.5 * (min(p[d] for p in c) + max(p[d] for p in c)) for d in range(2))

    pads.sort(key=lambda p: (round(center(p)[1], 3), center(p)[0]))
    for e, (t, tags, _) in enumerate(elements):
        if t in TRIANGLES and tags[0] == attr["exterior_boundary"]:
            tags[0] = EXTERIOR
    for a, part in [(GROUND, ground)] + [(PAD + k, p) for k, p in enumerate(pads, 1)]:
        for e in part:
            elements[e][1][0] = a

    # Region: tetrahedra with centroid in the box; fixed attribute numbers for the configs.
    box = list(zip(args.box[0::2], args.box[1::2]))
    material = {attr["substrate"]: "substrate", attr["vacuum"]: "vacuum"}
    count = collections.Counter()
    for etype, tags, nd in elements:
        if etype in TETRAHEDRA:
            centroid = [sum(nodes[v][d] for v in nd[:4]) / 4.0 for d in range(3)]
            inside = all(lo <= c <= hi for c, (lo, hi) in zip(centroid, box))
            tags[0] = (REGION if inside else ENVIRONMENT)[material[tags[0]]]
            count[tags[0]] += 1

    names = ['2 %d "exterior_boundary"' % EXTERIOR, '2 %d "ground"' % GROUND]
    names += ['2 %d "pad_%d"' % (PAD + k, k) for k in range(1, len(pads) + 1)]
    names += ['3 %d "%s"' % (a, m) for m, a in REGION.items()]
    names += ['3 %d "%s_environment"' % (a, m) for m, a in ENVIRONMENT.items()]
    write_msh2(args.output, names, nodes, elements)
    region = count[1] + count[2]
    print(
        "Wrote %s: %d pads, region %d tets (%.1f%%), environment %d tets"
        % (args.output, len(pads), region, 100.0 * region / sum(count.values()),
           count[10] + count[11])
    )
    inside = [k for k, p in enumerate(pads, 1)
              if all(lo <= c <= hi for c, (lo, hi) in zip(center(p), box))]
    print("  pads in the region: %s" % inside)


if __name__ == "__main__":
    main()
