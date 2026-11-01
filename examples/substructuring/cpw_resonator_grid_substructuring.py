# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Prepare the CPW resonator grid mesh (cpw_resonator_grid.jl) for substructuring (see
cpw_resonator_grid*.json).

The DeviceLayout mesh has no region/environment split. This script assigns the tetrahedra
whose centroid lies in a box around one resonator (by default resonator 6, in the middle row)
to the region (substrate 1, vacuum 2) and the rest to the environment (substrate 10, vacuum
11). The boundaries keep their attributes: 3 = exterior boundary, 4 = metal, 5, 6, ... = the
lumped ports. The region tetrahedra also get their own geometric (elementary) volumes, so that
viewers that show the mesh by geometric entity, like Gmsh, display the region apart from the
environment.

Usage (from this directory):

    python3 cpw_resonator_grid_substructuring.py [--input mesh/cpw_resonator_grid.msh2]
        [--output mesh/cpw_resonator_grid_substructuring.msh2]
        [--box XMIN XMAX YMIN YMAX ZMIN ZMAX]
"""

import argparse
import collections
import os

from remesh_region import TETRAHEDRA, read_msh2, write_msh2

REGION = {"substrate": 1, "vacuum": 2}
ENVIRONMENT = {"substrate": 10, "vacuum": 11}


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--input", default=os.path.join(here, "mesh", "cpw_resonator_grid.msh2")
    )
    parser.add_argument(
        "--output",
        default=os.path.join(here, "mesh", "cpw_resonator_grid_substructuring.msh2"),
    )
    parser.add_argument(
        "--box",
        type=float,
        nargs=6,
        default=[-900.0, 0.0, 650.0, 1400.0, -250.0, 250.0],
        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX"),
        help="region box, in mesh units (um)",
    )
    args = parser.parse_args()

    names, nodes, elements = read_msh2(args.input)
    attr = {name.split('"')[1]: int(name.split()[1]) for name in names}
    material = {attr["substrate"]: "substrate", attr["vacuum"]: "vacuum"}
    box = list(zip(args.box[0::2], args.box[1::2]))
    offset = 1 + max(tags[1] for etype, tags, _ in elements if etype in TETRAHEDRA)
    count = collections.Counter()
    for etype, tags, nd in elements:
        if etype in TETRAHEDRA:
            centroid = [sum(nodes[v][d] for v in nd[:4]) / 4.0 for d in range(3)]
            inside = all(lo <= c <= hi for c, (lo, hi) in zip(centroid, box))
            tags[0] = (REGION if inside else ENVIRONMENT)[material[tags[0]]]
            if inside:
                tags[1] += offset
            count[tags[0]] += 1

    names = [n for n in names if n.startswith("2 ")]
    names += ['3 %d "%s"' % (a, m) for m, a in REGION.items()]
    names += ['3 %d "%s_environment"' % (a, m) for m, a in ENVIRONMENT.items()]
    write_msh2(args.output, names, nodes, elements)
    print(
        "Wrote %s: region %d tets, environment %d tets"
        % (args.output, count[1] + count[2], count[10] + count[11])
    )


if __name__ == "__main__":
    main()
