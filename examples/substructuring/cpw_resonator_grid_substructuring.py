# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""
Prepare the CPW resonator grid mesh (cpw_resonator_grid.jl) for substructuring, and write the
configurations of the substructured runs from that of the regular one (cpw_resonator_grid.json).

The DeviceLayout mesh has no region/environment split. This script assigns the tetrahedra
whose centroid lies in a box around one resonator (by default resonator 6, in the middle row)
to the region (substrate 1, vacuum 2) and the rest to the environment (substrate 10, vacuum
11). The boundaries keep their attributes: 3 = exterior boundary, 4 = metal, 5, 6, ... = the
lumped ports. The region tetrahedra also get their own geometric (elementary) volumes, so that
viewers that show the mesh by geometric entity, like Gmsh, display the region apart from the
environment.

The configurations (cpw_resonator_grid_<run>.json) differ from the regular one by the
substructuring of the region and the environment, the output directory, and for each run:

  offline           condense the environment, save the model
  online            reuse the model: only the region is solved
  online_region     the same on the mesh of the region, which the offline run writes
  redesign          the region's substrate permittivity 0.4% higher
  offline_adaptive  sample a rational model of the environment over 15-21.5 GHz
  online_adaptive   the region against that model, reduced as in the adaptive fast
                    frequency sweep

The region is solved by FGMRES (whatever "KSPType"), so the configurations leave it out.

Usage (from this directory):

    python3 cpw_resonator_grid_substructuring.py [--input mesh/cpw_resonator_grid.msh2]
        [--output mesh/cpw_resonator_grid_substructuring.msh2]
        [--box XMIN XMAX YMIN YMAX ZMIN ZMAX]
"""

import argparse
import collections
import copy
import json
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
    write_configs(here)


def write_configs(here):
    with open(os.path.join(here, "cpw_resonator_grid.json")) as f:
        regular = json.load(f)

    def adaptive(config):
        driven = config["Solver"]["Driven"]
        driven["Samples"][0].update({"MinFreq": 15.0, "MaxFreq": 21.5})
        driven.update({"AdaptiveTol": 1.0e-8, "AdaptiveMaxSamples": 100})

    def region_mesh(config):
        config["Model"]["Mesh"] = "postpro/cpw_resonator_grid_offline/environment_region.mesh"

    def redesign(config):
        # The region's substrate as a material of its own, with permittivity 0.4% higher.
        materials = config["Domains"]["Materials"]
        substrate = next(m for m in materials if REGION["substrate"] in m["Attributes"])
        substrate["Attributes"].remove(REGION["substrate"])
        region = copy.deepcopy(substrate)
        region["Attributes"] = [REGION["substrate"]]
        region["Permittivity"] = [round(1.004 * e, 6) for e in substrate["Permittivity"]]
        materials.append(region)

    runs = [
        ("offline", None, None),
        ("online", "offline", None),
        ("online_region", "offline", region_mesh),
        ("redesign", "offline", redesign),
        ("offline_adaptive", None, adaptive),
        ("online_adaptive", "offline_adaptive", adaptive),
    ]
    for run, model, edit in runs:
        config = copy.deepcopy(regular)
        config["Problem"]["Output"] = "postpro/cpw_resonator_grid_" + run
        del config["Solver"]["Linear"]["KSPType"]
        substructuring = {
            "Region": {"Attributes": sorted(REGION.values())},
            "Environment": {"Attributes": sorted(ENVIRONMENT.values())},
        }
        if model:
            substructuring["ModelFile"] = (
                "postpro/cpw_resonator_grid_%s/environment.model" % model
            )
        config["Solver"]["Substructuring"] = substructuring
        if edit:
            edit(config)
        path = os.path.join(here, "cpw_resonator_grid_%s.json" % run)
        with open(path, "w") as f:
            json.dump(config, f, indent=2)
            f.write("\n")
        print("Wrote %s" % path)


if __name__ == "__main__":
    main()
