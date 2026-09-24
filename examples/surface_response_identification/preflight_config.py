#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Write a geometry-only surface-response preflight configuration for any Palace mesh.

The preflight (palace --surface-response-preflight CONFIG) needs: the metal attributes as
Ground (PEC) and Terminal boundaries (Electrostatic problem, so that Terminal indices give
conductor identity), typed SA / MS / MA dielectric interfaces with AutomaticEdges and the
library radius as the last EdgeDistances value, and the response-correction block naming the
library and the target interfaces. Nothing here depends on a particular device.

    python3 -m surface_response_identification.preflight_config --mesh M.msh2 \\
        --ground 5 6 7 --terminal 9 --sa 8 --library lib.json --output cfg.json \\
        [--ms 5 6 7 9] [--ma ...] [--uniform-levels 1] [--no-crack] [--l0 1e-6] [--radius 2]
"""

import argparse
import json
import os


def preflight_config(
    mesh,
    ground,
    terminals,
    sa,
    library,
    output_directory,
    ms=None,
    ma=None,
    radius=2.0,
    l0=1.0e-6,
    uniform_levels=0,
    crack=True,
    substrate_attributes=(1,),
    vacuum_attributes=(2,),
    substrate_permittivity=11.45,
    layer_thickness=0.002,
    sa_permittivity=4.0,
    ms_permittivity=11.47,
    ma_permittivity=10.0,
    order=1,
):
    metal = sorted(set(ground) | {a for group in terminals for a in group})
    ms = sorted(ms) if ms is not None else metal
    ma = sorted(ma) if ma is not None else metal
    # The classifier checks the interface layer thickness / permittivity of every target
    # against the library's Fabrication.InterfaceLayers to 1e-10: take them from the library.
    with open(library) as source:
        library_data = json.load(source)
    layers = library_data.get("Fabrication", {}).get("InterfaceLayers", {})
    radius = library_data.get("MatchingRadius", radius)
    thickness = {k: layers.get(k, {}).get("Thickness", layer_thickness) for k in ("SA", "MS", "MA")}
    permittivity = {
        "SA": layers.get("SA", {}).get("Permittivity", sa_permittivity),
        "MS": layers.get("MS", {}).get("Permittivity", ms_permittivity),
        "MA": layers.get("MA", {}).get("Permittivity", ma_permittivity),
    }
    substrate_permittivity = library_data.get("Fabrication", {}).get("SubstratePermittivity", substrate_permittivity)

    def dielectric(index, attributes, kind):
        return {
            "Index": index,
            "Attributes": list(attributes),
            "Type": kind,
            "Thickness": thickness[kind],
            "Permittivity": permittivity[kind],
            "LossTan": 1.0e-3,
            "AutomaticEdges": True,
            "LocalizeEdgeEnergy": False,
            "SaveLocalEdgeEnergy": False,
            "EdgeDistances": [radius],
        }

    dielectrics = []
    targets = []
    if sa:
        dielectrics.append(dielectric(1, sa, "SA"))
        targets.append(1)
    dielectrics.append(dielectric(2, ms, "MS"))
    dielectrics.append(dielectric(3, ma, "MA"))
    targets += [2, 3]
    boundaries = {"Ground": {"Attributes": sorted(ground)}, "Postprocessing": {"Dielectric": dielectrics}}
    if terminals:
        boundaries["Terminal"] = [{"Index": i + 1, "Attributes": sorted(group)} for i, group in enumerate(terminals)]
    model = {"Mesh": os.path.abspath(mesh), "L0": l0, "Refinement": {"MaxIts": 0, "UniformLevels": uniform_levels, "SerialUniformLevels": 0}}
    if not crack:
        model["CrackInternalBoundaryElements"] = False
    return {
        "Problem": {"Type": "Electrostatic", "Verbose": 1, "Output": os.path.abspath(output_directory), "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": model,
        "Domains": {
            "Materials": [
                {"Attributes": list(substrate_attributes), "Permittivity": substrate_permittivity},
                {"Attributes": list(vacuum_attributes), "Permittivity": 1.0},
            ]
        },
        "Boundaries": boundaries,
        "Solver": {
            "Order": order,
            "Electrostatic": {
                "Save": 0,
                "ResponseCorrection": {"Library": os.path.abspath(library), "TargetInterfaces": targets, "UnmatchedPolicy": "Warn"},
            },
            "Linear": {"Type": "BoomerAMG", "KSPType": "CG", "Tol": 1.0e-8, "MaxIts": 100},
        },
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--ground", type=int, nargs="+", required=True, help="grounded metal attributes (conductor 0)")
    parser.add_argument("--terminal", type=int, nargs="+", action="append", default=[], help="attributes of one separately identified conductor (repeatable)")
    parser.add_argument("--sa", type=int, nargs="*", default=[], help="substrate-air interface attributes (none: MS/MA only)")
    parser.add_argument("--ms", type=int, nargs="*", help="MS attributes (default: all metal)")
    parser.add_argument("--ma", type=int, nargs="*", help="MA attributes (default: all metal)")
    parser.add_argument("--substrate", type=int, nargs="+", default=[1])
    parser.add_argument("--vacuum", type=int, nargs="+", default=[2])
    parser.add_argument("--library", required=True)
    parser.add_argument("--output", required=True, help="configuration file to write")
    parser.add_argument("--postpro", help="Problem.Output (default: <output dir>/postpro)")
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--l0", type=float, default=1.0e-6)
    parser.add_argument("--uniform-levels", type=int, default=0)
    parser.add_argument("--no-crack", action="store_true", help="Model.CrackInternalBoundaryElements false")
    args = parser.parse_args(argv)
    output = os.path.abspath(args.output)
    postpro = args.postpro or os.path.join(os.path.dirname(output), "postpro")
    config = preflight_config(
        args.mesh, args.ground, args.terminal, args.sa, args.library, postpro, ms=args.ms, ma=args.ma, radius=args.radius, l0=args.l0,
        uniform_levels=args.uniform_levels, crack=not args.no_crack, substrate_attributes=args.substrate, vacuum_attributes=args.vacuum,
    )
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as target:
        json.dump(config, target, indent=2)
        target.write("\n")
    print(output)


if __name__ == "__main__":
    main()
