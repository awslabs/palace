#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare the conductor-tagged thin transmon for consistent electrostatic correction."""

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--postpro", type=Path, required=True)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--substrate-permittivity", type=float, default=11.45)
    parser.add_argument(
        "--source", type=Path, default=ROOT / "transmon_surface_coarse.json"
    )
    parser.add_argument(
        "--correction-mode",
        choices=("PostprocessOnly", "SelfConsistent", "Both"),
        default="Both",
    )
    args = parser.parse_args()
    if args.order <= 0 or args.substrate_permittivity <= 0.0:
        parser.error("order and substrate permittivity must be positive")

    config = json.loads(args.source.expanduser().resolve().read_text())
    config["Problem"] = {
        "Type": "Electrostatic",
        "Verbose": 1,
        "Output": str(args.postpro.expanduser().resolve()),
        "OutputFormats": {"Paraview": False, "GridFunction": False},
    }
    config["Model"]["Mesh"] = str(args.mesh.expanduser().resolve())
    config["Model"]["Refinement"] = {
        "MaxIts": 0,
        "UniformLevels": 0,
        "SerialUniformLevels": 0,
    }
    config["Domains"]["Materials"] = [
        {"Attributes": [1], "Permittivity": args.substrate_permittivity},
        {"Attributes": [2], "Permittivity": 1.0},
    ]
    postprocessing = config["Boundaries"]["Postprocessing"]
    for dielectric in postprocessing["Dielectric"]:
        dielectric["EdgeDistances"] = [2.0]
        dielectric["LocalizeEdgeEnergy"] = False
        dielectric["SaveLocalEdgeEnergy"] = False
        if dielectric["Type"] in ("MS", "MA"):
            dielectric["Attributes"] = [5, 6, 7, 9]
        if dielectric["Type"] == "MS":
            dielectric["Permittivity"] = args.substrate_permittivity
    config["Boundaries"] = {
        # Port patches 6 and 7 replace part of the metal sheet and must remain in the
        # grounded conductor union. Attribute 9 is the separately tagged transmon island.
        "Ground": {"Attributes": [5, 6, 7]},
        "Terminal": [{"Index": 1, "Attributes": [9]}],
        "Postprocessing": postprocessing,
    }
    config["Solver"] = {
        "Order": args.order,
        "Electrostatic": {
            "Save": 1,
            "ResponseCorrection": {
                "Library": str(args.library.expanduser().resolve()),
                "TargetInterfaces": [1, 2, 3],
                "UnmatchedPolicy": "Error",
                "CorrectionMode": args.correction_mode,
                "SolveTol": 1.0e-6,
            },
        },
        "Linear": {
            "Type": "BoomerAMG",
            "KSPType": "CG",
            "Tol": 1.0e-10,
            "MaxIts": 1000,
            "EstimatorTol": 1.0e-1,
            "EstimatorMG": True,
        },
    }
    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(config, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
