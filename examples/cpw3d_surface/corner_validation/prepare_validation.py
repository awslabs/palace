#!/usr/bin/env python3

"""Prepare corrected-thin and fabricated-reference corner validation runs."""

import argparse
import json
from pathlib import Path


INTERFACES = {
    "SA": (4.0, 2.0e-3),
    "MS": (11.47, 3.0e-4),
    "MA": (10.0, 3.0e-2),
}


def dielectric(index, attributes, interface_type, radius, automatic):
    permittivity, loss_tangent = INTERFACES[interface_type]
    result = {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": 2.0e-3,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
    }
    if automatic:
        result.update(
            {
                "AutomaticEdges": True,
                "LocalizeEdgeEnergy": True,
                "EdgeExcludeAttributes": [1],
                "EdgeDistances": [
                    0.125 * radius,
                    0.25 * radius,
                    0.5 * radius,
                    0.75 * radius,
                    radius,
                ],
            }
        )
    return result


def config(output, mesh, order, amr_iterations, fabricated, library, radius):
    interfaces = (
        [
            dielectric(1, [3], "SA", radius, True),
            dielectric(2, [2], "MS", radius, True),
            dielectric(3, [4], "MA", radius, True),
        ]
        if fabricated
        else [
            dielectric(1, [3], "SA", radius, True),
            dielectric(2, [2], "MS", radius, True),
            dielectric(3, [2], "MA", radius, True),
        ]
    )
    boundaries = {
        "Ground": {"Attributes": [1]},
        "Terminal": [
            {
                "Index": 1,
                "Attributes": [2, 4] if fabricated else [2],
            }
        ],
        "Postprocessing": {"Dielectric": interfaces},
    }
    solver = {
        "Order": order,
        "Electrostatic": {"Save": 0},
        "Linear": {
            "Type": "BoomerAMG",
            "KSPType": "CG",
            "Tol": 1.0e-10,
            "MaxIts": 1000,
            "EstimatorTol": 1.0e-6 if amr_iterations else 1.0e-1,
            "EstimatorMaxIts": 500 if amr_iterations else 20,
            "EstimatorMG": True,
        },
    }
    if not fabricated:
        solver["Electrostatic"]["ResponseCorrection"] = {
            "Library": str(library),
            "TargetInterfaces": [1, 2, 3],
            "UnmatchedPolicy": "Error",
        }
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 2,
            "Output": str(output),
            "OutputFormats": {"Paraview": False, "GridFunction": False},
        },
        "Model": {
            "Mesh": str(mesh),
            "L0": 1.0e-6,
            "Refinement": {"Tol": 1.0e-12, "MaxIts": amr_iterations},
        },
        "Domains": {
            "Materials": [
                {"Attributes": [1], "Permittivity": 11.47},
                {"Attributes": [2], "Permittivity": 1.0},
            ],
            "Postprocessing": {
                "Energy": [
                    {"Index": 1, "Attributes": [1]},
                    {"Index": 2, "Attributes": [2]},
                ]
            },
        },
        "Boundaries": boundaries,
        "Solver": solver,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path, required=True)
    parser.add_argument("--fabricated-mesh", type=Path, required=True)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--amr-iterations", type=int, default=0)
    args = parser.parse_args()
    if args.order < 1:
        parser.error("--order must be positive")
    if args.amr_iterations < 0:
        parser.error("--amr-iterations must be nonnegative")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    library = args.library.expanduser().resolve()
    thin_mesh = args.thin_mesh.expanduser().resolve()
    fabricated_mesh = args.fabricated_mesh.expanduser().resolve()
    for path in (library, thin_mesh, fabricated_mesh):
        if not path.is_file():
            raise FileNotFoundError(path)
    radius = float(json.loads(library.read_text())["MatchingRadius"])
    if radius <= 0.0:
        raise ValueError("The process library MatchingRadius must be positive")

    for name, mesh, fabricated in (
        ("thin-corrected", thin_mesh, False),
        ("fabricated-reference", fabricated_mesh, True),
    ):
        data = config(
            output / "postpro" / name,
            mesh,
            args.order,
            args.amr_iterations,
            fabricated,
            library,
            radius,
        )
        path = output / f"{name}.json"
        path.write_text(json.dumps(data, indent=2) + "\n")
        print(path)


if __name__ == "__main__":
    main()
