#!/usr/bin/env python3

"""Prepare smooth driven-Maxwell rounded-corner validation runs."""

import argparse
import json
import sys
from pathlib import Path

IDENTIFICATION_TOOLS = Path(__file__).resolve().parents[2] / "surface_response_identification"
if str(IDENTIFICATION_TOOLS) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION_TOOLS))
from crack_policy import CRACK_POLICY_DEFAULT, add_crack_policy_argument, apply_crack_policy  # noqa: E402


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
                "EdgeExcludeAttributes": [1, 5],
                "EdgeDistances": [radius],
                "LocalizeEdgeEnergy": False,
            }
        )
    return result


def config(output, mesh, order, frequency, fabricated, library, radius, crack_policy=CRACK_POLICY_DEFAULT):
    interfaces = (
        [
            dielectric(1, [3], "SA", radius, False),
            dielectric(2, [2], "MS", radius, False),
            dielectric(3, [4], "MA", radius, False),
        ]
        if fabricated
        else [
            dielectric(1, [3], "SA", radius, True),
            dielectric(2, [2], "MS", radius, True),
            dielectric(3, [2], "MA", radius, True),
        ]
    )
    solver = {
        "Order": order,
        "Driven": {
            "Samples": [{"Type": "Point", "Freq": [frequency], "SaveStep": 0}]
        },
        "Linear": {
            "Type": "STRUMPACK",
            "KSPType": "GMRES",
            "Tol": 1.0e-8,
            "MaxIts": 500,
            "EstimatorTol": 1.0e-1,
            "EstimatorMG": True,
        },
    }
    if not fabricated and library is not None:
        solver["SurfaceResponseCorrection"] = {
            "Library": str(library),
            "TargetInterfaces": [1, 2, 3],
            "UnmatchedPolicy": "Error",
        }
    return {
        "Problem": {
            "Type": "Driven",
            "Verbose": 2,
            "Output": str(output),
            "OutputFormats": {"Paraview": False, "GridFunction": False},
        },
        "Model": apply_crack_policy({
            "Mesh": str(mesh),
            "L0": 1.0e-6,
            "Refinement": {"MaxIts": 0},
        }, crack_policy),
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
        "Boundaries": {
            "PEC": {"Attributes": [1, 2, 4] if fabricated else [1, 2]},
            "LumpedPort": [
                {
                    "Index": 1,
                    "Attributes": [5],
                    "Direction": "+X",
                    "R": 50.0,
                    "Excitation": 1,
                }
            ],
            "Postprocessing": {"Dielectric": interfaces},
        },
        "Solver": solver,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path, required=True)
    parser.add_argument("--fabricated-mesh", type=Path, required=True)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--frequency", type=float, default=50.0)
    add_crack_policy_argument(parser)
    args = parser.parse_args()
    if args.order < 1:
        parser.error("--order must be positive")
    if args.frequency <= 0.0:
        parser.error("--frequency must be positive")

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

    cases = (
        ("thin-raw", thin_mesh, False, None),
        ("thin-corrected", thin_mesh, False, library),
        ("fabricated-reference", fabricated_mesh, True, None),
    )
    for name, mesh, fabricated, correction_library in cases:
        data = config(
            output / "postpro" / name,
            mesh,
            args.order,
            args.frequency,
            fabricated,
            correction_library,
            radius,
            crack_policy=args.crack_policy,
        )
        path = output / f"{name}.json"
        path.write_text(json.dumps(data, indent=2) + "\n")
        print(path)


if __name__ == "__main__":
    main()
