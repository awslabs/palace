#!/usr/bin/env python3

"""Prepare smooth driven-Maxwell rounded-corner validation runs."""

import argparse
import json
import sys
from pathlib import Path

IDENTIFICATION_TOOLS = Path(__file__).resolve().parents[2] / "surface_response_identification"
if str(IDENTIFICATION_TOOLS) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION_TOOLS))
from crack_policy import add_crack_policy_argument, apply_crack_policy  # noqa: E402
import developer_form  # noqa: E402

# The v1 default path (decision 2379 (3); CUT-PLAN step 5a; developer_form.py): by default this device-run writer
# writes the V1 FORM, uncracked along the perfect conductors (Palace refuses the global crack key there without the
# developer switch); --developer writes the RECORDED DEVELOPER FORM (the records' cracked mesh, crack policy legacy,
# and the corrected block with its default of record "CorrectionMode": "Both" and "Developer": true); --crack-policy
# legacy needs --developer. crack_policy.py (fingerprinted into the coupon cache keys) and the coupon generators keep
# their legacy default. Values from the v1 form are not comparable with records made before this cut until the
# benchmark's collocated twins are read.
DEVICE_RUN_CRACK_POLICY = developer_form.V1_CRACK_POLICY


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


def config(output, mesh, order, frequency, fabricated, library, radius, crack_policy=DEVICE_RUN_CRACK_POLICY,
           developer=False):
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
        recorded = {
            "Library": str(library),
            "TargetInterfaces": [1, 2, 3],
            "UnmatchedPolicy": "Error",
        }
        solver["SurfaceResponseCorrection"] = (
            developer_form.developer_block(recorded) if developer else developer_form.v1_block(recorded)
        )
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
    developer_form.add_developer_argument(parser)
    parser.set_defaults(crack_policy=None)  # the form's crack policy
    args = parser.parse_args()
    try:
        args.crack_policy = developer_form.device_run_crack_policy(args.crack_policy, args.developer)
    except ValueError as error:
        parser.error(str(error))
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
            developer=args.developer,
        )
        path = output / f"{name}.json"
        path.write_text(json.dumps(data, indent=2) + "\n")
        print(path)
    print(f"prepare_maxwell_validation.py: the {developer_form.form_name(args.developer)}", file=sys.stderr)


if __name__ == "__main__":
    main()
