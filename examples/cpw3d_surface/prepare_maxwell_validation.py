#!/usr/bin/env python3

"""Prepare compact Maxwell CPW surface-response validation configurations."""

import argparse
import copy
import json
import sys
from pathlib import Path

IDENTIFICATION_TOOLS = Path(__file__).resolve().parents[1] / "surface_response_identification"
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


THIN_INTERFACES = (
    (1, [3], "SA", 4.0, 2.0e-3),
    (2, [1, 2], "MS", 11.47, 3.0e-4),
    (3, [1, 2], "MA", 10.0, 3.0e-2),
)

FABRICATED_INTERFACES = (
    (1, [5], "SA", 4.0, 2.0e-3),
    (2, [1, 3], "MS", 11.47, 3.0e-4),
    (3, [2, 4], "MA", 10.0, 3.0e-2),
)


def dielectric_entries(
    fabricated, automatic_edges, edge_exclude_attributes=()
):
    entries = []
    interfaces = FABRICATED_INTERFACES if fabricated else THIN_INTERFACES
    for index, attributes, interface_type, permittivity, loss_tangent in interfaces:
        entry = {
            "Index": index,
            "Attributes": attributes,
            "Type": interface_type,
            "Thickness": 2.0e-3,
            "Permittivity": permittivity,
            "LossTan": loss_tangent,
        }
        if automatic_edges:
            entry["AutomaticEdges"] = True
            if edge_exclude_attributes:
                entry["EdgeExcludeAttributes"] = list(edge_exclude_attributes)
            entry["EdgeDistances"] = [2.0]
        entries.append(entry)
    return entries


def base_config(
    mesh,
    output,
    fabricated,
    automatic_edges,
    order,
    problem,
    eigenmode_target,
    crack_policy=DEVICE_RUN_CRACK_POLICY,
):
    eigenmode = problem == "eigenmode"
    port_attributes = (7, 8) if fabricated else (5, 6)
    pec_attributes = [1, 2, 3, 4, 6] if fabricated else [1, 2, 4]
    boundaries = {
        "PEC": {"Attributes": pec_attributes},
        "Postprocessing": {
            "Dielectric": dielectric_entries(
                fabricated,
                automatic_edges,
                port_attributes if eigenmode else (),
            )
        },
    }
    if not eigenmode:
        boundaries["WavePort"] = [
            {
                "Index": 1,
                "Attributes": [port_attributes[0]],
                "Mode": 1,
                "Excitation": 1,
            },
            {"Index": 2, "Attributes": [port_attributes[1]], "Mode": 1},
        ]
    solver = {
        "Order": order,
        "Linear": {
            "Type": "STRUMPACK",
            "KSPType": "GMRES",
            "Tol": 1.0e-8,
            "MaxIts": 500,
            "EstimatorTol": 1.0e-1,
            "EstimatorMG": True,
        },
    }
    if eigenmode:
        solver["Eigenmode"] = {
            "N": 1,
            "Target": eigenmode_target,
            "Tol": 1.0e-6,
            "Save": 1,
        }
    else:
        solver["Driven"] = {
            "Samples": [{"Type": "Point", "Freq": [5.0], "SaveStep": 0}]
        }
    return {
        "Problem": {
            "Type": "Eigenmode" if eigenmode else "Driven",
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
        "Boundaries": boundaries,
        "Solver": solver,
    }


def write_config(path, config):
    path.write_text(json.dumps(config, indent=2) + "\n")
    print(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--library",
        type=Path,
        required=True,
        help="Fabrication-process response library visible to Palace",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--length", type=int, choices=(50, 200), default=50)
    parser.add_argument(
        "--problem",
        choices=("driven", "eigenmode"),
        default="driven",
        help="Maxwell validation problem to prepare",
    )
    parser.add_argument(
        "--eigenmode-target",
        type=float,
        default=100.0,
        help="Eigenmode target frequency in GHz",
    )
    parser.add_argument(
        "--mesh-directory",
        type=Path,
        default=Path(__file__).resolve().parent / "mesh",
    )
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
    if args.eigenmode_target < 0.0:
        parser.error("--eigenmode-target must be nonnegative")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    mesh_directory = args.mesh_directory.expanduser().resolve()
    suffix = "" if args.length == 50 else f"_L{args.length}"
    thin_mesh = mesh_directory / f"cpw3d_surface_maxwell_thin{suffix}.msh"
    fabricated_mesh = (
        mesh_directory / f"cpw3d_surface_maxwell_fabricated{suffix}.msh"
    )
    for mesh in (thin_mesh, fabricated_mesh):
        if not mesh.is_file():
            raise FileNotFoundError(
                f"Missing {mesh}; run mesh/generate_maxwell_validation_meshes.jl"
            )

    baseline = base_config(
        thin_mesh,
        output / "thin_raw",
        False,
        True,
        args.order,
        args.problem,
        args.eigenmode_target,
        crack_policy=args.crack_policy,
    )
    corrected = copy.deepcopy(baseline)
    corrected["Problem"]["Output"] = str(output / "thin_corrected")
    recorded = {
        "Library": str(args.library),
        "TargetInterfaces": [1, 2, 3],
        "UnmatchedPolicy": "Error",
    }
    corrected["Solver"]["SurfaceResponseCorrection"] = (
        developer_form.developer_block(recorded) if args.developer else developer_form.v1_block(recorded)
    )
    fabricated = base_config(
        fabricated_mesh,
        output / "fabricated_reference",
        True,
        False,
        args.order,
        args.problem,
        args.eigenmode_target,
        crack_policy=args.crack_policy,
    )

    write_config(output / "thin_raw.json", baseline)
    write_config(output / "thin_corrected.json", corrected)
    write_config(output / "fabricated_reference.json", fabricated)
    print(f"prepare_maxwell_validation.py: the {developer_form.form_name(args.developer)}", file=sys.stderr)


if __name__ == "__main__":
    main()
