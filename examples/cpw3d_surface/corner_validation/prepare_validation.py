#!/usr/bin/env python3

"""Prepare corrected-thin and fabricated-reference corner validation runs."""

import argparse
import json
from pathlib import Path


# Defaults of the prototype process; a library with a Fabrication record (the version-2 coupon
# libraries) overrides the permittivities and thicknesses so that the corrected run passes the
# library's interface-layer validation (loss tangents are not part of the record).
INTERFACES = {
    "SA": (4.0, 2.0e-3, 2.0e-3),
    "MS": (11.47, 3.0e-4, 2.0e-3),
    "MA": (10.0, 3.0e-2, 2.0e-3),
}
SUBSTRATE_PERMITTIVITY = 11.47


def process_from_library(library):
    """(interfaces, substrate permittivity) of the library's Fabrication record, else the
    prototype defaults."""
    fabrication = json.loads(library.read_text()).get("Fabrication") or {}
    interfaces = dict(INTERFACES)
    for name, layer in (fabrication.get("InterfaceLayers") or {}).items():
        if name in interfaces:
            permittivity, loss_tangent, thickness = interfaces[name]
            interfaces[name] = (
                float(layer.get("Permittivity", permittivity)),
                loss_tangent,
                float(layer.get("Thickness", thickness)),
            )
    substrate = float(fabrication.get("SubstratePermittivity", SUBSTRATE_PERMITTIVITY))
    return interfaces, substrate


# Physical groups of mesh_corner_validation.jl / mesh_polygon_device.jl: 2D 1 outer, 3 SA,
# 7 outer_truncation (aperture); fabricated 2 MS / 4 MA, thin 2 thin_metal.
OUTER_ATTRIBUTES = [1, 7]
SA_ATTRIBUTE = 3


def dielectric(
    index,
    attributes,
    interface_type,
    radius,
    fabricated,
    edge_elements_per_radius,
    interfaces=INTERFACES,
):
    """Interface record with its edge lines at `radius`: the thin sheet's metal edge is the
    automatic perimeter (the device rule); the fabricated slab has no one-sided metal edge (every
    edge of MS bottom / MA sidewalls / MA top is a fold between metal faces, so `AutomaticEdges`
    finds nothing), and its edge lines are named as the perimeter of the SA surface minus the
    outer box, as the fabricated corner coupon does (generate_corner_response.py)."""
    permittivity, loss_tangent, thickness = interfaces[interface_type]
    result = {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": thickness,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "LocalizeEdgeEnergy": False,
        "EdgeExcludeAttributes": OUTER_ATTRIBUTES,
        "EdgeDistances": [radius],
    }
    if fabricated:
        result["EdgeAttributes"] = [SA_ATTRIBUTE]
    else:
        result["AutomaticEdges"] = True
    if edge_elements_per_radius:
        result["EdgeRefinement"] = {
            "Radius": radius,
            "ElementsPerRadius": edge_elements_per_radius,
            "OuterRadiusFactor": 1.0,
            "CoreIndicatorWeight": 0.0,
        }
    return result


def config(
    output,
    mesh,
    order,
    amr_iterations,
    fabricated,
    library,
    radius,
    edge_elements_per_radius,
    process=(INTERFACES, SUBSTRATE_PERMITTIVITY),
):
    refinement = 0 if fabricated else edge_elements_per_radius
    layers, substrate_permittivity = process
    interfaces = (
        [
            dielectric(1, [SA_ATTRIBUTE], "SA", radius, True, refinement, layers),
            dielectric(2, [2], "MS", radius, True, refinement, layers),
            dielectric(3, [4], "MA", radius, True, refinement, layers),
        ]
        if fabricated
        else [
            dielectric(1, [SA_ATTRIBUTE], "SA", radius, False, refinement, layers),
            dielectric(2, [2], "MS", radius, False, refinement, layers),
            dielectric(3, [2], "MA", radius, False, refinement, layers),
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
                {"Attributes": [1], "Permittivity": substrate_permittivity},
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
    parser.add_argument("--edge-elements-per-radius", type=int, default=1)
    args = parser.parse_args()
    if args.order < 1:
        parser.error("--order must be positive")
    if args.amr_iterations < 0:
        parser.error("--amr-iterations must be nonnegative")
    if args.edge_elements_per_radius < 1:
        parser.error("--edge-elements-per-radius must be positive")

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
    process = process_from_library(library)

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
            args.edge_elements_per_radius,
            process,
        )
        path = output / f"{name}.json"
        path.write_text(json.dumps(data, indent=2) + "\n")
        print(path)


if __name__ == "__main__":
    main()
