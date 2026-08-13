#!/usr/bin/env python3

"""Generate contour hat bases and Palace configs for one isolated-edge coupon."""

import argparse
import json
from pathlib import Path

import numpy as np


INTERFACE_PROPERTIES = {
    "SA": (4.0, 0.002),
    "MS": (11.47, 0.0003),
    "MA": (10.0, 0.03),
}


def contour_point(distance, radius):
    width = 2.0 * radius
    perimeter = 4.0 * width
    distance %= perimeter
    if distance < width:
        return np.array([-radius + distance, -radius, 0.0])
    distance -= width
    if distance < width:
        return np.array([radius, -radius + distance, 0.0])
    distance -= width
    if distance < width:
        return np.array([radius - distance, radius, 0.0])
    distance -= width
    return np.array([-radius, radius - distance, 0.0])


def write_bases(output, radius, metal_thickness, basis_size, samples):
    if basis_size <= 0:
        raise ValueError("basis size must be positive")
    if metal_thickness < 0.0 or metal_thickness >= 2.0 * radius:
        raise ValueError("metal thickness must be nonnegative and smaller than 2R")
    perimeter = 8.0 * radius
    spacing = perimeter / basis_size

    # The matching contour is traversed counterclockwise from its lower-left
    # corner. It meets the lower and upper fabricated-metal surfaces while
    # descending the left side. Insert both junctions and constrain the metal
    # interval so the thin and fabricated coupons have one common trace space.
    lower_junction = 7.0 * radius
    upper_junction = lower_junction - metal_thickness
    offset = lower_junction % spacing
    uniform_knots = (offset + spacing * np.arange(basis_size)) % perimeter
    knot_distances = np.unique(
        np.concatenate((uniform_knots, [upper_junction, lower_junction]))
    )
    knot_distances.sort()
    tolerance = 1.0e-12 * perimeter
    constrained = np.asarray(
        [
            upper_junction - tolerance
            <= distance
            <= lower_junction + tolerance
            for distance in knot_distances
        ]
    )
    free_indices = [index for index, fixed in enumerate(constrained) if not fixed]
    distances = np.unique(
        np.concatenate(
            (
                np.linspace(0.0, perimeter, samples, endpoint=False),
                knot_distances,
            )
        )
    )
    points = np.asarray([contour_point(distance, radius) for distance in distances])
    knots = np.asarray(
        [contour_point(distance, radius) for distance in knot_distances]
    )
    np.savetxt(
        output / "basis_points.csv",
        knots[free_indices],
        delimiter=",",
        header="x,y,z",
        comments="",
        fmt="%.16e",
    )

    for stale_path in output.glob("basis_hat*.csv"):
        stale_path.unlink()

    def hat_values(index):
        previous = knot_distances[index - 1]
        current = knot_distances[index]
        following = knot_distances[(index + 1) % len(knot_distances)]
        left_span = (current - previous) % perimeter
        right_span = (following - current) % perimeter
        backward = (current - distances) % perimeter
        forward = (distances - current) % perimeter
        left = np.where(
            backward <= left_span + tolerance, 1.0 - backward / left_span, 0.0
        )
        right = np.where(
            forward <= right_span + tolerance, 1.0 - forward / right_span, 0.0
        )
        return np.maximum(left, right)

    paths = []
    for output_index, knot_index in enumerate(free_indices, start=1):
        path = output / f"basis_hat{output_index:03d}.csv"
        np.savetxt(
            path,
            np.column_stack((points, hat_values(knot_index))),
            delimiter=",",
            header="x,y,z,V",
            comments="",
            fmt="%.16e",
        )
        paths.append(path)
    return paths


def write_heldout(output, traces, radius, metal_thickness):
    basis_points = np.atleast_2d(
        np.loadtxt(output / "basis_points.csv", delimiter=",", skiprows=1)
    )
    samples = [
        np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
        for path in traces
    ]
    coordinates = samples[0][:, :3]
    if any(
        sample.shape != samples[0].shape
        or not np.allclose(sample[:, :3], coordinates, rtol=0.0, atol=1.0e-14)
        for sample in samples[1:]
    ):
        raise ValueError("Contour basis traces do not share one sampling grid")

    def potential(points):
        x = points[:, 0] / radius
        y = points[:, 1] / radius
        distance_to_cut = np.hypot(
            points[:, 0] + radius,
            np.maximum.reduce(
                (
                    -points[:, 1],
                    points[:, 1] - metal_thickness,
                    np.zeros(len(points)),
                )
            ),
        )
        coordinate = np.clip(distance_to_cut / (radius / 3.0), 0.0, 1.0)
        cutoff = coordinate * coordinate * (3.0 - 2.0 * coordinate)
        return cutoff * (
            0.35
            + 0.20 * x
            - 0.15 * y
            + 0.08 * x * y
            + 0.06 * y * y
        )

    coefficients = potential(basis_points)
    values = potential(coordinates)
    trace = output / "heldout_trace.csv"
    np.savetxt(
        trace,
        np.column_stack((coordinates, values)),
        delimiter=",",
        header="x,y,z,V",
        comments="",
        fmt="%.16e",
    )
    np.savetxt(
        output / "heldout_coefficients.csv",
        coefficients,
        delimiter=",",
        header="coefficient_V",
        comments="",
        fmt="%.16e",
    )
    return trace


def dielectric(index, attributes, interface_type, thickness, permittivity):
    _, loss_tangent = INTERFACE_PROPERTIES[interface_type]
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": thickness,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "EdgeAttributes": [2],
        "EdgeExcludeAttributes": [1],
        "EdgeDistances": [0.2],
        "LocalizeEdgeEnergy": True,
        "SaveLocalEdgeEnergy": False,
        "EdgeFrameNormal": [0.0, 1.0, 0.0],
    }


def make_config(
    output,
    name,
    mesh,
    traces,
    fabricated,
    order,
    coupon_depth,
    substrate_permittivity,
    interface_layers,
):
    if fabricated:
        ground = [2, 3, 4]
        interfaces = [
            dielectric(1, [5, 6], "SA", *interface_layers["SA"]),
            dielectric(2, [2], "MS", *interface_layers["MS"]),
            dielectric(3, [3, 4], "MA", *interface_layers["MA"]),
        ]
    else:
        ground = [2]
        interfaces = [
            dielectric(1, [3], "SA", *interface_layers["SA"]),
            dielectric(2, [2], "MS", *interface_layers["MS"]),
            dielectric(3, [2], "MA", *interface_layers["MA"]),
        ]
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
        },
        "Model": {"Mesh": str(mesh), "L0": 1.0e-6, "Lc": coupon_depth},
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
        "Boundaries": {
            "Ground": {"Attributes": ground},
            "PrescribedPotential": [
                {
                    "Index": index,
                    "Attributes": [1],
                    "DataFile": str(trace),
                }
                for index, trace in enumerate(traces, start=1)
            ],
            "Postprocessing": {"Dielectric": interfaces},
        },
        "Solver": {
            "Order": order,
            "Device": "CPU",
            "Electrostatic": {
                "Save": 0,
                "ResponseMatrix": True,
                "AggregateResponseMatrix": True,
            },
            "Linear": {
                "Type": "BoomerAMG",
                "KSPType": "CG",
                "Tol": 1.0e-10,
                "MaxIts": 1000,
                "EstimatorTol": 1.0e-2,
                "EstimatorMaxIts": 20,
                "EstimatorMG": True,
            },
        },
    }


def write_library(
    output,
    name,
    radius,
    coupon_depth,
    metal_thickness,
    substrate_permittivity,
    interface_layers,
):
    library = {
        "Version": 3,
        "Name": name,
        "MatchingRadius": radius,
        "CouponDepth": coupon_depth,
        "Fabrication": {
            "LengthUnit": "um",
            "MetalThickness": metal_thickness,
            "SubstratePermittivity": substrate_permittivity,
            "InterfaceLayers": {
                interface_type: {
                    "Thickness": layer[0],
                    "Permittivity": layer[1],
                }
                for interface_type, layer in interface_layers.items()
            },
        },
        "Models": [
            {
                "Name": "isolated-edge",
                "Topology": "IsolatedEdge",
                "FabricatedMatrix":
                    "postpro/edge_fabricated/domain-response-matrix.csv",
                "ThinMatrix": "postpro/edge_thin/domain-response-matrix.csv",
                "FabricatedSurfaceMatrix":
                    "postpro/edge_fabricated/surface-response-matrix.csv",
                "ThinSurfaceMatrix":
                    "postpro/edge_thin/surface-response-matrix.csv",
                "BasisPoints": "basis_points.csv",
                "Interfaces": [
                    {"Type": "SA", "Coupon": 1},
                    {"Type": "MS", "Coupon": 2},
                    {"Type": "MA", "Coupon": 3},
                ],
            }
        ],
    }
    path = output / "process-library.json"
    path.write_text(json.dumps(library, indent=2) + "\n")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--basis-size", type=int, default=96)
    parser.add_argument("--samples", type=int, default=1200)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--coupon-depth", type=float, default=1055.0)
    parser.add_argument("--substrate-permittivity", type=float, default=11.47)
    parser.add_argument("--sa-thickness", type=float, default=0.002)
    parser.add_argument("--sa-permittivity", type=float, default=4.0)
    parser.add_argument("--ms-thickness", type=float, default=0.002)
    parser.add_argument("--ms-permittivity", type=float, default=11.47)
    parser.add_argument("--ma-thickness", type=float, default=0.002)
    parser.add_argument("--ma-permittivity", type=float, default=10.0)
    parser.add_argument("--library-name", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path)
    parser.add_argument("--fabricated-mesh", type=Path)
    args = parser.parse_args()
    material_values = (
        args.substrate_permittivity,
        args.sa_thickness,
        args.sa_permittivity,
        args.ms_thickness,
        args.ms_permittivity,
        args.ma_thickness,
        args.ma_permittivity,
    )
    if any(value <= 0.0 for value in material_values):
        parser.error("substrate and interface-layer properties must be positive")
    interface_layers = {
        "SA": (args.sa_thickness, args.sa_permittivity),
        "MS": (args.ms_thickness, args.ms_permittivity),
        "MA": (args.ma_thickness, args.ma_permittivity),
    }

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    thin_mesh = (
        args.thin_mesh.expanduser().resolve()
        if args.thin_mesh
        else output / "edge_thin.msh"
    )
    fabricated_mesh = (
        args.fabricated_mesh.expanduser().resolve()
        if args.fabricated_mesh
        else output / "edge_fabricated.msh"
    )
    traces = write_bases(
        output,
        args.radius,
        args.metal_thickness,
        args.basis_size,
        args.samples,
    )
    heldout_trace = write_heldout(
        output, traces, args.radius, args.metal_thickness
    )
    for name, mesh, fabricated in (
        ("edge_thin", thin_mesh, False),
        ("edge_fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            traces,
            fabricated,
            args.order,
            args.coupon_depth,
            args.substrate_permittivity,
            interface_layers,
        )
        path = output / f"{name}.json"
        path.write_text(json.dumps(config, indent=2) + "\n")
        print(path)
        heldout_name = f"heldout_{name}"
        heldout = make_config(
            output,
            heldout_name,
            mesh,
            [heldout_trace],
            fabricated,
            args.order,
            args.coupon_depth,
            args.substrate_permittivity,
            interface_layers,
        )
        heldout["Solver"]["Electrostatic"]["ResponseMatrix"] = False
        heldout["Solver"]["Electrostatic"]["AggregateResponseMatrix"] = False
        heldout_path = output / f"{heldout_name}.json"
        heldout_path.write_text(json.dumps(heldout, indent=2) + "\n")
        print(heldout_path)
    print(output / "basis_points.csv")
    print(
        write_library(
            output,
            args.library_name,
            args.radius,
            args.coupon_depth,
            args.metal_thickness,
            args.substrate_permittivity,
            interface_layers,
        )
    )


if __name__ == "__main__":
    main()
