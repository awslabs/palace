#!/usr/bin/env python3

"""Generate contour hat bases and Palace configs for a coupled edge-pair coupon."""

import argparse
import json
from pathlib import Path

import numpy as np


INTERFACE_PROPERTIES = {
    "SA": (4.0, 0.002),
    "MS": (11.47, 0.0003),
    "MA": (10.0, 0.03),
}


def contour_point(distance, half_width, radius):
    width = 2.0 * half_width
    height = 2.0 * radius
    perimeter = 2.0 * (width + height)
    distance %= perimeter
    if distance < width:
        return np.array([-half_width + distance, -radius, 0.0])
    distance -= width
    if distance < height:
        return np.array([half_width, -radius + distance, 0.0])
    distance -= height
    if distance < width:
        return np.array([half_width - distance, radius, 0.0])
    distance -= width
    return np.array([-half_width, radius - distance, 0.0])


def write_bases(
    output,
    separation,
    radius,
    metal_thickness,
    basis_size,
    samples,
    different_conductors,
    strip,
):
    if basis_size % 2:
        raise ValueError("The paired-edge contour basis size must be even")
    if strip and different_conductors:
        raise ValueError("A physical metal strip cannot use different conductors")
    half_width = 0.5 * separation + radius
    perimeter = 4.0 * (half_width + radius)
    spacing = perimeter / basis_size
    if metal_thickness < 0.0 or metal_thickness >= 2.0 * radius:
        raise ValueError("metal thickness must be nonnegative and smaller than 2R")
    if strip:
        # A central strip does not touch the matching contour, so every contour
        # coefficient is free and no fabricated-metal junctions need insertion.
        knot_distances = spacing * np.arange(basis_size)
        constrained = np.zeros(basis_size, dtype=bool)
    else:
        # Align knots with the lower thin-metal junctions and insert the upper
        # fabricated-metal junctions. Constraining both ends of each conductor cut
        # gives the thin and fabricated coupons the same trace space.
        right_lower = 2.0 * half_width + radius
        right_upper = right_lower + metal_thickness
        left_lower = right_lower + 0.5 * perimeter
        left_upper = left_lower - metal_thickness
        offset = right_lower % spacing
        uniform_knots = (offset + spacing * np.arange(basis_size)) % perimeter
        knot_distances = np.unique(
            np.concatenate(
                (uniform_knots, [right_lower, right_upper, left_upper, left_lower])
            )
        )
        knot_distances.sort()

        tolerance = 1.0e-12 * perimeter

        def in_interval(value, start, end):
            return start - tolerance <= value <= end + tolerance

        constrained = np.asarray(
            [
                in_interval(distance, right_lower, right_upper)
                or in_interval(distance, left_upper, left_lower)
                for distance in knot_distances
            ]
        )
    distances = np.unique(
        np.concatenate(
            (
                np.linspace(0.0, perimeter, samples, endpoint=False),
                knot_distances,
            )
        )
    )
    points = np.asarray(
        [contour_point(distance, half_width, radius) for distance in distances]
    )
    knots = np.asarray(
        [contour_point(distance, half_width, radius) for distance in knot_distances]
    )

    tolerance = 1.0e-12 * perimeter

    def in_interval(value, start, end):
        return start - tolerance <= value <= end + tolerance

    free_indices = [index for index, fixed in enumerate(constrained) if not fixed]
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
    for output_index, index in enumerate(free_indices, start=1):
        values = hat_values(index)
        path = output / f"basis_hat{output_index:03d}.csv"
        np.savetxt(
            path,
            np.column_stack((points, values)),
            delimiter=",",
            header="x,y,z,V",
            comments="",
            fmt="%.16e",
        )
        paths.append(path)

    conductor_trace = None
    open_contour_paths = []
    if different_conductors:
        knot_values = np.asarray(
            [
                1.0 if in_interval(distance, right_lower, right_upper) else 0.0
                for distance in knot_distances
            ]
        )
        values = np.zeros_like(distances)
        for index, value in enumerate(knot_values):
            if value:
                values += value * hat_values(index)
        conductor_trace = output / "basis_conductor_state.csv"
        np.savetxt(
            conductor_trace,
            np.column_stack((points, values)),
            delimiter=",",
            header="x,y,z,V",
            comments="",
            fmt="%.16e",
        )
        output_indices = {
            knot_index: output_index
            for output_index, knot_index in enumerate(free_indices, start=1)
        }
        lower = [
            index
            for index in free_indices
            if knot_distances[index] < right_lower - tolerance
            or knot_distances[index] > left_lower + tolerance
        ]
        lower.sort(
            key=lambda index: (knot_distances[index] - left_lower) % perimeter
        )
        upper = [
            index
            for index in free_indices
            if right_upper + tolerance < knot_distances[index]
            and knot_distances[index] < left_upper - tolerance
        ]
        upper.sort(key=lambda index: knot_distances[index], reverse=True)
        if len(lower) + len(upper) != len(free_indices):
            raise ValueError("Unable to partition free knots into open contour paths")
        open_contour_paths = [
            {
                "Indices": [output_indices[index] for index in contour],
                "StartConductor": 1,
                "EndConductor": 2,
            }
            for contour in (lower, upper)
        ]
    return paths, conductor_trace, open_contour_paths


def write_heldout(
    output,
    traces,
    conductor_trace,
    separation,
    radius,
    metal_thickness,
    strip,
):
    basis_points = np.atleast_2d(
        np.loadtxt(output / "basis_points.csv", delimiter=",", skiprows=1)
    )
    paths = list(traces)
    if conductor_trace is not None:
        paths.append(conductor_trace)
    samples = [
        np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
        for path in paths
    ]
    coordinates = samples[0][:, :3]
    if any(
        sample.shape != samples[0].shape
        or not np.allclose(sample[:, :3], coordinates, rtol=0.0, atol=1.0e-14)
        for sample in samples[1:]
    ):
        raise ValueError("Contour basis traces do not share one sampling grid")

    half_width = 0.5 * separation + radius

    def free_potential(points):
        x = points[:, 0] / radius
        y = points[:, 1] / radius
        if strip:
            cutoff = np.ones(len(points))
        else:
            vertical_distance = np.maximum.reduce(
                (
                    -points[:, 1],
                    points[:, 1] - metal_thickness,
                    np.zeros(len(points)),
                )
            )
            left = np.hypot(points[:, 0] + half_width, vertical_distance)
            right = np.hypot(points[:, 0] - half_width, vertical_distance)
            coordinate = np.clip(
                np.minimum(left, right) / (radius / 3.0), 0.0, 1.0
            )
            cutoff = coordinate * coordinate * (3.0 - 2.0 * coordinate)
        return cutoff * (
            0.35
            + 0.20 * x
            - 0.15 * y
            + 0.08 * x * y
            + 0.06 * y * y
        )

    coefficients = list(free_potential(basis_points))
    values = free_potential(coordinates)
    if conductor_trace is not None:
        conductor_coefficient = 0.17
        coefficients.append(conductor_coefficient)
        values += conductor_coefficient * samples[-1][:, 3]
    coefficients = np.asarray(coefficients)
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


def dielectric(
    index, attributes, interface_type, edge_attributes, thickness, permittivity
):
    _, loss_tangent = INTERFACE_PROPERTIES[interface_type]
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": thickness,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "EdgeAttributes": edge_attributes,
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
    conductor_trace,
    fabricated,
    different_conductors,
    order,
    coupon_depth,
    substrate_permittivity,
    interface_layers,
):
    if fabricated:
        ground = [2, 4, 5, 7, 8, 9] if different_conductors else [2, 4, 5]
        edge_attributes = [2, 7] if different_conductors else [2]
        interfaces = [
            dielectric(
                1, [3, 6], "SA", edge_attributes, *interface_layers["SA"]
            ),
            dielectric(
                2,
                [2, 7] if different_conductors else [2],
                "MS",
                edge_attributes,
                *interface_layers["MS"],
            ),
            dielectric(
                3,
                [4, 5, 8, 9] if different_conductors else [4, 5],
                "MA",
                edge_attributes,
                *interface_layers["MA"],
            ),
        ]
    else:
        ground = [2, 7] if different_conductors else [2]
        edge_attributes = [2, 7] if different_conductors else [2]
        interfaces = [
            dielectric(1, [3], "SA", edge_attributes, *interface_layers["SA"]),
            dielectric(
                2,
                edge_attributes,
                "MS",
                edge_attributes,
                *interface_layers["MS"],
            ),
            dielectric(
                3,
                edge_attributes,
                "MA",
                edge_attributes,
                *interface_layers["MA"],
            ),
        ]
    potentials = [
        {
            "Index": index,
            "Attributes": [1],
            "DataFile": str(trace),
        }
        for index, trace in enumerate(traces, start=1)
    ]
    if different_conductors:
        potentials.append(
            {
                "Index": len(potentials) + 1,
                "Attributes": [1],
                "TerminalAttributes": [7, 8, 9] if fabricated else [7],
                "DataFile": str(conductor_trace),
            }
        )
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
            "PrescribedPotential": potentials,
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
    separation,
    separation_tolerance,
    radius,
    coupon_depth,
    different_conductors,
    strip,
    open_contour_paths,
    metal_thickness,
    substrate_permittivity,
    interface_layers,
):
    half_gap = 0.5 * separation
    topology = (
        "SameConductorStrip"
        if strip
        else "DifferentConductorGap"
        if different_conductors
        else "SameConductorGap"
    )
    topology_name = (
        "same-conductor-strip"
        if strip
        else "different-conductor-gap"
        if different_conductors
        else "same-conductor-gap"
    )
    model = {
        "Name": f"{topology_name}-{separation:g}um",
        "Topology": topology,
        "Separation": separation,
        "SeparationTolerance": separation_tolerance,
        "FabricatedMatrix": "postpro/edge_pair_fabricated/domain-response-matrix.csv",
        "ThinMatrix": "postpro/edge_pair_thin/domain-response-matrix.csv",
        "FabricatedSurfaceMatrix":
            "postpro/edge_pair_fabricated/surface-response-matrix.csv",
        "ThinSurfaceMatrix": "postpro/edge_pair_thin/surface-response-matrix.csv",
        "BasisPoints": "basis_points.csv",
        "Interfaces": [
            {"Type": "SA", "Coupon": 1},
            {"Type": "MS", "Coupon": 2},
            {"Type": "MA", "Coupon": 3},
        ],
    }
    reference_offset = half_gap + 0.5 * radius
    if different_conductors:
        model["ConductorReferences"] = [
            [-reference_offset, 0.0, 0.0],
            [reference_offset, 0.0, 0.0],
        ]
        model["OpenContourPaths"] = open_contour_paths
    elif strip:
        model["Reference"] = [0.0, 0.0, 0.0]
    else:
        model["Reference"] = [-reference_offset, 0.0, 0.0]
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
        "Models": [model],
    }
    path = output / "process-library.json"
    path.write_text(json.dumps(library, indent=2) + "\n")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--separation", "--cutout-width", dest="separation", type=float, required=True
    )
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
    parser.add_argument("--separation-tolerance", type=float, default=1.0e-3)
    parser.add_argument(
        "--library-name",
        default="100nm-metal-50nm-overetch-paired-edge-prototype",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path)
    parser.add_argument("--fabricated-mesh", type=Path)
    topology = parser.add_mutually_exclusive_group()
    topology.add_argument("--different-conductors", action="store_true")
    topology.add_argument("--strip", action="store_true")
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

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    thin_mesh = (
        args.thin_mesh.resolve()
        if args.thin_mesh
        else output / "edge_pair_thin.msh"
    )
    fabricated_mesh = (
        args.fabricated_mesh.resolve()
        if args.fabricated_mesh
        else output / "edge_pair_fabricated.msh"
    )
    traces, conductor_trace, open_contour_paths = write_bases(
        output,
        args.separation,
        args.radius,
        args.metal_thickness,
        args.basis_size,
        args.samples,
        args.different_conductors,
        args.strip,
    )
    heldout_trace = write_heldout(
        output,
        traces,
        conductor_trace,
        args.separation,
        args.radius,
        args.metal_thickness,
        args.strip,
    )
    for name, mesh, fabricated in (
        ("edge_pair_thin", thin_mesh, False),
        ("edge_pair_fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            traces,
            conductor_trace,
            fabricated,
            args.different_conductors,
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
            conductor_trace,
            fabricated,
            args.different_conductors,
            args.order,
            args.coupon_depth,
            args.substrate_permittivity,
            interface_layers,
        )
        potential = {
            "Index": 1,
            "Attributes": [1],
            "DataFile": str(heldout_trace),
        }
        if args.different_conductors:
            potential["TerminalAttributes"] = (
                [7, 8, 9] if fabricated else [7]
            )
        heldout["Boundaries"]["PrescribedPotential"] = [potential]
        heldout["Solver"]["Electrostatic"]["ResponseMatrix"] = False
        heldout_path = output / f"{heldout_name}.json"
        heldout_path.write_text(json.dumps(heldout, indent=2) + "\n")
        print(heldout_path)
    print(output / "basis_points.csv")
    print(
        write_library(
            output,
            args.library_name,
            args.separation,
            args.separation_tolerance,
            args.radius,
            args.coupon_depth,
            args.different_conductors,
            args.strip,
            open_contour_paths,
            args.metal_thickness,
            args.substrate_permittivity,
            interface_layers,
        )
    )


if __name__ == "__main__":
    main()
