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


def write_bases(output, cutout_width, radius, basis_size, samples):
    if basis_size % 2:
        raise ValueError("The paired-edge contour basis size must be even")
    half_width = 0.5 * cutout_width + radius
    perimeter = 4.0 * (half_width + radius)
    distances = np.linspace(0.0, perimeter, samples, endpoint=False)
    points = np.asarray(
        [contour_point(distance, half_width, radius) for distance in distances]
    )
    # Align a knot with each point where the matching contour meets a grounded
    # metal region. The two intersections are half a perimeter apart, so an
    # even periodic basis resolves both. Missing these zero-potential junctions
    # creates a spurious singular response in an otherwise smooth trace fit.
    spacing = perimeter / basis_size
    right_metal_intersection = 2.0 * half_width + radius
    offset = right_metal_intersection % spacing
    knot_distances = offset + spacing * np.arange(basis_size)
    knots = np.asarray(
        [contour_point(distance, half_width, radius) for distance in knot_distances]
    )
    np.savetxt(
        output / "basis_points.csv",
        knots,
        delimiter=",",
        header="x,y,z",
        comments="",
        fmt="%.16e",
    )

    phase = basis_size * (distances - offset) / perimeter
    paths = []
    for index in range(basis_size):
        delta = np.abs(
            (phase - index + 0.5 * basis_size) % basis_size - 0.5 * basis_size
        )
        values = np.maximum(1.0 - delta, 0.0)
        path = output / f"basis_hat{index + 1:03d}.csv"
        np.savetxt(
            path,
            np.column_stack((points, values)),
            delimiter=",",
            header="x,y,z,V",
            comments="",
            fmt="%.16e",
        )
        paths.append(path)
    return paths


def dielectric(index, attributes, interface_type):
    permittivity, loss_tangent = INTERFACE_PROPERTIES[interface_type]
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": 0.002,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "EdgeAttributes": [2],
        "EdgeExcludeAttributes": [1],
        "EdgeDistances": [0.2],
        "LocalizeEdgeEnergy": True,
        "EdgeFrameNormal": [0.0, 1.0, 0.0],
    }


def make_config(output, name, mesh, traces, fabricated):
    if fabricated:
        ground = [2, 4, 5]
        interfaces = [
            dielectric(1, [3, 6], "SA"),
            dielectric(2, [2], "MS"),
            dielectric(3, [4, 5], "MA"),
        ]
    else:
        ground = [2]
        interfaces = [
            dielectric(1, [3], "SA"),
            dielectric(2, [2], "MS"),
            dielectric(3, [2], "MA"),
        ]
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
        },
        "Model": {"Mesh": str(mesh), "L0": 1.0e-6, "Lc": 1055.0},
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
            "Order": 2,
            "Device": "CPU",
            "Electrostatic": {"Save": 0, "ResponseMatrix": True},
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cutout-width", type=float, required=True)
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--basis-size", type=int, default=96)
    parser.add_argument("--samples", type=int, default=1200)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path)
    parser.add_argument("--fabricated-mesh", type=Path)
    args = parser.parse_args()

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
    traces = write_bases(
        output,
        args.cutout_width,
        args.radius,
        args.basis_size,
        args.samples,
    )
    for name, mesh, fabricated in (
        ("edge_pair_thin", thin_mesh, False),
        ("edge_pair_fabricated", fabricated_mesh, True),
    ):
        config = make_config(output, name, mesh, traces, fabricated)
        path = output / f"{name}.json"
        path.write_text(json.dumps(config, indent=2) + "\n")
        print(path)
    print(output / "basis_points.csv")


if __name__ == "__main__":
    main()
