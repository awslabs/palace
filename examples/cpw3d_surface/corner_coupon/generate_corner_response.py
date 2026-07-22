#!/usr/bin/env python3

"""Generate a spatial trace basis and Palace configs for a 3D corner coupon."""

import argparse
import json
from pathlib import Path

import numpy as np


INTERFACES = {
    "SA": (4.0, 2.0e-3),
    "MS": (11.47, 3.0e-4),
    "MA": (10.0, 3.0e-2),
}


def square_ring(half_width, z, size):
    if size < 8 or size % 8:
        raise ValueError("ring size must be a multiple of eight and at least eight")
    perimeter_coordinate = np.arange(size, dtype=float) * (8.0 * half_width / size)
    # Start at (-half_width, 0) so the Maxwell anchor path lies in the gap
    # instead of crossing the metal quadrant.
    perimeter_coordinate = (perimeter_coordinate + 7.0 * half_width) % (
        8.0 * half_width
    )
    points = []
    for coordinate in perimeter_coordinate:
        if coordinate < 2.0 * half_width:
            x = -half_width + coordinate
            y = -half_width
        elif coordinate < 4.0 * half_width:
            x = half_width
            y = -half_width + coordinate - 2.0 * half_width
        elif coordinate < 6.0 * half_width:
            x = half_width - coordinate + 4.0 * half_width
            y = half_width
        else:
            x = -half_width
            y = half_width - coordinate + 6.0 * half_width
        points.append((x, y, z))
    points = np.asarray(points)
    start = np.argmin(np.linalg.norm(points - np.array([-half_width, 0.0, z]), axis=1))
    return np.roll(points, -start, axis=0)


def connect_rings(triangles, first, second, size):
    for index in range(size):
        next_index = (index + 1) % size
        triangles.append((first + index, first + next_index, second + next_index))
        triangles.append((first + index, second + next_index, second + index))


def cap_ring(triangles, offset, size, reverse):
    for index in range(1, size - 1):
        triangle = (offset, offset + index, offset + index + 1)
        triangles.append(tuple(reversed(triangle)) if reverse else triangle)


def build_surface(
    radius,
    ring_size,
    metal_thickness,
    overetch_depth,
    cap_centers=False,
):
    # Keep the trace triangulation conforming to every fabrication plane that
    # reaches the matching surface. The coupon mesh resolves these intersections
    # so the narrow trace hats across the process zone have active boundary DOFs.
    levels = sorted(
        {
            -radius,
            -radius / 3.0,
            -overetch_depth,
            0.0,
            metal_thickness,
            radius / 3.0,
            radius,
        }
    )
    rings = [square_ring(radius, level, ring_size) for level in levels]
    rings.append(square_ring(radius / 3.0, radius, ring_size))
    rings.append(square_ring(radius / 3.0, -radius, ring_size))
    points = np.vstack(rings)

    triangles = []
    for ring in range(len(levels) - 1):
        connect_rings(
            triangles, ring * ring_size, (ring + 1) * ring_size, ring_size
        )
    top_inner = len(levels) * ring_size
    bottom_inner = (len(levels) + 1) * ring_size
    connect_rings(triangles, (len(levels) - 1) * ring_size, top_inner, ring_size)
    connect_rings(triangles, bottom_inner, 0, ring_size)
    if cap_centers:
        top_center = len(points)
        bottom_center = top_center + 1
        points = np.vstack(
            (points, (0.0, 0.0, radius), (0.0, 0.0, -radius))
        )
        for index in range(ring_size):
            next_index = (index + 1) % ring_size
            triangles.append(
                (top_inner + index, top_inner + next_index, top_center)
            )
            triangles.append(
                (bottom_inner + next_index, bottom_inner + index, bottom_center)
            )
    else:
        cap_ring(triangles, top_inner, ring_size, False)
        cap_ring(triangles, bottom_inner, ring_size, True)
    return points, np.asarray(triangles, dtype=int), [ring_size] * len(rings)


def write_basis(output, points, triangles):
    np.savetxt(
        output / "basis-points.csv",
        points,
        delimiter=",",
        header="x,y,z",
        comments="",
        fmt="%.16e",
    )
    trace_directory = output / "traces"
    trace_directory.mkdir(parents=True, exist_ok=True)
    for stale_trace in trace_directory.glob("basis-*.csv"):
        stale_trace.unlink()
    paths = []
    for basis in range(len(points)):
        path = trace_directory / f"basis-{basis + 1:03d}.csv"
        rows = []
        for triangle_index, triangle in enumerate(triangles, start=1):
            for vertex in triangle:
                rows.append(
                    (
                        *points[vertex],
                        1.0 if vertex == basis else 0.0,
                        triangle_index,
                    )
                )
        np.savetxt(
            path,
            np.asarray(rows),
            delimiter=",",
            header="x,y,z,V,triangle",
            comments="",
            fmt=("%.16e", "%.16e", "%.16e", "%.16e", "%d"),
        )
        paths.append(path)
    return paths


def write_surface_trace(path, points, triangles, values):
    rows = []
    for triangle_index, triangle in enumerate(triangles, start=1):
        for vertex in triangle:
            rows.append((*points[vertex], values[vertex], triangle_index))
    np.savetxt(
        path,
        np.asarray(rows),
        delimiter=",",
        header="x,y,z,V,triangle",
        comments="",
        fmt=("%.16e", "%.16e", "%.16e", "%.16e", "%d"),
    )


def heldout_potential(points, radius):
    x = points[:, 0] / radius
    y = points[:, 1] / radius
    z = points[:, 2] / radius
    return 0.35 + 0.20 * x - 0.15 * y + 0.10 * z + 0.08 * x * y + 0.06 * z * z


def dielectric(index, attributes, interface_type, radius):
    permittivity, loss_tangent = INTERFACES[interface_type]
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": 2.0e-3,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "LocalizeEdgeEnergy": True,
        "AutomaticEdges": True,
        "EdgeExcludeAttributes": [1],
        "EdgeFrameNormal": [0.0, 0.0, 1.0],
        "EdgeDistances": [radius],
    }


def make_config(output, name, mesh, traces, radius, order, fabricated):
    interfaces = (
        [
            dielectric(1, [3], "SA", radius),
            dielectric(2, [2], "MS", radius),
            dielectric(3, [4], "MA", radius),
        ]
        if fabricated
        else [
            dielectric(1, [3], "SA", radius),
            dielectric(2, [2], "MS", radius),
            dielectric(3, [2], "MA", radius),
        ]
    )
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
            "OutputFormats": {"Paraview": False, "GridFunction": False},
        },
        "Model": {
            "Mesh": str(mesh),
            "L0": 1.0e-6,
            "Refinement": {"MaxIts": 0},
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
        "Boundaries": {
            "Ground": {"Attributes": [2, 4] if fabricated else [2]},
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
            "Electrostatic": {"Save": 0, "ResponseMatrix": True},
            "Linear": {
                "Type": "BoomerAMG",
                "KSPType": "CG",
                "Tol": 1.0e-10,
                "MaxIts": 1000,
                "EstimatorTol": 1.0e-1,
                "EstimatorMaxIts": 20,
                "EstimatorMG": True,
            },
        },
    }


def write_library(output, radius, corner_radius, contour_groups, topology):
    topology_name = f"{topology.capitalize()}Corner"
    model_name = f"{topology}-corner-90deg"
    if corner_radius > 0.0:
        model_name += f"-r{corner_radius:g}um"
    reference = (
        [corner_radius, corner_radius, 0.0]
        if topology == "convex" and corner_radius > 0.0
        else [0.0, 0.0, 0.0]
    )
    corner_radius_tolerance = max(
        0.02 * corner_radius,
        1.0e-3 * radius,
    )
    library = {
        "Version": 1,
        "Name": (
            "100nm-metal-50nm-overetch-10nm-rounding-"
            f"{topology}-corner-r{corner_radius:g}um-prototype"
        ),
        "MatchingRadius": radius,
        "Models": [
            {
                "Name": model_name,
                "Topology": topology_name,
                "Angle": 90.0,
                "AngleTolerance": 2.0,
                "CornerRadius": corner_radius,
                "CornerRadiusTolerance": corner_radius_tolerance,
                "Reference": reference,
                "FabricatedMatrix": "postpro/fabricated/domain-response-matrix.csv",
                "ThinMatrix": "postpro/thin/domain-response-matrix.csv",
                "FabricatedSurfaceMatrix":
                    "postpro/fabricated/surface-response-matrix-aggregate.csv",
                "ThinSurfaceMatrix":
                    "postpro/thin/surface-response-matrix-aggregate.csv",
                "BasisPoints": "basis-points.csv",
                "ContourGroups": contour_groups,
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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path, required=True)
    parser.add_argument("--fabricated-mesh", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--corner-radius", type=float, default=0.0)
    parser.add_argument("--ring-size", type=int, default=8)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--overetch-depth", type=float, default=0.05)
    parser.add_argument(
        "--topology", choices=("convex", "concave"), default="convex"
    )
    args = parser.parse_args()
    if args.radius <= 0.0:
        parser.error("--radius must be positive")
    if not 0.0 <= args.corner_radius < args.radius:
        parser.error("--corner-radius must lie in [0, radius)")
    if args.order < 1:
        parser.error("--order must be positive")
    if not 0.0 < args.metal_thickness < args.radius:
        parser.error("--metal-thickness must lie between zero and the radius")
    if not 0.0 <= args.overetch_depth < args.radius:
        parser.error("--overetch-depth must be nonnegative and smaller than the radius")

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    thin_mesh = args.thin_mesh.expanduser().resolve()
    fabricated_mesh = args.fabricated_mesh.expanduser().resolve()
    for mesh in (thin_mesh, fabricated_mesh):
        if not mesh.is_file():
            raise FileNotFoundError(mesh)

    points, triangles, contour_groups = build_surface(
        args.radius,
        args.ring_size,
        args.metal_thickness,
        args.overetch_depth,
    )
    traces = write_basis(output, points, triangles)
    for name, mesh, fabricated in (
        ("thin", thin_mesh, False),
        ("fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            traces,
            args.radius,
            args.order,
            fabricated,
        )
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")
    library = write_library(
        output,
        args.radius,
        args.corner_radius,
        contour_groups,
        args.topology,
    )

    fine_points, fine_triangles, _ = build_surface(
        args.radius,
        max(16, 4 * args.ring_size),
        args.metal_thickness,
        args.overetch_depth,
        cap_centers=True,
    )
    heldout_trace = output / "heldout-trace.csv"
    write_surface_trace(
        heldout_trace,
        fine_points,
        fine_triangles,
        heldout_potential(fine_points, args.radius),
    )
    np.savetxt(
        output / "heldout-coefficients.csv",
        heldout_potential(points, args.radius),
        delimiter=",",
        header="coefficient_V",
        comments="",
        fmt="%.16e",
    )
    for name, mesh, fabricated in (
        ("heldout-thin", thin_mesh, False),
        ("heldout-fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            [heldout_trace],
            args.radius,
            args.order,
            fabricated,
        )
        config["Solver"]["Electrostatic"]["ResponseMatrix"] = False
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")

    print(f"Generated {len(points)} spatial basis traces")
    print(f"ContourGroups: {contour_groups}")
    print(output / "thin.json")
    print(output / "fabricated.json")
    print(library)
    print(output / "heldout-thin.json")
    print(output / "heldout-fabricated.json")


if __name__ == "__main__":
    main()
