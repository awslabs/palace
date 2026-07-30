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


def metal_band_cutoff(points, radius, metal_thickness):
    """Smoothly suppress a trace on both thin and fabricated PEC cuts."""
    transition = radius / 3.0
    distance = np.maximum(-points[:, 2], points[:, 2] - metal_thickness)
    coordinate = np.clip(distance / transition, 0.0, 1.0)
    return coordinate * coordinate * (3.0 - 2.0 * coordinate)


def heldout_potential(points, radius, metal_thickness):
    x = points[:, 0] / radius
    y = points[:, 1] / radius
    z = points[:, 2] / radius
    potential = (
        0.35
        + 0.20 * x
        - 0.15 * y
        + 0.10 * z
        + 0.08 * x * y
        + 0.06 * z * z
    )

    # Both coupons have conductor cuts in the matching surface at z = 0, while
    # the fabricated cut extends to z = metal_thickness. Use one smooth trace
    # which is compatible with both cuts. This avoids an order-dependent
    # Dirichlet jump where the matching and grounded boundaries meet.
    return metal_band_cutoff(points, radius, metal_thickness) * potential


def convergence_probe_potentials(points, radius, metal_thickness):
    """Return a small low-order trace space for response convergence tests."""
    x = points[:, 0] / radius
    y = points[:, 1] / radius
    z = points[:, 2] / radius
    cutoff = metal_band_cutoff(points, radius, metal_thickness)
    probes = (
        ("common", "cutoff", np.ones_like(x)),
        ("x-linear", "cutoff*x/R", x),
        ("y-linear", "cutoff*y/R", y),
        ("z-linear", "cutoff*z/R", z),
        ("xy-mixed", "cutoff*x*y/R^2", x * y),
        ("z-quadratic", "cutoff*z^2/R^2", z * z),
    )
    values = np.column_stack([cutoff * probe[2] for probe in probes])
    if np.linalg.matrix_rank(values) != len(probes):
        raise ValueError("corner convergence probes are linearly dependent")
    return [
        (name, expression, values[:, index])
        for index, (name, expression, _) in enumerate(probes)
    ]


def corner_frame(angle_degrees):
    angle = np.deg2rad(angle_degrees)
    if not 0.0 < angle < np.pi:
        raise ValueError("corner angle must lie strictly between zero and 180 degrees")
    first_normal = np.array([0.0, 1.0])
    second_normal = np.array([np.sin(angle), -np.cos(angle)])
    bisector = np.array([np.cos(0.5 * angle), np.sin(0.5 * angle)])
    return angle, first_normal, second_normal, bisector


def corner_center(angle_degrees, corner_radius):
    angle, _, _, bisector = corner_frame(angle_degrees)
    return corner_radius * bisector / np.sin(0.5 * angle)


def metal_footprint_mask(
    points,
    radius,
    angle_degrees,
    corner_radius,
    offset,
    topology,
):
    tolerance = 1.0e-12 * radius
    _, first_normal, second_normal, _ = corner_frame(angle_degrees)
    coordinates = points[:, :2]
    first_distance = coordinates @ first_normal
    second_distance = coordinates @ second_normal
    in_wedge = (first_distance >= offset - tolerance) & (
        second_distance >= offset - tolerance
    )
    if corner_radius > 0.0:
        arc_radius = corner_radius - offset
        if arc_radius <= 0.0:
            raise ValueError(
                "corner offset must be smaller than the plan-view corner radius"
            )
        center = corner_center(angle_degrees, corner_radius)
        rounded_wedge = (
            (first_distance >= corner_radius - tolerance)
            | (second_distance >= corner_radius - tolerance)
            | (
                np.sum((coordinates - center) ** 2, axis=1)
                <= arc_radius**2 + tolerance**2
            )
        )
        in_wedge &= rounded_wedge
    if topology == "convex":
        in_metal = in_wedge
    else:
        on_corner_boundary = (
            (np.abs(first_distance - offset) <= tolerance)
            & (second_distance >= corner_radius - tolerance)
        ) | (
            (np.abs(second_distance - offset) <= tolerance)
            & (first_distance >= corner_radius - tolerance)
        )
        in_metal = ~in_wedge | on_corner_boundary
    return in_metal


def thin_metal_mask(points, radius, angle_degrees, corner_radius, topology):
    tolerance = 1.0e-12 * radius
    return (np.abs(points[:, 2]) <= tolerance) & metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        0.0,
        topology,
    )


def pec_trace_mask(
    points,
    radius,
    angle_degrees,
    corner_radius,
    metal_thickness,
    sidewall_angle,
    topology,
):
    tolerance = 1.0e-12 * radius
    bottom = (np.abs(points[:, 2]) <= tolerance) & metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        0.0,
        topology,
    )
    pullback = metal_thickness / np.tan(np.deg2rad(sidewall_angle))
    top_offset = pullback if topology == "convex" else -pullback
    top_footprint = metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        top_offset,
        topology,
    )
    swept_footprint = top_footprint | metal_footprint_mask(
        points,
        radius,
        angle_degrees,
        corner_radius,
        0.0,
        topology,
    )
    top = (
        np.abs(points[:, 2] - metal_thickness) <= tolerance
    ) & swept_footprint
    return bottom | top


def dielectric(
    index, attributes, interface_type, radius, thickness, permittivity
):
    _, loss_tangent = INTERFACES[interface_type]
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": thickness,
        "Permittivity": permittivity,
        "LossTan": loss_tangent,
        "LocalizeEdgeEnergy": True,
        "SaveLocalEdgeEnergy": False,
        "AutomaticEdges": True,
        "EdgeExcludeAttributes": [1],
        "EdgeFrameNormal": [0.0, 0.0, 1.0],
        "EdgeDistances": [radius],
    }


def make_config(
    output,
    name,
    mesh,
    traces,
    radius,
    order,
    fabricated,
    substrate_permittivity,
    interface_layers,
):
    interfaces = (
        [
            dielectric(1, [3], "SA", radius, *interface_layers["SA"]),
            dielectric(2, [2], "MS", radius, *interface_layers["MS"]),
            dielectric(3, [4], "MA", radius, *interface_layers["MA"]),
        ]
        if fabricated
        else [
            dielectric(1, [3], "SA", radius, *interface_layers["SA"]),
            dielectric(2, [2], "MS", radius, *interface_layers["MS"]),
            dielectric(3, [2], "MA", radius, *interface_layers["MA"]),
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
                # These generated coupon solves never adapt the mesh.
                "EstimatorTol": 5.0e-1,
                "EstimatorMaxIts": 5,
                "EstimatorMG": True,
            },
        },
    }


def write_library(
    output,
    radius,
    angle_degrees,
    corner_radius,
    contour_groups,
    zero_trace_indices,
    topology,
    metal_thickness,
    overetch_depth,
    sidewall_angle,
    top_rounding,
    trench_rounding,
    substrate_permittivity,
    interface_layers,
):
    topology_name = f"{topology.capitalize()}Corner"
    model_name = f"{topology}-corner-{angle_degrees:g}deg"
    if corner_radius > 0.0:
        model_name += f"-r{corner_radius:g}um"
    reference = (
        [*corner_center(angle_degrees, corner_radius), 0.0]
        if topology == "convex" and corner_radius > 0.0
        else [0.0, 0.0, 0.0]
    )
    corner_radius_tolerance = max(
        0.02 * corner_radius,
        1.0e-3 * radius,
    )
    model = {
        "Name": model_name,
        "Topology": topology_name,
        "Angle": angle_degrees,
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
    if zero_trace_indices:
        model["ZeroTraceIndices"] = zero_trace_indices
    library = {
        "Version": 3,
        "Name": (
            "100nm-metal-50nm-overetch-10nm-rounding-"
            f"{topology}-corner-r{corner_radius:g}um-prototype"
        ),
        "MatchingRadius": radius,
        "Fabrication": {
            "LengthUnit": "um",
            "MetalThickness": metal_thickness,
            "OveretchDepth": overetch_depth,
            "SidewallAngle": sidewall_angle,
            "TopRounding": top_rounding,
            "TrenchRounding": trench_rounding,
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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path, required=True)
    parser.add_argument("--fabricated-mesh", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--angle", type=float, default=90.0)
    parser.add_argument("--corner-radius", type=float, default=0.0)
    parser.add_argument("--ring-size", type=int, default=8)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--overetch-depth", type=float, default=0.05)
    parser.add_argument("--sidewall-angle", type=float, default=80.0)
    parser.add_argument("--top-rounding", type=float, default=0.01)
    parser.add_argument("--trench-rounding", type=float, default=0.01)
    parser.add_argument("--substrate-permittivity", type=float, default=11.47)
    parser.add_argument("--sa-thickness", type=float, default=0.002)
    parser.add_argument("--sa-permittivity", type=float, default=4.0)
    parser.add_argument("--ms-thickness", type=float, default=0.002)
    parser.add_argument("--ms-permittivity", type=float, default=11.47)
    parser.add_argument("--ma-thickness", type=float, default=0.002)
    parser.add_argument("--ma-permittivity", type=float, default=10.0)
    parser.add_argument(
        "--topology", choices=("convex", "concave"), default="convex"
    )
    args = parser.parse_args()
    if args.radius <= 0.0:
        parser.error("--radius must be positive")
    if not 0.0 < args.angle < 180.0:
        parser.error("--angle must lie strictly between zero and 180 degrees")
    if not 0.0 <= args.corner_radius < args.radius:
        parser.error("--corner-radius must lie in [0, radius)")
    tangent_distance = args.corner_radius / np.tan(
        0.5 * np.deg2rad(args.angle)
    )
    if args.corner_radius > 0.0 and tangent_distance >= args.radius:
        parser.error(
            "rounded-corner tangency points must lie inside the matching box"
        )
    if args.order < 1:
        parser.error("--order must be positive")
    if not 0.0 < args.metal_thickness < args.radius:
        parser.error("--metal-thickness must lie between zero and the radius")
    if not 0.0 <= args.overetch_depth < args.radius:
        parser.error("--overetch-depth must be nonnegative and smaller than the radius")
    if not 0.0 < args.sidewall_angle <= 90.0:
        parser.error("--sidewall-angle must lie in (0, 90]")
    if not 0.0 <= args.top_rounding < args.metal_thickness:
        parser.error(
            "--top-rounding must be nonnegative and smaller than the metal thickness"
        )
    if not 0.0 <= args.trench_rounding <= args.overetch_depth:
        parser.error(
            "--trench-rounding must lie between zero and the overetch depth"
        )
    pullback = args.metal_thickness / np.tan(
        np.deg2rad(args.sidewall_angle)
    )
    if (
        args.corner_radius > 0.0
        and args.topology == "convex"
        and pullback >= args.corner_radius
    ):
        parser.error(
            "--sidewall-angle gives a pullback larger than the corner radius"
        )
    trench_pullback = args.overetch_depth / np.tan(
        np.deg2rad(args.sidewall_angle)
    )
    if (
        args.corner_radius > 0.0
        and args.topology == "concave"
        and trench_pullback >= args.corner_radius
    ):
        parser.error(
            "--sidewall-angle gives a trench pullback larger than the corner radius"
        )
    material_values = (
        args.substrate_permittivity,
        args.sa_thickness,
        args.sa_permittivity,
        args.ms_thickness,
        args.ms_permittivity,
        args.ma_thickness,
        args.ma_permittivity,
    )
    if any(not np.isfinite(value) or value <= 0.0 for value in material_values):
        parser.error("substrate and interface-layer properties must be finite and positive")
    interface_layers = {
        "SA": (args.sa_thickness, args.sa_permittivity),
        "MS": (args.ms_thickness, args.ms_permittivity),
        "MA": (args.ma_thickness, args.ma_permittivity),
    }

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
    zero_trace_indices = (
        np.flatnonzero(
            pec_trace_mask(
                points,
                args.radius,
                args.angle,
                args.corner_radius,
                args.metal_thickness,
                args.sidewall_angle,
                args.topology,
            )
        )
        + 1
    ).tolist()
    if not zero_trace_indices:
        raise ValueError("corner coupon has no PEC-constrained trace knots")
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
            args.substrate_permittivity,
            interface_layers,
        )
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")
    library = write_library(
        output,
        args.radius,
        args.angle,
        args.corner_radius,
        contour_groups,
        zero_trace_indices,
        args.topology,
        args.metal_thickness,
        args.overetch_depth,
        args.sidewall_angle,
        args.top_rounding,
        args.trench_rounding,
        args.substrate_permittivity,
        interface_layers,
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
        heldout_potential(
            fine_points,
            args.radius,
            args.metal_thickness,
        ),
    )
    heldout_coefficients = heldout_potential(
        points,
        args.radius,
        args.metal_thickness,
    )
    np.savetxt(
        output / "heldout-coefficients.csv",
        heldout_coefficients,
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
            args.substrate_permittivity,
            interface_layers,
        )
        config["Solver"]["Electrostatic"]["ResponseMatrix"] = False
        config["Solver"]["Electrostatic"]["AggregateResponseMatrix"] = False
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")

    probe_directory = output / "convergence-probes"
    probe_directory.mkdir(parents=True, exist_ok=True)
    for stale_probe in probe_directory.glob("probe-*.csv"):
        stale_probe.unlink()
    convergence_probes = convergence_probe_potentials(
        fine_points,
        args.radius,
        args.metal_thickness,
    )
    probe_paths = []
    probe_metadata = []
    for index, (name, expression, values) in enumerate(
        convergence_probes, start=1
    ):
        path = probe_directory / f"probe-{index:02d}-{name}.csv"
        write_surface_trace(path, fine_points, fine_triangles, values)
        probe_paths.append(path)
        probe_metadata.append(
            {
                "Index": index,
                "Name": name,
                "Expression": expression,
            }
        )
    (output / "probe-manifest.json").write_text(
        json.dumps({"Version": 1, "Probes": probe_metadata}, indent=2) + "\n"
    )
    for name, mesh, fabricated in (
        ("probe-thin", thin_mesh, False),
        ("probe-fabricated", fabricated_mesh, True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            probe_paths,
            args.radius,
            args.order,
            fabricated,
            args.substrate_permittivity,
            interface_layers,
        )
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")

    print(f"Generated {len(points)} spatial basis traces")
    print(f"Generated {len(convergence_probes)} convergence probe traces")
    print(f"ContourGroups: {contour_groups}")
    if zero_trace_indices:
        print(f"ZeroTraceIndices: {zero_trace_indices}")
    print(output / "thin.json")
    print(output / "fabricated.json")
    print(library)
    print(output / "heldout-thin.json")
    print(output / "heldout-fabricated.json")
    print(output / "probe-thin.json")
    print(output / "probe-fabricated.json")


if __name__ == "__main__":
    main()
