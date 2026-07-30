#!/usr/bin/env python3

"""Generate response inputs for an arbitrary parallel multi-edge coupon."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


INTERFACE_PROPERTIES = {
    "SA": (4.0, 0.002),
    "MS": (11.47, 0.0003),
    "MA": (10.0, 0.03),
}


def load_edges(path):
    data = json.loads(path.read_text())
    edges = data.get("Edges")
    if not isinstance(edges, list) or len(edges) < 3:
        raise ValueError("A parallel cluster specification needs at least three Edges")
    offsets = np.asarray([float(edge["Offset"]) for edge in edges])
    directions = np.asarray([int(edge["GapDirection"]) for edge in edges])
    conductors = np.asarray([int(edge["Conductor"]) for edge in edges])
    if (
        offsets[0] != 0.0
        or np.any(np.diff(offsets) <= 0.0)
        or np.any(np.abs(directions) != 1)
    ):
        raise ValueError("Invalid parallel edge offsets or gap directions")
    labels = []
    for conductor in conductors:
        if conductor not in labels:
            labels.append(conductor)
        if conductor != labels.index(conductor) + 1:
            raise ValueError("Conductor labels must be canonical")
    return data, offsets, directions, conductors


def metal_intervals(offsets, directions, conductors, radius):
    bounds = np.concatenate(([offsets[0] - radius], offsets, [offsets[-1] + radius]))
    intervals = []
    for index, (left, right) in enumerate(zip(bounds, bounds[1:])):
        if index == 0:
            metal = directions[0] == 1
            conductor = conductors[0]
        elif index == len(bounds) - 2:
            metal = directions[-1] == -1
            conductor = conductors[-1]
        else:
            metal = directions[index - 1] == -1
            if metal != (directions[index] == 1):
                raise ValueError(
                    "Adjacent gap directions do not define one material interval"
                )
            conductor = conductors[index]
            if metal and conductors[index - 1] != conductor:
                raise ValueError(
                    "Both boundaries of a metal interval must use one conductor"
                )
        if metal:
            intervals.append((float(left), float(right), int(conductor)))
    return intervals


def contour_point(distance, xmin, xmax, radius):
    width = xmax - xmin
    perimeter = 2.0 * (width + 2.0 * radius)
    distance %= perimeter
    if distance < width:
        return np.array([xmin + distance, -radius, 0.0])
    distance -= width
    if distance < 2.0 * radius:
        return np.array([xmax, -radius + distance, 0.0])
    distance -= 2.0 * radius
    if distance < width:
        return np.array([xmax - distance, radius, 0.0])
    distance -= width
    return np.array([xmin, radius - distance, 0.0])


def write_bases(
    output,
    offsets,
    directions,
    conductors,
    radius,
    metal_thickness,
    basis_size,
    samples,
):
    xmin = offsets[0] - radius
    xmax = offsets[-1] + radius
    width = xmax - xmin
    perimeter = 2.0 * (width + 2.0 * radius)
    spacing = perimeter / basis_size
    right_lower = width + radius
    right_upper = right_lower + metal_thickness
    left_lower = 2.0 * width + 3.0 * radius
    left_upper = left_lower - metal_thickness
    cuts = []
    if directions[0] == 1:
        cuts.append((left_upper, left_lower, int(conductors[0]), xmin))
    if directions[-1] == -1:
        cuts.append((right_lower, right_upper, int(conductors[-1]), xmax))

    offset = cuts[0][1] % spacing if cuts else 0.0
    uniform = (offset + spacing * np.arange(basis_size)) % perimeter
    endpoints = [value for cut in cuts for value in cut[:2]]
    knot_distances = np.unique(np.concatenate((uniform, endpoints)))
    knot_distances.sort()
    tolerance = 1.0e-12 * perimeter

    def cut_conductor(distance):
        for begin, end, conductor, _ in cuts:
            if begin - tolerance <= distance <= end + tolerance:
                return conductor
        return 0

    constrained_conductors = np.asarray(
        [cut_conductor(distance) for distance in knot_distances]
    )
    free_indices = np.flatnonzero(constrained_conductors == 0)
    distances = np.unique(
        np.concatenate(
            (
                np.linspace(0.0, perimeter, samples, endpoint=False),
                knot_distances,
            )
        )
    )
    points = np.asarray(
        [contour_point(distance, xmin, xmax, radius) for distance in distances]
    )
    knots = np.asarray(
        [contour_point(distance, xmin, xmax, radius) for distance in knot_distances]
    )
    np.savetxt(
        output / "basis_points.csv",
        knots[free_indices],
        delimiter=",",
        header="x,y,z",
        comments="",
        fmt="%.16e",
    )

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

    traces = []
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
        traces.append(path)

    conductor_count = int(max(conductors))
    conductor_traces = []
    for conductor in range(2, conductor_count + 1):
        values = np.zeros(len(distances))
        for index in np.flatnonzero(constrained_conductors == conductor):
            values += hat_values(index)
        path = output / f"basis_conductor_{conductor:03d}.csv"
        np.savetxt(
            path,
            np.column_stack((points, values)),
            delimiter=",",
            header="x,y,z,V",
            comments="",
            fmt="%.16e",
        )
        conductor_traces.append(path)

    intervals = metal_intervals(offsets, directions, conductors, radius)
    references = []
    portals = []
    for conductor in range(1, conductor_count + 1):
        candidates = [
            (right - left, 0.5 * (left + right))
            for left, right, label in intervals
            if label == conductor
        ]
        if not candidates:
            raise ValueError(f"Conductor {conductor} has no metal interval")
        _, reference_x = max(candidates)
        references.append([reference_x, 0.0, 0.0])
        edge = int(np.flatnonzero(conductors == conductor)[0])
        portal = width + 2.0 * radius + (xmax - offsets[edge])
        portals.append((portal % perimeter, conductor))
    portals.sort()
    output_indices = {
        knot_index: output_index
        for output_index, knot_index in enumerate(free_indices, start=1)
    }
    open_paths = []
    for index, (begin, conductor) in enumerate(portals):
        end, next_conductor = portals[(index + 1) % len(portals)]
        span = (end - begin) % perimeter
        arc = [
            knot_index
            for knot_index in free_indices
            if 0.0 <= (knot_distances[knot_index] - begin) % perimeter < span
        ]
        arc.sort(key=lambda knot_index: (knot_distances[knot_index] - begin) % perimeter)
        if not arc:
            raise ValueError("Trace basis is too coarse to connect conductor portals")
        open_paths.append(
            {
                "Indices": [output_indices[knot_index] for knot_index in arc],
                "StartConductor": conductor,
                "EndConductor": next_conductor,
            }
        )
    return (
        traces,
        conductor_traces,
        references,
        open_paths,
        points,
        cuts,
    )


def write_heldout(
    output,
    traces,
    conductor_traces,
    points,
    cuts,
    radius,
    metal_thickness,
):
    basis_points = np.atleast_2d(
        np.loadtxt(output / "basis_points.csv", delimiter=",", skiprows=1)
    )

    def free_potential(coordinates, polynomial):
        x = coordinates[:, 0] / radius
        y = coordinates[:, 1] / radius
        if cuts:
            distances = []
            for _, _, _, boundary_x in cuts:
                vertical = np.maximum.reduce(
                    (
                        -coordinates[:, 1],
                        coordinates[:, 1] - metal_thickness,
                        np.zeros(len(coordinates)),
                    )
                )
                distances.append(np.hypot(coordinates[:, 0] - boundary_x, vertical))
            distance = np.min(np.column_stack(distances), axis=1)
            coordinate = np.clip(distance / (radius / 3.0), 0.0, 1.0)
            cutoff = coordinate * coordinate * (3.0 - 2.0 * coordinate)
        else:
            cutoff = np.ones(len(coordinates))
        return cutoff * polynomial(x, y)

    heldout_polynomial = lambda x, y: (
            0.35
            + 0.20 * x
            - 0.15 * y
            + 0.08 * x * y
            + 0.06 * y * y
    )

    coefficients = list(free_potential(basis_points, heldout_polynomial))
    values = free_potential(points, heldout_polynomial)
    for path in conductor_traces:
        coefficient = 1.0
        coefficients.append(coefficient)
        trace = np.atleast_2d(np.loadtxt(path, delimiter=",", skiprows=1))
        values += coefficient * trace[:, 3]
    trace_path = output / "heldout_trace.csv"
    np.savetxt(
        trace_path,
        np.column_stack((points, values)),
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

    probe_definitions = (
        ("common", "cutoff", lambda x, y: np.ones_like(x)),
        ("x-linear", "cutoff*x/R", lambda x, y: x),
        ("y-linear", "cutoff*y/R", lambda x, y: y),
        ("xy-mixed", "cutoff*x*y/R^2", lambda x, y: x * y),
        ("x-quadratic", "cutoff*x^2/R^2", lambda x, y: x * x),
        ("y-quadratic", "cutoff*y^2/R^2", lambda x, y: y * y),
    )
    probe_paths = []
    probe_manifest = []
    for index, (name, expression, polynomial) in enumerate(
        probe_definitions, start=1
    ):
        path = output / f"probe_{index:02d}.csv"
        np.savetxt(
            path,
            np.column_stack((points, free_potential(points, polynomial))),
            delimiter=",",
            header="x,y,z,V",
            comments="",
            fmt="%.16e",
        )
        probe_paths.append(path)
        probe_manifest.append({"Name": name, "Expression": expression})
    for conductor in range(2, len(conductor_traces) + 2):
        probe_manifest.append(
            {
                "Name": f"conductor-{conductor}",
                "Expression": f"V_{conductor}-V_1",
            }
        )
    (output / "probe-manifest.json").write_text(
        json.dumps({"Version": 1, "Probes": probe_manifest}, indent=2) + "\n"
    )
    return trace_path, probe_paths


def dielectric(index, attributes, interface_type, edge_attributes, thickness, permittivity):
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
    conductor_traces,
    fabricated,
    order,
    coupon_depth,
    substrate_permittivity,
    interface_layers,
    conductor_count,
):
    ms = [10 + conductor for conductor in range(1, conductor_count + 1)]
    ma = [100 + conductor for conductor in range(1, conductor_count + 1)]
    ground = ms + ma if fabricated else ms
    interfaces = [
        dielectric(1, [3], "SA", ms, *interface_layers["SA"]),
        dielectric(2, ms, "MS", ms, *interface_layers["MS"]),
        dielectric(3, ma if fabricated else ms, "MA", ms, *interface_layers["MA"]),
    ]
    potentials = [
        {"Index": index, "Attributes": [1], "DataFile": str(trace)}
        for index, trace in enumerate(traces, start=1)
    ]
    for conductor, trace in enumerate(conductor_traces, start=2):
        attributes = [10 + conductor]
        if fabricated:
            attributes.append(100 + conductor)
        potentials.append(
            {
                "Index": len(potentials) + 1,
                "Attributes": [1],
                "TerminalAttributes": attributes,
                "DataFile": str(trace),
            }
        )
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
            "OutputFormats": {"Paraview": False, "GridFunction": False},
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
                "EstimatorTol": 5.0e-1,
                "EstimatorMaxIts": 5,
                "EstimatorMG": True,
            },
        },
    }


def write_library(
    output,
    name,
    edges,
    references,
    open_paths,
    radius,
    coupon_depth,
    metal_thickness,
    overetch_depth,
    sidewall_angle,
    top_rounding,
    trench_rounding,
    substrate_permittivity,
    interface_layers,
    model_name,
    edge_offset_tolerance,
):
    model = {
        "Name": model_name,
        "Topology": "ParallelEdgeCluster",
        "EdgeOffsetTolerance": edge_offset_tolerance,
        "Edges": edges,
        "ConductorReferences": references,
        "OpenContourPaths": open_paths,
        "FabricatedMatrix": "postpro/cluster_fabricated/domain-response-matrix.csv",
        "ThinMatrix": "postpro/cluster_thin/domain-response-matrix.csv",
        "FabricatedSurfaceMatrix":
            "postpro/cluster_fabricated/surface-response-matrix.csv",
        "ThinSurfaceMatrix": "postpro/cluster_thin/surface-response-matrix.csv",
        "BasisPoints": "basis_points.csv",
        "Interfaces": [
            {"Type": "SA", "Coupon": 1},
            {"Type": "MS", "Coupon": 2},
            {"Type": "MA", "Coupon": 3},
        ],
        "CouponDepth": coupon_depth,
    }
    library = {
        "Version": 3,
        "Name": name,
        "MatchingRadius": radius,
        "Fabrication": {
            "LengthUnit": "um",
            "MetalThickness": metal_thickness,
            "OveretchDepth": overetch_depth,
            "SidewallAngleDegrees": sidewall_angle,
            "TopRoundingRadius": top_rounding,
            "BottomRoundingRadius": trench_rounding,
            "SubstratePermittivity": substrate_permittivity,
            "InterfaceLayers": {
                interface: {
                    "Thickness": values[0],
                    "Permittivity": values[1],
                }
                for interface, values in interface_layers.items()
            },
        },
        "Models": [model],
    }
    path = output / "process-library.json"
    path.write_text(json.dumps(library, indent=2) + "\n")
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("spec", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path, required=True)
    parser.add_argument("--fabricated-mesh", type=Path, required=True)
    parser.add_argument("--radius", type=float, default=2.0)
    parser.add_argument("--metal-thickness", type=float, default=0.1)
    parser.add_argument("--overetch-depth", type=float, default=0.05)
    parser.add_argument("--sidewall-angle", type=float, default=80.0)
    parser.add_argument("--top-rounding", type=float, default=0.01)
    parser.add_argument("--trench-rounding", type=float, default=0.01)
    parser.add_argument("--basis-size", type=int, default=96)
    parser.add_argument("--samples", type=int, default=1200)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--coupon-depth", type=float, default=1055.0)
    parser.add_argument("--edge-offset-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--library-name", required=True)
    parser.add_argument("--model-name")
    parser.add_argument("--substrate-permittivity", type=float, default=11.47)
    for interface, (permittivity, thickness) in INTERFACE_PROPERTIES.items():
        parser.add_argument(f"--{interface.lower()}-thickness", type=float, default=thickness)
        parser.add_argument(
            f"--{interface.lower()}-permittivity",
            type=float,
            default=permittivity,
        )
    args = parser.parse_args()
    if args.basis_size <= 0 or args.samples < args.basis_size:
        parser.error("invalid basis size or sample count")
    if args.edge_offset_tolerance < 0.0:
        parser.error("--edge-offset-tolerance must be nonnegative")

    spec_path = args.spec.expanduser().resolve()
    data, offsets, directions, conductors = load_edges(spec_path)
    if args.model_name:
        model_name = args.model_name
    else:
        signature = hashlib.sha256(
            json.dumps(data["Edges"], sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()[:12]
        model_name = f"parallel-edge-cluster-{signature}"
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    (
        traces,
        conductor_traces,
        references,
        open_paths,
        points,
        cuts,
    ) = write_bases(
        output,
        offsets,
        directions,
        conductors,
        args.radius,
        args.metal_thickness,
        args.basis_size,
        args.samples,
    )
    heldout_trace, probe_traces = write_heldout(
        output,
        traces,
        conductor_traces,
        points,
        cuts,
        args.radius,
        args.metal_thickness,
    )
    interface_layers = {
        interface: (
            getattr(args, f"{interface.lower()}_thickness"),
            getattr(args, f"{interface.lower()}_permittivity"),
        )
        for interface in INTERFACE_PROPERTIES
    }
    for name, mesh, fabricated in (
        ("cluster_thin", args.thin_mesh.expanduser().resolve(), False),
        ("cluster_fabricated", args.fabricated_mesh.expanduser().resolve(), True),
    ):
        config = make_config(
            output,
            name,
            mesh,
            traces,
            conductor_traces,
            fabricated,
            args.order,
            args.coupon_depth,
            args.substrate_permittivity,
            interface_layers,
            int(max(conductors)),
        )
        path = output / f"{name}.json"
        path.write_text(json.dumps(config, indent=2) + "\n")
        heldout_name = f"heldout_{name}"
        heldout = make_config(
            output,
            heldout_name,
            mesh,
            [heldout_trace],
            conductor_traces,
            fabricated,
            args.order,
            args.coupon_depth,
            args.substrate_permittivity,
            interface_layers,
            int(max(conductors)),
        )
        potential = {
            "Index": 1,
            "Attributes": [1],
            "DataFile": str(heldout_trace),
        }
        terminal_attributes = []
        for conductor in range(2, int(max(conductors)) + 1):
            terminal_attributes.append(10 + conductor)
            if fabricated:
                terminal_attributes.append(100 + conductor)
        if terminal_attributes:
            potential["TerminalAttributes"] = terminal_attributes
        heldout["Boundaries"]["PrescribedPotential"] = [potential]
        heldout["Solver"]["Electrostatic"]["ResponseMatrix"] = False
        heldout_path = output / f"{heldout_name}.json"
        heldout_path.write_text(json.dumps(heldout, indent=2) + "\n")
        probe_name = f"probe-{'fabricated' if fabricated else 'thin'}"
        probe = make_config(
            output,
            probe_name,
            mesh,
            probe_traces,
            conductor_traces,
            fabricated,
            args.order,
            args.coupon_depth,
            args.substrate_permittivity,
            interface_layers,
            int(max(conductors)),
        )
        probe_path = output / f"{probe_name}.json"
        probe_path.write_text(json.dumps(probe, indent=2) + "\n")
        print(path)
        print(heldout_path)
        print(probe_path)
    print(
        write_library(
            output,
            args.library_name,
            data["Edges"],
            references,
            open_paths,
            args.radius,
            args.coupon_depth,
            args.metal_thickness,
            args.overetch_depth,
            args.sidewall_angle,
            args.top_rounding,
            args.trench_rounding,
            args.substrate_permittivity,
            interface_layers,
            model_name,
            args.edge_offset_tolerance,
        )
    )


if __name__ == "__main__":
    main()
