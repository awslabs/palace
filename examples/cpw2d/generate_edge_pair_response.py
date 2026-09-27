#!/usr/bin/env python3

"""Generate contour hat bases and Palace configs for a coupled edge-pair coupon.

Curved (axisymmetric) pair coupon: --axisymmetric-radius RHO places the same cross-section
in the (r, z) half-plane with the INNER edge at r = rho (the pair occupies r in
[rho, rho + separation]; mesh_edge_pair_coupon.jl with the same --axisymmetric-radius /
--convexity), so the coupon represents the concentric circular edges of an annular gap /
strip. --convexity names the curvature of the canonical left edge (the model's first edge
e1, conductor 1) relative to its gap direction as the identification records it: convex =
e1's metal inside the bend (e1 is the inner edge of a gap, the outer edge of a strip),
concave = the mirrored placement. The library record then carries Topology
Curved<straight topology>, Kappa = R / rho, Convexity, AxisymmetricRadius and CouponDepth =
2 pi (rho + separation / 2), the full-revolution length of the pair's centreline (the
FEATURES pair measure: the mean of the two sides).
"""

import argparse
import json
import math
from pathlib import Path

import numpy as np


class PairPlacement:
    """Map from the canonical pair cross-section (edges at -+ separation / 2) to the mesh:
    x -> rho + separation / 2 + sign x with sign +1 when e1 is the inner edge (a convex gap
    edge or a concave strip edge), -1 otherwise; the identity for a straight coupon (see
    generate_edge_response.CouponPlacement)."""

    def __init__(self, axisymmetric_radius=0.0, convexity="convex", separation=0.0,
                 strip=False):
        if axisymmetric_radius < 0.0:
            raise ValueError("axisymmetric radius must be nonnegative")
        if convexity not in ("convex", "concave"):
            raise ValueError("convexity must be convex or concave")
        self.rho = axisymmetric_radius
        self.convexity = convexity
        # e1's gap points toward +x (inner -> outward) for a gap, toward -x for a strip.
        self.e1_inner = (convexity == "convex") != bool(strip)
        self.sign = 1.0 if self.e1_inner else -1.0
        self.half_gap = 0.5 * separation

    @property
    def axisymmetric(self):
        return self.rho > 0.0

    @property
    def centreline_radius(self):
        return self.rho + self.half_gap

    def to_mesh(self, points):
        points = np.array(points, dtype=float, copy=True)
        if self.axisymmetric:
            points[..., 0] = self.centreline_radius + self.sign * points[..., 0]
        return points

    def to_canonical(self, points):
        points = np.array(points, dtype=float, copy=True)
        if self.axisymmetric:
            points[..., 0] = self.sign * (points[..., 0] - self.centreline_radius)
        return points


PLACEMENT = PairPlacement()


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
    # The hat traces (Palace DataFile) are in mesh coordinates; the library's basis points
    # stay in the canonical coupon frame (edges at -+ separation / 2), which the matcher
    # maps onto a device pair through the patch frame.
    points = PLACEMENT.to_mesh(
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
    # The right conductor of a different-conductor gap is a terminal at this potential in
    # the held-out solve; the polynomial blends to the potential of the conductor at each
    # cut (zero at the ground) so the trace is a continuous field across the cut (see
    # generate_edge_cluster_response.write_heldout).
    conductor_coefficient = 0.17 if conductor_trace is not None else 0.0

    def free_potential(points):
        points = PLACEMENT.to_canonical(points)
        x = points[:, 0] / radius
        y = points[:, 1] / radius
        polynomial = (
            0.35
            + 0.20 * x
            - 0.15 * y
            + 0.08 * x * y
            + 0.06 * y * y
        )
        if strip:
            return polynomial
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
        if conductor_coefficient:
            targets = np.where(right < left, conductor_coefficient, 0.0)
            return cutoff * polynomial + (1.0 - cutoff) * targets
        return cutoff * polynomial

    # basis_points.csv is canonical, the trace samples are in mesh coordinates.
    coefficients = list(free_potential(PLACEMENT.to_mesh(basis_points)))
    values = free_potential(coordinates)
    if conductor_trace is not None:
        coefficients.append(conductor_coefficient)
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


# See generate_edge_response.EDGE_DISTANCES (the same rule; decision 66 part C).
DEFAULT_EDGE_DISTANCES = [0.2]
EDGE_DISTANCES = list(DEFAULT_EDGE_DISTANCES)


def shells_requested():
    return EDGE_DISTANCES != DEFAULT_EDGE_DISTANCES


def ma_edge_attributes(foot, sidewall):
    """See generate_edge_response.ma_edge_attributes: the fabricated MA's edge points are
    the sidewall endpoints (bottom and top metal edge of every sidewall) under
    --edge-distances shells (per-edge rows: AggregateResponseMatrix false), the foot
    corners otherwise."""
    return sidewall if shells_requested() else foot


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
        "EdgeDistances": sorted(EDGE_DISTANCES),
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
                ma_edge_attributes(edge_attributes, [5, 9] if different_conductors else [5]),
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
    model = {"Mesh": str(mesh), "L0": 1.0e-6, "Lc": coupon_depth}
    if PLACEMENT.axisymmetric:
        model["Axisymmetric"] = True
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
        },
        "Model": model,
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
                "AggregateResponseMatrix": not shells_requested(),
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


def curved_model_fields(topology, topology_name, separation, radius):
    """Library record fields of a curved (axisymmetric) pair coupon: concentric edges with
    the inner one at rho have Kappa = R / rho, the response is a full-revolution energy so
    the CouponDepth is the centreline length 2 pi (rho + separation / 2), Convexity is the
    curvature of the model's first edge e1 relative to its gap direction (convex = metal
    inside the bend: e1 inner for a gap, e1 outer for a strip); a straight coupon is the
    kappa = 0 anchor and keeps its straight record."""
    if not PLACEMENT.axisymmetric:
        return {}
    kappa = radius / PLACEMENT.rho
    return {
        "Name": f"curved-{topology_name}-{separation:g}um-{PLACEMENT.convexity}-"
                f"kappa{kappa:.6g}",
        "Topology": "Curved" + topology,
        "Kappa": kappa,
        "Convexity": PLACEMENT.convexity.capitalize(),
        "AxisymmetricRadius": PLACEMENT.rho,
        "CouponDepth": 2.0 * math.pi * PLACEMENT.centreline_radius,
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
    model.update(curved_model_fields(topology, topology_name, separation, radius))
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
        "TraceLiftVersion": 2,
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
    parser.add_argument("--edge-distances", type=float, nargs="+", default=None,
                        help="localized-energy radii (um) of every interface (default 0.2 alone; the largest must be 0.2)")
    parser.add_argument("--thin-mesh", type=Path)
    parser.add_argument("--fabricated-mesh", type=Path)
    topology = parser.add_mutually_exclusive_group()
    topology.add_argument("--different-conductors", action="store_true")
    topology.add_argument("--strip", action="store_true")
    parser.add_argument("--axisymmetric-radius", type=float, default=0.0,
                        help="curved coupon: inner edge radius rho (um) in the (r, z) half-plane "
                             "(0 = straight); the mesh must be placed the same way")
    parser.add_argument("--convexity", choices=("convex", "concave"), default="convex",
                        help="curved coupon: curvature of the model's first edge e1 relative to "
                             "its gap direction (convex = metal inside the bend: e1 inner for a "
                             "gap, outer for a strip; concave = the mirrored placement)")
    args = parser.parse_args()
    global PLACEMENT
    PLACEMENT = PairPlacement(args.axisymmetric_radius, args.convexity, args.separation,
                              args.strip)
    if PLACEMENT.axisymmetric and PLACEMENT.rho <= args.radius:
        parser.error("--axisymmetric-radius must exceed the coupon radius")
    if args.edge_distances is not None:
        distances = sorted(set(args.edge_distances))
        if not distances or any(d <= 0.0 for d in distances) or max(distances) != DEFAULT_EDGE_DISTANCES[0]:
            parser.error("--edge-distances must be positive and end at the coupon radius 0.2")
        EDGE_DISTANCES[:] = distances
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
        heldout["Solver"]["Electrostatic"]["AggregateResponseMatrix"] = False
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
