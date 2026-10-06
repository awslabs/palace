#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The corner coupon's trace RESOLVABILITY gate (block (b) family 4 round 2, supervisor
decision 328, 2026-10-05).

The coupon's electrostatic solve prescribes every trace hat on the matching box by NODAL
interpolation onto the boundary degrees of freedom of the order-p H1 space
(laplaceoperator.cpp TracePotentialCoefficient / ProjectBdrCoefficient): a hat whose
support holds no boundary node is a zero response row (electrostaticsolver.cpp "Prescribed
trace has no active boundary degrees of freedom"), refused by finalize_corner_response.py's
positive-definiteness check; a hat whose support holds one or two nodes is a response of the
wrong magnitude that nothing refuses (the concave 48.75-degree coupon of decision 325: two
of the five inner free hats of the z = -R cap ring, 7.4 nm apart on a 300-nm far mesh, had
no active node in the fabricated p4 solve; the same hats of the 45 / 52.5 / 56.25 / 60
nodes read 10-70x scattered energies). This gate counts, for every free hat of the trace
basis and every coupon mesh, the boundary nodes of the order-p mesh at which the hat is
nonzero, exactly as the solver interpolates them (the closed Gauss-Lobatto lattice of the
H1 triangle element on every triangle of the matching surface; the box faces are planar, so
the order-2 geometry is affine there), and fails closed when any free hat has fewer than
the required count. It reads the mesh files (Gmsh MSH 2.2, binary or ASCII) and the
generator's trace mesh (trace-vertices.csv / trace-triangles.csv), and records the counts
per hat (trace-resolvability.json).
"""

import argparse
import csv
import json
import struct
import sys
from pathlib import Path

import numpy as np

# The required number of active boundary nodes per free hat (decision 328 (1)), scaled by p^2
# (the node count of a fixed region of an order-p surface mesh). CALIBRATION (the read-only
# audit of the published and built corner caches, family4/round2/audit, 2026-10-05, all at
# p3 unless stated; the five inner free hats of the two cap rings, |x|, |y| <= R / 3 at
# z = +-R, are the sparsest of the concave nodes): the fabricated domain diagonal of those hats
# scatters by a factor max / min 2.6-59 across the five hats at concave 60 degrees (1-3 active
# nodes; the 48.75-degree p4 coupon had 0 on two hats: zero rows), 1.6-1.9 at 75 (17-24 nodes
# at 75-82.5), 1.2-1.6 at 90-105 (28-44 nodes) and settles at <= 1.4 (the genuine angular
# variation) from 120 degrees on (52-56 nodes), as on every convex node whose cap hats hold
# >= 88. The threshold sits between the two regimes: 5 p^2 (45 at p3, 80 at p4). A corner
# coupon mesh sized at the knot gaps (mesh_corner_coupon.jl --trace-mesh, gradation 0.5)
# gives the cap hats 270-290 nodes at p4 (3.4x the gate) and every other free hat more.
MINIMUM_ACTIVE_NODES_PER_ORDER_SQUARED = 5
# A hat value below this at a node is not an active degree of freedom (the solver's zero
# test is the right-hand side norm against 100 eps; the barycentric evaluation here is exact
# to rounding, so this only excludes nodes on the support's boundary).
ACTIVE_VALUE = 1.0e-12
# Nodes closer than this (over R) coincide (the lattice of adjacent triangles shares edge
# nodes).
NODE_COINCIDENCE_OVER_R = 1.0e-9


def required_active_nodes(order):
    return MINIMUM_ACTIVE_NODES_PER_ORDER_SQUARED * int(order) ** 2


def gauss_lobatto_closed_points(order):
    """The closed Gauss-Lobatto points on [0, 1] of the order-p H1 element (MFEM
    poly1d.ClosedPoints with BasisType::GaussLobatto): 0, the roots of P'_p, 1."""
    if order < 1:
        raise ValueError("the element order must be positive")
    if order == 1:
        return np.array([0.0, 1.0])
    legendre = np.polynomial.legendre.Legendre.basis(order).deriv()
    interior = np.sort(np.real(legendre.roots()))
    return np.concatenate(([0.0], 0.5 * (1.0 + interior), [1.0]))


def lagrange_triangle_lattice(order):
    """The reference-triangle nodes (x, y) of MFEM's H1_TriangleElement of that order:
    vertices, cp[i] along each edge, interior (cp[i] / w, cp[j] / w) with w = cp[i] + cp[j] +
    cp[p - i - j] (fe_h1.cpp). The set is symmetric under vertex permutations."""
    cp = gauss_lobatto_closed_points(order)
    p = int(order)
    points = [(cp[0], cp[0]), (cp[p], cp[0]), (cp[0], cp[p])]
    for i in range(1, p):
        points.append((cp[i], cp[0]))
    for i in range(1, p):
        points.append((cp[p - i], cp[i]))
    for i in range(1, p):
        points.append((cp[0], cp[p - i]))
    for j in range(1, p):
        for i in range(1, p - j):
            w = cp[i] + cp[j] + cp[p - i - j]
            points.append((cp[i] / w, cp[j] / w))
    return np.asarray(points)


# Gmsh MSH 2.2 element types: (node count, is a triangle).
MSH_TRIANGLES = {2: 3, 9: 6, 20: 9, 21: 10, 22: 12, 23: 15, 24: 15, 25: 21}
MSH_NODE_COUNTS = {
    1: 2, 2: 3, 3: 4, 4: 4, 5: 8, 6: 6, 7: 5, 8: 3, 9: 6, 10: 9, 11: 10, 12: 27, 13: 18,
    14: 14, 15: 1, 16: 8, 17: 20, 18: 15, 19: 13, 20: 9, 21: 10, 22: 12, 23: 15, 24: 15,
    25: 21, 26: 4, 27: 5, 28: 6, 29: 20, 30: 35, 31: 56,
}


def read_msh2_boundary_triangles(path, physical_tag):
    """The node coordinates and the triangles (node indices, corner nodes first as Gmsh
    orders them) of the MSH 2.2 file's triangles carrying the physical tag."""
    with Path(path).open("rb") as stream:
        data = stream.read()
    header_end = data.index(b"$EndMeshFormat")
    header = data[: header_end].decode().split("\n")
    version, file_type, data_size = header[1].split()
    if not version.startswith("2."):
        raise ValueError(f"{path}: MSH version {version} is not 2.x")
    binary = int(file_type) == 1
    if binary:
        if int(data_size) != 8:
            raise ValueError(f"{path}: unsupported binary data size {data_size}")
        one = data.index(b"\n", data.index(version.encode())) + 1
        (endian_probe,) = struct.unpack("<i", data[one : one + 4])
        endian = "<" if endian_probe == 1 else ">"
    nodes_start = data.index(b"$Nodes") + len(b"$Nodes\n")
    count_end = data.index(b"\n", nodes_start)
    node_count = int(data[nodes_start:count_end])
    coordinates = {}
    cursor = count_end + 1
    if binary:
        record = struct.Struct(endian + "iddd")
        for _ in range(node_count):
            tag, x, y, z = record.unpack_from(data, cursor)
            coordinates[tag] = (x, y, z)
            cursor += record.size
    else:
        end = data.index(b"$EndNodes", cursor)
        for line in data[cursor:end].decode().split("\n"):
            if line.strip():
                tag, x, y, z = line.split()
                coordinates[int(tag)] = (float(x), float(y), float(z))
    elements_start = data.index(b"$Elements") + len(b"$Elements\n")
    count_end = data.index(b"\n", elements_start)
    element_count = int(data[elements_start:count_end])
    cursor = count_end + 1
    triangles = []
    if binary:
        block_header = struct.Struct(endian + "iii")
        read = 0
        while read < element_count:
            element_type, block_count, tag_count = block_header.unpack_from(data, cursor)
            cursor += block_header.size
            node_count = MSH_NODE_COUNTS[element_type]
            record = struct.Struct(endian + "i" * (1 + tag_count + node_count))
            for _ in range(block_count):
                values = record.unpack_from(data, cursor)
                cursor += record.size
                if element_type in MSH_TRIANGLES and values[1] == physical_tag:
                    triangles.append(values[1 + tag_count :])
            read += block_count
    else:
        end = data.index(b"$EndElements", cursor)
        for line in data[cursor:end].decode().split("\n"):
            if not line.strip():
                continue
            values = [int(v) for v in line.split()]
            element_type, tag_count = values[1], values[2]
            if element_type in MSH_TRIANGLES and values[3] == physical_tag:
                triangles.append(tuple(values[3 + tag_count :]))
    if not triangles:
        raise ValueError(f"{path}: no triangles carry the physical tag {physical_tag}")
    return coordinates, triangles


def boundary_lattice_nodes(coordinates, triangles, order, radius, planarity_tolerance=1.0e-9):
    """The order-p boundary node coordinates of the triangles (each triangle's lattice mapped
    affinely by its three corner nodes; the triangles must be straight-sided: every
    higher-order node of the file is checked against the affine map), de-duplicated."""
    lattice = lagrange_triangle_lattice(order)
    points = []
    for triangle in triangles:
        a, b, c = (np.asarray(coordinates[triangle[k]]) for k in range(3))
        for extra in triangle[3:]:
            # Gmsh's higher-order triangle nodes: edge nodes in order, then face nodes; the
            # affine map holds on a planar face with straight edges, which the box is.
            point = np.asarray(coordinates[extra])
            ab, ac = b - a, c - a
            normal = np.cross(ab, ac)
            if abs(np.dot(point - a, normal)) > planarity_tolerance * radius * np.linalg.norm(normal):
                raise ValueError("a matching-surface triangle is not planar")
        points.append(a + np.outer(lattice[:, 0], b - a) + np.outer(lattice[:, 1], c - a))
    points = np.concatenate(points)
    keys = np.round(points / (NODE_COINCIDENCE_OVER_R * radius)).astype(np.int64)
    _, unique = np.unique(keys, axis=0, return_index=True)
    return points[np.sort(unique)]


def read_trace_mesh(directory):
    """The generator's trace mesh: vertex points, per-vertex basis weights (a knot: {knot:
    1}; a slave: {parent_a: weight_a, parent_b: 1 - weight_a}) and triangles (0-based)."""
    directory = Path(directory)
    vertices = []
    weights = []
    with (directory / "trace-vertices.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            vertices.append((float(row["x"]), float(row["y"]), float(row["z"])))
            basis = int(row["basis"])
            if basis > 0:
                weights.append({basis - 1: 1.0})
            else:
                weight_a = float(row["weight_a"])
                weights.append(
                    {int(row["parent_a"]) - 1: weight_a, int(row["parent_b"]) - 1: 1.0 - weight_a}
                )
    triangles = []
    with (directory / "trace-triangles.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            triangles.append(
                (int(row["vertex_i"]) - 1, int(row["vertex_j"]) - 1, int(row["vertex_k"]) - 1)
            )
    return np.asarray(vertices), weights, triangles


def hat_active_counts(nodes, vertices, weights, triangles, basis_size, radius):
    """For every basis hat the number of nodes at which it exceeds ACTIVE_VALUE (nodes
    located in the trace triangles by barycentric coordinates; a node on a shared edge is
    counted once, with the continuous value) and the maximum value sampled."""
    counts = np.zeros(basis_size, dtype=int)
    maxima = np.zeros(basis_size)
    inside_tolerance = 1.0e-9
    assigned = np.zeros(len(nodes), dtype=bool)
    for triangle in triangles:
        a, b, c = (vertices[v] for v in triangle)
        ab, ac = b - a, c - a
        normal = np.cross(ab, ac)
        area2 = np.linalg.norm(normal)
        if area2 <= 0.0:
            continue
        normal = normal / area2
        candidates = np.flatnonzero(~assigned)
        if candidates.size == 0:
            break
        relative = nodes[candidates] - a
        off_plane = np.abs(relative @ normal) <= inside_tolerance * radius
        candidates = candidates[off_plane]
        relative = relative[off_plane]
        # Barycentric coordinates in the triangle's plane.
        d00, d01, d11 = ab @ ab, ab @ ac, ac @ ac
        d20, d21 = relative @ ab, relative @ ac
        denominator = d00 * d11 - d01 * d01
        v = (d11 * d20 - d01 * d21) / denominator
        w = (d00 * d21 - d01 * d20) / denominator
        u = 1.0 - v - w
        inside = (u >= -inside_tolerance) & (v >= -inside_tolerance) & (w >= -inside_tolerance)
        candidates, u, v, w = candidates[inside], u[inside], v[inside], w[inside]
        assigned[candidates] = True
        for vertex, barycentric in zip(triangle, (u, v, w)):
            for basis, weight in weights[vertex].items():
                values = weight * barycentric
                active = values > ACTIVE_VALUE
                counts[basis] += int(np.count_nonzero(active))
                if values.size:
                    maxima[basis] = max(maxima[basis], float(values.max()))
    return counts, maxima, int(np.count_nonzero(assigned)), int(len(nodes))


def audit_mesh(mesh, trace_directory, order, radius, zero_trace_indices, physical_tag=1):
    coordinates, triangles = read_msh2_boundary_triangles(mesh, physical_tag)
    nodes = boundary_lattice_nodes(coordinates, triangles, order, radius)
    vertices, weights, trace_triangles = read_trace_mesh(trace_directory)
    basis_size = max(max(w) for w in weights) + 1
    counts, maxima, located, total = hat_active_counts(
        nodes, vertices, weights, trace_triangles, basis_size, radius
    )
    zero = np.zeros(basis_size, dtype=bool)
    zero[[index - 1 for index in zero_trace_indices]] = True
    required = required_active_nodes(order)
    failing = [int(k) + 1 for k in np.flatnonzero((~zero) & (counts < required))]
    return {
        "Mesh": str(mesh),
        "Order": int(order),
        "MatchingSurfaceTriangles": len(triangles),
        "BoundaryNodes": total,
        "BoundaryNodesOnTraceSurface": located,
        "RequiredActiveNodes": required,
        "ActiveNodes": counts.tolist(),
        "MaximumHatValue": maxima.tolist(),
        "MinimumActiveNodesFreeHat": int(counts[~zero].min()),
        "FailingFreeHats": failing,
        "Passed": not failing,
    }


def read_zero_trace_indices(trace_directory):
    library = json.loads((Path(trace_directory) / "process-library.json").read_text())
    (model,) = library["Models"]
    return [int(index) for index in model["ZeroTraceIndices"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("directory", type=Path, help="the coupon cache (trace mesh + library)")
    parser.add_argument("--mesh", action="append", required=True, type=Path,
                        help="a coupon mesh to audit (repeatable: fabricated, thin, h-refined)")
    parser.add_argument("--order", type=int, required=True, help="the solve order p")
    parser.add_argument("--radius", type=float, required=True, help="the matching radius R")
    parser.add_argument("--report", type=Path, default=None,
                        help="the JSON report (default DIRECTORY/trace-resolvability.json)")
    parser.add_argument("--audit-only", action="store_true",
                        help="report without failing (the read-only audit of published caches)")
    args = parser.parse_args()
    zero = read_zero_trace_indices(args.directory)
    audits = [
        audit_mesh(mesh, args.directory, args.order, args.radius, zero) for mesh in args.mesh
    ]
    passed = all(audit["Passed"] for audit in audits)
    report = {
        "Version": 1,
        "Gate": "TraceResolvability",
        "MinimumActiveNodesPerOrderSquared": MINIMUM_ACTIVE_NODES_PER_ORDER_SQUARED,
        "ZeroTraceIndices": zero,
        "Meshes": audits,
        "Passed": passed,
    }
    report_path = args.report or (args.directory / "trace-resolvability.json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    for audit in audits:
        print(
            f"{Path(audit['Mesh']).name}: order {audit['Order']}, "
            f"{audit['MatchingSurfaceTriangles']} matching-surface triangles, "
            f"{audit['BoundaryNodes']} boundary nodes ({audit['BoundaryNodesOnTraceSurface']} "
            f"on the trace surface); free hats hold >= {audit['MinimumActiveNodesFreeHat']} "
            f"active nodes (required {audit['RequiredActiveNodes']}); "
            + ("PASS" if audit["Passed"] else f"FAIL on hats {audit['FailingFreeHats']}")
        )
    if not passed and not args.audit_only:
        print("trace resolvability gate FAILED", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
