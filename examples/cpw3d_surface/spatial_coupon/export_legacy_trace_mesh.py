#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Export explicit mortar trace connectivity from retained spatial-coupon traces."""

import argparse
import csv
from pathlib import Path


def coordinate(row):
    return tuple(float(row[key]) for key in ("x", "y", "z"))


def solve(matrix, rhs):
    augmented = [list(row) + [value] for row, value in zip(matrix, rhs)]
    size = len(rhs)
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1.0e-24:
            raise ValueError("legacy trace coordinate transform is singular")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                value - scale * source
                for value, source in zip(augmented[row], augmented[column])
            ]
    return [augmented[row][-1] for row in range(size)]


def affine_transform(source, target):
    normal = [[0.0] * 4 for _ in range(4)]
    rhs = [[0.0] * 4 for _ in range(3)]
    for raw, canonical in zip(source, target):
        row = (*raw, 1.0)
        for i in range(4):
            for j in range(4):
                normal[i][j] += row[i] * row[j]
            for d in range(3):
                rhs[d][i] += row[i] * canonical[d]
    coefficients = [solve(normal, values) for values in rhs]

    def transform(point):
        row = (*point, 1.0)
        return tuple(sum(coefficients[d][i] * row[i] for i in range(4)) for d in range(3))

    maximum_error = max(
        max(abs(a - b) for a, b in zip(transform(raw), canonical))
        for raw, canonical in zip(source, target)
    )
    if maximum_error > 1.0e-8:
        raise ValueError(f"legacy trace affine-fit error is {maximum_error:.3e}")
    return transform


def trace_values(path):
    values = {}
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            point = coordinate(row)
            value = float(row["V"])
            previous = values.setdefault(point, value)
            if abs(previous - value) > 1.0e-12:
                raise ValueError(f"inconsistent trace value at {point} in {path}")
    return values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--generation", type=Path, required=True)
    parser.add_argument("--basis-points", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--zero-trace-index", type=int, action="append", default=[])
    parser.add_argument("--conductor-count", type=int, default=1)
    args = parser.parse_args()

    root = args.generation.expanduser().resolve()
    basis_path = (args.basis_points or root / "basis-points.csv").expanduser().resolve()
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)

    triangles = {}
    vertex_order = []
    vertex_lookup = {}
    with (root / "zero-trace.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            point = coordinate(row)
            vertex = vertex_lookup.get(point)
            if vertex is None:
                vertex = len(vertex_order)
                vertex_lookup[point] = vertex
                vertex_order.append(point)
            triangles.setdefault(int(row["triangle"]), []).append(vertex)
    if any(len(vertices) != 3 for vertices in triangles.values()):
        raise ValueError("legacy zero trace does not contain triangular connectivity")

    canonical_basis = []
    with basis_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            canonical_basis.append(coordinate(row))
    raw_basis = []
    basis_by_vertex = {}
    for basis in range(1, len(canonical_basis) + 1):
        path = root / "traces" / f"basis-{basis:04d}.csv"
        if not path.is_file():
            # Older corner generation uses three-digit names.
            path = root / "traces" / f"basis-{basis:03d}.csv"
        values = trace_values(path)
        active = [point for point, value in values.items() if abs(value - 1.0) <= 1.0e-12]
        if len(active) != 1 or active[0] not in vertex_lookup:
            raise ValueError(f"cannot identify legacy trace vertex for basis {basis}")
        raw_basis.append(active[0])
        basis_by_vertex[vertex_lookup[active[0]]] = basis

    transform = affine_transform(raw_basis, canonical_basis)
    conductor_values = {}
    for conductor in range(2, args.conductor_count + 1):
        path = root / f"probe-conductor-{conductor}.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        conductor_values[conductor] = trace_values(path)
    constrained_basis = set(args.zero_trace_index)

    with (output / "trace-vertices.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("vertex", "x", "y", "z", "basis", "conductor"))
        for index, raw in enumerate(vertex_order, start=1):
            basis = basis_by_vertex.get(index - 1, 0)
            conductor = 1 if basis in constrained_basis else 0
            if basis == 0:
                conductor = 1
                for candidate, values in conductor_values.items():
                    if abs(values.get(raw, 0.0) - 1.0) <= 1.0e-12:
                        conductor = candidate
                        break
            writer.writerow((index, *[f"{value:.16e}" for value in transform(raw)], basis, conductor))

    with (output / "trace-triangles.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("triangle", "vertex_i", "vertex_j", "vertex_k"))
        for index in sorted(triangles):
            writer.writerow((index, *[vertex + 1 for vertex in triangles[index]]))

    print(output / "trace-vertices.csv")
    print(output / "trace-triangles.csv")


if __name__ == "__main__":
    main()
