#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Generate a periodic, extruded CPW mesh with zero-thickness metal sheets."""

import argparse
import itertools
import math
from pathlib import Path


def uniform_interval(start, end, count):
    return [start + (end - start) * i / count for i in range(count + 1)]


def cosine_interval(start, end, count):
    return [
        start + 0.5 * (end - start) * (1.0 - math.cos(math.pi * i / count))
        for i in range(count + 1)
    ]


def geometric_interval(start, end, count, first_step):
    length = end - start
    if count == 1:
        return [start, end]
    if not (0.0 < first_step < length / count):
        raise ValueError("The first geometric interval must be smaller than the mean")

    lower, upper = 1.0, 4.0
    for _ in range(100):
        ratio = 0.5 * (lower + upper)
        total = first_step * (ratio**count - 1.0) / (ratio - 1.0)
        if total < length:
            lower = ratio
        else:
            upper = ratio
    ratio = 0.5 * (lower + upper)
    weights = [first_step * ratio**i for i in range(count)]
    scale = length / sum(weights)
    points = [start]
    for weight in weights:
        points.append(points[-1] + scale * weight)
    points[-1] = end
    return points


def merge_intervals(*intervals):
    points = []
    for interval in intervals:
        points.extend(interval if not points else interval[1:])
    return points


def mirror_positive(points):
    if points[0] != 0.0:
        raise ValueError("Positive coordinate list must begin at zero")
    return [-value for value in reversed(points[1:])] + points


def subdivide_intervals(points, factor):
    if not isinstance(factor, int) or factor < 1:
        raise ValueError("The interval subdivision factor must be a positive integer")
    return [
        points[i] + (points[i + 1] - points[i]) * j / factor
        for i in range(len(points) - 1)
        for j in range(factor)
    ] + [points[-1]]


def determinant(vertices, coordinates):
    x0 = coordinates[vertices[0]]
    columns = [
        tuple(coordinates[vertices[column]][d] - x0[d] for d in range(3))
        for column in range(1, 4)
    ]
    return (
        columns[0][0]
        * (columns[1][1] * columns[2][2] - columns[1][2] * columns[2][1])
        - columns[1][0]
        * (columns[0][1] * columns[2][2] - columns[0][2] * columns[2][1])
        + columns[2][0]
        * (columns[0][1] * columns[1][2] - columns[0][2] * columns[1][1])
    )


def generate_mesh(
    output,
    strip_half_width=5.0,
    gap=3.0,
    ground_outer=80.0,
    box_half_width=160.0,
    box_half_height=160.0,
    length=10.0,
    strip_half_intervals=2,
    gap_intervals=2,
    ground_intervals=8,
    outer_intervals=4,
    height_intervals=10,
    longitudinal_intervals=3,
    transverse_subdivision_factor=1,
):
    interval_counts = {
        "strip_half_intervals": strip_half_intervals,
        "gap_intervals": gap_intervals,
        "ground_intervals": ground_intervals,
        "outer_intervals": outer_intervals,
        "height_intervals": height_intervals,
        "longitudinal_intervals": longitudinal_intervals,
    }
    invalid_counts = [
        name for name, count in interval_counts.items() if not isinstance(count, int) or count < 1
    ]
    if invalid_counts:
        raise ValueError(
            "Every interval count must be a positive integer: "
            + ", ".join(invalid_counts)
        )
    if (
        not isinstance(transverse_subdivision_factor, int)
        or transverse_subdivision_factor < 1
    ):
        raise ValueError("The transverse subdivision factor must be a positive integer")

    ground_inner = strip_half_width + gap
    if not (
        0.0 < strip_half_width < ground_inner < ground_outer < box_half_width
        and box_half_height > 0.0
        and length > 0.0
    ):
        raise ValueError("Invalid CPW dimensions")

    x_positive = merge_intervals(
        uniform_interval(0.0, strip_half_width, strip_half_intervals),
        uniform_interval(strip_half_width, ground_inner, gap_intervals),
        cosine_interval(ground_inner, ground_outer, ground_intervals),
        geometric_interval(
            ground_outer,
            box_half_width,
            outer_intervals,
            min(gap / 2.0, (box_half_width - ground_outer) / (2 * outer_intervals)),
        ),
    )
    y_positive = geometric_interval(
        0.0,
        box_half_height,
        height_intervals,
        min(gap / 2.0, box_half_height / (2 * height_intervals)),
    )
    x_values = mirror_positive(x_positive)
    y_values = mirror_positive(y_positive)
    x_values = subdivide_intervals(x_values, transverse_subdivision_factor)
    y_values = subdivide_intervals(y_values, transverse_subdivision_factor)
    z_values = uniform_interval(0.0, length, longitudinal_intervals)

    nx, ny, nz = len(x_values), len(y_values), len(z_values)

    def vertex(i, j, k):
        return (k * ny + j) * nx + i

    coordinates = [
        (x_values[i], y_values[j], z_values[k])
        for k in range(nz)
        for j in range(ny)
        for i in range(nx)
    ]

    elements = []
    permutations = tuple(itertools.permutations(range(3)))
    for k in range(nz - 1):
        for j in range(ny - 1):
            for i in range(nx - 1):
                cell_origin = (i, j, k)
                material = 1 if 0.5 * (y_values[j] + y_values[j + 1]) < 0.0 else 2
                for permutation in permutations:
                    index = list(cell_origin)
                    tetrahedron = [vertex(*index)]
                    for axis in permutation:
                        index[axis] += 1
                        tetrahedron.append(vertex(*index))
                    if determinant(tetrahedron, coordinates) < 0.0:
                        tetrahedron[1], tetrahedron[2] = (
                            tetrahedron[2],
                            tetrahedron[1],
                        )
                    elements.append((material, tuple(tetrahedron)))

    faces = {}
    for element, (_, tetrahedron) in enumerate(elements):
        for omitted in range(4):
            face = tuple(
                sorted(tetrahedron[local] for local in range(4) if local != omitted)
            )
            faces.setdefault(face, []).append(element)

    tolerance = 1.0e-12 * max(box_half_width, box_half_height, length)
    boundaries = []
    conductor_counts = {10: 0, 11: 0}
    for face, adjacent in faces.items():
        xyz = [coordinates[node] for node in face]

        def on_plane(axis, value):
            return all(abs(point[axis] - value) <= tolerance for point in xyz)

        attribute = None
        if len(adjacent) == 1:
            if on_plane(2, 0.0):
                attribute = 20
            elif on_plane(2, length):
                attribute = 21
            elif (
                on_plane(0, -box_half_width)
                or on_plane(0, box_half_width)
                or on_plane(1, -box_half_height)
                or on_plane(1, box_half_height)
            ):
                attribute = 12
            else:
                raise RuntimeError(f"Unclassified external face: {face}")
        elif len(adjacent) == 2 and on_plane(1, 0.0):
            xmin = min(point[0] for point in xyz)
            xmax = max(point[0] for point in xyz)
            if xmin >= -strip_half_width - tolerance and xmax <= strip_half_width + tolerance:
                attribute = 10
            elif (
                xmin >= ground_inner - tolerance and xmax <= ground_outer + tolerance
            ) or (
                xmin >= -ground_outer - tolerance
                and xmax <= -ground_inner + tolerance
            ):
                attribute = 11

        if attribute is not None:
            boundaries.append((attribute, face))
            if attribute in conductor_counts:
                conductor_counts[attribute] += 1

    expected_trace_faces = (
        4
        * strip_half_intervals
        * transverse_subdivision_factor
        * longitudinal_intervals
    )
    expected_ground_faces = (
        4 * ground_intervals * transverse_subdivision_factor * longitudinal_intervals
    )
    if conductor_counts != {10: expected_trace_faces, 11: expected_ground_faces}:
        raise RuntimeError(
            f"Unexpected conductor face counts: {conductor_counts}, "
            f"expected trace={expected_trace_faces}, ground={expected_ground_faces}"
        )

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="ascii") as mesh:
        mesh.write("MFEM mesh v1.0\n\n")
        mesh.write("dimension\n3\n\n")
        mesh.write(f"elements\n{len(elements)}\n")
        for attribute, tetrahedron in elements:
            mesh.write(f"{attribute} 4 {' '.join(map(str, tetrahedron))}\n")
        mesh.write(f"\nboundary\n{len(boundaries)}\n")
        for attribute, face in boundaries:
            mesh.write(f"{attribute} 2 {' '.join(map(str, face))}\n")
        mesh.write(f"\nvertices\n{len(coordinates)}\n3\n")
        for point in coordinates:
            mesh.write(" ".join(f"{value:.17g}" for value in point) + "\n")

    print(f"Wrote {output}")
    print(
        f"  vertices={len(coordinates)}, tetrahedra={len(elements)}, "
        f"boundary triangles={len(boundaries)}"
    )
    print(
        f"  trace triangles={conductor_counts[10]}, "
        f"ground triangles={conductor_counts[11]}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("cpw3d_singular.mesh"),
    )
    parser.add_argument("--strip-half-width", type=float, default=5.0)
    parser.add_argument("--gap", type=float, default=3.0)
    parser.add_argument("--ground-outer", type=float, default=80.0)
    parser.add_argument("--box-half-width", type=float, default=160.0)
    parser.add_argument("--box-half-height", type=float, default=160.0)
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--strip-half-intervals", type=int, default=2)
    parser.add_argument("--gap-intervals", type=int, default=2)
    parser.add_argument("--longitudinal-intervals", type=int, default=3)
    parser.add_argument("--ground-intervals", type=int, default=8)
    parser.add_argument("--outer-intervals", type=int, default=4)
    parser.add_argument("--height-intervals", type=int, default=10)
    parser.add_argument("--transverse-subdivision-factor", type=int, default=1)
    args = parser.parse_args()
    generate_mesh(
        args.output,
        strip_half_width=args.strip_half_width,
        gap=args.gap,
        ground_outer=args.ground_outer,
        box_half_width=args.box_half_width,
        box_half_height=args.box_half_height,
        length=args.length,
        strip_half_intervals=args.strip_half_intervals,
        gap_intervals=args.gap_intervals,
        outer_intervals=args.outer_intervals,
        longitudinal_intervals=args.longitudinal_intervals,
        ground_intervals=args.ground_intervals,
        height_intervals=args.height_intervals,
        transverse_subdivision_factor=args.transverse_subdivision_factor,
    )


if __name__ == "__main__":
    main()
