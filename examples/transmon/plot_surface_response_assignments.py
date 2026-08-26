#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Plot chip-plane metal geometry colored by assigned surface-response model family."""

import argparse
import csv
import json
import subprocess
import tempfile
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Circle, Patch


def ascii_mesh(path, gmsh, temporary):
    with path.open("rb") as stream:
        header = stream.read(32)
    if b"2.2 0 8" in header:
        return path
    output = Path(temporary) / "mesh-ascii.msh"
    subprocess.run(
        [
            str(gmsh),
            str(path),
            "-save",
            "-format",
            "msh2",
            "-setnumber",
            "Mesh.Binary",
            "0",
            "-o",
            str(output),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    return output


def chip_plane_triangles(path, attributes, chip_z, tolerance):
    lines = path.read_text().splitlines()
    node_section = lines.index("$Nodes")
    node_count = int(lines[node_section + 1])
    nodes = {}
    for line in lines[node_section + 2 : node_section + 2 + node_count]:
        values = line.split()
        nodes[int(values[0])] = tuple(float(value) for value in values[1:4])

    element_section = lines.index("$Elements")
    element_count = int(lines[element_section + 1])
    triangles = []
    for line in lines[element_section + 2 : element_section + 2 + element_count]:
        values = [int(value) for value in line.split()]
        element_type = values[1]
        tag_count = values[2]
        tags = values[3 : 3 + tag_count]
        element_nodes = values[3 + tag_count :]
        if element_type not in (2, 9) or not tags or tags[0] not in attributes:
            continue
        corners = element_nodes[:3]
        points = [nodes[node] for node in corners]
        if all(abs(point[2] - chip_z) <= tolerance for point in points):
            triangles.append([(point[0], point[1]) for point in points])
    return triangles


def perimeter_segments(triangles):
    counts = Counter()
    points = {}
    for triangle in triangles:
        for first, second in zip(triangle, triangle[1:] + triangle[:1]):
            a = tuple(round(value, 9) for value in first)
            b = tuple(round(value, 9) for value in second)
            edge = tuple(sorted((a, b)))
            counts[edge] += 1
            points[a] = first
            points[b] = second
    return [
        [points[first], points[second]]
        for (first, second), count in counts.items()
        if count == 1
    ]


def read_patches(path):
    with path.open(newline="") as stream:
        return [
            {key.strip(): float(value) for key, value in row.items()}
            for row in csv.DictReader(stream, skipinitialspace=True)
        ]


def family(topology):
    value = topology.lower()
    if "spatial" in value:
        return "spatial"
    if "corner" in value or "endpoint" in value or "junction" in value:
        return "corner"
    return "straight"


def support_footprint(patch, support_points):
    """Transform local SupportPoints corners into a chip-plane polygon in microns."""
    origin = [patch[f"origin {axis} (m)"] * 1.0e6 for axis in "xyz"]
    axes = [
        [patch[f"axis {name} {axis}"] for axis in "xyz"] for name in ("u", "v", "w")
    ]
    corners = []
    for local in support_points:
        point = list(origin)
        for axis in range(3):
            for direction in range(3):
                point[direction] += local[axis] * axes[axis][direction]
        corners.append(point)
    xs = [corner[0] for corner in corners]
    ys = [corner[1] for corner in corners]
    return [
        (min(xs), min(ys)),
        (max(xs), min(ys)),
        (max(xs), max(ys)),
        (min(xs), max(ys)),
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--patches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gmsh", type=Path, default=Path("gmsh"))
    parser.add_argument("--metal-attributes", type=int, nargs="+", default=[5, 6, 7, 9])
    parser.add_argument(
        "--chip-z", type=float, default=0.0, help="Chip plane in mesh units"
    )
    parser.add_argument("--matching-radius", type=float, default=2.0, help="Microns")
    parser.add_argument(
        "--library",
        type=Path,
        help="Process library used by the run; enables exact SupportPoints footprints",
    )
    parser.add_argument("--xlim", type=float, nargs=2)
    parser.add_argument("--ylim", type=float, nargs=2)
    parser.add_argument(
        "--label-models",
        action="store_true",
        help="Annotate corner and spatial patches with runtime model indices",
    )
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as temporary:
        mesh = ascii_mesh(args.mesh.expanduser().resolve(), args.gmsh, temporary)
        triangles = chip_plane_triangles(
            mesh, set(args.metal_attributes), args.chip_z, 1.0e-8
        )
    if not triangles:
        raise ValueError("No selected chip-plane metal triangles were found")

    metadata = json.loads(args.metadata.read_text())
    catalog = {
        int(model["Index"]): model
        for model in metadata["SurfaceResponse"]["ModelCatalog"]
    }
    patches = read_patches(args.patches)
    support_by_name = {}
    if args.library:
        library = json.loads(args.library.expanduser().read_text())
        for model in library.get("Models", []):
            if model.get("SupportPoints"):
                support_by_name[model["Name"]] = model["SupportPoints"]
    colors = {"straight": "#2878b5", "corner": "#f28e2b", "spatial": "#d62728"}

    figure, axis = plt.subplots(figsize=(12, 10))
    axis.add_collection(
        PolyCollection(triangles, facecolor="#e6e6e6", edgecolor="none", rasterized=True)
    )
    axis.add_collection(
        LineCollection(perimeter_segments(triangles), colors="#333333", linewidths=0.35)
    )

    grouped = {name: [] for name in colors}
    labels = []
    footprints = []
    for patch in patches:
        model_index = int(patch["model"])
        model = catalog[model_index]
        assigned_family = family(model["Topology"])
        point = (1.0e6 * patch["origin x (m)"], 1.0e6 * patch["origin y (m)"])
        grouped[assigned_family].append(point)
        support = support_by_name.get(model["Name"])
        if assigned_family == "spatial" and support:
            footprints.append(support_footprint(patch, support))
        if args.label_models and assigned_family != "straight":
            labels.append((*point, model_index))
    if grouped["straight"]:
        x, y = zip(*grouped["straight"])
        axis.scatter(x, y, s=1.0, color=colors["straight"], alpha=0.65, linewidths=0)
    for name in ("corner", "spatial"):
        for x, y in grouped[name]:
            if name == "spatial" and footprints:
                continue
            axis.add_patch(
                Circle(
                    (x, y),
                    args.matching_radius,
                    facecolor=colors[name],
                    edgecolor="none",
                    alpha=0.28,
                )
            )
        if grouped[name]:
            x, y = zip(*grouped[name])
            axis.scatter(x, y, s=7, color=colors[name], linewidths=0)
    if footprints:
        axis.add_collection(
            PolyCollection(
                footprints,
                facecolor=colors["spatial"],
                edgecolor=colors["spatial"],
                alpha=0.20,
                linewidths=1.0,
            )
        )
    for x, y, model in labels:
        axis.annotate(str(model), (x, y), xytext=(3, 3), textcoords="offset points", fontsize=6)

    axis.autoscale()
    if args.xlim:
        axis.set_xlim(args.xlim)
    if args.ylim:
        axis.set_ylim(args.ylim)
    axis.set_aspect("equal")
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_title("Single-transmon surface-response assignment at the chip plane")
    axis.legend(
        handles=[
            Patch(facecolor=colors["straight"], label="Straight-edge quadrature"),
            Patch(facecolor=colors["corner"], label="Corner neighborhood"),
            Patch(
                facecolor=colors["spatial"],
                label=(
                    "Spatial edge-cluster support"
                    if footprints
                    else "Spatial edge-cluster neighborhood"
                ),
            ),
        ],
        loc="best",
    )
    axis.text(
        0.01,
        0.01,
        "  ".join(f"{name}: {len(values)} patches" for name, values in grouped.items()),
        transform=axis.transAxes,
        fontsize=8,
        va="bottom",
    )
    figure.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=250)
    print(args.output)


if __name__ == "__main__":
    main()
