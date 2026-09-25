#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Streaming census of a large MSH 2.2 mesh (binary or ASCII), and clipped metal-triangle extraction.

Memory-maps the file (never `read()`s it) and streams the element section in chunks so a
chip-scale mesh (10^7 tetrahedra, 1 GB) is inventoried in a few hundred MB: per (dimension, physical attribute) the element count, the
physical name, the bounding box and, for surface attributes, the area-weighted unit normal
(|n_z| = 1 for a metal plane, 0 for walls) — what the preflight configuration needs to place
each chip's metal with its own EdgeFrameNormal.

    python3 -m surface_response_identification.mesh_census MESH --output census.json
    python3 -m surface_response_identification.mesh_census MESH --output tris.npz \\
        --extract 10 17 25 --window X0 X1 Y0 Y1 [--window ...]     # xy corners (float32), mean z and attribute of the
                                                                    # triangles of these attributes
                                                                    # inside any window (plot fill)
"""
import argparse
import hashlib
import json
import mmap
import struct

import numpy as np

from .msh2 import CORNER_NODES, ELEMENT_DIMENSION, NODES_PER_TYPE, TRIANGLE_TYPES, ascii_element_lines

CHUNK = 1 << 20


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 24), b""):
            digest.update(block)
    return digest.hexdigest()


def _section(data, name):
    start = data.find(f"${name}\n".encode()) + len(name) + 2
    end = data.find(f"$End{name}".encode(), start)
    return start, end


def open_mesh(path):
    """(mmap, node coordinates indexed by tag (dense: tag -> row), physical names, element source).

    Binary and ASCII MSH 2.2 are both streamed: the element source is a list of binary blocks
    (type, count, tag count, offset) or the (start, end) byte range of the ASCII element lines."""
    source = open(path, "rb")
    data = mmap.mmap(source.fileno(), 0, access=mmap.ACCESS_READ)
    fmt_start, fmt_end = _section(data, "MeshFormat")
    version, file_type, data_size = data[fmt_start:fmt_end].split()[:3]
    if version != b"2.2" or file_type not in (b"0", b"1") or (file_type == b"1" and data_size != b"8"):
        raise ValueError(f"{path}: MSH 2.2 required (found {version} {file_type} {data_size})")
    binary = file_type == b"1"
    names = {}
    if data.find(b"$PhysicalNames\n") >= 0:
        start, end = _section(data, "PhysicalNames")
        lines = data[start:end].decode().splitlines()
        for line in lines[1 : 1 + int(lines[0])]:
            dim, tag, name = line.split(maxsplit=2)
            names[(int(dim), int(tag))] = name.strip().strip('"')
    start, end = _section(data, "Nodes")
    count_end = data.find(b"\n", start)
    node_count = int(data[start:count_end])
    if binary:
        record = np.dtype([("tag", "<i4"), ("xyz", "<f8", (3,))])
        nodes = np.frombuffer(data, dtype=record, count=node_count, offset=count_end + 1)
        tags, xyz = nodes["tag"], nodes["xyz"]
    else:
        rows = np.fromstring(data[count_end + 1 : end], dtype=np.float64, sep=" ").reshape(node_count, 4)
        tags, xyz = rows[:, 0].astype(np.int64), rows[:, 1:4]
    coordinates = np.full((int(tags.max()) + 1, 3), np.nan)
    coordinates[tags] = xyz
    start, end = _section(data, "Elements")
    count_end = data.find(b"\n", start)
    element_count = int(data[start:count_end])
    if not binary:
        return data, coordinates, names, ("ascii", count_end + 1, end, element_count)
    blocks = []
    position = count_end + 1
    seen = 0
    while seen < element_count:
        element_type, following, tag_count = struct.unpack_from("<3i", data, position)
        position += 12
        width = 1 + tag_count + NODES_PER_TYPE[element_type]
        blocks.append((element_type, following, tag_count, position))
        position += following * width * 4
        seen += following
    return data, coordinates, names, ("binary", blocks)


def iter_elements(data, elements, dimensions=(0, 1, 2, 3)):
    """Yield (element_type, physical tags (m,), corner node tags (m, k)) in chunks."""
    if elements[0] == "ascii":
        _, start, end, _ = elements
        for width, rows in ascii_element_lines(data, start, end):
            for element_type in np.unique(rows[:, 1]):
                element_type = int(element_type)
                if ELEMENT_DIMENSION[element_type] not in dimensions:
                    continue
                select = rows[:, 1] == element_type
                block = rows[select]
                tag_count = int(block[0, 2])
                if not (block[:, 2] == tag_count).all() or width != 3 + tag_count + NODES_PER_TYPE[element_type]:
                    raise ValueError(f"inconsistent ASCII element lines of type {element_type} and width {width}")
                physical = block[:, 3] if tag_count else np.zeros(len(block), dtype=np.int64)
                yield element_type, physical, block[:, 3 + tag_count : 3 + tag_count + CORNER_NODES[element_type]]
        return
    for element_type, following, tag_count, position in elements[1]:
        if ELEMENT_DIMENSION[element_type] not in dimensions:
            continue
        width = 1 + tag_count + NODES_PER_TYPE[element_type]
        corners = CORNER_NODES[element_type]
        for first in range(0, following, CHUNK):
            count = min(CHUNK, following - first)
            block = np.frombuffer(data, dtype="<i4", count=count * width, offset=position + first * width * 4).reshape(count, width)
            physical = block[:, 1] if tag_count else np.zeros(count, dtype=np.int32)
            yield element_type, physical, block[:, 1 + tag_count : 1 + tag_count + corners]


def census(path):
    data, coordinates, names, elements = open_mesh(path)
    entries = {}
    for element_type, physical, corners in iter_elements(data, elements):
        dimension = ELEMENT_DIMENSION[element_type]
        for tag in np.unique(physical):
            select = physical == tag
            xyz = coordinates[corners[select]]  # (m, k, 3)
            entry = entries.setdefault((dimension, int(tag)), {
                "Dimension": dimension, "Attribute": int(tag), "Name": names.get((dimension, int(tag))),
                "Count": 0, "Types": {}, "Min": [np.inf] * 3, "Max": [-np.inf] * 3, "Area": 0.0, "NormalSum": [0.0] * 3,
            })
            entry["Count"] += int(select.sum())
            entry["Types"][str(element_type)] = entry["Types"].get(str(element_type), 0) + int(select.sum())
            entry["Min"] = np.minimum(entry["Min"], xyz.reshape(-1, 3).min(0)).tolist()
            entry["Max"] = np.maximum(entry["Max"], xyz.reshape(-1, 3).max(0)).tolist()
            if element_type in TRIANGLE_TYPES:
                normal = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
                area = 0.5 * np.linalg.norm(normal, axis=1)
                entry["Area"] += float(area.sum())
                # Orientation-free: sum |n| components so opposite orientations do not cancel.
                entry["NormalSum"] = (np.asarray(entry["NormalSum"]) + 0.5 * np.abs(normal).sum(0)).tolist()
    for entry in entries.values():
        if entry["Area"] > 0:
            entry["MeanAbsNormal"] = (np.asarray(entry["NormalSum"]) / entry["Area"]).tolist()
        del entry["NormalSum"]
        if entry["Area"] == 0:
            del entry["Area"]
    node_tags = int(np.isfinite(coordinates[:, 0]).sum())
    finite = coordinates[np.isfinite(coordinates[:, 0])]
    return {
        "Mesh": path, "SHA256": sha256(path), "Bytes": int(len(data)), "Nodes": node_tags,
        "BoundingBox": {"Min": finite.min(0).tolist(), "Max": finite.max(0).tolist()},
        "Elements": sorted(entries.values(), key=lambda e: (e["Dimension"], e["Attribute"])),
    }


def extract_triangles(path, attributes, windows):
    data, coordinates, names, elements = open_mesh(path)
    kept, tags, heights = [], [], []
    for element_type, physical, corners in iter_elements(data, elements, dimensions=(2,)):
        if element_type not in TRIANGLE_TYPES:
            continue
        select = np.isin(physical, attributes)
        if not select.any():
            continue
        xyz = coordinates[corners[select]]
        xy = xyz[:, :, :2]
        inside = np.zeros(len(xy), dtype=bool)
        for x0, x1, y0, y1 in windows:
            inside |= (xy[:, :, 0].max(1) >= x0) & (xy[:, :, 0].min(1) <= x1) & (xy[:, :, 1].max(1) >= y0) & (xy[:, :, 1].min(1) <= y1)
        kept.append(xy[inside].astype(np.float32))
        tags.append(physical[select][inside].astype(np.int32))
        heights.append(xyz[inside][:, :, 2].mean(1))  # mean z: the metal plane of the triangle (plots per plane)
    triangles = np.concatenate(kept) if kept else np.zeros((0, 3, 2), np.float32)
    return triangles, (np.concatenate(tags) if tags else np.zeros(0, np.int32)), (np.concatenate(heights) if heights else np.zeros(0))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mesh")
    parser.add_argument("--output", required=True)
    parser.add_argument("--extract", type=int, nargs="+", help="surface attributes whose triangles are written (npz: xy, attribute)")
    parser.add_argument("--window", type=float, nargs=4, action="append", default=[], metavar=("X0", "X1", "Y0", "Y1"),
                        help="keep only triangles overlapping one of these xy windows (default: all)")
    args = parser.parse_args(argv)
    if args.extract:
        windows = args.window or [[-np.inf, np.inf, -np.inf, np.inf]]
        triangles, tags, heights = extract_triangles(args.mesh, args.extract, windows)
        np.savez_compressed(args.output, xy=triangles, attribute=tags, z=heights, windows=np.asarray(windows, dtype=np.float64))
        print(json.dumps({"Triangles": int(len(triangles)), "Output": args.output}))
    else:
        result = census(args.mesh)
        with open(args.output, "w") as target:
            json.dump(result, target, indent=1)
        for entry in result["Elements"]:
            print(entry["Dimension"], entry["Attribute"], entry["Name"], entry["Count"],
                  [round(v, 3) for v in entry["Min"]], [round(v, 3) for v in entry["Max"]],
                  [round(v, 3) for v in entry.get("MeanAbsNormal", [])])


if __name__ == "__main__":
    main()
