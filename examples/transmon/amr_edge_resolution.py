#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Near-edge element sizes of an adapted device mesh written by Palace (SaveAdaptMesh).

Reads an "MFEM NC mesh v1.0/v1.1" tetrahedral mesh (the nonconforming format Palace writes
after every AMR refinement; a conforming "MFEM mesh v1.0" is not supported), extracts the
metal-edge segments from the root boundary faces of the metal attributes (edges belonging to
exactly one metal face) and reports, for the leaf tetrahedra whose centroid lies within each
requested distance band of a metal edge, the count and the minimum / percentiles of the
longest element edge. The result is one JSON object (--json) and a one-line summary, so a
refined mesh never has to leave the machine that solved it. Coordinates are in mesh units
(micrometres for the Palace device configurations, L0 = 1e-6).
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np

TETRAHEDRON = 4
TRIANGLE = 2
TET_EDGES = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])


def read_nc_mesh(path):
    """Leaf tetrahedra (node ids), boundary triangles (attribute, node ids) and node coordinates."""
    lines = Path(path).read_text().split("\n")
    if not lines or not lines[0].startswith("MFEM NC mesh v1."):
        raise ValueError(f"{path}: not an MFEM NC mesh (first line {lines[0][:40]!r})")
    scaled = lines[0].startswith("MFEM NC mesh v1.1")
    index = {}
    for number, line in enumerate(lines):
        if line in ("dimension", "elements", "boundary", "vertex_parents", "coordinates", "root_state", "rank", "nodes"):
            index[line] = number
    if int(lines[index["dimension"] + 1]) != 3:
        raise ValueError(f"{path}: only 3D meshes are supported")

    start = index["elements"] + 1
    count = int(lines[start])
    leaves = []
    for line in lines[start + 1:start + 1 + count]:
        fields = line.split()
        # rank attr geom ref_type nodes/children; "rank attr -1" marks an unused element.
        if len(fields) < 4 or fields[3] != "0":
            continue
        if int(fields[2]) != TETRAHEDRON:
            raise ValueError(f"{path}: leaf element geometry {fields[2]} is not a tetrahedron")
        leaves.append([int(v) for v in fields[4:8]])
    leaves = np.array(leaves, dtype=np.int64).reshape(-1, 4)

    boundary = []
    if "boundary" in index:
        start = index["boundary"] + 1
        count = int(lines[start])
        for line in lines[start + 1:start + 1 + count]:
            fields = line.split()
            if int(fields[1]) == TRIANGLE:
                boundary.append([int(fields[0])] + [int(v) for v in fields[2:5]])
    boundary = np.array(boundary, dtype=np.int64).reshape(-1, 4)

    if "coordinates" in index:
        start = index["coordinates"] + 1
        top_count = int(lines[start])
        space_dim = int(lines[start + 1])
        top = np.array([[float(v) for v in lines[start + 2 + i].split()[:space_dim]] for i in range(top_count)])
    elif "nodes" in index:
        # A mesh with a curvature GridFunction (Palace sets curvature order 1): the H1 P1 nodes
        # list every vertex; MFEM numbers the top-level (original) vertices first, in node order,
        # and refined vertices after them in a space-filling-curve order that is not recorded.
        # Only the top-level block is used; refined vertices come from vertex_parents.
        start = index["nodes"] + 1
        header = {}
        while lines[start].strip() and (":" in lines[start] or lines[start].strip() == "FiniteElementSpace"):
            key, _, value = lines[start].partition(":")
            header[key.strip()] = value.strip()
            start += 1
        if header.get("FiniteElementCollection") != "H1_3D_P1" or header.get("VDim") != "3":
            raise ValueError(f"{path}: mesh nodes are not a linear H1 vector GridFunction ({header})")
        while not lines[start].strip():
            start += 1
        values = []
        for line in lines[start:]:
            if not line.strip() or line.startswith("mfem_mesh_end"):
                break
            values.extend(float(v) for v in line.split())
        values = np.array(values)
        vertex_count = len(values) // 3
        nodes = values.reshape(vertex_count, 3) if header.get("Ordering") == "1" else values.reshape(3, vertex_count).T
        parent_count = int(lines[index["vertex_parents"] + 1]) if "vertex_parents" in index else 0
        top_count = vertex_count - parent_count
        top = nodes[:top_count]
    else:
        raise ValueError(f"{path}: neither a coordinates nor a nodes section")

    parents = np.zeros((0, 3), dtype=np.int64)
    scales = np.zeros(0)
    if "vertex_parents" in index:
        start = index["vertex_parents"] + 1
        count = int(lines[start])
        rows = [line.split() for line in lines[start + 1:start + 1 + count]]
        parents = np.array([[int(r[0]), int(r[1]), int(r[2])] for r in rows], dtype=np.int64)
        scales = np.array([float(r[3]) if scaled and len(r) > 3 else 0.5 for r in rows])

    node_count = int(max(top_count, leaves.max() + 1 if leaves.size else 0, parents[:, 0].max() + 1 if parents.size else 0))
    if parents.size and parents[:, 0].min() < top_count:
        raise ValueError(f"{path}: a top-level vertex has parents; the vertex numbering assumption failed")
    coordinates = np.full((node_count, 3), np.nan)
    coordinates[:top_count] = top
    # A refined vertex is the (scaled) point between its two parents; parents may themselves be
    # refined vertices, so resolve in passes until every listed vertex has coordinates.
    pending = np.ones(len(parents), dtype=bool)
    while pending.any():
        rows = parents[pending]
        ready = ~np.isnan(coordinates[rows[:, 1], 0]) & ~np.isnan(coordinates[rows[:, 2], 0])
        if not ready.any():
            raise ValueError(f"{path}: unresolvable vertex parents")
        s = scales[pending][ready][:, None]
        coordinates[rows[ready, 0]] = (1.0 - s) * coordinates[rows[ready, 1]] + s * coordinates[rows[ready, 2]]
        indices = np.flatnonzero(pending)
        pending[indices[ready]] = False
    return leaves, boundary, coordinates


def metal_edge_segments(boundary, metal_attributes, coordinates):
    """Mesh edges shared by a metal boundary triangle and a non-metal boundary triangle (the sheet
    boundary; Palace cracks the interior metal sheet, so a two-sided count cannot separate its edges),
    excluding edges on a face of the mesh bounding box (a sheet cut by the domain truncation)."""
    is_metal = np.isin(boundary[:, 0], list(metal_attributes))
    if not is_metal.any() or is_metal.all():
        raise ValueError("the boundary must carry metal and non-metal attributes")

    def edge_set(triangles):
        edges = np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]])
        edges.sort(axis=1)
        return np.unique(edges, axis=0)

    metal_edges = edge_set(boundary[is_metal, 1:4])
    other_edges = edge_set(boundary[~is_metal, 1:4])
    both = np.concatenate([metal_edges, other_edges])
    unique, counts = np.unique(both, axis=0, return_counts=True)
    segments = unique[counts == 2]
    finite = coordinates[~np.isnan(coordinates[:, 0])]
    lo, hi = finite.min(axis=0), finite.max(axis=0)
    scale = np.maximum(hi - lo, 1.0e-300)
    a, b = coordinates[segments[:, 0]], coordinates[segments[:, 1]]
    on_face = np.zeros(len(segments), dtype=bool)
    for bound in (lo, hi):
        both_ends = (np.abs(a - bound) <= 1.0e-9 * scale) & (np.abs(b - bound) <= 1.0e-9 * scale)
        on_face |= both_ends.any(axis=1)
    return segments[~on_face], segments[on_face]


def point_segment_distances(points, a, b):
    """Distance from every point to the nearest of the segments a[k]-b[k] (dense, chunked)."""
    result = np.full(len(points), np.inf)
    ab = b - a
    ab2 = np.einsum("ij,ij->i", ab, ab)
    step = max(1, 2_000_000 // max(1, len(a)))
    for start in range(0, len(points), step):
        p = points[start:start + step]
        ap = p[:, None, :] - a[None, :, :]
        t = np.clip(np.einsum("ikj,kj->ik", ap, ab) / ab2[None, :], 0.0, 1.0)
        d = ap - t[:, :, None] * ab[None, :, :]
        result[start:start + step] = np.sqrt(np.einsum("ikj,ikj->ik", d, d).min(axis=1))
    return result


def near_edge_distances(centroids, segments, coordinates, cutoff):
    """Centroid-to-metal-edge distance for centroids within cutoff (inf elsewhere), bucketed in the sheet plane."""
    a = coordinates[segments[:, 0]]
    b = coordinates[segments[:, 1]]
    result = np.full(len(centroids), np.inf)
    lo = np.minimum(a, b).min(axis=0) - cutoff
    hi = np.maximum(a, b).max(axis=0) + cutoff
    candidate = np.all((centroids >= lo) & (centroids <= hi), axis=1)
    cell = 4.0 * cutoff
    seg_lo = np.floor((np.minimum(a, b)[:, :2] - cutoff) / cell).astype(np.int64)
    seg_hi = np.floor((np.maximum(a, b)[:, :2] + cutoff) / cell).astype(np.int64)
    buckets = {}
    for k in range(len(a)):
        for i in range(seg_lo[k, 0], seg_hi[k, 0] + 1):
            for j in range(seg_lo[k, 1], seg_hi[k, 1] + 1):
                buckets.setdefault((i, j), []).append(k)
    indices = np.flatnonzero(candidate)
    cells = np.floor(centroids[indices, :2] / cell).astype(np.int64)
    keys = cells[:, 0] * 2_000_003 + cells[:, 1]
    order = np.argsort(keys, kind="stable")
    keys = keys[order]
    indices = indices[order]
    boundaries = np.flatnonzero(np.diff(keys)) + 1
    for group in np.split(np.arange(len(keys)), boundaries):
        if group.size == 0:
            continue
        i, j = cells[order[group[0]]]
        segs = buckets.get((int(i), int(j)))
        if not segs:
            continue
        pts = indices[group]
        result[pts] = point_segment_distances(centroids[pts], a[segs], b[segs])
    result[result > cutoff] = np.inf
    return result


def band_statistics(sizes, distances, radius):
    inside = distances <= radius
    if not inside.any():
        return {"Radius": radius, "Elements": 0}
    s = sizes[inside]
    return {"Radius": radius, "Elements": int(inside.sum()), "MinLongestEdge": float(s.min()),
            "P05LongestEdge": float(np.percentile(s, 5)), "MedianLongestEdge": float(np.median(s)),
            "MaxLongestEdge": float(s.max())}


def analyze(path, metal_attributes, bands, cutoff):
    leaves, boundary, coordinates = read_nc_mesh(path)
    segments, truncation = metal_edge_segments(boundary, metal_attributes, coordinates)
    if len(segments) == 0:
        raise ValueError(f"{path}: no metal edge segments away from the bounding box")
    seg_lengths = np.linalg.norm(coordinates[segments[:, 0]] - coordinates[segments[:, 1]], axis=1)
    truncation_length = float(np.linalg.norm(coordinates[truncation[:, 0]] - coordinates[truncation[:, 1]], axis=1).sum())
    vertices = coordinates[leaves]
    edge_vectors = vertices[:, TET_EDGES[:, 0], :] - vertices[:, TET_EDGES[:, 1], :]
    edge_lengths = np.linalg.norm(edge_vectors, axis=2)
    longest = edge_lengths.max(axis=1)
    shortest = edge_lengths.min(axis=1)
    centroids = vertices.mean(axis=1)
    distances = near_edge_distances(centroids, segments, coordinates, cutoff)
    touching = distances <= longest / 2.0
    result = {
        "Mesh": str(path), "LeafElements": int(len(leaves)), "MetalAttributes": sorted(metal_attributes),
        "MetalEdgeSegments": int(len(segments)), "MetalEdgeLength": float(seg_lengths.sum()),
        "TruncationSegments": int(len(truncation)), "TruncationLength": truncation_length,
        "GlobalMinLongestEdge": float(longest.min()), "GlobalMinShortestEdge": float(shortest.min()),
        "Bands": [band_statistics(longest, distances, r) for r in bands],
        "TouchingEdge": {"Elements": int(touching.sum()),
                         "MinLongestEdge": float(longest[touching].min()) if touching.any() else None,
                         "MedianLongestEdge": float(np.median(longest[touching])) if touching.any() else None},
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", type=Path, required=True, help="MFEM NC mesh written by Palace SaveAdaptMesh")
    parser.add_argument("--metal-attributes", type=int, nargs="+", default=[5, 6, 7, 9])
    parser.add_argument("--band", type=float, nargs="+", default=[2.0, 0.5, 0.1],
                        help="centroid-to-edge distance bands (mesh units)")
    parser.add_argument("--json", type=Path, help="write the result object here")
    args = parser.parse_args()
    result = analyze(args.mesh, set(args.metal_attributes), sorted(args.band, reverse=True), max(args.band))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=1) + "\n")
    bands = "; ".join(f"d<={b['Radius']:g}: n {b['Elements']}" + (f" min {b['MinLongestEdge']:.4g} med {b['MedianLongestEdge']:.4g}" if b["Elements"] else "")
                      for b in result["Bands"])
    print(f"{args.mesh.name}: {result['LeafElements']} leaves, {result['MetalEdgeSegments']} edge segments "
          f"({result['MetalEdgeLength']:.1f}), global min longest edge {result['GlobalMinLongestEdge']:.4g}; {bands}; "
          f"touching edge: n {result['TouchingEdge']['Elements']} min {result['TouchingEdge']['MinLongestEdge']} med {result['TouchingEdge']['MedianLongestEdge']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
