#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Geometric location of every prescribed-potential source (hat apex on the box cut
surface) of a coupon, from its trace files and plan-view boundary alone
(gallery-physics-06b locate_sources.py, case-independent).

* box = bounding box of the trace vertices; the z levels are the distinct apex
  heights (bottom face = the lowest, top face = the highest, metal-top level = the
  highest intermediate level, the others substrate / trench levels) - a role
  assignment that holds for ONE process layer with its metal ABOVE the plane (Nz = +1):
  the case's bound signature layers (Pz, Nz) are checked first and a downward layer or
  a second layer stops the location with a recorded scope guard
  (UnsupportedSourceGeometry: DownwardLayers / MultipleLayers) until the roles are
  assigned per layer band;
* metal loops from plan-view-boundary.csv (vertices in loop order): an edge whose two
  endpoints lie on the same box side is a box-boundary (continuation) edge, every
  other edge is a physical metal edge; the metal cross-section meets the cut along
  the continuation edges, a metal EDGE meets the cut at the endpoints of the physical
  edges that lie on the box boundary (junctions);
* the adjacent substrate surface at the cut is the trench (3100+s) unless a bound
  retained-etch.csv says the point is outside every etched loop (un-etched 3000+s).

usage: locate_sources.py --traces DIR --boundary plan-view-boundary.csv --signature mesh-signature.csv
       --out source-locations.csv [--retained-etch retained-etch.csv] [--geometry-out source-geometry.json]
"""
import argparse
import csv
import glob
import json
import math
from pathlib import Path

TOL = 1e-6
NUDGE = 1e-3
# The layer classes the z-level role assignment below does not cover (the ids are the
# recipe scope classes of mesh_stage_contract.RECIPE_SCOPE_SUPPORTED_CLASSES).
SOURCE_GEOMETRY_GUARDS = {
    "DownwardLayers": "a process layer with Nz = -1 puts the metal below its plane: the highest intermediate apex "
                      "level is not the metal top",
    "MultipleLayers": "two or more process layers stack their metal-top / trench levels: the two highest intermediate "
                      "apex levels are not one layer's cross-section"}


class UnsupportedSourceGeometry(ValueError):
    """The source geometry of a case outside the single upward layer the role assignment
    covers (recorded as a ScopeGuard stop, never a silent misclassification)."""

    def __init__(self, guard, layers):
        self.guard = guard
        self.layers = layers
        super().__init__(f"ScopeGuard[{guard}]: {SOURCE_GEOMETRY_GUARDS[guard]} (signature layers (Pz, Nz) {layers}); "
                         f"locate_sources assigns z-level roles for one upward layer only")


def signature_layers(rows):
    """The distinct (Pz, Nz) layers of mesh-signature rows (Nz defaults to +1)."""
    return sorted({(float(row["Pz"]), int(float(row.get("Nz", 1) or 1))) for row in rows})


def check_layers(layers):
    """Fail closed on the layer classes the role assignment does not cover."""
    layers = [(float(plane), int(sign)) for plane, sign in layers]
    if len(layers) > 1:
        raise UnsupportedSourceGeometry("MultipleLayers", layers)
    if any(sign < 0 for _, sign in layers):
        raise UnsupportedSourceGeometry("DownwardLayers", layers)
    return layers


def read_signature_rows(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def read_loops(boundary_path):
    loops = {}
    with open(boundary_path, newline="") as stream:
        for row in csv.DictReader(stream):
            loops.setdefault(int(row["Loop"]), []).append(
                (float(row["X"]), float(row["Y"]), row.get("Class", ""), int(row.get("Conductor", 0) or 0)))
    return loops


def read_etch_loops(path):
    loops = {}
    with open(path, newline="") as stream:
        for row in csv.DictReader(stream):
            loops.setdefault(int(row["Loop"]), []).append((float(row["X"]), float(row["Y"])))
    return loops


def trace_apex(path):
    """(index, (x, y, z)) of the V = 1 vertex of a basis trace file."""
    index = int(Path(path).stem.split("-")[1])
    with open(path, newline="") as stream:
        for row in csv.DictReader(stream):
            if abs(float(row["V"]) - 1.0) < 1e-9:
                return index, (float(row["x"]), float(row["y"]), float(row["z"]))
    raise ValueError(f"{path}: no apex (V = 1) vertex")


def trace_extent(paths):
    xs, ys = [], []
    for path in paths:
        with open(path, newline="") as stream:
            for row in csv.DictReader(stream):
                xs.append(float(row["x"]))
                ys.append(float(row["y"]))
    return (min(xs), max(xs), min(ys), max(ys))


def point_on_segment(p, a, b, tol=TOL):
    (px, py), (ax, ay), (bx, by) = p, a, b
    dx, dy = bx - ax, by - ay
    length2 = dx * dx + dy * dy
    if length2 == 0:
        return math.hypot(px - ax, py - ay) < tol
    t = ((px - ax) * dx + (py - ay) * dy) / length2
    if t < -tol or t > 1 + tol:
        return False
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy)) < tol


def inside(polygon, x, y):
    n = len(polygon)
    crossing = False
    for k in range(n):
        x0, y0 = polygon[k][:2]
        x1, y1 = polygon[(k + 1) % n][:2]
        if (y0 > y) != (y1 > y) and x < (x1 - x0) * (y - y0) / (y1 - y0) + x0:
            crossing = not crossing
    return crossing


def locate(trace_paths, loops, *, etch_loops=None, layers):
    """Rows (one per source) and the geometry summary; `layers` = the case's signature
    layers [(Pz, Nz)] (check_layers: one upward layer, else UnsupportedSourceGeometry)."""
    layers = check_layers(layers)
    box = trace_extent(trace_paths)

    def on_side(x, y):
        sides = []
        if abs(x - box[0]) < TOL:
            sides.append(f"x={box[0]:g}")
        if abs(x - box[1]) < TOL:
            sides.append(f"x={box[1]:g}")
        if abs(y - box[2]) < TOL:
            sides.append(f"y={box[2]:g}")
        if abs(y - box[3]) < TOL:
            sides.append(f"y={box[3]:g}")
        return sides

    def nudged(x, y):
        xi = x + (NUDGE if abs(x - box[0]) < TOL else -NUDGE if abs(x - box[1]) < TOL else 0.0)
        yi = y + (NUDGE if abs(y - box[2]) < TOL else -NUDGE if abs(y - box[3]) < TOL else 0.0)
        return xi, yi

    def metal_at(x, y):
        xi, yi = nudged(x, y)
        for vertices in loops.values():
            if inside(vertices, xi, yi):
                return vertices[0][3]
        return 0

    def etched(x, y):
        if etch_loops is None:
            return True
        xi, yi = nudged(x, y)
        return any(inside(polygon, xi, yi) for polygon in etch_loops.values())

    physical_edges, continuation_edges, junctions = [], [], []
    for vertices in loops.values():
        n = len(vertices)
        for k in range(n):
            (x0, y0, _, conductor), (x1, y1, _, _) = vertices[k], vertices[(k + 1) % n]
            shared = set(on_side(x0, y0)) & set(on_side(x1, y1))
            (continuation_edges if shared else physical_edges).append(((x0, y0), (x1, y1), conductor))
    for a, b, conductor in physical_edges:
        for point in (a, b):
            if on_side(*point):
                junctions.append((point, conductor))
    apexes = dict(trace_apex(path) for path in trace_paths)
    levels = sorted({round(z, 6) for _, _, z in apexes.values()})
    if len(levels) < 2:
        raise ValueError(f"the apexes span a single z level {levels}: no bottom / top face")
    bottom, top = levels[0], levels[-1]
    intermediate = levels[1:-1]
    metal_top = intermediate[-1] if intermediate else None
    # The metal cross-section spans the substrate surface to the metal top: the two
    # highest intermediate levels (the trench floor below them is substrate).
    cross_section_levels = set(intermediate[-2:])
    z_level = {}
    for level in levels:
        if level == bottom:
            z_level[level] = f"bottom face z={level:g}"
        elif level == top:
            z_level[level] = f"top face z={level:g}"
        elif level == metal_top:
            z_level[level] = f"metal-top level z={level:g}"
        else:
            z_level[level] = f"substrate/trench level z={level:g}"
    rows = []
    for index in sorted(apexes):
        x, y, z = apexes[index]
        level = round(z, 6)
        sides = on_side(x, y)
        lateral = ("box vertical-edge/corner" if len(sides) == 2 else f"side face {sides[0]}" if sides else "face interior")
        on_face = level in (bottom, top)
        corner3d = len(sides) == 2 and on_face
        cross_section = (any(point_on_segment((x, y), a, b) for a, b, _ in continuation_edges)
                         and level in cross_section_levels)
        junction_here = bool(sides) and any(math.hypot(p[0] - x, p[1] - y) < TOL for p, _ in junctions)
        junction_distance = min((math.hypot(p[0] - x, p[1] - y) for p, _ in junctions), default=math.inf) if sides else math.inf
        if not sides:
            adjacent = ""
        elif cross_section:
            conductor = metal_at(x, y)
            adjacent = f"metal {5000 + conductor}/{6000 + conductor} cross-section"
        elif not on_face:
            adjacent = "3100 trench" if etched(x, y) else "3000 un-etched"
        else:
            adjacent = ""
        rows.append({"index": index, "x": x, "y": y, "z": z, "z_level": z_level[level],
                     "z_role": ("bottom" if level == bottom else "top" if level == top else
                                "metal-top" if level == metal_top else "substrate"),
                     "lateral": lateral, "box_corner_3d": int(corner3d), "adjacent_surface_at_z0": adjacent,
                     "metal_meets_cut_here": int(cross_section), "metal_edge_junction_here": int(junction_here),
                     "junction_distance_um": f"{junction_distance:.4f}" if math.isfinite(junction_distance) else "",
                     "near_trench_cut_junction": int(bool(sides) and not on_face and level != metal_top)})
    geometry = {"sources": len(rows), "box": box, "z_levels": [z_level[level] for level in levels],
                "physical_metal_edges": physical_edges, "continuation_edges": continuation_edges,
                "metal_edge_cut_junctions": junctions,
                "sources_at_junctions": [r["index"] for r in rows if r["metal_edge_junction_here"]],
                "sources_on_metal_cross_section": [r["index"] for r in rows if r["metal_meets_cut_here"]],
                "box_corner_3d_sources": [r["index"] for r in rows if r["box_corner_3d"]],
                "face_interior_sources": [r["index"] for r in rows if r["lateral"] == "face interior"],
                "per_z_level": {z_level[level]: sum(1 for r in rows if r["z_level"] == z_level[level]) for level in levels},
                "retained_etch_bound": etch_loops is not None, "layers": layers,
                "layer_rule": "one upward layer (Nz = +1): metal-top = the highest intermediate apex level"}
    return rows, geometry


def write_locations(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def locate_directory(traces_dir, boundary_path, out_path, *, signature_path, retained_etch=None, geometry_out=None):
    trace_paths = sorted(glob.glob(str(Path(traces_dir) / "basis-*.csv")))
    if not trace_paths:
        raise ValueError(f"no basis-*.csv trace under {traces_dir}")
    rows, geometry = locate(trace_paths, read_loops(boundary_path),
                            etch_loops=read_etch_loops(retained_etch) if retained_etch else None,
                            layers=signature_layers(read_signature_rows(signature_path)))
    write_locations(rows, out_path)
    if geometry_out:
        Path(geometry_out).write_text(json.dumps(geometry, indent=1) + "\n")
    return rows, geometry


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--traces", required=True)
    parser.add_argument("--boundary", required=True)
    parser.add_argument("--signature", required=True, help="the case's mesh-signature.csv (layer check: Pz, Nz)")
    parser.add_argument("--out", required=True)
    parser.add_argument("--retained-etch")
    parser.add_argument("--geometry-out")
    args = parser.parse_args(argv)
    _, geometry = locate_directory(args.traces, args.boundary, args.out, signature_path=args.signature,
                                   retained_etch=args.retained_etch, geometry_out=args.geometry_out)
    print(json.dumps({key: value for key, value in geometry.items()
                      if key not in ("physical_metal_edges", "continuation_edges")}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
