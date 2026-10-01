#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Plan-view metal polygons of a validation window FROM THE CHIP MESH (lane W of decision 184, task 3).

    python3 -m surface_response_identification.window_polygons --triangles sct002-tris.npz --chip sct002.json \\
        --extract extracts/sct002-S1p.json --e0 e0/sct002-S1p.json --window S1p 290 765 -200 80 \\
        --output polygons/sct002-S1p.json [--verification polygons/sct002-S1p-verification.json]

Input: the window's mesh triangles (``mesh_census.py --extract ... --nodes``: node tags, float64 xyz, attribute, for
every surface attribute meeting the window grown by the halo), the chip description (``--chip``: planes with their
metal-sheet attributes, gap / bump attributes, process, substrate thicknesses), the window extract (``window_extract.py``:
the chip's perimeter loops, bodies and bump-joined conductors) and the window's E0 result (``window_query.py`` on the
extract: the Excitation = the quoted terminal under the rule of decision 180).

Per plane the metal = the triangles of the plane's metal-sheet attributes plus the bump faces lying in the plane
(horizontal faces of the bump attribute: the footprints). Connected components (shared mesh edges) are identified with
the chip's conductors by matching their boundary edges to the extract's segment table (segment -> loop -> body ->
chip conductor). Termination rule (decision 180 / plan n2): every triangle is clipped to the window box; the terminal
conductor of an ``OpenTerminatedTrace`` excitation is clipped to the box shrunk by ``Excitation.OpenSetback`` (the metal
stops 3 R + margin inside the wall); every other NON-ground conductor touching a wall is shorted to ground there by an
explicit metal bridge: the gap (dielectric) triangles in the strip of ``--bridge-width`` along that wall between the
conductor and its neighbouring metal become metal (the strip lines split every triangle so shared edges stay exact).
The union of the metal pieces is taken by edge counting (an edge of exactly one piece is a boundary), the boundary
edges are chained into loops (outer loops counter-clockwise, holes clockwise -> holes are emitted with the OUTER
orientation as lane F's mesher expects), and the result is written in the fabricated-window mesher's polygon-set
schema (``examples/cpw3d_surface/window_validation/README.md``, Version 1): ``ground`` for every grounded /
bridged / non-excited conductor, one label for the terminal, ``Bumps`` = the footprints with their conductor,
``Terminals``. ``WindowPolygons`` holds the extra bookkeeping (chip conductor per polygon, bridges, setback, heights).

Verification (``--verification``): every metal boundary edge not on a wall and outside the bridge / setback strips
must lie on a segment of the chip manifest (the extract) on the same plane within ``--tolerance``; every extract
segment entirely inside the box (uncut) must be covered by boundary edges over its whole length. Unmatched length,
uncovered segments and conductor conflicts are reported; ``Passed`` iff all three are zero.
"""

import argparse
import collections
import json
import math
import sys

import numpy as np

from . import window_query as WQ

DECIMALS = 6  # vertex keys: 1e-6 um


def key(p):
    return (round(float(p[0]), DECIMALS), round(float(p[1]), DECIMALS))


def clip_half_plane(polygon, axis, value, keep_greater):
    """Sutherland-Hodgman against x (axis 0) or y (axis 1) = value, keeping the side >= value (or <= value).
    Intersections are computed from the lexicographically smaller endpoint so shared edges clip identically."""
    out = []
    n = len(polygon)
    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        da, db = a[axis] - value, b[axis] - value
        ina = da >= -1e-12 if keep_greater else da <= 1e-12
        inb = db >= -1e-12 if keep_greater else db <= 1e-12
        if ina:
            out.append(a)
        if ina != inb and abs(da) > 1e-12 and abs(db) > 1e-12:
            out.append(intersection(a, b, axis, value))
        elif ina != inb and (abs(da) <= 1e-12 or abs(db) <= 1e-12):
            pass  # the crossing point is an endpoint already emitted (a) or to be emitted (b on the next step)
    return dedupe(out)


def intersection(a, b, axis, value):
    lo, hi = (a, b) if a < b else (b, a)
    t = (value - lo[axis]) / (hi[axis] - lo[axis])
    p = [lo[0] + t * (hi[0] - lo[0]), lo[1] + t * (hi[1] - lo[1])]
    p[axis] = value
    return (round(p[0], DECIMALS), round(p[1], DECIMALS))


def dedupe(polygon):
    out = []
    for p in polygon:
        p = key(p)
        if not out or out[-1] != p:
            out.append(p)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out


def split_by_line(polygon, axis, value):
    """Both sides of the line (pieces with < 3 vertices dropped)."""
    pieces = []
    for keep_greater in (True, False):
        part = clip_half_plane(polygon, axis, value, keep_greater)
        if len(part) >= 3 and abs(signed_area(part)) > 1e-14:
            pieces.append(part)
    return pieces


def signed_area(polygon):
    a = 0.0
    n = len(polygon)
    for i in range(n):
        x0, y0 = polygon[i]
        x1, y1 = polygon[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def centroid(polygon):
    xs = [p[0] for p in polygon]
    ys = [p[1] for p in polygon]
    return sum(xs) / len(xs), sum(ys) / len(ys)


def point_in_polygon(x, y, polygon):
    return WQ._point_in_polygon(x, y, polygon)


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, a):
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def components(nodes):
    """Connected components of triangles (node tag triples) sharing an edge."""
    uf = UnionFind(len(nodes))
    owner = {}
    for t, tri in enumerate(nodes):
        for i in range(3):
            e = (int(min(tri[i], tri[(i + 1) % 3])), int(max(tri[i], tri[(i + 1) % 3])))
            if e in owner:
                uf.union(owner[e], t)
            else:
                owner[e] = t
    return [uf.find(t) for t in range(len(nodes))]


def boundary_edges(nodes, xy):
    """Directed boundary edges (a, b) of the unclipped triangle set (edges of exactly one triangle), as point pairs."""
    count = collections.Counter()
    for tri in nodes:
        for i in range(3):
            count[(int(min(tri[i], tri[(i + 1) % 3])), int(max(tri[i], tri[(i + 1) % 3])))] += 1
    edges = []
    for t, tri in enumerate(nodes):
        for i in range(3):
            a, b = int(tri[i]), int(tri[(i + 1) % 3])
            if count[(min(a, b), max(a, b))] == 1:
                edges.append((t, tuple(xy[t][i]), tuple(xy[t][(i + 1) % 3])))
    return edges


class SegmentIndex:
    """Extract segments of one plane in a grid for point -> segment lookup."""

    def __init__(self, segments, indices, cell=10.0):
        self.segments = segments
        self.cell = cell
        self.grid = collections.defaultdict(list)
        for i in indices:
            (x0, y0, _), (x1, y1, _) = segments[i]["Key"]
            for gx in range(int(math.floor(min(x0, x1) / cell)), int(math.floor(max(x0, x1) / cell)) + 1):
                for gy in range(int(math.floor(min(y0, y1) / cell)), int(math.floor(max(y0, y1) / cell)) + 1):
                    self.grid[(gx, gy)].append(i)

    def on_segment(self, a, b, tolerance):
        """Index of a segment containing both a and b (collinear, within its extent), else None."""
        mx, my = 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])
        gx, gy = int(math.floor(mx / self.cell)), int(math.floor(my / self.cell))
        best = None
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for i in self.grid.get((gx + dx, gy + dy), []):
                    (x0, y0, _), (x1, y1, _) = self.segments[i]["Key"]
                    if WQ._segment_distance(a[0], a[1], (x0, y0), (x1, y1)) <= tolerance and \
                            WQ._segment_distance(b[0], b[1], (x0, y0), (x1, y1)) <= tolerance:
                        return i
        return best


def on_wall(a, b, box, tol=1e-9):
    for axis, value in ((0, box[0]), (0, box[1]), (1, box[2]), (1, box[3])):
        if abs(a[axis] - value) <= tol and abs(b[axis] - value) <= tol:
            return (axis, value)
    return None


def chain_loops(edges):
    """Closed loops from directed edges (a -> b); a vertex with several outgoing edges takes any unused one (a
    pinch point makes one loop of two touching loops: reported by the caller through the self-touch count)."""
    outgoing = collections.defaultdict(list)
    for a, b in edges:
        outgoing[a].append(b)
    used = collections.Counter()
    loops, open_chains = [], 0
    for a, b in edges:
        if used[(a, b)]:
            continue
        loop = [a]
        current = a
        nxt = b
        used[(a, b)] += 1
        while nxt != a:
            loop.append(nxt)
            candidates = [c for c in outgoing[nxt] if not used[(nxt, c)]]
            if not candidates:
                open_chains += 1
                break
            current, nxt = nxt, candidates[0]
            used[(current, nxt)] += 1
        else:
            loops.append(loop)
    return loops, open_chains


def polygons_with_holes(loops):
    """Outer loops (positive area) with their holes (negative area, inside the smallest enclosing outer loop)."""
    outers = [(signed_area(l), l) for l in loops if signed_area(l) > 0]
    holes = [(signed_area(l), l) for l in loops if signed_area(l) < 0]
    outers.sort(key=lambda t: t[0])  # smallest first
    result = [{"Outer": l, "Holes": [], "Area": a} for a, l in outers]
    unplaced = []
    for a, h in holes:
        x, y = h[0]
        # a hole vertex lies on the hole itself: test an interior point (the midpoint of a tiny inward offset is
        # fragile), so test all vertices by majority against each outer
        placed = False
        for rec in result:
            inside = sum(1 for p in h if point_in_polygon(p[0], p[1], rec["Outer"]) or any(k == p for k in rec["Outer"]))
            if inside >= max(1, (len(h) + 1) // 2):
                rec["Holes"].append(h)
                rec["Area"] += a
                placed = True
                break
        if not placed:
            unplaced.append(h)
    return result, unplaced


def wall_intervals(pieces, box, tol=1e-9):
    """Per wall: [(lo, hi, conductor root)] metal coverage from the piece edges lying on the wall."""
    intervals = collections.defaultdict(list)
    for piece, conductor in pieces:
        n = len(piece)
        for i in range(n):
            a, b = piece[i], piece[(i + 1) % n]
            w = on_wall(a, b, box, tol)
            if w is None:
                continue
            axis = 1 - w[0]
            intervals[w].append((min(a[axis], b[axis]), max(a[axis], b[axis]), conductor))
    merged = {}
    for w, items in intervals.items():
        items.sort()
        out = []
        for lo, hi, c in items:
            if out and lo <= out[-1][1] + tol and out[-1][2] == c:
                out[-1] = (out[-1][0], max(out[-1][1], hi), c)
            else:
                out.append((lo, hi, c))
        merged[w] = out
    return merged


def bridge_rectangles(intervals, box, width, grounded, terminal, warnings):
    """Bridge rectangles [x0, x1, y0, y1] in the gaps next to every non-ground non-terminal wall interval."""
    rectangles = []
    for (axis, value), items in intervals.items():
        lo_wall, hi_wall = (box[2], box[3]) if axis == 0 else (box[0], box[1])
        for k, (lo, hi, c) in enumerate(items):
            if c in grounded or c == terminal or c is None:
                continue
            for side in (-1, 1):
                if side < 0:
                    neighbour = items[k - 1] if k > 0 else None
                    g0, g1 = (neighbour[1] if neighbour else lo_wall), lo
                else:
                    neighbour = items[k + 1] if k + 1 < len(items) else None
                    g0, g1 = hi, (neighbour[0] if neighbour else hi_wall)
                if neighbour is not None and neighbour[2] == terminal:
                    warnings.append(f"conductor {c} is adjacent to the terminal along wall {axis}={value}: not bridged on that side")
                    continue
                if neighbour is None:
                    warnings.append(f"conductor {c} has no metal neighbour towards the box corner along wall {axis}={value}: bridge to the corner")
                if g1 - g0 <= 1e-9:
                    continue
                inner = value + width if value == (box[0] if axis == 0 else box[2]) else value - width
                if axis == 0:
                    rectangles.append([min(value, inner), max(value, inner), g0, g1, c])
                else:
                    rectangles.append([g0, g1, min(value, inner), max(value, inner), c])
    return rectangles


def inside_rectangle(p, r, tol=1e-9):
    return r[0] - tol <= p[0] <= r[1] + tol and r[2] - tol <= p[1] <= r[3] + tol


def bump_columns(footprints, warnings, pairing=0.5):
    """One bump object per column: the per-plane footprints (the bump faces drawn in each plane) paired by centroid
    within ``pairing`` um; the first plane's footprint is emitted, ``Planes`` lists the planes it was seen in and the
    conductors must agree (lane F's mesher needs the footprint on metal of the same conductor on both planes)."""
    columns = []
    for b in footprints:
        cx, cy = centroid(b["Footprint"])
        for col in columns:
            if math.dist((cx, cy), col["Centroid"]) <= pairing and b["Plane"] not in col["Planes"]:
                col["Planes"].append(b["Plane"])
                if b["Conductor"] != col["Conductor"]:
                    warnings.append(f"bump at ({cx:.2f}, {cy:.2f}) lands on {col['Conductor']} in {col['Planes'][0]} and on "
                                    f"{b['Conductor']} in {b['Plane']}")
                    col["ConductorConflict"] = True
                break
        else:
            columns.append({"Conductor": b["Conductor"], "Footprint": b["Footprint"], "Planes": [b["Plane"]], "Centroid": [cx, cy],
                            "ChipConductor": b["ChipConductor"], "Clipped": b["Clipped"], "Vertices": len(b["Footprint"])})
    return columns


def export(args):
    data = np.load(args.triangles, allow_pickle=False)
    nodes, xyz, attribute = data["nodes"], data["xyz"], data["attribute"]
    chip = json.load(open(args.chip))
    extract = json.load(open(args.extract))
    e0 = json.load(open(args.e0))["Windows"][args.window[0]]
    box = tuple(float(v) for v in args.window[1:])
    grown = (box[0] - args.halo, box[1] + args.halo, box[2] - args.halo, box[3] + args.halo)
    meets = (xyz[:, :, 0].max(1) >= grown[0]) & (xyz[:, :, 0].min(1) <= grown[1]) & \
            (xyz[:, :, 1].max(1) >= grown[2]) & (xyz[:, :, 1].min(1) <= grown[3])
    nodes, xyz, attribute = nodes[meets], xyz[meets], attribute[meets]
    loops_model = WQ.PerimeterLoops.from_extract(extract)
    bodies = loops_model.bodies
    segments = extract["Identification"]["Segments"]
    radius = float(extract["Identification"]["MatchingRadius"])
    excitation = e0["Excitation"]
    terminal = excitation["Conductor"]
    setback = float(excitation.get("OpenSetback") or 0.0) if excitation["Kind"] == "OpenTerminatedTrace" else 0.0
    bridge_width = args.bridge_width if args.bridge_width is not None else 3.0 * radius
    tolerance = args.tolerance
    warnings = []

    def conductor_of_body(body):
        rec = bodies.get(body)
        if rec is None:
            return None, True
        return rec["Conductor"], rec["ConductorGround"] or rec["Kind"] == "Ground"

    planes_out, bookkeeping, verification = [], [], {"Planes": {}, "Passed": True}
    bump_attrs = set(chip.get("Bump", []))
    bumps_out = []
    terminal_label = None
    for plane in chip["Planes"]:
        z = float(plane["SurfaceZ"])
        metal_attrs = set(plane["Attributes"])
        in_plane = np.all(np.abs(xyz[:, :, 2] - z) < 1e-6, axis=1)
        is_metal = in_plane & (np.isin(attribute, list(metal_attrs)) | np.isin(attribute, list(bump_attrs)))
        is_gap = in_plane & ~is_metal & ~np.isin(attribute, list(chip.get("Exterior", [])))
        is_bump_face = in_plane & np.isin(attribute, list(bump_attrs))
        metal_index = np.nonzero(is_metal)[0]
        xy_all = xyz[:, :, :2]
        # orient every triangle counter-clockwise so the union boundary is oriented
        tri_xy = []
        for t in range(len(nodes)):
            poly = [key(p) for p in xy_all[t]]
            tri_xy.append(poly if signed_area(poly) > 0 else poly[::-1])
        # components of the metal, identified with the chip's conductors through their boundary edges
        comp_of = components(nodes[metal_index])
        seg_plane = [i for i, s in enumerate(segments) if abs(float(s["Key"][0][2]) - z) < 1e-6]
        index = SegmentIndex(segments, seg_plane)
        comp_votes = collections.defaultdict(collections.Counter)
        for t_local, a, b in boundary_edges(nodes[metal_index], xy_all[metal_index]):
            i = index.on_segment(a, b, tolerance)
            if i is None:
                continue
            loop_index = loops_model.segment_loop.get(i)
            if loop_index is None:
                continue
            body = loops_model.loops[loop_index]["Body"]
            comp_votes[comp_of[t_local]][body] += 1
        comp_conductor, comp_ground = {}, {}
        for c in set(comp_of):
            votes = comp_votes.get(c)
            if not votes:
                comp_conductor[c], comp_ground[c] = None, True
                warnings.append(f"plane {plane['Name']}: metal component {c} has no boundary edge on a chip segment: treated as ground")
                continue
            conductors = collections.Counter()
            grounded = {}
            for body, n in votes.items():
                root, g = conductor_of_body(body)
                conductors[root] += n
                grounded[root] = g
            if len(conductors) > 1:
                warnings.append(f"plane {plane['Name']}: metal component {c} matches several conductors {dict(conductors)}: majority taken")
            root = conductors.most_common(1)[0][0]
            comp_conductor[c], comp_ground[c] = root, grounded[root]
        # clip: the box (the setback box for the terminal trace), then the bridge splitting lines
        full = [(tri_xy[t], comp_conductor[comp_of[k]]) for k, t in enumerate(metal_index)]
        clipped_full = []
        for poly, c in full:
            p = poly
            for axis, value, greater in ((0, box[0], True), (0, box[1], False), (1, box[2], True), (1, box[3], False)):
                p = clip_half_plane(p, axis, value, greater)
                if len(p) < 3:
                    break
            if len(p) >= 3 and abs(signed_area(p)) > 1e-14:
                clipped_full.append((p, c))
        intervals = wall_intervals(clipped_full, box)
        grounded_roots = {comp_conductor[k] for k, g in comp_ground.items() if g and comp_conductor[k] is not None}
        rectangles = bridge_rectangles(intervals, box, bridge_width, grounded_roots, terminal, warnings) if args.bridges else []
        setback_box = (box[0] + setback, box[1] - setback, box[2] + setback, box[3] - setback)
        pieces = []  # (polygon, conductor, metal?)
        pool = [(tri_xy[t], comp_conductor[comp_of[k]], True) for k, t in enumerate(metal_index)]
        pool += [(tri_xy[t], None, False) for t in np.nonzero(is_gap)[0]]
        for poly, c, metal in pool:
            b = setback_box if (metal and c == terminal and setback > 0.0) else box
            p = poly
            for axis, value, greater in ((0, b[0], True), (0, b[1], False), (1, b[2], True), (1, b[3], False)):
                p = clip_half_plane(p, axis, value, greater)
                if len(p) < 3:
                    break
            if len(p) < 3 or abs(signed_area(p)) <= 1e-14:
                continue
            parts = [p]
            for r in rectangles:
                for axis, value in ((0, r[0]), (0, r[1]), (1, r[2]), (1, r[3])):
                    parts = [q for part in parts for q in split_by_line(part, axis, value)]
            for q in parts:
                if metal:
                    pieces.append((q, c, True))
                else:
                    cx, cy = centroid(q)
                    hit = [r for r in rectangles if inside_rectangle((cx, cy), r)]
                    if hit:
                        pieces.append((q, hit[0][4], True))  # bridged gap: metal of the bridged conductor (-> ground)
        # union boundary of the metal pieces
        count = collections.Counter()
        directed = []
        piece_conductor = {}
        for q, c, metal in pieces:
            n = len(q)
            for i in range(n):
                a, b = q[i], q[(i + 1) % n]
                count[(min(a, b), max(a, b))] += 1
                directed.append((a, b, c))
        boundary = [(a, b, c) for a, b, c in directed if count[(min(a, b), max(a, b))] == 1]
        loops, open_chains = chain_loops([(a, b) for a, b, c in boundary])
        if open_chains:
            warnings.append(f"plane {plane['Name']}: {open_chains} boundary chains did not close")
        edge_conductor = {(a, b): c for a, b, c in boundary}
        polys, unplaced = polygons_with_holes(loops)
        if unplaced:
            warnings.append(f"plane {plane['Name']}: {len(unplaced)} hole loops without an enclosing outer loop")
        plane_polygons = []
        for rec in polys:
            votes = collections.Counter(edge_conductor.get((rec["Outer"][i], rec["Outer"][(i + 1) % len(rec["Outer"])]))
                                        for i in range(len(rec["Outer"])))
            root = votes.most_common(1)[0][0]
            is_terminal = terminal is not None and root == terminal
            label = "ground"
            if is_terminal:
                label = f"{'island' if excitation['Kind'] == 'Island' else 'trace'}_{'_'.join(str(b) for b in excitation['Bodies'])}"
                terminal_label = label
            plane_polygons.append({"Conductor": label, "Outer": [list(p) for p in rec["Outer"]],
                                   "Holes": [[list(p) for p in h[::-1]] for h in rec["Holes"]]})  # holes oriented like the outer
            bookkeeping.append({"Plane": plane["Name"], "ChipConductor": root, "Label": label, "Area": rec["Area"],
                                "Vertices": len(rec["Outer"]), "Holes": len(rec["Holes"]),
                                "OnWall": any(on_wall(rec["Outer"][i], rec["Outer"][(i + 1) % len(rec["Outer"])], box) for i in range(len(rec["Outer"])))})
        planes_out.append({"Name": plane["Name"], "SurfaceZ": z, "Facing": plane["Facing"],
                           "SubstrateThickness": plane["SubstrateThickness"], "Polygons": plane_polygons})
        # bump footprints in this plane: union boundary of the bump faces clipped to the box
        bump_pieces = []
        for t in np.nonzero(is_bump_face)[0]:
            p = tri_xy[t]
            for axis, value, greater in ((0, box[0], True), (0, box[1], False), (1, box[2], True), (1, box[3], False)):
                p = clip_half_plane(p, axis, value, greater)
                if len(p) < 3:
                    break
            if len(p) >= 3 and abs(signed_area(p)) > 1e-14:
                bump_pieces.append(p)
        bcount = collections.Counter()
        bdirected = []
        for q in bump_pieces:
            for i in range(len(q)):
                a, b = q[i], q[(i + 1) % len(q)]
                bcount[(min(a, b), max(a, b))] += 1
                bdirected.append((a, b))
        bloops, _ = chain_loops([(a, b) for a, b in bdirected if bcount[(min(a, b), max(a, b))] == 1])
        for l in bloops:
            if signed_area(l) <= 0:
                continue
            cx, cy = centroid(l)
            root = None
            for q, c, metal in pieces:
                if metal and c is not None and point_in_polygon(cx, cy, q):
                    root = c
                    break
            label = "ground" if root is None or root != terminal else terminal_label
            bumps_out.append({"Conductor": label, "Footprint": [list(p) for p in l], "Plane": plane["Name"], "ChipConductor": root,
                              "Clipped": any(on_wall(l[i], l[(i + 1) % len(l)], box) for i in range(len(l)))})
        # verification against the chip manifest (the extract): boundary edges off the walls and outside the
        # bridge / setback strips vs the segments; uncut segments covered
        strips = [r[:4] for r in rectangles]
        if setback > 0.0:
            strips += [[box[0], box[0] + setback, box[2], box[3]], [box[1] - setback, box[1], box[2], box[3]],
                       [box[0], box[1], box[2], box[2] + setback], [box[0], box[1], box[3] - setback, box[3]]]
        covered = collections.defaultdict(float)
        unmatched_length, matched_length, excluded_length, unmatched_examples = 0.0, 0.0, 0.0, []
        for a, b, c in boundary:
            length = math.dist(a, b)
            if on_wall(a, b, box):
                continue
            m = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]))
            if any(inside_rectangle(m, r) for r in strips):
                excluded_length += length
                continue
            i = index.on_segment(a, b, tolerance)
            if i is None:
                unmatched_length += length
                if len(unmatched_examples) < 20:
                    unmatched_examples.append([list(a), list(b)])
            else:
                matched_length += length
                covered[i] += length
        uncovered = []
        uncut = 0
        for i in seg_plane:
            s = segments[i]
            p0, p1 = s["Key"][0], s["Key"][1]
            if not (box[0] <= min(p0[0], p1[0]) and max(p0[0], p1[0]) <= box[1] and box[2] <= min(p0[1], p1[1]) and max(p0[1], p1[1]) <= box[3]):
                continue
            if s.get("Exclusion", {}).get("Class") in ("NonManifold", "NonPlanar"):
                continue  # bump outlines / airbridge feet are not metal-dielectric edges
            # the expected coverage = the segment's length outside the bridge / setback strips (union of t-intervals)
            length = float(s["Length"])
            in_strips = []
            for r in strips:
                t = WQ.clip_interval(p0, p1, r)
                if t is not None and t[1] > t[0]:
                    in_strips.append(t)
            in_strips.sort()
            excluded = 0.0
            cursor = 0.0
            for t0, t1 in in_strips:
                t0 = max(t0, cursor)
                if t1 > t0:
                    excluded += (t1 - t0) * length
                    cursor = t1
            expected = length - excluded
            if expected <= 1e-6:
                continue
            uncut += 1
            if covered.get(i, 0.0) < expected - 1e-6:
                uncovered.append({"Segment": i, "Length": length, "Expected": expected, "Covered": covered.get(i, 0.0), "Key": s["Key"],
                                  "Exclusion": s.get("Exclusion", {}).get("Class")})
        verification["Planes"][plane["Name"]] = {
            "MetalTriangles": int(len(metal_index)), "GapTriangles": int(is_gap.sum()), "Components": len(set(comp_of)),
            "Polygons": len(plane_polygons), "BoundaryEdges": len(boundary), "MatchedLength": matched_length,
            "UnmatchedLength": unmatched_length, "UnmatchedExamples": unmatched_examples,
            "ExcludedLength (bridges / setback)": excluded_length, "UncutSegments": uncut, "UncoveredSegments": uncovered,
            "Bridges": rectangles, "Setback": setback}
        if unmatched_length > 1e-6 or uncovered:
            verification["Passed"] = False
    if terminal is not None and terminal_label is None:
        warnings.append(f"terminal conductor {terminal} of the E0 excitation produced no polygon")
        verification["Passed"] = False
    columns = bump_columns(bumps_out, warnings)
    out = {"Version": 1, "Name": args.window[0], "Box": {"X": [box[0], box[1]], "Y": [box[2], box[3]]},
           "Process": chip.get("Process", {"MetalThickness": 0.1, "Overetch": 0.05}),
           "Planes": [{k: v for k, v in p.items()} for p in planes_out],
           "Bumps": [{"Conductor": b["Conductor"], "Footprint": b["Footprint"]} for b in columns],
           "Vacuum": chip.get("Vacuum", {"Below": 0.0, "Above": 0.0}),
           "Terminals": [terminal_label] if terminal_label else [],
           "WindowPolygons": {"Chip": chip.get("Name"), "Mesh": str(data["mesh"]) if "mesh" in data else None,
                              "Extract": args.extract, "E0": args.e0, "Excitation": excitation, "TerminalChipConductor": terminal,
                              "OpenSetback": setback, "BridgeWidth": bridge_width, "Polygons": bookkeeping,
                              "BumpFootprints": [{k: v for k, v in b.items() if k != "Footprint"} for b in columns],
                              "BumpHeight": chip.get("BumpHeight"), "Warnings": warnings}}
    verification["Warnings"] = warnings
    verification["Passed"] = verification["Passed"] and not any("did not close" in w or "several conductors" in w or "lands on" in w
                                                               for w in warnings)
    return out, verification


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--triangles", required=True, help="npz of mesh_census.py --extract ... --nodes")
    parser.add_argument("--chip", required=True, help="chip description JSON (Planes, Bump, Exterior, Process, Vacuum)")
    parser.add_argument("--extract", required=True, help="window extract of window_extract.py")
    parser.add_argument("--e0", required=True, help="window_query result on the extract (Excitation)")
    parser.add_argument("--window", nargs=5, required=True, metavar=("NAME", "X0", "X1", "Y0", "Y1"))
    parser.add_argument("--bridge-width", type=float, default=None, help="metal bridge width along the wall (default 3 R)")
    parser.add_argument("--no-bridges", dest="bridges", action="store_false")
    parser.add_argument("--tolerance", type=float, default=1e-4, help="segment match tolerance (um)")
    parser.add_argument("--halo", type=float, default=60.0, help="triangles meeting the box grown by this halo are considered (um)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--verification", default=None)
    args = parser.parse_args(argv)
    out, verification = export(args)
    with open(args.output, "w") as fh:
        json.dump(out, fh, separators=(",", ":"))
    if args.verification:
        with open(args.verification, "w") as fh:
            json.dump(verification, fh, indent=1)
    summary = {p: {k: v for k, v in rec.items() if k in ("MetalTriangles", "Components", "Polygons", "MatchedLength", "UnmatchedLength",
                                                          "UncutSegments", "Setback")} | {"Uncovered": len(rec["UncoveredSegments"]), "Bridges": len(rec["Bridges"])}
               for p, rec in verification["Planes"].items()}
    print(json.dumps({"Window": args.window[0], "Passed": verification["Passed"], "Terminals": out["Terminals"], "Bumps": len(out["Bumps"]),
                      "Planes": summary, "Warnings": verification["Warnings"][:10]}, indent=None))
    return 0 if verification["Passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
