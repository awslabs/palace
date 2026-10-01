#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Window extracts of a version-2 identification manifest (lane W of decision 184).

    python3 -m surface_response_identification.window_extract \\
        --manifest postpro/surface-response-requirements.json --chip sct002 \\
        --window S1p 290 765 -200 80 [--window ...] [--halo 60] --output-dir extracts/

A chip manifest is 200-400 MB of JSON; every window tool reads it ONCE here (in the background) and works
on the extracts afterwards. For every window the extract is a version-2 mini-manifest (readable by
``window_query``, ``segment_identity`` and ``plot_identification``) holding

* the segments with a portion inside the box grown by ``--halo`` (default 60 um: the plan's exclusion margin
  plus a bump footprint), renumbered, their feature portions re-indexed; ``WindowExtract.SegmentChipIndex``
  maps every extract segment to its chip index;
* the features with a kept portion (portions on dropped segments removed; ``ChipPortionCount`` records the
  chip's count), the whole Arcs / Exclusions / Conventions blocks, MatchingRadius, GeometryDigest;
* ``WindowExtract.Loops``: the chip's perimeter loops (``window_query.PerimeterLoops`` on the WHOLE chip:
  plane, closed, metal side by corner votes, body, parent, depth, the full polygon) whose bounding box meets the
  grown box — the loop / body model of the chip, not of the truncated extract;
* ``WindowExtract.Footprints``: the bump footprints (closed loops of NonManifold segments, chip-wide) whose
  bounding box meets the grown box, each with the metal body it lands on (the innermost closed perimeter
  loop containing its centroid on that plane; the plane's unbounded ground when none) and its partner
  footprint on the other plane (same centroid within 0.5 um; the bump height = the plane difference);
* ``WindowExtract.Bodies``: every body referenced by a kept loop or footprint with Kind (Ground / Island /
  OpenChain), plane and its CHIP conductor = the union of bodies joined by bump footprints anywhere on the
  chip (decision 184 M1: a bump joins the L1 body and the L2 body it lands on into one conductor).

Coordinates are the manifest's (absolute chip coordinates).
"""

import argparse
import collections
import json
import math
import resource
import sys
import time

from . import window_query as WQ

FOOTPRINT_MATCH_TOLERANCE = 0.5  # um between the centroids of the two planes' footprints of one bump


def polygon_area_centroid(polygon):
    area2, cx, cy = 0.0, 0.0, 0.0
    n = len(polygon)
    for i in range(n):
        x0, y0 = polygon[i]
        x1, y1 = polygon[(i + 1) % n]
        cross = x0 * y1 - x1 * y0
        area2 += cross
        cx += (x0 + x1) * cross
        cy += (y0 + y1) * cross
    if abs(area2) < 1e-18:
        xs = [p[0] for p in polygon]
        ys = [p[1] for p in polygon]
        return 0.0, (sum(xs) / n, sum(ys) / n)
    return 0.5 * area2, (cx / (3.0 * area2), cy / (3.0 * area2))


def trace_components(segments, indices):
    """Connected components of the given segments (shared endpoints, one plane each) as ordered polylines.
    Returns dicts Plane / Segments / Polygon / Closed / Area / BBox."""
    by_plane = collections.defaultdict(dict)
    for i in indices:
        z = round(float(segments[i]["Key"][0][2]), 3)
        for p in segments[i]["Key"]:
            by_plane[z].setdefault(WQ.point_key(p), []).append(i)
    seen = set()
    out = []
    for i in indices:
        if i in seen:
            continue
        z = round(float(segments[i]["Key"][0][2]), 3)
        adjacency = by_plane[z]
        stack, component = [i], []
        seen.add(i)
        while stack:
            j = stack.pop()
            component.append(j)
            for p in segments[j]["Key"]:
                for k in adjacency[WQ.point_key(p)]:
                    if k not in seen:
                        seen.add(k)
                        stack.append(k)
        polygon, closed = order_component(segments, component, adjacency)
        area, _ = polygon_area_centroid(polygon)
        xs = [p[0] for p in polygon]
        ys = [p[1] for p in polygon]
        out.append({"Plane": z, "Segments": component, "Polygon": polygon, "Closed": closed, "Area": area,
                    "BBox": [min(xs), max(xs), min(ys), max(ys)]})
    return out


def order_component(segments, component, adjacency):
    """Ordered polyline through a component (the traversal of ``window_query.PerimeterLoops._order``)."""
    comp = set(component)
    start = component[0]
    k0, k1 = WQ.point_key(segments[start]["Key"][0]), WQ.point_key(segments[start]["Key"][1])
    polygon = [k0[:2], k1[:2]]
    used = {start}
    current = k1
    closed = False
    while True:
        nxt = None
        for j in adjacency[current]:
            if j in comp and j not in used:
                nxt = j
                break
        if nxt is None:
            break
        used.add(nxt)
        a, b = WQ.point_key(segments[nxt]["Key"][0]), WQ.point_key(segments[nxt]["Key"][1])
        other = b if a == current else a
        if other == k0 and len(used) == len(comp):
            closed = True
            break
        polygon.append(other[:2])
        current = other
    if len(used) != len(comp):
        closed = False
    return polygon, closed


def innermost_loop(loops, plane, x, y):
    """Index of the smallest-area closed loop on the plane containing (x, y), else None."""
    best, best_area = None, None
    for idx, loop in loops.items():
        if loop["Plane"] != plane or not loop["Closed"]:
            continue
        bb = loop["BBox"]
        if x < bb[0] or x > bb[1] or y < bb[2] or y > bb[3]:
            continue
        area = abs(loop["Area"])
        if best_area is not None and area >= best_area:
            continue
        if WQ._point_in_polygon(x, y, loop["Polygon"]):
            best, best_area = idx, area
    return best


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, a):
        self.parent.setdefault(a, a)
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def body_kind(loops, body_index):
    if body_index < 0:
        return "Ground"
    body = loops[body_index]
    if not body["Closed"]:
        return "OpenChain"
    return "Ground" if body["Parent"] is None else "Island"


def footprint_body(loops, ground_bodies, plane, x, y):
    """Body of the metal under (x, y) on the plane: the innermost closed perimeter loop's body when the metal
    is inside it; a hole -> None (no metal under the footprint: reported); no loop -> the unbounded ground."""
    idx = innermost_loop(loops, plane, x, y)
    if idx is None:
        return ground_bodies.get(("Ground", plane), -1 - len(ground_bodies)), "UnboundedGround"
    loop = loops[idx]
    if loop["MetalInside"]:
        return loop["Body"], "Loop"
    return None, "Hole"


def chip_footprints(segments, loops, ground_bodies):
    """Every closed NonManifold loop of the chip with the body it lands on and its partner on the other plane."""
    indices = [i for i, s in enumerate(segments) if s.get("Exclusion", {}).get("Class") == "NonManifold"]
    footprints = []
    for comp in trace_components(segments, indices):
        area, (cx, cy) = polygon_area_centroid(comp["Polygon"])
        body, how = footprint_body(loops, ground_bodies, comp["Plane"], cx, cy) if comp["Closed"] else (None, "Open")
        footprints.append({"Plane": comp["Plane"], "Polygon": comp["Polygon"], "Closed": comp["Closed"], "Area": abs(area),
                           "Centroid": [cx, cy], "BBox": comp["BBox"], "Segments": comp["Segments"], "Body": body,
                           "BodyBy": how, "Partner": None, "Height": None})
    # partners: same centroid on another plane
    grid = collections.defaultdict(list)
    for k, fp in enumerate(footprints):
        if fp["Closed"]:
            grid[(round(fp["Centroid"][0] / FOOTPRINT_MATCH_TOLERANCE), round(fp["Centroid"][1] / FOOTPRINT_MATCH_TOLERANCE))].append(k)
    for k, fp in enumerate(footprints):
        if not fp["Closed"] or fp["Partner"] is not None:
            continue
        gx, gy = round(fp["Centroid"][0] / FOOTPRINT_MATCH_TOLERANCE), round(fp["Centroid"][1] / FOOTPRINT_MATCH_TOLERANCE)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for m in grid.get((gx + dx, gy + dy), []):
                    other = footprints[m]
                    if m == k or other["Plane"] == fp["Plane"] or other["Partner"] is not None:
                        continue
                    if math.dist(fp["Centroid"], other["Centroid"]) <= FOOTPRINT_MATCH_TOLERANCE:
                        fp["Partner"], other["Partner"] = m, k
                        h = abs(other["Plane"] - fp["Plane"])
                        fp["Height"] = other["Height"] = h
                        break
                if fp["Partner"] is not None:
                    break
            if fp["Partner"] is not None:
                break
    return footprints


def chip_conductors(loops, footprints):
    """Chip-wide conductor of every body: bodies joined by a bump (both footprints of a pair landing on metal)."""
    uf = UnionFind()
    for idx, loop in loops.items():
        uf.find(loop["Body"])
    joins = 0
    for k, fp in enumerate(footprints):
        m = fp["Partner"]
        if m is None or m < k or fp["Body"] is None or footprints[m]["Body"] is None:
            continue
        if uf.find(fp["Body"]) != uf.find(footprints[m]["Body"]):
            joins += 1
        uf.union(fp["Body"], footprints[m]["Body"])
    return uf, joins


def conductor_facts(loops, footprints, conductors):
    """Per chip conductor root: the number of bodies joined into it and whether one of them is a ground (a plane's
    outer boundary or the unbounded ground) - what a window needs to know about the parts of a conductor it does
    not see."""
    facts = {}
    bodies = {loop["Body"] for loop in loops.values()} | {fp["Body"] for fp in footprints if fp["Body"] is not None}
    for body in bodies:
        root = conductors.find(body)
        rec = facts.setdefault(root, {"Bodies": 0, "Ground": False})
        rec["Bodies"] += 1
        rec["Ground"] = rec["Ground"] or body_kind(loops, body) == "Ground"
    return facts


def bbox_meets(bb, box):
    return not (bb[1] < box[0] or bb[0] > box[1] or bb[3] < box[2] or bb[2] > box[3])


def extract_window(manifest, loops_model, footprints, conductors, facts, chip, name, box, halo, note):
    ident = manifest["Identification"]
    segments = ident["Segments"]
    grown = (box[0] - halo, box[1] + halo, box[2] - halo, box[3] + halo)
    keep = []
    for i, s in enumerate(segments):
        if WQ.clip_interval(s["Key"][0], s["Key"][1], grown) is not None:
            keep.append(i)
    index_map = {old: new for new, old in enumerate(keep)}
    new_segments = []
    feature_ids = set()
    for old in keep:
        s = dict(segments[old])
        for s0, s1, fid in s.get("Portions", []):
            feature_ids.add(fid)
        new_segments.append(s)
    features = []
    for f in ident["Features"]:
        if f["Id"] not in feature_ids:
            continue
        g = dict(f)
        portions = f.get("Portions", [])
        g["ChipPortionCount"] = len(portions)
        g["Portions"] = [[index_map[p[0]], p[1], p[2]] for p in portions if p[0] in index_map]
        if "Sides" in f and len(f["Sides"]) == len(portions):
            g["Sides"] = [side for side, p in zip(f["Sides"], portions) if p[0] in index_map]
        features.append(g)
    loops = loops_model.loops
    kept_loops = []
    bodies = {}
    for idx, loop in loops.items():
        if not bbox_meets(loop["BBox"], grown):
            continue
        kept_loops.append({"Index": idx, "Plane": loop["Plane"], "Closed": loop["Closed"], "Area": loop["Area"],
                           "BBox": loop["BBox"], "Depth": loop["Depth"], "Parent": loop["Parent"], "Body": loop["Body"],
                           "MetalInside": loop["MetalInside"], "CornerVotes": loop["CornerVotes"], "Polygon": loop["Polygon"],
                           "Segments": [index_map.get(j, -1) for j in loop["Segments"]]})
        bodies.setdefault(loop["Body"], None)
    kept_footprints = []
    for k, fp in enumerate(footprints):
        if not bbox_meets(fp["BBox"], grown):
            continue
        rec = dict(fp)
        rec["Index"] = k
        rec["Segments"] = [index_map.get(j, -1) for j in fp["Segments"]]
        rec["Inside"] = (box[0] <= fp["BBox"][0] and fp["BBox"][1] <= box[1] and box[2] <= fp["BBox"][2] and fp["BBox"][3] <= box[3])
        rec["Straddles"] = bbox_meets(fp["BBox"], box) and not rec["Inside"]
        kept_footprints.append(rec)
        if fp["Body"] is not None:
            bodies.setdefault(fp["Body"], None)
        if fp["Partner"] is not None and footprints[fp["Partner"]]["Body"] is not None:
            bodies.setdefault(footprints[fp["Partner"]]["Body"], None)
    body_records = {}
    for b in sorted(bodies):
        root = conductors.find(b)
        rec = {"Kind": body_kind(loops, b), "Conductor": root, "ConductorBodies": facts[root]["Bodies"],
               "ConductorGround": facts[root]["Ground"]}
        if b >= 0:
            loop = loops[b]
            rec.update({"Plane": loop["Plane"], "Area": abs(loop["Area"]), "BBox": loop["BBox"], "Closed": loop["Closed"]})
        else:
            plane = [p for (kind, p), v in loops_model.ground_bodies.items() if v == b]
            rec.update({"Plane": plane[0] if plane else None, "Area": None, "BBox": None, "Closed": None})
        body_records[str(b)] = rec
    out = {"Version": manifest.get("Version", 2), "Complete": False,
           "Identification": {"Version": ident.get("Version", 2), "MatchingRadius": ident["MatchingRadius"],
                              "Segments": new_segments, "Features": features, "Arcs": ident.get("Arcs", []),
                              "Vertices": [], "Exclusions": ident.get("Exclusions", []),
                              "Conventions": ident.get("Conventions", {}), "Totals": {},
                              "GeometryDigest": ident.get("GeometryDigest"), "Note": note},
           "WindowExtract": {"Chip": chip, "Window": name, "Box": list(box), "Halo": halo,
                             "SegmentChipIndex": keep, "Loops": kept_loops, "Footprints": kept_footprints,
                             "Bodies": body_records,
                             "GroundBodies": {str(p): v for (k, p), v in loops_model.ground_bodies.items()}}}
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--chip", required=True)
    parser.add_argument("--window", nargs=5, action="append", metavar=("NAME", "X0", "X1", "Y0", "Y1"), required=True)
    parser.add_argument("--halo", type=float, default=60.0)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--summary", default=None, help="JSON summary (loop / footprint / conductor counts)")
    args = parser.parse_args(argv)
    t0 = time.time()

    def log(msg):
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0 ** 2 if sys.platform == "darwin" else 1024.0)
        print(f"[{time.time() - t0:7.1f} s, maxrss {rss:6.0f} MB] {msg}", flush=True)

    manifest = json.load(open(args.manifest))
    log(f"manifest read: {len(manifest['Identification']['Segments'])} segments, {len(manifest['Identification']['Features'])} features")
    ident = manifest["Identification"]
    loops_model = WQ.PerimeterLoops(ident["Segments"], ident["Features"], float(ident["MatchingRadius"]))
    log(f"perimeter loops: {len(loops_model.loops)} ({sum(1 for l in loops_model.loops.values() if l['Closed'])} closed), "
        f"{loops_model.metal_side_checks} corner votes, {len(loops_model.metal_side_disagreements)} disagreements")
    footprints = chip_footprints(ident["Segments"], loops_model.loops, loops_model.ground_bodies)
    paired = sum(1 for fp in footprints if fp["Partner"] is not None)
    log(f"bump footprints: {len(footprints)} NonManifold loops ({sum(1 for f in footprints if f['Closed'])} closed, {paired} paired across planes, "
        f"{sum(1 for f in footprints if f['BodyBy'] == 'Hole')} in a hole)")
    conductors, joins = chip_conductors(loops_model.loops, footprints)
    facts = conductor_facts(loops_model.loops, footprints, conductors)
    roots = set(facts)
    log(f"chip conductors: {len(roots)} after {joins} bump joins of {len({l['Body'] for l in loops_model.loops.values()})} bodies")
    import os
    os.makedirs(args.output_dir, exist_ok=True)
    summary = {"Chip": args.chip, "Manifest": args.manifest, "Loops": len(loops_model.loops),
               "Footprints": len(footprints), "FootprintsPaired": paired, "BumpJoins": joins, "ChipConductors": len(roots),
               "Windows": {}}
    for w in args.window:
        name, box = w[0], tuple(float(v) for v in w[1:])
        note = f"window extract {args.chip} {name} box {list(box)} halo {args.halo} of {args.manifest}"
        out = extract_window(manifest, loops_model, footprints, conductors, facts, args.chip, name, box, args.halo, note)
        path = os.path.join(args.output_dir, f"{args.chip}-{name}.json")
        with open(path, "w") as fh:
            json.dump(out, fh, separators=(",", ":"))
        we = out["WindowExtract"]
        summary["Windows"][name] = {"Box": list(box), "Segments": len(out["Identification"]["Segments"]),
                                    "Features": len(out["Identification"]["Features"]), "Loops": len(we["Loops"]),
                                    "Footprints": len(we["Footprints"]), "Bodies": len(we["Bodies"]), "Path": path}
        log(f"{name}: {summary['Windows'][name]}")
    if args.summary:
        with open(args.summary, "w") as fh:
            json.dump(summary, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
