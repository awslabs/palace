#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Physical-consistency check of a version-2 identification manifest (decision 78).

The partition / determinism gates are bookkeeping: they cannot see a feature whose class
contradicts the geometry around it. This check samples every longitudinal feature portion and
asks the geometry directly, with a uniform grid over the perimeter segments (chip scale):

* `IsolatedEdge` / `CurvedEdge` portions: the length whose sample points face another
  (non-excluded) perimeter edge within 2R — closest point q of a segment with |q - p| < 2R and
  ACROSS the chain (|cos(q - p, tangent)| < 0.5). An isolated edge with metal across within the
  interaction radius is not isolated (the narrow-CPW misclassification of DS-SCT-001).
* Pair and stack portions (`SameConductorGap`, `DifferentConductorGap`, `SameConductorStrip`,
  `ParallelEdgeCluster`, the curved classes): the length whose sample points face a THIRD edge
  within 2R — a segment that belongs to none of the feature's sides. A pair with a third edge in
  reach is a multi-edge cross-section (ground - gap - trace - gap - ground narrower than 2R).
Distances are three-dimensional (the two metal planes of a flip chip do not face each other
when they are more than 2R apart); sites and boxes are in the plan view (x, y).

Hard gates (decision 82(2), 2026-09-25): the facing length is a gate, not a diagnostic — it
must be zero apart from the RECORDED exclusions, each reported with its length:
* `AtExactly2R`: the facing distance is exactly 2R on the decision grid (the strict rule: no
  interaction);
* `ClusterNeighbour`: the facing portion belongs to a `SpatialEdgeCluster` — the cluster region
  is the union of radius-R balls around the event cores (its claim radius R, decision 82(1)), so
  metal between R and 2R of a cluster's edge faces it without being inside it;
* `VertexNeighbour`: the facing portion belongs to a vertex feature (a corner window, a rounded
  corner's arc, an endpoint or junction window): the vertex describes its own neighbourhood;
* `ThroughVertex`: the sample and its facing point lie within 2R of one vertex feature (the
  arms of a corner or junction meet through it: design (b) 1, `ThroughVertexZoneOverR`);
* `SelfNeighbourhood`: the facing point lies on the sample's own chain less than pi R of arc
  length away (the bottom of a U of radius >= R faces itself within 2R without a fold:
  `SelfPairNeighbourhoodOverR`).

Output: per class the sampled length, the facing length and its fraction; the flagged samples
grouped into sites (samples within `--site-radius` of each other, default 25 R) ranked by
flagged length, each with its bounding box, feature ids and the closest facing distance — the
plot windows for a visual review.

    python3 -m surface_response_identification.facing_check MANIFEST --output facing.json [--spacing 0.5]
"""
import argparse
import json
import math
from collections import defaultdict

import numpy as np

ISOLATED_CLASSES = ("IsolatedEdge", "CurvedEdge")
PAIR_CLASSES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip",
                "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip",
                "ParallelEdgeCluster", "CurvedParallelEdgeCluster")
VERTEX_CLASSES = ("ConvexCorner", "ConcaveCorner", "Endpoint", "Junction")
EXCLUSION_CLASSES = ("AtExactly2R", "ClusterNeighbour", "VertexNeighbour", "ThroughVertex", "SelfNeighbourhood")
THROUGH_VERTEX_ZONE_OVER_R = 2.0
# A chain's points closer than pi R along the chain are within 2R of each other along any bend
# of radius >= R (Schur): the local neighbourhood that never pairs with itself
# (Identification.Conventions SelfPairNeighbourhoodOverR).
SELF_PAIR_NEIGHBOURHOOD_OVER_R = math.pi


class SegmentGrid:
    """Uniform grid over 2D segments: candidates for a point are the segments crossing the 3 x 3
    cells around it, so a query radius up to the cell size is exact."""

    def __init__(self, p0, p1, cell):
        p0, p1 = p0[:, :2], p1[:, :2]  # plan view
        self.p0, self.p1, self.cell = p0, p1, cell
        self.origin = np.minimum(p0, p1).min(0)
        lower = np.floor((np.minimum(p0, p1) - self.origin) / cell).astype(np.int64)
        upper = np.floor((np.maximum(p0, p1) - self.origin) / cell).astype(np.int64)
        keys, values = [], []
        for i, (lo, hi) in enumerate(zip(lower, upper)):
            for cx in range(lo[0], hi[0] + 1):
                for cy in range(lo[1], hi[1] + 1):
                    keys.append((cx << 32) + cy)
                    values.append(i)
        keys = np.asarray(keys, dtype=np.int64)
        values = np.asarray(values, dtype=np.int64)
        order = np.argsort(keys, kind="stable")
        self.keys, self.values = keys[order], values[order]
        self.unique, self.starts = np.unique(self.keys, return_index=True)
        self.ends = np.append(self.starts[1:], len(self.keys))

    def candidates(self, points):
        """(sample index, segment index) pairs for the 3 x 3 cells around every point."""
        base = np.floor((points[:, :2] - self.origin) / self.cell).astype(np.int64)
        sample_indices, segment_indices = [], []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                keys = ((base[:, 0] + dx) << 32) + (base[:, 1] + dy)
                position = np.searchsorted(self.unique, keys)
                position = np.minimum(position, len(self.unique) - 1)
                found = self.unique[position] == keys
                if not found.any():
                    continue
                starts, ends = self.starts[position[found]], self.ends[position[found]]
                counts = ends - starts
                sample_indices.append(np.repeat(np.nonzero(found)[0], counts))
                offsets = np.arange(counts.sum()) - np.repeat(np.cumsum(counts) - counts, counts)
                segment_indices.append(self.values[np.repeat(starts, counts) + offsets])
        if not sample_indices:
            return np.zeros(0, np.int64), np.zeros(0, np.int64)
        return np.concatenate(sample_indices), np.concatenate(segment_indices)


def load(manifest_path):
    if isinstance(manifest_path, dict):
        manifest = manifest_path
    else:
        with open(manifest_path) as source:
            manifest = json.load(source)
    identification = manifest["Identification"]
    segments = identification["Segments"]
    # Full 3D coordinates: the distances are measured in space, so the metal planes of a
    # flip chip (L1 at z = 0, L2 at z = 4.8 um: 2.4 R apart) never face each other; the grid
    # and the site boxes use the plan-view (x, y) projection.
    p0 = np.array([s["Key"][0][:3] for s in segments], dtype=float)
    p1 = np.array([s["Key"][1][:3] for s in segments], dtype=float)
    length = np.array([s["Length"] for s in segments], dtype=float)
    excluded = np.array(["Exclusion" in s for s in segments])
    radius = float(identification.get("MatchingRadius") or manifest.get("MatchingRadius") or 2.0)
    return identification, p0, p1, length, excluded, radius


# The classifier's decision grid (Identification.Conventions LengthQuantumOverR): an edge
# interacts iff its distance is below 2R on the quantized grid, strict less — a CPW gap or
# trace of exactly 2R (DS-SCT-001 / DS-SCT-002: 4 um at R = 2 um) does not interact and is
# not a defect; its length is reported separately (AtExactly2R).
LENGTH_QUANTUM_OVER_R = 1.0e-8


def facing_samples(grid, p0, p1, excluded, points, tangents, owners, own_mask, radius):
    """For every sample point: the closest facing (across) non-excluded segment whose distance is
    at most 2R (within the decision quantum) and neither the sample's own segment nor in own_mask
    (the feature's own sides); -1 when none; the foot on that segment. The caller separates the
    strictly interacting samples (quantized distance < 2R) from those at exactly 2R."""
    sample_indices, segment_indices = grid.candidates(points)
    a, b = p0[segment_indices], p1[segment_indices]
    ab = b - a
    p = points[sample_indices]
    t = np.clip(((p - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-30), 0.0, 1.0)
    foot = a + t[:, None] * ab
    d = foot - p
    r = np.sqrt((d * d).sum(1))
    cosine = np.abs((d * tangents[sample_indices]).sum(1)) / np.maximum(r, 1e-30)
    ok = ((r <= 2 * radius + 0.5 * LENGTH_QUANTUM_OVER_R * radius) & (r > 1e-9) & (cosine < 0.5) & ~excluded[segment_indices]
          & (segment_indices != owners[sample_indices]) & ~own_mask[segment_indices])
    best = np.full(len(points), -1, dtype=np.int64)
    best_r = np.full(len(points), np.inf)
    best_foot = np.zeros((len(points), 3))
    if ok.any():
        si, sj, ri, fi = sample_indices[ok], segment_indices[ok], r[ok], foot[ok]
        order = np.lexsort((ri, si))
        si, sj, ri, fi = si[order], sj[order], ri[order], fi[order]
        first = np.concatenate([[True], si[1:] != si[:-1]])
        best[si[first]] = sj[first]
        best_r[si[first]] = ri[first]
        best_foot[si[first]] = fi[first]
    return best, best_r, best_foot


def portion_samples(p0, p1, length, portions, spacing):
    """Sample points, tangents, weights (length per sample) and segment ids of feature portions."""
    points, tangents, weights, owners = [], [], [], []
    for seg, s0, s1 in portions:
        total = length[seg]
        if total <= 0 or s1 <= s0:
            continue
        tangent = (p1[seg] - p0[seg]) / total
        # Samples at the midpoints of `count` equal sub-intervals (never on a portion end,
        # where the claim changes hands).
        count = max(1, int(math.ceil((s1 - s0) / spacing)))
        s = s0 + (s1 - s0) * (np.arange(count) + 0.5) / count
        points.append(p0[seg] + tangent[None, :] * s[:, None])
        tangents.append(np.repeat(tangent[None, :], count, 0))
        weights.append(np.full(count, (s1 - s0) / count))
        owners.append(np.full(count, seg, dtype=np.int64))
    if not points:
        return None
    return np.concatenate(points), np.concatenate(tangents), np.concatenate(weights), np.concatenate(owners)


def group_sites(flags, site_radius):
    """Union-find of flagged samples within site_radius (grid hashing); sites sorted by length."""
    if not flags:
        return []
    xy = np.array([f["Point"][:2] for f in flags])
    cell = site_radius
    cells = defaultdict(list)
    for i, (cx, cy) in enumerate(np.floor(xy / cell).astype(np.int64)):
        cells[(int(cx), int(cy))].append(i)
    parent = list(range(len(flags)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for (cx, cy), members in cells.items():
        neighbours = [m for dx in (-1, 0, 1) for dy in (-1, 0, 1) for m in cells.get((cx + dx, cy + dy), [])]
        for i in members:
            for j in neighbours:
                if j > i and np.hypot(*(xy[i] - xy[j])) <= site_radius:
                    parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(len(flags)):
        groups[find(i)].append(i)
    sites = []
    for members in groups.values():
        rows = [flags[i] for i in members]
        pts = xy[members]
        closest = min(rows, key=lambda r: r["Distance"])
        sites.append({
            "MinDistancePoint": closest["Point"][:2],
            "Length": float(sum(r["Weight"] for r in rows)),
            "Samples": len(rows),
            "Classes": sorted({r["Class"] for r in rows}),
            "Features": sorted({r["Feature"] for r in rows}),
            "FacingSegments": sorted({r["FacingSegment"] for r in rows}),
            "MinDistance": float(min(r["Distance"] for r in rows)),
            "Box": [float(pts[:, 0].min()), float(pts[:, 0].max()), float(pts[:, 1].min()), float(pts[:, 1].max())],
            "Center": [float(pts[:, 0].mean()), float(pts[:, 1].mean())],
        })
    sites.sort(key=lambda s: -s["Length"])
    return sites


class ClaimLookup:
    """Feature claiming a point of a segment (portion by run parameter) and the vertex-feature
    points (corner / endpoint / junction vertices, the virtual corners of rounded corners)."""

    def __init__(self, identification, p0, p1, length):
        self.feature_type = {f["Id"]: f["Type"] for f in identification["Features"]}
        self.portions = {}
        for i, s in enumerate(identification["Segments"]):
            if "Portions" in s:
                self.portions[i] = [(float(a), float(b), int(f)) for a, b, f in s["Portions"]]
        self.p0, self.p1, self.length = p0, p1, length
        self.chain_of = {i: s.get("Chain") for i, s in enumerate(identification["Segments"])}
        self.chain_position, self.chain_length, self.chain_closed = self._chain_positions(identification["Segments"])
        # Vertex-feature points: the feature frame origins (the vertex, or the virtual corner
        # of a rounded corner; the vertex table carries mesh vertex indices, not coordinates).
        points = [f["Frame"]["Origin"][:3] for f in identification["Features"] if f["Type"] in VERTEX_CLASSES and f.get("Frame")]
        self.vertex_points = np.array(points, dtype=float).reshape(-1, 3)
        self.tree = None
        if len(self.vertex_points):
            try:
                from scipy.spatial import cKDTree
                self.tree = cKDTree(self.vertex_points)
            except ImportError:
                self.tree = None

    def _chain_positions(self, segments):
        """Arc-length position of every segment's Key[0] along its chain (segments chained by
        their endpoints, from an end for an open chain, from the first segment for a loop),
        the chain lengths and whether the chain is closed."""
        by_chain = defaultdict(list)
        for i, s in enumerate(segments):
            if s.get("Chain") is not None:
                by_chain[s["Chain"]].append(i)
        position, lengths, closed = {}, {}, {}

        def key(p):
            return tuple(round(float(v), 6) for v in p[:3])

        for chain, members in by_chain.items():
            ends = defaultdict(list)
            for i in members:
                ends[key(segments[i]["Key"][0])].append(i)
                ends[key(segments[i]["Key"][1])].append(i)
            starts = [k for k, v in ends.items() if len(v) == 1]
            closed[chain] = not starts
            current = key(segments[members[0]]["Key"][0]) if not starts else min(starts)
            seen = set()
            x = 0.0
            while True:
                nxt = [i for i in ends.get(current, []) if i not in seen]
                if not nxt:
                    break
                i = nxt[0]
                seen.add(i)
                a, b = key(segments[i]["Key"][0]), key(segments[i]["Key"][1])
                forward = a == current
                position[i] = (x, forward)
                x += float(segments[i]["Length"])
                current = b if forward else a
            for i in members:
                if i not in position:
                    position[i] = (x, True)  # disconnected remainder (not expected)
                    x += float(segments[i]["Length"])
            lengths[chain] = x
        return position, lengths, closed

    def arc_distance(self, segment_a, point_a, segment_b, point_b):
        """Arc length along the common chain between two points, or None on different chains."""
        ca, cb = self.chain_of.get(int(segment_a)), self.chain_of.get(int(segment_b))
        if ca is None or ca != cb:
            return None

        def at(seg, point):
            x, forward = self.chain_position[seg]
            s = float(np.linalg.norm(np.asarray(point) - self.p0[seg]))
            return x + (s if forward else self.length[seg] - s)

        d = abs(at(int(segment_a), point_a) - at(int(segment_b), point_b))
        if self.chain_closed.get(ca):
            d = min(d, self.chain_length[ca] - d)
        return d

    def feature_at(self, segment, point):
        """Feature id claiming the point of the segment (by its run parameter), or None."""
        portions = self.portions.get(int(segment))
        if not portions:
            return None
        s = float(np.linalg.norm(np.asarray(point) - self.p0[segment]))
        best = None
        for a, b, f in portions:
            if a - 1.0e-9 <= s <= b + 1.0e-9:
                return f
            if best is None or min(abs(s - a), abs(s - b)) < best[0]:
                best = (min(abs(s - a), abs(s - b)), f)
        return best[1] if best else None

    def through_vertex(self, point, facing_point, reach):
        """A vertex-feature point within reach of both points."""
        if self.tree is None:
            if not len(self.vertex_points):
                return False
            d_a = np.linalg.norm(self.vertex_points - np.asarray(point), axis=1)
            d_b = np.linalg.norm(self.vertex_points - np.asarray(facing_point), axis=1)
            return bool(np.any((d_a <= reach) & (d_b <= reach)))
        near = self.tree.query_ball_point(np.asarray(point), reach)
        if not near:
            return False
        d_b = np.linalg.norm(self.vertex_points[near] - np.asarray(facing_point), axis=1)
        return bool(np.any(d_b <= reach))


def classify_flag(lookup, flag, radius):
    """Recorded exclusion class of a flagged sample, or None (a genuine defect)."""
    seg = flag["FacingSegment"]
    foot = np.asarray(flag["FacingPoint"])
    feature = lookup.feature_at(seg, foot)
    ftype = lookup.feature_type.get(feature) if feature is not None else None
    if ftype == "SpatialEdgeCluster":
        return "ClusterNeighbour"
    if ftype in VERTEX_CLASSES:
        return "VertexNeighbour"
    if lookup.through_vertex(flag["Point"], foot, THROUGH_VERTEX_ZONE_OVER_R * radius):
        return "ThroughVertex"
    arc = lookup.arc_distance(flag["Segment"], flag["Point"], seg, foot)
    if arc is not None and arc < SELF_PAIR_NEIGHBOURHOOD_OVER_R * radius:
        return "SelfNeighbourhood"
    return None


def facing_check(manifest_path, spacing=0.5, site_radius_over_r=25.0, batch=200000):
    identification, p0, p1, length, excluded, radius = load(manifest_path)
    features = identification["Features"]
    grid = SegmentGrid(p0, p1, cell=2.0 * radius)
    lookup = ClaimLookup(identification, p0, p1, length)
    segments_of_chain = defaultdict(list)
    for i, s in enumerate(identification["Segments"]):
        if "Chain" in s:
            segments_of_chain[s["Chain"]].append(i)
    chain_of = {i: s.get("Chain") for i, s in enumerate(identification["Segments"])}

    def own_chain_segments(feature):
        chains = {chain_of[seg] for seg, _, _ in feature["Portions"]}
        return [i for c in chains if c is not None for i in segments_of_chain[c]]
    own_mask = np.zeros(len(p0), dtype=bool)
    totals = defaultdict(float)
    facing = defaultdict(float)
    at_2r = defaultdict(float)  # facing length at exactly 2R on the decision grid (not interacting)
    excluded_length = defaultdict(lambda: defaultdict(float))  # class -> exclusion -> length
    unexcluded = defaultdict(float)
    counts = defaultdict(int)
    quantum = LENGTH_QUANTUM_OVER_R * radius
    flagged_features = defaultdict(set)
    flags = []
    for feature in features:
        ftype = feature["Type"]
        if ftype not in ISOLATED_CLASSES + PAIR_CLASSES:
            continue
        sampled = portion_samples(p0, p1, length, feature.get("Portions", []), spacing)
        if sampled is None:
            continue
        points, tangents, weights, owners = sampled
        counts[ftype] += 1
        totals[ftype] += float(weights.sum())
        # A pair's / stack's own sides are its member EDGES (the chains of its portions): the
        # same chain's portions claimed by a neighbouring feature (the curved section of the
        # same pair, the next stack of a taper) are not a third edge.
        own_segments = own_chain_segments(feature) if ftype in PAIR_CLASSES else []
        own_mask[own_segments] = True
        for first in range(0, len(points), batch):
            sl = slice(first, first + batch)
            best, best_r, best_foot = facing_samples(grid, p0, p1, excluded, points[sl], tangents[sl], owners[sl], own_mask, radius)
            found = best >= 0
            # Strict less than 2R on the classifier's quantized grid; the rest is exactly 2R.
            hit = found & (np.round(best_r / quantum) < np.round(2.0 * radius / quantum))
            at_2r[ftype] += float(weights[sl][found & ~hit].sum())
            excluded_length[ftype]["AtExactly2R"] += float(weights[sl][found & ~hit].sum())
            if not hit.any():
                continue
            facing[ftype] += float(weights[sl][hit].sum())
            for k in np.nonzero(hit)[0]:
                flag = {"Class": ftype, "Feature": feature["Id"], "Segment": int(owners[sl][k]), "FacingSegment": int(best[k]),
                        "Distance": float(best_r[k]), "Weight": float(weights[sl][k]), "Point": [float(v) for v in points[sl][k]],
                        "FacingPoint": [float(v) for v in best_foot[k]]}
                flag["Exclusion"] = classify_flag(lookup, flag, radius)
                if flag["Exclusion"] is None:
                    unexcluded[ftype] += flag["Weight"]
                    flagged_features[ftype].add(feature["Id"])
                    flags.append(flag)
                else:
                    excluded_length[ftype][flag["Exclusion"]] += flag["Weight"]
        own_mask[own_segments] = False
    sites = group_sites(flags, site_radius_over_r * radius)
    histogram = defaultdict(float)
    for f in flags:
        histogram[f"{round(f['Distance'] / radius * 4) / 4:.2f}R"] += f["Weight"]
    by_class = {
        cls: {"Features": counts[cls], "SampledLength": totals[cls], "FacingLength": facing.get(cls, 0.0),
              "FacingFraction": facing.get(cls, 0.0) / totals[cls] if totals[cls] else 0.0,
              "AtExactly2RLength": at_2r.get(cls, 0.0),
              "ExcludedLength": {k: excluded_length[cls][k] for k in EXCLUSION_CLASSES if excluded_length[cls].get(k)},
              "UnexcludedFacingLength": unexcluded.get(cls, 0.0),
              "FlaggedFeatures": len(flagged_features.get(cls, ()))}
        for cls in sorted(totals)
    }
    isolated_total = sum(totals[c] for c in ISOLATED_CLASSES)
    isolated_facing = sum(facing.get(c, 0.0) for c in ISOLATED_CLASSES)
    pair_total = sum(totals[c] for c in PAIR_CLASSES)
    pair_facing = sum(facing.get(c, 0.0) for c in PAIR_CLASSES)
    def excluded_by(classes):
        out = defaultdict(float)
        for c in classes:
            for k, v in excluded_length[c].items():
                out[k] += v
        return {k: out[k] for k in EXCLUSION_CLASSES if out.get(k)}
    isolated_unexcluded = sum(unexcluded.get(c, 0.0) for c in ISOLATED_CLASSES)
    pair_unexcluded = sum(unexcluded.get(c, 0.0) for c in PAIR_CLASSES)
    return {
        "Manifest": manifest_path if isinstance(manifest_path, str) else None, "MatchingRadius": radius, "Spacing": spacing, "SiteRadius": site_radius_over_r * radius,
        "Segments": len(p0), "Features": len(features),
        "ByClass": by_class,
        "Isolated": {"SampledLength": isolated_total, "FacingLength": isolated_facing, "FacingFraction": isolated_facing / isolated_total if isolated_total else 0.0,
                     "AtExactly2RLength": sum(at_2r.get(c, 0.0) for c in ISOLATED_CLASSES),
                     "ExcludedLength": excluded_by(ISOLATED_CLASSES), "UnexcludedFacingLength": isolated_unexcluded},
        "Pairs": {"SampledLength": pair_total, "ThirdEdgeLength": pair_facing, "ThirdEdgeFraction": pair_facing / pair_total if pair_total else 0.0,
                  "AtExactly2RLength": sum(at_2r.get(c, 0.0) for c in PAIR_CLASSES),
                  "ExcludedLength": excluded_by(PAIR_CLASSES), "UnexcludedFacingLength": pair_unexcluded},
        "Gates": {"IsolatedFacing": isolated_unexcluded == 0.0, "PairThirdEdge": pair_unexcluded == 0.0},
        "KnifeEdgeRule": f"strict less than 2R on the grid of {LENGTH_QUANTUM_OVER_R:g} R (the classifier's rule); AtExactly2RLength is the facing length at exactly 2R, not flagged",
        "Exclusions": {"AtExactly2R": "facing at exactly 2R on the decision grid (no interaction under the strict rule)",
                       "ClusterNeighbour": "the facing portion belongs to a SpatialEdgeCluster (claim radius R around the event cores, decision 82(1))",
                       "VertexNeighbour": "the facing portion belongs to a vertex feature (corner window, rounded-corner arc, endpoint / junction window)",
                       "ThroughVertex": f"sample and facing point within {THROUGH_VERTEX_ZONE_OVER_R:g} R of one vertex feature (arms meeting through it, design (b) 1)",
                       "SelfNeighbourhood": f"facing point on the sample's own chain less than {SELF_PAIR_NEIGHBOURHOOD_OVER_R:.6g} R of arc length away (the local neighbourhood of a bend, SelfPairNeighbourhoodOverR)"},
        "DistanceHistogram": dict(sorted(histogram.items())),
        "Sites": sites,
        "FlaggedSamples": len(flags),
        "Flags": flags,
    }


def facing_gates(manifest, spacing=0.5):
    """The two hard gates of the audit (decision 82(2)): no isolated / curved edge portion faces
    another edge of its plane within 2R and no pair / stack side faces a third edge within 2R,
    apart from the recorded exclusions (reported with their lengths)."""
    result = facing_check(manifest, spacing=spacing)
    result.pop("Flags", None)
    gates = []
    for name, key, block in (("A8-isolated-facing", "IsolatedFacing", result["Isolated"]), ("A8-pair-third-edge", "PairThirdEdge", result["Pairs"])):
        gates.append({"Gate": name, "Status": "PASS" if result["Gates"][key] else "FAIL",
                      "Detail": {"UnexcludedFacingLength": block["UnexcludedFacingLength"], "FacingLength": block.get("FacingLength", block.get("ThirdEdgeLength")),
                                 "SampledLength": block["SampledLength"], "ExcludedLength": block["ExcludedLength"], "AtExactly2RLength": block["AtExactly2RLength"],
                                 "Sites": [{k: s[k] for k in ("Length", "Classes", "Features", "MinDistance", "Box")} for s in result["Sites"][:5]],
                                 "Basis": "facing_check.py: strict < 2R (3D) on the decision grid; exclusions AtExactly2R / ClusterNeighbour / VertexNeighbour / ThroughVertex recorded"}})
    return gates, result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("manifest")
    parser.add_argument("--output", required=True)
    parser.add_argument("--spacing", type=float, default=0.5, help="sample spacing along the portions (manifest units)")
    parser.add_argument("--site-radius", type=float, default=25.0, help="flagged samples closer than this (in R) form one site")
    parser.add_argument("--flags", help="also write every flagged sample (JSON lines) here")
    args = parser.parse_args(argv)
    result = facing_check(args.manifest, spacing=args.spacing, site_radius_over_r=args.site_radius)
    flags = result.pop("Flags")
    if args.flags:
        with open(args.flags, "w") as target:
            for flag in flags:
                target.write(json.dumps(flag) + "\n")
    with open(args.output, "w") as target:
        json.dump(result, target, indent=1)
    for cls, row in result["ByClass"].items():
        print(f"{cls}: {row['Features']} features, {row['SampledLength']:.1f} um; {row['FacingLength']:.1f} um ({100 * row['FacingFraction']:.1f}%) "
              f"{'face another edge' if cls in ISOLATED_CLASSES else 'face a third edge'} within 2R, of which {row['UnexcludedFacingLength']:.1f} um outside the recorded exclusions "
              f"({row['FlaggedFeatures']} features); excluded {({k: round(v, 1) for k, v in row['ExcludedLength'].items()})}; "
              f"{row['AtExactly2RLength']:.1f} um at exactly 2R (not flagged)")
    print(f"isolated: {result['Isolated']['UnexcludedFacingLength']:.1f} um unexcluded of {result['Isolated']['FacingLength']:.1f} facing / {result['Isolated']['SampledLength']:.1f} um; "
          f"pairs / stacks: {result['Pairs']['UnexcludedFacingLength']:.1f} um unexcluded of {result['Pairs']['ThirdEdgeLength']:.1f} with a third edge / {result['Pairs']['SampledLength']:.1f} um; "
          f"gates isolated {'PASS' if result['Gates']['IsolatedFacing'] else 'FAIL'}, pairs {'PASS' if result['Gates']['PairThirdEdge'] else 'FAIL'}; "
          f"{len(result['Sites'])} unexcluded sites; top: {[round(s['Length'], 1) for s in result['Sites'][:10]]}")


if __name__ == "__main__":
    main()
