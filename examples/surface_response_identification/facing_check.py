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
* Pair portions (`SameConductorGap`, `DifferentConductorGap`, `SameConductorStrip`, the curved
  pairs): the length whose sample points face a THIRD edge within 2R — a segment that belongs to
  neither side of the pair. A pair with a third edge in reach is a multi-edge cross-section
  (ground - gap - trace - gap - ground narrower than 2R).

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
                "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip")


class SegmentGrid:
    """Uniform grid over 2D segments: candidates for a point are the segments crossing the 3 x 3
    cells around it, so a query radius up to the cell size is exact."""

    def __init__(self, p0, p1, cell):
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
        base = np.floor((points - self.origin) / self.cell).astype(np.int64)
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
    with open(manifest_path) as source:
        manifest = json.load(source)
    identification = manifest["Identification"]
    segments = identification["Segments"]
    p0 = np.array([s["Key"][0][:2] for s in segments], dtype=float)
    p1 = np.array([s["Key"][1][:2] for s in segments], dtype=float)
    length = np.array([s["Length"] for s in segments], dtype=float)
    excluded = np.array(["Exclusion" in s for s in segments])
    radius = float(identification.get("MatchingRadius") or manifest.get("MatchingRadius") or 2.0)
    return identification, p0, p1, length, excluded, radius


def facing_samples(grid, p0, p1, excluded, points, tangents, owners, own_mask, radius):
    """For every sample point: the closest facing (across) non-excluded segment within 2R that is
    neither the sample's own segment nor in own_mask (the feature's own sides); -1 when none."""
    sample_indices, segment_indices = grid.candidates(points)
    a, b = p0[segment_indices], p1[segment_indices]
    ab = b - a
    p = points[sample_indices]
    t = np.clip(((p - a) * ab).sum(1) / np.maximum((ab * ab).sum(1), 1e-30), 0.0, 1.0)
    d = a + t[:, None] * ab - p
    r = np.sqrt((d * d).sum(1))
    cosine = np.abs((d * tangents[sample_indices]).sum(1)) / np.maximum(r, 1e-30)
    ok = ((r < 2 * radius) & (r > 1e-9) & (cosine < 0.5) & ~excluded[segment_indices]
          & (segment_indices != owners[sample_indices]) & ~own_mask[segment_indices])
    best = np.full(len(points), -1, dtype=np.int64)
    best_r = np.full(len(points), np.inf)
    if ok.any():
        si, sj, ri = sample_indices[ok], segment_indices[ok], r[ok]
        order = np.lexsort((ri, si))
        si, sj, ri = si[order], sj[order], ri[order]
        first = np.concatenate([[True], si[1:] != si[:-1]])
        best[si[first]] = sj[first]
        best_r[si[first]] = ri[first]
    return best, best_r


def portion_samples(p0, p1, length, portions, spacing):
    """Sample points, tangents, weights (length per sample) and segment ids of feature portions."""
    points, tangents, weights, owners = [], [], [], []
    for seg, s0, s1 in portions:
        total = length[seg]
        if total <= 0 or s1 <= s0:
            continue
        tangent = (p1[seg] - p0[seg]) / total
        count = max(2, int(math.ceil((s1 - s0) / spacing)) + 1)
        s = np.linspace(s0, s1, count)
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
    xy = np.array([f["Point"] for f in flags])
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
            "MinDistancePoint": closest["Point"],
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


def facing_check(manifest_path, spacing=0.5, site_radius_over_r=25.0, batch=200000):
    identification, p0, p1, length, excluded, radius = load(manifest_path)
    features = identification["Features"]
    grid = SegmentGrid(p0, p1, cell=2.0 * radius)
    own_mask = np.zeros(len(p0), dtype=bool)
    totals = defaultdict(float)
    facing = defaultdict(float)
    counts = defaultdict(int)
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
        own_segments = [seg for seg, _, _ in feature["Portions"]] if ftype in PAIR_CLASSES else []
        own_mask[own_segments] = True
        for first in range(0, len(points), batch):
            sl = slice(first, first + batch)
            best, best_r = facing_samples(grid, p0, p1, excluded, points[sl], tangents[sl], owners[sl], own_mask, radius)
            hit = best >= 0
            if not hit.any():
                continue
            facing[ftype] += float(weights[sl][hit].sum())
            flagged_features[ftype].add(feature["Id"])
            for k in np.nonzero(hit)[0]:
                flags.append({"Class": ftype, "Feature": feature["Id"], "Segment": int(owners[sl][k]), "FacingSegment": int(best[k]),
                              "Distance": float(best_r[k]), "Weight": float(weights[sl][k]), "Point": [float(v) for v in points[sl][k]]})
        own_mask[own_segments] = False
    sites = group_sites(flags, site_radius_over_r * radius)
    histogram = defaultdict(float)
    for f in flags:
        histogram[f"{round(f['Distance'] / radius * 4) / 4:.2f}R"] += f["Weight"]
    by_class = {
        cls: {"Features": counts[cls], "SampledLength": totals[cls], "FacingLength": facing.get(cls, 0.0),
              "FacingFraction": facing.get(cls, 0.0) / totals[cls] if totals[cls] else 0.0,
              "FlaggedFeatures": len(flagged_features.get(cls, ()))}
        for cls in sorted(totals)
    }
    isolated_total = sum(totals[c] for c in ISOLATED_CLASSES)
    isolated_facing = sum(facing.get(c, 0.0) for c in ISOLATED_CLASSES)
    pair_total = sum(totals[c] for c in PAIR_CLASSES)
    pair_facing = sum(facing.get(c, 0.0) for c in PAIR_CLASSES)
    return {
        "Manifest": manifest_path, "MatchingRadius": radius, "Spacing": spacing, "SiteRadius": site_radius_over_r * radius,
        "Segments": len(p0), "Features": len(features),
        "ByClass": by_class,
        "Isolated": {"SampledLength": isolated_total, "FacingLength": isolated_facing, "FacingFraction": isolated_facing / isolated_total if isolated_total else 0.0},
        "Pairs": {"SampledLength": pair_total, "ThirdEdgeLength": pair_facing, "ThirdEdgeFraction": pair_facing / pair_total if pair_total else 0.0},
        "DistanceHistogram": dict(sorted(histogram.items())),
        "Sites": sites,
        "FlaggedSamples": len(flags),
        "Flags": flags,
    }


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
              f"{'face another edge' if cls in ISOLATED_CLASSES else 'face a third edge'} within 2R ({row['FlaggedFeatures']} features)")
    print(f"isolated: {result['Isolated']['FacingLength']:.1f} / {result['Isolated']['SampledLength']:.1f} um facing; "
          f"pairs: {result['Pairs']['ThirdEdgeLength']:.1f} / {result['Pairs']['SampledLength']:.1f} um with a third edge; "
          f"{len(result['Sites'])} sites; top: {[round(s['Length'], 1) for s in result['Sites'][:10]]}")


if __name__ == "__main__":
    main()
