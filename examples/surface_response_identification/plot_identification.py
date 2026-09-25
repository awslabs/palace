#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Plan-view plots of a version-2 identification manifest for visual review (chip scale).

Every feature's perimeter portions are drawn colored by class; excluded segments dashed gray;
vertex features as markers; spatial / parallel clusters boxed and labeled. The metal fill comes
from `mesh_census.py --extract` (triangles of the metal attributes clipped to the plot windows,
`--metal tris.npz`); without it the perimeter outline alone is drawn. Metal planes are told apart
by the z coordinate of the segments (one overview per plane).

Windows: `--zoom X0 X1 Y0 Y1 [--zoom ...]` explicit; `--auto` adds (a) the largest clusters,
(b) one example of every class, (c) the top facing-check sites (`--facing facing.json`, from
`facing_check.py`); `--regions regions.json` adds named windows ({"name": [x0, x1, y0, y1], ...}).
Every figure is listed in `<prefix>-figures.json` (file, kind, window, what it shows).

    python3 -m surface_response_identification.plot_identification MANIFEST OUT_PREFIX \\
        [--metal tris.npz] [--facing facing.json] [--auto] [--regions regions.json] [--zoom ...]
"""
import argparse
import json
import os
import textwrap
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402

COLORS = {
    "IsolatedEdge": "#222222",
    "CurvedEdge": "#8e44ad",
    "SameConductorGap": "#1f77b4",
    "DifferentConductorGap": "#17becf",
    "SameConductorStrip": "#2ca02c",
    "CurvedSameConductorStrip": "#98df8a",
    "CurvedSameConductorGap": "#aec7e8",
    "CurvedDifferentConductorGap": "#9edae5",
    "ParallelEdgeCluster": "#ff7f0e",
    "SpatialEdgeCluster": "#d62728",
    "ConvexCorner": "#e377c2",
    "ConcaveCorner": "#bcbd22",
    "RoundedConvexCorner": "#e377c2",
    "RoundedConcaveCorner": "#bcbd22",
    "Endpoint": "#7f7f7f",
    "Junction": "#8c564b",
}
WIDTH = {"SpatialEdgeCluster": 2.5, "ParallelEdgeCluster": 2.0}
VERTEX_CLASSES = ("ConvexCorner", "ConcaveCorner", "RoundedConvexCorner", "RoundedConcaveCorner", "Endpoint", "Junction")
CLUSTER_CLASSES = ("SpatialEdgeCluster", "ParallelEdgeCluster")


def load(manifest_path):
    with open(manifest_path) as source:
        manifest = json.load(source)
    ident = manifest["Identification"]
    return ident["Segments"], ident["Features"], ident


def portion_xy(segments, seg, s0, s1):
    (x0, y0, _), (x1, y1, _) = segments[seg]["Key"]
    length = segments[seg]["Length"]
    t0, t1 = (s0 / length, s1 / length) if length > 0 else (0.0, 1.0)
    return ((x0 + t0 * (x1 - x0), y0 + t0 * (y1 - y0)), (x0 + t1 * (x1 - x0), y0 + t1 * (y1 - y0)))


def feature_box(segments, feature):
    pts = [p for seg, s0, s1 in feature.get("Portions", []) for p in portion_xy(segments, seg, s0, s1)]
    if not pts and "Frame" in feature:
        x, y = feature["Frame"]["Origin"][:2]
        pts = [(x, y)]
    if not pts:
        return None
    xs, ys = zip(*pts)
    return [min(xs), max(xs), min(ys), max(ys)]


def feature_plane(segments, feature):
    for seg, _, _ in feature.get("Portions", []):
        return round(segments[seg]["Key"][0][2], 6)
    if "Frame" in feature:
        return round(feature["Frame"]["Origin"][2], 6)
    return None


def segment_planes(segments):
    return sorted({round(s["Key"][0][2], 6) for s in segments})


def load_metal(path):
    if not path:
        return np.zeros((0, 3, 2)), np.zeros(0, dtype=int)
    data = np.load(path)
    return data["xy"], data["attribute"]


def draw(ax, segments, features, tris, tri_tags, view=None, plane=None, label_clusters=True, lw_scale=1.0,
         show_excluded=True, metal_palette=None, facing_sites=None):
    def in_view(box):
        return not view or not (box[1] < view[0] or box[0] > view[1] or box[3] < view[2] or box[2] > view[3])

    if len(tris):
        sel = np.ones(len(tris), dtype=bool)
        if view:
            sel = (tris[:, :, 0].max(1) >= view[0]) & (tris[:, :, 0].min(1) <= view[1]) & (tris[:, :, 1].max(1) >= view[2]) & (tris[:, :, 1].min(1) <= view[3])
        palette = metal_palette or {}
        face = [palette.get(int(t), "#d9d9d9") for t in tri_tags[sel]]
        ax.add_collection(PolyCollection(tris[sel], facecolors=face, edgecolors="none", zorder=0))
    by_type = defaultdict(list)
    for f in features:
        if plane is not None and feature_plane(segments, f) != plane:
            continue
        for seg, s0, s1 in f.get("Portions", []):
            line = portion_xy(segments, seg, s0, s1)
            if in_view([min(line[0][0], line[1][0]), max(line[0][0], line[1][0]), min(line[0][1], line[1][1]), max(line[0][1], line[1][1])]):
                by_type[f["Type"]].append(line)
    for ftype, lines in sorted(by_type.items()):
        ax.add_collection(LineCollection(lines, colors=COLORS.get(ftype, "#999999"),
                                         linewidths=WIDTH.get(ftype, 1.2) * lw_scale, zorder=2, label=f"{ftype} ({len(lines)})"))
    if show_excluded:
        excl = defaultdict(list)
        for s in segments:
            if "Exclusion" in s and (plane is None or round(s["Key"][0][2], 6) == plane):
                line = tuple(tuple(p[:2]) for p in s["Key"])
                if in_view([min(line[0][0], line[1][0]), max(line[0][0], line[1][0]), min(line[0][1], line[1][1]), max(line[0][1], line[1][1])]):
                    excl[s["Exclusion"] if isinstance(s["Exclusion"], str) else s["Exclusion"].get("Class", "Excluded")].append(line)
        for cls, lines in sorted(excl.items()):
            ax.add_collection(LineCollection(lines, colors="#999999", linewidths=0.8 * lw_scale, linestyles="dashed", zorder=1,
                                             label=f"Excluded {cls} ({len(lines)})"))
    for f in features:
        if f["Type"] in VERTEX_CLASSES and "Frame" in f and (plane is None or feature_plane(segments, f) == plane):
            x, y = f["Frame"]["Origin"][:2]
            if in_view([x, x, y, y]):
                ax.plot(x, y, marker="o", ms=4 * lw_scale, color=COLORS.get(f["Type"], "k"), zorder=3)
    if label_clusters:
        for f in features:
            if f["Type"] not in CLUSTER_CLASSES or (plane is not None and feature_plane(segments, f) != plane):
                continue
            box = feature_box(segments, f)
            if not box or not in_view(box):
                continue
            ax.add_patch(plt.Rectangle((box[0], box[2]), box[1] - box[0], box[3] - box[2], fill=False, ec="#d62728", lw=0.8, ls=":", zorder=4))
            ax.annotate(f"C{f['Id']}: {len(f['Portions'])} pieces, {f['Length']:.0f} um", (box[1], box[3]), fontsize=7, color="#d62728", zorder=5)
    if facing_sites:
        for k, site in enumerate(facing_sites):
            box = site["Box"]
            if not in_view(box):
                continue
            ax.add_patch(plt.Rectangle((box[0] - 1, box[2] - 1), box[1] - box[0] + 2, box[3] - box[2] + 2, fill=False, ec="#ff00ff", lw=1.0, ls="--", zorder=4))
            ax.annotate(f"F{k + 1}: {site['Length']:.0f} um {'/'.join(site['Classes'])}", (box[0], box[3] + 1), fontsize=7, color="#ff00ff", zorder=5)
    if view:
        ax.set_xlim(view[0], view[1])
        ax.set_ylim(view[2], view[3])
    else:
        ax.autoscale()
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")


def window_around(box, size):
    cx, cy = 0.5 * (box[0] + box[1]), 0.5 * (box[2] + box[3])
    half = max(size, box[1] - box[0], box[3] - box[2]) / 2 * 1.15
    return [cx - half, cx + half, cy - half, cy + half]


def auto_windows(segments, features, facing, zoom_size, cluster_count=5, site_count=10):
    windows = []
    clusters = sorted((f for f in features if f["Type"] in CLUSTER_CLASSES), key=lambda f: -f["Length"])
    for k, f in enumerate(clusters[:cluster_count]):
        box = feature_box(segments, f)
        windows.append({"Kind": "cluster", "Label": f"largest cluster {k + 1}: C{f['Id']} {f['Type']} {len(f['Portions'])} pieces {f['Length']:.1f} um",
                        "Window": window_around(box, zoom_size), "Feature": f["Id"]})
    seen = set()
    for f in sorted(features, key=lambda f: -f.get("Length", 0.0)):
        if f["Type"] in seen or f["Type"] in CLUSTER_CLASSES:
            continue
        box = feature_box(segments, f)
        if box is None:
            continue
        seen.add(f["Type"])
        windows.append({"Kind": "class-example", "Label": f"example of {f['Type']}: feature {f['Id']} {f.get('Length', 0.0):.1f} um",
                        "Window": window_around(box, zoom_size), "Feature": f["Id"]})
    if facing:
        for k, site in enumerate(facing.get("Sites", [])[:site_count]):
            box = site["Box"]
            label = f"facing site {k + 1}: {site['Length']:.1f} um of {'/'.join(site['Classes'])} facing within 2R (min {site['MinDistance']:.2f} um; features {site['Features'][:6]})"
            if max(box[1] - box[0], box[3] - box[2]) > 4 * zoom_size:
                # A site spanning a whole route: zoom on its closest-approach point.
                x, y = site.get("MinDistancePoint", site["Center"])
                windows.append({"Kind": "facing-site", "Label": label + f"; detail at the closest approach ({x:.1f}, {y:.1f}); site box x[{box[0]:.0f},{box[1]:.0f}] y[{box[2]:.0f},{box[3]:.0f}]",
                                "Window": window_around([x, x, y, y], 2 * zoom_size), "Site": k})
            else:
                windows.append({"Kind": "facing-site", "Label": label, "Window": window_around(box, zoom_size), "Site": k})
    return windows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest")
    ap.add_argument("out_prefix")
    ap.add_argument("--metal", help="npz of metal triangles (mesh_census.py --extract)")
    ap.add_argument("--facing", help="facing_check.py JSON: sites are boxed and zoomed")
    ap.add_argument("--regions", help="JSON {name: [x0, x1, y0, y1]} of named windows")
    ap.add_argument("--zoom", type=float, nargs=4, action="append", default=[], metavar=("X0", "X1", "Y0", "Y1"))
    ap.add_argument("--auto", action="store_true", help="largest clusters, one example per class, facing sites")
    ap.add_argument("--zoom-size", type=float, default=60.0, help="minimum side of an automatic zoom window (um)")
    ap.add_argument("--clusters", type=int, default=5)
    ap.add_argument("--sites", type=int, default=10)
    ap.add_argument("--dpi", type=int, default=170)
    ap.add_argument("--title", default="")
    args = ap.parse_args(argv)
    segments, features, ident = load(args.manifest)
    tris, tags = load_metal(args.metal)
    facing = json.load(open(args.facing)) if args.facing else None
    os.makedirs(os.path.dirname(os.path.abspath(args.out_prefix)), exist_ok=True)
    index = []
    planes = segment_planes(segments)
    for plane in planes:
        fig, ax = plt.subplots(figsize=(18, 11))
        plane_tris = tris
        draw(ax, segments, features, plane_tris, tags, plane=plane, label_clusters=True, facing_sites=facing["Sites"][: args.sites] if facing else None)
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc="upper right", fontsize=7)
        ax.set_title(f"{args.title} metal plane z = {plane:g} um: features by class (clusters boxed red, facing sites magenta)")
        fig.tight_layout()
        path = f"{args.out_prefix}-overview-z{plane:g}.png"
        fig.savefig(path, dpi=args.dpi)
        plt.close(fig)
        index.append({"File": os.path.basename(path), "Kind": "overview", "Plane": plane, "Label": f"chip overview of the metal plane at z = {plane:g} um"})
    windows = [{"Kind": "zoom", "Label": f"zoom x[{v[0]:.0f},{v[1]:.0f}] y[{v[2]:.0f},{v[3]:.0f}]", "Window": v} for v in args.zoom]
    if args.auto:
        windows += auto_windows(segments, features, facing, args.zoom_size, args.clusters, args.sites)
    if args.regions:
        with open(args.regions) as source:
            for name, window in json.load(source).items():
                windows.append({"Kind": "region", "Label": name, "Window": window})
    for k, entry in enumerate(windows):
        view = entry["Window"]
        fig, ax = plt.subplots(figsize=(11, 9))
        draw(ax, segments, features, tris, tags, view=view, lw_scale=1.6, facing_sites=facing["Sites"][: args.sites] if facing else None)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, loc="best", fontsize=7)
        ax.set_title("\n".join(textwrap.wrap(f"{entry['Kind']} {k + 1}: {entry['Label']}", 120)) + f"\nx[{view[0]:.1f},{view[1]:.1f}] y[{view[2]:.1f},{view[3]:.1f}] um", fontsize=8)
        fig.tight_layout()
        path = f"{args.out_prefix}-{entry['Kind']}-{k + 1:02d}.png"
        fig.savefig(path, dpi=args.dpi)
        plt.close(fig)
        index.append({"File": os.path.basename(path), **entry})
    with open(f"{args.out_prefix}-figures.json", "w") as target:
        json.dump(index, target, indent=1)
    print(json.dumps({"Figures": len(index), "Index": f"{args.out_prefix}-figures.json"}))


if __name__ == "__main__":
    main()
