#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Placement audit (gate A10, geometry only): every matched cluster / corner / stack / pair
model, mapped through the patch frame the dry run actually used, lands on the feature's
claimed perimeter portions within the signature parameter tolerance (1e-3 R).

    python3 -m surface_response_identification.placement_check \\
        --manifest postpro/surface-response-requirements.json [--patches ...csv] \\
        [--library process-library.json] [--output out.json]

The library's geometry is the model's own description (in library length units, scaled by
R_mesh / R_library):

* SpatialEdgeCluster: the model ``Edges`` (Point, GapDirection, ProcessNormal, Interval;
  tangent = gap x normal), or the chorded Signature portions of a Signature-only model
  (signature_library.cluster_plan_view_edges). Both directions are checked: every model
  edge endpoint lies on the feature's claimed geometry (straight portions as segments, arc
  portions as their fitted circle within the claimed angular range), and every claimed
  straight-portion endpoint / arc end lies on the mapped model polyline.
* ConvexCorner / ConcaveCorner: arms from the patch origin along +u and at the model
  ``Angle`` counterclockwise in (u, v) (the corner coupon's frame, generate_corner_response
  corner_frame); a rounded corner's claimed arc chords lie on the fillet circle of the
  model ``CornerRadius`` (centre CornerRadius x bisector / sin(Angle / 2)).
* ParallelEdgeCluster (stack): per longitudinal patch, origin + Offset_k u lies on the
  claimed portions of side k (sides in increasing lateral offset, reversed for chirality
  -1, as the patch construction orders them).
* SameConductorGap / DifferentConductorGap / SameConductorStrip: origin -+ Separation / 2 u
  lies on side 0 / 1 (swapped for chirality -1).

A model whose placement cannot be evaluated from its record (a curved pair, a topology the
gate does not cover) is listed, never silently passed. This gate is the one a misplaced
coupon fails while every coverage gate passes (the 2394fdb0c failure class: a version-2
cluster model canonicalised to another frame; a chorded arc model re-canonicalised from its
chords).
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import signature_library  # noqa: E402

# Identification.Conventions SignatureParameterToleranceOverR: two signatures agree when their
# lengths agree within this fraction of R; a placed model must agree with its feature at least
# as well.
SIGNATURE_TOLERANCE_OVER_R = 1.0e-3

CLUSTER_TYPES = ("SpatialEdgeCluster",)
CORNER_TYPES = ("ConvexCorner", "ConcaveCorner")
STACK_TYPES = ("ParallelEdgeCluster",)
PAIR_TYPES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip")


def gate(name, passed, detail, evaluable=True):
    status = "PASS" if passed else ("FAIL" if evaluable else "NOT-EVALUABLE")
    return {"Gate": name, "Status": status, "Detail": detail}


def load_patches(path):
    import csv
    rows = []
    with open(path, newline="") as source:
        for row in csv.DictReader(source):
            entry = dict(row)
            for key in ("Patch", "Feature", "ModelIndex", "Segment"):
                entry[key] = int(row[key])
            for key in ("Weight", "ModelWeight", "QuadratureWeight", "SideFactor", "CouponDepth", "S0", "S1"):
                entry[key] = float(row[key])
            entry["Origin"] = np.asarray([float(row[k]) for k in ("OriginX", "OriginY", "OriginZ")])
            entry["Axes"] = np.asarray([[float(row[f"Axis{a}{d}"]) for d in "XYZ"] for a in "UVW"])
            rows.append(entry)
    return rows


def _segment_distance(point, a, b):
    d = b - a
    length2 = float(np.dot(d, d))
    s = float(np.dot(point - a, d)) / length2 if length2 > 0.0 else 0.0
    s = min(max(s, 0.0), 1.0)
    return float(np.linalg.norm(point - (a + s * d)))


class ClaimedGeometry:
    """The claimed perimeter portions of one feature in mesh coordinates: straight
    sub-segments (per mesh segment) and, per fitted arc, the circle with the angular range
    of the claimed chords (the plane of the arc is the process plane of the feature)."""

    def __init__(self, feature, segments, arcs, normal):
        self.normal = np.asarray(normal, dtype=float)
        self.normal /= np.linalg.norm(self.normal)
        self.straight = []          # (p0, p1, side)
        self.arc_ranges = {}        # arc index -> dict(center, radius, u, v, angles, side)
        self.points = []            # claimed straight endpoints and arc ends (for the reverse check)
        sides = feature.get("Sides") or [0] * len(feature["Portions"])
        for (seg, s0, s1), side in zip(feature["Portions"], sides):
            segment = segments[int(seg)]
            key = np.asarray(segment["Key"], dtype=float)
            length = float(segment["Length"])
            p0 = key[0] + (float(s0) / length) * (key[1] - key[0])
            p1 = key[0] + (float(s1) / length) * (key[1] - key[0])
            if "Arc" not in segment:
                self.straight.append((p0, p1, int(side)))
                self.points.append(p0)
                self.points.append(p1)
                continue
            index = int(segment["Arc"])
            arc = arcs[index]
            entry = self.arc_ranges.get(index)
            if entry is None:
                center = np.asarray(arc["Center"], dtype=float)
                # In-plane basis of the arc: u toward the first claimed point, v = n x u.
                u = p0 - center
                u -= float(np.dot(u, self.normal)) * self.normal
                u /= np.linalg.norm(u)
                v = np.cross(self.normal, u)
                entry = {"center": center, "radius": float(arc["Radius"]), "u": u, "v": v, "angles": [], "side": int(side), "ends": {}}
                self.arc_ranges[index] = entry
            for q in (p0, p1):
                angle = self._angle(entry, q)
                entry["angles"].append(angle)
                entry["ends"][round(angle, 9)] = q
        for entry in self.arc_ranges.values():
            entry["min"], entry["max"] = min(entry["angles"]), max(entry["angles"])
            # The arc's two ends (extreme claimed angles) belong to the reverse check.
            self.points.append(entry["ends"][round(entry["min"], 9)])
            self.points.append(entry["ends"][round(entry["max"], 9)])

    @staticmethod
    def _angle(entry, q):
        r = q - entry["center"]
        return math.atan2(float(np.dot(r, entry["v"])), float(np.dot(r, entry["u"])))

    def _arc_distance(self, entry, point):
        r = point - entry["center"]
        out_of_plane = float(np.dot(r, self.normal))
        in_plane = r - out_of_plane * self.normal
        radial = abs(float(np.linalg.norm(in_plane)) - entry["radius"])
        angle = self._angle(entry, point)
        excess = max(0.0, entry["min"] - angle, angle - entry["max"])
        return math.hypot(radial, excess * entry["radius"], out_of_plane)

    def distance(self, point, side=None):
        """Distance from a point to the claimed geometry (optionally of one side)."""
        best = math.inf
        for p0, p1, s in self.straight:
            if side is None or s == side:
                best = min(best, _segment_distance(point, p0, p1))
        for entry in self.arc_ranges.values():
            if side is None or entry["side"] == side:
                best = min(best, self._arc_distance(entry, point))
        return best

    def sides(self):
        labels = {s for _, _, s in self.straight} | {e["side"] for e in self.arc_ranges.values()}
        return sorted(labels)


def model_cluster_edges(model, radius_library):
    """The model's Edges as (P0, P1) pairs in library units (tangent = gap x normal), from
    the stored Edges or the chorded Signature portions of a Signature-only model."""
    edges = model.get("Edges")
    if edges:
        pairs = []
        for edge in edges:
            point = np.asarray(edge["Point"], dtype=float)
            gap = np.asarray(edge["GapDirection"], dtype=float)
            normal = np.asarray(edge["ProcessNormal"], dtype=float)
            tangent = np.cross(gap, normal)
            tangent /= np.linalg.norm(tangent)
            a, b = edge["Interval"]
            pairs.append((point + float(a) * tangent, point + float(b) * tangent))
        return pairs, "Edges"
    signature = model.get("Signature")
    if signature and signature.get("Type") == "SpatialEdgeCluster":
        pairs = []
        for edge in signature_library.cluster_plan_view_edges(signature, radius_library):
            pairs.append((np.asarray([edge["P0"][0], edge["P0"][1], 0.0]), np.asarray([edge["P1"][0], edge["P1"][1], 0.0])))
        return pairs, "Signature"
    return None, None


def placement_gates(identification, patches, library, radius):
    """Gates A10 on the placement of every matched cluster / corner / stack / pair model;
    returns (gates, summary)."""
    tol = SIGNATURE_TOLERANCE_OVER_R * radius
    radius_library = float(library["MatchingRadius"])
    scale = radius / radius_library
    models = {m["Name"]: m for m in library["Models"]}
    features = {f["Id"]: f for f in identification["Features"]}
    segments = identification["Segments"]
    arcs = identification.get("Arcs", [])
    rows_by_feature = defaultdict(list)
    for row in patches:
        rows_by_feature[row["Feature"]].append(row)

    results = {"Clusters": [], "Corners": [], "Stacks": [], "Pairs": []}
    not_evaluable = []
    for feature_id, rows in sorted(rows_by_feature.items()):
        feature = features.get(feature_id)
        if feature is None or feature["Match"]["Status"] != "Matched":
            continue
        kind = feature["Type"]
        model_name = rows[0]["Model"]
        model = models.get(model_name)
        if model is None:
            # The dry run wrote the model's name; a library without it cannot be audited.
            not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name, "Reason": "model not in the library"})
            continue
        normal = np.asarray(feature["Frame"]["Axes"][2], dtype=float)
        claimed = ClaimedGeometry(feature, segments, arcs, normal)
        entry = {"Feature": feature_id, "Type": kind, "Model": model_name, "WorstOverR": 0.0, "Checks": 0, "Defects": []}

        def record(point_label, deviation, mapped=None):
            entry["Checks"] += 1
            entry["WorstOverR"] = max(entry["WorstOverR"], deviation / radius)
            if deviation > tol and len(entry["Defects"]) < 8:
                entry["Defects"].append({"Point": point_label, "DeviationOverR": deviation / radius, "Mapped": None if mapped is None else [float(x) for x in mapped]})

        if kind in CLUSTER_TYPES:
            pairs, source = model_cluster_edges(model, radius_library)
            if pairs is None:
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name, "Reason": "no Edges and no Signature"})
                continue
            if len(rows) != 1:
                entry["Defects"].append({"Point": "patch count", "DeviationOverR": math.inf, "Mapped": None, "Count": len(rows)})
            row = rows[0]
            origin, axes = row["Origin"], row["Axes"]
            mapped_pairs = []
            for k, (a, b) in enumerate(pairs):
                ma = origin + (scale * a) @ axes
                mb = origin + (scale * b) @ axes
                mapped_pairs.append((ma, mb))
                for label, q in ((f"edge {k} begin", ma), (f"edge {k} end", mb)):
                    record(label, claimed.distance(q), q)
            # Reverse: every claimed straight endpoint / arc end on the mapped model polyline.
            for j, q in enumerate(claimed.points):
                record(f"claimed point {j}", min(_segment_distance(q, ma, mb) for ma, mb in mapped_pairs), q)
            entry["EdgeSource"] = source
            entry["ModelEdges"] = len(pairs)
            results["Clusters"].append(entry)
        elif kind in CORNER_TYPES:
            angle = math.radians(float(model.get("Angle", feature["Signature"].get("AngleDegrees", 0.0))))
            corner_radius = float(model.get("CornerRadius", 0.0)) * scale
            if len(rows) != 1:
                entry["Defects"].append({"Point": "patch count", "DeviationOverR": math.inf, "Mapped": None, "Count": len(rows)})
            row = rows[0]
            origin, axes = row["Origin"], row["Axes"]
            rays = [np.asarray([1.0, 0.0]), np.asarray([math.cos(angle), math.sin(angle)])]
            bisector = np.asarray([math.cos(0.5 * angle), math.sin(0.5 * angle)])
            center = corner_radius * bisector / math.sin(0.5 * angle) if corner_radius > 0.0 else None
            tangent_length = corner_radius / math.tan(0.5 * angle) if corner_radius > 0.0 else 0.0

            def local(q):
                return axes @ (q - origin)

            def ray_distance(p):
                best = math.inf
                for ray in rays:
                    t = float(np.dot(p[:2], ray))
                    perpendicular = float(np.linalg.norm(p[:2] - t * ray))
                    along = max(0.0, tangent_length - t)
                    best = min(best, math.hypot(perpendicular, along, p[2]))
                return best

            for p0, p1, _ in claimed.straight:
                for label, q in (("straight begin", p0), ("straight end", p1)):
                    record(label, ray_distance(local(q)), q)
            for index, arc in claimed.arc_ranges.items():
                if center is None:
                    entry["Defects"].append({"Point": f"arc {index}", "DeviationOverR": math.inf, "Mapped": None, "Reason": "sharp model on a rounded feature"})
                    entry["Checks"] += 1
                    continue
                for q in (arc["ends"][round(arc["min"], 9)], arc["ends"][round(arc["max"], 9)]):
                    p = local(q)
                    record(f"arc {index} end", math.hypot(abs(float(np.linalg.norm(p[:2] - center)) - corner_radius), p[2]), q)
            if corner_radius == 0.0:
                # The sharp corner vertex is the patch origin.
                record("corner vertex", min(float(np.linalg.norm(q - origin)) for q in claimed.points), origin)
            entry["Angle"] = math.degrees(angle)
            entry["CornerRadius"] = corner_radius
            results["Corners"].append(entry)
        elif kind in STACK_TYPES:
            edges = model.get("Edges")
            if edges and "Offset" in edges[0]:
                offsets = [float(e["Offset"]) * scale for e in edges]
            elif model.get("Signature") and "Edges" in model["Signature"]:
                offsets = [float(e["OffsetOverR"]) * radius for e in model["Signature"]["Edges"]]
            else:
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name, "Reason": "no Edges and no Signature"})
                continue
            sides = claimed.sides()
            chirality = int(feature.get("Chirality", 1))
            ordered = sides if chirality >= 0 else list(reversed(sides))
            if len(ordered) != len(offsets):
                entry["Defects"].append({"Point": "side count", "DeviationOverR": math.inf, "Mapped": None, "Sides": len(ordered), "ModelEdges": len(offsets)})
                entry["Checks"] += 1
            else:
                for row in rows:
                    origin, axes = row["Origin"], row["Axes"]
                    for k, offset in enumerate(offsets):
                        q = origin + offset * axes[0]
                        record(f"patch {row['Patch']} edge {k}", claimed.distance(q, side=ordered[k]), q)
            entry["Patches"] = len(rows)
            results["Stacks"].append(entry)
        elif kind in PAIR_TYPES:
            if "Separation" in model:
                separation = float(model["Separation"]) * scale
            elif model.get("Signature") and "SeparationOverR" in model["Signature"]:
                separation = float(model["Signature"]["SeparationOverR"]) * radius
            else:
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name, "Reason": "no Separation and no Signature"})
                continue
            sides = claimed.sides()
            chirality = int(feature.get("Chirality", 1))
            ordered = sides if chirality >= 0 else list(reversed(sides))
            if len(ordered) != 2:
                entry["Defects"].append({"Point": "side count", "DeviationOverR": math.inf, "Mapped": None, "Sides": len(ordered)})
                entry["Checks"] += 1
            else:
                for row in rows:
                    origin, axes = row["Origin"], row["Axes"]
                    for k, sign in enumerate((-1.0, 1.0)):
                        q = origin + sign * 0.5 * separation * axes[0]
                        record(f"patch {row['Patch']} side {k}", claimed.distance(q, side=ordered[k]), q)
            entry["Patches"] = len(rows)
            results["Pairs"].append(entry)
        else:
            continue  # isolated / curved edges: the patch is the quadrature point itself

    gates = []
    summary = {"ToleranceOverR": SIGNATURE_TOLERANCE_OVER_R, "Scale": scale, "NotEvaluable": not_evaluable}
    for name, key in (("A10-placement-clusters", "Clusters"), ("A10-placement-corners", "Corners"), ("A10-placement-stacks", "Stacks"), ("A10-placement-pairs", "Pairs")):
        entries = results[key]
        defects = [e for e in entries if e["Defects"]]
        worst = max((e["WorstOverR"] for e in entries), default=None)
        detail = {
            "Features": len(entries),
            "Checks": sum(e["Checks"] for e in entries),
            "WorstDeviationOverR": worst,
            "ToleranceOverR": SIGNATURE_TOLERANCE_OVER_R,
            "Models": sorted({e["Model"] for e in entries}),
            "FeaturesWithDefects": len(defects),
            "Examples": [{k: v for k, v in e.items() if k != "Checks"} for e in defects[:6]],
            "Basis": "model geometry (library units x R_mesh / R_library) mapped through the dry-run patch frame lies on the feature's claimed portions (straight segments; arcs on their fitted circle within the claimed range) within the signature parameter tolerance",
        }
        gates.append(gate(name, not defects, detail, evaluable=bool(entries)))
        summary[key] = {"Features": len(entries), "WorstDeviationOverR": worst, "FeaturesWithDefects": len(defects)}
    if not_evaluable:
        gates.append(gate("A10-placement-evaluable", False, {"NotEvaluable": not_evaluable[:20], "Count": len(not_evaluable)}))
    return gates, summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", required=True, help="surface-response-requirements.json (version 2)")
    parser.add_argument("--patches", help="surface-response-patches.csv (default: next to the manifest)")
    parser.add_argument("--library", help="process-library.json (default: the manifest's Library.Path)")
    parser.add_argument("--output", help="write the gates and summary as JSON")
    args = parser.parse_args(argv)
    with open(args.manifest) as source:
        manifest = json.load(source)
    identification = manifest.get("Identification")
    if not identification:
        raise SystemExit("the manifest carries no version-2 Identification")
    patches_path = args.patches or os.path.join(os.path.dirname(os.path.abspath(args.manifest)), "surface-response-patches.csv")
    library_path = args.library or (manifest.get("Library") or {}).get("Path")
    if not library_path:
        raise SystemExit("no library: pass --library")
    with open(library_path) as source:
        library = json.load(source)
    gates, summary = placement_gates(identification, load_patches(patches_path), library, float(identification["MatchingRadius"]))
    result = {"Manifest": os.path.abspath(args.manifest), "Patches": os.path.abspath(patches_path), "Library": os.path.abspath(library_path), "Gates": gates, "Summary": summary,
              "Passed": all(g["Status"] == "PASS" for g in gates if g["Status"] != "NOT-EVALUABLE")}
    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        with open(args.output, "w") as out:
            json.dump(result, out, indent=1, default=str)
    for g in gates:
        detail = g["Detail"]
        print(f"{g['Gate']}: {g['Status']} (features {detail.get('Features', detail.get('Count'))}, worst {detail.get('WorstDeviationOverR')} R)")
    return 0 if result["Passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
