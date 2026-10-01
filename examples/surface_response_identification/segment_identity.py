#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Segment-level interior identity of two identification manifests (validation plan (b)-4, E1; the first
chip-scale A5 test): the window's thin mesh vs the chip, compared over the comparison region.

    python3 -m surface_response_identification.segment_identity \\
        --a chip/surface-response-requirements.json --b window/surface-response-requirements.json \\
        --region X0 X1 Y0 Y1 [--offset dx dy dz] [--spacing 0.25] --output out.json [--markdown out.md]

Both manifests are read as piecewise classifications of their perimeter polylines in absolute coordinates
(``--offset`` is added to B's coordinates): every assigned portion is (Type, signature Hash, topology key,
parameters, feature) and every excluded segment / portion is (Exclusion class). A's pieces inside the region
are sampled every ``--spacing`` R (default 0.25 R); for each sample the nearest piece of B on the same plane
within ``--tolerance`` R (default 1e-3 R, the signature parameter tolerance) is found and the two readings are
compared; then the same from B to A, so lost and gained features are both seen. Every sample is classified:

* ``Identical``: same Type and Hash (or the same exclusion class);
* ``A``: a straight-segment class (IsolatedEdge, sharp corners, straight pairs with ExactParameters,
  SpatialEdgeCluster portions on straight segments) that differs -> a DEFECT (class (A) requires identity);
* ``B``: a mean-separation class (ParallelEdgeCluster; any pair whose ExactParameters is false) with the same
  topology key but other parameters -> listed with the parameter deltas (the E0 prediction decides);
* ``C``: the re-mesh class (CurvedEdge / Curved* pairs, rounded corners, any portion on a fitted-arc segment
  or a cluster arc portion) with the same Type and, for clusters, the same EdgeCount -> allowed (residual 1);
* ``D``: everything else — a changed Type, an assigned portion facing an exclusion (or the reverse), a changed
  cluster EdgeCount, no B perimeter within tolerance (lost / gained perimeter) -> a DEFECT.

A sample on a re-meshed curve (CurvedEdge / Curved* pair, rounded corner, fitted-arc chord, cluster arc portion) whose
chord has no other-side perimeter within ``--tolerance`` is retried within ``--curve-tolerance`` R (default 0.5 R: the
other mesh's chords lie off this chord by up to the sagitta); a same-Type piece found there is class ``C`` (displaced
chord), anything else stays ``D``.

Cut-arc exclusions (VALIDATION-PLAN (h)-9, decision 214 (iii)): ``--exclude E0.json WINDOW`` reads the window's
``E1CutArcExclusions`` written by ``window_query.py`` (the chip arcs cut by a wall / the setback with the reading the
identification's arc rule predicts for the window: Collapse, DroppedJoint, None, and the pieces to exclude). A sample within
``--tolerance`` of an excluded piece (A: the chip's own chords / arms; B: the window's identical chords) is class ``Excluded``
(length reported per direction and per arc, NEVER a defect); its would-be class is kept per arc, and every A or D sample that
falls on a cut arc's segments OUTSIDE the excluded pieces is a ``PredictionMismatch`` (the rule said the arc reads identically
there) -> reported, and a defect like any A / D.

Acceptance: 0 length in classes A and D. The output lists per class the length, the count of samples and up
to 50 examples with coordinates, feature ids, both readings and the parameter deltas.
"""

import argparse
import collections
import json
import math
import sys

CLASS_A_TYPES = ("IsolatedEdge", "ConvexCorner", "ConcaveCorner", "Endpoint", "Junction", "SameConductorStrip",
                 "SameConductorGap", "DifferentConductorGap", "SpatialEdgeCluster")
CLASS_B_TYPES = ("ParallelEdgeCluster",)
CLASS_C_TYPES = ("CurvedEdge", "CurvedSameConductorStrip", "CurvedSameConductorGap", "CurvedDifferentConductorGap",
                 "CurvedParallelEdgeCluster")
PARAMETER_KEYS = ("SeparationOverR", "OffsetOverR", "RadiusOverR", "AngleDegrees", "CornerRadiusOverR", "P", "Arc")


def topology_key(signature):
    """Signature with its continuous parameters removed (the plan's topology key of decision 85(1))."""

    def strip(value):
        if isinstance(value, dict):
            return {k: strip(v) for k, v in value.items() if k not in PARAMETER_KEYS}
        if isinstance(value, list):
            return [strip(v) for v in value]
        return value

    return json.dumps(strip(signature), sort_keys=True, separators=(",", ":"))


def parameters(signature):
    found = {}

    def walk(value, prefix):
        if isinstance(value, dict):
            for k, v in value.items():
                if k in PARAMETER_KEYS and not isinstance(v, (dict, list)):
                    found[prefix + k] = v
                else:
                    walk(v, prefix + k + ".")
        elif isinstance(value, list):
            for i, v in enumerate(value):
                walk(v, prefix + f"[{i}].")

    walk(signature, "")
    return found


class Pieces:
    """Classified perimeter pieces of one manifest with a grid index for nearest-piece queries."""

    def __init__(self, manifest, offset=(0.0, 0.0, 0.0), cell=None):
        ident = manifest["Identification"]
        self.radius = float(ident["MatchingRadius"])
        self.features = {f["Id"]: f for f in ident["Features"]}
        self.pieces = []  # (a, b, plane, reading)
        self.offset = offset
        segments = ident["Segments"]
        exclusions = ident.get("Exclusions", [])
        self._load(segments, exclusions)
        self._index(cell)

    @classmethod
    def from_exclusions(cls, cut_arcs, radius):
        """Index of the cut-arc exclusion pieces (``window_query`` ``E1CutArcExclusions``), reading = the arc record."""
        self = cls.__new__(cls)
        self.radius = radius
        self.features = {}
        self.offset = (0.0, 0.0, 0.0)
        self.pieces = []
        for arc in cut_arcs.get("Arcs", []):
            for a, b in arc.get("Pieces", []):
                self.pieces.append((list(a), list(b), round(a[2], 3), {"Kind": "Excluded", "Arc": arc["Arc"], "Prediction": arc["Prediction"]}))
        self._index(None)
        return self

    def _load(self, segments, exclusions):
        offset = self.offset
        for i, s in enumerate(segments):
            p0 = [s["Key"][0][k] + offset[k] for k in range(3)]
            p1 = [s["Key"][1][k] + offset[k] for k in range(3)]
            length = float(s["Length"])
            plane = round(p0[2], 3)
            arc = s.get("Arc")
            if "Exclusion" in s:
                self.pieces.append((p0, p1, plane, {"Kind": "Exclusion", "Class": s["Exclusion"]["Class"], "Segment": i}))
                continue
            for s0, s1, fid in s.get("Portions", []):
                a, b = self._sub(p0, p1, length, s0, s1)
                f = self.features[fid]
                sig = f.get("Signature", {})
                reading = {"Kind": "Feature", "Type": f["Type"], "Hash": f["Hash"], "Feature": fid,
                           "Topology": topology_key(sig), "Parameters": parameters(sig),
                           "EdgeCount": sig.get("EdgeCount"), "Exact": bool(f.get("ExactParameters", True)),
                           "OnArc": arc is not None or self._portion_on_arc(f, fid, i),
                           "Rounded": float(sig.get("CornerRadiusOverR", 0.0) or 0.0) > 0.0, "Segment": i}
                self.pieces.append((a, b, plane, reading))
            for s0, s1, ex in s.get("ExcludedPortions", []):
                a, b = self._sub(p0, p1, length, s0, s1)
                self.pieces.append((a, b, plane, {"Kind": "Exclusion", "Class": exclusions[ex]["Class"], "Segment": i}))

    def _index(self, cell):
        self.cell = cell or 2.0 * self.radius
        self.grid = collections.defaultdict(list)
        for idx, (a, b, plane, _) in enumerate(self.pieces):
            for cx, cy in self._cells(a, b):
                self.grid[(plane, cx, cy)].append(idx)

    @staticmethod
    def _sub(p0, p1, length, s0, s1):
        t0, t1 = (s0 / length, s1 / length) if length > 0 else (0.0, 1.0)
        return ([p0[k] + t0 * (p1[k] - p0[k]) for k in range(3)], [p0[k] + t1 * (p1[k] - p0[k]) for k in range(3)])

    @staticmethod
    def _portion_on_arc(feature, fid, segment):
        sig = feature.get("Signature", {})
        if feature["Type"] != "SpatialEdgeCluster":
            return False
        return any("Arc" in p for p in sig.get("Portions", []))

    def _cells(self, a, b):
        xs = sorted((a[0], b[0]))
        ys = sorted((a[1], b[1]))
        c0, c1 = int(math.floor(xs[0] / self.cell)), int(math.floor(xs[1] / self.cell))
        d0, d1 = int(math.floor(ys[0] / self.cell)), int(math.floor(ys[1] / self.cell))
        return [(cx, cy) for cx in range(c0, c1 + 1) for cy in range(d0, d1 + 1)]

    def nearest(self, x, y, plane, tolerance):
        cx, cy = int(math.floor(x / self.cell)), int(math.floor(y / self.cell))
        best, best_d = None, tolerance
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for idx in self.grid.get((plane, cx + dx, cy + dy), ()):
                    a, b, _, reading = self.pieces[idx]
                    d = _segment_distance(x, y, a, b)
                    if d <= best_d:
                        best, best_d = reading, d
        return best, best_d


def _segment_distance(px, py, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    l2 = dx * dx + dy * dy
    if l2 <= 0.0:
        return math.hypot(px - a[0], py - a[1])
    t = max(0.0, min(1.0, ((px - a[0]) * dx + (py - a[1]) * dy) / l2))
    return math.hypot(px - a[0] - t * dx, py - a[1] - t * dy)


def classify(ra, rb):
    """Difference class of two readings at one sample (see the module docstring)."""
    if rb is None:
        return "D", "no perimeter within tolerance on the other side"
    if ra["Kind"] == "Exclusion" or rb["Kind"] == "Exclusion":
        if ra["Kind"] == rb["Kind"] and ra["Class"] == rb["Class"]:
            return "Identical", ""
        return "D", "exclusion vs assignment or a different exclusion class"
    if ra["Type"] == rb["Type"] and ra["Hash"] == rb["Hash"]:
        return "Identical", ""
    if ra["Type"] != rb["Type"]:
        return "D", "type changed"
    if ra["Type"] == "SpatialEdgeCluster" and ra["EdgeCount"] != rb["EdgeCount"]:
        return "D", "cluster edge count changed"
    curved = ra["Type"] in CLASS_C_TYPES or ra["Rounded"] or rb["Rounded"] or ra["OnArc"] or rb["OnArc"]
    if curved:
        return "C", "re-mesh class (curve / arc portion)"
    if ra["Type"] in CLASS_B_TYPES or not (ra["Exact"] and rb["Exact"]):
        if ra["Topology"] == rb["Topology"]:
            return "B", "mean-separation class: same topology, other parameters"
        return "D", "mean-separation class with a changed topology key"
    if ra["Type"] in CLASS_A_TYPES:
        return "A", "straight-segment class differs"
    return "D", "unclassified difference"


def parameter_deltas(ra, rb):
    if ra.get("Kind") != "Feature" or rb is None or rb.get("Kind") != "Feature":
        return {}
    pa, pb = ra["Parameters"], rb["Parameters"]
    return {k: [pa.get(k), pb.get(k)] for k in sorted(set(pa) | set(pb)) if pa.get(k) != pb.get(k)}


def is_curved(reading):
    return reading["Kind"] == "Feature" and (reading["Type"] in CLASS_C_TYPES or reading["OnArc"] or reading["Rounded"])


def compare(pieces_a, pieces_b, region, spacing_over_r=0.25, tolerance_over_r=1e-3, direction="A->B", examples=50,
            curve_tolerance_over_r=0.5, excluded=None, arc_of_segment=None, per_arc=None):
    """``excluded``: a ``Pieces.from_exclusions`` index (samples within the tolerance of its pieces are class Excluded);
    ``arc_of_segment``: chip segment index -> cut-arc id, for the prediction check (A samples by their own segment, B samples by
    the nearest A piece's segment); ``per_arc``: dict filled with the excluded length, would-be classes and mismatches per arc."""
    radius = pieces_a.radius
    spacing = spacing_over_r * radius
    tolerance = tolerance_over_r * radius
    curve_tolerance = curve_tolerance_over_r * radius
    x0, x1, y0, y1 = region
    totals = collections.defaultdict(lambda: {"Length": 0.0, "Samples": 0, "Examples": []})
    pairs = collections.Counter()
    arc_of_segment = arc_of_segment or {}
    for a, b, plane, reading in pieces_a.pieces:
        length = math.dist(a[:2], b[:2])
        if length <= 0.0:
            continue
        n = max(1, int(math.ceil(length / spacing)))
        h = length / n
        for k in range(n):
            t = (k + 0.5) / n
            x, y = a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])
            if not (x0 <= x <= x1 and y0 <= y <= y1):
                continue
            other, d = pieces_b.nearest(x, y, plane, tolerance)
            displaced = False
            if other is None and is_curved(reading):
                # re-meshed curve: the other mesh's chords lie off this chord by up to the sagitta
                other, d = pieces_b.nearest(x, y, plane, curve_tolerance)
                displaced = other is not None
            cls, why = classify(reading, other)
            if displaced:
                if cls in ("Identical", "C", "B"):
                    cls, why = "C", f"re-mesh class: chord displaced by {d / radius:.4f} R (same Type)"
                else:
                    why += f" (curve piece; nearest other-side perimeter {d / radius:.4f} R away)"
            chip_segment = reading.get("Segment") if direction == "A->B" else (other or {}).get("Segment")
            arc = arc_of_segment.get(chip_segment)
            hit = excluded.nearest(x, y, plane, tolerance)[0] if excluded is not None else None
            if hit is not None:
                record = per_arc.setdefault(hit["Arc"], {"Prediction": hit["Prediction"], "Excluded": collections.Counter(),
                                                         "WouldBe": collections.Counter(), "Mismatch": collections.Counter()})
                record["Excluded"][direction] += h
                record["WouldBe"][cls] += h
                cls, why = "Excluded", f"cut arc {hit['Arc']} ({hit['Prediction']}): would be {cls}"
            elif arc is not None and per_arc is not None:
                record = per_arc.setdefault(arc, {"Prediction": None, "Excluded": collections.Counter(), "WouldBe": collections.Counter(),
                                                  "Mismatch": collections.Counter()})
                if cls in ("A", "D"):
                    record["Mismatch"][direction] += h
                    why += f" [PredictionMismatch: on cut arc {arc} outside its excluded pieces]"
            entry = totals[cls]
            entry["Length"] += h
            entry["Samples"] += 1
            label = (reading.get("Type") or reading.get("Class"), (other or {}).get("Type") or (other or {}).get("Class"))
            pairs[(cls,) + label] += 1
            if cls != "Identical" and len(entry["Examples"]) < examples:
                entry["Examples"].append({"Point": [round(x, 6), round(y, 6), plane], "Why": why, "Direction": direction,
                                          "A": _brief(reading), "B": _brief(other), "Deltas": parameter_deltas(reading, other)})
    return totals, pairs


def _brief(reading):
    if reading is None:
        return None
    if reading["Kind"] == "Exclusion":
        return {"Exclusion": reading["Class"], "Segment": reading["Segment"]}
    return {"Type": reading["Type"], "Hash": reading["Hash"][:12], "Feature": reading["Feature"], "EdgeCount": reading["EdgeCount"],
            "Exact": reading["Exact"], "OnArc": reading["OnArc"], "Segment": reading["Segment"]}


def run(manifest_a, manifest_b, region, offset=(0.0, 0.0, 0.0), spacing_over_r=0.25, tolerance_over_r=1e-3, examples=50,
        curve_tolerance_over_r=0.5, exclusions=None):
    """``exclusions``: the window's ``E1CutArcExclusions`` record of ``window_query`` (chip = A), or None."""
    pa = Pieces(manifest_a)
    pb = Pieces(manifest_b, offset)
    excluded, arc_of_segment, per_arc = None, {}, {}
    if exclusions is not None:
        excluded = Pieces.from_exclusions(exclusions, pa.radius)
        for arc in exclusions.get("Arcs", []):
            for i in arc.get("Segments", []):
                arc_of_segment[i] = arc["Arc"]
            per_arc[arc["Arc"]] = {"Prediction": arc["Prediction"], "Excluded": collections.Counter(), "WouldBe": collections.Counter(),
                                   "Mismatch": collections.Counter()}
    fwd, pairs_fwd = compare(pa, pb, region, spacing_over_r, tolerance_over_r, "A->B", examples, curve_tolerance_over_r,
                             excluded, arc_of_segment, per_arc)
    rev, pairs_rev = compare(pb, pa, region, spacing_over_r, tolerance_over_r, "B->A", examples, curve_tolerance_over_r,
                             excluded, arc_of_segment, per_arc)
    result = {"Region": list(region), "Offset": list(offset), "SpacingOverR": spacing_over_r, "ToleranceOverR": tolerance_over_r,
              "CurveToleranceOverR": curve_tolerance_over_r,
              "MatchingRadius": [pa.radius, pb.radius], "Forward": dict(fwd), "Reverse": dict(rev),
              "PairsForward": [[list(k), v] for k, v in sorted(pairs_fwd.items(), key=lambda kv: [str(x) for x in kv[0]])],
              "PairsReverse": [[list(k), v] for k, v in sorted(pairs_rev.items(), key=lambda kv: [str(x) for x in kv[0]])]}
    defect = sum(fwd.get(c, {"Length": 0.0})["Length"] for c in ("A", "D")) + sum(rev.get(c, {"Length": 0.0})["Length"] for c in ("A", "D"))
    result["DefectLength"] = defect
    result["Verdict"] = "PASS" if defect == 0.0 else "FAIL"
    result["Summary"] = {d: {c: {"Length": v["Length"], "Samples": v["Samples"]} for c, v in t.items()} for d, t in (("A->B", fwd), ("B->A", rev))}
    if exclusions is not None:
        arcs = []
        for arc_id, record in sorted(per_arc.items()):
            arcs.append({"Arc": arc_id, "Prediction": record["Prediction"], "Excluded": dict(record["Excluded"]),
                         "WouldBe": dict(record["WouldBe"]), "MismatchLength": sum(record["Mismatch"].values())})
        result["CutArcs"] = {"Rule": exclusions.get("Rule"), "ExcludedLength": {d: fwd.get("Excluded", {"Length": 0.0})["Length"] if d == "A->B"
                                                                               else rev.get("Excluded", {"Length": 0.0})["Length"] for d in ("A->B", "B->A")},
                             "Arcs": arcs, "PredictionMismatchLength": sum(a["MismatchLength"] for a in arcs),
                             "PredictionMismatches": [a["Arc"] for a in arcs if a["MismatchLength"] > 0.0]}
    return result


def markdown(result):
    lines = [f"Segment identity over {result['Region']}: **{result['Verdict']}** (class A + D length {result['DefectLength']:.6f})",
             "| direction | class | length | samples |", "|---|---|---|---|"]
    for direction, table in result["Summary"].items():
        for cls in ("Identical", "A", "B", "C", "D", "Excluded"):
            if cls in table:
                lines.append(f"| {direction} | {cls} | {table[cls]['Length']:.4f} | {table[cls]['Samples']} |")
    if "CutArcs" in result:
        ca = result["CutArcs"]
        lines.append(f"\nCut-arc exclusions ((h)-9): excluded A->B {ca['ExcludedLength']['A->B']:.4f} / B->A {ca['ExcludedLength']['B->A']:.4f} um; "
                     f"prediction mismatches {len(ca['PredictionMismatches'])} ({ca['PredictionMismatchLength']:.4f} um)")
        lines += ["| arc | prediction | excluded A->B / B->A | would-be classes | mismatch um |", "|---|---|---|---|---|"]
        for a in ca["Arcs"]:
            lines.append(f"| {a['Arc']} | {a['Prediction']} | {a['Excluded'].get('A->B', 0.0):.3f} / {a['Excluded'].get('B->A', 0.0):.3f} | "
                         f"{ {k: round(v, 3) for k, v in a['WouldBe'].items()} } | {a['MismatchLength']:.3f} |")
    for direction, key in (("A->B", "Forward"), ("B->A", "Reverse")):
        for cls in ("A", "D", "B", "C"):
            entry = result[key].get(cls)
            if not entry or not entry["Examples"]:
                continue
            lines.append(f"\n{direction} class {cls} examples ({len(entry['Examples'])} of {entry['Samples']} samples):")
            for e in entry["Examples"][:10]:
                lines.append(f"- {e['Point']}: {e['Why']}; A {json.dumps(e['A'])}; B {json.dumps(e['B'])}; deltas {json.dumps(e['Deltas'])}")
    return "\n".join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--a", required=True, help="reference manifest (the chip)")
    parser.add_argument("--b", required=True, help="manifest under test (the window)")
    parser.add_argument("--region", nargs=4, type=float, required=True, metavar=("X0", "X1", "Y0", "Y1"))
    parser.add_argument("--offset", nargs=3, type=float, default=(0.0, 0.0, 0.0), help="added to B's coordinates")
    parser.add_argument("--spacing", type=float, default=0.25, help="sample spacing / R")
    parser.add_argument("--tolerance", type=float, default=1e-3, help="perimeter match tolerance / R")
    parser.add_argument("--curve-tolerance", type=float, default=0.5,
                        help="perimeter match tolerance / R for re-meshed curve pieces (chords displaced by up to the sagitta)")
    parser.add_argument("--examples", type=int, default=50)
    parser.add_argument("--exclude", nargs=2, default=None, metavar=("E0_JSON", "WINDOW"),
                        help="window_query E0 output and window name: its E1CutArcExclusions pieces are class Excluded ((h)-9)")
    parser.add_argument("--output", required=True)
    parser.add_argument("--markdown", default=None)
    args = parser.parse_args(argv)
    exclusions = json.load(open(args.exclude[0]))["Windows"][args.exclude[1]]["E1CutArcExclusions"] if args.exclude else None
    result = run(json.load(open(args.a)), json.load(open(args.b)), tuple(args.region), tuple(args.offset), args.spacing,
                 args.tolerance, args.examples, args.curve_tolerance, exclusions)
    result.update({"A": args.a, "B": args.b, "Exclude": list(args.exclude) if args.exclude else None})
    with open(args.output, "w") as out:
        json.dump(result, out, indent=1)
    text = markdown(result)
    if args.markdown:
        with open(args.markdown, "w") as out:
            out.write(text)
    print(text)
    return 0 if result["Verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
