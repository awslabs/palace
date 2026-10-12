#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Placement audit (gate A10, geometry only): every matched cluster / corner / stack / pair
model, mapped through the patch frame the dry run actually used, lands on the feature's
claimed perimeter portions within the signature parameter tolerance (1e-3 R; for a cluster
model the identification's arc-fit tolerance plus the quantum terms, see
cluster_tolerance_over_r).

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
  straight-portion endpoint / arc end lies on the mapped model polyline. The claimed angular
  range is parameterised as the identification parameterises it (Identifier::RunPoint): a
  claim cut at parameter s of a mesh segment lies on the circle at the angle the collinear
  RUN of the segment (its on-circle joints, RunIndex) maps s to linearly - not at the angle
  interpolated between the segment's own keys, which on a design-polygon arc sub-divided by
  the mesh are interior chord vertices up to the chord sagitta inside the circle (decision
  710: 0c94 on C3's arc 42 read 1.37e-3 R of angular excess with the model at the
  identification's claim end to 3.3e-7 R).
* ConvexCorner / ConcaveCorner: arms from the patch origin along +u and at the model
  ``Angle`` counterclockwise in (u, v) (the corner coupon's frame, generate_corner_response
  corner_frame); a rounded corner's claimed arc chords lie on the fillet circle of the
  model ``CornerRadius`` (centre CornerRadius x bisector / sin(Angle / 2)).
* ParallelEdgeCluster (stack): per longitudinal patch, origin + Offset_k u lies on the
  claimed portions of side k (sides in increasing lateral offset, reversed for chirality
  -1, as the patch construction orders them); arc segments read as their mesh chords (the
  patch sits on the mesh) with the chord's sagitta L^2 / (8 r) as an allowance (the design
  offsets are the concentric circles' separations), and the pair rule's own 5 % of the
  offset (PairSeparationToleranceRelative: the admissible variation of a pair's separation).
* SameConductorGap / DifferentConductorGap / SameConductorStrip and their Curved* classes:
  origin -+ Separation / 2 u lies on side 0 / 1 (swapped for chirality -1), chords and the
  5 % likewise. A curvature-family runtime model ``<anchor>@<convexity>-kappa<k>-<rule>``
  resolves to its anchor (same basis, separation and references); an angle-interpolated
  corner ``<base>@corner-angle<deg>-<rule>`` to its base coupon at the interpolated angle
  (the arms checked are the feature's); a first-order node (a
  straight-like bend splits its quadrature between the anchor and the family node at
  kappa = 1 / StraightBendRadiusOverR) is a library model. For every patch of such a curved
  coupon the model's first edge e1 must sit on the side of the coupon's recorded Convexity:
  the inner circle of a gap / the outer circle of a strip for Convex (read on the claimed
  fitted arc under e1; a chord without a fitted arc is not checked).
* A MIRROR-FORMED cluster placed on its REAL half by its requirement contract (decision 557
  (4); the feature's ``Mirror.Contract``, the operator's ``MirrorFormedPlacement`` record):
  the model's REAL edges (the placement record's ``EdgeWeights`` 1.0 = ``RealEdges``, cross-
  checked against the ENTRY stamp ``Edges[].Weight``) are checked as above; its IMAGE edges
  (weight 0.0, never applied to a cell) lie beyond the mirror plane where the device has no
  perimeter by construction and are LABELLED by name (``MirrorFormed.ImageEdges``), not
  checked. The context gate classifies every Context entry end by the contract's mirror
  planes (``Diagnostics.MirrorBand.Planes``: Inside = Offset - n . q, the image side is
  Inside < -tolerance): real-side ends are checked, image-side ends labelled by name
  (``ImageEnds``). A mirror-formed placement whose record lacks the weights or the planes is
  listed, never silently passed (decision 616 MINOR-8; v3-f E-4).

A model whose placement cannot be evaluated from its record (a topology the gate does not
cover, a runtime model the library does not resolve) is listed, never silently passed. This gate is the one a misplaced
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
# Identification.Conventions ArcFitToleranceOverR (kArcFitToleranceOverRadius): an arc's joint
# vertices lie on the fitted circle within this fraction of R, strictly. The identification
# serialises an arc portion's end ON the circle at the junction vertex's angle; the generator's
# ArcJoint snap (cluster_signature_geometry ARC_JOINT_SNAP_OVER_R = this tolerance +
# ARC_JOINT_SNAP_QUANTA quanta, read inclusive to INCLUSIVE_READ_QUANTA) puts the model's end
# on the vertex: a correctly placed arc end reads the vertex's radial fit residual eps
# (decision 708: a596 8.97e-4 R), bounded by the fit tolerance + the snap admission.
ARC_FIT_TOLERANCE_OVER_R = 1.0e-3
# Identification.Conventions SignatureLengthQuantumOverR: the signature grid; a claimed end and a
# model end are each rounded to it (a point rounding is at most sqrt(2) / 2 quantum).
SIGNATURE_QUANTUM_OVER_R = 1.0e-6
ARC_JOINT_SNAP_QUANTA = 2.0
INCLUSIVE_READ_QUANTA = 0.5
SERIALISATION_ROUNDINGS = 2
# Identification.Conventions DirectionQuantum (kDirectionQuantum): two consecutive mesh
# segments of a chain are one collinear run when their tangents' cosine is not below
# 1 - DirectionQuantum on the quantised grid (Identifier::BuildRunIndex).
DIRECTION_QUANTUM = 1.0e-12
# Identification.Conventions PairSeparationToleranceRelative: the bent-pair rule admits a pair
# whose sampled separation varies by this fraction (a slow taper is one pair at its mean
# separation), so a pair / stack model's edge lines may sit this far from the mesh edges by
# the rule's own definition; a taper wider than that about its mean fails the placement gate
# by construction (the coupon at the mean separation is placed where the local separation
# differs by more) - a recorded limitation of "one separation per pair", reported as such.
PAIR_SEPARATION_TOLERANCE = 0.05

CLUSTER_TYPES = ("SpatialEdgeCluster",)
CORNER_TYPES = ("ConvexCorner", "ConcaveCorner")
STACK_TYPES = ("ParallelEdgeCluster",)
PAIR_TYPES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip",
              "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip")
STRIP_TYPES = ("SameConductorStrip", "CurvedSameConductorStrip")


def resolve_model(models, name):
    """The library model a dry-run patch refers to: the model of that name; for a
    curvature-family runtime model `<anchor>@<convexity>-kappa<k>-<rule>` its straight anchor
    (the blend shares the anchor's basis, separation and conductor references) together with
    the blend's convexity; for an angle-interpolated corner `<base>@corner-angle<deg>-<rule>`
    (USER decision 121 (C)) a copy of its base coupon with `Angle` = the interpolated angle
    (the blend is placed in the FEATURE's frame at the feature's angle: the arms of the
    placement check are the interpolated corner's, not the base coupon's). Returns
    (model, convexity or None)."""
    model = models.get(name)
    if model is not None:
        return model, model.get("Convexity")
    if "@" in name:
        anchor, tail = name.rsplit("@", 1)
        model = models.get(anchor)
        if model is not None and "-kappa" in tail:
            convexity = tail.split("-kappa", 1)[0]
            return model, convexity.capitalize()
        if model is not None and tail.startswith("corner-angle"):
            angle = float(tail[len("corner-angle"):].rsplit("-", 1)[0])
            interpolated = dict(model)
            interpolated["Angle"] = angle
            interpolated["CornerRadius"] = 0.0
            return interpolated, model.get("Convexity")
    return None, None


def pair_side_convexity(claimed, row, separation, strip):
    """Convexity of the model's first edge e1 as placed by a pair patch, read on the claimed
    fitted arc under e1 (metal inside the bend = Convex: e1 on the inner circle of a gap or
    the outer circle of a strip); None when e1 lies on no fitted arc (a straight chord)."""
    origin, axes = row["Origin"], row["Axes"]
    e1 = origin - 0.5 * separation * axes[0]
    e2 = origin + 0.5 * separation * axes[0]
    best = None
    for entry in claimed.arc_ranges.values():
        d = abs(float(np.linalg.norm(e1 - entry["center"])) - entry["radius"])
        if best is None or d < best[0]:
            best = (d, entry)
    if best is None or best[0] > PAIR_SEPARATION_TOLERANCE * separation:
        return None
    center = best[1]["center"]
    e1_inner = float(np.linalg.norm(e1 - center)) < float(np.linalg.norm(e2 - center))
    return "Convex" if e1_inner != strip else "Concave"


def gate(name, passed, detail, evaluable=True):
    status = "PASS" if passed else ("FAIL" if evaluable else "NOT-EVALUABLE")
    return {"Gate": name, "Status": status, "Detail": detail}


MIRROR_FORMED_RULE = ("decision 557 (4) / impl-B5 CONTRACT.md section 4: a mirror-formed cluster is placed ONCE on its REAL half; the "
                      "model's real edges (EdgeWeights 1.0) are checked on the claimed portions, its image edges (0.0, never "
                      "applied to a cell) lie beyond the mirror plane by construction and are labelled, not checked; Context "
                      "entry ends on the image side of the contract's planes are labelled, not checked (v3-f E-4)")


def mirror_formed_placement(identification, feature, model, row):
    """The real / image split of a mirror-formed cluster placed on its real half by its
    requirement contract (the feature's Mirror.Contract): (split, reason). split = {RealEdges,
    ImageEdges (model Edge indices), Planes (indices into Diagnostics.MirrorBand.Planes)} or
    None; reason names why the split cannot be read when the feature IS such a placement
    (None for an ordinary feature). The weights are the operator's MirrorFormedPlacement patch
    record (Feature + Patch), cross-checked against the model's ENTRY stamp Edges[].Weight;
    the stamp alone when the record is absent (an older manifest). The stamp is read only
    when EVERY model edge carries a Weight (the C++ record stamps one weight per edge): a
    partial stamp is not evaluable, by name."""
    contract = (feature.get("Mirror") or {}).get("Contract")
    if not isinstance(contract, dict):
        return None, None
    edges = model.get("Edges") or []
    stamped_edges = sum(1 for edge in edges if "Weight" in edge)
    if 0 < stamped_edges < len(edges):
        return None, f"the model's ENTRY stamp Edges[].Weight is partial ({stamped_edges} of {len(edges)} edges stamped): every edge stamped or none"
    stamped = bool(edges) and stamped_edges == len(edges)
    stamp = [float(edge["Weight"]) for edge in edges] if stamped else []
    placement = None
    patches = ((identification.get("Diagnostics") or {}).get("MirrorFormedPlacement") or {}).get("Patches") or []
    for entry in patches:
        if int(entry.get("Feature", -1)) == int(feature["Id"]) and (row is None or int(entry.get("Patch", -1)) == int(row["Patch"])):
            placement = entry
            break
    if placement is not None:
        weights = [float(w) for w in placement["EdgeWeights"]]
        if stamped and weights != stamp:
            return None, "the MirrorFormedPlacement EdgeWeights differ from the model's ENTRY stamp Edges[].Weight"
    elif stamped:
        weights = stamp
    else:
        return None, "mirror-formed placement without a MirrorFormedPlacement record or an ENTRY stamp (Edges[].Weight)"
    if model.get("Edges") and len(weights) != len(model["Edges"]):
        return None, f"{len(weights)} edge weights for {len(model['Edges'])} model Edges"
    real = [k for k, w in enumerate(weights) if w > 0.0]
    image = [k for k, w in enumerate(weights) if w <= 0.0]
    planes = [int(k) for k in (contract.get("Planes") or (feature.get("Mirror") or {}).get("Planes") or [])]
    return {"RealEdges": real, "ImageEdges": image, "Planes": planes, "Weights": weights,
            "RealPortions": contract.get("RealPortions"), "Record": "MirrorFormedPlacement" if placement is not None else "EntryStamp"}, None


def mirror_plane_inside(identification, plane_indices):
    """Signed distance into the domain (>= 0 real side) of a mesh point from the nearest of the
    given mirror planes (surfaceresponsemirror.hpp MirrorPlane::Inside: Offset - n . p), read
    from Diagnostics.MirrorBand.Planes; None when the record lacks a named plane."""
    planes = ((identification.get("Diagnostics") or {}).get("MirrorBand") or {}).get("Planes") or []
    selected = []
    for k in plane_indices:
        if not 0 <= k < len(planes) or "Normal" not in planes[k] or "Offset" not in planes[k]:
            return None
        selected.append((np.asarray(planes[k]["Normal"], dtype=float), float(planes[k]["Offset"])))
    if not selected:
        return None

    def inside(point):
        return min(offset - float(np.dot(normal, point)) for normal, offset in selected)

    return inside


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


def _arc_chord_distance(point, a, b, center, radius):
    """Distance from a point to the arc of the fitted circle (centre, radius) that the chord
    a-b subtends (surfaceresponseoperator ArcChordDistance; block (b) DESIGN section 2 (a)):
    the point's in-plane projection inside the chord's angular interval reads the distance to
    the circle (radial residual + out-of-plane offset), outside it the distance to the nearer
    chord end; a chord collinear with the centre falls back to the straight distance."""
    ra, rb = a - center, b - center
    n = np.cross(ra, rb)
    n_norm = float(np.linalg.norm(n))
    if n_norm <= 1.0e-12 * max(float(np.linalg.norm(ra)) * float(np.linalg.norm(rb)), 1.0e-300):
        return _segment_distance(point, a, b)
    n = n / n_norm
    rq = point - center
    out_of_plane = float(np.dot(rq, n))
    in_plane = rq - out_of_plane * n
    within = float(np.dot(np.cross(ra, in_plane), n)) >= 0.0 and float(np.dot(np.cross(in_plane, rb), n)) >= 0.0
    if not within:
        return min(float(np.linalg.norm(point - a)), float(np.linalg.norm(point - b)))
    return math.hypot(abs(float(np.linalg.norm(in_plane)) - radius), out_of_plane)


def device_perimeter_distance(identification, point):
    """Distance from a point (mesh units) to the device perimeter of a manifest Identification
    (surfaceresponseoperator DevicePerimeterDistance): the Segments' Keys as straight chords, a
    segment on a fitted arc (Segments[].Arc) read on its circle's arc."""
    point = np.asarray(point, dtype=float)
    arcs = identification.get("Arcs", [])
    best = math.inf
    for segment in identification["Segments"]:
        key = np.asarray(segment["Key"], dtype=float)
        if "Arc" in segment and 0 <= int(segment["Arc"]) < len(arcs):
            arc = arcs[int(segment["Arc"])]
            best = min(best, _arc_chord_distance(point, key[0], key[1], np.asarray(arc["Center"], dtype=float), float(arc["Radius"])))
        else:
            best = min(best, _segment_distance(point, key[0], key[1]))
    return best


def context_gate(identification, patches, library, radius):
    """Gate A10-context (decision 282 section 3, arc-aware per block (b) DESIGN section 2
    (a)): every Context entry END of every placed contract-3 cluster model (the model's
    Signature placed by the dry-run patch frame) lies on the device perimeter within the
    cluster tolerance (cluster_tolerance_over_r: the fit tolerance + the quantum terms);
    returns (gate, entries)."""
    base_tolerance_over_r, base_terms = cluster_tolerance_over_r(identification)
    radius_library = float(library["MatchingRadius"])
    scale = radius / radius_library
    models = {m["Name"]: m for m in library["Models"]}
    features = {f["Id"]: f for f in identification["Features"]}
    entries = []
    for row in patches:
        feature = features.get(row["Feature"])
        if feature is None or feature["Type"] not in CLUSTER_TYPES or feature["Match"]["Status"] != "Matched":
            continue
        model, _ = resolve_model(models, row["Model"])
        signature = (model or {}).get("Signature") or {}
        if "Box" not in signature or (feature.get("Match") or {}).get("LegacyContract"):
            continue
        tolerance_over_r, _ = cluster_tolerance_over_r(identification, feature)
        tol = tolerance_over_r * radius
        origin, axes = row["Origin"], row["Axes"]
        entry = {"Feature": row["Feature"], "Model": row["Model"], "Checks": 0, "WorstOverR": 0.0, "ToleranceOverR": tolerance_over_r, "Defects": [], "ArcEnds": 0}
        split, reason = mirror_formed_placement(identification, feature, model or {}, row)
        inside = None
        if split is not None or reason is not None:
            entry["ImageEnds"] = []
            inside = mirror_plane_inside(identification, split["Planes"]) if split is not None else None
            if inside is None:
                entry["Defects"].append({"Point": "mirror planes", "DeviationOverR": math.inf, "Mapped": None,
                                         "Reason": reason or "no Diagnostics.MirrorBand.Planes record for the contract's Planes "
                                                             f"{(split or {}).get('Planes')}"})
                entry["Checks"] += 1
        for index, piece in enumerate(signature.get("Context", [])):
            P = [float(v) for v in piece["P"]]
            for label, (x, y) in (("begin", (P[0], P[1])), ("end", (P[2], P[3]))):
                q = origin + (scale * np.asarray([x * radius_library, y * radius_library, 0.0])) @ axes
                if inside is not None and inside(q) < -tol:
                    # The image side of the contract's mirror plane: the device has no
                    # perimeter there by construction; labelled, not checked.
                    entry["ImageEnds"].append(f"context {index} {label}")
                    continue
                distance = device_perimeter_distance(identification, q)
                entry["Checks"] += 1
                entry["WorstOverR"] = max(entry["WorstOverR"], distance / radius)
                if "Arc" in piece:
                    entry["ArcEnds"] += 1
                if distance > tol and len(entry["Defects"]) < 8:
                    entry["Defects"].append({"Point": f"context {index} {label}", "DeviationOverR": distance / radius, "Mapped": [float(v) for v in q]})
        entries.append(entry)
    defects = [e for e in entries if e["Defects"]]
    detail = {"Features": len(entries), "Checks": sum(e["Checks"] for e in entries), "ArcEnds": sum(e["ArcEnds"] for e in entries),
              "WorstDeviationOverR": max((e["WorstOverR"] for e in entries), default=None), "ToleranceOverR": base_tolerance_over_r,
              "ToleranceTerms": base_terms, "ToleranceRule": CLUSTER_TOLERANCE_RULE,
              "QuantumNearMatchedFeatures": sum(1 for e in entries if e["ToleranceOverR"] > base_tolerance_over_r),
              "MaxToleranceOverR": max((e["ToleranceOverR"] for e in entries), default=base_tolerance_over_r),
              "FeaturesWithDefects": len(defects),
              "MirrorFormedFeatures": sum(1 for e in entries if "ImageEnds" in e), "ImageEnds": sum(len(e.get("ImageEnds", [])) for e in entries),
              "Examples": [{k: v for k, v in e.items() if k != "Checks"} for e in defects[:6]],
              "Basis": "every Context entry end of a placed contract-3 cluster model lies on the device perimeter (segments as chords, fitted-arc segments on their circle's arc) within the cluster tolerance (the fit tolerance + the quantum terms): the C++ A10-extended check, mirrored; a mirror-formed placement's ends on the image side of its contract's planes are labelled (ImageEnds), not checked"}
    return gate("A10-context", not defects, detail, evaluable=bool(entries)), entries


def cluster_tolerance_over_r(identification, feature=None):
    """The A10 clusters / context tolerance of a placed cluster model in units of R (decisions
    708 / 710: the identification's fit tolerance plus the quantum terms): the arc-fit
    tolerance (a correctly placed arc end reads its junction vertex's radial fit residual,
    strictly below it) + the generator's ArcJoint snap admission beyond it (2 quanta, read
    inclusive to half a quantum) + the two serialisation roundings (the claimed end and the
    model end on the signature grid, sqrt(2) / 2 quantum each) + the quantum near-match term
    (sqrt(2) x Match.QuantumNearMatch.MaxDeltaQuanta quanta: the model's numbers differ from
    the feature's by that many quanta per coordinate; 0 for an exact match). The constants
    are the manifest's Identification.Conventions where recorded. Returns (tolerance, terms)."""
    conventions = identification.get("Conventions") or {}
    fit = float(conventions.get("ArcFitToleranceOverR", ARC_FIT_TOLERANCE_OVER_R))
    quantum = float(conventions.get("SignatureLengthQuantumOverR", SIGNATURE_QUANTUM_OVER_R))
    near_match = ((feature or {}).get("Match") or {}).get("QuantumNearMatch") or {}
    delta_quanta = float(near_match.get("MaxDeltaQuanta", 0.0))
    terms = {"FitToleranceOverR": fit,
             "JointSnapAdmissionOverR": (ARC_JOINT_SNAP_QUANTA + INCLUSIVE_READ_QUANTA) * quantum,
             "SerialisationRoundingOverR": SERIALISATION_ROUNDINGS * (math.sqrt(2.0) / 2.0) * quantum,
             "QuantumNearMatchOverR": math.sqrt(2.0) * delta_quanta * quantum,
             "QuantumNearMatchMaxDeltaQuanta": delta_quanta}
    return terms["FitToleranceOverR"] + terms["JointSnapAdmissionOverR"] + terms["SerialisationRoundingOverR"] + terms["QuantumNearMatchOverR"], terms


CLUSTER_TOLERANCE_RULE = ("decisions 708 / 710 (the chord-vs-arc class 611 / 653 / 655 / 663 / 708-709): the identification's arc-fit tolerance "
                          "(a placed arc end reads its junction vertex's radial fit residual, strictly below it) + the generator's ArcJoint snap "
                          "admission (2 quanta + the half-quantum inclusive read) + the two serialisation roundings (sqrt(2) quanta) + sqrt(2) x "
                          "MaxDeltaQuanta quanta for a quantum near-match; the claimed arc range is parameterised by the RUN's on-circle joints "
                          "(Identifier::RunPoint), not the mesh segment's keys")


class RunIndex:
    """The collinear RUNS of the fitted-arc mesh segments of an Identification, as the
    identification builds them (Identifier::BuildRunIndex / BuildRunArcGeometry): the segments
    of one chain on one arc, linked through shared keys and split where consecutive tangents
    are not collinear on the DirectionQuantum grid; a run's parameter s (the chord length from
    its start, every segment's [t0, t1] scaled so that the sum of the segment lengths is the
    chord) maps linearly onto the angle between the run's two end joints on the circle
    (RunPoint). A run is confined to one arc, as BuildRunArcGeometry requires of a run it
    parameterises on the circle. A segment on no run of more than one member is its own run
    (the mesh-resolved arc: the run convention and the segment-key convention coincide)."""

    def __init__(self, identification, radius):
        conventions = identification.get("Conventions") or {}
        self.direction_quantum = float(conventions.get("DirectionQuantum", DIRECTION_QUANTUM))
        self.segments = identification["Segments"]
        self.grid = 1.0e-9 * radius  # shared keys are the same mesh vertex: equal to roundoff
        self.runs = {}  # segment index -> {"members": {segment: (forward, t0, t1)}, "start", "end", "length"}
        groups = defaultdict(list)
        for index, segment in enumerate(self.segments):
            if "Arc" in segment:
                groups[(segment.get("Chain"), int(segment["Arc"]))].append(index)
        for members in groups.values():
            for path in self._paths(members):
                for run in self._split_collinear(path):
                    self._register(run)

    def _key(self, point):
        return tuple(int(round(float(x) / self.grid)) for x in point)

    def _paths(self, members):
        """The simple paths of the member segments through their shared keys, each as an
        ordered list of (segment, forward)."""
        ends = defaultdict(list)
        for index in members:
            for end in (0, 1):
                ends[self._key(self.segments[index]["Key"][end])].append(index)
        unused = set(members)
        paths = []
        while unused:
            start = None
            for index in sorted(unused):
                for end in (0, 1):
                    key = self._key(self.segments[index]["Key"][end])
                    if sum(1 for other in ends[key] if other in unused) == 1:
                        start = key
                        break
                if start is not None:
                    break
            if start is None:  # a closed loop of arc segments: start anywhere
                start = self._key(self.segments[min(unused)]["Key"][0])
            current, path = start, []
            while True:
                following = [other for other in ends[current] if other in unused]
                if not following:
                    break
                index = following[0]
                unused.discard(index)
                forward = self._key(self.segments[index]["Key"][0]) == current
                path.append((index, forward))
                current = self._key(self.segments[index]["Key"][1 if forward else 0])
            paths.append(path)
        return paths

    def _direction(self, index, forward):
        key = np.asarray(self.segments[index]["Key"], dtype=float)
        d = key[1] - key[0] if forward else key[0] - key[1]
        return d / np.linalg.norm(d)

    def _split_collinear(self, path):
        runs, run = [], [path[0]]
        for previous, current in zip(path, path[1:]):
            cosine = float(np.dot(self._direction(*previous), self._direction(*current)))
            if round(cosine / self.direction_quantum) >= round((1.0 - self.direction_quantum) / self.direction_quantum):
                run.append(current)
            else:
                runs.append(run)
                run = [current]
        runs.append(run)
        return runs

    def _register(self, run):
        first_index, first_forward = run[0]
        last_index, last_forward = run[-1]
        start = np.asarray(self.segments[first_index]["Key"][0 if first_forward else 1], dtype=float)
        end = np.asarray(self.segments[last_index]["Key"][1 if last_forward else 0], dtype=float)
        length = float(np.linalg.norm(end - start))
        total = sum(float(self.segments[index]["Length"]) for index, _ in run)
        ratio = length / total if total > 0.0 else 0.0
        record = {"members": {}, "start": start, "end": end, "length": length}
        t = 0.0
        for index, forward in run:
            segment_length = float(self.segments[index]["Length"])
            record["members"][index] = (forward, t * ratio, (t + segment_length) * ratio)
            t += segment_length
        for index, _ in run:
            self.runs[index] = record

    def run_segments(self, index):
        """The segment indices of the run holding segment ``index`` (in run order)."""
        run = self.runs.get(index)
        return sorted(run["members"], key=lambda k: run["members"][k][1]) if run else [index]

    def arc_point(self, index, s, center, radius, normal):
        """The point of the fitted circle (center, radius; plane normal) at a claim cut s along
        mesh segment ``index`` from its Key[0], parameterised as the identification does
        (RunPoint: the run parameter maps linearly onto the angle between the run's end
        joints); None when the segment is on no run or the run is degenerate (no chord, no
        sweep), where the caller reads the segment's own keys."""
        run = self.runs.get(index)
        if run is None or run["length"] <= 0.0:
            return None
        forward, t0, t1 = run["members"][index]
        segment_length = float(self.segments[index]["Length"])
        fraction = (s / segment_length if forward else 1.0 - s / segment_length) if segment_length > 0.0 else 0.0
        s_run = t0 + fraction * (t1 - t0)
        d0 = run["start"] - center
        d0 = d0 - float(np.dot(d0, normal)) * normal
        d1 = run["end"] - center
        d1 = d1 - float(np.dot(d1, normal)) * normal
        if float(np.linalg.norm(d0)) <= 1.0e-12 * radius or float(np.linalg.norm(d1)) <= 1.0e-12 * radius:
            return None
        u = d0 / np.linalg.norm(d0)
        v = np.cross(normal, u)
        sweep = math.atan2(float(np.dot(d1, v)), float(np.dot(d1, u)))
        if abs(sweep) <= 1.0e-12:
            return None
        angle = sweep * s_run / run["length"]
        return center + radius * (math.cos(angle) * u + math.sin(angle) * v)


class ClaimedGeometry:
    """The claimed perimeter portions of one feature in mesh coordinates: straight
    sub-segments (per mesh segment) and, per fitted arc, the circle with the angular range
    of the claimed chords (the plane of the arc is the process plane of the feature); a claim
    cut inside an arc run is placed on the circle by the run convention (RunIndex)."""

    def __init__(self, feature, segments, arcs, normal, runs=None):
        self.normal = np.asarray(normal, dtype=float)
        self.normal /= np.linalg.norm(self.normal)
        self.straight = []          # (p0, p1, side)
        self.chords = []            # (p0, p1, side): the claimed mesh chords of arc segments
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
            self.chords.append((p0, p1, int(side)))
            index = int(segment["Arc"])
            arc = arcs[index]
            entry = self.arc_ranges.get(index)
            if entry is None:
                center = np.asarray(arc["Center"], dtype=float)
                # In-plane basis of the arc: u toward the first claimed point, v = n x u.
                u = key[0] - center
                u -= float(np.dot(u, self.normal)) * self.normal
                u /= np.linalg.norm(u)
                v = np.cross(self.normal, u)
                entry = {"center": center, "radius": float(arc["Radius"]), "u": u, "v": v, "angles": [], "side": int(side), "ends": {}}
                self.arc_ranges[index] = entry
            # A claim cut inside a chord is a point of the ARC (the identification maps the
            # run parameter linearly onto the angle), not of the chord: the claimed end is
            # the circle point at the angle the RUN of collinear mesh segments between two
            # on-circle joints maps the cut to (Identifier::RunPoint; decision 710), read
            # between the segment's own keys only where the segment is on no run (a chord's
            # interior lies up to its sagitta inside the circle).
            a0, a1 = self._angle(entry, key[0]), self._angle(entry, key[1])
            sweep = math.remainder(a1 - a0, 2.0 * math.pi)
            for s in (float(s0), float(s1)):
                q = runs.arc_point(int(seg), s, entry["center"], entry["radius"], self.normal) if runs is not None else None
                if q is None:
                    angle = a0 + sweep * s / length
                    q = entry["center"] + entry["radius"] * (math.cos(angle) * entry["u"] + math.sin(angle) * entry["v"])
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

    def distance(self, point, side=None, chords=False):
        """Distance from a point to the claimed geometry (optionally of one side). Arc
        segments are read on their fitted circle (a cluster model is built on the design
        circle) or, with ``chords``, on the claimed mesh chords (a longitudinal patch sits on
        the mesh: its origin is a quadrature point of the chord, up to the sagitta inside the
        circle)."""
        best = math.inf
        for p0, p1, s in self.straight:
            if side is None or s == side:
                best = min(best, _segment_distance(point, p0, p1))
        if chords:
            for p0, p1, s in self.chords:
                if side is None or s == side:
                    best = min(best, _segment_distance(point, p0, p1))
            return best
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
    signature_tol = SIGNATURE_TOLERANCE_OVER_R * radius
    cluster_tolerance_base, cluster_terms = cluster_tolerance_over_r(identification)
    radius_library = float(library["MatchingRadius"])
    scale = radius / radius_library
    models = {m["Name"]: m for m in library["Models"]}
    features = {f["Id"]: f for f in identification["Features"]}
    segments = identification["Segments"]
    arcs = identification.get("Arcs", [])
    runs = RunIndex(identification, radius)
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
        # A feature's patches may refer to several runtime models (a straight-like bend
        # splits its quadrature between the anchor and the family's first-order node); every
        # one must resolve to a library model (a curvature-family blend to its anchor).
        model_names = sorted({r["Model"] for r in rows})
        model_name = model_names[0]
        resolved = {name: resolve_model(models, name) for name in model_names}
        missing = [name for name, (m, _) in resolved.items() if m is None]
        if missing:
            # The dry run wrote the model's name; a library without it cannot be audited.
            not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": ", ".join(missing), "Reason": "model not in the library"})
            continue
        model = resolved[model_name][0]
        normal = np.asarray(feature["Frame"]["Axes"][2], dtype=float)
        claimed = ClaimedGeometry(feature, segments, arcs, normal, runs)
        entry = {"Feature": feature_id, "Type": kind, "Model": model_name, "WorstOverR": 0.0, "Checks": 0, "Defects": []}
        # A cluster model's tolerance is the identification's fit tolerance + the quantum
        # terms (cluster_tolerance_over_r); every other model's the signature tolerance.
        if kind in CLUSTER_TYPES:
            entry["ToleranceOverR"], _ = cluster_tolerance_over_r(identification, feature)
            tol = entry["ToleranceOverR"] * radius
        else:
            tol = signature_tol

        def record(point_label, deviation, mapped=None, allowance=0.0):
            entry["Checks"] += 1
            entry["WorstOverR"] = max(entry["WorstOverR"], deviation / radius)
            if deviation > tol + allowance and len(entry["Defects"]) < 8:
                entry["Defects"].append({"Point": point_label, "DeviationOverR": deviation / radius, "AllowanceOverR": allowance / radius, "Mapped": None if mapped is None else [float(x) for x in mapped]})

        def chord_sagitta(row):
            """A longitudinal patch on a mesh chord of a fitted arc: the chord lies up to its
            sagitta L^2 / (8 r) inside the design circle, and the concentric sides' chords are
            offset by the same order; the model's design offsets are checked to within it."""
            if row["Segment"] < 0:
                return 0.0
            segment = segments[row["Segment"]]
            if "Arc" not in segment:
                return 0.0
            length = float(segment["Length"])
            return length * length / (8.0 * float(arcs[int(segment["Arc"])]["Radius"]))

        if kind in CLUSTER_TYPES:
            pairs, source = model_cluster_edges(model, radius_library)
            if pairs is None:
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name, "Reason": "no Edges and no Signature"})
                continue
            if len(rows) != 1:
                entry["Defects"].append({"Point": "patch count", "DeviationOverR": math.inf, "Mapped": None, "Count": len(rows)})
            row = rows[0]
            split, reason = mirror_formed_placement(identification, feature, model, row)
            if reason is not None:
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name, "Reason": f"mirror-formed placement: {reason}"})
                continue
            if split is not None and source != "Edges":
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name,
                                      "Reason": "mirror-formed placement by a Signature-only model: no Edges to split real / image"})
                continue
            if split is not None and len(split["Weights"]) != len(pairs):
                not_evaluable.append({"Feature": feature_id, "Type": kind, "Model": model_name,
                                      "Reason": f"mirror-formed placement: {len(split['Weights'])} edge weights for {len(pairs)} model edges"})
                continue
            origin, axes = row["Origin"], row["Axes"]
            mapped_pairs = []
            image_edges = []
            for k, (a, b) in enumerate(pairs):
                ma = origin + (scale * a) @ axes
                mb = origin + (scale * b) @ axes
                if split is not None and k in split["ImageEdges"]:
                    # An image edge (weight 0.0, never applied to a cell) lies beyond the
                    # mirror plane where the device has no perimeter: labelled, not checked.
                    image_edges.append({"Edge": k, "Mapped": [[float(x) for x in ma], [float(x) for x in mb]]})
                    continue
                mapped_pairs.append((ma, mb))
                for label, q in ((f"edge {k} begin", ma), (f"edge {k} end", mb)):
                    record(label, claimed.distance(q), q)
            # Reverse: every claimed straight endpoint / arc end on the mapped model polyline
            # (the real edges' for a mirror-formed placement: the claims are its real portions).
            if not mapped_pairs:
                entry["Defects"].append({"Point": "real edges", "DeviationOverR": math.inf, "Mapped": None, "Reason": "no real edge to place"})
                entry["Checks"] += 1
            for j, q in enumerate(claimed.points):
                record(f"claimed point {j}", min((_segment_distance(q, ma, mb) for ma, mb in mapped_pairs), default=math.inf), q)
            entry["EdgeSource"] = source
            entry["ModelEdges"] = len(pairs)
            if split is not None:
                entry["MirrorFormed"] = {"Record": split["Record"], "RealEdges": split["RealEdges"], "ImageEdges": image_edges,
                                         "RealPortions": split["RealPortions"], "Planes": split["Planes"], "Rule": MIRROR_FORMED_RULE}
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
                        record(f"patch {row['Patch']} edge {k}", claimed.distance(q, side=ordered[k], chords=True), q,
                               max(chord_sagitta(row), PAIR_SEPARATION_TOLERANCE * offset - tol))
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
            for name, (other, _) in resolved.items():
                if abs(float(other.get("Separation", separation / scale)) * scale - separation) > tol:
                    entry["Defects"].append({"Point": "model separations", "DeviationOverR": math.inf, "Mapped": None, "Models": model_names})
                    entry["Checks"] += 1
            sides = claimed.sides()
            chirality = int(feature.get("Chirality", 1))
            ordered = sides if chirality >= 0 else list(reversed(sides))
            strip = kind in STRIP_TYPES
            convexity_checks = 0
            convexity_not_evaluable = 0
            if len(ordered) != 2:
                entry["Defects"].append({"Point": "side count", "DeviationOverR": math.inf, "Mapped": None, "Sides": len(ordered)})
                entry["Checks"] += 1
            else:
                for row in rows:
                    origin, axes = row["Origin"], row["Axes"]
                    for k, sign in enumerate((-1.0, 1.0)):
                        q = origin + sign * 0.5 * separation * axes[0]
                        record(f"patch {row['Patch']} side {k}", claimed.distance(q, side=ordered[k], chords=True), q,
                               max(chord_sagitta(row), PAIR_SEPARATION_TOLERANCE * 0.5 * separation - tol))
                    # Curved coupon on a bent pair (a curvature family blend or a first-order
                    # node): the model's e1 must sit on the side of its recorded convexity
                    # (inner circle of a gap / outer circle of a strip for Convex).
                    expected = resolved[row["Model"]][1]
                    if expected in ("Convex", "Concave"):
                        placed = pair_side_convexity(claimed, row, separation, strip)
                        if placed is not None:
                            convexity_checks += 1
                            entry["Checks"] += 1
                            if placed != expected and len(entry["Defects"]) < 8:
                                entry["Defects"].append({"Point": f"patch {row['Patch']} convexity", "DeviationOverR": math.inf, "Mapped": [float(x) for x in origin], "Model": row["Model"], "Expected": expected, "Placed": placed})
                        else:
                            # e1 is not on a fitted arc within 5 % of the separation (a chorded
                            # bend without an arc: the straight-like first-order patches of a
                            # polyline meander): no convexity check there (review J/C4 m6).
                            convexity_not_evaluable += 1
            entry["Patches"] = len(rows)
            entry["Models"] = model_names
            entry["ConvexityChecks"] = convexity_checks
            entry["ConvexityNotEvaluable"] = convexity_not_evaluable
            results["Pairs"].append(entry)
        else:
            continue  # isolated / curved edges: the patch is the quadrature point itself

    gates = []
    summary = {"ToleranceOverR": SIGNATURE_TOLERANCE_OVER_R, "ClusterToleranceOverR": cluster_tolerance_base, "ClusterToleranceTerms": cluster_terms,
               "Scale": scale, "NotEvaluable": not_evaluable}
    for name, key in (("A10-placement-clusters", "Clusters"), ("A10-placement-corners", "Corners"), ("A10-placement-stacks", "Stacks"), ("A10-placement-pairs", "Pairs")):
        entries = results[key]
        defects = [e for e in entries if e["Defects"]]
        worst = max((e["WorstOverR"] for e in entries), default=None)
        detail = {
            "Features": len(entries),
            "Checks": sum(e["Checks"] for e in entries),
            "WorstDeviationOverR": worst,
            "ToleranceOverR": cluster_tolerance_base if key == "Clusters" else SIGNATURE_TOLERANCE_OVER_R,
            "Models": sorted({m for e in entries for m in e.get("Models", [e["Model"]])}),
            "FeaturesWithDefects": len(defects),
            "ConvexityChecks": sum(e.get("ConvexityChecks", 0) for e in entries),
            "ConvexityNotEvaluable": sum(e.get("ConvexityNotEvaluable", 0) for e in entries),
            "MirrorFormedFeatures": sum(1 for e in entries if "MirrorFormed" in e),
            "ImageEdges": sum(len(e["MirrorFormed"]["ImageEdges"]) for e in entries if "MirrorFormed" in e),
            "Examples": [{k: v for k, v in e.items() if k != "Checks"} for e in defects[:6]],
            "Basis": "model geometry (library units x R_mesh / R_library) mapped through the dry-run patch frame lies on the feature's claimed portions (straight segments; arcs on their fitted circle within the claimed range, the range parameterised by the run's on-circle joints as the identification does) within the tolerance: for a cluster model the fit tolerance + the quantum terms (ToleranceTerms), else the signature parameter tolerance; stack / pair patches on the mesh chords within the chord sagitta and the pair rule's 5 % of the offset (a slow taper wider than 5 % about its mean fails by construction: the coupon at the mean separation); a mirror-formed cluster's real edges only, its image edges labelled (MirrorFormed.ImageEdges)",
        }
        if key == "Clusters":
            detail["ToleranceTerms"] = cluster_terms
            detail["ToleranceRule"] = CLUSTER_TOLERANCE_RULE
            detail["QuantumNearMatchedFeatures"] = sum(1 for e in entries if e["ToleranceOverR"] > cluster_tolerance_base)
            detail["MaxToleranceOverR"] = max((e["ToleranceOverR"] for e in entries), default=cluster_tolerance_base)
        gates.append(gate(name, not defects, detail, evaluable=bool(entries)))
        summary[key] = {"Features": len(entries), "WorstDeviationOverR": worst, "FeaturesWithDefects": len(defects)}
    if not_evaluable:
        gates.append(gate("A10-placement-evaluable", False, {"NotEvaluable": not_evaluable[:20], "Count": len(not_evaluable)}))
    context, context_entries = context_gate(identification, patches, library, radius)
    gates.append(context)
    summary["Context"] = {"Features": len(context_entries), "WorstDeviationOverR": context["Detail"]["WorstDeviationOverR"],
                          "FeaturesWithDefects": sum(1 for e in context_entries if e["Defects"])}
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
