#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Geometry-identification audit: gates A1/A2 on a preflight manifest against the mesh.

    python3 -m surface_response_identification.audit --mesh M.msh2 --config C.json \\
        --manifest postpro/surface-response-requirements.json [--log palace.log] \\
        [--compare OTHER-requirements.json] [--output-prefix out/name]

Writes <prefix>.json and <prefix>.md; exit status 1 when any gate fails. The perimeter E is
recomputed from the mesh (perimeter.py); the manifest supplies the aggregate assignment.
With a version-2 manifest (top-level "Identification": per-segment assignment table, vertex
table, feature signatures, exclusions; palace/models/SURFACE-RESPONSE-IDENTIFICATION.md) the
gates are exact: every segment key is matched against the mesh perimeter, the portions must
cover every assigned segment exactly once, every corner / endpoint / junction vertex must carry
exactly one feature, and the excluded classes must be recorded with their length. Version-1
manifests fall back to the aggregate reading, where gates the manifest cannot support are
NOT-EVALUABLE and count as failures (decision 73(4)).

Patch dry run (phase 4): when ``surface-response-patches.csv`` lies next to the manifest (or is
given with --patches) the gates A7 check the Features-driven patch construction: every portion
of every matched feature is integrated by exactly one longitudinal quadrature (weights summing
to one), every vertex / cluster feature carries one patch, the patched feature set equals the
matched manifest features, and no patch touches an unmatched feature or an excluded segment.
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import manifest as M  # noqa: E402
from . import perimeter as P  # noqa: E402
from .msh2 import read_msh2  # noqa: E402

TURN_BINS = [0.0, 1.0e-6, 1.0, 2.0, 5.0, 10.0, 20.0, 29.999999, 30.000001, 45.0, 60.0, 89.999999, 90.000001, 120.0, 150.0, 180.0]


def gate(name, passed, detail, evaluable=True):
    status = "PASS" if passed else ("FAIL" if evaluable else "NOT-EVALUABLE")
    return {"Gate": name, "Status": status, "Detail": detail}


def perimeter_census(perimeter, radius, targets):
    by_kind = Counter(e.kind for e in perimeter.edges)
    length_by_kind = {k: perimeter.length(k) for k in by_kind}
    length_by_kind["CROSS_LAYER"] = perimeter.length("CROSS_LAYER")
    targeted = [e for e in perimeter.edges if e.kind == "PHYSICAL" and any(i[0] in targets for i in e.interfaces)]
    untargeted = [e for e in perimeter.edges if e.kind == "PHYSICAL" and not any(i[0] in targets for i in e.interfaces)]
    vertex_kinds = Counter(v.physical_kind for v in perimeter.vertices if v.physical_kind)
    turns = np.array([v.turn_degrees for v in perimeter.vertices if v.physical_kind in ("REGULAR", "CORNER") and v.turn_degrees is not None])
    histogram = np.histogram(turns, bins=TURN_BINS)[0].tolist() if turns.size else [0] * (len(TURN_BINS) - 1)
    corners = []
    excluded_vertices = []
    point_contacts = 0
    cuts = 0
    # The arc rule (decision 82(3)): corner vertices inside a fitted arc are absorbed (a
    # rounded corner or a bend describes them), not corners of their own.
    arcs = P.arc_groups(perimeter, radius)
    absorbed = {v for arc in arcs for v in arc["AbsorbedCorners"]}
    for index, v in enumerate(perimeter.vertices):
        if index in absorbed:
            continue
        if v.physical_kind in ("CORNER", "ENDPOINT", "JUNCTION"):
            if len({perimeter.edges[e].component for e in v.edges if perimeter.edges[e].kind not in P.CUT_KINDS}) > 1:
                point_contacts += 1
        if v.physical_kind == "ENDPOINT" and any(perimeter.edges[e].kind in P.CUT_KINDS for e in v.edges):
            # A physical chain ending on the truncation boundary or at a port face is a cut
            # (TruncationCut / PortCut), not a layout feature (decision 82(5)).
            cuts += 1
            continue
        if v.physical_kind in ("CORNER", "ENDPOINT", "JUNCTION"):
            # Every incident edge is an excluded class (embedded / non-planar metal) or the
            # vertex lies within 2R of off-plane metal: an excluded vertex of the manifest
            # (listed with the corners for the geometric comparison, not a feature vertex).
            excluded = index in perimeter.excluded_vertices or not any(perimeter.edges[e].kind == "PHYSICAL" for e in v.edges)
            if excluded:
                excluded_vertices.append([round(float(x), 9) for x in v.point])
            corners.append(
                {
                    "Point": [round(float(x), 9) for x in v.point],
                    "Kind": v.physical_kind,
                    "TurnDegrees": None if v.turn_degrees is None else round(v.turn_degrees, 6),
                    "InteriorAngleDegrees": None if v.turn_degrees is None else round(180.0 - v.turn_degrees, 6),
                    "Convex": v.convex,
                    "Degree": len([e for e in v.edges if perimeter.edges[e].kind == "PHYSICAL"]),
                    "Excluded": excluded,
                }
            )
    # Sub-threshold turns: vertices the classifier treats as straight although the layout
    # turns (polyline arcs); bend radius from the two incident segment lengths.
    smooth_turns = []
    for index, v in enumerate(perimeter.vertices):
        if v.physical_kind == "REGULAR" and v.turn_degrees is not None and v.turn_degrees > 1.0e-6:
            h = np.mean([perimeter.edges[e].length for e in v.edges if perimeter.edges[e].kind == "PHYSICAL"])
            bend_radius = h / (2.0 * math.sin(math.radians(v.turn_degrees) / 2.0))
            smooth_turns.append((v.turn_degrees, float(bend_radius)))
    smooth = np.array(smooth_turns) if smooth_turns else np.zeros((0, 2))
    runs = arcs
    interactions = P.edge_interactions(perimeter, radius)
    distances = np.array([i[2] for i in interactions]) if interactions else np.zeros(0)
    quantum = 1.0e-8 * radius
    knife = {
        "AtR": int(np.sum(np.abs(distances - radius) <= 1.0e-6 * radius)),
        "At2R": int(np.sum(np.abs(distances - 2.0 * radius) <= 1.0e-6 * radius)),
        "WithinHalfQuantumOfR": int(np.sum(np.abs(distances - radius) <= 0.5 * quantum)),
        "WithinHalfQuantumOf2R": int(np.sum(np.abs(distances - 2.0 * radius) <= 0.5 * quantum)),
    }
    parallel = [i for i in interactions if i[3] >= 1.0 - 1.0e-8]
    nonparallel = [i for i in interactions if i[3] < 1.0 - 1.0e-8]
    separations = Counter(round(i[2], 6) for i in parallel)
    return {
        "ProcessNormal": [round(float(x), 12) for x in perimeter.process_normal],
        "Planes": [round(p, 9) for p in perimeter.planes],
        "PlaneAreas": perimeter.plane_areas,
        "PrimaryPlane": perimeter.primary_plane,
        "NonplanarFaceArea": perimeter.nonplanar_face_area,
        "EdgesByClass": dict(by_kind),
        "LengthByClass": length_by_kind,
        "TargetedPhysicalEdges": len(targeted),
        "TargetedPhysicalLength": float(sum(e.length - e.cross_layer_length for e in targeted)),
        "TargetedCrossLayerLength": float(sum(e.cross_layer_length for e in targeted)),
        "MetalComponents": perimeter.components,
        "ExcludedVertices": len(excluded_vertices),
        "ExcludedVertexPoints": excluded_vertices[:50],
        "PointContacts": point_contacts,
        "UntargetedPhysicalEdges": len(untargeted),
        "UntargetedPhysicalLength": float(sum(e.length for e in untargeted)),
        "PhysicalChains": perimeter.chains,
        "VerticesByKind": dict(vertex_kinds),
        "TruncationCuts": cuts,
        "FeatureVertices": sum(1 for c in corners if not c["Excluded"]),
        "RoundedRuns": {"Runs": len(runs), "RoundedCorners": sum(1 for r in runs if r["Rounded"]), "Bends": sum(1 for r in runs if not r["Rounded"]), "AbsorbedCorners": len(absorbed), "Detail": runs[:50]},
        "TurnHistogram": {"Bins": TURN_BINS, "Counts": histogram},
        "SubThresholdTurns": {
            "Count": int(smooth.shape[0]),
            "MaxTurnDegrees": float(smooth[:, 0].max()) if smooth.size else 0.0,
            "MinBendRadius": float(smooth[:, 1].min()) if smooth.size else None,
            "MaxBendRadius": float(smooth[:, 1].max()) if smooth.size else None,
            "BendRadiusOverR": {
                "Min": float(smooth[:, 1].min() / radius) if smooth.size else None,
                "Median": float(np.median(smooth[:, 1]) / radius) if smooth.size else None,
            },
        },
        "Corners": corners,
        "Interactions": {
            "PairsWithin2R": len(interactions),
            "ParallelPairs": len(parallel),
            "NonparallelPairs": len(nonparallel),
            "ParallelSeparations": {str(k): v for k, v in sorted(separations.items())},
            "KnifeEdge": knife,
            "DistanceHistogram": np.histogram(distances, bins=[0, 0.5 * radius, radius * (1 - 1e-6), radius * (1 + 1e-6), 1.5 * radius, 2 * radius * (1 - 1e-6), 2 * radius * (1 + 1e-6), 2.01 * radius])[0].tolist() if distances.size else [],
        },
        "Conductors": {str(k): v for k, v in Counter(e.conductors for e in targeted).items()},
        "InterfaceSignatures": {str(k): v for k, v in Counter(e.interfaces for e in targeted).items()},
    }


def corner_clusters(perimeter, radius):
    """Single-linkage groups of non-regular physical vertices within 2R; centroid distances."""
    points = [v.point for v in perimeter.vertices if v.physical_kind in ("CORNER", "ENDPOINT", "JUNCTION")]
    n = len(points)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    # Candidate pairs from a k-d tree (a slightly larger radius), then the former exact test:
    # the single-linkage groups do not depend on the pair order.
    if n:
        tree = cKDTree(np.array(points))
        for i, j in sorted(tree.query_pairs(2.0 * radius * (1.0 + 1.0e-6))):
            if np.linalg.norm(points[i] - points[j]) <= 2.0 * radius * (1.0 + 1.0e-9):
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(n):
        groups[find(i)].append(i)
    clusters = []
    for members in groups.values():
        pts = np.array([points[m] for m in members])
        clusters.append({"Vertices": len(members), "Centroid": [round(float(x), 9) for x in pts.mean(axis=0)], "Diameter": float(max((np.linalg.norm(a - b) for a in pts for b in pts), default=0.0))})
    multi = [c for c in clusters if c["Vertices"] > 1]
    min_center_distance = math.inf
    if len(multi) > 1:
        # The closest pair of centroids: every centroid's nearest neighbour from a k-d tree
        # (the closest pair is among them), the value from the former norm.
        centroids = np.array([c["Centroid"] for c in multi])
        _, nearest = cKDTree(centroids).query(centroids, k=2)
        for i, j in enumerate(nearest[:, 1]):
            min_center_distance = min(min_center_distance, float(np.linalg.norm(np.array(multi[i]["Centroid"]) - np.array(multi[int(j)]["Centroid"]))))
    return {"Clusters": clusters, "MultiVertexClusters": len(multi), "MinMultiClusterCenterDistance": None if math.isinf(min_center_distance) else min_center_distance}


def _segment_key(p0, p1, decimals=7):
    a = tuple(round(float(x), decimals) for x in p0)
    b = tuple(round(float(x), decimals) for x in p1)
    return (a, b) if a <= b else (b, a)


# Audit edge kind -> manifest exclusion class (FOLD edges are NonPlanar exclusions).
EXCLUSION_CLASS_OF_AUDIT_KIND = {"NONPLANAR": "NonPlanar", "FOLD": "NonPlanar", "CROSS_LAYER": "CrossLayer", "NONMANIFOLD": "NonManifold", "EMBEDDED": "UndeterminedProcessSide", "BOX": "SimulationBoundary", "PORT": "Port"}
VERTEX_FEATURE_TYPES = ("ConvexCorner", "ConcaveCorner", "Endpoint", "Junction")


def identification_gates(identification, perimeter, census, radius, targets, compare=None):
    """Gates A1 / A2 from the version-2 contract; returns (gates, summary)."""
    gates = []
    tol = max(1.0e-9 * radius, 1.0e-12)
    segments = identification["Segments"]
    features = {f["Id"]: f for f in identification["Features"]}

    # Perimeter agreement by canonical segment key (mesh units on both sides).
    audit_keys = {}
    for i, e in enumerate(perimeter.edges):
        p0, p1 = perimeter.edge_points(e)
        audit_keys[_segment_key(p0, p1)] = i
    manifest_keys = {_segment_key(s["Key"][0], s["Key"][1]): j for j, s in enumerate(segments)}
    missing_in_manifest = [k for k in audit_keys if k not in manifest_keys]
    missing_in_audit = [k for k in manifest_keys if k not in audit_keys]
    # Keys that straddle a rounding boundary of the manifest's length grid: match the
    # leftovers by distance (1e-6 R on both endpoints).
    if missing_in_manifest and missing_in_audit:
        key_tolerance = 1.0e-6 * radius
        still_missing = []
        # Candidates by the first endpoint's grid cell (the former scan of every leftover
        # manifest key; the first match in list order wins, matched keys are consumed).
        audit_first = np.array([o[0] for o in missing_in_audit], dtype=float)
        leftover_grid = P.BoxGrid(audit_first, audit_first, 8.0 * radius)
        consumed = set()
        for k in missing_in_manifest:
            hit = None
            for index in leftover_grid.query(np.asarray(k[0], dtype=float), np.asarray(k[0], dtype=float), key_tolerance):
                index = int(index)
                if index in consumed:
                    continue
                o = missing_in_audit[index]
                if all(abs(a - b) <= key_tolerance for pa, pb in zip(k, o) for a, b in zip(pa, pb)):
                    hit = index
                    break
            if hit is None:
                still_missing.append(k)
            else:
                consumed.add(hit)
        missing_in_audit = [o for index, o in enumerate(missing_in_audit) if index not in consumed]
        missing_in_manifest = still_missing
    # The classifier works on Palace's mesh (cracking bisects the elements next to
    # under-resolved internal sheets; second-order edges are sampled at their mid-nodes), so
    # a manifest segment may be a piece of an audit edge: accept a leftover manifest segment
    # whose endpoints both lie on one audit edge, and a leftover audit edge covered by such
    # pieces; the perimeter lengths must then agree.
    key_tolerance = 1.0e-6 * radius
    covered_audit = set()
    split_pieces = 0
    still_missing_in_audit = []
    # Candidates: the leftover audit edges whose (tolerance-enlarged) boxes hold the first
    # endpoint, in list order (the former scan of every leftover audit edge; first match wins).
    if missing_in_manifest:
        edge_lower = np.array([np.minimum(o[0], o[1]) for o in missing_in_manifest], dtype=float)
        edge_upper = np.array([np.maximum(o[0], o[1]) for o in missing_in_manifest], dtype=float)
        edge_grid = P.BoxGrid(edge_lower, edge_upper, 8.0 * radius)
    for k in missing_in_audit:
        a, b = np.array(k[0]), np.array(k[1])
        hit = None
        for index in (edge_grid.query(a, a, 10.0 * key_tolerance) if missing_in_manifest else ()):
            o = missing_in_manifest[int(index)]
            p0, p1 = np.array(o[0]), np.array(o[1])
            d = p1 - p0
            ok = True
            for q in (a, b):
                t = float((q - p0) @ d) / float(d @ d)
                if t < -1.0e-9 or t > 1.0 + 1.0e-9 or np.linalg.norm(q - (p0 + t * d)) > key_tolerance:
                    ok = False
                    break
            if ok:
                hit = o
                break
        if hit is None:
            still_missing_in_audit.append(k)
        else:
            covered_audit.add(hit)
            split_pieces += 1
    missing_in_audit = still_missing_in_audit
    missing_in_manifest = [k for k in missing_in_manifest if k not in covered_audit]
    # Segments the classifier never sees are reported by kind, separately from PHYSICAL
    # disagreements.
    unseen_by_kind = Counter(perimeter.edges[audit_keys[k]].kind for k in missing_in_manifest)
    manifest_length = sum(float(s["Length"]) for s in segments)
    gates.append(
        gate(
            "perimeter-agreement",
            not missing_in_manifest and not missing_in_audit and abs(manifest_length - perimeter.length()) <= 1.0e-6 * max(perimeter.length(), radius),
            {
                "AuditSegments": len(perimeter.edges),
                "ManifestSegments": len(segments),
                "ManifestPerimeterLength": manifest_length,
                "AuditPerimeterLength": perimeter.length(),
                "SplitPieces": split_pieces,
                "UnseenByClassifier": dict(unseen_by_kind),
                "UnseenLength": float(sum(perimeter.edges[audit_keys[k]].length for k in missing_in_manifest)),
                "NotInMesh": len(missing_in_audit),
            },
        )
    )

    # Partition: every non-excluded segment covered exactly once by its portions.
    assigned_length = 0.0
    excluded_length = 0.0
    defects = []
    covered_segments = 0
    exclusion_records = identification["Exclusions"]
    for j, s in enumerate(segments):
        length = float(s["Length"])
        if "Exclusion" in s:
            excluded_length += length
            continue
        # Feature portions and analytically excluded portions (CrossLayer zones) tile the
        # segment together.
        portions = sorted((float(a), float(b), int(f)) for a, b, f in s.get("Portions", []))
        excluded_portions = sorted((float(a), float(b), -1 - int(x)) for a, b, x in s.get("ExcludedPortions", []))
        if not portions and not excluded_portions:
            defects.append({"Segment": j, "Defect": "no portion", "Length": length})
            continue
        covered = 0.0
        cursor = 0.0
        ok = True
        for a, b, f in sorted(portions + excluded_portions):
            if a < cursor - tol or b <= a:
                ok = False
            if f >= 0 and f not in features:
                ok = False
            if f < 0 and not (0 <= -1 - f < len(exclusion_records)):
                ok = False
            covered += b - a
            cursor = max(cursor, b)
        if not ok or abs(covered - length) > 1.0e-6 * max(length, radius) or abs(cursor - length) > 1.0e-6 * max(length, radius):
            defects.append({"Segment": j, "Defect": "gap or overlap", "Length": length, "Covered": covered, "Portions": (portions + excluded_portions)[:8]})
        assigned_length += sum(b - a for a, b, f in portions)
        excluded_length += sum(b - a for a, b, f in excluded_portions)
        covered_segments += 1
    targeted_length = census["TargetedPhysicalLength"]
    gates.append(
        gate(
            "A1-length-partition",
            abs(assigned_length + excluded_length - manifest_length) <= 1.0e-6 * max(manifest_length, radius) and abs(assigned_length - targeted_length) <= 1.0e-6 * max(targeted_length, radius),
            {"AssignedLength": assigned_length, "ExcludedLength": excluded_length, "ManifestPerimeterLength": manifest_length, "AuditTargetedPhysicalLength": targeted_length, "Deficit": targeted_length - assigned_length},
        )
    )
    unassigned = [j for j, s in enumerate(segments) if "Exclusion" not in s and not s.get("Portions") and not s.get("ExcludedPortions")]
    gates.append(gate("A1-count-partition", not unassigned and covered_segments + sum(1 for s in segments if "Exclusion" in s) == len(segments), {"Segments": len(segments), "Covered": covered_segments, "Excluded": sum(1 for s in segments if "Exclusion" in s), "Unassigned": len(unassigned)}))
    gates.append(gate("A1-multiplicity", not defects, {"Defects": len(defects), "Examples": defects[:10], "Basis": "per-segment portions sorted and contiguous, each portion in exactly one feature"}))

    # Feature portions must agree with the segment table (the two views of one assignment).
    feature_lengths = defaultdict(float)
    for s in segments:
        for a, b, f in s.get("Portions", []):
            feature_lengths[int(f)] += float(b) - float(a)
    length_mismatch = [f for f in features.values() if abs(feature_lengths.get(f["Id"], 0.0) - float(f["Length"])) > 1.0e-6 * max(float(f["Length"]), radius)]
    gates.append(gate("A1-weights", not length_mismatch, {"FeatureLengthMismatches": len(length_mismatch), "Basis": "every feature's Length equals the sum of its segment portions; one model per feature (weight 1)"}))

    # Vertex census: every corner / endpoint / junction (not a cut) and every rounded corner is
    # exactly one feature or cluster member.
    manifest_vertices = identification["Vertices"]
    typed = Counter(v["Type"] for v in manifest_vertices)
    unassigned_vertices = [v for v in manifest_vertices if v["Type"] in VERTEX_FEATURE_TYPES + ("RoundedCorner",) and v.get("Feature", -1) < 0]
    manifest_feature_vertices = sum(typed.get(t, 0) for t in VERTEX_FEATURE_TYPES)
    manifest_rounded = typed.get("RoundedCorner", 0)
    audit_feature_vertices = census["FeatureVertices"]
    audit_rounded = census["RoundedRuns"]["RoundedCorners"]
    manifest_angles = Counter()
    for v in manifest_vertices:
        if v["Type"] in ("ConvexCorner", "ConcaveCorner"):
            manifest_angles[round(180.0 - float(v["TurnDegrees"]), 6)] += 1
    audit_angles = Counter(round(180.0 - c["TurnDegrees"], 6) for c in census["Corners"] if c["Kind"] == "CORNER" and not c["Excluded"])
    # Excluded vertices (every incident run excluded, or within 2R of off-plane metal) are
    # accounted on both sides; point contacts are reported, never silent.
    manifest_point_contacts = sum(1 for v in manifest_vertices if v.get("PointContact"))
    gates.append(
        gate(
            "A1-vertex-census",
            not unassigned_vertices
            and manifest_feature_vertices == audit_feature_vertices
            and manifest_rounded == audit_rounded
            and manifest_angles == audit_angles
            and typed.get("Excluded", 0) == census["ExcludedVertices"]
            and manifest_point_contacts == census["PointContacts"],
            {
                "AuditCornerEndpointJunction": audit_feature_vertices,
                "ManifestCornerEndpointJunction": manifest_feature_vertices,
                "AuditRoundedCorners": audit_rounded,
                "ManifestRoundedCorners": manifest_rounded,
                "AuditExcludedVertices": census["ExcludedVertices"],
                "ManifestExcludedVertices": typed.get("Excluded", 0),
                "AuditPointContacts": census["PointContacts"],
                "ManifestPointContacts": manifest_point_contacts,
                "ManifestByType": dict(typed),
                "UnassignedVertices": len(unassigned_vertices),
                "AuditInteriorAngles": {str(k): v for k, v in sorted(audit_angles.items())},
                "ManifestCornerAngles": {str(k): v for k, v in sorted(manifest_angles.items())},
            },
        )
    )

    # Exclusions: every excluded class the audit finds in the mesh must be recorded with its
    # length (classes the classifier never sees are reported as the gap they are).
    recorded = defaultdict(float)
    for e in identification["Exclusions"]:
        recorded[e["Class"]] += float(e["Length"])
    exclusion_detail = {}
    exclusions_ok = True
    audit_by_class = defaultdict(float)
    for kind, cls in EXCLUSION_CLASS_OF_AUDIT_KIND.items():
        audit_by_class[cls] += census["LengthByClass"].get(kind, 0.0)
    for cls, audit_length in sorted(audit_by_class.items()):
        exclusion_detail[cls] = {"AuditLength": audit_length, "RecordedLength": recorded.get(cls, 0.0)}
        if abs(recorded.get(cls, 0.0) - audit_length) > 1.0e-6 * max(audit_length, radius):
            exclusions_ok = False
    exclusion_detail["Recorded"] = dict(recorded)
    exclusion_detail["UnseenByClassifier"] = dict(unseen_by_kind)
    gates.append(gate("A1-exclusions-recorded", exclusions_ok, exclusion_detail))

    # A2: cluster regions are disjoint by construction (merged when their cores are within 2R);
    # check the observable consequence: no two clusters share a segment portion, and report
    # the minimum distance between the frames' origins.
    clusters = [f for f in features.values() if f["Type"] == "SpatialEdgeCluster"]
    shared = 0
    owner = {}
    for f in clusters:
        for seg, a, b in f["Portions"]:
            for g, (c, d) in owner.get(seg, []):
                if min(b, d) - max(a, c) > tol and g != f["Id"]:
                    shared += 1
            owner.setdefault(seg, []).append((f["Id"], (a, b)))
    origins = [np.array(f["Frame"]["Origin"]) for f in clusters]
    min_origin_distance = min((float(np.linalg.norm(a - b)) for i, a in enumerate(origins) for b in origins[i + 1 :]), default=None)
    gates.append(gate("A2-cluster-balls", shared == 0, {"Clusters": len(clusters), "SharedPortions": shared, "MinOriginDistance": min_origin_distance, "TwoR": 2.0 * radius, "Basis": "cluster regions are the connected union of radius-R balls (merged below 2R); disjointness checked on the claimed portions"}))

    summary = {
        "Features": Counter(f["Type"] for f in features.values()),
        "FeatureLengths": {t: sum(float(f["Length"]) for f in features.values() if f["Type"] == t) for t in sorted({f["Type"] for f in features.values()})},
        "Matched": {f"{t}/{st}": n for (t, st), n in sorted(Counter((f["Type"], f["Match"]["Status"]) for f in features.values()).items())},
        "DistinctHashes": len({f["Hash"] for f in features.values()}),
        "Clusters": [{"Length": f["Length"], "EdgeCount": f["Signature"].get("EdgeCount"), "Vertices": len(f["Signature"].get("Vertices", [])), "Hash": f["Hash"][:12], "Chirality": f["Chirality"], "Status": f["Match"]["Status"]} for f in clusters],
        "Exclusions": identification["Exclusions"],
        "Totals": identification["Totals"],
        "GeometryDigest": identification["GeometryDigest"],
        "Conventions": identification.get("Conventions"),
    }
    return gates, summary


PATCH_LONGITUDINAL_TYPES = ("IsolatedEdge", "CurvedEdge", "SameConductorGap", "DifferentConductorGap", "SameConductorStrip", "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip", "ParallelEdgeCluster")
PATCH_SINGLE_TYPES = VERTEX_FEATURE_TYPES + ("SpatialEdgeCluster",)


def load_patches(path):
    """Rows of surface-response-patches.csv with numeric fields converted."""
    rows = []
    with open(path, newline="") as source:
        for row in csv.DictReader(source):
            entry = dict(row)
            for key in ("Patch", "Feature", "ModelIndex", "Segment"):
                entry[key] = int(row[key])
            for key in ("Weight", "ModelWeight", "QuadratureWeight", "SideFactor", "CouponDepth", "S0", "S1"):
                entry[key] = float(row[key])
            rows.append(entry)
    return rows


def patch_gates(identification, patches, radius):
    """Gates A7 on the patch dry run; returns (gates, summary)."""
    gates = []
    tol = 1.0e-6 * radius
    features = {f["Id"]: f for f in identification["Features"]}
    segments = identification["Segments"]
    matched = {f["Id"] for f in features.values() if f["Match"]["Status"] == "Matched"}
    rows_by_feature = defaultdict(list)
    for row in patches:
        rows_by_feature[row["Feature"]].append(row)
    patched = set(rows_by_feature)

    # A7-patch-features: the patched set is exactly the matched set.
    unknown = sorted(f for f in patched if f not in features)
    gates.append(
        gate(
            "A7-patch-features",
            patched == matched and not unknown,
            {"MatchedFeatures": len(matched), "PatchedFeatures": len(patched), "MatchedNotPatched": sorted(matched - patched)[:20], "PatchedNotMatched": sorted(patched - matched)[:20], "UnknownFeatureIds": unknown[:20], "Patches": len(patches)},
        )
    )

    # A7-patch-exclusions: no patch on an unmatched feature, an excluded segment or an
    # analytically excluded portion.
    violations = []
    for row in patches:
        if row["Feature"] not in matched:
            violations.append({"Patch": row["Patch"], "Defect": "unmatched feature", "Feature": row["Feature"]})
            continue
        if row["Segment"] >= 0:
            if row["Segment"] >= len(segments):
                violations.append({"Patch": row["Patch"], "Defect": "unknown segment", "Segment": row["Segment"]})
                continue
            segment = segments[row["Segment"]]
            if "Exclusion" in segment:
                violations.append({"Patch": row["Patch"], "Defect": "excluded segment", "Segment": row["Segment"], "Class": segment["Exclusion"]["Class"]})
            for a, b, _ in segment.get("ExcludedPortions", []):
                if min(row["S1"], float(b)) - max(row["S0"], float(a)) > tol:
                    violations.append({"Patch": row["Patch"], "Defect": "excluded portion", "Segment": row["Segment"]})
    gates.append(gate("A7-patch-exclusions", not violations, {"Violations": len(violations), "Examples": violations[:10]}))

    # A7-patch-coverage: the distinct intervals of a longitudinal feature are exactly its
    # manifest portions (each portion once); a vertex / cluster feature carries patches
    # without a portion.
    coverage_defects = []
    covered_length = 0.0
    weight_defects = []
    for feature_id in sorted(matched):
        feature = features[feature_id]
        rows = rows_by_feature.get(feature_id, [])
        if not rows:
            continue  # reported by A7-patch-features
        portions = sorted((int(seg), float(a), float(b)) for seg, a, b in feature["Portions"])
        groups = defaultdict(list)
        for row in rows:
            groups[(row["Segment"], round(row["S0"], 9), round(row["S1"], 9))].append(row)
        if feature["Type"] in PATCH_LONGITUDINAL_TYPES:
            intervals = sorted((seg, a, b) for (seg, a, b) in groups)
            if len(intervals) != len(portions) or any(i[0] != p[0] or abs(i[1] - p[1]) > tol or abs(i[2] - p[2]) > tol for i, p in zip(intervals, portions)):
                coverage_defects.append({"Feature": feature_id, "Type": feature["Type"], "Portions": len(portions), "PatchIntervals": len(intervals), "Examples": [(i, p) for i, p in zip(intervals, portions) if i[0] != p[0] or abs(i[1] - p[1]) > tol or abs(i[2] - p[2]) > tol][:4]})
            else:
                covered_length += sum(b - a for _, a, b in intervals)
            # Side factor: 1 / number of sides of the feature (pairs 1/2, clusters 1/n;
            # the manifest's Sides labels, both sides of a pair may share one chain).
            sides = set(feature.get("Sides", [])) or {segments[seg]["Chain"] for seg, _, _ in portions}
            expected_side = 1.0 / len(sides) if feature["Type"] not in ("IsolatedEdge", "CurvedEdge") else 1.0
            for key, group in groups.items():
                quadrature = sum(r["QuadratureWeight"] * r["ModelWeight"] for r in group)
                if abs(quadrature - 1.0) > 1.0e-9:
                    weight_defects.append({"Feature": feature_id, "Interval": key, "Defect": "quadrature x model weights do not sum to 1", "Sum": quadrature})
                for r in group:
                    expected = r["ModelWeight"] * r["QuadratureWeight"] * (r["S1"] - r["S0"]) * r["SideFactor"] / r["CouponDepth"] if r["CouponDepth"] > 0 else float("nan")
                    # S0 / S1 / CouponDepth are written on the manifest's length grid (<= 1e-10 R):
                    # the formula is checked to that grid on the portion length.
                    grid = r["ModelWeight"] * r["QuadratureWeight"] * r["SideFactor"] * 1.0e-9 * radius / r["CouponDepth"] if r["CouponDepth"] > 0 else 0.0
                    if not (r["CouponDepth"] > 0) or abs(r["Weight"] - expected) > 1.0e-9 * abs(expected) + grid:
                        weight_defects.append({"Feature": feature_id, "Patch": r["Patch"], "Defect": "weight formula", "Weight": r["Weight"], "Expected": expected})
                    if abs(r["SideFactor"] - expected_side) > 1.0e-12:
                        weight_defects.append({"Feature": feature_id, "Patch": r["Patch"], "Defect": "side factor", "SideFactor": r["SideFactor"], "Expected": expected_side})
        elif feature["Type"] in PATCH_SINGLE_TYPES:
            if any(r["Segment"] >= 0 for r in rows):
                coverage_defects.append({"Feature": feature_id, "Type": feature["Type"], "Defect": "vertex / cluster patch with a portion"})
            else:
                covered_length += sum(b - a for _, a, b in portions)
            model_weight = sum(r["ModelWeight"] for r in rows)
            if abs(model_weight - 1.0) > 1.0e-9:
                weight_defects.append({"Feature": feature_id, "Defect": "model weights do not sum to 1", "Sum": model_weight})
            for r in rows:
                if abs(r["Weight"] - r["ModelWeight"]) > 1.0e-12 or r["CouponDepth"] != 0.0:
                    weight_defects.append({"Feature": feature_id, "Patch": r["Patch"], "Defect": "vertex weight is not the model weight", "Weight": r["Weight"], "ModelWeight": r["ModelWeight"]})
        else:
            coverage_defects.append({"Feature": feature_id, "Type": feature["Type"], "Defect": "type without a patch construction"})
    matched_length = sum(float(features[f]["Length"]) for f in matched)
    assigned_length = float(identification["Totals"]["AssignedLength"])
    gates.append(
        gate(
            "A7-patch-coverage",
            not coverage_defects and abs(covered_length - matched_length) <= 1.0e-6 * max(matched_length, radius),
            {"Defects": len(coverage_defects), "Examples": coverage_defects[:10], "CoveredLength": covered_length, "MatchedLength": matched_length, "AssignedLength": assigned_length, "CoveredFractionOfAssigned": covered_length / assigned_length if assigned_length else None, "Basis": "every portion of a matched longitudinal feature is one quadrature interval; a vertex / cluster feature is one patch"},
        )
    )
    gates.append(gate("A7-patch-weights", not weight_defects, {"Defects": len(weight_defects), "Examples": weight_defects[:10], "Basis": "per interval sum(quadrature x model weight) = 1; weight = model x quadrature x length x side factor / coupon depth; side factor = 1 / chains of a pair or parallel cluster"}))
    summary = {
        "Patches": len(patches),
        "PatchesByTopology": dict(Counter(r["Topology"] for r in patches)),
        "PatchedFeatures": len(patched),
        "MatchedFeatures": len(matched),
        "CoveredLength": covered_length,
        "MatchedLength": matched_length,
        "AssignedLength": assigned_length,
        "CoveredFractionOfAssigned": covered_length / assigned_length if assigned_length else None,
        "Models": sorted({r["Model"] for r in patches}),
    }
    return gates, summary


def progress(message, started=[None]):
    """Section line with the wall time since the first call (stderr, flushed): a chip-scale
    audit runs for tens of minutes."""
    import time
    now = time.time()
    if started[0] is None:
        started[0] = now
    print(f"[audit {now - started[0]:8.1f} s] {message}", file=sys.stderr, flush=True)


def run_audit(args):
    with open(args.config) as source:
        config = json.load(source)
    progress(f"reading mesh {args.mesh}")
    mesh = read_msh2(args.mesh)
    progress(f"mesh read: {len(mesh.coordinates)} nodes; loading manifest")
    manifest = M.load_manifest(args.manifest)
    summary = M.summarize(manifest)
    radius = args.radius or summary["MatchingRadius"]
    if not radius:
        raise SystemExit("matching radius unknown: pass --radius")
    process_normal = None
    for entry in config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", []):
        if "EdgeFrameNormal" in entry:
            process_normal = entry["EdgeFrameNormal"]
    progress("extracting the perimeter from the mesh")
    perimeter = P.extract_perimeter(mesh, config, process_normal=process_normal, corner_tolerance_degrees=args.corner_tolerance, radius=radius)
    progress(f"perimeter: {len(perimeter.edges)} edges, {len(perimeter.vertices)} vertices, {perimeter.chains} chains; census (edge interactions, rounded runs)")
    targets = set(P.target_interfaces(config)) or {i for i, _ in P.interface_attributes(config)}
    census = perimeter_census(perimeter, radius, targets)
    progress("corner clusters")
    clusters = corner_clusters(perimeter, radius)
    progress(f"{len(clusters['Clusters'])} corner clusters; gates")
    log = M.parse_palace_log(args.log) if args.log else None

    gates = []
    statistics = summary["Statistics"]
    identification = manifest.get("Identification")
    v2_summary = None
    patch_summary = None
    if identification:
        gates, v2_summary = identification_gates(identification, perimeter, census, radius, targets)
        patches_path = getattr(args, "patches", None) or os.path.join(os.path.dirname(os.path.abspath(args.manifest)), "surface-response-patches.csv")
        if os.path.exists(patches_path):
            patch_gate_list, patch_summary = patch_gates(identification, load_patches(patches_path), radius)
            gates.extend(patch_gate_list)
            patch_summary["Path"] = os.path.abspath(patches_path)
    progress(f"{len(gates)} gates evaluated; summary")
    audit_segments = census["EdgesByClass"].get("PHYSICAL", 0) + census["EdgesByClass"].get("TRUNCATION", 0)
    if "MetalSegments" in statistics and not identification:
        bisected = bool(log and log.get("Bisection"))
        gates.append(
            gate(
                "perimeter-agreement",
                statistics["MetalSegments"] == audit_segments or bisected,
                {"ClassifierMetalSegments": statistics["MetalSegments"], "AuditPhysicalPlusTruncation": audit_segments, "LocalBisection": log.get("Bisection") if log else "unknown (no log)"},
            )
        )
        if "PhysicalChains" in statistics:
            gates.append(gate("chain-agreement", statistics["PhysicalChains"] == perimeter.chains, {"ClassifierPhysicalChains": statistics["PhysicalChains"], "AuditChains": perimeter.chains}))

    length = census["TargetedPhysicalLength"]
    if identification:
        deficit = length - float(identification["Totals"]["AssignedLength"])
        excluded_length = sum(census["LengthByClass"].get(k, 0.0) for k in ("NONPLANAR", "FOLD", "CROSS_LAYER", "NONMANIFOLD", "EMBEDDED", "BOX", "PORT"))
    else:
        assigned = summary["TranslationalLength"]
        deficit = length - assigned
        tolerance = max(1.0e-10 * radius * max(1, census["TargetedPhysicalEdges"]), 1.0e-9 * length)
        gates.append(
            gate(
                "A1-length-partition",
                abs(deficit) <= tolerance,
                {"TargetedPhysicalLength": length, "AssignedTranslationalLength": assigned, "Deficit": deficit, "DeficitFraction": deficit / length if length else None, "Tolerance": tolerance},
            )
        )
        omitted = log["OmittedSegments"] if log else None
        count_detail = {"TargetedPhysicalEdges": census["TargetedPhysicalEdges"], "AssignedTranslationalCount": summary["TranslationalCount"], "OmittedFromLog": omitted}
        if omitted is None:
            count_ok = summary["TranslationalCount"] == census["TargetedPhysicalEdges"]
        else:
            count_ok = summary["TranslationalCount"] + omitted == census["TargetedPhysicalEdges"] and omitted == 0
            count_detail["Reconciled"] = summary["TranslationalCount"] + omitted == census["TargetedPhysicalEdges"]
        gates.append(gate("A1-count-partition", count_ok, count_detail))
        gates.append(
            gate(
                "A1-multiplicity",
                count_ok and abs(deficit) <= tolerance,
                {"Basis": "aggregate proxy (count and length balance); the manifest carries no per-segment assignment, so an equal-length double count and gap cannot be separated", "ContractGap": "no per-segment assignment / feature positions in surface-response-requirements.json"},
                evaluable=False if not (count_ok and abs(deficit) <= tolerance) else True,
            )
        )
        gates.append(gate("A1-weights", not summary["WeightDefects"], {"Defects": summary["WeightDefects"]}))

        audit_vertices = census["FeatureVertices"] + census["RoundedRuns"]["RoundedCorners"]
        manifest_vertices = summary["VertexFeatureCount"]
        # The log's unmatched vertices are the Missing vertex records of the manifest (no model),
        # not an extra class; the residual is what the manifest does not enumerate at all.
        unmatched = log.get("UnmatchedVertices", 0) if log else 0
        residual = audit_vertices - manifest_vertices
        corner_angles = Counter(round(180.0 - c["TurnDegrees"], 6) for c in census["Corners"] if c["Kind"] == "CORNER" and not c["Excluded"])
        manifest_angles = Counter()
        for r in manifest["Requirements"]:
            if r["Topology"] in ("ConvexCorner", "ConcaveCorner"):
                manifest_angles[round(float(r["Geometry"].get("AngleDegrees", float("nan"))), 6)] += int(r["Count"])
        # A positive residual with spatial clusters present may be corners absorbed into the
        # clusters (the manifest does not enumerate them: contract gap) or silently dropped
        # vertices; without clusters it is a silent drop; a negative residual is a double count.
        cluster_records = len(summary["Clusters"])
        gates.append(
            gate(
                "A1-vertex-census",
                residual == 0,
                {
                    "AuditCornerEndpointJunction": audit_vertices,
                    "AuditByKind": {k: census["VerticesByKind"].get(k, 0) for k in ("CORNER", "ENDPOINT", "JUNCTION")},
                    "AuditTruncationCuts": census["TruncationCuts"],
                    "AuditRoundedCorners": census["RoundedRuns"]["RoundedCorners"],
                    "ManifestVertexFeatures": manifest_vertices,
                    "UnmatchedVerticesFromLog": unmatched,
                    "Residual": residual,
                    "ManifestClusterRecords": cluster_records,
                    "ResidualMeaning": "audit vertices minus manifest vertex records: > 0 with clusters = absorbed into clusters or dropped (not enumerated by the contract), > 0 without clusters = dropped, < 0 = double counted",
                    "AuditInteriorAngles": {str(k): v for k, v in sorted(corner_angles.items())},
                    "ManifestCornerAngles": {str(k): v for k, v in sorted(manifest_angles.items())},
                    "CornerTurnToleranceDegrees": args.corner_tolerance,
                },
                evaluable=not (residual > 0 and cluster_records > 0),
            )
        )
        excluded_length = sum(census["LengthByClass"].get(k, 0.0) for k in ("NONPLANAR", "FOLD", "CROSS_LAYER", "NONMANIFOLD", "EMBEDDED", "BOX", "PORT"))
        gates.append(
            gate(
                "A1-exclusions-recorded",
                excluded_length == 0.0,
                {
                    "NonplanarLength": census["LengthByClass"].get("NONPLANAR", 0.0),
                    "CrossLayerLength": census["LengthByClass"].get("CROSS_LAYER", 0.0),
                    "NonmanifoldLength": census["LengthByClass"].get("NONMANIFOLD", 0.0),
                    "UntargetedPhysicalLength": census["UntargetedPhysicalLength"],
                    "Meaning": "decision 73(3): excluded classes must be reported by the identification as coverage gaps; the manifest has no such record, so any nonzero excluded length is a silent omission",
                },
            )
        )
        a2_ok = clusters["MinMultiClusterCenterDistance"] is None or clusters["MinMultiClusterCenterDistance"] >= 2.0 * radius
        gates.append(
            gate(
                "A2-cluster-balls",
                a2_ok,
                {
                    "Basis": "audit-side single-linkage vertex clusters (<= 2R); the manifest carries no cluster positions",
                    "MultiVertexClusters": clusters["MultiVertexClusters"],
                    "ManifestClusterRecords": len(summary["Clusters"]),
                    "MinMultiClusterCenterDistance": clusters["MinMultiClusterCenterDistance"],
                    "TwoR": 2.0 * radius,
                },
            )
        )

    cluster_interval_length = 0.0
    for r in manifest["Requirements"]:
        if r["Topology"] in M.CLUSTER:
            for e in r["Geometry"].get("Edges", []):
                interval = e.get("Interval") if isinstance(e, dict) else None
                if interval and len(interval) == 2:
                    cluster_interval_length += (float(interval[1]) - float(interval[0])) * int(r["Count"])
    gap = {
        "Unit": manifest.get("LengthUnit", "mesh"),
        "ClusterIntervalLength": cluster_interval_length,
        "ClusterIntervalMeaning": "sum of the cluster records' edge intervals (the coupon description clipped to the matching ball), not an assignment of perimeter segments",
        "TotalTargetedLength": length,
        "OmittedLength": deficit,
        "OmittedFraction": deficit / length if length else None,
        "MissingModelLength": sum(v["MissingLength"] for k, v in summary["ByTopology"].items() if k in M.TRANSLATIONAL),
        "MissingModelFraction": (sum(v["MissingLength"] for k, v in summary["ByTopology"].items() if k in M.TRANSLATIONAL) / length) if length else None,
        "ExcludedLength": excluded_length,
        "UntargetedLength": census["UntargetedPhysicalLength"],
        "EnergyBound": "not available: the isolated-edge per-length defect needs the library model matrices (CSV), absent from the local library JSON",
        "OmittedByClassFromLog": log,
    }

    result = {
        "Inputs": {"Mesh": os.path.abspath(args.mesh), "Config": os.path.abspath(args.config), "Manifest": os.path.abspath(args.manifest), "Log": os.path.abspath(args.log) if args.log else None, "L0": config.get("Model", {}).get("L0")},
        "MatchingRadius": radius,
        "Perimeter": census,
        "CornerClusters": {k: v for k, v in clusters.items() if k != "Clusters"},
        "Manifest": {k: v for k, v in summary.items() if k != "Statistics"},
        "ClassifierStatistics": statistics,
        "Digest": {"Full": M.canonical_digest(manifest), "GeometryOnly": M.canonical_digest(manifest, geometry_only_counts=True), "Geometry": identification["GeometryDigest"] if identification else None},
        "Identification": v2_summary,
        "Patches": patch_summary,
        "Gates": gates,
        "GapBound": gap,
    }
    if args.compare:
        other = M.load_manifest(args.compare)
        result["Compare"] = {
            "Other": os.path.abspath(args.compare),
            "OtherDigest": {"Full": M.canonical_digest(other), "GeometryOnly": M.canonical_digest(other, geometry_only_counts=True)},
            "Diff": M.diff_manifests(manifest, other),
            "DiffGeometryOnly": M.diff_manifests(manifest, other, geometry_only_counts=True),
        }
        if identification and other.get("Identification"):
            same = identification["GeometryDigest"] == other["Identification"]["GeometryDigest"]
            mine = Counter((f["Type"], f["Hash"]) for f in identification["Features"])
            theirs = Counter((f["Type"], f["Hash"]) for f in other["Identification"]["Features"])
            # Lengths per signature agree to a tolerance (they are not hashed: roundoff).
            lengths_here, lengths_there = defaultdict(list), defaultdict(list)
            for f in identification["Features"]:
                lengths_here[f["Hash"]].append(float(f["Length"]))
            for f in other["Identification"]["Features"]:
                lengths_there[f["Hash"]].append(float(f["Length"]))
            length_mismatch = []
            for h in set(lengths_here) | set(lengths_there):
                a, b = sorted(lengths_here.get(h, [])), sorted(lengths_there.get(h, []))
                if len(a) != len(b) or any(abs(x - y) > 1.0e-6 * max(abs(x), radius) for x, y in zip(a, b)):
                    length_mismatch.append({"Hash": h[:12], "Here": a[:4], "There": b[:4]})
            result["Compare"]["GeometryDigestIdentical"] = same
            result["Compare"]["FeatureDiff"] = {"OnlyHere": sorted(f"{t}:{h[:12]}x{n}" for (t, h), n in (mine - theirs).items()), "OnlyThere": sorted(f"{t}:{h[:12]}x{n}" for (t, h), n in (theirs - mine).items()), "LengthMismatch": length_mismatch}
            gates.append(gate("A3/A5-set-identity", same and not length_mismatch, {"GeometryDigestIdentical": same, "FeatureDiff": result["Compare"]["FeatureDiff"], "V1Identical": result["Compare"]["Diff"]["Identical"]}))
        else:
            gates.append(gate("A3/A5-set-identity", result["Compare"]["Diff"]["Identical"], {"Added": len(result["Compare"]["Diff"]["Added"]), "Removed": len(result["Compare"]["Diff"]["Removed"]), "Changed": len(result["Compare"]["Diff"]["Changed"]), "GeometryOnlyIdentical": result["Compare"]["DiffGeometryOnly"]["Identical"]}))
    result["Passed"] = all(g["Status"] == "PASS" for g in gates)
    return result


def render_markdown(result):
    lines = [f"# Identification audit: {os.path.basename(result['Inputs']['Manifest'])}", ""]
    lines.append(f"Mesh `{result['Inputs']['Mesh']}`; R = {result['MatchingRadius']} (mesh units, L0 = {result['Inputs']['L0']}); overall **{'PASS' if result['Passed'] else 'FAIL'}**.")
    lines.append("")
    lines.append("| Gate | Status | Detail |")
    lines.append("|---|---|---|")
    for g in result["Gates"]:
        detail = json.dumps({k: v for k, v in g["Detail"].items() if k not in ("AuditInteriorAngles", "ManifestCornerAngles", "Meaning", "Basis", "ResidualMeaning", "ContractGap")}, default=str)
        lines.append(f"| {g['Gate']} | {g['Status']} | {detail} |")
    lines.append("")
    census = result["Perimeter"]
    lines.append("## Perimeter E (audit)")
    lines.append(f"- process normal {census['ProcessNormal']}, planes {census['Planes']} (areas {[round(a, 3) for a in census['PlaneAreas']]}), primary {census['PrimaryPlane']}")
    lines.append(f"- edges by class {census['EdgesByClass']}; length by class {{{', '.join(f'{k}: {v:.6f}' for k, v in census['LengthByClass'].items())}}}")
    lines.append(f"- targeted physical: {census['TargetedPhysicalEdges']} edges, L = {census['TargetedPhysicalLength']:.6f}; untargeted physical {census['UntargetedPhysicalEdges']} edges / {census['UntargetedPhysicalLength']:.6f}")
    lines.append(f"- chains {census['PhysicalChains']}; vertices {census['VerticesByKind']}")
    lines.append(f"- turn histogram (deg bins {census['TurnHistogram']['Bins']}): {census['TurnHistogram']['Counts']}")
    st = census["SubThresholdTurns"]
    lines.append(f"- sub-threshold turns (treated as straight by the 30 deg rule): {st['Count']}, max turn {st['MaxTurnDegrees']:.4f} deg, bend radius {st['MinBendRadius']} .. {st['MaxBendRadius']} (R_bend/R min {st['BendRadiusOverR']['Min']}, median {st['BendRadiusOverR']['Median']})")
    it = census["Interactions"]
    lines.append(f"- edge pairs within 2R: {it['PairsWithin2R']} (parallel {it['ParallelPairs']}, nonparallel {it['NonparallelPairs']}); parallel separations {it['ParallelSeparations']}; knife-edge {it['KnifeEdge']}")
    lines.append(f"- conductors on targeted edges {census['Conductors']}; interface signatures {census['InterfaceSignatures']}")
    lines.append("")
    if result.get("Identification"):
        ident = result["Identification"]
        lines.append("## Identification (manifest version 2)")
        lines.append(f"- features {dict(ident['Features'])}; lengths {{{', '.join(f'{k}: {v:.6f}' for k, v in ident['FeatureLengths'].items())}}}")
        lines.append(f"- matched {ident['Matched']}; distinct hashes {ident['DistinctHashes']}")
        lines.append(f"- totals {ident['Totals']}; exclusions {[(e['Class'], e['Count'], e['Length']) for e in ident['Exclusions']]}")
        lines.append(f"- clusters {ident['Clusters']}")
        lines.append(f"- geometry digest `{ident['GeometryDigest'][:16]}`")
        lines.append("")
    if result.get("Patches"):
        pt = result["Patches"]
        lines.append("## Patch dry run (surface-response-patches.csv)")
        lines.append(f"- patches {pt['Patches']} by topology {pt['PatchesByTopology']}; patched features {pt['PatchedFeatures']} of {pt['MatchedFeatures']} matched")
        lines.append(f"- covered length {pt['CoveredLength']:.6f} of assigned {pt['AssignedLength']:.6f} ({100.0 * (pt['CoveredFractionOfAssigned'] or 0):.4f} %); models {pt['Models'][:12]}")
        lines.append("")
    lines.append("## Manifest")
    for k, v in result["Manifest"]["ByTopology"].items():
        lines.append(f"- {k}: records {v['Records']}, count {v['Count']}, length {v['TotalEdgeLength']:.6f}, missing {v['Missing']} ({v['MissingLength']:.6f})")
    lines.append(f"- clusters: {result['Manifest']['Clusters']}")
    lines.append(f"- digest full `{result['Digest']['Full'][:16]}` geometry-only `{result['Digest']['GeometryOnly'][:16]}`; classifier statistics {result['ClassifierStatistics']}")
    lines.append("")
    gap = result["GapBound"]
    lines.append("## Gap bound (B1, length only)")
    lines.append(f"- omitted (unassigned) length {gap['OmittedLength']:.6f} = {100.0 * (gap['OmittedFraction'] or 0):.4f} % of L; missing-model length {gap['MissingModelLength']:.6f} = {100.0 * (gap['MissingModelFraction'] or 0):.4f} %; excluded classes {gap['ExcludedLength']:.6f}; untargeted {gap['UntargetedLength']:.6f}")
    lines.append(f"- {gap['EnergyBound']}")
    if gap["OmittedByClassFromLog"]:
        lines.append(f"- classifier log: {json.dumps(gap['OmittedByClassFromLog'])}")
    if "Compare" in result:
        lines.append("")
        lines.append("## Comparison")
        c = result["Compare"]
        lines.append(f"- other `{c['Other']}` digest full `{c['OtherDigest']['Full'][:16]}` geometry-only `{c['OtherDigest']['GeometryOnly'][:16]}`")
        lines.append(f"- diff: added {len(c['Diff']['Added'])}, removed {len(c['Diff']['Removed'])}, changed {len(c['Diff']['Changed'])}; geometry-only identical {c['DiffGeometryOnly']['Identical']}")
        for label in ("Added", "Removed", "Changed"):
            for entry in c["Diff"][label][:20]:
                lines.append(f"  - {label}: {json.dumps(entry, default=str)[:400]}")
    lines.append("")
    lines.append("## Corner / endpoint / junction vertices (audit)")
    for c in result["Perimeter"]["Corners"][:200]:
        lines.append(f"- {c['Kind']} at {c['Point']} turn {c['TurnDegrees']} interior {c['InteriorAngleDegrees']} convex {c['Convex']} degree {c['Degree']}")
    if len(result["Perimeter"]["Corners"]) > 200:
        lines.append(f"- ... {len(result['Perimeter']['Corners']) - 200} more in the JSON")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--config", required=True, help="the Palace configuration the preflight ran with")
    parser.add_argument("--manifest", required=True, help="surface-response-requirements.json")
    parser.add_argument("--log", help="palace stdout of the preflight (omission counts)")
    parser.add_argument("--compare", help="a second manifest to diff against (A3 / A5)")
    parser.add_argument("--patches", help="surface-response-patches.csv of the patch dry run (default: next to the manifest)")
    parser.add_argument("--radius", type=float, help="matching radius in mesh units (default: manifest)")
    parser.add_argument("--corner-tolerance", type=float, default=P.CORNER_ANGLE_TOLERANCE_DEGREES, help="turn (deg) above which a vertex is a corner (classifier: 30)")
    parser.add_argument("--output-prefix", help="write <prefix>.json and <prefix>.md")
    args = parser.parse_args(argv)
    result = run_audit(args)
    text = render_markdown(result)
    if args.output_prefix:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_prefix)), exist_ok=True)
        with open(args.output_prefix + ".json", "w") as out:
            json.dump(result, out, indent=1, default=str)
        with open(args.output_prefix + ".md", "w") as out:
            out.write(text)
    print(text if not args.output_prefix else "\n".join(text.splitlines()[:len(result["Gates"]) + 6]))
    return 0 if result["Passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
