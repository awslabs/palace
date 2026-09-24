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
"""

import argparse
import json
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

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
    targeted = [e for e in perimeter.edges if e.kind == "PHYSICAL" and any(i[0] in targets for i in e.interfaces)]
    untargeted = [e for e in perimeter.edges if e.kind == "PHYSICAL" and not any(i[0] in targets for i in e.interfaces)]
    vertex_kinds = Counter(v.physical_kind for v in perimeter.vertices if v.physical_kind)
    turns = np.array([v.turn_degrees for v in perimeter.vertices if v.physical_kind in ("REGULAR", "CORNER") and v.turn_degrees is not None])
    histogram = np.histogram(turns, bins=TURN_BINS)[0].tolist() if turns.size else [0] * (len(TURN_BINS) - 1)
    corners = []
    cuts = 0
    for index, v in enumerate(perimeter.vertices):
        if v.physical_kind == "ENDPOINT" and any(perimeter.edges[e].kind == "TRUNCATION" for e in v.edges):
            # A physical chain ending on the truncation boundary is a simulation cut, not a
            # layout feature (the classifier makes no Endpoint feature there).
            cuts += 1
            continue
        if v.physical_kind in ("CORNER", "ENDPOINT", "JUNCTION"):
            corners.append(
                {
                    "Point": [round(float(x), 9) for x in v.point],
                    "Kind": v.physical_kind,
                    "TurnDegrees": None if v.turn_degrees is None else round(v.turn_degrees, 6),
                    "InteriorAngleDegrees": None if v.turn_degrees is None else round(180.0 - v.turn_degrees, 6),
                    "Convex": v.convex,
                    "Degree": len([e for e in v.edges if perimeter.edges[e].kind == "PHYSICAL"]),
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
    runs = P.rounded_runs(perimeter, radius)
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
        "TargetedPhysicalLength": float(sum(e.length for e in targeted)),
        "UntargetedPhysicalEdges": len(untargeted),
        "UntargetedPhysicalLength": float(sum(e.length for e in untargeted)),
        "PhysicalChains": perimeter.chains,
        "VerticesByKind": dict(vertex_kinds),
        "TruncationCuts": cuts,
        "FeatureVertices": len(corners),
        "RoundedRuns": {"Runs": len(runs), "RoundedCorners": sum(1 for r in runs if r["Rounded"]), "Detail": runs[:50]},
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

    for i in range(n):
        for j in range(i + 1, n):
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
    for i in range(len(multi)):
        for j in range(i + 1, len(multi)):
            min_center_distance = min(min_center_distance, float(np.linalg.norm(np.array(multi[i]["Centroid"]) - np.array(multi[j]["Centroid"]))))
    return {"Clusters": clusters, "MultiVertexClusters": len(multi), "MinMultiClusterCenterDistance": None if math.isinf(min_center_distance) else min_center_distance}


def _segment_key(p0, p1, decimals=7):
    a = tuple(round(float(x), decimals) for x in p0)
    b = tuple(round(float(x), decimals) for x in p1)
    return (a, b) if a <= b else (b, a)


EXCLUSION_CLASS_OF_AUDIT_KIND = {"NONPLANAR": "NonPlanar", "CROSS_LAYER": "CrossLayer", "NONMANIFOLD": "NonManifold"}
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
        for k in missing_in_manifest:
            hit = next((o for o in missing_in_audit if all(abs(a - b) <= key_tolerance for pa, pb in zip(k, o) for a, b in zip(pa, pb))), None)
            if hit is None:
                still_missing.append(k)
            else:
                missing_in_audit.remove(hit)
        missing_in_manifest = still_missing
    # Segments the classifier never sees (sheets whose faces cancel in the odd-incidence
    # perimeter) are the audit's NONPLANAR / CROSS_LAYER / NONMANIFOLD edges in most meshes;
    # report them separately from PHYSICAL disagreements.
    unseen_by_kind = Counter(perimeter.edges[audit_keys[k]].kind for k in missing_in_manifest)
    manifest_length = sum(float(s["Length"]) for s in segments)
    gates.append(
        gate(
            "perimeter-agreement",
            unseen_by_kind.get("PHYSICAL", 0) == 0 and unseen_by_kind.get("TRUNCATION", 0) == 0 and not missing_in_audit,
            {
                "AuditSegments": len(perimeter.edges),
                "ManifestSegments": len(segments),
                "ManifestPerimeterLength": manifest_length,
                "AuditPerimeterLength": perimeter.length(),
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
    for j, s in enumerate(segments):
        length = float(s["Length"])
        if "Exclusion" in s:
            excluded_length += length
            continue
        portions = sorted((float(a), float(b), int(f)) for a, b, f in s.get("Portions", []))
        if not portions:
            defects.append({"Segment": j, "Defect": "no portion", "Length": length})
            continue
        covered = 0.0
        cursor = 0.0
        ok = True
        for a, b, f in portions:
            if a < cursor - tol or b <= a:
                ok = False
            if f not in features:
                ok = False
            covered += b - a
            cursor = max(cursor, b)
        if not ok or abs(covered - length) > 1.0e-6 * max(length, radius) or abs(cursor - length) > 1.0e-6 * max(length, radius):
            defects.append({"Segment": j, "Defect": "gap or overlap", "Length": length, "Covered": covered, "Portions": portions[:8]})
        assigned_length += covered
        covered_segments += 1
    targeted_length = census["TargetedPhysicalLength"]
    gates.append(
        gate(
            "A1-length-partition",
            abs(assigned_length + excluded_length - manifest_length) <= 1.0e-6 * max(manifest_length, radius) and abs(assigned_length - targeted_length) <= 1.0e-6 * max(targeted_length, radius),
            {"AssignedLength": assigned_length, "ExcludedLength": excluded_length, "ManifestPerimeterLength": manifest_length, "AuditTargetedPhysicalLength": targeted_length, "Deficit": targeted_length - assigned_length},
        )
    )
    unassigned = [j for j, s in enumerate(segments) if "Exclusion" not in s and not s.get("Portions")]
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
    audit_angles = Counter(round(180.0 - c["TurnDegrees"], 6) for c in census["Corners"] if c["Kind"] == "CORNER")
    gates.append(
        gate(
            "A1-vertex-census",
            not unassigned_vertices and manifest_feature_vertices == audit_feature_vertices and manifest_rounded == audit_rounded and manifest_angles == audit_angles and typed.get("Excluded", 0) == 0,
            {
                "AuditCornerEndpointJunction": audit_feature_vertices,
                "ManifestCornerEndpointJunction": manifest_feature_vertices,
                "AuditRoundedCorners": audit_rounded,
                "ManifestRoundedCorners": manifest_rounded,
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
    for kind, cls in EXCLUSION_CLASS_OF_AUDIT_KIND.items():
        audit_length = census["LengthByClass"].get(kind, 0.0)
        exclusion_detail[cls] = {"AuditLength": audit_length, "RecordedLength": recorded.get(cls, 0.0)}
        if audit_length > tol and recorded.get(cls, 0.0) < audit_length - 1.0e-6 * audit_length:
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


def run_audit(args):
    with open(args.config) as source:
        config = json.load(source)
    mesh = read_msh2(args.mesh)
    manifest = M.load_manifest(args.manifest)
    summary = M.summarize(manifest)
    radius = args.radius or summary["MatchingRadius"]
    if not radius:
        raise SystemExit("matching radius unknown: pass --radius")
    process_normal = None
    for entry in config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", []):
        if "EdgeFrameNormal" in entry:
            process_normal = entry["EdgeFrameNormal"]
    perimeter = P.extract_perimeter(mesh, config, process_normal=process_normal, corner_tolerance_degrees=args.corner_tolerance)
    targets = set(P.target_interfaces(config)) or {i for i, _ in P.interface_attributes(config)}
    census = perimeter_census(perimeter, radius, targets)
    clusters = corner_clusters(perimeter, radius)
    log = M.parse_palace_log(args.log) if args.log else None

    gates = []
    statistics = summary["Statistics"]
    identification = manifest.get("Identification")
    v2_summary = None
    if identification:
        gates, v2_summary = identification_gates(identification, perimeter, census, radius, targets)
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
        excluded_length = sum(census["LengthByClass"].get(k, 0.0) for k in ("NONPLANAR", "CROSS_LAYER", "NONMANIFOLD"))
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
        corner_angles = Counter(round(180.0 - c["TurnDegrees"], 6) for c in census["Corners"] if c["Kind"] == "CORNER")
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
        excluded_length = sum(census["LengthByClass"].get(k, 0.0) for k in ("NONPLANAR", "CROSS_LAYER", "NONMANIFOLD"))
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
