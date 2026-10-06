#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation helpers for frozen coupon mesh semantics.

The contract is deliberately independent of numeric label conventions.  Numeric
attributes are data; roles, materials, adjacency, corners, and protected
supports are all frozen before a mesh audit is run.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np


REQUIRED_ROLES = ("Signature", "Boundary", "Mask", "Process", "SemanticContract")


def _unique_nonempty(values, name):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} must be a nonempty list")
    if len(set(values)) != len(values):
        raise ValueError(f"{name} must contain unique values")
    return values


def validate_semantic_contract(data):
    if not isinstance(data, dict) or data.get("Version") != 1:
        raise ValueError("Unsupported semantic mesh contract")
    volumes = data.get("VolumeMaterials")
    if not isinstance(volumes, list) or not volumes:
        raise ValueError("VolumeMaterials must be nonempty")
    volume_attributes = []
    for item in volumes:
        if (not isinstance(item, dict) or not isinstance(item.get("Attribute"), int) or
                item["Attribute"] <= 0 or not isinstance(item.get("Material"), str) or
                not item["Material"].strip()):
            raise ValueError("Each volume material needs a positive attribute and a name")
        volume_attributes.append(item["Attribute"])
    _unique_nonempty(volume_attributes, "volume attributes")

    boundaries = data.get("BoundaryLabels")
    if not isinstance(boundaries, list) or not boundaries:
        raise ValueError("BoundaryLabels must be nonempty")
    boundary_attributes = []
    protected_roles = set()
    material_set = set(volume_attributes)
    for item in boundaries:
        adjacent = item.get("AdjacentMaterials") if isinstance(item, dict) else None
        allowed = item.get("AdjacentMaterialSets", [adjacent]) if isinstance(item, dict) else None
        if (not isinstance(item, dict) or not isinstance(item.get("Attribute"), int) or
                item["Attribute"] <= 0 or not isinstance(item.get("Role"), str) or
                not item["Role"].strip() or not isinstance(adjacent, list) or
                not adjacent or len(set(adjacent)) != len(adjacent) or
                not set(adjacent) <= material_set or not isinstance(allowed, list) or
                not allowed or any(not isinstance(values, list) or not values or
                    not set(values) <= material_set for values in allowed) or
                set().union(*(set(values) for values in allowed)) != set(adjacent)):
            raise ValueError("Each boundary label needs a role and valid adjacency sets")
        boundary_attributes.append(item["Attribute"])
        if item.get("Protected") is True:
            protected_roles.add(item["Role"])
    _unique_nonempty(boundary_attributes, "boundary attributes")

    corners = data.get("SemanticCorners")
    if (not isinstance(corners, list) or not corners or
            any(not isinstance(point, list) or len(point) != 3 or
                any(not isinstance(value, (int, float)) or not math.isfinite(value)
                    for value in point)
                for point in corners)):
        raise ValueError("SemanticCorners must contain at least one 3D point")
    supports = _unique_nonempty(data.get("ProtectedSupports"), "ProtectedSupports")
    if not all(isinstance(value, str) and value.strip() for value in supports):
        raise ValueError("ProtectedSupports must contain names")
    if not protected_roles or not protected_roles <= set(supports):
        raise ValueError("Protected boundary roles must be named protected supports")
    if data.get("UnmatchedPolicy") != "Error":
        raise ValueError("Semantic contract must use UnmatchedPolicy = Error")
    metric_roles = _unique_nonempty(data.get("MetricSurfaceRoles"), "MetricSurfaceRoles")
    cut_roles = _unique_nonempty(data.get("CutSurfaceRoles"), "CutSurfaceRoles")
    boundary_roles = {item["Role"] for item in boundaries}
    if not set(metric_roles) <= boundary_roles or not set(cut_roles) <= boundary_roles:
        raise ValueError("Metric and cut surface roles must name boundary roles")
    topology = data.get("FeatureTopology")
    if (not isinstance(topology, dict) or
            not isinstance(topology.get("PhysicalFeatureCount"), int) or
            topology["PhysicalFeatureCount"] <= 0 or
            not isinstance(topology.get("CADSubdivisionCount"), int) or
            topology["CADSubdivisionCount"] < 0 or
            not isinstance(topology.get("BoundaryPhysicalVertexCount"), int) or
            topology["BoundaryPhysicalVertexCount"] < 0 or
            not isinstance(topology.get("BoundaryContinuationVertexCount"), int) or
            topology["BoundaryContinuationVertexCount"] < 0):
        raise ValueError("FeatureTopology must contain nonnegative explicit topology counts")
    for name in ("CADSubdivisionEndpoints", "CutEndpoints"):
        points = topology.get(name)
        if (not isinstance(points, list) or
                any(not isinstance(point, list) or len(point) != 3 or
                    any(not isinstance(value, (int, float)) or not math.isfinite(value)
                        for value in point) for point in points)):
            raise ValueError(f"FeatureTopology {name} must contain finite 3D points")
    return data


def load_semantic_contract(path):
    return validate_semantic_contract(json.loads(Path(path).read_text()))


def _canonical_segment(row):
    point = np.array([float(row[name]) for name in ("Px", "Py", "Pz")])
    tangent = np.array([float(row[name]) for name in ("Tx", "Ty", "Tz")])
    norm = np.linalg.norm(tangent)
    if not np.isfinite(norm) or norm <= 0:
        raise ValueError("Signature has an invalid tangent")
    tangent /= norm
    first = point + float(row["S0"]) * tangent
    last = point + float(row["S1"]) * tangent
    pivot = int(np.argmax(np.abs(tangent)))
    if tangent[pivot] < 0:
        tangent *= -1
        first, last = last, first
    offset = point - np.dot(point, tangent) * tangent
    lo, hi = sorted((float(np.dot(first, tangent)), float(np.dot(last, tangent))))
    key = (int(row["Slot"]), int(row["Conductor"]),
           *np.round(tangent, 10), *np.round(offset, 10))
    return key, tangent, offset, lo, hi, first, last


def box_face_cut_end(previous_point, previous_class, point, point_class, following_point):
    """Whether a plan-view boundary vertex is a BOX-FACE CUT END rather than a semantic
    corner (block (b) design A2 / A6, supervisor decision 320).

    A vertex carries the class of its OUTGOING side, so a Physical vertex whose incoming
    side is a Continuation (box) side has a single metal side there.  Such a vertex is a
    cut end - no corner ball, the tube ends on the face - when its metal side is not
    exactly perpendicular to the face: the exact-arithmetic test on the quantised
    canonical coordinates is that the side's face-parallel coordinate differs between its
    two ends (theta > 0).  An exactly perpendicular side (theta == 0: every rectilinear
    coupon) keeps the legacy convention bitwise and the vertex stays a semantic corner.
    Two Physical sides meeting at a vertex are a corner at any angle."""
    if point_class != "Physical" or previous_class != "Continuation":
        return False
    constant = [c for c in range(2) if previous_point[c] == point[c]]
    if len(constant) != 1:
        raise ValueError("a Continuation side must run along one box face "
                         f"({previous_point} -> {point})")
    face_parallel = 1 - constant[0]
    return point[face_parallel] != following_point[face_parallel]


BOX_FACE_CUT_END_RULE = ("supervisor decision 320: a Physical vertex whose incoming side is a Continuation "
                         "(box) side has a single metal side there; it is a box-face cut end, not a "
                         "semantic corner, unless that side is exactly perpendicular to the face (its "
                         "face-parallel coordinate equal at both ends, exact arithmetic), in which case "
                         "the legacy rectilinear convention keeps it a corner bitwise")


def invariant_corner(previous_point, point, following_point):
    """Whether a semantic corner is INVARIANT (mesher design round 2 F5-A, supervisor
    decisions 351 / 358 / 363): the two plan-view boundary sides meeting at it are not
    exactly perpendicular.  Exact arithmetic on the quantised canonical coordinates: an
    axis-aligned side has an exactly zero component, so the dot product of the two side
    vectors (pointing away from the corner) is exactly 0.0 at every rectilinear corner and at
    the theta-0 box vertex of decision 320 (its metal side perpendicular to the box side);
    such a LEGACY corner keeps the vertex-0 corner measure, MaximumCornerAspect and the 3.8
    target bitwise.  An invariant corner is optimized on and judged by kappa_reg (the
    condition number of the affine map from the regular tetrahedron, order-invariant),
    descended to the fixed goal 3.8 and judged against the manifest's CornerShapeGate =
    min(E_pop, 5.0); a BridgingSliver candidate (a corner-incident cell whose four vertices
    all lie on the kink's two sidewalls, at least one strictly on each, in a wedge of obtuse
    opening) still above the gate after the descent triggers the corner-local reconnection
    pass (supervisor decision 365: the measure is the verdict, the predicate the trigger);
    the mesher (mesh_spatial_coupon.semantic_corner_kinds) evaluates the same predicate and
    fails closed on a contract disagreeing with it."""
    a = (previous_point[0] - point[0], previous_point[1] - point[1])
    b = (following_point[0] - point[0], following_point[1] - point[1])
    return a[0] * b[0] + a[1] * b[1] != 0.0


INVARIANT_CORNER_RULE = ("mesher design round 2 F5-A (supervisor decisions 351 / 358 / 363 / 365): a semantic "
                         "corner whose two plan-view boundary sides have a non-zero dot product on their "
                         "quantised coordinates (exact arithmetic; every rectilinear corner and the theta-0 box "
                         "vertex read exactly 0 and stay LEGACY, bitwise) is INVARIANT: its corner-incident seed "
                         "cells are optimized on and judged by kappa_reg, the condition number of the affine "
                         "map from the regular tetrahedron (order-invariant), descended to the fixed goal 3.8 "
                         "and judged against the manifest CornerShapeGate = min(E_pop, 5.0) (E_pop the kappa_reg "
                         "envelope of the (F)-qualified 90-degree corners); a BridgingSliver candidate (a "
                         "corner-incident cell whose four vertices all lie on the kink's two sidewalls, at least "
                         "one strictly on each, in a wedge of obtuse opening) still above the gate after the "
                         "descent triggers the corner-local reconnection pass (supervisor decision 365: the "
                         "measure is the verdict, the predicate the trigger); the mesher evaluates the same "
                         "predicate and a disagreeing contract fails closed")


def boundary_arc_tags(rows):
    """Per plan-view boundary row (file order) the arc tag of its OUTGOING side ({ArcId, ArcCx,
    ArcCy, ArcR, ArcSign} or None) and the joint record of its vertex ((turn or None, smooth)
    or None), from the ARC columns of block (b) design A1 (4) (generate_spatial_response.
    ARC_BOUNDARY_COLUMNS); every row None on a boundary without the columns (a legacy coupon)."""
    if not rows or "ArcId" not in rows[0]:
        return [None] * len(rows), [None] * len(rows)
    arcs, joints = [], []
    for row in rows:
        if row.get("ArcId") in (None, ""):
            arcs.append(None)
        else:
            arcs.append({"ArcId": int(row["ArcId"]), "ArcCx": float(row["ArcCx"]), "ArcCy": float(row["ArcCy"]),
                         "ArcR": float(row["ArcR"]), "ArcSign": int(row["ArcSign"])})
        if row.get("JointSmooth") in (None, ""):
            joints.append(None)
        else:
            turn = row.get("JointTurn")
            joints.append((None if turn in (None, "") else float(turn), int(row["JointSmooth"]) == 1))
    return arcs, joints


ARC_VERTEX_MATCH_TOLERANCE = 1.0e-6
ARC_VERTEX_RULE = ("block (b) design 1.2 (1) / A3 (1) (decision 303): a plan-view vertex strictly inside one arc "
                   "(both adjacent sides chords of the same ArcId: ArcInterior) is never a semantic corner; a vertex "
                   "ending an arc (ArcJoint) is a semantic corner only when its signature turn exceeds "
                   "JUNCTION_TANGENT_ANGLE = 1e-4 rad (JointSmooth 0); a smooth joint shares one tube section")


def arc_vertex_class(arcs, index):
    """'Interior' (both adjacent sides of one ArcId), 'Joint' (one adjacent arc side, or two of
    different arcs) or None (no adjacent arc side) for vertex `index` of a loop whose outgoing
    side tags are `arcs` (loop order)."""
    outgoing = arcs[index]
    incoming = arcs[index - 1]
    if outgoing is None and incoming is None:
        return None
    if outgoing is not None and incoming is not None and outgoing["ArcId"] == incoming["ArcId"]:
        return "Interior"
    return "Joint"


def boundary_semantic_corners(rows):
    """Semantic corners, box-face cut ends and invariant corners of plan-view boundary rows
    (Loop / Vertex / Class / X / Y / Plane, file order): the Physical vertices that are not
    box-face cut ends (box_face_cut_end), each at its Plane height, in file order, the
    excluded cut ends likewise, and the corners among them whose sides are not exactly
    perpendicular (invariant_corner).  Arc vertices (the ARC columns, design A1 (4)): an
    ArcInterior vertex and a smooth ArcJoint are never corners (ARC_VERTEX_RULE); a kinked
    ArcJoint is a corner like a Physical vertex (boundary_arc_vertices lists the excluded arc
    vertices)."""
    corners, cut_ends, invariant, _ = _boundary_vertex_classes(rows)
    return corners, cut_ends, invariant


def boundary_arc_vertices(rows):
    """The arc vertices boundary_semantic_corners excludes: {"Interior": [...], "SmoothJoints":
    [...]} at their Plane heights (file order); both empty on a boundary without arcs."""
    return _boundary_vertex_classes(rows)[3]


def _boundary_vertex_classes(rows):
    loops = {}
    for row, arc, joint in zip(rows, *boundary_arc_tags(rows)):
        loops.setdefault(row["Loop"], []).append((row, arc, joint))
    corners = []
    cut_ends = []
    invariant = []
    arc_vertices = {"Interior": [], "SmoothJoints": []}
    for loop in loops.values():
        points = [(float(row["X"]), float(row["Y"])) for row, _, _ in loop]
        classes = [row["Class"] for row, _, _ in loop]
        arcs = [arc for _, arc, _ in loop]
        n = len(loop)
        for index, (row, _, joint) in enumerate(loop):
            point = [points[index][0], points[index][1], float(row["Plane"])]
            vertex_class = arc_vertex_class(arcs, index)
            if vertex_class == "Interior":
                arc_vertices["Interior"].append(point)
                continue
            if vertex_class == "Joint" and joint is not None and joint[1]:
                arc_vertices["SmoothJoints"].append(point)
                continue
            if classes[index] != "Physical":
                continue
            if box_face_cut_end(points[index - 1], classes[index - 1], points[index],
                                classes[index], points[(index + 1) % n]):
                cut_ends.append(point)
            else:
                corners.append(point)
                if invariant_corner(points[index - 1], points[index], points[(index + 1) % n]):
                    invariant.append(point)
    return corners, cut_ends, invariant, arc_vertices


def invariant_corners(contract):
    """The contract's recorded invariant corners (Derivation.InvariantCorners.Points; an
    empty list when the record is absent: every corner legacy)."""
    derivation = contract.get("Derivation") if isinstance(contract, dict) else None
    record = derivation.get("InvariantCorners") if isinstance(derivation, dict) else None
    if record is None:
        return []
    points = record.get("Points") if isinstance(record, dict) else None
    if (not isinstance(points, list) or
            any(not isinstance(point, list) or len(point) != 3 or
                any(not isinstance(value, (int, float)) or not math.isfinite(value) for value in point)
                for point in points)):
        raise ValueError("Derivation.InvariantCorners.Points must contain finite 3D points")
    corners = contract.get("SemanticCorners", [])
    if any(all(np.linalg.norm(np.asarray(point) - np.asarray(corner)) > 1e-9 for corner in corners)
           for point in points):
        raise ValueError("Derivation.InvariantCorners.Points must be semantic corners of the contract")
    return points


def derive_feature_topology(signature_path, boundary_path, semantic_corners, tolerance=1e-9):
    with Path(signature_path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"Slot", "Conductor", "Px", "Py", "Pz", "Tx", "Ty", "Tz", "S0", "S1"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Signature lacks finite oriented segment columns")
    with Path(boundary_path).open(newline="") as stream:
        boundary = list(csv.DictReader(stream))
    if not boundary or "Class" not in boundary[0]:
        raise ValueError("Boundary contract lacks Physical/Continuation classes")
    classes = [row["Class"] for row in boundary]
    if any(value not in ("Physical", "Continuation") for value in classes):
        raise ValueError("Boundary contract has an unknown vertex class")
    # Arc chords (block (b) design 1.2 (1) / (5): "FeatureTopology counts arcs"): a signature
    # row whose ends are the two ends of a tagged chord of the plan-view boundary belongs to
    # that arc; one arc is ONE physical feature, its further chord rows are CAD subdivisions
    # of it, and the chord vertices strictly inside the arc are ArcInteriorEndpoints (the
    # tube is smooth there: neither a cut nor a subdivision end).  A smooth arc joint is a
    # CAD-subdivision end (one continuous metal edge, two CAD pieces); a kinked joint is a
    # semantic corner (semantic_corners) and a box-face arc end a cut end, as for lines.
    arc_tags, _ = boundary_arc_tags(boundary)
    chord_arc = {}
    if any(tag is not None for tag in arc_tags):
        loops = {}
        for row, tag in zip(boundary, arc_tags):
            loops.setdefault(row["Loop"], []).append((row, tag))
        for loop in loops.values():
            points = [np.asarray([float(row["X"]), float(row["Y"]), float(row["Plane"])]) for row, _ in loop]
            for index, (_, tag) in enumerate(loop):
                if tag is not None:
                    chord_arc[len(chord_arc)] = (tag["ArcId"], points[index], points[(index + 1) % len(loop)])
    arc_vertices = boundary_arc_vertices(boundary)
    interior = [np.asarray(point, dtype=float) for point in arc_vertices["Interior"]]
    smooth = [np.asarray(point, dtype=float) for point in arc_vertices["SmoothJoints"]]
    # The chord vertices are computed points quantised to the boundary's 1e-9 R grid (the
    # straight ends lie on the 1e-6 R signature grid and reproduce exactly): an arc vertex is
    # matched within ARC_VERTEX_MATCH_TOLERANCE (a nanometre; no two plan-view feature points
    # lie that close).
    arc_tolerance = max(tolerance, ARC_VERTEX_MATCH_TOLERANCE)

    def chord_of(first, last):
        for arc_id, a, b in chord_arc.values():
            if ((np.linalg.norm(a - first) <= arc_tolerance and np.linalg.norm(b - last) <= arc_tolerance) or
                    (np.linalg.norm(a - last) <= arc_tolerance and np.linalg.norm(b - first) <= arc_tolerance)):
                return arc_id
        return None

    grouped = {}
    endpoints = []
    arc_features = set()
    for row in rows:
        key, tangent, offset, lo, hi, first, last = _canonical_segment(row)
        arc_id = chord_of(first, last) if chord_arc else None
        if arc_id is not None:
            arc_features.add(arc_id)
        else:
            grouped.setdefault(key, []).append((lo, hi, tangent, offset))
        endpoints.extend((first, last))
    features = len(arc_features)
    subdivisions = [point.tolist() for point in smooth]
    for segments in grouped.values():
        segments.sort(key=lambda item: (item[0], item[1]))
        components = []
        for lo, hi, tangent, offset in segments:
            if components and lo <= components[-1][1] + tolerance:
                if abs(lo - components[-1][1]) <= tolerance:
                    subdivisions.append((offset + lo * tangent).tolist())
                components[-1] = (components[-1][0], max(components[-1][1], hi))
            else:
                components.append((lo, hi))
        features += len(components)
    corner_array = np.asarray(semantic_corners, dtype=float).reshape(-1, 3)
    cuts = []
    arc_interior = []
    for point in endpoints:
        if len(corner_array) and np.linalg.norm(corner_array - point, axis=1).min() <= tolerance:
            continue
        if any(np.linalg.norm(np.asarray(other) - point) <= tolerance for other in subdivisions):
            continue
        if any(np.linalg.norm(other - point) <= arc_tolerance for other in smooth):
            continue
        if any(np.linalg.norm(other - point) <= arc_tolerance for other in interior):
            if not any(np.linalg.norm(np.asarray(other) - point) <= tolerance for other in arc_interior):
                arc_interior.append(point.tolist())
            continue
        if not any(np.linalg.norm(np.asarray(other) - point) <= tolerance for other in cuts):
            cuts.append(point.tolist())
    canonical = lambda points: sorted([[float(value) for value in point] for point in points])
    topology = {"PhysicalFeatureCount": features,
                "CADSubdivisionCount": len(rows) - features,
                "BoundaryPhysicalVertexCount": classes.count("Physical"),
                "BoundaryContinuationVertexCount": classes.count("Continuation"),
                "CADSubdivisionEndpoints": canonical(subdivisions),
                "CutEndpoints": canonical(cuts)}
    if arc_features:
        # Recorded only where arcs exist, so every straight contract is unchanged.
        topology["ArcFeatureCount"] = len(arc_features)
        topology["ArcInteriorEndpoints"] = canonical(arc_interior)
        topology["ArcVertexRule"] = ARC_VERTEX_RULE
    return topology


def validate_feature_topology(contract, signature_path, boundary_path):
    actual = derive_feature_topology(signature_path, boundary_path,
                                     contract["SemanticCorners"])
    if contract["FeatureTopology"] != actual:
        raise ValueError("FeatureTopology differs from finite segments and boundary classes")
    return actual


def volume_attributes(contract):
    return {item["Attribute"] for item in contract["VolumeMaterials"]}


def boundary_attributes(contract):
    return {item["Attribute"] for item in contract["BoundaryLabels"]}


def boundary_adjacency(contract):
    return {item["Attribute"]: [set(values) for values in
            item.get("AdjacentMaterialSets", [item["AdjacentMaterials"]])]
            for item in contract["BoundaryLabels"]}


def metric_surface_attributes(contract):
    roles = set(contract["MetricSurfaceRoles"])
    return {item["Attribute"] for item in contract["BoundaryLabels"]
            if item["Role"] in roles}


def cut_surface_attributes(contract):
    roles = set(contract["CutSurfaceRoles"])
    return {item["Attribute"] for item in contract["BoundaryLabels"]
            if item["Role"] in roles}


def material_interface_attributes(contract):
    """Boundary labels whose triangles separate two volume materials.

    These are the dielectric interfaces (the etched trench floor and walls and
    the un-etched substrate-vacuum plane): an adjacency set with two materials
    means both materials meet across the same triangle.  Cut-surface roles are
    excluded even when their adjacency lists both materials, because each cut
    triangle touches one material (`AdjacentMaterialSets` [[1], [2]]); conductor
    labels touch one material and are excluded by the same rule.
    """
    cut = cut_surface_attributes(contract)
    return {attribute for attribute, sets in boundary_adjacency(contract).items()
            if attribute not in cut and any(len(values) >= 2 for values in sets)}


def simple_sharp_contract():
    """Explicit compatibility contract for the historical one-slot scout."""
    return validate_semantic_contract({
        "Version": 1,
        "VolumeMaterials": [
            {"Attribute": 1, "Material": "substrate"},
            {"Attribute": 2, "Material": "vacuum"},
        ],
        "BoundaryLabels": [
            {"Attribute": 1, "Role": "matching-surface", "AdjacentMaterials": [1, 2],
             "AdjacentMaterialSets": [[1], [2]], "Protected": True},
            {"Attribute": 3000, "Role": "substrate-vacuum", "AdjacentMaterials": [1, 2],
             "Protected": True},
            {"Attribute": 3100, "Role": "matching", "AdjacentMaterials": [1, 2],
             "Protected": True},
            {"Attribute": 5001, "Role": "metal-substrate", "AdjacentMaterials": [1],
             "Protected": True},
            {"Attribute": 6001, "Role": "metal-air", "AdjacentMaterials": [2],
             "Protected": True},
        ],
        "SemanticCorners": [[0.0, 0.0, 0.0]],
        "ProtectedSupports": ["matching-surface", "substrate-vacuum", "matching",
                              "metal-substrate", "metal-air"],
        "MetricSurfaceRoles": ["metal-substrate", "metal-air"],
        "CutSurfaceRoles": ["matching-surface"],
        "UnmatchedPolicy": "Error",
        "FeatureTopology": {
            "PhysicalFeatureCount": 1,
            "CADSubdivisionCount": 0,
            "BoundaryPhysicalVertexCount": 1,
            "BoundaryContinuationVertexCount": 1,
            "CADSubdivisionEndpoints": [],
            "CutEndpoints": [],
        },
    })
