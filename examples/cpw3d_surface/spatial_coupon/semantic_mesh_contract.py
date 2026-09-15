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


def derive_feature_topology(signature_path, boundary_path, semantic_corners, tolerance=1e-9):
    with Path(signature_path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"Slot", "Conductor", "Px", "Py", "Pz", "Tx", "Ty", "Tz", "S0", "S1"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Signature lacks finite oriented segment columns")
    grouped = {}
    endpoints = []
    for row in rows:
        key, tangent, offset, lo, hi, first, last = _canonical_segment(row)
        grouped.setdefault(key, []).append((lo, hi, tangent, offset))
        endpoints.extend((first, last))
    features = 0
    subdivisions = []
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
    for point in endpoints:
        if len(corner_array) and np.linalg.norm(corner_array - point, axis=1).min() <= tolerance:
            continue
        if any(np.linalg.norm(np.asarray(other) - point) <= tolerance for other in subdivisions):
            continue
        if not any(np.linalg.norm(np.asarray(other) - point) <= tolerance for other in cuts):
            cuts.append(point.tolist())
    with Path(boundary_path).open(newline="") as stream:
        boundary = list(csv.DictReader(stream))
    if not boundary or "Class" not in boundary[0]:
        raise ValueError("Boundary contract lacks Physical/Continuation classes")
    classes = [row["Class"] for row in boundary]
    if any(value not in ("Physical", "Continuation") for value in classes):
        raise ValueError("Boundary contract has an unknown vertex class")
    canonical = lambda points: sorted([[float(value) for value in point] for point in points])
    return {"PhysicalFeatureCount": features,
            "CADSubdivisionCount": len(rows) - features,
            "BoundaryPhysicalVertexCount": classes.count("Physical"),
            "BoundaryContinuationVertexCount": classes.count("Continuation"),
            "CADSubdivisionEndpoints": canonical(subdivisions),
            "CutEndpoints": canonical(cuts)}


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
