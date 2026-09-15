#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Validation helpers for frozen coupon mesh semantics.

The contract is deliberately independent of numeric label conventions.  Numeric
attributes are data; roles, materials, adjacency, corners, and protected
supports are all frozen before a mesh audit is run.
"""
import json
import math
from pathlib import Path


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
        if (not isinstance(item, dict) or not isinstance(item.get("Attribute"), int) or
                item["Attribute"] <= 0 or not isinstance(item.get("Role"), str) or
                not item["Role"].strip() or not isinstance(adjacent, list) or
                not adjacent or len(set(adjacent)) != len(adjacent) or
                not set(adjacent) <= material_set):
            raise ValueError("Each boundary label needs a role and nonempty valid adjacency")
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
    return data


def load_semantic_contract(path):
    return validate_semantic_contract(json.loads(Path(path).read_text()))


def volume_attributes(contract):
    return {item["Attribute"] for item in contract["VolumeMaterials"]}


def boundary_attributes(contract):
    return {item["Attribute"] for item in contract["BoundaryLabels"]}


def boundary_adjacency(contract):
    return {item["Attribute"]: set(item["AdjacentMaterials"])
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
            {"Attribute": 1, "Material": "vacuum"},
            {"Attribute": 2, "Material": "substrate"},
        ],
        "BoundaryLabels": [
            {"Attribute": 1, "Role": "outer", "AdjacentMaterials": [1],
             "Protected": True},
            {"Attribute": 3100, "Role": "matching", "AdjacentMaterials": [1, 2],
             "Protected": True},
            {"Attribute": 5001, "Role": "metal-substrate", "AdjacentMaterials": [1],
             "Protected": True},
            {"Attribute": 6001, "Role": "metal-air", "AdjacentMaterials": [2],
             "Protected": True},
        ],
        "SemanticCorners": [[0.0, 0.0, 0.0]],
        "ProtectedSupports": ["outer", "matching", "metal-substrate", "metal-air"],
        "MetricSurfaceRoles": ["metal-substrate", "metal-air"],
        "CutSurfaceRoles": ["outer"],
        "UnmatchedPolicy": "Error",
    })
