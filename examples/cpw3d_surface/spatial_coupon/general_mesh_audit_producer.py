#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Produce mesh-bound records for the geometry-independence evidence assembler."""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import sys
import tomllib

import meshio
import numpy as np

from audit_edge_metric_mesh import analyze, blocks, directional_widths
from general_mesh_manifest import canonical_sha256, sha256
from mesh_array_io import read_mesh
from mesh_stage_contract import (CANONICAL_STAGE_ORDER, PLACEMENT_STAGE_ORDER, STAGE_ORDER,
                                 sha256 as stage_sha256, validate_stage_dag)
from semantic_mesh_contract import load_semantic_contract


KINDS = ("bounded-run", "mesh-topology-quality", "mesh-complexity",
         "mesh-invariants", "variant-transform")


def _producer():
    path = Path(__file__).resolve()
    return {"Name": path.name, "SHA256": sha256(path)}


def _base(kind, case, variant, mesh, input_hashes, transform, command):
    return {"Version": 1, "Kind": kind, "CaseId": case, "Variant": variant,
            "MeshSHA256": sha256(mesh), "InputSHA256": input_hashes,
            "Transform": transform, "TransformSHA256": canonical_sha256(transform),
            "Command": list(command),
            "Environment": {name: os.environ.get(name, "") for name in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "JULIA_NUM_THREADS")},
            "Producer": _producer()}


def _global_diagonal_bands(mesh, physical_segments, normal_size):
    """Find long, narrow short-edge bands on one planar labeled support.

    A connected set spanning several orthogonal supports is not a geometric
    band: its bounding-box diagonal has no source meaning.  Components are
    therefore formed independently on the exact labeled planes.  The short
    scale is also bounded by twice the audited transverse target, so ordinary
    coarse-surface triangulation cannot masquerade as propagated trace sizing.
    """
    triangles, labels = blocks(mesh, "triangle")
    xyz = mesh.points[triangles]
    normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    pivot = np.argmax(abs(normals), axis=1)
    normals *= np.sign(normals[np.arange(len(normals)), pivot])[:, None]
    planes = np.round(np.column_stack(
        (normals, np.einsum("ij,ij->i", normals, xyz[:, 0]))), 8)
    planes[planes == 0] = 0.0
    _, patch = np.unique(np.column_stack((labels, planes)), axis=0, return_inverse=True)
    owners_by_edge, all_lengths = {}, []
    for owner, triangle in enumerate(triangles):
        for a, b in ((triangle[0], triangle[1]), (triangle[1], triangle[2]),
                     (triangle[2], triangle[0])):
            edge = tuple(sorted((int(a), int(b))))
            owners_by_edge.setdefault(edge, []).append(owner)
            all_lengths.append(np.linalg.norm(mesh.points[edge[0]] - mesh.points[edge[1]]))
    median = float(np.median(all_lengths))
    threshold = min(.6 * median, 2.0 * float(normal_size))
    short_by_patch = {}
    threshold_tolerance = 64.0 * np.finfo(float).eps * max(median, threshold, 1.0)
    for edge, owners in owners_by_edge.items():
        length = np.linalg.norm(mesh.points[edge[0]] - mesh.points[edge[1]])
        if (len(owners) == 2 and patch[owners[0]] == patch[owners[1]] and
                length <= threshold + threshold_tolerance):
            short_by_patch.setdefault(int(patch[owners[0]]), []).append(edge)
    segments = np.asarray(physical_segments, dtype=float).reshape(-1, 2, 3)
    bands, components = 0, []
    for patch_id, short in short_by_patch.items():
        adjacency = {}
        for first, last in short:
            adjacency.setdefault(first, set()).add(last)
            adjacency.setdefault(last, set()).add(first)
        patch_vertices = np.unique(triangles[patch == patch_id])
        surface_diameter = max(float(np.linalg.norm(
            np.ptp(mesh.points[patch_vertices], axis=0))), 1e-300)
        visited = set()
        for start in adjacency:
            if start in visited:
                continue
            stack = [start]; component = set()
            while stack:
                vertex = stack.pop()
                if vertex in component:
                    continue
                component.add(vertex); visited.add(vertex)
                stack.extend(adjacency.get(vertex, ()))
            points = mesh.points[list(component)]
            centered = points - points.mean(axis=0)
            _, singular, axes = np.linalg.svd(centered, full_matrices=False)
            projection = centered @ axes[0]
            first, last = int(np.argmin(projection)), int(np.argmax(projection))
            span = float(projection[last] - projection[first])
            if span <= .5 * surface_diameter:
                continue
            width = float(singular[1] / math.sqrt(len(points))) if len(singular) > 1 else 0.0
            line_like = width <= 2.0 * threshold
            direction = axes[0]
            alignment = 0.0
            for segment in segments:
                tangent = segment[1] - segment[0]
                tangent /= np.linalg.norm(tangent)
                alignment = max(alignment, float(abs(np.dot(direction, tangent))))
            aligned = alignment > 1 - 1e-6
            bands += int(line_like and not aligned)
            owner = int(np.flatnonzero(patch == patch_id)[0])
            components.append({
                "Attribute": int(labels[owner]), "Plane": planes[owner].tolist(),
                "Endpoints": [points[first].tolist(), points[last].tolist()],
                "Span": span, "RMSWidth": width, "PhysicalSegmentAlignment": alignment,
                "LineLike": line_like, "AlignedWithPhysicalSegment": aligned,
                "Vertices": len(component)})
    return {"GlobalDiagonalBands": bands,
            "ShortInternalEdges": sum(map(len, short_by_patch.values())),
            "MedianSurfaceEdgeLength": median, "ShortEdgeThreshold": threshold,
            "SizeFieldEvidence": {"NormalSize": float(normal_size),
                                  "MaximumPropagatedShortEdge": 2.0 * float(normal_size)},
            "LongShortEdgeComponents": components}


def _tetra_quality(mesh):
    tetrahedra, _ = blocks(mesh, "tetra")
    xyz = mesh.points[tetrahedra]
    jacobian = np.stack((xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0],
                         xyz[:, 3] - xyz[:, 0]), axis=2)
    determinant = np.linalg.det(jacobian)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    condition = singular[:, 0] / singular[:, -1]
    lengths = np.linalg.norm(jacobian, axis=1)
    scaled = np.abs(determinant) / np.prod(lengths, axis=1)
    if (not len(scaled) or not np.all(np.isfinite(condition)) or
            np.any(singular[:, -1] <= 0)):
        raise ValueError("Invalid tetrahedral Jacobian audit")
    quantiles = (0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0)
    return {"Samples": len(tetrahedra),
            "PositiveOrientation": bool(np.all(determinant > 0.0)),
            "MinimumScaledJacobian": float(scaled.min()),
            "MaximumJacobianCondition": float(condition.max()),
            "ScaledJacobianQuantiles": np.quantile(scaled, quantiles).tolist(),
            "JacobianConditionQuantiles": np.quantile(condition, quantiles).tolist()}


def _point_aspects(mesh, points):
    tetrahedra, _ = blocks(mesh, "tetra")
    xyz = mesh.points[tetrahedra]
    centers = xyz.mean(axis=1)
    jacobian = np.stack((xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0],
                         xyz[:, 3] - xyz[:, 0]), axis=2)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    aspects = singular[:, 0] / singular[:, -1]
    result = []
    for point in np.asarray(points, dtype=float).reshape(-1, 3):
        incident = np.flatnonzero(np.any(np.linalg.norm(xyz - point, axis=2) <= 1e-10, axis=1))
        selected = incident if len(incident) else [int(np.argmin(np.linalg.norm(centers-point, axis=1)))]
        result.append({"Point": point.tolist(), "MaximumAspect": float(aspects[selected].max()),
                       "Cells": int(len(selected))})
    return result


def _point_to_segments(points, segments):
    if not len(points) or not len(segments):
        return math.inf
    result = 0.0
    first, vector = segments[:, 0], segments[:, 1] - segments[:, 0]
    denominator = np.einsum("ij,ij->i", vector, vector)
    if np.any(denominator <= 0):
        raise ValueError("Degenerate protected-footprint boundary segment")
    for start in range(0, len(points), 512):
        delta = points[start:start + 512, None] - first[None]
        parameter = np.clip(np.einsum("ijk,jk->ij", delta, vector) / denominator, 0, 1)
        distance = np.linalg.norm(delta - parameter[..., None] * vector, axis=2)
        result = max(result, float(np.min(distance, axis=1).max()))
    return result


def _patch_edges(triangles):
    raw_edges = np.sort(
        triangles[:, ((0, 1), (1, 2), (2, 0))].reshape(-1, 2), axis=1)
    edges, count = np.unique(raw_edges, axis=0, return_counts=True)
    if np.any(count > 2):
        raise ValueError("Nonmanifold protected planar patch")
    return edges[count == 1]


def _split_pinch_vertices(points, triangles):
    """Give every wedge of a pinched boundary vertex its own vertex copy.

    Two sub-regions of one label may touch at a single vertex (a device footprint
    touching the coupon side, or a per-triangle slot labeling artifact).  Such a
    pinch vertex is a legitimate planar topology whose boundary has even degree
    2k: the incident triangles form k edge-connected fans, each with exactly two
    boundary edges.  Splitting the vertex per fan lets the boundary be traced as
    closed loops without changing the segment geometry.  Odd-degree (open)
    boundaries remain errors.  Returns the split points, triangles and the pinch
    coordinates.
    """
    boundary = _patch_edges(triangles)
    degree = {}
    for first, last in boundary:
        degree[int(first)] = degree.get(int(first), 0) + 1
        degree[int(last)] = degree.get(int(last), 0) + 1
    pinched = sorted(vertex for vertex, value in degree.items() if value > 2)
    if not pinched:
        return points, triangles, []
    triangles = triangles.copy()
    points = list(points)
    boundary_set = {tuple(map(int, edge)) for edge in boundary}
    coordinates = []
    for vertex in pinched:
        incident = np.flatnonzero(np.any(triangles == vertex, axis=1))
        # Fans: triangles incident to the vertex, connected through the interior
        # edges that contain the vertex.
        remaining = set(map(int, incident))
        fans = []
        while remaining:
            seed = remaining.pop(); fan = {seed}; stack = [seed]
            while stack:
                owner = stack.pop()
                others = [int(v) for v in triangles[owner] if v != vertex]
                for other in others:
                    edge = tuple(sorted((vertex, other)))
                    if edge in boundary_set:
                        continue
                    for candidate in list(remaining):
                        if other in triangles[candidate]:
                            remaining.remove(candidate); fan.add(candidate)
                            stack.append(candidate)
            fans.append(sorted(fan))
        if len(fans) * 2 != degree[vertex]:
            raise ValueError("Protected planar patch has an open or branched boundary")
        coordinates.append([float(value) for value in points[vertex]])
        for fan in fans[1:]:
            copy = len(points); points.append(points[vertex])
            for owner in fan:
                triangles[owner][triangles[owner] == vertex] = copy
    return np.asarray(points), triangles, coordinates


def _normalized_footprint_boundary(xyz, pinches=None):
    """Return subdivision-independent boundary segments and component/hole topology.

    Pinch vertices are resolved per wedge (see _split_pinch_vertices); their
    coordinates are appended to `pinches` when a list is given and they stay
    segment endpoints, so the normalized segment set is unchanged."""
    points, inverse = np.unique(xyz.reshape(-1, 3), axis=0, return_inverse=True)
    triangles = inverse.reshape(-1, 3)
    points, triangles, pinch_coordinates = _split_pinch_vertices(points, triangles)
    if pinches is not None:
        pinches.extend(pinch_coordinates)
    pinch_vertices = {index for index in range(len(points))
                      if any(np.array_equal(points[index], coordinate)
                             for coordinate in pinch_coordinates)}
    boundary = _patch_edges(triangles)
    adjacency = {}
    for first, last in boundary:
        adjacency.setdefault(int(first), set()).add(int(last))
        adjacency.setdefault(int(last), set()).add(int(first))
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError("Protected planar patch has an open or branched boundary")

    # Count triangle components and associate each closed boundary loop with its
    # owning component.  A connected planar component has one exterior loop;
    # every additional loop is a hole.  This distinguishes disconnected disks
    # from an annulus even when their aggregate boundary geometry is similar.
    triangle_neighbors = [set() for _ in triangles]
    owners = {}
    for owner, triangle in enumerate(triangles):
        for first, last in ((triangle[0], triangle[1]), (triangle[1], triangle[2]),
                            (triangle[2], triangle[0])):
            owners.setdefault(tuple(sorted((int(first), int(last)))), []).append(owner)
    for edge_owners in owners.values():
        if len(edge_owners) == 2:
            first, last = edge_owners
            triangle_neighbors[first].add(last)
            triangle_neighbors[last].add(first)
    triangle_component = np.full(len(triangles), -1, dtype=int)
    component_count = 0
    for start in range(len(triangles)):
        if triangle_component[start] >= 0:
            continue
        stack = [start]
        while stack:
            owner = stack.pop()
            if triangle_component[owner] >= 0:
                continue
            triangle_component[owner] = component_count
            stack.extend(triangle_neighbors[owner])
        component_count += 1

    edge_owner = {edge: values[0] for edge, values in owners.items()
                  if len(values) == 1}
    loops_by_component = [0] * component_count
    unvisited = {tuple(map(int, edge)) for edge in boundary}
    while unvisited:
        first, last = min(unvisited)
        start, previous, current = first, first, last
        unvisited.remove((first, last))
        component = int(triangle_component[edge_owner[(first, last)]])
        while current != start:
            choices = adjacency[current] - {previous}
            if len(choices) != 1:
                raise ValueError("Protected planar boundary is not a closed cycle")
            following = next(iter(choices))
            edge = tuple(sorted((current, following)))
            if edge not in unvisited:
                raise ValueError("Protected planar boundary cycle is inconsistent")
            if int(triangle_component[edge_owner[edge]]) != component:
                raise ValueError("Protected boundary loop crosses surface components")
            unvisited.remove(edge)
            previous, current = current, following
        loops_by_component[component] += 1

    # Remove degree-two collinear subdivision vertices.  The resulting segment
    # set represents the PL boundary connectivity rather than only its vertices.
    tolerance = 64.0 * np.finfo(float).eps
    changed = True
    while changed:
        changed = False
        for vertex in sorted(adjacency):
            if (vertex not in adjacency or len(adjacency[vertex]) != 2 or
                    vertex in pinch_vertices):
                continue
            first, last = sorted(adjacency[vertex])
            left, right = points[first] - points[vertex], points[last] - points[vertex]
            scale = max(np.linalg.norm(left) * np.linalg.norm(right), 1e-300)
            if (np.linalg.norm(np.cross(left, right)) <= tolerance * scale and
                    np.dot(left, right) < 0):
                adjacency[first].remove(vertex); adjacency[first].add(last)
                adjacency[last].remove(vertex); adjacency[last].add(first)
                del adjacency[vertex]
                changed = True
                break
    normalized = []
    for first, neighbors in adjacency.items():
        for last in neighbors:
            if first < last:
                normalized.append([points[first], points[last]])
    topology = {"Components": component_count,
                "BoundaryLoops": int(sum(loops_by_component)),
                "Holes": int(sum(value - 1 for value in loops_by_component)),
                "LoopsPerComponent": sorted(map(int, loops_by_component))}
    return np.asarray(normalized, dtype=float).reshape(-1, 2, 3), topology


def _footprint_boundary_comparison(left_triangles, right_triangles, diagnostics=None):
    """Compare normalized PL segment sets and planar component/hole topology.

    `diagnostics`, when a dict, receives the pinch vertices of both patches."""
    left_pinches, right_pinches = [], []
    left, left_topology = _normalized_footprint_boundary(left_triangles, left_pinches)
    right, right_topology = _normalized_footprint_boundary(right_triangles, right_pinches)
    if diagnostics is not None:
        diagnostics["PinchVertices"] = {
            "Reference": len(left_pinches), "Candidate": len(right_pinches),
            "ReferenceCoordinates": left_pinches, "CandidateCoordinates": right_pinches}
    # Endpoints alone discard connectivity.  Including each normalized segment
    # midpoint detects same-vertex rewiring while remaining invariant to boundary
    # subdivision (which normalization removes).
    def witnesses(segments):
        return np.concatenate((segments[:, 0], segments[:, 1], segments.mean(axis=1)))
    distance = max(_point_to_segments(witnesses(left), right),
                   _point_to_segments(witnesses(right), left))
    return distance, left_topology, right_topology


def _footprint_boundary_distance(left_triangles, right_triangles):
    """Compatibility helper returning normalized segment-set distance."""
    return _footprint_boundary_comparison(left_triangles, right_triangles)[0]


def _protected_surface_report(reference, candidate, contract):
    before, _ = analyze(reference, contract, require_material_names=True)
    after, _ = analyze(candidate, contract, require_material_names=True)
    protected_attributes = {item["Attribute"]: item["Role"]
                            for item in contract["BoundaryLabels"]
                            if item.get("Protected") is True}

    def patches(mesh, report):
        triangles, labels = blocks(mesh, "triangle")
        result = {}
        for key, area in report["PlanarPatchAreas"].items():
            values = tuple(float(value) for value in key.split())
            attribute, plane = int(values[0]), values[1:]
            if attribute not in protected_attributes:
                continue
            xyz = mesh.points[triangles[labels == attribute]]
            cross = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
            normals = cross / np.linalg.norm(cross, axis=1)[:, None]
            pivot = np.argmax(abs(normals), axis=1)
            normals *= np.sign(normals[np.arange(len(normals)), pivot])[:, None]
            triangle_planes = np.round(np.column_stack(
                (normals, np.einsum("ij,ij->i", normals, xyz[:, 0]))), 8)
            triangle_planes[triangle_planes == 0] = 0.0
            selected = xyz[np.all(triangle_planes == np.asarray(plane), axis=1)]
            result[(protected_attributes[attribute], plane)] = (float(area), selected)
        return result

    left, right = patches(reference, before), patches(candidate, after)
    if left.keys() != right.keys():
        return {"Actual": sorted(protected_attributes.values()),
                "PlaneSupportsMatch": False, "MaximumRelativeMeasureError": math.inf,
                "MaximumSupportVertexDistance": math.inf, "PatchCount": len(right)}
    area_errors, distances, topology_matches, by_patch = [], [], [], {}
    for key in left:
        left_area, left_triangles = left[key]
        right_area, right_triangles = right[key]
        area_error = abs(right_area - left_area) / max(abs(left_area), 1e-300)
        diagnostics = {}
        distance, left_topology, right_topology = _footprint_boundary_comparison(
            left_triangles, right_triangles, diagnostics)
        topology_match = left_topology == right_topology
        area_errors.append(float(area_error)); distances.append(float(distance))
        topology_matches.append(topology_match)
        by_patch[f"{key[0]} {' '.join(map(str, key[1]))}"] = {
            "RelativeMeasureError": float(area_error),
            "SupportVertexDistance": float(distance),
            "TopologyMatches": topology_match,
            "ReferenceTopology": left_topology,
            "CandidateTopology": right_topology,
            # Diagnostic only: pinched boundary vertices resolved per wedge.
            "PinchVertices": diagnostics["PinchVertices"]}
    return {"Actual": sorted(protected_attributes.values()), "PlaneSupportsMatch": True,
            "Comparison": "normalized-boundary-segment-sets-and-component-hole-topology",
            "TopologyMatches": all(topology_matches),
            "MaximumRelativeMeasureError": max(area_errors, default=0.0),
            "MaximumSupportVertexDistance": max(distances, default=0.0),
            "PatchCount": len(right), "ByPlaneSupport": by_patch}


def _expected_response_owner_attributes(contract):
    cut_roles = set(contract["CutSurfaceRoles"])
    return sorted(item["Attribute"] for item in contract["BoundaryLabels"]
                  if item["Role"] not in cut_roles)


def _ownership_report(path, quadrature_path, contract):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"attribute", "elements", "ambiguous_fraction", "unresolved_elements",
                "unresolved_fraction", "quadrature_rule", "quadrature_order",
                "quadrature_points", "quadrature_whole_measure",
                "quadrature_owned_measure", "quadrature_relative_closure",
                "quadrature_closure_tolerance", "quadrature_unmatched",
                "quadrature_overlaps", "quadrature_positive_weights"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Ownership partition report is incomplete")
    fields = ("quadrature_rule", "quadrature_order", "quadrature_points",
              "quadrature_whole_measure", "quadrature_owned_measure",
              "quadrature_relative_closure", "quadrature_closure_tolerance",
              "quadrature_unmatched", "quadrature_overlaps",
              "quadrature_positive_weights")
    if any(row[name] != rows[0][name] for row in rows[1:] for name in fields):
        raise ValueError("Response-ownership quadrature summary is inconsistent")
    summary = rows[0]
    closure = float(summary["quadrature_relative_closure"])
    tolerance = float(summary["quadrature_closure_tolerance"])
    unmatched = int(summary["quadrature_unmatched"])
    overlaps = int(summary["quadrature_overlaps"])
    positive = bool(int(summary["quadrature_positive_weights"]))
    with Path(quadrature_path).open(newline="") as stream:
        quadrature = list(csv.DictReader(stream))
    if (not quadrature or set(quadrature[0]) != {"attribute", "measure"} or
            any(set(row) != {"attribute", "measure"} for row in quadrature)):
        raise ValueError("Quadrature ownership partition is incomplete")
    attributes = [int(row["attribute"]) for row in quadrature]
    measures = [float(row["measure"]) for row in quadrature]
    expected = _expected_response_owner_attributes(contract)
    if (attributes != expected or len(attributes) != len(set(attributes)) or
            any(not math.isfinite(value) or value <= 0 for value in measures)):
        raise ValueError("Quadrature owner labels/measures differ from semantic contract")
    quadrature_owned = math.fsum(measures)
    whole = float(summary["quadrature_whole_measure"])
    owned = float(summary["quadrature_owned_measure"])
    partition_closure = abs(quadrature_owned - whole) / max(abs(whole), 1e-300)
    if (not math.isfinite(whole) or whole <= 0 or not math.isfinite(owned) or owned <= 0 or
            abs(quadrature_owned - owned) > 1e-12 * max(abs(owned), 1e-300) or
            partition_closure > 1e-12):
        raise ValueError("Quadrature owner partition does not close")
    return {
        "PhysicalSurfaceCoverage": {
            "InterfaceElements": sum(int(row["elements"]) for row in rows),
            "PartitionRows": len(rows), "Complete": True},
        "ResponseOwnership": {
            "UnmatchedPolicy": "Error", "QuadratureRule": summary["quadrature_rule"],
            "QuadratureOrder": int(summary["quadrature_order"]),
            "PositiveWeights": positive, "EvaluatedPoints": int(summary["quadrature_points"]),
            "UnmatchedPoints": unmatched, "OverlappingPoints": overlaps,
            "WholeMeasure": float(summary["quadrature_whole_measure"]),
            "OwnedMeasure": float(summary["quadrature_owned_measure"]),
            "RelativeClosureError": closure, "ClosureTolerance": tolerance,
            "ExpectedOwnerAttributes": expected, "OwnerAttributes": attributes,
            "OwnerMeasures": measures, "OwnerPartitionRelativeClosure": partition_closure,
            "NoDuplicateOrMissingOwners": attributes == expected,
            "Exhaustive": (positive and unmatched == 0 and overlaps == 0 and
                           closure <= tolerance and partition_closure <= 1e-12)},
        "WholeElementAmbiguityDiagnostics": {
            "AmbiguousRows": sum(float(row["ambiguous_fraction"]) > 0 for row in rows),
            "UnresolvedElements": sum(int(row["unresolved_elements"]) for row in rows),
            "AuthoritativeForResponseOwnership": False}}


def _signature_segments(path):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"Px", "Py", "Pz", "Tx", "Ty", "Tz", "S0", "S1"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Signature lacks segment geometry")
    segments = []
    for row in rows:
        point = np.array([float(row[name]) for name in ("Px", "Py", "Pz")])
        tangent = np.array([float(row[name]) for name in ("Tx", "Ty", "Tz")])
        ends = [point + float(row[name]) * tangent for name in ("S0", "S1")]
        segments.append(np.concatenate(ends).tolist())
    return segments


def topology_record(base, mesh_path, contract_path, recipe_path, process_path,
                    signature_path, reference_mesh_path, ownership_report_path,
                    ownership_quadrature_path, corner_tolerance=1e-8):
    mesh = read_mesh(mesh_path)
    contract = load_semantic_contract(contract_path)
    report, _ = analyze(mesh, contract, require_material_names=True)
    points = np.asarray(mesh.points)
    matrix = np.asarray(base["Transform"], dtype=float).reshape(4, 4)
    corners = np.asarray(contract["SemanticCorners"], dtype=float)
    transformed_corners = (np.column_stack((corners, np.ones(len(corners)))) @
                           matrix.T)[:, :3]
    for corner in transformed_corners:
        if np.linalg.norm(points - corner, axis=1).min() > corner_tolerance:
            raise ValueError("Semantic corner is absent from the audited mesh")
    recipe = json.loads(Path(recipe_path).read_text())
    process = tomllib.loads(Path(process_path).read_text())
    transformed_recipe = dict(recipe)
    normal = recipe["NormalSizeOverThickness"] * process["MetalThickness"]
    transformed_recipe["NormalSize"] = normal
    transformed_recipe["TangentialSize"] = normal * recipe["TangentialSizeOverNormalSize"]
    transformed_recipe["TruePhysicalCorners"] = contract["SemanticCorners"]
    segments = np.asarray(_signature_segments(signature_path), dtype=float).reshape(-1, 2, 3)
    homogeneous = np.concatenate((segments, np.ones((*segments.shape[:2], 1))), axis=2)
    transformed_recipe["PhysicalSegments"] = (homogeneous @ matrix.T)[..., :3].reshape(-1, 6).tolist()
    corners = np.asarray(transformed_recipe["TruePhysicalCorners"], dtype=float).reshape(-1, 3)
    if len(corners):
        homogeneous_corners = np.column_stack((corners, np.ones(len(corners))))
        transformed_recipe["TruePhysicalCorners"] = (homogeneous_corners @ matrix.T)[:, :3].tolist()
    widths = directional_widths(mesh, transformed_recipe)
    samples = next((value for value in widths.values() if value["Cells"]), None)
    if samples is None:
        raise ValueError("No directional-width samples")
    percentiles = samples["WidthsTangentialTransverse1Transverse2"]
    _, material_attributes = blocks(mesh, "tetra")
    _, boundary_attributes = blocks(mesh, "triangle")
    actual_boundary_attributes = sorted(int(value)
                                        for value in np.unique(boundary_attributes))
    topology = contract["FeatureTopology"]
    subdivision_points = np.asarray(topology["CADSubdivisionEndpoints"], dtype=float).reshape(-1, 3)
    cut_points = np.asarray(topology["CutEndpoints"], dtype=float).reshape(-1, 3)
    def transformed(points):
        if not len(points):
            return points
        return (np.column_stack((points, np.ones(len(points)))) @ matrix.T)[:, :3]
    reference_mesh = read_mesh(reference_mesh_path)
    transformed_reference = meshio.Mesh(
        (np.column_stack((reference_mesh.points, np.ones(len(reference_mesh.points)))) @
         matrix.T)[:, :3],
        [(cell.type, cell.data.copy()) for cell in reference_mesh.cells],
        point_data=reference_mesh.point_data, cell_data=reference_mesh.cell_data,
        field_data=reference_mesh.field_data)
    actual = {"ActualVolumeMaterials": [
                  {"Attribute": int(value),
                   "Material": report["PhysicalVolumeNames"][int(value)]}
                  for value in sorted(np.unique(material_attributes))],
              "ActualBoundaryAttributes": actual_boundary_attributes,
              "ActualAdjacency": {str(value): report["BoundaryAdjacency"][value]
                                  for value in actual_boundary_attributes},
              "OwnershipClosure": _ownership_report(
                  ownership_report_path, ownership_quadrature_path, contract),
              "ActualSemanticCorners": transformed_corners.tolist(),
              "CornerNeighborhoods": _point_aspects(mesh, transformed_corners),
              "SubdivisionNeighborhoods": _point_aspects(mesh, transformed(subdivision_points)),
              "CutNeighborhoods": _point_aspects(mesh, transformed(cut_points)),
              "ProtectedSurfaces": _protected_surface_report(
                  transformed_reference, mesh, contract),
              "AchievedAnisotropy": {"Samples": samples["Cells"],
                                     "NormalTarget": transformed_recipe["NormalSize"],
                                     "TangentialP50": percentiles[1][0],
                                     "Transverse1P90": percentiles[2][1],
                                     "Transverse2P90": percentiles[2][2]},
              "TraceDiagonal": _global_diagonal_bands(
                  mesh, transformed_recipe["PhysicalSegments"],
                  transformed_recipe["NormalSize"]),
              "MeshQuality": _tetra_quality(mesh)}
    base["Measurements"] = actual
    base["Topology"] = report
    base["ReferenceMeshSHA256"] = sha256(reference_mesh_path)
    base["OwnershipReportSHA256"] = sha256(ownership_report_path)
    base["OwnershipQuadratureSHA256"] = sha256(ownership_quadrature_path)
    return base


def _simplex_h1_dofs(mesh, order):
    tetrahedra, _ = blocks(mesh, "tetra")
    edges = np.sort(tetrahedra[:, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))]
                    .reshape(-1, 2), axis=1)
    faces = np.sort(tetrahedra[:, ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))]
                    .reshape(-1, 3), axis=1)
    vertices = len(np.unique(tetrahedra))
    edge_count, face_count = len(np.unique(edges, axis=0)), len(np.unique(faces, axis=0))
    return (vertices + max(order - 1, 0) * edge_count +
            max((order - 1) * (order - 2) // 2, 0) * face_count +
            max((order - 1) * (order - 2) * (order - 3) // 6, 0) * len(tetrahedra))


def complexity_record(base, mesh_path, contract_path, recipe_path):
    mesh = read_mesh(mesh_path)
    recipe = json.loads(Path(recipe_path).read_text())
    order = int(recipe["GeometryOrder"])
    topology = load_semantic_contract(contract_path)["FeatureTopology"]
    base["Measurements"] = {"Complexity": {
        "H1DOFs": int(_simplex_h1_dofs(mesh, order)),
        "FeatureCount": topology["PhysicalFeatureCount"],
        "CADSubdivisionCount": topology["CADSubdivisionCount"]}}
    return base


def _mesh_invariants(mesh):
    tetrahedra, materials = blocks(mesh, "tetra")
    xyz = mesh.points[tetrahedra]
    volume = np.abs(np.linalg.det(np.stack((xyz[:, 1] - xyz[:, 0],
        xyz[:, 2] - xyz[:, 0], xyz[:, 3] - xyz[:, 0]), axis=2))) / 6
    triangles, labels = blocks(mesh, "triangle")
    face = mesh.points[triangles]
    area = np.linalg.norm(np.cross(face[:, 1] - face[:, 0],
                                   face[:, 2] - face[:, 0]), axis=1) / 2
    values = {f"Volume:{int(attr)}": float(volume[materials == attr].sum())
              for attr in np.unique(materials)}
    values.update({f"Area:{int(attr)}": float(area[labels == attr].sum())
                   for attr in np.unique(labels)})
    return values


def invariants_record(base, mesh_path):
    base["Measurements"] = {"ComparisonInvariants": _mesh_invariants(read_mesh(mesh_path))}
    return base


def _inverse_transformed_mesh(mesh, matrix):
    result = meshio.Mesh(mesh.points.copy(), [(cell.type, cell.data.copy()) for cell in mesh.cells],
                         point_data=mesh.point_data, cell_data=mesh.cell_data,
                         field_data=mesh.field_data)
    result.points = (np.asarray(mesh.points) - matrix[:3, 3]) @ matrix[:3, :3]
    return result


def _physical_covariance_report(identity, transformed, contract, matrix):
    normalized = _inverse_transformed_mesh(transformed, matrix)
    left, _ = analyze(identity, contract, require_material_names=True)
    right, _ = analyze(normalized, contract, require_material_names=True)
    left_invariants = _mesh_invariants(identity)
    right_invariants = _mesh_invariants(normalized)
    left_quality, right_quality = _tetra_quality(identity), _tetra_quality(normalized)
    deterministic = {"PointCountEqual": len(identity.points) == len(normalized.points)}
    deterministic["TopologyEqual"] = all(
        np.array_equal(blocks(identity, kind)[0], blocks(normalized, kind)[0]) and
        np.array_equal(blocks(identity, kind)[1], blocks(normalized, kind)[1])
        for kind in ("triangle", "tetra"))
    deterministic["CoordinateMaximumError"] = None
    if deterministic["PointCountEqual"]:
        deterministic["CoordinateMaximumError"] = float(np.max(
            np.linalg.norm(identity.points - normalized.points, axis=1)))
    return {
        "ComparisonFrame": "SourceLocal",
        "LabelsMaterialsAdjacencyMatch": (
            left["PhysicalVolumeNames"] == right["PhysicalVolumeNames"] and
            left["BoundaryAdjacency"] == right["BoundaryAdjacency"] and
            set(blocks(identity, "triangle")[1]) == set(blocks(normalized, "triangle")[1]) and
            set(blocks(identity, "tetra")[1]) == set(blocks(normalized, "tetra")[1])),
        "ProtectedSurfaces": _protected_surface_report(identity, normalized, contract),
        "ReferenceInvariants": left_invariants,
        "TransformedInvariants": right_invariants,
        "ReferenceQuality": left_quality,
        "TransformedQuality": right_quality,
        "ReferencePoints": len(identity.points), "TransformedPoints": len(normalized.points),
        "ReferenceElements": len(blocks(identity, "tetra")[0]),
        "TransformedElements": len(blocks(normalized, "tetra")[0]),
        # Diagnostic only: physical acceptance never depends on this subsection.
        "DeterministicTopologyDiagnostic": deterministic,
    }


def variant_record(base, mesh_path, identity_mesh_path, identity_seed_mesh_path,
                   transform, contract_path, stage_reports, tolerance=1e-10):
    mesh = read_mesh(mesh_path)
    matrix = np.asarray(transform, dtype=float).reshape(4, 4)
    reports, _ = validate_stage_dag(stage_reports, mesh_path)
    publication = reports["proper-rigid-publication"]
    receipt = json.loads(Path(publication["Artifacts"]["transform-receipt"]["Path"]).read_text())
    canonical_path = Path(publication["Inputs"]["canonical-candidate-mesh"]["Path"])
    canonical = read_mesh(canonical_path)
    homogeneous = np.column_stack((canonical.points, np.ones(len(canonical.points))))
    expected = (homogeneous @ matrix.T)[:, :3]
    if mesh.points.shape != expected.shape:
        raise ValueError("Rigidly published point count differs from canonical candidate")
    error = float(np.max(np.linalg.norm(mesh.points - expected, axis=1)))
    if error > tolerance or receipt.get("MaximumCoordinateError", math.inf) > tolerance:
        raise ValueError("Final candidate does not apply the declared proper rigid transform")
    for kind in ("triangle", "tetra"):
        cells, refs = blocks(mesh, kind)
        canonical_cells, canonical_refs = blocks(canonical, kind)
        if not np.array_equal(cells, canonical_cells) or not np.array_equal(refs, canonical_refs):
            raise ValueError("Rigid publication changed exact connectivity or labels")
    contract = load_semantic_contract(contract_path)
    base["IdentityMeshPath"] = str(canonical_path.resolve())
    base["IdentityMeshSHA256"] = sha256(canonical_path)
    # Retain the field name for normalized-schema compatibility; the canonical
    # source-local seed is shared rather than rebuilt per placement.
    canonical_seed = Path(reports["seed-generation"]["Artifacts"]["seed-mesh"]["Path"])
    base["IdentitySeedMeshPath"] = str(canonical_seed.resolve())
    base["IdentitySeedMeshSHA256"] = sha256(canonical_seed)
    base["TransformMaximumCoordinateError"] = error
    base["TransformVerified"] = True
    base["TransformReceiptSHA256"] = sha256(
        publication["Artifacts"]["transform-receipt"]["Path"])
    base["Measurements"] = {"PhysicalCovariance": _physical_covariance_report(
        canonical, mesh, contract, matrix)}
    return base


def bounded_record(base, mesh_path, stage_reports):
    reports, digests = validate_stage_dag(stage_reports, mesh_path)
    canonical = [reports[name] for name in CANONICAL_STAGE_ORDER]
    placement = [reports[name] for name in PLACEMENT_STAGE_ORDER]
    canonical_resources = {
        "Seconds": sum(report["Seconds"] for report in canonical),
        "PeakRSSGiB": max(report["PeakProcessTreeRSSBytes"] for report in canonical) / 2**30}
    placement_resources = {
        "Seconds": sum(report["Seconds"] for report in placement),
        "PeakRSSGiB": max(report["PeakProcessTreeRSSBytes"] for report in placement) / 2**30}
    base["Measurements"] = {"Resources": {
        "ExitCode": 0,
        "Seconds": canonical_resources["Seconds"] + placement_resources["Seconds"],
        "PeakRSSGiB": max(canonical_resources["PeakRSSGiB"],
                           placement_resources["PeakRSSGiB"]),
        "Elements": len(blocks(read_mesh(mesh_path), "tetra")[0]),
        "CanonicalBuild": canonical_resources,
        "PlacementPublication": placement_resources}}
    base["BoundedStages"] = reports
    base["BoundedStageRecords"] = [
        {"Stage": stage, "Path": str(Path(stage_reports[stage]).resolve()),
         "SHA256": stage_sha256(stage_reports[stage])} for stage in STAGE_ORDER]
    base["StageRecordSHA256"] = sorted(digests)
    return base


def produce(kind, case, variant, mesh, inputs_path, transform_path, output, *,
            contract=None, recipe=None, process=None, signature=None, identity_mesh=None,
            identity_seed_mesh=None, stage_reports=None, command=None):
    if kind not in KINDS:
        raise ValueError("Unknown audit kind")
    inputs = json.loads(Path(inputs_path).read_text())
    transform = json.loads(Path(transform_path).read_text())
    base = _base(kind, case, variant, Path(mesh), inputs, transform,
                 command or [Path(__file__).name, kind])
    if kind == "bounded-run":
        record = bounded_record(base, mesh, stage_reports)
    elif kind == "mesh-topology-quality":
        reports, _ = validate_stage_dag(stage_reports, mesh)
        reference = reports["seed-generation"]["Artifacts"]["seed-mesh"]["Path"]
        publication = reports["proper-rigid-publication"]["Artifacts"]
        ownership = publication["ownership-partition"]["Path"]
        quadrature = publication["ownership-quadrature-partition"]["Path"]
        record = topology_record(base, mesh, contract, recipe, process, signature,
                                 reference, ownership, quadrature)
        record["Dependencies"] = {"SemanticContract": sha256(contract),
                                  "MeshRecipe": sha256(recipe), "Process": sha256(process),
                                  "Signature": sha256(signature)}
    elif kind == "mesh-complexity":
        record = complexity_record(base, mesh, contract, recipe)
        record["Dependencies"] = {"MeshRecipe": sha256(recipe),
                                  "SemanticContract": sha256(contract)}
    elif kind == "mesh-invariants":
        record = invariants_record(base, mesh)
    else:
        record = variant_record(base, mesh, identity_mesh, identity_seed_mesh,
                                transform, contract, stage_reports)
    output = Path(output)
    if output.exists():
        raise ValueError("Audit output must be fresh")
    output.write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=KINDS)
    parser.add_argument("case"); parser.add_argument("variant")
    parser.add_argument("mesh", type=Path); parser.add_argument("inputs", type=Path)
    parser.add_argument("transform", type=Path); parser.add_argument("output", type=Path)
    parser.add_argument("--contract", type=Path); parser.add_argument("--recipe", type=Path)
    parser.add_argument("--process", type=Path); parser.add_argument("--signature", type=Path)
    parser.add_argument("--identity-mesh", type=Path)
    parser.add_argument("--identity-seed-mesh", type=Path)
    parser.add_argument("--stage-report", action="append", default=[], metavar="STAGE=PATH")
    args = parser.parse_args()
    stage_reports = {}
    for value in args.stage_report:
        if "=" not in value:
            parser.error("--stage-report must be STAGE=PATH")
        stage, value = value.split("=", 1)
        if stage in stage_reports:
            parser.error("duplicate stage report")
        stage_reports[stage] = Path(value)
    produce(args.kind, args.case, args.variant, args.mesh, args.inputs, args.transform,
            args.output, contract=args.contract, recipe=args.recipe, process=args.process,
            signature=args.signature, identity_mesh=args.identity_mesh,
            identity_seed_mesh=args.identity_seed_mesh,
            stage_reports=stage_reports, command=sys.argv)


if __name__ == "__main__":
    main()
