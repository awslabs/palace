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

import numpy as np

from audit_edge_metric_mesh import analyze, blocks, directional_widths
from general_mesh_manifest import canonical_sha256, sha256
from mesh_array_io import read_mesh
from mesh_stage_contract import STAGE_ORDER, sha256 as stage_sha256, validate_stage_dag
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


def _global_diagonal_bands(mesh, physical_segments):
    triangles, labels = blocks(mesh, "triangle")
    xyz = mesh.points[triangles]
    normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    owners_by_edge = {}
    all_lengths = []
    for owner, triangle in enumerate(triangles):
        for a, b in ((triangle[0], triangle[1]), (triangle[1], triangle[2]),
                     (triangle[2], triangle[0])):
            edge = tuple(sorted((int(a), int(b))))
            owners_by_edge.setdefault(edge, []).append(owner)
            all_lengths.append(np.linalg.norm(mesh.points[edge[0]] - mesh.points[edge[1]]))
    median = float(np.median(all_lengths))
    short = []
    for edge, owners in owners_by_edge.items():
        length = np.linalg.norm(mesh.points[edge[0]] - mesh.points[edge[1]])
        if (len(owners) == 2 and labels[owners[0]] == labels[owners[1]] and
                abs(np.dot(normals[owners[0]], normals[owners[1]])) > 1 - 1e-10 and
                length < .6 * median):
            short.append(edge)
    adjacency = {}
    for first, last in short:
        adjacency.setdefault(first, set()).add(last)
        adjacency.setdefault(last, set()).add(first)
    segments = np.asarray(physical_segments, dtype=float).reshape(-1, 2, 3)
    bands = 0
    visited = set()
    surface_diameter = max(float(np.linalg.norm(np.ptp(mesh.points[triangles], axis=(0, 1)))),
                           1e-300)
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
        span = float(np.linalg.norm(np.ptp(points, axis=0)))
        if span <= .5 * surface_diameter:
            continue
        direction = points[np.argmax(np.linalg.norm(points - points[0], axis=1))] - points[0]
        direction /= np.linalg.norm(direction)
        aligned = False
        for segment in segments:
            tangent = segment[1] - segment[0]
            tangent /= np.linalg.norm(tangent)
            if abs(np.dot(direction, tangent)) > 1 - 1e-6:
                aligned = True
                break
        bands += int(not aligned)
    return {"GlobalDiagonalBands": bands, "ShortInternalEdges": len(short),
            "MedianSurfaceEdgeLength": median}


def _tetra_quality(mesh):
    tetrahedra, _ = blocks(mesh, "tetra")
    xyz = mesh.points[tetrahedra]
    jacobian = np.stack((xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0],
                         xyz[:, 3] - xyz[:, 0]), axis=2)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    condition = singular[:, 0] / singular[:, -1]
    lengths = np.linalg.norm(jacobian, axis=1)
    scaled = np.abs(np.linalg.det(jacobian)) / np.prod(lengths, axis=1)
    if (not len(scaled) or not np.all(np.isfinite(condition)) or
            np.any(singular[:, -1] <= 0)):
        raise ValueError("Invalid tetrahedral Jacobian audit")
    return {"Samples": len(tetrahedra),
            "MinimumScaledJacobian": float(scaled.min()),
            "MaximumJacobianCondition": float(condition.max())}


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


def _protected_surface_report(reference, candidate, contract):
    before, _ = analyze(reference, contract, require_material_names=True)
    after, _ = analyze(candidate, contract, require_material_names=True)
    protected_attributes = {item["Attribute"]: item["Role"]
                            for item in contract["BoundaryLabels"]
                            if item.get("Protected") is True}
    def areas(report):
        result = {attribute: 0.0 for attribute in protected_attributes}
        for key, value in report["PlanarPatchAreas"].items():
            attribute = int(float(key.split()[0]))
            if attribute in result:
                result[attribute] += value
        return result
    left, right = areas(before), areas(after)
    errors = {protected_attributes[key]: abs(right[key] - value) / max(abs(value), 1e-300)
              for key, value in left.items()}
    return {"Actual": sorted(protected_attributes.values()),
            "MaximumRelativeMeasureError": max(errors.values(), default=0.0),
            "RelativeMeasureErrorByRole": errors}


def _ownership_report(path):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"attribute", "elements", "ambiguous_fraction", "unresolved_elements",
                "unresolved_fraction"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Ownership partition report is incomplete")
    unresolved = sum(int(row["unresolved_elements"]) for row in rows)
    ambiguous = sum(float(row["ambiguous_fraction"]) > 0 for row in rows)
    return {"UnmatchedPolicy": "Error", "Unmatched": unresolved,
            "Overlaps": ambiguous, "Exhaustive": unresolved == 0 and ambiguous == 0,
            "PartitionRows": len(rows)}


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
                    corner_tolerance=1e-8):
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
    actual = {"ActualVolumeMaterials": [
                  {"Attribute": int(value),
                   "Material": report["PhysicalVolumeNames"][int(value)]}
                  for value in sorted(np.unique(material_attributes))],
              "ActualBoundaryAttributes": actual_boundary_attributes,
              "ActualAdjacency": {str(value): report["BoundaryAdjacency"][value]
                                  for value in actual_boundary_attributes},
              "OwnershipClosure": _ownership_report(ownership_report_path),
              "ActualSemanticCorners": transformed_corners.tolist(),
              "CornerNeighborhoods": _point_aspects(mesh, transformed_corners),
              "SubdivisionNeighborhoods": _point_aspects(mesh, transformed(subdivision_points)),
              "CutNeighborhoods": _point_aspects(mesh, transformed(cut_points)),
              "ProtectedSurfaces": _protected_surface_report(
                  read_mesh(reference_mesh_path), mesh, contract),
              "AchievedAnisotropy": {"Samples": samples["Cells"],
                                     "NormalTarget": transformed_recipe["NormalSize"],
                                     "TangentialP50": percentiles[1][0],
                                     "Transverse1P90": percentiles[2][1],
                                     "Transverse2P90": percentiles[2][2]},
              "TraceDiagonal": _global_diagonal_bands(
                  mesh, transformed_recipe["PhysicalSegments"]),
              "MeshQuality": _tetra_quality(mesh)}
    base["Measurements"] = actual
    base["Topology"] = report
    base["ReferenceMeshSHA256"] = sha256(reference_mesh_path)
    base["OwnershipReportSHA256"] = sha256(ownership_report_path)
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


def invariants_record(base, mesh_path):
    mesh = read_mesh(mesh_path)
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
    base["Measurements"] = {"ComparisonInvariants": values}
    return base


def variant_record(base, mesh_path, identity_mesh_path, transform, tolerance=1e-10):
    mesh, identity = read_mesh(mesh_path), read_mesh(identity_mesh_path)
    matrix = np.asarray(transform, dtype=float).reshape(4, 4)
    homogeneous = np.column_stack((identity.points, np.ones(len(identity.points))))
    expected = (homogeneous @ matrix.T)[:, :3]
    if mesh.points.shape != expected.shape:
        raise ValueError("Variant and identity point counts differ")
    error = float(np.max(np.linalg.norm(mesh.points - expected, axis=1)))
    if error > tolerance:
        raise ValueError("Candidate mesh does not apply the declared transform")
    for kind in ("triangle", "tetra"):
        cells, refs = blocks(mesh, kind)
        identity_cells, identity_refs = blocks(identity, kind)
        if not np.array_equal(cells, identity_cells) or not np.array_equal(refs, identity_refs):
            raise ValueError("Variant topology/labels differ from identity")
    base["IdentityMeshSHA256"] = sha256(identity_mesh_path)
    base["TransformMaximumCoordinateError"] = error
    base["TransformVerified"] = True
    base["Measurements"] = {}
    return base


def bounded_record(base, mesh_path, stage_reports):
    reports, digests = validate_stage_dag(stage_reports, mesh_path)
    base["Measurements"] = {"Resources": {
        "ExitCode": 0,
        "Seconds": sum(report["Seconds"] for report in reports.values()),
        "PeakRSSGiB": max(report["PeakProcessTreeRSSBytes"] for report in reports.values()) / 2**30,
        "Elements": len(blocks(read_mesh(mesh_path), "tetra")[0])}}
    base["BoundedStages"] = reports
    base["BoundedStageRecords"] = [
        {"Stage": stage, "Path": str(Path(stage_reports[stage]).resolve()),
         "SHA256": stage_sha256(stage_reports[stage])} for stage in STAGE_ORDER]
    base["StageRecordSHA256"] = sorted(digests)
    return base


def produce(kind, case, variant, mesh, inputs_path, transform_path, output, *,
            contract=None, recipe=None, process=None, signature=None, identity_mesh=None,
            stage_reports=None, command=None):
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
        ownership = reports["final-gmsh-publication"]["Artifacts"][
            "ownership-partition"]["Path"]
        record = topology_record(base, mesh, contract, recipe, process, signature,
                                 reference, ownership)
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
        record = variant_record(base, mesh, identity_mesh, transform)
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
            stage_reports=stage_reports, command=sys.argv)


if __name__ == "__main__":
    main()
