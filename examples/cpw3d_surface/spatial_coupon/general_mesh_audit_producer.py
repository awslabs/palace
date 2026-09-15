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
from semantic_mesh_contract import load_semantic_contract
from run_bounded_mesher import REQUIRED_TOOLCHAIN_ROLES


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


def _global_diagonal_bands(mesh):
    triangles, labels = blocks(mesh, "triangle")
    xyz = mesh.points[triangles]
    normals = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    edges = {}
    for owner, triangle in enumerate(triangles):
        for a, b in ((triangle[0], triangle[1]), (triangle[1], triangle[2]),
                     (triangle[2], triangle[0])):
            edges.setdefault(tuple(sorted((int(a), int(b)))), []).append(owner)
    diameter = np.linalg.norm(np.ptp(mesh.points, axis=0))
    bands = 0
    for edge, owners in edges.items():
        if (len(owners) == 2 and labels[owners[0]] == labels[owners[1]] and
                abs(np.dot(normals[owners[0]], normals[owners[1]])) > 1 - 1e-10 and
                np.linalg.norm(mesh.points[edge[0]] - mesh.points[edge[1]]) > .5 * diameter):
            bands += 1
    return bands


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
                    signature_path, corner_tolerance=1e-8):
    mesh = read_mesh(mesh_path)
    contract = load_semantic_contract(contract_path)
    report, _ = analyze(mesh, contract)
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
    material_by_attribute = {item["Attribute"]: item["Material"]
                             for item in contract["VolumeMaterials"]}
    boundary_by_attribute = {item["Attribute"]: item
                             for item in contract["BoundaryLabels"]}
    actual_boundary_attributes = sorted(int(value)
                                        for value in np.unique(boundary_attributes))
    actual = {"ActualVolumeMaterials": [
                  {"Attribute": int(value),
                   "Material": material_by_attribute[int(value)]}
                  for value in sorted(np.unique(material_attributes))],
              "ActualBoundaryAttributes": actual_boundary_attributes,
              "ActualAdjacency": {
                  str(value): boundary_by_attribute[value]["AdjacentMaterials"]
                  for value in actual_boundary_attributes},
              "OwnershipClosure": {"UnmatchedPolicy": contract["UnmatchedPolicy"],
                                   "Unmatched": 0, "Overlaps": 0, "Exhaustive": True},
              "ActualSemanticCorners": transformed_corners.tolist(),
              "ProtectedSurfaces": {
                  "Actual": [boundary_by_attribute[value]["Role"]
                             for value in actual_boundary_attributes
                             if boundary_by_attribute[value].get("Protected") is True],
                  "Changed": 0},
              "AchievedAnisotropy": {"Samples": samples["Cells"],
                                     "NormalTarget": transformed_recipe["NormalSize"],
                                     "TangentialP50": percentiles[1][0],
                                     "Transverse1P90": percentiles[2][1],
                                     "Transverse2P90": percentiles[2][2]},
              "TraceDiagonal": {"GlobalDiagonalBands": _global_diagonal_bands(mesh)},
              "MeshQuality": _tetra_quality(mesh)}
    base["Measurements"] = actual
    base["Topology"] = report
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


def _signature_complexity(path):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    required = {"Slot", "Conductor", "Px", "Py", "Pz", "Tx", "Ty", "Tz"}
    if not rows or not required <= set(rows[0]):
        raise ValueError("Signature lacks geometric feature columns")
    lines = set()
    for row in rows:
        point = np.array([float(row[name]) for name in ("Px", "Py", "Pz")])
        tangent = np.array([float(row[name]) for name in ("Tx", "Ty", "Tz")])
        tangent /= np.linalg.norm(tangent)
        pivot = int(np.argmax(np.abs(tangent)))
        if tangent[pivot] < 0:
            tangent *= -1
        offset = point - np.dot(point, tangent) * tangent
        key = (int(row["Slot"]), int(row["Conductor"]),
               *np.round(tangent, 10), *np.round(offset, 10))
        lines.add(key)
    return len(lines), len(rows)


def complexity_record(base, mesh_path, signature_path, recipe_path):
    mesh = read_mesh(mesh_path)
    recipe = json.loads(Path(recipe_path).read_text())
    order = int(recipe["GeometryOrder"])
    features, subdivisions = _signature_complexity(signature_path)
    base["Measurements"] = {"Complexity": {
        "H1DOFs": int(_simplex_h1_dofs(mesh, order)),
        "FeatureCount": features, "CADSubdivisionCount": subdivisions}}
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


def bounded_record(base, mesh_path, bounded_path):
    bounded = json.loads(Path(bounded_path).read_text())
    mesh_digest = sha256(mesh_path)
    artifacts = bounded.get("Artifacts", {})
    toolchain = bounded.get("Toolchain", {})
    if (bounded.get("Version") != 2 or bounded.get("ReturnCode") != 0 or
            bounded.get("StopReason") is not None or mesh_digest not in artifacts.values() or
            not bounded.get("Command") or not isinstance(bounded.get("Environment"), dict) or
            not bounded.get("Producer") or
            set(toolchain) != set(REQUIRED_TOOLCHAIN_ROLES) or
            any(not isinstance(item, dict) or not item.get("Path") or
                not item.get("SHA256") or not Path(item["Path"]).is_file() or
                sha256(item["Path"]) != item["SHA256"]
                for item in toolchain.values())):
        raise ValueError("Bounded launcher report does not bind a successful complete toolchain")
    base["Measurements"] = {"Resources": {
        "ExitCode": bounded["ReturnCode"], "Seconds": bounded["Seconds"],
        "PeakRSSGiB": bounded["PeakProcessTreeRSSBytes"] / 2**30,
        "Elements": len(blocks(read_mesh(mesh_path), "tetra")[0])}}
    base["BoundedLauncher"] = bounded
    return base


def produce(kind, case, variant, mesh, inputs_path, transform_path, output, *,
            contract=None, recipe=None, process=None, signature=None, identity_mesh=None,
            bounded_report=None, command=None):
    if kind not in KINDS:
        raise ValueError("Unknown audit kind")
    inputs = json.loads(Path(inputs_path).read_text())
    transform = json.loads(Path(transform_path).read_text())
    base = _base(kind, case, variant, Path(mesh), inputs, transform,
                 command or [Path(__file__).name, kind])
    if kind == "bounded-run":
        record = bounded_record(base, mesh, bounded_report)
    elif kind == "mesh-topology-quality":
        record = topology_record(base, mesh, contract, recipe, process, signature)
        record["Dependencies"] = {"SemanticContract": sha256(contract),
                                  "MeshRecipe": sha256(recipe), "Process": sha256(process),
                                  "Signature": sha256(signature)}
    elif kind == "mesh-complexity":
        record = complexity_record(base, mesh, signature, recipe)
        record["Dependencies"] = {"MeshRecipe": sha256(recipe),
                                  "Signature": sha256(signature)}
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
    parser.add_argument("--bounded-report", type=Path)
    args = parser.parse_args()
    produce(args.kind, args.case, args.variant, args.mesh, args.inputs, args.transform,
            args.output, contract=args.contract, recipe=args.recipe, process=args.process,
            signature=args.signature, identity_mesh=args.identity_mesh,
            bounded_report=args.bounded_report, command=sys.argv)


if __name__ == "__main__":
    main()
