#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Transform a coupon's coordinate-bearing source and semantic contracts.

The existing CSV plan-view formats are intentionally local-frame formats.  This
producer companion emits a 3-D support record for arbitrary proper rigid
transforms while preserving every non-coordinate label and material field.  The
transformed semantic contract is consumed directly by metric preparation.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def _matrix_product(left, right):
    return [[sum(left[i][k] * right[k][j] for k in range(3)) for j in range(3)]
            for i in range(3)]


def validate_rigid_transform(values, tolerance=1e-12):
    if (not isinstance(values, list) or len(values) != 16 or
            any(not isinstance(value, (int, float)) or not math.isfinite(value)
                for value in values)):
        raise ValueError("Rigid transform must contain 16 finite row-major values")
    matrix = [[float(values[4 * row + column]) for column in range(4)]
              for row in range(4)]
    if any(abs(matrix[3][column] - expected) > tolerance
           for column, expected in enumerate((0.0, 0.0, 0.0, 1.0))):
        raise ValueError("Rigid transform must have homogeneous last row [0, 0, 0, 1]")
    rotation = [row[:3] for row in matrix[:3]]
    transpose = [[rotation[j][i] for j in range(3)] for i in range(3)]
    gram = _matrix_product(transpose, rotation)
    if any(abs(gram[i][j] - (1.0 if i == j else 0.0)) > tolerance
           for i in range(3) for j in range(3)):
        raise ValueError("Rigid transform rotation must be orthogonal")
    determinant = (
        rotation[0][0] * (rotation[1][1] * rotation[2][2] - rotation[1][2] * rotation[2][1])
        - rotation[0][1] * (rotation[1][0] * rotation[2][2] - rotation[1][2] * rotation[2][0])
        + rotation[0][2] * (rotation[1][0] * rotation[2][1] - rotation[1][1] * rotation[2][0]))
    if abs(determinant - 1.0) > tolerance:
        raise ValueError("Rigid transform must preserve orientation")
    return matrix


def read_transform(path):
    data = json.loads(Path(path).read_text())
    if isinstance(data, dict):
        data = data.get("Transform")
    return validate_rigid_transform(data)


def transform_point(matrix, point):
    return [sum(matrix[row][column] * point[column] for column in range(3)) +
            matrix[row][3] for row in range(3)]


def transform_vector(matrix, vector):
    return [sum(matrix[row][column] * vector[column] for column in range(3))
            for row in range(3)]


def transform_semantic_contract(contract, matrix):
    result = json.loads(json.dumps(contract))
    result["SemanticCorners"] = [transform_point(matrix, point)
                                   for point in result["SemanticCorners"]]
    topology = result["FeatureTopology"]
    for name in ("CADSubdivisionEndpoints", "CutEndpoints"):
        topology[name] = [transform_point(matrix, point) for point in topology[name]]
    result["RigidTransform"] = [value for row in matrix for value in row]
    return result


def _read_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def transformed_supports(source, matrix, *, signature=None, boundary=None, mask=None):
    source = Path(source)
    signature = Path(signature or source / "mesh-signature.csv")
    boundary = Path(boundary or source / "plan-view-boundary.csv")
    mask = Path(mask or source / "plan-view-mask.csv")
    edges = []
    for row in _read_rows(signature):
        point = [float(row[name]) for name in ("Px", "Py", "Pz")]
        gap = [float(row[name]) for name in ("Gx", "Gy", "Gz")]
        tangent = [float(row[name]) for name in ("Tx", "Ty", "Tz")]
        normal = [0.0, 0.0, float(row["Nz"])]
        edges.append({
            "Index": int(row["Index"]), "Slot": int(row["Slot"]),
            "Conductor": int(row["Conductor"]),
            "Point": transform_point(matrix, point),
            "GapDirection": transform_vector(matrix, gap),
            "TangentDirection": transform_vector(matrix, tangent),
            "ProcessNormalDirection": transform_vector(matrix, normal),
            "Interval": [float(row["S0"]), float(row["S1"])],
            "VertexArm": int(row.get("VertexArm", 0)),
        })
    boundary = [{
        "Loop": int(row["Loop"]), "Vertex": int(row["Vertex"]),
        "Conductor": int(row["Conductor"]), "Hole": int(row["Hole"]),
        "Class": row["Class"],
        "Point": transform_point(matrix, [float(row["X"]), float(row["Y"]),
                                            float(row["Plane"])]),
    } for row in _read_rows(boundary)]
    mask = [{
        "Facet": int(row["Facet"]), "Conductor": int(row["Conductor"]),
        "Point": transform_point(matrix, [float(row["X"]), float(row["Y"]),
                                            float(row["Plane"])]),
    } for row in _read_rows(mask)]
    return {"Version": 1, "CoordinateSystem": "TransformedGlobal3D",
            "RigidTransform": [value for row in matrix for value in row],
            "Edges": edges, "BoundaryVertices": boundary, "MaskVertices": mask}


def write_transformed_contracts(source, transform_path, semantic_output, supports_output, *,
                                semantic_input=None, signature=None, boundary=None,
                                mask=None):
    source = Path(source)
    matrix = read_transform(transform_path)
    semantic_input = Path(semantic_input or source / "semantic-contract.json")
    contract = transform_semantic_contract(json.loads(semantic_input.read_text()), matrix)
    contract["SourceSemanticContractSHA256"] = hashlib.sha256(
        semantic_input.read_bytes()).hexdigest()
    Path(semantic_output).write_text(json.dumps(contract, indent=2) + "\n")
    Path(supports_output).write_text(json.dumps(transformed_supports(
        source, matrix, signature=signature, boundary=boundary, mask=mask), indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("transform", type=Path)
    parser.add_argument("semantic_output", type=Path)
    parser.add_argument("supports_output", type=Path)
    parser.add_argument("--semantic-input", type=Path)
    parser.add_argument("--signature", type=Path)
    parser.add_argument("--boundary", type=Path)
    parser.add_argument("--mask", type=Path)
    args = parser.parse_args()
    for output in (args.semantic_output, args.supports_output):
        if output.exists():
            parser.error(f"refuse to overwrite output: {output}")
    write_transformed_contracts(
        args.source, args.transform, args.semantic_output, args.supports_output,
        semantic_input=args.semantic_input, signature=args.signature,
        boundary=args.boundary, mask=args.mask)


if __name__ == "__main__":
    main()
