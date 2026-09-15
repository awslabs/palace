#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Publish a fresh Gmsh 2.2 mesh through a bound proper rigid transform."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess

import meshio
import numpy as np

from general_mesh_audit_producer import _mesh_invariants, _ownership_report, _tetra_quality
from transform_coupon_source_contract import (read_transform, transform_semantic_contract,
                                                transformed_supports)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _line(data, start):
    end = data.find(b"\n", start)
    if end < 0:
        raise ValueError("Truncated Gmsh section")
    raw = data[start:end]
    if raw.endswith(b"\r"):
        raw = raw[:-1]
    return raw, end + 1


def _node_section(data):
    marker = b"$Nodes"
    starts = []
    offset = 0
    while True:
        found = data.find(marker, offset)
        if found < 0:
            break
        if (found == 0 or data[found - 1:found] == b"\n"):
            starts.append(found)
        offset = found + len(marker)
    if len(starts) != 1:
        raise ValueError("Gmsh 2.2 mesh must have exactly one $Nodes section")
    start = starts[0]
    marker_line, cursor = _line(data, start)
    if marker_line != marker:
        raise ValueError("Malformed $Nodes marker")
    count_line, records = _line(data, cursor)
    try:
        count = int(count_line)
    except ValueError as error:
        raise ValueError("Invalid Gmsh node count") from error
    mesh_format = data.find(b"$MeshFormat")
    if mesh_format < 0:
        raise ValueError("Missing Gmsh mesh format")
    _, format_cursor = _line(data, mesh_format)
    format_line, format_cursor = _line(data, format_cursor)
    fields = format_line.split()
    if len(fields) != 3 or fields[0] != b"2.2" or fields[2] != b"8":
        raise ValueError("Rigid publication requires Gmsh 2.2 with 8-byte reals")
    binary = fields[1] == b"1"
    if fields[1] not in (b"0", b"1"):
        raise ValueError("Invalid Gmsh binary flag")
    if binary:
        endian_line, after_endian = _line(data, format_cursor + 4)
        endian_bytes = data[format_cursor:format_cursor + 4]
        if endian_line:
            raise ValueError("Malformed binary Gmsh endian marker")
        if struct.unpack("<i", endian_bytes)[0] == 1:
            endian = "<"
        elif struct.unpack(">i", endian_bytes)[0] == 1:
            endian = ">"
        else:
            raise ValueError("Invalid binary Gmsh endian marker")
        width = struct.calcsize(endian + "i3d")
        end = records + count * width
        if end > len(data):
            raise ValueError("Truncated binary node records")
        newline, after = _line(data, end)
        if newline:
            raise ValueError("Binary node records are not newline terminated")
        end_marker, finish = _line(data, after)
        if end_marker != b"$EndNodes":
            raise ValueError("Missing $EndNodes marker")
        tags, points = [], []
        for index in range(count):
            tag, x, y, z = struct.unpack_from(endian + "i3d", data, records + index * width)
            tags.append(tag); points.append((x, y, z))
        return {"binary": True, "endian": endian, "start": records, "end": end,
                "finish": finish, "count": count, "tags": tags,
                "points": np.asarray(points, dtype=float)}
    tags, points, lines = [], [], []
    cursor = records
    for _ in range(count):
        raw, following = _line(data, cursor)
        fields = raw.split()
        if len(fields) != 4:
            raise ValueError("ASCII Gmsh node record must contain tag and three coordinates")
        tags.append(int(fields[0])); points.append([float(value) for value in fields[1:]])
        lines.append((cursor, following)); cursor = following
    end_marker, finish = _line(data, cursor)
    if end_marker != b"$EndNodes":
        raise ValueError("Missing $EndNodes marker")
    return {"binary": False, "start": records, "end": cursor, "finish": finish,
            "count": count, "tags": tags, "points": np.asarray(points, dtype=float),
            "lines": lines}


def transform_gmsh22(source, output, matrix):
    source, output = Path(source).resolve(), Path(output).resolve()
    if source == output or not source.is_file() or output.exists():
        raise ValueError("Rigid publication requires distinct existing input and fresh output")
    data = source.read_bytes()
    section = _node_section(data)
    points = section["points"]
    transformed = points @ np.asarray(matrix, dtype=float)[:3, :3].T + np.asarray(matrix)[:3, 3]
    if not np.all(np.isfinite(transformed)):
        raise ValueError("Rigid transform produced nonfinite coordinates")
    if section["binary"]:
        records = bytearray(data[section["start"]:section["end"]])
        width = struct.calcsize(section["endian"] + "i3d")
        for index, (tag, point) in enumerate(zip(section["tags"], transformed)):
            struct.pack_into(section["endian"] + "i3d", records, index * width,
                             tag, *map(float, point))
        published = data[:section["start"]] + bytes(records) + data[section["end"]:]
    else:
        newline = b"\r\n" if b"\r\n" in data[section["start"]:section["end"]] else b"\n"
        records = b"".join((f"{tag} " + " ".join(format(value, '.17g') for value in point))
                           .encode() + newline
                           for tag, point in zip(section["tags"], transformed))
        published = data[:section["start"]] + records + data[section["end"]:]
    output.write_bytes(published)
    actual = _node_section(published)
    error = float(np.max(np.linalg.norm(actual["points"] - transformed, axis=1)))
    return section, actual, error


def _equal_data(left, right):
    if left.keys() != right.keys():
        return False
    return all(np.array_equal(np.asarray(left[key]), np.asarray(right[key])) for key in left)


def _exact_mesh_structure(left, right):
    return (len(left.cells) == len(right.cells) and
            all(a.type == b.type and np.array_equal(a.data, b.data)
                for a, b in zip(left.cells, right.cells)) and
            _equal_data(left.point_data, right.point_data) and
            left.cell_data.keys() == right.cell_data.keys() and
            all(len(left.cell_data[key]) == len(right.cell_data[key]) and
                all(np.array_equal(a, b) for a, b in zip(left.cell_data[key], right.cell_data[key]))
                for key in left.cell_data) and
            _equal_data(left.field_data, right.field_data))


def publish(canonical_mesh, transform_path, output_mesh, receipt_path, *,
            semantic_input, signature, boundary, mask, process, canonical_build_record,
            transformed_semantic, transformed_supports_path, ownership, ownership_quadrature,
            ownership_runtime, ownership_auditor, kind="fabricated", tolerance=1e-12):
    outputs = [Path(value) for value in (output_mesh, receipt_path, transformed_semantic,
                                         transformed_supports_path, ownership,
                                         ownership_quadrature)]
    if any(path.exists() for path in outputs):
        raise ValueError("Every rigid-publication output must be fresh")
    source_inputs = {"source-semantic-contract": Path(semantic_input),
                     "source-signature": Path(signature), "source-boundary": Path(boundary),
                     "source-mask": Path(mask), "source-process": Path(process)}
    if any(not path.is_file() for path in source_inputs.values()):
        raise ValueError("Every explicit source input must exist")
    matrix = read_transform(transform_path)
    canonical_build = json.loads(Path(canonical_build_record).read_text())
    canonical_digest = sha256(canonical_mesh)
    artifacts = canonical_build.get("CanonicalArtifacts", {})
    bound_candidate = artifacts.get("canonical-candidate-mesh", {})
    if bound_candidate.get("SHA256") != canonical_digest:
        raise ValueError("Canonical candidate differs from the bound canonical build")
    before_section, after_section, coordinate_error = transform_gmsh22(
        canonical_mesh, output_mesh, matrix)
    if before_section["tags"] != after_section["tags"]:
        raise ValueError("Rigid publication changed node tags")
    left, right = meshio.read(canonical_mesh), meshio.read(output_mesh)
    if not _exact_mesh_structure(left, right):
        raise ValueError("Rigid publication changed cell blocks, connectivity, labels, or data")
    expected = left.points @ np.asarray(matrix)[:3, :3].T + np.asarray(matrix)[:3, 3]
    coordinate_error = max(coordinate_error,
                           float(np.max(np.linalg.norm(right.points - expected, axis=1))))
    if coordinate_error > tolerance:
        raise ValueError("Rigid publication coordinate error exceeds the frozen tolerance")
    identity = np.array_equal(np.asarray(matrix), np.eye(4))
    output_digest = sha256(output_mesh)
    if not identity and output_digest == canonical_digest:
        raise ValueError("Nonidentity rigid publication reused canonical mesh bytes")
    quality = _tetra_quality(right)
    if quality["PositiveOrientation"] is not True:
        raise ValueError("Rigid publication changed tetrahedron orientation")
    before_invariants, after_invariants = _mesh_invariants(left), _mesh_invariants(right)
    measure_error = max((abs(after_invariants[key] - value) / max(abs(value), 1e-300)
                         for key, value in before_invariants.items()), default=0.0)
    if measure_error > 1e-11:
        raise ValueError("Rigid publication changed a material or physical measure")

    transform_digest = sha256(transform_path)
    semantic_digest = sha256(semantic_input)
    semantic = transform_semantic_contract(json.loads(Path(semantic_input).read_text()), matrix)
    semantic["SourceSemanticContractSHA256"] = semantic_digest
    semantic["CanonicalTransformSHA256"] = transform_digest
    supports = transformed_supports(
        Path(signature).parent, matrix, signature=signature, boundary=boundary, mask=mask,
        transform_sha256=transform_digest, semantic_sha256=semantic_digest)
    Path(transformed_semantic).write_text(json.dumps(semantic, indent=2) + "\n")
    Path(transformed_supports_path).write_text(json.dumps(supports, indent=2) + "\n")

    transform_csv = Path(str(receipt_path) + ".transform.csv")
    if transform_csv.exists():
        raise ValueError("Rigid-publication private transform path must be fresh")
    transform_csv.write_text(",".join(format(value, ".17g")
                                      for row in matrix for value in row) + "\n")
    command = [str(Path(ownership_runtime).resolve()), "--startup-file=no",
               str(Path(ownership_auditor).resolve()), kind,
               str(Path(output_mesh).resolve()), str(transform_csv.resolve()),
               str(Path(ownership).resolve()),
               "--process", str(Path(process).resolve()),
               "--signature", str(Path(signature).resolve()),
               "--boundary", str(Path(boundary).resolve())]
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"Transformed ownership audit failed with status {result.returncode}")
    generated_quadrature = Path(str(ownership) + ".quadrature.csv")
    if generated_quadrature.resolve() != Path(ownership_quadrature).resolve():
        if Path(ownership_quadrature).exists() or not generated_quadrature.is_file():
            raise ValueError("Ownership quadrature output path is inconsistent")
        generated_quadrature.replace(ownership_quadrature)
    contract = json.loads(Path(semantic_input).read_text())
    ownership_summary = _ownership_report(ownership, ownership_quadrature, contract)
    if ownership_summary["ResponseOwnership"]["Exhaustive"] is not True:
        raise ValueError("Transformed response ownership is not exhaustive")
    receipt = {
        "Version": 2,
        "CanonicalBuildId": canonical_build.get("CanonicalBuildId"),
        "CanonicalBuildSHA256": canonical_build.get("CanonicalBuildSHA256"),
        "CanonicalMeshSHA256": canonical_digest,
        "OutputMeshSHA256": output_digest,
        "TransformSHA256": transform_digest,
        "SourceInputSHA256": {name: sha256(path) for name, path in source_inputs.items()},
        "Transform": [value for row in matrix for value in row],
        "Identity": bool(identity),
        "MaximumCoordinateError": coordinate_error,
        "CoordinateTolerance": tolerance,
        "NodeCount": len(before_section["tags"]),
        "NodeTagsExact": True,
        "CellBlocksConnectivityLabelsAndDataExact": True,
        "PositiveOrientation": True,
        "MaximumRelativeMeasureError": measure_error,
        "MeshQuality": quality,
        "TransformedSemanticSHA256": sha256(transformed_semantic),
        "TransformedSupportsSHA256": sha256(transformed_supports_path),
        "OwnershipSHA256": sha256(ownership),
        "OwnershipQuadratureSHA256": sha256(ownership_quadrature),
        "OwnershipCommand": command,
        "OwnershipRuntimeSHA256": sha256(ownership_runtime),
        "OwnershipAuditorSHA256": sha256(ownership_auditor),
    }
    receipt["ReceiptPayloadSHA256"] = hashlib.sha256(json.dumps(
        receipt, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    Path(receipt_path).write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("canonical_mesh", type=Path)
    parser.add_argument("transform", type=Path)
    parser.add_argument("output_mesh", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--semantic-input", type=Path, required=True)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--boundary", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--process", type=Path, required=True)
    parser.add_argument("--canonical-build-record", type=Path, required=True)
    parser.add_argument("--transformed-semantic", type=Path, required=True)
    parser.add_argument("--transformed-supports", type=Path, required=True)
    parser.add_argument("--ownership", type=Path, required=True)
    parser.add_argument("--ownership-quadrature", type=Path, required=True)
    parser.add_argument("--ownership-runtime", type=Path, required=True)
    parser.add_argument("--ownership-auditor", type=Path, required=True)
    parser.add_argument("--kind", choices=("thin", "fabricated"), default="fabricated")
    parser.add_argument("--tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    try:
        publish(args.canonical_mesh, args.transform, args.output_mesh, args.receipt,
                semantic_input=args.semantic_input,
                signature=args.signature, boundary=args.boundary, mask=args.mask,
                process=args.process,
                canonical_build_record=args.canonical_build_record,
                transformed_semantic=args.transformed_semantic,
                transformed_supports_path=args.transformed_supports,
                ownership=args.ownership, ownership_quadrature=args.ownership_quadrature,
                ownership_runtime=args.ownership_runtime,
                ownership_auditor=args.ownership_auditor, kind=args.kind,
                tolerance=args.tolerance)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
