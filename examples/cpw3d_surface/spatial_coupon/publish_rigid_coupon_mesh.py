#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Publish a fresh Gmsh 2.2 mesh through a bound proper rigid transform, then split its
MA surfaces into the per-ring radial shells (supervisor decision 61a: every production
coupon carries the label-only shells of relabel_radial_ma_shells.py; the ring radii are the
gmsh-build census's tube ring set bound through the canonical build record, the shells are
classified in source-local coordinates through the inverse rigid map, and the census
<output>.radial-shells.json is bound by SHA-256 in the transform receipt (RadialShells)).
The transformed ownership audit runs on the shelled mesh (the auditor's family rule
re-derives the parent of every shell element and keeps its shell prefix)."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import tomllib

import meshio
import numpy as np

from general_mesh_audit_producer import _mesh_invariants, _ownership_report, _volume_quality
from mixed_mesh import parent_label_view
import relabel_radial_ma_shells as radial_shells
from transform_coupon_source_contract import (read_transform, transform_semantic_contract,
                                                transformed_supports)

RADIAL_SHELLS_SUFFIX = ".radial-shells.json"
RADIAL_SHELLS_KIND = "radial-ma-shells"


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


def _exact_mesh_structure(left, right, *, cell_data_keys=None):
    """Identical cell blocks, connectivity, point data, cell data and field data;
    `cell_data_keys` restricts the cell data compared and skips the field data (the
    radial-shell comparison: the physical labels of the parent view, not the shells'
    own elementary tags and names)."""
    keys = left.cell_data.keys() if cell_data_keys is None else cell_data_keys
    return (len(left.cells) == len(right.cells) and
            all(a.type == b.type and np.array_equal(a.data, b.data)
                for a, b in zip(left.cells, right.cells)) and
            _equal_data(left.point_data, right.point_data) and
            (cell_data_keys is not None or left.cell_data.keys() == right.cell_data.keys()) and
            all(key in right.cell_data and len(left.cell_data[key]) == len(right.cell_data[key]) and
                all(np.array_equal(a, b) for a, b in zip(left.cell_data[key], right.cell_data[key]))
                for key in keys) and
            (cell_data_keys is not None or _equal_data(left.field_data, right.field_data)))


def _bound_canonical_artifact(canonical_build, name):
    """The path of a canonical artifact bound in the canonical build record, its
    SHA-256 verified; None when the record binds no artifact of that name."""
    item = canonical_build.get("CanonicalArtifacts", {}).get(name)
    if item is None:
        return None
    path = Path(item["Path"])
    if not path.is_file() or sha256(path) != item.get("SHA256"):
        raise ValueError(f"Canonical artifact {name} is missing or differs from the bound canonical build")
    return path


def apply_radial_shells(output_mesh, census_path, *, canonical_build, matrix, signature, boundary, process,
                        parent_digest, canonical_digest):
    """Split the MA surfaces of the published (parent-labeled) mesh into the per-ring
    radial shells in place and write the census; returns the receipt record.  The ring
    radii are the gmsh-build census's PrismTubes.Section.RingRadii (bound through the
    canonical build record); a canonical build without a gmsh-build census (the legacy
    MMG pipeline) publishes no shell and records why."""
    census_path = Path(census_path)
    if census_path.exists():
        raise ValueError("Radial-shell census output must be fresh")
    build_census_path = _bound_canonical_artifact(canonical_build, "build-census")
    if build_census_path is None:
        return {"Applied": False, "Reason": "the canonical build binds no gmsh-build census (legacy pipeline): "
                                            "no tube ring set, the MA surfaces keep their parent labels"}
    build_census = json.loads(build_census_path.read_text())
    radii = [float(v) for v in (build_census.get("PrismTubes") or {}).get("Section", {}).get("RingRadii", [])]
    if not radii:
        raise ValueError("The gmsh-build census records no prism-tube ring set (PrismTubes.Section.RingRadii): "
                         "the radial MA shells need the production tubes")
    partition_path = _bound_canonical_artifact(canonical_build, "canonical-ownership-partition")
    if partition_path is None:
        raise ValueError("The canonical build binds no ownership partition: the shell areas cannot be checked")
    thickness = float(tomllib.loads(Path(process).read_text())["MetalThickness"])
    lines = radial_shells.metal_edge_lines(radial_shells.read_csv_rows(boundary), radial_shells.read_csv_rows(signature),
                                           thickness)
    rotation = np.asarray(matrix, dtype=float)[:3, :3]
    translation = np.asarray(matrix, dtype=float)[:3, 3]

    def source_coordinates(points):
        return (np.asarray(points, dtype=float) - translation) @ rotation

    data = Path(output_mesh).read_bytes()
    out, mesh, shells, relabeled, ma_labels = radial_shells.relabel(data, lines=lines, radii=radii,
                                                                    source_coordinates=source_coordinates)
    label_only = radial_shells.assert_label_only(data, out, mesh, relabeled)
    parents = radial_shells.parent_area_check(shells, ma_labels, partition_path)
    closure = radial_shells.closure_check(shells)
    Path(output_mesh).write_bytes(out)
    census = radial_shells.shell_census(
        kind=RADIAL_SHELLS_KIND, case_id=None, mesh_path=output_mesh, mesh_bytes=out,
        parent_mesh={"Path": None, "SHA256": parent_digest, "CanonicalMeshSHA256": canonical_digest,
                     "Rule": "the rigidly published mesh before the shell relabel (parent MA labels); not written"},
        radii=radii, thickness=thickness, lines=lines, shells=shells, parents=parents, label_only=label_only,
        closure=closure, certificate=None,
        extra={"Placement": {"Transform": [float(v) for row in np.asarray(matrix, dtype=float) for v in row],
                             "Rule": "shell of an element = ring interval of its centroid distance to the nearest metal "
                                     "edge line, both taken in source-local coordinates through the inverse rigid map"},
               "BuildCensus": {"Path": str(build_census_path), "SHA256": sha256(build_census_path)},
               "CanonicalOwnershipPartition": {"Path": str(partition_path), "SHA256": sha256(partition_path)}})
    census_path.write_text(json.dumps(census, indent=2) + "\n")
    return {"Applied": True, "Path": str(census_path.resolve()), "SHA256": sha256(census_path), "Kind": RADIAL_SHELLS_KIND,
            "Tool": census["Tool"], "ToolSHA256": census["ToolSHA256"], "ShellCount": len(census["Shells"]),
            "MAParents": ma_labels, "RingRadii": radii, "LabelOnly": label_only,
            "StraddlingFraction": closure["StraddlingFraction"], "RelativeClosure": closure["RelativeClosure"],
            "ParentAreaMaximumRelativeDifference": max(item["RelativeDifference"] for item in parents.values()),
            "Rule": "label-only: the $Nodes block and every element apart from the (physical, elementary) pair of the MA "
                    "elements are byte-identical to the parent-labeled publication (LabelOnly); the shells sum to the "
                    "canonical ownership partition's parent areas; the ownership audit below runs on the shelled mesh"}


def publish(canonical_mesh, transform_path, output_mesh, receipt_path, *,
            semantic_input, signature, boundary, mask, process, canonical_build_record,
            transformed_semantic, transformed_supports_path, ownership, ownership_quadrature,
            ownership_runtime, ownership_auditor, kind="fabricated", tolerance=1e-12):
    radial_shells_path = Path(str(output_mesh) + RADIAL_SHELLS_SUFFIX)
    outputs = [Path(value) for value in (output_mesh, receipt_path, transformed_semantic,
                                         transformed_supports_path, ownership,
                                         ownership_quadrature, radial_shells_path)]
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
    parent_digest = sha256(output_mesh)
    if not identity and parent_digest == canonical_digest:
        raise ValueError("Nonidentity rigid publication reused canonical mesh bytes")
    quality = _volume_quality(right)
    if quality["PositiveOrientation"] is not True:
        raise ValueError("Rigid publication changed volume element orientation")
    before_invariants, after_invariants = _mesh_invariants(left), _mesh_invariants(right)
    measure_error = max((abs(after_invariants[key] - value) / max(abs(value), 1e-300)
                         for key, value in before_invariants.items()), default=0.0)
    if measure_error > 1e-11:
        raise ValueError("Rigid publication changed a material or physical measure")
    # The per-ring radial MA shells (decision 61a), label-only on the published bytes.
    shells = apply_radial_shells(output_mesh, radial_shells_path, canonical_build=canonical_build, matrix=matrix,
                                 signature=signature, boundary=boundary, process=process,
                                 parent_digest=parent_digest, canonical_digest=canonical_digest)
    output_digest = sha256(output_mesh)
    if shells["Applied"]:
        shelled = parent_label_view(meshio.read(output_mesh))
        if not _exact_mesh_structure(right, shelled, cell_data_keys=("gmsh:physical",)):
            raise ValueError("Radial shell relabel changed the mesh beyond the MA shell labels")
        if _mesh_invariants(shelled) != after_invariants:
            raise ValueError("Radial shell relabel changed a material or physical measure")

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
        "ParentLabeledMeshSHA256": parent_digest,
        "RadialShells": shells,
        "TransformSHA256": transform_digest,
        "SourceInputSHA256": {name: sha256(path) for name, path in source_inputs.items()},
        "Transform": [value for row in matrix for value in row],
        "Identity": bool(identity),
        "MaximumCoordinateError": coordinate_error,
        "CoordinateTolerance": tolerance,
        "NodeCount": len(before_section["tags"]),
        "NodeTagsExact": True,
        "CellBlocksConnectivityLabelsAndDataExact": True,
        "CellBlocksConnectivityLabelsAndDataExactRule": "the rigidly published mesh before the radial shell relabel "
                                                        "(ParentLabeledMeshSHA256) equals the canonical mesh in every "
                                                        "cell block, connectivity, label and data; the output differs "
                                                        "from it in the MA shell labels only (RadialShells.LabelOnly)",
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
