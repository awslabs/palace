#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Bind producer measurements and artifacts to one manifest matrix entry.

This tool does not manufacture expected values.  Expected labels, materials,
corners, protected supports, and adjacency remain exclusively in the frozen
semantic contract and are read only by the gate.
"""
import argparse
import json
from pathlib import Path

from general_mesh_manifest import canonical_sha256, sha256, validate_manifest


def _binding(path, output):
    path = path.resolve()
    try:
        display = str(path.relative_to(output.parent.resolve()))
    except ValueError:
        display = str(path)
    return {"Path": display, "SHA256": sha256(path)}


def normalize(manifest_path, case_id, variant_id, raw_path, mesh_path, audit_paths, output):
    manifest_path, raw_path, mesh_path, output = map(Path,
        (manifest_path, raw_path, mesh_path, output))
    manifest = json.loads(manifest_path.read_text())
    repository, tools, matrix = validate_manifest(manifest, manifest_path.resolve())
    if (case_id, variant_id) not in matrix:
        raise ValueError("Case/variant is not in the required manifest matrix")
    case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
    variant = next(item for item in case["Variants"] if item["Id"] == variant_id)
    source = case["Source"]
    directory = Path(source["Directory"])
    directory = directory if directory.is_absolute() else repository / directory
    inputs = {role: item["SHA256"] for role, item in source["Files"].items()}
    for role, item in source["Files"].items():
        repository_name = item.get("RepositoryPath")
        path = Path(repository_name or item["Name"])
        if not path.is_absolute():
            path = (repository if repository_name else directory) / path
        if not path.is_file() or sha256(path) != item["SHA256"]:
            raise ValueError(f"Immutable {role} input changed")
    raw = json.loads(raw_path.read_text())
    if raw.get("CaseId") != case_id or raw.get("Variant") != variant_id:
        raise ValueError("Producer record names the wrong case or variant")
    required_actual = ("ActualVolumeMaterials", "ActualBoundaryAttributes", "ActualAdjacency",
                       "OwnershipClosure", "ActualSemanticCorners", "ProtectedSurfaces", "AchievedAnisotropy",
                       "TraceDiagonal", "Resources", "MeshQuality", "Complexity",
                       "ComparisonInvariants")
    if any(name not in raw for name in required_actual):
        raise ValueError("Producer record is incomplete")
    if output.exists():
        raise ValueError("Evidence output must be fresh")
    audit_paths = [Path(path) for path in audit_paths]
    if not audit_paths:
        raise ValueError("At least one independent producer audit is required")
    evidence = {"Version": 2, "CaseId": case_id, "Variant": variant_id,
                "TransformSHA256": canonical_sha256(variant["Transform"]),
                "InputSHA256": inputs, "ProcessSHA256": inputs["Process"],
                "SemanticContractSHA256": inputs["SemanticContract"],
                "RecipeSHA256": inputs["MeshRecipe"], "ToolSHA256": tools,
                "Mesh": _binding(mesh_path, output), "AuditRecords": []}
    for path in audit_paths:
        item = _binding(path, output)
        item["Kind"] = "producer-audit"
        evidence["AuditRecords"].append(item)
    evidence.update({name: raw[name] for name in required_actual})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(evidence, indent=2) + "\n")
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("case")
    parser.add_argument("variant")
    parser.add_argument("raw", type=Path)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit-record", type=Path, action="append", required=True)
    args = parser.parse_args()
    normalize(args.manifest, args.case, args.variant, args.raw, args.mesh,
              args.audit_record, args.output)


if __name__ == "__main__":
    main()
