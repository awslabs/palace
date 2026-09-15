#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Assemble separately produced, mesh-bound generality audit records."""
import argparse
import json
from pathlib import Path

from general_mesh_manifest import canonical_sha256, sha256, validate_manifest
from mesh_stage_contract import STAGE_ORDER, validate_stage_dag


REQUIRED_AUDITS = ("bounded-run", "mesh-topology-quality", "mesh-complexity",
                   "mesh-invariants", "variant-transform")


def _binding(path, output, **extra):
    path = path.resolve()
    try:
        display = str(path.relative_to(output.parent.resolve()))
    except ValueError:
        display = str(path)
    return {"Path": display, "SHA256": sha256(path), **extra}


def _source_paths(repository, case):
    source = case["Source"]
    directory = Path(source["Directory"])
    directory = directory if directory.is_absolute() else repository / directory
    paths = {}
    for role, item in source["Files"].items():
        repository_name = item.get("RepositoryPath")
        path = Path(repository_name or item["Name"])
        if not path.is_absolute():
            path = (repository if repository_name else directory) / path
        if not path.is_file() or sha256(path) != item["SHA256"]:
            raise ValueError(f"Immutable {role} input changed")
        paths[role] = path
    return paths


def normalize(manifest_path, case_id, variant_id, mesh_path, audit_paths, output):
    manifest_path, mesh_path, output = map(Path, (manifest_path, mesh_path, output))
    manifest = json.loads(manifest_path.read_text())
    repository, tools, matrix = validate_manifest(manifest, manifest_path.resolve())
    if (case_id, variant_id) not in matrix:
        raise ValueError("Case/variant is not in the required manifest matrix")
    case = next(item for item in manifest["Cases"] if item["Id"] == case_id)
    variant = next(item for item in case["Variants"] if item["Id"] == variant_id)
    paths = _source_paths(repository, case)
    inputs = {role: item["SHA256"] for role, item in case["Source"]["Files"].items()}
    mesh_digest = sha256(mesh_path)
    transform_digest = canonical_sha256(variant["Transform"])
    if output.exists():
        raise ValueError("Evidence output must be fresh")
    if set(audit_paths) != set(REQUIRED_AUDITS):
        raise ValueError("Exactly one record of every required audit kind is required")

    evidence = {"Version": 3, "CaseId": case_id, "Variant": variant_id,
                "TransformSHA256": transform_digest, "InputSHA256": inputs,
                "ProcessSHA256": inputs["Process"],
                "SemanticContractSHA256": inputs["SemanticContract"],
                "RecipeSHA256": inputs["MeshRecipe"], "ToolSHA256": tools,
                "Mesh": _binding(mesh_path, output), "AuditRecords": []}
    measurements = {}
    record_hashes = {mesh_digest}
    expected_dependencies = {
        "bounded-run": {},
        "mesh-topology-quality": {
            role: inputs[role]
            for role in ("SemanticContract", "MeshRecipe", "Process", "Signature")},
        "mesh-complexity": {
            role: inputs[role] for role in ("MeshRecipe", "SemanticContract")},
        "mesh-invariants": {},
        "variant-transform": {},
    }
    records_by_kind = {}
    for kind in REQUIRED_AUDITS:
        path = Path(audit_paths[kind])
        record = json.loads(path.read_text())
        records_by_kind[kind] = record
        producer = record.get("Producer", {})
        expected = {"Version": 1, "Kind": kind, "CaseId": case_id,
                    "Variant": variant_id, "MeshSHA256": mesh_digest,
                    "InputSHA256": inputs, "TransformSHA256": transform_digest}
        if any(record.get(key) != value for key, value in expected.items()):
            raise ValueError(f"{kind} record has stale or mismatched bindings")
        if record.get("Transform") != variant["Transform"]:
            raise ValueError(f"{kind} record does not bind the exact transform")
        dependencies = record.get("Dependencies", {})
        if dependencies != expected_dependencies[kind]:
            raise ValueError(f"{kind} record dependencies differ from frozen inputs")
        if (producer.get("Name") not in tools or
                tools[producer["Name"]] != producer.get("SHA256") or
                not isinstance(record.get("Command"), list) or not record["Command"] or
                not isinstance(record.get("Environment"), dict)):
            raise ValueError(f"{kind} producer/command/environment is not frozen")
        digest = sha256(path)
        if digest in record_hashes:
            raise ValueError("Audit records must be content-distinct")
        record_hashes.add(digest)
        for section, value in record.get("Measurements", {}).items():
            if section in measurements:
                raise ValueError(f"Measurement section produced twice: {section}")
            measurements[section] = value
        evidence["AuditRecords"].append(_binding(path, output, Kind=kind))
        if kind == "variant-transform":
            if record.get("TransformVerified") is not True:
                raise ValueError("Variant transform was not verified")
            evidence["IdentityMeshSHA256"] = record.get("IdentityMeshSHA256")
            evidence["TransformMaximumCoordinateError"] = record.get(
                "TransformMaximumCoordinateError")
        elif kind == "bounded-run":
            stage_items = record.get("BoundedStageRecords")
            if (not isinstance(stage_items, list) or len(stage_items) != len(STAGE_ORDER) or
                    {item.get("Stage") for item in stage_items} != set(STAGE_ORDER)):
                raise ValueError("Bounded stage records are incomplete")
            stage_paths = {item["Stage"]: Path(item["Path"]) for item in stage_items}
            reports, stage_digests = validate_stage_dag(stage_paths, mesh_path)
            if (sorted(stage_digests) != record.get("StageRecordSHA256") or
                    reports != record.get("BoundedStages")):
                raise ValueError("Bounded stage DAG differs from its producer record")
    bounded = records_by_kind["bounded-run"]["BoundedStages"]
    topology = records_by_kind["mesh-topology-quality"]
    if (topology.get("ReferenceMeshSHA256") !=
            bounded["seed-generation"]["Artifacts"]["seed-mesh"]["SHA256"] or
            topology.get("OwnershipReportSHA256") !=
            bounded["final-gmsh-publication"]["Artifacts"]["ownership-partition"]["SHA256"]):
        raise ValueError("Topology measurements do not bind the staged reference/ownership outputs")
    required_sections = {"Resources", "ActualVolumeMaterials", "ActualBoundaryAttributes",
                         "ActualAdjacency", "OwnershipClosure", "ActualSemanticCorners",
                         "CornerNeighborhoods", "SubdivisionNeighborhoods", "CutNeighborhoods",
                         "ProtectedSurfaces", "AchievedAnisotropy", "TraceDiagonal",
                         "MeshQuality", "Complexity", "ComparisonInvariants"}
    if set(measurements) != required_sections:
        raise ValueError("Producer records do not supply the exact measurement schema")
    evidence.update(measurements)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(evidence, indent=2) + "\n")
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path); parser.add_argument("case")
    parser.add_argument("variant"); parser.add_argument("mesh", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit-record", action="append", required=True,
                        metavar="KIND=PATH")
    args = parser.parse_args()
    records = {}
    for value in args.audit_record:
        if "=" not in value:
            parser.error("--audit-record must be KIND=PATH")
        kind, path = value.split("=", 1)
        if kind in records:
            parser.error("duplicate audit kind")
        records[kind] = Path(path)
    normalize(args.manifest, args.case, args.variant, args.mesh, records, args.output)


if __name__ == "__main__":
    main()
