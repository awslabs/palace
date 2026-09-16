#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Independent per-entry verification of one manifest case.

Every manifest variant of the case is judged with the unchanged production
functions (manifest validation, immutable-input hashes, evidence gates, mesh
readability, bound audit/stage records, exact canonical-build reuse) and the
frozen physical-covariance comparison is evaluated whenever both compared
variants have evidence - even when a variant fails - so one report lists every
failure of the case.  The report is written whether or not the case passed;
the exit status is nonzero unless everything passed.  Scope: a single case and
its manifest variants; no matrix, physics or release qualification.
"""
import argparse
import csv
import json
from pathlib import Path
import signal
import sys

from canonical_mesh_build import same_canonical_build
from general_mesh_manifest import (_check_artifact, _finite_number,
                                   _physical_comparison_failures, _validate_bound_records,
                                   _validate_mesh, audit_manifest_evidence, canonical_sha256,
                                   sha256, validate_manifest)
from semantic_mesh_contract import load_semantic_contract, validate_feature_topology

PRODUCTION_FUNCTIONS = ["validate_manifest", "validate_feature_topology",
                        "audit_manifest_evidence", "_check_artifact", "_validate_mesh",
                        "_validate_bound_records", "same_canonical_build",
                        "_physical_comparison_failures"]
_EXPECTED_ERRORS = (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError)


def _immutable_inputs(manifest_path, repository, case):
    """Hash-check every immutable input of the case; returns (hashes, paths)."""
    directory = Path(case["Source"]["Directory"])
    directory = directory if directory.is_absolute() else repository / directory
    hashes, paths = {}, {}
    for role, entry in case["Source"]["Files"].items():
        repository_name = entry.get("RepositoryPath")
        path = Path(repository_name or entry["Name"])
        if not path.is_absolute():
            path = (repository if repository_name else directory) / path
        if not path.is_file() or sha256(path) != entry["SHA256"]:
            raise ValueError(f"immutable {role} hash mismatch")
        hashes[role], paths[role] = entry["SHA256"], path
    contract = load_semantic_contract(paths["SemanticContract"])
    validate_feature_topology(contract, paths["Signature"], paths["Boundary"])
    with paths[case["Source"]["SignatureRole"]].open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or not set(case["Source"]["SignatureColumns"]).issubset(rows[0]):
        raise ValueError("signature preflight failed")
    return hashes, paths, contract


def verify_variant(case, variant, evidence_path, evidence, contract, hashes, paths, tools,
                   manifest, shared):
    """Judge one variant; `shared` accumulates the cross-variant uniqueness and
    canonical-reuse state (used meshes, used record digests, canonical builds).
    Returns the entry record with its gate failures; raises on binding errors."""
    binding = {"CaseId": case["Id"], "Variant": variant["Id"], "Transform": variant["Transform"],
               "TransformSHA256": canonical_sha256(variant["Transform"]),
               "InputSHA256": hashes, "ToolSHA256": tools,
               "StageToolSHA256": manifest["StageToolSHA256"], "Gates": manifest["Gates"]}
    failures = audit_manifest_evidence(evidence, manifest["Gates"], contract, binding)
    mesh_path = _check_artifact(evidence_path.parent, evidence.get("Mesh"), "audited mesh")
    _validate_mesh(mesh_path, contract)
    if evidence["Mesh"]["SHA256"] in shared["meshes"]:
        raise ValueError("audited meshes must be content-distinct per matrix entry")
    shared["meshes"].add(evidence["Mesh"]["SHA256"])
    variant_digests, canonical_digests, canonical_record = _validate_bound_records(
        evidence_path, evidence, binding, paths)
    if shared["variant_digests"] & variant_digests:
        raise ValueError("variant audit/placement records must be content-distinct")
    shared["variant_digests"].update(variant_digests)
    build_id = canonical_record["CanonicalBuildId"]
    previous = shared["canonical"].get(build_id)
    if previous is None:
        if any(canonical_digests & item[1] for item in shared["canonical"].values()):
            raise ValueError("canonical stages were reused under a different cache key")
        shared["canonical"][build_id] = (canonical_record, canonical_digests)
    elif (not same_canonical_build(previous[0], canonical_record) or
          previous[1] != canonical_digests):
        raise ValueError("shared canonical stages require exact cache key and hashes")
    protected = evidence["ProtectedSurfaces"]
    return {"Evidence": {"Path": str(evidence_path), "SHA256": sha256(evidence_path)},
            "Mesh": {"Path": str(mesh_path), "SHA256": sha256(mesh_path)},
            "GateFailures": failures, "Error": None, "Passed": not failures,
            "CanonicalBuildId": build_id,
            "CanonicalBuildSHA256": canonical_record["CanonicalBuildSHA256"],
            "CanonicalStageReportSHA256": sorted(canonical_digests),
            "VariantRecordSHA256": sorted(variant_digests),
            "Resources": evidence["Resources"], "MeshQuality": evidence["MeshQuality"],
            "ProtectedSurfaces": {key: protected.get(key) for key in
                                  ("MaximumRelativeMeasureError",
                                   "MaximumSupportVertexDistance", "PatchCount",
                                   "CoplanarSlotUnions")},
            "CornerNeighborhoods": evidence["CornerNeighborhoods"],
            "OwnershipClosure": evidence["OwnershipClosure"]["ResponseOwnership"],
            "TraceDiagonal": {key: value for key, value in evidence["TraceDiagonal"].items()
                              if key != "LongShortEdgeComponents"}}


def covariance_failures(manifest, case, evidence_by_variant):
    """The frozen physical-covariance comparison of the case, or why it could not run."""
    comparison = case["TransformComparison"]
    reference = evidence_by_variant.get(comparison["Reference"])
    transformed = evidence_by_variant.get(comparison["Transformed"])
    if reference is None or transformed is None:
        return ["missing transform evidence"], None
    identity_digest = reference["Mesh"]["SHA256"]
    coordinate_error = transformed.get("TransformMaximumCoordinateError")
    failures = []
    if (reference.get("IdentityMeshSHA256") != identity_digest or
            transformed.get("IdentityMeshSHA256") != identity_digest or
            transformed.get("IdentitySeedMeshSHA256") != reference.get("IdentitySeedMeshSHA256") or
            not _finite_number(coordinate_error, nonnegative=True) or
            coordinate_error > manifest["Gates"]["CornerTolerance"]):
        failures.append("exact source-seed covariance")
    failures.extend(_physical_comparison_failures(reference, transformed, comparison))
    return failures, coordinate_error


def verify_case(manifest_path, audit_root, case_id):
    """Verify every variant of `case_id` and the case covariance; never stops at the
    first failure.  Returns the report (Passed is true only with no failure)."""
    manifest_path, audit_root = Path(manifest_path).resolve(), Path(audit_root).resolve()
    manifest = json.loads(manifest_path.read_text())
    repository, tools, matrix = validate_manifest(manifest, manifest_path)
    case = next((item for item in manifest["Cases"] if item["Id"] == case_id), None)
    if case is None:
        raise ValueError(f"case {case_id} is not in the manifest")
    hashes, paths, contract = _immutable_inputs(manifest_path, repository, case)
    shared = {"meshes": set(), "variant_digests": set(), "canonical": {}}
    entries, evidence_by_variant, failures = {}, {}, []
    for variant in case["Variants"]:
        variant_id = variant["Id"]
        evidence_path = audit_root / f"{case_id}--{variant_id}.json"
        entry = {"Evidence": {"Path": str(evidence_path)}, "GateFailures": [], "Error": None,
                 "Passed": False}
        try:
            if (case_id, variant_id) not in matrix:
                raise ValueError("required entry is absent from manifest")
            evidence = json.loads(evidence_path.read_text())
            evidence_by_variant[variant_id] = evidence
            entry = verify_variant(case, variant, evidence_path, evidence, contract, hashes,
                                   paths, tools, manifest, shared)
        except _EXPECTED_ERRORS as error:
            entry["Error"] = str(error)
        entries[variant_id] = entry
        failures.extend(f"{variant_id}: gate failure: {name}" for name in entry["GateFailures"])
        if entry["Error"] is not None:
            failures.append(f"{variant_id}: {entry['Error']}")
    comparison_failures, coordinate_error = covariance_failures(manifest, case,
                                                                evidence_by_variant)
    failures.extend(f"covariance: {name}" for name in comparison_failures)
    reuse_failures = []
    if len(shared["canonical"]) != 1:
        reuse_failures.append("placements did not share exactly one canonical build")
    failures.extend(f"canonical reuse: {name}" for name in reuse_failures)
    return {"Version": 3,
            "Scope": "single case, all manifest variants and their covariance comparison; "
                     "every variant is judged and every failure reported; no matrix, "
                     "physics, or release qualification",
            "CaseId": case_id, "Passed": not failures, "Failures": failures,
            "Manifest": {"Path": str(manifest_path), "SHA256": sha256(manifest_path)},
            "InputSHA256": hashes, "ToolSHA256": tools,
            "StageToolSHA256": manifest["StageToolSHA256"],
            "Entries": entries,
            "SharedCanonicalBuildId": (next(iter(shared["canonical"]))
                                       if len(shared["canonical"]) == 1 else None),
            "CanonicalBuildIds": sorted(shared["canonical"]),
            "CanonicalReuseFailures": reuse_failures,
            "TransformComparisonFailures": comparison_failures,
            "TransformMaximumCoordinateError": coordinate_error,
            "ProductionFunctions": PRODUCTION_FUNCTIONS}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("audit_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("case_id")
    parser.add_argument("--timeout-seconds", type=int, default=1800,
                        help="wall-clock bound on the whole verification (default 1800)")
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError("verification output must be fresh")
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(
        TimeoutError(f"verification exceeded {args.timeout_seconds} seconds")))
    signal.alarm(args.timeout_seconds)
    report = verify_case(args.manifest, args.audit_root, args.case_id)
    signal.alarm(0)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in
                      ("Passed", "Failures", "SharedCanonicalBuildId",
                       "TransformComparisonFailures")}, indent=1))
    return 0 if report["Passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
