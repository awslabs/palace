#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Immutable identity and cache contract for one source-local canonical mesh build."""
import hashlib
import json
from pathlib import Path

from mesh_stage_contract import (CANONICAL_STAGE_ORDER, LEGACY_MMG_PIPELINE,
                                 canonical_artifact_roles, canonical_stage_order,
                                 canonical_tool_roles, pipeline_of, pipeline_of_tool_roles,
                                 sha256, validate_stage_report)


CANONICAL_SOURCE_ROLES = ("Signature", "Boundary", "Mask", "Process",
                          "SemanticContract", "MeshRecipe")
# Exact artifact and tool schemas are derived from the pipeline's canonical stages so
# the build hash binds every canonical output and every runtime-resolved tool.  The
# pipeline is identified by the exact canonical tool roles of the cache key (the
# Gmsh-only production pipeline's gmsh-build/mesher, or the legacy MMG stages); the
# module constants name the legacy schema.
CANONICAL_ARTIFACT_ROLES = canonical_artifact_roles(LEGACY_MMG_PIPELINE)
CANONICAL_TOOL_ROLES = canonical_tool_roles(LEGACY_MMG_PIPELINE)


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _digests(value, description, expected_names=None):
    if (not isinstance(value, dict) or not value or
            any(not isinstance(name, str) or not name or
                not isinstance(digest, str) or len(digest) != 64
                for name, digest in value.items())):
        raise ValueError(f"{description} must contain named SHA-256 digests")
    if expected_names is not None and set(value) != set(expected_names):
        raise ValueError(f"{description} names differ from the exact canonical schema")
    return dict(sorted(value.items()))


def cache_key_pipeline(canonical_tool_sha256):
    """The pipeline whose exact canonical tool roles a cache key (or tool dict) names."""
    if not isinstance(canonical_tool_sha256, dict):
        raise ValueError("Canonical tools must contain named SHA-256 digests")
    return pipeline_of_tool_roles(canonical_tool_sha256)


def canonical_cache_key(source_sha256, gates, canonical_tool_sha256):
    sources = _digests(source_sha256, "Canonical sources")
    if not set(CANONICAL_SOURCE_ROLES).issubset(sources):
        raise ValueError("Canonical sources omit a required immutable source role")
    pipeline = cache_key_pipeline(canonical_tool_sha256)
    tools = _digests(canonical_tool_sha256, "Canonical tools", canonical_tool_roles(pipeline))
    if not isinstance(gates, dict) or not gates:
        raise ValueError("Canonical qualification gates must be frozen")
    # Process, semantic and recipe hashes remain named explicitly even though all
    # immutable source roles are included. This makes cache-review failures legible.
    return {
        "Version": 2,
        "SourceSHA256": sources,
        "ProcessSHA256": sources["Process"],
        "SemanticContractSHA256": sources["SemanticContract"],
        "RecipeSHA256": sources["MeshRecipe"],
        "QualificationGates": gates,
        "CanonicalToolSHA256": tools,
    }


def build_record(source_sha256, gates, canonical_tool_sha256, artifacts,
                 stage_report_sha256):
    key = canonical_cache_key(source_sha256, gates, canonical_tool_sha256)
    pipeline = cache_key_pipeline(canonical_tool_sha256)
    if not isinstance(artifacts, dict) or set(artifacts) != canonical_artifact_roles(pipeline):
        raise ValueError("Canonical artifacts differ from the exact canonical-stage output schema")
    artifact_bindings = {}
    for name, item in artifacts.items():
        if (not isinstance(item, dict) or not item.get("Path") or
                not isinstance(item.get("SHA256"), str) or len(item["SHA256"]) != 64):
            raise ValueError("Canonical artifacts require paths and SHA-256 digests")
        artifact_bindings[name] = {"Path": str(Path(item["Path"]).resolve()),
                                   "SHA256": item["SHA256"]}
    record = {
        "Version": 2,
        "CanonicalBuildId": canonical_sha256(key),
        "CanonicalCacheKey": key,
        "CanonicalArtifacts": dict(sorted(artifact_bindings.items())),
        "CanonicalStageReportSHA256": _digests(stage_report_sha256, "Canonical stage reports",
                                               canonical_stage_order(pipeline)),
    }
    record["CanonicalBuildSHA256"] = canonical_sha256(record)
    return record


def build_record_from_stage_reports(source_sha256, gates, canonical_tool_sha256,
                                    report_paths):
    """Derive every canonical artifact and report hash from validated stage reports."""
    pipeline = pipeline_of(report_paths, canonical_only=True)
    if cache_key_pipeline(canonical_tool_sha256) != pipeline:
        raise ValueError("Canonical tool roles differ from the stage reports' pipeline")
    artifacts, report_sha256 = {}, {}
    for stage in canonical_stage_order(pipeline):
        path = Path(report_paths[stage])
        report = validate_stage_report(json.loads(path.read_text()), stage, pipeline=pipeline)
        for role, item in report["Tools"].items():
            if canonical_tool_sha256.get(f"{stage}/{role}") != item["SHA256"]:
                raise ValueError(f"Canonical stage tool differs from the frozen tool: {stage}/{role}")
        artifacts.update(report["Artifacts"])
        report_sha256[stage] = sha256(path)
    return build_record(source_sha256, gates, canonical_tool_sha256, artifacts, report_sha256)


def validate_build_record(record, source_sha256, gates, canonical_tool_sha256,
                          *, check_files=True):
    expected_key = canonical_cache_key(source_sha256, gates, canonical_tool_sha256)
    if (record.get("Version") != 2 or
            record.get("CanonicalCacheKey") != expected_key or
            record.get("CanonicalBuildId") != canonical_sha256(expected_key)):
        raise ValueError("Canonical build cache key differs from source, gates, or tools")
    payload = dict(record)
    claimed = payload.pop("CanonicalBuildSHA256", None)
    if claimed != canonical_sha256(payload):
        raise ValueError("Canonical build hash differs from its immutable payload")
    pipeline = cache_key_pipeline(canonical_tool_sha256)
    artifacts = record.get("CanonicalArtifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != canonical_artifact_roles(pipeline):
        raise ValueError("Canonical build artifacts differ from the exact canonical-stage output schema")
    _digests(record.get("CanonicalStageReportSHA256"), "Canonical stage reports",
             canonical_stage_order(pipeline))
    if check_files:
        for name, item in artifacts.items():
            if (not isinstance(item, dict) or not Path(item.get("Path", "")).is_file() or
                    sha256(item["Path"]) != item.get("SHA256")):
                raise ValueError(f"Canonical artifact changed: {name}")
    return record


def same_canonical_build(left, right):
    """True only for an exact cache-key, canonical-artifact and stage-report hash match."""
    return (left.get("CanonicalBuildId") == right.get("CanonicalBuildId") and
            left.get("CanonicalBuildSHA256") == right.get("CanonicalBuildSHA256") and
            left.get("CanonicalCacheKey") == right.get("CanonicalCacheKey") and
            left.get("CanonicalStageReportSHA256") == right.get("CanonicalStageReportSHA256") and
            {name: item.get("SHA256") for name, item in
             left.get("CanonicalArtifacts", {}).items()} ==
            {name: item.get("SHA256") for name, item in
             right.get("CanonicalArtifacts", {}).items()})


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_hashes", type=Path)
    parser.add_argument("gates", type=Path)
    parser.add_argument("canonical_tool_hashes", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--stage-report", action="append", required=True,
                        metavar="STAGE=PATH")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("canonical build record must be fresh")
    reports = {}
    for value in args.stage_report:
        if "=" not in value:
            parser.error("--stage-report must be STAGE=PATH")
        stage, raw_path = value.split("=", 1); path = Path(raw_path).resolve()
        if stage in reports or not path.is_file():
            parser.error("canonical stage reports must be unique and exist")
        reports[stage] = path
    try:
        record = build_record_from_stage_reports(
            json.loads(args.source_hashes.read_text()), json.loads(args.gates.read_text()),
            json.loads(args.canonical_tool_hashes.read_text()), reports)
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        parser.error(str(error))
    args.output.write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
