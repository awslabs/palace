#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Immutable identity and cache contract for one source-local canonical mesh build."""
import hashlib
import json
from pathlib import Path


CANONICAL_SOURCE_ROLES = ("Signature", "Boundary", "Mask", "Process",
                          "SemanticContract", "MeshRecipe")


def canonical_sha256(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _digests(value, description):
    if (not isinstance(value, dict) or not value or
            any(not isinstance(name, str) or not name or
                not isinstance(digest, str) or len(digest) != 64
                for name, digest in value.items())):
        raise ValueError(f"{description} must contain named SHA-256 digests")
    return dict(sorted(value.items()))


def canonical_cache_key(source_sha256, gates, canonical_tool_sha256):
    sources = _digests(source_sha256, "Canonical sources")
    if set(sources) != set(CANONICAL_SOURCE_ROLES):
        raise ValueError("Canonical sources do not match the immutable source contract")
    tools = _digests(canonical_tool_sha256, "Canonical tools")
    if not isinstance(gates, dict) or not gates:
        raise ValueError("Canonical qualification gates must be frozen")
    # Process, semantic and recipe hashes remain named explicitly even though all
    # immutable source roles are included. This makes cache-review failures legible.
    return {
        "Version": 1,
        "SourceSHA256": sources,
        "ProcessSHA256": sources["Process"],
        "SemanticContractSHA256": sources["SemanticContract"],
        "RecipeSHA256": sources["MeshRecipe"],
        "QualificationGates": gates,
        "CanonicalToolSHA256": tools,
    }


def build_record(source_sha256, gates, canonical_tool_sha256, artifacts):
    key = canonical_cache_key(source_sha256, gates, canonical_tool_sha256)
    artifact_bindings = {}
    for name, item in artifacts.items():
        if (not isinstance(item, dict) or not item.get("Path") or
                not isinstance(item.get("SHA256"), str) or len(item["SHA256"]) != 64):
            raise ValueError("Canonical artifacts require paths and SHA-256 digests")
        artifact_bindings[name] = {"Path": str(Path(item["Path"]).resolve()),
                                   "SHA256": item["SHA256"]}
    if not artifact_bindings or "candidate-mesh" not in artifact_bindings:
        raise ValueError("Canonical build must bind its published candidate")
    record = {
        "Version": 1,
        "CanonicalBuildId": canonical_sha256(key),
        "CanonicalCacheKey": key,
        "CanonicalArtifacts": dict(sorted(artifact_bindings.items())),
    }
    record["CanonicalBuildSHA256"] = canonical_sha256(record)
    return record


def validate_build_record(record, source_sha256, gates, canonical_tool_sha256,
                          *, check_files=True):
    expected_key = canonical_cache_key(source_sha256, gates, canonical_tool_sha256)
    if (record.get("Version") != 1 or
            record.get("CanonicalCacheKey") != expected_key or
            record.get("CanonicalBuildId") != canonical_sha256(expected_key)):
        raise ValueError("Canonical build cache key differs from source, gates, or tools")
    payload = dict(record)
    claimed = payload.pop("CanonicalBuildSHA256", None)
    if claimed != canonical_sha256(payload):
        raise ValueError("Canonical build hash differs from its immutable payload")
    artifacts = record.get("CanonicalArtifacts")
    if not isinstance(artifacts, dict) or "candidate-mesh" not in artifacts:
        raise ValueError("Canonical build has no candidate artifact")
    if check_files:
        from mesh_stage_contract import sha256
        for name, item in artifacts.items():
            if (not isinstance(item, dict) or not Path(item.get("Path", "")).is_file() or
                    sha256(item["Path"]) != item.get("SHA256")):
                raise ValueError(f"Canonical artifact changed: {name}")
    return record


def same_canonical_build(left, right):
    """True only for an exact cache-key and canonical-artifact-hash match."""
    return (left.get("CanonicalBuildId") == right.get("CanonicalBuildId") and
            left.get("CanonicalBuildSHA256") == right.get("CanonicalBuildSHA256") and
            left.get("CanonicalCacheKey") == right.get("CanonicalCacheKey") and
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
    parser.add_argument("--artifact", action="append", required=True, metavar="NAME=PATH")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("canonical build record must be fresh")
    artifacts = {}
    for value in args.artifact:
        if "=" not in value:
            parser.error("--artifact must be NAME=PATH")
        name, raw_path = value.split("=", 1); path = Path(raw_path).resolve()
        if not name or name in artifacts or not path.is_file():
            parser.error("canonical artifact names must be unique and paths must exist")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        artifacts[name] = {"Path": str(path), "SHA256": digest}
    try:
        record = build_record(json.loads(args.source_hashes.read_text()),
                              json.loads(args.gates.read_text()),
                              json.loads(args.canonical_tool_hashes.read_text()), artifacts)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    args.output.write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
