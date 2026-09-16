#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Refreeze the repository-tool SHA-256 digests of the production suite manifest and
mirror its `Tools` / `StageToolSHA256` into the labeled calibration manifest.

Only tools that live in the repository are recomputed: every `Tools` entry and the
stage-tool roles listed in STAGE_REPOSITORY_TOOLS.  Runtimes, the reviewed adapter
and the MMG library are machine-bound identities and are never touched here.  The
calibration manifest must carry exactly the production tool digests (asserted by
`test_general_mesh_manifest.py`), so both manifests are always refrozen together.
`--check` reports stale digests or a stale mirror without writing (exit 1).
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PRODUCTION_MANIFEST = HERE / "geometry-independence-suite.json"
CALIBRATION_MANIFEST = HERE / "geometry-independence-calibration-ma.json"
# Stage tool roles whose identity is a repository file (relative to this directory).
STAGE_REPOSITORY_TOOLS = {
    ("canonical-source-validation", "source-validator"): "transform_coupon_source_contract.py",
    ("seed-generation", "mesher"): "mesh_spatial_coupon.jl",
    ("metric-preparation", "metric-preparer"): "prepare_edge_metric_scout.py",
    ("native-adaptation-mmg", "adaptation-wrapper"): "run_native_mmg_adaptation.py",
    ("label-restoration", "label-restorer"): "restore_planar_metric_mesh.py",
    ("canonical-gmsh-publication", "publisher"): "relabel_frozen_interface_mesh.jl",
    ("proper-rigid-publication", "rigid-publisher"): "publish_rigid_coupon_mesh.py",
    ("proper-rigid-publication", "ownership-auditor"): "audit_rigid_coupon_ownership.jl",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def refrozen_production(manifest, manifest_path):
    """The production manifest with current repository-tool digests and the list of
    `(name, old, new)` changes."""
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    changes = []
    for tool in manifest["Tools"]:
        current = sha256(repository / tool["Path"])
        if current != tool["SHA256"]:
            changes.append((tool["Name"], tool["SHA256"], current))
            tool["SHA256"] = current
    for (stage, role), name in STAGE_REPOSITORY_TOOLS.items():
        current = sha256(manifest_path.parent / name)
        if manifest["StageToolSHA256"][stage][role] != current:
            changes.append((f"{stage}/{role}", manifest["StageToolSHA256"][stage][role], current))
            manifest["StageToolSHA256"][stage][role] = current
    return manifest, changes


def refreeze(production_path, calibration_path, *, check_only):
    """Refreeze both manifests (or only report); returns (changes, mirror_was_stale)."""
    production = json.loads(production_path.read_text())
    calibration = json.loads(calibration_path.read_text())
    production, changes = refrozen_production(production, production_path)
    mirror_stale = (calibration["Tools"] != production["Tools"] or
                    calibration["StageToolSHA256"] != production["StageToolSHA256"])
    if not check_only and (changes or mirror_stale):
        calibration["Tools"] = production["Tools"]
        calibration["StageToolSHA256"] = production["StageToolSHA256"]
        production_path.write_text(json.dumps(production, indent=2) + "\n")
        calibration_path.write_text(json.dumps(calibration, indent=2) + "\n")
    return changes, mirror_stale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production", type=Path, default=PRODUCTION_MANIFEST)
    parser.add_argument("--calibration", type=Path, default=CALIBRATION_MANIFEST)
    parser.add_argument("--check", action="store_true",
                        help="report stale digests without writing; exit 1 when stale")
    args = parser.parse_args()
    changes, mirror_stale = refreeze(args.production, args.calibration, check_only=args.check)
    for name, old, new in changes:
        print(f"{'stale' if args.check else 'refroze'} {name}: {old[:12]} -> {new[:12]}")
    if mirror_stale:
        print("calibration mirror " + ("stale" if args.check else "updated"))
    if not changes and not mirror_stale:
        print("frozen hashes current; calibration mirror equal")
    return 1 if args.check and (changes or mirror_stale) else 0


if __name__ == "__main__":
    sys.exit(main())
