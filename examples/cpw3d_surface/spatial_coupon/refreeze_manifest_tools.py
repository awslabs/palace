#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Refreeze the repository-tool SHA-256 digests of the production suite manifest and
mirror its `Tools` / shared `StageToolSHA256` stages into the labeled calibration
manifests.

The production manifest freezes the Gmsh-only pipeline's stages (supervisor decision
38: canonical-source-validation, gmsh-build, canonical-gmsh-publication,
proper-rigid-publication); the MA calibration manifest keeps the legacy MMG stages
(seed-generation, metric-preparation, native-adaptation-mmg, label-restoration) and
shares the other three; the sizing calibration manifest (supervisor decision 41)
freezes the Gmsh-only stages and shares all four.  Each manifest is refrozen over the stages it freezes; the
shared stages and `Tools` are mirrored from production.  Only tools that live in the
repository are recomputed: every `Tools` entry and the stage-tool roles listed in
STAGE_REPOSITORY_TOOLS.  Runtimes and the MMG library are machine-bound identities
and are never touched here.  The reviewed native adapter is
machine-bound too (a build of `adapt_edge_metric.cpp`, never committed): its
`native-adaptation-mmg/adapter-mmg` digest changes only when an explicit
`--adapter-mmg PATH` names the executable of a recorded build
(`testdata/adapter-build.json`: command, compiler, source/exe/dylib SHA-256, rpath)
whose source digest is the repository's `adapt_edge_metric.cpp`.  The Julia runtime
roles (JULIA_RUNTIME_ROLES) are machine-bound likewise and change only through an
explicit `--julia-runtime PATH` naming the launcher executable actually run.  Every calibration
manifest must carry exactly the production tool digests (asserted by
`test_general_mesh_manifest.py`), so the manifests are always refrozen together.
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
# The labeled calibration manifests: the legacy MMG edge-layer study and the
# Gmsh-only sizing calibration (supervisor decision 41; shares every stage).
CALIBRATION_MANIFESTS = (CALIBRATION_MANIFEST, HERE / "geometry-independence-calibration-sizing.json")
# Stage tool roles whose identity is a repository file (relative to this directory).
STAGE_REPOSITORY_TOOLS = {
    ("canonical-source-validation", "source-validator"): "transform_coupon_source_contract.py",
    ("seed-generation", "mesher"): "mesh_spatial_coupon.jl",
    ("gmsh-build", "mesher"): "mesh_spatial_coupon.jl",
    ("metric-preparation", "metric-preparer"): "prepare_edge_metric_scout.py",
    ("native-adaptation-mmg", "adaptation-wrapper"): "run_native_mmg_adaptation.py",
    ("label-restoration", "label-restorer"): "restore_planar_metric_mesh.py",
    ("canonical-gmsh-publication", "publisher"): "relabel_frozen_interface_mesh.jl",
    ("proper-rigid-publication", "rigid-publisher"): "publish_rigid_coupon_mesh.py",
    ("proper-rigid-publication", "ownership-auditor"): "audit_rigid_coupon_ownership.jl",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# Stage tool roles whose identity is the Julia launcher executable of this machine.
JULIA_RUNTIME_ROLES = (("seed-generation", "runtime"), ("gmsh-build", "runtime"),
                       ("canonical-gmsh-publication", "runtime"),
                       ("proper-rigid-publication", "ownership-runtime"))
ADAPTER_SOURCE = "adapt_edge_metric.cpp"
ADAPTER_BUILD_RECORD = "testdata/adapter-build.json"


def adapter_digest(adapter, manifest_path):
    """SHA-256 of the reviewed adapter executable, accepted only when the recorded
    build (`testdata/adapter-build.json`) names this executable digest and the
    repository's current `adapt_edge_metric.cpp` digest."""
    adapter = Path(adapter).resolve()
    if not adapter.is_file():
        raise ValueError(f"Adapter executable does not exist: {adapter}")
    digest = sha256(adapter)
    record = json.loads((manifest_path.parent / ADAPTER_BUILD_RECORD).read_text())
    if (record.get("ExecutableSHA256") != digest or
            record.get("SourceSHA256") != sha256(manifest_path.parent / ADAPTER_SOURCE)):
        raise ValueError("Adapter executable or adapt_edge_metric.cpp differ from the recorded build")
    return digest


def refrozen_production(manifest, manifest_path, adapter=None, julia_runtime=None):
    """The manifest with current repository-tool digests over the stages it freezes
    (the recorded adapter digest when `adapter` is given and the manifest freezes the
    MMG stage, the Julia launcher digest when `julia_runtime` is given) and the list
    of `(name, old, new)` changes."""
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    stage_tools = manifest["StageToolSHA256"]
    changes = []
    if julia_runtime is not None:
        launcher = Path(julia_runtime).resolve()
        if not launcher.is_file():
            raise ValueError(f"Julia runtime executable does not exist: {launcher}")
        current = sha256(launcher)
        for stage, role in JULIA_RUNTIME_ROLES:
            if stage in stage_tools and stage_tools[stage][role] != current:
                changes.append((f"{stage}/{role}", stage_tools[stage][role], current))
                stage_tools[stage][role] = current
    if adapter is not None and "native-adaptation-mmg" in stage_tools:
        current = adapter_digest(adapter, manifest_path)
        role = stage_tools["native-adaptation-mmg"]
        if role["adapter-mmg"] != current:
            changes.append(("native-adaptation-mmg/adapter-mmg", role["adapter-mmg"], current))
            role["adapter-mmg"] = current
    for tool in manifest["Tools"]:
        current = sha256(repository / tool["Path"])
        if current != tool["SHA256"]:
            changes.append((tool["Name"], tool["SHA256"], current))
            tool["SHA256"] = current
    for (stage, role), name in STAGE_REPOSITORY_TOOLS.items():
        if stage not in stage_tools:
            continue
        current = sha256(manifest_path.parent / name)
        if stage_tools[stage][role] != current:
            changes.append((f"{stage}/{role}", stage_tools[stage][role], current))
            stage_tools[stage][role] = current
    return manifest, changes


def shared_stages(production, calibration):
    return sorted(set(production["StageToolSHA256"]) & set(calibration["StageToolSHA256"]))


def refreeze(production_path, *calibration_paths, check_only, adapter=None, julia_runtime=None):
    """Refreeze the production manifest and every calibration manifest (or only
    report); returns (changes, mirror_was_stale).  A calibration manifest's own stages
    (the legacy stages) are refrozen from the repository; its `Tools` and the stages
    shared with production are mirrored from production."""
    production = json.loads(production_path.read_text())
    calibrations = [json.loads(path.read_text()) for path in calibration_paths]
    production, changes = refrozen_production(production, production_path, adapter, julia_runtime)
    mirror_stale = False
    for calibration_path, calibration in zip(calibration_paths, calibrations):
        calibration, calibration_changes = refrozen_production(calibration, calibration_path, adapter,
                                                               julia_runtime)
        # The calibration manifest's own stage roles (its legacy stages); its `Tools`
        # and shared stages are mirrored below.
        shared = shared_stages(production, calibration)
        changes += [(f"{calibration_path.name} {name}", old, new) for name, old, new in calibration_changes
                    if "/" in name and name.split("/")[0] not in shared]
        mirror_stale |= (calibration["Tools"] != production["Tools"] or
                         any(calibration["StageToolSHA256"][stage] != production["StageToolSHA256"][stage]
                             for stage in shared))
    if not check_only and (changes or mirror_stale):
        production_path.write_text(json.dumps(production, indent=2) + "\n")
        for calibration_path, calibration in zip(calibration_paths, calibrations):
            calibration["Tools"] = production["Tools"]
            for stage in shared_stages(production, calibration):
                calibration["StageToolSHA256"][stage] = production["StageToolSHA256"][stage]
            calibration_path.write_text(json.dumps(calibration, indent=2) + "\n")
    return changes, mirror_stale


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production", type=Path, default=PRODUCTION_MANIFEST)
    parser.add_argument("--calibration", type=Path, nargs="+", default=list(CALIBRATION_MANIFESTS),
                        help="labeled calibration manifests mirrored from production")
    parser.add_argument("--check", action="store_true",
                        help="report stale digests without writing; exit 1 when stale")
    parser.add_argument("--adapter-mmg", type=Path,
                        help="executable of the recorded adapter build whose digest becomes "
                             "the native-adaptation-mmg/adapter-mmg stage tool")
    parser.add_argument("--julia-runtime", type=Path,
                        help="Julia launcher executable of this machine whose digest becomes "
                             "the Julia runtime stage tools (seed, publisher, ownership)")
    args = parser.parse_args()
    try:
        changes, mirror_stale = refreeze(args.production, *args.calibration,
                                         check_only=args.check, adapter=args.adapter_mmg,
                                         julia_runtime=args.julia_runtime)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        parser.error(str(error))
    for name, old, new in changes:
        print(f"{'stale' if args.check else 'refroze'} {name}: {old[:12]} -> {new[:12]}")
    if mirror_stale:
        print("calibration mirror " + ("stale" if args.check else "updated"))
    if not changes and not mirror_stale:
        print("frozen hashes current; calibration mirror equal")
    return 1 if args.check and (changes or mirror_stale) else 0


if __name__ == "__main__":
    sys.exit(main())
