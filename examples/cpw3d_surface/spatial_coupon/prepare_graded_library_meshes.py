#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Freeze retained spatial-coupon contracts and emit bounded graded-mesh commands."""

import argparse
import copy
import hashlib
import json
import os
import pathlib
import shlex
import shutil


def digest(path):
    result = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(value, indent=2) + "\n")
    temporary.replace(path)


def surface_attributes(config):
    attributes = set(config["Boundaries"]["Ground"]["Attributes"])
    for source in config["Boundaries"]["PrescribedPotential"]:
        attributes.update(source.get("TerminalAttributes", []))
    for item in config["Boundaries"].get("Postprocessing", {}).get("Dielectric", []):
        attributes.update(item["Attributes"])
        attributes.update(item.get("EdgeAttributes", []))
    attributes.difference_update(config["Boundaries"]["PrescribedPotential"][0]["Attributes"])
    return sorted(attributes)


def source_contract(config):
    sources = config["Boundaries"]["PrescribedPotential"]
    files = []
    for source in sources:
        path = pathlib.Path(source["DataFile"])
        files.append(
            {
                "Index": source["Index"],
                "Attributes": source["Attributes"],
                "TerminalAttributes": source.get("TerminalAttributes", []),
                "Path": str(path),
                "Bytes": path.stat().st_size,
                "SHA256": digest(path),
            }
        )
    return {
        "Order": config["Solver"]["Order"],
        "Linear": copy.deepcopy(config["Solver"]["Linear"]),
        "SourceCount": len(sources),
        "Sources": files,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference_manifest", type=pathlib.Path)
    parser.add_argument("output", type=pathlib.Path)
    parser.add_argument("--tools", type=pathlib.Path, required=True)
    parser.add_argument("--julia-project", type=pathlib.Path, required=True)
    parser.add_argument("--profile", choices=("selected", "thin", "fabricated"), default="selected")
    parser.add_argument(
        "--matching-trace", choices=("none", "levels", "sides", "all"), required=True,
        help="Explicit experimental policy; none is diagnostic-only. Every policy needs trace/accuracy gates.",
    )
    arguments = parser.parse_args()
    output = arguments.output.resolve()
    manifest_path = output / "mesh-plan.json"
    manifest_path.exists() and parser.error(f"Refuse to overwrite {manifest_path}")
    reference = json.load(open(arguments.reference_manifest))
    tools = arguments.tools.resolve()
    inputs_root = output / "meshes" / "inputs"
    candidates = output / "meshes" / "candidates"
    inputs_root.mkdir(parents=True, exist_ok=True)
    candidates.mkdir(parents=True, exist_ok=True)
    process = {
        "Units": "um",
        "Radius": 2.0,
        "MetalThickness": 0.1,
        "Overetch": 0.05,
        "SidewallAngle": 90.0,
        "TopRounding": 0.0,
        "TrenchRounding": 0.0,
    }
    records = []
    cases = [case for case in reference["Cases"] if case["Model"].startswith("spatialedgecluster_")]
    for case in cases:
        key = case["Key"]
        kind = case["Kind"]
        if arguments.profile != "selected" and kind != arguments.profile:
            continue
        config_path = pathlib.Path(case["SourceConfig"])
        config = json.load(open(config_path))
        library_contract = {}
        if case.get("SourceLibrary"):
            library_path = pathlib.Path(case["SourceLibrary"])
            model = next(item for item in json.load(open(library_path))["Models"]
                         if item["Name"] == case["Model"])
            library_contract = {
                "Path": str(library_path), "SHA256": digest(library_path),
                "Model": copy.deepcopy(model),
                "ZeroTraceIndices": model.get("ZeroTraceIndices", []),
            }
        model_key = key.split("-")[0] + "-" + case["Model"]
        root = inputs_root / model_key
        root.mkdir(parents=True, exist_ok=True)
        origin = config_path.parent
        geometry_hashes = {}
        for name in ("mesh-signature.csv", "plan-view-mask.csv", "plan-view-boundary.csv"):
            source = origin / name
            if not source.is_file():
                raise RuntimeError(f"Missing retained geometry input {source}")
            destination = root / name
            if not destination.exists():
                shutil.copy2(source, destination)
            if digest(destination) != digest(source):
                raise RuntimeError(f"Staged geometry differs: {name}")
            geometry_hashes[name] = digest(source)
        process_path = root / "process.toml"
        process_text = "\n".join(
            f'{name} = "{value}"' if isinstance(value, str) else f"{name} = {str(value).lower()}"
            for name, value in process.items()
        ) + "\n"
        if process_path.exists() and process_path.read_text() != process_text:
            raise RuntimeError(f"Process contract differs at {process_path}")
        process_path.write_text(process_text)
        first_trace = pathlib.Path(config["Boundaries"]["PrescribedPotential"][0]["DataFile"])
        trace_path = root / "matching-trace.csv"
        if not trace_path.exists():
            shutil.copy2(first_trace, trace_path)
        if digest(trace_path) != digest(first_trace):
            raise RuntimeError("Trace copy differs")
        expected = root / f"expected-{kind}.csv"
        lines = ["dimension,attribute,measure", "2,1,0"]
        lines.extend(f"2,{attribute},0" for attribute in surface_attributes(config))
        lines.extend(("3,1,0", "3,2,0"))
        expected.write_text("\n".join(lines) + "\n")
        fine = 0.002 if kind == "thin" else 0.0005
        study = root / f"volume-{kind}.toml"
        study.write_text(
            "MaxElements = 4000000\nMaxNodes = 1200000\nOptimizeVolume = true\n\n"
            "[[Variants]]\nName = \"selected\"\nMinimumSize = 0.002\n"
            "NearGrowth = 0.5\nFarGrowth = 2.0\nTransitionDistance = 0.03\nMaximumSize = 0.5\n"
        )
        prefix = candidates / key
        geometry = pathlib.Path(str(prefix) + "-geometry.msh")
        family = pathlib.Path(str(prefix) + "-family.msh")
        frozen = pathlib.Path(str(prefix) + "-family-selected.msh")
        final = pathlib.Path(str(prefix) + ".msh")
        common = [
            "julia",
            "--startup-file=no",
            f"--project={arguments.julia_project.resolve()}",
            str(tools / "mesh_graded_tet_experiment.jl"),
            str(root),
            kind,
        ]
        # Even the no-CAD-line diagnostic validates the immutable trace geometry.
        trace_arguments = ["--matching-trace", str(trace_path)]
        geometry_command = [
            *common,
            str(geometry),
            "1.0",
            str(fine),
            "0.5",
            "--process",
            str(process_path),
            "--geometry-only",
            *trace_arguments,
        ]
        volume_command = [
            *common,
            str(family),
            "1.0",
            str(fine),
            "0.5",
            "--process",
            str(process_path),
            "--volume-study",
            str(study),
            "--reference-measures",
            str(geometry) + ".cad-measures.csv",
            *trace_arguments,
        ]
        relabel_command = [
            "julia",
            "--startup-file=no",
            f"--project={arguments.julia_project.resolve()}",
            str(tools / "relabel_frozen_interface_mesh.jl"),
            str(root),
            kind,
            str(frozen),
            str(final),
            "--expected-measures",
            str(expected),
        ]
        records.append(
            {
                "Key": key,
                "Model": case["Model"],
                "Kind": kind,
                "RetainedConfig": str(config_path),
                "RetainedConfigSHA256": digest(config_path),
                "RetainedMesh": case["Mesh"],
                "RetainedMeshSHA256": reference["Inputs"][case["Mesh"]]["SHA256"],
                "GeometryInputs": geometry_hashes,
                "Process": process,
                "MatchingTraceMode": arguments.matching_trace,
                "MatchingTraceSHA256": digest(trace_path),
                "FineSurfaceSizeMicrometres": fine,
                "MinimumVolumeSizeMicrometres": 0.002,
                "ExpectedInterfaceAttributes": surface_attributes(config),
                "SourceContract": source_contract(config),
                "LibraryContract": library_contract,
                "GeometryCommand": geometry_command,
                "VolumeCommand": volume_command,
                "RelabelCommand": relabel_command,
                "GeometryOutput": str(geometry),
                "FrozenOutput": str(frozen),
                "MeshOutput": str(final),
            }
        )
    plan = {
        "Version": 1,
        "Status": "PreparedNotGenerated",
        "ReferenceCampaign": str(arguments.reference_manifest.resolve()),
        "ReferenceCampaignSHA256": digest(arguments.reference_manifest),
        "Tools": str(tools),
        "JuliaProject": str(arguments.julia_project.resolve()),
        "MatchingTraceMode": arguments.matching_trace,
        "TraceSizing": {"TraceSize": 0.0, "TraceSizeScope": "off",
                        "TraceRelativeSize": 0.0, "TraceSurfaceGrowth": 4.0},
        "MeshFormat": "Gmsh 2.2 binary",
        "Scope": "Six retained spatial models; original sources, orders, and solver settings are immutable.",
        "Cases": records,
    }
    write_json(manifest_path, plan)
    runner = output / "generate-meshes.sh"
    environment = {
        "JULIA_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "TET_GEOMETRY_ORDER": "1",
        "TET_SURFACE_ALGORITHM": "5",
        "TET_ALGORITHM3D": "10",
        "TET_HXT_QUALITY": "0.1",
        "TET_TRACE_CONSTRAINT_MODE": arguments.matching_trace,
        "TET_TRACE_SIZE": "0",
        "TET_TRACE_SIZE_SCOPE": "off",
        "TET_TRACE_RELATIVE_SIZE": "0",
        "TET_TRACE_SURFACE_GROWTH": "4",
        "TET_VERBOSITY": "3",
    }
    lines = ["#!/bin/bash", "set -euo pipefail"]
    lines.extend(f"export {name}={shlex.quote(value)}" for name, value in environment.items())
    lines.append('only="${1:-all}"')
    for record in records:
        lines.append(f'if [[ "$only" == all || "$only" == {shlex.quote(record["Key"])} ]]; then')
        for stage, command, seconds in (
            ("geometry", record["GeometryCommand"], 900),
            ("volume", record["VolumeCommand"], 7200),
            ("relabel", record["RelabelCommand"], 1800),
        ):
            log = candidates / f'{record["Key"]}-{stage}.log'
            bounded = [
                "python3",
                str(tools / "run_bounded_mesher.py"),
                "--seconds",
                str(seconds),
                "--memory-gib",
                "110",
                "--log",
                str(log),
                "--",
                *command,
            ]
            lines.append("  " + shlex.join(bounded))
        lines.append("fi")
    runner.write_text("\n".join(lines) + "\n")
    runner.chmod(0o755)
    print(json.dumps({"Cases": len(records), "Manifest": str(manifest_path), "Runner": str(runner)}, indent=2))


if __name__ == "__main__":
    main()
