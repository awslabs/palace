#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Hash-bound, read-only launcher for experimental archived electrostatic diagnostics.

The request contains the entire source-ID/SHA256 table, fixed ZeroTraceIndices,
and dense coefficient rows. No nonzero coefficient is dropped. Does not solve the PDE.
Use the normal Palace EstimatorTol/EstimatorMaxIts/EstimatorMG controls in the config.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys

HEADER = struct.Struct("=QIIqqqq")
MAGIC = 0x50414C5253503031
CONFLICTS = ("PALACE_RESPONSE_ARCHIVE_ONLY", "PALACE_RESPONSE_REDUCE_ONLY",
             "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS", "PALACE_RESPONSE_BLOCK_SIZE")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    def reject(value):
        raise ValueError(f"nonfinite JSON number: {value}")
    return json.loads(Path(path).read_text(), parse_constant=reject)


def write(path, data):
    with Path(path).open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def validate_request(request, sources):
    ids = [s["Index"] for s in sources]
    if ids != sorted(set(ids)) or any(type(i) is not int or i <= 0 for i in ids):
        raise ValueError("configuration source IDs must be unique, positive and sorted")
    if request["Version"] != 1 or request["SourceIds"] != ids:
        raise ValueError("request SourceIds must equal the complete configuration source table")
    if request["SourceSHA256"] != [digest(s["DataFile"]) for s in sources]:
        raise ValueError("source SHA256 table mismatch")
    zero = request["ZeroTraceIndices"]
    if len(set(zero)) != len(zero) or not set(zero) <= set(ids):
        raise ValueError("invalid fixed ZeroTraceIndices")
    names = set()
    used = set()
    for row in request["Excitations"]:
        name, coefficients = row["Name"], row["Coefficients"]
        if not isinstance(name, str) or not name or name in names or len(coefficients) != len(ids):
            raise ValueError("excitation names must be unique and coefficient rows complete")
        names.add(name)
        for index, coefficient in zip(ids, coefficients):
            if type(coefficient) not in (int, float) or not math.isfinite(coefficient):
                raise ValueError("nonfinite/invalid coefficient")
            if coefficient != 0:
                if index in zero:
                    raise ValueError("nonzero coefficient on constrained source")
                used.add(index)
    if not names:
        raise ValueError("no excitations")
    for name in ("RelResidualTol", "AbsResidualTol", "BCAbsTolV", "MinEnergyJ", "MinCancellationRatio"):
        value = request["Validation"][name]
        if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
            raise ValueError(f"invalid explicit validation tolerance: {name}")
    if request["Validation"]["MinEnergyJ"] <= 0 or not 0 < request["Validation"]["MinCancellationRatio"] < 1:
        raise ValueError("MinEnergyJ must be positive and MinCancellationRatio must be in (0,1)")
    if type(request["WriteElementIndicators"]) is not bool:
        raise ValueError("WriteElementIndicators must be boolean")
    if "SurfaceQuadratureExtras" in request:
        extras = request["SurfaceQuadratureExtras"]
        if (not isinstance(extras, list) or not extras or extras[0] != 0 or
                any(type(x) is not int or x < 0 or x > 12 for x in extras) or
                extras != sorted(set(extras))):
            raise ValueError("SurfaceQuadratureExtras must be sorted unique integers in [0,12], starting at0")
    return used


def validate_basis_contract(request, sources):
    binding = request["BasisContract"]
    path = Path(binding["Path"]).resolve()
    if digest(path) != binding["SHA256"]:
        raise ValueError("basis contract SHA256 mismatch")
    contract = load(path)
    if contract["Sources"] != len(sources) or contract["ZeroTraceIndices"] != request["ZeroTraceIndices"]:
        raise ValueError("basis contract source count or fixed ZeroTraceIndices mismatch")
    # The preserved coupon contract encodes original source indices in these basenames.
    # Bind every association, including sources unused by the requested excitations.
    expected = []
    for filename, sha in contract["OutputSourceSHA256"].items():
        name = Path(filename).name
        match = re.fullmatch(r"basis-([0-9]{4})\.csv", name)
        if not match or int(match[1]) <= 0:
            raise ValueError("basis contract requires indexed basis-NNNN.csv filenames")
        expected.append((int(match[1]), name, sha))
    expected.sort()
    if len({entry[0] for entry in expected}) != len(expected):
        raise ValueError("ambiguous basis contract source indices")
    actual = [(s["Index"], Path(s["DataFile"]).name, digest(s["DataFile"])) for s in sources]
    if actual != expected:
        raise ValueError("basis contract ordered source index/name/hash mapping mismatch")
    return path


def validate_supported_inputs(config):
    model = config["Model"]
    if model.get("ExportPrerefinedMesh", False):
        raise ValueError("archive estimation rejects ExportPrerefinedMesh")
    if model.get("Partitioning", ""):
        raise ValueError("archive estimation does not support Partitioning")
    for surface in config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", []):
        if surface.get("OwnershipDataFile", ""):
            raise ValueError("archive estimation does not support OwnershipDataFile")


def archive_inputs(directory, used, ranks):
    paths, widths = [], {}
    for index in sorted(used):
        for rank in range(ranks):
            path = directory / f"source-{index:06d}-rank-{rank:06d}-V.bin"
            with path.open("rb") as stream:
                data = stream.read(HEADER.size)
            if len(data) != HEADER.size:
                raise ValueError(f"truncated archive header: {path}")
            magic, version, field, source, owner, size, width = HEADER.unpack(data)
            if (magic, version, field, source, owner, size) != (MAGIC, 1, 1, index, rank, ranks):
                raise ValueError(f"archive header mismatch: {path}")
            if width < 0 or width > 2**31 - 1 or path.stat().st_size != HEADER.size + 8 * width:
                raise ValueError(f"archive payload size mismatch: {path}")
            if rank in widths and widths[rank] != width:
                raise ValueError(f"inconsistent per-rank DOF count: {path}")
            widths[rank] = width
            paths.append(path)
    return paths, widths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True, help="new run directory, must not exist")
    parser.add_argument("--ranks", type=int, required=True)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--seconds", type=float, required=True)
    parser.add_argument("--memory-gib", type=float, required=True)
    args = parser.parse_args()
    for name in CONFLICTS:
        if name in os.environ:
            raise ValueError(f"conflicting archive flag: {name}")
    if args.ranks <= 0 or not all(math.isfinite(v) and v > 0 for v in (args.seconds, args.memory_gib)):
        raise ValueError("positive finite resource limits required")
    config_path, request_path = args.config.resolve(), args.request.resolve()
    binary, archive, output = args.binary.resolve(), args.archive.resolve(), args.output.resolve()
    if output.exists() or not archive.is_dir():
        raise ValueError("output must be new and archive directory must exist")
    if digest(binary) != args.binary_sha256:
        raise ValueError("binary SHA256 mismatch")
    config, request = load(config_path), load(request_path)
    validate_supported_inputs(config)
    sources = config["Boundaries"]["PrescribedPotential"]
    # Like the command-line solver, relative input paths refer to the invocation cwd.
    config["Model"]["Mesh"] = str(Path(config["Model"]["Mesh"]).resolve())
    for source in sources:
        source["DataFile"] = str(Path(source["DataFile"]).resolve())
    used = validate_request(request, sources)
    basis_contract = validate_basis_contract(request, sources)
    archive_paths, widths = archive_inputs(archive, used, args.ranks)
    inputs = [config_path, request_path, basis_contract, binary, Path(config["Model"]["Mesh"]),
              *[Path(s["DataFile"]) for s in sources], *archive_paths]
    for path in [archive, *inputs]:
        path = path.resolve()
        if path == output or path in output.parents or output in path.parents:
            raise ValueError(f"output overlaps input: {path}")
    hashes = {str(path): digest(path) for path in inputs}
    output.mkdir(parents=True, exist_ok=False)
    config["Problem"]["Output"] = str(output / "postpro")
    write(output / "config.json", config)
    write(output / "request.json", request)
    manifest = {"Status": "incomplete", "InputSHA256": hashes,
                "RunConfigSHA256": digest(output / "config.json"),
                "RequestSHA256": digest(output / "request.json"),
                "LauncherSHA256": digest(__file__), "Ranks": args.ranks,
                "ArchiveRankWidths": widths, "UsedSources": sorted(used),
                "SourceCount": len(sources), "ArchiveInputMutation": "not-yet-checked"}
    write(output / "provenance-before.json", manifest)
    env = dict(os.environ, PALACE_RESPONSE_ESTIMATE_ONLY="1",
               PALACE_RESPONSE_ESTIMATE_REQUEST=str(output / "request.json"),
               PALACE_RESPONSE_ARCHIVE_DIR=str(archive), PYTHONDONTWRITEBYTECODE="1")
    command = [args.mpiexec, "-n", str(args.ranks), str(binary), str(output / "config.json")]
    bounded = Path(__file__).with_name("run_bounded_mesher.py")
    try:
        result = subprocess.run([sys.executable, str(bounded), "--seconds", str(args.seconds),
                                 "--memory-gib", str(args.memory_gib), "--log", str(output / "solver.log"),
                                 "--", *command], env=env)
        manifest["ReturnCode"] = result.returncode
        manifest["Status"] = "complete" if result.returncode == 0 else "failed"
    finally:
        changed = [path for path, old in hashes.items() if not Path(path).exists() or digest(path) != old]
        manifest["ChangedInputs"] = changed
        manifest["OutputSHA256"] = {str(p.relative_to(output)): digest(p)
                                     for p in sorted(output.rglob('*')) if p.is_file()}
        manifest["ArchiveInputMutation"] = "detected" if changed else "none"
        write(output / "provenance-after.json", manifest)
    if changed:
        raise RuntimeError(f"input mutation detected: {changed}")
    raise SystemExit(result.returncode)


if __name__ == "__main__":
    main()
