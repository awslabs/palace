#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Configure + build the palace executable from a frozen CMake source tree on a Linux node.

Replaces the frozen-compile-list replay (build_linux_*.py: 102 recorded compile commands of one
historical configure) with a real configure of the frozen tree against the existing read-only
Linux dependency prefix, so any source change (new files, new CMake options, new schema files)
builds. The configure arguments are the ones recorded in the dependency prefix's own
`palace-cmake/tmp/palace-cfgcmd.txt` (ExternalProject step of the superbuild: compilers, ARM
Performance Libraries BLAS / LAPACK, every PALACE_WITH_* option, Release, static), with the
source and build directories re-pointed and `PALACE_BUILD_EXTERNAL_DEPS=OFF` (no Catch2
download; the unit tests are not built). Provenance: SHA-256 of every frozen source (checked
against `manifest.json` before and after), of the dependency prefix headers and libraries, the
compilers, cmake and mpirun; the configure command; the compiler / MPI versions; the resulting
binary's SHA-256 and `ldd`.

    python3 build_from_cmake_tree.py --source SRC --build BUILD --deps PREFIX [--jobs 96]

Runs inside a PBS job (job-linux-build.pbs): the login node has too few cores for the build.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def save(directory, name, data):
    with (directory / name).open("x") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def recorded_configure_arguments(deps):
    """The superbuild's recorded ExternalProject configure command for palace (a cmake list)."""
    text = (deps / "palace-cmake/tmp/palace-cfgcmd.txt").read_text()
    match = re.search(r"cmd='(.*)'", text, re.S)
    if not match:
        raise SystemExit(f"no configure command in {deps}/palace-cmake/tmp/palace-cfgcmd.txt")
    arguments = match.group(1).replace("$<SEMICOLON>", "\x00").split(";")
    return [a.replace("\x00", ";") for a in arguments]


def dependency_inputs(deps, tools):
    inputs = set(tools)
    inputs.update(p for p in (deps / "include").rglob("*") if p.is_file())
    inputs.update(p for p in (deps / "lib").glob("*") if p.is_file())
    inputs.update(p for p in (deps / "lib64").glob("*") if p.is_file() and (deps / "lib64").is_dir())
    inputs.add(deps / "palace-cmake/tmp/palace-cfgcmd.txt")
    return sorted(inputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="frozen tree (freeze_source_tree.py output) with manifest.json")
    parser.add_argument("--build", type=Path, required=True, help="new build root (must not exist)")
    parser.add_argument("--deps", type=Path, default=Path("/data/home/simlap/palace_builds/palace/build_c8g_main"),
                        help="read-only dependency prefix of the recorded superbuild")
    parser.add_argument("--cmake", type=Path, default=Path("/opt/cmake/bin/cmake"))
    parser.add_argument("--mpirun", type=Path, default=Path("/opt/openmpi/bin/mpirun"))
    parser.add_argument("--jobs", type=int, default=96)
    parser.add_argument("--binary-root", type=Path, help="where palace-<HEAD7>-<sha12>.bin is placed (default: parent of --build)")
    args = parser.parse_args()
    source = args.source.resolve()
    build = args.build.resolve()
    deps = args.deps.resolve()
    if build.exists():
        raise SystemExit(f"{build} exists: refusing to rebuild in place (every build is a new directory)")
    manifest = json.loads((source / "manifest.json").read_text())
    mismatched = [p for p, h in manifest["Files"].items() if not (source / p).is_file() or sha256(source / p) != h]
    if mismatched:
        raise SystemExit(f"frozen source does not match its manifest: {mismatched[:5]} ...")

    recorded = recorded_configure_arguments(deps)
    cmake = str(args.cmake)
    cxx = next(a.split("=", 1)[1] for a in recorded if a.startswith("-DCMAKE_CXX_COMPILER="))
    fortran = next(a.split("=", 1)[1] for a in recorded if a.startswith("-DCMAKE_Fortran_COMPILER="))
    describe = manifest.get("Describe", manifest["HEAD"][:12])
    git_flags = f'-DPALACE_GIT_COMMIT -DPALACE_GIT_COMMIT_ID=\\"{describe}-frozen\\"'
    configure = [cmake, str(source / "palace")]
    for argument in recorded[2:]:
        if argument.startswith("-DCMAKE_INSTALL_PREFIX="):
            argument = f"-DCMAKE_INSTALL_PREFIX={build}"
        elif argument.startswith("-DCMAKE_CXX_FLAGS="):
            # Only the git identity define (the frozen tree has no .git): codegen flags unchanged.
            argument = f"-DCMAKE_CXX_FLAGS={git_flags}"
        elif argument.startswith("-DPALACE_BUILD_EXTERNAL_DEPS="):
            argument = "-DPALACE_BUILD_EXTERNAL_DEPS=OFF"
        configure.append(argument)
    build.mkdir(parents=True)
    palace_build = build / "palace-build"
    palace_build.mkdir()

    tools = [Path(cxx), Path(fortran), args.cmake, args.mpirun, Path(__file__).resolve()]
    inputs = dependency_inputs(deps, tools)
    before = {str(p): sha256(p) for p in inputs}
    save(build, "provenance-before.json", {
        "SourceManifest": manifest, "SourceRoot": str(source), "DependencyPrefix": str(deps),
        "InputSHA256": before,
        "CompilerVersion": subprocess.check_output([cxx, "--version"], text=True),
        "FortranVersion": subprocess.check_output([fortran, "--version"], text=True),
        "CMakeVersion": subprocess.check_output([cmake, "--version"], text=True),
        "MPI": subprocess.check_output([str(args.mpirun), "--version"], text=True),
        "Environment": {k: os.environ.get(k) for k in ["PATH", "LD_LIBRARY_PATH", "LOADEDMODULES", "PBS_JOBID", "HOSTNAME"]},
        "Hostname": subprocess.check_output(["hostname"], text=True).strip(),
        "ConfigureCommand": configure, "RecordedConfigureCommand": recorded, "Jobs": args.jobs,
    })
    status = "incomplete"
    timing = {}
    try:
        started = time.time()
        print(shlex.join(configure), flush=True)
        with (build / "configure.log").open("w") as log:
            subprocess.run(configure, cwd=palace_build, check=True, stdout=log, stderr=subprocess.STDOUT)
        timing["ConfigureSeconds"] = time.time() - started
        started = time.time()
        build_command = [cmake, "--build", str(palace_build), "--target", "palace", "-j", str(args.jobs)]
        print(shlex.join(build_command), flush=True)
        with (build / "build.log").open("w") as log:
            subprocess.run(build_command, check=True, stdout=log, stderr=subprocess.STDOUT)
        timing["BuildSeconds"] = time.time() - started
        candidates = sorted(palace_build.glob("palace-*.bin")) + [palace_build / "palace"]
        executable = next(p for p in candidates if p.is_file() and os.access(p, os.X_OK))
        digest = sha256(executable)
        binary_root = (args.binary_root or build.parent).resolve()
        binary = binary_root / f"palace-{manifest['HEAD'][:7]}-{digest[:12]}.bin"
        with binary.open("xb") as stream:
            stream.write(executable.read_bytes())
        binary.chmod(0o555)
        ldd = subprocess.check_output(["ldd", str(binary)], text=True)
        dynamic = {t: sha256(t) for t in ldd.split() if t.startswith("/") and Path(t).is_file()}
        compile_commands = json.loads((palace_build / "compile_commands.json").read_text())
        save(build, "binary.json", {
            "Path": str(binary), "SHA256": digest, "BuiltFrom": str(executable), "HEAD": manifest["HEAD"],
            "Describe": describe, "Ldd": ldd, "DynamicDependenciesSHA256": dynamic,
            "TranslationUnits": len(compile_commands), "Timing": timing,
        })
        status = "complete"
    finally:
        changed = [str(p) for p in inputs if not p.exists() or sha256(p) != before[str(p)]]
        changed_sources = [p for p, h in manifest["Files"].items() if sha256(source / p) != h]
        save(build, "provenance-after.json", {"Status": status, "ChangedInputs": changed,
                                              "ChangedSources": changed_sources, "Timing": timing})
        if changed or changed_sources:
            raise SystemExit(f"inputs changed during the build: {changed[:3]} {changed_sources[:3]}")
    print(json.dumps(json.loads((build / "binary.json").read_text()) | {"Ldd": "..."}, indent=1))


if __name__ == "__main__":
    main()
