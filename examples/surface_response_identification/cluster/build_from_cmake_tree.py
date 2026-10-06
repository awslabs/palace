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
        [--mfem-patch extern/patch/mfem/<diff> ...] [--mfem-test TAG --mfem-test-ranks 4]

The vendored MFEM patches the frozen tree carries (`SRC/extern/patch/mfem/*.diff`, the set
`cmake/ExternalMFEM.cmake` applies in a superbuild) are ALWAYS applied: MFEM is rebuilt with
them on top of the dependency prefix's already patched MFEM checkout (`PREFIX/extern/mfem`,
copied to `BUILD/mfem-src`, the patches applied with `git apply`), configured with the prefix's
recorded MFEM options (`PREFIX/extern/mfem-cmake/tmp/mfem-cfgcmd.txt`) re-pointed to install
into BUILD; the palace configure then uses `MFEM_DIR=BUILD` and the build fails unless palace
resolved the rebuilt library. `--mfem-patch` names the set explicitly (repeatable); a list that
omits a frozen patch is refused, so no binary can silently link the prefix's MFEM without a
patch its source tree depends on. `binary.json` records the applied patch set (`MFEM.Patches`,
empty when the frozen tree carries none) and the MFEM library palace linked. `--mfem-test`
builds MFEM's `punit_tests` and runs the given Catch2 filter under `mpirun -n RANKS` (the build
fails if the tests fail); the output is `BUILD/mfem-punit-tests.log` and the summary is recorded
in `binary.json`.

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


def recorded_configure_arguments(deps, step="palace-cmake/tmp/palace-cfgcmd.txt"):
    """The superbuild's recorded ExternalProject configure command for a step (a cmake list)."""
    text = (deps / step).read_text()
    match = re.search(r"cmd='(.*)'", text, re.S)
    if not match:
        raise SystemExit(f"no configure command in {deps}/{step}")
    arguments = match.group(1).replace("$<SEMICOLON>", "\x00").split(";")
    return [a.replace("\x00", ";") for a in arguments]


def run_logged(command, log_path, cwd=None):
    print(shlex.join(command), flush=True)
    with log_path.open("w") as log:
        subprocess.run(command, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)


MFEM_PATCH_DIRECTORY = "extern/patch/mfem"


def frozen_mfem_patches(source):
    """The vendored MFEM patches the frozen tree carries (what cmake/ExternalMFEM.cmake applies)."""
    return sorted(str(p.relative_to(source)) for p in (source / MFEM_PATCH_DIRECTORY).glob("*.diff"))


def resolve_mfem_patches(source, manifest_files, requested):
    """The MFEM patches to apply: the frozen tree's vendored set, or an explicit list that covers it.

    A binary built from a tree that carries extern/patch/mfem/*.diff must not link the dependency
    prefix's MFEM without them (a superbuild of the same tree would apply them), so an empty or
    partial --mfem-patch list never silently drops a frozen patch; every applied patch must be in
    the hash-bound manifest.
    """
    frozen = frozen_mfem_patches(source)
    patches = [str(Path(p)) for p in requested] if requested else frozen
    omitted = [p for p in frozen if p not in patches]
    if omitted:
        raise SystemExit(f"--mfem-patch omits vendored MFEM patches of the frozen tree {omitted}: "
                         "the binary would link an MFEM without them (pass them too)")
    unlisted = [p for p in patches if p not in manifest_files]
    if unlisted:
        raise SystemExit(f"MFEM patches not in the frozen manifest: {unlisted}")
    return patches


def rebuild_mfem(source, build, deps, patches, cmake, jobs, timing):
    """Copy the prefix's MFEM checkout, apply the frozen vendored patches, build + install into BUILD."""
    mfem_source = build / "mfem-src"
    mfem_build = build / "mfem-build"
    started = time.time()
    shutil.copytree(deps / "extern/mfem", mfem_source, symlinks=True)
    applied = []
    for relative in patches:
        patch = source / relative
        subprocess.run(["git", "apply", "--check", str(patch)], cwd=mfem_source, check=True)
        subprocess.run(["git", "apply", str(patch)], cwd=mfem_source, check=True)
        applied.append({"Patch": relative, "SHA256": sha256(patch)})
    diffstat = subprocess.check_output(["git", "diff", "--stat"], cwd=mfem_source, text=True)
    recorded = recorded_configure_arguments(deps, "extern/mfem-cmake/tmp/mfem-cfgcmd.txt")
    configure = [cmake, str(mfem_source)]
    for argument in recorded[2:]:
        if argument.startswith("-DCMAKE_INSTALL_PREFIX="):
            argument = f"-DCMAKE_INSTALL_PREFIX={build}"
        configure.append(argument)
    mfem_build.mkdir()
    run_logged(configure, build / "mfem-configure.log", cwd=mfem_build)
    run_logged([cmake, "--build", str(mfem_build), "--target", "mfem", "-j", str(jobs)], build / "mfem-build.log")
    run_logged([cmake, "--install", str(mfem_build)], build / "mfem-install.log")
    timing["MFEMBuildSeconds"] = time.time() - started
    return {"Source": str(mfem_source), "CopiedFrom": str(deps / "extern/mfem"), "Patches": applied,
            "DiffStat": diffstat, "ConfigureCommand": configure, "Library": str(build / "lib/libmfem.a"),
            "LibrarySHA256": sha256(build / "lib/libmfem.a")}


def run_mfem_unit_tests(build, mpirun, ranks, test_filter, cmake, jobs, timing):
    """Build MFEM's punit_tests and run one Catch2 filter from its data-relative working directory."""
    mfem_build = build / "mfem-build"
    started = time.time()
    # punit_tests reads meshes from ../../data of its working directory (copy_data fills BUILD/mfem-build/data).
    run_logged([cmake, "--build", str(mfem_build), "--target", "punit_tests", "copy_data", "-j", str(jobs)],
               build / "mfem-build-punit-tests.log")
    command = [str(mpirun), "-n", str(ranks), str(mfem_build / "tests/unit/punit_tests"), test_filter]
    print(shlex.join(command), flush=True)
    log_path = build / "mfem-punit-tests.log"
    with log_path.open("w") as log:
        completed = subprocess.run(command, cwd=mfem_build / "tests/unit", stdout=log, stderr=subprocess.STDOUT)
    timing["MFEMUnitTestSeconds"] = time.time() - started
    # Catch2 colours its summary lines with ANSI escapes even when piped.
    plain = re.sub(r"\x1b\[[0-9;]*m", "", log_path.read_text())
    summary = [line for line in plain.splitlines()
               if line.startswith(("test cases:", "assertions:", "All tests passed"))]
    result = {"Command": command, "Ranks": ranks, "Filter": test_filter, "ReturnCode": completed.returncode,
              "Summary": summary, "Log": str(log_path)}
    if completed.returncode != 0:
        raise SystemExit(f"MFEM unit tests failed (rc {completed.returncode}): {summary} see {log_path}")
    return result


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
    parser.add_argument("--mfem-patch", action="append", default=[], metavar="RELATIVE_DIFF",
                        help="frozen-tree path of a vendored MFEM patch to apply on top of the prefix's MFEM "
                             "(default: every extern/patch/mfem/*.diff of the frozen tree; a list that omits one is refused)")
    parser.add_argument("--mfem-test", metavar="CATCH2_FILTER",
                        help="run MFEM punit_tests with this filter (requires a rebuilt MFEM, i.e. at least one patch)")
    parser.add_argument("--mfem-test-ranks", type=int, default=4)
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
    patches = resolve_mfem_patches(source, manifest["Files"], args.mfem_patch)
    if args.mfem_test and not patches:
        parser.error("--mfem-test needs a rebuilt MFEM: the frozen tree carries no extern/patch/mfem/*.diff and no --mfem-patch was given")
    print(f"MFEM patches: {patches or 'none (the dependency prefix MFEM is used as is)'}", flush=True)

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
        elif argument.startswith("-DMFEM_DIR=") and patches:
            argument = f"-DMFEM_DIR={build}"
        elif argument.startswith("-DCMAKE_PREFIX_PATH=") and patches:
            # find_library(MFEM_LIBRARY ... HINTS ${MFEM_DIR}/lib) searches CMAKE_PREFIX_PATH before HINTS,
            # so the dependency prefix's unpatched libmfem.a would win: this prefix goes first.
            argument = f"-DCMAKE_PREFIX_PATH={build};" + argument.split("=", 1)[1]
        configure.append(argument)
    if patches:
        configure.append(f"-DMFEM_LIBRARY={build}/lib/libmfem.a")
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
        "MFEMPatches": patches, "RequestedMFEMPatches": args.mfem_patch,
    })
    status = "incomplete"
    timing = {}
    # Always recorded: the applied patch set (empty = the prefix's MFEM as is) and the library palace linked.
    mfem = {"Patches": [], "Library": None}
    try:
        if patches:
            mfem = rebuild_mfem(source, build, deps, patches, cmake, args.jobs, timing)
            if args.mfem_test:
                mfem["UnitTests"] = run_mfem_unit_tests(build, args.mpirun, args.mfem_test_ranks, args.mfem_test,
                                                        cmake, args.jobs, timing)
        started = time.time()
        print(shlex.join(configure), flush=True)
        with (build / "configure.log").open("w") as log:
            subprocess.run(configure, cwd=palace_build, check=True, stdout=log, stderr=subprocess.STDOUT)
        timing["ConfigureSeconds"] = time.time() - started
        cache = (palace_build / "CMakeCache.txt").read_text()
        found = re.search(r"^MFEM_LIBRARY:\w+=(.*)$", cache, re.M)
        if patches:
            expected = str(build / "lib/libmfem.a")
            if not found or found.group(1) != expected:
                raise SystemExit(f"palace configured against {found and found.group(1)}, not the rebuilt {expected}")
        mfem["PalaceMFEMLibrary"] = found.group(1) if found else None
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
            "TranslationUnits": len(compile_commands), "Timing": timing, "MFEM": mfem,
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
