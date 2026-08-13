#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run the fixed-mesh transmon surface-response preflight benchmark."""

import argparse
import difflib
import hashlib
import json
import os
import platform
import shutil
import signal
import statistics
import subprocess
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
TRANSMON = HERE.parent
ROOT = TRANSMON.parent.parent
DEFAULT_BASELINE = HERE / "surface-response-preflight-baseline.json"
BENCHMARK_NAME = "transmon-surface-response-preflight-fixed-mesh"
REPO_URI = "repo://examples/transmon/benchmark/transmon_surface_process_seed.json"


def sha256_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_manifest(manifest, expected_library=None):
    """Remove environment-specific and observational fields from a preflight manifest."""
    result = json.loads(json.dumps(manifest))
    observed = Path(result["Library"]["Path"])
    if not observed.is_absolute():
        raise ValueError("Preflight manifest library path must be absolute")
    observed = observed.resolve(strict=True)
    if expected_library is not None and observed != expected_library.resolve(strict=True):
        raise ValueError(
            f"Preflight used {observed}, expected {expected_library.resolve(strict=True)}"
        )
    try:
        observed.relative_to(ROOT.resolve(strict=True))
    except ValueError as error:
        raise ValueError("Preflight library must be inside the repository") from error
    result["Library"]["Path"] = REPO_URI
    # Observability is intentionally excluded from the numerical/geometry correctness
    # baseline so Phase 0 can add counters without changing the canonical workload hash.
    result.pop("Statistics", None)
    return result


def canonical_json(value, pretty=False):
    if pretty:
        return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def canonical_sha256(manifest, expected_library=None):
    payload = canonical_json(canonical_manifest(manifest, expected_library))
    return hashlib.sha256(payload.encode()).hexdigest()


def manifest_summary(manifest):
    by_topology = {}
    for requirement in manifest["Requirements"]:
        entry = by_topology.setdefault(
            requirement["Topology"],
            {"Records": 0, "Count": 0, "TotalEdgeLength": 0.0},
        )
        entry["Records"] += 1
        entry["Count"] += int(requirement["Count"])
        entry["TotalEdgeLength"] += float(requirement["TotalEdgeLength"])
    return {
        "RequirementRecords": len(manifest["Requirements"]),
        "Summary": manifest["Summary"],
        "ByTopology": dict(sorted(by_topology.items())),
    }


def fixture(path, relative):
    return {"Path": f"repo://{relative}", "Sha256": sha256_file(path)}


def baseline_document(manifest, source, mesh, library):
    canonical = canonical_manifest(manifest, library)
    return {
        "SchemaVersion": 1,
        "Benchmark": BENCHMARK_NAME,
        "Fixtures": {
            "ConfigTemplate": fixture(
                source, "examples/transmon/transmon_surface_coarse.json"
            ),
            "Mesh": fixture(
                mesh, "examples/transmon/mesh/transmon_surface_p1.msh2"
            ),
            "ProcessLibrary": fixture(
                library,
                "examples/transmon/benchmark/transmon_surface_process_seed.json",
            ),
        },
        "CanonicalManifestSha256": hashlib.sha256(
            canonical_json(canonical).encode()
        ).hexdigest(),
        "CanonicalManifest": canonical,
        **manifest_summary(manifest),
    }


def prepare_config(source, library, output, postpro):
    """Prepare a bounded preflight: no AMR, edge refinement, or field output."""
    config = json.loads(source.read_text())
    mesh = Path(config["Model"]["Mesh"])
    if not mesh.is_absolute():
        mesh = source.parent / mesh
    mesh = mesh.resolve(strict=True)
    config["Model"]["Mesh"] = str(mesh)
    refinement = config["Model"]["Refinement"] = {
        "MaxIts": 0,
        "UniformLevels": 0,
        "SerialUniformLevels": 0,
        "SaveAdaptIterations": False,
        "SaveAdaptMesh": False,
    }
    assert refinement["MaxIts"] == 0
    config["Problem"]["Output"] = str(postpro.resolve())
    config["Problem"]["OutputFormats"] = {
        "Paraview": False,
        "GridFunction": False,
    }
    config["Solver"]["Eigenmode"]["Save"] = 0
    matching_radius = float(json.loads(library.read_text())["MatchingRadius"])
    for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        if "EdgeRefinement" in dielectric:
            raise ValueError("The fixed-mesh benchmark forbids EdgeRefinement")
        dielectric["EdgeDistances"] = [matching_radius]
        dielectric["LocalizeEdgeEnergy"] = False
        dielectric["SaveLocalEdgeEnergy"] = False
    config["Solver"]["SurfaceResponseCorrection"] = {
        "Library": str(library.resolve(strict=True)),
        "TargetInterfaces": [1, 2, 3],
        "UnmatchedPolicy": "Warn",
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(config, indent=2) + "\n")
    return config, mesh


def audit_resolved_config(path, expected_ranks, mesh, library):
    config = json.loads(path.read_text())
    refinement = config["Model"]["Refinement"]
    if refinement["MaxIts"] != 0 or refinement["UniformLevels"] != 0:
        raise RuntimeError("Resolved benchmark configuration enables refinement")
    for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        if "EdgeRefinement" in dielectric:
            raise RuntimeError("Resolved benchmark configuration enables EdgeRefinement")
        if dielectric["LocalizeEdgeEnergy"] or dielectric["SaveLocalEdgeEnergy"]:
            raise RuntimeError("Resolved benchmark configuration enables local output")
    if Path(config["Model"]["Mesh"]).resolve() != mesh:
        raise RuntimeError("Resolved benchmark mesh changed")
    response = config["Solver"]["SurfaceResponseCorrection"]
    if Path(response["Library"]).resolve() != library:
        raise RuntimeError("Resolved benchmark process library changed")
    output = config["Problem"]["OutputFormats"]
    if output["Paraview"] or output["GridFunction"]:
        raise RuntimeError("Resolved benchmark configuration enables field output")
    # The requested MPI size is recorded in palace.json, not the resolved configuration.
    if expected_ranks < 1:
        raise RuntimeError("Invalid expected MPI rank count")


def extract_metrics(metadata):
    required = (
        "ElapsedTime",
        "PeakMemoryMegabytes",
        "PeakNodeMemoryMegabytes",
        "Problem",
    )
    missing = [key for key in required if key not in metadata]
    if missing:
        raise RuntimeError(f"Preflight metadata is missing {missing}")
    return {key: metadata[key] for key in required}


def run_process(command, cwd, environment, timeout, stdout_path, stderr_path):
    start = time.perf_counter()
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        process = subprocess.Popen(
            command,
            cwd=cwd,
            env=environment,
            text=True,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        try:
            return_code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            raise TimeoutError(f"Benchmark timed out after {timeout} seconds")
    return return_code, time.perf_counter() - start


def build_command(launcher, launcher_args, numproc_flag, ranks, palace, config):
    return [
        str(launcher),
        *launcher_args,
        numproc_flag,
        str(ranks),
        str(palace),
        "--serial",
        "--surface-response-preflight",
        str(config),
    ]


def run_one(
    palace,
    launcher,
    launcher_args,
    numproc_flag,
    source,
    library,
    root,
    ranks,
    repeat,
    timeout,
):
    run_dir = root / f"np{ranks}-run{repeat}"
    postpro = run_dir / "postpro"
    config_path = run_dir / "config.json"
    run_dir.mkdir(parents=True, exist_ok=False)
    _, mesh = prepare_config(source, library, config_path, postpro)
    command = build_command(
        launcher, launcher_args, numproc_flag, ranks, palace, config_path
    )
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"
    return_code, wall_seconds = run_process(
        command,
        ROOT,
        environment,
        timeout,
        run_dir / "stdout.log",
        run_dir / "stderr.log",
    )
    if return_code:
        raise RuntimeError(
            f"Preflight failed for {ranks} ranks (see {run_dir / 'stderr.log'})"
        )
    manifest_path = postpro / "surface-response-requirements.json"
    resolved_path = postpro / "config_resolved.json"
    metadata_path = postpro / "palace.json"
    audit_resolved_config(resolved_path, ranks, mesh, library)
    metadata = json.loads(metadata_path.read_text())
    if int(metadata["Problem"]["MPISize"]) != ranks:
        raise RuntimeError(
            f"Requested {ranks} MPI ranks, Palace recorded "
            f"{metadata['Problem']['MPISize']}"
        )
    manifest = json.loads(manifest_path.read_text())
    canonical = canonical_manifest(manifest, library)
    (run_dir / "canonical-manifest.json").write_text(canonical_json(canonical, True))
    return {
        "Ranks": ranks,
        "Repeat": repeat,
        "WallSeconds": wall_seconds,
        "CanonicalManifestSha256": hashlib.sha256(
            canonical_json(canonical).encode()
        ).hexdigest(),
        "CanonicalManifest": canonical,
        "Manifest": manifest,
        "Metrics": extract_metrics(metadata),
        "RunDirectory": str(run_dir.relative_to(root)),
    }


def verify_fixture_hashes(baseline, source, mesh, library):
    observed = {
        "ConfigTemplate": sha256_file(source),
        "Mesh": sha256_file(mesh),
        "ProcessLibrary": sha256_file(library),
    }
    for name, digest in observed.items():
        expected = baseline["Fixtures"][name]["Sha256"]
        if digest != expected:
            raise RuntimeError(
                f"{name} fixture changed ({digest} != {expected}); review a new baseline"
            )


def check_baseline(run, baseline):
    observed = run["CanonicalManifest"]
    expected = baseline["CanonicalManifest"]
    if observed == expected:
        return
    difference = "".join(
        difflib.unified_diff(
            canonical_json(expected, True).splitlines(keepends=True),
            canonical_json(observed, True).splitlines(keepends=True),
            fromfile="checked-in baseline",
            tofile="observed manifest",
        )
    )
    raise RuntimeError("Preflight manifest changed:\n" + difference)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--palace", type=Path, default=ROOT / "build/bin/palace")
    parser.add_argument(
        "--config", type=Path, default=TRANSMON / "transmon_surface_coarse.json"
    )
    parser.add_argument(
        "--library",
        type=Path,
        default=HERE / "transmon_surface_process_seed.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--launcher",
        type=Path,
        default=Path(shutil.which("mpiexec") or shutil.which("mpirun") or "mpiexec"),
    )
    parser.add_argument("--launcher-arg", action="append", default=[])
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout-seconds", type=float, default=300.0)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--write-baseline-candidate", type=Path)
    args = parser.parse_args()
    if args.repeats < 1 or any(rank < 1 for rank in args.ranks):
        parser.error("ranks and repeats must be positive")

    source = args.config.expanduser().resolve(strict=True)
    library = args.library.expanduser().resolve(strict=True)
    palace = args.palace.expanduser().resolve(strict=True)
    launcher = args.launcher.expanduser().resolve(strict=True)
    output = args.output.expanduser().resolve()
    if output.exists():
        parser.error(f"output directory already exists: {output}")
    output.mkdir(parents=True)
    source_config = json.loads(source.read_text())
    mesh = Path(source_config["Model"]["Mesh"])
    if not mesh.is_absolute():
        mesh = source.parent / mesh
    mesh = mesh.resolve(strict=True)

    checked_in = None
    if not args.write_baseline_candidate:
        checked_in = json.loads(args.baseline.resolve(strict=True).read_text())
        verify_fixture_hashes(checked_in, source, mesh, library)
    elif args.write_baseline_candidate.resolve() == DEFAULT_BASELINE.resolve():
        parser.error("write a candidate outside the checked-in baseline path")
    elif sorted(args.ranks) != [1, 2, 4]:
        parser.error("baseline candidates require exactly --ranks 1 2 4")

    version = subprocess.run(
        [str(palace), "--serial", "--version"],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=True,
    ).stdout.strip()
    runs = []
    for ranks in args.ranks:
        for repeat in range(1, args.repeats + 1):
            print(f"Running transmon preflight: ranks={ranks}, repeat={repeat}")
            runs.append(
                run_one(
                    palace,
                    launcher,
                    args.launcher_arg,
                    args.numproc_flag,
                    source,
                    library,
                    output,
                    ranks,
                    repeat,
                    args.timeout_seconds,
                )
            )

    first = runs[0]["CanonicalManifest"]
    for run in runs[1:]:
        if run["CanonicalManifest"] != first:
            raise RuntimeError(
                f"Rank/repeat-dependent manifest for ranks={run['Ranks']}, "
                f"repeat={run['Repeat']}"
            )
    baseline = baseline_document(runs[0]["Manifest"], source, mesh, library)
    if args.write_baseline_candidate:
        candidate = args.write_baseline_candidate.expanduser().resolve()
        candidate.parent.mkdir(parents=True, exist_ok=True)
        candidate.write_text(canonical_json(baseline, True))
        print(f"Wrote baseline candidate: {candidate}")
    else:
        for run in runs:
            check_baseline(run, checked_in)

    results = {
        "SchemaVersion": 1,
        "Benchmark": BENCHMARK_NAME,
        "Host": {
            "Platform": platform.platform(),
            "Processor": platform.processor(),
            "Python": platform.python_version(),
        },
        "Palace": {
            "Path": str(palace),
            "Sha256": sha256_file(palace),
            "Version": version,
        },
        "Launcher": {
            "Path": str(launcher),
            "Arguments": args.launcher_arg,
            "NumProcFlag": args.numproc_flag,
        },
        "Environment": {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"},
        "Baseline": {
            "CanonicalManifestSha256": baseline["CanonicalManifestSha256"],
            "RequirementRecords": baseline["RequirementRecords"],
            "Summary": baseline["Summary"],
            "ByTopology": baseline["ByTopology"],
        },
        "Ranks": {},
    }
    for ranks in args.ranks:
        rank_runs = [run for run in runs if run["Ranks"] == ranks]
        results["Ranks"][str(ranks)] = {
            "WallSeconds": [run["WallSeconds"] for run in rank_runs],
            "MedianWallSeconds": statistics.median(
                run["WallSeconds"] for run in rank_runs
            ),
            "Metrics": [run["Metrics"] for run in rank_runs],
            "RunDirectories": [run["RunDirectory"] for run in rank_runs],
        }
    result_path = output / "benchmark-results.json"
    result_path.write_text(canonical_json(results, True))
    print(f"Benchmark results: {result_path}")


if __name__ == "__main__":
    main()
