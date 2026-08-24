#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run the fixed-mesh p=1 corrected single-transmon eigenmode benchmark."""

import argparse
import csv
import difflib
import json
import math
import os
import platform
import shutil
import statistics
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
BENCHMARK_ROOT = HERE.parent
TRANSMON = BENCHMARK_ROOT.parent
ROOT = TRANSMON.parent.parent
sys.path.insert(0, str(BENCHMARK_ROOT))
import run_surface_response_preflight as common  # noqa: E402


NAME = "transmon-surface-response-corrected-eigenmode-p1"
DEFAULT_LIBRARY = HERE / "library/process-library.json"
DEFAULT_BASELINE = HERE / "corrected-eigenmode-baseline.json"
SNAPSHOT_FILES = (
    "eig.csv",
    "surface-Q.csv",
    "surface-Q-corrected.csv",
    "surface-response-confidence.csv",
)
INTEGER_COLUMNS = {
    "m",
    "exc",
    "interface",
    "boundary-law parameters verified",
    "confidence pass",
    "self-consistent confidence pass",
}


def library_files(library):
    directory = library.parent
    return {
        str(path.relative_to(directory)): common.sha256_file(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def fixture_document(source, mesh, library):
    generation = HERE / "generation.json"
    return {
        "ConfigTemplate": common.fixture(
            source, "examples/transmon/transmon_surface_coarse.json"
        ),
        "Mesh": common.fixture(
            mesh, "examples/transmon/mesh/transmon_surface_p1.msh2"
        ),
        "ProcessLibrary": {
            "Path": "repo://examples/transmon/benchmark/corrected-p1/library/process-library.json",
            "Files": library_files(library),
        },
        "Generation": common.fixture(
            generation,
            "examples/transmon/benchmark/corrected-p1/generation.json",
        ),
    }


def prepare_config(source, library, output, postpro):
    config = json.loads(source.read_text())
    mesh = Path(config["Model"]["Mesh"])
    if not mesh.is_absolute():
        mesh = source.parent / mesh
    mesh = mesh.resolve(strict=True)
    config["Problem"]["Verbose"] = 1
    config["Problem"]["Output"] = str(postpro.resolve())
    config["Problem"]["OutputFormats"] = {
        "Paraview": False,
        "GridFunction": False,
    }
    config["Model"]["Mesh"] = str(mesh)
    config["Model"]["Refinement"] = {
        "MaxIts": 0,
        "UniformLevels": 0,
        "SerialUniformLevels": 0,
        "SaveAdaptIterations": False,
        "SaveAdaptMesh": False,
    }
    config["Solver"]["Order"] = 1
    config["Solver"]["Eigenmode"]["N"] = 1
    config["Solver"]["Eigenmode"]["Save"] = 0
    matching_radius = float(json.loads(library.read_text())["MatchingRadius"])
    for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        if "EdgeRefinement" in dielectric:
            raise ValueError("The corrected p=1 benchmark forbids EdgeRefinement")
        dielectric["EdgeDistances"] = [matching_radius]
        dielectric["LocalizeEdgeEnergy"] = False
        dielectric["SaveLocalEdgeEnergy"] = False
    config["Solver"]["SurfaceResponseCorrection"] = {
        "Library": str(library.resolve(strict=True)),
        "TargetInterfaces": [1, 2, 3],
        "UnmatchedPolicy": "Warn",
        "SolveTol": 1.0e-6,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(config, indent=2) + "\n")
    return mesh


def read_csv(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        result = []
        for row in reader:
            result.append({key.strip(): float(value) for key, value in row.items()})
        return result


def output_snapshot(postpro):
    return {name: read_csv(postpro / name) for name in SNAPSHOT_FILES}


def workload(metadata):
    response = metadata["SurfaceResponse"]
    correction = response["Correction"]
    interpolation = response["Interpolation"]
    return {
        "MeshElements": int(metadata["Problem"]["MeshElements"]),
        "DegreesOfFreedom": int(metadata["Problem"]["DegreesOfFreedom"]),
        "Models": int(correction["Models"]),
        "Patches": int(correction["Patches"]),
        "TraceCoefficients": int(correction["TraceCoefficients"]),
        "ContourLines": int(correction["ContourLines"]["Total"]),
        "PointQueries": int(interpolation["PointQueries"]["Total"]),
        "StencilNonzeros": int(interpolation["StencilNonzeros"]["Total"]),
    }


def compare_snapshots(expected, observed, relative_tolerance, absolute_tolerance):
    failures = []
    if expected.keys() != observed.keys():
        return ["snapshot file set changed"]
    for filename in expected:
        if len(expected[filename]) != len(observed[filename]):
            failures.append(f"{filename}: row count changed")
            continue
        for row_index, (expected_row, observed_row) in enumerate(
            zip(expected[filename], observed[filename])
        ):
            if expected_row.keys() != observed_row.keys():
                failures.append(f"{filename}:{row_index}: columns changed")
                continue
            for column in expected_row:
                first = expected_row[column]
                second = observed_row[column]
                if column in INTEGER_COLUMNS:
                    matches = first == second
                else:
                    matches = math.isclose(
                        first,
                        second,
                        rel_tol=relative_tolerance,
                        abs_tol=absolute_tolerance,
                    )
                if not matches:
                    failures.append(
                        f"{filename}:{row_index}:{column}: {second} != {first}"
                    )
    return failures


def verify_fixtures(expected, source, mesh, library):
    observed = fixture_document(source, mesh, library)
    if observed != expected:
        difference = "".join(
            difflib.unified_diff(
                common.canonical_json(expected, True).splitlines(keepends=True),
                common.canonical_json(observed, True).splitlines(keepends=True),
                fromfile="checked-in fixtures",
                tofile="observed fixtures",
            )
        )
        raise RuntimeError("Corrected benchmark fixtures changed:\n" + difference)


def audit_run(config_path, postpro, ranks, mesh, library):
    config = json.loads(config_path.read_text())
    if config["Model"]["Refinement"]["MaxIts"] != 0:
        raise RuntimeError("Corrected benchmark enables AMR")
    if config["Solver"]["Order"] != 1:
        raise RuntimeError("Corrected benchmark must use p=1")
    if Path(config["Model"]["Mesh"]).resolve() != mesh:
        raise RuntimeError("Corrected benchmark mesh changed")
    if Path(config["Solver"]["SurfaceResponseCorrection"]["Library"]).resolve() != library:
        raise RuntimeError("Corrected benchmark library changed")
    for dielectric in config["Boundaries"]["Postprocessing"]["Dielectric"]:
        if "EdgeRefinement" in dielectric:
            raise RuntimeError("Corrected benchmark enables EdgeRefinement")
    metadata = json.loads((postpro / "palace.json").read_text())
    if int(metadata["Problem"]["MPISize"]) != ranks:
        raise RuntimeError(
            f"Requested {ranks} ranks, Palace recorded {metadata['Problem']['MPISize']}"
        )
    return metadata


def run_one(
    palace,
    launcher,
    launcher_args,
    numproc_flag,
    source,
    library,
    output,
    ranks,
    repeat,
    timeout,
):
    run_dir = output / f"np{ranks}-run{repeat}"
    run_dir.mkdir(parents=True, exist_ok=False)
    postpro = run_dir / "postpro"
    config_path = run_dir / "config.json"
    mesh = prepare_config(source, library, config_path, postpro)
    command = [
        str(launcher),
        *launcher_args,
        numproc_flag,
        str(ranks),
        str(palace),
        "--serial",
        str(config_path),
    ]
    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = "1"
    environment["OMP_DYNAMIC"] = "FALSE"
    return_code, wall_seconds = common.run_process(
        command,
        ROOT,
        environment,
        timeout,
        run_dir / "stdout.log",
        run_dir / "stderr.log",
    )
    if return_code:
        raise RuntimeError(f"Corrected benchmark failed; see {run_dir}")
    metadata = audit_run(config_path, postpro, ranks, mesh, library)
    return {
        "Ranks": ranks,
        "Repeat": repeat,
        "WallSeconds": wall_seconds,
        "Workload": workload(metadata),
        "Snapshot": output_snapshot(postpro),
        "Metrics": common.extract_metrics(metadata),
        "CorrectedSolver": metadata["SurfaceResponse"].get("CorrectedSolver", {}),
        "RunDirectory": str(run_dir.relative_to(output)),
    }


def baseline_document(runs, source, mesh, library):
    by_rank = {}
    for run in runs:
        if run["Repeat"] == 1:
            by_rank[str(run["Ranks"])] = {
                "Workload": run["Workload"],
                "Snapshot": run["Snapshot"],
            }
    return {
        "SchemaVersion": 1,
        "Benchmark": NAME,
        "Qualified": False,
        "FEMOrder": 1,
        "Fixtures": fixture_document(source, mesh, library),
        "RelativeTolerance": 5.0e-4,
        "AbsoluteTolerance": 1.0e-20,
        "Ranks": by_rank,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--palace", type=Path, default=ROOT / "build/bin/palace")
    parser.add_argument(
        "--launcher",
        type=Path,
        default=Path(shutil.which("mpiexec") or shutil.which("mpirun") or "mpiexec"),
    )
    parser.add_argument("--launcher-arg", action="append", default=[])
    parser.add_argument("--numproc-flag", default="-n")
    parser.add_argument(
        "--config", type=Path, default=TRANSMON / "transmon_surface_coarse.json"
    )
    parser.add_argument("--library", type=Path, default=DEFAULT_LIBRARY)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 4])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--timeout-seconds", type=float, default=600.0)
    parser.add_argument("--write-baseline-candidate", type=Path)
    args = parser.parse_args()
    if args.repeats < 1 or any(rank < 1 for rank in args.ranks):
        parser.error("ranks and repeats must be positive")

    palace = args.palace.expanduser().resolve(strict=True)
    launcher = args.launcher.expanduser().resolve(strict=True)
    source = args.config.expanduser().resolve(strict=True)
    library = args.library.expanduser().resolve(strict=True)
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
    if args.write_baseline_candidate:
        if sorted(args.ranks) != [1, 2, 4]:
            parser.error("baseline candidates require exactly --ranks 1 2 4")
        if args.write_baseline_candidate.resolve() == DEFAULT_BASELINE.resolve():
            parser.error("write a candidate outside the checked-in baseline path")
    else:
        checked_in = json.loads(args.baseline.resolve(strict=True).read_text())
        verify_fixtures(checked_in["Fixtures"], source, mesh, library)

    runs = []
    for ranks in args.ranks:
        for repeat in range(1, args.repeats + 1):
            print(f"Running corrected p=1 transmon: ranks={ranks}, repeat={repeat}")
            run = run_one(
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
            runs.append(run)
            if checked_in:
                expected = checked_in["Ranks"].get(str(ranks))
                if expected is None:
                    raise RuntimeError(f"No checked-in baseline for {ranks} ranks")
                if run["Workload"] != expected["Workload"]:
                    raise RuntimeError(
                        f"Workload changed for {ranks} ranks:\n"
                        + json.dumps(
                            {
                                "Expected": expected["Workload"],
                                "Observed": run["Workload"],
                            },
                            indent=2,
                        )
                    )
                failures = compare_snapshots(
                    expected["Snapshot"],
                    run["Snapshot"],
                    checked_in["RelativeTolerance"],
                    checked_in["AbsoluteTolerance"],
                )
                if failures:
                    raise RuntimeError(
                        "Numerical baseline changed:\n" + "\n".join(failures[:50])
                    )

    candidate = baseline_document(runs, source, mesh, library)
    if args.write_baseline_candidate:
        destination = args.write_baseline_candidate.expanduser().resolve()
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(common.canonical_json(candidate, True))
        print(f"Wrote baseline candidate: {destination}")

    results = {
        "SchemaVersion": 1,
        "Benchmark": NAME,
        "Qualified": False,
        "Host": {
            "Platform": platform.platform(),
            "Processor": platform.processor(),
            "Python": platform.python_version(),
        },
        "Environment": {"OMP_NUM_THREADS": "1", "OMP_DYNAMIC": "FALSE"},
        "Runs": [
            {
                key: value
                for key, value in run.items()
                if key not in {"Snapshot"}
            }
            for run in runs
        ],
        "SummaryByRanks": {},
    }
    for ranks in args.ranks:
        samples = [run["WallSeconds"] for run in runs if run["Ranks"] == ranks]
        results["SummaryByRanks"][str(ranks)] = {
            "WallSeconds": samples,
            "MedianWallSeconds": statistics.median(samples),
        }
    (output / "benchmark-results.json").write_text(
        common.canonical_json(results, True)
    )
    print(f"Benchmark results: {output / 'benchmark-results.json'}")


if __name__ == "__main__":
    main()
