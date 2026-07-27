#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Generate and validate the controlled 3D zero-thickness CPW study."""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path

import validate
from mesh.generate_mesh import generate_mesh


ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / "study" / "manifest.json"

BASE_GEOMETRY = {
    "strip_half_width": 5.0,
    "gap": 3.0,
    "ground_outer": 80.0,
    "box_half_width": 160.0,
    "box_half_height": 160.0,
    "length": 10.0,
    "strip_half_intervals": 2,
    "gap_intervals": 2,
    "ground_intervals": 8,
    "outer_intervals": 4,
    "height_intervals": 10,
    "longitudinal_intervals": 3,
    "transverse_subdivision_factor": 1,
}

MESH_CASES = {
    "fixed_coarse": {},
    "fixed_medium": {
        "transverse_subdivision_factor": 2,
    },
    "fixed_fine": {
        "transverse_subdivision_factor": 4,
    },
    "large_box_fixed_ground": {
        "box_half_width": 320.0,
        "box_half_height": 320.0,
        "outer_intervals": 6,
        "height_intervals": 12,
    },
    "large_box_scaled_ground": {
        "ground_outer": 160.0,
        "box_half_width": 320.0,
        "box_half_height": 320.0,
        "ground_intervals": 12,
        "outer_intervals": 5,
        "height_intervals": 12,
    },
    "large_box_scaled_ground_medium": {
        "ground_outer": 160.0,
        "box_half_width": 320.0,
        "box_half_height": 320.0,
        "ground_intervals": 12,
        "outer_intervals": 5,
        "height_intervals": 12,
        "transverse_subdivision_factor": 2,
    },
    "xlarge_box_scaled_ground": {
        "ground_outer": 320.0,
        "box_half_width": 640.0,
        "box_half_height": 640.0,
        "ground_intervals": 17,
        "outer_intervals": 6,
        "height_intervals": 14,
    },
}

RUN_METHODS = {
    "fixed_coarse": ("standard", "singular_s1", "singular_s2", "singular_s3"),
    "fixed_medium": ("standard", "singular_s1", "singular_s2"),
    "fixed_fine": ("standard", "singular_s1"),
    "large_box_fixed_ground": ("singular_s1",),
    "large_box_scaled_ground": ("singular_s1",),
    "large_box_scaled_ground_medium": ("singular_s1",),
    "xlarge_box_scaled_ground": ("singular_s1",),
}

MESH_PARENTS = {
    "fixed_medium": "fixed_coarse",
    "fixed_fine": "fixed_medium",
    "large_box_scaled_ground_medium": "large_box_scaled_ground",
}


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_mesh_vertices(path):
    lines = path.read_text(encoding="ascii").splitlines()
    marker = lines.index("vertices")
    count = int(lines[marker + 1])
    dimension = int(lines[marker + 2])
    if dimension != 3:
        raise ValueError(f"Expected a 3D mesh: {path}")
    return [
        tuple(float(value) for value in line.split())
        for line in lines[marker + 3 : marker + 3 + count]
    ]


def verify_nested_vertices(coarse_path, fine_path):
    coarse = read_mesh_vertices(coarse_path)
    fine = read_mesh_vertices(fine_path)
    digits = 12
    fine_vertices = {
        tuple(round(value, digits) for value in point) for point in fine
    }
    unmatched = [
        point
        for point in coarse
        if tuple(round(value, digits) for value in point) not in fine_vertices
    ]
    if unmatched:
        raise ValueError(
            f"{len(unmatched)} vertices from {coarse_path} are absent in {fine_path}"
        )


def expected_mesh_counts(parameters):
    positive_x = sum(
        parameters[name]
        for name in (
            "strip_half_intervals",
            "gap_intervals",
            "ground_intervals",
            "outer_intervals",
        )
    )
    factor = parameters["transverse_subdivision_factor"]
    x_cells = 2 * positive_x * factor
    y_cells = 2 * parameters["height_intervals"] * factor
    z_cells = parameters["longitudinal_intervals"]
    return {
        "x_cells": x_cells,
        "y_cells": y_cells,
        "z_cells": z_cells,
        "tetrahedra": 6 * x_cells * y_cells * z_cells,
    }


def configure_run(base_standard, base_singular, case, method, mesh_path, abs_tol, rel_tol):
    singular = method.startswith("singular_s")
    config = copy.deepcopy(base_singular if singular else base_standard)
    run = f"{case}_{method}"
    config["Problem"]["Output"] = f"postpro/study/{run}"
    config["Problem"]["OutputFormats"]["Paraview"] = False
    config["Model"]["Mesh"] = mesh_path.relative_to(ROOT).as_posix()
    config["Solver"]["Electrostatic"]["Save"] = 0
    if singular:
        singular_order = int(method.removeprefix("singular_s"))
        options = config["Solver"]["SingularElements"]
        options["Order"] = singular_order
        options["AbsTol"] = abs_tol
        options["RelTol"] = rel_tol
    return run, config


def generate_study(abs_tol, rel_tol):
    mesh_dir = ROOT / "study" / "mesh"
    config_dir = ROOT / "study" / "config"
    mesh_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)

    with (ROOT / "cpw3d_standard.json").open(encoding="ascii") as stream:
        base_standard = json.load(stream)
    with (ROOT / "cpw3d_singular.json").open(encoding="ascii") as stream:
        base_singular = json.load(stream)

    manifest = {
        "version": 1,
        "length_unit_m": 1.0e-6,
        "reference": {
            "type": "infinite-domain conformal-map CPW",
            "strip_half_width": BASE_GEOMETRY["strip_half_width"],
            "ground_inner": BASE_GEOMETRY["strip_half_width"]
            + BASE_GEOMETRY["gap"],
            "substrate_permittivity": 11.45,
        },
        "quadrature": {
            "absolute_tolerance": abs_tol,
            "relative_tolerance": rel_tol,
            "maximum_subdivisions": base_singular["Solver"]["SingularElements"][
                "MaxSubdivisions"
            ],
        },
        "meshes": {},
        "runs": {},
    }

    for case, overrides in MESH_CASES.items():
        parameters = BASE_GEOMETRY | overrides
        mesh_path = mesh_dir / f"{case}.mesh"
        generate_mesh(mesh_path, **parameters)
        manifest["meshes"][case] = {
            "file": mesh_path.relative_to(ROOT).as_posix(),
            "sha256": file_sha256(mesh_path),
            "parameters_mesh_units": parameters,
            "counts": expected_mesh_counts(parameters),
        }
        if case in MESH_PARENTS:
            parent = MESH_PARENTS[case]
            parent_path = ROOT / manifest["meshes"][parent]["file"]
            verify_nested_vertices(parent_path, mesh_path)
            manifest["meshes"][case]["nested_parent"] = parent

        for method in RUN_METHODS[case]:
            run, config = configure_run(
                base_standard,
                base_singular,
                case,
                method,
                mesh_path,
                abs_tol,
                rel_tol,
            )
            config_path = config_dir / f"{run}.json"
            with config_path.open("w", encoding="ascii") as stream:
                json.dump(config, stream, indent=2)
                stream.write("\n")
            manifest["runs"][run] = {
                "mesh": case,
                "method": method,
                "config": config_path.relative_to(ROOT).as_posix(),
                "config_sha256": file_sha256(config_path),
                "output": config["Problem"]["Output"],
                "command_from_example_directory": [
                    "../../build/bin/palace",
                    "-np",
                    "1",
                    config_path.relative_to(ROOT).as_posix(),
                ],
            }

    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with MANIFEST.open("w", encoding="ascii") as stream:
        json.dump(manifest, stream, indent=2)
        stream.write("\n")
    print(f"Wrote {MANIFEST}")
    print(f"  meshes={len(manifest['meshes'])}, runs={len(manifest['runs'])}")


def validate_study(run_names):
    with MANIFEST.open(encoding="ascii") as stream:
        manifest = json.load(stream)
    reference_data = manifest["reference"]
    reference = validate.reference_capacitance_per_length(
        reference_data["strip_half_width"],
        reference_data["ground_inner"],
        reference_data["substrate_permittivity"],
    )
    length = BASE_GEOMETRY["length"] * manifest["length_unit_m"]
    selected = run_names or list(manifest["runs"])
    print(f"Infinite-domain reference C' = {reference:.12e} F/m")
    for run in selected:
        if run not in manifest["runs"]:
            raise ValueError(f"Unknown study run: {run}")
        postpro = ROOT / manifest["runs"][run]["output"]
        if not (postpro / "palace.json").is_file():
            print(f"{run}: missing output {postpro.relative_to(ROOT)}")
            continue
        validate.report_result(postpro, reference, length)
    if not run_names:
        summarize_study(manifest, reference, length)


def read_run_result(manifest, run, reference, length):
    if run not in manifest["runs"]:
        raise ValueError(f"Study manifest is missing required run: {run}")
    postpro = ROOT / manifest["runs"][run]["output"]
    required = (postpro / "palace.json", postpro / "terminal-C.csv")
    missing = [path.relative_to(ROOT) for path in required if not path.is_file()]
    if missing:
        raise ValueError(f"Study run {run} is missing output: {missing}")

    capacitance = validate.read_capacitance(postpro / "terminal-C.csv")
    metadata = validate.read_metadata(postpro)
    per_length = capacitance / length
    linear_solver = metadata["LinearSolver"]
    return {
        "run": run,
        "mesh_elements": int(metadata["Problem"]["MeshElements"]),
        "capacitance_per_length": per_length,
        "relative_error": per_length / reference - 1.0,
        "iterations": int(linear_solver["TotalIts"]),
        "converged": linear_solver.get("Converged"),
        "failed_solves": linear_solver.get("FailedSolves"),
    }


def richardson_summary(results, reference, refinement_ratio=2.0):
    if len(results) != 3:
        raise ValueError("Richardson extrapolation requires exactly three levels")
    values = [result["capacitance_per_length"] for result in results]
    coarse_difference = values[0] - values[1]
    fine_difference = values[1] - values[2]
    if (
        not all(math.isfinite(value) for value in values)
        or not math.isfinite(refinement_ratio)
        or refinement_ratio <= 1.0
        or coarse_difference * fine_difference <= 0.0
    ):
        raise ValueError("Richardson sequence is not finite and monotone")
    difference_ratio = abs(coarse_difference / fine_difference)
    if difference_ratio <= 1.0:
        raise ValueError("Richardson sequence does not show asymptotic convergence")
    observed_order = math.log(difference_ratio) / math.log(refinement_ratio)
    extrapolated = values[2] + (values[2] - values[1]) / (
        refinement_ratio**observed_order - 1.0
    )
    return {
        "observed_order": observed_order,
        "extrapolated": extrapolated,
        "relative_error": extrapolated / reference - 1.0,
    }


def convergence_label(result):
    if result["converged"] is True:
        return "yes"
    if result["converged"] is False:
        return f"NO ({result['failed_solves']} failed)"
    return "not recorded"


def summarize_study(manifest, reference, length):
    print("\nControlled study summary")
    sequences = {
        "standard p=1": (
            "fixed_coarse_standard",
            "fixed_medium_standard",
            "fixed_fine_standard",
        ),
        "singular s=1": (
            "fixed_coarse_singular_s1",
            "fixed_medium_singular_s1",
            "fixed_fine_singular_s1",
        ),
    }
    for label, runs in sequences.items():
        results = [
            read_run_result(manifest, run, reference, length) for run in runs
        ]
        richardson = richardson_summary(results, reference)
        counts = ", ".join(f"{result['mesh_elements']:,}" for result in results)
        errors = ", ".join(
            f"{abs(result['relative_error']):.4%}" for result in results
        )
        print(f"  {label}:")
        print(f"    tetrahedra = {counts}")
        print(f"    absolute reference errors = {errors}")
        print(
            f"    observed order = {richardson['observed_order']:.6f}, "
            f"Richardson C' = {richardson['extrapolated']:.12e} F/m, "
            f"reference error = {richardson['relative_error']:+.6%}"
        )

    order_runs = (
        "fixed_coarse_singular_s1",
        "fixed_coarse_singular_s2",
        "fixed_coarse_singular_s3",
    )
    print("  fixed coarse singular-order study:")
    for run in order_runs:
        result = read_run_result(manifest, run, reference, length)
        method = manifest["runs"][run]["method"]
        print(
            f"    {method}: error = {result['relative_error']:+.6%}, "
            f"iterations = {result['iterations']}, "
            f"converged = {convergence_label(result)}"
        )

    boundary_pairs = (
        (
            "double box at fixed ground",
            "fixed_coarse_singular_s1",
            "large_box_fixed_ground_singular_s1",
        ),
        (
            "extend ground in doubled box",
            "large_box_fixed_ground_singular_s1",
            "large_box_scaled_ground_singular_s1",
        ),
        (
            "double box and ground again",
            "large_box_scaled_ground_singular_s1",
            "xlarge_box_scaled_ground_singular_s1",
        ),
    )
    print("  finite-boundary sensitivity:")
    for label, baseline_run, changed_run in boundary_pairs:
        baseline = read_run_result(manifest, baseline_run, reference, length)
        changed = read_run_result(manifest, changed_run, reference, length)
        relative_change = (
            changed["capacitance_per_length"]
            / baseline["capacitance_per_length"]
            - 1.0
        )
        print(f"    {label}: {relative_change:+.6%}")


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate_parser = subparsers.add_parser("generate")
    generate_parser.add_argument("--abs-tol", type=float, default=1.0e-5)
    generate_parser.add_argument("--rel-tol", type=float, default=1.0e-5)
    validate_parser = subparsers.add_parser("validate")
    validate_parser.add_argument("runs", nargs="*")
    subparsers.add_parser("summarize")
    args = parser.parse_args()

    if args.command == "generate":
        if args.abs_tol <= 0.0 or args.rel_tol <= 0.0:
            raise ValueError("Study quadrature tolerances must be positive")
        generate_study(args.abs_tol, args.rel_tol)
    elif args.command == "validate":
        validate_study(args.runs)
    else:
        with MANIFEST.open(encoding="ascii") as stream:
            manifest = json.load(stream)
        reference_data = manifest["reference"]
        reference = validate.reference_capacitance_per_length(
            reference_data["strip_half_width"],
            reference_data["ground_inner"],
            reference_data["substrate_permittivity"],
        )
        length = BASE_GEOMETRY["length"] * manifest["length_unit_m"]
        print(f"Infinite-domain reference C' = {reference:.12e} F/m")
        summarize_study(manifest, reference, length)


if __name__ == "__main__":
    main()
