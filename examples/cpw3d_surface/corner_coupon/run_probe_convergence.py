#!/usr/bin/env python3

"""Run low-cost coupon probes across FEM orders and enforce convergence gates."""

import argparse
import json
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np

import compare_probe_convergence as comparison


KINDS = ("thin", "fabricated")
OUTPUT_FILES = ("domain-response-matrix.csv", "surface-response-matrix.csv")


def run(command):
    print("+ " + shlex.join(str(value) for value in command), flush=True)
    subprocess.run([str(value) for value in command], check=True)


def prepare_case(calibration, output, name, order):
    case = output / name
    case.mkdir(parents=True, exist_ok=True)
    for filename in ("probe-manifest.json", "process-library.json"):
        shutil.copy2(calibration / filename, case / filename)

    configs = []
    for kind in KINDS:
        source = calibration / f"probe-{kind}.json"
        config = json.loads(source.read_text())
        config["Solver"]["Order"] = order
        config["Problem"]["Output"] = str(case / "postpro" / f"probe-{kind}")
        destination = case / source.name
        destination.write_text(json.dumps(config, indent=2) + "\n")
        configs.append(destination)
    return case, configs


def complete(case):
    return all(
        (case / "postpro" / f"probe-{kind}" / filename).is_file()
        for kind in KINDS
        for filename in OUTPUT_FILES
    )


def configured_order(calibration):
    orders = {
        int(
            json.loads(
                (calibration / f"probe-{kind}.json").read_text()
            )["Solver"]["Order"]
        )
        for kind in KINDS
    }
    if len(orders) != 1:
        raise ValueError(f"{calibration} probe configs use inconsistent FEM orders")
    return orders.pop()


def metric(current, previous):
    norm = np.linalg.norm(current)
    matrix_change = np.linalg.norm(current - previous) / max(
        norm, np.finfo(float).tiny
    )
    return {
        "MatrixChangePercent": 100.0 * matrix_change,
        "WorstEnergyChangePercent":
            100.0 * comparison.worst_energy_change(current, previous),
    }


def compare_cases(previous, current, args, enforce):
    results = []
    previous_matrices = {
        (kind, quantity): matrix
        for kind, quantity, matrix in comparison.matrix_rows(previous)
    }
    passed = True
    for kind, quantity, matrix in comparison.matrix_rows(current):
        values = metric(matrix, previous_matrices[(kind, quantity)])
        convergence_quantity = (
            kind == "fabricated" or (kind == "defect" and quantity == "domain")
        )
        check = enforce and convergence_quantity
        limits = {}
        if check and kind == "fabricated":
            limits["MatrixChangePercent"] = args.max_fabricated_matrix_change
            limits["WorstEnergyChangePercent"] = args.max_fabricated_energy_change
        elif check and kind == "defect":
            limits["MatrixChangePercent"] = args.max_domain_defect_change

        failures = []
        for name, limit in limits.items():
            value = values[name]
            if not np.isfinite(value) or value > limit:
                failures.append(f"{name}={value:.6g}% > {limit:.6g}%")
        passed = passed and not failures
        results.append(
            {
                "Kind": kind,
                "Quantity": quantity,
                **values,
                "ConvergenceQuantity": convergence_quantity,
                "Checked": check,
                "Passed": not failures,
                "Failures": failures,
            }
        )
    return passed, results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("calibration", type=Path, nargs="?")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--palace", type=Path)
    parser.add_argument("--orders", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument(
        "--case",
        action="append",
        help=(
            "Ordered mesh-resolution case as NAME=CALIBRATION. Completed source "
            "cases are reused directly."
        ),
    )
    parser.add_argument("--fixed-order", type=int)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--max-fabricated-matrix-change", type=float, default=5.0)
    parser.add_argument("--max-fabricated-energy-change", type=float, default=10.0)
    parser.add_argument("--max-domain-defect-change", type=float, default=5.0)
    parser.add_argument(
        "--gate-all-transitions",
        action="store_true",
        help=(
            "Apply acceptance gates to every adjacent-order comparison; by default "
            "only the highest-order transition is gated"
        ),
    )
    args = parser.parse_args()

    if args.case and args.calibration is not None:
        parser.error("CALIBRATION and --case cannot be used together")
    if not args.case and args.calibration is None:
        parser.error("CALIBRATION is required without --case")
    if args.case and args.fixed_order is None:
        parser.error("--fixed-order is required with --case")
    if not args.case and args.fixed_order is not None:
        parser.error("--fixed-order requires --case")

    orders = sorted(set(args.orders))
    if not args.case and (len(orders) < 2 or any(order < 1 for order in orders)):
        parser.error("--orders requires at least two positive FEM orders")
    if args.fixed_order is not None and args.fixed_order < 1:
        parser.error("--fixed-order must be positive")
    if args.ranks < 1:
        parser.error("--ranks must be positive")

    required = (
        "probe-thin.json",
        "probe-fabricated.json",
        "probe-manifest.json",
        "process-library.json",
    )
    if args.case:
        levels = []
        for specification in args.case:
            if "=" not in specification:
                parser.error(f'--case "{specification}" must use NAME=CALIBRATION')
            name, value = specification.split("=", 1)
            if not name:
                parser.error("--case name cannot be empty")
            levels.append(
                (name, Path(value).expanduser().resolve(), args.fixed_order)
            )
        if len(levels) < 2:
            parser.error("--case requires at least two mesh-resolution levels")
        if len({name for name, _, _ in levels}) != len(levels):
            parser.error("--case names must be unique")
        study = "MeshResolution"
    else:
        calibration = args.calibration.expanduser().resolve()
        levels = [
            (f"p{order}", calibration, order)
            for order in orders
        ]
        study = "FEMOrder"

    for _, calibration, _ in levels:
        missing = [
            name for name in required if not (calibration / name).is_file()
        ]
        if missing:
            raise FileNotFoundError(
                f"{calibration} is missing probe input: {', '.join(missing)}"
            )
    if (
        not args.prepare_only
        and args.palace is None
        and not (args.case and all(complete(level[1]) for level in levels))
    ):
        parser.error(
            "--palace is required unless --prepare-only is used or every "
            "--case is complete"
        )

    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    cases = []
    for name, calibration, order in levels:
        if args.case and complete(calibration):
            actual_order = configured_order(calibration)
            if actual_order != order:
                raise ValueError(
                    f"{calibration} is complete at p={actual_order}, expected p={order}"
                )
            print(f"Reusing completed {name} probe responses in {calibration}")
            cases.append(calibration)
            continue
        case, configs = prepare_case(calibration, output, name, order)
        cases.append(case)
        if args.prepare_only:
            continue
        if complete(case) and not args.force:
            print(f"Reusing completed {name} probe responses")
            continue
        for config in configs:
            command = [args.palace]
            command.extend(["--serial"] if args.ranks == 1 else ["-np", args.ranks])
            command.append(config)
            run(command)
        if not complete(case):
            raise RuntimeError(f"Palace did not complete probe responses in {case}")

    if args.prepare_only:
        print(f"Prepared {len(cases)} {study} cases under {output}")
        return

    loaded = [
        comparison.load_case(f"{name}={case}")
        for (name, _, _), case in zip(levels, cases)
    ]
    reference_metadata = loaded[0]["metadata"]
    for case in loaded[1:]:
        if case["metadata"] != reference_metadata:
            raise ValueError(
                f"{case['root']} uses incompatible probe or coupon metadata"
            )
    comparisons = []
    passed = True
    for index, (previous, current) in enumerate(zip(loaded, loaded[1:])):
        gate_applied = args.gate_all_transitions or index == len(loaded) - 2
        step_passed, results = compare_cases(
            previous, current, args, gate_applied
        )
        passed = passed and step_passed
        comparisons.append(
            {
                "From": previous["name"],
                "To": current["name"],
                "GateApplied": gate_applied,
                "Passed": step_passed,
                "Metrics": results,
            }
        )

    report = {
        "Version": 1,
        "Study": study,
        "Calibration": (
            str(levels[0][1]) if study == "FEMOrder" else None
        ),
        "Orders": orders if study == "FEMOrder" else None,
        "Cases": (
            [
                {
                    "Name": name,
                    "Calibration": str(calibration),
                    "Order": order,
                }
                for name, calibration, order in levels
            ]
            if study == "MeshResolution"
            else None
        ),
        "ThresholdsPercent": {
            "FabricatedMatrixChange": args.max_fabricated_matrix_change,
            "FabricatedWorstEnergyChange": args.max_fabricated_energy_change,
            "DomainDefectMatrixChange": args.max_domain_defect_change,
        },
        "GateAllTransitions": args.gate_all_transitions,
        "Passed": passed,
        "Comparisons": comparisons,
    }
    destination = output / "probe-convergence.json"
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(destination)
    for step in comparisons:
        status = (
            ("PASS" if step["Passed"] else "FAIL")
            if step["GateApplied"]
            else "DIAGNOSTIC"
        )
        print(f"{step['From']} -> {step['To']}: {status}")
        for row in step["Metrics"]:
            if row["ConvergenceQuantity"]:
                print(
                    f"  {row['Kind']}/{row['Quantity']}: "
                    f"matrix={row['MatrixChangePercent']:.3f}%, "
                    f"worst-energy={row['WorstEnergyChangePercent']:.3f}%"
                )
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
