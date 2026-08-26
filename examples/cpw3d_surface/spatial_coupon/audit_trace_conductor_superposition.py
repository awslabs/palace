#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Audit exact trace/conductor superposition for a generated spatial coupon."""

import argparse
import csv
import json
import math
import shlex
import subprocess
from pathlib import Path

import numpy as np


def load(path):
    return json.loads(path.read_text())


def clean(row):
    return {key.strip(): value for key, value in row.items()}


def read_row(path):
    with path.open(newline="") as stream:
        return {
            key.strip(): float(value)
            for key, value in next(csv.DictReader(stream, skipinitialspace=True)).items()
        }


def read_domain_matrix(path, size):
    matrix = np.zeros((size, size))
    seen = set()
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream, skipinitialspace=True):
            row = clean(raw)
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            if not 0 <= i <= j < size or (i, j) in seen:
                raise ValueError(f"{path} contains an invalid basis pair")
            seen.add((i, j))
            matrix[i, j] = matrix[j, i] = float(row["Q_ij (J)"])
    if len(seen) != size * (size + 1) // 2:
        raise ValueError(f"{path} does not contain a complete {size} by {size} matrix")
    return matrix


def read_surface_matrices(path, size):
    unique = {}
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream, skipinitialspace=True):
            row = clean(raw)
            interface = int(float(row["interface"]))
            edge = int(float(row["edge"]))
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            key = (interface, edge, i, j)
            value = float(row["Q_total_ij (J)"])
            previous = unique.setdefault(key, value)
            scale = max(abs(previous), abs(value), np.finfo(float).tiny)
            if abs(previous - value) > 1.0e-10 * scale:
                raise ValueError(f"{path} contains inconsistent repeated energies")

    matrices = {}
    pairs = {}
    for (interface, edge, i, j), value in unique.items():
        if not 0 <= i <= j < size:
            raise ValueError(f"{path} contains an invalid basis pair")
        matrix = matrices.setdefault(interface, np.zeros((size, size)))
        matrix[i, j] += value
        if i != j:
            matrix[j, i] += value
        pairs.setdefault((interface, edge), set()).add((i, j))
    expected = {(i, j) for i in range(size) for j in range(i, size)}
    if any(found != expected for found in pairs.values()):
        raise ValueError(f"{path} contains an incomplete per-edge matrix")
    return matrices


def write_negative_trace(source, destination):
    with source.open(newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        fields = reader.fieldnames
    if not fields or "V" not in fields:
        raise ValueError(f"{source} has no V column")
    with destination.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            row["V"] = f"{-float(row['V']):.16e}"
            writer.writerow(row)


def palace_command(palace, ranks, config):
    return [
        str(palace),
        *( ["--serial"] if ranks == 1 else ["-np", str(ranks)]),
        str(config),
    ]


def run(command, cwd):
    print("+ " + shlex.join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def relative_error(predicted, direct, domain_scale):
    scale = max(abs(predicted), abs(direct), np.finfo(float).tiny)
    if scale <= 1.0e-14 * abs(domain_scale):
        return 0.0
    return abs(predicted - direct) / scale


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coupon", type=Path)
    parser.add_argument("--kind", choices=("thin", "fabricated"), default="fabricated")
    parser.add_argument("--palace", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--max-relative-error", type=float, default=1.0e-6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.ranks <= 0 or args.max_relative_error <= 0.0:
        parser.error("ranks and maximum relative error must be positive")
    if args.execute and args.palace is None:
        parser.error("--execute requires --palace")

    coupon = args.coupon.expanduser().resolve()
    output = (
        args.output.expanduser().resolve()
        if args.output
        else coupon / f"trace-conductor-audit-{args.kind}"
    )
    output.mkdir(parents=True, exist_ok=True)

    production = load(coupon / f"spatial_{args.kind}.json")
    heldout = load(coupon / f"heldout_spatial_{args.kind}.json")
    production_potentials = production["Boundaries"]["PrescribedPotential"]
    terminal_potentials = [
        entry for entry in production_potentials if entry.get("TerminalAttributes")
    ]
    if not terminal_potentials:
        raise ValueError("coupon has no conductor-state response basis")

    heldout_trace = dict(heldout["Boundaries"]["PrescribedPotential"][0])
    heldout_trace.pop("TerminalAttributes", None)
    heldout_trace["Index"] = 1
    matrix_potentials = [heldout_trace]
    for index, source in enumerate(terminal_potentials, start=2):
        source = dict(source)
        source["Index"] = index
        matrix_potentials.append(source)

    matrix = json.loads(json.dumps(production))
    matrix["Problem"]["Output"] = str(output / "matrix")
    matrix["Problem"]["Verbose"] = 0
    matrix["Boundaries"]["PrescribedPotential"] = matrix_potentials
    matrix["Solver"]["Electrostatic"]["ResponseMatrix"] = True
    matrix["Solver"]["Electrostatic"]["AggregateResponseMatrix"] = True
    matrix_path = output / "matrix.json"
    matrix_path.write_text(json.dumps(matrix, indent=2) + "\n")

    combined = json.loads(json.dumps(heldout))
    combined["Problem"]["Output"] = str(output / "combined")
    combined["Problem"]["Verbose"] = 0
    combined_path = output / "combined.json"
    combined_path.write_text(json.dumps(combined, indent=2) + "\n")

    negative_trace = output / "negative-heldout-trace.csv"
    write_negative_trace(coupon / "heldout-trace.csv", negative_trace)
    difference = json.loads(json.dumps(heldout))
    difference["Problem"]["Output"] = str(output / "difference")
    difference["Problem"]["Verbose"] = 0
    difference["Boundaries"]["PrescribedPotential"][0]["DataFile"] = str(
        negative_trace
    )
    difference_path = output / "difference.json"
    difference_path.write_text(json.dumps(difference, indent=2) + "\n")

    if not args.execute:
        print(matrix_path)
        print(combined_path)
        print(difference_path)
        return

    expected = (
        output / "matrix" / "domain-response-matrix.csv",
        output / "combined" / "domain-E.csv",
        output / "difference" / "domain-E.csv",
    )
    if args.force or not all(path.is_file() for path in expected):
        for path in (matrix_path, combined_path, difference_path):
            run(palace_command(args.palace, args.ranks, path), output)

    size = 1 + len(terminal_potentials)
    plus = np.ones(size)
    minus = np.ones(size)
    minus[0] = -1.0
    matrix_root = output / "matrix"
    domain_matrix = read_domain_matrix(
        matrix_root / "domain-response-matrix.csv", size
    )
    surface_matrices = read_surface_matrices(
        matrix_root / "surface-response-matrix.csv", size
    )

    results = []
    worst = 0.0
    for name, coefficients in (("combined", plus), ("difference", minus)):
        root = output / name
        domain_row = read_row(root / "domain-E.csv")
        direct_domain = domain_row["E_elec (J)"]
        predicted_domain = float(coefficients @ domain_matrix @ coefficients)
        error = relative_error(predicted_domain, direct_domain, direct_domain)
        worst = max(worst, error)
        results.append(
            {
                "Excitation": name,
                "Quantity": "domain",
                "DirectEnergy": direct_domain,
                "PredictedEnergy": predicted_domain,
                "RelativeError": error,
            }
        )
        surface_row = read_row(root / "surface-Q.csv")
        for interface, response in sorted(surface_matrices.items()):
            direct = surface_row.get(f"p_surf[{interface}]", 0.0) * direct_domain
            predicted = float(coefficients @ response @ coefficients)
            error = relative_error(predicted, direct, direct_domain)
            worst = max(worst, error)
            results.append(
                {
                    "Excitation": name,
                    "Quantity": f"interface-{interface}",
                    "DirectEnergy": direct,
                    "PredictedEnergy": predicted,
                    "RelativeError": error,
                }
            )

    report = {
        "Version": 1,
        "Coupon": str(coupon),
        "Kind": args.kind,
        "BasisSize": size,
        "MaximumRelativeError": worst,
        "Limit": args.max_relative_error,
        "Passed": bool(math.isfinite(worst) and worst <= args.max_relative_error),
        "Results": results,
    }
    report_path = output / "trace-conductor-superposition.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(report_path)
    if not report["Passed"]:
        raise SystemExit(
            f"trace/conductor superposition error {worst:.3e} exceeds "
            f"{args.max_relative_error:.3e}"
        )


if __name__ == "__main__":
    main()
