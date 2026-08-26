#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Replay an exported device patch trace directly on thin/fabricated coupons."""

import argparse
import csv
import json
import shlex
import subprocess
from pathlib import Path

import numpy as np


def read_rows(path):
    with path.open(newline="") as stream:
        return [
            {key.strip(): float(value) for key, value in row.items()}
            for row in csv.DictReader(stream, skipinitialspace=True)
        ]


def read_row(path):
    return read_rows(path)[0]


def read_domain_matrix(path, size):
    matrix = np.zeros((size, size))
    for row in read_rows(path):
        i = int(row["basis_i"]) - 1
        j = int(row["basis_j"]) - 1
        matrix[i, j] = matrix[j, i] = row["Q_ij (J)"]
    return matrix


def read_surface_matrices(path, size):
    unique = {}
    for row in read_rows(path):
        key = (
            int(row["interface"]),
            int(row["edge"]),
            int(row["basis_i"]) - 1,
            int(row["basis_j"]) - 1,
        )
        unique.setdefault(key, row["Q_total_ij (J)"])
    matrices = {}
    for (interface, _edge, i, j), value in unique.items():
        matrix = matrices.setdefault(interface, np.zeros((size, size)))
        matrix[i, j] += value
        if i != j:
            matrix[j, i] += value
    return matrices


def device_trace(path, source, patch):
    coefficients = {}
    for row in read_rows(path):
        if int(row["i"]) == source and int(row["patch"]) == patch:
            coefficients[int(row["coefficient"])] = row["value (V)"]
    if not coefficients:
        raise ValueError(f"No exported trace for source {source} patch {patch}")
    return np.asarray(
        [coefficients[index] for index in sorted(coefficients)]
    )


def write_trace(source, destination, values):
    rows = read_rows(source)
    if len(rows) % 3 != 0:
        raise ValueError(f"{source} is not a triangle trace file")
    with destination.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["x", "y", "z", "V", "triangle"])
        for row in rows:
            writer.writerow(
                [
                    f"{row['x']:.16e}",
                    f"{row['y']:.16e}",
                    f"{row['z']:.16e}",
                    f"{values[(row['x'], row['y'], row['z'])]:.16e}",
                    int(row["triangle"]),
                ]
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coupon", type=Path, help="Generated spatial coupon directory")
    parser.add_argument("--device-traces", type=Path, required=True)
    parser.add_argument("--source", type=int, required=True)
    parser.add_argument("--patch", type=int, required=True)
    parser.add_argument("--palace", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    coupon = args.coupon.expanduser().resolve()
    output = (
        args.output.expanduser().resolve()
        if args.output
        else coupon / f"device-trace-replay-s{args.source}-p{args.patch}"
    )
    output.mkdir(parents=True, exist_ok=True)

    coefficients = device_trace(
        args.device_traces.expanduser().resolve(), args.source, args.patch
    )
    basis_points = np.loadtxt(coupon / "basis-points.csv", delimiter=",", skiprows=1)
    library = json.loads((coupon / "process-library.json").read_text())
    model = library["Models"][0]
    contour_size = len(np.atleast_2d(basis_points))
    conductor_count = len(model.get("ConductorReferences", [None]))
    if len(coefficients) != contour_size + conductor_count - 1:
        raise ValueError(
            f"Exported trace has {len(coefficients)} coefficients; coupon expects "
            f"{contour_size} contour knots plus conductor states"
        )

    # Build the imposed matching-surface potential from the coupon's own zero trace
    # geometry and the exported nodal coefficients. Trace knots are matched by exact
    # canonical coordinates.
    zero_rows = read_rows(coupon / "zero-trace.csv")
    values = {(row["x"], row["y"], row["z"]): 0.0 for row in zero_rows}
    canonical = np.atleast_2d(basis_points)
    frame_points = {}
    for index, point in enumerate(canonical):
        frame_points[tuple(round(value, 12) for value in point)] = coefficients[index]
    # Coupon trace files store local (already canonical) coordinates.
    matched = 0
    for key in list(values):
        rounded = tuple(round(value, 12) for value in key)
        if rounded in frame_points:
            values[key] = frame_points[rounded]
            matched += 1
    if matched != len(frame_points):
        raise ValueError(
            f"Matched {matched} of {len(frame_points)} exported trace knots to coupon "
            "surface nodes; coordinate frames disagree"
        )

    trace_path = output / "device-trace.csv"
    write_trace(coupon / "zero-trace.csv", trace_path, values)

    # The direct replay imposes the exported contour trace with every conductor
    # grounded, so the matrix prediction must use zero conductor-state coefficients.
    # Conductor-state and cross contributions are validated separately by the
    # trace/conductor superposition audit; by linearity the two suffice together.
    imposed = coefficients.copy()
    imposed[contour_size:] = 0.0

    report = {"Version": 1, "Coupon": str(coupon), "Source": args.source,
              "Patch": args.patch, "Results": []}
    for kind in ("thin", "fabricated"):
        config = json.loads((coupon / f"heldout_spatial_{kind}.json").read_text())
        config["Problem"]["Output"] = str(output / kind)
        config["Problem"]["Verbose"] = 0
        potential = config["Boundaries"]["PrescribedPotential"][0]
        potential["DataFile"] = str(trace_path)
        # Ground every conductor so the direct excitation is the pure contour trace.
        potential.pop("TerminalAttributes", None)
        config_path = output / f"{kind}.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n")
        command = [
            str(args.palace),
            *( ["--serial"] if args.ranks == 1 else ["-np", str(args.ranks)]),
            str(config_path),
        ]
        print("+ " + shlex.join(command), flush=True)
        subprocess.run(command, cwd=output, check=True)

        direct_domain = read_row(output / kind / "domain-E.csv")["E_elec (J)"]
        surface = read_row(output / kind / "surface-Q.csv")
        matrix = read_domain_matrix(
            coupon / "postpro" / f"spatial_{kind}" / "domain-response-matrix.csv",
            len(coefficients),
        )
        surfaces = read_surface_matrices(
            coupon / "postpro" / f"spatial_{kind}" / "surface-response-matrix.csv",
            len(coefficients),
        )
        entry = {
            "Kind": kind,
            "DirectDomainEnergy": direct_domain,
            "PredictedDomainEnergy": float(imposed @ matrix @ imposed),
            "Interfaces": {},
        }
        for interface, response in sorted(surfaces.items()):
            entry["Interfaces"][str(interface)] = {
                "DirectEnergy": surface.get(f"p_surf[{interface}]", 0.0) * direct_domain,
                "PredictedEnergy": float(imposed @ response @ imposed),
            }
        report["Results"].append(entry)

    conductor_states = coefficients[contour_size:]
    report["ConductorStates"] = conductor_states.tolist()
    report_path = output / "device-trace-replay.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(report_path)


if __name__ == "__main__":
    main()
