#!/usr/bin/env python3

"""Validate and compact paired thin/fabricated corner response matrices."""

import argparse
import csv
from pathlib import Path

import numpy as np


def clean_header(row):
    return [entry.strip() for entry in row]


def read_domain_matrix(path, size):
    with path.open(newline="") as stream:
        rows = csv.reader(stream)
        header = clean_header(next(rows))
        columns = {name: index for index, name in enumerate(header)}
        matrix = np.zeros((size, size))
        count = 0
        for row in rows:
            i = int(float(row[columns["basis_i"]])) - 1
            j = int(float(row[columns["basis_j"]])) - 1
            value = float(row[columns["Q_ij (J)"]])
            matrix[i, j] = matrix[j, i] = value
            count += 1
    expected = size * (size + 1) // 2
    if count != expected:
        raise ValueError(f"{path} has {count} entries; expected {expected}")
    return matrix


def aggregate_localized_surface_matrix(source, destination, size):
    values = {}
    with source.open(newline="") as stream:
        rows = csv.reader(stream)
        header = clean_header(next(rows))
        columns = {name: index for index, name in enumerate(header)}
        for row in rows:
            key = (
                int(float(row[columns["interface"]])),
                int(float(row[columns["basis_i"]])),
                int(float(row[columns["basis_j"]])),
            )
            values[key] = values.get(key, 0.0) + float(row[columns["Q_ij (J)"]])

    interfaces = sorted({key[0] for key in values})
    expected = size * (size + 1) // 2
    for interface in interfaces:
        count = sum(key[0] == interface for key in values)
        if count != expected:
            raise ValueError(
                f"{source} interface {interface} has {count} entries; "
                f"expected {expected}"
            )

    with destination.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(
            ["interface", "edge", "basis_i", "basis_j", "Q_total_ij (J)"]
        )
        for key in sorted(values):
            interface, i, j = key
            writer.writerow((interface, 1, i, j, f"{values[key]:+.16e}"))

    matrices = {}
    for (interface, i, j), value in values.items():
        matrix = matrices.setdefault(interface, np.zeros((size, size)))
        matrix[i - 1, j - 1] = matrix[j - 1, i - 1] = value
    return matrices


def check_positive(name, matrix, strict):
    eigenvalues = np.linalg.eigvalsh(matrix)
    scale = max(abs(eigenvalues[-1]), 1.0e-300)
    tolerance = 1.0e-9 * scale
    if eigenvalues[0] < -tolerance or (strict and eigenvalues[0] <= tolerance):
        qualifier = "positive definite" if strict else "positive semidefinite"
        raise ValueError(
            f"{name} is not {qualifier}: "
            f"lambda_min={eigenvalues[0]:.6e}, lambda_max={eigenvalues[-1]:.6e}"
        )
    condition = eigenvalues[-1] / eigenvalues[0] if strict else np.inf
    return eigenvalues[0], eigenvalues[-1], condition


def read_single_row(path):
    with path.open(newline="") as stream:
        rows = csv.reader(stream)
        header = clean_header(next(rows))
        row = next(rows)
    return {name: float(row[index]) for index, name in enumerate(header)}


def read_interface_rows(path):
    with path.open(newline="") as stream:
        rows = csv.reader(stream)
        header = clean_header(next(rows))
        result = {}
        for row in rows:
            values = {
                name: float(row[index]) for index, name in enumerate(header)
            }
            result[int(values["interface"])] = values
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("calibration", type=Path)
    args = parser.parse_args()
    root = args.calibration.expanduser().resolve()

    basis = np.loadtxt(root / "basis-points.csv", delimiter=",", skiprows=1)
    size = len(np.atleast_2d(basis))
    results = {}
    for kind in ("thin", "fabricated"):
        postpro = root / "postpro" / kind
        domain_path = postpro / "domain-response-matrix.csv"
        surface_path = postpro / "surface-response-matrix.csv"
        if not domain_path.is_file() or not surface_path.is_file():
            raise FileNotFoundError(f"Incomplete {kind} response output in {postpro}")

        domain = read_domain_matrix(domain_path, size)
        compact_path = postpro / "surface-response-matrix-aggregate.csv"
        surfaces = aggregate_localized_surface_matrix(
            surface_path, compact_path, size
        )
        results[kind] = (domain, surfaces)

        minimum, maximum, condition = check_positive(
            f"{kind} domain response", domain, True
        )
        print(
            f"{kind}: basis={size}, domain eigenvalue range "
            f"[{minimum:.6e}, {maximum:.6e}], condition={condition:.3f}"
        )
        for interface, matrix in sorted(surfaces.items()):
            minimum, maximum, _ = check_positive(
                f"{kind} interface {interface} response", matrix, False
            )
            print(
                f"  interface {interface}: eigenvalue range "
                f"[{minimum:.6e}, {maximum:.6e}]"
            )
        print(f"  compact surface matrix: {compact_path}")

    thin_domain, thin_surfaces = results["thin"]
    fabricated_domain, fabricated_surfaces = results["fabricated"]
    domain_defect = np.linalg.norm(fabricated_domain - thin_domain) / np.linalg.norm(
        thin_domain
    )
    print(f"relative domain defect norm: {domain_defect:.6e}")
    for interface in sorted(thin_surfaces):
        defect = np.linalg.norm(
            fabricated_surfaces[interface] - thin_surfaces[interface]
        ) / np.linalg.norm(thin_surfaces[interface])
        print(f"relative interface {interface} defect norm: {defect:.6e}")

    coefficients_path = root / "heldout-coefficients.csv"
    if coefficients_path.is_file():
        coefficients = np.loadtxt(
            coefficients_path, delimiter=",", skiprows=1
        )
        for kind in ("thin", "fabricated"):
            postpro = root / "postpro" / f"heldout-{kind}"
            domain_path = postpro / "domain-E.csv"
            surface_path = postpro / "surface-Q.csv"
            if not domain_path.is_file() or not surface_path.is_file():
                print(f"held-out {kind}: not run")
                continue
            direct_domain = read_single_row(domain_path)["E_elec (J)"]
            surface_row = read_single_row(surface_path)
            edge_rows = read_interface_rows(postpro / "surface-Q-edge.csv")
            predicted_domain = coefficients @ results[kind][0] @ coefficients
            domain_error = predicted_domain / direct_domain - 1.0
            print(f"held-out {kind} domain error: {domain_error:+.6e}")
            for interface, matrix in sorted(results[kind][1].items()):
                direct_surface = surface_row[f"p_surf[{interface}]"] * direct_domain
                direct_surface -= edge_rows[interface]["E_out (J)"]
                predicted_surface = coefficients @ matrix @ coefficients
                error = predicted_surface / direct_surface - 1.0
                print(
                    f"  interface {interface} held-out error: {error:+.6e}"
                )


if __name__ == "__main__":
    main()
