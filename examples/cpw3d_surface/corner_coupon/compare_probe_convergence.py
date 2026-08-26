#!/usr/bin/env python3

"""Compare small projected corner-coupon response matrices across FEM orders."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def clean_row(row):
    return {key.strip(): value for key, value in row.items()}


def read_matrix(path, size, value_column):
    matrix = np.zeros((size, size))
    seen = set()
    with path.open(newline="") as stream:
        for input_row in csv.DictReader(stream, skipinitialspace=True):
            row = clean_row(input_row)
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            if not 0 <= i < size or not 0 <= j < size or i > j:
                raise ValueError(f"{path} contains invalid basis pair ({i + 1}, {j + 1})")
            key = (i, j)
            if key in seen:
                raise ValueError(f"{path} contains duplicate basis pair ({i + 1}, {j + 1})")
            seen.add(key)
            matrix[i, j] = matrix[j, i] = float(row[value_column])
    expected = size * (size + 1) // 2
    if len(seen) != expected:
        raise ValueError(f"{path} has {len(seen)} entries; expected {expected}")
    return matrix


def read_surface_matrices(path, size):
    matrices = {}
    seen = {}
    with path.open(newline="") as stream:
        for input_row in csv.DictReader(stream, skipinitialspace=True):
            row = clean_row(input_row)
            interface = int(float(row["interface"]))
            edge = int(float(row["edge"]))
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            if not 0 <= i < size or not 0 <= j < size or i > j:
                raise ValueError(f"{path} contains invalid basis pair ({i + 1}, {j + 1})")
            key = (edge, i, j)
            interface_seen = seen.setdefault(interface, set())
            if key in interface_seen:
                raise ValueError(
                    f"{path} contains duplicate interface/edge/basis entry "
                    f"({interface}, {edge}, {i + 1}, {j + 1})"
                )
            interface_seen.add(key)
            matrix = matrices.setdefault(interface, np.zeros((size, size)))
            value = float(row["Q_total_ij (J)"])
            matrix[i, j] += value
            if i != j:
                matrix[j, i] += value

    expected_pairs = {(i, j) for i in range(size) for j in range(i, size)}
    for interface, entries in seen.items():
        pairs = {(i, j) for _, i, j in entries}
        if pairs != expected_pairs:
            raise ValueError(
                f"{path} interface {interface} does not contain every basis pair"
            )
    return matrices


def response_metadata(root):
    manifest = json.loads((root / "probe-manifest.json").read_text())
    library = json.loads((root / "process-library.json").read_text())
    model = library["Models"][0]
    metadata = {
        "MatchingRadius": library["MatchingRadius"],
        "Fabrication": library.get("Fabrication"),
        "Topology": model["Topology"],
        "Angle": model.get("Angle"),
        "CornerRadius": model.get("CornerRadius", 0.0),
        "Probes": manifest["Probes"],
    }
    interface_names = {
        int(interface["Coupon"]): interface["Type"]
        for interface in model.get("Interfaces", [])
    }
    return (
        metadata,
        [probe["Name"] for probe in manifest["Probes"]],
        interface_names,
    )


def load_case(specification):
    if "=" not in specification:
        raise ValueError(f'Case "{specification}" must use NAME=CALIBRATION')
    name, path = specification.split("=", 1)
    if not name:
        raise ValueError("Case name cannot be empty")
    root = Path(path).expanduser().resolve()
    metadata, probes, interface_names = response_metadata(root)
    size = len(probes)
    if size == 0:
        raise ValueError(f"{root} defines no convergence probes")

    responses = {}
    for kind in ("thin", "fabricated"):
        postpro = root / "postpro" / f"probe-{kind}"
        domain_path = postpro / "domain-response-matrix.csv"
        surface_path = postpro / "surface-response-matrix.csv"
        if not domain_path.is_file() or not surface_path.is_file():
            raise FileNotFoundError(f"Incomplete probe response output in {postpro}")
        responses[kind] = {
            "domain": read_matrix(domain_path, size, "Q_ij (J)"),
            "surfaces": read_surface_matrices(surface_path, size),
        }
        # Fabrication can remove one requested local interface entirely. Treat an absent
        # coupon interface as an exact zero matrix, matching coverage-library padding.
        for interface in interface_names:
            responses[kind]["surfaces"].setdefault(
                interface, np.zeros((size, size))
            )
        if not responses[kind]["surfaces"].keys() <= interface_names.keys():
            raise ValueError(
                f"{surface_path} contains an unconfigured coupon interface"
            )
    return {
        "name": name,
        "root": root,
        "metadata": metadata,
        "interface_names": interface_names,
        "responses": responses,
    }


def worst_energy_change(current, previous):
    """Maximum relative quadratic-form change, when current is positive definite."""
    eigenvalues, eigenvectors = np.linalg.eigh(current)
    scale = max(np.max(np.abs(eigenvalues)), np.finfo(float).tiny)
    if eigenvalues[0] <= 1.0e-10 * scale:
        return float("nan")
    inverse_sqrt = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    relative_previous = inverse_sqrt @ previous @ inverse_sqrt
    return np.max(np.abs(np.linalg.eigvalsh(relative_previous) - 1.0))


def matrix_rows(case):
    thin = case["responses"]["thin"]
    fabricated = case["responses"]["fabricated"]
    interface_names = case["interface_names"]
    quantities = {"domain": (thin["domain"], fabricated["domain"])}
    interfaces = thin["surfaces"].keys() | fabricated["surfaces"].keys()
    for interface in sorted(interfaces):
        if interface not in thin["surfaces"] or interface not in fabricated["surfaces"]:
            raise ValueError(
                f"{case['root']} has inconsistent thin/fabricated interface matrices"
            )
        quantities[interface_names.get(interface, f"interface-{interface}")] = (
            thin["surfaces"][interface],
            fabricated["surfaces"][interface],
        )
    for quantity, (thin_matrix, fabricated_matrix) in quantities.items():
        yield "thin", quantity, thin_matrix
        yield "fabricated", quantity, fabricated_matrix
        if quantity == "domain":
            yield "defect", quantity, fabricated_matrix - thin_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        action="append",
        required=True,
        help="Ordered coupon level as NAME=CALIBRATION",
    )
    args = parser.parse_args()
    cases = [load_case(specification) for specification in args.case]
    reference_metadata = cases[0]["metadata"]
    for case in cases[1:]:
        if case["metadata"] != reference_metadata:
            raise ValueError(f"{case['root']} uses incompatible probe or coupon metadata")

    print(
        "case,kind,quantity,matrix_norm_J,matrix_change_pct,"
        "worst_energy_change_pct,min_eigenvalue_J,max_eigenvalue_J"
    )
    previous = {}
    for case in cases:
        for kind, quantity, matrix in matrix_rows(case):
            key = (kind, quantity)
            old = previous.get(key)
            norm = np.linalg.norm(matrix)
            change = (
                float("nan")
                if old is None
                else np.linalg.norm(matrix - old) / max(norm, np.finfo(float).tiny)
            )
            energy_change = (
                float("nan")
                if old is None
                else worst_energy_change(matrix, old)
            )
            eigenvalues = np.linalg.eigvalsh(matrix)
            print(
                f"{case['name']},{kind},{quantity},{norm:.12e},"
                f"{100.0 * change:+.6f},{100.0 * energy_change:+.6f},"
                f"{eigenvalues[0]:+.12e},{eigenvalues[-1]:+.12e}"
            )
            previous[key] = matrix


if __name__ == "__main__":
    main()
