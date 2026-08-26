#!/usr/bin/env python3

"""Qualify a straight-edge response coupon with convergence and held-out tests."""

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


def clean(row):
    return {key.strip(): value for key, value in row.items()}


def read_matrix(path, size, value_column):
    matrix = np.zeros((size, size))
    seen = set()
    with path.open(newline="") as stream:
        for input_row in csv.DictReader(stream, skipinitialspace=True):
            row = clean(input_row)
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            if not 0 <= i <= j < size or (i, j) in seen:
                raise ValueError(f"{path} contains an invalid or duplicate basis pair")
            seen.add((i, j))
            matrix[i, j] = matrix[j, i] = float(row[value_column])
    if len(seen) != size * (size + 1) // 2:
        raise ValueError(f"{path} does not contain a complete response matrix")
    return matrix


def read_surface_matrices(path, size):
    entries = {}
    with path.open(newline="") as stream:
        for input_row in csv.DictReader(stream, skipinitialspace=True):
            row = clean(input_row)
            interface = int(float(row["interface"]))
            edge = int(float(row["edge"]))
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            if not 0 <= i <= j < size:
                raise ValueError(f"{path} contains an invalid basis pair")
            key = (interface, edge, i, j)
            value = float(row["Q_total_ij (J)"])
            old = entries.setdefault(key, value)
            scale = max(abs(old), abs(value), np.finfo(float).tiny)
            if abs(old - value) > 1.0e-10 * scale:
                raise ValueError(f"{path} contains inconsistent repeated total energies")

    matrices = {}
    edge_pairs = {}
    for (interface, edge, i, j), value in entries.items():
        matrix = matrices.setdefault(interface, np.zeros((size, size)))
        matrix[i, j] += value
        if i != j:
            matrix[j, i] += value
        edge_pairs.setdefault((interface, edge), set()).add((i, j))
    expected = {(i, j) for i in range(size) for j in range(i, size)}
    if any(pairs != expected for pairs in edge_pairs.values()):
        raise ValueError(f"{path} contains an incomplete per-edge response matrix")
    return matrices


def read_one_row(path):
    with path.open(newline="") as stream:
        row = next(csv.DictReader(stream, skipinitialspace=True))
    return {key.strip(): float(value) for key, value in row.items()}


def library_model(root):
    library = json.loads((root / "process-library.json").read_text())
    if len(library.get("Models", [])) != 1:
        raise ValueError(f"{root} must contain exactly one coupon model")
    return library, library["Models"][0]


def resolve(root, value):
    path = Path(value)
    return path if path.is_absolute() else root / path


def active_indices(model, size):
    constrained = [int(index) - 1 for index in model.get("ZeroTraceIndices", [])]
    if any(index < 0 or index >= size for index in constrained):
        raise ValueError("ZeroTraceIndices contains an invalid basis index")
    if len(set(constrained)) != len(constrained):
        raise ValueError("ZeroTraceIndices contains duplicate basis indices")
    active = np.delete(np.arange(size), constrained)
    if not len(active):
        raise ValueError("Response model has no free trace basis functions")
    return active


def restrict(matrix, active):
    return matrix[np.ix_(active, active)]


def load_response(root):
    library, model = library_model(root)
    coefficients = np.atleast_1d(
        np.loadtxt(root / "heldout_coefficients.csv", delimiter=",", skiprows=1)
    )
    size = len(coefficients)
    responses = {}
    for kind in ("Thin", "Fabricated"):
        domain_path = resolve(root, model[f"{kind}Matrix"])
        surface_path = resolve(root, model[f"{kind}SurfaceMatrix"])
        responses[kind.lower()] = {
            "domain": read_matrix(domain_path, size, "Q_ij (J)"),
            "surfaces": read_surface_matrices(surface_path, size),
            "postpro": domain_path.parent,
        }
    return {
        "root": root,
        "library": library,
        "model": model,
        "coefficients": coefficients,
        "active": active_indices(model, size),
        "responses": responses,
    }


def minimum_eigenvalue_check(matrix, strict):
    eigenvalues = np.linalg.eigvalsh(matrix)
    scale = max(abs(eigenvalues[-1]), np.finfo(float).tiny)
    tolerance = 1.0e-9 * scale
    passed = eigenvalues[0] > tolerance if strict else eigenvalues[0] >= -tolerance
    return {
        "MinimumEigenvalue": eigenvalues[0],
        "MaximumEigenvalue": eigenvalues[-1],
        "Passed": bool(passed),
    }


def worst_energy_change(current, previous):
    eigenvalues, eigenvectors = np.linalg.eigh(current)
    scale = max(np.max(np.abs(eigenvalues)), np.finfo(float).tiny)
    active = eigenvalues > 1.0e-10 * scale
    if not np.any(active):
        return float("nan")
    basis = eigenvectors[:, active] / np.sqrt(eigenvalues[active])
    relative = basis.T @ previous @ basis
    return np.max(np.abs(np.linalg.eigvalsh(relative) - 1.0))


def convergence_check(name, current, previous, matrix_limit, energy_limit=None):
    matrix_change = np.linalg.norm(current - previous) / max(
        np.linalg.norm(current), np.finfo(float).tiny
    )
    energy_change = worst_energy_change(current, previous)
    passed = math.isfinite(matrix_change) and 100.0 * matrix_change <= matrix_limit
    if energy_limit is not None:
        passed = (
            passed
            and math.isfinite(energy_change)
            and 100.0 * energy_change <= energy_limit
        )
    return {
        "Quantity": name,
        "MatrixChangePercent": 100.0 * matrix_change,
        "WorstEnergyChangePercent": 100.0 * energy_change,
        "MatrixLimitPercent": matrix_limit,
        "EnergyLimitPercent": energy_limit,
        "Passed": bool(passed),
    }


def heldout_checks(response, maximum_error):
    coefficients = response["coefficients"]
    checks = []
    for kind, data in response["responses"].items():
        postpro = data["postpro"].parent / f"heldout_{data['postpro'].name}"
        domain_path = postpro / "domain-E.csv"
        surface_path = postpro / "surface-Q.csv"
        if not domain_path.is_file() or not surface_path.is_file():
            raise FileNotFoundError(f"Incomplete held-out solve in {postpro}")
        direct_domain = read_one_row(domain_path)["E_elec (J)"]
        predicted_domain = coefficients @ data["domain"] @ coefficients
        quantities = [("domain", direct_domain, predicted_domain)]
        surface_row = read_one_row(surface_path)
        for interface, matrix in sorted(data["surfaces"].items()):
            direct = surface_row.get(f"p_surf[{interface}]", 0.0) * direct_domain
            predicted = coefficients @ matrix @ coefficients
            quantities.append((f"interface-{interface}", direct, predicted))
        for quantity, direct, predicted in quantities:
            scale = max(abs(direct), abs(predicted), np.finfo(float).tiny)
            if max(abs(direct), abs(predicted)) <= 1.0e-14 * abs(direct_domain):
                error = 0.0
            elif abs(direct) <= 1.0e-14 * scale:
                error = math.copysign(math.inf, predicted)
            else:
                error = 100.0 * (predicted / direct - 1.0)
            checks.append(
                {
                    "Kind": kind,
                    "Quantity": quantity,
                    "DirectEnergy": direct,
                    "PredictedEnergy": predicted,
                    "RelativeErrorPercent": error,
                    "Passed": bool(
                        math.isfinite(error) and abs(error) <= maximum_error
                    ),
                }
            )
    return checks


def same_model(first, second):
    keys = (
        "Topology",
        "Separation",
        "Edges",
        "BoundaryCondition",
        "Interfaces",
        "ContourGroups",
        "ZeroTraceIndices",
    )
    return all(first.get(key) == second.get(key) for key in keys)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("calibration", type=Path)
    parser.add_argument(
        "--previous",
        type=Path,
        help=(
            "Optional lower-resolution calibration. Without it, qualification "
            "checks matrix definiteness and held-out accuracy only."
        ),
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-heldout-error", type=float, default=10.0)
    parser.add_argument("--max-fabricated-matrix-change", type=float, default=5.0)
    parser.add_argument("--max-fabricated-energy-change", type=float, default=10.0)
    parser.add_argument("--max-domain-defect-change", type=float, default=5.0)
    args = parser.parse_args()
    positive = (
        args.max_heldout_error,
        args.max_fabricated_matrix_change,
        args.max_fabricated_energy_change,
        args.max_domain_defect_change,
    )
    if any(value <= 0.0 for value in positive):
        parser.error("qualification thresholds must be positive")

    root = args.calibration.expanduser().resolve()
    current = load_response(root)
    previous_root = (
        args.previous.expanduser().resolve() if args.previous else None
    )
    previous = load_response(previous_root) if previous_root else None
    if previous:
        if not same_model(current["model"], previous["model"]):
            raise ValueError(
                "Current and previous calibrations describe different coupons"
            )
        if len(current["coefficients"]) != len(previous["coefficients"]):
            raise ValueError(
                "Current and previous coupons use different trace spaces"
            )

    matrix_checks = []
    active = current["active"]
    for kind, response in current["responses"].items():
        matrix_checks.append(
            {
                "Kind": kind,
                "Quantity": "domain",
                **minimum_eigenvalue_check(
                    restrict(response["domain"], active), kind == "fabricated"
                ),
            }
        )
        for interface, matrix in sorted(response["surfaces"].items()):
            matrix_checks.append(
                {
                    "Kind": kind,
                    "Quantity": f"interface-{interface}",
                    **minimum_eigenvalue_check(restrict(matrix, active), False),
                }
            )

    convergence = []
    if previous:
        fabricated = current["responses"]["fabricated"]
        old_fabricated = previous["responses"]["fabricated"]
        convergence.append(
            convergence_check(
                "fabricated-domain",
                restrict(fabricated["domain"], active),
                restrict(old_fabricated["domain"], active),
                args.max_fabricated_matrix_change,
                args.max_fabricated_energy_change,
            )
        )
        for interface, matrix in sorted(fabricated["surfaces"].items()):
            convergence.append(
                convergence_check(
                    f"fabricated-interface-{interface}",
                    restrict(matrix, active),
                    restrict(old_fabricated["surfaces"][interface], active),
                    args.max_fabricated_matrix_change,
                )
            )
        defect = (
            current["responses"]["fabricated"]["domain"]
            - current["responses"]["thin"]["domain"]
        )
        old_defect = (
            previous["responses"]["fabricated"]["domain"]
            - previous["responses"]["thin"]["domain"]
        )
        convergence.append(
            convergence_check(
                "domain-defect",
                restrict(defect, active),
                restrict(old_defect, active),
                args.max_domain_defect_change,
            )
        )
    heldout = heldout_checks(current, args.max_heldout_error)
    passed = all(
        check["Passed"] for check in matrix_checks + convergence + heldout
    )
    report = {
        "Version": 1,
        "Calibration": str(root),
        "PreviousCalibration": (
            str(previous_root) if previous_root is not None else None
        ),
        "Model": {
            "Name": current["model"].get("Name"),
            "Topology": current["model"].get("Topology"),
            "Separation": current["model"].get("Separation"),
        },
        "MatrixChecks": matrix_checks,
        "ConvergenceChecks": convergence,
        "HeldoutChecks": heldout,
        "Passed": passed,
    }
    destination = (
        args.output.expanduser().resolve()
        if args.output
        else root / "qualification.json"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(destination)
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
