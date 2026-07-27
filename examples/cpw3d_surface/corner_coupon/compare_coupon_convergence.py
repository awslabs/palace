#!/usr/bin/env python3

"""Compare response matrices and held-out energies across corner coupon levels."""

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def clean_row(row):
    return {key.strip(): value for key, value in row.items()}


def read_row(path):
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, skipinitialspace=True)
        try:
            return {
                key: float(value) for key, value in clean_row(next(rows)).items()
            }
        except StopIteration as exc:
            raise ValueError(f"{path} has no data rows") from exc


def read_interface_rows(path):
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, skipinitialspace=True)
        result = {}
        for input_row in rows:
            row = clean_row(input_row)
            result[int(float(row["interface"]))] = {
                key: float(value)
                for key, value in row.items()
            }
        return result


def read_domain_matrix(path, size):
    matrix = np.zeros((size, size))
    count = 0
    with path.open(newline="") as stream:
        for input_row in csv.DictReader(stream, skipinitialspace=True):
            row = clean_row(input_row)
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            matrix[i, j] = matrix[j, i] = float(row["Q_ij (J)"])
            count += 1
    expected = size * (size + 1) // 2
    if count != expected:
        raise ValueError(f"{path} has {count} entries; expected {expected}")
    return matrix


def read_surface_matrices(path, size):
    matrices = {}
    counts = {}
    with path.open(newline="") as stream:
        for input_row in csv.DictReader(stream, skipinitialspace=True):
            row = clean_row(input_row)
            interface = int(float(row["interface"]))
            i = int(float(row["basis_i"])) - 1
            j = int(float(row["basis_j"])) - 1
            value = float(row["Q_total_ij (J)"])
            matrix = matrices.setdefault(interface, np.zeros((size, size)))
            matrix[i, j] += value
            if i != j:
                matrix[j, i] += value
            counts[interface] = counts.get(interface, 0) + 1
    expected = size * (size + 1) // 2
    for interface, count in counts.items():
        if count != expected:
            raise ValueError(
                f"{path} interface {interface} has {count} entries; "
                f"expected {expected}"
            )
    return matrices


def parse_case(specification):
    if "=" not in specification:
        raise ValueError(f'Case "{specification}" must use NAME=CALIBRATION')
    name, path = specification.split("=", 1)
    if not name:
        raise ValueError("Case name cannot be empty")
    root = Path(path).expanduser().resolve()
    required = (
        root / "basis-points.csv",
        root / "heldout-coefficients.csv",
        root / "process-library.json",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing corner coupon convergence input: " + ", ".join(missing)
        )
    return name, root


def compatible_metadata(library):
    model = library["Models"][0]
    return {
        "MatchingRadius": library["MatchingRadius"],
        "Fabrication": library.get("Fabrication"),
        "Topology": model["Topology"],
        "Angle": model.get("Angle"),
        "CornerRadius": model.get("CornerRadius", 0.0),
        "ContourGroups": model["ContourGroups"],
        "ZeroTraceIndices": model.get("ZeroTraceIndices", []),
        "Interfaces": model.get("Interfaces", []),
    }


def load_case(name, root):
    basis = np.atleast_2d(
        np.loadtxt(root / "basis-points.csv", delimiter=",", skiprows=1)
    )
    coefficients = np.atleast_1d(
        np.loadtxt(
            root / "heldout-coefficients.csv",
            delimiter=",",
            skiprows=1,
        )
    )
    if len(coefficients) != len(basis):
        raise ValueError(
            f"{root} has {len(basis)} basis points but "
            f"{len(coefficients)} held-out coefficients"
        )
    library = json.loads((root / "process-library.json").read_text())
    model = library["Models"][0]
    zero_trace_indices = np.asarray(
        model.get("ZeroTraceIndices", []),
        dtype=int,
    ) - 1
    if np.any(zero_trace_indices < 0) or np.any(zero_trace_indices >= len(basis)):
        raise ValueError(f"{root} has an invalid ZeroTraceIndices basis index")
    active_indices = np.delete(np.arange(len(basis)), zero_trace_indices)
    if not len(active_indices):
        raise ValueError(f"{root} has no free trace coefficients")
    if np.any(np.abs(coefficients[zero_trace_indices]) > 1.0e-14):
        raise ValueError(
            f"{root} held-out excitation is nonzero at a PEC-constrained trace knot"
        )
    interface_names = {
        int(interface["Coupon"]): interface["Type"]
        for interface in model.get("Interfaces", [])
    }

    responses = {}
    for kind in ("thin", "fabricated"):
        postpro = root / "postpro" / kind
        domain_path = postpro / "domain-response-matrix.csv"
        surface_path = postpro / "surface-response-matrix-aggregate.csv"
        if domain_path.is_file() != surface_path.is_file():
            raise FileNotFoundError(
                f"{postpro} has only one of the domain and surface response matrices"
            )
        domain = (
            read_domain_matrix(domain_path, len(basis))
            if domain_path.is_file()
            else None
        )
        surfaces = (
            read_surface_matrices(surface_path, len(basis))
            if surface_path.is_file()
            else {}
        )
        heldout = root / "postpro" / f"heldout-{kind}"
        direct_domain = float("nan")
        direct_surfaces = {}
        if (heldout / "domain-E.csv").is_file():
            direct_domain = read_row(heldout / "domain-E.csv")["E_elec (J)"]
        if (
            np.isfinite(direct_domain)
            and (heldout / "surface-Q.csv").is_file()
            and (heldout / "surface-Q-edge.csv").is_file()
        ):
            surface_row = read_row(heldout / "surface-Q.csv")
            edge_rows = read_interface_rows(heldout / "surface-Q-edge.csv")
            for interface in edge_rows:
                direct_surfaces[interface] = (
                    surface_row[f"p_surf[{interface}]"] * direct_domain
                    - edge_rows[interface]["E_out (J)"]
                )
        responses[kind] = {
            "domain": domain,
            "surfaces": surfaces,
            "direct_domain": direct_domain,
            "direct_surfaces": direct_surfaces,
        }
    return {
        "name": name,
        "root": root,
        "basis": basis,
        "coefficients": coefficients,
        "metadata": compatible_metadata(library),
        "active_indices": active_indices,
        "interface_names": interface_names,
        "responses": responses,
    }


def relative_change(value, previous):
    if previous is None or previous == 0.0:
        return float("nan")
    return value / previous - 1.0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        action="append",
        required=True,
        help="Ordered coupon level as NAME=CALIBRATION",
    )
    parser.add_argument(
        "--allow-missing-fabrication-metadata",
        action="store_true",
        help=(
            "Allow a legacy case without Fabrication metadata to be compared "
            "with an otherwise compatible version-3 library"
        ),
    )
    args = parser.parse_args()
    cases = [load_case(*parse_case(specification)) for specification in args.case]

    reference = cases[0]
    scale = max(np.max(np.abs(reference["basis"])), 1.0)
    for case in cases[1:]:
        reference_metadata = reference["metadata"].copy()
        case_metadata = case["metadata"].copy()
        if (
            args.allow_missing_fabrication_metadata
            and (
                reference_metadata["Fabrication"] is None
                or case_metadata["Fabrication"] is None
            )
        ):
            reference_metadata.pop("Fabrication")
            case_metadata.pop("Fabrication")
        if case_metadata != reference_metadata:
            raise ValueError(
                f"{case['root']} has incompatible process or model metadata"
            )
        if case["basis"].shape != reference["basis"].shape or not np.allclose(
            case["basis"], reference["basis"], rtol=0.0, atol=1.0e-12 * scale
        ):
            raise ValueError(f"{case['root']} uses a different response basis")
        if not np.allclose(
            case["coefficients"],
            reference["coefficients"],
            rtol=0.0,
            atol=1.0e-14,
        ):
            raise ValueError(
                f"{case['root']} uses different held-out coefficients"
            )
        if not np.array_equal(case["active_indices"], reference["active_indices"]):
            raise ValueError(
                f"{case['root']} uses different PEC-constrained trace knots"
            )

    print(
        "case,kind,quantity,matrix_norm_J,matrix_change_pct,"
        "heldout_matrix_energy_J,heldout_matrix_change_pct,"
        "heldout_direct_energy_J,heldout_direct_change_pct,"
        "heldout_basis_error_pct"
    )
    previous = {}
    coefficients = reference["coefficients"]
    active = np.ix_(
        reference["active_indices"],
        reference["active_indices"],
    )
    for case in cases:
        for kind, response in case["responses"].items():
            quantities = [
                (
                    "domain",
                    response["domain"],
                    response["direct_domain"],
                )
            ]
            quantities.extend(
                (
                    case["interface_names"].get(
                        interface, f"interface-{interface}"
                    ),
                    response["surfaces"].get(interface),
                    response["direct_surfaces"].get(interface, float("nan")),
                )
                for interface in sorted(
                    response["surfaces"].keys()
                    | response["direct_surfaces"].keys()
                )
            )
            for quantity, matrix, direct in quantities:
                key = (kind, quantity)
                old = previous.get(key)
                reduced_matrix = (
                    matrix[active] if matrix is not None else None
                )
                matrix_norm = (
                    np.linalg.norm(reduced_matrix)
                    if reduced_matrix is not None
                    else float("nan")
                )
                matrix_energy = (
                    coefficients @ matrix @ coefficients
                    if matrix is not None
                    else float("nan")
                )
                matrix_change = (
                    float("nan")
                    if (
                        matrix is None
                        or old is None
                        or old["matrix"] is None
                    )
                    else np.linalg.norm(reduced_matrix - old["matrix"])
                    / max(matrix_norm, np.finfo(float).tiny)
                )
                energy_change = relative_change(
                    matrix_energy,
                    None if old is None else old["matrix_energy"],
                )
                direct_change = relative_change(
                    direct,
                    None if old is None else old["direct"],
                )
                basis_error = (
                    matrix_energy / direct - 1.0
                    if np.isfinite(direct) and direct != 0.0
                    else float("nan")
                )
                print(
                    f"{case['name']},{kind},{quantity},{matrix_norm:.12e},"
                    f"{100.0 * matrix_change:+.6f},"
                    f"{matrix_energy:.12e},{100.0 * energy_change:+.6f},"
                    f"{direct:.12e},{100.0 * direct_change:+.6f},"
                    f"{100.0 * basis_error:+.6f}"
                )
                previous[key] = {
                    "matrix": reduced_matrix,
                    "matrix_energy": matrix_energy,
                    "direct": direct,
                }


if __name__ == "__main__":
    main()
