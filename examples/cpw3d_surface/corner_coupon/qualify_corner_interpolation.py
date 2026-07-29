#!/usr/bin/env python3

"""Qualify interpolation between two rounded-corner response coupons."""

import argparse
import json
import math
from pathlib import Path

import numpy as np

import compare_probe_convergence as response


KINDS = ("thin", "fabricated")


def model_metadata(library):
    if len(library.get("Models", [])) != 1:
        raise ValueError("Corner interpolation inputs require exactly one model")
    model = library["Models"][0]
    return {
        "MatchingRadius": float(library["MatchingRadius"]),
        "Fabrication": library.get("Fabrication"),
        "Topology": model["Topology"],
        "Angle": float(model["Angle"]),
        "ContourGroups": model["ContourGroups"],
        "Interfaces": model.get("Interfaces", []),
        "BoundaryCondition": model.get("BoundaryCondition", "PEC"),
    }


def load_case(name, root):
    root = root.expanduser().resolve()
    library = json.loads((root / "process-library.json").read_text())
    model = library["Models"][0]
    radius = float(model.get("CornerRadius", 0.0))
    if model.get("Topology") not in ("ConvexCorner", "ConcaveCorner") or radius <= 0.0:
        raise ValueError(f"{root} is not a positive-radius corner calibration")

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
        raise ValueError(f"{root} has incompatible basis and held-out trace sizes")

    zero = np.asarray(model.get("ZeroTraceIndices", []), dtype=int) - 1
    if np.any(zero < 0) or np.any(zero >= len(basis)):
        raise ValueError(f"{root} has an invalid ZeroTraceIndices entry")
    active = np.delete(np.arange(len(basis)), zero)
    if not len(active):
        raise ValueError(f"{root} has no free response trace coefficients")

    interface_names = {
        int(interface["Coupon"]): interface["Type"]
        for interface in model.get("Interfaces", [])
    }
    responses = {}
    for kind in KINDS:
        postpro = root / "postpro" / kind
        domain = response.read_matrix(
            postpro / "domain-response-matrix.csv",
            len(basis),
            "Q_ij (J)",
        )
        surfaces = response.read_surface_matrices(
            postpro / "surface-response-matrix-aggregate.csv",
            len(basis),
        )
        if surfaces.keys() != interface_names.keys():
            raise ValueError(f"{root} has inconsistent interface response matrices")
        responses[kind] = {"domain": domain, "surfaces": surfaces}

    return {
        "Name": name,
        "Root": root,
        "ModelName": model["Name"],
        "Radius": radius,
        "Metadata": model_metadata(library),
        "Basis": basis,
        "Coefficients": coefficients,
        "Active": active,
        "Responses": responses,
        "InterfaceNames": interface_names,
    }


def relative_error(value, reference):
    if abs(reference) <= np.finfo(float).tiny:
        return float("nan")
    return 100.0 * (value / reference - 1.0)


def matrix_error(value, reference):
    scale = max(np.linalg.norm(reference), np.finfo(float).tiny)
    return 100.0 * np.linalg.norm(value - reference) / scale


def fixed_flux(case, coefficients=None):
    fabricated = case["Responses"]["fabricated"]["domain"]
    thin = case["Responses"]["thin"]["domain"]
    active = case["Active"]
    if coefficients is None:
        coefficients = case["Coefficients"]
    result = np.zeros(len(case["Basis"]))
    result[active] = np.linalg.solve(
        fabricated[np.ix_(active, active)],
        thin[np.ix_(active, active)] @ coefficients[active],
    )
    return result


def qualify_cases(lower, heldout, upper, max_matrix_error, max_energy_error):
    if not lower["Radius"] < heldout["Radius"] < upper["Radius"]:
        raise ValueError("Held-out corner radius must lie strictly inside the bracket")
    if lower["Metadata"] != heldout["Metadata"] or upper["Metadata"] != heldout["Metadata"]:
        raise ValueError("Corner interpolation inputs use incompatible process metadata")

    scale = max(np.max(np.abs(heldout["Basis"])), 1.0)
    for case in (lower, upper):
        if case["Basis"].shape != heldout["Basis"].shape or not np.allclose(
            case["Basis"],
            heldout["Basis"],
            rtol=0.0,
            atol=1.0e-12 * scale,
        ):
            raise ValueError("Corner interpolation inputs use different response bases")

    upper_weight = (
        (heldout["Radius"] - lower["Radius"])
        / (upper["Radius"] - lower["Radius"])
    )
    weights = (1.0 - upper_weight, upper_weight)
    endpoints = (lower, upper)
    checks = []

    def append_check(kind, quantity, metric, value, limit, checked):
        passed = math.isfinite(value) and (not checked or abs(value) <= limit)
        checks.append(
            {
                "Kind": kind,
                "Quantity": quantity,
                "Metric": metric,
                "ValuePercent": value,
                "LimitPercent": limit if checked else None,
                "Checked": checked,
                "Passed": passed,
            }
        )

    for kind in KINDS:
        exact = heldout["Responses"][kind]
        interpolated_domain = sum(
            weight * case["Responses"][kind]["domain"]
            for weight, case in zip(weights, endpoints)
        )
        append_check(
            kind,
            "domain",
            "MatrixError",
            matrix_error(interpolated_domain, exact["domain"]),
            max_matrix_error,
            True,
        )
        for interface, exact_surface in exact["surfaces"].items():
            interpolated_surface = sum(
                weight * case["Responses"][kind]["surfaces"][interface]
                for weight, case in zip(weights, endpoints)
            )
            append_check(
                kind,
                heldout["InterfaceNames"][interface],
                "MatrixError",
                matrix_error(interpolated_surface, exact_surface),
                max_matrix_error,
                kind == "fabricated",
            )

    exact_defect = (
        heldout["Responses"]["fabricated"]["domain"]
        - heldout["Responses"]["thin"]["domain"]
    )
    interpolated_defect = sum(
        weight
        * (
            case["Responses"]["fabricated"]["domain"]
            - case["Responses"]["thin"]["domain"]
        )
        for weight, case in zip(weights, endpoints)
    )
    append_check(
        "defect",
        "domain",
        "MatrixError",
        matrix_error(interpolated_defect, exact_defect),
        max_matrix_error,
        True,
    )
    coefficients = heldout["Coefficients"]
    append_check(
        "defect",
        "domain",
        "FixedTraceEnergyError",
        relative_error(
            coefficients @ interpolated_defect @ coefficients,
            coefficients @ exact_defect @ coefficients,
        ),
        max_energy_error,
        True,
    )

    exact_flux = fixed_flux(heldout)
    endpoint_flux = [fixed_flux(case, coefficients) for case in endpoints]
    for interface, exact_surface in heldout["Responses"]["fabricated"][
        "surfaces"
    ].items():
        interpolated_surface = sum(
            weight * case["Responses"]["fabricated"]["surfaces"][interface]
            for weight, case in zip(weights, endpoints)
        )
        append_check(
            "fabricated",
            heldout["InterfaceNames"][interface],
            "FixedTraceEnergyError",
            relative_error(
                coefficients @ interpolated_surface @ coefficients,
                coefficients @ exact_surface @ coefficients,
            ),
            max_energy_error,
            True,
        )
        interpolated_flux_energy = sum(
            weight
            * flux
            @ case["Responses"]["fabricated"]["surfaces"][interface]
            @ flux
            for weight, case, flux in zip(weights, endpoints, endpoint_flux)
        )
        append_check(
            "fabricated",
            heldout["InterfaceNames"][interface],
            "FixedFluxEnergyError",
            relative_error(
                interpolated_flux_energy,
                exact_flux @ exact_surface @ exact_flux,
            ),
            max_energy_error,
            True,
        )

    passed = all(check["Passed"] for check in checks)
    maximum_matrix_error = max(
        abs(check["ValuePercent"])
        for check in checks
        if check["Checked"] and check["Metric"] == "MatrixError"
    )
    maximum_energy_error = max(
        abs(check["ValuePercent"])
        for check in checks
        if check["Checked"] and check["Metric"] != "MatrixError"
    )
    return {
        "Version": 1,
        "Study": "CornerRadiusInterpolation",
        "Cases": {
            "Lower": str(lower["Root"]),
            "Heldout": str(heldout["Root"]),
            "Upper": str(upper["Root"]),
        },
        "Radii": {
            "Lower": lower["Radius"],
            "Heldout": heldout["Radius"],
            "Upper": upper["Radius"],
        },
        "Weights": {
            "Lower": weights[0],
            "Upper": weights[1],
        },
        "ThresholdsPercent": {
            "MatrixError": max_matrix_error,
            "EnergyError": max_energy_error,
        },
        "Checks": checks,
        "MaximumMatrixErrorPercent": maximum_matrix_error,
        "MaximumEnergyErrorPercent": maximum_energy_error,
        "LibraryRecord": {
            "LowerModel": lower["ModelName"],
            "UpperModel": upper["ModelName"],
            "Qualification": {
                "Method": "HeldOutCoupon",
                "Passed": passed,
                "HeldoutRadius": heldout["Radius"],
                "MaximumMatrixErrorPercent": maximum_matrix_error,
                "MaximumEnergyErrorPercent": maximum_energy_error,
            },
        },
        "Passed": passed,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--lower", type=Path, required=True)
    parser.add_argument("--heldout", type=Path, required=True)
    parser.add_argument("--upper", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-matrix-error", type=float, default=5.0)
    parser.add_argument("--max-energy-error", type=float, default=10.0)
    args = parser.parse_args()
    if args.max_matrix_error <= 0.0 or args.max_energy_error <= 0.0:
        parser.error("Interpolation error thresholds must be positive")

    report = qualify_cases(
        load_case("lower", args.lower),
        load_case("heldout", args.heldout),
        load_case("upper", args.upper),
        args.max_matrix_error,
        args.max_energy_error,
    )
    destination = args.output.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n")
    print(destination)
    print(
        "corner radius interpolation: "
        f"{'PASS' if report['Passed'] else 'FAIL'} "
        f"(matrix={report['MaximumMatrixErrorPercent']:.3f}%, "
        f"energy={report['MaximumEnergyErrorPercent']:.3f}%)"
    )
    if not report["Passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
