#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Compare an extruded CPW capacitance with the infinite-domain reference."""

import argparse
import csv
import json
import math
from pathlib import Path


VACUUM_PERMITTIVITY = 8.8541878128e-12


def complete_elliptic_integral_first_kind(modulus):
    if not 0.0 <= modulus < 1.0:
        raise ValueError("Elliptic modulus must satisfy 0 <= k < 1")
    arithmetic = 1.0
    geometric = math.sqrt(1.0 - modulus * modulus)
    while abs(arithmetic - geometric) > 4.0 * math.ulp(arithmetic):
        arithmetic, geometric = (
            0.5 * (arithmetic + geometric),
            math.sqrt(arithmetic * geometric),
        )
    return math.pi / (2.0 * arithmetic)


def reference_capacitance_per_length(strip_half_width, ground_inner, substrate_epsilon):
    modulus = strip_half_width / ground_inner
    complementary_modulus = math.sqrt(1.0 - modulus * modulus)
    effective_epsilon = 0.5 * (substrate_epsilon + 1.0)
    return (
        4.0
        * VACUUM_PERMITTIVITY
        * effective_epsilon
        * complete_elliptic_integral_first_kind(modulus)
        / complete_elliptic_integral_first_kind(complementary_modulus)
    )


def read_capacitance(path):
    with path.open(newline="", encoding="ascii") as stream:
        rows = list(csv.reader(stream))
    if len(rows) != 2 or len(rows[1]) != 2:
        raise ValueError(f"Unexpected one-terminal capacitance table: {path}")
    return float(rows[1][1])


def read_table(path):
    with path.open(newline="", encoding="ascii") as stream:
        reader = csv.reader(stream)
        rows = list(reader)
    if len(rows) < 2:
        raise ValueError(f"Empty table: {path}")
    header = [entry.strip() for entry in rows[0]]
    if len(header) != len(set(header)):
        raise ValueError(f"Duplicate columns in table: {path}")
    return [
        {name: value.strip() for name, value in zip(header, row, strict=True)}
        for row in rows[1:]
    ]


def read_metadata(postpro):
    path = postpro / "palace.json"
    with path.open(encoding="ascii") as stream:
        return json.load(stream)


def read_energy_consistency(postpro, capacitance):
    row = read_table(postpro / "domain-E.csv")[0]
    energy = float(row["E_elec (J)"])
    if not math.isfinite(energy) or energy <= 0.0:
        raise ValueError(f"Invalid electrostatic energy in {postpro / 'domain-E.csv'}")
    voltage_row = read_table(postpro / "terminal-V.csv")[0]
    voltage = float(voltage_row["V_inc[i] (V)"])
    if not math.isfinite(voltage) or voltage == 0.0:
        raise ValueError(f"Invalid terminal voltage in {postpro / 'terminal-V.csv'}")
    return energy, voltage, abs(2.0 * energy / (capacitance * voltage * voltage) - 1.0)


def read_slope_summary(postpro):
    path = postpro / "singular-edge-slopes.csv"
    if not path.is_file():
        return None
    rows = read_table(path)
    valid = [row for row in rows if int(row["fit_valid"]) != 0]
    if not valid:
        raise ValueError(f"No valid singular edge-slope fits in {path}")
    errors = [
        float(row["fitted_slope"]) - float(row["expected_slope"]) for row in valid
    ]
    r_squared = [float(row["R_squared"]) for row in valid]
    values = errors + r_squared
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"Nonfinite singular edge-slope diagnostic in {path}")
    return {
        "count": len(valid),
        "total": len(rows),
        "rms_error": math.sqrt(sum(error * error for error in errors) / len(errors)),
        "maximum_error": max(abs(error) for error in errors),
        "minimum_r_squared": min(r_squared),
    }


def report_result(postpro, reference, length):
    capacitance = read_capacitance(postpro / "terminal-C.csv")
    per_length = capacitance / length
    relative_error = per_length / reference - 1.0
    energy, voltage, energy_discrepancy = read_energy_consistency(postpro, capacitance)
    metadata = read_metadata(postpro)
    problem = metadata["Problem"]
    iterations = metadata["LinearSolver"]["TotalIts"]
    converged = metadata["LinearSolver"].get("Converged")
    operator_time = metadata["ElapsedTime"]["Durations"]["OperatorConstruction"]
    singular = metadata.get("SingularElements")
    standard_dofs = int(problem["DegreesOfFreedom"])
    enrichment_dofs = (
        int(singular["H1EnrichmentDegreesOfFreedom"]) if singular is not None else 0
    )
    combined_dofs = standard_dofs + enrichment_dofs

    print(
        f"{postpro}: C = {capacitance:.12e} F, "
        f"C' = {per_length:.12e} F/m, error = {relative_error:+.6%}"
    )
    print(
        f"  mesh = {problem['MeshElements']} tets, "
        f"DOFs = {combined_dofs} ({standard_dofs} standard"
        f"{f' + {enrichment_dofs} enrichment' if enrichment_dofs else ''}), "
        f"iterations = {iterations}, "
        f"operator = {operator_time:.3f} s, "
        f"|2 E/(C V^2) - 1| = {energy_discrepancy:.3e} "
        f"(E = {energy:.12e} J, V = {voltage:.6e} V)"
    )
    if converged is False:
        print(
            "  LINEAR SOLVE DID NOT CONVERGE: "
            f"{metadata['LinearSolver']['FailedSolves']} failed solve(s)"
        )

    slopes = read_slope_summary(postpro)
    if singular is None:
        if slopes is not None:
            raise ValueError(f"Slope output exists without singular metadata: {postpro}")
        return
    if slopes is None:
        raise ValueError(f"Singular metadata exists without slope output: {postpro}")

    quadrature = singular["Quadrature"]
    scaling = singular["StiffnessDiagonalScaling"]
    print(
        f"  singular order = {singular['SingularOrder']}, "
        f"H1 enrichment DOFs = {singular['H1EnrichmentDegreesOfFreedom']}, "
        f"leaves = {quadrature['TotalLeafCount']}, "
        f"max depth = {quadrature['MaximumDepth']}, "
        f"diagonal spread = {scaling['CombinedSpread']:.3e}"
    )
    consistency = singular.get("EnergyConsistency")
    if consistency is not None:
        if len(consistency) != 1:
            raise ValueError(
                f"Expected one singular energy-consistency record in {postpro}"
            )
        certificate = consistency[0]
        recorded_discrepancy = certificate[
            "DirectToAlgebraicRelativeDiscrepancy"
        ]
        if not math.isclose(
            energy_discrepancy, recorded_discrepancy, rel_tol=1.0e-8, abs_tol=1.0e-12
        ):
            raise ValueError(
                "Rounded CSV energy discrepancy disagrees with the full-precision "
                f"certificate in {postpro}: {energy_discrepancy} versus "
                f"{recorded_discrepancy}"
            )
        print(
            "  energy certificate = "
            f"{recorded_discrepancy:.3e} discrepancy <= "
            f"{certificate['CertifiedRelativeErrorBound']:.3e} bound "
            f"(integration {certificate['IntegrationRelativeErrorBound']:.3e}, "
            f"assembly {certificate['AssemblyRelativeErrorBound']:.3e}, "
            f"floating point {certificate['FloatingPointRelativeAllowance']:.3e})"
        )
    print(
        f"  slopes = {slopes['count']}/{slopes['total']} valid, "
        f"RMS error = {slopes['rms_error']:.3e}, "
        f"max error = {slopes['maximum_error']:.3e}, "
        f"minimum R^2 = {slopes['minimum_r_squared']:.6f}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("postpro", nargs="+", type=Path)
    parser.add_argument("--strip-half-width", type=float, default=5.0)
    parser.add_argument("--ground-inner", type=float, default=8.0)
    parser.add_argument("--substrate-epsilon", type=float, default=11.45)
    parser.add_argument("--length", type=float, default=10.0)
    parser.add_argument("--length-unit", type=float, default=1.0e-6)
    args = parser.parse_args()

    reference = reference_capacitance_per_length(
        args.strip_half_width, args.ground_inner, args.substrate_epsilon
    )
    print(f"Infinite-domain reference C' = {reference:.12e} F/m")
    for postpro in args.postpro:
        report_result(
            postpro,
            reference,
            args.length * args.length_unit,
        )


if __name__ == "__main__":
    main()
