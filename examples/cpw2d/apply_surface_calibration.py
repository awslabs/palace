#!/usr/bin/env python3

import argparse
import csv
import math
import re
from pathlib import Path


PARTICIPATION_RE = re.compile(r"p_surf\[(\d+)\]$")
QUALITY_FACTOR_RE = re.compile(r"Q_surf\[(\d+)\]$")


def read_table(path):
    with path.open(newline="") as stream:
        rows = list(csv.reader(stream))
    if len(rows) < 2:
        raise RuntimeError(f"No data rows found in {path}")
    header = [cell.strip() for cell in rows[0]]
    data = [
        dict(zip(header, (cell.strip() for cell in row)))
        for row in rows[1:]
        if any(cell.strip() for cell in row)
    ]
    return header, data


def read_surface_participations(directory):
    path = directory / "surface-Q.csv"
    header, rows = read_table(path)
    state_column = header[0]
    states = {}
    for row in rows:
        participations = {}
        quality_factors = {}
        for label in header:
            participation_match = PARTICIPATION_RE.fullmatch(label)
            quality_factor_match = QUALITY_FACTOR_RE.fullmatch(label)
            if participation_match:
                participations[int(participation_match.group(1))] = float(row[label])
            elif quality_factor_match:
                quality_factors[int(quality_factor_match.group(1))] = float(row[label])
        if not participations:
            raise RuntimeError(f"No participation columns found in {path}")
        if participations.keys() != quality_factors.keys():
            raise RuntimeError(
                f"Participation and quality-factor columns differ in {path}"
            )
        states[row[state_column]] = (participations, quality_factors)
    return state_column, states


def read_edge_participations(directory):
    path = directory / "surface-Q-edge.csv"
    header, rows = read_table(path)
    state_column = header[0]
    states = {}
    for row in rows:
        if int(row["exc"]) != 0:
            continue
        values = states.setdefault(row[state_column], {})
        key = (int(row["interface"]), float(row["R (m)"]))
        values[key] = (float(row["p_out"]), float(row["p_ann"]))
    if not states:
        raise RuntimeError(f"No edge data found in {path}")
    return state_column, states


def read_calibration(path):
    _, rows = read_table(path)
    required = {
        "interface",
        "interface_index",
        "radius_um",
        "proxy_interface",
        "outside_scale",
        "transfer_coefficient",
    }
    missing = required - rows[0].keys()
    if missing:
        raise RuntimeError(
            f"Calibration table {path} is missing columns: {sorted(missing)}"
        )
    return rows


def matching_edge_key(edge_values, interface, radius):
    candidates = [key for key in edge_values if key[0] == interface]
    if not candidates:
        raise RuntimeError(f"No edge diagnostics found for interface {interface}")
    key = min(candidates, key=lambda candidate: abs(candidate[1] - radius))
    if not math.isclose(key[1], radius, rel_tol=1.0e-8, abs_tol=1.0e-15):
        return None
    return key


def apply_calibration(calibration, input_directory):
    state_column, surface_states = read_surface_participations(input_directory)
    edge_state_column, edge_states = read_edge_participations(input_directory)
    if state_column != edge_state_column:
        raise RuntimeError(
            f"State columns differ between surface tables: "
            f"{state_column!r} and {edge_state_column!r}"
        )
    name_to_index = {
        row["interface"]: int(row["interface_index"]) for row in calibration
    }

    results = []
    for state, (raw, raw_quality) in surface_states.items():
        if state not in edge_states:
            raise RuntimeError(f"No edge diagnostics found for {state_column}={state}")
        edge = edge_states[state]
        state_results = []
        for coefficients in calibration:
            interface = int(coefficients["interface_index"])
            proxy_name = coefficients["proxy_interface"]
            if proxy_name not in name_to_index:
                raise RuntimeError(
                    f"Proxy interface {proxy_name!r} is absent from the calibration table"
                )
            proxy_interface = name_to_index[proxy_name]
            radius_um = float(coefficients["radius_um"])
            radius = 1.0e-6 * radius_um
            interface_key = matching_edge_key(edge, interface, radius)
            proxy_key = matching_edge_key(edge, proxy_interface, radius)
            if interface_key is None or proxy_key is None:
                continue
            participation_outside, _ = edge[interface_key]
            _, proxy_annulus = edge[proxy_key]
            outside_scale = float(coefficients["outside_scale"])
            transfer_coefficient = float(coefficients["transfer_coefficient"])
            corrected = (
                outside_scale * participation_outside
                + transfer_coefficient * proxy_annulus
            )

            raw_participation = raw[interface]
            raw_q = raw_quality[interface]
            loss_tangent = 1.0 / (raw_participation * raw_q)
            corrected_q = 1.0 / (corrected * loss_tangent)
            state_results.append(
                {
                    "state_column": state_column,
                    "state": state,
                    "interface": coefficients["interface"],
                    "interface_index": interface,
                    "radius_um": radius_um,
                    "p_raw": raw_participation,
                    "p_outside": participation_outside,
                    "proxy_interface": proxy_name,
                    "p_proxy_annulus": proxy_annulus,
                    "outside_scale": outside_scale,
                    "transfer_coefficient": transfer_coefficient,
                    "p_corrected": corrected,
                    "Q_raw": raw_q,
                    "Q_corrected": corrected_q,
                }
            )
        missing_interfaces = {
            int(row["interface_index"]) for row in calibration
        } - {row["interface_index"] for row in state_results}
        if missing_interfaces:
            raise RuntimeError(
                f"No matching calibration radii found for {state_column}={state}, "
                f"interfaces: {sorted(missing_interfaces)}"
            )
        results.extend(state_results)
    return results


def write_results(path, results):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)


def print_results(results):
    print(
        f"{'state':>8} {'interface':>9} {'R (um)':>8} {'p_raw':>14} "
        f"{'p_corrected':>14} {'Q_corrected':>14}"
    )
    for row in results:
        print(
            f"{row['state']:>8} {row['interface']:>9} {row['radius_um']:8.3f} "
            f"{row['p_raw']:14.6e} {row['p_corrected']:14.6e} "
            f"{row['Q_corrected']:14.6e}"
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Apply fabrication-stack edge coefficients to a Palace thin-metal "
            "surface-participation result."
        )
    )
    parser.add_argument(
        "--calibration",
        type=Path,
        required=True,
        help="Calibration CSV produced by surface_calibration.py.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Palace output directory containing surface-Q CSV files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Corrected participation CSV to write.",
    )
    args = parser.parse_args()

    results = apply_calibration(read_calibration(args.calibration), args.input)
    write_results(args.output, results)
    print_results(results)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
