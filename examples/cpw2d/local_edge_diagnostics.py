#!/usr/bin/env python3

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

from local_edge_correction import (
    TOPOLOGY_COLUMNS,
    TOPOLOGY_ID_COLUMNS,
    VERTEX_ENERGY_COLUMNS,
    edge_model_confidence,
    fit_singular_amplitude,
    fit_surface_annulus,
)


PARTICIPATION_RE = re.compile(r"p_surf\[(\d+)\](?: \((\d+)\))?$")


def read_table(path):
    with path.open(newline="") as stream:
        reader = csv.reader(stream)
        header = [cell.strip() for cell in next(reader)]
        rows = [
            dict(zip(header, (cell.strip() for cell in row)))
            for row in reader
            if any(cell.strip() for cell in row)
        ]
    if not rows:
        raise RuntimeError(f"No data rows found in {path}")
    return header, rows


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Validate localized Palace edge diagnostics and write per-edge annular "
            "participation densities."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Palace output directory containing surface-Q-edge CSV files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output CSV path (default: INPUT/surface-Q-edge-line.csv).",
    )
    parser.add_argument(
        "--tolerance",
        type=float,
        default=1.0e-8,
        help="Maximum relative mismatch allowed between aggregate and localized annuli.",
    )
    parser.add_argument(
        "--fit-count",
        type=int,
        default=3,
        help="Number of smallest radii used to fit each local singular amplitude.",
    )
    parser.add_argument(
        "--fit-residual-scale",
        type=float,
        default=0.1,
        help="Relative RMS at which edge-model confidence is reduced by one half.",
    )
    args = parser.parse_args()
    if args.fit_count < 3:
        parser.error("--fit-count must be at least 3")
    if args.fit_residual_scale < 0.0:
        parser.error("--fit-residual-scale must be nonnegative")
    output = args.output or args.input / "surface-Q-edge-line.csv"

    aggregate_header, aggregate_rows = read_table(
        args.input / "surface-Q-edge.csv"
    )
    surface_header, surface_rows = read_table(args.input / "surface-Q.csv")
    local_header, local_rows = read_table(
        args.input / "surface-Q-edge-local.csv"
    )
    topology_id_columns = TOPOLOGY_ID_COLUMNS & set(local_header)
    if topology_id_columns and not TOPOLOGY_COLUMNS <= set(local_header):
        raise RuntimeError(
            "Localized table has incomplete automatic edge topology columns; "
            f"missing={sorted(TOPOLOGY_COLUMNS - set(local_header))}"
        )
    has_topology = TOPOLOGY_COLUMNS <= set(local_header)
    vertex_energy_columns = VERTEX_ENERGY_COLUMNS & set(local_header)
    if vertex_energy_columns and vertex_energy_columns != VERTEX_ENERGY_COLUMNS:
        raise RuntimeError(
            "Localized table has incomplete non-regular vertex energy columns; "
            f"missing={sorted(VERTEX_ENERGY_COLUMNS - vertex_energy_columns)}"
        )
    has_vertex_energies = vertex_energy_columns == VERTEX_ENERGY_COLUMNS
    state_column = aggregate_header[0]
    if surface_header[0] != state_column or local_header[0] != state_column:
        raise RuntimeError(
            f"State columns differ between surface tables: {state_column!r}, "
            f"{surface_header[0]!r}, and {local_header[0]!r}"
        )

    excitations = defaultdict(set)
    for row in aggregate_rows:
        excitations[(row[state_column], int(row["interface"]))].add(int(row["exc"]))

    surface = {}
    for row in surface_rows:
        for label in surface_header:
            match = PARTICIPATION_RE.fullmatch(label)
            if not match:
                continue
            interface = int(match.group(1))
            if match.group(2) is None:
                candidates = excitations[(row[state_column], interface)]
                if len(candidates) != 1:
                    raise RuntimeError(
                        f"Cannot infer excitation for {state_column}={row[state_column]}, "
                        f"interface={interface}"
                    )
                excitation = next(iter(candidates))
            else:
                excitation = int(match.group(2))
            surface[(row[state_column], excitation, interface)] = float(row[label])

    aggregate = {}
    for row in aggregate_rows:
        key = (
            row[state_column],
            int(row["exc"]),
            int(row["interface"]),
            float(row["R (m)"]),
        )
        aggregate[key] = (float(row["p_out"]), float(row["p_ann"]))

    has_total = "p_total" in local_header
    has_inside = "p_in" in local_header
    localized_total = defaultdict(float)
    localized_inside = defaultdict(float)
    localized_annulus = defaultdict(float)
    local_edge_values = defaultdict(dict)
    for row in local_rows:
        edge_key = (
            row[state_column],
            int(row["exc"]),
            int(row["interface"]),
            int(row["edge"]),
        )
        radius = float(row["R (m)"])
        local_edge_values[edge_key][radius] = {
            "p_ann": float(row["p_ann"]),
            "p_bulk_ann": float(row["p_bulk_ann"]),
        }
    local_bulk_fits = {
        key: fit_singular_amplitude(values, args.fit_count)
        for key, values in local_edge_values.items()
    }
    local_surface_fits = {
        key: fit_surface_annulus(values, args.fit_count)
        for key, values in local_edge_values.items()
    }

    reduced = []
    for row in local_rows:
        key = (
            row[state_column],
            int(row["exc"]),
            int(row["interface"]),
            float(row["R (m)"]),
        )
        if has_total:
            localized_total[key] += float(row["p_total"])
        if has_inside:
            localized_inside[key] += float(row["p_in"])
        localized_annulus[key] += float(row["p_ann"])
        p0 = [float(row[f"{axis}0 (m)"]) for axis in "xyz"]
        p1 = [float(row[f"{axis}1 (m)"]) for axis in "xyz"]
        length = float(row["L_edge (m)"])
        total_participation = float(row["p_total"]) if has_total else math.nan
        inside_participation = float(row["p_in"]) if has_inside else math.nan
        participation = float(row["p_ann"])
        bulk_participation = float(row["p_bulk_ann"])
        radius = float(row["R (m)"])
        edge_key = (*key[:3], int(row["edge"]))
        amplitude, bulk_residual, _ = local_bulk_fits[edge_key]
        _, smooth_density, surface_residual, _ = local_surface_fits[edge_key]
        singular_fraction, model_confidence = edge_model_confidence(
            amplitude,
            bulk_participation / radius,
            bulk_residual,
            surface_residual,
            args.fit_residual_scale,
        )
        topology = {}
        if has_topology:
            vertex_types = (int(row["v0_type"]), int(row["v1_type"]))
            topology = {
                "automatic": int(row["automatic"]),
                "component": int(row["component"]),
                "chain": int(row["chain"]),
                "v0_type": vertex_types[0],
                "v1_type": vertex_types[1],
                "touches_nonregular_vertex": int(
                    any(vertex_type != 0 for vertex_type in vertex_types)
                ),
            }
        vertex_energies = {}
        if has_vertex_energies:
            vertex_energies = {
                "p_vertex_in": float(row["p_vertex_in"]),
                "p_bulk_vertex_ann": float(row["p_bulk_vertex_ann"]),
            }
        reduced.append(
            {
                state_column: row[state_column],
                "exc": int(row["exc"]),
                "interface": int(row["interface"]),
                "edge": int(row["edge"]),
                **topology,
                "R (m)": float(row["R (m)"]),
                "x_mid (m)": 0.5 * (p0[0] + p1[0]),
                "y_mid (m)": 0.5 * (p0[1] + p1[1]),
                "z_mid (m)": 0.5 * (p0[2] + p1[2]),
                "L_edge (m)": length,
                "p_total": total_participation,
                "p_total_per_m": (
                    total_participation / length if length > 0.0 else math.nan
                ),
                "p_in": inside_participation,
                "p_in_per_m": (
                    inside_participation / length if length > 0.0 else math.nan
                ),
                **vertex_energies,
                "p_ann": participation,
                "p_ann_per_m": participation / length if length > 0.0 else math.nan,
                "p_bulk_ann": bulk_participation,
                "p_bulk_ann_over_R (1/m)": bulk_participation / radius,
                "p_bulk_ann_over_RL (1/m^2)": (
                    bulk_participation / (radius * length)
                    if length > 0.0
                    else math.nan
                ),
                "singular_amplitude_per_m": amplitude,
                "smooth_surface_density_per_m": smooth_density,
                "bulk_fit_relative_rms": bulk_residual,
                "surface_fit_relative_rms": surface_residual,
                "fit_relative_rms": max(bulk_residual, surface_residual),
                "singular_fraction": singular_fraction,
                "model_confidence": model_confidence,
            }
        )

    if aggregate.keys() != localized_annulus.keys():
        missing = aggregate.keys() - localized_annulus.keys()
        extra = localized_annulus.keys() - aggregate.keys()
        raise RuntimeError(
            f"Aggregate and localized keys differ; missing={sorted(missing)}, "
            f"extra={sorted(extra)}"
        )

    worst_annulus = max(
        (
            abs(localized_annulus[key] - value[1]) / max(abs(value[1]), 1.0e-300),
            key,
            value[1],
            localized_annulus[key],
        )
        for key, value in aggregate.items()
    )
    worst_inside = None
    if has_inside:
        worst_inside = max(
            (
                abs(localized_inside[key] - (surface[key[:3]] - value[0]))
                / max(abs(surface[key[:3]] - value[0]), 1.0e-300),
                key,
                surface[key[:3]] - value[0],
                localized_inside[key],
            )
            for key, value in aggregate.items()
        )
    worst_total = None
    if has_total:
        worst_total = max(
            (
                abs(localized_total[key] - surface[key[:3]])
                / max(abs(surface[key[:3]]), 1.0e-300),
                key,
                surface[key[:3]],
                localized_total[key],
            )
            for key in aggregate
        )
    worst = max(
        candidate
        for candidate in (worst_annulus, worst_inside, worst_total)
        if candidate is not None
    )
    if worst[0] > args.tolerance:
        raise RuntimeError(
            f"Localized edge energies do not conserve the aggregate result: mismatch "
            f"{worst[0]:.3e} for {worst[1]}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=reduced[0].keys())
        writer.writeheader()
        writer.writerows(reduced)

    edge_count = len({(row["interface"], row["edge"]) for row in reduced})
    print(
        f"Wrote {len(reduced)} rows for {edge_count} interface-edge pairs to {output}"
    )
    print(
        f"Worst annulus conservation mismatch: {worst_annulus[0]:.3e}; "
        + (
            f"worst inside mismatch: {worst_inside[0]:.3e}"
            if worst_inside
            else "inside conservation unavailable (legacy table without p_in)"
        )
        + (
            f"; worst total mismatch: {worst_total[0]:.3e}"
            if worst_total
            else "; total conservation unavailable (legacy table without p_total)"
        )
    )


if __name__ == "__main__":
    main()
