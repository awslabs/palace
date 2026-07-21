#!/usr/bin/env python3

import argparse
import csv
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path


PARTICIPATION_RE = re.compile(r"p_surf\[(\d+)\](?: \((\d+)\))?$")
QUALITY_FACTOR_RE = re.compile(r"Q_surf\[(\d+)\](?: \((\d+)\))?$")
TOPOLOGY_ID_COLUMNS = {
    "automatic",
    "component",
    "chain",
    "v0_type",
    "v1_type",
}
TOPOLOGY_COLUMNS = TOPOLOGY_ID_COLUMNS | {"L_edge (m)"}
VERTEX_ENERGY_COLUMNS = {"p_vertex_in", "p_bulk_vertex_ann"}
REGULAR_VERTEX = 0


def canonical_state(value):
    try:
        return float(value)
    except ValueError:
        return value


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


def read_local_edges(directory):
    path = directory / "surface-Q-edge-local.csv"
    header, rows = read_table(path)
    required = {
        "exc",
        "interface",
        "edge",
        "R (m)",
        "p_total",
        "p_in",
        "p_ann",
        "p_bulk_ann",
    }
    missing = required - set(header)
    if missing:
        raise RuntimeError(
            f"{path} is missing columns {sorted(missing)}; rerun with a Palace build "
            "which writes total and localized inside energy"
        )
    topology_id_columns = TOPOLOGY_ID_COLUMNS & set(header)
    if topology_id_columns and not TOPOLOGY_COLUMNS <= set(header):
        raise RuntimeError(
            f"{path} has incomplete automatic edge topology columns; missing "
            f"{sorted(TOPOLOGY_COLUMNS - set(header))}"
        )
    has_topology = TOPOLOGY_COLUMNS <= set(header)
    vertex_energy_columns = VERTEX_ENERGY_COLUMNS & set(header)
    if vertex_energy_columns and vertex_energy_columns != VERTEX_ENERGY_COLUMNS:
        raise RuntimeError(
            f"{path} has incomplete non-regular vertex energy columns; missing "
            f"{sorted(VERTEX_ENERGY_COLUMNS - vertex_energy_columns)}"
        )
    has_vertex_energies = vertex_energy_columns == VERTEX_ENERGY_COLUMNS

    state_column = header[0]
    values = {}
    excitations = defaultdict(set)
    for row in rows:
        state = canonical_state(row[state_column])
        excitation = int(row["exc"])
        interface = int(row["interface"])
        edge = int(row["edge"])
        radius = float(row["R (m)"])
        key = (state, excitation, interface, edge, radius)
        value = {
            "p_total": float(row["p_total"]),
            "p_in": float(row["p_in"]),
            "p_ann": float(row["p_ann"]),
            "p_bulk_ann": float(row["p_bulk_ann"]),
        }
        if has_vertex_energies:
            value.update(
                {
                    "p_vertex_in": float(row["p_vertex_in"]),
                    "p_bulk_vertex_ann": float(row["p_bulk_vertex_ann"]),
                }
            )
        if has_topology:
            value.update(
                {
                    "automatic": int(row["automatic"]),
                    "component": int(row["component"]),
                    "chain": int(row["chain"]),
                    "vertex_types": (int(row["v0_type"]), int(row["v1_type"])),
                    "edge_length": float(row["L_edge (m)"]),
                }
            )
        values[key] = value
        excitations[(state, interface)].add(excitation)
    return state_column, values, excitations


def read_surface_states(directory, excitations):
    path = directory / "surface-Q.csv"
    header, rows = read_table(path)
    state_column = header[0]
    values = {}
    for row in rows:
        state = canonical_state(row[state_column])
        quality_factors = {}
        for label in header:
            match = QUALITY_FACTOR_RE.fullmatch(label)
            if match:
                quality_factors[(int(match.group(1)), match.group(2))] = float(
                    row[label]
                )
        for label in header:
            match = PARTICIPATION_RE.fullmatch(label)
            if not match:
                continue
            interface = int(match.group(1))
            excitation_text = match.group(2)
            if excitation_text is None:
                candidates = excitations[(state, interface)]
                if len(candidates) != 1:
                    raise RuntimeError(
                        f"Cannot infer excitation for {state_column}={state}, "
                        f"interface={interface} in {path}"
                    )
                excitation = next(iter(candidates))
            else:
                excitation = int(excitation_text)
            quality_key = (interface, excitation_text)
            quality = quality_factors.get(quality_key, math.inf)
            values[(state, excitation, interface)] = {
                "p": float(row[label]),
                "Q": quality,
            }
    if not values:
        raise RuntimeError(f"No surface participation columns found in {path}")
    return state_column, values


def group_local_edges(values):
    groups = defaultdict(dict)
    for (state, excitation, interface, edge, radius), value in values.items():
        if value.get("automatic", 0) and value.get("chain", -1) >= 0:
            group = ("chain", value["component"], value["chain"])
        else:
            group = ("edge", edge)
        radius_values = groups[(state, excitation, interface, group)]
        endpoint_count = sum(
            vertex_type != REGULAR_VERTEX
            for vertex_type in value.get("vertex_types", ())
        )
        edge_length = value.get("edge_length", 0.0)
        endpoint_fraction = (
            min(1.0, radius * endpoint_count / edge_length)
            if edge_length > 0.0
            else 0.0
        )
        vertex_risk_p_inside = value.get(
            "p_vertex_in", endpoint_fraction * value["p_in"]
        )
        vertex_risk_p_bulk_ann = value.get(
            "p_bulk_vertex_ann", endpoint_fraction * value["p_bulk_ann"]
        )
        if radius not in radius_values:
            radius_values[radius] = {
                **value,
                "segment_count": 1,
                "nonregular_endpoint_count": endpoint_count,
                "vertex_risk_p_inside": vertex_risk_p_inside,
                "vertex_risk_p_bulk_ann": vertex_risk_p_bulk_ann,
            }
            continue
        pooled = radius_values[radius]
        for name, entry in value.items():
            if name.startswith("p_"):
                pooled[name] += entry
        pooled["segment_count"] += 1
        pooled["edge_length"] = pooled.get("edge_length", 0.0) + value.get(
            "edge_length", 0.0
        )
        pooled["nonregular_endpoint_count"] += endpoint_count
        pooled["vertex_risk_p_inside"] += vertex_risk_p_inside
        pooled["vertex_risk_p_bulk_ann"] += vertex_risk_p_bulk_ann
    return groups


def summarize_edge_topology(edge_values, radius):
    segment_count = 0
    automatic_chain_count = 0
    nonregular_chain_count = 0
    total_length = 0.0
    unsupported_length = 0.0
    total_inside = 0.0
    unsupported_inside = 0.0
    total_bulk = 0.0
    unsupported_bulk = 0.0
    for values in edge_values.values():
        value = values[radius]
        segment_count += value.get("segment_count", 1)
        if not value.get("automatic", 0):
            continue
        automatic_chain_count += 1
        length = value.get("edge_length", 0.0)
        endpoint_count = value.get("nonregular_endpoint_count", 0)
        total_length += length
        unsupported_length += min(length, radius * endpoint_count)
        total_inside += value["p_in"]
        unsupported_inside += value.get("vertex_risk_p_inside", 0.0)
        total_bulk += value["p_bulk_ann"]
        unsupported_bulk += value.get("vertex_risk_p_bulk_ann", 0.0)
        nonregular_chain_count += endpoint_count > 0
    unmodeled_vertex_length_fraction = (
        min(1.0, unsupported_length / total_length) if total_length > 0.0 else 0.0
    )
    unmodeled_vertex_surface_fraction = (
        min(1.0, unsupported_inside / total_inside) if total_inside > 0.0 else 0.0
    )
    unmodeled_vertex_bulk_fraction = (
        min(1.0, unsupported_bulk / total_bulk) if total_bulk > 0.0 else 0.0
    )
    return {
        "segment_count": segment_count,
        "automatic_chain_count": automatic_chain_count,
        "nonregular_chain_count": nonregular_chain_count,
        "unmodeled_vertex_length_fraction": unmodeled_vertex_length_fraction,
        "unmodeled_vertex_surface_fraction": unmodeled_vertex_surface_fraction,
        "unmodeled_vertex_bulk_fraction": unmodeled_vertex_bulk_fraction,
        "unmodeled_vertex_fraction": max(
            unmodeled_vertex_length_fraction,
            unmodeled_vertex_surface_fraction,
            unmodeled_vertex_bulk_fraction,
        ),
    }


def select_single_calibration_state(groups, directory):
    states = {(key[0], key[1]) for key in groups}
    if len(states) != 1:
        raise RuntimeError(
            f"Expected one state and excitation in calibration directory {directory}, "
            f"found {sorted(states)}"
        )
    return next(iter(states))


def fit_singular_amplitude(radius_values, fit_count):
    points = sorted(
        (radius, value["p_bulk_ann"] / radius)
        for radius, value in radius_values.items()
        if radius > 0.0 and value["p_bulk_ann"] > 0.0
    )
    if len(points) < fit_count:
        return 0.0, math.inf, []
    points = points[:fit_count]
    scale = points[-1][0]
    x = [radius / scale for radius, _ in points]
    y = [value for _, value in points]
    x_mean = statistics.fmean(x)
    y_mean = statistics.fmean(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0.0:
        return 0.0, math.inf, points
    slope = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x, y)
    ) / denominator
    intercept = y_mean - slope * x_mean
    residual = math.sqrt(
        sum((y_value - (intercept + slope * x_value)) ** 2
            for x_value, y_value in zip(x, y))
        / max(sum(y_value * y_value for y_value in y), 1.0e-300)
    )
    return max(0.0, intercept), residual, points


def fit_surface_annulus(radius_values, fit_count):
    points = sorted(
        (radius, value["p_ann"])
        for radius, value in radius_values.items()
        if radius > 0.0 and value["p_ann"] >= 0.0
    )
    if len(points) < fit_count:
        return 0.0, 0.0, math.inf, []
    points = points[:fit_count]
    scale = points[-1][0]
    x = [radius / scale for radius, _ in points]
    y = [value for _, value in points]
    x_mean = statistics.fmean(x)
    y_mean = statistics.fmean(y)
    denominator = sum((value - x_mean) ** 2 for value in x)
    if denominator == 0.0:
        return 0.0, 0.0, math.inf, points
    scaled_slope = sum(
        (x_value - x_mean) * (y_value - y_mean)
        for x_value, y_value in zip(x, y)
    ) / denominator
    intercept = y_mean - scaled_slope * x_mean
    residual = math.sqrt(
        sum(
            (y_value - (intercept + scaled_slope * x_value)) ** 2
            for x_value, y_value in zip(x, y)
        )
        / max(sum(y_value * y_value for y_value in y), 1.0e-300)
    )
    return max(0.0, intercept), scaled_slope / scale, residual, points


def edge_model_confidence(
    amplitude, bulk_amplitude, bulk_residual, surface_residual, residual_scale
):
    singular_fraction = min(
        1.0, max(0.0, amplitude) / max(bulk_amplitude, 1.0e-300)
    )
    residual = max(bulk_residual, surface_residual)
    if residual_scale == 0.0:
        fit_confidence = 1.0
    elif math.isfinite(residual):
        fit_confidence = 1.0 / (1.0 + (residual / residual_scale) ** 2)
    else:
        fit_confidence = 0.0
    return singular_fraction, fit_confidence


def nearest_radius(values, requested):
    if not values:
        raise RuntimeError("No matching edge-diagnostic radii are available")
    radius = min(values, key=lambda candidate: abs(candidate - requested))
    if not math.isclose(radius, requested, rel_tol=1.0e-8, abs_tol=1.0e-15):
        raise RuntimeError(
            f"No edge diagnostic at R={requested:.6e} m; nearest is {radius:.6e} m"
        )
    return radius


def write_rows(path, rows):
    if not rows:
        raise RuntimeError(f"No rows to write to {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def read_resolved_configs(directory):
    candidates = list(directory.glob("*resolved.json"))
    candidates.extend(directory.glob("palace*.json"))
    suffix = directory.name.removeprefix("iteration")
    if suffix.isdigit():
        candidates.extend(directory.parent.glob("*resolved.json"))
        candidates.extend(directory.parent.glob("palace*.json"))
    configs = {}
    for path in sorted(set(candidates)):
        try:
            with path.open() as stream:
                configs[path] = json.load(stream)
        except (OSError, json.JSONDecodeError):
            continue
    return configs


def read_fem_order(directory):
    orders = {}
    for path, config in read_resolved_configs(directory).items():
        try:
            orders[path] = int(config["Solver"]["Order"])
        except (KeyError, ValueError):
            continue
    unique_orders = set(orders.values())
    if len(unique_orders) > 1:
        details = ", ".join(f"{path.name}: p={order}" for path, order in orders.items())
        raise RuntimeError(f"Conflicting FEM orders in {directory}: {details}")
    return next(iter(unique_orders), None)


def read_interface_metadata(directory):
    metadata_by_path = {}
    for path, config in read_resolved_configs(directory).items():
        try:
            length_scale = float(config["Model"]["L0"])
            entries = config["Boundaries"]["Postprocessing"]["Dielectric"]
            metadata_by_path[path] = {
                int(entry["Index"]): {
                    "interface_type": str(entry["Type"]),
                    "interface_thickness_m": length_scale
                    * float(entry["Thickness"]),
                    "interface_permittivity": float(entry["Permittivity"]),
                    "flux_recovery": bool(entry.get("FluxRecovery", False)),
                }
                for entry in entries
            }
        except (KeyError, TypeError, ValueError):
            continue
    if not metadata_by_path:
        return {}
    metadata_sets = list(metadata_by_path.values())
    reference = metadata_sets[0]
    if any(metadata != reference for metadata in metadata_sets[1:]):
        paths = ", ".join(path.name for path in metadata_by_path)
        raise RuntimeError(
            f"Conflicting interface metadata in {directory}: {paths}"
        )
    return reference


def validate_interface_metadata(
    calibration_metadata, target_metadata, target_interface, calibration_interface
):
    if calibration_interface not in calibration_metadata:
        raise RuntimeError(
            f"Calibration table has no metadata for interface {calibration_interface}"
        )
    if target_interface not in target_metadata:
        raise RuntimeError(
            f"Resolved target config has no metadata for interface {target_interface}"
        )
    calibration = calibration_metadata[calibration_interface]
    target = target_metadata[target_interface]
    mismatches = []
    if target["interface_type"] != calibration["interface_type"]:
        mismatches.append(
            f"type {target['interface_type']} != {calibration['interface_type']}"
        )
    if target["flux_recovery"] != calibration["flux_recovery"]:
        mismatches.append(
            f"flux recovery {target['flux_recovery']} != "
            f"{calibration['flux_recovery']}"
        )
    for key, label in (
        ("interface_thickness_m", "thickness"),
        ("interface_permittivity", "permittivity"),
    ):
        if not math.isclose(
            target[key], calibration[key], rel_tol=1.0e-9, abs_tol=1.0e-18
        ):
            mismatches.append(f"{label} {target[key]} != {calibration[key]}")
    if mismatches:
        raise RuntimeError(
            f"Target interface {target_interface} does not match calibration "
            f"interface {calibration_interface}: {', '.join(mismatches)}"
        )


def relative_spread(values):
    if len(values) < 2:
        return math.nan
    return (max(values) - min(values)) / max(abs(values[-1]), 1.0e-300)


def calibrate(
    thin_directory,
    fabricated_directory,
    output,
    fit_count,
    include_history=True,
):
    thin_order = read_fem_order(thin_directory)
    fabricated_order = read_fem_order(fabricated_directory)
    if (
        thin_order is not None
        and fabricated_order is not None
        and thin_order != fabricated_order
    ):
        raise RuntimeError(
            f"Calibration FEM orders do not match: thin p={thin_order}, "
            f"fabricated p={fabricated_order}"
    )
    fem_order = thin_order if thin_order is not None else fabricated_order
    thin_metadata = read_interface_metadata(thin_directory)
    fabricated_metadata = read_interface_metadata(fabricated_directory)
    if not thin_metadata or not fabricated_metadata:
        raise RuntimeError(
            "Calibration output directories must contain resolved configs with "
            "interface dielectric metadata"
        )
    if thin_metadata != fabricated_metadata:
        raise RuntimeError(
            "Thin and fabricated calibration interface metadata do not match"
        )

    _, thin_local, thin_excitations = read_local_edges(thin_directory)
    _, fabricated_local, fabricated_excitations = read_local_edges(
        fabricated_directory
    )
    _, thin_surface = read_surface_states(thin_directory, thin_excitations)
    _, fabricated_surface = read_surface_states(
        fabricated_directory, fabricated_excitations
    )
    thin_groups = group_local_edges(thin_local)
    fabricated_groups = group_local_edges(fabricated_local)
    thin_state = select_single_calibration_state(thin_groups, thin_directory)
    fabricated_state = select_single_calibration_state(
        fabricated_groups, fabricated_directory
    )

    thin_interfaces = {key[2] for key in thin_groups if key[:2] == thin_state}
    fabricated_interfaces = {
        key[2] for key in fabricated_groups if key[:2] == fabricated_state
    }
    if thin_interfaces != fabricated_interfaces:
        raise RuntimeError(
            "Thin and fabricated calibration interface indices do not match"
        )

    rows = []
    for interface in sorted(thin_interfaces):
        thin_edges = {
            key[3]: values
            for key, values in thin_groups.items()
            if key[:3] == (*thin_state, interface)
        }
        fabricated_edges = {
            key[3]: values
            for key, values in fabricated_groups.items()
            if key[:3] == (*fabricated_state, interface)
        }
        bulk_fits = {
            edge: fit_singular_amplitude(values, fit_count)
            for edge, values in thin_edges.items()
        }
        surface_fits = {
            edge: fit_surface_annulus(values, fit_count)
            for edge, values in thin_edges.items()
        }
        amplitude = sum(value[0] for value in bulk_fits.values())
        if amplitude <= 0.0:
            raise RuntimeError(
                f"No positive singular amplitude for calibration interface {interface}"
            )
        bulk_fit_residual = sum(
            value[0] * value[1]
            for value in bulk_fits.values()
            if math.isfinite(value[1])
        ) / amplitude
        assigned_total = sum(
            next(iter(values.values()))["p_total"] for values in thin_edges.values()
        )
        surface_fit_residual = sum(
            next(iter(thin_edges[edge].values()))["p_total"] * value[2]
            for edge, value in surface_fits.items()
            if math.isfinite(value[2])
        ) / max(assigned_total, 1.0e-300)
        radii = sorted(set.intersection(*(set(values) for values in thin_edges.values())))
        fabricated_radii = set.intersection(
            *(set(values) for values in fabricated_edges.values())
        )
        radii = [radius for radius in radii if radius in fabricated_radii]
        for radius in radii:
            metadata = thin_metadata[interface]
            thin_inside = sum(
                values[radius]["p_in"] for values in thin_edges.values()
            )
            fabricated_inside = sum(
                values[radius]["p_in"] for values in fabricated_edges.values()
            )
            thin_outside = thin_surface[(*thin_state, interface)]["p"] - thin_inside
            fabricated_outside = (
                fabricated_surface[(*fabricated_state, interface)]["p"]
                - fabricated_inside
            )
            thin_smooth_inside = sum(
                max(0.0, surface_fits[edge][1] * radius)
                for edge in thin_edges
            )
            thin_bulk_amplitude = sum(
                values[radius]["p_bulk_ann"] / radius
                for values in thin_edges.values()
            )
            singular_fraction = min(
                1.0, amplitude / max(thin_bulk_amplitude, 1.0e-300)
            )
            rows.append(
                {
                    "interface": interface,
                    "radius_um": 1.0e6 * radius,
                    "fem_order": fem_order if fem_order is not None else "",
                    **metadata,
                    "inside_singular_coefficient_m": (
                        fabricated_inside - thin_smooth_inside
                    )
                    / amplitude,
                    "outside_delta_coefficient_m": (
                        fabricated_outside - thin_outside
                    )
                    / amplitude,
                    "fabricated_p_in": fabricated_inside,
                    "thin_smooth_p_in": thin_smooth_inside,
                    "singular_amplitude_per_m": amplitude,
                    "calibration_singular_fraction": singular_fraction,
                    "bulk_fit_relative_rms": bulk_fit_residual,
                    "surface_fit_relative_rms": surface_fit_residual,
                    "fit_count": fit_count,
                    "edge_count": sum(
                        value[0] > 0.0 for value in bulk_fits.values()
                    ),
                }
            )
    for row in rows:
        row["inside_coefficient_amr_spread"] = math.nan
        row["outside_delta_amr_spread"] = math.nan
    if include_history:
        thin_iterations = {
            int(path.name.removeprefix("iteration")): path
            for path in thin_directory.glob("iteration*")
            if path.is_dir() and path.name.removeprefix("iteration").isdigit()
        }
        fabricated_iterations = {
            int(path.name.removeprefix("iteration")): path
            for path in fabricated_directory.glob("iteration*")
            if path.is_dir() and path.name.removeprefix("iteration").isdigit()
        }
        common_iterations = sorted(thin_iterations.keys() & fabricated_iterations.keys())
        history = []
        for iteration in common_iterations[-(fit_count - 1) :]:
            history.append(
                calibrate(
                    thin_iterations[iteration],
                    fabricated_iterations[iteration],
                    None,
                    fit_count,
                    include_history=False,
                )
            )
        history_maps = [
            {(row["interface"], row["radius_um"]): row for row in history_rows}
            for history_rows in history
        ]
        for row in rows:
            key = (row["interface"], row["radius_um"])
            matching = [
                history_map[key] for history_map in history_maps if key in history_map
            ]
            matching.append(row)
            row["inside_coefficient_amr_spread"] = relative_spread(
                [entry["inside_singular_coefficient_m"] for entry in matching]
            )
            row["outside_delta_amr_spread"] = relative_spread(
                [entry["outside_delta_coefficient_m"] for entry in matching]
            )
    if output is not None:
        write_rows(output, rows)
    return rows


def read_calibration(path):
    _, rows = read_table(path)
    required = {
        "interface",
        "radius_um",
        "fem_order",
        "interface_type",
        "interface_thickness_m",
        "interface_permittivity",
        "flux_recovery",
        "inside_singular_coefficient_m",
        "outside_delta_coefficient_m",
        "fit_count",
        "inside_coefficient_amr_spread",
        "outside_delta_amr_spread",
    }
    missing = required - rows[0].keys()
    if missing:
        raise RuntimeError(f"{path} is missing calibration columns {sorted(missing)}")
    return [
        {
            **row,
            "interface": int(row["interface"]),
            "radius_um": float(row["radius_um"]),
            "fem_order": (
                int(row["fem_order"]) if row.get("fem_order", "") else None
            ),
            "interface_type": row["interface_type"],
            "interface_thickness_m": float(row["interface_thickness_m"]),
            "interface_permittivity": float(row["interface_permittivity"]),
            "flux_recovery": row["flux_recovery"].lower() == "true",
            "inside_singular_coefficient_m": float(
                row["inside_singular_coefficient_m"]
            ),
            "outside_delta_coefficient_m": float(
                row["outside_delta_coefficient_m"]
            ),
            "fit_count": int(row["fit_count"]),
            "inside_coefficient_amr_spread": float(
                row["inside_coefficient_amr_spread"]
            ),
            "outside_delta_amr_spread": float(
                row["outside_delta_amr_spread"]
            ),
        }
        for row in rows
    ]


def parse_interface_map(text):
    try:
        target, calibration = (int(value) for value in text.split("=", 1))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Expected TARGET=CALIBRATION interface indices, got {text!r}"
        ) from exc
    if target <= 0 or calibration <= 0:
        raise argparse.ArgumentTypeError("Interface indices must be positive")
    return target, calibration


def selected_radii(available, requested_um):
    if not requested_um:
        return sorted(available)
    result = []
    for radius_um in requested_um:
        result.append(nearest_radius(available, 1.0e-6 * radius_um))
    return sorted(set(result))


def calculate_radius_rows(
    calibration,
    input_directory,
    interface_maps,
    requested_radii_um,
    fit_count,
    fit_residual_scale,
    reference_directory,
):
    coefficients = {
        (row["interface"], 1.0e-6 * row["radius_um"]): row for row in calibration
    }
    calibration_interfaces = {row["interface"] for row in calibration}
    configured_map = dict(interface_maps)

    state_column, local, excitations = read_local_edges(input_directory)
    surface_state_column, surface = read_surface_states(input_directory, excitations)
    if surface_state_column != state_column:
        raise RuntimeError(
            f"State columns differ between surface tables: {state_column!r} and "
            f"{surface_state_column!r}"
        )
    reference = {}
    if reference_directory is not None:
        _, reference = read_surface_states(reference_directory, excitations)

    groups = group_local_edges(local)
    interface_groups = defaultdict(dict)
    for (state, excitation, interface, edge), values in groups.items():
        interface_groups[(state, excitation, interface)][edge] = values

    radius_rows = []
    for key, edge_values in sorted(interface_groups.items()):
        state, excitation, interface = key
        calibration_interface = configured_map.get(interface, interface)
        if calibration_interface not in calibration_interfaces:
            raise RuntimeError(
                f"No calibration mapping for target interface {interface}; use "
                f"--interface-map {interface}=CALIBRATION"
            )
        raw = surface[key]["p"]
        raw_quality = surface[key]["Q"]
        bulk_fits = {
            edge: fit_singular_amplitude(values, fit_count)
            for edge, values in edge_values.items()
        }
        surface_fits = {
            edge: fit_surface_annulus(values, fit_count)
            for edge, values in edge_values.items()
        }
        available = set.intersection(*(set(values) for values in edge_values.values()))
        available = {
            radius
            for radius in available
            if any(
                interface_key == calibration_interface
                and math.isclose(
                    coefficient_radius,
                    radius,
                    rel_tol=1.0e-8,
                    abs_tol=1.0e-15,
                )
                for interface_key, coefficient_radius in coefficients
            )
        }
        if not available:
            raise RuntimeError(
                f"No common diagnostic and calibration radii for target interface "
                f"{interface} mapped to calibration interface {calibration_interface}"
            )
        for radius in selected_radii(available, requested_radii_um):
            coefficient_key = min(
                (
                    candidate
                    for candidate in coefficients
                    if candidate[0] == calibration_interface
                ),
                key=lambda candidate: abs(candidate[1] - radius),
            )
            coefficient = coefficients[coefficient_key]
            corrected_inside = 0.0
            total_assigned = 0.0
            corrected_outside = 0.0
            smooth_inside_sum = 0.0
            modeled_singular_sum = 0.0
            outside_delta_sum = 0.0
            amplitude_sum = 0.0
            confidence_weighted_assigned_sum = 0.0
            bulk_amplitude_sum = 0.0
            weighted_bulk_residual = 0.0
            weighted_surface_residual = 0.0
            failed_bulk_fit = False
            failed_surface_fit = False
            active_edges = 0
            sampled_edges = 0
            for edge, values in edge_values.items():
                value = values[radius]
                assigned = value["p_total"]
                inside = value["p_in"]
                outside = assigned - inside
                if outside < -1.0e-12 * max(abs(assigned), abs(inside), 1.0e-300):
                    raise RuntimeError(
                        f"Local inside participation exceeds assigned participation "
                        f"for interface {interface}, edge {edge}, R={radius:.6e} m"
                    )
                outside = max(0.0, outside)
                bulk_amplitude = value["p_bulk_ann"] / radius
                amplitude, bulk_residual, bulk_points = bulk_fits[edge]
                _, smooth_density, surface_residual, surface_points = surface_fits[edge]
                singular_fraction, fit_confidence = edge_model_confidence(
                    amplitude,
                    bulk_amplitude,
                    bulk_residual,
                    surface_residual,
                    fit_residual_scale,
                )
                smooth_inside = max(0.0, smooth_density * radius)
                modeled_singular = (
                    coefficient["inside_singular_coefficient_m"] * amplitude
                )
                outside_delta = (
                    coefficient["outside_delta_coefficient_m"] * amplitude
                )
                total_assigned += assigned
                corrected_inside += smooth_inside + modeled_singular
                corrected_outside += outside + outside_delta
                smooth_inside_sum += smooth_inside
                modeled_singular_sum += modeled_singular
                outside_delta_sum += outside_delta
                amplitude_sum += amplitude
                confidence_weighted_assigned_sum += fit_confidence * assigned
                bulk_amplitude_sum += bulk_amplitude
                sampled_edges += bool(bulk_points and surface_points)
                if math.isfinite(bulk_residual):
                    weighted_bulk_residual += assigned * bulk_residual
                elif assigned > 0.0:
                    failed_bulk_fit = True
                if math.isfinite(surface_residual):
                    weighted_surface_residual += assigned * surface_residual
                elif assigned > 0.0:
                    failed_surface_fit = True
                if (
                    amplitude > 0.0
                    and math.isfinite(bulk_residual)
                    and math.isfinite(surface_residual)
                ):
                    active_edges += 1
            fit_confidence_summary = min(
                1.0,
                confidence_weighted_assigned_sum
                / max(total_assigned, 1.0e-300),
            )
            bulk_fit_summary = (
                math.inf
                if failed_bulk_fit
                else weighted_bulk_residual / total_assigned
                if total_assigned > 0.0
                else math.inf
            )
            surface_fit_summary = (
                math.inf
                if failed_surface_fit
                else weighted_surface_residual / total_assigned
                if total_assigned > 0.0
                else math.inf
            )
            if not math.isclose(
                total_assigned, raw, rel_tol=1.0e-8, abs_tol=1.0e-15
            ):
                raise RuntimeError(
                    f"Local assigned participations sum to {total_assigned:.16e}, "
                    f"but interface {interface} participation is {raw:.16e}"
                )
            corrected = corrected_outside + corrected_inside
            topology = summarize_edge_topology(edge_values, radius)
            singular_fraction = min(
                1.0,
                amplitude_sum / max(bulk_amplitude_sum, 1.0e-300),
            )
            singular_identifiability = singular_fraction * fit_confidence_summary
            term_scale = (
                abs(corrected_outside - outside_delta_sum)
                + abs(smooth_inside_sum)
                + abs(modeled_singular_sum)
                + abs(outside_delta_sum)
            )
            smooth_fraction = abs(smooth_inside_sum) / max(term_scale, 1.0e-300)
            modeled_fraction = (
                abs(modeled_singular_sum) + abs(outside_delta_sum)
            ) / max(term_scale, 1.0e-300)
            model_risk = (
                smooth_fraction * (1.0 - fit_confidence_summary)
                + modeled_fraction * (1.0 - singular_identifiability)
            )
            model_confidence = (1.0 - min(1.0, model_risk)) * (
                1.0 - topology["unmodeled_vertex_fraction"]
            )
            if raw > 0.0 and math.isfinite(raw_quality) and raw_quality > 0.0:
                loss_tangent = 1.0 / (raw * raw_quality)
                corrected_quality = (
                    math.inf if corrected <= 0.0 else 1.0 / (corrected * loss_tangent)
                )
            else:
                corrected_quality = math.inf
            reference_value = reference.get(key, {}).get("p", math.nan)
            radius_rows.append(
                {
                    "method": "radius",
                    "state_column": state_column,
                    "state": state,
                    "exc": excitation,
                    "interface": interface,
                    "calibration_interface": calibration_interface,
                    "radius_um": 1.0e6 * radius,
                    "p_raw": raw,
                    "p_corrected": corrected,
                    "Q_corrected": corrected_quality,
                    "p_out": corrected_outside - outside_delta_sum,
                    "p_smooth_in": smooth_inside_sum,
                    "p_modeled_singular": modeled_singular_sum,
                    "p_outside_delta": outside_delta_sum,
                    "singular_fraction": singular_fraction,
                    "singular_identifiability": singular_identifiability,
                    "fit_confidence": fit_confidence_summary,
                    "modeled_fraction": modeled_fraction,
                    "model_confidence": model_confidence,
                    "fit_relative_rms": max(
                        bulk_fit_summary, surface_fit_summary
                    ),
                    "bulk_fit_relative_rms": bulk_fit_summary,
                    "surface_fit_relative_rms": surface_fit_summary,
                    "active_edges": active_edges,
                    "sampled_edges": sampled_edges,
                    "edge_count": len(edge_values),
                    **topology,
                    "sampled_edge_fraction": sampled_edges / len(edge_values),
                    "calibration_inside_amr_spread": coefficient[
                        "inside_coefficient_amr_spread"
                    ],
                    "calibration_outside_amr_spread": coefficient[
                        "outside_delta_amr_spread"
                    ],
                    "selected": False,
                    "window_min_um": math.nan,
                    "window_max_um": math.nan,
                    "radius_spread": math.nan,
                    "amr_spread": math.nan,
                    "relative_uncertainty": math.nan,
                    "p_reference": reference_value,
                    "relative_error": (
                        (corrected - reference_value) / reference_value
                        if reference_value > 0.0
                        else math.nan
                    ),
                }
            )
    return radius_rows


def summarize_window(candidates, amr_spread=math.nan):
    candidates = sorted(candidates, key=lambda row: row["radius_um"])
    corrected_values = [row["p_corrected"] for row in candidates]
    radii = [row["radius_um"] for row in candidates]
    corrected = statistics.fmean(corrected_values)
    method = "window-mean" if len(corrected_values) > 1 else "radius"
    template = candidates[0]
    reference_value = template["p_reference"]
    quality = template["Q_corrected"]
    if math.isfinite(quality) and corrected > 0.0:
        quality *= template["p_corrected"] / corrected
    radius_spread = (
        (max(corrected_values) - min(corrected_values))
        / max(abs(corrected), 1.0e-300)
    )
    uncertainty_components = [
        0.5 * radius_spread,
        amr_spread,
        max(row["calibration_inside_amr_spread"] for row in candidates),
        max(row["calibration_outside_amr_spread"] for row in candidates),
        statistics.median(row["fit_relative_rms"] for row in candidates),
        max(1.0 - row["model_confidence"] for row in candidates),
    ]
    relative_uncertainty = max(
        value for value in uncertainty_components if math.isfinite(value)
    )
    return {
        **template,
        "method": method,
        "radius_um": math.nan,
        "p_corrected": corrected,
        "Q_corrected": quality,
        "p_out": statistics.fmean(row["p_out"] for row in candidates),
        "p_smooth_in": statistics.fmean(
            row["p_smooth_in"] for row in candidates
        ),
        "p_modeled_singular": statistics.fmean(
            row["p_modeled_singular"] for row in candidates
        ),
        "p_outside_delta": statistics.fmean(
            row["p_outside_delta"] for row in candidates
        ),
        "singular_fraction": statistics.median(
            row["singular_fraction"] for row in candidates
        ),
        "singular_identifiability": statistics.median(
            row["singular_identifiability"] for row in candidates
        ),
        "fit_confidence": statistics.median(
            row["fit_confidence"] for row in candidates
        ),
        "modeled_fraction": statistics.median(
            row["modeled_fraction"] for row in candidates
        ),
        "model_confidence": statistics.median(
            row["model_confidence"] for row in candidates
        ),
        "unmodeled_vertex_fraction": max(
            row.get("unmodeled_vertex_fraction", 0.0) for row in candidates
        ),
        "unmodeled_vertex_length_fraction": max(
            row.get("unmodeled_vertex_length_fraction", 0.0)
            for row in candidates
        ),
        "unmodeled_vertex_surface_fraction": max(
            row.get("unmodeled_vertex_surface_fraction", 0.0)
            for row in candidates
        ),
        "unmodeled_vertex_bulk_fraction": max(
            row.get("unmodeled_vertex_bulk_fraction", 0.0)
            for row in candidates
        ),
        "fit_relative_rms": statistics.median(
            row["fit_relative_rms"] for row in candidates
        ),
        "bulk_fit_relative_rms": statistics.median(
            row["bulk_fit_relative_rms"] for row in candidates
        ),
        "surface_fit_relative_rms": statistics.median(
            row["surface_fit_relative_rms"] for row in candidates
        ),
        "sampled_edges": min(row["sampled_edges"] for row in candidates),
        "sampled_edge_fraction": min(
            row["sampled_edge_fraction"] for row in candidates
        ),
        "calibration_inside_amr_spread": max(
            row["calibration_inside_amr_spread"] for row in candidates
        ),
        "calibration_outside_amr_spread": max(
            row["calibration_outside_amr_spread"] for row in candidates
        ),
        "selected": True,
        "window_min_um": min(radii),
        "window_max_um": max(radii),
        "radius_spread": radius_spread,
        "amr_spread": amr_spread,
        "relative_uncertainty": relative_uncertainty,
        "relative_error": (
            (corrected - reference_value) / reference_value
            if reference_value > 0.0
            else math.nan
        ),
    }


def group_radius_rows(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["state"], row["exc"], row["interface"])].append(row)
    for candidates in groups.values():
        candidates.sort(key=lambda row: row["radius_um"])
    return groups


def matching_window(rows, radii):
    by_radius = {row["radius_um"]: row for row in rows}
    return [by_radius[radius] for radius in radii]


def choose_window(evaluated, tolerance):
    stable = [
        entry
        for entry in evaluated
        if math.isfinite(entry[1]) and entry[1] <= tolerance
    ]
    if stable:
        return stable[0]
    if any(math.isfinite(entry[1]) for entry in evaluated):
        return min(
            evaluated,
            key=lambda entry: (
                entry[1] if math.isfinite(entry[1]) else math.inf,
                entry[0][0]["radius_um"],
            ),
        )
    return evaluated[0]


def apply(
    calibration_path,
    input_directory,
    output,
    interface_maps,
    requested_radii_um,
    fit_count,
    fit_residual_scale,
    reference_directory,
    amr_samples=3,
    amr_tolerance=0.02,
    allow_order_mismatch=False,
):
    calibration = read_calibration(calibration_path)
    calibration_orders = {
        row["fem_order"] for row in calibration if row["fem_order"] is not None
    }
    if len(calibration_orders) > 1:
        raise RuntimeError(
            f"Calibration table contains multiple FEM orders: "
            f"{sorted(calibration_orders)}"
        )
    calibration_fit_counts = {row["fit_count"] for row in calibration}
    if len(calibration_fit_counts) != 1:
        raise RuntimeError(
            f"Calibration table contains multiple fit counts: "
            f"{sorted(calibration_fit_counts)}"
        )
    calibration_fit_count = next(iter(calibration_fit_counts))
    if fit_count is None:
        fit_count = calibration_fit_count
    elif fit_count != calibration_fit_count:
        raise RuntimeError(
            f"Calibration uses {calibration_fit_count} fit radii, but application "
            f"requested {fit_count}; recalibrate with the desired fit count"
        )
    calibration_metadata = {}
    for row in calibration:
        metadata = {
            key: row[key]
            for key in (
                "interface_type",
                "interface_thickness_m",
                "interface_permittivity",
                "flux_recovery",
            )
        }
        interface = row["interface"]
        if interface in calibration_metadata and calibration_metadata[interface] != metadata:
            raise RuntimeError(
                f"Calibration table contains conflicting metadata for interface "
                f"{interface}"
            )
        calibration_metadata[interface] = metadata
    target_order = read_fem_order(input_directory)
    if calibration_orders and target_order is not None:
        calibration_order = next(iter(calibration_orders))
        if calibration_order != target_order and not allow_order_mismatch:
            raise RuntimeError(
                f"Calibration uses FEM order p={calibration_order}, but target uses "
                f"p={target_order}; use a matched calibration or explicitly pass "
                "--allow-order-mismatch"
            )
    target_metadata = read_interface_metadata(input_directory)
    if not target_metadata:
        raise RuntimeError(
            "Target output directory must contain a resolved config with interface "
            "dielectric metadata"
        )
    configured_map = dict(interface_maps)
    for target_interface in target_metadata:
        calibration_interface = configured_map.get(target_interface, target_interface)
        if calibration_interface in calibration_metadata:
            validate_interface_metadata(
                calibration_metadata,
                target_metadata,
                target_interface,
                calibration_interface,
            )
    radius_rows = calculate_radius_rows(
        calibration,
        input_directory,
        interface_maps,
        requested_radii_um,
        fit_count,
        fit_residual_scale,
        reference_directory,
    )
    final_groups = group_radius_rows(radius_rows)

    history_groups = []
    if amr_samples > 1:
        iteration_directories = []
        for candidate in input_directory.glob("iteration*"):
            suffix = candidate.name.removeprefix("iteration")
            if suffix.isdigit():
                iteration_directories.append((int(suffix), candidate))
        iteration_directories.sort()
        for _, directory in iteration_directories[-(amr_samples - 1) :]:
            try:
                history_rows = calculate_radius_rows(
                    calibration,
                    directory,
                    interface_maps,
                    requested_radii_um,
                    fit_count,
                    fit_residual_scale,
                    None,
                )
            except RuntimeError:
                history_groups = []
                break
            history_groups.append(group_radius_rows(history_rows))
    history_groups.append(final_groups)

    rows = list(radius_rows)
    for key, candidates in sorted(final_groups.items()):
        if requested_radii_um:
            windows = [candidates]
        elif len(candidates) >= fit_count:
            windows = [
                candidates[start : start + fit_count]
                for start in range(len(candidates) - fit_count + 1)
            ]
        else:
            windows = [candidates]

        evaluated = []
        for window in windows:
            radii = [row["radius_um"] for row in window]
            history_values = []
            for groups in history_groups:
                if key not in groups:
                    continue
                try:
                    history_window = matching_window(groups[key], radii)
                except KeyError:
                    continue
                history_values.append(
                    summarize_window(history_window)["p_corrected"]
                )
            if len(history_values) > 1:
                amr_spread = (
                    max(history_values) - min(history_values)
                ) / max(abs(history_values[-1]), 1.0e-300)
            else:
                amr_spread = math.nan
            evaluated.append((window, amr_spread))

        selected, amr_spread = choose_window(evaluated, amr_tolerance)
        for row in selected:
            row["selected"] = True
        rows.append(summarize_window(selected, amr_spread))
    write_rows(output, rows)
    return rows


def print_calibration(rows):
    print(
        "interface  p  R (um)  C_inside (m)  C_outside (m)  singular fraction  "
        "bulk RMS  surface RMS  inside AMR spread  outside AMR spread"
    )
    for row in rows:
        print(
            f"{row['interface']:9d} {str(row['fem_order']):>2} "
            f"{row['radius_um']:7.2f} "
            f"{row['inside_singular_coefficient_m']:12.4e} "
            f"{row['outside_delta_coefficient_m']:13.4e} "
            f"{row['calibration_singular_fraction']:17.3f} "
            f"{row['bulk_fit_relative_rms']:8.3f} "
            f"{row['surface_fit_relative_rms']:11.3f} "
            f"{row['inside_coefficient_amr_spread']:17.2%} "
            f"{row['outside_delta_amr_spread']:18.2%}"
        )


def print_application(rows, fit_residual_scale, amr_tolerance):
    summaries = [row for row in rows if row["method"] != "radius"]
    print(
        "state  exc  interface  p_raw       p_corrected  model confidence  "
        "fit RMS  sampled edges  R window (um)  uncertainty  error"
    )
    for row in summaries:
        error = (
            f"{row['relative_error']:+7.2%}"
            if math.isfinite(row["relative_error"])
            else "    n/a"
        )
        print(
            f"{row['state']:>5} {row['exc']:4d} {row['interface']:10d} "
            f"{row['p_raw']:11.4e} {row['p_corrected']:11.4e} "
            f"{row['model_confidence']:16.3f} "
            f"{row['fit_relative_rms']:7.3f} "
            f"{row['sampled_edges']:5d}/{row['edge_count']:<5d} "
            f"{row['window_min_um']:4.2f}-{row['window_max_um']:4.2f} "
            f"{row['relative_uncertainty']:11.2%} "
            f"{error}"
        )
        warnings = []
        if (
            fit_residual_scale > 0.0
            and row["fit_relative_rms"] > fit_residual_scale
        ):
            warnings.append("poor local edge fit")
        if row["radius_spread"] > 0.1:
            warnings.append("large radius-window spread")
        if row["model_confidence"] < 0.9:
            warnings.append("low correction confidence")
        if row.get("unmodeled_vertex_fraction", 0.0) > 0.0:
            warnings.append("corner/endpoint neighborhoods use a straight-edge model")
        if row["calibration_inside_amr_spread"] > 0.05:
            warnings.append("calibration inside coefficient is AMR-unstable")
        if row["calibration_outside_amr_spread"] > 0.05:
            warnings.append("calibration outside delta is AMR-unstable")
        if math.isfinite(row["amr_spread"]):
            if row["amr_spread"] > amr_tolerance:
                warnings.append("AMR-unstable")
        else:
            warnings.append("AMR convergence unchecked")
        if row["relative_uncertainty"] > 0.1:
            warnings.append("large estimated uncertainty")
        if warnings:
            print(f"  warning: interface {row['interface']}: {', '.join(warnings)}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate and apply a local fabrication-aware correction to Palace "
            "thin-metal surface participations."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    calibrate_parser = subparsers.add_parser(
        "calibrate", help="Create fabrication-process coefficients from 2D results."
    )
    calibrate_parser.add_argument("--thin", type=Path, required=True)
    calibrate_parser.add_argument("--fabricated", type=Path, required=True)
    calibrate_parser.add_argument("--output", type=Path, required=True)
    calibrate_parser.add_argument("--fit-count", type=int, default=3)

    apply_parser = subparsers.add_parser(
        "apply", help="Apply process coefficients to a localized thin-metal result."
    )
    apply_parser.add_argument("--calibration", type=Path, required=True)
    apply_parser.add_argument("--input", type=Path, required=True)
    apply_parser.add_argument("--output", type=Path, required=True)
    apply_parser.add_argument("--reference", type=Path)
    apply_parser.add_argument(
        "--interface-map",
        action="append",
        type=parse_interface_map,
        default=[],
        metavar="TARGET=CALIBRATION",
    )
    apply_parser.add_argument(
        "--radius",
        action="append",
        type=float,
        default=[],
        metavar="R_UM",
        help="Use this radius in the window mean; repeat to select multiple radii.",
    )
    apply_parser.add_argument(
        "--fit-count",
        type=int,
        help="Number of smallest radii used by the edge fit; defaults to and must "
        "match the value recorded during calibration.",
    )
    apply_parser.add_argument(
        "--fit-residual-scale",
        type=float,
        default=0.1,
        help="Relative RMS at which edge-model confidence is reduced by one half.",
    )
    apply_parser.add_argument(
        "--amr-samples",
        type=int,
        default=3,
        help="Number of final AMR states used to select a stable radius window.",
    )
    apply_parser.add_argument(
        "--amr-tolerance",
        type=float,
        default=0.02,
        help="Maximum relative AMR spread accepted for the smallest radius window.",
    )
    apply_parser.add_argument(
        "--allow-order-mismatch",
        action="store_true",
        help="Apply despite a detected mismatch between calibration and target order.",
    )

    args = parser.parse_args()
    if args.command == "calibrate":
        if args.fit_count < 3:
            parser.error("--fit-count must be at least 3")
        rows = calibrate(args.thin, args.fabricated, args.output, args.fit_count)
        print_calibration(rows)
    else:
        if args.fit_count is not None and args.fit_count < 3:
            parser.error("--fit-count must be at least 3")
        if args.fit_residual_scale < 0.0:
            parser.error("--fit-residual-scale must be nonnegative")
        if args.amr_samples < 1:
            parser.error("--amr-samples must be positive")
        if args.amr_tolerance < 0.0:
            parser.error("--amr-tolerance must be nonnegative")
        rows = apply(
            args.calibration,
            args.input,
            args.output,
            args.interface_map,
            args.radius,
            args.fit_count,
            args.fit_residual_scale,
            args.reference,
            args.amr_samples,
            args.amr_tolerance,
            args.allow_order_mismatch,
        )
        print_application(rows, args.fit_residual_scale, args.amr_tolerance)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
