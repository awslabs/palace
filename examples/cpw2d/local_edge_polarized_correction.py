#!/usr/bin/env python3

import argparse
import math
import statistics
from collections import defaultdict
from pathlib import Path

import local_edge_correction as common


CHANNELS = (
    "top_n",
    "top_m",
    "top_l",
    "bottom_n",
    "bottom_m",
    "bottom_l",
)
SURFACE_COMPONENTS = ("n", "t")
DESCRIPTOR_CHANNELS = {
    ("MA", "n"): ("top_n",),
    ("MS", "n"): ("bottom_n",),
    ("SA", "n"): ("top_n",),
    ("SA", "t"): ("top_m",),
}


def active_components(interface_type):
    return tuple(
        component
        for component in SURFACE_COMPONENTS
        if (interface_type, component) in DESCRIPTOR_CHANNELS
    )


def read_local_edges(directory):
    path = directory / "surface-Q-edge-local.csv"
    header, rows = common.read_table(path)
    required = {
        "exc",
        "interface",
        "edge",
        "R (m)",
        "p_total",
        "p_in",
        "p_ann",
        "p_bulk_ann",
        "p_total_n",
        "p_total_t",
        "p_in_n",
        "p_in_t",
        "p_ann_n",
        "p_ann_t",
        *(f"p_bulk_{channel}_ann" for channel in CHANNELS),
    }
    missing = required - set(header)
    if missing:
        raise RuntimeError(
            f"{path} is missing polarized columns {sorted(missing)}; rerun with "
            "EdgeFrameNormal and a Palace build with polarized edge diagnostics"
        )
    topology_id_columns = common.TOPOLOGY_ID_COLUMNS & set(header)
    if topology_id_columns and not common.TOPOLOGY_COLUMNS <= set(header):
        raise RuntimeError(
            f"{path} has incomplete automatic edge topology columns; missing "
            f"{sorted(common.TOPOLOGY_COLUMNS - set(header))}"
        )
    has_topology = common.TOPOLOGY_COLUMNS <= set(header)
    vertex_energy_columns = common.VERTEX_ENERGY_COLUMNS & set(header)
    if vertex_energy_columns and vertex_energy_columns != common.VERTEX_ENERGY_COLUMNS:
        raise RuntimeError(
            f"{path} has incomplete non-regular vertex energy columns; missing "
            f"{sorted(common.VERTEX_ENERGY_COLUMNS - vertex_energy_columns)}"
        )
    has_vertex_energies = vertex_energy_columns == common.VERTEX_ENERGY_COLUMNS

    state_column = header[0]
    values = {}
    excitations = defaultdict(set)
    for row in rows:
        state = common.canonical_state(row[state_column])
        excitation = int(row["exc"])
        interface = int(row["interface"])
        edge = int(row["edge"])
        radius = float(row["R (m)"])
        value = {
            name: float(row[name])
            for name in (
                "p_total",
                "p_in",
                "p_ann",
                "p_bulk_ann",
                "p_total_n",
                "p_total_t",
                "p_in_n",
                "p_in_t",
                "p_ann_n",
                "p_ann_t",
                *(f"p_bulk_{channel}_ann" for channel in CHANNELS),
            )
        }
        if has_vertex_energies:
            value.update(
                {
                    "p_vertex_in": float(row["p_vertex_in"]),
                    "p_bulk_vertex_ann": float(row["p_bulk_vertex_ann"]),
                }
            )
        for total, parts in (
            ("p_total", ("p_total_n", "p_total_t")),
            ("p_in", ("p_in_n", "p_in_t")),
            ("p_ann", ("p_ann_n", "p_ann_t")),
            (
                "p_bulk_ann",
                tuple(f"p_bulk_{channel}_ann" for channel in CHANNELS),
            ),
        ):
            part_sum = sum(value[part] for part in parts)
            if not math.isclose(
                value[total], part_sum, rel_tol=5.0e-8, abs_tol=1.0e-18
            ):
                raise RuntimeError(
                    f"Polarized columns in {path} do not conserve {total} for "
                    f"interface {interface}, edge {edge}, R={radius:.6e} m: "
                    f"{part_sum:.16e} != {value[total]:.16e}"
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
        key = (state, excitation, interface, edge, radius)
        values[key] = value
        excitations[(state, interface)].add(excitation)
    return state_column, values, excitations


def fit_bulk_channel(radius_values, channel, fit_count):
    values = {
        radius: {"p_bulk_ann": value[f"p_bulk_{channel}_ann"]}
        for radius, value in radius_values.items()
    }
    return common.fit_singular_amplitude(values, fit_count)


def fit_descriptor(radius_values, interface_type, component, fit_count):
    channels = DESCRIPTOR_CHANNELS[(interface_type, component)]
    fits = {
        channel: fit_bulk_channel(radius_values, channel, fit_count)
        for channel in channels
    }
    amplitude = sum(fit[0] for fit in fits.values())
    if amplitude <= 0.0:
        return 0.0, math.inf, fits
    residual = sum(
        fit[0] * fit[1]
        for fit in fits.values()
        if math.isfinite(fit[1])
    ) / amplitude
    if any(fit[0] > 0.0 and not math.isfinite(fit[1]) for fit in fits.values()):
        residual = math.inf
    return amplitude, residual, fits


def fit_surface_component(radius_values, component, fit_count):
    values = {
        radius: {"p_ann": value[f"p_ann_{component}"]}
        for radius, value in radius_values.items()
    }
    return common.fit_surface_annulus(values, fit_count)


def check_assigned_total(local_edges, surface, state, interface, label):
    assigned = sum(
        next(iter(values.values()))["p_total"]
        for key, values in local_edges.items()
        if key[:3] == (*state, interface)
    )
    expected = surface[(*state, interface)]["p"]
    if not math.isclose(assigned, expected, rel_tol=1.0e-8, abs_tol=1.0e-15):
        raise RuntimeError(
            f"{label} local assigned participation for interface {interface} is "
            f"{assigned:.16e}, but surface-Q.csv reports {expected:.16e}"
        )


def calibrate(thin_directory, fabricated_directory, output, fit_count):
    thin_order = common.read_fem_order(thin_directory)
    fabricated_order = common.read_fem_order(fabricated_directory)
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
    thin_metadata = common.read_interface_metadata(thin_directory)
    fabricated_metadata = common.read_interface_metadata(fabricated_directory)
    if not thin_metadata or not fabricated_metadata:
        raise RuntimeError(
            "Calibration directories must contain resolved configs with interface "
            "dielectric metadata"
        )
    if thin_metadata != fabricated_metadata:
        raise RuntimeError(
            "Thin and fabricated calibration interface metadata do not match"
        )

    _, thin_local, thin_excitations = read_local_edges(thin_directory)
    _, fabricated_local, fabricated_excitations = read_local_edges(
        fabricated_directory
    )
    _, thin_surface = common.read_surface_states(thin_directory, thin_excitations)
    _, fabricated_surface = common.read_surface_states(
        fabricated_directory, fabricated_excitations
    )
    thin_groups = common.group_local_edges(thin_local)
    fabricated_groups = common.group_local_edges(fabricated_local)
    thin_state = common.select_single_calibration_state(thin_groups, thin_directory)
    fabricated_state = common.select_single_calibration_state(
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
        interface_type = thin_metadata[interface]["interface_type"]
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
        check_assigned_total(
            thin_groups, thin_surface, thin_state, interface, "Thin calibration"
        )
        check_assigned_total(
            fabricated_groups,
            fabricated_surface,
            fabricated_state,
            interface,
            "Fabricated calibration",
        )

        descriptor_fits = {
            (edge, component): fit_descriptor(
                values, interface_type, component, fit_count
            )
            for edge, values in thin_edges.items()
            for component in active_components(interface_type)
        }
        surface_fits = {
            (edge, component): fit_surface_component(
                values, component, fit_count
            )
            for edge, values in thin_edges.items()
            for component in active_components(interface_type)
        }
        channel_amplitudes = {
            channel: sum(
                fit_bulk_channel(values, channel, fit_count)[0]
                for values in thin_edges.values()
            )
            for channel in CHANNELS
        }
        radii = sorted(
            set.intersection(*(set(values) for values in thin_edges.values()))
            & set.intersection(*(set(values) for values in fabricated_edges.values()))
        )
        for component in active_components(interface_type):
            amplitude = sum(
                descriptor_fits[(edge, component)][0] for edge in thin_edges
            )
            if amplitude <= 0.0:
                raise RuntimeError(
                    f"No positive {component}-descriptor amplitude for calibration "
                    f"interface {interface} ({interface_type})"
                )
            bulk_residual = sum(
                descriptor_fits[(edge, component)][0]
                * descriptor_fits[(edge, component)][1]
                for edge in thin_edges
                if math.isfinite(descriptor_fits[(edge, component)][1])
            ) / amplitude
            assigned_component = sum(
                next(iter(values.values()))[f"p_total_{component}"]
                for values in thin_edges.values()
            )
            surface_residual = sum(
                next(iter(thin_edges[edge].values()))[f"p_total_{component}"]
                * surface_fits[(edge, component)][2]
                for edge in thin_edges
                if math.isfinite(surface_fits[(edge, component)][2])
            ) / max(assigned_component, 1.0e-300)
            for radius in radii:
                thin_inside = sum(
                    values[radius][f"p_in_{component}"]
                    for values in thin_edges.values()
                )
                fabricated_inside = sum(
                    values[radius][f"p_in_{component}"]
                    for values in fabricated_edges.values()
                )
                thin_total = sum(
                    values[radius][f"p_total_{component}"]
                    for values in thin_edges.values()
                )
                fabricated_total = sum(
                    values[radius][f"p_total_{component}"]
                    for values in fabricated_edges.values()
                )
                thin_smooth_inside = sum(
                    max(0.0, surface_fits[(edge, component)][1] * radius)
                    for edge in thin_edges
                )
                rows.append(
                    {
                        "interface": interface,
                        "surface_component": component,
                        "descriptor_channels": "+".join(
                            DESCRIPTOR_CHANNELS[(interface_type, component)]
                        ),
                        "radius_um": 1.0e6 * radius,
                        "fem_order": fem_order if fem_order is not None else "",
                        **thin_metadata[interface],
                        "inside_singular_coefficient_m": (
                            fabricated_inside - thin_smooth_inside
                        )
                        / amplitude,
                        "outside_delta_coefficient_m": (
                            (fabricated_total - fabricated_inside)
                            - (thin_total - thin_inside)
                        )
                        / amplitude,
                        "fabricated_p_in_component": fabricated_inside,
                        "thin_smooth_p_in_component": thin_smooth_inside,
                        "descriptor_amplitude_per_m": amplitude,
                        **{
                            f"calibration_{channel}_per_m": value
                            for channel, value in channel_amplitudes.items()
                        },
                        "bulk_fit_relative_rms": bulk_residual,
                        "surface_fit_relative_rms": surface_residual,
                        "fit_count": fit_count,
                    }
                )
    common.write_rows(output, rows)
    return rows


def read_calibration(path):
    _, rows = common.read_table(path)
    required = {
        "interface",
        "surface_component",
        "descriptor_channels",
        "radius_um",
        "fem_order",
        "interface_type",
        "interface_thickness_m",
        "interface_permittivity",
        "flux_recovery",
        "inside_singular_coefficient_m",
        "outside_delta_coefficient_m",
        "fit_count",
        *(f"calibration_{channel}_per_m" for channel in CHANNELS),
    }
    missing = required - rows[0].keys()
    if missing:
        raise RuntimeError(
            f"{path} is missing polarized calibration columns {sorted(missing)}"
        )
    return [
        {
            **row,
            "interface": int(row["interface"]),
            "radius_um": float(row["radius_um"]),
            "fem_order": (
                int(row["fem_order"]) if row.get("fem_order", "") else None
            ),
            "interface_thickness_m": float(row["interface_thickness_m"]),
            "interface_permittivity": float(row["interface_permittivity"]),
            "flux_recovery": row["flux_recovery"].lower() == "true",
            "inside_singular_coefficient_m": float(
                row["inside_singular_coefficient_m"]
            ),
            "outside_delta_coefficient_m": float(
                row["outside_delta_coefficient_m"]
            ),
            **{
                f"calibration_{channel}_per_m": float(
                    row[f"calibration_{channel}_per_m"]
                )
                for channel in CHANNELS
            },
            "fit_count": int(row["fit_count"]),
        }
        for row in rows
    ]


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
        (
            row["interface"],
            1.0e-6 * row["radius_um"],
            row["surface_component"],
        ): row
        for row in calibration
    }
    calibration_interfaces = {row["interface"] for row in calibration}
    configured_map = dict(interface_maps)
    target_metadata = common.read_interface_metadata(input_directory)

    state_column, local, excitations = read_local_edges(input_directory)
    surface_state_column, surface = common.read_surface_states(
        input_directory, excitations
    )
    if surface_state_column != state_column:
        raise RuntimeError(
            f"State columns differ between surface tables: {state_column!r} and "
            f"{surface_state_column!r}"
        )
    reference = {}
    if reference_directory is not None:
        _, reference = common.read_surface_states(reference_directory, excitations)

    groups = common.group_local_edges(local)
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
        calibration_row = next(
            row
            for row in calibration
            if row["interface"] == calibration_interface
        )
        calibration_channel_amplitudes = {
            channel: calibration_row[f"calibration_{channel}_per_m"]
            for channel in CHANNELS
        }
        interface_type = target_metadata[interface]["interface_type"]
        components = active_components(interface_type)
        raw = surface[key]["p"]
        raw_quality = surface[key]["Q"]
        descriptor_fits = {
            (edge, component): fit_descriptor(
                values, interface_type, component, fit_count
            )
            for edge, values in edge_values.items()
            for component in components
        }
        surface_fits = {
            (edge, component): fit_surface_component(
                values, component, fit_count
            )
            for edge, values in edge_values.items()
            for component in components
        }
        channel_fits = {
            (edge, channel): fit_bulk_channel(values, channel, fit_count)
            for edge, values in edge_values.items()
            for channel in CHANNELS
        }
        available = set.intersection(*(set(values) for values in edge_values.values()))
        available = {
            radius
            for radius in available
            if all(
                any(
                    candidate_interface == calibration_interface
                    and candidate_component == component
                    and math.isclose(
                        candidate_radius,
                        radius,
                        rel_tol=1.0e-8,
                        abs_tol=1.0e-15,
                    )
                    for (
                        candidate_interface,
                        candidate_radius,
                        candidate_component,
                    ) in coefficients
                )
                for component in components
            )
        }
        if not available:
            raise RuntimeError(
                f"No common diagnostic and polarized calibration radii for target "
                f"interface {interface}"
            )

        for radius in common.selected_radii(available, requested_radii_um):
            result = {
                "method": "radius",
                "state_column": state_column,
                "state": state,
                "exc": excitation,
                "interface": interface,
                "calibration_interface": calibration_interface,
                "interface_type": interface_type,
                "radius_um": 1.0e6 * radius,
                "p_raw": raw,
            }
            total_assigned = 0.0
            corrected = 0.0
            corrected_outside = 0.0
            smooth_inside_sum = 0.0
            modeled_singular_sum = 0.0
            outside_delta_sum = 0.0
            amplitude_total = 0.0
            bulk_amplitude_total = 0.0
            confidence_weighted_assigned_sum = 0.0
            weighted_bulk_residual = 0.0
            weighted_surface_residual = 0.0
            failed_bulk_fit = False
            failed_surface_fit = False
            for component in components:
                coefficient_key = min(
                    (
                        candidate
                        for candidate in coefficients
                        if candidate[0] == calibration_interface
                        and candidate[2] == component
                    ),
                    key=lambda candidate: abs(candidate[1] - radius),
                )
                coefficient = coefficients[coefficient_key]
                component_assigned = 0.0
                component_outside = 0.0
                component_smooth = 0.0
                component_modeled = 0.0
                component_outside_delta = 0.0
                component_amplitude = 0.0
                component_bulk_amplitude = 0.0
                for edge, values in edge_values.items():
                    value = values[radius]
                    assigned = value[f"p_total_{component}"]
                    inside = value[f"p_in_{component}"]
                    outside = assigned - inside
                    if outside < -1.0e-12 * max(
                        abs(assigned), abs(inside), 1.0e-300
                    ):
                        raise RuntimeError(
                            f"Local {component} inside participation exceeds assigned "
                            f"participation for interface {interface}, edge {edge}"
                        )
                    amplitude, residual, _ = descriptor_fits[(edge, component)]
                    _, smooth_density, surface_residual, _ = surface_fits[
                        (edge, component)
                    ]
                    bulk_amplitude = sum(
                        value[f"p_bulk_{channel}_ann"]
                        for channel in DESCRIPTOR_CHANNELS[
                            (interface_type, component)
                        ]
                    ) / radius
                    _, fit_confidence = common.edge_model_confidence(
                        amplitude,
                        bulk_amplitude,
                        residual,
                        surface_residual,
                        fit_residual_scale,
                    )
                    component_assigned += assigned
                    component_outside += max(0.0, outside)
                    component_smooth += max(0.0, smooth_density * radius)
                    component_modeled += (
                        coefficient["inside_singular_coefficient_m"] * amplitude
                    )
                    component_outside_delta += (
                        coefficient["outside_delta_coefficient_m"] * amplitude
                    )
                    component_amplitude += amplitude
                    component_bulk_amplitude += bulk_amplitude
                    confidence_weighted_assigned_sum += fit_confidence * assigned
                    if math.isfinite(residual):
                        weighted_bulk_residual += assigned * residual
                    elif assigned > 0.0:
                        failed_bulk_fit = True
                    if math.isfinite(surface_residual):
                        weighted_surface_residual += assigned * surface_residual
                    elif assigned > 0.0:
                        failed_surface_fit = True
                component_corrected = (
                    component_outside
                    + component_smooth
                    + component_modeled
                    + component_outside_delta
                )
                result.update(
                    {
                        f"p_raw_{component}": component_assigned,
                        f"p_out_{component}": component_outside,
                        f"p_smooth_in_{component}": component_smooth,
                        f"p_modeled_{component}": component_modeled,
                        f"p_outside_delta_{component}": component_outside_delta,
                        f"p_corrected_{component}": component_corrected,
                        f"descriptor_amplitude_{component}_per_m": component_amplitude,
                    }
                )
                total_assigned += component_assigned
                corrected += component_corrected
                corrected_outside += component_outside
                smooth_inside_sum += component_smooth
                modeled_singular_sum += component_modeled
                outside_delta_sum += component_outside_delta
                amplitude_total += component_amplitude
                bulk_amplitude_total += component_bulk_amplitude

            if not math.isclose(
                total_assigned, raw, rel_tol=1.0e-8, abs_tol=1.0e-15
            ):
                raise RuntimeError(
                    f"Polarized local assigned participations sum to "
                    f"{total_assigned:.16e}, but interface {interface} participation "
                    f"is {raw:.16e}"
                )
            channel_amplitudes = {
                channel: sum(
                    channel_fits[(edge, channel)][0] for edge in edge_values
                )
                for channel in CHANNELS
            }
            for channel, amplitude in channel_amplitudes.items():
                result[f"descriptor_{channel}_per_m"] = amplitude
            target_channel_total = sum(channel_amplitudes.values())
            calibration_channel_total = sum(
                calibration_channel_amplitudes.values()
            )
            if target_channel_total > 0.0 and calibration_channel_total > 0.0:
                descriptor_mismatch = 0.5 * sum(
                    abs(
                        channel_amplitudes[channel] / target_channel_total
                        - calibration_channel_amplitudes[channel]
                        / calibration_channel_total
                    )
                    for channel in CHANNELS
                )
            else:
                descriptor_mismatch = float(
                    (target_channel_total > 0.0)
                    != (calibration_channel_total > 0.0)
                )
            if interface_type == "SA":
                tangential_amplitude = (
                    channel_amplitudes["top_m"] + channel_amplitudes["top_l"]
                )
                unmodeled_fraction = channel_amplitudes["top_l"] / max(
                    tangential_amplitude, 1.0e-300
                )
            else:
                unmodeled_fraction = 0.0
            topology = common.summarize_edge_topology(edge_values, radius)
            fit_confidence = min(
                1.0,
                confidence_weighted_assigned_sum
                / max(total_assigned, 1.0e-300),
            )
            bulk_fit_residual = (
                math.inf
                if failed_bulk_fit
                else weighted_bulk_residual / total_assigned
                if total_assigned > 0.0
                else math.inf
            )
            surface_fit_residual = (
                math.inf
                if failed_surface_fit
                else weighted_surface_residual / total_assigned
                if total_assigned > 0.0
                else math.inf
            )
            fit_residual = max(bulk_fit_residual, surface_fit_residual)
            singular_fraction = min(
                1.0,
                amplitude_total / max(bulk_amplitude_total, 1.0e-300),
            )
            singular_identifiability = singular_fraction * fit_confidence
            term_scale = (
                abs(corrected_outside)
                + abs(smooth_inside_sum)
                + abs(modeled_singular_sum)
                + abs(outside_delta_sum)
            )
            smooth_fraction = abs(smooth_inside_sum) / max(
                term_scale, 1.0e-300
            )
            modeled_fraction = (
                abs(modeled_singular_sum) + abs(outside_delta_sum)
            ) / max(term_scale, 1.0e-300)
            model_risk = (
                smooth_fraction * (1.0 - fit_confidence)
                + modeled_fraction * (1.0 - singular_identifiability)
            )
            model_confidence = (
                1.0 - min(1.0, model_risk)
            ) * (1.0 - unmodeled_fraction) * (1.0 - descriptor_mismatch) * (
                1.0 - topology["unmodeled_vertex_fraction"]
            )
            if raw > 0.0 and math.isfinite(raw_quality) and raw_quality > 0.0:
                loss_tangent = 1.0 / (raw * raw_quality)
                corrected_quality = (
                    math.inf
                    if corrected <= 0.0
                    else 1.0 / (corrected * loss_tangent)
                )
            else:
                corrected_quality = math.inf
            reference_value = reference.get(key, {}).get("p", math.nan)
            result.update(
                {
                    "p_corrected": corrected,
                    "Q_corrected": corrected_quality,
                    "fit_relative_rms": fit_residual,
                    "bulk_fit_relative_rms": bulk_fit_residual,
                    "surface_fit_relative_rms": surface_fit_residual,
                    "singular_fraction": singular_fraction,
                    "singular_identifiability": singular_identifiability,
                    "fit_confidence": fit_confidence,
                    "descriptor_mismatch": descriptor_mismatch,
                    "unmodeled_longitudinal_fraction": unmodeled_fraction,
                    "model_confidence": model_confidence,
                    "modeled_fraction": modeled_fraction,
                    "edge_count": len(edge_values),
                    **topology,
                    "selected": False,
                    "window_min_um": math.nan,
                    "window_max_um": math.nan,
                    "radius_spread": math.nan,
                    "p_reference": reference_value,
                    "relative_error": (
                        (corrected - reference_value) / reference_value
                        if reference_value > 0.0
                        else math.nan
                    ),
                }
            )
            radius_rows.append(result)
    return radius_rows


def summarize_window(candidates):
    candidates = sorted(candidates, key=lambda row: row["radius_um"])
    result = dict(candidates[0])
    result["method"] = "window-mean" if len(candidates) > 1 else "radius"
    result["radius_um"] = math.nan
    result["selected"] = True
    result["window_min_um"] = candidates[0]["radius_um"]
    result["window_max_um"] = candidates[-1]["radius_um"]
    for key in candidates[0]:
        if (
            key.startswith("p_")
            and key not in {"p_raw", "p_reference"}
            or key.startswith("descriptor_")
        ):
            values = [row[key] for row in candidates]
            if all(isinstance(value, (int, float)) for value in values):
                result[key] = statistics.fmean(values)
    corrected_values = [row["p_corrected"] for row in candidates]
    result["radius_spread"] = (
        max(corrected_values) - min(corrected_values)
    ) / max(abs(result["p_corrected"]), 1.0e-300)
    result["fit_relative_rms"] = statistics.median(
        row["fit_relative_rms"] for row in candidates
    )
    result["bulk_fit_relative_rms"] = statistics.median(
        row["bulk_fit_relative_rms"] for row in candidates
    )
    result["surface_fit_relative_rms"] = statistics.median(
        row["surface_fit_relative_rms"] for row in candidates
    )
    result["singular_fraction"] = statistics.median(
        row["singular_fraction"] for row in candidates
    )
    result["singular_identifiability"] = statistics.median(
        row["singular_identifiability"] for row in candidates
    )
    result["fit_confidence"] = statistics.median(
        row["fit_confidence"] for row in candidates
    )
    result["descriptor_mismatch"] = statistics.median(
        row["descriptor_mismatch"] for row in candidates
    )
    result["modeled_fraction"] = statistics.median(
        row["modeled_fraction"] for row in candidates
    )
    result["unmodeled_longitudinal_fraction"] = statistics.median(
        row["unmodeled_longitudinal_fraction"] for row in candidates
    )
    result["unmodeled_vertex_fraction"] = max(
        row.get("unmodeled_vertex_fraction", 0.0) for row in candidates
    )
    for name in (
        "unmodeled_vertex_length_fraction",
        "unmodeled_vertex_surface_fraction",
        "unmodeled_vertex_bulk_fraction",
    ):
        result[name] = max(row.get(name, 0.0) for row in candidates)
    result["model_confidence"] = statistics.median(
        row["model_confidence"] for row in candidates
    )
    result["relative_error"] = (
        (result["p_corrected"] - result["p_reference"]) / result["p_reference"]
        if result["p_reference"] > 0.0
        else math.nan
    )
    return result


def apply(
    calibration_path,
    input_directory,
    output,
    interface_maps,
    requested_radii_um,
    fit_count,
    fit_residual_scale,
    reference_directory,
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
            f"requested {fit_count}"
        )
    target_order = common.read_fem_order(input_directory)
    if calibration_orders and target_order is not None:
        calibration_order = next(iter(calibration_orders))
        if calibration_order != target_order and not allow_order_mismatch:
            raise RuntimeError(
                f"Calibration uses FEM order p={calibration_order}, but target uses "
                f"p={target_order}; pass --allow-order-mismatch to override"
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
        calibration_metadata[row["interface"]] = metadata
    target_metadata = common.read_interface_metadata(input_directory)
    if not target_metadata:
        raise RuntimeError(
            "Target output directory must contain a resolved config with interface "
            "dielectric metadata"
        )
    configured_map = dict(interface_maps)
    for target_interface in target_metadata:
        calibration_interface = configured_map.get(
            target_interface, target_interface
        )
        if calibration_interface in calibration_metadata:
            common.validate_interface_metadata(
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
    rows = list(radius_rows)
    for candidates in common.group_radius_rows(radius_rows).values():
        for row in candidates:
            row["selected"] = True
        rows.append(summarize_window(candidates))
    common.write_rows(output, rows)
    return rows


def print_calibration(rows):
    print(
        "interface  type  component  descriptor  R (um)  C_inside (m)  "
        "C_outside (m)  bulk RMS  surface RMS"
    )
    for row in rows:
        print(
            f"{row['interface']:9d}  {row['interface_type']:>4}  "
            f"{row['surface_component']:>9}  "
            f"{row['descriptor_channels']:>10}  {row['radius_um']:6.2f}  "
            f"{row['inside_singular_coefficient_m']:12.4e}  "
            f"{row['outside_delta_coefficient_m']:13.4e}  "
            f"{row['bulk_fit_relative_rms']:8.3f}  "
            f"{row['surface_fit_relative_rms']:11.3f}"
        )


def print_application(rows):
    summaries = [row for row in rows if row["method"] != "radius"]
    print(
        "state  interface  type  p_raw       p_corrected  confidence  "
        "mismatch  longitudinal  R window (um)  error"
    )
    for row in summaries:
        error = (
            f"{row['relative_error']:+7.2%}"
            if math.isfinite(row["relative_error"])
            else "    n/a"
        )
        print(
            f"{row['state']:>5}  {row['interface']:9d}  "
            f"{row['interface_type']:>4}  {row['p_raw']:11.4e}  "
            f"{row['p_corrected']:11.4e}  {row['model_confidence']:10.3f}  "
            f"{row['descriptor_mismatch']:8.2%}  "
            f"{row['unmodeled_longitudinal_fraction']:12.2%}  "
            f"{row['window_min_um']:4.2f}-{row['window_max_um']:4.2f}  {error}"
        )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Calibrate and apply the experimental two-sided, polarization-aware "
            "thin-metal edge correction."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    calibrate_parser = subparsers.add_parser("calibrate")
    calibrate_parser.add_argument("--thin", type=Path, required=True)
    calibrate_parser.add_argument("--fabricated", type=Path, required=True)
    calibrate_parser.add_argument("--output", type=Path, required=True)
    calibrate_parser.add_argument("--fit-count", type=int, default=3)

    apply_parser = subparsers.add_parser("apply")
    apply_parser.add_argument("--calibration", type=Path, required=True)
    apply_parser.add_argument("--input", type=Path, required=True)
    apply_parser.add_argument("--output", type=Path, required=True)
    apply_parser.add_argument("--reference", type=Path)
    apply_parser.add_argument(
        "--interface-map",
        action="append",
        type=common.parse_interface_map,
        default=[],
        metavar="TARGET=CALIBRATION",
    )
    apply_parser.add_argument(
        "--radius",
        action="append",
        type=float,
        default=[],
        metavar="R_UM",
    )
    apply_parser.add_argument("--fit-count", type=int)
    apply_parser.add_argument("--fit-residual-scale", type=float, default=0.1)
    apply_parser.add_argument("--allow-order-mismatch", action="store_true")

    args = parser.parse_args()
    if args.command == "calibrate":
        if args.fit_count < 3:
            parser.error("--fit-count must be at least 3")
        rows = calibrate(
            args.thin, args.fabricated, args.output, args.fit_count
        )
        print_calibration(rows)
    else:
        if args.fit_count is not None and args.fit_count < 3:
            parser.error("--fit-count must be at least 3")
        if args.fit_residual_scale < 0.0:
            parser.error("--fit-residual-scale must be nonnegative")
        rows = apply(
            args.calibration,
            args.input,
            args.output,
            args.interface_map,
            args.radius,
            args.fit_count,
            args.fit_residual_scale,
            args.reference,
            args.allow_order_mismatch,
        )
        print_application(rows)
    print(f"\nWrote {args.output}")


if __name__ == "__main__":
    main()
