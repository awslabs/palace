#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


INTERFACE_NAMES = {1: "SA", 2: "MS", 3: "MA"}
PARTICIPATION_RE = re.compile(r"p_surf\[(\d+)\]$")


def read_rows(path):
    with path.open(newline="") as stream:
        rows = list(csv.reader(stream))
    if len(rows) < 2:
        raise RuntimeError(f"No data rows found in {path}")
    header = [cell.strip() for cell in rows[0]]
    return [
        dict(zip(header, (cell.strip() for cell in row)))
        for row in rows[1:]
        if any(cell.strip() for cell in row)
    ]


def read_raw_participations(directory):
    path = directory / "surface-Q.csv"
    rows = read_rows(path)
    row = rows[-1]
    values = {}
    for label, value in row.items():
        match = PARTICIPATION_RE.fullmatch(label)
        if match:
            values[int(match.group(1))] = float(value)
    if not values:
        raise RuntimeError(f"No participation columns found in {path}")
    return values


def read_edge_participations(directory):
    path = directory / "surface-Q-edge.csv"
    rows = read_rows(path)
    last_iteration = float(rows[-1]["i"])
    values = {}
    for row in rows:
        if float(row["i"]) != last_iteration or int(row["exc"]) != 0:
            continue
        key = (int(row["interface"]), float(row["R (m)"]))
        values[key] = (float(row["p_out"]), float(row["p_ann"]))
    if not values:
        raise RuntimeError(f"No final-iteration edge data found in {path}")
    return values


def calculate(thin_directory, fabricated_directory, proxy_interface):
    thin_raw = read_raw_participations(thin_directory)
    fabricated_raw = read_raw_participations(fabricated_directory)
    thin_edge = read_edge_participations(thin_directory)
    fabricated_edge = read_edge_participations(fabricated_directory)
    if thin_edge.keys() != fabricated_edge.keys():
        raise RuntimeError("Thin and fabricated edge tables do not have matching entries")

    results = []
    for interface, radius in sorted(thin_edge):
        thin_out, thin_annulus = thin_edge[(interface, radius)]
        fabricated_out, fabricated_annulus = fabricated_edge[(interface, radius)]
        _, proxy_thin_annulus = thin_edge[(proxy_interface, radius)]
        if fabricated_annulus <= 0.0:
            raise RuntimeError(
                f"Nonpositive fabricated annulus participation for interface "
                f"{interface}, R={radius}"
            )
        fabricated_inside = fabricated_raw[interface] - fabricated_out
        fabricated_shape_factor = fabricated_inside / fabricated_annulus
        annulus_scale = fabricated_annulus / thin_annulus
        outside_scale = fabricated_out / thin_out
        transfer_coefficient = fabricated_inside / proxy_thin_annulus
        corrected = (
            outside_scale * thin_out
            + transfer_coefficient * proxy_thin_annulus
        )
        reference = fabricated_raw[interface]
        results.append(
            {
                "interface": INTERFACE_NAMES.get(interface, str(interface)),
                "interface_index": interface,
                "radius_um": 1.0e6 * radius,
                "thin_raw": thin_raw[interface],
                "fabricated_raw": reference,
                "thin_outside": thin_out,
                "thin_annulus": thin_annulus,
                "fabricated_outside": fabricated_out,
                "fabricated_annulus": fabricated_annulus,
                "fabricated_inside": fabricated_inside,
                "fabricated_shape_factor": fabricated_shape_factor,
                "annulus_scale": annulus_scale,
                "proxy_interface": INTERFACE_NAMES.get(
                    proxy_interface, str(proxy_interface)
                ),
                "proxy_thin_annulus": proxy_thin_annulus,
                "outside_scale": outside_scale,
                "transfer_coefficient": transfer_coefficient,
                "corrected": corrected,
                "raw_relative_error": (thin_raw[interface] - reference) / reference,
                "relative_error": (corrected - reference) / reference,
            }
        )
    return results


def apply_calibration(calibration, thin_directory, fabricated_directory):
    thin_raw = read_raw_participations(thin_directory)
    fabricated_raw = read_raw_participations(fabricated_directory)
    thin_edge = read_edge_participations(thin_directory)
    proxy_interface = calibration[0]["proxy_interface"]
    proxy_index = next(
        index for index, name in INTERFACE_NAMES.items() if name == proxy_interface
    )

    results = []
    for coefficients in calibration:
        interface = coefficients["interface_index"]
        radius_um = coefficients["radius_um"]
        radius = 1.0e-6 * radius_um
        key = min(
            (key for key in thin_edge if key[0] == interface),
            key=lambda key: abs(key[1] - radius),
        )
        proxy_key = min(
            (key for key in thin_edge if key[0] == proxy_index),
            key=lambda key: abs(key[1] - radius),
        )
        thin_out, _ = thin_edge[key]
        _, proxy_annulus = thin_edge[proxy_key]
        corrected = (
            coefficients["outside_scale"] * thin_out
            + coefficients["transfer_coefficient"] * proxy_annulus
        )
        reference = fabricated_raw[interface]
        results.append(
            {
                "interface": INTERFACE_NAMES.get(interface, str(interface)),
                "interface_index": interface,
                "radius_um": radius_um,
                "thin_raw": thin_raw[interface],
                "fabricated_raw": reference,
                "thin_outside": thin_out,
                "proxy_thin_annulus": proxy_annulus,
                "outside_scale": coefficients["outside_scale"],
                "transfer_coefficient": coefficients["transfer_coefficient"],
                "corrected": corrected,
                "raw_relative_error": (thin_raw[interface] - reference) / reference,
                "relative_error": (corrected - reference) / reference,
            }
        )
    return results


def calculate_convergence(thin_directory, calibration):
    coefficients = {
        (row["interface_index"], row["radius_um"]): row["transfer_coefficient"]
        for row in calibration
    }
    outside_scales = {
        (row["interface_index"], row["radius_um"]): row["outside_scale"]
        for row in calibration
    }
    references = {
        row["interface_index"]: row["fabricated_raw"] for row in calibration
    }
    proxy_interface = calibration[0]["proxy_interface"]
    proxy_index = next(
        index for index, name in INTERFACE_NAMES.items() if name == proxy_interface
    )

    iteration_directories = []
    for directory in thin_directory.glob("iteration*"):
        match = re.fullmatch(r"iteration(\d+)", directory.name)
        if match:
            iteration_directories.append((int(match.group(1)), directory))
    iteration_directories.sort()
    iteration_directories.append((len(iteration_directories) + 1, thin_directory))

    rows = []
    for iteration, directory in iteration_directories:
        raw = read_raw_participations(directory)
        edge = read_edge_participations(directory)
        for interface, radius in sorted(edge):
            radius_um = 1.0e6 * radius
            thin_out, _ = edge[(interface, radius)]
            _, proxy_annulus = edge[(proxy_index, radius)]
            corrected = (
                outside_scales[(interface, radius_um)] * thin_out
                + coefficients[(interface, radius_um)] * proxy_annulus
            )
            reference = references[interface]
            rows.append(
                {
                    "iteration": iteration,
                    "interface": INTERFACE_NAMES.get(interface, str(interface)),
                    "interface_index": interface,
                    "radius_um": radius_um,
                    "thin_raw": raw[interface],
                    "corrected": corrected,
                    "fabricated_reference": reference,
                    "raw_relative_error": (raw[interface] - reference) / reference,
                    "corrected_relative_error": (corrected - reference) / reference,
                }
            )
    return rows


def write_results(path, results):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)


def print_convergence_summary(rows, requested_radius):
    radii = sorted({row["radius_um"] for row in rows})
    radius = min(radii, key=lambda value: abs(value - requested_radius))
    print(f"\nAMR spread at R = {radius:g} um:")
    print(f"{'interface':>9} {'raw spread':>12} {'corrected spread':>18}")
    interfaces = sorted({row["interface_index"] for row in rows})
    for interface in interfaces:
        selected = [
            row
            for row in rows
            if row["interface_index"] == interface and row["radius_um"] == radius
        ]
        reference = selected[0]["fabricated_reference"]
        raw = [row["thin_raw"] for row in selected]
        corrected = [row["corrected"] for row in selected]
        print(
            f"{INTERFACE_NAMES.get(interface, str(interface)):>9} "
            f"{(max(raw) - min(raw)) / reference:12.3%} "
            f"{(max(corrected) - min(corrected)) / reference:18.3%}"
        )


def print_results(results, title):
    print(title)
    print(
        f"{'interface':>9} {'R (um)':>8} {'S_out':>8} {'C_inner':>10} "
        f"{'p_corrected':>15} {'raw error':>10} {'corr error':>10}"
    )
    for row in results:
        print(
            f"{row['interface']:>9} {row['radius_um']:8.3f} "
            f"{row['outside_scale']:8.4f} "
            f"{row['transfer_coefficient']:10.5f} {row['corrected']:15.7e} "
            f"{row['raw_relative_error']:10.3%} {row['relative_error']:10.3%}"
        )


def main():
    parser = argparse.ArgumentParser(
        description="Calibrate thin-metal surface participation against a fabricated CPW."
    )
    parser.add_argument(
        "--thin",
        type=Path,
        default=Path("postpro/surface_calibration_thin"),
        help="Thin-metal Palace output directory.",
    )
    parser.add_argument(
        "--fabricated",
        type=Path,
        default=Path("postpro/surface_calibration_fabricated"),
        help="Fabrication-resolved Palace output directory.",
    )
    parser.add_argument(
        "--target-thin",
        type=Path,
        default=Path("postpro/surface_validation_thin"),
        help="Thin-metal output directory for cross-geometry validation.",
    )
    parser.add_argument(
        "--target-fabricated",
        type=Path,
        default=Path("postpro/surface_validation_fabricated"),
        help="Fabrication-resolved reference for cross-geometry validation.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("postpro/surface-calibration.csv"),
        help="Output CSV path.",
    )
    parser.add_argument(
        "--convergence-output",
        type=Path,
        default=Path("postpro/surface-calibration-convergence.csv"),
        help="Output CSV path for the thin-metal AMR history.",
    )
    parser.add_argument(
        "--validation-output",
        type=Path,
        default=Path("postpro/surface-calibration-validation.csv"),
        help="Output CSV path for cross-geometry validation.",
    )
    parser.add_argument(
        "--proxy-interface",
        type=int,
        choices=sorted(INTERFACE_NAMES),
        default=1,
        help="Thin interface index used as the shared edge-amplitude proxy.",
    )
    parser.add_argument(
        "--summary-radius",
        type=float,
        default=0.2,
        help="Matching radius in um used for the printed AMR-spread summary.",
    )
    args = parser.parse_args()

    results = calculate(args.thin, args.fabricated, args.proxy_interface)
    validation = apply_calibration(
        results, args.target_thin, args.target_fabricated
    )
    convergence = calculate_convergence(args.thin, results)
    write_results(args.output, results)
    write_results(args.convergence_output, convergence)
    write_results(args.validation_output, validation)
    print_results(validation, "Transferred 20/12 um CPW validation:")
    print_convergence_summary(convergence, args.summary_radius)
    print(f"\nWrote {args.output}")
    print(f"Wrote {args.convergence_output}")
    print(f"Wrote {args.validation_output}")


if __name__ == "__main__":
    main()
