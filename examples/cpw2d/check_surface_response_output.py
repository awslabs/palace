#!/usr/bin/env python3

"""Check that response correction preserves raw output and adds valid corrected data."""

import argparse
import csv
import math
from pathlib import Path


RAW_FILES = (
    "domain-E.csv",
    "surface-Q.csv",
    "terminal-C.csv",
    "terminal-Cinv.csv",
    "terminal-Cm.csv",
    "terminal-V.csv",
)


def read_rows(path):
    with path.open(newline="") as stream:
        return [
            {key.strip(): float(value) for key, value in row.items()}
            for row in csv.DictReader(stream, skipinitialspace=True)
        ]


def close(left, right, tolerance):
    if left == right:
        return True
    if not math.isfinite(left) or not math.isfinite(right):
        return False
    return abs(left - right) <= tolerance * max(abs(left), abs(right), 1.0e-300)


def corrected_groups(rows):
    groups = []
    current_key = None
    for row in rows:
        state_field = "f (GHz)" if "f (GHz)" in row else "m"
        key = (state_field, row[state_field], int(row["exc"]))
        if key != current_key:
            groups.append([])
            current_key = key
        groups[-1].append(row)
    return groups


def check_nonnegative(value, field, row_index, interface):
    if not math.isfinite(value) or value < 0.0:
        raise RuntimeError(
            f"Invalid {field} in row {row_index}, interface {interface}"
        )


def check(corrected_output, baseline_output, tolerance):
    for name in RAW_FILES:
        corrected_path = corrected_output / name
        baseline_path = baseline_output / name
        if corrected_path.exists() != baseline_path.exists():
            raise RuntimeError(f"Raw output presence differs for {name}")
        if corrected_path.exists() and corrected_path.read_bytes() != baseline_path.read_bytes():
            raise RuntimeError(f"Raw output differs for {name}")

    raw_rows = read_rows(corrected_output / "surface-Q.csv")
    corrected_rows = read_rows(corrected_output / "surface-Q-corrected.csv")
    groups = corrected_groups(corrected_rows)
    if len(raw_rows) != len(groups):
        raise RuntimeError(
            "Raw and corrected surface tables have different state/excitation counts"
        )

    for row_index, (raw, group) in enumerate(zip(raw_rows, groups), start=1):
        raw_indices = sorted(
            int(key.removeprefix("p_surf[").removesuffix("]"))
            for key in raw
            if key.startswith("p_surf[")
        )
        corrected_by_interface = {
            int(corrected["interface"]): corrected for corrected in group
        }
        if len(corrected_by_interface) != len(group):
            raise RuntimeError(
                f"Corrected row group {row_index} repeats an interface"
            )
        if set(raw_indices) != set(corrected_by_interface):
            raise RuntimeError(
                f"Raw and corrected interfaces differ in row group {row_index}"
            )

        state_field = "f (GHz)" if "f (GHz)" in raw else "m"
        corrected_state_field = "f (GHz)" if "f (GHz)" in group[0] else "m"
        if (
            state_field != corrected_state_field
            or not close(raw[state_field], group[0][corrected_state_field], tolerance)
        ):
            raise RuntimeError(
                f"Raw and corrected state values differ in row group {row_index}"
            )

        norm_raw = group[0]["E_norm raw (J)"]
        norm_corrected = group[0]["E_norm corrected (J)"]
        fixed_flux_fields = (
            "E_norm corrected fixed-flux (J)",
            "E_surf corrected fixed-flux (J)",
            "p_surf corrected fixed-flux",
            "Q_surf corrected fixed-flux",
            "trace closure spread",
        )
        fixed_flux_columns = [
            field in group[0] for field in fixed_flux_fields
        ]
        if any(fixed_flux_columns) and not all(fixed_flux_columns):
            raise RuntimeError(
                f"Incomplete fixed-flux output in row group {row_index}"
            )
        has_fixed_flux = all(fixed_flux_columns)
        norm_fixed_flux = (
            group[0]["E_norm corrected fixed-flux (J)"]
            if has_fixed_flux
            else None
        )
        if not math.isfinite(norm_raw) or norm_raw <= 0.0:
            raise RuntimeError(
                f"Invalid raw normalization energy in row group {row_index}"
            )
        if not math.isfinite(norm_corrected) or norm_corrected <= 0.0:
            raise RuntimeError(
                f"Invalid corrected normalization energy in row group {row_index}"
            )
        if has_fixed_flux and (
            not math.isfinite(norm_fixed_flux) or norm_fixed_flux <= 0.0
        ):
            raise RuntimeError(
                f"Invalid fixed-flux normalization energy in row group {row_index}"
            )

        for interface in raw_indices:
            corrected = corrected_by_interface[interface]
            raw_value = raw[f"p_surf[{interface}]"]
            copied_value = corrected["p_surf raw"]
            if not close(raw_value, copied_value, tolerance):
                raise RuntimeError(
                    f"Raw participation differs in row {row_index}, "
                    f"interface {interface}"
                )
            raw_quality = raw[f"Q_surf[{interface}]"]
            if not close(raw_quality, corrected["Q_surf raw"], tolerance):
                raise RuntimeError(
                    f"Raw quality factor differs in row {row_index}, "
                    f"interface {interface}"
                )
            if (
                not close(norm_raw, corrected["E_norm raw (J)"], tolerance)
                or not close(
                    norm_corrected,
                    corrected["E_norm corrected (J)"],
                    tolerance,
                )
                or (
                    has_fixed_flux
                    and not close(
                        norm_fixed_flux,
                        corrected["E_norm corrected fixed-flux (J)"],
                        tolerance,
                    )
                )
            ):
                raise RuntimeError(
                    f"Normalization energy differs within row group {row_index}"
                )

            energy_raw = corrected["E_surf raw (J)"]
            energy_corrected = corrected["E_surf corrected (J)"]
            participation_corrected = corrected["p_surf corrected"]
            check_nonnegative(energy_raw, "E_surf raw", row_index, interface)
            check_nonnegative(
                energy_corrected, "E_surf corrected", row_index, interface
            )
            check_nonnegative(
                copied_value, "p_surf raw", row_index, interface
            )
            check_nonnegative(
                participation_corrected,
                "p_surf corrected",
                row_index,
                interface,
            )
            if not close(copied_value, energy_raw / norm_raw, tolerance):
                raise RuntimeError(
                    f"Raw energy and participation disagree in row {row_index}, "
                    f"interface {interface}"
                )
            if not close(
                participation_corrected,
                energy_corrected / norm_corrected,
                tolerance,
            ):
                raise RuntimeError(
                    f"Corrected energy and participation disagree in row {row_index}, "
                    f"interface {interface}"
                )

            quality_corrected = corrected["Q_surf corrected"]
            if participation_corrected > 0.0:
                if not math.isfinite(quality_corrected) or quality_corrected <= 0.0:
                    raise RuntimeError(
                        f"Invalid Q_surf corrected in row {row_index}, "
                        f"interface {interface}"
                    )
            elif quality_corrected <= 0.0:
                raise RuntimeError(
                    f"Invalid zero-participation quality factor in row {row_index}, "
                    f"interface {interface}"
                )

            if has_fixed_flux:
                energy_fixed_flux = corrected[
                    "E_surf corrected fixed-flux (J)"
                ]
                participation_fixed_flux = corrected[
                    "p_surf corrected fixed-flux"
                ]
                quality_fixed_flux = corrected[
                    "Q_surf corrected fixed-flux"
                ]
                closure_spread = corrected["trace closure spread"]
                check_nonnegative(
                    energy_fixed_flux,
                    "E_surf corrected fixed-flux",
                    row_index,
                    interface,
                )
                check_nonnegative(
                    participation_fixed_flux,
                    "p_surf corrected fixed-flux",
                    row_index,
                    interface,
                )
                if not close(
                    participation_fixed_flux,
                    energy_fixed_flux / norm_fixed_flux,
                    tolerance,
                ):
                    raise RuntimeError(
                        "Fixed-flux energy and participation disagree in row "
                        f"{row_index}, interface {interface}"
                    )
                if participation_fixed_flux > 0.0:
                    if (
                        not math.isfinite(quality_fixed_flux)
                        or quality_fixed_flux <= 0.0
                    ):
                        raise RuntimeError(
                            f"Invalid Q_surf corrected fixed-flux in row "
                            f"{row_index}, interface {interface}"
                        )
                elif quality_fixed_flux <= 0.0:
                    raise RuntimeError(
                        f"Invalid zero-participation fixed-flux quality factor in "
                        f"row {row_index}, interface {interface}"
                    )
                if (
                    not math.isfinite(closure_spread)
                    or closure_spread < 0.0
                    or closure_spread > 1.0
                ):
                    raise RuntimeError(
                        f"Invalid trace closure spread in row {row_index}, "
                        f"interface {interface}"
                    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("corrected_output", type=Path)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-10)
    args = parser.parse_args()
    check(args.corrected_output, args.baseline_output, args.tolerance)
    print("Raw outputs are unchanged and corrected surface output is valid.")


if __name__ == "__main__":
    main()
