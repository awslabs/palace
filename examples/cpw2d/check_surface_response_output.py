#!/usr/bin/env python3

"""Check that response correction preserves raw output and adds valid corrected data.

Explicit marks (decisions 1908 / 1914; ../surface_response_identification/coverage_marks.py): on a MarkNotApplied tree (any
signal: a twin / exposure / mark column of a postpro table, or config_resolved.json's Coverage.Policy) the E / p cells this
check reads are read through their mark columns ('coverage mark E_norm', 'coverage mark E_surf'): an unmarked E / p cell that is
not a finite number, or a missing mark column, is a REFUSAL (exit 5); a row with a mark 1 is NOT-APPLIED (exit 3: the corrected
output is not checked, never 'valid'). Q_surf keeps this check's own rule (+inf only where p is 0). A Warn / pre-policy tree is
checked as before, except that a text token in a read corrected cell is a REFUSAL. Decision 1934: a mark token is read one
way ('1' / '0' exactly, kept as its text); a table naming a mark / twin / exposure column twice is a REFUSAL.
"""

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "surface_response_identification"))
import coverage_marks as cm  # noqa: E402

# The explicit-mark rule's module (decision 1914): this check is of record only with this file (sha256-16 asserted).
COVERAGE_MARKS_SHA256_16 = "3feaa623f05d4598"
if hashlib.sha256(Path(cm.__file__).read_bytes()).hexdigest()[:16] != COVERAGE_MARKS_SHA256_16:
    raise SystemExit(f"{cm.__file__}: sha256-16 != {COVERAGE_MARKS_SHA256_16} (this check reads marks with that file only)")


RAW_FILES = (
    "domain-E.csv",
    "surface-Q.csv",
    "terminal-C.csv",
    "terminal-Cinv.csv",
    "terminal-Cm.csv",
    "terminal-V.csv",
)


# The suffix of the MarkNotApplied twin columns (kCoveragePartialColumnSuffix, surfaceresponsecoverage.hpp).
PARTIAL_COLUMN_SUFFIX = " partial (uncovered part raw)"
NOT_APPLIED_EXIT = 3


class CorrectionNotApplied(RuntimeError):
    """NOT-APPLIED (decisions 1847 E-1 / 1914): a corrected row is explicitly marked by Coverage.Policy MarkNotApplied (its
    'coverage mark E_norm' or 'coverage mark E_surf' is 1): the corrected output is not checked and never reads as valid."""


# The corrected E / p cells this check reads, with the mark column that states each one's mark (decision 1914).
CHECKED_CORRECTED_CELLS = (
    ("E_norm corrected (J)", cm.MARK_E_NORM),
    ("E_norm corrected fixed-flux (J)", cm.MARK_E_NORM),
    ("E_surf corrected (J)", cm.MARK_E_SURF),
    ("p_surf corrected", cm.MARK_E_SURF),
    ("E_surf corrected fixed-flux (J)", cm.MARK_E_SURF),
    ("p_surf corrected fixed-flux", cm.MARK_E_SURF),
)


def check_marks(rows, marked_tree, path):
    """The explicit-mark reading of the corrected rows (decision 1914): every checked E / p cell through coverage_marks.read_cell
    (a Refusal first, on any row); then CorrectionNotApplied at the first row whose mark is 1. A Warn / pre-policy tree: a text
    token in a checked cell is refused; nothing else changes."""
    for row_index, row in enumerate(rows, start=1):
        where = f"{path} row {row_index} interface {row.get('interface')}"
        for column, mark in CHECKED_CORRECTED_CELLS:
            if column in row:
                cm.read_cell(row, column, mark, marked_tree, where)
    if not marked_tree:
        return
    for row_index, row in enumerate(rows, start=1):
        where = f"{path} row {row_index} interface {row.get('interface')}"
        for mark in (cm.MARK_E_NORM, cm.MARK_E_SURF):
            if cm.mark_is_set(row, mark, where):
                raise CorrectionNotApplied(
                    f"NOT-APPLIED: corrected row {row_index}, interface {row.get('interface')} is marked not applied by the "
                    f"coverage policy ('{mark}' = 1; the values are in the '...{PARTIAL_COLUMN_SUFFIX}' columns): the corrected "
                    "output is not checked"
                )


def marked_tree_signals(output):
    """The MarkNotApplied signals of a postpro directory (coverage_marks.tree_signals with its config_resolved.json)."""
    config_path = output / "config_resolved.json"
    config = json.loads(config_path.read_text()) if config_path.is_file() else None
    return cm.tree_signals(str(output), config)


def _cell(value):
    """A number, or the text of a text column (the coverage Type / status columns)."""
    try:
        return float(value)
    except ValueError:
        return value.strip()


def read_rows(path):
    """The rows of a CSV (values parsed by _cell; a MARK column's token kept as its text: the one reading of
    coverage_marks.mark_token). A header naming a mark / twin / exposure column twice is a Refusal (decision 1934)."""
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        cm.refuse_duplicate_columns([name.strip() for name in reader.fieldnames or []], str(path))
        return [
            {
                key.strip(): (value if cm.is_mark_header(key.strip()) else _cell(value))
                for key, value in row.items()
            }
            for row in reader
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
    check_marks(corrected_rows, bool(marked_tree_signals(corrected_output)), corrected_output / "surface-Q-corrected.csv")
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
        # The fixed-flux columns are parked (the v1 cut): written only with the developer switch.
        fixed_flux_fields = (
            "E_norm corrected fixed-flux (J)",
            "E_surf corrected fixed-flux (J)",
            "p_surf corrected fixed-flux",
            "Q_surf corrected fixed-flux",
        )
        fixed_flux_columns = [
            field in group[0] for field in fixed_flux_fields
        ]
        if any(fixed_flux_columns) and not all(fixed_flux_columns):
            raise RuntimeError(
                f"Incomplete fixed-flux output in row group {row_index}"
            )
        has_fixed_flux = all(fixed_flux_columns)
        has_closure_spread = "trace closure spread" in group[0]
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
            if has_closure_spread:
                closure_spread = corrected["trace closure spread"]
                if (
                    not math.isfinite(closure_spread)
                    or closure_spread < 0.0
                    or closure_spread > 1.0
                ):
                    raise RuntimeError(
                        f"Invalid trace closure spread in row {row_index}, "
                        f"interface {interface}"
                    )


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("corrected_output", type=Path)
    parser.add_argument("baseline_output", type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-10)
    args = parser.parse_args(argv)
    try:
        check(args.corrected_output, args.baseline_output, args.tolerance)
    except CorrectionNotApplied as error:
        print(error)
        return NOT_APPLIED_EXIT
    except cm.Refusal as error:
        print(error)
        return cm.REFUSAL_EXIT
    print("Raw outputs are unchanged and corrected surface output is valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
