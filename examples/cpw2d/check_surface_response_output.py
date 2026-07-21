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
    return abs(left - right) <= tolerance * max(abs(left), abs(right), 1.0e-300)


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
    if len(raw_rows) != len(corrected_rows):
        raise RuntimeError("Raw and corrected surface tables have different row counts")

    for row_index, (raw, corrected) in enumerate(zip(raw_rows, corrected_rows), start=1):
        raw_indices = sorted(
            int(key.removeprefix("p_surf[").removesuffix("]"))
            for key in raw
            if key.startswith("p_surf[")
        )
        for interface in raw_indices:
            raw_value = raw[f"p_surf[{interface}]"]
            copied_value = corrected[f"p_surf raw[{interface}]"]
            if not close(raw_value, copied_value, tolerance):
                raise RuntimeError(
                    f"Raw participation differs in row {row_index}, "
                    f"interface {interface}"
                )
            corrected_fields = (
                f"E_surf corrected[{interface}] (J)",
                f"p_surf corrected[{interface}]",
                f"Q_surf corrected[{interface}]",
            )
            for field in corrected_fields:
                value = corrected[field]
                if not math.isfinite(value) or value <= 0.0:
                    raise RuntimeError(
                        f"Invalid {field} in row {row_index}, interface {interface}"
                    )
        domain = corrected["E_elec corrected (J)"]
        if not math.isfinite(domain) or domain <= 0.0:
            raise RuntimeError(f"Invalid corrected domain energy in row {row_index}")


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
