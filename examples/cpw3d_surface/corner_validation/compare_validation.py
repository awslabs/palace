#!/usr/bin/env python3

"""Compare corrected thin-metal corner validation against the fabricated mesh."""

import argparse
import csv
from pathlib import Path


INTERFACES = {1: "SA", 2: "MS", 3: "MA"}


def read_row(path):
    with path.open(newline="") as stream:
        row = next(csv.DictReader(stream, skipinitialspace=True))
    return {key.strip(): float(value) for key, value in row.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--corrected", type=Path, required=True)
    parser.add_argument("--fabricated", type=Path, required=True)
    args = parser.parse_args()
    corrected = read_row(
        args.corrected.expanduser().resolve() / "surface-Q-corrected.csv"
    )
    fabricated = read_row(
        args.fabricated.expanduser().resolve() / "surface-Q.csv"
    )

    print(
        "interface,raw,fixed_trace,fixed_flux,self_consistent,fabricated,"
        "raw_error_pct,fixed_trace_error_pct,fixed_flux_error_pct,"
        "self_consistent_error_pct,closure_spread_pct"
    )
    for index, name in INTERFACES.items():
        raw = corrected[f"p_surf raw[{index}]"]
        fixed_trace = corrected[f"p_surf postprocessed fixed-trace[{index}]"]
        fixed_flux = corrected[f"p_surf postprocessed fixed-flux[{index}]"]
        self_consistent = corrected[f"p_surf corrected[{index}]"]
        reference = fabricated[f"p_surf[{index}]"]
        closure_spread = 100.0 * corrected[f"trace closure spread[{index}]"]
        errors = [
            100.0 * (value / reference - 1.0)
            for value in (raw, fixed_trace, fixed_flux, self_consistent)
        ]
        print(
            f"{name},{raw:.12e},{fixed_trace:.12e},{fixed_flux:.12e},"
            f"{self_consistent:.12e},{reference:.12e},"
            + ",".join(f"{error:+.3f}" for error in errors)
            + f",{closure_spread:.3f}"
        )


if __name__ == "__main__":
    main()
