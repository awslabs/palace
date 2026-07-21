#!/usr/bin/env python3

"""Compare compact driven-CPW raw, corrected, and fabricated-reference outputs."""

import argparse
import csv
from pathlib import Path
import sys


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cpw2d"))
from check_surface_response_output import check  # noqa: E402


INTERFACE_NAMES = {1: "SA", 2: "MS", 3: "MA"}


def read_rows(path):
    with path.open(newline="") as stream:
        return [
            {key.strip(): float(value) for key, value in row.items()}
            for row in csv.DictReader(stream, skipinitialspace=True)
        ]


def relative_error(value, reference):
    return 100.0 * (value / reference - 1.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--corrected", type=Path, required=True)
    parser.add_argument("--fabricated", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1.0e-10)
    args = parser.parse_args()

    check(args.corrected, args.baseline, args.tolerance)

    corrected_rows = read_rows(args.corrected / "surface-Q-corrected.csv")
    corrected = {int(row["interface"]): row for row in corrected_rows}
    fabricated_surface = read_rows(args.fabricated / "surface-Q.csv")[0]
    fabricated_domain = read_rows(args.fabricated / "domain-E.csv")[0]["E_elec (J)"]

    norm_raw = corrected_rows[0]["E_norm raw (J)"]
    norm_corrected = corrected_rows[0]["E_norm corrected (J)"]
    print("Normalization energy (J)")
    print(f"  thin raw:        {norm_raw:.12e}")
    print(
        f"  thin corrected:  {norm_corrected:.12e}"
        f" ({relative_error(norm_corrected, fabricated_domain):+.3f}%)"
    )
    print(f"  fabricated:      {fabricated_domain:.12e}")
    print()
    print(
        "Interface  p_raw          p_corrected    p_fabricated"
        "   raw error   corrected error"
    )
    for interface, name in INTERFACE_NAMES.items():
        row = corrected[interface]
        p_raw = row["p_surf raw"]
        p_corrected = row["p_surf corrected"]
        p_fabricated = fabricated_surface[f"p_surf[{interface}]"]
        print(
            f"{name:>9}  {p_raw:.6e}  {p_corrected:.6e}  {p_fabricated:.6e}"
            f"  {relative_error(p_raw, p_fabricated):+9.3f}%"
            f"  {relative_error(p_corrected, p_fabricated):+15.3f}%"
        )
        fabricated_energy = p_fabricated * fabricated_domain
        print(
            f"{'energy':>9}  {row['E_surf raw (J)']:.6e}"
            f"  {row['E_surf corrected (J)']:.6e}  {fabricated_energy:.6e}"
            f"  {relative_error(row['E_surf raw (J)'], fabricated_energy):+9.3f}%"
            f"  {relative_error(row['E_surf corrected (J)'], fabricated_energy):+15.3f}%"
        )


if __name__ == "__main__":
    main()
