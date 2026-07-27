#!/usr/bin/env python3

"""Compare compact Maxwell CPW thin and fabricated-reference outputs."""

import argparse
import csv
import math
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
    parser.add_argument(
        "--mode",
        type=int,
        default=1,
        help="Eigenmode index to compare; ignored for driven output",
    )
    args = parser.parse_args()
    if args.mode < 1:
        parser.error("--mode must be positive")

    check(args.corrected, args.baseline, args.tolerance)

    corrected_rows = read_rows(args.corrected / "surface-Q-corrected.csv")
    if corrected_rows and "m" in corrected_rows[0]:
        corrected_rows = [
            row for row in corrected_rows if int(row["m"]) == args.mode
        ]
        if not corrected_rows:
            raise RuntimeError(f"Corrected output has no mode {args.mode}")
    corrected = {int(row["interface"]): row for row in corrected_rows}
    fabricated_surface_rows = read_rows(args.fabricated / "surface-Q.csv")
    fabricated_domain_rows = read_rows(args.fabricated / "domain-E.csv")
    if fabricated_surface_rows and "m" in fabricated_surface_rows[0]:
        fabricated_surface_rows = [
            row for row in fabricated_surface_rows if int(row["m"]) == args.mode
        ]
        fabricated_domain_rows = [
            row for row in fabricated_domain_rows if int(row["m"]) == args.mode
        ]
        if not fabricated_surface_rows or not fabricated_domain_rows:
            raise RuntimeError(f"Fabricated output has no mode {args.mode}")
        print(f"Eigenmode {args.mode}")
    fabricated_surface = fabricated_surface_rows[0]
    fabricated_domain = fabricated_domain_rows[0]["E_elec (J)"]

    norm_raw = corrected_rows[0]["E_norm raw (J)"]
    norm_fixed_trace = corrected_rows[0]["E_norm corrected (J)"]
    norm_fixed_flux = corrected_rows[0]["E_norm corrected fixed-flux (J)"]
    norm_self_consistent = corrected_rows[0].get("E_norm self-consistent (J)")
    has_self_consistent = (
        norm_self_consistent is not None
        and math.isfinite(norm_self_consistent)
    )
    print("Normalization energy (J)")
    print(f"  thin raw:          {norm_raw:.12e}")
    print(
        f"  fixed trace:       {norm_fixed_trace:.12e}"
        f" ({relative_error(norm_fixed_trace, fabricated_domain):+.3f}%)"
    )
    print(
        f"  fixed flux:        {norm_fixed_flux:.12e}"
        f" ({relative_error(norm_fixed_flux, fabricated_domain):+.3f}%)"
    )
    if has_self_consistent:
        print(
            f"  self-consistent:   {norm_self_consistent:.12e}"
            f" ({relative_error(norm_self_consistent, fabricated_domain):+.3f}%)"
        )
    print(f"  fabricated:        {fabricated_domain:.12e}")
    print()
    print(
        "Interface  p_raw          p_fixed_trace  p_fixed_flux  "
        + (" p_self_consist" if has_self_consistent else "")
        + "  p_fabricated   raw error   trace error   flux error"
        + ("   self error" if has_self_consistent else "")
        + "   closure spread"
    )
    for interface, name in INTERFACE_NAMES.items():
        row = corrected[interface]
        p_raw = row["p_surf raw"]
        p_fixed_trace = row["p_surf corrected"]
        p_fixed_flux = row["p_surf corrected fixed-flux"]
        p_self_consistent = row.get("p_surf self-consistent")
        p_fabricated = fabricated_surface[f"p_surf[{interface}]"]
        fabricated_energy = p_fabricated * fabricated_domain
        self_columns = ""
        self_error = ""
        self_energy = ""
        self_energy_error = ""
        if has_self_consistent:
            self_columns = f"  {p_self_consistent:.6e}"
            self_error = (
                f"  {relative_error(p_self_consistent, p_fabricated):+10.3f}%"
            )
            self_energy_value = row["E_surf self-consistent (J)"]
            self_energy = f"  {self_energy_value:.6e}"
            self_energy_error = (
                f"  {relative_error(self_energy_value, fabricated_energy):+10.3f}%"
            )
        print(
            f"{name:>9}  {p_raw:.6e}  {p_fixed_trace:.6e}  {p_fixed_flux:.6e}"
            f"{self_columns}  {p_fabricated:.6e}"
            f"  {relative_error(p_raw, p_fabricated):+9.3f}%"
            f"  {relative_error(p_fixed_trace, p_fabricated):+11.3f}%"
            f"  {relative_error(p_fixed_flux, p_fabricated):+10.3f}%"
            f"{self_error}"
            f"  {100.0 * row['trace closure spread']:+13.3f}%"
        )
        print(
            f"{'energy':>9}  {row['E_surf raw (J)']:.6e}"
            f"  {row['E_surf corrected (J)']:.6e}"
            f"  {row['E_surf corrected fixed-flux (J)']:.6e}"
            f"{self_energy}"
            f"  {fabricated_energy:.6e}"
            f"  {relative_error(row['E_surf raw (J)'], fabricated_energy):+9.3f}%"
            f"  {relative_error(row['E_surf corrected (J)'], fabricated_energy):+15.3f}%"
            "  "
            f"{relative_error(row['E_surf corrected fixed-flux (J)'], fabricated_energy):+10.3f}%"
            f"{self_energy_error}"
        )

    confidence_rows = read_rows(
        args.corrected / "surface-response-confidence.csv"
    )
    if confidence_rows and "m" in confidence_rows[0]:
        confidence_rows = [
            row for row in confidence_rows if int(row["m"]) == args.mode
        ]
    if not confidence_rows:
        raise RuntimeError("Corrected output has no matching confidence row")
    confidence = confidence_rows[-1]
    print()
    print(
        "Closure diagnostic:"
        f" max spread = {100.0 * confidence['max trace closure spread']:.3f}%,"
        f" confidence pass = {int(confidence['confidence pass'])}"
    )
    if has_self_consistent:
        print(
            "Self-consistent diagnostic:"
            " f ="
            f" {confidence['Re{f self-consistent} (GHz)']:.6f}"
            f"{confidence['Im{f self-consistent} (GHz)']:+.6e}i GHz,"
            f" M-overlap = {confidence['self-consistent M-overlap']:.9f},"
            " confidence pass ="
            f" {int(confidence['self-consistent confidence pass'])}"
        )


if __name__ == "__main__":
    main()
