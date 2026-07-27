#!/usr/bin/env python3

"""Summarize order/mesh convergence for corrected corner validation cases."""

import argparse
import csv
from pathlib import Path


INTERFACES = {1: "SA", 2: "MS", 3: "MA"}
CORRECTED_COLUMNS = {
    "raw": "p_surf raw[{index}]",
    "fixed_trace": "p_surf postprocessed fixed-trace[{index}]",
    "fixed_flux": "p_surf postprocessed fixed-flux[{index}]",
    "self_consistent": "p_surf corrected[{index}]",
}


def read_row(path):
    with path.open(newline="") as stream:
        rows = csv.DictReader(stream, skipinitialspace=True)
        try:
            row = next(rows)
        except StopIteration as exc:
            raise ValueError(f"{path} has no data rows") from exc
    return {key.strip(): float(value) for key, value in row.items()}


def output_path(path, filename):
    if not path:
        return None
    result = Path(path).expanduser().resolve()
    return result / filename if result.is_dir() else result


def parse_case(specification):
    if "=" not in specification:
        raise ValueError(f'Case "{specification}" must use NAME=PATH')
    name, paths = specification.split("=", 1)
    if not name:
        raise ValueError("Case name cannot be empty")
    if "," in paths:
        corrected, fabricated = paths.split(",", 1)
        corrected_path = output_path(corrected, "surface-Q-corrected.csv")
        fabricated_path = output_path(fabricated, "surface-Q.csv")
    else:
        root = Path(paths).expanduser().resolve()
        corrected_path = (
            root / "postpro" / "thin-corrected" / "surface-Q-corrected.csv"
        )
        fabricated_path = (
            root / "postpro" / "fabricated-reference" / "surface-Q.csv"
        )
    corrected = (
        read_row(corrected_path)
        if corrected_path is not None and corrected_path.is_file()
        else None
    )
    fabricated = (
        read_row(fabricated_path)
        if fabricated_path is not None and fabricated_path.is_file()
        else None
    )
    if corrected is None and fabricated is None:
        raise FileNotFoundError(
            f'Case "{name}" has no corrected or fabricated convergence output'
        )
    return name, corrected, fabricated


def relative_change(value, previous):
    if previous is None or previous == 0.0:
        return float("nan")
    return 100.0 * (value / previous - 1.0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        action="append",
        required=True,
        help=(
            "Ordered level as NAME=validation-root or "
            "NAME=corrected-output,fabricated-output"
        ),
    )
    args = parser.parse_args()
    cases = [parse_case(specification) for specification in args.case]

    print(
        "case,interface,quantity,value,reference,error_pct,"
        "change_from_previous_pct"
    )
    previous = {}
    for case_name, corrected, fabricated in cases:
        for index, interface_name in INTERFACES.items():
            reference = (
                fabricated[f"p_surf[{index}]"]
                if fabricated is not None
                else float("nan")
            )
            if fabricated is not None:
                quantity = "fabricated"
                change = relative_change(
                    reference, previous.get((index, quantity))
                )
                print(
                    f"{case_name},{interface_name},{quantity},"
                    f"{reference:.12e},{reference:.12e},+0.000000,"
                    f"{change:+.6f}"
                )
                previous[(index, quantity)] = reference
            if corrected is None:
                continue
            for quantity, column in CORRECTED_COLUMNS.items():
                value = corrected[column.format(index=index)]
                error = (
                    100.0 * (value / reference - 1.0)
                    if fabricated is not None
                    else float("nan")
                )
                change = relative_change(value, previous.get((index, quantity)))
                print(
                    f"{case_name},{interface_name},{quantity},{value:.12e},"
                    f"{reference:.12e},{error:+.6f},{change:+.6f}"
                )
                previous[(index, quantity)] = value


if __name__ == "__main__":
    main()
