#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Entry-wise comparison of two reducer outputs (domain-response-matrix.csv and
surface-response-matrix.csv) of the same coupon: the acceptance of the per-coupon source
split (decision 61b) - the matrices of a coupon reduced from the union of N worker
jobs' archives must equal a single job's to roundoff.  Rows are matched by their key
columns (basis_i, basis_j; interface, edge, basis_i, basis_j), every numeric column is
compared, and the record carries the number of entries, the exact matches, the largest
absolute and relative differences per column (and per source = basis_i) and the verdict
against --tolerance: the
largest absolute difference of a column relative to the column's largest entry (default
1e-9, roundoff of a 12-significant-digit print; the per-entry relative differences are
reported - a near-zero cross term may differ by more at the same absolute roundoff).
The CSVs print 12 significant digits, so equal values are bit-for-bit equal as text.

usage: compare_split_matrices.py DIR_A DIR_B [--tolerance REL] [--out PATH]
"""
import argparse
import csv
import json
import math
from pathlib import Path

MATRICES = {"domain": ("domain-response-matrix.csv", ("basis_i", "basis_j")),
            "surface": ("surface-response-matrix.csv", ("interface", "edge", "basis_i", "basis_j"))}
DEFAULT_RELATIVE_TOLERANCE = 1e-9


def read_rows(path):
    with open(path, newline="") as stream:
        reader = csv.reader(stream)
        header = [name.strip() for name in next(reader)]
        rows = [[cell.strip() for cell in row] for row in reader if row]
    return header, rows


def key_of(header, row, keys):
    columns = {name.split(" ")[0]: index for index, name in enumerate(header)}
    return tuple(int(float(row[columns[key]])) for key in keys)


def compare_matrix(path_a, path_b, keys, tolerance):
    header_a, rows_a = read_rows(path_a)
    header_b, rows_b = read_rows(path_b)
    if header_a != header_b:
        raise ValueError(f"headers differ: {header_a} vs {header_b}")
    key_columns = {header_a.index(next(name for name in header_a if name.split(" ")[0] == key)) for key in keys}
    value_columns = [index for index in range(len(header_a)) if index not in key_columns]
    by_key_b = {key_of(header_b, row, keys): row for row in rows_b}
    if len(by_key_b) != len(rows_b):
        raise ValueError(f"{path_b}: repeated row keys")
    if {key_of(header_a, row, keys) for row in rows_a} != set(by_key_b):
        raise ValueError(f"{path_a} and {path_b} hold different row keys")
    columns = {}
    for index in value_columns:
        name = header_a[index]
        exact = 0
        max_abs = max_rel = 0.0
        worst = None
        scale = max((abs(float(row[index])) for row in rows_a), default=0.0)
        for row in rows_a:
            other = by_key_b[key_of(header_a, row, keys)]
            if row[index] == other[index]:
                exact += 1
            a, b = float(row[index]), float(other[index])
            difference = abs(a - b)
            relative = difference / max(abs(a), abs(b)) if max(abs(a), abs(b)) > 0.0 else 0.0
            if difference > max_abs:
                max_abs = difference
            if relative > max_rel:
                max_rel = relative
                worst = {"Key": key_of(header_a, row, keys), "A": a, "B": b}
        per_source = {}
        for row in rows_a:
            other = by_key_b[key_of(header_a, row, keys)]
            a, b = float(row[index]), float(other[index])
            relative = abs(a - b) / max(abs(a), abs(b)) if max(abs(a), abs(b)) > 0.0 else 0.0
            source = key_of(header_a, row, keys)[keys.index("basis_i")]
            per_source[source] = max(per_source.get(source, 0.0), relative)
        columns[name] = {"Entries": len(rows_a), "ExactText": exact, "MaxAbsoluteDifference": max_abs,
                         "MaxRelativeDifference": max_rel, "MaxRelativeToLargestEntry": (max_abs / scale if scale > 0.0 else 0.0),
                         "LargestEntry": scale, "Worst": worst,
                         "PerSourceMaxRelativeDifference": {str(source): per_source[source] for source in sorted(per_source)},
                         "WithinTolerance": bool((max_abs / scale if scale > 0.0 else 0.0) <= tolerance or not math.isfinite(tolerance))}
    return {"A": str(path_a), "B": str(path_b), "Rows": len(rows_a), "Columns": columns,
            "WithinTolerance": all(column["WithinTolerance"] for column in columns.values())}


def compare(directory_a, directory_b, *, tolerance=DEFAULT_RELATIVE_TOLERANCE):
    record = {"Tolerance": tolerance, "Rule": ("rows matched by key columns, every numeric column compared; equal when "
                                                "the largest relative entry difference is within the tolerance"),
              "Matrices": {}}
    for kind, (name, keys) in MATRICES.items():
        record["Matrices"][kind] = compare_matrix(Path(directory_a) / name, Path(directory_b) / name, keys, tolerance)
    record["Equal"] = all(matrix["WithinTolerance"] for matrix in record["Matrices"].values())
    record["MaxRelativeDifference"] = max(column["MaxRelativeDifference"] for matrix in record["Matrices"].values()
                                          for column in matrix["Columns"].values())
    record["MaxRelativeToLargestEntry"] = max(column["MaxRelativeToLargestEntry"] for matrix in record["Matrices"].values()
                                              for column in matrix["Columns"].values())
    record["ExactTextFraction"] = (sum(column["ExactText"] for matrix in record["Matrices"].values() for column in matrix["Columns"].values())
                                   / max(1, sum(column["Entries"] for matrix in record["Matrices"].values() for column in matrix["Columns"].values())))
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory_a", type=Path)
    parser.add_argument("directory_b", type=Path)
    parser.add_argument("--tolerance", type=float, default=DEFAULT_RELATIVE_TOLERANCE)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    record = compare(args.directory_a, args.directory_b, tolerance=args.tolerance)
    if args.out:
        args.out.write_text(json.dumps(record, indent=2) + "\n")
    for kind, matrix in record["Matrices"].items():
        for name, column in matrix["Columns"].items():
            print(f"{kind} {name}: {column['Entries']} entries, {column['ExactText']} bit-for-bit, max |rel| "
                  f"{column['MaxRelativeDifference']:.3e}, max |abs| {column['MaxAbsoluteDifference']:.3e} "
                  f"({column['MaxRelativeToLargestEntry']:.3e} of the largest entry)")
    print(f"EQUAL {record['Equal']} (max difference relative to the largest entry {record['MaxRelativeToLargestEntry']:.3e}, "
          f"tolerance {record['Tolerance']:.1e}; max per-entry relative difference {record['MaxRelativeDifference']:.3e})")
    return 0 if record["Equal"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
