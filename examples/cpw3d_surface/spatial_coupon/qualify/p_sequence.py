#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""p-sequence controls of a run vs its reference (gallery-physics-06b p_sequence.py, the
control sources, orders and directories as arguments).

For every control source and observable (domain energy E_ii, whole-interface energies
Q_total MA / MS / SA, SA normal / tangential, participations p_X = Q_total_X / E_ii) the
table lists the values at the low / main / high orders with the increments
d_low = (main - low) / |high|, d_high = (high - main) / |high|, the contraction ratio
r = d_high / d_low and, where 0 < r < 1, the Aitken limit p_inf = high + (high - main)
r / (1 - r) (geometric continuation; optimistic for edge-singular functionals), next to
the reference value and the signed distances (value - ref) / |ref|.

usage: p_sequence.py --main DIR --main-order P --config RUN_CONFIG [--low DIR --low-order P]
       [--high DIR --high-order P] [--reference DIR] --control I ... --out-md PATH --out-json PATH
"""
import argparse
import csv
import json
import math
from pathlib import Path

INTERFACE_TYPES = ("MA", "MS", "SA")
OBSERVABLES = ["E", "Q_MA", "Q_MS", "Q_SA", "Q_SA_normal", "Q_SA_tangential", "p_MA", "p_MS", "p_SA"]
# The sharp-edge extrapolated MA (ma_tail.py) of a radial-shell run: Q_MA_sharp = Q_MA +
# the tail inside the innermost ring, p_MA_sharp = Q_MA_sharp / E; the tail is measured
# at every order from that order's own shells.
SHARP_OBSERVABLES = ["Q_MA_tail", "Q_MA_sharp", "p_MA_sharp"]
GATED_OBSERVABLES = ("E", "p_MA", "p_MS", "p_SA")
# The p-sequence gate of a radial-shell run evaluates p_MA_sharp in place of p_MA (decision 61a).
SHARP_GATED_OBSERVABLES = ("E", "p_MA_sharp", "p_MS", "p_SA")


def rows(path):
    with open(path, newline="") as stream:
        return [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(stream)]


def observables(directory, interface_types, ma_tails=None):
    """Source -> observables of a reducer directory; `interface_types` = the run
    config's interface index -> type map (a participation of a type sums every
    interface of that type); `ma_tails` (ma_tail.tails of the directory) adds the
    sharp-edge MA observables."""
    if interface_types is None:
        raise ValueError("the interface index -> type map of the run config is required (never assumed)")
    names = {int(index): str(name) for index, name in interface_types.items()}
    directory = Path(directory)
    out = {}
    for row in rows(directory / "domain-response-matrix.csv"):
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        if i == j:
            out.setdefault(i, {})["E"] = float(row["Q_ij (J)"])
    # The whole-interface diagonal row per (interface, source): the smallest edge /
    # largest radius group (compare_matrices.diagonal_surface_keys), then summed per type.
    diagonal = {}
    for row in rows(directory / "surface-response-matrix.csv"):
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        if i != j:
            continue
        interface = int(float(row["interface"]))
        if interface not in names:
            raise ValueError(f"surface matrix interface {interface} is not in the config's interface map {names}")
        group = (int(float(row["edge"])), -float(row["R (m)"]))
        current = diagonal.get((interface, i))
        if current is None or group < current[0]:
            diagonal[(interface, i)] = (group, row)
    for (interface, i), (_, row) in sorted(diagonal.items()):
        name = names[interface]
        record = out.setdefault(i, {})
        record[f"Q_{name}"] = record.get(f"Q_{name}", 0.0) + float(row["Q_total_ij (J)"])
        if name == "SA":
            record["Q_SA_normal"] = record.get("Q_SA_normal", 0.0) + float(row["Q_total_ij normal (J)"])
            record["Q_SA_tangential"] = record.get("Q_SA_tangential", 0.0) + float(row["Q_total_ij tangential (J)"])
    for record in out.values():
        for name in INTERFACE_TYPES:
            if f"Q_{name}" in record and record.get("E"):
                record[f"p_{name}"] = record[f"Q_{name}"] / record["E"]
    if ma_tails is not None:
        for i, tail in ma_tails.items():
            record = out.get(i)
            if record is not None and "Q_MA" in record and record.get("E"):
                record["Q_MA_tail"] = tail["Q_MA_tail"]
                record["Q_MA_sharp"] = record["Q_MA"] + tail["Q_MA_tail"]
                record["p_MA_sharp"] = record["Q_MA_sharp"] / record["E"]
    return out


def fmt(x, digits=4):
    return "n/a" if x is None or not (isinstance(x, float) and math.isfinite(x)) else f"{x:.{digits}e}"


def pc(x):
    return "n/a" if x is None or not (isinstance(x, float) and math.isfinite(x)) else f"{100 * x:+.2f}%"


def rel(a, b):
    if a is None or b is None:
        return None
    return (b - a) / abs(a) if a else math.nan


def sequence(low, main, high):
    """d_low, d_high relative to |high|, contraction ratio and Aitken limit (None where undefined)."""
    if high is None or main is None:
        return {"d_low": None, "d_high": None, "r": None, "p_inf": None}
    d_high = (high - main) / abs(high) if high else math.nan
    if low is None:
        return {"d_low": None, "d_high": d_high, "r": None, "p_inf": None}
    d_low = (main - low) / abs(high) if high else math.nan
    r = d_high / d_low if d_low else None
    p_inf = high + (high - main) * r / (1 - r) if r is not None and 0 < r < 1 else None
    return {"d_low": d_low, "d_high": d_high, "r": r, "p_inf": p_inf}


def p_sequence(runs, controls, interface_types, reference_interface_types=None, ma_tails=None, reference_ma_side=None):
    """`runs` = {"low": dir or None, "main": dir, "high": dir or None, "ref": dir or None};
    returns {control: {observable: {values, seq, vs_ref}}}.  The reference matrices are
    labeled by `reference_interface_types` (default: the run's map).  `ma_tails` =
    {run key: ma_tail.tails of that directory} adds the sharp-edge MA observables
    (SHARP_OBSERVABLES); the reference's p_MA_sharp is `reference_ma_side` (ma_tail.
    reference_side: extrapolated, modelled or raw as recorded there)."""
    ma_tails = ma_tails or {}
    data = {key: observables(path, reference_interface_types if key == "ref" and reference_interface_types is not None
                             else interface_types, ma_tails.get(key))
            for key, path in runs.items() if path is not None}
    if "ref" in data and reference_ma_side is not None and ma_tails:
        for i, side in reference_ma_side["PerSource"].items():
            if i in data["ref"] and "p_MA" in data["ref"][i]:
                data["ref"][i]["p_MA_sharp"] = side["p_MA_sharp"]
                data["ref"][i]["Q_MA_sharp"] = side["p_MA_sharp"] * data["ref"][i]["E"]
    names = OBSERVABLES + (SHARP_OBSERVABLES if ma_tails else [])
    summary = {}
    for i in controls:
        summary[i] = {}
        for name in names:
            values = {key: data[key].get(i, {}).get(name) for key in ("low", "main", "high", "ref") if key in data}
            seq = sequence(values.get("low"), values.get("main"), values.get("high"))
            vs_ref = {}
            if "ref" in data:
                vs_ref = {key: rel(values["ref"], values.get(key)) for key in ("low", "main", "high")}
                vs_ref["p_inf"] = rel(values["ref"], seq["p_inf"])
            summary[i][name] = {"values": values, "seq": seq, "vs_ref": vs_ref}
    return summary


def markdown_report(summary, orders, title):
    labels = {key: f"p{orders[key]}" for key in orders}
    lines = [f"# {title}", "",
             f"d_low = ({labels.get('main')} - {labels.get('low')})/|{labels.get('high')}|, d_high = ({labels.get('high')} - "
             f"{labels.get('main')})/|{labels.get('high')}|, r = d_high/d_low, p_inf = Aitken limit (only for 0 < r < 1); "
             '"vs ref" = (value - ref)/|ref|.', ""]
    for i, by_observable in summary.items():
        lines += [f"## source {i}", "",
                  f"| observable | {labels.get('low', 'low')} | {labels.get('main', 'main')} | {labels.get('high', 'high')} | d_low | d_high | r | p_inf | ref | low vs ref | main vs ref | high vs ref | p_inf vs ref |",
                  "|---|" + "---:|" * 12]
        for name, record in by_observable.items():
            values, seq, vs_ref = record["values"], record["seq"], record["vs_ref"]
            lines.append("| " + " | ".join([name, fmt(values.get("low")), fmt(values.get("main")), fmt(values.get("high")),
                                            pc(seq["d_low"]), pc(seq["d_high"]), fmt(seq["r"], 2), fmt(seq["p_inf"]),
                                            fmt(values.get("ref")), pc(vs_ref.get("low")), pc(vs_ref.get("main")),
                                            pc(vs_ref.get("high")), pc(vs_ref.get("p_inf"))]) + " |")
        lines.append("")
    sharp = any("p_MA_sharp" in by_observable for by_observable in summary.values())
    lines += ["## Compact: p_MA, p_MS, E, p_SA at the controls (signed, relative to the reference)", "",
              "| source | quantity | low vs ref | main vs ref | high vs ref | p_inf vs ref | d_low | d_high | r |",
              "|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name in ("p_MA",) + (("p_MA_sharp",) if sharp else ()) + ("p_MS", "E", "p_SA"):
        for i, by_observable in summary.items():
            record = by_observable[name]
            lines.append("| " + " | ".join([str(i), name, pc(record["vs_ref"].get("low")), pc(record["vs_ref"].get("main")),
                                            pc(record["vs_ref"].get("high")), pc(record["vs_ref"].get("p_inf")),
                                            pc(record["seq"]["d_low"]), pc(record["seq"]["d_high"]), fmt(record["seq"]["r"], 2)]) + " |")
    lines.append("")
    return "\n".join(lines) + "\n"


def write_p_sequence(runs, orders, controls, out_md, out_json, *, title, interface_types, reference_interface_types=None,
                     ma_tails=None, reference_ma_side=None):
    summary = p_sequence(runs, controls, interface_types, reference_interface_types, ma_tails=ma_tails,
                         reference_ma_side=reference_ma_side)
    Path(out_md).write_text(markdown_report(summary, orders, title))
    Path(out_json).write_text(json.dumps({"Orders": orders, "Controls": list(controls),
                                          "Sources": {str(i): record for i, record in summary.items()}}, indent=2) + "\n")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--main", required=True)
    parser.add_argument("--main-order", type=int, required=True)
    parser.add_argument("--low")
    parser.add_argument("--low-order", type=int)
    parser.add_argument("--high")
    parser.add_argument("--high-order", type=int)
    parser.add_argument("--reference")
    parser.add_argument("--control", type=int, action="append", required=True)
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--title", default="p-sequence controls")
    parser.add_argument("--config", required=True, help="the run's Palace config (interface index -> type map)")
    args = parser.parse_args(argv)
    interface_types = {int(entry["Index"]): entry["Type"]
                       for entry in json.loads(Path(args.config).read_text())["Boundaries"]["Postprocessing"]["Dielectric"]}
    runs = {"low": args.low, "main": args.main, "high": args.high, "ref": args.reference}
    orders = {key: order for key, order in (("low", args.low_order), ("main", args.main_order), ("high", args.high_order))
              if order is not None}
    write_p_sequence(runs, orders, args.control, args.out_md, args.out_json, title=args.title, interface_types=interface_types)
    print(Path(args.out_md).read_text().splitlines()[-40:])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
