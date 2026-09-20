#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Key-source table of a qualification (four-edge-physics-11 key_sources.py without its
hand-picked source lists): per offset quantity the free within-1 / 2 / 5% counts, the
worst and second-worst free source, the largest |offset| sources and - when other
pairings are given (e.g. the run vs an earlier mesh) - the largest movers.  Sources
named explicitly with --source are listed first.

usage: key_sources.py --out PATH [--zero-trace I ...] [--source I ...] [--top N]
       label=comparison.json [label=comparison.json ...]   (the first pairing is the main one)
"""
import argparse
import json
import math
from pathlib import Path

OFFSET_KEYS = ("E_rel", "p_SA_rel", "p_MS_rel", "p_MA_rel")
DEFAULT_TOP = 8


def pc(x):
    return "n/a" if x is None or not (isinstance(x, float) and math.isfinite(x)) else f"{100 * x:+.2f}%"


def key_sources(comparisons, zero_trace, *, sources=(), top=DEFAULT_TOP):
    """`comparisons` = {label: PerSource} (the first is the main pairing); returns the
    record {quantity: {Counts per label, Worst, Largest, Movers per other label, Named}}."""
    main_label = next(iter(comparisons))
    main = comparisons[main_label]
    zero_trace = set(zero_trace)
    free = sorted(int(i) for i in main if int(i) not in zero_trace)
    record = {"MainPairing": main_label, "Free": free, "Named": list(sources), "Quantities": {}}
    for key in OFFSET_KEYS:
        entry = {"Counts": {}, "Named": {}, "Largest": [], "Movers": {}}
        for label, data in comparisons.items():
            values = {i: data[str(i)][key] for i in free if str(i) in data and key in data[str(i)]
                      and math.isfinite(data[str(i)][key])}
            if not values:
                continue
            ranked = sorted(values, key=lambda i: -abs(values[i]))
            entry["Counts"][label] = {"n": len(values),
                                      "within_1pct": sum(abs(v) < 0.01 for v in values.values()),
                                      "within_2pct": sum(abs(v) < 0.02 for v in values.values()),
                                      "within_5pct": sum(abs(v) < 0.05 for v in values.values()),
                                      "Worst": [{"Source": i, "Offset": values[i]} for i in ranked[:2]]}
            if label == main_label:
                entry["Largest"] = [{"Source": i, "Offset": values[i]} for i in ranked[:top]]
            else:
                entry["Movers"][label] = [{"Source": i, "Offset": values[i]} for i in ranked[:top]]
        for i in sources:
            entry["Named"][str(i)] = {label: data.get(str(i), {}).get(key) for label, data in comparisons.items()}
        record["Quantities"][key] = entry
    return record


def markdown_report(record):
    lines = [f"# Key sources of {record['MainPairing']} (free view: {len(record['Free'])} sources)", ""]
    for key, entry in record["Quantities"].items():
        lines += [f"## {key}", ""]
        if record["Named"]:
            labels = list(next(iter(entry["Named"].values())).keys()) if entry["Named"] else []
            lines += ["| source | " + " | ".join(labels) + " |", "|---:|" + "---:|" * len(labels)]
            for i in record["Named"]:
                lines.append(f"| {i} | " + " | ".join(pc(entry['Named'][str(i)][label]) for label in labels) + " |")
            lines.append("")
        for label, counts in entry["Counts"].items():
            worst = ", ".join(f"{pc(w['Offset'])} at {w['Source']}" for w in counts["Worst"])
            lines.append(f"- {label}: within 1% {counts['within_1pct']}, within 2% {counts['within_2pct']}, within 5% "
                         f"{counts['within_5pct']} of {counts['n']}; worst {worst}")
        lines.append("- largest |offset|: " + ", ".join(f"{item['Source']} {pc(item['Offset'])}" for item in entry["Largest"]))
        for label, movers in entry["Movers"].items():
            lines.append(f"- largest movers {label}: " + ", ".join(f"{item['Source']} {pc(item['Offset'])}" for item in movers))
        lines.append("")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, required=True, help="Markdown output (a .json sibling is written too)")
    parser.add_argument("--zero-trace", type=int, action="append", default=[])
    parser.add_argument("--source", type=int, action="append", default=[])
    parser.add_argument("--top", type=int, default=DEFAULT_TOP)
    parser.add_argument("comparisons", nargs="+", metavar="label=comparison.json")
    args = parser.parse_args(argv)
    comparisons = {item.split("=", 1)[0]: json.loads(Path(item.split("=", 1)[1]).read_text())["PerSource"]
                   for item in args.comparisons}
    record = key_sources(comparisons, args.zero_trace, sources=args.source, top=args.top)
    args.out.write_text(markdown_report(record))
    args.out.with_suffix(".json").write_text(json.dumps(record, indent=2) + "\n")
    print(args.out.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
