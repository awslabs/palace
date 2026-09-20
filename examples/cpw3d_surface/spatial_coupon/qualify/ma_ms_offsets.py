#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""MA / MS / SA participation and energy offsets of a run vs its reference (and any
other pairing), all sources, with the magnitude-weighted view of p_MA (weights =
reference p_MA; weighted mean and weighted median of the signed relative offset; the
same on the strongest-reference-MA sources) and the within-1/2/5/10% counts
(gallery-physics-06b ma_ms_offsets.py, the pairings and the reference directory as
arguments).

usage: ma_ms_offsets.py --reference-dir DIR --out-md PATH --out-json PATH [--zero-trace I ...]
       [--strongest N] label=comparison.json [label=comparison.json ...]
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics

OFFSET_KEYS = ["p_MA_rel", "p_MS_rel", "p_SA_rel", "E_rel"]
DEFAULT_STRONGEST = 20


def percentile(values, q):
    values = sorted(values)
    if not values:
        return math.nan
    k = (len(values) - 1) * q
    low, high = math.floor(k), math.ceil(k)
    return values[low] + (values[high] - values[low]) * (k - low)


def distribution(values, keys=None):
    pairs = [(v, k) for v, k in zip(values, keys or [None] * len(values)) if isinstance(v, float) and math.isfinite(v)]
    values = [v for v, _ in pairs]
    if not values:
        return {"n": 0}
    worst = max(pairs, key=lambda pair: abs(pair[0]))
    return {"n": len(values), "signed_median": statistics.median(values), "p10": percentile(values, 0.1),
            "p90": percentile(values, 0.9), "min": min(values), "max": max(values),
            "abs_median": statistics.median(abs(v) for v in values), "abs_max": abs(worst[0]), "worst_source": worst[1],
            "negative": sum(v < 0 for v in values),
            "within_1pct": sum(abs(v) < 0.01 for v in values), "within_2pct": sum(abs(v) < 0.02 for v in values),
            "within_5pct": sum(abs(v) < 0.05 for v in values), "within_10pct": sum(abs(v) < 0.10 for v in values)}


def weighted(offsets, weights):
    """Weighted mean (= relative offset of the weight-aggregated p_MA when the weights are
    the reference p_MA) and weighted median of the signed offsets."""
    pairs = sorted((o, w) for o, w in zip(offsets, weights) if math.isfinite(o))
    total = sum(w for _, w in pairs)
    mean = sum(o * w for o, w in pairs) / total
    accumulated = 0.0
    median = pairs[-1][0]
    for offset, weight in pairs:
        accumulated += weight
        if accumulated >= total / 2:
            median = offset
            break
    return {"n": len(pairs), "weighted_mean": mean, "weighted_median": median, "weight_total": total}


def reference_p_ma(reference_dir):
    """Source -> reference p_MA = Q_total_ii(MA) / E_ii."""
    reference_dir = Path(reference_dir)
    energies, ma = {}, {}
    with (reference_dir / "domain-response-matrix.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            row = {k.strip(): v.strip() for k, v in row.items()}
            if int(float(row["basis_i"])) == int(float(row["basis_j"])):
                energies[int(float(row["basis_i"]))] = float(row["Q_ij (J)"])
    with (reference_dir / "surface-response-matrix.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            row = {k.strip(): v.strip() for k, v in row.items()}
            if int(float(row["interface"])) == 1 and int(float(row["basis_i"])) == int(float(row["basis_j"])):
                ma[int(float(row["basis_i"]))] = float(row["Q_total_ij (J)"])
    return {i: ma[i] / energies[i] for i in ma if energies.get(i)}


def strongest_sources(ref_pma, free, count):
    """The `count` strongest reference-MA free sources (every free source when fewer)."""
    return sorted(free, key=lambda i: -ref_pma[i])[:count]


def pc(x):
    return "n/a" if x is None or not (isinstance(x, float) and math.isfinite(x)) else f"{100 * x:+.2f}%"


def offsets(comparisons, main_pairing, zero_trace, ref_pma, *, strongest=DEFAULT_STRONGEST):
    """The distributions / weighted view record; `comparisons` = {pairing: PerSource}."""
    per_source = comparisons[main_pairing]
    sources = sorted(int(i) for i in per_source)
    zero_trace = set(zero_trace)
    free = [i for i in sources if i not in zero_trace and i in ref_pma]
    strong = strongest_sources(ref_pma, free, strongest)
    out = {"MainPairing": main_pairing, "Distributions": {}, "Weighted": {}, "StrongestSources": strong,
           "StrongestCount": strongest, "Free": free, "ZeroTrace": sorted(zero_trace & set(sources)),
           "ReferencePMA": {str(i): ref_pma[i] for i in free}}
    for pairing, data in comparisons.items():
        out["Distributions"][pairing] = {}
        for key in OFFSET_KEYS:
            for view, selection in (("free", free), ("raw", sources), (f"strongest{strongest}", strong)):
                selection = [i for i in selection if str(i) in data]
                out["Distributions"][pairing][f"{key}:{view}"] = distribution([data[str(i)].get(key) for i in selection], selection)
        if "p_MA_rel" in next(iter(data.values()), {}):
            out["Weighted"][pairing] = {}
            for view, selection in (("free", free), (f"strongest{strongest}", strong)):
                selection = [i for i in selection if str(i) in data and "p_MA_rel" in data[str(i)]]
                if not selection:
                    continue
                values = [data[str(i)]["p_MA_rel"] for i in selection]
                record = weighted(values, [ref_pma[i] for i in selection])
                dist = distribution(values, selection)
                record.update({"unweighted_median": dist["signed_median"], "p10": dist["p10"], "p90": dist["p90"],
                               "within_1pct": dist["within_1pct"], "within_2pct": dist["within_2pct"]})
                out["Weighted"][pairing][view] = record
    return out


def markdown_report(record, comparisons, title):
    free, strong, ref_pma = record["Free"], record["StrongestSources"], record["ReferencePMA"]
    lines = [f"# {title}", "",
             f"Sources: {len(free) + len(record['ZeroTrace'])} raw, {len(free)} free (contract ZeroTraceIndices excluded). "
             "Relative offset = (run - reference)/|reference|.", "",
             "## Signed distributions of the per-source (diagonal) relative offsets", "",
             "| pairing | quantity | view | n | signed median | p10 | p90 | min | max | |median| | |max| (source) | negative | <1% | <2% | <5% | <10% |",
             "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for pairing, by_key in record["Distributions"].items():
        for key_view, d in by_key.items():
            if d["n"]:
                key, view = key_view.split(":")
                lines.append(f"| {pairing} | {key} | {view} | {d['n']} | {pc(d['signed_median'])} | {pc(d['p10'])} | {pc(d['p90'])} | "
                             f"{pc(d['min'])} | {pc(d['max'])} | {pc(d['abs_median'])} | {pc(d['abs_max'])} ({d['worst_source']}) | "
                             f"{d['negative']} | {d['within_1pct']} | {d['within_2pct']} | {d['within_5pct']} | {d['within_10pct']} |")
    if strong:
        lines += ["", f"## Magnitude-weighted view of p_MA (weights = reference p_MA; free view and the {len(strong)} strongest-MA sources)", "",
                  f"Strongest-MA sources (reference p_MA descending): {strong}; reference p_MA range "
                  f"{ref_pma[str(strong[0])]:.3e} .. {ref_pma[str(strong[-1])]:.3e}; weakest free source "
                  f"{min(ref_pma[str(i)] for i in free):.3e}.", "",
                  "| pairing | view | n | weighted mean | weighted median | unweighted median | p10 | p90 | <1% | <2% |",
                  "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for pairing, by_view in record["Weighted"].items():
            for view, w in by_view.items():
                lines.append(f"| {pairing} | {view} | {w['n']} | {pc(w['weighted_mean'])} | {pc(w['weighted_median'])} | "
                             f"{pc(w['unweighted_median'])} | {pc(w['p10'])} | {pc(w['p90'])} | {w['within_1pct']} | {w['within_2pct']} |")
    lines += ["", "## Per-source table (free view)", "",
              "| source | ref p_MA | " + " | ".join(f"{p} {q}" for p in comparisons for q in OFFSET_KEYS) + " |",
              "|---:|---:|" + "---:|" * (len(comparisons) * len(OFFSET_KEYS))]
    for i in free:
        lines.append(f"| {i} | {ref_pma[str(i)]:.2e} | " + " | ".join(
            pc(comparisons[p][str(i)].get(q)) if str(i) in comparisons[p] else "n/a" for p in comparisons for q in OFFSET_KEYS) + " |")
    return "\n".join(lines) + "\n"


def write_offsets(comparisons, main_pairing, zero_trace, reference_dir, out_md, out_json, *, title,
                  strongest=DEFAULT_STRONGEST):
    ref_pma = reference_p_ma(reference_dir)
    record = offsets(comparisons, main_pairing, zero_trace, ref_pma, strongest=strongest)
    Path(out_md).write_text(markdown_report(record, comparisons, title))
    Path(out_json).write_text(json.dumps(record, indent=2) + "\n")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference-dir", required=True, help="reducer directory of the reference matrices")
    parser.add_argument("--out-md", required=True)
    parser.add_argument("--out-json", required=True)
    parser.add_argument("--zero-trace", type=int, action="append", default=[])
    parser.add_argument("--strongest", type=int, default=DEFAULT_STRONGEST)
    parser.add_argument("--title", default="Offsets vs the reference")
    parser.add_argument("comparisons", nargs="+", metavar="label=comparison.json",
                        help="the first pairing is the main one (run vs reference)")
    args = parser.parse_args(argv)
    comparisons = {item.split("=", 1)[0]: json.loads(Path(item.split("=", 1)[1]).read_text())["PerSource"]
                   for item in args.comparisons}
    record = write_offsets(comparisons, next(iter(comparisons)), args.zero_trace, args.reference_dir, args.out_md,
                           args.out_json, title=args.title, strongest=args.strongest)
    print(json.dumps(record["Weighted"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
