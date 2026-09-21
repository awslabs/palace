#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Compare Palace response matrices (domain + surface) of a run against a reference
(four-edge-physics-11 compare_matrices.py as an importable module).

Common (basis_i, basis_j) pairs only.  Emits a long CSV of every compared entry, a
Markdown summary (diagonal / participation / off-diagonal views) and a JSON summary
with the per-source relative offsets the class statistics and gates are evaluated
on.  The surface matrix's interface indices are labeled by the run config's
Postprocessing.Dielectric map (never assumed; a participation p_X sums every interface
of type X).  Zero-trace sources (diagonal domain energy of exactly 0 in either input, plus the
contract's ZeroTraceIndices: basis knots PEC-constrained in the library whose 3D
coupon energies are nonetheless nonzero) are reported in the raw view and excluded
from the free view.

usage: compare_matrices.py --reference DIR --run DIR --out-prefix PATH --config RUN_CONFIG [--zero-trace I ...]
       [--contract basis-contract.json] [--locations source-locations.csv]
       [--reference-label TEXT] [--run-label TEXT] [--max-entry-rows N]
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics

SURFACE_QUANTITIES = ["Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)", "Q_total_ij (J)",
                      "Q_total_ij normal (J)", "Q_total_ij tangential (J)"]
INTERFACE_TYPES = ("MA", "MS", "SA")
OFFSET_KEYS = ("E_rel", "p_MA_rel", "p_MS_rel", "p_SA_rel")


def interface_names_of(interface_types):
    """Interface index -> type from a config's map ({index: type}, case_inputs.interface_types);
    None (unknown) fails closed: the labels of the surface matrix rows are never assumed."""
    if interface_types is None:
        raise ValueError("the interface index -> type map of the run config is required (never assumed)")
    names = {int(index): str(name) for index, name in interface_types.items()}
    unknown = sorted(set(names.values()) - set(INTERFACE_TYPES))
    if unknown:
        raise ValueError(f"interface types {unknown} outside {INTERFACE_TYPES}")
    return names


def type_energy(surface, diagonal_keys, names, interface_type, i, quantity="Q_total_ij (J)"):
    """The whole-interface diagonal energy of source i summed over every interface of
    the type (None when the type has no interface entry for i)."""
    values = [surface[diagonal_keys[(interface, i)]][quantity] for interface, name in names.items()
              if name == interface_type and (interface, i) in diagonal_keys]
    return sum(values) if values else None


def rows(path):
    with open(path, newline="") as stream:
        return [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(stream)]


def load(directory):
    """(domain {(i, j): Q_ij}, surface {(interface, edge, R, i, j): {quantity: value}})."""
    directory = Path(directory)
    domain = {}
    for row in rows(directory / "domain-response-matrix.csv"):
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        domain[(min(i, j), max(i, j))] = float(row["Q_ij (J)"])
    surface = {}
    for row in rows(directory / "surface-response-matrix.csv"):
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        key = (int(float(row["interface"])), int(float(row["edge"])), float(row["R (m)"]), min(i, j), max(i, j))
        surface[key] = {quantity: float(row[quantity]) for quantity in SURFACE_QUANTITIES}
    return domain, surface


def diagonal_surface_keys(surface):
    """(interface, i) -> the whole-interface diagonal entry key (the smallest edge /
    largest radius group of the interface), independent of the reference's edge / radius
    vocabulary."""
    keys = {}
    for key in surface:
        interface, edge, radius, i, j = key
        if i != j:
            continue
        current = keys.get((interface, i))
        if current is None or (edge, -radius) < (current[1], -current[2]):
            keys[(interface, i)] = key
    return keys


def rel(a, b):
    return (b - a) / abs(a) if a != 0 else (math.nan if b == 0 else math.inf)


def summarize(values):
    values = sorted(abs(v) for v in values if math.isfinite(v))
    if not values:
        return {"n": 0}

    def pct(q):
        return values[min(len(values) - 1, int(math.ceil(q * len(values))) - 1)] if len(values) > 1 else values[0]
    return {"n": len(values), "median": statistics.median(values), "p90": pct(0.9), "p99": pct(0.99), "max": values[-1]}


def fmt(x):
    return "n/a" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.3e}"


def compare(reference_dir, run_dir, *, zero_trace_indices=(), locations=None, interface_types=None,
            reference_interface_types=None):
    """Every compared entry (out_rows), the per-source offsets and the source lists;
    `interface_types` = the run config's interface index -> type map,
    `reference_interface_types` the reference config's (default: the run's; a radial-
    shell run labels many MA interfaces the reference labels as one - the per-type sums
    compare, the per-entry rows cover the common indices only)."""
    names = interface_names_of(interface_types)
    reference_names = interface_names_of(reference_interface_types) if reference_interface_types is not None else names
    ref_domain, ref_surface = load(reference_dir)
    run_domain, run_surface = load(run_dir)
    unlabeled = sorted({key[0] for key in run_surface} - set(names))
    if unlabeled:
        raise ValueError(f"surface matrix interfaces {unlabeled} are not in the config's interface map {names}")
    unlabeled = sorted({key[0] for key in ref_surface} - set(reference_names))
    if unlabeled:
        raise ValueError(f"reference surface matrix interfaces {unlabeled} are not in the reference config's interface map "
                         f"{reference_names}")
    common_names = {index: name for index, name in names.items() if reference_names.get(index) == name}
    common = sorted(set(ref_domain) & set(run_domain))
    sources = sorted({i for i, _ in common} | {j for _, j in common})
    contract_zero = {int(i) for i in zero_trace_indices}
    zero_trace = sorted(i for i in sources
                        if ref_domain.get((i, i), 0.0) == 0.0 or run_domain.get((i, i), 0.0) == 0.0 or i in contract_zero)
    free = [i for i in sources if i not in zero_trace]
    locations = locations or {}
    diagonal_keys = diagonal_surface_keys(ref_surface)
    run_diagonal_keys = diagonal_surface_keys(run_surface)
    out_rows = []

    def add(kind, interface, quantity, i, j, r, v):
        scale = None
        if kind == "domain":
            dii, djj = ref_domain[(i, i)], ref_domain[(j, j)]
            scale = math.sqrt(abs(dii * djj)) if dii and djj else None
        else:
            key_i = diagonal_keys.get((interface, i))
            key_j = diagonal_keys.get((interface, j))
            if key_i is not None and key_j is not None:
                qi, qj = ref_surface[key_i][quantity], ref_surface[key_j][quantity]
                scale = math.sqrt(abs(qi * qj)) if qi and qj else None
        out_rows.append({"kind": kind, "interface": interface or "", "name": names.get(interface, "domain"),
                         "quantity": quantity, "basis_i": i, "basis_j": j, "diagonal": int(i == j),
                         "zero_trace": int(i in zero_trace or j in zero_trace), "reference": r, "run": v,
                         "rel_diff": rel(r, v), "scaled_diff": (v - r) / scale if scale else math.nan})

    for (i, j) in common:
        add("domain", None, "Q_ij (J)", i, j, ref_domain[(i, j)], run_domain[(i, j)])
    for key in sorted(set(ref_surface) & set(run_surface)):
        interface, edge, radius, i, j = key
        if interface not in common_names:
            continue
        for quantity in SURFACE_QUANTITIES:
            add("surface", interface, quantity, i, j, ref_surface[key][quantity], run_surface[key][quantity])
    per_source = {}
    for i in sources:
        er, ev = ref_domain[(i, i)], run_domain[(i, i)]
        per_source[i] = {"E_rel": rel(er, ev), "E_ref": er, "E_run": ev}
        for name in INTERFACE_TYPES:
            qr = type_energy(ref_surface, diagonal_keys, reference_names, name, i)
            qv = type_energy(run_surface, run_diagonal_keys, names, name, i)
            if qr is not None and qv is not None and er and ev:
                pr = qr / er
                pv = qv / ev
                per_source[i][f"p_{name}_rel"] = rel(pr, pv)
                per_source[i][f"p_{name}_ref"] = pr
                per_source[i][f"p_{name}_run"] = pv
    return {"Rows": out_rows, "PerSource": per_source, "Sources": sources, "ZeroTrace": zero_trace, "Free": free,
            "ContractZeroTrace": sorted(i for i in contract_zero if i in sources), "Locations": locations,
            "ReferenceSurface": ref_surface, "RunSurface": run_surface, "InterfaceNames": common_names,
            "RunInterfaceNames": names, "ReferenceInterfaceNames": reference_names}


def summary_record(comparison):
    """The JSON summary: source lists, distribution summaries and the per-source offsets."""
    out_rows, per_source = comparison["Rows"], comparison["PerSource"]
    return {"Sources": comparison["Sources"], "ZeroTrace": comparison["ZeroTrace"], "Free": comparison["Free"],
            "Rows": len(out_rows),
            "DomainDiagonalMaxRel": summarize([r["rel_diff"] for r in out_rows if r["kind"] == "domain" and r["diagonal"]]).get("max"),
            "DomainDiagonalRelFree": summarize([r["rel_diff"] for r in out_rows if r["kind"] == "domain" and r["diagonal"] and not r["zero_trace"]]),
            "DomainDiagonalRelRaw": summarize([r["rel_diff"] for r in out_rows if r["kind"] == "domain" and r["diagonal"]]),
            "DomainOffDiagonalMaxScaled": summarize([r["scaled_diff"] for r in out_rows if r["kind"] == "domain" and not r["diagonal"]]).get("max"),
            "DomainOffDiagonalScaledFree": summarize([r["scaled_diff"] for r in out_rows if r["kind"] == "domain" and not r["diagonal"] and not r["zero_trace"]]),
            "PerSource": {str(i): {key: value for key, value in record.items() if key.endswith("_rel")}
                          for i, record in per_source.items()}}


def markdown_report(comparison, *, reference_label, run_label, max_entry_rows=0):
    out_rows, per_source = comparison["Rows"], comparison["PerSource"]
    sources, zero_trace, free, locations = (comparison["Sources"], comparison["ZeroTrace"], comparison["Free"],
                                            comparison["Locations"])
    contract_zero = comparison["ContractZeroTrace"]
    names = comparison["InterfaceNames"]
    common_pairs = sum(1 for r in out_rows if r["kind"] == "domain")
    lines = [f"# Response-matrix comparison: {run_label} vs {reference_label}", "",
             f"Common sources: {len(sources)} ({sources[0]}..{sources[-1]}; {common_pairs} common domain pairs). "
             f"Zero-trace sources (excluded from the free view): {zero_trace or 'none'}"
             + (f" - {len(contract_zero)} of them from the contract ZeroTraceIndices (PEC-constrained library knots; "
                f"their 3D energies are nonzero and are shown in the raw/diagonal tables)." if contract_zero else "."),
             "Relative difference = (run - reference)/|reference|; scaled difference = (run - reference)/sqrt(|Q_ii Q_jj|) of the reference.", ""]
    lines += ["## Per-source (diagonal) domain energy and participations (raw view: all common sources; ZT = contract zero-trace knot)", "",
              "| source | ZT | location | E_ref (J) | E_run (J) | rel diff | "
              + " | ".join(f"p_{n} ref | p_{n} run | rel diff" for n in INTERFACE_TYPES) + " |",
              "|---:|---|---|---:|---:|---:|" + "---:|" * 9]
    for i in sources:
        record = per_source[i]
        loc = locations.get(i)
        loc_text = f"({loc['x']}, {loc['y']}, {loc['z']}) {loc['lateral']}; {loc['adjacent_surface_at_z0']}" if loc else ""
        cells = [str(i), "ZT" if i in zero_trace else "", loc_text, fmt(record["E_ref"]), fmt(record["E_run"]), fmt(record["E_rel"])]
        for name in INTERFACE_TYPES:
            if f"p_{name}_rel" in record:
                cells += [fmt(record[f"p_{name}_ref"]), fmt(record[f"p_{name}_run"]), fmt(record[f"p_{name}_rel"])]
            else:
                cells += ["n/a"] * 3
        lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "## Distribution of per-source relative differences (|rel diff|)", "",
              "| quantity | view | n | median | p90 | p99 | max | source at max |", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for q in OFFSET_KEYS:
        for view, selection in (("raw", sources), ("free", free)):
            values = {i: per_source[i][q] for i in selection if q in per_source[i] and math.isfinite(per_source[i][q])}
            stats = summarize(values.values())
            worst = max(values, key=lambda k: abs(values[k])) if values else ""
            lines.append(f"| {q} | {view} | {stats['n']} | {fmt(stats.get('median'))} | {fmt(stats.get('p90'))} | "
                         f"{fmt(stats.get('p99'))} | {fmt(stats.get('max'))} | {worst} |")
    lines += ["", "## Sources with the largest deviations (raw view; sorted by |rel diff| of the domain energy)", "",
              "| source | ZT | E rel diff | p_MA rel | p_MS rel | p_SA rel | apex (x, y, z) | z-level | lateral | adjacent surface |",
              "|---:|---|---:|---:|---:|---:|---|---|---|---|"]
    ranked = sorted(sources, key=lambda k: -abs(per_source[k]["E_rel"]) if math.isfinite(per_source[k]["E_rel"]) else 0)[:15]
    for i in ranked:
        loc = locations.get(i, {})
        ps = per_source[i]
        lines.append(f"| {i} | {'ZT' if i in zero_trace else ''} | {fmt(ps['E_rel'])} | {fmt(ps.get('p_MA_rel'))} | "
                     f"{fmt(ps.get('p_MS_rel'))} | {fmt(ps.get('p_SA_rel'))} | ({loc.get('x')}, {loc.get('y')}, {loc.get('z')}) | "
                     f"{loc.get('z_level', '')} | {loc.get('lateral', '')} | {loc.get('adjacent_surface_at_z0', '')} |")
    lines += ["", "Participation p = Q_total_ii(interface)/E_ii(domain) (whole-interface energy incl. the effective layer normalization used by Palace).", ""]
    for view in ("raw", "free"):
        lines += [f"## Diagonal (per-source) surface energies, all surface quantities ({view} view)", "",
                  "| interface | quantity | n | median rel diff | p90 | p99 | max rel diff | source at max | ref magnitude range |",
                  "|---|---|---:|---:|---:|---:|---:|---:|---|"]
        for interface in names:
            for q in SURFACE_QUANTITIES:
                selection = [r for r in out_rows if r["kind"] == "surface" and r["interface"] == interface and r["quantity"] == q
                             and r["diagonal"] and (view == "raw" or not r["zero_trace"])]
                stats = summarize([r["rel_diff"] for r in selection])
                finite = [r for r in selection if math.isfinite(r["rel_diff"])]
                worst = max(finite, key=lambda r: abs(r["rel_diff"]))["basis_i"] if finite else ""
                magnitudes = [abs(r["reference"]) for r in selection]
                lines.append(f"| {names[interface]} ({interface}) | {q} | {stats['n']} | {fmt(stats.get('median'))} | "
                             f"{fmt(stats.get('p90'))} | {fmt(stats.get('p99'))} | {fmt(stats.get('max'))} | {worst} | "
                             f"{fmt(min(magnitudes)) if magnitudes else 'n/a'} .. {fmt(max(magnitudes)) if magnitudes else 'n/a'} |")
        lines.append("")
    for view in ("raw", "free"):
        lines += [f"## Off-diagonal (cross) terms ({view} view)", "",
                  "| kind | quantity | n | median rel diff | p90 rel | max rel diff | median scaled diff | p90 scaled | p99 scaled | max scaled diff | pair at max scaled |",
                  "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
        for kind, interface in [("domain", None)] + [("surface", k) for k in names]:
            for q in (["Q_ij (J)"] if kind == "domain" else SURFACE_QUANTITIES):
                selection = [r for r in out_rows if r["kind"] == kind and (r["interface"] or None) == interface and r["quantity"] == q
                             and not r["diagonal"] and (view == "raw" or not r["zero_trace"])]
                stats = summarize([r["rel_diff"] for r in selection])
                scaled = summarize([r["scaled_diff"] for r in selection])
                finite = [r for r in selection if math.isfinite(r["scaled_diff"])]
                worst = max(finite, key=lambda r: abs(r["scaled_diff"])) if finite else None
                label = "domain" if kind == "domain" else f"{names[interface]} ({interface})"
                lines.append(f"| {label} | {q} | {stats['n']} | {fmt(stats.get('median'))} | {fmt(stats.get('p90'))} | {fmt(stats.get('max'))} | "
                             f"{fmt(scaled.get('median'))} | {fmt(scaled.get('p90'))} | {fmt(scaled.get('p99'))} | {fmt(scaled.get('max'))} | "
                             f"{(worst['basis_i'], worst['basis_j']) if worst else ''} |")
        lines.append("")
    shown = out_rows if max_entry_rows <= 0 else sorted(
        out_rows, key=lambda r: -abs(r["scaled_diff"]) if math.isfinite(r["scaled_diff"]) else 0)[:max_entry_rows]
    lines += ["", f"## {'Every compared entry' if max_entry_rows <= 0 else f'{len(shown)} entries with the largest |scaled diff| (full list in the CSV)'}", "",
              "| kind | quantity | i | j | ZT | reference | run | rel diff | scaled diff |", "|---|---|---:|---:|---|---:|---:|---:|---:|"]
    for r in shown:
        lines.append(f"| {r['name']} | {r['quantity']} | {r['basis_i']} | {r['basis_j']} | {'ZT' if r['zero_trace'] else ''} | "
                     f"{fmt(r['reference'])} | {fmt(r['run'])} | {fmt(r['rel_diff'])} | {fmt(r['scaled_diff'])} |")
    return "\n".join(lines) + "\n"


def read_locations(path):
    return {int(row["index"]): row for row in rows(path)} if path else {}


def write_comparison(reference_dir, run_dir, out_prefix, *, zero_trace_indices=(), locations=None,
                     reference_label="reference", run_label="run", max_entry_rows=0, interface_types=None,
                     reference_interface_types=None):
    """CSV + Markdown + JSON at out_prefix; returns the JSON summary."""
    comparison = compare(reference_dir, run_dir, zero_trace_indices=zero_trace_indices, locations=locations,
                         interface_types=interface_types, reference_interface_types=reference_interface_types)
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    with open(f"{out_prefix}.csv", "w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(comparison["Rows"][0]))
        writer.writeheader()
        writer.writerows(comparison["Rows"])
    Path(f"{out_prefix}.md").write_text(markdown_report(comparison, reference_label=reference_label, run_label=run_label,
                                                        max_entry_rows=max_entry_rows))
    summary = summary_record(comparison)
    Path(f"{out_prefix}.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--reference-label", default="reference")
    parser.add_argument("--run-label", default="run")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--contract", help="basis-contract.json; its ZeroTraceIndices are excluded from the free view")
    parser.add_argument("--zero-trace", type=int, action="append", default=[], help="zero-trace source index (repeatable)")
    parser.add_argument("--locations", help="source-locations.csv from locate_sources.py")
    parser.add_argument("--max-entry-rows", type=int, default=0, help="cap the per-entry Markdown table (0 = all rows)")
    parser.add_argument("--config", required=True, help="the run's Palace config: its Postprocessing.Dielectric entries "
                                                        "label the surface matrix interfaces (index -> type)")
    args = parser.parse_args(argv)
    interface_types = {int(entry["Index"]): entry["Type"]
                       for entry in json.loads(Path(args.config).read_text())["Boundaries"]["Postprocessing"]["Dielectric"]}
    zero = list(args.zero_trace)
    if args.contract:
        zero += [int(i) for i in json.loads(Path(args.contract).read_text()).get("ZeroTraceIndices", [])]
    summary = write_comparison(args.reference, args.run, args.out_prefix, zero_trace_indices=zero,
                               locations=read_locations(args.locations), reference_label=args.reference_label,
                               run_label=args.run_label, max_entry_rows=args.max_entry_rows, interface_types=interface_types)
    print(json.dumps({key: value for key, value in summary.items() if key != "PerSource"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
