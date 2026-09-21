#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Radial MA profile of a radial-shell relabel run (supervisor decision 56): MA per
tube-ring shell per source, the power law of the sharp-edge tail and the extrapolated
remainder inside the innermost ring.

Inputs: a `coupon-library qualify` case root of a case relabeled by
relabel_radial_ma_shells.py (library-qualification.json: Inputs.RadialShells.Interfaces =
interface index -> shell; results/main/<prefix>-p<order>/reducer matrices), the
reference reducer (strongest sources by reference p_MA) and, optionally, the production
run of the parent mesh (the reproduction check: the shells summed per source vs the
parent's whole-MA interface).

Per source, order and edge kind (top = the sharp 90-degree metal edge at the top of the
sidewall, bottom = the metal / trench edge at the process plane), the shell energies
Q_k = Q_total_ii of ring k are the MA integrand integrated over the ring [r_(k-1), r_k)
along the edges.  A power law integrand f(r) = c r^alpha gives Q_k = c L (r_k^(1+alpha) -
r_(k-1)^(1+alpha)) / (1 + alpha) with L the edge length carrying the ring; alpha and c
are fitted by least squares on log Q_k over resolved rings (the innermost ring holds the
singularity a polynomial element cannot represent): over rings 2..K (the whole resolved
profile) and over rings 2..4 (the inner rings, 0.25-3.75 nm, where the edge asymptotics
hold before the metal thickness and the strip width bend the profile); the slope
uncertainty is the standard error of the fit.  The remainder inside the innermost ring is
c L r_1^(1+alpha) / (1 + alpha) - Q_1 (the model's ring-1 energy minus the resolved one)
for each fit, and for the theoretical alpha = -2/3 anchored on ring 2 alone (c from Q_2:
the most local extrapolation, Q_1 model = Q_2 r_1^(1/3) / (r_2^(1/3) - r_1^(1/3))) and
fitted over rings 2..K.  The extrapolated sharp-edge MA is Q_MA + remainder(top) +
remainder(bottom); the deficit is the remainder over Q_MA.  A profile is a clean power law
when the log residual RMS of the rings-2..K fit is below 0.05 and every local slope
between neighbouring fitted rings lies within 0.25 of alpha.

usage: radial_ma_profile.py --record library-qualification.json --case CASE_ID [--production DIR --production-prefix P]
       --out-json PATH --out-md PATH [--strongest N] [--detail I ...]
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import compare_matrices  # noqa: E402
import ma_ms_offsets  # noqa: E402

THEORETICAL_ALPHA = -2.0 / 3.0
CLEAN_RESIDUAL_RMS = 0.05
CLEAN_LOCAL_SLOPE_TOLERANCE = 0.25
DEFAULT_STRONGEST = 20
KINDS = ("top", "bottom")
# Estimator name -> (first ring, last ring or None = K, fixed alpha or None = free).
ESTIMATORS = {"Fit2-K": (2, None, None), "Fit2-4": (2, 4, None),
              "Theory@2": (2, 2, THEORETICAL_ALPHA), "Theory2-K": (2, None, THEORETICAL_ALPHA)}
PRIMARY_FIT = "Fit2-K"


def rows(path):
    with open(path, newline="") as stream:
        return [{key.strip(): value.strip() for key, value in row.items()} for row in csv.DictReader(stream)]


def diagonal_energies(reducer):
    """(E {i: Q_ii}, surface {interface: {i: Q_total_ii}}) of a reducer directory
    (whole-interface diagonal rows: compare_matrices.diagonal_surface_keys)."""
    energies = {}
    for row in rows(Path(reducer) / "domain-response-matrix.csv"):
        i, j = int(float(row["basis_i"])), int(float(row["basis_j"]))
        if i == j:
            energies[i] = float(row["Q_ij (J)"])
    _, surface = compare_matrices.load(reducer)
    keys = compare_matrices.diagonal_surface_keys(surface)
    diagonal = {}
    for (interface, i), key in keys.items():
        diagonal.setdefault(interface, {})[i] = surface[key]["Q_total_ij (J)"]
    return energies, diagonal


def ring_energy_factor(alpha, inner, outer):
    """(r_k^(1+alpha) - r_(k-1)^(1+alpha)) / (1 + alpha) for alpha != -1."""
    return (outer ** (1.0 + alpha) - inner ** (1.0 + alpha)) / (1.0 + alpha)


def fit_power_law(shells, alpha=None):
    """Least squares of log Q_k = log(c L) + log g_alpha(r_(k-1), r_k) over `shells`
    [(inner, outer, Q)]; alpha free (1-D minimisation of the log residual, standard
    error from the local linearisation) or fixed.  Returns {alpha, alpha_se, log_cL,
    residual_rms, n}."""
    inner = np.asarray([s[0] for s in shells])
    outer = np.asarray([s[1] for s in shells])
    q = np.asarray([s[2] for s in shells])
    if len(q) < (2 if alpha is None else 1) or np.any(q <= 0.0):
        return None
    logq = np.log(q)

    def residuals(a):
        g = np.log(ring_energy_factor(a, inner, outer))
        log_cl = float(np.mean(logq - g))
        return logq - g - log_cl, log_cl

    if alpha is None:
        grid = np.linspace(-0.98, 1.5, 2481)
        costs = np.asarray([float(np.sum(residuals(a)[0] ** 2)) for a in grid])
        best = grid[int(np.argmin(costs))]
        # Golden-section refinement around the grid minimum.
        lo, hi = best - 0.002, best + 0.002
        for _ in range(60):
            m1, m2 = lo + (hi - lo) * 0.381966, hi - (hi - lo) * 0.381966
            if np.sum(residuals(m1)[0] ** 2) < np.sum(residuals(m2)[0] ** 2):
                hi = m2
            else:
                lo = m1
        alpha_hat = 0.5 * (lo + hi)
        free = True
    else:
        alpha_hat = float(alpha)
        free = False
    residual, log_cl = residuals(alpha_hat)
    n = len(q)
    rms = float(math.sqrt(np.mean(residual ** 2)))
    alpha_se = None
    if free and n > 2:
        # d(log g)/d alpha at alpha_hat, centred (the intercept absorbs the mean).
        h = 1e-6
        derivative = (np.log(ring_energy_factor(alpha_hat + h, inner, outer)) - np.log(ring_energy_factor(alpha_hat - h, inner, outer))) / (2 * h)
        derivative = derivative - derivative.mean()
        sigma2 = float(np.sum(residual ** 2)) / (n - 2)
        alpha_se = float(math.sqrt(sigma2 / float(np.sum(derivative ** 2)))) if np.sum(derivative ** 2) > 0 else None
    return {"alpha": float(alpha_hat), "alpha_se": alpha_se, "log_cL": float(log_cl), "residual_rms": rms, "n": int(n),
            "free": free}


def local_slopes(shells):
    """Slopes of log density (Q / (r_k - r_(k-1))) vs log of the interval midpoint between
    neighbouring rings."""
    out = []
    for (i0, o0, q0), (i1, o1, q1) in zip(shells, shells[1:]):
        if q0 <= 0 or q1 <= 0:
            out.append(None)
            continue
        d0, d1 = q0 / (o0 - i0), q1 / (o1 - i1)
        out.append(float(math.log(d1 / d0) / math.log(((o1 + i1) / 2) / ((o0 + i0) / 2))))
    return out


def profile_of(shells_by_ring, *, estimators=ESTIMATORS, primary=PRIMARY_FIT):
    """The fits and remainders of one (source, order, kind) profile; `shells_by_ring` =
    {ring: (inner, outer, Q)} for the rings 1..K.  Estimates[name] = {Alpha, AlphaSE,
    ResidualRMS, Rings, Ring1Model, Ring1Resolved, Remainder, Ring1ResolvedOverModel}."""
    rings = sorted(shells_by_ring)
    ordered = [shells_by_ring[k] for k in rings]
    inner1, outer1, q1 = shells_by_ring[rings[0]]
    result = {"Rings": rings, "Q": [s[2] for s in ordered], "Density": [s[2] / (s[1] - s[0]) for s in ordered],
              "LocalSlopes": local_slopes(ordered), "Estimates": {}}
    for name, (first, last, alpha) in estimators.items():
        last = rings[-1] if last is None else last
        selected = [shells_by_ring[k] for k in rings if first <= k <= last]
        model = fit_power_law(selected, alpha=alpha)
        if model is None or model["alpha"] <= -1.0:
            result["Estimates"][name] = None
            continue
        ring1_model = math.exp(model["log_cL"]) * ring_energy_factor(model["alpha"], 0.0, outer1)
        result["Estimates"][name] = {"Alpha": model["alpha"], "AlphaSE": model["alpha_se"], "ResidualRMS": model["residual_rms"],
                                     "Rings": [first, last], "Ring1Model": ring1_model, "Ring1Resolved": q1,
                                     "Remainder": ring1_model - q1, "Ring1ResolvedOverModel": q1 / ring1_model if ring1_model else None}
    fit = result["Estimates"].get(primary)
    first = estimators[primary][0]
    # The slopes between the fitted rings only (slope index k - 1 joins rings k and k + 1).
    slopes = [s for s in result["LocalSlopes"][max(0, first - 1):] if s is not None]
    result["CleanPowerLaw"] = bool(fit is not None and fit["ResidualRMS"] < CLEAN_RESIDUAL_RMS and
                                   all(abs(s - fit["Alpha"]) <= CLEAN_LOCAL_SLOPE_TOLERANCE for s in slopes))
    return result


def analyze(record_path, case_id, *, production=None, production_prefix=None, strongest=DEFAULT_STRONGEST):
    record = json.loads(Path(record_path).read_text())
    case = next(item for item in record["Cases"] if item["Case"] == case_id)
    shell_map = {int(index): value for index, value in case["Inputs"]["RadialShells"]["Interfaces"].items()}
    interface_types = {int(index): name for index, name in case["Inputs"]["Interfaces"].items()}
    ma_interfaces = sorted(index for index, name in interface_types.items() if name == "MA")
    if set(ma_interfaces) != set(shell_map):
        raise ValueError(f"the MA interfaces {ma_interfaces} are not the shell interfaces {sorted(shell_map)}")
    zero_trace = set(case["Sources"]["ZeroTrace"])
    results = Path(case["Root"]) / "results" / "main"
    mains = [item for item in case["Stages"] if item["Role"] == "main"]
    reference = case["Reference"]
    reference_dir = Path(reference["Results"])
    reference_config = json.loads(Path(reference["Config"]).read_text())
    reference_types = {int(e["Index"]): e["Type"] for e in reference_config["Boundaries"]["Postprocessing"]["Dielectric"]}
    ref_pma = ma_ms_offsets.reference_p_ma(reference_dir, reference_types)
    out = {"Case": case_id, "Record": str(record_path), "TheoreticalAlpha": THEORETICAL_ALPHA,
           "Estimators": {name: list(value) for name, value in ESTIMATORS.items()}, "PrimaryFit": PRIMARY_FIT,
           "CleanRule": {"ResidualRMS": CLEAN_RESIDUAL_RMS, "LocalSlopeTolerance": CLEAN_LOCAL_SLOPE_TOLERANCE},
           "Shells": {str(k): v for k, v in shell_map.items()}, "Orders": {}, "Reproduction": None}
    ring_shells = {kind: sorted((index, value) for index, value in shell_map.items() if value["Kind"] == kind)
                   for kind in KINDS}
    far = [index for index, value in shell_map.items() if value["Kind"] == "far"]
    sources_all = None
    for item in mains:
        reducer = results / item["Prefix"] / "reducer"
        energies, diagonal = diagonal_energies(reducer)
        sources = sorted(i for i in energies if i not in zero_trace and energies[i] > 0)
        sources_all = sources
        strongest_sources = ma_ms_offsets.strongest_sources(ref_pma, [i for i in sources if i in ref_pma], strongest)
        per_source = {}
        for i in sources:
            total = sum(diagonal[index][i] for index in ma_interfaces)
            entry = {"E": energies[i], "Q_MA": total, "p_MA": total / energies[i], "Far": sum(diagonal[index][i] for index in far),
                     "Kinds": {}}
            remainders = {name: 0.0 for name in ESTIMATORS}
            top_remainders = {name: 0.0 for name in ESTIMATORS}
            for kind in KINDS:
                by_ring = {value["Ring"]: (value["InnerRadius"], value["OuterRadius"], diagonal[index][i])
                           for index, value in ring_shells[kind]}
                profile = profile_of(by_ring)
                profile["Share"] = sum(s[2] for s in by_ring.values()) / total if total else None
                entry["Kinds"][kind] = profile
                for name in remainders:
                    if profile["Estimates"].get(name):
                        remainders[name] += profile["Estimates"][name]["Remainder"]
                        if kind == "top":
                            top_remainders[name] += profile["Estimates"][name]["Remainder"]
            entry["Remainder"] = remainders
            entry["ExtrapolatedQ_MA"] = {name: total + value for name, value in remainders.items()}
            entry["Deficit"] = {name: value / total for name, value in remainders.items()}
            entry["DeficitTop"] = {name: value / total for name, value in top_remainders.items()}
            entry["Strongest"] = i in strongest_sources
            per_source[str(i)] = entry

        def alpha_of(i, kind, name):
            estimate = per_source[str(i)]["Kinds"][kind]["Estimates"].get(name)
            return estimate["Alpha"] if estimate else None

        strong = [i for i in strongest_sources]
        summary = {"Alpha": {}, "Deficit": {}, "DeficitTop": {}, "CleanPowerLaw": {}}
        for name in ESTIMATORS:
            for kind in KINDS:
                values = [alpha_of(i, kind, name) for i in sources if alpha_of(i, kind, name) is not None]
                summary["Alpha"][f"{kind}:{name}"] = {
                    "Median": statistics.median(values), "Quartiles": [float(np.percentile(values, q)) for q in (25, 75)],
                    "StrongestMedian": statistics.median([alpha_of(i, kind, name) for i in strong if alpha_of(i, kind, name) is not None]),
                    "SEMedian": statistics.median([per_source[str(i)]["Kinds"][kind]["Estimates"][name]["AlphaSE"] or 0.0 for i in sources])}
            for key in ("Deficit", "DeficitTop"):
                values = [per_source[str(i)][key][name] for i in sources]
                summary[key][name] = {"Median": statistics.median(values),
                                      "Quartiles": [float(np.percentile(values, q)) for q in (25, 75)],
                                      "StrongestMedian": statistics.median([per_source[str(i)][key][name] for i in strong]),
                                      "At": {str(i): per_source[str(i)][key][name] for i in (53, 58) if str(i) in per_source}}
        summary["CleanPowerLaw"] = {kind: sum(per_source[str(i)]["Kinds"][kind]["CleanPowerLaw"] for i in sources) for kind in KINDS}
        summary["TopShareMedian"] = statistics.median(per_source[str(i)]["Kinds"]["top"]["Share"] for i in sources)
        summary["BottomShareMedian"] = statistics.median(per_source[str(i)]["Kinds"]["bottom"]["Share"] for i in sources)
        summary["FarShareMedian"] = statistics.median(per_source[str(i)]["Far"] / per_source[str(i)]["Q_MA"] for i in sources)
        summary["Ring1ShareMedian"] = {kind: statistics.median(per_source[str(i)]["Kinds"][kind]["Q"][0] / per_source[str(i)]["Q_MA"]
                                                               for i in sources) for kind in KINDS}
        summary["LocalSlopeMedians"] = {kind: [statistics.median(per_source[str(i)]["Kinds"][kind]["LocalSlopes"][k] for i in sources)
                                               for k in range(len(per_source[str(sources[0])]["Kinds"][kind]["LocalSlopes"]))]
                                        for kind in KINDS}
        out["Orders"][f"p{item['Order']}"] = {"Prefix": item["Prefix"], "Sources": sources, "Strongest": strongest_sources,
                                             "PerSource": per_source, "Summary": summary}
    out["PStep"] = p_step(out)
    if production is not None:
        out["Reproduction"] = reproduction(out, results, mains, Path(production), production_prefix, ma_interfaces, interface_types)
    return out


def p_step(out):
    """Ring-by-ring p-step between the two lowest and highest main orders: median over the
    sources of Q_k(high) / Q_k(low) - 1 per kind and ring, and of the total MA step."""
    orders = sorted(out["Orders"], key=lambda name: int(name[1:]))
    if len(orders) < 2:
        return None
    low, high = out["Orders"][orders[0]], out["Orders"][orders[-1]]
    sources = [str(i) for i in low["Sources"] if str(i) in high["PerSource"]]
    steps = {kind: [statistics.median(high["PerSource"][i]["Kinds"][kind]["Q"][k] / low["PerSource"][i]["Kinds"][kind]["Q"][k] - 1.0
                                      for i in sources)
                    for k in range(len(low["PerSource"][sources[0]]["Kinds"][kind]["Q"]))] for kind in KINDS}
    steps["far"] = statistics.median(high["PerSource"][i]["Far"] / low["PerSource"][i]["Far"] - 1.0 for i in sources)
    steps["Q_MA"] = statistics.median(high["PerSource"][i]["Q_MA"] / low["PerSource"][i]["Q_MA"] - 1.0 for i in sources)
    steps["Ring1ShareOfStep"] = {kind: statistics.median(
        (high["PerSource"][i]["Kinds"][kind]["Q"][0] - low["PerSource"][i]["Kinds"][kind]["Q"][0]) /
        (high["PerSource"][i]["Q_MA"] - low["PerSource"][i]["Q_MA"]) for i in sources
        if high["PerSource"][i]["Q_MA"] != low["PerSource"][i]["Q_MA"]) for kind in KINDS}
    return {"Orders": [orders[0], orders[-1]], "MedianRelativeStep": steps}


def reproduction(out, results, mains, production, production_prefix, ma_interfaces, interface_types):
    """Per order: the shells summed per source vs the production run's whole-MA interface
    (same mesh geometry, same Tol): max relative differences of E, Q_MA and Q_MS."""
    checks = {}
    for item in mains:
        prefix = f"{production_prefix}-p{item['Order']}"
        production_reducer = production / prefix / "reducer"
        if not (production_reducer / "surface-response-matrix.csv").is_file():
            checks[f"p{item['Order']}"] = {"Available": False}
            continue
        run_e, run_surface = diagonal_energies(results / item["Prefix"] / "reducer")
        prod_e, prod_surface = diagonal_energies(production_reducer)
        ms_run = [index for index, name in interface_types.items() if name == "MS"]
        prod_types = None
        config = production / prefix / "reducer.json"
        if config.is_file():
            prod_types = {int(e["Index"]): e["Type"] for e in json.loads(config.read_text())["Boundaries"]["Postprocessing"]["Dielectric"]}
        prod_ma = [k for k, v in (prod_types or {1: "MA", 2: "MS"}).items() if v == "MA"]
        prod_ms = [k for k, v in (prod_types or {1: "MA", 2: "MS"}).items() if v == "MS"]
        common = sorted(set(run_e) & set(prod_e))
        differences = {"E": [], "Q_MA": [], "Q_MS": []}
        for i in common:
            if prod_e[i] == 0.0:
                continue
            differences["E"].append(abs(run_e[i] - prod_e[i]) / abs(prod_e[i]))
            qm = sum(prod_surface[k][i] for k in prod_ma)
            qs = sum(run_surface[k][i] for k in ma_interfaces)
            differences["Q_MA"].append(abs(qs - qm) / abs(qm) if qm else 0.0)
            ms_p = sum(prod_surface[k][i] for k in prod_ms)
            ms_r = sum(run_surface[k][i] for k in ms_run)
            differences["Q_MS"].append(abs(ms_r - ms_p) / abs(ms_p) if ms_p else 0.0)
        checks[f"p{item['Order']}"] = {"Available": True, "Production": str(production_reducer), "Sources": len(common),
                                      "MaxRelative": {key: max(values) for key, values in differences.items()},
                                      "MedianRelative": {key: statistics.median(values) for key, values in differences.items()}}
    return checks


def pc(x, digits=2):
    return "n/a" if x is None or not math.isfinite(x) else f"{100 * x:+.{digits}f}%"


def markdown(out, *, detail):
    lines = [f"# Radial MA profile of {out['Case']}", ""]
    for order, block in out["Orders"].items():
        summary = block["Summary"]
        n = len(block["Sources"])
        lines += [f"## {order}: {n} free sources (strongest-{len(block['Strongest'])} by reference p_MA)",
                  f"MA shares (median): top-edge rings {100 * summary['TopShareMedian']:.1f}%, bottom-edge rings "
                  f"{100 * summary['BottomShareMedian']:.1f}%, far {100 * summary['FarShareMedian']:.1f}%; ring 1 (0-0.25 nm) top "
                  f"{100 * summary['Ring1ShareMedian']['top']:.1f}%, bottom {100 * summary['Ring1ShareMedian']['bottom']:.1f}%",
                  "local slopes (median over sources), rings 1-2 .. 6-7: top " +
                  ", ".join(f"{v:+.2f}" for v in summary["LocalSlopeMedians"]["top"]) + "; bottom " +
                  ", ".join(f"{v:+.2f}" for v in summary["LocalSlopeMedians"]["bottom"]), "",
                  "| estimator | alpha top median (quartiles) | se | strongest | alpha bottom | deficit median (quartiles) | strongest | 53 | 58 | top-edge part |",
                  "|---|---|---:|---:|---:|---|---:|---:|---:|---:|"]
        for name in out["Estimators"]:
            a = summary["Alpha"][f"top:{name}"]
            b = summary["Alpha"][f"bottom:{name}"]
            d = summary["Deficit"][name]
            t = summary["DeficitTop"][name]
            lines.append(f"| {name} | {a['Median']:+.3f} ({a['Quartiles'][0]:+.3f} / {a['Quartiles'][1]:+.3f}) | {a['SEMedian']:.3f} | "
                         f"{a['StrongestMedian']:+.3f} | {b['Median']:+.3f} | {pc(d['Median'])} ({pc(d['Quartiles'][0])} / {pc(d['Quartiles'][1])}) | "
                         f"{pc(d['StrongestMedian'])} | {pc(d['At'].get('53'))} | {pc(d['At'].get('58'))} | {pc(t['Median'])} |")
        lines += [f"clean power law (Fit2-K residual RMS < {out['CleanRule']['ResidualRMS']}, local slopes within "
                  f"{out['CleanRule']['LocalSlopeTolerance']}): top {summary['CleanPowerLaw']['top']}/{n}, bottom {summary['CleanPowerLaw']['bottom']}/{n}", "",
                  "| source | p_MA | top / bottom / far share | top Q_k / Q_MA, rings 1..7 | slopes 2-3, 3-4, 4-5 | alpha Fit2-K (se) | alpha Fit2-4 | clean | deficit Fit2-K / Fit2-4 / Theory@2 | ring1 resolved / Theory@2 model |",
                  "|---:|---:|---|---|---|---|---:|---|---|---:|"]
        for i in detail:
            entry = block["PerSource"].get(str(i))
            if entry is None:
                continue
            top, bottom = entry["Kinds"]["top"], entry["Kinds"]["bottom"]
            e = top["Estimates"]
            lines.append(f"| {i}{'*' if entry['Strongest'] else ''} | {entry['p_MA']:.3e} | {100 * top['Share']:.0f} / {100 * bottom['Share']:.0f} / "
                         f"{100 * entry['Far'] / entry['Q_MA']:.0f}% | " + " ".join(f"{q / entry['Q_MA']:.3f}" for q in top["Q"]) + " | "
                         + " ".join(f"{v:+.2f}" for v in top["LocalSlopes"][1:4]) + f" | {e['Fit2-K']['Alpha']:+.3f} ({e['Fit2-K']['AlphaSE'] or 0:.3f}) | "
                         f"{e['Fit2-4']['Alpha']:+.3f} | {'yes' if top['CleanPowerLaw'] else 'no'} | "
                         f"{pc(entry['Deficit']['Fit2-K'])} / {pc(entry['Deficit']['Fit2-4'])} / {pc(entry['Deficit']['Theory@2'])} | "
                         f"{e['Theory@2']['Ring1ResolvedOverModel']:.3f} |")
        lines.append("")
    if out.get("PStep"):
        step = out["PStep"]["MedianRelativeStep"]
        lines += [f"## p-step {out['PStep']['Orders'][0]} -> {out['PStep']['Orders'][1]} per ring (median relative change of Q_k)",
                  "top rings 1..7: " + ", ".join(pc(v) for v in step["top"]) + "; bottom rings 1..7: " + ", ".join(pc(v) for v in step["bottom"])
                  + f"; far {pc(step['far'])}; Q_MA {pc(step['Q_MA'])}; ring-1 share of the Q_MA step: top {pc(step['Ring1ShareOfStep']['top'])}, "
                  f"bottom {pc(step['Ring1ShareOfStep']['bottom'])}", ""]
    if out.get("Reproduction"):
        lines.append("## Reproduction of the production run (shells summed per source)")
        for order, check in out["Reproduction"].items():
            if check.get("Available"):
                m = check["MaxRelative"]
                lines.append(f"{order}: {check['Sources']} sources, max |rel| E {m['E']:.2e}, Q_MA {m['Q_MA']:.2e}, Q_MS {m['Q_MS']:.2e}")
            else:
                lines.append(f"{order}: no production stage")
        lines.append("")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--record", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--production", type=Path, help="results/main of the parent mesh's production run")
    parser.add_argument("--production-prefix", help="stage prefix of the production run (default: the parent case id)")
    parser.add_argument("--strongest", type=int, default=DEFAULT_STRONGEST)
    parser.add_argument("--detail", type=int, action="append", default=[], help="sources tabulated (default: the strongest)")
    parser.add_argument("--out-json", type=Path, required=True)
    parser.add_argument("--out-md", type=Path, required=True)
    args = parser.parse_args(argv)
    out = analyze(args.record, args.case, production=args.production, production_prefix=args.production_prefix,
                  strongest=args.strongest)
    detail = args.detail or sorted(next(iter(out["Orders"].values()))["Strongest"])
    args.out_json.write_text(json.dumps(out, indent=2) + "\n")
    args.out_md.write_text(markdown(out, detail=detail))
    print(args.out_md.read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
