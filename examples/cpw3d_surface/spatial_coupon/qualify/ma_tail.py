#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Sharp-edge MA extrapolation of a radial-shell run (user decision 60(1), supervisor
decision 61a): per source, MA_raw = the sum of the per-ring MA shells, MA_tail = the
remainder of the sharp-edge law inside the innermost ring, MA_sharp = MA_raw + MA_tail.

The tail is the `Consistent` estimator of radial_ma_profile (the decision-56 evidence,
qualify/radial-ma-20260921/RADIAL-MA.md): the top edge (the sharp 90-degree metal edge)
follows the r^(-2/3) law anchored on ring 2 (Theory@2: Q_1 model = Q_2 r_1^(1/3) /
(r_2^(1/3) - r_1^(1/3)), the remainder = model - resolved ring 1), the bottom edge (the
metal / trench edge at the process plane) its own law fitted over rings 2-4 (Fit2-4).
The ring-1 factor (resolved / model) is measured on THIS run at every order, never
carried over from another run.  Per source the record keeps alpha of the top edge
fitted over rings 2-4 with its standard error, the ring-1 factor, MA_raw / MA_tail /
MA_sharp, and the estimator spread (the Fit2-4 and Theory@2 deficits next to the
`Consistent` one: the honest uncertainty of the extrapolation).

A reference is extrapolated by the same rule only when its ring / edge sizing is
recorded (a shelled run: its own shell map); a reference without recorded sizing (the
graded_v2 tet references) is 'reference unextrapolated': its sharp-edge deficit is
MODELLED from its edge size via the eps^(1/3) law - the deficit of an edge integral cut
off at eps scales as eps^(1/3) (decision 55 / HE-CHECK), so per source
deficit_ref = deficit_run(Consistent) x (eps_ref / eps_run)^(1/3) with eps_run the run's
innermost ring radius - and the comparison reports raw-vs-raw next to sharp-vs-modelled.
"""
import csv
import json
import math
from pathlib import Path
import statistics
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import radial_ma_profile as profile  # noqa: E402

ESTIMATOR = profile.HEADLINE_ESTIMATOR
ALPHA_ESTIMATOR = "Fit2-4"
RING1_ESTIMATOR = "Theory@2"
SPREAD_ESTIMATORS = ("Fit2-4", "Theory@2", ESTIMATOR)
EPS_LAW_EXPONENT = 1.0 / 3.0
RULE = ("MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` "
        "estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the "
        "bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy "
        "minus the resolved ring 1; the ring-1 factor is measured on this run at each order")
REFERENCE_MODEL_RULE = ("reference unextrapolated (no ring / edge sizing recorded): its sharp-edge deficit is modelled per "
                        "source as deficit_run(Consistent) x (eps_ref / eps_run)^(1/3) from its edge size eps_ref and the run's "
                        "innermost ring radius eps_run (the eps^(1/3) law of the edge MA inside a cutoff, decision 55)")


def shell_map_of(radial_shells_inputs):
    """{interface index: shell record} of a library-qualification Inputs.RadialShells."""
    return {int(index): value for index, value in radial_shells_inputs["Interfaces"].items()}


def tails(reducer, shell_map, interface_types):
    """Per source (positive domain energy) the extrapolation record of one reducer
    directory: {E, Q_MA_raw, Q_MA_tail, Q_MA_sharp, p_MA_raw, p_MA_sharp, Alpha, AlphaSE,
    Ring1Factor, Deficit {estimator: remainder / Q_MA_raw}, Kinds {top / bottom:
    {Share, Remainder, Estimates}}}."""
    types = {int(index): str(name) for index, name in interface_types.items()}
    ma_interfaces = sorted(index for index, name in types.items() if name == "MA")
    if set(ma_interfaces) != set(shell_map):
        raise ValueError(f"the MA interfaces {ma_interfaces} are not the shell interfaces {sorted(shell_map)}")
    ring_shells = {kind: sorted((index, value) for index, value in shell_map.items() if value["Kind"] == kind)
                   for kind in profile.KINDS}
    energies, diagonal = profile.diagonal_energies(reducer)
    out = {}
    for i in sorted(i for i in energies if energies[i] > 0):
        raw = math.fsum(diagonal[index][i] for index in ma_interfaces)
        entry = {"E": energies[i], "Q_MA_raw": raw, "Kinds": {}}
        remainders = {name: 0.0 for name in SPREAD_ESTIMATORS}
        for kind in profile.KINDS:
            by_ring = {value["Ring"]: (value["InnerRadius"], value["OuterRadius"], diagonal[index][i])
                       for index, value in ring_shells[kind]}
            if not by_ring:
                continue
            fitted = profile.profile_of(by_ring)
            estimates = {name: fitted["Estimates"].get(name) for name in ("Fit2-4", "Theory@2")}
            consistent = fitted["Estimates"].get(profile.kind_estimator(ESTIMATOR, kind))
            kind_record = {"Share": (sum(s[2] for s in by_ring.values()) / raw) if raw else None,
                           "Remainder": consistent["Remainder"] if consistent else 0.0,
                           "Estimator": profile.kind_estimator(ESTIMATOR, kind), "Estimates": estimates}
            entry["Kinds"][kind] = kind_record
            for name in SPREAD_ESTIMATORS:
                estimate = fitted["Estimates"].get(profile.kind_estimator(name, kind))
                if estimate:
                    remainders[name] += estimate["Remainder"]
        top = entry["Kinds"].get("top", {}).get("Estimates", {})
        alpha = top.get(ALPHA_ESTIMATOR)
        ring1 = top.get(RING1_ESTIMATOR)
        tail = remainders[ESTIMATOR]
        entry.update({"Q_MA_tail": tail, "Q_MA_sharp": raw + tail,
                      "p_MA_raw": raw / energies[i], "p_MA_sharp": (raw + tail) / energies[i],
                      "Alpha": alpha["Alpha"] if alpha else None, "AlphaSE": alpha["AlphaSE"] if alpha else None,
                      "Ring1Factor": ring1["Ring1ResolvedOverModel"] if ring1 else None,
                      "Deficit": {name: (value / raw if raw else None) for name, value in remainders.items()}})
        out[i] = entry
    return out


def _median(values):
    values = [v for v in values if v is not None and math.isfinite(v)]
    return statistics.median(values) if values else None


def _quartiles(values):
    values = [v for v in values if v is not None and math.isfinite(v)]
    return [float(np.percentile(values, q)) for q in (25, 75)] if values else None


def _range(values):
    values = [v for v in values if v is not None and math.isfinite(v)]
    return [min(values), max(values)] if values else None


def summary(per_source, sources, *, strongest=(), at=(53, 58)):
    """Distributions over `sources`: alpha (median, quartiles, SE median), the ring-1
    factor (median, range), the deficit of every spread estimator (median, quartiles,
    strongest-median, at the named sources) and the tail share of the p-MA."""
    sources = [i for i in sources if i in per_source]
    strong = [i for i in strongest if i in per_source]
    deficits = {}
    for name in SPREAD_ESTIMATORS:
        values = [per_source[i]["Deficit"][name] for i in sources]
        deficits[name] = {"Median": _median(values), "Quartiles": _quartiles(values),
                          "StrongestMedian": _median([per_source[i]["Deficit"][name] for i in strong]),
                          "At": {str(i): per_source[i]["Deficit"][name] for i in at if i in per_source}}
    alphas = [per_source[i]["Alpha"] for i in sources]
    return {"Sources": len(sources), "Estimator": ESTIMATOR,
            "Alpha": {"Estimator": ALPHA_ESTIMATOR, "Edge": "top", "Median": _median(alphas), "Quartiles": _quartiles(alphas),
                      "SEMedian": _median([per_source[i]["AlphaSE"] for i in sources]),
                      "StrongestMedian": _median([per_source[i]["Alpha"] for i in strong]),
                      "At": {str(i): [per_source[i]["Alpha"], per_source[i]["AlphaSE"]] for i in at if i in per_source},
                      "Theoretical": profile.THEORETICAL_ALPHA},
            "Ring1Factor": {"Estimator": RING1_ESTIMATOR, "Edge": "top",
                            "Median": _median([per_source[i]["Ring1Factor"] for i in sources]),
                            "Range": _range([per_source[i]["Ring1Factor"] for i in sources])},
            "Deficit": deficits,
            "EstimatorSpread": {"Median": [deficits[name]["Median"] for name in ("Fit2-4", ESTIMATOR)],
                                "At": {str(i): [deficits[name]["At"].get(str(i)) for name in ("Fit2-4", ESTIMATOR)] for i in at},
                                "Rule": "Fit2-4 and Consistent (Theory@2 top) deficits: the realised uncertainty of the tail"}}


def modelled_reference(reference_p_ma, run_tails, *, reference_edge_size, run_inner_radius):
    """The reference's modelled sharp-edge p_MA per source (eps^(1/3) law from its edge
    size, the run's Consistent deficit); `reference_edge_size` and `run_inner_radius` in
    the same unit.  Returns {i: {p_MA_raw, p_MA_modelled, ModelledDeficit}}."""
    factor = (reference_edge_size / run_inner_radius) ** EPS_LAW_EXPONENT
    out = {}
    for i, p_raw in reference_p_ma.items():
        tail = run_tails.get(i)
        if tail is None or tail["Deficit"][ESTIMATOR] is None:
            continue
        deficit = tail["Deficit"][ESTIMATOR] * factor
        out[i] = {"p_MA_raw": p_raw, "p_MA_modelled": p_raw * (1.0 + deficit), "ModelledDeficit": deficit}
    return out


def edge_size_text(size, length_unit):
    """The edge size in nm when the unit is um (the coupon unit), else as given."""
    return f"{size * 1000.0:g} nm" if length_unit == "um" else f"{size:g} {length_unit}"


def reference_side(reference_p_ma, run_tails, *, reference_tails=None, reference_edge_size=None, run_inner_radius=None,
                   length_unit="um"):
    """The reference's sharp-edge p_MA per source and the rule applied: its own
    extrapolation when `reference_tails` (a shelled reference) is given, the eps^(1/3)
    model when its edge size is given, else the raw value annotated unmodelled."""
    if reference_tails is not None:
        return {"Rule": "reference extrapolated by the same rule (its own per-ring shells)", "Extrapolated": True,
                "Modelled": False,
                "PerSource": {i: {"p_MA_raw": reference_p_ma[i], "p_MA_sharp": reference_tails[i]["p_MA_sharp"]}
                              for i in reference_p_ma if i in reference_tails}}
    if reference_edge_size is not None:
        modelled = modelled_reference(reference_p_ma, run_tails, reference_edge_size=reference_edge_size,
                                      run_inner_radius=run_inner_radius)
        deficits = [entry["ModelledDeficit"] for entry in modelled.values()]
        median = _median(deficits)
        return {"Rule": REFERENCE_MODEL_RULE, "Extrapolated": False, "Modelled": True,
                "ReferenceEdgeSize": reference_edge_size, "RunInnerRadius": run_inner_radius, "LengthUnit": length_unit,
                "EpsLawExponent": EPS_LAW_EXPONENT, "ModelledDeficitMedian": median,
                "ModelledDeficitQuartiles": _quartiles(deficits),
                "Annotation": (f"reference unextrapolated; sharp-edge deficit modelled at {100 * median:.2f}% (median) from its "
                               f"{edge_size_text(reference_edge_size, length_unit)} edge size via the eps^(1/3) law")
                              if median is not None else None,
                "PerSource": {i: {"p_MA_raw": entry["p_MA_raw"], "p_MA_sharp": entry["p_MA_modelled"],
                                  "ModelledDeficit": entry["ModelledDeficit"]} for i, entry in modelled.items()}}
    return {"Rule": "reference unextrapolated and unmodelled (no ring / edge sizing recorded, no edge size given): its raw "
                    "p_MA stands in for the sharp value", "Extrapolated": False, "Modelled": False,
            "Annotation": "reference unextrapolated; sharp-edge deficit not modelled (no edge size given)",
            "PerSource": {i: {"p_MA_raw": value, "p_MA_sharp": value} for i, value in reference_p_ma.items()}}


def sharp_offsets(per_source, run_tails, reference):
    """Add p_MA_sharp_run / p_MA_sharp_ref / p_MA_sharp_rel (and the tail fields) to a
    compare_matrices per-source record in place; `reference` = reference_side output."""
    side = reference["PerSource"]
    for i, record in per_source.items():
        tail = run_tails.get(i)
        ref = side.get(i)
        if tail is None or ref is None or not record.get("E_run"):
            continue
        record["p_MA_sharp_run"] = tail["p_MA_sharp"]
        record["p_MA_sharp_ref"] = ref["p_MA_sharp"]
        record["p_MA_sharp_rel"] = (tail["p_MA_sharp"] - ref["p_MA_sharp"]) / abs(ref["p_MA_sharp"]) if ref["p_MA_sharp"] else math.nan
        record["MA_tail_run"] = tail["Q_MA_tail"]
        record["MA_deficit_run"] = tail["Deficit"][ESTIMATOR]
    return per_source


def csv_rows(per_source):
    return [{"source": i, "E": e["E"], "Q_MA_raw": e["Q_MA_raw"], "Q_MA_tail": e["Q_MA_tail"], "Q_MA_sharp": e["Q_MA_sharp"],
             "p_MA_raw": e["p_MA_raw"], "p_MA_sharp": e["p_MA_sharp"], "deficit_consistent": e["Deficit"][ESTIMATOR],
             "deficit_fit2_4": e["Deficit"]["Fit2-4"], "deficit_theory_2": e["Deficit"]["Theory@2"],
             "alpha_top_fit2_4": e["Alpha"], "alpha_se": e["AlphaSE"], "ring1_factor_top": e["Ring1Factor"],
             "top_share": (e["Kinds"].get("top") or {}).get("Share"), "bottom_share": (e["Kinds"].get("bottom") or {}).get("Share")}
            for i, e in sorted(per_source.items())]


def write_csv(path, per_source):
    rows = csv_rows(per_source)
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]) if rows else ["source"])
        writer.writeheader()
        writer.writerows(rows)


def pc(x, digits=2):
    return "n/a" if x is None or not math.isfinite(x) else f"{100 * x:+.{digits}f}%"


def f3(x):
    """A fitted value to three decimals; n/a when the fit has no sample (a run without a
    reference has no strongest-source set)."""
    return "n/a" if x is None or not math.isfinite(x) else f"{x:.3f}"


def markdown(record, title):
    lines = [f"# {title}", "", record["Rule"], ""]
    for order, block in record["Orders"].items():
        s = block["Summary"]
        alpha, ring1 = s["Alpha"], s["Ring1Factor"]
        q = alpha["Quartiles"] or [None, None]
        lines += [f"## {order} ({block['Prefix']}; {s['Sources']} sources)", "",
                  f"- alpha (top edge, rings 2-4): median {f3(alpha['Median'])} (quartiles {f3(q[0])} / {f3(q[1])}), se median "
                  f"{f3(alpha['SEMedian'])}; strongest-{len(block['Strongest'])} median {f3(alpha['StrongestMedian'])}; theory {f3(alpha['Theoretical'])}"
                  if alpha["Median"] is not None else "- alpha: no top-edge fit",
                  f"- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median {f3(ring1['Median'])}, range "
                  f"{f3((ring1['Range'] or [None, None])[0])}-{f3((ring1['Range'] or [None, None])[1])}"
                  if ring1["Median"] is not None else "- ring-1 factor: n/a",
                  "", "| deficit estimator | median (quartiles) | strongest median | " + " | ".join(f"at {i}" for i in s["Deficit"][ESTIMATOR]["At"]) + " |",
                  "|---|---|---:|" + "---:|" * len(s["Deficit"][ESTIMATOR]["At"])]
        for name, d in s["Deficit"].items():
            qq = d["Quartiles"] or [None, None]
            lines.append(f"| {name}{' (tail)' if name == ESTIMATOR else ''} | {pc(d['Median'])} ({pc(qq[0])} / {pc(qq[1])}) | "
                         f"{pc(d['StrongestMedian'])} | " + " | ".join(pc(v) for v in d["At"].values()) + " |")
        lines.append("")
    reference = record.get("Reference")
    if reference:
        lines += ["## Reference", "", f"- {reference['Rule']}", f"- {reference.get('Annotation') or 'extrapolated'}", ""]
    return "\n".join(lines) + "\n"


def write_record(path_json, path_md, record, title):
    Path(path_json).write_text(json.dumps(record, indent=2) + "\n")
    Path(path_md).write_text(markdown(record, title))
