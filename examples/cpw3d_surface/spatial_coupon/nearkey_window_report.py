#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The WINDOW side of the near-key reuse (DESIGN v2 section 4.4, phase 1 - no C++): the reused share per
Type and the reuse error budget beside the D4 figures, from a run's `surface-response-model-energy.csv`
(the fixed-trace evaluation 0, every source summed), its `palace.json` ModelCatalog (model index -> Name),
its `config_resolved.json` (interface index -> class) and the library JSON (the reused models' recorded
bounds).

    share_T  = SUM_reused E_T(model) / E_ref,T
    budget_T = SUM_reused Bound_T(model) x E_T(model) / E_ref,T          (points of the window; a bound: signs may cancel)
    MA       : the budget of record is on the MA_sharp bound (decision 422): Bound_MA_sharp(model) x E_MA(model) / E_ref,MA
               (the raw-MA share, conservative by the factor (1 + t_donor) ~ 1.02); the raw-MA budget beside as information

E_ref,T = the reference energies given (--reference-energies JSON {SA, MS, MA} in J, the window's reference
of record); without them the run's own fixed-trace totals (the sum of every model row + nothing else is NOT
the window total, so the share is then reported against the summed MODEL energies and labelled so). A
window whose reuse budget exceeds 1 point on any Type is flagged ReuseBudgetLarge (the USER's call whether
it gates; DESIGN 4.4).

    nearkey_window_report.py --library LIB.json --postpro DIR [--reference-energies REF.json] --output REPORT.json
"""
import argparse
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(HERE / "qualify") not in sys.path:
    sys.path.insert(0, str(HERE / "qualify"))
import nearkey_detection as detection  # noqa: E402
import spatial_qualification  # noqa: E402

TYPES = ("SA", "MS", "MA")
EVALUATION_FIXED_TRACE = 0
REUSE_BUDGET_LARGE_POINTS = 1.0


class WindowReportError(ValueError):
    """An input the report cannot read (fail closed)."""


def rows(path):
    with open(path, newline="") as handle:
        return [{k.strip(): v.strip() for k, v in r.items()} for r in csv.DictReader(handle)]


def model_energies(postpro, classes, evaluation=EVALUATION_FIXED_TRACE):
    """{model index: {Type: energy (J)}} of the fixed-trace evaluation, every source summed; the
    interface columns `fabricated surface energy[k] (J)` mapped to their class by `classes`."""
    table = rows(Path(postpro) / "surface-response-model-energy.csv")
    if not table:
        raise WindowReportError(f"{postpro}: empty surface-response-model-energy.csv")
    columns = {}
    for name in table[0]:
        if name.startswith("fabricated surface energy["):
            k = int(name[len("fabricated surface energy["):].split("]")[0])
            columns[k] = name
    energies = {}
    for r in table:
        if int(float(r["evaluation"])) != evaluation:
            continue
        index = int(float(r["model"]))
        entry = energies.setdefault(index, {T: 0.0 for T in TYPES})
        for k, name in columns.items():
            cls = classes.get(k)
            if cls in TYPES:
                entry[cls] += float(r[name])
    return energies


def model_names(postpro):
    """{model index: Name} from palace.json SurfaceResponse.ModelCatalog."""
    record = json.loads((Path(postpro) / "palace.json").read_text())
    catalog = record.get("SurfaceResponse", {}).get("ModelCatalog")
    if not catalog:
        raise WindowReportError(f"{postpro}: palace.json carries no SurfaceResponse.ModelCatalog")
    return {int(m["Index"]): m["Name"] for m in catalog}


def interface_classes(postpro):
    config = json.loads((Path(postpro) / "config_resolved.json").read_text())
    return spatial_qualification.interface_classes(config)


def reused_models(library):
    return [m for m in library.get("Models", []) if m.get("QualificationStatus") == detection.STATUS_REUSED]


def window_report(library, postpro, reference_energies=None):
    """The reused share / budget record of one run (see the module docstring)."""
    classes = interface_classes(postpro)
    energies = model_energies(postpro, classes)
    names = model_names(postpro)
    reused = {m["Name"]: m for m in reused_models(library)}
    totals = {T: sum(e[T] for e in energies.values()) for T in TYPES}
    if reference_energies is not None:
        for T in TYPES:
            if T not in reference_energies or float(reference_energies[T]) <= 0:
                raise WindowReportError(f"--reference-energies needs a positive {T} energy (J)")
        reference = {T: float(reference_energies[T]) for T in TYPES}
        convention = "A: E_T(model) / E_ref,T with the reference energies given (the window's reference of record)"
    else:
        reference = totals
        convention = "summed-model: E_T(model) / SUM_models E_T (no reference energies given; NOT the window total)"
    per_model = []
    share = {T: 0.0 for T in TYPES}
    budget = {T: 0.0 for T in TYPES}
    budget["MA_sharp"] = 0.0
    for index, name in sorted(names.items()):
        if name not in reused:
            continue
        model = reused[name]
        prediction = model.get("PredictedReuseError")
        if not prediction:
            raise WindowReportError(f"reused model {name} carries no PredictedReuseError record")
        e = energies.get(index, {T: 0.0 for T in TYPES})
        entry = {"Index": index, "Model": name, "ReuseMode": model.get("ReuseMode"), "Donor": (model.get("ReusedFrom") or {}).get("Donor"),
                 "W": (model.get("NearKey") or {}).get("W"), "Energy": e, "Share": {}, "Bound": {}, "BudgetPoints": {}}
        for T in TYPES:
            s = e[T] / reference[T] if reference[T] > 0 else 0.0
            entry["Share"][T] = s
            entry["Bound"][T] = prediction[T]["Bound"]
            entry["BudgetPoints"][T] = 100.0 * prediction[T]["Bound"] * s
            share[T] += s
            budget[T] += entry["BudgetPoints"][T]
        entry["Bound"]["MA_sharp"] = prediction["MA_sharp"]["Bound"]
        entry["BudgetPoints"]["MA_sharp"] = 100.0 * prediction["MA_sharp"]["Bound"] * entry["Share"]["MA"]
        budget["MA_sharp"] += entry["BudgetPoints"]["MA_sharp"]
        per_model.append(entry)
    flagged = [T for T in ("SA", "MS", "MA_sharp") if budget[T] > REUSE_BUDGET_LARGE_POINTS]
    return {"Postpro": str(postpro), "Library": library.get("Name"), "RuleVersion": (library.get("NearKeyReuse") or {}).get("RuleVersion"),
            "Evaluation": "fixed trace (0)", "InterfaceClasses": classes, "ShareConvention": convention,
            "ReferenceEnergies": reference, "ModelTotals": totals, "ReusedModels": per_model,
            "ReusedShare": share, "ReuseBudgetPoints": budget,
            "BudgetOfRecord": {"SA": budget["SA"], "MS": budget["MS"], "MA_sharp": budget["MA_sharp"], "MA_raw_information": budget["MA"]},
            "ReuseBudgetLarge": flagged, "Counts": {"Reused": len(per_model), "Models": len(names)},
            "Rule": "DESIGN v2 4.4 phase 1: share_T = SUM_reused E_T / E_ref,T; budget_T = SUM_reused Bound_T x E_T / E_ref,T (points; a bound, "
                    "signs may cancel); MA of record on the MA_sharp bound x the raw-MA share (conservative by (1 + t_donor)); "
                    f"ReuseBudgetLarge above {REUSE_BUDGET_LARGE_POINTS} point"}


def markdown(report):
    lines = ["| reused model | mode | donor | W % | share SA / MS / MA % | bound SA / MS / MA / MA_sharp % | budget SA / MS / MA_sharp (MA raw) points |",
             "|---|---|---|---|---|---|---|"]
    for m in report["ReusedModels"]:
        lines.append(f"| {m['Model']} | {m['ReuseMode']} | {m['Donor']} | {100 * (m['W'] or 0):+.3f} | "
                     + " / ".join(f"{100 * m['Share'][T]:.2f}" for T in TYPES) + " | "
                     + " / ".join(f"{100 * m['Bound'][T]:.3f}" for T in TYPES) + f" / {100 * m['Bound']['MA_sharp']:.3f} | "
                     + " / ".join(f"{m['BudgetPoints'][T]:.3f}" for T in ("SA", "MS", "MA_sharp")) + f" ({m['BudgetPoints']['MA']:.3f}) |")
    b = report["ReuseBudgetPoints"]
    s = report["ReusedShare"]
    lines.append(f"| **window** | | | | **{100 * s['SA']:.2f} / {100 * s['MS']:.2f} / {100 * s['MA']:.2f}** | | "
                 f"**{b['SA']:.3f} / {b['MS']:.3f} / {b['MA_sharp']:.3f} ({b['MA']:.3f})** |")
    lines.append(f"\nShare convention: {report['ShareConvention']}. ReuseBudgetLarge: {report['ReuseBudgetLarge'] or 'none'}.")
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--library", required=True, help="the process library the run used (its reused models carry the bounds)")
    parser.add_argument("--postpro", required=True, help="the run's postpro directory")
    parser.add_argument("--reference-energies", help="JSON {SA, MS, MA} reference interface energies (J) of the window")
    parser.add_argument("--output", required=True, help="the report JSON (a .md beside it)")
    args = parser.parse_args(argv)
    try:
        library = json.loads(Path(args.library).read_text())
        reference = json.loads(Path(args.reference_energies).read_text()) if args.reference_energies else None
        report = window_report(library, args.postpro, reference)
    except (WindowReportError, OSError, json.JSONDecodeError, KeyError) as error:
        print(f"NEARKEY_WINDOW_REPORT_FAILED: {error}", file=sys.stderr)
        return 1
    Path(args.output).write_text(json.dumps(report, indent=1) + "\n")
    Path(args.output).with_suffix(".md").write_text(markdown(report))
    print(markdown(report))
    return 0


if __name__ == "__main__":
    sys.exit(main())
