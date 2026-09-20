#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Machine-readable qualification gates (qualification-gates.json) evaluated on a
compare_matrices summary, the classify_sources classes, the reference p_MA and the
p-sequence controls.  Every statement records its measured numbers next to the
threshold; the verdict is Passed / Failed / PendingQualification (no reference:
only the p-sequence controls are evaluated, never Passed).

usage: gates.py --comparison run-vs-reference.json --classes source-classes.csv
       --reference-dir DIR --p-sequence p-sequence.json --reference-order P --out PATH
       [--gates qualification-gates.json] [--zero-trace I ...]
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from ma_ms_offsets import reference_p_ma, strongest_sources, weighted  # noqa: E402
from p_sequence import GATED_OBSERVABLES  # noqa: E402

GATES_FILE = HERE / "qualification-gates.json"
VERDICT_PASSED, VERDICT_FAILED, VERDICT_PENDING = "Passed", "Failed", "PendingQualification"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_gates(path=GATES_FILE):
    return json.loads(Path(path).read_text()), sha256(path)


def free_offsets(per_source, key, free):
    return {i: per_source[str(i)][key] for i in free
            if str(i) in per_source and key in per_source[str(i)] and math.isfinite(per_source[str(i)][key])}


def counts_of(values):
    return {"n": len(values), "within_1pct": sum(abs(v) < 0.01 for v in values.values()),
            "within_2pct": sum(abs(v) < 0.02 for v in values.values()),
            "within_5pct": sum(abs(v) < 0.05 for v in values.values()),
            "signed_median": statistics.median(values.values()) if values else None,
            "abs_median": statistics.median(abs(v) for v in values.values()) if values else None,
            "worst_source": max(values, key=lambda i: abs(values[i])) if values else None,
            "worst_offset": (values[max(values, key=lambda i: abs(values[i]))] if values else None)}


def evaluate_energy(gate, per_source, classes, free, wide_classes):
    wide = [i for i in free if classes.get(i) in wide_classes]
    values = free_offsets(per_source, "E_rel", wide)
    all_free = free_offsets(per_source, "E_rel", free)
    bound = gate["MaximumAbsoluteRelativeOffset"]
    failing = sorted(i for i, v in values.items() if abs(v) >= bound)
    return {"Statement": gate["Statement"], "Passed": not failing and bool(values),
            "WideClasses": list(wide_classes), "WideSources": len(values), "WithinBound": len(values) - len(failing),
            "Bound": bound, "FailingSources": failing, "Wide": counts_of(values), "AllFree": counts_of(all_free)}


def evaluate_ma(gate, per_source, free, ref_pma):
    values = free_offsets(per_source, "p_MA_rel", free)
    counts = counts_of(values)
    strong = strongest_sources(ref_pma, [i for i in values if i in ref_pma], gate["StrongestCount"])
    strong_values = {i: values[i] for i in strong}
    weight = weighted(list(values.values()), [ref_pma[i] for i in values]) if values else {"weighted_mean": None}
    checks = {"SignedMedianWithinBound": counts["signed_median"] is not None and abs(counts["signed_median"]) < gate["MaximumAbsoluteSignedMedian"],
              "WeightedMeanWithinBound": weight["weighted_mean"] is not None and abs(weight["weighted_mean"]) < gate["MaximumAbsoluteWeightedMean"],
              "StrongestWithinBound": bool(strong_values) and all(abs(v) < gate["MaximumAbsoluteRelativeOffsetStrongest"] for v in strong_values.values()),
              "FreeWithinBound": bool(values) and all(abs(v) < gate["MaximumAbsoluteRelativeOffsetFree"] for v in values.values())}
    return {"Statement": gate["Statement"], "Passed": all(checks.values()), "Checks": checks, "Free": counts,
            "WeightedMean": weight["weighted_mean"], "StrongestSources": strong, "Strongest": counts_of(strong_values),
            "StrongestFailing": sorted(i for i, v in strong_values.items() if abs(v) >= gate["MaximumAbsoluteRelativeOffsetStrongest"]),
            "FreeFailing": sorted(i for i, v in values.items() if abs(v) >= gate["MaximumAbsoluteRelativeOffsetFree"]),
            "Bounds": {key: gate[key] for key in ("MaximumAbsoluteSignedMedian", "MaximumAbsoluteWeightedMean",
                                                  "MaximumAbsoluteRelativeOffsetStrongest", "MaximumAbsoluteRelativeOffsetFree")}}


def evaluate_ms(gate, per_source, free):
    values = free_offsets(per_source, "p_MS_rel", free)
    counts = counts_of(values)
    checks = {"SignedMedianWithinBound": counts["signed_median"] is not None and abs(counts["signed_median"]) < gate["MaximumAbsoluteSignedMedian"],
              "FreeWithinBound": bool(values) and all(abs(v) < gate["MaximumAbsoluteRelativeOffsetFree"] for v in values.values())}
    return {"Statement": gate["Statement"], "Passed": all(checks.values()), "Checks": checks, "Free": counts,
            "FreeFailing": sorted(i for i, v in values.items() if abs(v) >= gate["MaximumAbsoluteRelativeOffsetFree"]),
            "Bounds": {key: gate[key] for key in ("MaximumAbsoluteSignedMedian", "MaximumAbsoluteRelativeOffsetFree")}}


def evaluate_sa(gate, per_source, free):
    values = free_offsets(per_source, "p_SA_rel", free)
    counts = counts_of(values)
    n = counts["n"]
    fraction2 = counts["within_2pct"] / n if n else 0.0
    fraction5 = counts["within_5pct"] / n if n else 0.0
    checks = {"FractionWithin2Percent": fraction2 >= gate["MinimumFractionWithin2Percent"] and n > 0,
              "FractionWithin5Percent": fraction5 >= gate["MinimumFractionWithin5Percent"] and n > 0}
    return {"Statement": gate["Statement"], "Passed": all(checks.values()), "Checks": checks, "Free": counts,
            "FractionWithin2Percent": fraction2, "FractionWithin5Percent": fraction5,
            "Bounds": {key: gate[key] for key in ("MinimumFractionWithin2Percent", "MinimumFractionWithin5Percent")}}


def evaluate_p_sequence(gate, p_sequence_summary):
    """`p_sequence_summary` = p_sequence.p_sequence output {control: {observable: {seq}}};
    the step of every gated observable towards the higher order (d_high = (high - main)
    / |high|) must be within the bound; d_low and the contraction ratio are reported."""
    controls = {}
    failing = []
    for control, by_observable in p_sequence_summary.items():
        record = {}
        for name in GATED_OBSERVABLES:
            seq = by_observable[name]["seq"]
            bound = gate["MaximumAbsoluteEnergyStep"] if name == "E" else gate["MaximumAbsoluteParticipationStep"]
            step = seq.get("d_high")
            ok = step is not None and math.isfinite(step) and abs(step) < bound
            record[name] = {"StepToHigherOrder": step, "StepFromLowerOrder": seq.get("d_low"), "Bound": bound,
                            "Passed": ok, "r": seq.get("r")}
            if not ok:
                failing.append({"Control": int(control), "Observable": name, "StepToHigherOrder": step, "Bound": bound})
        controls[str(control)] = record
    return {"Statement": gate["Statement"], "Passed": not failing and bool(controls), "Controls": controls, "Failing": failing}


def evaluate(gates, *, comparison=None, classes=None, ref_pma=None, p_sequence_summary=None, reference_order=None,
             free=None, gates_sha256=None):
    """The gate record.  Without `comparison` (no reference matrices) only the p-sequence
    controls are evaluated and the verdict is PendingQualification."""
    record = {"GatesFile": str(GATES_FILE), "GatesSHA256": gates_sha256, "FreeView": gates["FreeView"],
              "ReferenceAnchor": f"vs p{reference_order} anchor" if reference_order is not None else None,
              "Gates": {}}
    table = gates["Gates"]
    if comparison is not None:
        per_source = comparison["PerSource"]
        free = list(free if free is not None else comparison["Free"])
        record["Free"] = free
        record["Gates"]["E"] = evaluate_energy(table["E"], per_source, classes, free, gates["WideClasses"])
        record["Gates"]["p_MA"] = evaluate_ma(table["p_MA"], per_source, free, ref_pma)
        record["Gates"]["p_MS"] = evaluate_ms(table["p_MS"], per_source, free)
        record["Gates"]["p_SA"] = evaluate_sa(table["p_SA"], per_source, free)
    if p_sequence_summary is not None:
        record["Gates"]["PSequenceControls"] = evaluate_p_sequence(table["PSequenceControls"], p_sequence_summary)
    passed = {name: gate["Passed"] for name, gate in record["Gates"].items()}
    record["GatesPassed"] = passed
    if comparison is None:
        record["Verdict"] = VERDICT_PENDING
        record["Reason"] = ("no reference matrices: only the p-sequence controls were evaluated"
                            + ("" if all(passed.values()) else f"; failing {[k for k, v in passed.items() if not v]}"))
    elif all(passed.values()) and passed:
        record["Verdict"] = VERDICT_PASSED
        record["Reason"] = "every gate passed"
    else:
        record["Verdict"] = VERDICT_FAILED
        record["Reason"] = f"failing gates {[k for k, v in passed.items() if not v]}"
    return record


def read_classes(path):
    with open(path, newline="") as stream:
        return {int(row["index"]): row["class"] for row in csv.DictReader(stream)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--comparison", help="compare_matrices JSON of the main run vs the reference (omit when no reference)")
    parser.add_argument("--classes", help="source-classes.csv of classify_sources.py")
    parser.add_argument("--reference-dir", help="reducer directory of the reference matrices")
    parser.add_argument("--p-sequence", help="p_sequence.py JSON")
    parser.add_argument("--reference-order", type=int)
    parser.add_argument("--gates", type=Path, default=GATES_FILE)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    gates, digest = load_gates(args.gates)
    comparison = json.loads(Path(args.comparison).read_text()) if args.comparison else None
    p_sequence_summary = None
    if args.p_sequence:
        loaded = json.loads(Path(args.p_sequence).read_text())
        p_sequence_summary = {int(k): v for k, v in loaded["Sources"].items()}
    record = evaluate(gates, comparison=comparison, classes=read_classes(args.classes) if args.classes else None,
                      ref_pma=reference_p_ma(args.reference_dir) if args.reference_dir else None,
                      p_sequence_summary=p_sequence_summary, reference_order=args.reference_order, gates_sha256=digest)
    args.out.write_text(json.dumps(record, indent=2) + "\n")
    print(record["Verdict"], record["Reason"])
    return 0 if record["Verdict"] == VERDICT_PASSED else 1


if __name__ == "__main__":
    raise SystemExit(main())
