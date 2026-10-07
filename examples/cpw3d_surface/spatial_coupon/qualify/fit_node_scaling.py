#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The cost model's NodeScaling calibration (decision 457 (1) / (3)): how a stage's node-used
memory and wall time divide between the nodes of a multi-node job, measured on ONE stored
coupon run at 1 node (the run_stages statuses of record) and at 2 nodes (the multi-node
runner's status: per-node sampled peaks, Palace per-node peaks).

Memory, per stage kind (Worker / Reducer from the response stages; LocalEdge = the ordinary
Palace stage, measured on the dense-twin run through Palace's per-node peak report): the
per-node peak at 2 nodes against the one-node peak gives the REPLICATED fraction r of
  per_node(N) = used(1) x (r + (1 - r) / N)      ->      r = 2 x per_node(2) / used(1) - 1,
clamped to [0, 1] (r < 0 = better than even division, planned as 0; the largest r over the
stages of a kind is kept: fail-closed).  Time, per part: the measured 1 -> 2 node speedup
s and the exponent log2(s) of T(N) = T(1) / N^exponent (WorkerPerSource from the per-source
timings, WorkerNonSource from the worker wall minus its sources' seconds, ReducerSetup /
ReducerReduction from the reducer's elapsed-time report, LocalEdge from the ordinary stage's
Palace total); estimate_stages caps every exponent at 1.

usage: fit_node_scaling.py --one-node-status STATUS.json ... --two-node-status STATUS.json
       [--ordinary-one-node-log palace.log --ordinary-two-node-stage NAME] --out node-scaling.json
"""
import argparse
import json
import math
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from refit_cost_model import REDUCTION_TIMERS, memory_gb, parse_elapsed_time_report  # noqa: E402
from run_stages import parse_log  # noqa: E402

GIB = 2**30
EXPONENT_RULE = ("T(N) = T(1) / N^Exponent with Exponent = log2(the measured 1 -> 2 node speedup of the part); estimate_stages caps "
                 "it at 1 (no superlinear speedup is budgeted); a negative exponent (slower at 2 nodes) is kept")
MEMORY_RULE = ("per-node node-used GiB at N nodes = the one-node figure x (ReplicatedFraction + (1 - ReplicatedFraction) / N); "
               "ReplicatedFraction = 2 x per-node peak at 2 nodes / the one-node peak - 1, clamped to [0, 1], the largest over the "
               "measured stages of the kind")


def stage_kind(name):
    if name.endswith("-reducer"):
        return "Reducer"
    if name.endswith("-worker") or re.search(r"-worker-block\d+$", name):
        return "Worker"
    return "LocalEdge"


def stage_index(statuses):
    index = {}
    for status in statuses:
        for stage in status["Stages"]:
            index[stage["Name"]] = stage
    return index


def worker_parts(stage):
    timings = stage["Parsed"].get("SourceTiming") or []
    total = sum(t["TotalSeconds"] for t in timings)
    return {"Sources": len(timings), "PerSourceSeconds": total / len(timings) if timings else None,
            "NonSourceSeconds": max(stage["WallSeconds"] - total, 0.0), "WallSeconds": stage["WallSeconds"],
            "MeanPCGIterations": sum(t["Iterations"] for t in timings) / len(timings) if timings else None}


def reducer_parts(stage):
    timers = parse_elapsed_time_report(stage["Parsed"]["ElapsedTimeReport"])
    reduction = sum(timers.get(name, 0.0) for name in REDUCTION_TIMERS)
    wall = stage["WallSeconds"]
    fraction = reduction / sum(timers.values())
    return {"ReductionSeconds": wall * fraction, "SetupSeconds": wall * (1.0 - fraction), "WallSeconds": wall}


def per_node_peak_bytes(stage):
    """The largest per-node sampled peak of a stage (a multi-node runner stage carries the
    per-node table; a one-node stage the node's figure)."""
    per_node = stage.get("NodePeakUsedBytesSampledPerNode")
    if per_node:
        return max(per_node.values()), len(per_node)
    return stage["NodePeakUsedBytesSampled"], 1


def replicated_fraction(used_one, per_node_two, nodes_two=2):
    raw = nodes_two * per_node_two / used_one - 1.0
    raw = raw * (1.0 / (nodes_two - 1.0))   # the general-N inversion of r + (1 - r) / N
    return {"Raw": raw, "ReplicatedFraction": min(1.0, max(0.0, raw))}


def exponent(speedup):
    return math.log2(speedup) if speedup > 0 else None


PALACE_GB_PER_GIB = 1.0737   # the cost model's PalaceGBPerGiB (Palace prints GB as 1e9 bytes)


def fit(one_node_statuses, two_node_status, *, ordinary_one_node_log=None, ordinary_two_node_stage=None,
        palace_gb_per_gib=PALACE_GB_PER_GIB):
    one = stage_index(one_node_statuses)
    two = stage_index([two_node_status])
    nodes_two = int(two_node_status.get("Nodes") and len(two_node_status["Nodes"]) or 2)
    memory = {"Worker": [], "Reducer": [], "LocalEdge": []}
    times = {"WorkerPerSource": [], "WorkerNonSource": [], "ReducerSetup": [], "ReducerReduction": [], "LocalEdge": []}
    matched = []
    for name, stage_two in two.items():
        stage_one = one.get(name)
        if stage_one is None or stage_two["State"] != "complete" or stage_one["State"] != "complete":
            continue
        kind = stage_kind(name)
        used_one, _ = per_node_peak_bytes(stage_one)
        per_node_two, node_count = per_node_peak_bytes(stage_two)
        fraction = replicated_fraction(used_one, per_node_two, node_count)
        record = {"Stage": name, "Kind": kind, "OneNodeUsedGiB": used_one / GIB, "TwoNodePerNodeUsedGiB": per_node_two / GIB,
                  "TwoNodePerNode": {host: value / GIB for host, value in (stage_two.get("NodePeakUsedBytesSampledPerNode") or {}).items()},
                  "PalacePeakOneNode": stage_one["Parsed"].get("PalacePeakMemory"), "PalacePeakTwoNodes": stage_two["Parsed"].get("PalacePeakMemory"),
                  **fraction}
        memory[kind].append(record)
        if kind == "Worker":
            parts_one, parts_two = worker_parts(stage_one), worker_parts(stage_two)
            record["Time"] = {"OneNode": parts_one, "TwoNodes": parts_two}
            if parts_one["PerSourceSeconds"] and parts_two["PerSourceSeconds"]:
                times["WorkerPerSource"].append((parts_one["PerSourceSeconds"] / parts_two["PerSourceSeconds"], name))
            if parts_two["NonSourceSeconds"] > 0:
                times["WorkerNonSource"].append((parts_one["NonSourceSeconds"] / parts_two["NonSourceSeconds"], name))
        elif kind == "Reducer":
            parts_one, parts_two = reducer_parts(stage_one), reducer_parts(stage_two)
            record["Time"] = {"OneNode": parts_one, "TwoNodes": parts_two}
            times["ReducerSetup"].append((parts_one["SetupSeconds"] / parts_two["SetupSeconds"], name))
            times["ReducerReduction"].append((parts_one["ReductionSeconds"] / parts_two["ReductionSeconds"], name))
        matched.append(record)
    ordinary = None
    if ordinary_one_node_log is not None and ordinary_two_node_stage is not None:
        parsed_one = parse_log(Path(ordinary_one_node_log).read_text(errors="replace"))
        stage_two = two[ordinary_two_node_stage]
        peak_one = memory_gb(parsed_one["PalacePeakMemory"]["Total"])
        peak_two_max = memory_gb(stage_two["Parsed"]["PalacePeakMemory"]["Max"])
        peak_two_total = memory_gb(stage_two["Parsed"]["PalacePeakMemory"]["Total"])
        palace_based = replicated_fraction(peak_one, peak_two_max, nodes_two)
        # The planner compares NODE-USED figures, whose per-node overhead (OS, MPI, the runner)
        # does not divide: the one-node node-used peak of a run outside the runner (no sample) is
        # reconstructed as its Palace peak in GiB + the measured 2-node per-node overhead
        # (node used - Palace per-node peak); the larger of the two fractions is kept.
        per_node_used = stage_two.get("NodePeakUsedBytesSampledPerNode") or {}
        node_used_based = None
        if per_node_used:
            used_two_gib = max(per_node_used.values()) / GIB
            overhead_gib = used_two_gib - peak_two_max / palace_gb_per_gib
            used_one_gib = peak_one / palace_gb_per_gib + overhead_gib
            node_used_based = {**replicated_fraction(used_one_gib, used_two_gib, nodes_two), "OneNodeUsedGiBReconstructed": used_one_gib,
                               "TwoNodePerNodeUsedGiB": used_two_gib, "PerNodeOverheadGiB": overhead_gib}
        chosen = node_used_based if node_used_based and node_used_based["ReplicatedFraction"] > palace_based["ReplicatedFraction"] else palace_based
        ordinary = {"Stage": ordinary_two_node_stage, "Kind": "LocalEdge",
                    "Basis": ("the larger of the Palace per-node peak fraction and the node-used fraction with the one-node node-used "
                              "reconstructed (Palace peak / PalaceGBPerGiB + the 2-node per-node overhead); the one-node run of record ran "
                              "outside the runner: no node-used sample"),
                    "OneNodePalacePeakGB": peak_one, "TwoNodePalacePeakMaxGB": peak_two_max, "TwoNodePalacePeakTotalGB": peak_two_total,
                    "PalaceBased": palace_based, "NodeUsedBased": node_used_based,
                    "TwoNodePerNodeUsedGiB": {host: value / GIB for host, value in per_node_used.items()},
                    "OneNodePalaceTotalSeconds": parsed_one["PalaceTotalSeconds"], "TwoNodePalaceTotalSeconds": stage_two["Parsed"]["PalaceTotalSeconds"],
                    "OneNodePCG": parsed_one.get("PCG"), "TwoNodePCG": stage_two["Parsed"].get("PCG"),
                    "Raw": chosen["Raw"], "ReplicatedFraction": chosen["ReplicatedFraction"]}
        memory["LocalEdge"].append(ordinary)
        times["LocalEdge"].append((parsed_one["PalaceTotalSeconds"] / stage_two["Parsed"]["PalaceTotalSeconds"], ordinary_two_node_stage))
    block = {"MeasuredNodes": nodes_two, "Rule": {"Memory": MEMORY_RULE, "Time": EXPONENT_RULE}, "Memory": {}, "Time": {}}
    for kind, records in memory.items():
        if records:
            block["Memory"][kind] = {"ReplicatedFraction": round(max(r["ReplicatedFraction"] for r in records), 4),
                                     "SetBy": max(records, key=lambda r: r["ReplicatedFraction"])["Stage"],
                                     "Raw": [round(r["Raw"], 4) for r in records]}
    if "LocalEdge" not in block["Memory"] and "Worker" in block["Memory"]:
        block["Memory"]["LocalEdge"] = {**block["Memory"]["Worker"], "Rule": "carried from Worker: the ordinary stage was not measured"}
    for part, values in times.items():
        if values:
            # The smallest measured speedup of the part (fail-closed on the time side).
            speedup, name = min(values)
            block["Time"][part] = {"Speedup": round(speedup, 4), "Exponent": round(exponent(speedup), 4), "SetBy": name,
                                   "Measured": [{"Speedup": round(s, 4), "Stage": n} for s, n in values]}
    if "LocalEdge" not in block["Time"]:
        block["Time"]["LocalEdge"] = {"Speedup": 1.0, "Exponent": 0.0, "Rule": "not measured: no speedup assumed"}
    block["Measured"] = {"Stages": matched, "Ordinary": ordinary,
                         "TwoNodeJob": {"PBSJobID": two_node_status.get("PBSJobID"), "Hosts": two_node_status.get("Nodes"),
                                        "Preflight": {host: node for host, node in (two_node_status.get("Preflight", {}).get("Nodes") or {}).items()}},
                         "OneNodeJobs": [status.get("PBSJobID") for status in one_node_statuses]}
    return block


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--one-node-status", type=Path, action="append", required=True, help="run_stages status.json of the one-node run of record (repeatable)")
    parser.add_argument("--two-node-status", type=Path, required=True, help="the multi-node runner's status.json")
    parser.add_argument("--ordinary-one-node-log", type=Path, help="palace.log of the ordinary (dense-twin) stage at one node")
    parser.add_argument("--ordinary-two-node-stage", help="the ordinary stage's name in --two-node-status")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    block = fit([json.loads(path.read_text()) for path in args.one_node_status], json.loads(args.two_node_status.read_text()),
                ordinary_one_node_log=args.ordinary_one_node_log, ordinary_two_node_stage=args.ordinary_two_node_stage)
    args.out.write_text(json.dumps(block, indent=2) + "\n")
    print(json.dumps({"Memory": block["Memory"], "Time": {part: {k: v for k, v in value.items() if k != "Measured"} for part, value in block["Time"].items()}}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
