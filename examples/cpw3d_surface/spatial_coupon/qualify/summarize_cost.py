#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Per-source cost metrics, node-hours and full-coupon extrapolation from a
run_stages.py status.json (four-edge-physics-11 summarize_cost.py; the source count,
node count and block size are arguments / read from the status).

Per worker + reducer stage pair: elements, H1 / ND / RT, per-source solve and total
seconds, PCG iterations, seconds per PCG iteration, worker / reducer Palace and wall
seconds, reducer pairs and block pairs, node peak used GiB (sampled), Palace peak, the
max single-rank RSS; the full-coupon estimate at `full_sources` sources; node-hours =
wall seconds x nodes / 3600 (the "node-h per coupon" of the suite table).

usage: summarize_cost.py STATUS.json --out PATH [--nodes N] [--block-size B]
"""
import argparse
import json
from pathlib import Path
import statistics

GIB = 2**30
DEFAULT_BLOCK_SIZE = 6


def stage_cost(stages, prefix, *, full_sources, nodes, block_size=DEFAULT_BLOCK_SIZE):
    worker, reducer = stages[f"{prefix}-worker"], stages[f"{prefix}-reducer"]
    parsed_worker, parsed_reducer = worker["Parsed"], reducer["Parsed"]
    timings = parsed_worker["SourceTiming"]
    per_source_total = [t["TotalSeconds"] for t in timings]
    per_source_solve = [t["SolveSeconds"] for t in timings]
    iterations = [t["Iterations"] for t in timings]
    n = len(timings)
    worker_non_source = parsed_worker["PalaceTotalSeconds"] - sum(per_source_total)
    seconds_per_iteration = sum(per_source_solve) / sum(iterations) if sum(iterations) else None
    pairs = n * (n + 1) // 2
    blocks = -(-n // block_size)
    block_pairs = blocks * (blocks + 1) // 2
    # Reducer setup is everything except the block-pair work; approximate it by the
    # worker's non-source time (same mesh, same operator / linear setup).  The remainder
    # is attributed to the block pairs actually reduced.
    reducer_setup = min(worker_non_source, parsed_reducer["PalaceTotalSeconds"])
    reducer_pair_seconds = max(parsed_reducer["PalaceTotalSeconds"] - reducer_setup, 0.0)
    full_pairs = full_sources * (full_sources + 1) // 2
    mean_total = statistics.mean(per_source_total) if per_source_total else None
    estimate = None
    if mean_total is not None:
        estimate = {"WorkerSeconds": worker_non_source + full_sources * mean_total,
                    "ReducerSecondsPairScaled": reducer_setup + reducer_pair_seconds * full_pairs / pairs}
        estimate["TotalSeconds"] = estimate["WorkerSeconds"] + estimate["ReducerSecondsPairScaled"]
    wall = worker["WallSeconds"] + reducer["WallSeconds"]
    return {"Elements": parsed_worker.get("Elements"), "H1": parsed_worker.get("H1"), "ND": parsed_worker.get("ND"),
            "RT": parsed_worker.get("RT"), "Order": parsed_worker.get("Order"),
            "Sources": [t["Index"] for t in timings], "PCGIterations": iterations,
            "MeanPCGIterations": statistics.mean(iterations) if iterations else None,
            "MaxPCGIterations": max(iterations) if iterations else None,
            "PerSourceSolveSeconds": per_source_solve, "PerSourceTotalSeconds": per_source_total,
            "MeanPerSourceTotalSeconds": mean_total, "SecondsPerPCGIteration": seconds_per_iteration,
            "WorkerPalaceTotalSeconds": parsed_worker["PalaceTotalSeconds"], "WorkerWallSeconds": worker["WallSeconds"],
            "WorkerNonSourceSeconds": worker_non_source,
            "ReducerPalaceTotalSeconds": parsed_reducer["PalaceTotalSeconds"], "ReducerWallSeconds": reducer["WallSeconds"],
            "ReducerPairs": pairs, "ReducerBlockPairs": block_pairs, "ReducerPairSecondsEstimate": reducer_pair_seconds,
            "StageWallSeconds": wall, "NodeHours": wall * nodes / 3600.0,
            "Memory": {"WorkerNodePeakUsedGiBSampled": (worker.get("NodePeakUsedBytesSampled") or 0) / GIB,
                       "ReducerNodePeakUsedGiBSampled": (reducer.get("NodePeakUsedBytesSampled") or 0) / GIB,
                       "WorkerPalacePeakTotal": parsed_worker.get("PalacePeakMemory", {}).get("Total"),
                       "ReducerPalacePeakTotal": parsed_reducer.get("PalacePeakMemory", {}).get("Total"),
                       "WorkerMaxSingleRankRSSGiB": (worker.get("MaxSingleProcessRSSBytes") or 0) / GIB,
                       "ReducerMaxSingleRankRSSGiB": (reducer.get("MaxSingleProcessRSSBytes") or 0) / GIB},
            f"FullCouponEstimate{full_sources}Sources": estimate}


def merge_split_statuses(statuses, jobs):
    """One coupon-level status from the per-job statuses of a split coupon (decision
    61b): `statuses` = job name -> run_stages.py status, `jobs` = the coupon's job
    records (Name, Kind, StageNames).  Every main stage's block workers
    (<prefix>-worker-block<k>) merge into one <prefix>-worker record - wall seconds and
    Palace totals summed (the node time actually spent), the per-source timings and PCG
    counts concatenated in block order, peaks as the largest block's - so stage_cost
    reads the split exactly as a single job's stage; the reducer and the first job's
    control / local-edge stages are copied; TotalSeconds = the sum over the jobs (node
    seconds), and Jobs carries every job's own total."""
    merged = {"Version": 2, "Split": True, "Jobs": {}, "Stages": [], "PBSJobID": None, "Host": None,
              "StartUTC": None, "EndUTC": None, "TotalSeconds": 0.0}
    blocks = {}
    for job in jobs:
        status = statuses[job["Name"]]
        merged["Jobs"][job["Name"]] = {"Kind": job["Kind"], "PBSJobID": status.get("PBSJobID"), "Host": status.get("Host"),
                                       "StartUTC": status.get("StartUTC"), "EndUTC": status.get("EndUTC"),
                                       "TotalSeconds": status.get("TotalSeconds"), "State": status.get("State")}
        merged["TotalSeconds"] += status.get("TotalSeconds") or 0.0
        merged["StartUTC"] = min(filter(None, (merged["StartUTC"], status.get("StartUTC"))), default=None)
        merged["EndUTC"] = max(filter(None, (merged["EndUTC"], status.get("EndUTC"))), default=None)
        for stage in status["Stages"]:
            name = stage["Name"]
            if "-worker-block" in name:
                prefix, block = name.rsplit("-worker-block", 1)
                blocks.setdefault(prefix, []).append((int(block), stage))
            else:
                merged["Stages"].append(stage)
    for prefix, items in blocks.items():
        items.sort(key=lambda pair: pair[0])
        stages = [stage for _, stage in items]
        parsed = {"SourceTiming": [], "PCG": [], "Iterations": [], "Nonconvergence": [], "PairsProgress": [], "BlockPairsProgress": []}
        for stage in stages:
            for key in parsed:
                parsed[key].extend(stage["Parsed"].get(key) or [])
        first = stages[0]["Parsed"]
        for key in ("Order", "H1", "ND", "RT", "Elements"):
            if key in first:
                parsed[key] = first[key]
        parsed["PalaceTotalSeconds"] = sum(stage["Parsed"].get("PalaceTotalSeconds") or 0.0 for stage in stages)
        peaks = [stage["Parsed"].get("PalacePeakMemory") for stage in stages if stage["Parsed"].get("PalacePeakMemory")]
        if peaks:
            parsed["PalacePeakMemory"] = max(peaks, key=lambda peak: _memory_bytes(peak.get("Total")))
        parsed["Blocks"] = [{"Block": block, "Name": stage["Name"], "WallSeconds": stage.get("WallSeconds"),
                             "Sources": [t["Index"] for t in stage["Parsed"].get("SourceTiming") or []],
                             "PalaceTotalSeconds": stage["Parsed"].get("PalaceTotalSeconds"),
                             "PalacePeakTotal": (stage["Parsed"].get("PalacePeakMemory") or {}).get("Total")}
                            for block, stage in items]
        merged["Stages"].append({"Name": f"{prefix}-worker", "State": "complete" if all(s["State"] == "complete" for s in stages) else "incomplete",
                                 "WallSeconds": sum(s.get("WallSeconds") or 0.0 for s in stages),
                                 "NodePeakUsedBytesSampled": max((s.get("NodePeakUsedBytesSampled") or 0) for s in stages),
                                 "MaxSingleProcessRSSBytes": max((s.get("MaxSingleProcessRSSBytes") or 0) for s in stages),
                                 "Blocks": len(stages), "Parsed": parsed})
    merged["State"] = "complete" if all(s["State"] == "complete" for s in merged["Stages"]) else "incomplete"
    return merged


def _memory_bytes(text):
    """'123.4G' -> bytes (the Palace peak-memory report unit suffixes)."""
    if not text:
        return 0.0
    units = {"K": 2**10, "M": 2**20, "G": 2**30, "T": 2**40}
    text = str(text).strip()
    if text[-1] in units:
        return float(text[:-1]) * units[text[-1]]
    return float(text)


def summarize(status, *, full_sources, nodes, block_size=DEFAULT_BLOCK_SIZE):
    stages = {stage["Name"]: stage for stage in status["Stages"]}
    summary = {"PBSJobID": status.get("PBSJobID"), "Host": status.get("Host"), "StartUTC": status.get("StartUTC"),
               "EndUTC": status.get("EndUTC"), "JobTotalSeconds": status.get("TotalSeconds"), "Nodes": nodes,
               "JobNodeHours": (status.get("TotalSeconds") or 0.0) * nodes / 3600.0, "Stages": {}}
    if status.get("Jobs"):
        summary["Jobs"] = {name: {**job, "NodeHours": (job.get("TotalSeconds") or 0.0) * nodes / 3600.0}
                           for name, job in status["Jobs"].items()}
        summary["JobTotalSecondsRule"] = "the sum of the split jobs' runner totals (node seconds); per job under Jobs"
    for prefix in sorted({name[:-len("-worker")] for name in stages if name.endswith("-worker")}):
        if all(stages.get(f"{prefix}-{kind}", {}).get("State") == "complete" for kind in ("worker", "reducer")):
            summary["Stages"][prefix] = stage_cost(stages, prefix, full_sources=full_sources, nodes=nodes, block_size=block_size)
        else:
            summary["Stages"][prefix] = {"State": {kind: stages.get(f"{prefix}-{kind}", {}).get("State") for kind in ("worker", "reducer")}}
    for name, stage in stages.items():
        if not name.endswith(("-worker", "-reducer")):
            summary["Stages"][name] = {"State": stage.get("State"), "WallSeconds": stage.get("WallSeconds"),
                                       "NodeHours": (stage.get("WallSeconds") or 0.0) * nodes / 3600.0,
                                       "H1": stage.get("Parsed", {}).get("H1"), "PCG": stage.get("Parsed", {}).get("PCG"),
                                       "PalacePeakTotal": stage.get("Parsed", {}).get("PalacePeakMemory", {}).get("Total")}
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("status", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--full-sources", type=int, required=True, help="the coupon's source count (extrapolation target)")
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--block-size", type=int, default=DEFAULT_BLOCK_SIZE)
    args = parser.parse_args(argv)
    summary = summarize(json.loads(args.status.read_text()), full_sources=args.full_sources, nodes=args.nodes,
                        block_size=args.block_size)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "Stages"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
