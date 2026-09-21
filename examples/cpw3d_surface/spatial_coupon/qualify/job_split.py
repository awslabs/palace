# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The per-coupon SOURCE SPLIT of `coupon-library qualify` (supervisor decision 61b):
a coupon's main-order source set is partitioned into N contiguous blocks run as N
independent worker jobs - the same mesh and configs apart from the PrescribedPotential
subset - every block writing into the stage's ONE archive directory (Palace archives one
file per source, rank and field: the union is the directory), and the main-order
reducer runs once, in its own job, on that union after every worker job completed
(the reduction is linear in the archive union: the matrices equal a single job's to
roundoff).  The p-sequence controls and the local-edge stage stay in the first worker
job; when they leave no room for a source block there, the first job carries them alone
(recorded as ControlsJob "separate").  N = 1 is today's single job, byte-identical.

Policy {Mode: speed | frugal | fixed, MaxJobs, WalltimeSeconds, FixedJobs}, every
job estimated from the recorded cost model at the largest PCG factor with the
preflight-and-margin factor (estimate_stages.estimate's fail-closed rule):
  frugal - the smallest N whose every job fits the walltime;
  speed  - the N minimizing the estimated critical path (the longest worker job, then
           the reducer job); the smallest such N on a tie, so it never uses more jobs
           than shorten the path;
  fixed  - N = FixedJobs (fail closed when a job does not fit).
N is bounded by MaxJobs (the run's concurrency, --max-jobs), the cluster's user job cap
and the source count.  The blocks are balanced on the estimate: the first job's block
is shortened by the controls' and local-edge share so every worker job ends together.
A coupon fails closed only when even the maximal split does not fit (or the reducer
job / the first job's fixed stages alone do not fit any job).
"""
import math

MODES = ("speed", "frugal", "fixed")
DEFAULT_MODE = "frugal"
SPLIT_RULE = ("N contiguous source blocks as N independent worker jobs at every main order (job 1 also runs the "
              "p-sequence controls and the local-edge stage), every block archived into the stage's one archive "
              "directory; one reducer job per coupon on the archive union after every worker job completed; N = 1 "
              "is the single job of the recorded campaigns; estimates per job = the cost model at the largest PCG "
              "factor x PreflightAndMarginFactor + PreflightSeconds; a job fits when that is below the walltime")


def normalize_policy(mode, *, max_jobs, walltime_seconds, fixed_jobs=None, user_job_cap=None, origin=None):
    if mode not in MODES:
        raise ValueError(f"job policy mode {mode!r} is not one of {MODES}")
    if mode == "fixed" and (not isinstance(fixed_jobs, int) or isinstance(fixed_jobs, bool) or fixed_jobs < 1):
        raise ValueError("job policy 'fixed' needs FixedJobs >= 1 (--fixed-jobs)")
    if not isinstance(max_jobs, int) or max_jobs < 1:
        raise ValueError("job policy needs MaxJobs >= 1")
    return {"Mode": mode, "MaxJobs": int(max_jobs), "WalltimeSeconds": float(walltime_seconds),
            "FixedJobs": int(fixed_jobs) if mode == "fixed" else None,
            "UserJobCap": user_job_cap, "Origin": origin or "command line", "Rule": SPLIT_RULE}


def balanced_blocks(source_count, jobs, per_source, extras):
    """Contiguous block sizes (job 1 first) balancing s * per_source (+ extras on job 1)."""
    if jobs <= 1:
        return [source_count]
    if per_source <= 0.0:
        base, rest = divmod(source_count, jobs)
        return [base + (1 if k < rest else 0) for k in range(jobs)]
    target = (source_count * per_source + extras) / jobs
    first = int(math.floor((target - extras) / per_source))
    first = max(0, min(source_count, first))
    candidates = []
    for s1 in {first, min(source_count, first + 1)}:
        rest = source_count - s1
        base, remainder = divmod(rest, jobs - 1)
        others = [base + (1 if k < remainder else 0) for k in range(jobs - 1)]
        longest = max(s1 * per_source + extras, max(others) * per_source)
        candidates.append((longest, -s1, [s1] + others))
    return min(candidates)[2]


def stage_seconds(estimate_stage, factor):
    return estimate_stage["ByPCGFactor"][factor]["StageSecondsEstimate"]


def plan_split(*, indices, layout, estimate, policy, model, profile):
    """The split record for a coupon: `indices` = the main-order source indices in
    config order, `layout` = qualify_library.stage_layout, `estimate` =
    estimate_stages.estimate over that layout (its Stages carry the per-source and
    non-source worker seconds of every main stage)."""
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    worst = max(factors, key=float)
    margin, preflight = model["PreflightAndMarginFactor"], model["PreflightSeconds"]
    walltime = policy["WalltimeSeconds"]
    mains = [item for item in layout if item["Role"] == "main"]
    fixed_stages = [item for item in layout if item["Role"] != "main"]
    source_count = len(indices)

    def job_seconds(block_size, *, with_fixed, factor):
        seconds = 0.0
        for item in mains:
            stage = estimate["Stages"][item["EstimateKey"]]
            if block_size > 0:
                seconds += stage["WorkerNonSourceSecondsEstimate"] + block_size * stage["ByPCGFactor"][factor]["PerSourceSecondsEstimate"]
        if with_fixed:
            seconds += sum(stage_seconds(estimate["Stages"][item["EstimateKey"]], factor) for item in fixed_stages)
        return seconds

    def reducer_seconds():
        return sum(estimate["Stages"][item["EstimateKey"]]["ReducerSecondsEstimate"] for item in mains)

    def with_margin(seconds):
        return seconds * margin + preflight

    per_source_worst = sum(estimate["Stages"][item["EstimateKey"]]["ByPCGFactor"][worst]["PerSourceSecondsEstimate"]
                           for item in mains)
    extras_worst = job_seconds(0, with_fixed=True, factor=worst)

    def candidate(jobs):
        if jobs == 1:
            single = {factor: with_margin(job_seconds(source_count, with_fixed=True, factor=factor) + reducer_seconds())
                      for factor in factors}
            return {"N": 1, "Blocks": [list(indices)], "ControlsJob": "worker-1",
                    "Jobs": [{"Name": "single", "Kind": "single", "Block": 1, "Sources": list(indices),
                              "SecondsEstimateWithPreflightAndMargin": single, "Fits": single[worst] < walltime}],
                    "CriticalPathEstimateSeconds": single, "NodeSecondsEstimate": single,
                    "Fits": single[worst] < walltime}
        sizes = balanced_blocks(source_count, jobs, per_source_worst, extras_worst)
        blocks, start = [], 0
        for size in sizes:
            blocks.append(list(indices[start:start + size]))
            start += size
        job_records = []
        for k, block in enumerate(blocks, start=1):
            seconds = {factor: with_margin(job_seconds(len(block), with_fixed=(k == 1), factor=factor)) for factor in factors}
            job_records.append({"Name": f"worker-{k}", "Kind": "worker", "Block": k, "Sources": block,
                                "SecondsEstimateWithPreflightAndMargin": seconds, "Fits": seconds[worst] < walltime})
        reducer = {factor: with_margin(reducer_seconds()) for factor in factors}
        job_records.append({"Name": "reducer", "Kind": "reducer", "Block": None, "Sources": list(indices),
                            "SecondsEstimateWithPreflightAndMargin": reducer, "Fits": reducer[worst] < walltime})
        critical = {factor: max(job["SecondsEstimateWithPreflightAndMargin"][factor] for job in job_records[:-1]) + reducer[factor]
                    for factor in factors}
        node_seconds = {factor: sum(job["SecondsEstimateWithPreflightAndMargin"][factor] for job in job_records) for factor in factors}
        return {"N": jobs, "Blocks": blocks, "ControlsJob": "worker-1" if blocks[0] else "separate",
                "Jobs": job_records, "CriticalPathEstimateSeconds": critical, "NodeSecondsEstimate": node_seconds,
                "Fits": all(job["Fits"] for job in job_records)}

    largest = max(1, min(policy["MaxJobs"], policy.get("UserJobCap") or policy["MaxJobs"], source_count))
    if policy["Mode"] == "fixed":
        if policy["FixedJobs"] > max(largest, 1):
            raise ValueError(f"job policy fixed {policy['FixedJobs']} exceeds the largest split {largest} "
                             f"(MaxJobs {policy['MaxJobs']}, user job cap {policy.get('UserJobCap')}, {source_count} sources)")
        wanted = [policy["FixedJobs"]]
    else:
        wanted = list(range(1, largest + 1))
    candidates = [candidate(jobs) for jobs in wanted]
    summary = [{"N": item["N"], "Fits": item["Fits"], "ControlsJob": item["ControlsJob"],
                "LongestJobSeconds": max(job["SecondsEstimateWithPreflightAndMargin"][worst] for job in item["Jobs"]),
                "CriticalPathSeconds": item["CriticalPathEstimateSeconds"][worst],
                "NodeSeconds": item["NodeSecondsEstimate"][worst]} for item in candidates]
    fitting = [item for item in candidates if item["Fits"]]
    chosen = None
    if fitting:
        if policy["Mode"] == "frugal":
            chosen = min(fitting, key=lambda item: item["N"])
        elif policy["Mode"] == "speed":
            chosen = min(fitting, key=lambda item: (item["CriticalPathEstimateSeconds"][worst], item["N"]))
        else:
            chosen = fitting[0]
    record = {"Policy": policy, "Rule": SPLIT_RULE, "WorstPCGFactor": worst, "WalltimeSeconds": walltime,
              "LargestSplit": largest, "Candidates": summary, "SourceCount": source_count}
    if chosen is None:
        longest = min(summary, key=lambda item: item["LongestJobSeconds"])
        record.update({"Fits": False, "N": None, "Blocks": None, "Jobs": None, "ControlsJob": None,
                       "CriticalPathEstimateSeconds": None, "NodeSecondsEstimate": None,
                       "Decision": (f"does NOT fit: even the maximal split N = {max(wanted)} leaves a job of "
                                    f"{longest['LongestJobSeconds'] / 60:.0f} min at {worst}x PCG (+{100 * (margin - 1):.0f}% and preflight) "
                                    f"against the {walltime / 3600:.1f} h walltime (policy {policy['Mode']}, MaxJobs {policy['MaxJobs']}): "
                                    f"fail closed, not submitted")})
        return record
    record.update({key: chosen[key] for key in ("N", "Blocks", "Jobs", "ControlsJob", "CriticalPathEstimateSeconds",
                                                  "NodeSecondsEstimate", "Fits")})
    record["Decision"] = (f"policy {policy['Mode']}: N = {chosen['N']} of at most {largest} "
                          f"(blocks {[len(block) for block in chosen['Blocks']]} sources; controls + local-edge in "
                          f"{chosen['ControlsJob']}); longest job {max(job['SecondsEstimateWithPreflightAndMargin'][worst] for job in chosen['Jobs']) / 60:.0f} min, "
                          f"critical path {chosen['CriticalPathEstimateSeconds'][worst] / 60:.0f} min, node time "
                          f"{chosen['NodeSecondsEstimate'][worst] / 3600:.2f} h at {worst}x PCG (+{100 * (margin - 1):.0f}% and preflight) "
                          f"against the {walltime / 3600:.1f} h walltime")
    return record
