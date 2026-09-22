#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Refit qualify/cost-model.json from the measured rates of a completed
`coupon-library qualify` run (decision 64a: the device library re-run with the frozen
executable 170439c4... and the streaming reducer at block size 48).

Every rate of the model is expressed on ONE reference mesh (--reference-case: its H1
entity counts from the build record give the closed-form H1 of every stage, the check
estimate_stages.load_cost_model applies) and is the LARGEST value the run's coupons
imply once scaled to that mesh by the exact H1 ratio - the fail-closed side: the
estimate at the measured PCG counts is at least the measured wall of every stage of
every coupon of the run (the self-check the record carries, per coupon and stage).
Per worker + reducer stage (p3 / p5 controls at 8 sources, the p4 main at every
source): seconds per PCG iteration, mean / max PCG counts, the per-source seconds
outside the solve, the worker's non-source wall (wall - sources x mean per-source
seconds: launch overhead included), the reducer's setup wall and its block-pair
seconds inverted through the estimator's evaluation / Gram split at the reference
source count, Palace peaks, node used GiB (the non-Palace remainder taken as the
largest measured) and the archive GB per source; the local-edge stage's solve and
non-solve wall from Palace's elapsed-time report.  ReducerEvaluationFraction and
ReducerResidentFieldGBPerMillionH1 cannot be re-measured from a run at one block
size and keep the previous values (the resident-field growth for b > 48 is an upper
bound under the streaming Gram).  The previous model is kept as a file and bound by
digest under Previous.

usage: refit_cost_model.py --qualification library-qualification.json [--build-record library-build.json]
       --previous cost-model.json --previous-kept PATH [--reference-case CASE] --out cost-model.json --record refit.json
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import estimate_stages  # noqa: E402
from mixed_mesh import h1_dofs_from_counts  # noqa: E402

GIB = 2**30
MEMORY_UNITS = {"K": 2**10, "M": 2**20, "G": 2**30, "T": 2**40}
LINEAR_SOLVE_TIMERS = ("Linear Solve", "Setup", "Preconditioner", "Coarse Solve")
REDUCTION_TIMERS = ("Archive Reduction", "Archive Read", "Sample Evaluation", "Gram Assembly", "Domain Gram")
KEPT_KEYS = ("ReducerEvaluationFraction", "ReducerEvaluationFractionRule", "ReducerResidentFieldGBPerMillionH1",
             "ReducerResidentFieldGBRule", "PCGFactors", "PreflightAndMarginFactor", "PreflightSeconds", "NodeFitFraction",
             "PalaceGBPerGiB")


class RefitError(ValueError):
    """The run's records do not support a refit."""


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recorded_path(path):
    """A record path relative to the spatial_coupon directory when it lies inside it (the
    committed records), else absolute: the refit of committed records is reproducible."""
    path = Path(path).resolve()
    try:
        return str(path.relative_to(HERE.parent))
    except ValueError:
        return str(path)


def round_up(value, decimals):
    """Rounded towards +infinity at `decimals` (a recorded rate never rounds below the measurement)."""
    scale = 10 ** decimals
    return math.ceil(value * scale - 1e-9) / scale


def memory_gb(text):
    """Palace's peak-memory figure '48.2G' -> 48.2, the model's PalacePeakGB unit (the
    figure as Palace prints it, the convention of the physics-11 model; PalaceGBPerGiB
    converts it to the node's GiB); K / M / T suffixes are scaled to G."""
    text = str(text).strip()
    unit = text[-1]
    if unit in MEMORY_UNITS:
        return float(text[:-1]) * MEMORY_UNITS[unit] / MEMORY_UNITS["G"]
    return float(text)


def parse_elapsed_time_report(text):
    """Palace's elapsed-time report -> {timer name: average seconds} (the third column;
    the indented sub-timers are additive with their parent line, as the Total shows)."""
    timers = {}
    for line in text.splitlines():
        match = re.match(r"^\s*([A-Za-z][A-Za-z ]*?)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s*$", line)
        if match and match.group(1).strip() not in ("Total",):
            timers[match.group(1).strip()] = float(match.group(4))
    if not timers:
        raise RefitError("no timers in the elapsed-time report")
    return timers


def linear_solve_seconds(timers):
    return sum(timers.get(name, 0.0) for name in LINEAR_SOLVE_TIMERS)


def archive_gb(deletion):
    """The archive GB of a case's main stage from its ArchiveDeletion record (du -sh
    sizes before deletion, one line per archive: the main stage's is the largest; None
    when the run recorded no sizes)."""
    sizes = (deletion or {}).get("SizesBeforeDeletion") or ""
    values = []
    for line in sizes.splitlines():
        size = line.split("\t")[0].strip()
        if size:
            values.append(memory_gb(size))
    return max(values) if values else None


def resolve(record_dir, case, name, recorded):
    """A record's path: the copied one next to the qualification record, else the recorded one."""
    copied = record_dir / case / name
    if copied.exists():
        return copied
    if recorded and Path(recorded).exists():
        return Path(recorded)
    raise RefitError(f"{case}: {name} is neither at {copied} nor at the recorded {recorded}")


def load_statuses(record_dir, case_record):
    """The run_stages status of every job of the coupon (job name -> status): copied under
    <case>/results/<job>/status.json, else the run's results tree."""
    case = case_record["Case"]
    root = Path(case_record["Root"])
    statuses = {}
    for job in case_record["Jobs"]:
        name = job["Name"]
        relative = Path(job["RemoteDirectory"]).relative_to(case_record["Remote"]["Case"] + "/main")
        candidates = [record_dir / case / "results" / name / "status.json", root / "results" / "main" / relative / "status.json"]
        path = next((path for path in candidates if path.exists()), None)
        if path is None:
            raise RefitError(f"{case}: no status.json of job {name} at {candidates}")
        statuses[name] = json.loads(path.read_text())
    return statuses


def stage_index(statuses):
    """Stage name -> (job name, stage record) over the coupon's jobs."""
    index = {}
    for job, status in statuses.items():
        for stage in status["Stages"]:
            index[stage["Name"]] = (job, stage)
    return index


def stage_measurement(prefix, cost_stage, stages):
    """The measured quantities of one worker + reducer stage: the merged per-source rates
    of summarize_cost, the worker's non-source wall per JOB (the largest block's wall minus
    its sources' total seconds: one launch per worker job) and the reducer's wall split by
    Palace's elapsed-time report into the setup and the per-source reduction work."""
    n = len(cost_stage["Sources"])
    mean_total = cost_stage["MeanPerSourceTotalSeconds"]
    blocks = []
    for name, (job, stage) in stages.items():
        if name == f"{prefix}-worker" or name.startswith(f"{prefix}-worker-block"):
            timings = stage["Parsed"].get("SourceTiming") or []
            total = sum(t["TotalSeconds"] for t in timings)
            blocks.append({"Name": name, "Job": job, "Sources": len(timings), "WallSeconds": stage["WallSeconds"],
                           "NonSourceSeconds": max(stage["WallSeconds"] - total, 0.0),
                           "MeanPerSourceSeconds": total / len(timings) if timings else None})
    if not blocks:
        raise RefitError(f"{prefix}: no worker stage in the statuses")
    # Contiguous source blocks differ in their PCG counts: the per-source seconds the
    # planner budgets every block with are the largest block mean (= the coupon mean for
    # a single worker job), so every worker job of the run is covered.
    block_mean = max(block["MeanPerSourceSeconds"] for block in blocks if block["MeanPerSourceSeconds"] is not None)
    reducer_job, reducer = stages[f"{prefix}-reducer"]
    timers = parse_elapsed_time_report(reducer["Parsed"]["ElapsedTimeReport"])
    reduction_fraction = sum(timers.get(name, 0.0) for name in REDUCTION_TIMERS) / sum(timers.values())
    reducer_wall = reducer["WallSeconds"]
    pair = reducer_wall * reduction_fraction
    return {"H1": cost_stage["H1"], "Sources": n, "Order": cost_stage["Order"],
            "SecondsPerPCGIteration": cost_stage["SecondsPerPCGIteration"],
            "MeanPCGIterations": cost_stage["MeanPCGIterations"], "MaxPCGIterations": cost_stage["MaxPCGIterations"],
            "MeanPerSourceSeconds": mean_total, "LargestBlockMeanPerSourceSeconds": block_mean,
            "PerSourceNonSolveSeconds": block_mean - cost_stage["SecondsPerPCGIteration"] * cost_stage["MeanPCGIterations"],
            "WorkerWallSeconds": cost_stage["WorkerWallSeconds"], "WorkerBlocks": blocks,
            "WorkerNonSourceSeconds": max(block["NonSourceSeconds"] for block in blocks),
            "ReducerWallSeconds": reducer_wall, "ReducerReductionFraction": reduction_fraction,
            "ReducerPairSeconds": pair, "ReducerSetupSeconds": reducer_wall - pair, "ReducerJob": reducer_job,
            "ReducerBlockSize": cost_stage["ReducerBlockSize"], "StageWallSeconds": cost_stage["StageWallSeconds"],
            "WorkerPalacePeakGB": memory_gb(cost_stage["Memory"]["WorkerPalacePeakTotal"]),
            "ReducerPalacePeakGB": memory_gb(cost_stage["Memory"]["ReducerPalacePeakTotal"]),
            "WorkerNodeUsedGiB": cost_stage["Memory"]["WorkerNodePeakUsedGiBSampled"],
            "ReducerNodeUsedGiB": cost_stage["Memory"]["ReducerNodePeakUsedGiBSampled"]}


def local_edge_measurement(name, stages):
    job, stage = stages[name]
    parsed = stage["Parsed"]
    timers = parse_elapsed_time_report(parsed["ElapsedTimeReport"])
    solve = linear_solve_seconds(timers)
    wall = stage["WallSeconds"]
    return {"Name": name, "Job": job, "Order": parsed["Order"], "H1": parsed["H1"], "Sources": len(parsed["Iterations"]),
            "WallSeconds": wall, "PalaceTotalSeconds": parsed["PalaceTotalSeconds"],
            "LinearSolveSeconds": solve, "NonSolveSeconds": max(wall - solve, 0.0),
            "PalacePeakGB": memory_gb(parsed["PalacePeakMemory"]["Total"]),
            "NodeUsedGiB": (stage.get("NodePeakUsedBytesSampled") or 0) / GIB, "PCG": parsed.get("PCG")}


def measure_case(record_dir, case_record, build_case):
    case = case_record["Case"]
    cost_summary = json.loads(resolve(record_dir, case, "cost-summary.json", case_record["Cost"]["Path"]).read_text())
    statuses = load_statuses(record_dir, case_record)
    stages = stage_index(statuses)
    counts = build_case["H1"]["EntityCounts"]
    measured_stages = {}
    local = None
    for item in case_record["Stages"]:
        if item["Kind"] == "response":
            cost_stage = cost_summary["Stages"].get(item["Prefix"])
            if not cost_stage or "H1" not in cost_stage:
                raise RefitError(f"{case}: stage {item['Prefix']} has no complete cost record")
            measured = stage_measurement(item["Prefix"], cost_stage, stages)
            measured["Role"] = item["Role"]
            measured_stages[f"p{item['Order']}"] = measured
        elif item["Kind"] == "local-edge":
            measured = local = local_edge_measurement(item["Prefix"], stages)
        else:
            continue
        closed = h1_dofs_from_counts(counts, item["Order"])
        if closed != measured["H1"]:
            raise RefitError(f"{case}: closed-form H1 {closed} at p{item['Order']} differs from the run's {measured['H1']}")
    if local is None:
        raise RefitError(f"{case}: no local-edge stage")
    jobs = {}
    for name, status in statuses.items():
        job = case_record["Cost"]["Jobs"][name]
        jobs[name] = {"PBSJobID": job["PBSJobID"], "ActualSeconds": job["ActualSeconds"],
                      "EstimateSecondsWithPreflightAndMargin": job["EstimateSecondsWithPreflightAndMargin"],
                      "Stages": [{"Name": stage["Name"], "WallSeconds": stage["WallSeconds"],
                                  "Sources": len(stage["Parsed"].get("SourceTiming") or [])} for stage in status["Stages"]]}
    return {"Case": case, "EntityCounts": counts, "MeshSHA256": case_record["Mesh"]["SHA256"],
            "Elements": build_case["Elements"], "Stages": measured_stages, "LocalEdge": local,
            "ArchiveGB": archive_gb(case_record.get("ArchiveDeletion")), "Jobs": jobs}


def _largest(values):
    """(value, case) of the largest scaled measurement."""
    return max(values, key=lambda pair: pair[0])


def refit_stage(order, measurements, reference, previous):
    """The model's stage record at `order` on the reference mesh: every rate the largest
    the coupons imply once scaled to the reference H1 (ratio = H1_ref / H1_c)."""
    key = f"p{order}"
    ref_h1 = h1_dofs_from_counts(reference["EntityCounts"], order)
    ref_sources = reference["Stages"][key]["Sources"]
    block_size = {m["Stages"][key]["ReducerBlockSize"] for m in measurements}
    if len(block_size) != 1:
        raise RefitError(f"{key}: the coupons ran the reducer at different block sizes {sorted(block_size)}")
    block_size = block_size.pop()
    fraction = previous["ReducerEvaluationFraction"]
    ref_evaluations = estimate_stages.source_evaluations(ref_sources, block_size)
    ref_pairs = ref_sources * (ref_sources + 1) // 2
    gb_per_gib = previous["PalaceGBPerGiB"]
    scaled = {name: [] for name in ("SecondsPerPCGIteration", "PerSourceNonSolveSeconds", "WorkerNonSourceSeconds",
                                    "ReducerSetupSeconds", "ReducerPairSeconds", "WorkerPalacePeakGB", "ReducerPalacePeakGB",
                                    "WorkerNodeOverheadGiB", "ReducerNodeOverheadGiB", "ArchiveGBPerSource")}
    pcg_mean, pcg_max = [], []
    for m in measurements:
        s = m["Stages"][key]
        ratio = ref_h1 / s["H1"]
        for name in ("SecondsPerPCGIteration", "PerSourceNonSolveSeconds", "WorkerNonSourceSeconds", "ReducerSetupSeconds",
                     "WorkerPalacePeakGB", "ReducerPalacePeakGB"):
            scaled[name].append((s[name] * ratio, m["Case"]))
        evaluations = estimate_stages.source_evaluations(s["Sources"], block_size)
        pairs = s["Sources"] * (s["Sources"] + 1) // 2
        scale = fraction * evaluations / ref_evaluations + (1.0 - fraction) * pairs / ref_pairs
        scaled["ReducerPairSeconds"].append((s["ReducerPairSeconds"] * ratio / scale, m["Case"]))
        scaled["WorkerNodeOverheadGiB"].append((s["WorkerNodeUsedGiB"] - s["WorkerPalacePeakGB"] / gb_per_gib, m["Case"]))
        scaled["ReducerNodeOverheadGiB"].append((s["ReducerNodeUsedGiB"] - s["ReducerPalacePeakGB"] / gb_per_gib, m["Case"]))
        if s["Role"] == "main" and m["ArchiveGB"] is not None:
            scaled["ArchiveGBPerSource"].append((m["ArchiveGB"] * ratio / s["Sources"], m["Case"]))
        pcg_mean.append((s["MeanPCGIterations"], m["Case"]))
        pcg_max.append((s["MaxPCGIterations"], m["Case"]))
    chosen = {name: _largest(values) for name, values in scaled.items() if values}
    chosen["MeanPCGIterations"] = _largest(pcg_mean)
    chosen["MaxPCGIterations"] = _largest(pcg_max)
    spi, mean_pcg = chosen["SecondsPerPCGIteration"][0], chosen["MeanPCGIterations"][0]
    worker_peak, reducer_peak = chosen["WorkerPalacePeakGB"][0], chosen["ReducerPalacePeakGB"][0]
    reducer_pair = chosen["ReducerPairSeconds"][0]
    if "ArchiveGBPerSource" in chosen:
        archive = chosen["ArchiveGBPerSource"][0] * ref_sources
    else:
        # A control stage: the previous model's archive rate scaled to the reference (the
        # control archives are deleted with the main's, their sizes not recorded apart).
        prev = previous["Stages"][key]
        archive = prev["ArchiveGB"] * ref_h1 / prev["H1"] * ref_sources / prev["Sources"]
        chosen["ArchiveGBPerSource"] = (archive / ref_sources, "previous model (scaled)")
    stage = {"H1": ref_h1, "Sources": ref_sources, "SecondsPerPCGIteration": round_up(spi, 5),
             "MeanPCGIterations": round_up(mean_pcg, 3), "MaxPCGIterations": int(chosen["MaxPCGIterations"][0]),
             "MeanPerSourceSeconds": round_up(chosen["PerSourceNonSolveSeconds"][0] + spi * mean_pcg, 3),
             "WorkerNonSourceSeconds": round_up(chosen["WorkerNonSourceSeconds"][0], 2),
             "ReducerPalaceSeconds": round_up(chosen["ReducerSetupSeconds"][0] + reducer_pair, 2),
             "ReducerPairSeconds": round_up(reducer_pair, 2),
             "WorkerPalacePeakGB": round_up(worker_peak, 1), "ReducerPalacePeakGB": round_up(reducer_peak, 1),
             "WorkerNodeUsedGiB": round_up(chosen["WorkerNodeOverheadGiB"][0] + worker_peak / gb_per_gib, 2),
             "ReducerNodeUsedGiB": round_up(chosen["ReducerNodeOverheadGiB"][0] + reducer_peak / gb_per_gib, 2),
             "ArchiveGB": round_up(archive, 2)}
    provenance = {name: {"Value": value, "SetBy": case} for name, (value, case) in chosen.items()}
    return stage, provenance, block_size


def refit_local_edge(measurements, reference):
    order = reference["LocalEdge"]["Order"]
    ref_h1 = h1_dofs_from_counts(reference["EntityCounts"], order)
    sources = {m["LocalEdge"]["Sources"] for m in measurements}
    if len(sources) != 1:
        raise RefitError(f"the local-edge stages ran on different control counts {sorted(sources)}")
    scaled = {name: [] for name in ("NonSolveSeconds", "LinearSolveSeconds", "PalacePeakGB", "NodeUsedGiB")}
    for m in measurements:
        local = m["LocalEdge"]
        if local["Order"] != order:
            raise RefitError(f"{m['Case']}: local-edge order {local['Order']} differs from the reference's {order}")
        ratio = ref_h1 / local["H1"]
        for name in scaled:
            scaled[name].append((local[name] * ratio, m["Case"]))
    chosen = {name: _largest(values) for name, values in scaled.items()}
    solve, non_solve = chosen["LinearSolveSeconds"][0], chosen["NonSolveSeconds"][0]
    stage = {"Order": order, "H1": ref_h1, "Sources": sources.pop(), "PalacePeakGB": round_up(chosen["PalacePeakGB"][0], 1),
             "NodeUsedGiB": round_up(chosen["NodeUsedGiB"][0], 1), "PalaceTotalSeconds": round_up(solve + non_solve, 1),
             "LinearSolveSeconds": round_up(solve, 1), "NonSolveSeconds": round_up(non_solve, 1)}
    return stage, {name: {"Value": value, "SetBy": case} for name, (value, case) in chosen.items()}


def self_check(model, measurements, block_size, profile):
    """estimate_stages on every coupon of the run under `model`, read per JOB as the split
    planner does: a worker job costs the non-source estimate plus its sources x the
    per-source estimate for every block it ran, plus the control and local-edge stages it
    carried; the reducer job the reducer estimate.  Per job: the 1.0x estimate over the
    sum of its measured stage walls (>= 1 required), the worst-PCG-factor estimate with
    preflight and margin over the job's runner total (the recorded ActualOverEstimate
    inverted).  Per coupon: the single-job estimate and whether it fits."""
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    worst = max(factors, key=float)
    margin, preflight = model["PreflightAndMarginFactor"], model["PreflightSeconds"]
    out = {}
    for m in measurements:
        stages = [(f"{order}-{s['Sources']}", int(order[1:]), s["Sources"]) for order, s in m["Stages"].items()]
        local = m["LocalEdge"]
        estimate = estimate_stages.estimate(m["EntityCounts"], stages, model=model, profile=profile, block_size=block_size,
                                            local_edge=(local["Name"], local["Order"], local["Sources"]))
        by_prefix = {}
        for order, s in m["Stages"].items():
            by_prefix[order] = estimate["Stages"][f"{order}-{s['Sources']}"]
        prefixes = {order: (f"{m['Case']}-p{order[1:]}" if s["Role"] == "main" else f"{m['Case']}-p{order[1:]}-control")
                    for order, s in m["Stages"].items()}

        def stage_estimate(stage, factor):
            name, sources = stage["Name"], stage["Sources"]
            if name == local["Name"]:
                return estimate["Stages"][name]["ByPCGFactor"][factor]["StageSecondsEstimate"]
            for order, prefix in prefixes.items():
                est = by_prefix[order]
                if name == f"{prefix}-reducer":
                    return est["ReducerSecondsEstimate"]
                if name == f"{prefix}-worker" or name.startswith(f"{prefix}-worker-block"):
                    return est["WorkerNonSourceSecondsEstimate"] + sources * est["ByPCGFactor"][factor]["PerSourceSecondsEstimate"]
            raise RefitError(f"{m['Case']}: stage {name} is not in the estimate")

        jobs = {}
        for name, job in m["Jobs"].items():
            actual_stages = sum(stage["WallSeconds"] for stage in job["Stages"])
            est_1x = sum(stage_estimate(stage, "1.0") for stage in job["Stages"])
            est_worst = sum(stage_estimate(stage, worst) for stage in job["Stages"]) * margin + preflight
            jobs[name] = {"Stages": [stage["Name"] for stage in job["Stages"]], "ActualStageWallSeconds": actual_stages,
                          "ActualJobSeconds": job["ActualSeconds"], "Estimate1xSeconds": est_1x,
                          "Estimate1xOverStageWall": est_1x / actual_stages,
                          "EstimateWorstWithMarginSeconds": est_worst,
                          "EstimateWorstWithMarginOverJob": est_worst / job["ActualSeconds"] if job["ActualSeconds"] else None,
                          "PerStage": {stage["Name"]: stage_estimate(stage, "1.0") / stage["WallSeconds"] for stage in job["Stages"]}}
        out[m["Case"]] = {"Jobs": jobs,
                          "MinimumJobEstimate1xOverStageWall": min(job["Estimate1xOverStageWall"] for job in jobs.values()),
                          "MinimumStageEstimate1xOverWall": min(ratio for job in jobs.values() for ratio in job["PerStage"].values()),
                          "SingleJob": {"FitsOneJob": estimate["FitsOneJob"],
                                        "JobSecondsEstimateByPCGFactor": estimate["JobSecondsEstimateByPCGFactor"],
                                        "JobSecondsEstimateWithPreflightAndMargin": estimate["JobSecondsEstimateWithPreflightAndMargin"],
                                        "MaxPalacePeakGBEstimate": estimate["MaxPalacePeakGBEstimate"]},
                          "ActualNodeSeconds": sum(job["ActualSeconds"] for job in m["Jobs"].values()),
                          "EstimateWorstWithMarginOverActualNodeSeconds": (
                              sum(job["EstimateWorstWithMarginSeconds"] for job in jobs.values())
                              / sum(job["ActualSeconds"] for job in m["Jobs"].values()))}
    out["MinimumStageEstimateOverActual"] = min(case["MinimumStageEstimate1xOverWall"] for case in out.values())
    out["Conservative"] = out["MinimumStageEstimateOverActual"] >= 1.0 - 1e-9
    return out


def refit(qualification_path, *, build_record_path=None, previous_path, previous_kept, reference_case=None, profile=None):
    qualification_path = Path(qualification_path)
    record_dir = qualification_path.parent
    qualification = json.loads(qualification_path.read_text())
    build_path = Path(build_record_path) if build_record_path else record_dir / "library-build.json"
    if not build_path.exists():
        build_path = Path(qualification["BuildRecord"]["Path"])
    build = json.loads(build_path.read_text())
    build_cases = {case["Case"]: case for case in build["Cases"]}
    previous = json.loads(Path(previous_path).read_text())
    previous_kept = Path(previous_kept)
    if sha256(previous_kept) != sha256(previous_path):
        raise RefitError(f"the kept previous model {previous_kept} is not byte-identical to {previous_path}")
    profile = profile or json.loads(estimate_stages.CLUSTER_PROFILE.read_text())
    measurements = [measure_case(record_dir, case, build_cases[case["Case"]]) for case in qualification["Cases"]
                    if case.get("Cost") and case["Cost"].get("Jobs")]
    if not measurements:
        raise RefitError("no coupon of the run has a complete cost record")
    if reference_case is None:
        reference = max(measurements, key=lambda m: max(s["H1"] for s in m["Stages"].values()))
    else:
        reference = next((m for m in measurements if m["Case"] == reference_case), None)
        if reference is None:
            raise RefitError(f"--reference-case {reference_case} is not a measured coupon of the run")
    orders = sorted({order for m in measurements for order in m["Stages"]})
    stages, provenance, block_sizes = {}, {}, set()
    for order in orders:
        with_stage = [m for m in measurements if order in m["Stages"]]
        stages[order], provenance[order], block_size = refit_stage(int(order[1:]), with_stage, reference, previous)
        block_sizes.add(block_size)
    if len(block_sizes) != 1:
        raise RefitError(f"the stages ran the reducer at different block sizes {sorted(block_sizes)}")
    block_size = block_sizes.pop()
    local_edge, provenance["LocalEdge"] = refit_local_edge(measurements, reference)
    executables = sorted({(case.get("FrozenExecutable") or {}).get("SHA256") for case in qualification["Cases"]
                          if case.get("Cost") and case["Cost"].get("Jobs")})
    if len(executables) != 1 or executables[0] is None:
        raise RefitError(f"the coupons ran different frozen executables {executables}")
    executable = executables[0]
    pbs = sorted(job["PBSJobID"].split(".")[0] for m in measurements for job in m["Jobs"].values() if job.get("PBSJobID"))
    model = {"Version": 2,
             "Copyright": previous["Copyright"], "SPDX-License-Identifier": previous["SPDX-License-Identifier"],
             "Purpose": ("Measured per-stage cost rates of the Palace response worker / reducer / ordinary local-edge stages on "
                         f"one r8g.48xlarge node (192 ranks, frozen executable {executable}, streaming one-pass Gram "
                         f"reducer at block size {block_size}) refit by refit_cost_model.py from the "
                         f"{len(measurements)}-coupon device library qualification {qualification_path.name} (tool commit "
                         f"{qualification.get('ToolCommit')}, PBS {', '.join(pbs)}). Every rate is on the reference mesh "
                         f"{reference['Case']} and is the LARGEST the run's coupons imply once scaled to it by the exact H1 "
                         "ratio (fail-closed: the 1x estimate is at least the measured wall of every stage of every coupon; "
                         "SelfCheck). estimate_stages.py scales them to another hybrid prism-tube coupon mesh by the exact H1 DOF "
                         "ratio with 1x / 1.5x / 2x the measured mean PCG counts."),
             "MeasuredMesh": {"Label": (f"{reference['Case']} identity (device library 2026-09-22; "
                                        f"{reference['Elements']['Tetrahedron']:,} tets + {reference['Elements']['Prism']:,} prisms + "
                                        f"{reference['Elements']['Pyramid']:,} pyramids, {reference['EntityCounts']['Vertices']:,} nodes)"),
                              "SHA256": reference["MeshSHA256"], "EntityCounts": reference["EntityCounts"]},
             "MeasuredBlockSize": block_size}
    for key in KEPT_KEYS:
        model[key] = previous[key]
    model["ReducerEvaluationFractionRule"] = (previous["ReducerEvaluationFractionRule"]
                                             + f" KEPT by the {qualification_path.name} refit: one block size ({block_size}) cannot "
                                             "separate the evaluation and Gram parts; ReducerPairSeconds is inverted through this "
                                             "split at the reference source count so the 1x reducer estimate covers every coupon of the run")
    model["ReducerResidentFieldGBRule"] = (previous["ReducerResidentFieldGBRule"]
                                          + f" KEPT by the {qualification_path.name} refit as an upper bound: the streaming one-pass "
                                          "Gram executable keeps N x (4 Q_local + L_local) x 8 bytes resident whatever b, so any "
                                          f"b > {block_size} adds at most this much (not re-measured: every reducer ran at b = {block_size})")
    model["Stages"] = stages
    model["LocalEdge"] = local_edge
    model["RateRule"] = ("per stage and quantity the largest value over the run's coupons after scaling to the reference H1 "
                         "(the per-source non-solve seconds from the largest worker-block mean per-source seconds of the coupon, so "
                         "every block of a split coupon is covered; the worker's non-source seconds per worker job) "
                         "(times, peaks: x H1_ref / H1_c; node used GiB: the largest non-Palace remainder plus the reference peak; "
                         "the reducer block-pair seconds: divided by the estimator's evaluation / Gram scale of the coupon relative to "
                         "the reference source count; per-source archive GB: x H1 ratio / sources); the worker's non-source and the "
                         "reducer's setup seconds are wall-based (wall - sources x mean per-source seconds / wall - block-pair "
                         "seconds) so launch overhead is inside the estimate; MeanPerSourceSeconds = the largest per-source "
                         "non-solve seconds + the largest seconds per PCG iteration x the largest mean PCG count")
    model["Provenance"] = {"Qualification": {"Path": recorded_path(qualification_path), "SHA256": sha256(qualification_path),
                                             "ToolCommit": qualification.get("ToolCommit"), "Root": qualification.get("Root")},
                           "BuildRecord": {"Path": recorded_path(build_path), "SHA256": sha256(build_path),
                                           "Commit": build["Library"].get("Commit")},
                           "PBSJobs": pbs, "FrozenExecutable": executable, "ReferenceCase": reference["Case"],
                           "Coupons": {m["Case"]: {"EntityCounts": m["EntityCounts"], "MeshSHA256": m["MeshSHA256"],
                                                   "Stages": m["Stages"], "LocalEdge": m["LocalEdge"], "ArchiveGB": m["ArchiveGB"],
                                                   "Jobs": m["Jobs"]} for m in measurements},
                           "SetBy": provenance}
    model["Previous"] = {"Path": previous_kept.name, "SHA256": sha256(previous_kept),
                         "Provenance": ("the model every qualify run up to this refit used: " + previous["Purpose"])}
    check = self_check({**model, "Path": "refit", "SHA256": None}, measurements, block_size, profile)
    check_previous = self_check({**previous, "Path": str(previous_path), "SHA256": sha256(previous_path)}, measurements,
                                block_size, profile)
    model["SelfCheck"] = {"Rule": ("estimate_stages on every coupon of the run under this model, read per job as the split planner "
                                   "does: at the measured PCG counts every stage's estimate over its measured wall is >= 1 "
                                   "(MinimumStageEstimateOverActual) and so is every job's; the jobs at the worst PCG factor with "
                                   "preflight and margin over the coupon's measured node seconds is the conservatism ratio "
                                   "(EstimateWorstWithMarginOverActual); FitsOneJob is the single-job decision under this model; "
                                   "PreviousModel is the same check under the previous model"),
                          "MinimumStageEstimateOverActual": check["MinimumStageEstimateOverActual"],
                          "Conservative": check["Conservative"],
                          "EstimateWorstWithMarginOverActual": {case: check[case]["EstimateWorstWithMarginOverActualNodeSeconds"]
                                                                for case in check if case in model["Provenance"]["Coupons"]},
                          "MinimumJobEstimate1xOverStageWall": {case: check[case]["MinimumJobEstimate1xOverStageWall"]
                                                                for case in check if case in model["Provenance"]["Coupons"]},
                          "FitsOneJob": {case: check[case]["SingleJob"]["FitsOneJob"]
                                         for case in check if case in model["Provenance"]["Coupons"]},
                          "PreviousModel": {"MinimumStageEstimateOverActual": check_previous["MinimumStageEstimateOverActual"],
                                            "EstimateWorstWithMarginOverActual": {
                                                case: check_previous[case]["EstimateWorstWithMarginOverActualNodeSeconds"]
                                                for case in check_previous if case in model["Provenance"]["Coupons"]},
                                            "FitsOneJob": {case: check_previous[case]["SingleJob"]["FitsOneJob"]
                                                           for case in check_previous if case in model["Provenance"]["Coupons"]}}}
    if not check["Conservative"]:
        raise RefitError(f"the refit is not conservative: minimum stage estimate / actual {check['MinimumStageEstimateOverActual']:.3f}")
    record = {"Model": {key: model[key] for key in ("Version", "MeasuredMesh", "MeasuredBlockSize", "Stages", "LocalEdge", "Previous")},
              "SelfCheck": check, "SelfCheckPreviousModel": check_previous, "SetBy": provenance}
    return model, record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--qualification", type=Path, required=True, help="library-qualification.json of the completed run")
    parser.add_argument("--build-record", type=Path, help="library-build.json (default: next to the qualification, else its recorded path)")
    parser.add_argument("--previous", type=Path, default=estimate_stages.COST_MODEL, help="the model to replace")
    parser.add_argument("--previous-kept", type=Path, required=True, help="the byte-identical copy of --previous that stays in the tree")
    parser.add_argument("--reference-case", help="the coupon whose mesh the rates are expressed on (default: the largest H1)")
    parser.add_argument("--cluster-profile", type=Path, default=estimate_stages.CLUSTER_PROFILE)
    parser.add_argument("--out", type=Path, required=True, help="the refit cost model")
    parser.add_argument("--record", type=Path, required=True, help="the refit record (measurements, maxima, self-check)")
    args = parser.parse_args(argv)
    model, record = refit(args.qualification, build_record_path=args.build_record, previous_path=args.previous,
                          previous_kept=args.previous_kept, reference_case=args.reference_case,
                          profile=json.loads(args.cluster_profile.read_text()))
    args.out.write_text(json.dumps(model, indent=2) + "\n")
    args.record.write_text(json.dumps(record, indent=2) + "\n")
    estimate_stages.load_cost_model(args.out)
    check = model["SelfCheck"]
    print(f"refit {args.out}: reference {model['Provenance']['ReferenceCase']}, block size {model['MeasuredBlockSize']}, "
          f"minimum stage estimate / actual {check['MinimumStageEstimateOverActual']:.3f}; job estimate (worst PCG + margin) / "
          "actual: " + ", ".join(f"{case} {ratio:.2f}" for case, ratio in check["EstimateWorstWithMarginOverActual"].items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
