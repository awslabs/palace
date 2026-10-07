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
digest under Previous.  The ReducerPeakStreaming block (USER decision 2026-09-22 (B):
the reducer peak as a node-used line of the streaming executable, calibrated on measured
node-used peaks, not on the Palace peaks this refit reads) is carried from
--reducer-peak-streaming (default: the --previous model when it carries one) and must
name the run's frozen executable.

usage: refit_cost_model.py --qualification library-qualification.json [--build-record library-build.json]
       --previous cost-model.json --previous-kept PATH [--reducer-peak-streaming MODEL] [--reference-case CASE]
       --out cost-model.json --record refit.json
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


# A run's H1 may exceed the mesh entities' closed form by the cracked interior-boundary
# vertices of a thin coupon (~1 %); anything larger is another mesh.
H1_CLOSED_FORM_TOLERANCE = 0.02


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
        candidates = [record_dir / case / "results" / name / "status.json", root / "results" / "main" / relative / "status.json",
                      record_dir / Path(root).name / "results" / "main" / relative / "status.json"]
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
                           "MeanPerSourceSeconds": total / len(timings) if timings else None,
                           "MeanPCGIterations": sum(t["Iterations"] for t in timings) / len(timings) if timings else None,
                           "NodeUsedGiB": (stage.get("NodePeakUsedBytesSampled") or 0) / GIB,
                           "PalacePeakGB": memory_gb((stage["Parsed"].get("PalacePeakMemory") or {}).get("Total", "0G"))})
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
            "ReducerTimers": timers, "ReducerStageName": f"{prefix}-reducer",
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
            # A thin coupon's cracked interior boundary duplicates vertices at run time (a
            # few 0.1 % more DOFs than the mesh entities' closed form the planner estimates
            # with): the rates are scaled by the closed form the planner uses, the run's count
            # recorded; a larger difference is a different mesh.
            if abs(closed - measured["H1"]) > H1_CLOSED_FORM_TOLERANCE * closed:
                raise RefitError(f"{case}: closed-form H1 {closed} at p{item['Order']} differs from the run's {measured['H1']}")
            measured["RunH1"] = measured["H1"]
            measured["H1"] = closed
            measured["H1Rule"] = ("the closed form of the mesh entities (the planner's count); RunH1 is Palace's count of the run "
                                  "(the cracked interior boundary of a thin coupon duplicates vertices)")
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


# The per-block PCG / memory figures and the reducer timers were added to stage_measurement
# for the measured refit (decision 457 (2)); the decision-64a refit's provenance keeps its
# recorded shape (the committed device model is its byte-for-byte output).
LEGACY_BLOCK_KEYS = ("Name", "Job", "Sources", "WallSeconds", "NonSourceSeconds", "MeanPerSourceSeconds")
MEASURED_STAGE_KEYS = ("ReducerTimers", "ReducerStageName")


def legacy_stage_view(stage):
    view = {key: value for key, value in stage.items() if key not in MEASURED_STAGE_KEYS}
    view["WorkerBlocks"] = [{key: block[key] for key in LEGACY_BLOCK_KEYS} for block in stage["WorkerBlocks"]]
    return view


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


def refit(qualification_path, *, build_record_path=None, previous_path, previous_kept, streaming_path=None,
          reference_case=None, profile=None):
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
    streaming_source = json.loads(Path(streaming_path).read_text()) if streaming_path else previous
    streaming = streaming_source.get("ReducerPeakStreaming")
    if streaming_path and streaming is None:
        raise RefitError(f"--reducer-peak-streaming {streaming_path} carries no ReducerPeakStreaming block")
    if streaming is not None and streaming.get("Executable") != executable:
        raise RefitError(f"the ReducerPeakStreaming calibration names executable {streaming.get('Executable')}, "
                         f"the run's frozen executable is {executable}")
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
    if streaming is not None:
        # A node-used calibration of the streaming executable (measured node peaks, not the
        # Palace peaks this refit scales): carried as stated.
        model["ReducerPeakStreaming"] = {**streaming,
                                         "CarriedRule": (f"carried by the {qualification_path.name} refit unchanged: the block is a "
                                                         "node-used calibration of the streaming executable (USER decision "
                                                         "2026-09-22 (B), library-device-thin-01's reducers), not re-derivable from "
                                                         "the Palace peaks this refit reads; its Executable is checked against the "
                                                         "run's frozen executable; estimate_stages applies it to every reducer "
                                                         "stage of this model")}
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
                                                   "Stages": {order: legacy_stage_view(stage) for order, stage in m["Stages"].items()},
                                                   "LocalEdge": m["LocalEdge"], "ArchiveGB": m["ArchiveGB"],
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
    parser.add_argument("--qualification", type=Path, required=True, action="append",
                        help="library-qualification.json of a completed run (repeatable with --measured)")
    parser.add_argument("--build-record", type=Path, help="library-build.json (default: next to the qualification, else its recorded path)")
    parser.add_argument("--previous", type=Path, default=estimate_stages.COST_MODEL, help="the model to replace")
    parser.add_argument("--previous-kept", type=Path, required=True, help="the byte-identical copy of --previous that stays in the tree")
    parser.add_argument("--reducer-peak-streaming", type=Path,
                        help="the model whose ReducerPeakStreaming block is carried (default: --previous when it has one)")
    parser.add_argument("--reference-case", help="the coupon whose mesh the rates are expressed on (default: the largest H1)")
    parser.add_argument("--cluster-profile", type=Path, default=estimate_stages.CLUSTER_PROFILE)
    parser.add_argument("--measured", action="store_true",
                        help="the Version-3 refit on several runs (refit_measured: mean rates + per-stage safety factors, the reducer "
                             "line refit, the cap table; decision 457 (2) / 458)")
    parser.add_argument("--label", help="the label of the measured runs in the model's text (with --measured)")
    parser.add_argument("--out", type=Path, required=True, help="the refit cost model")
    parser.add_argument("--record", type=Path, required=True, help="the refit record (measurements, maxima, self-check)")
    args = parser.parse_args(argv)
    profile = json.loads(args.cluster_profile.read_text())
    if args.measured:
        model, record = refit_measured(args.qualification, previous_path=args.previous, previous_kept=args.previous_kept,
                                       reference_case=args.reference_case, profile=profile, label=args.label)
    else:
        if len(args.qualification) != 1:
            parser.error("one --qualification without --measured")
        model, record = refit(args.qualification[0], build_record_path=args.build_record, previous_path=args.previous,
                              previous_kept=args.previous_kept, streaming_path=args.reducer_peak_streaming,
                              reference_case=args.reference_case, profile=profile)
    args.out.write_text(json.dumps(model, indent=2) + "\n")
    args.record.write_text(json.dumps(record, indent=2) + "\n")
    estimate_stages.load_cost_model(args.out)
    check = model["SelfCheck"]
    print(f"refit {args.out}: reference {model['Provenance']['ReferenceCase']}, block size {model['MeasuredBlockSize']}, "
          f"minimum stage estimate / actual {check['MinimumStageEstimateOverActual']:.3f}; job estimate (worst PCG + margin) / "
          "actual: " + ", ".join(f"{case} {ratio:.2f}" for case, ratio in check["EstimateWorstWithMarginOverActual"].items()))
    if args.measured:
        caps = check["CapCheck"]
        print(f"cap check: minimum new cap / measured wall {caps['MinimumMarginNewCapOverWall']:.2f}; pessimism this model "
              f"{check['Pessimism']['ThisModel']}; previous {check['Pessimism']['PreviousModel']}")
    return 0




# ---------------------------------------------------------------------------------------
# Version-3 refit on several measured runs (decision 457 (2) / 458): mean rates with a
# per-stage safety factor from the measured residuals, the reducer node-used line refit,
# the old-vs-new cap table of every stored coupon and stage.
# ---------------------------------------------------------------------------------------
MODEL_VERSION_MEASURED = 3
# The measured policy (decision 457 (2)): the largest PCG factor covers the largest measured
# coupon-mean PCG ratio with this headroom; the runner overhead measured <= 0.7 % of the
# stage walls and <= 2.5 s of preflight on 4.7-9 M-element coupons - the margin keeps 10 %
# and 60 s for a 1 GB mesh with 1,000 traces.
PCG_HEADROOM = 1.2
MEASURED_MARGIN_FACTOR = 1.10
MEASURED_PREFLIGHT_SECONDS = 60
RATE_NAMES = ("SecondsPerPCGIteration", "PerSourceNonSolveSeconds", "WorkerNonSourceSeconds", "ReducerSetupSeconds",
              "ReducerPairSeconds")
PEAK_NAMES = ("WorkerPalacePeakGB", "ReducerPalacePeakGB", "WorkerNodeOverheadGiB", "ReducerNodeOverheadGiB")


def statistics(values):
    ordered = sorted(values)
    n = len(ordered)
    if n == 0:
        return {"Count": 0}
    median = ordered[n // 2] if n % 2 else 0.5 * (ordered[n // 2 - 1] + ordered[n // 2])
    return {"Count": n, "Min": ordered[0], "Median": median, "Mean": sum(ordered) / n, "Max": ordered[-1]}


def solve_normal_equations(rows, targets):
    """Least squares of targets ~ rows x coefficients (rows = list of equal-length lists) by
    the normal equations (Gaussian elimination with partial pivoting; stdlib only)."""
    k = len(rows[0])
    a = [[sum(row[i] * row[j] for row in rows) for j in range(k)] + [sum(row[i] * t for row, t in zip(rows, targets))] for i in range(k)]
    for column in range(k):
        pivot = max(range(column, k), key=lambda r: abs(a[r][column]))
        a[column], a[pivot] = a[pivot], a[column]
        if abs(a[column][column]) < 1e-300:
            raise RefitError("the reducer node-used fit is singular")
        for r in range(k):
            if r != column:
                f = a[r][column] / a[column][column]
                a[r] = [x - f * y for x, y in zip(a[r], a[column])]
    return [a[i][k] / a[i][i] for i in range(k)]


def fit_reducer_node_used(points):
    """NodeUsedGiB = a + (b + c x Sources) x H1 / 1e6 over the measured reducer stages
    (`points` = [(H1, Sources, NodeUsedGiB, label), ...]); a negative coefficient is fixed at
    0 and the rest refit (the line stays physical: a baseline, a per-DOF mesh / space term and
    the per-source resident fields).  The safety factor is the largest measured / fitted
    ratio (>= 1), so the line covers every measured reducer."""
    free = [0, 1, 2]
    coefficients = [0.0, 0.0, 0.0]
    while free:
        rows = [[[1.0, h1 / 1e6, sources * h1 / 1e6][i] for i in free] for h1, sources, _, _ in points]
        solved = solve_normal_equations(rows, [used for _, _, used, _ in points])
        coefficients = [0.0, 0.0, 0.0]
        for index, value in zip(free, solved):
            coefficients[index] = value
        negative = [index for index in free if coefficients[index] < 0.0]
        if not negative:
            break
        free = [index for index in free if index not in negative]
    a, b, c = coefficients

    def fitted(h1, sources):
        return a + (b + c * sources) * h1 / 1e6
    ratios = [(used / fitted(h1, sources), label) for h1, sources, used, label in points]
    safety = max(1.0, max(ratio for ratio, _ in ratios))
    return {"NodeBaselineGiB": round_up(a, 3), "NodeUsedGiBPerMillionH1": round_up(b, 5), "NodeUsedGiBPerMillionH1PerSource": round_up(c, 6),
            "SafetyFactor": round_up(safety, 3), "SafetyFactorSetBy": max(ratios)[1],
            "Residuals": statistics([ratio for ratio, _ in ratios]),
            "Measured": [{"Label": label, "H1": h1, "Sources": sources, "NodeUsedGiB": used, "FittedGiB": fitted(h1, sources),
                          "Ratio": used / fitted(h1, sources)} for h1, sources, used, label in points]}


def build_case_from_estimate(record_dir, case_record):
    """The build-record view of a case (H1 entity counts, element counts) from its recorded
    stage estimate when the run's library-build.json is not at hand (the stage-2 lanes keep
    it under the registration root of another machine)."""
    path = record_dir / Path(case_record["Root"]).name / "preflight" / "stage-estimate.json"
    if not path.exists():
        raise RefitError(f"{case_record['Case']}: neither a build record nor a stage estimate at {path}")
    mesh = json.loads(path.read_text())["Mesh"]
    counts = mesh["EntityCounts"]
    return {"Case": case_record["Case"], "H1": {"EntityCounts": counts},
            "Elements": {"Tetrahedron": counts["Tetrahedra"], "Prism": counts["Prisms"], "Pyramid": counts["Pyramids"],
                         "Total": mesh["VolumeElements"]}}


def mean_rate_stage(order, measurements, reference, previous):
    """The Version-3 stage record at `order`: every time rate the MEAN of the scaled
    measurements (memory: the largest, fail-closed), the PCG counts the mean / max over
    the coupons; the time safety factors from the residuals are set by stage_safety."""
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
    scaled = {name: [] for name in RATE_NAMES + PEAK_NAMES + ("ArchiveGBPerSource",)}
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
    chosen = {}
    for name in RATE_NAMES:
        values = [value for value, _ in scaled[name]]
        chosen[name] = {"Value": sum(values) / len(values), "Rule": "mean", "Statistics": statistics(values),
                        "Largest": list(_largest(scaled[name]))}
    for name in PEAK_NAMES + (("ArchiveGBPerSource",) if scaled["ArchiveGBPerSource"] else ()):
        value, case = _largest(scaled[name])
        chosen[name] = {"Value": value, "Rule": "largest", "SetBy": case, "Statistics": statistics([v for v, _ in scaled[name]])}
    chosen["MeanPCGIterations"] = {"Value": sum(v for v, _ in pcg_mean) / len(pcg_mean), "Rule": "mean of the coupon means",
                                   "Statistics": statistics([v for v, _ in pcg_mean])}
    chosen["MaxPCGIterations"] = {"Value": _largest(pcg_max)[0], "Rule": "largest", "SetBy": _largest(pcg_max)[1]}
    spi, mean_pcg = chosen["SecondsPerPCGIteration"]["Value"], chosen["MeanPCGIterations"]["Value"]
    worker_peak, reducer_peak = chosen["WorkerPalacePeakGB"]["Value"], chosen["ReducerPalacePeakGB"]["Value"]
    reducer_pair = chosen["ReducerPairSeconds"]["Value"]
    if "ArchiveGBPerSource" in chosen:
        archive = chosen["ArchiveGBPerSource"]["Value"] * ref_sources
    else:
        prev = previous["Stages"][key]
        archive = prev["ArchiveGB"] * ref_h1 / prev["H1"] * ref_sources / prev["Sources"]
        chosen["ArchiveGBPerSource"] = {"Value": archive / ref_sources, "Rule": "previous model (scaled)"}
    stage = {"H1": ref_h1, "Sources": ref_sources, "SecondsPerPCGIteration": round_up(spi, 5),
             "MeanPCGIterations": round_up(mean_pcg, 3), "MaxPCGIterations": int(chosen["MaxPCGIterations"]["Value"]),
             "MeanPerSourceSeconds": round_up(chosen["PerSourceNonSolveSeconds"]["Value"] + spi * mean_pcg, 3),
             "WorkerNonSourceSeconds": round_up(chosen["WorkerNonSourceSeconds"]["Value"], 2),
             "ReducerPalaceSeconds": round_up(chosen["ReducerSetupSeconds"]["Value"] + reducer_pair, 2),
             "ReducerPairSeconds": round_up(reducer_pair, 2),
             "WorkerPalacePeakGB": round_up(worker_peak, 1), "ReducerPalacePeakGB": round_up(reducer_peak, 1),
             "WorkerNodeUsedGiB": round_up(chosen["WorkerNodeOverheadGiB"]["Value"] + worker_peak / gb_per_gib, 2),
             "ReducerNodeUsedGiB": round_up(chosen["ReducerNodeOverheadGiB"]["Value"] + reducer_peak / gb_per_gib, 2),
             "ArchiveGB": round_up(archive, 2)}
    return stage, chosen, block_size


def stage_safety(order, stage, measurements, previous, block_size):
    """The per-stage time safety factors of a Version-3 stage: the largest measured wall /
    estimate at the MEAN rates and the coupon's OWN mean PCG count over every measured worker
    block and reducer (so the PCG spread stays with the PCGFactors), floored at 1; the residual
    distributions recorded.  Memory residuals (measured / estimate) are recorded for the record."""
    key = f"p{order}"
    fraction = previous["ReducerEvaluationFraction"]
    gb_per_gib = previous["PalaceGBPerGiB"]
    ref_evaluations = estimate_stages.source_evaluations(stage["Sources"], block_size)
    ref_pairs = stage["Sources"] * (stage["Sources"] + 1) // 2
    worker, reducer, worker_memory, reducer_memory = [], [], [], []
    per_source_non_solve = stage["MeanPerSourceSeconds"] - stage["SecondsPerPCGIteration"] * stage["MeanPCGIterations"]
    for m in measurements:
        s = m["Stages"][key]
        ratio = s["H1"] / stage["H1"]
        for block in s["WorkerBlocks"]:
            if not block["Sources"]:
                continue
            pcg = block["MeanPCGIterations"] if block["MeanPCGIterations"] is not None else s["MeanPCGIterations"]
            estimate = (stage["WorkerNonSourceSeconds"] + block["Sources"] * (per_source_non_solve + stage["SecondsPerPCGIteration"] * pcg)) * ratio
            worker.append((block["WallSeconds"] / estimate, f"{m['Case']}:{block['Name']}"))
            if block["NodeUsedGiB"]:
                worker_memory.append((block["NodeUsedGiB"] / (stage["WorkerNodeUsedGiB"] - stage["WorkerPalacePeakGB"] / gb_per_gib
                                                              + stage["WorkerPalacePeakGB"] * ratio / gb_per_gib), f"{m['Case']}:{block['Name']}"))
        evaluations = estimate_stages.source_evaluations(s["Sources"], block_size)
        pairs = s["Sources"] * (s["Sources"] + 1) // 2
        scale = fraction * evaluations / ref_evaluations + (1.0 - fraction) * pairs / ref_pairs
        estimate = ((stage["ReducerPalaceSeconds"] - stage["ReducerPairSeconds"]) + stage["ReducerPairSeconds"] * scale) * ratio
        reducer.append((s["ReducerWallSeconds"] / estimate, f"{m['Case']}:{s['ReducerStageName']}"))
        reducer_memory.append((s["ReducerNodeUsedGiB"], f"{m['Case']}"))
    return {"SafetyFactor": {"Worker": round_up(max(1.0, max(r for r, _ in worker)), 3),
                             "Reducer": round_up(max(1.0, max(r for r, _ in reducer)), 3)},
            "SafetyFactorSetBy": {"Worker": max(worker)[1], "Reducer": max(reducer)[1]},
            "Residuals": {"Worker": statistics([r for r, _ in worker]), "Reducer": statistics([r for r, _ in reducer]),
                          "WorkerNodeUsedGiB": statistics([r for r, _ in worker_memory])},
            "ResidualRule": ("measured wall / the estimate at the model's mean rates and the block's own mean PCG count, per "
                             "measured worker block and reducer; SafetyFactor = the largest (>= 1): the 1x estimate at a coupon's "
                             "own PCG count covers every measured stage; the PCG spread is the PCGFactors' business")}


def local_edge_mean(measurements, reference):
    """The Version-3 local-edge record: mean solve / non-solve rates, the largest peaks, the
    safety factor from the residuals."""
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
    solve = sum(v for v, _ in scaled["LinearSolveSeconds"]) / len(measurements)
    non_solve = sum(v for v, _ in scaled["NonSolveSeconds"]) / len(measurements)
    peak, node_used = _largest(scaled["PalacePeakGB"]), _largest(scaled["NodeUsedGiB"])
    residuals = []
    for m in measurements:
        local = m["LocalEdge"]
        ratio = local["H1"] / ref_h1
        residuals.append((local["WallSeconds"] / ((non_solve + solve) * ratio), m["Case"]))
    stage = {"Order": order, "H1": ref_h1, "Sources": sources.pop(), "PalacePeakGB": round_up(peak[0], 1),
             "NodeUsedGiB": round_up(node_used[0], 1), "PalaceTotalSeconds": round_up(solve + non_solve, 1),
             "LinearSolveSeconds": round_up(solve, 1), "NonSolveSeconds": round_up(non_solve, 1),
             "SafetyFactor": {"LocalEdge": round_up(max(1.0, max(r for r, _ in residuals)), 3)},
             "SafetyFactorSetBy": max(residuals)[1], "Residuals": {"LocalEdge": statistics([r for r, _ in residuals])}}
    provenance = {"LinearSolveSeconds": {"Value": solve, "Rule": "mean", "Statistics": statistics([v for v, _ in scaled["LinearSolveSeconds"]])},
                  "NonSolveSeconds": {"Value": non_solve, "Rule": "mean", "Statistics": statistics([v for v, _ in scaled["NonSolveSeconds"]])},
                  "PalacePeakGB": {"Value": peak[0], "Rule": "largest", "SetBy": peak[1]},
                  "NodeUsedGiB": {"Value": node_used[0], "Rule": "largest", "SetBy": node_used[1]}}
    return stage, provenance


def cap_check(model, measurements, block_size, profile, qualification_cases):
    """The old-vs-new cap table (decision 458): per measured coupon and stage the measured
    wall, the cap of record (the run's plan, library-qualification Plan.Caps) and the cap
    under `model` (build_plan.stage_caps on the new estimate; a block worker at its block's
    source count), the margin new cap / wall; the minimum margin must be >= 1."""
    import build_plan
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    deadline = profile["DeadlineSeconds"] - profile["DeadlineMarginSeconds"]
    rows = []
    for m in measurements:
        case_record = qualification_cases[m["Case"]]
        caps_of_record = (case_record.get("Plan") or {}).get("Caps") or {}
        stages = [(f"{order}-{s['Sources']}", int(order[1:]), s["Sources"]) for order, s in m["Stages"].items()]
        local = m["LocalEdge"]
        estimate = estimate_stages.estimate(m["EntityCounts"], stages, model=model, profile=profile, block_size=block_size,
                                            local_edge=(local["Name"], local["Order"], local["Sources"]))
        for order, s in m["Stages"].items():
            est = estimate["Stages"][f"{order}-{s['Sources']}"]
            for block in s["WorkerBlocks"]:
                if not block["Sources"]:
                    continue
                view = build_plan.block_estimate(est, block["Sources"], factors) if block["Sources"] != s["Sources"] else est
                cap, minimum = build_plan.stage_caps(view, "worker", factors, deadline)
                rows.append({"Case": m["Case"], "Stage": block["Name"], "Kind": "worker", "Sources": block["Sources"],
                             "MeasuredWallSeconds": block["WallSeconds"], "CapOfRecordSeconds": (caps_of_record.get(block["Name"]) or [None])[0],
                             "NewCapSeconds": cap, "NewMinimumSeconds": minimum})
            cap, minimum = build_plan.stage_caps(est, "reducer", factors, deadline)
            rows.append({"Case": m["Case"], "Stage": s["ReducerStageName"], "Kind": "reducer", "Sources": s["Sources"],
                         "MeasuredWallSeconds": s["ReducerWallSeconds"], "CapOfRecordSeconds": (caps_of_record.get(s["ReducerStageName"]) or [None])[0],
                         "NewCapSeconds": cap, "NewMinimumSeconds": minimum})
        cap, minimum = build_plan.stage_caps(estimate["Stages"][local["Name"]], "local-edge", factors, deadline)
        rows.append({"Case": m["Case"], "Stage": local["Name"], "Kind": "local-edge", "Sources": local["Sources"],
                     "MeasuredWallSeconds": local["WallSeconds"], "CapOfRecordSeconds": (caps_of_record.get(local["Name"]) or [None])[0],
                     "NewCapSeconds": cap, "NewMinimumSeconds": minimum})
    for row in rows:
        row["MarginNewCapOverWall"] = row["NewCapSeconds"] / row["MeasuredWallSeconds"] if row["MeasuredWallSeconds"] else None
        row["NewCapOverCapOfRecord"] = (row["NewCapSeconds"] / row["CapOfRecordSeconds"]) if row["CapOfRecordSeconds"] else None
    margins = [row["MarginNewCapOverWall"] for row in rows if row["MarginNewCapOverWall"] is not None]
    below = [row for row in rows if row["MarginNewCapOverWall"] is not None and row["MarginNewCapOverWall"] < 1.0]
    return {"Rule": ("per measured coupon and stage: the measured wall, the cap of record (the run's plan) and the cap under the new "
                     "model (CapSeconds = 2 x the estimate at the largest PCG factor, a block worker at its block's source count); "
                     "a stage whose new cap is below its measured wall fails the refit (decision 458)"),
            "Rows": rows, "MinimumMarginNewCapOverWall": min(margins) if margins else None,
            "MarginStatistics": statistics(margins), "StagesBelowNewCap": [row["Stage"] for row in below],
            "NewCapOverCapOfRecord": statistics([row["NewCapOverCapOfRecord"] for row in rows if row["NewCapOverCapOfRecord"]])}


def refit_measured(qualification_paths, *, previous_path, previous_kept, reference_case=None, profile=None, label=None):
    """The Version-3 model from several completed runs (decision 457 (2) / 458): every time
    rate the mean of the scaled measurements with a per-stage SafetyFactor from the residuals
    (estimate_stages applies it), the memory peaks the largest (fail-closed), the reducer
    node-used line (ReducerPeakStreaming) refit on every measured reducer with its safety
    factor from the residuals, the previous model kept and bound by digest, the self-check,
    the pessimism of both models and the old-vs-new cap table."""
    paths = [Path(path) for path in qualification_paths]
    previous = json.loads(Path(previous_path).read_text())
    previous_kept = Path(previous_kept)
    if sha256(previous_kept) != sha256(previous_path):
        raise RefitError(f"the kept previous model {previous_kept} is not byte-identical to {previous_path}")
    profile = profile or json.loads(estimate_stages.CLUSTER_PROFILE.read_text())
    measurements, cases, runs, executables = [], {}, [], {}
    for qualification_path in paths:
        qualification = json.loads(qualification_path.read_text())
        record_dir = qualification_path.parent
        build_path = record_dir / "library-build.json"
        if not build_path.exists():
            build_path = Path(qualification["BuildRecord"]["Path"])
        build_cases = {}
        if build_path.exists():
            build_cases = {case["Case"]: case for case in json.loads(build_path.read_text())["Cases"]}
        for case in qualification["Cases"]:
            if not (case.get("Cost") and case["Cost"].get("Jobs")):
                continue
            if case["Case"] in cases:
                raise RefitError(f"{case['Case']} is measured by two runs")
            build_case = build_cases.get(case["Case"]) or build_case_from_estimate(record_dir, case)
            measurement = measure_case(record_dir, case, build_case)
            measurement["Run"] = recorded_path(qualification_path)
            measurements.append(measurement)
            cases[case["Case"]] = case
            executable = (case.get("FrozenExecutable") or {}).get("SHA256")
            executables.setdefault(executable, []).append(case["Case"])
        runs.append({"Path": recorded_path(qualification_path), "SHA256": sha256(qualification_path),
                     "ToolCommit": qualification.get("ToolCommit"), "Root": qualification.get("Root"),
                     "BuildRecord": ({"Path": recorded_path(build_path), "SHA256": sha256(build_path)} if build_path.exists() else
                                     {"Path": str(build_path), "Absent": "entity counts read from each case's recorded stage estimate"})})
    if not measurements:
        raise RefitError("no coupon of the runs has a complete cost record")
    if reference_case is None:
        reference = max(measurements, key=lambda m: max(s["H1"] for s in m["Stages"].values()))
    else:
        reference = next((m for m in measurements if m["Case"] == reference_case), None)
        if reference is None:
            raise RefitError(f"--reference-case {reference_case} is not a measured coupon of the runs")
    orders = sorted({order for m in measurements for order in m["Stages"]})
    stages, provenance, block_sizes = {}, {}, set()
    for order in orders:
        with_stage = [m for m in measurements if order in m["Stages"]]
        stages[order], provenance[order], block_size = mean_rate_stage(int(order[1:]), with_stage, reference, previous)
        block_sizes.add(block_size)
    if len(block_sizes) != 1:
        raise RefitError(f"the stages ran the reducer at different block sizes {sorted(block_sizes)}")
    block_size = block_sizes.pop()
    for order in orders:
        with_stage = [m for m in measurements if order in m["Stages"]]
        stages[order].update(stage_safety(int(order[1:]), stages[order], with_stage, previous, block_size))
    local_edge, provenance["LocalEdge"] = local_edge_mean(measurements, reference)
    reducer_points = [(s["H1"], s["Sources"], s["ReducerNodeUsedGiB"], f"{m['Case']}:{order}")
                      for m in measurements for order, s in m["Stages"].items()]
    streaming_fit = fit_reducer_node_used(reducer_points)
    pbs = sorted(job["PBSJobID"].split(".")[0] for m in measurements for job in m["Jobs"].values() if job.get("PBSJobID"))
    label = label or "stage-2"
    model = {"Version": MODEL_VERSION_MEASURED,
             "Copyright": previous["Copyright"], "SPDX-License-Identifier": previous["SPDX-License-Identifier"],
             "Purpose": ("Measured per-stage cost rates of the Palace response worker / reducer / ordinary local-edge stages on one "
                         f"node (192 ranks, frozen executables {sorted(executables)}, streaming one-pass Gram reducer at block size "
                         f"{block_size}) refit by refit_cost_model.refit_measured from {len(measurements)} measured coupons of "
                         f"{len(paths)} qualify runs ({label}; PBS {pbs[0]}..{pbs[-1]}). Every time rate is the MEAN of the "
                         f"measurements scaled to the reference mesh {reference['Case']} by the exact H1 ratio and every stage carries a "
                         "SafetyFactor = the largest measured wall / mean-rate estimate at the coupon's own PCG count (decision 458: "
                         "a per-stage factor from the residuals, not a constant); memory peaks are the largest scaled measurement; the "
                         "reducer node-used line is refit on every measured reducer with its own residual safety factor. "
                         "estimate_stages.py scales the rates to another coupon mesh by the exact H1 ratio with 1x / 1.5x / 2x the "
                         "mean PCG counts and applies the safety factors; the cap of every measured stage of record is above its "
                         "measured wall (CapCheck)."),
             "MeasuredMesh": {"Label": (f"{reference['Case']} identity ({label}; "
                                        f"{reference['Elements']['Tetrahedron']:,} tets + {reference['Elements']['Prism']:,} prisms + "
                                        f"{reference['Elements']['Pyramid']:,} pyramids, {reference['EntityCounts']['Vertices']:,} nodes)"),
                              "SHA256": reference["MeshSHA256"], "EntityCounts": reference["EntityCounts"]},
             "MeasuredBlockSize": block_size}
    for key in KEPT_KEYS:
        model[key] = previous[key]
    # The policy factors are recalibrated on the measured runs too (decision 457 (2)): the
    # PCG spread that matters for a stage's wall is the spread of the coupon MEAN counts
    # (the single-source maxima the previous 2x covered average out over a stage), and
    # the runner's preflight / between-stage overhead is measured per job.
    pcg_ratios = {order: max(m["Stages"][order]["MeanPCGIterations"] for m in measurements if order in m["Stages"])
                  / stages[order]["MeanPCGIterations"] for order in orders}
    largest_pcg_ratio = max(pcg_ratios.values())
    worst_factor = max(1.5, round_up(PCG_HEADROOM * largest_pcg_ratio, 1))
    model["PCGFactors"] = [1.0, round_up((1.0 + worst_factor) / 2.0, 1), worst_factor]   # one decimal: the factors key ByPCGFactor as "%.1f"
    overheads = [job["ActualSeconds"] / sum(stage["WallSeconds"] for stage in job["Stages"]) for m in measurements for job in m["Jobs"].values()
                 if job.get("ActualSeconds") and job["Stages"]]
    model["PreflightAndMarginFactor"] = MEASURED_MARGIN_FACTOR
    model["PreflightSeconds"] = MEASURED_PREFLIGHT_SECONDS
    model["PolicyRule"] = {
        "PCGFactors": (f"the largest factor = max(1.5, {PCG_HEADROOM} x the largest measured coupon-mean PCG count over the model mean, "
                       f"rounded up to 0.1) = {worst_factor} (measured ratios per order {json.dumps({k: round(v, 3) for k, v in pcg_ratios.items()})}); "
                       f"the middle factor halfway rounded up to 0.1; previous {previous['PCGFactors']} (set when the single-source maximum, not the stage "
                       "mean, was the reference)"),
        "PreflightAndMarginFactor": (f"{MEASURED_MARGIN_FACTOR}: the measured runner overhead (job total / the sum of its stage walls) is "
                                     f"{statistics(overheads)} over {len(overheads)} jobs; previous {previous['PreflightAndMarginFactor']}"),
        "PreflightSeconds": (f"{MEASURED_PREFLIGHT_SECONDS}: the measured preflight (mesh + trace digests, node checks) is "
                             f"{round(max(job['ActualSeconds'] - sum(stage['WallSeconds'] for stage in job['Stages']) for m in measurements for job in m['Jobs'].values() if job.get('ActualSeconds') and job['Stages']), 1)} s "
                             f"at most on 4.7-9 M-element coupons; previous {previous['PreflightSeconds']}"),
        "MeasuredOverheads": statistics(overheads)}
    model["ReducerEvaluationFractionRule"] = (previous["ReducerEvaluationFractionRule"].split(" KEPT by")[0]
                                             + f" KEPT by the {label} refit: one block size ({block_size}) cannot separate the "
                                             "evaluation and Gram parts")
    model["ReducerResidentFieldGBRule"] = (previous["ReducerResidentFieldGBRule"].split(" KEPT by")[0]
                                          + f" KEPT by the {label} refit as an upper bound (every reducer ran at b = {block_size})")
    model["ReducerPeakStreaming"] = {
        "Executable": sorted(executables),
        **{key: streaming_fit[key] for key in ("NodeBaselineGiB", "NodeUsedGiBPerMillionH1", "NodeUsedGiBPerMillionH1PerSource", "SafetyFactor")},
        "Rule": ("ReducerPalacePeakGBEstimate of every reducer stage: SafetyFactor x PalaceGBPerGiB x (NodeBaselineGiB + "
                 "(NodeUsedGiBPerMillionH1 + NodeUsedGiBPerMillionH1PerSource x sources) x H1 / 1e6), a node-used line fit by least "
                 f"squares on the {len(reducer_points)} measured reducer node-used peaks of the {label} runs (every order, 8-source "
                 "controls and the main stages); SafetyFactor = the largest measured / fitted ratio (decision 458), so the line covers "
                 "every measured reducer; the per-source term is the resident archived fields (~8 bytes per H1 DOF per source)"),
        "SafetyFactorSetBy": streaming_fit["SafetyFactorSetBy"], "Residuals": streaming_fit["Residuals"],
        "Measured": streaming_fit["Measured"],
        "Previous": {"Executable": previous["ReducerPeakStreaming"].get("Executable"),
                     **{key: previous["ReducerPeakStreaming"][key] for key in ("NodeBaselineGiB", "NodeUsedGiBPerMillionH1",
                                                                                "NodeUsedGiBPerMillionH1PerSource", "SafetyFactor")},
                     "Rule": "the line of the previous model (library-device-thin-01's reducers, USER decision 2026-09-22 (B))"}}
    model["Stages"] = stages
    model["LocalEdge"] = local_edge
    model["RateRule"] = ("per stage and time quantity the MEAN over the measured coupons after scaling to the reference H1 (times: x "
                        "H1_ref / H1_c; the reducer block-pair seconds divided by the estimator's evaluation / Gram scale of the coupon "
                        "relative to the reference source count), the PCG counts the mean of the coupon means (MaxPCGIterations the "
                        "largest); memory: the largest scaled Palace peak, the largest non-Palace node remainder plus the reference "
                        "peak, the largest per-source archive GB; SafetyFactor per stage (Worker / Reducer / LocalEdge) = the largest "
                        "measured wall / mean-rate estimate at the coupon's own PCG count (>= 1), applied by estimate_stages to the "
                        "stage's times; the worker's non-source and the reducer's setup seconds are wall-based so launch overhead is "
                        "inside the estimate")
    model["Provenance"] = {"Runs": runs, "PBSJobs": pbs, "FrozenExecutables": executables, "ReferenceCase": reference["Case"],
                           "Coupons": {m["Case"]: {"Run": m["Run"], "EntityCounts": m["EntityCounts"], "MeshSHA256": m["MeshSHA256"],
                                                   "Stages": {order: {key: value for key, value in s.items() if key != "ReducerTimers"}
                                                              for order, s in m["Stages"].items()},
                                                   "LocalEdge": m["LocalEdge"], "ArchiveGB": m["ArchiveGB"], "Jobs": m["Jobs"]}
                                       for m in measurements},
                           "SetBy": provenance}
    model["Previous"] = {"Path": previous_kept.name, "SHA256": sha256(previous_kept),
                         "Provenance": ("the model every qualify run up to this refit used: " + previous["Purpose"])}
    check = self_check({**model, "Path": "refit", "SHA256": None}, measurements, block_size, profile)
    check_previous = self_check({**previous, "Path": str(previous_path), "SHA256": sha256(previous_path)}, measurements,
                                block_size, profile)
    caps = cap_check({**model, "Path": "refit", "SHA256": None}, measurements, block_size, profile, cases)
    pessimism = {case: check[case]["EstimateWorstWithMarginOverActualNodeSeconds"] for case in check if case in cases}
    pessimism_previous = {case: check_previous[case]["EstimateWorstWithMarginOverActualNodeSeconds"] for case in check_previous if case in cases}
    model["SelfCheck"] = {"Rule": ("estimate_stages on every coupon of the runs under this model, read per job as the split planner "
                                   "does: at the measured PCG counts every job's estimate over its measured stage walls "
                                   "(MinimumJobEstimate1xOverStageWall) and the stage estimates over the walls "
                                   "(MinimumStageEstimateOverActual: below 1 only where a coupon's own PCG count exceeds the mean - "
                                   "the per-stage SafetyFactor is set at the coupon's own count, see ResidualRule); the jobs at the "
                                   "worst PCG factor with preflight and margin over the coupon's measured node seconds is the "
                                   "PESSIMISM (EstimateWorstWithMarginOverActual: the stated pessimism of this model and of the "
                                   "previous one); CapCheck: no measured stage of record above its new cap"),
                          "MinimumStageEstimateOverActual": check["MinimumStageEstimateOverActual"],
                          "Conservative": check["Conservative"],
                          "ConservativeAtOwnPCG": all(
                              stages[order]["Residuals"][kind]["Max"] <= stages[order]["SafetyFactor"][kind] + 1e-9
                              for order in orders for kind in ("Worker", "Reducer"))
                          and local_edge["Residuals"]["LocalEdge"]["Max"] <= local_edge["SafetyFactor"]["LocalEdge"] + 1e-9,
                          "ConservativeAtOwnPCGRule": ("every measured worker block, reducer and local-edge stage: measured wall <= the "
                                                       "estimate at the mean rates and its own PCG count x the stage's SafetyFactor"),
                          "EstimateWorstWithMarginOverActual": pessimism,
                          "Pessimism": {"ThisModel": statistics(list(pessimism.values())),
                                        "PreviousModel": statistics(list(pessimism_previous.values()))},
                          "MinimumJobEstimate1xOverStageWall": {case: check[case]["MinimumJobEstimate1xOverStageWall"] for case in pessimism},
                          "FitsOneJob": {case: check[case]["SingleJob"]["FitsOneJob"] for case in pessimism},
                          "PreviousModel": {"MinimumStageEstimateOverActual": check_previous["MinimumStageEstimateOverActual"],
                                            "EstimateWorstWithMarginOverActual": pessimism_previous,
                                            "FitsOneJob": {case: check_previous[case]["SingleJob"]["FitsOneJob"] for case in pessimism_previous}},
                          "CapCheck": {key: caps[key] for key in ("Rule", "MinimumMarginNewCapOverWall", "MarginStatistics",
                                                                  "StagesBelowNewCap", "NewCapOverCapOfRecord")}}
    if caps["StagesBelowNewCap"]:
        raise RefitError(f"the refit fails: {len(caps['StagesBelowNewCap'])} measured stages are above their new cap "
                         f"(minimum margin {caps['MinimumMarginNewCapOverWall']:.3f}): {caps['StagesBelowNewCap'][:6]}")
    record = {"Model": {key: model[key] for key in ("Version", "MeasuredMesh", "MeasuredBlockSize", "Stages", "LocalEdge",
                                                     "ReducerPeakStreaming", "Previous")},
              "SelfCheck": check, "SelfCheckPreviousModel": check_previous, "CapCheck": caps, "SetBy": provenance}
    return model, record


if __name__ == "__main__":
    raise SystemExit(main())
