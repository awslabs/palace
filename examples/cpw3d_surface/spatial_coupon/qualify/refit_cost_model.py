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


def load_status(record_dir, case_record):
    """The run_stages status carrying the local-edge stage (the controls job's) and every
    stage's elapsed-time report: copied under <case>/results/<job>/status.json, else the
    run's results tree."""
    case = case_record["Case"]
    jobs = case_record["Cost"]["Jobs"]
    candidates = [record_dir / case / "results" / name / "status.json" for name in jobs]
    root = Path(case_record["Root"])
    candidates.append(root / "results" / "main" / "status.json")
    for name in jobs:
        candidates.append(root / "results" / "main" / name / "status.json")
    statuses = [json.loads(path.read_text()) for path in candidates if path.exists()]
    if not statuses:
        raise RefitError(f"{case}: no run_stages status.json under {record_dir / case / 'results'} or {root / 'results'}")
    return statuses


def stage_measurement(cost_stage, *, wall_based=True):
    """The measured quantities of one worker + reducer stage of summarize_cost."""
    n = len(cost_stage["Sources"])
    mean_total = cost_stage["MeanPerSourceTotalSeconds"]
    worker_wall = cost_stage["WorkerWallSeconds"]
    reducer_wall = cost_stage["ReducerWallSeconds"]
    pair = cost_stage["ReducerPairSecondsEstimate"]
    return {"H1": cost_stage["H1"], "Sources": n, "Order": cost_stage["Order"],
            "SecondsPerPCGIteration": cost_stage["SecondsPerPCGIteration"],
            "MeanPCGIterations": cost_stage["MeanPCGIterations"], "MaxPCGIterations": cost_stage["MaxPCGIterations"],
            "MeanPerSourceSeconds": mean_total,
            "PerSourceNonSolveSeconds": mean_total - cost_stage["SecondsPerPCGIteration"] * cost_stage["MeanPCGIterations"],
            "WorkerWallSeconds": worker_wall, "WorkerNonSourceSeconds": max(worker_wall - n * mean_total, 0.0),
            "ReducerWallSeconds": reducer_wall, "ReducerPairSeconds": pair,
            "ReducerSetupSeconds": max(reducer_wall - pair, 0.0),
            "ReducerBlockSize": cost_stage["ReducerBlockSize"], "StageWallSeconds": cost_stage["StageWallSeconds"],
            "WorkerPalacePeakGB": memory_gb(cost_stage["Memory"]["WorkerPalacePeakTotal"]),
            "ReducerPalacePeakGB": memory_gb(cost_stage["Memory"]["ReducerPalacePeakTotal"]),
            "WorkerNodeUsedGiB": cost_stage["Memory"]["WorkerNodePeakUsedGiBSampled"],
            "ReducerNodeUsedGiB": cost_stage["Memory"]["ReducerNodePeakUsedGiBSampled"]}


def local_edge_measurement(statuses, cost_summary):
    name, stage = next(((name, stage) for name, stage in cost_summary["Stages"].items() if name.endswith("-local-edge")), (None, None))
    if stage is None:
        raise RefitError("no local-edge stage in the cost summary")
    parsed = None
    node_used = None
    for status in statuses:
        for item in status["Stages"]:
            if item["Name"] == name:
                parsed = item["Parsed"]
                node_used = (item.get("NodePeakUsedBytesSampled") or 0) / GIB
    if parsed is None:
        raise RefitError(f"{name}: not in any status.json")
    timers = parse_elapsed_time_report(parsed["ElapsedTimeReport"])
    solve = linear_solve_seconds(timers)
    wall = stage["WallSeconds"]
    return {"Name": name, "Order": parsed["Order"], "H1": parsed["H1"], "Sources": len(parsed["Iterations"]),
            "WallSeconds": wall, "PalaceTotalSeconds": parsed["PalaceTotalSeconds"],
            "LinearSolveSeconds": solve, "NonSolveSeconds": max(wall - solve, 0.0),
            "PalacePeakGB": memory_gb(parsed["PalacePeakMemory"]["Total"]), "NodeUsedGiB": node_used,
            "PCG": parsed.get("PCG")}


def measure_case(record_dir, case_record, build_case):
    case = case_record["Case"]
    cost_summary = json.loads(resolve(record_dir, case, "cost-summary.json", case_record["Cost"]["Path"]).read_text())
    statuses = load_status(record_dir, case_record)
    counts = build_case["H1"]["EntityCounts"]
    stages = {}
    for item in case_record["Stages"]:
        if item["Kind"] != "response":
            continue
        cost_stage = cost_summary["Stages"].get(item["Prefix"])
        if not cost_stage or "H1" not in cost_stage:
            raise RefitError(f"{case}: stage {item['Prefix']} has no complete cost record")
        measured = stage_measurement(cost_stage)
        closed = h1_dofs_from_counts(counts, item["Order"])
        if closed != measured["H1"]:
            raise RefitError(f"{case}: closed-form H1 {closed} at p{item['Order']} differs from the run's {measured['H1']}")
        measured["Role"] = item["Role"]
        stages[f"p{item['Order']}"] = measured
    local = local_edge_measurement(statuses, cost_summary)
    if h1_dofs_from_counts(counts, local["Order"]) != local["H1"]:
        raise RefitError(f"{case}: closed-form H1 at the local-edge order differs from the run's {local['H1']}")
    return {"Case": case, "EntityCounts": counts, "MeshSHA256": case_record["Mesh"]["SHA256"],
            "Elements": build_case["Elements"], "Stages": stages, "LocalEdge": local,
            "ArchiveGB": archive_gb(case_record.get("ArchiveDeletion")),
            "Jobs": {name: {"PBSJobID": job["PBSJobID"], "ActualSeconds": job["ActualSeconds"],
                            "EstimateSecondsWithPreflightAndMargin": job["EstimateSecondsWithPreflightAndMargin"]}
                     for name, job in case_record["Cost"]["Jobs"].items()}}


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
    """estimate_stages on every coupon of the run under `model`: the 1.0x estimate of
    every stage against its measured wall (>= 1 required), the job at the worst PCG
    factor with preflight and margin against the sum of the measured stage walls."""
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    worst = max(factors, key=float)
    out = {}
    for m in measurements:
        stages = [(f"{order}-{s['Sources']}", int(order[1:]), s["Sources"]) for order, s in m["Stages"].items()]
        local = m["LocalEdge"]
        estimate = estimate_stages.estimate(m["EntityCounts"], stages, model=model, profile=profile, block_size=block_size,
                                            local_edge=(local["Name"], local["Order"], local["Sources"]))
        per_stage = {}
        actual_total = 0.0
        for order, s in m["Stages"].items():
            est = estimate["Stages"][f"{order}-{s['Sources']}"]
            worker = est["ByPCGFactor"]["1.0"]["WorkerSecondsEstimate"]
            reducer = est["ReducerSecondsEstimate"]
            per_stage[order] = {"WorkerEstimateOverActual": worker / s["WorkerWallSeconds"],
                                "ReducerEstimateOverActual": reducer / s["ReducerWallSeconds"],
                                "StageEstimateOverActual": (worker + reducer) / s["StageWallSeconds"],
                                "ActualStageWallSeconds": s["StageWallSeconds"], "Estimate1xSeconds": worker + reducer}
            actual_total += s["StageWallSeconds"]
        est_local = estimate["Stages"][local["Name"]]["ByPCGFactor"]["1.0"]["StageSecondsEstimate"]
        per_stage["local-edge"] = {"StageEstimateOverActual": est_local / local["WallSeconds"],
                                   "ActualStageWallSeconds": local["WallSeconds"], "Estimate1xSeconds": est_local}
        actual_total += local["WallSeconds"]
        out[m["Case"]] = {"Stages": per_stage,
                          "ActualStageWallSecondsTotal": actual_total,
                          "Estimate1xOverActual": estimate["JobSecondsEstimateByPCGFactor"]["1.0"] / actual_total,
                          "EstimateWorstWithMarginOverActual": estimate["JobSecondsEstimateWithPreflightAndMargin"][worst] / actual_total,
                          "FitsOneJob": estimate["FitsOneJob"],
                          "JobSecondsEstimateWithPreflightAndMargin": estimate["JobSecondsEstimateWithPreflightAndMargin"],
                          "MaxPalacePeakGBEstimate": estimate["MaxPalacePeakGBEstimate"]}
    out["MinimumStageEstimateOverActual"] = min(stage["StageEstimateOverActual"]
                                                for case in out.values() for stage in case["Stages"].values())
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
    model["SelfCheck"] = {"Rule": ("estimate_stages on every coupon of the run under this model at the measured PCG counts: every "
                                   "stage's estimate over its measured wall is >= 1 (MinimumStageEstimateOverActual), the job at the "
                                   "worst PCG factor with preflight and margin over the sum of the measured stage walls is the "
                                   "conservatism ratio; PreviousModel is the same check under the previous model"),
                          "MinimumStageEstimateOverActual": check["MinimumStageEstimateOverActual"],
                          "Conservative": check["Conservative"],
                          "EstimateWorstWithMarginOverActual": {case: check[case]["EstimateWorstWithMarginOverActual"]
                                                                for case in check if case in model["Provenance"]["Coupons"]},
                          "Estimate1xOverActual": {case: check[case]["Estimate1xOverActual"]
                                                   for case in check if case in model["Provenance"]["Coupons"]},
                          "PreviousModel": {"MinimumStageEstimateOverActual": check_previous["MinimumStageEstimateOverActual"],
                                            "EstimateWorstWithMarginOverActual": {
                                                case: check_previous[case]["EstimateWorstWithMarginOverActual"]
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
