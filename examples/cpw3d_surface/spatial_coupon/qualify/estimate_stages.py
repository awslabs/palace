#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""DOF, memory and time estimate of a coupon's physics stages (four-edge-physics-11 /
gallery-physics-06b estimate_stages.py, parametrized by the cost model and the mesh).

H1 counts are Palace's exact closed form on the mesh's entity counts
(mixed_mesh.h1_dofs_from_counts; the counts come from library-build.json `H1.EntityCounts`
or are computed from the mesh).  Every stage is scaled from the cost model's measured
rates (per PCG iteration, per source, reducer setup + per pair, Palace peak, node used
GiB, archive GB) by the H1 ratio, with 1x / 1.5x / 2x the measured mean PCG counts;
the reducer splits into a setup part (~ the worker's non-source time) and a block-pair
part: the fraction ReducerEvaluationFraction of the measured pair seconds is the field
read + evaluation work (N x ceil(N / b) source evaluations at block size b, decision 62(1):
the measured rates ran at MeasuredBlockSize), the remainder the Gram work scaled with the
number of source pairs.  The decision is fail-closed: the job
fits when the worst-PCG total with the preflight and margin factor is below the walltime
and the largest Palace peak plus the runner's headroom stays under NodeFitFraction of
the node.  A Version-3 model (refit_cost_model.refit_measured) carries per-stage
SafetyFactor figures the times are multiplied by.

Multi-node stages (decision 457): plan_nodes attaches to every stage the minimum node
count whose estimated per-node node-used peak fits an instance's admission guard
(the model's NodeScaling calibration: a replicated fraction + the distributed remainder
/ N), 1 when the one-node rules of record admit the stage, None (fail closed) above the
profile's MaximumNodesPerJob; scale_to_nodes divides every time part by N^exponent (the
measured 1 -> 2 node speedups) for the node assignment the driver chose.

usage: estimate_stages.py --entity-counts JSON --stage ORDER:SOURCES ... [--local-edge ORDER:SOURCES]
       [--reducer-block-size B] [--cost-model PATH] [--cluster-profile PATH] [--out PATH]
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from mixed_mesh import h1_dofs_from_counts, h1_entity_counts  # noqa: E402
from build_plan import DEFAULT_REDUCER_BLOCK_SIZE, GIB, largest_node_gib  # noqa: E402

COST_MODEL = HERE / "cost-model.json"
# The model every qualify run up to the decision-64a refit used (four-edge-physics-11 V-a,
# b28 executable, reducer block size 6): kept for the recorded estimates it produced;
# refit_cost_model.py binds it by digest under the device model's Previous.
PREVIOUS_COST_MODEL = HERE / "cost-model-physics11-b28.json"
# The decision-64a device-library refit (2026-09-22, the largest-rate rule) every qualify run
# of stage 2 planned with (groups A / B / C, pair 5): kept for the plans of record it produced
# (the bitwise proof of decision 458 replays them under it) and bound by digest under the
# current model's Previous.
DEVICE_COST_MODEL = HERE / "cost-model-device-20260922.json"
CLUSTER_PROFILE = HERE / "cluster-profile.json"


def load_cost_model(path=COST_MODEL):
    model = json.loads(Path(path).read_text())
    counts = model["MeasuredMesh"]["EntityCounts"]
    check = {}
    for name, stage in model["Stages"].items():
        closed = h1_dofs_from_counts(counts, int(name[1:]))
        check[name] = {"ClosedForm": closed, "Measured": stage["H1"]}
        if closed != stage["H1"]:
            raise ValueError(f"cost model {path}: closed-form H1 {closed} at {name} differs from the measured "
                             f"{stage['H1']} - the model file is inconsistent")
    model["ClosedFormCheck"] = check
    model["Path"] = str(Path(path).resolve())
    model["SHA256"] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return model


def entity_counts_of_mesh(mesh_path):
    from mesh_array_io import read_mesh
    return h1_entity_counts(read_mesh(mesh_path))


def blocks_of(sources, block_size):
    return -(-sources // block_size)


def block_pairs(sources, block_size):
    blocks = blocks_of(sources, block_size)
    return blocks * (blocks + 1) // 2


def source_evaluations(sources, block_size):
    """Archived source fields the reducer reads and evaluates: every source once per block
    pair its block takes part in (ceil(N / b) block pairs per block)."""
    return sources * blocks_of(sources, block_size)


def streaming_reducer_peak_gb(model, h1, sources):
    """The reducer peak estimate (GB) under the streaming one-pass executable
    (cost-model ReducerPeakStreaming.Rule): SafetyFactor x PalaceGBPerGiB x the node-used
    line (baseline + (per-million-H1 + per-source-per-million-H1 x sources) x H1 / 1e6) -
    conservative by construction against every measured node-used peak of the calibration."""
    streaming = model["ReducerPeakStreaming"]
    node_used_gib = (streaming["NodeBaselineGiB"]
                     + (streaming["NodeUsedGiBPerMillionH1"] + streaming["NodeUsedGiBPerMillionH1PerSource"] * sources) * h1 / 1e6)
    return streaming["SafetyFactor"] * model["PalaceGBPerGiB"] * node_used_gib


def safety_factor(measured, kind):
    """The per-stage time safety factor of a Version-3 model (refit_cost_model.refit_measured:
    the largest measured wall / mean-rate estimate at the coupon's own PCG count over the
    measured runs; `kind` Worker / Reducer / LocalEdge); 1 for a model without one."""
    factors = measured.get("SafetyFactor")
    if factors is None:
        return 1.0
    return float(factors[kind]) if isinstance(factors, dict) else float(factors)


def estimate_stage(model, order, sources, counts, block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    """One worker + reducer stage at `order` on `sources` sources, the reducer at
    PALACE_RESPONSE_BLOCK_SIZE `block_size`; the times carry the model's per-stage safety
    factors (Version 3) when present."""
    measured = model["Stages"].get(f"p{order}")
    if measured is None:
        raise ValueError(f"the cost model has no measured stage at order {order} (orders {sorted(model['Stages'])})")
    if int(block_size) < 1:
        raise ValueError(f"the reducer block size must be an integer >= 1, not {block_size}")
    ratio = h1_dofs_from_counts(counts, order) / measured["H1"]
    pairs = sources * (sources + 1) // 2
    measured_pairs = measured["Sources"] * (measured["Sources"] + 1) // 2
    evaluations = source_evaluations(sources, block_size)
    measured_evaluations = source_evaluations(measured["Sources"], model["MeasuredBlockSize"])
    evaluation_fraction = model["ReducerEvaluationFraction"]
    worker_safety, reducer_safety = safety_factor(measured, "Worker"), safety_factor(measured, "Reducer")
    solve = measured["SecondsPerPCGIteration"] * measured["MeanPCGIterations"] * ratio * worker_safety
    other = (measured["MeanPerSourceSeconds"] - measured["SecondsPerPCGIteration"] * measured["MeanPCGIterations"]) * ratio * worker_safety
    reducer_setup = (measured["ReducerPalaceSeconds"] - measured["ReducerPairSeconds"]) * ratio * reducer_safety
    reducer_evaluation = measured["ReducerPairSeconds"] * ratio * evaluation_fraction * evaluations / measured_evaluations * reducer_safety
    reducer_gram = measured["ReducerPairSeconds"] * ratio * (1.0 - evaluation_fraction) * pairs / measured_pairs * reducer_safety
    reducer = reducer_setup + reducer_evaluation + reducer_gram
    gb_per_gib = model["PalaceGBPerGiB"]
    h1 = h1_dofs_from_counts(counts, order)
    # 2b archived source fields are resident per block pair: the peak grows from the
    # measured block size by the measured per-field cost (decision 62(1)).
    resident_fields_gb = (max(2 * int(block_size) - 2 * model["MeasuredBlockSize"], 0)
                          * model["ReducerResidentFieldGBPerMillionH1"] * h1 / 1e6)
    previous_reducer_peak = measured["ReducerPalacePeakGB"] * ratio + resident_fields_gb
    streaming = model.get("ReducerPeakStreaming")
    reducer_peak = streaming_reducer_peak_gb(model, h1, sources) if streaming else previous_reducer_peak
    stage = {"Order": order, "Sources": sources, "H1Estimate": h1,
             "DOFRatioVsMeasured": ratio, "ReducerPairs": pairs, "ReducerBlockSize": int(block_size),
             "ReducerBlockPairs": block_pairs(sources, block_size), "ReducerSourceEvaluations": evaluations,
             "ReducerSecondsEstimate": reducer,
             "ReducerSecondsEstimateParts": {"Setup": reducer_setup, "Evaluation": reducer_evaluation, "Gram": reducer_gram},
             "SafetyFactor": {"Worker": worker_safety, "Reducer": reducer_safety},
             "WorkerNonSourceSecondsEstimate": measured["WorkerNonSourceSeconds"] * ratio * worker_safety,
             "WorkerPalacePeakGBEstimate": measured["WorkerPalacePeakGB"] * ratio,
             "ReducerPalacePeakGBEstimate": reducer_peak,
             "ReducerPalacePeakGBEstimatePrevious": previous_reducer_peak,
             "ReducerPalacePeakGBRule": ("ReducerPeakStreaming (a node-used figure at the streaming executable; the Previous term "
                                         "recorded)" if streaming else "Previous: measured peak x H1 ratio + resident fields"),
             "ReducerResidentFieldsGBEstimate": resident_fields_gb,
             "ArchiveGBEstimate": measured["ArchiveGB"] * ratio * sources / measured["Sources"],
             "ByPCGFactor": {}}
    stage["NodeUsedGiBEstimateWorker"] = (measured["WorkerNodeUsedGiB"] - measured["WorkerPalacePeakGB"] / gb_per_gib
                                          + stage["WorkerPalacePeakGBEstimate"] / gb_per_gib)
    # The streaming estimate already is a conservative node-used figure (baseline included).
    stage["NodeUsedGiBEstimateReducer"] = (reducer_peak / gb_per_gib if streaming else
                                           measured["ReducerNodeUsedGiB"] - measured["ReducerPalacePeakGB"] / gb_per_gib
                                           + stage["ReducerPalacePeakGBEstimate"] / gb_per_gib)
    for factor in model["PCGFactors"]:
        worker = stage["WorkerNonSourceSecondsEstimate"] + sources * (other + solve * factor)
        stage["ByPCGFactor"][f"{factor:.1f}"] = {"MeanPCGIterations": measured["MeanPCGIterations"] * factor,
                                                 "PerSourceSecondsEstimate": other + solve * factor,
                                                 "WorkerSecondsEstimate": worker,
                                                 "StageSecondsEstimate": worker + reducer}
    return stage


def estimate_local_edge(model, order, sources, counts):
    measured = model["LocalEdge"]
    if int(measured["Order"]) != int(order):
        raise ValueError(f"the cost model measured the local-edge stage at order {measured['Order']}, not {order}")
    ratio = h1_dofs_from_counts(counts, order) / measured["H1"]
    gb_per_gib = model["PalaceGBPerGiB"]
    stage = {"Order": order, "Sources": sources, "H1Estimate": h1_dofs_from_counts(counts, order),
             "DOFRatioVsMeasured": ratio, "PalacePeakGBEstimate": measured["PalacePeakGB"] * ratio,
             "NodeUsedGiBEstimate": (measured["NodeUsedGiB"] - measured["PalacePeakGB"] / gb_per_gib
                                     + measured["PalacePeakGB"] * ratio / gb_per_gib),
             "ByPCGFactor": {}}
    # The measured solve time scales with the number of controls; the setup / estimator
    # part is taken as fixed (identical to the measurement at its control count).
    per_source = sources / measured["Sources"]
    safety = safety_factor(measured, "LocalEdge")
    stage["SafetyFactor"] = safety
    for factor in model["PCGFactors"]:
        stage["ByPCGFactor"][f"{factor:.1f}"] = {
            "StageSecondsEstimate": (measured["NonSolveSeconds"] + measured["LinearSolveSeconds"] * factor * per_source) * ratio * safety}
    return stage


def estimate(counts, stages, *, local_edge=None, model=None, profile=None, block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    """The estimate record: `stages` = [(name, order, sources), ...] worker + reducer
    stages at reducer block size `block_size`, `local_edge` = (name, order, sources) or None."""
    model = model or load_cost_model()
    profile = profile or json.loads(CLUSTER_PROFILE.read_text())
    node_gib, walltime = largest_node_gib(profile), profile["WalltimeSeconds"]
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    totals = {factor: 0.0 for factor in factors}
    peak_gb = 0.0
    out = {"Mesh": {"EntityCounts": counts, "VolumeElements": counts["Tetrahedra"] + counts["Prisms"] + counts["Pyramids"],
                    "H1ByOrder": {f"p{p}": h1_dofs_from_counts(counts, p) for p in (1, 2, 3, 4, 5)}},
           "CostModel": {"Path": model["Path"], "SHA256": model["SHA256"], "MeasuredMesh": model["MeasuredMesh"]["SHA256"],
                         "ClosedFormCheck": model.get("ClosedFormCheck"), "MeasuredBlockSize": model["MeasuredBlockSize"],
                         "ReducerEvaluationFraction": model["ReducerEvaluationFraction"],
                         "ReducerResidentFieldGBPerMillionH1": model["ReducerResidentFieldGBPerMillionH1"],
                         "PalaceGBPerGiB": model["PalaceGBPerGiB"]},
           "ReducerBlockSize": int(block_size),
           "Stages": {}}
    for name, order, sources in stages:
        stage = estimate_stage(model, order, sources, counts, block_size)
        stage["FitsNode"] = max(stage["NodeUsedGiBEstimateWorker"], stage["NodeUsedGiBEstimateReducer"]) < model["NodeFitFraction"] * node_gib
        for factor in factors:
            totals[factor] += stage["ByPCGFactor"][factor]["StageSecondsEstimate"]
        peak_gb = max(peak_gb, stage["WorkerPalacePeakGBEstimate"], stage["ReducerPalacePeakGBEstimate"])
        out["Stages"][name] = stage
    if local_edge is not None:
        name, order, sources = local_edge
        stage = estimate_local_edge(model, order, sources, counts)
        stage["FitsNode"] = stage["NodeUsedGiBEstimate"] < model["NodeFitFraction"] * node_gib
        for factor in factors:
            totals[factor] += stage["ByPCGFactor"][factor]["StageSecondsEstimate"]
        peak_gb = max(peak_gb, stage["PalacePeakGBEstimate"])
        out["Stages"][name] = stage
    out["JobSecondsEstimateByPCGFactor"] = dict(totals)
    out["JobSecondsEstimateWithPreflightAndMargin"] = {
        factor: total * model["PreflightAndMarginFactor"] + model["PreflightSeconds"] for factor, total in totals.items()}
    out["MaxPalacePeakGBEstimate"] = peak_gb
    out["NodeGiB"] = node_gib
    out["WalltimeSeconds"] = walltime
    worst = max(factors, key=float)
    job = out["JobSecondsEstimateWithPreflightAndMargin"][worst]
    fits_time = job < walltime
    fits_memory = peak_gb / model["PalaceGBPerGiB"] + 60 < model["NodeFitFraction"] * node_gib
    out["FitsOneJob"] = bool(fits_time and fits_memory)
    out["Decision"] = (f"every stage fits ONE job of {walltime / 3600:.0f} h: runner total {totals[factors[0]] / 60:.0f} min at the "
                       f"measured PCG counts, {totals[worst] / 60:.0f} min at {worst}x (+{100 * (model['PreflightAndMarginFactor'] - 1):.0f}% "
                       f"and preflight: {job / 60:.0f} min); largest Palace peak {peak_gb:.0f} GB of {node_gib:.0f} GiB"
                       if out["FitsOneJob"] else
                       f"does NOT fit one job of {walltime / 3600:.0f} h at {worst}x PCG ({job / 60:.0f} min; largest Palace peak "
                       f"{peak_gb:.0f} GB of {node_gib:.0f} GiB): fail closed, not submitted")
    plan_nodes(out, model, profile)
    return out


# ---------------------------------------------------------------------------------------
# Multi-node stages (decision 457): the node count of a stage is the minimum N whose
# estimated per-node node-used peak fits an instance's admission guard; a stage that fits
# one node under the one-node rules of record keeps them (its figures are untouched).
# ---------------------------------------------------------------------------------------
NODE_SCALING_RULE = ("cost-model NodeScaling (measured at 1 and 2 nodes on one stored coupon, decision 457): per-node node-used "
                     "GiB at N nodes = the one-node figure x (ReplicatedFraction + (1 - ReplicatedFraction) / N) per stage kind "
                     "(Worker / Reducer / LocalEdge); every time part scales as T(1) / N^Exponent with the exponent = log2 of the "
                     "measured 1 -> 2 node speedup of that part (WorkerPerSource, WorkerNonSource, ReducerSetup, ReducerReduction, "
                     "LocalEdge), never above 1; a stage fits N nodes when its largest per-node figure is <= MemoryFitFraction x "
                     "MemoryGiB of an instance x (1 - PerNodeGuardMargin) (decision 463; the admission guard MinimumMemAvailableBytes "
                     "is checked on every node); NodesRequired is "
                     "the minimum such N in 2..MaximumNodesPerJob (1 when the one-node rules of record admit the stage); a stage "
                     "above MaximumNodesPerJob fails closed")
ONE_NODE_HEADROOM_GIB = 60
MEMORY_KINDS = {"Worker": "Worker", "Reducer": "Reducer", "LocalEdge": "LocalEdge"}
# Decision 463: the estimated per-node peak must stay under the admission guard by a margin
# (the measured p5 dense twin at 2 nodes peaked at 460.6 GiB against the 460.8 GiB m8g guard);
# the model's NodeScaling.PerNodeGuardMargin (from the measured residuals, >= 0.10) overrides.
DEFAULT_PER_NODE_GUARD_MARGIN = 0.10


def per_node_guard_margin(model):
    scaling = (model or {}).get("NodeScaling") or {}
    return float(scaling.get("PerNodeGuardMargin", DEFAULT_PER_NODE_GUARD_MARGIN))


def one_node_fits(model, profile, node_used_gib, palace_peak_gb):
    """The one-node rules of record for one stage: the node-used estimate under NodeFitFraction
    of the largest node (estimate's FitsNode) and the Palace peak plus the runner's headroom
    under the same bound (estimate's FitsOneJob memory term)."""
    bound = model["NodeFitFraction"] * largest_node_gib(profile)
    return bool(node_used_gib < bound and palace_peak_gb / model["PalaceGBPerGiB"] + ONE_NODE_HEADROOM_GIB < bound)


def replicated_fraction(model, kind):
    scaling = model.get("NodeScaling")
    if scaling is None:
        raise ValueError("the cost model carries no NodeScaling calibration: a stage that does not fit one node cannot be planned")
    return float(scaling["Memory"][MEMORY_KINDS[kind]]["ReplicatedFraction"])


def per_node_used_gib(model, kind, node_used_one, nodes):
    """The per-node node-used GiB of a `kind` (Worker / Reducer / LocalEdge) stage at `nodes`
    nodes from its one-node figure (NODE_SCALING_RULE)."""
    nodes = int(nodes)
    if nodes < 1:
        raise ValueError(f"the node count must be >= 1, not {nodes}")
    if nodes == 1:
        return float(node_used_one)
    fraction = replicated_fraction(model, kind)
    return float(node_used_one) * (fraction + (1.0 - fraction) / nodes)


def time_exponent(model, part):
    scaling = model.get("NodeScaling")
    if scaling is None:
        raise ValueError("the cost model carries no NodeScaling calibration: a stage that does not fit one node cannot be planned")
    return min(float(scaling["Time"][part]["Exponent"]), 1.0)


def scaled_seconds(model, part, seconds_one, nodes):
    return float(seconds_one) if int(nodes) == 1 else float(seconds_one) / float(nodes) ** time_exponent(model, part)


def select_instance_for_nodes(profile, per_node_gib, nodes, margin=DEFAULT_PER_NODE_GUARD_MARGIN):
    """The first instance of Instances whose admission guard (MemoryFitFraction x MemoryGiB)
    holds `per_node_gib` (the largest estimated per-node node-used GiB of the stages a job
    runs at `nodes` nodes) with the guard margin (decision 463: per node <= guard x (1 - margin));
    None when no instance holds it."""
    fraction = float(profile["MemoryFitFraction"])
    for item in profile["Instances"]:
        guard_gib = fraction * item["MemoryGiB"]
        if per_node_gib <= guard_gib * (1.0 - margin):
            return {"Type": item["Type"], "vCPUs": item["vCPUs"], "MemoryGiB": item["MemoryGiB"], "NodeGiB": item["NodeGiB"],
                    "Nodes": int(nodes), "PerNodeUsedGiBEstimate": float(per_node_gib), "MemoryFitFraction": fraction, "Fits": True,
                    "PerNodeGuardMargin": margin, "PerNodeGuardGiB": guard_gib * (1.0 - margin),
                    "MinimumMemAvailableBytes": int(fraction * item["MemoryGiB"] * GIB),
                    "Candidates": [entry["Type"] for entry in profile["Instances"]], "Rule": profile["MultiNodeRule"]}
    return None


def stage_memory_kinds(stage):
    """(kind, node-used GiB, Palace peak GB) of every executable of a stage record."""
    if "NodeUsedGiBEstimateWorker" in stage:
        return [("Worker", stage["NodeUsedGiBEstimateWorker"], stage["WorkerPalacePeakGBEstimate"]),
                ("Reducer", stage["NodeUsedGiBEstimateReducer"], stage["ReducerPalacePeakGBEstimate"])]
    return [("LocalEdge", stage["NodeUsedGiBEstimate"], stage["PalacePeakGBEstimate"])]


def nodes_required(model, profile, stage):
    """The node plan of one stage record: {NodesRequired, OneNodeFits, Instance (at NodesRequired),
    PerNodeUsedGiB (by N), MaximumNodesPerJob}; NodesRequired None when even MaximumNodesPerJob
    nodes do not hold the stage (or the model has no NodeScaling calibration: Reason)."""
    kinds = stage_memory_kinds(stage)
    maximum = int(profile.get("MaximumNodesPerJob", 1))
    record = {"OneNodeFits": all(one_node_fits(model, profile, used, peak) for _, used, peak in kinds),
              "MaximumNodesPerJob": maximum, "PerNodeUsedGiB": {}, "Rule": NODE_SCALING_RULE}
    if record["OneNodeFits"]:
        record.update(NodesRequired=1, Instance=None)
        return record
    try:
        for nodes in range(2, maximum + 1):
            per_node = {kind: per_node_used_gib(model, kind, used, nodes) for kind, used, _ in kinds}
            record["PerNodeUsedGiB"][str(nodes)] = per_node
            instance = select_instance_for_nodes(profile, max(per_node.values()), nodes, per_node_guard_margin(model))
            if instance is not None:
                record.update(NodesRequired=nodes, Instance=instance)
                return record
    except ValueError as error:
        record.update(NodesRequired=None, Instance=None, Reason=str(error))
        return record
    largest = max((per_node_used_gib(model, kind, used, maximum) for kind, used, _ in kinds), default=0.0)
    margin = per_node_guard_margin(model)
    record.update(NodesRequired=None, Instance=None,
                  Reason=(f"the per-node node-used estimate at MaximumNodesPerJob = {maximum} nodes ({largest:.0f} GiB) exceeds every "
                          f"instance's admission guard less the {100 * margin:.0f} % margin "
                          f"({[round(profile['MemoryFitFraction'] * item['MemoryGiB'] * (1 - margin)) for item in profile['Instances']]} GiB)"))
    return record


def plan_nodes(estimate, model, profile):
    """Attach the node plan to an estimate record: every stage's NodePlan (nodes_required) and
    the summary Nodes = {MaximumNodesPerJob, Required: {stage: N}, MultiNode: any N > 1, Fits:
    every stage holds <= MaximumNodesPerJob nodes, Decision}; a coupon whose every stage fits one
    node is reported as such and its figures are untouched."""
    required = {}
    for name, stage in estimate["Stages"].items():
        stage["NodePlan"] = nodes_required(model, profile, stage)
        required[name] = stage["NodePlan"]["NodesRequired"]
    fits = all(nodes is not None for nodes in required.values())
    multi = fits and any(nodes > 1 for nodes in required.values())
    if not fits:
        failed = {name: estimate["Stages"][name]["NodePlan"].get("Reason") for name, nodes in required.items() if nodes is None}
        decision = (f"does NOT fit {profile.get('MaximumNodesPerJob', 1)} nodes (MaximumNodesPerJob): {failed}: fail closed, not submitted")
    elif multi:
        decision = ("multi-node stages (decision 457): nodes required per stage " + ", ".join(
            f"{name} {nodes}" + (f" ({estimate['Stages'][name]['NodePlan']['Instance']['Type']}, per node "
                                 f"{estimate['Stages'][name]['NodePlan']['Instance']['PerNodeUsedGiBEstimate']:.0f} GiB)" if nodes > 1 else "")
            for name, nodes in required.items()))
    else:
        decision = "every stage fits one node under the one-node rules of record"
    estimate["Nodes"] = {"MaximumNodesPerJob": int(profile.get("MaximumNodesPerJob", 1)), "Required": required, "MultiNode": multi,
                         "Fits": fits, "Decision": decision, "PerNodeGuardMargin": per_node_guard_margin(model), "Rule": NODE_SCALING_RULE}
    return estimate["Nodes"]


def scale_stage_to_nodes(model, stage, nodes):
    """A copy of a stage record at `nodes` nodes: every time part divided by N^Exponent of its
    NodeScaling part, the per-node memory figures added (PerNodeUsedGiBEstimate*,
    PerNodePalacePeakGB*), the one-node figures kept under OneNode; `nodes` = 1 returns a copy."""
    nodes = int(nodes)
    scaled = json.loads(json.dumps(stage))
    scaled["Nodes"] = nodes
    if nodes == 1:
        return scaled
    scaled["OneNode"] = {key: stage[key] for key in stage if key not in ("NodePlan",)}
    if "NodeUsedGiBEstimateWorker" in stage:
        parts = stage["ReducerSecondsEstimateParts"]
        setup = scaled_seconds(model, "ReducerSetup", parts["Setup"], nodes)
        evaluation = scaled_seconds(model, "ReducerReduction", parts["Evaluation"], nodes)
        gram = scaled_seconds(model, "ReducerReduction", parts["Gram"], nodes)
        scaled["ReducerSecondsEstimateParts"] = {"Setup": setup, "Evaluation": evaluation, "Gram": gram}
        scaled["ReducerSecondsEstimate"] = setup + evaluation + gram
        scaled["WorkerNonSourceSecondsEstimate"] = scaled_seconds(model, "WorkerNonSource", stage["WorkerNonSourceSecondsEstimate"], nodes)
        sources = stage["Sources"]
        for factor, figures in stage["ByPCGFactor"].items():
            per_source = scaled_seconds(model, "WorkerPerSource", figures["PerSourceSecondsEstimate"], nodes)
            worker = scaled["WorkerNonSourceSecondsEstimate"] + sources * per_source
            scaled["ByPCGFactor"][factor] = {"MeanPCGIterations": figures["MeanPCGIterations"], "PerSourceSecondsEstimate": per_source,
                                             "WorkerSecondsEstimate": worker, "StageSecondsEstimate": worker + scaled["ReducerSecondsEstimate"]}
        for kind, suffix in (("Worker", "Worker"), ("Reducer", "Reducer")):
            factor = replicated_fraction(model, kind) + (1.0 - replicated_fraction(model, kind)) / nodes
            scaled[f"PerNodeUsedGiBEstimate{suffix}"] = per_node_used_gib(model, kind, stage[f"NodeUsedGiBEstimate{suffix}"], nodes)
            scaled[f"PerNodePalacePeakGBEstimate{suffix}"] = stage[f"{suffix}PalacePeakGBEstimate"] * factor
    else:
        for factor, figures in stage["ByPCGFactor"].items():
            scaled["ByPCGFactor"][factor] = {"StageSecondsEstimate": scaled_seconds(model, "LocalEdge", figures["StageSecondsEstimate"], nodes)}
        factor = replicated_fraction(model, "LocalEdge") + (1.0 - replicated_fraction(model, "LocalEdge")) / nodes
        scaled["PerNodeUsedGiBEstimate"] = per_node_used_gib(model, "LocalEdge", stage["NodeUsedGiBEstimate"], nodes)
        scaled["PerNodePalacePeakGBEstimate"] = stage["PalacePeakGBEstimate"] * factor
    return scaled


def stage_per_node_used_gib(stage):
    """The largest per-node node-used GiB a (scaled) stage record carries (the one-node figure
    of an unscaled stage)."""
    if "NodeUsedGiBEstimateWorker" in stage:
        return max(stage.get("PerNodeUsedGiBEstimateWorker", stage["NodeUsedGiBEstimateWorker"]),
                   stage.get("PerNodeUsedGiBEstimateReducer", stage["NodeUsedGiBEstimateReducer"]))
    return stage.get("PerNodeUsedGiBEstimate", stage["NodeUsedGiBEstimate"])


def scale_to_nodes(estimate, assignment, model, profile):
    """The estimate record of a multi-node coupon: `assignment` = stage name -> nodes (the
    main stages share one count, the control / local-edge group another: build_plan /
    job_split read them per job); every stage scaled by scale_stage_to_nodes, the job totals
    recomputed, FitsOneJob = the time term alone (the memory term is the per-node fit the
    assignment was chosen for), the Decision stating the nodes per stage."""
    scaled = json.loads(json.dumps(estimate))
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    totals = {factor: 0.0 for factor in factors}
    walltime = profile["WalltimeSeconds"]
    for name, stage in estimate["Stages"].items():
        nodes = int(assignment[name])
        if stage["NodePlan"]["NodesRequired"] is None or nodes < stage["NodePlan"]["NodesRequired"]:
            raise ValueError(f"stage {name} needs {stage['NodePlan']['NodesRequired']} nodes, assigned {nodes}")
        if nodes > int(profile["MaximumNodesPerJob"]):
            raise ValueError(f"stage {name} assigned {nodes} nodes above MaximumNodesPerJob {profile['MaximumNodesPerJob']}")
        record = scale_stage_to_nodes(model, stage, nodes)
        instance = (select_instance_for_nodes(profile, stage_per_node_used_gib(record), nodes, per_node_guard_margin(model))
                    if nodes > 1 else None)
        if nodes > 1 and instance is None:
            raise ValueError(f"stage {name} at {nodes} nodes fits no instance's admission guard")
        record["Instance"] = instance
        for factor in factors:
            totals[factor] += record["ByPCGFactor"][factor]["StageSecondsEstimate"]
        scaled["Stages"][name] = record
    scaled["JobSecondsEstimateByPCGFactor"] = dict(totals)
    scaled["JobSecondsEstimateWithPreflightAndMargin"] = {
        factor: total * model["PreflightAndMarginFactor"] + model["PreflightSeconds"] for factor, total in totals.items()}
    worst = max(factors, key=float)
    job = scaled["JobSecondsEstimateWithPreflightAndMargin"][worst]
    scaled["FitsOneJob"] = bool(job < walltime)
    scaled["Nodes"]["Assigned"] = {name: int(nodes) for name, nodes in assignment.items()}
    scaled["Nodes"]["NodeSecondsEstimateByPCGFactor"] = {
        factor: sum(scaled["Stages"][name]["ByPCGFactor"][factor]["StageSecondsEstimate"] * int(assignment[name]) for name in scaled["Stages"])
        for factor in factors}
    per_stage = ", ".join(f"{name} on {int(assignment[name])} node(s)"
                          + (f" ({scaled['Stages'][name]['Instance']['Type']}, per node "
                             f"{scaled['Stages'][name]['Instance']['PerNodeUsedGiBEstimate']:.0f} GiB)" if int(assignment[name]) > 1 else "")
                          for name in scaled["Stages"])
    scaled["Decision"] = (f"multi-node coupon (decision 457): {per_stage}; runner total {totals[factors[0]] / 60:.0f} min at the measured "
                          f"PCG counts, {totals[worst] / 60:.0f} min at {worst}x (+{100 * (model['PreflightAndMarginFactor'] - 1):.0f}% and "
                          f"preflight: {job / 60:.0f} min) in one job of {walltime / 3600:.0f} h: "
                          + ("fits" if scaled["FitsOneJob"] else "does NOT fit one job (the job split decides)"))
    return scaled


def parse_stage(text):
    order, sources = text.split(":")
    return int(order), int(sources)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--entity-counts", type=Path, help="JSON with the H1 entity counts (library-build.json H1.EntityCounts)")
    parser.add_argument("--mesh", type=Path, help="mesh to count the entities of (when no --entity-counts)")
    parser.add_argument("--stage", action="append", default=[], metavar="ORDER:SOURCES", help="worker + reducer stage (repeatable)")
    parser.add_argument("--local-edge", metavar="ORDER:SOURCES", help="the ordinary-path local-edge stage")
    parser.add_argument("--reducer-block-size", type=int, default=DEFAULT_REDUCER_BLOCK_SIZE,
                        help=f"PALACE_RESPONSE_BLOCK_SIZE of the reducer stages (default {DEFAULT_REDUCER_BLOCK_SIZE})")
    parser.add_argument("--cost-model", type=Path, default=COST_MODEL)
    parser.add_argument("--cluster-profile", type=Path, default=CLUSTER_PROFILE)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args(argv)
    if (args.entity_counts is None) == (args.mesh is None):
        parser.error("exactly one of --entity-counts / --mesh")
    counts = json.loads(args.entity_counts.read_text()) if args.entity_counts else entity_counts_of_mesh(args.mesh)
    stages = [(f"p{order}-{sources}", order, sources) for order, sources in map(parse_stage, args.stage)]
    local = None
    if args.local_edge:
        order, sources = parse_stage(args.local_edge)
        local = (f"local-edge-p{order}-{sources}", order, sources)
    record = estimate(counts, stages, local_edge=local, model=load_cost_model(args.cost_model),
                      profile=json.loads(args.cluster_profile.read_text()), block_size=args.reducer_block_size)
    text = json.dumps(record, indent=2) + "\n"
    if args.out:
        args.out.write_text(text)
    print(record["Decision"])
    return 0 if record["FitsOneJob"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
