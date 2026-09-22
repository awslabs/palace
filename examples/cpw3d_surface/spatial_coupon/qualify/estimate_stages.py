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
fits when the 2x-PCG total with the preflight and margin factor is below the walltime
and the largest Palace peak plus the runner's headroom stays under NodeFitFraction of
the node.

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
from build_plan import DEFAULT_REDUCER_BLOCK_SIZE  # noqa: E402

COST_MODEL = HERE / "cost-model.json"
# The model every qualify run up to the decision-64a refit used (four-edge-physics-11 V-a,
# b28 executable, reducer block size 6): kept for the recorded estimates it produced;
# refit_cost_model.py binds it by digest under the current model's Previous.
PREVIOUS_COST_MODEL = HERE / "cost-model-physics11-b28.json"
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


def estimate_stage(model, order, sources, counts, block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    """One worker + reducer stage at `order` on `sources` sources, the reducer at
    PALACE_RESPONSE_BLOCK_SIZE `block_size`."""
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
    solve = measured["SecondsPerPCGIteration"] * measured["MeanPCGIterations"] * ratio
    other = (measured["MeanPerSourceSeconds"] - measured["SecondsPerPCGIteration"] * measured["MeanPCGIterations"]) * ratio
    reducer_setup = (measured["ReducerPalaceSeconds"] - measured["ReducerPairSeconds"]) * ratio
    reducer_evaluation = measured["ReducerPairSeconds"] * ratio * evaluation_fraction * evaluations / measured_evaluations
    reducer_gram = measured["ReducerPairSeconds"] * ratio * (1.0 - evaluation_fraction) * pairs / measured_pairs
    reducer = reducer_setup + reducer_evaluation + reducer_gram
    gb_per_gib = model["PalaceGBPerGiB"]
    h1 = h1_dofs_from_counts(counts, order)
    # 2b archived source fields are resident per block pair: the peak grows from the
    # measured block size by the measured per-field cost (decision 62(1)).
    resident_fields_gb = (max(2 * int(block_size) - 2 * model["MeasuredBlockSize"], 0)
                          * model["ReducerResidentFieldGBPerMillionH1"] * h1 / 1e6)
    stage = {"Order": order, "Sources": sources, "H1Estimate": h1,
             "DOFRatioVsMeasured": ratio, "ReducerPairs": pairs, "ReducerBlockSize": int(block_size),
             "ReducerBlockPairs": block_pairs(sources, block_size), "ReducerSourceEvaluations": evaluations,
             "ReducerSecondsEstimate": reducer,
             "ReducerSecondsEstimateParts": {"Setup": reducer_setup, "Evaluation": reducer_evaluation, "Gram": reducer_gram},
             "WorkerNonSourceSecondsEstimate": measured["WorkerNonSourceSeconds"] * ratio,
             "WorkerPalacePeakGBEstimate": measured["WorkerPalacePeakGB"] * ratio,
             "ReducerPalacePeakGBEstimate": measured["ReducerPalacePeakGB"] * ratio + resident_fields_gb,
             "ReducerResidentFieldsGBEstimate": resident_fields_gb,
             "ArchiveGBEstimate": measured["ArchiveGB"] * ratio * sources / measured["Sources"],
             "ByPCGFactor": {}}
    stage["NodeUsedGiBEstimateWorker"] = (measured["WorkerNodeUsedGiB"] - measured["WorkerPalacePeakGB"] / gb_per_gib
                                          + stage["WorkerPalacePeakGBEstimate"] / gb_per_gib)
    stage["NodeUsedGiBEstimateReducer"] = (measured["ReducerNodeUsedGiB"] - measured["ReducerPalacePeakGB"] / gb_per_gib
                                           + stage["ReducerPalacePeakGBEstimate"] / gb_per_gib)
    for factor in model["PCGFactors"]:
        worker = measured["WorkerNonSourceSeconds"] * ratio + sources * (other + solve * factor)
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
    for factor in model["PCGFactors"]:
        stage["ByPCGFactor"][f"{factor:.1f}"] = {
            "StageSecondsEstimate": (measured["NonSolveSeconds"] + measured["LinearSolveSeconds"] * factor * per_source) * ratio}
    return stage


def estimate(counts, stages, *, local_edge=None, model=None, profile=None, block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    """The estimate record: `stages` = [(name, order, sources), ...] worker + reducer
    stages at reducer block size `block_size`, `local_edge` = (name, order, sources) or None."""
    model = model or load_cost_model()
    profile = profile or json.loads(CLUSTER_PROFILE.read_text())
    node_gib, walltime = profile["NodeGiB"], profile["WalltimeSeconds"]
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    totals = {factor: 0.0 for factor in factors}
    peak_gb = 0.0
    out = {"Mesh": {"EntityCounts": counts, "VolumeElements": counts["Tetrahedra"] + counts["Prisms"] + counts["Pyramids"],
                    "H1ByOrder": {f"p{p}": h1_dofs_from_counts(counts, p) for p in (1, 2, 3, 4, 5)}},
           "CostModel": {"Path": model["Path"], "SHA256": model["SHA256"], "MeasuredMesh": model["MeasuredMesh"]["SHA256"],
                         "ClosedFormCheck": model.get("ClosedFormCheck"), "MeasuredBlockSize": model["MeasuredBlockSize"],
                         "ReducerEvaluationFraction": model["ReducerEvaluationFraction"],
                         "ReducerResidentFieldGBPerMillionH1": model["ReducerResidentFieldGBPerMillionH1"]},
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
    return out


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
