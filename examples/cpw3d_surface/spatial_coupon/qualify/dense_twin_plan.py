#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The (F) dense twin as a (multi-node) Palace job (decision 457 (4)): the dense held-out
traces of a coupon solved on its identity mesh at the library and control orders
(`coupon_library.py spatial-qualify traces` writes the configs: every trace a
PrescribedPotential excitation of ONE ordinary Palace run per twin and order, the response
matrix off), planned as a run_stages plan - one ordinary stage per run, configs / mesh /
traces pinned, caps from the dense-twin cost model - on the minimum node count whose
per-node peak fits an instance's admission guard (estimate_stages NodeScaling, the ordinary
stage's LocalEdge kind), with the job script of build_plan (select = N, hostfile mpirun).

The dense-twin cost model (dense-twin-model.json, `calibrate`): per order the LARGEST
measured Palace peak per million H1 and the solve seconds per million H1 per trace / the
non-solve seconds per million H1 of the measured dense-twin runs (the pair-5 fab / thin twins
at p4 / p5, PBS 57628 / 57629), scaled to another coupon by the exact closed-form H1.

usage: dense_twin_plan.py plan --dense-dir DIR --run NAME ... --entity-counts JSON --binary-sha256 SHA
           --remote-root ROOT --job-dir REMOTE_DIR --out LOCAL_DIR [--nodes N] [--cost-model PATH]
           [--dense-twin-model PATH] [--cluster-profile PATH] [--walltime-seconds S]
       dense_twin_plan.py calibrate --log NAME=palace.log ... --executable SHA --pbs-jobs IDS --out dense-twin-model.json
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import build_plan  # noqa: E402
import estimate_stages  # noqa: E402
from mixed_mesh import h1_dofs_from_counts  # noqa: E402
from refit_cost_model import linear_solve_seconds, memory_gb, parse_elapsed_time_report, round_up  # noqa: E402
from run_stages import parse_log  # noqa: E402

DENSE_TWIN_MODEL = HERE / "dense-twin-model.json"
STAGE_RULE = ("an ordinary Palace run (no PALACE_RESPONSE_* variable: the per-source postprocessing of surface-Q*.csv / domain-E.csv) "
              "of the dense traces at one order; memory and time scale from the dense-twin model by the exact closed-form H1 of the "
              "mesh at that order; the node count is the LocalEdge kind of the cost model's NodeScaling (the same ordinary stage)")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def calibrate(logs, *, executable, pbs_jobs):
    """The dense-twin model from measured runs: `logs` = {name: palace.log path}; per order the
    largest Palace GB per million H1, solve seconds per million H1 per trace and non-solve
    seconds per million H1 (fail-closed: the largest over the runs of the order)."""
    runs, by_order = {}, {}
    for name, path in logs.items():
        parsed = parse_log(Path(path).read_text(errors="replace"))
        timers = parse_elapsed_time_report(parsed["ElapsedTimeReport"])
        solve = linear_solve_seconds(timers)
        traces = len(parsed["Iterations"])
        h1 = parsed["H1"]
        run = {"Order": parsed["Order"], "H1": h1, "Traces": traces, "PCG": parsed["PCG"], "PalacePeakGB": memory_gb(parsed["PalacePeakMemory"]["Total"]),
               "PalaceTotalSeconds": parsed["PalaceTotalSeconds"], "LinearSolveSeconds": solve,
               "PalaceGBPerMillionH1": memory_gb(parsed["PalacePeakMemory"]["Total"]) / (h1 / 1e6),
               "SolveSecondsPerMillionH1PerTrace": solve / traces / (h1 / 1e6),
               "NonSolveSecondsPerMillionH1": (parsed["PalaceTotalSeconds"] - solve) / (h1 / 1e6),
               "MeanPCGIterations": sum(parsed["PCG"]) / len(parsed["PCG"]) if parsed["PCG"] else None}
        runs[name] = run
        by_order.setdefault(f"p{parsed['Order']}", []).append((name, run))
    orders = {}
    for order, items in sorted(by_order.items()):
        largest = {key: max(items, key=lambda item: item[1][key]) for key in
                   ("PalaceGBPerMillionH1", "SolveSecondsPerMillionH1PerTrace", "NonSolveSecondsPerMillionH1", "MeanPCGIterations")}
        orders[order] = {key: round_up(largest[key][1][key], 5) for key in largest}
        orders[order]["SetBy"] = {key: largest[key][0] for key in largest}
        orders[order]["Runs"] = [name for name, _ in items]
    return {"Version": 1, "Copyright": "Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.",
            "SPDX-License-Identifier": "Apache-2.0",
            "Purpose": ("Measured rates of the (F) dense-twin runs (one ordinary Palace run of the dense held-out traces per twin and "
                        f"order; executable {executable}, one node x 192 ranks, PBS {pbs_jobs}): per order the largest Palace peak per "
                        "million H1, solve seconds per million H1 per trace and non-solve seconds per million H1 over the measured runs; "
                        "dense_twin_plan.py scales them to another coupon mesh by the exact closed-form H1 at that order (decision 457 (4))"),
            "Executable": executable, "PBSJobs": pbs_jobs, "Orders": orders, "Runs": runs, "Rule": STAGE_RULE}


def estimate_run(model, dense_model, counts, order, traces, *, profile):
    """One dense-twin run's estimate: Palace peak, node-used GiB (the cost model's LocalEdge
    node overhead added), seconds at the cost model's PCG factors, the one-node fit and the
    nodes required (estimate_stages.nodes_required on a LocalEdge-shaped stage record)."""
    rates = dense_model["Orders"].get(f"p{order}")
    if rates is None:
        raise ValueError(f"the dense-twin model has no measured order p{order} (orders {sorted(dense_model['Orders'])})")
    h1 = h1_dofs_from_counts(counts, order)
    million = h1 / 1e6
    gb_per_gib = model["PalaceGBPerGiB"]
    local = model["LocalEdge"]
    peak_gb = rates["PalaceGBPerMillionH1"] * million
    overhead_gib = local["NodeUsedGiB"] - local["PalacePeakGB"] / gb_per_gib
    stage = {"Order": order, "Sources": traces, "H1Estimate": h1, "PalacePeakGBEstimate": peak_gb,
             "NodeUsedGiBEstimate": overhead_gib + peak_gb / gb_per_gib, "ByPCGFactor": {}, "Rule": STAGE_RULE}
    for factor in model["PCGFactors"]:
        stage["ByPCGFactor"][f"{factor:.1f}"] = {"StageSecondsEstimate": (rates["NonSolveSecondsPerMillionH1"]
                                                                          + rates["SolveSecondsPerMillionH1PerTrace"] * traces * factor) * million}
    stage["NodePlan"] = estimate_stages.nodes_required(model, profile, stage)
    return stage


def plan_dense_twins(*, dense_dir, runs, counts, binary_sha256, remote_root, job_dir, model, dense_model, profile,
                     nodes=None, walltime_seconds=None, case_id=None):
    """The run_stages plan and job script of the dense-twin runs `runs` of `dense_dir` (its
    dense-traces.json names the configs): one ordinary stage per run in order, the node
    count = the largest NodesRequired over the runs (or `nodes`), the instance by the largest
    per-node figure, every config / mesh / trace pinned."""
    dense_dir = Path(dense_dir)
    traces_record = json.loads((dense_dir / "dense-traces.json").read_text())
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    worst = max(factors, key=float)
    deadline = profile["DeadlineSeconds"] - profile["DeadlineMarginSeconds"]
    stages, estimates, pinned = [], {}, {}
    required = 1
    for name in runs:
        config_path = Path(traces_record["Configs"][name])
        config = json.loads(config_path.read_text())
        order = int(config["Solver"]["Order"])
        traces = len(config["Boundaries"]["PrescribedPotential"])
        estimate = estimate_run(model, dense_model, counts, order, traces, profile=profile)
        if estimate["NodePlan"]["NodesRequired"] is None:
            raise ValueError(f"{name}: {estimate['NodePlan'].get('Reason')}: fail closed")
        required = max(required, estimate["NodePlan"]["NodesRequired"])
        estimates[name] = estimate
        pinned[str(config_path)] = sha256(config_path)
        pinned.setdefault(config["Model"]["Mesh"], sha256(config["Model"]["Mesh"]) if Path(config["Model"]["Mesh"]).exists() else None)
        for trace in config["Boundaries"]["PrescribedPotential"]:
            if trace.get("DataFile"):
                pinned.setdefault(trace["DataFile"], sha256(trace["DataFile"]) if Path(trace["DataFile"]).exists() else None)
    job_nodes = int(nodes) if nodes else required
    if job_nodes < required:
        raise ValueError(f"--nodes {job_nodes} is below the {required} nodes the runs require: fail closed")
    scaled = {}
    for name, estimate in estimates.items():
        scaled[name] = estimate_stages.scale_stage_to_nodes(model, estimate, job_nodes)
        high = scaled[name]["ByPCGFactor"][worst]["StageSecondsEstimate"]
        low = scaled[name]["ByPCGFactor"][min(factors, key=float)]["StageSecondsEstimate"]
        cap = min(build_plan.round_up(build_plan.CAP_FACTOR * high), deadline)
        config_path = traces_record["Configs"][name]
        stages.append({"Name": f"{name}-dense", "Config": config_path, "Environment": {}, "CapSeconds": cap,
                       "MinimumSeconds": min(build_plan.round_up(low), cap), "Requires": []})
    per_node = max(estimate_stages.stage_per_node_used_gib(stage) for stage in scaled.values())
    peak_gb = max(stage["PalacePeakGBEstimate"] for stage in estimates.values())
    if job_nodes > 1:
        instance = estimate_stages.select_instance_for_nodes(profile, per_node, job_nodes, estimate_stages.per_node_guard_margin(model))
        if instance is None:
            raise ValueError(f"no instance's admission guard holds {per_node:.0f} GiB per node at {job_nodes} nodes: fail closed")
        instance["EstimatedPalacePeakGB"] = peak_gb
    else:
        instance = build_plan.select_instance(profile, peak_gb, model["PalaceGBPerGiB"])
    missing = [path for path, digest in pinned.items() if digest is None]
    total = {factor: sum(stage["ByPCGFactor"][factor]["StageSecondsEstimate"] for stage in scaled.values()) for factor in factors}
    walltime = int(walltime_seconds or profile["WalltimeSeconds"])
    plan = {"Version": build_plan.PLAN_VERSION, "Case": case_id or traces_record.get("Model"), "Job": "dense-twins", "JobKind": "dense-twin",
            "Purpose": (f"decision 457 (4): the (F) dense twins {list(runs)} of {dense_dir} as ONE Palace job on {job_nodes} node(s) "
                        f"({instance['Type']}); {STAGE_RULE}"),
            "Ranks": profile["Ranks"], "DeadlineSeconds": min(profile["DeadlineSeconds"], walltime - 900),
            "DeadlineMarginSeconds": profile["DeadlineMarginSeconds"],
            "Instance": instance, "MinimumMemAvailableBytes": instance["MinimumMemAvailableBytes"],
            "Binary": f"{remote_root}/{profile['BinaryPattern'].format(sha256=binary_sha256)}", "BinarySHA256": binary_sha256,
            "MPIExec": f"{remote_root}/{profile['MPIExecWrapper']}", "PinnedSHA256": pinned, "Stages": stages,
            "DenseTwinModel": {"Path": str(DENSE_TWIN_MODEL), "Executable": dense_model["Executable"]},
            "Estimate": {"Runs": scaled, "Nodes": job_nodes, "NodesRequired": required, "PerNodeUsedGiB": per_node,
                         "JobSecondsEstimateByPCGFactor": total,
                         "JobSecondsEstimateWithPreflightAndMargin": {factor: seconds * model["PreflightAndMarginFactor"] + model["PreflightSeconds"]
                                                                      for factor, seconds in total.items()}},
            "CapRule": (f"CapSeconds = {build_plan.CAP_FACTOR:g} x the run estimate at {worst}x the measured PCG counts rounded up to "
                        f"{build_plan.CAP_ROUNDING_SECONDS} s and bounded by the deadline; MinimumSeconds = the estimate at 1.0x rounded up"),
            "UnpinnedInputs": missing}
    if job_nodes > 1:
        plan.update(build_plan.multi_node_fields(profile, job_nodes, instance))
        plan["Purpose"] = f"{plan['Purpose']}; {build_plan.MULTI_NODE_REDUCTION_NOTE}"
    script = build_plan.render_job_script(profile=profile, remote_root=remote_root, remote_case_root=job_dir,
                                          runner=f"{job_dir}/run_stages.py", job_name=f"{profile['JobNamePrefix']}-dense-{(case_id or 'twins')}"[:64],
                                          walltime_seconds=walltime, instance_type=instance["Type"], job_directory=job_dir, nodes=job_nodes)
    return plan, script


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    cal = sub.add_parser("calibrate")
    cal.add_argument("--log", action="append", required=True, metavar="NAME=PALACE.LOG")
    cal.add_argument("--executable", required=True)
    cal.add_argument("--pbs-jobs", required=True)
    cal.add_argument("--out", type=Path, default=DENSE_TWIN_MODEL)
    pl = sub.add_parser("plan")
    pl.add_argument("--dense-dir", type=Path, required=True)
    pl.add_argument("--run", action="append", required=True)
    pl.add_argument("--entity-counts", type=Path, required=True)
    pl.add_argument("--binary-sha256", required=True)
    pl.add_argument("--remote-root", required=True)
    pl.add_argument("--job-dir", required=True, help="the remote directory the plan runs in (plan.json, job.pbs, run_stages.py copied there)")
    pl.add_argument("--out", type=Path, required=True)
    pl.add_argument("--nodes", type=int)
    pl.add_argument("--case")
    pl.add_argument("--walltime-seconds", type=int)
    pl.add_argument("--cost-model", type=Path, default=estimate_stages.COST_MODEL)
    pl.add_argument("--dense-twin-model", type=Path, default=DENSE_TWIN_MODEL)
    pl.add_argument("--cluster-profile", type=Path, default=estimate_stages.CLUSTER_PROFILE)
    args = parser.parse_args(argv)
    if args.command == "calibrate":
        logs = dict(item.split("=", 1) for item in args.log)
        model = calibrate(logs, executable=args.executable, pbs_jobs=args.pbs_jobs)
        args.out.write_text(json.dumps(model, indent=2) + "\n")
        print(json.dumps(model["Orders"], indent=1))
        return 0
    model = estimate_stages.load_cost_model(args.cost_model)
    profile = json.loads(args.cluster_profile.read_text())
    plan, script = plan_dense_twins(dense_dir=args.dense_dir, runs=args.run, counts=json.loads(args.entity_counts.read_text()),
                                    binary_sha256=args.binary_sha256, remote_root=args.remote_root, job_dir=args.job_dir, model=model,
                                    dense_model=json.loads(args.dense_twin_model.read_text()), profile=profile, nodes=args.nodes,
                                    walltime_seconds=args.walltime_seconds, case_id=args.case)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    (args.out / "job.pbs").write_text(script)
    print(json.dumps({"Nodes": plan.get("Nodes", 1), "Instance": plan["Instance"]["Type"], "Ranks": plan["Ranks"],
                      "Stages": [(stage["Name"], stage["CapSeconds"]) for stage in plan["Stages"]],
                      "PerNodeUsedGiB": plan["Estimate"]["PerNodeUsedGiB"], "UnpinnedInputs": plan["UnpinnedInputs"]}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
