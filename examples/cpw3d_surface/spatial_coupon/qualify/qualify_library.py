#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""`coupon-library qualify`: physics qualification of the coupons of a library-build.json
against their graded_v2 references (supervisor decision 48, second command).

Per passed coupon of the build record:
 1. inputs: the run config, the source traces and the zero-trace knots are derived from
    the CASE ITSELF (case_inputs.py: process-library.json materials / interface layers /
    interface index -> type map, the bound trace basis regenerated in the mesh frame,
    the identity mesh's $PhysicalNames, the manifest recipe's PhysicsRun Order and
    Linear.Tol); the identity mesh of the build record is hash-verified; the case's
    signature layers must be one upward layer (locate_sources.check_layers: a
    DownwardLayers / MultipleLayers case stops with a recorded ScopeGuard).  A reference
    campaign directory (--reference) is matched BY CONTENT (the frozen basis-contract.json
    digest) to `inputs-<key>/` and `case-<key>-fabricated/reducer/` (the reference
    matrices, when the reference completed); the reference's own config must equal the
    derived one apart from Model.Mesh, Problem.Output, the DataFile directory,
    Solver.Order and Solver.Linear.Tol (fail closed otherwise); `--reference none` runs
    the coupon on its own inputs (verdict PendingQualification or Failed, never Passed).
    A coupon with per-ring radial MA shells (every production build since decision 61a:
    build record `RadialShells`, the placement stage's census; a calibration relabel:
    `Relabel`) derives its base config against the parent MA labels - the reference
    comparison view - and runs the shell-expanded config (case_inputs.expand_radial_shells:
    one MA interface per shell), the reference matrices labeled by the reference's own
    interface map; a labeled Calibration.PhysicsRun.LinearTol deviation replaces the
    recipe tolerance;
 2. sources: every PrescribedPotential source of the derived config; the control
    sources by geometric class (locate_sources / classify_sources: one per class in
    priority order, cycling) unless --control-source names them;
 3. configs: worker / reducer at every --orders order on all sources, at every
    --controls order on the controls, the ordinary-path local-edge stage at the main
    order on the controls (build_configs.py); the remote paths follow the layout
    <remote root>/<run>/<case>/{mesh, inputs/traces, main/<prefix>-p<order>[-control]};
 4. estimate (estimate_stages.py, the cost model scaled by the exact H1 counts of the
    build record's entity counts) and the JOB POLICY (job_split.py, decision 61b): the
    main-order source set is split into N contiguous blocks run as N independent worker
    jobs archiving into the stage's one archive directory, then one reducer job on the
    union (N = 1: the single job of the recorded campaigns, byte-identical); the policy
    {speed | frugal | fixed} from --job-policy / --fixed-jobs, else the manifest's
    PhysicsRun.JobPolicy default, else frugal; fail closed only when even the maximal
    split (<= --max-jobs, the user job cap, the source count) does not fit the walltime
    or the Palace peak exceeds the node fraction; one plan.json per job with
    estimate-derived caps and pinned SHA-256 of the mesh, its configs and every trace;
    job.pbs from the cluster profile;
 5. --dry-run stops here (plans / configs / estimates / qualification-gates.json written,
    nothing contacted); otherwise up to --max-jobs of the run's jobs are queued / running
    at once (each qsub under the user job cap, recorded and counted at submission): a
    coupon's worker jobs, then - every worker status.json complete and the archive union
    counted (sources x ranks potential files) - its reducer job; every active job is
    polled read-only once per interval, and a coupon whose last job left the queue is
    fetched (never the archives), hash-verified CSV by CSV against the remote digests,
    matrix-validated (run_graded_library_case.validate_matrix: complete, symmetric,
    nonnegative), its remote archives deleted (recorded) and analyzed while the other
    jobs run; the next ready job takes the freed slot;
 6. qualification: comparisons vs the reference and between the p levels, class
    statistics, MA / MS / SA offsets, p-sequence controls, key sources, cost; the
    frozen gate table (qualification-gates.json, digest recorded) -> Passed / Failed, or
    PendingQualification when there are no reference matrices and every p-sequence
    control passed (never Passed; a failed control is Failed).  The main orders are the
    recipe's PhysicsRun order, then --orders, then the reference's own order when it
    differs; the same-order comparison is gated (the others are informational); a
    participation of an interface the coupon does not postprocess is NotApplicable
    (gates.py).  The MA of a shelled coupon is the SHARP-EDGE value (user decision 60(1),
    ma_tail.py): per source and order MA_raw = the shell sum, MA_tail = the `Consistent`
    remainder inside the innermost ring (top edge -2/3 anchored on ring 2 with the
    ring-1 factor measured on this run, bottom edge its own rings-2..4 law), MA_sharp =
    MA_raw + MA_tail, with alpha (rings 2-4) and its standard error, the ring-1 factor
    and the estimator spread recorded (comparison/ma-tail.{json,md}, ma-sharp-<stage>.csv);
    the p_MA gate and the p-sequence MA control evaluate p_MA_sharp (the gate table's
    Quantity / MAObservable) against the reference's sharp value - extrapolated by the same
    rule when its ring / edge sizing is recorded, else modelled from --reference-edge-size-nm
    by the eps^(1/3) law and annotated 'reference unextrapolated' (raw-vs-raw recorded
    next to sharp-vs-modelled); a coupon without shells is gated on the raw p_MA and says so.

Records: ROOT/library-qualification.json (per coupon: verdict per gate and class
offsets, PCG, node-h, x the reference cost, the job policy and split - N, blocks, per-job
estimate / actual / node-h, the coupon's critical path; library totals: node-h summed
over every job, critical-path wall clock from first submission to last fetch, jobs vs
the cap, splits, coupons stopped and why),
ROOT/qualification-gates.json (the table used), ROOT/process-library.json (the
Palace-facing Version-3 library, LIBRARY_RULE: the sources' header, every model's files
copied under ROOT/models/<slug>/, the per-source MA_raw / MA_sharp of a shelled coupon and
the qualification bound; LibraryQualified only when Passed; ThinMatrix null + NotLoadable
until the thin coupon exists), ROOT/process-library-preflight.json (the geometry-only
variant for palace --surface-response-preflight, PREFLIGHT_RULE).
"""
import argparse
import calendar
import csv
import hashlib
import json
import math
import re
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent
for path in (str(HERE), str(TOOLS)):
    if path not in sys.path:
        sys.path.insert(0, path)
import build_configs  # noqa: E402
import build_plan  # noqa: E402
import case_inputs  # noqa: E402
import classify_sources  # noqa: E402
import compare_matrices  # noqa: E402
import estimate_stages  # noqa: E402
import gates as gate_evaluation  # noqa: E402
import job_split  # noqa: E402
import key_sources  # noqa: E402
import locate_sources  # noqa: E402
import ma_ms_offsets  # noqa: E402
import ma_tail  # noqa: E402
import p_sequence  # noqa: E402
import reference_campaign  # noqa: E402
import remote as remote_side  # noqa: E402
import summarize_cost  # noqa: E402
from run_graded_library_case import validate_matrix  # noqa: E402

LIBRARY_QUALIFICATION_RECORD = "library-qualification.json"
PROCESS_LIBRARY_RECORD = "process-library.json"
GATES_COPY = "qualification-gates.json"
QUALIFICATION_VERSION = 1
DEFAULT_CONTROL_COUNT = 8
DEFAULT_MONITOR_INTERVAL = 90
DEFAULT_MONITOR_POLLS = 240
STATUS_QUALIFIED, STATUS_PENDING, STATUS_FAILED, STATUS_PLANNED, STATUS_SKIPPED, STATUS_UNSUPPORTED = (
    "qualified", "pending-qualification", "failed", "planned", "skipped", "unsupported-class")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n")


def parse_orders(text):
    orders = []
    for item in text.split(","):
        item = item.strip().lower()
        if not item:
            continue
        if not item.startswith("p") or not item[1:].isdigit():
            raise argparse.ArgumentTypeError(f"orders are p<N>, not {item!r}")
        orders.append(int(item[1:]))
    if len(set(orders)) != len(orders):
        raise argparse.ArgumentTypeError(f"repeated order in {text!r}")
    return orders


def parse_remote(text):
    if text is None:
        return None
    if ":" not in text:
        raise argparse.ArgumentTypeError("--remote expects HOST:ROOT")
    host, root = text.split(":", 1)
    if not host or not root.startswith("/"):
        raise argparse.ArgumentTypeError("--remote expects HOST:/absolute/root")
    return {"Host": host, "Root": root.rstrip("/")}


class CaseStop(Exception):
    """A per-coupon fail-closed stop (recorded as StoppedBy; the other coupons continue)."""

    def __init__(self, kind, message, **fields):
        super().__init__(message)
        self.record = {"Kind": kind, "Message": message, **fields}


def manifest_case(manifest, case_id):
    case = next((item for item in manifest["Cases"] if item["Id"] == case_id), None)
    if case is None:
        raise CaseStop("Manifest", f"case {case_id} is not in the manifest {manifest.get('Path')}")
    return case


def source_directory(manifest_path, manifest, case):
    directory = Path(case["Source"]["Directory"])
    if directory.is_absolute():
        return directory
    return (Path(manifest_path).parent / manifest["RepositoryRoot"]).resolve() / directory


ORDERS_RULE = ("main stages at the recipe's PhysicsRun order (the library order: cost coupon, local-edge stage, "
               "p-sequence main), then every --orders order, then the reference's own order when it differs (the "
               "same-order comparison is the gated one)")


def main_orders_of(recipe_order, orders, reference_order):
    """The main orders of a coupon: the recipe order, --orders, the reference order."""
    main_orders = [recipe_order]
    for order in list(orders or []) + ([reference_order] if reference_order is not None else []):
        if order not in main_orders:
            main_orders.append(order)
    return main_orders


def gated_order(main_orders, reference_order):
    """The gated main order: the reference order when the run solves it, else the first."""
    return reference_order if reference_order in main_orders else main_orders[0]


def stage_layout(prefix, main_orders, control_orders, source_count, control_count):
    """The run order: full stages at every main order, the control stages from the highest
    order down (the step to the higher order is secured before the cheaper low control, as
    the recorded campaigns did), the local-edge stage at the first main order."""
    layout = []
    control_orders = sorted(control_orders, reverse=True)
    for order in main_orders:
        layout.append({"Prefix": f"{prefix}-p{order}", "Kind": "response", "Order": order, "Sources": source_count,
                       "EstimateKey": f"p{order}-{source_count}", "Role": "main"})
    for order in control_orders:
        layout.append({"Prefix": f"{prefix}-p{order}-control", "Kind": "response", "Order": order, "Sources": control_count,
                       "EstimateKey": f"p{order}-control-{control_count}", "Role": "control"})
    layout.append({"Prefix": f"{prefix}-p{main_orders[0]}-local-edge", "Kind": "local-edge", "Order": main_orders[0],
                   "Sources": control_count, "EstimateKey": f"local-edge-p{main_orders[0]}-{control_count}", "Role": "local-edge"})
    return layout


def prepare_case(case_record, *, manifest_path, manifest, args, root, remote, profile, cost_model, gates, gates_digest):
    """Steps 1-4 of one coupon: the per-case record with the plan written locally."""
    case_id = case_record["Case"]
    case_root = root / case_id
    case_root.mkdir(parents=True, exist_ok=True)
    record = {"Case": case_id, "Status": None, "StoppedBy": None, "BuildStatus": case_record["Status"],
              "BuildPassed": case_record["Passed"], "Mesh": None, "Reference": None, "Sources": None, "Controls": None,
              "StagePrefix": None, "Stages": None, "Estimate": None, "Plan": None, "Configs": None,
              "Root": str(case_root), "Remote": None}
    if not case_record["Passed"]:
        raise CaseStop("Build", f"the build record marks {case_id} {case_record['Status']} (Passed false): not qualified")
    case = manifest_case(manifest, case_id)
    files = case["Source"]["Files"]
    if "BasisContract" not in files:
        raise CaseStop("Manifest", f"{case_id} binds no trace basis (no basis-contract.json): no sources to qualify")
    # 1. mesh (hash-verified now), the case's own run inputs, the reference (by content).
    identity = case_record["Variants"].get("identity")
    if not identity or not Path(identity["Path"]).is_file():
        raise CaseStop("Mesh", f"the build record's identity mesh of {case_id} is missing ({identity})")
    actual = sha256(identity["Path"])
    if actual != identity["SHA256"]:
        raise CaseStop("Mesh", f"identity mesh SHA256 {actual} differs from the build record {identity['SHA256']}")
    record["Mesh"] = {"Local": identity["Path"], "SHA256": identity["SHA256"], "Verified": True}
    directory = source_directory(manifest_path, manifest, case)
    # The per-ring MA shells: a production build records them under RadialShells (the
    # placement stage's census, decision 61a), a calibration relabel under Relabel.
    radial_shells = None
    shell_binding = case_record.get("RadialShells") or case_record.get("Relabel")
    if shell_binding:
        shells = shell_binding.get("Shells") or {}
        if not shells.get("Path") or not Path(shells["Path"]).is_file() or sha256(shells["Path"]) != shells.get("SHA256"):
            raise CaseStop("Mesh", f"the build record's radial-shell census of {case_id} is missing or its SHA256 differs ({shells})")
        radial_shells = json.loads(Path(shells["Path"]).read_text())
        if radial_shells.get("Mesh", {}).get("SHA256") != identity["SHA256"]:
            raise CaseStop("Mesh", f"the radial-shell census binds mesh {radial_shells.get('Mesh', {}).get('SHA256')}, the build "
                                   f"record {identity['SHA256']}")
        record["Mesh"]["RadialShells"] = {"Kind": shell_binding.get("Kind"), "Parent": shell_binding.get("Parent"),
                                          "ParentMesh": shell_binding.get("ParentMesh"), "Shells": shells,
                                          "ShellCount": len(radial_shells["Shells"]), "RingRadii": radial_shells.get("RingRadii"),
                                          "Binding": "RadialShells" if case_record.get("RadialShells") else "Relabel"}
    try:
        physics_run = case_inputs.physics_run_parameters(manifest, case)
        signature_layers = locate_sources.signature_layers(
            locate_sources.read_signature_rows(directory / files["Signature"]["Name"]))
        locate_sources.check_layers(signature_layers)
    except locate_sources.UnsupportedSourceGeometry as guard:
        raise CaseStop("ScopeGuard", str(guard), Guard=guard.guard, Layers=[list(layer) for layer in guard.layers],
                       Rule="locate_sources assigns the z-level roles of the sources for one upward process layer; "
                            "a DownwardLayers / MultipleLayers case is refused until the roles are assigned per layer band")
    except case_inputs.CaseInputError as error:
        raise CaseStop("Manifest", str(error))
    prefix = args.stage_prefix or case_id
    run_name = root.name
    remote_root = remote["Root"] if remote else "<remote-root>"
    remote_run = f"{remote_root}/{run_name}"
    remote_case = f"{remote_run}/{case_id}"
    remote_mesh = f"{remote_case}/mesh/identity-{identity['SHA256'][:12]}.msh"
    remote_traces = f"{remote_case}/inputs/traces"
    try:
        run_config, inputs = case_inputs.derive(case, directory, mesh_path=identity["Path"], physics_run=physics_run,
                                                out_dir=case_root / "inputs", mesh=remote_mesh, output_root=f"{remote_case}/main",
                                                traces_remote=remote_traces, radial_shells=radial_shells)
    except (case_inputs.CaseInputError, ValueError) as error:
        raise CaseStop("Inputs", f"the run inputs of {case_id} cannot be derived from its sources: {error}")
    write_json(case_root / "inputs" / "run-config.json", run_config)
    # A radial-shell relabel compares its BASE config (parent labels, one MA entry) with
    # the reference; the run config carries the shell entries (case_inputs.expand_radial_shells).
    comparison_config = inputs.pop("BaseConfig", None)
    if comparison_config is not None:
        write_json(case_root / "inputs" / "base-run-config.json", comparison_config)
        inputs["BaseConfigPath"] = str(case_root / "inputs" / "base-run-config.json")
        inputs["BaseConfigSHA256"] = sha256(case_root / "inputs" / "base-run-config.json")
    else:
        comparison_config = run_config
    record["Inputs"] = {**inputs, "Config": str(case_root / "inputs" / "run-config.json"),
                        "ConfigSHA256": sha256(case_root / "inputs" / "run-config.json")}
    interface_types = {int(index): name for index, name in inputs["Interfaces"].items()}
    interfaces = inputs["InterfaceTypes"]
    reference = None
    reference_interface_types = None
    if args.reference is not None:
        try:
            reference = reference_campaign.resolve(args.reference, files["BasisContract"]["SHA256"])
        except reference_campaign.ReferenceError as error:
            raise CaseStop("Reference", str(error))
        if reference is None:
            raise CaseStop("Reference", f"no reference campaign inputs bind the basis contract of {case_id} "
                                        f"({files['BasisContract']['SHA256'][:12]}...) under {args.reference}: pass "
                                        f"--reference none to run the coupon on its own inputs")
        reference_config = json.loads(Path(reference["Config"]).read_text())
        reference_interface_types = case_inputs.interface_types(reference_config)
        differences = case_inputs.config_differences(comparison_config, reference_config, ignore_solver=("Order", "Linear.Tol"))
        if differences:
            raise CaseStop("Reference", f"the reference config {reference['Config']} differs from the config derived from "
                                        f"the case's own sources outside Mesh / Output / DataFile / Order / Tol: {differences[:12]}",
                           Differences=differences)
        reference["ConfigEqualsDerived"] = {"ApartFrom": list(case_inputs.PATH_FIELDS) + ["Solver.Order", "Solver.Linear.Tol"],
                                            "ComparedConfig": ("the base config (parent labels) of the radial-shell relabel"
                                                               if comparison_config is not run_config else "the run config"),
                                            "ReferenceOrder": reference["ReferenceOrder"], "ReferenceLinearTol": reference["LinearTol"],
                                            "RunOrder": physics_run["Order"], "RunLinearTol": physics_run["LinearTol"]}
    record["Reference"] = reference if reference is not None else {
        "Rule": "--reference none: the coupon runs on its own inputs; without reference matrices the verdict is "
                "PendingQualification (every p-sequence control passed) or Failed, never Passed"}
    reference_order = reference["ReferenceOrder"] if reference else None
    # 2. sources and controls.
    sources = inputs["Sources"]
    indices = [source["Index"] for source in sources]
    terminals = [source["Index"] for source in sources if source["Terminal"]]
    zero_trace = inputs["ZeroTraceIndices"]
    trace_dir = Path(inputs["Traces"]["Directory"])
    comparison_dir = case_root / "comparison"
    comparison_dir.mkdir(exist_ok=True)
    try:
        rows, geometry = locate_sources.locate_directory(
            trace_dir, inputs["PlanViewBoundary"], comparison_dir / "source-locations.csv", signature_path=inputs["Signature"],
            retained_etch=inputs["RetainedEtch"], geometry_out=comparison_dir / "source-geometry.json")
    except locate_sources.UnsupportedSourceGeometry as guard:
        raise CaseStop("ScopeGuard", str(guard), Guard=guard.guard)
    locations = {row["index"]: {key: str(value) for key, value in row.items()} for row in rows}
    classes = classify_sources.classify_all(locations, zero_trace, terminals, thresholds=classify_sources.thresholds_of_gates(gates))
    free = [i for i in indices if i not in set(zero_trace) and i in locations]
    if args.control_source:
        controls = sorted(args.control_source)
        unknown = [i for i in controls if i not in indices]
        if unknown:
            raise CaseStop("Controls", f"--control-source {unknown} are not sources of {case_id}")
        control_rule = "explicit --control-source"
    else:
        controls = classify_sources.choose_controls(classes, args.control_count, free)
        control_rule = f"{args.control_count} by class (classify_sources.choose_controls: one per class in priority order, cycling)"
    record["Sources"] = {"Count": len(indices), "Indices": indices, "ZeroTrace": zero_trace, "Free": len(free),
                         "Terminals": terminals, "Classes": {str(i): name for i, (name, _) in classes.items()},
                         "Geometry": {key: geometry[key] for key in ("box", "z_levels", "per_z_level", "layers")}}
    record["Controls"] = {"Indices": controls, "Rule": control_rule,
                          "Classes": {str(i): classes[i][0] for i in controls if i in classes}}
    # 3. stage layout, remote layout, configs.
    main_orders = main_orders_of(physics_run["Order"], args.orders, reference_order)
    control_orders = list(args.controls)
    record["Orders"] = {"Main": [f"p{order}" for order in main_orders], "Controls": [f"p{order}" for order in control_orders],
                        "RecipeOrder": f"p{physics_run['Order']}", "RecipeLinearTol": physics_run["LinearTol"],
                        "ReferenceOrder": f"p{reference_order}" if reference_order is not None else None,
                        "Gated": f"p{gated_order(main_orders, reference_order)}", "Rule": ORDERS_RULE}
    layout = stage_layout(prefix, main_orders, control_orders, len(indices), len(controls))
    record["StagePrefix"] = prefix
    record["Stages"] = layout
    # 4. estimate (the single-job view, then the split under the job policy: fail closed
    # only when even the maximal split does not fit - decision 61b).
    counts = (case_record.get("H1") or {}).get("EntityCounts")
    counts_origin = "library-build.json H1.EntityCounts"
    if counts is None:
        counts = estimate_stages.entity_counts_of_mesh(identity["Path"])
        counts_origin = "computed from the identity mesh (the build record carries no entity counts)"
    stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in layout if item["Kind"] == "response"]
    local_edge = next((item for item in layout if item["Kind"] == "local-edge"), None)
    reducer_block_size = reducer_block_size_of(args, physics_run)
    frozen_binary = frozen_binary_of(args, physics_run)
    try:
        estimate = estimate_stages.estimate(counts, stages, model=cost_model, profile=profile,
                                            local_edge=(local_edge["EstimateKey"], local_edge["Order"], local_edge["Sources"]),
                                            block_size=reducer_block_size["Value"])
    except ValueError as error:
        raise CaseStop("Estimate", str(error))
    estimate["EntityCountsOrigin"] = counts_origin
    estimate["ReducerBlockSizeOrigin"] = reducer_block_size["Origin"]
    policy = job_policy_of(args, physics_run, profile)
    try:
        split = job_split.plan_split(indices=indices, layout=layout, estimate=estimate, policy=policy, model=cost_model,
                                     profile=profile)
    except ValueError as error:
        raise CaseStop("JobPolicy", str(error), Policy=policy)
    estimate["Split"] = split
    write_json(case_root / "preflight" / "stage-estimate.json", estimate)
    record["JobPolicy"] = policy
    record["ReducerBlockSize"] = reducer_block_size
    record["FrozenExecutable"] = frozen_binary
    record["Estimate"] = {"Path": str(case_root / "preflight" / "stage-estimate.json"), "FitsOneJob": estimate["FitsOneJob"],
                         "Decision": estimate["Decision"], "H1ByOrder": estimate["Mesh"]["H1ByOrder"],
                         "JobSecondsEstimateByPCGFactor": estimate["JobSecondsEstimateByPCGFactor"],
                         "JobSecondsEstimateWithPreflightAndMargin": estimate["JobSecondsEstimateWithPreflightAndMargin"],
                         "MaxPalacePeakGBEstimate": estimate["MaxPalacePeakGBEstimate"],
                         "ReducerBlockSize": estimate["ReducerBlockSize"],
                         "MainStageNodeHoursEstimate": {
                             factor: value["StageSecondsEstimate"] * profile["Nodes"] / 3600.0
                             for factor, value in estimate["Stages"][layout[0]["EstimateKey"]]["ByPCGFactor"].items()}}
    record["Split"] = {key: split[key] for key in ("N", "Fits", "Decision", "ControlsJob", "LargestSplit", "Candidates",
                                                    "CriticalPathEstimateSeconds", "NodeSecondsEstimate", "WorstPCGFactor")}
    record["Split"]["Blocks"] = [len(block) for block in split["Blocks"]] if split["Blocks"] else None
    record["Split"]["Rule"] = job_split.SPLIT_RULE
    fits_memory = estimate["MaxPalacePeakGBEstimate"] / cost_model["PalaceGBPerGiB"] + 60 < cost_model["NodeFitFraction"] * profile["NodeGiB"]
    if not fits_memory:
        raise CaseStop("Estimate", estimate["Decision"], Estimate=record["Estimate"])
    if not split["Fits"]:
        raise CaseStop("Estimate", split["Decision"], Estimate=record["Estimate"], Split=record["Split"])
    # 3. configs (the split's blocks at every main order) and one plan / job script per job.
    config_digests = {}
    local_configs = {}
    for item in layout:
        directory_ = case_root / "main" / item["Prefix"]
        if item["Role"] == "main" and split["N"] > 1:
            config_digests[item["Prefix"]] = build_configs.write_split_stage(
                directory_, run_config, remote_mesh, f"{remote_case}/main/{item['Prefix']}", split["Blocks"], remote_traces,
                order=item["Order"])
        elif item["Kind"] == "response":
            subset = indices if item["Role"] == "main" else controls
            worker, reducer = build_configs.derive(run_config, remote_mesh, f"{remote_case}/main/{item['Prefix']}",
                                                   subset, remote_traces, order=item["Order"])
            config_digests[item["Prefix"]] = build_configs.write_stage(directory_, worker, reducer)
        else:
            config, _ = build_configs.derive(run_config, remote_mesh, f"{remote_case}/main/{item['Prefix']}", controls,
                                             remote_traces, order=item["Order"], save_local_edge_energy=True)
            config_digests[item["Prefix"]] = build_configs.write_local_edge_stage(
                directory_, config, f"{remote_case}/main/{item['Prefix']}/output")
        local_configs[item["Prefix"]] = str(directory_)
    record["Configs"] = {"Directories": local_configs, "SHA256": config_digests,
                         "DerivedFrom": "the case's own sources (case_inputs.derive) and the recipe PhysicsRun",
                         "RunConfig": record["Inputs"]["Config"], "RunConfigSHA256": record["Inputs"]["ConfigSHA256"],
                         "ReferenceConfig": reference["Config"] if reference else None,
                         "ReferenceConfigSHA256": reference["ConfigSHA256"] if reference else None,
                         "ReferenceConfigRole": reference["ConfigRole"] if reference else None,
                         "ReferenceOrder": reference_order, "LinearTol": physics_run["LinearTol"], "Order": physics_run["Order"]}
    trace_pins = {f"{remote_traces}/{source['Name']}": source["SHA256"] for source in sources}
    binary = f"{remote_root}/{profile['BinaryPattern'].format(sha256=frozen_binary['SHA256'])}"
    mpiexec = f"{remote_root}/{profile['MPIExecWrapper']}"
    purpose = (f"coupon-library qualify of {case_id} (identity mesh {identity['SHA256']}, {len(indices)} sources, "
               f"{len(zero_trace)} contract zero-trace knots; config derived from the case's sources at the recipe Order "
               f"{physics_run['Order']}, Linear.Tol {physics_run['LinearTol']}) "
               + (f"against the graded_v2 reference inputs {reference['Key']} (reference Order {reference['ReferenceOrder']}, "
                  f"Linear.Tol {reference['LinearTol']})" if reference else "without a reference (--reference none)")
               + f"; stages {[item['Prefix'] for item in layout]}; controls {controls} ({control_rule}); {estimate['Decision']}")
    factors = [f"{factor:.1f}" for factor in cost_model["PCGFactors"]]
    jobs = []
    if split["N"] == 1:
        plan = build_plan.build_plan(case_id=case_id, remote_case_root=remote_case,
                                     mesh={"Remote": remote_mesh, "SHA256": identity["SHA256"], "Local": identity["Path"]},
                                     stage_layout=layout, estimate=estimate, config_digests=config_digests, trace_pins=trace_pins,
                                     profile=profile, binary=binary, binary_sha256=frozen_binary["SHA256"], mpiexec=mpiexec,
                                     purpose=purpose, factors=factors, reducer_block_size=reducer_block_size["Value"])
        write_json(case_root / "main" / "plan.json", plan)
        job_script = build_plan.render_job_script(profile=profile, remote_root=remote_root, remote_case_root=remote_case,
                                                  runner=f"{remote_run}/run_stages.py",
                                                  job_name=f"{profile['JobNamePrefix']}-{case_id}"[:64],
                                                  walltime_seconds=profile["WalltimeSeconds"])
        (case_root / "main" / "job.pbs").write_text(job_script)
        jobs.append({"Name": "single", "Kind": "single", "Block": 1, "Sources": indices, "Requires": [],
                     "Directory": str(case_root / "main"), "RemoteDirectory": f"{remote_case}/main",
                     "SubmissionRecord": str(case_root / "submission.json"), "StageNames": [stage["Name"] for stage in plan["Stages"]],
                     "Estimate": split["Jobs"][0]["SecondsEstimateWithPreflightAndMargin"], "Plan": str(case_root / "main" / "plan.json")})
    else:
        for split_job in split["Jobs"]:
            name = split_job["Name"]
            directory_ = case_root / "main" / "jobs" / name
            remote_directory = f"{remote_case}/main/jobs/{name}"
            plan = build_plan.build_job_plan(case_id=case_id, job_name=name, remote_case_root=remote_case,
                                             mesh={"Remote": remote_mesh, "SHA256": identity["SHA256"], "Local": identity["Path"]},
                                             stage_layout=layout, split_job=split_job, estimate=estimate, config_digests=config_digests,
                                             trace_pins=trace_pins, profile=profile, binary=binary,
                                             binary_sha256=frozen_binary["SHA256"], mpiexec=mpiexec,
                                             purpose=f"{purpose}; job {name} of the split {split['Decision']}", factors=factors,
                                             reducer_block_size=reducer_block_size["Value"])
            write_json(directory_ / "plan.json", plan)
            job_script = build_plan.render_job_script(profile=profile, remote_root=remote_root, remote_case_root=remote_case,
                                                      runner=f"{remote_run}/run_stages.py",
                                                      job_name=f"{profile['JobNamePrefix']}-{case_id}-{name}"[:64],
                                                      walltime_seconds=profile["WalltimeSeconds"], job_directory=remote_directory)
            (directory_ / "job.pbs").write_text(job_script)
            jobs.append({"Name": name, "Kind": split_job["Kind"], "Block": split_job["Block"], "Sources": split_job["Sources"],
                         "Requires": ([job["Name"] for job in split["Jobs"] if job["Kind"] == "worker"]
                                      if split_job["Kind"] == "reducer" else []),
                         "Directory": str(directory_), "RemoteDirectory": remote_directory,
                         "SubmissionRecord": str(case_root / f"submission-{name}.json"), "StageNames": plan["StageNames"],
                         "Estimate": split_job["SecondsEstimateWithPreflightAndMargin"], "Plan": str(directory_ / "plan.json")})
    record["Plan"] = {"Path": jobs[0]["Plan"] if split["N"] == 1 else str(case_root / "main" / "jobs"),
                      "Pins": len(plan["PinnedSHA256"]),
                      "StageNames": [name for job in jobs for name in job["StageNames"]],
                      "Caps": {stage["Name"]: [stage["CapSeconds"], stage["MinimumSeconds"]]
                               for job in jobs for stage in json.loads(Path(job["Plan"]).read_text())["Stages"]},
                      "JobScript": str(case_root / "main" / "job.pbs") if split["N"] == 1 else None}
    # The job records are live: submission, monitor, status and times are written into
    # them as the run proceeds (every checkpoint carries them).
    record["Jobs"] = jobs
    record["Remote"] = {"Root": remote_root, "Run": remote_run, "Case": remote_case, "Mesh": remote_mesh,
                        "Traces": remote_traces, "Binary": binary, "MPIExec": mpiexec}
    record["Status"] = STATUS_PLANNED
    return record, {"jobs": jobs, "split": split, "sources": sources, "locations": locations, "classes": classes,
                    "layout": layout, "reference": reference, "reference_order": reference_order, "interfaces": interfaces,
                    "interface_types": interface_types, "reference_interface_types": reference_interface_types,
                    "controls": controls, "zero_trace": zero_trace, "identity": identity,
                    "radial_shells": inputs.get("RadialShells"),
                    "reference_edge_size": (getattr(args, "reference_edge_size_nm", None) / 1000.0
                                            if getattr(args, "reference_edge_size_nm", None) else None)}


def reducer_block_size_of(args, physics_run):
    """The reducer's PALACE_RESPONSE_BLOCK_SIZE: --reducer-block-size, else the manifest's
    PhysicsRun.ReducerBlockSize, else build_plan.DEFAULT_REDUCER_BLOCK_SIZE (decision 62(1))."""
    if args.reducer_block_size is not None:
        value, origin = args.reducer_block_size, "--reducer-block-size"
    elif physics_run.get("ReducerBlockSize") is not None:
        value, origin = physics_run["ReducerBlockSize"], "manifest ProductionRecipe.PhysicsRun.ReducerBlockSize"
    else:
        value, origin = build_plan.DEFAULT_REDUCER_BLOCK_SIZE, "built-in default build_plan.DEFAULT_REDUCER_BLOCK_SIZE"
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise CaseStop("ReducerBlockSize", f"the reducer block size must be an integer >= 1, not {value!r} ({origin})")
    return {"Value": value, "Origin": origin, "Rule": build_plan.REDUCER_BLOCK_SIZE_RULE}


def frozen_binary_of(args, physics_run):
    """The frozen Palace executable's SHA-256: --frozen-binary-sha256, else the manifest's
    PhysicsRun.FrozenExecutable, else build_plan.DEFAULT_FROZEN_BINARY_SHA256 (decision 63)."""
    if args.frozen_binary_sha256 is not None:
        value, origin = args.frozen_binary_sha256, "--frozen-binary-sha256"
    elif physics_run.get("FrozenExecutableSHA256") is not None:
        value, origin = physics_run["FrozenExecutableSHA256"], "manifest ProductionRecipe.PhysicsRun.FrozenExecutable"
    else:
        value, origin = build_plan.DEFAULT_FROZEN_BINARY_SHA256, "built-in default build_plan.DEFAULT_FROZEN_BINARY_SHA256"
    if not isinstance(value, str) or not value:
        raise CaseStop("FrozenExecutable", f"the frozen executable digest must be a non-empty string, not {value!r} ({origin})")
    return {"SHA256": value, "Origin": origin, "Rule": build_plan.FROZEN_BINARY_RULE}


def job_policy_of(args, physics_run, profile):
    """The coupon's job policy: --job-policy (and --fixed-jobs), else the manifest's
    PhysicsRun.JobPolicy default, else frugal; MaxJobs = --max-jobs, the walltime and
    the user job cap from the cluster profile."""
    manifest_default = physics_run.get("JobPolicy") or {}
    mode = args.job_policy or manifest_default.get("Mode") or job_split.DEFAULT_MODE
    fixed = args.fixed_jobs if args.fixed_jobs is not None else manifest_default.get("FixedJobs")
    origin = ("--job-policy" if args.job_policy else
              "manifest ProductionRecipe.PhysicsRun.JobPolicy" if manifest_default.get("Mode") else
              f"built-in default {job_split.DEFAULT_MODE}")
    try:
        return job_split.normalize_policy(mode, max_jobs=args.max_jobs, walltime_seconds=profile["WalltimeSeconds"],
                                          fixed_jobs=fixed, user_job_cap=profile["UserJobCap"], origin=origin)
    except ValueError as error:
        raise CaseStop("JobPolicy", str(error))


def upload_case(record, context, *, remote, profile):
    """Mesh, traces (every DataFile) and the main/ directory to the remote case root."""
    case_root = Path(record["Root"])
    staging = case_root / "upload"
    if staging.exists():
        shutil.rmtree(staging)
    (staging / "mesh").mkdir(parents=True)
    (staging / "inputs" / "traces").mkdir(parents=True)
    shutil.copyfile(context["identity"]["Path"], staging / "mesh" / Path(record["Remote"]["Mesh"]).name)
    for source in context["sources"]:
        shutil.copyfile(source["Path"], staging / "inputs" / "traces" / source["Name"])
    shutil.copytree(case_root / "main", staging / "main")
    commands = [remote_side.upload(remote["Host"], staging, record["Remote"]["Case"]),
                remote_side.upload_file(remote["Host"], HERE / "run_stages.py", record["Remote"]["Run"] + "/run_stages.py")]
    shutil.rmtree(staging)
    return {"Commands": commands, "UTC": remote_side.utc()}


def wait(seconds, slice_seconds=10):
    """Sleep in short slices against the wall clock (one long sleep of an idle driver was
    observed not to return on macOS; the slices keep the driver's timers short)."""
    end = time.time() + seconds
    while True:
        remaining = end - time.time()
        if remaining <= 0:
            return
        time.sleep(min(slice_seconds, remaining))


def plans_text(paths):
    """The plans of a coupon's jobs concatenated in path order (the --resume identity check
    reads the recorded plans by glob and the derived ones from the job records: one order)."""
    return "\n".join(path.read_text() for path in sorted(paths, key=str))


def job_plans_text(record):
    return plans_text(Path(job["Plan"]) for job in record["Jobs"])


def resume_submissions(record, context, plans_before):
    """Under --resume: adopt the submissions a previous driver recorded for this coupon's
    jobs (<case>/submission.json or <case>/submission-<job>.json) when the plans are
    byte-identical to the ones just derived; returns the adopted jobs (their Submission
    set) - the jobs without a record stay pending."""
    case_root = Path(record["Root"])
    if plans_before is None:
        return []
    if plans_before != job_plans_text(record):
        raise CaseStop("Resume", f"the plans derived now differ from the plans the recorded submissions under {case_root} "
                                 f"ran: not resumable (inputs or tooling changed)")
    adopted = []
    for job in context["jobs"]:
        submission_path = Path(job["SubmissionRecord"])
        if not submission_path.is_file():
            continue
        submission = json.loads(submission_path.read_text())
        job["Submission"] = submission
        job["Monitor"] = {"Polls": 0, "LastJobState": None, "LastStages": None, "Resumed": True}
        job["SubmittedAt"] = calendar.timegm(time.strptime(submission["UTC"], "%Y-%m-%dT%H:%M:%SZ"))
        job["Upload"] = {"Resumed": True, "SubmissionRecord": str(submission_path)}
        if job["Kind"] == "single":
            record["Submission"] = submission
            record["Monitor"] = job["Monitor"]
            record["SubmittedAt"] = job["SubmittedAt"]
        adopted.append(job)
    if adopted:
        record["Upload"] = {"Resumed": True, "SubmissionRecords": [job["SubmissionRecord"] for job in adopted]}
        record.setdefault("Monitor", {"Polls": 0, "LastJobState": None, "LastStages": None, "Resumed": True})
    return adopted


def submit_job(record, context, job, *, remote, profile, log):
    """Step 5a: upload the coupon (once, before its first job) and qsub one job under the
    user job cap (recorded per job)."""
    remote_case = record["Remote"]["Case"]
    if not record.get("Upload"):
        record["Upload"] = upload_case(record, context, remote=remote, profile=profile)
    try:
        submission = remote_side.submit(remote["Host"], profile["PBSBin"], f"{job['RemoteDirectory']}/job.pbs",
                                        job["RemoteDirectory"], job_cap=profile["UserJobCap"])
    except RuntimeError as error:
        raise CaseStop("JobBudget", str(error), Job=job["Name"])
    job["Submission"] = submission
    job["Monitor"] = {"Polls": 0, "LastJobState": None, "LastStages": None}
    job["SubmittedAt"] = time.time()
    write_json(job["SubmissionRecord"], submission)
    if job["Kind"] == "single":
        record["Submission"] = submission
        record["Monitor"] = job["Monitor"]
    log(f"{record['Case']}: submitted {job['Name']} {submission['Job']} ({submission['UserJobsBefore']} user jobs before, "
        f"cap {profile['UserJobCap']}; {remote_case})")
    return submission


def poll_job(record, job, *, remote, profile, log):
    """One read-only poll of a submitted job; True when it has left the queue (a held job -
    PBS H, the SOCA capacity retry - is still in the queue: remote.IN_QUEUE_STATES)."""
    submission = job["Submission"]
    poll = remote_side.poll(remote["Host"], profile["PBSBin"], submission["Job"], f"{job['RemoteDirectory']}/status.json")
    job["Monitor"]["Polls"] += 1
    job["Monitor"]["LastPollUTC"] = poll["UTC"]
    if not poll.get("Reachable", True):
        # A failed ssh round trip says nothing about the job: the poll counts against the
        # budget, the job stays active and the failure is recorded (never "left the queue").
        job["Monitor"]["TransportFailures"] = job["Monitor"].get("TransportFailures", 0) + 1
        log(f"== {poll['UTC']} {record['Case']} {job['Name']} job {submission['Job']} unreachable (ssh rc {poll['SSHReturnCode']}; "
            f"transport failure {job['Monitor']['TransportFailures']}, job kept active)")
        return False
    stages = ([(s["Name"], s["State"], round(s.get("WallSeconds", 0))) for s in poll["Status"]["Stages"]] if poll["Status"] else None)
    job["Monitor"]["LastJobState"] = poll["JobState"]
    job["Monitor"]["LastStages"] = stages
    held = re.search(r"(comment = .*|error_message = .*)", poll["QStat"]) if poll["JobState"] == "H" else None
    log(f"== {poll['UTC']} {record['Case']} {job['Name']} job {submission['Job']} state {poll['JobState']} stages {stages}"
        + (f" (held, still queued: {held.group(1).strip()})" if held else ""))
    return not remote_side.in_queue(poll["JobState"])


def poll_case(record, *, remote, profile, log):
    """The single-job poll of a coupon (the recorded campaigns' driver interface)."""
    job = {"Name": "single", "Submission": record["Submission"], "Monitor": record["Monitor"],
           "RemoteDirectory": f"{record['Remote']['Case']}/main"}
    return poll_job(record, job, remote=remote, profile=profile, log=log)


def status_failures(status):
    """The incomplete stages and PCG non-convergences of a runner status (empty = complete)."""
    incomplete = [stage["Name"] for stage in status["Stages"] if stage["State"] != "complete"]
    nonconvergence = {stage["Name"]: stage.get("Parsed", {}).get("Nonconvergence") for stage in status["Stages"]
                      if stage.get("Parsed", {}).get("Nonconvergence")}
    return incomplete, nonconvergence


def complete_worker_job(record, context, job, *, remote, profile):
    """A worker job of a split coupon left the queue: its status.json (read now, read-only)
    must be complete - every stage complete, no PCG non-convergence - else the coupon
    stops (fail closed; the coupon's other jobs are left to finish on their own and are
    recorded)."""
    status = remote_side.read_json(remote["Host"], f"{job['RemoteDirectory']}/status.json")
    job["FinishedAt"] = time.time()
    if status is None:
        raise CaseStop("Stages", f"{job['Name']} job {job['Submission']['Job']} left the queue without a status.json "
                                 f"({job['RemoteDirectory']})", Job=job["Name"], Submission=job["Submission"])
    incomplete, nonconvergence = status_failures(status)
    job["Status"] = {"State": status["State"], "TotalSeconds": status.get("TotalSeconds"), "Host": status.get("Host"),
                     "Stages": {stage["Name"]: (stage["State"], stage.get("WallSeconds")) for stage in status["Stages"]}}
    if status["State"] != "complete" or incomplete or nonconvergence:
        raise CaseStop("Stages", f"{job['Name']} job {job['Submission']['Job']}: runner state {status['State']}: incomplete "
                                 f"{incomplete}, PCG non-convergence {nonconvergence}", Job=job["Name"], Submission=job["Submission"])


def verify_archive_union(record, context, *, remote, profile):
    """Every worker job completed: the archive directory of every main stage holds one
    potential file per (source, rank) - the union the reducer job reduces (recorded)."""
    union = {}
    expected = len(context["sources"]) * profile["Ranks"]
    for item in context["layout"]:
        if item["Role"] != "main":
            continue
        directory = f"{record['Remote']['Case']}/main/{item['Prefix']}/archive"
        found = remote_side.count_archive_potentials(remote["Host"], directory)
        union[item["Prefix"]] = {"Directory": directory, "Expected": expected, "Found": found, "OK": found == expected,
                                 "Rule": "sources x ranks files source-*-rank-*-V.bin (one potential per source and rank)"}
    record["ArchiveUnion"] = {"UTC": remote_side.utc(), "Stages": union}
    write_json(Path(record["Root"]) / "archive-union.json", record["ArchiveUnion"])
    short = {prefix: (value["Found"], value["Expected"]) for prefix, value in union.items() if not value["OK"]}
    if short:
        raise CaseStop("ArchiveUnion", f"the archive union is incomplete (found, expected): {short}", ArchiveUnion=union)


def finish_case(record, context, *, remote, profile):
    """Step 5b after the coupon's last job left the queue: fetch (never the archives),
    hash-verify every CSV against the remote, validate every matrix, delete the remote
    archives (recorded).  Returns the local results directory (results/main)."""
    case_root = Path(record["Root"])
    remote_case = record["Remote"]["Case"]
    jobs = context["jobs"]
    results = case_root / "results"
    results.mkdir(exist_ok=True)
    try:
        fetch_command = remote_side.fetch(remote["Host"], f"{remote_case}/main", results / "main")
    except subprocess.CalledProcessError as error:
        raise CaseStop("Fetch", f"rsync of {remote_case}/main returned {error.returncode} (transport failure; the job's "
                                f"results stay on the remote: run again with --resume)",
                       Submission=[job.get("Submission") for job in jobs], Command=error.cmd)
    record["Fetch"] = {"Command": fetch_command, "UTC": remote_side.utc(),
                       "QStatHistory": {job["Name"]: remote_side.qstat_history(remote["Host"], profile["PBSBin"], job["Submission"]["Job"])
                                        for job in jobs}}
    for job in jobs:
        suffix = "" if job["Kind"] == "single" else f"-{job['Name']}"
        (results / f"qstat-xf{suffix}.txt").write_text(record["Fetch"]["QStatHistory"][job["Name"]])
    if len(jobs) == 1:
        record["Fetch"]["QStatHistory"] = record["Fetch"]["QStatHistory"]["single"]
    statuses = {}
    for job in jobs:
        status_path = results / "main" / Path(job["RemoteDirectory"]).relative_to(f"{remote_case}/main") / "status.json"
        if not status_path.is_file():
            raise CaseStop("Fetch", f"no status.json fetched for {job['Name']} ({status_path})", Submission=job.get("Submission"))
        status = json.loads(status_path.read_text())
        incomplete, nonconvergence = status_failures(status)
        if status["State"] != "complete" or incomplete or nonconvergence:
            raise CaseStop("Stages", f"{job['Name']}: runner state {status['State']}: incomplete {incomplete}, PCG non-convergence "
                                     f"{nonconvergence}", Submission=job.get("Submission"), StatusPath=str(status_path))
        statuses[job["Name"]] = status
    # Hash-verify every fetched CSV against the remote digests, validate every matrix.
    csv_local = sorted(path for path in (results / "main").rglob("*.csv"))
    remote_paths = [f"{remote_case}/main/{path.relative_to(results / 'main')}" for path in csv_local]
    remote_digests = remote_side.remote_sha256(remote["Host"], remote_paths)
    verification = {}
    mismatches = []
    for local, remote_path in zip(csv_local, remote_paths):
        local_digest = sha256(local)
        ok = remote_digests.get(remote_path) == local_digest
        verification[str(local.relative_to(results))] = {"Local": local_digest, "Remote": remote_digests.get(remote_path), "OK": ok}
        if not ok:
            mismatches.append(str(local))
    record["ResultDigests"] = verification
    write_json(case_root / "result-csv-sha256.json", verification)
    if mismatches:
        raise CaseStop("Verification", f"fetched CSVs differ from the remote digests: {mismatches}")
    matrices = {}
    for item in context["layout"]:
        if item["Kind"] != "response":
            continue
        reducer = results / "main" / item["Prefix"] / "reducer"
        subset = [source["Index"] for source in context["sources"]] if item["Role"] == "main" else context["controls"]
        try:
            matrices[item["Prefix"]] = {"domain": validate_matrix(str(reducer / "domain-response-matrix.csv"), subset, "domain"),
                                        "surface": validate_matrix(str(reducer / "surface-response-matrix.csv"), subset, "surface")}
        except (ValueError, FileNotFoundError) as error:
            raise CaseStop("MatrixValidation", f"{item['Prefix']}: {error}")
    record["MatrixValidation"] = {prefix: {kind: {"Rows": value["Rows"], "BasisSize": value["BasisSize"]}
                                           for kind, value in kinds.items()} for prefix, kinds in matrices.items()}
    write_json(case_root / "matrix-validation.json", matrices)
    archives = [f"{remote_case}/main/{item['Prefix']}/archive" for item in context["layout"] if item["Kind"] == "response"]
    record["ArchiveDeletion"] = remote_side.delete_archives(remote["Host"], archives)
    write_json(case_root / "remote-archive-deletion.json", record["ArchiveDeletion"])
    context["statuses"] = statuses
    return results


def analyze_case(record, context, results, *, gates, gates_digest, profile):
    """Step 6: comparisons, classes, offsets, p-sequence, key sources, cost, gates."""
    case_root = Path(record["Root"])
    comparison_dir = case_root / "comparison"
    reference = context["reference"]
    reference_order, interface_types = context["reference_order"], context["interface_types"]
    # The reference labels its surface matrix by its own config (one MA entry); a radial-
    # shell run labels its own by the shell entries - both sum per type.
    reference_interface_types = context.get("reference_interface_types") or interface_types
    layout, controls, zero_trace = context["layout"], context["controls"], context["zero_trace"]
    locations, classes = context["locations"], context["classes"]
    mains = [item for item in layout if item["Role"] == "main"]
    main = mains[0]
    gated = next(item for item in mains if item["Order"] == gated_order([item["Order"] for item in mains], reference_order))
    main_dir = results / "main" / main["Prefix"] / "reducer"
    record["MainReducer"] = str(main_dir)
    control_dirs = {item["Order"]: results / "main" / item["Prefix"] / "reducer" for item in layout if item["Role"] == "control"}
    reference_dir = Path(reference["Results"]) if reference and reference.get("Results") else None
    comparisons = {}
    labels = {}
    if reference_dir is not None:
        for item in mains:
            labels[f"{item['Prefix']}-vs-reference"] = (reference_dir, results / "main" / item["Prefix"] / "reducer",
                                                       f"reference p{reference_order}", f"{item['Prefix']} (all sources)")
        for order, directory in control_dirs.items():
            labels[f"{main['Prefix'].rsplit('-p', 1)[0]}-p{order}-control-vs-reference"] = (
                reference_dir, directory, f"reference p{reference_order}", f"p{order} control ({len(controls)} sources)")
    for order, directory in control_dirs.items():
        if order > main["Order"]:
            labels[f"{main['Prefix']}-vs-p{order}-control"] = (main_dir, directory, main["Prefix"], f"p{order} control")
        else:
            labels[f"p{order}-control-vs-{main['Prefix']}"] = (directory, main_dir, f"p{order} control", main["Prefix"])
    # The sharp-edge MA extrapolation (decision 61a): the tail of every reducer directory
    # from its own shells; the reference's side extrapolated / modelled / raw as recorded.
    shell_map = ma_tail.shell_map_of(context["radial_shells"]) if context.get("radial_shells") else None
    tails_cache = {}

    def tails_of(directory):
        if shell_map is None:
            return None
        key = str(directory)
        if key not in tails_cache:
            tails_cache[key] = ma_tail.tails(directory, shell_map, interface_types)
        return tails_cache[key]

    gated_dir = results / "main" / gated["Prefix"] / "reducer"
    reference_ma_side = None
    if reference_dir is not None and shell_map is not None:
        reference_ma_side = ma_tail.reference_side(
            ma_ms_offsets.reference_p_ma(reference_dir, reference_interface_types), tails_of(gated_dir),
            reference_edge_size=context.get("reference_edge_size"), run_inner_radius=context["radial_shells"]["RingRadii"][0])
    for label, (ref, run, ref_label, run_label) in labels.items():
        against_reference = ref == reference_dir
        comparisons[label] = compare_matrices.write_comparison(
            ref, run, comparison_dir / label, zero_trace_indices=zero_trace, locations=locations, reference_label=ref_label,
            run_label=run_label, max_entry_rows=200, interface_types=interface_types,
            reference_interface_types=(reference_interface_types if against_reference else interface_types),
            run_ma_tails=(tails_of(run) if against_reference else None),
            reference_ma_side=(reference_ma_side if against_reference else None))
    class_stats = classify_sources.write_class_report(classes, locations, list(comparisons.items()),
                                                      comparison_dir / "source-classes.csv", comparison_dir / "class-statistics.md")
    main_label = f"{gated['Prefix']}-vs-reference"
    informational = [f"{item['Prefix']}-vs-reference" for item in mains if item is not gated and reference_dir is not None]
    orders_of = {"main": main["Order"]}
    runs = {"main": main_dir, "low": None, "high": None, "ref": reference_dir}
    lower = [order for order in control_dirs if order < main["Order"]]
    higher = [order for order in control_dirs if order > main["Order"]]
    if lower:
        orders_of["low"] = max(lower)
        runs["low"] = control_dirs[max(lower)]
    if higher:
        orders_of["high"] = min(higher)
        runs["high"] = control_dirs[min(higher)]
    p_sequence_summary = p_sequence.write_p_sequence(runs, orders_of, controls, comparison_dir / "p-sequence-controls.md",
                                                     comparison_dir / "p-sequence-controls.json",
                                                     title=f"p-sequence controls of {record['Case']}", interface_types=interface_types,
                                                     reference_interface_types=reference_interface_types,
                                                     ma_tails=({key: tails_of(path) for key, path in runs.items()
                                                                if path is not None and key != "ref"} if shell_map else None),
                                                     reference_ma_side=reference_ma_side)
    ma_tail_record = None
    if shell_map is not None:
        ref_pma = ma_ms_offsets.reference_p_ma(reference_dir, reference_interface_types) if reference_dir else {}
        ma_tail_record = {"Rule": ma_tail.RULE, "Estimator": ma_tail.ESTIMATOR, "RingRadii": context["radial_shells"]["RingRadii"],
                          "GatedOrder": f"p{gated['Order']}", "Orders": {}, "Reference": reference_ma_side}
        for item in mains:
            directory = results / "main" / item["Prefix"] / "reducer"
            per_source = tails_of(directory)
            free = [i for i in per_source if i not in set(zero_trace)]
            strongest = ma_ms_offsets.strongest_sources(ref_pma, [i for i in free if i in ref_pma],
                                                        gates["Gates"]["p_MA"]["StrongestCount"]) if ref_pma else []
            csv_path = comparison_dir / f"ma-sharp-{item['Prefix']}.csv"
            ma_tail.write_csv(csv_path, per_source)
            ma_tail_record["Orders"][f"p{item['Order']}"] = {
                "Prefix": item["Prefix"], "Strongest": strongest, "CSV": str(csv_path),
                "Summary": ma_tail.summary(per_source, free, strongest=strongest),
                "PerSource": {str(i): {key: per_source[i][key] for key in ("Q_MA_raw", "Q_MA_tail", "Q_MA_sharp", "p_MA_raw",
                                                                            "p_MA_sharp", "Alpha", "AlphaSE", "Ring1Factor")}
                              | {"Deficit": per_source[i]["Deficit"]} for i in sorted(per_source)}}
        ma_tail.write_record(comparison_dir / "ma-tail.json", comparison_dir / "ma-tail.md", ma_tail_record,
                             f"Sharp-edge MA extrapolation of {record['Case']}")
        ma_tail_record["Path"] = str(comparison_dir / "ma-tail.json")
    record["MATail"] = ({key: value for key, value in ma_tail_record.items() if key != "Orders"}
                        | {"Orders": {order: {key: value for key, value in block.items() if key != "PerSource"}
                                      for order, block in ma_tail_record["Orders"].items()}}
                        if ma_tail_record else {"Applied": False, "Reason": "no radial MA shells on this run: the MA is raw"})
    context["ma_tail"] = ma_tail_record
    offsets = None
    keys = None
    if reference_dir is not None:
        per_source = {label: summary["PerSource"] for label, summary in comparisons.items()}
        offsets = ma_ms_offsets.write_offsets(per_source, main_label, zero_trace, reference_dir,
                                              comparison_dir / "ma-ms-offsets.md", comparison_dir / "ma-ms-offsets.json",
                                              title=f"Offsets of {record['Case']} vs the reference",
                                              strongest=gates["Gates"]["p_MA"]["StrongestCount"],
                                              interface_types=reference_interface_types)
        keys = key_sources.key_sources(per_source, zero_trace)
        (comparison_dir / "key-sources.md").write_text(key_sources.markdown_report(keys))
        write_json(comparison_dir / "key-sources.json", keys)
    jobs = context["jobs"]
    statuses = context.get("statuses") or {
        job["Name"]: json.loads((results / "main" / Path(job["RemoteDirectory"]).relative_to(record["Remote"]["Case"] + "/main")
                                 / "status.json").read_text()) for job in jobs}
    status = statuses["single"] if len(jobs) == 1 and jobs[0]["Kind"] == "single" else summarize_cost.merge_split_statuses(statuses, jobs)
    cost = summarize_cost.summarize(status, full_sources=len(context["sources"]), nodes=profile["Nodes"])
    worst = record["Split"]["WorstPCGFactor"]
    cost["Jobs"] = {job["Name"]: {"Kind": job["Kind"], "Block": job["Block"], "SourceCount": len(job["Sources"]),
                                  "PBSJobID": (job.get("Submission") or {}).get("Job"),
                                  "SubmittedUTC": (job.get("Submission") or {}).get("UTC"),
                                  "EstimateSecondsWithPreflightAndMargin": job["Estimate"][worst],
                                  "ActualSeconds": statuses[job["Name"]].get("TotalSeconds"),
                                  "ActualOverEstimate": ((statuses[job["Name"]].get("TotalSeconds") or 0.0) / job["Estimate"][worst]
                                                         if job["Estimate"][worst] else None),
                                  "NodeHours": (statuses[job["Name"]].get("TotalSeconds") or 0.0) * profile["Nodes"] / 3600.0,
                                  "StartUTC": statuses[job["Name"]].get("StartUTC"), "EndUTC": statuses[job["Name"]].get("EndUTC"),
                                  "Host": statuses[job["Name"]].get("Host")}
                    for job in jobs}
    submitted = [job["SubmittedAt"] for job in jobs if job.get("SubmittedAt")]
    cost["CriticalPathSeconds"] = (record["FetchedAt"] - min(submitted)) if submitted and record.get("FetchedAt") else None
    cost["CriticalPathRule"] = ("first submission of the coupon's jobs to its fetch (a split coupon: the worker jobs in "
                                "parallel, then the reducer job's queue wait and run)")
    cost["Split"] = {"N": record["Split"]["N"], "Blocks": record["Split"]["Blocks"], "ControlsJob": record["Split"]["ControlsJob"],
                     "Policy": record["JobPolicy"]["Mode"],
                     "CriticalPathEstimateSeconds": record["Split"]["CriticalPathEstimateSeconds"][worst],
                     "NodeSecondsEstimate": record["Split"]["NodeSecondsEstimate"][worst]}

    # Every analysis output goes under the case root; the results directory is read only
    # (it may be a recorded campaign's tree).
    write_json(case_root / "cost-summary.json", cost)
    reference_cost = reference_node_hours(reference, profile)
    main_cost = cost["Stages"].get(main["Prefix"], {})
    gate_record = gate_evaluation.evaluate(
        gates, comparison=comparisons.get(main_label), classes={i: name for i, (name, _) in classes.items()},
        ref_pma=ma_ms_offsets.reference_p_ma(reference_dir, reference_interface_types) if reference_dir else None,
        p_sequence_summary=p_sequence_summary, reference_order=reference_order, gates_sha256=gates_digest,
        interfaces=context["interfaces"], gated_order=gated["Order"])
    write_json(case_root / "qualification.json", gate_record)
    record["Qualification"] = {"Verdict": gate_record["Verdict"], "Reason": gate_record["Reason"],
                               "GatesPassed": gate_record["GatesPassed"], "NotApplicable": gate_record["NotApplicable"],
                               "ReferenceAnchor": gate_record["ReferenceAnchor"],
                               "GatedStage": gated["Prefix"], "GatedOrder": f"p{gated['Order']}", "GatedComparison": main_label,
                               "Path": str(case_root / "qualification.json"),
                               "MAQuantity": (gate_record["Gates"].get("p_MA") or {}).get("Quantity"),
                               "MASharp": ({"Summary": ma_tail_record["Orders"][f"p{gated['Order']}"]["Summary"],
                                            "Reference": {key: value for key, value in (reference_ma_side or {}).items()
                                                          if key != "PerSource"} or None}
                                           if ma_tail_record else None),
                               "ClassStatistics": class_stats.get(main_label),
                               "Offsets": (offsets["Distributions"].get(main_label) if offsets else None),
                               "WeightedPMA": (offsets["Weighted"].get(main_label) if offsets else None),
                               "Informational": {label: {"ClassStatistics": class_stats.get(label),
                                                         "Offsets": (offsets["Distributions"].get(label) if offsets else None),
                                                         "WeightedPMA": (offsets["Weighted"].get(label) if offsets else None)}
                                                 for label in informational}}
    stage_keys = ("H1", "Order", "MeanPCGIterations", "MaxPCGIterations", "SecondsPerPCGIteration", "WorkerWallSeconds",
                  "ReducerWallSeconds", "StageWallSeconds", "NodeHours")
    record["Cost"] = {"MainStage": {key: main_cost.get(key) for key in stage_keys},
                      "MainStages": {item["Prefix"]: {key: cost["Stages"].get(item["Prefix"], {}).get(key) for key in stage_keys}
                                     for item in mains},
                      "JobNodeHours": cost["JobNodeHours"], "JobTotalSeconds": cost["JobTotalSeconds"],
                      "Jobs": cost["Jobs"], "Split": cost["Split"], "CriticalPathSeconds": cost["CriticalPathSeconds"],
                      "CriticalPathRule": cost["CriticalPathRule"],
                      "ReferenceNodeHours": reference_cost,
                      "MainStageOverReference": (main_cost.get("NodeHours") / reference_cost
                                                 if reference_cost and main_cost.get("NodeHours") else None),
                      "Path": str(case_root / "cost-summary.json")}
    record["Status"] = {gate_evaluation.VERDICT_PASSED: STATUS_QUALIFIED, gate_evaluation.VERDICT_PENDING: STATUS_PENDING,
                        gate_evaluation.VERDICT_FAILED: STATUS_FAILED}[gate_record["Verdict"]]
    return gate_record


def reference_node_hours(reference, profile):
    """Worker + reducer wall of the reference campaign case (its status.json) in node-hours."""
    if not reference or not reference.get("Results"):
        return None
    status_path = Path(reference["Results"]).parent / "status.json"
    if not status_path.is_file():
        return None
    status = json.loads(status_path.read_text())
    seconds = 0.0
    nodes = 1
    for stage in status.get("Stages", []):
        seconds += float(stage.get("WallSeconds") or 0.0)
        command = stage.get("Command") or []
        if "-n" in command:
            ranks = int(command[command.index("-n") + 1])
            nodes = max(nodes, -(-ranks // profile["Ranks"]))
    return seconds * nodes / 3600.0 if seconds else None


# The library written for Palace (Version 3): the header of the cases' source
# process-library.json files, every model's files copied under ROOT/models/<slug>/.
PROCESS_LIBRARY_PREFLIGHT_RECORD = "process-library-preflight.json"
LIBRARY_HEADER_EXCLUDED = ("Name", "Models")
MODEL_FILE_NAMES = {"FabricatedMatrix": "fabricated-domain-response-matrix.csv",
                    "FabricatedSurfaceMatrix": "fabricated-surface-response-matrix.csv",
                    "ThinMatrix": "thin-domain-response-matrix.csv",
                    "ThinSurfaceMatrix": "thin-surface-response-matrix.csv",
                    "BasisPoints": "basis-points.csv"}
TRACE_MESH_FILE_NAMES = {"Vertices": "trace-vertices.csv", "Triangles": "trace-triangles.csv"}
THIN_FIELDS = ("ThinMatrix", "ThinSurfaceMatrix")
REDUCER_FILE_NAMES = {"FabricatedMatrix": "domain-response-matrix.csv", "FabricatedSurfaceMatrix": "surface-response-matrix.csv"}
LIBRARY_RULE = ("Version-3 library: the header (MatchingRadius, Fabrication.InterfaceLayers, ...) is the cases' source "
                "process-library.json header (identical across the cases, else the writer stops); every model's matrices, "
                "basis points and trace mesh are copies under models/<slug>/ (paths relative to this file); an entry is "
                "LibraryQualified only with the verdict Passed; PendingQualification entries carry their matrices and "
                "p-sequence controls but no accuracy statement; a model without thin matrices has ThinMatrix null and a "
                "NotLoadable reason (Palace requires the thin response of every model: this file loads only when "
                "NotLoadable is null for every model; process-library-preflight.json is the geometry-only variant)")
PREFLIGHT_RULE = ("geometry-only variant for palace --surface-response-preflight: Palace reads ThinMatrix unconditionally "
                  "as a string, so a NotLoadable model points ThinMatrix / ThinSurfaceMatrix at the copy paths its thin "
                  "matrices will take under models/<slug>/ (absent files: a solve fails closed at file open); never a "
                  "solve input")
NOT_LOADABLE_THIN = "no thin response matrices: {missing} not produced (the source library's {sources} do not exist)"
NOT_LOADABLE_SUPPORT = "no SupportPoints (Palace requires the matching-volume corners of a TraceLiftVersion >= 2 model)"
BASIS_POINTS_FROM_TRACE = ("the producer's basis-points.csv is absent: derived from TraceMesh.Vertices (the rows with basis > 0 in "
                           "basis order, x / y / z in the library frame - the producer's own definition, byte-identical)")


def model_slug(name):
    return re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")


def library_header(library, path):
    """The top-level keys of a source process-library.json apart from Name / Models,
    checked for what Palace requires of a Version-3 library."""
    if library.get("Version") != 3:
        raise ValueError(f"{path}: the source process library is Version {library.get('Version')!r}, not 3")
    if not isinstance(library.get("MatchingRadius"), (int, float)) or library["MatchingRadius"] <= 0.0:
        raise ValueError(f"{path}: the source process library has no positive MatchingRadius")
    if not isinstance((library.get("Fabrication") or {}).get("InterfaceLayers"), dict):
        raise ValueError(f"{path}: the source process library has no Fabrication.InterfaceLayers")
    return {key: value for key, value in library.items() if key not in LIBRARY_HEADER_EXCLUDED}


def same_metadata(first, second):
    """Structural equality with numbers compared to 1e-12 relative."""
    if isinstance(first, dict) and isinstance(second, dict):
        return first.keys() == second.keys() and all(same_metadata(first[key], second[key]) for key in first)
    if isinstance(first, list) and isinstance(second, list):
        return len(first) == len(second) and all(same_metadata(a, b) for a, b in zip(first, second))
    numbers = (int, float)
    if isinstance(first, numbers) and isinstance(second, numbers) and not isinstance(first, bool) and not isinstance(second, bool):
        return math.isclose(float(first), float(second), rel_tol=1.0e-12, abs_tol=0.0)
    return first == second


def consistent_header(headers):
    """One header for the library; the writer stops when two sources disagree on any key."""
    if not headers:
        return None
    first_path, first = headers[0]
    for path, header in headers[1:]:
        for key in sorted(first.keys() | header.keys()):
            if key not in first or key not in header or not same_metadata(first[key], header[key]):
                raise ValueError(f"the source process libraries disagree on {key}: {first_path} has {first.get(key)!r}, "
                                 f"{path} has {header.get(key)!r}")
    return dict(first)


def resolve_model_file(value, directories):
    """An absolute path, or the first of `directories` holding the relative path; None when absent."""
    path = Path(value)
    if path.is_absolute():
        return path if path.is_file() else None
    for directory in directories:
        candidate = Path(directory) / path
        if candidate.is_file():
            return candidate.resolve()
    return None


def write_basis_points_from_trace_vertices(trace_vertices, destination):
    """basis-points.csv as generate_spatial_response writes it: the trace vertices with a
    basis index, in basis order, x / y / z (%.16e, header x,y,z)."""
    with open(trace_vertices, newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if int(row["basis"]) > 0]
    rows.sort(key=lambda row: int(row["basis"]))
    if [int(row["basis"]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError(f"{trace_vertices}: the basis indices are not 1..N")
    with open(destination, "w") as stream:
        stream.write("x,y,z\n")
        for row in rows:
            stream.write(",".join(f"{float(row[name]):.16e}" for name in ("x", "y", "z")) + "\n")


def copy_model_files(model, *, root, directories):
    """Copies the model's matrices, basis points and trace mesh under ROOT/models/<slug>/,
    rewriting the fields to paths relative to ROOT.  A required file that does not exist
    stops the writer (BasisPoints is derived from the trace vertices when the producer's
    file is absent); a missing thin matrix becomes null and is returned as missing."""
    relative_directory = Path("models") / model_slug(model["Name"])
    destination = Path(root) / relative_directory
    destination.mkdir(parents=True, exist_ok=True)
    missing = []
    sources = []
    for field, filename in MODEL_FILE_NAMES.items():
        value = model.get(field)
        source = resolve_model_file(value, directories) if value is not None else None
        if source is None:
            if field in THIN_FIELDS:
                missing.append(field)
                sources.append(str(value))
                model[field] = None
                continue
            if field == "BasisPoints" and "TraceMesh" in model:
                trace_vertices = resolve_model_file(model["TraceMesh"]["Vertices"], directories)
                if trace_vertices is not None:
                    write_basis_points_from_trace_vertices(trace_vertices, destination / filename)
                    model[field] = str(relative_directory / filename)
                    model["BasisPointsRule"] = BASIS_POINTS_FROM_TRACE
                    continue
            raise FileNotFoundError(f"{model['Name']} {field} {value!r} does not exist under {[str(d) for d in directories]}")
        shutil.copy2(source, destination / filename)
        model[field] = str(relative_directory / filename)
    if "TraceMesh" in model:
        trace_mesh = dict(model["TraceMesh"])
        for field, filename in TRACE_MESH_FILE_NAMES.items():
            source = resolve_model_file(trace_mesh[field], directories)
            if source is None:
                raise FileNotFoundError(f"{model['Name']} TraceMesh.{field} {trace_mesh[field]!r} does not exist under "
                                        f"{[str(d) for d in directories]}")
            shutil.copy2(source, destination / filename)
            trace_mesh[field] = str(relative_directory / filename)
        model["TraceMesh"] = trace_mesh
    reasons = []
    if missing:
        reasons.append(NOT_LOADABLE_THIN.format(missing=" / ".join(missing), sources=" / ".join(sources)))
    if not model.get("SupportPoints"):
        reasons.append(NOT_LOADABLE_SUPPORT)
    model["NotLoadable"] = {"Reason": "; ".join(reasons), "Missing": missing} if reasons else None
    return model


def preflight_process_library(library):
    """The geometry-only variant of a written library (PREFLIGHT_RULE)."""
    preflight = json.loads(json.dumps(library))
    preflight["PreflightOnly"] = True
    preflight["PreflightRule"] = PREFLIGHT_RULE
    for model in preflight["Models"]:
        if model.get("NotLoadable") is None:
            continue
        relative_directory = Path("models") / model_slug(model["Name"])
        for field in model["NotLoadable"]["Missing"]:
            model[field] = str(relative_directory / MODEL_FILE_NAMES[field])
    return preflight


def process_library_entries(records, contexts, *, manifest_path, manifest, root, merge_into=None):
    """The process-library entries (LIBRARY_RULE): every qualified / pending coupon's model
    (from the case's own process-library.json) with the fetched response matrices copied
    under ROOT/models/<slug>/ and the qualification bound; LibraryQualified only when
    Passed.  `merge_into` = a previous qualify run's process-library.json (read only): its
    models this run did not qualify are kept ahead of this run's (their files copied from
    the previous root, the previous case root or the model's source directory), a model of
    the same Name is replaced, and MergedFrom records the file, its digest and the names
    kept / replaced.  The header is the source libraries' (consistent_header)."""
    models = []
    headers = []
    for record in records:
        context = contexts.get(record["Case"])
        if context is None or record.get("Qualification") is None:
            continue
        case = manifest_case(manifest, record["Case"])
        directory = source_directory(manifest_path, manifest, case)
        library_path = directory / case["Source"]["Files"]["ProcessLibrary"]["Name"]
        library = json.loads(library_path.read_text())
        headers.append((library_path, library_header(library, library_path)))
        model = dict(library["Models"][0])
        main = next(item for item in context["layout"] if item["Role"] == "main")
        reducer = Path(record.get("MainReducer") or Path(record["Root"]) / "results" / "main" / main["Prefix"] / "reducer")
        for field, filename in REDUCER_FILE_NAMES.items():
            model[field] = str(reducer / filename)
        model["CouponMesh"] = {"Path": record["Mesh"]["Local"], "SHA256": record["Mesh"]["SHA256"]}
        model["Qualification"] = {"Verdict": record["Qualification"]["Verdict"], "Record": record["Qualification"]["Path"],
                                  "ReferenceAnchor": record["Qualification"]["ReferenceAnchor"], "Order": main["Order"]}
        model["LibraryQualified"] = record["Qualification"]["Verdict"] == gate_evaluation.VERDICT_PASSED
        model["SourceProcessLibrary"] = {"Path": str(library_path), "SHA256": sha256(library_path)}
        tail = context.get("ma_tail")
        order_block = (tail or {}).get("Orders", {}).get(f"p{main['Order']}")
        model["MA"] = ({"Rule": tail["Rule"], "Order": main["Order"], "RingRadii": tail["RingRadii"],
                        "Record": str(Path(tail["Path"]).relative_to(root)),
                        "MA_raw": {i: entry["Q_MA_raw"] for i, entry in order_block["PerSource"].items()},
                        "MA_sharp": {i: entry["Q_MA_sharp"] for i, entry in order_block["PerSource"].items()},
                        "p_MA_raw": {i: entry["p_MA_raw"] for i, entry in order_block["PerSource"].items()},
                        "p_MA_sharp": {i: entry["p_MA_sharp"] for i, entry in order_block["PerSource"].items()},
                        "Summary": order_block["Summary"]}
                       if order_block else {"Rule": "raw MA only: the run carries no radial MA shells", "MA_sharp": None})
        models.append(copy_model_files(model, root=root, directories=[directory]))
    merged = None
    if merge_into is not None:
        merge_into = Path(merge_into)
        previous = json.loads(merge_into.read_text())
        names = {model["Name"] for model in models}
        kept = []
        for previous_model in previous["Models"]:
            if previous_model["Name"] in names:
                continue
            model = dict(previous_model)
            source_library = Path(model["SourceProcessLibrary"]["Path"])
            headers.append((source_library, library_header(json.loads(source_library.read_text()), source_library)))
            qualification_record = Path(model["Qualification"]["Record"])
            kept.append(copy_model_files(model, root=root, directories=[merge_into.parent, qualification_record.parents[1],
                                                                        source_library.parent]))
        merged = {"Path": str(merge_into), "SHA256": sha256(merge_into), "Root": previous.get("Root"),
                  "Kept": [model["Name"] for model in kept],
                  "Replaced": [model["Name"] for model in previous["Models"] if model["Name"] in names],
                  "Rule": "models of the previous library this run did not qualify are kept ahead of this run's (their "
                          "files copied into this root); a model of the same Name is replaced by this run's; the previous "
                          "file is not modified"}
        models = kept + models
    header = consistent_header(headers) or {}
    not_loadable = [model["Name"] for model in models if model["NotLoadable"] is not None]
    return {**header, "Name": "coupon-library", "Command": "coupon-library qualify", "Root": str(root),
            "Rule": LIBRARY_RULE,
            "HeaderSources": [{"Path": str(path), "SHA256": sha256(path)} for path, _ in headers],
            "Loadable": {"Palace": not not_loadable, "Models": len(models) - len(not_loadable), "NotLoadable": not_loadable,
                         "Preflight": PROCESS_LIBRARY_PREFLIGHT_RECORD},
            "MergedFrom": merged, "Models": models}


def library_totals(records, *, args, remote, profile, wall_seconds, first_submission, last_fetch, jobs, cap):
    node_hours = sum((record.get("Cost") or {}).get("JobNodeHours") or 0.0 for record in records)
    return {"CouponsAttempted": len(records),
            "CouponsQualified": sum(record["Status"] == STATUS_QUALIFIED for record in records),
            "CouponsPending": sum(record["Status"] == STATUS_PENDING for record in records),
            "CouponsFailed": sum(record["Status"] == STATUS_FAILED for record in records),
            "CouponsPlanned": sum(record["Status"] == STATUS_PLANNED for record in records),
            "CouponsSkipped": sum(record["Status"] == STATUS_SKIPPED for record in records),
            "CouponsUnsupported": sum(record["Status"] == STATUS_UNSUPPORTED for record in records),
            "StoppedCoupons": {record["Case"]: record["StoppedBy"] for record in records if record["StoppedBy"]},
            "NodeHours": node_hours,
            "CriticalPathSeconds": ((last_fetch - first_submission) if first_submission and last_fetch else None),
            "CriticalPathRule": "first submission to the last fetch of this run (the jobs of different coupons overlap up to MaxJobs)",
            "JobWallSeconds": {record["Case"]: (record.get("Cost") or {}).get("JobTotalSeconds") for record in records
                               if record.get("Cost")},
            "WallClockSeconds": wall_seconds, "JobsSubmitted": jobs, "UserJobCap": cap,
            "MaxJobs": args.max_jobs, "JobsRule": "every qsub of this run is counted at submission (a stop after the "
                                                   "submission keeps its job counted); a split coupon counts N worker jobs "
                                                   "and its reducer job",
            "JobPolicy": {"Mode": args.job_policy or "manifest default or frugal", "FixedJobs": args.fixed_jobs,
                          "MaxJobs": args.max_jobs, "Rule": job_split.SPLIT_RULE},
            "ReducerBlockSize": {"CommandLine": args.reducer_block_size,
                                 "PerCase": {record["Case"]: record["ReducerBlockSize"] for record in records
                                             if record.get("ReducerBlockSize")},
                                 "Rule": build_plan.REDUCER_BLOCK_SIZE_RULE},
            "Splits": {record["Case"]: {"N": record["Split"]["N"], "Blocks": record["Split"]["Blocks"],
                                        "Policy": record["JobPolicy"]["Mode"],
                                        "CriticalPathSeconds": (record.get("Cost") or {}).get("CriticalPathSeconds"),
                                        "NodeHours": (record.get("Cost") or {}).get("JobNodeHours")}
                       for record in records if record.get("Split")},
            "DryRun": args.dry_run, "Remote": remote,
            "Orders": [f"p{order}" for order in (args.orders or [])], "OrdersRule": ORDERS_RULE,
            "Controls": [f"p{order}" for order in args.controls],
            "FrozenBinarySHA256": {"CommandLine": args.frozen_binary_sha256,
                                   "PerCase": {record["Case"]: record["FrozenExecutable"] for record in records
                                               if record.get("FrozenExecutable")},
                                   "Default": build_plan.DEFAULT_FROZEN_BINARY_SHA256, "Rule": build_plan.FROZEN_BINARY_RULE},
            "ClusterProfile": profile["Name"]}


def tool_commit():
    """The repository commit of the qualify tooling (None outside a git checkout)."""
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True,
                                       stderr=subprocess.DEVNULL).strip()
    except (subprocess.CalledProcessError, OSError):
        return None


def log_line(message):
    print(message, flush=True)


def run_qualify(args, *, log=log_line):
    build_path = Path(args.build_record).resolve()
    build = json.loads(build_path.read_text())
    manifest_path = Path(build["Library"]["Manifest"]["Path"])
    manifest = json.loads(manifest_path.read_text())
    manifest["Path"] = str(manifest_path)
    if sha256(manifest_path) != build["Library"]["Manifest"]["SHA256"]:
        raise ValueError(f"the manifest {manifest_path} changed since the build record ({build['Library']['Manifest']['SHA256'][:12]}...)")
    profile = json.loads(Path(args.cluster_profile).read_text())
    cost_model = estimate_stages.load_cost_model(args.cost_model)
    gates, gates_digest = gate_evaluation.load_gates(args.gates)
    remote = parse_remote(args.remote)
    if not args.dry_run and remote is None:
        raise ValueError("--remote HOST:ROOT is required unless --dry-run")
    if args.max_jobs > profile["UserJobCap"]:
        raise ValueError(f"--max-jobs {args.max_jobs} exceeds the cluster profile's user job cap {profile['UserJobCap']}")
    root = Path(args.root) if args.root else Path(
        f"/tmp/coupon-library-qualify-{build['Library']['Commit']}-{time.strftime('%Y%m%d-%H%M%S')}")
    root.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.gates, root / GATES_COPY)
    selected = build["Cases"]
    if args.case:
        by_id = {case["Case"]: case for case in build["Cases"]}
        unknown = [case_id for case_id in args.case if case_id not in by_id]
        if unknown:
            raise ValueError(f"cases {unknown} are not in the build record")
        selected = [by_id[case_id] for case_id in args.case]
    start = time.time()
    records, contexts = [], {}
    partial = root / LIBRARY_QUALIFICATION_RECORD

    def stop(record, case_record, exception):
        if record is None:
            record = {"Case": case_record["Case"], "Root": str(root / case_record["Case"])}
        record["Status"] = (STATUS_UNSUPPORTED if exception.record["Kind"] == "ScopeGuard" else
                            STATUS_SKIPPED if exception.record["Kind"] in ("Build", "Manifest", "Reference") else STATUS_FAILED)
        record["StoppedBy"] = exception.record
        log(f"{record['Case']}: {record['Status']} - {exception.record['Kind']}: {exception.record['Message']}")
        return record

    def checkpoint():
        write_json(partial, {"Version": QUALIFICATION_VERSION, "Command": "coupon-library qualify", "Partial": True, "Cases": records})

    # Steps 1-4 for every coupon (a stop records the coupon; the others continue).
    def recorded_plans(case_id):
        case_root = root / case_id
        plans = list(case_root.glob("main/plan.json")) + list(case_root.glob("main/jobs/*/plan.json"))
        return plans_text(plans) if plans else None
    plans_before = {case["Case"]: recorded_plans(case["Case"]) for case in selected} if args.resume else {}
    pending = []
    for case_record in selected:
        record = None
        try:
            record, context = prepare_case(case_record, manifest_path=manifest_path, manifest=manifest, args=args, root=root,
                                           remote=remote, profile=profile, cost_model=cost_model, gates=gates,
                                           gates_digest=gates_digest)
            contexts[record["Case"]] = context
            log(f"{record['Case']}: planned {record['Plan']['StageNames']} ({record['Split']['Decision']})")
            if not args.dry_run:
                pending.append((record, context))
        except CaseStop as exception:
            record = stop(record, case_record, exception)
        records.append(record)
        checkpoint()
    # Step 5: up to --max-jobs of this run's jobs in the queue at once (a coupon's worker
    # jobs, then its reducer job once every worker completed and the archive union is
    # counted; a single-job coupon is one job), a read-only poll of every active job per
    # interval, fetch / verify / analyze each coupon as its last job leaves the queue, the
    # next ready job submitted into the freed slot.  --resume adopts the submissions a
    # previous driver of this root recorded (<case>/submission*.json) instead of
    # uploading and submitting again.
    jobs = 0
    first_submission = last_fetch = None
    active = []      # (record, context, job) submitted and in the queue
    waiting = []     # (record, context, job) not yet submitted

    def drop_coupon(record, context, stopped_job=None):
        """A coupon stopped: its unsubmitted jobs are dropped; its other submitted jobs
        are left to finish on their own (the driver never qdels) and recorded."""
        nonlocal waiting
        waiting = [item for item in waiting if item[0] is not record]
        left = [(job["Name"], job["Submission"]["Job"]) for job in context["jobs"]
                if job is not stopped_job and job.get("Submission") and job.get("FinishedAt") is None]
        if left and record.get("StoppedBy") is not None:
            record["StoppedBy"]["JobsLeftRunning"] = left
            record["StoppedBy"]["JobsLeftRunningRule"] = ("the coupon's other jobs finish on their own (the driver never qdels); "
                                                          "their archives stay under the remote case until deleted by hand")

    for record, context in pending:
        adopted = []
        if args.resume:
            try:
                adopted = resume_submissions(record, context, plans_before.get(record["Case"]))
            except CaseStop as exception:
                stop(record, None, exception)
                checkpoint()
                continue
        for job in context["jobs"]:
            if job in adopted:
                active.append((record, context, job))
                jobs += 1
                first_submission = min(first_submission or job["SubmittedAt"], job["SubmittedAt"])
                log(f"{record['Case']}: resumed {job['Name']} job {job['Submission']['Job']} submitted {job['Submission']['UTC']}")
            else:
                waiting.append((record, context, job))
    checkpoint()

    def ready(record, context, job):
        done = {j["Name"] for j in context["jobs"] if j.get("Status", {}).get("State") == "complete"}
        return all(name in done for name in job["Requires"])

    while waiting or active:
        for item in [item for item in waiting if ready(*item) and len(active) < args.max_jobs]:
            if len(active) >= args.max_jobs:
                break
            record, context, job = item
            waiting.remove(item)
            try:
                if job["Kind"] == "reducer" and record.get("ArchiveUnion") is None:
                    verify_archive_union(record, context, remote=remote, profile=profile)
                    log(f"{record['Case']}: archive union complete "
                        f"{[(k, v['Found']) for k, v in record['ArchiveUnion']['Stages'].items()]}")
                submit_job(record, context, job, remote=remote, profile=profile, log=log)
            except CaseStop as exception:
                stop(record, None, exception)
                # A qsub that went through before the stop is a job of this run (counted).
                if job.get("Submission"):
                    jobs += 1
                drop_coupon(record, context, job)
                checkpoint()
                continue
            jobs += 1
            first_submission = first_submission or job["SubmittedAt"]
            if job["Kind"] == "single":
                record["SubmittedAt"] = job["SubmittedAt"]
            active.append((record, context, job))
            checkpoint()
        if not active:
            if waiting:
                # Nothing active and nothing ready: the waiting jobs' requirements can never
                # complete (their coupons stopped) - drop them.
                for record, context, job in list(waiting):
                    if record.get("StoppedBy") is None:
                        stop(record, None, CaseStop("Scheduler", f"{job['Name']} waits on {job['Requires']} that never completed"))
                    drop_coupon(record, context)
                checkpoint()
            break
        wait(args.monitor_interval)
        still_active = []
        for record, context, job in active:
            if record.get("StoppedBy") is not None:
                continue
            try:
                done = poll_job(record, job, remote=remote, profile=profile, log=log)
                if not done and job["Monitor"]["Polls"] >= args.monitor_polls:
                    raise CaseStop("Monitor", f"{job['Name']} job {job['Submission']['Job']} still {job['Monitor']['LastJobState']} after "
                                              f"{job['Monitor']['Polls']} polls: fetch later with the recorded job id",
                                   Submission=job["Submission"])
                if not done:
                    still_active.append((record, context, job))
                    continue
                job["FinishedAt"] = time.time()
                if job["Kind"] == "worker":
                    complete_worker_job(record, context, job, remote=remote, profile=profile)
                    log(f"{record['Case']}: {job['Name']} complete ({job['Status']['TotalSeconds']:.0f} s)")
                    continue
                results = finish_case(record, context, remote=remote, profile=profile)
                record["FetchedAt"] = time.time()
                last_fetch = record["FetchedAt"]
                analyze_case(record, context, results, gates=gates, gates_digest=gates_digest, profile=profile)
                log(f"{record['Case']}: {record['Qualification']['Verdict']} ({record['Qualification']['Reason']})")
            except CaseStop as exception:
                stop(record, None, exception)
                drop_coupon(record, context, job)
            checkpoint()
        active = [item for item in still_active if item[0].get("StoppedBy") is None]
    record = {"Version": QUALIFICATION_VERSION, "Command": "coupon-library qualify", "Root": str(root),
              "ToolCommit": tool_commit(),
              "BuildRecord": {"Path": str(build_path), "SHA256": sha256(build_path), "Commit": build["Library"]["Commit"]},
              "Gates": {"Path": str(root / GATES_COPY), "SHA256": gates_digest},
              "CostModel": {"Path": str(args.cost_model), "SHA256": sha256(args.cost_model)},
              "ClusterProfile": {"Path": str(args.cluster_profile), "SHA256": sha256(args.cluster_profile)},
              "Cases": records,
              "Library": library_totals(records, args=args, remote=remote, profile=profile, wall_seconds=time.time() - start,
                                        first_submission=first_submission, last_fetch=last_fetch, jobs=jobs,
                                        cap=profile["UserJobCap"])}
    write_json(root / LIBRARY_QUALIFICATION_RECORD, record)
    library = process_library_entries(records, contexts, manifest_path=manifest_path, manifest=manifest, root=root,
                                      merge_into=args.merge_into)
    write_json(root / PROCESS_LIBRARY_RECORD, library)
    write_json(root / PROCESS_LIBRARY_PREFLIGHT_RECORD, preflight_process_library(library))
    return record


def add_arguments(parser):
    parser.add_argument("--build-record", type=Path, required=True, help="library-build.json of coupon-library build")
    parser.add_argument("--reference", type=lambda text: None if text.lower() == "none" else Path(text), required=True,
                        help="graded_v2 reference campaign directory (inputs-<key>/, case-<key>-fabricated/) or 'none' "
                             "(the coupon runs on its own inputs; verdict PendingQualification / Failed)")
    parser.add_argument("--remote", help="HOST:ROOT of the frozen executable and MPI wrapper (required unless --dry-run)")
    parser.add_argument("--orders", type=parse_orders, default=[], help="additional main orders, full stages (the recipe's "
                                                                       "PhysicsRun order is always the first main order)")
    parser.add_argument("--controls", type=parse_orders, default=[3, 5], help="control orders (default p3,p5)")
    parser.add_argument("--control-count", type=int, default=DEFAULT_CONTROL_COUNT, help="controls chosen by class (default 8)")
    parser.add_argument("--control-source", type=int, action="append", help="explicit control source (repeatable; overrides the class choice)")
    parser.add_argument("--max-jobs", type=int, default=40, help="at most this many of the run's jobs queued / running at once "
                                                                  "(a coupon's split is bounded by it too; every qsub is counted "
                                                                  "against the cluster profile's user cap)")
    parser.add_argument("--job-policy", choices=job_split.MODES, default=None,
                        help="per-coupon source split policy (decision 61b): speed = the N <= --max-jobs minimizing the "
                             "estimated critical path, frugal = the fewest jobs that fit the walltime, fixed = --fixed-jobs; "
                             "default: the manifest's ProductionRecipe.PhysicsRun.JobPolicy, else frugal")
    parser.add_argument("--fixed-jobs", type=int, default=None, help="N of --job-policy fixed")
    parser.add_argument("--reducer-block-size", type=int, default=None,
                        help="PALACE_RESPONSE_BLOCK_SIZE of every reducer stage (decision 62(1)); default: the manifest's "
                             f"ProductionRecipe.PhysicsRun.ReducerBlockSize, else {build_plan.DEFAULT_REDUCER_BLOCK_SIZE}")
    parser.add_argument("--frozen-binary-sha256", default=None,
                        help="SHA-256 of the frozen Palace executable under ROOT (decision 63); default: the manifest's "
                             f"ProductionRecipe.PhysicsRun.FrozenExecutable, else {build_plan.DEFAULT_FROZEN_BINARY_SHA256[:12]}...")
    parser.add_argument("--stage-prefix", help="stage name prefix (default: the case id)")
    parser.add_argument("--case", action="append", default=None, help="case id (repeatable; default: every case of the build record)")
    parser.add_argument("--root", type=Path, help="local run root (default /tmp/coupon-library-qualify-<commit>-<ts>)")
    parser.add_argument("--reference-edge-size-nm", type=float, default=None,
                        help="edge size (nm) of a reference whose ring / edge sizing is not recorded (the graded_v2 tet "
                             "references: 0.5): its sharp-edge deficit is modelled per source from the run's by the eps^(1/3) "
                             "law and recorded 'reference unextrapolated' (ma_tail.py); omitted = the raw reference stands in")
    parser.add_argument("--merge-into", type=Path, default=None,
                        help="a previous qualify run's process-library.json (read only): ROOT/process-library.json holds "
                             "its models this run did not qualify plus this run's (the same Name replaced), MergedFrom recorded")
    parser.add_argument("--dry-run", action="store_true", help="write plans / configs / estimates / gates; contact nothing")
    parser.add_argument("--resume", action="store_true", help="adopt the job ids a previous driver of this --root recorded "
                                                               "(<case>/submission.json, plan byte-identical) instead of "
                                                               "submitting again; monitor / fetch / analyze from there")
    parser.add_argument("--monitor-interval", type=int, default=DEFAULT_MONITOR_INTERVAL, help="seconds between read-only polls")
    parser.add_argument("--monitor-polls", type=int, default=DEFAULT_MONITOR_POLLS)
    parser.add_argument("--cluster-profile", type=Path, default=HERE / "cluster-profile.json")
    parser.add_argument("--cost-model", type=Path, default=estimate_stages.COST_MODEL)
    parser.add_argument("--gates", type=Path, default=gate_evaluation.GATES_FILE)


def run_from_args(args):
    record = run_qualify(args)
    totals = record["Library"]
    for case in record["Cases"]:
        stopped = case.get("StoppedBy")
        verdict = (case.get("Qualification") or {}).get("Verdict")
        print(f"{case['Case']}: {case['Status']}" + (f" verdict {verdict}" if verdict else "")
              + (f" stopped-by {stopped['Kind']}: {stopped['Message']}" if stopped else ""), flush=True)
    print(f"LIBRARY attempted {totals['CouponsAttempted']} qualified {totals['CouponsQualified']} pending {totals['CouponsPending']} "
          f"failed {totals['CouponsFailed']} planned {totals['CouponsPlanned']} skipped {totals['CouponsSkipped']} "
          f"unsupported {totals['CouponsUnsupported']} "
          f"node-h {totals['NodeHours']:.3f} jobs {totals['JobsSubmitted']}/{totals['MaxJobs']} (cap {totals['UserJobCap']}); "
          f"record {Path(record['Root']) / LIBRARY_QUALIFICATION_RECORD}")
    if args.dry_run:
        return 0 if totals["CouponsPlanned"] == totals["CouponsAttempted"] else 1
    return 0 if totals["CouponsQualified"] == totals["CouponsAttempted"] else 1


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    return run_from_args(parser.parse_args(argv))


if __name__ == "__main__":
    sys.exit(main())
