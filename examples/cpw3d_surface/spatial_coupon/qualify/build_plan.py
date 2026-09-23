#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The bounded-runner plan of one coupon (four-edge-physics-11 / gallery-physics-06b
build_plan.py, parametrized): pinned SHA-256 of the mesh, every stage config and every
trace file; the stages in order - the main-order worker / reducer on every source,
then for each control order the worker / reducer on the control sources, then the
ordinary-path local-edge stage (main order, controls, SaveLocalEdgeEnergy); the frozen
executable and MPI wrapper the runner must find under the remote root.

Stage caps are estimate-derived (the rule of the gallery runs): CapSeconds = 2 x the
stage's estimate at the largest PCG factor, rounded up to CAP_ROUNDING_SECONDS and
bounded by the deadline; MinimumSeconds = the estimate at the measured PCG counts,
rounded up (a stage is not started with less time than its 1x estimate).  The
reducer's cap covers the reducer estimate alone (the worker has its own).
"""
import math

CAP_ROUNDING_SECONDS = 300
CAP_FACTOR = 2.0
WORKER_ENVIRONMENT = {"PALACE_RESPONSE_ARCHIVE_ONLY": "1", "PALACE_RESPONSE_SOURCE_TIMING": "1"}
# The reducer reads, differentiates and evaluates every archived source field once per
# block pair of PALACE_RESPONSE_BLOCK_SIZE sources (N x ceil(N / b) field evaluations for
# N sources): b = 6 (physics-01..-13, gallery, the decision-58 library run) evaluated
# every source 32 times at 191 sources.  Decision 62(1) raises the recorded default to
# 48; the memory rationale is REDUCER_BLOCK_SIZE_RULE (the plan's MinimumMemAvailableBytes
# admission guard is unchanged).
DEFAULT_REDUCER_BLOCK_SIZE = 48
REDUCER_BLOCK_SIZE_RULE = ("PALACE_RESPONSE_BLOCK_SIZE = b: the reducer keeps 2b archived source fields resident per block pair "
                           "(each ~14 MB per rank at p4: V, E and their grid functions with ghosts, ~2.7 GB node-wide) and "
                           "evaluates every source ceil(N / b) times; b = 48 -> 96 resident fields, estimated 350-450 GB of "
                           "the r8g.48xlarge's 1,485 GiB (b = 6 measured 194 GB); the per-entry sums are identical in "
                           "identical sample order, so the matrices are independent of b (decision 62(1) acceptance: "
                           "bit-identical CSVs at b = 6 and b = 48). With the streaming one-pass Gram executable "
                           "(decision 62(4), SHA-256 170439c4...) the reducer evaluates every source once whatever b and "
                           "keeps N x (4 Q_local + L_local) x 8 bytes resident, so b = N (a single block) is permitted "
                           "(decision 63 rule; PBS 46718: p4 78 sources at b = 48 peak 33.3 GB vs 91.2 GB with b28); the "
                           "default stays 48")
PLAN_VERSION = 3
# The frozen Palace executable every qualify stage runs (decision 63): the streaming
# one-pass Gram build 170439c4... (decision 62(4), PBS 46717 build of the b1e7e9e9d
# source freeze b3103728...) replaced the b28f089a... executable of every earlier
# campaign (physics-01..-13, gallery, the decision-58 library run, the decision-61
# acceptances); --frozen-binary-sha256 overrides, the manifest's
# ProductionRecipe.PhysicsRun.FrozenExecutable records the library default.
DEFAULT_FROZEN_BINARY_SHA256 = "170439c4a9fc5d5ce329310812055be5fb83a4a7f288024b57b3b83551cbe70b"
PREVIOUS_FROZEN_BINARY_SHA256 = "b28f089ae12c25863493566b2b8ca11af2c8ffb0e273e7aa67a2b42046eacf27"
FROZEN_BINARY_RULE = ("the frozen palace-archive-estimate-<sha256>.bin under the remote root runs every stage; the plan pins "
                      "the digest and run_stages verifies it in the preflight; 170439c4... (streaming one-pass Gram, decision "
                      "62(4)) reproduces b28f089a... at p3 bit for bit and at p4 to 8.4e-13 per entry (PBS 46718); "
                      "--frozen-binary-sha256 overrides the default and the origin is recorded per coupon")


GIB = 1024 ** 3


def largest_node_gib(profile):
    """The OS-visible memory of the profile's largest instance (the fallback): the node the
    whole-coupon fit checks are made against."""
    return max(instance["NodeGiB"] for instance in profile["Instances"])


def select_instance(profile, peak_gb, gb_per_gib):
    """The instance a job runs on (profile InstanceRule): the first of Instances whose
    MemoryFitFraction x MemoryGiB holds the job's largest estimated Palace peak, else the
    last (largest) one with Fits false (the caller's whole-coupon fit check stops the case);
    the record carries the admission guard MinimumMemAvailableBytes = the fraction of the
    chosen instance's memory."""
    fraction = float(profile["MemoryFitFraction"])
    peak_gib = float(peak_gb) / float(gb_per_gib)
    instances = profile["Instances"]
    chosen = next((item for item in instances if peak_gib <= fraction * item["MemoryGiB"]), instances[-1])
    return {"Type": chosen["Type"], "vCPUs": chosen["vCPUs"], "MemoryGiB": chosen["MemoryGiB"], "NodeGiB": chosen["NodeGiB"],
            "EstimatedPalacePeakGB": float(peak_gb), "EstimatedPalacePeakGiB": peak_gib, "MemoryFitFraction": fraction,
            "Fits": peak_gib <= fraction * chosen["MemoryGiB"],
            "MinimumMemAvailableBytes": int(fraction * chosen["MemoryGiB"] * GIB),
            "Candidates": [item["Type"] for item in instances], "Rule": profile["InstanceRule"]}


def reducer_environment(block_size):
    block_size = int(block_size)
    if block_size < 1:
        raise ValueError(f"the reducer block size must be an integer >= 1, not {block_size}")
    return {"PALACE_RESPONSE_REDUCE_ONLY": "1", "PALACE_RESPONSE_BLOCK_SIZE": str(block_size)}


REDUCER_ENVIRONMENT = reducer_environment(DEFAULT_REDUCER_BLOCK_SIZE)


def round_up(seconds, rounding=CAP_ROUNDING_SECONDS):
    return int(math.ceil(max(seconds, 0.0) / rounding) * rounding)


def stage_caps(estimate_stage, kind, factors, deadline):
    """(CapSeconds, MinimumSeconds) of a worker / reducer / local-edge stage."""
    lowest, highest = min(factors, key=float), max(factors, key=float)
    by_factor = estimate_stage["ByPCGFactor"]
    if kind == "worker":
        low, high = by_factor[lowest]["WorkerSecondsEstimate"], by_factor[highest]["WorkerSecondsEstimate"]
    elif kind == "reducer":
        low = high = estimate_stage["ReducerSecondsEstimate"]
    else:
        low, high = by_factor[lowest]["StageSecondsEstimate"], by_factor[highest]["StageSecondsEstimate"]
    cap = min(round_up(CAP_FACTOR * high), deadline)
    minimum = min(round_up(low), cap)
    return cap, minimum


def stage(remote_case_root, prefix, kind, cap, minimum, reducer_block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    environment = {"PALACE_RESPONSE_ARCHIVE_DIR": f"{remote_case_root}/main/{prefix}/archive"}
    environment.update(WORKER_ENVIRONMENT if kind == "worker" else reducer_environment(reducer_block_size))
    return {"Name": f"{prefix}-{kind}", "Config": f"{remote_case_root}/main/{prefix}/{kind}.json",
            "Environment": environment, "CapSeconds": cap, "MinimumSeconds": minimum,
            "Requires": [] if kind == "worker" else [f"{prefix}-worker"]}


def ordinary_stage(remote_case_root, prefix, cap, minimum):
    # No PALACE_RESPONSE_* variables: run_stages.py strips them from the inherited
    # environment, so Palace runs its ordinary per-source postprocessing (surface-Q*.csv
    # including the per-segment surface-Q-edge-local.csv) and the in-memory matrix.
    return {"Name": prefix, "Config": f"{remote_case_root}/main/{prefix}/config.json", "Environment": {},
            "CapSeconds": cap, "MinimumSeconds": minimum, "Requires": []}


def build_plan(*, case_id, remote_case_root, mesh, stage_layout, estimate, config_digests, trace_pins,
               profile, binary, binary_sha256, mpiexec, purpose, factors, reducer_block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    """The plan dict.  `stage_layout` = [{"Prefix", "Kind": "response" | "local-edge",
    "Order", "Sources", "EstimateKey"}, ...] in run order; `config_digests` = stage
    prefix -> {file name -> sha256}; `trace_pins` = remote trace path -> sha256;
    `reducer_block_size` = PALACE_RESPONSE_BLOCK_SIZE of every reducer stage."""
    pinned = {mesh["Remote"]: mesh["SHA256"]}
    for prefix, digests in config_digests.items():
        for name, digest in digests.items():
            pinned[f"{remote_case_root}/main/{prefix}/{name}"] = digest
    pinned.update(trace_pins)
    deadline = profile["DeadlineSeconds"] - profile["DeadlineMarginSeconds"]
    stages = []
    for item in stage_layout:
        est = estimate["Stages"][item["EstimateKey"]]
        if item["Kind"] == "response":
            for kind in ("worker", "reducer"):
                cap, minimum = stage_caps(est, kind, factors, deadline)
                stages.append(stage(remote_case_root, item["Prefix"], kind, cap, minimum, reducer_block_size))
        else:
            cap, minimum = stage_caps(est, "local-edge", factors, deadline)
            stages.append(ordinary_stage(remote_case_root, item["Prefix"], cap, minimum))
    instance = select_instance(profile, estimate["MaxPalacePeakGBEstimate"], estimate["CostModel"]["PalaceGBPerGiB"])
    return {"Version": PLAN_VERSION, "Case": case_id, "Purpose": purpose,
            "Ranks": profile["Ranks"], "DeadlineSeconds": profile["DeadlineSeconds"],
            "DeadlineMarginSeconds": profile["DeadlineMarginSeconds"],
            "Instance": instance, "MinimumMemAvailableBytes": instance["MinimumMemAvailableBytes"],
            "Binary": binary, "BinarySHA256": binary_sha256, "MPIExec": mpiexec,
            "PinnedSHA256": pinned, "Stages": stages,
            "ReducerBlockSize": int(reducer_block_size), "ReducerBlockSizeRule": REDUCER_BLOCK_SIZE_RULE,
            "MeshSHA256": mesh["SHA256"], "MeshRemote": mesh["Remote"], "MeshLocal": mesh["Local"],
            "CapRule": (f"CapSeconds = {CAP_FACTOR:g} x the stage estimate at {max(factors, key=float)}x the measured PCG "
                        f"counts rounded up to {CAP_ROUNDING_SECONDS} s and bounded by the deadline; MinimumSeconds = "
                        f"the estimate at {min(factors, key=float)}x rounded up"),
            "EstimateSummary": {"JobSecondsEstimateByPCGFactor": estimate["JobSecondsEstimateByPCGFactor"],
                                "JobSecondsEstimateWithPreflightAndMargin": estimate["JobSecondsEstimateWithPreflightAndMargin"],
                                "MaxPalacePeakGBEstimate": estimate["MaxPalacePeakGBEstimate"],
                                "FitsOneJob": estimate["FitsOneJob"], "Decision": estimate["Decision"]}}


def block_stage(remote_case_root, prefix, block, cap, minimum):
    """The worker stage of source block `block` of a split main stage (decision 61b): its
    own config and Problem.Output, the stage's ONE archive directory (the union)."""
    environment = {"PALACE_RESPONSE_ARCHIVE_DIR": f"{remote_case_root}/main/{prefix}/archive"}
    environment.update(WORKER_ENVIRONMENT)
    return {"Name": f"{prefix}-worker-block{block}", "Config": f"{remote_case_root}/main/{prefix}/worker-block{block}.json",
            "Environment": environment, "CapSeconds": cap, "MinimumSeconds": minimum, "Requires": []}


def block_estimate(estimate_stage, block_size, factors):
    """The estimate view of one source block of a main stage (stage_caps' shape)."""
    return {"ByPCGFactor": {factor: {"WorkerSecondsEstimate": (estimate_stage["WorkerNonSourceSecondsEstimate"]
                                                               + block_size * estimate_stage["ByPCGFactor"][factor]["PerSourceSecondsEstimate"])}
                            for factor in factors},
            "ReducerSecondsEstimate": estimate_stage["ReducerSecondsEstimate"]}


def build_job_plan(*, case_id, job_name, remote_case_root, mesh, stage_layout, split_job, estimate, config_digests,
                   trace_pins, profile, binary, binary_sha256, mpiexec, purpose, factors,
                   reducer_block_size=DEFAULT_REDUCER_BLOCK_SIZE):
    """The plan of one job of a split coupon (decision 61b).  `split_job` = a job record
    of job_split.plan_split (Kind worker: the main-stage worker of its block - job 1 also
    the control stages and the local-edge stage; Kind reducer: the main-stage reducers,
    Requires empty - the driver submits it after every worker job completed and the
    archive union is counted).  Pins: the mesh, the configs of the stages this job runs
    and every trace (the reducer config names them all)."""
    kind, block = split_job["Kind"], split_job["Block"]
    stage_names = []
    deadline = profile["DeadlineSeconds"] - profile["DeadlineMarginSeconds"]
    stages = []
    peak_gb = 0.0   # the largest estimated Palace peak of the stages THIS job runs (its instance)
    pinned = {mesh["Remote"]: mesh["SHA256"]}
    for item in stage_layout:
        est = estimate["Stages"][item["EstimateKey"]]
        if item["Role"] == "main":
            if kind == "worker" and split_job["Sources"]:
                cap, minimum = stage_caps(block_estimate(est, len(split_job["Sources"]), factors), "worker", factors, deadline)
                stages.append(block_stage(remote_case_root, item["Prefix"], block, cap, minimum))
                peak_gb = max(peak_gb, est["WorkerPalacePeakGBEstimate"])
                pinned[f"{remote_case_root}/main/{item['Prefix']}/worker-block{block}.json"] = config_digests[item["Prefix"]][f"worker-block{block}.json"]
            elif kind == "reducer":
                cap, minimum = stage_caps(est, "reducer", factors, deadline)
                reducer = stage(remote_case_root, item["Prefix"], "reducer", cap, minimum, reducer_block_size)
                reducer["Requires"] = []
                stages.append(reducer)
                peak_gb = max(peak_gb, est["ReducerPalacePeakGBEstimate"])
                pinned[f"{remote_case_root}/main/{item['Prefix']}/reducer.json"] = config_digests[item["Prefix"]]["reducer.json"]
        elif kind == "worker" and block == 1:
            for name, digest in config_digests[item["Prefix"]].items():
                pinned[f"{remote_case_root}/main/{item['Prefix']}/{name}"] = digest
            if item["Kind"] == "response":
                for stage_kind in ("worker", "reducer"):
                    cap, minimum = stage_caps(est, stage_kind, factors, deadline)
                    stages.append(stage(remote_case_root, item["Prefix"], stage_kind, cap, minimum, reducer_block_size))
                peak_gb = max(peak_gb, est["WorkerPalacePeakGBEstimate"], est["ReducerPalacePeakGBEstimate"])
            else:
                cap, minimum = stage_caps(est, "local-edge", factors, deadline)
                stages.append(ordinary_stage(remote_case_root, item["Prefix"], cap, minimum))
                peak_gb = max(peak_gb, est["PalacePeakGBEstimate"])
    pinned.update(trace_pins)
    stage_names = [item["Name"] for item in stages]
    instance = select_instance(profile, peak_gb, estimate["CostModel"]["PalaceGBPerGiB"])
    return {"Version": PLAN_VERSION, "Case": case_id, "Job": job_name, "JobKind": kind, "Block": block,
            "BlockSources": split_job["Sources"] if kind == "worker" else None, "Purpose": purpose,
            "Ranks": profile["Ranks"], "DeadlineSeconds": profile["DeadlineSeconds"],
            "DeadlineMarginSeconds": profile["DeadlineMarginSeconds"],
            "Instance": instance, "MinimumMemAvailableBytes": instance["MinimumMemAvailableBytes"],
            "Binary": binary, "BinarySHA256": binary_sha256, "MPIExec": mpiexec,
            "PinnedSHA256": pinned, "Stages": stages, "StageNames": stage_names,
            "ReducerBlockSize": int(reducer_block_size), "ReducerBlockSizeRule": REDUCER_BLOCK_SIZE_RULE,
            "MeshSHA256": mesh["SHA256"], "MeshRemote": mesh["Remote"], "MeshLocal": mesh["Local"],
            "CapRule": (f"CapSeconds = {CAP_FACTOR:g} x the stage estimate at {max(factors, key=float)}x the measured PCG "
                        f"counts rounded up to {CAP_ROUNDING_SECONDS} s and bounded by the deadline; MinimumSeconds = "
                        f"the estimate at {min(factors, key=float)}x rounded up; a block worker is estimated at its "
                        f"block's source count"),
            "EstimateSummary": {"JobSecondsEstimateWithPreflightAndMargin": split_job["SecondsEstimateWithPreflightAndMargin"],
                                "Fits": split_job["Fits"]}}


def render_job_script(*, profile, remote_root, remote_case_root, runner, job_name, walltime_seconds, instance_type,
                      job_directory=None):
    """The PBS job script of one coupon job (one exclusive node of `instance_type`, the
    plan's chosen instance; the runner reads the plan of `job_directory`, default
    <case>/main - the single job of a coupon)."""
    if instance_type not in {item["Type"] for item in profile["Instances"]}:
        raise ValueError(f"{instance_type!r} is not an instance of the cluster profile")
    job_directory = job_directory or f"{remote_case_root}/main"
    hours, rest = divmod(int(walltime_seconds), 3600)
    minutes, seconds = divmod(rest, 60)
    modules = " ".join(profile["Modules"])
    lines = ["#!/bin/bash",
             "# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.",
             "# SPDX-License-Identifier: Apache-2.0",
             f"#PBS -N {job_name}",
             f"#PBS -q {profile['Queue']}",
             f"#PBS -P {profile['Project']}",
             "#PBS -r n",
             f"#PBS -l {profile['SelectResources']}",
             f"#PBS -l {profile['PlaceResources']}",
             f"#PBS -l instance_type={instance_type}",
             f"#PBS -l {profile['ExtraResources']}",
             f"#PBS -l walltime={hours:02d}:{minutes:02d}:{seconds:02d}",
             "#PBS -j oe",
             f"#PBS -o {job_directory}/pbs.log",
             "set -euo pipefail",
             "unset PYTHONOPTIMIZE",
             "source /etc/profile.d/modules.sh",
             f"module load {modules}",
             "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ARMPL_NUM_THREADS=1 "
             "VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1",
             'for name in ${!PALACE_RESPONSE_@}; do unset "$name"; done',
             f"R={remote_root}",
             f"D={job_directory}",
             "trap 'code=$?; printf \"{\\\"ExitCode\\\":%d,\\\"JobID\\\":\\\"%s\\\",\\\"UTC\\\":\\\"%s\\\"}\\n\" \"$code\" "
             "\"${PBS_JOBID:-}\" \"$(date -u +%FT%TZ)\" > \"$D/pbs-status.json\"' EXIT",
             'mkdir "$D/tmp"; export TMPDIR="$D/tmp"',
             'date -u +%FT%TZ; hostname; cat "$PBS_NODEFILE"',
             f'python3 "{runner}" "$D/plan.json"']
    return "\n".join(lines) + "\n"
