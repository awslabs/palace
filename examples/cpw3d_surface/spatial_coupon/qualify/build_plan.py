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
REDUCER_ENVIRONMENT = {"PALACE_RESPONSE_REDUCE_ONLY": "1", "PALACE_RESPONSE_BLOCK_SIZE": "6"}
PLAN_VERSION = 2


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


def stage(remote_case_root, prefix, kind, cap, minimum):
    environment = {"PALACE_RESPONSE_ARCHIVE_DIR": f"{remote_case_root}/main/{prefix}/archive"}
    environment.update(WORKER_ENVIRONMENT if kind == "worker" else REDUCER_ENVIRONMENT)
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
               profile, binary, binary_sha256, mpiexec, purpose, factors):
    """The plan dict.  `stage_layout` = [{"Prefix", "Kind": "response" | "local-edge",
    "Order", "Sources", "EstimateKey"}, ...] in run order; `config_digests` = stage
    prefix -> {file name -> sha256}; `trace_pins` = remote trace path -> sha256."""
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
                stages.append(stage(remote_case_root, item["Prefix"], kind, cap, minimum))
        else:
            cap, minimum = stage_caps(est, "local-edge", factors, deadline)
            stages.append(ordinary_stage(remote_case_root, item["Prefix"], cap, minimum))
    return {"Version": PLAN_VERSION, "Case": case_id, "Purpose": purpose,
            "Ranks": profile["Ranks"], "DeadlineSeconds": profile["DeadlineSeconds"],
            "DeadlineMarginSeconds": profile["DeadlineMarginSeconds"],
            "MinimumMemAvailableBytes": profile["MinimumMemAvailableBytes"],
            "Binary": binary, "BinarySHA256": binary_sha256, "MPIExec": mpiexec,
            "PinnedSHA256": pinned, "Stages": stages,
            "MeshSHA256": mesh["SHA256"], "MeshRemote": mesh["Remote"], "MeshLocal": mesh["Local"],
            "CapRule": (f"CapSeconds = {CAP_FACTOR:g} x the stage estimate at {max(factors, key=float)}x the measured PCG "
                        f"counts rounded up to {CAP_ROUNDING_SECONDS} s and bounded by the deadline; MinimumSeconds = "
                        f"the estimate at {min(factors, key=float)}x rounded up"),
            "EstimateSummary": {"JobSecondsEstimateByPCGFactor": estimate["JobSecondsEstimateByPCGFactor"],
                                "JobSecondsEstimateWithPreflightAndMargin": estimate["JobSecondsEstimateWithPreflightAndMargin"],
                                "MaxPalacePeakGBEstimate": estimate["MaxPalacePeakGBEstimate"],
                                "FitsOneJob": estimate["FitsOneJob"], "Decision": estimate["Decision"]}}


def render_job_script(*, profile, remote_root, remote_case_root, runner, job_name, walltime_seconds):
    """The PBS job script of one coupon (one exclusive node; the runner reads the plan)."""
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
             f"#PBS -l instance_type={profile['InstanceType']}",
             f"#PBS -l {profile['ExtraResources']}",
             f"#PBS -l walltime={hours:02d}:{minutes:02d}:{seconds:02d}",
             "#PBS -j oe",
             f"#PBS -o {remote_case_root}/main/pbs.log",
             "set -euo pipefail",
             "unset PYTHONOPTIMIZE",
             "source /etc/profile.d/modules.sh",
             f"module load {modules}",
             "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 ARMPL_NUM_THREADS=1 "
             "VECLIB_MAXIMUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1",
             'for name in ${!PALACE_RESPONSE_@}; do unset "$name"; done',
             f"R={remote_root}",
             f"D={remote_case_root}/main",
             "trap 'code=$?; printf \"{\\\"ExitCode\\\":%d,\\\"JobID\\\":\\\"%s\\\",\\\"UTC\\\":\\\"%s\\\"}\\n\" \"$code\" "
             "\"${PBS_JOBID:-}\" \"$(date -u +%FT%TZ)\" > \"$D/pbs-status.json\"' EXIT",
             'mkdir "$D/tmp"; export TMPDIR="$D/tmp"',
             'date -u +%FT%TZ; hostname; cat "$PBS_NODEFILE"',
             f'python3 "{runner}" "$D/plan.json"']
    return "\n".join(lines) + "\n"
