#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The decision-62(4) acceptance job plan (perf-gram-01): on the perf-reducer-01 archives of
two-edge-8dd4bc70f183 (PBS 46685: the p3 control, 8 sources, and the p4 stage, 78 sources,
whose b28 reductions at block size 48 are recorded) the NEW executable reduces the same
archives at block size 48; then a fresh p4 worker archive is reduced by the OLD (b28) and
the NEW executable in the same job (the realistic timing pair). Every stage names its
executable (run_stages.py stage Binary / BinarySHA256).

usage: plan_gram_acceptance.py NEW_SHA256 OUT_DIR
"""
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import build_plan  # noqa: E402

ROOT = "/data/home/simlap/coupon_accuracy_assessment_20260913"
OLD_SHA = "b28f089ae12c25863493566b2b8ca11af2c8ffb0e273e7aa67a2b42046eacf27"
PREVIOUS = f"{ROOT}/perf-reducer-01/two-edge-8dd4bc70f183/main"
RUN = f"{ROOT}/perf-gram-01"
CASE = f"{RUN}/two-edge-8dd4bc70f183"
JOB = f"{CASE}/main"


def main():
    new_sha, out = sys.argv[1], Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    profile = json.loads((Path(__file__).resolve().parents[2] / "cluster-profile.json").read_text())
    recorded = json.load(open("/tmp/perf-reducer-01/fetch/recorded-plan-worker-1.json"))
    pinned = {k: v for k, v in recorded["PinnedSHA256"].items() if k.endswith(".msh") or "/inputs/traces/" in k}
    binaries = {"old": (f"{ROOT}/palace-archive-estimate-{OLD_SHA}.bin", OLD_SHA),
                "new": (f"{ROOT}/palace-archive-estimate-{new_sha}.bin", new_sha)}
    stages = []

    def config(src, name, output):
        d = json.load(open(src)); d["Problem"]["Output"] = output
        text = json.dumps(d, indent=2) + "\n"
        (out / name).write_text(text)
        pinned[f"{JOB}/{name}"] = hashlib.sha256(text.encode()).hexdigest()
        return f"{JOB}/{name}"

    def reducer(name, src, archive, which, cap):
        path, sha = binaries[which]
        stages.append({"Name": name, "Config": config(src, f"{name}.json", f"{JOB}/{name}"),
                       "Environment": {"PALACE_RESPONSE_ARCHIVE_DIR": archive, **build_plan.reducer_environment(48)},
                       "CapSeconds": cap, "MinimumSeconds": 60, "Requires": [], "Binary": path, "BinarySHA256": sha})

    # The recorded archives (PBS 46685) reduced by the new executable.
    reducer("p3-control-reducer-new", "/tmp/perf-reducer-01/fetch/p3-reducer.json", f"{PREVIOUS}/p3-control/archive", "new", 600)
    reducer("p4-reducer-new", "/tmp/perf-reducer-01/fetch/p4-reducer.json", f"{PREVIOUS}/p4/archive", "new", 900)
    # A fresh p4 worker archive reduced by both executables.
    archive = f"{JOB}/p4-fresh/archive"
    stages.append({"Name": "p4-fresh-worker", "Config": config("/tmp/perf-reducer-01/fetch/p4-reducer.json", "p4-fresh-worker.json", f"{JOB}/p4-fresh/worker"),
                   "Environment": {"PALACE_RESPONSE_ARCHIVE_DIR": archive, **build_plan.WORKER_ENVIRONMENT},
                   "CapSeconds": 1800, "MinimumSeconds": 60, "Requires": [], "Binary": binaries["old"][0], "BinarySHA256": OLD_SHA})
    for which in ("old", "new"):
        reducer(f"p4-fresh-reducer-{which}", "/tmp/perf-reducer-01/fetch/p4-reducer.json", archive, which, 900)
        stages[-1]["Requires"] = ["p4-fresh-worker"]
    plan = {"Version": build_plan.PLAN_VERSION, "Case": "two-edge-8dd4bc70f183", "Job": "perf-gram-01", "JobKind": "acceptance",
            "Purpose": "decision 62(4) acceptance: the NEW streaming-Gram executable reduces the perf-reducer-01 archives (PBS 46685; "
                       "the b28 reductions at block size 48 are recorded there) and, on a fresh p4 worker archive, the OLD (b28) and "
                       "NEW executables reduce the same archive at block size 48 (roundoff-identical matrices <= 1e-12 relative, "
                       "reducer wall and peak memory old vs new)",
            "Ranks": profile["Ranks"], "DeadlineSeconds": 3300, "DeadlineMarginSeconds": 120,
            "MinimumMemAvailableBytes": profile["MinimumMemAvailableBytes"],
            "Binary": binaries["old"][0], "BinarySHA256": OLD_SHA, "MPIExec": recorded["MPIExec"],
            "NewBinary": binaries["new"][0], "NewBinarySHA256": new_sha,
            "PinnedSHA256": pinned, "Stages": stages, "StageNames": [s["Name"] for s in stages],
            "ReducerBlockSize": 48, "ReducerBlockSizeRule": build_plan.REDUCER_BLOCK_SIZE_RULE,
            "MeshSHA256": recorded["MeshSHA256"], "MeshRemote": recorded["MeshRemote"], "MeshLocal": recorded["MeshLocal"]}
    (out / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    (out / "job.pbs").write_text(build_plan.render_job_script(profile=profile, remote_root=ROOT, remote_case_root=CASE,
                                                             runner=f"{RUN}/run_stages.py", job_name="coupon-perf62-4-two-edge",
                                                             walltime_seconds=3600, job_directory=JOB))
    print(json.dumps([(s["Name"], s["BinarySHA256"][:8], s["Environment"].get("PALACE_RESPONSE_BLOCK_SIZE")) for s in stages]))


if __name__ == "__main__":
    main()
