#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Replay the plans of record of a completed `coupon-library qualify` run with the CURRENT
planning code and compare them byte for byte (decision 457 (1): one-node plans are
bitwise unchanged for every coupon that fits one node; decision 458: the proof is made
under the KEPT cost model the run planned with, digest-checked).

Per coupon of the run's library-qualification.json: the stage estimate is recomputed
from the recorded entity counts, reducer block size and layout under --cost-model (whose
digest must equal the recorded estimate's CostModel.SHA256), the split under the recorded
job policy, and every job's plan.json and job.pbs are rebuilt from the recorded inputs
(mesh, pins, purpose, binary, MPI wrapper) and compared with the recorded files.  The
record lists, per coupon and job, Identical and the first differing keys.

usage: replay_plans.py --qualification library-qualification.json --cost-model PATH [--local-root DIR]
       [--cluster-profile PATH] --out record.json
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
import job_split  # noqa: E402


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def local_case_root(case, qualification_path, local_root):
    """The local directory of a case root recorded on another machine: the recorded
    root's last component under the qualification record's directory (or --local-root)."""
    base = Path(local_root) if local_root else qualification_path.parent
    return base / Path(case["Root"]).name


def local_path(recorded, case, case_root):
    """A recorded absolute path under the case root mapped to the local case root."""
    relative = Path(recorded).relative_to(case["Root"])
    return case_root / relative


def first_differences(a, b, prefix="", limit=12):
    """The first differing JSON paths between two values."""
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                out.append(f"{prefix}/{key}: {'missing in replay' if key not in a else 'missing in record'}")
            else:
                out.extend(first_differences(a[key], b[key], f"{prefix}/{key}", limit))
            if len(out) >= limit:
                return out[:limit]
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{prefix}: length {len(a)} vs {len(b)}")
        for index, (x, y) in enumerate(zip(a, b)):
            out.extend(first_differences(x, y, f"{prefix}[{index}]", limit))
            if len(out) >= limit:
                return out[:limit]
    elif a != b:
        out.append(f"{prefix}: {str(a)[:80]!r} vs {str(b)[:80]!r}")
    return out[:limit]


def replay_case(case, qualification_path, *, model, profile, local_root):
    case_root = local_case_root(case, qualification_path, local_root)
    record = {"Case": case["Case"], "Status": case["Status"], "LocalRoot": str(case_root), "Jobs": {}}
    estimate_path = case_root / "preflight" / "stage-estimate.json"
    if not estimate_path.is_file():
        record["Skipped"] = f"no stage estimate at {estimate_path}"
        return record
    recorded_estimate = json.loads(estimate_path.read_text())
    if recorded_estimate["CostModel"]["SHA256"] != model["SHA256"]:
        record["Skipped"] = (f"the recorded estimate planned with cost model {recorded_estimate['CostModel']['SHA256'][:12]}, "
                             f"not the replay's {model['SHA256'][:12]}")
        return record
    layout = case["Stages"]
    counts = recorded_estimate["Mesh"]["EntityCounts"]
    block_size = recorded_estimate["ReducerBlockSize"]
    stages = [(item["EstimateKey"], item["Order"], item["Sources"]) for item in layout if item["Kind"] == "response"]
    local_edge = next(item for item in layout if item["Kind"] == "local-edge")
    estimate = estimate_stages.estimate(counts, stages, model=model, profile=profile, block_size=block_size,
                                        local_edge=(local_edge["EstimateKey"], local_edge["Order"], local_edge["Sources"]))
    record["Estimate"] = {"DecisionIdentical": estimate["Decision"] == recorded_estimate["Decision"],
                          "FitsOneJob": estimate["FitsOneJob"], "RecordedFitsOneJob": recorded_estimate["FitsOneJob"],
                          "NodesRequired": estimate["Nodes"]["Required"], "MultiNode": estimate["Nodes"]["MultiNode"]}
    policy = dict(case["JobPolicy"])
    indices = case["Sources"]["Indices"]
    split = job_split.plan_split(indices=indices, layout=layout, estimate=estimate, policy=policy, model=model, profile=profile)
    recorded_split = recorded_estimate["Split"]
    record["Split"] = {"N": split["N"], "RecordedN": recorded_split["N"],
                       "BlocksIdentical": split["Blocks"] == recorded_split["Blocks"],
                       "JobsIdentical": split["Jobs"] == recorded_split["Jobs"]}
    factors = [f"{factor:.1f}" for factor in model["PCGFactors"]]
    jobs_by_name = {job["Name"]: job for job in (split["Jobs"] or [])}
    all_identical = record["Estimate"]["DecisionIdentical"] and record["Split"]["BlocksIdentical"]
    for job in case["Jobs"]:
        plan_path = local_path(job["Plan"], case, case_root)
        job_record = {"Plan": str(plan_path)}
        record["Jobs"][job["Name"]] = job_record
        if not plan_path.is_file():
            job_record["Skipped"] = "plan.json not stored locally"
            all_identical = False
            continue
        stored_text = plan_path.read_text()
        stored = json.loads(stored_text)
        remote_case = stored["MeshRemote"].rsplit("/mesh/", 1)[0]
        remote_root = stored["Binary"].rsplit("/", 1)[0]
        remote_run = remote_case.rsplit("/", 1)[0]
        mesh = {"Remote": stored["MeshRemote"], "SHA256": stored["MeshSHA256"], "Local": stored["MeshLocal"]}
        config_digests, trace_pins = {}, {}
        for path, digest in stored["PinnedSHA256"].items():
            if path.startswith(f"{remote_case}/main/"):
                prefix, name = path[len(f"{remote_case}/main/"):].split("/", 1)
                config_digests.setdefault(prefix, {})[name] = digest
            elif path.startswith(f"{remote_case}/inputs/traces/"):
                trace_pins[path] = digest
        common = dict(case_id=case["Case"], remote_case_root=remote_case, mesh=mesh, stage_layout=layout, estimate=estimate,
                      config_digests=config_digests, trace_pins=trace_pins, profile=profile, binary=stored["Binary"],
                      binary_sha256=stored["BinarySHA256"], mpiexec=stored["MPIExec"], purpose=stored["Purpose"], factors=factors,
                      reducer_block_size=block_size)
        if job["Kind"] == "single":
            replayed = build_plan.build_plan(**common)
            job_directory = None
            job_name = f"{profile['JobNamePrefix']}-{case['Case']}"[:64]
        else:
            split_job = jobs_by_name.get(job["Name"])
            if split_job is None:
                job_record["Skipped"] = f"the replayed split has no job {job['Name']}"
                all_identical = False
                continue
            replayed = build_plan.build_job_plan(job_name=job["Name"], split_job=split_job, **common)
            job_directory = f"{remote_case}/main/jobs/{job['Name']}"
            job_name = f"{profile['JobNamePrefix']}-{case['Case']}-{job['Name']}"[:64]
        replayed_text = json.dumps(replayed, indent=2) + "\n"
        job_record["PlanIdentical"] = replayed_text == stored_text
        if not job_record["PlanIdentical"]:
            job_record["PlanDifferences"] = first_differences(replayed, stored)
        script_path = plan_path.parent / "job.pbs"
        if script_path.is_file():
            script = build_plan.render_job_script(profile=profile, remote_root=remote_root, remote_case_root=remote_case,
                                                  runner=f"{remote_run}/run_stages.py", job_name=job_name,
                                                  walltime_seconds=profile["WalltimeSeconds"], instance_type=replayed["Instance"]["Type"],
                                                  job_directory=job_directory, nodes=replayed.get("Nodes", 1))
            job_record["ScriptIdentical"] = script == script_path.read_text()
        else:
            job_record["ScriptIdentical"] = None
        all_identical = all_identical and job_record["PlanIdentical"] and job_record["ScriptIdentical"] is not False
    record["Identical"] = bool(all_identical)
    return record


def replay(qualification_path, *, cost_model_path, profile_path=None, local_root=None):
    qualification_path = Path(qualification_path)
    qualification = json.loads(qualification_path.read_text())
    model = estimate_stages.load_cost_model(cost_model_path)
    profile = json.loads(Path(profile_path or estimate_stages.CLUSTER_PROFILE).read_text())
    cases = [replay_case(case, qualification_path, model=model, profile=profile, local_root=local_root)
             for case in qualification["Cases"] if case.get("Jobs")]
    compared = [case for case in cases if "Skipped" not in case]
    return {"Qualification": str(qualification_path), "QualificationSHA256": sha256(qualification_path),
            "CostModel": {"Path": model["Path"], "SHA256": model["SHA256"]}, "ClusterProfile": str(profile_path or estimate_stages.CLUSTER_PROFILE),
            "Cases": cases, "Compared": len(compared), "Identical": sum(1 for case in compared if case["Identical"]),
            "Plans": sum(len(case["Jobs"]) for case in compared),
            "PlansIdentical": sum(1 for case in compared for job in case["Jobs"].values() if job.get("PlanIdentical")),
            "ScriptsIdentical": sum(1 for case in compared for job in case["Jobs"].values() if job.get("ScriptIdentical")),
            "AllIdentical": bool(compared) and all(case["Identical"] for case in compared)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--qualification", type=Path, required=True, action="append",
                        help="library-qualification.json of a completed run (repeatable)")
    parser.add_argument("--cost-model", type=Path, required=True, help="the cost model the run planned with (digest-checked)")
    parser.add_argument("--cluster-profile", type=Path, default=estimate_stages.CLUSTER_PROFILE)
    parser.add_argument("--local-root", type=Path, help="directory holding the case roots (default: next to each qualification record)")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    records = [replay(path, cost_model_path=args.cost_model, profile_path=args.cluster_profile, local_root=args.local_root)
               for path in args.qualification]
    out = {"Runs": records, "Compared": sum(r["Compared"] for r in records), "Identical": sum(r["Identical"] for r in records),
           "Plans": sum(r["Plans"] for r in records), "PlansIdentical": sum(r["PlansIdentical"] for r in records),
           "ScriptsIdentical": sum(r["ScriptsIdentical"] for r in records),
           "AllIdentical": all(r["AllIdentical"] for r in records) and any(r["Compared"] for r in records)}
    args.out.write_text(json.dumps(out, indent=2) + "\n")
    print(f"replayed {out['Compared']} coupons, {out['Plans']} plans: {out['PlansIdentical']} plans and {out['ScriptsIdentical']} scripts "
          f"byte-identical; coupons identical {out['Identical']} / {out['Compared']}")
    return 0 if out["AllIdentical"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
