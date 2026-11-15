# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""A FAKE paladin CLI for test_paladin_backend.py (Y-SPEC s2 TESTS): it answers `palace
submit-detailed`, `palace describe-run`, `palace cancel`, `palace my-jobs` and `--version` with
envelopes of the shapes paladin 0.2.16 prints (--json; read from real rows), keeps its runs in
the JSON file named by FAKE_PALADIN_STATE and appends every argv to FAKE_PALADIN_LOG.

The state's "Plan" maps a run label to its scripted end:
  completed (default) | failed (Palace failed: no time report, runner exit 1) |
  failed-upload (Palace completed, the upload harness failed: PREREG-ESSURR-NOTE-1's case) |
  oom (out_of_memory) | hang (running until cancelled) | refuse (the submit is refused, no row) |
  transport (every submit: an envelope ok false + a submission_failed row without a job id) |
  transport-once (the first submit only) | short-archive (an archive writer leaves one file out) |
  drift (paladin ran a configuration with one more changed leaf).
The state's "Lookups" is a list of answers for the next `palace my-jobs` calls, one consumed per
call (decision 1820 (1); the review r1799's R1 / R2 / R3): empty (rc 0, no rows) | omitted (rc 0,
the first row dropped, the warning rows_omitted: the installed client's form) | nokey (rc 0, data
without 'jobs') | fail (rc 1, ok false).  An empty list or none = the true listing.
A run is 'queued' at the submit, 'running' at its first describe-run (started = now) and ends
at the next one (completed = started + 600 s), when its run_dir is written: log.out, console.log,
palace_config.json (paladin's rewrite of Mesh / Output / Device / SaveAdaptMesh) and postpro/
(palace.json, palace_config_resolved.json; the two response matrices of a reducer or a
response-matrix solve); an archive writer writes source-NNNNNN-rank-NNNNNN-V.bin files.
"""
import copy
import datetime
import json
import os
import sys

LOG = """Running with {ranks} MPI processes
H1 (p = 2): 1,234, ND (p = 2): 5,678, RT (p = 2): 9,012
PCG solver converged in 7 iterations
{tail}"""
TIME_REPORT = """Elapsed Time Report (s)   Min.   Max.   Avg.
-------------------------------------------------
Total                     10.0   10.0   10.0
"""


def envelope(command, data, ok=True, error=None, warnings=()):
    return {"ok": ok, "command": command, "schema_version": "1", "data": data, "error": error,
            "warnings": list(warnings)}


def option(argv, name):
    return argv[argv.index(name) + 1] if name in argv else None


def metadata_of(argv):
    pairs = [argv[i + 1] for i, part in enumerate(argv) if part == "-m"]
    return dict(pair.split("=", 1) for pair in pairs)


def write_outputs(run, end):
    run_dir = run["run_dir"]
    config = run["config"]
    ranks = run["ranks"]
    os.makedirs(os.path.join(run_dir, "postpro"), exist_ok=True)
    finished = end in ("completed", "failed-upload", "short-archive", "drift")
    with open(os.path.join(run_dir, "log.out"), "w") as f:
        f.write(LOG.format(ranks=ranks, tail=TIME_REPORT if finished else "MPI_ABORT\n"))
    with open(os.path.join(run_dir, "console.log"), "w") as f:
        f.write("Palace SIF:     dev/fake/neoverse_v1.sif\nPaladin addon:  paladin:0.2.16\n")
        if end == "failed-upload":
            f.write("ERROR qdash_das.palace.runner.process_harness: harness upload final tick failed for run\n"
                    "Traceback (most recent call last):\n"
                    "pyarrow.lib.ArrowTypeError: (\"Expected bytes, got a 'float' object\")\n")
        f.write("Runner finished (exit 0)\n" if finished else "Runner finished (exit 1)\n")
    ran = copy.deepcopy(config)
    ran.setdefault("Model", {})["Mesh"] = "mesh.msh"
    ran["Model"].setdefault("Refinement", {})["SaveAdaptMesh"] = True
    ran.setdefault("Problem", {})["Output"] = "postpro"
    ran.setdefault("Solver", {})["Device"] = "CPU"
    if end == "drift":
        ran["Solver"]["Order"] = 99
    with open(os.path.join(run_dir, "palace_config.json"), "w") as f:
        json.dump(ran, f)
    if not finished:
        return
    electrostatic = config.get("Solver", {}).get("Electrostatic", {})
    archive = electrostatic.get("ResponseArchive", {})
    with open(os.path.join(run_dir, "postpro", "palace.json"), "w") as f:
        json.dump({"Fake": True}, f)
    with open(os.path.join(run_dir, "postpro", "palace_config_resolved.json"), "w") as f:
        json.dump(ran, f)
    if archive.get("Mode") == "Reduce" or (not archive.get("Mode") and electrostatic.get("ResponseMatrix")):
        for name in ("domain-response-matrix.csv", "surface-response-matrix.csv"):
            with open(os.path.join(run_dir, "postpro", name), "w") as f:
                f.write(f"i,j,value\n1,1,{len(name)}.0\n")
    if archive.get("Mode") == "Archive":
        os.makedirs(archive["Directory"], exist_ok=True)
        sources = [entry["Index"] for entry in config.get("Boundaries", {}).get("PrescribedPotential", [])]
        for number, source in enumerate(sources):
            for rank in range(ranks):
                if end == "short-archive" and number == 0 and rank == 0:
                    continue
                with open(os.path.join(archive["Directory"], f"source-{source:06d}-rank-{rank:06d}-V.bin"), "wb") as f:
                    f.write(bytes([1]) * (64 + rank))


def main(argv):
    state_path = os.environ["FAKE_PALADIN_STATE"]
    with open(os.environ["FAKE_PALADIN_LOG"], "a") as f:
        f.write(json.dumps(argv) + "\n")
    if argv == ["--version"]:
        print("0.2.16")
        return 0
    with open(state_path) as f:
        state = json.load(f)
    runs = state.setdefault("Runs", {})
    command = " ".join(argv[:2])
    now = datetime.datetime.now(datetime.timezone.utc)
    if command == "palace submit-detailed":
        label = option(argv, "--label")
        end = state.get("Plan", {}).get(label, "completed")
        if end == "refuse":
            print(json.dumps(envelope(command, None, ok=False, error="refused by the fake")))
            return 2
        run_id = state.get("NextID", 1000)
        state["NextID"] = run_id + 1
        transport = end == "transport" or (end == "transport-once" and not state.get("Transported", {}).get(label))
        with open(argv[2]) as f:
            config = json.load(f)
        runs[str(run_id)] = {"label": label, "status": "submission_failed" if transport else "queued", "started": None,
                             "completed": None, "nodes": int(option(argv, "--nodes")),
                             "ranks": int(option(argv, "--ranks-per-node")), "instance_type": option(argv, "--instance-type"),
                             "job_id": None if transport else str(run_id + 5000), "metadata": metadata_of(argv),
                             "run_dir": os.path.join(state["Root"], f"{label}_{run_id}"), "config": config,
                             "end": "completed" if end == "transport-once" else end}
        if transport:
            state.setdefault("Transported", {})[label] = True
            with open(state_path, "w") as f:
                json.dump(state, f)
            print(json.dumps(envelope(command, None, ok=False, error="partial_submit: 1 failed transport")))
            return 1
        os.makedirs(runs[str(run_id)]["run_dir"])
        data = {"jobs": [{"run_id": run_id, "job_id": str(run_id + 5000), "label": label, "status": "queued",
                          "nodes": runs[str(run_id)]["nodes"]}]}
    elif command == "palace describe-run":
        run = runs[argv[2]]
        if run["status"] == "queued":
            run["status"], run["started"] = "running", now.isoformat()
        elif run["status"] == "running" and run["end"] != "hang":
            run["status"] = {"failed": "failed", "failed-upload": "failed", "oom": "out_of_memory"}.get(run["end"], "completed")
            started = datetime.datetime.fromisoformat(run["started"])
            run["completed"] = (started + datetime.timedelta(seconds=600)).isoformat()
            write_outputs(run, run["end"])
        data = {"run_id": int(argv[2]), "label": run["label"], "status": run["status"], "started": run["started"],
                "completed": run["completed"],
                "block": {"backend_data": {"run_dir": run["run_dir"], "ranks_per_node_resolved": run["ranks"]},
                          "user_metadata": run["metadata"]},
                "output_artifacts": ([{"artifact_id": "A1", "label": "log", "type": "FILE_UPLOAD"},
                                      {"artifact_id": "A2", "label": "palace", "type": "DICT"}]
                                     if run["status"] == "completed" else [])}
    elif command == "palace cancel":
        run = runs[argv[2]]
        run["status"], run["completed"] = "canceled", now.isoformat()
        os.makedirs(run["run_dir"], exist_ok=True)
        data = {"results": [{"run_id": int(argv[2]), "status": "canceled", "terminal": True, "ok": True, "error": None}]}
    elif command == "palace my-jobs":
        wanted = metadata_of(argv)
        data = {"jobs": [{"run_id": int(run_id), "job_id": run["job_id"], "label": run["label"], "status": run["status"],
                          "nodes": run["nodes"], "instance_type": run["instance_type"], "metadata": run["metadata"],
                          "started": run["started"], "completed": run["completed"]}
                         for run_id, run in sorted(runs.items(), reverse=True)
                         if all(run["metadata"].get(key) == value for key, value in wanted.items())]}
        answer = state.get("Lookups", []).pop(0) if state.get("Lookups") else None
        warnings = []
        if answer == "empty":
            data = {"jobs": []}
        elif answer == "omitted":
            data = {"jobs": data["jobs"][1:]}
            warnings = [{"code": "rows_omitted", "message": "Skipping PalaceRun in listing (the fake)"}]
        elif answer == "nokey":
            data = {"count": 0}
        elif answer == "fail":
            with open(state_path, "w") as f:
                json.dump(state, f)
            print(json.dumps(envelope(command, None, ok=False, error="the fake's failed listing")))
            return 1
        with open(state_path, "w") as f:
            json.dump(state, f)
        print(json.dumps(envelope(command, data, warnings=warnings)))
        return 0
    else:
        print(json.dumps(envelope(command, None, ok=False, error="unknown to the fake")))
        return 2
    with open(state_path, "w") as f:
        json.dump(state, f)
    print(json.dumps(envelope(command, data)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
