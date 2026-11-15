# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""A FAKE paladin CLI for test_paladin_backend.py (Y-SPEC s2 TESTS): it answers
`palace submit-detailed`, `palace describe-run` and `palace cancel` with envelopes of the
shapes the installed paladin prints (--json), keeps its runs in the JSON file named by
FAKE_PALADIN_STATE and appends every argv to FAKE_PALADIN_LOG (one JSON line per call).

The state's "Plan" maps a run label to its scripted end: "completed" (the default),
"failed", "hang" (stays running until cancelled) or "refuse" (the submit fails).  A run is
'queued' at the submit, 'running' at its first describe-run (started = now) and ends at the
next one; its run_dir (under the state's "Root") holds a log.out and a postpro/ table.
"""
import datetime
import json
import os
import sys


def envelope(command, data, ok=True, error=None):
    return {"ok": ok, "command": command, "schema_version": "1", "data": data,
            "error": error, "warnings": []}


def option(argv, name):
    return argv[argv.index(name) + 1] if name in argv else None


def main(argv):
    state_path = os.environ["FAKE_PALADIN_STATE"]
    with open(os.environ["FAKE_PALADIN_LOG"], "a") as f:
        f.write(json.dumps(argv) + "\n")
    with open(state_path) as f:
        state = json.load(f)
    runs = state.setdefault("Runs", {})
    command = " ".join(argv[:2])
    now = datetime.datetime.now(datetime.timezone.utc)
    if command == "palace submit-detailed":
        label = option(argv, "--label")
        if state.get("Plan", {}).get(label) == "refuse":
            print(json.dumps(envelope(command, None, ok=False, error="refused by the fake")))
            return 2
        run_id = state.get("NextID", 1000)
        state["NextID"] = run_id + 1
        run_dir = os.path.join(state["Root"], f"{label}_{run_id}")
        os.makedirs(os.path.join(run_dir, "postpro"))
        with open(argv[2]) as f:
            config = json.load(f)
        with open(os.path.join(run_dir, "log.out"), "w") as f:
            f.write("H1 (p = 2): 1,234, ND (p = 2): 5,678, RT (p = 2): 9,012\n"
                    "PCG solver converged in 7 iterations\n")
        with open(os.path.join(run_dir, "postpro", "domain-response-matrix.csv"), "w") as f:
            f.write("i,j,value\n1,1,2.0\n")
        with open(os.path.join(run_dir, "postpro", "run-config.json"), "w") as f:
            json.dump(config, f)
        runs[str(run_id)] = {"label": label, "status": "queued", "started": None,
                             "completed": None, "nodes": int(option(argv, "--nodes")),
                             "run_dir": run_dir,
                             "end": state.get("Plan", {}).get(label, "completed")}
        data = {"jobs": [{"run_id": run_id, "job_id": str(run_id + 5000), "label": label,
                          "status": "queued", "nodes": runs[str(run_id)]["nodes"]}]}
    elif command == "palace describe-run":
        run = runs[argv[2]]
        if run["status"] == "queued":
            run["status"], run["started"] = "running", now.isoformat()
        elif run["status"] == "running" and run["end"] != "hang":
            run["status"] = run["end"]
            started = datetime.datetime.fromisoformat(run["started"])
            run["completed"] = (started + datetime.timedelta(seconds=600)).isoformat()
        data = {"run_id": int(argv[2]), "label": run["label"], "status": run["status"],
                "started": run["started"], "completed": run["completed"],
                "block": {"backend_data": {"run_dir": run["run_dir"]}}}
    elif command == "palace cancel":
        run = runs[argv[2]]
        run["status"], run["completed"] = "canceled", now.isoformat()
        data = {"results": [{"run_id": int(argv[2]), "status": "canceled", "terminal": True,
                             "ok": True, "error": None}]}
    else:
        print(json.dumps(envelope(command, None, ok=False, error="unknown to the fake")))
        return 2
    with open(state_path, "w") as f:
        json.dump(state, f)
    print(json.dumps(envelope(command, data)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
