# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form Y-b of `coupon-library qualify` (lane chainport; USER decision 1674 (2); Y-SPEC s2),
FIRST PIECE: a run_stages plan as paladin runs, ONE RUN PER STAGE.  Paladin runs one
configuration file per run and passes no environment variable and no command-line flag to
Palace, so each stage's PALACE_RESPONSE_* variables become the stage configuration's
Solver.Electrostatic.ResponseArchive keys (the Palace change Y-a, d8a603ae0a).  --dry-run
writes the run configurations and prints the submit commands; it contacts nothing.  --submit
drives the planned runs in order through the paladin CLI: submit-detailed (the run_id is
recorded in status.json before anything else; a stage with a run_id is never submitted
again), the wait (describe-run polled until the run ends; past the stage's CapSeconds from
its start the run is cancelled), and the fetch of <run_dir>/postpro and log.out by rsync.
A failed or cancelled run stops the drive; nothing is resubmitted.

  python3 paladin_backend.py --dry-run PLAN --selector dev-simlapointe-sp-pin-<sha10>
      --lane LANE --out DIR [--nodes N] [--map REMOTE_PREFIX=LOCAL_PREFIX]
  python3 paladin_backend.py --submit PLAN ... [--paladin-cli PATH] [--remote-host HOST]
      [--poll-seconds S]
"""
import argparse
import copy
import datetime
import hashlib
import json
import os
import subprocess
import sys
import time

# The stage variables of build_plan.py and their keys.  Every other variable is refused.
MODE_VARIABLES = {
    "PALACE_RESPONSE_ARCHIVE_ONLY": "Archive",
    "PALACE_RESPONSE_REDUCE_ONLY": "Reduce",
    "PALACE_RESPONSE_ESTIMATE_ONLY": "Estimate",
}
PATH_VARIABLES = {
    "PALACE_RESPONSE_ARCHIVE_DIR": "Directory",
    "PALACE_RESPONSE_ESTIMATE_REQUEST": "EstimateRequest",
}
SWITCH_VARIABLES = {
    "PALACE_RESPONSE_SOURCE_TIMING": "SourceTiming",
    "PALACE_RESPONSE_RECYCLE_INITIAL_GUESS": "RecycleInitialGuess",
}
INTEGER_VARIABLES = {"PALACE_RESPONSE_BLOCK_SIZE": "BlockSize"}


def response_archive_keys(environment):
    """The ResponseArchive block equivalent to a stage's Environment ({} = none)."""
    keys = {}
    for name, value in sorted(environment.items()):
        if name in MODE_VARIABLES:
            if value != "1":
                raise ValueError(f"{name}={value!r}: a plan sets a mode variable to '1' only")
            if "Mode" in keys:
                raise ValueError(f"{name}: more than one response-archive mode")
            keys["Mode"] = MODE_VARIABLES[name]
        elif name in PATH_VARIABLES:
            if not os.path.isabs(value):
                # Worker and reducer are different paladin runs, each in its own run_dir.
                raise ValueError(f"{name}={value!r}: a paladin run needs an absolute path")
            keys[PATH_VARIABLES[name]] = value
        elif name in SWITCH_VARIABLES:
            if value != "1":
                raise ValueError(f"{name}={value!r}: a plan sets a switch to '1' only")
            keys[SWITCH_VARIABLES[name]] = True
        elif name in INTEGER_VARIABLES:
            number = int(value)
            if number < 1:
                raise ValueError(f"{name}={value!r} must be positive")
            keys[INTEGER_VARIABLES[name]] = number
        else:
            raise ValueError(f"{name}: no configuration key carries this variable")
    return keys


def run_config(config, environment):
    """The stage configuration with its Environment as ResponseArchive keys."""
    keys = response_archive_keys(environment)
    out = copy.deepcopy(config)
    electrostatic = out.setdefault("Solver", {}).setdefault("Electrostatic", {})
    if "ResponseArchive" in electrostatic:
        raise ValueError("the stage configuration already has a ResponseArchive block")
    if keys:
        electrostatic["ResponseArchive"] = keys
    return out


def map_path(path, mapping):
    for remote, local in mapping:
        if path == remote or path.startswith(remote.rstrip("/") + "/"):
            return local.rstrip("/") + path[len(remote.rstrip("/")):]
    return path


def plan_runs(plan, *, plan_sha256, selector, lane, nodes, out, mapping=()):
    """The paladin runs of a plan (one per stage), their run configurations written to out."""
    ranks = int(plan["Ranks"])
    if nodes < 1 or ranks % nodes:
        raise ValueError(f"{ranks} ranks do not split over {nodes} nodes")
    instance = plan["Instance"]["Type"] if isinstance(plan["Instance"], dict) else plan["Instance"]
    case = plan["Case"]
    names = [stage["Name"] for stage in plan["Stages"]]
    runs = []
    for stage in plan["Stages"]:
        for required in stage.get("Requires", []):
            if required not in names[: names.index(stage["Name"])]:
                raise ValueError(f"{stage['Name']} requires {required}, not an earlier stage")
        with open(map_path(stage["Config"], mapping)) as f:
            config = run_config(json.load(f), stage.get("Environment", {}))
        path = os.path.join(out, stage["Name"] + ".json")
        with open(path, "w") as f:
            json.dump(config, f, indent=2)
        submit = ["paladin", "palace", "submit-detailed", path,
                  "--instance-type", instance, "--nodes", str(nodes),
                  "--ranks-per-node", str(ranks // nodes), "--palace", selector,
                  "--label", stage["Name"],
                  "-m", "programme=spr", "-m", f"lane={lane}", "-m", f"case={case}",
                  "-m", f"stage={stage['Name']}", "-m", f"plan={plan_sha256[:16]}",
                  "--comment", f"{case} {stage['Name']} (coupon-library qualify, paladin form)",
                  "--json"]
        runs.append({"Stage": stage["Name"], "RunConfig": path,
                     "ResponseArchive": config["Solver"]["Electrostatic"].get(
                         "ResponseArchive", {}),
                     "Mesh": config.get("Model", {}).get("Mesh"),
                     "Requires": list(stage.get("Requires", [])),
                     "CapSeconds": stage.get("CapSeconds"),
                     "Nodes": nodes, "RanksPerNode": ranks // nodes, "Submit": submit})
    return runs


TERMINAL = {"completed", "failed", "canceled", "cancelled"}


class PaladinError(RuntimeError):
    """A paladin call that failed, or a run that did not complete: the drive stops."""


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc)


def paladin_call(cli, arguments):
    """Run `cli arguments` (a --json command); the envelope's data, or PaladinError."""
    completed = subprocess.run([cli] + list(arguments), capture_output=True, text=True)
    try:
        envelope = json.loads(completed.stdout)
    except ValueError:
        raise PaladinError(f"{' '.join(arguments[:2])}: no JSON envelope "
                           f"(rc {completed.returncode})") from None
    if completed.returncode != 0 or not envelope.get("ok"):
        raise PaladinError(f"{' '.join(arguments[:2])}: {envelope.get('error')!r} "
                           f"(rc {completed.returncode})")
    return envelope["data"]


def write_status(path, status):
    temporary = path + ".tmp"
    with open(temporary, "w") as f:
        json.dump(status, f, indent=1)
    os.replace(temporary, path)


def wait_run(cli, run_id, cap_seconds, *, poll_seconds, now=utc_now, sleep=time.sleep):
    """Poll describe-run until the run ends; past cap_seconds from its start, cancel it.
    Returns (the last describe-run data, whether this wait cancelled the run)."""
    cancelled = False
    while True:
        data = paladin_call(cli, ["palace", "describe-run", str(run_id), "--json"])
        if str(data.get("status", "")).lower() in TERMINAL:
            return data, cancelled
        started = data.get("started")
        if cap_seconds and started and not cancelled:
            running = (now() - datetime.datetime.fromisoformat(started)).total_seconds()
            if running > float(cap_seconds):
                paladin_call(cli, ["palace", "cancel", str(run_id), "--json"])
                cancelled = True
        sleep(poll_seconds)


def fetch_run(run_dir, local, *, remote_host, rsync="rsync"):
    """rsync <run_dir>/postpro and <run_dir>/log.out into local.  The archives are not in
    the run_dir (their Directory is an absolute path under the remote case root)."""
    os.makedirs(local, exist_ok=True)
    source = f"{remote_host}:{run_dir}" if remote_host else run_dir
    transport = ["-e", "ssh -o BatchMode=yes -o ConnectTimeout=20"] if remote_host else []
    for part, target in (("postpro/", os.path.join(local, "postpro") + "/"),
                         ("log.out", os.path.join(local, "log.out"))):
        command = [rsync, "-a"] + transport + [source.rstrip("/") + "/" + part, target]
        if subprocess.run(command, capture_output=True).returncode != 0:
            raise PaladinError(f"the fetch of {part} from {run_dir} failed")


def log_summary(path):
    """run_stages.parse_log on a fetched log.out (one parser), without its long reports."""
    import run_stages  # the same qualify/ directory

    with open(path, errors="replace") as f:
        parsed = run_stages.parse_log(f.read())
    for long_entry in ("ElapsedTimeReport", "PeakMemoryReport", "Iterations"):
        parsed.pop(long_entry, None)
    return parsed


def drive(runs, *, status_path, case, cli, out, remote_host, poll_seconds, rsync="rsync",
          now=utc_now, sleep=time.sleep):
    """Submit, wait for and fetch the runs IN ORDER.  status.json records each stage's
    RunID right after its submit; a stage with a RunID is never submitted again."""
    status = {"Backend": "paladin", "Case": case, "Stages": {}}
    if os.path.isfile(status_path):
        with open(status_path) as f:
            status = json.load(f)
    stages = status["Stages"]
    by_name = {run["Stage"]: run for run in runs}
    for run in runs:
        name = run["Stage"]
        record = stages.setdefault(name, {})
        if record.get("State") == "completed" and record.get("Fetched"):
            continue
        for required in run["Requires"]:
            if stages.get(required, {}).get("State") != "completed":
                raise PaladinError(f"{name} requires {required}, which has not completed")
            other = by_name[required]
            directory = run["ResponseArchive"].get("Directory")
            if (directory and directory == other["ResponseArchive"].get("Directory") and
                    (run["Nodes"], run["RanksPerNode"]) !=
                    (other["Nodes"], other["RanksPerNode"])):
                raise PaladinError(f"{name} reads the archive of {required} on another "
                                   "partition (the archive is per rank)")
        if "RunID" not in record:
            data = paladin_call(cli, run["Submit"][1:])
            record["RunID"] = int(data["jobs"][0]["run_id"])
            record["Submitted"] = now().isoformat()
            write_status(status_path, status)
        data, cancelled = wait_run(cli, record["RunID"], run.get("CapSeconds"),
                                   poll_seconds=poll_seconds, now=now, sleep=sleep)
        backend_data = data.get("block", {}).get("backend_data", {})
        record.update({"State": "timed-out" if cancelled else str(data["status"]).lower(),
                       "RunDir": backend_data.get("run_dir"),
                       "Started": data.get("started"), "Completed": data.get("completed")})
        if record["Started"] and record["Completed"]:
            seconds = (datetime.datetime.fromisoformat(record["Completed"]) -
                       datetime.datetime.fromisoformat(record["Started"])).total_seconds()
            record["NodeHours"] = run["Nodes"] * seconds / 3600.0
        write_status(status_path, status)
        if record["State"] != "completed":
            raise PaladinError(f"{name}: run {record['RunID']} ended {record['State']}; "
                               "nothing is resubmitted")
        local = os.path.join(out, name)
        fetch_run(record["RunDir"], local, remote_host=remote_host, rsync=rsync)
        record["Fetched"] = local
        record["Log"] = log_summary(os.path.join(local, "log.out"))
        write_status(status_path, status)
    return status

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("plan")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--dry-run", action="store_true")
    mode.add_argument("--submit", action="store_true")
    parser.add_argument("--selector", required=True)
    parser.add_argument("--lane", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--map", action="append", default=[])
    parser.add_argument("--paladin-cli", default="paladin")
    parser.add_argument("--remote-host", default="soca-green-job")
    parser.add_argument("--poll-seconds", type=float, default=60.0)
    args = parser.parse_args(argv)
    with open(args.plan, "rb") as f:
        data = f.read()
    mapping = [tuple(entry.split("=", 1)) for entry in args.map]
    os.makedirs(args.out, exist_ok=True)
    runs = plan_runs(json.loads(data), plan_sha256=hashlib.sha256(data).hexdigest(),
                     selector=args.selector, lane=args.lane, nodes=args.nodes,
                     out=args.out, mapping=mapping)
    with open(os.path.join(args.out, "paladin-runs.json"), "w") as f:
        json.dump(runs, f, indent=1)
    if args.submit:
        status_path = os.path.join(args.out, "status.json")
        code = 0
        try:
            drive(runs, status_path=status_path, case=json.loads(data)["Case"],
                  cli=args.paladin_cli, out=args.out, remote_host=args.remote_host or None,
                  poll_seconds=args.poll_seconds)
        except PaladinError as error:
            print(f"STOP: {error}", file=sys.stderr)
            code = 3
        if os.path.isfile(status_path):
            with open(status_path) as f:
                stages = json.load(f)["Stages"]
            for name, record in stages.items():
                print(json.dumps({"Stage": name, "RunID": record.get("RunID"),
                                  "State": record.get("State"),
                                  "NodeHours": record.get("NodeHours")}))
        return code
    for run in runs:
        print(json.dumps({"Stage": run["Stage"], "ResponseArchive": run["ResponseArchive"],
                          "Requires": run["Requires"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
