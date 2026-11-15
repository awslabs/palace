# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Form Y-b of `coupon-library qualify` (lane chainport; USER decision 1674 (2); Y-SPEC s2),
FIRST PIECE: a run_stages plan as paladin runs, ONE RUN PER STAGE.  Paladin runs one
configuration file per run and passes no environment variable and no command-line flag to
Palace, so each stage's PALACE_RESPONSE_* variables become the stage configuration's
Solver.Electrostatic.ResponseArchive keys (the Palace change Y-a, d8a603ae0a).  This piece
PLANS ONLY: --dry-run writes the run configurations and prints the submit commands; it
contacts nothing.  Submit / wait / cancel / status / fetch / resume are the next pieces.

  python3 paladin_backend.py --dry-run PLAN --selector dev-simlapointe-sp-pin-<sha10>
      --lane LANE --out DIR [--nodes N] [--map REMOTE_PREFIX=LOCAL_PREFIX]
"""
import argparse
import copy
import hashlib
import json
import os
import sys

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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("plan")
    parser.add_argument("--dry-run", action="store_true", required=True)
    parser.add_argument("--selector", required=True)
    parser.add_argument("--lane", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--map", action="append", default=[])
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
    for run in runs:
        print(json.dumps({"Stage": run["Stage"], "ResponseArchive": run["ResponseArchive"],
                          "Requires": run["Requires"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
