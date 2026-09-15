#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run one bounded meshing stage with immutable input/output/tool bindings."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time

from mesh_stage_contract import (STAGE_INPUTS, STAGE_OUTPUTS, STAGE_TOOLS, binding,
                                 sha256)


def tree_rss(root):
    result = subprocess.run(
        ["ps", "-axo", "pid=,ppid=,rss="], capture_output=True, text=True,
        check=True, timeout=3,
    )
    rows = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
    children = {root}
    while True:
        added = {pid for pid, parent, _ in rows if parent in children} - children
        if not added:
            break
        children.update(added)
    return sum(rss * 1024 for pid, _, rss in rows if pid in children)


def stop(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        process.wait(timeout=3)
    except subprocess.TimeoutExpired:
        pass
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def _named_paths(values, option):
    result = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{option} must be NAME=PATH in staged mode")
        name, value = value.split("=", 1)
        if not name or name in result:
            raise ValueError(f"{option} names must be nonempty and unique")
        result[name] = Path(value).resolve()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--memory-gib", type=float, default=8)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--stage", choices=tuple(STAGE_TOOLS))
    parser.add_argument("--input", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--artifact", action="append", default=[], metavar="NAME=PATH")
    parser.add_argument("--tool", action="append", default=[], metavar="ROLE=PATH")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or not all(math.isfinite(x) and x > 0
                              for x in (args.seconds, args.memory_gib)):
        parser.error("finite positive time/memory limits and a command are required")
    try:
        inputs = _named_paths(args.input, "--input")
        tools = _named_paths(args.tool, "--tool")
        artifacts = (_named_paths(args.artifact, "--artifact") if args.stage else
                     {str(Path(value).resolve()): Path(value).resolve()
                      for value in args.artifact})
    except ValueError as error:
        parser.error(str(error))
    if args.stage is None:
        if inputs or tools or any("=" in value for value in args.artifact):
            parser.error("named inputs, artifacts, and tools require --stage")
    else:
        if set(inputs) != STAGE_INPUTS[args.stage] or set(artifacts) != STAGE_OUTPUTS[args.stage]:
            parser.error("stage input/output names do not match the frozen stage contract")
        if set(tools) != STAGE_TOOLS[args.stage] or any(not path.is_file()
                                                         for path in tools.values()):
            parser.error("stage tool roles must name every required existing tool")
        command_paths = {Path(value).resolve() for value in command
                         if isinstance(value, str) and Path(value).is_file()}
        if not set(tools.values()) <= command_paths:
            parser.error("every frozen stage tool must appear in the executed command")
        if any(not path.is_file() for path in inputs.values()):
            parser.error("every declared stage input must exist before execution")
    if args.log.exists():
        parser.error("refusing to overwrite an experiment log")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", JULIA_NUM_THREADS="1")
    input_bindings = {name: binding(path) for name, path in inputs.items()}
    start = time.monotonic()
    peak = 0
    reason = None
    with args.log.open("w") as log:
        process = subprocess.Popen(command, env=env, stdout=log, stderr=subprocess.STDOUT,
                                   start_new_session=True)
        try:
            while process.poll() is None:
                peak = max(peak, tree_rss(process.pid))
                if peak > args.memory_gib * 2**30:
                    reason = "memory_limit"
                    break
                if time.monotonic() - start >= args.seconds:
                    reason = "timeout"
                    break
                time.sleep(0.05)
        finally:
            if process.poll() is None:
                stop(process)
        code = process.wait()
    artifact_bindings = {}
    if code == 0 and reason is None:
        for name, artifact in artifacts.items():
            if not artifact.is_file():
                parser.error(f"declared artifact was not produced: {artifact}")
            artifact_bindings[name] = binding(artifact)
    report = {
        "Version": 3 if args.stage else 2,
        "Command": command,
        "Environment": {name: env[name] for name in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "JULIA_NUM_THREADS")},
        "Producer": {"Name": Path(__file__).name,
                     "SHA256": sha256(Path(__file__).resolve())},
        "Seconds": time.monotonic() - start,
        "PeakProcessTreeRSSBytes": peak,
        "ReturnCode": code,
        "StopReason": reason,
        "Limits": {"Seconds": args.seconds, "MemoryGiB": args.memory_gib},
        "Artifacts": artifact_bindings,
    }
    if args.stage:
        report.update({"Stage": args.stage, "Inputs": input_bindings,
                       "Tools": {role: binding(path) for role, path in tools.items()}})
    args.log.with_suffix(args.log.suffix + ".json").write_text(
        json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(124 if reason == "timeout" else 137 if reason else code)


if __name__ == "__main__":
    main()
