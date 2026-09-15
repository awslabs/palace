#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run a meshing experiment with a process-group timeout and conservative RSS cap."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time


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
    # The launcher can exit before a descendant that ignores SIGTERM. Always
    # terminate the remaining process group, not only a still-running launcher.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


REQUIRED_TOOLCHAIN_ROLES = ("runtime", "mesher", "adaptor", "mmg")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--memory-gib", type=float, default=8)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--artifact", type=Path, action="append", default=[])
    parser.add_argument("--require-complete-toolchain", action="store_true")
    parser.add_argument("--tool", action="append", default=[], metavar="ROLE=PATH")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or not all(math.isfinite(x) and x > 0
                              for x in (args.seconds, args.memory_gib)):
        parser.error("finite positive time/memory limits and a command are required")
    tool_paths = {}
    for value in args.tool:
        if "=" not in value:
            parser.error("--tool must be ROLE=PATH")
        role, value = value.split("=", 1)
        path = Path(value).resolve()
        if role in tool_paths or role not in REQUIRED_TOOLCHAIN_ROLES or not path.is_file():
            parser.error("tool roles must be unique required roles naming existing files")
        tool_paths[role] = path
    if args.require_complete_toolchain:
        if set(tool_paths) != set(REQUIRED_TOOLCHAIN_ROLES):
            parser.error("runtime, mesher, adaptor, and mmg tool digests are required")
        command_files = {Path(value).resolve() for value in command[:2]
                         if Path(value).is_file()}
        if not {tool_paths["runtime"], tool_paths["mesher"]} <= command_files:
            parser.error("the command must directly invoke the frozen runtime and mesher")
    elif tool_paths:
        parser.error("--tool requires --require-complete-toolchain")
    if args.log.exists():
        parser.error("refusing to overwrite an experiment log")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", JULIA_NUM_THREADS="1")
    start = time.monotonic()
    peak = 0
    reason = None
    with args.log.open("w") as log:
        process = subprocess.Popen(
            command, env=env, stdout=log, stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            while process.poll() is None:
                peak = max(peak, tree_rss(process.pid))
                if peak > args.memory_gib * 2**30:
                    reason = "memory_limit"
                    break
                if time.monotonic() - start >= args.seconds:
                    reason = "timeout"
                    break
                time.sleep(0.5)
        finally:
            if process.poll() is None:
                stop(process)
        code = process.wait()
    artifacts = {}
    if code == 0 and reason is None:
        for artifact in args.artifact:
            if not artifact.is_file():
                parser.error(f"declared artifact was not produced: {artifact}")
            artifacts[str(artifact.resolve())] = sha256(artifact)
    report = {
        "Version": 2, "Command": command,
        "Environment": {name: env[name] for name in
                        ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "JULIA_NUM_THREADS")},
        "Producer": {"Name": Path(__file__).name,
                     "SHA256": sha256(Path(__file__).resolve())},
        "Seconds": time.monotonic() - start,
        "PeakProcessTreeRSSBytes": peak, "ReturnCode": code, "StopReason": reason,
        "Limits": {"Seconds": args.seconds, "MemoryGiB": args.memory_gib},
        "Toolchain": {role: {"Path": str(path), "SHA256": sha256(path)}
                      for role, path in tool_paths.items()},
        "Artifacts": artifacts,
    }
    args.log.with_suffix(args.log.suffix + ".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(124 if reason == "timeout" else 137 if reason else code)


if __name__ == "__main__":
    main()
