#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run a meshing experiment with a process-group timeout and conservative RSS cap."""
import argparse
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180)
    parser.add_argument("--memory-gib", type=float, default=8)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or not all(math.isfinite(x) and x>0 for x in (args.seconds,args.memory_gib)):
        parser.error("finite positive time/memory limits and a command are required")
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
    report = {
        "Command": command, "Seconds": time.monotonic() - start,
        "PeakProcessTreeRSSBytes": peak, "ReturnCode": code, "StopReason": reason,
        "Limits": {"Seconds": args.seconds, "MemoryGiB": args.memory_gib},
    }
    args.log.with_suffix(args.log.suffix + ".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    raise SystemExit(124 if reason == "timeout" else 137 if reason else code)


if __name__ == "__main__":
    main()
