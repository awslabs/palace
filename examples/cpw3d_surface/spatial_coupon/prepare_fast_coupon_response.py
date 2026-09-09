#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare an exact, memory-bounded, two-command coupon response build."""

import argparse
import json
import pathlib
import shlex
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent
PREPARE = ROOT / "prepare_blockwise_response.py"


def load(path):
    return json.loads(path.read_text())


def discover(root):
    candidates = [
        ("thin", root / "spatial_thin.json"),
        ("fabricated", root / "spatial_fabricated.json"),
        ("thin", root / "thin.json"),
        ("fabricated", root / "fabricated.json"),
    ]
    result = []
    seen = set()
    for name, path in candidates:
        if path.is_file() and name not in seen:
            result.append((name, path.resolve()))
            seen.add(name)
    if {name for name, _ in result} != {"thin", "fabricated"}:
        raise ValueError(
            f"{root} must contain thin/fabricated or spatial_thin/spatial_fabricated configs"
        )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coupon", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--ranks", type=int, default=6)
    parser.add_argument("--reduction-block-size", type=int, default=6)
    parser.add_argument("--local", action="store_true")
    parser.add_argument("--keep-archives", action="store_true")
    args = parser.parse_args()
    if args.ranks <= 0 or args.reduction_block_size <= 0:
        parser.error("ranks and reduction block size must be positive")

    coupon = args.coupon.expanduser().resolve()
    output = (
        args.output.expanduser().resolve()
        if args.output
        else coupon / "fast-response"
    )
    output.mkdir(parents=True, exist_ok=True)
    cases = discover(coupon)
    manifest_cases = []
    worker_commands = []
    reducer_commands = []
    install_commands = []
    for name, config in cases:
        case_output = output / name
        command = [
            "python3",
            str(PREPARE),
            str(config),
            "--output",
            str(case_output),
            "--execution-mode",
            "streaming",
            "--ranks",
            str(args.ranks),
            "--reduction-block-size",
            str(args.reduction_block_size),
        ]
        if args.local:
            command.append("--local")
        subprocess.run(command, check=True)
        source = load(config)
        destination = pathlib.Path(source["Problem"]["Output"])
        if not destination.is_absolute():
            destination = (config.parent / destination).resolve()
        worker_commands.append(shlex.quote(str(case_output / "run-blocks.sh")))
        reducer_commands.append(shlex.quote(str(case_output / "run-reducer.sh")))
        reducer_output = case_output / "postpro" / "reducer"
        install_commands.extend(
            [
                f"mkdir -p {shlex.quote(str(destination))}",
                "cp "
                f"{shlex.quote(str(reducer_output / 'domain-response-matrix.csv'))} "
                f"{shlex.quote(str(destination / 'domain-response-matrix.csv'))}",
                "cp "
                f"{shlex.quote(str(reducer_output / 'surface-response-matrix.csv'))} "
                f"{shlex.quote(str(destination / 'surface-response-matrix.csv'))}",
                "cp "
                f"{shlex.quote(str(reducer_output / 'surface-response-matrix.csv'))} "
                f"{shlex.quote(str(destination / 'surface-response-matrix-aggregate.csv'))}",
            ]
        )
        manifest_cases.append(
            {
                "Name": name,
                "Config": str(config),
                "WorkDirectory": str(case_output),
                "Destination": str(destination),
            }
        )

    script = """#!/bin/bash
set -euo pipefail
: "${PALACE:?Set PALACE to the Palace executable}"
MPIEXEC=${MPIEXEC:-mpirun}
"""
    script += "\n# Stream exact basis fields while reusing one operator/preconditioner.\n"
    script += "\n".join(worker_commands) + "\n"
    script += "\n# Reduce bounded source blocks into the exact dense matrices.\n"
    script += "\n".join(reducer_commands) + "\n"
    script += "\n# Install matrices at the paths referenced by process-library.json.\n"
    script += "\n".join(install_commands) + "\n"
    if not args.keep_archives:
        script += "\n# Remove exact volume-field archives after successful installation.\n"
        script += "\n".join(
            f"rm -rf {shlex.quote(str(output / name / 'archive'))}"
            for name, _ in cases
        ) + "\n"
    run = output / "run-all.sh"
    run.write_text(script)
    run.chmod(0o755)
    manifest = {
        "Version": 1,
        "Coupon": str(coupon),
        "Output": str(output),
        "Ranks": args.ranks,
        "ReductionBlockSize": args.reduction_block_size,
        "Local": args.local,
        "Exact": True,
        "KeepArchives": args.keep_archives,
        "Cases": manifest_cases,
        "RunScript": str(run),
    }
    path = output / "fast-response-manifest.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(path)


if __name__ == "__main__":
    main()
