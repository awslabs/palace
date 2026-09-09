#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare exact source-parallel response archives and a blockwise reducer config."""

import argparse
import copy
import json
import pathlib
import shlex


def load(path):
    return json.loads(path.read_text())


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def response_sources(config):
    sources = config.get("Boundaries", {}).get("PrescribedPotential")
    if not isinstance(sources, list) or not sources:
        raise ValueError("config has no Boundaries.PrescribedPotential sources")
    indices = [int(source["Index"]) for source in sources]
    if len(set(indices)) != len(indices) or any(index <= 0 for index in indices):
        raise ValueError("prescribed-potential source indices must be unique and positive")
    return sources


def set_output(config, output):
    config.setdefault("Problem", {})["Output"] = str(output)


def set_response_matrix(config, enabled):
    electrostatic = config.setdefault("Solver", {}).setdefault("Electrostatic", {})
    electrostatic["ResponseMatrix"] = bool(enabled)
    electrostatic["AggregateResponseMatrix"] = bool(enabled)
    electrostatic["Save"] = 0


def shell_command(
    config,
    archive,
    ranks,
    reduce=False,
    block_size=None,
    archive_only=False,
    recycle_initial_guess=False,
    local=False,
):
    environment = [f"PALACE_RESPONSE_ARCHIVE_DIR={shlex.quote(str(archive))}"]
    if reduce:
        environment.append("PALACE_RESPONSE_REDUCE_ONLY=1")
        environment.append(f"PALACE_RESPONSE_BLOCK_SIZE={block_size}")
    if archive_only:
        environment.append("PALACE_RESPONSE_ARCHIVE_ONLY=1")
    if recycle_initial_guess:
        environment.append("PALACE_RESPONSE_RECYCLE_INITIAL_GUESS=1")
    prefix = " ".join(environment)
    hostfile = "" if local else ' --hostfile "$PBS_NODEFILE"'
    return (
        f"{prefix} \"$MPIEXEC\" -n {ranks}{hostfile} "
        f"\"$PALACE\" {shlex.quote(str(config))}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    parser.add_argument("--archive", type=pathlib.Path)
    parser.add_argument("--sources-per-block", type=int, default=3)
    parser.add_argument("--reduction-block-size", type=int, default=3)
    parser.add_argument("--ranks", type=int, default=6)
    parser.add_argument(
        "--execution-mode", choices=("streaming", "blocks"), default="streaming"
    )
    parser.add_argument("--local", action="store_true")
    parser.add_argument("--recycle-initial-guess", action="store_true")
    args = parser.parse_args()
    if args.sources_per_block <= 0 or args.reduction_block_size <= 0 or args.ranks <= 0:
        parser.error("block sizes and ranks must be positive")

    config_path = args.config.expanduser().resolve()
    output = args.output.expanduser().resolve()
    archive = (
        args.archive.expanduser().resolve()
        if args.archive
        else output / "archive"
    )
    output.mkdir(parents=True, exist_ok=True)
    archive.mkdir(parents=True, exist_ok=True)
    config = load(config_path)
    sources = response_sources(config)

    blocks = []
    block_commands = []
    source_blocks = (
        [sources]
        if args.execution_mode == "streaming"
        else [
            sources[begin : begin + args.sources_per_block]
            for begin in range(0, len(sources), args.sources_per_block)
        ]
    )
    for number, selected in enumerate(source_blocks):
        block = copy.deepcopy(config)
        block["Boundaries"]["PrescribedPotential"] = selected
        # Keep response-matrix postprocessing enabled so Palace recovers and archives the
        # electric flux. Streaming archive-only mode skips the in-process matrix; block
        # mode may emit an incidental within-block matrix. The reducer emits the complete
        # cross-block matrix in either case.
        set_response_matrix(block, True)
        block_path = output / "blocks" / f"block-{number:04d}.json"
        set_output(block, output / "postpro" / f"block-{number:04d}")
        write(block_path, block)
        indices = [int(source["Index"]) for source in selected]
        blocks.append({"Block": number, "Config": str(block_path), "Sources": indices})
        block_commands.append(
            shell_command(
                block_path,
                archive,
                args.ranks,
                archive_only=args.execution_mode == "streaming",
                recycle_initial_guess=(
                    args.execution_mode == "streaming"
                    and args.recycle_initial_guess
                ),
                local=args.local,
            )
        )

    reducer = copy.deepcopy(config)
    set_response_matrix(reducer, True)
    reducer_path = output / "reducer.json"
    set_output(reducer, output / "postpro" / "reducer")
    write(reducer_path, reducer)
    reducer_command = shell_command(
        reducer_path,
        archive,
        args.ranks,
        reduce=True,
        block_size=args.reduction_block_size,
        local=args.local,
    )

    preamble = """#!/bin/bash
set -euo pipefail
: "${PALACE:?Set PALACE to the Palace executable}"
MPIEXEC=${MPIEXEC:-mpirun}
"""
    if not args.local:
        preamble += ': "${PBS_NODEFILE:?Run in a PBS allocation or use --local}"\n'
    (output / "block-commands.txt").write_text("\n".join(block_commands) + "\n")
    (output / "run-blocks.sh").write_text(
        preamble + "\n" + "\n".join(block_commands) + "\n"
    )
    (output / "run-reducer.sh").write_text(preamble + "\n" + reducer_command + "\n")
    (output / "run-blocks.sh").chmod(0o755)
    (output / "run-reducer.sh").chmod(0o755)

    manifest = {
        "Version": 1,
        "SourceConfig": str(config_path),
        "Archive": str(archive),
        "SourceCount": len(sources),
        "ExecutionMode": args.execution_mode,
        "Local": args.local,
        "RecycleInitialGuess": (
            args.execution_mode == "streaming"
            and args.recycle_initial_guess
        ),
        "BlockCount": len(blocks),
        "SourcesPerBlock": args.sources_per_block,
        "ReductionBlockSize": args.reduction_block_size,
        "Ranks": args.ranks,
        "Blocks": blocks,
        "ReducerConfig": str(reducer_path),
        "EstimatedArchiveFieldFiles": 2 * len(sources) * args.ranks,
    }
    write(output / "blockwise-manifest.json", manifest)
    print(output / "blockwise-manifest.json")


if __name__ == "__main__":
    main()
