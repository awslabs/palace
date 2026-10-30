#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The Julia executable a `julia` command resolves to.

The `julia` on PATH of a juliaup install is the launcher (`julialauncher`), a
dispatcher that reads the channel configuration and `exec`s the installed version's
binary; juliaup replaces the launcher file on its own schedule (self-update), while
the binary it dispatches to changes only when another Julia version is installed
and selected.  The machine-bound Julia runtime identity of the stage manifests
(`StageToolSHA256` runtime roles; supervisor decision 1014, OI-JULIAUP-SELF-UPDATE-PIN)
is therefore the digest of the resolved binary, never of the launcher: every Julia
stage runs the resolved binary as argv[0] and records it as its runtime tool, and
`refreeze_manifest_tools.py --julia-runtime` freezes the same binary.  A launcher
given to the resolver is run once to report its own executable and version.
"""
from pathlib import Path
import subprocess

# Printed by the Julia process itself: its executable (`Base.julia_cmd()[1]`) and version.
JULIA_RUNTIME_QUERY = 'print(joinpath(Sys.BINDIR, Base.julia_exename()), "\\n", VERSION)'


def resolve_julia_runtime(julia):
    """`(executable, version)` of the Julia binary that running `julia` executes:
    the resolved path of `joinpath(Sys.BINDIR, Base.julia_exename())` and the
    `VERSION` string, both reported by `julia` itself.  `julia` may be a launcher,
    a symlink or the binary (which resolves to itself)."""
    julia = Path(julia)
    if not julia.is_file():
        raise ValueError(f"Julia runtime executable does not exist: {julia}")
    try:
        result = subprocess.run([str(julia), "--startup-file=no", "-e", JULIA_RUNTIME_QUERY],
                                capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValueError(f"Julia runtime {julia} did not report its executable: {error}") from error
    lines = result.stdout.splitlines()
    if result.returncode != 0 or len(lines) != 2 or not lines[0] or not lines[1]:
        raise ValueError(f"Julia runtime {julia} did not report its executable and version "
                         f"(exit {result.returncode}): {result.stderr.strip()}")
    executable = Path(lines[0]).resolve()
    if not executable.is_file():
        raise ValueError(f"Julia runtime {julia} reports a missing executable: {executable}")
    return executable, lines[1]
