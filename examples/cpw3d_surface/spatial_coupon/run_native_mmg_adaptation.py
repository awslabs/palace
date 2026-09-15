#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Invoke the reviewed native MMG adapter with the metric recipe's bound hmax.

The adapter's dynamically loaded MMG3D library is resolved from its link table and
rpath entries, so the library actually executed is bound by path and SHA-256.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys


MMG_LIBRARY_NAME = re.compile(r"libmmg3d(?:[.-][0-9]+)*\.(?:dylib|so)(?:\.[0-9]+)*")
# Dynamic-loader search overrides change which library the adapter loads; they are
# removed from the adapter environment so the link/rpath resolution is exact.
LOADER_SEARCH_OVERRIDES = ("DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH",
                           "DYLD_INSERT_LIBRARIES", "DYLD_FRAMEWORK_PATH",
                           "LD_LIBRARY_PATH", "LD_PRELOAD", "LD_AUDIT")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tool_output(command):
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as error:
        raise ValueError(f"Link inspection tool is unavailable: {command[0]}") from error
    if result.returncode != 0:
        raise ValueError(f"Adapter is not an inspectable native executable: {result.stderr.strip()}")
    return result.stdout


def linked_libraries(executable):
    """Return (dependency references, rpath entries) from the executable's link table."""
    executable = str(Path(executable).resolve())
    if sys.platform == "darwin":
        references = []
        for line in _tool_output(["otool", "-L", executable]).splitlines()[1:]:
            line = line.strip()
            if line:
                references.append(line.split(" (compatibility", 1)[0].strip())
        rpaths, pending = [], False
        for line in _tool_output(["otool", "-l", executable]).splitlines():
            fields = line.split()
            if fields[:2] == ["cmd", "LC_RPATH"]:
                pending = True
            elif pending and fields[:1] == ["path"]:
                rpaths.append(line.split("path", 1)[1].rsplit(" (offset", 1)[0].strip())
                pending = False
        return references, rpaths
    references, rpaths = [], []
    for line in _tool_output(["readelf", "-d", executable]).splitlines():
        match = re.search(r"\((NEEDED|RUNPATH|RPATH)\)\s+.*\[(.*)\]\s*$", line)
        if match is None:
            continue
        if match.group(1) == "NEEDED":
            references.append(match.group(2))
        else:
            rpaths.extend(entry for entry in match.group(2).split(":") if entry)
    return references, rpaths


def _search_directories(rpaths, executable_directory):
    directories = []
    for entry in rpaths:
        for token in ("@loader_path", "@executable_path", "$ORIGIN", "${ORIGIN}"):
            entry = entry.replace(token, str(executable_directory))
        directories.append(Path(entry))
    return directories


def resolve_mmg_library(adapter):
    """Resolve the single MMG3D library the adapter will load, failing closed."""
    adapter = Path(adapter).resolve()
    references, rpaths = linked_libraries(adapter)
    candidates = [reference for reference in references
                  if MMG_LIBRARY_NAME.fullmatch(Path(reference).name)]
    if len(candidates) != 1:
        raise ValueError("Adapter must link exactly one MMG3D shared library")
    reference = candidates[0]
    search = _search_directories(rpaths, adapter.parent)
    if reference.startswith("@rpath/"):
        located = [directory / reference[len("@rpath/"):] for directory in search]
    elif reference.startswith(("@loader_path/", "@executable_path/")):
        located = [adapter.parent / reference.split("/", 1)[1]]
    elif Path(reference).is_absolute():
        located = [Path(reference)]
    elif "/" not in reference:
        located = [directory / reference for directory in search]
    else:
        raise ValueError(f"MMG library reference is not resolvable from link/rpath: {reference}")
    existing = [path for path in located if path.is_file()]
    if not existing:
        raise ValueError(f"MMG library reference does not resolve to a file: {reference}")
    resolved = existing[0].resolve()
    return {"LinkReference": reference, "Path": str(resolved), "SHA256": sha256(resolved)}


def effective_far_size(recipe):
    policy = recipe.get("FarFieldBudgetPolicy")
    if not isinstance(policy, dict) or policy.get("Name") != "seed-fraction-far-field-v1":
        raise ValueError("Metric recipe has no reviewed far-field budget policy")
    hmax = policy.get("EffectiveFarSize")
    if (not isinstance(hmax, (int, float)) or isinstance(hmax, bool) or
            not math.isfinite(hmax) or hmax <= 0):
        raise ValueError("Metric recipe has an invalid effective far size")
    if recipe.get("FarSize") != hmax:
        raise ValueError("Recorded far-field policy was not consumed by the metric recipe")
    requested = policy.get("RequestedFarSize")
    pressure = policy.get("Pressure")
    if (not isinstance(requested, (int, float)) or not math.isfinite(requested) or
            requested <= 0 or not isinstance(pressure, (int, float)) or
            not math.isfinite(pressure) or pressure < 1 or
            not math.isclose(hmax, requested * pressure, rel_tol=1e-14, abs_tol=0.0)):
        raise ValueError("Metric recipe far-field policy is internally inconsistent")
    return float(hmax)


def run(adapter, seed, metric, pins, recipe_path, output, receipt, *, mmg_library, hmin,
        hgrad, mode, fixed_triangles=None, hausd=1e-8):
    paths = [Path(value).resolve() for value in
             (adapter, seed, metric, pins, recipe_path, output, receipt, mmg_library)]
    adapter, seed, metric, pins, recipe_path, output, receipt, mmg_library = paths
    if any(not path.is_file() for path in (adapter, seed, metric, pins, recipe_path,
                                           mmg_library)):
        raise ValueError("Every native-adaptation input, adapter, and MMG library must exist")
    if output.exists() or receipt.exists():
        raise ValueError("Native-adaptation outputs must be fresh")
    if (not all(math.isfinite(value) for value in (hmin, hgrad, hausd)) or
            hmin <= 0 or (hgrad != -1 and hgrad <= 1) or hausd <= 0):
        raise ValueError("Invalid native-adaptation controls")
    removed = sorted(name for name in LOADER_SEARCH_OVERRIDES if name in os.environ)
    environment = {name: value for name, value in os.environ.items()
                   if name not in LOADER_SEARCH_OVERRIDES}
    library = resolve_mmg_library(adapter)
    if library["Path"] != str(mmg_library):
        raise ValueError("Adapter link/rpath resolves a different MMG library than --mmg-library")
    adapter_sha256 = sha256(adapter)
    recipe = json.loads(recipe_path.read_text())
    hmax = effective_far_size(recipe)
    if hmax < hmin:
        raise ValueError("Policy hmax is smaller than hmin")
    command = [str(adapter), str(seed), str(metric), str(pins), str(output),
               format(hmin, ".17g"), format(hmax, ".17g"), format(hgrad, ".17g"), mode]
    if mode in ("freeze-matching", "freeze-matching-no-move", "freeze-selected"):
        if fixed_triangles is None or not Path(fixed_triangles).is_file():
            raise ValueError("Selected-surface adaptation requires fixed triangles")
        command.extend([str(Path(fixed_triangles).resolve()), format(hausd, ".17g")])
    elif fixed_triangles is not None:
        raise ValueError("Fixed triangles are only valid for a selected-surface mode")
    result = subprocess.run(command, check=False, env=environment)
    if result.returncode != 0 or not output.is_file():
        raise RuntimeError(f"Native MMG adapter failed with status {result.returncode}")
    if sha256(mmg_library) != library["SHA256"] or sha256(adapter) != adapter_sha256:
        raise RuntimeError("Adapter or MMG library changed during native adaptation")
    record = {
        "Version": 2,
        "RecipeSHA256": sha256(recipe_path),
        "MetricSHA256": sha256(metric),
        "AdapterSHA256": adapter_sha256,
        "MMGLibraryPath": library["Path"],
        "MMGLibrarySHA256": library["SHA256"],
        "MMGLibraryLinkReference": library["LinkReference"],
        "LoaderSearchOverridesRemoved": removed,
        "EffectiveFarSize": hmax,
        "HmaxArgument": command[6],
        "Command": command,
        "OutputSHA256": sha256(output),
    }
    receipt.write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=Path)
    parser.add_argument("metric", type=Path)
    parser.add_argument("pins", type=Path)
    parser.add_argument("recipe", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("receipt", type=Path)
    parser.add_argument("--adapter", type=Path, required=True)
    parser.add_argument("--mmg-library", type=Path, required=True,
                        help="MMG3D shared library the adapter resolves through link/rpath")
    parser.add_argument("--hmin", type=float, required=True)
    parser.add_argument("--hgrad", type=float, required=True)
    parser.add_argument("--mode", default="freeze-selected",
                        choices=("adapt", "no-move", "freeze-surface", "freeze-matching",
                                 "freeze-matching-no-move", "freeze-selected"))
    parser.add_argument("--fixed-triangles", type=Path)
    parser.add_argument("--hausd", type=float, default=1e-8)
    args = parser.parse_args()
    try:
        run(args.adapter, args.seed, args.metric, args.pins, args.recipe, args.output,
            args.receipt, mmg_library=args.mmg_library, hmin=args.hmin, hgrad=args.hgrad,
            mode=args.mode, fixed_triangles=args.fixed_triangles, hausd=args.hausd)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
