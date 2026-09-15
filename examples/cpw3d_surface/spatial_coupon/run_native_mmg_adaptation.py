#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Invoke the reviewed native MMG adapter with the metric recipe's bound hmax."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def run(adapter, seed, metric, pins, recipe_path, output, receipt, *, hmin, hgrad,
        mode, fixed_triangles=None, hausd=1e-8):
    paths = [Path(value).resolve() for value in
             (adapter, seed, metric, pins, recipe_path, output, receipt)]
    adapter, seed, metric, pins, recipe_path, output, receipt = paths
    if any(not path.is_file() for path in (adapter, seed, metric, pins, recipe_path)):
        raise ValueError("Every native-adaptation input and adapter must exist")
    if output.exists() or receipt.exists():
        raise ValueError("Native-adaptation outputs must be fresh")
    if (not all(math.isfinite(value) for value in (hmin, hgrad, hausd)) or
            hmin <= 0 or (hgrad != -1 and hgrad <= 1) or hausd <= 0):
        raise ValueError("Invalid native-adaptation controls")
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
    result = subprocess.run(command, check=False)
    if result.returncode != 0 or not output.is_file():
        raise RuntimeError(f"Native MMG adapter failed with status {result.returncode}")
    record = {
        "Version": 1,
        "RecipeSHA256": sha256(recipe_path),
        "MetricSHA256": sha256(metric),
        "AdapterSHA256": sha256(adapter),
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
            args.receipt, hmin=args.hmin, hgrad=args.hgrad, mode=args.mode,
            fixed_triangles=args.fixed_triangles, hausd=args.hausd)
    except (OSError, ValueError, RuntimeError, json.JSONDecodeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
