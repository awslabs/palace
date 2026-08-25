#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Run direct held-out solves and qualify every generated spatial coupon."""

import argparse
import concurrent.futures
import json
import shlex
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
QUALIFIER = ROOT / "qualify_surface_response.py"


def run(command, cwd, check=True):
    command = [str(value) for value in command]
    print("+ " + shlex.join(command) + f" (cwd={cwd})", flush=True)
    return subprocess.run(command, cwd=cwd, check=check).returncode


def palace_command(palace, ranks, config):
    return [palace, *( ["--serial"] if ranks == 1 else ["-np", ranks]), config]


def complete(root):
    names = ("domain-E.csv", "surface-Q.csv", "surface-Q-edge.csv")
    return all(
        (root / "postpro" / f"heldout_spatial_{kind}" / name).is_file()
        for kind in ("thin", "fabricated")
        for name in names
    )


def qualify(root, args):
    try:
        if not complete(root) or args.force:
            for kind in ("thin", "fabricated"):
                run(
                    palace_command(
                        args.palace,
                        args.ranks,
                        root / f"heldout_spatial_{kind}.json",
                    ),
                    root,
                )
        report = root / "heldout-qualification.json"
        code = run(
            [
                sys.executable,
                QUALIFIER,
                root,
                "--output",
                report,
                "--max-heldout-error",
                args.max_heldout_error,
            ],
            root,
            check=False,
        )
        data = json.loads(report.read_text())
        return {
            "Coupon": str(root),
            "Report": str(report),
            "Passed": code == 0 and bool(data.get("Passed")),
            "Failure": None,
        }
    except Exception as error:
        return {
            "Coupon": str(root),
            "Report": None,
            "Passed": False,
            "Failure": f"{type(error).__name__}: {error}",
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("generation", type=Path)
    parser.add_argument("--palace", type=Path, required=True)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--max-heldout-error", type=float, default=10.0)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.ranks <= 0 or args.jobs <= 0 or args.max_heldout_error <= 0.0:
        parser.error("ranks, jobs, and maximum held-out error must be positive")

    generation = args.generation.expanduser().resolve()
    coupons = sorted(
        root.parent
        for root in (generation / "cache").glob("spatial-*/process-library.json")
    )
    if not coupons:
        raise FileNotFoundError(f"No generated spatial coupons found under {generation}")
    if args.jobs == 1:
        results = [qualify(coupon, args) for coupon in coupons]
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as executor:
            results = list(executor.map(lambda coupon: qualify(coupon, args), coupons))

    summary = {
        "Version": 1,
        "Generation": str(generation),
        "CouponCount": len(results),
        "Passed": all(result["Passed"] for result in results),
        "Results": results,
    }
    output = (
        args.output.expanduser().resolve()
        if args.output
        else generation / "spatial-heldout-qualification.json"
    )
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(output)
    if not summary["Passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
