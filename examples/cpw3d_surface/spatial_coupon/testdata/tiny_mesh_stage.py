#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Deterministic stage executable for bounded mesh-DAG regression tests."""
import argparse
import hashlib
import json
from pathlib import Path

import meshio


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rewrite(source, output, binary):
    meshio.write(output, meshio.read(source), file_format="gmsh22", binary=binary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("metric", "adapt", "restore", "publish"))
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--metric", type=Path)
    parser.add_argument("--mmg-seed", type=Path)
    parser.add_argument("--pins", type=Path)
    parser.add_argument("--recipe", type=Path)
    parser.add_argument("--fixed-triangles", type=Path)
    parser.add_argument("--ownership", type=Path)
    args = parser.parse_args()
    if args.stage == "metric":
        if (args.mmg_seed is None or args.pins is None or args.recipe is None or
                args.fixed_triangles is None):
            parser.error("metric requires MMG seed, pins, fixed triangles, and recipe")
        args.output.write_text(json.dumps({"SeedSHA256": digest(args.source),
                                          "Metric": [1.0, 0.0, 1.0]}) + "\n")
        rewrite(args.source, args.mmg_seed, False)
        args.pins.write_text("1\n")
        args.fixed_triangles.write_text("1\n")
        args.recipe.write_text(json.dumps({"SeedSHA256": digest(args.source)}) + "\n")
    elif args.stage == "adapt":
        if (args.metric is None or not args.metric.is_file() or
                args.pins is None or not args.pins.is_file() or
                args.fixed_triangles is None or not args.fixed_triangles.is_file()):
            parser.error("adapt requires metric, pins, and fixed triangles")
        rewrite(args.source, args.output, True)
    elif args.stage == "restore":
        if args.recipe is None or not args.recipe.is_file():
            parser.error("restore requires --recipe")
        rewrite(args.source, args.output, False)
    else:
        if args.ownership is None:
            parser.error("publish requires --ownership")
        rewrite(args.source, args.output, True)
        args.ownership.write_text(
            "attribute,elements,area,ambiguous_area,ambiguous_fraction,"
            "unresolved_elements,unresolved_area,unresolved_fraction\n"
            "1,7,1,0,0,0,0,0\n2,3,1,0,0,0,0,0\n3,1,1,0,0,0,0,0\n")


if __name__ == "__main__":
    main()
