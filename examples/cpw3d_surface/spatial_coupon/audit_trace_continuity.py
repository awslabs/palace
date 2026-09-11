#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Check triangle trace values at all geometric vertices, including T-junctions.

Matching-face meshes may have collinear boundary nodes omitted by a cap fan.
Checking shared vertex IDs alone misses the resulting side/cap trace jumps.
This diagnostic never modifies sources or excludes any source from generation.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def read_trace(path):
    triangles = defaultdict(list)
    with Path(path).open(newline="") as stream:
        for row in csv.DictReader(stream):
            triangles[int(row["triangle"])].append(
                (tuple(float(row[k]) for k in ("x", "y", "z")), float(row["V"]))
            )
    if not triangles or any(len(tri) != 3 for tri in triangles.values()):
        raise ValueError("Trace must contain three vertices per triangle")
    return dict(triangles)


def audit(triangles, tolerance=1e-10):
    values = defaultdict(list)
    for tri in triangles.values():
        for point, value in tri:
            values[point].append(value)
    points = np.asarray(list(values))
    vertex_values = np.asarray([entries[0] for entries in values.values()])
    mismatch = max(max(entries) - min(entries) for entries in values.values())
    worst = None
    violations = 0
    for index, tri in sorted(triangles.items()):
        xyz = np.asarray([p for p, _ in tri])
        a, b, c = xyz
        ab, ac = b - a, c - a
        normal = np.cross(ab, ac)
        length = np.linalg.norm(normal)
        if length <= 0:
            raise ValueError("Degenerate trace triangle")
        distance = np.abs((points - a) @ normal) / length
        candidates = np.flatnonzero(distance <= tolerance)
        # SVD least squares avoids cancellation for thin cap triangles.
        uv = np.linalg.lstsq(np.column_stack((ab, ac)), (points[candidates] - a).T, rcond=None)[0].T
        barycentric = np.column_stack((1 - uv.sum(axis=1), uv))
        inside = (barycentric.min(axis=1) >= -tolerance) & (barycentric.max(axis=1) <= 1 + tolerance)
        for row in np.flatnonzero(inside):
            vertex = candidates[row]
            value = barycentric[row] @ np.asarray([v for _, v in tri])
            error = float(abs(vertex_values[vertex] - value))
            if error > tolerance:
                violations += 1
                if worst is None or error > worst["AbsoluteJump"]:
                    worst = {"Point": points[vertex].tolist(), "Triangle": index,
                             "VertexValue": float(vertex_values[vertex]),
                             "TriangleValue": float(value), "AbsoluteJump": error}
    return {"ContinuousAtAllVertices": violations == 0 and mismatch <= tolerance,
            "Violations": violations, "SharedVertexValueMismatch": mismatch,
            "Worst": worst}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--model-library", type=Path)
    parser.add_argument("--model")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refuse to overwrite trace audit")
    config = json.loads(args.config.read_text())
    zero = set()
    if args.model_library:
        model = next(m for m in json.loads(args.model_library.read_text())["Models"] if m["Name"] == args.model)
        zero = set(model.get("ZeroTraceIndices", []))
    results = []
    for source in config["Boundaries"]["PrescribedPotential"]:
        path = Path(source["DataFile"])
        if not path.is_absolute():
            path = args.config.parent / path
        with path.open() as stream:
            header = stream.readline()
        if "triangle" not in header:
            raise ValueError("This diagnostic requires explicit triangle traces")
        result = audit(read_trace(path))
        results.append({"Index": source["Index"], "DataFile": str(path),
                        "RuntimeZeroTraceIndex": source["Index"] in zero, **result})
    report = {"Config": str(args.config), "Scope": "Declared triangle trace continuity; not a PDE or matrix accuracy test",
              "Sources": results, "DiscontinuousSources": sum(not r["ContinuousAtAllVertices"] for r in results),
              "DiscontinuousRuntimeActiveSources": sum(not r["ContinuousAtAllVertices"] and not r["RuntimeZeroTraceIndex"] for r in results)}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "Sources"}, indent=2))


if __name__ == "__main__":
    main()
