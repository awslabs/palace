#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Create separate continuous-cap trace inputs, preserving all retained nodal values.

This is a source-definition correction, NOT a byte-identical source migration.
It never overwrites the retained configuration, traces, mesh, or response library.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from audit_trace_continuity import audit, read_trace
from generate_spatial_response import cap_ring


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def repaired_topology(traces):
    points = sorted({p for tri in traces.values() for p, _ in tri})
    array = np.asarray(points)
    index = {point: i for i, point in enumerate(points)}
    zmin, zmax = array[:, 2].min(), array[:, 2].max()
    # Keep every original side triangle and its orientation. Only replace caps.
    triangles = [tuple(index[p] for p, _ in tri) for tri in traces.values()
                 if not (all(p[2] == zmin for p, _ in tri) or
                         all(p[2] == zmax for p, _ in tri))]
    for z, reverse in ((zmin, True), (zmax, False)):
        vertices = np.flatnonzero(array[:, 2] == z)
        xy = array[vertices, :2]
        low, high = xy.min(axis=0), xy.max(axis=0)
        tolerance = 1e-10 * max(1., np.max(high - low))
        if not np.all(np.any((np.abs(xy - low) <= tolerance) |
                             (np.abs(xy - high) <= tolerance), axis=1)):
            raise ValueError("Expected a rectangular cap ring without interior vertices")
        center = xy.mean(axis=0)
        vertices = vertices[np.argsort(np.arctan2(xy[:, 1] - center[1], xy[:, 0] - center[0]))]
        cap = []
        cap_ring(cap, array[vertices], 0, len(vertices), reverse)
        triangles.extend(tuple(int(vertices[i]) for i in tri) for tri in cap)
    return points, triangles


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refuse to reuse a trace-repair output directory")
    config = json.loads(args.config.read_text())
    records = []
    points = triangles = original_topology = None
    generated = []
    for source in config["Boundaries"]["PrescribedPotential"]:
        path = Path(source["DataFile"])
        if not path.is_absolute():
            path = args.config.parent / path
        traces = read_trace(path)
        topology = [[p for p, _ in tri] for tri in traces.values()]
        if original_topology is None:
            original_topology = topology
            points, triangles = repaired_topology(traces)
        elif topology != original_topology:
            raise ValueError("Source trace topologies differ")
        values = {}
        for tri in traces.values():
            for point, value in tri:
                if point in values and abs(values[point] - value) > 1e-12:
                    raise ValueError("Conflicting values at a retained node")
                values[point] = value
        repaired = {i + 1: [(points[v], values[points[v]]) for v in tri]
                    for i, tri in enumerate(triangles)}
        result = audit(repaired)
        if not result["ContinuousAtAllVertices"]:
            raise ValueError(f"Corrected source {source['Index']} remains discontinuous: {result}")
        generated.append((source, repaired))
        records.append({"Index": source["Index"], "Original": str(path),
                        "OriginalSHA256": sha(path), "ContinuousAtAllVertices": True})
    args.output.mkdir(parents=True)
    for (source, repaired), record in zip(generated, records):
        path = args.output / f"source-{source['Index']:04d}.csv"
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(("x", "y", "z", "V", "triangle"))
            for index, tri in repaired.items():
                for point, value in tri:
                    writer.writerow((*[f"{x:.16e}" for x in point], f"{value:.16e}", index))
        record.update(Corrected=str(path.resolve()), CorrectedSHA256=sha(path))
        source["DataFile"] = str(path.resolve())
    (args.output / "config.json").write_text(json.dumps(config, indent=2) + "\n")
    report = {"OriginalConfig": str(args.config.resolve()), "OriginalConfigSHA256": sha(args.config),
              "SourceDefinitionChanged": True, "BasisNodesAndAllNodalValuesPreserved": True,
              "SourceIndicesAndCountPreserved": True, "MeshUnchanged": True,
              "SolverOrderAndSettingsUnchanged": True, "SourceCount": len(records),
              "Vertices": len(points), "TrianglesBefore": len(original_topology),
              "TrianglesAfter": len(triangles), "Sources": records}
    (args.output / "trace-repair.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "Sources"}, indent=2))


if __name__ == "__main__":
    main()
