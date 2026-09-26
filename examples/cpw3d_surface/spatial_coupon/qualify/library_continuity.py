#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Library continuity gate (decision 82(4), 2026-09-25): the strict interaction rule (two edges
interact iff their separation is below 2R) is only consistent if a pair / stack model at a
separation of 2R responds like two isolated edges. For every model of a process library whose
edges are all at least 2R (1 - band) apart from their neighbours (a pair at the threshold, a
stack whose consecutive separations reach 2R) the per-edge response — the sum of the diagonal
domain-response entries Q_ii of the basis points assigned to the edge, per unit edge length —
must equal the isolated-edge model's per-unit-length response within the recorded tolerance
(qualification-gates.json, LibraryContinuity). A library without a pair / stack at the
threshold gets the verdict NotApplicable (nothing to compare); a library without an
isolated-edge model cannot be checked (NotApplicable, reason recorded).

Model data (process-library.json): `Edges` [{Point, GapDirection, Interval, ...}] in the
model's frame (R = the library's MatchingRadius in the same units), the domain response
matrix (`FabricatedMatrix`, CSV basis_i, basis_j, Q_ij) and `BasisPoints` (CSV x, y, z of every
basis function, in the model's frame) resolved relative to the library file's directory or its
`Root`. Basis points are assigned to the nearest edge line (the edge's Point along its
longitudinal axis = the process normal x gap direction).

    python3 library_continuity.py process-library.json [--gates qualification-gates.json]
"""
import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
GATES_FILE = os.path.join(HERE, "qualification-gates.json")
ISOLATED_TOPOLOGIES = ("IsolatedEdge",)
LONGITUDINAL_TOPOLOGIES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip",
                           "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip",
                           "ParallelEdgeCluster", "CurvedParallelEdgeCluster", "SpatialEdgeCluster")


def read_matrix_diagonal(path):
    """Diagonal Q_ii (1-based basis index -> value) of a domain-response-matrix.csv."""
    diagonal = {}
    with open(path) as source:
        reader = csv.reader(source)
        header = next(reader)
        if len(header) < 3:
            raise ValueError(f"{path}: expected basis_i, basis_j, Q_ij")
        for row in reader:
            if not row or not row[0].strip():
                continue
            i, j, q = int(round(float(row[0]))), int(round(float(row[1]))), float(row[2])
            if i == j:
                diagonal[i] = q
    return diagonal


def read_basis_points(path):
    with open(path) as source:
        reader = csv.reader(source)
        next(reader)
        return [tuple(float(v) for v in row[:3]) for row in reader if row and row[0].strip()]


def resolve(library_path, library, value):
    """A model file path relative to the library's directory or its Root."""
    if value is None:
        return None
    if os.path.isabs(value) and os.path.exists(value):
        return value
    for base in (os.path.dirname(os.path.abspath(library_path)), library.get("Root")):
        if base and os.path.exists(os.path.join(base, value)):
            return os.path.join(base, value)
    return None


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def edge_frames(edges):
    """Per edge: point, unit longitudinal axis (process normal x gap direction), length."""
    frames = []
    for edge in edges:
        p = tuple(float(v) for v in edge["Point"][:3])
        g = tuple(float(v) for v in edge["GapDirection"][:3])
        n = tuple(float(v) for v in edge.get("ProcessNormal", [0.0, 0.0, 1.0])[:3])
        t = cross(n, g)
        norm = math.sqrt(sum(v * v for v in t)) or 1.0
        t = tuple(v / norm for v in t)
        interval = edge.get("Interval")
        length = float(interval[1]) - float(interval[0]) if interval else None
        frames.append({"Point": p, "Axis": t, "Gap": g, "Length": length})
    return frames


def distance_to_edge_line(point, frame):
    d = tuple(point[k] - frame["Point"][k] for k in range(3))
    along = sum(d[k] * frame["Axis"][k] for k in range(3))
    perpendicular = tuple(d[k] - along * frame["Axis"][k] for k in range(3))
    return math.sqrt(sum(v * v for v in perpendicular))


def per_edge_response(model, library_path, library):
    """Per edge (in the model's Edges order): summed diagonal response of the basis points
    nearest to the edge line, divided by the edge length; None when the files are missing."""
    matrix_path = resolve(library_path, library, model.get("FabricatedMatrix") or model.get("Matrix"))
    basis_path = resolve(library_path, library, model.get("BasisPoints"))
    if matrix_path is None or basis_path is None or not model.get("Edges"):
        return None
    diagonal = read_matrix_diagonal(matrix_path)
    points = read_basis_points(basis_path)
    frames = edge_frames(model["Edges"])
    totals = [0.0] * len(frames)
    for index, point in enumerate(points, start=1):
        nearest = min(range(len(frames)), key=lambda k: distance_to_edge_line(point, frames[k]))
        totals[nearest] += diagonal.get(index, 0.0)
    responses = []
    for total, frame in zip(totals, frames):
        if not frame["Length"] or frame["Length"] <= 0.0:
            return None
        responses.append(total / frame["Length"])
    return responses


def consecutive_separations(model):
    """Separations between consecutive edges along the lateral axis of a 2-edge or stack model
    (the edges sorted by their offset along the first edge's gap direction)."""
    frames = edge_frames(model["Edges"])
    if len(frames) < 2:
        return []
    g = frames[0]["Gap"]
    offsets = sorted(sum((f["Point"][k] - frames[0]["Point"][k]) * g[k] for k in range(3)) for f in frames)
    return [b - a for a, b in zip(offsets, offsets[1:])]


def matching_radius(library, library_path):
    """The library's MatchingRadius, or that of the models' source process libraries (a
    qualify-written process-library.json keeps the header of its sources per model)."""
    if "MatchingRadius" in library:
        return float(library["MatchingRadius"])
    radii = set()
    for model in library.get("Models", []):
        source = model.get("SourceProcessLibrary") or {}
        if "MatchingRadius" in source:
            radii.add(float(source["MatchingRadius"]))
        else:
            path = resolve(library_path, library, source.get("Path"))
            if path:
                with open(path) as handle:
                    header = json.load(handle)
                if "MatchingRadius" in header:
                    radii.add(float(header["MatchingRadius"]))
    if len(radii) != 1:
        raise ValueError(f"MatchingRadius unknown or inconsistent across the models' sources: {sorted(radii)}")
    return radii.pop()


def evaluate(library_path, gates=None, library=None):
    """The continuity gate record for a process library."""
    if library is None:
        with open(library_path) as source:
            library = json.load(source)
    if gates is None:
        with open(GATES_FILE) as source:
            gates = json.load(source)
    table = gates["Gates"]["LibraryContinuity"]
    tolerance = float(table["MaximumRelativeOffset"])
    band = float(table["ThresholdBandOverR"])
    record = {"Gate": "LibraryContinuity", "Statement": table["Statement"], "MatchingRadius": None,
              "MaximumRelativeOffset": tolerance, "ThresholdBandOverR": band, "Isolated": None, "Models": [], "Verdict": None}
    try:
        radius = matching_radius(library, library_path)
    except (ValueError, OSError) as exception:
        record["Verdict"] = "NotApplicable"
        record["Reason"] = str(exception)
        return record
    record["MatchingRadius"] = radius
    isolated = [m for m in library.get("Models", []) if m.get("Topology") in ISOLATED_TOPOLOGIES and m.get("Edges")]
    isolated_response = None
    for model in isolated:
        response = per_edge_response(model, library_path, library)
        if response:
            isolated_response = response[0]
            record["Isolated"] = {"Model": model["Name"], "ResponsePerLength": isolated_response}
            break
    candidates = []
    for model in library.get("Models", []):
        if model.get("Topology") not in LONGITUDINAL_TOPOLOGIES or len(model.get("Edges", [])) < 2:
            continue
        separations = consecutive_separations(model)
        if separations and all(s >= 2.0 * radius * (1.0 - band) for s in separations):
            candidates.append((model, separations))
    if not candidates:
        record["Verdict"] = "NotApplicable"
        record["Reason"] = "no pair / stack model with every consecutive separation at or above 2R (1 - band)"
        return record
    if isolated_response is None:
        record["Verdict"] = "NotApplicable"
        record["Reason"] = "no IsolatedEdge model with a readable matrix / basis points to compare against"
        record["Candidates"] = [m["Name"] for m, _ in candidates]
        return record
    worst = 0.0
    for model, separations in candidates:
        response = per_edge_response(model, library_path, library)
        entry = {"Model": model["Name"], "Topology": model["Topology"], "SeparationsOverR": [s / radius for s in separations]}
        if response is None:
            entry["Status"] = "Unreadable"
            record["Models"].append(entry)
            continue
        offsets = [(r - isolated_response) / abs(isolated_response) if isolated_response else math.inf for r in response]
        entry["PerEdgeResponse"] = response
        entry["RelativeOffsets"] = offsets
        entry["Status"] = "PASS" if all(abs(o) <= tolerance for o in offsets) else "FAIL"
        worst = max(worst, max(abs(o) for o in offsets))
        record["Models"].append(entry)
    statuses = {m["Status"] for m in record["Models"]}
    record["WorstRelativeOffset"] = worst
    record["Verdict"] = "Failed" if "FAIL" in statuses else ("NotApplicable" if statuses <= {"Unreadable"} else "Passed")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("library")
    parser.add_argument("--gates", default=GATES_FILE)
    parser.add_argument("--output", help="write the gate record here (default: stdout)")
    args = parser.parse_args(argv)
    with open(args.gates) as source:
        gates = json.load(source)
    record = evaluate(args.library, gates)
    text = json.dumps(record, indent=1)
    if args.output:
        with open(args.output, "w") as target:
            target.write(text + "\n")
    print(text if not args.output else f"{record['Verdict']}: {args.output}")
    return 0 if record["Verdict"] != "Failed" else 1


if __name__ == "__main__":
    sys.exit(main())
