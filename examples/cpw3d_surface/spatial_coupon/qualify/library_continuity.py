#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Library continuity gate (decision 82(4), 2026-09-25; per-basis-function normalisation by USER
decision 117(3), 2026-09-28): the strict interaction rule (two edges interact iff their
separation is below 2R) is only consistent if a pair / stack model at a separation of 2R
responds like two isolated edges. For every model of a process library whose edges are all at
least 2R (1 - band) apart from their neighbours (a pair at the threshold, a stack whose
consecutive separations reach 2R) the response of every basis function assigned to an edge is
compared with the isolated-edge model's basis function at the SAME position relative to the
edge (matched in the edge's local frame: offset toward the gap, height above the process
plane, longitudinal offset; within MatchToleranceOverSpacing x the isolated hat spacing). Two
gates per edge (USER decision 2026-09-28 on the decision-117 review MAJOR-2): (a) the
energy-weighted per-basis-function offset (sum of the matched candidate Q_ii over the sum of
their isolated partners' Q_ii, minus one) within MaximumRelativeOffset (1 %), AND (b) every
matched basis function's own offset within PerBasisFunctionLimit (5 %; a single hat carries
the mesh noise of two independently meshed coupons, +-1-2 %, so the per-hat limit is looser
than the aggregate); the verdict fails if any matched hat exceeds (b). An edge with ZERO
matched basis functions is a FAILURE (a basis that matches nothing cannot pass silently; it
was NotApplicable before the decision). The per-hat offsets (median, worst) are recorded.
Basis functions within CornerExclusionOverR x R of a corner of the isolated
coupon's box (|u| and |v| both above R - c) are not compared: the isolated box's truncation
changes their support (the hats next to a box corner differ by 10-60 % between an isolated
and a stack coupon on the same physics), and candidate hats without an isolated partner at
their position are reported as unmatched, never gated. The recorded former definition (sum
of Q_ii over the hats nearest an edge divided by the edge length, an artefact of the hat
count per edge: -23 / -51 / -24 % at the 1.985 R stack while the per-hat response agrees
within 1 %) is evaluated next to it as PreviousDefinition and does not decide the verdict.
A library without a pair / stack at the threshold gets the verdict NotApplicable (nothing to
compare); a library without an isolated-edge model cannot be checked (NotApplicable, reason
recorded).

Model data (process-library.json): `Edges` [{Point, GapDirection, Interval, ...}] in the
model's frame (R = the library's MatchingRadius in the same units), the domain response
matrix (`FabricatedMatrix`, CSV basis_i, basis_j, Q_ij) and `BasisPoints` (CSV x, y, z of every
basis function, in the model's frame) resolved relative to the library file's directory or its
`Root`. Basis points are assigned to the nearest edge line (the edge's Point along its
longitudinal axis = the process normal x gap direction). The two-dimensional cross-section
coupons (isolated edge, pairs by Separation, stacks by their offset Edges) carry no Point /
Interval geometry: their edge lines are lifted from the 2D frame (`planar_edges`), the edge
length is the model's CouponDepth (Palace's 2D normalisation, one depth per library).

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


def local_coordinates(point, frame):
    """A basis point in the edge's local frame: (u, v, w) = offset toward the gap, height
    along the process normal, longitudinal offset from the edge's Point."""
    d = tuple(point[k] - frame["Point"][k] for k in range(3))
    g = frame["Gap"]
    norm = math.sqrt(sum(v * v for v in g)) or 1.0
    u = sum(d[k] * g[k] for k in range(3)) / norm
    w = sum(d[k] * frame["Axis"][k] for k in range(3))
    n = cross(frame["Gap"], frame["Axis"])
    n_norm = math.sqrt(sum(v * v for v in n)) or 1.0
    v = sum(d[k] * n[k] for k in range(3)) / n_norm
    return (u, v, w)


def hat_spacing(points):
    """The smallest distance between two distinct basis points (the hat spacing of a
    perimeter basis); None with fewer than two distinct points."""
    best = None
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            d = math.dist(points[i], points[j])
            if d > 0.0 and (best is None or d < best):
                best = d
    return best


# The two-dimensional cross-section coupons (examples/cpw2d: mesh_edge_coupon.jl,
# mesh_edge_pair_coupon.jl, mesh_edge_cluster_coupon.jl) share one frame: the process
# plane is y = 0 with the process normal +y, the edges run along z (the translation-invariant
# direction) and lie at x = 0 (isolated), x = -+ separation / 2 (pairs) and x = the offsets
# (stacks); every basis point has z = 0 and a per-edge response is per CouponDepth of the
# translation-invariant direction (Palace's 2D normalisation; every model of one library uses
# the same depth, so the per-length comparison is consistent).
PLANAR_PROCESS_NORMAL = (0.0, 1.0, 0.0)


def planar_edges(model, radius):
    """The edge list of a two-dimensional model without stored Edges (isolated / pair), or a
    stack's offset edges given as {Offset, GapDirection, Conductor} lifted to the planar
    frame; None when the model is not a 2D cross-section coupon."""
    depth = model.get("CouponDepth")
    if depth is None:
        return None
    interval = [-0.5 * float(depth), 0.5 * float(depth)]

    def edge(x, gap, conductor=1):
        return {"Point": [x, 0.0, 0.0], "GapDirection": [float(gap), 0.0, 0.0], "ProcessNormal": list(PLANAR_PROCESS_NORMAL),
                "Interval": interval, "Conductor": conductor}

    topology = model.get("Topology")
    edges = model.get("Edges")
    if edges and all("Point" in e for e in edges):
        return None   # a three-dimensional model: its Edges are complete
    if topology in ISOLATED_TOPOLOGIES and not edges:
        return [edge(0.0, 1.0)]
    if topology in ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip") and not edges:
        separation = model.get("Separation")
        if separation is None and model.get("Signature"):
            separation = float(model["Signature"]["SeparationOverR"]) * radius
        if separation is None:
            return None
        strip = topology == "SameConductorStrip"
        return [edge(-0.5 * float(separation), -1.0 if strip else 1.0), edge(0.5 * float(separation), 1.0 if strip else -1.0,
                                                                              2 if topology == "DifferentConductorGap" else 1)]
    if topology in ("ParallelEdgeCluster",) and edges and all("Offset" in e for e in edges):
        return [edge(float(e["Offset"]), float(e["GapDirection"]), int(e.get("Conductor", 1))) for e in edges]
    return None


def model_edges(model, radius):
    """A model's edges with Point / GapDirection / ProcessNormal / Interval: its own Edges
    (three-dimensional models) or the planar lift of a 2D cross-section coupon."""
    planar = planar_edges(model, radius)
    if planar is not None:
        return planar
    edges = model.get("Edges") or []
    return edges if edges and all("Point" in e for e in edges) else []


def model_basis(model, library_path, library, radius=None):
    """A model's diagonal response, basis points and edge frames; None when the files are
    missing or the model carries no edges."""
    matrix_path = resolve(library_path, library, model.get("FabricatedMatrix") or model.get("Matrix"))
    basis_path = resolve(library_path, library, model.get("BasisPoints"))
    edges = model_edges(model, radius)
    if matrix_path is None or basis_path is None or not edges:
        return None
    return {"Diagonal": read_matrix_diagonal(matrix_path), "Points": read_basis_points(basis_path),
            "Frames": edge_frames(edges)}


def assign_to_edges(basis):
    """Per basis point (1-based index): the index of the nearest edge line."""
    frames = basis["Frames"]
    return {index: min(range(len(frames)), key=lambda k: distance_to_edge_line(point, frames[k]))
            for index, point in enumerate(basis["Points"], start=1)}


def isolated_reference(basis):
    """The isolated-edge model's basis functions in its edge's local frame: a list of
    (local coordinates, Q_ii) and the hat spacing."""
    frame = basis["Frames"][0]
    table = [(local_coordinates(point, frame), basis["Diagonal"].get(index, 0.0))
             for index, point in enumerate(basis["Points"], start=1)]
    return {"Table": table, "Spacing": hat_spacing(basis["Points"])}


def in_corner_region(local, radius, exclusion):
    """Whether a basis function lies within `exclusion` of a corner of the isolated coupon's
    +-R box in the cross-section (both |u| and |v| above R - exclusion)."""
    u, v, _ = local
    return abs(u) > radius - exclusion and abs(v) > radius - exclusion


def per_basis_function_comparison(model, basis, reference, radius, match_tolerance, corner_exclusion):
    """Per edge (in the model's Edges order): the candidate's basis functions matched to the
    isolated model's at the same local position, the energy-weighted matched offset (gate
    statistic), the per-hat offsets (median, worst) and the counts of corner-excluded and
    unmatched basis functions."""
    spacing = reference["Spacing"]
    if spacing is None:
        return None
    tolerance = match_tolerance * spacing
    exclusion = corner_exclusion * radius
    assignment = assign_to_edges(basis)
    edges = [{"Matched": 0, "CornerExcluded": 0, "Unmatched": 0, "CandidateEnergy": 0.0, "IsolatedEnergy": 0.0,
              "Offsets": []} for _ in basis["Frames"]]
    for index, point in enumerate(basis["Points"], start=1):
        k = assignment[index]
        local = local_coordinates(point, basis["Frames"][k])
        entry = edges[k]
        if in_corner_region(local, radius, exclusion):
            entry["CornerExcluded"] += 1
            continue
        partner = min(reference["Table"], key=lambda row: math.dist(row[0], local))
        if math.dist(partner[0], local) > tolerance:
            entry["Unmatched"] += 1
            continue
        q = basis["Diagonal"].get(index, 0.0)
        entry["Matched"] += 1
        entry["CandidateEnergy"] += q
        entry["IsolatedEnergy"] += partner[1]
        entry["Offsets"].append((q - partner[1]) / abs(partner[1]) if partner[1] else math.inf)
    result = []
    for entry in edges:
        offsets = sorted(entry.pop("Offsets"))
        if entry["Matched"] and entry["IsolatedEnergy"]:
            entry["MatchedEnergyOffset"] = (entry["CandidateEnergy"] - entry["IsolatedEnergy"]) / abs(entry["IsolatedEnergy"])
            middle = len(offsets) // 2
            entry["MedianPerHatOffset"] = offsets[middle] if len(offsets) % 2 else 0.5 * (offsets[middle - 1] + offsets[middle])
            entry["WorstPerHatOffset"] = max(offsets, key=abs)
        else:
            entry["MatchedEnergyOffset"] = None
        result.append(entry)
    return result


def per_edge_response(model, library_path, library, radius=None):
    """The former gate quantity, per edge (in the model's Edges order): summed diagonal
    response of the basis points nearest to the edge line, divided by the edge length;
    None when the files are missing."""
    basis = model_basis(model, library_path, library, radius)
    if basis is None:
        return None
    diagonal, frames = basis["Diagonal"], basis["Frames"]
    totals = [0.0] * len(frames)
    for index, nearest in assign_to_edges(basis).items():
        totals[nearest] += diagonal.get(index, 0.0)
    responses = []
    for total, frame in zip(totals, frames):
        if not frame["Length"] or frame["Length"] <= 0.0:
            return None
        responses.append(total / frame["Length"])
    return responses


def consecutive_separations(model, radius=None):
    """Separations between consecutive edges along the lateral axis of a 2-edge or stack model
    (the edges sorted by their offset along the first edge's gap direction)."""
    frames = edge_frames(model_edges(model, radius))
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
    match_tolerance = float(table["MatchToleranceOverSpacing"])
    corner_exclusion = float(table["CornerExclusionOverR"])
    per_hat_limit = float(table["PerBasisFunctionLimit"])
    record = {"Gate": "LibraryContinuity", "Statement": table["Statement"], "Quantity": table["Quantity"],
              "PreviousQuantity": table["PreviousQuantity"], "MatchingRadius": None,
              "MaximumRelativeOffset": tolerance, "PerBasisFunctionLimit": per_hat_limit, "ThresholdBandOverR": band,
              "MatchToleranceOverSpacing": match_tolerance,
              "CornerExclusionOverR": corner_exclusion, "Isolated": None, "Models": [], "Verdict": None}
    try:
        radius = matching_radius(library, library_path)
    except (ValueError, OSError) as exception:
        record["Verdict"] = "NotApplicable"
        record["Reason"] = str(exception)
        return record
    record["MatchingRadius"] = radius
    isolated = [m for m in library.get("Models", []) if m.get("Topology") in ISOLATED_TOPOLOGIES and model_edges(m, radius)]
    reference = None
    isolated_response = None
    for model in isolated:
        basis = model_basis(model, library_path, library, radius)
        response = per_edge_response(model, library_path, library, radius)
        if basis is None or not response:
            continue
        reference = isolated_reference(basis)
        if reference["Spacing"] is None:
            reference = None
            continue
        isolated_response = response[0]
        record["Isolated"] = {"Model": model["Name"], "BasisFunctions": len(basis["Points"]), "HatSpacing": reference["Spacing"],
                              "ResponsePerLength": isolated_response}
        break
    candidates = []
    for model in library.get("Models", []):
        if model.get("Topology") not in LONGITUDINAL_TOPOLOGIES or len(model_edges(model, radius)) < 2:
            continue
        separations = consecutive_separations(model, radius)
        if separations and all(s >= 2.0 * radius * (1.0 - band) for s in separations):
            candidates.append((model, separations))
    if not candidates:
        record["Verdict"] = "NotApplicable"
        record["Reason"] = "no pair / stack model with every consecutive separation at or above 2R (1 - band)"
        return record
    if reference is None:
        record["Verdict"] = "NotApplicable"
        record["Reason"] = "no IsolatedEdge model with a readable matrix / basis points to compare against"
        record["Candidates"] = [m["Name"] for m, _ in candidates]
        return record
    worst = 0.0
    worst_hat = 0.0
    for model, separations in candidates:
        basis = model_basis(model, library_path, library, radius)
        entry = {"Model": model["Name"], "Topology": model["Topology"], "SeparationsOverR": [s / radius for s in separations]}
        comparison = None if basis is None else per_basis_function_comparison(model, basis, reference, radius, match_tolerance,
                                                                              corner_exclusion)
        if comparison is None:
            entry["Status"] = "Unreadable"
            record["Models"].append(entry)
            continue
        entry["PerEdge"] = comparison
        offsets = [edge["MatchedEnergyOffset"] for edge in comparison]
        entry["RelativeOffsets"] = offsets
        entry["WorstPerHatOffsets"] = [edge.get("WorstPerHatOffset") for edge in comparison]
        if any(o is None for o in offsets):
            # An edge with no matched basis function fails: nothing of it was compared.
            entry["Status"] = "FAIL"
            entry["Reason"] = "an edge has no basis function matched to an isolated basis function at its position (zero matched hats: Failed, not NotApplicable)"
        else:
            aggregate_pass = all(abs(o) <= tolerance for o in offsets)
            hats = [h for h in entry["WorstPerHatOffsets"] if h is not None]
            per_hat_pass = all(abs(h) <= per_hat_limit for h in hats)
            entry["AggregateStatus"] = "PASS" if aggregate_pass else "FAIL"
            entry["PerBasisFunctionStatus"] = "PASS" if per_hat_pass else "FAIL"
            entry["Status"] = "PASS" if aggregate_pass and per_hat_pass else "FAIL"
            worst = max(worst, max(abs(o) for o in offsets))
            worst_hat = max([worst_hat] + [abs(h) for h in hats])
        response = per_edge_response(model, library_path, library, radius)
        if response is not None and isolated_response:
            former = [(r - isolated_response) / abs(isolated_response) for r in response]
            entry["PreviousDefinition"] = {"PerEdgeResponse": response, "RelativeOffsets": former,
                                           "Status": "PASS" if all(abs(o) <= tolerance for o in former) else "FAIL"}
        record["Models"].append(entry)
    statuses = {m["Status"] for m in record["Models"]}
    record["WorstRelativeOffset"] = worst
    record["WorstPerBasisFunctionOffset"] = worst_hat
    record["Verdict"] = "Failed" if "FAIL" in statuses else ("NotApplicable" if not statuses & {"PASS"} else "Passed")
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
