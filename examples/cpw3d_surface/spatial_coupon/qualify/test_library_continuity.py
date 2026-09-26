# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Synthetic library data for the continuity gate (decision 82(4)): an isolated-edge model and
pair / stack models at separations below, at and above 2R with per-edge responses equal to,
or off from, the isolated edge's. Run: python3 -m pytest test_library_continuity.py"""
import json
import os
import tempfile

import library_continuity as LC

R = 2.1


def write_model(root, name, topology, edge_offsets, per_edge_total, length=4.0, points_per_edge=3):
    """A model whose edges lie at the given lateral offsets (x) along z in [-L/2, L/2], with
    `points_per_edge` basis points per edge on the edge line and diagonal entries summing to
    per_edge_total[k] for edge k."""
    directory = os.path.join(root, name)
    os.makedirs(directory, exist_ok=True)
    basis, rows = [], []
    index = 1
    for k, x in enumerate(edge_offsets):
        for m in range(points_per_edge):
            z = -0.5 * length + length * (m + 0.5) / points_per_edge
            basis.append((x, 0.0, z))
            rows.append((index, index, per_edge_total[k] / points_per_edge))
            # An off-diagonal entry (ignored by the gate).
            if index > 1:
                rows.append((index, index - 1, -0.1 * per_edge_total[k]))
            index += 1
    with open(os.path.join(directory, "basis-points.csv"), "w") as target:
        target.write("x,y,z\n")
        for p in basis:
            target.write(f"{p[0]!r},{p[1]!r},{p[2]!r}\n")
    with open(os.path.join(directory, "domain-response-matrix.csv"), "w") as target:
        target.write("  basis_i,  basis_j,                   Q_ij (J)\n")
        for i, j, q in rows:
            target.write(f" {i:.2e}, {j:.2e}, {q:+.12e}\n")
    edges = []
    for k, x in enumerate(edge_offsets):
        sign = -1.0 if k % 2 == 0 else 1.0  # alternating gap directions (a gap between edges 0 and 1)
        edges.append({"Point": [x, 0.0, 0.0], "GapDirection": [sign, 0.0, 0.0], "ProcessNormal": [0.0, 1.0, 0.0],
                      "Interval": [-0.5 * length, 0.5 * length], "Conductor": 1, "BoundaryCondition": {"Type": "PEC"}})
    return {"Name": name, "Topology": topology, "Edges": edges,
            "FabricatedMatrix": os.path.join(name, "domain-response-matrix.csv"), "BasisPoints": os.path.join(name, "basis-points.csv")}


def library_with(root, models):
    library = {"Version": 1, "Name": "synthetic", "MatchingRadius": R, "Models": models}
    path = os.path.join(root, "process-library.json")
    with open(path, "w") as target:
        json.dump(library, target)
    return path, library


def gates():
    with open(LC.GATES_FILE) as source:
        return json.load(source)


def test_pair_at_2r_equal_to_isolated_passes():
    with tempfile.TemporaryDirectory() as root:
        iso = write_model(root, "isolated", "IsolatedEdge", [0.0], [8.0])
        pair = write_model(root, "pair-2R", "SameConductorGap", [0.0, 2.0 * R], [8.0, 8.0])
        below = write_model(root, "pair-1R", "SameConductorGap", [0.0, R], [12.0, 12.0])  # interacting: not gated
        path, library = library_with(root, [iso, pair, below])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "Passed"
        assert [m["Model"] for m in record["Models"]] == ["pair-2R"]
        assert record["Models"][0]["Status"] == "PASS"
        assert abs(record["Isolated"]["ResponsePerLength"] - 2.0) < 1e-12


def test_pair_at_2r_off_by_more_than_tolerance_fails():
    with tempfile.TemporaryDirectory() as root:
        iso = write_model(root, "isolated", "IsolatedEdge", [0.0], [8.0])
        pair = write_model(root, "pair-2R", "SameConductorGap", [0.0, 2.0 * R], [8.0, 8.4])  # 5 % on one edge
        path, library = library_with(root, [iso, pair])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "Failed"
        assert record["Models"][0]["Status"] == "FAIL"
        assert abs(record["WorstRelativeOffset"] - 0.05) < 1e-9


def test_stack_with_one_narrow_separation_is_not_gated():
    with tempfile.TemporaryDirectory() as root:
        iso = write_model(root, "isolated", "IsolatedEdge", [0.0], [8.0])
        stack = write_model(root, "stack", "ParallelEdgeCluster", [0.0, 2.0 * R, 3.0 * R], [8.0, 8.0, 20.0])
        path, library = library_with(root, [iso, stack])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "NotApplicable"
        wide = write_model(root, "stack-wide", "ParallelEdgeCluster", [0.0, 2.0 * R, 4.1 * R], [8.0, 8.0, 8.0])
        path, library = library_with(root, [iso, stack, wide])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "Passed"
        assert [m["Model"] for m in record["Models"]] == ["stack-wide"]


def test_band_admits_a_separation_just_below_2r():
    with tempfile.TemporaryDirectory() as root:
        iso = write_model(root, "isolated", "IsolatedEdge", [0.0], [8.0])
        band = gates()["Gates"]["LibraryContinuity"]["ThresholdBandOverR"]
        pair = write_model(root, "pair-band", "SameConductorGap", [0.0, 2.0 * R * (1.0 - 0.5 * band)], [8.0, 8.0])
        path, library = library_with(root, [iso, pair])
        assert LC.evaluate(path, gates(), library)["Verdict"] == "Passed"


def test_without_isolated_model_is_not_applicable():
    with tempfile.TemporaryDirectory() as root:
        pair = write_model(root, "pair-2R", "SameConductorGap", [0.0, 2.0 * R], [8.0, 8.0])
        path, library = library_with(root, [pair])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "NotApplicable"
        assert record["Candidates"] == ["pair-2R"]
