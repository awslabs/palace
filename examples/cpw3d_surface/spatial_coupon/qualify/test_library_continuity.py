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


def write_planar_model(root, name, topology, edge_x, per_edge_total, depth=1055.0, points_per_edge=3, **fields):
    """A two-dimensional cross-section model (no Point / Interval Edges): basis points in the
    plane z = 0 next to each edge line x = edge_x[k], the geometry given as the 2D builders
    write it (Separation for pairs, offset Edges for stacks, nothing for the isolated edge)."""
    directory = os.path.join(root, name)
    os.makedirs(directory, exist_ok=True)
    basis, rows, index = [], [], 1
    for k, x in enumerate(edge_x):
        for m in range(points_per_edge):
            basis.append((x + 0.01 * (m - 1), 0.3, 0.0))
            rows.append((index, index, per_edge_total[k] / points_per_edge))
            index += 1
    with open(os.path.join(directory, "basis-points.csv"), "w") as target:
        target.write("x,y,z\n")
        for p in basis:
            target.write(f"{p[0]!r},{p[1]!r},{p[2]!r}\n")
    with open(os.path.join(directory, "domain-response-matrix.csv"), "w") as target:
        target.write("  basis_i,  basis_j,                   Q_ij (J)\n")
        for i, j, q in rows:
            target.write(f" {i:.2e}, {j:.2e}, {q:+.12e}\n")
    return {"Name": name, "Topology": topology, "CouponDepth": depth, **fields,
            "FabricatedMatrix": os.path.join(name, "domain-response-matrix.csv"), "BasisPoints": os.path.join(name, "basis-points.csv")}


def test_planar_cross_section_models_are_gated_through_their_2d_frame():
    """The 2D isolated / pair / stack coupons carry no Point / Interval Edges: their edge lines
    are lifted from the 2D frame (x = 0; -+ s / 2; the offsets) with the CouponDepth as the
    edge length, so the recorded gate applies to real cross-section libraries."""
    with tempfile.TemporaryDirectory() as root:
        iso = write_planar_model(root, "isolated-edge", "IsolatedEdge", [0.0], [8.0])
        pair = write_planar_model(root, "strip-2R", "SameConductorStrip", [-R, R], [8.0, 8.0], Separation=2.0 * R)
        narrow = write_planar_model(root, "strip-1R", "SameConductorStrip", [-0.5 * R, 0.5 * R], [12.0, 12.0], Separation=R)
        stack = write_planar_model(root, "stack-2R", "ParallelEdgeCluster", [0.0, 2.0 * R, 4.0 * R], [8.0, 8.0, 8.16],
                                   Edges=[{"Offset": 0.0, "GapDirection": 1, "Conductor": 1}, {"Offset": 2.0 * R, "GapDirection": -1, "Conductor": 2},
                                          {"Offset": 4.0 * R, "GapDirection": 1, "Conductor": 2}])
        path, library = library_with(root, [iso, pair, narrow, stack])
        record = LC.evaluate(path, gates(), library)
        assert [m["Model"] for m in record["Models"]] == ["strip-2R", "stack-2R"]
        assert abs(record["Isolated"]["ResponsePerLength"] - 8.0 / 1055.0) < 1e-15
        assert record["Models"][0]["Status"] == "PASS"
        assert record["Models"][1]["Status"] == "FAIL"   # 2 % on the third edge
        assert record["Verdict"] == "Failed"
        assert abs(record["WorstRelativeOffset"] - 0.02) < 1e-9


def write_point_model(root, name, topology, points, diagonal, depth=1055.0, **fields):
    """A planar model with explicit basis points (x, y) and diagonal entries."""
    directory = os.path.join(root, name)
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "basis-points.csv"), "w") as target:
        target.write("x,y,z\n")
        for x, y in points:
            target.write(f"{x!r},{y!r},0.0\n")
    with open(os.path.join(directory, "domain-response-matrix.csv"), "w") as target:
        target.write("  basis_i,  basis_j,                   Q_ij (J)\n")
        for index, q in enumerate(diagonal, start=1):
            target.write(f" {index:.2e}, {index:.2e}, {q:+.12e}\n")
    return {"Name": name, "Topology": topology, "CouponDepth": depth, **fields,
            "FabricatedMatrix": os.path.join(name, "domain-response-matrix.csv"), "BasisPoints": os.path.join(name, "basis-points.csv")}


def box_top(x0, gap, spacing, count):
    """Basis points along the top of a +-R box around the edge at x0 (u toward the gap)."""
    return [(x0 + gap * (-R + spacing * m), R) for m in range(count)]


def test_hat_count_artefact_passes_per_basis_function_and_records_the_former_failure():
    """Decision 117(3): a threshold pair whose edges carry FEWER basis functions than the
    isolated coupon (its box is shared between the edges) but whose basis functions at
    matched positions carry the isolated energies passes; the former per-length sum
    (a hat-count ratio) is recorded as a failure it no longer decides."""
    with tempfile.TemporaryDirectory() as root:
        spacing = 0.125 * R
        iso_points = box_top(0.0, 1.0, spacing, 17)                          # u = -R .. +R
        iso_q = [1.0, 1.5] + [2.0] * 13 + [1.5, 1.0]                          # corner hats carry less
        iso = write_point_model(root, "isolated-edge", "IsolatedEdge", iso_points, iso_q)
        # The pair at exactly 2R shares the box: each edge keeps the hats on its own side up
        # to the midpoint (u < +R for the first edge; +R is its neighbour's -R).
        pair_points = box_top(-R, 1.0, spacing, 17)[:-1] + box_top(R, -1.0, spacing, 17)[:-1]
        pair_q = ([1.0, 1.5] + [2.0] * 13 + [1.5]) * 2
        pair = write_point_model(root, "gap-2R", "SameConductorGap", pair_points, pair_q, Separation=2.0 * R)
        path, library = library_with(root, [iso, pair])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "Passed"
        entry = record["Models"][0]
        assert entry["Status"] == "PASS"
        assert all(abs(o) < 1e-12 for o in entry["RelativeOffsets"])
        for edge in entry["PerEdge"]:
            # u = -R, -0.875 R and +0.875 R lie in the 0.25 R corner region of the isolated box.
            assert edge["Matched"] == 13 and edge["CornerExcluded"] == 3 and edge["Unmatched"] == 0
            assert abs(edge["MedianPerHatOffset"]) < 1e-12 and abs(edge["WorstPerHatOffset"]) < 1e-12
        former = entry["PreviousDefinition"]
        assert former["Status"] == "FAIL"
        assert all(abs(o - (30.0 - 31.0) / 31.0) < 1e-12 for o in former["RelativeOffsets"])
        assert record["Isolated"]["BasisFunctions"] == 17
        assert abs(record["Isolated"]["HatSpacing"] - spacing) < 1e-12


def test_corner_region_hats_are_excluded_and_unmatched_hats_are_reported():
    with tempfile.TemporaryDirectory() as root:
        spacing = 0.125 * R
        iso = write_point_model(root, "isolated-edge", "IsolatedEdge", box_top(0.0, 1.0, spacing, 17), [1.0, 1.5] + [2.0] * 13 + [1.5, 1.0])
        # The strip's corner-region hats are far off (the box-truncation effect of the isolated
        # coupon), one hat sits at a position the isolated coupon does not have (unmatched);
        # the comparable hats agree.
        # The strip's edges sit at -+R with the gaps outward; the shared midpoint hat (x = 0)
        # is written once and belongs to the first edge.
        first = box_top(-R, -1.0, spacing, 17)[:-1] + [(-R - 1.5 * R, 0.0)]
        second = box_top(R, 1.0, spacing, 17)[1:-1]
        edge_q = [1.6, 2.6] + [2.0] * 13 + [7.0]
        strip = write_point_model(root, "strip-2R", "SameConductorStrip", first + second, edge_q + [9.9] + edge_q[1:], Separation=2.0 * R)
        path, library = library_with(root, [iso, strip])
        record = LC.evaluate(path, gates(), library)
        entry = record["Models"][0]
        assert record["Verdict"] == "Passed" and entry["Status"] == "PASS"
        assert entry["PerEdge"][0]["CornerExcluded"] == 3 and entry["PerEdge"][0]["Unmatched"] == 1
        assert entry["PerEdge"][0]["Matched"] == 13 and entry["PerEdge"][1]["Matched"] == 13
        assert entry["PerEdge"][1]["CornerExcluded"] == 2 and entry["PerEdge"][1]["Unmatched"] == 0
        # One comparable hat off by 20 % (1.5 % of the edge's matched energy) fails.
        moved = [1.6, 2.6] + [2.0] * 12 + [2.4] + [7.0]
        strip = write_point_model(root, "strip-2R", "SameConductorStrip", first + second, moved + [9.9] + edge_q[1:], Separation=2.0 * R)
        path, library = library_with(root, [iso, strip])
        record = LC.evaluate(path, gates(), library)
        entry = record["Models"][0]
        assert record["Verdict"] == "Failed" and entry["Status"] == "FAIL"
        assert abs(entry["RelativeOffsets"][0] - 0.4 / 26.0) < 1e-12
        assert abs(entry["PerEdge"][0]["WorstPerHatOffset"] - 0.2) < 1e-12
        assert abs(entry["RelativeOffsets"][1]) < 1e-12
        assert abs(record["WorstRelativeOffset"] - 0.4 / 26.0) < 1e-12


def test_candidate_without_matched_hats_is_not_comparable():
    with tempfile.TemporaryDirectory() as root:
        spacing = 0.125 * R
        iso = write_point_model(root, "isolated-edge", "IsolatedEdge", box_top(0.0, 1.0, spacing, 17), [2.0] * 17)
        # Hats half a box above the isolated coupon's: no partner within half a spacing.
        points = [(x, y + 0.5 * R) for x, y in box_top(-R, 1.0, spacing, 17)[:-1] + box_top(R, -1.0, spacing, 17)[:-1]]
        pair = write_point_model(root, "gap-2R", "SameConductorGap", points, [2.0] * 32, Separation=2.0 * R)
        path, library = library_with(root, [iso, pair])
        record = LC.evaluate(path, gates(), library)
        assert record["Verdict"] == "NotApplicable"
        assert record["Models"][0]["Status"] == "NotComparable"
        assert all(edge["Matched"] == 0 for edge in record["Models"][0]["PerEdge"])
