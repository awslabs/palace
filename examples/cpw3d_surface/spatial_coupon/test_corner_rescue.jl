# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# Round-4 lane F fix-up, supervisor decision 635 (the 56e7 diagnostic, PBS 60933): the CORNER
# RESCUE of optimize_required_region! - an invariant corner still above CornerShapeGate after
# its kappa_reg descent and bridging-sliver rounds (today's fail-closed path) takes the F5-B
# reconnection with the trigger widened to every corner-incident cell above the gate. The
# fixtures are frozen (displacement ratio 0: no vertex can move, the stand-in for the
# guard-locked sub-nm wedge whose every direction the edge floor rejects), so the descent
# changes nothing and only the reconnection can act.

# Three cells around the interior edge corner-b: two healthy wedges and one flat 170-degree
# wedge (kappa_reg ~ 9.9 > the gate 5.0) that the 3-2 flip (the edge removal of corner-b,
# ring of 3) replaces by two cells of kappa_reg ~ 2.4 (the edge's ends over the ring
# triangle, which the edge pierces).
function locked_wedge()
    s = 0.01
    points = s .* hcat([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.6, 0.0, 0.5],
                       [-0.1, 0.7, 0.5], [-0.1, -0.7, 0.5])
    tetrahedra = [(1, 2, 3, 4), (1, 2, 4, 5), (1, 2, 5, 3)]
    return points, tetrahedra
end

const WEDGE_SIDES = [(walls=[[1.0, 0.0], [cosd(150.0), sind(150.0)]], metal=[cosd(75.0), sind(75.0)])]

function run_locked(points, tetrahedra, triangles; gate, displacement_ratio=0.0)
    spans = Tuple{Vector{Float64}, Vector{Float64}}[]
    return optimize_required_region!(copy(points), copy(tetrahedra), triangles, [(0.0, 0.0, 0.0)],
                                     0.1, 0.025, spans, 0.0, 2.0, 0.0, 0.05, 4.0, 0.01, 1000.0,
                                     displacement_ratio, 1.0e-9;
                                     corner_kinds=[:invariant], corner_sides=WEDGE_SIDES,
                                     corner_shape_gate=gate)
end

@testset "the frozen wedge fails the gate under the descent alone; the rescue's edge removal brings it under" begin
    points, tetrahedra = locked_wedge()
    kappas = [tetrahedron_regular_condition([points[:, i] for i in cell]) for cell in tetrahedra]
    @test kappas[2] > 5.0 && all(kappas[[1, 3]] .< 2.0)
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0.0 for cell in tetrahedra)
    record, moved, _, reconnection = run_locked(points, tetrahedra, NTuple{3, Int}[]; gate=5.0)
    row = record["CornerMeasures"][1]
    @test row["Kind"] == "Invariant" && row["Passed"] && row["After"] <= 5.0
    # The descent moved nothing (frozen): AfterDescent is the wedge's kappa_reg, the bridging
    # path found no candidate (no wall triangle), so today's path would have failed closed.
    @test isempty(moved) && row["Moves"] == 0
    @test row["AfterDescent"] == kappas[2] && isempty(row["Reconnections"])
    @test row["BridgingSlivers"] == Dict("Before" => 0, "After" => 0, "AboveGateAfter" => 0)
    rescue = row["CornerRescue"]
    @test rescue["Before"] == kappas[2] && rescue["After"] == row["After"] < 2.5
    @test rescue["Rounds"] == 1 && rescue["CellsAboveGateAfter"] == 0
    @test length(rescue["Reconnections"]) == 1
    removal = rescue["Reconnections"][1]
    @test removal["Kind"] == "edge-removal" && removal["Round"] == 1
    @test sort(removal["Element"]) == [1, 2]
    @test removal["Before"] == kappas[2] && removal["After"] == row["After"]
    @test removal["ReplacedCells"] == 3 && removal["AddedCells"] == 2
    @test occursin("EVERY corner-incident cell", rescue["Rule"]) && occursin("decision 635", rescue["Rule"])
    @test record["CornerRescues"] == 1 && record["CornerRescueReconnections"] == 1
    @test record["CornerRescueRule"] == CORNER_RESCUE_RULE == rescue["Rule"]
    @test record["CornerReconnections"] == 0
    @test record["ReconnectionReplacedCells"] == 3 && record["ReconnectionAddedCells"] == 2
    @test record["CornerAspectsAfter"] == [row["After"]] && record["CornerAspectsBefore"] == [kappas[2]]
    @test length(reconnection.replaced) == 3 && length(reconnection.added) == 2
    # The two new cells are the edge's ends (the corner, b) over the ring triangle (3, 4, 5).
    @test sort([only(v for v in cell if v in (1, 2)) for (_, cell) in reconnection.added]) == [1, 2]
    @test all(issubset((3, 4, 5), cell) for (_, cell) in reconnection.added)
end

@testset "a wedge whose faces and edges are all surface has no reconnection: still fails closed" begin
    points, tetrahedra = locked_wedge()
    triangles = [(1, 2, 4), (1, 2, 5), (1, 4, 5), (2, 4, 5)]
    @test_throws ErrorException run_locked(points, tetrahedra, triangles; gate=5.0)
    try
        run_locked(points, tetrahedra, triangles; gate=5.0)
    catch err
        @test occursin("exceeds CornerShapeGate", err.msg)
    end
end

@testset "the rescue runs only on the fail-closed path: a corner under the gate carries no rescue record" begin
    points, tetrahedra = locked_wedge()
    record, _ = run_locked(points, tetrahedra, NTuple{3, Int}[]; gate=12.0)
    row = record["CornerMeasures"][1]
    @test row["Passed"] && row["After"] > 5.0 && !haskey(row, "CornerRescue")
    @test !haskey(record, "CornerRescues") && !haskey(record, "CornerRescueRule") &&
          !haskey(record, "CornerRescueReconnections")
    @test record["ReconnectionReplacedCells"] == 0 && record["ReconnectionAddedCells"] == 0
    # The production-size fan of test_invariant_corners.jl (a corner that descends under the
    # target by moves alone): no rescue, the record keys as before.
    ring = 6; radius = 0.025
    fan = zeros(3, ring + 2)
    for i in 1:ring
        angle = deg2rad(150.0) * (i - 1) / (ring - 1)
        fan[:, 1 + i] = [radius * cos(angle), radius * sin(angle), 0.0]
    end
    fan[:, ring + 2] = [0.4 * radius * cosd(75.0), 0.4 * radius * sind(75.0), 0.3 * radius]
    cells = [(1, 1 + i, 2 + i, ring + 2) for i in 1:(ring - 1)]
    faces = [(1, 1 + i, 2 + i) for i in 1:(ring - 1)]
    record, _ = run_locked(fan, cells, faces; gate=5.0, displacement_ratio=0.75)
    row = record["CornerMeasures"][1]
    @test row["Passed"] && row["Moves"] > 0 && !haskey(row, "CornerRescue")
    @test !haskey(record, "CornerRescues")
end
