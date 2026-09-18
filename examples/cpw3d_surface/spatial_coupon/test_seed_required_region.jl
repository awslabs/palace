# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# The seed-side required-region optimization (supervisor decision 30) on small
# hand-built tetrahedral fans: the region rule, the movement bases, the bounded
# coordinate descent and its floors.

@testset "required region rule matches the metric stage" begin
    points = hcat([0.01, 0.01, 0.01], [0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.05],
                  [0.2, 0.0, 0.0], [0.3, 0.0, 0.0], [0.2, 0.1, 0.0], [0.2, 0.0, 0.1],
                  [5.0, 0.03, 0.0], [5.2, 0.03, 0.0], [5.1, 0.5, 0.0], [5.1, 0.03, 0.5],
                  [5.0, 0.04, 0.0], [5.2, 0.04, 0.0], [5.1, 0.5, 0.0], [5.1, 0.04, 0.5])
    tetrahedra = [(1, 2, 3, 4), (5, 6, 7, 8), (9, 10, 11, 12), (13, 14, 15, 16)]
    spans = [([1.0, 0.0, 0.0], [9.0, 0.0, 0.0])]
    reach = 0.028 * 1.05 + 0.004
    required, distance = required_region_cells(points, tetrahedra, [(0.0, 0.0, 0.0)], 0.1,
                                               spans, reach)
    @test collect(required) == [true, false, true, false]
    @test distance[9] ≈ 0.03 && distance[13] ≈ 0.04
    @test span_point_distance([5.0, -2.0, 0.05], [9.9, -2.0, 0.0], [0.2, -2.0, 0.0]) ≈ 0.05
    @test span_point_distance([12.0, -2.0, 0.0], [9.9, -2.0, 0.0], [0.2, -2.0, 0.0]) ≈ 2.1
    none, _ = required_region_cells(points, tetrahedra, [(0.0, 0.0, 0.0)], 0.1,
                                    Tuple{Vector{Float64}, Vector{Float64}}[], 0.0)
    @test collect(none) == [true, false, false, false]
end

@testset "surface movement bases: plane, line, fixed" begin
    points = hcat([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
                  [1.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0])
    # Node 1 lies on z = 0, y = 0 and x = 0 (three planes): fixed. Node 2 lies on
    # z = 0 and y = 0: moves along x. Node 5 lies on z = 0 only: moves in the plane.
    triangles = [(1, 2, 5), (1, 5, 3), (1, 2, 6), (1, 6, 4), (1, 3, 7), (1, 7, 4)]
    bases = surface_movement_bases(points, triangles)
    @test size(bases[1], 2) == 0
    @test size(bases[2], 2) == 1 && abs(abs(bases[2][1, 1]) - 1.0) < 1e-12
    @test size(bases[5], 2) == 2 && all(abs.(bases[5][3, :]) .< 1e-12)
    @test !haskey(bases, 8)
end

# A fan of tetrahedra around a fixed corner at the origin: ring vertices on the
# plane z = 0 at radius r, one interior apex above (r/2, 0) whose height sets the
# quality (the corner's three edges become coplanar as the height vanishes).
function corner_fan(apex_height; ring=6, radius=0.025)
    points = zeros(3, ring + 2)
    for i in 1:ring
        angle = 2pi * (i - 1) / ring
        points[:, 1 + i] = [radius * cos(angle), radius * sin(angle), 0.0]
    end
    points[:, ring + 2] = [0.5 * radius, 0.0, apex_height]
    tetrahedra = [(1, 1 + i, 1 + mod(i, ring) + 1, ring + 2) for i in 1:ring]
    triangles = [(1, 1 + i, 1 + mod(i, ring) + 1) for i in 1:ring]
    return points, tetrahedra, triangles
end

@testset "bounded descent lifts a needle apex within the bound and floors" begin
    points, tetrahedra, triangles = corner_fan(0.0002)
    original = copy(points)
    incident = [Int[] for _ in axes(points, 2)]
    for (k, cell) in enumerate(tetrahedra), i in cell
        push!(incident[i], k)
    end
    bases = surface_movement_bases(points, triangles)
    # The corner (1) and the ring lie on z = 0; the corner is additionally pinned by
    # making it fixed (three planes are not present here, so fix it explicitly).
    bases[1] = zeros(3, 0)
    scaled = [tetrahedron_scaled_jacobian([points[:, i] for i in cell]) for cell in tetrahedra]
    @test all(<(0.02), scaled)
    floors = min.(scaled, 0.02)
    bounds = fill(0.75 * 0.025, size(points, 2))
    targets = collect(eachindex(tetrahedra))
    value, moves = optimize_seed_cells!(points, original, tetrahedra, incident, bases, floors,
                                        bounds, targets, :scaled, 0.02)
    @test moves > 0
    @test value >= 0.02
    displacement = [norm(points[:, i] .- original[:, i]) for i in axes(points, 2)]
    @test all(displacement .<= 0.75 * 0.025 * (1 + 1e-12))
    # Only the apex (interior) moved; the ring stayed on z = 0 and the corner is fixed.
    @test displacement[1] == 0.0
    @test all(abs.(points[3, 2:7]) .< 1e-15)
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
    # The aspect objective on the same fan: the maximum incident aspect decreases.
    points2, tetrahedra2, triangles2 = corner_fan(0.002)
    original2 = copy(points2)
    bases2 = surface_movement_bases(points2, triangles2); bases2[1] = zeros(3, 0)
    scaled2 = [tetrahedron_scaled_jacobian([points2[:, i] for i in cell]) for cell in tetrahedra2]
    before = maximum(tetrahedron_aspect([points2[:, i] for i in cell]) for cell in tetrahedra2)
    aspect, _ = optimize_seed_cells!(points2, original2, tetrahedra2, incident, bases2,
                                     min.(scaled2, 0.02), bounds, targets, :aspect, 3.8)
    @test aspect < before
    @test aspect <= 3.8 + 1e-6
    # The guards of decision 34: with an edge floor above every edge of the apex
    # no move may shorten one of its edges (min(original, floor) = original), and
    # with condition ceilings at the cells' own values no touched cell may end
    # with a higher Jacobian condition; the scaled objective is still reached.
    points3, tetrahedra3, triangles3 = corner_fan(0.0002)
    original3 = copy(points3)
    bases3 = surface_movement_bases(points3, triangles3); bases3[1] = zeros(3, 0)
    scaled3 = [tetrahedron_scaled_jacobian([points3[:, i] for i in cell]) for cell in tetrahedra3]
    condition3 = [tetrahedron_aspect([points3[:, i] for i in cell]) for cell in tetrahedra3]
    apex = size(points3, 2)
    lengths_before = [norm(points3[:, j] .- points3[:, apex]) for j in 1:(apex - 1)]
    value3, moves3 = optimize_seed_cells!(points3, original3, tetrahedra3, incident, bases3,
                                          min.(scaled3, 0.02), bounds, targets, :scaled, 0.02;
                                          condition_ceilings=condition3,
                                          edge_floors=fill(1.0, size(points3, 2)))
    @test moves3 > 0 && value3 >= 0.02
    lengths_after = [norm(points3[:, j] .- points3[:, apex]) for j in 1:(apex - 1)]
    @test all(lengths_after .>= lengths_before .* (1 - 1e-9))
    @test all(tetrahedron_aspect([points3[:, i] for i in cell]) <= ceiling * (1 + 1e-9)
              for (cell, ceiling) in zip(tetrahedra3, condition3))
    # A floor below every edge never binds: the unguarded result is reproduced.
    points4, tetrahedra4, triangles4 = corner_fan(0.0002)
    bases4 = surface_movement_bases(points4, triangles4); bases4[1] = zeros(3, 0)
    optimize_seed_cells!(points4, copy(points4), tetrahedra4, incident, bases4,
                         min.(scaled3, 0.02), bounds, targets, :scaled, 0.02;
                         edge_floors=fill(1.0e-6, size(points4, 2)))
    @test points4 == points
end

@testset "descent directions cover every basis combination" begin
    @test length(descent_directions(Matrix{Float64}(I, 3, 3))) == 26
    @test length(descent_directions(Matrix{Float64}(I, 3, 3)[:, 1:2])) == 8
    @test length(descent_directions(Matrix{Float64}(I, 3, 3)[:, 1:1])) == 2
    @test all(abs(norm(d) - 1.0) < 1e-12 for d in descent_directions(Matrix{Float64}(I, 3, 3)))
end

# A corner fan plus a layer needle whose repair carries its free apex across the
# layer reach: the record counts the required set of the FINAL positions.
function fan_with_layer_needle()
    points, tetrahedra, triangles = corner_fan(0.0125)
    n = size(points, 2)
    layer = hcat([1.5, 0.0, 0.030], [1.55, 0.0, 0.030], [1.5, 0.02, 0.030],
                 [1.525, 0.007, 0.0302],
                 [1.5, 0.0, 0.06], [1.55, 0.0, 0.06], [1.5, 0.02, 0.06])
    points = hcat(points, layer)
    needle = (n + 1, n + 2, n + 3, n + 4)
    neighbor = (n + 5, n + 7, n + 6, n + 4)
    return points, vcat(tetrahedra, [needle, neighbor]), triangles, needle, neighbor
end

@testset "required set is recomputed on the moved positions and gated there" begin
    points, tetrahedra, triangles, needle, neighbor = fan_with_layer_needle()
    spans = [([1.0, 0.0, 0.0], [2.0, 0.0, 0.0])]
    edge_size, thickness, zigzag = 0.004, 0.028, 0.05
    reach = thickness * (1.0 + zigzag) + edge_size
    corners = [(0.0, 0.0, 0.0)]
    before, _ = required_region_cells(points, tetrahedra, corners, 0.1, spans, reach)
    # The needle (vertices within reach) and its neighbor (sharing the apex) are
    # required before any move; the neighbor's other vertices are beyond reach.
    @test before[length(tetrahedra) - 1] && before[length(tetrahedra)]
    @test tetrahedron_scaled_jacobian([points[:, i] for i in needle]) < 0.02
    record, moved = optimize_required_region!(points, tetrahedra, triangles, corners, 0.1,
                                              0.025, spans, edge_size, 2.0, thickness, zigzag,
                                              4.0, 0.01, 1000.0, 0.75, 1e-9)
    after, _ = required_region_cells(points, tetrahedra, corners, 0.1, spans, reach)
    @test record["RequiredTetrahedra"] == count(after)
    @test record["RequiredTetrahedraBeforeMoves"] == count(before)
    @test record["RequiredCellsBelowGateAfter"] == 0
    @test record["LayerRequiredReach"] == reach
    @test tetrahedron_scaled_jacobian([points[:, i] for i in needle]) >= 0.02
    @test !isempty(moved)
    # The apex left the reach: the neighbor is no longer required and the record
    # says so (a second recomputation ran because the set changed).
    if !after[length(tetrahedra)]
        @test record["RequiredTetrahedra"] == count(before) - 1
        @test record["RequiredSetRecomputations"] == 2
    else
        @test record["RequiredSetRecomputations"] == 1
    end
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
    # Without any layer the set is the corner ball only and never changes.
    points0, tetrahedra0, triangles0 = corner_fan(0.0125)
    record0, _ = optimize_required_region!(points0, tetrahedra0, triangles0, corners, 0.1, 0.025,
                                           Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0,
                                           0.0, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9)
    @test record0["RequiredTetrahedra"] == length(tetrahedra0) == record0["RequiredTetrahedraBeforeMoves"]
    @test record0["LayerRequiredReach"] === nothing
    @test record0["RequiredSetRecomputations"] == 1
end

# The edge-layer quality rule (decision 32): a 1.1 nm x 2 nm x 250 nm layer sliver
# on the span whose scaled Jacobian is far below the gate by construction (every
# vertex within 2 nm of the span, so the bounded repair moves at most ~2 nm; no
# edge below EdgeSize, so nothing is collapsed).
function fan_with_layer_slab(height)
    points, tetrahedra, triangles = corner_fan(0.0125)
    n = size(points, 2)
    slab = hcat([1.5, 0.0, 0.0], [1.75, 0.0, 0.0], [1.75, 0.002, 0.0], [1.75, 0.002, height])
    return hcat(points, slab), vcat(tetrahedra, [(n + 1, n + 2, n + 3, n + 4)]), triangles
end

function positively_oriented(points, cell)
    return tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 ? cell :
           (cell[1], cell[3], cell[2], cell[4])
end

@testset "edge layer quality rule gates orientation and edge aspect, not the scaled Jacobian" begin
    spans = [([1.0, 0.0, 0.0], [2.0, 0.0, 0.0])]
    corners = [(0.0, 0.0, 0.0)]
    edge_size, thickness, zigzag = 0.001, 0.031, 0.05
    points, tetrahedra, triangles = fan_with_layer_slab(0.0011)
    slab = tetrahedra[end]
    @test tetrahedron_scaled_jacobian([points[:, i] for i in slab]) < 0.01
    @test tetrahedron_edge_aspect([points[:, i] for i in slab]) > 200.0
    # Without the rule the slab is a scaled-Jacobian-gated required cell and the
    # bounded repair (0.75 x ~2 nm) cannot lift it to the gate: fails closed.
    @test_throws ErrorException optimize_required_region!(
        copy(points), copy(tetrahedra), triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9)
    # With the rule the slab is gated by orientation and edge aspect: the aspect
    # repair lifts the apex within the bound and the scaled Jacobian is a diagnostic.
    moved_points = copy(points); moved_cells = copy(tetrahedra)
    record, moved, collapse = optimize_required_region!(
        moved_points, moved_cells, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9; edge_layer_maximum_aspect=200.0)
    rule = record["EdgeLayerQualityRule"]
    @test rule["LayerCells"] == 1 && record["SubSizeEdgeCollapse"]["CollapsedVertices"] == 0
    @test isempty(collapse.cells_removed) && isempty(collapse.cells_remapped) &&
          moved_cells == tetrahedra
    @test rule["MaximumEdgeAspect"] == 200.0 && rule["EdgeAspectTarget"] == 190.0
    @test rule["MaximumEdgeAspectBefore"] > 200.0
    @test rule["MaximumEdgeAspectAfter"] <= 200.0 && rule["CellsAboveBoundAfter"] == 0
    @test rule["MinimumScaledJacobian"] < 0.01
    @test rule["ScaledJacobianRoundoffFloor"] == 1.0e-12
    @test record["ScaledJacobianGateCells"] == record["RequiredTetrahedra"] - 1
    @test record["RequiredCellsBelowGateAfter"] == 0
    @test tetrahedron_edge_aspect([moved_points[:, i] for i in slab]) <= 200.0
    @test tetrahedron_scaled_jacobian([moved_points[:, i] for i in slab]) > 0.0
    @test !isempty(moved)
    # A bound the bounded repair cannot reach fails closed.
    @test_throws ErrorException optimize_required_region!(
        copy(points), copy(tetrahedra), triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9; edge_layer_maximum_aspect=1.5)
    # The rule needs a seeded layer.
    points0, tetrahedra0, triangles0 = corner_fan(0.0125)
    @test_throws ErrorException optimize_required_region!(
        points0, tetrahedra0, triangles0, corners, 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9;
        edge_layer_maximum_aspect=100.0)
    # Edge aspect: unit right tetrahedron sqrt(6); a 50 x 50 x 1 slab ~ 70.7.
    @test tetrahedron_edge_aspect([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]) ≈ sqrt(6)
    @test tetrahedron_edge_aspect([[0.0, 0, 0], [50.0, 0, 0], [0, 50.0, 0], [0, 0, 1.0]]) ≈ 70.7389567 atol=1e-6
    @test tetrahedron_edge_aspect([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [1.0, 1.0, 0]]) == Inf
end

# A seed volume vertex 0.28 nm from a row node (below EdgeSize = the adapter hmin)
# inside the layer: collapsed onto that node; the cell containing both vanishes,
# the other incident cell is remapped, and the survivors satisfy the rule. The
# collapse runs with and without the layer quality rule (decision 34: the aspect-4
# layers do not use the rule); without it the cavity must also keep the replaced
# cells' minimum scaled Jacobian (the guard the flat cell trivially allows).
@testset "sub-EdgeSize seed vertices inside the layer are collapsed onto their neighbour" begin
    spans = [([1.0, 0.0, 0.0], [2.0, 0.0, 0.0])]
    corners = [(0.0, 0.0, 0.0)]
    edge_size, thickness, zigzag = 0.001, 0.031, 0.05
    points, tetrahedra, triangles = corner_fan(0.0125)
    n = size(points, 2)
    extra = hcat([1.5, 0.0, 0.0], [1.55, 0.0, 0.0], [1.55, 0.002, 0.0], [1.5, 0.0002, 0.0002],
                 [1.55, 0.002, 0.003])
    points = hcat(points, extra)
    w, p2, p3, v, q = n + 1, n + 2, n + 3, n + 4, n + 5
    flat = positively_oriented(points, (w, p2, p3, v))
    other = positively_oriented(points, (v, p2, p3, q))
    tetrahedra = vcat(tetrahedra, [flat, other])
    # The slab base is a face (row nodes are surface vertices).
    triangles = vcat(triangles, [(w, p2, p3)])
    cells_before = length(tetrahedra)
    @test tetrahedron_edge_aspect([points[:, i] for i in flat]) > 200.0
    record, moved, collapse = optimize_required_region!(
        points, tetrahedra, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9; edge_layer_maximum_aspect=200.0)
    rule = record["EdgeLayerQualityRule"]
    collapsed = record["SubSizeEdgeCollapse"]
    @test collapsed["CollapseSize"] == edge_size && collapsed["CollapseThreshold"] == 0.5 * edge_size
    @test collapsed["ScaledJacobianGate"] == 0.0
    @test collapsed["CollapsedVertices"] == 1 && collapsed["CollapsedCells"] == 1 &&
          collapsed["RemappedCells"] == 1 && collapsed["CollapsedInteriorVertices"] == 1
    @test collapsed["CollapsedTriangles"] == 0 && collapsed["RemappedTriangles"] == 0
    @test collapsed["ShortestCollapsedEdge"] ≈ norm(extra[:, 4] .- extra[:, 1])
    @test collapsed["Collapses"][1]["Position"] == extra[:, 4]
    @test collapsed["Collapses"][1]["TargetPosition"] == extra[:, 1]
    @test collapsed["Collapses"][1]["Surface"] == false
    @test collapse.cells_removed == [cells_before - 1]   # the flat cell (first appended) vanished
    @test length(tetrahedra) == cells_before - 1
    @test !(v in Iterators.flatten(tetrahedra))   # the collapsed vertex is orphaned
    @test rule["LayerCells"] == 1 && rule["CellsAboveBoundAfter"] == 0
    @test rule["MaximumEdgeAspectBeforeCollapse"] > 200.0
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
    # The remapped cell is the other cell with v replaced by w.
    @test tetrahedra[end] == ntuple(i -> other[i] == v ? w : other[i], 4)
    # Without the rule the same collapse happens under the scaled-Jacobian gate
    # guard (the other cell is below the gate before and after) and no rule record
    # is written.
    points2, tetrahedra2, _ = corner_fan(0.0125)
    points2 = hcat(points2, extra); tetrahedra2 = vcat(tetrahedra2, [flat, other])
    record2, _, collapse2 = optimize_required_region!(
        points2, tetrahedra2, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9)
    @test record2["EdgeLayerQualityRule"] === nothing
    @test record2["SubSizeEdgeCollapse"]["ScaledJacobianGate"] == 0.01
    @test record2["SubSizeEdgeCollapse"]["CollapsedVertices"] == 1
    @test collapse2.cells_removed == [cells_before - 1] && length(tetrahedra2) == cells_before - 1
    @test record2["SubSizeEdgeCollapse"]["Collapses"][1]["CavityMaximumJacobianCondition"] <
          record2["SubSizeEdgeCollapse"]["Collapses"][1]["ReplacedMaximumJacobianCondition"]
    @test record2["RequiredCellsAboveConditionAfter"] == 0
    @test record2["RequiredMaximumJacobianConditionAfter"] <= 1000.0
    @test record2["MaximumJacobianCondition"] == 1000.0
end

# The metal top face z = 0.1 next to the ridge y = 0 carrying an EdgeSize 1 nm
# layer: ridge nodes 12 nm apart, the first row at 1 nm, the second at 3 nm, and a
# Gmsh face vertex 0.045 nm from the first-row node (the located EL1c cells:
# 0.02-0.24 nm edges between a face vertex and a 1 nm row node, condition up to
# 5116; supervisor decisions 33/34). Vacuum above and metal below the face.
function layer_face_with_sliver()
    points, tetrahedra, triangles = corner_fan(0.0125)
    n = size(points, 2)
    face = hcat([5.0, 0.001, 0.1], [5.00004, 0.00102, 0.1], [5.012, 0.0, 0.1], [4.988, 0.0, 0.1],
                [5.0, 0.003, 0.1], [5.012, 0.003, 0.1], [4.988, 0.003, 0.1],
                [5.0, 0.0015, 0.1025], [5.0, 0.0015, 0.099])
    points = hcat(points, face)
    w, v, a, b, c, d, e, q, p = n .+ (1:9)
    face_triangles = [(b, a, w), (a, v, w), (a, d, v), (d, c, v), (c, w, v), (c, e, w), (e, b, w)]
    cells = vcat([positively_oriented(points, (t..., q)) for t in face_triangles],
                 [positively_oriented(points, (t..., p)) for t in face_triangles])
    entities = vcat(fill(1, length(triangles)), fill(7, length(face_triangles)))
    return (points, vcat(tetrahedra, cells), vcat(triangles, face_triangles), entities,
            (w=w, v=v, a=a, b=b, c=c, d=d, e=e, q=q, p=p))
end

@testset "layer-surface sub-EdgeSize vertices are collapsed along their face (EL1c cells)" begin
    spans = [([1.0, 0.0, 0.1], [9.0, 0.0, 0.1])]
    corners = [(0.0, 0.0, 0.0)]
    edge_size, thickness, zigzag = 0.001, 0.031, 0.05
    points, tetrahedra, triangles, entities, node = layer_face_with_sliver()
    cells_before = length(tetrahedra); triangles_before = length(triangles)
    sliver = [k for (k, cell) in enumerate(tetrahedra) if node.v in cell && node.w in cell]
    @test length(sliver) == 4
    @test maximum(tetrahedron_aspect([points[:, i] for i in tetrahedra[k]]) for k in sliver) > 1000.0
    # The first row is a seeded line (1D support) on the face; the ridge nodes lie
    # on the ridge line and would carry the sidewall too.
    lines = [(node.w, node.c)]; line_entities = [11]
    # Without the collapse the required region fails the condition gate (the
    # slivers are scaled-Jacobian-gated required cells no bounded move repairs).
    record, moved, collapse = optimize_required_region!(
        points, tetrahedra, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 1000.0, 0.75, 1e-9;
        triangle_entities=entities, lines=lines, line_entities=line_entities)
    collapsed = record["SubSizeEdgeCollapse"]
    @test collapsed["CollapsedVertices"] == 1 && collapsed["CollapsedSurfaceVertices"] == 1
    row = collapsed["Collapses"][1]
    @test row["Surface"] == true
    @test row["Position"] ≈ [5.00004, 0.00102, 0.1] && row["TargetPosition"] ≈ [5.0, 0.001, 0.1]
    @test row["EdgeLength"] ≈ norm([0.00004, 0.00002, 0.0])
    @test row["ReplacedMaximumJacobianCondition"] > 1000.0 && row["CavityMaximumJacobianCondition"] < 1000.0
    @test row["CavityMinimumScaledJacobian"] >= min(row["ReplacedMinimumScaledJacobian"], 0.01)
    @test row["CavityMaximumEdgeAspect"] < row["ReplacedMaximumEdgeAspect"]
    # Four cells (both sides of the two triangles on the v-w edge) vanished, the
    # other four cells of v were remapped onto w; two face triangles vanished and
    # two were remapped, keeping their normals; the row line is untouched.
    @test length(collapse.cells_removed) == 4 && length(collapse.cells_remapped) == 4
    @test length(tetrahedra) == cells_before - 4
    @test length(collapse.triangles_removed) == 2 && length(collapse.triangles_remapped) == 2
    @test length(triangles) == triangles_before - 2
    @test isempty(collapse.lines_removed) && isempty(collapse.lines_remapped) && lines == [(node.w, node.c)]
    @test !(node.v in Iterators.flatten(tetrahedra)) && !(node.v in Iterators.flatten(triangles))
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
    @test all(abs(points[3, i] - 0.1) < 1e-15 for t in triangles[(end - 4):end] for i in t)
    for t in triangles[(end - 4):end]
        normal = cross(points[:, t[2]] .- points[:, t[1]], points[:, t[3]] .- points[:, t[1]])
        @test normal[3] > 0.0
    end
    @test record["RequiredCellsAboveConditionAfter"] == 0
    @test record["RequiredMaximumJacobianConditionAfter"] <= 1000.0
    @test record["RequiredCellsBelowGateAfter"] == 0
    # The row node w is never merged onto the face vertex (its line support is not
    # v's) nor onto another face vertex: with v excluded from the candidates, w
    # keeps its sub-size edge and nothing is collapsed.
    points3, tetrahedra3, triangles3, entities3, node3 = layer_face_with_sliver()
    candidate = trues(size(points3, 2)); candidate[node3.v] = false
    supports = vertex_supports(size(points3, 2), triangles3, entities3, lines, line_entities)
    result = collapse_short_edges!(points3, tetrahedra3, triangles3, [(node3.w, node3.c)],
                                   candidate, falses(size(points3, 2)), supports, 0.5 * edge_size;
                                   scaled_jacobian_gate=0.01)
    @test isempty(result.records) && length(tetrahedra3) == cells_before
    # A fixed (CAD point) vertex is never collapsed either (with the row line kept,
    # the row node has no admissible target and nothing is collapsed).
    points4, tetrahedra4, triangles4, entities4, node4 = layer_face_with_sliver()
    fixed = falses(size(points4, 2)); fixed[node4.v] = true
    result4 = collapse_short_edges!(points4, tetrahedra4, triangles4, [(node4.w, node4.c)],
                                    trues(size(points4, 2)), fixed,
                                    vertex_supports(size(points4, 2), triangles4, entities4,
                                                    [(node4.w, node4.c)], line_entities),
                                    0.5 * edge_size; scaled_jacobian_gate=0.01)
    @test isempty(result4.records)
    # Without the row line the row node itself may move along the face onto a
    # face neighbour (its face support is preserved), never onto the fixed vertex.
    points6, tetrahedra6, triangles6, entities6, node6 = layer_face_with_sliver()
    result6 = collapse_short_edges!(points6, tetrahedra6, triangles6, NTuple{2, Int}[],
                                    trues(size(points6, 2)), fixed,
                                    vertex_supports(size(points6, 2), triangles6, entities6,
                                                    NTuple{2, Int}[], Int[]),
                                    0.5 * edge_size; scaled_jacobian_gate=0.01)
    @test length(result6.records) == 1 && result6.records[1]["Vertex"] == node6.w
    @test result6.records[1]["Target"] != node6.v
    @test abs(result6.records[1]["TargetPosition"][3] - 0.1) < 1e-15
    # Without a layer or a corner grading nothing is a candidate and the record
    # says so.
    points0, tetrahedra0, triangles0 = corner_fan(0.0125)
    record0, _, collapse0 = optimize_required_region!(
        points0, tetrahedra0, triangles0, corners, 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, zigzag, 4.0, 0.01, 1000.0, 0.75,
        1e-9)
    @test record0["SubSizeEdgeCollapse"]["CollapseSize"] === nothing
    @test record0["SubSizeEdgeCollapse"]["CandidateVertices"] == 0
    @test record0["SubSizeEdgeCollapse"]["CollapsedVertices"] == 0
    # The condition gate fails closed (a bound below the fan's own condition).
    points5, tetrahedra5, triangles5 = corner_fan(0.0125)
    @test_throws ErrorException optimize_required_region!(
        points5, tetrahedra5, triangles5, corners, 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, zigzag, 4.0, 0.01, 1.5, 0.75,
        1e-9)
end

# The CornerSize 1 nm probe (decision 33 amendment): at the semantic corner
# (2, 0, 0) where the etch sidewall x = 2 (z < 0), the substrate top z = 0 and
# the metal sidewall y = 0 meet, a Gmsh surface vertex on the sidewall x = 2 at
# (2, 0.57, -0.34) nm sits 0.41 nm from the ridge node (2, 0.79, 0) nm; the cells
# on that edge made the corner aspect 4.81 > 4. The vertex is merged along its
# plane onto the ridge node; the corner itself (a CAD point in three planes) and
# the ridge nodes are never moved.
@testset "corner-ball sub-CornerSize surface vertices are collapsed along their plane" begin
    corner = (2.0, 0.0, 0.0)
    xyz = hcat([2.0, 0.0, 0.0], [2.0, 0.00079, 0.0], [2.0, 0.00057, -0.00034], [2.0, 0.0, -0.00075],
               [2.0, 0.0008, -0.0008], [2.00078, 0.0, 0.0], [2.0008, 0.0008, 0.0],
               [2.0008, 0.0, -0.0008], [2.0004, 0.0004, -0.0004])
    o, w, v, b, e, a, f, g, q = 1:9
    # Surface triangles per plane with their entity: x = 2 (1), z = 0 (2), y = 0 (3).
    faces = [((o, w, v), 1), ((o, v, b), 1), ((w, e, v), 1), ((v, e, b), 1),
             ((o, a, w), 2), ((a, f, w), 2), ((o, a, b), 3), ((a, g, b), 3)]
    triangles = [t for (t, _) in faces]; entities = [entity for (_, entity) in faces]
    tetrahedra = [positively_oriented(xyz, (t..., q)) for t in triangles]
    points = copy(xyz)
    fixed = falses(9); fixed[o] = true
    grading = CornerGrading(0.001, 2.0, 0.025, 0.1)
    corner_cells(cells) = [k for (k, cell) in enumerate(cells) if o in cell]
    aspect_before = maximum(tetrahedron_aspect([points[:, i] for i in tetrahedra[k]])
                            for k in corner_cells(tetrahedra))
    @test norm(xyz[:, v] .- xyz[:, w]) < 0.001
    record, moved, collapse = optimize_required_region!(
        points, tetrahedra, triangles, [corner], 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, 0.05, 4.0, 0.01, 1000.0, 0.75,
        1e-9; corner_grading=grading, triangle_entities=entities, fixed=fixed)
    collapsed = record["SubSizeEdgeCollapse"]
    @test collapsed["CollapseSize"] == 0.001 && collapsed["CollapseThreshold"] == 0.0005
    @test collapsed["CandidateVertices"] == 9
    @test collapsed["CollapsedVertices"] >= 1 && collapsed["CollapsedSurfaceVertices"] >= 1
    rows = collapsed["Collapses"]
    first = rows[1]
    @test first["Surface"] == true && first["Position"] ≈ xyz[:, v]
    @test first["TargetPosition"] ≈ xyz[:, w]
    # Every collapsed surface vertex stayed on its own plane: the target lies in
    # every plane of the vertex (here x = 2), and the corner never moved.
    @test all(abs(row["TargetPosition"][1] - 2.0) < 1e-12 for row in rows if row["Surface"])
    @test points[:, o] == xyz[:, o] && !(o in moved)
    @test all(points[:, i] == xyz[:, i] for i in (w, a, b))
    @test !(v in Iterators.flatten(tetrahedra)) && !(v in Iterators.flatten(triangles))
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
    @test record["CornerAspectsAfter"][1] <= 4.0
    @test record["CornerAspectsAfter"][1] <= aspect_before
    @test record["RequiredCellsAboveConditionAfter"] == 0
    # The same vertex without its plane's supports on the target (the ridge node
    # labeled as another face) is not collapsed onto it: supports must include.
    points2 = copy(xyz); tetrahedra2 = [positively_oriented(xyz, (t..., q)) for t in triangles]
    triangles2 = copy(triangles)
    supports = vertex_supports(9, triangles2, entities, NTuple{2, Int}[], Int[])
    delete!(supports[w], (2, 1))
    candidate = falses(9); candidate[v] = true
    result = collapse_short_edges!(points2, tetrahedra2, triangles2, NTuple{2, Int}[], candidate,
                                   fixed, supports, 0.0005; scaled_jacobian_gate=0.01)
    @test all(row["TargetPosition"] != xyz[:, w] for row in result.records)
end

# The Gmsh connectivity mutation of the seed collapse (apply_seed_cell_collapse!)
# round-trips through the seed's writer: a free interior vertex of a tiny box mesh
# is collapsed onto a neighbour, the vanished cells are removed, the remapped cells
# replaced, and the written MSH 2.2 file carries exactly the used points and the
# census element count (no orphaned vertex, no lost or duplicated cell).
@testset "Gmsh seed collapse round-trip: written points == used points, cells == census" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("collapse-round-trip")
        box = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [box], 1, "volume")
        gmsh.model.addPhysicalGroup(2, [tag for (_, tag) in gmsh.model.getEntities(2)], 2, "wall")
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.2)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.2)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        points = reshape(copy(coordinates), 3, :)
        index = Dict(tag => i for (i, tag) in enumerate(node_tags))
        tetrahedra, tags, entities = gmsh_volume_cells(index)
        triangles, _, _ = gmsh_entity_cells(index, 2, Val(3))
        surface = falses(size(points, 2))
        for triangle in triangles, i in triangle
            surface[i] = true
        end
        incident = [Int[] for _ in axes(points, 2)]
        for (k, cell) in enumerate(tetrahedra), i in cell
            push!(incident[i], k)
        end
        interior = [i for i in axes(points, 2) if !surface[i] && !isempty(incident[i])]
        @test !isempty(interior)
        # Collapse the interior vertex v onto its nearest interior neighbour w (a
        # free-free edge, as collapse_short_layer_edges! produces): cells with both
        # vertices vanish, the other cells of v are remapped onto w.
        v = interior[1]
        neighbours = unique(j for k in incident[v] for j in tetrahedra[k] if j != v && !surface[j])
        @test !isempty(neighbours)
        w = neighbours[argmin([norm(points[:, j] .- points[:, v]) for j in neighbours])]
        removed = sort!([k for k in incident[v] if w in tetrahedra[k]])
        remapped = Dict(k => ntuple(i -> tetrahedra[k][i] == v ? w : tetrahedra[k][i], 4)
                        for k in incident[v] if !(w in tetrahedra[k]))
        @test !isempty(removed) && !isempty(remapped)
        expected_cells = length(tetrahedra) - length(removed)
        for (k, cell) in remapped
            tetrahedra[k] = cell
        end
        deleteat!(tetrahedra, removed)
        @test length(tetrahedra) == expected_cells
        @test !(v in Iterators.flatten(tetrahedra))
        used = Set(i for cell in tetrahedra for i in cell)
        expected_keys = Set(sort!([round.(points[:, i]; digits=9) for i in cell]) for cell in tetrahedra)
        model_count = apply_seed_cell_collapse!(node_tags, tags, entities, removed, remapped)
        @test model_count == expected_cells
        path = joinpath(mktempdir(), "collapse-round-trip.msh")
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(path)
        gmsh.model.remove()
        gmsh.open(path)
        written_tags, written_coordinates, _ = gmsh.model.mesh.getNodes()
        written = reshape(written_coordinates, 3, :)
        written_index = Dict(tag => i for (i, tag) in enumerate(written_tags))
        _, written_element_tags, written_nodes = gmsh.model.mesh.getElements(3)
        written_cells = Vector{NTuple{4, Int}}()
        for (block_tags, block) in zip(written_element_tags, written_nodes)
            for start in 1:4:length(block)
                push!(written_cells, ntuple(i -> written_index[block[start + i - 1]], 4))
            end
        end
        @test length(written_cells) == expected_cells
        # Every written point is used by a volume cell (the orphaned v is dropped)
        # and the written point count is the used point count.
        written_used = Set(i for cell in written_cells for i in cell)
        @test length(written_used) == size(written, 2) == length(used)
        @test all(norm(written[:, i] .- points[:, v]) > 1e-9 for i in axes(written, 2))
        written_keys = Set(sort!([round.(written[:, i]; digits=9) for i in cell])
                           for cell in written_cells)
        @test written_keys == expected_keys
    finally
        gmsh.finalize()
    end
end

# The surface collapse round-trips through the Gmsh model too: a face-interior
# vertex of the tiny box is merged along its face (collapse_short_edges! with the
# face entities as supports, a collapse size above every edge so the vertex is a
# candidate), the vanished/remapped tetrahedra, triangles and lines are applied per
# dimension, and the written mesh keeps the wall area, the census counts and only
# the used points.
@testset "Gmsh surface collapse round-trip: triangles remapped along the face, area kept" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("surface-collapse-round-trip")
        box = gmsh.model.occ.addBox(0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(3, [box], 1, "volume")
        gmsh.model.addPhysicalGroup(2, [tag for (_, tag) in gmsh.model.getEntities(2)], 2, "wall")
        gmsh.option.setNumber("Mesh.MeshSizeMin", 0.25)
        gmsh.option.setNumber("Mesh.MeshSizeMax", 0.25)
        gmsh.option.setNumber("Mesh.Algorithm3D", 1)
        gmsh.model.mesh.generate(3)
        node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
        points = reshape(copy(coordinates), 3, :)
        index = Dict(tag => i for (i, tag) in enumerate(node_tags))
        tetrahedra, tags, entities = gmsh_volume_cells(index)
        triangles, triangle_tags, triangle_entities = gmsh_entity_cells(index, 2, Val(3))
        lines, line_tags, line_entities = gmsh_entity_cells(index, 1, Val(2))
        fixed = gmsh_point_nodes(index, size(points, 2))
        @test count(fixed) == 8
        supports = vertex_supports(size(points, 2), triangles, triangle_entities, lines,
                                   line_entities)
        # A face-interior vertex: one face support, no line, not a CAD point.
        interior = [i for i in axes(points, 2)
                    if length(supports[i]) == 1 && first(supports[i])[1] == 2 && !fixed[i]]
        @test !isempty(interior)
        candidate = falses(size(points, 2)); candidate[interior[1]] = true
        v = interior[1]
        face = first(supports[v])[2]
        # Make the vertex a sub-size insertion: slide it along its face halfway to
        # its nearest face neighbour (the model node follows).
        neighbours = unique(j for t in triangles for j in t if v in t && j != v && (2, face) in supports[j])
        w = neighbours[argmin([norm(points[:, j] .- points[:, v]) for j in neighbours])]
        points[:, v] = points[:, w] .+ 0.5 .* (points[:, v] .- points[:, w])
        gmsh.model.mesh.setNode(node_tags[v], points[:, v], Float64[])
        @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
        cells_before = length(tetrahedra); triangles_before = length(triangles)
        result = collapse_short_edges!(points, tetrahedra, triangles, lines, candidate, fixed,
                                       supports, 0.2; scaled_jacobian_gate=0.01)
        @test length(result.records) == 1
        row = result.records[1]
        @test row["Surface"] && row["Vertex"] == v
        @test (2, face) in supports[row["Target"]]
        @test length(tetrahedra) == cells_before - length(result.cells_removed)
        @test length(triangles) == triangles_before - length(result.triangles_removed)
        @test length(result.triangles_removed) == 2 && !isempty(result.triangles_remapped)
        @test !(v in Iterators.flatten(tetrahedra)) && !(v in Iterators.flatten(triangles))
        expected_cells = length(tetrahedra); expected_triangles = length(triangles)
        @test apply_seed_cell_collapse!(node_tags, tags, entities, result.cells_removed,
                                        result.cells_remapped) == expected_cells
        @test apply_seed_cell_collapse!(node_tags, triangle_tags, triangle_entities,
                                        result.triangles_removed, result.triangles_remapped;
                                        dimension=2) == expected_triangles
        @test apply_seed_cell_collapse!(node_tags, line_tags, line_entities, result.lines_removed,
                                        result.lines_remapped; dimension=1) == length(lines)
        path = joinpath(mktempdir(), "surface-collapse-round-trip.msh")
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(path)
        gmsh.model.remove()
        gmsh.open(path)
        written_tags, written_coordinates, _ = gmsh.model.mesh.getNodes()
        written = reshape(written_coordinates, 3, :)
        written_index = Dict(tag => i for (i, tag) in enumerate(written_tags))
        _, cell_tags, cell_nodes = gmsh.model.mesh.getElements(3)
        _, face_tags, face_nodes = gmsh.model.mesh.getElements(2)
        @test sum(length(block) for block in cell_tags) == expected_cells
        @test sum(length(block) for block in face_tags) == expected_triangles
        area = 0.0
        for block in face_nodes, start in 1:3:length(block)
            a, b, c = (written[:, written_index[block[start + i]]] for i in 0:2)
            area += 0.5 * norm(cross(b .- a, c .- a))
        end
        @test area ≈ 6.0 atol=1e-12
        used = Set(written_index[node] for block in cell_nodes for node in block)
        @test length(used) == size(written, 2)
        @test all(norm(written[:, i] .- points[:, v]) > 1e-9 for i in axes(written, 2))
    finally
        gmsh.finalize()
    end
end

# Corner grading (supervisor decision 33): the shell law from CornerSize, its Gmsh
# field expression, the ridge size law through a graded ball, the optimizer's
# corner-graded local bounds and the un-layered edge length census.
@testset "corner grading shells, field expression and local bounds" begin
    grading = CornerGrading(0.004, 2.0, 0.025, 0.1)
    @test corner_shell_radii(grading) ≈ [0.004, 0.012, 0.028, 0.1]
    @test corner_grading_reach(grading) ≈ 0.021
    @test [corner_ball_size(grading, d) for d in (0.0, 0.0039, 0.004, 0.0119, 0.012, 0.0279, 0.028, 0.09)] ≈
          [0.004, 0.004, 0.008, 0.008, 0.016, 0.016, 0.025, 0.025]
    plain = CornerGrading(0.0, 2.0, 0.025, 0.1)
    @test corner_shell_radii(plain) == [0.1] && corner_grading_reach(plain) == 0.0
    @test all(corner_ball_size(plain, d) == 0.025 for d in (0.0, 0.05, 0.1))
    # The field expression is the staircase inside the ball (step(x) = 1 for x >= 0)
    # and the band grading slope beyond the radius; without grading it is lc_fine.
    expression = corner_size_expression("F2", grading, 0.16, 0.2)
    @test expression == "min(0.16,min(0.025,0.004+0.004*step(F2-0.004)+0.008*step(F2-0.012)+" *
                        "0.009000000000000001*step(F2-0.028))+(0.16-0.025)*max(F2-0.1,0)/0.2)"
    @test corner_size_expression("F2", plain, 0.16, 0.2) == "min(0.16,0.025+(0.16-0.025)*max(F2-0.1,0)/0.2)"
    # Ridge size through the ball: the shells inside, lc_fine at the radius, the
    # band slope to lc_tangent outside.
    corners = [(0.0, 0.0, 0.0)]
    @test corner_curve_size([0.002, 0.0, 0.0], corners, grading, 0.05, 0.675) ≈ 0.004
    @test corner_curve_size([0.02, 0.0, 0.0], corners, grading, 0.05, 0.675) ≈ 0.016
    @test corner_curve_size([0.1, 0.0, 0.0], corners, grading, 0.05, 0.675) ≈ 0.025
    @test corner_curve_size([0.137037, 0.0, 0.0], corners, grading, 0.05, 0.675) ≈ 0.05 atol=1e-6
    @test corner_curve_size([0.05, 0.0, 0.0], corners, plain, 0.05, 0.675) ≈ 0.025
    # Optimizer bounds follow the shells: the needle apex of the fan (12.5 nm from
    # the corner, third shell of size 16 nm) may move 0.75 x 16 nm instead of
    # 0.75 x lc_fine; the census records the CornerSize.
    points, tetrahedra, triangles = corner_fan(0.0002)
    record, moved, _ = optimize_required_region!(
        points, tetrahedra, triangles, corners, 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, 0.05, 4.0, 0.01, 1000.0, 0.75, 1e-9;
        corner_grading=grading)
    @test record["CornerSize"] == 0.004
    @test record["MaximumDisplacement"] > 0.0
    @test record["MaximumDisplacement"] <= 0.75 * 0.016 * (1 + 1e-9)
    @test record["MaximumDisplacementOverBound"] <= 1.0 + 1e-9
    graded_displacement = record["MaximumDisplacement"]
    points2, tetrahedra2, triangles2 = corner_fan(0.0002)
    record2, _, _ = optimize_required_region!(
        points2, tetrahedra2, triangles2, corners, 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, 0.05, 4.0, 0.01, 1000.0, 0.75, 1e-9)
    @test record2["CornerSize"] == 0.0
    @test record2["MaximumDisplacement"] >= graded_displacement
    # Un-layered edge length per corner: the distance from a corner to the nearest
    # span end of each layered curve ending at it (curves not ending there ignored).
    curves = [Dict{String, Any}("Start" => [0.1, 0.0, 0.0], "End" => [9.9, 0.0, 0.0],
                                "CurveEnds" => [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]),
              Dict{String, Any}("Start" => [0.0, 0.12, 0.0], "End" => [0.0, 7.9, 0.0],
                                "CurveEnds" => [[0.0, 0.0, 0.0], [0.0, 8.0, 0.0]])]
    rows = unlayered_edge_length_per_corner([(0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (5.0, 5.0, 0.0)],
                                            curves, 0.1)
    @test rows[1]["LayeredEdges"] == 2 && rows[1]["UnlayeredLengths"] ≈ [0.1, 0.12]
    @test rows[1]["Maximum"] ≈ 0.12 && rows[1]["Minimum"] ≈ 0.1
    @test rows[2]["LayeredEdges"] == 1 && rows[2]["Maximum"] ≈ 0.1
    @test rows[3]["LayeredEdges"] == 0 && rows[3]["Maximum"] === nothing
end
