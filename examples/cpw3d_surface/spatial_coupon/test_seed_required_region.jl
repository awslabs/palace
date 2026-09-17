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
                                              4.0, 0.01, 0.75, 1e-9)
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
                                           0.0, zigzag, 4.0, 0.01, 0.75, 1e-9)
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
        thickness, zigzag, 4.0, 0.01, 0.75, 1e-9)
    # With the rule the slab is gated by orientation and edge aspect: the aspect
    # repair lifts the apex within the bound and the scaled Jacobian is a diagnostic.
    moved_points = copy(points); moved_cells = copy(tetrahedra)
    record, moved, removed, remapped = optimize_required_region!(
        moved_points, moved_cells, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 0.75, 1e-9; edge_layer_maximum_aspect=200.0)
    rule = record["EdgeLayerQualityRule"]
    @test rule["LayerCells"] == 1 && rule["CollapsedVertices"] == 0
    @test isempty(removed) && isempty(remapped) && moved_cells == tetrahedra
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
        thickness, zigzag, 4.0, 0.01, 0.75, 1e-9; edge_layer_maximum_aspect=1.5)
    # The rule needs a seeded layer.
    points0, tetrahedra0, triangles0 = corner_fan(0.0125)
    @test_throws ErrorException optimize_required_region!(
        points0, tetrahedra0, triangles0, corners, 0.1, 0.025,
        Tuple{Vector{Float64}, Vector{Float64}}[], 0.0, 2.0, 0.0, zigzag, 4.0, 0.01, 0.75, 1e-9;
        edge_layer_maximum_aspect=100.0)
    # Edge aspect: unit right tetrahedron sqrt(6); a 50 x 50 x 1 slab ~ 70.7.
    @test tetrahedron_edge_aspect([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]) ≈ sqrt(6)
    @test tetrahedron_edge_aspect([[0.0, 0, 0], [50.0, 0, 0], [0, 50.0, 0], [0, 0, 1.0]]) ≈ 70.7389567 atol=1e-6
    @test tetrahedron_edge_aspect([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [1.0, 1.0, 0]]) == Inf
end

# A seed volume vertex 0.28 nm from a row node (below EdgeSize = the adapter hmin)
# inside the layer: collapsed onto that node; the cell containing both vanishes,
# the other incident cell is remapped, and the survivors satisfy the rule.
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
    # The slab base is a face (row nodes are surface vertices, never collapsed).
    triangles = vcat(triangles, [(w, p2, p3)])
    cells_before = length(tetrahedra)
    @test tetrahedron_edge_aspect([points[:, i] for i in flat]) > 200.0
    record, moved, removed, remapped = optimize_required_region!(
        points, tetrahedra, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 0.75, 1e-9; edge_layer_maximum_aspect=200.0)
    rule = record["EdgeLayerQualityRule"]
    @test rule["CollapsedVertices"] == 1 && rule["CollapsedCells"] == 1 && rule["RemappedCells"] == 1
    @test rule["ShortestCollapsedEdge"] ≈ norm(extra[:, 4] .- extra[:, 1])
    @test removed == [cells_before - 1]           # the flat cell (first appended) vanished
    @test length(tetrahedra) == cells_before - 1
    @test !(v in Iterators.flatten(tetrahedra))   # the collapsed vertex is orphaned
    @test rule["LayerCells"] == 1 && rule["CellsAboveBoundAfter"] == 0
    @test all(tetrahedron_scaled_jacobian([points[:, i] for i in cell]) > 0 for cell in tetrahedra)
    # The remapped cell is the other cell with v replaced by w.
    @test tetrahedra[end] == ntuple(i -> other[i] == v ? w : other[i], 4)
    # Without the rule nothing is collapsed (the scaled-Jacobian repair moves the
    # free vertex instead) and no rule record is written.
    points2, tetrahedra2, _ = corner_fan(0.0125)
    points2 = hcat(points2, extra); tetrahedra2 = vcat(tetrahedra2, [flat, other])
    record2, _, removed2, remapped2 = optimize_required_region!(
        points2, tetrahedra2, triangles, corners, 0.1, 0.025, spans, edge_size, 2.0,
        thickness, zigzag, 4.0, 0.01, 0.75, 1e-9)
    @test record2["EdgeLayerQualityRule"] === nothing
    @test isempty(removed2) && isempty(remapped2) && length(tetrahedra2) == cells_before
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
        triangles = gmsh_linear_cells(index, 2, Val(3))
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
