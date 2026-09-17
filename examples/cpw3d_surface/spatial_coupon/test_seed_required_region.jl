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
