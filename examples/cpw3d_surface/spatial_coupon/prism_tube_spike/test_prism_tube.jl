# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Targeted tests of the prism-tube SPIKE: cross-section topology and orientation,
# tube frames, and the tiny hybrid build (prisms + explicit pyramids + tets) as a
# smoke test with its element counts and label areas.
#
# usage: julia --project=test/examples test_prism_tube.jl

using Test
include(joinpath(@__DIR__, "mesh_tiny_hybrid.jl"))

@testset "Cross-section rings and rays" begin
    section = TubeSection(0.00025, 2.0, 7, [-90.0 + 30.0 * j for j in 0:9], fill(2, 9))
    @test ring_sizes(section) ≈ [0.00025 * 2.0^k for k in 0:6]
    @test tube_radius(section) ≈ 0.03175
    @test section_node_count(section) == 1 + 7 * 10
    @test material_groups(section) == [(0, 8, 2)]
    @test cad_rays(section) == [0, 9]
    bottom = TubeSection(0.00025, 2.0, 7, [180.0 + 30.0 * j for j in 0:9], vcat(fill(1, 3), fill(2, 6)))
    @test material_groups(bottom) == [(0, 2, 1), (3, 8, 2)]
    @test cad_rays(bottom) == [0, 3, 9]
    # Every sector triangle is counterclockwise in (u, w) and the triangles tile
    # the sector: 1 fan + 2 per annulus.
    uw = section_coordinates(section)
    for j in 0:8
        triangles = sector_triangles(section, j)
        @test length(triangles) == 1 + 2 * 6
        for (a, b, c) in triangles
            area = 0.5 * ((uw[1, b] - uw[1, a]) * (uw[2, c] - uw[2, a]) -
                          (uw[2, b] - uw[2, a]) * (uw[1, c] - uw[1, a]))
            @test area > 0.0
        end
    end
    sector_area = sum(0.5 * ((uw[1, b] - uw[1, a]) * (uw[2, c] - uw[2, a]) -
                             (uw[2, b] - uw[2, a]) * (uw[1, c] - uw[1, a]))
                      for j in 0:8 for (a, b, c) in sector_triangles(section, j))
    @test sector_area ≈ 9 * 0.5 * tube_radius(section)^2 * sin(deg2rad(30.0))
    @test_throws ErrorException TubeSection(0.0, 2.0, 7, [0.0, 90.0], [2])
    @test_throws ErrorException TubeSection(0.001, 2.0, 7, [0.0, 90.0], [2, 2])
end

@testset "Tube frame" begin
    tube = EdgeTube([2.0, 0.0, 0.1], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 0.1, 7.9, 0.05)
    @test tube.layers == 156
    @test tube.e ≈ [1.0, 0.0, 0.0]
    @test tube_point(tube, 0.0, 0.0, tube.s_start) ≈ [2.1, 0.0, 0.1]
    @test tube_station(tube, tube.layers) ≈ 7.9
    # A length that is not a multiple of the spacing takes one more layer; the
    # achieved spacing (<= the requested one) is recorded by tube_spacing.
    uneven = EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 1.03, 0.05)
    @test uneven.layers == 21 && tube_spacing(uneven) ≈ 1.03 / 21
    @test_throws ErrorException EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 1.03, 0.0)
end

@testset "Tiny hybrid build" begin
    directory = mktempdir()
    census = tiny_hybrid_mesh("hybrid", joinpath(directory, "tiny.msh");
                              census_path=joinpath(directory, "census.json"))
    @test census["Prism"]["Count"] == 2 * 40 * 9 * 13
    @test census["Pyramid"]["Count"] == 2 * 40 * 9
    @test census["Prism"]["NonPositive"] == 0
    @test census["Pyramid"]["NonPositive"] == 0
    @test census["Tetrahedron"]["NonPositive"] == 0
    @test census["DuplicateVolumeElementsRemoved"] == 0
    @test census["Areas"]["6001"]["Area"] ≈ 2.2
    @test census["Areas"]["5001"]["Area"] ≈ 2.0
    @test census["Areas"]["3100"]["Area"] ≈ 1.2
    @test census["Areas"]["6001"]["Elements"]["Quadrangle"] == 3 * 7 * 40
    @test census["Areas"]["5001"]["Elements"]["Quadrangle"] == 7 * 40
    @test census["Areas"]["3100"]["Elements"]["Quadrangle"] == 7 * 40
    @test isfile(joinpath(directory, "tiny.msh"))
    @test isfile(joinpath(directory, "census.json"))
end
