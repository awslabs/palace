# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

const TOL = FOOTPRINT_COLLINEAR_TOLERANCE

rotate(points, angle, shift) = [
    (cos(angle) * p[1] - sin(angle) * p[2] + shift[1],
     sin(angle) * p[1] + cos(angle) * p[2] + shift[2]) for p in points]

@testset "Near-collinear duplicate vertices are merged within the tolerance" begin
    # A rectangle whose top edge carries the four-edge device footprint's
    # near-collinear chain: turns of 2e-8, 8e-7 and 1.6e-6 radians (CSV precision)
    # plus an exact duplicate vertex. One CAD face per genuine facet remains.
    chain = [(8.812013816937261, 5.923428529140843), (8.512016587819119, 5.924717596130073),
             (8.212019383521778, 5.926006417356859), (7.912022126551145, 5.927295712420648)]
    points = [(10.0, 0.0), (10.0, 5.922139468088279), (9.11201104699479, 5.922139468088279),
              chain..., chain[end], (7.0, 5.93), (7.0, 0.0)]
    simplified, record = simplify_footprint_polygon(points, TOL)
    @test record["OriginalVertices"] == length(points)
    @test record["RemovedVertexCount"] == length(points) - record["Vertices"]
    @test record["RemovedVertexIndices"] ⊆ 3:8
    # Exactly one of the two coincident vertices (7, 8) is removed.
    @test count(in(record["RemovedVertexIndices"]), (7, 8)) == 1
    @test record["MaximumRelativeDeviation"] <= TOL
    @test record["MaximumDeviation"] <= TOL * record["MaximumDeviationLocalScale"]
    # Every removed vertex lies within the tolerance of the simplified polygon.
    for index in record["RemovedVertexIndices"]
        distance = minimum(point_segment_distance_2d(points[index], simplified[i],
                                                     simplified[mod1(i + 1, length(simplified))])
                           for i in eachindex(simplified))
        @test distance <= TOL * maximum(hypot((simplified[i] .- simplified[mod1(i + 1, length(simplified))])...)
                                        for i in eachindex(simplified))
    end
    # Remaining consecutive edges are never within the tolerance of each other:
    # the metric stage classifies every kept wall junction as a genuine dihedral.
    for i in eachindex(simplified)
        a = simplified[mod1(i - 1, length(simplified))]; b = simplified[i]
        c = simplified[mod1(i + 1, length(simplified))]
        u = (b[1] - a[1], b[2] - a[2]); v = (c[1] - b[1], c[2] - b[2])
        @test abs(cross2d(u, v)) / (hypot(u...) * hypot(v...)) > TOL
    end
    # The whole chain lies within the tolerance of one edge from vertex 3 to the
    # chain's end (the individual 1.6e-6 turn does not make a facet on its own).
    @test record["RemovedVertexIndices"] == [4, 5, 6, 7]
    @test record["Vertices"] == length(simplified)
    @test simplified[1] == points[1] && simplified[end] == points[end]
end

@testset "A genuine bend is kept and exact polygons are untouched" begin
    square = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
    simplified, record = simplify_footprint_polygon(square, TOL)
    @test simplified == square
    @test record["RemovedVertexCount"] == 0 && isempty(record["RemovedVertexIndices"])
    @test record["MaximumDeviation"] == 0.0
    # A 5-degree bend in the middle of the top edge is a facet, not noise.
    bend = 0.5 * tand(5.0)
    bent = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.0 + bend), (0.0, 1.0)]
    simplified, record = simplify_footprint_polygon(bent, TOL)
    @test simplified == bent
    @test record["RemovedVertexCount"] == 0
    # A bend just above the tolerance stays; one just below it merges.
    for (angle, kept) in ((4.0 * TOL, true), (0.5 * TOL, false))
        offset = 0.5 * tan(angle)
        polygon = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.5, 1.0 + offset), (0.0, 1.0)]
        simplified, record = simplify_footprint_polygon(polygon, TOL)
        @test (record["RemovedVertexCount"] == 0) == kept
        kept || @test record["RemovedVertexIndices"] == [4]
    end
    # A triangle is never reduced; fewer than three vertices is an error.
    triangle = [(0.0, 0.0), (1.0, 0.0), (0.5, 1e-9)]
    @test simplify_footprint_polygon(triangle, TOL)[1] == triangle
    @test_throws ErrorException simplify_footprint_polygon([(0.0, 0.0), (1.0, 0.0)], TOL)
    @test_throws ErrorException simplify_footprint_polygon(square, 0.0)
end

@testset "Simplification is covariant under in-plane rigid motion" begin
    chain = [(8.812013816937261, 5.923428529140843), (8.512016587819119, 5.924717596130073),
             (8.212019383521778, 5.926006417356859), (7.912022126551145, 5.927295712420648)]
    points = [(10.0, 0.0), (10.0, 5.922139468088279), (9.11201104699479, 5.922139468088279),
              chain..., (7.0, 5.93), (7.0, 0.0)]
    _, record = simplify_footprint_polygon(points, TOL)
    @test record["RemovedVertexCount"] > 0
    for (angle, shift) in ((0.63, (1.2, -0.7)), (-2.1, (0.0, 3.0)), (pi / 2, (-5.0, 5.0)))
        moved = rotate(points, angle, shift)
        simplified, moved_record = simplify_footprint_polygon(moved, TOL)
        @test moved_record["RemovedVertexIndices"] == record["RemovedVertexIndices"]
        @test moved_record["Vertices"] == record["Vertices"]
        @test isapprox(moved_record["MaximumDeviation"], record["MaximumDeviation"];
                       rtol=1e-6, atol=1e-15)
        expected = rotate([points[i] for i in eachindex(points)
                           if !(i in record["RemovedVertexIndices"])], angle, shift)
        @test all(isapprox(collect(a), collect(b); atol=1e-12)
                  for (a, b) in zip(simplified, expected))
    end
end
