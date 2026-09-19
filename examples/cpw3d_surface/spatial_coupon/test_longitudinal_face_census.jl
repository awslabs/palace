# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# A ridge-to-ridge face (sidewall) whose longitudinal curves run along x from a
# semantic corner at the origin. Sizes are the coupon recipe's; the face height is
# arbitrary and the census must not depend on it.
const LC_FINE = 0.025
const LC_TANGENT = 0.1
const LC_FAR = 0.16
const SLOPE = (LC_FAR - LC_FINE) / 0.2
const RADIUS = 0.1
const CORNER = (0.0, 0.0, 0.0)

function sidewall_face(length, height)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.model.add("sidewall")
    surface = gmsh.model.occ.addRectangle(0.0, 0.0, 0.0, length, height)
    gmsh.model.occ.synchronize()
    longitudinal = Int32[]
    for (_, curve) in gmsh.model.getBoundary([(2, surface)], false, false, false)
        lower, upper = gmsh.model.getParametrizationBounds(1, curve)
        derivative = gmsh.model.getDerivative(1, curve, [0.5 * (lower[1] + upper[1])])
        abs(derivative[1]) / norm(derivative[1:3]) > 1.0 - 1.0e-9 && push!(longitudinal, curve)
    end
    field = gmsh.model.mesh.field
    field.add("AttractorAnisoCurve", 1)
    field.setNumbers(1, "CurvesList", Float64.(longitudinal))
    field.setNumber(1, "DistMin", 0.0); field.setNumber(1, "DistMax", 0.2)
    field.setNumber(1, "SizeMinNormal", LC_FINE); field.setNumber(1, "SizeMaxNormal", LC_FAR)
    field.setNumber(1, "SizeMinTangent", LC_TANGENT); field.setNumber(1, "SizeMaxTangent", LC_FAR)
    field.setNumber(1, "Sampling", 100)
    field.add("MathEval", 2)
    field.setString(2, "F", "min(F1,$(LC_FAR))")
    field.setAsBackgroundMesh(2)
    for (name, value) in [("Mesh.MeshSizeMin", LC_FINE), ("Mesh.MeshSizeMax", LC_FAR),
                          ("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0),
                          ("Mesh.MeshSizeFromCurvature", 0), ("Mesh.MinimumCurvePoints", 3)]
        gmsh.option.setNumber(name, value)
    end
    return sort!(longitudinal)
end

function census_after_generation(longitudinal)
    gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    gmsh.model.mesh.generate(2)
    reach = corner_law_reach(RADIUS, LC_FINE, LC_TANGENT, SLOPE)
    rows = longitudinal_face_census(longitudinal, [CORNER], reach)
    gmsh.finalize()
    return rows
end

@testset "Misaligned ridge rows leave full-height triangles the census detects" begin
    longitudinal = sidewall_face(8.0, 0.1)
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    # Uniform rows with different counts: the ridge nodes drift out of alignment.
    for (curve, intervals) in zip(longitudinal, (80, 83))
        lower, upper = gmsh.model.getParametrizationBounds(1, curve)
        parameters = collect(range(lower[1], upper[1]; length=intervals + 1))[2:(end - 1)]
        add_explicit_curve_mesh!(curve, parameters, gmsh.model.getValue(1, curve, parameters),
                                 next_node, point_nodes)
    end
    rows = census_after_generation(longitudinal)
    @test length(rows) == 1
    row = rows[1]
    @test row["LongitudinalCurves"] == 2
    @test row["FullHeightTriangles"] > 0
    @test row["FullHeightTrianglesAwayFromCorners"] > 0
    @test row["InteriorNodesAwayFromCorners"] < row["Triangles"] ÷ 8
    @test length(row["InteriorNodeHistogramAlongEdge"]) == LONGITUDINAL_FACE_HISTOGRAM_BINS
    @test sum(row["InteriorNodeHistogramAlongEdge"]) == row["InteriorNodesAwayFromCorners"]
end

@testset "Grid-preserving corner law keeps the interior row and no full-height triangles" begin
    longitudinal = sidewall_face(8.0, 0.1)
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    reach = corner_law_reach(RADIUS, LC_FINE, LC_TANGENT, SLOPE)
    for curve in longitudinal
        placed = corner_isotropic_curve_nodes(curve, [CORNER], CornerGrading(0.0, 2.0, LC_FINE, RADIUS),
                                              LC_TANGENT, SLOPE)
        @test placed !== nothing
        parameters, coordinates = placed
        xyz = reshape(coordinates, 3, :)
        # Beyond the first grid node past the law's reach (within one grid interval of
        # it) every node sits on the transfinite lc_tangent grid.
        lower, upper = gmsh.model.getParametrizationBounds(1, curve)
        intervals = ceil(Int, gmsh.model.occ.getMass(1, curve) / LC_TANGENT)
        grid = collect(range(lower[1], upper[1]; length=intervals + 1))
        for (parameter, i) in zip(parameters, axes(xyz, 2))
            norm(xyz[:, i] .- CORNER) > reach + LC_TANGENT || continue
            @test minimum(abs.(grid .- parameter)) <= 1.0e-12 * (upper[1] - lower[1])
        end
        # Inside the ball the spacing is the isotropic size, never above it (the gap
        # is equidistributed with ceil(integral) intervals, decision 43).
        inside = sort!([xyz[1, i] for i in axes(xyz, 2) if norm(xyz[:, i] .- CORNER) <= RADIUS])
        if length(inside) >= 2
            @test all(0.8 * LC_FINE .<= diff(inside) .<= LC_FINE * (1.0 + 1.0e-9))
        end
        add_explicit_curve_mesh!(curve, parameters, coordinates, next_node, point_nodes)
    end
    rows = census_after_generation(longitudinal)
    @test length(rows) == 1
    row = rows[1]
    @test row["FullHeightTriangles"] == 0
    @test row["FullHeightTrianglesAwayFromCorners"] == 0
    @test row["InteriorNodesAwayFromCorners"] > row["Triangles"] ÷ 8
    @test all(>(0), row["InteriorNodeHistogramAlongEdge"][2:end])
end

@testset "Composed curve law grades a dip strictly inside a grid interval (decision 43)" begin
    longitudinal = sidewall_face(8.0, 0.1)
    curve = longitudinal[1]
    lower, upper = gmsh.model.getParametrizationBounds(1, curve)
    intervals = ceil(Int, gmsh.model.occ.getMass(1, curve) / LC_TANGENT)
    grid = collect(range(lower[1], upper[1]; length=intervals + 1))
    grid_x = [gmsh.model.getValue(1, curve, [g])[1] for g in grid]
    # A trace-like dip: 20 nm at x = 3.05 growing at slope 4 back to the spacing
    # within 20 nm, entirely inside the grid interval [3.0, 3.1]; the limiter
    # (slope 0.5 at growth 2) widens the graded region to 3.05 +/- 0.16 = within
    # the neighbouring intervals [2.8, 3.3], whose grid nodes keep the law < spacing.
    dip = 3.05
    law(point) = min(LC_TANGENT, 0.02 + 4.0 * abs(point[1] - dip))
    grading = CornerGrading(0.0, 2.0, LC_FINE, RADIUS)
    far_corner = [(100.0, 0.0, 0.0)]
    @test composed_curve_nodes(curve, LC_TANGENT, p -> LC_TANGENT, 2.0, far_corner, grading) === nothing
    placed = composed_curve_nodes(curve, LC_TANGENT, law, 2.0, far_corner, grading)
    @test placed !== nothing
    parameters, coordinates, record = placed
    xyz = reshape(coordinates, 3, :)
    x = sort!(collect(xyz[1, :]))
    @test issorted(parameters)
    # Every grid node outside the limiter's reach stays a node; the grid nodes at
    # 2.9 / 3.0 / 3.1 / 3.2 lie under the limited law and are replaced (five grid
    # intervals graded).
    for g in grid_x[2:(end - 1)]
        if abs(g - dip) > 0.16 + 1.0e-9
            @test minimum(abs.(x .- g)) <= 1.0e-9
        end
    end
    nodes = vcat(grid_x[1], x, grid_x[end])
    spacings = diff(nodes)
    @test record["GridIntervalsKept"] == intervals - 5
    @test count(3.0 < xi < 3.1 for xi in x) >= 3
    @test all(s -> s <= LC_TANGENT * (1.0 + 1.0e-9), spacings)
    @test minimum(spacings) <= 0.03
    # Each interval is no longer than the (unlimited) law at its midpoint, and the
    # neighbour ratio stays within the growth cap along the whole curve.
    midpoints = 0.5 .* (nodes[1:(end - 1)] .+ nodes[2:end])
    @test all(spacings[i] <= law((midpoints[i], 0.0, 0.0)) * (1.0 + 1.0e-9) for i in eachindex(spacings))
    ratios = spacings[2:end] ./ spacings[1:(end - 1)]
    @test max(maximum(ratios), 1.0 / minimum(ratios)) <= 2.0 * (1.0 + 1.0e-9)
    @test record["Graded"] && record["Spacing"] == LC_TANGENT && record["GridIntervals"] == intervals
    @test record["NodeSpacing"]["Minimum"] == minimum(spacings)
    @test record["NodeSpacing"]["Maximum"] <= LC_TANGENT * (1.0 + 1.0e-9)
    @test record["AchievedOverPrescribed"]["Maximum"] <= 1.0 + 1.0e-8
    @test record["PrescribedMinimum"] <= 0.03
    gmsh.finalize()
end

@testset "Curves out of the law's reach keep the transfinite spacing" begin
    sidewall_face(8.0, 0.1)
    far_corner = (100.0, 0.0, 0.0)
    for (_, curve) in gmsh.model.getEntities(1)
        @test corner_isotropic_curve_nodes(curve, [far_corner],
                                           CornerGrading(0.0, 2.0, LC_FINE, RADIUS), LC_TANGENT,
                                           SLOPE) === nothing
    end
    gmsh.finalize()
end
