# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# The Gmsh-only prism tube recipe of mesh_spatial_coupon.jl (supervisor decision
# 38): the mixed-element quality census computed by the mesher itself (Gmsh's
# minSJ is the high-order mapping Jacobian, identically 1 for order-1 elements),
# its fail-closed gates, the tube ring rule, the recorded extrusion spacing, the
# straight metal edge segments with their corner clearance and the explicit
# volume size laws.

# A discrete volume carrying hand-built linear cells of every type.
function install_cells!(points, cells)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("mixed")
    volume = gmsh.model.addDiscreteEntity(3)
    gmsh.model.mesh.addNodes(3, volume, collect(1:size(points, 2)), vec(points))
    for (type, connectivity) in cells
        gmsh.model.mesh.addElementsByType(volume, type, Int[], vec(hcat(connectivity...)))
    end
    return volume
end

@testset "mixed volume quality: corner scaled Jacobian, condition, orientation" begin
    # A regular tetrahedron, a sliver (SJ ~ 6e-4), a right prism, a right pyramid.
    points = hcat([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
                  [2.0, 0.0, 0.0], [3.0, 0.0, 0.0], [2.0, 1.0, 0.0], [2.5, 0.5, 0.0006],
                  [4.0, 0.0, 0.0], [5.0, 0.0, 0.0], [4.0, 1.0, 0.0],
                  [4.0, 0.0, 1.0], [5.0, 0.0, 1.0], [4.0, 1.0, 1.0],
                  [6.0, 0.0, 0.0], [7.0, 0.0, 0.0], [7.0, 1.0, 0.0], [6.0, 1.0, 0.0], [6.5, 0.5, 1.0])
    cells = [(4, [[1, 2, 3, 4], [5, 6, 7, 8]]), (6, [[9, 10, 11, 12, 13, 14]]),
             (7, [[15, 16, 17, 18, 19]])]
    install_cells!(points, cells)
    quality = mixed_volume_quality()
    gmsh.finalize()
    @test quality["Total"] == 4
    tets = quality["Tetrahedron"]
    @test tets["Count"] == 2 && tets["PositiveOrientation"]
    # The sliver: |det| / (|e1| |e2| |e3|) of the vertex-0 frame (production definition).
    e1 = [1.0, 0.0, 0.0]; e2 = [0.0, 1.0, 0.0]; e3 = [0.5, 0.5, 0.0006]
    sliver = abs(det(hcat(e1, e2, e3))) / (norm(e1) * norm(e2) * norm(e3))
    @test tets["MinimumScaledJacobian"] ≈ sliver
    @test sliver < 0.01 && tets["CellsBelowScaledJacobian0.01"] == 1
    @test tets["MaximumJacobianCondition"] > 1000.0 && tets["CellsAboveCondition1000"] == 1
    prisms = quality["Prism"]
    @test prisms["Count"] == 1 && prisms["PositiveOrientation"]
    # The right-angle corners have an orthonormal frame; at the two acute corners of
    # the right-triangle section the frame is ((-1,1,0),(-1,0,0),(0,0,1)) up to symmetry:
    # scaled Jacobian 1/sqrt(2) and condition phi^2.
    @test prisms["MinimumScaledJacobian"] ≈ 1.0 / sqrt(2.0)
    @test prisms["MaximumJacobianCondition"] ≈ ((1.0 + sqrt(5.0)) / 2.0)^2
    @test prisms["ScaledJacobianQuantiles"][end] ≈ 1.0 / sqrt(2.0)
    pyramids = quality["Pyramid"]
    @test pyramids["Count"] == 1 && pyramids["PositiveOrientation"]
    @test 0.0 < pyramids["MinimumScaledJacobian"] < 1.0
    @test isfinite(pyramids["MaximumJacobianCondition"]) && pyramids["MaximumJacobianCondition"] > 1.0
    # Every corner frame of the linear cells appears once.
    @test length(VOLUME_CORNER_FRAMES[4]) == 1 && length(VOLUME_CORNER_FRAMES[6]) == 6 &&
          length(VOLUME_CORNER_FRAMES[7]) == 4
    @test_throws ErrorException gate_mixed_volume_quality(quality, 0.01, 1000.0)
    quality["Tetrahedron"]["MinimumScaledJacobian"] = 0.02
    quality["Tetrahedron"]["MaximumJacobianCondition"] = 10.0
    @test gate_mixed_volume_quality(quality, 0.01, 1000.0) === nothing
    quality["Prism"]["MaximumJacobianCondition"] = 2000.0
    @test_throws ErrorException gate_mixed_volume_quality(quality, 0.01, 1000.0)
end

@testset "mixed volume quality: an inverted prism fails orientation" begin
    points = hcat([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
                  [0.0, 0.0, -1.0], [1.0, 0.0, -1.0], [0.0, 1.0, -1.0])
    install_cells!(points, [(6, [[1, 2, 3, 4, 5, 6]])])
    quality = mixed_volume_quality()
    gmsh.finalize()
    @test !quality["Prism"]["PositiveOrientation"] && quality["Prism"]["NonpositiveCells"] == 1
    @test_throws ErrorException gate_mixed_volume_quality(quality, 0.01, 1000.0)
end

@testset "tube ring rule and extrusion spacing" begin
    # r_K + h_K <= bound: 0.25 nm rings of ratio 2 up to K = 7 (31.75 + 16 nm) fit 50 nm,
    # K = 8 (63.75 + 32 nm) does not.
    @test tube_ring_count(0.00025, 2.0, 0.05) == 7
    @test tube_ring_count(0.00025, 2.0, 0.1) == 8
    @test tube_ring_count(0.001, 2.0, 0.0025) == 1
    @test_throws ErrorException tube_ring_count(0.002, 2.0, 0.0025)
    tube = EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 1.0, 0.05)
    @test tube.layers == 20 && tube_spacing(tube) ≈ 0.05
    @test tube.stations ≈ collect(0.0:0.05:1.0) && tube_station(tube, 20) == 1.0
    @test tube_station(tube, 2.5) ≈ 0.125
    tube = EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 1.01, 0.05)
    @test tube.layers == 21 && tube_spacing(tube) ≈ 1.01 / 21 && tube_spacing(tube) <= 0.05
    @test_throws ErrorException EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0],
                                         0.0, 1.0, 0.0)
end

@testset "tube layers follow the size field along the axis (decision 40)" begin
    tube = EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 1.0, 0.05)
    # A uniform field at the spacing gives uniform layers within the spacing.
    stations, positions, sizes = graded_tube_stations(0.0, 0.99, s -> 0.05, 0.05, 2.0)
    @test length(stations) == 21 && stations ≈ collect(range(0.0, 0.99; length=21))
    @test positions[1] == 0.0 && positions[end] == 0.99 && all(sizes .<= 0.05)
    # A field above the spacing is capped at the spacing; the layers never exceed it.
    stations, _, _ = graded_tube_stations(0.0, 1.0, s -> 1.0, 0.05, 2.0)
    @test length(stations) in (21, 22) && maximum(diff(stations)) < 0.05
    @test maximum(diff(stations)) ≈ 0.05 rtol = 0.05
    # A fine end: the corner-ball staircase 0.25 -> 25 nm within 0.1 from s = 0,
    # then the spacing. The first layer is the inner size, the layers grow by at
    # most the growth ratio, the far layers are the spacing, and the tube length
    # is preserved exactly.
    grading = CornerGrading(0.00025, 2.0, 0.025, 0.1)
    fine_end(s) = min(0.05, corner_ball_size(grading, s) + 0.675 * max(s - 0.1, 0.0))
    stations, positions, sizes = graded_tube_stations(0.0, 1.0, fine_end, 0.05, 2.0)
    thickness = diff(stations)
    graded = EdgeTube(tube, stations)
    @test graded.layers == length(thickness) && graded.stations[end] == 1.0
    # (a layer equidistributed in 1 / size is the size at its midpoint, 1.3 x the
    # inner size where the field grows at the limited slope)
    @test 0.5 * 0.00025 <= thickness[1] <= 2.0 * 0.00025
    @test all(thickness .> 0.0) && all(thickness .<= 0.05)
    ratios = thickness[2:end] ./ thickness[1:(end - 1)]
    @test maximum(ratios) <= 2.0 && minimum(ratios) >= 0.5
    @test maximum(ratios) > 1.5                      # the geometric grading is used
    @test thickness[end] ≈ 0.05 rtol = 0.05          # the far end is at the spacing
    @test count(>(0.045), thickness) >= 15           # most of the tube is at the spacing
    @test graded.layers > 20 + 7                     # the fine end adds the graded layers
    # The layer thickness follows the prescribed size along the whole axis.
    statistics = tube_layer_statistics(graded, positions, sizes)
    @test statistics["Minimum"] ≈ thickness[1]
    @test statistics["Maximum"] ≈ 0.05 rtol = 0.05
    @test statistics["AtStart"] ≈ thickness[1] && statistics["AtEnd"] ≈ thickness[end]
    @test statistics["PrescribedAtStart"] ≈ 0.00025 && statistics["PrescribedAtEnd"] ≈ 0.05
    @test 0.7 <= statistics["AchievedOverPrescribed"]["Minimum"]
    @test statistics["AchievedOverPrescribed"]["Maximum"] <= 1.3
    @test statistics["MaximumNeighbourRatio"] <= 2.0
    # The mesh follows the stations: the pyramid apex of layer i sits mid-layer.
    @test tube_station(graded, 0) == 0.0 && tube_station(graded, graded.layers) == 1.0
    @test tube_station(graded, 0.5) ≈ 0.5 * thickness[1]
    # Negative: uniform layers where the field prescribes finer. The uniform first
    # layer spans sizes down to the inner size (200 x too thick); every graded layer
    # is within 1.5 x the finest field over its span (the limited slope 1/2 bounds
    # the field variation across a layer to a quarter of the layer).
    @test 0.05 / minimum(fine_end.(range(0.0, 0.05; length=201))) >= 200.0
    @test all(thickness[i] <= 1.5 * minimum(fine_end.(range(stations[i], stations[i + 1]; length=9)))
              for i in eachindex(thickness))
    # A field that steps down sharply is gradient-limited so that the neighbour
    # ratio stays within the growth; the limiter reach is (growth - 1) / growth.
    step_field(s) = s < 0.5 ? 0.05 : 0.001
    stations, positions, sizes = graded_tube_stations(0.0, 1.0, step_field, 0.05, 2.0)
    thickness = diff(stations)
    @test maximum(thickness[2:end] ./ thickness[1:(end - 1)]) <= 2.0
    @test minimum(thickness[2:end] ./ thickness[1:(end - 1)]) >= 0.5
    @test minimum(thickness) ≈ 0.001 rtol = 0.15
    # A surface end (decision 41): the trace rule prescribes 0.0217 on the cut
    # surface at s = 0 growing at FarGrowth 0.5 into the volume. Equidistribution
    # puts the first layer at its midpoint size (1.3 x the surface value); with
    # `surface_start` the first layer is the surface value, the second layer within
    # the growth, the rest equidistributed and the length preserved.
    surface_field(s) = min(0.05, 0.0217 + 0.5 * s)
    midpoint, _, _ = graded_tube_stations(0.0, 1.0, surface_field, 0.05, 2.0)
    @test 1.2 * 0.0217 <= diff(midpoint)[1] <= 1.4 * 0.0217
    surfaced, positions, sizes = graded_tube_stations(0.0, 1.0, surface_field, 0.05, 2.0;
                                                      surface_start=true)
    thickness = diff(surfaced)
    @test thickness[1] ≈ 0.0217 atol = 1.0e-12
    @test 1.0 < thickness[2] / thickness[1] <= 2.0
    @test all(thickness .> 0.0) && all(thickness .< 0.05) && surfaced[end] == 1.0
    @test length(surfaced) in (length(midpoint), length(midpoint) + 1)
    statistics = tube_layer_statistics(EdgeTube(tube, surfaced), positions, sizes)
    @test statistics["AtStart"] <= statistics["PrescribedAtStart"] ≈ 0.0217
    @test statistics["AchievedOverPrescribed"]["Minimum"] < 0.9      # the surface layer
    @test statistics["MaximumNeighbourRatio"] <= 2.0
    # A surface end where the field dips below the surface value within the first
    # layer: the layer is at most the field over its span (0.04 at the surface, 0.02
    # at s = 0.04: the iteration 0.04 -> 0.02 stops at 0.02, the field over
    # [0, 0.02] being >= 0.03); both ends on a surface; the end layer at the far end.
    dip_field(s) = min(0.05, 0.02 + 0.5 * abs(s - 0.04))
    dipped, positions, sizes = graded_tube_stations(0.0, 1.0, dip_field, 0.05, 2.0;
                                                    surface_start=true)
    @test diff(dipped)[1] ≈ 0.02 rtol = 0.05                 # (sampled field minimum)
    @test diff(dipped)[1] <= minimum(dip_field.(range(0.0, diff(dipped)[1]; length=201)))
    both, _, _ = graded_tube_stations(0.0, 1.0, s -> min(0.05, 0.0217 + 0.5 * min(s, 1.0 - s)),
                                      0.05, 2.0; surface_start=true, surface_end=true)
    @test diff(both)[1] ≈ 0.0217 atol = 1.0e-12
    @test diff(both)[end] ≈ 0.0217 atol = 1.0e-12
    @test maximum(diff(both)[2:end] ./ diff(both)[1:(end - 1)]) <= 2.0
    @test_throws ErrorException graded_tube_stations(0.0, 0.03, surface_field, 0.05, 2.0;
                                                     surface_start=true, surface_end=true)
    # Without the flags the stations are unchanged by the surface option.
    @test graded_tube_stations(0.0, 1.0, surface_field, 0.05, 2.0)[1] == midpoint
    # Guards: a non-positive field, a non-increasing station list, wrong ends.
    @test_throws ErrorException graded_tube_stations(0.0, 1.0, s -> 0.0, 0.05, 2.0)
    @test_throws ErrorException graded_tube_stations(0.0, 1.0, s -> 0.05, 0.05, 1.0)
    @test_throws ErrorException EdgeTube(tube, [0.0, 0.5, 0.5, 1.0])
    @test_throws ErrorException EdgeTube(tube, [0.0, 0.5, 0.9])
    # The axis size law composes the corner law, the trace rule and the band laws
    # without the tube rule (which reads NormalSize on the axis).
    record = prepare_tube_band_sizing!([([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])],
                                       [([1.0, -1.0, 0.0], [1.0, 1.0, 0.0])],
                                       0.04, 0.025, 0.16, 0.5, [(0.0, 0.0, 0.0)], 0.1)
    @test record["BandSegments"] == 1
    @test tube_band_size(0.5, 0.0, 0.0, 1.0) ≈ 0.025            # the tube rule on the axis
    @test feature_band_size(0.5, 0.0, 0.0, 1.0) ≈ 0.16          # band + corner exterior only
    @test tube_axis_size((0.5, 0.0, 0.0), 0.05, [(0.0, 0.0, 0.0)], grading, 0.675) ≈ 0.05
    @test tube_axis_size((1.0, 0.0, 0.0), 0.05, [(0.0, 0.0, 0.0)], grading, 0.675) ≈ 0.025
    @test tube_axis_size((0.98, 0.0, 0.0), 0.05, [(0.0, 0.0, 0.0)], grading, 0.675) ≈ 0.045
    @test tube_axis_size((0.0, 0.0, 0.0), 0.05, [(0.0, 0.0, 0.0)], grading, 0.675) ≈ 0.00025
    @test tube_axis_size((0.05, 0.0, 0.0), 0.05, [(0.0, 0.0, 0.0)], grading, 0.675) ≈ 0.025
    @test tube_axis_size((0.12, 0.0, 0.0), 0.05, [(0.0, 0.0, 0.0)], grading, 0.675) ≈
          min(0.05, 0.025 + 0.5 * 0.02)                          # corner exterior law
    empty!(TUBE_AXIS_SEGMENTS); empty!(BAND_SEGMENTS); empty!(CORNER_EXTERIOR_POINTS)
end

@testset "metal edge segments: outward normals, corner clearance, box continuation" begin
    # An L-shaped conductor whose two arms end on the outer box x = 10 and y = 8;
    # the semantic corners are the loop vertices inside the box.
    points = [(10.0, 0.0), (2.0, 0.0), (0.0, -2.0), (0.0, 8.0), (-1.0, 8.0), (-1.0, -3.0),
              (10.0, -3.0)]
    loop = (conductor=1, plane=0.0, hole=false, points=points,
            classes=fill("Physical", length(points)))
    corners = [(2.0, 0.0, 0.0), (0.0, -2.0, 0.0), (-1.0, -3.0, 0.0)]
    lower = [-1.0, -3.0]; upper = [10.0, 8.0]
    clearance(angle) = 0.03 / tan(0.5 * angle) + 0.016
    segments = metal_edge_segments([loop], corners, clearance, lower, upper, 1.0e-9)
    # Box sides are skipped: (0,8)->(-1,8) on y = 8, (-1,8)->(-1,-3) on x = -1,
    # (-1,-3)->(10,-3) on y = -3, (10,-3)->(10,0) on x = 10.
    @test length(segments) == 3
    first = segments[1]
    @test first.start == [10.0, 0.0] && first.stop == [2.0, 0.0]
    @test first.normal ≈ [0.0, 1.0]      # away from the metal below the edge
    @test first.s_start == 0.0            # (10, 0) is a box continuation vertex
    angle = acos(dot([1.0, 0.0], [-1.0, -1.0] ./ sqrt(2.0)))
    @test first.s_end ≈ 8.0 - clearance(angle)
    @test first.corner_angles[1] ≈ pi && first.corner_angles[2] ≈ angle
    second = segments[2]
    @test second.normal ≈ [-1.0, 1.0] ./ sqrt(2.0)  # the metal lies below the diagonal
    @test second.s_start ≈ clearance(angle)
    third = segments[3]
    @test third.start == [0.0, -2.0] && third.stop == [0.0, 8.0] && third.normal ≈ [1.0, 0.0]
    @test third.s_end ≈ 10.0             # (0, 8) lies on the box
    # A loop vertex that is neither a corner nor on the box is rejected.
    @test_throws ErrorException metal_edge_segments([loop], corners[2:end], clearance, lower,
                                                    upper, 1.0e-9)
end

@testset "interior conductor loops (holes): inward tube normals, tube counts, facing width (decision 48)" begin
    # An outer square conductor with a square hole; every vertex is a semantic corner.
    outer = [(-0.8, -0.8), (0.8, -0.8), (0.8, 0.8), (-0.8, 0.8)]
    inner = [(-0.3, -0.3), (0.3, -0.3), (0.3, 0.3), (-0.3, 0.3)]
    loops = [(conductor=1, plane=0.0, hole=false, points=outer, classes=fill("Physical", 4)),
             (conductor=1, plane=0.0, hole=true, points=inner, classes=fill("Physical", 4))]
    corners = [(p[1], p[2], 0.0) for p in vcat(outer, inner)]
    lower = [-2.3, -2.3]; upper = [2.3, 2.3]
    clearance(angle) = 0.03 / tan(0.5 * angle) + 0.016
    segments = metal_edge_segments(loops, corners, clearance, lower, upper, 1.0e-9)
    # Two tubes per straight side of every loop: 8 sides -> 16 tubes, the same count
    # rule as an exterior loop (metal_loop_records agrees).
    @test length(segments) == 8
    @test 2 * length(segments) ==
          2 * sum(record["Sides"] for record in metal_loop_records(loops, lower, upper, 1.0e-9))
    for segment in segments
        midpoint = 0.5 .* (segment.start .+ segment.stop)
        # The normal points away from the metal: outward for the exterior loop, into the
        # hole (towards its centre) for the hole loop.
        towards_centre = dot(segment.normal, -midpoint) > 0.0
        @test segment.hole == towards_centre
        @test abs(norm(segment.normal) - 1.0) <= 1.0e-12
        # Every hole corner is a right angle between two tube edges: the same clearance
        # as an exterior right-angle corner.
        @test all(angle -> isapprox(angle, pi / 2), segment.corner_angles)
        @test segment.s_start == clearance(pi / 2) && segment.s_end == segment.span - clearance(pi / 2)
    end
    # A hole side missing from a device etch footprint fails closed.
    footprint = [(conductor=1, plane=0.0, hole=true, points=[(-0.3, -0.3), (0.3, -0.3), (0.3, 0.2), (-0.3, 0.2)],
                  classes=fill("Physical", 4))]
    hole_sides = [segment for segment in segments if segment.hole]
    @test assert_etch_carries_edge(footprint, hole_sides[1], 1.0e-9)
    message = try assert_etch_carries_edge(footprint, hole_sides[3], 1.0e-9); "" catch e; e.msg end
    @test occursin("ScopeGuard[FootprintWithoutEdge]", message)
    # The facing width of a hole is the distance between its non-adjacent sides.
    @test hole_facing_width(inner, 1.0e-9) ≈ 0.6
    @test hole_facing_width([(0.0, 0.0), (1.0, 0.0), (1.0, 0.05), (0.0, 0.05)], 1.0e-9) ≈ 0.05
    @test segment_segment_distance_2d((0.0, 0.0), (1.0, 0.0), (0.5, -1.0), (0.5, 1.0)) == 0.0
    @test segment_segment_distance_2d((0.0, 0.0), (1.0, 0.0), (2.0, 1.0), (3.0, 1.0)) ≈ sqrt(2.0)
end

# A single-conductor, single-slot coupon whose counterclockwise plan-view loop is
# `points`, written as the frozen inputs (signature, boundary, mask, semantic
# contract) of the production mesher, with the process normal `sign`.
function write_loop_inputs(directory, points, sign; plane=0.0)
    n = length(points)
    open(joinpath(directory, "signature.csv"), "w") do io
        println(io, "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm")
        for i in 1:n
            a = points[i]; b = points[i % n + 1]
            span = hypot(b[1] - a[1], b[2] - a[2])
            t = ((b[1] - a[1]) / span, (b[2] - a[2]) / span)
            println(io, join([i, 0, 1, 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), plane, t[2], -t[1], 0,
                              t[1], t[2], 0, sign, -span / 2, span / 2, 0], ","))
        end
    end
    open(joinpath(directory, "boundary.csv"), "w") do io
        println(io, "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y")
        for (i, point) in enumerate(points)
            println(io, join([1, i, 1, plane, 0, "Physical", point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for point in points
            println(io, join([1, 1, plane, point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "semantic.json"), "w") do io
        write_json(io, Dict{String, Any}(
            "Version" => 1, "SemanticCorners" => [[p[1], p[2], plane] for p in points]))
    end
    return joinpath(directory, "signature.csv"), joinpath(directory, "boundary.csv"),
           joinpath(directory, "mask.csv"), joinpath(directory, "semantic.json")
end

write_strip_inputs(directory, sign; x=0.6, y=0.2, plane=0.0) =
    write_loop_inputs(directory, [(-x, -y), (x, -y), (x, y), (-x, y)], sign; plane=plane)

# Build the loop coupon under coarse prism-tube options at Radius 0.5; returns the
# census and the mesh path.
function build_strip_coupon(directory, sign; points=nothing, stem=sign > 0 ? "up" : "down")
    signature, boundary, mask, semantic = points === nothing ?
        write_strip_inputs(directory, sign) : write_loop_inputs(directory, points, sign)
    mesh = joinpath(directory, "coupon-$stem.msh")
    census = joinpath(directory, "census-$stem.json")
    generate_spatial_coupon(; signature=signature, mask=mask, boundary=boundary, fabricated=true,
                            filename=mesh, radius=0.5, metal_thickness=0.1, overetch=0.05,
                            sidewall_angle=90.0, top_rounding=0.0, trench_rounding=0.0,
                            lc_fine=0.05, lc_tangent=0.1, lc_far=0.3, max_nodes=2_000_000,
                            max_elements=2_000_000, semantic_contract=semantic,
                            corner_isotropy_radius=0.1, corner_census=census, edge_size=0.01,
                            edge_growth_ratio=2.0, corner_size=0.01, prism_tubes=true,
                            far_growth=0.5, maximum_corner_aspect=4.0,
                            minimum_scaled_jacobian=0.01, maximum_jacobian_condition=1000.0,
                            quality_displacement_over_normal=0.75)
    return parse_json(read(census, String)), mesh
end

# Node coordinates and per-type element node tags of a written mesh.
function read_mesh_cells(path)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(path)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = Dict(tag => coordinates[(3i - 2):(3i)] for (i, tag) in enumerate(node_tags))
    types, _, element_nodes = gmsh.model.mesh.getElements(3)
    cells = Dict(Int(type) => reshape(block, GMSH_LINEAR_VOLUME_TYPES[Int(type)][2], :)
                 for (type, block) in zip(types, element_nodes))
    gmsh.finalize()
    return points, cells
end

@testset "downward process layers: tube frame, layer gap, box padding, mirror covariance (decision 48)" begin
    guard_message(f) = try f(); "" catch e; e.msg end
    # The coupon box pads each layer's plane by Overetch on the substrate side (-Nz) and
    # MetalThickness on the metal side (+Nz).
    row(z, sign) = (slot=0, conductor=1, point=(0.0, 0.0, z), gap=(0.0, -1.0, 0.0),
                    tangent=(1.0, 0.0, 0.0), interval=(-1.0, 1.0), normal_sign=Float64(sign),
                    vertex_arm=false)
    up_lower, up_upper = coupon_bounds([row(0.0, 1)], 2.0, 0.1, 0.05)
    @test up_lower[3] ≈ -2.05 && up_upper[3] ≈ 2.1
    down_lower, down_upper = coupon_bounds([row(0.0, -1)], 2.0, 0.1, 0.05)
    @test down_lower[3] ≈ -2.1 && down_upper[3] ≈ 2.05
    opposed_lower, opposed_upper = coupon_bounds([row(0.0, 1), row(0.6, -1)], 2.0, 0.1, 0.05)
    @test opposed_lower[3] ≈ -2.05 && opposed_upper[3] ≈ 2.65
    @test "DownwardLayers" in RECIPE_SCOPE_SUPPORTED_CLASSES &&
          any(guard -> guard[1] == "NarrowLayerGap", RECIPE_SCOPE_GUARDS)
    # The tube frame of a downward layer: b = -z, the top tube on plane - thickness, the
    # extrusion sense mirrored; a narrow gap between facing layers fails closed.
    points = [(-0.6, -0.2), (0.6, -0.2), (0.6, 0.2), (-0.6, 0.2)]
    loop(plane) = (conductor=1, plane=plane, hole=false, points=points, classes=fill("Physical", 4))
    corners(plane) = [(p[1], p[2], plane) for p in points]
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("frames")
    lower = [-2.1, -1.7]; upper = [2.1, 1.7]
    _, _, up_tubes, _, description = build_edge_tubes!(
        gmsh.model.occ, [(plane=0.0, sign=1, edges=nothing)], [loop(0.0)], nothing, corners(0.0),
        0.01, 2.0, 30.0, 0.1, 0.05, 0.1, 0.1, 0.05, lower, upper, 1.0e-9)
    _, _, down_tubes, _, _ = build_edge_tubes!(
        gmsh.model.occ, [(plane=0.0, sign=-1, edges=nothing)], [loop(0.0)], nothing, corners(0.0),
        0.01, 2.0, 30.0, 0.1, 0.05, 0.1, 0.1, 0.05, lower, upper, 1.0e-9)
    @test length(up_tubes) == length(down_tubes) == 8
    for ((up, up_section), (down, down_section)) in zip(up_tubes, down_tubes)
        @test up.b == [0.0, 0.0, 1.0] && down.b == [0.0, 0.0, -1.0]
        @test up.n == down.n && up.e ≈ -down.e
        @test up.origin[3] ≈ -down.origin[3]                # plane + thickness <-> plane - thickness
        @test up_section.materials == down_section.materials && up_section.angles == down_section.angles
        # The tube covers the same edge interval in space: its end points are mirrored.
        ends(t) = sort([tube_point(t, 0.0, 0.0, t.s_start)[1:2], tube_point(t, 0.0, 0.0, t.s_end)[1:2]])
        @test ends(up) ≈ ends(down)
    end
    @test occursin("tube frame", lowercase(description["TubeFrameRule"])) || haskey(description, "TubeFrameRule")
    facing_reach = description["FacingReach"]
    @test facing_reach ≈ description["Radius"] + description["PyramidHeight"] + 2 * 0.05
    # Facing layers: the gap between the metal faces must exceed twice the reach.
    narrow = [(plane=0.0, sign=1, edges=nothing), (plane=0.2 + 2 * facing_reach, sign=-1, edges=nothing)]
    message = guard_message(() -> build_edge_tubes!(
        gmsh.model.occ, narrow, [loop(0.0), loop(narrow[2].plane)], nothing,
        vcat(corners(0.0), corners(narrow[2].plane)), 0.01, 2.0, 30.0, 0.1, 0.05, 0.1, 0.1, 0.05,
        lower, upper, 1.0e-9))
    @test occursin("ScopeGuard[NarrowLayerGap]", message)
    wide = [(plane=0.0, sign=1, edges=nothing), (plane=0.2 + 2 * facing_reach + 1.0e-3, sign=-1, edges=nothing)]
    _, _, wide_tubes, _, _ = build_edge_tubes!(
        gmsh.model.occ, wide, [loop(0.0), loop(wide[2].plane)], nothing,
        vcat(corners(0.0), corners(wide[2].plane)), 0.01, 2.0, 30.0, 0.1, 0.05, 0.1, 0.1, 0.05,
        lower, upper, 1.0e-9)
    @test length(wide_tubes) == 16
    gmsh.finalize()
    # Mirror covariance: a downward-only coupon is the z-reflection of the upward one.
    mktempdir() do directory
        up_census, up_mesh = build_strip_coupon(directory, 1)
        down_census, down_mesh = build_strip_coupon(directory, -1)
        up_rows = up_census["PrismTubes"]["Tubes"]; down_rows = down_census["PrismTubes"]["Tubes"]
        @test all(row["Layer"] == 1 for row in up_rows) && all(row["Layer"] == -1 for row in down_rows)
        @test all(row["Origin"][3] ≈ (row["Edge"] == "top" ? 0.1 : 0.0) for row in up_rows)
        @test all(row["Origin"][3] ≈ (row["Edge"] == "top" ? -0.1 : 0.0) for row in down_rows)
        @test up_census["Scope"]["ExhibitedClasses"] == ["ExteriorLoops"]
        @test down_census["Scope"]["ExhibitedClasses"] == ["DownwardLayers", "ExteriorLoops"]
        @test up_census["CouponBox"]["Lower"][3] ≈ -down_census["CouponBox"]["Upper"][3]
        @test up_census["CouponBox"]["Upper"][3] ≈ -down_census["CouponBox"]["Lower"][3]
        # Everything the recipe determines mirrors exactly: the interface areas (CAD),
        # the tube prisms and pyramids (installed explicitly) and their node sets; every
        # gate passes on both. Gmsh's tetrahedral meshing of the mirrored CAD is not
        # reflection-covariant (the Delaunay kernel depends on orientation), so the tet
        # count and its quality extremes are reported, not asserted equal.
        @test [row["Attribute"] for row in up_census["InterfaceAreas"]] ==
              [row["Attribute"] for row in down_census["InterfaceAreas"]]
        @test [row["Area"] for row in up_census["InterfaceAreas"]] ≈
              [row["Area"] for row in down_census["InterfaceAreas"]]
        up_quality = up_census["PrismTubes"]["Quality"]; down_quality = down_census["PrismTubes"]["Quality"]
        for kind in ("Tetrahedron", "Prism", "Pyramid")
            @test up_quality[kind]["PositiveOrientation"] && down_quality[kind]["PositiveOrientation"]
            @test max(up_quality[kind]["MaximumJacobianCondition"],
                      down_quality[kind]["MaximumJacobianCondition"]) <= 1000.0
        end
        @test min(up_quality["Tetrahedron"]["MinimumScaledJacobian"],
                  down_quality["Tetrahedron"]["MinimumScaledJacobian"]) >= 0.01
        @test all(<=(4.0), up_census["SeedQualityOptimization"]["CornerAspectsAfter"])
        @test all(<=(4.0), down_census["SeedQualityOptimization"]["CornerAspectsAfter"])
        @info "mirror covariance: tetrahedra" up=up_quality["Tetrahedron"]["Count"] down=down_quality["Tetrahedron"]["Count"]
        # The tube rows mirror: same edge intervals and layer counts, origin z negated,
        # extrusion sense reversed (e = n x b with b = (0, 0, Nz)).
        @test length(up_rows) == length(down_rows) == 8
        for (up_row, down_row) in zip(up_rows, down_rows)
            @test up_row["Layers"] == down_row["Layers"] && up_row["Edge"] == down_row["Edge"]
            @test up_row["Length"] ≈ down_row["Length"] && up_row["Start"] ≈ down_row["Start"]
            @test up_row["Extrusion"] ≈ -down_row["Extrusion"] && up_row["Normal"] ≈ down_row["Normal"]
            @test up_row["Origin"][3] ≈ -down_row["Origin"][3]
        end
        up_points, up_cells = read_mesh_cells(up_mesh)
        down_points, down_cells = read_mesh_cells(down_mesh)
        for type in (6, 7)
            @test size(up_cells[type]) == size(down_cells[type])
            @test up_quality[GMSH_LINEAR_VOLUME_TYPES[type][1]]["Count"] == size(up_cells[type], 2)
            # The tube node sets are mirrored up to the axial layer placement: the layer
            # stations follow the axis size field sampled adaptively from s_start
            # (graded_tube_stations, decision 40), so reversing the extrusion sense moves a
            # station by less than one sampling step TangentialSize /
            # TUBE_LAYER_SAMPLES_PER_SIZE (the same direction dependence two oppositely
            # traversed exterior edges have); the cross-section (rings, rays, cap ends) is
            # exact.
            up_nodes = [up_points[tag] for tag in unique(vec(up_cells[type]))]
            down_nodes = [down_points[tag] .* [1.0, 1.0, -1.0] for tag in unique(vec(down_cells[type]))]
            @test length(up_nodes) == length(down_nodes)
            deviation = maximum(minimum(norm(a .- b) for b in down_nodes) for a in up_nodes)
            @info "mirror covariance: tube node deviation" type deviation
            @test deviation <= 0.1 / TUBE_LAYER_SAMPLES_PER_SIZE
            ends(points, cells, flip) = sort([points[tag] .* [1.0, 1.0, flip] for tag in unique(vec(cells[type]))
                                              if abs(abs(points[tag][1]) - 0.55) <= 1.0e-9 ||
                                                 abs(abs(points[tag][2]) - 0.15) <= 1.0e-9];
                                             by=xyz -> round.(xyz, digits=9))
            up_ends = ends(up_points, up_cells, 1.0); down_ends = ends(down_points, down_cells, -1.0)
            @test !isempty(up_ends) && length(up_ends) == length(down_ends)
            @test maximum(norm(a .- b) for (a, b) in zip(up_ends, down_ends)) <= 1.0e-9
        end
    end
end

@testset "un-etched plane of a Radius-0.5 coupon is labeled 3000 (decisions 48 / 49)" begin
    # The producer-default collar (3 x Radius) of an L-shaped conductor leaves the box
    # corner opposite its notch un-etched: box [-1.5, 2.7] x [-1.5, 2.5] minus the
    # mitered collar polygon = the notch [1.8, 2.7] x [1.8, 2.5] = 0.63 um^2 at the layer
    # plane. Before decision 49 the mesher compared the surface's z-range with the plane
    # at 1e-7 x Radius = 5e-8, below the 1e-7 padding of the OCC bounding box, so this
    # plane was labeled 3100 (etched trench) for every coupon with Radius < 1 um.
    l_shape = [(0.0, 0.0), (1.2, 0.0), (1.2, 0.3), (0.3, 0.3), (0.3, 1.0), (0.0, 1.0)]
    mktempdir() do directory
        census, _ = build_strip_coupon(directory, 1; points=l_shape, stem="l-shape")
        areas = Dict(row["Attribute"] => row["Area"] for row in census["InterfaceAreas"])
        @test haskey(areas, 3000) && haskey(areas, 3100)
        @test isapprox(areas[3000], 0.63; atol=1.0e-9)
        # The trench: the collar floor (16.17 minus the metal 0.57) plus the trench walls
        # under the metal edges (perimeter 4.4) and along the notch (0.9 + 0.7; the other
        # collar sides lie on the box), both x Overetch 0.05.
        @test isapprox(areas[3100], 16.17 - 0.57 + (4.4 + 1.6) * 0.05; atol=1.0e-9)
        @test census["CouponBox"]["Lower"][1:2] ≈ [-1.5, -1.5] && census["CouponBox"]["Upper"][1:2] ≈ [2.7, 2.5]
        @test all(polygon["Construction"] == "MiterOffset" for polygon in census["FootprintPolygons"])
        @test census["FootprintSimplification"]["CollarUnionPolygons"] == 0
    end
end

@testset "producer-default collar of facing sides closer than twice the collar: union footprint (decision 54a)" begin
    guard_message(f) = try f(); "" catch e; e.msg end
    tolerance = 1.0e-9
    # A metal band y in [0, 10] with a notch x in (-4, 4), y in (0, 6) open to the gap
    # below; box [-10, 10]^2. Under the 6 um collar the two notch walls face each other
    # across 8 < 12 um, so the miter polygon folds back through the notch (the wall at
    # x = -4 offsets to x = 2, the wall at x = 4 to x = -2) and self-intersects.
    notched(width) = [(-10.0, 0.0), (-width / 2, 0.0), (-width / 2, 6.0), (width / 2, 6.0),
                      (width / 2, 0.0), (10.0, 0.0), (10.0, 10.0), (-10.0, 10.0)]
    classes = vcat(fill("Physical", 5), fill("Continuation", 3))
    loop(width) = (conductor=1, plane=0.0, hole=false, points=notched(width), classes=classes)
    box = ([-10.0, -10.0], [10.0, 10.0])
    narrow = offset_loop_points(loop(8.0), -6.0, tolerance)
    @test polygon_is_simple(notched(8.0), tolerance) && !polygon_is_simple(narrow, tolerance)
    @test !polygon_is_simple([(0.0, 0.0), (1.0, 0.0), (1.0, 0.0), (0.0, 1.0)], tolerance)   # zero side
    @test !polygon_is_simple([(0.0, 0.0), (2.0, 0.0), (1.0, 0.0), (1.0, 1.0)], tolerance)   # fold-back
    points, construction = collar_loop_points(loop(8.0), -6.0, box, tolerance)
    @test construction == "CollarUnion"
    # The union of the loop, the side rectangles and the convex-corner kites clipped to
    # the box is the half-plane y >= -6 of the box: the notch is etched throughout.
    simplified, _ = simplify_footprint_polygon(points, FOOTPRINT_COLLINEAR_TOLERANCE)
    @test simplified == [(-10.0, -6.0), (10.0, -6.0), (10.0, 10.0), (-10.0, 10.0)]
    @test polygon_area2(points) > 0.0
    # A notch wider than twice the collar keeps the miter polygon (the same region).
    wide, wide_construction = collar_loop_points(loop(16.0), -6.0, box, tolerance)
    @test wide_construction == "MiterOffset" && wide == offset_loop_points(loop(16.0), -6.0, tolerance)
    @test polygon_is_simple(wide, tolerance)
    # Zero offset is the loop itself; an inward (positive) self-intersecting offset has
    # no union form; the union needs the box.
    @test collar_loop_points(loop(8.0), 0.0, nothing, tolerance) == (notched(8.0), "MiterOffset")
    @test occursin("Inward offset", guard_message(() -> collar_loop_points(
        (conductor=1, plane=0.0, hole=false, points=notched(8.0), classes=fill("Physical", 8)),
        3.0, box, tolerance)))
    @test occursin("coupon box", guard_message(() -> collar_loop_points(loop(8.0), -6.0, nothing, tolerance)))
    # The pieces are counterclockwise convex polygons inside the box: the loop, five side
    # rectangles (the three box sides are Continuation) and the two convex notch-opening
    # kites (the box junctions and the loop's convex box corners give degenerate kites).
    pieces = collar_pieces(loop(8.0), -6.0, narrow, box..., tolerance)
    @test length(pieces) == 1 + 5 + 2
    @test all(polygon_area2(piece) > 0.0 for piece in pieces)
    @test all(-10.0 - tolerance <= p[d] <= 10.0 + tolerance for piece in pieces for p in piece for d in 1:2)
    @test [(-4.0, 0.0), (-4.0, -6.0), (2.0, -6.0), (2.0, 0.0)] in pieces
    # A keyhole (entry 8 wide, chamber 32 wide and 14 tall under a 6 um collar) leaves an
    # un-etched island inside the chamber: the collar region is not one polygon.
    keyhole = [(-20.0, 0.0), (-4.0, 0.0), (-4.0, 4.0), (-16.0, 4.0), (-16.0, 18.0), (16.0, 18.0),
               (16.0, 4.0), (4.0, 4.0), (4.0, 0.0), (20.0, 0.0), (20.0, 20.0), (-20.0, 20.0)]
    island = (conductor=1, plane=0.0, hole=false, points=keyhole,
              classes=vcat(fill("Physical", 10), fill("Continuation", 2)))
    message = guard_message(() -> collar_loop_points(island, -6.0, ([-20.0, -10.0], [20.0, 20.0]), tolerance))
    @test occursin("ScopeGuard[FootprintTopology]", message) && occursin("island", message)
    @test any(guard -> guard[1] == "FootprintTopology" && guard[2] == "build", RECIPE_SCOPE_GUARDS)
    # The failure mode of the device 5-edge coupon ("tube tool ... has 3 volume
    # descendants") on a Radius-0.5 coupon: a notch exactly as wide as the 1.5 um collar,
    # so the miter side of each notch wall lands on the facing wall, the axis of that
    # wall's tubes, and the OCC face built from the self-intersecting polygon split the
    # bottom tube's vacuum sectors. With the union footprint every tube is matched, the
    # collar covers the whole box (every un-notched side is within 1.5 um of the box) and
    # the census records the construction.
    notch = [(-2.0, -0.2), (2.0, -0.2), (2.0, 0.8), (0.75, 0.8), (0.75, 0.0), (-0.75, 0.0),
             (-0.75, 0.8), (-2.0, 0.8)]
    mktempdir() do directory
        census, _ = build_strip_coupon(directory, 1; points=notch, stem="notch")
        polygons = census["FootprintPolygons"]
        @test length(polygons) == 1 && polygons[1]["Construction"] == "CollarUnion"
        @test census["FootprintSimplification"]["CollarUnionPolygons"] == 1
        lower = census["CouponBox"]["Lower"]; upper = census["CouponBox"]["Upper"]
        @test Set(Tuple.(polygons[1]["Points"])) ==
              Set([(lower[1], lower[2]), (upper[1], lower[2]), (upper[1], upper[2]), (lower[1], upper[2])])
        # 8 sides x (top tube + bottom substrate + bottom vacuum sectors) = 24 tube volumes,
        # every one a single fragment descendant.
        @test length(census["PrismTubes"]["Tubes"]) == 16
        @test length(census["PrismTubes"]["Volumes"]) == 24
        areas = Dict(row["Attribute"] => row["Area"] for row in census["InterfaceAreas"])
        @test !haskey(areas, 3000)
        @test census["PrismTubes"]["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
    end
end

@testset "recipe scope: guard ids, exhibited classes, metal loop records (decision 48)" begin
    guard_message(f) = try f(); "" catch e; e.msg end
    points = [(10.0, 0.0), (2.0, 0.0), (0.0, -2.0), (0.0, 8.0), (-1.0, 8.0), (-1.0, -3.0),
              (10.0, -3.0)]
    loop = (conductor=1, plane=0.0, hole=false, points=points,
            classes=fill("Physical", length(points)))
    corners = [(2.0, 0.0, 0.0), (0.0, -2.0, 0.0), (-1.0, -3.0, 0.0)]
    lower = [-1.0, -3.0]; upper = [10.0, 8.0]
    clearance(angle) = 0.03 / tan(0.5 * angle) + 0.016
    # Every guard id is unique, has a statement and a detection origin.
    ids = [guard[1] for guard in RECIPE_SCOPE_GUARDS]
    @test allunique(ids) && isempty(intersect(ids, RECIPE_SCOPE_SUPPORTED_CLASSES))
    @test all(guard[2] in ("inputs", "build") && !isempty(guard[3]) for guard in RECIPE_SCOPE_GUARDS)
    @test_throws ErrorException scope_error("NotAGuard", "")
    # The guard messages carry ScopeGuard[<id>].
    hole = (conductor=1, plane=0.0, hole=true, points=points, classes=loop.classes)
    @test occursin("ScopeGuard[FreeEdgeEnds]",
                   guard_message(() -> metal_edge_segments([loop], corners[2:end], clearance, lower, upper, 1.0e-9)))
    @test occursin("ScopeGuard[ShortEdges]",
                   guard_message(() -> metal_edge_segments([loop], corners, angle -> 5.0, lower, upper, 1.0e-9)))
    @test occursin("ScopeGuard[NarrowTransverseBound]", guard_message(() -> tube_ring_count(0.01, 2.0, 0.005)))
    segment = (start=[2.0, 0.0], stop=[5.0, 0.0])
    @test occursin("ScopeGuard[FootprintWithoutEdge]",
                   guard_message(() -> assert_etch_carries_edge([loop], segment, 1.0e-9)))
    # Exhibited classes from the frozen inputs: an upward exterior loop shows only
    # supported classes; a hole, a downward layer, rounding, a sloped sidewall, thin
    # metal and a missing trench show their guarded class.
    edges = [(slot=0, conductor=1), (slot=0, conductor=1)]
    layers = [(plane=0.0, sign=1)]
    @test exhibited_scope_classes(edges, [loop], layers, true, 90.0, 0.0, 0.0, 0.03, false, false) ==
          ["ExteriorLoops"]
    continued = (loop..., classes=vcat("Continuation", fill("Physical", length(points) - 1)))
    @test exhibited_scope_classes([(slot=0, conductor=1), (slot=1, conductor=2)], [continued],
                                  [(plane=0.0, sign=1), (plane=0.6, sign=1)], true, 90.0, 0.0, 0.0,
                                  0.03, true, true) ==
          ["ContinuationVertices", "DeviceFootprint", "ExteriorLoops", "MultipleConductors",
           "MultipleLayers", "MultipleSlots", "TraceBasis"]
    @test exhibited_scope_classes(edges, [loop, hole], layers, true, 90.0, 0.0, 0.0, 0.03, false, false) ==
          ["ExteriorLoops", "HoleLoops"]
    @test "HoleLoops" in RECIPE_SCOPE_SUPPORTED_CLASSES && "NarrowHoles" in ids
    @test exhibited_scope_classes(edges, [loop], [(plane=0.6, sign=-1)], true, 90.0, 0.0, 0.0, 0.03,
                                  false, false) == ["DownwardLayers", "ExteriorLoops"]
    @test "DownwardLayers" in RECIPE_SCOPE_SUPPORTED_CLASSES
    @test exhibited_scope_classes(edges, [loop], layers, false, 80.0, 0.005, 0.001, 0.0, false, false) ==
          ["ExteriorLoops", "NoTrench", "SlopedSidewalls", "ThinMetal", "TopRounding", "TrenchRounding"]
    # The metal loop records count the sides not on the box: the L-shaped loop has
    # three (its four other sides lie on the box faces).
    records = metal_loop_records([loop, hole], lower, upper, 1.0e-9)
    @test [record["Sides"] for record in records] == [3, 3]
    @test [record["Hole"] for record in records] == [false, true]
    @test records[1]["Vertices"] == 7 && records[1]["Conductor"] == 1 && records[1]["Plane"] == 0.0
    scope = recipe_scope_record(["ExteriorLoops"], [loop], lower, upper, 1.0e-9)
    @test scope["Recipe"] == "prism-tubes" && scope["GuardedClasses"] == ids
    @test scope["SupportedClasses"] == RECIPE_SCOPE_SUPPORTED_CLASSES
    @test [guard["Id"] for guard in scope["Guards"]] == ids
    @test 2 * sum(record["Sides"] for record in scope["MetalLoops"]) == 6
end

@testset "tube and band volume size laws" begin
    @test segment_point_distance(0.5, 1.0, 0.0, (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)) ≈ 1.0
    @test segment_point_distance(2.0, 0.0, 0.0, (0.0, 0.0, 0.0, 1.0, 0.0, 0.0)) ≈ 1.0
    record = prepare_tube_band_sizing!([([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])],
                                       [([0.0, 5.0, 0.0], [1.0, 5.0, 0.0])],
                                       0.04, 0.025, 0.16, 0.5)
    @test record["TubeAxisSegments"] == 1 && record["BandSegments"] == 1
    @test record["ProtectedDistance"] ≈ 0.05
    # Tube law: NormalSize up to the tube surface offset, then FarGrowth to FarSize.
    @test tube_band_size(0.5, 0.01, 0.0, 1.0) ≈ 0.025
    @test tube_band_size(0.5, 0.04, 0.0, 1.0) ≈ 0.025
    @test tube_band_size(0.5, 0.14, 0.0, 1.0) ≈ 0.025 + 0.5 * 0.1
    @test tube_band_size(0.5, 2.0, 0.0, 1.0) ≈ 0.16
    # The background size is never exceeded.
    @test tube_band_size(0.5, 2.0, 0.0, 0.1) ≈ 0.1
    # Band law: NormalSize + r inside the protected distance, then FarGrowth.
    @test tube_band_size(0.5, 5.0, 0.0, 1.0) ≈ 0.025
    @test tube_band_size(0.5, 5.0, 0.03, 1.0) ≈ 0.055
    @test tube_band_size(0.5, 5.0, 0.15, 1.0) ≈ 0.025 + 0.05 + 0.5 * 0.1
    @test tube_band_size(0.5, 5.0, 3.0, 1.0) ≈ 0.16
    @test band_law_size(0.03, 0.025, 0.16, 0.5) ≈ 0.055
    # Corner-ball exterior law (decision 39b): NormalSize at the ball radius, then
    # FarGrowth, up to FarSize; the interior is left to the background law.
    record = prepare_tube_band_sizing!([([0.0, 0.0, 0.0], [1.0, 0.0, 0.0])], [],
                                       0.04, 0.025, 0.16, 0.5, [(0.0, 10.0, 0.0)], 0.1)
    @test record["CornerExteriorPoints"] == 1 && record["CornerExteriorGrowth"] ≈ 0.5
    @test record["CornerIsotropyRadius"] ≈ 0.1 && record["TraceBasisVolumeGrowth"] ≈ 0.5
    @test tube_band_size(0.0, 10.05, 0.0, 1.0) ≈ 0.025
    @test tube_band_size(0.0, 10.1, 0.0, 1.0) ≈ 0.025
    @test tube_band_size(0.0, 10.3, 0.0, 1.0) ≈ 0.025 + 0.5 * 0.2
    @test tube_band_size(0.0, 11.0, 0.0, 1.0) ≈ 0.16
    @test corner_exterior_size(0.2, 0.025, 0.16, 0.5) ≈ 0.075
    empty!(TUBE_AXIS_SEGMENTS); empty!(BAND_SEGMENTS); empty!(CORNER_EXTERIOR_POINTS)
    @test tube_band_size(0.5, 5.0, 0.0, 1.0) ≈ 1.0
end

@testset "junction band: first-layer transverse size and footprint sides" begin
    # A junction line along x on the plane y = 0 (the cut face z = 0 carries the
    # matching label 1): two surface triangles with an edge on the line and heights
    # 0.025 (the band law) and 0.05 (the ridge grid), one triangle touching the
    # line at a vertex only (height 0.03), one triangle away from the line; two
    # tetrahedra on the line with apex heights 0.02 and 0.05.
    points = hcat([0.0, 0.0, 0.0], [0.025, 0.0, 0.0], [0.05, 0.0, 0.0],
                  [0.0125, 0.025, 0.0], [0.0375, 0.05, 0.0], [0.025, 0.03, 0.0],
                  [0.2, 0.2, 0.0], [0.3, 0.2, 0.0], [0.25, 0.3, 0.0],
                  [0.0125, 0.01, 0.02], [0.0375, 0.01, 0.05])
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("junction")
    surface = gmsh.model.addDiscreteEntity(2)
    volume = gmsh.model.addDiscreteEntity(3)
    gmsh.model.mesh.addNodes(2, surface, collect(1:size(points, 2)), vec(points))
    gmsh.model.mesh.addElementsByType(surface, 2, Int[],
                                      vec(hcat([1, 2, 4], [2, 3, 5], [2, 6, 3], [7, 8, 9])))
    gmsh.model.mesh.addElementsByType(volume, 4, Int[], vec(hcat([1, 2, 4, 10], [2, 3, 5, 11])))
    gmsh.model.addPhysicalGroup(2, [surface], 1)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    xyz = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    segments = [([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0])]
    surface_layer = first_layer_transverse_statistics(xyz, index, segments, 0.025, 2, 1)
    @test surface_layer["Elements"] == 3 && surface_layer["Prescribed"] == 0.025
    @test surface_layer["TransverseP50"] ≈ 0.03 && surface_layer["TransverseMaximum"] ≈ 0.05
    @test surface_layer["TransverseP10"] ≈ 0.025
    @test surface_layer["AchievedOverPrescribedP50"] ≈ 1.2
    @test surface_layer["AchievedOverPrescribedP90"] ≈ 2.0
    tet_layer = first_layer_transverse_statistics(xyz, index, segments, 0.025, 3, 0)
    # The tetrahedra's largest node distances to the line y = z = 0 are 0.025 (the
    # surface node) and sqrt(0.01^2 + 0.05^2) (the apex).
    @test tet_layer["Elements"] == 2
    @test tet_layer["TransverseP50"] ≈ 0.025 && tet_layer["TransverseMaximum"] ≈ hypot(0.01, 0.05)
    # A line no element touches, and a line whose segment ends before the elements
    # (the node on the line must lie within the segment span).
    @test first_layer_transverse_statistics(xyz, index, [([0.0, 1.0, 0.0], [1.0, 1.0, 0.0])],
                                            0.025, 2, 1)["Elements"] == 0
    @test first_layer_transverse_statistics(xyz, index, [([-1.0, 0.0, 0.0], [-0.5, 0.0, 0.0])],
                                            0.025, 2, 1)["Elements"] == 0
    # The band tetrahedron statistic reports the segment count.
    band = band_tetrahedron_statistics(xyz, index, segments, 0.025)
    @test band["Cells"] == 2 && band["Segments"] == 1
    # Shell statistics against a law: both tetrahedra have their centroid within
    # the first NormalSize shell of the line (0.009 and 0.015); a constant law of
    # 0.025 gives the ratio mean edge / 0.025; an empty shell reports nothing.
    shells = shell_size_statistics(xyz, index, c -> hypot(c[2], c[3]), 0.0, [0.025, 0.05],
                                   (c, d) -> 0.025)
    @test length(shells) == 2 && shells[1]["Cells"] == 2 && shells[2]["Cells"] == 0
    @test shells[1]["Lower"] == 0.0 && shells[1]["Upper"] == 0.025 && shells[2]["Lower"] == 0.025
    mean_edges = sort([sum(norm(points[:, i] .- points[:, j]) for (i, j) in
                           ((a, b) for a in cell for b in cell if a < b)) / 6
                       for cell in ([1, 2, 4, 10], [2, 3, 5, 11])])
    @test shells[1]["MeanEdgeP50"] ≈ mean_edges[1] && shells[1]["MeanEdgeP90"] ≈ mean_edges[2]
    @test shells[1]["AchievedOverPrescribedP50"] ≈ mean_edges[1] / 0.025
    @test shells[1]["LongestEdgeMaximum"] ≈ maximum(norm(points[:, i] .- points[:, j])
                                                    for cell in ([1, 2, 4, 10], [2, 3, 5, 11])
                                                    for i in cell for j in cell)
    @test shells[2]["MeanEdgeP50"] === nothing
    gmsh.finalize()
    # Footprint sides on the outer box are separated from the interior (feature) sides.
    polygons = [Dict{String, Any}("Plane" => -0.05,
                                  "Points" => [[-6.0, -8.0], [-6.0, 8.0], [2.0, 8.0], [2.0, -8.0]])]
    interior, on_box = footprint_polygon_segments(polygons, [-6.0, -8.0, -1.0], [10.0, 8.0, 1.0], 1.0e-9)
    @test length(on_box) == 3 && length(interior) == 1
    @test interior[1] == ([2.0, 8.0, -0.05], [2.0, -8.0, -0.05])
end
