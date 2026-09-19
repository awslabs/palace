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
    hole = (conductor=1, plane=0.0, hole=true, points=points, classes=loop.classes)
    @test_throws ErrorException metal_edge_segments([hole], corners, clearance, lower, upper, 1.0e-9)
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
