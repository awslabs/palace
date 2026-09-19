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
    tube = EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 1.01, 0.05)
    @test tube.layers == 21 && tube_spacing(tube) ≈ 1.01 / 21 && tube_spacing(tube) <= 0.05
    @test_throws ErrorException EdgeTube([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0],
                                         0.0, 1.0, 0.0)
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
    empty!(TUBE_AXIS_SEGMENTS); empty!(BAND_SEGMENTS)
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
    gmsh.finalize()
    # Footprint sides on the outer box are separated from the interior (feature) sides.
    polygons = [Dict{String, Any}("Plane" => -0.05,
                                  "Points" => [[-6.0, -8.0], [-6.0, 8.0], [2.0, 8.0], [2.0, -8.0]])]
    interior, on_box = footprint_polygon_segments(polygons, [-6.0, -8.0, -1.0], [10.0, 8.0, 1.0], 1.0e-9)
    @test length(on_box) == 3 && length(interior) == 1
    @test interior[1] == ([2.0, 8.0, -0.05], [2.0, -8.0, -0.05])
end
