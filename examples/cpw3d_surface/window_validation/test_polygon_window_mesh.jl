# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests of the fabricated window mesher:
#   julia --project=. test_polygon_window_mesh.jl
# Schema checks, the z-level stack of the transmon reference (29 levels of the recorded
# r50 / t5 manifest), input refusals (overlapping polygons, a bump off metal), and the tiny
# synthetic two-level window meshed at r 0.2 um / t 5 um: every attribute area and volume
# against its analytic value, then validate_window_mesh.jl (adjacency rules, positivity);
# the boundary-layer Thickness margin (decision 188) and its explicit exact-sum option; the
# cross-plane reconciliation of coincident two-plane edges (snap rule, order independence,
# tie-breaking, the same-plane refusal).

using Test
using JSON

include(joinpath(@__DIR__, "PolygonWindowMesh.jl"))
using .PolygonWindowMesh
include(joinpath(@__DIR__, "synthetic_two_level_window.jl"))

const PWM = PolygonWindowMesh

@testset "polygon set schema" begin
    spec = read_polygon_set(synthetic_two_level_window())
    @test length(spec.planes) == 2
    @test spec.planes[1].facing == 1 && spec.planes[2].facing == -1
    @test spec.terminals == ["trace_l1", "trace_l2"]
    @test PWM.attribute_table(spec) ==
          Dict("ground" => (4, 5), "trace_l1" => (7, 8), "trace_l2" => (10, 11))
    names = Dict(name => (d, a) for (d, a, name) in PWM.physical_names(spec))
    @test names["trace_l2_substrate"] == (2, 11) && names["substrate_backside"] == (2, 9)
    bad = synthetic_two_level_window()
    bad["Terminals"] = ["trace_l1"]
    @test_throws ErrorException read_polygon_set(bad)
    swapped = synthetic_two_level_window()
    swapped["Planes"][2]["Facing"] = "up"
    @test_throws ErrorException read_polygon_set(swapped)
end

@testset "z levels of the transmon reference" begin
    spec = read_polygon_set(
        Dict(
            "Box" => Dict("X" => [-2000.0, 2000.0], "Y" => [-712.0, 2988.0]),
            "Planes" => [
                Dict(
                    "Name" => "L1",
                    "SurfaceZ" => 0.0,
                    "Facing" => "up",
                    "SubstrateThickness" => 525.00000015,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [
                                [-2000, -712],
                                [2000, -712],
                                [2000, 2988],
                                [-2000, 2988]
                            ]
                        )
                    ]
                )
            ],
            "Vacuum" => Dict("Below" => 475.0, "Above" => 1000.00000015)
        )
    )
    stack = PWM.z_levels(spec, 2, 1)
    recorded = [
        -1000.00000015,
        -525.00000015,
        -100.0,
        -50.0,
        -20.0,
        -10.0,
        -5.0,
        -2.0,
        -1.0,
        -0.5,
        -0.2,
        -0.1,
        -0.05,
        0.0,
        0.05,
        0.1,
        0.15,
        0.2,
        0.3,
        0.5,
        1.0,
        2.0,
        5.0,
        10.0,
        20.0,
        50.0,
        100.0,
        200.0,
        1000.00000015
    ]
    @test length(stack.levels) == 29
    @test all(isapprox.(stack.levels, recorded; atol=1.0e-9))
    @test stack.backsides == [-525.00000015]
    two = read_polygon_set(synthetic_two_level_window())
    gap = PWM.z_levels(two, 2, 1)
    @test gap.levels[1] == -20.0 && gap.levels[end] == 24.8
    @test isempty(gap.backsides)
    # Both metal bands and trench bands are levels; the gap is graded from both sides.
    for z in (-0.05, 0.0, 0.05, 0.1, 4.7, 4.75, 4.8, 4.85, 2.4)
        @test any(isapprox.(gap.levels, z; atol=1.0e-9))
    end
    @test PWM.material(two, PWM.PartitionClass(["ground", "ground"], 1), 2.4) == 0
    @test PWM.material(two, PWM.PartitionClass(["ground", "ground"], 0), 2.4) == 2
    @test PWM.material(two, PWM.PartitionClass(["", "ground"], 0), -0.025) == 2
    @test PWM.material(two, PWM.PartitionClass(["ground", ""], 0), -0.025) == 1
    @test PWM.material(two, PWM.PartitionClass(["", ""], 0), 4.825) == 2
    @test PWM.material(two, PWM.PartitionClass(["", "ground"], 0), 4.825) == 1
    @test PWM.material(two, PWM.PartitionClass(["", "ground"], 0), 4.75) == 0
end

# Two planes whose ground edges are nominally coincident along y = 10 (different vertex
# samplings) and share a chamfer sampled as one chord (L1) versus two chords (L2, the middle
# vertex 0.03 um off the chord); L2's vertices are listed in a different order.
function coincident_edge_set(; l2_rotation=0, l2_polygon_order=[1, 2])
    l2_ground = [
        [0.0, 0.0],
        [40.0, 0.0],
        [40.0, 10.0],
        [25.0, 10.0],
        [15.0, 10.0],
        [10.0, 10.0],
        [5.0, 12.5 + 0.03],
        [0.0, 15.0]
    ]
    l2_ground = circshift(l2_ground, l2_rotation)
    l2_polygons = [
        Dict("Conductor" => "ground", "Outer" => l2_ground),
        Dict(
            "Conductor" => "island",
            "Outer" => [[10.0, 20.0], [30.0, 20.0], [30.0, 25.0], [10.0, 25.0]]
        )
    ][l2_polygon_order]
    return Dict(
        "Version" => 1,
        "Name" => "coincident",
        "MatchingRadius" => 1.9,
        "Box" => Dict("X" => [0.0, 40.0], "Y" => [0.0, 30.0]),
        "Planes" => [
            Dict(
                "Name" => "L1",
                "SurfaceZ" => 0.0,
                "Facing" => "up",
                "SubstrateThickness" => 20.0,
                "Polygons" => [
                    Dict(
                        "Conductor" => "ground",
                        "Outer" => [
                            [0.0, 0.0],
                            [40.0, 0.0],
                            [40.0, 10.0],
                            [30.0, 10.0],
                            [20.0, 10.0 + 0.02],
                            [10.0, 10.0],
                            [0.0, 15.0]
                        ]
                    )
                ]
            ),
            Dict(
                "Name" => "L2",
                "SurfaceZ" => 4.8,
                "Facing" => "down",
                "SubstrateThickness" => 20.0,
                "Polygons" => l2_polygons
            )
        ]
    )
end

@testset "cross-plane reconciliation" begin
    spec = read_polygon_set(coincident_edge_set())
    delta, rule = PWM.cross_plane_snap_distance(spec, NaN)
    @test delta ≈ 0.0475 && rule == "0.025 x MatchingRadius"
    @test PWM.cross_plane_snap_distance(spec, 0.1) == (0.1, "override")
    @test PWM.cross_plane_snap_distance(
        read_polygon_set(synthetic_two_level_window()),
        NaN
    )[1] ≈ 0.0475
    reconciled, report = reconcile_planes(spec, delta)
    l1 = reconciled.planes[1].polygons[1].outer
    l2 = reconciled.planes[2].polygons[1].outer
    # The L1 vertex (20, 10.02) is fixed (the lower plane never moves); L2's (25, 10) and
    # (15, 10) snap onto the L1 segments through it (0.01 um off the line y = 10) and L2's
    # chamfer midpoint snaps onto L1's chord (its foot, 0.027 um away); then L1 receives
    # L2's three snapped vertices exactly and L2 receives L1's (30, 10) and (20, 10.02)
    # within delta, so the run (40, 10) -> (0, 15) is one identical point sequence.
    @test report["applied"] && report["moved_vertices"] == 3
    @test report["max_displacement_um"] ≈ 0.03 / hypot(2.0, 1.0) * 2.0 atol = 1.0e-6
    @test report["inserted_vertices"] == Dict("L1" => 3, "L2" => 2)
    run(points) = points[findfirst(==((40.0, 10.0)), points):end]
    @test run(l1) == run(l2)
    @test length(run(l1)) == 8 && (30.0, 10.0) in run(l2) && (20.0, 10.02) in run(l2)
    @test all(p -> abs(p[2] - 10.0) <= 0.02 + 1.0e-12, run(l2)[2:6])
    @test report["coincident_segments"] == 7
    @test report["coincident_run_length_um"] ≈ 30.0 + hypot(10.0, 5.0) rtol = 1.0e-3
    # Order independence: L2's loop rotated and its polygons permuted give the same chains.
    for (rotation, order) in ((3, [1, 2]), (5, [2, 1]))
        other, other_report = reconcile_planes(
            read_polygon_set(
                coincident_edge_set(; l2_rotation=rotation, l2_polygon_order=order)
            ),
            delta
        )
        @test Set(other.planes[1].polygons[1].outer) == Set(l1)
        ground = other.planes[2].polygons[order[1] == 1 ? 1 : 2].outer
        @test Set(ground) == Set(l2)
        @test other_report["moved_vertices"] == 3 &&
              other_report["coincident_segments"] == 7
    end
    # Tie: an L2 vertex equidistant from two L1 vertices goes to the smaller coordinates.
    tie = read_polygon_set(
        Dict(
            "Version" => 1,
            "MatchingRadius" => 1.9,
            "Box" => Dict("X" => [0.0, 10.0], "Y" => [0.0, 10.0]),
            "Planes" => [
                Dict(
                    "Name" => "L1",
                    "SurfaceZ" => 0.0,
                    "Facing" => "up",
                    "SubstrateThickness" => 5.0,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0.0, 0.0], [10.0, 0.0], [10.0, 2.0], [0.0, 2.0]]
                        ),
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0.0, 2.04], [10.0, 2.04], [10.0, 4.0], [0.0, 4.0]]
                        )
                    ]
                ),
                Dict(
                    "Name" => "L2",
                    "SurfaceZ" => 4.8,
                    "Facing" => "down",
                    "SubstrateThickness" => 5.0,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0.0, 2.02], [10.0, 2.02], [10.0, 6.0], [0.0, 6.0]]
                        )
                    ]
                )
            ]
        )
    )
    # The same-plane slot of 0.04 um is refused at delta 0.0475 ...
    @test_throws ErrorException reconcile_planes(tie, 0.0475)
    # ... and with a smaller override the equidistant L2 vertices take the lower L1 edge.
    snapped, _ = reconcile_planes(tie, 0.03)
    @test snapped.planes[2].polygons[1].outer[1:2] == [(0.0, 2.0), (10.0, 2.0)]
    # A two-plane set without MatchingRadius and without an override is refused.
    bare = coincident_edge_set()
    delete!(bare, "MatchingRadius")
    @test_throws ErrorException PWM.cross_plane_snap_distance(read_polygon_set(bare), NaN)
end

@testset "input refusals" begin
    overlapping = synthetic_two_level_window()
    push!(
        overlapping["Planes"][1]["Polygons"],
        Dict(
            "Conductor" => "extra",
            "Outer" => rectangle(50.0, 60.0, 5.0, 15.0),
            "Holes" => []
        )
    )
    overlapping["Terminals"] = ["extra", "trace_l1", "trace_l2"]
    @test_throws ErrorException mesh_polygon_window(
        read_polygon_set(overlapping),
        0.2,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
    off_metal = synthetic_two_level_window()
    off_metal["Bumps"][1]["Footprint"] = regular_polygon(42.0, 25.0, 2.0, 8)
    @test_throws ErrorException mesh_polygon_window(
        read_polygon_set(off_metal),
        0.2,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
end

@testset "synthetic two-level window" begin
    directory = mktempdir()
    output = joinpath(directory, "synthetic_r200.msh2")
    spec = read_polygon_set(synthetic_two_level_window())
    manifest = mesh_polygon_window(spec, 0.2, 5.0, output; verbose=false)
    @test manifest["nonpositive"] == 0
    # Decision 188: the band Thickness carries the 1e-6 margin by default (3 layers of
    # first height 0.2: exact sum 1.4); the exact-sum mode is explicit and recorded.
    @test manifest["radial_layers"] == 3
    @test manifest["radial_band_thickness_mode"] == "geometric_sum_x_1p000001"
    @test manifest["radial_band_thickness_um"] ≈ 1.4 * (1.0 + 1.0e-6) rtol = 1.0e-12
    exact = mesh_polygon_window(
        spec,
        0.2,
        5.0,
        joinpath(directory, "synthetic_r200_exact.msh2");
        verbose=false,
        plan_only=true,
        exact_band_thickness=true
    )
    @test exact["radial_band_thickness_mode"] == "exact_geometric_sum"
    @test exact["radial_band_thickness_um"] ≈ 1.4 rtol = 1.0e-12
    @test exact["plan_perimeter_edges"] == manifest["plan_perimeter_edges"]
    @test manifest["first_layer_normal_outliers_above_1p5x_target"] == 0
    @test manifest["first_layer_normal_height_um"]["maximum"] <= 0.22
    @test manifest["perimeter_tangent_length_um"]["maximum"] <= 5.0 + 1.0e-9
    @test manifest["conductor_components"]["L1"]["trace_l1"]["components"] == 1
    @test manifest["conductor_components"]["L2"]["ground"]["components"] == 1
    # Analytic areas: box 80 x 60, gap holes 42 x 16, traces 30 x 4, bump 16-gon of radius 4,
    # metal 0.1, overetch 0.05, gap 4.8 (bump column 4.6), substrates 20 each side.
    box, hole, trace = 4800.0, 672.0, 120.0
    ground = box - hole
    bump_area = 0.5 * 16 * 16.0 * sin(2pi / 16)
    bump_perimeter = 16 * 8.0 * sin(pi / 16)
    hole_perimeter, trace_perimeter = 2 * (42.0 + 16.0), 2 * (30.0 + 4.0)
    areas = manifest["surface_area_um2"]
    @test areas["5"] ≈ 2ground rtol = 1.0e-9
    @test areas["4"] ≈
          2 * (ground - bump_area) + 2 * hole_perimeter * 0.1 + bump_perimeter * 4.6 rtol =
        1.0e-9
    @test areas["7"] ≈ trace + trace_perimeter * 0.1 rtol = 1.0e-9
    @test areas["8"] ≈ trace rtol = 1.0e-9
    @test areas["10"] ≈ trace + trace_perimeter * 0.1 rtol = 1.0e-9
    @test areas["11"] ≈ trace rtol = 1.0e-9
    @test areas["6"] ≈ 2 * (hole - trace) + 2 * (hole_perimeter + trace_perimeter) * 0.05 rtol =
        1.0e-9
    @test areas["3"] ≈ 2 * (80.0 + 60.0) * (44.8 - 0.2) + 2box rtol = 1.0e-9
    @test !haskey(areas, "9")
    volumes = manifest["volume_um3"]
    @test volumes["1"] ≈ 2 * (box * 20.0 - (hole - trace) * 0.05) rtol = 1.0e-9
    @test volumes["2"] ≈
          4.6box + 2 * 0.1 * (hole - trace) + 2 * 0.05 * (hole - trace) - bump_area * 4.6 rtol =
        1.0e-9
    validator = joinpath(@__DIR__, "validate_window_mesh.jl")
    run(
        pipeline(
            `$(Base.julia_cmd()) --project=$(@__DIR__) $validator $output`;
            stdout=devnull
        )
    )
    validation = JSON.parsefile(replace(output, r"\.msh2$" => ".validation.json"))
    @test all(iszero, values(validation["surface_adjacency_errors"]))
    @test validation["mesh_quality"]["nonpositive"] == 0
    @test validation["nodes"] == manifest["nodes"]
    @test validation["surface_area_um2"]["4"] ≈ areas["4"] rtol = 1.0e-9
    @test Set(keys(validation["surface_adjacency_errors"])) == Set([
        "exterior_boundary",
        "ground_air",
        "ground_substrate",
        "substrate_air",
        "trace_l1_air",
        "trace_l1_substrate",
        "trace_l2_air",
        "trace_l2_substrate"
    ])
end
