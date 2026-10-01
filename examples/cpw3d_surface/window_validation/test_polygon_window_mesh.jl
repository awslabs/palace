# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests of the fabricated window mesher:
#   julia --project=. test_polygon_window_mesh.jl
# Schema checks, the z-level stack of the transmon reference (29 levels of the recorded
# r50 / t5 manifest), input refusals (overlapping polygons, a bump off metal), and the tiny
# synthetic two-level window meshed at r 0.2 um / t 5 um: every attribute area and volume
# against its analytic value, then validate_window_mesh.jl (adjacency rules, positivity);
# the boundary-layer Thickness margin (decision 188) and its explicit exact-sum option.

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
