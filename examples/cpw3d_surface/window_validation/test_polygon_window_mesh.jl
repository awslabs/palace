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
# tie-breaking, the same-plane refusal); the own structured band (facing distance per segment
# and side with other bodies' corners, row termination, fans, inward corners and T junctions,
# wall ends, the collision refusal) and its first-layer / area exactness on slot windows.

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
    @test delta ≈ 0.095 && rule == "0.05 x MatchingRadius"
    @test PWM.cross_plane_snap_distance(spec, 0.1) == (0.1, "override")
    @test PWM.cross_plane_snap_distance(
        read_polygon_set(synthetic_two_level_window()),
        NaN
    )[1] ≈ 0.095
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
    # The same-plane slot of 0.04 um is refused at delta 0.095 ...
    @test_throws ErrorException reconcile_planes(tie, 0.095)
    # ... and with a smaller override the equidistant L2 vertices take the lower L1 edge.
    snapped, _ = reconcile_planes(tie, 0.03)
    @test snapped.planes[2].polygons[1].outer[1:2] == [(0.0, 2.0), (10.0, 2.0)]
    # A two-plane set without MatchingRadius and without an override is refused.
    bare = coincident_edge_set()
    delete!(bare, "MatchingRadius")
    @test_throws ErrorException PWM.cross_plane_snap_distance(read_polygon_set(bare), NaN)
    # No window-cut sliver configuration here: every rule counter is zero.
    rules = report["window_cut_sliver_rules"]
    @test rules["wall_constrained_snaps"] == 0 &&
          rules["wall_insertions_refused"] == 0 &&
          rules["single_segment_insertions"] == 0 &&
          rules["merged_consecutive_vertices"] == 0 &&
          rules["collapsed_spikes"] == 0
end

# Two planes for the window-cut sliver rules (decision 306). The C3 class: an L2 finger whose
# top edge y = 10.08 lies within delta of L1's ground edge y = 10 (both finger corners snap onto
# it) while L1 carries the vertex P = (19.93, 10) on that run, 0.07 um from the snapped finger
# corner B = (20, 10): P lies ON the snapped top edge and 0.07 um from the finger's (slightly
# tilted) side, so an insertion into both segments would give the spike P, B, P.
function finger_spike_set()
    return Dict(
        "Version" => 1,
        "Name" => "finger-spike",
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
                            [19.93, 10.0],
                            [0.0, 10.0]
                        ]
                    )
                ]
            ),
            Dict(
                "Name" => "L2",
                "SurfaceZ" => 4.8,
                "Facing" => "down",
                "SubstrateThickness" => 20.0,
                "Polygons" => [
                    Dict(
                        "Conductor" => "finger",
                        "Outer" => [
                            [15.0, 5.0],
                            [19.98, 5.0],
                            [20.0, 10.08],
                            [15.0, 10.08]
                        ]
                    )
                ]
            )
        ]
    )
end

# The C1 class: L2 ground covers the top-left box corner along the top wall y = 30 (a
# window-cut line); L1 ground is a wedge whose edge meets the wall 0.24 um from the corner with
# an 18-nm first chord to the chip node (0.25, 29.985), 15 nm below the wall and within delta
# of L2's wall segment (an insertion would bend the cut line into a 2 r x 1.5 r notch at r10);
# L2's wall vertex (0.3, 30) is nearer to that off-wall node (0.052 um) than to L1's wall
# vertex (0.24, 30) (0.06 um): the wall rule keeps it on the wall.
function wall_notch_set()
    return Dict(
        "Version" => 1,
        "Name" => "wall-notch",
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
                            [0.24, 30.0],
                            [0.25, 29.985],
                            [1.43, 28.28],
                            [4.0, 28.28],
                            [4.0, 30.0]
                        ]
                    )
                ]
            ),
            Dict(
                "Name" => "L2",
                "SurfaceZ" => 4.8,
                "Facing" => "down",
                "SubstrateThickness" => 20.0,
                "Polygons" => [
                    Dict(
                        "Conductor" => "ground",
                        "Outer" => [
                            [0.0, 20.0],
                            [40.0, 20.0],
                            [40.0, 30.0],
                            [0.3, 30.0],
                            [0.0, 30.0]
                        ]
                    )
                ]
            )
        ]
    )
end

@testset "window-cut sliver rules (decision 306)" begin
    tolerance = PWM.PLAN_NODE_MERGE_TOLERANCE_UM
    # clean_loop: consecutive vertices within the identity quantum merge onto the first, a
    # spike collapses onto its neighbours, a self-touching loop is refused with the point.
    loop = PWM.Point2[(0.0, 0.0), (10.0, 0.0), (10.0, 1.0e-7), (10.0, 10.0), (0.0, 10.0)]
    cleaned, merged, spikes = PWM.clean_loop(loop, tolerance, "test")
    @test cleaned == [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)] &&
          merged == 1 &&
          spikes == 0
    spike = PWM.Point2[
        (0.0, 0.0),
        (10.0, 0.0),
        (5.0, 10.0),
        (4.0, 10.0),
        (5.0, 10.0),
        (0.0, 10.0)
    ]
    cleaned, merged, spikes = PWM.clean_loop(spike, tolerance, "test")
    @test cleaned == [(0.0, 0.0), (10.0, 0.0), (5.0, 10.0), (0.0, 10.0)] &&
          merged == 0 &&
          spikes == 1
    # A spike whose return point is within the quantum of the departure point, then a merge.
    nearly = PWM.Point2[
        (0.0, 0.0),
        (10.0, 0.0),
        (5.0, 10.0),
        (4.0, 10.0),
        (5.0 + 1.0e-7, 10.0),
        (0.0, 10.0)
    ]
    cleaned, merged, spikes = PWM.clean_loop(nearly, tolerance, "test")
    @test cleaned == [(0.0, 0.0), (10.0, 0.0), (5.0, 10.0), (0.0, 10.0)] && spikes == 1
    eight =
        PWM.Point2[(0.0, 0.0), (2.0, 0.0), (1.0, 1.0), (2.0, 2.0), (0.0, 2.0), (1.0, 1.0)]
    message = try
        PWM.clean_loop(eight, tolerance, "Plane L2 loop 1")
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("Plane L2 loop 1 visits (1.0, 1.0) twice", message)

    # C3 class: P is inserted into the nearest segment only (the snapped top edge), never into
    # the finger's side as well; the loops are simple and the finger corner stays at B.
    spec = read_polygon_set(finger_spike_set())
    delta, _ = PWM.cross_plane_snap_distance(spec, NaN)
    reconciled, report = reconcile_planes(spec, delta)
    rules = report["window_cut_sliver_rules"]
    finger = reconciled.planes[2].polygons[1].outer
    # The snapped finger corner is the segment foot (20 + 4e-15, 10): compared to 1e-12.
    same_points(a, b) =
        length(a) == length(b) && all(
            isapprox(p[1], q[1]; atol=1.0e-12) && isapprox(p[2], q[2]; atol=1.0e-12) for
            (p, q) in zip(a, b)
        )
    @test same_points(
        finger,
        [(15.0, 5.0), (19.98, 5.0), (20.0, 10.0), (19.93, 10.0), (15.0, 10.0)]
    )
    @test length(unique(finger)) == length(finger)
    @test same_points(
        reconciled.planes[1].polygons[1].outer,
        [
            (0.0, 0.0),
            (40.0, 0.0),
            (40.0, 10.0),
            (20.0, 10.0),
            (19.93, 10.0),
            (15.0, 10.0),
            (0.0, 10.0)
        ]
    )
    @test report["moved_vertices"] == 2 && report["max_displacement_um"] ≈ 0.08
    @test rules["single_segment_insertions"] == 1 &&
          rules["collapsed_spikes"] == 0 &&
          rules["merged_consecutive_vertices"] == 0 &&
          rules["wall_insertions_refused"] == 0
    @test report["inserted_vertices"] == Dict("L1" => 2, "L2" => 1)
    manifest = mesh_polygon_window(
        spec,
        0.01,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
    @test manifest["cross_plane_reconciliation"]["window_cut_sliver_rules"]["single_segment_insertions"] ==
          1
    @test manifest["first_layer_normal_outliers_above_1p5x_target"] == 0

    # C1 class: the off-wall chip node is not inserted into L2's wall segment, L2's wall
    # vertex snaps along the wall onto L1's wall vertex, and the plan (18-nm chord 0.24 um from
    # the box corner, at r10) meshes.
    spec = read_polygon_set(wall_notch_set())
    reconciled, report = reconcile_planes(spec, delta)
    rules = report["window_cut_sliver_rules"]
    l2 = reconciled.planes[2].polygons[1].outer
    @test all(p -> p[2] == 30.0, filter(p -> p[2] > 29.0, l2))
    @test (0.24, 30.0) in l2 && (4.0, 30.0) in l2 && !((0.3, 30.0) in l2)
    @test !((0.25, 29.985) in l2)
    @test rules["wall_insertions_refused"] == 1 &&
          rules["wall_constrained_snaps"] == 1 &&
          rules["single_segment_insertions"] == 0
    @test report["moved_vertices"] == 1 && report["max_displacement_um"] ≈ 0.06
    manifest = mesh_polygon_window(
        spec,
        0.01,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
    @test manifest["cross_plane_reconciliation"]["window_cut_sliver_rules"]["wall_insertions_refused"] ==
          1
    @test manifest["perimeter_tangent_length_um"]["minimum"] ≈ hypot(0.01, 0.015) atol =
        1.0e-9
    @test manifest["first_layer_normal_outliers_above_1p5x_target"] == 0
    # Wall ends: the wedge chain (both sides) and L2's edge y = 20 (both sides), two each.
    @test manifest["band"]["wall_end_columns"] == 8
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
    # The own band: exactly 3 rows of first height 0.2 (band 1.4) by construction.
    @test manifest["radial_layers"] == 3
    @test manifest["radial_band_thickness_mode"] == "own_structured_band"
    @test manifest["radial_band_thickness_um"] ≈ 1.4 rtol = 1.0e-12
    # Inside the 16-gon bump (radius 4, chords 1.56 um) the facing rule sees the next-but-one
    # chord as a front 1.56 um away: 2 rows on the bump's inner band, 3 everywhere else.
    @test manifest["band"]["mode"] == "own" && manifest["band"]["applied_minimum_rows"] == 2
    @test manifest["band"]["per_side_segment_length_um_by_rows"]["2"] ≈
          16 * 8.0 * sin(pi / 16) rtol = 1.0e-9
    @test manifest["first_layer_normal_height_um"]["minimum"] ≈ 0.2 rtol = 1.0e-9
    @test manifest["first_layer_normal_height_um"]["maximum"] ≈ 0.2 rtol = 1.0e-9
    # Gmsh mode (decision 188): the band Thickness carries the 1e-6 margin by default; the
    # exact-sum mode is explicit and recorded.
    margin = mesh_polygon_window(
        spec,
        0.2,
        5.0,
        joinpath(directory, "synthetic_r200_margin.msh2");
        verbose=false,
        plan_only=true,
        band_mode=:gmsh
    )
    @test margin["radial_band_thickness_mode"] == "geometric_sum_x_1p000001"
    @test margin["radial_band_thickness_um"] ≈ 1.4 * (1.0 + 1.0e-6) rtol = 1.0e-12
    exact = mesh_polygon_window(
        spec,
        0.2,
        5.0,
        joinpath(directory, "synthetic_r200_exact.msh2");
        verbose=false,
        plan_only=true,
        band_mode=:gmsh,
        exact_band_thickness=true
    )
    @test exact["radial_band_thickness_mode"] == "exact_geometric_sum"
    @test exact["radial_band_thickness_um"] ≈ 1.4 rtol = 1.0e-12
    @test exact["plan_perimeter_edges"] == manifest["plan_perimeter_edges"]
    @test margin["plan_perimeter_edges"] == manifest["plan_perimeter_edges"]
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

# ---------------------------------------------------------------------------------------------
# Own structured band (structured_band.jl).

single_plane_set(name, outer; box_x=[0.0, 20.0], box_y=[0.0, 80.0]) = Dict(
    "Version" => 1,
    "Name" => name,
    "Box" => Dict("X" => box_x, "Y" => box_y),
    "Planes" => [
        Dict(
            "Name" => "L1",
            "SurfaceZ" => 0.0,
            "Facing" => "up",
            "SubstrateThickness" => 20.0,
            "Polygons" => [Dict("Conductor" => "ground", "Outer" => outer)]
        )
    ],
    "Vacuum" => Dict("Below" => 0.0, "Above" => 20.0)
)

@testset "facing distance per segment and side" begin
    # An L1 edge along y = 0 (points 1-2); a 2-um finger tip at y = 0.75 over x in [5, 7] with
    # perpendicular sides (points 3-6); a collinear continuation (points 7-8); a curve sharing
    # point 2 (points 2-9).
    curves = [
        PWM.PlanCurve(1, (0.0, 0.0), (10.0, 0.0), (1, 2), true),
        PWM.PlanCurve(2, (5.0, 0.75), (7.0, 0.75), (3, 4), true),
        PWM.PlanCurve(3, (5.0, 0.75), (5.0, 5.0), (3, 5), true),
        PWM.PlanCurve(4, (7.0, 0.75), (7.0, 5.0), (4, 6), true),
        PWM.PlanCurve(5, (12.0, 0.0), (20.0, 0.0), (7, 8), true),
        PWM.PlanCurve(6, (10.0, 0.0), (10.0, 0.3), (2, 9), true)
    ]
    d(q1, q2, normal) = PWM.segment_facing_distance(q1, q2, normal, curves[1], curves, 10.0)
    # The segment ending at (5, 0) sees the finger corner (review M1: the corner shadow).
    @test d((4.0, 0.0), (5.0, 0.0), (0.0, 1.0)) ≈ 0.75
    # The segment before it sees the corner diagonally (continuous, not Inf).
    @test d((3.0, 0.0), (4.0, 0.0), (0.0, 1.0)) ≈ hypot(1.0, 0.75)
    # Under the tip: the exact segment distance 0.75 (per segment, no sampling).
    @test d((5.5, 0.0), (6.5, 0.0), (0.0, 1.0)) ≈ 0.75
    # The other side of the edge faces nothing (per side, m8). The last segment sees the
    # finger's far side (2.136 um), not the corner-adjacent curve at point 2 (distance 0) nor
    # the collinear continuation (2.0 um): both excluded.
    @test d((5.0, 0.0), (6.0, 0.0), (0.0, -1.0)) == Inf
    @test d((9.0, 0.0), (10.0, 0.0), (0.0, 1.0)) ≈ hypot(2.0, 0.75)
    @test PWM.band_cap_rows(0.75, 0.01, 7) == 4
    @test PWM.band_cap_rows(hypot(1.0, 0.75), 0.01, 7) == 5
    @test PWM.band_cap_rows(Inf, 0.01, 7) == 7
    @test PWM.band_cap_rows(0.01, 0.01, 7) == 1
end

# A straight open chain with its interior above, wall ends perpendicular at both ends.
function straight_chain(bases)
    n = length(bases)
    return PWM.MetalChain(
        false,
        bases,
        fill(1, n - 1),
        fill(false, n),
        (0.0, 1.0),
        (0.0, 1.0),
        1,
        2
    )
end

# `distances`: the facing distance per segment (default: nothing faces the chain).
function band_of(
    chain,
    segment_rows,
    heights;
    fan_turn=deg2rad(90.0),
    distances=fill(Inf, length(segment_rows))
)
    statistics = PWM.BandStatistics()
    columns, _, node_rows =
        PWM.chain_columns(chain, segment_rows, distances, heights, fan_turn, statistics)
    triangles = NTuple{3, PWM.Point2}[]
    for c = 1:(length(columns) - 1)
        PWM.push_band_elements!(triangles, columns[c], columns[c + 1], statistics, "test")
    end
    chain.closed &&
        PWM.push_band_elements!(triangles, columns[end], columns[1], statistics, "test")
    return columns, node_rows, triangles, statistics
end

near(p, q; atol=1.0e-12) = hypot(p[1] - q[1], p[2] - q[2]) <= atol

@testset "row termination by collapse" begin
    heights = PWM.band_heights(0.01, 2.0, 7)
    @test heights ≈ [0.01, 0.03, 0.07, 0.15, 0.31, 0.63, 1.27]
    chain = straight_chain([(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (15.0, 0.0)])
    columns, node_rows, triangles, statistics = band_of(chain, [7, 4, 4], heights)
    @test node_rows == [7, 4, 4, 4]
    @test length(columns) == 4 && length(columns[1].nodes) == 7
    @test near(columns[1].nodes[7], (0.0, 1.27))
    @test near(columns[2].nodes[end], (5.0, 0.15))
    # 4 quads between every pair (12), the three excess rows of the first column collapsed.
    @test statistics.quads == 12 && statistics.collapse_triangles == 3
    @test length(triangles) == 2 * 12 + 3
    @test all(t -> PWM.orient(t...) > 0.0, triangles)
    @test statistics.quad_min_abs_sin ≈ 1.0
    # The band top is one edge per column pair; the step is inside the band.
    tops = [c.nodes[end] for c in columns]
    @test all(near.(tops, [(0.0, 1.27), (5.0, 0.15), (10.0, 0.15), (15.0, 0.15)]))
    # A step in the other direction collapses onto the shorter (first) column.
    _, _, _, reversed = band_of(chain, [4, 4, 7], heights)
    @test reversed.collapse_triangles == 3 && reversed.quads == 12
end

@testset "outward corners: scaled bisector and fans" begin
    heights = PWM.band_heights(0.05, 2.0, 5)
    # The band around a 10 x 10 metal square (the gap partition walks the hole clockwise).
    bases = PWM.Point2[]
    for (a, b) in (
        ((0.0, 0.0), (0.0, 10.0)),
        ((0.0, 10.0), (10.0, 10.0)),
        ((10.0, 10.0), (10.0, 0.0)),
        ((10.0, 0.0), (0.0, 0.0))
    )
        push!(bases, a, ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2))
    end
    chain = PWM.MetalChain(
        true,
        bases,
        [1, 1, 2, 2, 3, 3, 4, 4],
        [isodd(i) for i = 1:8],
        nothing,
        nothing,
        0,
        0
    )
    columns, node_rows, triangles, statistics = band_of(chain, fill(5, 8), heights)
    # Default threshold 90: the right-angle corners are scaled bisectors, exactly h_k from
    # both edge lines (the recorded corner treatment), no fan.
    @test statistics.fans == 0 && statistics.mitre_corners == 4
    @test length(columns) == 8
    @test near(columns[1].nodes[end], (-1.55, -1.55))
    @test statistics.quad_min_abs_sin ≈ sin(pi / 4) atol = 1.0e-9
    @test all(t -> PWM.orient(t...) > 0.0, triangles)
    @test all(node_rows .== 5)
    # Threshold 45: fans of 7 columns at every corner, row 1 triangles, exact radii.
    columns, _, triangles, fanned =
        band_of(chain, fill(5, 8), heights; fan_turn=deg2rad(45.0))
    @test fanned.fans == 4 && fanned.fan_triangles == 4 * 6
    @test length(columns) == 4 * 7 + 4
    @test all(hypot(c.nodes[end]...) ≈ 1.55 for c in columns[1:7])
    @test near(columns[1].nodes[end], (0.0, -1.55))
    @test near(columns[7].nodes[end], (-1.55, 0.0))
    @test fanned.quad_min_abs_sin > 0.97
    @test all(t -> PWM.orient(t...) > 0.0, triangles)
    @test PWM.check_band_collisions(triangles, 5.0, "test") > 0
    # The length cap of the scaled column (review M1): a corner whose adjacent segments face
    # a front at d keeps only the rows with h_k sqrt 2 <= 0.4 d. At d = 2.0 (r 0.05) the
    # segment rule allows 4 rows (h_4 = 0.75) but the sqrt 2-scaled column only 3 (0.495);
    # the straight columns keep their 4 rows and the step collapses.
    d = 2.0
    @test PWM.band_cap_rows(d, 0.05, 5) == 4
    columns, node_rows, triangles, capped =
        band_of(chain, fill(4, 8), heights; distances=fill(d, 8))
    @test node_rows == [3, 4, 3, 4, 3, 4, 3, 4]
    @test capped.scale_capped_columns == 4 && capped.scale_clamped_columns == 0
    @test capped.max_mitre_scale ≈ sqrt(2.0) && capped.max_inward_scale == 1.0
    @test hypot(columns[1].nodes[end]...) ≈ 0.35 * sqrt(2.0)
    @test hypot(columns[1].nodes[end]...) <= 0.4 * d
    @test capped.collapse_triangles == 8 && all(t -> PWM.orient(t...) > 0.0, triangles)
    # Fans have scale 1 and are not capped by the length rule.
    _, fan_rows, _, fan_capped =
        band_of(chain, fill(4, 8), heights; fan_turn=deg2rad(45.0), distances=fill(d, 8))
    @test all(fan_rows .== 4) && fan_capped.scale_capped_columns == 0
end

@testset "fan threshold tolerance" begin
    # Review M2: a right angle whose 1e-9-rounded coordinates turn by 90.0000001 deg stays a
    # mitre (tolerance 1e-6 deg); turns beyond the tolerance fan.
    heights = PWM.band_heights(0.05, 2.0, 5)
    function corner_chain(turn_deg)
        # Edge 1 along +x (partition above), edge 2 turning right by turn_deg down to y = 0.
        phi = deg2rad(turn_deg)
        t2 = (cos(-phi), sin(-phi))
        bases = [(0.0, 10.0), (5.0, 10.0), (10.0, 10.0)]
        push!(bases, (10.0 + 5.0 * t2[1], 10.0 + 5.0 * t2[2]))
        push!(bases, (10.0 + 10.0 * t2[1], 10.0 + 10.0 * t2[2]))
        return PWM.MetalChain(
            false,
            bases,
            [1, 1, 2, 2],
            [false, false, true, false, false],
            (0.0, 1.0),
            (1.0, 0.0),
            1,
            2
        )
    end
    @test PWM.FAN_TURN_TOLERANCE_DEG == 1.0e-6
    for turn_deg in (90.0, 90.0000001, 90.0 + 0.9e-6)
        _, _, _, statistics = band_of(corner_chain(turn_deg), fill(5, 4), heights)
        @test statistics.fans == 0 && statistics.mitre_corners == 1
    end
    for turn_deg in (90.0 + 1.1e-6, 90.00001, 90.01)
        _, _, _, statistics = band_of(corner_chain(turn_deg), fill(5, 4), heights)
        @test statistics.fans == 1 && statistics.mitre_corners == 0
    end
    # The threshold itself still moves with the option.
    _, _, _, lowered =
        band_of(corner_chain(90.0), fill(5, 4), heights; fan_turn=deg2rad(89.0))
    @test lowered.fans == 1
end

@testset "inward corners, T junctions and wall ends" begin
    heights = PWM.band_heights(0.01, 2.0, 7)
    # A 90-degree left turn with 0.5-um segments: the clearance rule h_k <= 0.5 s tan(45)
    # = 0.25 keeps 4 rows at the corner and its neighbours, 7 rows farther away.
    bases = [
        (0.0, 0.0),
        (5.0, 0.0),
        (10.0, 0.0),
        (10.5, 0.0),
        (10.5, 0.5),
        (10.5, 5.5),
        (10.5, 10.5)
    ]
    chain = PWM.MetalChain(
        false,
        bases,
        [1, 1, 1, 2, 2, 2],
        [false, false, false, true, false, false, false],
        (0.0, 1.0),
        (-1.0, 0.0),
        1,
        2
    )
    columns, node_rows, triangles, statistics = band_of(chain, fill(7, 6), heights)
    @test statistics.inward_corners == 1
    @test node_rows[4] == 4 && node_rows[3] == 4 && node_rows[5] == 4
    @test node_rows[2] == 7 && node_rows[6] == 7
    @test statistics.clearance_capped_columns == 3 &&
          statistics.clearance_clamped_columns == 0
    # The corner column is the bisector scaled by sqrt 2: exactly h_k from both edges.
    @test near(columns[4].nodes[end], (10.5 - 0.15, 0.15))
    @test all(t -> PWM.orient(t...) > 0.0, triangles)
    @test PWM.check_band_collisions(triangles, 5.0, "test") > 0
    # A straight-through T node (two collinear curves meeting): the exact normal, no cap.
    straight = PWM.MetalChain(
        false,
        [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0)],
        [1, 2],
        [false, true, false],
        (0.0, 1.0),
        (0.0, 1.0),
        1,
        2
    )
    columns, node_rows, _, straight_statistics = band_of(straight, [7, 7], heights)
    @test node_rows == [7, 7, 7] && straight_statistics.inward_corners == 0
    @test near(columns[2].nodes[end], (5.0, 1.27))
    # An oblique wall end: the end column runs along the wall, scaled by 1 / sin(theta).
    # (walls x = 0 and x = 16, the partition above the edge: the end columns run up the walls)
    oblique = PWM.MetalChain(
        false,
        [(0.0, 0.0), (8.0, 6.0), (16.0, 12.0)],
        [1, 1],
        [false, false, false],
        (0.0, 1.0),
        (0.0, 1.0),
        1,
        2
    )
    columns, node_rows, triangles, oblique_statistics = band_of(oblique, [7, 7], heights)
    @test oblique_statistics.wall_end_columns == 2 && node_rows == [7, 7, 7]
    @test near(columns[1].nodes[end], (0.0, 1.27 / 0.8))
    @test near(columns[3].nodes[end], (16.0, 12.0 + 1.27 / 0.8))
    @test near(columns[2].nodes[end], (8.0 - 1.27 * 0.6, 6.0 + 1.27 * 0.8))
    @test all(t -> PWM.orient(t...) > 0.0, triangles)
    @test oblique_statistics.max_wall_end_scale ≈ 1.25
    # The wall-end column (scale 1.25) facing a front at 2.0 um keeps h_k x 1.25 <= 0.8: 6
    # rows (0.7875) instead of 7; the straight column keeps 7 (review M1).
    _, capped_rows, _, capped = band_of(oblique, [7, 7], heights; distances=[2.0, 2.0])
    @test capped_rows == [6, 7, 6] && capped.scale_capped_columns == 2
    # The inward corner above with a front at 0.3 um on the band side: the sqrt 2-scaled
    # corner column keeps h_k sqrt 2 <= 0.12 (3 rows, 0.099), tighter than the clearance's 4.
    _, inward_rows, _, inward_capped =
        band_of(chain, fill(7, 6), heights; distances=fill(0.3, 6))
    @test inward_rows[4] == 3 && inward_rows[3] == 4
    @test inward_capped.max_inward_scale ≈ sqrt(2.0)
    # The mitre threshold cannot exceed 120 degrees (the scaled column would exceed 2 bands).
    @test_throws ErrorException mesh_polygon_window(
        read_polygon_set(single_plane_set("x", [[0, 0], [20, 0], [20, 10], [0, 10]])),
        0.05,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true,
        fan_turn_angle_deg=130.0
    )
end

@testset "own band on slot windows" begin
    # The 1-um strip Gmsh could not mesh at r10: 7 rows outside, 5 rows (0.31 um) inside.
    strip = read_polygon_set(
        single_plane_set(
            "strip-1um",
            [
                [0, 0],
                [9.5, 0],
                [9.5, 70],
                [10.5, 70],
                [10.5, 0],
                [20, 0],
                [20, 80],
                [0, 80]
            ]
        )
    )
    manifest = mesh_polygon_window(
        strip,
        0.01,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
    band = manifest["band"]
    @test band["mode"] == "own" &&
          manifest["radial_band_thickness_mode"] == "own_structured_band"
    @test band["per_side_segment_length_um_by_rows"] == Dict("7" => 142.0, "5" => 140.0)
    @test band["wall_end_columns"] == 4 && band["fans"] == 0 && band["inward_corners"] == 2
    @test band["segments_below_2p5r"] == 0 && band["clearance_clamped_columns"] == 0
    @test manifest["plan_boundary_layer_end_points"] == 4
    @test manifest["first_layer_normal_height_um"]["minimum"] ≈ 0.01 rtol = 1.0e-9
    @test manifest["first_layer_normal_height_um"]["maximum"] ≈ 0.01 rtol = 1.0e-9
    @test occursin("per straight run", manifest["band_cap"]["rule_statistic"])
    # The same window with the Gmsh band (the recorded path) still runs at r50.
    gmsh_manifest = mesh_polygon_window(
        strip,
        0.05,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true,
        band_mode=:gmsh
    )
    @test gmsh_manifest["band"]["mode"] == "gmsh"
    @test gmsh_manifest["radial_band_thickness_mode"] == "geometric_sum_x_1p000001"
    @test gmsh_manifest["plan_perimeter_edges"] == manifest["plan_perimeter_edges"]
    # A 0.15-um slot at r 0.1: the first rows (never dropped) collide -> refused with the
    # location, not meshed.
    collision = read_polygon_set(
        single_plane_set(
            "slot-collision",
            [
                [0, 0],
                [9.5, 0],
                [9.5, 70],
                [9.65, 70],
                [9.65, 0],
                [20, 0],
                [20, 80],
                [0, 80]
            ]
        )
    )
    message = try
        mesh_polygon_window(
            collision,
            0.1,
            5.0,
            tempname() * ".msh2";
            verbose=false,
            plan_only=true
        )
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("collide", message) ||
          occursin("collision", message) ||
          occursin("no region", message)
    @test occursin("partition", message)
    # An oblique wall end and a chamfer: exact plan areas, first-layer heights exactly r.
    chamfer = read_polygon_set(
        single_plane_set(
            "wall-chamfer",
            [[0, 0], [40, 0], [40, 10], [16.25, 10], [7.123, 10.712], [0, 12.347]];
            box_x=[0.0, 40.0],
            box_y=[0.0, 30.0]
        )
    )
    chamfer_manifest = mesh_polygon_window(
        chamfer,
        0.05,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
    @test chamfer_manifest["band"]["wall_end_columns"] == 4
    @test chamfer_manifest["first_layer_normal_height_um"]["minimum"] ≈ 0.05 rtol = 1.0e-9
    @test chamfer_manifest["first_layer_normal_height_um"]["maximum"] ≈ 0.05 rtol = 1.0e-9
end

# Two convex metal corners facing each other diagonally at corner-to-corner distance D (a
# ground square and an island square in a 20 x 20 box): the adjacent segments see each other
# at exactly D, so the segment rule alone let the two sqrt 2-scaled mitre tips reach 1.13 D
# (review M1).
function diagonal_corner_set(D)
    e = D / sqrt(2.0)
    return Dict(
        "Version" => 1,
        "Name" => "diagonal-corners",
        "Box" => Dict("X" => [0.0, 20.0], "Y" => [0.0, 20.0]),
        "Planes" => [
            Dict(
                "Name" => "L1",
                "SurfaceZ" => 0.0,
                "Facing" => "up",
                "SubstrateThickness" => 20.0,
                "Polygons" => [
                    Dict(
                        "Conductor" => "ground",
                        "Outer" => [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]]
                    ),
                    Dict(
                        "Conductor" => "island",
                        "Outer" => [
                            [10.0 + e, 10.0 + e],
                            [20.0, 10.0 + e],
                            [20.0, 20.0],
                            [10.0 + e, 20.0]
                        ]
                    )
                ]
            )
        ],
        "Vacuum" => Dict("Below" => 0.0, "Above" => 20.0)
    )
end

@testset "diagonally facing convex corners (M1)" begin
    # The reviewer's refusal windows D in [2.5 h_k, 2.83 h_k): r50 0.9 / 2.0 um, r10 1.6 / 3.2
    # um were refused by the collision backstop before the length cap; they mesh now, with the
    # two corner columns shortened and every other column as before.
    for (radial_um, D) in ((0.05, 0.9), (0.05, 2.0), (0.01, 1.6), (0.01, 3.2))
        manifest = mesh_polygon_window(
            read_polygon_set(diagonal_corner_set(D)),
            radial_um,
            5.0,
            tempname() * ".msh2";
            verbose=false,
            plan_only=true
        )
        band = manifest["band"]
        @test band["scale_capped_columns"] == 2 && band["scale_clamped_columns"] == 0
        @test band["mitre_outward_corners"] == 2 && band["fans"] == 0
        @test band["max_mitre_scale"] ≈ sqrt(2.0)
        @test manifest["first_layer_normal_height_um"]["minimum"] ≈ radial_um rtol = 1.0e-9
        @test manifest["first_layer_normal_height_um"]["maximum"] ≈ radial_um rtol = 1.0e-9
    end
    # Outside the windows (D = 1.5 um at r50: 3 rows, 0.35 sqrt 2 = 0.495 <= 0.6) nothing is
    # capped.
    manifest = mesh_polygon_window(
        read_polygon_set(diagonal_corner_set(1.5)),
        0.05,
        5.0,
        tempname() * ".msh2";
        verbose=false,
        plan_only=true
    )
    @test manifest["band"]["scale_capped_columns"] == 0
end

@testset "wall-touching wedge (writer constraint)" begin
    # SCHEMA: a metal vertex on the box wall without an edge along the wall makes the gap
    # partition touch itself at that vertex (four boundary curves): refused with the location.
    wedge =
        read_polygon_set(single_plane_set("wedge", [[5.0, 0.0], [9.0, 4.0], [1.0, 4.0]]))
    message = try
        mesh_polygon_window(
            wedge,
            0.05,
            5.0,
            tempname() * ".msh2";
            verbose=false,
            plan_only=true
        )
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("touches itself", message)
end
