# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# R2 "AcuteCornerRetry" (supervisor decisions 636 (3) / 638 / 639 / 641; round4/corner-surface
# DESIGN.md section 11): the RETRY-ON-FAIL wrapper of generate_spatial_coupon, its trigger
# predicate, the (R-a) corner-line remap judged on the FINAL configuration, the dihedral floor
# f(phi), the fail-closed texts, and - on the production-size synthetic tips of
# test_face_end_tubes.jl (write_tip_inputs) - the retry on the 18.75-degree tip that fails today,
# the identity by construction of a passing tip, and the geometric fail-closed of a 15-degree tip.

const TIP_INPUTS = read(joinpath(@__DIR__, "test_face_end_tubes.jl"), String)
let source = TIP_INPUTS
    start = findfirst("function write_tip_inputs(", source)
    stop = findfirst("function build_tip_coupon(", source)
    include_string(Main, source[start[1]:(stop[1] - 1)])
end

const FIXTURE_OPENING = deg2rad(18.75)
const FIXTURE_OVERETCH = 0.05
const FIXTURE_CORNER_SIZE = 0.00025
const FIXTURE_FLOOR = SEED_COLLAPSE_SIZE_FRACTION * FIXTURE_CORNER_SIZE   # 0.125 nm
const FIXTURE_TOLERANCE = 1.0e-9

fixture_sides(phi=FIXTURE_OPENING) =
    (walls=[[cos(phi / 2), sin(phi / 2)], [cos(phi / 2), -sin(phi / 2)]], metal=[1.0, 0.0])

# A synthetic acute-tip star (um): the corner O, the corner line below (c1 .. c3) and above
# (t), the z = 0 ring a0 / b0 on the lines A / B (r 0.233 nm < r_min 0.384 nm), the walls'
# first ring at depth (9 on wall A at r 0.15, 10 on wall B at r 0.12), the second ring at r
# 0.5 nm (outside r_min), a z = 0 face node and one metal-sidewall node above z = 0 per wall.
function fixture_star(; sidewall_between=false, degenerate=false, duplicate=false)
    phi = FIXTURE_OPENING
    wall_1 = [cos(phi / 2), sin(phi / 2), 0.0]
    wall_2 = [cos(phi / 2), -sin(phi / 2), 0.0]
    nm = 1.0e-3
    at(direction, r, z) = r * nm .* direction .+ [0.0, 0.0, z * nm]
    points = [
        [0.0, 0.0, 0.0],                 # 1  O (CAD point)
        [0.0, 0.0, -0.24nm],             # 2  c1
        [0.0, 0.0, -0.52nm],             # 3  c2
        [0.0, 0.0, -0.9nm],              # 4  c3
        at(wall_1, 0.233, 0.0),          # 5  a0 (line A)
        at(wall_2, 0.233, 0.0),          # 6  b0 (line B)
        at(wall_1, 0.5, 0.0),            # 7  a1 (line A)
        at(wall_2, 0.5, 0.0),            # 8  b1 (line B)
        at(wall_1, 0.15, -0.16),         # 9  wall A first ring -> c1
        at(wall_2, 0.12, -0.11),         # 10 wall B first ring -> O
        at(wall_1, 0.5, -0.45),          # 11 wall A second ring
        at(wall_2, 0.5, -0.45),          # 12 wall B second ring
        at(wall_1, 0.5, -0.9),           # 13
        at(wall_2, 0.5, -0.9),           # 14
        [0.6nm, 0.0, 0.0],               # 15 z = 0 face interior
        at(wall_1, 0.3, 0.3),            # 16 metal sidewall A
        at(wall_2, 0.3, 0.3),            # 17 metal sidewall B
        [0.0, 0.0, 0.3nm],               # 18 corner line above z = 0
        at(wall_1, 0.15, 0.05)           # 19 a sidewall node between a0 and O (the 28-degree class)
    ]
    P = hcat(points...)
    fixed = falses(size(P, 2)); fixed[1] = true
    # Planar faces oriented to a fixed normal each (the orientation only has to be consistent).
    function oriented(triangle, normal)
        n = cross(P[:, triangle[2]] .- P[:, triangle[1]], P[:, triangle[3]] .- P[:, triangle[1]])
        return dot(n, normal) >= 0.0 ? triangle : (triangle[1], triangle[3], triangle[2])
    end
    normal_1 = cross(wall_1, [0.0, 0.0, 1.0]); normal_2 = cross(wall_2, [0.0, 0.0, 1.0]); up = [0.0, 0.0, 1.0]
    triangles = NTuple{3, Int}[]; entities = Int[]
    function face!(entity, normal, list)
        for triangle in list
            push!(triangles, oriented(triangle, normal)); push!(entities, entity)
        end
    end
    face!(201, normal_1, [(1, 5, 9), (5, 7, 9), (7, 11, 9), (1, 9, 2), (2, 9, 11), (2, 11, 3), (3, 11, 13), (3, 13, 4)])
    face!(202, normal_2, [(1, 6, 10), (6, 8, 10), (8, 12, 10), (1, 10, 2), (2, 10, 12), (2, 12, 3), (3, 12, 14), (3, 14, 4)])
    if duplicate
        face!(203, up, [(1, 6, 5), (5, 7, 15), (6, 7, 15), (6, 8, 15)])
    else
        face!(203, up, [(1, 6, 5), (5, 6, 15), (5, 15, 7), (6, 8, 15)])
    end
    if sidewall_between
        face!(301, normal_1, [(1, 5, 19), (5, 16, 19), (19, 16, 18), (1, 19, 18), (5, 7, 16)])
    else
        face!(301, normal_1, [(1, 5, 16), (5, 7, 16), (1, 16, 18)])
    end
    face!(302, normal_2, [(1, 6, 17), (6, 8, 17), (1, 17, 18)])
    degenerate && face!(201, normal_1, [(9, 3, 4)])
    lines = [(1, 5), (5, 7), (1, 6), (6, 8), (1, 2), (2, 3), (3, 4), (1, 18)]
    line_entities = [101, 101, 102, 102, 103, 103, 103, 104]
    surface = falses(size(P, 2))
    for triangle in triangles, i in triangle
        surface[i] = true
    end
    for line in lines, i in line
        surface[i] = true
    end
    candidate_allowed = [surface[i] && !fixed[i] for i in axes(P, 2)]
    return (points=P, triangles=triangles, entities=entities, lines=lines,
            line_entities=line_entities, candidate_allowed=candidate_allowed)
end

remap_fixture(star; non_simplex=falses(size(star.points, 2)), sides=fixture_sides()) =
    acute_corner_ring_remap(star.points, (0.0, 0.0, 0.0), sides, FIXTURE_FLOOR, FIXTURE_CORNER_SIZE,
                            FIXTURE_OVERETCH, FIXTURE_TOLERANCE, star.candidate_allowed, non_simplex,
                            star.triangles, star.entities, star.lines, star.line_entities)

@testset "the dihedral floor f(phi) reproduces the design table" begin
    for (degrees, expected) in ((15.0, 5.371), (17.0, 4.731), (18.75, 4.283), (20.0, 4.010),
                                (22.5, 3.555), (28.0, 2.836), (45.0, 1.707), (90.0, 1.414))
        @test isapprox(kappa_reg_dihedral_floor(deg2rad(degrees)), expected; atol=1.0e-3)
    end
    # The acute-tip gate (open item R2-G): f crosses CornerShapeGate 5.0 between 16.0 and 16.2
    # degrees (5.031 / 4.968; the design text's rounded "16.3") and the target 3.8 between 21.0
    # and 21.2 degrees (3.815 / 3.778); the regular tetrahedron's dihedral reaches 1.
    @test kappa_reg_dihedral_floor(deg2rad(16.0)) > 5.0 > kappa_reg_dihedral_floor(deg2rad(16.2))
    @test kappa_reg_dihedral_floor(deg2rad(21.0)) > 3.8 > kappa_reg_dihedral_floor(deg2rad(21.2))
    @test isapprox(kappa_reg_dihedral_floor(deg2rad(16.5)), 4.877; atol=1.0e-3)
    @test isapprox(kappa_reg_dihedral_floor(acos(1 / 3)), 1.0; atol=1.0e-3)
    @test_throws ErrorException kappa_reg_dihedral_floor(0.0)
    @test_throws ErrorException kappa_reg_dihedral_floor(pi)
end

@testset "the locked-class trigger predicate (decision 639 (1) (b))" begin
    row(; kind="Invariant", phi=18.75, metal_sign=1.0, sides=true) = Dict{String, Any}(
        "Kind" => kind, "Point" => [0.0, 0.0, 0.0], "After" => 5.43, "Gate" => 5.0,
        "Sides" => sides ? Dict{String, Any}(
            "Walls" => [[cos(deg2rad(phi) / 2), sin(deg2rad(phi) / 2)],
                        [cos(deg2rad(phi) / 2), -sin(deg2rad(phi) / 2)]],
            "Metal" => [metal_sign, 0.0]) : nothing)
    @test acute_corner_locked_class(row(), true)
    @test acute_corner_locked_class(row(; phi=28.0), true)             # 2 sin(14 deg) = 0.484 < 0.5
    @test !acute_corner_locked_class(row(; phi=28.96), true)           # 2 sin(14.48 deg) = 0.5001
    @test !acute_corner_locked_class(row(; phi=52.816), true)          # 448693's sharpest corner
    @test !acute_corner_locked_class(row(), false)                     # the thin kind is out
    @test !acute_corner_locked_class(row(; kind="Legacy"), true)
    @test !acute_corner_locked_class(row(; sides=false), true)
    @test !acute_corner_locked_class(row(; metal_sign=-1.0), true)     # the reflex metal sector
    phi, convex = corner_metal_opening(fixture_sides(deg2rad(18.75)))
    @test isapprox(rad2deg(phi), 18.75; atol=1.0e-9) && convex
    # The wrapper's trigger: ONLY the CornerShapeGate failure, every row in the class, else nothing.
    only_corner = SeedGateFailure(["Seed invariant semantic-corner kappa_reg exceeds CornerShapeGate"], [row()])
    entries = acute_corner_retry_entries(only_corner, true)
    @test entries == [Dict{String, Any}("Point" => [0.0, 0.0, 0.0], "FirstAttemptAfter" => 5.43)]
    @test acute_corner_retry_entries(only_corner, false) === nothing
    @test acute_corner_retry_entries(SeedGateFailure(["corner", "Seed required region keeps 1 cells below"], [row()]), true) === nothing
    @test acute_corner_retry_entries(SeedGateFailure(["Seed required region keeps 1 cells below"], Dict{String, Any}[]), true) === nothing
    @test acute_corner_retry_entries(SeedGateFailure(["corner"], [row(), row(; phi=52.816)]), true) === nothing
    @test sprint(showerror, SeedGateFailure(["a", "b"], Dict{String, Any}[])) == "a; b"
    @test occursin("SIMULTANEOUSLY", ACUTE_CORNER_RETRY_RULE) && occursin("decision 363", ACUTE_CORNER_RETRY_RULE)
end

@testset "the (R-a) remap on the synthetic acute-tip star: candidates, corner-line targets, the final configuration" begin
    star = fixture_star()
    result = remap_fixture(star)
    @test isempty(result.findings)
    @test sort(result.candidates) == [5, 6, 9, 10]
    @test result.map == Dict(5 => 1, 6 => 1, 9 => 2, 10 => 1)        # z = 0 ring -> O; wall A node -> c1; wall B node -> O
    @test sort(result.corner_line) == [1, 2, 3, 4]                    # the corner line within [-Overetch, 0] only
    record = result.record
    @test isapprox(record["OpeningDegrees"], 18.75; atol=1.0e-9)
    @test isapprox(record["FirstRingRadius"], FIXTURE_FLOOR / (2 * sin(FIXTURE_OPENING / 2)); atol=1.0e-15)
    @test isapprox(record["FirstRingRadius"], 0.3837e-3; atol=1.0e-6)
    @test record["Floor"] == FIXTURE_FLOOR
    @test isapprox(record["CrossWedgeFirstRing"], 2 * FIXTURE_CORNER_SIZE * sin(FIXTURE_OPENING / 2); atol=1.0e-15)
    @test isapprox(record["DihedralFloor"], 4.283; atol=1.0e-3)
    @test record["TrianglesRemoved"] == 10 && record["TrianglesRewritten"] == 8
    @test record["LinesRemoved"] == 2 && record["LinesRewritten"] == 2
    @test record["InvalidRemaps"] == 0 && record["CornerLineNodes"] == 4
    @test length(record["Candidates"]) == 4 && record["Candidates"][1]["TargetZ"] == 0.0
    @test [c["Position"][3] for c in record["Candidates"]] == sort([c["Position"][3] for c in record["Candidates"]]; rev=true)
    kinds = Dict((row["Dimension"], row["Entity"]) => row["Kind"] for row in record["Entities"])
    @test kinds == Dict((1, 101) => "ZeroPlane", (1, 102) => "ZeroPlane", (2, 201) => "PedestalWedge",
                        (2, 202) => "PedestalWedge", (2, 203) => "ZeroPlane", (2, 301) => "MetalSidewall",
                        (2, 302) => "MetalSidewall")
    by_entity = Dict(row["Entity"] => row for row in record["Entities"])
    @test by_entity[201]["TrianglesRemoved"] == 3 && by_entity[201]["TrianglesRewritten"] == 2
    @test by_entity[202]["TrianglesRemoved"] == 3 && by_entity[202]["TrianglesRewritten"] == 2
    @test by_entity[203]["TrianglesRemoved"] == 2 && by_entity[203]["TrianglesRewritten"] == 2
    @test by_entity[301]["TrianglesRemoved"] == 1 && by_entity[301]["TrianglesRewritten"] == 1
    @test by_entity[101]["LinesRemoved"] == 1 && by_entity[101]["LinesRewritten"] == 1
    # Every rewritten triangle keeps its normal; the vanished ones are exactly those with a repeated vertex.
    @test count(new -> new === nothing, values(result.triangles)) == 10
    for (k, new) in result.triangles
        new === nothing && continue
        old = star.triangles[k]
        n_old = cross(star.points[:, old[2]] .- star.points[:, old[1]], star.points[:, old[3]] .- star.points[:, old[1]])
        n_new = cross(star.points[:, new[2]] .- star.points[:, new[1]], star.points[:, new[3]] .- star.points[:, new[1]])
        @test dot(n_old, n_new) > 0.0
    end
    @test result.lines == Dict(1 => nothing, 2 => (1, 7), 3 => nothing, 4 => (1, 8))
    # The floor guarantee on the cleaned wedge: the surviving second ring at 0.5 nm on both walls is
    # 2 x 0.5 sin(phi / 2) = 0.163 nm apart, above the 0.125-nm floor.
    cleaned = record["CleanedWedge"]
    @test isapprox(cleaned["ShortestCrossWedgePair"], 2 * 0.5e-3 * sin(FIXTURE_OPENING / 2); atol=1.0e-12)
    @test cleaned["ShortestCrossWedgePair"] >= FIXTURE_FLOOR && cleaned["Nodes"] == 7
    @test occursin("corner-line node nearest in z", record["Rule"]) || occursin("CORNER-LINE node nearest in z", record["Rule"])
end

@testset "the final-configuration validity: reversed / degenerate / duplicate / non-simplex fail closed" begin
    reversed = remap_fixture(fixture_star(; sidewall_between=true))
    @test length(reversed.findings) == 1 && occursin("reversed triangle", reversed.findings[1]) &&
          occursin("(entity 301)", reversed.findings[1])
    @test reversed.record["InvalidRemaps"] == 1
    degenerate = remap_fixture(fixture_star(; degenerate=true))
    @test any(f -> occursin("degenerate triangle", f) && occursin("[9, 3, 4] -> [2, 3, 4]", f), degenerate.findings)
    duplicate = remap_fixture(fixture_star(; duplicate=true))
    @test any(f -> occursin("duplicate triangles", f), duplicate.findings)
    star = fixture_star()
    non_simplex = falses(size(star.points, 2)); non_simplex[5] = true
    tube = remap_fixture(star; non_simplex=non_simplex)
    @test tube.findings == ["non-simplex element references the mapped node 5"]
    # Not the locked class / the reflex sector: the remap refuses.
    @test_throws ErrorException remap_fixture(star; sides=fixture_sides(deg2rad(60.0)))
    @test_throws ErrorException remap_fixture(star; sides=(walls=fixture_sides().walls, metal=[-1.0, 0.0]))
    # The second-failure text: the gate failure of record first, then the f(phi) statement.
    again = SeedGateFailure(["Seed invariant semantic-corner kappa_reg exceeds CornerShapeGate after optimization: 5.6 > 5.0"],
                            [Dict{String, Any}("Point" => [0.0, 0.0, 0.0], "After" => 5.6, "Gate" => 5.0)])
    entries = [Dict{String, Any}("Point" => [0.0, 0.0, 0.0], "FirstAttemptAfter" => 6.96,
                                 "OpeningDegrees" => 15.0, "DihedralFloor" => 5.371)]
    message = acute_corner_retry_failure_message(again, entries)
    @test startswith(message, "Seed invariant semantic-corner kappa_reg exceeds CornerShapeGate after optimization: 5.6 > 5.0; acute-corner retry (decision 639) failed again")
    @test occursin("first attempt 6.96, retry 5.6, dihedral floor f(15.0 deg) = 5.371 > CornerShapeGate 5.0: the gate is unattainable at this tip by geometry", message)
    entries[1]["DihedralFloor"] = 4.283; entries[1]["OpeningDegrees"] = 18.75
    @test !occursin("unattainable", acute_corner_retry_failure_message(again, entries))
end

# The production-size synthetic tip (test_invariant_corner_builds.jl's build_production_tip option
# set) through generate_spatial_coupon_once (today's path) or the retry wrapper.
function build_tip(directory, phi, theta; fabricated, once=false)
    inputs = write_tip_inputs(directory, phi, theta; radius=0.5, side_length=0.6, side_length_b=0.6)
    census = joinpath(directory, "census.json")
    mesh = joinpath(directory, "coupon.msh")
    edge = fabricated ? 0.00025 : 0.002
    build = once ? generate_spatial_coupon_once : generate_spatial_coupon
    build(; signature=inputs.signature, mask=inputs.mask, boundary=inputs.boundary, fabricated=fabricated,
          filename=mesh, radius=0.5, metal_thickness=0.1, overetch=0.05, sidewall_angle=90.0,
          top_rounding=0.0, trench_rounding=0.0, lc_fine=0.025, lc_tangent=0.05, lc_far=0.152,
          max_nodes=4_000_000, max_elements=4_000_000, mesh_order=1, semantic_contract=inputs.semantic,
          corner_isotropy_radius=0.1, corner_census=census, edge_size=edge, edge_growth_ratio=2.0,
          corner_size=edge, prism_tubes=true, far_growth=0.5, maximum_corner_aspect=4.0,
          minimum_scaled_jacobian=0.01, maximum_jacobian_condition=1000.0,
          quality_displacement_over_normal=0.75, corner_shape_gate=5.0)
    return read(mesh), parse_json(read(census, String))
end

# The census without its per-run leaves: the wall clocks (every `Timings` sub-record) and the
# input path of the semantic contract (a per-directory path).
strip_volatile(value::Dict) =
    Dict(k => strip_volatile(v) for (k, v) in value if !(k in ("Timings", "SemanticContract")))
strip_volatile(value::Vector) = [strip_volatile(v) for v in value]
strip_volatile(value) = value

function has_key(value::Dict, name)
    haskey(value, name) && return true
    return any(has_key(v, name) for v in values(value))
end
has_key(value::Vector, name) = any(has_key(v, name) for v in value)
has_key(value, name) = false

@testset "a passing tip (22.5 degrees, fabricated) is byte-identical through the wrapper and carries no retry key" begin
    mktempdir() do directory
        mesh_once, census_once = build_tip(mkpath(joinpath(directory, "once")), 22.5, 10.0; fabricated=true, once=true)
        mesh_wrapped, census_wrapped = build_tip(mkpath(joinpath(directory, "wrapped")), 22.5, 10.0; fabricated=true)
        @test mesh_once == mesh_wrapped
        @test strip_volatile(census_once) == strip_volatile(census_wrapped)
        @test !has_key(census_wrapped, "AcuteCornerRetry") && !has_key(census_wrapped, "AcuteCornerRetries")
        row = census_wrapped["SeedQualityOptimization"]["CornerMeasures"][1]
        @test row["Passed"] && row["After"] <= 4.3
    end
end

@testset "the 18.75-degree tip (theta 10) fails today and is rebuilt by the retry; the record is on the retried row only" begin
    mktempdir() do directory
        failure = try
            build_tip(mkpath(joinpath(directory, "once")), 18.75, 10.0; fabricated=true, once=true)
            nothing
        catch caught
            caught
        end
        gmsh.isInitialized() == 1 && gmsh.finalize()
        @test failure isa SeedGateFailure
        @test length(failure.failures) == 1 && length(failure.invariant_failed) == 1
        @test startswith(sprint(showerror, failure),
                         "Seed invariant semantic-corner kappa_reg exceeds CornerShapeGate after optimization: 5.44")
        first_after = failure.invariant_failed[1]["After"]
        @test 5.0 < first_after < 5.5
        @test acute_corner_locked_class(failure.invariant_failed[1], true)
        @test acute_corner_retry_entries(failure, false) === nothing
        _, census = build_tip(mkpath(joinpath(directory, "wrapped")), 18.75, 10.0; fabricated=true)
        optimization = census["SeedQualityOptimization"]
        @test optimization["AcuteCornerRetries"] == [[0.0, 0.0, 0.0]]
        @test optimization["AcuteCornerRetryRule"] == ACUTE_CORNER_RETRY_RULE
        row = optimization["CornerMeasures"][1]
        @test row["Passed"] && row["After"] <= 5.0 && row["Kind"] == "Invariant"
        retry = row["AcuteCornerRetry"]
        @test retry["FirstAttemptAfter"] == first_after
        @test isapprox(retry["OpeningDegrees"], 18.75; atol=1.0e-6)
        @test retry["InvalidRemaps"] == 0 && length(retry["Candidates"]) == 6
        @test isapprox(retry["FirstRingRadius"], 0.3837e-3; atol=1.0e-6) && retry["Floor"] == 0.000125
        @test isapprox(retry["DihedralFloor"], 4.283; atol=1.0e-3)
        @test retry["TrianglesRewritten"] > 0 && retry["TrianglesRemoved"] > 0 && retry["LinesRemoved"] == 2
        @test count(e -> e["Kind"] == "MetalSidewall", retry["Entities"]) == 2
        @test retry["CleanedWedge"]["ShortestCrossWedgePair"] >= 0.000125
        # The floor is satisfiable after the cleanup: the ball's shortest edge is at or above it.
        @test census["Corners"][1]["EdgeMinimum"] >= 0.000125
        @test census["Corners"][1]["Measure"]["Value"] == row["After"]
    end
end

@testset "a 15-degree tip is retried and fails closed again with the dihedral-floor statement (open item R2-G)" begin
    mktempdir() do directory
        failure = try
            build_tip(directory, 15.0, 10.0; fabricated=true)
            nothing
        catch caught
            caught
        end
        gmsh.isInitialized() == 1 && gmsh.finalize()
        @test failure isa ErrorException
        message = sprint(showerror, failure)
        @test startswith(message, "Seed invariant semantic-corner kappa_reg exceeds CornerShapeGate after optimization:")
        @test occursin(r"acute-corner retry \(decision 639\) failed again at \[0\.0, 0\.0, 0\.0\]: first attempt 7\.\d+, retry 5\.\d+", message)
        @test occursin(r"dihedral floor f\(1[45]\.\d+ deg\) = 5\.371", message)
        @test occursin("> CornerShapeGate 5.0: the gate is unattainable at this tip by geometry", message)
        @test occursin("\"InvalidRemaps\": 0", message) && occursin("AcuteCornerRetry records [ {", message)
    end
end
