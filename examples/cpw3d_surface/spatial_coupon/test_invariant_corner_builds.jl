# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# Mesher design round 2 F5-A (supervisor decisions 351 / 358 / 363 / 365) on the family-5
# production-size synthetic tip family (the V7-c / test_face_end_tubes write_tip_inputs
# geometry at the production corner-shell structure: R 0.5 box, 0.6-um sides, NormalSize
# 0.025, TangentialSize 0.05, FarSize 0.152, EdgeSize = CornerSize 0.25 nm fabricated / 2 nm
# thin): the invariant corner measure kappa_reg replaces the vertex-0 lottery at the
# non-perpendicular corners, the legacy 90-degree corner is bitwise untouched, and a
# contract predating the rule fails closed.

const TIP_INPUTS = read(joinpath(@__DIR__, "test_face_end_tubes.jl"), String)
# Only the fixture writer of test_face_end_tubes.jl (its testsets are not run here).
let source = TIP_INPUTS
    start = findfirst("function write_tip_inputs(", source)
    stop = findfirst("function build_tip_coupon(", source)
    include_string(Main, source[start[1]:(stop[1] - 1)])
end

function build_production_tip(
    directory,
    phi,
    theta;
    fabricated,
    corner_shape_gate=5.0,
    legacy_contract=false
)
    inputs = write_tip_inputs(
        directory,
        phi,
        theta;
        radius=0.5,
        side_length=0.6,
        side_length_b=0.6
    )
    if legacy_contract
        # The contract as derived before the rule: SemanticCorners only.
        open(inputs.semantic, "w") do io
            return write_json(
                io,
                Dict{String, Any}(
                    "Version" => 1,
                    "SemanticCorners" => [[c[1], c[2], c[3]] for c in inputs.corners]
                )
            )
        end
    end
    census = joinpath(directory, "census.json")
    mesh = joinpath(directory, "coupon.msh")
    edge = fabricated ? 0.00025 : 0.002
    options = corner_shape_gate > 0.0 ? (; corner_shape_gate) : (;)
    generate_spatial_coupon(;
        signature=inputs.signature,
        mask=inputs.mask,
        boundary=inputs.boundary,
        fabricated=fabricated,
        filename=mesh,
        radius=0.5,
        metal_thickness=0.1,
        overetch=0.05,
        sidewall_angle=90.0,
        top_rounding=0.0,
        trench_rounding=0.0,
        lc_fine=0.025,
        lc_tangent=0.05,
        lc_far=0.152,
        max_nodes=4_000_000,
        max_elements=4_000_000,
        mesh_order=1,
        semantic_contract=inputs.semantic,
        corner_isotropy_radius=0.1,
        corner_census=census,
        edge_size=edge,
        edge_growth_ratio=2.0,
        corner_size=edge,
        prism_tubes=true,
        far_growth=0.5,
        maximum_corner_aspect=4.0,
        minimum_scaled_jacobian=0.01,
        maximum_jacobian_condition=1000.0,
        quality_displacement_over_normal=0.75,
        options...
    )
    return parse_json(read(census, String)), inputs
end

corner_rows(census) = census["SeedQualityOptimization"]["CornerMeasures"]

@testset "22.5-degree tips (O1 / a260b5b3d9d0 class): kappa_reg after-values, legacy gate would fail" begin
    mktempdir() do directory
        for (fabricated, bound) in ((false, 3.5), (true, 4.3))
            census, inputs = build_production_tip(
                joinpath(mkpath(joinpath(directory, fabricated ? "fab" : "thin"))),
                22.5,
                10.0;
                fabricated=fabricated
            )
            @test length(inputs.corners) == 1 && inputs.invariant == [[0.0, 0.0, 0.0]]
            @test census["SemanticCornerKinds"] == ["Invariant"]
            @test census["InvariantCorners"]["Points"] == [[0.0, 0.0, 0.0]] &&
                  census["InvariantCorners"]["CornerShapeGate"] == 5.0 &&
                  census["InvariantCorners"]["Target"] == 3.8
            row = corner_rows(census)[1]
            @test row["Kind"] == "Invariant" && row["Measure"] == "RegularCondition"
            @test row["Target"] == 3.8 && row["Gate"] == 5.0 && row["Passed"]
            # B11's after-value clauses on the synthetic twins of the production tips.
            @test row["After"] <= bound
            @test row["Before"] > row["After"]
            # The vertex-0 measure of the same optimized cells reads above the legacy gate
            # 4.0 (the node-order lottery the invariant measure retires): the legacy rule
            # would have failed this corner.
            @test row["Information"]["KappaV0Max"] > 4.0 || !fabricated
            @test row["Information"]["KappaRegMax"] == row["After"]
            @test row["BridgingSlivers"]["AboveGateAfter"] == 0
            measure = census["Corners"][1]["Measure"]
            @test measure["Kind"] == "Invariant" && measure["Name"] == "RegularCondition"
            @test isapprox(measure["Value"], row["After"]; atol=1.0e-9) &&
                  measure["Gate"] == 5.0
            @test census["SeedQualityOptimization"]["CornerAspectsAfter"] == [row["After"]]
        end
    end
end

@testset "a contract without the InvariantCorners record, or no gate, fails closed" begin
    mktempdir() do directory
        @test_throws ErrorException build_production_tip(
            directory,
            22.5,
            10.0;
            fabricated=false,
            legacy_contract=true
        )
        @test_throws ErrorException build_production_tip(
            directory,
            22.5,
            10.0;
            fabricated=false,
            corner_shape_gate=0.0
        )
    end
end

@testset "135-degree fabricated kink (the a260b5b3d9d0 concave 225) passes; the 150 / 151 flat-kink slivers are reconnected (F5-B)" begin
    mktempdir() do directory
        census, _ = build_production_tip(
            joinpath(mkpath(joinpath(directory, "k135"))),
            135.0,
            10.0;
            fabricated=true
        )
        row = corner_rows(census)[1]
        @test row["Kind"] == "Invariant" && row["Passed"] && row["After"] <= 3.5
        @test isempty(row["Reconnections"])
        # The flat-kink sliver (family-5 root cause B): kappa_reg stays in the identified band
        # [5.2, 5.6] under the measure alone (AfterDescent), above the cap 5.0, as a
        # bridging-sliver candidate above the gate - the reconnection trigger (decision 365);
        # the edge removal of its bridging edge (the 3-2 flip) and the second descent bring
        # the corner under 0.95 x CornerShapeGate (the design's 6 / 6 selection bar) with no
        # candidate left (B12).
        for (phi, theta) in ((150.0, 10.0), (151.0, 20.0))
            census, _ = build_production_tip(
                joinpath(mkpath(joinpath(directory, "k$(phi)-$(theta)"))),
                phi,
                theta;
                fabricated=true
            )
            row = corner_rows(census)[1]
            @test row["Kind"] == "Invariant" && row["Passed"]
            @test 5.2 <= row["AfterDescent"] <= 5.6
            @test row["After"] <= 0.95 * 5.0
            @test length(row["Reconnections"]) == 1
            reconnection = row["Reconnections"][1]
            @test reconnection["Kind"] == "edge-removal"
            @test reconnection["After"] < reconnection["Before"]
            # Exact records (decision 392 MINOR-7): the one trigger cell is the corner's worst
            # cell after the descent, so Before is AfterDescent to the digit, the removed
            # bridging edge is a 2-vertex element and the ring of 3 cells became 2.
            @test reconnection["Before"] == row["AfterDescent"] &&
                  reconnection["Round"] == 1
            @test length(reconnection["Element"]) == 2 &&
                  reconnection["Element"][1] != reconnection["Element"][2]
            @test reconnection["ReplacedCells"] == 3 && reconnection["AddedCells"] == 2
            @test row["BridgingSlivers"]["Before"] >= 1 &&
                  row["BridgingSlivers"]["AboveGateAfter"] == 0
            optimization = census["SeedQualityOptimization"]
            @test optimization["CornerReconnections"] == 1
            @test optimization["ReconnectionReplacedCells"] == 3 &&
                  optimization["ReconnectionAddedCells"] == 2
            rule = optimization["CornerReconnectionRule"]
            @test occursin("2-3 flip", rule) && occursin("edge removal", rule)
            @test occursin("trigger cell's", rule) && occursin("reused being skipped", rule)
        end
    end
end

@testset "SEAM: the tip bisector above the minimum opening, the unrefined seams below it (decision 368)" begin
    mktempdir() do directory
        # phi_split at the gate 5.0: the split-fan floor f(phi / 2) = 5.0 / 1.2.
        minimum_opening = rad2deg(tip_bisector_minimum_opening(5.0))
        @test 31.0 < minimum_opening < 32.0
        @test sector_regular_floor(deg2rad(11.25)) ≈ 5.8619 atol = 1.0e-3
        @test sector_regular_floor(deg2rad(60.0)) ≈ 1.0 atol = 1.0e-12
        @test sector_regular_floor(deg2rad(minimum_opening / 2)) ≈
              5.0 / TIP_BISECTOR_DESCENT_EXCESS atol = 1.0e-9
        # 45-degree thin tip (O4 class): the bisector to the half-width station 0.065 um, seam
        # count 0, kappa_reg within the gate.
        census, _ = build_production_tip(
            joinpath(mkpath(joinpath(directory, "t45"))),
            45.0,
            10.0;
            fabricated=false
        )
        bisectors = census["TipBisectors"]
        @test bisectors["Count"] == 1 && bisectors["UnrefinedTips"] == 0
        curve = bisectors["Curves"][1]
        @test curve["OpeningDegrees"] ≈ 45.0 atol = 1.0e-6
        @test curve["StationRule"] == "HalfWidthReachesCornerLaw"
        @test isapprox(curve["Length"], 0.0654; atol=1.0e-4)
        @test curve["Nodes"] >= 5 && length(curve["EmbeddedIn"]) == 1
        @test census["ThinSheetSeams"]["Count"] == 0 &&
              census["ThinSheetSeams"]["UnrefinedCrackSeams"] === nothing
        @test corner_rows(census)[1]["Passed"] && corner_rows(census)[1]["After"] <= 5.0
        # 22.5-degree thin tip (O1 / a260b5b3d9d0 class): below the minimum opening, no
        # bisector; the fan's pinched seams are recorded as UnrefinedCrackSeams, all within the
        # tip's tube station; the kappa_reg after-value is the unsplit one (B11 <= 3.5).
        census, _ = build_production_tip(
            joinpath(mkpath(joinpath(directory, "t22"))),
            22.5,
            10.0;
            fabricated=false
        )
        @test census["TipBisectors"]["Count"] == 0 &&
              census["TipBisectors"]["UnrefinedTips"] == 1
        seams = census["ThinSheetSeams"]
        @test seams["Count"] >= 1 && seams["UnrefinedCrackSeams"]["Count"] == seams["Count"]
        @test seams["UnrefinedCrackSeams"]["RefineCrackElements"] === false
        @test seams["UnrefinedCrackSeams"]["Tips"][1]["OpeningDegrees"] ≈ 22.5 atol = 1.0e-6
        @test corner_rows(census)[1]["After"] <= 3.5
    end
end

@testset "legacy 90-degree corners: the vertex-0 path, Information recorded, kinds Legacy" begin
    mktempdir() do directory
        census, inputs = build_production_tip(
            directory,
            90.0,
            0.0;
            fabricated=true,
            corner_shape_gate=0.0
        )
        @test length(inputs.corners) == 2 && isempty(inputs.invariant)
        @test census["SemanticCornerKinds"] == ["Legacy", "Legacy"]
        @test isempty(census["InvariantCorners"]["Points"]) &&
              census["InvariantCorners"]["CornerShapeGate"] === nothing
        for row in corner_rows(census)
            @test row["Kind"] == "Legacy" && row["Measure"] == "VertexFrameCondition"
            @test row["Target"] == 3.8 && row["Gate"] == 4.0 && row["After"] <= 4.0
            @test row["Information"]["KappaRegMax"] > 0.0 &&
                  row["Information"]["EtaMin"] > 0.0
        end
        @test all(row["Measure"]["Kind"] == "Legacy" for row in census["Corners"])
    end
end
