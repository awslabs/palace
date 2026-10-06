# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# Mesher design round 2 F5-A (supervisor decisions 351 / 358 / 363): the order-invariant
# corner measure kappa_reg, the exact dot-product corner-kind predicate, the contract
# check, the bridging-sliver predicate and the gated required-region optimization.

const REGULAR = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.5, sqrt(3.0) / 2.0, 0.0],
    [0.5, sqrt(3.0) / 6.0, sqrt(2.0 / 3.0)]
]
const TRIRECTANGULAR = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]

permutations4() =
    [p for p in Iterators.product(1:4, 1:4, 1:4, 1:4) if length(unique(p)) == 4]

@testset "kappa_reg: regular 1, trirectangular 2, order-invariant, factor-2 bounds" begin
    @test tetrahedron_regular_condition(REGULAR) ≈ 1.0 atol = 1.0e-12
    @test tetrahedron_regular_condition(TRIRECTANGULAR) ≈ 2.0 atol = 1.0e-12
    @test tetrahedron_aspect(TRIRECTANGULAR) ≈ 1.0 atol = 1.0e-12
    @test tetrahedron_aspect(REGULAR) ≈ 2.0 atol = 1.0e-12
    # The O1 thin 22.5-degree tip cell type: the corner first, two ring vertices at 2.04 nm
    # on the sheet edges 22.5 degrees apart, the fourth vertex 2.6 nm above the sheet - the
    # vertex-0 measure reads by first vertex, kappa_reg does not.
    tip = [
        [0.0, 0.0, 0.0],
        [0.00204, 0.0, 0.0],
        [0.00204 * cosd(22.5), 0.00204 * sind(22.5), 0.0],
        [0.0012, 0.0004, 0.0026]
    ]
    frames = Set{Float64}()
    regular = Set{Float64}()
    for p in permutations4()
        cell = tip[collect(p)]
        push!(frames, round(tetrahedron_aspect(cell); digits=6))
        push!(regular, round(tetrahedron_regular_condition(cell); digits=9))
    end
    @test length(frames) > 1
    @test length(regular) == 1
    @test tetrahedron_aspect(tip) >= cotd(11.25) * (1.0 - 1.0e-9)
    low, high = tetrahedron_frame_conditions(tip)
    @test low <= tetrahedron_aspect(tip) <= high
    kappa_reg = tetrahedron_regular_condition(tip)
    for p in permutations4()
        v0 = tetrahedron_aspect(tip[collect(p)])
        @test v0 / 2.0 <= kappa_reg * (1.0 + 1.0e-12) &&
              kappa_reg <= 2.0 * v0 * (1.0 + 1.0e-12)
    end
    @test tetrahedron_mean_ratio(REGULAR) ≈ 1.0 atol = 1.0e-12
    @test 0.8 < tetrahedron_mean_ratio(TRIRECTANGULAR) < 0.9
    @test tetrahedron_mean_ratio(tip) < 0.7
    @test INVARIANT_CORNER_TARGET == 3.8 == 0.95 * 4.0
end

# A plan-view loop as read_boundary returns it.
loop(points; plane=0.0, classes=fill("Physical", length(points))) =
    (conductor=1, plane=plane, hole=false, points=points, classes=classes)

# The plan-view quantum of the production coupons (Radius 1.9 um): 1.9e-9.
const TEST_QUANTUM = PLAN_VIEW_QUANTUM_OVER_RADIUS * 1.9

@testset "corner kinds: the exact dot-product predicate" begin
    # A rectilinear polygon with quantised (many-digit) coordinates: every corner legacy.
    square = loop([
        (2.0515421000000003, 0.6218871),
        (2.0515421000000003, -1.0296138000000001),
        (7.2265417, -1.0296138000000001),
        (7.2265417, 0.6218871)
    ])
    corners = [(p[1], p[2], 0.0) for p in square.points]
    kinds, sides, dots = semantic_corner_kinds(corners, [square], 1.0e-8, TEST_QUANTUM)
    @test all(kind -> kind === :legacy, kinds) && all(==(0), dots)
    @test all(norm(wall) ≈ 1.0 for side in sides for wall in side.walls)
    # The metal direction points into the (counter-clockwise) square: at the first vertex
    # (top left, walls down and right) it is the diagonal (1, -1) / sqrt 2.
    @test sides[1].metal ≈ [1.0, -1.0] ./ sqrt(2.0)
    @test in_metal_sector(sides[1], [1.0, -0.5]) && !in_metal_sector(sides[1], [-1.0, 1.0])
    # A 22.5-degree tip and a 150-degree kink are invariant; the box vertex where a
    # perpendicular metal side meets the box side (theta 0, decision 320) is legacy.
    tip = loop(
        [(0.0, 0.0), (0.6, 0.0), (0.6, 0.1), (0.6 * cosd(22.5), 0.6 * sind(22.5))];
        classes=["Physical", "Continuation", "Physical", "Physical"]
    )
    kinds, sides, dots = semantic_corner_kinds(
        [(0.0, 0.0, 0.0), (0.6, 0.0, 0.0)],
        [tip],
        1.0e-8,
        TEST_QUANTUM
    )
    @test kinds == [:invariant, :legacy]
    @test dots[1] > 0 && dots[2] == 0
    @test sort(sides[1].walls; by=first) ≈ [[cosd(22.5), sind(22.5)], [1.0, 0.0]]
    # The metal of the tip is the 22.5-degree wedge between its sides.
    @test sides[1].metal ≈ [cosd(11.25), sind(11.25)]
    @test in_metal_sector(sides[1], [cosd(10.0), sind(10.0)])
    @test !in_metal_sector(sides[1], [cosd(100.0), sind(100.0)]) &&
          !in_metal_sector(sides[1], [cosd(-30.0), sind(-30.0)])
    kink =
        loop([(0.0, 0.0), (0.6, 0.0), (0.6, 0.3), (0.6 * cosd(150.0), 0.6 * sind(150.0))])
    kinds, sides, dots =
        semantic_corner_kinds([(0.0, 0.0, 0.0)], [kink], 1.0e-8, TEST_QUANTUM)
    @test kinds == [:invariant] && dots[1] < 0
    @test sides[1].metal ≈ [cosd(75.0), sind(75.0)]
    # A hole loop (the metal outside): the metal direction flips to the 210-degree side.
    hole = (conductor=1, plane=0.0, hole=true, points=kink.points, classes=kink.classes)
    _, hole_sides, _ =
        semantic_corner_kinds([(0.0, 0.0, 0.0)], [hole], 1.0e-8, TEST_QUANTUM)
    @test hole_sides[1].metal ≈ -[cosd(75.0), sind(75.0)]
    # A corner absent from the loops, or present twice, fails closed; the plane matters.
    @test_throws ErrorException semantic_corner_kinds(
        [(5.0, 5.0, 0.0)],
        [square],
        1.0e-8,
        TEST_QUANTUM
    )
    @test_throws ErrorException semantic_corner_kinds(
        [(0.0, 0.0, 1.0)],
        [kink],
        1.0e-8,
        TEST_QUANTUM
    )
    @test_throws ErrorException semantic_corner_kinds(
        [(0.0, 0.0, 0.0)],
        [kink, kink],
        1.0e-8,
        TEST_QUANTUM
    )
end

# The seven semantic corners of the SCT loop end 1b26671c9080 (sct002-S1p, Radius 1.9;
# merge-preparation registration PBS 57005) with their two plan-view loop neighbours, as
# the plan-view boundary CSV spells them: the two 135-degree corners, the box corner, a
# rigidly ROTATED perpendicular corner (both sides tilted 9.5e-7 rad: (2.0000008, -0.0000019)
# and (0.0000019, 2.0000008), exactly 1000 quanta off the axes) and three corners whose
# sides are tilted by real ~1e-6-rad amounts (supervisor decision 416's erratum to 410).
const LOOP_END_CORNERS = [
    (
        before=(10.8001909, -14.6336328),
        point=(14.6001928, -10.833638500000001),
        after=(14.600213700000001, 11.1663608),
        kind=:invariant
    ),
    (
        before=(14.6001928, -10.833638500000001),
        point=(14.600213700000001, 11.1663608),
        after=(10.8002023, 14.966381700000001),
        kind=:invariant
    ),
    (
        before=(-19.730892, 14.966381700000001),
        point=(-19.730892, 1.1663834),
        after=(-6.8997949, 1.1663834),
        kind=:legacy
    ),
    (
        before=(-6.899793, 3.1663823),
        point=(-8.8997938, 3.1663842),
        after=(-8.8997919, 5.166385),
        kind=:legacy
    ),
    (
        before=(-8.8997938, 3.1663842),
        point=(-8.8997919, 5.166385),
        after=(-7.8997915, 5.1663831),
        kind=:invariant
    ),
    (
        before=(-8.8997919, 5.166385),
        point=(-7.8997915, 5.1663831),
        after=(-7.8997915, 6.1663835),
        kind=:invariant
    ),
    (
        before=(-7.8998029, -5.8336156),
        point=(-7.899801, -4.8336171000000006),
        after=(-19.730892, -4.8336171000000006),
        kind=:invariant
    )
]

@testset "corner kinds: exact on the quantum counts, rotation-independent (decision 416)" begin
    # Every coordinate of the loop end is an integer count of its 1.9e-9 quantum.
    for corner in LOOP_END_CORNERS,
        point in (corner.before, corner.point, corner.after),
        value in point

        @test abs(value / TEST_QUANTUM - round(value / TEST_QUANTUM)) < 1.0e-5
    end
    for corner in LOOP_END_CORNERS
        triangle = loop([corner.before, corner.point, corner.after])
        kinds, _, dots = semantic_corner_kinds(
            [(corner.point[1], corner.point[2], 0.0)],
            [triangle],
            1.0e-8,
            TEST_QUANTUM
        )
        @test kinds == [corner.kind]
        @test (dots[1] == 0) == (corner.kind === :legacy)
    end
    # The rotated perpendicular corner: the float dot product of the float side vectors is
    # NOT zero (-1.776e-15, the coordinates' rounding), the integer one is exactly 0.
    rotated = LOOP_END_CORNERS[4]
    a = collect(rotated.before) .- collect(rotated.point)
    b = collect(rotated.after) .- collect(rotated.point)
    @test a[1] * b[1] + a[2] * b[2] != 0.0
    @test plan_view_quantum_count.(a, TEST_QUANTUM) == [1052632000, -1000]
    @test plan_view_quantum_count.(b, TEST_QUANTUM) == [1000, 1052632000]
    # The same corner rotated rigidly to other angles (integer rotations on the quantum
    # grid: 3-4-5 and 5-12-13 triangles of 1e8 quanta) stays legacy; the real tilts read
    # 526316000 x 1000 (twice) and 6226890000 x 1000 quanta squared.
    for (c, s) in ((3, 4), (4, 3), (-4, 3), (5, 12))
        q = 100_000_000
        p = (0.38, -0.76)
        before = p .+ TEST_QUANTUM .* (c * q, s * q)
        after = p .+ TEST_QUANTUM .* (-s * q, c * q)
        kinds, _, dots = semantic_corner_kinds(
            [(p[1], p[2], 0.0)],
            [loop([before, p, after])],
            1.0e-8,
            TEST_QUANTUM
        )
        @test kinds == [:legacy] && dots[1] == 0
    end
    tilted = LOOP_END_CORNERS[5]
    _, _, dots = semantic_corner_kinds(
        [(tilted.point[1], tilted.point[2], 0.0)],
        [loop([tilted.before, tilted.point, tilted.after])],
        1.0e-8,
        TEST_QUANTUM
    )
    @test dots[1] == 526316000 * 1000
    _, _, dots = semantic_corner_kinds(
        [(LOOP_END_CORNERS[7].point[1], LOOP_END_CORNERS[7].point[2], 0.0)],
        [
            loop([
                LOOP_END_CORNERS[7].before,
                LOOP_END_CORNERS[7].point,
                LOOP_END_CORNERS[7].after
            ])
        ],
        1.0e-8,
        TEST_QUANTUM
    )
    @test dots[1] == 6226890000 * 1000
end

@testset "contract record of the invariant corners: present iff the mesher finds them" begin
    mktempdir() do directory
        corners = [(0.0, 0.0, 0.0), (0.6, 0.0, 0.0)]
        kinds = [:invariant, :legacy]
        path = joinpath(directory, "semantic.json")
        write(
            path,
            """{"Version": 1, "SemanticCorners": [[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]]}"""
        )
        recorded = read_invariant_corners(path, copy(IDENTITY_RIGID_TRANSFORM))
        @test isempty(recorded)
        # A contract without the record while the mesher finds an invariant corner: regenerate.
        @test_throws ErrorException check_invariant_corner_contract(
            corners,
            kinds,
            recorded,
            1.0e-8
        )
        @test check_invariant_corner_contract(
            corners,
            [:legacy, :legacy],
            recorded,
            1.0e-8
        ) == 0
        write(
            path,
            """{"Version": 1, "SemanticCorners": [[0.0, 0.0, 0.0], [0.6, 0.0, 0.0]],
                "Derivation": {"InvariantCorners": {"Rule": "r", "Points": [[0.0, 0.0, 0.0]]}}}"""
        )
        recorded = read_invariant_corners(path, copy(IDENTITY_RIGID_TRANSFORM))
        @test recorded == [(0.0, 0.0, 0.0)]
        @test check_invariant_corner_contract(corners, kinds, recorded, 1.0e-8) == 1
        # Recorded but perpendicular by the mesher: a disagreement fails closed as well.
        @test_throws ErrorException check_invariant_corner_contract(
            corners,
            [:legacy, :legacy],
            recorded,
            1.0e-8
        )
    end
end

# A 150-degree fabricated kink at the origin: wall 1 along +x, wall 2 at 150 degrees,
# both vertical (z in [-0.05, 0]); the S4 sliver (the corner, a vertex on the corner
# line, one on wall 1, the apex on wall 2) against a healthy cell with an interior apex.
@testset "bridging sliver predicate: all-surface cell spanning both sidewalls" begin
    d1 = [1.0, 0.0];
    d2 = [cosd(150.0), sind(150.0)]
    points = hcat(
        [0.0, 0.0, 0.0],                      # 1 the corner
        [0.0, 0.0, -0.00026],                 # 2 on the corner line (both walls)
        [0.0001, 0.0, -0.00013],              # 3 on wall 1
        [0.000179 * d2[1], 0.000179 * d2[2], -0.0001],  # 4 on wall 2
        [0.0002, 0.0, 0.0],                   # 5 on wall 1 (the MS plane edge)
        [0.0, 0.0, -0.0005],                  # 6 on the corner line, lower
        [0.0003 * d2[1], 0.0003 * d2[2], -0.0003],   # 7 on wall 2
        [0.00005, 0.00015, -0.0001],          # 8 interior (dielectric wedge)
        [0.0002, 0.0, -0.0004]
    )               # 9 on wall 1
    triangles = [
        (1, 2, 3),
        (2, 3, 9),
        (3, 5, 9),
        (1, 3, 5),  # wall 1
        (1, 2, 4),
        (2, 4, 7),
        (4, 7, 6),
        (2, 6, 7)
    ]  # wall 2
    tetrahedra = [
        (1, 2, 3, 4),   # the S4-type sliver: three on wall 1, the apex on wall 2
        (1, 3, 5, 8),   # wall 1 + an interior vertex: not bridging
        (1, 2, 3, 9),   # wall 1 only, all surface: not bridging
        (2, 3, 7, 8),   # both walls but an interior vertex: not bridging
        (1, 3, 4, 8)
    ]   # both walls, interior apex: not bridging
    surface = [i <= 7 || i == 9 for i = 1:9]
    # The metal fills the 210-degree sector (a concave metal corner): the sliver sits in
    # the 150-degree vacuum wedge between the walls.
    sides = (walls=[d1, d2], metal=(-[cosd(75.0), sind(75.0)]))
    center = [0.0, 0.0, 0.0]
    masks = sidewall_vertex_masks(points, triangles, center, sides, 0.1)
    @test masks[1] == 0x03 && masks[2] == 0x03 && masks[6] == 0x03
    @test masks[3] == 0x01 && masks[5] == 0x01 && masks[9] == 0x01
    @test masks[4] == 0x02 && masks[7] == 0x02
    @test !haskey(masks, 8)
    @test bridging_sliver_cells(points, tetrahedra, 1:5, masks, center, sides) == [1]
    # The same cells with the metal on the 150-degree side (a convex kink of 150 degrees,
    # the synthetic family / S4): the wedge containing the cell is still obtuse (150) and
    # the candidate stands; the dielectric 210-degree wedge is obtuse as well.
    convex = (walls=[d1, d2], metal=[cosd(75.0), sind(75.0)])
    @test bridging_sliver_cells(points, tetrahedra, 1:5, masks, center, convex) == [1]
    @test wedge_opening(convex, [cosd(75.0), sind(75.0)]) ≈ deg2rad(150.0)
    @test wedge_opening(convex, [cosd(-75.0), sind(-75.0)]) ≈ deg2rad(210.0)
    # A 22.5-degree fabricated tip: the same combinatorial cell type in the acute wedge under
    # the metal is a healthy fan cell, not a candidate (decision 365).
    tip = (walls=[[1.0, 0.0], [cosd(22.5), sind(22.5)]], metal=[cosd(11.25), sind(11.25)])
    tip_points = hcat(
        [0.0, 0.0, 0.0],
        [0.0, 0.0, -0.00025],
        [0.00025, 0.0, -0.0001],
        [0.00025 * cosd(22.5), 0.00025 * sind(22.5), -0.0001]
    )
    tip_masks = sidewall_vertex_masks(tip_points, [(1, 2, 3), (1, 2, 4)], center, tip, 0.1)
    @test tip_masks == Dict(1 => 0x03, 2 => 0x03, 3 => 0x01, 4 => 0x02)
    @test isempty(
        bridging_sliver_cells(tip_points, [(1, 2, 3, 4)], 1:1, tip_masks, center, tip)
    )
    @test wedge_opening(tip, [cosd(11.25), sind(11.25)]) ≈ deg2rad(22.5)
    # A cell with a vertex on a third surface (the MS plane, vertex 5 at z = 0 off the walls
    # ... here vertex 8 made a surface vertex) is not a candidate either.
    ms_masks = copy(masks)
    @test isempty(
        bridging_sliver_cells(points, [(1, 2, 3, 8)], 1:1, ms_masks, center, sides)
    )
    # A horizontal triangle (the MS plane) is not a wall; a thin coupon's kink has no
    # vertical triangle at all and hence no wall.
    flat = sidewall_vertex_masks(points, [(1, 5, 8)], center, sides, 0.1)
    @test isempty(flat)
    @test isempty(bridging_sliver_cells(points, tetrahedra, 1:5, flat, center, sides))
    @test tetrahedron_regular_condition([points[:, i] for i in tetrahedra[1]]) > 4.0
end

# A fan of tetrahedra at a corner whose three CAD edges span a 150-degree kink in the
# plane z = 0 (the metal edges) and the vertical: the corner cells need the invariant
# gate once the corner is declared invariant.
function kink_fan(; phi=150.0, ring=6, radius=0.025)
    points = zeros(3, ring + 2)
    for i = 1:ring
        angle = deg2rad(phi) * (i - 1) / (ring - 1)
        points[:, 1 + i] = [radius * cos(angle), radius * sin(angle), 0.0]
    end
    points[:, ring + 2] =
        [0.4 * radius * cosd(phi / 2), 0.4 * radius * sind(phi / 2), 0.3 * radius]
    tetrahedra = [(1, 1 + i, 2 + i, ring + 2) for i = 1:(ring - 1)]
    triangles = [(1, 1 + i, 2 + i) for i = 1:(ring - 1)]
    return points, tetrahedra, triangles
end

@testset "invariant corner in optimize_required_region!: gate required, kappa_reg judged" begin
    points, tetrahedra, triangles = kink_fan()
    corners = [(0.0, 0.0, 0.0)]
    sides =
        [(walls=[[1.0, 0.0], [cosd(150.0), sind(150.0)]], metal=[cosd(75.0), sind(75.0)])]
    spans = Tuple{Vector{Float64}, Vector{Float64}}[]
    run(kinds; gate=0.0) = optimize_required_region!(
        copy(points),
        copy(tetrahedra),
        triangles,
        corners,
        0.1,
        0.025,
        spans,
        0.0,
        2.0,
        0.0,
        0.05,
        4.0,
        0.01,
        1000.0,
        0.75,
        1.0e-9;
        corner_kinds=kinds,
        corner_sides=kinds == [:invariant] ? sides : [nothing],
        corner_shape_gate=gate
    )
    # Legacy: the production record, Kind Legacy, the vertex-0 measure.
    record, _ = run([:legacy])
    legacy = record["CornerMeasures"][1]
    @test legacy["Kind"] == "Legacy" && legacy["Measure"] == "VertexFrameCondition"
    @test legacy["Target"] == 3.8 && legacy["Gate"] == 4.0 && legacy["Passed"]
    @test record["InvariantCorners"] == 0 && record["CornerShapeGate"] === nothing
    @test haskey(legacy["Information"], "KappaRegMax") &&
          legacy["Information"]["EtaMin"] > 0.0
    @test record["CornerAspectsAfter"] == [legacy["After"]]
    # Invariant without a gate: fail closed.
    @test_throws ErrorException run([:invariant])
    record, _ = run([:invariant]; gate=5.0)
    invariant = record["CornerMeasures"][1]
    @test invariant["Kind"] == "Invariant" && invariant["Measure"] == "RegularCondition"
    @test invariant["Target"] == 3.8 && invariant["Gate"] == 5.0 && invariant["Passed"]
    @test invariant["BridgingSlivers"] ==
          Dict("Before" => 0, "After" => 0, "AboveGateAfter" => 0)
    @test invariant["After"] <= 5.0 && invariant["After"] == record["CornerAspectsAfter"][1]
    @test record["InvariantCorners"] == 1 && record["CornerShapeGate"] == 5.0
    @test record["InvariantCornerTarget"] == 3.8
    @test occursin("regular tetrahedron", record["InvariantCornerRule"])
    # A gate below the achieved value fails closed with the invariant message.
    if invariant["After"] > 1.5
        @test_throws ErrorException run([:invariant]; gate=1.5)
    end
end
