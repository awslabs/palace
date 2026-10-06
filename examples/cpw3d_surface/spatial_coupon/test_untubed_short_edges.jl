# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# Short metal edges against the corner clearances (block (b) design AMENDMENT 1 A9
# family 3; supervisor decision 304): a straight metal side whose tube interval is
# below the tube's inner ring size carries no tube; its corner balls mesh it at
# CornerSize; the census Scope records UntubedEdges[] {Side, Span, Clearances,
# Interval, Corners, CoveredByBalls} and TubeCount drops by TubesPerSide per untubed
# side. Before this change the side failed closed as ScopeGuard[ShortEdges].

guard_message(f) =
    try
        f()
        ""
    catch e
        e isa ErrorException ? e.msg : sprint(showerror, e)
    end

@testset "metal edge segments: the untubed short-side class and its knife edge" begin
    # A metal finger of width w between two 90-degree corners ((2, -w/2) -> (2, w/2) is
    # the tip side); the metal fills x <= 0 and the finger [0, 2] x [-w/2, w/2]; the box
    # [-1, 3] x [-1, 1] closes the loop. Clearance c(90) = 0.05 at both tip corners, so
    # the tip's interval is w - 0.1 against the inner ring size 0.01.
    lower = [-1.0, -1.0]
    upper = [3.0, 1.0]
    clearance(angle) = 0.03 / tan(0.5 * angle) + max(0.02, 1.25 * 0.01 / sin(0.5 * angle))
    @test clearance(pi / 2) ≈ 0.05
    function finger_segments(w; edge_size=0.01, corner_radius=0.1, clearance=clearance)
        points = [
            (lower[1], lower[2]),
            (0.0, lower[2]),
            (0.0, -w / 2),
            (2.0, -w / 2),
            (2.0, w / 2),
            (0.0, w / 2),
            (0.0, upper[2]),
            (lower[1], upper[2])
        ]
        classes = [
            "Continuation",
            "Physical",
            "Physical",
            "Physical",
            "Physical",
            "Physical",
            "Continuation",
            "Continuation"
        ]
        loop = (conductor=1, plane=0.0, hole=false, points=points, classes=classes)
        # Decision 320: (0, lower.y) is a Physical-class theta-0 box vertex (a legacy corner).
        corners = [(p[1], p[2], 0.0) for p in points[2:6]]
        return metal_edge_segments(
            [loop],
            corners,
            clearance,
            lower,
            upper,
            1.0e-9;
            edge_size=edge_size,
            corner_radius=corner_radius
        )
    end
    # Above the knife edge: interval 0.02 >= 0.01 -> a one-layer tube of two inner rings.
    wide = finger_segments(0.12)
    @test length(wide) == 5 && all(!segment.untubed for segment in wide)
    tip = wide[3]
    @test tip.start == [2.0, -0.06] && tip.stop == [2.0, 0.06]
    @test tip.s_start ≈ 0.05 && tip.s_end ≈ 0.07
    # Exactly at the knife edge (binary-exact numbers: clearance 0.25, EdgeSize 2^-6, width
    # 0.515625): interval == EdgeSize carries a tube; one ulp of width below it does not.
    exact = finger_segments(0.515625; edge_size=0.015625, clearance=angle -> 0.25)
    @test !exact[3].untubed && exact[3].s_end - exact[3].s_start == 0.015625
    @test finger_segments(prevfloat(0.515625); edge_size=0.015625, clearance=angle -> 0.25)[3].untubed
    # Below it: interval 0.009 -> untubed; the clearances and the span are kept as derived,
    # the two 0.1 balls cover the 0.109 span (CoveredByBalls).
    below = finger_segments(0.109)
    @test [segment.untubed for segment in below] == [false, false, true, false, false]
    @test below[3].s_start ≈ 0.05 && below[3].span - below[3].s_end ≈ 0.05
    @test below[3].corners == (true, true) && below[3].covered_by_balls
    # A negative interval (the former ShortEdges refusal) is untubed too, and the message of
    # the former guard is gone: no exception.
    narrow = finger_segments(0.08)
    @test narrow[3].untubed && narrow[3].s_end - narrow[3].s_start ≈ -0.02
    @test narrow[3].covered_by_balls
    @test guard_message(() -> finger_segments(0.08)) == ""
    # Coverage follows the balls that exist: with a 0.05 ball radius the two balls reach
    # 0.1 < 0.109 -> CoveredByBalls false (the middle is graded by the two-sided law).
    @test !finger_segments(0.109; corner_radius=0.05)[3].covered_by_balls
    # The finger's long sides (span 2 at any width) are unchanged bitwise by the class.
    for k in (2, 4)
        @test wide[k].s_start == below[k].s_start && wide[k].s_end == below[k].s_end
        @test !below[k].untubed
    end
    # Without an inner ring size (edge_size 0) an exactly empty interval is a degenerate
    # tube and fails closed; a negative one is the untubed class.
    @test occursin(
        "leaves no tube interval",
        guard_message(() -> finger_segments(0.1; edge_size=0.0))
    )
    @test finger_segments(0.08; edge_size=0.0)[3].untubed
end

@testset "an untubed side between two acute corners leaves an uncovered middle (CoveredByBalls false)" begin
    # An isosceles metal island with 20-degree base angles: the base of span 0.45 lies
    # between two acute corners whose clearances (0.03 / tan 10 + 1.25 x 0.01 / sin 10 =
    # 0.242 each) exceed the 0.1 ball radius: untubed, and 2 x 0.1 < 0.45.
    base = 0.45
    height = 0.5 * base * tand(20.0)
    points = [(0.0, 0.0), (base, 0.0), (0.5 * base, height)]
    loop = (conductor=1, plane=0.0, hole=false, points=points, classes=fill("Physical", 3))
    corners = [(p[1], p[2], 0.0) for p in points]
    clearance(angle) = 0.03 / tan(0.5 * angle) + max(0.02, 1.25 * 0.01 / sin(0.5 * angle))
    segments = metal_edge_segments(
        [loop],
        corners,
        clearance,
        [-2.0, -2.0],
        [2.0, 2.0],
        1.0e-9;
        edge_size=0.01,
        corner_radius=0.1
    )
    island_base = segments[1]
    @test island_base.span ≈ base && island_base.corners == (true, true)
    @test all(angle -> isapprox(angle, deg2rad(20.0)), island_base.corner_angles)
    @test island_base.s_start ≈ clearance(deg2rad(20.0)) && island_base.s_start > 0.1
    @test island_base.untubed && !island_base.covered_by_balls
end

# A metal finger of width `width` (um) protruding 2 um from a half-plane of metal into
# the coupon: the finger's tip side is the short side between two 90-degree corners.
# The box is the mesher's own (the two finger rows padded by Radius 0.5); the boundary
# classes are the plan-view builder's (a vertex carries its outgoing side's class) and
# the contract corners follow decision 320 (the theta-0 Physical-class box vertex at
# the bottom of the x = 0 side is a legacy corner).
function write_finger_inputs(directory, width; plane=0.0)
    length_x = 2.0
    rows = [
        (
            point=(0.5 * length_x, sign * 0.5 * width, plane),
            tangent=(1.0, 0.0, 0.0),
            gap=(0.0, sign, 0.0),
            interval=(-0.5 * length_x, 0.5 * length_x),
            normal_sign=1.0,
            vertex_arm=false,
            slot=0,
            conductor=1
        ) for sign in (-1.0, 1.0)
    ]
    lower, upper = row_coupon_bounds(rows, 0.5, 0.1, 0.05)
    polygon = [
        (lower[1], lower[2]),
        (0.0, lower[2]),
        (0.0, -0.5 * width),
        (length_x, -0.5 * width),
        (length_x, 0.5 * width),
        (0.0, 0.5 * width),
        (0.0, upper[2]),
        (lower[1], upper[2])
    ]
    classes = [
        "Continuation",
        "Physical",
        "Physical",
        "Physical",
        "Physical",
        "Physical",
        "Continuation",
        "Continuation"
    ]
    corners = [(p[1], p[2], plane) for p in polygon[2:6]]
    open(joinpath(directory, "signature.csv"), "w") do io
        println(io, "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm")
        for (i, row) in enumerate(rows)
            println(
                io,
                join(
                    [
                        i,
                        0,
                        1,
                        row.point...,
                        row.gap...,
                        row.tangent...,
                        1,
                        row.interval...,
                        0
                    ],
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "boundary.csv"), "w") do io
        println(io, "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y")
        for (i, point) in enumerate(polygon)
            println(io, join([1, i, 1, plane, 0, classes[i], point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for point in polygon
            println(io, join([1, 1, plane, point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "semantic.json"), "w") do io
        return write_json(
            io,
            Dict{String, Any}(
                "Version" => 1,
                "SemanticCorners" => [[c[1], c[2], c[3]] for c in corners]
            )
        )
    end
    return (
        signature=joinpath(directory, "signature.csv"),
        boundary=joinpath(directory, "boundary.csv"),
        mask=joinpath(directory, "mask.csv"),
        semantic=joinpath(directory, "semantic.json"),
        polygon=polygon,
        corners=corners,
        lower=lower,
        upper=upper,
        tip=((length_x, -0.5 * width), (length_x, 0.5 * width))
    )
end

function build_finger_coupon(
    directory,
    width;
    fabricated=true,
    stem="finger",
    overetch=0.05
)
    inputs = write_finger_inputs(directory, width)
    mesh = joinpath(directory, "coupon-$stem.msh")
    census = joinpath(directory, "census-$stem.json")
    generate_spatial_coupon(;
        signature=inputs.signature,
        mask=inputs.mask,
        boundary=inputs.boundary,
        fabricated=fabricated,
        filename=mesh,
        radius=0.5,
        metal_thickness=0.1,
        overetch=overetch,
        sidewall_angle=90.0,
        top_rounding=0.0,
        trench_rounding=0.0,
        lc_fine=0.05,
        lc_tangent=0.1,
        lc_far=0.3,
        max_nodes=2_000_000,
        max_elements=2_000_000,
        semantic_contract=inputs.semantic,
        corner_isotropy_radius=0.1,
        corner_census=census,
        edge_size=0.01,
        edge_growth_ratio=2.0,
        corner_size=0.01,
        prism_tubes=true,
        far_growth=0.5,
        maximum_corner_aspect=4.0,
        minimum_scaled_jacobian=0.01,
        maximum_jacobian_condition=1000.0,
        quality_displacement_over_normal=0.75
    )
    return parse_json(read(census, String)), mesh, inputs
end

# Mesh nodes within `reach` of `point`.
function nodes_near(path, point, reach)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(path)
    _, coordinates, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()
    xyz = reshape(coordinates, 3, :)
    return [xyz[:, i] for i = 1:size(xyz, 2) if norm(xyz[:, i] .- point) <= reach]
end

untubed_tip(census, inputs) = [
    row for row in census["Scope"]["UntubedEdges"] if (
        row["Side"]["Start"] == collect(inputs.tip[1]) &&
        row["Side"]["Stop"] == collect(inputs.tip[2])
    ) || (
        row["Side"]["Start"] == collect(inputs.tip[2]) &&
        row["Side"]["Stop"] == collect(inputs.tip[1])
    )
]

@testset "synthetic 0.09-um finger: the tip side is untubed (thin sheet tube; fabricated top + bottom tubes); a 0.12 finger keeps a one-layer tube" begin
    mktempdir() do directory
        # Tubes at these sizes (EdgeSize 0.01, ratio 2): the thin sheet tube's bound is the
        # 0.1 ball radius and the fabricated one's min(Overetch 0.05, thickness / 2, 0.1) ->
        # both 2 rings (0.01, 0.02): R 0.03, h_K 0.02, h_pyr 0.01, c(90) = 0.05; the 0.09 tip
        # leaves -0.01 < 0.01 -> untubed, covered by the two 0.1 balls (2 x 0.1 >= 0.09).
        thin, thin_mesh, inputs =
            build_finger_coupon(directory, 0.09; fabricated=false, stem="thin")
        scope = thin["Scope"]
        @test "UntubedShortEdges" in scope["SupportedClasses"]
        @test !("ShortEdges" in scope["GuardedClasses"])
        @test occursin("minus the untubed short sides", scope["Rule"])
        @test occursin("UntubedEdges", scope["UntubedShortEdgeRule"])
        @test length(scope["UntubedEdges"]) == 1
        record = untubed_tip(thin, inputs)[1]
        @test record["Span"] ≈ 0.09 && record["Clearances"] ≈ [0.05, 0.05]
        @test record["Interval"] ≈ -0.01 && record["Corners"] == [true, true]
        @test record["CoveredByBalls"] === true
        @test record["Side"]["Plane"] == 0.0 &&
              record["Side"]["Conductor"] == 1 &&
              record["Side"]["Layer"] == 1
        tubes = thin["PrismTubes"]
        # 5 straight sides, one untubed -> 4 sheet tubes.
        @test sum(loop["Sides"] for loop in scope["MetalLoops"]) == 5
        @test tubes["TubeCount"] == 4 && length(tubes["Tubes"]) == 4
        # No tube runs along the tip (x = 2, extrusion along y).
        @test !any(
            abs(row["Extrusion"][2]) > 0.5 && abs(row["Origin"][1] - 2.0) <= 1.0e-9 for
            row in tubes["Tubes"]
        )
        @test occursin("untubed short sides", tubes["Rule"])
        # The tip is meshed by the corner balls at CornerSize: the nodes on the tip line are
        # spaced below the ball's innermost shells, none farther than the ball radius apart.
        tip_nodes = [
            p for p in nodes_near(thin_mesh, [2.0, 0.0, 0.0], 0.06) if
            abs(p[1] - 2.0) <= 1.0e-9 && abs(p[3]) <= 1.0e-9
        ]
        ys = sort([p[2] for p in tip_nodes])
        @test length(ys) >= 5 && ys[1] ≈ -0.045 && ys[end] ≈ 0.045
        @test maximum(diff(ys)) <= 0.03 + 1.0e-9
        @test tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        @test tubes["Quality"]["Prism"]["PositiveOrientation"]
        # Fabricated at the same tube sizes: the tip's top AND bottom tube are dropped
        # (TubeCount 10 -> 8), the tip edges at z = 0 and z = 0.1 meshed by the balls.
        fab, fab_mesh, _ = build_finger_coupon(directory, 0.09; fabricated=true, stem="fab")
        @test length(fab["Scope"]["UntubedEdges"]) == 1 &&
              untubed_tip(fab, inputs)[1]["Interval"] ≈ -0.01
        @test fab["PrismTubes"]["TubeCount"] == 8 && length(fab["PrismTubes"]["Tubes"]) == 8
        @test fab["PrismTubes"]["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        @test fab["PrismTubes"]["Quality"]["Prism"]["PositiveOrientation"] &&
              fab["PrismTubes"]["Quality"]["Pyramid"]["PositiveOrientation"]
        # (the ball's shells 0.01 / 0.03 / 0.07: no node interval on the tip line exceeds the
        # third shell size 0.04)
        for z in (0.0, 0.1)
            tip_nodes = [
                p for p in nodes_near(fab_mesh, [2.0, 0.0, z], 0.06) if
                abs(p[1] - 2.0) <= 1.0e-9 && abs(p[3] - z) <= 1.0e-9
            ]
            ys = sort([p[2] for p in tip_nodes])
            @test length(ys) >= 4
            @test ys[1] ≈ -0.045 && ys[end] ≈ 0.045
            @test maximum(diff(ys)) <= 0.04 + 1.0e-9
        end
        # The control above the knife edge: a 0.12 finger leaves 0.02 >= 0.01 -> the tip
        # carries a one-layer sheet tube of 0.02 (two inner rings in length; the S1p v3
        # fingers' 15-22-nm thin tubes are this case): nothing untubed, TubeCount 5.
        control, _, _ =
            build_finger_coupon(directory, 0.12; fabricated=false, stem="control")
        @test control["Scope"]["UntubedEdges"] == []
        @test control["PrismTubes"]["TubeCount"] == 5
        short_rows = [
            row for
            row in control["PrismTubes"]["Tubes"] if abs(row["Length"] - 0.02) <= 1.0e-9
        ]
        @test length(short_rows) == 1 && short_rows[1]["Layers"] == 1
    end
end

@testset "NarrowMetal (decision 347): tubed sides facing each other across a metal strip narrower than 2 x the envelope" begin
    lower = [-1.0, -1.0]
    upper = [3.0, 1.0]
    clearance(angle) = 0.03 / tan(0.5 * angle) + max(0.02, 1.25 * 0.01 / sin(0.5 * angle))
    # One finger of width w: its two long sides face each other across w of metal (their
    # tube intervals run from 0.05 past the root to 0.05 before the tip).
    function finger(w)
        points = [
            (lower[1], lower[2]),
            (0.0, lower[2]),
            (0.0, -w / 2),
            (2.0, -w / 2),
            (2.0, w / 2),
            (0.0, w / 2),
            (0.0, upper[2]),
            (lower[1], upper[2])
        ]
        classes = [
            "Continuation",
            "Physical",
            "Physical",
            "Physical",
            "Physical",
            "Physical",
            "Continuation",
            "Continuation"
        ]
        loop = (conductor=1, plane=0.0, hole=false, points=points, classes=classes)
        return metal_edge_segments(
            [loop],
            [(p[1], p[2], 0.0) for p in points[2:6]],
            clearance,
            lower,
            upper,
            1.0e-9;
            edge_size=0.01,
            corner_radius=0.1
        )
    end
    @test metal_facing_width(finger(0.07), 0.03, 1.0e-9) ≈ 0.07
    message = guard_message(() -> metal_facing_width(finger(0.07), 0.04, 1.0e-9))
    @test occursin("ScopeGuard[NarrowMetal]", message) &&
          occursin("0.07", message) &&
          occursin("0.04", message)
    @test metal_facing_width(finger(0.09), 0.04, 1.0e-9) ≈ 0.09
    # Two 0.2 fingers separated by a 0.05 dielectric slot: the slot sides face each other
    # across DIELECTRIC (their outward normals point at each other) and are not a metal
    # strip; the fingers' own widths (0.2) are the facing metal widths.
    slot = 0.05
    points = [
        (lower[1], lower[2]),
        (0.0, lower[2]),
        (0.0, -slot / 2 - 0.2),
        (2.0, -slot / 2 - 0.2),
        (2.0, -slot / 2),
        (0.0, -slot / 2),
        (0.0, slot / 2),
        (2.0, slot / 2),
        (2.0, slot / 2 + 0.2),
        (0.0, slot / 2 + 0.2),
        (0.0, upper[2]),
        (lower[1], upper[2])
    ]
    classes = vcat("Continuation", fill("Physical", 9), "Continuation", "Continuation")
    loop = (conductor=1, plane=0.0, hole=false, points=points, classes=classes)
    two = metal_edge_segments(
        [loop],
        [(p[1], p[2], 0.0) for p in points[2:10]],
        clearance,
        lower,
        upper,
        1.0e-9;
        edge_size=0.01,
        corner_radius=0.1
    )
    # (the 0.05 slot root between its two concave corners is an untubed short side)
    @test length(two) == 9 && count(segment.untubed for segment in two) == 1
    @test two[5].untubed && two[5].span ≈ slot
    @test metal_facing_width(two, 0.04, 1.0e-9) ≈ 0.2
    # An untubed side faces nothing: the 0.09 finger's tip is untubed, its long sides keep
    # the 0.09 facing width; a side list without a facing pair reads Inf.
    @test metal_facing_width(finger(0.09)[[1, 5]], 0.04, 1.0e-9) == Inf
    @test "NarrowMetal" in [guard[1] for guard in RECIPE_SCOPE_GUARDS]
    # The full thin build of a 0.07 finger at these sizes (sheet tube R 0.03 + h_pyr 0.01 =
    # 0.04 envelope: 2 x 0.04 > 0.07) stops at the guard before any CAD tube is built; the
    # 0.09 finger of the testset above builds.
    mktempdir() do directory
        message = guard_message(
            () -> build_finger_coupon(directory, 0.07; fabricated=false, stem="narrow")
        )
        @test occursin("ScopeGuard[NarrowMetal]", message)
        @test occursin("across 0.07 of metal", message)
        envelope = match(r"tube envelope ([0-9.e-]+)", message)
        @test envelope !== nothing && parse(Float64, envelope[1]) ≈ 0.04
    end
end
