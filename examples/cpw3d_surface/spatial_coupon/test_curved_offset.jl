# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Exact curved plan-view offsets and the collar island rule (supervisor decision 246).
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

const TOLERANCE = 1.0e-9

exterior(points; classes=fill("Physical", length(points))) =
    (conductor=1, plane=0.0, hole=false, points=points, classes=classes)
hole_loop(points; classes=fill("Physical", length(points))) =
    (conductor=1, plane=0.0, hole=true, points=points, classes=classes)

# Counterclockwise rounded rectangle [x0, x1] x [y0, y1] with quarter arcs of radius r
# (`chords` chords each); the straight sides join consecutive arcs tangentially.
function rounded_rectangle(x0, y0, x1, y1, r; chords=9, perturbation=0.0)
    centers = (
        (x1 - r, y0 + r, -pi / 2),
        (x1 - r, y1 - r, 0.0),
        (x0 + r, y1 - r, pi / 2),
        (x0 + r, y0 + r, pi)
    )
    points = NTuple{2, Float64}[]
    for (index, (cx, cy, a0)) in enumerate(centers), k = 0:chords
        a = a0 + k * (pi / 2) / chords + (index == 1 ? perturbation : 0.0)
        push!(points, (cx + r * cos(a), cy + r * sin(a)))
    end
    return points
end

# Counterclockwise circle (one closed arc run).
circle(cx, cy, r; chords=36) = [
    (cx + r * cos(2pi * k / chords), cy + r * sin(2pi * k / chords)) for k = 0:(chords - 1)
]

guard_message(f) =
    try
        f()
        ""
    catch e
        e.msg
    end

distance_to_circle(p, center, r) = abs(hypot(p[1] - center[1], p[2] - center[2]) - r)

# The reference (pre-decision-246) straight implementations, verbatim, for the
# bitwise straight-coupon comparison.
function reference_offset_loop_points(loop, distance, tolerance)
    abs(distance) <= tolerance && return loop.points
    metal_side = loop_orientation(loop.points) * (loop.hole ? -1.0 : 1.0)
    shifted = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for index in eachindex(loop.points)
        first = loop.points[index]
        second = loop.points[mod1(index + 1, length(loop.points))]
        direction = (second[1] - first[1], second[2] - first[2])
        segment_length = hypot(direction...)
        shift = loop.classes[index] == "Physical" ? distance : 0.0
        normal = (
            -metal_side * direction[2] / segment_length,
            metal_side * direction[1] / segment_length
        )
        push!(
            shifted,
            ((first[1] + shift * normal[1], first[2] + shift * normal[2]), direction)
        )
    end
    points = NTuple{2, Float64}[]
    for index in eachindex(shifted)
        previous = shifted[mod1(index - 1, length(shifted))]
        current = shifted[index]
        denominator = cross2d(previous[2], current[2])
        offset = (current[1][1] - previous[1][1], current[1][2] - previous[1][2])
        coordinate = cross2d(offset, current[2]) / denominator
        push!(
            points,
            (
                previous[1][1] + coordinate * previous[2][1],
                previous[1][2] + coordinate * previous[2][2]
            )
        )
    end
    return points
end

function reference_offset_hole_points(loop, distance, tolerance)
    distance >= -tolerance && return reference_offset_loop_points(loop, distance, tolerance)
    points = copy(loop.points)
    orientation = loop_orientation(points)
    clipped = copy(points)
    for i in eachindex(points)
        loop.classes[i] == "Physical" || continue
        a, b = points[i], points[mod1(i + 1, length(points))]
        direction = (b[1] - a[1], b[2] - a[2]);
        edge_length = hypot(direction...)
        normal = (
            -orientation * direction[2] / edge_length,
            orientation * direction[1] / edge_length
        )
        signed(p) = normal[1] * (p[1] - a[1]) + normal[2] * (p[2] - a[2]) + distance
        result = NTuple{2, Float64}[]
        isempty(clipped) && return result
        for j in eachindex(clipped)
            p, q = clipped[j], clipped[mod1(j + 1, Base.length(clipped))]
            dp, dq = signed(p), signed(q)
            dp >= -tolerance && push!(result, p)
            if (dp >= -tolerance) != (dq >= -tolerance)
                t = dp / (dp - dq)
                push!(result, (p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2])))
            end
        end
        clipped = result
    end
    cleaned = NTuple{2, Float64}[]
    for p in clipped
        (
            isempty(cleaned) ||
            hypot(p[1] - cleaned[end][1], p[2] - cleaned[end][2]) > tolerance
        ) && push!(cleaned, p)
    end
    if Base.length(cleaned) > 1 &&
       hypot(cleaned[1][1] - cleaned[end][1], cleaned[1][2] - cleaned[end][2]) <= tolerance
        pop!(cleaned)
    end
    Base.length(cleaned) >= 3 || return NTuple{2, Float64}[]
    area2 = sum(
        cross2d(cleaned[i], cleaned[mod1(i + 1, Base.length(cleaned))]) for
        i in eachindex(cleaned)
    )
    return abs(area2) > tolerance^2 ? cleaned : NTuple{2, Float64}[]
end

@testset "annulus: whole-circle loops offset concentrically, a collapsed circular hole vanishes" begin
    disc = exterior(circle(1.0, -2.0, 2.0))
    for distance in (-0.5, 0.75)
        offset = offset_loop_points(disc, distance, TOLERANCE)
        @test length(offset) == 36
        @test all(
            distance_to_circle(p, (1.0, -2.0), 2.0 - distance) <= 1.0e-12 for p in offset
        )
        runs = circular_arc_runs(offset, TOLERANCE)
        @test length(runs) == 1 &&
              isapprox(runs[1].radius, 2.0 - distance; atol=1.0e-12) &&
              length(runs[1].edge_indices) == 36
    end
    @test occursin(
        "collapses the circular",
        guard_message(() -> offset_loop_points(disc, 2.0, TOLERANCE))
    )
    ring = hole_loop(circle(1.0, -2.0, 1.0))
    shrunk = offset_hole_points(ring, -0.3, TOLERANCE)
    @test length(shrunk) == 36 &&
          all(distance_to_circle(p, (1.0, -2.0), 0.7) <= 1.0e-12 for p in shrunk)
    @test isapprox(circular_arc_runs(shrunk, TOLERANCE)[1].radius, 0.7; atol=1.0e-12)
    grown = offset_hole_points(ring, 0.2, TOLERANCE)
    @test all(distance_to_circle(p, (1.0, -2.0), 1.2) <= 1.0e-12 for p in grown)
    @test isempty(offset_hole_points(ring, -1.0, TOLERANCE))
    @test isempty(offset_hole_points(ring, -1.5, TOLERANCE))
    # An annulus lofted as the mask would be: the disc's collar and the hole's shrink
    # are the same concentric construction.
    @test all(
        distance_to_circle(p, (1.0, -2.0), 2.3) <= 1.0e-12 for
        p in offset_loop_points(disc, -0.3, TOLERANCE)
    )
end

@testset "quarter-arc corners: convex arcs grow by the collar, straight sides shift, joints stay continuous" begin
    r = 1.0
    loop = exterior(rounded_rectangle(-4.0, -3.0, 4.0, 3.0, r))
    for distance in (-1.5, 0.4)
        offset = offset_loop_points(loop, distance, TOLERANCE)
        @test length(offset) == length(loop.points)
        runs = circular_arc_runs(offset, TOLERANCE)
        @test length(runs) == 4
        for run in runs
            @test isapprox(run.radius, r - distance; atol=1.0e-9)
            @test any(
                hypot(run.center[1] - c[1], run.center[2] - c[2]) <= 1.0e-9 for
                c in ((3.0, -2.0), (3.0, 2.0), (-3.0, 2.0), (-3.0, -2.0))
            )
            @test isapprox(run.angle, pi / 2; atol=1.0e-9)
            # The run's end points are exact points of the offset circle: the junction
            # with the tangent straight side is continuous.
            for index in run.point_indices
                @test distance_to_circle(offset[index], run.center, r - distance) <= 1.0e-9
            end
        end
        # The straight sides are the originals shifted by the offset.
        @test any(isapprox(p[2], 3.0 - distance; atol=1.0e-12) for p in offset)
        @test any(isapprox(p[1], -4.0 + distance; atol=1.0e-12) for p in offset)
        @test polygon_is_simple(offset, TOLERANCE)
        points, construction =
            collar_loop_points(loop, distance, ([-10.0, -10.0], [10.0, 10.0]), TOLERANCE)
        @test construction == "MiterOffset" && points == offset
    end
    # physical_segments (the rounding primitives) see the offset arcs.
    primitives = physical_segments([loop], -1.5, TOLERANCE)
    @test count(p -> p.kind == :arc, primitives) == 4 &&
          count(p -> p.kind == :line, primitives) == 4
    @test all(isapprox(p.radius, 2.5; atol=1.0e-9) for p in primitives if p.kind == :arc)
end

@testset "degenerate inner arc (r < offset) collapses onto its neighbours' miter" begin
    # An L-shaped conductor whose inner (concave-metal) corner carries a fillet of
    # radius 0.5: metal = [0, 4] x [0, 1] union [0, 1] x [0, 4].
    r = 0.5
    chords = 9
    # The fillet runs from (1 + r, 1) to (1, 1 + r) about (1 + r, 1 + r), clockwise:
    # a concave-metal arc of the counterclockwise loop (0,0) -> (4,0) -> (4,1) ->
    # fillet -> (1,4) -> (0,4).
    fillet = [
        (1.0 + r - r * cos(a), 1.0 + r - r * sin(a)) for
        a in range(pi / 2, 0.0; length=chords + 1)
    ]
    points = vcat([(0.0, 0.0), (4.0, 0.0), (4.0, 1.0)], fillet, [(1.0, 4.0), (0.0, 4.0)])
    loop = exterior(points)
    runs = circular_arc_runs(points, TOLERANCE)
    @test length(runs) == 1 && isapprox(runs[1].radius, r; atol=1.0e-12)
    sharp =
        exterior([(0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (1.0, 1.0), (1.0, 4.0), (0.0, 4.0)])
    # Collar 1.5 > r: the fillet collapses; the polygon is the sharp L's miter polygon.
    collapsed = offset_loop_points(loop, -1.5, TOLERANCE)
    @test length(collapsed) == 6
    reference = reference_offset_loop_points(sharp, -1.5, TOLERANCE)
    @test all(
        minimum(hypot(p[1] - q[1], p[2] - q[2]) for q in reference) <= 1.0e-12 for
        p in collapsed
    )
    @test any(
        isapprox(p[1], 2.5; atol=1.0e-12) && isapprox(p[2], 2.5; atol=1.0e-12) for
        p in collapsed
    )
    # Collar 0.3 < r: the fillet shrinks concentrically to 0.2.
    shrunk = offset_loop_points(loop, -0.3, TOLERANCE)
    @test length(shrunk) == length(points)
    run = only(circular_arc_runs(shrunk, TOLERANCE))
    @test isapprox(run.radius, r - 0.3; atol=1.0e-9) &&
          hypot(run.center[1] - 1.5, run.center[2] - 1.5) <= 1.0e-9
    # Exactly at r the arc is degenerate too (radius <= tolerance): collapsed.
    @test length(offset_loop_points(loop, -r, TOLERANCE)) == 6
    # An inward offset larger than a convex fillet collapses it onto the sharp corner.
    rounded = exterior(rounded_rectangle(-4.0, -3.0, 4.0, 3.0, 0.3))
    inward = offset_loop_points(rounded, 0.5, TOLERANCE)
    @test length(inward) == 4
    @test Set(round.(p; digits=9) for p in inward) ==
          Set([(3.5, -2.5), (3.5, 2.5), (-3.5, 2.5), (-3.5, -2.5)])
end

@testset "tangent arc-line joints under the identification's tangency noise" begin
    r = 3.0
    for perturbation in (0.0, 1.0e-6, -2.0e-6)
        loop = exterior(
            rounded_rectangle(
                -10.0,
                -8.0,
                10.0,
                8.0,
                r;
                chords=18,
                perturbation=perturbation
            )
        )
        runs = circular_arc_runs(loop.points, TOLERANCE)
        @test length(runs) == 4
        offset = curved_offset_loop(loop, -5.7, runs, TOLERANCE)
        @test !offset.bridged
        for junction in offset.junctions
            before = offset.items[junction.before]
            after = offset.items[junction.after]
            line = before.kind == :line ? before : after
            arc = before.kind == :arc ? before : after
            # The junction point lies on the shifted line exactly and on the offset
            # circle to the tangency noise squared.
            u = unit2d(line.direction)
            @test abs(
                cross2d(
                    u,
                    (
                        junction.point[1] - line.shifted_start[1],
                        junction.point[2] - line.shifted_start[2]
                    )
                )
            ) <= 1.0e-12
            @test distance_to_circle(junction.point, arc.center, arc.radius) <= 1.0e-9
            @test !junction.convex
            @test hypot(
                before.shifted_stop[1] - after.shifted_start[1],
                before.shifted_stop[2] - after.shifted_start[2]
            ) <= 1.0e-4
        end
        # polygon_wire would rebuild four exact arcs of radius r + 5.7.
        rebuilt = circular_arc_runs(offset.points, TOLERANCE)
        @test length(rebuilt) == 4 &&
              all(isapprox(run.radius, r + 5.7; atol=1.0e-9) for run in rebuilt)
        @test polygon_is_simple(offset.points, TOLERANCE)
    end
    # Two arcs meeting at a kink (a lens of two quarter arcs) fail closed.
    lens = vcat(
        [(2.0 * cos(a), 2.0 * sin(a)) for a in range(0.0, pi / 2; length=10)],
        [(2.0 + 2.0 * cos(a), 2.0 + 2.0 * sin(a)) for a in range(pi, 3pi / 2; length=10)][2:(end - 1)]
    )
    @test length(circular_arc_runs(lens, TOLERANCE)) == 2
    @test occursin(
        "kink",
        guard_message(() -> offset_loop_points(exterior(lens), -0.5, TOLERANCE))
    )
end

@testset "overlapping collars across a 2-um channel with a rounded end: union, collapsed arc bridged" begin
    box = ([-10.0, -10.0], [10.0, 10.0])
    radius = 1.9
    collar = 3radius
    # Metal y >= 0 of the box, cut by a channel of width `width` from the box side
    # y = 0 up to a semicircular end (centre (0, 6 - width / 2)).
    function channelled(width)
        half = width / 2
        cy = 8.0 - half
        # Over the top from the left wall's end to the right wall's end (clockwise
        # about the centre: a concave-metal arc of the counterclockwise loop).
        chain = [(-half * cos(a), cy + half * sin(a)) for a in range(0.0, pi; length=19)]
        points = vcat(
            [(-10.0, 0.0), (-half, 0.0)],
            chain,
            [(half, 0.0), (10.0, 0.0), (10.0, 10.0), (-10.0, 10.0)]
        )
        classes = vcat(fill("Physical", length(points) - 3), fill("Continuation", 3))
        return exterior(points; classes=classes)
    end
    narrow = channelled(2.0)
    @test loop_orientation(narrow.points) == 1.0
    runs = circular_arc_runs(narrow.points, TOLERANCE)
    @test length(runs) == 1 &&
          isapprox(runs[1].radius, 1.0; atol=1.0e-9) &&
          isapprox(runs[1].angle, pi; atol=1.0e-9)
    offset = curved_offset_loop(narrow, -collar, runs, TOLERANCE)
    @test offset.bridged && count(item -> item.collapsed, offset.items) == 1
    # A bridged offset is never a simple offset: the non-collar callers
    # (offset_loop_points: retained-mask / metal lofts, rounding primitives) fail
    # closed; only collar_loop_points routes it to the union.
    message = guard_message(() -> offset_loop_points(narrow, -collar, TOLERANCE))
    @test occursin("bridges a collapsed circular arc", message) &&
          occursin("not a simple offset", message)
    @test offset_loop_points(narrow, -0.3, TOLERANCE) ==
          curved_offset_loop(narrow, -0.3, runs, TOLERANCE).points
    absorbed = Dict{String, Any}[]
    rule = collar_island_rule(narrow, -collar, box, radius, TOLERANCE)
    points, construction = collar_loop_points(
        narrow,
        -collar,
        box,
        TOLERANCE;
        island_rule=rule,
        absorbed=absorbed
    )
    @test construction == "CollarUnion" && isempty(absorbed)
    simplified, _ = simplify_footprint_polygon(points, FOOTPRINT_COLLINEAR_TOLERANCE)
    # The channel is etched throughout; the bottom sides' collar reaches y = -5.7.
    @test simplified == [(-10.0, -collar), (10.0, -collar), (10.0, 10.0), (-10.0, 10.0)]
    # Every channel point is inside the union; the metal is too.
    @test all(
        point_in_polygon(p, points, TOLERANCE) for
        p in ((0.0, 3.0), (0.9, 7.5), (0.0, 7.99))
    )
    # A channel wider than twice the collar keeps a slot: the semicircle (radius 6 >
    # collar) shrinks concentrically to 0.3 about its centre (0, 2), the walls shift to
    # x = -/+ 0.3 and the miter polygon is simple.
    wide = channelled(12.0)
    wide_points, wide_construction = collar_loop_points(wide, -collar, box, TOLERANCE)
    @test wide_construction == "MiterOffset"
    run = only(circular_arc_runs(wide_points, TOLERANCE))
    @test isapprox(run.radius, 6.0 - collar; atol=1.0e-9) &&
          hypot(run.center[1], run.center[2] - 2.0) <= 1.0e-9
    @test any(
        isapprox(p[1], -0.3; atol=1.0e-9) && isapprox(p[2], -collar; atol=1.0e-9) for
        p in wide_points
    )
    @test any(
        isapprox(p[1], 0.3; atol=1.0e-9) && isapprox(p[2], 2.0; atol=1.0e-9) for
        p in wide_points
    )
end

@testset "straight loops are bitwise unchanged (reference implementation)" begin
    notched(width) = [
        (-10.0, 0.0),
        (-width / 2, 0.0),
        (-width / 2, 6.0),
        (width / 2, 6.0),
        (width / 2, 0.0),
        (10.0, 0.0),
        (10.0, 10.0),
        (-10.0, 10.0)
    ]
    classes = vcat(fill("Physical", 5), fill("Continuation", 3))
    l_shape = [(0.0, 0.0), (1.2, 0.0), (1.2, 0.3), (0.3, 0.3), (0.3, 1.0), (0.0, 1.0)]
    keyhole = [
        (-20.0, 0.0),
        (-4.0, 0.0),
        (-4.0, 4.0),
        (-16.0, 4.0),
        (-16.0, 18.0),
        (16.0, 18.0),
        (16.0, 4.0),
        (4.0, 4.0),
        (4.0, 0.0),
        (20.0, 0.0),
        (20.0, 20.0),
        (-20.0, 20.0)
    ]
    taper = [(0.0, 0.0), (5.0, -0.7), (9.0, 0.4), (10.0, 3.0), (6.0, 4.5), (1.0, 3.9)]
    loops = [
        exterior(notched(8.0); classes=classes),
        exterior(notched(16.0); classes=classes),
        exterior(l_shape),
        exterior(keyhole; classes=vcat(fill("Physical", 10), fill("Continuation", 2))),
        exterior(taper)
    ]
    for loop in loops, distance in (-6.0, -1.5, -0.05, 0.02, 0.1)
        @test offset_loop_points(loop, distance, TOLERANCE) ==
              reference_offset_loop_points(loop, distance, TOLERANCE)
    end
    hexagon = [(2.0 * cos(a), 2.0 * sin(a)) for a in range(0.0, 2pi; length=7)[1:6]]
    holes = [
        hole_loop(hexagon),
        hole_loop([(0.0, 0.0), (3.0, 0.0), (3.0, 1.0), (0.0, 1.0)]),
        hole_loop(taper)
    ]
    for hole in holes, distance in (-0.4, -0.2, -0.05, 0.1, 0.3)
        @test offset_hole_points(hole, distance, TOLERANCE) ==
              reference_offset_hole_points(hole, distance, TOLERANCE)
    end
    # The collar union of the straight notched loop (decision 54a) is untouched.
    box = ([-10.0, -10.0], [10.0, 10.0])
    points, construction = collar_loop_points(loops[1], -6.0, box, TOLERANCE)
    @test construction == "CollarUnion"
    simplified, _ = simplify_footprint_polygon(points, FOOTPRINT_COLLINEAR_TOLERANCE)
    @test simplified == [(-10.0, -6.0), (10.0, -6.0), (10.0, 10.0), (-10.0, 10.0)]
    @test length(
        collar_pieces(
            loops[1],
            -6.0,
            offset_loop_points(loops[1], -6.0, TOLERANCE),
            box...,
            TOLERANCE
        )
    ) == 1 + 5 + 2
end

@testset "island rule (decision 246 B): the loop end's island is absorbed and recorded" begin
    loop = only(
        read_boundary(
            joinpath(
                @__DIR__,
                "testdata",
                "loop-end-5ed91f8890c0",
                "plan-view-boundary.csv"
            )
        )
    )
    radius = 1.9
    tolerance = 1.0e-7 * radius
    xs = first.(loop.points);
    ys = last.(loop.points)
    box = ((minimum(xs), minimum(ys), 0.0), (maximum(xs), maximum(ys), 1.0))
    runs = circular_arc_runs(loop.points, tolerance)
    @test length(runs) == 16
    offset = curved_offset_loop(loop, -3radius, runs, tolerance)
    @test !offset.bridged && count(item -> item.collapsed, offset.items) == 10
    @test count(item -> item.kind == :arc && !item.collapsed, offset.items) == 6
    @test all(
        isapprox(item.radius, item.original_radius + 3radius; atol=1.0e-9) for
        item in offset.items if item.kind == :arc && !item.collapsed
    )
    # Without the rule the island fails closed as before.
    message = guard_message(() -> collar_loop_points(loop, -3radius, box, tolerance))
    @test occursin("ScopeGuard[FootprintTopology]", message) && occursin("island", message)
    absorbed = Dict{String, Any}[]
    rule = collar_island_rule(loop, -3radius, box, radius, tolerance)
    @test rule.cap == 0.05 * radius && rule.collar == 3radius
    points, construction = collar_loop_points(
        loop,
        -3radius,
        box,
        tolerance;
        island_rule=rule,
        absorbed=absorbed
    )
    @test construction == "CollarUnion"
    simplified, _ = simplify_footprint_polygon(points, FOOTPRINT_COLLINEAR_TOLERANCE)
    @test simplified == [
        (box[1][1], box[1][2]),
        (box[2][1], box[1][2]),
        (box[2][1], box[2][2]),
        (box[1][1], box[2][2])
    ]
    @test length(absorbed) == 1
    island = absorbed[1]
    @test island["Vertices"] >= 3
    @test 0.005 <= island["Area"] <= 0.02
    @test 0.03 <= island["MaximumExcess"] <= island["MaximumExcessBound"] <= rule.cap
    @test island["MaximumExcessBound"] <= 0.025 * radius
    @test island["ExcessCap"] == rule.cap && island["Collar"] == 3radius
    at = island["MaximumExcessPoint"]
    @test hypot(at[1] - 0.86, at[2] + 0.19) <= 0.1
    @test all(abs(p[1] - 0.86) <= 0.1 && abs(p[2] + 0.2) <= 0.3 for p in island["Points"])
    # The census record flows through the footprint record.
    record = footprint_record(
        1,
        0.0,
        false,
        simplified,
        Dict{String, Any}(),
        "CollarUnion";
        absorbed_islands=absorbed
    )
    @test record["AbsorbedIslands"] == absorbed
    @test !haskey(
        footprint_record(1, 0.0, false, simplified, Dict{String, Any}(), "MiterOffset"),
        "AbsorbedIslands"
    )
end

@testset "island rule: an island needing 0.06 R is refused, 0.03 R absorbed, a box-face notch untouched" begin
    radius = 2.0
    collar = 3radius
    box = ([-20.0, -10.0], [20.0, 20.0])
    # A keyhole: entry 8 wide, chamber (2 x collar + 2 x excess) wide and 14 tall under
    # the collar: the island is 2 x excess wide and 2 tall at the chamber's centre.
    function keyhole(excess)
        half = collar + excess
        points = [
            (-20.0, 0.0),
            (-4.0, 0.0),
            (-4.0, 4.0),
            (-half, 4.0),
            (-half, 18.0),
            (half, 18.0),
            (half, 4.0),
            (4.0, 4.0),
            (4.0, 0.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (-20.0, 20.0)
        ]
        return exterior(points; classes=vcat(fill("Physical", 10), fill("Continuation", 2)))
    end
    for excess in (0.03 * radius, 0.0499 * radius)
        loop = keyhole(excess)
        absorbed = Dict{String, Any}[]
        points, construction = collar_loop_points(
            loop,
            -collar,
            box,
            TOLERANCE;
            island_rule=collar_island_rule(loop, -collar, box, radius, TOLERANCE),
            absorbed=absorbed
        )
        @test construction == "CollarUnion" && length(absorbed) == 1
        island = absorbed[1]
        @test isapprox(island["Area"], 2excess * 2.0; atol=1.0e-9)
        @test isapprox(island["MaximumExcessBound"], excess; atol=1.0e-9)
        @test isapprox(island["MaximumExcess"], excess; atol=1.0e-3 * excess)
        simplified, _ = simplify_footprint_polygon(points, FOOTPRINT_COLLINEAR_TOLERANCE)
        @test simplified == [(-20.0, -6.0), (20.0, -6.0), (20.0, 20.0), (-20.0, 20.0)]
    end
    refused = keyhole(0.06 * radius)
    message = guard_message(
        () -> collar_loop_points(
            refused,
            -collar,
            box,
            TOLERANCE;
            island_rule=collar_island_rule(refused, -collar, box, radius, TOLERANCE),
            absorbed=Dict{String, Any}[]
        )
    )
    @test occursin("ScopeGuard[FootprintTopology]", message) &&
          occursin("above the admitted", message)
    # The un-etched notch at the box face (decision 229): two leads separated by a gap
    # of 2 x collar + 2 x 0.03 R opening at the top box face stay separated by the
    # notch, which is part of the outer boundary, not an island; a 2-um slot in the
    # bottom side routes the collar to the union construction.
    excess = 0.03 * radius
    half = collar + excess
    notch = exterior(
        [
            (-20.0, 0.0),
            (-1.0, 0.0),
            (-1.0, 3.0),
            (1.0, 3.0),
            (1.0, 0.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (half, 20.0),
            (half, 8.0),
            (-half, 8.0),
            (-half, 20.0),
            (-20.0, 20.0)
        ];
        classes=vcat(
            fill("Physical", 6),
            [
                "Continuation",
                "Physical",
                "Physical",
                "Physical",
                "Continuation",
                "Continuation"
            ]
        )
    )
    @test !polygon_is_simple(offset_loop_points(notch, -collar, TOLERANCE), TOLERANCE)
    absorbed = Dict{String, Any}[]
    points, construction = collar_loop_points(
        notch,
        -collar,
        box,
        TOLERANCE;
        island_rule=collar_island_rule(notch, -collar, box, radius, TOLERANCE),
        absorbed=absorbed
    )
    @test construction == "CollarUnion" && isempty(absorbed)
    simplified, _ = simplify_footprint_polygon(points, FOOTPRINT_COLLINEAR_TOLERANCE)
    @test Set(round.(p; digits=9) for p in simplified) == Set([
        (-20.0, -6.0),
        (20.0, -6.0),
        (20.0, 20.0),
        (excess, 20.0),
        (excess, 14.0),
        (-excess, 14.0),
        (-excess, 20.0),
        (-20.0, 20.0)
    ])
    @test !point_in_polygon((0.0, 17.0), points, TOLERANCE)
    @test point_in_polygon((0.0, 1.5), points, TOLERANCE)
    # The same notch without the slot is a simple miter polygon: the notch stays too.
    plain = exterior(
        [
            (-20.0, 0.0),
            (20.0, 0.0),
            (20.0, 20.0),
            (half, 20.0),
            (half, 8.0),
            (-half, 8.0),
            (-half, 20.0),
            (-20.0, 20.0)
        ];
        classes=[
            "Physical",
            "Continuation",
            "Physical",
            "Physical",
            "Physical",
            "Continuation",
            "Continuation",
            "Physical"
        ]
    )
    plain_points, plain_construction = collar_loop_points(
        plain,
        -collar,
        box,
        TOLERANCE;
        island_rule=collar_island_rule(plain, -collar, box, radius, TOLERANCE),
        absorbed=absorbed
    )
    @test plain_construction == "MiterOffset" && isempty(absorbed)
    @test !point_in_polygon((0.0, 17.0), plain_points, TOLERANCE)
end

@testset "island rule: a kite-bounded island fails closed when the measured excess exceeds the cap" begin
    radius = 2.0
    collar = 3radius
    excess = 0.03 * radius
    box = ([-20.0, -10.0], [20.0, 30.0])
    # A stepped keyhole: entry 8 wide, a narrow chamber (2 x collar + 2 x excess)
    # wide from y = 4 to 11, then a wide chamber (16 wide) up to the top wall y = 22.
    # The step corners (-/+ (collar + excess), 11) are convex metal corners whose
    # miter kites (squares of side collar) bound the island laterally above y = 11:
    # the island is x in -/+ excess, y in (10, 16), bound = excess (half its width),
    # but its points are farther than the collar from every metal point (the step
    # corners, the top wall, the far walls): the exact excess reaches ~1.17 = 0.58 R.
    half = collar + excess
    stepped = exterior(
        [
            (-20.0, 0.0),
            (-4.0, 0.0),
            (-4.0, 4.0),
            (-half, 4.0),
            (-half, 11.0),
            (-8.0, 11.0),
            (-8.0, 22.0),
            (8.0, 22.0),
            (8.0, 11.0),
            (half, 11.0),
            (half, 4.0),
            (4.0, 4.0),
            (4.0, 0.0),
            (20.0, 0.0),
            (20.0, 30.0),
            (-20.0, 30.0)
        ];
        classes=vcat(fill("Physical", 13), fill("Continuation", 3))
    )
    rule = collar_island_rule(stepped, -collar, box, radius, TOLERANCE)
    island = [(-excess, 10.0), (excess, 10.0), (excess, 16.0), (-excess, 16.0)]
    @test isapprox(island_excess_bound(island), excess; atol=1.0e-12)
    @test island_excess_bound(island) <= rule.cap
    measured, at = island_maximum_excess(island, rule, TOLERANCE)
    @test 1.1 <= measured <= 1.2 && measured > rule.cap
    @test abs(at[1]) <= 1.0e-3 && 14.5 <= at[2] <= 15.2
    # The bound alone would admit the island; the gate also needs the measurement.
    message = guard_message(
        () -> collar_loop_points(
            stepped,
            -collar,
            box,
            TOLERANCE;
            island_rule=rule,
            absorbed=Dict{String, Any}[]
        )
    )
    @test occursin("ScopeGuard[FootprintTopology]", message) &&
          occursin("above the admitted", message) &&
          occursin("the measurement exceeds the bound", message)
    # The same stepped chamber with its top wall at y = 16.1 leaves the island
    # x in -/+ excess, y in (10, 10.1), still bounded below by the entry corners'
    # kites, but its measured excess (the narrow walls, excess; the top wall, 0.05)
    # agrees with the cap: absorbed and recorded with both measures within the cap.
    low = exterior(
        [p[2] == 22.0 ? (p[1], 16.1) : p for p in stepped.points];
        classes=stepped.classes
    )
    absorbed = Dict{String, Any}[]
    points, construction = collar_loop_points(
        low,
        -collar,
        box,
        TOLERANCE;
        island_rule=collar_island_rule(low, -collar, box, radius, TOLERANCE),
        absorbed=absorbed
    )
    @test construction == "CollarUnion" && length(absorbed) == 1
    record = absorbed[1]
    @test isapprox(record["Area"], 2excess * 0.1; atol=1.0e-9)
    @test isapprox(record["MaximumExcessBound"], 0.05; atol=1.0e-9)
    @test isapprox(record["MaximumExcess"], excess; atol=1.0e-3 * excess)
    @test record["MaximumExcess"] <= record["ExcessCap"] &&
          record["MaximumExcessBound"] <= record["ExcessCap"]
end

@testset "hole shrink of a rounded-rectangle hole: concentric arcs, collapse, vanishing" begin
    r = 0.5
    hole = hole_loop(rounded_rectangle(-3.0, -1.0, 3.0, 1.0, r))
    shrunk = offset_hole_points(hole, -0.2, TOLERANCE)
    runs = circular_arc_runs(shrunk, TOLERANCE)
    @test length(runs) == 4 &&
          all(isapprox(run.radius, r - 0.2; atol=1.0e-9) for run in runs)
    @test all(abs(p[1]) <= 2.8 + 1.0e-9 && abs(p[2]) <= 0.8 + 1.0e-9 for p in shrunk)
    @test any(isapprox(p[2], 0.8; atol=1.0e-9) for p in shrunk)
    # Collar >= r: the arcs collapse onto the straight sides' corners.
    sharp = offset_hole_points(hole, -0.6, TOLERANCE)
    @test length(sharp) == 4 &&
          Set(round.(p; digits=9) for p in sharp) ==
          Set([(-2.4, -0.4), (2.4, -0.4), (2.4, 0.4), (-2.4, 0.4)])
    @test isempty(offset_hole_points(hole, -1.0, TOLERANCE))
    @test isempty(offset_hole_points(hole, -1.2, TOLERANCE))
    grown = offset_hole_points(hole, 0.3, TOLERANCE)
    @test all(
        isapprox(run.radius, r + 0.3; atol=1.0e-9) for
        run in circular_arc_runs(grown, TOLERANCE)
    )
end

@testset "round 3 class (5) interim 5B (decision 510; part M 2.2): a concave arc leaving the box by less than the collar fails closed at ScopeGuard[CollarFaceEnd]" begin
    # A metal slab above the box face y = 0 with a concave circular bite (the metal OUTSIDE the
    # circle of radius 10 centred at (0, -c)) whose two ends lie on the face: the 5.7 collar,
    # offset away from the metal, SHRINKS the circle to 4.3; it meets the face line only when
    # the circle protrudes into the box by rho - c > 5.7 (the 32dc558f4810 mechanism: a CPW
    # ground's bend leaving the box by 2.18 um < 3 R).
    rho = 10.0
    collar = 5.7
    function bitten(c; chords=12)
        half = sqrt(rho^2 - c^2)
        a0 = atan(c, -half)                 # the left face intersection, over the top, to the right
        a1 = atan(c, half)
        arc = [(rho * cos(a), -c + rho * sin(a)) for a in range(a0, a1; length=chords + 1)]
        arc[1] = (-half, 0.0)
        arc[end] = (half, 0.0)
        points = vcat([(-20.0, 0.0)], arc, [(20.0, 0.0), (20.0, 20.0), (-20.0, 20.0)])
        n = length(points)
        classes = [i == 1 || i >= n - 3 ? "Continuation" : "Physical" for i = 1:n]
        return exterior(points; classes=classes)
    end
    for c in (6.0, 4.4)                     # protrusions 4.0 and 5.6 < 5.7: the guard fires
        loop = bitten(c)
        runs = circular_arc_runs(loop.points, TOLERANCE)
        @test length(runs) == 1 && isapprox(runs[1].radius, rho; atol=1.0e-9)
        message = guard_message(() -> curved_offset_loop(loop, -collar, runs, TOLERANCE))
        @test occursin("ScopeGuard[CollarFaceEnd]", message) &&
              occursin("box-face line", message) &&
              occursin("collar width $collar", message)
        @test isapprox(
            parse(Float64, match(r"line ([0-9.e+-]+) from its centre", message)[1]),
            c;
            atol=1.0e-9
        )
        @test isapprox(
            parse(Float64, match(r"rho - h_face = ([0-9.e+-]+)", message)[1]),
            rho - c;
            atol=1.0e-9
        )
    end
    # Protrusion 7 > 5.7: the shrunk circle meets the face line; both junctions resolve ON the
    # concentric offset circle of radius rho - collar.
    loop = bitten(3.0)
    runs = circular_arc_runs(loop.points, TOLERANCE)
    offset = curved_offset_loop(loop, -collar, runs, TOLERANCE)
    @test !offset.bridged && all(j.point !== nothing for j in offset.junctions)
    shrunk = only(item for item in offset.items if item.kind == :arc)
    @test isapprox(shrunk.radius, rho - collar; atol=1.0e-9) &&
          isapprox(shrunk.center[2], -3.0; atol=1.0e-9)
    arc_index = findfirst(item -> item.kind == :arc, offset.items)
    arc_junctions =
        [j for j in offset.junctions if j.before == arc_index || j.after == arc_index]
    @test length(arc_junctions) == 2 && all(
        distance_to_circle(j.point, shrunk.center, shrunk.radius) <= 1.0e-9 &&
        abs(j.point[2]) <= 1.0e-9 for j in arc_junctions
    )
    # 510 MINOR-8 (b): a line-arc CORNER kink whose shifted (Physical) line misses the shrunk
    # concave circle is named by the same guard as the untested corner case: a 1-wide slot
    # from the face into a circular chamber of radius 6 (its walls shifted by 5.7 miss the
    # 0.3 circle).
    chamber = 6.0
    y_wall = 10.0 - sqrt(chamber^2 - 0.25)
    a_left = atan(y_wall - 10.0, -0.5)
    a_right = atan(y_wall - 10.0, 0.5)
    sweep = mod(a_left - a_right, 2pi)       # clockwise the long way round, over the top
    chords = 24
    bite = [
        (
            chamber * cos(a_left - sweep * k / chords),
            10.0 + chamber * sin(a_left - sweep * k / chords)
        ) for k = 0:chords
    ]
    bite[1] = (-0.5, y_wall)
    bite[end] = (0.5, y_wall)
    slot_points = vcat(
        [(-20.0, 0.0), (-0.5, 0.0)],
        bite,
        [(0.5, 0.0), (20.0, 0.0), (20.0, 30.0), (-20.0, 30.0)]
    )
    ns = length(slot_points)
    slot_classes = [i == 1 || i >= ns - 3 ? "Continuation" : "Physical" for i = 1:ns]
    slot = exterior(slot_points; classes=slot_classes)
    slot_runs = circular_arc_runs(slot.points, TOLERANCE)
    @test length(slot_runs) == 1 && isapprox(slot_runs[1].radius, chamber; atol=1.0e-9)
    message = guard_message(() -> curved_offset_loop(slot, -collar, slot_runs, TOLERANCE))
    @test occursin("ScopeGuard[CollarFaceEnd]", message) &&
          occursin("corner kink whose offsets do not meet (untested)", message)
end
