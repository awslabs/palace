# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Our own structured boundary-layer band (`band_mode = :own`, supervisor decisions 196 / 198;
# design: structured-band-20261002/DESIGN.md). Included by PolygonWindowMesh.jl.
#
# Gmsh's BoundaryLayer field takes one Thickness per curve set and terminates a layer only
# at a corner, so the per-point band cap of decision 191 (band <= 0.4 x the local facing
# distance, first rows unchanged) cannot be realised with it, and its extrusion fails where two
# fronts overlap (the S1p junction region). Here the band quads are built by this code and
# Gmsh only triangulates the remaining region of every partition:
#
#   * a column at every transfinite node of every metal curve of the partition's loops: nodes
#     p + h_k * scale * u, h_k = r (2^k - 1), k = 1..n (the loop is walked with the partition
#     on the LEFT: outer loops counter-clockwise, holes clockwise; u points inward);
#   * rows per 1D SEGMENT from the facing distance of that segment (the exact distance to every
#     other plan curve clipped to the half-plane on the band side; curves sharing an endpoint
#     and collinear continuations excluded; other bodies' corners count), applied to both end
#     columns of the segment; a SCALED column (mitre, inward, wall end) is further capped by
#     its LENGTH: h_k * scale <= 0.4 x the smaller facing distance of its two adjacent
#     segments (review M1: the mitre tip reaches sqrt 2 h along the diagonal, so two convex
#     corners facing each other diagonally at D would otherwise reach 1.13 D together); where
#     neighbouring columns differ in rows the excess rows of the taller column collapse to
#     triangles on the shorter column's top node (the band top is one edge per node pair, Gmsh
#     never sees the step);
#   * corners: a right turn (OUTWARD: the band wraps around a convex metal corner) up to the
#     fan threshold and every left turn (INWARD: a concave metal corner, a T / X junction of
#     the cross-plane snap) get the bisector column scaled by 1 / cos(turn / 2), so the row-k
#     node is exactly h_k from both edge lines — the recorded transmon generator listed no
#     FanPointsList and its first-layer heights (0.995-1.005 r) are consistent with this
#     column; a right turn above the threshold (compared with the angular tolerance
#     FAN_TURN_TOLERANCE_DEG, so an exact right angle and one whose 1e-9-rounded coordinates
#     turn by 90.0000000x deg are both mitred) gets a fan of FAN_COLUMNS columns; at an inward
#     corner of interior angle theta every column at arc distance s along the two adjacent
#     curves (the corner column: s = the smaller adjacent spacing) keeps h_k <= 0.5 s tan(theta
#     / 2) so neighbouring columns never cross (the HANDOFF's half-distance rule);
#   * a metal curve ending on the window wall: the end column runs ALONG the wall into the
#     partition (scale 1 / sin(theta), theta the partition's interior angle at the wall), the
#     wall as the second edge of the clearance rule;
#   * backstop: every band quad strictly convex, every band triangle positive, no two band
#     edges of the partition intersecting (grid-accelerated exact test), every band node
#     inside the box, every end column shorter than its wall curve — any failure is a
#     deterministic refusal naming the partition and the location;
#   * the remaining region (band tops + wall curves + end columns) is a fresh built-in-kernel
#     Gmsh model meshed with Algorithm 6 (transfinite 1 segment per band edge), its triangles
#     joined to the band triangles; band + region area = the partition's OCC mass within
#     REGION_AREA_TOLERANCE (relative).

const FAN_COLUMNS = 7
const DEFAULT_FAN_TURN_ANGLE_DEG = 90.0
# A turn fans only if it exceeds the threshold by more than this tolerance: right angles whose
# coordinates carry the writer's 1e-9 rounding turn by 90.0000000x deg and must stay mitres.
const FAN_TURN_TOLERANCE_DEG = 1.0e-6
# Above this turn the mitre column 1 / cos(turn / 2) would exceed 2 x the band: must be a fan.
# The limit bounds OUTWARD corners only; inward mitres (1 / cos(turn / 2)) and wall-end
# columns (1 / sin theta) scale without limit (geometrically forced), their LENGTH is capped
# by the facing distance and the maximum scale is recorded in the manifest.
const MAX_FAN_TURN_ANGLE_DEG = 120.0
# Band + region area versus the partition's OCC mass. The construction closes to ~1e-11; the
# transmon's exported box (1.5e-7 um outside its ground outline, merged by OCC) leaves 8e-11.
const REGION_AREA_TOLERANCE = 1.0e-9
const CORNER_CLEARANCE_FRACTION = 0.5
const COLLINEAR_TOLERANCE_UM = 1.0e-7
const HALF_PLANE_TOLERANCE_UM = 1.0e-9
const BAND_NODE_KEY_UM = 1.0e-7

struct PlanCurve
    tag::Int32
    a::Point2
    b::Point2
    points::NTuple{2, Int32}
    metal::Bool
end

function point_coordinates(tag::Int32)
    xyz = gmsh.model.get_value(0, tag, Float64[])
    return (xyz[1], xyz[2])
end

# Every straight curve of the fragmented plan with its end point tags and the metal flag.
function plan_curves(surfaces, metal_edge_curves)
    curves = PlanCurve[]
    for (_, tag) in unique(gmsh.model.get_boundary(surfaces, false, false, false))
        tag = abs(tag)
        gmsh.model.get_type(1, tag) == "Line" ||
            error("Plan curve $tag is not a straight line")
        points =
            [abs(p) for (_, p) in gmsh.model.get_boundary([(1, tag)], false, false, false)]
        length(points) == 2 || error("Plan curve $tag has $(length(points)) end points")
        push!(
            curves,
            PlanCurve(
                tag,
                point_coordinates(points[1]),
                point_coordinates(points[2]),
                (points[1], points[2]),
                tag in metal_edge_curves
            )
        )
    end
    sort!(curves; by=c -> c.tag)
    return curves
end

# Transfinite nodes of a metal curve, from points[1] to points[2]: m = ceil(L / t - 1e-6)
# segments (the relative guard keeps an exact multiple of t from gaining a segment by
# rounding; identical to the exporter's chord count on every recorded edge).
function curve_transfinite_nodes(curve::PlanCurve, tangential_um::Float64)
    length_um = point_distance(curve.a, curve.b)
    m = max(1, ceil(Int, length_um / tangential_um - 1.0e-6))
    return [
        (
            curve.a[1] + (curve.b[1] - curve.a[1]) * i / m,
            curve.a[2] + (curve.b[2] - curve.a[2]) * i / m
        ) for i = 0:m
    ]
end

# A loop of a partition: vertex point tags and the curve (index into the plan curves) joining
# vertex i to vertex i + 1, walked with the partition's interior on the left.
struct PartitionLoop
    vertices::Vector{Int32}
    vertex_xy::Vector{Point2}
    curves::Vector{Int}
end

struct PartitionGeometry
    surface::Int32
    area::Float64
    loops::Vector{PartitionLoop} # the outer loop first
end

function partition_geometry(surface::Int32, curves::Vector{PlanCurve}, box)
    index_of = Dict(curve.tag => k for (k, curve) in enumerate(curves))
    members = [
        index_of[abs(tag)] for
        (_, tag) in unique(gmsh.model.get_boundary([(2, surface)], false, false, false))
    ]
    by_point = Dict{Int32, Vector{Int}}()
    for k in members, p in curves[k].points
        push!(get!(by_point, p, Int[]), k)
    end
    for (p, incident) in by_point
        length(incident) == 2 || error(
            "Partition $surface touches itself at point $p ($(point_coordinates(p))): " *
            "$(length(incident)) incident boundary curves"
        )
    end
    for k in members
        curves[k].metal ||
            on_box_wall(box, gmsh.model.get_bounding_box(1, curves[k].tag)) ||
            error(
                "Non-metal plan curve $(curves[k].tag) of partition $surface is off the wall"
            )
    end
    visited = Set{Int}()
    loops = PartitionLoop[]
    for start in members
        start in visited && continue
        vertices = Int32[]
        loop_curves = Int[]
        current, point = start, curves[start].points[1]
        while true
            push!(visited, current)
            push!(loop_curves, current)
            push!(vertices, point)
            c = curves[current]
            point = c.points[1] == point ? c.points[2] : c.points[1]
            next = only(filter(!=(current), by_point[point]))
            next == start && break
            current = next
        end
        push!(loops, PartitionLoop(vertices, point_coordinates.(vertices), loop_curves))
    end
    areas = [polygon_signed_area(loop.vertex_xy) for loop in loops]
    outer = argmax(abs.(areas))
    oriented = PartitionLoop[]
    for (k, loop) in enumerate(loops)
        want_ccw = k == outer
        if (areas[k] > 0.0) != want_ccw
            loop = PartitionLoop(
                reverse(loop.vertices),
                reverse(loop.vertex_xy),
                circshift(reverse(loop.curves), -1)
            )
        end
        k == outer ? pushfirst!(oriented, loop) : push!(oriented, loop)
    end
    return PartitionGeometry(surface, gmsh.model.occ.get_mass(2, surface), oriented)
end

# ---------------------------------------------------------------------------------------------
# Facing distance per 1D segment and side (review M1 / M2, per side m8).

function orient(a::Point2, b::Point2, c::Point2)
    return (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
end

function segments_cross(a::Point2, b::Point2, c::Point2, d::Point2)
    d1, d2 = orient(c, d, a), orient(c, d, b)
    d3, d4 = orient(a, b, c), orient(a, b, d)
    return d1 * d2 < 0.0 && d3 * d4 < 0.0
end

function segment_segment_distance(a::Point2, b::Point2, c::Point2, d::Point2)
    segments_cross(a, b, c, d) && return 0.0
    return min(
        segment_projection(a, c, d)[1],
        segment_projection(b, c, d)[1],
        segment_projection(c, a, b)[1],
        segment_projection(d, a, b)[1]
    )
end

# The part of segment cd on the open side {x : (x - origin) . normal > tolerance}, or nothing.
function clip_to_half_plane(c::Point2, d::Point2, origin::Point2, normal::Point2)
    sc = (c[1] - origin[1]) * normal[1] + (c[2] - origin[2]) * normal[2]
    sd = (d[1] - origin[1]) * normal[1] + (d[2] - origin[2]) * normal[2]
    tolerance = HALF_PLANE_TOLERANCE_UM
    sc <= tolerance && sd <= tolerance && return nothing
    sc > tolerance && sd > tolerance && return (c, d)
    t = (tolerance - sc) / (sd - sc)
    x = (c[1] + t * (d[1] - c[1]), c[2] + t * (d[2] - c[2]))
    return sc > tolerance ? (c, x) : (x, d)
end

function line_distance(p::Point2, a::Point2, b::Point2)
    return abs(orient(a, b, p)) / point_distance(a, b)
end

"""
    segment_facing_distance(q1, q2, normal, curve, curves, reach) -> distance

Exact distance from the 1D segment `q1 q2` of `curve` to the nearest other plan curve (metal
of every plane and bump, or box wall) lying on the band side (the open half-plane of `normal`
through the segment). Curves sharing an end point with `curve` (corner-adjacent) and curves
collinear with it (continuations of the run) are excluded; everything else counts in full,
including other bodies' corners. Curves farther than `reach` (bounding boxes) are skipped.
"""
function segment_facing_distance(
    q1::Point2,
    q2::Point2,
    normal::Point2,
    curve::PlanCurve,
    curves::Vector{PlanCurve},
    reach::Float64
)
    best = Inf
    xmin, xmax = minmax(q1[1], q2[1])
    ymin, ymax = minmax(q1[2], q2[2])
    for other in curves
        other.tag == curve.tag && continue
        (other.points[1] in curve.points || other.points[2] in curve.points) && continue
        min(other.a[1], other.b[1]) > xmax + reach && continue
        max(other.a[1], other.b[1]) < xmin - reach && continue
        min(other.a[2], other.b[2]) > ymax + reach && continue
        max(other.a[2], other.b[2]) < ymin - reach && continue
        line_distance(other.a, curve.a, curve.b) <= COLLINEAR_TOLERANCE_UM &&
            line_distance(other.b, curve.a, curve.b) <= COLLINEAR_TOLERANCE_UM &&
            continue
        clipped = clip_to_half_plane(other.a, other.b, q1, normal)
        clipped === nothing && continue
        best = min(best, segment_segment_distance(q1, q2, clipped[1], clipped[2]))
    end
    return best
end

band_heights(radial_um::Float64, growth::Float64, rows::Int) =
    [radial_um * (growth^k - 1.0) / (growth - 1.0) for k = 1:rows]

# Largest row count whose height stays within `limit` (at least one row: the first row is
# never dropped; the caller counts the clamped columns and the backstop decides).
function rows_within(heights::Vector{Float64}, limit::Float64)
    rows = count(<=(limit), heights)
    return max(1, rows), rows == 0
end

# ---------------------------------------------------------------------------------------------
# Band construction per partition.

struct BandColumn
    base::Point2
    nodes::Vector{Point2}
end

mutable struct BandStatistics
    quads::Int
    collapse_triangles::Int
    fan_triangles::Int
    fans::Int
    mitre_corners::Int
    inward_corners::Int
    clearance_capped_columns::Int
    clearance_clamped_columns::Int
    scale_capped_columns::Int # scaled columns shortened by the length cap (M1)
    scale_clamped_columns::Int # scaled columns whose first row already exceeds the cap
    wall_end_columns::Int
    segments_below_2p5r::Int
    quad_min_abs_sin::Float64
    triangle_min_angle_deg::Float64
    max_mitre_scale::Float64 # outward (mitre) corners: 1 / cos(turn / 2)
    max_inward_scale::Float64 # inward corners: 1 / cos(turn / 2)
    max_wall_end_scale::Float64 # wall ends: 1 / sin(theta)
    length_um_by_rows::Dict{Int, Float64}
end
BandStatistics() = BandStatistics(
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    1.0,
    180.0,
    1.0,
    1.0,
    1.0,
    Dict{Int, Float64}()
)

unit(v::Point2) = (h=hypot(v[1], v[2]); (v[1] / h, v[2] / h))
left_normal(t::Point2) = (-t[2], t[1])
rotate(v::Point2, angle::Float64) =
    (v[1] * cos(angle) - v[2] * sin(angle), v[1] * sin(angle) + v[2] * cos(angle))

# One metal chain of a loop: its base nodes in walking order (a closed chain repeats none),
# the curve of every segment, and the wall directions at the ends of an open chain.
struct MetalChain
    closed::Bool
    bases::Vector{Point2}
    segment_curve::Vector{Int} # curve index per segment (bases[i] -> bases[i + 1])
    corner::Vector{Bool} # base i is a joint between two curves
    wall_start::Union{Nothing, Point2} # unit direction along the wall into the partition
    wall_end::Union{Nothing, Point2}
    wall_start_vertex::Int32 # vertex tags of the chain ends (0 for a closed chain)
    wall_end_vertex::Int32
end

function loop_metal_chains(
    loop::PartitionLoop,
    curves::Vector{PlanCurve},
    nodes_by_curve::Dict{Int, Vector{Point2}}
)
    m = length(loop.curves)
    vertex_xy = loop.vertex_xy
    walked_nodes(i) = begin
        k = loop.curves[i]
        nodes = nodes_by_curve[k]
        curves[k].points[1] == loop.vertices[i] ? nodes : reverse(nodes)
    end
    is_metal = [curves[k].metal for k in loop.curves]
    chains = MetalChain[]
    if all(is_metal)
        bases = Point2[]
        segment_curve = Int[]
        corner = Bool[]
        for i = 1:m
            nodes = walked_nodes(i)
            for j = 1:(length(nodes) - 1)
                push!(bases, nodes[j])
                push!(segment_curve, loop.curves[i])
                push!(corner, j == 1)
            end
        end
        push!(
            chains,
            MetalChain(true, bases, segment_curve, corner, nothing, nothing, 0, 0)
        )
        return chains
    end
    any(is_metal) || return chains
    # Start right after a wall curve so every open chain is walked whole.
    first_wall = findfirst(!, is_metal)
    order = [mod1(first_wall + s, m) for s = 1:m]
    i = 1
    while i <= m
        if !is_metal[order[i]]
            i += 1
            continue
        end
        start = i
        while i <= m && is_metal[order[i]]
            i += 1
        end
        run = order[start:(i - 1)]
        bases = Point2[]
        segment_curve = Int[]
        corner = Bool[]
        for (r, c) in enumerate(run)
            nodes = walked_nodes(c)
            for j = 1:(length(nodes) - 1)
                push!(bases, nodes[j])
                push!(segment_curve, loop.curves[c])
                push!(corner, j == 1 && r > 1)
            end
            r == length(run) && push!(bases, nodes[end])
        end
        push!(corner, false)
        # Wall directions: at the start back along the preceding wall curve, at the end
        # forward along the following wall curve.
        start_vertex = run[1]
        previous_vertex = mod1(start_vertex - 1, m)
        end_vertex = mod1(run[end] + 1, m)
        next_vertex = mod1(end_vertex + 1, m)
        wall_start = unit((
            vertex_xy[previous_vertex][1] - vertex_xy[start_vertex][1],
            vertex_xy[previous_vertex][2] - vertex_xy[start_vertex][2]
        ))
        wall_end = unit((
            vertex_xy[next_vertex][1] - vertex_xy[end_vertex][1],
            vertex_xy[next_vertex][2] - vertex_xy[end_vertex][2]
        ))
        push!(
            chains,
            MetalChain(
                false,
                bases,
                segment_curve,
                corner,
                wall_start,
                wall_end,
                loop.vertices[start_vertex],
                loop.vertices[end_vertex]
            )
        )
    end
    return chains
end

# Rows per segment of a chain from the facing rule (both end columns of a segment), and the
# facing distance of every segment (the length cap of the scaled columns).
function chain_segment_rows(
    chain::MetalChain,
    curves::Vector{PlanCurve},
    heights::Vector{Float64},
    radial_um::Float64,
    statistics::BandStatistics
)
    rows_max = length(heights)
    reach = heights[end] / BAND_CAP_FRACTION + 1.0
    n_segments = length(chain.segment_curve)
    rows = Vector{Int}(undef, n_segments)
    distances = Vector{Float64}(undef, n_segments)
    for s = 1:n_segments
        q1 = chain.bases[s]
        q2 =
            chain.closed ? chain.bases[mod1(s + 1, length(chain.bases))] :
            chain.bases[s + 1]
        normal = left_normal(unit((q2[1] - q1[1], q2[2] - q1[2])))
        d = segment_facing_distance(
            q1,
            q2,
            normal,
            curves[chain.segment_curve[s]],
            curves,
            reach
        )
        rows[s] = band_cap_rows(d, radial_um, rows_max)
        distances[s] = d
        d < 2.5radial_um && (statistics.segments_below_2p5r += 1)
        statistics.length_um_by_rows[rows[s]] =
            get(statistics.length_um_by_rows, rows[s], 0.0) + point_distance(q1, q2)
    end
    return rows, distances
end

# Columns of a chain (corners, fans, wall ends, the clearance caps, the length cap of the
# scaled columns) and the per-column rows.
function chain_columns(
    chain::MetalChain,
    segment_rows::Vector{Int},
    segment_distances::Vector{Float64},
    heights::Vector{Float64},
    fan_turn::Float64,
    statistics::BandStatistics
)
    n = length(chain.bases)
    n_segments = length(chain.segment_curve)
    segment_direction(s) = begin
        q1 = chain.bases[s]
        q2 = chain.closed ? chain.bases[mod1(s + 1, n)] : chain.bases[s + 1]
        unit((q2[1] - q1[1], q2[2] - q1[2]))
    end
    spacing(s) = begin
        q1 = chain.bases[s]
        q2 = chain.closed ? chain.bases[mod1(s + 1, n)] : chain.bases[s + 1]
        point_distance(q1, q2)
    end
    # Node rows: the minimum of the adjacent segments; node facing distance likewise (the
    # length cap of a scaled column).
    node_rows = Vector{Int}(undef, n)
    node_distance = Vector{Float64}(undef, n)
    for i = 1:n
        before = chain.closed ? mod1(i - 1, n) : i - 1
        after = chain.closed ? i : (i <= n_segments ? i : 0)
        adjacent = filter(>=(1), (before, after))
        node_rows[i] = minimum(segment_rows[s] for s in adjacent)
        node_distance[i] = minimum(segment_distances[s] for s in adjacent)
    end
    # Clearance caps along the two curves adjacent to every inward corner and wall end:
    # h_k <= 0.5 s tan(theta / 2) at arc distance s from the corner.
    function apply_clearance!(
        i_corner,
        theta,
        spacing_at_corner,
        backward_curve,
        forward_curve
    )
        limit_at(s) = CORNER_CLEARANCE_FRACTION * s * tan(theta / 2)
        cap!(i, s) = begin
            rows, clamped = rows_within(heights, limit_at(s))
            if rows < node_rows[i]
                node_rows[i] = rows
                statistics.clearance_capped_columns += 1
            end
            clamped && (statistics.clearance_clamped_columns += 1)
        end
        cap!(i_corner, spacing_at_corner)
        # Backward along the preceding curve, forward along the following curve, while the
        # bound can still bind.
        for (direction, wanted) in ((-1, backward_curve), (1, forward_curve))
            wanted === nothing && continue
            s_arc, i, steps = 0.0, i_corner, 0
            while steps < n
                segment = direction == 1 ? i : (chain.closed ? mod1(i - 1, n) : i - 1)
                (1 <= segment <= n_segments) || break
                chain.segment_curve[segment] == wanted || break
                j = chain.closed ? mod1(i + direction, n) : i + direction
                (1 <= j <= n) || break
                s_arc += spacing(segment)
                limit_at(s_arc) >= heights[end] && break
                cap!(j, s_arc)
                i, steps = j, steps + 1
            end
        end
    end
    turns = zeros(n)
    for i = 1:n
        chain.corner[i] || continue
        before = chain.closed ? mod1(i - 1, n) : i - 1
        t1, t2 = segment_direction(before), segment_direction(i)
        turns[i] = atan(t1[1] * t2[2] - t1[2] * t2[1], t1[1] * t2[1] + t1[2] * t2[2])
        if turns[i] > 0.0
            statistics.inward_corners += 1
            apply_clearance!(
                i,
                pi - turns[i],
                min(spacing(before), spacing(i)),
                chain.segment_curve[before],
                chain.segment_curve[i]
            )
        end
    end
    if !chain.closed
        for (i, wall, segment) in
            ((1, chain.wall_start, 1), (n, chain.wall_end, n_segments))
            t = segment_direction(segment)
            i == 1 || (t = (-t[1], -t[2])) # the metal direction away from the wall end
            theta = acos(clamp(t[1] * wall[1] + t[2] * wall[2], -1.0, 1.0))
            apply_clearance!(
                i,
                theta,
                spacing(segment),
                i == 1 ? nothing : chain.segment_curve[segment],
                i == 1 ? chain.segment_curve[segment] : nothing
            )
            statistics.wall_end_columns += 1
        end
    end
    # Columns. A scaled column (scale > 1: mitre, inward, wall end) is capped by its LENGTH,
    # h_k * scale <= BAND_CAP_FRACTION x the node's facing distance (review M1), so the whole
    # band over a segment stays within 0.4 x that segment's facing distance.
    columns = BandColumn[]
    column_node = Int[] # chain base index of every column
    cap_scaled!(i, scale) = begin
        scale > 1.0 || return nothing
        rows, clamped =
            rows_within(heights .* scale, BAND_CAP_FRACTION * node_distance[i])
        if rows < node_rows[i]
            node_rows[i] = rows
            statistics.scale_capped_columns += 1
        end
        clamped && (statistics.scale_clamped_columns += 1)
        return nothing
    end
    push_column!(i, u, scale) = begin
        p = chain.bases[i]
        push!(
            columns,
            BandColumn(
                p,
                [
                    (p[1] + h * scale * u[1], p[2] + h * scale * u[2]) for
                    h in heights[1:node_rows[i]]
                ]
            )
        )
        push!(column_node, i)
    end
    for i = 1:n
        if !chain.closed && (i == 1 || i == n)
            segment = i == 1 ? 1 : n_segments
            wall = i == 1 ? chain.wall_start : chain.wall_end
            t = segment_direction(segment)
            i == 1 || (t = (-t[1], -t[2]))
            sin_theta = abs(t[1] * wall[2] - t[2] * wall[1])
            sin_theta > 1.0e-9 ||
                error("Metal curve at $(chain.bases[i]) runs along the window wall")
            statistics.max_wall_end_scale =
                max(statistics.max_wall_end_scale, 1.0 / sin_theta)
            cap_scaled!(i, 1.0 / sin_theta)
            push_column!(i, wall, 1.0 / sin_theta)
            continue
        end
        before = chain.closed ? mod1(i - 1, n) : i - 1
        t1, t2 = segment_direction(before), segment_direction(i)
        n1, n2 = left_normal(t1), left_normal(t2)
        turn = turns[i]
        if turn < -(fan_turn + deg2rad(FAN_TURN_TOLERANCE_DEG))
            statistics.fans += 1
            for j = 0:(FAN_COLUMNS - 1)
                push_column!(i, rotate(n1, turn * j / (FAN_COLUMNS - 1)), 1.0)
            end
        else
            scale = 1.0 / cos(turn / 2)
            if chain.corner[i] && turn < 0.0
                statistics.mitre_corners += 1
                statistics.max_mitre_scale = max(statistics.max_mitre_scale, scale)
            elseif turn > 0.0
                statistics.max_inward_scale = max(statistics.max_inward_scale, scale)
            end
            cap_scaled!(i, scale)
            push_column!(i, unit((n1[1] + n2[1], n1[2] + n2[2])), scale)
        end
    end
    return columns, column_node, node_rows
end

function quad_min_abs_sin(a::Point2, b::Point2, c::Point2, d::Point2)
    best = 1.0
    corners = (a, b, c, d)
    for i = 1:4
        p, q, s = corners[mod1(i - 1, 4)], corners[i], corners[mod1(i + 1, 4)]
        u = unit((p[1] - q[1], p[2] - q[2]))
        v = unit((s[1] - q[1], s[2] - q[2]))
        best = min(best, abs(u[1] * v[2] - u[2] * v[1]))
    end
    return best
end

function triangle_min_angle_deg(a::Point2, b::Point2, c::Point2)
    best = 180.0
    corners = (a, b, c)
    for i = 1:3
        p, q, s = corners[mod1(i - 1, 3)], corners[i], corners[mod1(i + 1, 3)]
        u = unit((p[1] - q[1], p[2] - q[2]))
        v = unit((s[1] - q[1], s[2] - q[2]))
        best = min(best, acosd(clamp(u[1] * v[1] + u[2] * v[2], -1.0, 1.0)))
    end
    return best
end

# Elements between two consecutive columns (CCW), with the convexity / positivity checks.
function push_band_elements!(
    triangles::Vector{NTuple{3, Point2}},
    left::BandColumn,
    right::BandColumn,
    statistics::BandStatistics,
    where::String
)
    function push_triangle!(a, b, c, kind)
        orient(a, b, c) > 0.0 ||
            error("Band $kind triangle at $a is not positive ($where); the fronts collide")
        statistics.triangle_min_angle_deg =
            min(statistics.triangle_min_angle_deg, triangle_min_angle_deg(a, b, c))
        return push!(triangles, (a, b, c))
    end
    function push_quad!(a, b, c, d)
        (
            orient(a, b, c) > 0.0 &&
            orient(b, c, d) > 0.0 &&
            orient(c, d, a) > 0.0 &&
            orient(d, a, b) > 0.0
        ) || error("Band quad at $a is not strictly convex ($where); the fronts collide")
        statistics.quad_min_abs_sin =
            min(statistics.quad_min_abs_sin, quad_min_abs_sin(a, b, c, d))
        statistics.quads += 1
        return push!(triangles, (a, b, c), (a, c, d))
    end
    L, R = left.nodes, right.nodes
    if left.base == right.base
        # Fan sector: row 1 is a triangle on the corner, rows >= 2 quads.
        length(L) == length(R) || error("Fan columns with different rows ($where)")
        push_triangle!(left.base, R[1], L[1], "fan")
        statistics.fan_triangles += 1
        for k = 2:length(L)
            push_quad!(L[k - 1], R[k - 1], R[k], L[k])
        end
        return nothing
    end
    m = min(length(L), length(R))
    push_quad!(left.base, right.base, R[1], L[1])
    for k = 2:m
        push_quad!(L[k - 1], R[k - 1], R[k], L[k])
    end
    for k = (m + 1):length(L)
        push_triangle!(L[k - 1], R[m], L[k], "collapse")
        statistics.collapse_triangles += 1
    end
    for k = (m + 1):length(R)
        push_triangle!(R[k - 1], R[k], L[m], "collapse")
        statistics.collapse_triangles += 1
    end
    return nothing
end

# Exact segment intersection incl. touching (an end point on the other segment's interior).
function segments_touch(a::Point2, b::Point2, c::Point2, d::Point2)
    segments_cross(a, b, c, d) && return true
    scale = max(point_distance(a, b), point_distance(c, d))
    tolerance = 1.0e-9 * scale
    on(p, u, v) = begin
        abs(orient(u, v, p)) <= tolerance * point_distance(u, v) || return false
        t =
            ((p[1] - u[1]) * (v[1] - u[1]) + (p[2] - u[2]) * (v[2] - u[2])) /
            point_distance(u, v)^2
        return -1.0e-9 < t < 1.0 + 1.0e-9
    end
    return on(a, c, d) || on(b, c, d) || on(c, a, b) || on(d, a, b)
end

# Backstop: no two band edges of the partition (not sharing a node) intersect or touch.
function check_band_collisions(
    triangles::Vector{NTuple{3, Point2}},
    cell::Float64,
    where::String
)
    index_of = Dict{NTuple{2, Int64}, Int}()
    xy = Point2[]
    key(p) = (round(Int64, p[1] / BAND_NODE_KEY_UM), round(Int64, p[2] / BAND_NODE_KEY_UM))
    node(p) = get!(index_of, key(p)) do
        push!(xy, p)
        return length(xy)
    end
    edges = Set{Tuple{Int, Int}}()
    for (a, b, c) in triangles
        i, j, k = node(a), node(b), node(c)
        push!(edges, minmax(i, j), minmax(j, k), minmax(k, i))
    end
    edge_list = collect(edges)
    grid = Dict{Tuple{Int64, Int64}, Vector{Int}}()
    for (e, (i, j)) in enumerate(edge_list)
        p, q = xy[i], xy[j]
        for cx = floor(Int64, min(p[1], q[1]) / cell):floor(Int64, max(p[1], q[1]) / cell),
            cy = floor(Int64, min(p[2], q[2]) / cell):floor(Int64, max(p[2], q[2]) / cell)

            push!(get!(grid, (cx, cy), Int[]), e)
        end
    end
    for members in values(grid)
        for u = 1:length(members), v = (u + 1):length(members)
            e, f = members[u], members[v]
            i, j = edge_list[e]
            k, l = edge_list[f]
            (i == k || i == l || j == k || j == l) && continue
            segments_touch(xy[i], xy[j], xy[k], xy[l]) && error(
                "Band collision in $where: edge $(xy[i]) - $(xy[j]) meets edge " *
                "$(xy[k]) - $(xy[l])"
            )
        end
    end
    return length(edge_list)
end

# ---------------------------------------------------------------------------------------------
# The remaining region of a partition: band tops, end columns and wall curves, meshed by Gmsh.

struct RegionEdge
    a::Point2
    b::Point2
    transfinite::Bool # a band edge (1 segment) or a wall curve (Gmsh's size)
end

function region_loop_edges(
    loop::PartitionLoop,
    chains::Vector{MetalChain},
    chain_tops::Vector{Vector{Point2}},
    curves::Vector{PlanCurve}
)
    m = length(loop.curves)
    vertex_xy = loop.vertex_xy
    is_metal = [curves[k].metal for k in loop.curves]
    edges = RegionEdge[]
    if all(is_metal)
        length(chains) == 1 || error("A closed metal loop must carry one chain")
        top = chain_tops[1]
        for i in eachindex(top)
            push!(edges, RegionEdge(top[i], top[mod1(i + 1, length(top))], true))
        end
        return edges
    end
    # Replacement of a metal end vertex by its end column's top (on the wall).
    end_top = Dict{Int32, Point2}()
    for (chain, top) in zip(chains, chain_tops)
        end_top[chain.wall_start_vertex] = top[1]
        end_top[chain.wall_end_vertex] = top[end]
    end
    chain_by_start = Dict(chain.wall_start_vertex => c for (c, chain) in enumerate(chains))
    first_wall = findfirst(!, is_metal)
    order = [mod1(first_wall + s, m) for s = 1:m]
    k = 1
    while k <= m
        i = order[k]
        v, w = loop.vertices[i], loop.vertices[mod1(i + 1, m)]
        if is_metal[i]
            c = get(chain_by_start, v, 0)
            c > 0 || error("Metal chain start not found at vertex $v")
            top = chain_tops[c]
            for j = 1:(length(top) - 1)
                push!(edges, RegionEdge(top[j], top[j + 1], true))
            end
            # Skip the chain's remaining curves.
            while k <= m && is_metal[order[k]]
                k += 1
            end
            continue
        end
        # A straight wall run: consecutive collinear wall curves (the wall is fragmented by
        # the other plane's wall vertices, which carry no metal) between two metal ends or box
        # corners become ONE region edge, so an end column may span several fragments.
        run_start = i
        direction = unit((
            vertex_xy[mod1(i + 1, m)][1] - vertex_xy[i][1],
            vertex_xy[mod1(i + 1, m)][2] - vertex_xy[i][2]
        ))
        while k < m
            next_index = order[k + 1]
            is_metal[next_index] && break
            haskey(end_top, loop.vertices[next_index]) && break
            next_direction = unit((
                vertex_xy[mod1(next_index + 1, m)][1] - vertex_xy[next_index][1],
                vertex_xy[mod1(next_index + 1, m)][2] - vertex_xy[next_index][2]
            ))
            abs(direction[1] * next_direction[2] - direction[2] * next_direction[1]) <=
            1.0e-9 || break
            k += 1
            i = next_index
        end
        w = loop.vertices[mod1(i + 1, m)]
        p = get(end_top, v, vertex_xy[run_start])
        q = get(end_top, w, vertex_xy[mod1(i + 1, m)])
        remaining = (q[1] - p[1]) * direction[1] + (q[2] - p[2]) * direction[2]
        remaining > 1.0e-9 || error(
            "Band end columns overlap along the wall run $(vertex_xy[run_start]) - " *
            "$(vertex_xy[mod1(i + 1, m)]) (remaining $remaining um)"
        )
        push!(edges, RegionEdge(p, q, false))
        k += 1
    end
    return edges
end

function mesh_region(
    partition::PartitionGeometry,
    loops_edges::Vector{Vector{RegionEdge}},
    class::Int,
    raw_triangles::Vector{Tuple{NTuple{3, Point2}, Int}},
    plan_model::String
)
    gmsh.model.add("region_$(partition.surface)")
    geo = gmsh.model.geo
    point_tag = Dict{NTuple{2, Int64}, Int32}()
    key(p) = (round(Int64, p[1] / BAND_NODE_KEY_UM), round(Int64, p[2] / BAND_NODE_KEY_UM))
    point(p) = get!(point_tag, key(p)) do
        return geo.add_point(p[1], p[2], 0.0)
    end
    transfinite = Int32[]
    loop_tags = Int32[]
    region_area = 0.0
    for edges in loops_edges
        lines = Int32[]
        for edge in edges
            a, b = point(edge.a), point(edge.b)
            a != b || error(
                "Degenerate region edge at $(edge.a) of partition $(partition.surface)"
            )
            line = geo.add_line(a, b)
            push!(lines, line)
            edge.transfinite && push!(transfinite, line)
            region_area += 0.5 * (edge.a[1] * edge.b[2] - edge.b[1] * edge.a[2])
        end
        push!(loop_tags, geo.add_curve_loop(lines))
    end
    region_area > 0.0 || error(
        "Partition $(partition.surface): the band leaves no region (signed area $region_area)"
    )
    surface = geo.add_plane_surface(loop_tags)
    geo.synchronize()
    for line in transfinite
        gmsh.model.mesh.set_transfinite_curve(line, 2)
    end
    gmsh.model.mesh.generate(2)
    node_tags, node_coordinates, _ = gmsh.model.mesh.get_nodes()
    coordinate_by_tag = Dict{UInt64, Point2}()
    for (index, tag) in enumerate(node_tags)
        coordinate_by_tag[tag] =
            (node_coordinates[3index - 2], node_coordinates[3index - 1])
    end
    area = 0.0
    element_types, _, nodes_by_type = gmsh.model.mesh.get_elements(2, surface)
    for (element_type, element_nodes) in zip(element_types, nodes_by_type)
        _, dimension, _, node_count, _, primary =
            gmsh.model.mesh.get_element_properties(element_type)
        dimension == 2 || continue
        primary in (3, 4) || error("Expected triangular or quadrilateral region elements")
        for offset = 0:node_count:(length(element_nodes) - node_count)
            corners = [coordinate_by_tag[element_nodes[offset + i]] for i = 1:primary]
            split =
                primary == 3 ? ((corners[1], corners[2], corners[3]),) :
                ((corners[1], corners[2], corners[3]), (corners[1], corners[3], corners[4]))
            for triangle in split
                area += triangle_area(triangle...)
                push!(raw_triangles, (triangle, class))
            end
        end
    end
    gmsh.model.remove()
    gmsh.model.set_current(plan_model)
    return area, region_area
end

"""
    mesh_partition_band(partition, curves, nodes_by_curve, heights, radial_um, fan_turn,
                        box, class, raw_triangles, statistics, plan_model, collision_cell_um,
                        chains_out) -> (band area, region area, minimum rows, band triangles)

Build the structured band of one partition (its own elements), refuse collisions, mesh the
remaining region with Gmsh, append every triangle (band and region) to `raw_triangles` with
the partition's class, record every chain's columns in `chains_out` (closed flag, class,
columns; the graded sweep's input), and check band + region = the OCC area.
"""
function mesh_partition_band(
    partition::PartitionGeometry,
    curves::Vector{PlanCurve},
    nodes_by_curve::Dict{Int, Vector{Point2}},
    heights::Vector{Float64},
    radial_um::Float64,
    fan_turn::Float64,
    box,
    class::Int,
    raw_triangles::Vector{Tuple{NTuple{3, Point2}, Int}},
    statistics::BandStatistics,
    plan_model::String,
    collision_cell_um::Float64,
    chains_out::Vector{Tuple{Bool, Int, Vector{BandColumn}}}
)
    where = "partition $(partition.surface)"
    band = NTuple{3, Point2}[]
    loops_edges = Vector{RegionEdge}[]
    minimum_rows = length(heights)
    for loop in partition.loops
        chains = loop_metal_chains(loop, curves, nodes_by_curve)
        tops = Vector{Point2}[]
        for chain in chains
            segment_rows, segment_distances =
                chain_segment_rows(chain, curves, heights, radial_um, statistics)
            columns, _, node_rows = chain_columns(
                chain,
                segment_rows,
                segment_distances,
                heights,
                fan_turn,
                statistics
            )
            minimum_rows = min(minimum_rows, minimum(node_rows))
            for c = 1:(length(columns) - 1)
                push_band_elements!(band, columns[c], columns[c + 1], statistics, where)
            end
            chain.closed &&
                push_band_elements!(band, columns[end], columns[1], statistics, where)
            push!(tops, [column.nodes[end] for column in columns])
            push!(chains_out, (chain.closed, class, columns))
        end
        push!(loops_edges, region_loop_edges(loop, chains, tops, curves))
    end
    for (a, b, c) in band, p in (a, b, c)
        (
            box[1] - 1.0e-9 <= p[1] <= box[2] + 1.0e-9 &&
            box[3] - 1.0e-9 <= p[2] <= box[4] + 1.0e-9
        ) || error("Band node $p of $where lies outside the box")
    end
    isempty(band) || check_band_collisions(band, collision_cell_um, where)
    band_area = sum(triangle_area(t...) for t in band; init=0.0)
    for triangle in band
        push!(raw_triangles, (triangle, class))
    end
    region_mesh_area, region_polygon_area =
        mesh_region(partition, loops_edges, class, raw_triangles, plan_model)
    total = band_area + region_mesh_area
    abs(total - partition.area) <= REGION_AREA_TOLERANCE * max(1.0, partition.area) ||
        error(
            "$where: band $band_area + region $region_mesh_area (polygon $region_polygon_area) " *
            "= $total differs from the OCC area $(partition.area)"
        )
    return band_area, region_mesh_area, minimum_rows, length(band)
end
