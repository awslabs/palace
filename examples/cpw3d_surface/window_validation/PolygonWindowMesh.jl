# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabricated (finite-metal, overetched) reference mesher for validation WINDOWS given as a
# POLYGON SET (validation plan section (c), USER decision 184): the single transmon's
# shared-plan tensor-sweep generator generalised from the DeviceLayout geometry to
#
#   * an arbitrary number of conductors per plane, each labelled; the label `ground` is shared
#     by every plane and every bump, every other label is a terminal with its own metal-air /
#     metal-substrate attribute pair;
#   * one or two metal levels (planes): a plane facing `up` carries its metal on top of a
#     substrate below its surface z (the transmon chip), a plane facing `down` carries its
#     metal below a substrate above its surface (the flip-chip L2). Both planes' metal edges
#     constrain ONE plan mesh, swept through explicit z levels with each plane's metal band,
#     overetch band and the graded gap between the levels;
#   * bumps: plan footprints extruded as excluded metal columns from the L1 metal top to the
#     L2 metal bottom, joining the two grounds (their footprint must lie on metal of BOTH
#     planes);
#   * a truncated box: the plan rectangle, each plane's substrate thickness and optional
#     vacuum beyond the outermost substrate backsides.
#
# Process geometry as the transmon reference: metal thickness 0.1 um with vertical sidewalls,
# 50-nm overetch recess of the exposed substrate, metal volume excluded from the domain, the
# complete conforming substrate-air interface incl. the overetch step, separate metal-air /
# metal-substrate shells per conductor; SA / MS / MA participations are surface integrals on
# those shells (no meshed 2-nm layers). Resolution: edge-normal boundary layers of first
# height r (geometric growth 2, band ~1.55 um) on every metal edge of every plane and bump,
# tangent spacing t along the edges, r-resolved metal and overetch bands in z.
#
# Method (as generate_shared_plan_mesh.jl): the plan rectangle is fragmented by every polygon
# (all planes) and bump footprint; each resulting partition is classified by the fragment map
# (conductor per plane, bump); each partition gets a one-sided boundary-layer band on its
# metal-edge curves (transfinite tangent spacing) — by default OUR OWN structured band
# (structured_band.jl, `band_mode=:own`: per-segment / per-side band cap, deterministic, Gmsh
# triangulates only the remaining region; supervisor decisions 191 / 196 / 198), or Gmsh's
# BoundaryLayer field on the partition's copy (`band_mode=:gmsh`, the recorded generator's
# path) — the partitions' coincident nodes are welded by coordinate, the plan is swept through
# the z levels, prisms are split conformally into tetrahedra, the metal is omitted, unused
# nodes compacted and the physical attributes written directly to an ASCII MSH2 file with a
# JSON manifest (counts, per-attribute areas / volumes, first-layer heights, z levels).
#
# Deliberate deviation from the recorded transmon generator in Gmsh mode (supervisor decision
# 188): the BoundaryLayer Thickness is the geometric sum r (2^n - 1) times (1 + 1e-6). The
# recorded generator passes the exact sum, so floating-point rounding decides whether a
# column gets n or n - 1 rows (~8 % of the columns short at r10, ~20 % at r50, also in the
# recorded meshes), and on Linux that mix fails Gmsh's edge recovery at r10. With the margin
# every column has exactly n rows on every platform. `exact_band_thickness=true`
# (`--exact-band-thickness`) reproduces the recorded meshes.
#
# Polygon-set JSON (micrometres; `Version` 1):
#   {"Name": ..., "Box": {"X": [x0, x1], "Y": [y0, y1]},
#    "Process": {"MetalThickness": 0.1, "Overetch": 0.05},
#    "Planes": [{"Name": "L1", "SurfaceZ": 0.0, "Facing": "up", "SubstrateThickness": 525.0,
#                "Polygons": [{"Conductor": "ground", "Outer": [[x, y], ...],
#                              "Holes": [[[x, y], ...], ...]}, ...]}, ...],
#    "Bumps": [{"Conductor": "ground", "Footprint": [[x, y], ...]}],
#    "Vacuum": {"Below": 0.0, "Above": 0.0},
#    "Terminals": ["island", ...] (optional attribute order; default sorted labels)}
# Polygons of one plane must not overlap; vertices may lie on the box wall (a ground plane
# reaching the wall); wall edges are not metal edges (no boundary layer).
#
# Attributes: 3D 1 substrate, 2 vacuum; 2D 3 exterior_boundary, 4 ground_air,
# 5 ground_substrate, 6 substrate_air, 7 <terminal 1>_air, 8 <terminal 1>_substrate,
# 9 substrate_backside (a substrate face against vacuum beyond the box: SA diagnostic),
# 10 / 11 <terminal 2>_air / _substrate, ... (the transmon reference's table for one
# terminal). The manifest records the table.

module PolygonWindowMesh

import Gmsh: gmsh
using JSON
using Printf
using SHA
using Statistics

export read_polygon_set, mesh_polygon_window, plan_partitions, reconcile_planes

const Point2 = NTuple{2, Float64}

struct Polygon
    conductor::String
    outer::Vector{Point2}
    holes::Vector{Vector{Point2}}
end

struct Plane
    name::String
    surface_z::Float64
    facing::Int # +1: metal above the surface (substrate below); -1: metal below (flip-chip)
    substrate_thickness::Float64
    polygons::Vector{Polygon}
end

struct Bump
    conductor::String
    footprint::Vector{Point2}
end

struct PolygonSet
    name::String
    box::NTuple{4, Float64} # xmin, xmax, ymin, ymax
    metal_thickness::Float64
    overetch::Float64
    planes::Vector{Plane}
    bumps::Vector{Bump}
    vacuum_below::Float64
    vacuum_above::Float64
    terminals::Vector{String}
    matching_radius::Float64 # the identification's R (um); NaN when the set has none
end

const GROUND = "ground"

point_list(values) = [(Float64(p[1]), Float64(p[2])) for p in values]

function polygon_signed_area(points::Vector{Point2})
    n = length(points)
    return 0.5 * sum(
        points[i][1] * points[mod1(i + 1, n)][2] - points[mod1(i + 1, n)][1] * points[i][2]
        for i = 1:n
    )
end

function read_polygon_set(data::AbstractDict)
    get(data, "Version", 1) == 1 || error("Unsupported polygon-set version")
    box = data["Box"]
    x, y = box["X"], box["Y"]
    x[1] < x[2] && y[1] < y[2] || error("Degenerate box")
    process = get(data, "Process", Dict())
    planes = Plane[]
    labels = Set{String}()
    for plane in data["Planes"]
        facing = plane["Facing"]
        facing in ("up", "down") || error("Plane facing must be up or down")
        polygons = Polygon[]
        for polygon in plane["Polygons"]
            outer = point_list(polygon["Outer"])
            length(outer) >= 3 || error("Polygon outer loop needs 3 vertices")
            holes = [point_list(h) for h in get(polygon, "Holes", [])]
            all(h -> length(h) >= 3, holes) || error("Polygon hole needs 3 vertices")
            label = String(polygon["Conductor"])
            label == GROUND || push!(labels, label)
            push!(polygons, Polygon(label, outer, holes))
        end
        isempty(polygons) && error("Plane $(plane["Name"]) has no polygons")
        push!(
            planes,
            Plane(
                String(plane["Name"]),
                Float64(plane["SurfaceZ"]),
                facing == "up" ? 1 : -1,
                Float64(plane["SubstrateThickness"]),
                polygons
            )
        )
    end
    1 <= length(planes) <= 2 || error("One or two planes are supported")
    sort!(planes; by=p -> p.surface_z)
    if length(planes) == 2
        planes[1].facing == 1 && planes[2].facing == -1 || error(
            "Two planes must be a lower plane facing up and an upper plane facing down"
        )
        planes[2].surface_z - planes[1].surface_z >
        2 * Float64(get(process, "MetalThickness", 0.1)) ||
            error("Planes too close for two metal bands")
    end
    bumps = Bump[]
    for bump in get(data, "Bumps", [])
        footprint = point_list(bump["Footprint"])
        length(footprint) >= 3 || error("Bump footprint needs 3 vertices")
        label = String(get(bump, "Conductor", GROUND))
        label == GROUND || push!(labels, label)
        push!(bumps, Bump(label, footprint))
    end
    isempty(bumps) || length(planes) == 2 || error("Bumps need two planes")
    vacuum = get(data, "Vacuum", Dict())
    terminals = String.(get(data, "Terminals", sort!(collect(labels))))
    Set(terminals) == labels || error("Terminals must list every non-ground conductor once")
    length(unique(terminals)) == length(terminals) || error("Duplicate terminal label")
    return PolygonSet(
        String(get(data, "Name", "window")),
        (Float64(x[1]), Float64(x[2]), Float64(y[1]), Float64(y[2])),
        Float64(get(process, "MetalThickness", 0.1)),
        Float64(get(process, "Overetch", 0.05)),
        planes,
        bumps,
        Float64(get(vacuum, "Below", 0.0)),
        Float64(get(vacuum, "Above", 0.0)),
        terminals,
        haskey(data, "MatchingRadius") ? Float64(data["MatchingRadius"]) : NaN
    )
end

read_polygon_set(path::AbstractString) = read_polygon_set(JSON.parsefile(path))

# ---------------------------------------------------------------------------------------------
# Cross-plane reconciliation (supervisor decision on the S1p trial, 2026-10-01): edges of the
# two planes that are nominally coincident in plan (the aligned ground edges and rounded
# corners of both chips) arrive with different vertex samplings (chip-mesh nodes, chords), and
# the fragmentation then produces sliver partitions and sub-0.1-um pieces that break the
# boundary layers. Before meshing, every vertex of the UPPER plane within delta of a vertex
# (first) or else a segment of the lower plane is moved onto it (ties by coordinates, so the
# result does not depend on the input order); then every upper vertex lying exactly on a lower
# segment is inserted into the lower chain, and every lower vertex within delta of an upper
# segment is inserted into the upper chain (bending that edge by at most delta): coincident
# runs become identical point sequences. The lower plane's geometry never moves.
# delta = CROSS_PLANE_SNAP_FRACTION x MatchingRadius = 0.05 R (0.095 um at R = 1.9 um: the
# identification's own joint-noise resolution; cross-plane edges are separated vertically by the
# chip gap, so smaller plan offsets between them are physically irrelevant. Raised from 0.025 R
# after the S1p chamfers left a 0.029-R needle, supervisor decision 2026-10-01); an explicit
# `cross_plane_snap_um` override is recorded in the manifest. Within ONE plane nothing is
# snapped: two polygons closer than delta (a sub-delta slot or touching metal) are refused and
# the writer resolves them. The reference therefore differs from the thin (chip-mesh) geometry
# by at most delta on the reconciled runs; the manifest records the moved vertices, the maximum
# displacement, the inserted vertices and the coincident run length.
#
# Window-cut sliver rules (supervisor decision 306; the stage-2 C1 / C3 windows): the d210
# extraction cuts the chip polygons at the box walls, so chip-mesh nodes sit within delta of
# the walls and of each other's chains, and the reconciliation above produced two degenerate
# loops that no mesher can take:
#   * a vertex within delta of TWO consecutive upper segments (near their shared vertex) was
#     inserted into both, giving the zero-area spike P, B, P (C3: an L1 ground vertex 0.073 um
#     from a snapped L2 finger corner; Gmsh "Curve loop is not closed"). Rule: a vertex is
#     inserted into at most ONE segment per loop, the nearest (ties by the foot's
#     coordinates), and never into a loop it already belongs to;
#   * a lower vertex 0.015 um below the top wall was inserted into the upper plane's WALL
#     segment, bending the window-cut line into a 0.022 x 0.015-um notch partition (C1: a 2 r x
#     1.5 r sliver whose wall-end columns collide). Rule: a segment on a box wall (a window-cut
#     line, not a metal edge) is never bent — it takes only vertices lying on the wall — and a
#     vertex on a wall snaps only along that wall (to lower vertices on it or lower segments
#     along it), so the planes' outlines stay on the cut lines;
#   * after the snaps and insertions every loop is cleaned deterministically: consecutive
#     vertices within the plan-node identity quantum PLAN_NODE_MERGE_TOLERANCE_UM (the
#     coordinate welding tolerance of `mesh_plan`: such vertices would be ONE plan node and a
#     zero-length OCC edge) are merged onto the first, a spike (a vertex whose two neighbours
#     coincide) collapses onto its neighbours, repeated until stable; a loop that still visits
#     a point twice (self-touching) is refused naming the point.
# None of the rules carries a case constant (delta and the identity quantum are the existing
# ones); every count is recorded in the manifest's `cross_plane_reconciliation`, and where no
# such configuration occurs the result is identical to the rules above (S1p / C2p regression).

const CROSS_PLANE_SNAP_FRACTION = 0.05
const ON_SEGMENT_TOLERANCE_UM = 1.0e-9
# Plan nodes closer than this are welded into one node by `mesh_plan` (coordinate keys), so two
# polygon vertices this close are one vertex.
const PLAN_NODE_MERGE_TOLERANCE_UM = 1.0e-6

point_distance(p::Point2, q::Point2) = hypot(p[1] - q[1], p[2] - q[2])

# Distance from p to the segment ab, the clamped parameter and the foot point.
function segment_projection(p::Point2, a::Point2, b::Point2)
    abx, aby = b[1] - a[1], b[2] - a[2]
    length2 = abx^2 + aby^2
    t =
        length2 == 0.0 ? 0.0 :
        clamp(((p[1] - a[1]) * abx + (p[2] - a[2]) * aby) / length2, 0.0, 1.0)
    foot = (a[1] + t * abx, a[2] + t * aby)
    return point_distance(p, foot), t, foot
end

plane_loops(plane::Plane) =
    [loop for polygon in plane.polygons for loop in (polygon.outer, polygon.holes...)]

loop_segments(loop::Vector{Point2}) =
    [(loop[i], loop[mod1(i + 1, length(loop))]) for i in eachindex(loop)]

# Polygons of one plane closer than delta (a sub-delta slot or touching metal): refused.
function check_same_plane_separation(plane::Plane, delta::Float64)
    for (i, a) in enumerate(plane.polygons), (j, b) in enumerate(plane.polygons)
        i == j && continue
        for loop in (a.outer, a.holes...), v in loop
            for other in (b.outer, b.holes...), (p, q) in loop_segments(other)
                distance, _, _ = segment_projection(v, p, q)
                distance < delta && error(
                    "Plane $(plane.name): polygons $i ($(a.conductor)) and $j " *
                    "($(b.conductor)) are $distance um apart at $v, closer than the " *
                    "cross-plane snap distance $delta um (a sub-delta slot or touching " *
                    "metal; the writer must resolve it)"
                )
            end
        end
    end
end

# The box walls a point lies on exactly (the window-cut lines): 1 / 2 x = xmin / xmax, 3 / 4
# y = ymin / ymax; a box corner lies on two.
function point_walls(p::Point2, box)
    walls = Int[]
    p[1] == box[1] && push!(walls, 1)
    p[1] == box[2] && push!(walls, 2)
    p[2] == box[3] && push!(walls, 3)
    p[2] == box[4] && push!(walls, 4)
    return walls
end

# The walls a segment lies along (both ends on the same wall).
segment_walls(a::Point2, b::Point2, box) =
    intersect(point_walls(a, box), point_walls(b, box))

# Move v onto the nearest target vertex within delta, else onto the nearest target segment
# within delta (equal distances: the lexicographically smallest target point); else keep v.
function snap_point(v::Point2, vertices::Vector{Point2}, segments, delta::Float64)
    best, best_distance = v, Inf
    for q in vertices
        d = point_distance(v, q)
        d <= delta || continue
        if d < best_distance - ON_SEGMENT_TOLERANCE_UM ||
           (abs(d - best_distance) <= ON_SEGMENT_TOLERANCE_UM && q < best)
            best, best_distance = q, d
        end
    end
    isfinite(best_distance) && return best
    for (a, b) in segments
        d, _, foot = segment_projection(v, a, b)
        d <= delta || continue
        if d < best_distance - ON_SEGMENT_TOLERANCE_UM ||
           (abs(d - best_distance) <= ON_SEGMENT_TOLERANCE_UM && foot < best)
            best, best_distance = foot, d
        end
    end
    return best
end

# Vertices on a box wall move only along that wall: the candidates are the lower vertices on
# the wall and the lower segments along it (their feet lie on the wall). Returns the snapped
# point and whether the wall constraint changed the unconstrained result.
function snap_wall_point(
    v::Point2,
    vertices::Vector{Point2},
    segments,
    delta::Float64,
    box,
    walls::Vector{Int}
)
    on_same_wall(q) = any(in(walls), point_walls(q, box))
    wall_vertices = filter(on_same_wall, vertices)
    wall_segments =
        [(a, b) for (a, b) in segments if any(in(walls), segment_walls(a, b, box))]
    w = snap_point(v, wall_vertices, wall_segments, delta)
    return w, w != snap_point(v, vertices, segments, delta)
end

# Consecutive vertices within `tolerance` are one vertex (the first kept); returns the loop and
# the number of dropped vertices.
function drop_repeated_vertices(loop::Vector{Point2}, tolerance::Float64=0.0)
    cleaned = Point2[]
    dropped = 0
    for p in loop
        if !isempty(cleaned) && point_distance(p, cleaned[end]) <= tolerance
            dropped += 1
            continue
        end
        push!(cleaned, p)
    end
    while length(cleaned) > 1 && point_distance(cleaned[end], cleaned[1]) <= tolerance
        pop!(cleaned)
        dropped += 1
    end
    return cleaned, dropped
end

# Insert every point of `points` that is not yet a vertex of `loop` into the NEAREST segment of
# the loop within `tolerance` (foot strictly inside the segment; equal distances: the
# lexicographically smallest foot) — one segment per point, so a point near a shared vertex of
# two segments never appears twice in the loop. A segment lying on a box wall (a window-cut
# line) takes only points on the wall (within ON_SEGMENT_TOLERANCE_UM): the wall is never bent.
# Returns the loop, the number of inserted points, the largest distance by which an inserted
# point bends its segment, the number of points that had more than one admissible segment and
# the number of points kept off a wall segment they were within `tolerance` of.
function insert_points_on_segments(
    loop::Vector{Point2},
    points::Vector{Point2},
    tolerance::Float64,
    box
)
    segments = loop_segments(loop)
    wall_segment = [!isempty(segment_walls(a, b, box)) for (a, b) in segments]
    hits = [Tuple{Float64, Float64, Point2}[] for _ in segments]
    vertex_set = Set(loop)
    multiple, wall_refused = 0, 0
    for p in points
        p in vertex_set && continue
        best, best_distance, best_t, best_foot = 0, Inf, 0.0, p
        admissible, refused = 0, false
        for (s, (a, b)) in enumerate(segments)
            d, t, foot = segment_projection(p, a, b)
            d <= tolerance && 0.0 < t < 1.0 || continue
            if wall_segment[s] && d > ON_SEGMENT_TOLERANCE_UM
                refused = true
                continue
            end
            admissible += 1
            if d < best_distance - ON_SEGMENT_TOLERANCE_UM ||
               (abs(d - best_distance) <= ON_SEGMENT_TOLERANCE_UM && foot < best_foot)
                best, best_distance, best_t, best_foot = s, d, t, foot
            end
        end
        refused && (wall_refused += 1)
        best == 0 && continue
        admissible > 1 && (multiple += 1)
        push!(hits[best], (best_t, best_distance, p))
    end
    result = Point2[]
    inserted, max_bend = 0, 0.0
    for (s, (a, _)) in enumerate(segments)
        push!(result, a)
        for (_, d, p) in sort!(hits[s])
            p == result[end] && continue
            push!(result, p)
            inserted += 1
            max_bend = max(max_bend, d)
        end
    end
    cleaned, _ = drop_repeated_vertices(result)
    return cleaned, inserted, max_bend, multiple, wall_refused
end

# Deterministic clean-up of a reconciled loop: consecutive vertices within `tolerance` merged
# onto the first, a spike (a vertex whose two neighbours coincide within `tolerance`: a
# zero-area excursion) collapsed onto its neighbours, repeated until stable. Returns the loop,
# the merged count and the collapsed-spike count; a loop that still visits a point twice is
# refused with the point (`where` names the plane and loop).
function clean_loop(loop::Vector{Point2}, tolerance::Float64, where::String)
    merged, spikes = 0, 0
    while true
        loop, dropped = drop_repeated_vertices(loop, tolerance)
        merged += dropped
        n = length(loop)
        spike =
            n >= 3 ?
            findfirst(
                i ->
                    point_distance(loop[mod1(i - 1, n)], loop[mod1(i + 1, n)]) <= tolerance,
                1:n
            ) : nothing
        spike === nothing && break
        deleteat!(loop, sort!(unique([spike, mod1(spike + 1, n)])))
        spikes += 1
    end
    seen = Dict{NTuple{2, Int64}, Point2}()
    for p in loop
        key = (round(Int64, p[1] / tolerance), round(Int64, p[2] / tolerance))
        haskey(seen, key) && error(
            "$where visits $p twice after the cross-plane reconciliation (a self-touching " *
            "loop; the writer must resolve it)"
        )
        seen[key] = p
    end
    return loop, merged, spikes
end

function rebuild_plane(plane::Plane, loops::Vector{Vector{Point2}})
    polygons = Polygon[]
    k = 0
    for polygon in plane.polygons
        outer = loops[k + 1]
        holes = loops[(k + 2):(k + 1 + length(polygon.holes))]
        k += 1 + length(polygon.holes)
        all(l -> length(l) >= 3, (outer, holes...)) ||
            error("Polygon of plane $(plane.name) degenerated by the reconciliation")
        push!(polygons, Polygon(polygon.conductor, outer, holes))
    end
    return Plane(
        plane.name,
        plane.surface_z,
        plane.facing,
        plane.substrate_thickness,
        polygons
    )
end

"""
    reconcile_planes(spec, delta) -> (reconciled spec, report)

Snap the upper plane's vertices within `delta` onto the lower plane (vertices first, then
segments; a vertex on a box wall only along that wall), insert the cross vertices on both planes
(exactly into the lower chains, within `delta` into the upper chains; one segment per vertex,
wall segments never bent), clean every loop (merged vertices, collapsed spikes; a self-touching
loop refused), refuse same-plane polygons closer than `delta`. A single-plane set is returned
unchanged.
"""
function reconcile_planes(spec::PolygonSet, delta::Float64)
    report = Dict{String, Any}("applied" => false)
    length(spec.planes) == 2 || return spec, report
    delta > 0.0 || error("Cross-plane snap distance must be positive")
    lower, upper = spec.planes
    box = spec.box
    for plane in spec.planes
        check_same_plane_separation(plane, delta)
    end
    lower_loops = [copy(loop) for loop in plane_loops(lower)]
    upper_loops = [copy(loop) for loop in plane_loops(upper)]
    lower_vertices = reduce(vcat, lower_loops)
    lower_segments = reduce(vcat, loop_segments.(lower_loops))
    moved, max_displacement, wall_constrained = 0, 0.0, 0
    for (k, loop) in enumerate(upper_loops)
        snapped = Point2[]
        for v in loop
            walls = point_walls(v, box)
            w = if isempty(walls)
                snap_point(v, lower_vertices, lower_segments, delta)
            else
                w, constrained = snap_wall_point(
                    v,
                    lower_vertices,
                    lower_segments,
                    delta,
                    box,
                    walls
                )
                constrained && (wall_constrained += 1)
                w
            end
            if w != v
                moved += 1
                max_displacement = max(max_displacement, point_distance(v, w))
            end
            push!(snapped, w)
        end
        upper_loops[k], _ = drop_repeated_vertices(snapped)
    end
    upper_vertices = reduce(vcat, upper_loops)
    inserted_lower, inserted_upper = 0, 0
    single_segment, wall_refused = 0, 0
    for (k, loop) in enumerate(lower_loops)
        lower_loops[k], n, _, multiple, refused =
            insert_points_on_segments(loop, upper_vertices, ON_SEGMENT_TOLERANCE_UM, box)
        inserted_lower += n
        single_segment += multiple
        wall_refused += refused
    end
    lower_vertices = reduce(vcat, lower_loops)
    for (k, loop) in enumerate(upper_loops)
        upper_loops[k], n, bend, multiple, refused =
            insert_points_on_segments(loop, lower_vertices, delta, box)
        inserted_upper += n
        max_displacement = max(max_displacement, bend)
        single_segment += multiple
        wall_refused += refused
    end
    merged, spikes = 0, 0
    for (plane, loops) in ((lower, lower_loops), (upper, upper_loops)),
        k in eachindex(loops)

        loops[k], m, s = clean_loop(
            loops[k],
            PLAN_NODE_MERGE_TOLERANCE_UM,
            "Plane $(plane.name) loop $k"
        )
        merged += m
        spikes += s
    end
    # Coincident metal runs (segments on the box wall are not metal edges and not counted).
    on_wall((a, b)) = !isempty(segment_walls(a, b, box))
    lower_set = Set(minmax(a, b) for (a, b) in reduce(vcat, loop_segments.(lower_loops)))
    coincident = [
        (a, b) for (a, b) in reduce(vcat, loop_segments.(upper_loops)) if
        minmax(a, b) in lower_set && !on_wall((a, b))
    ]
    report = Dict{String, Any}(
        "applied" => true,
        "delta_um" => delta,
        "moved_vertices" => moved,
        "max_displacement_um" => max_displacement,
        "inserted_vertices" =>
            Dict(lower.name => inserted_lower, upper.name => inserted_upper),
        "coincident_segments" => length(coincident),
        "coincident_run_length_um" =>
            sum(point_distance(a, b) for (a, b) in coincident; init=0.0),
        "window_cut_sliver_rules" => Dict{String, Any}(
            "rules" =>
                "a vertex on a box wall snaps only along that wall; a wall segment " *
                "(window-cut line) takes only vertices on the wall; a vertex is inserted " *
                "into at most one segment per loop (the nearest); consecutive vertices " *
                "within the plan-node identity quantum merged, spikes collapsed",
            "wall_constrained_snaps" => wall_constrained,
            "wall_insertions_refused" => wall_refused,
            "single_segment_insertions" => single_segment,
            "merge_tolerance_um" => PLAN_NODE_MERGE_TOLERANCE_UM,
            "merged_consecutive_vertices" => merged,
            "collapsed_spikes" => spikes
        )
    )
    reconciled = PolygonSet(
        spec.name,
        spec.box,
        spec.metal_thickness,
        spec.overetch,
        [rebuild_plane(lower, lower_loops), rebuild_plane(upper, upper_loops)],
        spec.bumps,
        spec.vacuum_below,
        spec.vacuum_above,
        spec.terminals,
        spec.matching_radius
    )
    return reconciled, report
end

# The snap distance of a set: the rule 0.05 R, or an explicit override (recorded).
function cross_plane_snap_distance(spec::PolygonSet, override::Float64)
    length(spec.planes) == 2 || return NaN, "none"
    isnan(override) || return override, "override"
    isnan(spec.matching_radius) && error(
        "Two planes need MatchingRadius (um) in the polygon set (delta = " *
        "$(CROSS_PLANE_SNAP_FRACTION) R) or an explicit cross_plane_snap_um"
    )
    return CROSS_PLANE_SNAP_FRACTION * spec.matching_radius,
    "$(CROSS_PLANE_SNAP_FRACTION) x MatchingRadius"
end

# Attribute table: 4 / 5 ground, 7 / 8 the first terminal, 10 / 11 the second, ...
function attribute_table(spec::PolygonSet)
    table = Dict{String, Tuple{Int, Int}}(GROUND => (4, 5))
    for (index, label) in enumerate(spec.terminals)
        table[label] = index == 1 ? (7, 8) : (10 + 2 * (index - 2), 11 + 2 * (index - 2))
    end
    return table
end

function physical_names(spec::PolygonSet)
    table = attribute_table(spec)
    names = [
        (3, 1, "substrate"),
        (3, 2, "vacuum"),
        (2, 3, "exterior_boundary"),
        (2, 4, "ground_air"),
        (2, 5, "ground_substrate"),
        (2, 6, "substrate_air"),
        (2, 9, "substrate_backside")
    ]
    for label in spec.terminals
        air, substrate = table[label]
        push!(names, (2, air, label * "_air"), (2, substrate, label * "_substrate"))
    end
    return sort!(names; by=entry -> (entry[1], entry[2]))
end

# ---------------------------------------------------------------------------------------------
# Plan geometry: the rectangle fragmented by every polygon of every plane and every bump.

# Gmsh's OCC plane surface expects every hole wire with the SAME orientation as the outer
# wire (it reverses the holes itself), so every loop is passed counterclockwise.
function add_polygon_surface(occ, outer::Vector{Point2}, holes::Vector{Vector{Point2}})
    function loop(points)
        polygon_signed_area(points) < 0.0 && (points = reverse(points))
        tags = [occ.add_point(p[1], p[2], 0.0) for p in points]
        lines = [
            occ.add_line(tags[i], tags[mod1(i + 1, length(tags))]) for i in eachindex(tags)
        ]
        return occ.add_curve_loop(lines)
    end
    return occ.add_plane_surface(vcat(loop(outer), [loop(h) for h in holes]))
end

# Partition classification: conductor label per plane ("" = none) and bump index (0 = none).
struct PartitionClass
    conductors::Vector{String}
    bump::Int
end
Base.:(==)(a::PartitionClass, b::PartitionClass) =
    a.conductors == b.conductors && a.bump == b.bump
Base.hash(c::PartitionClass, h::UInt) = hash((c.conductors, c.bump), h)

function on_box_wall(box, bbox; tolerance=1.0e-5)
    return (abs(bbox[1] - box[1]) <= tolerance && abs(bbox[4] - box[1]) <= tolerance) ||
           (abs(bbox[1] - box[2]) <= tolerance && abs(bbox[4] - box[2]) <= tolerance) ||
           (abs(bbox[2] - box[3]) <= tolerance && abs(bbox[5] - box[3]) <= tolerance) ||
           (abs(bbox[2] - box[4]) <= tolerance && abs(bbox[5] - box[4]) <= tolerance)
end

"""
    plan_partitions(spec) -> (surfaces, classes, metal_edge_curves)

Build the fragmented plan in the current Gmsh model. Returns the partition surfaces
(dimtags), their classes, and the set of curve tags that are metal edges (a curve whose two
adjacent partitions differ in conductor on some plane or in bump membership).
"""
function plan_partitions(spec::PolygonSet)
    occ = gmsh.model.occ
    box = spec.box
    rectangle = occ.add_rectangle(box[1], box[3], 0.0, box[2] - box[1], box[4] - box[3])
    tools = Tuple{Int32, Int32}[]
    tool_meaning = Tuple{Int, Int, String}[] # (plane index or 0, bump index or 0, label)
    for (plane_index, plane) in enumerate(spec.planes), polygon in plane.polygons
        push!(tools, (2, add_polygon_surface(occ, polygon.outer, polygon.holes)))
        push!(tool_meaning, (plane_index, 0, polygon.conductor))
    end
    for (bump_index, bump) in enumerate(spec.bumps)
        push!(tools, (2, add_polygon_surface(occ, bump.footprint, Vector{Point2}[])))
        push!(tool_meaning, (0, bump_index, bump.conductor))
    end
    fragments, fragment_map = occ.fragment([(2, rectangle)], tools)
    occ.synchronize()
    all(entity -> entity[1] == 2, fragments) ||
        error("Plan fragmentation produced non-surface entities")
    rectangle_fragments = Set(fragment_map[1])
    surfaces = sort!(collect(rectangle_fragments))
    n_planes = length(spec.planes)
    classes = Dict(surface => PartitionClass(fill("", n_planes), 0) for surface in surfaces)
    for (tool_index, (plane_index, bump_index, label)) in enumerate(tool_meaning)
        for fragment in fragment_map[tool_index + 1]
            haskey(classes, fragment) || error(
                "A polygon of $(label) leaves the plan rectangle (fragment $fragment)"
            )
            class = classes[fragment]
            if plane_index > 0
                isempty(class.conductors[plane_index]) || error(
                    "Polygons overlap on plane $(spec.planes[plane_index].name) " *
                    "($(class.conductors[plane_index]) and $label)"
                )
                class.conductors[plane_index] = label
            else
                class.bump == 0 || error("Bump footprints overlap")
                classes[fragment] = PartitionClass(class.conductors, bump_index)
            end
        end
    end
    for surface in surfaces
        class = classes[surface]
        if class.bump > 0
            all(!isempty, class.conductors) ||
                error("Bump $(class.bump) footprint is not on metal of both planes")
        end
    end
    # Metal edges: curves between partitions of different class, strictly inside the box.
    metal_edge_curves = Set{Int32}()
    for (_, tag) in unique(gmsh.model.get_boundary(surfaces, false, false, false))
        tag = abs(tag)
        upward, _ = gmsh.model.get_adjacencies(1, tag)
        owners = [(2, Int32(s)) for s in upward if haskey(classes, (2, Int32(s)))]
        if length(owners) == 1
            on_box_wall(box, gmsh.model.get_bounding_box(1, tag)) ||
                error("Plan curve $tag bounds one partition but is not on the box wall")
        elseif length(owners) == 2
            classes[owners[1]] == classes[owners[2]] || push!(metal_edge_curves, tag)
        else
            error("Plan curve $tag bounds $(length(owners)) partitions")
        end
    end
    isempty(metal_edge_curves) && error("No metal edges found")
    return surfaces, classes, metal_edge_curves
end

# Curves are matched between the original and its copies by geometry.
function curve_signature(dimension, tag)
    entity = (dimension, abs(tag))
    bbox = gmsh.model.get_bounding_box(entity...)
    length_um = gmsh.model.occ.get_mass(entity...)
    bounds = gmsh.model.get_parametrization_bounds(entity...)
    mid = gmsh.model.get_value(entity..., [0.5 * (bounds[1][1] + bounds[2][1])])
    return Tuple(round.(collect((bbox..., length_um, mid[1], mid[2])); digits=6))
end

function triangle_area(a::Point2, b::Point2, c::Point2)
    return 0.5 * abs((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
end

function metric_summary(values)
    return Dict(
        "minimum" => minimum(values),
        "p01" => quantile(values, 0.01),
        "median" => median(values),
        "p99" => quantile(values, 0.99),
        "maximum" => maximum(values)
    )
end

# ---------------------------------------------------------------------------------------------
# z levels

const SUBSTRATE_SIDE_OFFSETS_UM = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
const VACUUM_SIDE_OFFSETS_UM =
    [0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]

struct ZStack
    levels::Vector{Float64}
    z_bottom::Float64
    z_top::Float64
    backsides::Vector{Float64} # substrate faces against vacuum beyond the box
    gap_midpoint::Float64 # two planes: the level halfway between the metal tops (else NaN)
end

const STEP_FACE_Z_GRADINGS = (:legacy, :mirrored)

"""
    step_face_z_levels(spec, radial_um) -> Vector{Float64}

The `step_face_z_grading = :mirrored` levels (D1 of the reference-quality lane, supervisor
decision 406): the plan band's row heights r (2^k - 1), k = 1, 2, ..., mirrored into z beyond
BOTH faces of every plane's fabricated step — metal top + r, 3r, 7r, ... and trench bottom
- r, 3r, 7r, ... — while they stay strictly inside the first fixed cell beyond that face (the
first of `VACUUM_SIDE_OFFSETS_UM` above the surface, of `SUBSTRATE_SIDE_OFFSETS_UM` below),
inside the plane's substrate and, with two planes, on the plane's own side of the gap
midpoint. The legacy stack (`:legacy`, the recorded stage-1 / S2.1 reference family)
resolves z in r only inside the step, so the air cell at the metal top corner and the
substrate cell under the trench foot are r x 50 nm at every r and only p resolves the corner
singularities there (the measured MA / SA under-read of the window references). The graded
sweep thins these levels away from the edge like any level (they are not interfaces).
"""
function step_face_z_levels(spec::PolygonSet, radial_um::Float64)
    radial_um > 0.0 || error("step_face_z_levels needs a positive radial size")
    m, o = spec.metal_thickness, spec.overetch
    planes = spec.planes
    midpoint =
        length(planes) == 2 ?
        0.5 * sum(plane.surface_z + plane.facing * m for plane in planes) : NaN
    levels = Float64[]
    for plane in planes
        s, f = plane.surface_z, plane.facing
        faces = (
            # (face z, outward sign, room to the first fixed level beyond the face)
            (s + f * m, f, VACUUM_SIDE_OFFSETS_UM[1] - m),
            (s - f * o, -f, min(SUBSTRATE_SIDE_OFFSETS_UM[1], plane.substrate_thickness) - o)
        )
        for (face, sign, room) in faces
            k = 1
            while true
                offset = radial_um * (2.0^k - 1.0)
                offset < room - 1.0e-9 || break
                z = face + sign * offset
                # Two planes: the vacuum side stops at the gap midpoint like the fixed offsets.
                (isnan(midpoint) || sign != f || f * (midpoint - z) > 0.0) || break
                push!(levels, z)
                k += 1
            end
        end
    end
    return sort!(levels)
end

"""
    z_levels(spec, metal_layers, trench_layers; step_face_levels=Float64[]) -> ZStack

The production z stack: the r-resolved trench and metal bands of every plane, the fixed
substrate- and vacuum-side offsets, the backsides, the box ends and, with two planes, the gap
midpoint; `step_face_levels` (`step_face_z_levels`) are merged in as plain levels.
"""
function z_levels(
    spec::PolygonSet,
    metal_layers::Int,
    trench_layers::Int;
    step_face_levels::Vector{Float64}=Float64[]
)
    m, o = spec.metal_thickness, spec.overetch
    planes = spec.planes
    levels = copy(step_face_levels)
    backsides = Float64[]
    metal_tops = Float64[]
    midpoint = NaN
    for plane in planes
        s, f = plane.surface_z, plane.facing
        # The r-resolved trench and metal bands, then the graded substrate side.
        append!(levels, range(s - f * o, s; length=trench_layers + 1))
        append!(levels, range(s, s + f * m; length=metal_layers + 1))
        push!(metal_tops, s + f * m)
        for offset in SUBSTRATE_SIDE_OFFSETS_UM
            offset < plane.substrate_thickness || break
            push!(levels, s - f * offset)
        end
        push!(levels, s - f * plane.substrate_thickness)
    end
    if length(planes) == 2
        # Gap between the metal tops: each plane grades toward the midpoint.
        midpoint = 0.5 * (metal_tops[1] + metal_tops[2])
        for plane in planes
            s, f = plane.surface_z, plane.facing
            for offset in VACUUM_SIDE_OFFSETS_UM
                z = s + f * offset
                f * (midpoint - z) > 0.0 || break
                push!(levels, z)
            end
        end
        push!(levels, midpoint)
        lower, upper = planes
        z_bottom = lower.surface_z - lower.substrate_thickness - spec.vacuum_below
        z_top = upper.surface_z + upper.substrate_thickness + spec.vacuum_above
        spec.vacuum_below > 0.0 &&
            push!(backsides, lower.surface_z - lower.substrate_thickness)
        spec.vacuum_above > 0.0 &&
            push!(backsides, upper.surface_z + upper.substrate_thickness)
    else
        plane = planes[1]
        s, f = plane.surface_z, plane.facing
        for offset in VACUUM_SIDE_OFFSETS_UM
            push!(levels, s + f * offset)
        end
        backside = s - f * plane.substrate_thickness
        if f == 1
            spec.vacuum_above > m || error("Vacuum.Above must exceed the metal thickness")
            z_bottom = backside - spec.vacuum_below
            z_top = s + spec.vacuum_above
            spec.vacuum_below > 0.0 && push!(backsides, backside)
        else
            spec.vacuum_below > m || error("Vacuum.Below must exceed the metal thickness")
            z_bottom = s - spec.vacuum_below
            z_top = backside + spec.vacuum_above
            spec.vacuum_above > 0.0 && push!(backsides, backside)
        end
    end
    push!(levels, z_bottom, z_top)
    filter!(z -> z_bottom - 1.0e-9 <= z <= z_top + 1.0e-9, levels)
    sort!(levels)
    merged = Float64[levels[1]]
    for z in levels[2:end]
        z - merged[end] > 1.0e-9 && push!(merged, z)
    end
    return ZStack(merged, z_bottom, z_top, backsides, midpoint)
end

# Material of a prism at height z for a partition class: 1 substrate, 2 vacuum, 0 excluded.
function material(spec::PolygonSet, class::PartitionClass, z::Float64)
    m, o = spec.metal_thickness, spec.overetch
    for (k, plane) in enumerate(spec.planes)
        s, f = plane.surface_z, plane.facing
        d = f * (z - s) # signed distance from the surface toward the metal side
        metal = !isempty(class.conductors[k])
        if -plane.substrate_thickness < d < -o
            return Int8(1)
        elseif -o < d < 0.0
            return metal ? Int8(1) : Int8(2)
        elseif 0.0 < d < m
            metal && return Int8(0)
        end
    end
    if class.bump > 0
        lower, upper = spec.planes
        lower.surface_z + m < z < upper.surface_z - m && return Int8(0)
    end
    return Int8(2)
end

# ---------------------------------------------------------------------------------------------

# Relative margin on the BoundaryLayer Thickness (decision 188): every column gets exactly
# `radial_layers` rows instead of n or n - 1 by rounding of the exact geometric sum.
const BAND_THICKNESS_MARGIN = 1.0e-6

# ---------------------------------------------------------------------------------------------
# Local band cap (supervisor decision on the S1p trial, 2026-10-01): Gmsh's 2D boundary-layer
# extrusion fails ("Edge not recovered" on the layer front, or inverted quads) where the fronts
# of two facing metal edges overlap — a 1-um slot at r10 with the 1.27-um band, the SCT junction
# region (2-um fingers, L1 / L2 edges 0.75 um apart in plan). Resolution rule, not a failure
# fallback: at every point of a metal edge the band thickness is at most BAND_CAP_FRACTION times
# the local facing distance = the distance to the nearest NON-ADJACENT metal edge of any plane
# in the shared plan (curves sharing an endpoint are corners, not facing fronts). The first rows
# r, 2r, 4r, ... stay; only the outer graded rows are dropped: rows(d) = floor(log2(1 + 0.4 d /
# r)) capped at the global row count. Gmsh's field takes one Thickness per curve and terminates
# a layer only at a corner (the end column runs along the adjacent curve), so the row count is
# the minimum along every straight run of collinear curves and changes only at corners; the
# curves of a partition are grouped into one BoundaryLayer field per row count, each with its
# layer ends declared. Applied symmetrically (the rule is per curve, both facing fronts see the
# same distance). Manifest: capped length, the minimum band, length per row count.
const BAND_CAP_FRACTION = 0.4
const BAND_CAP_SAMPLE_SPACING_UM = 0.5
const BAND_CAP_COLLINEARITY = 1.0e-3 # |sin| of the angle between curves of one straight run

struct MetalCurve
    tag::Int32
    a::Point2
    b::Point2
    points::NTuple{2, Int32} # endpoint tags (adjacency)
end

function metal_curves(metal_edge_curves)
    curves = MetalCurve[]
    for tag in sort!(collect(metal_edge_curves))
        gmsh.model.get_type(1, tag) == "Line" ||
            error("Plan curve $tag is not a straight line")
        points =
            [abs(p) for (_, p) in gmsh.model.get_boundary([(1, tag)], false, false, false)]
        length(points) == 2 || error("Plan curve $tag has $(length(points)) end points")
        bounds = gmsh.model.get_parametrization_bounds(1, tag)
        a = gmsh.model.get_value(1, tag, [bounds[1][1]])
        b = gmsh.model.get_value(1, tag, [bounds[2][1]])
        push!(curves, MetalCurve(tag, (a[1], a[2]), (b[1], b[2]), (points[1], points[2])))
    end
    return curves
end

# Distance from p (on curve `index`) to the nearest metal curve that faces it: p's foot lies
# inside the other curve, and the other curve is not on the line through p (a piece of the
# same straight run, or a collinear edge ahead: the layers are normal to the run, never meet).
function facing_distance(p::Point2, index::Int, curves::Vector{MetalCurve})
    best = Inf
    for (j, other) in enumerate(curves)
        j == index && continue
        dx, dy = other.b[1] - other.a[1], other.b[2] - other.a[2]
        line_distance =
            abs((p[1] - other.a[1]) * dy - (p[2] - other.a[2]) * dx) / hypot(dx, dy)
        line_distance <= ON_SEGMENT_TOLERANCE_UM && continue
        d, t, _ = segment_projection(p, other.a, other.b)
        # Facing means the foot of p lies inside the other curve (a corner-adjacent curve at
        # an obtuse or right angle is not facing; a V-shaped sliver's other side is).
        1.0e-7 < t < 1.0 - 1.0e-7 || continue
        best = min(best, d)
    end
    return best
end

band_cap_rows(distance::Float64, radial_um::Float64, rows_max::Int) =
    isfinite(distance) ?
    clamp(floor(Int, log2(1.0 + BAND_CAP_FRACTION * distance / radial_um)), 1, rows_max) :
    rows_max

"""
    band_cap_rows_by_curve(metal_edge_curves, radial_um, rows_max) -> rows by curve signature

Row count of every fragmented metal curve: the minimum over its samples of the facing-distance
rule, then the minimum along every straight run (collinear curves joined at a point of degree
two), because Gmsh terminates a layer of one thickness only at a corner — the end column runs
along the adjacent curve — so the row count may change at corners (fans) but not along a run.
"""
function band_cap_rows_by_curve(metal_edge_curves, radial_um::Float64, rows_max::Int)
    curves = metal_curves(metal_edge_curves)
    rows = Int[]
    for (index, curve) in enumerate(curves)
        length_um = point_distance(curve.a, curve.b)
        n = max(1, ceil(Int, length_um / BAND_CAP_SAMPLE_SPACING_UM))
        at(t) = (
            curve.a[1] + t * (curve.b[1] - curve.a[1]),
            curve.a[2] + t * (curve.b[2] - curve.a[2])
        )
        push!(
            rows,
            minimum(
                band_cap_rows(
                    facing_distance(at(k / n), index, curves),
                    radial_um,
                    rows_max
                ) for k = 0:n
            )
        )
    end
    # Straight runs: union of collinear curve pairs meeting at a point of degree two.
    by_point = Dict{Int32, Vector{Int}}()
    for (index, curve) in enumerate(curves), point in curve.points
        push!(get!(by_point, point, Int[]), index)
    end
    parent = collect(1:length(curves))
    find(i) = parent[i] == i ? i : (parent[i] = find(parent[i]))
    for (_, members) in by_point
        length(members) == 2 || continue
        u, v = curves[members[1]], curves[members[2]]
        du = (u.b[1] - u.a[1], u.b[2] - u.a[2])
        dv = (v.b[1] - v.a[1], v.b[2] - v.a[2])
        cross = du[1] * dv[2] - du[2] * dv[1]
        abs(cross) <= BAND_CAP_COLLINEARITY * hypot(du...) * hypot(dv...) || continue
        parent[find(members[1])] = find(members[2])
    end
    run_rows = Dict{Int, Int}()
    for index in eachindex(curves)
        root = find(index)
        run_rows[root] = min(get(run_rows, root, rows_max), rows[index])
    end
    return Dict{Any, Int}(
        curve_signature(1, curve.tag) => run_rows[find(index)] for
        (index, curve) in enumerate(curves)
    )
end

include("structured_band.jl")

# The region mesh size (Gmsh `Mesh.MeshSizeMax`); also the largest plan size the graded sweep's
# interior stacks are thinned to.
const REGION_MESH_SIZE_MAX_UM = 30.0
# Region size grading of the graded sweep's plan (decision 275): the region mesh size grows from
# the band's station size t at the band tops with this slope (size per unit distance,
# dimensionless) up to REGION_MESH_SIZE_MAX_UM; the Gmsh Distance field's sampling per curve.
const REGION_GRADING_SLOPE = 1.0
const REGION_GRADING_DISTANCE_SAMPLING = 20

# A metal chain of the own band for the graded sweep: its columns as plan indices.
struct GradedChain
    closed::Bool
    class::Int
    plan_rows::Vector{Vector{Int}} # per column: plan index of the base, then of rows 1..n
end

struct PlanMesh
    xy::Vector{Point2}
    triangles::Vector{NTuple{3, Int}}
    triangle_class::Vector{Int} # index into classes
    classes::Vector{PartitionClass}
    band_triangle::BitVector # a triangle of the own band (false: region / Gmsh mode)
    chains::Vector{GradedChain} # the own band's chains (empty in Gmsh mode)
end

function mesh_plan(
    spec::PolygonSet,
    radial_um::Float64,
    tangential_um::Float64;
    radial_growth::Float64=2.0,
    radial_band_um::Float64=1.55,
    exact_band_thickness::Bool=false,
    band_cap::Symbol=:none,
    band_mode::Symbol=:own,
    fan_turn_angle_deg::Float64=DEFAULT_FAN_TURN_ANGLE_DEG,
    region_grading::Bool=false,
    vertex_column_grading::Bool=false,
    vertex_grading_min_turn_deg::Float64=DEFAULT_VERTEX_GRADING_MIN_TURN_DEG,
    verbose::Bool=true
)
    band_cap in (:none, :partition, :curve) ||
        error("band_cap must be :none, :partition or :curve")
    band_mode in (:own, :gmsh) || error("band_mode must be :own or :gmsh")
    region_grading && band_mode != :own && error("region_grading needs band_mode own")
    vertex_column_grading &&
        band_mode != :own &&
        error("vertex_column_grading needs band_mode own")
    0.0 <= fan_turn_angle_deg <= MAX_FAN_TURN_ANGLE_DEG || error(
        "fan_turn_angle_deg must lie in [0, $(MAX_FAN_TURN_ANGLE_DEG)]: above it the mitre " *
        "column 1 / cos(turn / 2) exceeds twice the band"
    )
    band_mode == :gmsh ||
        band_cap == :none ||
        error("band_cap applies to band_mode gmsh only (own mode always applies the rule)")
    radial_layers = max(1, round(Int, log2(1.0 + radial_band_um / radial_um)))
    plan_model = gmsh.model.get_current()
    surfaces, class_by_surface, metal_edge_curves = plan_partitions(spec)
    rows_by_signature = band_cap_rows_by_curve(metal_edge_curves, radial_um, radial_layers)
    metal_signatures = Set(curve_signature(1, tag) for tag in metal_edge_curves)
    class_list = unique(class_by_surface[s] for s in surfaces)
    class_index = Dict(class => i for (i, class) in enumerate(class_list))

    # Own mode: the band is exactly radial_layers rows of h_k = r (2^k - 1) by construction.
    band_thickness(rows) =
        radial_um * (radial_growth^rows - 1.0) / (radial_growth - 1.0) *
        (exact_band_thickness || band_mode == :own ? 1.0 : 1.0 + BAND_THICKNESS_MARGIN)
    radial_thickness = band_thickness(radial_layers)
    # Band-cap record, PER STRAIGHT RUN (the minimum along collinear curves = what Gmsh's
    # `:curve` mode can apply; the own band applies the rule per 1D segment and side, recorded
    # under `band`): metal-edge length per row count (every curve once).
    length_by_rows = Dict{Int, Float64}()
    for curve in metal_curves(metal_edge_curves)
        rows = rows_by_signature[curve_signature(1, curve.tag)]
        length_by_rows[rows] =
            get(length_by_rows, rows, 0.0) + point_distance(curve.a, curve.b)
    end
    capped_length = sum(l for (rows, l) in length_by_rows if rows < radial_layers; init=0.0)
    band_cap_record = Dict{String, Any}(
        "mode" => string(band_cap),
        "band_mode" => string(band_mode),
        "fraction" => BAND_CAP_FRACTION,
        "rows_max" => radial_layers,
        "rule_statistic" =>
            "per straight run (minimum along collinear curves; Gmsh " *
            ":curve mode), both sides; the own band applies the rule per " *
            "1D segment and side: see `band`",
        "rule_capped_length_um" => capped_length,
        "rule_minimum_rows" => minimum(keys(length_by_rows)),
        "rule_minimum_band_um" => band_thickness(minimum(keys(length_by_rows))),
        "rule_length_um_by_rows" => Dict(string(k) => v for (k, v) in length_by_rows)
    )
    for option in (
        "General.NumThreads",
        "Mesh.MaxNumThreads1D",
        "Mesh.MaxNumThreads2D",
        "Mesh.MaxNumThreads3D"
    )
        gmsh.option.set_number(option, 1)
    end
    gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.set_number("Mesh.MeshSizeMin", radial_um)
    gmsh.option.set_number("Mesh.MeshSizeMax", REGION_MESH_SIZE_MAX_UM)
    gmsh.option.set_number("Mesh.BoundaryLayerFanElements", 7)
    gmsh.option.set_number("Mesh.Algorithm", 6)
    verbose && println(
        "Plan partitions: ",
        length(surfaces),
        ", metal-edge curves: ",
        length(metal_edge_curves),
        ", boundary layer first=",
        radial_um,
        " um, layers=",
        radial_layers,
        ", thickness=",
        radial_thickness,
        " um (",
        band_mode == :own ? "own structured band" :
        exact_band_thickness ? "exact geometric sum" : "geometric sum x (1 + 1e-6)",
        "); band cap (per run): ",
        band_cap_record
    )

    raw_triangles = Tuple{NTuple{3, Point2}, Int}[]
    band_flags = BitVector()
    raw_chains = Tuple{Bool, Int, Vector{BandColumn}}[]
    boundary_layer_ends = 0
    partition_rows = fill(radial_layers, length(surfaces))
    band_record = Dict{String, Any}("mode" => string(band_mode))
    if band_mode == :own
        heights = band_heights(radial_um, radial_growth, radial_layers)
        curves = plan_curves(surfaces, metal_edge_curves)
        0.0 <= vertex_grading_min_turn_deg < 180.0 ||
            error("vertex_grading_min_turn_deg must lie in [0, 180)")
        vertex_angles =
            vertex_column_grading ?
            plan_vertex_angles(curves, spec.box, deg2rad(vertex_grading_min_turn_deg)) :
            Dict{Int32, Float64}()
        ladder_stations = 0
        graded_curve_ends = 0
        nodes_by_curve = Dict{Int, Vector{Point2}}()
        for (k, curve) in enumerate(curves)
            curve.metal || continue
            if vertex_column_grading
                phi_a = get(vertex_angles, curve.points[1], NaN)
                phi_b = get(vertex_angles, curve.points[2], NaN)
                nodes_by_curve[k], added = vertex_graded_curve_nodes(
                    curve,
                    tangential_um,
                    radial_um,
                    phi_a,
                    phi_b
                )
                ladder_stations += added
                graded_curve_ends += !isnan(phi_a) + !isnan(phi_b)
            else
                nodes_by_curve[k] = curve_transfinite_nodes(curve, tangential_um)
            end
        end
        statistics = BandStatistics()
        band_triangles = 0
        for (index, surface) in enumerate(surfaces)
            partition = partition_geometry(surface[2], curves, spec.box)
            class = class_index[class_by_surface[surface]]
            band_area, region_area, minimum_rows, triangles = mesh_partition_band(
                partition,
                curves,
                nodes_by_curve,
                heights,
                radial_um,
                deg2rad(fan_turn_angle_deg),
                spec.box,
                class,
                raw_triangles,
                statistics,
                plan_model,
                max(2.0 * heights[end], tangential_um),
                raw_chains,
                region_grading ? tangential_um : 0.0
            )
            append!(band_flags, trues(triangles))
            append!(band_flags, falses(length(raw_triangles) - length(band_flags)))
            partition_rows[index] = minimum_rows
            band_triangles += triangles
            verbose && println(
                "  partition ",
                index,
                "/",
                length(surfaces),
                " surface ",
                surface[2],
                " class ",
                class_list[class].conductors,
                " bump ",
                class_list[class].bump,
                " band area ",
                band_area,
                " region area ",
                region_area,
                " loops ",
                length(partition.loops),
                " minimum rows ",
                minimum_rows
            )
        end
        boundary_layer_ends = statistics.wall_end_columns
        merge!(
            band_record,
            Dict{String, Any}(
                "fan_turn_angle_deg" => fan_turn_angle_deg,
                "fan_turn_tolerance_deg" => FAN_TURN_TOLERANCE_DEG,
                "fan_columns" => FAN_COLUMNS,
                "corner_clearance_fraction" => CORNER_CLEARANCE_FRACTION,
                "rows_max" => radial_layers,
                "heights_um" => heights,
                "band_triangles" => band_triangles,
                "quads" => statistics.quads,
                "collapse_triangles" => statistics.collapse_triangles,
                "fan_triangles" => statistics.fan_triangles,
                "fans" => statistics.fans,
                "mitre_outward_corners" => statistics.mitre_corners,
                "inward_corners" => statistics.inward_corners,
                "clearance_capped_columns" => statistics.clearance_capped_columns,
                "clearance_clamped_columns" => statistics.clearance_clamped_columns,
                "scale_capped_columns" => statistics.scale_capped_columns,
                "scale_clamped_columns" => statistics.scale_clamped_columns,
                "max_mitre_scale" => statistics.max_mitre_scale,
                "max_inward_scale" => statistics.max_inward_scale,
                "max_wall_end_scale" => statistics.max_wall_end_scale,
                "region_area_tolerance" => REGION_AREA_TOLERANCE,
                "region_grading" =>
                    region_grading ?
                    Dict{String, Any}(
                        "size_min_um" => tangential_um,
                        "size_max_um" => REGION_MESH_SIZE_MAX_UM,
                        "slope" => REGION_GRADING_SLOPE,
                        "distance_max_um" =>
                            (REGION_MESH_SIZE_MAX_UM - tangential_um) /
                            REGION_GRADING_SLOPE
                    ) : false,
                "wall_end_columns" => statistics.wall_end_columns,
                "vertex_column_grading" =>
                    vertex_column_grading ?
                    Dict{String, Any}(
                        "min_turn_deg" => vertex_grading_min_turn_deg,
                        "vertices" => length(vertex_angles),
                        "graded_curve_ends" => graded_curve_ends,
                        "ladder_stations" => ladder_stations,
                        "min_wedge_angle_deg" =>
                            isempty(vertex_angles) ? nothing :
                            rad2deg(minimum(values(vertex_angles))),
                        "first_station_rule" => "s_1 = r (2^j - 1) >= 2 r / tan(phi / 2)"
                    ) : false,
                "segments_below_2p5r" => statistics.segments_below_2p5r,
                "quad_min_abs_sin" => statistics.quad_min_abs_sin,
                "triangle_min_angle_deg" => statistics.triangle_min_angle_deg,
                "per_side_segment_length_um_by_rows" =>
                    Dict(string(k) => v for (k, v) in statistics.length_um_by_rows),
                "per_side_capped_length_um" => sum(
                    l for (rows, l) in statistics.length_um_by_rows if rows < radial_layers;
                    init=0.0
                )
            )
        )
        verbose && println("Own band: ", band_record)
    else
        # One-sided boundary layers need every source curve to bound a single surface: copy each
        # partition and mesh the copies one at a time (coincident copies never share a Delaunay
        # problem); the copies' interface nodes are welded by coordinate afterwards.
        copies = Tuple{Int32, Int32}[]
        copy_class = Int[]
        for surface in surfaces
            append!(copies, gmsh.model.occ.copy([surface]))
            push!(copy_class, class_index[class_by_surface[surface]])
        end
        gmsh.model.occ.synchronize()
        gmsh.model.set_visibility(gmsh.model.get_entities(), 0, true)
        gmsh.option.set_number("Mesh.MeshOnlyVisible", 1)
        copy_curves = unique(gmsh.model.get_boundary(copies, false, false, false))
        metal_copy_curves = Int32[]
        for (_, tag) in copy_curves
            if curve_signature(1, tag) in metal_signatures
                push!(metal_copy_curves, abs(tag))
                length_um = gmsh.model.occ.get_mass(1, abs(tag))
                segments = max(1, ceil(Int, length_um / tangential_um - 1.0e-6))
                gmsh.model.mesh.set_transfinite_curve(abs(tag), segments + 1)
            end
        end
        metal_copy_set = Set(metal_copy_curves)
        length(metal_copy_curves) == 2 * length(metal_edge_curves) || error(
            "Metal-edge curve copies: $(length(metal_copy_curves)) for " *
            "$(length(metal_edge_curves)) curves"
        )
        for (copy_index, surface) in enumerate(copies)
            gmsh.model.set_visibility(copies, 0, true)
            gmsh.model.set_visibility([surface], 1, true)
            sources = [
                Float64(abs(tag)) for
                (_, tag) in gmsh.model.get_boundary([surface], false, false, false) if
                abs(tag) in metal_copy_set
            ]
            # One BoundaryLayer field per row count (the band cap), each with its layer end points
            # (a metal edge ending on the window wall, a junction with another row count).
            fields = Int32[]
            if !isempty(sources)
                ends = boundary_layer_end_points(sources)
                boundary_layer_ends += length(ends)
                curve_rows =
                    [rows_by_signature[curve_signature(1, Int32(tag))] for tag in sources]
                if band_cap == :none
                    fill!(curve_rows, radial_layers)
                elseif band_cap == :partition
                    fill!(curve_rows, minimum(curve_rows))
                end
                by_rows = Dict{Int, Vector{Float64}}()
                for (tag, rows) in zip(sources, curve_rows)
                    push!(get!(by_rows, rows, Float64[]), tag)
                end
                partition_rows[copy_index] = minimum(curve_rows)
                for rows in sort!(collect(keys(by_rows)))
                    group = by_rows[rows]
                    field = gmsh.model.mesh.field.add("BoundaryLayer")
                    push!(fields, field)
                    gmsh.model.mesh.field.set_numbers(field, "CurvesList", group)
                    # The field's own layer ends: the wall ends and the junctions with the
                    # curves of another row count (Gmsh needs both declared).
                    group_ends = boundary_layer_end_points(group)
                    isempty(group_ends) ||
                        gmsh.model.mesh.field.set_numbers(field, "PointsList", group_ends)
                    gmsh.model.mesh.field.set_number(field, "Size", radial_um)
                    gmsh.model.mesh.field.set_number(field, "Ratio", radial_growth)
                    gmsh.model.mesh.field.set_number(
                        field,
                        "Thickness",
                        band_thickness(rows)
                    )
                    gmsh.model.mesh.field.set_number(field, "Quads", 1)
                    gmsh.model.mesh.field.set_number(field, "IntersectMetrics", 1)
                    gmsh.model.mesh.field.set_as_boundary_layer(field)
                end
            end
            gmsh.model.mesh.generate(2)
            node_tags, node_coordinates, _ = gmsh.model.mesh.get_nodes()
            coordinate_by_tag = Dict{UInt64, Point2}()
            for (index, tag) in enumerate(node_tags)
                abs(node_coordinates[3index]) <= 1.0e-6 || error("Plan node is not on z=0")
                coordinate_by_tag[tag] =
                    (node_coordinates[3index - 2], node_coordinates[3index - 1])
            end
            area = 0.0
            element_types, _, nodes_by_type = gmsh.model.mesh.get_elements(2, surface[2])
            for (element_type, element_nodes) in zip(element_types, nodes_by_type)
                _, dimension, _, node_count, _, primary =
                    gmsh.model.mesh.get_element_properties(element_type)
                dimension == 2 || continue
                primary in (3, 4) ||
                    error("Expected triangular or quadrilateral plan elements")
                for offset = 0:node_count:(length(element_nodes) - node_count)
                    corners =
                        [coordinate_by_tag[element_nodes[offset + i]] for i = 1:primary]
                    split =
                        primary == 3 ? ((corners[1], corners[2], corners[3]),) :
                        (
                            (corners[1], corners[2], corners[3]),
                            (corners[1], corners[3], corners[4])
                        )
                    for triangle in split
                        area += triangle_area(triangle...)
                        push!(raw_triangles, (triangle, copy_class[copy_index]))
                    end
                end
            end
            occ_area = gmsh.model.occ.get_mass(surface...)
            abs(area - occ_area) <= 1.0e-6 * max(1.0, occ_area) || error(
                "Partition $(surface[2]) mesh area $area differs from OCC area $occ_area"
            )
            verbose && println(
                "  partition ",
                copy_index,
                "/",
                length(copies),
                " surface ",
                surface[2],
                " class ",
                class_list[copy_class[copy_index]].conductors,
                " bump ",
                class_list[copy_class[copy_index]].bump,
                " area ",
                area,
                " sources ",
                length(sources)
            )
            for field in fields
                gmsh.model.mesh.field.remove(field)
            end
            gmsh.model.mesh.clear(copies)
        end
    end # band_mode
    isempty(raw_triangles) && error("No plan triangles extracted")
    band_mode == :own || append!(band_flags, falses(length(raw_triangles)))
    length(band_flags) == length(raw_triangles) || error("Band flags out of step")
    band_cap_record["applied_rows_by_partition"] = partition_rows
    band_cap_record["applied_minimum_rows"] = minimum(partition_rows)
    band_record["applied_rows_by_partition"] = partition_rows
    band_record["applied_minimum_rows"] = minimum(partition_rows)

    merge_tolerance = PLAN_NODE_MERGE_TOLERANCE_UM
    coordinate_index = Dict{NTuple{2, Int64}, Int}()
    xy = Point2[]
    function node_index(p::Point2)
        key = (round(Int64, p[1] / merge_tolerance), round(Int64, p[2] / merge_tolerance))
        index = get(coordinate_index, key, 0)
        if index == 0
            push!(xy, p)
            index = length(xy)
            coordinate_index[key] = index
        else
            q = xy[index]
            hypot(p[1] - q[1], p[2] - q[2]) <= 2merge_tolerance ||
                error("Plan-node merge tolerance collision")
        end
        return index
    end
    triangles = NTuple{3, Int}[]
    triangle_class = Int[]
    for (triangle, class) in raw_triangles
        push!(
            triangles,
            (node_index(triangle[1]), node_index(triangle[2]), node_index(triangle[3]))
        )
        push!(triangle_class, class)
    end
    # The own band's chains as plan indices (every column node is a band triangle corner).
    function welded_index(p::Point2)
        key = (round(Int64, p[1] / merge_tolerance), round(Int64, p[2] / merge_tolerance))
        index = get(coordinate_index, key, 0)
        index > 0 || error("Band column node $p is not a plan node")
        return index
    end
    chains = [
        GradedChain(
            closed,
            class,
            [
                [welded_index(column.base); welded_index.(column.nodes)] for
                column in columns
            ]
        ) for (closed, class, columns) in raw_chains
    ]
    return PlanMesh(xy, triangles, triangle_class, class_list, band_flags, chains),
    radial_layers,
    radial_thickness,
    boundary_layer_ends,
    band_cap_record,
    band_record
end

# Points where a metal edge of the partition ends against a non-metal boundary curve (the
# window wall cutting a conductor): the end of exactly one source curve. Declared to the
# BoundaryLayer field as `PointsList`, Gmsh terminates the layer there with a column of the
# layer's own heights along the wall; undeclared, it reuses the wall's far-field 1D nodes for
# that column and produces inverted quads tens of micrometres long. A closed metal loop (the
# transmon, the synthetic window) has no such point.
function boundary_layer_end_points(sources::Vector{Float64})
    counts = Dict{Int32, Int}()
    for tag in sources,
        (_, point) in gmsh.model.get_boundary([(1, Int32(tag))], false, false, false)

        counts[abs(point)] = get(counts, abs(point), 0) + 1
    end
    return sort!([Float64(point) for (point, count) in counts if count == 1])
end

# Plan edge -> owning triangles; metal perimeter edges per plane / bump with their conductor.
struct PlanTopology
    edge_incidence::Dict{Tuple{Int, Int}, Vector{Int}}
    open_edges::Int
    plane_edge_conductor::Vector{Dict{Tuple{Int, Int}, String}} # per plane
    bump_edge_conductor::Dict{Tuple{Int, Int}, String}
end

function plan_topology(spec::PolygonSet, plan::PlanMesh)
    edge_incidence = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (index, t) in enumerate(plan.triangles)
        for (a, b) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
            push!(get!(edge_incidence, a < b ? (a, b) : (b, a), Int[]), index)
        end
    end
    any(owners -> length(owners) > 2, values(edge_incidence)) &&
        error("Plan mesh has nonmanifold edges")
    box = spec.box
    function on_wall(edge)
        a, b = plan.xy[edge[1]], plan.xy[edge[2]]
        tolerance = 1.0e-5
        return (abs(a[1] - box[1]) <= tolerance && abs(b[1] - box[1]) <= tolerance) ||
               (abs(a[1] - box[2]) <= tolerance && abs(b[1] - box[2]) <= tolerance) ||
               (abs(a[2] - box[3]) <= tolerance && abs(b[2] - box[3]) <= tolerance) ||
               (abs(a[2] - box[4]) <= tolerance && abs(b[2] - box[4]) <= tolerance)
    end
    open_edges = [edge for (edge, owners) in edge_incidence if length(owners) == 1]
    unexpected = count(!on_wall, open_edges)
    unexpected == 0 || error("Plan mesh has $unexpected open edges off the box wall")
    plane_edge_conductor = [Dict{Tuple{Int, Int}, String}() for _ in spec.planes]
    bump_edge_conductor = Dict{Tuple{Int, Int}, String}()
    for (edge, owners) in edge_incidence
        length(owners) == 2 || continue
        c1 = plan.classes[plan.triangle_class[owners[1]]]
        c2 = plan.classes[plan.triangle_class[owners[2]]]
        for k in eachindex(spec.planes)
            a, b = c1.conductors[k], c2.conductors[k]
            a == b && continue
            (isempty(a) || isempty(b)) ||
                error("Conductors $a and $b touch on plane $(spec.planes[k].name)")
            plane_edge_conductor[k][edge] = isempty(a) ? b : a
        end
        if c1.bump != c2.bump
            (c1.bump == 0 || c2.bump == 0) || error("Bumps touch")
            bump_edge_conductor[edge] = spec.bumps[max(c1.bump, c2.bump)].conductor
        end
    end
    return PlanTopology(
        edge_incidence,
        length(open_edges),
        plane_edge_conductor,
        bump_edge_conductor
    )
end

function edge_metrics(plan::PlanMesh, topology::PlanTopology, radial_um, tangential_um)
    metal_edges = Set{Tuple{Int, Int}}()
    for edges in topology.plane_edge_conductor
        union!(metal_edges, keys(edges))
    end
    union!(metal_edges, keys(topology.bump_edge_conductor))
    tangents = Float64[]
    heights = Float64[]
    outliers = 0
    for edge in metal_edges
        a, b = plan.xy[edge[1]], plan.xy[edge[2]]
        tangent = hypot(b[1] - a[1], b[2] - a[2])
        push!(tangents, tangent)
        for index in topology.edge_incidence[edge]
            t = plan.triangles[index]
            third = only(setdiff((t[1], t[2], t[3]), edge))
            height = 2 * triangle_area(a, b, plan.xy[third]) / tangent
            push!(heights, height)
            height > 1.5radial_um && (outliers += 1)
        end
    end
    tangent_summary = metric_summary(tangents)
    height_summary = metric_summary(heights)
    tangent_summary["maximum"] <= 1.01tangential_um ||
        error("Metal-edge tangent target was not enforced: $tangent_summary")
    height_summary["maximum"] <= 1.1radial_um ||
        error("Metal-edge first-layer normal target was not enforced: $height_summary")
    height_summary["minimum"] >= 0.9radial_um || error(
        "Metal-edge first-layer normal height below 0.9 x the target (compressed fronts): " *
        "$height_summary"
    )
    return length(metal_edges), tangent_summary, height_summary, outliers
end

# Connected components of every conductor label per plane (a diagnostic: an island split in
# two pieces or a ground split by the window is reported, not refused).
function conductor_components(spec::PolygonSet, plan::PlanMesh, topology::PlanTopology)
    result = Dict{String, Any}()
    for (k, plane) in enumerate(spec.planes)
        parent = collect(1:length(plan.triangles))
        function root(i)
            while parent[i] != i
                parent[i] = parent[parent[i]]
                i = parent[i]
            end
            return i
        end
        for owners in values(topology.edge_incidence)
            length(owners) == 2 || continue
            a = plan.classes[plan.triangle_class[owners[1]]].conductors[k]
            b = plan.classes[plan.triangle_class[owners[2]]].conductors[k]
            (a == b && !isempty(a)) && (parent[root(owners[1])] = root(owners[2]))
        end
        areas = Dict{String, Dict{Int, Float64}}()
        for (index, t) in enumerate(plan.triangles)
            label = plan.classes[plan.triangle_class[index]].conductors[k]
            isempty(label) && continue
            component = get!(areas, label, Dict{Int, Float64}())
            component[root(index)] =
                get(component, root(index), 0.0) +
                triangle_area(plan.xy[t[1]], plan.xy[t[2]], plan.xy[t[3]])
        end
        result[plane.name] = Dict(
            label => Dict(
                "components" => length(components),
                "area_um2" => sum(values(components)),
                "component_areas_um2" => sort!(collect(values(components)); rev=true)
            ) for (label, components) in areas
        )
    end
    return result
end

# ---------------------------------------------------------------------------------------------

include("graded_sweep.jl")

signed_volume(p1, p2, p3, p4) =
    (
        (p2[1] - p1[1]) *
        ((p3[2] - p1[2]) * (p4[3] - p1[3]) - (p3[3] - p1[3]) * (p4[2] - p1[2])) -
        (p2[2] - p1[2]) *
        ((p3[1] - p1[1]) * (p4[3] - p1[3]) - (p3[3] - p1[3]) * (p4[1] - p1[1])) +
        (p2[3] - p1[3]) *
        ((p3[1] - p1[1]) * (p4[2] - p1[2]) - (p3[2] - p1[2]) * (p4[1] - p1[1]))
    ) / 6

function face_area(p1, p2, p3)
    ux, uy, uz = p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]
    vx, vy, vz = p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3]
    return 0.5 * hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
end

"""
    mesh_polygon_window(spec, radial_um, tangential_um, output; verbose=true,
                        plan_only=false, band_mode=:own, fan_turn_angle_deg=90.0,
                        exact_band_thickness=false, cross_plane_snap_um=NaN,
                        band_cap=:none, sweep=:tensor, alpha=1.0, beta=3.0,
                        region_grading=(sweep == :graded), region_ring=true,
                        region_z_grading=false, step_face_z_grading=:legacy,
                        vertex_column_grading=false,
                        vertex_grading_min_turn_deg=30.0) -> manifest

Generate the fabricated reference mesh of a polygon set and write `output` (ASCII MSH2) with
its JSON manifest next to it. `band_mode=:own` (default) builds the structured boundary-layer
band itself with the per-segment, per-side band cap (structured_band.jl); `:gmsh` uses Gmsh's
BoundaryLayer field as the recorded transmon generator did. `fan_turn_angle_deg`: outward
corners turning more than this (by more than `FAN_TURN_TOLERANCE_DEG`) get a fan of columns,
the others the scaled bisector column (90: the recorded corner treatment). `sweep=:tensor`
(default) sweeps every plan triangle through every z level (the recorded family);
`sweep=:graded` (own band only) sweeps the graded (n, z) cross-section of graded_sweep.jl along
the band columns with `alpha` (row k's z spacing >= alpha x its width) and `beta` (row k ends
beta h_k beyond the fabricated step), dimensionless. Graded sweep only (decision 275, both on by
default): `region_grading` grades the Gmsh region's plan size from the band's station size t at
the band tops with the slope `REGION_GRADING_SLOPE` up to `REGION_MESH_SIZE_MAX_UM` (the tensor
sweep's plan is the recorded family's and cannot be graded); `region_ring` puts the region
nodes adjacent to the band on the band's outermost stack Z_K instead of their plan-size ladder
step; `region_z_grading` (decision 276 option (ii)) puts every region node on the geometric
stack grown from the fabricated steps' faces (`region_z_graded_levels`). Both sweeps
(reference-quality lane, decision 406; both default to the recorded family):
`step_face_z_grading=:mirrored` (D1) adds the band row heights r (2^k - 1) as z levels beyond
both faces of every plane's fabricated step (`step_face_z_levels`; `:legacy` leaves the first
cell beyond each face at 50 nm); `vertex_column_grading=true` (D3, own band only) grades the
band's column stations along every metal curve toward its plan vertices (joints turning by at
least `vertex_grading_min_turn_deg`, junctions always) with the same r (2^k - 1) ladder
(`vertex_graded_curve_nodes`). Gmsh mode only: `exact_band_thickness=true` passes
the exact geometric sum as the boundary-layer Thickness (the recorded transmon generator's
formula; see the header); `band_cap` applies the local band cap rule (`:none`: the rule is only
recorded; `:partition`: one band per partition, the minimum over its curves; `:curve`: one
BoundaryLayer field per row count — Gmsh fails on many of its transitions).
`cross_plane_snap_um` overrides the cross-plane snap distance 0.05 x MatchingRadius of a
two-plane set (see `reconcile_planes`).
"""
function mesh_polygon_window(
    spec::PolygonSet,
    radial_um::Float64,
    tangential_um::Float64,
    output::AbstractString;
    verbose::Bool=true,
    plan_only::Bool=false,
    exact_band_thickness::Bool=false,
    cross_plane_snap_um::Float64=NaN,
    band_cap::Symbol=:none,
    band_mode::Symbol=:own,
    fan_turn_angle_deg::Float64=DEFAULT_FAN_TURN_ANGLE_DEG,
    sweep::Symbol=:tensor,
    alpha::Float64=DEFAULT_GRADED_ALPHA,
    beta::Float64=DEFAULT_GRADED_BETA,
    region_grading::Bool=(sweep == :graded),
    region_ring::Bool=true,
    region_z_grading::Bool=false,
    step_face_z_grading::Symbol=:legacy,
    vertex_column_grading::Bool=false,
    vertex_grading_min_turn_deg::Float64=DEFAULT_VERTEX_GRADING_MIN_TURN_DEG
)
    output = abspath(output)
    sweep in (:tensor, :graded) || error("sweep must be :tensor or :graded")
    step_face_z_grading in STEP_FACE_Z_GRADINGS ||
        error("step_face_z_grading must be one of $(STEP_FACE_Z_GRADINGS)")
    !vertex_column_grading ||
        band_mode == :own ||
        error("vertex_column_grading needs the own structured band (band_mode own)")
    sweep == :tensor ||
        band_mode == :own ||
        error("The graded sweep needs the own structured band (band_mode own)")
    sweep == :graded ||
        !region_grading ||
        error(
            "region_grading belongs to the graded sweep (the tensor sweep's plan is fixed)"
        )
    radial_growth = 2.0
    snap_delta, snap_rule = cross_plane_snap_distance(spec, cross_plane_snap_um)
    spec, reconciliation = reconcile_planes(spec, snap_delta)
    reconciliation["rule"] = snap_rule
    reconciliation["matching_radius_um"] =
        isnan(spec.matching_radius) ? nothing : spec.matching_radius
    verbose &&
        reconciliation["applied"] &&
        println("Cross-plane reconciliation: ", reconciliation)
    metal_layers = max(2, ceil(Int, spec.metal_thickness / radial_um - 1.0e-9))
    trench_layers = max(1, ceil(Int, spec.overetch / radial_um - 1.0e-9))
    gmsh.initialize()
    plan,
    radial_layers,
    radial_thickness,
    boundary_layer_ends,
    band_cap_record,
    band_record = try
        gmsh.option.set_number("General.Verbosity", 2)
        gmsh.model.add(spec.name)
        mesh_plan(
            spec,
            radial_um,
            tangential_um;
            exact_band_thickness=exact_band_thickness,
            band_cap=band_cap,
            band_mode=band_mode,
            fan_turn_angle_deg=fan_turn_angle_deg,
            region_grading=region_grading,
            vertex_column_grading=vertex_column_grading,
            vertex_grading_min_turn_deg=vertex_grading_min_turn_deg,
            verbose=verbose
        )
    finally
        gmsh.finalize()
    end
    topology = plan_topology(spec, plan)
    box = spec.box
    plan_area = sum(
        triangle_area(plan.xy[t[1]], plan.xy[t[2]], plan.xy[t[3]]) for t in plan.triangles
    )
    expected_area = (box[2] - box[1]) * (box[4] - box[3])
    abs(plan_area - expected_area) / expected_area < 1.0e-8 ||
        error("Plan area mismatch: $plan_area vs $expected_area")
    metal_edge_count, tangent_summary, height_summary, height_outliers =
        edge_metrics(plan, topology, radial_um, tangential_um)
    components = conductor_components(spec, plan, topology)
    verbose && println(
        "Plan nodes: ",
        length(plan.xy),
        ", triangles: ",
        length(plan.triangles),
        ", metal perimeter edges: ",
        metal_edge_count,
        "\nTangent (um): ",
        tangent_summary,
        "\nFirst-layer normal height (um): ",
        height_summary,
        ", outliers > 1.5 r: ",
        height_outliers,
        "\nConductor components: ",
        components
    )
    step_face_levels =
        step_face_z_grading == :mirrored ? step_face_z_levels(spec, radial_um) : Float64[]
    stack = z_levels(spec, metal_layers, trench_layers; step_face_levels=step_face_levels)
    zs = stack.levels
    manifest = Dict{String, Any}(
        "mesh" => output,
        "polygon_set" => spec.name,
        "radial_target_um" => radial_um,
        "tangential_target_um" => tangential_um,
        "radial_growth" => radial_growth,
        "radial_layers" => radial_layers,
        "radial_band_thickness_um" => radial_thickness,
        "radial_band_thickness_mode" =>
            band_mode == :own ? "own_structured_band" :
            exact_band_thickness ? "exact_geometric_sum" : "geometric_sum_x_1p000001",
        "metal_thickness_um" => spec.metal_thickness,
        "overetch_um" => spec.overetch,
        "metal_layers" => metal_layers,
        "trench_layers" => trench_layers,
        "planes" => [
            Dict(
                "name" => p.name,
                "surface_z_um" => p.surface_z,
                "facing" => p.facing == 1 ? "up" : "down",
                "substrate_thickness_um" => p.substrate_thickness
            ) for p in spec.planes
        ],
        "bumps" => length(spec.bumps),
        "plan_nodes" => length(plan.xy),
        "plan_triangles" => length(plan.triangles),
        "plan_partition_classes" => length(plan.classes),
        "cross_plane_reconciliation" => reconciliation,
        "plan_open_outer_edges" => topology.open_edges,
        "plan_perimeter_edges" => metal_edge_count,
        "plan_boundary_layer_end_points" => boundary_layer_ends,
        "band_cap" => band_cap_record,
        "band" => band_record,
        "perimeter_tangent_length_um" => tangent_summary,
        "first_layer_normal_height_um" => height_summary,
        "first_layer_normal_outliers_above_1p5x_target" => height_outliers,
        "conductor_components" => components,
        "z_bounds_um" => Dict("bottom" => stack.z_bottom, "top" => stack.z_top),
        "substrate_backsides_z_um" => stack.backsides,
        "z_levels" => zs,
        "sweep" => string(sweep),
        "step_face_z_grading" => string(step_face_z_grading),
        "step_face_z_levels_um" => step_face_levels,
        "attributes" => Dict(
            name => Dict("dimension" => d, "attribute" => a) for
            (d, a, name) in physical_names(spec)
        )
    )
    plan_only && return manifest

    # Nodes: every plan node at every z level, (level - 1) x plan nodes + plan index; the
    # unused ones (inside the metal, thinned away by the graded sweep) are compacted below.
    n_plan = length(plan.xy)
    n_levels = length(zs)
    nodes = Vector{NTuple{3, Float64}}(undef, n_plan * n_levels)
    for level = 1:n_levels, index = 1:n_plan
        nodes[(level - 1) * n_plan + index] =
            (plan.xy[index][1], plan.xy[index][2], zs[level])
    end
    tetrahedra = NTuple{4, Int32}[]
    tetrahedron_attribute = Int8[]
    tetrahedron_class = Int16[] # partition class per tetrahedron (horizontal face tags)
    if sweep == :graded
        tetrahedra, tetrahedron_attribute, tetrahedron_class, graded_record =
            graded_sweep_elements(
                spec,
                plan,
                topology,
                stack,
                radial_um,
                radial_growth,
                radial_layers,
                alpha,
                beta;
                region_ring=region_ring,
                region_z_grading=region_z_grading,
                verbose=verbose
            )
        manifest["graded_sweep"] = graded_record
    else
        # Tensor sweep: every plan triangle x every z interval is a prism of one material (or
        # excluded).
        material_table = [
            material(spec, class, 0.5 * (zs[i] + zs[i + 1])) for
            class in plan.classes, i = 1:(n_levels - 1)
        ]
        sizehint!(tetrahedra, length(plan.triangles) * (n_levels - 1) * 3)
        function push_prism!(a1, a2, a3, b1, b2, b3, attribute, class)
            # Conforming split: the diagonal choice on every quad side follows global node order.
            if a2 < a1 && a2 <= a3
                a1, a2, a3 = a2, a3, a1
                b1, b2, b3 = b2, b3, b1
            elseif a3 < a1 && a3 <= a2
                a1, a2, a3 = a3, a1, a2
                b1, b2, b3 = b3, b1, b2
            end
            if a2 < a3
                push!(tetrahedra, (a1, a2, a3, b3), (a1, a2, b3, b2), (a1, b2, b3, b1))
            else
                push!(tetrahedra, (a1, a2, a3, b2), (a1, a3, b3, b2), (a1, b2, b3, b1))
            end
            push!(tetrahedron_attribute, attribute, attribute, attribute)
            push!(tetrahedron_class, class, class, class)
            return nothing
        end
        for (index, t) in enumerate(plan.triangles)
            class = plan.triangle_class[index]
            for interval = 1:(n_levels - 1)
                attribute = material_table[class, interval]
                attribute == 0 && continue
                lower = Int32((interval - 1) * n_plan)
                upper = Int32(interval * n_plan)
                push_prism!(
                    lower + t[1],
                    lower + t[2],
                    lower + t[3],
                    upper + t[1],
                    upper + t[2],
                    upper + t[3],
                    attribute,
                    Int16(class)
                )
            end
        end
    end
    for index in eachindex(tetrahedra)
        t = tetrahedra[index]
        volume = signed_volume(nodes[t[1]], nodes[t[2]], nodes[t[3]], nodes[t[4]])
        if volume < 0
            tetrahedra[index] = (t[1], t[3], t[2], t[4])
            volume = -volume
        end
        volume > 0 || error("Degenerate tetrahedron $index")
    end
    verbose && println("z levels: ", n_levels, ", tetrahedra: ", length(tetrahedra))

    # Faces: count of owning tetrahedra, the sum of their attributes and the first owner.
    face_data = Dict{NTuple{3, Int32}, Tuple{Int8, Int8, Int32}}()
    sizehint!(face_data, 2 * length(tetrahedra) + length(plan.triangles) * 4)
    for (index, t) in enumerate(tetrahedra)
        attribute = tetrahedron_attribute[index]
        for face in
            ((t[1], t[2], t[3]), (t[1], t[2], t[4]), (t[1], t[3], t[4]), (t[2], t[3], t[4]))
            a, b, c = face
            a > b && ((a, b) = (b, a))
            b > c && ((b, c) = (c, b))
            a > b && ((a, b) = (b, a))
            previous = get(face_data, (a, b, c), (Int8(0), Int8(0), Int32(index)))
            face_data[(a, b, c)] =
                (previous[1] + Int8(1), previous[2] + attribute, previous[3])
        end
    end

    table = attribute_table(spec)
    plane_by_z = Dict{Float64, Tuple{Int, Symbol}}() # metal top / metal bottom z per plane
    for (k, plane) in enumerate(spec.planes)
        plane_by_z[plane.surface_z] = (k, :substrate)
        plane_by_z[plane.surface_z + plane.facing * spec.metal_thickness] = (k, :air)
    end
    level_of_z = Dict(z => i for (i, z) in enumerate(zs))
    base_index(node) = mod(Int(node) - 1, n_plan) + 1
    level_index(node) = (Int(node) - 1) ÷ n_plan + 1
    lower_metal_top =
        spec.planes[1].surface_z + spec.planes[1].facing * spec.metal_thickness
    upper_metal_bottom =
        length(spec.planes) == 2 ?
        spec.planes[2].surface_z + spec.planes[2].facing * spec.metal_thickness : NaN

    surface_elements = Tuple{Int, NTuple{3, Int32}}[]
    for (face, (count, attribute_sum, owner)) in face_data
        if count == 1
            levels = (level_index(face[1]), level_index(face[2]), level_index(face[3]))
            coordinates = (nodes[face[1]], nodes[face[2]], nodes[face[3]])
            attribute = 0
            if levels[1] == levels[2] == levels[3]
                z = zs[levels[1]]
                if levels[1] == 1 || levels[1] == n_levels
                    attribute = 3
                else
                    haskey(plane_by_z, z) || error("Unclassified horizontal face at z=$z")
                    k, side = plane_by_z[z]
                    # The owning tetrahedron's partition (the metal face of its plan class).
                    class = plan.classes[tetrahedron_class[owner]]
                    label = class.conductors[k]
                    isempty(label) &&
                        error("Horizontal metal face without conductor at z=$z")
                    attribute = side == :air ? table[label][1] : table[label][2]
                end
            else
                plan_nodes = unique(base_index.(collect(face)))
                if all(abs(p[1] - box[1]) <= 1.0e-6 for p in coordinates) ||
                   all(abs(p[1] - box[2]) <= 1.0e-6 for p in coordinates) ||
                   all(abs(p[2] - box[3]) <= 1.0e-6 for p in coordinates) ||
                   all(abs(p[2] - box[4]) <= 1.0e-6 for p in coordinates)
                    # A box wall face (incl. the graded sweep's cross-section triangles at a
                    # wall-end column, which span up to three plan nodes and two levels).
                    attribute = 3
                else
                    length(plan_nodes) == 2 || error("Non-vertical boundary face")
                    edge = Tuple(sort(plan_nodes))
                    zmid = sum(p[3] for p in coordinates) / 3
                    label = ""
                    for (k, plane) in enumerate(spec.planes)
                        d = plane.facing * (zmid - plane.surface_z)
                        if 0.0 < d < spec.metal_thickness
                            label = get(topology.plane_edge_conductor[k], edge, "")
                            isempty(label) &&
                                error("Sidewall face off the metal perimeter at z=$zmid")
                        end
                    end
                    if isempty(label) && lower_metal_top < zmid < upper_metal_bottom
                        label = get(topology.bump_edge_conductor, edge, "")
                        isempty(label) &&
                            error("Bump sidewall face off a bump perimeter at z=$zmid")
                    end
                    isempty(label) && error("Unclassified vertical face at z=$zmid")
                    attribute = table[label][1]
                end
            end
            push!(surface_elements, (attribute, face))
        elseif count == 2 && attribute_sum == 3
            levels = (level_index(face[1]), level_index(face[2]), level_index(face[3]))
            z = levels[1] == levels[2] == levels[3] ? zs[levels[1]] : NaN
            attribute = any(b -> abs(z - b) <= 1.0e-9, stack.backsides) ? 9 : 6
            push!(surface_elements, (attribute, face))
        end
    end
    face_data = nothing
    GC.gc()

    # Compact the nodes inside the excluded metal.
    uncompacted = length(nodes)
    used = falses(length(nodes))
    for t in tetrahedra, node in t
        used[node] = true
    end
    remap = zeros(Int32, length(nodes))
    compacted = NTuple{3, Float64}[]
    sizehint!(compacted, count(used))
    for index in eachindex(nodes)
        if used[index]
            push!(compacted, nodes[index])
            remap[index] = Int32(length(compacted))
        end
    end
    tetrahedra = [(remap[t[1]], remap[t[2]], remap[t[3]], remap[t[4]]) for t in tetrahedra]
    surface_elements = [
        (attribute, (remap[f[1]], remap[f[2]], remap[f[3]])) for
        (attribute, f) in surface_elements
    ]
    nodes = compacted

    names = physical_names(spec)
    valid_attributes = Set(a for (d, a, _) in names if d == 2)
    all(element -> element[1] in valid_attributes, surface_elements) ||
        error("Surface element with an attribute outside the table")
    surface_counts = Dict{String, Int}()
    surface_areas = Dict{String, Float64}()
    for (attribute, f) in surface_elements
        key = string(attribute)
        surface_counts[key] = get(surface_counts, key, 0) + 1
        surface_areas[key] =
            get(surface_areas, key, 0.0) + face_area(nodes[f[1]], nodes[f[2]], nodes[f[3]])
    end
    volume_counts = Dict{String, Int}()
    volumes = Dict{String, Float64}()
    for (index, t) in enumerate(tetrahedra)
        key = string(tetrahedron_attribute[index])
        volume_counts[key] = get(volume_counts, key, 0) + 1
        volumes[key] =
            get(volumes, key, 0.0) +
            signed_volume(nodes[t[1]], nodes[t[2]], nodes[t[3]], nodes[t[4]])
    end
    for (d, a, name) in names
        d == 2 &&
            a != 9 &&
            !haskey(surface_counts, string(a)) &&
            error("Physical surface $name (attribute $a) is empty")
    end
    if sweep == :graded
        # The 3D backstop on the written mesh: positive tetrahedra (above) whose volume per
        # material equals the plan's analytic volume, so no cell overlaps or is missing.
        expected_volume = manifest["graded_sweep"]["expected_volume_um3"]
        keys(volumes) == keys(expected_volume) ||
            error("Graded sweep: materials $(keys(volumes)) vs $(keys(expected_volume))")
        for (key, expected) in expected_volume
            abs(volumes[key] - expected) <= VOLUME_CLOSURE_TOLERANCE * expected || error(
                "Graded sweep: the written volume of material $key, $(volumes[key]) um^3, " *
                "differs from the plan's $expected um^3"
            )
        end
    end

    mkpath(dirname(output))
    open(output, "w") do stream
        print(stream, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n")
        written_names = filter(
            entry -> entry[1] == 3 || haskey(surface_counts, string(entry[2])),
            names
        )
        print(stream, "\$PhysicalNames\n$(length(written_names))\n")
        for (dimension, tag, name) in written_names
            print(stream, "$dimension $tag \"$name\"\n")
        end
        print(stream, "\$EndPhysicalNames\n")
        print(stream, "\$Nodes\n$(length(nodes))\n")
        for (index, point) in enumerate(nodes)
            @printf(stream, "%d %.16g %.16g %.16g\n", index, point[1], point[2], point[3])
        end
        print(stream, "\$EndNodes\n")
        print(stream, "\$Elements\n$(length(surface_elements) + length(tetrahedra))\n")
        element = 0
        for (attribute, f) in surface_elements
            element += 1
            print(stream, "$element 2 2 $attribute $attribute $(f[1]) $(f[2]) $(f[3])\n")
        end
        for (index, t) in enumerate(tetrahedra)
            element += 1
            attribute = Int(tetrahedron_attribute[index])
            print(
                stream,
                "$element 4 2 $attribute $attribute $(t[1]) $(t[2]) $(t[3]) $(t[4])\n"
            )
        end
        return print(stream, "\$EndElements\n")
    end
    merge!(
        manifest,
        Dict(
            "bytes" => filesize(output),
            "sha256" => bytes2hex(open(sha256, output)),
            "uncompacted_nodes" => uncompacted,
            "removed_unused_nodes" => uncompacted - length(nodes),
            "nodes" => length(nodes),
            "tetrahedra" => length(tetrahedra),
            "nonpositive" => 0,
            "volume_attribute_counts" => volume_counts,
            "volume_um3" => volumes,
            "surface_triangles" => length(surface_elements),
            "surface_attribute_counts" => surface_counts,
            "surface_area_um2" => surface_areas
        )
    )
    open(replace(output, r"\.msh2$" => ".json"), "w") do stream
        JSON.print(stream, manifest, 2)
        return println(stream)
    end
    verbose && println(
        "Saved ",
        output,
        ": nodes ",
        length(nodes),
        ", tetrahedra ",
        length(tetrahedra),
        ", surface triangles ",
        length(surface_elements),
        "\nSurface counts: ",
        surface_counts,
        "\nSurface areas (um^2): ",
        surface_areas,
        "\nVolumes (um^3): ",
        volumes
    )
    return manifest
end

end # module
