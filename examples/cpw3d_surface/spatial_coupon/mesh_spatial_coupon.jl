# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabrication-resolved local coupon for endpoint, junction, and exact spatial clusters.

import Gmsh: gmsh
using DelimitedFiles
using LinearAlgebra
using SHA
include(joinpath(@__DIR__, "prism_edge_tubes.jl"))

function signature_integer(value, name)
    value isa Real && isfinite(value) && value==round(value) ||
        error("Spatial signature $name must be an integer")
    return Int(value)
end

function read_edges(path)
    data, header = readdlm(path, ',', header=true)
    data = ndims(data) == 1 ? reshape(data, 1, :) : data
    names = vec(String.(header))
    columns = Dict(name => index for (index, name) in enumerate(names))
    required = (
        "Slot",
        "Conductor",
        "Px",
        "Py",
        "Pz",
        "Gx",
        "Gy",
        "Gz",
        "Tx",
        "Ty",
        "Tz",
        "Nz",
        "S0",
        "S1",
        "VertexArm"
    )
    all(haskey(columns, name) for name in required) ||
        error("Spatial signature is missing required columns")
    edges = NamedTuple[]
    for row in axes(data, 1)
        push!(
            edges,
            (
                slot=signature_integer(data[row, columns["Slot"]], "Slot"),
                conductor=signature_integer(data[row, columns["Conductor"]], "Conductor"),
                point=(
                    Float64(data[row, columns["Px"]]),
                    Float64(data[row, columns["Py"]]),
                    Float64(data[row, columns["Pz"]])
                ),
                gap=(
                    Float64(data[row, columns["Gx"]]),
                    Float64(data[row, columns["Gy"]]),
                    Float64(data[row, columns["Gz"]])
                ),
                tangent=(
                    Float64(data[row, columns["Tx"]]),
                    Float64(data[row, columns["Ty"]]),
                    Float64(data[row, columns["Tz"]])
                ),
                normal_sign=Float64(data[row, columns["Nz"]]),
                interval=(
                    Float64(data[row, columns["S0"]]),
                    Float64(data[row, columns["S1"]])
                ),
                vertex_arm=Bool(signature_integer(data[row, columns["VertexArm"]], "VertexArm")),
                # Device-plan coupons (decision 282 rule B3): a context row is a device edge
                # of the plan that is not a claim (a continuation chain piece, chain = true,
                # or a foreign edge); its attributes and owner lookup are those of any edge,
                # the coupon box comes from the process library's SupportBox.
                context=haskey(columns, "Context") ?
                        Bool(signature_integer(data[row, columns["Context"]], "Context")) : false,
                chain=haskey(columns, "Chain") ?
                      Bool(signature_integer(data[row, columns["Chain"]], "Chain")) : false
            )
        )
    end
    isempty(edges) && error("Spatial signature contains no edges")
    all(0 <= edge.slot < 10 && 0 < edge.conductor < 100 for edge in edges) ||
        error("Spatial signature requires slots 0:9 and conductor labels 1:99")
    all(all(isfinite, (edge.point...,edge.gap...,edge.tangent...,edge.interval...,
                      edge.normal_sign)) && edge.interval[1]<=edge.interval[2] &&
        (edge.vertex_arm || edge.interval[1]<edge.interval[2]) for edge in edges) ||
        error("Spatial signature contains non-finite values or invalid intervals")
    all(abs(norm(edge.gap)-1)<1e-6 && abs(norm(edge.tangent)-1)<1e-6 &&
        abs(dot(edge.gap,edge.tangent))<1e-6 for edge in edges) ||
        error("Spatial signature requires orthonormal gap/tangent frames")
    all(
        abs(edge.gap[3]) < 1.0e-8 &&
            abs(edge.tangent[3]) < 1.0e-8 &&
            abs(abs(edge.normal_sign)-1) < 1.0e-8 for edge in edges
    ) || error("Spatial coupon requires parallel or antiparallel process planes")
    edges = [merge(edge,(normal_sign=sign(edge.normal_sign),)) for edge in edges]
    register_edge_chains!(edges)
    return edges
end

# Coupon box rule for CAD-subdivided edges (supervisor decision 47): collinear rows of
# one metal edge that touch end to end (same slot, conductor, plane, tangent and gap,
# no vertex arm) form a chain, and the box extension is decided on the chain's union
# interval instead of each row's own - the subdivision of a straight edge by the
# source CAD must not change the coupon box. The union is expressed in each member
# row's own edge coordinate (EDGE_CHAIN_UNIONS, filled by read_edges), so
# extended_interval keeps its signature and every call site (box, strips, ownership)
# sees one rule. Rows that are not chained (every registered real case: the probe of
# 2026-09-19 found touching chains in the one-edge-cad-subdivided fixture only) keep
# the single-row rule unchanged.
const EDGE_CHAIN_UNIONS = Dict{Any, Tuple{Float64, Float64}}()
const EDGE_CHAIN_RECORDS = Dict{String, Any}[]
const COUPON_BOX_RULE =
    "a row is extended by 2 x Radius at an end reaching Radius from its reference point " *
    "(vertex arms at their far end); collinear touching rows of one metal edge (same slot, " *
    "conductor, plane, tangent, gap; no vertex arm) form a chain judged on the union interval " *
    "about its midpoint, so the source CAD's subdivision of a straight edge never changes the " *
    "coupon box (decision 47; EdgeChains lists every chain with its union length and rows); " *
    "the box is padded by Radius around the extended rows and, per process layer sign Nz, " *
    "by Overetch on the substrate side (-Nz) and MetalThickness on the metal side (+Nz) of " *
    "every row's plane (decision 48)"

edge_chain_key(edge) = (edge.slot, edge.conductor, edge.point, edge.tangent, edge.gap, edge.interval)

function register_edge_chains!(edges)
    empty!(EDGE_CHAIN_UNIONS); empty!(EDGE_CHAIN_RECORDS)
    n = length(edges)
    tolerance = 1.0e-9 * max(1.0, maximum(maximum(abs, edge.point) for edge in edges))
    endpoints(i) = [add(edges[i].point, scale(s, edges[i].tangent)) for s in edges[i].interval]
    same_line(i, j) = edges[i].slot == edges[j].slot && edges[i].conductor == edges[j].conductor &&
        !edges[i].vertex_arm && !edges[j].vertex_arm &&
        all(abs(edges[i].tangent[d] - edges[j].tangent[d]) <= tolerance &&
            abs(edges[i].gap[d] - edges[j].gap[d]) <= tolerance for d in 1:3) &&
        abs(edges[i].point[3] - edges[j].point[3]) <= tolerance &&
        begin
            delta = (edges[j].point[1] - edges[i].point[1], edges[j].point[2] - edges[i].point[2], 0.0)
            along = sum(delta[d] * edges[i].tangent[d] for d in 1:3)
            sqrt(sum((delta[d] - along * edges[i].tangent[d])^2 for d in 1:3)) <= tolerance
        end
    touching(chain, j) = any(sqrt(sum((a[d] - b[d])^2 for d in 1:3)) <= tolerance
                             for i in chain for a in endpoints(i) for b in endpoints(j))
    used = falses(n)
    for i in 1:n
        (used[i] || edges[i].vertex_arm) && continue
        chain = [i]; used[i] = true
        changed = true
        while changed
            changed = false
            for j in 1:n
                (used[j] || !same_line(i, j) || !touching(chain, j)) && continue
                push!(chain, j); used[j] = true; changed = true
            end
        end
        length(chain) > 1 || continue
        origin, tangent = edges[chain[1]].point, edges[chain[1]].tangent
        coordinate(k, s) = sum((edges[k].point[d] - origin[d]) * tangent[d] for d in 1:3) + s
        u0 = minimum(coordinate(k, edges[k].interval[1]) for k in chain)
        u1 = maximum(coordinate(k, edges[k].interval[2]) for k in chain)
        for k in chain
            shift = coordinate(k, 0.0)
            EDGE_CHAIN_UNIONS[edge_chain_key(edges[k])] = (u0 - shift, u1 - shift)
        end
        push!(EDGE_CHAIN_RECORDS, Dict{String, Any}(
            "Rows" => sort(chain), "Slot" => edges[chain[1]].slot, "Conductor" => edges[chain[1]].conductor,
            "UnionLength" => u1 - u0,
            "Union" => [collect(add(origin, scale(u0, tangent))), collect(add(origin, scale(u1, tangent)))]))
    end
    return edges
end

function read_mask(path)
    path === nothing && return NamedTuple[]
    data, header = readdlm(path, ',', header=true)
    data = ndims(data) == 1 ? reshape(data, 1, :) : data
    names = vec(String.(header))
    columns = Dict(name => index for (index, name) in enumerate(names))
    required = ("Facet", "Conductor", "Plane", "X", "Y")
    all(haskey(columns, name) for name in required) ||
        error("Plan-view mask is missing required columns")
    facets = NamedTuple[]
    for facet_index in
        sort!(unique(Int(round(data[row, columns["Facet"]])) for row in axes(data, 1)))
        rows = [
            row for row in axes(data, 1) if
            Int(round(data[row, columns["Facet"]])) == facet_index
        ]
        conductor = Int(round(data[first(rows), columns["Conductor"]]))
        plane = Float64(data[first(rows), columns["Plane"]])
        points = [
            (Float64(data[row, columns["X"]]), Float64(data[row, columns["Y"]])) for
            row in rows
        ]
        all(Int(round(data[row, columns["Conductor"]])) == conductor for row in rows) &&
        all(Float64(data[row, columns["Plane"]]) == plane for row in rows) ||
            error("Inconsistent plan-view mask facet $facet_index")
        conductor > 0 &&
        isfinite(plane) &&
        length(points) >= 3 &&
        all(all(isfinite, point) for point in points) ||
            error("Invalid plan-view mask facet $facet_index")
        push!(facets, (conductor=conductor, plane=plane, points=points))
    end
    isempty(facets) && error("Plan-view mask contains no facets")
    return facets
end

function read_boundary(path)
    path === nothing && return NamedTuple[]
    data, header = readdlm(path, ',', header=true)
    data = ndims(data) == 1 ? reshape(data, 1, :) : data
    names = vec(String.(header))
    columns = Dict(name => index for (index, name) in enumerate(names))
    required = ("Loop", "Vertex", "Conductor", "Plane", "Hole", "Class", "X", "Y")
    all(haskey(columns, name) for name in required) ||
        error("Plan-view boundary is missing required columns")
    loops = NamedTuple[]
    for loop_index in
        sort!(unique(Int(round(data[row, columns["Loop"]])) for row in axes(data, 1)))
        rows = sort!(
            [
                row for row in axes(data, 1) if
                Int(round(data[row, columns["Loop"]])) == loop_index
            ];
            by=row -> Int(round(data[row, columns["Vertex"]]))
        )
        conductor = Int(round(data[first(rows), columns["Conductor"]]))
        plane = Float64(data[first(rows), columns["Plane"]])
        hole = Bool(round(Int, data[first(rows), columns["Hole"]]))
        points = [
            (Float64(data[row, columns["X"]]), Float64(data[row, columns["Y"]])) for
            row in rows
        ]
        classes = [String(data[row, columns["Class"]]) for row in rows]
        all(Int(round(data[row, columns["Conductor"]])) == conductor for row in rows) &&
        all(Float64(data[row, columns["Plane"]]) == plane for row in rows) &&
        all(Bool(round(Int, data[row, columns["Hole"]])) == hole for row in rows) ||
            error("Inconsistent plan-view boundary loop $loop_index")
        conductor > 0 &&
        isfinite(plane) &&
        length(points) >= 3 &&
        all(all(isfinite, point) for point in points) &&
        all(value in ("Physical", "Continuation") for value in classes) ||
            error("Invalid plan-view boundary loop $loop_index")
        arcs, joints = read_boundary_arc_tags(data, columns, rows, loop_index)
        push!(
            loops,
            (conductor=conductor, plane=plane, hole=hole, points=points, classes=classes,
             arcs=arcs, joints=joints)
        )
    end
    isempty(loops) && error("Plan-view boundary contains no loops")
    return loops
end

# The arc columns of the plan-view boundary (block (b) design A1 (4); generate_spatial_response
# ARC_BOUNDARY_COLUMNS): per vertex row the tag of its OUTGOING side - (id, centre, radius,
# sign) of the rebuilt circle the side is a chord of, or nothing on a straight side - and the
# joint record of the vertex - (turn, smooth) at an arc end, nothing elsewhere. A boundary
# without the columns (every legacy coupon) carries nothing: the untagged path is taken bitwise.
const ARC_BOUNDARY_COLUMNS = ("ArcId", "ArcCx", "ArcCy", "ArcR", "ArcSign", "JointTurn", "JointSmooth")

blank_cell(value) = value === "" || (value isa AbstractString && isempty(strip(value)))

function read_boundary_arc_tags(data, columns, rows, loop_index)
    arcs = Vector{Union{Nothing, NamedTuple}}(nothing, length(rows))
    joints = Vector{Union{Nothing, NamedTuple}}(nothing, length(rows))
    any(haskey(columns, name) for name in ARC_BOUNDARY_COLUMNS) || return arcs, joints
    all(haskey(columns, name) for name in ARC_BOUNDARY_COLUMNS) ||
        error("Plan-view boundary carries a partial arc column set")
    cell(row, name) = data[row, columns[name]]
    for (k, row) in enumerate(rows)
        if !blank_cell(cell(row, "ArcId"))
            id = Int(round(Float64(cell(row, "ArcId"))))
            centre = (Float64(cell(row, "ArcCx")), Float64(cell(row, "ArcCy")))
            radius = Float64(cell(row, "ArcR"))
            sign = Int(round(Float64(cell(row, "ArcSign"))))
            id > 0 && all(isfinite, centre) && isfinite(radius) && radius > 0.0 && sign in (-1, 1) ||
                error("Invalid arc tag on plan-view boundary loop $loop_index row $k")
            arcs[k] = (id=id, centre=centre, radius=radius, sign=sign)
        end
        if !blank_cell(cell(row, "JointSmooth"))
            smooth = Int(round(Float64(cell(row, "JointSmooth"))))
            smooth in (0, 1) || error("Invalid joint tag on plan-view boundary loop $loop_index row $k")
            turn = blank_cell(cell(row, "JointTurn")) ? nothing : Float64(cell(row, "JointTurn"))
            joints[k] = (turn=turn, smooth=smooth == 1)
        end
    end
    return arcs, joints
end

# Whether a loop (read_boundary) carries arc tags.
loop_has_arcs(loop) = haskey(loop, :arcs) && any(arc !== nothing for arc in loop.arcs)

# The arc runs of a TAGGED loop (design A1 (4)): one run per ArcId over its consecutive
# chord sides (the loop rotated so that no run straddles its end), in the format of
# circular_arc_runs (centre, radius, point_indices, edge_indices, orientation, angle) plus
# `id`, `sign` and `sweep` (the signed angular travel of the run). The fit is SEEDED by the
# tags: the circle through the run's first, middle and last vertex must agree with the tagged
# centre / radius and every tagged vertex must lie on it within the fit tolerance of
# fitted_arc_run (max(64 tol, 2e-7 rho)); any disagreement fails closed - one source of truth
# (the signature), one circle (the tagged one) for the CAD. A one-chord run is admitted (no
# four-chord minimum).
function tagged_arc_runs(loop, tolerance)
    points = loop.points
    n = length(points)
    runs = NamedTuple[]
    loop_has_arcs(loop) || return runs
    ids = [arc === nothing ? 0 : arc.id for arc in loop.arcs]
    # Start at a straight side (or at a run start) so that no run straddles the end.
    start = findfirst(i -> ids[i] != 0 && ids[mod1(i - 1, n)] != ids[i], 1:n)
    start === nothing && error("a plan-view loop made of one closed arc is not supported by the tube recipe")
    seen = Set{Int}()
    for step in 1:n
        index = mod1(start + step - 1, n)
        id = ids[index]
        id == 0 && continue
        ids[mod1(index - 1, n)] == id && continue    # not the first chord of its run
        id in seen && error("arc $id of the plan-view boundary is not one consecutive run of chords")
        push!(seen, id)
        edge_indices = Int[]
        k = index
        while ids[k] == id && length(edge_indices) < n
            push!(edge_indices, k)
            k = mod1(k + 1, n)
        end
        point_indices = vcat(edge_indices, [mod1(edge_indices[end] + 1, n)])
        tag = loop.arcs[index]
        all(loop.arcs[e] !== nothing && loop.arcs[e].centre == tag.centre && loop.arcs[e].radius == tag.radius &&
            loop.arcs[e].sign == tag.sign for e in edge_indices) ||
            error("arc $id of the plan-view boundary carries two circles")
        circle = (center=tag.centre, radius=tag.radius)
        fit_tolerance = arc_fit_tolerance(tag.radius, tolerance)
        fit = circle_through(points[point_indices[1]], points[point_indices[cld(length(point_indices), 2)]],
                             points[point_indices[end]], tolerance)
        if fit === nothing
            # Two vertices (one chord): the tagged circle must pass through both.
            length(point_indices) == 2 ||
                error("arc $id of the plan-view boundary: its vertices are not concyclic")
        else
            hypot(fit.center[1] - tag.centre[1], fit.center[2] - tag.centre[2]) <= fit_tolerance &&
            abs(fit.radius - tag.radius) <= fit_tolerance ||
                error("arc $id of the plan-view boundary: the fitted circle (centre $(fit.center), radius " *
                      "$(fit.radius)) disagrees with the tagged circle (centre $(tag.centre), radius " *
                      "$(tag.radius)) beyond $fit_tolerance")
        end
        run = fitted_arc_run(points, point_indices, edge_indices, circle, tolerance)
        run === nothing &&
            error("arc $id of the plan-view boundary: its vertices do not lie on the tagged circle within " *
                  "$fit_tolerance, or do not travel monotonically around it")
        push!(runs, (run..., id=id, sign=tag.sign,
                     sweep=run_sweep(points, run.point_indices, run.center, run.orientation)))
    end
    return runs
end

# The signed angular travel of a run from its first to its last vertex about its centre
# (unwrapped along the orientation; a closed run reads a full turn): a function of the two
# END vertices alone, so the equal-fraction split points of polygon_wire and the arc tubes
# do not depend on the chord count (design A7 MINOR-7).
function run_sweep(points, point_indices, center, orientation)
    first = points[point_indices[1]]
    last = points[point_indices[end]]
    theta_first = atan(first[2] - center[2], first[1] - center[1])
    theta_last = atan(last[2] - center[2], last[1] - center[1])
    first == last && return orientation * 2.0 * pi
    travel = mod(orientation * (theta_last - theta_first), 2.0 * pi)
    return orientation * travel
end

# The number of <= 90-degree parts an arc of signed sweep `sweep` is split into (OCC circle
# arcs and revolves stay below a half turn); the parts are equal angular fractions, so no
# part is ever a sliver (a 90.000001-degree arc has two parts of 45.0000005 degrees).
arc_part_count(sweep) = max(1, ceil(Int, abs(sweep) / (0.5 * pi) - 1.0e-9))

# The split angles theta_0 = the first vertex's angle, ..., theta_parts = the last vertex's.
function arc_split_angles(points, point_indices, center, sweep)
    first = points[point_indices[1]]
    theta_first = atan(first[2] - center[2], first[1] - center[1])
    parts = arc_part_count(sweep)
    return [theta_first + sweep * k / parts for k in 0:parts]
end

# The arc runs of a loop: the tagged runs when the boundary carries arc tags, else the
# untagged fit (circular_arc_runs: bitwise for every legacy coupon).
function loop_arc_runs(loop, tolerance)
    loop_has_arcs(loop) && return tagged_arc_runs(loop, tolerance)
    return circular_arc_runs(loop.points, tolerance)
end

add(a, b) = ntuple(i -> a[i] + b[i], 3)
scale(value, vector) = ntuple(i -> value * vector[i], 3)

const IDENTITY_RIGID_TRANSFORM = Float64[
    1 0 0 0;
    0 1 0 0;
    0 0 1 0;
    0 0 0 1
]

function rigid_transform(values::AbstractVector{<:Real})
    length(values) == 16 || error("Rigid transform must contain 16 row-major values")
    matrix = Matrix(reshape(Float64.(values), 4, 4)')
    all(isfinite, matrix) || error("Rigid transform contains a non-finite value")
    isapprox(matrix[4, :], [0.0, 0.0, 0.0, 1.0]; atol=1.0e-12, rtol=0.0) ||
        error("Rigid transform must have homogeneous last row [0, 0, 0, 1]")
    rotation = matrix[1:3, 1:3]
    isapprox(rotation' * rotation, Matrix{Float64}(I, 3, 3); atol=1.0e-12, rtol=0.0) ||
        error("Rigid transform rotation must be orthogonal")
    isapprox(det(rotation), 1.0; atol=1.0e-12, rtol=0.0) ||
        error("Rigid transform must preserve orientation")
    return matrix
end

function parse_rigid_transform(value::String)
    fields = split(value, ',')
    length(fields) == 16 || error("Rigid transform must be 16 comma-separated values")
    return rigid_transform(parse.(Float64, fields))
end

function transform_point(matrix, point)
    value = matrix * [point[1], point[2], point[3], 1.0]
    return (value[1], value[2], value[3])
end

function transform_vector(matrix, vector)
    value = matrix[1:3, 1:3] * collect(vector)
    return (value[1], value[2], value[3])
end

function inverse_transform_point(matrix, point)
    rotation = matrix[1:3, 1:3]
    value = rotation' * (collect(point) - matrix[1:3, 4])
    return (value[1], value[2], value[3])
end

function transform_edge_contract(edge, matrix)
    return merge(edge, (
        point=transform_point(matrix, edge.point),
        gap=transform_vector(matrix, edge.gap),
        tangent=transform_vector(matrix, edge.tangent),
        process_normal=transform_vector(matrix, (0.0, 0.0, edge.normal_sign))
    ))
end

# Minimal strict JSON reader/writer: the test Julia project has no JSON package,
# and the seed must consume the frozen semantic contract exactly as written.
function parse_json(text::AbstractString)
    characters = collect(text)
    position = Ref(1)
    skip_space() = while position[] <= length(characters) && isspace(characters[position[]])
        position[] += 1
    end
    function current()
        position[] <= length(characters) || error("Unexpected end of JSON")
        return characters[position[]]
    end
    function expect(token)
        stop = position[] + length(token) - 1
        stop <= length(characters) && String(characters[position[]:stop]) == token ||
            error("Invalid JSON near position $(position[])")
        position[] = stop + 1
    end
    function parse_string()
        expect("\"")
        buffer = IOBuffer()
        while true
            position[] <= length(characters) || error("Unterminated JSON string")
            character = characters[position[]]
            position[] += 1
            character == '"' && return String(take!(buffer))
            if character == '\\'
                position[] <= length(characters) || error("Unterminated JSON escape")
                escape = characters[position[]]
                position[] += 1
                if escape == 'u'
                    stop = position[] + 3
                    stop <= length(characters) || error("Invalid JSON unicode escape")
                    write(buffer, Char(parse(UInt32, String(characters[position[]:stop]); base=16)))
                    position[] = stop + 1
                else
                    mapped = Dict('"' => '"', '\\' => '\\', '/' => '/', 'b' => '\b',
                                  'f' => '\f', 'n' => '\n', 'r' => '\r', 't' => '\t')
                    haskey(mapped, escape) || error("Invalid JSON escape")
                    write(buffer, mapped[escape])
                end
            else
                write(buffer, character)
            end
        end
    end
    function parse_number()
        start = position[]
        while position[] <= length(characters) && characters[position[]] in "+-0123456789.eE"
            position[] += 1
        end
        token = String(characters[start:(position[] - 1)])
        isempty(token) && error("Invalid JSON near position $start")
        value = tryparse(Int, token)
        value === nothing || return value
        value = tryparse(Float64, token)
        value === nothing && error("Invalid JSON number $token")
        return value
    end
    function parse_value()
        skip_space()
        character = current()
        if character == '{'
            position[] += 1
            result = Dict{String, Any}()
            skip_space()
            if current() == '}'
                position[] += 1
                return result
            end
            while true
                skip_space()
                key = parse_string()
                skip_space()
                expect(":")
                result[key] = parse_value()
                skip_space()
                current() == ',' && (position[] += 1; continue)
                expect("}")
                return result
            end
        elseif character == '['
            position[] += 1
            result = Any[]
            skip_space()
            if current() == ']'
                position[] += 1
                return result
            end
            while true
                push!(result, parse_value())
                skip_space()
                current() == ',' && (position[] += 1; continue)
                expect("]")
                return result
            end
        elseif character == '"'
            return parse_string()
        elseif character == 't'
            expect("true"); return true
        elseif character == 'f'
            expect("false"); return false
        elseif character == 'n'
            expect("null"); return nothing
        else
            return parse_number()
        end
    end
    value = parse_value()
    skip_space()
    position[] > length(characters) || error("Trailing characters after JSON value")
    return value
end

function write_json(io::IO, value, indent::Int=0)
    pad = repeat(" ", indent)
    if value isa AbstractDict
        keys_sorted = sort!(collect(keys(value)))
        println(io, "{")
        for (index, key) in enumerate(keys_sorted)
            print(io, pad, "  \"", key, "\": ")
            write_json(io, value[key], indent + 2)
            println(io, index < length(keys_sorted) ? "," : "")
        end
        print(io, pad, "}")
    elseif value isa AbstractVector || value isa Tuple
        if isempty(value)
            print(io, "[]")
        else
            println(io, "[")
            for (index, item) in enumerate(value)
                print(io, pad, "  ")
                write_json(io, item, indent + 2)
                println(io, index < length(value) ? "," : "")
            end
            print(io, pad, "]")
        end
    elseif value isa AbstractString
        print(io, "\"", escape_string(value), "\"")
    elseif value isa Bool
        print(io, value ? "true" : "false")
    elseif value isa Integer
        print(io, value)
    elseif value isa Real
        isfinite(value) || error("JSON cannot record a non-finite number")
        print(io, Float64(value))
    elseif value === nothing
        print(io, "null")
    else
        error("Unsupported JSON value of type $(typeof(value))")
    end
end

function json_point(value)
    value isa AbstractVector && length(value) == 3 && all(x -> x isa Real && isfinite(x), value) ||
        error("Semantic contract corners must be finite 3D points")
    return (Float64(value[1]), Float64(value[2]), Float64(value[3]))
end

# Contract semantic corners pulled back into the seed's source-local frame.
function read_semantic_corners(path, transform)
    contract = parse_json(read(path, String))
    contract isa AbstractDict && haskey(contract, "SemanticCorners") ||
        error("Semantic contract lacks SemanticCorners")
    corners = contract["SemanticCorners"]
    corners isa AbstractVector && !isempty(corners) ||
        error("Semantic contract must list at least one semantic corner")
    placement = haskey(contract, "RigidTransform") ?
        rigid_transform(Float64.(contract["RigidTransform"])) : copy(IDENTITY_RIGID_TRANSFORM)
    isapprox(placement, transform; atol=1.0e-12, rtol=0.0) ||
        error("Semantic contract placement differs from the seed rigid transform")
    return [inverse_transform_point(transform, json_point(corner)) for corner in corners]
end

# Model point tags coincident with every semantic corner; fail closed otherwise.
function semantic_corner_points(corners, tolerance)
    entities = gmsh.model.getEntities(0)
    coordinates = Dict(tag => gmsh.model.getValue(dim, tag, Float64[]) for (dim, tag) in entities)
    tags = Int32[]
    for corner in corners
        matches = [tag for (tag, xyz) in coordinates if norm(xyz .- collect(corner)) <= tolerance]
        isempty(matches) && error("Semantic corner $(corner) is absent from the seed CAD")
        append!(tags, matches)
    end
    return sort!(unique(tags))
end

# Isotropic size law inside a corner ball (supervisor decision 33): with a
# CornerSize the size grows geometrically from CornerSize at the corner point by
# the recorded ratio to lc_fine (NormalSize) in shells - shell k (k = 1, 2, ...)
# has size CornerSize ratio^(k-1) and ends at the cumulative radius
# CornerSize (ratio^k - 1) / (ratio - 1), exactly the edge layer's rows with
# CornerSize for EdgeSize (4 nm: sizes 4/8/16 nm to 4/12/28 nm, lc_fine beyond) -
# and is lc_fine from the last shell (the grading reach (lc_fine - CornerSize) /
# (ratio - 1)) to the ball radius; without a CornerSize (0, production) the ball
# is uniformly lc_fine. The seed uses the shell (staircase) form, so the ridge
# nodes through a corner fall on the shell radii and the corner cells have the
# shell sizes; the metric stage prescribes the continuous form CornerSize +
# (ratio - 1) d of the same law (edge_volume_metric.corner_ball_size), which the
# shells never exceed.
struct CornerGrading
    corner_size::Float64
    ratio::Float64
    lc_fine::Float64
    radius::Float64
end

corner_grading_reach(grading::CornerGrading) =
    grading.corner_size > 0.0 ? (grading.lc_fine - grading.corner_size) / (grading.ratio - 1.0) : 0.0

# Radii of the shell boundaries inside the ball: the cumulative geometric
# offsets below lc_fine (edge_layer_row_offsets with CornerSize), then the radius.
corner_shell_radii(grading::CornerGrading) =
    grading.corner_size > 0.0 ?
        vcat(edge_layer_row_offsets(grading.corner_size, grading.ratio, grading.lc_fine),
             grading.radius) : [grading.radius]

function corner_ball_size(grading::CornerGrading, distance)
    grading.corner_size > 0.0 || return grading.lc_fine
    size = grading.corner_size
    for shell_radius in edge_layer_row_offsets(grading.corner_size, grading.ratio, grading.lc_fine)
        distance >= shell_radius || break
        size *= grading.ratio
    end
    return min(grading.lc_fine, size)
end

function corner_size_expression(distance, grading::CornerGrading, lc_far, transition_width)
    # The corner law inside the ball (lc_fine uniformly, or the shells of the
    # geometric grading from CornerSize: Gmsh MathEval step(x) is 1 for x >= 0),
    # then the same linear grading slope as the process edge band up to the far
    # size.
    lc_fine, radius = grading.lc_fine, grading.radius
    inside = "$(lc_fine)"
    if grading.corner_size > 0.0
        inside = "$(grading.corner_size)"
        size = grading.corner_size
        for shell_radius in edge_layer_row_offsets(grading.corner_size, grading.ratio, lc_fine)
            next_size = min(lc_fine, size * grading.ratio)
            inside *= "+$(next_size - size)*step($(distance)-$(shell_radius))"
            size = next_size
        end
        inside = "min($(lc_fine),$(inside))"
    end
    return "min($(lc_far),$(inside)+($(lc_far)-$(lc_fine))*" *
           "max($(distance)-$(radius),0)/$(transition_width))"
end

# Tangential size of a longitudinal curve through the corner balls: the corner
# law inside a ball, the process-band grading slope up to lc_tangent outside.
function corner_curve_size(point, corners, grading::CornerGrading, lc_tangent, slope)
    distance = minimum(norm(point .- corner) for corner in corners)
    return min(lc_tangent, corner_ball_size(grading, distance) +
                           slope * max(distance - grading.radius, 0.0))
end

# Samples per spacing of the chord polyline that tabulates the arclength of a
# longitudinal curve span; a resolution constant, not a mesh target.
const CURVE_SIZE_SAMPLES_PER_FINE_LENGTH = 16

# Interior node parameters of a longitudinal curve under the corner law alone
# (the pre-decision-43 rule of the spike and the face-census tests): the
# composed rule without the trace and band laws. Returns nothing when the curve
# is out of reach of every ball.
function corner_isotropic_curve_nodes(curve, corners, grading::CornerGrading, lc_tangent, slope)
    placed = composed_curve_nodes(
        curve, lc_tangent, point -> corner_curve_size(point, corners, grading, lc_tangent, slope),
        grading.ratio, corners, grading)
    return placed === nothing ? nothing : (placed[1], placed[2])
end

# Parameters in (from, to) where the curve crosses a corner ball boundary
# (distance to the nearest corner == radius), located by sampling and bisection.
# A crossing closer than half the ball's boundary size (lc_fine) to a gap end is
# not a node: a CAD vertex on or next to the ball boundary (the top of a
# metal-thickness vertical corner edge when the thickness equals the radius) would
# otherwise get a node at roundoff distance and degenerate cells.
function ball_boundary_parameters(curve, from, to, corners, grading::CornerGrading)
    radius = grading.radius
    samples = 256
    parameters = collect(range(from, to; length=samples + 1))
    xyz = reshape(gmsh.model.getValue(1, curve, parameters), 3, :)
    excess(i) = minimum(norm(xyz[:, i] .- collect(corner)) for corner in corners) - radius
    excess_at(t) = minimum(norm(gmsh.model.getValue(1, curve, [t]) .- collect(corner))
                           for corner in corners) - radius
    crossings = Float64[]
    for i in 1:samples
        left, right = excess(i), excess(i + 1)
        (left == 0.0 || sign(left) == sign(right)) && continue
        lo, hi = parameters[i], parameters[i + 1]
        for _ in 1:60
            mid = 0.5 * (lo + hi)
            if sign(excess_at(mid)) == sign(left)
                lo = mid
            else
                hi = mid
            end
        end
        crossing = 0.5 * (lo + hi)
        from < crossing < to || continue
        point = gmsh.model.getValue(1, curve, [crossing])
        separation = min(norm(point .- xyz[:, 1]), norm(point .- xyz[:, end]))
        separation >= 0.5 * grading.lc_fine && push!(crossings, crossing)
    end
    return crossings
end

# Samples of the composed law per explicit-spacing grid interval when deciding
# whether the interval keeps the grid; a resolution constant, not a mesh target.
const CURVE_LAW_SAMPLES_PER_INTERVAL = 8

const CURVE_SPACING_RULE =
    "every explicitly 1D-meshed longitudinal curve (band curves at NormalSize, un-tubed " *
    "metal ridge parts at TangentialSize) follows min(its explicit Spacing, the composed " *
    "size field on the curve: corner-ball law of the graded points, trace-basis rule, band " *
    "rule, corner exterior rule), sampled along the curve and gradient-limited to " *
    "(GrowthRatio - 1) / GrowthRatio as on the tube axes (graded_tube_stations, the " *
    "decision-40 rule); the Spacing grid is kept on every grid interval where the limited " *
    "law equals the Spacing (ridge alignment across faces), each remaining gap between kept " *
    "grid nodes (split at the corner ball boundaries) is equidistributed in the arclength " *
    "integral of the reciprocal limited law with ceil(integral) intervals, so every node " *
    "interval is <= the size it spans; AchievedOverPrescribed = each node interval over the " *
    "limited law at its midpoint - decision 43"

# Chord-polyline arclength table of the curve span [from, to] at
# CURVE_SIZE_SAMPLES_PER_FINE_LENGTH samples per spacing (exact for the straight
# feature curves): parameter_at(s) and arclength_at(t) by linear interpolation.
function curve_arclength_table(curve, from, to, spacing)
    to > from || error("Longitudinal curve span is not increasing")
    chord = norm(gmsh.model.getValue(1, curve, [to]) .- gmsh.model.getValue(1, curve, [from]))
    samples = max(16, ceil(Int, CURVE_SIZE_SAMPLES_PER_FINE_LENGTH * chord / spacing))
    parameters = collect(range(from, to; length=samples + 1))
    xyz = reshape(gmsh.model.getValue(1, curve, parameters), 3, :)
    arclength = zeros(samples + 1)
    for i in 1:samples
        arclength[i + 1] = arclength[i] + norm(xyz[:, i + 1] .- xyz[:, i])
    end
    arclength[end] > 0.0 || error("Longitudinal curve span has no length")
    function parameter_at(s)
        i = clamp(searchsortedlast(arclength, s), 1, samples)
        fraction = (s - arclength[i]) / (arclength[i + 1] - arclength[i])
        return parameters[i] + fraction * (parameters[i + 1] - parameters[i])
    end
    function arclength_at(t)
        i = clamp(searchsortedlast(parameters, t), 1, samples)
        fraction = (t - parameters[i]) / (parameters[i + 1] - parameters[i])
        return arclength[i] + fraction * (arclength[i + 1] - arclength[i])
    end
    return (; parameters, arclength, parameter_at, arclength_at)
end

# Linear interpolation of the sampled (positions, sizes) law.
function sampled_size_at(positions, sizes, s)
    i = clamp(searchsortedlast(positions, s), 1, length(positions) - 1)
    fraction = (s - positions[i]) / (positions[i + 1] - positions[i])
    return sizes[i] + fraction * (sizes[i + 1] - sizes[i])
end

# Interior node parameters of an explicitly 1D-meshed longitudinal curve under the
# composed size law (CURVE_SPACING_RULE): `size_at(point)` is the law already
# capped at `spacing`, `growth` the layer growth ratio, `corners`/`grading` the
# graded points and corner grading for the ball boundary nodes (nothing without
# corner isotropy). Returns nothing when the limited law is the spacing on every
# grid interval (the caller keeps the transfinite grid), otherwise (parameters,
# coordinates, record) with the per-curve spacing statistics.
function composed_curve_nodes(curve, spacing, size_at, growth, corners, grading)
    lower, upper = gmsh.model.getParametrizationBounds(1, curve)
    curve_length = gmsh.model.occ.getMass(1, curve)
    intervals = max(1, ceil(Int, curve_length / spacing))
    grid = collect(range(lower[1], upper[1]; length=intervals + 1))
    table = curve_arclength_table(curve, lower[1], upper[1], spacing)
    law_at(s) = size_at(Tuple(gmsh.model.getValue(1, curve, [table.parameter_at(s)])))
    # The whole-curve sampled and gradient-limited law (graded_tube_stations
    # returns its adaptive samples), so the neighbour ratio bound holds across
    # kept grid intervals and graded gaps alike.
    _, positions, sizes = graded_tube_stations(0.0, table.arclength[end], law_at, spacing, growth)
    limited(s) = sampled_size_at(positions, sizes, s)
    threshold = spacing * (1.0 - 1.0e-9)
    grid_s = [table.arclength_at(t) for t in grid]
    # A grid interval keeps the grid when the limited law is the spacing at its ends
    # and every sample inside; a grid node stays a node when the law is the spacing
    # there (a dip strictly inside an interval is graded between its two grid nodes).
    at_spacing = [limited(grid_s[i]) >= threshold && limited(grid_s[i + 1]) >= threshold &&
                  all(sizes[j] >= threshold for j in searchsortedfirst(positions, grid_s[i]):
                                                     searchsortedlast(positions, grid_s[i + 1]))
                  for i in 1:intervals]
    all(at_spacing) && return nothing
    anchors = [i for i in 1:(intervals + 1)
               if i == 1 || i == intervals + 1 || limited(grid_s[i]) >= threshold]
    interior = Float64[]
    for (a, b) in zip(anchors[1:(end - 1)], anchors[2:end])
        if !(b == a + 1 && at_spacing[a])
            bounds = corners !== nothing && grading.corner_size > 0.0 ?
                ball_boundary_parameters(curve, grid[a], grid[b], corners, grading) : Float64[]
            stops = vcat(grid[a], bounds, grid[b])
            for (from, to) in zip(stops[1:(end - 1)], stops[2:end])
                stations, _, _ = graded_tube_stations(table.arclength_at(from), table.arclength_at(to),
                                                      limited, spacing, growth)
                append!(interior, table.parameter_at(s) for s in stations[2:(end - 1)])
                to < grid[b] && push!(interior, to)
            end
        end
        b <= intervals && push!(interior, grid[b])
    end
    coordinates = gmsh.model.getValue(1, curve, interior)
    node_s = vcat(0.0, [table.arclength_at(t) for t in interior], table.arclength[end])
    achieved = diff(node_s)
    all(achieved .> 0.0) || error("Longitudinal curve $curve received a coincident node")
    prescribed = [limited(0.5 * (node_s[i] + node_s[i + 1])) for i in eachindex(achieved)]
    ratios = achieved ./ prescribed
    record = Dict{String, Any}(
        "Length" => curve_length, "Spacing" => spacing, "Graded" => true,
        "GridIntervals" => intervals, "GridIntervalsKept" => count(at_spacing),
        "InteriorNodes" => length(interior),
        "NodeSpacing" => Dict{String, Any}("Minimum" => minimum(achieved),
                                           "P50" => sorted_median(achieved),
                                           "Maximum" => maximum(achieved)),
        "PrescribedMinimum" => minimum(sizes),
        "AchievedOverPrescribed" => Dict{String, Any}("Minimum" => minimum(ratios),
                                                     "P50" => sorted_median(ratios),
                                                     "Maximum" => maximum(ratios)))
    return interior, coordinates, record
end

# Install an explicit curve mesh (endpoints, interior nodes, line elements) so
# that Mesh.MeshOnlyEmpty keeps it when the remaining entities are generated.
function add_explicit_curve_mesh!(curve, parameters, coordinates, next_node, point_nodes)
    lower, upper = gmsh.model.getParametrizationBounds(1, curve)
    ends = gmsh.model.getValue(1, curve, [lower[1], upper[1]])
    _, points = gmsh.model.getAdjacencies(1, curve)
    length(points) == 2 || error("Longitudinal curve $curve is not bounded by two points")
    function endpoint_node(target)
        tag = points[argmin([norm(gmsh.model.getValue(0, point, Float64[]) .- target)
                             for point in points])]
        return get!(point_nodes, tag) do
            next_node[] += 1
            gmsh.model.mesh.addNodes(0, tag, [next_node[]], gmsh.model.getValue(0, tag, Float64[]))
            next_node[]
        end
    end
    first_node = endpoint_node(ends[1:3])
    last_node = endpoint_node(ends[4:6])
    interior = collect((next_node[] + 1):(next_node[] + length(parameters)))
    next_node[] += length(parameters)
    isempty(interior) || gmsh.model.mesh.addNodes(1, curve, interior, coordinates, parameters)
    sequence = [first_node; interior; last_node]
    connectivity = Int[]
    for i in 1:(length(sequence) - 1)
        push!(connectivity, sequence[i], sequence[i + 1])
    end
    gmsh.model.mesh.addElementsByType(curve, 1, Int[], connectivity)
    return length(interior)
end

# Metal surface families (thin-film metal, metal-substrate, metal-vacuum): the
# feature curves bounding one of them are metal edges.
const METAL_SURFACE_FAMILIES = (4, 5, 6)

# Geometric transverse edge layer: the layer k (k = 1, 2, ...) has size
# edge_size * ratio^(k - 1) and every layer strictly finer than the normal size
# lc_fine is seeded, so the rows lie at the cumulative distances
# edge_size * (ratio^k - 1) / (ratio - 1) from the edge.
function edge_layer_row_offsets(edge_size, ratio, lc_fine)
    offsets = Float64[]
    size = edge_size
    while size < lc_fine
        push!(offsets, isempty(offsets) ? size : offsets[end] + size)
        size *= ratio
    end
    return offsets
end

# Unit normal of a planar face from its boundary: the ridge tangent and the
# first boundary point off the ridge line.
function planar_face_normal(face, origin, tangent, tolerance)
    for (curve_dim, curve) in gmsh.model.getBoundary([(2, face)], false, false, false)
        curve_dim == 1 || continue
        lower, upper = gmsh.model.getParametrizationBounds(1, curve)
        for parameter in range(lower[1], upper[1]; length=5)
            point = gmsh.model.getValue(1, curve, [parameter])
            offset = point .- origin
            lateral = offset .- dot(offset, tangent) .* tangent
            norm(lateral) > tolerance || continue
            normal = cross(tangent, lateral)
            return normal ./ norm(normal)
        end
    end
    return error("Face $face has no boundary point off the ridge line")
end

# Smallest power of two n with spacing / n <= target (at least 1).
function tangential_subdivision(spacing, target)
    target > 0.0 || error("tangential subdivision target must be positive")
    n = 1
    while spacing / n > target
        n *= 2
    end
    return n
end

# Explicit node rows parallel to a metal ridge on one adjacent planar face: one
# row per geometric layer, each row's nodes on the ridge span grid subdivided by
# the row's nested power of two (aligned rows keep the strips 2D-mesher
# independent and their edges Delaunay). Every second node of a row is set
# EDGE_LAYER_ROW_ZIGZAG of the row offset farther from the ridge: perfectly
# aligned parallel rows form rectangles whose cocircular corners are
# Delaunay-degenerate, and the boundary recovery was measured to leave
# zero-volume tets in the face plane there; the zigzag makes every quad a
# non-cyclic right trapezoid with a unique triangulation. A row is placed only
# while the next layer still fits inside the face, so rows of different ridges
# or a face narrower than the layer never collide. Returns the rows' record or
# nothing when no row fits.
const EDGE_LAYER_ROW_ZIGZAG = 0.05

function edge_layer_row_nodes(span, offset, inward, subdivision)
    columns = Vector{Float64}[]
    for i in 1:(size(span, 2) - 1), j in 0:(subdivision - 1)
        i == 1 && j == 0 && continue
        push!(columns, span[:, i] .+ (j / subdivision) .* (span[:, i + 1] .- span[:, i]))
    end
    interior = isempty(columns) ? zeros(3, 0) : reduce(hcat, columns)
    scale = [offset * (1.0 + EDGE_LAYER_ROW_ZIGZAG * (i % 2)) for i in axes(interior, 2)]
    return interior .+ inward * scale'
end

function add_edge_layer_rows!(occ, face, ridge_xyz, offsets, tolerance)
    size(ridge_xyz, 2) >= 2 || return nothing
    origin = ridge_xyz[:, 1]
    tangent = ridge_xyz[:, end] .- origin
    tangent ./= norm(tangent)
    normal = planar_face_normal(face, origin, tangent, tolerance)
    inward = cross(normal, tangent)
    middle = ridge_xyz[:, (size(ridge_xyz, 2) + 1) ÷ 2]
    probe = 2 * offsets[1]
    if gmsh.model.isInside(2, face, middle .+ probe .* inward, false) == 0
        gmsh.model.isInside(2, face, middle .- probe .* inward, false) == 0 && return nothing
        inward = -inward
    end
    rows = Int32[]
    for offset in offsets
        fits = all(gmsh.model.isInside(2, face, ridge_xyz[:, i] .+ (2offset) .* inward,
                                       false) == 1
                   for i in (1, (size(ridge_xyz, 2) + 1) ÷ 2, size(ridge_xyz, 2)))
        fits || break
        first_point = occ.addPoint((ridge_xyz[:, 1] .+ offset .* inward)...)
        last_point = occ.addPoint((ridge_xyz[:, end] .+ offset .* inward)...)
        push!(rows, occ.addLine(first_point, last_point))
    end
    isempty(rows) && return nothing
    return Dict{String, Any}("Face" => Int(face), "Rows" => rows, "Inward" => collect(inward))
end

# Mesh the embedded rows explicitly with the subdivided span grid positions
# offset into the face; the rows were created in the same order as the offsets.
function mesh_edge_layer_rows!(record, span, offsets, subdivisions, next_node, point_nodes)
    inward = record["Inward"]
    for (row, offset, subdivision) in zip(record["Rows"], offsets, subdivisions)
        lower, upper = gmsh.model.getParametrizationBounds(1, row)
        start = gmsh.model.getValue(1, row, [lower[1]])
        first_point = span[:, 1] .+ offset .* inward
        last_point = span[:, end] .+ offset .* inward
        length_along = norm(last_point .- first_point)
        forward = norm(start .- first_point) <= norm(start .- last_point)
        interior = edge_layer_row_nodes(span, offset, inward, subdivision)
        fractions = [dot(interior[:, i] .- first_point, last_point .- first_point) /
                     length_along^2 for i in axes(interior, 2)]
        parameters = [forward ? lower[1] + f * (upper[1] - lower[1]) :
                                upper[1] - f * (upper[1] - lower[1]) for f in fractions]
        add_explicit_curve_mesh!(row, parameters, vec(interior), next_node, point_nodes)
    end
    return length(record["Rows"])
end

function tetrahedron_aspect(xyz)
    jacobian = hcat(xyz[2] .- xyz[1], xyz[3] .- xyz[1], xyz[4] .- xyz[1])
    singular = svdvals(jacobian)
    return singular[1] / singular[end]
end

# ---------------------------------------------------------------------------
# Order-invariant corner cell measures (mesher design round 2 F5-A, supervisor
# decisions 351 / 358 / 363). tetrahedron_aspect is MFEM's linear-tet Jacobian
# condition against the RIGHT-ANGLED reference tetrahedron in the STORED node order:
# the same cell reads differently by first vertex (O1 thin 22.5-degree tip: 5.03 /
# 4.63 / 4.61 / 9.10), a corner-first cell spanning a tip of opening phi is floored
# at cot(phi / 2), and the solver reorients the cells at load, so the measure is a
# node-order lottery at non-perpendicular corners. kappa_reg is the condition number
# of the affine map from the REGULAR tetrahedron (J W^-1, W the regular tetrahedron's
# edge matrix): the regular tetrahedron's relabellings are orthogonal maps, so the
# value is the same in every node order; kappa_reg = 1 for the regular tetrahedron,
# 2 for the trirectangular one, and lies in [kappa_v0 / 2, 2 kappa_v0] for every
# frame (kappa(W) = 2). It is the objective and the verdict at INVARIANT corners
# (semantic_corner_kinds); LEGACY (exactly perpendicular) corners keep
# tetrahedron_aspect bitwise. The frame extremes and the mean ratio are recorded
# as information at every corner (FIX-UPS F.2 / F.6).
const REGULAR_TETRAHEDRON_EDGES_INVERSE =
    inv(hcat([1.0, 0.0, 0.0], [0.5, sqrt(3.0) / 2.0, 0.0], [0.5, sqrt(3.0) / 6.0, sqrt(2.0 / 3.0)]))

function tetrahedron_regular_condition(xyz)
    jacobian = hcat(xyz[2] .- xyz[1], xyz[3] .- xyz[1], xyz[4] .- xyz[1]) *
               REGULAR_TETRAHEDRON_EDGES_INVERSE
    singular = svdvals(jacobian)
    return singular[1] / singular[end]
end

# (minimum, maximum) of tetrahedron_aspect over the four first-vertex choices.
function tetrahedron_frame_conditions(xyz)
    values = [tetrahedron_aspect(xyz[[i; setdiff(1:4, i)]]) for i in 1:4]
    return minimum(values), maximum(values)
end

# Mean ratio 12 (3 V)^(2/3) / sum of the squared edge lengths: 1 for the regular
# tetrahedron, 0 for a flat one (information).
function tetrahedron_mean_ratio(xyz)
    a = xyz[2] .- xyz[1]; b = xyz[3] .- xyz[1]; c = xyz[4] .- xyz[1]
    volume = abs(dot(a, cross(b, c))) / 6.0
    squared = sum(sum((xyz[i] .- xyz[j]) .^ 2) for i in 1:4 for j in (i + 1):4)
    return 12.0 * (3.0 * volume)^(2.0 / 3.0) / squared
end

# The optimizer goal at an invariant corner: the constant 3.8 (= 0.95 x 4.0, the
# number every legacy corner descends to), decoupled from the verdict CornerShapeGate
# (FIX-UPS F.2: the after-values are then the E3 production predictions by
# construction; a goal tied to the gate would stop the descent at shapes the
# qualified population never produced).
const INVARIANT_CORNER_TARGET = 3.8

const INVARIANT_CORNER_RULE =
    "a semantic corner is LEGACY when the two plan-view boundary sides meeting at it have " *
    "an exactly zero dot product in exact integer arithmetic on their quantum counts (every " *
    "coordinate divided by the plan-view quantum 1e-9 R and rounded to the generator's " *
    "integer; every rectilinear corner, the theta-0 box vertex of decision 320 - its metal " *
    "side exactly perpendicular to the box side - and every rigidly rotated perpendicular " *
    "corner, whatever its rotation): it keeps the vertex-0 Jacobian condition " *
    "tetrahedron_aspect, MaximumCornerAspect and the 0.95 target bitwise; every other corner " *
    "is INVARIANT: its corner-incident seed cells are optimized on and judged by kappa_reg, " *
    "the condition number of the affine map from the regular tetrahedron (order-invariant), " *
    "descended to the fixed goal 3.8 and judged against CornerShapeGate = min(E_pop, 5.0) " *
    "(E_pop the kappa_reg envelope of the (F)-qualified 90-degree corners); a BridgingSliver " *
    "candidate (a corner-incident cell whose four vertices all lie on the kink's two sidewalls, " *
    "at least one strictly on each, in a wedge of obtuse opening) still above the gate after " *
    "the descent triggers the corner-local reconnection pass (supervisor decision 365: the " *
    "measure is the verdict, the predicate the trigger); the contract " *
    "(derive_semantic_contract: Derivation.InvariantCorners) and the mesher evaluate the " *
    "same predicate on the same quantum and a disagreement fails closed (supervisor " *
    "decisions 351 / 358 / 363 / 416)"

# The plan-view boundary's coordinate quantum per unit coupon radius: the generator
# (generate_spatial_response.plan_view_boundary_loops) snaps every loop vertex to the
# 1e-9 R grid and writes the float k x 1e-9 R (semantic_mesh_contract.
# PLAN_VIEW_QUANTUM_OVER_RADIUS is the same number; both tools count in 1e-9 x --radius).
const PLAN_VIEW_QUANTUM_OVER_RADIUS = 1.0e-9

# The integer count of `quantum` nearest to a plan-view coordinate (the generator's
# rounding: half away from zero), exact for every generated coordinate - the float
# k x quantum lies within 1e-6 quanta of k. Int128: two sides of ~1e10 quanta each
# multiply beyond Int64.
function plan_view_quantum_count(value, quantum)
    scaled = value / quantum
    return scaled >= 0.0 ? floor(Int128, scaled + 0.5) : ceil(Int128, scaled - 0.5)
end

# The two plan-view boundary sides meeting at every semantic corner, as unit 2D
# directions pointing away from the corner, the unit direction into the metal between
# them (the kink's metal sector: inside an exterior loop, outside a hole), and the
# corner's kind (:legacy / :invariant) by the exact dot-product predicate above, evaluated
# in integer arithmetic on the quantum counts of the three vertices (`quantum` = 1e-9 R,
# the generator's grid): exactly 0 at every perpendicular corner, axis-aligned or rigidly
# rotated. (The float dot product of the float side vectors is exact only for axis-aligned
# sides: the loop end's rotated perpendicular corner read -1.776e-15 and was classed
# invariant by both tools; supervisor decision 416.) The returned dots are those integers.
# Every corner must be exactly one boundary vertex of its plane (fail closed otherwise).
# `corners` and `loops` share the seed frame.
function semantic_corner_kinds(corners, loops, tolerance, quantum)
    kinds = Symbol[]; sides = NamedTuple[]; dots = Int128[]
    for corner in corners
        matches = NamedTuple[]
        for loop in loops
            abs(loop.plane - corner[3]) <= tolerance || continue
            n = length(loop.points)
            for i in 1:n
                p = loop.points[i]
                hypot(p[1] - corner[1], p[2] - corner[2]) <= tolerance || continue
                push!(matches, (before=loop.points[mod1(i - 1, n)],
                                after=loop.points[mod1(i + 1, n)], point=p, loop=loop))
            end
        end
        length(matches) == 1 ||
            error("Semantic corner $(corner) matches $(length(matches)) plan-view boundary " *
                  "vertices of its plane (exactly one is required)")
        before, after, p, loop = matches[1].before, matches[1].after, matches[1].point, matches[1].loop
        a = [before[1] - p[1], before[2] - p[2]]
        b = [after[1] - p[1], after[2] - p[2]]
        all(norm(v) > 0.0 for v in (a, b)) ||
            error("Degenerate plan-view boundary side at the semantic corner $(corner)")
        quantised(point) =
            (plan_view_quantum_count(point[1], quantum), plan_view_quantum_count(point[2], quantum))
        qp = quantised(p)
        qa = quantised(before) .- qp
        qb = quantised(after) .- qp
        product = qa[1] * qb[1] + qa[2] * qb[2]
        push!(kinds, product == 0 ? :legacy : :invariant)
        wall_1, wall_2 = a ./ norm(a), b ./ norm(b)
        # The metal direction: the bisector (or its perpendicular at a straight joint),
        # oriented by a probe just inside the loop (inside an exterior loop = metal).
        bisector = wall_1 .+ wall_2
        metal = norm(bisector) > 1.0e-9 ? bisector ./ norm(bisector) : [-wall_1[2], wall_1[1]]
        probe_distance = 1.0e-3 * min(norm(a), norm(b))
        probe = (p[1] + probe_distance * metal[1], p[2] + probe_distance * metal[2])
        point_in_polygon(probe, loop.points, tolerance) == !loop.hole || (metal = -metal)
        probe = (p[1] + probe_distance * metal[1], p[2] + probe_distance * metal[2])
        point_in_polygon(probe, loop.points, tolerance) == !loop.hole ||
            error("Unable to orient the metal side of the semantic corner $(corner)")
        push!(sides, (walls=[wall_1, wall_2], metal=metal))
        push!(dots, product)
    end
    return kinds, sides, dots
end

# Whether the plan-view direction `u` (from the corner) lies in the metal sector of a
# corner: the sector from wall 1 to wall 2 that contains the metal direction.
function in_metal_sector(sides, u)
    angle(v) = mod(atan(v[2], v[1]) - atan(sides.walls[1][2], sides.walls[1][1]), 2.0 * pi)
    span = angle(sides.walls[2])
    return (angle(u) <= span) == (angle(sides.metal) <= span)
end

# The opening angle (radians) of the wedge between the two sidewalls that contains the
# plan-view direction `u`: the metal sector's angle or its complement to 2 pi.
function wedge_opening(sides, u)
    angle(v) = mod(atan(v[2], v[1]) - atan(sides.walls[1][2], sides.walls[1][1]), 2.0 * pi)
    span = angle(sides.walls[2])
    metal_span = angle(sides.metal) <= span ? span : 2.0 * pi - span
    return in_metal_sector(sides, u) ? metal_span : 2.0 * pi - metal_span
end

# The contract's recorded invariant corners (Derivation.InvariantCorners.Points,
# pulled back into the seed frame; an empty list when the record is absent).
function read_invariant_corners(path, transform)
    contract = parse_json(read(path, String))
    derivation = get(contract, "Derivation", nothing)
    derivation isa AbstractDict && haskey(derivation, "InvariantCorners") || return NTuple{3, Float64}[]
    record = derivation["InvariantCorners"]
    record isa AbstractDict && record["Points"] isa AbstractVector ||
        error("Semantic contract Derivation.InvariantCorners must record Points")
    return [inverse_transform_point(transform, json_point(point)) for point in record["Points"]]
end

# The mesher's corner kinds against the contract's record: the set of invariant
# corners must agree exactly (a contract derived before the rule, or by a tool
# disagreeing with the mesher, fails closed: regenerate it).
function check_invariant_corner_contract(corners, kinds, recorded, tolerance)
    found = [collect(corner) for (corner, kind) in zip(corners, kinds) if kind === :invariant]
    matched(point, list) = any(norm(collect(point) .- other) <= tolerance for other in list)
    missing_in_contract = [point for point in found if !matched(point, [collect(r) for r in recorded])]
    missing_in_mesher = [point for point in recorded if !matched(point, found)]
    isempty(missing_in_contract) && isempty(missing_in_mesher) ||
        error("the semantic contract's invariant corners disagree with the plan-view boundary: " *
              "invariant by the mesher but not recorded $(missing_in_contract), recorded but " *
              "perpendicular by the mesher $(missing_in_mesher): regenerate the contract " *
              "(derive_semantic_contract records Derivation.InvariantCorners; $(INVARIANT_CORNER_RULE))")
    return length(found)
end

# A sidewall of an invariant corner: the vertical CAD plane through the corner along
# one of its sides. A mesh vertex lies on wall i when a VERTICAL surface triangle
# containing it lies in that plane on the side's half (a thin coupon has no vertical
# triangle at a kink, so no wall and no bridging sliver there; a box-face triangle is
# vertical but never in a metal side's plane). Returns per vertex a 2-bit wall mask.
const SIDEWALL_PLANE_TOLERANCE = 1.0e-8

function sidewall_vertex_masks(points, triangles, center, sides, radius)
    directions = sides.walls
    normals = [[-d[2], d[1]] for d in directions]
    in_wall(p, w) = abs((p[1] - center[1]) * normals[w][1] + (p[2] - center[2]) * normals[w][2]) <=
                    SIDEWALL_PLANE_TOLERANCE &&
                    (p[1] - center[1]) * directions[w][1] + (p[2] - center[2]) * directions[w][2] >=
                    -SIDEWALL_PLANE_TOLERANCE
    # The vertices of the vertical triangles lying in a wall plane near the corner ...
    wall_vertices = Set{Int}()
    for triangle in triangles
        xyz = [points[:, i] for i in triangle]
        any(norm(p .- center) <= radius for p in xyz) || continue
        n = cross(xyz[2] .- xyz[1], xyz[3] .- xyz[1])
        length_n = norm(n)
        length_n > 0.0 && abs(n[3]) <= 1.0e-6 * length_n || continue
        any(all(in_wall(p, w) for p in xyz) for w in eachindex(directions)) || continue
        union!(wall_vertices, triangle)
    end
    # ... each masked by the wall planes it lies in (a corner-line vertex lies in both,
    # whichever wall's triangles happen to contain it).
    masks = Dict{Int, UInt8}()
    for i in wall_vertices
        mask = 0x00
        for w in eachindex(directions)
            in_wall(points[:, i], w) && (mask |= UInt8(1 << (w - 1)))
        end
        masks[i] = mask
    end
    return masks
end

# FIX-UPS F.1 (b) as ruled by supervisor decision 365: a corner-incident cell is a
# BRIDGING SLIVER candidate when all four of its vertices lie on the kink's two
# SIDEWALLS only (none interior, none on a third CAD surface), at least one strictly on
# wall 1 (off wall 2) and one strictly on wall 2, in a wedge whose opening at the corner
# is OBTUSE (> 90 degrees exactly) - the S4 150-degree fab cell: the corner and a
# corner-line vertex on both walls, one vertex on each wall, in the 150-degree wedge
# under the metal (kappa_reg 5.28-6.16). The same combinatorial type is the ideal
# trirectangular cell at 90 degrees and an unavoidable healthy wedge cell at a
# fabricated tip (the substrate fan under a 22.5-degree tip, kappa_reg 3.4-4.2), so the
# predicate is NOT a verdict: the verdict is kappa_reg <= CornerShapeGate, and the
# predicate restricted to cells still above the gate after the kappa_reg descent is the
# TRIGGER of the corner-local reconnection pass (F5-B).
function bridging_sliver_cells(points, tetrahedra, cells, masks, center, sides)
    bridging = Int[]
    for k in cells
        cell = tetrahedra[k]
        all(get(masks, i, 0x00) != 0x00 for i in cell) || continue
        wall_1 = any(get(masks, i, 0x00) == 0x01 for i in cell)
        wall_2 = any(get(masks, i, 0x00) == 0x02 for i in cell)
        wall_1 && wall_2 || continue
        centroid = sum(points[:, i] for i in cell) ./ 4
        wedge_opening(sides, [centroid[1] - center[1], centroid[2] - center[2]]) > pi / 2 || continue
        push!(bridging, k)
    end
    return bridging
end

# Edge lengths and cell aspects of the linear seed inside each semantic corner
# ball, and per shell of the corner law (corner_shell_radii: the geometric shell
# boundaries, then the radius) the ball edges by midpoint distance, their length
# percentiles against the law's size at the shell, and the cells by centroid.
function corner_shell_census(points, edges, cells, center, grading::CornerGrading)
    boundaries = vcat(0.0, corner_shell_radii(grading))
    rows = Dict{String, Any}[]
    for (inner, outer) in zip(boundaries[1:(end - 1)], boundaries[2:end])
        midpoint(a, b) = norm(0.5 .* (points[:, a] .+ points[:, b]) .- center)
        lengths = sort!([norm(points[:, a] .- points[:, b]) for (a, b) in edges
                         if inner < midpoint(a, b) <= outer])
        centroid_count = count(cell -> inner < norm(sum(points[:, i] for i in cell) ./ 4 .- center) <= outer,
                               cells)
        target = corner_ball_size(grading, inner)
        push!(rows, Dict{String, Any}(
            "InnerRadius" => inner, "OuterRadius" => outer, "TargetSize" => target,
            "Edges" => length(lengths), "Cells" => centroid_count,
            "EdgeP50" => isempty(lengths) ? nothing : sorted_median(lengths),
            "EdgeP90" => isempty(lengths) ? nothing :
                         lengths[clamp(ceil(Int, 0.9 * length(lengths)), 1, length(lengths))],
            "EdgeMaximum" => isempty(lengths) ? nothing : lengths[end],
            "EdgesOverSqrt2TargetSize" => count(>(sqrt(2.0) * target), lengths)))
    end
    return rows
end

function seed_corner_census(corners, grading::CornerGrading, tolerance;
                            corner_kinds=fill(:legacy, length(corners)),
                            maximum_corner_aspect=0.0, corner_shape_gate=0.0)
    radius, isotropic_size = grading.radius, grading.lc_fine
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    types, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
    tetrahedra = Vector{NTuple{4, Int}}()
    for (type, tags, block) in zip(types, element_tags, element_nodes)
        isempty(tags) && continue
        # Tetrahedral blocks only (a tube's prisms and pyramids are not corner-ball
        # cells); only the vertex nodes define the linear cells (high-order nodes follow).
        gmsh.model.mesh.getElementProperties(type)[1] |> startswith("Tetrahedron") || continue
        nodes_per_element = length(block) ÷ length(tags)
        for start in 1:nodes_per_element:length(block)
            push!(tetrahedra, ntuple(i -> index[block[start + i - 1]], 4))
        end
    end
    threshold = sqrt(2.0) * isotropic_size
    rows = Dict{String, Any}[]
    for (k, corner) in enumerate(corners)
        center = collect(corner)
        inside = [norm(points[:, i] .- center) <= radius for i in axes(points, 2)]
        at_corner = [norm(points[:, i] .- center) <= tolerance for i in axes(points, 2)]
        ball = [cell for cell in tetrahedra if all(inside[i] for i in cell)]
        incident = [cell for cell in tetrahedra if any(at_corner[i] for i in cell)]
        edges = Set{Tuple{Int, Int}}()
        for cell in ball, i in 1:4, j in (i + 1):4
            push!(edges, (min(cell[i], cell[j]), max(cell[i], cell[j])))
        end
        lengths = sort!([norm(points[:, a] .- points[:, b]) for (a, b) in edges])
        aspects = [tetrahedron_aspect([points[:, i] for i in cell]) for cell in incident]
        ring = unique([i for cell in incident for i in cell if !at_corner[i]])
        ring_radii = [norm(points[:, i] .- center) for i in ring]
        # Design round 2 F5-A: the corner's kind and the measure of its verdict
        # (kappa_reg at an invariant corner, the vertex-0 condition at a legacy one), with
        # the order-invariant measures of every corner as information (FIX-UPS F.6).
        invariant = corner_kinds[k] === :invariant
        cells = [[points[:, i] for i in cell] for cell in incident]
        frames = [tetrahedron_frame_conditions(xyz) for xyz in cells]
        kappa_reg = [tetrahedron_regular_condition(xyz) for xyz in cells]
        measure = isempty(cells) ? nothing : Dict{String, Any}(
            "Kind" => invariant ? "Invariant" : "Legacy",
            "Name" => invariant ? "RegularCondition" : "VertexFrameCondition",
            "Value" => invariant ? maximum(kappa_reg) : maximum(aspects),
            "Gate" => invariant ? (corner_shape_gate > 0.0 ? corner_shape_gate : nothing) :
                      (maximum_corner_aspect > 0.0 ? maximum_corner_aspect : nothing),
            "KappaRegMax" => maximum(kappa_reg), "KappaV0Max" => maximum(aspects),
            "KappaMinMax" => maximum(first(f) for f in frames),
            "KappaMaxMax" => maximum(last(f) for f in frames),
            "EtaMin" => minimum(tetrahedron_mean_ratio(xyz) for xyz in cells))
        push!(rows, Dict{String, Any}(
            "Corner" => k - 1, "Point" => collect(corner), "Measure" => measure,
            "BallCells" => length(ball), "BallEdges" => length(lengths),
            "EdgeMinimum" => isempty(lengths) ? nothing : lengths[1],
            "EdgeMedian" => isempty(lengths) ? nothing : sorted_median(lengths),
            "EdgeMaximum" => isempty(lengths) ? nothing : lengths[end],
            "EdgesOverSqrt2IsotropicSize" => count(>(threshold), lengths),
            "FractionOverSqrt2IsotropicSize" =>
                isempty(lengths) ? nothing : count(>(threshold), lengths) / length(lengths),
            "IncidentCells" => length(incident),
            "IncidentMaximumAspect" => isempty(aspects) ? nothing : maximum(aspects),
            "RingRadiusMinimum" => isempty(ring_radii) ? nothing : minimum(ring_radii),
            "RingRadiusMaximum" => isempty(ring_radii) ? nothing : maximum(ring_radii),
            "RingEdgesOverSqrt2IsotropicSize" => count(>(threshold), ring_radii),
            "Shells" => corner_shell_census(points, edges, ball, center, grading)))
    end
    return rows
end

# ---------------------------------------------------------------------------
# Seed-side quality optimization of the MMG required region (supervisor decision
# 30). The metric stage marks the seed cells inside the semantic corner balls
# (centroid within CornerIsotropyRadius) and the edge layer (a vertex within
# LayerThickness x (1 + RowZigzag) + EdgeSize of a span) as MMG required
# tetrahedra, so their quality after adaptation is exactly the seed's: the seed
# must satisfy the corner-aspect and scaled-Jacobian gates on them itself. The
# repair is the restorer's bounded rule on the seed: vertices move within
# ratio x the local prescribed size (lc_fine, EdgeSize-graded in the layer),
# surface vertices stay on their planes (null space of the incident triangle
# normals; three independent planes fix a vertex), every touched cell keeps a
# positive orientation and a scaled Jacobian of at least min(original, target).

function tetrahedron_scaled_jacobian(xyz)
    a = xyz[2] .- xyz[1]; b = xyz[3] .- xyz[1]; c = xyz[4] .- xyz[1]
    return dot(a, cross(b, c)) / (norm(a) * norm(b) * norm(c))
end

# Longest edge over shortest height (altitude): the edge-layer quality rule's
# aspect (edge_volume_metric.tetrahedron_edge_aspect; the shortest height is
# |det| / (2 x largest face area)). Infinite for a flat cell.
function tetrahedron_edge_aspect(xyz)
    a = xyz[2] .- xyz[1]; b = xyz[3] .- xyz[1]; c = xyz[4] .- xyz[1]
    determinant = abs(dot(a, cross(b, c)))
    longest = maximum(norm(xyz[i] .- xyz[j]) for i in 1:4 for j in (i + 1):4)
    area = maximum(0.5 * norm(cross(xyz[q] .- xyz[p], xyz[r] .- xyz[p]))
                   for (p, q, r) in ((2, 3, 4), (1, 3, 4), (1, 2, 4), (1, 2, 3)))
    return determinant > 0.0 ? longest * 2.0 * area / determinant : Inf
end

# Scaled Jacobian below which a layer cell is flat to roundoff
# (edge_volume_metric.EDGE_LAYER_ORIENTATION_FLOOR).
const EDGE_LAYER_ORIENTATION_FLOOR = 1.0e-12

# 3D point-to-segment distance (the plan-view point_segment_distance below is 2D).
function span_point_distance(point, start, stop)
    v = stop .- start
    span_length = norm(v)
    span_length > 0.0 || error("Degenerate edge layer span")
    axial = clamp(dot(point .- start, v) / span_length, 0.0, span_length)
    return norm(point .- start .- (axial / span_length) .* v)
end

# The metric stage's REQUIRED_TETRAHEDRA_RULE on the seed (same regions, same
# reach); the third result is the layer part (EDGE_LAYER_CELL_RULE: a vertex
# within the reach of a span).
function required_region_cells(points, tetrahedra, corners, radius, spans, reach)
    required = falses(length(tetrahedra))
    in_layer = falses(length(tetrahedra))
    span_distance = fill(Inf, size(points, 2))
    if !isempty(spans)
        for i in axes(points, 2), (start, stop) in spans
            span_distance[i] = min(span_distance[i], span_point_distance(points[:, i], start, stop))
        end
    end
    for (k, cell) in enumerate(tetrahedra)
        centroid = sum(points[:, i] for i in cell) ./ 4
        in_layer[k] = !isempty(spans) && any(span_distance[i] <= reach for i in cell)
        if in_layer[k] || any(norm(centroid .- collect(corner)) <= radius for corner in corners)
            required[k] = true
        end
    end
    return required, span_distance, in_layer
end

# Movement basis of every surface vertex: null space of its triangle normals and,
# for a vertex on a CAD curve (a line element), of the plane orthogonal to the
# curve: the vertex stays on the curve.  A curve between two coplanar faces (the
# thin sheet edge between the metal sheet and the SA plane, decision 66) is
# invisible to the triangle normals alone, and a sheet-edge vertex moved in-plane
# would change the metal footprint.
function surface_movement_bases(points, triangles, lines=NTuple{2, Int}[])
    normals = Dict{Int, Vector{Vector{Float64}}}()
    for triangle in triangles
        a, b, c = (points[:, i] for i in triangle)
        n = cross(b .- a, c .- a)
        length_n = norm(n)
        length_n > 0.0 || continue
        for i in triangle
            push!(get!(normals, i, Vector{Float64}[]), n ./ length_n)
        end
    end
    for line in lines
        a, b = (points[:, i] for i in line)
        d = b .- a
        length_d = norm(d)
        length_d > 0.0 || continue
        complement = nullspace(reshape(d ./ length_d, 1, 3))
        for i in line, column in eachcol(complement)
            push!(get!(normals, i, Vector{Float64}[]), collect(column))
        end
    end
    bases = Dict{Int, Matrix{Float64}}()
    for (node, rows) in normals
        decomposition = svd(reduce(hcat, rows)'; full=true)
        rank = count(>(1.0e-6), decomposition.S)
        bases[node] = decomposition.V[:, (rank + 1):end]
    end
    return bases
end

# The aspect objective descends on the p-norm of the target aspects (a smooth
# proxy of the maximum that keeps descending where several cells tie for the
# worst; the gate/target are always judged on the true maximum).
const SEED_ASPECT_PROXY_POWER = 32

# Every nonzero {-1, 0, 1} combination of the basis columns, normalized: 26
# directions for a free vertex, 8 in a plane, 2 on a line.
function descent_directions(basis)
    k = size(basis, 2)
    directions = Vector{Float64}[]
    for code in 0:(3^k - 1)
        coefficients = [Float64(((code ÷ 3^(j - 1)) % 3) - 1) for j in 1:k]
        all(iszero, coefficients) && continue
        direction = basis * coefficients
        push!(directions, direction ./ norm(direction))
    end
    return directions
end

# Greedy bounded coordinate descent on one target set: minimize the maximum
# aspect (:aspect, the Jacobian condition number; :edge_aspect, the layer rule's
# longest edge over shortest height; both through their p-norm proxy) or raise
# the minimum scaled Jacobian (:scaled) of the target cells until the goal is met
# on the true maximum/minimum. Every touched cell keeps its scaled-Jacobian floor
# and, with `ceilings` (per cell, Inf = none), its edge-aspect ceiling: a repair
# may never flatten a neighbouring layer cell; with `condition_ceilings` (per
# cell) its Jacobian-condition ceiling, and with `edge_floors` (per vertex, the
# sub-size collapse threshold) no edge of the moved vertex may become shorter
# than min(its original length, the floor): a repair may never make the
# sub-size/high-condition cells the collapse removes (measured on EL1c: the
# descent moved face vertices to 0.013-0.24 nm of a row or ridge node, 43-56
# cells above condition 1000; supervisor decision 34). Returns (achieved true
# value, accepted moves).
function optimize_seed_cells!(points, original, tetrahedra, incident, bases, floors, bounds,
                              targets, objective::Symbol, goal; ceilings=nothing,
                              condition_ceilings=nothing, edge_floors=nothing)
    function cell_xyz(cell)
        return [points[:, i] for i in cell]
    end
    # :regular_condition is the invariant corner objective (kappa_reg, design round 2
    # F5-A) through the same p-norm proxy as :aspect.
    minimizing = objective in (:aspect, :edge_aspect, :regular_condition)
    aspect_of = objective === :edge_aspect ? tetrahedron_edge_aspect :
                objective === :regular_condition ? tetrahedron_regular_condition : tetrahedron_aspect
    function value()
        if minimizing
            aspects = [aspect_of(cell_xyz(tetrahedra[k])) for k in targets]
            scale = maximum(aspects)
            return scale * sum((aspects ./ scale) .^ SEED_ASPECT_PROXY_POWER)^(1 / SEED_ASPECT_PROXY_POWER)
        end
        return -minimum(tetrahedron_scaled_jacobian(cell_xyz(tetrahedra[k])) for k in targets)
    end
    function achieved()
        if minimizing
            return maximum(aspect_of(cell_xyz(tetrahedra[k])) for k in targets)
        end
        return -value()
    end
    vertices = unique(i for k in targets for i in tetrahedra[k])
    movable = [(i, get(bases, i, Matrix{Float64}(I, 3, 3))) for i in vertices]
    filter!(pair -> size(pair[2], 2) > 0, movable)
    current = value()
    moves = 0
    for _ in 1:60
        (minimizing ? achieved() <= goal : achieved() >= goal) && break
        improved = false
        for (vertex, basis) in movable
            cells = incident[vertex]
            bound = bounds[vertex]
            neighbours = unique(j for k in cells for j in tetrahedra[k] if j != vertex)
            edge_floor = edge_floors === nothing ? 0.0 : edge_floors[vertex]
            edge_limits = [min(norm(points[:, j] .- points[:, vertex]), edge_floor)
                           for j in neighbours]
            for direction in descent_directions(basis)
                step = bound / 2
                while step > bound / 128
                    previous = points[:, vertex]
                    points[:, vertex] = previous .+ step .* direction
                    accepted = norm(points[:, vertex] .- original[:, vertex]) <=
                               bound * (1.0 + 1.0e-12)
                    if accepted && edge_floor > 0.0
                        accepted = all(norm(points[:, j] .- points[:, vertex]) >=
                                       limit * (1.0 - 1.0e-9)
                                       for (j, limit) in zip(neighbours, edge_limits))
                    end
                    if accepted
                        for k in cells
                            xyz = cell_xyz(tetrahedra[k])
                            scaled = tetrahedron_scaled_jacobian(xyz)
                            if !(scaled > 0.0) || scaled < floors[k] - 1.0e-9 ||
                               (ceilings !== nothing && isfinite(ceilings[k]) &&
                                tetrahedron_edge_aspect(xyz) > ceilings[k] * (1.0 + 1.0e-9)) ||
                               (condition_ceilings !== nothing && isfinite(condition_ceilings[k]) &&
                                tetrahedron_aspect(xyz) > condition_ceilings[k] * (1.0 + 1.0e-9))
                                accepted = false
                                break
                            end
                        end
                    end
                    if accepted
                        candidate = value()
                        if candidate < current - 1.0e-12
                            current = candidate
                            improved = true
                            moves += 1
                            break
                        end
                    end
                    points[:, vertex] = previous
                    step /= 2
                end
            end
        end
        improved || break
    end
    return achieved(), moves
end

# ---------------------------------------------------------------------------
# F5-B (mesher design round 2 section 1.4 (b) / FIX-UPS F.1 (b); supervisor decisions
# 363 / 365): the corner-local reconnection pass for the fabricated flat-kink sliver
# (family-5 root cause B: a corner-incident cell with its four vertices on the two
# sidewalls, kappa_reg 5.28-6.16 at 150-151 degrees, which the bounded descent cannot
# repair - the vertices are confined to the walls). Trigger = a bridging-sliver
# candidate still above CornerShapeGate after the kappa_reg descent (decision 365).
# Each trigger cell takes the better of two operations: the 2-3 flip across one of its
# INTERIOR faces (a face that is no surface triangle, shared with a cell of the same
# volume entity; the union of the two cells convex so the three new cells are positively
# oriented) and the EDGE REMOVAL of one of its interior edges (the closed ring of n >= 3
# cells of one entity around the edge re-triangulated by every fan: 2 (n - 2) cells; the
# 3-2 flip for n = 3 - the operation every measured S4-type sliver used, whose two
# interior faces have non-convex unions so no 2-3 flip applies), each by the smallest
# resulting maximum kappa_reg (the edge removal preferred when its quality is <= the
# flip's), which must improve on the trigger cell; every new cell keeps the TRIGGER
# cell's scaled-Jacobian floor (min(its original value, 2 x MinimumScaledJacobian)) and
# its Jacobian-condition ceiling (max(its original value, JacobianConditionTarget));
# then the descent runs again. Boundary triangles are never flipped; the moves of the
# descent stay bounded as before. Recorded CornerReconnections per corner.
const CORNER_RECONNECTION_RULE =
    "a bridging-sliver candidate above CornerShapeGate after the kappa_reg descent (a " *
    "corner-incident cell with its four vertices on the two sidewalls, one strictly on each, " *
    "in an obtuse wedge) is replaced by the better of its 2-3 flip across one of its interior " *
    "faces (no surface triangle; the neighbour in the same volume entity; the three new cells " *
    "positively oriented) and its edge removal of one of its interior edges (no surface edge; " *
    "the closed ring of n >= 3 cells of one volume entity around the edge re-triangulated by " *
    "every fan into 2 (n - 2) positively oriented cells; the 3-2 flip for n = 3), each by the " *
    "smallest resulting maximum kappa_reg (the edge removal preferred when its quality is <= " *
    "the flip's), which must improve on the cell; every new cell keeps the trigger cell's " *
    "floor min(its original scaled Jacobian, 2 x MinimumScaledJacobian) and its ceiling " *
    "max(its original Jacobian condition, JacobianConditionTarget); the descent then runs " *
    "again; up to CornerReconnectionRounds rounds per corner, a trigger slot an earlier " *
    "reconnection of the round reused being skipped until the next round (mesher design " *
    "round 2 F5-B, decisions 363 / 365 / 392)"
const CORNER_RECONNECTION_ROUNDS = 4

function oriented_tetrahedron(points, cell)
    xyz = [points[:, i] for i in cell]
    volume = dot(xyz[2] .- xyz[1], cross(xyz[3] .- xyz[1], xyz[4] .- xyz[1]))
    volume == 0.0 && return nothing
    return volume > 0.0 ? cell : (cell[1], cell[2], cell[4], cell[3])
end

# The best 2-3 flip of cell k: (neighbour j, the three new cells, their maximum
# kappa_reg) or nothing.
function best_two_three_flip(points, tetrahedra, incident, surface_faces, cell_entity, k,
                             scaled_floor, condition_ceiling)
    cell = tetrahedra[k]
    best = nothing
    for opposite in 1:4
        face = Tuple(sort!([cell[i] for i in 1:4 if i != opposite]))
        face in surface_faces && continue
        apex_k = cell[opposite]
        neighbours = [j for j in incident[face[1]]
                      if j != k && all(v in tetrahedra[j] for v in face)]
        length(neighbours) == 1 || continue
        j = neighbours[1]
        cell_entity(j) == cell_entity(k) || continue
        apex_j = only(v for v in tetrahedra[j] if !(v in face))
        candidates = NTuple{4, Int}[]
        valid = true
        for (a, b) in ((face[1], face[2]), (face[2], face[3]), (face[3], face[1]))
            oriented = oriented_tetrahedron(points, (apex_k, apex_j, a, b))
            oriented === nothing && (valid = false; break)
            push!(candidates, oriented)
        end
        valid || continue
        # The three cells fill the union only when the segment apex_k - apex_j pierces
        # the face: the signed volumes of the three cells then have one sign and sum to
        # the pair's volume.
        volume(c) = (xyz = [points[:, i] for i in c];
                     dot(xyz[2] .- xyz[1], cross(xyz[3] .- xyz[1], xyz[4] .- xyz[1])) / 6.0)
        pair = abs(volume(cell)) + abs(volume(tetrahedra[j]))
        total = sum(abs(volume(c)) for c in candidates)
        quality = maximum(tetrahedron_regular_condition([points[:, i] for i in c])
                          for c in candidates)
        scaled = minimum(tetrahedron_scaled_jacobian([points[:, i] for i in c])
                         for c in candidates)
        condition = maximum(tetrahedron_aspect([points[:, i] for i in c]) for c in candidates)
        abs(total - pair) <= 1.0e-9 * pair || continue
        scaled >= scaled_floor * (1.0 - 1.0e-9) || continue
        condition <= condition_ceiling * (1.0 + 1.0e-9) || continue
        if best === nothing || quality < best.quality
            best = (neighbour=j, cells=candidates, quality=quality, face=face)
        end
    end
    return best
end

# The best edge removal of cell k: for an interior edge (no surface edge) with n >= 3
# cells around it in one closed ring of the same volume entity, the ring polygon is
# re-triangulated (every fan of the ring) and each ring triangle gives the two cells
# with the edge's ends (the 3-2 flip for n = 3, the 4-4 flip for n = 4, ...); the
# candidate is valid when every new cell is positively oriented and the volumes sum to
# the ring's. Returns (ring cells, new cells, quality, edge) or nothing.
function best_edge_removal(points, tetrahedra, incident, surface_edges, cell_entity, k,
                           scaled_floor, condition_ceiling)
    cell = tetrahedra[k]
    best = nothing
    volume(c) = (xyz = [points[:, i] for i in c];
                 dot(xyz[2] .- xyz[1], cross(xyz[3] .- xyz[1], xyz[4] .- xyz[1])) / 6.0)
    for a in 1:4, b in (a + 1):4
        u, v = minmax(cell[a], cell[b])
        (u, v) in surface_edges && continue
        ring = [j for j in incident[u] if v in tetrahedra[j]]
        length(ring) >= 3 || continue
        all(cell_entity(j) == cell_entity(k) for j in ring) || continue
        others = Dict(j => Tuple(w for w in tetrahedra[j] if w != u && w != v) for j in ring)
        all(length(others[j]) == 2 for j in ring) || continue
        # Chain the ring: consecutive cells share one ring vertex; a closed single cycle.
        cycle = Int[]
        current = ring[1]
        vertex = others[current][1]
        visited = Set{Int}()
        closed = true
        for _ in 1:length(ring)
            push!(cycle, vertex); push!(visited, current)
            vertex = others[current][1] == vertex ? others[current][2] : others[current][1]
            next = [j for j in ring if !(j in visited) && vertex in others[j]]
            if isempty(next)
                closed = (length(visited) == length(ring)) && vertex == cycle[1]
                break
            end
            length(next) == 1 || (closed = false; break)
            current = next[1]
        end
        closed && length(cycle) == length(ring) && length(unique(cycle)) == length(ring) || continue
        n = length(cycle)
        ring_volume = sum(abs(volume(tetrahedra[j])) for j in ring)
        for pivot in 1:n
            triangles = [(cycle[pivot], cycle[mod1(pivot + i, n)], cycle[mod1(pivot + i + 1, n)])
                         for i in 1:(n - 2)]
            candidates = NTuple{4, Int}[]
            valid = true
            for (wa, wb, wc) in triangles, apex in (u, v)
                oriented = oriented_tetrahedron(points, (apex, wa, wb, wc))
                oriented === nothing && (valid = false; break)
                push!(candidates, oriented)
            end
            valid || continue
            total = sum(abs(volume(c)) for c in candidates)
            quality = maximum(tetrahedron_regular_condition([points[:, i] for i in c])
                              for c in candidates)
            scaled = minimum(tetrahedron_scaled_jacobian([points[:, i] for i in c])
                             for c in candidates)
            condition = maximum(tetrahedron_aspect([points[:, i] for i in c])
                                for c in candidates)
            abs(total - ring_volume) <= 1.0e-9 * ring_volume || continue
            scaled >= scaled_floor * (1.0 - 1.0e-9) || continue
            condition <= condition_ceiling * (1.0 + 1.0e-9) || continue
            if best === nothing || quality < best.quality
                best = (ring=ring, cells=candidates, quality=quality, edge=(u, v))
            end
        end
    end
    return best
end

# Below-target cells of a target set grouped into vertex-sharing components.
function below_target_components(cells, tetrahedra, incident, value_of, target)
    bad = Set(k for k in cells if value_of(k) < target)
    components = Vector{Vector{Int}}()
    while !isempty(bad)
        first = pop!(bad)
        component = [first]; stack = [first]
        while !isempty(stack)
            k = pop!(stack)
            for i in tetrahedra[k], j in incident[i]
                if j in bad
                    delete!(bad, j); push!(component, j); push!(stack, j)
                end
            end
        end
        push!(components, sort!(component))
    end
    return components
end

# Tolerance (um) within which a collapse target counts as lying in a plane of the
# collapsed surface vertex (its triangle normals are exact to roundoff on the
# planar CAD faces and on the seeded rows).
const SEED_COLLAPSE_PLANE_TOLERANCE = 1.0e-8

# An edge is a sub-size insertion when shorter than this fraction of the smallest
# prescribed size: the seeded structure itself has edges at the size (ridge to
# first row = EdgeSize, to roundoff) and Gmsh fills between the rows with edges
# down to 0.98 x the size, while the measured artefacts (a Gmsh volume vertex
# 0.3 nm from a 1 nm row node; face vertices 0.02-0.24 nm from a 1 nm row node;
# a 0.41 nm edge at a 1 nm corner shell; 0.48 nm at a 4 nm layer) are all below
# 0.41 x the size. A threshold at the size collapsed 14,667 legitimate layer face
# vertices and 2,677 row nodes on the EL1c seed.
const SEED_COLLAPSE_SIZE_FRACTION = 0.5

# Distinct unit normals of the triangles incident to every surface vertex (the
# planes the vertex lies in).
function surface_vertex_normals(points, triangles)
    normals = Dict{Int, Vector{Vector{Float64}}}()
    for triangle in triangles
        a, b, c = (points[:, i] for i in triangle)
        n = cross(b .- a, c .- a)
        length_n = norm(n)
        length_n > 0.0 || continue
        n = n ./ length_n
        for i in triangle
            rows = get!(normals, i, Vector{Float64}[])
            any(abs(dot(row, n)) > 1.0 - 1.0e-9 for row in rows) || push!(rows, n)
        end
    end
    return normals
end

# Collapse the seed's own sub-size insertions inside the required region: a
# candidate vertex (within the layer reach or, with corner grading, inside a
# corner ball) whose shortest incident edge is below `threshold` -
# SEED_COLLAPSE_SIZE_FRACTION x the smallest prescribed size (EdgeSize, the
# adapter hmin, or CornerSize), so such an edge is below the prescribed metric
# everywhere - is merged onto one of its neighbours. A free interior vertex (in no surface triangle) may merge onto any
# neighbour; a surface vertex only along its own surface: onto a surface
# neighbour lying in every plane of the vertex (within SEED_COLLAPSE_PLANE_TOLERANCE)
# whose line/triangle supports (`supports`: the (dimension, entity) pairs of the
# incident line and triangle elements) include the vertex's own, so a vertex on
# one face moves within that face, a vertex on a ridge or a seeded row along it,
# and the remapped triangles keep their orientation; `fixed` vertices (CAD
# points, the semantic corners among them) and vertices in three independent
# planes are never collapsed. Among the valid cavities (every remapped cell
# positively oriented above EDGE_LAYER_ORIENTATION_FLOOR, every remapped triangle
# keeping its normal) the one with the smallest maximum edge aspect is taken,
# accepted only when the cavity is no worse than the cells it replaces in maximum
# edge aspect (the restorer's corner-ball collapse rule, decision 15a) and in
# maximum Jacobian condition, better in one of them, and keeps every cell at or
# above `scaled_jacobian_gate` (the gate judging the cells; 0 under the layer
# quality rule) unless a replaced cell was already below it.
# The bounded descent then repairs the cavity further. Bounded moves alone cannot
# repair the cells around such a vertex (measured: a Gmsh volume vertex 0.3 nm
# from a row node and 0.04 nm under the face plane, 13 cells at aspect 1150,
# every direction inverting a cell or staying above 300; the collapse onto the
# row node gives a worst cell of 235; a Gmsh face vertex 0.02-0.24 nm from a 1 nm
# row node on the metal faces, 56 layer cells above Jacobian condition 1000 -
# supervisor decision 34). Mutates `tetrahedra`, `triangles` and `lines` in place
# (elements containing both vertices are deleted, the others remapped) and
# returns a NamedTuple: per dimension the deleted original indices and
# Dict(original index => remapped element), the collapse records (position,
# target, surface flag, edge length, replaced/cavity quality) and the shortest
# collapsed edge.
function collapse_short_edges!(points, tetrahedra, triangles, lines, candidate, fixed,
                               supports, threshold; scaled_jacobian_gate::Float64)
    surface = falses(size(points, 2))
    for triangle in triangles, i in triangle
        surface[i] = true
    end
    normals = surface_vertex_normals(points, triangles)
    incident = [Int[] for _ in axes(points, 2)]
    for (k, cell) in enumerate(tetrahedra), i in cell
        push!(incident[i], k)
    end
    incident_triangles = [Int[] for _ in axes(points, 2)]
    for (t, triangle) in enumerate(triangles), i in triangle
        push!(incident_triangles[i], t)
    end
    incident_lines = [Int[] for _ in axes(points, 2)]
    for (l, line) in enumerate(lines), i in line
        push!(incident_lines[i], l)
    end
    deleted = Set{Int}(); remapped = Dict{Int, NTuple{4, Int}}()
    deleted_triangles = Set{Int}(); remapped_triangles = Dict{Int, NTuple{3, Int}}()
    deleted_lines = Set{Int}(); remapped_lines = Dict{Int, NTuple{2, Int}}()
    current(k) = get(remapped, k, tetrahedra[k])
    current_triangle(t) = get(remapped_triangles, t, triangles[t])
    current_line(l) = get(remapped_lines, l, lines[l])
    triangle_normal(triangle) = cross(points[:, triangle[2]] .- points[:, triangle[1]],
                                      points[:, triangle[3]] .- points[:, triangle[1]])
    records = Dict{String, Any}[]; shortest = Inf
    # The least constrained vertices go first (a free Gmsh face vertex is merged
    # onto its row node before the row node is considered): by support count.
    candidates = sort!([i for i in axes(points, 2)
                        if candidate[i] && !fixed[i] && !isempty(incident[i]) &&
                           (!surface[i] || length(get(normals, i, Vector{Float64}[])) < 3)];
                       by=i -> (length(supports[i]), i))
    for v in candidates
        cells = [k for k in incident[v] if !(k in deleted)]
        isempty(cells) && continue
        neighbours = unique(j for k in cells for j in current(k) if j != v)
        lengths = [norm(points[:, j] .- points[:, v]) for j in neighbours]
        shortest_length = minimum(lengths)
        shortest_length < threshold || continue
        faces = [t for t in incident_triangles[v] if !(t in deleted_triangles)]
        curves = [l for l in incident_lines[v] if !(l in deleted_lines)]
        replaced_xyz = [[points[:, i] for i in current(k)] for k in cells]
        replaced_worst = maximum(tetrahedron_edge_aspect(xyz) for xyz in replaced_xyz)
        replaced_scaled = minimum(tetrahedron_scaled_jacobian(xyz) for xyz in replaced_xyz)
        replaced_condition = maximum(tetrahedron_aspect(xyz) for xyz in replaced_xyz)
        best = nothing; best_aspect = Inf; best_scaled = 0.0; best_condition = Inf; w = 0
        best_faces = Dict{Int, NTuple{3, Int}}()
        for target in neighbours
            if surface[v]
                # Along the vertex's own surface: a surface target in every plane of
                # v carrying every support of v.
                surface[target] || continue
                fixed[target] && continue
                all(abs(dot(n, points[:, target] .- points[:, v])) <= SEED_COLLAPSE_PLANE_TOLERANCE
                    for n in normals[v]) || continue
                issubset(supports[v], supports[target]) || continue
            end
            candidate_cells = Dict{Int, NTuple{4, Int}}()
            valid = true; worst = 0.0; least = Inf; condition = 0.0
            for k in cells
                cell = current(k)
                target in cell && continue
                replaced = ntuple(i -> cell[i] == v ? target : cell[i], 4)
                xyz = [points[:, i] for i in replaced]
                scaled = tetrahedron_scaled_jacobian(xyz)
                if !(scaled > EDGE_LAYER_ORIENTATION_FLOOR)
                    valid = false
                    break
                end
                worst = max(worst, tetrahedron_edge_aspect(xyz))
                least = min(least, scaled)
                condition = max(condition, tetrahedron_aspect(xyz))
                candidate_cells[k] = replaced
            end
            valid || continue
            candidate_faces = Dict{Int, NTuple{3, Int}}()
            for t in faces
                triangle = current_triangle(t)
                target in triangle && continue
                replaced = ntuple(i -> triangle[i] == v ? target : triangle[i], 3)
                before = triangle_normal(triangle); after = triangle_normal(replaced)
                if !(norm(after) > 0.0 && dot(before, after) > 0.0)
                    valid = false
                    break
                end
                candidate_faces[t] = replaced
            end
            valid || continue
            # A collapse must leave remapped cells covering the cavity, no worse than
            # the cells it replaces in edge aspect and Jacobian condition and better
            # in one of them (the vertex's worst cell is among the deleted ones: a
            # remapped cell is the replaced cell shifted by the collapsed edge, so
            # without an improvement the comparison would be roundoff), and, under
            # the scaled-Jacobian gate, no cell below the gate where none was.
            isempty(candidate_cells) && continue
            worst <= replaced_worst * (1.0 + 1.0e-9) || continue
            condition <= replaced_condition * (1.0 + 1.0e-9) || continue
            (worst < replaced_worst * (1.0 - 1.0e-9) ||
             condition < replaced_condition * (1.0 - 1.0e-9)) || continue
            (replaced_scaled < scaled_jacobian_gate ||
             least >= scaled_jacobian_gate * (1.0 - 1.0e-9)) || continue
            if worst < best_aspect
                best, best_aspect, best_scaled, best_condition, w = candidate_cells, worst, least,
                                                                    condition, target
                best_faces = candidate_faces
            end
        end
        best === nothing && continue
        for k in cells
            if w in current(k)
                push!(deleted, k)
            else
                remapped[k] = best[k]
                push!(incident[w], k)
            end
        end
        for t in faces
            if w in current_triangle(t)
                push!(deleted_triangles, t)
            else
                remapped_triangles[t] = best_faces[t]
                push!(incident_triangles[w], t)
            end
        end
        for l in curves
            line = current_line(l)
            if w in line
                push!(deleted_lines, l)
            else
                remapped_lines[l] = ntuple(i -> line[i] == v ? w : line[i], 2)
                push!(incident_lines[w], l)
            end
        end
        shortest = min(shortest, shortest_length)
        push!(records, Dict{String, Any}(
            "Vertex" => v, "Position" => points[:, v], "Target" => w,
            "TargetPosition" => points[:, w], "Surface" => surface[v],
            "EdgeLength" => shortest_length, "Cells" => length(cells),
            "Triangles" => length(faces), "Lines" => length(curves),
            "ReplacedMaximumEdgeAspect" => replaced_worst,
            "CavityMaximumEdgeAspect" => best_aspect,
            "ReplacedMinimumScaledJacobian" => replaced_scaled,
            "CavityMinimumScaledJacobian" => best_scaled,
            "ReplacedMaximumJacobianCondition" => replaced_condition,
            "CavityMaximumJacobianCondition" => best_condition))
    end
    function apply!(elements, deleted_set, remapped_map)
        for k in keys(remapped_map)
            k in deleted_set && delete!(remapped_map, k)
        end
        for (k, element) in remapped_map
            elements[k] = element
        end
        local removed_indices = sort!(collect(deleted_set))
        deleteat!(elements, removed_indices)
        return removed_indices
    end
    removed = apply!(tetrahedra, deleted, remapped)
    removed_triangles = apply!(triangles, deleted_triangles, remapped_triangles)
    removed_lines = apply!(lines, deleted_lines, remapped_lines)
    return (cells_removed=removed, cells_remapped=remapped,
            triangles_removed=removed_triangles, triangles_remapped=remapped_triangles,
            lines_removed=removed_lines, lines_remapped=remapped_lines,
            records=records, shortest=shortest)
end

# The (dimension, entity) supports of every vertex from the line and triangle
# elements it belongs to (Set per vertex; empty for a volume vertex).
function vertex_supports(point_count, triangles, triangle_entities, lines, line_entities)
    supports = [Set{Tuple{Int, Int}}() for _ in 1:point_count]
    for (triangle, entity) in zip(triangles, triangle_entities), i in triangle
        push!(supports[i], (2, Int(entity)))
    end
    for (line, entity) in zip(lines, line_entities), i in line
        push!(supports[i], (1, Int(entity)))
    end
    return supports
end

# Optimize the required region of a linear tetrahedral mesh in place and gate it.
# The required set is the metric stage's rule evaluated on the FINAL vertex
# positions: the moves can carry vertices across the corner-ball radius or the
# layer reach, so after the first pass the set is recomputed; a changed set gets
# one more pass and is recomputed again, and the gates are judged on that final
# set, which is the set the metric stage will list (the census count must equal
# the recipe count). With edge_layer_maximum_aspect > 0 (the layer quality rule,
# supervisor decision 32, calibration only) the layer cells (a vertex within the
# reach) are gated by positive orientation above EDGE_LAYER_ORIENTATION_FLOOR and
# by longest edge / shortest height <= the bound (repaired by bounded descent on
# that aspect, target 0.95 x the bound) and the scaled-Jacobian gate judges the
# other required cells; with 0 every required cell is gated by the scaled
# Jacobian. The scaled-Jacobian-gated cells are also bounded by the Jacobian
# condition number maximum_jacobian_condition (the restorer's/manifest gate
# MaximumJacobianCondition; supervisor decision 34), reported and not repaired.
# Before the moves the sub-size edges of the required region are collapsed
# (collapse_short_edges!; `triangle_entities`/`lines`/`line_entities` carry the
# surface supports and `fixed` the CAD points). Fails closed when a corner exceeds
# maximum_corner_aspect, a scaled-Jacobian-gated cell stays below
# minimum_scaled_jacobian or above maximum_jacobian_condition, or a layer cell
# fails the rule. Returns (census record, indices of the moved vertices, the
# collapse NamedTuple of collapse_short_edges!).
function optimize_required_region!(points, tetrahedra, triangles, corners, radius, lc_fine,
                                   spans, edge_size, growth_ratio, layer_thickness, row_zigzag,
                                   maximum_corner_aspect, minimum_scaled_jacobian,
                                   maximum_jacobian_condition, displacement_ratio, tolerance;
                                   edge_layer_maximum_aspect=0.0,
                                   corner_grading=CornerGrading(0.0, growth_ratio, lc_fine, radius),
                                   triangle_entities=zeros(Int, length(triangles)),
                                   lines=NTuple{2, Int}[], line_entities=Int[],
                                   fixed=falses(size(points, 2)),
                                   corner_kinds=fill(:legacy, length(corners)),
                                   corner_sides=fill(nothing, length(corners)),
                                   corner_shape_gate=0.0,
                                   cell_entities=zeros(Int, length(tetrahedra)))
    original = copy(points)
    original_cell_count = length(tetrahedra)
    length(cell_entities) == original_cell_count || error("Every seed cell needs its volume entity")
    length(corner_kinds) == length(corners) == length(corner_sides) ||
        error("Every semantic corner needs a kind and its sides")
    invariant_corners = count(kind -> kind === :invariant, corner_kinds)
    invariant_corners == 0 || (isfinite(corner_shape_gate) && corner_shape_gate > 1.0) ||
        error("$(invariant_corners) invariant (non-perpendicular) semantic corners need " *
              "--corner-shape-gate (CornerShapeGate = min(E_pop, 5.0) above 1; $(INVARIANT_CORNER_RULE))")
    reach = isempty(spans) ? 0.0 : layer_thickness * (1.0 + row_zigzag) + edge_size
    layer_rule = edge_layer_maximum_aspect > 0.0
    layer_rule && isempty(spans) &&
        error("The edge layer quality rule needs a seeded edge layer")
    isfinite(maximum_jacobian_condition) && maximum_jacobian_condition > 1.0 ||
        error("The required-region Jacobian condition bound must be a finite number above 1")
    required, span_distance, in_layer = required_region_cells(points, tetrahedra, corners,
                                                              radius, spans, reach)
    aspect_before_collapse = !layer_rule ? 0.0 :
        maximum((tetrahedron_edge_aspect([points[:, i] for i in tetrahedra[k]])
                 for k in findall(in_layer)); init=0.0)
    # The collapse size is the smallest prescribed size (EdgeSize inside the layer,
    # CornerSize inside a graded ball); its candidates are the vertices of those
    # regions. Without a layer or a corner grading nothing is below lc_fine by
    # construction and nothing is collapsed.
    corner_distance = [minimum((norm(points[:, i] .- collect(corner)) for corner in corners);
                               init=Inf) for i in axes(points, 2)]
    collapse_sizes = filter(>(0.0), [isempty(spans) ? 0.0 : edge_size, corner_grading.corner_size])
    collapse_size = isempty(collapse_sizes) ? 0.0 : minimum(collapse_sizes)
    collapse_threshold = SEED_COLLAPSE_SIZE_FRACTION * collapse_size
    candidate = [(!isempty(spans) && span_distance[i] <= reach) ||
                 (corner_grading.corner_size > 0.0 && corner_distance[i] <= radius)
                 for i in axes(points, 2)]
    collapse = collapse_short_edges!(points, tetrahedra, triangles, lines, candidate, fixed,
                                     vertex_supports(size(points, 2), triangles, triangle_entities,
                                                     lines, line_entities),
                                     collapse_threshold;
                                     scaled_jacobian_gate=layer_rule ? 0.0 : minimum_scaled_jacobian)
    if !isempty(collapse.records)
        required, span_distance, in_layer = required_region_cells(points, tetrahedra, corners,
                                                                  radius, spans, reach)
    end
    incident = [Int[] for _ in axes(points, 2)]
    for (k, cell) in enumerate(tetrahedra), i in cell
        push!(incident[i], k)
    end
    # F5-B bookkeeping: the original index (Gmsh tag slot) of every current cell (0 for a
    # cell created by a reconnection), its volume entity, the cells a reconnection
    # replaced (original indices) and added (current index => (entity, cell)); a cell the
    # collapse remapped keeps its fresh Gmsh element and is never reconnected.
    origin = setdiff(1:original_cell_count, collapse.cells_removed)
    length(origin) == length(tetrahedra) || error("Seed collapse bookkeeping lost a cell")
    entity_of = [cell_entities[o] for o in origin]
    remapped_by_collapse = Set(keys(collapse.cells_remapped))
    reconnection_replaced = Int[]
    reconnection_added = Dict{Int, Tuple{Int, NTuple{4, Int}}}()
    surface_faces = Set(Tuple(sort!(collect(triangle))) for triangle in triangles)
    surface_edges = Set(minmax(triangle[i], triangle[mod1(i + 1, 3)]) for triangle in triangles for i in 1:3)
    dead = Set{Int}()
    bases = surface_movement_bases(points, triangles, lines)
    scaled_target = 2.0 * minimum_scaled_jacobian
    corner_target = 0.95 * maximum_corner_aspect
    original_scaled = [tetrahedron_scaled_jacobian([points[:, i] for i in cell])
                       for cell in tetrahedra]
    all(>(0.0), original_scaled) || error("Seed contains a nonpositive tetrahedron")
    floors = min.(original_scaled, scaled_target)
    # Under the layer rule a layer cell's guards are the rule's: orientation above
    # the roundoff floor and max(original edge aspect, target) as its ceiling
    # through every pass (corner, scaled-Jacobian and aspect repairs alike); its
    # vertex-0 scaled Jacobian is not a quality measure there (a layer cell's value
    # depends on which vertex is first, and a floor of min(original, target) was
    # measured to block every aspect repair).
    ceilings = fill(Inf, length(tetrahedra))
    # Every touched cell also keeps max(original, target) Jacobian condition (the
    # layer cells under the rule excepted: the condition gate judges the cells
    # outside the layer there) and, inside the collapse region, the edges of a
    # moved vertex stay at or above min(original, the collapse threshold).
    condition_target = 0.95 * maximum_jacobian_condition
    condition_ceilings = [max(tetrahedron_aspect([points[:, i] for i in cell]), condition_target)
                          for cell in tetrahedra]
    edge_floors = [candidate[i] ? collapse_threshold : 0.0 for i in axes(points, 2)]
    if layer_rule
        for k in findall(in_layer)
            floors[k] = EDGE_LAYER_ORIENTATION_FLOOR
            ceilings[k] = max(tetrahedron_edge_aspect([points[:, i] for i in tetrahedra[k]]),
                              0.95 * edge_layer_maximum_aspect)
            condition_ceilings[k] = Inf
        end
    end
    # Local prescribed size: lc_fine, EdgeSize-graded inside the layer reach,
    # CornerSize-graded inside a graded corner ball (the bounds are set once, from
    # the original positions).
    bounds = [displacement_ratio * min(lc_fine,
              isempty(spans) ? lc_fine : edge_size + (growth_ratio - 1.0) * span_distance[i],
              corner_distance[i] <= radius ? corner_ball_size(corner_grading, corner_distance[i]) :
                                             lc_fine)
              for i in axes(points, 2)]
    # Per corner: the measure of its kind (legacy: tetrahedron_aspect, target 0.95 x
    # MaximumCornerAspect, verdict MaximumCornerAspect - bitwise the production path;
    # invariant: kappa_reg, goal INVARIANT_CORNER_TARGET, verdict corner_shape_gate), plus
    # the frame extremes / mean ratio and the bridging-sliver candidates (decision 365:
    # the reconnection trigger = a candidate still above the gate) as information.
    corner_before = Float64[]; corner_after = Float64[]; corner_moves = Int[]
    corner_measures = Dict{String, Any}[]
    corner_incident(center) = [k for (k, cell) in enumerate(tetrahedra)
                               if !(k in dead) && any(norm(points[:, i] .- center) <= tolerance for i in cell)]
    function corner_information(targets)
        cells = [[points[:, i] for i in tetrahedra[k]] for k in targets]
        frames = [tetrahedron_frame_conditions(xyz) for xyz in cells]
        return Dict{String, Any}(
            "KappaRegMax" => maximum(tetrahedron_regular_condition(xyz) for xyz in cells),
            "KappaV0Max" => maximum(tetrahedron_aspect(xyz) for xyz in cells),
            "KappaMinMax" => maximum(first(f) for f in frames),
            "KappaMaxMax" => maximum(last(f) for f in frames),
            "EtaMin" => minimum(tetrahedron_mean_ratio(xyz) for xyz in cells),
            "Cells" => length(targets))
    end
    function bridging_slivers(center, sides, targets)
        sides === nothing && return Int[]
        masks = sidewall_vertex_masks(points, triangles, center, sides, radius)
        return bridging_sliver_cells(points, tetrahedra, targets, masks, center, sides)
    end
    above_gate(cells, gate) = [k for k in cells
                               if tetrahedron_regular_condition([points[:, i] for i in tetrahedra[k]]) > gate]
    for (corner, kind, sides) in zip(corners, corner_kinds, corner_sides)
        center = collect(corner)
        targets = corner_incident(center)
        isempty(targets) && error("Semantic corner is absent from the seed volume mesh")
        invariant = kind === :invariant
        measure_of = invariant ? tetrahedron_regular_condition : tetrahedron_aspect
        objective = invariant ? :regular_condition : :aspect
        target = invariant ? INVARIANT_CORNER_TARGET : corner_target
        gate = invariant ? corner_shape_gate : maximum_corner_aspect
        before = maximum(measure_of([points[:, i] for i in tetrahedra[k]]) for k in targets)
        information_before = corner_information(targets)
        bridging_before = invariant ? bridging_slivers(center, sides, targets) : Int[]
        push!(corner_before, before)
        if before <= target
            after, moves = before, 0
        else
            after, moves = optimize_seed_cells!(points, original, tetrahedra, incident, bases,
                                                floors, bounds, targets, objective, target;
                                                ceilings=ceilings,
                                                condition_ceilings=condition_ceilings,
                                                edge_floors=edge_floors)
        end
        bridging_after = invariant ? bridging_slivers(center, sides, targets) : Int[]
        # Decision 365: the reconnection trigger = a bridging-sliver candidate still above
        # the gate after the kappa_reg descent (none at a legacy corner). F5-B: each trigger
        # cell takes its best 2-3 flip or edge removal, then the descent runs again, for up
        # to CORNER_RECONNECTION_ROUNDS rounds while triggers remain and a flip was found.
        trigger = invariant ? above_gate(bridging_after, gate) : Int[]
        reconnections = Dict{String, Any}[]
        after_descent = after
        for round in 1:CORNER_RECONNECTION_ROUNDS
            isempty(trigger) && break
            flips = 0
            # The trigger cells as they were when the round started (decision 392 MINOR-7): an
            # earlier reconnection of this round may reuse a later trigger's slot for a NEW
            # cell (the replaced slots take the new cells); such a slot is skipped so every
            # Element / Before record names a trigger cell - the slot's cell, if still a
            # candidate above the gate after the descent, triggers again next round.
            trigger_cells = [tetrahedra[k] for k in trigger]
            for (slot, k) in enumerate(trigger)
                k in dead && continue
                tetrahedra[k] == trigger_cells[slot] || continue
                origin[k] in remapped_by_collapse && continue
                before_flip = tetrahedron_regular_condition([points[:, i] for i in tetrahedra[k]])
                flip = best_two_three_flip(points, tetrahedra, incident, surface_faces,
                                           p -> entity_of[p], k, floors[k], condition_ceilings[k])
                removal = best_edge_removal(points, tetrahedra, incident, surface_edges,
                                            p -> entity_of[p], k, floors[k], condition_ceilings[k])
                flip !== nothing && origin[flip.neighbour] in remapped_by_collapse && (flip = nothing)
                removal !== nothing && any(origin[j] in remapped_by_collapse for j in removal.ring) &&
                    (removal = nothing)
                flip === nothing && removal === nothing && continue
                use_removal = removal !== nothing &&
                              (flip === nothing || removal.quality <= flip.quality)
                quality = use_removal ? removal.quality : flip.quality
                quality < before_flip || continue
                replaced = use_removal ? removal.ring : [k, flip.neighbour]
                new_cells = use_removal ? removal.cells : flip.cells
                entity = entity_of[k]
                for p in replaced
                    for i in tetrahedra[p]
                        filter!(!=(p), incident[i])
                    end
                    origin[p] > 0 && push!(reconnection_replaced, origin[p])
                    origin[p] = 0
                    delete!(reconnection_added, p)
                end
                # The new cells take the replaced slots; extra slots die, extra cells append.
                slots = Int[]
                for (index, new_cell) in enumerate(new_cells)
                    if index <= length(replaced)
                        p = replaced[index]
                        tetrahedra[p] = new_cell
                        scaled_new = tetrahedron_scaled_jacobian([points[:, i] for i in new_cell])
                        floors[p] = min(scaled_new, scaled_target); ceilings[p] = Inf
                        condition_ceilings[p] = max(tetrahedron_aspect([points[:, i] for i in new_cell]),
                                                    condition_target)
                    else
                        push!(tetrahedra, new_cell)
                        p = length(tetrahedra)
                        push!(origin, 0); push!(entity_of, entity)
                        scaled_new = tetrahedron_scaled_jacobian([points[:, i] for i in new_cell])
                        push!(floors, min(scaled_new, scaled_target)); push!(ceilings, Inf)
                        push!(condition_ceilings,
                              max(tetrahedron_aspect([points[:, i] for i in new_cell]), condition_target))
                    end
                    push!(slots, p)
                    reconnection_added[p] = (entity, new_cell)
                    for i in new_cell
                        push!(incident[i], p)
                    end
                end
                for p in replaced[(length(new_cells) + 1):end]
                    push!(dead, p)
                end
                flips += 1
                push!(reconnections, Dict{String, Any}(
                    "Round" => round, "Kind" => use_removal ? "edge-removal" : "2-3",
                    "Element" => use_removal ? collect(removal.edge) : collect(flip.face),
                    "Before" => before_flip, "After" => quality,
                    "ReplacedCells" => length(replaced), "AddedCells" => length(new_cells)))
            end
            flips == 0 && break
            targets = corner_incident(center)
            after, descent_moves = optimize_seed_cells!(points, original, tetrahedra, incident,
                                                        bases, floors, bounds, targets, objective,
                                                        target; ceilings=ceilings,
                                                        condition_ceilings=condition_ceilings,
                                                        edge_floors=edge_floors)
            moves += descent_moves
            bridging_after = bridging_slivers(center, sides, targets)
            trigger = above_gate(bridging_after, gate)
        end
        push!(corner_after, after); push!(corner_moves, moves)
        push!(corner_measures, Dict{String, Any}(
            "Point" => center, "Kind" => invariant ? "Invariant" : "Legacy",
            "Measure" => invariant ? "RegularCondition" : "VertexFrameCondition",
            "Before" => before, "After" => after, "Target" => target, "Gate" => gate,
            "Moves" => moves,
            "Sides" => sides === nothing ? nothing :
                       Dict{String, Any}("Walls" => [collect(w) for w in sides.walls],
                                         "Metal" => collect(sides.metal)),
            "BridgingSlivers" => Dict{String, Any}("Before" => length(bridging_before),
                                                   "After" => length(bridging_after),
                                                   "AboveGateAfter" => length(trigger)),
            "AfterDescent" => after_descent,
            "Reconnections" => reconnections,
            "InformationBefore" => information_before,
            "Information" => corner_information(targets),
            "Passed" => after <= gate))
    end
    if !isempty(reconnection_added)
        # The reconnections changed the cell list: the dead slots (an edge removal replaces
        # n cells by 2 (n - 2)) are compacted out of every per-cell array, the incidence
        # rebuilt and the required set recomputed on the new list.
        if !isempty(dead)
            keep = [p for p in eachindex(tetrahedra) if !(p in dead)]
            remap = Dict(p => q for (q, p) in enumerate(keep))
            tetrahedra_kept = tetrahedra[keep]
            empty!(tetrahedra); append!(tetrahedra, tetrahedra_kept)
            floors = floors[keep]; ceilings = ceilings[keep]
            condition_ceilings = condition_ceilings[keep]
            origin = origin[keep]; entity_of = entity_of[keep]
            reconnection_added = Dict(remap[p] => value for (p, value) in reconnection_added if haskey(remap, p))
            empty!(dead)
            for list in incident
                empty!(list)
            end
            for (k, cell) in enumerate(tetrahedra), i in cell
                push!(incident[i], k)
            end
        end
        required, span_distance, in_layer = required_region_cells(points, tetrahedra, corners,
                                                                  radius, spans, reach)
    end
    scaled_of(k) = tetrahedron_scaled_jacobian([points[:, i] for i in tetrahedra[k]])
    edge_aspect_of(k) = tetrahedron_edge_aspect([points[:, i] for i in tetrahedra[k]])
    condition_of(k) = tetrahedron_aspect([points[:, i] for i in tetrahedra[k]])
    required_cells = findall(required)
    required_before_moves = length(required_cells)
    # Under the layer rule the scaled-Jacobian gate judges the required cells
    # outside the layer; the layer cells are judged by the rule.
    scaled_gated(cells, layer) = layer_rule ? [k for k in cells if !layer[k]] : cells
    layer_cells = layer_rule ? findall(in_layer) : Int[]
    gated_cells = scaled_gated(required_cells, in_layer)
    required_before = minimum(scaled_of(k) for k in gated_cells; init=Inf)
    below_before = count(k -> scaled_of(k) < scaled_target, gated_cells)
    condition_before = maximum(condition_of(k) for k in gated_cells; init=0.0)
    layer_aspect_target = 0.95 * edge_layer_maximum_aspect
    layer_aspect_before = maximum(edge_aspect_of(k) for k in layer_cells; init=0.0)
    layer_above_before = count(k -> edge_aspect_of(k) > layer_aspect_target, layer_cells)
    layer_minimum_before = minimum(scaled_of(k) for k in layer_cells; init=Inf)
    # Components of below-target scaled-Jacobian-gated cells (and, under the rule,
    # of layer cells above the aspect target) sharing a vertex are repaired
    # together; then the required set is recomputed on the moved positions and a
    # changed set gets one more pass (see above).
    repair_moves = 0; components = 0; recomputations = 0
    layer_repair_moves = 0; layer_components = 0
    for pass in 1:2
        for component in below_target_components(gated_cells, tetrahedra, incident,
                                                 scaled_of, scaled_target)
            _, moves = optimize_seed_cells!(points, original, tetrahedra, incident, bases,
                                            floors, bounds, component, :scaled, scaled_target;
                                            ceilings=ceilings, condition_ceilings=condition_ceilings,
                                            edge_floors=edge_floors)
            repair_moves += moves; components += 1
        end
        for component in below_target_components(layer_cells, tetrahedra, incident,
                                                 k -> -edge_aspect_of(k), -layer_aspect_target)
            _, moves = optimize_seed_cells!(points, original, tetrahedra, incident, bases,
                                            floors, bounds, component, :edge_aspect,
                                            layer_aspect_target; ceilings=ceilings,
                                            condition_ceilings=condition_ceilings,
                                            edge_floors=edge_floors)
            layer_repair_moves += moves; layer_components += 1
        end
        recomputed, _, recomputed_layer = required_region_cells(points, tetrahedra, corners,
                                                                radius, spans, reach)
        recomputed_cells = findall(recomputed)
        recomputations += 1
        changed = recomputed_cells != required_cells ||
                  (layer_rule && findall(recomputed_layer) != layer_cells)
        required_cells = recomputed_cells
        in_layer = recomputed_layer
        layer_cells = layer_rule ? findall(in_layer) : Int[]
        gated_cells = scaled_gated(required_cells, in_layer)
        changed || break
    end
    required_after = minimum(scaled_of(k) for k in gated_cells; init=Inf)
    below_after = count(k -> scaled_of(k) < scaled_target, gated_cells)
    below_gate = count(k -> scaled_of(k) < minimum_scaled_jacobian, gated_cells)
    condition_after = maximum(condition_of(k) for k in gated_cells; init=0.0)
    above_condition = count(k -> condition_of(k) > maximum_jacobian_condition, gated_cells)
    layer_aspect_after = maximum(edge_aspect_of(k) for k in layer_cells; init=0.0)
    layer_above_after = count(k -> edge_aspect_of(k) > layer_aspect_target, layer_cells)
    layer_above_bound = count(k -> edge_aspect_of(k) > edge_layer_maximum_aspect, layer_cells)
    layer_minimum_after = minimum(scaled_of(k) for k in layer_cells; init=Inf)
    layer_flat = count(k -> !(scaled_of(k) > EDGE_LAYER_ORIENTATION_FLOOR), layer_cells)
    displacement = [norm(points[:, i] .- original[:, i]) for i in axes(points, 2)]
    moved = findall(>(0.0), displacement)
    all(displacement[i] <= bounds[i] * (1.0 + 1.0e-12) for i in moved) ||
        error("Seed quality optimization exceeded its displacement bound")
    layer_record = !layer_rule ? nothing : Dict{String, Any}(
        "Rule" => "inside the recorded edge layer (a vertex within LayerRequiredReach of a " *
                  "span) a cell passes when its scaled Jacobian exceeds the roundoff floor " *
                  "(positive orientation) and its longest edge over its shortest height is " *
                  "at most MaximumEdgeAspect (repaired by bounded descent on that aspect, " *
                  "target 0.95 x the bound); the scaled-Jacobian gate judges the required " *
                  "cells outside the layer; the layer scaled Jacobian is a diagnostic " *
                  "(supervisor decision 32, calibration only)",
        "MaximumEdgeAspect" => edge_layer_maximum_aspect,
        "EdgeAspectTarget" => layer_aspect_target,
        "ScaledJacobianRoundoffFloor" => EDGE_LAYER_ORIENTATION_FLOOR,
        "LayerCells" => length(layer_cells),
        "MaximumEdgeAspectBeforeCollapse" => aspect_before_collapse,
        "MaximumEdgeAspectBefore" => layer_aspect_before,
        "MaximumEdgeAspectAfter" => layer_aspect_after,
        "CellsAboveTargetBefore" => layer_above_before,
        "CellsAboveTargetAfter" => layer_above_after,
        "CellsAboveBoundAfter" => layer_above_bound,
        "CellsBelowRoundoffFloorAfter" => layer_flat,
        "MinimumScaledJacobianBefore" => layer_minimum_before,
        "MinimumScaledJacobian" => layer_minimum_after,
        "RepairComponents" => layer_components, "RepairMoves" => layer_repair_moves)
    record = Dict{String, Any}(
        "Rule" => "seed cells inside the corner balls (centroid within CornerIsotropyRadius) " *
                  "and the edge layer (a vertex within LayerRequiredReach = LayerThickness x " *
                  "(1 + RowZigzag) + EdgeSize of a span) are MMG required tetrahedra whose " *
                  "quality after adaptation is the seed's; the seed repairs them by bounded " *
                  "coordinate descent (surface vertices in their planes, " *
                  "DisplacementBoundOverNormal x the local prescribed size, touched cells " *
                  "keep min(original, target) scaled Jacobian), recomputes the set on the " *
                  "moved positions (one more pass if it changed) and fails closed on the " *
                  "gates judged on the final set, which the metric stage lists; under the " *
                  "EdgeLayerQualityRule the scaled-Jacobian statistics cover the required " *
                  "cells outside the layer (ScaledJacobianGateCells)",
        "EdgeLayerQualityRule" => layer_record,
        "SubSizeEdgeCollapse" => Dict{String, Any}(
            "Rule" => "before the moves, a required-region seed vertex (within the layer reach " *
                      "or, with corner grading, inside a corner ball) with an incident edge " *
                      "shorter than CollapseThreshold = ThresholdOverSize x CollapseSize (the " *
                      "smallest prescribed size: EdgeSize, the adapter hmin, or CornerSize) is " *
                      "merged onto the neighbour whose cavity " *
                      "(every remapped cell positively oriented above the roundoff floor, every " *
                      "remapped triangle keeping its normal) has the smallest maximum edge " *
                      "aspect, no worse than the cells it replaces in maximum edge aspect and " *
                      "maximum Jacobian condition, better in one of them, and with no cell " *
                      "below ScaledJacobianGate (0: none, the layer quality rule) unless a " *
                      "replaced cell already was; a free interior " *
                      "vertex onto any neighbour, a surface vertex only along its own surface " *
                      "(a surface neighbour in every plane of the vertex carrying every " *
                      "line/triangle support of the vertex); CAD points and vertices in three " *
                      "planes are never collapsed (supervisor decisions 32 and 34)",
            "CollapseSize" => collapse_size > 0.0 ? collapse_size : nothing,
            "ThresholdOverSize" => SEED_COLLAPSE_SIZE_FRACTION,
            "CollapseThreshold" => collapse_size > 0.0 ? collapse_threshold : nothing,
            "PlaneTolerance" => SEED_COLLAPSE_PLANE_TOLERANCE,
            "ScaledJacobianGate" => layer_rule ? 0.0 : minimum_scaled_jacobian,
            "CandidateVertices" => count(candidate),
            "CollapsedVertices" => length(collapse.records),
            "CollapsedInteriorVertices" => count(row -> !row["Surface"], collapse.records),
            "CollapsedSurfaceVertices" => count(row -> row["Surface"], collapse.records),
            "CollapsedCells" => length(collapse.cells_removed),
            "RemappedCells" => length(collapse.cells_remapped),
            "CollapsedTriangles" => length(collapse.triangles_removed),
            "RemappedTriangles" => length(collapse.triangles_remapped),
            "CollapsedLines" => length(collapse.lines_removed),
            "RemappedLines" => length(collapse.lines_remapped),
            "ShortestCollapsedEdge" => isfinite(collapse.shortest) ? collapse.shortest : nothing,
            "Collapses" => [Dict{String, Any}(name => row[name] for name in
                                              ("Position", "TargetPosition", "Surface", "EdgeLength",
                                               "Cells", "Triangles", "Lines",
                                               "ReplacedMaximumEdgeAspect", "CavityMaximumEdgeAspect",
                                               "ReplacedMinimumScaledJacobian",
                                               "CavityMinimumScaledJacobian",
                                               "ReplacedMaximumJacobianCondition",
                                               "CavityMaximumJacobianCondition"))
                            for row in collapse.records]),
        "ScaledJacobianGateCells" => length(gated_cells),
        "MaximumCornerAspect" => maximum_corner_aspect, "CornerAspectTarget" => corner_target,
        "InvariantCornerRule" => INVARIANT_CORNER_RULE,
        "InvariantCornerTarget" => INVARIANT_CORNER_TARGET,
        "CornerReconnectionRule" => CORNER_RECONNECTION_RULE,
        "CornerReconnectionRounds" => CORNER_RECONNECTION_ROUNDS,
        "CornerReconnections" => sum(length(row["Reconnections"]) for row in corner_measures; init=0),
        "ReconnectionReplacedCells" => length(reconnection_replaced),
        "ReconnectionAddedCells" => length(reconnection_added),
        "CornerShapeGate" => corner_shape_gate > 0.0 ? corner_shape_gate : nothing,
        "InvariantCorners" => invariant_corners,
        "CornerMeasures" => corner_measures,
        "MinimumScaledJacobian" => minimum_scaled_jacobian, "ScaledJacobianTarget" => scaled_target,
        "MaximumJacobianCondition" => maximum_jacobian_condition,
        "JacobianConditionRule" => "the Jacobian condition number (largest over smallest " *
                                   "singular value of the edge-vector Jacobian, the audit " *
                                   "producer's MaximumJacobianCondition) of every " *
                                   "scaled-Jacobian-gated required cell is at most " *
                                   "MaximumJacobianCondition; not an objective: every " *
                                   "touched cell keeps max(original, JacobianConditionTarget) " *
                                   "through the moves and a moved vertex in the collapse " *
                                   "region keeps its edges at or above min(original, " *
                                   "CollapseThreshold)",
        "JacobianConditionTarget" => condition_target,
        "RequiredMaximumJacobianConditionBefore" => condition_before,
        "RequiredMaximumJacobianConditionAfter" => condition_after,
        "RequiredCellsAboveConditionAfter" => above_condition,
        "DisplacementBoundOverNormal" => displacement_ratio,
        "LocalSizeRule" => "min(NormalSize, EdgeSize + (GrowthRatio - 1) x span distance " *
                           "inside the layer, the corner law inside a corner ball)",
        "CornerSize" => corner_grading.corner_size,
        "LayerRequiredReach" => isempty(spans) ? nothing : reach,
        "RequiredTetrahedra" => length(required_cells),
        "RequiredTetrahedraBeforeMoves" => required_before_moves,
        "RequiredSetRecomputations" => recomputations,
        "CornerAspectsBefore" => corner_before, "CornerAspectsAfter" => corner_after,
        "CornerMoves" => corner_moves,
        "RequiredMinimumScaledJacobianBefore" => required_before,
        "RequiredMinimumScaledJacobianAfter" => required_after,
        "RequiredCellsBelowTargetBefore" => below_before,
        "RequiredCellsBelowTargetAfter" => below_after,
        "RequiredCellsBelowGateAfter" => below_gate,
        "RepairComponents" => components, "RepairMoves" => repair_moves,
        "MovedVertices" => length(moved),
        "MaximumDisplacement" => isempty(moved) ? 0.0 : maximum(displacement[moved]),
        "MaximumDisplacementOverBound" =>
            isempty(moved) ? 0.0 : maximum(displacement[i] / bounds[i] for i in moved))
    for row in corner_measures
        println("Seed corner $(row["Kind"]) $(row["Point"]): $(row["Measure"]) $(row["Before"]) -> " *
                "$(row["After"]) (target $(row["Target"]), gate $(row["Gate"]), moves $(row["Moves"])), " *
                "bridging slivers $(row["BridgingSlivers"]["Before"]) -> " *
                "$(row["BridgingSlivers"]["After"]) (above the gate $(row["BridgingSlivers"]["AboveGateAfter"])), " *
                "reconnections $(length(row["Reconnections"])) (after the descent $(row["AfterDescent"])), " *
                "kappa_reg $(row["Information"]["KappaRegMax"]) " *
                "kappa_v0 $(row["Information"]["KappaV0Max"]) kappa_min $(row["Information"]["KappaMinMax"]) " *
                "kappa_max $(row["Information"]["KappaMaxMax"]) eta $(row["Information"]["EtaMin"])")
    end
    println("Seed required-region optimization: corners $(corner_before) -> $(corner_after), " *
            "required cells $(required_before_moves) -> $(length(required_cells)) " *
            "($(recomputations) recomputations) min scaled Jacobian $(required_before) -> " *
            "$(required_after) (below target $(below_before) -> $(below_after)), max " *
            "Jacobian condition $(condition_before) -> $(condition_after) (above the bound " *
            "$(maximum_jacobian_condition): $(above_condition)), moved vertices $(length(moved))")
    println("Seed sub-size edge collapse: collapse size $(collapse_size), threshold " *
            "$(collapse_threshold), " *
            "$(count(candidate)) candidate vertices, collapsed $(length(collapse.records)) " *
            "($(count(row -> row["Surface"], collapse.records)) surface, " *
            "$(length(collapse.cells_removed)) cells removed, " *
            "$(length(collapse.cells_remapped)) remapped, " *
            "$(length(collapse.triangles_removed)) triangles removed, " *
            "$(length(collapse.triangles_remapped)) remapped, " *
            "$(length(collapse.lines_removed)) lines removed, " *
            "$(length(collapse.lines_remapped)) remapped), shortest collapsed edge " *
            "$(collapse.shortest)")
    for row in collapse.records
        println("  collapsed $(row["Surface"] ? "surface" : "interior") vertex at " *
                "$(row["Position"]) onto $(row["TargetPosition"]) (edge $(row["EdgeLength"]), " *
                "$(row["Cells"]) cells): max edge aspect $(row["ReplacedMaximumEdgeAspect"]) -> " *
                "$(row["CavityMaximumEdgeAspect"]), min scaled Jacobian " *
                "$(row["ReplacedMinimumScaledJacobian"]) -> $(row["CavityMinimumScaledJacobian"]), " *
                "max condition $(row["ReplacedMaximumJacobianCondition"]) -> " *
                "$(row["CavityMaximumJacobianCondition"])")
    end
    # Failing cells are located for the report: centroid, nearest-corner distance,
    # span distance, edge lengths and vertices (position, surface flag, moved flag).
    surface = falses(size(points, 2))
    for triangle in triangles, i in triangle
        surface[i] = true
    end
    moved_set = Set(moved)
    function describe(label, k, value)
        cell = tetrahedra[k]
        centroid = sum(points[:, i] for i in cell) ./ 4
        lengths = [norm(points[:, cell[i]] .- points[:, cell[j]]) for i in 1:4 for j in (i + 1):4]
        println("  $label cell $k: $(value), centroid $(centroid), " *
                "corner distance $(minimum(norm(centroid .- collect(c)) for c in corners)), " *
                "span distance $(minimum(span_distance[i] for i in cell)), " *
                "surface vertices $(count(surface[i] for i in cell)), edges $(sort(lengths)), " *
                "vertices $([(points[:, i], surface[i], i in moved_set) for i in cell])")
    end
    failures = String[]
    legacy_failed = [row["After"] for row in corner_measures
                     if row["Kind"] == "Legacy" && row["After"] > row["Gate"]]
    invariant_failed = [row for row in corner_measures
                        if row["Kind"] == "Invariant" && row["After"] > row["Gate"]]
    for row in corner_measures
        row["Passed"] && continue
        center = row["Point"]
        incident = corner_incident(center)
        measure = row["Kind"] == "Invariant" ? tetrahedron_regular_condition : tetrahedron_aspect
        measure_of(k) = measure([points[:, i] for i in tetrahedra[k]])
        for k in sort(incident; by=measure_of, rev=true)[1:min(6, end)]
            describe("corner $(center) $(row["Measure"])", k, "$(row["Measure"]) $(measure_of(k))")
        end
        if row["BridgingSlivers"]["After"] > 0
            sides = (walls=row["Sides"]["Walls"], metal=row["Sides"]["Metal"])
            masks = sidewall_vertex_masks(points, triangles, center, sides, radius)
            for k in bridging_sliver_cells(points, tetrahedra, incident, masks, center, sides)
                describe("corner $(center) BridgingSliver", k,
                         "kappa_reg $(tetrahedron_regular_condition([points[:, i] for i in tetrahedra[k]])) " *
                         "walls $([get(masks, i, 0x00) for i in tetrahedra[k]])")
            end
        end
    end
    isempty(legacy_failed) ||
        push!(failures, "Seed semantic-corner aspect exceeds the gate after optimization: " *
                        "$(maximum(legacy_failed)) > $(maximum_corner_aspect)")
    isempty(invariant_failed) ||
        push!(failures, "Seed invariant semantic-corner kappa_reg exceeds CornerShapeGate after " *
                        "optimization: $(maximum(row["After"] for row in invariant_failed)) > " *
                        "$(corner_shape_gate) at $([row["Point"] for row in invariant_failed]) " *
                        "(bridging-sliver cells above the gate: " *
                        "$([row["BridgingSlivers"]["AboveGateAfter"] for row in invariant_failed]))")
    if below_gate > 0
        failing = sort([k for k in gated_cells if scaled_of(k) < minimum_scaled_jacobian];
                       by=scaled_of)
        for k in failing[1:min(10, end)]
            describe("below-gate", k, "scaled Jacobian $(scaled_of(k))")
        end
        push!(failures, "Seed required region keeps $(below_gate) cells below the " *
                        "scaled-Jacobian gate $(minimum_scaled_jacobian) after optimization " *
                        "(minimum $(required_after))")
    end
    if above_condition > 0
        failing = sort([k for k in gated_cells if condition_of(k) > maximum_jacobian_condition];
                       by=condition_of, rev=true)
        for k in failing[1:min(10, end)]
            describe("above-condition", k, "Jacobian condition $(condition_of(k))")
        end
        push!(failures, "Seed required region keeps $(above_condition) cells above the " *
                        "Jacobian condition gate $(maximum_jacobian_condition) after " *
                        "optimization (maximum $(condition_after))")
    end
    isempty(failures) || error(join(failures, "; "))
    if layer_rule
        println("Seed edge-layer quality rule: $(length(layer_cells)) layer cells, " *
                "maximum edge aspect $(aspect_before_collapse) -> $(layer_aspect_before) -> " *
                "$(layer_aspect_after) (bound $(edge_layer_maximum_aspect)), minimum scaled " *
                "Jacobian (diagnostic) $(layer_minimum_after)")
        layer_flat == 0 ||
            error("Seed edge layer keeps $(layer_flat) cells flat to roundoff (scaled Jacobian " *
                  "<= $(EDGE_LAYER_ORIENTATION_FLOOR))")
        layer_above_bound == 0 ||
            error("Seed edge layer keeps $(layer_above_bound) cells above the edge aspect bound " *
                  "$(edge_layer_maximum_aspect) after optimization (maximum $(layer_aspect_after))")
    end
    reconnection = (replaced=reconnection_replaced,
                    added=[(entity, cell) for (entity, cell) in values(reconnection_added)])
    return record, moved, collapse, reconnection
end

# The linear cells of one dimension with their Gmsh element tags and entity tags
# (only the vertex nodes define a linear cell; high-order nodes follow).
# Only the linear simplex blocks of the dimension are read (prisms, pyramids and
# quadrangles of a tube mesh are not simplices and stay untouched).
function gmsh_entity_cells(index, dimension, ::Val{width}) where {width}
    cells = Vector{NTuple{width, Int}}(); tags = UInt[]; entities = Int[]
    for (_, entity) in gmsh.model.getEntities(dimension)
        types, element_tags, element_nodes = gmsh.model.mesh.getElements(dimension, entity)
        for (type, block_tags, block) in zip(types, element_tags, element_nodes)
            isempty(block_tags) && continue
            Int(type) == GMSH_LINEAR_ELEMENT_TYPE[dimension] || continue
            nodes_per_element = length(block) ÷ length(block_tags)
            nodes_per_element >= width || continue
            for (n, start) in enumerate(1:nodes_per_element:length(block))
                push!(cells, ntuple(i -> index[block[start + i - 1]], Val(width)))
                push!(tags, block_tags[n]); push!(entities, entity)
            end
        end
    end
    return cells, tags, entities
end

# The volume tetrahedra with their Gmsh element tags and volume entity tags.
gmsh_volume_cells(index) = gmsh_entity_cells(index, 3, Val(4))

# The nodes of the CAD points (never collapsed; the semantic corners among them).
function gmsh_point_nodes(index, point_count)
    fixed = falses(point_count)
    for (_, point) in gmsh.model.getEntities(0)
        tags, _, _ = gmsh.model.mesh.getNodes(0, point)
        for tag in tags
            haskey(index, tag) && (fixed[index[tag]] = true)
        end
    end
    return fixed
end

function optimize_required_seed_region!(corners, radius, lc_fine, spans, edge_size,
                                        growth_ratio, layer_thickness, row_zigzag,
                                        maximum_corner_aspect, minimum_scaled_jacobian,
                                        maximum_jacobian_condition, displacement_ratio, tolerance,
                                        edge_layer_maximum_aspect,
                                        corner_grading=CornerGrading(0.0, growth_ratio, lc_fine,
                                                                     radius);
                                        fixed_node_tags=UInt[],
                                        corner_kinds=fill(:legacy, length(corners)),
                                        corner_sides=fill(nothing, length(corners)),
                                        corner_shape_gate=0.0)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(copy(coordinates), 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    tetrahedra, tags, entities = gmsh_volume_cells(index)
    triangles, triangle_tags, triangle_entities = gmsh_entity_cells(index, 2, Val(3))
    lines, line_tags, line_entities = gmsh_entity_cells(index, 1, Val(2))
    # The CAD point nodes and every node of a non-simplex cell (the prism tubes and
    # their pyramids) stay where they are.
    fixed = gmsh_point_nodes(index, size(points, 2))
    for tag in fixed_node_tags
        fixed[index[tag]] = true
    end
    record, moved, collapse, reconnection = optimize_required_region!(
        points, tetrahedra, triangles, corners, radius, lc_fine, spans, edge_size, growth_ratio,
        layer_thickness, row_zigzag, maximum_corner_aspect, minimum_scaled_jacobian,
        maximum_jacobian_condition, displacement_ratio, tolerance;
        edge_layer_maximum_aspect=edge_layer_maximum_aspect, corner_grading=corner_grading,
        triangle_entities=triangle_entities, lines=lines, line_entities=line_entities,
        fixed=fixed, corner_kinds=corner_kinds, corner_sides=corner_sides,
        corner_shape_gate=corner_shape_gate, cell_entities=Int.(entities))
    for i in moved
        gmsh.model.mesh.setNode(node_tags[i], points[:, i], Float64[])
    end
    apply_seed_cell_collapse!(node_tags, tags, entities, collapse.cells_removed,
                              collapse.cells_remapped)
    apply_seed_cell_reconnection!(node_tags, tags, entities, reconnection.replaced,
                                  reconnection.added)
    apply_seed_cell_collapse!(node_tags, triangle_tags, triangle_entities,
                              collapse.triangles_removed, collapse.triangles_remapped; dimension=2)
    apply_seed_cell_collapse!(node_tags, line_tags, line_entities, collapse.lines_removed,
                              collapse.lines_remapped; dimension=1)
    return record
end

# Gmsh element type of the linear cell of each dimension (line, triangle, tetrahedron).
const GMSH_LINEAR_ELEMENT_TYPE = Dict(1 => 1, 2 => 2, 3 => 4)

# Apply a seed collapse (collapse_short_edges!) of one dimension to the Gmsh
# model: the vanished elements (`removed`, original indices) are deleted and the
# remapped ones (original index => new vertex-index cell) are replaced in their
# entity with fresh element tags; the orphaned vertices are dropped by the writer
# (Mesh.SaveAll is off), so the written points are exactly the used points.
# Returns the number of elements of that dimension the model carries afterwards.
function apply_seed_cell_collapse!(node_tags, tags, entities, removed, remapped; dimension=3)
    element_type = GMSH_LINEAR_ELEMENT_TYPE[dimension]
    width = dimension + 1
    if !isempty(removed) || !isempty(remapped)
        by_entity = Dict{Int, Vector{UInt}}()
        for k in vcat(removed, collect(keys(remapped)))
            push!(get!(by_entity, entities[k], UInt[]), tags[k])
        end
        for (entity, element_tags) in by_entity
            gmsh.model.mesh.removeElements(dimension, entity, element_tags)
        end
        next_tag = gmsh.model.mesh.getMaxElementTag()
        added = Dict{Int, Vector{Int}}()
        for (k, cell) in remapped
            append!(get!(added, entities[k], Int[]), [Int(node_tags[i]) for i in cell])
        end
        for (entity, nodes) in added
            count = length(nodes) ÷ width
            gmsh.model.mesh.addElementsByType(entity, element_type,
                                              collect((next_tag + 1):(next_tag + count)), nodes)
            next_tag += count
        end
    end
    _, element_tags, _ = gmsh.model.mesh.getElements(dimension)
    return sum(length(block) for block in element_tags; init=0)
end

# Apply a corner reconnection (F5-B) to the Gmsh model: the replaced cells (original
# indices, untouched by the collapse) are deleted and the added cells created in their
# volume entity with fresh element tags.
function apply_seed_cell_reconnection!(node_tags, tags, entities, replaced, added)
    isempty(replaced) && isempty(added) && return 0
    by_entity = Dict{Int, Vector{UInt}}()
    for k in replaced
        push!(get!(by_entity, entities[k], UInt[]), tags[k])
    end
    for (entity, element_tags) in by_entity
        gmsh.model.mesh.removeElements(3, entity, element_tags)
    end
    next_tag = gmsh.model.mesh.getMaxElementTag()
    nodes_by_entity = Dict{Int, Vector{Int}}()
    for (entity, cell) in added
        append!(get!(nodes_by_entity, entity, Int[]), [Int(node_tags[i]) for i in cell])
    end
    for (entity, nodes) in nodes_by_entity
        count = length(nodes) ÷ 4
        gmsh.model.mesh.addElementsByType(entity, GMSH_LINEAR_ELEMENT_TYPE[3],
                                          collect((next_tag + 1):(next_tag + count)), nodes)
        next_tag += count
    end
    return length(added)
end

function sorted_median(values)
    sorted = sort(values)
    n = length(sorted)
    return isodd(n) ? sorted[(n + 1) ÷ 2] : 0.5 * (sorted[n ÷ 2] + sorted[n ÷ 2 + 1])
end

# Area of every physical surface label of the linear seed (source-local frame),
# so the etched footprint is asserted from the mesh rather than assumed from
# the producer's inputs. Reported, not gated.
function interface_areas()
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    rows = Dict{String, Any}[]
    for (dim, attribute) in gmsh.model.getPhysicalGroups(2)
        area = 0.0
        triangles = 0
        quadrangles = 0
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, attribute)
            types, element_tags, element_nodes = gmsh.model.mesh.getElements(2, entity)
            for (type, tags, block) in zip(types, element_tags, element_nodes)
                isempty(tags) && continue
                nodes_per_element = length(block) ÷ length(tags)
                name = gmsh.model.mesh.getElementProperties(type)[1]
                if startswith(name, "Triangle")
                    for start in 1:nodes_per_element:length(block)
                        a, b, c = (points[:, index[block[start + i]]] for i in 0:2)
                        area += 0.5 * norm(cross(b .- a, c .- a))
                        triangles += 1
                    end
                elseif startswith(name, "Quadrilateral")
                    for start in 1:nodes_per_element:length(block)
                        a, b, c, d = (points[:, index[block[start + i]]] for i in 0:3)
                        area += 0.5 * norm(cross(b .- a, c .- a)) + 0.5 * norm(cross(c .- a, d .- a))
                        quadrangles += 1
                    end
                else
                    error("Unsupported surface element $name in the interface areas")
                end
            end
        end
        push!(rows, Dict{String, Any}("Attribute" => Int(attribute),
                                      "Name" => gmsh.model.getPhysicalName(dim, attribute),
                                      "Triangles" => triangles, "Quadrangles" => quadrangles,
                                      "Area" => area))
    end
    sort!(rows; by=row -> row["Attribute"])
    return rows
end

# Bins of the along-edge histogram of a longitudinal face's interior nodes; a
# census resolution constant, not a mesh target.
const LONGITUDINAL_FACE_HISTOGRAM_BINS = 10

# Per semantic corner, the un-layered metal-edge length: for every layered curve
# whose CAD endpoint is the corner, the distance from the corner to the span end
# nearest to it (the edge between them carries no layer row), and the maximum.
function unlayered_edge_length_per_corner(corners, layer_curves, radius)
    rows = Dict{String, Any}[]
    for (k, corner) in enumerate(corners)
        center = collect(corner)
        lengths = Float64[]
        for row in layer_curves
            ends = (Vector{Float64}(row["Start"]), Vector{Float64}(row["End"]))
            for (curve_end, span_end) in zip(row["CurveEnds"], ends)
                norm(Vector{Float64}(curve_end) .- center) <= 1.0e-9 * radius || continue
                push!(lengths, norm(span_end .- center))
            end
        end
        push!(rows, Dict{String, Any}(
            "Corner" => k - 1, "Point" => center, "LayeredEdges" => length(lengths),
            "UnlayeredLengths" => lengths,
            "Maximum" => isempty(lengths) ? nothing : maximum(lengths),
            "Minimum" => isempty(lengths) ? nothing : minimum(lengths)))
    end
    return rows
end

# Distance from a semantic corner beyond which the corner size law is lc_tangent,
# so the seed is meant to be unchanged there.
function corner_law_reach(radius, lc_fine, lc_tangent, slope)
    return lc_tangent > lc_fine ? radius + (lc_tangent - lc_fine) / slope : radius
end

# Interior-node and full-height-triangle census of every surface bounded by at
# least two longitudinal feature curves (sidewalls and other ridge-to-ridge
# faces). The interior row of such a face is alignment-fragile: ridge rows that
# are misaligned across the face leave triangles whose vertices all lie on the
# bounding curves and span two longitudinal curves ("full-height" triangles)
# with no interior node. Nodes and triangles farther than `reach` from every
# semantic corner are counted separately, since the seed is meant to be
# unchanged there. Reported, not gated.
function longitudinal_face_census(longitudinal_curves, corners, reach)
    longitudinal = Set{Int32}(longitudinal_curves)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    away(i) = all(norm(points[:, i] .- corner) > reach for corner in corners)
    rows = Dict{String, Any}[]
    for (_, surface) in gmsh.model.getEntities(2)
        curves = Int32[curve for (dim, curve) in
                       gmsh.model.getBoundary([(2, surface)], false, false, false)]
        face_longitudinal = [curve for curve in curves if curve in longitudinal]
        length(face_longitudinal) >= 2 || continue
        curve_nodes = Dict{Int, Set{Int32}}()
        boundary = Set{Int}()
        for curve in curves
            tags, _, _ = gmsh.model.mesh.getNodes(1, curve, true)
            for tag in tags
                node = index[tag]
                push!(boundary, node)
                curve in longitudinal &&
                    push!(get!(curve_nodes, node, Set{Int32}()), curve)
            end
        end
        types, element_tags, element_nodes = gmsh.model.mesh.getElements(2, surface)
        triangles = Vector{NTuple{3, Int}}()
        for (type, tags, block) in zip(types, element_tags, element_nodes)
            isempty(tags) && continue
            startswith(gmsh.model.mesh.getElementProperties(type)[1], "Triangle") || continue
            nodes_per_element = length(block) ÷ length(tags)
            for start in 1:nodes_per_element:length(block)
                push!(triangles, ntuple(i -> index[block[start + i - 1]], 3))
            end
        end
        face_nodes = unique([node for triangle in triangles for node in triangle])
        interior = [node for node in face_nodes if !(node in boundary)]
        interior_away = [node for node in interior if away(node)]
        function full_height(triangle)
            all(node in boundary for node in triangle) || return false
            spanned = union((get(curve_nodes, node, Set{Int32}()) for node in triangle)...)
            length(spanned) >= 2 || return false
            return !any(all(curve in get(curve_nodes, node, Set{Int32}()) for node in triangle)
                        for curve in spanned)
        end
        full = [triangle for triangle in triangles if full_height(triangle)]
        # Along-edge histogram of the interior nodes away from the corners, over
        # the face's extent along its first longitudinal curve's direction.
        lower, upper = gmsh.model.getParametrizationBounds(1, face_longitudinal[1])
        derivative = gmsh.model.getDerivative(1, face_longitudinal[1],
                                              [0.5 * (lower[1] + upper[1])])
        direction = derivative[1:3] ./ norm(derivative[1:3])
        along = [dot(points[:, node], direction) for node in face_nodes]
        start, stop = extrema(along)
        histogram = zeros(Int, LONGITUDINAL_FACE_HISTOGRAM_BINS)
        for node in interior_away
            fraction = (dot(points[:, node], direction) - start) / max(stop - start, eps())
            histogram[clamp(floor(Int, fraction * LONGITUDINAL_FACE_HISTOGRAM_BINS) + 1, 1,
                            LONGITUDINAL_FACE_HISTOGRAM_BINS)] += 1
        end
        push!(rows, Dict{String, Any}(
            "Surface" => Int(surface),
            "PhysicalGroups" => sort!(Int.(gmsh.model.getPhysicalGroupsForEntity(2, surface))),
            "LongitudinalCurves" => length(face_longitudinal),
            "Triangles" => length(triangles),
            "InteriorNodes" => length(interior),
            "InteriorNodesAwayFromCorners" => length(interior_away),
            "InteriorNodeHistogramAlongEdge" => histogram,
            "FullHeightTriangles" => length(full),
            "FullHeightTrianglesAwayFromCorners" =>
                count(triangle -> all(away(node) for node in triangle), full)))
    end
    sort!(rows; by=row -> row["Surface"])
    return rows
end

function extended_interval(edge, radius)
    first, second = edge.interval
    extension = 2radius
    union = get(EDGE_CHAIN_UNIONS, edge_chain_key(edge), nothing)
    if union !== nothing
        # Decision 47: a chained row is extended where it carries the chain's outer end,
        # and only when the chain's union reaches Radius from its own midpoint.
        u0, u1 = union
        tolerance = 1.0e-10radius
        if (u1 - u0) / 2 >= radius - tolerance
            abs(first - u0) <= tolerance && (first -= extension)
            abs(second - u1) <= tolerance && (second += extension)
        end
    elseif edge.vertex_arm
        if abs(first) <= 1.0e-10radius
            second += extension
        elseif abs(second) <= 1.0e-10radius
            first -= extension
        end
    else
        tolerance = 1.0e-10radius
        first <= -radius + tolerance && (first -= extension)
        second >= radius - tolerance && (second += extension)
    end
    return first, second
end

function strip_points(edge, radius, side, shift=0.0)
    first, second = extended_interval(edge, radius)
    width = 3radius
    p0 = add(edge.point, scale(first, edge.tangent))
    p1 = add(edge.point, scale(second, edge.tangent))
    q0 = add(p0, scale(side * shift, edge.gap))
    q1 = add(p1, scale(side * shift, edge.gap))
    r1 = add(p1, scale(side * width, edge.gap))
    r0 = add(p0, scale(side * width, edge.gap))
    return (q0, q1, r1, r0)
end

function convex_polygons_overlap(first, second, tolerance)
    for polygon in (first, second)
        for index in eachindex(polygon)
            start = polygon[index]
            stop = polygon[mod1(index + 1, length(polygon))]
            delta = (stop[1] - start[1], stop[2] - start[2])
            axis = (-delta[2], delta[1])
            norm = hypot(axis...)
            norm <= tolerance && continue
            axis = (axis[1] / norm, axis[2] / norm)
            first_projection = [point[1] * axis[1] + point[2] * axis[2] for point in first]
            second_projection =
                [point[1] * axis[1] + point[2] * axis[2] for point in second]
            if maximum(first_projection) <= minimum(second_projection) + tolerance ||
               maximum(second_projection) <= minimum(first_projection) + tolerance
                return false
            end
        end
    end
    return true
end

function validate_plan_view_geometry(edges, radius, tolerance, facets)
    !isempty(facets) && return
    polygons = [strip_points(edge, radius, -1.0) for edge in edges]
    for first in eachindex(edges), second = (first + 1):length(edges)
        a = edges[first]
        b = edges[second]
        same_layer =
            abs(a.point[3] - b.point[3]) <= tolerance && a.normal_sign == b.normal_sign
        if same_layer &&
           a.conductor != b.conductor &&
           convex_polygons_overlap(polygons[first], polygons[second], tolerance)
            error(
                "The edge-only spatial signature reconstructs overlapping " *
                "plan-view metal for conductors $(a.conductor) and " *
                "$(b.conductor). Exact coupon generation requires additional " *
                "plan-view conductor boundaries."
            )
        end
    end
end

function point_in_polygon(point, polygon, tolerance)
    inside = false
    previous = polygon[end]
    for current in polygon
        edge = (current[1] - previous[1], current[2] - previous[2])
        relative = (point[1] - previous[1], point[2] - previous[2])
        cross = edge[1] * relative[2] - edge[2] * relative[1]
        projection = relative[1] * edge[1] + relative[2] * edge[2]
        edge_norm_squared = edge[1]^2 + edge[2]^2
        if abs(cross) <= tolerance * max(hypot(edge...), 1.0) &&
           -tolerance <= projection <= edge_norm_squared + tolerance
            return true
        end
        if (previous[2] > point[2]) != (current[2] > point[2])
            intersection = previous[1] + (point[2] - previous[2]) * edge[1] / edge[2]
            intersection > point[1] && (inside = !inside)
        end
        previous = current
    end
    return inside
end

# The circular segments between the chords of the tagged arcs and their circles, registered
# from the plan-view boundary (register_mask_arc_segments!): the mask facets are the chorded
# polygons, the metal is bounded by the exact arcs, so a point inside a circular segment of
# conductor / plane flips the facet verdict (design 1.2 (1): one source of truth for the
# metal; a point of the exact metal wall of a convex arc lies outside every chord facet by
# up to the sagitta). Rows: (conductor, plane, chord start, chord end, outward normal of the
# chord away from the centre, centre, radius).
const MASK_ARC_SEGMENTS = NamedTuple[]

function register_mask_arc_segments!(loops, tolerance)
    empty!(MASK_ARC_SEGMENTS)
    for loop in loops
        loop_has_arcs(loop) || continue
        for run in tagged_arc_runs(loop, tolerance)
            for k in 1:(length(run.point_indices) - 1)
                a = loop.points[run.point_indices[k]]
                b = loop.points[run.point_indices[k + 1]]
                mid = (0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]))
                outward = (mid[1] - run.center[1], mid[2] - run.center[2])
                scale = hypot(outward...)
                scale > tolerance || continue
                push!(MASK_ARC_SEGMENTS, (conductor=loop.conductor, plane=loop.plane, a=a, b=b,
                                          normal=(outward[1] / scale, outward[2] / scale),
                                          centre=run.center, radius=run.radius, sign=run.sign))
            end
        end
    end
    return length(MASK_ARC_SEGMENTS)
end

# The exact-arc verdict for a plan-view point near a registered arc of conductor / plane:
# `true` on the arc itself (within the tolerance of the circle, over the chord's extent: the
# metal boundary belongs to the metal), `true` / `false` in a circular segment (strictly
# inside the circle, beyond the chord away from the centre) when the dielectric lies
# outside (sign +1: the segment is metal) / inside the circle (sign -1: dielectric); nothing
# away from every registered arc (the facet verdict stands).
function point_in_mask_arc_segment(point, conductor, plane, tolerance)
    for row in MASK_ARC_SEGMENTS
        row.conductor == conductor && abs(row.plane - plane) <= tolerance || continue
        distance = hypot(point[1] - row.centre[1], point[2] - row.centre[2])
        distance <= row.radius + tolerance || continue
        beyond = (point[1] - row.a[1]) * row.normal[1] + (point[2] - row.a[2]) * row.normal[2]
        beyond > -tolerance || continue
        # Within the chord's extent along the chord direction.
        d = (row.b[1] - row.a[1], row.b[2] - row.a[2])
        t = ((point[1] - row.a[1]) * d[1] + (point[2] - row.a[2]) * d[2]) / (d[1]^2 + d[2]^2)
        -1.0e-9 <= t <= 1.0 + 1.0e-9 || continue
        abs(distance - row.radius) <= tolerance && return true
        beyond > tolerance || continue
        return row.sign == 1
    end
    return nothing
end

function point_in_mask(facets, point, conductor, plane, tolerance)
    inside = any(
        facet.conductor == conductor &&
        abs(facet.plane - plane) <= tolerance &&
        point_in_polygon((point[1], point[2]), facet.points, tolerance) for facet in facets
    )
    isempty(MASK_ARC_SEGMENTS) && return inside
    verdict = point_in_mask_arc_segment(point, conductor, plane, tolerance)
    return verdict === nothing ? inside : verdict
end

function circle_through(first, second, third, tolerance)
    ax = second[1] - first[1]
    ay = second[2] - first[2]
    bx = third[1] - first[1]
    by = third[2] - first[2]
    determinant = 2.0 * (ax * by - ay * bx)
    scale = max(hypot(ax, ay), hypot(bx, by), 1.0)
    abs(determinant) > tolerance * scale || return nothing
    a2 = ax^2 + ay^2
    b2 = bx^2 + by^2
    center = (
        first[1] + (by * a2 - ay * b2) / determinant,
        first[2] + (ax * b2 - bx * a2) / determinant
    )
    radius = hypot(first[1] - center[1], first[2] - center[2])
    radius > tolerance || return nothing
    return (center=center, radius=radius)
end

function fitted_arc_run(points, point_indices, edge_indices, circle, tolerance)
    circle === nothing && return nothing
    radial = [
        (points[index][1] - circle.center[1], points[index][2] - circle.center[2]) for
        index in point_indices
    ]
    angle_steps = [
        atan(
            cross2d(radial[index], radial[index + 1]),
            radial[index][1] * radial[index + 1][1] +
            radial[index][2] * radial[index + 1][2]
        ) for index = 1:(length(radial) - 1)
    ]
    orientation = sign(sum(angle_steps))
    orientation != 0.0 &&
    all(sign(angle) == orientation for angle in angle_steps if angle != 0.0) ||
        return nothing
    residual = maximum(
        abs(
            hypot(
                points[index][1] - circle.center[1],
                points[index][2] - circle.center[2]
            ) - circle.radius
        ) for index in point_indices
    )
    residual <= max(64tolerance, 2.0e-7 * circle.radius) || return nothing
    return (
        center=circle.center,
        radius=circle.radius,
        point_indices=point_indices,
        edge_indices=edge_indices,
        orientation=orientation,
        angle=sum(abs, angle_steps)
    )
end

function circular_arc_runs(points, tolerance)
    length(points) >= 4 || return NamedTuple[]
    circles = [
        circle_through(
            points[index],
            points[mod1(index + 1, length(points))],
            points[mod1(index + 2, length(points))],
            tolerance
        ) for index in eachindex(points)
    ]
    compatible(first, second) =
        first !== nothing &&
        second !== nothing &&
        hypot(first.center[1] - second.center[1], first.center[2] - second.center[2]) <=
        max(32tolerance, 1.0e-7 * max(first.radius, second.radius)) &&
        abs(first.radius - second.radius) <=
        max(32tolerance, 1.0e-7 * max(first.radius, second.radius))

    compatible_pairs = [
        compatible(circles[index], circles[mod1(index + 1, length(points))]) for
        index in eachindex(points)
    ]
    any(compatible_pairs) || return NamedTuple[]
    if all(compatible_pairs)
        # Four rectangle vertices are also co-circular. A process-generated smooth
        # closed curve has many more samples, so reconstruct only sufficiently dense
        # loops and retain ordinary low-sided polygons exactly.
        length(points) >= 8 || return NamedTuple[]
        point_indices = vcat(collect(eachindex(points)), firstindex(points))
        run = fitted_arc_run(
            points,
            point_indices,
            collect(eachindex(points)),
            circles[firstindex(points)],
            tolerance
        )
        return run === nothing ? NamedTuple[] : [run]
    end

    runs = NamedTuple[]
    for seed in eachindex(points)
        compatible_pairs[seed] || continue
        compatible_pairs[mod1(seed - 1, length(points))] && continue
        pair_count = 1
        while pair_count < length(points) &&
            compatible_pairs[mod1(seed + pair_count, length(points))]
            pair_count += 1
        end
        triple_count = pair_count + 1
        edge_count = triple_count + 1
        # Four edges is the smallest useful smooth reconstruction. Requiring this
        # rejects accidental co-circular closure vertices without affecting the
        # process-rounded chains, which are exported at substantially higher resolution.
        edge_count >= 4 || continue
        point_indices = [mod1(seed + step, length(points)) for step = 0:(triple_count + 1)]
        edge_indices = [mod1(seed + step, length(points)) for step = 0:triple_count]
        circle = circle_through(
            points[first(point_indices)],
            points[point_indices[cld(length(point_indices), 2)]],
            points[last(point_indices)],
            tolerance
        )
        run = fitted_arc_run(points, point_indices, edge_indices, circle, tolerance)
        run === nothing || push!(runs, run)
    end
    return runs
end

# The OCC wire of a plan-view polygon: straight sides as lines, every recognised circular
# arc run as exact OCC circle arcs. `runs` = the loop's runs (loop_arc_runs: the TAGGED runs
# of a tagged metal loop, else the untagged fit; an offset loop has none to pass and is
# fitted here). A run is split into arc_part_count(sweep) parts of EQUAL angular fraction
# at new points on the circle (not at chord vertices): the split angles depend on the run's
# two end vertices alone, so the metal loft, the trench wall under it (the same end
# vertices, the concentric circle) and the arc tubes (ArcTube parts: the same rule) are
# split at the same radial planes and the fragment keeps every tube entity whole, whatever
# the chord count (design 1.2 (3), A7 MINOR-7). A straight loop has no run: bitwise.
function polygon_wire(occ, points, z; runs=nothing)
    tolerance =
        1.0e-9 * max(
            maximum(point[1] for point in points) - minimum(point[1] for point in points),
            maximum(point[2] for point in points) - minimum(point[2] for point in points),
            1.0
        )
    runs = runs === nothing ? circular_arc_runs(points, tolerance) : runs
    projected = collect(points)
    for run in runs, index in run.point_indices[2:(end - 1)]
        radial = (points[index][1] - run.center[1], points[index][2] - run.center[2])
        scale = run.radius / hypot(radial...)
        projected[index] =
            (run.center[1] + scale * radial[1], run.center[2] + scale * radial[2])
    end
    # OCC points are created where a curve needs them: every vertex of a straight side and
    # the two ends of a run (a chord vertex strictly inside a run would be a dangling CAD
    # point that Gmsh would still carry as a node). A straight loop creates every vertex.
    tags = Dict{Int, Int32}()
    point_tag(index) = get!(tags, index) do
        occ.addPoint(projected[index][1], projected[index][2], z)
    end
    curve_for_edge = Dict{Int, Int32}()
    covered = falses(length(points))
    for run in runs
        any(covered[run.edge_indices]) && continue
        sweep = run_sweep(points, run.point_indices, run.center, run.orientation)
        angles = arc_split_angles(points, run.point_indices, run.center, sweep)
        parts = length(angles) - 1
        parts <= length(run.edge_indices) ||
            error("a circular arc run of $(length(run.edge_indices)) chords cannot carry $parts parts")
        center = occ.addPoint(run.center[1], run.center[2], z)
        split_tags = vcat(
            [point_tag(run.point_indices[1])],
            [occ.addPoint(run.center[1] + run.radius * cos(angles[k]),
                          run.center[2] + run.radius * sin(angles[k]), z) for k in 2:parts],
            [point_tag(run.point_indices[end])]
        )
        for k in 1:parts
            curve_for_edge[run.edge_indices[k]] =
                occ.addCircleArc(split_tags[k], center, split_tags[k + 1])
        end
        covered[run.edge_indices] .= true
    end
    curves = Int32[]
    for index in eachindex(points)
        if haskey(curve_for_edge, index)
            push!(curves, curve_for_edge[index])
        elseif !covered[index]
            push!(curves, occ.addLine(point_tag(index), point_tag(mod1(index + 1, length(points)))))
        end
    end
    return occ.addWire(curves)
end

cross2d(first, second) = first[1] * second[2] - first[2] * second[1]

function loop_orientation(points)
    area2=sum(points[i][1]*points[mod1(i+1,length(points))][2]-
              points[mod1(i+1,length(points))][1]*points[i][2] for i in eachindex(points))
    area2!=0 || error("Degenerate plan-view loop")
    return sign(area2)
end

# The plan-view loop offset by `distance` along the metal-side normal of every
# Physical side (`distance` < 0: away from the metal; Continuation sides, on the
# box, never move). A straight loop is the miter polygon of its shifted sides; a
# loop with recognised circular arcs (circular_arc_runs) is offset exactly by
# curved_offset_loop: concentric arcs, tangent joints, collapsed degenerate arcs. A
# bridged curved offset (a collapsed arc whose neighbours' offsets do not meet) is
# never a simple offset: it fails closed here; only collar_loop_points, which
# routes it to the collar union, accepts it.
function offset_loop_points(loop, distance, tolerance)
    abs(distance)<=tolerance && return loop.points
    runs = loop_arc_runs(loop, tolerance)
    if !isempty(runs)
        offset = curved_offset_loop(loop, distance, runs, tolerance)
        offset.bridged &&
            error("Offset $distance of the plan-view loop of conductor $(loop.conductor) on " *
                  "plane $(loop.plane) bridges a collapsed circular arc whose neighbours' " *
                  "offsets do not meet: not a simple offset (only the collar union takes it)")
        return offset.points
    end
    metal_side=loop_orientation(loop.points)*(loop.hole ? -1.0 : 1.0)
    shifted = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for index in eachindex(loop.points)
        first = loop.points[index]
        second = loop.points[mod1(index + 1, length(loop.points))]
        direction = (second[1] - first[1], second[2] - first[2])
        segment_length = hypot(direction...)
        segment_length > tolerance ||
            error("Plan-view boundary contains a zero-length segment")
        shift = loop.classes[index] == "Physical" ? distance : 0.0
        normal = (-metal_side*direction[2] / segment_length,
                   metal_side*direction[1] / segment_length)
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
        abs(denominator) > tolerance ||
            error("Plan-view taper has a singular boundary vertex")
        offset = (current[1][1] - previous[1][1], current[1][2] - previous[1][2])
        coordinate = cross2d(offset, current[2]) / denominator
        point = (
            previous[1][1] + coordinate * previous[2][1],
            previous[1][2] + coordinate * previous[2][2]
        )
        hypot(point[1] - loop.points[index][1], point[2] - loop.points[index][2]) <=
        8.0 * max(abs(distance), tolerance) ||
            error("Plan-view taper produces an unresolved miter")
        push!(points, point)
    end
    return points
end

# ---------------------------------------------------------------------------
# Exact offsets of curved plan-view loops (supervisor decision 246).
#
# A loop with recognised circular arcs is a sequence of items in loop order: a
# straight side (one edge) or an arc run (its consecutive chord edges). Under an
# offset `distance` (< 0: away from the metal) a Physical straight side shifts
# along its metal-side normal as in offset_loop_points; a Physical arc becomes the
# concentric arc of radius r - s x distance, s = +1 for a convex-metal arc (the
# metal inside its circle, the arc bulging away from the metal) and -1 for a
# concave-metal arc, keeping the run's chord vertices as exact points of the new
# circle (polygon_wire rebuilds the OCC arc from them, like the mask's own arcs).
# A concave arc smaller than the offset (r - |distance| <= tolerance) is
# degenerate: it collapses to the junction of its two neighbours, the exact
# boundary of the dilated region there. Junctions: two lines meet at their miter;
# a line and an arc meet where the shifted line meets the new circle (the foot of
# the centre on the line when they are tangent, else the intersection nearest the
# arc's shifted end); two arcs must stay tangent (their shifted ends agree). A
# collapsed arc whose neighbours' offsets do not meet (antiparallel channel walls,
# an unresolved miter) is bridged by the two shifted ends and the polygon is marked
# `bridged`: it is never a simple offset and goes to the collar union.

# Tolerance of the arc fit of circular_arc_runs (fitted_arc_run), at radius `r`.
arc_fit_tolerance(r, tolerance) = max(64tolerance, 2.0e-7 * r)

# A junction whose two items' directions of travel differ by at most this angle
# (radians) is tangent: the identification's joint noise (CSV rounding of the arc
# vertices, ~1e-6) never makes a corner of a rounded joint, and a genuine corner
# turns by far more. Its two shifted ends (apart by |distance| x the angle) snap to
# the junction point; a larger turn is a corner (convex: a miter kite in the union).
const JUNCTION_TANGENT_ANGLE = 1.0e-4

perp2d(v) = (-v[2], v[1])

unit2d(v) = (v[1] / hypot(v...), v[2] / hypot(v...))

# Items of the loop in loop order, rotated so that no arc run straddles the end of
# the sequence: (kind, edges, run) with `edges` the 1-based edge indices (edge i
# joins points[i] and points[i + 1]) and `run` the index into `runs` (0: a line).
function loop_edge_items(points, runs)
    n = length(points)
    run_of_edge = zeros(Int, n)
    for (k, run) in enumerate(runs), index in run.edge_indices
        run_of_edge[index] == 0 ||
            error("Plan-view boundary edge $index belongs to two circular arcs")
        run_of_edge[index] = k
    end
    first_edge_of_run = Dict(k => first(run.edge_indices) for (k, run) in enumerate(runs))
    start = findfirst(i -> run_of_edge[i] == 0 || first_edge_of_run[run_of_edge[i]] == i, 1:n)
    start === nothing && error("Plan-view boundary has no arc start")
    items = NamedTuple[]
    index = start
    for _ in 1:n
        edge = mod1(index, n)
        k = run_of_edge[edge]
        if k == 0
            push!(items, (kind=:line, edges=[edge], run=0))
            index += 1
        else
            edges = copy(runs[k].edge_indices)
            edges[1] == edge ||
                error("Circular arc of the plan-view boundary starts inside another item")
            push!(items, (kind=:arc, edges=edges, run=k))
            index += length(edges)
        end
        mod1(index, n) == start && break
    end
    sum(length(item.edges) for item in items) == n ||
        error("Plan-view boundary items do not cover the loop once")
    return items
end

# The offset geometry of one item (see curved_offset_loop): `start` / `stop` are
# the item's vertices, `shifted_start` / `shifted_stop` their offsets along the
# item, `tangent_start` / `tangent_stop` the directions of travel there and, for
# an arc, `interior` the offset chord vertices between them.
function offset_item_geometry(item, loop, runs, distance, metal_side, tolerance)
    points = loop.points
    n = length(points)
    start = points[item.edges[1]]
    stop = points[mod1(item.edges[end] + 1, n)]
    classes = unique(loop.classes[item.edges])
    length(classes) == 1 || error("A circular arc of the plan-view boundary of conductor " *
                                   "$(loop.conductor) mixes Physical and Continuation chords")
    shift = classes[1] == "Physical" ? distance : 0.0
    if item.kind == :line
        direction = (stop[1] - start[1], stop[2] - start[2])
        segment_length = hypot(direction...)
        segment_length > tolerance ||
            error("Plan-view boundary contains a zero-length segment")
        normal = (-metal_side * direction[2] / segment_length,
                  metal_side * direction[1] / segment_length)
        offset = (shift * normal[1], shift * normal[2])
        return (kind=:line, edges=item.edges, start=start, stop=stop,
                shifted_start=(start[1] + offset[1], start[2] + offset[2]),
                shifted_stop=(stop[1] + offset[1], stop[2] + offset[2]),
                direction=direction, tangent_start=direction, tangent_stop=direction,
                collapsed=false, interior=NTuple{2, Float64}[], shift=offset)
    end
    run = runs[item.run]
    convex_sign = run.orientation * metal_side
    radius = run.radius - convex_sign * shift
    scaled(p) = (run.center[1] + radius * (p[1] - run.center[1]) /
                                 hypot(p[1] - run.center[1], p[2] - run.center[2]),
                 run.center[2] + radius * (p[2] - run.center[2]) /
                                 hypot(p[1] - run.center[1], p[2] - run.center[2]))
    tangent(p) = run.orientation .* perp2d((p[1] - run.center[1], p[2] - run.center[2]))
    collapsed = radius <= tolerance
    interior = collapsed ? NTuple{2, Float64}[] :
               [scaled(points[edge]) for edge in item.edges[2:end]]
    return (kind=:arc, edges=item.edges, start=start, stop=stop,
            shifted_start=collapsed ? start : scaled(start),
            shifted_stop=collapsed ? stop : scaled(stop),
            center=run.center, original_radius=run.radius, radius=radius,
            orientation=run.orientation, convex_sign=convex_sign,
            tangent_start=tangent(start), tangent_stop=tangent(stop),
            collapsed=collapsed, interior=interior, angle=run.angle)
end

# Where the shifted line (`point`, `direction`) meets the offset circle of `arc`:
# the foot of the centre when they are tangent (the arc-fit tolerance), else the
# intersection nearest `near`; nothing when the line misses the circle.
function line_circle_junction(point, direction, arc, near, tolerance)
    u = unit2d(direction)
    w = (arc.center[1] - point[1], arc.center[2] - point[2])
    h = cross2d(u, w)
    along = w[1] * u[1] + w[2] * u[2]
    foot = (point[1] + along * u[1], point[2] + along * u[2])
    slack = arc_fit_tolerance(arc.radius, tolerance)
    abs(h) > arc.radius + slack && return nothing
    arc.radius - abs(h) <= slack && return foot
    half = sqrt(arc.radius^2 - h^2)
    candidates = ((foot[1] + half * u[1], foot[2] + half * u[2]),
                  (foot[1] - half * u[1], foot[2] - half * u[2]))
    return argmin(c -> hypot(c[1] - near[1], c[2] - near[2]), candidates)
end

# The junction point of two consecutive non-collapsed items (nothing when their
# offsets do not meet).
function offset_junction(before, after, distance, tolerance)
    if before.kind == :line && after.kind == :line
        denominator = cross2d(before.direction, after.direction)
        abs(denominator) > tolerance * hypot(before.direction...) * hypot(after.direction...) ||
            return nothing
        offset = (after.shifted_start[1] - before.shifted_start[1],
                  after.shifted_start[2] - before.shifted_start[2])
        coordinate = cross2d(offset, after.direction) / denominator
        return (before.shifted_start[1] + coordinate * before.direction[1],
                before.shifted_start[2] + coordinate * before.direction[2])
    elseif before.kind == :line
        return line_circle_junction(before.shifted_start, before.direction, after,
                                    after.shifted_start, tolerance)
    elseif after.kind == :line
        return line_circle_junction(after.shifted_start, after.direction, before,
                                    before.shifted_stop, tolerance)
    end
    gap = hypot(before.shifted_stop[1] - after.shifted_start[1],
                before.shifted_stop[2] - after.shifted_start[2])
    gap <= arc_fit_tolerance(max(before.radius, after.radius), tolerance) ||
        error("Two circular arcs of the plan-view boundary meet at a kink (their offsets " *
              "by $distance differ by $gap): only tangent arc-arc joints are supported")
    return (0.5 * (before.shifted_stop[1] + after.shifted_start[1]),
            0.5 * (before.shifted_stop[2] + after.shifted_start[2]))
end

# The offset polygon of a curved loop: `points` (the miter polygon with exact arc
# vertices), `bridged` (a collapsed arc's neighbours did not meet: not a simple
# offset), `items` (offset_item_geometry per item) and `junctions`, one per pair of
# consecutive non-collapsed items: (before, after) item indices, `vertex` (the
# loop vertex, or the two vertices of the collapsed arcs between them), `point`
# (nothing when bridged), `tangent` (the items join within JUNCTION_TANGENT_ANGLE)
# and `convex` (the metal turns convexly there).
function curved_offset_loop(loop, distance, runs, tolerance)
    points = loop.points
    n = length(points)
    metal_side = loop_orientation(points) * (loop.hole ? -1.0 : 1.0)
    if length(runs) == 1 && length(runs[1].edge_indices) == n
        # A whole circle: concentric, or collapsed (no neighbour to collapse onto).
        run = runs[1]
        classes = unique(loop.classes)
        length(classes) == 1 ||
            error("A circular plan-view loop mixes Physical and Continuation chords")
        radius = run.radius - run.orientation * metal_side *
                              (classes[1] == "Physical" ? distance : 0.0)
        radius > tolerance ||
            error("Offset $distance collapses the circular plan-view loop of conductor " *
                  "$(loop.conductor) on plane $(loop.plane) (radius $(run.radius))")
        offset = [(run.center[1] + radius * (p[1] - run.center[1]) /
                                   hypot(p[1] - run.center[1], p[2] - run.center[2]),
                   run.center[2] + radius * (p[2] - run.center[2]) /
                                   hypot(p[1] - run.center[1], p[2] - run.center[2]))
                  for p in points]
        return (points=offset, bridged=false, items=NamedTuple[], junctions=NamedTuple[],
                metal_side=metal_side)
    end
    items = [offset_item_geometry(item, loop, runs, distance, metal_side, tolerance)
             for item in loop_edge_items(points, runs)]
    live = [k for k in eachindex(items) if !items[k].collapsed]
    isempty(live) && error("Offset $distance collapses every arc of the plan-view loop of " *
                           "conductor $(loop.conductor) on plane $(loop.plane)")
    offset = NTuple{2, Float64}[]
    junctions = NamedTuple[]
    bridged = false
    for (position, k) in enumerate(live)
        before = items[k]
        next = live[mod1(position + 1, length(live))]
        after = items[next]
        skipped = Int[]
        j = mod1(k + 1, length(items))
        while j != next
            push!(skipped, j)
            j = mod1(j + 1, length(items))
        end
        collapsed_between = !isempty(skipped)
        vertices = collapsed_between ? [before.stop, after.start] : [before.stop]
        point = offset_junction(before, after, distance, tolerance)
        resolved = point !== nothing && any(
            hypot(point[1] - v[1], point[2] - v[2]) <= 8.0 * max(abs(distance), tolerance)
            for v in vertices)
        if !resolved && !collapsed_between
            if point === nothing && before.kind == :line && after.kind == :line
                error("Plan-view taper has a singular boundary vertex at $(before.stop)")
            elseif point === nothing
                error("Plan-view taper offset by $distance: the shifted straight side misses " *
                      "the offset circle of the circular arc it meets at $(before.stop) (a " *
                      "line-arc kink whose offsets do not meet)")
            end
            error("Plan-view taper produces an unresolved miter at $(before.stop)")
        end
        turn = metal_side * atan(cross2d(before.tangent_stop, after.tangent_start),
                                 before.tangent_stop[1] * after.tangent_start[1] +
                                 before.tangent_stop[2] * after.tangent_start[2])
        tangent = !collapsed_between && abs(turn) <= JUNCTION_TANGENT_ANGLE
        convex = !collapsed_between && turn > JUNCTION_TANGENT_ANGLE
        append!(offset, before.interior)
        if resolved
            push!(offset, point)
        else
            bridged = true
            push!(offset, before.shifted_stop)
            push!(offset, after.shifted_start)
        end
        push!(junctions, (before=k, after=next, vertices=vertices,
                          point=resolved ? point : nothing, tangent=tangent, convex=convex,
                          collapsed=skipped))
    end
    # The cycle [interior_1, J_12, interior_2, ..., interior_m, J_m1] is the offset
    # polygon in loop order (J_m1 closes onto interior_1).
    return (points=offset, bridged=bridged, items=items, junctions=junctions,
            metal_side=metal_side)
end

# A closed plan-view polygon is simple when no side is degenerate (shorter than
# `tolerance`, or folding straight back onto its predecessor) and no two
# non-adjacent sides come within `tolerance` of each other.
function polygon_is_simple(points, tolerance)
    n = length(points)
    n >= 3 || return false
    side(i) = (points[mod1(i, n)], points[mod1(i + 1, n)])
    for i in 1:n
        a, b = side(i)
        u = (b[1] - a[1], b[2] - a[2])
        hypot(u...) > tolerance || return false
        c, d = side(i + 1)
        v = (d[1] - c[1], d[2] - c[2])
        u[1] * v[1] + u[2] * v[2] < 0.0 && abs(cross2d(u, v)) <= tolerance * hypot(u...) &&
            return false
    end
    for i in 1:n, j in (i + 2):n
        (i == 1 && j == n) && continue
        segment_segment_distance_2d(side(i)..., side(j)...) <= tolerance && return false
    end
    return true
end

# Consecutive duplicate vertices (within `tolerance`) removed, the closing
# duplicate included.
function dedupe_polygon_points(points, tolerance)
    result = NTuple{2, Float64}[]
    for point in points
        (isempty(result) || hypot(point[1] - result[end][1], point[2] - result[end][2]) > tolerance) &&
            push!(result, (Float64(point[1]), Float64(point[2])))
    end
    while length(result) > 1 &&
          hypot(result[1][1] - result[end][1], result[1][2] - result[end][2]) <= tolerance
        pop!(result)
    end
    return result
end

# Sutherland-Hodgman clip of a polygon to the plan-view rectangle of the coupon box
# (vertices on a box side within `tolerance` are kept; cut vertices lie exactly on
# the side). Returns the deduplicated result, possibly empty.
function clip_polygon_to_box(points, lower, upper, tolerance)
    result = collect(points)
    for (axis, bound, sign) in ((1, lower[1], 1.0), (1, upper[1], -1.0),
                                (2, lower[2], 1.0), (2, upper[2], -1.0))
        inside(p) = sign * (p[axis] - bound) >= -tolerance
        output = NTuple{2, Float64}[]
        m = length(result)
        for k in 1:m
            p = result[k]
            q = result[mod1(k + 1, m)]
            inside(p) && push!(output, (Float64(p[1]), Float64(p[2])))
            if inside(p) != inside(q)
                t = (bound - p[axis]) / (q[axis] - p[axis])
                push!(output, axis == 1 ? (bound, p[2] + t * (q[2] - p[2])) :
                                          (p[1] + t * (q[1] - p[1]), bound))
            end
        end
        result = output
    end
    return dedupe_polygon_points(result, tolerance)
end

polygon_area2(points) = sum(cross2d(points[i], points[mod1(i + 1, length(points))])
                            for i in eachindex(points); init=0.0)

# The convex pieces whose union is the etch collar of an exterior loop under the
# miter rule of offset_loop_points (`distance` < 0: away from the metal): the loop
# itself, one rectangle per shifted (Physical) side and, at every convex metal
# corner, the kite between the two rectangle ends and the miter point (`miter`
# is the miter polygon; a right-angle box junction gives a degenerate kite). Every
# piece is clipped to the coupon box: the region outside the box is never lofted.
# A loop with circular arcs takes curved_collar_pieces.
function collar_pieces(loop, distance, miter, lower, upper, tolerance)
    points = loop.points
    n = length(points)
    runs = loop_arc_runs(loop, tolerance)
    isempty(runs) || return curved_collar_pieces(
        loop, distance, curved_offset_loop(loop, distance, runs, tolerance), lower, upper,
        tolerance)
    metal_side = loop_orientation(points) * (loop.hole ? -1.0 : 1.0)
    directions = Vector{NTuple{2, Float64}}(undef, n)
    shifts = Vector{NTuple{2, Float64}}(undef, n)
    for i in 1:n
        a = points[i]
        b = points[mod1(i + 1, n)]
        direction = (b[1] - a[1], b[2] - a[2])
        segment_length = hypot(direction...)
        shift = loop.classes[i] == "Physical" ? distance : 0.0
        directions[i] = direction
        shifts[i] = (-metal_side * direction[2] / segment_length * shift,
                     metal_side * direction[1] / segment_length * shift)
    end
    pieces = [[(Float64(p[1]), Float64(p[2])) for p in points]]
    for i in 1:n
        hypot(shifts[i]...) > tolerance || continue
        a = points[i]
        b = points[mod1(i + 1, n)]
        push!(pieces, [(a[1], a[2]), (b[1], b[2]), (b[1] + shifts[i][1], b[2] + shifts[i][2]),
                       (a[1] + shifts[i][1], a[2] + shifts[i][2])])
    end
    for i in 1:n
        h = mod1(i - 1, n)
        metal_side * cross2d(directions[h], directions[i]) >
        tolerance * hypot(directions[h]...) * hypot(directions[i]...) || continue
        p = points[i]
        push!(pieces, [(p[1], p[2]), (p[1] + shifts[h][1], p[2] + shifts[h][2]),
                       (miter[i][1], miter[i][2]), (p[1] + shifts[i][1], p[2] + shifts[i][2])])
    end
    clipped = Vector{NTuple{2, Float64}}[]
    for piece in pieces
        piece = clip_polygon_to_box(dedupe_polygon_points(piece, tolerance), lower, upper, tolerance)
        length(piece) >= 3 && abs(polygon_area2(piece)) > tolerance^2 || continue
        push!(clipped, polygon_area2(piece) > 0.0 ? piece : reverse(piece))
    end
    return clipped
end

# The pieces whose union is the etch collar of a curved exterior loop (`offset` =
# curved_offset_loop(loop, distance, ...), `distance` < 0): the loop, one rectangle
# per shifted straight side, one annular sector per offset arc (between the arc's
# chord vertices and their concentric offsets: the exact dilation of the arc), one
# fan (the arc's chord vertices and its centre: the whole sector is within the
# offset of a collapsed arc) per collapsed arc and the kite of every convex
# junction whose two shifted ends differ. Neighbouring pieces share their junction
# point exactly: a continuous (tangent) junction snaps both ends to it, a fan's
# apex is where its two neighbours' end normals meet (the centre, to the arc fit's
# noise), so that the union sees collinear shared sides, never a sliver.
function curved_collar_pieces(loop, distance, offset, lower, upper, tolerance)
    points = loop.points
    n = length(points)
    items = offset.items
    isempty(items) && error("The collar union of a circular plan-view loop is its concentric offset")
    start_corner = [item.shifted_start for item in items]
    stop_corner = [item.shifted_stop for item in items]
    kites = Vector{NTuple{2, Float64}}[]
    apex = Dict{Int, NTuple{2, Float64}}()
    for junction in offset.junctions
        before = items[junction.before]
        after = items[junction.after]
        if junction.point !== nothing && junction.tangent
            stop_corner[junction.before] = junction.point
            start_corner[junction.after] = junction.point
        elseif junction.point !== nothing && junction.convex
            push!(kites, [junction.vertices[1], before.shifted_stop, junction.point,
                          after.shifted_start])
        end
        if length(junction.collapsed) == 1
            k = junction.collapsed[1]
            arc = items[k]
            u = (before.shifted_stop[1] - before.stop[1], before.shifted_stop[2] - before.stop[2])
            v = (after.shifted_start[1] - after.start[1], after.shifted_start[2] - after.start[2])
            denominator = cross2d(u, v)
            if abs(denominator) > tolerance * hypot(u...) * hypot(v...)
                w = (arc.stop[1] - arc.start[1], arc.stop[2] - arc.start[2])
                t = cross2d(w, v) / denominator
                candidate = (arc.start[1] + t * u[1], arc.start[2] + t * u[2])
                hypot(candidate[1] - arc.center[1], candidate[2] - arc.center[2]) <=
                1.0e-3 * arc.original_radius && (apex[k] = candidate)
            end
        end
    end
    pieces = [[(Float64(p[1]), Float64(p[2])) for p in points]]
    for (k, item) in enumerate(items)
        chords = [points[edge] for edge in item.edges[2:end]]
        if item.kind == :line
            hypot(item.shift...) > tolerance || continue
            push!(pieces, [item.start, item.stop, stop_corner[k], start_corner[k]])
        elseif item.collapsed
            push!(pieces, vcat([item.start], chords, [item.stop, get(apex, k, item.center)]))
        else
            abs(item.radius - item.original_radius) > tolerance || continue
            push!(pieces, vcat([item.start], chords, [item.stop, stop_corner[k]],
                               reverse(item.interior), [start_corner[k]]))
        end
    end
    append!(pieces, kites)
    clipped = Vector{NTuple{2, Float64}}[]
    for piece in pieces
        piece = clip_polygon_to_box(dedupe_polygon_points(piece, tolerance), lower, upper, tolerance)
        length(piece) >= 3 && abs(polygon_area2(piece)) > tolerance^2 || continue
        push!(clipped, polygon_area2(piece) > 0.0 ? piece : reverse(piece))
    end
    return clipped
end

# Outer boundary of the union of simple polygons: every piece edge is split where
# any other piece edge meets it (crossings, touching ends and collinear overlaps),
# a sub-segment is on the union boundary when exactly one of its two sides is
# inside some piece (probed a quarter of the shortest sub-segment away), and the
# boundary sub-segments are chained with the union on their left, starting at the
# lexicographically smallest vertex. Fails closed (ScopeGuard FootprintTopology)
# when the boundary is not one simple loop: a vertex with two outgoing boundary
# sub-segments (the region touches itself) or sub-segments left over (holes) —
# unless `island_rule` admits them (collar_island_rule / absorb_collar_islands:
# the records of the absorbed islands are pushed to `absorbed`).
function polygon_union_boundary(pieces, tolerance; island_rule=nothing, absorbed=nothing)
    edges = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for piece in pieces, k in eachindex(piece)
        push!(edges, (piece[k], piece[mod1(k + 1, length(piece))]))
    end
    subsegments = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for (index, (a, b)) in enumerate(edges)
        u = (b[1] - a[1], b[2] - a[2])
        len = hypot(u...)
        parameters = [0.0, 1.0]
        for (other, (c, d)) in enumerate(edges)
            other == index && continue
            v = (d[1] - c[1], d[2] - c[2])
            w = (c[1] - a[1], c[2] - a[2])
            denominator = cross2d(u, v)
            if abs(denominator) <= tolerance * len * hypot(v...)
                abs(cross2d(u, w)) <= tolerance * len || continue
                for point in (c, d)
                    t = ((point[1] - a[1]) * u[1] + (point[2] - a[2]) * u[2]) / len^2
                    -tolerance / len < t < 1.0 + tolerance / len &&
                        push!(parameters, clamp(t, 0.0, 1.0))
                end
            else
                t = cross2d(w, v) / denominator
                s = cross2d(w, u) / denominator
                -tolerance / len <= t <= 1.0 + tolerance / len &&
                    -tolerance / hypot(v...) <= s <= 1.0 + tolerance / hypot(v...) &&
                    push!(parameters, clamp(t, 0.0, 1.0))
            end
        end
        sort!(parameters)
        merged = [parameters[1]]
        for t in parameters
            t - merged[end] > tolerance / len && push!(merged, t)
        end
        merged[end] < 1.0 && (merged[end] = 1.0)
        at(t) = (a[1] + t * u[1], a[2] + t * u[2])
        for k in 1:(length(merged) - 1)
            push!(subsegments, (at(merged[k]), at(merged[k + 1])))
        end
    end
    probe = 0.25 * minimum(hypot(b[1] - a[1], b[2] - a[2]) for (a, b) in subsegments)
    inside_union(p) = any(point_in_polygon(p, piece, 1.0e-3 * probe) for piece in pieces)
    boundary = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for (a, b) in subsegments
        u = (b[1] - a[1], b[2] - a[2])
        len = hypot(u...)
        normal = (-u[2] / len, u[1] / len)
        mid = (0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]))
        left = inside_union((mid[1] + probe * normal[1], mid[2] + probe * normal[2]))
        right = inside_union((mid[1] - probe * normal[1], mid[2] - probe * normal[2]))
        left == right && continue
        segment = left ? (a, b) : (b, a)
        same(p, q) = hypot(p[1] - q[1], p[2] - q[2]) <= tolerance
        any(same(s[1], segment[1]) && same(s[2], segment[2]) for s in boundary) && continue
        any(same(s[1], segment[2]) && same(s[2], segment[1]) for s in boundary) &&
            error("Collar union boundary sub-segment $segment has the union on both sides")
        push!(boundary, segment)
    end
    vertices = NTuple{2, Float64}[]
    function vertex_index(p)
        i = findfirst(v -> hypot(v[1] - p[1], v[2] - p[2]) <= tolerance, vertices)
        i === nothing || return i
        push!(vertices, p)
        return length(vertices)
    end
    starts = [vertex_index(a) for (a, _) in boundary]
    stops = [vertex_index(b) for (_, b) in boundary]
    outgoing = [findall(==(v), starts) for v in eachindex(vertices)]
    incoming = [findall(==(v), stops) for v in eachindex(vertices)]
    for v in eachindex(vertices)
        length(outgoing[v]) == 1 && length(incoming[v]) == 1 ||
            scope_error("FootprintTopology",
                        "the collar region touches itself at $(vertices[v]) " *
                        "($(length(outgoing[v])) outgoing / $(length(incoming[v])) incoming " *
                        "boundary sub-segments)")
    end
    start = argmin(vertices)
    loop = [start]
    used = falses(length(boundary))
    current = start
    while true
        k = outgoing[current][1]
        used[k] = true
        current = stops[k]
        current == start && break
        push!(loop, current)
    end
    all(used) || absorb_collar_islands!(absorbed, island_rule, vertices, starts, stops,
                                        outgoing, used, tolerance)
    return [vertices[i] for i in loop]
end

# Supervisor decision 246 (B): the 3R collar stands for the real overetch, which
# removes ALL exposed substrate, so an un-etched region created only by the collar
# geometry is an artefact. An un-etched region bounded entirely by collar
# boundaries (not touching the box face) whose every point lies within
# COLLAR_ISLAND_EXCESS_CAP_OVER_RADIUS x Radius beyond the collar is absorbed into
# the etched collar and recorded (count, area, maximum excess); a larger island, or
# one touching the box face, fails closed as before (ScopeGuard FootprintTopology).
const COLLAR_ISLAND_EXCESS_CAP_OVER_RADIUS = 0.05

# The island rule of a collar union: the coupon box, the metal boundary primitives
# of the loop (physical_segments at zero offset: exact arcs and lines), the collar
# width and the admitted excess.
function collar_island_rule(loop, distance, box, radius, tolerance)
    return (lower=box[1], upper=box[2], primitives=physical_segments([loop], 0.0, tolerance),
            collar=abs(distance), cap=COLLAR_ISLAND_EXCESS_CAP_OVER_RADIUS * radius)
end

# Distance from `point` to the metal boundary described by `primitives`.
metal_distance(point, primitives, tolerance) =
    minimum(point_primitive_distance(point, primitive, tolerance) for primitive in primitives)

# Half the smallest width of the island polygon over its edge normals: an upper
# bound of its inradius, hence of the excess of any of its points beyond the
# CONSTRUCTED collar boundary (every boundary point is on that boundary, at the
# collar distance from the metal along the rectangles and arcs but farther along a
# convex-corner miter kite, and the metal distance is 1-Lipschitz). The exact metal
# distance is measured by island_maximum_excess; the gate needs both within the cap.
function island_excess_bound(points)
    n = length(points)
    bound = Inf
    for i in 1:n
        a = points[i]
        b = points[mod1(i + 1, n)]
        normal = perp2d((b[1] - a[1], b[2] - a[2]))
        hypot(normal...) > 0.0 || continue
        normal = unit2d(normal)
        projections = [p[1] * normal[1] + p[2] * normal[2] for p in points]
        bound = min(bound, 0.5 * (maximum(projections) - minimum(projections)))
    end
    return bound
end

# The largest excess of an island point beyond the collar found deterministically:
# the best of the island's vertices and of a 32 x 32 grid over its bounding box,
# refined by a pattern search (8 directions, steps from a quarter of the island's
# extent down to 1e-6 of the collar) inside the island. The rigorous bound is
# island_excess_bound; this value is the recorded estimate.
function island_maximum_excess(points, rule, tolerance)
    excess(p) = metal_distance(p, rule.primitives, tolerance) - rule.collar
    inside(p) = point_in_polygon(p, points, tolerance)
    lower = (minimum(p[1] for p in points), minimum(p[2] for p in points))
    upper = (maximum(p[1] for p in points), maximum(p[2] for p in points))
    best, best_point = -Inf, points[1]
    candidates = vcat(points, [(lower[1] + (i + 0.5) / 32 * (upper[1] - lower[1]),
                                lower[2] + (j + 0.5) / 32 * (upper[2] - lower[2]))
                               for i in 0:31 for j in 0:31])
    for candidate in candidates
        inside(candidate) || continue
        value = excess(candidate)
        if value > best
            best, best_point = value, candidate
        end
    end
    step = 0.25 * max(upper[1] - lower[1], upper[2] - lower[2])
    directions = [(cos(k * pi / 4), sin(k * pi / 4)) for k in 0:7]
    while step > 1.0e-6 * rule.collar
        improved = true
        while improved
            improved = false
            for direction in directions
                candidate = (best_point[1] + step * direction[1], best_point[2] + step * direction[2])
                inside(candidate) || continue
                value = excess(candidate)
                if value > best
                    best, best_point, improved = value, candidate, true
                end
            end
        end
        step *= 0.5
    end
    return best, best_point
end

# The boundary sub-segments left over by the outer loop of polygon_union_boundary
# are the boundaries of un-etched islands. Each is chained (every vertex has one
# outgoing and one incoming sub-segment), judged by the island rule and recorded;
# any island touching the box face or exceeding the admitted excess fails closed.
function absorb_collar_islands!(absorbed, rule, vertices, starts, stops, outgoing, used,
                                tolerance)
    rule === nothing &&
        scope_error("FootprintTopology",
                    "the collar region encloses $(count(!, used)) boundary sub-segments " *
                    "of an un-etched island inside a dielectric gap")
    islands = Vector{NTuple{2, Float64}}[]
    while !all(used)
        start = starts[findfirst(!, used)]
        island = [vertices[start]]
        current = start
        while true
            k = outgoing[current][1]
            used[k] && error("Collar island boundary revisits $(vertices[current])")
            used[k] = true
            current = stops[k]
            current == start && break
            push!(island, vertices[current])
        end
        push!(islands, island)
    end
    on_box(p) = any(abs(p[d] - rule.lower[d]) <= tolerance || abs(p[d] - rule.upper[d]) <= tolerance
                    for d in 1:2)
    for island in islands
        any(on_box, island) &&
            scope_error("FootprintTopology",
                        "the collar region leaves an un-etched region of $(length(island)) " *
                        "boundary sub-segments touching the box face (the island rule " *
                        "applies inside the box only)")
        bound = island_excess_bound(island)
        excess, at = island_maximum_excess(island, rule, tolerance)
        # The bound measures beyond the CONSTRUCTED collar boundary (which carries
        # the convex-corner miter kites), the measurement the exact metal distance:
        # an island bounded by a kite edge can exceed the cap by measurement while
        # its bound passes. Both must pass; they fail closed when they disagree.
        bound <= rule.cap && excess <= rule.cap ||
            scope_error("FootprintTopology",
                        "the collar region encloses an un-etched island of " *
                        "$(length(island)) boundary sub-segments near $at whose excess beyond " *
                        "the collar $(rule.collar) is bounded by $bound (half its smallest " *
                        "width) and measured $excess (exact metal distance), above the " *
                        "admitted $(rule.cap)" *
                        (excess > bound ? " (the measurement exceeds the bound: the island " *
                                          "is bounded by a convex-corner miter kite, beyond " *
                                          "the round offset)" : ""))
        absorbed === nothing || push!(absorbed, Dict{String, Any}(
            "Vertices" => length(island), "Area" => 0.5 * abs(polygon_area2(island)),
            "Points" => [collect(p) for p in island],
            "MaximumExcess" => excess, "MaximumExcessPoint" => collect(at),
            "MaximumExcessBound" => bound, "ExcessCap" => rule.cap, "Collar" => rule.collar))
    end
    return absorbed
end

# The plan-view polygon of an exterior loop offset by `distance` for a loft: the
# miter polygon of offset_loop_points while it is simple (construction
# "MiterOffset"). When two Physical sides face each other across a dielectric gap
# narrower than twice the offset (the device 5-edge coupon: 6 um notches against
# the 6 um producer-default collar), the miter polygon folds back through the gap
# and self-intersects; the OCC face built from it acquires arbitrary seams that
# split the tube volumes of the facing edges ("tube tool ... has 3 volume
# descendants"). The collar region is then assembled as the union of
# collar_pieces clipped to the coupon box `box` = (lower, upper) and its outer
# boundary is traced (construction "CollarUnion"); the two constructions describe
# the same region wherever the miter polygon is simple. Inward offsets
# (`distance` > 0, sloped-sidewall pullbacks) have no union form and fail closed.
# A curved loop (circular_arc_runs) is offset exactly (curved_offset_loop /
# curved_collar_pieces); its bridged offset polygon is never a simple offset. The
# union's un-etched islands fail closed unless `island_rule` (collar_island_rule)
# admits them; `absorbed` collects their records.
function collar_loop_points(loop, distance, box, tolerance; island_rule=nothing,
                            absorbed=nothing)
    runs = abs(distance) <= tolerance ? NamedTuple[] : loop_arc_runs(loop, tolerance)
    offset = isempty(runs) ? nothing : curved_offset_loop(loop, distance, runs, tolerance)
    miter = offset === nothing ? offset_loop_points(loop, distance, tolerance) : offset.points
    (offset === nothing || !offset.bridged) && polygon_is_simple(miter, tolerance) &&
        return miter, "MiterOffset"
    distance < 0.0 ||
        error("Inward offset $distance of the plan-view loop of conductor $(loop.conductor) " *
              "on plane $(loop.plane) self-intersects")
    box === nothing && error("The collar union of a self-intersecting offset needs the coupon box")
    pieces = offset === nothing ? collar_pieces(loop, distance, miter, box[1], box[2], tolerance) :
             curved_collar_pieces(loop, distance, offset, box[1], box[2], tolerance)
    return polygon_union_boundary(pieces, tolerance; island_rule=island_rule,
                                  absorbed=absorbed), "CollarUnion"
end

# Two consecutive etch-footprint edges are one edge when every vertex between
# their outer endpoints lies within this fraction of the merged edge's length from
# the merged edge. It is the same dimensionless roundoff-scale bound as
# edge_volume_metric.COPLANAR_TOLERANCE (1e-6, the sine of the dihedral between
# the two trench-wall faces the edges would create): a merged vertex never leaves
# a near-coplanar sliver face for the metric stage to see, and a kept vertex
# always bends the wall by more than the tolerance, so it is a genuine facet.
const FOOTPRINT_COLLINEAR_TOLERANCE = 1.0e-6

function point_segment_distance_2d(point, first, second)
    direction = (second[1] - first[1], second[2] - first[2])
    span = direction[1]^2 + direction[2]^2
    span > 0.0 || return hypot(point[1] - first[1], point[2] - first[2])
    parameter = clamp(((point[1] - first[1]) * direction[1] +
                       (point[2] - first[2]) * direction[2]) / span, 0.0, 1.0)
    return hypot(point[1] - first[1] - parameter * direction[1],
                 point[2] - first[2] - parameter * direction[2])
end

# Merge consecutive collinear edges of a closed plan-view footprint polygon before
# any CAD face is created from it: one CAD face per genuine facet. A vertex is
# removed when it, and every original vertex already merged into its two edges,
# deviates from the merged edge by at most `tolerance` times that edge's length;
# the vertex of smallest relative deviation is removed first until none qualifies.
# Returns the simplified points and the record (removed 1-based original vertex
# indices, maximum deviation and its local scale). Fails closed if the result is
# not a polygon or exceeds the tolerance.
function simplify_footprint_polygon(points, tolerance)
    length(points) >= 3 || error("Footprint polygon needs at least three vertices")
    tolerance > 0.0 || error("Footprint collinearity tolerance must be positive")
    kept = collect(eachindex(points))
    removed = Int[]
    maximum_deviation = 0.0
    maximum_scale = 0.0
    maximum_relative = 0.0
    function merged_span(k)
        # Original vertices strictly between the kept neighbours of kept[k].
        before = kept[mod1(k - 1, length(kept))]
        after = kept[mod1(k + 1, length(kept))]
        span = Int[]
        index = mod1(before + 1, length(points))
        while index != after
            push!(span, index)
            index = mod1(index + 1, length(points))
        end
        return before, after, span
    end
    while length(kept) > 3
        best = 0
        best_relative = Inf
        best_deviation = 0.0
        best_scale = 0.0
        for k in eachindex(kept)
            before, after, span = merged_span(k)
            scale = hypot(points[after][1] - points[before][1],
                          points[after][2] - points[before][2])
            scale > 0.0 || continue
            deviation = maximum(point_segment_distance_2d(points[index], points[before],
                                                          points[after]) for index in span)
            relative = deviation / scale
            if relative <= tolerance && relative < best_relative
                best, best_relative, best_deviation, best_scale = k, relative, deviation, scale
            end
        end
        best == 0 && break
        push!(removed, kept[best])
        deleteat!(kept, best)
        if best_relative >= maximum_relative
            maximum_deviation, maximum_scale, maximum_relative =
                best_deviation, best_scale, best_relative
        end
    end
    sort!(removed)
    simplified = [points[index] for index in kept]
    maximum_relative <= tolerance ||
        error("Footprint simplification exceeded its collinearity tolerance")
    abs(sum(cross2d(simplified[i], simplified[mod1(i + 1, length(simplified))])
            for i in eachindex(simplified))) > 0.0 ||
        error("Footprint simplification produced a degenerate polygon")
    record = Dict{String, Any}(
        "OriginalVertices" => length(points), "Vertices" => length(simplified),
        "RemovedVertexCount" => length(removed), "RemovedVertexIndices" => removed,
        "MaximumDeviation" => maximum_deviation,
        "MaximumDeviationLocalScale" => maximum_scale,
        "MaximumRelativeDeviation" => maximum_relative, "Tolerance" => tolerance)
    return simplified, record
end

# `construction` names how the polygon was built before simplification:
# "MiterOffset" / "CollarUnion" (collar_loop_points), "HoleOffset"
# (offset_hole_points) or "EdgeStrip" (loft_strip); a CollarUnion polygon records
# the un-etched islands absorbed by the island rule (absorb_collar_islands!).
function footprint_record(conductor, plane, hole, points, record, construction;
                          absorbed_islands=nothing)
    footprint = Dict{String, Any}(
        "Conductor" => conductor, "Plane" => plane, "Hole" => hole,
        "Points" => [collect(point) for point in points], "Simplification" => record,
        "Construction" => construction)
    construction == "CollarUnion" &&
        (footprint["AbsorbedIslands"] = absorbed_islands === nothing ? Dict{String, Any}[] :
                                        absorbed_islands)
    return footprint
end

# Simplify the bottom and top polygons of a footprint loft; when `footprint` is a
# vector, the loft must be prismatic (one polygon) and the polygon is recorded.
function simplified_loft_polygons(bottom_points, top_points, footprint, conductor, plane, hole,
                                  construction; absorbed_islands=nothing)
    bottom_points, bottom_record =
        simplify_footprint_polygon(bottom_points, FOOTPRINT_COLLINEAR_TOLERANCE)
    top_points, _ = simplify_footprint_polygon(top_points, FOOTPRINT_COLLINEAR_TOLERANCE)
    if footprint !== nothing
        bottom_points == top_points ||
            error("Footprint recording requires a prismatic (vertical-wall) loft")
        push!(footprint, footprint_record(conductor, plane, hole, bottom_points, bottom_record,
                                          construction; absorbed_islands=absorbed_islands))
    end
    return bottom_points, top_points
end

# A prismatic loft (equal bottom and top polygons) whose wire carries circular arcs is
# EXTRUDED, so its walls are analytic OCC cylinders: the fragment then recognises the arc
# tubes' revolved sidewall rays (the same cylinders) as the same faces and merges them
# (design 1.2 (3)); a ThruSections loft turns the arcs into BSpline walls that OCC does not
# merge with a cylinder (a dangling coincident face Gmsh cannot mesh). A straight polygon
# (and a non-prismatic loft) keeps the ThruSections construction bitwise.
# The TAGGED runs of a loop for the wire of `points` when `points` are the loop's own vertices
# (an offset-0 loft, a sheet face: the same vertices, so the run indices hold) and the loop
# carries arc tags; nothing otherwise (polygon_wire then fits the points itself). One source of
# truth for the split angles of every surface coincident with an arc tube: the tagged runs of
# the metal loop are read by the tubes, the metal loft, the trench wall under the metal edge
# and the thin sheet face alike (the untagged fit of the 1e-9 R quantised chord vertices
# carries ~1e-7-degree noise, which flips the part count of an exactly 90-degree arc).
function loop_wire_runs(loop, points, tolerance)
    loop_has_arcs(loop) && points == loop.points || return nothing
    return tagged_arc_runs(loop, tolerance)
end

function loft_polygon(occ, bottom_points, top_points, z0, z1; runs=nothing)
    if bottom_points == top_points
        tolerance = 1.0e-9 * max(
            maximum(point[1] for point in bottom_points) - minimum(point[1] for point in bottom_points),
            maximum(point[2] for point in bottom_points) - minimum(point[2] for point in bottom_points),
            1.0)
        if runs !== nothing || !isempty(circular_arc_runs(bottom_points, tolerance))
            face = occ.addPlaneSurface([polygon_wire(occ, bottom_points, z0; runs=runs)])
            entities = occ.extrude([(2, face)], 0.0, 0.0, z1 - z0)
            volumes = [(dim, tag) for (dim, tag) in entities if dim == 3]
            length(volumes) == 1 || error("Plan-view mask extrusion produced $(length(volumes)) volumes")
            return volumes
        end
    end
    bottom = polygon_wire(occ, bottom_points, z0)
    top = polygon_wire(occ, top_points, z1)
    entities = occ.addThruSections([bottom, top], -1, true, false, -1, "C0")
    volumes = [(dim, tag) for (dim, tag) in entities if dim == 3]
    isempty(volumes) && error("Plan-view mask loft produced no volume")
    return volumes
end

# A hole offset by `distance`: grown (`distance` >= 0, the metal side pulled back)
# by offset_loop_points; shrunk (`distance` < 0: the collar of its metal sides)
# by clipping the hole polygon to the inward offset of every Physical side. A
# straight side clips by its shifted half-plane; a circular arc (convex, like the
# hole) by the chords of its concentric offset arc of radius r - |distance| (its
# chord vertices are exact points of that circle, which polygon_wire rebuilds), and
# an arc smaller than the offset clips nothing: its neighbours' half-planes imply
# every tangent half-plane of a collapsed arc. The result may be empty.
function offset_hole_points(loop,distance,tolerance)
    distance>=-tolerance && return offset_loop_points(loop,distance,tolerance)
    points=copy(loop.points)
    orientation=loop_orientation(points)
    # Convex hole erosion is an intersection of inward-offset half-planes. It
    # may vanish: intersecting offset lines alone would incorrectly reopen a
    # reflected polygon after collapse. Nonconvex topology changes fail closed.
    for i in eachindex(points)
        a,b,c=points[i],points[mod1(i+1,length(points))],points[mod1(i+2,length(points))]
        orientation*cross2d((b[1]-a[1],b[2]-a[2]),(c[1]-b[1],c[2]-b[2]))>=-tolerance ||
            error("Shrinking a nonconvex fabrication hole requires topology-aware offset support")
    end
    runs=loop_arc_runs(loop,tolerance)
    sides=Tuple{NTuple{2,Float64},NTuple{2,Float64},Float64}[]
    arc_edges=falses(length(points))
    for run in runs
        classes=unique(loop.classes[run.edge_indices])
        length(classes)==1 ||
            error("A circular arc of the fabrication hole of conductor $(loop.conductor) mixes " *
                  "Physical and Continuation chords")
        arc_edges[run.edge_indices].=true
        classes[1]=="Physical" || continue
        radius=run.radius+run.orientation*orientation*distance
        if radius<=tolerance
            # A collapsed whole-circle hole vanishes; a collapsed arc is implied by its
            # neighbours' half-planes.
            length(run.edge_indices)==length(points) && return NTuple{2,Float64}[]
            continue
        end
        scaled(p)=(run.center[1]+radius*(p[1]-run.center[1])/hypot(p[1]-run.center[1],p[2]-run.center[2]),
                   run.center[2]+radius*(p[2]-run.center[2])/hypot(p[1]-run.center[1],p[2]-run.center[2]))
        for (first,second) in zip(run.point_indices[1:(end-1)],run.point_indices[2:end])
            push!(sides,(scaled(points[first]),scaled(points[second]),0.0))
        end
    end
    for i in eachindex(points)
        loop.classes[i]=="Physical" && !arc_edges[i] || continue
        push!(sides,(points[i],points[mod1(i+1,length(points))],distance))
    end
    clipped=copy(points)
    for (a,b,shift) in sides
        direction=(b[1]-a[1],b[2]-a[2]);edge_length=hypot(direction...)
        normal=(-orientation*direction[2]/edge_length,orientation*direction[1]/edge_length)
        signed(p)=normal[1]*(p[1]-a[1])+normal[2]*(p[2]-a[2])+shift
        result=NTuple{2,Float64}[]
        isempty(clipped) && return result
        for j in eachindex(clipped)
            p,q=clipped[j],clipped[mod1(j+1,Base.length(clipped))]
            dp,dq=signed(p),signed(q)
            dp>=-tolerance && push!(result,p)
            if (dp>=-tolerance)!=(dq>=-tolerance)
                t=dp/(dp-dq)
                push!(result,(p[1]+t*(q[1]-p[1]),p[2]+t*(q[2]-p[2])))
            end
        end
        clipped=result
    end
    cleaned=NTuple{2,Float64}[]
    for p in clipped
        (isempty(cleaned) || hypot(p[1]-cleaned[end][1],p[2]-cleaned[end][2])>tolerance) && push!(cleaned,p)
    end
    if Base.length(cleaned)>1 && hypot(cleaned[1][1]-cleaned[end][1],cleaned[1][2]-cleaned[end][2])<=tolerance
        pop!(cleaned)
    end
    Base.length(cleaned)>=3 || return NTuple{2,Float64}[]
    area2=sum(cross2d(cleaned[i],cleaned[mod1(i+1,Base.length(cleaned))]) for i in eachindex(cleaned))
    return abs(area2)>tolerance^2 ? cleaned : NTuple{2,Float64}[]
end

# `simplify` merges collinear polygon edges (etch footprints: one CAD face per
# genuine facet); `footprint` additionally records every simplified polygon; `box`
# = (lower, upper) is the coupon box the collar union of a self-intersecting
# exterior offset is clipped to (collar_loop_points); `island_radius` (the coupon
# Radius) enables the collar island rule (collar_island_rule) of that union.
function loft_mask_offsets(occ, loops, z0, z1, bottom_offset, top_offset, tolerance;
                           simplify=false, footprint=nothing, box=nothing, island_radius=nothing)
    footprint === nothing || simplify || error("Footprint recording requires simplification")
    outers = [loop for loop in loops if !loop.hole]
    holes = [loop for loop in loops if loop.hole]
    isempty(outers) && error("Plan-view mask has no exterior loop")
    result = Tuple{Int32, Int32}[]
    hole_owners=zeros(Int,length(holes))
    for outer in outers
        island_rule(offset) = island_radius === nothing || box === nothing ? nothing :
                              collar_island_rule(outer, offset, box, island_radius, tolerance)
        bottom_islands = Dict{String, Any}[]
        top_islands = Dict{String, Any}[]
        bottom_points, bottom_construction =
            collar_loop_points(outer, bottom_offset, box, tolerance;
                               island_rule=island_rule(bottom_offset), absorbed=bottom_islands)
        top_points, top_construction =
            collar_loop_points(outer, top_offset, box, tolerance;
                               island_rule=island_rule(top_offset), absorbed=top_islands)
        if simplify
            bottom_construction == top_construction ||
                error("Footprint loft mixes the $bottom_construction and $top_construction " *
                      "constructions")
            bottom_points, top_points = simplified_loft_polygons(
                bottom_points, top_points, footprint, outer.conductor, z0, false,
                bottom_construction; absorbed_islands=bottom_islands)
        end
        volume = loft_polygon(occ, bottom_points, top_points, z0, z1;
                              runs=loop_wire_runs(outer, bottom_points, tolerance))
        cutters = Tuple{Int32, Int32}[]
        for (index,hole) in enumerate(holes)
            hole.conductor==outer.conductor && abs(hole.plane-outer.plane)<=tolerance || continue
            point_in_polygon(hole.points[1], outer.points, tolerance) || continue
            hole_owners[index]+=1
            bottom_hole=offset_hole_points(hole,bottom_offset,tolerance)
            top_hole=offset_hole_points(hole,top_offset,tolerance)
            isempty(bottom_hole) && isempty(top_hole) && continue
            isempty(bottom_hole)==isempty(top_hole) ||
                error("Fabrication hole collapses across loft height; unsupported topology change")
            if simplify
                bottom_hole, top_hole = simplified_loft_polygons(
                    bottom_hole, top_hole, footprint, hole.conductor, z0, true, "HoleOffset")
            end
            append!(
                cutters,
                loft_polygon(occ,bottom_hole,top_hole,z0,z1; runs=loop_wire_runs(hole, bottom_hole, tolerance))
            )
        end
        if !isempty(cutters)
            volume, _ = occ.cut(volume, cutters)
            volume = [(dim, tag) for (dim, tag) in volume if dim == 3]
        end
        append!(result, volume)
    end
    all(==(1),hole_owners) || error("Every fabrication hole must belong to exactly one exterior conductor/layer loop")
    return fuse_all(occ, result)
end

function loft_mask(occ, loops, z0, z1, pullback, tolerance; simplify=false, footprint=nothing)
    return loft_mask_offsets(occ, loops, z0, z1, 0.0, pullback, tolerance;
                             simplify=simplify, footprint=footprint)
end

# The expanded collar polygons are the producer-default etch footprint: they are
# simplified and recorded like a device footprint. The retained metal mask is not.
function boundary_strips(occ, loops, radius, z0, z1, pullback, tolerance; footprint=nothing,
                         box=nothing)
    expanded_volumes = Tuple{Int32, Int32}[]
    retained_volumes = Tuple{Int32, Int32}[]
    width = 3radius
    for conductor in sort!(unique(loop.conductor for loop in loops))
        conductor_loops = [loop for loop in loops if loop.conductor == conductor]
        append!(expanded_volumes,
                loft_mask_offsets(occ, conductor_loops, z0, z1, -width, -width, tolerance;
                                  simplify=true, footprint=footprint, box=box,
                                  island_radius=radius))
        append!(retained_volumes,
                loft_mask_offsets(occ, conductor_loops, z0, z1, 0.0, -pullback, tolerance))
    end
    # A conductor's etch collar must not remove substrate beneath another conductor.
    # Subtract the complete retained mask after combining the expanded collars.
    expanded = fuse_all(occ, expanded_volumes)
    retained = fuse_all(occ, retained_volumes)
    strip, _ = occ.cut(expanded, retained)
    strip = [(dim, tag) for (dim, tag) in strip if dim == 3]
    isempty(strip) && error("Classified boundary strip produced no volume")
    return fuse_all(occ, strip)
end

function loft_strip(occ, edge, radius, side, z0, z1, pullback; footprint=nothing)
    # Producer-default per-edge trench strips are etch footprint polygons too
    # (recorded when `footprint` is given); strips are rectangles, so the
    # simplification is a no-op that keeps one rule for every lofted polygon.
    bottom_points, top_points = simplified_loft_polygons(
        collect(strip_points(edge, radius, side)),
        collect(strip_points(edge, radius, side, pullback)), footprint, edge.conductor, z0,
        false, "EdgeStrip")
    bottom = polygon_wire(occ, bottom_points, z0)
    top = polygon_wire(occ, top_points, z1)
    entities = occ.addThruSections([bottom, top], -1, true, false, -1, "C0")
    volumes = [(dim, tag) for (dim, tag) in entities if dim == 3]
    isempty(volumes) && error("Spatial strip loft produced no volume")
    return volumes
end

function extruded_strip(occ, edge, radius, side, z0, dz)
    surface = occ.addPlaneSurface([polygon_wire(occ, strip_points(edge, radius, side), z0)])
    return [
        (dim, tag) for (dim, tag) in occ.extrude([(2, surface)], 0.0, 0.0, dz) if dim == 3
    ]
end

function planar_mask_surfaces(occ, loops, tolerance)
    surfaces = Tuple{Int32,Int32}[]
    holes=[loop for loop in loops if loop.hole]
    owners=zeros(Int,length(holes))
    for outer in loops
        outer.hole && continue
        wires = [polygon_wire(occ, outer.points, outer.plane; runs=loop_wire_runs(outer, outer.points, tolerance))]
        for (i,hole) in enumerate(holes)
            hole.conductor == outer.conductor &&
                abs(hole.plane - outer.plane) <= tolerance || continue
            point_in_polygon(hole.points[1], outer.points, tolerance) || continue
            owners[i]+=1
            push!(wires, polygon_wire(occ, hole.points, hole.plane; runs=loop_wire_runs(hole, hole.points, tolerance)))
        end
        push!(surfaces, (2, occ.addPlaneSurface(wires)))
    end
    all(==(1),owners) || error("Every mask hole must belong to exactly one exterior loop of its conductor/layer")
    isempty(surfaces) && error("Planar mask has no exterior surfaces")
    return surfaces
end

function mask_prism(occ, facets, z0, dz)
    volumes = Tuple{Int32, Int32}[]
    for facet in facets
        surface = occ.addPlaneSurface([polygon_wire(occ, facet.points, z0)])
        append!(
            volumes,
            [
                (dim, tag) for
                (dim, tag) in occ.extrude([(2, surface)], 0.0, 0.0, dz) if dim == 3
            ]
        )
    end
    return fuse_all(occ, volumes)
end

function apply_plan_view_mask(occ, volumes, facets, lower, upper)
    isempty(facets) && return volumes
    mask = mask_prism(occ, facets, lower[3], upper[3] - lower[3])
    isempty(mask) && error("Plan-view mask produced no volume")
    result, _ = occ.intersect(volumes, mask)
    result = [(dim, tag) for (dim, tag) in result if dim == 3]
    isempty(result) && error("Plan-view mask removed all edge-strip metal")
    return result
end

function fuse_all(occ, volumes)
    isempty(volumes) && return Tuple{Int32, Int32}[]
    result = [volumes[1]]
    for volume in volumes[2:end]
        result, _ = occ.fuse(result, [volume])
    end
    return result
end

function boundary_curves(volumes)
    gmsh.model.occ.synchronize()
    surfaces = [
        entity for
        entity in gmsh.model.getBoundary(volumes, false, false, false) if entity[1] == 2
    ]
    return unique(
        tag for
        (dim, tag) in gmsh.model.getBoundary(surfaces, false, false, false) if dim == 1
    )
end

function curve_lies_on_plane(curve,z,tolerance)
    lower,upper=gmsh.model.getParametrizationBounds(1,curve)
    parameters=[lower[1],(lower[1]+upper[1])/2,upper[1]]
    values=gmsh.model.getValue(1,curve,parameters)
    # OCC bounding boxes have scale-independent padding: test the curve itself.
    return all(abs(values[3i]-z)<=tolerance for i in 1:3)
end

function fillet_plane_edges(occ, volumes, radius, z, tolerance)
    radius <= 0.0 && return volumes
    curves = [curve for curve in boundary_curves(volumes) if
              curve_lies_on_plane(curve,z,tolerance)]
    isempty(curves) && error("Requested rounding found no curves on process plane $z")
    rounded = occ.fillet(Int32[tag for (dim, tag) in volumes if dim == 3], curves, [radius])
    result = [(dim, tag) for (dim, tag) in rounded if dim == 3]
    isempty(result) && error("Requested rounding produced no solid")
    return result
end

function point_segment_distance(point, first, second)
    direction = (second[1] - first[1], second[2] - first[2])
    length_squared = direction[1]^2 + direction[2]^2
    length_squared > 0.0 || return hypot(point[1] - first[1], point[2] - first[2])
    coordinate = clamp(
        ((point[1] - first[1]) * direction[1] + (point[2] - first[2]) * direction[2]) /
        length_squared,
        0.0,
        1.0
    )
    closest = (first[1] + coordinate * direction[1], first[2] + coordinate * direction[2])
    return hypot(point[1] - closest[1], point[2] - closest[2])
end

function physical_segments(loops, offset, tolerance)
    primitives = NamedTuple[]
    for loop in loops
        points = offset_loop_points(loop, offset, tolerance)
        covered = falses(length(points))
        for run in circular_arc_runs(points, tolerance)
            all(loop.classes[index] == "Physical" for index in run.edge_indices) || continue
            push!(
                primitives,
                (
                    kind=:arc,
                    center=run.center,
                    radius=run.radius,
                    first=points[first(run.point_indices)],
                    last=points[last(run.point_indices)],
                    orientation=run.orientation,
                    angle=run.angle
                )
            )
            covered[run.edge_indices] .= true
        end
        for index in eachindex(points)
            loop.classes[index] == "Physical" && !covered[index] || continue
            push!(
                primitives,
                (
                    kind=:line,
                    first=points[index],
                    last=points[mod1(index + 1, length(points))]
                )
            )
        end
    end
    return primitives
end

function directed_angle(first, second, orientation)
    angle =
        orientation *
        atan(cross2d(first, second), first[1] * second[1] + first[2] * second[2])
    return mod(angle, 2pi)
end

function point_primitive_distance(point, primitive, tolerance)
    primitive.kind == :line &&
        return point_segment_distance(point, primitive.first, primitive.last)
    radial = (point[1] - primitive.center[1], point[2] - primitive.center[2])
    start =
        (primitive.first[1] - primitive.center[1], primitive.first[2] - primitive.center[2])
    angle = directed_angle(start, radial, primitive.orientation)
    angle_tolerance = tolerance / max(primitive.radius, tolerance)
    if angle <= primitive.angle + angle_tolerance
        return abs(hypot(radial...) - primitive.radius)
    end
    return min(
        hypot(point[1] - primitive.first[1], point[2] - primitive.first[2]),
        hypot(point[1] - primitive.last[1], point[2] - primitive.last[2])
    )
end

function point_on_curve(curve)
    lower, upper = gmsh.model.getParametrizationBounds(1, curve)
    length(lower) == 1 && length(upper) == 1 ||
        error("Unexpected curve parametrization for curve $curve")
    value = gmsh.model.getValue(1, curve, [(lower[1] + upper[1]) / 2])
    return (value[1], value[2])
end

function fillet_physical_edges(occ, volumes, radius, z, primitives, tolerance)
    radius <= 0.0 && return volumes
    gmsh.model.occ.synchronize()
    curves = Int32[]
    for curve in boundary_curves(volumes)
        curve_lies_on_plane(curve,z,tolerance) || continue
        point = point_on_curve(curve)
        any(
            point_primitive_distance(point, primitive, tolerance) <= 10tolerance for
            primitive in primitives
        ) && push!(curves, curve)
    end
    isempty(curves) && error("Requested physical-edge rounding found no curves on plane $z")
    rounded = occ.fillet(Int32[tag for (dim, tag) in volumes if dim == 3], curves, [radius])
    result = [(dim, tag) for (dim, tag) in rounded if dim == 3]
    isempty(result) && error("Requested physical-edge rounding produced no solid")
    return result
end

# The SupportBox [x0, y0, x1, y1] of the single model of a process library (a device-plan
# coupon, decision 282: the signature's grown box in the canonical frame = the mesh frame),
# or nothing for a legacy model.
function read_support_box(process_library_path)
    process_library_path === nothing && return nothing
    library = parse_json(read(process_library_path, String))
    models = get(library, "Models", Any[])
    length(models) == 1 || error("The process library must bind exactly one model")
    haskey(models[1], "SupportBox") || return nothing
    box = ntuple(d -> Float64(models[1]["SupportBox"][d]), 4)
    all(isfinite, box) && box[3] > box[1] && box[4] > box[2] ||
        error("The process library's SupportBox is not a box")
    return box
end

function coupon_bounds(edges, radius, metal_thickness, overetch; support_box=nothing)
    if support_box !== nothing
        # The z extent from the rows (the plane, the metal and the trench), the plane box
        # the signature's.
        lower, upper = row_coupon_bounds(edges, radius, metal_thickness, overetch)
        return (support_box[1], support_box[2], lower[3]),
               (support_box[3], support_box[4], upper[3])
    end
    any(get(edge, :context, false) for edge in edges) &&
        error("Device-plan coupon rows (Context) need the process library's SupportBox")
    return row_coupon_bounds(edges, radius, metal_thickness, overetch)
end

function row_coupon_bounds(edges, radius, metal_thickness, overetch)
    points = NTuple{3, Float64}[]
    for edge in edges
        first, second = extended_interval(edge, radius)
        for coordinate in (first, second)
            boundary = add(edge.point, scale(coordinate, edge.tangent))
            for side in (-1.0, 1.0)
                # The final radius of box padding places the transverse matching
                # boundary 2R from the physical edge. The metal strip extends 3R
                # inward, so it is truncated without an artificial back edge.
                push!(points, add(boundary, scale(side * radius, edge.gap)))
            end
        end
    end
    lower = ntuple(index -> minimum(point[index] for point in points) - radius, 3)
    upper = ntuple(index -> maximum(point[index] for point in points) + radius, 3)
    # Vertical padding per process layer sign: the trench (Overetch) lies on the
    # substrate side of the plane, -Nz, and the metal (MetalThickness) on +Nz.
    lower = (
        lower[1],
        lower[2],
        min(lower[3], minimum(edge.point[3] - radius -
                              (edge.normal_sign > 0 ? overetch : metal_thickness)
                              for edge in edges))
    )
    upper = (
        upper[1],
        upper[2],
        max(upper[3], maximum(edge.point[3] + radius +
                              (edge.normal_sign > 0 ? metal_thickness : overetch)
                              for edge in edges))
    )
    return lower, upper
end

function layer_groups(edges, tolerance)
    groups = Dict{Tuple{Int, Int}, Vector{eltype(edges)}}()
    for edge in edges
        plane = round(Int, edge.point[3] / tolerance)
        key = (plane, Int(edge.normal_sign))
        push!(get!(groups, key, eltype(edges)[]), edge)
    end
    layers = [
        (
            plane = sum(edge.point[3] for edge in group) / length(group),
            sign  = key[2],
            edges = group
        ) for (key, group) in groups
    ]
    sort!(layers, by=layer -> layer.plane)
    for first in eachindex(layers), second = (first + 1):length(layers)
        a = layers[first]
        b = layers[second]
        overlap =
            (a.sign > 0 && b.plane < a.plane) ||
            (a.sign < 0 && b.plane > a.plane) ||
            (b.sign > 0 && a.plane < b.plane) ||
            (b.sign < 0 && a.plane > b.plane)
        overlap && error("Spatial coupon substrate half-spaces overlap")
    end
    return layers
end

# The process layer whose band [plane - Nz x Overetch, plane + Nz x MetalThickness]
# contains the z-range [zmin, zmax] of a fabricated surface (OCC bounding boxes carry
# their own tolerance: the box tolerance is used); fail closed when none or several do
# (layer_groups rejects overlapping half-spaces, so bands are disjoint).
function surface_process_layer(layers, zmin, zmax, metal_thickness, overetch, tolerance)
    matches = [layer for layer in layers
               if min(layer.plane - layer.sign * overetch, layer.plane + layer.sign * metal_thickness) -
                  tolerance <= zmin &&
                  zmax <= max(layer.plane - layer.sign * overetch,
                              layer.plane + layer.sign * metal_thickness) + tolerance]
    length(matches) == 1 ||
        error("Fabricated surface spanning z in [$zmin, $zmax] lies in $(length(matches)) process " *
              "layer bands")
    return only(matches)
end

function on_outer_box(bounds, lower, upper, tolerance)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    return (
        (abs(xmin - lower[1]) < tolerance && abs(xmax - lower[1]) < tolerance) ||
        (abs(xmin - upper[1]) < tolerance && abs(xmax - upper[1]) < tolerance) ||
        (abs(ymin - lower[2]) < tolerance && abs(ymax - lower[2]) < tolerance) ||
        (abs(ymin - upper[2]) < tolerance && abs(ymax - upper[2]) < tolerance) ||
        (abs(zmin - lower[3]) < tolerance && abs(zmax - lower[3]) < tolerance) ||
        (abs(zmin - upper[3]) < tolerance && abs(zmax - upper[3]) < tolerance)
    )
end

function point_on_surface(tag)
    center = gmsh.model.occ.getCenterOfMass(2, tag)
    coordinate = collect(center)
    # The centre of mass of a CURVED face (the cylindrical wall of an arc side: a ThruSections
    # BSpline) lies off the surface by rho (1 - sinc(dtheta / 2)) although isInside projects
    # it onto the face and says yes; the centre is taken only when it lies ON the surface
    # (a planar face: bitwise as before), else the probes below find a true surface point.
    if gmsh.model.isInside(2, tag, coordinate) > 0
        projected = gmsh.model.getValue(2, tag, gmsh.model.getParametrization(2, tag, coordinate))
        scale = max(1.0, maximum(abs, coordinate))
        norm(collect(projected) .- coordinate) <= 1.0e-9 * scale && return center
    end

    # A trimmed annulus or a thin ribbon can miss every point of a uniform UV
    # grid. Probe inward from real boundary curves at a scale derived from the
    # face area/perimeter, testing both orientations instead of assuming winding.
    curves = [curve for (dim,curve) in gmsh.model.getBoundary([(2,tag)],false,false,false)
              if dim==1]
    perimeter = sum(gmsh.model.occ.getMass(1,curve) for curve in curves)
    area = gmsh.model.occ.getMass(2,tag)
    area>0 && perimeter>0 || error("Degenerate trimmed surface $tag")
    width = area/perimeter
    for curve in curves
        lo,hi = gmsh.model.getParametrizationBounds(1,curve)
        for fraction in (.5,.25,.75)
            parameter = [lo[1]+fraction*(hi[1]-lo[1])]
            boundary = gmsh.model.getValue(1,curve,parameter)
            tangent = gmsh.model.getDerivative(1,curve,parameter)
            norm(tangent)>0 || continue
            uv = gmsh.model.getParametrization(2,tag,boundary)
            normal = gmsh.model.getNormal(tag,uv)
            inward = cross(normal,tangent)
            norm(inward)>0 || continue
            inward/=norm(inward)
            for factor in (.25,.1,.5,.01), side in (-1.,1.)
                distance = factor*width
                candidate = boundary+side*distance*inward
                parameter2 = gmsh.model.getParametrization(2,tag,candidate)
                gmsh.model.isInside(2,tag,parameter2,true)>0 || continue
                value = gmsh.model.getValue(2,tag,parameter2)
                norm(value-candidate)<=0.5distance && norm(value-boundary)>=0.5distance || continue
                return Tuple(value)
            end
        end
    end

    lower, upper = gmsh.model.getParametrizationBounds(2, tag)
    length(lower) == 2 && length(upper) == 2 ||
        error("Unexpected surface parametrization for surface $tag")
    for samples in (5, 11, 21)
        for i = 1:samples, j = 1:samples
            parameter = [
                lower[1] + (i - 0.5) * (upper[1] - lower[1]) / samples,
                lower[2] + (j - 0.5) * (upper[2] - lower[2]) / samples
            ]
            if gmsh.model.isInside(2, tag, parameter, true) > 0
                value = gmsh.model.getValue(2, tag, parameter)
                return (value[1], value[2], value[3])
            end
        end
    end
    return error("Unable to find a point on trimmed surface $tag")
end

function segment_distance(edge, point, radius)
    first, second = extended_interval(edge, radius)
    delta = (point[1] - edge.point[1], point[2] - edge.point[2], point[3] - edge.point[3])
    tangent_norm2=sum(value^2 for value in edge.tangent)
    tangent_norm2>0 || error("Zero-length edge tangent")
    coordinate = clamp(sum(delta[i]*edge.tangent[i] for i in 1:3)/tangent_norm2,
                       first,second)
    closest = add(edge.point, scale(coordinate, edge.tangent))
    return sqrt(sum((point[index] - closest[index])^2 for index = 1:3))
end

function nearest_edge(edges, point, radius)
    isempty(edges) && error("Cannot assign ownership without edges")
    distances=[segment_distance(edge,point,radius) for edge in edges]
    best=minimum(distances)
    # At a bisector prefer a stable physical label, not input row order. Slot ties
    # have identical output ownership and require no coordinate-dependent rule.
    tied=[i for i in eachindex(edges) if distances[i]<=best+1e-12max(radius,1.)]
    index=tied[argmin((edges[i].conductor,edges[i].slot) for i in tied)]
    return edges[index]
end

const METAL_SLOT_STRIDE = 100
function metal_surface_attribute(base, slot, conductor)
    0 <= slot < 10 || error("Metal interface slot must lie in [0, 10)")
    1 <= conductor < METAL_SLOT_STRIDE ||
        error("Metal conductor label must lie in [1, $METAL_SLOT_STRIDE)")
    return base + METAL_SLOT_STRIDE * slot + conductor
end

function nearest_metal_edge(edges, facets, point, radius, tolerance; surface=nothing)
    candidates = if isempty(facets)
        edges
    else
        [
            edge for edge in edges if
            point_in_mask(facets, point, edge.conductor, edge.point[3], tolerance)
        ]
    end
    isempty(candidates) && error(
        "Unable to assign a metal surface to a conductor mask (surface point $point" *
        (surface === nothing ? "" :
         "; surface $(surface[1]) of type $(gmsh.model.getType(2, surface[1])), bounding box " *
         "$(surface[2]), adjacent volumes $(surface[3])") * ")")
    return nearest_edge(candidates, point, radius)
end

function point_in_metal(edge, point, radius, tolerance, facets)
    abs(point[3] - edge.point[3]) <= tolerance || return false
    # An exact mask is the conductor geometry, not merely a clipping filter for a
    # finite collection of edge strips. Interior metal can extend beyond those strips.
    if !isempty(facets)
        return point_in_mask(facets, point, edge.conductor, edge.point[3], tolerance)
    end
    delta = (point[1] - edge.point[1], point[2] - edge.point[2], point[3] - edge.point[3])
    longitudinal = delta[1] * edge.tangent[1] + delta[2] * edge.tangent[2]
    transverse = delta[1] * edge.gap[1] + delta[2] * edge.gap[2]
    first, second = extended_interval(edge, radius)
    return first - tolerance <= longitudinal <= second + tolerance &&
           -3radius - tolerance <= transverse <= tolerance &&
           (
               isempty(facets) ||
               point_in_mask(facets, point, edge.conductor, edge.point[3], tolerance)
           )
end

function surface_family(attribute::Int)
    # Ignore conductor/slot suffixes when deciding whether a curve is a physical process
    # feature. Attributes 5001 and 5101, for example, are both MS surfaces; the curve
    # separating their bookkeeping patches is not a material or geometric edge.
    return div(attribute, 1000)
end

function coplanar_surfaces(surfaces::Vector{Int32}, tolerance::Float64)
    length(surfaces) >= 2 || return false
    all(gmsh.model.getType(2, surface) == "Plane" for surface in surfaces) || return false
    # OCC bounding boxes have geometric padding, so their nominally zero width
    # need not be below a nanometre-scale tolerance. Test the surfaces themselves.
    origin = collect(gmsh.model.occ.getCenterOfMass(2, first(surfaces)))
    uv = gmsh.model.getParametrization(2, first(surfaces), origin)
    normal = gmsh.model.getNormal(first(surfaces), uv)
    for surface in surfaces
        center = collect(gmsh.model.occ.getCenterOfMass(2, surface))
        abs(dot(normal, center-origin)) <= tolerance || return false
        lo, hi = gmsh.model.getParametrizationBounds(2, surface)
        for a in (.2, .5, .8), b in (.2, .5, .8)
            parameter = [lo[1]+a*(hi[1]-lo[1]), lo[2]+b*(hi[2]-lo[2])]
            point = gmsh.model.getValue(2, surface, parameter)
            direction = gmsh.model.getNormal(surface, parameter)
            abs(dot(normal, point-origin)) <= tolerance || return false
            abs(dot(normal, direction)) >= 1-1e-10 || return false
        end
    end
    return true
end

# Validate before selecting any CAD/sizing subset, including the diagnostic `none`.
# This is local triangle validity, not complete-box coverage or FEM trace accuracy.
function trace_triangle_areas(triangles)
    isempty(triangles) && error("Empty trace geometry")
    areas = Dict{Int,Float64}()
    for (index, triangle) in triangles
        index isa Integer && index > 0 || error("Trace triangle indices must be positive integers")
        length(triangle) == 3 || error("Trace triangle needs three vertices")
        all(p -> length(p) == 3 && all(isfinite, p), triangle) ||
            error("Trace vertices must have three finite coordinates")
        a, b, c = triangle
        ab = ntuple(d -> b[d] - a[d], 3)
        ac = ntuple(d -> c[d] - a[d], 3)
        area2 = norm(cross(collect(ab), collect(ac)))
        isfinite(area2) && area2 > 0 || error("Degenerate trace triangle")
        areas[index] = area2
    end
    return areas
end

function read_matching_trace_triangles(path)
    data, header = readdlm(path, ',', header=true)
    names = vec(String.(header)); columns = Dict(name=>i for (i,name) in enumerate(names))
    length(columns) == length(names) || error("Duplicate matching trace columns")
    all(haskey(columns,key) for key in ("x","y","z","triangle")) ||
        error("Matching trace must have x,y,z,triangle columns")
    triangles = Dict{Int,Vector{NTuple{3,Float64}}}()
    for row in axes(data,1)
        p = ntuple(d->Float64(data[row,columns[("x","y","z")[d]]]),3)
        index = signature_integer(data[row,columns["triangle"]], "trace triangle")
        push!(get!(triangles,index,NTuple{3,Float64}[]),p)
    end
    trace_triangle_areas(triangles)
    return triangles
end

function matching_trace_lines(occ, path, lower, upper, tolerance; mode="all")
    mode in ("all","sides","levels","none") || error("Unknown matching trace constraint mode")
    source = read_matching_trace_triangles(path)
    mode == "none" && return Tuple{Int32,Int32}[]
    triangles = Dict(index => [ntuple(d->abs(p[d]-lower[d])<tolerance ? lower[d] :
        abs(p[d]-upper[d])<tolerance ? upper[d] : p[d],3) for p in triangle]
        for (index,triangle) in source)
    trace_triangle_areas(triangles)
    lines=Tuple{Int32,Int32}[]
    if mode == "levels"
        levels=sort!(unique(p[3] for tri in Base.values(triangles) for p in tri))
        for z in levels
            lower[3]+tolerance<z<upper[3]-tolerance || continue
            ring=[(lower[1],lower[2],z),(upper[1],lower[2],z),
                  (upper[1],upper[2],z),(lower[1],upper[2],z)]
            tags=[occ.addPoint(p...) for p in ring]
            for i in 1:4
                push!(lines,(1,occ.addLine(tags[i],tags[mod1(i+1,4)])))
            end
        end
        println("Matching-trace level constraints: $(length(lines)) segments")
        return lines
    end
    mode in ("all","sides") || error("Unknown matching trace constraint mode")
    points=Dict{NTuple{3,Float64},Int32}()
    seen=Set{Tuple{NTuple{3,Float64},NTuple{3,Float64}}}()
    skipped=0
    for tri in Base.values(triangles)
        length(tri)==3 || error("Matching trace triangle needs three vertices")
        for (a,b) in ((tri[1],tri[2]),(tri[2],tri[3]),(tri[3],tri[1]))
            key = isless(a,b) ? (a,b) : (b,a)
            key in seen && continue
            push!(seen,key)
            dimensions = mode == "sides" ? (1:2) : (1:3)
            if !any((a[d]==lower[d] && b[d]==lower[d]) ||
                    (a[d]==upper[d] && b[d]==upper[d]) for d in dimensions)
                skipped+=1;continue
            end
            sum((a[d]-b[d])^2 for d in 1:3)>tolerance^2 || continue
            pa=get!(points,a) do;occ.addPoint(a...);end
            pb=get!(points,b) do;occ.addPoint(b...);end
            push!(lines,(1,occ.addLine(pa,pb)))
        end
    end
    println("Matching-trace constraints: boundary segments=$(length(lines)), off-box segments=$skipped")
    return lines
end

# ---------------------------------------------------------------------------
# Trace-basis cut-surface sizing.
#
# The bound trace basis (basis-contract.json, trace-vertices.csv,
# trace-triangles.csv) is the closed box triangulation whose vertex hats are the
# response sources; its vertices are stored in the process-library frame and
# mapped to the mesh frame exactly as the campaign producer does (local z = the
# process normal, local x = the gap direction of the first edge). The size rule
# has one parameter, TraceBasisSizeRatio (dimensionless): on the cut surface the
# element size must not exceed the ratio times the minimum altitude of the basis
# triangle containing the point. The hat of a basis vertex varies linearly over
# the whole triangle with gradient 1 / (its altitude: the distance from the
# vertex to the opposite edge), so the Dirichlet datum's variation scale in a
# triangle is its minimum altitude = 2 x area / longest edge (the altitude of the
# vertex opposite the longest edge); the shortest edge equals it for right
# slivers but overstates it by longest / shortest x sin(angle) for a needle (a
# 0.05 x 8.7 um triangle with an 11 nm altitude: decision 43). The support is
# resolved where the hat varies only when the whole triangle is discretized at
# that scale. Away from the surface the size grows with the process-band grading
# slope up to lc_far, so far from narrow hats the far size stays. It composes
# with the scalar background field through Gmsh's size callback, which needs a
# named function (closure trampolines are unsupported on aarch64), hence the
# module-level state.
const TRACE_BASIS_TRIANGLES = NTuple{9, Float64}[]  # narrow triangles (requested < lc_far)
const TRACE_BASIS_SIZES = Float64[]                 # ratio x minimum altitude of each
const TRACE_BASIS_APEXES = NTuple{3, Float64}[]     # endpoints of the shortest edge of each
# Report-only classification of a basis triangle as a needle: its minimum altitude
# below this fraction of its shortest edge (the shortest-edge proxy overstated the
# hat scale by more than 1 / 0.6; altitude / shortest edge is 1 for right slivers,
# 0.866 equilateral, 0.707 right isosceles). Not a mesh parameter: it never enters a
# size. Mirrors trace_basis.NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE (the single documented
# definition; the stage contract requires the census value to equal it).
const NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE = 0.6
const TRACE_BASIS_SIZE_MEASURE =
    "minimum altitude of the basis triangle (2 x area / longest edge = 1 / the largest " *
    "gradient of its three vertex hats); the shortest edge equals it for right slivers " *
    "(where TraceBasisSizeRatio 0.5 was calibrated, physics-11) and overstates it for " *
    "needles - decision 43"
const TRACE_BASIS_NEEDLE_RULE =
    "a basis triangle whose minimum altitude is below NeedleAltitudeOverShortestEdge x " *
    "its shortest edge (report-only: needles are a property of the reference basis " *
    "triangulation, not of the coupon mesh)"
const TRACE_BASIS_FAR = Ref(Inf)
const TRACE_BASIS_SLOPE = Ref(0.0)
const TRACE_BASIS_CALLBACK_ROOTS = Any[]

function process_frame(library, model_name)
    library isa AbstractDict && library["Models"] isa AbstractVector ||
        error("Process library lacks Models")
    models = [model for model in library["Models"] if get(model, "Name", nothing) == model_name]
    length(models) == 1 || error("Process library must contain the trace basis model exactly once")
    model = models[1]
    entries = get(model, "Topology", nothing) == "SpatialEdgeCluster" ? get(model, "Edges", nothing) :
              get(model, "Arms", nothing)
    entries isa AbstractVector && !isempty(entries) ||
        error("Process library model has no complete edge geometry")
    # A device-plan coupon (spatial-support contract v3, decision 282) is built in its
    # signature's canonical frame: the mesh frame is the identity.
    haskey(model, "SupportBox") && return Matrix{Float64}(I, 3, 3)
    normal = collect(Float64, entries[1]["ProcessNormal"])
    gap = collect(Float64, entries[1]["GapDirection"])
    length(normal) == 3 && length(gap) == 3 && all(isfinite, normal) && all(isfinite, gap) &&
        norm(normal) > 0 && norm(gap) > 0 || error("Process library edge frame vectors are invalid")
    normal ./= norm(normal); gap ./= norm(gap)
    abs(dot(normal, gap)) <= 1.0e-9 || error("ProcessNormal and GapDirection must be orthogonal")
    axis_y = cross(normal, gap); axis_y ./= norm(axis_y)
    return vcat(gap', axis_y', normal')
end

function read_trace_basis(contract_path, vertices_path, triangles_path, library_path;
                          tolerance=1.0e-8)
    contract = parse_json(read(contract_path, String))
    contract isa AbstractDict && get(contract, "Version", nothing) == 1 &&
        get(contract, "Model", nothing) isa AbstractString && get(contract, "Geometry", nothing) isa AbstractDict ||
        error("Unsupported trace basis contract")
    geometry = contract["Geometry"]
    get(geometry, "OffBoxTriangles", nothing) == 0 && get(geometry, "ClosedOrientedSurface", nothing) === true ||
        error("Unsupported trace basis contract")
    lower = ntuple(d -> Float64(geometry["Lower"][d]), 3)
    upper = ntuple(d -> Float64(geometry["Upper"][d]), 3)
    all(isfinite, lower) && all(isfinite, upper) && all(upper[d] > lower[d] for d in 1:3) ||
        error("Trace basis contract box is invalid")
    frame = process_frame(parse_json(read(library_path, String)), contract["Model"])
    vertex_data, vertex_header = readdlm(vertices_path, ',', header=true)
    vertex_columns = Dict(String(name) => i for (i, name) in enumerate(vec(vertex_header)))
    all(haskey(vertex_columns, key) for key in ("vertex", "x", "y", "z", "basis", "conductor")) ||
        error("Trace vertices must have vertex,x,y,z,basis,conductor columns")
    triangle_data, triangle_header = readdlm(triangles_path, ',', header=true)
    triangle_columns = Dict(String(name) => i for (i, name) in enumerate(vec(triangle_header)))
    all(haskey(triangle_columns, key) for key in ("triangle", "vertex_i", "vertex_j", "vertex_k")) ||
        error("Trace triangles must have triangle,vertex_i,vertex_j,vertex_k columns")
    vertex_count = size(vertex_data, 1); triangle_count = size(triangle_data, 1)
    [signature_integer(vertex_data[i, vertex_columns["vertex"]], "trace vertex") for i in 1:vertex_count] ==
        collect(1:vertex_count) || error("Trace vertices must be numbered contiguously from one")
    [signature_integer(triangle_data[i, triangle_columns["triangle"]], "trace triangle") for i in 1:triangle_count] ==
        collect(1:triangle_count) || error("Trace triangles must be numbered contiguously from one")
    vertex_count == geometry["Vertices"] && triangle_count == geometry["Triangles"] ||
        error("Trace basis files differ from the contract geometry counts")
    scale = maximum(upper[d] - lower[d] for d in 1:3)
    points = NTuple{3, Float64}[]
    for i in 1:vertex_count
        canonical = ntuple(d -> Float64(vertex_data[i, vertex_columns[("x", "y", "z")[d]]]), 3)
        all(isfinite, canonical) || error("Trace vertex coordinates must be finite")
        point = ntuple(r -> sum(frame[r, d] * canonical[d] for d in 1:3), 3)
        on_face = any(abs(point[d] - lower[d]) <= tolerance * scale ||
                      abs(point[d] - upper[d]) <= tolerance * scale for d in 1:3)
        inside = all(lower[d] - tolerance * scale <= point[d] <= upper[d] + tolerance * scale for d in 1:3)
        on_face && inside || error("Trace basis vertices are not on the contract box surface")
        push!(points, point)
    end
    triangles = NTuple{3, Int}[]
    for i in 1:triangle_count
        triangle = ntuple(k -> signature_integer(
            triangle_data[i, triangle_columns[("vertex_i", "vertex_j", "vertex_k")[k]]], "trace vertex index"), 3)
        all(1 <= v <= vertex_count for v in triangle) && length(unique(triangle)) == 3 ||
            error("Trace triangle connectivity is invalid")
        a, b, c = (points[v] for v in triangle)
        norm(cross(collect(b .- a), collect(c .- a))) > 0 || error("Degenerate trace basis triangle")
        any(all(abs(points[v][d] - lower[d]) <= tolerance * scale for v in triangle) ||
            all(abs(points[v][d] - upper[d]) <= tolerance * scale for v in triangle) for d in 1:3) ||
            error("Trace basis triangle is not in one box face plane")
        push!(triangles, triangle)
    end
    hashes = Dict{String, Any}(
        "BasisContract" => bytes2hex(sha256(read(contract_path))),
        "TraceVertices" => bytes2hex(sha256(read(vertices_path))),
        "TraceTriangles" => bytes2hex(sha256(read(triangles_path))),
        "ProcessLibrary" => bytes2hex(sha256(read(library_path))))
    return (; points, triangles, lower, upper, frame, model=String(contract["Model"]), hashes)
end

trace_basis_edge_lengths(points, triangle) = ntuple(
    k -> norm(collect(points[triangle[mod1(k + 1, 3)]] .- points[triangle[k]])), 3)

# Minimum altitude of a basis triangle: 2 x area / longest edge, the smallest
# distance from a vertex to its opposite edge = 1 / the largest hat gradient.
function trace_basis_minimum_altitude(points, triangle, lengths)
    a = collect(points[triangle[1]]); b = collect(points[triangle[2]]); c = collect(points[triangle[3]])
    doubled_area = norm(cross(b .- a, c .- a))
    doubled_area > 0.0 || error("Degenerate trace basis triangle")
    return doubled_area / maximum(lengths)
end

# Euclidean distance from a point to the triangle (a, b, c): the plane projection
# when it falls inside, otherwise the nearest edge.
function point_triangle_distance(p, a, b, c)
    ab = b .- a; ac = c .- a; d = p .- a
    daa = dot(ab, ab); dab = dot(ab, ac); dbb = dot(ac, ac)
    dpa = dot(d, ab); dpb = dot(d, ac)
    denominator = daa * dbb - dab * dab
    v = (dbb * dpa - dab * dpb) / denominator
    w = (daa * dpb - dab * dpa) / denominator
    if v >= 0.0 && w >= 0.0 && v + w <= 1.0
        n = cross(collect(ab), collect(ac))
        return abs(dot(d, n)) / norm(n)
    end
    best = Inf
    for (first, last) in ((a, b), (b, c), (c, a))
        vector = last .- first; delta = p .- first
        t = clamp(dot(delta, vector) / dot(vector, vector), 0.0, 1.0)
        best = min(best, norm(collect(delta .- t .* vector)))
    end
    return best
end

# Size prescribed by the narrow basis triangles at (x, y, z), never above `lc`.
function trace_basis_size(x, y, z, lc)
    size = min(lc, TRACE_BASIS_FAR[])
    slope = TRACE_BASIS_SLOPE[]
    @inbounds for index in eachindex(TRACE_BASIS_TRIANGLES)
        requested = TRACE_BASIS_SIZES[index]
        requested < size || continue
        t = TRACE_BASIS_TRIANGLES[index]
        a = (t[1], t[2], t[3]); b = (t[4], t[5], t[6]); c = (t[7], t[8], t[9])
        # Bounding-box distance lower bound before the exact distance.
        bound = 0.0
        for d in 1:3
            low = min(a[d], b[d], c[d]); high = max(a[d], b[d], c[d])
            excess = max(low - (x, y, z)[d], (x, y, z)[d] - high, 0.0)
            bound += excess * excess
        end
        requested + slope * sqrt(bound) < size || continue
        size = min(size, requested + slope * point_triangle_distance((x, y, z), a, b, c))
    end
    return size
end

# Segment-distance size laws of the prism tube recipe, evaluated with the trace
# rule in the size callback (exact distances; Gmsh's Distance field samples curves).
# Tube axes: NormalSize at the tube surface (axis distance TUBE_BAND_OFFSET[] =
# R + pyramid height) growing with FarGrowth to FarSize. Band curves (the feature
# curves outside the tubes: junction lines, footprint edges, un-tubed edge parts):
# the metric band law NormalSize + RadialGrowth x min(r, ProtectedDistance) +
# FarGrowth x max(r - ProtectedDistance, 0), ProtectedDistance = 2 x NormalSize.
const TUBE_AXIS_SEGMENTS = NTuple{6, Float64}[]
const BAND_SEGMENTS = NTuple{6, Float64}[]
const TUBE_BAND_OFFSET = Ref(0.0)
const TUBE_BAND_NORMAL = Ref(Inf)
const TUBE_BAND_FAR = Ref(Inf)
const TUBE_BAND_FAR_GROWTH = Ref(0.0)
const BAND_RADIAL_GROWTH = 1.0
const BAND_PROTECTED_DISTANCE_OVER_NORMAL = 2.0
# Coupon-scale size bound recorded in the census (SizeBounds; see the option checks).
const SIZE_BOUND_RULE =
    "TangentialSize = min(--lc-tangent, FarSize): FarSize = FarSizeOverRadius x Radius is " *
    "the coarsest size the coupon admits (its matching-surface size), so the along-edge " *
    "tube extrusion / ridge spacing never exceeds it; dimensionless in Radius x " *
    "FarSizeOverRadius / --lc-tangent, identity for every coupon with FarSize >= --lc-tangent; " *
    "NormalSize and EdgeSize = CornerSize are resolution sizes and are never bounded (a fine " *
    "size above FarSize fails closed)"
# Corner-ball exterior law (decision 39b): NormalSize (the ball boundary size) +
# FarGrowth x the distance beyond CornerIsotropyRadius from the nearest graded
# point (semantic corners and tube cap centres), up to FarSize.
const CORNER_EXTERIOR_POINTS = NTuple{3, Float64}[]
const CORNER_EXTERIOR_RADIUS = Ref(0.0)

function segment_point_distance(x, y, z, segment)
    ax, ay, az, bx, by, bz = segment
    vx, vy, vz = bx - ax, by - ay, bz - az
    dx, dy, dz = x - ax, y - ay, z - az
    t = clamp((dx * vx + dy * vy + dz * vz) / (vx * vx + vy * vy + vz * vz), 0.0, 1.0)
    return sqrt((dx - t * vx)^2 + (dy - t * vy)^2 + (dz - t * vz)^2)
end

function tube_band_size(x, y, z, lc)
    size = lc
    normal = TUBE_BAND_NORMAL[]
    far = TUBE_BAND_FAR[]
    growth = TUBE_BAND_FAR_GROWTH[]
    isempty(TUBE_AXIS_SEGMENTS) && isempty(BAND_SEGMENTS) && return size
    axis_distance = Inf
    @inbounds for segment in TUBE_AXIS_SEGMENTS
        axis_distance = min(axis_distance, segment_point_distance(x, y, z, segment))
    end
    if isfinite(axis_distance)
        size = min(size, min(far, normal + growth * max(axis_distance - TUBE_BAND_OFFSET[], 0.0)))
    end
    return feature_band_size(x, y, z, size)
end

# The band law of the feature curves outside the tubes and the corner-ball
# exterior law, never above `lc` (the callback laws without the tube rule).
function feature_band_size(x, y, z, lc)
    size = lc
    normal = TUBE_BAND_NORMAL[]
    far = TUBE_BAND_FAR[]
    growth = TUBE_BAND_FAR_GROWTH[]
    band_distance = Inf
    @inbounds for segment in BAND_SEGMENTS
        band_distance = min(band_distance, segment_point_distance(x, y, z, segment))
    end
    if isfinite(band_distance)
        size = min(size, band_law_size(band_distance, normal, far, growth))
    end
    corner_distance = Inf
    @inbounds for point in CORNER_EXTERIOR_POINTS
        corner_distance = min(corner_distance, sqrt((x - point[1])^2 + (y - point[2])^2 + (z - point[3])^2))
    end
    if isfinite(corner_distance)
        size = min(size, corner_exterior_size(corner_distance, normal, far, growth))
    end
    return size
end

# Size prescribed on a tube axis (supervisor decision 40): the composed size
# field of the callback evaluated on the edge, without the tube rule (which
# describes the tube's exterior and would read NormalSize on the axis): the
# minimum of TangentialSize, the corner-ball law of the semantic corners along
# the edge (corner_curve_size: the ball grading from CornerSize and the
# process-band slope beyond the radius), the trace-basis volume law, the band law
# of the feature curves outside the tubes and the corner-ball exterior law (whose
# graded points include the tube cap centres: NormalSize within
# CornerIsotropyRadius of a cap). The cap centres are NOT ball-law points of the
# axis: grading the prism layers to CornerSize at a cap makes the lateral pyramid
# faces CornerSize x outer-arc slivers by construction, and the corner-ball
# tetrahedra against them fail the scaled-Jacobian gate (four-edge: 38 cells below
# 0.01, minimum 0.0060). The tube's layer boundaries equidistribute the arclength
# integral of the reciprocal size (graded_tube_stations).
function tube_axis_size(point, lc_tangent, corners, grading::CornerGrading, slope)
    size = corner_curve_size(point, corners, grading, lc_tangent, slope)
    size = trace_basis_size(point[1], point[2], point[3], size)
    return feature_band_size(point[1], point[2], point[3], size)
end

const TUBE_LAYER_RULE =
    "layer thickness = the composed size field on the tube axis (TubeAxisSizeLaw), " *
    "gradient-limited along the axis to (GrowthRatio - 1) / GrowthRatio and " *
    "equidistributed in the arclength integral of its reciprocal with ceil(integral) " *
    "layers, so every layer is <= the size it spans (<= TangentialSize); at a tube end on " *
    "the outer box the end layer is the field at the surface (its minimum over the layer " *
    "span) - decision 41; consecutive layers differ by at most GrowthRatio " *
    "(MaximumNeighbourRatio, fail closed) - decision 40"
const TUBE_AXIS_SIZE_LAW =
    "min(TangentialSize, corner-ball law of the semantic corners along the edge [CornerSize " *
    "with GrowthRatio to NormalSize inside CornerIsotropyRadius, the process-band slope " *
    "beyond], trace-basis volume rule, band rule, corner exterior rule [NormalSize within " *
    "CornerIsotropyRadius of a semantic corner or tube cap centre, FarGrowth beyond]); " *
    "excluded: the tube rule (the tube's exterior, NormalSize on the axis) and the cap " *
    "centres as ball-law points (CornerSize layers at a cap make the pyramid faces slivers " *
    "and the corner-ball tetrahedra fail the scaled-Jacobian gate)"

# The band law at distance r from a feature curve outside the tubes.
function band_law_size(r, normal, far, growth)
    protected = BAND_PROTECTED_DISTANCE_OVER_NORMAL * normal
    return min(far, normal + BAND_RADIAL_GROWTH * min(r, protected) + growth * max(r - protected, 0.0))
end

# The corner-ball exterior law at distance d from a graded point.
corner_exterior_size(d, normal, far, growth) =
    min(far, normal + growth * max(d - CORNER_EXTERIOR_RADIUS[], 0.0))

function prepare_tube_band_sizing!(axis_segments, band_segments, offset, normal, far, far_growth,
                                   graded_points=NTuple{3, Float64}[], corner_radius=0.0)
    empty!(TUBE_AXIS_SEGMENTS); empty!(BAND_SEGMENTS); empty!(CORNER_EXTERIOR_POINTS)
    for point in graded_points
        push!(CORNER_EXTERIOR_POINTS, (point[1], point[2], point[3]))
    end
    CORNER_EXTERIOR_RADIUS[] = corner_radius
    for (a, b) in axis_segments
        push!(TUBE_AXIS_SEGMENTS, (a[1], a[2], a[3], b[1], b[2], b[3]))
    end
    for (a, b) in band_segments
        push!(BAND_SEGMENTS, (a[1], a[2], a[3], b[1], b[2], b[3]))
    end
    TUBE_BAND_OFFSET[] = offset; TUBE_BAND_NORMAL[] = normal; TUBE_BAND_FAR[] = far
    TUBE_BAND_FAR_GROWTH[] = far_growth
    return Dict{String, Any}(
        "TubeAxisSegments" => length(TUBE_AXIS_SEGMENTS),
        "TubeAxisLength" => sum(norm([s[4] - s[1], s[5] - s[2], s[6] - s[3]]) for s in TUBE_AXIS_SEGMENTS; init=0.0),
        "TubeSurfaceOffset" => offset, "NormalSize" => normal, "FarSize" => far,
        "FarGrowth" => far_growth,
        "TubeRule" => "size = min(FarSize, NormalSize + FarGrowth x max(distance to the tube " *
                      "axis - TubeSurfaceOffset, 0)); TubeSurfaceOffset = tube radius + pyramid height",
        "BandSegments" => length(BAND_SEGMENTS),
        "BandLength" => sum(norm([s[4] - s[1], s[5] - s[2], s[6] - s[3]]) for s in BAND_SEGMENTS; init=0.0),
        "RadialGrowth" => BAND_RADIAL_GROWTH,
        "ProtectedDistance" => BAND_PROTECTED_DISTANCE_OVER_NORMAL * normal,
        "BandRule" => "size = min(FarSize, NormalSize + RadialGrowth x min(r, ProtectedDistance) + " *
                      "FarGrowth x max(r - ProtectedDistance, 0)), r the distance to the nearest " *
                      "feature curve outside the tubes (junction lines, footprint edges, un-tubed " *
                      "metal edge parts); ProtectedDistance = 2 x NormalSize",
        "JunctionVolumeRule" => "the junction lines carry the band law throughout the volume " *
                                "(BandRule, exact segment distance in the size callback) and are " *
                                "1D-meshed at NormalSize (BandCurves) - decision 39c",
        "CornerExteriorPoints" => length(CORNER_EXTERIOR_POINTS),
        "CornerIsotropyRadius" => corner_radius,
        "CornerExteriorGrowth" => far_growth,
        "CornerExteriorRule" => "size = min(FarSize, NormalSize + FarGrowth x max(d - " *
                                "CornerIsotropyRadius, 0)), d the distance to the nearest graded " *
                                "point (semantic corner or tube cap centre); the ball interior is " *
                                "the background graded-point law - decision 39b",
        "TraceBasisVolumeGrowth" => far_growth,
        "TraceBasisVolumeRule" => "size = min(FarSize, TraceBasisSizeRatio x the minimum " *
                                  "altitude of the basis triangle + FarGrowth x the distance " *
                                  "to the basis triangle) throughout the volume " *
                                  "(TraceBasisSizing.GradingSlope = FarGrowth) - decisions 39a, 43",
        "Composition" => "min(background field [anisotropic curve attractor, graded-point law], " *
                         "trace rule, tube rule, band rule, corner exterior rule) in the Gmsh " *
                         "size callback")
end

function trace_basis_size_callback(dim, tag, x, y, z, lc, data)::Cdouble
    return tube_band_size(x, y, z, trace_basis_size(x, y, z, lc))
end

function install_trace_basis_callback!()
    callback = @cfunction(trace_basis_size_callback, Cdouble,
                          (Cint, Cint, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cvoid}))
    push!(TRACE_BASIS_CALLBACK_ROOTS, callback)
    ierr = Ref{Cint}()
    ccall((:gmshModelMeshSetSizeCallback, gmsh.lib), Cvoid,
          (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cint}), callback, C_NULL, ierr)
    ierr[] == 0 || error(gmsh.logger.getLastError())
    return nothing
end

# Loads the narrow basis triangles into the callback state; returns the census record.
function prepare_trace_basis_sizing!(basis, ratio, lc_far, slope)
    isfinite(ratio) && ratio > 0.0 || error("Trace basis size ratio must be a positive finite number")
    empty!(TRACE_BASIS_TRIANGLES); empty!(TRACE_BASIS_SIZES); empty!(TRACE_BASIS_APEXES)
    TRACE_BASIS_FAR[] = lc_far; TRACE_BASIS_SLOPE[] = slope
    requested = Float64[]
    altitudes = Float64[]
    needles = 0
    narrow_needles = 0
    edges = Dict{Tuple{Int, Int}, Float64}()
    for triangle in basis.triangles
        lengths = trace_basis_edge_lengths(basis.points, triangle)
        altitude = trace_basis_minimum_altitude(basis.points, triangle, lengths)
        push!(altitudes, altitude)
        push!(requested, ratio * altitude)
        if altitude < NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE * minimum(lengths)
            needles += 1
            requested[end] < lc_far && (narrow_needles += 1)
        end
        for k in 1:3
            key = minmax(triangle[k], triangle[mod1(k + 1, 3)])
            edges[key] = lengths[k]
        end
        if requested[end] < lc_far
            push!(TRACE_BASIS_TRIANGLES, ntuple(i -> basis.points[triangle[(i - 1) ÷ 3 + 1]][mod1(i, 3)], 9))
            push!(TRACE_BASIS_SIZES, requested[end])
            # The narrow hat's apexes: the endpoints of the shortest edge, the two
            # vertices opposite the two longest edges, whose hats are the steepest.
            shortest = argmin(lengths)
            for vertex in (triangle[shortest], triangle[mod1(shortest + 1, 3)])
                apex = Tuple(Float64.(basis.points[vertex]))
                apex in TRACE_BASIS_APEXES || push!(TRACE_BASIS_APEXES, apex)
            end
        end
    end
    edge_lengths = collect(values(edges))
    return Dict{String, Any}(
        "Ratio" => ratio,
        "RatioIsDimensionless" => true,
        "Rule" => "element size <= TraceBasisSizeRatio x the minimum altitude of the basis " *
                  "triangle nearest the point (per-triangle rule: the hat varies over the " *
                  "whole triangle) + GradingSlope x the distance to that triangle, up to " *
                  "FarSize, throughout the volume (GradingSlope = FarGrowth under the prism " *
                  "tube recipe, decision 39a; the process-band slope otherwise); the only " *
                  "size parameter is the dimensionless ratio, applied through the Gmsh size " *
                  "callback on top of the background field",
        "SizeMeasure" => TRACE_BASIS_SIZE_MEASURE,
        "NeedleRule" => TRACE_BASIS_NEEDLE_RULE,
        "NeedleAltitudeOverShortestEdge" => NEEDLE_ALTITUDE_OVER_SHORTEST_EDGE,
        "NeedleTriangles" => needles,
        "NeedleTrianglesBelowFarSize" => narrow_needles,
        "MinimumBasisAltitude" => minimum(altitudes),
        "Model" => basis.model,
        "Frame" => [collect(basis.frame[r, :]) for r in 1:3],
        "InputSHA256" => basis.hashes,
        "Lower" => collect(basis.lower), "Upper" => collect(basis.upper),
        "Vertices" => length(basis.points), "Triangles" => length(basis.triangles),
        "UniqueEdges" => length(edge_lengths),
        "MinimumBasisEdge" => minimum(edge_lengths),
        "BasisEdgesBelowFarSize" => count(<(lc_far), edge_lengths),
        "TrianglesBelowFarSize" => length(TRACE_BASIS_TRIANGLES),
        "MinimumRequestedSize" => minimum(requested),
        "FarSize" => lc_far, "GradingSlope" => slope,
        "Apexes" => length(TRACE_BASIS_APEXES),
        "MeshFrameTriangles" => [[collect(basis.points[v]) for v in triangle] for triangle in basis.triangles])
end

# ---------------------------------------------------------------------------
# Prism edge tubes (supervisor decisions 37 and 38): the Gmsh-only recipe
# replaces the tetrahedral edge layer by a localized prism tube around every
# straight metal edge (top and bottom edge of every metal segment). The tube
# cross-section has geometric rings of size EdgeSize x GrowthRatio^(k-1) (the
# corner-ball shells with CornerSize = EdgeSize have the same law), is extruded at
# the tangential spacing, and its lateral quadrangles carry explicit pyramids
# (prism_edge_tubes.jl). A tube ends at the outer box where the edge does, and
# TubeCornerClearance before a semantic corner: the largest tube-free clearance
# any two tubes meeting at that corner need, R / tan(phi / 2) for the in-plane
# angle phi between the edges, plus one outer ring size. The un-tubed edge part
# inside the corner ball is meshed with tetrahedra graded from EdgeSize at the
# tube cap centre (the cap centre is a graded point of the corner law) and from
# CornerSize at the corner, so the cap triangles meet cells of their own size.

# Recipe scope (supervisor decision 48): every fail-closed guard of the prism-tube
# recipe is a recorded statement. RECIPE_SCOPE_GUARDS lists the geometry classes
# the recipe does not build, keyed by a stable id that the guard's error message
# carries as "ScopeGuard[<id>]" (scope_error), so that the build drivers record an
# unsupported class distinctly from any other failure; DetectedFrom says whether
# the class is visible in the frozen inputs ("inputs") or only in a derived
# quantity during the build ("build"). RECIPE_SCOPE_SUPPORTED_CLASSES are the
# input classes the recipe builds; exhibited_scope_classes lists the classes an
# input exhibits from the frozen inputs alone, except UntubedShortEdges (block (b)
# design A9 family 3, decision 304): a metal side whose tube interval is shorter
# than the tube's inner ring size carries no tube (its corner balls mesh it at
# CornerSize), a class known only once the clearances are derived and exhibited
# by the Scope.UntubedEdges records. The census Scope block records them and
# mesh_stage_contract.py binds the same lists.
const RECIPE_SCOPE_RECIPE = "prism-tubes"
const RECIPE_SCOPE_SUPPORTED_CLASSES = [
    "ArcSides", "ContinuationVertices", "DeviceFootprint", "DownwardLayers", "ExteriorLoops",
    "HoleLoops", "MultipleConductors", "MultipleLayers", "MultipleSlots", "ThinMetal", "TraceBasis",
    "UntubedShortEdges"]
# Decision 391 MAJOR-2 (ii): the arc tube paths that no full build has exercised fail
# closed until the four synthetic full builds of the next mesher lane (arc face ends at
# 45 / 70 degrees, a kinked arc / line joint, a 5e-5-rad smooth joint, a 2e-4-rad corner
# joint) lift the guards. The tested range of an arc joint's turn is the loop end
# 1b26671c9080's: 32 smooth arc / line joints turning by 3.0e-8 .. 1.6e-6 rad on the
# signature (the mesher's post-snap tilts read <= 1.27e-6 rad on its built thin and
# measurement-only fabricated meshes) and the exact part splits of one arc (turn 0).
const ARC_JOINT_TURN_BOUND = 1.6e-6
const RECIPE_SCOPE_GUARDS = [
    ("ArcTubeRadiusVsCurvature", "build",
     "an arc metal side whose radius is below four times the tube envelope (Radius + " *
     "PyramidHeight): the revolved sections would fold (block (b) design 1.2 (3))"),
    ("ArcFaceEnds", "build",
     "an arc metal side with an end on the outer box (a box-face cut end at any tilt, an end " *
     "exactly perpendicular to the face or a box-vertex corner): no full build has exercised " *
     "an arc tube reaching a box face (supervisor decision 391 MAJOR-2 (ii); lifted by the " *
     "synthetic arc face-end builds at 45 / 70 degrees)"),
    ("ArcJointTilt", "build",
     "an arc metal side meeting another side (straight or arc) at a joint whose turn (the " *
     "angle between the arc's end tangent and the other side's direction) exceeds " *
     "$(ARC_JOINT_TURN_BOUND) rad, the loop end's tested range (smooth joints turning more, " *
     "kinked arc / line joints and arc corners): no full build has exercised a sheared or a " *
     "capped arc joint beyond it (supervisor decision 391 MAJOR-2 (ii); lifted by the " *
     "synthetic 5e-5-rad smooth, 2e-4-rad corner and kinked arc / line joint builds)"),
    ("TopRounding", "inputs",
     "rounded metal top edges (TopRounding > 0): the tube rings surround a sharp edge"),
    ("TrenchRounding", "inputs",
     "rounded trench edges (TrenchRounding > 0): the bottom tube splits its substrate and " *
     "vacuum sectors at a sharp trench wall"),
    ("SlopedSidewalls", "inputs",
     "sidewall angle below 90 degrees: the tube sections assume vertical metal faces"),
    ("NoTrench", "inputs",
     "Overetch 0: the bottom tube needs an etched trench on its substrate side"),
    ("ShallowTrench", "build",
     "trench shallower than the tube radius plus the pyramid height: the pyramids would " *
     "reach the trench floor"),
    ("NarrowTransverseBound", "build",
     "the tube inner size does not fit min(Overetch, MetalThickness / 2, " *
     "CornerIsotropyRadius): no ring fits"),
    ("NarrowHoles", "build",
     "a hole narrower than twice the tube reach (Radius + PyramidHeight + the band's " *
     "ProtectedDistance): the tubes facing each other across it would overlap"),
    ("NarrowLayerGap", "build",
     "a vacuum gap between the metal faces of an upward and a downward process layer " *
     "narrower than twice the tube reach: the tubes facing each other across it would overlap"),
    ("NarrowMetal", "build",
     "a metal strip narrower than twice the tube envelope (Radius + PyramidHeight) between " *
     "two tubed metal sides of one plane facing each other across the metal: the tube parts " *
     "over the metal would overlap (supervisor decision 347)"),
    ("SteepFaceCrossing", "build",
     "a tube end on a box face whose tilt lies beyond the validity ceiling of the capped end " *
     "block, 2 PyramidHeight |tan theta| >= the section's end-spacing cap lc_cap (the largest " *
     "axial spacing at which the section's own prism frames read <= 0.95 x the Jacobian-" *
     "condition ceiling): no sheared layer can keep its pyramid apex inside the box within " *
     "the cap (mesher design round 2 F2b 3.2, supervisor decisions 358 / 363 / 437)"),
    ("FreeEdgeEnds", "build",
     "a metal edge end that is neither a semantic corner nor on the outer box"),
    ("FootprintWithoutEdge", "build",
     "an explicit etch footprint with no side coincident with a metal edge"),
    ("FootprintTopology", "build",
     "a producer-default etch collar whose region is not one simple polygon (Physical sides " *
     "facing each other across a dielectric gap narrower than twice the collar enclose an " *
     "un-etched island larger than the island rule admits (COLLAR_ISLAND_EXCESS_CAP_OVER_RADIUS " *
     "x Radius beyond the collar) or touching the box face, or the region touches itself): " *
     "the collar loft takes one simple polygon")]
const RECIPE_SCOPE_RULE =
    "the prism-tube recipe builds every input whose classes are all in SupportedClasses; " *
    "an input exhibiting a class in GuardedClasses fails closed at the guard whose " *
    "error message carries ScopeGuard[<Id>] (Guards); ExhibitedClasses are the classes " *
    "of this input among both lists, from the frozen inputs (loops, layers, process), " *
    "UntubedShortEdges excepted (exhibited by a non-empty UntubedEdges: " *
    "UntubedShortEdgeRule); MetalLoops counts per plan-view loop the straight sides not on " *
    "the outer box plus the parts of its tagged arcs (ArcSides: one run of chords = one arc " *
    "side split into ceil(sweep / 90 degrees) parts), and TubeCount = TubesPerSide x their " *
    "sum minus the untubed short sides (a top and a bottom tube per side of every loop of a " *
    "fabricated coupon, one sheet tube per side of a thin coupon; decision 66)"
# Family 3 (design A9; decision 304): the threshold at interval == EdgeSize separates an
# untubed side from a one-layer tube of at least one inner ring in length (both
# valid); an interval <= 0 was the ShortEdges refusal before.
const UNTUBED_SHORT_EDGE_RULE =
    "a straight metal side whose tube interval Span - Clearances[1] - Clearances[2] is below " *
    "the tube's inner ring size EdgeSize carries no tube (an interval >= EdgeSize carries a " *
    "tube of at least one inner ring in length); its corner balls (CornerIsotropyRadius about " *
    "each semantic corner at its ends) mesh it with tetrahedra graded from CornerSize = " *
    "EdgeSize, CoveredByBalls when the balls reach over the whole span (CornerIsotropyRadius x " *
    "the corners at its ends >= Span), otherwise the uncovered middle follows the two-sided " *
    "corner law; recorded UntubedEdges[] {Side, Span, Clearances, Interval, Corners, " *
    "CoveredByBalls}"

scope_guard_statement(id) = RECIPE_SCOPE_GUARDS[findfirst(guard -> guard[1] == id, RECIPE_SCOPE_GUARDS)][3]

function scope_error(id, detail)
    any(guard -> guard[1] == id, RECIPE_SCOPE_GUARDS) || error("Unknown scope guard $id")
    error("ScopeGuard[$id]: $(scope_guard_statement(id)); $detail")
end

# The classes an input exhibits (sorted), from the plan-view loops, the signature
# layers and the process options, before anything is built.
function exhibited_scope_classes(edges, loops, layers, fabricated, sidewall_angle, top_rounding,
                                 trench_rounding, overetch, device_footprint, trace_basis_bound)
    classes = String[]
    any(!loop.hole for loop in loops) && push!(classes, "ExteriorLoops")
    any(loop.hole for loop in loops) && push!(classes, "HoleLoops")
    any(class == "Continuation" for loop in loops for class in loop.classes) &&
        push!(classes, "ContinuationVertices")
    length(unique(edge.slot for edge in edges)) > 1 && push!(classes, "MultipleSlots")
    length(unique(edge.conductor for edge in edges)) > 1 && push!(classes, "MultipleConductors")
    length(layers) > 1 && push!(classes, "MultipleLayers")
    any(layer.sign < 0 for layer in layers) && push!(classes, "DownwardLayers")
    device_footprint && push!(classes, "DeviceFootprint")
    trace_basis_bound && push!(classes, "TraceBasis")
    top_rounding > 0.0 && push!(classes, "TopRounding")
    trench_rounding > 0.0 && push!(classes, "TrenchRounding")
    sidewall_angle != 90.0 && push!(classes, "SlopedSidewalls")
    fabricated || push!(classes, "ThinMetal")
    overetch == 0.0 && push!(classes, "NoTrench")
    any(loop_has_arcs(loop) for loop in loops) && push!(classes, "ArcSides")
    return sort!(classes)
end

# A plan-view side lies on the outer box when both ends are on the same box face.
function side_on_box_face(p, q, lower, upper, tolerance)
    return any((abs(p[d] - lower[d]) <= tolerance && abs(q[d] - lower[d]) <= tolerance) ||
               (abs(p[d] - upper[d]) <= tolerance && abs(q[d] - upper[d]) <= tolerance)
               for d in 1:2)
end

# Per plan-view loop: the straight sides not on the outer box (each carries a top
# and a bottom tube) plus the parts of its tagged arcs (ArcSides, design 1.2 (2)); an arc
# loop records its arcs {ArcId, Centre, Radius, SweepDegrees, Parts, Chords, Sign}.
function metal_loop_records(loops, lower, upper, tolerance)
    records = Dict{String, Any}[]
    for (index, loop) in enumerate(loops)
        n = length(loop.points)
        runs = loop_has_arcs(loop) ? tagged_arc_runs(loop, tolerance) : NamedTuple[]
        arc_edges = falses(n)
        for run in runs
            arc_edges[run.edge_indices] .= true
        end
        sides = count(!arc_edges[i] && !side_on_box_face(loop.points[i], loop.points[i % n + 1], lower, upper,
                                                         tolerance) for i in 1:n)
        parts = sum(arc_part_count(run.sweep) for run in runs; init=0)
        record = Dict{String, Any}(
            "Loop" => index, "Conductor" => loop.conductor, "Plane" => loop.plane,
            "Hole" => loop.hole, "Vertices" => n, "Sides" => sides + parts)
        if !isempty(runs)
            record["StraightSides"] = sides
            record["ArcParts"] = parts
            record["Arcs"] = [Dict{String, Any}(
                "ArcId" => run.id, "Centre" => [run.center[1], run.center[2]], "Radius" => run.radius,
                "SweepDegrees" => rad2deg(run.sweep), "Parts" => arc_part_count(run.sweep),
                "Chords" => length(run.edge_indices), "Sign" => run.sign) for run in runs]
        end
        push!(records, record)
    end
    return records
end

function recipe_scope_record(exhibited, loops, lower, upper, tolerance, untubed_edges)
    return Dict{String, Any}(
        "Rule" => RECIPE_SCOPE_RULE, "Recipe" => RECIPE_SCOPE_RECIPE,
        "SupportedClasses" => copy(RECIPE_SCOPE_SUPPORTED_CLASSES),
        "GuardedClasses" => [guard[1] for guard in RECIPE_SCOPE_GUARDS],
        "Guards" => [Dict{String, Any}("Id" => id, "DetectedFrom" => detected,
                                       "Statement" => statement)
                     for (id, detected, statement) in RECIPE_SCOPE_GUARDS],
        "ExhibitedClasses" => copy(exhibited),
        "MetalLoops" => metal_loop_records(loops, lower, upper, tolerance),
        "UntubedShortEdgeRule" => UNTUBED_SHORT_EDGE_RULE,
        "UntubedEdges" => untubed_edges)
end

# The box face (d, side) a side end on the outer box is cut by: the face among those
# the end lies on that the side crosses most steeply (|direction[d]| largest). An end
# at a box corner whose side crosses both faces obliquely has no single face plane
# for its end block and fails closed.
function crossed_box_face(point, direction, lower, upper, tolerance)
    faces = Tuple{Int, Int}[]
    for d in 1:2
        abs(point[d] - lower[d]) <= tolerance && push!(faces, (d, 0))
        abs(point[d] - upper[d]) <= tolerance && push!(faces, (d, 1))
    end
    isempty(faces) && return nothing
    best = argmax([abs(direction[d]) for (d, _) in faces])
    for (k, (d, _)) in enumerate(faces)
        k == best && continue
        # At a box corner the side must lie IN the other face (direction[d] == 0).
        direction[d] == 0.0 ||
            error("metal edge end $point crosses the box at a box corner obliquely " *
                  "(direction $direction): no single face plane for its end block")
    end
    return faces[best]
end

box_face_name(d, side) = string(("x", "y")[d], side)

# Straight metal edges of the plan-view boundary loops: the Physical sides (the
# sides not lying on the outer box), with the horizontal normal pointing away
# from the metal and the tube interval shrunk by the corner clearance at semantic
# corners (0 at box continuation vertices). Returns rows with start, stop,
# direction, normal, span, s_start, s_end, plane, conductor, corner angles, the
# face ends, the legacy box-vertex corners and the untubed short-side class.
#
# Untubed short sides (design A9 family 3, decision 304): a side whose tube interval
# s_end - s_start is below `edge_size` (the tube's inner ring size) carries no tube -
# `untubed` true, s_start / s_end kept as derived (the interval may be negative), the
# side meshed by its corner balls at CornerSize; `covered_by_balls` says whether the
# balls of radius `corner_radius` about its semantic corners reach over the whole
# span (an untubed side between two acute corners can leave an uncovered middle,
# graded by the two-sided corner law). A side that is not untubed must keep a
# positive interval (with edge_size 0 an exactly empty interval fails closed).
#
# Box-face ends (design A2 / A6, supervisor decision 320): an end on the outer box
# where no other metal side of the plane meets the edge is a FACE END when the
# edge is not exactly perpendicular to the face - the exact-arithmetic test on the
# quantised canonical coordinates is that the edge's face-parallel coordinate
# differs between its two ends (theta > 0) - whatever the contract's vertex class:
# clearance 0, no ball, the sheared end block, contract class box-face cut end
# (derive_semantic_contract mirrors the rule; a contract still listing such an
# end as a corner fails closed here). An exactly perpendicular end (theta == 0,
# every rectilinear coupon) keeps the legacy convention bitwise: a Physical-class
# box vertex is a semantic corner with h_K + its ball (recorded
# LegacyBoxVertexCorner), a Continuation-class one has clearance 0. Two metal
# sides meeting at a box vertex are a corner at any theta. The theta -> 0 class
# boundary is a designed discontinuity between two valid treatments 3 R past the
# claims (decision 320; dropping the theta-0 box balls is a recorded follow-up).
function metal_edge_segments(loops, corners, clearance_of_angle, lower, upper, tolerance;
                             edge_size=0.0, corner_radius=0.0)
    segments = NamedTuple[]
    on_box(p) = any(abs(p[d] - lower[d]) <= tolerance || abs(p[d] - upper[d]) <= tolerance
                    for d in 1:2)
    is_corner(p, plane) = any(norm(collect(corner) .- [p[1], p[2], plane]) <= tolerance
                              for corner in corners)
    # Every side first - the straight sides and the ARC parts (design 1.2 (2) / A3: an arc
    # run of the tagged boundary is one ArcSide, split into arc_part_count equal angular
    # parts, every part a side) - so that the angle between the sides meeting at a corner is
    # known before the clearance is applied. Each side carries its two END descriptors
    # (point, the unit direction pointing AWAY from the side at that end) for the corner
    # angles; a smooth joint (a vertex whose JointSmooth tag is 1, or the split between two
    # parts of one arc) is never a corner: the two tubes share one cross-section there.
    sides = NamedTuple[]
    smooth_vertices = Set{Tuple{Float64, Float64, Float64}}()
    for loop in loops
        # The metal lies inside an exterior loop and outside a hole (loop.hole: the
        # polygon interior is dielectric), so the tube normal, which points away from
        # the metal, points out of an exterior loop and into a hole.
        metal_inside = !loop.hole
        n = length(loop.points)
        # Arc sides come from the TAGGED runs alone (design A1 (4)): an untagged loop (every
        # legacy coupon) is read as straight sides, bitwise.
        runs = loop_has_arcs(loop) ? tagged_arc_runs(loop, tolerance) : NamedTuple[]
        arc_edges = falses(n)
        for run in runs
            arc_edges[run.edge_indices] .= true
        end
        if haskey(loop, :joints)
            for (i, joint) in enumerate(loop.joints)
                joint === nothing || !joint.smooth || push!(smooth_vertices, (loop.points[i][1], loop.points[i][2], loop.plane))
            end
        end
        for i in 1:n
            arc_edges[i] && continue
            p = loop.points[i]
            q = loop.points[i % n + 1]
            side_on_box_face(p, q, lower, upper, tolerance) && continue
            direction = [q[1] - p[1], q[2] - p[2]]
            span = norm(direction)
            span > tolerance || error("Degenerate metal edge $p -> $q")
            direction ./= span
            normal = [direction[2], -direction[1]]
            midpoint = [0.5 * (p[1] + q[1]), 0.5 * (p[2] + q[2])]
            probe = midpoint .+ 1.0e-3 .* normal
            if point_in_polygon((probe[1], probe[2]), loop.points, tolerance) == metal_inside
                normal .*= -1.0
            end
            probe = midpoint .- 1.0e-3 .* normal
            point_in_polygon((probe[1], probe[2]), loop.points, tolerance) == metal_inside ||
                error("Unable to orient the metal edge $p -> $q")
            for (point, other) in ((p, q), (q, p))
                is_corner(point, loop.plane) || on_box(point) ||
                    (point[1], point[2], loop.plane) in smooth_vertices ||
                    scope_error("FreeEdgeEnds", "metal edge end $point of $p -> $q")
            end
            push!(sides, (kind=:straight, start=[p[1], p[2]], stop=[q[1], q[2]], direction=direction,
                          normal=normal, span=span, plane=loop.plane, conductor=loop.conductor,
                          hole=loop.hole,
                          start_corner=is_corner(p, loop.plane), stop_corner=is_corner(q, loop.plane),
                          away=(-direction, direction), tangents=(direction, direction),
                          arc=nothing))
        end
        for run in runs
            # The run's travel: from its first vertex to its last, sweep signed about the
            # centre; sigma = +1 when the dielectric lies outside the circle (convex metal).
            centre = [run.center[1], run.center[2]]
            rho = run.radius
            sweep = run.sweep
            abs(sweep) > 0.0 || error("arc $(run.id) has no sweep")
            p_first = loop.points[run.point_indices[1]]
            p_last = loop.points[run.point_indices[end]]
            theta_first = atan(p_first[2] - centre[2], p_first[1] - centre[1])
            travel = sign(sweep)
            # The metal side of the circle, probed on the radial line through the run's middle
            # chord VERTEX (a point of both the circle and the chord polygon: 1 nm inside the
            # circle lies inside the chord polygon of a convex arc, 1 nm outside lies in the
            # metal of a concave one; a probe at a chord midpoint would sit in the sagitta).
            vertex = loop.points[run.point_indices[cld(length(run.point_indices), 2)]]
            theta_mid = atan(vertex[2] - centre[2], vertex[1] - centre[1])
            probe_out = centre .+ (rho + 1.0e-3) .* [cos(theta_mid), sin(theta_mid)]
            probe_in = centre .+ (rho - 1.0e-3) .* [cos(theta_mid), sin(theta_mid)]
            outside_metal = point_in_polygon((probe_out[1], probe_out[2]), loop.points, tolerance) == metal_inside
            inside_metal = point_in_polygon((probe_in[1], probe_in[2]), loop.points, tolerance) == metal_inside
            outside_metal != inside_metal || error("Unable to orient the arc $(run.id) of conductor $(loop.conductor)")
            sigma = inside_metal ? 1.0 : -1.0
            Int(sigma) == run.sign ||
                error("arc $(run.id): the plan-view metal side (sigma $(Int(sigma))) disagrees with the tagged " *
                      "ArcSign $(run.sign)")
            chords = [([loop.points[run.point_indices[k]]...], [loop.points[run.point_indices[k + 1]]...])
                      for k in 1:(length(run.point_indices) - 1)]
            length(chords) >= 4 ||
                error("arc $(run.id) of conductor $(loop.conductor) has $(length(chords)) chords: the untagged " *
                      "arc fit of the trench footprint needs at least four (a short arc is not supported yet)")
            angles = arc_split_angles(loop.points, run.point_indices, run.center, sweep)
            parts = length(angles) - 1
            for k in 1:parts
                theta_a, theta_b = angles[k], angles[k + 1]
                a = k == 1 ? [p_first[1], p_first[2]] : centre .+ rho .* [cos(theta_a), sin(theta_a)]
                b = k == parts ? [p_last[1], p_last[2]] : centre .+ rho .* [cos(theta_b), sin(theta_b)]
                tangent(theta) = travel .* [-sin(theta), cos(theta)]
                part_sweep = theta_b - theta_a
                span = rho * abs(part_sweep)
                normal_a = sigma .* [cos(theta_a), sin(theta_a)]
                for (point, first_part, last_part) in ((a, k == 1, false), (b, false, k == parts))
                    (first_part || last_part) || continue
                    is_corner(point, loop.plane) || on_box(point) ||
                        (point[1], point[2], loop.plane) in smooth_vertices ||
                        scope_error("FreeEdgeEnds", "arc $(run.id) end $point")
                end
                push!(sides, (kind=:arc, start=a, stop=b, direction=tangent(theta_a), normal=normal_a,
                              span=span, plane=loop.plane, conductor=loop.conductor, hole=loop.hole,
                              start_corner=k == 1 && is_corner(a, loop.plane),
                              stop_corner=k == parts && is_corner(b, loop.plane),
                              away=(-tangent(theta_a), tangent(theta_b)),
                              tangents=(tangent(theta_a), tangent(theta_b)),
                              arc=(id=run.id, centre=centre, rho=rho, sigma=sigma, theta_start=theta_a,
                                   theta_end=theta_b, sweep=part_sweep, part=k, parts=parts,
                                   run_sweep=sweep, chords=chords)))
            end
        end
    end
    isempty(sides) && error("No straight metal edges found")
    # In-plane angle between two tube edges meeting at a corner: the smallest
    # angle between the directions pointing away from the corner (pi when the
    # corner has a single tube edge). An arc's away direction is its end tangent.
    function corner_angle(point, plane, own_direction_away)
        angle = Float64(pi)
        for side in sides
            abs(side.plane - plane) <= tolerance || continue
            for (end_point, away) in ((side.start, side.away[1]), (side.stop, side.away[2]))
                norm(end_point .- point) <= tolerance || continue
                dot(away, own_direction_away) >= 1.0 - 1.0e-12 && continue
                angle = min(angle, acos(clamp(dot(away, own_direction_away), -1.0, 1.0)))
            end
        end
        return angle
    end
    # Metal sides of the plane meeting at a point (the side itself included).
    function sides_at(point, plane)
        return count(abs(side.plane - plane) <= tolerance &&
                     (norm(side.start .- point) <= tolerance || norm(side.stop .- point) <= tolerance)
                     for side in sides)
    end
    # The face end of a side end on the box, or nothing: (face name, outward normal,
    # theta, axis, value) when the end has a single metal side and theta > 0 (the tangent
    # at the end for an arc).
    function face_end(point, side, tangent, segment_end)
        on_box(point) || return nothing
        sides_at(point, side.plane) == 1 || return nothing
        face = crossed_box_face(point, tangent, lower, upper, tolerance)
        face === nothing && return nothing
        d, box_side = face
        if side.kind == :straight
            # theta == 0 exactly: the face-parallel coordinate is the same at both ends.
            side.start[3 - d] == side.stop[3 - d] && return nothing
        else
            abs(tangent[d]) == 1.0 && return nothing
        end
        normal = [0.0, 0.0]
        normal[d] = box_side == 0 ? -1.0 : 1.0
        theta = acos(clamp(abs(tangent[d]), 0.0, 1.0))
        return (face=box_face_name(d, box_side), normal=normal, theta=theta, axis=d,
                value=box_side == 0 ? lower[d] : upper[d])
    end
    # A smooth joint at a side end: the OTHER side of the plane meeting the end there
    # (its index in `sides` and which of its ends), when the vertex is a smooth joint vertex
    # or the end is the split between two parts of one arc; nothing otherwise.
    function smooth_joint(side_index, point, plane)
        side = sides[side_index]
        split = side.kind == :arc && (
            (norm(side.start .- point) <= tolerance && side.arc.part > 1) ||
            (norm(side.stop .- point) <= tolerance && side.arc.part < side.arc.parts))
        split || (point[1], point[2], plane) in smooth_vertices || return nothing
        partners = Tuple{Int, Int}[]
        for (j, other) in enumerate(sides)
            j == side_index && continue
            abs(other.plane - plane) <= tolerance || continue
            norm(other.start .- point) <= tolerance && push!(partners, (j, 1))
            norm(other.stop .- point) <= tolerance && push!(partners, (j, 2))
        end
        length(partners) == 1 ||
            error("the smooth joint at $point of plane $plane meets $(length(partners)) other sides (one expected)")
        side.kind == :arc || sides[partners[1][1]].kind == :arc ||
            error("the smooth joint at $point joins two straight sides: the builder merges those (design A3 (1))")
        return partners[1]
    end
    # Decision 391 MAJOR-2 (ii): an arc side is built only in the tested configuration -
    # no end on the outer box, and every joint with another side of its plane turning by at
    # most ARC_JOINT_TURN_BOUND (a tangent joint turns by 0; the split between two parts of
    # one arc is exact). The turn is read on the sides' away directions (the arc's end
    # tangent against the other side's direction leaving the joint), by atan of the cross
    # and dot products (exact to the rounding of the directions, unlike acos near 1).
    function arc_end_guards(side)
        side.kind == :arc || return
        for (point, away) in ((side.start, side.away[1]), (side.stop, side.away[2]))
            on_box(point) && scope_error("ArcFaceEnds",
                                          "arc $(side.arc.id) part $(side.arc.part) of conductor " *
                                          "$(side.conductor) ends at $point on the outer box")
            for other in sides
                other === side && continue
                abs(other.plane - side.plane) <= tolerance || continue
                for (other_point, other_away) in ((other.start, other.away[1]), (other.stop, other.away[2]))
                    norm(other_point .- point) <= tolerance || continue
                    continuing = -other_away
                    turn = atan(abs(cross2d((away[1], away[2]), (continuing[1], continuing[2]))),
                                dot(away, continuing))
                    turn <= ARC_JOINT_TURN_BOUND ||
                        scope_error("ArcJointTilt",
                                    "arc $(side.arc.id) part $(side.arc.part) of conductor $(side.conductor) " *
                                    "meets a $(other.kind == :arc ? "arc" : "straight") side at $point " *
                                    "with a turn of $turn rad")
                end
            end
        end
    end
    for side in sides
        arc_end_guards(side)
    end
    for (index, side) in enumerate(sides)
        start_joint = smooth_joint(index, side.start, side.plane)
        stop_joint = smooth_joint(index, side.stop, side.plane)
        start_face = start_joint === nothing ? face_end(side.start, side, side.tangents[1], 1) : nothing
        stop_face = stop_joint === nothing ? face_end(side.stop, side, side.tangents[2], 2) : nothing
        for (point, face, corner) in ((side.start, start_face, side.start_corner),
                                      (side.stop, stop_face, side.stop_corner))
            face === nothing || !corner ||
                error("the semantic contract lists the box-face cut end $point (face " *
                      "$(face.face), tilt $(rad2deg(face.theta)) degrees) as a semantic corner: " *
                      "regenerate the contract (decision 320: a box-face end with a single " *
                      "metal side that is not exactly perpendicular to the face is a cut end)")
        end
        for (point, joint, corner) in ((side.start, start_joint, side.start_corner),
                                       (side.stop, stop_joint, side.stop_corner))
            joint === nothing || !corner ||
                error("the semantic contract lists the smooth joint $point as a semantic corner: " *
                      "regenerate the contract (design A3 (1): a joint turning by at most 1e-4 rad shares " *
                      "one tube section, no ball)")
        end
        start_corner = side.start_corner && start_face === nothing
        stop_corner = side.stop_corner && stop_face === nothing
        start_angle = start_corner ? corner_angle(side.start, side.plane, side.away[1]) : Float64(pi)
        stop_angle = stop_corner ? corner_angle(side.stop, side.plane, side.away[2]) : Float64(pi)
        s_start = start_corner ? clearance_of_angle(start_angle) : 0.0
        s_end = side.span - (stop_corner ? clearance_of_angle(stop_angle) : 0.0)
        untubed = s_end - s_start < edge_size
        untubed || s_end > s_start ||
            error("metal edge $(side.start) -> $(side.stop) of span $(side.span) leaves no tube " *
                  "interval against clearances $(s_start) and $(side.span - s_end)")
        untubed && side.kind == :arc &&
            error("arc $(side.arc.id) part $(side.arc.part) of span $(side.span) leaves a tube interval below " *
                  "the inner ring size against its clearances: an untubed arc is not supported")
        covered_by_balls = corner_radius * (start_corner + stop_corner) >= side.span
        # A theta-0 box vertex kept as a corner by the legacy convention.
        legacy_box_corners = (start_corner && on_box(side.start) && start_angle >= pi - 1.0e-9,
                              stop_corner && on_box(side.stop) && stop_angle >= pi - 1.0e-9)
        push!(segments, (side..., s_start=s_start, s_end=s_end,
                         corner_angles=(start_angle, stop_angle),
                         face_ends=(start_face, stop_face),
                         joints=(start_joint, stop_joint),
                         legacy_box_corners=legacy_box_corners,
                         corners=(start_corner, stop_corner),
                         untubed=untubed, covered_by_balls=covered_by_balls))
    end
    return segments
end

# Distance between two plan-view segments (0 when they intersect).
function segment_segment_distance_2d(a, b, c, d)
    orientation(p, q, r) = sign(cross2d((q[1] - p[1], q[2] - p[2]), (r[1] - p[1], r[2] - p[2])))
    if orientation(a, b, c) * orientation(a, b, d) < 0 && orientation(c, d, a) * orientation(c, d, b) < 0
        return 0.0
    end
    return min(point_segment_distance_2d(a, c, d), point_segment_distance_2d(b, c, d),
               point_segment_distance_2d(c, a, b), point_segment_distance_2d(d, a, b))
end

# The closest pair of points of two plan-view segments (p on a-b, q on c-d), among
# the projections of each endpoint onto the other segment.
function segment_closest_points_2d(a, b, c, d)
    function project(point, first, second)
        direction = (second[1] - first[1], second[2] - first[2])
        span = direction[1]^2 + direction[2]^2
        span > 0.0 || return (first[1], first[2])
        t = clamp(((point[1] - first[1]) * direction[1] + (point[2] - first[2]) * direction[2]) / span,
                  0.0, 1.0)
        return (first[1] + t * direction[1], first[2] + t * direction[2])
    end
    candidates = [((a[1], a[2]), project(a, c, d)), ((b[1], b[2]), project(b, c, d)),
                  (project(c, a, b), (c[1], c[2])), (project(d, a, b), (d[1], d[2]))]
    return candidates[argmin([hypot(p[1] - q[1], p[2] - q[2]) for (p, q) in candidates])]
end

# ---------------------------------------------------------------------------
# SEAM (mesher design round 2 section 5 / FIX-UPS F.7; supervisor decisions 350 / 351 /
# 358 / 363): the acute thin-tip pinched seam. Palace opens the thin sheet as a crack
# by duplicating its interior vertices; a sheet edge shared by two sheet triangles
# whose two vertices both lie on the sheet's FREE boundary (never duplicated) is a
# coarse crack edge that geodata.cpp refine_crack_elements must split conformally,
# which a prism / pyramid tube mesh refuses (the O1 thin and the stored O4 thin V7
# meshes: 6 / 1 such seams at their 22.5 / 45-degree tips). At a convex sheet tip of
# opening phi < 90 degrees the Delaunay fan at the tip is one triangle (its opposite
# edge joins the two sheet edges) and seams follow along the wedge while the local
# size exceeds the wedge width. The fix embeds a TIP BISECTOR curve in the sheet
# surface (an OCC line fragmented with the sheet: Gmsh records it as an embedded
# curve of the sheet face, so every triangle edge crossing the wedge meets a bisector
# node and every fan / wedge triangle has an interior, duplicable vertex) from the
# tip to the first station d_b where the half width d sin(phi / 2) reaches the corner
# law s(d) (CornerSize -> NormalSize inside the ball, FarGrowth beyond), capped at the
# tube start station (beyond it the sheet edges are tube axes inside the structured
# tubes and no free-sheet triangle can join them). Exact predicate: a semantic corner
# whose two sides have a positive quantised dot product (opening < 90 degrees) with
# the metal in the acute sector between them. Thin kind under the prism-tube recipe.
const TIP_BISECTOR_RULE =
    "a convex thin-sheet tip whose two plan-view sides have a positive quantised dot product " *
    "(opening phi < 90 degrees) with the metal in the acute sector between them carries an " *
    "embedded bisector curve in the sheet surface from the tip to the first station d_b with " *
    "d_b sin(phi / 2) >= s(d_b) (s the corner law: CornerSize graded to NormalSize inside " *
    "CornerIsotropyRadius, FarGrowth beyond), capped at the tube start station of the shorter " *
    "clearance (the station whose along-edge projection reaches the side's tube start; beyond " *
    "it the sheet edges are tube axes); the curve clears every tube envelope (fail closed " *
    "otherwise) and is meshed by the composed size field, so no sheet triangle edge joins the " *
    "two sheet edges: seam count 0 by conformity (ThinSheetSeams); supervisor decisions 350 / " *
    "351 / 358 / 363 (mesher design round 2, SEAM)"

const TIP_BISECTOR_STATION_STEP = 1.0e-4

# Supervisor decision 368: the bisector splits the tip fan into two phi / 2 sectors, so
# the corner cells on those sectors carry a kappa_reg FLOOR of their own - the 2D
# condition number of the affine map from the equilateral triangle onto the isoceles
# triangle of apex angle alpha = phi / 2 (singular-value interlacing: no tetrahedron on
# such a sector reads below it; exact: 5.862 at phi 22.5, 5.272 at 25, 4.385 at 30,
# reproduced by a free 3D minimisation to 4 digits). The bounded descent at the
# production sizes reaches the floor within TIP_BISECTOR_DESCENT_EXCESS (the synthetic
# production-size thin family 26-45 degrees, three exit tilts: measured / floor at most
# 1.19; seed-impl REPORT section 2). The bisector is therefore embedded only at tips of
# opening phi >= TipBisectorMinimumOpening = 2 alpha* with f(alpha*) = CornerShapeGate /
# TIP_BISECTOR_DESCENT_EXCESS (31.556 degrees at the gate 5.0: 31 / 32 degrees measured
# 4.85 / 4.82 with the bisector, 29 / 30 up to 5.41 / 5.04 - the law's crossing); a
# sharper convex tip keeps its fan and its pinched seams, which the mesher records as
# UnrefinedCrackSeams for the solve config (Model.RefineCrackElements false; FIX-UPS F.7:
# every DOF of the PEC sheet, the seam edge's included, is Dirichlet, so the coupling an
# unrefined seam leaves is between two constrained values).
const TIP_BISECTOR_DESCENT_EXCESS = 1.2

function sector_regular_floor(alpha)
    matrix = [1.0 (2.0 * cos(alpha) - 1.0) / sqrt(3.0); 0.0 2.0 * sin(alpha) / sqrt(3.0)]
    singular = svdvals(matrix)
    return singular[1] / singular[end]
end

function tip_bisector_minimum_opening(corner_shape_gate)
    target = corner_shape_gate / TIP_BISECTOR_DESCENT_EXCESS
    target > sector_regular_floor(pi / 4) || return Float64(pi)
    low, high = 1.0e-6, pi / 4
    for _ in 1:200
        middle = 0.5 * (low + high)
        sector_regular_floor(middle) > target ? (low = middle) : (high = middle)
    end
    return 2.0 * 0.5 * (low + high)
end

const UNREFINED_CRACK_SEAMS_RULE =
    "a convex thin tip sharper than TipBisectorMinimumOpening keeps its one-triangle fan and the " *
    "pinched seams along its wedge (every seam within the tip's tube start station, fail closed " *
    "otherwise); the solve config of such a thin coupon sets Model.RefineCrackElements false " *
    "(every DOF of the PEC sheet, the seam edge's included, is Dirichlet: the coupling an " *
    "unrefined seam leaves is between two constrained values, with no effect on the potential " *
    "or the per-side charge integrals) - supervisor decision 368 (3); the (F) qualification of " *
    "such a coupon is the acceptance"

# The bisector tools of every convex thin tip at or above the minimum opening (OCC
# lines, to be fragmented with the sheet), their records, and the sharper tips whose
# seams stay unrefined.
function tip_bisector_tools!(occ, corners, kinds, sides, segments, radius_envelope,
                             grading::CornerGrading, lc_far, far_growth, tolerance;
                             corner_shape_gate=0.0)
    tools = Tuple{Int32, Int32}[]
    records = Dict{String, Any}[]
    unrefined = Dict{String, Any}[]
    minimum_opening = tip_bisector_minimum_opening(corner_shape_gate)
    size_law(d) = d <= grading.radius ? corner_ball_size(grading, d) :
                  min(lc_far, grading.lc_fine + far_growth * (d - grading.radius))
    for (k, (corner, kind, side)) in enumerate(zip(corners, kinds, sides))
        kind === :invariant || continue
        a, b = side.walls
        product = a[1] * b[1] + a[2] * b[2]
        product > 0.0 || continue
        bisector = (a .+ b) ./ norm(a .+ b)
        in_metal_sector(side, bisector) || continue
        phi = acos(clamp(product, -1.0, 1.0))
        half_sine, half_cosine = sin(phi / 2), cos(phi / 2)
        point = (corner[1], corner[2])
        # The tube start stations of the two sides at this corner (the whole span of an
        # untubed side), projected onto the bisector.
        clearances = Float64[]
        for segment in segments
            abs(segment.plane - corner[3]) <= tolerance || continue
            if norm(segment.start .- collect(point)) <= tolerance
                push!(clearances, segment.untubed ? segment.span : segment.s_start)
            elseif norm(segment.stop .- collect(point)) <= tolerance
                push!(clearances, segment.untubed ? segment.span : segment.span - segment.s_end)
            end
        end
        length(clearances) == 2 ||
            error("Thin tip $(corner) has $(length(clearances)) metal sides in its plane (two expected)")
        tube_station = minimum(clearances) / half_cosine
        if phi < minimum_opening
            push!(unrefined, Dict{String, Any}(
                "Corner" => k - 1, "Point" => collect(corner), "OpeningDegrees" => rad2deg(phi),
                "TubeStation" => tube_station,
                "SplitFloor" => sector_regular_floor(phi / 2),
                "MinimumOpeningDegrees" => rad2deg(minimum_opening)))
            continue
        end
        half_width_station = nothing
        d = TIP_BISECTOR_STATION_STEP
        while d <= tube_station
            if d * half_sine >= size_law(d)
                half_width_station = d
                break
            end
            d += TIP_BISECTOR_STATION_STEP
        end
        length_b = half_width_station === nothing ? tube_station : half_width_station
        stop = (point[1] + length_b * bisector[1], point[2] + length_b * bisector[2])
        # The curve must clear every tube envelope of the plane (its own two tubes start
        # beyond the station by construction: the corner clearance R / tan(phi / 2) +
        # the envelope margin is where the two envelopes separate on the bisector).
        for segment in segments
            abs(segment.plane - corner[3]) <= tolerance && !segment.untubed || continue
            axis_start = segment.start .+ segment.s_start .* segment.direction
            axis_stop = segment.start .+ segment.s_end .* segment.direction
            distance = segment_segment_distance_2d(point, stop, axis_start, axis_stop)
            distance >= radius_envelope * (1.0 - 1.0e-9) ||
                error("the tip bisector curve of the thin tip $(corner) (length $(length_b)) comes " *
                      "within $(distance) of the tube axis $(segment.start) -> $(segment.stop) " *
                      "(envelope $(radius_envelope)): TipBisectorTube")
        end
        p0 = occ.addPoint(corner[1], corner[2], corner[3])
        p1 = occ.addPoint(stop[1], stop[2], corner[3])
        push!(tools, (Int32(1), occ.addLine(p0, p1)))
        push!(records, Dict{String, Any}(
            "Corner" => k - 1, "Point" => collect(corner), "OpeningDegrees" => rad2deg(phi),
            "Direction" => collect(bisector), "Length" => length_b, "End" => [stop[1], stop[2], corner[3]],
            "HalfWidthStation" => half_width_station, "TubeStation" => tube_station,
            "StationRule" => half_width_station === nothing ? "TubeStart" : "HalfWidthReachesCornerLaw",
            "SplitFloor" => sector_regular_floor(phi / 2),
            "MinimumOpeningDegrees" => rad2deg(minimum_opening),
            "Curves" => Int[], "Nodes" => 0, "EmbeddedIn" => Int[]))
    end
    return tools, records, unrefined
end

# After the fragment: the curve pieces of every bisector tool (domain_map entries) and,
# after the mesh, their node counts and the sheet surfaces they are embedded in (fail
# closed when a piece is embedded in no surface of the tip's plane).
function bind_tip_bisector_curves!(records, tool_map, tolerance)
    for (record, descendants) in zip(records, tool_map)
        curves = [Int(tag) for (dim, tag) in descendants if dim == 1]
        isempty(curves) && error("the tip bisector curve of $(record["Point"]) vanished in the fragment")
        record["Curves"] = curves
    end
    return records
end

function census_tip_bisector_curves!(records, tolerance)
    for record in records
        nodes = 0
        embedded = Int[]
        for curve in record["Curves"]
            tags, _, _ = gmsh.model.mesh.getNodes(1, curve, true)
            nodes = max(nodes, length(tags))
            for (dim, surface) in gmsh.model.getEntities(2)
                (1, curve) in [(Int(d), Int(t)) for (d, t) in gmsh.model.mesh.getEmbedded(2, surface)] || continue
                # The OCC bounding box carries the CAD's absolute tolerance (1e-7).
                _, _, zmin, _, _, zmax = gmsh.model.getBoundingBox(2, surface)
                plane_tolerance = max(tolerance, 1.0e-6)
                abs(zmin - record["Point"][3]) <= plane_tolerance &&
                    abs(zmax - record["Point"][3]) <= plane_tolerance ||
                    error("the tip bisector curve $(curve) of $(record["Point"]) is embedded in the " *
                          "surface $(surface) outside the sheet plane")
                push!(embedded, Int(surface))
            end
        end
        isempty(embedded) &&
            error("the tip bisector curve of $(record["Point"]) is embedded in no sheet surface")
        record["Nodes"] = nodes
        record["EmbeddedIn"] = sort!(unique!(embedded))
        nodes >= 2 || error("the tip bisector curve of $(record["Point"]) carries $(nodes) nodes")
    end
    return records
end

# The pinched-seam census of the thin sheet (Palace's crack opener, geodata.cpp
# refine_crack_elements): over the sheet faces (the 4000-family surfaces; triangles and
# quadrangles), an edge shared by two faces whose two vertices both lie on the sheet's
# FREE boundary (an edge of one face, off the outer box) is a pinched seam. Recorded
# ThinSheetSeams {Count, Edges}; a thin coupon fails closed on a non-zero count.
function thin_sheet_seam_census(sheet_surfaces, lower, upper, tolerance)
    faces = Vector{Vector{Int}}()
    for surface in sheet_surfaces
        types, _, blocks = gmsh.model.mesh.getElements(2, surface)
        for (type, block) in zip(types, blocks)
            name, _, _, nodes_per_element, _, _ = gmsh.model.mesh.getElementProperties(type)
            (startswith(name, "Triangle") || startswith(name, "Quadrilateral")) || continue
            corners = startswith(name, "Triangle") ? 3 : 4
            for start in 1:nodes_per_element:length(block)
                push!(faces, [Int(block[start + i - 1]) for i in 1:corners])
            end
        end
    end
    counts = Dict{Tuple{Int, Int}, Int}()
    for face in faces, i in eachindex(face)
        edge = minmax(face[i], face[mod1(i + 1, length(face))])
        counts[edge] = get(counts, edge, 0) + 1
    end
    coordinates = Dict{Int, Vector{Float64}}()
    for node in unique(reduce(vcat, faces; init=Int[]))
        coordinates[node] = gmsh.model.mesh.getNode(node)[1]
    end
    on_box(node) = any(abs(coordinates[node][d] - lower[d]) <= tolerance ||
                       abs(coordinates[node][d] - upper[d]) <= tolerance for d in 1:2)
    free = Set{Int}()
    for (edge, count) in counts
        count == 1 || continue
        for node in edge
            on_box(node) || push!(free, node)
        end
    end
    seams = [edge for (edge, count) in counts if count == 2 && all(in(free), edge)]
    return Dict{String, Any}(
        "Rule" => "a sheet edge shared by two sheet faces whose two vertices both lie on the " *
                  "sheet's free boundary (off the outer box) is a pinched seam: Palace's crack " *
                  "opener (geodata.cpp refine_crack_elements) must split it conformally, which " *
                  "a prism / pyramid mesh refuses; every thin coupon requires Count 0 (mesher " *
                  "design round 2 SEAM, decision 363)",
        "SheetSurfaces" => length(sheet_surfaces), "SheetFaces" => length(faces),
        "FreeBoundaryVertices" => length(free), "Count" => length(seams),
        "Edges" => [[coordinates[edge[1]], coordinates[edge[2]]] for edge in seams])
end

# Narrow metal strips (supervisor decision 347; family 6 of the mesher generality list):
# the tube of a metal side reaches Radius + PyramidHeight over the metal (the thin
# sheet tube's metal half, the fabricated top tube's quadrant above the top face and the
# bottom tube's quadrant under the metal), so two tubed sides of one plane whose tube
# intervals face each other ACROSS THE METAL (each lies on the metal side of the
# other's outward normal, at the closest points of the intervals) must be more than
# twice that envelope apart, like the sides of a hole (NarrowHoles). Returns the
# smallest facing width found (Inf without a facing pair); the guard fails closed
# with the measured width and the reach. Adjacent sides (sharing an end) meet at a
# corner, where the clearance rule applies; untubed short sides carry no tube.
function metal_facing_width(segments, envelope_radius, tolerance)
    width = Inf
    n = length(segments)
    for i in 1:n, j in (i + 1):n
        first, second = segments[i], segments[j]
        abs(first.plane - second.plane) <= tolerance || continue
        (first.untubed || second.untubed) && continue
        # Arc sides are not in the facing test (their interval is not a plan-view segment;
        # a narrow curved strip fails closed at the fragment, as every tube overlap does).
        (first.kind == :arc || second.kind == :arc) && continue
        shared = any(norm(p .- q) <= tolerance for p in (first.start, first.stop)
                     for q in (second.start, second.stop))
        shared && continue
        a = first.start .+ first.s_start .* first.direction
        b = first.start .+ first.s_end .* first.direction
        c = second.start .+ second.s_start .* second.direction
        d = second.start .+ second.s_end .* second.direction
        p, q = segment_closest_points_2d(a, b, c, d)
        separation = [q[1] - p[1], q[2] - p[2]]
        distance = norm(separation)
        distance > tolerance || continue
        # Across the metal: the other interval lies against each outward normal.
        dot(first.normal, separation) < -tolerance && dot(second.normal, separation) > tolerance ||
            continue
        if distance < width
            width = distance
            width <= 2.0 * envelope_radius &&
                scope_error("NarrowMetal", "the tubed metal sides $(first.start) -> $(first.stop) " *
                                           "and $(second.start) -> $(second.stop) on plane " *
                                           "$(first.plane) face each other across $distance of " *
                                           "metal against twice the tube envelope " *
                                           "$(envelope_radius) (Radius + PyramidHeight)")
        end
    end
    return width
end

# Width of a hole for its facing tubes: the smallest distance between two
# non-adjacent sides of the loop (adjacent sides meet at a corner, where the corner
# clearance rule applies).
function hole_facing_width(points, tolerance)
    n = length(points)
    n >= 4 || return Inf
    width = Inf
    for i in 1:n, j in (i + 2):n
        (i == 1 && j == n) && continue
        width = min(width, segment_segment_distance_2d(points[i], points[i % n + 1],
                                                       points[j], points[j % n + 1]))
    end
    return width
end

# The etch footprint must carry the metal edge (trench wall under the sidewall) so
# that the bottom tube's substrate sectors and vacuum sectors split exactly at the
# wall; the fragment then keeps every tube entity whole (match_tube_entities is the
# fail-closed check for producer-default footprints).
function assert_etch_carries_edge(etch_loops, segment, tolerance)
    function carried(start, stop)
        for loop in etch_loops
            n = length(loop.points)
            for i in 1:n
                p = loop.points[i]
                q = loop.points[i % n + 1]
                for (a, b) in ((p, q), (q, p))
                    norm([a[1], a[2]] .- start) <= tolerance &&
                        norm([b[1], b[2]] .- stop) <= tolerance && return true
                end
            end
        end
        return false
    end
    if haskey(segment, :arc) && segment.arc !== nothing
        # An arc side (design A7 MINOR-2): the footprint carries the whole metal ARC - every
        # chord of the tagged run (the trench wall under the arc is the same cylinder).
        for (a, b) in segment.arc.chords
            carried(a, b) || scope_error("FootprintWithoutEdge",
                                         "etch footprint has no side coincident with the chord $a -> $b of " *
                                         "the metal arc $(segment.arc.id)")
        end
        return true
    end
    carried(segment.start, segment.stop) && return true
    scope_error("FootprintWithoutEdge",
                "etch footprint has no side coincident with the metal edge $(segment.start) -> " *
                "$(segment.stop)")
end

# Ring count of the tube: the largest K with r_K + h_K <= the smallest transverse
# bound (trench depth below the bottom edge, metal thickness between the two
# edges, corner isotropy radius), so the pyramid apexes (h_K / 2 outside the
# tube) stay one outer ring size away from the trench floor and the two tubes of
# one sidewall never meet.
function tube_ring_count(edge_size, ratio, bound)
    rings = 0
    while true
        next = rings + 1
        radius = edge_size * (ratio^next - 1.0) / (ratio - 1.0)
        size = edge_size * ratio^(next - 1)
        radius + size <= bound || break
        rings = next
    end
    rings >= 1 || scope_error("NarrowTransverseBound",
                              "the tube inner size $edge_size does not fit the transverse bound $bound")
    return rings
end

const TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING = 0.5

# The OCC tube volumes of every straight metal edge of every upward process
# layer (top edge at plane + thickness, bottom edge at the plane), before the
# fragment. Returns the tool list, the per-volume records, the (tube, section)
# pairs, the edge segments and the section description for the census.
function build_edge_tubes!(occ, layers, loops, etch_loops, corners, edge_size, ratio,
                           sector_degrees, metal_thickness, overetch, corner_radius, lc_tangent,
                           lc_fine, lower, upper, tolerance; fabricated::Bool=true,
                           maximum_jacobian_condition::Float64=0.0)
    fabricated && (overetch > 0.0 || scope_error("NoTrench", "Overetch $overetch"))
    sectors = round(Int, 270.0 / sector_degrees)
    abs(sectors * sector_degrees - 270.0) <= 1.0e-9 || error("Tube sector angle must divide 270 degrees")
    per_quadrant = round(Int, 90.0 / sector_degrees)
    abs(per_quadrant * sector_degrees - 90.0) <= 1.0e-9 || error("Tube sector angle must divide 90 degrees")
    # Thin metal (decision 66): the sheet has no thickness and no trench, so the only
    # transverse bound is the corner isotropy radius (the facing rules below still hold).
    bound = fabricated ? min(overetch, 0.5 * metal_thickness, corner_radius) : corner_radius
    rings = tube_ring_count(edge_size, ratio, bound)
    # Top edge: vacuum from the sidewall ray (-90) through the outward normal (0)
    # and up (90) to the top-face ray (180). Bottom edge: substrate from the
    # metal bottom face ray (180) to the trench wall ray (270), vacuum from the
    # trench wall past the outward normal (360) to the sidewall ray (450).
    # Thin sheet edge: one full-turn section split by the sheet plane on the whole
    # diameter - substrate from the sheet's substrate side (180) down (270) to the SA
    # floor ray (360), vacuum from the SA floor up (450) to the sheet's vacuum side
    # (540 = 180: a closed section, the sheet ray is both bounding rays).
    top_section = TubeSection(edge_size, ratio, rings,
                              [-90.0 + sector_degrees * j for j in 0:sectors], fill(2, sectors))
    bottom_section = TubeSection(edge_size, ratio, rings,
                                 [180.0 + sector_degrees * j for j in 0:sectors],
                                 vcat(fill(1, per_quadrant), fill(2, sectors - per_quadrant)))
    full_turn = 4 * per_quadrant
    sheet_section = TubeSection(edge_size, ratio, rings,
                                [180.0 + sector_degrees * j for j in 0:full_turn],
                                vcat(fill(1, 2 * per_quadrant), fill(2, 2 * per_quadrant)))
    radius = tube_radius(top_section)
    outer_ring = ring_sizes(top_section)[end]
    pyramid_height = TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING * outer_ring
    # The face-end end-spacing cap of the section (design round 2 F2b 3.2): derived from
    # the section's own prism frames and the Jacobian-condition ceiling, needed only where
    # a face end exists (a coupon without one never reads it: regime I is bitwise and a
    # theta-0 end has no FaceEnd); the ceiling is required then (fail closed).
    face_end_cap = Ref{Union{Nothing, Float64}}(nothing)
    function face_end_spacing_cap_of_coupon()
        if face_end_cap[] === nothing
            maximum_jacobian_condition > 1.0 ||
                error("a face end needs the Jacobian-condition ceiling " *
                      "(--maximum-jacobian-condition) to derive its end-spacing cap")
            face_end_cap[] = face_end_spacing_cap(fabricated ? top_section : sheet_section,
                                                  maximum_jacobian_condition)
        end
        return face_end_cap[]
    end
    # The validity ceiling of the capped end block (design 3.2): beyond it no sheared layer
    # keeps its pyramid apex inside the box within the cap - fail closed.
    function guard_steep_face_crossing(face)
        cap = face_end_spacing_cap_of_coupon()
        apex = 2.0 * pyramid_height * abs(tan(face.theta))
        apex < cap ||
            scope_error("SteepFaceCrossing",
                        "the tube end on face $(face.face) at $(rad2deg(face.theta)) degrees " *
                        "needs a thinnest layer 2 PyramidHeight |tan theta| = $apex at or above " *
                        "the end-spacing cap $cap of the tube section (the validity ceiling is " *
                        "$(rad2deg(atan(cap / (2.0 * pyramid_height)))) degrees)")
        return cap
    end
    fabricated && (radius + pyramid_height < overetch ||
        scope_error("ShallowTrench", "tube radius $radius + pyramid height $pyramid_height " *
                                     "against Overetch $overetch"))
    # The corner clearance (design A3 (3), decision 303): the tube solids retract by
    # R / tan(phi / 2) so that their sidewalls meet at the corner bisector, plus the
    # larger of the outer ring h_K and the pyramid-envelope margin 1.25 h_pyr /
    # sin(phi / 2) (the apexes h_pyr outside the two tube surfaces need 2 h_pyr + a
    # h_pyr / 2 margin of plan-view separation 2 (c - R / tan(phi / 2)) sin(phi / 2)).
    # The max() is h_K exactly wherever sin(phi / 2) >= 0.625 (phi >= 77.4 degrees:
    # every rectilinear corner, a single tube edge at phi = pi), so those stay bitwise.
    clearance(angle) = (angle < pi - 1.0e-9 ? radius / tan(0.5 * angle) : 0.0) +
                       max(outer_ring, 1.25 * pyramid_height / sin(0.5 * angle))
    envelope_radius = radius + pyramid_height
    # Facing tubes (the sides of a hole carry tubes pointing into it): the hole must be
    # wider than two tube reaches, a reach being the tube radius, the pyramid height and
    # the band's protected distance, so that no two tube bands overlap across it.
    facing_reach = radius + pyramid_height + BAND_PROTECTED_DISTANCE_OVER_NORMAL * lc_fine
    for loop in loops
        loop.hole || continue
        width = hole_facing_width(loop.points, tolerance)
        width > 2.0 * facing_reach ||
            scope_error("NarrowHoles", "hole of conductor $(loop.conductor) on plane " *
                                       "$(loop.plane) is $width wide against twice the tube " *
                                       "reach $facing_reach")
    end
    tools = Tuple{Int32, Int32}[]
    records = TubeRecord[]
    tubes = Tuple{AbstractTube, TubeSection}[]
    segments = NamedTuple[]
    # The smooth joints (design A3 (2)) as (adopting tube index, its end, owning tube index,
    # its end) into `tubes`, resolved after every tube exists: an arc part owns the section
    # shared with a straight neighbour; the earlier part of one arc owns the split section.
    joints = Tuple{Int, Int, Int, Int}[]
    arc_tubes = 0
    shared_sections = 0
    # The untubed short sides (design A9 family 3): recorded, no tube, not in `segments`
    # (which stays in lock-step with `tubes`: TubesPerSide tubes per segment).
    untubed_edges = Dict{String, Any}[]
    # The narrowest metal strip between facing tubed sides (decision 347; Inf without one).
    metal_facing = Inf
    # Facing process layers (an upward layer below a downward one, the only pair
    # layer_groups admits): the vacuum gap between their metal top faces must exceed
    # two tube reaches, like the width of a hole.
    face_offset = fabricated ? metal_thickness : 0.0
    for lower_layer in layers, upper_layer in layers
        lower_layer.sign > 0 && upper_layer.sign < 0 || continue
        gap = (upper_layer.plane - face_offset) - (lower_layer.plane + face_offset)
        gap > 2.0 * facing_reach ||
            scope_error("NarrowLayerGap", "vacuum gap $gap between the metal faces of the layers " *
                                          "at z = $(lower_layer.plane) (Nz = 1) and z = " *
                                          "$(upper_layer.plane) (Nz = -1) against twice the tube " *
                                          "reach $facing_reach")
    end
    for layer in layers
        layer_loops = [loop for loop in loops if abs(loop.plane - layer.plane) <= tolerance]
        isempty(layer_loops) && error("Plan-view boundary is missing the tube layer $(layer.plane)")
        layer_segments = [(segment..., layer_sign=layer.sign) for segment in
                          metal_edge_segments(layer_loops, corners, clearance, lower, upper, tolerance;
                                              edge_size=edge_size, corner_radius=corner_radius)]
        if etch_loops !== nothing
            for segment in layer_segments
                assert_etch_carries_edge(etch_loops, segment, tolerance)
            end
        end
        metal_facing = min(metal_facing, metal_facing_width(layer_segments, envelope_radius, tolerance))
        # The tubes of this layer's segments in segment order (TubesPerSide per segment); the
        # joints are resolved on the indices into `tubes` once every segment has its tubes.
        layer_tube_index = Dict{Int, Int}()
        tubed_segments = Tuple{Int, NamedTuple}[]
        for (segment_index, segment) in enumerate(layer_segments)
            if segment.untubed
                push!(untubed_edges, untubed_edge_record(segment))
                continue
            end
            b = [0.0, 0.0, Float64(layer.sign)]
            placements = fabricated ?
                ((layer.plane + layer.sign * metal_thickness, top_section),
                 (layer.plane, bottom_section)) :
                ((layer.plane, sheet_section),)
            layer_tube_index[segment_index] = length(tubes)
            if segment.kind == :arc
                # An ARC part (design 1.2 (3)): the tube revolves about the arc centre; its
                # travel e = n x b runs with the loop (orientation -sigma b_z = travel) or
                # against it, as the straight tube's `along` does.
                arc = segment.arc
                envelope_radius <= 0.25 * arc.rho ||
                    scope_error("ArcTubeRadiusVsCurvature",
                                "arc $(arc.id) of radius $(arc.rho) against the tube envelope " *
                                "$(envelope_radius) (Radius + PyramidHeight)")
                orientation = -arc.sigma * layer.sign
                travel = sign(arc.sweep)
                along = orientation * travel
                theta0 = along > 0.0 ? arc.theta_start : arc.theta_end
                face_ends = FaceEnd[]
                for (segment_end, face) in enumerate(segment.face_ends)
                    face === nothing && continue
                    end_index = (along > 0.0) == (segment_end == 1) ? 0 : 1
                    normal = [face.normal[1], face.normal[2], 0.0]
                    push!(face_ends, FaceEnd(end_index, face.face, normal, face.theta, 0.0, 0.0,
                                             envelope_radius, pyramid_height, lc_tangent;
                                             spacing_cap=guard_steep_face_crossing(face),
                                             face_axis=face.axis, face_value=face.value))
                end
                for (z, section) in placements
                    s_start, s_end = along > 0.0 ? (segment.s_start, segment.s_end) :
                                     (segment.span - segment.s_end, segment.span - segment.s_start)
                    tube = ArcTube([arc.centre[1], arc.centre[2], z], arc.rho, arc.sigma, b, theta0,
                                   s_start, s_end, lc_tangent; face_ends=face_ends)
                    push!(tubes, (tube, section))
                    arc_tubes += 1
                    for (tool, group) in add_tube_volumes!(occ, tube, section; box=(lower, upper))
                        push!(records, TubeRecord(tube, section, group, tool))
                        push!(tools, tool)
                    end
                end
            else
                # The tube frame follows the layer: b is the process normal (Nz), so the
                # sections' "up" (towards the metal top face) and the extrusion sense
                # e = n x b mirror for a downward layer; the top tube sits on the metal top
                # face at plane + Nz x thickness, the bottom tube on the plane.
                n = [segment.normal[1], segment.normal[2], 0.0]
                e = cross(n, b)
                along = dot(e[1:2], segment.direction)
                abs(abs(along) - 1.0) <= 1.0e-12 || error("Tube frame is not aligned with the edge")
                # Face ends (design A2): the segment's start / stop face ends map onto the
                # tube's start (end index 0) / end (1) through the extrusion sense; the face
                # plane's slope in tube coordinates is kappa = (N . n, N . b) / (N . e).
                face_ends = FaceEnd[]
                for (segment_end, face) in enumerate(segment.face_ends)
                    face === nothing && continue
                    end_index = (along > 0.0) == (segment_end == 1) ? 0 : 1
                    normal = [face.normal[1], face.normal[2], 0.0]
                    along_normal = dot(normal, e)
                    abs(along_normal) > 0.0 || error("Tube axis lies in the box face $(face.face)")
                    push!(face_ends, FaceEnd(end_index, face.face, normal, face.theta,
                                             dot(normal, n) / along_normal, dot(normal, b) / along_normal,
                                             envelope_radius, pyramid_height, lc_tangent;
                                             spacing_cap=guard_steep_face_crossing(face),
                                             face_axis=face.axis, face_value=face.value))
                end
                for (z, section) in placements
                    # Smooth joints (design A3 (2)): the straight tube's end section at a joint
                    # with an arc part is the arc's radial section there (the arc's frame at the
                    # joint vertex; the tilt = the post-snap kink between the two tangents).
                    joint_ends = JointEnd[]
                    for (segment_end, joint) in enumerate(segment.joints)
                        joint === nothing && continue
                        other = layer_segments[joint[1]]
                        other.kind == :arc || continue
                        end_index = (along > 0.0) == (segment_end == 1) ? 0 : 1
                        theta = joint[2] == 1 ? other.arc.theta_start : other.arc.theta_end
                        radial = [cos(theta), sin(theta), 0.0]
                        n_arc = other.arc.sigma .* radial
                        origin = [other.arc.centre[1] + other.arc.rho * radial[1],
                                  other.arc.centre[2] + other.arc.rho * radial[2], z]
                        # The arc's travel tangent at the joint, pointed OUT of the straight tube.
                        tangent = [-radial[2], radial[1], 0.0]
                        outward = end_index == 1 ? e : -e
                        dot(tangent, outward) < 0.0 && (tangent = -tangent)
                        tilt = acos(clamp(dot(tangent, outward), -1.0, 1.0))
                        push!(joint_ends, JointEnd(end_index, origin, n_arc, b, tangent, tilt,
                                                   envelope_radius, lc_tangent))
                    end
                    tube = if along > 0.0
                        EdgeTube([segment.start[1], segment.start[2], z], n, b, segment.s_start,
                                 segment.s_end, lc_tangent; face_ends=face_ends, joints=joint_ends)
                    else
                        EdgeTube([segment.stop[1], segment.stop[2], z], n, b,
                                 segment.span - segment.s_end, segment.span - segment.s_start, lc_tangent;
                                 face_ends=face_ends, joints=joint_ends)
                    end
                    push!(tubes, (tube, section))
                    for (tool, group) in add_tube_volumes!(occ, tube, section; box=(lower, upper))
                        push!(records, TubeRecord(tube, section, group, tool))
                        push!(tools, tool)
                    end
                end
            end
            push!(segments, segment)
            push!(tubed_segments, (segment_index, segment))
        end
        # The joint table of this layer: for every smooth joint the ADOPTING tube (a straight
        # tube, or the later part of one arc) and the OWNING tube (the arc part), per
        # placement (top / bottom / sheet), with the tube ends through each tube's travel.
        tube_end(tube::AbstractTube, point, z) = begin
            start = tube_point(tube, 0.0, 0.0, tube.s_start)
            stop = tube_point(tube, 0.0, 0.0, tube.s_end)
            # The joint vertex is the CAD end of the tube interval (clearance 0 there).
            d_start = norm(start .- [point[1], point[2], z])
            d_stop = norm(stop .- [point[1], point[2], z])
            d_start <= tolerance && d_stop <= tolerance && error("a tube of zero length at a joint")
            d_start <= tolerance ? 0 : d_stop <= tolerance ? 1 :
                error("the tube does not end at the joint vertex $point (ends at $start / $stop)")
        end
        placements_count = fabricated ? 2 : 1
        for (segment_index, segment) in tubed_segments
            for (segment_end, joint) in enumerate(segment.joints)
                joint === nothing && continue
                other_index, _ = joint
                other = layer_segments[other_index]
                other.untubed && error("a smooth joint with an untubed side")
                point = segment_end == 1 ? segment.start : segment.stop
                # Owner: the arc part; between two parts of one arc, the earlier part.
                owner_is_other = other.kind == :arc && (segment.kind != :arc || other.arc.part < segment.arc.part)
                owner_is_other || continue      # recorded from the adopting side only
                for placement in 1:placements_count
                    adopting = layer_tube_index[segment_index] + placement
                    owning = layer_tube_index[other_index] + placement
                    z = tubes[adopting][1] isa ArcTube ? tubes[adopting][1].centre[3] : tubes[adopting][1].origin[3]
                    push!(joints, (adopting, tube_end(tubes[adopting][1], point, z),
                                   owning, tube_end(tubes[owning][1], point, z)))
                    shared_sections += 1
                end
            end
        end
    end
    description = Dict{String, Any}(
        "Rings" => rings, "RingSizes" => ring_sizes(top_section),
        "RingRadii" => copy(top_section.ring_radii), "Radius" => radius,
        "TransverseBound" => bound,
        "RingRule" => fabricated ?
                      "largest K with r_K + h_K <= min(Overetch, MetalThickness / 2, " *
                      "CornerIsotropyRadius)" :
                      "largest K with r_K + h_K <= CornerIsotropyRadius (thin sheet: no " *
                      "trench and no thickness bound the tube)",
        "Kind" => fabricated ? "fabricated" : "thin",
        "TubesPerSide" => fabricated ? 2 : 1,
        "TubesPerSideRule" => fabricated ?
                              "a top tube and a bottom tube per straight metal side" :
                              "one sheet tube per straight metal side at the process plane " *
                              "(decision 66: the thin sheet edge is one line where the metal " *
                              "sheet, the SA floor and both dielectrics meet)",
        # The thin coupon is recorded at its cutoff, never converged (decision 66): the
        # inner ring size is the MA / MS cutoff of every thin surface participation.
        "ThinCutoff" => fabricated ? nothing : edge_size,
        "SectorDegrees" => sector_degrees, "Sectors" => sectors,
        "PyramidHeight" => pyramid_height,
        "PyramidHeightOverOuterRing" => TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING,
        "CornerClearanceRule" => "R / tan(phi / 2) + max(h_K, 1.25 h_pyr / sin(phi / 2)) before " *
                                 "a semantic corner, phi the smallest in-plane angle between the " *
                                 "tube edges meeting there (h_K alone for a single tube edge; the " *
                                 "max() is h_K for phi >= 77.4 degrees, the pyramid-envelope margin " *
                                 "below); 0 at box continuation vertices and at face ends",
        "FaceEndRule" => "a tube end on a box face with a single metal side at the vertex and a " *
                         "tilt theta > 0 (the edge's face-parallel coordinate differs between its " *
                         "two ends: exact arithmetic, no tolerance) ends ON the face: the CAD " *
                         "solid is extruded over-long by (R + h_pyr) |tan theta| + TangentialSize " *
                         "and intersected with the coupon box before the fragment; the mesh ends " *
                         "with m sheared layers of axial spacing lc_end whose last station is the " *
                         "face plane (Tubes[].FaceEnds). Regime I (4 h_pyr |tan theta| <= lc_cap): " *
                         "lc_end = max(TangentialSize, 4 h_pyr |tan theta|), m = ceil(2 (R + h_pyr) " *
                         "|tan theta| / lc_end) (block (b) design A2 (4), bitwise); regime II (4 h_pyr " *
                         "|tan theta| > lc_cap): lc_end = lc_cap, m = max(ceil((R + h_pyr) |tan theta| / " *
                         "(lc_cap - 2 h_pyr |tan theta|)), the regime-I count at lc_cap), so the " *
                         "thinnest layer keeps the apex rule t_min >= 2 h_pyr |tan theta| (mesher " *
                         "design round 2 F2b 3.2); lc_cap = FaceEndSpacingCap (FaceEndSpacingCapRule); " *
                         "2 h_pyr |tan theta| >= lc_cap fails closed at ScopeGuard[SteepFaceCrossing]; " *
                         "theta == 0 keeps the unchanged perpendicular end (decisions 302 / 320)",
        "FaceEndSpacingCapRule" => "the largest axial spacing at which the section's own prism " *
                                   "corner frames (every sector triangle of the section at each " *
                                   "vertex, the axial edge orthogonal: singular values {lc, " *
                                   "sigma_1, sigma_2} of the frame) read a Jacobian condition " *
                                   "<= $(FACE_END_CONDITION_MARGIN) x MaximumJacobianCondition: " *
                                   "lc_cap = $(FACE_END_CONDITION_MARGIN) x MaximumJacobianCondition " *
                                   "x min sigma_2 (closed form; the production sections' smallest " *
                                   "sigma_2 is the second ring's (a, c, d) triangle at its ring-2 " *
                                   "vertex); null on a coupon without a face end",
        "FaceEndSpacingCap" => face_end_cap[],
        "BoxVertexRule" => "a box-face vertex with a single metal side is a cut end (FaceEnd, no " *
                           "ball, clearance 0) when theta > 0; at theta == 0 exactly the legacy " *
                           "convention holds bitwise: a Physical-class box vertex is a semantic " *
                           "corner with h_K + its ball (Tubes[].LegacyBoxVertexCorner), a " *
                           "Continuation-class one has clearance 0; two metal sides meeting at a " *
                           "box vertex are a corner at any theta (decision 320). An ARC end exactly " *
                           "perpendicular to a box face (|tangent[d]| == 1.0) is a cut end for the " *
                           "contract (box_face_cut_end reads the next chord vertex's face-parallel " *
                           "coordinate) and the legacy perpendicular end for the mesher (clearance " *
                           "0, no ball, no disagreement error since the contract lists no corner " *
                           "there): consistent in effect, recorded here (mesher review R2 MINOR-5, " *
                           "decision 391); every arc end on the outer box fails closed at " *
                           "ScopeGuard[ArcFaceEnds] until the next mesher lane lifts it",
        "EnvelopeRadius" => envelope_radius,
        "MetalFacingRule" => "two tubed metal sides of one plane whose tube intervals face each " *
                             "other across the metal (each on the metal side of the other's " *
                             "outward normal at their closest points) are more than 2 x " *
                             "EnvelopeRadius = 2 (Radius + PyramidHeight) apart, else the build " *
                             "fails closed at ScopeGuard[NarrowMetal] with the measured width " *
                             "(supervisor decision 347)",
        "MetalFacingWidth" => isfinite(metal_facing) ? metal_facing : nothing,
        "FacingReach" => facing_reach,
        "FacingRule" => "every hole is wider than 2 x (Radius + PyramidHeight + " *
                        "ProtectedDistance) between any two of its non-adjacent sides, and the " *
                        "vacuum gap between the metal top faces of an upward and a downward " *
                        "process layer exceeds the same 2 x reach, so the tubes facing each " *
                        "other across a hole or a layer gap keep disjoint bands",
        "TubeFrameRule" => "b = (0, 0, Nz) of the tube's process layer, n the horizontal " *
                           "normal away from the metal, e = n x b; the top tube lies on the " *
                           "metal top face at plane + Nz x MetalThickness, the bottom tube on " *
                           "the plane (Tubes[].Layer records Nz); a thin sheet tube lies on the plane",
        "Top" => fabricated ? Dict("Angles" => top_section.angles, "Materials" => top_section.materials) : nothing,
        "Bottom" => fabricated ? Dict("Angles" => bottom_section.angles, "Materials" => bottom_section.materials) : nothing,
        "Sheet" => fabricated ? nothing :
                   Dict("Angles" => sheet_section.angles, "Materials" => sheet_section.materials,
                        "Closed" => sheet_section.closed))
    if arc_tubes > 0
        # Recorded only where arcs exist (block (b) design 1.2 (3) / A3 (2)), so every straight
        # census is unchanged.
        description["ArcTubes"] = arc_tubes
        description["SharedSections"] = shared_sections
        description["ArcTubeRule"] = ARC_TUBE_RULE
        description["SmoothJointRule"] = SMOOTH_JOINT_RULE
        description["ArcCurvatureBound"] = 0.25
    end
    return tools, records, tubes, segments, description, untubed_edges, joints
end

const ARC_TUBE_RULE =
    "block (b) design 1.2 (3) (decision 303): the tube of an arc metal side is the revolve of " *
    "the section (the same rings, sectors and pyramids as a straight tube) about the vertical " *
    "axis through the tagged arc centre; the axis coordinate is the arc length on the edge " *
    "circle, the layers follow the composed size field along it (LayerRule); the sidewall " *
    "rays are the cylinder of radius rho about the centre - the exact metal / trench wall - " *
    "so the fragment splits the sectors at the wall; an arc is split into ceil(sweep / 90 " *
    "degrees) parts of equal angular fraction (the same radial planes as the metal loft and " *
    "the trench wall of polygon_wire), consecutive parts sharing one section; the guard " *
    "ScopeGuard[ArcTubeRadiusVsCurvature] requires Radius + PyramidHeight <= ArcCurvatureBound x rho"
const SMOOTH_JOINT_RULE =
    "block (b) design A3 (1)-(2) (decision 303): a joint whose signature turn is at most " *
    "1e-4 rad (JointSmooth) is no corner - no ball, no caps, clearance 0 - and the two tubes " *
    "share ONE cross-section owned by the arc (its radial plane through the joint vertex, " *
    "its frame): the straight tube's end nodes are the arc tube's nodes, its CAD solid is " *
    "extruded over-long and cut by that radial plane when the post-snap tilt is positive " *
    "(Tubes[].Joints: TiltRadians, PlaneCut); a kinked joint is a corner under " *
    "CornerClearanceRule (the arc tube retracts by the clearance along its arc)"

# The census record of an untubed short side (design A9 family 3): the side, its
# span, the two corner clearances (0 at a box vertex or face end), the tube interval
# they leave, which ends are semantic corners and whether their balls cover the span.
function untubed_edge_record(segment)
    return Dict{String, Any}(
        "Side" => Dict{String, Any}(
            "Start" => copy(segment.start), "Stop" => copy(segment.stop),
            "Plane" => segment.plane, "Conductor" => segment.conductor,
            "Hole" => segment.hole, "Layer" => segment.layer_sign),
        "Span" => segment.span,
        "Clearances" => [segment.s_start, segment.span - segment.s_end],
        "Interval" => segment.s_end - segment.s_start,
        "Corners" => collect(segment.corners),
        "CoveredByBalls" => segment.covered_by_balls)
end

# Gmsh element type of the linear volume cells and their vertex counts.
const GMSH_LINEAR_VOLUME_TYPES = Dict(4 => ("Tetrahedron", 4), 6 => ("Prism", 6), 7 => ("Pyramid", 5))

# Corner frames of the linear volume cells: at each listed vertex the three edges
# to its neighbours span the local Jacobian (the same frames as the Python
# mixed-mesh audits). Tetrahedra keep the production vertex-0 frame.
const VOLUME_CORNER_FRAMES = Dict(
    4 => [(1, 2, 3, 4)],
    6 => [(1, 2, 3, 4), (2, 3, 1, 5), (3, 1, 2, 6), (4, 6, 5, 1), (5, 4, 6, 2), (6, 5, 4, 3)],
    7 => [(1, 2, 4, 5), (2, 3, 1, 5), (3, 4, 2, 5), (4, 1, 3, 5)])

# Node tags of every non-simplex volume cell (the tube prisms and pyramids).
function non_simplex_volume_nodes()
    tags = UInt[]
    types, _, element_nodes = gmsh.model.mesh.getElements(3)
    for (type, block) in zip(types, element_nodes)
        Int(type) == GMSH_LINEAR_ELEMENT_TYPE[3] && continue
        append!(tags, block)
    end
    return unique!(tags)
end

# Per-type quality of every linear volume cell in the model: orientation (every
# corner Jacobian determinant positive), Jacobian condition (largest / smallest
# singular value of the corner edge matrix, maximum over the corners) and, for
# tetrahedra, the scaled Jacobian of the vertex-0 frame (production definition).
function mixed_volume_quality()
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    types, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
    result = Dict{String, Any}()
    total = 0
    for (type, tags, block) in zip(types, element_tags, element_nodes)
        isempty(tags) && continue
        haskey(GMSH_LINEAR_VOLUME_TYPES, Int(type)) ||
            error("Unsupported volume element type $type")
        name, width = GMSH_LINEAR_VOLUME_TYPES[Int(type)]
        cell_count = length(tags)
        total += cell_count
        frames = VOLUME_CORNER_FRAMES[Int(type)]
        condition = zeros(cell_count)
        scaled = fill(Inf, cell_count)
        determinant = fill(Inf, cell_count)
        for n in 1:cell_count
            start = (n - 1) * width
            for (v, a, b, c) in frames
                p0 = points[:, index[block[start + v]]]
                e1 = points[:, index[block[start + a]]] .- p0
                e2 = points[:, index[block[start + b]]] .- p0
                e3 = points[:, index[block[start + c]]] .- p0
                jacobian = hcat(e1, e2, e3)
                det_value = det(jacobian)
                determinant[n] = min(determinant[n], det_value)
                singular = svdvals(jacobian)
                condition[n] = max(condition[n], singular[1] / max(singular[end], 1.0e-300))
                scaled[n] = min(scaled[n], abs(det_value) / max(norm(e1) * norm(e2) * norm(e3), 1.0e-300))
            end
        end
        quantile(values, q) = (sorted = sort(values);
                               sorted[clamp(ceil(Int, q * length(sorted)), 1, length(sorted))])
        result[name] = Dict{String, Any}(
            "Count" => cell_count, "GmshType" => Int(type),
            "PositiveOrientation" => all(>(0.0), determinant),
            "NonpositiveCells" => count(<=(0.0), determinant),
            "MinimumScaledJacobian" => minimum(scaled),
            "ScaledJacobianQuantiles" => [quantile(scaled, q) for q in (0.0, 0.01, 0.5, 1.0)],
            "MaximumJacobianCondition" => maximum(condition),
            "JacobianConditionQuantiles" => [quantile(condition, q) for q in (0.0, 0.5, 0.99, 1.0)],
            "CellsBelowScaledJacobian0.01" => count(<(0.01), scaled),
            "CellsAboveCondition1000" => count(>(1000.0), condition))
    end
    result["Total"] = total
    return result
end

# Gate the per-type quality (fail closed): positive orientation and the Jacobian
# condition bound for every type, the scaled-Jacobian bound for tetrahedra.
function gate_mixed_volume_quality(quality, minimum_scaled_jacobian, maximum_jacobian_condition)
    failures = String[]
    for (name, record) in quality
        record isa Dict || continue
        record["PositiveOrientation"] ||
            push!(failures, "$name: $(record["NonpositiveCells"]) nonpositive cells")
        record["MaximumJacobianCondition"] <= maximum_jacobian_condition ||
            push!(failures, "$name: Jacobian condition $(record["MaximumJacobianCondition"]) > $maximum_jacobian_condition")
        name == "Tetrahedron" && record["MinimumScaledJacobian"] < minimum_scaled_jacobian &&
            push!(failures, "Tetrahedron: scaled Jacobian $(record["MinimumScaledJacobian"]) < $minimum_scaled_jacobian")
    end
    isempty(failures) || error("Mixed-element quality gates failed: " * join(failures, "; "))
    return nothing
end

# Edge-length statistics of the surface elements of one physical label restricted
# to elements whose centroid lies within `reach` of a segment set (or all of them
# when the set is nothing): the achieved cut-surface / band sizes.
function surface_size_statistics(points, index, attribute, segments, reach)
    lengths = Float64[]
    elements = 0
    for entity in gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
        types, element_tags, element_nodes = gmsh.model.mesh.getElements(2, entity)
        for (type, tags, block) in zip(types, element_tags, element_nodes)
            isempty(tags) && continue
            width = length(block) ÷ length(tags)
            for start in 0:width:(length(block) - 1)
                nodes = [points[:, index[block[start + i]]] for i in 1:width]
                if segments !== nothing
                    centroid = sum(nodes) ./ width
                    minimum(span_point_distance(centroid, a, b) for (a, b) in segments) <= reach ||
                        continue
                end
                elements += 1
                for i in 1:width
                    push!(lengths, norm(nodes[i] .- nodes[i % width + 1]))
                end
            end
        end
    end
    isempty(lengths) && return Dict{String, Any}("Elements" => 0)
    sort!(lengths)
    at(q) = lengths[clamp(ceil(Int, q * length(lengths)), 1, length(lengths))]
    return Dict{String, Any}("Elements" => elements, "Edges" => length(lengths),
                             "EdgeMinimum" => lengths[1], "EdgeP10" => at(0.1),
                             "EdgeP50" => at(0.5), "EdgeP90" => at(0.9),
                             "EdgeMaximum" => lengths[end])
end

# Achieved cut-surface size against the trace rule: for every narrow basis
# triangle (requested size below the far size) the mean edge length of the
# matching-surface elements whose centroid lies in it, over the requested size.
function trace_basis_achieved_sizes(points, index, matching_attribute, basis_record)
    triangles = TRACE_BASIS_TRIANGLES
    isempty(triangles) && return nothing
    sums = zeros(length(triangles)); counts = zeros(Int, length(triangles))
    for entity in gmsh.model.getEntitiesForPhysicalGroup(2, matching_attribute)
        types, element_tags, element_nodes = gmsh.model.mesh.getElements(2, entity)
        for (type, tags, block) in zip(types, element_tags, element_nodes)
            isempty(tags) && continue
            width = length(block) ÷ length(tags)
            for start in 0:width:(length(block) - 1)
                nodes = [points[:, index[block[start + i]]] for i in 1:width]
                centroid = sum(nodes) ./ width
                mean_edge = sum(norm(nodes[i] .- nodes[i % width + 1]) for i in 1:width) / width
                for (k, t) in enumerate(triangles)
                    a = (t[1], t[2], t[3]); b = (t[4], t[5], t[6]); c = (t[7], t[8], t[9])
                    point_triangle_distance(Tuple(centroid), a, b, c) <= 1.0e-9 || continue
                    sums[k] += mean_edge; counts[k] += 1
                    break
                end
            end
        end
    end
    ratios = [counts[k] > 0 ? sums[k] / counts[k] / TRACE_BASIS_SIZES[k] : NaN
              for k in eachindex(triangles)]
    covered = filter(isfinite, ratios)
    return Dict{String, Any}(
        "NarrowBasisTriangles" => length(triangles),
        "TrianglesWithMatchingElements" => length(covered),
        "AchievedOverRequestedMaximum" => isempty(covered) ? nothing : maximum(covered),
        "AchievedOverRequestedP50" => isempty(covered) ? nothing : sorted_median(covered),
        "Rule" => "mean edge length of the matching-surface elements whose centroid lies " *
                  "in a narrow basis triangle over the triangle's requested size " *
                  "(TraceBasisSizeRatio x its minimum altitude)")
end

# Size statistics of the tetrahedra whose centroid lies within `reach` of a
# segment set: mean edge length percentiles (the junction / footprint band sizes).
function band_tetrahedron_statistics(points, index, segments, reach)
    isempty(segments) && return Dict{String, Any}("Cells" => 0, "Segments" => 0)
    sizes = Float64[]
    types, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
    for (type, tags, block) in zip(types, element_tags, element_nodes)
        Int(type) == 4 || continue
        for start in 0:4:(length(block) - 1)
            nodes = [points[:, index[block[start + i]]] for i in 1:4]
            centroid = sum(nodes) ./ 4
            minimum(span_point_distance(centroid, a, b) for (a, b) in segments) <= reach || continue
            push!(sizes, sum(norm(nodes[i] .- nodes[j]) for i in 1:4 for j in (i + 1):4) / 6)
        end
    end
    isempty(sizes) && return Dict{String, Any}("Cells" => 0, "Segments" => length(segments))
    sort!(sizes)
    at(q) = sizes[clamp(ceil(Int, q * length(sizes)), 1, length(sizes))]
    return Dict{String, Any}("Cells" => length(sizes), "Segments" => length(segments),
                             "Reach" => reach,
                             "MeanEdgeP10" => at(0.1), "MeanEdgeP50" => at(0.5),
                             "MeanEdgeP90" => at(0.9), "MeanEdgeMaximum" => sizes[end])
end

# Achieved volume size against a prescribed law in distance shells (decision 39):
# every tetrahedron whose centroid lies at distance `distance(centroid)` within
# the last shell bound is binned into the shells (lower bound `lower`, then the
# `bounds`); per shell the mean and longest edge percentiles and the ratio of the
# mean edge to `law(centroid, distance)` are reported. Nothing is gated.
function shell_size_statistics(points, index, distance, lower, bounds, law)
    shells = [Dict{String, Any}("Lower" => a, "Upper" => b, "Mean" => Float64[],
                                "Longest" => Float64[], "Ratio" => Float64[])
              for (a, b) in zip(vcat(lower, bounds[1:(end - 1)]), bounds)]
    types, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
    for (type, tags, block) in zip(types, element_tags, element_nodes)
        Int(type) == GMSH_LINEAR_ELEMENT_TYPE[3] || continue
        for start in 0:4:(length(block) - 1)
            nodes = [points[:, index[block[start + i]]] for i in 1:4]
            centroid = sum(nodes) ./ 4
            d = distance(centroid)
            lower <= d <= bounds[end] || continue
            shell = shells[something(findfirst(>=(d), bounds))]
            lengths = [norm(nodes[i] .- nodes[j]) for i in 1:4 for j in (i + 1):4]
            push!(shell["Mean"], sum(lengths) / 6)
            push!(shell["Longest"], maximum(lengths))
            push!(shell["Ratio"], sum(lengths) / 6 / law(centroid, d))
        end
    end
    rows = Dict{String, Any}[]
    for shell in shells
        at(values, q) = isempty(values) ? nothing :
            sort(values)[clamp(ceil(Int, q * length(values)), 1, length(values))]
        push!(rows, Dict{String, Any}(
            "Lower" => shell["Lower"], "Upper" => shell["Upper"], "Cells" => length(shell["Mean"]),
            "MeanEdgeP50" => at(shell["Mean"], 0.5), "MeanEdgeP90" => at(shell["Mean"], 0.9),
            "LongestEdgeP50" => at(shell["Longest"], 0.5), "LongestEdgeP90" => at(shell["Longest"], 0.9),
            "LongestEdgeMaximum" => at(shell["Longest"], 1.0),
            "AchievedOverPrescribedP50" => at(shell["Ratio"], 0.5),
            "AchievedOverPrescribedP90" => at(shell["Ratio"], 0.9)))
    end
    return rows
end

# Distance of a point to the infinite line through a segment.
function line_point_distance(point, start, stop)
    v = stop .- start
    return norm(cross(point .- start, v)) / norm(v)
end

# Achieved first-layer transverse size against a segment set (the junction lines):
# for every element of one dimension (the matching-surface elements of a physical
# label, or every tetrahedron) with a node on a segment, the largest distance of
# its nodes to the line through that segment - the thickness of the first cell
# layer against the line, which the band law prescribes as `prescribed`
# (NormalSize at r = 0). Percentiles and achieved-over-prescribed ratios are
# reported; nothing is gated.
function first_layer_transverse_statistics(points, index, segments, prescribed, dim, attribute)
    on_line = 1.0e-6 * prescribed
    thickness = Float64[]
    entities = dim == 3 ? [-1] : gmsh.model.getEntitiesForPhysicalGroup(dim, attribute)
    for entity in entities
        types, element_tags, element_nodes = gmsh.model.mesh.getElements(dim, entity)
        for (type, tags, block) in zip(types, element_tags, element_nodes)
            isempty(tags) && continue
            dim == 3 && Int(type) != GMSH_LINEAR_ELEMENT_TYPE[3] && continue
            width = length(block) ÷ length(tags)
            for start in 0:width:(length(block) - 1)
                nodes = [points[:, index[block[start + i]]] for i in 1:width]
                nearest = 0
                for (k, (a, b)) in enumerate(segments)
                    any(span_point_distance(node, a, b) <= on_line for node in nodes) || continue
                    nearest = k
                    break
                end
                nearest == 0 && continue
                a, b = segments[nearest]
                push!(thickness, maximum(line_point_distance(node, a, b) for node in nodes))
            end
        end
    end
    isempty(thickness) && return Dict{String, Any}("Elements" => 0, "Prescribed" => prescribed)
    sort!(thickness)
    at(q) = thickness[clamp(ceil(Int, q * length(thickness)), 1, length(thickness))]
    return Dict{String, Any}("Elements" => length(thickness), "Prescribed" => prescribed,
                             "TransverseP10" => at(0.1), "TransverseP50" => at(0.5),
                             "TransverseP90" => at(0.9), "TransverseMaximum" => thickness[end],
                             "AchievedOverPrescribedP50" => at(0.5) / prescribed,
                             "AchievedOverPrescribedP90" => at(0.9) / prescribed)
end

# The segments of the footprint polygons split into the sides on the outer box
# (box sides, not feature curves) and the interior sides (feature curves carrying
# the band law).
function footprint_polygon_segments(footprint_polygons, lower, upper, tolerance)
    interior = Tuple{Vector{Float64}, Vector{Float64}}[]
    on_box = Tuple{Vector{Float64}, Vector{Float64}}[]
    for polygon in footprint_polygons
        polygon_points = polygon["Points"]
        plane = Float64(polygon["Plane"])
        for i in eachindex(polygon_points)
            a = polygon_points[i]; b = polygon_points[mod1(i + 1, length(polygon_points))]
            segment = ([a[1], a[2], plane], [b[1], b[2], plane])
            bounds = [min(a[1], b[1]), min(a[2], b[2]), plane, max(a[1], b[1]), max(a[2], b[2]), plane]
            push!(on_outer_box(bounds, lower, upper, tolerance) ? on_box : interior, segment)
        end
    end
    return interior, on_box
end

# Quality of the tetrahedra with a vertex within `reach` of a point (the cap
# regions: the tube cap triangles meet these cells).
function point_region_tetrahedra(points, index, center, reach)
    scaled = Float64[]; condition = Float64[]
    types, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
    for (type, tags, block) in zip(types, element_tags, element_nodes)
        Int(type) == GMSH_LINEAR_ELEMENT_TYPE[3] || continue
        for start in 0:4:(length(block) - 1)
            nodes = [points[:, index[block[start + i]]] for i in 1:4]
            any(norm(node .- center) <= reach for node in nodes) || continue
            push!(scaled, tetrahedron_scaled_jacobian(nodes))
            push!(condition, tetrahedron_aspect(nodes))
        end
    end
    return Dict{String, Any}("Point" => collect(center), "Reach" => reach, "Cells" => length(scaled),
                             "MinimumScaledJacobian" => isempty(scaled) ? nothing : minimum(scaled),
                             "MaximumJacobianCondition" => isempty(condition) ? nothing : maximum(condition))
end

# The prism tube census (source-local frame; reported, the gates are applied by
# gate_mixed_volume_quality and the element budget).
function prism_tube_census(tubes, segments, description, states, volume_census, face_meshes,
                           curve_meshes, timings, band_record, quality, graded_points,
                           semantic_corners, junction_curves, band_curves, footprint_polygons,
                           box, lc_fine, lc_tangent, lc_far, far_growth, trace_basis_record,
                           max_elements, element_count, layer_records)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    radius = description["Radius"]
    lower, upper, box_tolerance = box
    on_box(point) = on_outer_box(vcat(point, point), lower, upper, box_tolerance)
    tube_rows = Dict{String, Any}[]
    tubes_per_side = description["TubesPerSide"]
    for (k, (tube, section)) in enumerate(tubes)
        segment = segments[(k + tubes_per_side - 1) ÷ tubes_per_side]
        start_point = tube_point(tube, 0.0, 0.0, tube.s_start)
        end_point = tube_point(tube, 0.0, 0.0, tube.s_end)
        # An arc tube's Origin / Normal / Extrusion are its start frame (design 1.2 (5)).
        frame = tube isa ArcTube ? arc_frame(tube, tube.s_start) : (origin=tube.origin, n=tube.n, e=tube.e)
        push!(tube_rows, Dict{String, Any}(
            "Origin" => frame.origin, "Normal" => frame.n, "Extrusion" => frame.e,
            "Start" => tube.s_start, "End" => tube.s_end, "Length" => tube.s_end - tube.s_start,
            "StartPoint" => start_point, "EndPoint" => end_point,
            "EndsOnBox" => [on_box(start_point), on_box(end_point)],
            "Layers" => tube.layers, "Spacing" => tube_spacing(tube),
            "LayerThickness" => layer_records[k],
            "Materials" => section.materials, "Conductor" => segment.conductor,
            "Plane" => segment.plane, "Hole" => segment.hole, "Layer" => segment.layer_sign,
            "CornerAngles" => collect(segment.corner_angles),
            "Edge" => tubes_per_side == 1 ? "sheet" : isodd(k) ? "top" : "bottom"))
        # Face ends (design A2) and theta-0 legacy box-vertex corners (decision 320) are
        # recorded only where they occur, so a rectilinear census is unchanged.
        isempty(tube.face_ends) ||
            (tube_rows[end]["FaceEnds"] = [face_end_record(face_end) for face_end in tube.face_ends])
        if any(segment.legacy_box_corners)
            tube_rows[end]["LegacyBoxVertexCorner"] = [
                point for (point, legacy) in ((segment.start, segment.legacy_box_corners[1]),
                                              (segment.stop, segment.legacy_box_corners[2])) if legacy]
        end
        if tube isa ArcTube
            arc = segment.arc
            tube_rows[end]["Arc"] = Dict{String, Any}(
                "ArcId" => arc.id, "Centre" => [arc.centre[1], arc.centre[2]], "Radius" => arc.rho,
                "Sign" => Int(arc.sigma), "Part" => arc.part, "Parts" => arc.parts,
                "SweepDegrees" => rad2deg(abs(arc.sweep)), "RunSweepDegrees" => rad2deg(abs(arc.run_sweep)),
                "ThetaStartDegrees" => rad2deg(tube.theta0),
                "Orientation" => tube.orientation)
        end
        isempty(tube.joints) ||
            (tube_rows[end]["Joints"] = [joint_end_record(joint) for joint in tube.joints])
    end
    arc_rows = count(haskey(row, "Arc") for row in tube_rows)
    # The straight tubes' joint ends (Tubes[].Joints) and the arc part splits (an arc row whose
    # Part < Parts shares its end section with the next part) make up the shared sections.
    joint_rows = sum(length(row["Joints"]) for row in tube_rows if haskey(row, "Joints"); init=0)
    part_splits = count(haskey(row, "Arc") && row["Arc"]["Part"] < row["Arc"]["Parts"] for row in tube_rows)
    face_end_rows = [row["FaceEnds"] for row in tube_rows if haskey(row, "FaceEnds")]
    face_end_count = sum(length, face_end_rows; init=0)
    legacy_box_corners = sum(length(row["LegacyBoxVertexCorner"]) for row in tube_rows
                             if haskey(row, "LegacyBoxVertexCorner"); init=0)
    thicknesses = reduce(vcat, tube_layer_thicknesses(tube) for (tube, _) in tubes)
    spacings = [row["Spacing"] for row in tube_rows]
    neighbour_ratios = [row["LayerThickness"]["MaximumNeighbourRatio"] for row in tube_rows]
    achieved_over_prescribed = [row["LayerThickness"]["AchievedOverPrescribed"]["P50"]
                                for row in tube_rows]
    inner_arc = description["RingSizes"][1] * deg2rad(description["SectorDegrees"])
    caps = [point_region_tetrahedra(points, index, collect(point), radius)
            for point in graded_points[(length(semantic_corners) + 1):end]]
    cap_scaled = filter(!isnothing, [cap["MinimumScaledJacobian"] for cap in caps])
    cap_condition = filter(!isnothing, [cap["MaximumJacobianCondition"] for cap in caps])
    junction_segments = curve_segments(junction_curves)
    footprint_segments, footprint_box_segments =
        footprint_polygon_segments(footprint_polygons, box...)
    band_curve_segments = curve_segments(band_curves)
    # Achieved volume sizes against the three volume laws of decision 39, in
    # NormalSize shells: around the trace-basis apexes (within CornerIsotropyRadius),
    # outside the corner balls (semantic corners) and around the junction lines.
    corner_radius = band_record["CornerIsotropyRadius"]
    band_record["Achieved"] = Dict{String, Any}(
        "Rule" => "tetrahedra binned by centroid distance in NormalSize shells; the mean edge " *
                  "over the law evaluated at the centroid (AchievedOverPrescribed); " *
                  "LongestEdgeP50 is the physics-09 observable (median longest edge within " *
                  "CornerIsotropyRadius of the narrow-hat apexes)",
        "TraceApexes" => isempty(TRACE_BASIS_APEXES) ? nothing : Dict{String, Any}(
            "Apexes" => length(TRACE_BASIS_APEXES), "Reach" => corner_radius,
            "Law" => "TraceBasisSizing (volume rule)",
            "Shells" => shell_size_statistics(
                points, index,
                c -> minimum(norm(c .- collect(apex)) for apex in TRACE_BASIS_APEXES),
                0.0, [lc_fine, 2lc_fine, corner_radius],
                (c, d) -> trace_basis_size(c[1], c[2], c[3], lc_far))),
        "CornerExterior" => Dict{String, Any}(
            "Corners" => length(semantic_corners), "Law" => "CornerExteriorRule",
            "Shells" => shell_size_statistics(
                points, index,
                c -> minimum(norm(c .- collect(corner)) for corner in semantic_corners),
                corner_radius, corner_radius .+ [2lc_fine, 4lc_fine, 8lc_fine],
                (c, d) -> corner_exterior_size(d, lc_fine, lc_far, far_growth))),
        "JunctionLines" => Dict{String, Any}(
            "Segments" => length(junction_segments), "Law" => "BandRule",
            "Shells" => shell_size_statistics(
                points, index,
                c -> minimum(span_point_distance(c, a, b) for (a, b) in junction_segments),
                0.0, [lc_fine, 2lc_fine, 4lc_fine, 8lc_fine],
                (c, d) -> band_law_size(d, lc_fine, lc_far, far_growth))))
    return Dict{String, Any}(
        "Rule" => "every straight metal edge (top and bottom edge of every metal segment) " *
                  "carries a prism tube on its dielectric side, except the untubed short sides " *
                  "(Scope.UntubedEdges: a tube interval below the inner ring size; their corner " *
                  "balls mesh them at CornerSize): geometric rings of size " *
                  "EdgeSize x GrowthRatio^(k-1) (RingSizes) over the sectors (SectorDegrees), " *
                  "extruded along the edge in Layers whose thickness follows the composed size " *
                  "field on the axis (LayerRule, TubeAxisSizeLaw; Spacing = the largest layer " *
                  "<= TangentialSize), the lateral " *
                  "quadrangles closed by explicit pyramids (PyramidHeight); the tube ends at the " *
                  "outer box and CornerClearance before a semantic corner (CornerClearanceRule); " *
                  "the cap centre before a corner is a graded point of the corner law " *
                  "(CornerSize = EdgeSize, the same shells as the rings), so the un-tubed edge " *
                  "part and the cap region are tetrahedra graded from EdgeSize; no tetrahedral " *
                  "edge layer",
        "Section" => description,
        "InnerSize" => description["RingSizes"][1], "GrowthRatio" => description["RingSizes"][2] / description["RingSizes"][1],
        "TangentialSize" => lc_tangent, "NormalSize" => lc_fine, "FarSize" => lc_far,
        "FarGrowth" => far_growth,
        "Tubes" => tube_rows, "TubeCount" => length(tubes),
        "TotalTubeLength" => sum(row["Length"] for row in tube_rows),
        "Layers" => sum(row["Layers"] for row in tube_rows),
        "FaceEnds" => Dict{String, Any}(
            "Rule" => description["FaceEndRule"], "BoxVertexRule" => description["BoxVertexRule"],
            "Count" => face_end_count,
            "EndBlockLayers" => sum(record["Layers"] for records in face_end_rows for record in records;
                                    init=0),
            "LegacyBoxVertexCorners" => legacy_box_corners),
        "SpacingMinimum" => minimum(thicknesses), "SpacingMaximum" => maximum(thicknesses),
        # The largest face-end block spacing lc_end (0 without face ends): the bound the
        # tube design statement and the census validator add to TangentialSize (design A2).
        "FaceEndSpacingMaximum" => maximum([record["EndSpacing"] for records in face_end_rows
                                            for record in records]; init=0.0),
        # Arc tubes and smooth joints (design 1.2 (3) / A3 (2)), recorded only where they occur.
        "ArcTubes" => arc_rows == 0 ? nothing : Dict{String, Any}(
            "Rule" => description["ArcTubeRule"], "SmoothJointRule" => description["SmoothJointRule"],
            "Count" => arc_rows, "JointEnds" => joint_rows, "PartSplits" => part_splits,
            "SharedSections" => description["SharedSections"],
            "TotalArcLength" => sum(row["Length"] for row in tube_rows if haskey(row, "Arc"))),
        "LayerRule" => TUBE_LAYER_RULE, "TubeAxisSizeLaw" => TUBE_AXIS_SIZE_LAW,
        "LayerGrowthCap" => description["RingSizes"][2] / description["RingSizes"][1],
        "LayerThickness" => Dict{String, Any}(
            "Rule" => "over every layer of every tube: Minimum / P50 / Maximum thickness, the " *
                      "largest neighbour ratio, the per-tube P50 achieved-over-prescribed " *
                      "(layer thickness over the gradient-limited axis size at its midpoint) and " *
                      "the number of layers below TangentialSize / GrowthRatio (the layers the " *
                      "axis field refines beyond the tangential grid); " *
                      "per tube in Tubes[].LayerThickness (ends: AtStart / AtEnd against " *
                      "PrescribedAtStart / PrescribedAtEnd; at an end on the outer box, EndsOnBox, " *
                      "the end layer is the surface value: AtStart / AtEnd <= the prescribed size)",
            "Minimum" => minimum(thicknesses),
            "P50" => sort(thicknesses)[cld(length(thicknesses), 2)],
            "Maximum" => maximum(thicknesses),
            "MaximumNeighbourRatio" => maximum(neighbour_ratios),
            "AchievedOverPrescribedP50Range" => [minimum(achieved_over_prescribed),
                                                 maximum(achieved_over_prescribed)],
            "LayersBelowTangentialSizeOverGrowthRatio" =>
                count(<(lc_tangent / (description["RingSizes"][2] / description["RingSizes"][1])),
                      thicknesses)),
        "InnermostArc" => inner_arc,
        "MaximumPrismEdgeAspect" => maximum(spacings) / min(inner_arc, description["RingSizes"][1]),
        "Volumes" => volume_census,
        "Prisms" => sum(Int[row["Prisms"] for row in volume_census]),
        "Pyramids" => sum(Int[row["Pyramids"] for row in volume_census]),
        "Faces" => face_meshes, "Curves" => curve_meshes, "Timings" => timings,
        "SizeLaws" => band_record,
        "Quality" => quality,
        "CapRegions" => Dict{String, Any}(
            "Rule" => "tetrahedra with a vertex within the tube radius of a cap centre before a corner",
            "Caps" => length(caps),
            "MinimumScaledJacobian" => isempty(cap_scaled) ? nothing : minimum(cap_scaled),
            "MaximumJacobianCondition" => isempty(cap_condition) ? nothing : maximum(cap_condition),
            "Regions" => caps),
        "CutSurface" => Dict{String, Any}(
            "Rule" => "edge lengths of the matching-surface elements (label 1): all, within " *
                      "NormalSize of a junction line, and against the trace rule",
            "All" => surface_size_statistics(points, index, 1, nothing, 0.0),
            "NearJunctions" => surface_size_statistics(points, index, 1, junction_segments, lc_fine),
            "TraceRule" => trace_basis_record === nothing ? nothing :
                           trace_basis_achieved_sizes(points, index, 1, trace_basis_record)),
        "Bands" => Dict{String, Any}(
            "Rule" => "mean edge length of the tetrahedra whose centroid lies within NormalSize " *
                      "of the junction lines / the footprint edges (Footprint: the footprint " *
                      "sides not on the outer box, the feature curves; FootprintOnBox: the box " *
                      "sides, a far-field sample); FirstLayer: the achieved first-layer transverse " *
                      "size against the junction lines (matching-surface elements and tetrahedra " *
                      "with a node on a line) over the prescribed NormalSize (BandRule at r = 0)",
            "Prescribed" => lc_fine,
            "Junction" => band_tetrahedron_statistics(points, index, junction_segments, lc_fine),
            "JunctionFirstLayer" => Dict{String, Any}(
                "CutSurface" => first_layer_transverse_statistics(points, index, junction_segments,
                                                                  lc_fine, 2, 1),
                "Tetrahedra" => first_layer_transverse_statistics(points, index, junction_segments,
                                                                  lc_fine, 3, 0)),
            "Footprint" => band_tetrahedron_statistics(points, index, footprint_segments, lc_fine),
            "FootprintOnBox" => band_tetrahedron_statistics(points, index, footprint_box_segments,
                                                            lc_fine)),
        "BandCurves" => Dict{String, Any}(
            "Rule" => "the non-metal longitudinal feature curves (junction lines and footprint " *
                      "edges parallel to a metal edge) are 1D-meshed at Spacing = NormalSize, the " *
                      "band law on the line, composed with the corner law and, since decision 43, " *
                      "the whole size field along the curve (census CurveSpacing); the metal " *
                      "ridges keep the TangentialSize grid (they carry the tubes)",
            "Count" => length(band_curves),
            "TotalLength" => sum(norm(b .- a) for (a, b) in band_curve_segments; init=0.0),
            "Spacing" => lc_fine,
            "Segments" => [vcat(a, b) for (a, b) in band_curve_segments]),
        "FarFieldBudgetPolicy" => Dict{String, Any}(
            "Name" => "gmsh-only-fail-closed-cap",
            "Rule" => "the requested far size is used as is (Pressure 1); the element budget " *
                      "is a fail-closed gate on the generated mesh",
            "RequestedFarSize" => lc_far, "EffectiveFarSize" => lc_far, "Pressure" => 1.0,
            "Elements" => element_count, "MaximumElements" => max_elements,
            "ElementBudgetFraction" => element_count / max_elements))
end

# Straight segments (start, stop) of a curve list from the CAD endpoints.
# With `subdivide`, a CIRCLE curve (an arc side's ridge, its concentric trench / collar
# arcs) is sampled into chords whose sagitta stays below `sagitta`
# (CURVE_SEGMENT_SAGITTA_OVER_RADIUS x the arc radius by default), so the exact
# segment-distance size laws see the arc, not a chord; every other curve, and every curve
# without `subdivide` (the census records: one segment per curve), is its chord as before.
const CURVE_SEGMENT_SAGITTA_OVER_RADIUS = 1.0e-3

function curve_segments(curves; subdivide=false, sagitta=nothing)
    segments = Tuple{Vector{Float64}, Vector{Float64}}[]
    for curve in curves
        lower, upper = gmsh.model.getParametrizationBounds(1, curve)
        if !subdivide || gmsh.model.getType(1, curve) != "Circle"
            xyz = gmsh.model.getValue(1, curve, [lower[1], upper[1]])
            push!(segments, (collect(xyz[1:3]), collect(xyz[4:6])))
            continue
        end
        curve_length = gmsh.model.occ.getMass(1, curve)
        # OCC parametrises circles by the angle: the radius follows from the length.
        angle = max(upper[1] - lower[1], 1.0e-12)
        radius = curve_length / angle
        bound = sagitta === nothing ? CURVE_SEGMENT_SAGITTA_OVER_RADIUS * radius : sagitta
        step = 2.0 * acos(clamp(1.0 - bound / radius, -1.0, 1.0))
        count = max(1, ceil(Int, angle / max(step, 1.0e-9)))
        count = min(count, 4096)
        parameters = collect(range(lower[1], upper[1]; length=count + 1))
        xyz = reshape(gmsh.model.getValue(1, curve, parameters), 3, :)
        for i in 1:count
            push!(segments, (collect(xyz[:, i]), collect(xyz[:, i + 1])))
        end
    end
    return segments
end

function generate_spatial_coupon(;
    signature::String,
    mask::Union{Nothing, String}=nothing,
    boundary::Union{Nothing, String}=nothing,
    etch_boundary::Union{Nothing, String}=nothing,
    fabricated::Bool,
    radius::Float64          = 2.0,
    metal_thickness::Float64 = 0.1,
    overetch::Float64        = 0.05,
    sidewall_angle::Float64  = 80.0,
    top_rounding::Float64    = 0.01,
    trench_rounding::Float64 = 0.01,
    lc_fine::Float64         = 0.02,
    lc_tangent::Float64      = 0.0,
    lc_far::Float64          = 0.3,
    process_core_width::Float64 = 0.0,
    process_fine_width::Float64 = 0.0,
    process_grading_power::Float64 = 1.7,
    max_nodes::Int           = 500_000,
    max_elements::Int        = 2_000_000,
    mesh_order::Int          = 1,
    mesh_control::Union{Nothing, Function}=nothing,
    surface_constraints::Union{Nothing, Function}=nothing,
    mesh_postprocess::Union{Nothing, Function}=nothing,
    optimize_volume::Bool=true,
    matching_trace::Union{Nothing,String}=nothing,
    matching_trace_mode::String="all",
    geometry_only::Bool=false,
    transform::Matrix{Float64}=copy(IDENTITY_RIGID_TRANSFORM),
    semantic_contract::Union{Nothing, String}=nothing,
    corner_isotropy_radius::Float64=0.0,
    corner_census::Union{Nothing, String}=nothing,
    labels_only::Union{Nothing, String}=nothing,
    trace_basis_contract::Union{Nothing, String}=nothing,
    trace_vertices::Union{Nothing, String}=nothing,
    trace_triangles::Union{Nothing, String}=nothing,
    process_library::Union{Nothing, String}=nothing,
    trace_basis_size_ratio::Float64=1.0,
    edge_size::Float64=0.0,
    edge_growth_ratio::Float64=2.0,
    edge_layer_aspect::Float64=4.0,
    maximum_corner_aspect::Float64=0.0,
    minimum_scaled_jacobian::Float64=0.0,
    maximum_jacobian_condition::Float64=0.0,
    quality_displacement_over_normal::Float64=0.0,
    edge_layer_maximum_aspect::Float64=0.0,
    corner_size::Float64=0.0,
    prism_tubes::Bool=false,
    tube_sector_degrees::Float64=30.0,
    far_growth::Float64=0.0,
    corner_shape_gate::Float64=0.0,
    filename::String
)
    matching_trace_mode in ("all","sides","levels","none") ||
        error("Unknown matching trace constraint mode")
    radius > 0.0 || error("radius must be positive")
    metal_thickness > 0.0 || error("metal thickness must be positive")
    0.0 <= overetch < radius || error("overetch must lie in [0, radius)")
    0.0 < sidewall_angle <= 90.0 || error("sidewall angle must lie in (0, 90]")
    0.0 <= top_rounding < metal_thickness ||
        error("top rounding must be smaller than metal thickness")
    0.0 <= trench_rounding <= overetch || error("trench rounding must not exceed overetch")
    lc_fine > 0.0 || error("fine mesh size must be positive")
    lc_far >= lc_fine || error("far mesh size must not be smaller than fine mesh size")
    # Coupon-scale size bound: FarSize (FarSizeOverRadius x Radius) is the coarsest
    # size the coupon admits - the size at its matching surface - so the tangential
    # spacing (the along-edge tube extrusion / ridge grid, a coarsening size fixed by
    # the recipe in absolute units) is bounded by it: TangentialSize = min(--lc-tangent,
    # FarSize).  Dimensionless: the bound acts exactly when Radius < --lc-tangent /
    # FarSizeOverRadius and leaves every larger coupon unchanged.  The resolution sizes
    # (NormalSize, EdgeSize = CornerSize) are not bounded: a fine size above FarSize is
    # a contradictory recipe and fails closed above.  Requested and bound values are
    # recorded in the census (SizeBounds) and bound by the stage contract.
    lc_tangent >= 0.0 || error("tangential mesh size must be nonnegative")
    requested_lc_tangent = lc_tangent
    lc_tangent = lc_tangent == 0.0 ? 0.0 : min(lc_tangent, lc_far)
    (lc_tangent == 0.0 || lc_fine <= lc_tangent) ||
        error("tangential mesh size must not be smaller than the fine size")
    process_core_width = process_core_width > 0.0 ?
                         process_core_width :
                         max(2metal_thickness, 4overetch, 8lc_fine)
    process_fine_width >= 0.0 || error("process fine width must be nonnegative")
    process_core_width > process_fine_width ||
        error("process core width must exceed the fully fine region width")
    process_grading_power > 0.0 || error("process grading power must be positive")
    max_nodes > 0 || error("maximum node budget must be positive")
    max_elements > 0 || error("maximum element budget must be positive")
    transform = rigid_transform(vec(transform'))
    corner_isotropy_radius >= 0.0 || error("corner isotropy radius must be nonnegative")
    corner_isotropy = semantic_contract !== nothing
    (corner_isotropy == (corner_isotropy_radius > 0.0) == (corner_census !== nothing)) ||
        error("Corner isotropy requires --semantic-contract, --corner-isotropy-radius " *
              "and --corner-census together")
    corner_isotropy && corner_census == filename &&
        error("Corner census must not overwrite the mesh output")
    # The labels-only pass (decision 62(2): the registration probe) stops right after the
    # CAD-entity labelling, before any mesh generation; it needs the semantic contract
    # (the census fields it records) and writes neither the mesh nor the corner census
    # (its output may take the corner census path: it is the probe's census).
    labels_only === nothing || corner_isotropy ||
        error("--labels-only requires the semantic-corner isotropy options")
    labels_only !== nothing && labels_only == filename &&
        error("Labels-only output must not overwrite the mesh output")
    # The seed-side required-region gates (decision 30) come together and need the
    # corner balls; without them the seed is neither optimized nor gated.
    seed_quality_gates = maximum_corner_aspect > 0.0
    (seed_quality_gates == (minimum_scaled_jacobian > 0.0) ==
     (maximum_jacobian_condition > 0.0) == (quality_displacement_over_normal > 0.0)) ||
        error("--maximum-corner-aspect, --minimum-scaled-jacobian, " *
              "--maximum-jacobian-condition and --maximum-quality-displacement-over-normal " *
              "are required together")
    seed_quality_gates && !corner_isotropy &&
        error("Seed quality gates require the semantic-corner isotropy options")
    (!seed_quality_gates || (maximum_corner_aspect > 1.0 && minimum_scaled_jacobian < 1.0 &&
                             isfinite(maximum_jacobian_condition) &&
                             maximum_jacobian_condition > 1.0)) ||
        error("Seed quality gates must satisfy MaximumCornerAspect > 1, " *
              "MinimumScaledJacobian < 1 and a finite MaximumJacobianCondition > 1")
    # The edge-layer quality rule (decision 32, calibration only) needs the seed
    # gates and an edge layer to apply to.
    isfinite(edge_layer_maximum_aspect) && edge_layer_maximum_aspect >= 0.0 ||
        error("edge layer maximum aspect must be a nonnegative finite number")
    edge_layer_maximum_aspect == 0.0 || (edge_layer_maximum_aspect > 1.0 && seed_quality_gates &&
                                         edge_size > 0.0) ||
        error("--edge-layer-maximum-aspect requires a bound above 1, the seed quality gates " *
              "and --edge-size")
    # The seed honors the same isotropic corner ball the metric stage prescribes:
    # NormalSize (lc_fine) inside CornerIsotropyRadius around every contract
    # semantic corner, graded to the far size with the process-band slope.
    semantic_corners = corner_isotropy ? read_semantic_corners(semantic_contract, transform) :
                       NTuple{3, Float64}[]
    trace_basis_paths = (trace_basis_contract, trace_vertices, trace_triangles, process_library)
    trace_basis_bound = all(path -> path !== nothing, trace_basis_paths)
    trace_basis_bound || all(path -> path === nothing, trace_basis_paths) ||
        error("Trace basis sizing requires --trace-basis-contract, --trace-vertices, " *
              "--trace-triangles and --process-library together")
    isfinite(trace_basis_size_ratio) && trace_basis_size_ratio > 0.0 ||
        error("Trace basis size ratio must be a positive finite number")
    # The trace rule composes with the scalar corner-isotropy background field
    # through the size callback; it is recorded in the census, so it needs both.
    !trace_basis_bound || corner_isotropy ||
        error("Trace basis sizing requires the semantic contract corner isotropy and census")
    trace_basis = trace_basis_bound ? read_trace_basis(trace_basis_contract, trace_vertices,
                                                       trace_triangles, process_library) : nothing
    # Geometric transverse edge layer seeded on the metal-edge faces (EdgeSize at the
    # edge, ratio per layer, up to the normal size lc_fine); it needs the ridge grid
    # (lc_tangent) and is recorded in the census for the metric stage.
    edge_size >= 0.0 || error("edge size must be nonnegative")
    edge_layer = edge_size > 0.0
    !edge_layer || edge_size < lc_fine ||
        error("edge size must be smaller than the fine (normal) mesh size")
    isfinite(edge_growth_ratio) && edge_growth_ratio > 1.0 ||
        error("edge growth ratio must exceed 1")
    isfinite(edge_layer_aspect) && edge_layer_aspect >= 1.0 ||
        error("edge layer aspect must be at least 1")
    !edge_layer || (lc_tangent > 0.0 && corner_isotropy) ||
        error("The edge layer requires --lc-tangent and the semantic contract corner census")
    edge_layer_offsets = edge_layer ?
        edge_layer_row_offsets(edge_size, edge_growth_ratio, lc_fine) : Float64[]
    # Corner grading (decision 33): CornerSize at every semantic corner point growing
    # by the edge-layer ratio to lc_fine inside the corner ball (0 = uniform lc_fine,
    # the production ball); the grading must reach lc_fine inside the ball.
    isfinite(corner_size) && corner_size >= 0.0 || error("corner size must be nonnegative")
    corner_size == 0.0 || corner_isotropy ||
        error("--corner-size requires the semantic contract corner isotropy and census")
    corner_size == 0.0 || corner_size < lc_fine ||
        error("corner size must be smaller than the fine (normal) mesh size")
    corner_grading = CornerGrading(corner_size, edge_growth_ratio, lc_fine, corner_isotropy_radius)
    corner_size == 0.0 || corner_grading_reach(corner_grading) <= corner_isotropy_radius ||
        error("The corner grading must reach the normal size inside the corner ball")
    # Prism edge tubes (decision 38): EdgeSize is the inner ring size and the edge
    # growth ratio the ring ratio; the tetrahedral edge layer is not built. The
    # tubes need the tangential spacing, the corner balls with their grading, the
    # seed quality gates (applied per element type) and the far growth of the
    # band grading around the tubes.
    if prism_tubes
        edge_size > 0.0 || error("Prism tubes require --edge-size (the inner ring size)")
        lc_tangent > 0.0 || error("Prism tubes require --lc-tangent (the extrusion spacing)")
        corner_isotropy && corner_size > 0.0 ||
            error("Prism tubes require the semantic contract corner balls with --corner-size")
        corner_size == edge_size ||
            error("Prism tubes require CornerSize == EdgeSize: one graded law for the semantic " *
                  "corners and the tube cap centres")
        seed_quality_gates || error("Prism tubes require the seed quality gates")
        isfinite(far_growth) && far_growth > 0.0 ||
            error("Prism tubes require --far-growth > 0 (band grading around the tubes)")
        isfinite(tube_sector_degrees) && 0.0 < tube_sector_degrees <= 90.0 ||
            error("Tube sector angle must lie in (0, 90] degrees")
        sidewall_angle == 90.0 || scope_error("SlopedSidewalls", "SidewallAngle $sidewall_angle")
        top_rounding == 0.0 || scope_error("TopRounding", "TopRounding $top_rounding")
        trench_rounding == 0.0 || scope_error("TrenchRounding", "TrenchRounding $trench_rounding")
        matching_trace === nothing && surface_constraints === nothing ||
            error("Prism tubes cannot be combined with matching-trace or surface constraints")
        mesh_order == 1 || error("Prism tubes require a linear mesh")
        edge_layer_maximum_aspect == 0.0 || error("Prism tubes replace the edge layer rule")
        edge_layer = false
        edge_layer_offsets = Float64[]
    else
        far_growth == 0.0 || error("--far-growth belongs to the prism tube recipe")
    end

    edges = read_edges(signature)
    if length(unique(edge.slot for edge in edges)) > 1 &&
       mesh_postprocess === nothing && !geometry_only
        error("Multi-slot tetrahedral coupons require an explicit surface-partition postprocessor; whole-CAD-face labeling is not sufficient")
    end
    facets = read_mask(mask)
    boundary_loops = read_boundary(boundary)
    isempty(boundary_loops) ||
        !isempty(facets) ||
        error("A classified plan-view boundary requires the corresponding mask facets")
    lower, upper = coupon_bounds(edges, radius, metal_thickness, overetch;
                                 support_box=read_support_box(process_library))
    tolerance = 1.0e-7 * radius
    # The tagged arcs of the boundary (design A1 (4)): the mask tests see the exact arcs.
    arc_boundary = any(loop_has_arcs(loop) for loop in boundary_loops)
    register_mask_arc_segments!(boundary_loops, tolerance)
    if trace_basis !== nothing
        all(abs(trace_basis.lower[d] - lower[d]) <= tolerance &&
            abs(trace_basis.upper[d] - upper[d]) <= tolerance for d in 1:3) ||
            error("Trace basis box differs from the coupon box")
    end
    outer_tolerance = 1.0e-4 * radius
    validate_plan_view_geometry(edges, radius, tolerance, facets)
    # Design round 2 F5-A (decision 363): the kind of every semantic corner by the exact
    # dot-product predicate on its two plan-view boundary sides (INVARIANT_CORNER_RULE;
    # integer arithmetic on the generator's 1e-9 R quantum counts, decision 416), checked
    # against the contract's Derivation.InvariantCorners (fail closed on a disagreement);
    # invariant corners need the CornerShapeGate. Without a classified boundary (no
    # plan-view loops) every corner keeps the legacy convention.
    corner_shape_gate == 0.0 || (isfinite(corner_shape_gate) && corner_shape_gate > 1.0) ||
        error("--corner-shape-gate must be a finite number above 1 (0: none)")
    corner_kinds, corner_sides = if corner_isotropy && !isempty(boundary_loops)
        kinds, sides, _ = semantic_corner_kinds(semantic_corners, boundary_loops, tolerance,
                                                PLAN_VIEW_QUANTUM_OVER_RADIUS * radius)
        check_invariant_corner_contract(semantic_corners, kinds,
                                        read_invariant_corners(semantic_contract, transform),
                                        tolerance)
        kinds, sides
    else
        fill(:legacy, length(semantic_corners)), fill(nothing, length(semantic_corners))
    end
    invariant_corner_count = count(kind -> kind === :invariant, corner_kinds)
    invariant_corner_count == 0 || !seed_quality_gates || corner_shape_gate > 1.0 ||
        error("$(invariant_corner_count) invariant (non-perpendicular) semantic corners " *
              "$([collect(c) for (c, k) in zip(semantic_corners, corner_kinds) if k === :invariant]) " *
              "need --corner-shape-gate (CornerShapeGate = min(E_pop, 5.0); $(INVARIANT_CORNER_RULE))")
    println("Semantic corner kinds: $(length(semantic_corners) - invariant_corner_count) legacy " *
            "(exactly perpendicular), $(invariant_corner_count) invariant")
    layers = layer_groups(edges, tolerance)
    pullback_metal = metal_thickness / tan(deg2rad(sidewall_angle))
    etch_loops = etch_boundary === nothing ? nothing : read_boundary(etch_boundary)
    if etch_loops !== nothing
        # The bound footprint is the trench's plan view: a fabricated build lofts it and
        # needs sharp vertical geometry; a thin build (decision 66) etches nothing - it
        # checks the footprint carries every metal edge and records the bound file.
        !fabricated || (sidewall_angle == 90.0 && top_rounding == 0.0 && trench_rounding == 0.0) ||
            error("Explicit etch footprints require sharp vertical fabricated geometry")
        isempty(etch_loops) && error("Empty explicit etch footprint")
    end
    pullback_trench = overetch > 0.0 ? overetch / tan(deg2rad(sidewall_angle)) : 0.0
    scope_classes = exhibited_scope_classes(edges, boundary_loops, layers, fabricated,
                                            sidewall_angle, top_rounding, trench_rounding,
                                            overetch, etch_loops !== nothing,
                                            trace_basis !== nothing)

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("spatial_coupon_$(fabricated ? "fabricated" : "thin")")
    occ = gmsh.model.occ
    outer = (
        3,
        occ.addBox(
            lower[1],
            lower[2],
            lower[3],
            upper[1] - lower[1],
            upper[2] - lower[2],
            upper[3] - lower[3]
        )
    )

    substrates = Tuple{Int32, Int32}[]
    layer_substrates = Vector{Vector{Tuple{Int32, Int32}}}()
    # Every etch footprint polygon (device or producer default) as lofted, after
    # collinear-edge simplification; recorded in the census when it is written.
    footprint_polygons = corner_isotropy ? Dict{String, Any}[] : nothing
    for layer in layers
        slab = if layer.sign > 0
            [(
                3,
                occ.addBox(
                    lower[1],
                    lower[2],
                    lower[3],
                    upper[1] - lower[1],
                    upper[2] - lower[2],
                    layer.plane - lower[3]
                )
            )]
        else
            [(
                3,
                occ.addBox(
                    lower[1],
                    lower[2],
                    layer.plane,
                    upper[1] - lower[1],
                    upper[2] - lower[2],
                    upper[3] - layer.plane
                )
            )]
        end
        if fabricated && overetch > 0.0
            layer_loops = [
                loop for
                loop in boundary_loops if abs(loop.plane - layer.plane) <= tolerance
            ]
            trenches = if etch_loops !== nothing
                footprint = [loop for loop in etch_loops
                             if abs(loop.plane - layer.plane) <= tolerance]
                isempty(footprint) && error("Explicit etch footprint is missing a process layer")
                loft_mask(occ, footprint, layer.plane,
                          layer.plane - layer.sign * overetch, 0.0, tolerance;
                          simplify=true, footprint=footprint_polygons)
            elseif isempty(boundary_loops)
                result = Tuple{Int32, Int32}[]
                for edge in layer.edges
                    append!(
                        result,
                        loft_strip(
                            occ,
                            edge,
                            radius,
                            1.0,
                            layer.plane,
                            layer.plane - layer.sign * overetch,
                            pullback_trench;
                            footprint=footprint_polygons
                        )
                    )
                end
                fuse_all(occ, result)
            else
                isempty(layer_loops) && error(
                    "Plan-view boundary is missing fabrication layer $(layer.plane)"
                )
                boundary_strips(
                    occ,
                    layer_loops,
                    radius,
                    layer.plane,
                    layer.plane - layer.sign * overetch,
                    pullback_trench,
                    tolerance;
                    footprint=footprint_polygons,
                    box=(lower, upper)
                )
            end
            trench = if isempty(boundary_loops)
                fillet_plane_edges(
                    occ,
                    trenches,
                    trench_rounding,
                    layer.plane - layer.sign * overetch,
                    tolerance
                )
            else
                fillet_physical_edges(
                    occ,
                    trenches,
                    trench_rounding,
                    layer.plane - layer.sign * overetch,
                    physical_segments(layer_loops, -pullback_trench, tolerance),
                    tolerance
                )
            end
            slab, _ = occ.cut(slab, trench)
        end
        push!(layer_substrates, slab)
        append!(substrates, slab)
    end

    substrate_seed = Tuple{Int32, Int32}[]
    vacuum_seed = Tuple{Int32, Int32}[]
    domains = Tuple{Int32, Int32}[]
    tube_tools = Tuple{Int32, Int32}[]
    tube_records = TubeRecord[]
    tubes = Tuple{AbstractTube, TubeSection}[]
    tube_segments = NamedTuple[]
    tube_untubed_edges = Dict{String, Any}[]
    tube_joints = Tuple{Int, Int, Int, Int}[]
    tube_sections = Dict{String, Any}()
    tube_volumes = Int32[]
    tip_bisectors = Dict{String, Any}[]
    unrefined_tips = Dict{String, Any}[]
    if fabricated
        metal = Tuple{Int32, Int32}[]
        for layer in layers
            for conductor in sort!(unique(edge.conductor for edge in layer.edges))
                conductor_edges =
                    [edge for edge in layer.edges if edge.conductor == conductor]
                conductor_facets = [
                    facet for facet in facets if facet.conductor == conductor &&
                    abs(facet.plane - layer.plane) <= tolerance
                ]
                !isempty(facets) &&
                    isempty(conductor_facets) &&
                    error("Plan-view mask is missing conductor $conductor")
                conductor_loops = [
                    loop for loop in boundary_loops if loop.conductor == conductor &&
                    abs(loop.plane - layer.plane) <= tolerance
                ]
                conductor_metal = if isempty(boundary_loops)
                    result = Tuple{Int32, Int32}[]
                    for edge in conductor_edges
                        append!(
                            result,
                            loft_strip(
                                occ,
                                edge,
                                radius,
                                -1.0,
                                layer.plane,
                                layer.plane + layer.sign * metal_thickness,
                                pullback_metal
                            )
                        )
                    end
                    result = fuse_all(occ, result)
                    apply_plan_view_mask(occ, result, conductor_facets, lower, upper)
                else
                    isempty(conductor_loops) &&
                        error("Plan-view boundary is missing conductor $conductor")
                    loft_mask(
                        occ,
                        conductor_loops,
                        layer.plane,
                        layer.plane + layer.sign * metal_thickness,
                        pullback_metal,
                        tolerance
                    )
                end
                conductor_metal = if isempty(boundary_loops)
                    fillet_plane_edges(
                        occ,
                        conductor_metal,
                        top_rounding,
                        layer.plane + layer.sign * metal_thickness,
                        tolerance
                    )
                else
                    fillet_physical_edges(
                        occ,
                        conductor_metal,
                        top_rounding,
                        layer.plane + layer.sign * metal_thickness,
                        physical_segments(conductor_loops, pullback_metal, tolerance),
                        tolerance
                    )
                end
                append!(metal, conductor_metal)
            end
        end
        field, _ = occ.cut([outer], metal, -1, true, true)
        vacuum, _ = occ.cut(field, substrates, -1, true, false)
        objects = vcat(substrates, vacuum)
        if prism_tubes
            tube_tools, tube_records, tubes, tube_segments, tube_sections, tube_untubed_edges,
            tube_joints =
                build_edge_tubes!(occ, layers, boundary_loops, etch_loops, semantic_corners,
                                  edge_size, edge_growth_ratio, tube_sector_degrees,
                                  metal_thickness, overetch, corner_isotropy_radius,
                                  lc_tangent, lc_fine, lower, upper, tolerance;
                                  maximum_jacobian_condition=maximum_jacobian_condition)
        end
        domains, domain_map = occ.fragment(objects, tube_tools)
        substrate_seed = domain_map[1:length(substrates)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[(length(substrates) + 1):length(objects)] |> Iterators.flatten |> collect
        if prism_tubes
            tube_volumes = Int32[tube_volume_after_fragment(record, index, domain_map, length(objects))
                                 for (index, record) in enumerate(tube_records)]
            for (index, record) in enumerate(tube_records)
                material_seed = record.group[3] == 1 ? substrate_seed : vacuum_seed
                (3, tube_volumes[index]) in material_seed ||
                    error("Tube volume $(tube_volumes[index]) material differs from its section material")
            end
        end
    else
        vacuum, _ = occ.cut([outer], substrates, -1, true, false)
        tools = Tuple{Int32, Int32}[]
        depth = upper[3] - lower[3]
        # Thin metal under the prism tubes (decision 66): one sheet tube per straight
        # metal side at the process plane, fragmented together with the sheet faces.
        prism_tubes && isempty(boundary_loops) &&
            error("Prism tubes on thin metal require the classified plan-view boundary")
        for conductor in sort!(unique(edge.conductor for edge in edges))
            conductor_edges = [edge for edge in edges if edge.conductor == conductor]
            conductor_facets = [facet for facet in facets if facet.conductor == conductor]
            conductor_loops = [loop for loop in boundary_loops if loop.conductor == conductor]
            conductor_tools = if !isempty(conductor_loops)
                # Thin metal is a sheet, not a full-height volume partition. Supplying
                # surface tools avoids artificial vertical seams throughout the bulk.
                planar_mask_surfaces(occ, conductor_loops, tolerance)
            elseif !isempty(conductor_facets)
                faces = [(Int32(2), occ.addPlaneSurface([
                    polygon_wire(occ, facet.points, facet.plane)])) for facet in conductor_facets]
                fuse_all(occ, faces)
            else
                isempty(facets) || error("Plan-view mask is missing conductor $conductor")
                result = Tuple{Int32, Int32}[]
                for edge in conductor_edges
                    append!(result, extruded_strip(occ, edge, radius, -1.0, lower[3], depth))
                end
                fuse_all(occ, result)
            end
            append!(tools, conductor_tools)
        end
        sheet_tools = length(tools)
        if prism_tubes
            tube_tools, tube_records, tubes, tube_segments, tube_sections, tube_untubed_edges,
            tube_joints =
                build_edge_tubes!(occ, layers, boundary_loops, etch_loops, semantic_corners,
                                  edge_size, edge_growth_ratio, tube_sector_degrees,
                                  metal_thickness, overetch, corner_isotropy_radius,
                                  lc_tangent, lc_fine, lower, upper, tolerance; fabricated=false,
                                  maximum_jacobian_condition=maximum_jacobian_condition)
            append!(tools, tube_tools)
            # SEAM (design round 2): the tip bisector curves of the convex thin tips below
            # 90 degrees, fragmented with the sheet after the tube tools (the tube map
            # offsets are unchanged).
            bisector_tools, tip_bisectors, unrefined_tips = tip_bisector_tools!(
                occ, semantic_corners, corner_kinds, corner_sides, tube_segments,
                tube_sections["Radius"] + tube_sections["PyramidHeight"], corner_grading, lc_far,
                far_growth, tolerance; corner_shape_gate=corner_shape_gate)
            bisector_offset = length(tools)
            append!(tools, bisector_tools)
        end
        objects = vcat(substrates, vacuum)
        domains, domain_map = occ.fragment(objects, tools)
        isempty(tip_bisectors) ||
            bind_tip_bisector_curves!(tip_bisectors,
                                      domain_map[(length(objects) + bisector_offset + 1):(length(objects) + length(tools))],
                                      tolerance)
        substrate_seed = domain_map[1:length(substrates)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[(length(substrates) + 1):(length(substrates) + length(vacuum))] |>
            Iterators.flatten |>
            collect
        if prism_tubes
            # The tube tools follow the sheet tools in the fragment map.
            tube_volumes = Int32[tube_volume_after_fragment(record, index, domain_map,
                                                            length(objects) + sheet_tools)
                                 for (index, record) in enumerate(tube_records)]
            for (index, record) in enumerate(tube_records)
                material_seed = record.group[3] == 1 ? substrate_seed : vacuum_seed
                (3, tube_volumes[index]) in material_seed ||
                    error("Tube volume $(tube_volumes[index]) material differs from its section material")
            end
        end
    end
    surface_tools = surface_constraints === nothing ? Tuple{Int32,Int32}[] :
        surface_constraints(occ, boundary_loops, [(layer.plane,layer.sign) for layer in layers], tolerance)
    trace_tools = matching_trace === nothing ? Tuple{Int32,Int32}[] :
        matching_trace_lines(occ, matching_trace, lower, upper, tolerance;
                             mode=matching_trace_mode)
    append!(trace_tools,surface_tools)
    if !isempty(trace_tools)
        old_domains = [(dim,tag) for (dim,tag) in domains if dim==3]
        old_substrate = Set(substrate_seed); old_vacuum = Set(vacuum_seed)
        domains, trace_map = occ.fragment(old_domains, trace_tools)
        substrate_seed = Tuple{Int32,Int32}[]; vacuum_seed = Tuple{Int32,Int32}[]
        for (i,domain) in enumerate(old_domains)
            descendants = [(dim,tag) for (dim,tag) in trace_map[i] if dim==3]
            domain in old_substrate && append!(substrate_seed, descendants)
            domain in old_vacuum && append!(vacuum_seed, descendants)
        end
    end
    occ.synchronize()
    corner_point_tags = corner_isotropy ? semantic_corner_points(semantic_corners, tolerance) :
                        Int32[]

    domain_tags = Set(tag for (dim, tag) in domains if dim == 3)
    substrate_tags = sort!(
        unique(tag for (dim, tag) in substrate_seed if dim == 3 && tag in domain_tags)
    )
    vacuum_tags =
        sort!(unique(tag for (dim, tag) in vacuum_seed if dim == 3 && tag in domain_tags))
    isempty(substrate_tags) && error("No substrate volumes were generated")
    isempty(vacuum_tags) && error("No vacuum volumes were generated")
    substrate_set = Set(substrate_tags)
    vacuum_set = Set(vacuum_tags)

    matching = Int32[]
    boundary_groups = Dict{Int, Vector{Int32}}()
    # Material interfaces: interior surfaces with a substrate volume on one side
    # and a vacuum volume on the other (trench floor and walls, un-etched plane).
    interface_surfaces = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        adjacent_substrate = [volume for volume in up if volume in substrate_set]
        adjacent_vacuum = [volume for volume in up if volume in vacuum_set]
        isempty(adjacent_substrate) && isempty(adjacent_vacuum) && continue
        bounds = gmsh.model.getBoundingBox(dim, tag)
        if on_outer_box(bounds, lower, upper, outer_tolerance)
            push!(matching, tag)
            continue
        end

        # A face between two volumes of one material (a tube volume and the
        # volume around it) is an interior face.
        prism_tubes && length(up) == 2 && (isempty(adjacent_substrate) || isempty(adjacent_vacuum)) &&
            continue
        point = point_on_surface(tag)
        attribute = 0
        if fabricated
            # A fabricated surface belongs to the process layer whose band
            # [plane - Nz x Overetch, plane + Nz x MetalThickness] contains it, and is
            # owned by that layer's edges: with two layers (decision 48) the nearest edge
            # of another layer must not decide the un-etched / etched plane or the slot.
            _, _, zmin, _, _, zmax = bounds
            layer = surface_process_layer(layers, zmin, zmax, metal_thickness, overetch,
                                          outer_tolerance)
            layer_edges = layer.edges
            if !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                edge = nearest_edge(layer_edges, point, radius)
                # The un-etched plane is the interface lying in the layer plane; the
                # OCC bounding box carries the CAD's absolute tolerance (1e-7), so it is
                # compared with the box tolerance like every other bounding box here
                # (the source tolerance 1e-7 x Radius is below it for Radius < 1 um).
                attribute =
                    abs(zmin - layer.plane) < outer_tolerance &&
                    abs(zmax - layer.plane) < outer_tolerance ? 3000 + edge.slot :
                    3100 + edge.slot
                push!(interface_surfaces, tag)
            elseif !isempty(adjacent_substrate)
                owner = nearest_metal_edge(layer_edges, facets, point, radius, tolerance;
                                           surface=(tag, bounds, up))
                attribute = metal_surface_attribute(5000, owner.slot, owner.conductor)
            elseif !isempty(adjacent_vacuum)
                owner = nearest_metal_edge(layer_edges, facets, point, radius, tolerance;
                                           surface=(tag, bounds, up))
                attribute = metal_surface_attribute(6000, owner.slot, owner.conductor)
            end
        else
            edge = nearest_edge(edges, point, radius)
            metal_edges = [
                candidate for candidate in edges if
                point_in_metal(candidate, point, radius, tolerance, facets)
            ]
            if !isempty(metal_edges)
                owner = nearest_edge(metal_edges, point, radius)
                attribute = metal_surface_attribute(4000, owner.slot, owner.conductor)
            elseif !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                attribute = 3000 + edge.slot
                push!(interface_surfaces, tag)
            end
        end
        attribute > 0 && push!(get!(boundary_groups, attribute, Int32[]), tag)
    end
    isempty(matching) && error("No matching surface was generated")

    if !prism_tubes
        gmsh.model.addPhysicalGroup(3, substrate_tags, 1, "substrate")
        gmsh.model.addPhysicalGroup(3, vacuum_tags, 2, "vacuum")
    end
    gmsh.model.addPhysicalGroup(2, unique(matching), 1, "matching_surface")
    for (attribute, surfaces) in sort(collect(boundary_groups))
        unique!(surfaces)
        gmsh.model.addPhysicalGroup(2, surfaces, attribute, "surface_$attribute")
    end

    if labels_only !== nothing
        # Decision 62(2): the registration probe learns the label set here - every
        # physical surface group exists on the CAD entities after the last occ.fragment
        # (the trace-basis fragments included) - so the pass stops before any mesh
        # generation.  The record carries the fields the derivation and the registration
        # consume: InterfaceAreas[].Attribute / Name (no areas: nothing is meshed), the
        # recipe scope record and the semantic contract identity.
        ispath(labels_only) && error("Labels-only output already exists")
        # Step 4.3: the probe also matches every tube's structural entities against the
        # fragmented CAD (match_tube_entities, fail closed), so a fragment that splits a tube
        # entity - the first thing a production build would stop on - is caught at the
        # registration probe (the single-descendant check alone let the arc tubes through).
        matched_tube_entities = 0
        for (index, record) in enumerate(tube_records)
            matched_tube_entities += length(match_tube_entities(
                tube_volumes[index], record.tube, record.section, record.group, 1.0e-6 * radius))
        end
        label_rows = [Dict{String, Any}("Attribute" => 1, "Name" => "matching_surface",
                                        "CADSurfaces" => length(unique(matching)))]
        for (attribute, surfaces) in sort(collect(boundary_groups))
            push!(label_rows, Dict{String, Any}("Attribute" => attribute, "Name" => "surface_$attribute",
                                                "CADSurfaces" => length(surfaces)))
        end
        open(labels_only, "w") do stream
            write_json(stream, Dict{String, Any}(
                "Version" => 1, "Frame" => "SourceLocal", "LabelsOnly" => true,
                "Purpose" => "Labels-only pass of the production-option build (decision 62(2)): the " *
                             "physical surface groups assigned on the CAD entities, before any mesh " *
                             "generation; the registration probe's census (InterfaceAreas carries the " *
                             "label set, no areas - nothing is meshed)",
                "Scope" => prism_tubes ?
                    recipe_scope_record(scope_classes, boundary_loops, lower, upper, tolerance,
                                        tube_untubed_edges) :
                    nothing,
                "SemanticContract" => semantic_contract,
                "SemanticContractSHA256" => bytes2hex(sha256(read(semantic_contract))),
                "SemanticCorners" => [collect(corner) for corner in semantic_corners],
                "CouponBox" => Dict{String, Any}("Radius" => radius, "Lower" => collect(lower),
                                                 "Upper" => collect(upper)),
                # The face ends of the tubes (design A2; decision 320), known at the CAD
                # stage: the registration probe records them before any mesh exists.
                "PrismTubeEntitiesMatched" => matched_tube_entities,
                "PrismTubeFaceEnds" => prism_tubes ?
                    Dict{String, Any}(
                        "Count" => sum(length(tube.face_ends) for (tube, _) in tubes; init=0),
                        # The section's end-spacing cap lc_cap (design round 2 F2b; null without a face end).
                        "SpacingCap" => tube_sections["FaceEndSpacingCap"],
                        "Tubes" => [Dict{String, Any}(
                                        "Tube" => k,
                                        "Origin" => tube isa ArcTube ? arc_frame(tube, tube.s_start).origin : tube.origin,
                                        "Extrusion" => tube isa ArcTube ? arc_frame(tube, tube.s_start).e : tube.e,
                                        "FaceEnds" => [face_end_record(f) for f in tube.face_ends])
                                    for (k, (tube, _)) in enumerate(tubes) if !isempty(tube.face_ends)],
                        "LegacyBoxVertexCorners" => [point for segment in tube_segments
                                                     for (point, legacy) in ((segment.start, segment.legacy_box_corners[1]),
                                                                             (segment.stop, segment.legacy_box_corners[2]))
                                                     if legacy]) : nothing,
                "InterfaceAreas" => label_rows))
            println(stream)
        end
        println("Labels-only census: $labels_only")
        for row in label_rows
            println("  label $(row["Attribute"]) $(row["Name"]): CAD surfaces=$(row["CADSurfaces"])")
        end
        println("Spatial coupon labels only: fabricated=$fabricated, edges=$(length(edges)), " *
                "layers=$(length(layers)), labels=$(length(label_rows))")
        return gmsh.finalize()
    end

    # Build the refinement source from physical process edges, not every OCC fragment
    # curve. Exact plan-view masks and slot-partitioned interface surfaces can introduce
    # thousands of coplanar bookkeeping seams; treating those as fabrication features
    # makes the nanometer-scale size field cover most of a large coupon.
    curve_surfaces = Dict{Int32, Vector{Tuple{Int, Int32}}}()
    candidate_curves = Set{Int32}()
    for (attribute, surfaces) in boundary_groups
        for surface in surfaces
            for (curve_dim, curve) in
                gmsh.model.getBoundary([(2, surface)], false, false, false)
                curve_dim == 1 || continue
                on_outer_box(
                    gmsh.model.getBoundingBox(curve_dim, curve),
                    lower,
                    upper,
                    outer_tolerance
                ) && continue
                push!(candidate_curves, curve)
                push!(get!(curve_surfaces, curve, Tuple{Int, Int32}[]), (attribute, surface))
            end
        end
    end
    feature_curves = Int32[]
    discarded_seams = Int32[]
    for curve in candidate_curves
        records = unique(curve_surfaces[curve])
        attributes = unique(first(record) for record in records)
        surfaces = unique(last(record) for record in records)
        families = unique(surface_family(attribute) for attribute in attributes)
        artificial_seam = length(surfaces) >= 2 && length(families) == 1 &&
                          coplanar_surfaces(surfaces, tolerance)
        push!(artificial_seam ? discarded_seams : feature_curves, curve)
    end
    sort!(feature_curves)
    isempty(feature_curves) && error("No physical process-feature curves were generated")
    # The lines where a material interface meets the outer box (the cut/trench and
    # cut/un-etched-plane junctions) are feature curves with the same band as the
    # process edges: the physics pilot located every AMR mark within 0.3 um of
    # them. They are curves of interface surfaces lying on the outer box.
    junction_curves = Int32[]
    for surface in interface_surfaces
        for (curve_dim, curve) in gmsh.model.getBoundary([(2, surface)], false, false, false)
            curve_dim == 1 || continue
            on_outer_box(gmsh.model.getBoundingBox(curve_dim, curve), lower, upper,
                         outer_tolerance) || continue
            push!(junction_curves, curve)
        end
    end
    sort!(unique!(junction_curves))
    isempty(junction_curves) &&
        error("No cut-surface/material-interface junction curves were generated")
    junction_length = sum(gmsh.model.occ.getMass(1, curve) for curve in junction_curves)
    append!(feature_curves, junction_curves)
    sort!(unique!(feature_curves))
    # Tube entities: their curves are meshed explicitly (never placed or used as
    # attractors), the tube axes carry the band grading around the tubes, and the
    # cap centres before a corner are graded points of the corner law.
    tube_states = TubeMesh[]
    tube_curves = Set{Int32}()
    tube_axis_curves = Int32[]
    tube_cap_points = Int32[]
    tube_volume_groups = Vector{Tuple{Int32, Tuple{Int, Int, Int}, Dict}}[]
    # The ends of every tube that are smooth joints (design A3 (2)): no cap there, no graded
    # cap centre; the shared section is installed once (the owner's), adopted by the other.
    tube_joint_ends = [Set{Int}() for _ in tubes]
    if prism_tubes
        for (adopting, own_end, owning, owner_end) in tube_joints
            push!(tube_joint_ends[adopting], own_end)
            push!(tube_joint_ends[owning], owner_end)
        end
        for (k, (tube, section)) in enumerate(tubes)
            volumes = Tuple{Int32, Tuple{Int, Int, Int}, Dict}[]
            for (index, record) in enumerate(tube_records)
                record.tube === tube || continue
                matched = match_tube_entities(tube_volumes[index], tube, section, record.group,
                                              1.0e-6 * radius)
                push!(volumes, (tube_volumes[index], record.group, matched))
                for ((kind, id), tag) in matched
                    kind in (:edge_line, :outer_line, :cap_polygon, :cap_ray) && push!(tube_curves, tag)
                    kind == :edge_line && push!(tube_axis_curves, tag)
                    if kind == :edge_point
                        id[2] in tube_joint_ends[k] && continue
                        point = gmsh.model.getValue(0, tag, Float64[])
                        on_outer_box(vcat(point, point), lower, upper, outer_tolerance) ||
                            push!(tube_cap_points, tag)
                    end
                end
            end
            isempty(volumes) && error("A tube has no fragmented volume")
            push!(tube_volume_groups, volumes)
        end
        sort!(unique!(tube_axis_curves))
        sort!(unique!(tube_cap_points))
        isempty(tube_cap_points) && error("No tube ends before a semantic corner")
    end
    graded_points = copy(semantic_corners)
    graded_point_tags = copy(corner_point_tags)
    for tag in tube_cap_points
        push!(graded_points, Tuple(gmsh.model.getValue(0, tag, Float64[])))
        push!(graded_point_tags, tag)
    end
    longitudinal_curves = Int32[]
    corner_curves = Int32[]
    corner_grading_slope = (lc_far - lc_fine) / (process_core_width - process_fine_width)
    # Under the prism tube recipe the trace rule is a volume law growing with
    # FarGrowth from the cut surface (decision 39a); the legacy seed keeps the
    # process-band slope.
    trace_basis_record = trace_basis === nothing ? nothing :
        prepare_trace_basis_sizing!(trace_basis, trace_basis_size_ratio, lc_far,
                                    prism_tubes ? far_growth : corner_grading_slope)
    # The band grading around the tubes and the band law of the remaining feature
    # curves (junction lines, footprint edges, the un-tubed edge parts) are exact
    # segment-distance laws evaluated in the size callback with the trace rule.
    tube_band_record = prism_tubes ?
        prepare_tube_band_sizing!(curve_segments(tube_axis_curves; subdivide=true),
                                  curve_segments([curve for curve in feature_curves
                                                  if !(curve in tube_curves)]; subdivide=true),
                                  tube_sections["Radius"] + tube_sections["PyramidHeight"],
                                  lc_fine, lc_far, far_growth, graded_points,
                                  corner_isotropy_radius) :
        nothing
    # Tube layers follow the composed size field on the axis (decision 40): the
    # tubes are re-stationed with the laws prepared above before their mesh is
    # installed; the cross-section rings are unchanged. At an end on the outer box
    # the end layer is the size at the surface (decision 41: the narrow hats decay
    # from the cut surface at the trace rule's size there).
    tube_layer_records = Dict{String, Any}[]
    if prism_tubes
        for (k, (tube, section)) in enumerate(tubes)
            on_box(s) = on_outer_box(vcat(tube_point(tube, 0.0, 0.0, s), tube_point(tube, 0.0, 0.0, s)),
                                     lower, upper, outer_tolerance)
            axis_size(s) = tube_axis_size(tube_point(tube, 0.0, 0.0, s), lc_tangent, semantic_corners,
                                          corner_grading, corner_grading_slope)
            graded = if isempty(tube.face_ends)
                stations, axis_positions, axis_sizes = graded_tube_stations(
                    tube.s_start, tube.s_end, axis_size, lc_tangent, edge_growth_ratio;
                    surface_start=on_box(tube.s_start), surface_end=on_box(tube.s_end))
                tube isa ArcTube ? ArcTube(tube, stations) : EdgeTube(tube, stations)
            elseif tube isa ArcTube
                stations, fraction, block_end, axis_positions, axis_sizes = face_ended_tube_stations(
                    tube, axis_size, lc_tangent, edge_growth_ratio;
                    surface_start=on_box(tube.s_start), surface_end=on_box(tube.s_end))
                ArcTube(tube, stations; block_fraction=fraction, block_end=block_end)
            else
                # A face end replaces the surface layer of its end by the sheared end
                # block (design A2); the other end keeps the surface rule.
                stations, shear_u, shear_w, axis_positions, axis_sizes = face_ended_tube_stations(
                    tube, axis_size, lc_tangent, edge_growth_ratio;
                    surface_start=on_box(tube.s_start), surface_end=on_box(tube.s_end))
                EdgeTube(tube, stations; shear_u=shear_u, shear_w=shear_w)
            end
            tubes[k] = (graded, section)
            push!(tube_layer_records, tube_layer_statistics(graded, axis_positions, axis_sizes))
            push!(tube_states, TubeMesh(graded, section, tube_volume_groups[k];
                                        pyramid_height=tube_sections["PyramidHeight"]))
        end
        for (adopting, own_end, owning, owner_end) in tube_joints
            register_joint!(tube_states[adopting], own_end, tube_states[owning], owner_end)
        end
    end
    # Installation order (design A3 (2)): the owner of every shared section before the tube
    # adopting it - the arc tubes (parts in arc order) before the straight tubes.
    tube_install_order = vcat([i for (i, state) in enumerate(tube_states) if state.tube isa ArcTube],
                              [i for (i, state) in enumerate(tube_states) if !(state.tube isa ArcTube)])
    for (adopting, _, owning, _) in tube_joints
        findfirst(==(owning), tube_install_order) < findfirst(==(adopting), tube_install_order) ||
            error("a tube adopts a shared section from a tube installed after it")
    end
    tube_meshed_curves = Set{Int32}()
    tube_meshed_faces = Set{Int32}()
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    tube_curve_meshes = Vector{Dict{String, Any}}(undef, length(tube_states))
    for i in tube_install_order
        tube_curve_meshes[i] = install_tube_curves!(tube_states[i], next_node, point_nodes;
                                                    meshed=tube_meshed_curves)
    end
    # Interior ridge node parameters and coordinates (curve order) of every
    # longitudinal curve; the mesh is assigned after the edge layer curves are
    # known, because a layer ridge is subdivided inside its span.
    ridge_parameters = Dict{Int32, Vector{Float64}}()
    ridge_nodes = Dict{Int32, Matrix{Float64}}()
    transfinite_curves = Set{Int32}()
    # Band curves (Gmsh-only recipe): the non-metal longitudinal feature curves -
    # the cut-surface / material-interface junction lines and the footprint edges
    # parallel to a metal edge - are 1D-meshed at the band law's size on the line,
    # NormalSize (BandRule at r = 0; the growth away from the line is the size
    # callback's), composed with the corner law, so that the first cell layer
    # against the line is NormalSize transversally: on the lc_tangent ridge grid
    # the first layer was TangentialSize (the surface mesher cannot refine a curve's
    # nodes). The metal ridges keep the lc_tangent grid: they carry the tubes.
    # Every explicitly placed curve follows min(its spacing, the composed size
    # field on the curve) - CURVE_SPACING_RULE (decision 43): a basis sliver
    # crossing a band curve is resolved at the trace rule along the curve too.
    junction_set = Set(junction_curves)
    band_curves = Int32[]
    curve_spacing_records = Dict{String, Any}[]
    function curve_size_law(point, spacing)
        size = corner_isotropy ?
            corner_curve_size(point, graded_points, corner_grading, spacing, corner_grading_slope) :
            spacing
        trace_basis_record === nothing ||
            (size = trace_basis_size(point[1], point[2], point[3], size))
        tube_band_record === nothing ||
            (size = feature_band_size(point[1], point[2], point[3], size))
        return size
    end
    if lc_tangent > 0.0
        for curve in feature_curves
            curve in tube_curves && continue
            lower_parameter, upper_parameter =
                gmsh.model.getParametrizationBounds(1, curve)
            parameter = 0.5 * (lower_parameter[1] + upper_parameter[1])
            derivative = gmsh.model.getDerivative(1, curve, [parameter])
            tangent = derivative[1:3]
            tangent ./= norm(tangent)
            # A circle feature curve of an arc coupon (a metal arc ridge or a concentric
            # trench / collar arc) is longitudinal like the straight ridge of a straight
            # side (design 1.2 (5)); its tangent at the midpoint matches no chord row.
            arc_curve = arc_boundary && gmsh.model.getType(1, curve) != "Line"
            if arc_curve || any(abs(dot(tangent, edge.tangent)) >= 1.0 - 1.0e-6 for edge in edges)
                push!(longitudinal_curves, curve)
                band_curve = prism_tubes && (curve in junction_set ||
                    !any(surface_family(first(record)) in METAL_SURFACE_FAMILIES
                         for record in get(curve_surfaces, curve, Tuple{Int, Int32}[])))
                band_curve && push!(band_curves, curve)
                curve_spacing = band_curve ? lc_fine : lc_tangent
                placed = composed_curve_nodes(
                    curve, curve_spacing, point -> curve_size_law(point, curve_spacing),
                    edge_growth_ratio, corner_isotropy ? graded_points : nothing, corner_grading)
                curve_record = Dict{String, Any}(
                    "Curve" => Int(curve),
                    "Kind" => curve in junction_set ? "junction" : band_curve ? "band" : "metal",
                    "Segment" => vcat(gmsh.model.getValue(1, curve, [lower_parameter[1]]),
                                      gmsh.model.getValue(1, curve, [upper_parameter[1]])))
                if placed === nothing
                    curve_length = gmsh.model.occ.getMass(1, curve)
                    point_count = max(2, ceil(Int, curve_length / curve_spacing) + 1)
                    grid = collect(range(lower_parameter[1], upper_parameter[1];
                                         length=point_count))[2:(end - 1)]
                    ridge_parameters[curve] = grid
                    ridge_nodes[curve] = reshape(gmsh.model.getValue(1, curve, grid), 3, :)
                    push!(transfinite_curves, curve)
                    uniform = curve_length / (point_count - 1)
                    merge!(curve_record, Dict{String, Any}(
                        "Length" => curve_length, "Spacing" => curve_spacing, "Graded" => false,
                        "GridIntervals" => point_count - 1, "GridIntervalsKept" => point_count - 1,
                        "InteriorNodes" => length(grid),
                        "NodeSpacing" => Dict{String, Any}("Minimum" => uniform, "P50" => uniform,
                                                           "Maximum" => uniform),
                        "PrescribedMinimum" => curve_spacing,
                        "AchievedOverPrescribed" => Dict{String, Any}(
                            "Minimum" => uniform / curve_spacing, "P50" => uniform / curve_spacing,
                            "Maximum" => uniform / curve_spacing)))
                else
                    ridge_parameters[curve] = placed[1]
                    ridge_nodes[curve] = reshape(placed[2], 3, :)
                    merge!(curve_record, placed[3])
                end
                push!(curve_spacing_records, curve_record)
            end
        end
    end
    # Edge layer: on every face bounding a metal edge (a longitudinal feature
    # curve of a metal surface family, junction curves excluded), one embedded row
    # per geometric layer along the ridge span on the lc_tangent grid outside the
    # corner size law. Inside the span the ridge and the rows are refined
    # tangentially so that every layer cell keeps an aspect of at most
    # EdgeLayerAspect (a tetrahedron with one corner whose three edges are all
    # tangential has a scaled Jacobian of (hn/ht)^2, so the scaled-Jacobian gate
    # bounds the layer anisotropy): the ridge grid interval is subdivided by
    # TangentialSubdivision = the smallest power of two with
    # lc_tangent / n <= EdgeLayerAspect x EdgeSize, row k by the nested power of
    # two with lc_tangent / m <= EdgeLayerAspect x (its layer size), and the
    # ridge subdivision halves interval by interval towards the span ends (the
    # taper), where no row is seeded. Prisms would not have the aspect limit;
    # they were not pursued.
    edge_layer_curves = Dict{String, Any}[]
    edge_layer_records = Tuple{Dict{String, Any}, Matrix{Float64}}[]
    edge_layer_subdivision = 0
    edge_layer_row_subdivisions = Int[]
    edge_layer_taper = Int[]
    edge_layer_ridge_nodes = 0
    layer_to_ball = edge_layer && corner_size > 0.0
    if edge_layer
        edge_layer_subdivision = tangential_subdivision(lc_tangent, edge_layer_aspect * edge_size)
        edge_layer_row_subdivisions = [
            min(edge_layer_subdivision,
                tangential_subdivision(lc_tangent,
                                       edge_layer_aspect * edge_size * edge_growth_ratio^(k - 1)))
            for k in eachindex(edge_layer_offsets)]
        edge_layer_taper = layer_to_ball ? Int[] :
            [2^j for j in 1:(round(Int, log2(edge_layer_subdivision)) - 1)]
        for curve in longitudinal_curves
            curve in junction_set && continue
            records = get(curve_surfaces, curve, Tuple{Int, Int32}[])
            any(surface_family(first(record)) in METAL_SURFACE_FAMILIES
                for record in records) || continue
            xyz = ridge_nodes[curve]
            corner_distance = [minimum(norm(xyz[:, i] .- collect(corner))
                                       for corner in semantic_corners) for i in axes(xyz, 2)]
            # Without corner grading the rows lie on the lc_tangent grid outside the
            # corner size law, after the taper intervals; with corner grading
            # (decision 33) the span starts at the ball boundary node (the graded
            # ball is the transition: ~lc_fine ridge spacing at the radius against
            # the first row's subdivision), so no ridge interval is left without rows
            # outside the ball and no taper is needed.
            grid_indices = layer_to_ball ?
                findall(>=(corner_isotropy_radius * (1.0 - 1.0e-9)), corner_distance) :
                findall([corner_curve_size(Tuple(xyz[:, i]), semantic_corners, corner_grading,
                                           lc_tangent, corner_grading_slope) >= lc_tangent
                         for i in axes(xyz, 2)])
            taper = length(edge_layer_taper)
            length(grid_indices) >= 2taper + 2 || continue
            grid_indices == collect(grid_indices[1]:grid_indices[end]) ||
                error("Longitudinal curve $curve has a non-contiguous edge layer span")
            row_indices = grid_indices[(1 + taper):(end - taper)]
            span = xyz[:, row_indices]
            unlayered = (corner_distance[row_indices[1]], corner_distance[row_indices[end]])
            # The curve's CAD endpoints in span order (the span runs with the curve's
            # interior node order, from the lower parameter end).
            curve_bounds = gmsh.model.getParametrizationBounds(1, curve)
            curve_xyz = gmsh.model.getValue(1, curve, [curve_bounds[1][1], curve_bounds[2][1]])
            curve_ends = [collect(curve_xyz[1:3]), collect(curve_xyz[4:6])]
            faces = Dict{String, Any}[]
            for face in sort!(unique(last(record) for record in records))
                record = add_edge_layer_rows!(occ, face, span, edge_layer_offsets, tolerance)
                record === nothing && continue
                record["RowNodes"] = sum(size(edge_layer_row_nodes(span, offset, record["Inward"],
                                                                   subdivision), 2)
                                         for (offset, subdivision) in
                                             zip(edge_layer_offsets[1:length(record["Rows"])],
                                                 edge_layer_row_subdivisions))
                push!(faces, record)
                push!(edge_layer_records, (record, span))
            end
            isempty(faces) && continue
            # The ridge mesh: grid nodes everywhere, the span intervals subdivided
            # (taper at both ends, full subdivision under the rows).
            parameters = ridge_parameters[curve]
            refined = Float64[]
            for i in eachindex(parameters)
                push!(refined, parameters[i])
                i in grid_indices[1]:(grid_indices[end] - 1) || continue
                position = i - grid_indices[1] + 1
                from_end = grid_indices[end] - i
                n = position <= taper ? edge_layer_taper[position] :
                    from_end <= taper ? edge_layer_taper[from_end] : edge_layer_subdivision
                for j in 1:(n - 1)
                    push!(refined, parameters[i] + j / n * (parameters[i + 1] - parameters[i]))
                end
            end
            ridge_parameters[curve] = refined
            edge_layer_ridge_nodes += length(refined) - length(parameters)
            push!(edge_layer_curves, Dict{String, Any}(
                "Curve" => Int(curve), "Start" => collect(span[:, 1]),
                "End" => collect(span[:, end]), "SpanLength" => norm(span[:, end] .- span[:, 1]),
                "CurveLength" => gmsh.model.occ.getMass(1, curve),
                "GridNodes" => size(span, 2), "TaperIntervalsPerEnd" => taper,
                "CurveEnds" => curve_ends,
                "UnlayeredLengthAtEnds" => collect(unlayered),
                "RidgeNodesAdded" => length(refined) - length(parameters),
                "Faces" => [Dict{String, Any}("Face" => face["Face"],
                                              "Rows" => length(face["Rows"]),
                                              "RowNodes" => face["RowNodes"])
                            for face in faces]))
        end
        isempty(edge_layer_records) && error("No edge layer row fits on any metal-edge face")
        occ.synchronize()
    end
    layer_curve_set = Set{Int32}(Int32(row["Curve"]) for row in edge_layer_curves)
    for curve in longitudinal_curves
        if curve in transfinite_curves && !(curve in layer_curve_set)
            gmsh.model.mesh.setTransfiniteCurve(curve, length(ridge_parameters[curve]) + 2)
        else
            # Transfinite spacing cannot follow the corner ball or the layer
            # subdivision; place the curve nodes explicitly and keep them through
            # generation.
            parameters = ridge_parameters[curve]
            add_explicit_curve_mesh!(curve, parameters, gmsh.model.getValue(1, curve, parameters),
                                     next_node, point_nodes)
            push!(corner_curves, curve)
        end
    end
    if edge_layer
        for (record, span) in edge_layer_records
            gmsh.model.mesh.embed(1, record["Rows"], 2, record["Face"])
            mesh_edge_layer_rows!(record, span, edge_layer_offsets, edge_layer_row_subdivisions,
                                  next_node, point_nodes)
        end
        println("Edge layer: edge_size=$edge_size, growth_ratio=$edge_growth_ratio, " *
                "aspect=$edge_layer_aspect, layers=$(length(edge_layer_offsets)), " *
                "row_offsets=$(edge_layer_offsets), " *
                "tangential_subdivision=$edge_layer_subdivision, " *
                "row_subdivisions=$edge_layer_row_subdivisions, taper=$edge_layer_taper, " *
                "curves=$(length(edge_layer_curves)), " *
                "rows=$(sum(length(record["Rows"]) for (record, _) in edge_layer_records)), " *
                "row_nodes=$(sum(record["RowNodes"] for (record, _) in edge_layer_records)), " *
                "ridge_nodes_added=$edge_layer_ridge_nodes, " *
                "span_length=$(sum(row["SpanLength"] for row in edge_layer_curves))")
    end
    isempty(corner_curves) && !prism_tubes || gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    prism_tubes && gmsh.option.setNumber("Mesh.Renumber", 0)
    # The sizes below lc_fine a field can ask for are the trace rule's and the corner
    # grading's, so the Gmsh floor follows the smallest requested trace size
    # (recorded in the census) and CornerSize.
    mesh_size_minimum = trace_basis_record === nothing ? lc_fine :
        min(lc_fine, trace_basis_record["MinimumRequestedSize"])
    corner_size > 0.0 && (mesh_size_minimum = min(mesh_size_minimum, corner_size))
    prism_tubes && (mesh_size_minimum = min(mesh_size_minimum, edge_size))
    # The trace rule and the tube/band laws (prepared before the tube layers were
    # placed) compose with the background field in the size callback.
    if trace_basis_record !== nothing || tube_band_record !== nothing
        install_trace_basis_callback!()
    end
    if trace_basis_record !== nothing
        trace_basis_record["MeshSizeMinimum"] = mesh_size_minimum
        println("Trace basis sizing: ratio=$(trace_basis_size_ratio), triangles=" *
                "$(trace_basis_record["Triangles"]), below far=$(trace_basis_record["TrianglesBelowFarSize"]), " *
                "basis edges below far=$(trace_basis_record["BasisEdgesBelowFarSize"]), " *
                "minimum requested=$(trace_basis_record["MinimumRequestedSize"]), " *
                "mesh size minimum=$mesh_size_minimum")
    end
    println(
        "Spatial mesh features: candidates=$(length(candidate_curves)), " *
        "physical=$(length(feature_curves)), " *
        "junction=$(length(junction_curves)) (length $junction_length), " *
        "longitudinal=$(length(longitudinal_curves)), " *
        "band_curves=$(length(band_curves)), " *
        "corner_isotropic_longitudinal=$(length(corner_curves)), " *
        "discarded_coplanar_seams=$(length(discarded_seams)), " *
        "fine_width=$process_fine_width, core_width=$process_core_width, " *
        "grading_power=$process_grading_power"
    )
    if lc_tangent > 0.0
        # Gmsh's curve-attractor metric resolves the process cross-section normally while
        # keeping the smooth longitudinal edge direction coarse. The metric follows the
        # nearest curve tangent, so differently oriented cluster edges do not require a
        # single global frame; intersections naturally receive the stricter local metric.
        gmsh.model.mesh.field.add("AttractorAnisoCurve", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList",
                                         Float64.([curve for curve in feature_curves
                                                   if !(curve in tube_curves)]))
        gmsh.model.mesh.field.setNumber(1, "DistMin", process_fine_width)
        gmsh.model.mesh.field.setNumber(1, "DistMax", process_core_width)
        gmsh.model.mesh.field.setNumber(1, "SizeMinNormal", lc_fine)
        gmsh.model.mesh.field.setNumber(1, "SizeMaxNormal", lc_far)
        gmsh.model.mesh.field.setNumber(1, "SizeMinTangent", lc_tangent)
        gmsh.model.mesh.field.setNumber(1, "SizeMaxTangent", lc_far)
        gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
        background = 1
        if corner_isotropy
            # MathEval evaluates the attractor's own size ("F1") unchanged outside
            # the corner balls; Min/MinAniso wrappers would not (they re-derive the
            # anisotropic child's size), so the band prescription stays identical.
            gmsh.model.mesh.field.add("Distance", 2)
            gmsh.model.mesh.field.setNumbers(2, "PointsList", Float64.(graded_point_tags))
            gmsh.model.mesh.field.add("MathEval", 3)
            gmsh.model.mesh.field.setString(3, "F", "min(F1," * corner_size_expression(
                "F2", corner_grading, lc_far, process_core_width - process_fine_width) * ")")
            background = 3
        end
        gmsh.model.mesh.field.setAsBackgroundMesh(background)
    else
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(feature_curves))
        # Isotropic fallback with immediate power-law grading away from process edges.
        transition_width = process_core_width - process_fine_width
        distance_expression = "max(F1-$(process_fine_width),0)"
        size_expression =
            "min($(lc_far),$(lc_fine)+($(lc_far)-$(lc_fine))*" *
            "($(distance_expression)/$(transition_width))^$(process_grading_power))"
        if corner_isotropy
            gmsh.model.mesh.field.add("Distance", 3)
            gmsh.model.mesh.field.setNumbers(3, "PointsList", Float64.(corner_point_tags))
            size_expression = "min($(size_expression)," * corner_size_expression(
                "F3", corner_grading, lc_far, transition_width) * ")"
        end
        gmsh.model.mesh.field.add("MathEval", 2)
        gmsh.model.mesh.field.setString(2, "F", size_expression)
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
    end
    for (name, value) in [
        ("Mesh.MeshSizeMin", mesh_size_minimum),
        ("Mesh.MeshSizeMax", lc_far),
        ("Mesh.Algorithm3D", 1),
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MinimumCirclePoints", 24),
        ("Mesh.MinimumCurvePoints", 3),
        ("Mesh.MshFileVersion", 2.2),
        ("Mesh.Binary", 1)
    ]
        gmsh.option.setNumber(name, value)
    end
    if mesh_control !== nothing
        # All controls receive matching entity tags separately from physical interfaces.
        mesh_control(feature_curves, boundary_groups, matching, lower, upper)
    end
    if geometry_only
        return gmsh.finalize()
    end
    tube_timings = Dict{String, Float64}()
    tube_face_meshes = Dict{String, Any}[]
    tube_volume_census = Dict{String, Any}[]
    tube_discrete = Dict{Int, Vector{Int32}}()
    if prism_tubes
        # Three phases around the generator (prism_edge_tubes.jl): the tube curves
        # are meshed, Gmsh meshes the surfaces, the tube faces are replaced by the
        # explicit cap triangles / radial and pyramid faces, the OCC tube volumes
        # are removed and Gmsh meshes the remaining volumes with tetrahedra, then
        # discrete volumes receive the prisms and pyramids.
        # The revolved faces of the arc tubes (design 1.2 (3)) are periodic OCC surfaces that
        # Gmsh's surface mesher refuses with the explicit boundary nodes ("Impossible to mesh
        # periodic surface"); their meshes are replaced by the explicit prism faces anyway, so
        # on an arc coupon every tube face is hidden from the 2D pass (Mesh.MeshOnlyVisible)
        # and shown again before the tube faces are installed. A straight coupon keeps the
        # unchanged pass (its planar tube faces are meshed and replaced as before: bitwise).
        hidden_tube_faces = Int32[]
        if any(state.tube isa ArcTube for state in tube_states)
            hidden_tube_faces = unique(reduce(vcat, [reduce(vcat, values(state.faces); init=Int32[])
                                                      for state in tube_states]; init=Int32[]))
            gmsh.model.setVisibility([(2, face) for face in hidden_tube_faces], 0)
            gmsh.option.setNumber("Mesh.MeshOnlyVisible", 1)
        end
        tube_timings["Generate2D"] = @elapsed gmsh.model.mesh.generate(2)
        tube_timings["TubeFaces"] = @elapsed begin
            tube_face_meshes = Vector{Dict{String, Any}}(undef, length(tube_states))
            for i in tube_install_order
                tube_face_meshes[i] = install_tube_faces!(tube_states[i], next_node;
                                                          meshed=tube_meshed_faces)
            end
        end
        remove_tube_volumes!(tube_states)
        # The 3D pass repeats the surface pass for the faces it finds pending (the explicitly
        # meshed periodic faces among them), so the arc coupon's tube faces stay hidden
        # through it; their explicit meshes bound the tetrahedra as on a straight coupon.
        tube_timings["Generate3D"] = @elapsed gmsh.model.mesh.generate(3)
        if !isempty(hidden_tube_faces)
            gmsh.model.setVisibility([(2, face) for face in hidden_tube_faces], 1)
            gmsh.option.setNumber("Mesh.MeshOnlyVisible", 0)
        end
        before_duplicates = sum(length(tags) for tags in gmsh.model.mesh.getElements(3)[2]; init=0)
        gmsh.model.mesh.removeDuplicateElements([(3, tag) for (dim, tag) in gmsh.model.getEntities(3)])
        after_duplicates = sum(length(tags) for tags in gmsh.model.mesh.getElements(3)[2]; init=0)
        before_duplicates == after_duplicates ||
            error("Gmsh produced $(before_duplicates - after_duplicates) duplicate volume elements")
        tube_discrete, tube_volume_census = finalize_tube_volumes!(tube_states[tube_install_order])
        gmsh.model.addPhysicalGroup(3, vcat(setdiff(substrate_tags, tube_volumes),
                                            get(tube_discrete, 1, Int32[])), 1, "substrate")
        gmsh.model.addPhysicalGroup(3, vcat(setdiff(vacuum_tags, tube_volumes),
                                            get(tube_discrete, 2, Int32[])), 2, "vacuum")
    else
        gmsh.model.mesh.generate(3)
    end
    (trace_basis_record === nothing && tube_band_record === nothing) ||
        gmsh.model.mesh.removeSizeCallback()
    # Reject oversized linear meshes before allocating their high-order nodes.
    _, linear_tags, _ = gmsh.model.mesh.getElements(3)
    sum(length(tags) for tags in linear_tags) <= max_elements ||
        error("Linear spatial coupon exceeds element budget before order elevation")
    if lc_tangent == 0.0 && optimize_volume
        gmsh.model.mesh.optimize("Netgen")
    end
    seed_quality = seed_quality_gates ?
        optimize_required_seed_region!(
            semantic_corners, corner_isotropy_radius, lc_fine,
            [(Vector{Float64}(row["Start"]), Vector{Float64}(row["End"]))
             for row in edge_layer_curves],
            edge_size, edge_growth_ratio,
            isempty(edge_layer_offsets) ? 0.0 : edge_layer_offsets[end],
            EDGE_LAYER_ROW_ZIGZAG, maximum_corner_aspect, minimum_scaled_jacobian,
            maximum_jacobian_condition, quality_displacement_over_normal, tolerance,
            edge_layer_maximum_aspect,
            corner_grading; fixed_node_tags=prism_tubes ? non_simplex_volume_nodes() : UInt[],
            corner_kinds=corner_kinds, corner_sides=corner_sides,
            corner_shape_gate=corner_shape_gate) :
        nothing
    # Per-type quality of the mixed mesh, gated after the corner-ball optimization
    # (fail closed): positive orientation and Jacobian condition for every type,
    # scaled Jacobian for the tetrahedra.
    mixed_quality = prism_tubes ? mixed_volume_quality() : nothing
    prism_tubes && gate_mixed_volume_quality(mixed_quality, minimum_scaled_jacobian,
                                             maximum_jacobian_condition)
    census_rows = corner_isotropy ?
        seed_corner_census(semantic_corners, corner_grading, tolerance;
                           corner_kinds=corner_kinds, maximum_corner_aspect=maximum_corner_aspect,
                           corner_shape_gate=corner_shape_gate) :
        Dict{String, Any}[]
    corner_reach = corner_isotropy ?
        corner_law_reach(corner_isotropy_radius, lc_fine, lc_tangent, corner_grading_slope) :
        0.0
    # The ridge-to-ridge face census concerns the process edges' faces (metal
    # sidewalls), not the box faces the junction curves bound.
    face_rows = corner_isotropy ?
        longitudinal_face_census(setdiff(longitudinal_curves, junction_curves),
                                 semantic_corners, corner_reach) :
        Dict{String, Any}[]
    # SEAM (design round 2): the embedded tip bisector curves (node counts, the sheet
    # surfaces they are embedded in) and the pinched-seam census of the thin sheet, which
    # every thin coupon must pass with Count 0 (fail closed).
    isempty(tip_bisectors) || census_tip_bisector_curves!(tip_bisectors, tolerance)
    thin_sheet_seams = nothing
    if !fabricated && corner_isotropy
        sheet_surfaces = sort!(unique!(reduce(vcat, [surfaces for (attribute, surfaces) in boundary_groups
                                                     if 4000 <= attribute < 5000]; init=Int32[])))
        thin_sheet_seams = thin_sheet_seam_census(sheet_surfaces, lower, upper, tolerance)
        println("Thin sheet seam census: $(thin_sheet_seams["Count"]) pinched seams on " *
                "$(thin_sheet_seams["SheetFaces"]) sheet faces of $(length(sheet_surfaces)) surfaces; " *
                "tip bisectors $(length(tip_bisectors)), unrefined tips $(length(unrefined_tips))")
        for record in tip_bisectors
            println("  tip bisector at $(record["Point"]) ($(record["OpeningDegrees"]) degrees): length " *
                    "$(record["Length"]) ($(record["StationRule"])), $(record["Nodes"]) nodes, embedded in " *
                    "surfaces $(record["EmbeddedIn"])")
        end
        # Decision 368: the seams of a tip below the minimum opening stay (recorded for the
        # solve config); a seam away from every such tip is the fail-closed class.
        attributed = [edge for edge in thin_sheet_seams["Edges"]
                      if any(norm(0.5 .* (edge[1] .+ edge[2]) .- tip["Point"]) <= tip["TubeStation"] * (1.0 + 1.0e-6)
                             for tip in unrefined_tips)]
        unattributed = length(thin_sheet_seams["Edges"]) - length(attributed)
        thin_sheet_seams["UnrefinedCrackSeams"] = isempty(unrefined_tips) ? nothing : Dict{String, Any}(
            "Rule" => UNREFINED_CRACK_SEAMS_RULE, "Count" => length(attributed), "Edges" => attributed,
            "Tips" => unrefined_tips, "RefineCrackElements" => false)
        for tip in unrefined_tips
            println("  unrefined tip at $(tip["Point"]) ($(tip["OpeningDegrees"]) degrees < " *
                    "$(tip["MinimumOpeningDegrees"])): split floor $(tip["SplitFloor"]), seams within its tube " *
                    "station $(count(edge -> norm(0.5 .* (edge[1] .+ edge[2]) .- tip["Point"]) <= tip["TubeStation"] * (1.0 + 1.0e-6), attributed))")
        end
        unattributed == 0 ||
            error("ThinSheetSeams: $(unattributed) pinched seam edges on the thin sheet away from every " *
                  "unrefined tip (Palace's crack opener cannot split them on a prism / pyramid mesh): " *
                  "$([edge for edge in thin_sheet_seams["Edges"] if !(edge in attributed)][1:min(6, end)])")
    end
    gmsh.model.mesh.setOrder(mesh_order)
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    _, volume_element_tags, _ = gmsh.model.mesh.getElements(3)
    node_count = length(node_tags)
    element_count = sum(length(tags) for tags in volume_element_tags)
    println(
        "Spatial mesh budget: nodes=$node_count/$max_nodes, " *
        "volume_elements=$element_count/$max_elements"
    )
    node_count <= max_nodes ||
        error("Spatial coupon exceeds node budget: $node_count > $max_nodes")
    element_count <= max_elements ||
        error("Spatial coupon exceeds element budget: $element_count > $max_elements")
    element_count_after_generation = element_count
    if mesh_order > 1
        gmsh.model.mesh.optimize("HighOrderElastic", true, 20)
        gmsh.model.mesh.optimize("HighOrder", true, 20)
        _, element_tags, _ = gmsh.model.mesh.getElements(3)
        tags = reduce(vcat, element_tags; init=UInt64[])
        scaled_jacobians = gmsh.model.mesh.getElementQualities(tags, "minSJ")
        minimum_jacobian = minimum(scaled_jacobians)
        minimum_jacobian > 0.0 ||
            error("High-order spatial coupon contains a nonpositive Jacobian")
        println("High-order mesh minimum scaled Jacobian: $minimum_jacobian")
    end
    all_volume_tags=reduce(vcat,gmsh.model.mesh.getElements(3)[2];init=UInt64[])
    volume_qualities=gmsh.model.mesh.getElementQualities(all_volume_tags,"minSICN")
    minimum_signed_inverse_condition=minimum(volume_qualities)
    if !(minimum_signed_inverse_condition>1e-10)
        # Locate the degenerate cells before failing: their centroids tell which
        # construction (edge layer rows, corner ball, trace planes) produced them.
        degenerate=all_volume_tags[volume_qualities .<= 1e-10]
        for tag in degenerate[1:min(end, 8)]
            _,nodes,_,_=gmsh.model.mesh.getElement(tag)
            centroid=sum(gmsh.model.mesh.getNode(node)[1] for node in nodes) ./ length(nodes)
            println("Degenerate volume element $tag at $centroid")
        end
        error("Spatial coupon has invalid or near-singular elements: " *
              "minSICN=$minimum_signed_inverse_condition, count=$(length(degenerate))")
    end
    if mesh_postprocess !== nothing
        # Ownership contracts are planar source-local data.  Classify and certify
        # every interface element in that frame before applying the final rigid
        # transform; otherwise tilted transforms invalidate the XY/layer queries.
        mesh_postprocess(edges, boundary_loops, radius)
    end
    # The per-label areas are those of the labels the seed is written with: a
    # multi-slot postprocessor replaces the CAD-face groups by slot/conductor
    # groups, so the census is taken after it (areas are rigid-invariant, and the
    # source-local frame is kept by measuring before the placement transform).
    area_rows = corner_isotropy ? interface_areas() : Dict{String, Any}[]
    tube_record = prism_tubes ?
        prism_tube_census(tubes, tube_segments, tube_sections, tube_states, tube_volume_census,
                          tube_face_meshes, tube_curve_meshes, tube_timings, tube_band_record,
                          mixed_quality, graded_points, semantic_corners, junction_curves,
                          band_curves, footprint_polygons, (lower, upper, outer_tolerance),
                          lc_fine, lc_tangent, lc_far, far_growth,
                          trace_basis_record, max_elements, element_count_after_generation,
                          tube_layer_records) :
        nothing
    if transform != IDENTITY_RIGID_TRANSFORM
        connectivity = gmsh.model.mesh.getElements()
        gmsh.model.mesh.affineTransform(vec(transform'))
        connectivity == gmsh.model.mesh.getElements() ||
            error("Rigid transform changed mesh connectivity")
    end
    gmsh.write(filename)
    if corner_isotropy
        # Recorded seed-stage artifact, source-local frame, no pass/fail gate.
        ispath(corner_census) && error("Corner census output already exists")
        open(corner_census, "w") do stream
            write_json(stream, Dict{String, Any}(
                "Version" => 1, "Frame" => "SourceLocal",
                "Purpose" => "Seed corner-ball census, longitudinal-face census and interface areas; reported, not a qualification gate",
                "Scope" => prism_tubes ?
                    recipe_scope_record(scope_classes, boundary_loops, lower, upper, tolerance,
                                        tube_untubed_edges) :
                    nothing,
                "SemanticContract" => semantic_contract,
                "SemanticContractSHA256" => bytes2hex(sha256(read(semantic_contract))),
                "RigidTransform" => vec(transform'),
                "SemanticCorners" => [collect(corner) for corner in semantic_corners],
                "SemanticCornerKinds" => [kind === :invariant ? "Invariant" : "Legacy"
                                          for kind in corner_kinds],
                "InvariantCorners" => Dict{String, Any}(
                    "Rule" => INVARIANT_CORNER_RULE,
                    "Points" => [collect(corner) for (corner, kind) in
                                 zip(semantic_corners, corner_kinds) if kind === :invariant],
                    "CornerShapeGate" => corner_shape_gate > 0.0 ? corner_shape_gate : nothing,
                    "Target" => INVARIANT_CORNER_TARGET),
                "CornerIsotropyRadius" => corner_isotropy_radius,
                "IsotropicSize" => lc_fine,
                "Sqrt2IsotropicSize" => sqrt(2.0) * lc_fine,
                "FarSize" => lc_far,
                "CouponBox" => Dict{String, Any}(
                    "Rule" => COUPON_BOX_RULE, "Radius" => radius,
                    "Lower" => collect(lower), "Upper" => collect(upper),
                    "EdgeChains" => copy(EDGE_CHAIN_RECORDS),
                    "ChainedRows" => sum(Int[length(chain["Rows"]) for chain in EDGE_CHAIN_RECORDS]),
                    "ExtendedChains" => count(chain -> chain["UnionLength"] / 2 >= radius - 1.0e-10radius,
                                              EDGE_CHAIN_RECORDS)),
                "SizeBounds" => Dict{String, Any}(
                    "Rule" => SIZE_BOUND_RULE,
                    "FarSize" => lc_far,
                    "RequestedTangentialSize" => requested_lc_tangent,
                    "TangentialSize" => lc_tangent,
                    "TangentialSizeBoundByFarSize" => lc_tangent < requested_lc_tangent),
                "GradingTransitionWidth" => process_core_width - process_fine_width,
                "LongitudinalCurves" => length(longitudinal_curves),
                "CornerIsotropicLongitudinalCurves" => length(corner_curves),
                "CurveSpacing" => Dict{String, Any}(
                    "Rule" => CURVE_SPACING_RULE,
                    "GrowthRatio" => edge_growth_ratio,
                    "Count" => length(curve_spacing_records),
                    "GradedCurves" => count(record["Graded"] for record in curve_spacing_records),
                    "Curves" => curve_spacing_records),
                "JunctionCurves" => Dict{String, Any}(
                    "Count" => length(junction_curves),
                    "TotalLength" => junction_length,
                    # Source-local CAD endpoints of every junction curve (the Gmsh-only
                    # audit's junction band lines; a curved junction is recorded by its
                    # chord and flagged).
                    "Segments" => [vcat(a, b) for (a, b) in curve_segments(junction_curves)],
                    "CurvedCurves" => count(
                        abs(gmsh.model.occ.getMass(1, curve) - norm(b .- a)) >
                        1.0e-9 * max(norm(b .- a), 1.0)
                        for (curve, (a, b)) in zip(junction_curves, curve_segments(junction_curves))),
                    "Rule" => "curves of material-interface surfaces (substrate on one " *
                              "side, vacuum on the other) lying on the outer box: the " *
                              "cut-surface junctions of the trench floor/walls and the " *
                              "un-etched plane; feature curves with the process-edge band"),
                "CornerLawReach" => corner_reach,
                "CornerGrading" => Dict{String, Any}(
                    "CornerSize" => corner_size,
                    "GrowthRatio" => edge_growth_ratio,
                    "NormalSize" => lc_fine,
                    "Radius" => corner_isotropy_radius,
                    "Reach" => corner_grading_reach(corner_grading),
                    "ShellRadii" => corner_shell_radii(corner_grading),
                    "ShellSizes" => [corner_ball_size(corner_grading, inner)
                                     for inner in vcat(0.0, corner_shell_radii(corner_grading)[1:(end - 1)])],
                    "Rule" => corner_size > 0.0 ?
                        "inside every corner ball the isotropic size grows geometrically from " *
                        "CornerSize at the corner point in shells: shell k has size CornerSize x " *
                        "GrowthRatio^(k-1) and ends at the cumulative radius CornerSize " *
                        "(GrowthRatio^k - 1) / (GrowthRatio - 1) (ShellRadii, ShellSizes), " *
                        "NormalSize from the Reach to the ball radius; the Gmsh background " *
                        "field, the ridge node placement (nodes on the shell radii and on the " *
                        "ball boundary) and the seed optimizer's local bounds follow the shells, " *
                        "which never exceed the metric stage's continuous law min(NormalSize, " *
                        "CornerSize + (GrowthRatio - 1) x distance) (supervisor decision 33)" :
                        "uniform NormalSize inside every corner ball (no corner grading)"),
                "LongitudinalFaceHistogramBins" => LONGITUDINAL_FACE_HISTOGRAM_BINS,
                "LongitudinalFaces" => face_rows,
                "InterfaceAreaUnits" => "um^2",
                "EtchBoundary" => etch_boundary === nothing ? "producer-default" : etch_boundary,
                "EtchBoundarySHA256" => etch_boundary === nothing ? nothing :
                                        bytes2hex(sha256(read(etch_boundary))),
                "FootprintCollinearTolerance" => FOOTPRINT_COLLINEAR_TOLERANCE,
                "FootprintSimplification" => Dict{String, Any}(
                    "Rule" => "consecutive footprint edges are merged when every vertex " *
                              "between their outer endpoints lies within " *
                              "FootprintCollinearTolerance times the merged edge length " *
                              "of the merged edge, before CAD face creation",
                    "Polygons" => length(footprint_polygons),
                    "RemovedVertices" => sum(Int[polygon["Simplification"]["RemovedVertexCount"]
                                                 for polygon in footprint_polygons]),
                    "MaximumRelativeDeviation" => maximum(
                        Float64[polygon["Simplification"]["MaximumRelativeDeviation"]
                                for polygon in footprint_polygons]; init=0.0),
                    "ConstructionRule" => "an exterior loop's collar polygon is the miter " *
                                          "offset of its Physical sides (MiterOffset) while " *
                                          "that polygon is simple; a self-intersecting miter " *
                                          "polygon (facing sides closer than twice the collar) " *
                                          "is replaced by the outer boundary of the union of " *
                                          "the loop, the per-side collar rectangles and the " *
                                          "convex-corner miter kites, clipped to the coupon " *
                                          "box (CollarUnion); circular arcs of the loop are " *
                                          "offset exactly as concentric arcs (r + collar for a " *
                                          "convex-metal arc, r - collar for a concave one, " *
                                          "collapsed onto its neighbours' junction when " *
                                          "r <= collar), tangent joints staying continuous " *
                                          "(a junction turning by at most " *
                                          "$JUNCTION_TANGENT_ANGLE rad snaps its two shifted " *
                                          "ends, at most collar x that angle apart, to one " *
                                          "point: the one sub-nanometre tolerance of the exact " *
                                          "offsets), and the union " *
                                          "takes an annular sector per arc; where an offset " *
                                          "arc crosses another collar front on the union " *
                                          "boundary the crossing vertex is the chord-polyline " *
                                          "intersection (within the chord sagitta of the " *
                                          "exact arc) and the arc run is interrupted there",
                    "IslandRule" => "an un-etched region bounded entirely by collar " *
                                    "boundaries (not touching the box face) whose every " *
                                    "point lies within IslandExcessCap beyond the collar is " *
                                    "absorbed into the etched collar and recorded " *
                                    "(AbsorbedIslands of its polygon: area, maximum excess " *
                                    "found and its bound); the bound (half the island's " *
                                    "smallest width) measures beyond the constructed collar " *
                                    "boundary, kites included, the maximum excess found " *
                                    "measures the exact metal distance, and both must stay " *
                                    "within the cap; a larger island fails closed " *
                                    "(ScopeGuard FootprintTopology)",
                    "IslandExcessCapOverRadius" => COLLAR_ISLAND_EXCESS_CAP_OVER_RADIUS,
                    "IslandExcessCap" => COLLAR_ISLAND_EXCESS_CAP_OVER_RADIUS * radius,
                    "AbsorbedIslands" => sum(Int[length(get(polygon, "AbsorbedIslands", []))
                                                for polygon in footprint_polygons]),
                    "AbsorbedIslandArea" => sum(Float64[island["Area"]
                                                        for polygon in footprint_polygons
                                                        for island in get(polygon, "AbsorbedIslands", [])]),
                    "AbsorbedIslandMaximumExcess" => maximum(
                        Float64[island["MaximumExcess"] for polygon in footprint_polygons
                                for island in get(polygon, "AbsorbedIslands", [])]; init=0.0),
                    "AbsorbedIslandMaximumExcessBound" => maximum(
                        Float64[island["MaximumExcessBound"] for polygon in footprint_polygons
                                for island in get(polygon, "AbsorbedIslands", [])]; init=0.0),
                    "CollarUnionPolygons" => count(polygon["Construction"] == "CollarUnion"
                                                   for polygon in footprint_polygons)),
                "FootprintPolygons" => footprint_polygons,
                "PrismTubes" => tube_record,
                "EdgeLayer" => prism_tubes ? nothing : Dict{String, Any}(
                    "EdgeSize" => edge_size,
                    "GrowthRatio" => edge_growth_ratio,
                    "Aspect" => edge_layer_aspect,
                    "NormalSize" => lc_fine,
                    "TangentialSize" => lc_tangent,
                    "Layers" => length(edge_layer_offsets),
                    "RowOffsets" => edge_layer_offsets,
                    "LayerThickness" => isempty(edge_layer_offsets) ? 0.0 : edge_layer_offsets[end],
                    "TangentialSubdivision" => edge_layer_subdivision,
                    "RowSubdivisions" => edge_layer_row_subdivisions,
                    "RowTangentialSpacings" =>
                        [lc_tangent / n for n in edge_layer_row_subdivisions],
                    "TaperSubdivisions" => edge_layer_taper,
                    "CornerTaperOffset" => layer_to_ball ? corner_isotropy_radius :
                                           corner_reach + length(edge_layer_taper) * lc_tangent,
                    "LayerReachesCornerBall" => layer_to_ball,
                    "UnlayeredEdgeLengthPerCorner" => unlayered_edge_length_per_corner(
                        semantic_corners, edge_layer_curves, corner_isotropy_radius),
                    "RowZigzag" => EDGE_LAYER_ROW_ZIGZAG,
                    "Rule" => "on every face bounding a metal edge (longitudinal feature " *
                              "curve of a metal surface family, junction curves excluded) " *
                              "one embedded explicit node row per geometric layer of size " *
                              "EdgeSize x GrowthRatio^(k-1) below NormalSize, at the " *
                              "cumulative layer distance from the edge; inside the row span " *
                              "the ridge lc_tangent grid is subdivided by TangentialSubdivision " *
                              "(smallest power of two with spacing <= Aspect x EdgeSize) and " *
                              "row k by the nested power of two with spacing <= Aspect x its " *
                              "size, so every layer cell has aspect <= Aspect (a tetrahedron " *
                              "corner with three tangential edges has scaled Jacobian " *
                              "(hn/ht)^2); the ridge subdivision halves interval by interval " *
                              "(TaperSubdivisions) towards the span ends, where no row is " *
                              "seeded; rows start at the ridge grid outside the corner size " *
                              "law after the taper (CornerTaperOffset from a semantic corner) " *
                              "or, with corner grading (LayerReachesCornerBall), at the ridge " *
                              "node on the corner ball boundary with no taper, so the " *
                              "un-layered edge length per corner is the ball radius " *
                              "(UnlayeredEdgeLengthPerCorner: per semantic corner, the " *
                              "distances from the corner to the nearest span end of each " *
                              "adjacent layered edge and their maximum); " *
                              "every second row node is RowZigzag of the offset farther out (no " *
                              "Delaunay-degenerate rectangles); rows stop where the next " *
                              "layer leaves the face",
                    "Curves" => edge_layer_curves,
                    "TotalSpanLength" =>
                        sum(Float64[row["SpanLength"] for row in edge_layer_curves]),
                    "Rows" =>
                        sum(Int[length(record["Rows"]) for (record, _) in edge_layer_records]),
                    "RowNodes" =>
                        sum(Int[record["RowNodes"] for (record, _) in edge_layer_records]),
                    "RidgeNodesAdded" => edge_layer_ridge_nodes),
                "TraceBasisSizing" => trace_basis_record,
                "SeedQualityOptimization" => seed_quality,
                "TipBisectors" => Dict{String, Any}(
                    "Rule" => TIP_BISECTOR_RULE, "Count" => length(tip_bisectors),
                    "Curves" => tip_bisectors,
                    "MinimumOpeningDegrees" => corner_shape_gate > 0.0 ?
                                               rad2deg(tip_bisector_minimum_opening(corner_shape_gate)) :
                                               nothing,
                    "DescentExcess" => TIP_BISECTOR_DESCENT_EXCESS,
                    "SplitFloorRule" => "the 2D condition number of the affine map from the " *
                                        "equilateral triangle onto the isoceles triangle of apex " *
                                        "angle phi / 2 (singular-value interlacing: the kappa_reg " *
                                        "floor of every cell on a split fan sector); the bisector is " *
                                        "embedded for phi >= 2 alpha* with floor(alpha*) = " *
                                        "CornerShapeGate / DescentExcess (supervisor decision 368)",
                    "UnrefinedTips" => length(unrefined_tips)),
                "ThinSheetSeams" => thin_sheet_seams,
                "InterfaceAreas" => area_rows,
                "Corners" => census_rows))
            println(stream)
        end
        println("Seed corner census: $corner_census")
        for row in census_rows
            println("  corner $(row["Corner"]) $(row["Point"]): ball edges=$(row["BallEdges"]) " *
                    "min/median/max=$(row["EdgeMinimum"])/$(row["EdgeMedian"])/$(row["EdgeMaximum"]) " *
                    "over sqrt2*hn=$(row["EdgesOverSqrt2IsotropicSize"]) " *
                    "incident aspect=$(row["IncidentMaximumAspect"])")
        end
        for row in area_rows
            println("  interface $(row["Attribute"]) $(row["Name"]): triangles=$(row["Triangles"]) " *
                    "area=$(row["Area"]) um^2")
        end
        for polygon in footprint_polygons
            record = polygon["Simplification"]
            println("  footprint polygon conductor $(polygon["Conductor"]) plane $(polygon["Plane"]) " *
                    "hole=$(polygon["Hole"]): vertices $(record["OriginalVertices"]) -> " *
                    "$(record["Vertices"]), removed $(record["RemovedVertexIndices"]), " *
                    "max deviation $(record["MaximumDeviation"]) " *
                    "(relative $(record["MaximumRelativeDeviation"]))")
        end
        for row in face_rows
            println("  longitudinal face $(row["Surface"]) $(row["PhysicalGroups"]): " *
                    "triangles=$(row["Triangles"]) interior nodes=$(row["InteriorNodes"]) " *
                    "(away from corners $(row["InteriorNodesAwayFromCorners"])) " *
                    "full-height triangles=$(row["FullHeightTriangles"]) " *
                    "(away from corners $(row["FullHeightTrianglesAwayFromCorners"]))")
        end
    end
    metadata_path = filename * ".metadata.json"
    open(metadata_path, "w") do stream
        println(stream, "{")
        println(stream, "  \"Version\": 1,")
        println(stream, "  \"MetalSurfacePartition\": \"InterfaceSlotAndConductor\",")
        println(stream, "  \"NodeCount\": $node_count,")
        println(stream, "  \"VolumeElementCount\": $element_count,")
        println(stream, "  \"MinimumSignedInverseCondition\": $minimum_signed_inverse_condition,")
        println(stream, "  \"FineSize\": $lc_fine,")
        println(stream, "  \"TangentialSize\": $lc_tangent,")
        println(stream, "  \"FarSize\": $lc_far,")
        println(stream, "  \"ProcessCoreWidth\": $process_core_width,")
        println(stream, "  \"ProcessFineWidth\": $process_fine_width,")
        println(stream, "  \"ProcessGradingPower\": $process_grading_power,")
        println(stream, "  \"CornerIsotropyRadius\": $corner_isotropy_radius,")
        println(stream, "  \"CornerIsotropicSize\": $(corner_isotropy ? lc_fine : 0.0),")
        println(stream, "  \"EdgeSize\": $edge_size,")
        println(stream, "  \"EdgeGrowthRatio\": $edge_growth_ratio,")
        println(stream, "  \"EdgeLayerAspect\": $edge_layer_aspect,")
        println(stream, "  \"EdgeLayerMaximumAspect\": $edge_layer_maximum_aspect,")
        println(stream, "  \"EdgeLayerRows\": $(length(edge_layer_offsets)),")
        println(stream, "  \"PrismTubes\": $(prism_tubes ? "true" : "false"),")
        println(stream, "  \"FarGrowth\": $far_growth,")
        println(stream, "  \"SemanticCornerCount\": $(length(semantic_corners)),")
        println(stream, "  \"RigidTransform\": [$(join(vec(transform'), ", "))],")
        println(stream, "  \"Algorithm3D\": $(Int(gmsh.option.getNumber("Mesh.Algorithm3D"))),")
        println(stream, "  \"MeshOrder\": $mesh_order")
        println(stream, "}")
    end
    println("Spatial mesh metadata: $metadata_path")
    println(
        "Spatial coupon: fabricated=$fabricated, edges=$(length(edges)), " *
        "layers=$(length(layers)), file=$filename"
    )
    return gmsh.finalize()
end

function parse_options(args)
    length(args) >= 3 || error(
        "Usage: mesh_spatial_coupon.jl SIGNATURE.csv thin|fabricated " *
        "OUTPUT.msh [options]"
    )
    args[2] in ("thin", "fabricated") || error("Coupon kind must be thin or fabricated")
    options = Dict{String, Any}(
        "signature" => abspath(args[1]),
        "fabricated" => args[2] == "fabricated",
        "filename" => abspath(args[3])
    )
    names = Dict(
        "--mask" => ("mask", String),
        "--boundary" => ("boundary", String),
        "--radius" => ("radius", Float64),
        "--metal-thickness" => ("metal_thickness", Float64),
        "--overetch" => ("overetch", Float64),
        "--sidewall-angle" => ("sidewall_angle", Float64),
        "--top-radius" => ("top_rounding", Float64),
        "--bottom-radius" => ("trench_rounding", Float64),
        "--lc-fine" => ("lc_fine", Float64),
        "--lc-tangent" => ("lc_tangent", Float64),
        "--lc-far" => ("lc_far", Float64),
        "--process-core-width" => ("process_core_width", Float64),
        "--process-fine-width" => ("process_fine_width", Float64),
        "--process-grading-power" => ("process_grading_power", Float64),
        "--max-nodes" => ("max_nodes", Int),
        "--max-elements" => ("max_elements", Int),
        "--mesh-order" => ("mesh_order", Int),
        "--rigid-transform" => ("transform", Matrix{Float64}),
        "--interface-ownership-report" => ("interface_ownership_report", String),
        "--semantic-contract" => ("semantic_contract", String),
        "--corner-isotropy-radius" => ("corner_isotropy_radius", Float64),
        "--corner-census" => ("corner_census", String),
        "--labels-only" => ("labels_only", String),
        "--etch-boundary" => ("etch_boundary", String),
        "--trace-basis-contract" => ("trace_basis_contract", String),
        "--trace-vertices" => ("trace_vertices", String),
        "--trace-triangles" => ("trace_triangles", String),
        "--process-library" => ("process_library", String),
        "--trace-basis-size-ratio" => ("trace_basis_size_ratio", Float64),
        "--edge-size" => ("edge_size", Float64),
        "--edge-growth-ratio" => ("edge_growth_ratio", Float64),
        "--edge-layer-aspect" => ("edge_layer_aspect", Float64),
        "--maximum-corner-aspect" => ("maximum_corner_aspect", Float64),
        "--minimum-scaled-jacobian" => ("minimum_scaled_jacobian", Float64),
        "--maximum-jacobian-condition" => ("maximum_jacobian_condition", Float64),
        "--maximum-quality-displacement-over-normal" =>
            ("quality_displacement_over_normal", Float64),
        "--edge-layer-maximum-aspect" => ("edge_layer_maximum_aspect", Float64),
        "--corner-size" => ("corner_size", Float64),
        "--prism-tubes" => ("prism_tubes", Bool),
        "--tube-sector-degrees" => ("tube_sector_degrees", Float64),
        "--far-growth" => ("far_growth", Float64),
        "--corner-shape-gate" => ("corner_shape_gate", Float64)
    )
    index = 4
    while index <= length(args)
        flag = args[index]
        haskey(names, flag) || error("Unknown option: $flag")
        index < length(args) || error("Missing value for option: $flag")
        name, type = names[flag]
        options[name] = if type === String
            abspath(args[index + 1])
        elseif type === Matrix{Float64}
            parse_rigid_transform(args[index + 1])
        else
            parse(type, args[index + 1])
        end
        index += 2
    end
    return options
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_options(ARGS)
    postprocess = nothing
    if haskey(options, "interface_ownership_report")
        include(joinpath(@__DIR__, "label_interface_patches.jl"))
        report = options["interface_ownership_report"]
        fabricated = options["fabricated"]
        thickness = get(options, "metal_thickness", 0.1)
        overetch = get(options, "overetch", 0.05)
        postprocess = (edges, loops, radius) -> label_interface_patches(
            edges, loops, radius, report; fabricated=fabricated,
            metal_thickness=thickness, overetch=overetch)
    end
    generate_spatial_coupon(;
        signature       = options["signature"],
        mask            = get(options, "mask", nothing),
        boundary        = get(options, "boundary", nothing),
        fabricated      = options["fabricated"],
        filename        = options["filename"],
        radius          = get(options, "radius", 2.0),
        metal_thickness = get(options, "metal_thickness", 0.1),
        overetch        = get(options, "overetch", 0.05),
        sidewall_angle  = get(options, "sidewall_angle", 80.0),
        top_rounding    = get(options, "top_rounding", 0.01),
        trench_rounding = get(options, "trench_rounding", 0.01),
        lc_fine         = get(options, "lc_fine", 0.02),
        lc_tangent      = get(options, "lc_tangent", 0.0),
        lc_far          = get(options, "lc_far", 0.3),
        process_core_width = get(options, "process_core_width", 0.0),
        process_fine_width = get(options, "process_fine_width", 0.0),
        process_grading_power = get(options, "process_grading_power", 1.7),
        max_nodes       = get(options, "max_nodes", 500_000),
        max_elements    = get(options, "max_elements", 2_000_000),
        mesh_order      = get(options, "mesh_order", 1),
        mesh_postprocess = postprocess,
        transform       = get(options, "transform", copy(IDENTITY_RIGID_TRANSFORM)),
        semantic_contract = get(options, "semantic_contract", nothing),
        corner_isotropy_radius = get(options, "corner_isotropy_radius", 0.0),
        corner_census   = get(options, "corner_census", nothing),
        labels_only     = get(options, "labels_only", nothing),
        etch_boundary   = get(options, "etch_boundary", nothing),
        trace_basis_contract = get(options, "trace_basis_contract", nothing),
        trace_vertices  = get(options, "trace_vertices", nothing),
        trace_triangles = get(options, "trace_triangles", nothing),
        process_library = get(options, "process_library", nothing),
        trace_basis_size_ratio = get(options, "trace_basis_size_ratio", 1.0),
        edge_size = get(options, "edge_size", 0.0),
        edge_growth_ratio = get(options, "edge_growth_ratio", 2.0),
        edge_layer_aspect = get(options, "edge_layer_aspect", 4.0),
        maximum_corner_aspect = get(options, "maximum_corner_aspect", 0.0),
        minimum_scaled_jacobian = get(options, "minimum_scaled_jacobian", 0.0),
        maximum_jacobian_condition = get(options, "maximum_jacobian_condition", 0.0),
        quality_displacement_over_normal =
            get(options, "quality_displacement_over_normal", 0.0),
        edge_layer_maximum_aspect = get(options, "edge_layer_maximum_aspect", 0.0),
        corner_size = get(options, "corner_size", 0.0),
        prism_tubes = get(options, "prism_tubes", false),
        tube_sector_degrees = get(options, "tube_sector_degrees", 30.0),
        far_growth = get(options, "far_growth", 0.0),
        corner_shape_gate = get(options, "corner_shape_gate", 0.0)
    )
end
