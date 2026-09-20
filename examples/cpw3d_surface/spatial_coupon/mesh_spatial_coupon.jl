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
                vertex_arm=Bool(signature_integer(data[row, columns["VertexArm"]], "VertexArm"))
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
        push!(
            loops,
            (conductor=conductor, plane=plane, hole=hole, points=points, classes=classes)
        )
    end
    isempty(loops) && error("Plan-view boundary contains no loops")
    return loops
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

function seed_corner_census(corners, grading::CornerGrading, tolerance)
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
        push!(rows, Dict{String, Any}(
            "Corner" => k - 1, "Point" => collect(corner),
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

# Movement basis of every surface vertex: null space of its triangle normals.
function surface_movement_bases(points, triangles)
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
    aspect_of = objective === :edge_aspect ? tetrahedron_edge_aspect : tetrahedron_aspect
    function value()
        if objective === :aspect || objective === :edge_aspect
            aspects = [aspect_of(cell_xyz(tetrahedra[k])) for k in targets]
            scale = maximum(aspects)
            return scale * sum((aspects ./ scale) .^ SEED_ASPECT_PROXY_POWER)^(1 / SEED_ASPECT_PROXY_POWER)
        end
        return -minimum(tetrahedron_scaled_jacobian(cell_xyz(tetrahedra[k])) for k in targets)
    end
    function achieved()
        if objective === :aspect || objective === :edge_aspect
            return maximum(aspect_of(cell_xyz(tetrahedra[k])) for k in targets)
        end
        return -value()
    end
    minimizing = objective === :aspect || objective === :edge_aspect
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
                                   fixed=falses(size(points, 2)))
    original = copy(points)
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
    bases = surface_movement_bases(points, triangles)
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
    corner_before = Float64[]; corner_after = Float64[]; corner_moves = Int[]
    for corner in corners
        center = collect(corner)
        targets = [k for (k, cell) in enumerate(tetrahedra)
                   if any(norm(points[:, i] .- center) <= tolerance for i in cell)]
        isempty(targets) && error("Semantic corner is absent from the seed volume mesh")
        before = maximum(tetrahedron_aspect([points[:, i] for i in tetrahedra[k]]) for k in targets)
        push!(corner_before, before)
        if before <= corner_target
            push!(corner_after, before); push!(corner_moves, 0)
            continue
        end
        after, moves = optimize_seed_cells!(points, original, tetrahedra, incident, bases, floors,
                                            bounds, targets, :aspect, corner_target;
                                            ceilings=ceilings, condition_ceilings=condition_ceilings,
                                            edge_floors=edge_floors)
        push!(corner_after, after); push!(corner_moves, moves)
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
    if !all(<=(maximum_corner_aspect), corner_after)
        for (corner, after) in zip(corners, corner_after)
            after <= maximum_corner_aspect && continue
            center = collect(corner)
            incident = [k for (k, cell) in enumerate(tetrahedra)
                        if any(norm(points[:, i] .- center) <= tolerance for i in cell)]
            aspect_of(k) = tetrahedron_aspect([points[:, i] for i in tetrahedra[k]])
            for k in sort(incident; by=aspect_of, rev=true)[1:min(6, end)]
                describe("corner $(center) aspect", k, "aspect $(aspect_of(k))")
            end
        end
        push!(failures, "Seed semantic-corner aspect exceeds the gate after optimization: " *
                        "$(maximum(corner_after)) > $(maximum_corner_aspect)")
    end
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
    return record, moved, collapse
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
                                        fixed_node_tags=UInt[])
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
    record, moved, collapse = optimize_required_region!(
        points, tetrahedra, triangles, corners, radius, lc_fine, spans, edge_size, growth_ratio,
        layer_thickness, row_zigzag, maximum_corner_aspect, minimum_scaled_jacobian,
        maximum_jacobian_condition, displacement_ratio, tolerance;
        edge_layer_maximum_aspect=edge_layer_maximum_aspect, corner_grading=corner_grading,
        triangle_entities=triangle_entities, lines=lines, line_entities=line_entities,
        fixed=fixed)
    for i in moved
        gmsh.model.mesh.setNode(node_tags[i], points[:, i], Float64[])
    end
    apply_seed_cell_collapse!(node_tags, tags, entities, collapse.cells_removed,
                              collapse.cells_remapped)
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

function point_in_mask(facets, point, conductor, plane, tolerance)
    return any(
        facet.conductor == conductor &&
        abs(facet.plane - plane) <= tolerance &&
        point_in_polygon((point[1], point[2]), facet.points, tolerance) for facet in facets
    )
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

function polygon_wire(occ, points, z)
    tolerance =
        1.0e-9 * max(
            maximum(point[1] for point in points) - minimum(point[1] for point in points),
            maximum(point[2] for point in points) - minimum(point[2] for point in points),
            1.0
        )
    runs = circular_arc_runs(points, tolerance)
    projected = collect(points)
    for run in runs, index in run.point_indices[2:(end - 1)]
        radial = (points[index][1] - run.center[1], points[index][2] - run.center[2])
        scale = run.radius / hypot(radial...)
        projected[index] =
            (run.center[1] + scale * radial[1], run.center[2] + scale * radial[2])
    end
    tags = [occ.addPoint(point[1], point[2], z) for point in projected]
    curve_for_edge = Dict{Int, Int32}()
    covered = falses(length(points))
    for run in runs
        any(covered[run.edge_indices]) && continue
        parts = max(1, ceil(Int, run.angle / (0.5 * pi)))
        split = unique(round.(Int, range(1, length(run.point_indices), length=parts + 1)))
        center = occ.addPoint(run.center[1], run.center[2], z)
        for (first, second) in zip(split, split[2:end])
            curve_for_edge[run.edge_indices[first]] = occ.addCircleArc(
                tags[run.point_indices[first]],
                center,
                tags[run.point_indices[second]]
            )
        end
        covered[run.edge_indices] .= true
    end
    curves = Int32[]
    for index in eachindex(tags)
        if haskey(curve_for_edge, index)
            push!(curves, curve_for_edge[index])
        elseif !covered[index]
            push!(curves, occ.addLine(tags[index], tags[mod1(index + 1, length(tags))]))
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

function offset_loop_points(loop, distance, tolerance)
    abs(distance)<=tolerance && return loop.points
    isempty(circular_arc_runs(loop.points,tolerance)) ||
        error("Nonzero offsets of curved plan-view boundaries require exact curved-offset support")
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

function footprint_record(conductor, plane, hole, points, record)
    return Dict{String, Any}(
        "Conductor" => conductor, "Plane" => plane, "Hole" => hole,
        "Points" => [collect(point) for point in points], "Simplification" => record)
end

# Simplify the bottom and top polygons of a footprint loft; when `footprint` is a
# vector, the loft must be prismatic (one polygon) and the polygon is recorded.
function simplified_loft_polygons(bottom_points, top_points, footprint, conductor, plane, hole)
    bottom_points, bottom_record =
        simplify_footprint_polygon(bottom_points, FOOTPRINT_COLLINEAR_TOLERANCE)
    top_points, _ = simplify_footprint_polygon(top_points, FOOTPRINT_COLLINEAR_TOLERANCE)
    if footprint !== nothing
        bottom_points == top_points ||
            error("Footprint recording requires a prismatic (vertical-wall) loft")
        push!(footprint, footprint_record(conductor, plane, hole, bottom_points, bottom_record))
    end
    return bottom_points, top_points
end

function loft_polygon(occ, bottom_points, top_points, z0, z1)
    bottom = polygon_wire(occ, bottom_points, z0)
    top = polygon_wire(occ, top_points, z1)
    entities = occ.addThruSections([bottom, top], -1, true, false, -1, "C0")
    volumes = [(dim, tag) for (dim, tag) in entities if dim == 3]
    isempty(volumes) && error("Plan-view mask loft produced no volume")
    return volumes
end

function offset_hole_points(loop,distance,tolerance)
    distance>=-tolerance && return offset_loop_points(loop,distance,tolerance)
    isempty(circular_arc_runs(loop.points,tolerance)) ||
        error("Shrinking a curved fabrication hole requires exact curved-offset support")
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
    clipped=copy(points)
    for i in eachindex(points)
        loop.classes[i]=="Physical" || continue
        a,b=points[i],points[mod1(i+1,length(points))]
        direction=(b[1]-a[1],b[2]-a[2]);edge_length=hypot(direction...)
        normal=(-orientation*direction[2]/edge_length,orientation*direction[1]/edge_length)
        signed(p)=normal[1]*(p[1]-a[1])+normal[2]*(p[2]-a[2])+distance
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
# genuine facet); `footprint` additionally records every simplified polygon.
function loft_mask_offsets(occ, loops, z0, z1, bottom_offset, top_offset, tolerance;
                           simplify=false, footprint=nothing)
    footprint === nothing || simplify || error("Footprint recording requires simplification")
    outers = [loop for loop in loops if !loop.hole]
    holes = [loop for loop in loops if loop.hole]
    isempty(outers) && error("Plan-view mask has no exterior loop")
    result = Tuple{Int32, Int32}[]
    hole_owners=zeros(Int,length(holes))
    for outer in outers
        bottom_points = offset_loop_points(outer, bottom_offset, tolerance)
        top_points = offset_loop_points(outer, top_offset, tolerance)
        if simplify
            bottom_points, top_points = simplified_loft_polygons(
                bottom_points, top_points, footprint, outer.conductor, z0, false)
        end
        volume = loft_polygon(occ, bottom_points, top_points, z0, z1)
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
                    bottom_hole, top_hole, footprint, hole.conductor, z0, true)
            end
            append!(
                cutters,
                loft_polygon(occ,bottom_hole,top_hole,z0,z1)
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
function boundary_strips(occ, loops, radius, z0, z1, pullback, tolerance; footprint=nothing)
    expanded_volumes = Tuple{Int32, Int32}[]
    retained_volumes = Tuple{Int32, Int32}[]
    width = 3radius
    for conductor in sort!(unique(loop.conductor for loop in loops))
        conductor_loops = [loop for loop in loops if loop.conductor == conductor]
        append!(expanded_volumes,
                loft_mask_offsets(occ, conductor_loops, z0, z1, -width, -width, tolerance;
                                  simplify=true, footprint=footprint))
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
        false)
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
        wires = [polygon_wire(occ, outer.points, outer.plane)]
        for (i,hole) in enumerate(holes)
            hole.conductor == outer.conductor &&
                abs(hole.plane - outer.plane) <= tolerance || continue
            point_in_polygon(hole.points[1], outer.points, tolerance) || continue
            owners[i]+=1
            push!(wires, polygon_wire(occ, hole.points, hole.plane))
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

function coupon_bounds(edges, radius, metal_thickness, overetch)
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
    gmsh.model.isInside(2, tag, coordinate) > 0 && return center

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

function nearest_metal_edge(edges, facets, point, radius, tolerance)
    candidates = if isempty(facets)
        edges
    else
        [
            edge for edge in edges if
            point_in_mask(facets, point, edge.conductor, edge.point[3], tolerance)
        ]
    end
    isempty(candidates) && error("Unable to assign a metal surface to a conductor mask")
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
# input exhibits from the frozen inputs alone. The census Scope block records
# them and mesh_stage_contract.py binds the same lists.
const RECIPE_SCOPE_RECIPE = "prism-tubes"
const RECIPE_SCOPE_SUPPORTED_CLASSES = [
    "ContinuationVertices", "DeviceFootprint", "DownwardLayers", "ExteriorLoops", "HoleLoops",
    "MultipleConductors", "MultipleLayers", "MultipleSlots", "TraceBasis"]
const RECIPE_SCOPE_GUARDS = [
    ("TopRounding", "inputs",
     "rounded metal top edges (TopRounding > 0): the tube rings surround a sharp edge"),
    ("TrenchRounding", "inputs",
     "rounded trench edges (TrenchRounding > 0): the bottom tube splits its substrate and " *
     "vacuum sectors at a sharp trench wall"),
    ("SlopedSidewalls", "inputs",
     "sidewall angle below 90 degrees: the tube sections assume vertical metal faces"),
    ("ThinMetal", "inputs",
     "thin (non-fabricated) metal: no sidewall and no top / bottom edge pair to tube"),
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
    ("FreeEdgeEnds", "build",
     "a metal edge end that is neither a semantic corner nor on the outer box"),
    ("ShortEdges", "build",
     "a metal edge shorter than the corner clearances at its ends: no tube interval remains"),
    ("FootprintWithoutEdge", "build",
     "an explicit etch footprint with no side coincident with a metal edge")]
const RECIPE_SCOPE_RULE =
    "the prism-tube recipe builds every input whose classes are all in SupportedClasses; " *
    "an input exhibiting a class in GuardedClasses fails closed at the guard whose " *
    "error message carries ScopeGuard[<Id>] (Guards); ExhibitedClasses are the classes " *
    "of this input among both lists, from the frozen inputs (loops, layers, process); " *
    "MetalLoops counts per plan-view loop the straight sides not on the outer box, and " *
    "TubeCount = 2 x their sum (a top and a bottom tube per side of every loop)"

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
    return sort!(classes)
end

# A plan-view side lies on the outer box when both ends are on the same box face.
function side_on_box_face(p, q, lower, upper, tolerance)
    return any((abs(p[d] - lower[d]) <= tolerance && abs(q[d] - lower[d]) <= tolerance) ||
               (abs(p[d] - upper[d]) <= tolerance && abs(q[d] - upper[d]) <= tolerance)
               for d in 1:2)
end

# Per plan-view loop: the straight sides not on the outer box (each carries a top
# and a bottom tube).
function metal_loop_records(loops, lower, upper, tolerance)
    records = Dict{String, Any}[]
    for (index, loop) in enumerate(loops)
        n = length(loop.points)
        sides = count(!side_on_box_face(loop.points[i], loop.points[i % n + 1], lower, upper,
                                        tolerance) for i in 1:n)
        push!(records, Dict{String, Any}(
            "Loop" => index, "Conductor" => loop.conductor, "Plane" => loop.plane,
            "Hole" => loop.hole, "Vertices" => n, "Sides" => sides))
    end
    return records
end

function recipe_scope_record(exhibited, loops, lower, upper, tolerance)
    return Dict{String, Any}(
        "Rule" => RECIPE_SCOPE_RULE, "Recipe" => RECIPE_SCOPE_RECIPE,
        "SupportedClasses" => copy(RECIPE_SCOPE_SUPPORTED_CLASSES),
        "GuardedClasses" => [guard[1] for guard in RECIPE_SCOPE_GUARDS],
        "Guards" => [Dict{String, Any}("Id" => id, "DetectedFrom" => detected,
                                       "Statement" => statement)
                     for (id, detected, statement) in RECIPE_SCOPE_GUARDS],
        "ExhibitedClasses" => copy(exhibited),
        "MetalLoops" => metal_loop_records(loops, lower, upper, tolerance))
end

# Straight metal edges of the plan-view boundary loops: the Physical sides (the
# sides not lying on the outer box), with the horizontal normal pointing away
# from the metal and the tube interval shrunk by the corner clearance at semantic
# corners (0 at box continuation vertices). Returns rows with start, stop,
# direction, normal, span, s_start, s_end, plane, conductor, corner angles.
function metal_edge_segments(loops, corners, clearance_of_angle, lower, upper, tolerance)
    segments = NamedTuple[]
    on_box(p) = any(abs(p[d] - lower[d]) <= tolerance || abs(p[d] - upper[d]) <= tolerance
                    for d in 1:2)
    is_corner(p, plane) = any(norm(collect(corner) .- [p[1], p[2], plane]) <= tolerance
                              for corner in corners)
    # Every straight side first, so that the angle between the sides meeting at a
    # corner is known before the clearance is applied.
    sides = NamedTuple[]
    for loop in loops
        # The metal lies inside an exterior loop and outside a hole (loop.hole: the
        # polygon interior is dielectric), so the tube normal, which points away from
        # the metal, points out of an exterior loop and into a hole.
        metal_inside = !loop.hole
        n = length(loop.points)
        for i in 1:n
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
                    scope_error("FreeEdgeEnds", "metal edge end $point of $p -> $q")
            end
            push!(sides, (start=[p[1], p[2]], stop=[q[1], q[2]], direction=direction,
                          normal=normal, span=span, plane=loop.plane, conductor=loop.conductor,
                          hole=loop.hole,
                          start_corner=is_corner(p, loop.plane), stop_corner=is_corner(q, loop.plane)))
        end
    end
    isempty(sides) && error("No straight metal edges found")
    # In-plane angle between two tube edges meeting at a corner: the smallest
    # angle between the directions pointing away from the corner (pi when the
    # corner has a single tube edge).
    function corner_angle(point, plane, own_direction_away)
        angle = Float64(pi)
        for side in sides
            abs(side.plane - plane) <= tolerance || continue
            for (end_point, away) in ((side.start, side.direction), (side.stop, -side.direction))
                norm(end_point .- point) <= tolerance || continue
                dot(away, own_direction_away) >= 1.0 - 1.0e-12 && continue
                angle = min(angle, acos(clamp(dot(away, own_direction_away), -1.0, 1.0)))
            end
        end
        return angle
    end
    for side in sides
        start_angle = side.start_corner ? corner_angle(side.start, side.plane, side.direction) : Float64(pi)
        stop_angle = side.stop_corner ? corner_angle(side.stop, side.plane, -side.direction) : Float64(pi)
        s_start = side.start_corner ? clearance_of_angle(start_angle) : 0.0
        s_end = side.span - (side.stop_corner ? clearance_of_angle(stop_angle) : 0.0)
        s_end > s_start ||
            scope_error("ShortEdges", "metal edge $(side.start) -> $(side.stop) of span " *
                                      "$(side.span) against clearances $(s_start) and " *
                                      "$(side.span - s_end)")
        push!(segments, (side..., s_start=s_start, s_end=s_end,
                         corner_angles=(start_angle, stop_angle)))
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
    for loop in etch_loops
        n = length(loop.points)
        for i in 1:n
            p = loop.points[i]
            q = loop.points[i % n + 1]
            for (a, b) in ((p, q), (q, p))
                norm([a[1], a[2]] .- segment.start) <= tolerance &&
                    norm([b[1], b[2]] .- segment.stop) <= tolerance && return true
            end
        end
    end
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
                           lc_fine, lower, upper, tolerance)
    overetch > 0.0 || scope_error("NoTrench", "Overetch $overetch")
    sectors = round(Int, 270.0 / sector_degrees)
    abs(sectors * sector_degrees - 270.0) <= 1.0e-9 || error("Tube sector angle must divide 270 degrees")
    per_quadrant = round(Int, 90.0 / sector_degrees)
    abs(per_quadrant * sector_degrees - 90.0) <= 1.0e-9 || error("Tube sector angle must divide 90 degrees")
    bound = min(overetch, 0.5 * metal_thickness, corner_radius)
    rings = tube_ring_count(edge_size, ratio, bound)
    # Top edge: vacuum from the sidewall ray (-90) through the outward normal (0)
    # and up (90) to the top-face ray (180). Bottom edge: substrate from the
    # metal bottom face ray (180) to the trench wall ray (270), vacuum from the
    # trench wall past the outward normal (360) to the sidewall ray (450).
    top_section = TubeSection(edge_size, ratio, rings,
                              [-90.0 + sector_degrees * j for j in 0:sectors], fill(2, sectors))
    bottom_section = TubeSection(edge_size, ratio, rings,
                                 [180.0 + sector_degrees * j for j in 0:sectors],
                                 vcat(fill(1, per_quadrant), fill(2, sectors - per_quadrant)))
    radius = tube_radius(top_section)
    outer_ring = ring_sizes(top_section)[end]
    pyramid_height = TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING * outer_ring
    radius + pyramid_height < overetch ||
        scope_error("ShallowTrench", "tube radius $radius + pyramid height $pyramid_height " *
                                     "against Overetch $overetch")
    clearance(angle) = (angle < pi - 1.0e-9 ? radius / tan(0.5 * angle) : 0.0) + outer_ring
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
    tubes = Tuple{EdgeTube, TubeSection}[]
    segments = NamedTuple[]
    # Facing process layers (an upward layer below a downward one, the only pair
    # layer_groups admits): the vacuum gap between their metal top faces must exceed
    # two tube reaches, like the width of a hole.
    for lower_layer in layers, upper_layer in layers
        lower_layer.sign > 0 && upper_layer.sign < 0 || continue
        gap = (upper_layer.plane - metal_thickness) - (lower_layer.plane + metal_thickness)
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
                          metal_edge_segments(layer_loops, corners, clearance, lower, upper, tolerance)]
        if etch_loops !== nothing
            for segment in layer_segments
                assert_etch_carries_edge(etch_loops, segment, tolerance)
            end
        end
        for segment in layer_segments
            # The tube frame follows the layer: b is the process normal (Nz), so the
            # sections' "up" (towards the metal top face) and the extrusion sense
            # e = n x b mirror for a downward layer; the top tube sits on the metal top
            # face at plane + Nz x thickness, the bottom tube on the plane.
            n = [segment.normal[1], segment.normal[2], 0.0]
            b = [0.0, 0.0, Float64(layer.sign)]
            e = cross(n, b)
            along = dot(e[1:2], segment.direction)
            abs(abs(along) - 1.0) <= 1.0e-12 || error("Tube frame is not aligned with the edge")
            for (z, section) in ((layer.plane + layer.sign * metal_thickness, top_section),
                                 (layer.plane, bottom_section))
                tube = if along > 0.0
                    EdgeTube([segment.start[1], segment.start[2], z], n, b, segment.s_start,
                             segment.s_end, lc_tangent)
                else
                    EdgeTube([segment.stop[1], segment.stop[2], z], n, b,
                             segment.span - segment.s_end, segment.span - segment.s_start, lc_tangent)
                end
                push!(tubes, (tube, section))
                for (tool, group) in add_tube_volumes!(occ, tube, section)
                    push!(records, TubeRecord(tube, section, group, tool))
                    push!(tools, tool)
                end
            end
        end
        append!(segments, layer_segments)
    end
    description = Dict{String, Any}(
        "Rings" => rings, "RingSizes" => ring_sizes(top_section),
        "RingRadii" => copy(top_section.ring_radii), "Radius" => radius,
        "TransverseBound" => bound,
        "RingRule" => "largest K with r_K + h_K <= min(Overetch, MetalThickness / 2, " *
                      "CornerIsotropyRadius)",
        "SectorDegrees" => sector_degrees, "Sectors" => sectors,
        "PyramidHeight" => pyramid_height,
        "PyramidHeightOverOuterRing" => TUBE_PYRAMID_HEIGHT_OVER_OUTER_RING,
        "CornerClearanceRule" => "R / tan(phi / 2) + h_K before a semantic corner, phi the " *
                                 "smallest in-plane angle between the tube edges meeting there " *
                                 "(h_K alone for a single tube edge); 0 at box continuation vertices",
        "FacingReach" => facing_reach,
        "FacingRule" => "every hole is wider than 2 x (Radius + PyramidHeight + " *
                        "ProtectedDistance) between any two of its non-adjacent sides, and the " *
                        "vacuum gap between the metal top faces of an upward and a downward " *
                        "process layer exceeds the same 2 x reach, so the tubes facing each " *
                        "other across a hole or a layer gap keep disjoint bands",
        "TubeFrameRule" => "b = (0, 0, Nz) of the tube's process layer, n the horizontal " *
                           "normal away from the metal, e = n x b; the top tube lies on the " *
                           "metal top face at plane + Nz x MetalThickness, the bottom tube on " *
                           "the plane (Tubes[].Layer records Nz)",
        "Top" => Dict("Angles" => top_section.angles, "Materials" => top_section.materials),
        "Bottom" => Dict("Angles" => bottom_section.angles, "Materials" => bottom_section.materials))
    return tools, records, tubes, segments, description
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
    for (k, (tube, section)) in enumerate(tubes)
        segment = segments[(k + 1) ÷ 2]
        start_point = tube_point(tube, 0.0, 0.0, tube.s_start)
        end_point = tube_point(tube, 0.0, 0.0, tube.s_end)
        push!(tube_rows, Dict{String, Any}(
            "Origin" => tube.origin, "Normal" => tube.n, "Extrusion" => tube.e,
            "Start" => tube.s_start, "End" => tube.s_end, "Length" => tube.s_end - tube.s_start,
            "StartPoint" => start_point, "EndPoint" => end_point,
            "EndsOnBox" => [on_box(start_point), on_box(end_point)],
            "Layers" => tube.layers, "Spacing" => tube_spacing(tube),
            "LayerThickness" => layer_records[k],
            "Materials" => section.materials, "Conductor" => segment.conductor,
            "Plane" => segment.plane, "Hole" => segment.hole, "Layer" => segment.layer_sign,
            "CornerAngles" => collect(segment.corner_angles),
            "Edge" => isodd(k) ? "top" : "bottom"))
    end
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
                  "carries a prism tube on its dielectric side: geometric rings of size " *
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
        "SpacingMinimum" => minimum(thicknesses), "SpacingMaximum" => maximum(thicknesses),
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
function curve_segments(curves)
    segments = Tuple{Vector{Float64}, Vector{Float64}}[]
    for curve in curves
        lower, upper = gmsh.model.getParametrizationBounds(1, curve)
        xyz = gmsh.model.getValue(1, curve, [lower[1], upper[1]])
        push!(segments, (collect(xyz[1:3]), collect(xyz[4:6])))
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
        fabricated || scope_error("ThinMetal", "kind thin")
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
    lower, upper = coupon_bounds(edges, radius, metal_thickness, overetch)
    tolerance = 1.0e-7 * radius
    if trace_basis !== nothing
        all(abs(trace_basis.lower[d] - lower[d]) <= tolerance &&
            abs(trace_basis.upper[d] - upper[d]) <= tolerance for d in 1:3) ||
            error("Trace basis box differs from the coupon box")
    end
    outer_tolerance = 1.0e-4 * radius
    validate_plan_view_geometry(edges, radius, tolerance, facets)
    layers = layer_groups(edges, tolerance)
    pullback_metal = metal_thickness / tan(deg2rad(sidewall_angle))
    etch_loops = etch_boundary === nothing ? nothing : read_boundary(etch_boundary)
    if etch_loops !== nothing
        fabricated && sidewall_angle == 90.0 && top_rounding == 0.0 &&
            trench_rounding == 0.0 || error("Explicit etch footprints require sharp vertical fabricated geometry")
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
                    footprint=footprint_polygons
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
    tubes = Tuple{EdgeTube, TubeSection}[]
    tube_segments = NamedTuple[]
    tube_sections = Dict{String, Any}()
    tube_volumes = Int32[]
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
            tube_tools, tube_records, tubes, tube_segments, tube_sections =
                build_edge_tubes!(occ, layers, boundary_loops, etch_loops, semantic_corners,
                                  edge_size, edge_growth_ratio, tube_sector_degrees,
                                  metal_thickness, overetch, corner_isotropy_radius,
                                  lc_tangent, lc_fine, lower, upper, tolerance)
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
        domains, domain_map = occ.fragment(vcat(substrates, vacuum), tools)
        substrate_seed = domain_map[1:length(substrates)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[(length(substrates) + 1):(length(substrates) + length(vacuum))] |>
            Iterators.flatten |>
            collect
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
                owner = nearest_metal_edge(layer_edges, facets, point, radius, tolerance)
                attribute = metal_surface_attribute(5000, owner.slot, owner.conductor)
            elseif !isempty(adjacent_vacuum)
                owner = nearest_metal_edge(layer_edges, facets, point, radius, tolerance)
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
    if prism_tubes
        for (tube, section) in tubes
            volumes = Tuple{Int32, Tuple{Int, Int, Int}, Dict}[]
            for (index, record) in enumerate(tube_records)
                record.tube === tube || continue
                matched = match_tube_entities(tube_volumes[index], tube, section, record.group,
                                              1.0e-6 * radius)
                push!(volumes, (tube_volumes[index], record.group, matched))
                for ((kind, _), tag) in matched
                    kind in (:edge_line, :outer_line, :cap_polygon, :cap_ray) && push!(tube_curves, tag)
                    kind == :edge_line && push!(tube_axis_curves, tag)
                    if kind == :edge_point
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
        prepare_tube_band_sizing!(curve_segments(tube_axis_curves),
                                  curve_segments([curve for curve in feature_curves
                                                  if !(curve in tube_curves)]),
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
            stations, axis_positions, axis_sizes = graded_tube_stations(
                tube.s_start, tube.s_end,
                s -> tube_axis_size(tube_point(tube, 0.0, 0.0, s), lc_tangent, semantic_corners,
                                    corner_grading, corner_grading_slope),
                lc_tangent, edge_growth_ratio;
                surface_start=on_box(tube.s_start), surface_end=on_box(tube.s_end))
            graded = EdgeTube(tube, stations)
            tubes[k] = (graded, section)
            push!(tube_layer_records, tube_layer_statistics(graded, axis_positions, axis_sizes))
            push!(tube_states, TubeMesh(graded, section, tube_volume_groups[k];
                                        pyramid_height=tube_sections["PyramidHeight"]))
        end
    end
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    tube_curve_meshes = [install_tube_curves!(state, next_node, point_nodes) for state in tube_states]
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
            if any(abs(dot(tangent, edge.tangent)) >= 1.0 - 1.0e-6 for edge in edges)
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
        tube_timings["Generate2D"] = @elapsed gmsh.model.mesh.generate(2)
        tube_timings["TubeFaces"] = @elapsed begin
            tube_face_meshes = [install_tube_faces!(state, next_node) for state in tube_states]
        end
        remove_tube_volumes!(tube_states)
        tube_timings["Generate3D"] = @elapsed gmsh.model.mesh.generate(3)
        before_duplicates = sum(length(tags) for tags in gmsh.model.mesh.getElements(3)[2]; init=0)
        gmsh.model.mesh.removeDuplicateElements([(3, tag) for (dim, tag) in gmsh.model.getEntities(3)])
        after_duplicates = sum(length(tags) for tags in gmsh.model.mesh.getElements(3)[2]; init=0)
        before_duplicates == after_duplicates ||
            error("Gmsh produced $(before_duplicates - after_duplicates) duplicate volume elements")
        tube_discrete, tube_volume_census = finalize_tube_volumes!(tube_states)
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
            corner_grading; fixed_node_tags=prism_tubes ? non_simplex_volume_nodes() : UInt[]) :
        nothing
    # Per-type quality of the mixed mesh, gated after the corner-ball optimization
    # (fail closed): positive orientation and Jacobian condition for every type,
    # scaled Jacobian for the tetrahedra.
    mixed_quality = prism_tubes ? mixed_volume_quality() : nothing
    prism_tubes && gate_mixed_volume_quality(mixed_quality, minimum_scaled_jacobian,
                                             maximum_jacobian_condition)
    census_rows = corner_isotropy ?
        seed_corner_census(semantic_corners, corner_grading, tolerance) :
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
                    recipe_scope_record(scope_classes, boundary_loops, lower, upper, tolerance) :
                    nothing,
                "SemanticContract" => semantic_contract,
                "SemanticContractSHA256" => bytes2hex(sha256(read(semantic_contract))),
                "RigidTransform" => vec(transform'),
                "SemanticCorners" => [collect(corner) for corner in semantic_corners],
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
                                for polygon in footprint_polygons]; init=0.0)),
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
        "--far-growth" => ("far_growth", Float64)
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
        far_growth = get(options, "far_growth", 0.0)
    )
end
