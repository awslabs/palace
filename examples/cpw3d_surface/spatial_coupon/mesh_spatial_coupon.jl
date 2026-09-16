# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabrication-resolved local coupon for endpoint, junction, and exact spatial clusters.

import Gmsh: gmsh
using DelimitedFiles
using LinearAlgebra
using SHA

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
    return [merge(edge,(normal_sign=sign(edge.normal_sign),)) for edge in edges]
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

function corner_size_expression(distance, lc_fine, lc_far, radius, transition_width)
    # Isotropic fine size inside the corner ball, then the same linear grading
    # slope as the process edge band up to the far size.
    return "min($(lc_far),$(lc_fine)+($(lc_far)-$(lc_fine))*" *
           "max($(distance)-$(radius),0)/$(transition_width))"
end

# Tangential size of a longitudinal curve through the corner balls: lc_fine
# inside a ball, the process-band grading slope up to lc_tangent outside.
function corner_curve_size(point, corners, radius, lc_fine, lc_tangent, slope)
    distance = minimum(norm(point .- corner) for corner in corners)
    return min(lc_tangent, lc_fine + slope * max(distance - radius, 0.0))
end

# Samples per fine length for the size-weighted arclength quadrature that places
# longitudinal curve nodes; a resolution constant, not a mesh target.
const CURVE_SIZE_SAMPLES_PER_FINE_LENGTH = 16

# Interior node parameters of a longitudinal curve that reaches a corner ball. The
# transfinite lc_tangent grid is kept wherever the corner size law equals
# lc_tangent, so the ridge nodes of one process edge stay aligned across its
# faces away from the corners exactly as without the corner ball (misaligned
# ridge rows were measured to remove the sidewall interior nodes). Each gap
# between kept grid nodes, or between a curve end and the first kept grid node,
# is filled by size-weighted arclength quadrature of the law: lc_fine inside the
# ball, the process-band grading slope up to lc_tangent. Returns nothing when the
# curve is out of reach of every ball, so the caller keeps the ordinary
# transfinite spacing unchanged.
function corner_isotropic_curve_nodes(curve, corners, radius, lc_fine, lc_tangent, slope)
    lower, upper = gmsh.model.getParametrizationBounds(1, curve)
    intervals = max(1, ceil(Int, gmsh.model.occ.getMass(1, curve) / lc_tangent))
    grid = collect(range(lower[1], upper[1]; length=intervals + 1))
    xyz = reshape(gmsh.model.getValue(1, curve, grid), 3, :)
    at_tangent = [corner_curve_size(Tuple(xyz[:, i]), corners, radius, lc_fine, lc_tangent,
                                    slope) >= lc_tangent for i in axes(xyz, 2)]
    all(at_tangent) && return nothing
    anchors = [i for i in 1:(intervals + 1) if at_tangent[i] || i == 1 || i == intervals + 1]
    interior = Float64[]
    for (a, b) in zip(anchors[1:(end - 1)], anchors[2:end])
        if !(b == a + 1 && at_tangent[a] && at_tangent[b])
            append!(interior, graded_curve_parameters(curve, grid[a], grid[b], corners, radius,
                                                      lc_fine, lc_tangent, slope))
        end
        b <= intervals && push!(interior, grid[b])
    end
    return interior, gmsh.model.getValue(1, curve, interior)
end

# Interior parameters of the curve span [from, to] equidistributed in the
# arclength integral of 1 / corner size law.
function graded_curve_parameters(curve, from, to, corners, radius, lc_fine, lc_tangent, slope)
    samples = max(64, ceil(Int, CURVE_SIZE_SAMPLES_PER_FINE_LENGTH *
                              norm(gmsh.model.getValue(1, curve, [to]) .-
                                   gmsh.model.getValue(1, curve, [from])) / lc_fine))
    parameters = collect(range(from, to; length=samples + 1))
    xyz = reshape(gmsh.model.getValue(1, curve, parameters), 3, :)
    sizes = [corner_curve_size(Tuple(xyz[:, i]), corners, radius, lc_fine, lc_tangent, slope)
             for i in axes(xyz, 2)]
    cumulative = zeros(samples + 1)
    for i in 1:samples
        cumulative[i + 1] = cumulative[i] + norm(xyz[:, i + 1] .- xyz[:, i]) *
                                            0.5 * (1.0 / sizes[i] + 1.0 / sizes[i + 1])
    end
    intervals = max(1, round(Int, cumulative[end]))
    interior = Float64[]
    for k in 1:(intervals - 1)
        target = k * cumulative[end] / intervals
        i = clamp(searchsortedlast(cumulative, target), 1, samples)
        fraction = (target - cumulative[i]) / (cumulative[i + 1] - cumulative[i])
        push!(interior, parameters[i] + fraction * (parameters[i + 1] - parameters[i]))
    end
    return interior
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

function tetrahedron_aspect(xyz)
    jacobian = hcat(xyz[2] .- xyz[1], xyz[3] .- xyz[1], xyz[4] .- xyz[1])
    singular = svdvals(jacobian)
    return singular[1] / singular[end]
end

# Edge lengths and cell aspects of the linear seed inside each semantic corner ball.
function seed_corner_census(corners, radius, isotropic_size, tolerance)
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    points = reshape(coordinates, 3, :)
    index = Dict(tag => i for (i, tag) in enumerate(node_tags))
    _, element_tags, element_nodes = gmsh.model.mesh.getElements(3)
    tetrahedra = Vector{NTuple{4, Int}}()
    for (tags, block) in zip(element_tags, element_nodes)
        isempty(tags) && continue
        nodes_per_element = length(block) ÷ length(tags)
        # Only the vertex nodes define the linear cells (high-order nodes follow).
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
            "RingEdgesOverSqrt2IsotropicSize" => count(>(threshold), ring_radii)))
    end
    return rows
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
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, attribute)
            _, element_tags, element_nodes = gmsh.model.mesh.getElements(2, entity)
            for (tags, block) in zip(element_tags, element_nodes)
                isempty(tags) && continue
                nodes_per_element = length(block) ÷ length(tags)
                for start in 1:nodes_per_element:length(block)
                    a, b, c = (points[:, index[block[start + i]]] for i in 0:2)
                    area += 0.5 * norm(cross(b .- a, c .- a))
                    triangles += 1
                end
            end
        end
        push!(rows, Dict{String, Any}("Attribute" => Int(attribute),
                                      "Name" => gmsh.model.getPhysicalName(dim, attribute),
                                      "Triangles" => triangles, "Area" => area))
    end
    sort!(rows; by=row -> row["Attribute"])
    return rows
end

# Bins of the along-edge histogram of a longitudinal face's interior nodes; a
# census resolution constant, not a mesh target.
const LONGITUDINAL_FACE_HISTOGRAM_BINS = 10

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
        _, element_tags, element_nodes = gmsh.model.mesh.getElements(2, surface)
        triangles = Vector{NTuple{3, Int}}()
        for (tags, block) in zip(element_tags, element_nodes)
            isempty(tags) && continue
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
    if edge.vertex_arm
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

function loft_mask_offsets(occ, loops, z0, z1, bottom_offset, top_offset, tolerance)
    outers = [loop for loop in loops if !loop.hole]
    holes = [loop for loop in loops if loop.hole]
    isempty(outers) && error("Plan-view mask has no exterior loop")
    result = Tuple{Int32, Int32}[]
    hole_owners=zeros(Int,length(holes))
    for outer in outers
        volume = loft_polygon(
            occ,
            offset_loop_points(outer, bottom_offset, tolerance),
            offset_loop_points(outer, top_offset, tolerance),
            z0,
            z1
        )
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

function loft_mask(occ, loops, z0, z1, pullback, tolerance)
    return loft_mask_offsets(occ, loops, z0, z1, 0.0, pullback, tolerance)
end

function boundary_strips(occ, loops, radius, z0, z1, pullback, tolerance)
    expanded_volumes = Tuple{Int32, Int32}[]
    retained_volumes = Tuple{Int32, Int32}[]
    width = 3radius
    for conductor in sort!(unique(loop.conductor for loop in loops))
        conductor_loops = [loop for loop in loops if loop.conductor == conductor]
        append!(expanded_volumes,
                loft_mask_offsets(occ, conductor_loops, z0, z1, -width, -width, tolerance))
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

function loft_strip(occ, edge, radius, side, z0, z1, pullback)
    bottom = polygon_wire(occ, strip_points(edge, radius, side), z0)
    top = polygon_wire(occ, strip_points(edge, radius, side, pullback), z1)
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
    lower = (
        lower[1],
        lower[2],
        min(lower[3], minimum(edge.point[3] for edge in edges) - radius - overetch)
    )
    upper = (
        upper[1],
        upper[2],
        max(upper[3], maximum(edge.point[3] for edge in edges) + radius + metal_thickness)
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
    (lc_tangent == 0.0 || lc_fine <= lc_tangent <= lc_far) ||
        error("tangential mesh size must lie between fine and far sizes")
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
    # The seed honors the same isotropic corner ball the metric stage prescribes:
    # NormalSize (lc_fine) inside CornerIsotropyRadius around every contract
    # semantic corner, graded to the far size with the process-band slope.
    semantic_corners = corner_isotropy ? read_semantic_corners(semantic_contract, transform) :
                       NTuple{3, Float64}[]

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
                          layer.plane - layer.sign * overetch, 0.0, tolerance)
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
                            pullback_trench
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
                    tolerance
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
        domains, domain_map = occ.fragment(vcat(substrates, vacuum), [])
        substrate_seed = domain_map[1:length(substrates)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[(length(substrates) + 1):end] |> Iterators.flatten |> collect
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

        point = point_on_surface(tag)
        edge = nearest_edge(edges, point, radius)
        attribute = 0
        if fabricated
            if !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                _, _, zmin, _, _, zmax = bounds
                attribute =
                    abs(zmin - edge.point[3]) < tolerance &&
                    abs(zmax - edge.point[3]) < tolerance ? 3000 + edge.slot :
                    3100 + edge.slot
            elseif !isempty(adjacent_substrate)
                owner = nearest_metal_edge(edges, facets, point, radius, tolerance)
                attribute = metal_surface_attribute(5000, owner.slot, owner.conductor)
            elseif !isempty(adjacent_vacuum)
                owner = nearest_metal_edge(edges, facets, point, radius, tolerance)
                attribute = metal_surface_attribute(6000, owner.slot, owner.conductor)
            end
        else
            metal_edges = [
                candidate for candidate in edges if
                point_in_metal(candidate, point, radius, tolerance, facets)
            ]
            if !isempty(metal_edges)
                owner = nearest_edge(metal_edges, point, radius)
                attribute = metal_surface_attribute(4000, owner.slot, owner.conductor)
            elseif !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                attribute = 3000 + edge.slot
            end
        end
        attribute > 0 && push!(get!(boundary_groups, attribute, Int32[]), tag)
    end
    isempty(matching) && error("No matching surface was generated")

    gmsh.model.addPhysicalGroup(3, substrate_tags, 1, "substrate")
    gmsh.model.addPhysicalGroup(3, vacuum_tags, 2, "vacuum")
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
    longitudinal_curves = Int32[]
    corner_curves = Int32[]
    corner_grading_slope = (lc_far - lc_fine) / (process_core_width - process_fine_width)
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    if lc_tangent > 0.0
        for curve in feature_curves
            lower_parameter, upper_parameter =
                gmsh.model.getParametrizationBounds(1, curve)
            parameter = 0.5 * (lower_parameter[1] + upper_parameter[1])
            derivative = gmsh.model.getDerivative(1, curve, [parameter])
            tangent = derivative[1:3]
            tangent ./= norm(tangent)
            if any(abs(dot(tangent, edge.tangent)) >= 1.0 - 1.0e-6 for edge in edges)
                push!(longitudinal_curves, curve)
                placed = corner_isotropy ? corner_isotropic_curve_nodes(
                    curve, semantic_corners, corner_isotropy_radius, lc_fine, lc_tangent,
                    corner_grading_slope) : nothing
                if placed === nothing
                    curve_length = gmsh.model.occ.getMass(1, curve)
                    point_count = max(2, ceil(Int, curve_length / lc_tangent) + 1)
                    gmsh.model.mesh.setTransfiniteCurve(curve, point_count)
                else
                    # Transfinite spacing cannot follow the corner ball; place the
                    # curve nodes explicitly and keep them through generation.
                    add_explicit_curve_mesh!(curve, placed..., next_node, point_nodes)
                    push!(corner_curves, curve)
                end
            end
        end
    end
    isempty(corner_curves) || gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
    println(
        "Spatial mesh features: candidates=$(length(candidate_curves)), " *
        "physical=$(length(feature_curves)), " *
        "longitudinal=$(length(longitudinal_curves)), " *
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
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(feature_curves))
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
            gmsh.model.mesh.field.setNumbers(2, "PointsList", Float64.(corner_point_tags))
            gmsh.model.mesh.field.add("MathEval", 3)
            gmsh.model.mesh.field.setString(3, "F", "min(F1," * corner_size_expression(
                "F2", lc_fine, lc_far, corner_isotropy_radius,
                process_core_width - process_fine_width) * ")")
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
                "F3", lc_fine, lc_far, corner_isotropy_radius, transition_width) * ")"
        end
        gmsh.model.mesh.field.add("MathEval", 2)
        gmsh.model.mesh.field.setString(2, "F", size_expression)
        gmsh.model.mesh.field.setAsBackgroundMesh(2)
    end
    for (name, value) in [
        ("Mesh.MeshSizeMin", lc_fine),
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
    gmsh.model.mesh.generate(3)
    # Reject oversized linear meshes before allocating their high-order nodes.
    _, linear_tags, _ = gmsh.model.mesh.getElements(3)
    sum(length(tags) for tags in linear_tags) <= max_elements ||
        error("Linear spatial coupon exceeds element budget before order elevation")
    if lc_tangent == 0.0 && optimize_volume
        gmsh.model.mesh.optimize("Netgen")
    end
    census_rows = corner_isotropy ?
        seed_corner_census(semantic_corners, corner_isotropy_radius, lc_fine, tolerance) :
        Dict{String, Any}[]
    corner_reach = corner_isotropy ?
        corner_law_reach(corner_isotropy_radius, lc_fine, lc_tangent, corner_grading_slope) :
        0.0
    face_rows = corner_isotropy ?
        longitudinal_face_census(longitudinal_curves, semantic_corners, corner_reach) :
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
    minimum_signed_inverse_condition=minimum(
        gmsh.model.mesh.getElementQualities(all_volume_tags,"minSICN"))
    minimum_signed_inverse_condition>1e-10 ||
        error("Spatial coupon has invalid or near-singular elements: minSICN=$minimum_signed_inverse_condition")
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
                "Scope" => "Seed corner-ball census, longitudinal-face census and interface areas; reported, not a qualification gate",
                "SemanticContract" => semantic_contract,
                "SemanticContractSHA256" => bytes2hex(sha256(read(semantic_contract))),
                "RigidTransform" => vec(transform'),
                "SemanticCorners" => [collect(corner) for corner in semantic_corners],
                "CornerIsotropyRadius" => corner_isotropy_radius,
                "IsotropicSize" => lc_fine,
                "Sqrt2IsotropicSize" => sqrt(2.0) * lc_fine,
                "FarSize" => lc_far,
                "GradingTransitionWidth" => process_core_width - process_fine_width,
                "LongitudinalCurves" => length(longitudinal_curves),
                "CornerIsotropicLongitudinalCurves" => length(corner_curves),
                "CornerLawReach" => corner_reach,
                "LongitudinalFaceHistogramBins" => LONGITUDINAL_FACE_HISTOGRAM_BINS,
                "LongitudinalFaces" => face_rows,
                "InterfaceAreaUnits" => "um^2",
                "EtchBoundary" => etch_boundary === nothing ? "producer-default" : etch_boundary,
                "EtchBoundarySHA256" => etch_boundary === nothing ? nothing :
                                        bytes2hex(sha256(read(etch_boundary))),
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
        "--etch-boundary" => ("etch_boundary", String)
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
        etch_boundary   = get(options, "etch_boundary", nothing)
    )
end
