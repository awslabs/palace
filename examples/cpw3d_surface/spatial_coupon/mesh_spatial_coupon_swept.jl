# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Explicit swept-prism mesher for a translationally invariant spatial edge coupon.
#
# This is intentionally a separate executable from mesh_spatial_coupon.jl. The latter
# remains the experimental field-based tetrahedral mesher. This first swept milestone
# accepts one continuing straight edge and constructs the complete cross-section from
# explicit coordinates before sweeping it with independent longitudinal coordinates.

include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

mutable struct CrossSectionBuilder
    geo
    edge
    s0::Float64
    grading_width::Float64
    point_tags::Dict{NTuple{2, Float64}, Int32}
    point_coordinates::Dict{Int32, NTuple{2, Float64}}
    line_tags::Dict{Tuple{Int32, Int32}, Int32}
    line_orientations::Dict{Tuple{Int32, Int32}, Tuple{Int32, Int32}}
    arc_centers::Dict{Tuple{NTuple{2, Float64}, NTuple{2, Float64}}, NTuple{2, Float64}}
    surfaces::Vector{Tuple{Int32, Symbol}}
    near_feature_spacings::Vector{Float64}
end

coordinate_key(value) = round(Float64(value), digits=14)
point_key(u, w) = (coordinate_key(u), coordinate_key(w))
edge_key(first, second) = first < second ? (first, second) : (second, first)
coordinate_edge_key(first, second) = first < second ? (first, second) : (second, first)

function global_point(edge, u, s, w)
    sign = edge.normal_sign
    return (
        edge.point[1] + u * edge.gap[1] + s * edge.tangent[1],
        edge.point[2] + u * edge.gap[2] + s * edge.tangent[2],
        edge.point[3] + sign * w
    )
end

function add_cross_section_point!(builder, u, w)
    key = point_key(u, w)
    return get!(builder.point_tags, key) do
        point = global_point(builder.edge, u, builder.s0, w)
        tag = builder.geo.addPoint(point...)
        builder.point_coordinates[tag] = key
        return tag
    end
end

function add_arc_center!(builder, center)
    key = point_key(center...)
    return get!(builder.point_tags, key) do
        point = global_point(builder.edge, center[1], builder.s0, center[2])
        tag = builder.geo.addPoint(point...)
        builder.point_coordinates[tag] = key
        return tag
    end
end

function oriented_curve!(builder, first, second)
    p0 = add_cross_section_point!(builder, first...)
    p1 = add_cross_section_point!(builder, second...)
    key = edge_key(p0, p1)
    if !haskey(builder.line_tags, key)
        coordinate_key_pair = coordinate_edge_key(point_key(first...), point_key(second...))
        builder.line_tags[key] = if haskey(builder.arc_centers, coordinate_key_pair)
            center = add_arc_center!(builder, builder.arc_centers[coordinate_key_pair])
            builder.geo.addCircleArc(p0, center, p1)
        else
            builder.geo.addLine(p0, p1)
        end
        builder.line_orientations[key] = (p0, p1)
        builder.geo.mesh.setTransfiniteCurve(builder.line_tags[key], 2)
        return builder.line_tags[key]
    end
    curve = builder.line_tags[key]
    return builder.line_orientations[key] == (p0, p1) ? curve : -curve
end

function add_triangle!(builder, first, second, third, material)
    points = (first, second, third)
    area = 0.5 * (
        first[1] * (second[2] - third[2]) +
        second[1] * (third[2] - first[2]) +
        third[1] * (first[2] - second[2])
    )
    abs(area) > 1.0e-15 || error("Swept cross-section contains a degenerate triangle")
    if area < 0.0
        points = (first, third, second)
    end
    point_tags = [add_cross_section_point!(builder, point...) for point in points]
    curves = Int32[
        oriented_curve!(builder, points[1], points[2]),
        oriented_curve!(builder, points[2], points[3]),
        oriented_curve!(builder, points[3], points[1])
    ]
    loop = builder.geo.addCurveLoop(curves)
    surface = builder.geo.addPlaneSurface([loop])
    builder.geo.mesh.setTransfiniteSurface(surface, "Left", point_tags)
    push!(builder.surfaces, (surface, material))
    return surface
end

function graded_distances(width, lc_normal, lc_far, power, grading_width)
    width >= 0.0 || error("Graded-coordinate width must be nonnegative")
    width == 0.0 && return [0.0]
    width <= lc_normal && return [0.0, width]
    coordinates = [0.0, lc_normal]
    distance = lc_normal
    step = lc_normal
    while true
        fraction = min(distance / grading_width, 1.0)
        step = lc_normal + (lc_far - lc_normal) * fraction^power
        remaining = width - distance
        remaining <= step && break
        if remaining - step < 0.45 * min(step, lc_far)
            break
        end
        distance += step
        push!(coordinates, distance)
    end
    remaining = width - coordinates[end]
    if length(coordinates) >= 2 &&
       remaining < 0.45 * (coordinates[end] - coordinates[end - 1])
        previous = coordinates[end - 1]
        total = width - previous
        coordinates[end] = previous + 0.5 * total
    end
    push!(coordinates, width)
    return unique(coordinates)
end

function graded_interval(
    first,
    last,
    lc_normal,
    lc_far,
    power,
    grading_width;
    refine=:left
)
    last > first || error("Graded interval must have positive width")
    distances = graded_distances(
        last - first,
        lc_normal,
        lc_far,
        power,
        grading_width
    )
    if refine == :left
        return first .+ distances
    elseif refine == :right
        return last .- reverse(distances)
    end
    error("Unknown graded interval refinement side: $refine")
end

function anchored_interval(
    first,
    anchor,
    last,
    lc_normal,
    lc_far,
    power,
    grading_width
)
    first < anchor < last || error("Anchor must lie inside its interval")
    left = graded_interval(
        first,
        anchor,
        lc_normal,
        lc_far,
        power,
        grading_width;
        refine=:right
    )
    right = graded_interval(
        anchor,
        last,
        lc_normal,
        lc_far,
        power,
        grading_width;
        refine=:left
    )
    return vcat(left[1:(end - 1)], right)
end

function record_near_spacing!(builder, coordinates, side)
    length(coordinates) >= 2 || return
    spacing = side == :left ? coordinates[2] - coordinates[1] :
              coordinates[end] - coordinates[end - 1]
    push!(builder.near_feature_spacings, spacing)
end

function row_from_left!(builder, left, right, lc_normal, lc_far, power)
    values = graded_interval(
        left,
        right,
        lc_normal,
        lc_far,
        power,
        builder.grading_width;
        refine=:left
    )
    record_near_spacing!(builder, values, :left)
    return values
end

function row_from_right!(builder, left, right, lc_normal, lc_far, power)
    values = graded_interval(
        left,
        right,
        lc_normal,
        lc_far,
        power,
        builder.grading_width;
        refine=:right
    )
    record_near_spacing!(builder, values, :right)
    return values
end

function row_from_anchor!(builder, left, anchor, right, lc_normal, lc_far, power)
    values = anchored_interval(
        left,
        anchor,
        right,
        lc_normal,
        lc_far,
        power,
        builder.grading_width
    )
    index = findfirst(==(anchor), values)
    index === nothing && error("Unable to recover the explicit normal-coordinate anchor")
    index > 1 && push!(builder.near_feature_spacings, values[index] - values[index - 1])
    index < length(values) &&
        push!(builder.near_feature_spacings, values[index + 1] - values[index])
    return values
end

function zipper_strip!(builder, lower, upper, material)
    z0, x0 = lower
    z1, x1 = upper
    z1 > z0 || error("Swept cross-section rows must be strictly ordered")
    length(x0) >= 2 && length(x1) >= 2 ||
        error("Swept cross-section row contains fewer than two points")
    scale0 = x0[end] - x0[1]
    scale1 = x1[end] - x1[1]
    scale0 > 0.0 && scale1 > 0.0 || error("Swept cross-section row has zero width")
    i = 1
    j = 1
    while i < length(x0) || j < length(x1)
        if i == length(x0)
            add_triangle!(builder, (x0[i], z0), (x1[j + 1], z1), (x1[j], z1), material)
            j += 1
        elseif j == length(x1)
            add_triangle!(builder, (x0[i], z0), (x0[i + 1], z0), (x1[j], z1), material)
            i += 1
        else
            next0 = (x0[i + 1] - x0[1]) / scale0
            next1 = (x1[j + 1] - x1[1]) / scale1
            if next0 <= next1
                add_triangle!(
                    builder,
                    (x0[i], z0),
                    (x0[i + 1], z0),
                    (x1[j], z1),
                    material
                )
                i += 1
            else
                add_triangle!(
                    builder,
                    (x0[i], z0),
                    (x1[j + 1], z1),
                    (x1[j], z1),
                    material
                )
                j += 1
            end
        end
    end
    return
end

function add_region!(builder, rows, material)
    length(rows) >= 2 || return
    for (lower, upper) in zip(rows, rows[2:end])
        zipper_strip!(builder, lower, upper, material)
    end
    return
end

function uniform_segment(first, last, size)
    length = hypot(last[1] - first[1], last[2] - first[2])
    count = max(1, ceil(Int, length / size))
    return [
        (
            first[1] + index * (last[1] - first[1]) / count,
            first[2] + index * (last[2] - first[2]) / count
        ) for index = 0:count
    ]
end

function fillet_data(vertex, first_direction, second_direction, radius, lc_normal)
    radius <= 0.0 && return nothing
    d0 = first_direction ./ hypot(first_direction...)
    d1 = second_direction ./ hypot(second_direction...)
    angle = acos(clamp(d0[1] * d1[1] + d0[2] * d1[2], -1.0, 1.0))
    1.0e-8 < angle < pi - 1.0e-8 || error("Invalid straight-edge fillet angle")
    setback = radius / tan(0.5 * angle)
    first = (vertex[1] + setback * d0[1], vertex[2] + setback * d0[2])
    second = (vertex[1] + setback * d1[1], vertex[2] + setback * d1[2])
    bisector = (d0[1] + d1[1], d0[2] + d1[2])
    bisector_norm = hypot(bisector...)
    center_distance = radius / sin(0.5 * angle)
    center = (
        vertex[1] + center_distance * bisector[1] / bisector_norm,
        vertex[2] + center_distance * bisector[2] / bisector_norm
    )
    radial_first = (first[1] - center[1], first[2] - center[2])
    radial_second = (second[1] - center[1], second[2] - center[2])
    delta = atan(
        radial_first[1] * radial_second[2] - radial_first[2] * radial_second[1],
        radial_first[1] * radial_second[1] + radial_first[2] * radial_second[2]
    )
    count = max(1, ceil(Int, abs(delta) * radius / lc_normal))
    start_angle = atan(radial_first[2], radial_first[1])
    points = [
        (
            center[1] + radius * cos(start_angle + index * delta / count),
            center[2] + radius * sin(start_angle + index * delta / count)
        ) for index = 0:count
    ]
    return (first=first, second=second, center=center, points=points)
end

function register_arc!(builder, points, center)
    for (first, second) in zip(points, points[2:end])
        key = coordinate_edge_key(point_key(first...), point_key(second...))
        builder.arc_centers[key] = point_key(center...)
    end
    return
end

function unique_path(points)
    result = NTuple{2, Float64}[]
    for point in points
        value = (Float64(point[1]), Float64(point[2]))
        distance = if isempty(result)
            Inf
        else
            hypot(value[1] - result[end][1], value[2] - result[end][2])
        end
        distance > 1.0e-13 && push!(result, value)
    end
    return result
end

function fabricated_boundaries!(
    builder,
    metal_thickness,
    overetch,
    sidewall_angle,
    top_rounding,
    trench_rounding,
    lc_normal
)
    slope = cot(deg2rad(sidewall_angle))
    edge = (0.0, 0.0)

    trench_vertex = (overetch * slope, -overetch)
    trench_side_direction = (
        edge[1] - trench_vertex[1],
        edge[2] - trench_vertex[2]
    )
    trench_fillet = fillet_data(
        trench_vertex,
        (1.0, 0.0),
        trench_side_direction,
        trench_rounding,
        lc_normal
    )
    trench_bottom = trench_fillet === nothing ? trench_vertex : trench_fillet.first
    trench_side = trench_fillet === nothing ? trench_vertex : trench_fillet.second
    trench_path = NTuple{2, Float64}[trench_bottom]
    if trench_fillet !== nothing
        append!(trench_path, trench_fillet.points[2:end])
        register_arc!(builder, trench_fillet.points, trench_fillet.center)
    end
    side_points = uniform_segment(trench_side, edge, lc_normal)
    append!(trench_path, side_points[2:end])

    metal_vertex = (-metal_thickness * slope, metal_thickness)
    metal_side_direction = (
        edge[1] - metal_vertex[1],
        edge[2] - metal_vertex[2]
    )
    metal_fillet = fillet_data(
        metal_vertex,
        metal_side_direction,
        (-1.0, 0.0),
        top_rounding,
        lc_normal
    )
    metal_side = metal_fillet === nothing ? metal_vertex : metal_fillet.first
    metal_top = metal_fillet === nothing ? metal_vertex : metal_fillet.second
    metal_path = uniform_segment(edge, metal_side, lc_normal)
    if metal_fillet !== nothing
        append!(metal_path, metal_fillet.points[2:end])
        register_arc!(builder, metal_fillet.points, metal_fillet.center)
    end

    trench_path = unique_path(trench_path)
    metal_path = unique_path(metal_path)
    trench_is_monotone = all(
        trench_path[index + 1][2] > trench_path[index][2] for
        index = 1:(length(trench_path) - 1)
    )
    trench_is_monotone ||
        error("Trench boundary is not monotone in the fabrication coordinate")
    metal_is_monotone = all(
        metal_path[index + 1][2] > metal_path[index][2] for
        index = 1:(length(metal_path) - 1)
    )
    metal_is_monotone ||
        error("Metal boundary is not monotone in the fabrication coordinate")
    return trench_path, metal_path, trench_bottom[1], metal_top[1]
end

function build_thin_cross_section!(
    builder,
    umin,
    umax,
    wmin,
    wmax,
    lc_normal,
    lc_far,
    power
)
    x = row_from_anchor!(builder, umin, 0.0, umax, lc_normal, lc_far, power)
    lower_w = graded_interval(
        wmin,
        0.0,
        lc_normal,
        lc_far,
        power,
        builder.grading_width;
        refine=:right
    )
    upper_w = graded_interval(
        0.0,
        wmax,
        lc_normal,
        lc_far,
        power,
        builder.grading_width;
        refine=:left
    )
    record_near_spacing!(builder, lower_w, :right)
    record_near_spacing!(builder, upper_w, :left)
    add_region!(builder, [(w, x) for w in lower_w], :substrate)
    add_region!(builder, [(w, x) for w in upper_w], :vacuum)
    return
end

function build_fabricated_cross_section!(
    builder,
    umin,
    umax,
    wmin,
    wmax,
    metal_thickness,
    overetch,
    sidewall_angle,
    top_rounding,
    trench_rounding,
    lc_normal,
    lc_far,
    power
)
    trench, metal, trench_anchor, metal_anchor = fabricated_boundaries!(
        builder,
        metal_thickness,
        overetch,
        sidewall_angle,
        top_rounding,
        trench_rounding,
        lc_normal
    )

    lower_w = graded_interval(
        wmin,
        -overetch,
        lc_normal,
        lc_far,
        power,
        builder.grading_width;
        refine=:right
    )
    lower_rows = [
        (
            w,
            row_from_anchor!(
                builder,
                umin,
                trench_anchor,
                umax,
                lc_normal,
                lc_far,
                power
            )
        ) for w in lower_w
    ]
    add_region!(builder, lower_rows, :substrate)

    substrate_trench_rows = [
        (
            point[2],
            row_from_right!(builder, umin, point[1], lc_normal, lc_far, power)
        ) for point in trench
    ]
    vacuum_trench_rows = [
        (
            point[2],
            row_from_left!(builder, point[1], umax, lc_normal, lc_far, power)
        ) for point in trench
    ]
    add_region!(builder, substrate_trench_rows, :substrate)
    add_region!(builder, vacuum_trench_rows, :vacuum)

    vacuum_metal_rows = [
        (
            point[2],
            row_from_left!(builder, point[1], umax, lc_normal, lc_far, power)
        ) for point in metal
    ]
    add_region!(builder, vacuum_metal_rows, :vacuum)

    upper_w = graded_interval(
        metal_thickness,
        wmax,
        lc_normal,
        lc_far,
        power,
        builder.grading_width;
        refine=:left
    )
    upper_rows = [
        (
            w,
            row_from_anchor!(
                builder,
                umin,
                metal_anchor,
                umax,
                lc_normal,
                lc_far,
                power
            )
        ) for w in upper_w
    ]
    add_region!(builder, upper_rows, :vacuum)
    return
end

function local_coordinates(edge, point)
    delta = (
        point[1] - edge.point[1],
        point[2] - edge.point[2],
        point[3] - edge.point[3]
    )
    return (
        delta[1] * edge.gap[1] + delta[2] * edge.gap[2],
        delta[1] * edge.tangent[1] + delta[2] * edge.tangent[2],
        edge.normal_sign * delta[3]
    )
end

function surface_local_nodes(edge, tag)
    _, coordinates, _ = gmsh.model.mesh.getNodes(2, tag, true, false)
    isempty(coordinates) && error("Surface $tag has no mesh nodes")
    return [
        local_coordinates(
            edge,
            (
                coordinates[index],
                coordinates[index + 1],
                coordinates[index + 2]
            )
        ) for index = 1:3:length(coordinates)
    ]
end

function on_outer_surface(points, umin, umax, smin, smax, wmin, wmax, tolerance)
    for (coordinate, lower, upper) in (
        (1, umin, umax),
        (2, smin, smax),
        (3, wmin, wmax)
    )
        values = [point[coordinate] for point in points]
        (all(abs(value - lower) <= tolerance for value in values) ||
         all(abs(value - upper) <= tolerance for value in values)) && return true
    end
    return false
end

function prism_face_counts()
    counts = Dict{Tuple{Vararg{UInt64}}, Int}()
    volume_types, _, _ = gmsh.model.mesh.getElements(3)
    for element_type in volume_types
        name, _, _, node_count, _, primary_count =
            gmsh.model.mesh.getElementProperties(element_type)
        occursin("Prism", name) ||
            error("Swept mesh generated unsupported element type $name")
        primary_count == 6 || error("Unexpected primary-node count for $name")
        _, connectivity = gmsh.model.mesh.getElementsByType(element_type)
        for offset = 1:node_count:length(connectivity)
            nodes = connectivity[offset:(offset + primary_count - 1)]
            faces = (
                (nodes[1], nodes[2], nodes[3]),
                (nodes[4], nodes[5], nodes[6]),
                (nodes[1], nodes[2], nodes[5], nodes[4]),
                (nodes[2], nodes[3], nodes[6], nodes[5]),
                (nodes[3], nodes[1], nodes[4], nodes[6])
            )
            for face in faces
                key = Tuple(sort!(UInt64[face...]))
                counts[key] = get(counts, key, 0) + 1
            end
        end
    end
    nonmanifold = count(>(2), values(counts))
    nonmanifold == 0 || error("Swept mesh contains $nonmanifold nonmanifold faces")
    return (
        boundary=count(==(1), values(counts)),
        interior=count(==(2), values(counts)),
        nonmanifold=nonmanifold
    )
end

function json_array(values)
    return "[" * join(string.(values), ", ") * "]"
end

function write_swept_metadata(
    path;
    node_count,
    element_count,
    surface_attributes,
    volume_attributes,
    volume_element_types,
    edge_count,
    lc_normal,
    lc_tangent,
    lc_far,
    core_width,
    grading_power,
    longitudinal_coordinates,
    normal_spacings,
    minimum_jacobian,
    mesh_order,
    face_counts,
    cross_section_nodes,
    fabricated
)
    open(path, "w") do stream
        println(stream, "{")
        println(stream, "  \"Version\": 2,")
        println(stream, "  \"MeshingMode\": \"Swept\",")
        println(stream, "  \"SweptTopology\": \"StraightEdgePrisms\",")
        println(stream, "  \"Fabricated\": $(fabricated ? "true" : "false"),")
        println(stream, "  \"MetalSurfacePartition\": \"InterfaceSlotAndConductor\",")
        println(stream, "  \"NodeCount\": $node_count,")
        println(stream, "  \"VolumeElementCount\": $element_count,")
        println(stream, "  \"SweptVolumeElementCount\": $element_count,")
        println(stream, "  \"TransitionVolumeElementCount\": 0,")
        println(stream, "  \"SurfaceAttributes\": $(json_array(surface_attributes)),")
        println(stream, "  \"VolumeAttributes\": $(json_array(volume_attributes)),")
        println(
            stream,
            "  \"VolumeElementTypes\": [" *
            join(("\"" * value * "\"" for value in volume_element_types), ", ") *
            "],"
        )
        println(stream, "  \"InputEdgeCount\": $edge_count,")
        println(stream, "  \"FineSize\": $lc_normal,")
        println(stream, "  \"NormalSize\": $lc_normal,")
        println(stream, "  \"TangentialSize\": $lc_tangent,")
        println(stream, "  \"FarSize\": $lc_far,")
        println(stream, "  \"ProcessCoreWidth\": $core_width,")
        println(stream, "  \"ProcessGradingPower\": $grading_power,")
        println(
            stream,
            "  \"LongitudinalStationCount\": $(length(longitudinal_coordinates)),"
        )
        println(
            stream,
            "  \"LongitudinalCoordinates\": $(json_array(longitudinal_coordinates)),"
        )
        maximum_tangent = maximum(diff(longitudinal_coordinates))
        println(stream, "  \"MeasuredMaximumLongitudinalSpacing\": $maximum_tangent,")
        println(
            stream,
            "  \"MeasuredNearFeatureNormalSpacings\": " *
            "$(json_array(normal_spacings)),"
        )
        println(
            stream,
            "  \"MeasuredMaximumNearFeatureNormalSpacing\": " *
            "$(maximum(normal_spacings)),"
        )
        println(stream, "  \"CrossSectionNodeCount\": $cross_section_nodes,")
        println(stream, "  \"BoundaryFaceCount\": $(face_counts.boundary),")
        println(stream, "  \"InteriorFaceCount\": $(face_counts.interior),")
        println(stream, "  \"NonmanifoldFaceCount\": $(face_counts.nonmanifold),")
        println(stream, "  \"MinimumScaledJacobian\": $minimum_jacobian,")
        println(stream, "  \"MeshOrder\": $mesh_order")
        println(stream, "}")
    end
    return
end

function generate_spatial_coupon_swept(;
    signature::String,
    mask::Union{Nothing, String}=nothing,
    boundary::Union{Nothing, String}=nothing,
    fabricated::Bool,
    radius::Float64=2.0,
    metal_thickness::Float64=0.1,
    overetch::Float64=0.05,
    sidewall_angle::Float64=80.0,
    top_rounding::Float64=0.01,
    trench_rounding::Float64=0.01,
    lc_normal::Float64=0.02,
    lc_tangent::Float64=0.2,
    lc_far::Float64=0.3,
    process_core_width::Float64=0.0,
    process_grading_power::Float64=1.7,
    max_nodes::Int=500_000,
    max_elements::Int=2_000_000,
    mesh_order::Int=2,
    filename::String
)
    radius > 0.0 || error("radius must be positive")
    metal_thickness > 0.0 || error("metal thickness must be positive")
    0.0 <= overetch < radius || error("overetch must lie in [0, radius)")
    0.0 < sidewall_angle <= 90.0 || error("sidewall angle must lie in (0, 90]")
    0.0 <= top_rounding < metal_thickness ||
        error("top rounding must be smaller than metal thickness")
    0.0 <= trench_rounding <= overetch ||
        error("trench rounding must not exceed overetch")
    lc_normal > 0.0 || error("normal mesh size must be positive")
    lc_tangent > 0.0 || error("tangential mesh size must be positive")
    lc_far >= lc_normal || error("far mesh size must not be smaller than normal size")
    process_grading_power > 0.0 || error("process grading power must be positive")
    max_nodes > 0 || error("maximum node budget must be positive")
    max_elements > 0 || error("maximum element budget must be positive")
    mesh_order >= 2 || error("swept spatial coupons require quadratic mesh order or higher")
    mask === nothing ||
        error("The straight-edge swept milestone does not accept a plan-view mask")
    boundary === nothing ||
        error("The straight-edge swept milestone does not accept a classified boundary")

    edges = read_edges(signature)
    length(edges) == 1 || error(
        "The straight-edge swept milestone requires exactly one edge; " *
        "endpoint, corner, and interacting-edge transition cores are not implemented"
    )
    edge = only(edges)
    edge.vertex_arm && error(
        "The straight-edge swept milestone requires a continuing edge, not a vertex arm"
    )
    abs(edge.interval[1] - edge.interval[2]) > 0.0 || error("edge interval is empty")
    abs(dot(edge.gap, edge.tangent)) <= 1.0e-10 || error("edge frame is not orthogonal")

    first, second = extended_interval(edge, radius)
    smin = first - radius
    smax = second + radius
    umin = -2radius
    umax = 2radius
    wmin = -radius - overetch
    wmax = radius + metal_thickness
    process_core_width = process_core_width > 0.0 ?
                         process_core_width :
                         max(2metal_thickness, 4overetch, 8lc_normal)
    minimum_process_dimension = min(
        metal_thickness,
        overetch > 0.0 ? overetch : metal_thickness
    )
    if fabricated && lc_normal > minimum_process_dimension
        @warn(
            "Normal spacing exceeds at least one fabrication dimension",
            lc_normal,
            metal_thickness,
            overetch
        )
    end

    longitudinal_intervals = max(1, ceil(Int, (smax - smin) / lc_tangent))
    longitudinal_coordinates = collect(range(smin, smax, length=longitudinal_intervals + 1))

    initialized = false
    try
        gmsh.initialize()
        initialized = true
        gmsh.option.setNumber("General.Verbosity", 2)
        gmsh.model.add("spatial_coupon_swept_$(fabricated ? "fabricated" : "thin")")
        geo = gmsh.model.geo
        builder = CrossSectionBuilder(
            geo,
            edge,
            smin,
            process_core_width,
            Dict{NTuple{2, Float64}, Int32}(),
            Dict{Int32, NTuple{2, Float64}}(),
            Dict{Tuple{Int32, Int32}, Int32}(),
            Dict{Tuple{Int32, Int32}, Tuple{Int32, Int32}}(),
            Dict{Tuple{NTuple{2, Float64}, NTuple{2, Float64}}, NTuple{2, Float64}}(),
            Tuple{Int32, Symbol}[],
            Float64[]
        )

        if fabricated
            build_fabricated_cross_section!(
                builder,
                umin,
                umax,
                wmin,
                wmax,
                metal_thickness,
                overetch,
                sidewall_angle,
                top_rounding,
                trench_rounding,
                lc_normal,
                lc_far,
                process_grading_power
            )
        else
            build_thin_cross_section!(
                builder,
                umin,
                umax,
                wmin,
                wmax,
                lc_normal,
                lc_far,
                process_grading_power
            )
        end
        isempty(builder.surfaces) && error("Swept cross-section produced no surfaces")

        input_surfaces = [(2, surface) for (surface, _) in builder.surfaces]
        delta_s = smax - smin
        output = geo.extrude(
            input_surfaces,
            delta_s * edge.tangent[1],
            delta_s * edge.tangent[2],
            0.0,
            [longitudinal_intervals],
            [1.0],
            true
        )
        geo.synchronize()

        substrate_tags = Int32[]
        vacuum_tags = Int32[]
        position = 1
        for (_, material) in builder.surfaces
            position + 1 <= length(output) || error("Unable to decode swept volume output")
            output[position][1] == 2 || error("Swept output is missing its top surface")
            output[position + 1][1] == 3 || error("Swept output is missing its volume")
            volume = output[position + 1][2]
            push!(material == :substrate ? substrate_tags : vacuum_tags, volume)
            position += 5 # top, volume, and the three lateral surfaces of a triangle
        end
        position == length(output) + 1 || error("Unexpected swept extrusion output layout")
        unique!(substrate_tags)
        unique!(vacuum_tags)

        for (name, value) in [
            ("Mesh.MeshSizeMin", lc_normal),
            ("Mesh.MeshSizeMax", lc_far),
            ("Mesh.MeshSizeExtendFromBoundary", 0),
            ("Mesh.MeshSizeFromPoints", 0),
            ("Mesh.MeshSizeFromCurvature", 0),
            ("Mesh.MshFileVersion", 2.2),
            ("Mesh.Binary", 1)
        ]
            gmsh.option.setNumber(name, value)
        end
        gmsh.model.mesh.generate(3)
        gmsh.model.mesh.setOrder(mesh_order)
        mesh_order > 1 && gmsh.model.mesh.optimize("HighOrder", true, 20)

        substrate_set = Set(substrate_tags)
        vacuum_set = Set(vacuum_tags)
        matching = Int32[]
        boundary_groups = Dict{Int, Vector{Int32}}()
        tolerance = 1.0e-8 * radius
        for (dimension, tag) in gmsh.model.getEntities(2)
            up, _ = gmsh.model.getAdjacencies(dimension, tag)
            adjacent_substrate = [volume for volume in up if volume in substrate_set]
            adjacent_vacuum = [volume for volume in up if volume in vacuum_set]
            isempty(adjacent_substrate) && isempty(adjacent_vacuum) && continue
            points = surface_local_nodes(edge, tag)
            if on_outer_surface(points, umin, umax, smin, smax, wmin, wmax, tolerance)
                push!(matching, tag)
                continue
            end
            attribute = 0
            if fabricated
                if !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                    w_values = [point[3] for point in points]
                    attribute = all(abs(value) <= tolerance for value in w_values) ?
                                3000 + edge.slot : 3100 + edge.slot
                elseif length(unique(up)) == 1 && !isempty(adjacent_substrate)
                    attribute = metal_surface_attribute(5000, edge.slot, edge.conductor)
                elseif length(unique(up)) == 1 && !isempty(adjacent_vacuum)
                    attribute = metal_surface_attribute(6000, edge.slot, edge.conductor)
                end
            elseif !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                u = sum(point[1] for point in points) / length(points)
                attribute = u <= tolerance ?
                            metal_surface_attribute(4000, edge.slot, edge.conductor) :
                            3000 + edge.slot
            end
            attribute > 0 && push!(get!(boundary_groups, attribute, Int32[]), tag)
        end
        isempty(matching) && error("No matching surface was generated")

        gmsh.model.addPhysicalGroup(3, sort!(unique(substrate_tags)), 1, "substrate")
        gmsh.model.addPhysicalGroup(3, sort!(unique(vacuum_tags)), 2, "vacuum")
        gmsh.model.addPhysicalGroup(2, sort!(unique(matching)), 1, "matching_surface")
        for (attribute, surfaces) in sort(collect(boundary_groups))
            gmsh.model.addPhysicalGroup(
                2,
                sort!(unique(surfaces)),
                attribute,
                "surface_$attribute"
            )
        end

        node_tags, _, _ = gmsh.model.mesh.getNodes()
        _, volume_element_tags, _ = gmsh.model.mesh.getElements(3)
        node_count = length(node_tags)
        element_count = sum(length(tags) for tags in volume_element_tags)
        node_count <= max_nodes ||
            error("Swept spatial coupon exceeds node budget: $node_count > $max_nodes")
        element_count <= max_elements || error(
            "Swept spatial coupon exceeds element budget: $element_count > $max_elements"
        )
        element_tags = reduce(vcat, volume_element_tags; init=UInt64[])
        isempty(element_tags) && error("Swept spatial coupon contains no volume elements")
        scaled_jacobians = gmsh.model.mesh.getElementQualities(element_tags, "minSJ")
        all(isfinite, scaled_jacobians) || error("Swept mesh contains a nonfinite Jacobian")
        minimum_jacobian = minimum(scaled_jacobians)
        minimum_jacobian > 0.0 || error("Swept mesh contains a nonpositive Jacobian")
        face_counts = prism_face_counts()

        normal_spacings = [
            value for value in builder.near_feature_spacings if value > tolerance
        ]
        isempty(normal_spacings) && error("No near-feature normal spacing was measured")
        maximum(normal_spacings) <= 1.05 * lc_normal || error(
            "Measured normal spacing $(maximum(normal_spacings)) exceeds " *
            "requested value $lc_normal"
        )
        maximum(diff(longitudinal_coordinates)) <= 1.05 * lc_tangent || error(
            "Measured longitudinal spacing exceeds requested value $lc_tangent"
        )

        mkpath(dirname(filename))
        gmsh.write(filename)

        # Reopen the serialized artifact and independently verify the evidence that is
        # persisted in metadata. This catches MSH-format filtering or physical-group
        # mistakes that are invisible in the in-memory model.
        gmsh.clear()
        gmsh.open(filename)
        written_node_tags, _, _ = gmsh.model.mesh.getNodes()
        written_types, written_element_tags, _ = gmsh.model.mesh.getElements(3)
        written_node_count = length(written_node_tags)
        written_element_count = sum(length(tags) for tags in written_element_tags)
        # Circle-center construction points are CAD helpers, not mesh nodes attached to
        # physical elements, and MSH 2.2 intentionally omits them. The serialized count
        # is authoritative for budgets and metadata.
        written_node_count <= node_count || error(
            "Written swept mesh gained nodes: $written_node_count > $node_count"
        )
        written_node_count <= max_nodes || error(
            "Written swept spatial coupon exceeds node budget: " *
            "$written_node_count > $max_nodes"
        )
        written_element_count == element_count || error(
            "Written swept mesh element count changed: " *
            "$written_element_count != $element_count"
        )
        written_tags = reduce(vcat, written_element_tags; init=UInt64[])
        written_jacobians = gmsh.model.mesh.getElementQualities(written_tags, "minSJ")
        written_minimum_jacobian = minimum(written_jacobians)
        written_minimum_jacobian > 0.0 ||
            error("Serialized swept mesh contains a nonpositive Jacobian")
        isapprox(written_minimum_jacobian, minimum_jacobian; rtol=1.0e-12, atol=1.0e-14) ||
            error("Serialized swept mesh Jacobian does not match in-memory validation")
        written_face_counts = prism_face_counts()
        written_face_counts == face_counts ||
            error("Serialized swept mesh face topology changed after writing")
        surface_attributes = sort!([
            Int(tag) for (dimension, tag) in gmsh.model.getPhysicalGroups(2)
        ])
        volume_attributes = sort!([
            Int(tag) for (dimension, tag) in gmsh.model.getPhysicalGroups(3)
        ])
        volume_element_types = sort!([
            gmsh.model.mesh.getElementProperties(element_type)[1] for
            element_type in written_types
        ])

        write_swept_metadata(
            filename * ".metadata.json";
            node_count=written_node_count,
            element_count=written_element_count,
            surface_attributes=surface_attributes,
            volume_attributes=volume_attributes,
            volume_element_types=volume_element_types,
            edge_count=length(edges),
            lc_normal=lc_normal,
            lc_tangent=lc_tangent,
            lc_far=lc_far,
            core_width=process_core_width,
            grading_power=process_grading_power,
            longitudinal_coordinates=longitudinal_coordinates,
            normal_spacings=normal_spacings,
            minimum_jacobian=written_minimum_jacobian,
            mesh_order=mesh_order,
            face_counts=written_face_counts,
            cross_section_nodes=length(builder.point_coordinates),
            fabricated=fabricated
        )
        println(
            "Swept spatial coupon: fabricated=$fabricated, prisms=$element_count, " *
            "nodes=$node_count, minSJ=$minimum_jacobian, file=$filename"
        )
    catch
        isfile(filename) && rm(filename; force=true)
        isfile(filename * ".metadata.json") && rm(filename * ".metadata.json"; force=true)
        rethrow()
    finally
        initialized && gmsh.finalize()
    end
    return
end

function parse_swept_options(args)
    length(args) >= 3 || error(
        "Usage: mesh_spatial_coupon_swept.jl SIGNATURE.csv thin|fabricated " *
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
        "--lc-normal" => ("lc_normal", Float64),
        "--lc-fine" => ("lc_normal", Float64),
        "--lc-tangent" => ("lc_tangent", Float64),
        "--lc-far" => ("lc_far", Float64),
        "--process-core-width" => ("process_core_width", Float64),
        "--process-fine-width" => ("ignored_process_fine_width", Float64),
        "--process-grading-power" => ("process_grading_power", Float64),
        "--normal-growth-ratio" => ("normal_growth_ratio", Float64),
        "--max-nodes" => ("max_nodes", Int),
        "--max-elements" => ("max_elements", Int),
        "--mesh-order" => ("mesh_order", Int)
    )
    seen_normal = nothing
    index = 4
    while index <= length(args)
        flag = args[index]
        haskey(names, flag) || error("Unknown option: $flag")
        index < length(args) || error("Missing value for option: $flag")
        name, type = names[flag]
        value = type === String ? abspath(args[index + 1]) : parse(type, args[index + 1])
        if flag in ("--lc-normal", "--lc-fine")
            seen_normal === nothing || seen_normal == value ||
                error("--lc-normal and --lc-fine specify different values")
            seen_normal = value
        end
        options[name] = value
        index += 2
    end
    get(options, "ignored_process_fine_width", 0.0) == 0.0 || error(
        "--process-fine-width is not part of the explicit swept-coordinate model"
    )
    return options
end

include(joinpath(@__DIR__, "mesh_spatial_coupon_edge_cluster.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_swept_options(ARGS)
    edge_count = length(read_edges(options["signature"]))
    common = (
        signature=options["signature"],
        fabricated=options["fabricated"],
        filename=options["filename"],
        radius=get(options, "radius", 2.0),
        metal_thickness=get(options, "metal_thickness", 0.1),
        overetch=get(options, "overetch", 0.05),
        sidewall_angle=get(options, "sidewall_angle", 80.0),
        top_rounding=get(options, "top_rounding", 0.01),
        trench_rounding=get(options, "trench_rounding", 0.01),
        lc_normal=get(options, "lc_normal", 0.02),
        lc_tangent=get(options, "lc_tangent", 0.2),
        lc_far=get(options, "lc_far", 0.3),
        process_core_width=get(options, "process_core_width", 0.0),
        process_grading_power=get(options, "process_grading_power", 1.7),
        max_nodes=get(options, "max_nodes", 500_000),
        max_elements=get(options, "max_elements", 2_000_000),
        mesh_order=get(options, "mesh_order", 2)
    )
    mask = get(options, "mask", nothing)
    boundary = get(options, "boundary", nothing)
    if mask !== nothing || boundary !== nothing
        mask === nothing && error("Masked edge-cluster meshing requires --mask")
        boundary === nothing && error("Masked edge-cluster meshing requires --boundary")
        generate_masked_edge_cluster_coupon(;
            common...,
            mask=mask,
            boundary=boundary,
            normal_growth_ratio=get(options, "normal_growth_ratio", 1.4)
        )
    elseif edge_count == 1
        generate_spatial_coupon_swept(;
            common...,
            mask=nothing,
            boundary=nothing
        )
    else
        error("Spatial edge clusters with multiple edges require --mask and --boundary")
    end
end
