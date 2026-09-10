# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Masked spatial edge-cluster transition mesher. This file is included by
# mesh_spatial_coupon_swept.jl after its common parsing and metadata helpers.

function geometric_distances(width, first_spacing, ratio)
    width >= 0.0 || error("Geometric-coordinate width must be nonnegative")
    first_spacing > 0.0 || error("Geometric first spacing must be positive")
    ratio > 1.0 || error("Geometric growth ratio must exceed one")
    width == 0.0 && return [0.0]
    width <= first_spacing && return [0.0, width]
    result = [0.0]
    distance = 0.0
    spacing = first_spacing
    while distance + spacing < width
        distance += spacing
        push!(result, distance)
        spacing *= ratio
    end
    remainder = width - result[end]
    if length(result) > 1 && remainder < 0.45 * (result[end] - result[end - 1])
        previous = result[end - 1]
        total = width - previous
        result[end] = previous + 0.5 * total
    end
    push!(result, width)
    return unique(result)
end

function geometric_interval(first, last, spacing, ratio; refine=:left)
    last > first || error("Geometric interval must have positive width")
    distances = geometric_distances(last - first, spacing, ratio)
    if refine == :left
        return first .+ distances
    elseif refine == :right
        return last .- reverse(distances)
    elseif refine == :both
        middle = 0.5 * (first + last)
        left = geometric_interval(first, middle, spacing, ratio; refine=:left)
        right = geometric_interval(middle, last, spacing, ratio; refine=:right)
        return sort!(unique(vcat(left, right)))
    end
    error("Unknown geometric refinement side: $refine")
end

function merge_coordinate_blocks(blocks)
    values = Float64[]
    for block in blocks
        for value in block
            (isempty(values) || !isapprox(value, values[end]; atol=1.0e-13, rtol=0.0)) &&
                push!(values, value)
        end
    end
    return values
end

function clip_segment_to_box(first, second, lower, upper, tolerance)
    x0, y0 = first
    dx = second[1] - x0
    dy = second[2] - y0
    begin_parameter = 0.0
    end_parameter = 1.0
    for (p, q) in (
        (-dx, x0 - lower[1]),
        (dx, upper[1] - x0),
        (-dy, y0 - lower[2]),
        (dy, upper[2] - y0)
    )
        if abs(p) <= tolerance
            q >= -tolerance || return nothing
            continue
        end
        ratio = q / p
        if p < 0.0
            begin_parameter = max(begin_parameter, ratio)
        else
            end_parameter = min(end_parameter, ratio)
        end
        begin_parameter <= end_parameter + tolerance || return nothing
    end
    clipped_first = (x0 + begin_parameter * dx, y0 + begin_parameter * dy)
    clipped_second = (x0 + end_parameter * dx, y0 + end_parameter * dy)
    hypot(
        clipped_second[1] - clipped_first[1],
        clipped_second[2] - clipped_first[2]
    ) > tolerance || return nothing
    return clipped_first, clipped_second
end

function segment_key(first, second, tolerance)
    quantize(point) = (
        round(Int, point[1] / tolerance),
        round(Int, point[2] / tolerance)
    )
    a = quantize(first)
    b = quantize(second)
    return a < b ? (a, b) : (b, a)
end

function point_segment_projection(point, first, second)
    direction = (second[1] - first[1], second[2] - first[2])
    length_squared = direction[1]^2 + direction[2]^2
    length_squared > 0.0 || return (
        distance=hypot(point[1] - first[1], point[2] - first[2]),
        coordinate=0.0
    )
    coordinate = (
        (point[1] - first[1]) * direction[1] +
        (point[2] - first[2]) * direction[2]
    ) / length_squared
    closest = (
        first[1] + clamp(coordinate, 0.0, 1.0) * direction[1],
        first[2] + clamp(coordinate, 0.0, 1.0) * direction[2]
    )
    return (
        distance=hypot(point[1] - closest[1], point[2] - closest[2]),
        coordinate=coordinate
    )
end

function active_plan_segments(edges, radius)
    result = NamedTuple[]
    for edge in edges
        first, second = extended_interval(edge, radius)
        p0 = add(edge.point, scale(first, edge.tangent))
        p1 = add(edge.point, scale(second, edge.tangent))
        push!(
            result,
            (
                first=(p0[1], p0[2]),
                second=(p1[1], p1[2]),
                gap=(edge.gap[1], edge.gap[2]),
                edge=edge
            )
        )
    end
    return result
end

function nearest_active_segment(segments, point)
    index = argmin(
        point_segment_projection(point, segment.first, segment.second).distance for
        segment in segments
    )
    return segments[index]
end

function plan_point_in_mask(facets, point, conductor, plane, tolerance)
    return point_in_mask(
        facets,
        (point[1], point[2], plane),
        conductor,
        plane,
        tolerance
    )
end

function plan_conductor(facets, point, plane, tolerance)
    conductors = [
        conductor for conductor in sort!(unique(facet.conductor for facet in facets)) if
        plan_point_in_mask(facets, point, conductor, plane, tolerance)
    ]
    length(conductors) <= 1 || error(
        "Plan-view masks for conductors $(join(conductors, ", ")) overlap"
    )
    return isempty(conductors) ? 0 : only(conductors)
end

function classified_fabrication_primitives(loops, tolerance)
    primitives = NamedTuple[]
    for loop in loops
        for primitive in physical_segments([loop], 0.0, tolerance)
            push!(primitives, (; primitive..., conductor=loop.conductor))
        end
    end
    return primitives
end

function planar_mesh_with_transition_rows(
    edges,
    facets,
    loops,
    lower,
    upper,
    radius,
    lc_normal,
    lc_tangent,
    lc_far,
    core_width,
    growth_ratio,
    model_name
)
    tolerance = 1.0e-9 * radius
    gmsh.model.add(model_name)
    occ = gmsh.model.occ
    width = upper[1] - lower[1]
    height = upper[2] - lower[2]
    segments = active_plan_segments(edges, radius)
    distances = geometric_distances(core_width, lc_normal, growth_ratio)

    raw_segments = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for loop in loops
        for index in eachindex(loop.points)
            push!(
                raw_segments,
                (loop.points[index], loop.points[mod1(index + 1, length(loop.points))])
            )
        end
    end
    for segment in segments
        for distance in distances
            for side in (distance == 0.0 ? (0.0,) : (-distance, distance))
                shift = (side * segment.gap[1], side * segment.gap[2])
                first = (segment.first[1] + shift[1], segment.first[2] + shift[2])
                second = (segment.second[1] + shift[1], segment.second[2] + shift[2])
                clipped = clip_segment_to_box(
                    first,
                    second,
                    (lower[1], lower[2]),
                    (upper[1], upper[2]),
                    tolerance
                )
                clipped === nothing || push!(raw_segments, clipped)
            end
        end
    end

    tags = Tuple{Int32, Int32}[]
    seen = Set{Any}()
    for (first, second) in raw_segments
        key = segment_key(first, second, tolerance)
        key in seen && continue
        push!(seen, key)
        p0 = occ.addPoint(first[1], first[2], edges[1].point[3])
        p1 = occ.addPoint(second[1], second[2], edges[1].point[3])
        push!(tags, (1, occ.addLine(p0, p1)))
    end
    outer = (2, occ.addRectangle(lower[1], lower[2], edges[1].point[3], width, height))
    partition, _ = occ.fragment([outer], tags)
    plan_surfaces = unique(entity[2] for entity in partition if entity[1] == 2)
    isempty(plan_surfaces) && error("Edge-cluster constraints removed the plan surface")
    occ.synchronize()

    constrained_curves = Int32[]
    tangent_spacings = Float64[]
    active_tangent_spacings = Float64[]
    for (_, curve) in gmsh.model.getEntities(1)
        bounds = gmsh.model.getBoundingBox(1, curve)
        midpoint = (0.5 * (bounds[1] + bounds[4]), 0.5 * (bounds[2] + bounds[5]))
        curve_length = gmsh.model.occ.getMass(1, curve)
        near_active = any(
            begin
                projection = point_segment_projection(
                    midpoint,
                    segment.first,
                    segment.second
                )
                projection.distance <= core_width + tolerance &&
                    -tolerance <= projection.coordinate <= 1.0 + tolerance
            end for segment in segments
        )
        requested = near_active ? lc_tangent : lc_far
        point_count = max(2, ceil(Int, curve_length / requested) + 1)
        gmsh.model.mesh.setTransfiniteCurve(curve, point_count)
        spacing = curve_length / (point_count - 1)
        push!(constrained_curves, curve)
        push!(tangent_spacings, spacing)
        near_active && push!(active_tangent_spacings, spacing)
    end

    for (name, value) in [
        ("Mesh.MeshSizeMin", lc_normal),
        ("Mesh.MeshSizeMax", lc_far),
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.Algorithm", 6)
    ]
        gmsh.option.setNumber(name, value)
    end
    gmsh.model.mesh.generate(2)

    node_tags, coordinates, _ = gmsh.model.mesh.getNodes(-1, -1, false, false)
    coordinate_map = Dict{UInt64, NTuple{2, Float64}}()
    for (index, tag) in enumerate(node_tags)
        coordinate_map[tag] = (
            coordinates[3index - 2],
            coordinates[3index - 1]
        )
    end
    types, element_tags, connectivity = gmsh.model.mesh.getElements(2)
    triangles = NTuple{3, UInt64}[]
    for (element_type, tags_for_type, nodes_for_type) in
        zip(types, element_tags, connectivity)
        name, _, _, node_count, _, primary_count =
            gmsh.model.mesh.getElementProperties(element_type)
        name == "Triangle 3" ||
            error("Edge-cluster plan mesh generated unsupported element $name")
        node_count == 3 && primary_count == 3 ||
            error("Edge-cluster plan mesh must be linear before extrusion")
        for offset = 1:3:length(nodes_for_type)
            push!(triangles, Tuple(nodes_for_type[offset:(offset + 2)]))
        end
    end
    isempty(triangles) && error("Edge-cluster transition plan contains no triangles")
    plan_edges = Dict{Tuple{UInt64, UInt64}, Int}()
    for triangle in triangles
        for (first, second) in (
            (triangle[1], triangle[2]),
            (triangle[2], triangle[3]),
            (triangle[3], triangle[1])
        )
            key = first < second ? (first, second) : (second, first)
            plan_edges[key] = get(plan_edges, key, 0) + 1
        end
    end
    count(>(2), values(plan_edges)) == 0 ||
        error("Edge-cluster transition plan contains a nonmanifold edge")
    unmatched_internal_edges = 0
    for ((first, second), count) in plan_edges
        count == 1 || continue
        a = coordinate_map[first]
        b = coordinate_map[second]
        on_outer = (
            (abs(a[1] - lower[1]) <= tolerance && abs(b[1] - lower[1]) <= tolerance) ||
            (abs(a[1] - upper[1]) <= tolerance && abs(b[1] - upper[1]) <= tolerance) ||
            (abs(a[2] - lower[2]) <= tolerance && abs(b[2] - lower[2]) <= tolerance) ||
            (abs(a[2] - upper[2]) <= tolerance && abs(b[2] - upper[2]) <= tolerance)
        )
        on_outer || (unmatched_internal_edges += 1)
    end
    unmatched_internal_edges == 0 || error(
        "Edge-cluster transition plan contains $unmatched_internal_edges unmatched " *
        "internal edges"
    )
    return (
        node_tags=sort!(collect(keys(coordinate_map))),
        coordinates=coordinate_map,
        triangles=triangles,
        active_segments=segments,
        fabrication_primitives=classified_fabrication_primitives(loops, tolerance),
        normal_distances=distances,
        maximum_tangent_spacing=isempty(active_tangent_spacings) ?
                                0.0 :
                                maximum(active_tangent_spacings),
        maximum_far_spacing=isempty(tangent_spacings) ? 0.0 : maximum(tangent_spacings),
        embedded_curve_count=length(constrained_curves),
        plan_surface_count=length(plan_surfaces),
        unmatched_internal_edge_count=unmatched_internal_edges
    )
end

function edge_cluster_z_coordinates(
    plane,
    lower,
    upper,
    fabricated,
    metal_thickness,
    overetch,
    lc_normal,
    growth_ratio
)
    if fabricated
        return merge_coordinate_blocks([
            geometric_interval(
                lower,
                plane - overetch,
                lc_normal,
                growth_ratio;
                refine=:right
            ),
            geometric_interval(
                plane - overetch,
                plane,
                lc_normal,
                growth_ratio;
                refine=:both
            ),
            geometric_interval(
                plane,
                plane + metal_thickness,
                lc_normal,
                growth_ratio;
                refine=:both
            ),
            geometric_interval(
                plane + metal_thickness,
                upper,
                lc_normal,
                growth_ratio;
                refine=:left
            )
        ])
    end
    return merge_coordinate_blocks([
        geometric_interval(lower, plane, lc_normal, growth_ratio; refine=:right),
        geometric_interval(plane, upper, lc_normal, growth_ratio; refine=:left)
    ])
end

function triangle_centroid(triangle, coordinates)
    points = [coordinates[tag] for tag in triangle]
    return (
        sum(point[1] for point in points) / 3.0,
        sum(point[2] for point in points) / 3.0
    )
end

function triangle_phase(
    triangle,
    coordinates,
    facets,
    fabrication_primitives,
    plane,
    radius,
    tolerance;
    etch_full_gap=false
)
    point = triangle_centroid(triangle, coordinates)
    conductor = plan_conductor(facets, point, plane, tolerance)
    conductor > 0 && return (:metal, conductor)
    distances = [
        point_primitive_distance(point, primitive, tolerance) for
        primitive in fabrication_primitives
    ]
    nearest = argmin(distances)
    return etch_full_gap || distances[nearest] <= 3radius + tolerance ?
           (:trench, fabrication_primitives[nearest].conductor) :
           (:ordinary, 0)
end

function interval_material(fabricated, phase, zmid, plane, metal_thickness, overetch)
    if !fabricated
        return zmid < plane ? :substrate : :vacuum
    elseif zmid < plane - overetch
        return :substrate
    elseif zmid < plane
        return phase == :trench ? :vacuum : :substrate
    elseif zmid < plane + metal_thickness
        return phase == :metal ? :dummy : :vacuum
    end
    return :vacuum
end

function add_surface_elements!(surface_data, attribute, kind, connectivity)
    entry = get!(surface_data, attribute) do
        Dict(:triangles => UInt64[], :quads => UInt64[])
    end
    append!(entry[kind], connectivity)
    return
end

function edge_cluster_owner(edges, segments, point, conductor=0)
    candidates = [
        (edge=edge, segment=segment) for (edge, segment) in zip(edges, segments) if
        conductor == 0 || edge.conductor == conductor
    ]
    isempty(candidates) && error("Unable to assign edge-cluster surface owner")
    index = argmin(
        point_segment_projection(
            point,
            candidate.segment.first,
            candidate.segment.second
        ).distance for candidate in candidates
    )
    return candidates[index].edge
end

function horizontal_attribute(
    fabricated,
    z,
    phase,
    conductor,
    edges,
    segments,
    point,
    plane,
    metal_thickness,
    overetch,
    tolerance
)
    owner = edge_cluster_owner(edges, segments, point, conductor)
    if !fabricated && abs(z - plane) <= tolerance
        return phase == :metal ?
               metal_surface_attribute(4000, owner.slot, owner.conductor) :
               3000 + owner.slot
    elseif fabricated && abs(z - (plane - overetch)) <= tolerance && phase == :trench
        return 3100 + owner.slot
    elseif fabricated && abs(z - plane) <= tolerance
        if phase == :metal
            return metal_surface_attribute(5000, owner.slot, owner.conductor)
        elseif phase == :ordinary
            return 3000 + owner.slot
        end
    elseif fabricated && abs(z - (plane + metal_thickness)) <= tolerance && phase == :metal
        return metal_surface_attribute(6000, owner.slot, owner.conductor)
    end
    return 0
end

function vertical_attribute(
    fabricated,
    zmid,
    first_phase,
    second_phase,
    first_conductor,
    second_conductor,
    edges,
    segments,
    point,
    plane,
    metal_thickness,
    overetch
)
    fabricated || return 0
    conductor = first_conductor > 0 ? first_conductor : second_conductor
    owner = edge_cluster_owner(edges, segments, point, conductor)
    phases = Set((first_phase, second_phase))
    if plane - overetch < zmid < plane && :trench in phases
        return 3100 + owner.slot
    elseif plane < zmid < plane + metal_thickness && :metal in phases
        return metal_surface_attribute(6000, owner.slot, owner.conductor)
    end
    return 0
end

function add_discrete_edge_cluster_mesh!(
    plan,
    edges,
    facets,
    fabricated,
    z_coordinates,
    lower,
    upper,
    radius,
    metal_thickness,
    overetch,
    max_nodes,
    max_elements,
    mesh_order;
    etch_full_gap=false
)
    tolerance = 1.0e-8 * radius
    plane = edges[1].point[3]
    plan_tags = plan.node_tags
    node_index = Dict(tag => index for (index, tag) in enumerate(plan_tags))
    layer_count = length(z_coordinates)
    plan_node_count = length(plan_tags)
    node_tags = UInt64[]
    coordinates = Float64[]
    for (layer, z) in enumerate(z_coordinates)
        for tag in plan_tags
            point = plan.coordinates[tag]
            push!(node_tags, UInt64((layer - 1) * plan_node_count + node_index[tag]))
            append!(coordinates, (point[1], point[2], z))
        end
    end
    global_node(layer, tag) = UInt64((layer - 1) * plan_node_count + node_index[tag])

    phases = [
        triangle_phase(
            triangle,
            plan.coordinates,
            facets,
            plan.fabrication_primitives,
            plane,
            radius,
            tolerance;
            etch_full_gap=etch_full_gap
        ) for triangle in plan.triangles
    ]
    substrate_connectivity = UInt64[]
    vacuum_connectivity = UInt64[]
    substrate_elements = UInt64[]
    vacuum_elements = UInt64[]
    element_tag = UInt64(1)
    for layer = 1:(layer_count - 1)
        zmid = 0.5 * (z_coordinates[layer] + z_coordinates[layer + 1])
        for (triangle, phase_record) in zip(plan.triangles, phases)
            material = interval_material(
                fabricated,
                phase_record[1],
                zmid,
                plane,
                metal_thickness,
                overetch
            )
            material == :dummy && continue
            prism = UInt64[
                global_node(layer, triangle[1]),
                global_node(layer, triangle[2]),
                global_node(layer, triangle[3]),
                global_node(layer + 1, triangle[1]),
                global_node(layer + 1, triangle[2]),
                global_node(layer + 1, triangle[3])
            ]
            if material == :substrate
                push!(substrate_elements, element_tag)
                append!(substrate_connectivity, prism)
            else
                push!(vacuum_elements, element_tag)
                append!(vacuum_connectivity, prism)
            end
            element_tag += 1
        end
    end
    linear_element_count = length(substrate_elements) + length(vacuum_elements)
    linear_element_count <= max_elements || error(
        "Edge-cluster coupon exceeds linear element budget: " *
        "$linear_element_count > $max_elements"
    )

    # Do not register nodes inside removed metal cells. They belong to no field
    # element and Gmsh drops them on serialization, which otherwise changes the
    # linear-mesh node count (and wastes memory before writing).
    used_nodes=Set(substrate_connectivity)
    union!(used_nodes,vacuum_connectivity)
    keep=findall(tag->tag in used_nodes,node_tags)
    node_tags=node_tags[keep]
    used_coordinates=Float64[]
    sizehint!(used_coordinates,3length(keep))
    for i in keep, d in 1:3
        push!(used_coordinates,coordinates[3i-3+d])
    end
    coordinates=used_coordinates

    substrate_entity = gmsh.model.addDiscreteEntity(3, 1)
    vacuum_entity = gmsh.model.addDiscreteEntity(3, 2)
    gmsh.model.mesh.addNodes(3, substrate_entity, node_tags, coordinates)
    !isempty(substrate_elements) && gmsh.model.mesh.addElementsByType(
        substrate_entity,
        6,
        substrate_elements,
        substrate_connectivity
    )
    !isempty(vacuum_elements) && gmsh.model.mesh.addElementsByType(
        vacuum_entity,
        6,
        vacuum_elements,
        vacuum_connectivity
    )

    edge_neighbors = Dict{Tuple{UInt64, UInt64}, Vector{Int}}()
    for (index, triangle) in enumerate(plan.triangles)
        for (a, b) in (
            (triangle[1], triangle[2]),
            (triangle[2], triangle[3]),
            (triangle[3], triangle[1])
        )
            key = a < b ? (a, b) : (b, a)
            push!(get!(edge_neighbors, key, Int[]), index)
        end
    end

    surface_data = Dict{Int, Dict{Symbol, Vector{UInt64}}}()
    # Matching planes at the lower and upper coupon bounds.
    for (layer, z) in ((1, z_coordinates[1]), (layer_count, z_coordinates[end]))
        for (triangle, phase_record) in zip(plan.triangles, phases)
            adjacent_layer = layer == 1 ? 1 : layer_count - 1
            zmid = 0.5 * (
                z_coordinates[adjacent_layer] + z_coordinates[adjacent_layer + 1]
            )
            material = interval_material(
                fabricated,
                phase_record[1],
                zmid,
                plane,
                metal_thickness,
                overetch
            )
            material == :dummy && continue
            add_surface_elements!(
                surface_data,
                1,
                :triangles,
                UInt64[
                    global_node(layer, triangle[1]),
                    global_node(layer, triangle[2]),
                    global_node(layer, triangle[3])
                ]
            )
        end
    end

    # Horizontal process interfaces.
    for layer = 2:(layer_count - 1)
        z = z_coordinates[layer]
        for (triangle, phase_record) in zip(plan.triangles, phases)
            below = interval_material(
                fabricated,
                phase_record[1],
                0.5 * (z_coordinates[layer - 1] + z),
                plane,
                metal_thickness,
                overetch
            )
            above = interval_material(
                fabricated,
                phase_record[1],
                0.5 * (z + z_coordinates[layer + 1]),
                plane,
                metal_thickness,
                overetch
            )
            below == above && continue
            point = triangle_centroid(triangle, plan.coordinates)
            attribute = horizontal_attribute(
                fabricated,
                z,
                phase_record[1],
                phase_record[2],
                edges,
                plan.active_segments,
                point,
                plane,
                metal_thickness,
                overetch,
                tolerance
            )
            attribute > 0 || continue
            add_surface_elements!(
                surface_data,
                attribute,
                :triangles,
                UInt64[
                    global_node(layer, triangle[1]),
                    global_node(layer, triangle[2]),
                    global_node(layer, triangle[3])
                ]
            )
        end
    end

    # Vertical outer matching boundary and material-transition walls.
    for (edge_key, neighbors) in edge_neighbors
        first_tag, second_tag = edge_key
        midpoint = (
            0.5 * (plan.coordinates[first_tag][1] + plan.coordinates[second_tag][1]),
            0.5 * (plan.coordinates[first_tag][2] + plan.coordinates[second_tag][2])
        )
        boundary_edge = length(neighbors) == 1
        on_outer_box = (
            (abs(plan.coordinates[first_tag][1] - lower[1]) <= tolerance &&
             abs(plan.coordinates[second_tag][1] - lower[1]) <= tolerance) ||
            (abs(plan.coordinates[first_tag][1] - upper[1]) <= tolerance &&
             abs(plan.coordinates[second_tag][1] - upper[1]) <= tolerance) ||
            (abs(plan.coordinates[first_tag][2] - lower[2]) <= tolerance &&
             abs(plan.coordinates[second_tag][2] - lower[2]) <= tolerance) ||
            (abs(plan.coordinates[first_tag][2] - upper[2]) <= tolerance &&
             abs(plan.coordinates[second_tag][2] - upper[2]) <= tolerance)
        )
        for layer = 1:(layer_count - 1)
            zmid = 0.5 * (z_coordinates[layer] + z_coordinates[layer + 1])
            attribute = 0
            if boundary_edge && on_outer_box
                phase_record = phases[neighbors[1]]
                material = interval_material(
                    fabricated,
                    phase_record[1],
                    zmid,
                    plane,
                    metal_thickness,
                    overetch
                )
                material == :dummy || (attribute = 1)
            elseif length(neighbors) == 2 && phases[neighbors[1]] != phases[neighbors[2]]
                first_phase = phases[neighbors[1]]
                second_phase = phases[neighbors[2]]
                first_material = interval_material(
                    fabricated,
                    first_phase[1],
                    zmid,
                    plane,
                    metal_thickness,
                    overetch
                )
                second_material = interval_material(
                    fabricated,
                    second_phase[1],
                    zmid,
                    plane,
                    metal_thickness,
                    overetch
                )
                first_material == second_material && continue
                attribute = vertical_attribute(
                    fabricated,
                    zmid,
                    first_phase[1],
                    second_phase[1],
                    first_phase[2],
                    second_phase[2],
                    edges,
                    plan.active_segments,
                    midpoint,
                    plane,
                    metal_thickness,
                    overetch
                )
            end
            attribute > 0 || continue
            add_surface_elements!(
                surface_data,
                attribute,
                :quads,
                UInt64[
                    global_node(layer, first_tag),
                    global_node(layer, second_tag),
                    global_node(layer + 1, second_tag),
                    global_node(layer + 1, first_tag)
                ]
            )
        end
    end

    surface_entities = Dict{Int, Int32}()
    for (attribute, data) in sort(collect(surface_data); by=first)
        entity = gmsh.model.addDiscreteEntity(2)
        surface_entities[attribute] = entity
        if !isempty(data[:triangles])
            count = div(length(data[:triangles]), 3)
            tags = collect(element_tag:(element_tag + count - 1))
            gmsh.model.mesh.addElementsByType(entity, 2, tags, data[:triangles])
            element_tag += count
        end
        if !isempty(data[:quads])
            count = div(length(data[:quads]), 4)
            tags = collect(element_tag:(element_tag + count - 1))
            gmsh.model.mesh.addElementsByType(entity, 3, tags, data[:quads])
            element_tag += count
        end
    end

    gmsh.model.addPhysicalGroup(3, [substrate_entity], 1, "substrate")
    gmsh.model.addPhysicalGroup(3, [vacuum_entity], 2, "vacuum")
    for (attribute, entity) in sort(collect(surface_entities); by=first)
        name = attribute == 1 ? "matching_surface" : "surface_$attribute"
        gmsh.model.addPhysicalGroup(2, [entity], attribute, name)
    end

    gmsh.model.mesh.setOrder(mesh_order)
    node_count = length(gmsh.model.mesh.getNodes()[1])
    node_count <= max_nodes ||
        error("Edge-cluster coupon exceeds node budget: $node_count > $max_nodes")
    volume_types, volume_element_tags, _ = gmsh.model.mesh.getElements(3)
    volume_count = sum(length(tags) for tags in volume_element_tags)
    volume_count <= max_elements || error(
        "Edge-cluster coupon exceeds element budget: $volume_count > $max_elements"
    )
    all_tags = reduce(vcat, volume_element_tags; init=UInt64[])
    minimum_jacobian = minimum(gmsh.model.mesh.getElementQualities(all_tags, "minSJ"))
    minimum_jacobian > 0.0 || error("Edge-cluster coupon contains a nonpositive Jacobian")
    face_counts = prism_face_counts()
    return (
        node_count=node_count,
        volume_count=volume_count,
        minimum_jacobian=minimum_jacobian,
        surface_attributes=sort!(collect(keys(surface_entities))),
        volume_element_types=sort!([
            gmsh.model.mesh.getElementProperties(element_type)[1] for
            element_type in volume_types
        ]),
        plan_triangle_count=length(plan.triangles),
        z_station_count=layer_count,
        face_counts=face_counts
    )
end

function write_edge_cluster_metadata(
    path,
    evidence,
    plan,
    z_coordinates,
    lc_normal,
    lc_tangent,
    lc_far,
    core_width,
    growth_ratio,
    mesh_order,
    fabricated,
    edge_count
)
    open(path, "w") do stream
        println(stream, "{")
        println(stream, "  \"Version\": 2,")
        println(stream, "  \"MeshingMode\": \"ExplicitEdgeClusterTransition\",")
        println(stream, "  \"SweptTopology\": \"PlanarRowsAndPrisms\",")
        println(stream, "  \"Fabricated\": $(fabricated ? "true" : "false"),")
        println(stream, "  \"MetalSurfacePartition\": \"InterfaceSlotAndConductor\",")
        println(stream, "  \"NodeCount\": $(evidence.node_count),")
        println(stream, "  \"VolumeElementCount\": $(evidence.volume_count),")
        println(stream, "  \"SweptVolumeElementCount\": $(evidence.volume_count),")
        println(stream, "  \"TransitionVolumeElementCount\": $(evidence.volume_count),")
        println(stream, "  \"InputEdgeCount\": $edge_count,")
        println(stream, "  \"EtchFullGap\": $(hasproperty(evidence,:etch_full_gap) && evidence.etch_full_gap),")
        println(stream, "  \"FineSize\": $lc_normal,")
        println(stream, "  \"NormalSize\": $lc_normal,")
        println(stream, "  \"TangentialSize\": $lc_tangent,")
        println(stream, "  \"FarSize\": $lc_far,")
        println(stream, "  \"ProcessCoreWidth\": $core_width,")
        println(stream, "  \"NormalGrowthRatio\": $growth_ratio,")
        println(stream, "  \"NormalCoordinates\": $(json_array(plan.normal_distances)),")
        println(
            stream,
            "  \"MeasuredFirstNormalSpacing\": $(plan.normal_distances[2]),"
        )
        println(stream, "  \"NormalCoordinateCount\": $(length(plan.normal_distances)),")
        println(
            stream,
            "  \"MeasuredMaximumLongitudinalSpacing\": " *
            "$(plan.maximum_tangent_spacing),"
        )
        println(
            stream,
            "  \"MeasuredMaximumFarSpacing\": $(plan.maximum_far_spacing),"
        )
        println(stream, "  \"FabricationNormalCoordinates\": $(json_array(z_coordinates)),")
        println(stream, "  \"PlanTriangleCount\": $(evidence.plan_triangle_count),")
        println(stream, "  \"EmbeddedCurveCount\": $(plan.embedded_curve_count),")
        println(
            stream,
            "  \"UnmatchedInternalPlanEdgeCount\": " *
            "$(plan.unmatched_internal_edge_count),"
        )
        println(stream, "  \"BoundaryFaceCount\": $(evidence.face_counts.boundary),")
        println(stream, "  \"InteriorFaceCount\": $(evidence.face_counts.interior),")
        println(stream, "  \"NonmanifoldFaceCount\": $(evidence.face_counts.nonmanifold),")
        println(
            stream,
            "  \"SurfaceAttributes\": " *
            "$(json_array(evidence.surface_attributes)),"
        )
        println(stream, "  \"VolumeAttributes\": [1, 2],")
        println(
            stream,
            "  \"VolumeElementTypes\": [" *
            join(("\"" * value * "\"" for value in evidence.volume_element_types), ", ") *
            "],"
        )
        println(stream, "  \"MinimumScaledJacobian\": $(evidence.minimum_jacobian),")
        println(stream, "  \"MeshOrder\": $mesh_order")
        println(stream, "}")
    end
    return
end

function generate_masked_edge_cluster_coupon(;
    signature::String,
    mask::String,
    boundary::String,
    fabricated::Bool,
    radius::Float64,
    metal_thickness::Float64,
    overetch::Float64,
    sidewall_angle::Float64,
    top_rounding::Float64,
    trench_rounding::Float64,
    lc_normal::Float64,
    lc_tangent::Float64,
    lc_far::Float64,
    process_core_width::Float64,
    process_grading_power::Float64,
    normal_growth_ratio::Float64,
    max_nodes::Int,
    max_elements::Int,
    mesh_order::Int,
    filename::String,
    etch_full_gap::Bool=false
)
    abs(sidewall_angle - 90.0) <= 1.0e-12 || error(
        "Masked edge-cluster transition meshing requires 90 degree sidewalls"
    )
    top_rounding == 0.0 || error(
        "Masked edge-cluster transition meshing requires zero top rounding"
    )
    trench_rounding == 0.0 || error(
        "Masked edge-cluster transition meshing requires zero trench rounding"
    )
    normal_growth_ratio > 1.0 || error("normal growth ratio must exceed one")
    edges = read_edges(signature)
    !isempty(edges) || error(
        "Masked edge-cluster transition meshing requires at least one edge"
    )
    facets = read_mask(mask)
    loops = read_boundary(boundary)
    isempty(facets) && error("Masked edge-cluster meshing requires mask facets")
    isempty(loops) && error("Masked edge-cluster meshing requires boundary loops")
    all(!edge.vertex_arm for edge in edges) ||
        error("Masked edge-cluster meshing requires continuing spatial edges")
    edge_conductors = Set(edge.conductor for edge in edges)
    facet_conductors = Set(facet.conductor for facet in facets)
    loop_conductors = Set(loop.conductor for loop in loops)
    edge_conductors == facet_conductors == loop_conductors || error(
        "Edge, mask, and classified-boundary conductor sets do not agree"
    )
    all(abs(edge.point[3] - edges[1].point[3]) <= 1.0e-10radius for edge in edges) ||
        error("Masked edge-cluster edges must share one fabrication plane")
    all(abs(facet.plane - edges[1].point[3]) <= 1.0e-10radius for facet in facets) ||
        error("Mask facets do not lie on the edge-cluster fabrication plane")
    all(abs(loop.plane - edges[1].point[3]) <= 1.0e-10radius for loop in loops) ||
        error("Classified boundaries do not lie on the fabrication plane")
    all(edge.normal_sign == edges[1].normal_sign for edge in edges) ||
        error("Masked edge-cluster edges must share one fabrication normal")
    edges[1].normal_sign > 0.0 || error(
        "Masked edge-cluster transition meshing requires positive process normal"
    )

    lower, upper = coupon_bounds(edges, radius, metal_thickness, overetch)
    core_width = process_core_width > 0.0 ?
                 process_core_width :
                 max(2metal_thickness, 4overetch, 8lc_normal)
    core_width >= 8lc_normal || error(
        "Edge-cluster collar width must contain at least eight first-layer spacings"
    )
    initialized = false
    try
        gmsh.initialize()
        initialized = true
        gmsh.option.setNumber("General.Verbosity", 2)
        plan = planar_mesh_with_transition_rows(
            edges,
            facets,
            loops,
            lower,
            upper,
            radius,
            lc_normal,
            lc_tangent,
            lc_far,
            core_width,
            normal_growth_ratio,
            "spatial_coupon_edge_cluster_plan"
        )
        plane = edges[1].point[3]
        z_coordinates = edge_cluster_z_coordinates(
            plane,
            lower[3],
            upper[3],
            fabricated,
            metal_thickness,
            overetch,
            lc_normal,
            normal_growth_ratio
        )

        # Replace the temporary OCC plan model with one discrete 3D model.
        gmsh.clear()
        gmsh.model.add(
            "spatial_coupon_edge_cluster_$(fabricated ? "fabricated" : "thin")"
        )
        evidence = add_discrete_edge_cluster_mesh!(
            plan,
            edges,
            facets,
            fabricated,
            z_coordinates,
            lower,
            upper,
            radius,
            metal_thickness,
            overetch,
            max_nodes,
            max_elements,
            mesh_order;
            etch_full_gap=etch_full_gap
        )
        mkpath(dirname(filename))
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(filename)

        gmsh.clear()
        gmsh.open(filename)
        written_nodes = length(gmsh.model.mesh.getNodes()[1])
        written_types, written_tags, _ = gmsh.model.mesh.getElements(3)
        written_elements = sum(length(tags) for tags in written_tags)
        written_nodes == evidence.node_count || error(
            "Serialized edge-cluster node count changed: " *
            "$written_nodes != $(evidence.node_count)"
        )
        written_elements == evidence.volume_count || error(
            "Serialized edge-cluster element count changed: " *
            "$written_elements != $(evidence.volume_count)"
        )
        all_written_tags = reduce(vcat, written_tags; init=UInt64[])
        written_jacobian = minimum(
            gmsh.model.mesh.getElementQualities(all_written_tags, "minSJ")
        )
        written_jacobian > 0.0 ||
            error("Serialized edge-cluster mesh contains a nonpositive Jacobian")
        written_faces = prism_face_counts()
        written_faces == evidence.face_counts ||
            error("Serialized edge-cluster face topology changed after writing")
        written_attributes = sort!([
            Int(tag) for (_, tag) in gmsh.model.getPhysicalGroups(2)
        ])
        written_attributes == evidence.surface_attributes || error(
            "Serialized edge-cluster physical attributes changed after writing"
        )
        serialized_evidence = merge(
            evidence,
            (
                minimum_jacobian=written_jacobian,
                etch_full_gap=etch_full_gap,
                volume_element_types=sort!([
                    gmsh.model.mesh.getElementProperties(element_type)[1] for
                    element_type in written_types
                ]),
                face_counts=written_faces
            )
        )
        write_edge_cluster_metadata(
            filename * ".metadata.json",
            serialized_evidence,
            plan,
            z_coordinates,
            lc_normal,
            lc_tangent,
            lc_far,
            core_width,
            normal_growth_ratio,
            mesh_order,
            fabricated,
            length(edges)
        )
        println(
            "Edge-cluster transition coupon: edges=$(length(edges)), " *
            "fabricated=$fabricated, " *
            "nodes=$(serialized_evidence.node_count), " *
            "prisms=$(serialized_evidence.volume_count), " *
            "minSJ=$(serialized_evidence.minimum_jacobian), file=$filename"
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
