# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabrication-resolved local coupon for endpoint, junction, and exact spatial clusters.

import Gmsh: gmsh
using DelimitedFiles

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
                slot=Int(round(data[row, columns["Slot"]])),
                conductor=Int(round(data[row, columns["Conductor"]])),
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
                normal_sign=sign(Float64(data[row, columns["Nz"]])),
                interval=(
                    Float64(data[row, columns["S0"]]),
                    Float64(data[row, columns["S1"]])
                ),
                vertex_arm=Bool(round(Int, data[row, columns["VertexArm"]]))
            )
        )
    end
    isempty(edges) && error("Spatial signature contains no edges")
    all(edge.slot >= 0 && edge.conductor > 0 for edge in edges) ||
        error("Spatial signature has invalid slot or conductor labels")
    all(
        abs(edge.gap[3]) < 1.0e-8 &&
            abs(edge.tangent[3]) < 1.0e-8 &&
            abs(edge.normal_sign) == 1 for edge in edges
    ) || error("Spatial coupon requires parallel or antiparallel process planes")
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
            (
                conductor=conductor,
                plane=plane,
                hole=hole,
                points=points,
                classes=classes
            )
        )
    end
    isempty(loops) && error("Plan-view boundary contains no loops")
    return loops
end

add(a, b) = ntuple(i -> a[i] + b[i], 3)
scale(value, vector) = ntuple(i -> value * vector[i], 3)

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
        (
            points[index][1] - circle.center[1],
            points[index][2] - circle.center[2]
        ) for index in point_indices
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
        abs(hypot(points[index][1] - circle.center[1],
                  points[index][2] - circle.center[2]) - circle.radius) for
        index in point_indices
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
        hypot(
            first.center[1] - second.center[1],
            first.center[2] - second.center[2]
        ) <= max(32tolerance, 1.0e-7 * max(first.radius, second.radius)) &&
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
        point_indices = [
            mod1(seed + step, length(points)) for step = 0:(triple_count + 1)
        ]
        edge_indices = [
            mod1(seed + step, length(points)) for step = 0:triple_count
        ]
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
    tolerance = 1.0e-9 * max(
        maximum(point[1] for point in points) - minimum(point[1] for point in points),
        maximum(point[2] for point in points) - minimum(point[2] for point in points),
        1.0
    )
    runs = circular_arc_runs(points, tolerance)
    projected = collect(points)
    for run in runs, index in run.point_indices[2:(end - 1)]
        radial = (
            points[index][1] - run.center[1],
            points[index][2] - run.center[2]
        )
        scale = run.radius / hypot(radial...)
        projected[index] = (
            run.center[1] + scale * radial[1],
            run.center[2] + scale * radial[2]
        )
    end
    tags = [occ.addPoint(point[1], point[2], z) for point in projected]
    curve_for_edge = Dict{Int, Int32}()
    covered = falses(length(points))
    for run in runs
        any(covered[run.edge_indices]) && continue
        parts = max(1, ceil(Int, run.angle / (0.5 * pi)))
        split = unique(
            round.(
                Int,
                range(1, length(run.point_indices), length=parts + 1)
            )
        )
        center = occ.addPoint(run.center[1], run.center[2], z)
        for (first, second) in zip(split, split[2:end])
            curve_for_edge[run.edge_indices[first]] =
                occ.addCircleArc(
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

function offset_loop_points(loop, distance, tolerance)
    distance == 0.0 && return loop.points
    shifted = Tuple{NTuple{2, Float64}, NTuple{2, Float64}}[]
    for index in eachindex(loop.points)
        first = loop.points[index]
        second = loop.points[mod1(index + 1, length(loop.points))]
        direction = (second[1] - first[1], second[2] - first[2])
        segment_length = hypot(direction...)
        segment_length > tolerance ||
            error("Plan-view boundary contains a zero-length segment")
        shift = loop.classes[index] == "Physical" ? distance : 0.0
        normal = (-direction[2] / segment_length, direction[1] / segment_length)
        push!(
            shifted,
            (
                (first[1] + shift * normal[1], first[2] + shift * normal[2]),
                direction
            )
        )
    end
    points = NTuple{2, Float64}[]
    for index in eachindex(shifted)
        previous = shifted[mod1(index - 1, length(shifted))]
        current = shifted[index]
        denominator = cross2d(previous[2], current[2])
        abs(denominator) > tolerance ||
            error("Plan-view taper has a singular boundary vertex")
        offset = (
            current[1][1] - previous[1][1],
            current[1][2] - previous[1][2]
        )
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

function loft_mask_offsets(
    occ,
    loops,
    z0,
    z1,
    bottom_offset,
    top_offset,
    tolerance
)
    outers = [loop for loop in loops if !loop.hole]
    holes = [loop for loop in loops if loop.hole]
    isempty(outers) && error("Plan-view mask has no exterior loop")
    result = Tuple{Int32, Int32}[]
    for outer in outers
        volume = loft_polygon(
            occ,
            offset_loop_points(outer, bottom_offset, tolerance),
            offset_loop_points(outer, top_offset, tolerance),
            z0,
            z1
        )
        cutters = Tuple{Int32, Int32}[]
        for hole in holes
            point_in_polygon(hole.points[1], outer.points, tolerance) || continue
            append!(
                cutters,
                loft_polygon(
                    occ,
                    offset_loop_points(hole, bottom_offset, tolerance),
                    offset_loop_points(hole, top_offset, tolerance),
                    z0,
                    z1
                )
            )
        end
        if !isempty(cutters)
            volume, _ = occ.cut(volume, cutters)
            volume = [(dim, tag) for (dim, tag) in volume if dim == 3]
        end
        append!(result, volume)
    end
    return fuse_all(occ, result)
end

function loft_mask(occ, loops, z0, z1, pullback, tolerance)
    return loft_mask_offsets(occ, loops, z0, z1, 0.0, pullback, tolerance)
end

function boundary_strips(occ, loops, radius, z0, z1, pullback, tolerance)
    volumes = Tuple{Int32, Int32}[]
    width = 3radius
    for conductor in sort!(unique(loop.conductor for loop in loops))
        conductor_loops = [loop for loop in loops if loop.conductor == conductor]
        expanded = loft_mask_offsets(
            occ,
            conductor_loops,
            z0,
            z1,
            -width,
            -width,
            tolerance
        )
        retained = loft_mask_offsets(
            occ,
            conductor_loops,
            z0,
            z1,
            0.0,
            -pullback,
            tolerance
        )
        strip, _ = occ.cut(expanded, retained)
        strip = [(dim, tag) for (dim, tag) in strip if dim == 3]
        isempty(strip) && error("Classified boundary strip produced no volume")
        append!(volumes, strip)
    end
    return fuse_all(occ, volumes)
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

function fillet_plane_edges(occ, volumes, radius, z, tolerance)
    radius <= 0.0 && return volumes
    curves = Int32[]
    for curve in boundary_curves(volumes)
        _, _, zmin, _, _, zmax = gmsh.model.getBoundingBox(1, curve)
        if abs(zmin - z) < tolerance && abs(zmax - z) < tolerance
            push!(curves, curve)
        end
    end
    isempty(curves) && return volumes
    rounded = occ.fillet(Int32[tag for (dim, tag) in volumes if dim == 3], curves, [radius])
    result = [(dim, tag) for (dim, tag) in rounded if dim == 3]
    return isempty(result) ? volumes : result
end

function point_segment_distance(point, first, second)
    direction = (second[1] - first[1], second[2] - first[2])
    length_squared = direction[1]^2 + direction[2]^2
    length_squared > 0.0 || return hypot(point[1] - first[1], point[2] - first[2])
    coordinate = clamp(
        ((point[1] - first[1]) * direction[1] +
         (point[2] - first[2]) * direction[2]) / length_squared,
        0.0,
        1.0
    )
    closest = (
        first[1] + coordinate * direction[1],
        first[2] + coordinate * direction[2]
    )
    return hypot(point[1] - closest[1], point[2] - closest[2])
end

function physical_segments(loops, offset, tolerance)
    primitives = NamedTuple[]
    for loop in loops
        points = offset_loop_points(loop, offset, tolerance)
        covered = falses(length(points))
        for run in circular_arc_runs(points, tolerance)
            all(loop.classes[index] == "Physical" for index in run.edge_indices) ||
                continue
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
    angle = orientation * atan(cross2d(first, second), first[1] * second[1] +
                                                       first[2] * second[2])
    return mod(angle, 2pi)
end

function point_primitive_distance(point, primitive, tolerance)
    primitive.kind == :line &&
        return point_segment_distance(point, primitive.first, primitive.last)
    radial = (
        point[1] - primitive.center[1],
        point[2] - primitive.center[2]
    )
    start = (
        primitive.first[1] - primitive.center[1],
        primitive.first[2] - primitive.center[2]
    )
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
        _, _, zmin, _, _, zmax = gmsh.model.getBoundingBox(1, curve)
        abs(zmin - z) < tolerance && abs(zmax - z) < tolerance || continue
        point = point_on_curve(curve)
        any(
            point_primitive_distance(point, primitive, tolerance) <= 10tolerance for
            primitive in primitives
        ) && push!(curves, curve)
    end
    isempty(curves) && return volumes
    rounded = occ.fillet(Int32[tag for (dim, tag) in volumes if dim == 3], curves, [radius])
    result = [(dim, tag) for (dim, tag) in rounded if dim == 3]
    return isempty(result) ? volumes : result
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
    coordinate =
        clamp(delta[1] * edge.tangent[1] + delta[2] * edge.tangent[2], first, second)
    closest = add(edge.point, scale(coordinate, edge.tangent))
    return sqrt(sum((point[index] - closest[index])^2 for index = 1:3))
end

function nearest_edge(edges, point, radius)
    return edges[argmin(segment_distance(edge, point, radius) for edge in edges)]
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
    boxes = [gmsh.model.getBoundingBox(2, surface) for surface in surfaces]
    for axis in 1:3
        if all(abs(box[axis + 3] - box[axis]) <= tolerance for box in boxes)
            coordinate = boxes[1][axis]
            if all(abs(box[axis] - coordinate) <= tolerance for box in boxes)
                return true
            end
        end
    end
    return false
end

function generate_spatial_coupon(;
    signature::String,
    mask::Union{Nothing, String}=nothing,
    boundary::Union{Nothing, String}=nothing,
    fabricated::Bool,
    radius::Float64          = 2.0,
    metal_thickness::Float64 = 0.1,
    overetch::Float64        = 0.05,
    sidewall_angle::Float64  = 80.0,
    top_rounding::Float64    = 0.01,
    trench_rounding::Float64 = 0.01,
    lc_fine::Float64         = 0.02,
    lc_far::Float64          = 0.3,
    process_core_width::Float64 = 0.0,
    max_nodes::Int           = 500_000,
    max_elements::Int        = 2_000_000,
    mesh_order::Int          = 1,
    filename::String
)
    radius > 0.0 || error("radius must be positive")
    metal_thickness > 0.0 || error("metal thickness must be positive")
    0.0 <= overetch < radius || error("overetch must lie in [0, radius)")
    0.0 < sidewall_angle <= 90.0 || error("sidewall angle must lie in (0, 90]")
    0.0 <= top_rounding < metal_thickness ||
        error("top rounding must be smaller than metal thickness")
    0.0 <= trench_rounding <= overetch || error("trench rounding must not exceed overetch")
    lc_fine > 0.0 || error("fine mesh size must be positive")
    lc_far >= lc_fine || error("far mesh size must not be smaller than fine mesh size")
    process_core_width = process_core_width > 0.0 ?
                         process_core_width :
                         max(2metal_thickness, 4overetch, 8lc_fine)
    process_core_width > 2lc_fine ||
        error("process core width must exceed twice the fine mesh size")
    max_nodes > 0 || error("maximum node budget must be positive")
    max_elements > 0 || error("maximum element budget must be positive")

    edges = read_edges(signature)
    facets = read_mask(mask)
    boundary_loops = read_boundary(boundary)
    isempty(boundary_loops) || !isempty(facets) ||
        error("A classified plan-view boundary requires the corresponding mask facets")
    lower, upper = coupon_bounds(edges, radius, metal_thickness, overetch)
    tolerance = 1.0e-7 * radius
    outer_tolerance = 1.0e-4 * radius
    validate_plan_view_geometry(edges, radius, tolerance, facets)
    layers = layer_groups(edges, tolerance)
    pullback_metal = metal_thickness / tan(deg2rad(sidewall_angle))
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
                loop for loop in boundary_loops if
                abs(loop.plane - layer.plane) <= tolerance
            ]
            trenches = if isempty(boundary_loops)
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
                isempty(layer_loops) &&
                    error("Plan-view boundary is missing fabrication layer $(layer.plane)")
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
                    physical_segments(
                        layer_loops,
                        -pullback_trench,
                        tolerance
                    ),
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
                    loop for loop in boundary_loops if
                    loop.conductor == conductor &&
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
                    apply_plan_view_mask(
                        occ,
                        result,
                        conductor_facets,
                        lower,
                        upper
                    )
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
            conductor_tools = Tuple{Int32, Int32}[]
            for edge in conductor_edges
                append!(
                    conductor_tools,
                    extruded_strip(occ, edge, radius, -1.0, lower[3], depth)
                )
            end
            conductor_tools = fuse_all(occ, conductor_tools)
            conductor_facets = [facet for facet in facets if facet.conductor == conductor]
            !isempty(facets) &&
                isempty(conductor_facets) &&
                error("Plan-view mask is missing conductor $conductor")
            conductor_tools =
                apply_plan_view_mask(occ, conductor_tools, conductor_facets, lower, upper)
            append!(tools, conductor_tools)
        end
        domains, domain_map = occ.fragment(vcat(substrates, vacuum), tools)
        substrate_seed = domain_map[1:length(substrates)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[(length(substrates) + 1):(length(substrates) + length(vacuum))] |>
            Iterators.flatten |>
            collect
    end
    occ.synchronize()

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
    println(
        "Spatial mesh features: candidates=$(length(candidate_curves)), " *
        "physical=$(length(feature_curves)), " *
        "discarded_coplanar_seams=$(length(discarded_seams)), " *
        "core_width=$process_core_width"
    )
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(feature_curves))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2lc_fine)
    gmsh.model.mesh.field.setNumber(2, "DistMax", process_core_width)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    for (name, value) in [
        ("Mesh.MeshSizeMin", lc_fine),
        ("Mesh.MeshSizeMax", lc_far),
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
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
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
    gmsh.write(filename)
    metadata_path = filename * ".metadata.json"
    open(metadata_path, "w") do stream
        println(stream, "{")
        println(stream, "  \"Version\": 1,")
        println(stream, "  \"MetalSurfacePartition\": \"InterfaceSlotAndConductor\",")
        println(stream, "  \"NodeCount\": $node_count,")
        println(stream, "  \"VolumeElementCount\": $element_count,")
        println(stream, "  \"FineSize\": $lc_fine,")
        println(stream, "  \"FarSize\": $lc_far,")
        println(stream, "  \"ProcessCoreWidth\": $process_core_width,")
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
        "--lc-far" => ("lc_far", Float64),
        "--process-core-width" => ("process_core_width", Float64),
        "--max-nodes" => ("max_nodes", Int),
        "--max-elements" => ("max_elements", Int),
        "--mesh-order" => ("mesh_order", Int)
    )
    index = 4
    while index <= length(args)
        flag = args[index]
        haskey(names, flag) || error("Unknown option: $flag")
        index < length(args) || error("Missing value for option: $flag")
        name, type = names[flag]
        options[name] =
            type === String ? abspath(args[index + 1]) : parse(type, args[index + 1])
        index += 2
    end
    return options
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_options(ARGS)
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
        lc_far          = get(options, "lc_far", 0.3),
        process_core_width = get(options, "process_core_width", 0.0),
        max_nodes       = get(options, "max_nodes", 500_000),
        max_elements    = get(options, "max_elements", 2_000_000),
        mesh_order      = get(options, "mesh_order", 1)
    )
end
