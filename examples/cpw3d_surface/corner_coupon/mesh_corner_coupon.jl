# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Local three-dimensional coupon for a convex or concave metal corner.
#
# Coordinates are in microns. The corner is at the origin, the first metal-edge
# arm points along +x, the second is counterclockwise at the requested angle,
# and +z points from the substrate into vacuum.

import Gmsh: gmsh

function polygon_surface(occ, points)
    tags = [occ.addPoint(point...) for point in points]
    lines = [
        occ.addLine(tags[i], tags[mod1(i + 1, length(tags))])
        for i in eachindex(tags)
    ]
    return occ.addPlaneSurface([occ.addCurveLoop(lines)])
end

function corner_frame(angle_degrees)
    angle = deg2rad(angle_degrees)
    0.0 < angle < pi ||
        error("Corner angle must lie strictly between zero and 180 degrees")
    first = (1.0, 0.0)
    second = (cos(angle), sin(angle))
    first_normal = (0.0, 1.0)
    second_normal = (sin(angle), -cos(angle))
    bisector = (cos(0.5angle), sin(0.5angle))
    return (angle, first, second, first_normal, second_normal, bisector)
end

function rounded_wedge_wire(
    occ,
    extent,
    angle_degrees,
    corner_radius,
    offset,
    z,
)
    angle, first, second, first_normal, second_normal, bisector =
        corner_frame(angle_degrees)
    apex_distance = offset / sin(0.5angle)
    apex = (
        apex_distance * bisector[1],
        apex_distance * bisector[2],
    )
    outer_radius = 4.0extent + abs(apex_distance)
    outer_first = (
        apex[1] + outer_radius * first[1],
        apex[2] + outer_radius * first[2],
    )
    outer_second = (
        apex[1] + outer_radius * second[1],
        apex[2] + outer_radius * second[2],
    )
    apex_tag = occ.addPoint(apex[1], apex[2], z)
    outer_first_tag = occ.addPoint(outer_first[1], outer_first[2], z)
    outer_second_tag = occ.addPoint(outer_second[1], outer_second[2], z)
    outer_arc =
        occ.addCircleArc(outer_first_tag, apex_tag, outer_second_tag)

    if corner_radius == 0.0
        curves = [
            occ.addLine(apex_tag, outer_first_tag),
            outer_arc,
            occ.addLine(outer_second_tag, apex_tag),
        ]
        return occ.addWire(curves)
    end

    arc_radius = corner_radius - offset
    arc_radius > 0.0 ||
        error("Corner offset must be smaller than the plan-view corner radius")
    center_distance = corner_radius / sin(0.5angle)
    center = (
        center_distance * bisector[1],
        center_distance * bisector[2],
    )
    first_tangent = (
        center[1] - arc_radius * first_normal[1],
        center[2] - arc_radius * first_normal[2],
    )
    second_tangent = (
        center[1] - arc_radius * second_normal[1],
        center[2] - arc_radius * second_normal[2],
    )
    first_tangent_tag =
        occ.addPoint(first_tangent[1], first_tangent[2], z)
    second_tangent_tag =
        occ.addPoint(second_tangent[1], second_tangent[2], z)
    center_tag = occ.addPoint(center[1], center[2], z)
    curves = [
        occ.addLine(first_tangent_tag, outer_first_tag),
        outer_arc,
        occ.addLine(outer_second_tag, second_tangent_tag),
        occ.addCircleArc(second_tangent_tag, center_tag, first_tangent_tag),
    ]
    return occ.addWire(curves)
end

function tapered_rounded_wedge(
    occ,
    extent,
    angle_degrees,
    corner_radius,
    z0,
    height,
    offset_bottom,
    offset_top,
)
    bottom = rounded_wedge_wire(
        occ,
        extent,
        angle_degrees,
        corner_radius,
        offset_bottom,
        z0,
    )
    top = rounded_wedge_wire(
        occ,
        extent,
        angle_degrees,
        corner_radius,
        offset_top,
        z0 + height,
    )
    volumes = occ.addThruSections([bottom, top], -1, true, false, -1, "C0")
    result = [(dim, tag) for (dim, tag) in volumes if dim == 3]
    isempty(result) && error("Rounded-corner wedge loft did not create a volume")
    return result
end

function corner_coordinates(point, angle_degrees)
    _, _, _, first_normal, second_normal, _ =
        corner_frame(angle_degrees)
    return (
        first_normal[1] * point[1] + first_normal[2] * point[2],
        second_normal[1] * point[1] + second_normal[2] * point[2],
    )
end

function in_rounded_wedge(
    point,
    angle_degrees,
    corner_radius,
    offset,
    tolerance,
)
    first_distance, second_distance =
        corner_coordinates(point, angle_degrees)
    first_distance >= offset - tolerance &&
        second_distance >= offset - tolerance || return false
    corner_radius == 0.0 && return true

    angle, _, _, _, _, bisector = corner_frame(angle_degrees)
    center_distance = corner_radius / sin(0.5angle)
    center = (
        center_distance * bisector[1],
        center_distance * bisector[2],
    )
    arc_radius = corner_radius - offset
    return first_distance >= corner_radius - tolerance ||
           second_distance >= corner_radius - tolerance ||
           (point[1] - center[1])^2 + (point[2] - center[2])^2 <=
               arc_radius^2 + tolerance^2
end

function on_rounded_wedge_curve(
    point,
    angle_degrees,
    corner_radius,
    offset,
    tolerance,
)
    first_distance, second_distance =
        corner_coordinates(point, angle_degrees)
    if corner_radius == 0.0
        return (
            abs(first_distance - offset) <= tolerance ||
            abs(second_distance - offset) <= tolerance
        )
    end

    on_first =
        abs(first_distance - offset) <= tolerance &&
        second_distance >= corner_radius - tolerance
    on_second =
        abs(second_distance - offset) <= tolerance &&
        first_distance >= corner_radius - tolerance
    angle, _, _, _, _, bisector = corner_frame(angle_degrees)
    center_distance = corner_radius / sin(0.5angle)
    center = (
        center_distance * bisector[1],
        center_distance * bisector[2],
    )
    arc_radius = corner_radius - offset
    on_arc =
        abs(hypot(point[1] - center[1], point[2] - center[2]) - arc_radius) <=
            tolerance &&
        first_distance <= corner_radius + tolerance &&
        second_distance <= corner_radius + tolerance
    return on_first || on_second || on_arc
end

function boundary_curves(volumes)
    gmsh.model.occ.synchronize()
    surfaces = [
        entity for entity in gmsh.model.getBoundary(volumes, false, false, false)
        if entity[1] == 2
    ]
    return unique(tag for (dim, tag) in
                  gmsh.model.getBoundary(surfaces, false, false, false) if dim == 1)
end

function point_on_curve(curve)
    lower, upper = gmsh.model.getParametrizationBounds(1, curve)
    length(lower) == 1 && length(upper) == 1 ||
        error("Unexpected curve parametrization for curve $curve")
    value = gmsh.model.getValue(1, curve, [(lower[1] + upper[1]) / 2])
    return (value[1], value[2], value[3])
end

function fillet_volume_edges(occ, volumes, radius, selector, name)
    radius <= 0.0 && return volumes
    curves = Int32[]
    for curve in boundary_curves(volumes)
        bounds = gmsh.model.getBoundingBox(1, curve)
        selector(curve, bounds) && push!(curves, curve)
    end
    isempty(curves) && error("No $name curves were found for radius $radius")
    rounded = occ.fillet(
        Int32[tag for (dim, tag) in volumes if dim == 3],
        Int32.(curves),
        [radius],
    )
    result = [(dim, tag) for (dim, tag) in rounded if dim == 3]
    isempty(result) && error("Filleting $name removed every volume")
    return result
end

function on_outer_box(bounds, center, radius, tolerance)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    bounds_on_box =
        (abs(xmin + radius) < tolerance &&
         abs(xmax + radius) < tolerance) ||
        (abs(xmin - radius) < tolerance &&
         abs(xmax - radius) < tolerance) ||
        (abs(ymin + radius) < tolerance &&
         abs(ymax + radius) < tolerance) ||
        (abs(ymin - radius) < tolerance &&
         abs(ymax - radius) < tolerance) ||
        (abs(zmin + radius) < tolerance &&
         abs(zmax + radius) < tolerance) ||
        (abs(zmin - radius) < tolerance &&
         abs(zmax - radius) < tolerance)
    bounds_on_box && return true

    # OCC Boolean operations involving circular wedge closures can enlarge the
    # pre-mesh bounding box of an otherwise planar box face. Its centroid still
    # lies exactly on the matching box and is a more stable secondary test.
    spans = (xmax - xmin, ymax - ymin, zmax - zmin)
    return any(
        abs(abs(center[axis]) - radius) < tolerance &&
        spans[axis] < 1.0e-3radius for axis in 1:3
    )
end

function generate_corner_coupon(;
    topology::Symbol = :convex,
    fabricated::Bool = false,
    radius::Float64 = 2.0,
    angle_degrees::Float64 = 90.0,
    corner_radius::Float64 = 0.0,
    metal_thickness::Float64 = 0.1,
    overetch_depth::Float64 = 0.05,
    sidewall_angle::Float64 = 80.0,
    top_rounding::Float64 = 0.01,
    trench_rounding::Float64 = 0.01,
    lc_fine::Float64 = 0.02,
    lc_far::Float64 = 0.3,
    mesh_order::Int = 1,
    filename::String,
)
    radius > 0.0 || error("radius must be positive")
    topology in (:convex, :concave) ||
        error("topology must be :convex or :concave")
    convex = topology == :convex
    0.0 < angle_degrees < 180.0 ||
        error("angle_degrees must lie strictly between zero and 180")
    0.0 <= corner_radius < radius ||
        error("corner_radius must lie in [0, radius)")
    if corner_radius > 0.0
        tangent_distance =
            corner_radius / tan(0.5deg2rad(angle_degrees))
        tangent_distance < radius ||
            error("Rounded-corner tangency points must lie inside the matching box")
    end
    0.0 < sidewall_angle <= 90.0 ||
        error("sidewall_angle must be in (0, 90] degrees")
    if fabricated
        metal_thickness > 0.0 || error("metal_thickness must be positive")
        overetch_depth >= 0.0 || error("overetch_depth cannot be negative")
        0.0 <= top_rounding < metal_thickness ||
            error("top_rounding must be nonnegative and smaller than metal thickness")
        0.0 <= trench_rounding <= overetch_depth ||
            error("trench_rounding must lie between zero and overetch depth")
    end

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("$(topology)_corner_$(fabricated ? "fabricated" : "thin")")
    occ = gmsh.model.occ

    substrate_seed = Tuple{Int32,Int32}[]
    vacuum_seed = Tuple{Int32,Int32}[]
    if fabricated
        substrate_base =
            occ.addBox(-radius, -radius, -radius, 2radius, 2radius, radius)
        substrate = [(3, substrate_base)]
        if overetch_depth > 0.0
            trench_pullback =
                overetch_depth / tan(deg2rad(sidewall_angle))
            corner_radius > 0.0 && !convex &&
                corner_radius <= trench_pullback &&
                error("Trench pullback exceeds the concave corner radius")
            trench_slab = [(
                3,
                occ.addBox(
                    -radius,
                    -radius,
                    -overetch_depth,
                    2radius,
                    2radius,
                    overetch_depth,
                ),
            )]
            trench_wedge = tapered_rounded_wedge(
                occ,
                radius,
                angle_degrees,
                corner_radius,
                -overetch_depth,
                overetch_depth,
                convex ? -trench_pullback : trench_pullback,
                0.0,
            )
            notch = if convex
                result, _ = occ.cut(trench_slab, trench_wedge)
                result
            else
                result, _ = occ.intersect(trench_slab, trench_wedge)
                result
            end
            substrate, _ = occ.cut(substrate, notch)
            trench_offset = convex ? -trench_pullback : trench_pullback
            tolerance = 1.0e-6 * radius
            substrate = fillet_volume_edges(
                occ, substrate, trench_rounding,
                (curve, bounds) -> begin
                    xmin, ymin, zmin, xmax, ymax, zmax = bounds
                    on_floor = abs(zmin + overetch_depth) < tolerance &&
                               abs(zmax + overetch_depth) < tolerance
                    point = point_on_curve(curve)
                    on_trench = on_rounded_wedge_curve(
                        point,
                        angle_degrees,
                        corner_radius,
                        trench_offset,
                        tolerance,
                    )
                    on_floor && on_trench
                end,
                "trench-bottom",
            )
        end

        metal_pullback =
            metal_thickness / tan(deg2rad(sidewall_angle))
        corner_radius > 0.0 && convex &&
            corner_radius <= metal_pullback &&
            error("Sidewall pullback exceeds the convex corner radius")
        metal_slab = [(
            3,
            occ.addBox(
                -1.1radius,
                -1.1radius,
                0.0,
                2.2radius,
                2.2radius,
                metal_thickness,
            ),
        )]
        metal_offset = convex ? metal_pullback : -metal_pullback
        metal_wedge = tapered_rounded_wedge(
            occ,
            1.1radius,
            angle_degrees,
            corner_radius,
            0.0,
            metal_thickness,
            0.0,
            metal_offset,
        )
        metal = if convex
            result, _ = occ.intersect(metal_slab, metal_wedge)
            result
        else
            result, _ = occ.cut(metal_slab, metal_wedge)
            result
        end
        tolerance = 1.0e-6 * radius
        metal = fillet_volume_edges(
            occ, metal, top_rounding,
            (curve, bounds) -> begin
                xmin, ymin, zmin, xmax, ymax, zmax = bounds
                on_top = abs(zmin - metal_thickness) < tolerance &&
                         abs(zmax - metal_thickness) < tolerance
                point = point_on_curve(curve)
                on_side = on_rounded_wedge_curve(
                    point,
                    angle_degrees,
                    corner_radius,
                    metal_offset,
                    tolerance,
                )
                on_top && on_side
            end,
            "top-metal",
        )
        outer = occ.addBox(
            -radius, -radius, -radius, 2radius, 2radius, 2radius)
        field, _ = occ.cut([(3, outer)], metal, -1, true, true)
        vacuum, _ = occ.cut(field, substrate, -1, true, false)
        domains, domain_map = occ.fragment(vcat(substrate, vacuum), [])
        substrate_seed = domain_map[1:length(substrate)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[length(substrate) + 1:end] |> Iterators.flatten |> collect
    else
        substrate_box = (3, occ.addBox(
            -radius, -radius, -radius, 2radius, 2radius, radius))
        vacuum_box = (3, occ.addBox(
            -radius, -radius, 0.0, 2radius, 2radius, radius))
        square = polygon_surface(
            occ,
            [
                (-radius, -radius, 0.0),
                (radius, -radius, 0.0),
                (radius, radius, 0.0),
                (-radius, radius, 0.0),
            ],
        )
        wedge = occ.addPlaneSurface([
            rounded_wedge_wire(
                occ,
                radius,
                angle_degrees,
                corner_radius,
                0.0,
                0.0,
            ),
        ])
        footprints = if convex
            result, _ = occ.intersect([(2, square)], [(2, wedge)])
            result
        else
            result, _ = occ.cut([(2, square)], [(2, wedge)])
            result
        end
        footprint_surfaces =
            [(dim, tag) for (dim, tag) in footprints if dim == 2]
        isempty(footprint_surfaces) &&
            error("Corner footprint does not intersect the matching box")
        tools = Tuple{Int32,Int32}[]
        for footprint in footprint_surfaces
            lower = occ.extrude([footprint], 0.0, 0.0, -radius)
            upper = occ.extrude([footprint], 0.0, 0.0, radius)
            append!(
                tools,
                [
                    entity for entity in vcat(lower, upper)
                    if entity[1] == 3
                ],
            )
        end
        domains, domain_map =
            occ.fragment([substrate_box, vacuum_box], tools)
        append!(substrate_seed, domain_map[1])
        append!(vacuum_seed, domain_map[2])
    end
    occ.synchronize()

    domain_tags = Set(tag for (dim, tag) in domains if dim == 3)
    substrate_tags =
        sort!(unique(tag for (dim, tag) in substrate_seed
                     if dim == 3 && tag in domain_tags))
    vacuum_tags =
        sort!(unique(tag for (dim, tag) in vacuum_seed
                     if dim == 3 && tag in domain_tags))
    isempty(substrate_tags) && error("No substrate volumes were generated")
    isempty(vacuum_tags) && error("No vacuum volumes were generated")

    substrate_set = Set(substrate_tags)
    vacuum_set = Set(vacuum_tags)
    tolerance = 1.0e-6 * radius
    matching = Int32[]
    thin_metal = Int32[]
    ms = Int32[]
    ma = Int32[]
    sa = Int32[]

    for (dim, tag) in gmsh.model.getEntities(2)
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        adjacent_substrate = [volume for volume in up if volume in substrate_set]
        adjacent_vacuum = [volume for volume in up if volume in vacuum_set]
        isempty(adjacent_substrate) && isempty(adjacent_vacuum) && continue

        bounds = gmsh.model.getBoundingBox(dim, tag)
        center = gmsh.model.occ.getCenterOfMass(dim, tag)
        if on_outer_box(bounds, center, radius, tolerance)
            push!(matching, tag)
            continue
        end

        if fabricated
            if !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                push!(sa, tag)
            elseif !isempty(adjacent_substrate)
                push!(ms, tag)
            elseif !isempty(adjacent_vacuum)
                push!(ma, tag)
            end
        else
            xmin, ymin, zmin, xmax, ymax, zmax = bounds
            on_interface = abs(zmin) < tolerance && abs(zmax) < tolerance
            x, y = center[1], center[2]
            in_wedge = in_rounded_wedge(
                (x, y),
                angle_degrees,
                corner_radius,
                0.0,
                tolerance,
            )
            on_metal = convex ?
                in_wedge :
                !in_wedge
            if on_interface && on_metal
                push!(thin_metal, tag)
            elseif on_interface &&
                   !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                push!(sa, tag)
            end
        end
    end

    groups = [
        (3, substrate_tags, 1, "substrate"),
        (3, vacuum_tags, 2, "vacuum"),
        (2, matching, 1, "matching_surface"),
    ]
    if fabricated
        append!(groups, [
            (2, ms, 2, "MS"),
            (2, sa, 3, "SA"),
            (2, ma, 4, "MA"),
        ])
    else
        append!(groups, [
            (2, thin_metal, 2, "thin_metal"),
            (2, sa, 3, "SA"),
        ])
    end
    for (dim, entities, tag, name) in groups
        isempty(entities) && error("Empty physical group: $name")
        gmsh.model.addPhysicalGroup(dim, entities, tag, name)
    end

    feature_surfaces = fabricated ? vcat(ms, sa, ma) : vcat(thin_metal, sa)
    feature_curves = Int32[]
    for surface in feature_surfaces
        for (dim, curve) in gmsh.model.getBoundary(
            [(2, surface)], false, false, false)
            dim == 1 || continue
            # Curves where a fabrication plane intersects the matching surface
            # must be resolved as carefully as the interior process geometry.
            # The spatial response basis has knots on these planes.
            push!(feature_curves, curve)
        end
    end
    unique!(feature_curves)
    isempty(feature_curves) && error("No process-feature curves were generated")
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(feature_curves))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2lc_fine)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5radius)
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
        ("Mesh.Binary", 1),
    ]
        gmsh.option.setNumber(name, value)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.setOrder(mesh_order)
    if mesh_order > 1
        gmsh.model.mesh.optimize("HighOrderElastic", true, 20)
        gmsh.model.mesh.optimize("HighOrder", true, 20)
        _, element_tags, _ = gmsh.model.mesh.getElements(3)
        tags = reduce(vcat, element_tags; init = UInt64[])
        scaled_jacobians =
            gmsh.model.mesh.getElementQualities(tags, "minSJ")
        minimum_jacobian = minimum(scaled_jacobians)
        minimum_jacobian > 0.0 ||
            error("High-order corner coupon contains a nonpositive Jacobian")
        println("High-order mesh minimum scaled Jacobian: $minimum_jacobian")
    end
    gmsh.write(filename)

    println("Corner coupon: topology=$(topology), fabricated=$(fabricated), " *
            "angle=$(angle_degrees) deg, R=$(radius) um, " *
            "corner radius=$(corner_radius) um")
    for (dim, tag) in gmsh.model.getPhysicalGroups()
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  dim=$dim tag=$tag name=$name ($(length(entities)) entities)")
    end
    println("  nodes=$(length(gmsh.model.mesh.getNodes()[1]))")
    println("  file=$(filename)")
    gmsh.finalize()
end

function parse_command_line(args)
    length(args) >= 2 ||
        error("Usage: mesh_corner_coupon.jl " *
              "[convex-|concave-]thin|fabricated OUTPUT.msh [CORNER_RADIUS_UM] " *
              "[--radius VALUE ...]")
    kind = args[1]
    aliases = Dict(
        "thin" => (:convex, false),
        "fabricated" => (:convex, true),
        "convex-thin" => (:convex, false),
        "convex-fabricated" => (:convex, true),
        "concave-thin" => (:concave, false),
        "concave-fabricated" => (:concave, true),
    )
    haskey(aliases, kind) || error("Unknown coupon kind: $kind")
    topology, fabricated = aliases[kind]
    options = Dict{String,Any}(
        "topology" => topology,
        "fabricated" => fabricated,
        "filename" => abspath(args[2]),
        "angle_degrees" => 90.0,
        "corner_radius" => 0.0,
    )
    names = Dict(
        "--radius" => ("radius", Float64),
        "--angle" => ("angle_degrees", Float64),
        "--corner-radius" => ("corner_radius", Float64),
        "--metal-thickness" => ("metal_thickness", Float64),
        "--overetch" => ("overetch_depth", Float64),
        "--sidewall-angle" => ("sidewall_angle", Float64),
        "--top-radius" => ("top_rounding", Float64),
        "--bottom-radius" => ("trench_rounding", Float64),
        "--lc-fine" => ("lc_fine", Float64),
        "--lc-far" => ("lc_far", Float64),
        "--mesh-order" => ("mesh_order", Int),
    )
    index = 3
    if index <= length(args) && !startswith(args[index], "--")
        options["corner_radius"] = parse(Float64, args[index])
        index += 1
    end
    while index <= length(args)
        flag = args[index]
        haskey(names, flag) || error("Unknown option: $flag")
        index < length(args) || error("Missing value for option: $flag")
        name, type = names[flag]
        options[name] = parse(type, args[index + 1])
        index += 2
    end
    return options
end

if abspath(PROGRAM_FILE) == @__FILE__
    options = parse_command_line(ARGS)
    generate_corner_coupon(;
        topology = options["topology"],
        fabricated = options["fabricated"],
        filename = options["filename"],
        radius = get(options, "radius", 2.0),
        angle_degrees = options["angle_degrees"],
        corner_radius = options["corner_radius"],
        metal_thickness = get(options, "metal_thickness", 0.1),
        overetch_depth = get(options, "overetch_depth", 0.05),
        sidewall_angle = get(options, "sidewall_angle", 80.0),
        top_rounding = get(options, "top_rounding", 0.01),
        trench_rounding = get(options, "trench_rounding", 0.01),
        lc_fine = get(options, "lc_fine", 0.02),
        lc_far = get(options, "lc_far", 0.3),
        mesh_order = get(options, "mesh_order", 1),
    )
end
