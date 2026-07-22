# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Local three-dimensional coupon for a convex 90 degree metal corner.
#
# Coordinates are in microns. The corner is at the origin, the two metal-edge
# arms point along +x and +y, and +z points from the substrate into vacuum.

import Gmsh: gmsh

function polygon_surface(occ, points)
    tags = [occ.addPoint(point...) for point in points]
    lines = [
        occ.addLine(tags[i], tags[mod1(i + 1, length(tags))])
        for i in eachindex(tags)
    ]
    return occ.addPlaneSurface([occ.addCurveLoop(lines)])
end

function extruded_polygon(occ, points, direction)
    surface = polygon_surface(occ, points)
    entities = occ.extrude([(2, surface)], direction...)
    return only(tag for (dim, tag) in entities if dim == 3)
end

function rounded_rectangle_wire(occ, half_x, half_y, radius, z)
    0.0 <= radius < min(half_x, half_y) ||
        error("Rounded-rectangle radius must lie in [0, min(half_x, half_y))")
    if radius == 0.0
        points = [
            occ.addPoint(-half_x, -half_y, z),
            occ.addPoint(half_x, -half_y, z),
            occ.addPoint(half_x, half_y, z),
            occ.addPoint(-half_x, half_y, z),
        ]
        curves = [
            occ.addLine(points[i], points[mod1(i + 1, length(points))])
            for i in eachindex(points)
        ]
        return occ.addWire(curves)
    end

    points = [
        occ.addPoint(half_x - radius, -half_y, z),
        occ.addPoint(half_x, -half_y + radius, z),
        occ.addPoint(half_x, half_y - radius, z),
        occ.addPoint(half_x - radius, half_y, z),
        occ.addPoint(-half_x + radius, half_y, z),
        occ.addPoint(-half_x, half_y - radius, z),
        occ.addPoint(-half_x, -half_y + radius, z),
        occ.addPoint(-half_x + radius, -half_y, z),
    ]
    centers = [
        occ.addPoint(half_x - radius, -half_y + radius, z),
        occ.addPoint(half_x - radius, half_y - radius, z),
        occ.addPoint(-half_x + radius, half_y - radius, z),
        occ.addPoint(-half_x + radius, -half_y + radius, z),
    ]
    curves = [
        occ.addCircleArc(points[1], centers[1], points[2]),
        occ.addLine(points[2], points[3]),
        occ.addCircleArc(points[3], centers[2], points[4]),
        occ.addLine(points[4], points[5]),
        occ.addCircleArc(points[5], centers[3], points[6]),
        occ.addLine(points[6], points[7]),
        occ.addCircleArc(points[7], centers[4], points[8]),
        occ.addLine(points[8], points[1]),
    ]
    return occ.addWire(curves)
end

function rounded_rectangle_surface(occ, half_x, half_y, radius, z)
    return occ.addPlaneSurface([
        rounded_rectangle_wire(occ, half_x, half_y, radius, z),
    ])
end

function tapered_rounded_rectangle(
    occ, half_x, half_y, radius, z0, height, pullback_bottom, pullback_top)
    bottom_half_x = half_x - pullback_bottom
    bottom_half_y = half_y - pullback_bottom
    top_half_x = half_x - pullback_top
    top_half_y = half_y - pullback_top
    bottom_radius = radius - pullback_bottom
    top_radius = radius - pullback_top
    minimum_radius = min(bottom_radius, top_radius)
    minimum_radius >= 0.0 ||
        error("Sidewall pullback exceeds the plan-view corner radius")

    bottom = rounded_rectangle_wire(
        occ, bottom_half_x, bottom_half_y, bottom_radius, z0)
    top = rounded_rectangle_wire(
        occ, top_half_x, top_half_y, top_radius, z0 + height)
    volumes = occ.addThruSections([bottom, top], -1, true, false, -1, "C0")
    result = [(dim, tag) for (dim, tag) in volumes if dim == 3]
    isempty(result) && error("Rounded-rectangle loft did not create a volume")
    return result
end

function rounded_convex_corner_wire(occ, extent, radius, offset, z)
    offset < extent ||
        error("Rounded-corner offset must be smaller than the extent")
    radius >= offset ||
        error("Sidewall pullback exceeds the plan-view corner radius")
    if radius == offset
        points = [
            occ.addPoint(offset, offset, z),
            occ.addPoint(extent, offset, z),
            occ.addPoint(extent, extent, z),
            occ.addPoint(offset, extent, z),
        ]
        curves = [
            occ.addLine(points[i], points[mod1(i + 1, length(points))])
            for i in eachindex(points)
        ]
        return occ.addWire(curves)
    end

    points = [
        occ.addPoint(radius, offset, z),
        occ.addPoint(extent, offset, z),
        occ.addPoint(extent, extent, z),
        occ.addPoint(offset, extent, z),
        occ.addPoint(offset, radius, z),
    ]
    center = occ.addPoint(radius, radius, z)
    curves = [
        occ.addLine(points[1], points[2]),
        occ.addLine(points[2], points[3]),
        occ.addLine(points[3], points[4]),
        occ.addLine(points[4], points[5]),
        occ.addCircleArc(points[5], center, points[1]),
    ]
    return occ.addWire(curves)
end

function tapered_rounded_convex_corner(
    occ, extent, radius, z0, height, offset_bottom, offset_top)
    bottom = rounded_convex_corner_wire(
        occ, extent, radius, offset_bottom, z0)
    top = rounded_convex_corner_wire(
        occ, extent, radius, offset_top, z0 + height)
    volumes = occ.addThruSections([bottom, top], -1, true, false, -1, "C0")
    result = [(dim, tag) for (dim, tag) in volumes if dim == 3]
    isempty(result) && error("Rounded-corner loft did not create a volume")
    return result
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

function fillet_volume_edges(occ, volumes, radius, selector, name)
    radius <= 0.0 && return volumes
    curves = [
        curve for curve in boundary_curves(volumes)
        if selector(gmsh.model.getBoundingBox(1, curve))
    ]
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

function add_tapered_corner_metal(
    occ, radius, thickness, sidewall_angle, convex)
    extension = 0.1 * radius
    pullback = thickness / tan(deg2rad(sidewall_angle))

    x_profile = if convex
        [
            (0.0, -extension, 0.0),
            (radius + extension, -extension, 0.0),
            (radius + extension, -extension, thickness),
            (pullback, -extension, thickness),
        ]
    else
        [
            (-radius - extension, -extension, 0.0),
            (0.0, -extension, 0.0),
            (-pullback, -extension, thickness),
            (-radius - extension, -extension, thickness),
        ]
    end
    x_prism =
        extruded_polygon(occ, x_profile, (0.0, radius + 2extension, 0.0))

    y_profile = if convex
        [
            (-extension, 0.0, 0.0),
            (-extension, radius + extension, 0.0),
            (-extension, radius + extension, thickness),
            (-extension, pullback, thickness),
        ]
    else
        [
            (-extension, -radius - extension, 0.0),
            (-extension, 0.0, 0.0),
            (-extension, -pullback, thickness),
            (-extension, -radius - extension, thickness),
        ]
    end
    y_prism =
        extruded_polygon(occ, y_profile, (radius + 2extension, 0.0, 0.0))

    metal, _ = convex ?
        occ.intersect([(3, x_prism)], [(3, y_prism)]) :
        occ.fuse([(3, x_prism)], [(3, y_prism)])
    return metal
end

function add_overetch_notch(occ, radius, depth, sidewall_angle, convex)
    extension = 0.1 * radius
    pullback = depth / tan(deg2rad(sidewall_angle))

    x_profile = if convex
        [
            (-radius - extension, -radius - extension, 0.0),
            (0.0, -radius - extension, 0.0),
            (-pullback, -radius - extension, -depth),
            (-radius - extension, -radius - extension, -depth),
        ]
    else
        [
            (0.0, -radius - extension, 0.0),
            (radius + extension, -radius - extension, 0.0),
            (radius + extension, -radius - extension, -depth),
            (pullback, -radius - extension, -depth),
        ]
    end
    x_notch =
        extruded_polygon(occ, x_profile, (0.0, 2radius + 2extension, 0.0))

    y_profile = if convex
        [
            (-radius - extension, -radius - extension, 0.0),
            (-radius - extension, 0.0, 0.0),
            (-radius - extension, -pullback, -depth),
            (-radius - extension, -radius - extension, -depth),
        ]
    else
        [
            (-radius - extension, 0.0, 0.0),
            (-radius - extension, radius + extension, 0.0),
            (-radius - extension, radius + extension, -depth),
            (-radius - extension, pullback, -depth),
        ]
    end
    y_notch =
        extruded_polygon(occ, y_profile, (2radius + 2extension, 0.0, 0.0))

    notch, _ = convex ?
        occ.fuse([(3, x_notch)], [(3, y_notch)]) :
        occ.intersect([(3, x_notch)], [(3, y_notch)])
    return notch
end

function on_outer_box(bounds, radius, tolerance)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    return (abs(xmin + radius) < tolerance &&
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
end

function generate_corner_coupon(;
    topology::Symbol = :convex,
    fabricated::Bool = false,
    radius::Float64 = 2.0,
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
    0.0 <= corner_radius < radius ||
        error("corner_radius must lie in [0, radius)")
    !convex && corner_radius > 0.0 &&
        error("Rounded concave coupons are not implemented yet")
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
            notch = if corner_radius > 0.0
                extension = 0.1 * radius
                trench_slab = [(
                    3,
                    occ.addBox(
                        -radius, -radius, -overetch_depth,
                        2radius, 2radius, overetch_depth),
                )]
                pedestal = tapered_rounded_convex_corner(
                    occ, radius + extension, corner_radius,
                    -overetch_depth, overetch_depth, -trench_pullback, 0.0)
                result, _ = occ.cut(trench_slab, pedestal)
                result
            else
                add_overetch_notch(
                    occ, radius, overetch_depth, sidewall_angle, convex)
            end
            substrate, _ = occ.cut(substrate, notch)
            trench_location =
                (convex ? -1.0 : 1.0) *
                trench_pullback
            tolerance = 1.0e-6 * radius
            substrate = fillet_volume_edges(
                occ, substrate, trench_rounding,
                bounds -> begin
                    xmin, ymin, zmin, xmax, ymax, zmax = bounds
                    on_floor = abs(zmin + overetch_depth) < tolerance &&
                               abs(zmax + overetch_depth) < tolerance
                    on_trench = if corner_radius > 0.0
                        xmin <= corner_radius + tolerance &&
                        ymin <= corner_radius + tolerance &&
                        xmax >= trench_location - tolerance &&
                        ymax >= trench_location - tolerance
                    else
                        (abs(xmin - trench_location) < tolerance &&
                         abs(xmax - trench_location) < tolerance) ||
                        (abs(ymin - trench_location) < tolerance &&
                         abs(ymax - trench_location) < tolerance)
                    end
                    on_floor && on_trench
                end,
                "trench-bottom",
            )
        end

        metal = if corner_radius > 0.0
            tapered_rounded_convex_corner(
                occ, 1.1radius, corner_radius,
                0.0, metal_thickness, 0.0,
                metal_thickness / tan(deg2rad(sidewall_angle)))
        else
            add_tapered_corner_metal(
                occ, radius, metal_thickness, sidewall_angle, convex)
        end
        top_location =
            (convex ? 1.0 : -1.0) *
            metal_thickness / tan(deg2rad(sidewall_angle))
        tolerance = 1.0e-6 * radius
        metal = fillet_volume_edges(
            occ, metal, top_rounding,
            bounds -> begin
                xmin, ymin, zmin, xmax, ymax, zmax = bounds
                on_top = abs(zmin - metal_thickness) < tolerance &&
                         abs(zmax - metal_thickness) < tolerance
                on_side = if corner_radius > 0.0
                    xmin <= corner_radius + tolerance &&
                    ymin <= corner_radius + tolerance
                else
                    (abs(xmin - top_location) < tolerance &&
                     abs(xmax - top_location) < tolerance) ||
                    (abs(ymin - top_location) < tolerance &&
                     abs(ymax - top_location) < tolerance)
                end
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
        footprint = if convex && corner_radius > 0.0
            wire = rounded_convex_corner_wire(
                occ, radius, corner_radius, 0.0, 0.0)
            occ.addPlaneSurface([wire])
        else
            polygon_surface(
                occ,
                convex ?
                    [(0.0, 0.0, 0.0), (radius, 0.0, 0.0),
                     (radius, radius, 0.0), (0.0, radius, 0.0)] :
                    [(-radius, -radius, 0.0), (radius, -radius, 0.0),
                     (radius, 0.0, 0.0), (0.0, 0.0, 0.0),
                     (0.0, radius, 0.0), (-radius, radius, 0.0)],
            )
        end
        lower = occ.extrude([(2, footprint)], 0.0, 0.0, -radius)
        upper = occ.extrude([(2, footprint)], 0.0, 0.0, radius)
        tools = [
            entity for entity in vcat(lower, upper)
            if entity[1] == 3
        ]
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
        if on_outer_box(bounds, radius, tolerance)
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
            center = gmsh.model.occ.getCenterOfMass(dim, tag)
            x, y = center[1], center[2]
            in_quadrant = x >= -tolerance && y >= -tolerance
            rounded_convex_metal =
                x >= corner_radius - tolerance ||
                y >= corner_radius - tolerance ||
                (x - corner_radius)^2 + (y - corner_radius)^2 <=
                    corner_radius^2 + tolerance^2
            on_metal = convex ?
                in_quadrant && rounded_convex_metal :
                !in_quadrant
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
        ("Mesh.MshFileVersion", 2.2),
        ("Mesh.Binary", 1),
    ]
        gmsh.option.setNumber(name, value)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(mesh_order)
    mesh_order > 1 && gmsh.model.mesh.optimize("HighOrderElastic")
    gmsh.write(filename)

    println("Corner coupon: topology=$(topology), fabricated=$(fabricated), " *
            "R=$(radius) um, corner radius=$(corner_radius) um")
    for (dim, tag) in gmsh.model.getPhysicalGroups()
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  dim=$dim tag=$tag name=$name ($(length(entities)) entities)")
    end
    println("  nodes=$(length(gmsh.model.mesh.getNodes()[1]))")
    println("  file=$(filename)")
    gmsh.finalize()
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) in (2, 3) ||
        error("Usage: mesh_corner_coupon.jl " *
              "[convex-|concave-]thin|fabricated OUTPUT.msh [CORNER_RADIUS_UM]")
    kind = ARGS[1]
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
    corner_radius = length(ARGS) == 3 ? parse(Float64, ARGS[3]) : 0.0
    generate_corner_coupon(
        topology = topology,
        fabricated = fabricated,
        corner_radius = corner_radius,
        filename = abspath(ARGS[2]),
    )
end
