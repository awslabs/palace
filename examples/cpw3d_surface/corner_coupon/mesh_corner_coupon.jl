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
            notch = add_overetch_notch(
                occ, radius, overetch_depth, sidewall_angle, convex)
            substrate, _ = occ.cut(substrate, notch)
            trench_location =
                (convex ? -1.0 : 1.0) *
                overetch_depth / tan(deg2rad(sidewall_angle))
            tolerance = 1.0e-6 * radius
            substrate = fillet_volume_edges(
                occ, substrate, trench_rounding,
                bounds -> begin
                    xmin, ymin, zmin, xmax, ymax, zmax = bounds
                    on_floor = abs(zmin + overetch_depth) < tolerance &&
                               abs(zmax + overetch_depth) < tolerance
                    on_trench = (abs(xmin - trench_location) < tolerance &&
                                 abs(xmax - trench_location) < tolerance) ||
                                (abs(ymin - trench_location) < tolerance &&
                                 abs(ymax - trench_location) < tolerance)
                    on_floor && on_trench
                end,
                "trench-bottom",
            )
        end

        metal = add_tapered_corner_metal(
            occ, radius, metal_thickness, sidewall_angle, convex)
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
                on_side = (abs(xmin - top_location) < tolerance &&
                           abs(xmax - top_location) < tolerance) ||
                          (abs(ymin - top_location) < tolerance &&
                           abs(ymax - top_location) < tolerance)
                on_top && on_side
            end,
            "top-metal",
        )
        outer = occ.addBox(
            -radius, -radius, -radius, 2radius, 2radius, 2radius)
        field, _ = occ.cut([(3, outer)], metal, -1, true, false)
        vacuum, _ = occ.cut(field, substrate, -1, true, false)
        domains, domain_map = occ.fragment(vcat(substrate, vacuum), [])
        substrate_seed = domain_map[1:length(substrate)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[length(substrate) + 1:end] |> Iterators.flatten |> collect
    else
        # Splitting at x=0 and y=0 makes the metal quadrant and exposed
        # substrate interface distinct geometric surfaces.
        volumes = Tuple{Int32,Int32}[]
        material = Symbol[]
        for (z0, depth, name) in ((-radius, radius, :substrate),
                                  (0.0, radius, :vacuum))
            for x0 in (-radius, 0.0), y0 in (-radius, 0.0)
                push!(volumes, (3, occ.addBox(x0, y0, z0, radius, radius, depth)))
                push!(material, name)
            end
        end
        domains, domain_map = occ.fragment(volumes, [])
        for (index, mapped) in enumerate(domain_map)
            if material[index] == :substrate
                append!(substrate_seed, mapped)
            else
                append!(vacuum_seed, mapped)
            end
        end
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
            in_quadrant = xmin >= -tolerance && ymin >= -tolerance
            on_metal = convex ? in_quadrant : !in_quadrant
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
            "R=$(radius) um")
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
    length(ARGS) == 2 ||
        error("Usage: mesh_corner_coupon.jl " *
              "[convex-|concave-]thin|fabricated OUTPUT.msh")
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
    generate_corner_coupon(
        topology = topology,
        fabricated = fabricated,
        filename = abspath(ARGS[2]),
    )
end
