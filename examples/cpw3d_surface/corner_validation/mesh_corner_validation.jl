# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabrication-resolved and thin-metal validation meshes for either a finite
# rectangular island or a rectangular aperture in a metal sheet.

import Gmsh: gmsh

include(joinpath(@__DIR__, "..", "corner_coupon", "mesh_corner_coupon.jl"))

function tapered_rectangle(
    occ, half_x, half_y, z0, height, pullback_bottom, pullback_top)
    extension = 0.1 * max(half_x, half_y)
    x_profile = [
        (-half_x + pullback_bottom, -half_y - extension, z0),
        (half_x - pullback_bottom, -half_y - extension, z0),
        (half_x - pullback_top, -half_y - extension, z0 + height),
        (-half_x + pullback_top, -half_y - extension, z0 + height),
    ]
    x_prism = extruded_polygon(
        occ, x_profile, (0.0, 2half_y + 2extension, 0.0))
    y_profile = [
        (-half_x - extension, -half_y + pullback_bottom, z0),
        (-half_x - extension, half_y - pullback_bottom, z0),
        (-half_x - extension, half_y - pullback_top, z0 + height),
        (-half_x - extension, -half_y + pullback_top, z0 + height),
    ]
    y_prism = extruded_polygon(
        occ, y_profile, (2half_x + 2extension, 0.0, 0.0))
    volume, _ = occ.intersect([(3, x_prism)], [(3, y_prism)])
    return volume
end

function on_outer_box(bounds, half_box, depth, height, tolerance)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    return (abs(xmin + half_box) < tolerance &&
            abs(xmax + half_box) < tolerance) ||
           (abs(xmin - half_box) < tolerance &&
            abs(xmax - half_box) < tolerance) ||
           (abs(ymin + half_box) < tolerance &&
            abs(ymax + half_box) < tolerance) ||
           (abs(ymin - half_box) < tolerance &&
            abs(ymax - half_box) < tolerance) ||
           (abs(zmin + depth) < tolerance &&
            abs(zmax + depth) < tolerance) ||
           (abs(zmin - height) < tolerance &&
            abs(zmax - height) < tolerance)
end

function generate_corner_validation(;
    fabricated::Bool = false,
    aperture::Bool = false,
    island_width::Float64 = 8.0,
    island_height::Float64 = 6.0,
    corner_radius::Float64 = 0.0,
    half_box::Float64 = 12.0,
    substrate_depth::Float64 = 8.0,
    vacuum_height::Float64 = 8.0,
    metal_thickness::Float64 = 0.1,
    overetch_depth::Float64 = 0.05,
    sidewall_angle::Float64 = 80.0,
    top_rounding::Float64 = 0.01,
    trench_rounding::Float64 = 0.01,
    lc_fine::Float64 = fabricated ? 0.01 : 0.08,
    lc_far::Float64 = 1.5,
    mesh_order::Int = 1,
    maxwell_excitation::Bool = false,
    filename::String,
)
    half_x = 0.5 * island_width
    half_y = 0.5 * island_height
    half_x > 0.0 && half_y > 0.0 ||
        error("island dimensions must be positive")
    0.0 <= corner_radius < min(half_x, half_y) ||
        error("corner_radius must lie in [0, min(island_width, island_height) / 2)")
    half_box > max(half_x, half_y) ||
        error("half_box must contain the island")
    substrate_depth > overetch_depth ||
        error("substrate_depth must exceed overetch_depth")
    vacuum_height > metal_thickness ||
        error("vacuum_height must exceed metal_thickness")
    aperture && maxwell_excitation &&
        error("The Maxwell excitation is currently defined only for an island")

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add(
        "corner_validation_" * (aperture ? "aperture_" : "island_") *
        (fabricated ? "fabricated" : "thin"))
    occ = gmsh.model.occ

    substrate_seed = Tuple{Int32,Int32}[]
    vacuum_seed = Tuple{Int32,Int32}[]
    if fabricated
        angle = deg2rad(sidewall_angle)
        metal_pullback = metal_thickness / tan(angle)
        trench_pullback = overetch_depth / tan(angle)

        substrate_box = [(
            3,
            occ.addBox(
                -half_box, -half_box, -substrate_depth,
                2half_box, 2half_box, substrate_depth),
        )]
        notch = if aperture
            corner_radius > 0.0 ?
                tapered_rounded_rectangle(
                    occ, half_x, half_y, corner_radius,
                    -overetch_depth, overetch_depth, trench_pullback, 0.0) :
                tapered_rectangle(
                    occ, half_x, half_y, -overetch_depth, overetch_depth,
                    trench_pullback, 0.0)
        else
            trench_slab = [(
                3,
                occ.addBox(
                    -half_box, -half_box, -overetch_depth,
                    2half_box, 2half_box, overetch_depth),
            )]
            pedestal = corner_radius > 0.0 ?
                tapered_rounded_rectangle(
                    occ, half_x, half_y, corner_radius,
                    -overetch_depth, overetch_depth, -trench_pullback, 0.0) :
                tapered_rectangle(
                    occ, half_x, half_y, -overetch_depth, overetch_depth,
                    -trench_pullback, 0.0)
            occ.cut(trench_slab, pedestal)[1]
        end
        substrate, _ = occ.cut(substrate_box, notch)
        tolerance = 1.0e-6 * half_box
        substrate = fillet_volume_edges(
            occ, substrate, trench_rounding,
            bounds -> begin
                xmin, ymin, zmin, xmax, ymax, zmax = bounds
                on_floor = abs(zmin + overetch_depth) < tolerance &&
                           abs(zmax + overetch_depth) < tolerance
                on_pedestal =
                    maximum(abs.((xmin, ymin, xmax, ymax))) <
                    half_box - tolerance
                on_floor && on_pedestal
            end,
            "validation trench-bottom",
        )

        metal = if aperture
            metal_sheet = [(
                3,
                occ.addBox(
                    -half_box, -half_box, 0.0,
                    2half_box, 2half_box, metal_thickness),
            )]
            opening = corner_radius > 0.0 ?
                tapered_rounded_rectangle(
                    occ, half_x, half_y, corner_radius,
                    0.0, metal_thickness, 0.0, -metal_pullback) :
                tapered_rectangle(
                    occ, half_x, half_y, 0.0, metal_thickness,
                    0.0, -metal_pullback)
            occ.cut(metal_sheet, opening)[1]
        else
            corner_radius > 0.0 ?
                tapered_rounded_rectangle(
                    occ, half_x, half_y, corner_radius,
                    0.0, metal_thickness, 0.0, metal_pullback) :
                tapered_rectangle(
                    occ, half_x, half_y, 0.0, metal_thickness,
                    0.0, metal_pullback)
        end
        metal = fillet_volume_edges(
            occ, metal, top_rounding,
            bounds -> begin
                xmin, ymin, zmin, xmax, ymax, zmax = bounds
                on_top = abs(zmin - metal_thickness) < tolerance &&
                         abs(zmax - metal_thickness) < tolerance
                off_outer_box =
                    maximum(abs.((xmin, ymin, xmax, ymax))) <
                    half_box - tolerance
                on_top && (!aperture || off_outer_box)
            end,
            "validation top-metal",
        )
        outer = [(
            3,
            occ.addBox(
                -half_box, -half_box, -substrate_depth,
                2half_box, 2half_box, substrate_depth + vacuum_height),
        )]
        field, _ = occ.cut(outer, metal, -1, true, true)
        vacuum, _ = occ.cut(field, substrate, -1, true, false)
        domains, domain_map = occ.fragment(vcat(substrate, vacuum), [])
        substrate_seed = domain_map[1:length(substrate)] |> Iterators.flatten |> collect
        vacuum_seed =
            domain_map[length(substrate) + 1:end] |> Iterators.flatten |> collect
    else
        substrate_box = (3, occ.addBox(
            -half_box, -half_box, -substrate_depth,
            2half_box, 2half_box, substrate_depth))
        vacuum_box = (3, occ.addBox(
            -half_box, -half_box, 0.0,
            2half_box, 2half_box, vacuum_height))
        footprint =
            rounded_rectangle_surface(occ, half_x, half_y, corner_radius, 0.0)
        lower = occ.extrude(
            [(2, footprint)], 0.0, 0.0, -substrate_depth)
        upper = occ.extrude(
            [(2, footprint)], 0.0, 0.0, vacuum_height)
        tools = [
            entity for entity in vcat(lower, upper)
            if entity[1] == 3
        ]
        domains, domain_map =
            occ.fragment([substrate_box, vacuum_box], tools)
        append!(substrate_seed, domain_map[1])
        append!(vacuum_seed, domain_map[2])
    end

    lumped_port_seed = Tuple{Int32,Int32}[]
    sa_port_exclusion_seed = Tuple{Int32,Int32}[]
    port_width = 0.0
    if maxwell_excitation
        port_width = 0.5
        port = (
            2,
            occ.addRectangle(
                half_x, -0.5port_width, 0.0, half_box - half_x, port_width),
        )
        tools = [port]
        if fabricated
            trench_pullback =
                overetch_depth / tan(deg2rad(sidewall_angle))
            exclusion_start = half_x + trench_pullback
            push!(
                tools,
                (
                    2,
                    occ.addRectangle(
                        exclusion_start,
                        -0.5port_width,
                        -overetch_depth,
                        half_box - exclusion_start,
                        port_width,
                    ),
                ),
            )
        end
        old_domains = domains
        old_substrate = Set(
            tag for (dim, tag) in substrate_seed if dim == 3)
        old_vacuum = Set(
            tag for (dim, tag) in vacuum_seed if dim == 3)
        domains, domain_map =
            occ.fragment(old_domains, tools, -1, true, false)
        substrate_seed = Tuple{Int32,Int32}[]
        vacuum_seed = Tuple{Int32,Int32}[]
        for (index, (dim, tag)) in enumerate(old_domains)
            dim == 3 || continue
            tag in old_substrate && append!(substrate_seed, domain_map[index])
            tag in old_vacuum && append!(vacuum_seed, domain_map[index])
        end
        append!(lumped_port_seed, domain_map[length(old_domains) + 1])
        fabricated &&
            append!(
                sa_port_exclusion_seed,
                domain_map[length(old_domains) + 2],
            )
    end
    occ.synchronize()

    domain_tags = Set(tag for (dim, tag) in domains if dim == 3)
    substrate_tags = sort!(unique(
        tag for (dim, tag) in substrate_seed
        if dim == 3 && tag in domain_tags))
    vacuum_tags = sort!(unique(
        tag for (dim, tag) in vacuum_seed
        if dim == 3 && tag in domain_tags))
    substrate_set = Set(substrate_tags)
    vacuum_set = Set(vacuum_tags)
    lumped_port_set = Set(
        tag for (dim, tag) in lumped_port_seed if dim == 2)
    sa_port_exclusion_set = Set(
        tag for (dim, tag) in sa_port_exclusion_seed if dim == 2)
    tolerance = 1.0e-6 * half_box
    outer = Int32[]
    outer_truncation = Int32[]
    lumped_port = Int32[]
    sa_port_exclusion = Int32[]
    thin_metal = Int32[]
    ms = Int32[]
    sa = Int32[]
    ma = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        bounds = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bounds
        on_lumped_port =
            maxwell_excitation &&
            abs(zmin) < tolerance && abs(zmax) < tolerance &&
            xmin >= half_x - tolerance && xmax <= half_box + tolerance &&
            ymin >= -0.5port_width - tolerance &&
            ymax <= 0.5port_width + tolerance
        if tag in lumped_port_set || on_lumped_port
            push!(lumped_port, tag)
            continue
        elseif tag in sa_port_exclusion_set
            push!(sa_port_exclusion, tag)
            continue
        end
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        adjacent_substrate = [volume for volume in up if volume in substrate_set]
        adjacent_vacuum = [volume for volume in up if volume in vacuum_set]
        isempty(adjacent_substrate) && isempty(adjacent_vacuum) && continue
        if on_outer_box(
            bounds, half_box, substrate_depth, vacuum_height, tolerance)
            on_horizontal_outer =
                (abs(zmin + substrate_depth) < tolerance &&
                 abs(zmax + substrate_depth) < tolerance) ||
                (abs(zmin - vacuum_height) < tolerance &&
                 abs(zmax - vacuum_height) < tolerance)
            if aperture && !on_horizontal_outer
                push!(outer_truncation, tag)
            else
                push!(outer, tag)
            end
        elseif fabricated
            if !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                push!(sa, tag)
            elseif !isempty(adjacent_substrate)
                push!(ms, tag)
            elseif !isempty(adjacent_vacuum)
                push!(ma, tag)
            end
        else
            on_interface = abs(zmin) < tolerance && abs(zmax) < tolerance
            inside_footprint =
                xmin >= -half_x - tolerance &&
                xmax <= half_x + tolerance &&
                ymin >= -half_y - tolerance &&
                ymax <= half_y + tolerance
            on_thin_metal =
                on_interface && (aperture ? !inside_footprint : inside_footprint)
            if on_thin_metal
                push!(thin_metal, tag)
            elseif on_interface &&
                   !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                push!(sa, tag)
            end
        end
    end

    for surface in lumped_port
        up, _ = gmsh.model.getAdjacencies(2, surface)
        if isempty(up)
            length(vacuum_tags) == 1 ||
                error("Cannot identify the vacuum volume containing the lumped port")
            gmsh.model.mesh.embed(2, [surface], 3, only(vacuum_tags))
        end
    end

    groups = [
        (3, substrate_tags, 1, "substrate"),
        (3, vacuum_tags, 2, "vacuum"),
        (2, outer, 1, "outer"),
    ]
    aperture &&
        push!(
            groups,
            (2, outer_truncation, 7, "outer_truncation"),
        )
    maxwell_excitation &&
        push!(groups, (2, lumped_port, 5, "lumped_port"))
    fabricated && maxwell_excitation &&
        push!(
            groups,
            (2, sa_port_exclusion, 6, "SA_port_exclusion"),
        )
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

    feature_surfaces = fabricated ? vcat(ms, sa, ma) : thin_metal
    excluded_feature_curves = Set(
        curve for (dim, curve) in gmsh.model.getBoundary(
            [(2, surface) for surface in sa_port_exclusion],
            false,
            false,
            false,
        )
        if dim == 1
    )
    feature_curves = Int32[]
    for surface in feature_surfaces
        for (dim, curve) in gmsh.model.getBoundary(
            [(2, surface)], false, false, false)
            dim == 1 || continue
            curve in excluded_feature_curves && continue
            bounds = gmsh.model.getBoundingBox(dim, curve)
            on_outer_box(
                bounds, half_box, substrate_depth, vacuum_height, tolerance) ||
                push!(feature_curves, curve)
        end
    end
    unique!(feature_curves)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(feature_curves))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2lc_fine)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 2.0)
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

    println("Corner validation: geometry=$(aperture ? "aperture" : "island"), " *
            "fabricated=$(fabricated), " *
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

function main(args)
    length(args) >= 2 ||
        error("Usage: mesh_corner_validation.jl thin|fabricated OUTPUT.msh " *
              "[CORNER_RADIUS_UM] [island|aperture] [maxwell]")
    kind = args[1]
    kind in ("thin", "fabricated") || error("Unknown validation kind: $kind")
    corner_radius = 0.0
    radius_set = false
    aperture = false
    maxwell_excitation = false
    for argument in args[3:end]
        if argument == "island"
            aperture = false
        elseif argument == "aperture"
            aperture = true
        elseif argument == "maxwell"
            maxwell_excitation = true
        elseif !radius_set
            corner_radius = parse(Float64, argument)
            radius_set = true
        else
            error("Unknown validation argument: $argument")
        end
    end
    generate_corner_validation(
        fabricated = kind == "fabricated",
        aperture = aperture,
        corner_radius = corner_radius,
        maxwell_excitation = maxwell_excitation,
        filename = abspath(args[2]),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
