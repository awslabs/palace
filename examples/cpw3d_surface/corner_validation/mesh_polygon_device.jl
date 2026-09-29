# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabrication-resolved and thin-metal meshes of a polygonal metal island (convex corners) or a
# polygonal aperture in a metal sheet (concave corners) for the device-style check of the
# angle-interpolated corner family (USER decision 121 (C)): a device whose sharp corners lie at
# several NON-node interior angles, corrected on the thin mesh and compared with its fabricated
# reference (prepare_validation.py / compare_validation.py of this directory read the same
# physical groups as mesh_corner_validation.jl: 3D 1 substrate, 2 vacuum; 2D 1 outer,
# fabricated 2 MS / 3 SA / 4 MA, thin 2 thin_metal / 3 SA; aperture 7 outer_truncation).
#
# Polygon: a closed list of plan-view vertices (microns, counterclockwise) given on the command
# line as x1,y1 x2,y2 ..., or the named preset `heptagon` = the check's convex polygon with the
# interior angles 82 / 97 / 112 / 128 / 143 / 166 / 172 deg (one per interval of the family's
# nodes 75 / 90 / 105 / 120 / 135 / 150 / 165 / 180, the last two in the first-order regime).
# Fabrication: metal thickness, overetch depth, sidewall angle (90 = vertical), no rounding
# (the R 1.9 library's Fabrication record); a sloped sidewall offsets the polygon inward by
# thickness / tan(angle) at the top (loft between the two offset polygons).
#
#   julia mesh_polygon_device.jl thin|fabricated OUTPUT.msh island|aperture heptagon|x1,y1 x2,y2 ...
#       [--half-box 24] [--lc-fine 0.01] [--lc-far 2.0] [--order 1] [--metal 0.1] [--overetch 0.05]
#       [--sidewall 90]

import Gmsh: gmsh

function heptagon_points()
    # Interior angles 172 / 166 / 143 / 128 / 112 / 97 / 82 (sum 900 = 5 x 180): exterior turns
    # 8 / 14 / 37 / 52 / 68 / 83 / 98 deg (sum 360), the two largest last so that the closing
    # pair of edges is well conditioned. Five edge lengths are chosen (>= 8 um so the corners
    # are isolated events at R 1.9) and the last two close the polygon (24.3 / 21.9 um).
    turns = deg2rad.([8.0, 14.0, 37.0, 52.0, 68.0, 83.0, 98.0])
    chosen = [10.0, 9.0, 8.0, 9.0, 10.0]
    headings = cumsum(vcat(0.0, turns[1:(end - 1)]))
    directions = [(cos(h), sin(h)) for h in headings]
    # Closure: sum_i L_i d_i = 0 with L_6, L_7 unknown.
    rx = -sum(chosen[i] * directions[i][1] for i in 1:5)
    ry = -sum(chosen[i] * directions[i][2] for i in 1:5)
    a, b = directions[6], directions[7]
    det = a[1] * b[2] - a[2] * b[1]
    l6 = (rx * b[2] - ry * b[1]) / det
    l7 = (a[1] * ry - a[2] * rx) / det
    lengths = vcat(chosen, [l6, l7])
    all(lengths .> 0.0) || error("heptagon closure gives a non-positive edge length: $lengths")
    points = Tuple{Float64, Float64}[]
    x, y = 0.0, 0.0
    for i in 1:7
        push!(points, (x, y))
        x += lengths[i] * directions[i][1]
        y += lengths[i] * directions[i][2]
    end
    abs(x) < 1e-9 && abs(y) < 1e-9 || error("heptagon does not close: ($x, $y)")
    cx = sum(p[1] for p in points) / 7
    cy = sum(p[2] for p in points) / 7
    return [(p[1] - cx, p[2] - cy) for p in points]
end

function interior_angles_degrees(points)
    n = length(points)
    angles = Float64[]
    for i in 1:n
        p0, p1, p2 = points[mod1(i - 1, n)], points[i], points[mod1(i + 1, n)]
        u = (p0[1] - p1[1], p0[2] - p1[2])
        v = (p2[1] - p1[1], p2[2] - p1[2])
        c = (u[1] * v[1] + u[2] * v[2]) / (hypot(u...) * hypot(v...))
        cross = u[1] * v[2] - u[2] * v[1]
        a = acosd(clamp(c, -1.0, 1.0))
        # Counterclockwise polygon: a left turn (cross of the outgoing edges > 0) is convex.
        push!(angles, cross < 0.0 ? a : 360.0 - a)
    end
    return angles
end

function signed_area(points)
    n = length(points)
    return 0.5 * sum(
        points[i][1] * points[mod1(i + 1, n)][2] - points[mod1(i + 1, n)][1] * points[i][2] for
        i in 1:n
    )
end

# Inward offset of a simple polygon by `distance` (the intersection of the offset edge lines
# at every vertex; exact for offsets small against the edges, convex and concave vertices
# alike). A negative distance offsets outward.
function offset_polygon(points, distance)
    n = length(points)
    out = Tuple{Float64, Float64}[]
    for i in 1:n
        p0, p1, p2 = points[mod1(i - 1, n)], points[i], points[mod1(i + 1, n)]
        d0 = (p1[1] - p0[1], p1[2] - p0[2])
        d1 = (p2[1] - p1[1], p2[2] - p1[2])
        l0, l1 = hypot(d0...), hypot(d1...)
        # Left normals (inward for a counterclockwise polygon).
        n0 = (-d0[2] / l0, d0[1] / l0)
        n1 = (-d1[2] / l1, d1[1] / l1)
        # Offset lines: (x - p1) . n0 = distance and (x - p1) . n1 = distance.
        det = n0[1] * n1[2] - n0[2] * n1[1]
        if abs(det) < 1e-12
            push!(out, (p1[1] + distance * n0[1], p1[2] + distance * n0[2]))
        else
            c0 = n0[1] * p1[1] + n0[2] * p1[2] + distance
            c1 = n1[1] * p1[1] + n1[2] * p1[2] + distance
            push!(out, ((c0 * n1[2] - c1 * n0[2]) / det, (n0[1] * c1 - n1[1] * c0) / det))
        end
    end
    return out
end

function polygon_wire(occ, points, z)
    tags = [occ.addPoint(p[1], p[2], z) for p in points]
    lines = [occ.addLine(tags[i], tags[mod1(i + 1, length(tags))]) for i in eachindex(tags)]
    return occ.addCurveLoop(lines)
end

# Prism between the polygon offset inward by pullback_bottom at z0 and by pullback_top at
# z0 + height (a loft; a straight extrusion when both pullbacks agree).
function tapered_polygon(occ, points, z0, height, pullback_bottom, pullback_top)
    if pullback_bottom == pullback_top
        surface = occ.addPlaneSurface([polygon_wire(occ, offset_polygon(points, pullback_bottom), z0)])
        extruded = occ.extrude([(2, surface)], 0.0, 0.0, height)
        return [entity for entity in extruded if entity[1] == 3]
    end
    bottom = polygon_wire(occ, offset_polygon(points, pullback_bottom), z0)
    top = polygon_wire(occ, offset_polygon(points, pullback_top), z0 + height)
    lofted = occ.addThruSections([bottom, top], -1, true, true)
    return [entity for entity in lofted if entity[1] == 3]
end

function generate_polygon_device(;
    points::Vector{Tuple{Float64, Float64}},
    fabricated::Bool,
    aperture::Bool,
    half_box::Float64        = 24.0,
    substrate_depth::Float64 = 8.0,
    vacuum_height::Float64   = 8.0,
    metal_thickness::Float64 = 0.1,
    overetch_depth::Float64  = 0.05,
    sidewall_angle::Float64  = 90.0,
    lc_fine::Float64         = fabricated ? 0.01 : 0.08,
    lc_far::Float64          = 2.0,
    mesh_order::Int          = 1,
    filename::String
)
    polygon_area = signed_area(points)
    polygon_area > 0.0 || error("polygon vertices must be counterclockwise")
    extent = maximum(max(abs(p[1]), abs(p[2])) for p in points)
    half_box > extent + 2.0 || error("half_box must contain the polygon with a margin")
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(
        "polygon_device_" * (aperture ? "aperture_" : "island_") *
        (fabricated ? "fabricated" : "thin")
    )
    occ = gmsh.model.occ
    tolerance = 1.0e-6 * half_box
    substrate_seed = Tuple{Int32, Int32}[]
    vacuum_seed = Tuple{Int32, Int32}[]
    if fabricated
        angle = deg2rad(sidewall_angle)
        metal_pullback = metal_thickness / tan(angle)
        trench_pullback = overetch_depth / tan(angle)
        substrate_box = [(
            3,
            occ.addBox(-half_box, -half_box, -substrate_depth, 2half_box, 2half_box, substrate_depth)
        )]
        # Overetch: the substrate is removed by overetch_depth wherever there is no metal (the
        # metal footprint protects a pedestal under the island; the aperture is the notch).
        notch = if aperture
            tapered_polygon(occ, points, -overetch_depth, overetch_depth, trench_pullback, 0.0)
        else
            trench_slab = [(
                3,
                occ.addBox(-half_box, -half_box, -overetch_depth, 2half_box, 2half_box, overetch_depth)
            )]
            pedestal =
                tapered_polygon(occ, points, -overetch_depth, overetch_depth, -trench_pullback, 0.0)
            occ.cut(trench_slab, pedestal)[1]
        end
        substrate, _ = occ.cut(substrate_box, notch)
        metal = if aperture
            sheet = [(
                3,
                occ.addBox(-half_box, -half_box, 0.0, 2half_box, 2half_box, metal_thickness)
            )]
            opening = tapered_polygon(occ, points, 0.0, metal_thickness, 0.0, -metal_pullback)
            occ.cut(sheet, opening)[1]
        else
            tapered_polygon(occ, points, 0.0, metal_thickness, 0.0, metal_pullback)
        end
        outer = [(
            3,
            occ.addBox(
                -half_box,
                -half_box,
                -substrate_depth,
                2half_box,
                2half_box,
                substrate_depth + vacuum_height
            )
        )]
        field, _ = occ.cut(outer, metal, -1, true, true)
        vacuum, _ = occ.cut(field, substrate, -1, true, false)
        domains, domain_map = occ.fragment(vcat(substrate, vacuum), [])
        substrate_seed = domain_map[1:length(substrate)] |> Iterators.flatten |> collect
        vacuum_seed = domain_map[(length(substrate) + 1):end] |> Iterators.flatten |> collect
    else
        substrate_box =
            (3, occ.addBox(-half_box, -half_box, -substrate_depth, 2half_box, 2half_box, substrate_depth))
        vacuum_box =
            (3, occ.addBox(-half_box, -half_box, 0.0, 2half_box, 2half_box, vacuum_height))
        footprint = occ.addPlaneSurface([polygon_wire(occ, points, 0.0)])
        lower = occ.extrude([(2, footprint)], 0.0, 0.0, -substrate_depth)
        upper = occ.extrude([(2, footprint)], 0.0, 0.0, vacuum_height)
        tools = [entity for entity in vcat(lower, upper) if entity[1] == 3]
        domains, domain_map = occ.fragment([substrate_box, vacuum_box], tools)
        append!(substrate_seed, domain_map[1])
        append!(vacuum_seed, domain_map[2])
    end
    occ.synchronize()

    domain_tags = Set(tag for (dim, tag) in domains if dim == 3)
    substrate_tags =
        sort!(unique(tag for (dim, tag) in substrate_seed if dim == 3 && tag in domain_tags))
    vacuum_tags =
        sort!(unique(tag for (dim, tag) in vacuum_seed if dim == 3 && tag in domain_tags))
    substrate_set = Set(substrate_tags)
    vacuum_set = Set(vacuum_tags)
    on_outer_box(bounds) = begin
        xmin, ymin, zmin, xmax, ymax, zmax = bounds
        (abs(xmin + half_box) < tolerance && abs(xmax + half_box) < tolerance) ||
        (abs(xmin - half_box) < tolerance && abs(xmax - half_box) < tolerance) ||
        (abs(ymin + half_box) < tolerance && abs(ymax + half_box) < tolerance) ||
        (abs(ymin - half_box) < tolerance && abs(ymax - half_box) < tolerance) ||
        (abs(zmin + substrate_depth) < tolerance && abs(zmax + substrate_depth) < tolerance) ||
        (abs(zmin - vacuum_height) < tolerance && abs(zmax - vacuum_height) < tolerance)
    end
    outer = Int32[]
    outer_truncation = Int32[]
    thin_metal = Int32[]
    ms = Int32[]
    sa = Int32[]
    ma = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        bounds = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bounds
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        adjacent_substrate = [volume for volume in up if volume in substrate_set]
        adjacent_vacuum = [volume for volume in up if volume in vacuum_set]
        isempty(adjacent_substrate) && isempty(adjacent_vacuum) && continue
        if on_outer_box(bounds)
            on_horizontal_outer =
                (abs(zmin + substrate_depth) < tolerance && abs(zmax + substrate_depth) < tolerance) ||
                (abs(zmin - vacuum_height) < tolerance && abs(zmax - vacuum_height) < tolerance)
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
            on_interface || continue
            # The footprint surface is the one interface surface with the polygon's area (the
            # complement's centroid may lie inside the polygon, so the area decides).
            inside_footprint = abs(gmsh.model.occ.getMass(dim, tag) - polygon_area) < 1e-6 * polygon_area
            if aperture ? !inside_footprint : inside_footprint
                push!(thin_metal, tag)
            elseif !isempty(adjacent_substrate) && !isempty(adjacent_vacuum)
                push!(sa, tag)
            end
        end
    end

    groups = [(3, substrate_tags, 1, "substrate"), (3, vacuum_tags, 2, "vacuum"), (2, outer, 1, "outer")]
    aperture && push!(groups, (2, outer_truncation, 7, "outer_truncation"))
    if fabricated
        append!(groups, [(2, ms, 2, "MS"), (2, sa, 3, "SA"), (2, ma, 4, "MA")])
    else
        append!(groups, [(2, thin_metal, 2, "thin_metal"), (2, sa, 3, "SA")])
    end
    for (dim, entities, tag, name) in groups
        isempty(entities) && error("Empty physical group: $name")
        gmsh.model.addPhysicalGroup(dim, entities, tag, name)
    end

    # Mesh size: lc_fine within 2 lc_fine of the metal edges (the perimeter curves of the metal
    # surfaces off the outer box), lc_far beyond 2 um.
    feature_surfaces = fabricated ? vcat(ms, sa, ma) : thin_metal
    feature_curves = Int32[]
    for surface in feature_surfaces
        for (dim, curve) in gmsh.model.getBoundary([(2, surface)], false, false, false)
            dim == 1 || continue
            on_outer_box(gmsh.model.getBoundingBox(dim, curve)) || push!(feature_curves, curve)
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
        ("Mesh.Binary", 1)
    ]
        gmsh.option.setNumber(name, value)
    end
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(mesh_order)
    mesh_order > 1 && gmsh.model.mesh.optimize("HighOrderElastic")
    gmsh.write(filename)

    angles = interior_angles_degrees(points)
    println(
        "Polygon device: geometry=$(aperture ? "aperture" : "island"), fabricated=$(fabricated), " *
        "vertices=$(length(points)), interior angles (deg)=$(round.(angles; digits=3)), " *
        "metal=$(metal_thickness) overetch=$(overetch_depth) sidewall=$(sidewall_angle)"
    )
    for (dim, tag) in gmsh.model.getPhysicalGroups()
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  dim=$dim tag=$tag name=$name ($(length(entities)) entities)")
    end
    println("  nodes=$(length(gmsh.model.mesh.getNodes()[1]))")
    println("  file=$(filename)")
    return gmsh.finalize()
end

function main(args)
    length(args) >= 4 || error(
        "Usage: mesh_polygon_device.jl thin|fabricated OUTPUT.msh island|aperture " *
        "heptagon|x1,y1 x2,y2 ... [--half-box H] [--lc-fine L] [--lc-far L] [--order N] " *
        "[--metal T] [--overetch D] [--sidewall A]"
    )
    kind = args[1]
    kind in ("thin", "fabricated") || error("Unknown kind: $kind")
    geometry = args[3]
    geometry in ("island", "aperture") || error("Unknown geometry: $geometry")
    points = Tuple{Float64, Float64}[]
    options = Dict{String, Float64}()
    i = 4
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            i < length(args) || error("missing value for $a")
            options[a[3:end]] = parse(Float64, args[i + 1])
            i += 2
        elseif a == "heptagon"
            append!(points, heptagon_points())
            i += 1
        else
            x, y = split(a, ",")
            push!(points, (parse(Float64, x), parse(Float64, y)))
            i += 1
        end
    end
    length(points) >= 3 || error("a polygon needs at least three vertices")
    fabricated = kind == "fabricated"
    return generate_polygon_device(
        points=points,
        fabricated=fabricated,
        aperture=geometry == "aperture",
        half_box=get(options, "half-box", 24.0),
        lc_fine=get(options, "lc-fine", fabricated ? 0.01 : 0.08),
        lc_far=get(options, "lc-far", 2.0),
        mesh_order=Int(get(options, "order", 1.0)),
        metal_thickness=get(options, "metal", 0.1),
        overetch_depth=get(options, "overetch", 0.05),
        sidewall_angle=get(options, "sidewall", 90.0),
        filename=abspath(args[2])
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
