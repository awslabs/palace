# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Two-dimensional local coupon for a pair of nearby edges bounding a cutout.
# Dimensions are in microns. The two metal regions are equipotential, as they
# are for the two sides of a ground-plane cutout.

import Gmsh: gmsh

function fillet(c, p, n, radius)
    radius <= 0 && return (c, c, c)
    d1 = let v = (c[1] - p[1], c[2] - p[2])
        length = hypot(v...)
        (v[1] / length, v[2] / length)
    end
    d2 = let v = (n[1] - c[1], n[2] - c[2])
        length = hypot(v...)
        (v[1] / length, v[2] / length)
    end
    cosine = clamp(d1[1] * d2[1] + d1[2] * d2[2], -1.0, 1.0)
    tangent = radius / tan((pi - acos(cosine)) / 2)
    point_in = (c[1] - tangent * d1[1], c[2] - tangent * d1[2])
    point_out = (c[1] + tangent * d2[1], c[2] + tangent * d2[2])
    normal = (-d1[2], d1[1])
    center = (point_in[1] + radius * normal[1],
              point_in[2] + radius * normal[2])
    return (point_in, center, point_out)
end

function shape(occ, corners)
    count = length(corners)
    data = [
        fillet(corners[i][1], corners[mod1(i - 1, count)][1],
               corners[mod1(i + 1, count)][1], corners[i][2])
        for i in 1:count
    ]
    curves = Int32[]
    for i in 1:count
        j = mod1(i + 1, count)
        p1 = occ.addPoint(data[i][3][1], data[i][3][2], 0.0)
        p2 = occ.addPoint(data[j][1][1], data[j][1][2], 0.0)
        push!(curves, occ.addLine(p1, p2))
        if corners[j][2] > 0
            center = occ.addPoint(data[j][2][1], data[j][2][2], 0.0)
            p3 = occ.addPoint(data[j][3][1], data[j][3][2], 0.0)
            push!(curves, occ.addCircleArc(p2, center, p3))
        end
    end
    return occ.addPlaneSurface([occ.addCurveLoop(curves)])
end

function generate_edge_pair_coupon(;
    cutout_width::Float64,
    radius::Float64 = 2.0,
    fabricated::Bool = false,
    t_metal::Float64 = 0.1,
    overetch::Float64 = 0.05,
    sidewall_angle::Float64 = 80.0,
    r_top::Float64 = 0.01,
    r_bottom::Float64 = 0.01,
    lc_fine::Float64 = 0.002,
    lc_far::Float64 = 0.05,
    mesh_order::Int = 2,
    filename::String,
)
    cutout_width > 0 || error("cutout_width must be positive")
    radius > 0 || error("radius must be positive")
    half_gap = 0.5 * cutout_width
    half_width = half_gap + radius
    tolerance = 1.0e-5

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("edge_pair_coupon")
    occ = gmsh.model.occ

    if fabricated
        angle = deg2rad(sidewall_angle)
        metal_pullback = t_metal / tan(angle)
        trench_pullback = overetch / tan(angle)
        extension = 0.2

        outer = occ.addRectangle(-half_width, -radius, 0.0,
                                 2half_width, 2radius)
        left_metal = shape(occ, [
            ((-half_width - extension, 0.0), 0.0),
            ((-half_gap, 0.0), 0.0),
            ((-half_gap - metal_pullback, t_metal), r_top),
            ((-half_width - extension, t_metal), 0.0),
        ])
        right_metal = shape(occ, [
            ((half_gap, 0.0), 0.0),
            ((half_width + extension, 0.0), 0.0),
            ((half_width + extension, t_metal), 0.0),
            ((half_gap + metal_pullback, t_metal), r_top),
        ])
        substrate_base =
            occ.addRectangle(-half_width, -radius, 0.0, 2half_width, radius)
        trench = shape(occ, [
            ((-half_gap, 0.0), 0.0),
            ((-half_gap + trench_pullback, -overetch), r_bottom),
            ((half_gap - trench_pullback, -overetch), r_bottom),
            ((half_gap, 0.0), 0.0),
        ])
        occ.synchronize()

        substrate, _ = occ.cut([(2, substrate_base)], [(2, trench)])
        field, _ = occ.cut([(2, outer)], [(2, left_metal), (2, right_metal)])
        vacuum, _ = occ.cut(field, substrate, -1, true, false)
        occ.fragment(vcat(vacuum, substrate), [])
    else
        x_coordinates = (-half_width, -half_gap, half_gap, half_width)
        surfaces = Tuple{Int32,Int32}[]
        for i in 1:3
            width = x_coordinates[i + 1] - x_coordinates[i]
            push!(surfaces,
                  (2, occ.addRectangle(x_coordinates[i], -radius, 0.0,
                                       width, radius)))
            push!(surfaces,
                  (2, occ.addRectangle(x_coordinates[i], 0.0, 0.0,
                                       width, radius)))
        end
        occ.fragment(surfaces, [])
    end
    occ.synchronize()

    substrate_surfaces = Int32[]
    vacuum_surfaces = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        _, y, _ = occ.getCenterOfMass(dim, tag)
        push!(y < (fabricated ? -0.1 * overetch : 0.0) ?
              substrate_surfaces : vacuum_surfaces, tag)
    end

    outer_curves = Int32[]
    ms_curves = Int32[]
    ma_horizontal = Int32[]
    ma_side = Int32[]
    sa_side = Int32[]
    sa_floor = Int32[]
    model_surfaces = Set(vcat(substrate_surfaces, vacuum_surfaces))
    for (dim, tag) in gmsh.model.getEntities(1)
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        isempty([surface for surface in up if surface in model_surfaces]) && continue
        xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(dim, tag)
        xmid = 0.5 * (xmin + xmax)
        ymid = 0.5 * (ymin + ymax)
        horizontal = ymax - ymin < tolerance
        vertical = xmax - xmin < tolerance
        on_outer =
            (vertical && (abs(xmid + half_width) < tolerance ||
                          abs(xmid - half_width) < tolerance)) ||
            (horizontal && (abs(ymid + radius) < tolerance ||
                            abs(ymid - radius) < tolerance))
        if on_outer
            push!(outer_curves, tag)
        elseif !fabricated && horizontal && abs(ymid) < tolerance
            if xmax <= -half_gap + tolerance || xmin >= half_gap - tolerance
                push!(ms_curves, tag)
            elseif xmin >= -half_gap - tolerance && xmax <= half_gap + tolerance
                push!(sa_floor, tag)
            end
        elseif fabricated && horizontal && abs(ymid) < tolerance &&
               (xmax <= -half_gap + tolerance || xmin >= half_gap - tolerance)
            push!(ms_curves, tag)
        elseif fabricated && horizontal && abs(ymid - t_metal) < tolerance
            push!(ma_horizontal, tag)
        elseif fabricated && ymin >= -tolerance &&
               ymax <= t_metal + tolerance && !horizontal
            push!(ma_side, tag)
        elseif fabricated && ymin >= -overetch - tolerance &&
               ymax <= tolerance && !horizontal
            push!(sa_side, tag)
        elseif fabricated && horizontal &&
               abs(ymid + overetch) < tolerance &&
               xmin >= -half_gap - tolerance && xmax <= half_gap + tolerance
            push!(sa_floor, tag)
        end
    end

    groups = [
        (2, substrate_surfaces, 1, "substrate"),
        (2, vacuum_surfaces, 2, "vacuum"),
        (1, outer_curves, 1, "matching_contour"),
        (1, ms_curves, 2, fabricated ? "MS" : "thin_metal"),
        (1, sa_floor, 3, fabricated ? "SA_floor" : "SA"),
    ]
    if fabricated
        append!(groups, [
            (1, ma_horizontal, 4, "MA_horizontal"),
            (1, ma_side, 5, "MA_side"),
            (1, sa_side, 6, "SA_side"),
        ])
    end
    for (dim, entities, tag, name) in groups
        isempty(entities) && error("Empty physical group: $name")
        gmsh.model.addPhysicalGroup(dim, entities, tag, name)
    end

    features = fabricated ?
        vcat(ms_curves, ma_horizontal, ma_side, sa_side, sa_floor) :
        vcat(ms_curves, sa_floor)
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(features))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.02)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    for (name, value) in [
        ("Mesh.MeshSizeMin", lc_fine),
        ("Mesh.MeshSizeMax", lc_far),
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
    ]
        gmsh.option.setNumber(name, value)
    end

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.setOrder(mesh_order)
    mesh_order > 1 && gmsh.model.mesh.optimize("HighOrderElastic")
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    println("Pair coupon: cutout=$(cutout_width) um, R=$(radius) um, " *
            "fabricated=$(fabricated)")
    for (dim, tag) in gmsh.model.getPhysicalGroups()
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  dim=$dim tag=$tag name=$name ($(length(entities)) entities)")
    end
    println("  nodes=$(length(gmsh.model.mesh.getNodes()[1]))")
    println("  file=$filename")
    gmsh.finalize()
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 3 ||
        error("Usage: mesh_edge_pair_coupon.jl thin|fabricated CUTOUT_WIDTH OUTPUT.msh")
    kind = ARGS[1]
    kind in ("thin", "fabricated") || error("Unknown coupon kind: $kind")
    generate_edge_pair_coupon(
        cutout_width = parse(Float64, ARGS[2]),
        fabricated = kind == "fabricated",
        filename = abspath(ARGS[3]),
    )
end
