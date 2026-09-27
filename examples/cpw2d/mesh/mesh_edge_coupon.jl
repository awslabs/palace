# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Two-dimensional local coupon for one isolated thin-metal edge. Dimensions
# are in microns. The edge is at the origin, metal extends toward negative x,
# substrate occupies negative y, and the matching contour is a square of
# half-width R.
#
# Curved (axisymmetric) coupon: with `axisymmetric_radius` rho > 0 the same
# cross-section is placed in the (r, z) half-plane (x = r) with the edge at
# r = rho, for Palace `Model.Axisymmetric`: `convexity = :convex` keeps the
# metal inside the circle (r < rho, a disk edge), `:concave` mirrors the
# cross-section so the metal lies outside (r > rho, a hole edge). The contour
# is [rho - R, rho + R] x [-R, R]; every physical group keeps its meaning.

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

function generate_edge_coupon(;
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
    axisymmetric_radius::Float64 = 0.0,
    convexity::Symbol = :convex,
    filename::String,
)
    radius > 0 || error("radius must be positive")
    axisymmetric_radius >= 0 || error("axisymmetric_radius must be nonnegative")
    axisymmetric_radius == 0 || axisymmetric_radius > radius ||
        error("axisymmetric_radius must exceed the coupon radius (the contour must not cross the axis)")
    convexity in (:convex, :concave) || error("convexity must be :convex or :concave")
    # Map from the canonical cross-section (edge at x = 0, metal at x < 0) to
    # the mesh: X(x) = rho + sign x; the inverse U(x) classifies mesh entities.
    sign = convexity == :convex ? 1.0 : -1.0
    rho = axisymmetric_radius
    X(x) = rho + sign * x
    U(x) = sign * (x - rho)
    function rect(x, y, dx, dy)
        x1, x2 = minmax(X(x), X(x + dx))
        return occ.addRectangle(x1, y, 0.0, x2 - x1, dy)
    end
    function polygon(corners)
        mapped = [((X(c[1][1]), c[1][2]), c[2]) for c in corners]
        return shape(occ, sign > 0 ? mapped : reverse(mapped))
    end
    t_metal > 0 || error("t_metal must be positive")
    overetch > 0 || error("overetch must be positive")
    0 < sidewall_angle <= 90 ||
        error("sidewall_angle must be in (0, 90] degrees")
    r_top >= 0 || error("r_top must be nonnegative")
    r_bottom >= 0 || error("r_bottom must be nonnegative")
    lc_fine > 0 || error("lc_fine must be positive")
    lc_far >= lc_fine || error("lc_far must be at least lc_fine")
    mesh_order > 0 || error("mesh_order must be positive")
    tolerance = 1.0e-5

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("edge_coupon")
    occ = gmsh.model.occ

    if fabricated
        angle = deg2rad(sidewall_angle)
        metal_pullback = t_metal / tan(angle)
        trench_pullback = overetch / tan(angle)
        extension = 0.2
        outer = rect(-radius, -radius, 2radius, 2radius)
        metal = polygon([
            ((-radius - extension, 0.0), 0.0),
            ((0.0, 0.0), 0.0),
            ((-metal_pullback, t_metal), r_top),
            ((-radius - extension, t_metal), 0.0),
        ])
        substrate_base = rect(-radius, -radius, 2radius, radius)
        trench = polygon([
            ((0.0, 0.0), 0.0),
            ((trench_pullback, -overetch), r_bottom),
            ((radius + extension, -overetch), 0.0),
            ((radius + extension, 0.0), 0.0),
        ])
        occ.synchronize()
        substrate, _ = occ.cut([(2, substrate_base)], [(2, trench)])
        field, _ = occ.cut([(2, outer)], [(2, metal)])
        vacuum, _ = occ.cut(field, substrate, -1, true, false)
        occ.fragment(vcat(vacuum, substrate), [])
    else
        surfaces = [
            rect(-radius, -radius, radius, radius),
            rect(0.0, -radius, radius, radius),
            rect(-radius, 0.0, radius, radius),
            rect(0.0, 0.0, radius, radius),
        ]
        occ.fragment([(2, surface) for surface in surfaces], [])
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
        xmid = U(0.5 * (xmin + xmax))
        ymid = 0.5 * (ymin + ymax)
        horizontal = ymax - ymin < tolerance
        vertical = xmax - xmin < tolerance
        on_outer =
            (vertical && (abs(xmid + radius) < tolerance ||
                          abs(xmid - radius) < tolerance)) ||
            (horizontal && (abs(ymid + radius) < tolerance ||
                            abs(ymid - radius) < tolerance))
        if on_outer
            push!(outer_curves, tag)
        elseif !fabricated && horizontal && abs(ymid) < tolerance
            push!(xmid < 0.0 ? ms_curves : sa_floor, tag)
        elseif fabricated && horizontal && abs(ymid) < tolerance && xmid < 0.0
            push!(ms_curves, tag)
        elseif fabricated && horizontal && abs(ymid - t_metal) < tolerance
            push!(ma_horizontal, tag)
        elseif fabricated && ymin >= -tolerance &&
               ymax <= t_metal + tolerance && !horizontal
            push!(ma_side, tag)
        elseif fabricated && ymin >= -overetch - tolerance &&
               ymax <= tolerance && !horizontal
            push!(sa_side, tag)
        elseif fabricated && horizontal && abs(ymid + overetch) < tolerance &&
               xmid > 0.0
            push!(sa_floor, tag)
        end
    end

    groups = [
        (2, substrate_surfaces, 1, "substrate"),
        (2, vacuum_surfaces, 2, "vacuum"),
        (1, outer_curves, 1, "matching_contour"),
        (1, ms_curves, 2, fabricated ? "MS" : "thin_metal"),
        (1, sa_floor, fabricated ? 6 : 3, fabricated ? "SA_floor" : "SA"),
    ]
    if fabricated
        append!(groups, [
            (1, ma_horizontal, 3, "MA_horizontal"),
            (1, ma_side, 4, "MA_side"),
            (1, sa_side, 5, "SA_side"),
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
    gmsh.model.mesh.field.setNumber(2, "DistMax", min(0.5, 0.5 * radius))
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

    println("Edge coupon: R=$(radius) um, fabricated=$(fabricated)" *
            (rho > 0 ? ", axisymmetric rho=$(rho) um ($(convexity))" : ""))
    for (dim, tag) in gmsh.model.getPhysicalGroups()
        name = gmsh.model.getPhysicalName(dim, tag)
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  dim=$dim tag=$tag name=$name ($(length(entities)) entities)")
    end
    println("  nodes=$(length(gmsh.model.mesh.getNodes()[1]))")
    println("  file=$filename")
    gmsh.finalize()
end

function main(args)
    length(args) >= 2 ||
        error("Usage: mesh_edge_coupon.jl thin|fabricated OUTPUT.msh " *
              "[--radius R] [--metal-thickness T] [--overetch D] " *
              "[--sidewall-angle A] [--top-radius R] [--bottom-radius R] " *
              "[--lc-fine H] [--lc-far H] [--mesh-order P] " *
              "[--axisymmetric-radius RHO] [--convexity convex|concave]")
    kind = args[1]
    kind in ("thin", "fabricated") || error("Unknown coupon kind: $kind")
    options = Dict{String,String}()
    index = 3
    while index <= length(args)
        name = args[index]
        startswith(name, "--") || error("Expected an option, found: $name")
        index < length(args) || error("Missing value for option: $name")
        haskey(options, name) && error("Repeated option: $name")
        options[name] = args[index + 1]
        index += 2
    end
    allowed = Set([
        "--radius",
        "--metal-thickness",
        "--overetch",
        "--sidewall-angle",
        "--top-radius",
        "--bottom-radius",
        "--lc-fine",
        "--lc-far",
        "--mesh-order",
        "--axisymmetric-radius",
        "--convexity",
    ])
    unknown = setdiff(Set(keys(options)), allowed)
    isempty(unknown) || error("Unknown option(s): $(join(sort(collect(unknown)), ", "))")
    float_option(name, default) =
        haskey(options, name) ? parse(Float64, options[name]) : default
    int_option(name, default) =
        haskey(options, name) ? parse(Int, options[name]) : default

    generate_edge_coupon(
        radius = float_option("--radius", 2.0),
        fabricated = kind == "fabricated",
        t_metal = float_option("--metal-thickness", 0.1),
        overetch = float_option("--overetch", 0.05),
        sidewall_angle = float_option("--sidewall-angle", 80.0),
        r_top = float_option("--top-radius", 0.01),
        r_bottom = float_option("--bottom-radius", 0.01),
        lc_fine = float_option("--lc-fine", 0.002),
        lc_far = float_option("--lc-far", 0.05),
        mesh_order = int_option("--mesh-order", 2),
        axisymmetric_radius = float_option("--axisymmetric-radius", 0.0),
        convexity = Symbol(get(options, "--convexity", "convex")),
        filename = abspath(args[2]),
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
