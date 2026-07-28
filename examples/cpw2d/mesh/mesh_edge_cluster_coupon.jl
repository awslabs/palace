# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generic two-dimensional coupon for three or more parallel metal edges.

import Gmsh: gmsh
using DelimitedFiles

function read_edges(path)
    data, header = readdlm(path, ',', header = true)
    names = vec(String.(header))
    columns = Dict(name => index for (index, name) in enumerate(names))
    all(haskey(columns, name) for name in
        ("Offset", "GapDirection", "Conductor")) ||
        error("Cluster CSV must contain Offset, GapDirection, and Conductor")
    return [
        Dict(
            "Offset" => Float64(data[row, columns["Offset"]]),
            "GapDirection" => Int(data[row, columns["GapDirection"]]),
            "Conductor" => Int(data[row, columns["Conductor"]]),
        )
        for row in axes(data, 1)
    ]
end

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

function classify_intervals(edges, radius)
    length(edges) >= 3 || error("A parallel cluster requires at least three edges")
    offsets = Float64[edge["Offset"] for edge in edges]
    directions = Int[edge["GapDirection"] for edge in edges]
    conductors = Int[edge["Conductor"] for edge in edges]
    offsets[1] == 0.0 || error("The first edge offset must be zero")
    issorted(offsets) && all(diff(offsets) .> 0.0) ||
        error("Edge offsets must be strictly increasing")
    all(abs.(directions) .== 1) ||
        error("GapDirection must be +1 or -1")

    bounds = vcat(offsets[1] - radius, offsets, offsets[end] + radius)
    intervals = NamedTuple[]
    for i in 1:length(bounds)-1
        metal = if i == 1
            directions[1] == 1
        elseif i == length(bounds) - 1
            directions[end] == -1
        else
            right_of_left = directions[i - 1] == -1
            left_of_right = directions[i] == 1
            right_of_left == left_of_right ||
                error("Adjacent edge gap directions do not define one material interval")
            right_of_left
        end
        conductor = 0
        if metal
            conductor = i == 1 ? conductors[1] :
                        i == length(bounds) - 1 ? conductors[end] :
                        begin
                            conductors[i - 1] == conductors[i] ||
                                error("Both boundaries of a metal strip must use one conductor")
                            conductors[i]
                        end
        end
        push!(intervals, (
            left = bounds[i],
            right = bounds[i + 1],
            left_edge = i > 1,
            right_edge = i < length(bounds) - 1,
            metal = metal,
            conductor = conductor,
        ))
    end
    labels = unique(interval.conductor for interval in intervals if interval.metal)
    sort(labels) == collect(1:maximum(labels)) ||
        error("Conductor labels must be canonical and contiguous")
    return intervals
end

function interval_at(intervals, x, metal)
    tolerance = 1.0e-8
    candidates = [
        interval for interval in intervals
        if interval.metal == metal &&
           interval.left - tolerance <= x <= interval.right + tolerance
    ]
    isempty(candidates) && return nothing
    return candidates[argmin(
        abs(x - 0.5 * (interval.left + interval.right))
        for interval in candidates
    )]
end

function generate_edge_cluster_coupon(;
    spec::String,
    fabricated::Bool,
    radius::Float64 = 2.0,
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
    edges = read_edges(spec)
    intervals = classify_intervals(edges, radius)
    xmin = intervals[1].left
    xmax = intervals[end].right
    width = xmax - xmin
    tolerance = 1.0e-6 * max(radius, width)
    0.0 < sidewall_angle <= 90.0 ||
        error("sidewall_angle must lie in (0, 90]")

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("parallel_edge_cluster_$(fabricated ? "fabricated" : "thin")")
    occ = gmsh.model.occ

    if fabricated
        metal_pullback = t_metal / tan(deg2rad(sidewall_angle))
        trench_pullback = overetch / tan(deg2rad(sidewall_angle))
        outer = occ.addRectangle(xmin, -radius, 0.0, width, 2radius)
        metals = Tuple{Int32,Int32}[]
        for interval in intervals
            interval.metal || continue
            top_left = interval.left + (interval.left_edge ? metal_pullback : 0.0)
            top_right = interval.right - (interval.right_edge ? metal_pullback : 0.0)
            top_left < top_right ||
                error("Metal sidewall pullback closes a cluster strip")
            surface = shape(occ, [
                ((interval.left, 0.0), 0.0),
                ((interval.right, 0.0), 0.0),
                ((top_right, t_metal), interval.right_edge ? r_top : 0.0),
                ((top_left, t_metal), interval.left_edge ? r_top : 0.0),
            ])
            push!(metals, (2, surface))
        end
        substrate_base = occ.addRectangle(xmin, -radius, 0.0, width, radius)
        trenches = Tuple{Int32,Int32}[]
        if overetch > 0.0
            for interval in intervals
                interval.metal && continue
                bottom_left =
                    interval.left + (interval.left_edge ? trench_pullback : 0.0)
                bottom_right =
                    interval.right - (interval.right_edge ? trench_pullback : 0.0)
                bottom_left < bottom_right ||
                    error("Trench sidewall pullback closes a cluster gap")
                surface = shape(occ, [
                    ((interval.left, 0.0), 0.0),
                    ((bottom_left, -overetch), interval.left_edge ? r_bottom : 0.0),
                    ((bottom_right, -overetch), interval.right_edge ? r_bottom : 0.0),
                    ((interval.right, 0.0), 0.0),
                ])
                push!(trenches, (2, surface))
            end
        end
        substrate = isempty(trenches) ? [(2, substrate_base)] :
            first(occ.cut([(2, substrate_base)], trenches))
        field = first(occ.cut([(2, outer)], metals))
        vacuum = first(occ.cut(field, substrate, -1, true, false))
        occ.fragment(vcat(vacuum, substrate), [])
    else
        surfaces = Tuple{Int32,Int32}[]
        for interval in intervals
            push!(surfaces, (
                2,
                occ.addRectangle(
                    interval.left, -radius, 0.0,
                    interval.right - interval.left, radius),
            ))
            push!(surfaces, (
                2,
                occ.addRectangle(
                    interval.left, 0.0, 0.0,
                    interval.right - interval.left, radius),
            ))
        end
        occ.fragment(surfaces, [])
    end
    occ.synchronize()

    substrate_surfaces = Int32[]
    vacuum_surfaces = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        _, y, _ = occ.getCenterOfMass(dim, tag)
        push!(y < (fabricated ? -0.1 * max(overetch, t_metal) : 0.0) ?
              substrate_surfaces : vacuum_surfaces, tag)
    end
    substrate_set = Set(substrate_surfaces)
    vacuum_set = Set(vacuum_surfaces)

    outer_curves = Int32[]
    sa_curves = Int32[]
    ms_curves = Dict{Int,Vector{Int32}}()
    ma_curves = Dict{Int,Vector{Int32}}()
    for interval in intervals
        if interval.metal
            get!(ms_curves, interval.conductor, Int32[])
            get!(ma_curves, interval.conductor, Int32[])
        end
    end
    model_surfaces = Set(vcat(substrate_surfaces, vacuum_surfaces))
    for (dim, tag) in gmsh.model.getEntities(1)
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        adjacent = [surface for surface in up if surface in model_surfaces]
        isempty(adjacent) && continue
        xmin_curve, ymin, _, xmax_curve, ymax, _ =
            gmsh.model.getBoundingBox(dim, tag)
        xmid = 0.5 * (xmin_curve + xmax_curve)
        ymid = 0.5 * (ymin + ymax)
        horizontal = ymax - ymin < tolerance
        vertical = xmax_curve - xmin_curve < tolerance
        on_outer =
            (vertical && (abs(xmid - xmin) < tolerance ||
                          abs(xmid - xmax) < tolerance)) ||
            (horizontal && (abs(ymid + radius) < tolerance ||
                            abs(ymid - radius) < tolerance))
        if on_outer
            push!(outer_curves, tag)
            continue
        end

        adjacent_substrate = any(surface in substrate_set for surface in adjacent)
        adjacent_vacuum = any(surface in vacuum_set for surface in adjacent)
        if fabricated && adjacent_substrate && adjacent_vacuum
            push!(sa_curves, tag)
        elseif fabricated && ymin >= -tolerance &&
               ymax <= t_metal + tolerance
            interval = interval_at(intervals, xmid, true)
            if interval !== nothing
                target = horizontal && abs(ymid) < tolerance ?
                    ms_curves : ma_curves
                push!(target[interval.conductor], tag)
            end
        elseif !fabricated && horizontal && abs(ymid) < tolerance
            interval = interval_at(intervals, xmid, true)
            if interval === nothing
                push!(sa_curves, tag)
            else
                push!(ms_curves[interval.conductor], tag)
            end
        end
    end

    groups = [
        (2, substrate_surfaces, 1, "substrate"),
        (2, vacuum_surfaces, 2, "vacuum"),
        (1, outer_curves, 1, "matching_contour"),
        (1, sa_curves, 3, "SA"),
    ]
    for conductor in sort(collect(keys(ms_curves)))
        push!(groups, (
            1,
            ms_curves[conductor],
            10 + conductor,
            fabricated ? "MS_$conductor" : "thin_metal_$conductor",
        ))
        if fabricated
            push!(groups, (
                1,
                ma_curves[conductor],
                100 + conductor,
                "MA_$conductor",
            ))
        end
    end
    for (dim, entities, tag, name) in groups
        isempty(entities) && error("Empty physical group: $name")
        gmsh.model.addPhysicalGroup(dim, entities, tag, name)
    end

    features = vcat(
        sa_curves,
        collect(Iterators.flatten(values(ms_curves))),
        collect(Iterators.flatten(values(ma_curves))),
    )
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(features))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 10lc_fine)
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
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.setOrder(mesh_order)
    mesh_order > 1 && gmsh.model.mesh.optimize("HighOrderElastic")
    gmsh.write(filename)
    println("Parallel-edge cluster: fabricated=$fabricated, edges=$(length(edges))")
    println("  file=$filename")
    gmsh.finalize()
end

function parse_options(args)
    length(args) >= 3 ||
        error("Usage: mesh_edge_cluster_coupon.jl EDGES.csv thin|fabricated OUTPUT.msh [options]")
    options = Dict{String,Any}(
        "spec" => abspath(args[1]),
        "fabricated" => args[2] == "fabricated",
        "filename" => abspath(args[3]),
    )
    args[2] in ("thin", "fabricated") || error("Kind must be thin or fabricated")
    names = Dict(
        "--radius" => ("radius", Float64),
        "--metal-thickness" => ("t_metal", Float64),
        "--overetch" => ("overetch", Float64),
        "--sidewall-angle" => ("sidewall_angle", Float64),
        "--top-radius" => ("r_top", Float64),
        "--bottom-radius" => ("r_bottom", Float64),
        "--lc-fine" => ("lc_fine", Float64),
        "--lc-far" => ("lc_far", Float64),
        "--mesh-order" => ("mesh_order", Int),
    )
    index = 4
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
    options = parse_options(ARGS)
    generate_edge_cluster_coupon(;
        spec = options["spec"],
        fabricated = options["fabricated"],
        filename = options["filename"],
        radius = get(options, "radius", 2.0),
        t_metal = get(options, "t_metal", 0.1),
        overetch = get(options, "overetch", 0.05),
        sidewall_angle = get(options, "sidewall_angle", 80.0),
        r_top = get(options, "r_top", 0.01),
        r_bottom = get(options, "r_bottom", 0.01),
        lc_fine = get(options, "lc_fine", 0.002),
        lc_far = get(options, "lc_far", 0.05),
        mesh_order = get(options, "mesh_order", 2),
    )
end
