# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Synthetic thin-metal layouts for the geometry-identification stress tests (Gmsh / OCC).
#
# Reads the plain-text layout specification written by synthetic_layouts.py and writes one
# binary MSH 2.2 mesh per layout. Every layout is a substrate box (z < 0) under a vacuum box
# (z > 0) with zero-thickness metal sheets on the z = 0 process plane; optionally a facing
# metal sheet on a second plane z = z_layer (the cross-layer class of decision 73(3)) and
# vertical metal walls (the non-planar class). Physical groups (same convention as the
# transmon example): 1 substrate, 2 vacuum, 3 outer (every box face), 8 substrate_air (the
# non-metal part of the process plane), and the metal attributes named by the specification.
#
#     julia synthetic_layouts.jl SPEC.txt OUTPUT_DIRECTORY [NAME ...]
#
# Specification grammar (one statement per line, '#' comments):
#     layout NAME
#     box HALF_X HALF_Y DEPTH HEIGHT
#     size LC_FINE LC_FAR
#     polygon ATTRIBUTE            # a metal sheet on z = 0; the loop follows
#     hole                          # the next loop is a hole of the current polygon
#     v X Y                         # loop vertex
#     arc X Y CX CY                 # circular arc from the previous vertex to (X, Y), centre (CX, CY)
#     close                         # end of the loop
#     sheet ATTRIBUTE Z             # a metal polygon (following loop) on the plane z = Z
#     wall ATTRIBUTE X0 Y0 X1 Y1 H  # a vertical metal rectangle from (X0, Y0, 0) to (X1, Y1, H)
#     port ATTRIBUTE                # a non-metal polygon on z = 0 with its own physical group
#                                   # (a lumped port face; the loop follows)
#     end

import Gmsh: gmsh

mutable struct Loop
    points::Vector{Tuple{Float64, Float64}}
    arcs::Dict{Int, Tuple{Float64, Float64}}  # index of the end vertex -> centre
end
Loop() = Loop(Tuple{Float64, Float64}[], Dict{Int, Tuple{Float64, Float64}}())

mutable struct Sheet
    attribute::Int
    z::Float64
    loops::Vector{Loop}  # first = outer boundary, the rest holes
end

struct Wall
    attribute::Int
    x0::Float64
    y0::Float64
    x1::Float64
    y1::Float64
    height::Float64
end

mutable struct Layout
    name::String
    half_x::Float64
    half_y::Float64
    depth::Float64
    height::Float64
    lc_fine::Float64
    lc_far::Float64
    sheets::Vector{Sheet}
    walls::Vector{Wall}
    ports::Vector{Sheet}
end
Layout(name) = Layout(name, 30.0, 30.0, 15.0, 15.0, 1.0, 6.0, Sheet[], Wall[], Sheet[])

function parse_specification(path)
    layouts = Layout[]
    current = nothing
    sheet = nothing
    loop = nothing
    pending_hole = false
    for raw in eachline(path)
        line = strip(first(split(raw, '#'; limit=2)))
        isempty(line) && continue
        fields = split(line)
        keyword = fields[1]
        numbers = [parse(Float64, f) for f in fields[2:end] if tryparse(Float64, f) !== nothing]
        if keyword == "layout"
            current = Layout(String(fields[2]))
            push!(layouts, current)
        elseif keyword == "box"
            current.half_x, current.half_y, current.depth, current.height = numbers
        elseif keyword == "size"
            current.lc_fine, current.lc_far = numbers
        elseif keyword == "polygon" || keyword == "sheet" || keyword == "port"
            z = keyword == "sheet" ? numbers[2] : 0.0
            sheet = Sheet(Int(numbers[1]), z, Loop[])
            push!(keyword == "port" ? current.ports : current.sheets, sheet)
            loop = Loop()
            push!(sheet.loops, loop)
            pending_hole = false
        elseif keyword == "hole"
            loop = Loop()
            push!(sheet.loops, loop)
        elseif keyword == "v"
            push!(loop.points, (numbers[1], numbers[2]))
        elseif keyword == "arc"
            push!(loop.points, (numbers[1], numbers[2]))
            loop.arcs[length(loop.points)] = (numbers[3], numbers[4])
        elseif keyword == "close"
            loop = nothing
        elseif keyword == "wall"
            push!(current.walls, Wall(Int(numbers[1]), numbers[2:6]...))
        elseif keyword == "end"
            current = nothing
        else
            error("unknown statement: $line")
        end
    end
    return layouts
end

function add_loop(occ, loop, z)
    tags = Int32[]
    n = length(loop.points)
    for (x, y) in loop.points
        push!(tags, occ.addPoint(x, y, z))
    end
    curves = Int32[]
    for i in 1:n
        a = tags[i]
        b = tags[i == n ? 1 : i + 1]
        j = i == n ? 1 : i + 1
        if haskey(loop.arcs, j)
            cx, cy = loop.arcs[j]
            centre = occ.addPoint(cx, cy, z)
            push!(curves, occ.addCircleArc(a, centre, b))
        else
            push!(curves, occ.addLine(a, b))
        end
    end
    return occ.addCurveLoop(curves)
end

function add_sheet(occ, sheet)
    loops = [add_loop(occ, loop, sheet.z) for loop in sheet.loops]
    return occ.addPlaneSurface(loops)
end

function add_wall(occ, wall)
    p = [
        occ.addPoint(wall.x0, wall.y0, 0.0),
        occ.addPoint(wall.x1, wall.y1, 0.0),
        occ.addPoint(wall.x1, wall.y1, wall.height),
        occ.addPoint(wall.x0, wall.y0, wall.height)
    ]
    curves = [occ.addLine(p[i], p[i == 4 ? 1 : i + 1]) for i in 1:4]
    return occ.addPlaneSurface([occ.addCurveLoop(curves)])
end

function on_outer_box(bounds, layout, tolerance)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    return (abs(xmin + layout.half_x) < tolerance && abs(xmax + layout.half_x) < tolerance) ||
           (abs(xmin - layout.half_x) < tolerance && abs(xmax - layout.half_x) < tolerance) ||
           (abs(ymin + layout.half_y) < tolerance && abs(ymax + layout.half_y) < tolerance) ||
           (abs(ymin - layout.half_y) < tolerance && abs(ymax - layout.half_y) < tolerance) ||
           (abs(zmin + layout.depth) < tolerance && abs(zmax + layout.depth) < tolerance) ||
           (abs(zmin - layout.height) < tolerance && abs(zmax - layout.height) < tolerance)
end

function generate(layout, output_directory)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(layout.name)
    occ = gmsh.model.occ
    substrate = occ.addBox(-layout.half_x, -layout.half_y, -layout.depth, 2layout.half_x, 2layout.half_y, layout.depth)
    vacuum = occ.addBox(-layout.half_x, -layout.half_y, 0.0, 2layout.half_x, 2layout.half_y, layout.height)
    tools = Tuple{Int32, Int32}[]
    for sheet in layout.sheets
        push!(tools, (2, add_sheet(occ, sheet)))
    end
    for wall in layout.walls
        push!(tools, (2, add_wall(occ, wall)))
    end
    for port in layout.ports
        push!(tools, (2, add_sheet(occ, port)))
    end
    domains, domain_map = occ.fragment([(3, substrate), (3, vacuum)], tools)
    occ.synchronize()
    substrate_tags = [tag for (dim, tag) in domain_map[1] if dim == 3]
    vacuum_tags = [tag for (dim, tag) in domain_map[2] if dim == 3]

    tolerance = 1.0e-6 * max(layout.half_x, layout.half_y)
    # The metal faces are the fragment images of the tool sheets (domain_map entries after
    # the two boxes); every other face on the process plane is substrate_air.
    metal = Dict{Int, Vector{Int32}}()
    tool_attributes = vcat([sheet.attribute for sheet in layout.sheets], [wall.attribute for wall in layout.walls])
    metal_faces = Set{Int32}()
    for (k, attribute) in enumerate(tool_attributes)
        for (dim, tag) in domain_map[2 + k]
            dim == 2 || continue
            push!(get!(metal, attribute, Int32[]), tag)
            push!(metal_faces, tag)
        end
    end
    # Port faces: non-metal, their own physical group, not substrate_air.
    port_groups = Dict{Int, Vector{Int32}}()
    for (k, port) in enumerate(layout.ports)
        for (dim, tag) in domain_map[2 + length(tool_attributes) + k]
            dim == 2 || continue
            push!(get!(port_groups, port.attribute, Int32[]), tag)
            push!(metal_faces, tag)  # excluded from substrate_air below
        end
    end
    outer = Int32[]
    substrate_air = Int32[]
    for (dim, tag) in gmsh.model.getEntities(2)
        tag in metal_faces && continue
        bounds = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bounds
        if on_outer_box(bounds, layout, tolerance)
            push!(outer, tag)
        elseif abs(zmin) < tolerance && abs(zmax) < tolerance
            push!(substrate_air, tag)
        end
    end

    gmsh.model.addPhysicalGroup(3, substrate_tags, 1, "substrate")
    gmsh.model.addPhysicalGroup(3, vacuum_tags, 2, "vacuum")
    gmsh.model.addPhysicalGroup(2, outer, 3, "outer")
    isempty(substrate_air) || gmsh.model.addPhysicalGroup(2, substrate_air, 8, "substrate_air")
    for (attribute, tags) in sort!(collect(port_groups); by=first)
        gmsh.model.addPhysicalGroup(2, tags, attribute, "port_$(attribute)")
    end
    metal_curves = Int32[]
    for (attribute, tags) in sort!(collect(metal); by=first)
        gmsh.model.addPhysicalGroup(2, tags, attribute, "metal_$(attribute)")
        for (dim, curve) in gmsh.model.getBoundary([(2, tag) for tag in tags], false, false, false)
            dim == 1 || continue
            bounds = gmsh.model.getBoundingBox(dim, curve)
            on_outer_box(bounds, layout, tolerance) || push!(metal_curves, abs(curve))
        end
    end
    unique!(metal_curves)

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(metal_curves))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 200)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", layout.lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", layout.lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2layout.lc_fine)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 6layout.lc_fine)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeMin", layout.lc_fine)
    gmsh.option.setNumber("Mesh.MeshSizeMax", layout.lc_far)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.SaveAll", 0)
    gmsh.model.mesh.generate(3)
    path = joinpath(output_directory, layout.name * ".msh2")
    gmsh.write(path)
    nodes = length(gmsh.model.mesh.getNodes()[1])
    gmsh.finalize()
    println("$(layout.name) $(path) nodes=$(nodes) metal=$(join(sort(collect(keys(metal))), ",")) substrate_air=$(length(substrate_air))")
    return path
end

function main(args)
    length(args) >= 2 || error("usage: julia synthetic_layouts.jl SPEC OUTPUT_DIRECTORY [NAME ...]")
    layouts = parse_specification(args[1])
    mkpath(args[2])
    selected = Set(args[3:end])
    for layout in layouts
        (isempty(selected) || layout.name in selected) || continue
        generate(layout, args[2])
    end
end

main(ARGS)
