# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Thin window mesh of a validation window from its polygon set (lane W of USER decision 184, task 4): the
# surface-response device recipe of the chip meshes — zero-thickness metal sheets on the plane(s), the
# substrate(s) and the vacuum as boxes, bumps as void prisms between the planes with a metal shell, the
# chip's attribute numbers (so the chip's preflight configurations apply to the window) — built from
# EXACTLY the polygon set lane F's fabricated mesher reads (`window_polygons.py`, schema Version 1), the
# chip's chords kept (no snapping; decision 190 (c)).
#
#     julia --project=<env with Gmsh, JSON> mesh_thin_window.jl WINDOW.json CHIP.json OUTPUT.msh2 [LC_FINE LC_FAR]
#
# WINDOW.json: the polygon set (Box, Planes with Polygons / Conductor labels, Bumps, Vacuum, Terminals).
# CHIP.json: the chip description of window_polygons.py (per plane `Attributes[1]` = the ground metal
# sheet attribute, `Gap[1]` = the substrate-air attribute; `Bump[1]`, `Exterior[1]`; `Volumes` =
# {"Substrate": [per plane], "Vacuum": n}; `TerminalAttributeBase` b: the terminal of plane k gets b + k).
# Output: OUTPUT.msh2 (MSH 2.2 binary, order 1) and OUTPUT.json (attribute table, counts, sizes).
#
# Geometry: plane k facing up has its substrate box in [SurfaceZ - T, SurfaceZ], facing down in
# [SurfaceZ, SurfaceZ + T]; the vacuum fills the box between the substrates (two planes) or above the
# single plane up to `Vacuum.Above`; `Vacuum.Below` / `Above` > 0 add vacuum boxes beyond the backsides.
# Bumps: the footprint extruded from the lower plane to the upper plane is CUT from the vacuum (a void
# column; its sidewalls and both footprint faces carry the bump attribute, as the chip mesh draws them).
# Metal polygons (outer loop + holes) are plane surfaces fragmented into the volumes; their fragment
# images are the metal faces (label -> attribute), every other face on a plane is the plane's gap, every
# face on the outer box is the exterior boundary. Mesh size: LC_FINE on the metal edges (Distance /
# Threshold field) growing to LC_FAR (defaults 4 / 40 um).

import Gmsh: gmsh
import JSON
using SHA: sha256

function add_polygon(occ, polygon, z)
    loops = Int32[]
    for loop in vcat([polygon["Outer"]], get(polygon, "Holes", []))
        tags = Int32[occ.addPoint(Float64(p[1]), Float64(p[2]), z) for p in loop]
        n = length(tags)
        curves = Int32[occ.addLine(tags[i], tags[i == n ? 1 : i + 1]) for i = 1:n]
        push!(loops, occ.addCurveLoop(curves))
    end
    return occ.addPlaneSurface(loops)
end

function bbox_within(bounds, x0, x1, y0, y1, z0, z1, tol)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    return xmin >= x0 - tol &&
           xmax <= x1 + tol &&
           ymin >= y0 - tol &&
           ymax <= y1 + tol &&
           zmin >= z0 - tol &&
           zmax <= z1 + tol
end

function on_box_wall(bounds, x0, x1, y0, y1, z0, z1, tol)
    xmin, ymin, zmin, xmax, ymax, zmax = bounds
    return (abs(xmin - x0) < tol && abs(xmax - x0) < tol) ||
           (abs(xmin - x1) < tol && abs(xmax - x1) < tol) ||
           (abs(ymin - y0) < tol && abs(ymax - y0) < tol) ||
           (abs(ymin - y1) < tol && abs(ymax - y1) < tol) ||
           (abs(zmin - z0) < tol && abs(zmax - z0) < tol) ||
           (abs(zmin - z1) < tol && abs(zmax - z1) < tol)
end

function mesh_thin_window(window_path, chip_path, output, lc_fine, lc_far)
    spec = JSON.parsefile(window_path)
    chip = JSON.parsefile(chip_path)
    x0, x1 = Float64.(spec["Box"]["X"])
    y0, y1 = Float64.(spec["Box"]["Y"])
    planes = spec["Planes"]
    chip_planes = chip["Planes"]
    length(planes) == length(chip_planes) || error(
        "the polygon set has $(length(planes)) planes, the chip description $(length(chip_planes))"
    )
    vacuum = get(spec, "Vacuum", Dict("Below" => 0.0, "Above" => 0.0))
    below, above = Float64(get(vacuum, "Below", 0.0)), Float64(get(vacuum, "Above", 0.0))
    volumes = chip["Volumes"]
    base = Int(get(chip, "TerminalAttributeBase", 1000))
    terminals = String.(get(spec, "Terminals", String[]))

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add(String(get(spec, "Name", "window")))
    occ = gmsh.model.occ
    # substrate boxes; the vacuum between the planes (or above the single plane); extra vacuum beyond the backsides
    z_levels = Float64[]
    substrate_boxes = Int32[]
    substrate_span = Tuple{Float64, Float64}[]
    for plane in planes
        z = Float64(plane["SurfaceZ"])
        t = Float64(plane["SubstrateThickness"])
        push!(z_levels, z)
        lo, hi = plane["Facing"] == "up" ? (z - t, z) : (z, z + t)
        push!(substrate_span, (lo, hi))
        push!(substrate_boxes, occ.addBox(x0, y0, lo, x1 - x0, y1 - y0, hi - lo))
    end
    zmin_all = minimum(first.(substrate_span))
    zmax_all = maximum(last.(substrate_span))
    vacuum_boxes = Int32[]
    if length(planes) == 2
        z_low, z_high = minimum(z_levels), maximum(z_levels)
        push!(vacuum_boxes, occ.addBox(x0, y0, z_low, x1 - x0, y1 - y0, z_high - z_low))
    else
        above > 0.0 || error("a single plane needs Vacuum.Above > 0")
        z = z_levels[1]
        if planes[1]["Facing"] == "up"
            push!(vacuum_boxes, occ.addBox(x0, y0, z, x1 - x0, y1 - y0, above))
            zmax_all = z + above
            above = 0.0
        else
            push!(vacuum_boxes, occ.addBox(x0, y0, z - above, x1 - x0, y1 - y0, above))
            zmin_all = z - above
            above = 0.0
        end
    end
    if below > 0.0
        push!(vacuum_boxes, occ.addBox(x0, y0, zmin_all - below, x1 - x0, y1 - y0, below))
        zmin_all -= below
    end
    if above > 0.0
        push!(vacuum_boxes, occ.addBox(x0, y0, zmax_all, x1 - x0, y1 - y0, above))
        zmax_all += above
    end
    # bump columns cut from the vacuum between the planes
    bumps = get(spec, "Bumps", [])
    if !isempty(bumps)
        length(planes) == 2 || error("bumps need two planes")
        z_low, z_high = minimum(z_levels), maximum(z_levels)
        prisms = Tuple{Int32, Int32}[]
        for bump in bumps
            face = add_polygon(occ, Dict("Outer" => bump["Footprint"]), z_low)
            extruded = occ.extrude([(2, face)], 0.0, 0.0, z_high - z_low)
            for (dim, tag) in extruded
                dim == 3 && push!(prisms, (3, tag))
            end
        end
        cut, _ = occ.cut([(3, vacuum_boxes[1])], prisms)
        vacuum_boxes[1] = cut[1][2]
    end
    # metal polygons as tool sheets
    tools = Tuple{Int32, Int32}[]
    tool_labels = Tuple{Int, String}[]  # (plane index, conductor label)
    for (k, plane) in enumerate(planes)
        z = Float64(plane["SurfaceZ"])
        for polygon in plane["Polygons"]
            push!(tools, (2, add_polygon(occ, polygon, z)))
            push!(tool_labels, (k, String(polygon["Conductor"])))
        end
    end
    objects = vcat(
        [(Int32(3), t) for t in substrate_boxes],
        [(Int32(3), t) for t in vacuum_boxes]
    )
    _, domain_map = occ.fragment(objects, tools)
    occ.synchronize()
    n_objects = length(objects)
    substrate_tags = [Int32[] for _ in planes]
    for (k, _) in enumerate(planes)
        substrate_tags[k] = [tag for (dim, tag) in domain_map[k] if dim == 3]
    end
    vacuum_tags = Int32[]
    for j = 1:length(vacuum_boxes)
        append!(
            vacuum_tags,
            [tag for (dim, tag) in domain_map[length(planes) + j] if dim == 3]
        )
    end
    tol = 1.0e-6 * max(x1 - x0, y1 - y0, zmax_all - zmin_all)
    # metal faces per (plane, label); bump faces (planar footprint faces re-assigned, sidewalls by bounding box)
    metal = Dict{Tuple{Int, String}, Vector{Int32}}()
    assigned = Set{Int32}()
    for (i, (k, label)) in enumerate(tool_labels)
        for (dim, tag) in domain_map[n_objects + i]
            dim == 2 || continue
            push!(get!(metal, (k, label), Int32[]), tag)
            push!(assigned, tag)
        end
    end
    bump_faces = Int32[]
    if !isempty(bumps)
        z_low, z_high = minimum(z_levels), maximum(z_levels)
        footprints = [
            (
                minimum(Float64(p[1]) for p in b["Footprint"]),
                maximum(Float64(p[1]) for p in b["Footprint"]),
                minimum(Float64(p[2]) for p in b["Footprint"]),
                maximum(Float64(p[2]) for p in b["Footprint"])
            ) for b in bumps
        ]
        for (dim, tag) in gmsh.model.getEntities(2)
            bounds = gmsh.model.getBoundingBox(dim, tag)
            for (bx0, bx1, by0, by1) in footprints
                if bbox_within(bounds, bx0, bx1, by0, by1, z_low, z_high, tol)
                    push!(bump_faces, tag)
                    break
                end
            end
        end
        bump_set = Set(bump_faces)
        for (key, tags) in metal
            filter!(t -> !(t in bump_set), tags)
        end
        union!(assigned, bump_set)
    end
    exterior = Int32[]
    gap = [Int32[] for _ in planes]
    for (dim, tag) in gmsh.model.getEntities(2)
        tag in assigned && continue
        bounds = gmsh.model.getBoundingBox(dim, tag)
        if on_box_wall(bounds, x0, x1, y0, y1, zmin_all, zmax_all, tol)
            push!(exterior, tag)
            continue
        end
        for (k, z) in enumerate(z_levels)
            if abs(bounds[3] - z) < tol && abs(bounds[6] - z) < tol
                push!(gap[k], tag)
                break
            end
        end
    end
    # physical groups with the chip's attribute numbers
    attributes = Dict{String, Any}()
    for (k, plane) in enumerate(planes)
        a = Int(volumes["Substrate"][k])
        gmsh.model.addPhysicalGroup(
            3,
            substrate_tags[k],
            a,
            "substrate_$(lowercase(plane["Name"]))"
        )
        attributes["substrate_$(lowercase(plane["Name"]))"] = a
    end
    gmsh.model.addPhysicalGroup(3, vacuum_tags, Int(volumes["Vacuum"]), "vacuum")
    attributes["vacuum"] = Int(volumes["Vacuum"])
    gmsh.model.addPhysicalGroup(2, exterior, Int(chip["Exterior"][1]), "exterior_boundary")
    attributes["exterior_boundary"] = Int(chip["Exterior"][1])
    metal_curves = Int32[]
    for (k, plane) in enumerate(planes)
        name = plane["Name"]
        isempty(gap[k]) || gmsh.model.addPhysicalGroup(
            2,
            gap[k],
            Int(chip_planes[k]["Gap"][1]),
            "gap_$(name)"
        )
        attributes["gap_$(name)"] = Int(chip_planes[k]["Gap"][1])
        for ((kk, label), tags) in sort!(collect(metal); by=first)
            kk == k || continue
            isempty(tags) && continue
            attribute =
                label == "ground" ? Int(chip_planes[k]["Attributes"][1]) :
                base + k + 10 * (findfirst(==(label), terminals) - 1)
            gmsh.model.addPhysicalGroup(2, tags, attribute, "$(label)_$(name)")
            attributes["$(label)_$(name)"] = attribute
            for (dim, curve) in
                gmsh.model.getBoundary([(2, t) for t in tags], false, false, false)
                dim == 1 || continue
                bounds = gmsh.model.getBoundingBox(dim, curve)
                on_box_wall(bounds, x0, x1, y0, y1, zmin_all, zmax_all, tol) ||
                    push!(metal_curves, abs(curve))
            end
        end
    end
    if !isempty(bump_faces)
        gmsh.model.addPhysicalGroup(2, bump_faces, Int(chip["Bump"][1]), "bump_surface")
        attributes["bump_surface"] = Int(chip["Bump"][1])
        for (dim, curve) in
            gmsh.model.getBoundary([(2, t) for t in bump_faces], false, false, false)
            dim == 1 && push!(metal_curves, abs(curve))
        end
    end
    unique!(metal_curves)
    # mesh size: fine on the metal edges, coarse elsewhere
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(metal_curves))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 2lc_fine)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 10lc_fine)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc_fine)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
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
    gmsh.write(output)
    # manifest
    counts = Dict{String, Any}()
    for (dim, tag) in gmsh.model.getPhysicalGroups()
        name = gmsh.model.getPhysicalName(dim, tag)
        n = 0
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity)
            n += sum(length.(element_tags); init=0)
        end
        counts[name] = Dict("Dimension" => dim, "Attribute" => tag, "Elements" => n)
    end
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    manifest = Dict(
        "Window" => window_path,
        "Chip" => chip_path,
        "Output" => output,
        "SHA256" => bytes2hex(open(sha256, output)),
        "Bytes" => filesize(output),
        "Nodes" => length(node_tags),
        "Box" => Dict("X" => [x0, x1], "Y" => [y0, y1], "Z" => [zmin_all, zmax_all]),
        "Planes" => [
            Dict("Name" => p["Name"], "SurfaceZ" => p["SurfaceZ"], "Facing" => p["Facing"]) for p in planes
        ],
        "Bumps" => length(bumps),
        "Terminals" => terminals,
        "Attributes" => attributes,
        "Groups" => counts,
        "MeshSize" => Dict("Fine" => lc_fine, "Far" => lc_far),
        "GmshVersion" => gmsh.option.getString("General.Version"),
        "Recipe" => "thin sheets at the planes, void bump columns with a metal shell, order 1, MSH 2.2 binary"
    )
    open(replace(output, r"\.msh2$" => ".json"), "w") do stream
        JSON.print(stream, manifest, 1)
        return println(stream)
    end
    gmsh.finalize()
    return manifest
end

function main(args)
    length(args) >= 3 || error(
        "usage: mesh_thin_window.jl WINDOW.json CHIP.json OUTPUT.msh2 [LC_FINE LC_FAR]"
    )
    lc_fine = length(args) >= 4 ? parse(Float64, args[4]) : 4.0
    lc_far = length(args) >= 5 ? parse(Float64, args[5]) : 40.0
    manifest = mesh_thin_window(args[1], args[2], args[3], lc_fine, lc_far)
    println(
        "Saved ",
        args[3],
        ": nodes ",
        manifest["Nodes"],
        ", groups ",
        join(
            sort!([
                "$(k)=$(v["Attribute"]):$(v["Elements"])" for (k, v) in manifest["Groups"]
            ]),
            " "
        )
    )
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end
