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
# {"Substrate": [per plane], "Vacuum": n}; `TerminalAttributeBase` b: terminal i (1-based in `Terminals`)
# gets b + 10 (i - 1) + k on plane k and b + 10 (i - 1) + (number of planes + 1) for its bump shells).
# Output: OUTPUT.msh2 (MSH 2.2 binary, order 1) and OUTPUT.json (attribute table, counts, sizes); the
# attribute table lists exactly the physical groups the mesh carries (a plane that is one metal sheet
# over the whole window has no `gap_<plane>`).
# Bump shells are grouped per conductor (decision 321): the ground bumps keep the chip's bump attribute
# as `bump_surface`; the bumps of terminal `<label>` form `bump_<label>` (a bump whose conductor is
# neither `ground` nor a terminal is refused), so a terminal's bump disc is never grounded.
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
#
# Ring orientation (decision 1022, census finding F-1): gmsh's `occ.addPlaneSurface` assumes
# every hole wire arrives in the SAME orientation as the outer wire and reverses it itself;
# a hole given in the opposite orientation yields an invalid face (area = outer PLUS holes)
# that `occ.fragment` can leave as a FREE sheet (meshed twice: orphan + duplicate triangles,
# an MFEM abort at mesh load) whenever the hole contains another metal polygon. GDS does not
# enforce ring orientation, so every loop is normalised to the mesher's contract before it
# reaches gmsh: outers (and bump footprints) counter-clockwise, holes co-oriented with their
# outer (a loop is re-traversed backwards from its first vertex; the vertices are
# unchanged). Three gates then fail closed (a named error, no mesh written). Before the
# fragment, the FACE GATE: the OCC area of every metal face must equal the polygon's net
# area (outer minus holes); a hole outside its outer, nested holes, a self-intersecting ring
# or a hole in the opposite orientation make gmsh build a malformed face silently (the
# census CHR-012 set carried four ground rings declared as holes of an island). After the
# fragment, the CONFORMITY GATE: every image of a metal sheet must be shared by exactly two
# volumes (the substrate and the vacuum), except an image inside a bump footprint column,
# which is one-sided by construction (the void column); a free or otherwise one-sided image
# (e.g. a metal polygon reaching beyond the box wall) is an error. And the COINCIDENT-SHEET
# GATE: an image claimed by more than one metal sheet (coincident sheets = duplicate
# boundary triangles) is an error. The manifest records `ReorientedLoops` (outers / holes),
# `ToolImagesNotShared` and `ToolImagesShared` (both 0 for a written mesh).

import Gmsh: gmsh
import JSON
using SHA: sha256

# Twice the signed area of a closed loop (positive = counter-clockwise).
function signed_area2(loop)
    a = 0.0
    n = length(loop)
    for i = 1:n
        p, q = loop[i], loop[i == n ? 1 : i + 1]
        a += Float64(p[1]) * Float64(q[2]) - Float64(q[1]) * Float64(p[2])
    end
    return a
end

# The loop traversed in the opposite direction, starting from the same first vertex.
reverse_loop(loop) = vcat([loop[1]], reverse(loop[2:end]))

# The polygon's loops in the mesher's contract: the outer counter-clockwise, every hole
# co-oriented with it. Returns the loops and the counts of reoriented outers and holes.
# `normalise = false` hands the loops to gmsh as given (the behaviour before decision 1022;
# the pre-fragment face gate then refuses the invalid face: used by the tests).
function oriented_loops(polygon; normalise=true)
    outer = polygon["Outer"]
    reoriented_outers = 0
    if normalise && signed_area2(outer) < 0
        outer = reverse_loop(outer)
        reoriented_outers = 1
    end
    loops = Any[outer]
    reoriented_holes = 0
    for hole in get(polygon, "Holes", [])
        if normalise && signed_area2(hole) < 0
            hole = reverse_loop(hole)
            reoriented_holes += 1
        end
        push!(loops, hole)
    end
    return loops, reoriented_outers, reoriented_holes
end

# Adds the polygon as a plane surface at height z; returns the surface tag and the
# reoriented (outers, holes) counts.
function add_polygon(occ, polygon, z; normalise=true)
    loops, reoriented_outers, reoriented_holes = oriented_loops(polygon; normalise)
    loop_tags = Int32[]
    for loop in loops
        tags = Int32[occ.addPoint(Float64(p[1]), Float64(p[2]), z) for p in loop]
        n = length(tags)
        curves = Int32[occ.addLine(tags[i], tags[i == n ? 1 : i + 1]) for i = 1:n]
        push!(loop_tags, occ.addCurveLoop(curves))
    end
    return occ.addPlaneSurface(loop_tags), reoriented_outers, reoriented_holes
end

# The polygon's net area (outer minus holes) from the shoelace formula.
function polygon_net_area(polygon)
    return abs(signed_area2(polygon["Outer"]) / 2) -
           sum(abs(signed_area2(hole) / 2) for hole in get(polygon, "Holes", []); init=0.0)
end

# The pre-fragment face gate: the OCC area of the plane surface built for the polygon must
# equal the polygon's net area, else gmsh built a malformed face (a hole outside its outer,
# nested holes, a self-intersecting ring). Errors with the conductor named; `occ` must be
# synchronized.
function check_face_area(face, polygon, label)
    net = polygon_net_area(polygon)
    area = gmsh.model.occ.getMass(2, face)
    abs(area - net) <= 1.0e-6 * max(net, 1.0) || error(
        "malformed metal face for $(label): OCC area $(area) vs the polygon's net area " *
        "$(net) (outer minus holes): a hole outside its outer, nested holes or a " *
        "self-intersecting ring; no mesh written"
    )
    return area
end

# The fragment images claimed by more than one metal sheet (coincident sheets: both groups
# would carry the same triangles). Returns [(face tag, [conductor labels])].
function shared_tool_images(domain_map, n_objects, tool_labels)
    owners = Dict{Int32, Vector{String}}()
    for (i, (k, label)) in enumerate(tool_labels)
        for (dim, tag) in domain_map[n_objects + i]
            dim == 2 || continue
            push!(get!(owners, tag, String[]), "$(label) (plane $(k))")
        end
    end
    shared = [(tag, labels) for (tag, labels) in owners if length(labels) > 1]
    return sort!(shared; by=first)
end

# The conformity gate on the fragment images of the metal sheets (decision 1022): every
# image must be shared by exactly two volumes, except a one-sided image inside a bump
# footprint column (the bump columns are void: `footprint_columns` = [(x0, x1, y0, y1, z0,
# z1)]). Returns the images that fail.
function unshared_tool_images(domain_map, n_objects, tool_labels, footprint_columns, tol)
    unshared = Tuple{String, Int32, Int}[]  # (conductor label, face tag, adjacent volumes)
    for (i, (k, label)) in enumerate(tool_labels)
        for (dim, tag) in domain_map[n_objects + i]
            dim == 2 || continue
            volumes, _ = gmsh.model.getAdjacencies(2, tag)
            length(volumes) == 2 && continue
            if length(volumes) == 1
                bounds = gmsh.model.getBoundingBox(2, tag)
                any(
                    column -> bbox_within(bounds, column..., tol),
                    footprint_columns
                ) && continue
            end
            push!(unshared, ("$(label) (plane $(k))", tag, length(volumes)))
        end
    end
    return unshared
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

function mesh_thin_window(
    window_path,
    chip_path,
    output,
    lc_fine,
    lc_far;
    normalise_loops=true
)
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
    reoriented_outers = 0
    reoriented_holes = 0
    footprint_columns = NTuple{6, Float64}[]
    if !isempty(bumps)
        length(planes) == 2 || error("bumps need two planes")
        z_low, z_high = minimum(z_levels), maximum(z_levels)
        prisms = Tuple{Int32, Int32}[]
        for bump in bumps
            face, outers, holes = add_polygon(
                occ,
                Dict("Outer" => bump["Footprint"]),
                z_low;
                normalise=normalise_loops
            )
            reoriented_outers += outers
            reoriented_holes += holes
            extruded = occ.extrude([(2, face)], 0.0, 0.0, z_high - z_low)
            for (dim, tag) in extruded
                dim == 3 && push!(prisms, (3, tag))
            end
            push!(
                footprint_columns,
                (
                    minimum(Float64(p[1]) for p in bump["Footprint"]),
                    maximum(Float64(p[1]) for p in bump["Footprint"]),
                    minimum(Float64(p[2]) for p in bump["Footprint"]),
                    maximum(Float64(p[2]) for p in bump["Footprint"]),
                    z_low,
                    z_high
                )
            )
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
            face, outers, holes = add_polygon(occ, polygon, z; normalise=normalise_loops)
            reoriented_outers += outers
            reoriented_holes += holes
            push!(tools, (2, face))
            push!(tool_labels, (k, String(polygon["Conductor"])))
        end
    end
    # the pre-fragment face gate: every metal face has the polygon's net area
    occ.synchronize()
    polygons = [polygon for plane in planes for polygon in plane["Polygons"]]
    for (i, (dim, face)) in enumerate(tools)
        k, label = tool_labels[i]
        check_face_area(face, polygons[i], "$(label) (plane $(k))")
    end
    objects = vcat(
        [(Int32(3), t) for t in substrate_boxes],
        [(Int32(3), t) for t in vacuum_boxes]
    )
    _, domain_map = occ.fragment(objects, tools)
    occ.synchronize()
    n_objects = length(objects)
    tol = 1.0e-6 * max(x1 - x0, y1 - y0, zmax_all - zmin_all)
    # the conformity gate (decision 1022): a metal sheet image not shared by the substrate
    # and the vacuum would be meshed as a free sheet (orphan + duplicate triangles); fail
    # closed, no mesh
    unshared =
        unshared_tool_images(domain_map, n_objects, tool_labels, footprint_columns, tol)
    if !isempty(unshared)
        first_unshared = [
            "$(label) face $(tag) with $(n) adjacent volume(s)" for
            (label, tag, n) in unshared[1:min(end, 5)]
        ]
        error(
            "conformity gate: $(length(unshared)) metal sheet image(s) not shared by two " *
            "volumes after the fragment (a free sheet would be meshed twice); no mesh " *
            "written. First: " * join(first_unshared, "; ")
        )
    end
    shared = shared_tool_images(domain_map, n_objects, tool_labels)
    if !isempty(shared)
        first_shared = [
            "face $(tag) claimed by $(join(labels, ", "))" for
            (tag, labels) in shared[1:min(end, 5)]
        ]
        error(
            "conformity gate: $(length(shared)) fragment image(s) claimed by more than " *
            "one metal sheet (coincident sheets would carry duplicate triangles); no mesh " *
            "written. First: " * join(first_shared, "; ")
        )
    end
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
    bump_conductor = Dict{Int32, String}()  # bump face -> the conductor of its bump
    if !isempty(bumps)
        bump_labels = [String(get(b, "Conductor", "ground")) for b in bumps]
        for label in bump_labels
            label == "ground" ||
                label in terminals ||
                error(
                    "bump conductor $(label) is neither ground nor one of Terminals $(terminals)"
                )
        end
        for (dim, tag) in gmsh.model.getEntities(2)
            bounds = gmsh.model.getBoundingBox(dim, tag)
            for (i, column) in enumerate(footprint_columns)
                if bbox_within(bounds, column..., tol)
                    push!(bump_faces, tag)
                    bump_conductor[tag] = bump_labels[i]
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
        if !isempty(gap[k])
            gmsh.model.addPhysicalGroup(
                2,
                gap[k],
                Int(chip_planes[k]["Gap"][1]),
                "gap_$(name)"
            )
            attributes["gap_$(name)"] = Int(chip_planes[k]["Gap"][1])
        end
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
        # Ground bump shells keep the chip's bump attribute (the legacy single group); a terminal's
        # bump shells get their own group so the config writer puts them into that Terminal.
        ground_bump_faces = [t for t in bump_faces if bump_conductor[t] == "ground"]
        if !isempty(ground_bump_faces)
            gmsh.model.addPhysicalGroup(
                2,
                ground_bump_faces,
                Int(chip["Bump"][1]),
                "bump_surface"
            )
            attributes["bump_surface"] = Int(chip["Bump"][1])
        end
        for (i, label) in enumerate(terminals)
            faces = [t for t in bump_faces if bump_conductor[t] == label]
            isempty(faces) && continue
            attribute = base + 10 * (i - 1) + length(planes) + 1
            gmsh.model.addPhysicalGroup(2, faces, attribute, "bump_$(label)")
            attributes["bump_$(label)"] = attribute
        end
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
        "ReorientedLoops" =>
            Dict("Outers" => reoriented_outers, "Holes" => reoriented_holes),
        "ToolImagesNotShared" => length(unshared),
        "ToolImagesShared" => length(shared),
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
