# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Step 0 of the prism-tube SPIKE (decision 37): a small box with one straight metal
# ridge on an etched substrate, meshed either as a hybrid (prism edge tubes around
# the top and bottom metal edges, pyramids at the tube/tet interface, tets
# elsewhere) or as pure tets with the same size field. Gmsh 2.2 binary output for
# the local Palace electrostatic check.
#
# usage: julia --project=test/examples mesh_tiny_hybrid.jl hybrid|tets output.msh [--inner-size um]
#        [--spacing um] [--lc-edge um] [--pyramid-height um] [--census output.json]

include(joinpath(@__DIR__, "prism_tube.jl"))
include(joinpath(@__DIR__, "hybrid_census.jl"))

function tiny_hybrid_mesh(kind, filename; inner_size=0.00025, ratio=2.0, rings=7,
                          spacing=0.05, sector_degrees=30.0, lc_edge=0.025, lc_far=0.15,
                          pyramid_height=0.008, census_path=nothing)
    kind in ("hybrid", "tets") || error("kind must be hybrid or tets")
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 2)
    gmsh.model.add("tiny_$kind")
    occ = gmsh.model.occ
    box = (3, occ.addBox(-1.0, -1.0, -0.55, 2.0, 2.0, 1.15))
    slab = [(3, occ.addBox(-1.0, -1.0, -0.55, 2.0, 2.0, 0.55))]
    trench = [(3, occ.addBox(0.0, -1.0, -0.05, 0.5, 2.0, 0.05))]
    substrate, _ = occ.cut(slab, trench)
    metal = [(3, occ.addBox(-1.0, -1.0, 0.0, 1.0, 2.0, 0.1))]
    field, _ = occ.cut([box], metal)
    vacuum, _ = occ.cut(field, substrate, -1, true, false)
    objects = vcat(substrate, vacuum)
    sectors = round(Int, 270.0 / sector_degrees)
    top_section = TubeSection(inner_size, ratio, rings,
                              [-90.0 + sector_degrees * j for j in 0:sectors], fill(2, sectors))
    per_quadrant = round(Int, 90.0 / sector_degrees)
    bottom_section = TubeSection(inner_size, ratio, rings,
                                 [180.0 + sector_degrees * j for j in 0:sectors],
                                 vcat(fill(1, per_quadrant), fill(2, sectors - per_quadrant)))
    tubes = Tuple{EdgeTube, TubeSection}[]
    records = TubeRecord[]
    if kind == "hybrid"
        push!(tubes, (EdgeTube([0.0, 1.0, 0.1], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 2.0, spacing),
                      top_section))
        push!(tubes, (EdgeTube([0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.0, 2.0, spacing),
                      bottom_section))
        for (tube, section) in tubes
            for (tool, group) in add_tube_volumes!(occ, tube, section)
                push!(records, TubeRecord(tube, section, group, tool))
            end
        end
    end
    tools = [record.tool for record in records]
    domains, fragment_map = occ.fragment(objects, tools)
    occ.synchronize()
    substrate_tags = Int32[]
    vacuum_tags = Int32[]
    for (i, object) in enumerate(objects)
        descendants = [tag for (dim, tag) in fragment_map[i] if dim == 3]
        append!(object in substrate ? substrate_tags : vacuum_tags, descendants)
    end
    unique!(substrate_tags)
    unique!(vacuum_tags)
    tube_volumes = Int32[]
    for (index, record) in enumerate(records)
        push!(tube_volumes, tube_volume_after_fragment(record, index, fragment_map, length(objects)))
    end
    lower = (-1.0, -1.0, -0.55)
    upper = (1.0, 1.0, 0.6)
    # OCC bounding boxes carry a gap of ~1e-7; classify with a looser tolerance.
    tolerance = 1.0e-5
    on_box(bounds) = any(abs(bounds[d] - lower[d]) < tolerance && abs(bounds[d + 3] - lower[d]) < tolerance ||
                         abs(bounds[d] - upper[d]) < tolerance && abs(bounds[d + 3] - upper[d]) < tolerance
                         for d in 1:3)
    substrate_set = Set(substrate_tags)
    vacuum_set = Set(vacuum_tags)
    groups = Dict{Int, Vector{Int32}}()
    for (dim, tag) in gmsh.model.getEntities(2)
        up, _ = gmsh.model.getAdjacencies(dim, tag)
        bounds = gmsh.model.getBoundingBox(dim, tag)
        if on_box(bounds)
            # Top and bottom box faces (touching only vacuum / only substrate) are
            # the terminal (7); the lateral box faces, which the metal ridge meets,
            # stay the natural matching surface (1) so that no potential jump
            # (infinite energy) is prescribed along a metal/box contact line.
            horizontal = abs(bounds[3] - bounds[6]) < tolerance
            push!(get!(groups, horizontal ? 7 : 1, Int32[]), tag)
            continue
        end
        in_substrate = count(v -> v in substrate_set, up)
        in_vacuum = count(v -> v in vacuum_set, up)
        attribute = if in_substrate > 0 && in_vacuum > 0
            abs(bounds[3]) < tolerance && abs(bounds[6]) < tolerance ? 3000 : 3100
        elseif length(up) == 1 && in_substrate == 1
            5001
        elseif length(up) == 1 && in_vacuum == 1
            6001
        else
            0   # interior face between two volumes of one material
        end
        attribute > 0 && push!(get!(groups, attribute, Int32[]), tag)
    end
    for (attribute, surfaces) in sort(collect(groups))
        gmsh.model.addPhysicalGroup(2, surfaces, attribute,
                                    attribute == 1 ? "matching_surface" :
                                    attribute == 7 ? "terminal" : "surface_$attribute")
    end
    # Explicit tube meshes, phase 1 (points and curves).
    next_node = Ref(0)
    point_nodes = Dict{Int32, Int}()
    states = TubeMesh[]
    for (tube, section) in tubes
        volumes = Tuple{Int32, Tuple{Int, Int, Int}, Dict}[]
        for (index, record) in enumerate(records)
            record.tube === tube || continue
            volume = tube_volumes[index]
            push!(volumes, (volume, record.group,
                            match_tube_entities(volume, tube, section, record.group, 1.0e-6)))
        end
        state = TubeMesh(tube, section, volumes; pyramid_height=pyramid_height)
        install_tube_curves!(state, next_node, point_nodes)
        push!(states, state)
    end
    if !isempty(records)
        gmsh.option.setNumber("Mesh.MeshOnlyEmpty", 1)
        # generate() renumbers nodes by default; the explicit tube node tags must
        # stay valid across the 2D and 3D passes.
        gmsh.option.setNumber("Mesh.Renumber", 0)
    end
    # Size field: lc_edge at the two metal edges growing linearly to lc_far.
    edge_curves = Int32[]
    for (dim, tag) in gmsh.model.getEntities(1)
        bounds = gmsh.model.getBoundingBox(dim, tag)
        abs(bounds[1]) < tolerance && abs(bounds[4]) < tolerance &&
            (abs(bounds[3] - bounds[6]) < tolerance) && (abs(bounds[3]) < tolerance ||
                                                         abs(bounds[3] - 0.1) < tolerance) &&
            push!(edge_curves, tag)
    end
    isempty(edge_curves) && error("no metal edge curves found")
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(edge_curves))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 200)
    gmsh.model.mesh.field.add("MathEval", 2)
    gmsh.model.mesh.field.setString(2, "F", "min($lc_far,$lc_edge+($lc_far-$lc_edge)*F1/0.4)")
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    for (name, value) in [("Mesh.MeshSizeMin", min(lc_edge, inner_size)), ("Mesh.MeshSizeMax", lc_far),
                          ("Mesh.Algorithm3D", 1), ("Mesh.MeshSizeExtendFromBoundary", 0),
                          ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0),
                          ("Mesh.MshFileVersion", 2.2), ("Mesh.Binary", 1)]
        gmsh.option.setNumber(name, value)
    end
    elapsed = @elapsed begin
        gmsh.model.mesh.generate(2)
        tube_faces = [install_tube_faces!(state, next_node) for state in states]
        remove_tube_volumes!(states)
        gmsh.model.mesh.generate(3)
    end
    duplicates = remove_duplicate_volume_elements!()
    tube_discrete, tube_census = finalize_tube_volumes!(states)
    outside_substrate = setdiff(substrate_tags, tube_volumes)
    outside_vacuum = setdiff(vacuum_tags, tube_volumes)
    gmsh.model.addPhysicalGroup(3, vcat(outside_substrate, get(tube_discrete, 1, Int32[])), 1,
                                "substrate")
    gmsh.model.addPhysicalGroup(3, vcat(outside_vacuum, get(tube_discrete, 2, Int32[])), 2,
                                "vacuum")
    census = hybrid_volume_census()
    census["Areas"] = physical_surface_areas()
    census["Tubes"] = tube_census
    census["TubeFaces"] = tube_faces
    census["GenerateSeconds"] = elapsed
    census["DuplicateVolumeElementsRemoved"] = duplicates
    census["Kind"] = kind
    isempty(records) ||
        (census["Anisotropy"] = tube_anisotropy(tubes[1][1], tubes[1][2], pyramid_height))
    gmsh.write(filename)
    gmsh.finalize()
    if census_path !== nothing
        open(census_path, "w") do io
            write_census_json(io, census)
        end
    end
    return census
end

function parse_tiny_options(args)
    options = Dict{Symbol, Any}()
    index = 3
    while index <= length(args)
        if args[index] == "--inner-size"
            options[:inner_size] = parse(Float64, args[index + 1])
        elseif args[index] == "--spacing"
            options[:spacing] = parse(Float64, args[index + 1])
        elseif args[index] == "--lc-edge"
            options[:lc_edge] = parse(Float64, args[index + 1])
        elseif args[index] == "--pyramid-height"
            options[:pyramid_height] = parse(Float64, args[index + 1])
        elseif args[index] == "--census"
            options[:census_path] = abspath(args[index + 1])
        else
            error("unknown option $(args[index])")
        end
        index += 2
    end
    return options
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("usage: mesh_tiny_hybrid.jl hybrid|tets output.msh [--inner-size um] [--spacing um] [--lc-edge um] [--pyramid-height um] [--census json]")
    census = tiny_hybrid_mesh(ARGS[1], abspath(ARGS[2]); parse_tiny_options(ARGS)...)
    for (name, value) in sort(collect(census); by=first)
        name in ("Areas", "Tubes") && continue
        println(name, ": ", value)
    end
    println("Areas: ", census["Areas"])
end
