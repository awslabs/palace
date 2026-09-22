# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Recompute response ownership on a rigidly published final mesh. Classification
# is performed in immutable source-local coordinates through the inverse transform.
# Every consumed source file is an explicit, bound option; no directory is scanned.
using Gmsh: gmsh
using TOML
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))
include(joinpath(@__DIR__, "label_interface_patches.jl"))

# Snapshots of the mesh the audit must leave unchanged, as tag-sorted arrays per element
# type (decision 62 step 3, proposal 6): the same guarantee as a per-element dictionary -
# every element's type, tag, owner attribute and node list, independent of the order Gmsh
# returns them in - without one vector allocation per element (the 3M-element device
# coupons peaked near the 8 GiB stage bound with the dictionaries).
function sorted_by_tag(elements, connectivity, nnode)
    order = sortperm(elements)
    nodes = reshape(connectivity, Int(nnode), :)[:, order]
    return elements[order], vec(nodes)
end

function surface_assignment()
    by_type = Dict{Int,Tuple{Vector{UInt64},Vector{Int},Vector{UInt64}}}()
    for (_, attribute) in gmsh.model.getPhysicalGroups(2)
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
            types, tags, nodes = gmsh.model.mesh.getElements(2, entity)
            for (type, elements, connectivity) in zip(types, tags, nodes)
                record = get!(by_type, Int(type)) do
                    (UInt64[], Int[], UInt64[])
                end
                append!(record[1], elements)
                append!(record[2], fill(Int(attribute), length(elements)))
                append!(record[3], connectivity)
            end
        end
    end
    result = Tuple{Int,Vector{UInt64},Vector{Int},Vector{UInt64}}[]
    for type in sort!(collect(keys(by_type)))
        elements, attributes, connectivity = by_type[type]
        _, _, _, nnode, _, _ = gmsh.model.mesh.getElementProperties(type)
        order = sortperm(elements)
        sorted = elements[order]
        any(sorted[i] == sorted[i + 1] for i in 1:(length(sorted) - 1)) &&
            error("Surface element has multiple owners")
        push!(result, (type, sorted, attributes[order], vec(reshape(connectivity, Int(nnode), :)[:, order])))
    end
    all_tags = reduce(vcat, (r[2] for r in result); init=UInt64[])
    length(unique(all_tags)) == length(all_tags) || error("Surface element has multiple owners")
    return result
end

function all_elements(dimension)
    types, tags, nodes = gmsh.model.mesh.getElements(dimension)
    result = Tuple{Int,Vector{UInt64},Vector{UInt64}}[]
    for (type, elements, connectivity) in zip(types, tags, nodes)
        _, _, _, nnode, _, _ = gmsh.model.mesh.getElementProperties(type)
        sorted, sorted_nodes = sorted_by_tag(elements, connectivity, nnode)
        push!(result, (Int(type), sorted, sorted_nodes))
    end
    sort!(result; by=first)
    return result
end

const OWNERSHIP_SOURCE_OPTIONS = ("--process", "--signature", "--boundary")

function parse_arguments(args)
    usage = "thin|fabricated mesh transform.csv report.csv --process P --signature S --boundary B"
    length(args) >= 4 || error(usage)
    kind = args[1]
    kind in ("thin", "fabricated") || error("Kind must be thin or fabricated")
    mesh, transform_path, report = abspath.((args[2], args[3], args[4]))
    sources = Dict{String,String}()
    index = 5
    while index <= length(args)
        option = args[index]
        option in OWNERSHIP_SOURCE_OPTIONS || error("Unknown option $option; usage: $usage")
        index < length(args) || error("Missing value for option $option")
        haskey(sources, option) && error("Duplicate option $option")
        sources[option] = abspath(args[index + 1])
        index += 2
    end
    for option in OWNERSHIP_SOURCE_OPTIONS
        haskey(sources, option) || error("Missing required option $option")
        isfile(sources[option]) || error("Source input for $option does not exist")
    end
    return kind, mesh, transform_path, report, sources
end

function main()
    kind, mesh, transform_path, report, sources = parse_arguments(ARGS)
    all(isfile, (mesh, transform_path)) || error("Mesh or transform is missing")
    any(isfile, (report, report * ".quadrature.csv")) && error("Ownership outputs must be fresh")
    transform_values = only(readlines(transform_path))
    transform = parse_rigid_transform(transform_values)
    process = TOML.parsefile(sources["--process"])
    process["Units"] == "um" || error("Process units must be um")
    edges = read_edges(sources["--signature"])
    loops = read_boundary(sources["--boundary"])
    gmsh.initialize()
    try
        gmsh.open(mesh)
        before_nodes = gmsh.model.mesh.getNodes()[1:2]
        before_surface = surface_assignment()
        before_volume = all_elements(3)
        label_interface_patches(
            edges, loops, Float64(process["Radius"]), report;
            minimum_size = 0.0, fabricated = kind == "fabricated",
            metal_thickness = Float64(process["MetalThickness"]),
            overetch = Float64(process["Overetch"]),
            ownership_coordinates = point -> inverse_transform_point(transform, point))
        before_nodes == gmsh.model.mesh.getNodes()[1:2] || error("Ownership audit changed nodes")
        before_volume == all_elements(3) || error("Ownership audit changed volume connectivity")
        before_surface == surface_assignment() ||
            error("Final surface owner labels differ from inverse-local reconstruction")
    finally
        gmsh.isInitialized() != 0 && gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
