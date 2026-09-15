# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Recompute response ownership on a rigidly published final mesh. Classification
# is performed in immutable source-local coordinates through the inverse transform.
using Gmsh: gmsh
using TOML
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))
include(joinpath(@__DIR__, "label_interface_patches.jl"))

function surface_assignment()
    result = Dict{UInt64,Tuple{Int,Vector{UInt64}}}()
    for (_, attribute) in gmsh.model.getPhysicalGroups(2)
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
            types, tags, nodes = gmsh.model.mesh.getElements(2, entity)
            for (type, elements, connectivity) in zip(types, tags, nodes)
                _, _, _, nnode, _, _ = gmsh.model.mesh.getElementProperties(type)
                for (i, element) in enumerate(elements)
                    haskey(result, element) && error("Surface element has multiple owners")
                    result[element] = (Int(attribute), collect(connectivity[(i-1)*nnode+1:i*nnode]))
                end
            end
        end
    end
    return result
end

function all_elements(dimension)
    result = Dict{UInt64,Tuple{Int,Vector{UInt64}}}()
    types, tags, nodes = gmsh.model.mesh.getElements(dimension)
    for (type, elements, connectivity) in zip(types, tags, nodes)
        _, _, _, nnode, _, _ = gmsh.model.mesh.getElementProperties(type)
        for (i, element) in enumerate(elements)
            result[element] = (Int(type), collect(connectivity[(i-1)*nnode+1:i*nnode]))
        end
    end
    return result
end

function main()
    length(ARGS) == 5 || error("source thin|fabricated mesh transform.csv report.csv")
    source, mesh, transform_path, report = abspath.((ARGS[1], ARGS[3], ARGS[4], ARGS[5]))
    kind = ARGS[2]
    kind in ("thin", "fabricated") || error("Kind must be thin or fabricated")
    all(isfile, (mesh, transform_path)) || error("Mesh or transform is missing")
    any(isfile, (report, report * ".quadrature.csv")) && error("Ownership outputs must be fresh")
    transform_values = only(readlines(transform_path))
    transform = parse_rigid_transform(transform_values)
    process = TOML.parsefile(joinpath(source, "process.toml"))
    process["Units"] == "um" || error("Process units must be um")
    edges = read_edges(joinpath(source, "mesh-signature.csv"))
    loops = read_boundary(joinpath(source, "plan-view-boundary.csv"))
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
