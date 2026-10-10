# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate mesh from examples/pml_layered:
#   julia --project=.. -e 'include("mesh/mesh.jl"); generate_pml_layered_mesh()'

using Gmsh: gmsh

"""
    generate_pml_layered_mesh(; filename, w, b, L, d, h, order, verbose)

Generate a layered parallel-plate waveguide terminated by a PML: PEC plates at y = ±b/2,
natural (PMC) walls at x = 0 and x = w (a single element across), a substrate for y < 0
and vacuum for y > 0. The substrate/vacuum interface crosses the +z PML layer, so the PML
is only reflectionless if the stretch is the same in both PML regions.

Attributes:

  - Volume 1, 2: substrate, vacuum (physical region, 0 < z < L)
  - Volume 3, 4: substrate, vacuum (PML region, L < z < L + d)
  - Boundary 5: wave port, z = 0
  - Boundary 6: PEC plates, y = ±b/2
  - Boundary 7: PEC termination behind the PML, z = L + d
  - Boundary 8: PMC walls, x = 0 and x = w

Units are meters (Palace L0 = 1).
"""
function generate_pml_layered_mesh(;
    filename::AbstractString="mesh/pml_layered.msh",
    w::Real=0.002,
    b::Real=0.01,
    L::Real=0.04,
    d::Real=0.03,
    h::Real=0.00125,
    order::Integer=1,
    verbose::Integer=1
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("pml_layered")
    kernel = gmsh.model.occ
    tags = Int[]
    for (y0, y1) in [(-b / 2, 0.0), (0.0, b / 2)], (z0, z1) in [(0.0, L), (L, L + d)]
        push!(tags, kernel.addBox(0.0, y0, z0, w, y1 - y0, z1 - z0))
    end
    kernel.fragment([(3, t) for t in tags], [])
    kernel.synchronize()

    groups = Dict{String, Vector{Int}}()
    for (dim, tag) in gmsh.model.getEntities(3)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)
        sub = (ymin + ymax) / 2 < 0.0
        pml = (zmin + zmax) / 2 > L
        name = (pml ? "pml_" : "guide_") * (sub ? "sub" : "vac")
        push!(get!(groups, name, Int[]), tag)
    end
    for (k, name) in enumerate(["guide_sub", "guide_vac", "pml_sub", "pml_vac"])
        gmsh.model.addPhysicalGroup(3, groups[name], k, name)
    end

    port, plates, pec_end, walls = Int[], Int[], Int[], Int[]
    tol = 0.05 * h
    for (dim, tag) in
        gmsh.model.getBoundary([(3, t) for v in values(groups) for t in v], true, false)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
        if zmax < tol
            push!(port, tag)
        elseif zmin > L + d - tol
            push!(pec_end, tag)
        elseif ymax - ymin < tol
            push!(plates, tag)
        else
            push!(walls, tag)
        end
    end
    gmsh.model.addPhysicalGroup(2, port, 5, "port")
    gmsh.model.addPhysicalGroup(2, plates, 6, "pec_plates")
    gmsh.model.addPhysicalGroup(2, pec_end, 7, "pec_end")
    gmsh.model.addPhysicalGroup(2, walls, 8, "pmc_walls")

    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.option.setNumber("Mesh.MeshSizeMin", h)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, "..", filename))
    println("Wrote $filename with $(length(gmsh.model.mesh.getNodes()[1])) nodes")
    return gmsh.finalize()
end
