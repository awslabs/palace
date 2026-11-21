# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create the mesh of an n x n lattice of superconducting rings
(annular films with a flux hole each) in a box, for the magnetostatic substructuring example
(see flux_lattice_*.json). The elements in a box around the middle ring form the region
(domain attribute 1), the rest the environment (domain attribute 2).

## How to run

From this directory, run:

```bash
julia --project -e 'include("flux_lattice.jl"); generate_flux_lattice()'
```

which writes mesh/flux_lattice.msh.
=#

using Gmsh: gmsh

"""
    generate_flux_lattice(;
        n = 5,
        pitch = 10.0,
        R = 3.0,
        r = 1.0,
        box = (60.0, 60.0, 30.0),
        region = (10.0, 10.0, 10.0),
        mesh_size_fine = 0.25,
        mesh_size_medium = 0.6,
        mesh_size_coarse = 4.0,
        filename = "flux_lattice.msh",
        verbose = 1
    )

Generate the mesh of an `n` x `n` lattice of rings with spacing `pitch` on the plane z = 0,
centered in a box of size `box`. Each ring is an annular film of outer radius `R` and inner
radius `r`, with a hole surface filling the inner circle. The region is a box of size
`region` around the middle ring. Lengths are in μm.

Attributes: domains 1 (region) and 2 (environment); boundaries 1 (box walls), 10 + k (film
of ring k), h + k (hole of ring k) with h = max(40, 10 + n²), and the rings k = 1, 2, ...
numbered by rows of increasing y, then by increasing x.
"""
function generate_flux_lattice(;
    n=5,
    pitch=10.0,
    R=3.0,
    r=1.0,
    box=(60.0, 60.0, 30.0),
    region=(10.0, 10.0, 10.0),
    mesh_size_fine=0.25,
    mesh_size_medium=0.6,
    mesh_size_coarse=4.0,
    filename="flux_lattice.msh",
    verbose=1
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("flux_lattice")
    occ = gmsh.model.occ

    # Box, region box, and the films and holes on z = 0.
    outer = occ.addBox(-box[1] / 2, -box[2] / 2, -box[3] / 2, box...)
    inner = occ.addBox(-region[1] / 2, -region[2] / 2, -region[3] / 2, region...)
    films, holes, centers = Int[], Int[], Tuple{Float64, Float64}[]
    for j = 1:n, i = 1:n
        c = ((i - (n + 1) / 2) * pitch, (j - (n + 1) / 2) * pitch)
        disk = occ.addDisk(c[1], c[2], 0.0, R, R)
        hole = occ.addDisk(c[1], c[2], 0.0, r, r)
        film, _ = occ.cut([(2, disk)], [(2, hole)], -1, true, false)
        push!(films, film[1][2])
        push!(holes, hole)
        push!(centers, c)
    end
    _, map = occ.fragment(
        [(3, outer)],
        [(3, inner); [(2, s) for s in films]; [(2, s) for s in holes]]
    )
    occ.synchronize()

    # Region and environment volumes, and the surfaces of each film and hole.
    region_vol = [t for (d, t) in map[2] if d == 3]
    env_vol = setdiff([t for (d, t) in map[1] if d == 3], region_vol)
    film_tags = [[t for (d, t) in map[2 + k] if d == 2] for k = 1:(n * n)]
    hole_tags = [[t for (d, t) in map[2 + n * n + k] if d == 2] for k = 1:(n * n)]
    walls = [t for (d, t) in gmsh.model.getBoundary([(3, v) for v in env_vol]) if d == 2]
    walls = filter(walls) do t
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, abs(t))
        return any(abs.([x0, y0, z0, x1, y1, z1]) .≈ [box..., box...] ./ 2)
    end
    walls = unique(abs.(walls))

    gmsh.model.addPhysicalGroup(3, region_vol, 1, "region")
    gmsh.model.addPhysicalGroup(3, env_vol, 2, "environment")
    gmsh.model.addPhysicalGroup(2, walls, 1, "walls")
    hole_offset = max(40, 10 + n * n)
    for k = 1:(n * n)
        gmsh.model.addPhysicalGroup(2, film_tags[k], 10 + k, "film_$k")
        gmsh.model.addPhysicalGroup(2, hole_tags[k], hole_offset + k, "hole_$k")
    end

    # Mesh size: fine at the hole edges, medium at the film edges, coarse elsewhere.
    curves(tags) = unique([
        abs(t) for (d, t) in gmsh.model.getBoundary([(2, s) for s in tags]) if d == 1
    ])
    hole_curves = curves(vcat(hole_tags...))
    film_curves = setdiff(curves(vcat(film_tags...)), hole_curves)
    field(curves_list, size_min, dist_min, dist_max) = begin
        f = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f, "CurvesList", curves_list)
        g = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(g, "InField", f)
        gmsh.model.mesh.field.setNumber(g, "SizeMin", size_min)
        gmsh.model.mesh.field.setNumber(g, "SizeMax", mesh_size_coarse)
        gmsh.model.mesh.field.setNumber(g, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(g, "DistMax", dist_max)
        return g
    end
    fields = [
        field(hole_curves, mesh_size_fine, 0.2, 4.0),
        field(film_curves, mesh_size_medium, 0.5, 6.0)
    ]
    m = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(m, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(m)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)
    mesh_path = joinpath(@__DIR__, "mesh", filename)
    mkpath(dirname(mesh_path))
    gmsh.write(mesh_path)
    println("Generated mesh: $(mesh_path)")
    return gmsh.finalize()
end
