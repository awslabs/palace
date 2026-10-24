# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create the mesh of a square superconducting plate with an
n x n lattice of flux holes in a box, for the magnetostatic substructuring example (see
flux_sheet_*.json). Unlike flux_lattice.jl, all holes share one film. The elements in a box
around the middle hole form the region (domain attribute 1), the rest the environment
(domain attribute 2).

## How to run

From this directory, run:

```bash
julia --project -e 'include("flux_sheet.jl"); generate_flux_sheet()'
```

which writes mesh/flux_sheet.msh.
=#

using Gmsh: gmsh

"""
    generate_flux_sheet(;
        n = 5,
        pitch = 10.0,
        r = 1.0,
        plate = 50.0,
        box = (60.0, 60.0, 30.0),
        region = (10.0, 10.0, 10.0),
        mesh_size_fine = 0.25,
        mesh_size_film = 1.0,
        mesh_size_edge = 0.6,
        mesh_size_coarse = 4.0,
        filename = "flux_sheet.msh",
        verbose = 1
    )

Generate the mesh of a square plate of side `plate` on the plane z = 0 with an `n` x `n`
lattice of circular holes of radius `r` and spacing `pitch`, centered in a box of size
`box`. Each hole has a hole surface filling it. The region is a box of size `region`
around the middle hole. Lengths are in μm.

Attributes: domains 1 (region) and 2 (environment); boundaries 1 (box walls), 10 (the
plate), 40 + k (hole k), with the holes k = 1, 2, ... numbered by rows of increasing y,
then by increasing x.
"""
function generate_flux_sheet(;
    n=5,
    pitch=10.0,
    r=1.0,
    plate=50.0,
    box=(60.0, 60.0, 30.0),
    region=(10.0, 10.0, 10.0),
    mesh_size_fine=0.25,
    mesh_size_film=1.0,
    mesh_size_edge=0.6,
    mesh_size_coarse=4.0,
    filename="flux_sheet.msh",
    verbose=1
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("flux_sheet")
    occ = gmsh.model.occ

    # Box, region box, and the plate with its holes on z = 0.
    outer = occ.addBox(-box[1] / 2, -box[2] / 2, -box[3] / 2, box...)
    inner = occ.addBox(-region[1] / 2, -region[2] / 2, -region[3] / 2, region...)
    sheet = occ.addRectangle(-plate / 2, -plate / 2, 0.0, plate, plate)
    holes = Int[]
    for j = 1:n, i = 1:n
        c = ((i - (n + 1) / 2) * pitch, (j - (n + 1) / 2) * pitch)
        push!(holes, occ.addDisk(c[1], c[2], 0.0, r, r))
    end
    film, _ = occ.cut([(2, sheet)], [(2, h) for h in holes], -1, true, false)
    _, map = occ.fragment(
        [(3, outer)],
        [(3, inner); [(d, t) for (d, t) in film]; [(2, h) for h in holes]]
    )
    occ.synchronize()

    # Region and environment volumes, the plate surfaces (split by the region box), and
    # the surface of each hole.
    nf = length(film)
    region_vol = [t for (d, t) in map[2] if d == 3]
    env_vol = setdiff([t for (d, t) in map[1] if d == 3], region_vol)
    film_tags = unique(vcat([[t for (d, t) in map[2 + k] if d == 2] for k = 1:nf]...))
    hole_tags = [[t for (d, t) in map[2 + nf + k] if d == 2] for k = 1:(n * n)]
    walls = [t for (d, t) in gmsh.model.getBoundary([(3, v) for v in env_vol]) if d == 2]
    walls = filter(walls) do t
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, abs(t))
        return any(abs.([x0, y0, z0, x1, y1, z1]) .≈ [box..., box...] ./ 2)
    end
    walls = unique(abs.(walls))

    gmsh.model.addPhysicalGroup(3, region_vol, 1, "region")
    gmsh.model.addPhysicalGroup(3, env_vol, 2, "environment")
    gmsh.model.addPhysicalGroup(2, walls, 1, "walls")
    gmsh.model.addPhysicalGroup(2, film_tags, 10, "film")
    for k = 1:(n * n)
        gmsh.model.addPhysicalGroup(2, hole_tags[k], 40 + k, "hole_$k")
    end

    # Mesh size: fine at the hole edges, medium at the plate edge, the film size on the
    # plate, coarse elsewhere.
    curves(tags) = unique([
        abs(t) for (d, t) in gmsh.model.getBoundary([(2, s) for s in tags]) if d == 1
    ])
    hole_curves = curves(vcat(hole_tags...))
    edge_curves = setdiff(curves(film_tags), hole_curves)
    threshold(f, size_min, dist_min, dist_max) = begin
        g = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(g, "InField", f)
        gmsh.model.mesh.field.setNumber(g, "SizeMin", size_min)
        gmsh.model.mesh.field.setNumber(g, "SizeMax", mesh_size_coarse)
        gmsh.model.mesh.field.setNumber(g, "DistMin", dist_min)
        gmsh.model.mesh.field.setNumber(g, "DistMax", dist_max)
        return g
    end
    distance(key, tags) = begin
        f = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f, key, tags)
        return f
    end
    fields = [
        threshold(distance("CurvesList", hole_curves), mesh_size_fine, 0.2, 3.0),
        threshold(distance("CurvesList", edge_curves), mesh_size_edge, 0.5, 4.0),
        threshold(distance("SurfacesList", film_tags), mesh_size_film, 0.0, 8.0)
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
