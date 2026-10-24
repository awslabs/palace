# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate meshes for the two-sided London film tests (from the two_sided_london
# directory):
#   julia -e 'include("mesh/mesh.jl"); generate_pec_film_mesh(); generate_hex_washer_mesh()'

using Gmsh: gmsh

"""
    generate_pec_film_mesh(; filename, ...)

Generate a tetrahedral mesh of a film with a circular hole inside a PEC box. Three edges
of the film lie on the PEC walls and the fourth is free, so the film meets the PEC along
one connected arc and the hole fluxoid stays well defined.

# Attributes

  - Volume 1: the box interior, on both sides of the film
  - Boundary 2: PEC box walls
  - Boundary 8: the film, |x| ≤ L and -L ≤ y ≤ b at z = 0
  - Boundary 9: the hole cap, centered at y = y0 on z = 0

Units: μm (Palace L0 = 1e-6).
"""
function generate_pec_film_mesh(;
    filename::AbstractString="mesh/pec_film.msh",
    L::Real=5.0,
    b::Real=3.0,
    y0::Real=-1.0,
    r::Real=1.5,
    h_min::Real=0.8,
    h_max::Real=2.5,
    verbose::Integer=0
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("pec_film")
    occ = gmsh.model.occ
    box = occ.addBox(-L, -L, -L, 2L, 2L, 2L)
    sq = occ.addRectangle(-L, -L, 0, 2L, L + b)
    disk = occ.addDisk(0, y0, 0, r, r)
    occ.fragment([(3, box)], [(2, sq), (2, disk)])
    occ.synchronize()

    vols = [t for (d, t) in gmsh.model.getEntities(3)]
    film, hole, walls = Int[], Int[], Int[]
    for (d, t) in gmsh.model.getEntities(2)
        bb = gmsh.model.getBoundingBox(d, t)
        if abs(bb[3]) < 1e-6 && abs(bb[6]) < 1e-6
            push!((bb[4] - bb[1] < 2r + 0.1) ? hole : film, t)
        else
            push!(walls, t)
        end
    end
    gmsh.model.addPhysicalGroup(3, vols, 1)
    gmsh.model.addPhysicalGroup(2, walls, 2)
    gmsh.model.addPhysicalGroup(2, film, 8)
    gmsh.model.addPhysicalGroup(2, hole, 9)

    # Refine toward the film.
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", vcat(film, hole))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", h_min)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", h_max)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.2)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 3.0)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(joinpath(@__DIR__, basename(filename)))
    return gmsh.finalize()
end

"""
    generate_hex_washer_mesh(; filename, ...)

Generate a structured hexahedral mesh of a square washer centered in a PEC box, extruded
from the film plane so the mesh is mirror-symmetric about it.

# Attributes

  - Volume 1: the box interior, on both sides of the film
  - Boundary 2: PEC box walls
  - Boundary 8: the washer, a ≥ max(|x|, |y|) ≥ b at z = 0
  - Boundary 9: the hole cap, max(|x|, |y|) < b at z = 0

Units: μm (Palace L0 = 1e-6).
"""
function generate_hex_washer_mesh(;
    filename::AbstractString="mesh/hex_washer.msh",
    L::Real=5.0,
    a::Real=3.0,
    b::Real=1.0,
    n_xy::Integer=3,
    n_z::Integer=5,
    verbose::Integer=0
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("hex_washer")
    geo = gmsh.model.geo

    # A 5 × 5 grid of transfinite quadrilateral patches on z = 0, extruded to ±L.
    xs = [-L, -a, -b, b, a, L]
    p = [geo.addPoint(x, y, 0.0) for x in xs, y in xs]
    hx = [geo.addLine(p[i, j], p[i + 1, j]) for i = 1:5, j = 1:6]
    hy = [geo.addLine(p[i, j], p[i, j + 1]) for i = 1:6, j = 1:5]
    for t in vcat(vec(hx), vec(hy))
        geo.mesh.setTransfiniteCurve(t, n_xy + 1)
    end
    patches = Tuple{Int, Int}[]
    for i = 1:5, j = 1:5
        loop = geo.addCurveLoop([hx[i, j], hy[i + 1, j], -hx[i, j + 1], -hy[i, j]])
        s = geo.addPlaneSurface([loop])
        geo.mesh.setTransfiniteSurface(s)
        geo.mesh.setRecombine(2, s)
        push!(patches, (2, s))
    end
    geo.extrude(patches, 0, 0, L, [n_z], [], true)
    geo.extrude(patches, 0, 0, -L, [n_z], [], true)
    geo.synchronize()

    vols = [t for (d, t) in gmsh.model.getEntities(3)]
    film, hole, walls = Int[], Int[], Int[]
    for (d, t) in gmsh.model.getEntities(2)
        bb = gmsh.model.getBoundingBox(d, t)
        on(k, v) = abs(bb[k] - v) < 1e-8 && abs(bb[k + 3] - v) < 1e-8
        if on(3, 0.0)
            xm, ym = 0.5 * (bb[1] + bb[4]), 0.5 * (bb[2] + bb[5])
            if max(abs(xm), abs(ym)) < b
                push!(hole, t)
            elseif max(abs(xm), abs(ym)) < a
                push!(film, t)
            end
        elseif any(on(k, v) for k = 1:3, v in (-L, L))
            push!(walls, t)
        end
    end
    gmsh.model.addPhysicalGroup(3, vols, 1)
    gmsh.model.addPhysicalGroup(2, walls, 2)
    gmsh.model.addPhysicalGroup(2, film, 8)
    gmsh.model.addPhysicalGroup(2, hole, 9)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(joinpath(@__DIR__, basename(filename)))
    return gmsh.finalize()
end
