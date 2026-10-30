# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate the mesh for the interior far-field surface test (from the farfield_interior
# directory):
#   julia -e 'include("mesh/mesh.jl"); generate_farfield_interior_mesh()'

using Gmsh: gmsh

"""
    generate_farfield_interior_mesh(; filename, h, nb, na, verbose)

Generate a structured hexahedral mesh of the cube [-B, B]^3 (B = nb h / 2) with an
interior closed surface S, the cube [-a, a]^3 (a = na h / 2). nb and na are odd, so that
the origin, where the dipole sits, is the center of an element and no mesh line passes
through it.

# Attributes

  - Volume 4: inside S
  - Volume 1: between S and the outer boundary
  - Boundary 3: outer boundary (absorbing)
  - Boundary 5: interior surface S

Units are meters.
"""
function generate_farfield_interior_mesh(;
    filename::AbstractString="mesh/farfield_interior.msh",
    h::Real=0.1 / 9,
    nb::Integer=9,
    na::Integer=5,
    verbose::Integer=2
)
    @assert isodd(nb) && isodd(na) && na < nb
    B, a = nb * h / 2, na * h / 2
    tol = 1.0e-3 * h

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    if "farfield_interior" in gmsh.model.list()
        gmsh.model.setCurrent("farfield_interior")
        gmsh.model.remove()
    end
    gmsh.model.add("farfield_interior")

    cuts = [-B, -a, a, B]
    boxes = Dict{NTuple{3, Int}, Int}()
    for i = 1:3, j = 1:3, k = 1:3
        lo = (cuts[i], cuts[j], cuts[k])
        hi = (cuts[i + 1], cuts[j + 1], cuts[k + 1])
        boxes[(i, j, k)] = gmsh.model.occ.addBox(lo..., (hi .- lo)...)
    end
    gmsh.model.occ.fragment([(3, t) for t in values(boxes)], Tuple{Int32, Int32}[])
    gmsh.model.occ.synchronize()

    inner = boxes[(2, 2, 2)]
    outer_vols = [t for (ijk, t) in boxes if ijk != (2, 2, 2)]
    gmsh.model.addPhysicalGroup(3, [inner], 4, "inside_S")
    gmsh.model.addPhysicalGroup(3, outer_vols, 1, "vacuum")

    outer, surf = Int[], Int[]
    for (_, s) in gmsh.model.getEntities(2)
        bb = gmsh.model.getBoundingBox(2, s)
        lo, hi = bb[1:3], bb[4:6]
        flat = findfirst(d -> hi[d] - lo[d] < tol, 1:3)
        isnothing(flat) && continue
        c = 0.5 * (lo[flat] + hi[flat])
        others = [d for d = 1:3 if d != flat]
        if abs(abs(c) - B) < tol
            push!(outer, s)
        elseif abs(abs(c) - a) < tol &&
               all(lo[d] > -a - tol && hi[d] < a + tol for d in others)
            push!(surf, s)
        end
    end
    @assert length(surf) == 6
    gmsh.model.addPhysicalGroup(2, outer, 3, "outer")
    gmsh.model.addPhysicalGroup(2, surf, 5, "S")

    for (_, c) in gmsh.model.getEntities(1)
        bb = gmsh.model.getBoundingBox(1, c)
        len = maximum(bb[4:6] .- bb[1:3])
        gmsh.model.mesh.setTransfiniteCurve(c, round(Int, len / h) + 1)
    end
    for (_, s) in gmsh.model.getEntities(2)
        gmsh.model.mesh.setTransfiniteSurface(s)
        gmsh.model.mesh.setRecombine(2, s)
    end
    for (_, v) in gmsh.model.getEntities(3)
        gmsh.model.mesh.setTransfiniteVolume(v)
    end

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)
    gmsh.finalize()
    return nothing
end
