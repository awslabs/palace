# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

Generates a small rectangular waveguide mesh for validating PML reflection against
a wave-port S11 measurement.

## Geometry

    z=0      z=L       z=L+d
    |--------|----------|
    |  guide |   PML    |
    |  (1)   |   (2)    |
    |--------|----------|
    ↑
    wave port

- Guide cross-section: a × b = 0.1 × 0.05 m (standard rectangular waveguide style).
- Guide length:        L = 0.3 m (≈ 2.6 guided wavelengths at 3 GHz).
- PML thickness:       d = 0.1 m (≈ 0.87 guided wavelengths at 3 GHz).

## Physical groups (mesh attributes)

  Domain   1  — guide interior (physical region)
  Domain   2  — PML slab on the +z face
  Boundary 3  — wave port on z = 0 (where S11 is measured)
  Boundary 4  — PEC lateral walls (±x, ±y, full length including PML)
  Boundary 5  — PEC termination on the outer +z face of the PML

## Usage

```bash
julia -e 'include("mesh.jl"); generate_waveguide_mesh(; filename="waveguide.msh")'
```
=#

using Gmsh: gmsh

extract_tag(object) = extract_tag(only(object))
extract_tag(object::Tuple) = last(object)
extract_tag(object::Integer) = error("pass a tuple, not integer")

bbox(x::Tuple) = gmsh.model.occ.get_bounding_box(x...)
xmin(x) = bbox(x)[1];
ymin(x) = bbox(x)[2];
zmin(x) = bbox(x)[3]
xmax(x) = bbox(x)[4];
ymax(x) = bbox(x)[5];
zmax(x) = bbox(x)[6]

"""
    generate_waveguide_mesh(; filename, a, b, L, d, n_pml_slabs, mesh_size, verbose, gui)

Generate the rectangular waveguide + PML slab mesh.

# Arguments

  - filename    - output .msh filename (written in this directory)
  - a, b        - waveguide cross-section dimensions (m)
  - L           - guide length (m)
  - d           - total PML thickness (m)
  - n_pml_slabs - number of axis-aligned PML slabs along +z. Each slab becomes its
    own mesh attribute so Palace's per-attribute PML stretch tensor
    discretizes the continuous σ(z) at `n_pml_slabs` successive
    levels. More slabs ⇒ closer to the ideal continuously-varying
    σ(z) ⇒ lower residual |S11|. Attributes assigned in z-order,
    starting at 2 (slab closest to the physical region) through
    `1 + n_pml_slabs` (outermost slab).
  - mesh_size   - target element size (m)
  - verbose     - gmsh verbosity (0-5)
  - gui         - open gmsh GUI after mesh generation
"""
function generate_waveguide_mesh(;
    filename::AbstractString="waveguide.msh",
    a::Real=0.1,
    b::Real=0.05,
    L::Real=0.3,
    d::Real=0.1,
    n_pml_slabs::Integer=1,
    mesh_size::Real=0.015,
    verbose::Integer=3,
    gui::Bool=false
)
    @assert n_pml_slabs >= 1 "n_pml_slabs must be >= 1"
    gmsh.initialize()
    kernel = gmsh.model.occ
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "waveguide" in gmsh.model.list()
        gmsh.model.setCurrent("waveguide")
        gmsh.model.remove()
    end
    gmsh.model.add("waveguide")

    # Build the guide + N PML sub-slabs of equal thickness d/N along +z. Fragment
    # glues all of them so shared faces are conformal.
    slab_dz = d / n_pml_slabs
    guide = kernel.addBox(-a/2, -b/2, 0.0, a, b, L)
    pml_slab_tags = [
        kernel.addBox(-a/2, -b/2, L + k*slab_dz, a, b, slab_dz) for k = 0:(n_pml_slabs - 1)
    ]
    kernel.fragment([(3, guide)], [(3, tag) for tag in pml_slab_tags])
    kernel.synchronize()

    # Identify all domains by centroid z. The guide has z_mid < L; each PML slab k
    # has z_mid ≈ L + (k + 0.5) * slab_dz.
    all_3d = kernel.getEntities(3)
    @assert length(all_3d) == 1 + n_pml_slabs
    zmid(x) = 0.5 * (zmin(x) + zmax(x))
    guide_dimtag = filter(e -> zmid(e) < L, all_3d) |> only
    # Sort PML slabs by z from inner (closest to physical) to outer.
    pml_dimtags = sort(filter(e -> zmid(e) > L, all_3d); by=zmid)
    @assert length(pml_dimtags) == n_pml_slabs

    # Boundary faces.
    all_2d = kernel.getEntities(2)
    eps = mesh_size / 100

    # Wave port on z = 0: degenerate in z, zmax ≈ 0.
    port_faces = filter(e -> (zmax(e) - zmin(e) < eps) && abs(zmax(e)) < eps, all_2d)
    @assert length(port_faces) == 1

    # PEC termination on z = L + d.
    pec_end_faces =
        filter(e -> (zmax(e) - zmin(e) < eps) && abs(zmax(e) - (L + d)) < eps, all_2d)
    @assert length(pec_end_faces) == 1

    # Lateral PEC walls: all remaining faces that are degenerate in x or y.
    lateral_faces = filter(e -> begin
        dx = xmax(e) - xmin(e);
        dy = ymax(e) - ymin(e);
        dz = zmax(e) - zmin(e)
        # Faces on ±x or ±y walls have dx or dy ≈ 0.
        (dx < eps) || (dy < eps)
    end, all_2d)
    # Exclude any accidental internal interface faces (z-extent should be nontrivial).
    lateral_faces = filter(e -> (zmax(e) - zmin(e)) > eps, lateral_faces)

    # Create physical groups. Attributes come out in creation order:
    #   1         = guide
    #   2 .. 1+N  = PML slabs, inner-to-outer
    #   next      = port, pec_lateral, pec_end (boundary attributes)
    guide_attr = gmsh.model.addPhysicalGroup(3, [extract_tag(guide_dimtag)], -1, "guide")
    pml_slab_attrs = Int[]
    for (k, dt) in enumerate(pml_dimtags)
        push!(
            pml_slab_attrs,
            gmsh.model.addPhysicalGroup(3, [extract_tag(dt)], -1, "pml_slab_$k")
        )
    end
    port_attr = gmsh.model.addPhysicalGroup(2, extract_tag.(port_faces), -1, "port")
    pec_lateral_attr =
        gmsh.model.addPhysicalGroup(2, extract_tag.(lateral_faces), -1, "pec_lateral")
    pec_end_attr =
        gmsh.model.addPhysicalGroup(2, extract_tag.(pec_end_faces), -1, "pec_end")

    # Uniform mesh size.
    gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.CharacteristicLengthFromCurvature", 0)

    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    gmsh.model.mesh.generate(3)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\n=== Mesh generated: $filename ($n_pml_slabs PML slab(s)) ===")
    println("  guide:       attr $guide_attr, bbox z ∈ [0, $L]")
    for (k, attr) in enumerate(pml_slab_attrs)
        z0 = L + (k - 1) * slab_dz
        z1 = L + k * slab_dz
        println(
            "  pml_slab_$k:  attr $attr, bbox z ∈ [$(round(z0; digits=4)), $(round(z1; digits=4))]"
        )
    end
    println("  port:        attr $port_attr (2D, z = 0)")
    println("  pec_lateral: attr $pec_lateral_attr (2D, ±x, ±y walls, entire length)")
    println("  pec_end:     attr $pec_end_attr (2D, z = $(L+d))")
    println()

    if gui
        gmsh.fltk.run()
    end
    return gmsh.finalize()
end
