# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create the mesh of the line-fed rectangular microstrip patch
antenna of Sheen et al., IEEE Trans. Microw. Theory Tech. 38 (1990) 849-857.

A 12.45 mm x 16 mm patch sits on a 0.794 mm thick substrate of relative permittivity 2.2,
fed by a 2.46 mm wide microstrip line whose edge is 2.09 mm from the edge of the patch. The
ground plane is the bottom of the box. The feed line runs to the y = 0 face, which carries
a wave port around the line. The other faces cut through the substrate and are open.

The mesh has the following attributes:

  - Volume 1: substrate (0 < z < h)
  - Volume 2: air
  - Boundary 3: ground plane (z = 0)
  - Boundary 4: patch and feed line (z = h)
  - Boundary 5: wave port (y = 0, around the feed line)
  - Boundary 6: the rest of the y = 0 face
  - Boundary 7: the other faces (x = ±X, y = Y, z = Z)

Lengths are in mm (Palace "L0": 1.0e-3).

## Prerequisites

This script requires the Gmsh Julia package. Install the packages of the examples with

```bash
julia --project=examples -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From the repository root, run:
```bash
julia --project=examples -e 'include("examples/patch_antenna/mesh/mesh.jl"); generate_patch_antenna_mesh()'
```
To visualize the mesh in Gmsh's graphical interface, add the `gui=true` parameter.
=#

using Gmsh: gmsh

"""
    generate_patch_antenna_mesh(;
        filename::AbstractString="patch_antenna.msh",
        refinement::Real=1.0,
        margin::Real=6.0,
        port_width::Real=16.0,
        port_height::Real=5.0,
        size_open_edge::Real=1.0,
        dist_open_edge::Real=4.0,
        order::Integer=1,
        verbose::Integer=5,
        gui::Bool=false
    )

Generate the mesh of the patch antenna, written to `filename` in this directory.

# Arguments

  - filename - name of the mesh file
  - refinement - the element sizes are divided by this factor
  - margin - distance from the patch and the feed line to the open faces of the box, in mm
  - port_width - width of the wave port, centered on the feed line, in mm
  - port_height - height of the wave port above the ground plane, in mm
  - size_open_edge - element size at the edges between the open faces, in mm, graded to the
    size in the air over `dist_open_edge`
  - dist_open_edge - distance from the edges between the open faces over which the elements
    grow to the size in the air, in mm
  - order - polynomial order of the mesh geometry
  - verbose - Gmsh verbosity level
  - gui - whether to launch the Gmsh GUI
"""
function generate_patch_antenna_mesh(;
    filename::AbstractString="patch_antenna.msh",
    refinement::Real=1.0,
    margin::Real=6.0,
    port_width::Real=16.0,
    port_height::Real=5.0,
    size_open_edge::Real=1.0,
    dist_open_edge::Real=4.0,
    order::Integer=1,
    verbose::Integer=5,
    gui::Bool=false
)
    # Geometry of Sheen et al. (1990), Fig. 4, in mm.
    h = 0.794          # Substrate thickness
    L1 = 12.45         # Patch, along x (the edge with the feed line)
    L2 = 16.0          # Patch, along y (along the feed line)
    w = 2.46           # Feed line width
    offset = 2.09      # From the edge of the patch to the edge of the feed line
    line_length = 20.0 # From the port to the patch: 50 cells of 0.4 mm in the paper

    # The patch is centered on x = 0; the feed line runs from y = 0 to the patch.
    x_line = -L1 / 2 + offset
    y_patch = line_length
    X = L1 / 2 + margin
    Y = y_patch + L2 + margin
    Z = h + margin
    x_port = x_line + w / 2 - port_width / 2
    @assert -X < x_port && x_port + port_width < X && h < port_height < Z

    # Element sizes: at the edges of the metal, in the substrate and in the air far from the
    # metal, graded between 0.5 and 6 mm from the edges of the metal.
    size_edge = 0.5 / refinement
    size_substrate = 2.0 / refinement
    size_air = 4.0 / refinement

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    if "patch_antenna" in gmsh.model.list()
        gmsh.model.setCurrent("patch_antenna")
        gmsh.model.remove()
    end
    gmsh.model.add("patch_antenna")
    occ = gmsh.model.occ

    box = occ.addBox(-X, 0.0, 0.0, 2X, Y, Z)
    substrate = occ.addBox(-X, 0.0, 0.0, 2X, Y, h)
    patch = occ.addRectangle(-L1 / 2, y_patch, h, L1, L2)
    line = occ.addRectangle(x_line, 0.0, h, w, y_patch)
    port = occ.addRectangle(x_port, 0.0, 0.0, port_width, port_height)
    occ.rotate([(2, port)], x_port, 0.0, 0.0, 1.0, 0.0, 0.0, π / 2)
    occ.fragment([(3, box), (3, substrate)], [(2, patch), (2, line), (2, port)])
    occ.synchronize()

    # Identify the entities by their bounding boxes.
    ϵ = 1.0e-6
    function entities_in(dim, xmin, ymin, zmin, xmax, ymax, zmax)
        return last.(
            gmsh.model.getEntitiesInBoundingBox(
                xmin - ϵ,
                ymin - ϵ,
                zmin - ϵ,
                xmax + ϵ,
                ymax + ϵ,
                zmax + ϵ,
                dim
            )
        )
    end
    substrate_vols = entities_in(3, -X, 0.0, 0.0, X, Y, h)
    air_vols = setdiff(entities_in(3, -X, 0.0, 0.0, X, Y, Z), substrate_vols)
    ground = entities_in(2, -X, 0.0, 0.0, X, Y, 0.0)
    metal = vcat(
        entities_in(2, -L1 / 2, y_patch, h, L1 / 2, y_patch + L2, h),
        entities_in(2, x_line, 0.0, h, x_line + w, y_patch, h)
    )
    port_surfs = entities_in(2, x_port, 0.0, 0.0, x_port + port_width, 0.0, port_height)
    port_face = setdiff(entities_in(2, -X, 0.0, 0.0, X, 0.0, Z), port_surfs)
    open_faces = vcat(
        entities_in(2, -X, 0.0, 0.0, -X, Y, Z),
        entities_in(2, X, 0.0, 0.0, X, Y, Z),
        entities_in(2, -X, Y, 0.0, X, Y, Z),
        entities_in(2, -X, 0.0, Z, X, Y, Z)
    )
    @assert length(metal) == 2 && length(port_surfs) == 2

    gmsh.model.addPhysicalGroup(3, substrate_vols, 1, "substrate")
    gmsh.model.addPhysicalGroup(3, air_vols, 2, "air")
    gmsh.model.addPhysicalGroup(2, ground, 3, "ground")
    gmsh.model.addPhysicalGroup(2, metal, 4, "metal")
    gmsh.model.addPhysicalGroup(2, port_surfs, 5, "port")
    gmsh.model.addPhysicalGroup(2, port_face, 6, "port_face")
    gmsh.model.addPhysicalGroup(2, open_faces, 7, "open")

    # Mesh size: fine at the edges of the metal, graded to the substrate and air sizes.
    metal_curves =
        unique(last.(gmsh.model.getBoundary([(2, s) for s in metal], false, false)))
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", metal_curves)
    gmsh.model.mesh.field.setNumber(1, "Sampling", 200)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", size_edge)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", size_air)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.5)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 6.0)
    gmsh.model.mesh.field.add("Box", 3)
    gmsh.model.mesh.field.setNumber(3, "VIn", size_substrate)
    gmsh.model.mesh.field.setNumber(3, "VOut", size_air)
    gmsh.model.mesh.field.setNumber(3, "XMin", -X)
    gmsh.model.mesh.field.setNumber(3, "XMax", X)
    gmsh.model.mesh.field.setNumber(3, "YMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "YMax", Y)
    gmsh.model.mesh.field.setNumber(3, "ZMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "ZMax", h)
    # Fine along the edges between the open faces, where the PML extruded from them has its
    # edge regions, filled by the columns of the boundary elements next to the edges.
    open_edges = vcat(
        entities_in(1, -X, 0.0, Z, -X, Y, Z),
        entities_in(1, X, 0.0, Z, X, Y, Z),
        entities_in(1, -X, Y, 0.0, -X, Y, Z),
        entities_in(1, X, Y, 0.0, X, Y, Z),
        entities_in(1, -X, Y, Z, X, Y, Z)
    )
    gmsh.model.mesh.field.add("Distance", 5)
    gmsh.model.mesh.field.setNumbers(5, "CurvesList", open_edges)
    gmsh.model.mesh.field.setNumber(5, "Sampling", 200)
    gmsh.model.mesh.field.add("Threshold", 6)
    gmsh.model.mesh.field.setNumber(6, "InField", 5)
    gmsh.model.mesh.field.setNumber(6, "SizeMin", size_open_edge / refinement)
    gmsh.model.mesh.field.setNumber(6, "SizeMax", size_air)
    gmsh.model.mesh.field.setNumber(6, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(6, "DistMax", dist_open_edge)
    gmsh.model.mesh.field.add("Min", 4)
    gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3, 6])
    gmsh.model.mesh.field.setAsBackgroundMesh(4)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # Deterministic meshing algorithms.
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    n_tet = length(gmsh.model.mesh.getElementsByType(4)[1])
    println("\nFinished generating $(filename): $(n_tet) tetrahedra")
    println("Box: x = ±$(X), 0 < y < $(Y), 0 < z < $(Z) mm")
    println("Port: $(x_port) < x < $(x_port + port_width), z < $(port_height) mm")

    if gui
        gmsh.fltk.run()
    end
    return gmsh.finalize()
end
