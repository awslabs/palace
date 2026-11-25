# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Mesh generator for the Nb microstrip of Talanov et al., "Propagation of picosecond pulses
# on superconducting transmission line interconnects" (2022): a strip of width w and a
# ground plane, both of thickness d, separated by a SiO₂ layer of thickness s, with vacuum
# above, in a PEC box. Half of the cross-section is modeled, with the symmetry plane at
# x = 0. Lengths are in μm.
#
# Domain attributes: 1 Nb (only without `surface`), 2 SiO₂, 3 vacuum. Boundary attributes:
# 1 PEC box, 2 symmetry plane, 3 Nb surfaces (only with `surface`), 4 back face of the
# ground plane (only without `surface`), 5 and 6 the ports at z = 0 and z = l (only in 3D).

import Gmsh: gmsh

"""
    generate_superconducting_line_mesh(; kwargs...)

Generate the 2D cross-section (`l = 0`) or the 3D line of length `l` of the microstrip.
With `surface = true` the Nb films are not meshed, and their surfaces are boundaries for a
surface impedance model: the strip is a hole, and the ground plane is the bottom of the
SiO₂ layer.
"""
function generate_superconducting_line_mesh(;
    w::Real          = 4.0,    # Strip width
    d::Real          = 0.2,    # Nb thickness
    s::Real          = 0.15,   # SiO₂ thickness
    x_box::Real      = 6.0,    # Half width of the PEC box
    h_box::Real      = 4.0,    # Height of the PEC box above the ground plane
    l::Real          = 0.0,    # Line length (3D for l > 0)
    surface::Bool    = false,  # Model the Nb films as surfaces
    lc_film::Real    = surface ? 0.05 : 0.03,  # Mesh size at the Nb surfaces
    lc_far::Real     = 1.0,    # Mesh size far from the strip
    n_layers::Int    = 8,      # Number of element layers along the line (3D)
    order::Int       = 1,      # Mesh order (the geometry is straight-sided)
    filename::String = "superconducting_line.msh",
    verbose::Int     = 0
)
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("superconducting_line")
    kernel = gmsh.model.occ

    sio2 = kernel.addRectangle(0.0, 0.0, 0.0, x_box, s)
    vacuum = kernel.addRectangle(0.0, s, 0.0, x_box, h_box - s)
    strip = kernel.addRectangle(0.0, s, 0.0, w / 2, d)
    if surface
        vacuum, _ = kernel.cut([(2, vacuum)], [(2, strip)])
        kernel.fragment([(2, sio2); vacuum], [])
    else
        ground = kernel.addRectangle(0.0, -d, 0.0, x_box, d)
        kernel.fragment([(2, sio2), (2, vacuum), (2, strip), (2, ground)], [])
    end
    kernel.synchronize()

    # Mesh size: lc_film at the Nb surfaces, growing to lc_far.
    tol = 1e-6
    films = [
        tag for (_, tag) in gmsh.model.getEntities(1) if begin
            xmin, ymin, _, xmax, ymax, _ = gmsh.model.getBoundingBox(1, tag)
            ymax > -d - tol && ymin < s + d + tol && xmax > tol && xmin < x_box - tol
        end
    ]
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", films)
    gmsh.model.mesh.field.setNumber(1, "Sampling", 200)
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_film)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", d)
    gmsh.model.mesh.field.setNumber(2, "DistMax", w)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    # Extrude the cross-section along z with structured layers for a 3D line.
    dim = (l > 0.0) ? 3 : 2
    if dim == 3
        kernel.extrude(gmsh.model.getEntities(2), 0.0, 0.0, l, [n_layers], [], true)
        kernel.synchronize()
    end

    # Domains by their center, and boundaries by their position.
    in_strip(x, y) = x < w / 2 + tol && y > s - tol && y < s + d + tol
    domains = Dict(1 => Int[], 2 => Int[], 3 => Int[])
    for (_, tag) in gmsh.model.getEntities(dim)
        x, y, _ = kernel.getCenterOfMass(dim, tag)
        push!(domains[(y < 0.0 || in_strip(x, y)) ? 1 : ((y < s) ? 2 : 3)], tag)
    end
    boundaries = Dict(attr => Int[] for attr = 1:6)
    for (_, tag) in gmsh.model.getEntities(dim - 1)
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim - 1, tag)
        xc, yc = (xmin + xmax) / 2, (ymin + ymax) / 2
        attr =
            (dim == 3 && zmax < tol) ? 5 :
            (dim == 3 && zmin > l - tol) ? 6 :
            (xmax < tol) ? 2 :
            (abs(xmin - x_box) < tol || abs(ymax - h_box) < tol) ? 1 :
            (!surface && abs(ymax + d) < tol) ? 4 :
            (surface && (abs(yc) < tol || in_strip(xc, yc))) ? 3 : 0
        attr > 0 && push!(boundaries[attr], tag)
    end
    for (attr, tags) in domains
        !isempty(tags) && gmsh.model.addPhysicalGroup(dim, tags, attr)
    end
    for (attr, tags) in boundaries
        !isempty(tags) && gmsh.model.addPhysicalGroup(dim - 1, tags, attr)
    end

    gmsh.model.mesh.generate(dim)
    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(joinpath(@__DIR__, filename))
    if verbose > 0
        println("Wrote $(joinpath(@__DIR__, filename))")
    end
    return gmsh.finalize()
end
