# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generated with:
# julia -e 'include("mesh.jl"); generate_parallel_plate_mesh(filename="parallel_plate.msh")'

using Gmsh: gmsh

"""
    generate_parallel_plate_mesh(;
        filename::AbstractString,
        refinement::Integer = 0,
        order::Integer      = 1,
        width_mm::Real      = 7.5,    # X, transverse to port current
        height_mm::Real     = 1.0,    # Y, port current direction (plate separation)
        length_mm::Real     = 15.0,   # Z, propagation direction
        verbose::Integer    = 5,
        gui::Bool           = false
    )

Generate the mesh for the rational-impedance verification example: a section of a
parallel-plate (TEM) line. The two plates (y = 0, y = height) are PEC; the two side walls
(x = 0, x = width) are PMC, giving a clean TEM mode with uniform E_y. The two end faces
(z = 0, z = length) are uniform lumped-port surfaces (port 1 = source, port 2 = termination
under test).

Physical group tags (used as attributes in the JSON):
1 -> dielectric (vacuum volume)
2 -> conductor  (PEC plates,  y = 0 and y = height)
3 -> sidewall   (PMC walls,   x = 0 and x = width)
4 -> port1      (z = 0)
5 -> port2      (z = length)
"""
function generate_parallel_plate_mesh(;
    filename::AbstractString,
    refinement::Integer = 0,
    order::Integer      = 1,
    width_mm::Real      = 7.5,
    height_mm::Real     = 1.0,
    length_mm::Real     = 15.0,
    verbose::Integer    = 5,
    gui::Bool           = false
)
    @assert width_mm > 0 && height_mm > 0 && length_mm > 0
    @assert refinement >= 0
    @assert order > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "parallel_plate" in gmsh.model.list()
        gmsh.model.setCurrent("parallel_plate")
        gmsh.model.remove()
    end
    gmsh.model.add("parallel_plate")

    w = width_mm
    h = height_mm
    L = length_mm

    # Single hexahedral box [0,w] x [0,h] x [0,L].
    box = kernel.addBox(0.0, 0.0, 0.0, w, h, L)
    kernel.synchronize()

    # Target element size (a few elements across the plate gap; finer with refinement).
    lc = h / (2.0 * 2.0^refinement)

    # Identify the six faces by their center of mass.
    tol       = 1.0e-6 * max(w, h, L)
    conductor = Int[]   # y = 0, y = h  -> PEC
    sidewall  = Int[]   # x = 0, x = w  -> PMC
    port1     = Int[]   # z = 0
    port2     = Int[]   # z = L
    for (dim, tag) in gmsh.model.getEntities(2)
        c = gmsh.model.occ.getCenterOfMass(dim, tag)
        if abs(c[3] - 0.0) < tol
            push!(port1, tag)
        elseif abs(c[3] - L) < tol
            push!(port2, tag)
        elseif abs(c[2] - 0.0) < tol || abs(c[2] - h) < tol
            push!(conductor, tag)
        elseif abs(c[1] - 0.0) < tol || abs(c[1] - w) < tol
            push!(sidewall, tag)
        else
            error("Unclassified face at center $(c)!")
        end
    end
    @assert length(conductor) == 2 && length(sidewall) == 2
    @assert length(port1) == 1 && length(port2) == 1

    # Physical groups (explicit tags so the JSON attributes are deterministic).
    gmsh.model.addPhysicalGroup(3, [box], 1, "dielectric")
    gmsh.model.addPhysicalGroup(2, conductor, 2, "conductor")
    gmsh.model.addPhysicalGroup(2, sidewall, 3, "sidewall")
    gmsh.model.addPhysicalGroup(2, port1, 4, "port1")
    gmsh.model.addPhysicalGroup(2, port2, 5, "port2")

    # Uniform unstructured mesh.
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh.
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("  Dielectric: 1")
        println("  Conductor (PEC): 2")
        println("  Sidewall (PMC): 3")
        println("  Port 1 (source): 4")
        println("  Port 2 (termination): 5")
        println()
    end

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
