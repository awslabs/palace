# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate a coarse multi-ring mesh (N concentric wire rings, each with its own current
# terminal) for exercising the Short inactive-port operator cache with more than two ports.
#
# Generate with:
#   julia -e 'include("mesh/multiring_mesh.jl"); generate_multiring_mesh(filename="multiring.msh")'
#
# Physical group / attribute layout (matches the order below):
#   1            : domain
#   2            : farfield (PEC)
#   3            : rings (PEC, the wire ring surfaces)
#   4 .. 3+N     : terminal_1 .. terminal_N  (SurfaceCurrent ports, ordered inner->outer)

using Gmsh: gmsh
using LinearAlgebra

"""
    generate_multiring_mesh(;
        filename::AbstractString,
        radii::AbstractVector{<:Real}      = [10.0, 40.0, 100.0],
        wire_width                         = 1.0,
        rot_center::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        rot_axis::AbstractVector{<:Real}   = [0.0, 0.0, 1.0],
        rot_θ::Real                        = π / 6,
        verbose::Integer                   = 5,
        gui::Bool                          = false
    )

Generate a mesh with `length(radii)` concentric wire rings using Gmsh. Each ring carries a
small square terminal surface (the current port). The mesh is intentionally coarse so it is
cheap enough for a regression/profiling fixture. Flux-loop holes are omitted since this
fixture only exercises SurfaceCurrent ports.
"""
function generate_multiring_mesh(;
    filename::AbstractString,
    radii::AbstractVector{<:Real}      = [10.0, 40.0, 100.0],
    wire_width                         = 1.0,
    rot_center::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
    rot_axis::AbstractVector{<:Real}   = [0.0, 0.0, 1.0],
    rot_θ::Real                        = π / 6,
    verbose::Integer                   = 5,
    gui::Bool                          = false
)
    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "multiring" in gmsh.model.list()
        gmsh.model.setCurrent("multiring")
        gmsh.model.remove()
    end
    gmsh.model.add("multiring")

    outer_radius = maximum(radii)
    farfield_radius = 10.0 * outer_radius

    # Coarse mesh parameters (larger than the rings example to keep the fixture small/fast).
    l_ring = 4.0
    l_farfield = 400.0

    p0 = kernel.addPoint(0.0, 0.0, 0.0)
    h0 = 0.5 * wire_width

    # Build one wire ring (square cross-section torus surface) plus its terminal surface for
    # each requested radius. Returns (terminal_surface, ring_surface).
    function add_ring(radius)
        r1 = radius - h0
        r2 = radius + h0
        x1 = sqrt(r1^2 - h0^2)
        x2 = sqrt(r2^2 - h0^2)

        p1 = kernel.addPoint(x1, -h0, 0.0, l_ring)
        p2 = kernel.addPoint(x1, h0, 0.0, l_ring)
        p3 = kernel.addPoint(x2, h0, 0.0, l_ring)
        p4 = kernel.addPoint(x2, -h0, 0.0, l_ring)

        l1 = kernel.addLine(p1, p2)
        l2 = kernel.addLine(p2, p3)
        l3 = kernel.addLine(p3, p4)
        l4 = kernel.addLine(p4, p1)

        terminal_loop = kernel.addCurveLoop([l1, l2, l3, l4])
        terminal = kernel.addPlaneSurface([terminal_loop])

        p5 = kernel.addPoint(-r1, 0.0, 0.0, l_ring)
        p6 = kernel.addPoint(-r2, 0.0, 0.0, l_ring)

        a1 = kernel.addCircleArc(p2, p0, p5)
        a2 = kernel.addCircleArc(p5, p0, p1)
        a3 = kernel.addCircleArc(p3, p0, p6)
        a4 = kernel.addCircleArc(p6, p0, p4)

        ring_loop = kernel.addCurveLoop([a1, a2, -l4, -a4, -a3, -l2])
        ring = kernel.addPlaneSurface([ring_loop])

        return terminal, ring
    end

    terminals = Int[]
    rings = Int[]
    for r in radii
        terminal, ring = add_ring(r)
        push!(terminals, terminal)
        push!(rings, ring)
    end

    # External box.
    domain = kernel.addBox(
        -farfield_radius,
        -farfield_radius,
        -farfield_radius,
        2.0 * farfield_radius,
        2.0 * farfield_radius,
        2.0 * farfield_radius
    )

    # Rotate everything (matches the rings example so the current directions carry over).
    rot_axis ./= norm(rot_axis)
    kernel.rotate(
        kernel.getEntities(),
        rot_center[1],
        rot_center[2],
        rot_center[3],
        rot_axis[1],
        rot_axis[2],
        rot_axis[3],
        rot_θ
    )

    kernel.synchronize()

    # Physical groups (attribute numbering follows registration order).
    domain_group = gmsh.model.addPhysicalGroup(3, [domain], -1, "domain")

    _, farfield_boundaries = gmsh.model.getAdjacencies(3, domain)
    farfield_group = gmsh.model.addPhysicalGroup(2, farfield_boundaries, -1, "farfield")

    rings_group = gmsh.model.addPhysicalGroup(2, rings, -1, "rings")

    terminal_groups = Int[]
    for (i, terminal) in enumerate(terminals)
        push!(
            terminal_groups,
            gmsh.model.addPhysicalGroup(2, [terminal], -1, "terminal_$(i)")
        )
    end

    # Meshing.
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_ring)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", vcat(rings, terminals))
    gmsh.model.mesh.field.setNumber(1, "Power", 1.0)
    gmsh.model.mesh.field.setNumber(1, "DistMax", 6.0 * outer_radius)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    mesh_curves = last.(
        gmsh.model.getBoundary(
            [(2, x) for x in vcat(rings, terminals)],
            true,
            false,
            false
        )
    )

    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", mesh_curves)
    gmsh.model.mesh.field.setNumber(2, "Sampling", 20)

    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 2)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", l_ring)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "DistMax", 6.0 * outer_radius)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [1, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.model.mesh.embed(2, vcat(terminals, rings), 3, domain)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating multi-ring mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Farfield boundaries: ", farfield_group)
    println("Ring boundaries: ", rings_group)
    for (i, g) in enumerate(terminal_groups)
        println("Terminal $(i): ", g)
    end
    println()

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
