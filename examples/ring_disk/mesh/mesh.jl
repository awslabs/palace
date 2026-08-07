# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate example mesh with:
# julia -e 'include("mesh/mesh.jl"); generate_ring_disk_mesh(filename="ring_disk.msh")'

using Gmsh: gmsh
using LinearAlgebra

"""
    generate_ring_disk_mesh(;
        filename::AbstractString,
        wire_width                         = 1.0,
        ring_radius                        = 100.0,
        disk_radius                        = 50.0,
        hole_radius                        = 10.0,
        rot_center::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
        rot_axis::AbstractVector{<:Real}   = [0.0, 0.0, 1.0],
        rot_θ::Real                        = 0.0,
        mesh_order::Integer                = 1,
        verbose::Integer                   = 5,
        gui::Bool                          = false
    )

Generate a mesh for the ring and holed disk example using Gmsh. A circular ring carrying a
surface current encircles a concentric, coplanar metallic disk with a central circular hole.
The ring provides a current excitation and the hole in the disk provides a flux loop, so the
geometry supports current-only, flux-only, and mixed current/flux magnetostatic excitations.

# Arguments

  - filename - the filename to use for the generated mesh
  - wire_width - width of the ring
  - ring_radius - radius of the ring, measured to the center of the wire
  - disk_radius - outer radius of the central disk
  - hole_radius - radius of the hole at the center of the disk
  - rot_center - center of rotation
  - rot_axis - axis of rotation
  - rot_θ - angle of rotation about rot_axis, originating at rot_center. The default of zero
    keeps the geometry in the ``z = 0`` plane so that the flux loop direction is ``+Z``
  - mesh_order - polynomial order of the mesh geometry
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_ring_disk_mesh(;
    filename::AbstractString,
    wire_width                         = 1.0,
    ring_radius                        = 100.0,
    disk_radius                        = 50.0,
    hole_radius                        = 10.0,
    rot_center::AbstractVector{<:Real} = [0.0, 0.0, 0.0],
    rot_axis::AbstractVector{<:Real}   = [0.0, 0.0, 1.0],
    rot_θ::Real                        = 0.0,
    mesh_order::Integer                = 1,
    verbose::Integer                   = 5,
    gui::Bool                          = false
)
    @assert 0.0 < hole_radius < disk_radius
    @assert disk_radius < ring_radius - 0.5 * wire_width

    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "ring_disk" in gmsh.model.list()
        gmsh.model.setCurrent("ring_disk")
        gmsh.model.remove()
    end
    gmsh.model.add("ring_disk")

    # Geometry parameters (in μm)
    farfield_radius = 10.0 * ring_radius

    # Mesh parameters
    l_ring = 2.0
    l_hole = 2.0
    l_disk = 5.0
    l_farfield = 200.0

    # Origin
    p0 = kernel.addPoint(0.0, 0.0, 0.0)

    # Ring, constructed as a full annulus interrupted by a small radial terminal surface
    # which carries the surface current excitation.
    h0 = 0.5 * wire_width
    r1 = ring_radius - h0
    r2 = ring_radius + h0
    x1 = sqrt(r1^2 - h0^2)
    x2 = sqrt(r2^2 - h0^2)

    pr1 = kernel.addPoint(x1, -h0, 0.0, l_ring)
    pr2 = kernel.addPoint(x1, h0, 0.0, l_ring)
    pr3 = kernel.addPoint(x2, h0, 0.0, l_ring)
    pr4 = kernel.addPoint(x2, -h0, 0.0, l_ring)

    lr1 = kernel.addLine(pr1, pr2)
    lr2 = kernel.addLine(pr2, pr3)
    lr3 = kernel.addLine(pr3, pr4)
    lr4 = kernel.addLine(pr4, pr1)

    terminal_loop = kernel.addCurveLoop([lr1, lr2, lr3, lr4])
    terminal = kernel.addPlaneSurface([terminal_loop])

    pr5 = kernel.addPoint(-r1, 0.0, 0.0, l_ring)
    pr6 = kernel.addPoint(-r2, 0.0, 0.0, l_ring)

    ar1 = kernel.addCircleArc(pr2, p0, pr5)
    ar2 = kernel.addCircleArc(pr5, p0, pr1)
    ar3 = kernel.addCircleArc(pr3, p0, pr6)
    ar4 = kernel.addCircleArc(pr6, p0, pr4)

    ring_loop = kernel.addCurveLoop([ar1, ar2, -lr4, -ar4, -ar3, -lr2])
    ring = kernel.addPlaneSurface([ring_loop])

    # Disk with a central hole, concentric and coplanar with the ring. The disk surface is
    # the metal and the hole surface spans the opening through which flux is driven.
    disk_circle = kernel.addCircle(0.0, 0.0, 0.0, disk_radius)
    disk_outer_loop = kernel.addCurveLoop([disk_circle])
    hole_circle = kernel.addCircle(0.0, 0.0, 0.0, hole_radius)
    hole_loop = kernel.addCurveLoop([hole_circle])

    disk = kernel.addPlaneSurface([disk_outer_loop, hole_loop])
    hole = kernel.addPlaneSurface([hole_loop])

    # Auxiliary surface spanning the annular gap between the disk and the inner edge of the
    # ring. Together with the disk and hole surfaces it tiles the full interior of the ring,
    # which allows the flux linked by the ring to be postprocessed.
    ring_interior_loop = kernel.addCurveLoop([ar1, ar2, lr1])
    gap = kernel.addPlaneSurface([ring_interior_loop, disk_outer_loop])

    # Add external box
    domain = kernel.addBox(
        -farfield_radius,
        -farfield_radius,
        -farfield_radius,
        2.0 * farfield_radius,
        2.0 * farfield_radius,
        2.0 * farfield_radius
    )

    # Apply a rotation transformation to all entities in the model
    if !iszero(rot_θ)
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
    end

    kernel.synchronize()

    # Add physical groups
    domain_group = gmsh.model.addPhysicalGroup(3, [domain], -1, "domain")

    _, farfield_boundaries = gmsh.model.getAdjacencies(3, domain)
    farfield_group = gmsh.model.addPhysicalGroup(2, farfield_boundaries, -1, "farfield")

    ring_group = gmsh.model.addPhysicalGroup(2, [ring], -1, "ring")
    terminal_group = gmsh.model.addPhysicalGroup(2, [terminal], -1, "terminal")
    disk_group = gmsh.model.addPhysicalGroup(2, [disk], -1, "disk")
    hole_group = gmsh.model.addPhysicalGroup(2, [hole], -1, "hole")
    gap_group = gmsh.model.addPhysicalGroup(2, [gap], -1, "gap")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", min(l_ring, l_hole))
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", [ring, terminal, disk, hole, gap])
    gmsh.model.mesh.field.setNumber(1, "Power", 1.0)
    gmsh.model.mesh.field.setNumber(1, "DistMax", 6.0 * ring_radius)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    # Refine towards the ring wire and towards the edges of the hole and disk, which carry
    # the strongest field gradients. Every threshold must grade all the way up to
    # l_farfield: a Threshold field evaluates to SizeMax beyond DistMax rather than to no
    # constraint, so a field capped at a small SizeMax would constrain the entire far-field
    # volume once combined with the Min field below.
    ring_curves = last.(
        gmsh.model.getBoundary([(2, x) for x in [ring, terminal]], true, false, false)
    )

    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", ring_curves)
    gmsh.model.mesh.field.setNumber(2, "Sampling", 30)

    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 2)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", l_ring)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "DistMax", 6.0 * ring_radius)

    gmsh.model.mesh.field.add("Distance", 4)
    gmsh.model.mesh.field.setNumbers(4, "CurvesList", [hole_circle, disk_circle])
    gmsh.model.mesh.field.setNumber(4, "Sampling", 60)

    gmsh.model.mesh.field.add("Threshold", 5)
    gmsh.model.mesh.field.setNumber(5, "InField", 4)
    gmsh.model.mesh.field.setNumber(5, "SizeMin", l_hole)
    gmsh.model.mesh.field.setNumber(5, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(5, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(5, "DistMax", 6.0 * ring_radius)

    # Hold the disk and its interior at a moderate size so the flux through the hole is well
    # resolved across the whole metal surface, not only near its edges.
    gmsh.model.mesh.field.add("Constant", 6)
    gmsh.model.mesh.field.setNumber(6, "VIn", l_disk)
    gmsh.model.mesh.field.setNumber(6, "VOut", l_farfield)
    gmsh.model.mesh.field.setNumbers(6, "SurfacesList", [disk, hole])

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [1, 3, 5, 6])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.model.mesh.embed(2, [ring, terminal, disk, hole, gap], 3, domain)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.model.mesh.generate(3)
    if mesh_order > 1
        gmsh.model.mesh.setOrder(mesh_order)
    end

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Farfield boundaries: ", farfield_group)
    println("Ring: ", ring_group)
    println("Terminal: ", terminal_group)
    println("Disk: ", disk_group)
    println("Hole: ", hole_group)
    println("Gap: ", gap_group)
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
