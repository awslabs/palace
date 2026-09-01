# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Generate example mesh with:
# julia -e 'include("mesh/mesh.jl"); generate_microstrip_mesh(filename="microstrip.msh")'

using Gmsh: gmsh

"""
    generate_microstrip_mesh(;
        filename::AbstractString,
        strip_length::Real   = 100.0,
        plate_width::Real     = 40.0,
        separation::Real      = 2.0,
        terminal_length::Real = 4.0,
        l_conductor::Real     = 4.0,
        l_farfield::Real      = 60.0,
        farfield_pad::Real    = 200.0,
        verbose::Integer      = 5,
        gui::Bool             = false
    )

Generate a mesh for the microstrip example using Gmsh: a superconducting strip over an
equal-width ground plane, shorted at both ends into a closed single-turn current loop. A
short terminal patch on the top strip carries the current excitation. The conductor and
terminal surfaces are tagged separately so they can be assigned PEC (geometric-inductance
baseline) or Superconductor (finite London penetration depth) boundary conditions. The
parallel-plate limit `separation << plate_width` yields the analytic loop inductance
`L ≈ μ₀ * separation * strip_length / plate_width`.

# Arguments

  - filename - the filename to use for the generated mesh
  - strip_length - length of the strip and ground plane (x), μm
  - plate_width - width of the plates (y), μm
  - separation - gap between strip and ground plane (z), μm
  - terminal_length - length of the current-source patch on the strip (x), μm
  - l_conductor - target mesh size on the conductor surfaces
  - l_farfield - target mesh size at the farfield boundary
  - farfield_pad - padding from the loop to the farfield box, μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_microstrip_mesh(;
    filename::AbstractString,
    strip_length::Real   = 100.0,
    plate_width::Real     = 40.0,
    separation::Real      = 2.0,
    terminal_length::Real = 4.0,
    l_conductor::Real     = 4.0,
    l_farfield::Real      = 60.0,
    farfield_pad::Real    = 200.0,
    verbose::Integer      = 5,
    gui::Bool             = false
)
    kernel = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "microstrip" in gmsh.model.list()
        gmsh.model.setCurrent("microstrip")
        gmsh.model.remove()
    end
    gmsh.model.add("microstrip")

    # Geometry parameters (in μm)
    hl = 0.5 * strip_length
    hz = 0.5 * separation
    ht = 0.5 * terminal_length

    # Points: top strip (z = +hz), split into left / terminal / right along x
    a  = kernel.addPoint(-hl, 0.0, hz, l_conductor)
    t1 = kernel.addPoint(-ht, 0.0, hz, l_conductor)
    t2 = kernel.addPoint(ht, 0.0, hz, l_conductor)
    b  = kernel.addPoint(hl, 0.0, hz, l_conductor)
    d  = kernel.addPoint(-hl, plate_width, hz, l_conductor)
    t3 = kernel.addPoint(-ht, plate_width, hz, l_conductor)
    t4 = kernel.addPoint(ht, plate_width, hz, l_conductor)
    c  = kernel.addPoint(hl, plate_width, hz, l_conductor)

    # Points: ground plane (z = -hz)
    e = kernel.addPoint(-hl, 0.0, -hz, l_conductor)
    f = kernel.addPoint(hl, 0.0, -hz, l_conductor)
    g = kernel.addPoint(hl, plate_width, -hz, l_conductor)
    h = kernel.addPoint(-hl, plate_width, -hz, l_conductor)

    function quad(p1, p2, p3, p4)
        l1 = kernel.addLine(p1, p2)
        l2 = kernel.addLine(p2, p3)
        l3 = kernel.addLine(p3, p4)
        l4 = kernel.addLine(p4, p1)
        return kernel.addPlaneSurface([kernel.addCurveLoop([l1, l2, l3, l4])])
    end

    # Top strip (left / terminal / right), ground plane, and the two shorting end walls
    top_left  = quad(a, t1, t3, d)
    terminal  = quad(t1, t2, t4, t3)
    top_right = quad(t2, b, c, t4)
    bottom    = quad(e, f, g, h)
    wall_p    = quad(b, f, g, c)
    wall_m    = quad(a, e, h, d)

    conductor_surfaces = [top_left, top_right, bottom, wall_p, wall_m]
    embedded_surfaces  = vcat(conductor_surfaces, terminal)

    # Add external box
    domain = kernel.addBox(
        -hl - farfield_pad,
        -farfield_pad,
        -farfield_pad,
        strip_length + 2.0 * farfield_pad,
        plate_width + 2.0 * farfield_pad,
        separation + 2.0 * farfield_pad
    )

    # Imprint the conductor and terminal surfaces into the domain volume so the mesh
    # conforms to them (a plain embed of the closed conductor shell trips tetgen).
    kernel.fragment([(3, domain)], [(2, s) for s in embedded_surfaces])
    kernel.synchronize()

    # The fragment renumbers surfaces, so locate the physical groups by bounding box.
    volumes = [v[2] for v in gmsh.model.getEntities(3)]
    @assert length(volumes) == 1
    domain = volumes[1]

    tol = 1.0e-3
    surfs_in(x0, y0, z0, x1, y1, z1) = [
        s[2] for s in gmsh.model.getEntitiesInBoundingBox(
            x0 - tol, y0 - tol, z0 - tol, x1 + tol, y1 + tol, z1 + tol, 2
        )
    ]

    terminal_surfs = surfs_in(-ht, 0.0, hz, ht, plate_width, hz)
    top_surfs = setdiff(surfs_in(-hl, 0.0, hz, hl, plate_width, hz), terminal_surfs)
    bottom_surfs = surfs_in(-hl, 0.0, -hz, hl, plate_width, -hz)
    wall_p_surfs = surfs_in(hl, 0.0, -hz, hl, plate_width, hz)
    wall_m_surfs = surfs_in(-hl, 0.0, -hz, -hl, plate_width, hz)
    conductor_surfs = vcat(top_surfs, bottom_surfs, wall_p_surfs, wall_m_surfs)

    _, farfield_boundaries = gmsh.model.getAdjacencies(3, domain)

    # Add physical groups
    domain_group = gmsh.model.addPhysicalGroup(3, [domain], -1, "domain")
    farfield_group = gmsh.model.addPhysicalGroup(2, farfield_boundaries, -1, "farfield")
    conductor_group = gmsh.model.addPhysicalGroup(2, conductor_surfs, -1, "conductor")
    terminal_group = gmsh.model.addPhysicalGroup(2, terminal_surfs, -1, "terminal")

    # Generate mesh (fine on the conductor, coarsening toward the farfield)
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_conductor)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", conductor_surfs)
    gmsh.model.mesh.field.setNumber(1, "Sampling", 50)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", l_conductor)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(2, "DistMin", separation)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5 * farfield_pad)
    gmsh.model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(2)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    println("\nFinished generating mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("Farfield boundaries: ", farfield_group)
    println("Conductor: ", conductor_group)
    println("Terminal: ", terminal_group)
    println("Loop squares (perimeter / width): ",
            (2.0 * strip_length + 2.0 * separation) / plate_width)
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
