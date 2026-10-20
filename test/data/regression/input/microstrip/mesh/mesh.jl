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

Generate a mesh for the microstrip regression cases using Gmsh: two equal-width plates
(strip and ground plane) joined at both ends by shorting walls, in a box of air. A short
terminal patch on the top strip carries the current excitation of the magnetostatic cases.
The parallel-plate limit `separation << plate_width` yields the analytic loop inductance
`L ≈ μ₀ * separation * strip_length / plate_width`.

The planes of the shorting walls, x = ±strip_length/2, split the air beyond the ends from the
air around the line, and the gap between the plates is its own volume, so one mesh serves
two problems:

  - Magnetostatic (μm): every volume is air, and the plates and walls are conductors (PEC
    or Superconductor) closing a current loop driven at the terminal.
  - Driven with wave ports: the gap holds the substrate, and the air beyond the ends has no
    material, so it is removed when the mesh is loaded and the end planes become the ports.
    Their elements are oriented with the normal pointing away from the line, out of the
    remaining domain, as the wave ports require.

Volume attributes: 1 air around the line, 2 gap between the plates, 3 air beyond the ends.
Boundary attributes: 2 farfield box, 3 plates, 4 terminal, 5 and 6 the walls at x =
-strip_length/2 and x = +strip_length/2, and 7 and 8 the rest of those end planes.

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
    strip_length::Real    = 100.0,
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

    # External box, and the end planes x = ±hl through it
    x0, y0, z0 = -hl - farfield_pad, -farfield_pad, -farfield_pad
    x1, y1, z1 = hl + farfield_pad, plate_width + farfield_pad, separation + farfield_pad
    domain = kernel.addBox(x0, y0, z0, x1 - x0, y1 - y0, z1 - z0)
    function end_plane(x)
        p = [
            kernel.addPoint(x, y0, z0),
            kernel.addPoint(x, y1, z0),
            kernel.addPoint(x, y1, z1),
            kernel.addPoint(x, y0, z1)
        ]
        return quad(p...)
    end
    plane_m, plane_p = end_plane(-hl), end_plane(hl)

    # The gap between the plates is a volume, bounded by the plates and walls.
    gap = kernel.addBox(-hl, 0.0, -hz, strip_length, plate_width, separation)

    # Imprint everything into the domain volume so the mesh conforms to it.
    kernel.fragment(
        [(3, domain)],
        vcat(
            [(3, gap)],
            [
                (2, s) for s in
                [top_left, terminal, top_right, bottom, wall_p, wall_m, plane_m, plane_p]
            ]
        )
    )
    kernel.synchronize()

    # The fragment renumbers entities, so locate the physical groups by bounding box.
    tol = 1.0e-3
    in_box(dim, x0, y0, z0, x1, y1, z1) = [
        e[2] for e in gmsh.model.getEntitiesInBoundingBox(
            x0 - tol,
            y0 - tol,
            z0 - tol,
            x1 + tol,
            y1 + tol,
            z1 + tol,
            dim
        )
    ]
    volumes = [v[2] for v in gmsh.model.getEntities(3)]
    gap_vols = in_box(3, -hl, 0.0, -hz, hl, plate_width, hz)
    end_vols = vcat(in_box(3, x0, y0, z0, -hl, y1, z1), in_box(3, hl, y0, z0, x1, y1, z1))
    air_vols = setdiff(volumes, gap_vols, end_vols)
    @assert length(gap_vols) == 1 && length(end_vols) == 2 && length(air_vols) == 1

    terminal_surfs = in_box(2, -ht, 0.0, hz, ht, plate_width, hz)
    plate_surfs = setdiff(
        vcat(
            in_box(2, -hl, 0.0, hz, hl, plate_width, hz),
            in_box(2, -hl, 0.0, -hz, hl, plate_width, -hz)
        ),
        terminal_surfs
    )
    wall_m_surfs = in_box(2, -hl, 0.0, -hz, -hl, plate_width, hz)
    wall_p_surfs = in_box(2, hl, 0.0, -hz, hl, plate_width, hz)
    plane_m_surfs = setdiff(in_box(2, -hl, y0, z0, -hl, y1, z1), wall_m_surfs)
    plane_p_surfs = setdiff(in_box(2, hl, y0, z0, hl, y1, z1), wall_p_surfs)
    conductor_surfs = vcat(plate_surfs, wall_m_surfs, wall_p_surfs)
    farfield_surfs = [
        s[2] for s in gmsh.model.getEntities(2) if begin
            bb = gmsh.model.getBoundingBox(2, s[2])
            any(
                abs(bb[i] - v) < tol && abs(bb[i + 3] - v) < tol for
                (i, v) in ((1, x0), (1, x1), (2, y0), (2, y1), (3, z0), (3, z1))
            )
        end
    ]

    # Add physical groups
    gmsh.model.addPhysicalGroup(3, air_vols, 1, "air")
    gmsh.model.addPhysicalGroup(3, gap_vols, 2, "gap")
    gmsh.model.addPhysicalGroup(3, end_vols, 3, "end")
    gmsh.model.addPhysicalGroup(2, farfield_surfs, 2, "farfield")
    gmsh.model.addPhysicalGroup(2, plate_surfs, 3, "plates")
    gmsh.model.addPhysicalGroup(2, terminal_surfs, 4, "terminal")
    gmsh.model.addPhysicalGroup(2, wall_m_surfs, 5, "wall_m")
    gmsh.model.addPhysicalGroup(2, wall_p_surfs, 6, "wall_p")
    gmsh.model.addPhysicalGroup(2, plane_m_surfs, 7, "end_m")
    gmsh.model.addPhysicalGroup(2, plane_p_surfs, 8, "end_p")

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

    # Elements on the end planes are interior faces of the full mesh, which keep the
    # orientation written here. Orient them with the normal pointing away from the line,
    # which is outward once the air beyond the ends is removed.
    for (surfs, sign) in (
        (vcat(wall_m_surfs, plane_m_surfs), -1.0),
        (vcat(wall_p_surfs, plane_p_surfs), 1.0)
    )
        for s in surfs
            _, _, node_tags = gmsh.model.mesh.getElements(2, s)
            nodes_per_element =
                length(node_tags[1]) ÷ length(gmsh.model.mesh.getElements(2, s)[2][1])
            normal_x(k) = begin
                tags =
                    node_tags[1][(nodes_per_element * (k - 1) + 1):(nodes_per_element * (k - 1) + 3)]
                p = [gmsh.model.mesh.getNode(t)[1] for t in tags]
                u, v = p[2] - p[1], p[3] - p[1]
                u[2] * v[3] - u[3] * v[2]
            end
            if sign * normal_x(1) < 0.0
                gmsh.model.mesh.reverse([(2, s)])
                _, _, node_tags = gmsh.model.mesh.getElements(2, s)
            end
            n_elements = length(node_tags[1]) ÷ nodes_per_element
            @assert all(sign * normal_x(k) > 0.0 for k = 1:n_elements)
        end
    end

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    println("\nFinished generating mesh. Physical group tags:")
    println("Domains: 1 air, 2 gap, 3 air beyond the ends")
    println("Boundaries: 2 farfield, 3 plates, 4 terminal, 5/6 walls, 7/8 end planes")
    println(
        "Loop squares (perimeter / width): ",
        (2.0 * strip_length + 2.0 * separation) / plate_width
    )
    println()

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
