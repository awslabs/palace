# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    function generate_coplanar_waveguide_lumped_mesh(;
        refinement::Integer = 0,
        trace_width_μm::Real = 30.0,
        gap_width_μm::Real = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real = 800.0,
        substrate_height_μm::Real = 500.0,
        length_μm::Real = 4000.0,
        filename::AbstractString,
        verbose::Integer = 1,
    )

Generate a mesh for the coplanar waveguide with lumped ports using Gmsh

# Arguments

  - refinement - measure of how many elements to include, 0 is least
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - length_μm - length of the waveguides, in μm
  - filename - the filename to use for the generated mesh
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
"""
function generate_coplanar_waveguide_lumped_mesh(;
    refinement::Integer       = 0,
    trace_width_μm::Real      = 30.0,
    gap_width_μm::Real        = 18.0,
    separation_width_μm::Real = 200.0,
    ground_width_μm::Real     = 800.0,
    substrate_height_μm::Real = 500.0,
    length_μm::Real           = 4000.0,
    filename::AbstractString,
    verbose::Integer=1
)
    @assert refinement >= 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert separation_width_μm > 0
    @assert ground_width_μm > 0
    @assert substrate_height_μm > 0
    @assert length_μm > 0

    kernel = gmsh.model.occ

    gmsh.initialize()
    # gmsh.option.setNumber("General.Verbosity", verbose)

    # Add model
    if "cpw" in gmsh.model.list()
        gmsh.model.setCurrent("cpw")
        gmsh.model.remove()
    end
    gmsh.model.add("cpw")

    sep_dz = 1000.0
    sep_dy = 0.5 * sep_dz

    # Mesh parameters
    l_trace = 1.5 * trace_width_μm * 2^-refinement
    max_length = min(length_μm, ground_width_μm, substrate_height_μm)
    l_farfield = 2.0 * max_length * 2^-refinement

    # Chip pattern
    dy = 0.0
    g1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, ground_width_μm)
    dy += ground_width_μm
    n1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    t1 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    g2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, separation_width_μm)
    dy += separation_width_μm
    n3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    t2 = kernel.addRectangle(0.0, dy, 0.0, length_μm, trace_width_μm)
    dy += trace_width_μm
    n4 = kernel.addRectangle(0.0, dy, 0.0, length_μm, gap_width_μm)
    dy += gap_width_μm
    g3 = kernel.addRectangle(0.0, dy, 0.0, length_μm, ground_width_μm)
    dy += ground_width_μm

    # Substrate
    substrate =
        kernel.addBox(0.0, 0.0, -substrate_height_μm, length_μm, dy, substrate_height_μm)

    # Exterior box
    domain = kernel.addBox(
        -0.5 * sep_dy,
        -sep_dy,
        -sep_dz,
        length_μm + sep_dy,
        dy + 2.0 * sep_dy,
        2.0 * sep_dz
    )
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    dy = ground_width_μm
    p1a = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2a = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p1b = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p2b = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + separation_width_μm
    p3a = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p4a = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)
    dy += gap_width_μm + trace_width_μm
    p3b = kernel.addRectangle(0.0, dy, 0.0, gap_width_μm, gap_width_μm)
    p4b = kernel.addRectangle(length_μm - gap_width_μm, dy, 0.0, gap_width_μm, gap_width_μm)

    # Embedding
    geom_dimtags = filter(x -> x[1] == 2 || x[1] == 3, kernel.getEntities())
    _, geom_map = kernel.fragment(geom_dimtags, [])

    kernel.synchronize()

    # Add physical groups
    si_domain = geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)]
    @assert length(si_domain) == 1
    si_domain = last(first(si_domain))

    air_domain = filter(
        x -> x != (3, si_domain),
        geom_map[findfirst(x -> x == (3, domain), geom_dimtags)]
    )
    @assert length(air_domain) == 1
    air_domain = last(first(air_domain))

    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "si")
    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")

    metal =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x in [(2, g1), (2, g2), (2, g3), (2, t1), (2, t2)],
                        geom_dimtags
                    )]
                )
            )
        )

    metal_group = gmsh.model.addPhysicalGroup(2, metal, -1, "metal")

    port1a = last.(geom_map[findfirst(x -> x == (2, p1a), geom_dimtags)])
    port2a = last.(geom_map[findfirst(x -> x == (2, p2a), geom_dimtags)])
    port3a = last.(geom_map[findfirst(x -> x == (2, p3a), geom_dimtags)])
    port4a = last.(geom_map[findfirst(x -> x == (2, p4a), geom_dimtags)])
    port1b = last.(geom_map[findfirst(x -> x == (2, p1b), geom_dimtags)])
    port2b = last.(geom_map[findfirst(x -> x == (2, p2b), geom_dimtags)])
    port3b = last.(geom_map[findfirst(x -> x == (2, p3b), geom_dimtags)])
    port4b = last.(geom_map[findfirst(x -> x == (2, p4b), geom_dimtags)])

    port1a_group = gmsh.model.addPhysicalGroup(2, port1a, -1, "port1a")
    port2a_group = gmsh.model.addPhysicalGroup(2, port2a, -1, "port2a")
    port3a_group = gmsh.model.addPhysicalGroup(2, port3a, -1, "port3a")
    port4a_group = gmsh.model.addPhysicalGroup(2, port4a, -1, "port4a")
    port1b_group = gmsh.model.addPhysicalGroup(2, port1b, -1, "port1b")
    port2b_group = gmsh.model.addPhysicalGroup(2, port2b, -1, "port2b")
    port3b_group = gmsh.model.addPhysicalGroup(2, port3b, -1, "port3b")
    port4b_group = gmsh.model.addPhysicalGroup(2, port4b, -1, "port4b")

    gap =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x in [(2, n1), (2, n2), (2, n3), (2, n4)],
                        geom_dimtags
                    )]
                )
            )
        )
    filter!(
        x -> !(
            x in port1a ||
            x in port2a ||
            x in port3a ||
            x in port4a ||
            x in port1b ||
            x in port2b ||
            x in port3b ||
            x in port4b
        ),
        gap
    )

    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

    farfield =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x[1] == 2 && x[2] in domain_boundary,
                        geom_dimtags
                    )]
                )
            )
        )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", 0.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gap_points = filter(
        x -> x[1] == 0,
        gmsh.model.getBoundary([(2, z) for z in gap], false, true, true)
    )
    gap_curves =
        last.(
            filter(
                x -> x[1] == 1,
                gmsh.model.getBoundary([(2, z) for z in gap], false, false, false)
            )
        )
    gmsh.model.mesh.setSize(gap_points, l_trace)

    gmsh.model.mesh.field.add("Extend", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", gap)
    gmsh.model.mesh.field.setNumber(1, "Power", 1.0)
    gmsh.model.mesh.field.setNumber(1, "DistMax", sep_dz)
    gmsh.model.mesh.field.setNumber(1, "SizeMax", l_farfield)

    gmsh.model.mesh.field.add("Distance", 2)
    gmsh.model.mesh.field.setNumbers(2, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumber(2, "Sampling", 10)

    gmsh.model.mesh.field.add("Threshold", 3)
    gmsh.model.mesh.field.setNumber(3, "InField", 2)
    gmsh.model.mesh.field.setNumber(3, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(3, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(3, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(3, "DistMax", sep_dz)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [1, 3])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 8)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(1)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Si domain: ", si_domain_group)
        println("Air domain: ", air_domain_group)
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", metal_group)
        println("Trace negative boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println("Port 3: ", port3a_group, ", ", port3b_group)
        println("Port 4: ", port4a_group, ", ", port4b_group)
        println()
    end

    # Optionally launch GUI
    if "gui" in lowercase.(ARGS)
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
