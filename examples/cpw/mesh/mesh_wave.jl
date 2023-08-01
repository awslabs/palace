# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    function generate_coplanar_waveguide_wave_mesh(;
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

Generate a mesh for the coplanar waveguide with wave ports using Gmsh

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
function generate_coplanar_waveguide_wave_mesh(;
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
    gmsh.option.setNumber("General.Verbosity", verbose)

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
    max_length = max(length_μm, ground_width_μm, substrate_height_μm)
    l_farfield = 1.0 * max_length * 2^-refinement

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
    domain =
        kernel.addBox(0.0, -sep_dy, -sep_dz, length_μm, dy + 2.0 * sep_dy, 2.0 * sep_dz)
    _, domain_boundary = kernel.getSurfaceLoops(domain)
    @assert length(domain_boundary) == 1
    domain_boundary = first(domain_boundary)

    # Ports
    cy1 = ground_width_μm + gap_width_μm + 0.5 * trace_width_μm
    cy2 =
        cy1 +
        0.5 * trace_width_μm +
        gap_width_μm +
        separation_width_μm +
        gap_width_μm +
        0.5 * trace_width_μm
    dzp = trace_width_μm + 2.0 * gap_width_μm
    dyp = 2.0 * dzp
    let pa, pb, l
        pa = kernel.addPoint(0.0, cy1 - 0.5 * dyp, -0.5 * dzp)
        pb = kernel.addPoint(0.0, cy1 + 0.5 * dyp, -0.5 * dzp)
        l = kernel.addLine(pa, pb)
        global p1 = first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp)))[2]
    end
    let pa, pb, l
        pa = kernel.addPoint(0.0, cy2 - 0.5 * dyp, -0.5 * dzp)
        pb = kernel.addPoint(0.0, cy2 + 0.5 * dyp, -0.5 * dzp)
        l = kernel.addLine(pa, pb)
        global p3 = first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp)))[2]
    end
    let pa, pb, l
        pa = kernel.addPoint(length_μm, cy1 - 0.5 * dyp, -0.5 * dzp)
        pb = kernel.addPoint(length_μm, cy1 + 0.5 * dyp, -0.5 * dzp)
        l = kernel.addLine(pa, pb)
        global p2 = first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp)))[2]
    end
    let pa, pb, l
        pa = kernel.addPoint(length_μm, cy2 - 0.5 * dyp, -0.5 * dzp)
        pb = kernel.addPoint(length_μm, cy2 + 0.5 * dyp, -0.5 * dzp)
        l = kernel.addLine(pa, pb)
        global p4 = first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp)))[2]
    end
    let pa, pb, l
        pa = kernel.addPoint(0.0, -sep_dy, -sep_dz)
        pb = kernel.addPoint(0.0, dy + sep_dy, -sep_dz)
        l = kernel.addLine(pa, pb)
        global p5 =
            first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, 2.0 * sep_dz)))[2]
    end
    let pa, pb, l
        pa = kernel.addPoint(length_μm, -sep_dy, -sep_dz)
        pb = kernel.addPoint(length_μm, dy + sep_dy, -sep_dz)
        l = kernel.addLine(pa, pb)
        global p6 =
            first(filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, 2.0 * sep_dz)))[2]
    end

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

    metal_group = gmsh.model.addPhysicalGroup(2, metal, -1, "metal")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

    port1 = last.(geom_map[findfirst(x -> x == (2, p1), geom_dimtags)])
    port2 = last.(geom_map[findfirst(x -> x == (2, p2), geom_dimtags)])
    port3 = last.(geom_map[findfirst(x -> x == (2, p3), geom_dimtags)])
    port4 = last.(geom_map[findfirst(x -> x == (2, p4), geom_dimtags)])

    port1_group = gmsh.model.addPhysicalGroup(2, port1, -1, "port1")
    port2_group = gmsh.model.addPhysicalGroup(2, port2, -1, "port2")
    port3_group = gmsh.model.addPhysicalGroup(2, port3, -1, "port3")
    port4_group = gmsh.model.addPhysicalGroup(2, port4, -1, "port4")

    port5 = last.(geom_map[findfirst(x -> x == (2, p5), geom_dimtags)])
    port6 = last.(geom_map[findfirst(x -> x == (2, p6), geom_dimtags)])
    ends = vcat(port5, port6)
    filter!(x -> !(x in port1 || x in port2 || x in port3 || x in port4), ends)

    ends_group = gmsh.model.addPhysicalGroup(2, ends, -1, "ends")

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
    filter!(
        x -> !(x in port1 || x in port2 || x in port3 || x in port4 || x in ends),
        farfield
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
        println("End boundaries: ", ends_group)
        println("Metal boundaries: ", metal_group)
        println("Trace negative boundaries: ", gap_group)

        println("\nPorts:")
        println("Port 1: ", port1_group)
        println("Port 2: ", port2_group)
        println("Port 3: ", port3_group)
        println("Port 4: ", port4_group)
        println()
    end

    # Optionally launch GUI
    if "gui" in lowercase.(ARGS)
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
