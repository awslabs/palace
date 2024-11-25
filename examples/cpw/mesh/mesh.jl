# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

using Gmsh: gmsh

"""
    function generate_cpw_wave_mesh(;
        filename::AbstractString,
        refinement::Integer       = 0,
        order::Integer            = 1,
        trace_width_μm::Real      = 30.0,
        gap_width_μm::Real        = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real     = 800.0,
        substrate_height_μm::Real = 500.0,
        metal_height_μm::Real     = 0.0,
        remove_metal_vol::Bool    = true,
        length_μm::Real           = 4000.0,
        coax_ports::Bool          = false,
        verbose::Integer          = 5,
        gui::Bool                 = false
    )

Generate a mesh for the coplanar waveguide with wave ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - coax_ports - flag to use coaxial lumped ports instead of wave ports
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cpw_wave_mesh(;
    filename::AbstractString,
    refinement::Integer       = 0,
    order::Integer            = 1,
    trace_width_μm::Real      = 30.0,
    gap_width_μm::Real        = 18.0,
    separation_width_μm::Real = 200.0,
    ground_width_μm::Real     = 800.0,
    substrate_height_μm::Real = 500.0,
    metal_height_μm::Real     = 0.0,
    remove_metal_vol::Bool    = true,
    length_μm::Real           = 4000.0,
    coax_ports::Bool          = false,
    verbose::Integer          = 5,
    gui::Bool                 = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert separation_width_μm > 0
    @assert ground_width_μm > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
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
    l_trace = 1.5 * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * substrate_height_μm * (2.0^-refinement)

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

    # Metal thickness
    metal_boundary = [g1, g2, g3, t1, t2]
    metal_boundary_top = typeof(metal_boundary)(undef, 0)
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
            append!(metal_boundary_top, first(boundary))
        end
        filter!(x -> !(x in metal_boundary), metal_boundary_top)
    end

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
    if coax_ports
        ra = 0.5 * trace_width_μm
        rb = 0.5 * trace_width_μm + gap_width_μm
        let da, db
            da = kernel.addDisk(0.0, cy1, 0.0, ra, ra, -1, [1, 0, 0], [])
            db = kernel.addDisk(0.0, cy1, 0.0, rb, rb, -1, [1, 0, 0], [])
            global p1 = last(first(first(kernel.cut((2, db), (2, da)))))
        end
        let da, db
            da = kernel.addDisk(0.0, cy2, 0.0, ra, ra, -1, [1, 0, 0], [])
            db = kernel.addDisk(0.0, cy2, 0.0, rb, rb, -1, [1, 0, 0], [])
            global p3 = last(first(first(kernel.cut((2, db), (2, da)))))
        end
        let da, db
            da = kernel.addDisk(length_μm, cy1, 0.0, ra, ra, -1, [1, 0, 0], [])
            db = kernel.addDisk(length_μm, cy1, 0.0, rb, rb, -1, [1, 0, 0], [])
            global p2 = last(first(first(kernel.cut((2, db), (2, da)))))
        end
        let da, db
            da = kernel.addDisk(length_μm, cy2, 0.0, ra, ra, -1, [1, 0, 0], [])
            db = kernel.addDisk(length_μm, cy2, 0.0, rb, rb, -1, [1, 0, 0], [])
            global p4 = last(first(first(kernel.cut((2, db), (2, da)))))
        end
    else
        dyp1 = 0.5 * (cy2 - cy1)
        dyp2 = dyp1
        dzp1 = 0.5 * (sep_dz + substrate_height_μm)
        dzp2 = substrate_height_μm
        let pa, pb, l
            pa = kernel.addPoint(0.0, cy1 - dyp2, -dzp1)
            pb = kernel.addPoint(0.0, cy1 + dyp1, -dzp1)
            l = kernel.addLine(pa, pb)
            global p1 = first(
                filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp1 + dzp2))
            )[2]
        end
        let pa, pb, l
            pa = kernel.addPoint(0.0, cy2 - dyp1, -dzp1)
            pb = kernel.addPoint(0.0, cy2 + dyp2, -dzp1)
            l = kernel.addLine(pa, pb)
            global p3 = first(
                filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp1 + dzp2))
            )[2]
        end
        let pa, pb, l
            pa = kernel.addPoint(length_μm, cy1 - dyp2, -dzp1)
            pb = kernel.addPoint(length_μm, cy1 + dyp1, -dzp1)
            l = kernel.addLine(pa, pb)
            global p2 = first(
                filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp1 + dzp2))
            )[2]
        end
        let pa, pb, l
            pa = kernel.addPoint(length_μm, cy2 - dyp1, -dzp1)
            pb = kernel.addPoint(length_μm, cy2 + dyp2, -dzp1)
            l = kernel.addLine(pa, pb)
            global p4 = first(
                filter(x -> x[1] == 2, kernel.extrude([1, l], 0.0, 0.0, dzp1 + dzp2))
            )[2]
        end
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
    metal_domains =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
                )
            )
        )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in
            last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "si")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

    port1 = last.(geom_map[findfirst(x -> x == (2, p1), geom_dimtags)])
    port2 = last.(geom_map[findfirst(x -> x == (2, p2), geom_dimtags)])
    port3 = last.(geom_map[findfirst(x -> x == (2, p3), geom_dimtags)])
    port4 = last.(geom_map[findfirst(x -> x == (2, p4), geom_dimtags)])

    port1_group = gmsh.model.addPhysicalGroup(2, port1, -1, "port1")
    port2_group = gmsh.model.addPhysicalGroup(2, port2, -1, "port2")
    port3_group = gmsh.model.addPhysicalGroup(2, port3, -1, "port3")
    port4_group = gmsh.model.addPhysicalGroup(2, port4, -1, "port4")

    end1 = last.(geom_map[findfirst(x -> x == (2, p5), geom_dimtags)])
    end2 = last.(geom_map[findfirst(x -> x == (2, p6), geom_dimtags)])
    filter!(x -> !(x in port1 || x in port3), end1)
    filter!(x -> !(x in port2 || x in port4), end2)

    end1_group = gmsh.model.addPhysicalGroup(2, end1, -1, "end1")
    end2_group = gmsh.model.addPhysicalGroup(2, end2, -1, "end2")

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
        x -> !(
            x in port1 || x in port2 || x in port3 || x in port4 || x in end1 || x in end2
        ),
        farfield
    )

    farfield_group = gmsh.model.addPhysicalGroup(2, farfield, -1, "farfield")

    trace =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x[1] == 2 && x[2] in metal_boundary,
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
                        x -> x[1] == 2 && x[2] in [n1, n2, n3, n4],
                        geom_dimtags
                    )]
                )
            )
        )

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

    trace_top =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x[1] == 2 && x[2] in metal_boundary_top,
                        geom_dimtags
                    )]
                )
            )
        )
    filter!(
        x -> !(
            x in port1 || x in port2 || x in port3 || x in port4 || x in end1 || x in end2
        ),
        trace_top
    )

    trace_top_group = gmsh.model.addPhysicalGroup(2, trace_top, -1, "trace2")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gap_points =
        last.(
            filter(
                x -> x[1] == 0,
                gmsh.model.getBoundary([(2, z) for z in gap], false, true, true)
            )
        )
    gap_curves =
        last.(
            filter(
                x -> x[1] == 1,
                gmsh.model.getBoundary([(2, z) for z in gap], false, false, false)
            )
        )

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", gap_points)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", gap)
    gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.9 * sep_dz)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    for tag in Iterators.flatten((gap, trace, trace_top))
        gmsh.model.mesh.setAlgorithm(2, tag, 8)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Near-end boundaries: ", end1_group)
        println("Far-end boundaries: ", end2_group)
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nPorts:")
        println("Port 1: ", port1_group)
        println("Port 2: ", port2_group)
        println("Port 3: ", port3_group)
        println("Port 4: ", port4_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end

"""
    function generate_cpw_lumped_mesh(;
        filename::AbstractString,
        refinement::Integer       = 0,
        order::Integer            = 1,
        trace_width_μm::Real      = 30.0,
        gap_width_μm::Real        = 18.0,
        separation_width_μm::Real = 200.0,
        ground_width_μm::Real     = 800.0,
        substrate_height_μm::Real = 500.0,
        metal_height_μm::Real     = 0.0,
        remove_metal_vol::Bool    = true,
        length_μm::Real           = 4000.0,
        verbose::Integer          = 5,
        gui::Bool                 = false
    )

Generate a mesh for the coplanar waveguide with lumped ports using Gmsh

# Arguments

  - filename - the filename to use for the generated mesh
  - refinement - measure of how many elements to include, 0 is least
  - order - the polynomial order of the approximation, minimum 1
  - trace_width_μm - width of the coplanar waveguide trace, in μm
  - gap_width_μm - width of the coplanar waveguide gap, in μm
  - separation_width_μm - separation distance between the two waveguides, in μm
  - ground_width_μm - width of the ground plane, in μm
  - substrate_height_μm - height of the substrate, in μm
  - metal_height_μm - metal thickness, in μm
  - remove_metal_vol - for positive metal thickness, whether to remove the metal domain
  - length_μm - length of the waveguides, in μm
  - verbose - flag to dictate the level of print to REPL, passed to Gmsh
  - gui - whether to launch the Gmsh GUI on mesh generation
"""
function generate_cpw_lumped_mesh(;
    filename::AbstractString,
    refinement::Integer       = 0,
    order::Integer            = 1,
    trace_width_μm::Real      = 30.0,
    gap_width_μm::Real        = 18.0,
    separation_width_μm::Real = 200.0,
    ground_width_μm::Real     = 800.0,
    substrate_height_μm::Real = 500.0,
    metal_height_μm::Real     = 0.0,
    remove_metal_vol::Bool    = true,
    length_μm::Real           = 4000.0,
    verbose::Integer          = 5,
    gui::Bool                 = false
)
    @assert refinement >= 0
    @assert order > 0
    @assert trace_width_μm > 0
    @assert gap_width_μm > 0
    @assert separation_width_μm > 0
    @assert ground_width_μm > 0
    @assert substrate_height_μm > 0
    @assert metal_height_μm >= 0
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
    l_trace = 1.5 * trace_width_μm * (2.0^-refinement)
    l_farfield = 1.0 * substrate_height_μm * (2.0^-refinement)

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

    # Metal thickness
    metal_boundary = [g1, g2, g3, t1, t2]
    metal_boundary_top = typeof(metal_boundary)(undef, 0)
    metal = typeof(metal_boundary)(undef, 0)
    if metal_height_μm > 0
        metal_dimtags =
            kernel.extrude([(2, x) for x in metal_boundary], 0.0, 0.0, metal_height_μm)
        metal = [x[2] for x in filter(x -> x[1] == 3, metal_dimtags)]
        for domain in metal
            _, boundary = kernel.getSurfaceLoops(domain)
            @assert length(boundary) == 1
            append!(metal_boundary_top, first(boundary))
        end
        filter!(x -> !(x in metal_boundary), metal_boundary_top)
    end

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
    metal_domains =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(x -> x[1] == 3 && x[2] in metal, geom_dimtags)]
                )
            )
        )

    si_domain = last.(geom_map[findfirst(x -> x == (3, substrate), geom_dimtags)])
    @assert length(si_domain) == 1
    si_domain = first(si_domain)

    air_domain = last.(geom_map[findfirst(x -> x == (3, domain), geom_dimtags)])
    filter!(x -> !(x == si_domain || x in metal_domains), air_domain)
    @assert length(air_domain) == 1
    air_domain = first(air_domain)

    if length(metal_domains) > 0 && remove_metal_vol
        remove_dimtags = [(3, x) for x in metal_domains]
        for tag in
            last.(
            filter(
                x -> x[1] == 2,
                gmsh.model.getBoundary(
                    [(3, z) for z in metal_domains],
                    false,
                    false,
                    false
                )
            )
        )
            normal = gmsh.model.getNormal(tag, [0, 0])
            if abs(normal[1]) == 1.0
                push!(remove_dimtags, (2, tag))
            end
        end
        kernel.remove(remove_dimtags)
        kernel.synchronize()
        filter!.(x -> !(x in remove_dimtags), geom_map)
        empty!(metal_domains)
    end

    air_domain_group = gmsh.model.addPhysicalGroup(3, [air_domain], -1, "air")
    si_domain_group = gmsh.model.addPhysicalGroup(3, [si_domain], -1, "si")
    metal_domain_group = gmsh.model.addPhysicalGroup(3, metal_domains, -1, "metal")

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

    trace =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x[1] == 2 && x[2] in metal_boundary,
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
                        x -> x[1] == 2 && x[2] in [n1, n2, n3, n4],
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

    trace_group = gmsh.model.addPhysicalGroup(2, trace, -1, "trace")
    gap_group = gmsh.model.addPhysicalGroup(2, gap, -1, "gap")

    trace_top =
        last.(
            collect(
                Iterators.flatten(
                    geom_map[findall(
                        x -> x[1] == 2 && x[2] in metal_boundary_top,
                        geom_dimtags
                    )]
                )
            )
        )

    trace_top_group = gmsh.model.addPhysicalGroup(2, trace_top, -1, "trace2")

    # Generate mesh
    gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gap_points =
        last.(
            filter(
                x -> x[1] == 0,
                gmsh.model.getBoundary([(2, z) for z in gap], false, true, true)
            )
        )
    gap_curves =
        last.(
            filter(
                x -> x[1] == 1,
                gmsh.model.getBoundary([(2, z) for z in gap], false, false, false)
            )
        )

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", gap_points)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", gap_curves)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", gap)
    gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(length_μm / l_trace))

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(2, "DistMin", trace_width_μm)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.9 * sep_dz)

    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    for tag in Iterators.flatten((gap, trace, trace_top))
        gmsh.model.mesh.setAlgorithm(2, tag, 8)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)

    # Save mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print some information
    if verbose > 0
        println("\nFinished generating mesh. Physical group tags:")
        println("Air domain: ", air_domain_group)
        println("Si domain: ", si_domain_group)
        if length(metal_domains) > 0
            println("Metal domain: ", metal_domain_group)
        end
        println("Farfield boundaries: ", farfield_group)
        println("Metal boundaries: ", trace_group)
        if length(trace_top) > 0
            println("Extruded metal boundaries: ", trace_top_group)
        end
        println("Negative trace boundaries: ", gap_group)

        println("\nMultielement lumped ports:")
        println("Port 1: ", port1a_group, ", ", port1b_group)
        println("Port 2: ", port2a_group, ", ", port2b_group)
        println("Port 3: ", port3a_group, ", ", port3b_group)
        println("Port 4: ", port4a_group, ", ", port4b_group)
        println()
    end

    # Optionally launch GUI
    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
