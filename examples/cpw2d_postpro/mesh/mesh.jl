# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script creates a 2D CPW cross-section mesh for electrostatic analysis.
It uses Gmsh for geometry and meshing, then writes MFEM native format with explicit
internal boundary elements at the metal plane (z=0).

MFEM's Gmsh reader does not import 1D physical groups on internal edges for 2D meshes,
so we write the .mesh file directly to include internal boundary elements for the metal
traces and substrate-air interface.

## How to run

```bash
julia -e 'include("mesh.jl"); generate_cpw_cross_section(; filename="cpw2d_postpro.mesh")'
```
=#

using Gmsh: gmsh

"""
    generate_cpw_cross_section(;
        filename, trace_width_μm=30.0, gap_width_μm=18.0, ground_width_μm=800.0,
        substrate_height_μm=500.0, air_height_μm=500.0,
        lc_gap=1.5, lc_far=40.0, verbose=5, gui=false)

Generate a 2D CPW cross-section mesh in MFEM format. Output attributes:
Domain 1: air (z > 0)
Domain 2: substrate (z < 0)
Boundary 3: ground_left (z=0, y < g)
Boundary 4: ground_right (z=0, y > g+2s+w)
Boundary 5: trace (z=0, g+s < y < g+s+w)
Boundary 6: outer boundary
Boundary 7: SA interface (z=0 gaps, g < y < g+s and g+s+w < y < g+2s+w)
"""
function generate_cpw_cross_section(;
    filename::AbstractString,
    trace_width_μm::Real=30.0,
    gap_width_μm::Real=18.0,
    ground_width_μm::Real=800.0,
    substrate_height_μm::Real=500.0,
    air_height_μm::Real=500.0,
    lc_gap::Real=1.5,
    lc_far::Real=40.0,
    verbose::Integer=5,
    gui::Bool=false
)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "cpw_xsec" in gmsh.model.list()
        gmsh.model.setCurrent("cpw_xsec")
        gmsh.model.remove()
    end
    gmsh.model.add("cpw_xsec")

    w = trace_width_μm
    s = gap_width_μm
    g = ground_width_μm
    h_sub = substrate_height_μm
    h_air = air_height_μm
    W = 2g + 2s + w

    trace_y0 = g + s
    trace_y1 = g + s + w

    # Two rectangles: air and substrate
    air = kernel.addRectangle(0.0, 0.0, 0.0, W, h_air)
    sub = kernel.addRectangle(0.0, -h_sub, 0.0, W, h_sub)

    # Fragment for conformal mesh at z=0
    frag_result, _ = kernel.fragment([(2, air)], [(2, sub)])
    kernel.synchronize()

    # Mesh size control
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc_gap / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    # Refine near trace region
    # Find edges at z=0 near the trace
    all_1d = gmsh.model.getEntities(1)
    trace_edge_tags = Int[]
    for (dim, tag) in all_1d
        bb = gmsh.model.getBoundingBox(dim, tag)
        zmid = (bb[3] + bb[6]) / 2.0
        ymid = (bb[2] + bb[5]) / 2.0
        if abs(zmid) < 0.5 && ymid > g - s && ymid < g + 2s + w + s
            push!(trace_edge_tags, tag)
        end
    end

    if !isempty(trace_edge_tags)
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(trace_edge_tags))
        gmsh.model.mesh.field.setNumber(1, "Sampling", 200)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_gap)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
        gmsh.model.mesh.field.setNumber(2, "DistMin", s)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 10.0 * s)

        gmsh.model.mesh.field.setAsBackgroundMesh(2)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    end

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(1)  # Linear mesh for MFEM native format

    # Extract mesh data from Gmsh
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    # Build node index map (Gmsh tags are 1-based but may have gaps)
    node_map = Dict{Int, Int}()
    for (i, tag) in enumerate(node_tags)
        node_map[tag] = i
    end
    nv = length(node_tags)

    # Get 2D elements (triangles)
    all_2d_entities = gmsh.model.getEntities(2)
    elements = []
    for (dim, tag) in all_2d_entities
        bb = gmsh.model.getBoundingBox(dim, tag)
        zmid = (bb[3] + bb[6]) / 2.0
        attr = zmid >= -0.5 ? 1 : 2  # 1=air, 2=substrate

        elem_types, elem_tags_list, elem_node_tags = gmsh.model.mesh.getElements(dim, tag)
        for (etype, etags, enodes) in zip(elem_types, elem_tags_list, elem_node_tags)
            if etype == 2  # Triangle
                nnpe = 3
                for i = 1:length(etags)
                    n1 = node_map[enodes[(i - 1) * nnpe + 1]]
                    n2 = node_map[enodes[(i - 1) * nnpe + 2]]
                    n3 = node_map[enodes[(i - 1) * nnpe + 3]]
                    push!(elements, (attr, n1, n2, n3))
                end
            end
        end
    end
    ne = length(elements)

    # Build edge-to-element connectivity to find boundary and internal edges
    edge_elements = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (idx, (attr, n1, n2, n3)) in enumerate(elements)
        for (a, b) in [(n1, n2), (n2, n3), (n3, n1)]
            edge = (min(a, b), max(a, b))
            if !haskey(edge_elements, edge)
                edge_elements[edge] = Int[]
            end
            push!(edge_elements[edge], idx)
        end
    end

    # Classify edges
    boundary_elements = []
    tol = 0.5

    for ((n1, n2), elems) in edge_elements
        x1, y1 = node_coords[(node_map == nothing ? 0 : 0) + 1], 0.0  # placeholder
        # Get actual coordinates
        idx1 = findfirst(==(n1), [node_map[t] for t in node_tags])
        idx2 = findfirst(==(n2), [node_map[t] for t in node_tags])
        y1_coord = node_coords[3 * (idx1 - 1) + 1]  # x in Gmsh (= y in our convention)
        z1_coord = node_coords[3 * (idx1 - 1) + 2]  # y in Gmsh (= z in our convention)
        y2_coord = node_coords[3 * (idx2 - 1) + 1]
        z2_coord = node_coords[3 * (idx2 - 1) + 2]

        ymid = (y1_coord + y2_coord) / 2.0
        zmid = (z1_coord + z2_coord) / 2.0

        is_external = length(elems) == 1
        is_internal_interface =
            length(elems) == 2 && elements[elems[1]][1] != elements[elems[2]][1]  # Different domain attrs

        if is_external
            # External boundary
            if abs(zmid + h_sub) < tol ||
               abs(zmid - h_air) < tol ||
               abs(ymid) < tol ||
               abs(ymid - W) < tol
                push!(boundary_elements, (6, n1, n2))  # Outer boundary
            end
        elseif is_internal_interface && abs(zmid) < tol
            # Internal interface at z=0 between air and substrate
            if ymid < g + tol
                push!(boundary_elements, (3, n1, n2))  # Ground left
            elseif ymid < trace_y0 + tol
                push!(boundary_elements, (7, n1, n2))  # SA interface (left gap)
            elseif ymid < trace_y1 + tol
                push!(boundary_elements, (5, n1, n2))  # Trace
            elseif ymid < g + 2s + w + tol
                push!(boundary_elements, (7, n1, n2))  # SA interface (right gap)
            else
                push!(boundary_elements, (4, n1, n2))  # Ground right
            end
        end
    end
    nbe = length(boundary_elements)

    # Write MFEM mesh
    open(joinpath(@__DIR__, filename), "w") do f
        write(f, "MFEM mesh v1.0\n\n")
        write(f, "dimension\n2\n\n")

        write(f, "elements\n$ne\n")
        for (attr, n1, n2, n3) in elements
            write(f, "$attr 2 $(n1-1) $(n2-1) $(n3-1)\n")  # 0-indexed, type 2=triangle
        end

        write(f, "\nboundary\n$nbe\n")
        for (attr, n1, n2) in boundary_elements
            write(f, "$attr 1 $(n1-1) $(n2-1)\n")  # 0-indexed, type 1=segment
        end

        write(f, "\nvertices\n$nv\n2\n")
        for i = 1:nv
            # Find the original tag for this index
            orig_tag = node_tags[i]
            x = node_coords[3 * (i - 1) + 1]
            y = node_coords[3 * (i - 1) + 2]
            write(f, "$x $y\n")
        end
    end

    println("\nFinished generating CPW cross-section mesh.")
    println("  Vertices: $nv, Elements: $ne, Boundary elements: $nbe")
    println("  Domain attrs: 1=air, 2=substrate")
    println("  Boundary attrs: 3=ground_left, 4=ground_right, 5=trace, 6=outer, 7=SA")
    println()

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
