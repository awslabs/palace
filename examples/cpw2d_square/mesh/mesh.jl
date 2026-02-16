# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script creates a 2D CPW mesh that closely mirrors the 3D CPW example geometry
as an x-y slice at the metal plane. The lumped ports are square patches (gap_width ×
gap_width) in the gap region at each end of the CPW, matching the 3D port geometry.

The trace is cut from the domain (not meshed), and its boundary becomes PEC. The port
boundary condition is applied to the edge of each square patch closest to the trace.

In 2D, the trace must be cut (not kept as a meshed PEC domain) because a PEC region
fully isolates the two gaps — there is no third dimension for the field to wrap around.

## How to run

```bash
julia -e 'include("mesh.jl"); generate_cpw2d_square_mesh(; filename="cpw2d_square.msh")'
```
=#

using Gmsh: gmsh

"""
    generate_cpw2d_square_mesh(;
        filename::AbstractString,
        trace_width_μm::Real = 30.0,
        gap_width_μm::Real = 18.0,
        ground_width_μm::Real = 800.0,
        cpw_length_μm::Real = 4000.0,
        lc_gap::Real = 3.0,
        lc_far::Real = 60.0,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a 2D CPW mesh with square lumped port patches, mirroring a horizontal slice
of the 3D CPW example. The trace is boolean-cut from the domain so it is not meshed.

Layout (Y coordinates):
y=0          : bottom of domain
y=g          : bottom ground / bottom gap boundary
y=g+s        : bottom gap / trace boundary (PEC, from cut)
y=g+s+w      : trace / top gap boundary (PEC, from cut)
y=g+2s+w     : top gap / top ground boundary
y=2g+2s+w    : top of domain

Port patches are gap_width × gap_width squares at x=0..s and x=L-s..L in each gap.
The port edge is the horizontal side of each square closest to the trace.
"""
function generate_cpw2d_square_mesh(;
    filename::AbstractString,
    trace_width_μm::Real=30.0,
    gap_width_μm::Real=18.0,
    ground_width_μm::Real=800.0,
    cpw_length_μm::Real=4000.0,
    lc_gap::Real=3.0,
    lc_far::Real=60.0,
    verbose::Integer=5,
    gui::Bool=false
)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "cpw2d_square" in gmsh.model.list()
        gmsh.model.setCurrent("cpw2d_square")
        gmsh.model.remove()
    end
    gmsh.model.add("cpw2d_square")

    w = trace_width_μm
    s = gap_width_μm
    g = ground_width_μm
    L = cpw_length_μm
    H = 2g + 2s + w

    # 1. Create the full outer rectangle
    outer = kernel.addRectangle(0.0, 0.0, 0.0, L, H)

    # 2. Create the trace rectangle and CUT it from the domain.
    #    The trace boundary becomes PEC edges in the mesh.
    trace = kernel.addRectangle(0.0, g + s, 0.0, L, w)
    cut_result, _ = kernel.cut([(2, outer)], [(2, trace)])
    # The cut may produce multiple surfaces (top half and bottom half)
    domain_tags_from_cut = [t for (d, t) in cut_result if d == 2]
    @assert length(domain_tags_from_cut) >= 1

    # 3. Create port patch rectangles (s × s squares in each gap at each end)
    # Bottom gap: y from g to g+s
    p1a = kernel.addRectangle(0.0, g, 0.0, s, s)           # Port 1 bottom, left end
    p2a = kernel.addRectangle(L - s, g, 0.0, s, s)         # Port 2 bottom, right end
    # Top gap: y from g+s+w to g+2s+w
    p1b = kernel.addRectangle(0.0, g + s + w, 0.0, s, s)   # Port 1 top, left end
    p2b = kernel.addRectangle(L - s, g + s + w, 0.0, s, s) # Port 2 top, right end

    # 4. Fragment the domain with the port patches to create conformal boundaries
    frag_result, frag_map = kernel.fragment(
        [(2, t) for t in domain_tags_from_cut],
        [(2, p1a), (2, p2a), (2, p1b), (2, p2b)]
    )
    kernel.synchronize()

    # 5. Classify surfaces — everything is domain (trace was cut out)
    all_2d = gmsh.model.getEntities(2)
    domain_tags = [t for (d, t) in all_2d]

    # 6. Classify 1D curves
    all_1d = gmsh.model.getEntities(1)
    tol = 0.5

    trace_bot_y = g + s
    trace_top_y = g + s + w

    trace_pec_curves = Int[]    # Trace boundary from the cut
    outer_pec_curves = Int[]    # Top/bottom/left/right domain boundary
    port1_bot_curves = Int[]
    port1_top_curves = Int[]
    port2_bot_curves = Int[]
    port2_top_curves = Int[]
    other_curves = Int[]

    for (dim, tag) in all_1d
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        dx = xmax - xmin
        dy = ymax - ymin

        is_horizontal = dy < tol && dx > tol
        is_vertical = dx < tol && dy > tol
        is_short = dx < s + tol

        # Trace PEC: horizontal edges at trace_bot_y or trace_top_y
        # Port edges are SHORT horizontal segments at the trace boundary near x=0 or x=L
        if is_horizontal && abs(ymid - trace_bot_y) < tol
            if is_short && xmid < s + tol
                push!(port1_bot_curves, tag)
            elseif is_short && xmid > L - s - tol
                push!(port2_bot_curves, tag)
            else
                push!(trace_pec_curves, tag)
            end
        elseif is_horizontal && abs(ymid - trace_top_y) < tol
            if is_short && xmid < s + tol
                push!(port1_top_curves, tag)
            elseif is_short && xmid > L - s - tol
                push!(port2_top_curves, tag)
            else
                push!(trace_pec_curves, tag)
            end
            # Outer boundary
        elseif is_horizontal && (abs(ymid) < tol || abs(ymid - H) < tol)
            push!(outer_pec_curves, tag)
        elseif is_vertical && (abs(xmid) < tol || abs(xmid - L) < tol)
            push!(outer_pec_curves, tag)
        else
            push!(other_curves, tag)
        end
    end

    # Remaining interior curves (port patch edges, etc.) are left as natural BC
    all_pec = unique(vcat(trace_pec_curves, outer_pec_curves))

    # 7. Create physical groups
    domain_group = gmsh.model.addPhysicalGroup(2, domain_tags, -1, "domain")
    pec_group = gmsh.model.addPhysicalGroup(1, all_pec, -1, "pec")
    p1b_group = gmsh.model.addPhysicalGroup(1, unique(port1_bot_curves), -1, "port1_bot")
    p1t_group = gmsh.model.addPhysicalGroup(1, unique(port1_top_curves), -1, "port1_top")
    p2b_group = gmsh.model.addPhysicalGroup(1, unique(port2_bot_curves), -1, "port2_bot")
    p2t_group = gmsh.model.addPhysicalGroup(1, unique(port2_top_curves), -1, "port2_top")

    # 8. Mesh sizing
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc_gap / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    # Refine near trace edges
    trace_ref_curves = unique(
        vcat(
            trace_pec_curves,
            port1_bot_curves,
            port1_top_curves,
            port2_bot_curves,
            port2_top_curves
        )
    )
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(trace_ref_curves))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_gap)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", s)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 10.0 * s)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)

    # 9. Write
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating CPW 2D (square ports) mesh.")
    println("Physical group tags:")
    println("  Domain: ", domain_group)
    println("  PEC (trace + outer): ", pec_group)
    println("  Port 1 bottom: ", p1b_group)
    println("  Port 1 top: ", p1t_group)
    println("  Port 2 bottom: ", p2b_group)
    println("  Port 2 top: ", p2t_group)
    println()

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
