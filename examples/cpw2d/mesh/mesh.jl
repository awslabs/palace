# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses Gmsh to create a 2D mesh for a coplanar waveguide (CPW) in top-down
view, with two square lumped ports at each end.

The geometry consists of a center trace (PEC) running the length of the domain, flanked by
two gaps and ground planes (PEC boundaries at the top and bottom). The trace does not extend
to the domain ends, leaving space for the lumped ports.

Each port has two elements connecting the trace to the top and bottom ground planes.

## Prerequisites

```bash
julia -e 'using Pkg; Pkg.add(Pkg.PackageSpec(name="Gmsh", uuid="705231aa-382f-11e9-3f0c-b7cb4346fdeb"))'
```

## How to run

```bash
julia -e 'include("mesh.jl"); generate_cpw2d_mesh(; filename="cpw2d.msh")'
```
=#

using Gmsh: gmsh

"""
    generate_cpw2d_mesh(;
        filename::AbstractString,
        trace_w::Real = 30.0,
        gap_s::Real = 18.0,
        cpw_length::Real = 4000.0,
        port_inset::Real = 30.0,
        lc_gap::Real = 3.0,
        lc_far::Real = 30.0,
        verbose::Integer = 5,
        gui::Bool = false
    )

Generate a 2D CPW mesh with lumped ports.

# Arguments
  - trace_w - center trace width (μm)
  - gap_s - gap between trace and ground (μm)
  - cpw_length - total CPW length (μm)
  - port_inset - distance from domain edge to port / trace end (μm)
  - lc_gap - mesh size in the gap region (μm)
  - lc_far - mesh size far from the gap (μm)
"""
function generate_cpw2d_mesh(;
    filename::AbstractString,
    trace_w::Real=30.0,
    gap_s::Real=18.0,
    cpw_length::Real=4000.0,
    port_inset::Real=30.0,
    lc_gap::Real=3.0,
    lc_far::Real=30.0,
    verbose::Integer=5,
    gui::Bool=false
)
    kernel = gmsh.model.occ
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)

    if "cpw2d" in gmsh.model.list()
        gmsh.model.setCurrent("cpw2d")
        gmsh.model.remove()
    end
    gmsh.model.add("cpw2d")

    # Derived dimensions
    half_w = trace_w / 2.0
    half_h = half_w + gap_s  # half-height of domain

    # Create the outer rectangle (full domain)
    outer = kernel.addRectangle(0.0, -half_h, 0.0, cpw_length, 2.0 * half_h)

    # Create the center trace rectangle (will be cut from domain)
    trace = kernel.addRectangle(port_inset, -half_w, 0.0,
                                 cpw_length - 2.0 * port_inset, trace_w)

    # Cut the trace from the domain
    cut_result, _ = kernel.cut([(2, outer)], [(2, trace)])
    @assert length(cut_result) == 1 && cut_result[1][1] == 2
    domain_tag = cut_result[1][2]

    # Now we need to fragment the domain to create port edges.
    # Create small port rectangles that will split the gap at the port locations.
    # Port 1 (left): vertical line at x = port_inset in both gaps
    # Port 2 (right): vertical line at x = cpw_length - port_inset in both gaps
    # We use thin rectangles as "blades" to fragment the mesh.
    eps = 0.1  # thin blade width
    blade1 = kernel.addRectangle(port_inset - eps / 2, -half_h, 0.0, eps, gap_s)
    blade2 = kernel.addRectangle(port_inset - eps / 2, half_w, 0.0, eps, gap_s)
    blade3 = kernel.addRectangle(cpw_length - port_inset - eps / 2, -half_h, 0.0, eps, gap_s)
    blade4 = kernel.addRectangle(cpw_length - port_inset - eps / 2, half_w, 0.0, eps, gap_s)

    # Fragment the domain with the blades to create shared edges at the port locations
    frag_result, frag_map = kernel.fragment(
        [(2, domain_tag)],
        [(2, blade1), (2, blade2), (2, blade3), (2, blade4)]
    )

    kernel.synchronize()

    # Get all 2D surfaces (these form the domain)
    all_surfaces = gmsh.model.getEntities(2)
    domain_tags = [t for (d, t) in all_surfaces]

    # Get all 1D curves (boundaries)
    all_curves = gmsh.model.getEntities(1)

    # Classify curves by location
    tol = 1.0
    trace_curves = Int[]
    ground_curves = Int[]
    end_curves = Int[]
    port1_bot_curves = Int[]
    port1_top_curves = Int[]
    port2_bot_curves = Int[]
    port2_top_curves = Int[]
    other_curves = Int[]

    for (dim, tag) in all_curves
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        dx = xmax - xmin
        dy = ymax - ymin

        is_vertical = dx < tol && dy > tol
        is_horizontal = dy < tol && dx > tol

        if is_horizontal && abs(ymid - half_w) < tol && xmid > port_inset - tol &&
           xmid < cpw_length - port_inset + tol
            # Top edge of trace
            push!(trace_curves, tag)
        elseif is_horizontal && abs(ymid + half_w) < tol && xmid > port_inset - tol &&
               xmid < cpw_length - port_inset + tol
            # Bottom edge of trace
            push!(trace_curves, tag)
        elseif is_vertical && abs(xmid - port_inset) < tol &&
               ymid > port_inset - tol && abs(ymid) > half_w - tol
            # Vertical trace end at left (part of trace boundary)
            if abs(ymid) < half_w + tol
                push!(trace_curves, tag)
            end
        elseif is_vertical && abs(xmid - (cpw_length - port_inset)) < tol &&
               abs(ymid) < half_w + tol
            # Vertical trace end at right
            push!(trace_curves, tag)
        elseif is_horizontal && abs(ymid - half_h) < tol
            # Top boundary (ground)
            push!(ground_curves, tag)
        elseif is_horizontal && abs(ymid + half_h) < tol
            # Bottom boundary (ground)
            push!(ground_curves, tag)
        elseif is_vertical && abs(xmid) < tol
            # Left end
            push!(end_curves, tag)
        elseif is_vertical && abs(xmid - cpw_length) < tol
            # Right end
            push!(end_curves, tag)
        elseif is_vertical && abs(xmid - port_inset) < tol && ymid < -half_w + tol
            # Port 1 bottom gap
            push!(port1_bot_curves, tag)
        elseif is_vertical && abs(xmid - port_inset) < tol && ymid > half_w - tol
            # Port 1 top gap
            push!(port1_top_curves, tag)
        elseif is_vertical && abs(xmid - (cpw_length - port_inset)) < tol &&
               ymid < -half_w + tol
            # Port 2 bottom gap
            push!(port2_bot_curves, tag)
        elseif is_vertical && abs(xmid - (cpw_length - port_inset)) < tol &&
               ymid > half_w - tol
            # Port 2 top gap
            push!(port2_top_curves, tag)
        else
            push!(other_curves, tag)
        end
    end

    # Create physical groups
    domain_group = gmsh.model.addPhysicalGroup(2, domain_tags, -1, "domain")

    pec_trace_curves = unique(trace_curves)
    pec_trace_group = -1
    if !isempty(pec_trace_curves)
        pec_trace_group = gmsh.model.addPhysicalGroup(1, pec_trace_curves, -1, "pec_trace")
    end

    pec_ground_group = gmsh.model.addPhysicalGroup(1, unique(ground_curves), -1, "pec_ground")
    pec_ends_group = gmsh.model.addPhysicalGroup(1, unique(end_curves), -1, "pec_ends")

    port1_bot_group = gmsh.model.addPhysicalGroup(1, unique(port1_bot_curves), -1, "port1_bot")
    port1_top_group = gmsh.model.addPhysicalGroup(1, unique(port1_top_curves), -1, "port1_top")
    port2_bot_group = gmsh.model.addPhysicalGroup(1, unique(port2_bot_curves), -1, "port2_bot")
    port2_top_group = gmsh.model.addPhysicalGroup(1, unique(port2_top_curves), -1, "port2_top")

    # Assign any remaining boundary curves to PEC (catch-all for blade edges, etc.)
    remaining = setdiff(
        Set([t for (_, t) in all_curves]),
        Set(vcat(trace_curves, ground_curves, end_curves,
                 port1_bot_curves, port1_top_curves,
                 port2_bot_curves, port2_top_curves))
    )
    if !isempty(remaining)
        # These are internal blade fragment edges — add them to trace PEC
        if pec_trace_group > 0
            gmsh.model.removePhysicalGroups([(1, pec_trace_group)])
        end
        pec_trace_group = gmsh.model.addPhysicalGroup(
            1, unique(vcat(trace_curves, collect(remaining))), -1, "pec_trace")
    end

    # Mesh size control: fine in the gap, coarser away
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc_gap / 2.0)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    # Add mesh size field to refine near the trace
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(pec_trace_curves))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_gap)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", gap_s)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 5.0 * gap_s)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(2)

    # Write mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    println("\nFinished generating CPW 2D mesh. Physical group tags:")
    println("Domain: ", domain_group)
    println("PEC trace: ", pec_trace_group)
    println("PEC ground: ", pec_ground_group)
    println("PEC ends: ", pec_ends_group)
    println("Port 1 bottom: ", port1_bot_group)
    println("Port 1 top: ", port1_top_group)
    println("Port 2 bottom: ", port2_bot_group)
    println("Port 2 top: ", port2_top_group)
    println()

    if gui
        gmsh.fltk.run()
    end

    return gmsh.finalize()
end
