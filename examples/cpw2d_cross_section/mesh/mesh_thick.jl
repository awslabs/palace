# CPW 2D cross-section mesh generator with finite metal thickness.
# Generates a Gmsh mesh for a coplanar waveguide cross-section with:
# - Finite-thickness metal trace (PEC hole in the domain)
# - Substrate below, vacuum above
# - Ground planes on both sides
# - Physical groups for Palace boundary conditions and postprocessing
#
# Usage: julia mesh_cpw2d_thick.jl
#
# The mesh is 2D (XY plane), with X along the CPW width and Y vertical.
# The substrate surface is at y = 0. Metal sits on top of the substrate.

import Gmsh: gmsh

function generate_cpw2d_mesh(;
    # Geometry parameters (all in μm)
    w_trace::Float64     = 15.0,          # Trace width
    w_total_cpw::Float64 = 22.0,      # Total width of trace + 2 gaps (constant)
    w_ground::Float64    = 500.0,        # Ground plane width on each side
    t_metal::Float64     = 0.1,           # Metal thickness (must be > 0; use mesh_thin.jl for t=0)
    h_substrate::Float64 = 525.0,     # Substrate thickness
    h_vacuum::Float64    = 500.0,        # Vacuum region height above substrate

    # Mesh parameters
    lc_gap::Float64=0.3,            # Mesh size in the gap region
    lc_trace::Float64=1.0,          # Mesh size on the trace edges
    lc_far::Float64=50.0,           # Mesh size far from the trace
    mesh_order::Int=2,              # Mesh element order (1 or 2)

    # Output
    filename::String="cpw2d_thick.msh",
    verbose::Int=0
)
    # Derived dimensions
    @assert t_metal > 0 "Metal thickness must be > 0. Use mesh_thin.jl for zero-thickness PEC."
    w_gap = (w_total_cpw - w_trace) / 2.0
    @assert w_gap > 0 "Gap width must be positive! Increase w_total_cpw or decrease w_trace."
    w_box = 2.0 * (w_gap + w_ground) + w_trace  # Total simulation box width

    # Center the trace at x = w_box/2
    x_center = w_box / 2.0
    x_trace_left = x_center - w_trace / 2.0
    x_trace_right = x_center + w_trace / 2.0
    x_ground_left_inner = x_trace_left - w_gap    # Left ground inner edge
    x_ground_right_inner = x_trace_right + w_gap   # Right ground inner edge

    # Y coordinates: substrate from -h_substrate to 0, metal from 0 to t_metal
    y_sub_bot = -h_substrate
    y_sub_top = 0.0
    y_metal_bot = 0.0
    y_metal_top = t_metal
    y_vac_top = y_metal_top + h_vacuum

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("cpw2d_thick")

    kernel = gmsh.model.occ

    # Substrate rectangle: full width, below y=0
    sub = kernel.addRectangle(0.0, y_sub_bot, 0.0, w_box, h_substrate)

    # Vacuum rectangle: full width, above metal
    vac = kernel.addRectangle(0.0, y_metal_top, 0.0, w_box, h_vacuum)

    # Left gap fill (between left ground and trace, at metal height)
    gap_left = kernel.addRectangle(x_ground_left_inner, y_metal_bot, 0.0, w_gap, t_metal)

    # Right gap fill (between trace and right ground, at metal height)
    gap_right = kernel.addRectangle(x_trace_right, y_metal_bot, 0.0, w_gap, t_metal)

    # Fragment everything to create shared boundaries
    all_surfs = [(2, sub), (2, vac), (2, gap_left), (2, gap_right)]
    _, frag_map = kernel.fragment(all_surfs, [])
    kernel.synchronize()

    # === Classify entities by bounding box ===
    tol = 1e-6

    # Get all 2D surfaces (domains)
    all_surfaces = gmsh.model.getEntities(2)
    substrate_domains = Int[]
    vacuum_domains = Int[]
    gap_domains = Int[]  # Only for finite metal thickness

    for (dim, tag) in all_surfaces
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        ymid = (ymin + ymax) / 2.0

        if ymid < y_sub_top - tol
            push!(substrate_domains, tag)
        elseif ymin > y_metal_bot - tol && ymax < y_metal_top + tol
            push!(gap_domains, tag)
        else
            push!(vacuum_domains, tag)
        end
    end

    # Get all 1D curves (boundaries)
    all_curves = gmsh.model.getEntities(1)

    pec_curves = Int[]          # PEC boundaries (trace, ground planes)
    outer_curves = Int[]        # Outer box boundary (for PMC or ABC)
    voltage_left_curves = Int[] # Left gap voltage path
    voltage_right_curves = Int[]# Right gap voltage path
    current_curves = Int[]      # Current loop around trace

    for (dim, tag) in all_curves
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        dx = xmax - xmin
        dy = ymax - ymin
        is_horiz = dy < tol
        is_vert = dx < tol

        # Outer box boundaries
        if is_vert && (abs(xmin) < tol || abs(xmax - w_box) < tol)
            push!(outer_curves, tag)
            continue
        end
        if is_horiz && (abs(ymin - y_sub_bot) < tol || abs(ymax - y_vac_top) < tol)
            push!(outer_curves, tag)
            continue
        end

        # Finite metal: PEC is the trace hole boundary + ground surfaces.
        begin

            # Trace top/bottom horizontal edges
            if is_horiz && xmin > x_trace_left - tol && xmax < x_trace_right + tol
                if abs(ymid - y_metal_top) < tol || abs(ymid - y_metal_bot) < tol
                    push!(pec_curves, tag)
                    continue
                end
            end

            # Trace left/right vertical edges
            if is_vert && ymid > y_metal_bot - tol && ymid < y_metal_top + tol
                if abs(xmid - x_trace_left) < tol || abs(xmid - x_trace_right) < tol
                    push!(pec_curves, tag)
                    continue
                end
            end

            # Ground plane PEC: top surface of ground at y = t_metal
            # Left ground: x ∈ [0, x_ground_left_inner]
            if is_horiz && abs(ymid - y_metal_top) < tol
                if xmin < x_ground_left_inner - w_gap + tol &&
                   xmax < x_ground_left_inner + tol
                    push!(pec_curves, tag)
                    continue
                end
                # Right ground: x ∈ [x_ground_right_inner, w_box]
                if xmin > x_ground_right_inner - tol
                    push!(pec_curves, tag)
                    continue
                end
            end

            # Ground plane PEC at y = 0 (bottom of metal region)
            if is_horiz && abs(ymid - y_metal_bot) < tol
                if xmax < x_ground_left_inner + tol || xmin > x_ground_right_inner - tol
                    push!(pec_curves, tag)
                    continue
                end
            end

            # Ground side walls (vertical at inner edges)
            if is_vert && ymid > y_metal_bot - tol && ymid < y_metal_top + tol
                if abs(xmid - x_ground_left_inner) < tol ||
                   abs(xmid - x_ground_right_inner) < tol
                    push!(pec_curves, tag)
                    continue
                end
            end

            # Voltage path: horizontal edges in the gap at y = y_sub_top (substrate surface)
            # Right gap: x ∈ [x_trace_right, x_ground_right_inner] at y = 0
            if is_horiz && abs(ymid - y_sub_top) < tol
                if xmin > x_trace_right - tol && xmax < x_ground_right_inner + tol
                    push!(voltage_right_curves, tag)
                    continue
                end
            end

            # Current loop: domain interfaces forming a rectangle around the trace
            # These are at the gap domain boundaries (not PEC)
            # Left gap left edge (x = x_ground_left_inner)
            if is_vert &&
               abs(xmid - x_ground_left_inner) < tol &&
               ymin > y_metal_bot - tol &&
               ymax < y_metal_top + tol
                # Already classified as PEC above
            end
            # Bottom of gap domains (y = 0)
            if is_horiz && abs(ymid - y_metal_bot) < tol
                if xmin > x_ground_left_inner - tol && xmax < x_trace_left + tol
                    push!(current_curves, tag)
                    continue
                end
                if xmin > x_trace_right - tol && xmax < x_ground_right_inner + tol
                    push!(current_curves, tag)
                    continue
                end
            end
            # Top of gap domains (y = t_metal)
            if is_horiz && abs(ymid - y_metal_top) < tol
                if xmin > x_ground_left_inner - tol && xmax < x_trace_left + tol
                    push!(current_curves, tag)
                    continue
                end
                if xmin > x_trace_right - tol && xmax < x_ground_right_inner + tol
                    push!(current_curves, tag)
                    continue
                end
            end
            # Left side of gap region
            if is_vert &&
               abs(xmid - x_ground_left_inner) < tol &&
               ymin > y_metal_bot - tol &&
               ymax < y_metal_top + tol
                push!(current_curves, tag)
                continue
            end
            # Right side of gap region
            if is_vert &&
               abs(xmid - x_ground_right_inner) < tol &&
               ymin > y_metal_bot - tol &&
               ymax < y_metal_top + tol
                push!(current_curves, tag)
                continue
            end
        end
    end

    # === Physical groups ===
    domain_idx = 1
    bdr_idx = 1

    # Domains
    sub_attr = domain_idx
    gmsh.model.addPhysicalGroup(2, substrate_domains, sub_attr, "substrate")
    domain_idx += 1

    vac_attr = domain_idx
    gmsh.model.addPhysicalGroup(2, vacuum_domains, vac_attr, "vacuum")
    domain_idx += 1

    gap_attr = -1
    if !isempty(gap_domains)
        gap_attr = domain_idx
        gmsh.model.addPhysicalGroup(2, gap_domains, gap_attr, "gap_vacuum")
        domain_idx += 1
    end

    # Boundaries
    pec_attr = bdr_idx
    if !isempty(pec_curves)
        gmsh.model.addPhysicalGroup(1, pec_curves, pec_attr, "pec")
        bdr_idx += 1
    end

    outer_attr = bdr_idx
    if !isempty(outer_curves)
        gmsh.model.addPhysicalGroup(1, outer_curves, outer_attr, "outer")
        bdr_idx += 1
    end

    voltage_right_attr = bdr_idx
    if !isempty(voltage_right_curves)
        gmsh.model.addPhysicalGroup(
            1,
            voltage_right_curves,
            voltage_right_attr,
            "voltage_right"
        )
        bdr_idx += 1
    end

    voltage_left_attr = bdr_idx
    if !isempty(voltage_left_curves)
        gmsh.model.addPhysicalGroup(
            1,
            voltage_left_curves,
            voltage_left_attr,
            "voltage_left"
        )
        bdr_idx += 1
    end

    current_attr = bdr_idx
    if !isempty(current_curves)
        gmsh.model.addPhysicalGroup(1, current_curves, current_attr, "current_loop")
        bdr_idx += 1
    end

    # === Mesh size control ===
    # Distance field from trace edges
    trace_edge_curves = Int[]
    # Trace hole vertical edges (for distance-based mesh refinement)
    for (dim, tag) in all_curves
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        if (abs(xmid - x_trace_left) < tol || abs(xmid - x_trace_right) < tol) &&
           ymid > y_metal_bot - tol &&
           ymid < y_metal_top + tol
            push!(trace_edge_curves, tag)
        end
    end

    # Use all PEC curves for distance field if no specific trace edges found
    if isempty(trace_edge_curves)
        trace_edge_curves = pec_curves
    end

    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(trace_edge_curves))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_gap)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", w_gap)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 10.0 * w_gap)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeMin", lc_gap)
    gmsh.option.setNumber("Mesh.MeshSizeMax", lc_far)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # Generate mesh
    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(mesh_order)

    # Write mesh
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(joinpath(@__DIR__, filename))

    # Print summary
    println("=== CPW 2D Cross-Section Mesh ===")
    println("  Trace width:    $w_trace μm")
    println("  Gap width:      $w_gap μm")
    println("  Ground width:   $w_ground μm")
    println("  Metal thickness: $t_metal μm")
    println("  Box size:       $w_box × $(h_substrate + t_metal + h_vacuum) μm")
    println("  Substrate:      domain attr $sub_attr  (ε_r = 11.47)")
    println("  Vacuum:         domain attr $vac_attr")
    if gap_attr > 0
        println("  Gap vacuum:     domain attr $gap_attr")
    end
    println("  PEC:            bdr attr $pec_attr  ($(length(pec_curves)) curves)")
    println("  Outer box:      bdr attr $outer_attr  ($(length(outer_curves)) curves)")
    if !isempty(voltage_right_curves)
        println(
            "  Voltage (R):    bdr attr $voltage_right_attr  ($(length(voltage_right_curves)) curves)"
        )
    end
    if !isempty(voltage_left_curves)
        println(
            "  Voltage (L):    bdr attr $voltage_left_attr  ($(length(voltage_left_curves)) curves)"
        )
    end
    if !isempty(current_curves)
        println(
            "  Current loop:   bdr attr $current_attr  ($(length(current_curves)) curves)"
        )
    end
    println("  Mesh file:      $filename")

    # Return attribute info for config generation
    attrs = Dict(
        "substrate" => sub_attr,
        "vacuum" => vac_attr,
        "gap" => gap_attr,
        "pec" => pec_attr,
        "outer" => outer_attr,
        "voltage_right" => voltage_right_attr,
        "voltage_left" => !isempty(voltage_left_curves) ? voltage_left_attr : -1,
        "current" => !isempty(current_curves) ? current_attr : -1,
        "x_trace_right" => x_trace_right,
        "x_ground_right" => x_ground_right_inner,
        "y_surface" => y_sub_top
    )

    gmsh.finalize()
    return attrs
end

# Run with default parameters
attrs = generate_cpw2d_mesh()
