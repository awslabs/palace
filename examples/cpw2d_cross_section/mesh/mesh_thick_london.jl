# CPW 2D cross-section mesh with meshed metal regions for London equation.
# Unlike mesh_thick.jl (which treats metals as PEC holes), this script meshes
# the metal trace and ground planes as separate domains. The London equation
# is applied inside the metal domains via Palace's material properties.
#
# Usage: julia mesh_thick_london.jl
#
# The mesh is 2D (XY plane), with X along the CPW width and Y vertical.
# The substrate surface is at y = 0. Metal sits on top of the substrate.

import Gmsh: gmsh

function generate_cpw2d_london_mesh(;
    # Geometry parameters (all in μm)
    w_trace::Float64     = 15.0,         # Trace width
    w_total_cpw::Float64 = 22.0,         # Total width of trace + 2 gaps
    w_ground::Float64    = 100.0,        # Ground plane width on each side
    t_metal::Float64     = 0.1,          # Metal thickness
    h_substrate::Float64 = 525.0,        # Substrate thickness
    h_vacuum::Float64    = 500.0,        # Vacuum region height

    # Mesh parameters — finer near the gap to resolve the thin metal layer.
    lc_gap::Float64=0.1,            # Mesh size near the gap (~metal thickness, gives ~2 layers)
    lc_trace::Float64=1.0,          # Mesh size on the trace edges
    lc_far::Float64=50.0,           # Mesh size far from the trace
    mesh_order::Int=2,              # Mesh element order (1 or 2)

    # Output
    filename::String="cpw2d_thick_london.msh",
    verbose::Int=0
)
    # Derived dimensions
    @assert t_metal > 0 "Metal thickness must be > 0."
    w_gap = (w_total_cpw - w_trace) / 2.0
    @assert w_gap > 0 "Gap width must be positive!"
    w_box = 2.0 * (w_gap + w_ground) + w_trace  # Total simulation box width

    # Center the trace at x = w_box/2
    x_center = w_box / 2.0
    x_trace_left = x_center - w_trace / 2.0
    x_trace_right = x_center + w_trace / 2.0
    x_ground_left_inner = x_trace_left - w_gap
    x_ground_right_inner = x_trace_right + w_gap

    # Y coordinates
    y_sub_bot = -h_substrate
    y_sub_top = 0.0
    y_metal_bot = 0.0
    y_metal_top = t_metal
    y_vac_top = y_metal_top + h_vacuum

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("cpw2d_thick_london")

    kernel = gmsh.model.occ

    # Substrate: full width, below y=0
    sub = kernel.addRectangle(0.0, y_sub_bot, 0.0, w_box, h_substrate)

    # Vacuum: full width, above metal
    vac = kernel.addRectangle(0.0, y_metal_top, 0.0, w_box, h_vacuum)

    # Metal regions at y ∈ [0, t_metal]:
    # Left ground plane: x ∈ [0, x_ground_left_inner]
    metal_left_ground = kernel.addRectangle(0.0, y_metal_bot, 0.0,
                                            x_ground_left_inner, t_metal)

    # Left gap: x ∈ [x_ground_left_inner, x_trace_left]
    gap_left = kernel.addRectangle(x_ground_left_inner, y_metal_bot, 0.0,
                                   w_gap, t_metal)

    # Trace: x ∈ [x_trace_left, x_trace_right]
    metal_trace = kernel.addRectangle(x_trace_left, y_metal_bot, 0.0,
                                      w_trace, t_metal)

    # Right gap: x ∈ [x_trace_right, x_ground_right_inner]
    gap_right = kernel.addRectangle(x_trace_right, y_metal_bot, 0.0,
                                    w_gap, t_metal)

    # Right ground plane: x ∈ [x_ground_right_inner, w_box]
    metal_right_ground = kernel.addRectangle(x_ground_right_inner, y_metal_bot, 0.0,
                                             w_ground, t_metal)

    # Fragment everything to create shared boundaries
    all_surfs = [(2, sub), (2, vac), (2, metal_left_ground), (2, gap_left),
                 (2, metal_trace), (2, gap_right), (2, metal_right_ground)]
    _, frag_map = kernel.fragment(all_surfs, [])
    kernel.synchronize()

    # === Classify entities by bounding box ===
    tol = 1e-6

    all_surfaces = gmsh.model.getEntities(2)
    substrate_domains = Int[]
    vacuum_domains = Int[]
    gap_domains = Int[]
    metal_domains = Int[]

    for (dim, tag) in all_surfaces
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        ymid = (ymin + ymax) / 2.0
        xmid = (xmin + xmax) / 2.0

        if ymid < y_sub_top - tol
            # Below substrate surface → substrate
            push!(substrate_domains, tag)
        elseif ymin > y_metal_bot - tol && ymax < y_metal_top + tol
            # In the metal layer (y ∈ [0, t_metal])
            if xmid > x_ground_left_inner - tol && xmid < x_trace_left + tol
                # Left gap
                push!(gap_domains, tag)
            elseif xmid > x_trace_right - tol && xmid < x_ground_right_inner + tol
                # Right gap
                push!(gap_domains, tag)
            else
                # Metal region (trace, left ground, or right ground)
                push!(metal_domains, tag)
            end
        else
            # Above metal → vacuum
            push!(vacuum_domains, tag)
        end
    end

    # Get all 1D curves (boundaries)
    all_curves = gmsh.model.getEntities(1)

    outer_curves = Int[]           # Outer box boundary
    voltage_right_curves = Int[]   # Right gap voltage path at y=0
    current_curves = Int[]         # Current loop around trace

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

        # Voltage path: right gap bottom (y = 0), x ∈ [x_trace_right, x_ground_right_inner]
        if is_horiz && abs(ymid - y_sub_top) < tol &&
           xmin > x_trace_right - tol && xmax < x_ground_right_inner + tol
            push!(voltage_right_curves, tag)
            continue
        end

        # Current loop: edges of gap domains forming a rectangle around the trace
        # Bottom of gaps (y = 0)
        if is_horiz && abs(ymid - y_metal_bot) < tol
            if (xmin > x_ground_left_inner - tol && xmax < x_trace_left + tol) ||
               (xmin > x_trace_right - tol && xmax < x_ground_right_inner + tol)
                push!(current_curves, tag)
                continue
            end
        end
        # Top of gaps (y = t_metal)
        if is_horiz && abs(ymid - y_metal_top) < tol
            if (xmin > x_ground_left_inner - tol && xmax < x_trace_left + tol) ||
               (xmin > x_trace_right - tol && xmax < x_ground_right_inner + tol)
                push!(current_curves, tag)
                continue
            end
        end
        # Vertical edges of gaps (metal-gap interfaces)
        if is_vert && ymin > y_metal_bot - tol && ymax < y_metal_top + tol
            if abs(xmid - x_trace_left) < tol || abs(xmid - x_trace_right) < tol ||
               abs(xmid - x_ground_left_inner) < tol ||
               abs(xmid - x_ground_right_inner) < tol
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

    gap_attr = domain_idx
    gmsh.model.addPhysicalGroup(2, gap_domains, gap_attr, "gap_vacuum")
    domain_idx += 1

    metal_attr = domain_idx
    gmsh.model.addPhysicalGroup(2, metal_domains, metal_attr, "metal")
    domain_idx += 1

    # Boundaries (NO PEC — metals are meshed domains with London equation)
    outer_attr = bdr_idx
    gmsh.model.addPhysicalGroup(1, outer_curves, outer_attr, "outer")
    bdr_idx += 1

    voltage_right_attr = bdr_idx
    if !isempty(voltage_right_curves)
        gmsh.model.addPhysicalGroup(1, voltage_right_curves, voltage_right_attr,
                                    "voltage_right")
        bdr_idx += 1
    end

    current_attr = bdr_idx
    if !isempty(current_curves)
        gmsh.model.addPhysicalGroup(1, current_curves, current_attr, "current_loop")
        bdr_idx += 1
    end

    # === Mesh size control ===
    # Fine mesh inside metal (need to resolve London penetration depth)
    # Distance field from metal-gap interfaces (trace edges)
    trace_edge_curves = Int[]
    for (dim, tag) in all_curves
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        dx = xmax - xmin
        is_vert = dx < tol
        if is_vert && ymid > y_metal_bot - tol && ymid < y_metal_top + tol
            if abs(xmid - x_trace_left) < tol || abs(xmid - x_trace_right) < tol
                push!(trace_edge_curves, tag)
            end
        end
    end

    if isempty(trace_edge_curves)
        trace_edge_curves = current_curves
    end

    # Distance field from trace edges — same strategy as mesh_thick.jl.
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
    gmsh.write(joinpath(@__DIR__, filename))

    # Print summary
    println("=== CPW 2D Cross-Section Mesh (London) ===")
    println("  Trace width:    $w_trace μm")
    println("  Gap width:      $w_gap μm")
    println("  Ground width:   $w_ground μm")
    println("  Metal thickness: $t_metal μm")
    println("  Box size:       $w_box × $(h_substrate + t_metal + h_vacuum) μm")
    println("  Substrate:      domain attr $sub_attr")
    println("  Vacuum:         domain attr $vac_attr")
    println("  Gap vacuum:     domain attr $gap_attr")
    println("  Metal:          domain attr $metal_attr  ($(length(metal_domains)) surfaces)")
    println("  Outer box:      bdr attr $outer_attr")
    if !isempty(voltage_right_curves)
        println("  Voltage (R):    bdr attr $voltage_right_attr")
    end
    if !isempty(current_curves)
        println("  Current loop:   bdr attr $current_attr")
    end
    println("  Mesh file:      $filename")

    attrs = Dict(
        "substrate" => sub_attr,
        "vacuum" => vac_attr,
        "gap" => gap_attr,
        "metal" => metal_attr,
        "outer" => outer_attr,
        "voltage_right" => !isempty(voltage_right_curves) ? voltage_right_attr : -1,
        "current" => !isempty(current_curves) ? current_attr : -1,
    )

    gmsh.finalize()
    return attrs
end

# Run with default parameters
attrs = generate_cpw2d_london_mesh()
