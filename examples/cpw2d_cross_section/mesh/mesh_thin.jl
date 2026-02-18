# CPW 2D cross-section mesh with infinitely thin (zero-thickness) PEC metal.
# The trace and ground are 1D line segments at y=0 (substrate surface).
# The domain is split into substrate (below y=0) and vacuum (above y=0) regions.
#
# Usage: julia mesh_cpw2d_thin.jl

import Gmsh: gmsh

function generate_cpw2d_thin_mesh(;
    w_trace::Float64 = 15.0,
    w_total_cpw::Float64 = 22.0,
    w_ground::Float64 = 500.0,
    h_substrate::Float64 = 525.0,
    h_vacuum::Float64 = 500.0,
    lc_gap::Float64 = 0.1,
    lc_far::Float64 = 50.0,
    mesh_order::Int = 2,
    filename::String = "cpw2d_thin.msh",
    verbose::Int = 0
)
    w_gap = (w_total_cpw - w_trace) / 2.0
    w_box = 2.0 * (w_gap + w_ground) + w_trace
    x_center = w_box / 2.0
    x_trace_left = x_center - w_trace / 2.0
    x_trace_right = x_center + w_trace / 2.0
    x_gl = x_trace_left - w_gap   # Left ground inner edge
    x_gr = x_trace_right + w_gap  # Right ground inner edge

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("cpw2d_thin")

    # Build geometry from individual points and lines, not rectangles,
    # so the y=0 boundary is naturally split into trace, gap, and ground segments.

    # Points (counter-clockwise for substrate, then vacuum)
    # Substrate bottom: y = -h_substrate
    p1  = gmsh.model.occ.addPoint(0.0,     -h_substrate, 0.0)
    p2  = gmsh.model.occ.addPoint(w_box,   -h_substrate, 0.0)
    # y=0 level: split at ground/gap/trace boundaries
    p3  = gmsh.model.occ.addPoint(w_box,   0.0, 0.0)
    p4  = gmsh.model.occ.addPoint(x_gr,    0.0, 0.0)
    p5  = gmsh.model.occ.addPoint(x_trace_right, 0.0, 0.0)
    p6  = gmsh.model.occ.addPoint(x_trace_left,  0.0, 0.0)
    p7  = gmsh.model.occ.addPoint(x_gl,    0.0, 0.0)
    p8  = gmsh.model.occ.addPoint(0.0,     0.0, 0.0)
    # Vacuum top: y = h_vacuum
    p9  = gmsh.model.occ.addPoint(0.0,     h_vacuum, 0.0)
    p10 = gmsh.model.occ.addPoint(w_box,   h_vacuum, 0.0)

    # Lines — substrate boundary (clockwise looking from outside = CCW for interior)
    l_sub_bot   = gmsh.model.occ.addLine(p1, p2)   # Bottom
    l_sub_right = gmsh.model.occ.addLine(p2, p3)   # Right wall
    l_ground_r  = gmsh.model.occ.addLine(p3, p4)   # Right ground at y=0
    l_gap_r     = gmsh.model.occ.addLine(p4, p5)   # Right gap at y=0
    l_trace     = gmsh.model.occ.addLine(p5, p6)   # Trace at y=0
    l_gap_l     = gmsh.model.occ.addLine(p6, p7)   # Left gap at y=0
    l_ground_l  = gmsh.model.occ.addLine(p7, p8)   # Left ground at y=0
    l_sub_left  = gmsh.model.occ.addLine(p8, p1)   # Left wall

    cl_sub = gmsh.model.occ.addCurveLoop([l_sub_bot, l_sub_right, l_ground_r, l_gap_r,
                                           l_trace, l_gap_l, l_ground_l, l_sub_left])
    s_sub = gmsh.model.occ.addPlaneSurface([cl_sub])

    # Vacuum boundary
    l_vac_left  = gmsh.model.occ.addLine(p8, p9)   # Left wall (shared with substrate at p8)
    l_vac_top   = gmsh.model.occ.addLine(p9, p10)  # Top
    l_vac_right = gmsh.model.occ.addLine(p10, p3)  # Right wall (shared at p3)

    # Vacuum uses the same y=0 lines but reversed
    cl_vac = gmsh.model.occ.addCurveLoop([
        -l_ground_l, -l_gap_l, -l_trace, -l_gap_r, -l_ground_r,
        l_vac_right, -l_vac_top, -l_vac_left
    ])

    # Wait, OCC won't let me reverse lines this way. Let me use fragment instead.
    # Actually with OCC, let me just create the vacuum rectangle separately and fragment.

    gmsh.model.occ.remove([(1, l_vac_left), (1, l_vac_top), (1, l_vac_right),
                            (2, s_sub)], true)
    gmsh.model.occ.remove([(0, p9), (0, p10)], true)

    # Start over with a simpler approach: create substrate and vacuum as rectangles,
    # but with the y=0 boundary points embedded.

    gmsh.model.occ.synchronize()
    gmsh.clear()
    gmsh.model.add("cpw2d_thin")

    # Create 5 substrate rectangles side by side, and 5 vacuum rectangles on top.
    # This naturally splits y=0 into ground/gap/trace/gap/ground segments.

    # Substrate regions (y from -h_substrate to 0)
    sub_lg   = gmsh.model.occ.addRectangle(0.0,           -h_substrate, 0.0, x_gl,           h_substrate)
    sub_gapl = gmsh.model.occ.addRectangle(x_gl,          -h_substrate, 0.0, w_gap,          h_substrate)
    sub_tr   = gmsh.model.occ.addRectangle(x_trace_left,  -h_substrate, 0.0, w_trace,        h_substrate)
    sub_gapr = gmsh.model.occ.addRectangle(x_trace_right, -h_substrate, 0.0, w_gap,          h_substrate)
    sub_rg   = gmsh.model.occ.addRectangle(x_gr,          -h_substrate, 0.0, w_box - x_gr,   h_substrate)

    # Vacuum regions (y from 0 to h_vacuum)
    vac_lg   = gmsh.model.occ.addRectangle(0.0,           0.0, 0.0, x_gl,           h_vacuum)
    vac_gapl = gmsh.model.occ.addRectangle(x_gl,          0.0, 0.0, w_gap,          h_vacuum)
    vac_tr   = gmsh.model.occ.addRectangle(x_trace_left,  0.0, 0.0, w_trace,        h_vacuum)
    vac_gapr = gmsh.model.occ.addRectangle(x_trace_right, 0.0, 0.0, w_gap,          h_vacuum)
    vac_rg   = gmsh.model.occ.addRectangle(x_gr,          0.0, 0.0, w_box - x_gr,   h_vacuum)

    # Fragment to merge shared boundaries
    all_surfs = [(2, s) for s in [sub_lg, sub_gapl, sub_tr, sub_gapr, sub_rg,
                                   vac_lg, vac_gapl, vac_tr, vac_gapr, vac_rg]]
    _, frag_map = gmsh.model.occ.fragment(all_surfs, [])
    gmsh.model.occ.synchronize()

    # Classify entities
    tol = 1e-6
    all_surfaces = gmsh.model.getEntities(2)
    substrate_tags = Int[]
    vacuum_tags = Int[]
    for (dim, tag) in all_surfaces
        bb = gmsh.model.getBoundingBox(dim, tag)
        ymid = (bb[2] + bb[5]) / 2.0
        if ymid < -tol
            push!(substrate_tags, tag)
        else
            push!(vacuum_tags, tag)
        end
    end

    all_curves = gmsh.model.getEntities(1)
    pec_trace = Int[]
    pec_ground = Int[]
    gap_curves = Int[]
    outer_curves = Int[]

    for (dim, tag) in all_curves
        bb = gmsh.model.getBoundingBox(dim, tag)
        xmin, ymin, _, xmax, ymax, _ = bb
        xmid = (xmin + xmax) / 2.0
        ymid = (ymin + ymax) / 2.0
        dx = xmax - xmin
        dy = ymax - ymin
        is_horiz = dy < tol
        is_vert = dx < tol

        # Outer box
        if is_vert && (abs(xmin) < tol || abs(xmax - w_box) < tol)
            push!(outer_curves, tag); continue
        end
        if is_horiz && (abs(ymin + h_substrate) < tol || abs(ymax - h_vacuum) < tol)
            push!(outer_curves, tag); continue
        end

        # y=0 horizontal curves
        if is_horiz && abs(ymid) < tol
            # Trace
            if xmin > x_trace_left - tol && xmax < x_trace_right + tol
                push!(pec_trace, tag); continue
            end
            # Ground
            if xmax < x_gl + tol || xmin > x_gr - tol
                push!(pec_ground, tag); continue
            end
            # Gap
            if (xmin > x_gl - tol && xmax < x_trace_left + tol) ||
               (xmin > x_trace_right - tol && xmax < x_gr + tol)
                push!(gap_curves, tag); continue
            end
        end
    end

    # Physical groups
    sub_attr = 1
    gmsh.model.addPhysicalGroup(2, substrate_tags, sub_attr, "substrate")
    vac_attr = 2
    gmsh.model.addPhysicalGroup(2, vacuum_tags, vac_attr, "vacuum")

    bdr_idx = 1
    pec_trace_attr = bdr_idx
    gmsh.model.addPhysicalGroup(1, pec_trace, pec_trace_attr, "pec_trace")
    bdr_idx += 1

    pec_ground_attr = bdr_idx
    gmsh.model.addPhysicalGroup(1, pec_ground, pec_ground_attr, "pec_ground")
    bdr_idx += 1

    outer_attr = bdr_idx
    gmsh.model.addPhysicalGroup(1, outer_curves, outer_attr, "outer")
    bdr_idx += 1

    gap_attr = -1
    if !isempty(gap_curves)
        gap_attr = bdr_idx
        gmsh.model.addPhysicalGroup(1, gap_curves, gap_attr, "gap")
        bdr_idx += 1
    end

    # Mesh size control
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", Float64.(vcat(pec_trace, gap_curves)))
    gmsh.model.mesh.field.setNumber(1, "Sampling", 100)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_gap)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc_far)
    gmsh.model.mesh.field.setNumber(2, "DistMin", w_gap)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 10.0 * w_gap)

    gmsh.model.mesh.field.setAsBackgroundMesh(2)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.setOrder(mesh_order)

    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(joinpath(@__DIR__, filename))

    println("=== CPW 2D Thin Metal Mesh ===")
    println("  Trace: $w_trace μm, Gap: $w_gap μm, Ground: $w_ground μm")
    println("  Substrate: domain attr $sub_attr, Vacuum: domain attr $vac_attr")
    println("  PEC trace: bdr attr $pec_trace_attr ($(length(pec_trace)) curves)")
    println("  PEC ground: bdr attr $pec_ground_attr ($(length(pec_ground)) curves)")
    println("  Outer: bdr attr $outer_attr, Gap: bdr attr $gap_attr")

    gmsh.finalize()
end

generate_cpw2d_thin_mesh()
