# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Unified CPW 2D cross-section mesh generator.
#
# Supports: single-chip or flip-chip, thin or thick metal, tapered or vertical
# sidewalls, rounded or sharp corners, overetch or not, oxide layers or not.
# Optional outer vacuum below L1 substrate (and above L2 in flipchip), and a
# narrow L1 substrate pocket inside the flipchip cutout.
#
# Geometry (all dimensions in μm):
#   Single chip: substrate below y=0, vacuum above.
#   Flip-chip:   L1 substrate below y=0, vacuum gap y=0..gap, L2 substrate above y=gap.
#
# Thick metal:   trapezoidal cross-section with optional sidewall taper and corner rounding.
#                Overetch trenches in the substrate with optional bottom rounding.
# Thin metal:    zero-thickness PEC at the substrate-vacuum interface.
# Oxide layers:  2D annular rings around each PEC shape, fragmented into MS/MA/SA domains.
# Outer vacuum:  optional vacuum regions below L1 substrate (always) and above L2
#                substrate (flipchip only), creating L0_SA / L2_SA interfaces.
# L1 pocket:     optional narrow substrate pocket inside the flipchip cutout,
#                with independent width and depth.

import Gmsh: gmsh

function generate_cpw2d_mesh(;
    # --- Topology ---
    flipchip::Bool       = false,

    # --- CPW geometry (μm) ---
    w_trace::Float64     = 10.0,
    w_gap::Float64       = 6.0,
    w_ground::Float64    = 300.0,

    # --- Flip-chip only ---
    cow::Float64         = 50.0,     # L1 cutout width
    fc_gap::Float64      = 5.0,      # vacuum gap height
    bb_offset::Float64   = 200.0,    # bump bond offset from center
    bb_width::Float64    = 25.0,     # bump bond width

    # --- Substrate ---
    # Pool E onwards: L1 and L2 substrate thicknesses can be independent.
    # h_substrate, if provided, sets both layers (Pool A-D backwards-compat).
    h_substrate::Float64    = 525.0,
    h_substrate_L1::Float64 = -1.0,   # -1 = inherit from h_substrate
    h_substrate_L2::Float64 = -1.0,   # -1 = inherit from h_substrate
    h_vacuum::Float64       = 500.0,  # single-chip only: vacuum height above metal

    # --- Outer vacuum (backside / topside) ---
    h_vacuum_outer::Float64 = 0.0,   # vacuum below L1 substrate (and above L2 substrate for flipchip)

    # --- Flip-chip L1 pocket (narrower than cutout, deeper than overetch typically) ---
    pocket_width::Float64   = 0.0,   # 0 = no pocket
    pocket_depth::Float64   = 0.0,   # 0 = no pocket

    # --- Metal profile (per-layer) ---
    # Pool C onwards: L1 and L2 metal/overetch are independent. Single-chip
    # uses *_L1 only.
    t_metal_L1::Float64      = 0.0,   # L1 metal thickness; 0 = thin (PEC line)
    t_metal_L2::Float64      = 0.0,   # L2 metal thickness (flipchip only)
    overetch_L1::Float64     = 0.0,   # L1 overetch depth
    overetch_L2::Float64     = 0.0,   # L2 overetch depth (flipchip only)
    sidewall_angle::Float64  = 90.0,  # degrees from horizontal (90=vertical)
    r_top::Float64           = 0.0,   # fillet radius at vacuum-facing tips
    r_bot::Float64           = 0.0,   # fillet radius at overetch trench bottom

    # --- Oxide ---
    t_oxide::Float64     = 0.0,       # 0 = no oxide layers

    # --- Mesh ---
    lc_fine::Float64     = 0.03,
    lc_far::Float64      = 60.0,
    mesh_order::Int      = 2,

    # --- Output ---
    filename::String     = "cpw2d.msh",
    verbose::Int         = 0
)
    # For single-chip, t_metal_L2 / overetch_L2 are ignored.
    # Resolve per-layer substrate thicknesses (default to symmetric h_substrate).
    h_sub_L1 = h_substrate_L1 < 0 ? h_substrate : h_substrate_L1
    h_sub_L2 = h_substrate_L2 < 0 ? h_substrate : h_substrate_L2

    thick_L1 = t_metal_L1 > 0
    thick_L2 = flipchip && t_metal_L2 > 0
    thick    = thick_L1 || thick_L2
    if flipchip && (thick_L1 != thick_L2)
        error("Mixed thin/thick metal between L1 and L2 not supported "*
              "(t_metal_L1=$t_metal_L1, t_metal_L2=$t_metal_L2)")
    end
    has_overetch_L1 = overetch_L1 > 0
    has_overetch_L2 = flipchip && overetch_L2 > 0
    has_overetch    = has_overetch_L1 || has_overetch_L2
    has_oxide = t_oxide > 0
    has_taper = thick && sidewall_angle < 90.0 - 1e-6
    has_outer_vac = h_vacuum_outer > 0
    has_pocket = flipchip && pocket_width > 0 && pocket_depth > 0
    if has_pocket && pocket_width >= cow
        error("pocket_width ($pocket_width) must be smaller than cow ($cow)")
    end

    θ = deg2rad(sidewall_angle)
    tp_L1 = thick_L1 ? (sidewall_angle < 89.9 ? t_metal_L1 / tan(θ) : 0.0) : 0.0
    tp_L2 = thick_L2 ? (sidewall_angle < 89.9 ? t_metal_L2 / tan(θ) : 0.0) : 0.0
    to_L1 = has_overetch_L1 ? (sidewall_angle < 89.9 ? overetch_L1 / tan(θ) : 0.0) : 0.0
    to_L2 = has_overetch_L2 ? (sidewall_angle < 89.9 ? overetch_L2 / tan(θ) : 0.0) : 0.0
    oe_L1 = overetch_L1
    oe_L2 = overetch_L2
    tox = t_oxide

    # Box dimensions
    if flipchip
        gap = fc_gap
        w_box = 2.0 * (bb_offset + bb_width/2 + w_ground)
        ext = max(r_top, 1.0) + 1.0
    else
        gap = 0.0
        w_box = w_trace + 2*w_gap + 2*w_ground
        ext = max(r_top, 1.0) + 1.0
    end
    xc = w_box / 2.0

    # Key x-coordinates at substrate surface
    x_tl = xc - w_trace/2;  x_tr = xc + w_trace/2
    x_gl = x_tl - w_gap;    x_gr = x_tr + w_gap

    # Flip-chip specific
    if flipchip
        x_cl = xc - cow/2;  x_cr = xc + cow/2
        x_bll = xc - bb_offset - bb_width/2
        x_blr = xc - bb_offset + bb_width/2
        x_brl = xc + bb_offset - bb_width/2
        x_brr = xc + bb_offset + bb_width/2
        y2b = gap - t_metal_L2
        y2t = gap
        y2oe = gap + oe_L2
        x_pl = xc - pocket_width/2; x_pr = xc + pocket_width/2
    else
        x_cl = x_gl; x_cr = x_gr  # single chip: cutout = full CPW gap region
    end

    # Y coordinates
    # Interfaces between substrate and vacuum:
    #   y = -h_sub_L1            (L1 substrate / L0 outer vacuum)
    #   y = gap + h_sub_L2       (L2 substrate / outer vacuum; flipchip only)
    y_l1_bot = -h_sub_L1
    y_l2_top = flipchip ? gap + h_sub_L2 : t_metal_L1 + h_vacuum
    y_sub_bot = has_outer_vac ? y_l1_bot - h_vacuum_outer : y_l1_bot
    y_vac_top = (flipchip && has_outer_vac) ? y_l2_top + h_vacuum_outer : y_l2_top

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    gmsh.model.add("cpw2d")
    K = gmsh.model.occ

    # =====================================================================
    # Fillet and shape helpers (from mesh_flipchip_thick.jl)
    # =====================================================================
    function fillet(c, p, n, r)
        r <= 0 && return (c, c, c)
        d1 = let v=(c[1]-p[1],c[2]-p[2]); L=hypot(v...); (v[1]/L,v[2]/L) end
        d2 = let v=(n[1]-c[1],n[2]-c[2]); L=hypot(v...); (v[1]/L,v[2]/L) end
        ce = clamp(d1[1]*d2[1]+d1[2]*d2[2],-1.0,1.0)
        hi = (π-acos(ce))/2; t = r/tan(hi)
        Ti = (c[1]-t*d1[1], c[2]-t*d1[2])
        To = (c[1]+t*d2[1], c[2]+t*d2[2])
        nm = (-d1[2], d1[1])
        C = (Ti[1]+r*nm[1], Ti[2]+r*nm[2])
        return (Ti, C, To)
    end

    function shape(corners)
        nc = length(corners)
        fd = [fillet(corners[i][1], corners[mod1(i-1,nc)][1], corners[mod1(i+1,nc)][1],
                     corners[i][2]) for i in 1:nc]
        curves = Int32[]
        for i in 1:nc
            j = mod1(i+1,nc)
            p1 = K.addPoint(fd[i][3][1], fd[i][3][2], 0.0)
            p2 = K.addPoint(fd[j][1][1], fd[j][1][2], 0.0)
            push!(curves, K.addLine(p1, p2))
            if corners[j][2] > 0
                pc = K.addPoint(fd[j][2][1], fd[j][2][2], 0.0)
                p3 = K.addPoint(fd[j][3][1], fd[j][3][2], 0.0)
                push!(curves, K.addCircleArc(p2, pc, p3))
            end
        end
        return K.addPlaneSurface([K.addCurveLoop(curves)])
    end

    # =====================================================================
    # Oxide ring helpers (from mesh_flipchip_oxide_v3.jl)
    # =====================================================================
    function offset_corners(corners, d, sign; fixed_edges=Set{Int}())
        nc = length(corners)
        edge_n = Vector{Tuple{Float64,Float64}}(undef, nc)
        edge_t = Vector{Tuple{Float64,Float64}}(undef, nc)
        for i in 1:nc
            j = mod1(i+1, nc)
            dx = corners[j][1][1] - corners[i][1][1]
            dy = corners[j][1][2] - corners[i][1][2]
            L = hypot(dx, dy)
            edge_t[i] = (dx/L, dy/L)
            edge_n[i] = i in fixed_edges ? (0.0, 0.0) : (sign*dy/L, sign*-dx/L)
        end
        result = Vector{Tuple{Tuple{Float64,Float64},Float64}}(undef, nc)
        for i in 1:nc
            pt, r = corners[i]
            ei = mod1(i-1, nc); eo = i
            n_in = edge_n[ei]; n_out = edge_n[eo]
            in_fixed = ei in fixed_edges; out_fixed = eo in fixed_edges
            bx = n_in[1]+n_out[1]; by = n_in[2]+n_out[2]; blen = hypot(bx, by)
            if blen < 1e-12
                result[i] = (pt, r > 0 ? max(r+sign*d, 0.0) : 0.0); continue
            end
            if in_fixed != out_fixed
                # One edge fixed: slide along the fixed edge to maintain offset from non-fixed edge
                fe = in_fixed ? ei : eo
                nfe = in_fixed ? n_out : n_in
                ft = edge_t[fe]
                dot_nt = nfe[1]*ft[1] + nfe[2]*ft[2]
                if abs(dot_nt) > 1e-12
                    omag = d / dot_nt
                else
                    omag = 0.0
                end
                new_pt = (pt[1]+omag*ft[1], pt[2]+omag*ft[2])
            else
                ref = hypot(n_in...) > 0.5 ? n_in : n_out
                cos_half = max((bx*ref[1]+by*ref[2])/blen, 0.1)
                omag = d / cos_half
                new_pt = (pt[1]+omag*bx/blen, pt[2]+omag*by/blen)
            end
            new_r = r > 0 ? max(r+sign*d, 0.0) : 0.0
            result[i] = (new_pt, new_r)
        end
        return result
    end

    function shape_with_ref(corners, ref_corners)
        nc = length(corners)
        fd = [fillet(corners[i][1], corners[mod1(i-1,nc)][1], corners[mod1(i+1,nc)][1],
                     corners[i][2]) for i in 1:nc]
        fd_ref = [fillet(ref_corners[i][1], ref_corners[mod1(i-1,nc)][1],
                         ref_corners[mod1(i+1,nc)][1], ref_corners[i][2]) for i in 1:nc]
        for i in 1:nc
            r = corners[i][2]; r > 0 || continue
            C = fd_ref[i][2]
            Ti_old = fd[i][1]; To_old = fd[i][3]
            di = (Ti_old[1]-C[1], Ti_old[2]-C[2]); Li = hypot(di...)
            do_ = (To_old[1]-C[1], To_old[2]-C[2]); Lo = hypot(do_...)
            Li < 1e-12 && continue; Lo < 1e-12 && continue
            fd[i] = ((C[1]+r*di[1]/Li, C[2]+r*di[2]/Li), C,
                     (C[1]+r*do_[1]/Lo, C[2]+r*do_[2]/Lo))
        end
        curves = Int32[]
        for i in 1:nc
            j = mod1(i+1,nc)
            Ti_out = fd[i][3]; Tj_in = fd[j][1]
            d_line = hypot(Ti_out[1]-Tj_in[1], Ti_out[2]-Tj_in[2])
            p1 = K.addPoint(Ti_out[1], Ti_out[2], 0.0)
            p2 = K.addPoint(Tj_in[1], Tj_in[2], 0.0)
            d_line > 1e-10 && push!(curves, K.addLine(p1, p2))
            if corners[j][2] > 0
                Cj = fd[j][2]; Tj_out = fd[j][3]
                d_arc = hypot(Tj_in[1]-Tj_out[1], Tj_in[2]-Tj_out[2])
                if d_arc > 1e-10
                    pc = K.addPoint(Cj[1], Cj[2], 0.0)
                    p3 = K.addPoint(Tj_out[1], Tj_out[2], 0.0)
                    push!(curves, K.addCircleArc(p2, pc, p3))
                end
            end
        end
        return K.addPlaneSurface([K.addCurveLoop(curves)])
    end

    function make_ring(inner_corners, outer_corners; inward=false)
        inner_s = shape(inner_corners)
        outer_s = shape_with_ref(outer_corners, inner_corners)
        if inward
            ring, _ = K.cut([(2, inner_s)], [(2, outer_s)], -1, true, true)
        else
            ring, _ = K.cut([(2, outer_s)], [(2, inner_s)], -1, true, true)
        end
        return [t for (d, t) in ring]
    end

    # =====================================================================
    # Build geometry
    # =====================================================================

    # --- Domain rectangles ---
    # Outer vacuum regions (below L1 substrate, and above L2 substrate for flipchip).
    outer_vac_surfs = Tuple{Int32,Int32}[]
    if has_outer_vac
        ov_bot = K.addRectangle(0.0, y_sub_bot, 0.0, w_box, h_vacuum_outer)
        push!(outer_vac_surfs, (2, ov_bot))
        if flipchip
            ov_top = K.addRectangle(0.0, y_l2_top, 0.0, w_box, h_vacuum_outer)
            push!(outer_vac_surfs, (2, ov_top))
        end
    end

    if flipchip
        l1_sub = K.addRectangle(0.0, y_l1_bot, 0.0, w_box, h_sub_L1)
        vac = K.addRectangle(0.0, 0.0, 0.0, w_box, gap)
        l2_sub = K.addRectangle(0.0, gap, 0.0, w_box, h_sub_L2)
    else
        l1_sub = K.addRectangle(0.0, y_l1_bot, 0.0, w_box, h_sub_L1)
        vac = K.addRectangle(0.0, 0.0, 0.0, w_box, y_l2_top)
    end

    # Flipchip L1 pocket: rectangle cut out of L1 substrate, between y=0 and
    # y=-pocket_depth, centered on xc with width pocket_width. Treated as
    # a vacuum surface that replaces substrate in this region.
    pocket_surf = 0
    if has_pocket
        pocket_surf = K.addRectangle(x_pl, -pocket_depth, 0.0, pocket_width, pocket_depth)
    end

    # --- PEC shapes and boolean operations ---
    vac_result = Tuple{Int32,Int32}[]
    oe_surfs = Tuple{Int32,Int32}[]
    pec_corners = Dict{String, Vector{Tuple{Tuple{Float64,Float64},Float64}}}()

    if thick
        # === THICK METAL ===
        if flipchip
            # L2 trace
            c = [((x_tl, y2t), 0.0), ((x_tl+tp_L2, y2b), r_top),
                 ((x_tr-tp_L2, y2b), r_top), ((x_tr, y2t), 0.0)]
            pec_corners["L2_trace"] = c
            l2_trace = shape(c)

            # L2 ground left
            c = [((-ext, y2t), 0.0), ((-ext, y2b), 0.0),
                 ((x_gl-tp_L2, y2b), r_top), ((x_gl, y2t), 0.0)]
            pec_corners["L2_gl"] = c
            l2_gl = shape(c)

            # L2 ground right
            c = [((x_gr, y2t), 0.0), ((x_gr+tp_L2, y2b), r_top),
                 ((w_box+ext, y2b), 0.0), ((w_box+ext, y2t), 0.0)]
            pec_corners["L2_gr"] = c
            l2_gr = shape(c)

            # L1 ground left
            c = [((-ext, 0.0), 0.0), ((x_cl, 0.0), 0.0),
                 ((x_cl-tp_L1, t_metal_L1), r_top), ((-ext, t_metal_L1), 0.0)]
            pec_corners["L1_gl"] = c
            l1_gl = shape(c)

            # L1 ground right
            c = [((x_cr, 0.0), 0.0), ((w_box+ext, 0.0), 0.0),
                 ((w_box+ext, t_metal_L1), 0.0), ((x_cr+tp_L1, t_metal_L1), r_top)]
            pec_corners["L1_gr"] = c
            l1_gr = shape(c)

            # Bump bonds
            bb_l = K.addRectangle(x_bll, 0.0, 0.0, bb_width, gap)
            bb_r = K.addRectangle(x_brl, 0.0, 0.0, bb_width, gap)

            # Cut PEC from vacuum
            vr, _ = K.cut([(2, vac)],
                [(2, l2_trace), (2, l2_gl), (2, l2_gr),
                 (2, l1_gl), (2, l1_gr), (2, bb_l), (2, bb_r)])
            vac_result = vr

            # Overetch shapes (L1 and L2 are independent now).
            oe_surfs_list = Tuple{Int32,Int32}[]
            if has_overetch_L2
                c = [((x_tr, y2t), 0.0), ((x_gr, y2t), 0.0),
                     ((x_gr-to_L2, y2oe), r_bot), ((x_tr+to_L2, y2oe), r_bot)]
                pec_corners["L2_oe_r"] = c
                l2_oe_r = shape(c)
                c = [((x_gl, y2t), 0.0), ((x_tl, y2t), 0.0),
                     ((x_tl-to_L2, y2oe), r_bot), ((x_gl+to_L2, y2oe), r_bot)]
                pec_corners["L2_oe_l"] = c
                l2_oe_l = shape(c)
                push!(oe_surfs_list, (2, l2_oe_r))
                push!(oe_surfs_list, (2, l2_oe_l))
            end
            if has_overetch_L1
                c = [((x_cl+to_L1, -oe_L1), r_bot), ((x_cr-to_L1, -oe_L1), r_bot),
                     ((x_cr, 0.0), 0.0), ((x_cl, 0.0), 0.0)]
                pec_corners["L1_oe"] = c
                l1_oe = shape(c)
                push!(oe_surfs_list, (2, l1_oe))
            end
            oe_surfs = oe_surfs_list

            # Fragment overetch with substrates and any extras (outer vacuum, pocket)
            extras = Tuple{Int32,Int32}[]
            append!(extras, outer_vac_surfs)
            has_pocket && push!(extras, (2, pocket_surf))
            sub_oe = vcat(vac_result, [(2, l1_sub), (2, l2_sub)], oe_surfs, extras)
            _, _ = K.fragment(sub_oe, [])
            K.synchronize()

        else
            # Single chip, thick metal — L1 only.
            # Trace
            c = [((x_tl, 0.0), 0.0), ((x_tr, 0.0), 0.0),
                 ((x_tr-tp_L1, t_metal_L1), r_top), ((x_tl+tp_L1, t_metal_L1), r_top)]
            pec_corners["trace"] = c
            trace_surf = shape(c)

            # Ground left
            c = [((-ext, 0.0), 0.0), ((x_gl, 0.0), 0.0),
                 ((x_gl-tp_L1, t_metal_L1), r_top), ((-ext, t_metal_L1), 0.0)]
            pec_corners["ground_l"] = c
            gl_surf = shape(c)

            # Ground right
            c = [((x_gr, 0.0), 0.0), ((w_box+ext, 0.0), 0.0),
                 ((w_box+ext, t_metal_L1), 0.0), ((x_gr+tp_L1, t_metal_L1), r_top)]
            pec_corners["ground_r"] = c
            gr_surf = shape(c)

            # Cut PEC from vacuum
            vr, _ = K.cut([(2, vac)], [(2, trace_surf), (2, gl_surf), (2, gr_surf)])
            vac_result = vr

            # Overetch (CCW: TL → BL → BR → TR matching mesh_overetch.jl)
            if has_overetch_L1
                c_oe_r = [((x_tr, 0.0), 0.0),
                          ((x_tr+to_L1, -oe_L1), r_bot), ((x_gr-to_L1, -oe_L1), r_bot),
                          ((x_gr, 0.0), 0.0)]
                pec_corners["oe_r"] = c_oe_r
                oe_r = shape(c_oe_r)
                c_oe_l = [((x_gl, 0.0), 0.0),
                          ((x_gl+to_L1, -oe_L1), r_bot), ((x_tl-to_L1, -oe_L1), r_bot),
                          ((x_tl, 0.0), 0.0)]
                pec_corners["oe_l"] = c_oe_l
                oe_l = shape(c_oe_l)
                oe_surfs = [(2, oe_r), (2, oe_l)]
            end

            extras = Tuple{Int32,Int32}[]
            append!(extras, outer_vac_surfs)
            sub_oe = vcat(vac_result, [(2, l1_sub)], oe_surfs, extras)
            _, _ = K.fragment(sub_oe, [])
            K.synchronize()
        end

    else
        # === THIN METAL (PEC lines at interface) ===
        if flipchip
            # Bump bonds
            bb_l = K.addRectangle(x_bll, 0.0, 0.0, bb_width, gap)
            bb_r = K.addRectangle(x_brl, 0.0, 0.0, bb_width, gap)
            vr, _ = K.cut([(2, vac)], [(2, bb_l), (2, bb_r)])
            vac_result = vr

            # Splitting lines inside the gap
            split_x = sort(unique(filter(x -> x > 0 && x < w_box, [
                x_cl, x_cr, x_gl, x_tl, x_tr, x_gr])))
            split_lines = Tuple{Int,Int}[]
            for x in split_x
                p1 = K.addPoint(x, 0.0, 0.0)
                p2 = K.addPoint(x, gap, 0.0)
                push!(split_lines, (1, K.addLine(p1, p2)))
            end
            extras = Tuple{Int,Int}[]
            append!(extras, outer_vac_surfs)
            has_pocket && push!(extras, (2, pocket_surf))
            all_ents = vcat(vac_result, split_lines, [(2, l1_sub), (2, l2_sub)], extras)
            _, _ = K.fragment(all_ents, [])
            K.synchronize()
        else
            # Single chip thin: fragment with splitting lines
            split_x = sort(unique(filter(x -> x > 0 && x < w_box, [x_gl, x_tl, x_tr, x_gr])))
            split_lines = Tuple{Int,Int}[]
            for x in split_x
                p1 = K.addPoint(x, 0.0, 0.0)
                p2 = K.addPoint(x, y_l2_top, 0.0)
                push!(split_lines, (1, K.addLine(p1, p2)))
            end
            extras = Tuple{Int,Int}[]
            append!(extras, outer_vac_surfs)
            all_ents = vcat([(2, vac), (2, l1_sub)], split_lines, extras)
            _, _ = K.fragment(all_ents, [])
            K.synchronize()
        end
    end

    # =====================================================================
    # Oxide rings (if t_oxide > 0 and thick metal)
    # =====================================================================
    oxide_surfs = Tuple{Int, Symbol}[]
    all_oxide_tags = Set{Int}()

    if has_oxide && thick
        if flipchip
            # Metal rings (outward)
            for (nm, fe) in [("L2_trace", Set{Int}()),
                             ("L2_gl", Set([1])), ("L2_gr", Set([3])),
                             ("L1_gl", Set([4])), ("L1_gr", Set([2]))]
                haskey(pec_corners, nm) || continue
                # For ground planes: use box-clamped corners for the ring
                c = pec_corners[nm]
                ring_c = c
                if nm in ("L2_gl", "L2_gr", "L1_gl", "L1_gr")
                    ring_c = [((clamp(pt[1], 0.0, w_box), pt[2]), r) for (pt, r) in c]
                end
                lbl = startswith(nm, "L1") ? :metal_L1 : :metal_L2
                outer = offset_corners(ring_c, tox, +1; fixed_edges=fe)
                for s in make_ring(ring_c, outer)
                    push!(oxide_surfs, (s, lbl))
                end
            end
            # Overetch rings (inward, opening edge fixed)
            if has_overetch
                for (nm, fe) in [("L2_oe_r", Set([1])), ("L2_oe_l", Set([1])),
                                 ("L1_oe", Set([3]))]
                    haskey(pec_corners, nm) || continue
                    c = pec_corners[nm]
                    lbl = startswith(nm, "L1") ? :SA_L1 : :SA_L2
                    inner = offset_corners(c, tox, -1; fixed_edges=fe)
                    for s in make_ring(c, inner; inward=true)
                        push!(oxide_surfs, (s, lbl))
                    end
                end
            end
        else
            # Single chip metal rings
            for (nm, fe) in [("trace", Set{Int}()),
                             ("ground_l", Set([4])), ("ground_r", Set([2]))]
                haskey(pec_corners, nm) || continue
                c = pec_corners[nm]
                ring_c = c
                if nm in ("ground_l", "ground_r")
                    ring_c = [((clamp(pt[1], 0.0, w_box), pt[2]), r) for (pt, r) in c]
                end
                outer = offset_corners(ring_c, tox, +1; fixed_edges=fe)
                for s in make_ring(ring_c, outer)
                    push!(oxide_surfs, (s, :metal))
                end
            end
            if has_overetch
                for (nm, fe) in [("oe_r", Set([4])), ("oe_l", Set([4]))]
                    haskey(pec_corners, nm) || continue
                    c = pec_corners[nm]
                    inner = offset_corners(c, tox, -1; fixed_edges=fe)
                    for s in make_ring(c, inner; inward=true)
                        push!(oxide_surfs, (s, :SA))
                    end
                end
            else
                # SA oxide strips at y=0 in the CPW gaps (substrate-air interface)
                sa_l = K.addRectangle(x_gl, -tox, 0.0, x_tl - x_gl, tox)
                sa_r = K.addRectangle(x_tr, -tox, 0.0, x_gr - x_tr, tox)
                push!(oxide_surfs, (sa_l, :SA))
                push!(oxide_surfs, (sa_r, :SA))
            end
        end

        # Fragment oxide with existing domains
        if !isempty(oxide_surfs)
            existing = [(2, t) for (_, t) in gmsh.model.getEntities(2)]
            oxide_dt = [(2, s) for (s, _) in oxide_surfs]
            all_for_frag = vcat(oxide_dt, existing)
            _, frag_map = K.fragment(all_for_frag, [])
            K.synchronize()

            n_oxide = length(oxide_dt)
            oxide_frag_raw = Dict{Symbol, Set{Int}}()
            for i in 1:n_oxide
                typ = oxide_surfs[i][2]
                s = get!(oxide_frag_raw, typ, Set{Int}())
                for (d, t) in frag_map[i]; push!(s, t) end
            end
            for s in values(oxide_frag_raw); union!(all_oxide_tags, s) end
        end
    end

    # =====================================================================
    # Domain classification
    # =====================================================================
    tol = 1e-6
    l1_doms = Int[]; l2_doms = Int[]; gap_doms = Int[]

    # Threshold for "substrate vs vacuum" within the L1 pocket region: any
    # surface centered inside the pocket footprint (x∈[x_pl,x_pr]) and below
    # y=0 but above y=-pocket_depth is pocket vacuum, not substrate.
    function _in_pocket(xm, ym)
        has_pocket || return false
        return xm > x_pl - tol && xm < x_pr + tol &&
               ym < -tol && ym > -pocket_depth - tol
    end

    for (dim, tag) in gmsh.model.getEntities(2)
        tag in all_oxide_tags && continue
        bb = gmsh.model.getBoundingBox(dim, tag)
        xm = (bb[1]+bb[4])/2
        ym = (bb[2]+bb[5])/2

        # Outer vacuum (below L1 substrate, and above L2 substrate for flipchip)
        if has_outer_vac
            if ym < y_l1_bot - tol
                push!(gap_doms, tag); continue
            end
            if flipchip && ym > y_l2_top + tol
                push!(gap_doms, tag); continue
            end
        end

        # Pocket cavity inside L1 substrate → vacuum
        if _in_pocket(xm, ym)
            push!(gap_doms, tag); continue
        end

        if flipchip
            if ym < -oe_L1 - tol;        push!(l1_doms, tag)
            elseif ym > gap+oe_L2+tol;   push!(l2_doms, tag)
            else;                         push!(gap_doms, tag) end
        else
            if ym < -oe_L1 - tol;        push!(l1_doms, tag)
            else;                         push!(gap_doms, tag) end
        end
    end

    # Oxide domain classification
    l1_ms_ox=Int[]; l1_ma_ox=Int[]; l2_ms_ox=Int[]; l2_ma_ox=Int[]
    l1_sa_ox=Int[]; l2_sa_ox=Int[]; sc_ms_ox=Int[]; sc_ma_ox=Int[]; sc_sa_ox=Int[]

    if has_oxide && thick
        if flipchip
            l1_sa_set = get(oxide_frag_raw, :SA_L1, Set{Int}())
            l2_sa_set = get(oxide_frag_raw, :SA_L2, Set{Int}())
            for tag in get(oxide_frag_raw, :metal_L1, Set{Int}())
                tag in l1_sa_set && continue
                bb = gmsh.model.getBoundingBox(2, tag)
                ym = (bb[2]+bb[5])/2
                if ym < 0.0; push!(l1_ms_ox, tag) else push!(l1_ma_ox, tag) end
            end
            for tag in get(oxide_frag_raw, :metal_L2, Set{Int}())
                tag in l2_sa_set && continue
                bb = gmsh.model.getBoundingBox(2, tag)
                ym = (bb[2]+bb[5])/2
                if ym > gap; push!(l2_ms_ox, tag) else push!(l2_ma_ox, tag) end
            end
            l1_sa_ox = collect(l1_sa_set)
            l2_sa_ox = collect(l2_sa_set)
        else
            sa_set = get(oxide_frag_raw, :SA, Set{Int}())
            for tag in get(oxide_frag_raw, :metal, Set{Int}())
                tag in sa_set && continue
                bb = gmsh.model.getBoundingBox(2, tag)
                ym = (bb[2]+bb[5])/2
                if ym < 0.0; push!(sc_ms_ox, tag) else push!(sc_ma_ox, tag) end
            end
            sc_sa_ox = collect(sa_set)
        end
    end

    # =====================================================================
    # Boundary classification
    # =====================================================================
    model_surfs = Set(tag for (dim, tag) in gmsh.model.getEntities(2))
    # Substrate surface sets (to distinguish real substrate-vacuum interfaces
    # from internal vacuum-vacuum curves produced by fragmentation).
    substrate_surfs = Set(vcat(l1_doms, l2_doms))
    function _is_sub_vac(adj)
        length(adj) == 2 || return false
        n_sub = count(s -> s in substrate_surfs, adj)
        return n_sub == 1
    end
    mr = max(r_top, r_bot, tp_L1, tp_L2, to_L1, to_L2) + 0.01

    pec_trace = Int[]; pec_ground = Int[]; bb_pec = Int[]; outer_c = Int[]
    sa_curves = Int[]
    # Single-chip thick-metal split buckets. The MS faces provide the nominal
    # thin-sheet edge locations used by the surface-participation calibration.
    sc_trace_ms_curves = Int[]; sc_trace_ma_curves = Int[]
    sc_ground_ms_curves = Int[]; sc_ground_ma_curves = Int[]
    # Flipchip-specific split buckets (layer × role × interface-type × face
    # orientation / feature). Disjoint; together they cover all PEC +
    # substrate-air faces (minus bump bonds). Terminal BC = L2_trace_*;
    # Ground BC = *_ground_*.
    #
    # Metal-substrate (MS) contacts are single horizontal faces -> not split.
    l1_ground_ms_curves = Int[]
    l2_trace_ms_curves  = Int[]
    l2_ground_ms_curves = Int[]
    # Metal-air (MA): split into horizontal (top/bottom) faces and side walls.
    l1_ground_ma_horiz = Int[]; l1_ground_ma_side = Int[]
    l2_trace_ma_horiz  = Int[]; l2_trace_ma_side  = Int[]
    l2_ground_ma_horiz = Int[]; l2_ground_ma_side = Int[]
    # L1 substrate-air: split by feature.
    l1_sa_cutout         = Int[]   # flat cutout floor at y=0 (no-overetch case)
    l1_sa_pocket_side    = Int[]   # pocket vertical walls
    l1_sa_pocket_bottom  = Int[]   # pocket horizontal floor at y=-pocket_depth
    l1_sa_overetch_side  = Int[]   # overetch trench walls
    l1_sa_overetch_floor = Int[]   # overetch trench floor at y=-oe_L1
    # L2 substrate-air: split by feature (no pocket on L2).
    l2_sa_gap            = Int[]   # flat gap surface at y=gap (no-overetch case)
    l2_sa_overetch_side  = Int[]   # overetch trench walls
    l2_sa_overetch_floor = Int[]   # overetch trench floor at y=gap+oe_L2
    l0_sa_curves = Int[]    # L1 substrate / outer vacuum below (y=y_l1_bot)
    l3_sa_curves = Int[]    # L2 substrate / outer vacuum above (y=y_l2_top; flipchip only)

    # For thick metal with oxide: PEC = curve with exactly 1 adjacent model surface (oxide inner face)
    # For thick metal without oxide: classify by geometry (like mesh_flipchip_thick.jl)
    # For thin metal: classify by x-range at interface y

    for (dim, tag) in gmsh.model.getEntities(1)
        up, _ = gmsh.model.getAdjacencies(1, tag)
        adj = [s for s in up if s in model_surfs]
        isempty(adj) && continue

        bb = gmsh.model.getBoundingBox(dim, tag)
        xn,yn,_,xx,yx,_ = bb
        xm=(xn+xx)/2; ym=(yn+yx)/2; dx=xx-xn; dy=yx-yn
        ih = dy<tol; iv = dx<tol

        # Outer box
        if iv && (abs(xn)<tol || abs(xx-w_box)<tol); push!(outer_c,tag); continue end
        if ih && (abs(yn-y_sub_bot)<tol || abs(yx-y_vac_top)<tol); push!(outer_c,tag); continue end

        # L0_SA / L3_SA: horizontal substrate-vacuum interfaces with outer vacuum
        if has_outer_vac && ih && abs(ym - y_l1_bot) < tol
            push!(l0_sa_curves, tag); continue
        end
        if flipchip && has_outer_vac && ih && abs(ym - y_l2_top) < tol
            push!(l3_sa_curves, tag); continue
        end

        # Pocket walls (flipchip): substrate-vacuum interface at x=x_pl, x=x_pr
        # for y∈[-pocket_depth, 0], and at y=-pocket_depth for x∈[x_pl, x_pr].
        # Pocket is in L1 substrate, so classify as L1_SA — but only if this
        # curve actually separates substrate from vacuum. Inside the overetch
        # trench (y∈[-oe,0]) the pocket-wall locus is vacuum-vacuum and must be
        # skipped.
        if has_pocket
            in_pocket_y = yn > -pocket_depth - tol && yx < tol
            in_pocket_x = xn > x_pl - tol && xx < x_pr + tol
            is_pocket_vert = iv && in_pocket_y &&
                             (abs(xm - x_pl) < tol || abs(xm - x_pr) < tol)
            is_pocket_bot  = ih && in_pocket_x && abs(ym + pocket_depth) < tol
            if is_pocket_vert || is_pocket_bot
                if _is_sub_vac(adj)
                    if is_pocket_bot
                        push!(l1_sa_pocket_bottom, tag)
                    else
                        push!(l1_sa_pocket_side, tag)
                    end
                end
                continue
            end
        end

        if has_oxide && thick
            n_adj = length(adj)
            n_ox = count(s -> s in all_oxide_tags, adj)
            # PEC: exactly 1 adjacent surface, and it's oxide
            if n_adj == 1 && n_ox == 1
                if flipchip
                    in_trace_x = xn > x_tl - mr && xx < x_tr + mr
                    in_l2_y = yn > y2b - tol && yx < y2t + tol
                    if in_l2_y && in_trace_x
                        push!(pec_trace, tag)
                    else
                        push!(pec_ground, tag)
                    end
                else
                    in_trace_x = xn > x_tl - mr && xx < x_tr + mr
                    in_metal_y = yn > -tol && yx < t_metal_L1 + tol
                    if in_trace_x && in_metal_y
                        push!(pec_trace, tag)
                    else
                        push!(pec_ground, tag)
                    end
                end
                continue
            end
            # Bump bonds (flipchip): 1 adj, not oxide
            if flipchip && n_adj == 1 && n_ox == 0
                in_bbl = xn>x_bll-tol && xx<x_blr+tol
                in_bbr = xn>x_brl-tol && xx<x_brr+tol
                if (in_bbl||in_bbr) && yn>-tol && yx<gap+tol
                    push!(bb_pec,tag); continue
                end
            end

        elseif thick
            # Thick metal, no oxide: geometric classification.
            # PEC vs SA split: a curve is a metal PEC iff it's in direct contact
            # with a metal face (y in [0,t_metal] for L1 or [y2b,y2t] for L2).
            # An SA curve is a substrate↔vacuum interface with no metal touching.
            if flipchip
                in_bbl = xn>x_bll-tol && xx<x_blr+tol
                in_bbr = xn>x_brl-tol && xx<x_brr+tol
                in_bb = in_bbl || in_bbr

                in_trace_x = xn > x_tl - mr && xx < x_tr + mr
                in_lgap_x = xn > x_gl-mr && xx < x_tl+mr
                in_rgap_x = xn > x_tr-mr && xx < x_gr+mr
                in_gap_x = in_lgap_x || in_rgap_x
                in_cut_x = xn > x_cl-mr && xx < x_cr+mr

                # --- Bump bonds ---
                if in_bb && yn>-tol && yx<gap+tol
                    push!(bb_pec,tag); continue
                end

                # --- SA classification first (to avoid scooping by ground checks) ---
                # Only accept substrate↔vacuum curves (filters out internal
                # vacuum-vacuum segments from fragmentation, e.g. when a pocket
                # overlaps the overetch trench).
                if has_overetch_L1
                    # L1 overetch trench surfaces (inside cutout, y∈[-oe_L1,0]).
                    # Floor = horizontal at y=-oe_L1; walls = everything else.
                    if yn > -oe_L1-tol && yx < tol && in_cut_x && !(ih && abs(ym)<tol)
                        if _is_sub_vac(adj)
                            if ih && abs(ym + oe_L1) < tol
                                push!(l1_sa_overetch_floor, tag)
                            else
                                push!(l1_sa_overetch_side, tag)
                            end
                        end
                        continue
                    end
                end
                if has_overetch_L2
                    # L2 overetch trench surfaces (above L2, y∈[gap, gap+oe_L2]).
                    # Floor = horizontal at y=gap+oe_L2; walls = everything else.
                    if yn > y2t-tol && yx < y2oe+tol && in_gap_x &&
                       !(ih && abs(ym-gap)<tol)
                        if _is_sub_vac(adj)
                            if ih && abs(ym - y2oe) < tol
                                push!(l2_sa_overetch_floor, tag)
                            else
                                push!(l2_sa_overetch_side, tag)
                            end
                        end
                        continue
                    end
                end
                # Horizontal substrate↔vacuum interface at y=0 in the L1 cutout but
                # OUTSIDE the metal footprint (not under L1 grounds). With overetch
                # this layer is vacuum-vacuum and filtered out.
                if ih && abs(ym)<tol && in_cut_x && !in_bb
                    # _is_sub_vac filters out the vacuum-vacuum segment over a
                    # pocket opening (gap above, pocket cavity below): that y=0
                    # line is NOT a substrate-air interface.
                    if !has_overetch_L1 && _is_sub_vac(adj)
                        push!(l1_sa_cutout, tag)
                    end
                    continue
                end
                # Horizontal substrate↔vacuum at y=gap in CPW gap regions (L2 side).
                if ih && abs(ym-gap)<tol && in_gap_x && !in_bb
                    if !has_overetch_L2 && _is_sub_vac(adj)
                        push!(l2_sa_gap, tag)
                    end
                    continue
                end

                # --- Metal PEC surfaces ---
                # Classify each PEC curve by layer × role × face direction:
                #   L1 ground: y∈[0,t_metal]; MS=horizontal at y=0, MA=elsewhere
                #   L2 trace/ground: y∈[y2b,y2t]; MS=horizontal at y=y2t=gap,
                #                                 MA=elsewhere (bottom + sidewalls)

                # L2 trace. MS = horizontal contact at y=gap; MA horizontal =
                # bottom face at y=y2b; MA side = tapered/vertical sidewalls.
                if yn > y2b-tol && yx < y2t+tol && in_trace_x && !in_bb
                    if ih && abs(ym - gap) < tol
                        push!(l2_trace_ms_curves, tag)
                    elseif ih
                        push!(l2_trace_ma_horiz, tag)
                    else
                        push!(l2_trace_ma_side, tag)
                    end
                    continue
                end
                # L2 ground (metal at gap-t_metal..gap, outside trace, outside BB)
                if yn > y2b-tol && yx < y2t+tol && !in_trace_x && !in_bb
                    if ih && abs(ym - gap) < tol
                        push!(l2_ground_ms_curves, tag)
                    elseif ih
                        push!(l2_ground_ma_horiz, tag)
                    else
                        push!(l2_ground_ma_side, tag)
                    end
                    continue
                end
                # L1 ground (metal at 0..t_metal_L1), excluding y=0 inside cutout.
                # MS = horizontal contact at y=0; MA horizontal = top at y=t_metal_L1;
                # MA side = tapered/vertical sidewalls.
                if yn > -tol && yx < t_metal_L1+tol && !in_bb
                    if !(ih && abs(ym)<tol && in_cut_x)
                        if ih && abs(ym) < tol
                            push!(l1_ground_ms_curves, tag)
                        elseif ih
                            push!(l1_ground_ma_horiz, tag)
                        else
                            push!(l1_ground_ma_side, tag)
                        end
                        continue
                    end
                end
            else
                # Single chip thick, no oxide
                in_trace_x = xn > x_tl-mr && xx < x_tr+mr
                in_rgap = xn > x_tr-mr && xx < x_gr+mr
                in_lgap = xn > x_gl-mr && xx < x_tl+mr
                if yn > -tol && yx < t_metal_L1+tol && in_trace_x
                    if ih && abs(ym) < tol
                        push!(sc_trace_ms_curves, tag)
                    else
                        push!(sc_trace_ma_curves, tag)
                    end
                    continue
                end
                if ih && abs(ym)<tol && (in_rgap||in_lgap)
                    if has_overetch_L1
                        # y=0 in gap is vacuum-vacuum interface, not SA
                    else
                        push!(sa_curves, tag)
                    end
                    continue
                end
                if yn > -tol && yx < t_metal_L1+tol && !in_trace_x
                    if ih && abs(ym) < tol
                        push!(sc_ground_ms_curves, tag)
                    else
                        push!(sc_ground_ma_curves, tag)
                    end
                    continue
                end
                if has_overetch_L1
                    if yn > -oe_L1-tol && yx < tol && (in_rgap||in_lgap)
                        push!(sa_curves, tag); continue end
                end
            end

        else
            # Thin metal: classify by x-range at interface y
            if flipchip
                # y=0 (L1)
                if ih && abs(ym)<tol
                    in_bbl = xn>x_bll-tol && xx<x_blr+tol
                    in_bbr = xn>x_brl-tol && xx<x_brr+tol
                    in_cutout = xn>x_cl-tol && xx<x_cr+tol
                    if in_bbl||in_bbr; push!(bb_pec,tag)
                    elseif in_cutout; push!(sa_curves,tag)
                    else; push!(pec_ground,tag) end
                    continue
                end
                # y=gap (L2)
                if ih && abs(ym-gap)<tol
                    in_bbl = xn>x_bll-tol && xx<x_blr+tol
                    in_bbr = xn>x_brl-tol && xx<x_brr+tol
                    in_trace = xn>x_tl-tol && xx<x_tr+tol
                    in_gap_l = xn>x_gl-tol && xx<x_tl+tol
                    in_gap_r = xn>x_tr-tol && xx<x_gr+tol
                    if in_bbl||in_bbr; push!(bb_pec,tag)
                    elseif in_trace; push!(pec_trace,tag)
                    elseif in_gap_l||in_gap_r; push!(sa_curves,tag)
                    else; push!(pec_ground,tag) end
                    continue
                end
                # Vertical BB sides
                if iv && yn>-tol && yx<gap+tol
                    near_bb = any(abs(xm-bx)<tol for bx in [x_bll,x_blr,x_brl,x_brr])
                    if near_bb; push!(bb_pec,tag); continue end
                end
            else
                # Single chip thin: y=0 interface
                if ih && abs(ym)<tol
                    if xn>=x_tl-tol && xx<=x_tr+tol; push!(pec_trace,tag)
                    elseif xx<=x_gl+tol || xn>=x_gr-tol; push!(pec_ground,tag)
                    elseif (xn>=x_gl-tol && xx<=x_tl+tol) || (xn>=x_tr-tol && xx<=x_gr+tol)
                        push!(sa_curves,tag)
                    end
                    continue
                end
            end
        end
    end

    # =====================================================================
    # Physical groups
    # =====================================================================
    di = 1
    if flipchip
        gmsh.model.addPhysicalGroup(2, l1_doms, di, "L1_substrate"); di+=1
        gmsh.model.addPhysicalGroup(2, l2_doms, di, "L2_substrate"); di+=1
    else
        gmsh.model.addPhysicalGroup(2, l1_doms, di, "substrate"); di+=1
    end
    gmsh.model.addPhysicalGroup(2, gap_doms, di, "vacuum"); di+=1

    if has_oxide && thick
        if flipchip
            for (d, nm) in [(l1_ms_ox,"L1_MS_oxide"),(l1_ma_ox,"L1_MA_oxide"),
                            (l2_ms_ox,"L2_MS_oxide"),(l2_ma_ox,"L2_MA_oxide"),
                            (l1_sa_ox,"L1_SA_oxide"),(l2_sa_ox,"L2_SA_oxide")]
                !isempty(d) && (gmsh.model.addPhysicalGroup(2, d, di, nm); di+=1)
            end
        else
            for (d, nm) in [(sc_ms_ox,"MS_oxide"),(sc_ma_ox,"MA_oxide"),
                            (sc_sa_ox,"SA_oxide")]
                !isempty(d) && (gmsh.model.addPhysicalGroup(2, d, di, nm); di+=1)
            end
        end
    end

    bi = 1
    if !flipchip && thick && !has_oxide
        # Disjoint interface groups for the single-chip fabrication calibration.
        for (d, nm) in [
            (sc_trace_ms_curves,  "trace_MS"),
            (sc_trace_ma_curves,  "trace_MA"),
            (sc_ground_ms_curves, "ground_MS"),
            (sc_ground_ma_curves, "ground_MA"),
            (sa_curves,            "SA"),
        ]
            !isempty(d) && (gmsh.model.addPhysicalGroup(1, d, bi, nm); bi+=1)
        end
    else
        # Legacy groups used by thin metal and the oxide path.
        !isempty(pec_trace) && (gmsh.model.addPhysicalGroup(1,pec_trace,bi,"PEC_trace"); bi+=1)
        !isempty(pec_ground) && (gmsh.model.addPhysicalGroup(1,pec_ground,bi,"PEC_ground"); bi+=1)
        !isempty(bb_pec) && (gmsh.model.addPhysicalGroup(1,bb_pec,bi,"bump_bonds"); bi+=1)
        !isempty(sa_curves) && (gmsh.model.addPhysicalGroup(1,sa_curves,bi,"SA"); bi+=1)
    end
    # Flipchip path: fine-grained split (layer × role × interface × face/feature).
    # MS contacts are single horizontal faces; MA split into horiz/side; SA
    # split by feature (cutout/gap, pocket side/bottom, overetch side/floor).
    # Names are stable identifiers — downstream config maps groups BY NAME,
    # not by attribute number (numbering shifts with pocket/overetch presence).
    for (d, nm) in [
        (l1_ground_ms_curves,   "L1_ground_MS"),
        (l1_ground_ma_horiz,    "L1_ground_MA_horiz"),
        (l1_ground_ma_side,     "L1_ground_MA_side"),
        (l1_sa_cutout,          "L1_SA_cutout"),
        (l1_sa_pocket_side,     "L1_SA_pocket_side"),
        (l1_sa_pocket_bottom,   "L1_SA_pocket_bottom"),
        (l1_sa_overetch_side,   "L1_SA_overetch_side"),
        (l1_sa_overetch_floor,  "L1_SA_overetch_floor"),
        (l2_trace_ms_curves,    "L2_trace_MS"),
        (l2_trace_ma_horiz,     "L2_trace_MA_horiz"),
        (l2_trace_ma_side,      "L2_trace_MA_side"),
        (l2_ground_ms_curves,   "L2_ground_MS"),
        (l2_ground_ma_horiz,    "L2_ground_MA_horiz"),
        (l2_ground_ma_side,     "L2_ground_MA_side"),
        (l2_sa_gap,             "L2_SA_gap"),
        (l2_sa_overetch_side,   "L2_SA_overetch_side"),
        (l2_sa_overetch_floor,  "L2_SA_overetch_floor"),
        (l0_sa_curves,          "L0_SA"),
        (l3_sa_curves,          "L3_SA"),
        (outer_c,               "outer"),
    ]
        !isempty(d) && (gmsh.model.addPhysicalGroup(1, d, bi, nm); bi+=1)
    end

    # =====================================================================
    # Mesh size control
    # =====================================================================
    pts = Float64[]
    if thick
        if flipchip
            for (x,y) in [(x_tl+tp_L2,y2b),(x_tr-tp_L2,y2b),(x_gl-tp_L2,y2b),(x_gr+tp_L2,y2b),
                           (x_cl-tp_L1,t_metal_L1),(x_cr+tp_L1,t_metal_L1)]
                push!(pts, Float64(K.addPoint(x,y,0.0))) end
        else
            for (x,y) in [(x_tl+tp_L1,t_metal_L1),(x_tr-tp_L1,t_metal_L1),
                           (x_gl-tp_L1,t_metal_L1),(x_gr+tp_L1,t_metal_L1)]
                push!(pts, Float64(K.addPoint(x,y,0.0))) end
        end
    else
        if flipchip
            for (x,y) in [(x_tl,gap),(x_tr,gap),(x_gl,gap),(x_gr,gap),
                           (x_cl,0.0),(x_cr,0.0)]
                push!(pts, Float64(K.addPoint(x,y,0.0))) end
        else
            for x in [x_tl, x_tr, x_gl, x_gr]
                push!(pts, Float64(K.addPoint(x, 0.0, 0.0))) end
        end
    end
    K.synchronize()

    field_id = 1
    gmsh.model.mesh.field.add("Distance", field_id)
    gmsh.model.mesh.field.setNumbers(field_id, "PointsList", pts)
    field_id += 1
    gmsh.model.mesh.field.add("Threshold", field_id)
    gmsh.model.mesh.field.setNumber(field_id, "InField", field_id-1)
    gmsh.model.mesh.field.setNumber(field_id, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(field_id, "SizeMax", lc_far)
    dist_scale = flipchip ? gap : w_gap
    gmsh.model.mesh.field.setNumber(field_id, "DistMin",
        thick ? 5.0*max(tox, max(t_metal_L1, t_metal_L2)/50) : 0.1)
    gmsh.model.mesh.field.setNumber(field_id, "DistMax", 200.0*dist_scale)
    thresh_id = field_id
    field_id += 1

    fields_to_min = [Float64(thresh_id)]

    # Oxide curve refinement
    if has_oxide && thick
        oxide_curves = Int[]
        for tag in collect(all_oxide_tags)
            _, down = gmsh.model.getAdjacencies(2, tag)
            append!(oxide_curves, down)
        end
        unique!(oxide_curves)
        lc_oxide = 10.0 * tox
        gmsh.model.mesh.field.add("Distance", field_id)
        gmsh.model.mesh.field.setNumbers(field_id, "CurvesList", Float64.(oxide_curves))
        field_id += 1
        gmsh.model.mesh.field.add("Threshold", field_id)
        gmsh.model.mesh.field.setNumber(field_id, "InField", field_id-1)
        gmsh.model.mesh.field.setNumber(field_id, "SizeMin", lc_oxide)
        gmsh.model.mesh.field.setNumber(field_id, "SizeMax", lc_far)
        gmsh.model.mesh.field.setNumber(field_id, "DistMin", 100.0*tox)
        gmsh.model.mesh.field.setNumber(field_id, "DistMax", 40.0*dist_scale)
        push!(fields_to_min, Float64(field_id))
        field_id += 1
    end

    gmsh.model.mesh.field.add("Min", field_id)
    gmsh.model.mesh.field.setNumbers(field_id, "FieldsList", fields_to_min)
    gmsh.model.mesh.field.setAsBackgroundMesh(field_id)

    for (k,v) in [("Mesh.MeshSizeMin",lc_fine),("Mesh.MeshSizeMax",lc_far),
                   ("Mesh.MeshSizeExtendFromBoundary",0),("Mesh.MeshSizeFromPoints",0),
                   ("Mesh.MeshSizeFromCurvature",0)]
        gmsh.option.setNumber(k,v) end

    gmsh.model.mesh.generate(2)
    gmsh.model.mesh.optimize("Netgen")
    gmsh.model.mesh.setOrder(mesh_order)
    if mesh_order > 1
        gmsh.model.mesh.optimize("HighOrderElastic")
    end

    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    gmsh.option.setNumber("Mesh.Binary",1)
    gmsh.write(joinpath(@__DIR__,filename))

    # Summary
    println("=== CPW 2D Cross-Section Mesh ===")
    println("  Topology: $(flipchip ? "flip-chip" : "single chip")")
    if thick
        if flipchip
            println("  Metal: thick L1=$(t_metal_L1*1000)nm L2=$(t_metal_L2*1000)nm")
        else
            println("  Metal: thick L1=$(t_metal_L1*1000)nm")
        end
        println("  Sidewall: $(sidewall_angle)° r_top=$(r_top*1000)nm r_bot=$(r_bot*1000)nm")
        if flipchip
            (has_overetch_L1 || has_overetch_L2) && println(
                "  Overetch: L1=$(oe_L1*1000)nm L2=$(oe_L2*1000)nm")
        else
            has_overetch_L1 && println("  Overetch: L1=$(oe_L1*1000)nm")
        end
    else
        println("  Metal: thin")
    end
    has_oxide && println("  Oxide: $(tox*1000)nm")
    has_outer_vac && println("  Outer vacuum: $h_vacuum_outer μm")
    has_pocket && println("  L1 pocket: w=$(pocket_width)μm depth=$(pocket_depth)μm")
    println("  w_trace=$w_trace w_gap=$w_gap μm")
    flipchip && println("  cow=$cow gap=$fc_gap μm")
    for (dim, tag) in gmsh.model.getPhysicalGroups(2)
        name = gmsh.model.getPhysicalName(dim, tag)
        ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  Domain $tag: $name ($(length(ents)) surfs)")
    end
    for (dim, tag) in gmsh.model.getPhysicalGroups(1)
        name = gmsh.model.getPhysicalName(dim, tag)
        ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  Boundary $tag: $name ($(length(ents)) curves)")
    end
    println("  Nodes: $(length(gmsh.model.mesh.getNodes()[1]))")
    println("  File: $filename")

    gmsh.finalize()
end
