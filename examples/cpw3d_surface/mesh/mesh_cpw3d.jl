# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Unified CPW 2D/quasi-2D cross-section mesh generator.
#
# Supports: single-chip or flip-chip, thin or thick metal, tapered or vertical
# sidewalls, rounded or sharp corners, overetch or not, oxide layers or not.
# Optional extrusion in the +z direction for quasi-2D (3D) simulations.
#
# Geometry (all dimensions in μm):
#   Single chip: substrate below y=0, vacuum above.
#   Flip-chip:   L1 substrate below y=0, vacuum gap y=0..gap, L2 substrate above y=gap.
#
# Thick metal:   trapezoidal cross-section with optional sidewall taper and corner rounding.
#                Overetch trenches in the substrate with optional bottom rounding.
# Thin metal:    zero-thickness PEC at the substrate-vacuum interface.
# Oxide layers:  2D annular rings around each PEC shape, fragmented into MS/MA/SA domains.
# Extrusion:     if extrude_length > 0, the 2D cross-section is extruded in +z.
#                Physical groups become 3D volumes (domains) and 2D surfaces (boundaries).
#                The two end-cap faces at z=0 and z=extrude_length are added to "outer".

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
    h_substrate::Float64 = 525.0,
    h_vacuum::Float64    = 500.0,    # single-chip only: vacuum height above metal

    # --- Metal profile ---
    t_metal::Float64         = 0.0,   # 0 = thin (PEC line), >0 = thick
    overetch_depth::Float64  = 0.0,   # 0 = no overetch
    sidewall_angle::Float64  = 90.0,  # degrees from horizontal (90=vertical)
    r_top::Float64           = 0.0,   # fillet radius at vacuum-facing tips
    r_bot::Float64           = 0.0,   # fillet radius at overetch trench bottom

    # --- Oxide ---
    t_oxide::Float64     = 0.0,       # 0 = no oxide layers

    # --- Mesh ---
    lc_fine::Float64     = 0.03,
    lc_far::Float64      = 60.0,
    mesh_order::Int      = 2,

    # --- Extrusion (quasi-2D) ---
    extrude_length::Float64 = 0.0,    # 0 = pure 2D, >0 = extrude +z by this length (μm)
    extrude_nz::Int         = 0,      # number of layers along z; 0 = auto (based on lc_fine)

    # --- Output ---
    filename::String     = "cpw2d.msh",
    verbose::Int         = 0
)
    thick = t_metal > 0
    has_overetch = overetch_depth > 0
    has_oxide = t_oxide > 0
    has_taper = thick && sidewall_angle < 90.0 - 1e-6

    θ = deg2rad(sidewall_angle)
    tp = thick ? (sidewall_angle < 89.9 ? t_metal / tan(θ) : 0.0) : 0.0
    to = has_overetch ? (sidewall_angle < 89.9 ? overetch_depth / tan(θ) : 0.0) : 0.0
    oe = overetch_depth
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
        y2b = gap - t_metal
        y2t = gap
        y2oe = gap + oe
    else
        x_cl = x_gl; x_cr = x_gr  # single chip: cutout = full CPW gap region
    end

    # Y coordinates
    y_sub_bot = -h_substrate
    y_vac_top = flipchip ? gap + h_substrate : t_metal + h_vacuum

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
    if flipchip
        l1_sub = K.addRectangle(0.0, -h_substrate, 0.0, w_box, h_substrate)
        vac = K.addRectangle(0.0, 0.0, 0.0, w_box, gap)
        l2_sub = K.addRectangle(0.0, gap, 0.0, w_box, h_substrate)
    else
        l1_sub = K.addRectangle(0.0, y_sub_bot, 0.0, w_box, h_substrate)
        vac = K.addRectangle(0.0, 0.0, 0.0, w_box, y_vac_top)
    end

    # --- PEC shapes and boolean operations ---
    vac_result = Tuple{Int32,Int32}[]
    oe_surfs = Tuple{Int32,Int32}[]
    pec_corners = Dict{String, Vector{Tuple{Tuple{Float64,Float64},Float64}}}()

    if thick
        # === THICK METAL ===
        if flipchip
            # L2 trace
            c = [((x_tl, y2t), 0.0), ((x_tl+tp, y2b), r_top),
                 ((x_tr-tp, y2b), r_top), ((x_tr, y2t), 0.0)]
            pec_corners["L2_trace"] = c
            l2_trace = shape(c)

            # L2 ground left
            c = [((-ext, y2t), 0.0), ((-ext, y2b), 0.0),
                 ((x_gl-tp, y2b), r_top), ((x_gl, y2t), 0.0)]
            pec_corners["L2_gl"] = c
            l2_gl = shape(c)

            # L2 ground right
            c = [((x_gr, y2t), 0.0), ((x_gr+tp, y2b), r_top),
                 ((w_box+ext, y2b), 0.0), ((w_box+ext, y2t), 0.0)]
            pec_corners["L2_gr"] = c
            l2_gr = shape(c)

            # L1 ground left
            c = [((-ext, 0.0), 0.0), ((x_cl, 0.0), 0.0),
                 ((x_cl-tp, t_metal), r_top), ((-ext, t_metal), 0.0)]
            pec_corners["L1_gl"] = c
            l1_gl = shape(c)

            # L1 ground right
            c = [((x_cr, 0.0), 0.0), ((w_box+ext, 0.0), 0.0),
                 ((w_box+ext, t_metal), 0.0), ((x_cr+tp, t_metal), r_top)]
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

            # Overetch shapes
            if has_overetch
                c = [((x_tr, y2t), 0.0), ((x_gr, y2t), 0.0),
                     ((x_gr-to, y2oe), r_bot), ((x_tr+to, y2oe), r_bot)]
                pec_corners["L2_oe_r"] = c
                l2_oe_r = shape(c)
                c = [((x_gl, y2t), 0.0), ((x_tl, y2t), 0.0),
                     ((x_tl-to, y2oe), r_bot), ((x_gl+to, y2oe), r_bot)]
                pec_corners["L2_oe_l"] = c
                l2_oe_l = shape(c)
                c = [((x_cl+to, -oe), r_bot), ((x_cr-to, -oe), r_bot),
                     ((x_cr, 0.0), 0.0), ((x_cl, 0.0), 0.0)]
                pec_corners["L1_oe"] = c
                l1_oe = shape(c)
                oe_surfs = [(2, l1_oe), (2, l2_oe_r), (2, l2_oe_l)]
            end

            # Fragment overetch with substrates
            sub_oe = vcat(vac_result, [(2, l1_sub), (2, l2_sub)], oe_surfs)
            _, _ = K.fragment(sub_oe, [])
            K.synchronize()

        else
            # Single chip, thick metal
            # Trace
            c = [((x_tl, 0.0), 0.0), ((x_tr, 0.0), 0.0),
                 ((x_tr-tp, t_metal), r_top), ((x_tl+tp, t_metal), r_top)]
            pec_corners["trace"] = c
            trace_surf = shape(c)

            # Ground left
            c = [((-ext, 0.0), 0.0), ((x_gl, 0.0), 0.0),
                 ((x_gl-tp, t_metal), r_top), ((-ext, t_metal), 0.0)]
            pec_corners["ground_l"] = c
            gl_surf = shape(c)

            # Ground right
            c = [((x_gr, 0.0), 0.0), ((w_box+ext, 0.0), 0.0),
                 ((w_box+ext, t_metal), 0.0), ((x_gr+tp, t_metal), r_top)]
            pec_corners["ground_r"] = c
            gr_surf = shape(c)

            # Cut PEC from vacuum
            vr, _ = K.cut([(2, vac)], [(2, trace_surf), (2, gl_surf), (2, gr_surf)])
            vac_result = vr

            # Overetch (CCW: TL → BL → BR → TR matching mesh_overetch.jl)
            if has_overetch
                c_oe_r = [((x_tr, 0.0), 0.0),
                          ((x_tr+to, -oe), r_bot), ((x_gr-to, -oe), r_bot),
                          ((x_gr, 0.0), 0.0)]
                pec_corners["oe_r"] = c_oe_r
                oe_r = shape(c_oe_r)
                c_oe_l = [((x_gl, 0.0), 0.0),
                          ((x_gl+to, -oe), r_bot), ((x_tl-to, -oe), r_bot),
                          ((x_tl, 0.0), 0.0)]
                pec_corners["oe_l"] = c_oe_l
                oe_l = shape(c_oe_l)
                oe_surfs = [(2, oe_r), (2, oe_l)]
            end

            sub_oe = vcat(vac_result, [(2, l1_sub)], oe_surfs)
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
            all_ents = vcat(vac_result, split_lines, [(2, l1_sub), (2, l2_sub)])
            _, _ = K.fragment(all_ents, [])
            K.synchronize()
        else
            # Single chip thin: fragment with splitting lines
            split_x = sort(unique(filter(x -> x > 0 && x < w_box, [x_gl, x_tl, x_tr, x_gr])))
            split_lines = Tuple{Int,Int}[]
            for x in split_x
                p1 = K.addPoint(x, 0.0, 0.0)
                p2 = K.addPoint(x, y_vac_top, 0.0)
                push!(split_lines, (1, K.addLine(p1, p2)))
            end
            all_ents = vcat([(2, vac), (2, l1_sub)], split_lines)
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
    # Extrusion (quasi-2D)
    # =====================================================================
    # Do the extrusion here, before classification. Since the geometry is
    # z-invariant, (x,y) bounding boxes can classify both the 2D (dim=2) and
    # 3D (dim=3) entities with the same logic.
    extruded = extrude_length > 0
    # Map from 2D oxide surface tags to 3D oxide volume tags (for classification).
    oxide_2d_to_3d = Dict{Int,Int}()

    if extruded
        all_2d_entities = [(2, tag) for (_, tag) in gmsh.model.getEntities(2)]
        # If extrude_nz > 0, use structured extrusion: the 2D cross-section mesh
        # is swept in z with exactly nz layers, producing prism elements and
        # guaranteeing identical xy resolution at every z. Otherwise use
        # unstructured extrusion (3D mesher decides z resolution via the
        # cross-section size field, isotropically scaled).
        if extrude_nz > 0
            out = K.extrude(all_2d_entities, 0.0, 0.0, extrude_length,
                            [extrude_nz], [1.0], false)
        else
            out = K.extrude(all_2d_entities, 0.0, 0.0, extrude_length,
                            Int[], Float64[], false)
        end
        K.synchronize()

        # Decode: for each input surface, outputs appear in order
        #   (2, top_face), (3, volume), (2, lat_face_1), (2, lat_face_2), ...
        # Build a map from original 2D surface tag → 3D volume tag.
        pos = 1
        for (_, s) in all_2d_entities
            top = out[pos]; pos += 1
            vol = out[pos]; pos += 1
            @assert top[1] == 2 && vol[1] == 3
            if s in all_oxide_tags
                oxide_2d_to_3d[s] = vol[2]
            end
            # skip lateral faces
            _, bdr = gmsh.model.getAdjacencies(2, s)
            pos += length(bdr)
        end
        # After extrusion, oxide tags refer to volumes.
        if !isempty(all_oxide_tags)
            all_oxide_tags = Set{Int}(values(oxide_2d_to_3d))
        end
    end

    # Classification dimension
    dom_dim = extruded ? 3 : 2
    bdr_dim = extruded ? 2 : 1

    # =====================================================================
    # Domain classification (uses (x,y) bounding box, works for 2D or 3D)
    # =====================================================================
    tol = 1e-6
    l1_doms = Int[]; l2_doms = Int[]; gap_doms = Int[]

    for (dim, tag) in gmsh.model.getEntities(dom_dim)
        tag in all_oxide_tags && continue
        bb = gmsh.model.getBoundingBox(dim, tag)
        ym = (bb[2]+bb[5])/2
        if flipchip
            if ym < -oe - tol;        push!(l1_doms, tag)
            elseif ym > gap+oe+tol;   push!(l2_doms, tag)
            else;                      push!(gap_doms, tag) end
        else
            if ym < -oe - tol;        push!(l1_doms, tag)
            else;                      push!(gap_doms, tag) end
        end
    end

    # Oxide domain classification (remap 2D oxide tags to 3D if extruded)
    l1_ms_ox=Int[]; l1_ma_ox=Int[]; l2_ms_ox=Int[]; l2_ma_ox=Int[]
    l1_sa_ox=Int[]; l2_sa_ox=Int[]; sc_ms_ox=Int[]; sc_ma_ox=Int[]; sc_sa_ox=Int[]

    map_ox = extruded ? (t -> oxide_2d_to_3d[t]) : identity

    if has_oxide && thick
        if flipchip
            l1_sa_set = get(oxide_frag_raw, :SA_L1, Set{Int}())
            l2_sa_set = get(oxide_frag_raw, :SA_L2, Set{Int}())
            for tag in get(oxide_frag_raw, :metal_L1, Set{Int}())
                tag in l1_sa_set && continue
                bb = gmsh.model.getBoundingBox(2, tag)  # BB computed on 2D (pre-extrusion)
                ym = (bb[2]+bb[5])/2
                dst = ym < 0.0 ? l1_ms_ox : l1_ma_ox
                push!(dst, map_ox(tag))
            end
            for tag in get(oxide_frag_raw, :metal_L2, Set{Int}())
                tag in l2_sa_set && continue
                bb = gmsh.model.getBoundingBox(2, tag)
                ym = (bb[2]+bb[5])/2
                dst = ym > gap ? l2_ms_ox : l2_ma_ox
                push!(dst, map_ox(tag))
            end
            l1_sa_ox = [map_ox(t) for t in l1_sa_set]
            l2_sa_ox = [map_ox(t) for t in l2_sa_set]
        else
            sa_set = get(oxide_frag_raw, :SA, Set{Int}())
            for tag in get(oxide_frag_raw, :metal, Set{Int}())
                tag in sa_set && continue
                bb = gmsh.model.getBoundingBox(2, tag)
                ym = (bb[2]+bb[5])/2
                dst = ym < 0.0 ? sc_ms_ox : sc_ma_ox
                push!(dst, map_ox(tag))
            end
            sc_sa_ox = [map_ox(t) for t in sa_set]
        end
    end

    # =====================================================================
    # Boundary classification
    # =====================================================================
    model_doms = Set(tag for (dim, tag) in gmsh.model.getEntities(dom_dim))
    mr = max(r_top, r_bot, tp, to) + 0.01

    pec_trace = Int[]; pec_ground = Int[]; bb_pec = Int[]; outer_c = Int[]
    sa_curves = Int[]
    port_front = Int[]   # z = 0 end cap (waveport in 3D)
    port_back  = Int[]   # z = extrude_length end cap (waveport in 3D)

    # For thick metal with oxide: PEC = curve/face with exactly 1 adjacent
    # model domain (oxide inner face).
    # For thick metal without oxide: classify by geometry.
    # For thin metal: classify by x-range at interface y.
    # When extruded, lateral faces use their (x,y) extent which matches the
    # originating 1D curve; end-cap faces (z-extent ≈ 0) go to "outer".

    for (dim, tag) in gmsh.model.getEntities(bdr_dim)
        up, _ = gmsh.model.getAdjacencies(bdr_dim, tag)
        adj = [s for s in up if s in model_doms]
        isempty(adj) && continue

        bb = gmsh.model.getBoundingBox(dim, tag)
        xn,yn,zn,xx,yx,zx = bb
        xm=(xn+xx)/2; ym=(yn+yx)/2; dx=xx-xn; dy=yx-yn
        ih = dy<tol; iv = dx<tol

        # Extrusion end caps (z-extent ≈ 0): these are the z=0 and z=L faces,
        # which are copies of the 2D cross-section. Classified as separate
        # port_front (z=0) and port_back (z=extrude_length) boundaries for
        # use as waveport BCs in RF simulations.
        if extruded && (zx - zn) < tol
            zm = (zn + zx) / 2
            if abs(zm) < tol
                push!(port_front, tag)
            else
                push!(port_back, tag)
            end
            continue
        end

        # Outer box (lateral faces at cross-section boundary)
        if iv && (abs(xn)<tol || abs(xx-w_box)<tol); push!(outer_c,tag); continue end
        if ih && (abs(yn-y_sub_bot)<tol || abs(yx-y_vac_top)<tol); push!(outer_c,tag); continue end

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
                    in_metal_y = yn > -tol && yx < t_metal + tol
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
            # Thick metal, no oxide: geometric classification
            if flipchip
                # Bump bonds
                in_bbl = xn>x_bll-tol && xx<x_blr+tol
                in_bbr = xn>x_brl-tol && xx<x_brr+tol
                if (in_bbl||in_bbr) && yn>-tol && yx<gap+tol
                    push!(bb_pec,tag); continue end

                # L2 trace
                in_trace_x = xn > x_tl - mr && xx < x_tr + mr
                if yn > y2b-tol && yx < y2t+tol && in_trace_x && !(in_bbl||in_bbr)
                    push!(pec_trace, tag); continue end

                # L2 ground
                if yn > y2b-tol && yx < y2t+tol && !in_trace_x && !(in_bbl||in_bbr)
                    push!(pec_ground, tag); continue end

                # L1 ground
                if yn > -tol && yx < t_metal+tol && !(in_bbl||in_bbr)
                    push!(pec_ground, tag); continue end

                # SA (overetch)
                if has_overetch
                    in_cut_x = xn > x_cl-mr && xx < x_cr+mr
                    if yn > -oe-tol && yx < tol && in_cut_x
                        if !(ih && abs(ym)<tol); push!(sa_curves, tag); continue end
                    end
                    in_lgap = xn > x_gl-mr && xx < x_tl+mr
                    in_rgap = xn > x_tr-mr && xx < x_gr+mr
                    if yn > y2t-tol && yx < y2oe+tol && (in_lgap||in_rgap)
                        if !(ih && abs(ym-gap)<tol); push!(sa_curves, tag); continue end
                    end
                end
            else
                # Single chip thick, no oxide
                in_trace_x = xn > x_tl-mr && xx < x_tr+mr
                in_rgap = xn > x_tr-mr && xx < x_gr+mr
                in_lgap = xn > x_gl-mr && xx < x_tl+mr
                if yn > -tol && yx < t_metal+tol && in_trace_x
                    push!(pec_trace, tag); continue end
                if ih && abs(ym)<tol && (in_rgap||in_lgap)
                    if has_overetch
                        # y=0 in gap is vacuum-vacuum interface, not SA
                    else
                        push!(sa_curves, tag)
                    end
                    continue
                end
                if yn > -tol && yx < t_metal+tol && !in_trace_x
                    push!(pec_ground, tag); continue end
                if has_overetch
                    if yn > -oe-tol && yx < tol && (in_rgap||in_lgap)
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
        gmsh.model.addPhysicalGroup(dom_dim, l1_doms, di, "L1_substrate"); di+=1
        gmsh.model.addPhysicalGroup(dom_dim, l2_doms, di, "L2_substrate"); di+=1
    else
        gmsh.model.addPhysicalGroup(dom_dim, l1_doms, di, "substrate"); di+=1
    end
    gmsh.model.addPhysicalGroup(dom_dim, gap_doms, di, "vacuum"); di+=1

    if has_oxide && thick
        if flipchip
            for (d, nm) in [(l1_ms_ox,"L1_MS_oxide"),(l1_ma_ox,"L1_MA_oxide"),
                            (l2_ms_ox,"L2_MS_oxide"),(l2_ma_ox,"L2_MA_oxide"),
                            (l1_sa_ox,"L1_SA_oxide"),(l2_sa_ox,"L2_SA_oxide")]
                !isempty(d) && (gmsh.model.addPhysicalGroup(dom_dim, d, di, nm); di+=1)
            end
        else
            for (d, nm) in [(sc_ms_ox,"MS_oxide"),(sc_ma_ox,"MA_oxide"),
                            (sc_sa_ox,"SA_oxide")]
                !isempty(d) && (gmsh.model.addPhysicalGroup(dom_dim, d, di, nm); di+=1)
            end
        end
    end

    bi = 1
    !isempty(pec_trace) && (gmsh.model.addPhysicalGroup(bdr_dim,pec_trace,bi,"PEC_trace"); bi+=1)
    !isempty(pec_ground) && (gmsh.model.addPhysicalGroup(bdr_dim,pec_ground,bi,"PEC_ground"); bi+=1)
    !isempty(bb_pec) && (gmsh.model.addPhysicalGroup(bdr_dim,bb_pec,bi,"bump_bonds"); bi+=1)
    !isempty(sa_curves) && (gmsh.model.addPhysicalGroup(bdr_dim,sa_curves,bi,"SA"); bi+=1)
    !isempty(outer_c) && (gmsh.model.addPhysicalGroup(bdr_dim,outer_c,bi,"outer"); bi+=1)
    !isempty(port_front) && (gmsh.model.addPhysicalGroup(bdr_dim,port_front,bi,"port_front"); bi+=1)
    !isempty(port_back) && (gmsh.model.addPhysicalGroup(bdr_dim,port_back,bi,"port_back"); bi+=1)

    # =====================================================================
    # Mesh size control
    # =====================================================================
    # Metal edge locations (x, y). In 2D we use points at (x,y,0); in 3D we
    # use vertical lines spanning the extrusion so the distance field is
    # z-invariant and preserves the cross-section resolution at every z.
    edge_pts_xy = Tuple{Float64,Float64}[]
    if thick
        if flipchip
            append!(edge_pts_xy, [(x_tl+tp,y2b),(x_tr-tp,y2b),(x_gl-tp,y2b),(x_gr+tp,y2b),
                                  (x_cl-tp,t_metal),(x_cr+tp,t_metal)])
        else
            append!(edge_pts_xy, [(x_tl+tp,t_metal),(x_tr-tp,t_metal),
                                  (x_gl-tp,t_metal),(x_gr+tp,t_metal)])
        end
    else
        if flipchip
            append!(edge_pts_xy, [(x_tl,gap),(x_tr,gap),(x_gl,gap),(x_gr,gap),
                                  (x_cl,0.0),(x_cr,0.0)])
        else
            for x in [x_tl, x_tr, x_gl, x_gr]
                push!(edge_pts_xy, (x, 0.0))
            end
        end
    end

    pts = Float64[]      # Distance point list
    crvs = Float64[]     # Distance curve list (unstructured 3D only)
    use_curves = extruded && extrude_nz == 0
    if use_curves
        # Unstructured 3D: span the extrusion with vertical lines so the size
        # field is z-invariant.
        for (x, y) in edge_pts_xy
            p0 = K.addPoint(x, y, 0.0)
            p1 = K.addPoint(x, y, extrude_length)
            push!(crvs, Float64(K.addLine(p0, p1)))
        end
    else
        # 2D or structured 3D: the 2D cross-section mesh at z=0 determines
        # everything (structured extrusion replicates it). Points at z=0 only
        # guarantees an L-independent xy mesh.
        for (x, y) in edge_pts_xy
            push!(pts, Float64(K.addPoint(x, y, 0.0)))
        end
    end
    K.synchronize()

    field_id = 1
    gmsh.model.mesh.field.add("Distance", field_id)
    if use_curves
        gmsh.model.mesh.field.setNumbers(field_id, "CurvesList", crvs)
    else
        gmsh.model.mesh.field.setNumbers(field_id, "PointsList", pts)
    end
    field_id += 1
    gmsh.model.mesh.field.add("Threshold", field_id)
    gmsh.model.mesh.field.setNumber(field_id, "InField", field_id-1)
    gmsh.model.mesh.field.setNumber(field_id, "SizeMin", lc_fine)
    gmsh.model.mesh.field.setNumber(field_id, "SizeMax", lc_far)
    dist_scale = flipchip ? gap : w_gap
    gmsh.model.mesh.field.setNumber(field_id, "DistMin", thick ? 5.0*max(tox,t_metal/50) : 0.1)
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

    if extruded
        gmsh.model.mesh.generate(3)
    else
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Netgen")
    end
    gmsh.model.mesh.setOrder(mesh_order)
    if mesh_order > 1 && !extruded
        gmsh.model.mesh.optimize("HighOrderElastic")
    end

    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    gmsh.option.setNumber("Mesh.Binary",1)
    gmsh.write(joinpath(@__DIR__,filename))

    # Summary
    println("=== CPW $(extruded ? "Quasi-2D (3D)" : "2D") Cross-Section Mesh ===")
    println("  Topology: $(flipchip ? "flip-chip" : "single chip")")
    println("  Metal: $(thick ? "thick ($(t_metal*1000)nm)" : "thin")")
    if thick
        println("  Sidewall: $(sidewall_angle)° r_top=$(r_top*1000)nm r_bot=$(r_bot*1000)nm")
        has_overetch && println("  Overetch: $(oe*1000)nm")
    end
    has_oxide && println("  Oxide: $(tox*1000)nm")
    println("  w_trace=$w_trace w_gap=$w_gap μm")
    flipchip && println("  cow=$cow gap=$fc_gap μm")
    extruded && println("  Extrusion: $extrude_length μm in +z")
    dom_ent = extruded ? "vols" : "surfs"
    bdr_ent = extruded ? "faces" : "curves"
    for (dim, tag) in gmsh.model.getPhysicalGroups(dom_dim)
        name = gmsh.model.getPhysicalName(dim, tag)
        ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  Domain $tag: $name ($(length(ents)) $dom_ent)")
    end
    for (dim, tag) in gmsh.model.getPhysicalGroups(bdr_dim)
        name = gmsh.model.getPhysicalName(dim, tag)
        ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        println("  Boundary $tag: $name ($(length(ents)) $bdr_ent)")
    end
    println("  Nodes: $(length(gmsh.model.mesh.getNodes()[1]))")
    println("  File: $filename")

    gmsh.finalize()
end
