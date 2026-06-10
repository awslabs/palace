# Test EGADS sensitivity for an extruded LINE body via manual geometry_dot.
#
# extrude_dot! fails for straight-edge extrusions because OCC simplifies
# lateral EXTRUSION surfaces to PLANE. The workaround: set geometry_dot
# directly on each face, edge, and node of the extruded body.
#
# Test case: unit square [0,L]×[0,L] in XY plane, extruded along Z by H.
# Parameters: side length L and extrude distance H.
#
# Key insight: for PLANE surfaces, P(u,v) = origin + u*udir + v*vdir.
# The sensitivity is: dP/dp = dorigin/dp + u*(dudir/dp) + v*(dvdir/dp).
# Since udir and vdir are fixed-orientation unit vectors for a rigid plane,
# only dorigin/dp is nonzero for a translating plane.
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev test_extrude_dot_line.jl

using LinearAlgebra
using Printf
using EGADS

# ─── Build a square FACEBODY from 4 line segments ─────────────────────────

function make_square_face_body(ctx, L)
    c1 = [0.0, 0.0, 0.0]
    c2 = [L,   0.0, 0.0]
    c3 = [L,   L,   0.0]
    c4 = [0.0, L,   0.0]

    n1 = EGADS.make_node(ctx, c1...)
    n2 = EGADS.make_node(ctx, c2...)
    n3 = EGADS.make_node(ctx, c3...)
    n4 = EGADS.make_node(ctx, c4...)

    e1 = EGADS.make_edge_from_nodes(ctx, n1, n2)
    e2 = EGADS.make_edge_from_nodes(ctx, n2, n3)
    e3 = EGADS.make_edge_from_nodes(ctx, n3, n4)
    e4 = EGADS.make_edge_from_nodes(ctx, n4, n1)

    loop = EGADS.make_wire_loop(ctx, [e1, e2, e3, e4])
    face = EGADS.make_plane_face(ctx, loop)
    body = EGADS.make_face_body(ctx, face)

    return body
end

# ─── Analytical sensitivity for the extruded square prism ─────────────────
#
# For parameter L (side length scaling), a point at (x, y, z) on the prism
# has sensitivity dx/dL = (x/L, y/L, 0) — uniform radial scaling in XY.
#
# For parameter H (extrude distance), a point at (x, y, z) on the prism
# has sensitivity dx/dH = (0, 0, z/H) — linear interpolation in Z.

function analytical_dxdp_dL(x, y, z, L)
    return [x / L, y / L, 0.0]
end

function analytical_dxdp_dH(x, y, z, H)
    return [0.0, 0.0, z / H]
end

# ─── Set geometry_dot on all sub-topologies of the extruded body ──────────
#
# This is the manual alternative to extrude_dot!, which fails when OCC
# simplifies surfaces. We set dots on every node, edge, and face geometry
# of the extruded body directly.
#
# The key is computing the correct derivative for each geometry object:
# - Nodes: dx_dot = analytical sensitivity at the node position
# - Lines: origin_dot + direction_dot from endpoint derivatives
# - Planes: origin_dot (plane orientation is fixed for pure scaling/translation)

function set_body_dots_manual!(body, analytical_fn)
    # Set all node sensitivities first
    nodes = EGADS.get_body_topos(body, EGADS.OClass.NODE)
    for n in nodes
        nt = EGADS.get_topology(n)
        pos = Float64[nt.limits[1:3]...]
        pos_dot = Float64[analytical_fn(pos...)...]
        EGADS.set_geometry_dot!(n, EGADS.OClass.NODE, Int32(0),
                                nothing, pos, pos_dot)
    end

    # Set all edge sensitivities
    edges = EGADS.get_body_topos(body, EGADS.OClass.EDGE)
    for e in edges
        et = EGADS.get_topology(e)
        isnothing(et.geom) && continue

        geom_info = EGADS.get_geometry(et.geom)
        rvec = Float64[geom_info.rvec...]

        if geom_info.mtype == EGADS.MType.LINE
            # LINE: rvec = [origin(3), direction(3)]
            origin = rvec[1:3]
            direction = rvec[4:6]
            endpoint = origin + direction

            origin_dot = Float64[analytical_fn(origin...)...]
            endpoint_dot = Float64[analytical_fn(endpoint...)...]
            direction_dot = endpoint_dot - origin_dot

            rvec_dot = Float64[origin_dot..., direction_dot...]
        else
            # For non-line curves, set zero dot (conservative fallback)
            rvec_dot = zeros(length(rvec))
        end

        EGADS.set_geometry_dot!(et.geom, EGADS.OClass.CURVE, geom_info.mtype,
                                nothing, rvec, rvec_dot)

        # Edge t-range: for distance-parameterized lines, t ∈ [0, |dir|]
        range_info = EGADS.get_range(e)
        t0, t1 = range_info.range[1:2]

        if geom_info.mtype == EGADS.MType.LINE
            direction = rvec[4:6]
            direction_dot = rvec_dot[4:6]
            edge_len = norm(direction)
            edge_len_dot = edge_len > 0 ? dot(direction, direction_dot) / edge_len : 0.0
            EGADS.set_range_dot!(e, EGADS.OClass.EDGE,
                                 Float64[t0, t1], Float64[0.0, edge_len_dot])
        else
            EGADS.set_range_dot!(e, EGADS.OClass.EDGE,
                                 Float64[t0, t1], Float64[0.0, 0.0])
        end
    end

    # Set all face surface sensitivities
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    for f in faces
        ft = EGADS.get_topology(f)
        isnothing(ft.geom) && continue

        geom_info = EGADS.get_geometry(ft.geom)
        rvec = Float64[geom_info.rvec...]

        if geom_info.mtype == EGADS.MType.PLANE
            # PLANE: rvec = [origin(3), udir(3), vdir(3)]
            # P(u,v) = origin + u*udir + v*vdir
            # For a linear deformation f(x) = Ax:
            #   dP/dp = A*origin + u*(A*udir) + v*(A*vdir)
            # So: origin_dot = f(origin), udir_dot = f(udir_endpoint) - f(origin),
            #     vdir_dot = f(vdir_endpoint) - f(origin)
            # But since f is linear: udir_dot = A*udir = f(origin+udir) - f(origin)
            origin = rvec[1:3]
            udir = rvec[4:6]
            vdir = rvec[7:9]

            f_origin = Float64[analytical_fn(origin...)...]
            f_udir_end = Float64[analytical_fn((origin .+ udir)...)...]
            f_vdir_end = Float64[analytical_fn((origin .+ vdir)...)...]

            origin_dot = f_origin
            udir_dot = f_udir_end .- f_origin
            vdir_dot = f_vdir_end .- f_origin

            rvec_dot = Float64[origin_dot..., udir_dot..., vdir_dot...]
        else
            rvec_dot = zeros(length(rvec))
        end

        EGADS.set_geometry_dot!(ft.geom, EGADS.OClass.SURFACE, geom_info.mtype,
                                nothing, rvec, rvec_dot)
    end

    return nothing
end

# ─── Main test ────────────────────────────────────────────────────────────

function main()
    println("═══ EGADS Manual Sensitivity Test (extruded square prism) ═══")
    println()

    ctx = EGADS.EGADSContext()

    L = 1.0
    H = 2.0
    dir = [0.0, 0.0, 1.0]

    println("Building square face body (L=$L)...")
    src = make_square_face_body(ctx, L)
    println("  Source body created.")

    println("Extruding along Z by H=$H...")
    body = EGADS.extrude(src, H, dir)
    println("  Extruded body created.")

    tess = EGADS.tessellate(body; max_edge=0.2, sag=0.01, angle=15.0)
    n_tess = EGADS.tess_status(tess).npts
    println("  Tessellation: $n_tess points")

    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    println("  Faces: $(length(faces))")
    println()

    # ═══ Test 1: d/dL (side length sensitivity) ═══════════════════════════

    println("─── Test 1: d/dL (side length sensitivity) ───")

    set_body_dots_manual!(body, (x, y, z) -> analytical_dxdp_dL(x, y, z, L))

    max_err = 0.0
    n_pts = 0

    for fi in eachindex(faces)
        fmesh = EGADS.get_face_mesh(tess, body, fi)
        face_max_err = 0.0
        for j in 1:size(fmesh.xyz, 1)
            u, v = fmesh.uv[j, 1], fmesh.uv[j, 2]
            res = EGADS.evaluate_dot(faces[fi], [u, v], [0.0, 0.0])
            dx_dp = res.result_dot.xyz
            x, y, z = fmesh.xyz[j, :]

            analytical = analytical_dxdp_dL(x, y, z, L)
            err = norm(dx_dp - analytical)
            face_max_err = max(face_max_err, err)
            max_err = max(max_err, err)
            n_pts += 1
        end
        @printf("  Face %d: max error = %.2e\n", fi, face_max_err)
    end

    println("  Overall ($n_pts pts): max error = $max_err")
    pass1 = max_err < 1e-10
    println("  Test 1 (d/dL): $(pass1 ? "PASS ✓" : "FAIL ✗")")
    println()

    # ═══ Test 2: d/dH (extrude distance sensitivity) ═════════════════════

    println("─── Test 2: d/dH (extrude distance sensitivity) ───")

    set_body_dots_manual!(body, (x, y, z) -> analytical_dxdp_dH(x, y, z, H))

    max_err_dH = 0.0
    n_dH = 0

    for fi in eachindex(faces)
        fmesh = EGADS.get_face_mesh(tess, body, fi)
        face_max_err = 0.0
        for j in 1:size(fmesh.xyz, 1)
            u, v = fmesh.uv[j, 1], fmesh.uv[j, 2]
            res = EGADS.evaluate_dot(faces[fi], [u, v], [0.0, 0.0])
            dx_dp = res.result_dot.xyz
            x, y, z = fmesh.xyz[j, :]

            analytical = analytical_dxdp_dH(x, y, z, H)
            err = norm(dx_dp - analytical)
            face_max_err = max(face_max_err, err)
            max_err_dH = max(max_err_dH, err)
            n_dH += 1
        end
        @printf("  Face %d: max error = %.2e\n", fi, face_max_err)
    end

    println("  Overall ($n_dH pts): max error = $max_err_dH")
    pass2 = max_err_dH < 1e-10
    println("  Test 2 (d/dH): $(pass2 ? "PASS ✓" : "FAIL ✗")")
    println()

    # ═══ Test 3: FD validation (d/dL) ═════════════════════════════════════

    println("─── Test 3: Finite difference validation (d/dL) ───")

    h_fd = 1e-6
    L_pert = L + h_fd
    src_pert = make_square_face_body(ctx, L_pert)
    body_pert = EGADS.extrude(src_pert, H, dir)
    tess_pert = EGADS.tessellate(body_pert; max_edge=0.2, sag=0.01, angle=15.0)

    props0 = EGADS.tess_mass_props(tess)
    props1 = EGADS.tess_mass_props(tess_pert)

    vol0 = props0[1]
    vol1 = props1[1]
    dvol_dL_fd = (vol1 - vol0) / h_fd

    # Analytical: V = L² * H, dV/dL = 2*L*H
    dvol_dL_analytical = 2 * L * H

    fd_err = abs(dvol_dL_fd - dvol_dL_analytical) / abs(dvol_dL_analytical)
    @printf("  V₀ = %.6f, V₁ = %.6f\n", vol0, vol1)
    @printf("  dV/dL (FD):         %.10f\n", dvol_dL_fd)
    @printf("  dV/dL (analytical): %.10f\n", dvol_dL_analytical)
    @printf("  Relative error:     %.2e\n", fd_err)

    pass3 = fd_err < 1e-4
    println("  Test 3 (FD vol): $(pass3 ? "PASS ✓" : "FAIL ✗")")
    println()

    # ═══ Summary ══════════════════════════════════════════════════════════

    println("═══ Summary ═══")
    println("  Test 1 (d/dL evaluate_dot): $(pass1 ? "PASS ✓" : "FAIL ✗")")
    println("  Test 2 (d/dH evaluate_dot): $(pass2 ? "PASS ✓" : "FAIL ✗")")
    println("  Test 3 (FD volume):         $(pass3 ? "PASS ✓" : "FAIL ✗")")

    EGADS.close(ctx)

    all_pass = pass1 && pass2 && pass3
    all_pass || error("extrude_dot test FAILED")
    println("\n  All tests passed.")
end

main()
