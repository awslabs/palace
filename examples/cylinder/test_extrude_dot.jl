# Test EGADS extrude_dot! following the calling convention from egadsHLevel_dot.c.
#
# This test constructs a circle face body from low-level EGADS primitives,
# sets geometry derivatives on ALL sub-objects (circle curve, node, edge range,
# plane surface), extrudes it, then calls extrude_dot! and verifies evaluate_dot
# returns the correct dx/dp values.
#
# The key requirement (documented in the C examples) is that the SOURCE body
# must have geometry derivatives set on its sub-geometries BEFORE calling
# extrude_dot!. Without this, EG_extrude_dot returns EGADS_NODATA.
#
# Test case: circle of radius R in XY plane, extruded along Z by distance H.
# Two parameters tested:
#   d/dR: analytical dx/dp = (x/R, y/R, 0) on the lateral surface
#   d/dH: analytical dx/dp on top cap = (0, 0, 1) scaled by extrude_dot
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev test_extrude_dot.jl

using LinearAlgebra
using Printf
using EGADS

# ─── Build a circle FACEBODY (equivalent to makeCircleBody in C) ─────────

function make_circle_face_body(ctx, xcent, xax, yax, radius)
    # 1. Create circle curve geometry: rvec = [center(3), xaxis(3), yaxis(3), radius]
    circle = EGADS.make_circle(ctx, xcent, xax, yax, radius)

    # 2. Create edge on the circle (full 2π, periodic → ONENODE)
    edge = EGADS.make_edge(ctx, circle, 0.0, 2π)

    # 3. Create loop from the edge
    loop = EGADS.make_wire_loop(ctx, [edge])

    # 4. Create planar face from the loop
    face = EGADS.make_plane_face(ctx, loop)

    # 5. Wrap in a FACEBODY
    body = EGADS.make_face_body(ctx, face)

    return body
end

# ─── Set derivatives on circle FACEBODY (equivalent to setCircleBody_dot) ──

function set_circle_body_dot!(body, xcent, xcent_dot, xax, xax_dot, yax, yax_dot,
                               radius, radius_dot)
    # Traverse topology: FACEBODY → FACE → LOOP → EDGE → (NODE, CIRCLE)
    body_topo = EGADS.get_topology(body)
    face = body_topo.children[1]

    face_topo = EGADS.get_topology(face)
    plane_geom = face_topo.geom
    loop = face_topo.children[1]

    loop_topo = EGADS.get_topology(loop)
    edge = loop_topo.children[1]

    edge_topo = EGADS.get_topology(edge)
    circle_geom = edge_topo.geom
    node = edge_topo.children[1]

    # (a) Set edge t-range sensitivity (range [0, 2π] doesn't change)
    EGADS.set_range_dot!(edge, EGADS.OClass.EDGE,
                         Float64[0.0, 2π], Float64[0.0, 0.0])

    # (b) Set circle geometry sensitivity
    #     rvec = [center(3), xaxis(3), yaxis(3), radius]
    circle_data = Float64[xcent..., xax..., yax..., radius]
    circle_data_dot = Float64[xcent_dot..., xax_dot..., yax_dot..., radius_dot]
    EGADS.set_geometry_dot!(circle_geom, EGADS.OClass.CURVE, EGADS.MType.CIRCLE,
                            nothing, circle_data, circle_data_dot)

    # (c) Read back the axes (EGADS may orthonormalize)
    dot_info = EGADS.get_geometry_dot(circle_geom)
    rvec = dot_info.reals
    rvec_dot = dot_info.reals_dot
    dx = rvec[4:6]
    dy = rvec[7:9]
    dx_dot = rvec_dot[4:6]
    dy_dot = rvec_dot[7:9]

    # (d) Set node sensitivity: node position = center + dx * radius
    node_pos = xcent .+ dx .* radius
    node_pos_dot = xcent_dot .+ dx_dot .* radius .+ dx .* radius_dot
    EGADS.set_geometry_dot!(node, EGADS.OClass.NODE, Int32(0),
                            nothing, node_pos, node_pos_dot)

    # (e) Set plane geometry sensitivity
    #     rvec = [center(3), udir(3), vdir(3)]
    plane_data = Float64[xcent..., dx..., dy...]
    plane_data_dot = Float64[xcent_dot..., dx_dot..., dy_dot...]
    EGADS.set_geometry_dot!(plane_geom, EGADS.OClass.SURFACE, EGADS.MType.PLANE,
                            nothing, plane_data, plane_data_dot)

    return nothing
end

# ─── Main test ────────────────────────────────────────────────────────────

function main()
    println("═══ EGADS extrude_dot! Test ═══")
    println()

    ctx = EGADS.EGADSContext()

    R = 1.0
    H = 2.0
    xcent = [0.0, 0.0, 0.0]
    xax   = [1.0, 0.0, 0.0]
    yax   = [0.0, 1.0, 0.0]
    dir   = [0.0, 0.0, 1.0]

    # ─── Build circle face body ───────────────────────────────────────────

    println("Building circle face body (R=$R)...")
    src = make_circle_face_body(ctx, xcent, xax, yax, R)
    println("  Source body created.")

    # ─── Extrude ──────────────────────────────────────────────────────────

    println("Extruding along Z by H=$H...")
    body = EGADS.extrude(src, H, dir)
    println("  Extruded body created.")

    # ─── Tessellate for evaluation points ─────────────────────────────────

    tess = EGADS.tessellate(body; max_edge=0.2, sag=0.01, angle=15.0)
    n_tess = EGADS.tess_status(tess).npts
    println("  Tessellation: $n_tess points")

    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    println("  Faces: $(length(faces))")

    # ═══ Test 1: d/dR (radius sensitivity) ═══════════════════════════════

    println()
    println("─── Test 1: d/dR (radius sensitivity) ───")

    # Set source derivatives for d/dR:
    #   circle radius_dot = 1, everything else = 0
    set_circle_body_dot!(src, xcent, zeros(3), xax, zeros(3), yax, zeros(3), R, 1.0)

    # Verify source has geometry dots
    src_topo = EGADS.get_topology(src)
    src_face = src_topo.children[1]
    src_face_topo = EGADS.get_topology(src_face)
    src_loop = src_face_topo.children[1]
    src_loop_topo = EGADS.get_topology(src_loop)
    src_edge = src_loop_topo.children[1]
    src_edge_topo = EGADS.get_topology(src_edge)
    circle_geom = src_edge_topo.geom
    has_dot = EGADS.has_geometry_dot(circle_geom)
    println("  Source circle has_geometry_dot: $has_dot")

    # Call extrude_dot! with dist_dot=0, dir_dot=0 (only radius changes)
    println("  Calling extrude_dot!...")
    EGADS.extrude_dot!(body, src, H, 0.0, dir, zeros(3))
    println("  extrude_dot! succeeded!")

    # Verify body has geometry dots
    has_dot_body = EGADS.has_geometry_dot(body)
    println("  Extruded body has_geometry_dot: $has_dot_body")

    # Check evaluate_dot on each face
    max_err_lateral = 0.0
    n_lateral = 0
    max_err_cap = 0.0
    n_cap = 0

    for fi in eachindex(faces)
        fmesh = EGADS.get_face_mesh(tess, body, fi)
        topo = EGADS.get_topology(faces[fi])
        geom = EGADS.get_geometry(topo.geom)
        is_curved = geom.mtype != EGADS.MType.PLANE

        for j in 1:size(fmesh.xyz, 1)
            u, v = fmesh.uv[j, 1], fmesh.uv[j, 2]
            res = EGADS.evaluate_dot(faces[fi], [u, v], [0.0, 0.0])
            dx_dp = res.result_dot.xyz
            x, y, z = fmesh.xyz[j, :]

            if is_curved
                # Lateral surface: dx/dR = (x/R, y/R, 0)
                analytical = [x / R, y / R, 0.0]
                err = norm(dx_dp - analytical)
                max_err_lateral = max(max_err_lateral, err)
                n_lateral += 1
            else
                # Cap faces: dx/dR should be nonzero (caps scale with R)
                # Top cap: dx/dR = (x/R, y/R, 0), Bottom cap: same
                analytical = [x / R, y / R, 0.0]
                err = norm(dx_dp - analytical)
                max_err_cap = max(max_err_cap, err)
                n_cap += 1
            end
        end
    end

    println("  Lateral surface ($n_lateral pts): max error = $max_err_lateral")
    println("  Cap faces ($n_cap pts): max error = $max_err_cap")

    pass1 = max_err_lateral < 1e-10 && max_err_cap < 1e-10
    println("  Test 1 (d/dR): $(pass1 ? "PASS ✓" : "FAIL ✗")")

    # ═══ Test 2: d/dH (extrude distance sensitivity) ═════════════════════

    println()
    println("─── Test 2: d/dH (extrude distance sensitivity) ───")

    # Zero out source dots (radius is not changing for d/dH)
    EGADS.zero_geometry_dot!(src)
    set_circle_body_dot!(src, xcent, zeros(3), xax, zeros(3), yax, zeros(3), R, 0.0)

    # Call extrude_dot! with dist_dot=1 (H is changing)
    println("  Calling extrude_dot!...")
    EGADS.extrude_dot!(body, src, H, 1.0, dir, zeros(3))
    println("  extrude_dot! succeeded!")

    # For d/dH: the extrusion direction is Z, so dx/dH = (0, 0, z/H)
    # Points at z=0 don't move, points at z=H move by (0,0,1).
    # Interior points move linearly: (0, 0, z/H).
    max_err_dH = 0.0
    n_dH = 0

    for fi in eachindex(faces)
        fmesh = EGADS.get_face_mesh(tess, body, fi)
        for j in 1:size(fmesh.xyz, 1)
            u, v = fmesh.uv[j, 1], fmesh.uv[j, 2]
            res = EGADS.evaluate_dot(faces[fi], [u, v], [0.0, 0.0])
            dx_dp = res.result_dot.xyz
            x, y, z = fmesh.xyz[j, :]

            # Analytical: dx/dH = (0, 0, z/H)
            analytical = [0.0, 0.0, z / H]
            err = norm(dx_dp - analytical)
            max_err_dH = max(max_err_dH, err)
            n_dH += 1
        end
    end

    println("  All faces ($n_dH pts): max error = $max_err_dH")

    pass2 = max_err_dH < 1e-10
    println("  Test 2 (d/dH): $(pass2 ? "PASS ✓" : "FAIL ✗")")

    # ═══ Test 3: FD validation ════════════════════════════════════════════

    println()
    println("─── Test 3: Finite difference validation (d/dR) ───")

    h_fd = 1e-6
    R_pert = R + h_fd
    src_pert = make_circle_face_body(ctx, xcent, xax, yax, R_pert)
    body_pert = EGADS.extrude(src_pert, H, dir)
    tess_pert = EGADS.tessellate(body_pert; max_edge=0.2, sag=0.01, angle=15.0)

    # Use mapTessBody to match tessellation points between original and perturbed
    # For now, just compare tess mass properties as a simpler FD check
    props0 = EGADS.tess_mass_props(tess)
    props1 = EGADS.tess_mass_props(tess_pert)

    vol0 = props0[1]
    vol1 = props1[1]
    dvol_dR_fd = (vol1 - vol0) / h_fd

    # Analytical: V = π R² H, dV/dR = 2π R H
    dvol_dR_analytical = 2π * R * H

    fd_err = abs(dvol_dR_fd - dvol_dR_analytical) / abs(dvol_dR_analytical)
    @printf("  V₀ = %.6f, V₁ = %.6f\n", vol0, vol1)
    @printf("  dV/dR (FD):         %.10f\n", dvol_dR_fd)
    @printf("  dV/dR (analytical): %.10f\n", dvol_dR_analytical)
    @printf("  Relative error:     %.2e\n", fd_err)

    pass3 = fd_err < 1e-4
    println("  Test 3 (FD vol): $(pass3 ? "PASS ✓" : "FAIL ✗")")

    # ═══ Summary ══════════════════════════════════════════════════════════

    println()
    println("═══ Summary ═══")
    println("  Test 1 (d/dR evaluate_dot): $(pass1 ? "PASS ✓" : "FAIL ✗")")
    println("  Test 2 (d/dH evaluate_dot): $(pass2 ? "PASS ✓" : "FAIL ✗")")
    println("  Test 3 (FD volume):         $(pass3 ? "PASS ✓" : "FAIL ✗")")

    EGADS.close(ctx)

    all_pass = pass1 && pass2 && pass3
    all_pass || error("extrude_dot test FAILED")
    println("\n  All tests passed.")
end

main()
