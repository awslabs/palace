# Full EGADS sensitivity chain: EGADS surface dx/dp + Laplacian → df/dp.
#
# For each design parameter p, the chain computes df/dp via:
#   1. Set EGADS body derivatives for parameter p (make_solid_body_dot!)
#   2. Classify faces by evaluate_dot magnitude (multi-point sampling)
#   3. For each Palace boundary vertex on a "moving" face: get dx/dp from EGADS
#   4. Propagate dx/dp from boundary to interior via Laplacian
#   5. Contract: df/dp = Σᵢ (df/dxᵢ) · (dxᵢ/dp)
#
# Tests both R (radius) and H (height) parameters on a cylinder.
#
# Prerequisites:
#   - Palace must have been run with MeshSensitivity=true producing .mesh + .solb
#   - EGADS.jl must be available (via DeviceLayout.jl project or standalone)
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev \
#          examples/cylinder/chain_test_egads.jl [postpro_dir] [radius]

using LinearAlgebra
using SparseArrays
using NearestNeighbors: KDTree, nn

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))

using EGADS

# ─── Per-parameter face classification ───────────────────────────────────────

"""
    classify_sensitivity_faces(faces, n_samples=5) -> BitVector

Classify which EGADS faces have nonzero evaluate_dot for the currently-set
parameter. Samples evaluate_dot at multiple UV points per face (center + 4
corners of the UV range). A face is "moving" if ANY sample has ‖dx/dp‖ > tol.

This replaces the previous geometry-type-based classification (curved vs planar)
which is only correct for d/dR on a primitive cylinder.
"""
function classify_sensitivity_faces(faces; tol=1e-14)
    face_is_moving = falses(length(faces))
    for (fi, f) in enumerate(faces)
        range_info = EGADS.get_range(f)
        u0, u1, v0, v1 = range_info.range[1:4]
        um = (u0 + u1) / 2
        vm = (v0 + v1) / 2
        # Sample at 5 points: center + 4 corners (offset slightly from edges)
        du = (u1 - u0) * 0.1
        dv = (v1 - v0) * 0.1
        sample_uvs = [
            (um, vm),
            (u0 + du, v0 + dv),
            (u1 - du, v0 + dv),
            (u0 + du, v1 - dv),
            (u1 - du, v1 - dv),
        ]
        for (u, v) in sample_uvs
            res = EGADS.evaluate_dot(f, [u, v], [0.0, 0.0])
            if norm(res.result_dot.xyz) > tol
                face_is_moving[fi] = true
                break
            end
        end
    end
    return face_is_moving
end

# ─── Generic EGADS→Palace chain for one parameter ───────────────────────────

"""
    compute_egads_chain(body, faces, mesh, dfdx, is_hex_boundary;
                        proj_tol=1e-8) -> (dfdp, dxdp_full)

Compute df/dp for the currently-set parameter derivatives on `body`.

1. Classify faces by evaluate_dot magnitude (multi-point sampling)
2. For each Palace boundary vertex on a "moving" face: get dx/dp via inv_evaluate + evaluate_dot
3. Propagate from boundary to interior via Laplacian
4. Contract: df/dp = Σᵢ (df/dxᵢ) · (dxᵢ/dp)

Returns the scalar df/dp and the full Nx3 dxdp matrix (after Laplacian propagation).
"""
function compute_egads_chain(body, faces, mesh, dfdx, is_hex_boundary; proj_tol=1e-8)
    nv = size(mesh.vertices, 1)
    face_is_moving = classify_sensitivity_faces(faces)
    n_moving = count(face_is_moving)
    n_static = length(faces) - n_moving

    boundary_indices = findall(is_hex_boundary)
    egads_dxdp = zeros(Float64, nv, 3)
    is_sensitivity_boundary = falses(nv)
    n_on_moving = 0
    n_on_static = 0

    for vi in boundary_indices
        pt = mesh.vertices[vi, :]  # Both EGADS body and mesh are in cm
        found = false
        for fi in eachindex(faces)
            face_is_moving[fi] || continue
            uv, result = EGADS.inv_evaluate(faces[fi], pt)
            dist = norm(result[1:3] - pt)
            if dist < proj_tol
                res = EGADS.evaluate_dot(faces[fi], uv, [0.0, 0.0])
                egads_dxdp[vi, :] .= res.result_dot.xyz
                is_sensitivity_boundary[vi] = true
                found = true
                n_on_moving += 1
                break
            end
        end
        found || (n_on_static += 1)
    end

    propagated = propagate_boundary_sensitivity(mesh, egads_dxdp, is_sensitivity_boundary)
    dfdp = dot(vec(dfdx), vec(propagated))

    return (dfdp=dfdp, dxdp=propagated, n_moving=n_moving, n_static=n_static,
            n_on_moving=n_on_moving, n_on_static=n_on_static)
end

# ─── Main test ───────────────────────────────────────────────────────────────

function main_egads()
    postpro_dir =
        length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "postpro", "cavity_pec_sens_p1_ref")
    radius = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2.74  # cm

    mesh_file = joinpath(postpro_dir, "mesh_for_sensitivity.mesh")
    solb_file = joinpath(postpro_dir, "mesh_sensitivity.solb")
    eig_file = joinpath(postpro_dir, "eig.csv")

    println("═══ EGADS Multi-Parameter Sensitivity Chain Test ═══")
    println("  Postpro: $postpro_dir")
    println("  Radius:  $radius cm")
    println()

    # Read Palace output
    mesh = read_mesh(mesh_file)
    sol = read_solb(solb_file)
    freq = read_eigenfrequency(eig_file; mode=1)
    nv = size(mesh.vertices, 1)
    dfdx = sol.data
    is_hex_boundary = classify_boundary_vertices(mesh)
    println("  Palace mesh: $nv vertices, freq = $freq GHz")

    # ─── Build EGADS cylinder (in cm to match mesh coordinates) ─────────
    R = radius       # cm
    H = 2 * R        # height = diameter (H = 2R for this test case)

    ctx = EGADSContext()
    body = EGADS.make_cylinder(ctx, [0.0, 0.0, 0.0], [0.0, 0.0, H], R)
    faces = EGADS.get_body_topos(body, OClass.FACE)
    data = Float64[0.0, 0.0, 0.0, 0.0, 0.0, H, R]

    # ═══ Parameter 1: Radius (R) ═════════════════════════════════════════

    println("═══ Parameter 1: Radius (d/dR) ═══")
    data_dot_R = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    EGADS.make_solid_body_dot!(body, SolidType.CYLINDER, data, data_dot_R)

    result_R = compute_egads_chain(body, faces, mesh, dfdx, is_hex_boundary)
    nondim_R = result_R.dfdp * radius / freq

    println("  Moving faces: $(result_R.n_moving), static: $(result_R.n_static)")
    println("  Vertices on moving: $(result_R.n_on_moving), on static: $(result_R.n_on_static)")
    println("  df/dR         = $(result_R.dfdp) GHz/cm")
    println("  df/dR * R / f = $nondim_R")
    println("  Expected:       -1.0")
    println("  Error:          $(abs(nondim_R + 1.0))")

    r_pass = abs(nondim_R + 1.0) < 5e-3
    println("  Status:         $(r_pass ? "PASS ✓" : "FAIL ✗")")
    println()

    # ═══ Parameter 2: Height (H) ═════════════════════════════════════════

    println("═══ Parameter 2: Height (d/dH) ═══")
    data_dot_H = Float64[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    EGADS.make_solid_body_dot!(body, SolidType.CYLINDER, data, data_dot_H)

    result_H = compute_egads_chain(body, faces, mesh, dfdx, is_hex_boundary)
    # For TM₀₁₀: f = c·x₀₁/(2πR√εᵣ), independent of H. So df/dH ≈ 0.
    # The chain should produce near-zero through correct mechanism:
    # moving faces have evaluate_dot ≠ 0, but the contraction cancels.
    nondim_H = result_H.dfdp * H / freq

    println("  Moving faces: $(result_H.n_moving), static: $(result_H.n_static)")
    println("  Vertices on moving: $(result_H.n_on_moving), on static: $(result_H.n_on_static)")
    println("  df/dH         = $(result_H.dfdp) GHz/cm")
    println("  df/dH * H / f = $nondim_H")
    println("  Expected:       ~0 (TM₀₁₀ is H-independent)")

    # For H: the gradient should be near zero but not exactly zero due to
    # discrete approximation. Accept |df/dH * H / f| < 0.01
    h_pass = abs(nondim_H) < 0.01
    println("  Status:         $(h_pass ? "PASS ✓" : "FAIL ✗")")
    println()

    # ═══ Summary ═════════════════════════════════════════════════════════

    println("═══ Summary ═══")
    println("  R chain (df/dR*R/f ≈ -1.0):  $(r_pass ? "PASS ✓" : "FAIL ✗")")
    println("  H chain (df/dH*H/f ≈  0.0):  $(h_pass ? "PASS ✓" : "FAIL ✗")")

    close(ctx)

    all_pass = r_pass && h_pass
    all_pass || error("EGADS chain test FAILED")
    return println("\n  All tests passed.")
end

main_egads()
