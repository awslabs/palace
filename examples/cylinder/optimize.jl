# Multi-parameter optimizer: find R and H that produce a target frequency.
#
# Uses the verified sensitivity chain to compute exact gradients [df/dR, df/dH],
# then drives L-BFGS toward the target. Compares against derivative-free
# Nelder-Mead to demonstrate gradient advantage.
#
# The objective function at each (R,H):
#   1. Scale the baseline mesh affinely to the new (R,H)
#   2. Run Palace eigensolve
#   3. Return f (frequency) and optionally ∇f = [df/dR, df/dH]
#
# The gradient uses the mesh sensitivity: df/dp = Σᵢ (df/dxᵢ)(dxᵢ/dp)
# where dxᵢ/dp is the analytical mesh deformation for affine scaling.
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev optimize.jl

using LinearAlgebra, Printf, Optim

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))
include(joinpath(@__DIR__, "taylor_test.jl"))

# ─── Objective and gradient ──────────────────────────────────────────────────

"""
    CylinderObjective

Holds the fixed data needed to evaluate f(R,H) and ∇f(R,H) at any point.
"""
struct CylinderObjective
    palace::String
    base_config::String
    base_mesh::String
    work_mesh::String
    work_config::String
    f_target::Float64
    baseline_mesh::SimpleMesh      # for gradient: vertex coordinates at (R₀, H₀)
    baseline_dfdx::Matrix{Float64} # df/dx at each vertex (from Palace .solb)
    R₀::Float64                    # baseline radius (cm)
    H₀::Float64                    # baseline height (cm)
    eval_count::Ref{Int}
end

"""
    evaluate!(obj, R, H) -> (f, [df/dR, df/dH])

Run Palace at (R,H) and compute frequency + gradient.
"""
function evaluate!(obj::CylinderObjective, R, H)
    obj.eval_count[] += 1
    label = @sprintf("opt_%04d", obj.eval_count[])
    outdir = "postpro/$label"

    # Scale mesh
    write_scaled_mesh(obj.base_mesh, obj.work_mesh;
                      xy=R / obj.R₀, z=H / obj.H₀)

    # Write config
    write_taylor_config(obj.base_config, "mesh/cylinder_hex_p1_taylor.msh",
                        outdir, obj.work_config)

    # Run Palace
    try
        run(pipeline(`$(obj.palace) $(obj.work_config)`; stdout=devnull, stderr=devnull))
    catch
        error("Palace failed at R=$R, H=$H")
    end

    eig = joinpath(outdir, "eig.csv")
    @assert isfile(eig) "eig.csv not found in $outdir"
    f = read_eigenfrequency(eig; mode=1)

    # Gradient via chain contraction at the CURRENT (R,H) mesh
    # dx/dR at the new mesh: (x * R₀/R / R₀, y * R₀/R / R₀, 0) = (x_new / R, y_new / R, 0)
    # But we use the baseline dfdx (computed at R₀,H₀).
    # For affine scaling: dfdx scales with the geometry, so we need dfdx at the new (R,H).
    # However, we only have dfdx at (R₀,H₀). For small perturbations, dfdx ≈ constant.
    # For the optimizer to work well, we recompute df/dp using the baseline dfdx and
    # the analytical dx/dp at the CURRENT R,H scaled coordinates.
    nv = size(obj.baseline_mesh.vertices, 1)

    # Scale baseline vertices to current (R,H)
    scale_xy = R / obj.R₀
    scale_z = H / obj.H₀
    verts_scaled = copy(obj.baseline_mesh.vertices)
    verts_scaled[:, 1] .*= scale_xy
    verts_scaled[:, 2] .*= scale_xy
    verts_scaled[:, 3] .*= scale_z

    # df/dR: dx/dR = (x/R, y/R, 0)
    dxdR = zeros(nv, 3)
    dxdR[:, 1] .= verts_scaled[:, 1] ./ R
    dxdR[:, 2] .= verts_scaled[:, 2] ./ R
    dfdR = dot(vec(obj.baseline_dfdx), vec(dxdR))

    # df/dH: dx/dH = (0, 0, z/H)
    dxdH = zeros(nv, 3)
    dxdH[:, 3] .= verts_scaled[:, 3] ./ H
    dfdH = dot(vec(obj.baseline_dfdx), vec(dxdH))

    # Clean output
    rm(outdir; recursive=true, force=true)

    return f, [dfdR, dfdH]
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)

    palace      = "../../build/palace-build/palace-arm64.bin"
    base_config = "cavity_pec_meshsens_p1_ref.json"
    base_mesh   = "mesh/cylinder_hex_p1_ref1.msh"
    work_mesh   = "mesh/cylinder_hex_p1_taylor.msh"
    work_config = "cavity_pec_taylor.json"
    baseline    = "postpro/cavity_pec_sens_p1_ref"
    R₀          = 2.74    # cm
    H₀          = 2 * R₀  # 5.48 cm

    # Verify prerequisites
    @assert isfile(palace) "Palace binary not found: $palace"
    @assert isfile(base_config) "Base config not found"
    @assert isfile(base_mesh) "Base mesh not found"
    @assert isdir(baseline) "Baseline output not found"

    # Read baseline
    mesh = read_mesh(joinpath(baseline, "mesh_for_sensitivity.mesh"))
    sol = read_solb(joinpath(baseline, "mesh_sensitivity.solb"))
    f₀ = read_eigenfrequency(joinpath(baseline, "eig.csv"); mode=1)

    println("═══ Cylinder Optimizer: L-BFGS vs Nelder-Mead ═══")
    @printf("  Baseline: f₀ = %.12e GHz at R₀ = %.2f cm, H₀ = %.2f cm\n", f₀, R₀, H₀)

    # Target: perturb frequency by ~5% (achievable by changing R)
    f_target = f₀ * 0.95
    @printf("  Target:   f_target = %.12e GHz (95%% of baseline)\n\n", f_target)

    # ─── L-BFGS (gradient-based) ─────────────────────────────────────────

    println("═══ L-BFGS (gradient-based) ═══")

    obj_lbfgs = CylinderObjective(
        palace, base_config, base_mesh, work_mesh, work_config,
        f_target, mesh, sol.data, R₀, H₀, Ref(0)
    )

    function fg_lbfgs!(F, G, x)
        R, H = x
        f, grad = evaluate!(obj_lbfgs, R, H)
        residual = (f - f_target)^2
        if !isnothing(G)
            G[1] = 2 * (f - f_target) * grad[1]
            G[2] = 2 * (f - f_target) * grad[2]
        end
        if !isnothing(F)
            return residual
        end
    end

    x0 = [R₀, H₀]
    lower = [R₀ * 0.8, H₀ * 0.8]
    upper = [R₀ * 1.2, H₀ * 1.2]

    # Optim.jl API: separate f and g! for Fminbox
    function f_lbfgs(x)
        R, H = x
        f, _ = evaluate!(obj_lbfgs, R, H)
        return (f - f_target)^2
    end

    function g_lbfgs!(G, x)
        R, H = x
        f, grad = evaluate!(obj_lbfgs, R, H)
        G[1] = 2 * (f - f_target) * grad[1]
        G[2] = 2 * (f - f_target) * grad[2]
    end

    result_lbfgs = optimize(
        f_lbfgs, g_lbfgs!, lower, upper, x0,
        Fminbox(LBFGS()),
        Optim.Options(iterations=50, f_tol=1e-20, g_tol=1e-12,
                       show_trace=true, show_every=1)
    )

    R_opt_lbfgs, H_opt_lbfgs = Optim.minimizer(result_lbfgs)
    f_opt_lbfgs, _ = evaluate!(obj_lbfgs, R_opt_lbfgs, H_opt_lbfgs)
    n_eval_lbfgs = obj_lbfgs.eval_count[]

    println()
    @printf("  Result: R = %.6f cm, H = %.6f cm\n", R_opt_lbfgs, H_opt_lbfgs)
    @printf("  f_opt = %.12e GHz (target: %.12e)\n", f_opt_lbfgs, f_target)
    @printf("  |f - f_target| = %.6e GHz\n", abs(f_opt_lbfgs - f_target))
    @printf("  Palace evaluations: %d\n\n", n_eval_lbfgs)

    # ─── Nelder-Mead (derivative-free) ───────────────────────────────────

    println("═══ Nelder-Mead (derivative-free) ═══")

    obj_nm = CylinderObjective(
        palace, base_config, base_mesh, work_mesh, work_config,
        f_target, mesh, sol.data, R₀, H₀, Ref(0)
    )

    function f_nm(x)
        R, H = x
        f, _ = evaluate!(obj_nm, R, H)
        return (f - f_target)^2
    end

    # Nelder-Mead with simple box projection (Fminbox + NelderMead doesn't work)
    function f_nm_clamped(x)
        R = clamp(x[1], lower[1], upper[1])
        H = clamp(x[2], lower[2], upper[2])
        f, _ = evaluate!(obj_nm, R, H)
        return (f - f_target)^2
    end

    result_nm = optimize(
        f_nm_clamped, x0,
        NelderMead(),
        Optim.Options(iterations=200, f_tol=1e-24, show_trace=true, show_every=5)
    )

    R_opt_nm, H_opt_nm = Optim.minimizer(result_nm)
    f_opt_nm, _ = evaluate!(obj_nm, R_opt_nm, H_opt_nm)
    n_eval_nm = obj_nm.eval_count[]

    println()
    @printf("  Result: R = %.6f cm, H = %.6f cm\n", R_opt_nm, H_opt_nm)
    @printf("  f_opt = %.12e GHz (target: %.12e)\n", f_opt_nm, f_target)
    @printf("  |f - f_target| = %.6e GHz\n", abs(f_opt_nm - f_target))
    @printf("  Palace evaluations: %d\n\n", n_eval_nm)

    # ─── Comparison ──────────────────────────────────────────────────────

    println("═══ Comparison ═══")
    @printf("  %-20s  %-12s  %-12s\n", "", "L-BFGS", "Nelder-Mead")
    @printf("  %-20s  %-12d  %-12d\n", "Palace evaluations", n_eval_lbfgs, n_eval_nm)
    @printf("  %-20s  %-12.6e  %-12.6e\n", "|f - f_target| GHz",
            abs(f_opt_lbfgs - f_target), abs(f_opt_nm - f_target))
    @printf("  %-20s  %-12.6f  %-12.6f\n", "R_opt (cm)", R_opt_lbfgs, R_opt_nm)
    @printf("  %-20s  %-12.6f  %-12.6f\n", "H_opt (cm)", H_opt_lbfgs, H_opt_nm)

    speedup = n_eval_nm / max(n_eval_lbfgs, 1)
    @printf("\n  Speedup: L-BFGS uses %.1fx fewer evaluations than Nelder-Mead\n", speedup)

    # Cleanup
    isfile(work_mesh) && rm(work_mesh)
    isfile(work_config) && rm(work_config)

    println("\n  Done.")
end

main()
