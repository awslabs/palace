# Gradient-based optimization of transmon cavity eigenfrequency.
#
# Uses mesh-perturbation (same mesh, shifted vertices) instead of geometry-rebuild
# to eliminate TetGen re-meshing noise. This is the key advantage of the derivative
# pipeline: the gradient lives on a FIXED mesh, enabling clean convergence.
#
# Pipeline:
#   1. Build baseline EGADS geometry → TetGen mesh → Palace (with MeshSensitivity)
#   2. Compute dx/dp via EGADS provenance-based chain
#   3. At each optimizer step: perturb mesh by Δp * dx/dp, re-run Palace
#   4. Read new df/dx, recompute gradient for the next step
#
# Compares L-BFGS (gradient-based) vs Nelder-Mead (derivative-free).
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev transmon_optimize.jl

using LinearAlgebra, Printf, Optim
using EGADS
using DeviceLayout
using TetGen

const DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

include(joinpath(@__DIR__, "shared", "palace_io.jl"))
include(joinpath(@__DIR__, "transmon_sensitivity_test.jl"))

# ─── Write a perturbed mesh ─────────────────────────────────────────────────

"""Write a .mesh file with vertices shifted by `delta * dxdp`."""
function write_perturbed_mesh(filename, base_mesh, dxdp, delta)
    nv = size(base_mesh.vertices, 1)
    open(filename, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        println(io, "Vertices\n$nv")
        for vi in 1:nv
            x = base_mesh.vertices[vi, 1] + delta * dxdp[vi, 1]
            y = base_mesh.vertices[vi, 2] + delta * dxdp[vi, 2]
            z = base_mesh.vertices[vi, 3] + delta * dxdp[vi, 3]
            @printf(io, "%.15e %.15e %.15e 0\n", x, y, z)
        end
        _write_elements(io, "Tetrahedra", base_mesh.tetrahedra, base_mesh.tet_refs, 4)
        _write_elements(io, "Hexahedra", base_mesh.hexahedra, base_mesh.hex_refs, 8)
        _write_elements(io, "Triangles", base_mesh.triangles, base_mesh.tri_refs, 3)
        _write_elements(io, "Quadrilaterals", base_mesh.quadrilaterals, base_mesh.quad_refs, 4)
        println(io, "\nEnd")
    end
end

function _write_elements(io, keyword, elems, refs, ncols)
    ne = size(elems, 1)
    ne > 0 || return
    println(io, "\n$keyword\n$ne")
    for e in 1:ne
        for c in 1:ncols
            print(io, elems[e, c], " ")
        end
        println(io, refs[e])
    end
end

# ─── Evaluate at a perturbed mesh ───────────────────────────────────────────

"""
Run Palace on a mesh perturbed by delta from baseline.
Returns frequency (and optionally gradient via chain contraction at the new mesh).
"""
function eval_perturbed(palace_bin, base_mesh, dxdp_full, nfaces, delta;
                         label="eval", with_sensitivity=false)
    mesh_file = "mesh/transmon_$(label).mesh"
    config_file = "transmon_$(label).json"
    outdir = "postpro/transmon_$(label)"

    mkpath("mesh")
    write_perturbed_mesh(mesh_file, base_mesh, dxdp_full, delta)
    write_config(config_file, mesh_file, nfaces; output_dir=outdir)

    if !with_sensitivity
        config_str = read(config_file, String)
        config_str = replace(config_str, "\"MeshSensitivity\": true" => "\"MeshSensitivity\": false")
        config_str = replace(config_str, "\"Sensitivity\": true" => "\"Sensitivity\": false")
        write(config_file, config_str)
    end

    mkpath(outdir)
    run(pipeline(`$palace_bin $config_file`; stdout=devnull, stderr=devnull))
    freq = read_eigenfrequency(joinpath(outdir, "eig.csv"); mode=1)

    # If sensitivity was requested, read the new df/dx and compute updated gradient
    grad = nothing
    if with_sensitivity
        solb = joinpath(outdir, "mesh_sensitivity.solb")
        mesh_p = joinpath(outdir, "mesh_for_sensitivity.mesh")
        if isfile(solb) && isfile(mesh_p)
            sol = read_solb(solb)
            # Gradient at the perturbed point: df/dp = dot(dfdx_new, dxdp)
            # dxdp is the same shape direction — valid for small perturbations
            grad = dot(vec(sol.data), vec(dxdp_full[1:size(sol.data, 1), :]))
        end
    end

    rm(outdir; recursive=true, force=true)
    rm(config_file; force=true)
    rm(mesh_file; force=true)

    return freq, grad
end

# ─── Main ───────────────────────────────────────────────────────────────────

function optimize_main()
    cd(@__DIR__)

    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found"

    cap_length_0 = 1.0  # baseline cap_length (cm)

    println("═══ Transmon Cavity Optimizer: L-BFGS vs Nelder-Mead ═══\n")

    # ─── Step 1: Baseline ────────────────────────────────────────────────

    println("Step 1: Baseline geometry + sensitivity...")
    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)
    baseline = run_palace_eval(ctx, cap_length_0; label="opt_baseline", palace_bin=palace)

    # Read Palace mesh sensitivity
    solb_file = joinpath(baseline.outdir, "mesh_sensitivity.solb")
    mesh_palace_file = joinpath(baseline.outdir, "mesh_for_sensitivity.mesh")
    palace_mesh = read_mesh(mesh_palace_file)
    sol = read_solb(solb_file)
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv = size(palace_mesh.vertices, 1)

    # Compute dx/dp at baseline
    geo = build_geometry(ctx; cap_length=cap_length_0)
    chain = compute_chain(geo, geo.result, palace_mesh, dfdx, is_boundary;
                           cap_length=cap_length_0)
    EGADS.close(ctx)

    f0 = baseline.freq
    dfdp0 = chain.dfdp
    nfaces = baseline.minfo.nfaces

    @printf("  f₀ = %.6f GHz, df/dp = %.6e GHz/cm\n", f0, dfdp0)
    @printf("  Palace mesh: %d verts, %d boundary\n\n", nv, count(is_boundary))

    # Target: 95% of baseline
    f_target = f0 * 0.95
    @printf("  Target: f_target = %.6f GHz (95%%)\n\n", f_target)

    # Clean up baseline output (keep mesh data in memory)
    rm(baseline.outdir; recursive=true, force=true)
    rm(baseline.config_file; force=true)
    rm(baseline.mesh_file; force=true)

    # ─── Step 2: L-BFGS (gradient-based, mesh-perturbation) ─────────────

    println("═══ L-BFGS (gradient-based) ═══")
    n_eval_lbfgs = Ref(0)

    function f_lbfgs(x)
        delta = x[1]  # perturbation magnitude (cap_length = cap_length_0 + delta)
        n_eval_lbfgs[] += 1
        freq, _ = eval_perturbed(palace, palace_mesh, chain.dxdp_full, nfaces, delta;
                                  label=@sprintf("lbfgs_%04d", n_eval_lbfgs[]))
        residual = (freq - f_target)^2
        @printf("  [%d] Δp=%.4f  f=%.6f  |f-ft|=%.4e\n",
                n_eval_lbfgs[], delta, freq, abs(freq - f_target))
        return residual
    end

    function g_lbfgs!(G, x)
        delta = x[1]
        n_eval_lbfgs[] += 1
        freq, grad = eval_perturbed(palace, palace_mesh, chain.dxdp_full, nfaces, delta;
                                     label=@sprintf("lbfgs_%04d", n_eval_lbfgs[]),
                                     with_sensitivity=true)
        if isnothing(grad)
            G[1] = 2 * (freq - f_target) * dfdp0  # fallback to baseline gradient
        else
            G[1] = 2 * (freq - f_target) * grad
        end
        @printf("  [%d] Δp=%.4f  f=%.6f  grad=%.4e\n",
                n_eval_lbfgs[], delta, freq, isnothing(grad) ? dfdp0 : grad)
    end

    # Optimize in delta-space: cap_length = cap_length_0 + delta
    # Bounds: ±30% of cap_length_0 → delta ∈ [-0.3, 0.3]
    result_lbfgs = optimize(
        f_lbfgs, g_lbfgs!,
        [-0.3],   # lower
        [0.3],    # upper
        [0.0],    # start at baseline
        Fminbox(LBFGS()),
        Optim.Options(iterations=15, f_abstol=1e-16, g_tol=1e-10)
    )

    delta_opt = Optim.minimizer(result_lbfgs)[1]
    cl_opt_lbfgs = cap_length_0 + delta_opt
    n_eval_lbfgs_final = n_eval_lbfgs[]

    # Final evaluation at optimum
    n_eval_lbfgs[] += 1
    f_opt_lbfgs, _ = eval_perturbed(palace, palace_mesh, chain.dxdp_full, nfaces, delta_opt;
                                     label="lbfgs_final")

    @printf("\n  Result: cap_length = %.4f cm (Δp = %.4f)\n", cl_opt_lbfgs, delta_opt)
    @printf("  f_opt = %.6f GHz (target: %.6f)\n", f_opt_lbfgs, f_target)
    @printf("  |f - f_target| = %.6e GHz\n", abs(f_opt_lbfgs - f_target))
    @printf("  Palace evaluations: %d\n\n", n_eval_lbfgs[])

    # ─── Step 3: Nelder-Mead (derivative-free, mesh-perturbation) ────────

    println("═══ Nelder-Mead (derivative-free) ═══")
    n_eval_nm = Ref(0)

    function f_nm(x)
        delta = clamp(x[1], -0.3, 0.3)
        n_eval_nm[] += 1
        freq, _ = eval_perturbed(palace, palace_mesh, chain.dxdp_full, nfaces, delta;
                                  label=@sprintf("nm_%04d", n_eval_nm[]))
        residual = (freq - f_target)^2
        @printf("  [%d] Δp=%.4f  f=%.6f  |f-ft|=%.4e\n",
                n_eval_nm[], delta, freq, abs(freq - f_target))
        return residual
    end

    result_nm = optimize(
        f_nm, [0.0],
        NelderMead(),
        Optim.Options(iterations=30, f_abstol=1e-16)
    )

    delta_opt_nm = Optim.minimizer(result_nm)[1]
    cl_opt_nm = cap_length_0 + delta_opt_nm

    n_eval_nm[] += 1
    f_opt_nm, _ = eval_perturbed(palace, palace_mesh, chain.dxdp_full, nfaces, delta_opt_nm;
                                  label="nm_final")

    @printf("\n  Result: cap_length = %.4f cm (Δp = %.4f)\n", cl_opt_nm, delta_opt_nm)
    @printf("  f_opt = %.6f GHz (target: %.6f)\n", f_opt_nm, f_target)
    @printf("  |f - f_target| = %.6e GHz\n", abs(f_opt_nm - f_target))
    @printf("  Palace evaluations: %d\n\n", n_eval_nm[])

    # ─── Comparison ──────────────────────────────────────────────────────

    println("═══ Comparison ═══")
    @printf("  %-22s  %-14s  %-14s\n", "", "L-BFGS", "Nelder-Mead")
    @printf("  %-22s  %-14d  %-14d\n", "Palace evaluations", n_eval_lbfgs[], n_eval_nm[])
    @printf("  %-22s  %-14.6e  %-14.6e\n", "|f - f_target| GHz",
            abs(f_opt_lbfgs - f_target), abs(f_opt_nm - f_target))
    @printf("  %-22s  %-14.4f  %-14.4f\n", "cap_length_opt (cm)",
            cl_opt_lbfgs, cl_opt_nm)

    speedup = n_eval_nm[] / max(n_eval_lbfgs[], 1)
    @printf("\n  Speedup: L-BFGS uses %.1f× fewer evaluations\n", speedup)

    println("\n  Done.")
end

optimize_main()
