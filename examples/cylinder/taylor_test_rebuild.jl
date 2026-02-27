# Geometry rebuild Taylor test — validates the full geometry→mesh→solve pipeline.
#
# Unlike taylor_test.jl (affine mesh scaling), this script:
#   1. Rebuilds the cylinder geometry at each perturbed radius R₀+ε
#   2. Re-meshes via Gmsh (new mesh topology, different node count)
#   3. Runs Palace on the new mesh
#   4. Checks that |f(R₀+ε) - f(R₀) - ε·df/dR| = O(ε²)
#
# This catches bugs that affine scaling misses: mesh topology changes,
# tessellation quality at different parameters, Gmsh element distribution.
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev taylor_test_rebuild.jl

using LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))
include(joinpath(@__DIR__, "mesh", "mesh.jl"))

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)

    palace      = "../../build/palace-build/palace-arm64.bin"
    base_config = "cavity_pec_meshsens_p1_ref.json"
    baseline    = "postpro/cavity_pec_sens_p1_ref"
    R₀          = 2.74    # cm
    H₀          = 2 * R₀  # 5.48 cm (aspect_ratio=1.0 means height = 2*radius)

    @assert isfile(palace) "Palace binary not found: $(abspath(palace))"
    @assert isdir(baseline) "Baseline output not found: $baseline"

    # Read baseline results
    mesh = read_mesh(joinpath(baseline, "mesh_for_sensitivity.mesh"))
    sol = read_solb(joinpath(baseline, "mesh_sensitivity.solb"))
    f₀ = read_eigenfrequency(joinpath(baseline, "eig.csv"); mode=1)

    # Compute df/dR from the baseline sensitivity
    nv = size(mesh.vertices, 1)
    dfdx = sol.data
    dxdR = zeros(nv, 3)
    dxdR[:, 1] .= mesh.vertices[:, 1] ./ R₀
    dxdR[:, 2] .= mesh.vertices[:, 2] ./ R₀
    dfdR = dot(vec(dfdx), vec(dxdR))

    println("═══ Geometry Rebuild Taylor Test ═══")
    @printf("  f₀    = %.12e GHz\n", f₀)
    @printf("  R₀    = %.2f cm, H₀ = %.2f cm\n", R₀, H₀)
    @printf("  df/dR = %.12e GHz/cm  (df/dR·R/f = %.10f)\n", dfdR, dfdR * R₀ / f₀)
    println()

    # Taylor test: rebuild geometry at R₀ + ε for each ε
    epsilons = [1e-1, 1e-2, 1e-3, 1e-4]
    remainders = Float64[]
    freqs = Float64[]

    for ε in epsilons
        R_pert = R₀ + ε
        label = @sprintf("1e-%d", round(Int, -log10(ε)))
        mesh_file = "mesh/cylinder_hex_taylor_rebuild.msh"
        config_file = "cavity_pec_taylor_rebuild.json"
        outdir = "postpro/taylor_rebuild_$label"

        # Clean previous output
        isdir(outdir) && rm(outdir; recursive=true)

        # Rebuild mesh at perturbed radius
        @printf("  ε = %.0e: meshing (R=%.4f)... ", ε, R_pert)
        flush(stdout)
        generate_cylindrical_mesh(
            filename="cylinder_hex_taylor_rebuild.msh",
            refinement=2,
            order=1,
            mesh_type=2,
            radius=R_pert,
            aspect_ratio=1.0,
            verbose=0
        )

        # Write Palace config pointing to the new mesh
        text = read(base_config, String)
        text = replace(text, r"\"Mesh\":\s*\"[^\"]*\"" => "\"Mesh\": \"$mesh_file\"")
        text = replace(text, r"\"Output\":\s*\"[^\"]*\"" => "\"Output\": \"$outdir\"")
        text = replace(text, r"\"Sensitivity\":\s*true" => "\"Sensitivity\": false")
        text = replace(text, r"\"MeshSensitivity\":\s*true" => "\"MeshSensitivity\": false")
        write(config_file, text)

        # Run Palace
        print("solving... ")
        flush(stdout)
        t0 = time()
        try
            run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
        catch
            error("Palace failed at R=$R_pert (ε=$ε)")
        end
        dt = time() - t0

        eig = joinpath(outdir, "eig.csv")
        @assert isfile(eig) "eig.csv not found in $outdir"
        f_ε = read_eigenfrequency(eig; mode=1)

        r = abs(f_ε - f₀ - ε * dfdR)
        push!(remainders, r)
        push!(freqs, f_ε)
        @printf("f = %.12e, remainder = %.4e (%.1fs)\n", f_ε, r, dt)

        # Clean up
        isdir(outdir) && rm(outdir; recursive=true)
        isfile(config_file) && rm(config_file)
    end

    # Clean up mesh
    isfile("mesh/cylinder_hex_taylor_rebuild.msh") &&
        rm("mesh/cylinder_hex_taylor_rebuild.msh")

    # Convergence table
    println()
    println("─── R Convergence Table (geometry rebuild) ───")
    @printf("  %-10s  %-22s  %-14s  %s\n", "ε", "f [GHz]", "remainder", "rate")
    println("  " * "─"^58)

    rates = Float64[]
    for (k, ε) in enumerate(epsilons)
        if k == 1
            @printf("  %-10.0e  %+.12e  %.4e      —\n", ε, freqs[k], remainders[k])
        else
            rate = log10(remainders[k] / remainders[k-1]) /
                   log10(epsilons[k] / epsilons[k-1])
            push!(rates, rate)
            @printf("  %-10.0e  %+.12e  %.4e  %5.2f\n",
                    ε, freqs[k], remainders[k], rate)
        end
    end

    # Check rates
    n_good = count(r -> r > 1.8, rates)
    pass = n_good >= 2

    println()
    @printf("  Rates > 1.8: %d / %d  (need ≥ 2)\n", n_good, length(rates))
    println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")

    if !pass
        println("\n  WARNING: Geometry rebuild Taylor test did not achieve O(ε²) convergence.")
        println("  This may indicate mesh-topology-dependent errors in the sensitivity chain.")
    end

    println("\n  Done.")
    return pass
end

main()
