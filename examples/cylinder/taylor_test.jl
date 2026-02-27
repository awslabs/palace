# Taylor test for cylinder mesh sensitivity — 2 parameters (R, H).
#
# Perturbs the cylinder radius (R) and height (H) independently, re-solves
# Palace at each perturbation, and verifies the first-order Taylor remainder
# converges quadratically for both parameters:
#   |f(p+ε) - f(p) - ε·df/dp| = O(ε²)
#
# For TM₀₁₀ mode: df/dR·R/f = -1 (strong), df/dH ≈ 0 (mode is H-independent).
#
# Reads baseline f₀ and df/dx from existing Phase A2 output (no re-run needed).
# Usage: cd examples/cylinder && julia taylor_test.jl

using LinearAlgebra, Printf

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))

# ─── Binary Gmsh v2 mesh scaling ──────────────────────────────────────────

"""
    find_nodes_offset(buf) -> (pos, n_nodes)

Find the start of binary node data in a Gmsh v2 file.
Returns the byte offset of the first node and the node count.
"""
function find_nodes_offset(buf::Vector{UInt8})
    target = UInt8['$', 'N', 'o', 'd', 'e', 's', '\n']
    pos = 0
    for i in 1:length(buf) - length(target) + 1
        if @view(buf[i:i+length(target)-1]) == target
            pos = i + length(target)
            break
        end
    end
    pos == 0 && error("No \$Nodes section found")
    nl = findnext(==(UInt8('\n')), buf, pos)
    n_nodes = parse(Int, String(buf[pos:nl-1]))
    return nl + 1, n_nodes
end

read_f64(buf, pos) = reinterpret(Float64, buf[pos:pos+7])[1]

function write_f64!(buf, pos, val)
    buf[pos:pos+7] .= reinterpret(UInt8, Float64[val])
end

"""
    scale_mesh!(buf, n_nodes, pos; xy=1.0, z=1.0)

Scale node coordinates in-place. Each node: 4-byte tag + 3×8-byte doubles.
"""
function scale_mesh!(buf, n_nodes, pos; xy=1.0, z=1.0)
    for _ in 1:n_nodes
        pos += 4  # skip tag
        if xy != 1.0
            write_f64!(buf, pos,     read_f64(buf, pos)     * xy)
            write_f64!(buf, pos + 8, read_f64(buf, pos + 8) * xy)
        end
        if z != 1.0
            write_f64!(buf, pos + 16, read_f64(buf, pos + 16) * z)
        end
        pos += 24
    end
end

function write_scaled_mesh(src, dst; xy=1.0, z=1.0)
    buf = read(src)
    pos, n_nodes = find_nodes_offset(buf)
    scale_mesh!(buf, n_nodes, pos; xy=xy, z=z)
    write(dst, buf)
end

# ─── Config generation ─────────────────────────────────────────────────────

function write_taylor_config(base_config, mesh_path, output_dir, dest)
    text = read(base_config, String)
    text = replace(text, r"\"Mesh\":\s*\"[^\"]*\"" => "\"Mesh\": \"$mesh_path\"")
    text = replace(text, r"\"Output\":\s*\"[^\"]*\"" => "\"Output\": \"$output_dir\"")
    text = replace(text, r"\"Sensitivity\":\s*true" => "\"Sensitivity\": false")
    text = replace(text, r"\"MeshSensitivity\":\s*true" => "\"MeshSensitivity\": false")
    write(dest, text)
end

eps_label(ε) = @sprintf("1e-%d", round(Int, -log10(ε)))

# ─── Single-parameter Taylor test ─────────────────────────────────────────

"""
    run_taylor_test(; palace, base_config, base_mesh, work_mesh, work_config,
                      f₀, dfdp, epsilons, param_name, scale_fn, prefix)

Run a Taylor test for one parameter. Returns (pass, rates).
"""
function run_taylor_test(;
    palace, base_config, base_mesh, work_mesh, work_config,
    f₀, dfdp, epsilons, param_name, scale_fn, prefix
)
    # Clean old output
    for ε in epsilons
        dir = joinpath("postpro", "$(prefix)_$(eps_label(ε))")
        isdir(dir) && rm(dir; recursive=true)
    end

    remainders = Float64[]
    freqs      = Float64[]

    for ε in epsilons
        label  = eps_label(ε)
        outdir = "postpro/$(prefix)_$label"

        # Scale mesh
        scale_fn(base_mesh, work_mesh, ε)

        # Write config
        write_taylor_config(base_config, "mesh/cylinder_hex_p1_taylor.msh",
                            outdir, work_config)

        # Run Palace
        @printf("  ε = %.0e: running Palace... ", ε)
        try
            run(pipeline(`$palace $work_config`; stdout=devnull, stderr=devnull))
        catch
            error("Palace failed for $param_name perturbation ε=$ε")
        end

        eig = joinpath(outdir, "eig.csv")
        @assert isfile(eig) "eig.csv not found in $outdir"
        f_ε = read_eigenfrequency(eig; mode=1)

        r = abs(f_ε - f₀ - ε * dfdp)
        push!(remainders, r)
        push!(freqs, f_ε)
        @printf("f = %.12e, remainder = %.4e\n", f_ε, r)
    end

    # Convergence table
    println()
    println("─── $param_name Convergence Table ───")
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

    # Pass criterion depends on gradient magnitude.
    # Zero gradient: verify remainder drops below noise (gradient correctly = 0).
    # Nonzero gradient: verify O(ε²) convergence rate.
    if abs(dfdp) < 1e-10 && remainders[end] < 1e-12
        pass = true
        println()
        println("  Gradient ≈ 0: remainder at smallest ε = $(remainders[end]) (< 1e-12)")
        println("  Status: PASS ✓ (zero gradient confirmed)")
    else
        n_good = count(>(1.8), rates)
        pass   = n_good >= 2
        println()
        println("  Rates > 1.8: $n_good / $(length(rates))  (need ≥ 2)")
        println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    end
    return pass, rates
end

# ─── Main ──────────────────────────────────────────────────────────────────

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
    epsilons    = [1e-1, 1e-2, 1e-3, 1e-4]

    # Verify prerequisites
    @assert isfile(palace)      "Palace binary not found: $palace"
    @assert isfile(base_config) "Base config not found: $base_config"
    @assert isfile(base_mesh)   "Base mesh not found: $base_mesh"
    @assert isdir(baseline)     "Baseline output not found — run Phase A2 first"
    for f in ("eig.csv", "mesh_for_sensitivity.mesh", "mesh_sensitivity.solb")
        @assert isfile(joinpath(baseline, f)) "$f not found in $baseline"
    end

    # ─── Read baseline ───────────────────────────────────────────────────

    println("═══ Taylor Test: 2-Parameter Cylinder Mesh Sensitivity ═══")
    println()

    f₀   = read_eigenfrequency(joinpath(baseline, "eig.csv"); mode=1)
    mesh  = read_mesh(joinpath(baseline, "mesh_for_sensitivity.mesh"))
    sol   = read_solb(joinpath(baseline, "mesh_sensitivity.solb"))
    nv    = size(mesh.vertices, 1)

    # Chain contraction for R: dx/dR = (x/R, y/R, 0)
    dxdR = zeros(nv, 3)
    for i in 1:nv
        dxdR[i, 1] = mesh.vertices[i, 1] / R₀
        dxdR[i, 2] = mesh.vertices[i, 2] / R₀
    end
    dfdR = dot(vec(sol.data), vec(dxdR))

    # Chain contraction for H: dx/dH = (0, 0, z/H)
    dxdH = zeros(nv, 3)
    for i in 1:nv
        dxdH[i, 3] = mesh.vertices[i, 3] / H₀
    end
    dfdH = dot(vec(sol.data), vec(dxdH))

    @printf("  f₀      = %.12e GHz\n", f₀)
    @printf("  R₀      = %.2f cm,  H₀ = %.2f cm\n", R₀, H₀)
    @printf("  df/dR   = %.12e GHz/cm  (df/dR·R/f = %.10f, expect -1)\n", dfdR, dfdR*R₀/f₀)
    @printf("  df/dH   = %.12e GHz/cm  (expect ≈ 0 for TM₀₁₀)\n\n", dfdH)

    # ─── Taylor test for R ───────────────────────────────────────────────

    println("═══ Parameter 1: Radius (R) ═══")
    pass_R, _ = run_taylor_test(;
        palace, base_config, base_mesh, work_mesh, work_config,
        f₀, dfdp=dfdR, epsilons, param_name="R",
        scale_fn=(src, dst, ε) -> write_scaled_mesh(src, dst; xy=(R₀ + ε) / R₀),
        prefix="taylor_R",
    )

    # ─── Taylor test for H ───────────────────────────────────────────────

    println()
    println("═══ Parameter 2: Height (H) ═══")
    pass_H, _ = run_taylor_test(;
        palace, base_config, base_mesh, work_mesh, work_config,
        f₀, dfdp=dfdH, epsilons, param_name="H",
        scale_fn=(src, dst, ε) -> write_scaled_mesh(src, dst; z=(H₀ + ε) / H₀),
        prefix="taylor_H",
    )

    # ─── Summary ─────────────────────────────────────────────────────────

    println()
    println("═══ Summary ═══")
    println("  R (radius): $(pass_R ? "PASS ✓" : "FAIL ✗")")
    println("  H (height): $(pass_H ? "PASS ✓" : "FAIL ✗")")

    # Cleanup
    isfile(work_mesh) && rm(work_mesh)
    isfile(work_config) && rm(work_config)

    all_pass = pass_R && pass_H
    all_pass || error("Taylor test FAILED")
    return println("\n  All tests passed.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
