# Sensitivity chain verification: analytical dx/dR for cylinder TM₀₁₀ mode.
#
# Reads Palace mesh sensitivity output (df/dx at each vertex from .solb) and
# the companion mesh (vertex coordinates from .mesh), computes the analytical
# mesh deformation dx/dR = (x/R, y/R, 0) for a cylinder radius change, contracts
# to get df/dR, and verifies df/dR * R / f ≈ -1.0.
#
# This is a self-contained script — no DeviceLayout or EGADS dependencies.
# Usage: julia examples/cylinder/chain_test.jl [postpro_dir] [radius]

using LinearAlgebra, DelimitedFiles, SparseArrays

include(joinpath(@__DIR__, "..", "shared", "palace_io.jl"))

# ─── Main test ────────────────────────────────────────────────────────────────

function main()
    # Parse arguments
    postpro_dir =
        length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "postpro", "cavity_pec_sens")
    radius = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 2.74  # cm, matches mesh.jl default

    mesh_file = joinpath(postpro_dir, "mesh_for_sensitivity.mesh")
    solb_file = joinpath(postpro_dir, "mesh_sensitivity.solb")
    eig_file = joinpath(postpro_dir, "eig.csv")

    println("═══ Sensitivity Chain Test ═══")
    println("  Mesh file: $mesh_file")
    println("  Solb file: $solb_file")
    println("  Eig file:  $eig_file")
    println("  Radius:    $radius cm")
    println()

    # ─── Step 3: Analytical chain (no EGADS, no Laplacian) ────────────────

    # Read mesh and sensitivity
    @assert isfile(mesh_file) "Mesh file not found: $mesh_file\nRun Palace with companion .mesh output first."
    @assert isfile(solb_file) "Solb file not found: $solb_file"
    @assert isfile(eig_file) "Eigenfrequency file not found: $eig_file"

    mesh = read_mesh(mesh_file)
    sol = read_solb(solb_file)
    freq = read_eigenfrequency(eig_file; mode=1)

    nv = size(mesh.vertices, 1)
    println("  Vertices:    $nv")
    println("  Hexahedra:   $(size(mesh.hexahedra, 1))")
    println("  Tetrahedra:  $(size(mesh.tetrahedra, 1))")
    println("  Quads:       $(size(mesh.quadrilaterals, 1))")
    println("  Triangles:   $(size(mesh.triangles, 1))")
    println("  Sol npts:    $(sol.npts)")
    println("  Frequency:   $freq GHz")
    println()

    @assert nv == sol.npts "Vertex count mismatch: mesh=$nv, solb=$(sol.npts)"

    dfdx = sol.data  # Nx3 matrix of df/dx at each vertex

    # Analytical dx/dR for cylinder radius change: radial scaling (x/R, y/R, 0)
    # This is exact for ALL vertices (surface + interior) because the natural
    # deformation for a uniform radius change is uniform radial scaling.
    R = radius  # cm (mesh coordinates are in cm)
    dxdR_analytical = zeros(Float64, nv, 3)
    for i = 1:nv
        dxdR_analytical[i, 1] = mesh.vertices[i, 1] / R
        dxdR_analytical[i, 2] = mesh.vertices[i, 2] / R
        # dz/dR = 0 for pure radial scaling (height doesn't depend on radius)
    end

    # Contract: df/dR = Σᵢ dot(dfdx[i,:], dxdR[i,:])
    dfdR = dot(vec(dfdx), vec(dxdR_analytical))

    # Nondimensional check: df/dR * R / f should be -1 for TM₀₁₀
    # The eigenfrequency of TM₀₁₀ is f = c * x₀₁ / (2π R √εᵣ) where x₀₁ ≈ 2.405
    # So df/dR = -f/R, giving df/dR * R/f = -1.
    nondim = dfdR * R / freq

    println("─── Step 3: Analytical chain (all vertices) ───")
    println("  df/dR           = $dfdR GHz/cm")
    println("  R               = $R cm")
    println("  f               = $freq GHz")
    println("  df/dR * R / f   = $nondim")
    println("  Expected:         -1.0")
    println("  Error:            $(abs(nondim + 1.0))")

    step3_pass = abs(nondim + 1.0) < 1e-3
    println("  Status:           $(step3_pass ? "PASS ✓" : "FAIL ✗")")
    println()

    # ─── Step 4: Laplacian propagation verification ───────────────────────

    is_boundary = classify_boundary_vertices(mesh)
    n_boundary = count(is_boundary)
    n_interior = nv - n_boundary
    println("─── Step 4: Laplacian propagation ───")
    println("  Boundary vertices: $n_boundary")
    println("  Interior vertices: $n_interior")

    if n_interior > 0
        # Use analytical boundary values, propagate to interior via Laplacian
        boundary_dxdR = copy(dxdR_analytical)
        # Zero out interior values to simulate boundary-only input
        for i = 1:nv
            if !is_boundary[i]
                boundary_dxdR[i, :] .= 0.0
            end
        end

        propagated = propagate_boundary_sensitivity(mesh, boundary_dxdR, is_boundary)

        # Compare propagated interior values to analytical
        errors = Float64[]
        for i = 1:nv
            if !is_boundary[i]
                err = norm(propagated[i, :] - dxdR_analytical[i, :])
                push!(errors, err)
            end
        end

        l2_err = sqrt(sum(e^2 for e in errors) / length(errors))
        max_err = maximum(errors)

        println("  L2 error (interior): $l2_err")
        println("  Max error (interior): $max_err")

        println("  (Per-vertex errors expected for non-convex hex mesh)")

        # The real test: verify the chain contraction with propagated dx/dR
        dxdR_propagated = copy(propagated)
        dfdR_propagated = dot(vec(dfdx), vec(dxdR_propagated))
        nondim_propagated = dfdR_propagated * R / freq
        println("  df/dR * R / f (propagated): $nondim_propagated")
        println("  Error vs -1.0:               $(abs(nondim_propagated + 1.0))")

        # Pass criterion: chain contraction with Laplacian-propagated dx/dR ≈ -1.0
        step4_pass = abs(nondim_propagated + 1.0) < 1e-1
        println("  Status:              $(step4_pass ? "PASS ✓" : "FAIL ✗")")
    else
        println("  No interior vertices — Laplacian test skipped (all boundary)")
        step4_pass = true
    end
    println()

    # ─── Summary ──────────────────────────────────────────────────────────

    println("═══ Summary ═══")
    println("  Step 3 (analytical chain):     $(step3_pass ? "PASS" : "FAIL")")
    println("  Step 4 (Laplacian propagation): $(step4_pass ? "PASS" : "FAIL")")

    all_pass = step3_pass && step4_pass
    if !all_pass
        error("Chain test FAILED")
    end
    return println("\n  All tests passed.")
end

# Run main() when executed directly (not when included by another script)
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
