# End-to-end test: EGADS geometry → TetGen volume mesh → Palace eigensolve → sensitivity chain.
#
# Step 1: Build a box via EGADS (EGADSKernel would do this, but we start with raw EGADS)
# Step 2: Tessellate → surface triangles
# Step 3: TetGen → volume tets
# Step 4: Write .mesh file for Palace
# Step 5: Run Palace with MeshSensitivity
# Step 6: dx/dp from EGADS evaluate_dot → chain contraction
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev egads_pipeline_test.jl

using LinearAlgebra, Printf
using EGADS
using TetGen

# Shared Palace mesh/solb readers and boundary utilities
include(joinpath(@__DIR__, "shared", "palace_io.jl"))

# ─── Step 1: Build box geometry via EGADS ─────────────────────────────────

function build_box(ctx, Lx, Ly, Lz)
    body = EGADS.make_box(ctx, 0.0, 0.0, 0.0, Lx, Ly, Lz)
    data = Float64[0.0, 0.0, 0.0, Lx, Ly, Lz]
    return body, data
end

# ─── Step 2+3: Tessellate → TetGen → volume mesh ─────────────────────────

function make_volume_mesh(body; max_edge=0.5, sag=0.01, angle=15.0)
    # Tessellate the EGADS body to get surface triangles
    tess = EGADS.tessellate(body; max_edge=max_edge, sag=sag, angle=angle)
    status = EGADS.tess_status(tess)
    npts = status.npts
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)

    # Collect all surface vertices and triangles
    # Build global vertex list from tessellation
    pts = zeros(Float64, 3, npts)
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end

    # Collect surface triangles with face IDs as markers
    tris = Int32[]
    tri_markers = Int32[]
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        ntri = size(fm.tris, 1)
        for t in 1:ntri
            # Convert local triangle indices to global
            for k in 1:3
                ig = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])
                push!(tris, Int32(ig))
            end
            push!(tri_markers, Int32(fi))  # face index as boundary marker
        end
    end
    ntri_total = length(tri_markers)

    println("  Surface: $npts vertices, $ntri_total triangles, $(length(faces)) faces")

    # Build TetGen input
    input = RawTetGenIO{Float64}()
    input.pointlist = pts

    # Build facet matrix: 3×ntri (each column is a triangle)
    tri_matrix = reshape(tris, 3, ntri_total)
    facetlist!(input, tri_matrix)

    # Set triangle markers
    input.facetmarkerlist = tri_markers

    # Tetrahedralize: 'p' = PLC, 'q' = quality, 'a' = max volume
    output = tetrahedralize(input, "pq1.2a0.01")

    nvol = numberofpoints(output)
    ntet = numberoftetrahedra(output)
    nsurf = numberoftrifaces(output)
    println("  Volume:  $nvol vertices, $ntet tets, $nsurf surface tris")

    return output, tess, npts
end

# ─── Step 4: Write .mesh file ─────────────────────────────────────────────

function write_palace_mesh(filename, tetgen_out, n_surface_pts)
    pts = tetgen_out.pointlist
    nvol = size(pts, 2)

    # Tet connectivity (1-based)
    tets = tetgen_out.tetrahedronlist
    ntet = size(tets, 2)

    # Surface triangles with markers
    trifaces = tetgen_out.trifacelist
    trimarkers = tetgen_out.trifacemarkerlist
    ntri = size(trifaces, 2)

    open(filename, "w") do io
        println(io, "MeshVersionFormatted 2")
        println(io)
        println(io, "Dimension 3")
        println(io)

        # Vertices
        println(io, "Vertices")
        println(io, nvol)
        for i in 1:nvol
            @printf(io, "%.15e %.15e %.15e 0\n", pts[1, i], pts[2, i], pts[3, i])
        end
        println(io)

        # Tetrahedra (attribute 1 = volume material)
        println(io, "Tetrahedra")
        println(io, ntet)
        for i in 1:ntet
            @printf(io, "%d %d %d %d 1\n", tets[1, i], tets[2, i], tets[3, i], tets[4, i])
        end
        println(io)

        # Surface triangles with face markers as boundary attributes
        println(io, "Triangles")
        println(io, ntri)
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", trifaces[1, i], trifaces[2, i], trifaces[3, i], trimarkers[i])
        end
        println(io)

        println(io, "End")
    end
    println("  Wrote $filename ($nvol verts, $ntet tets, $ntri surface tris)")
end

# ─── Step 5+6: Palace config + run ───────────────────────────────────────

function write_palace_config(filename, mesh_path, Lx, Ly, Lz, nfaces)
    # PEC on all faces — simple resonant cavity
    face_attrs = collect(1:nfaces)

    config = """{
  "Problem": {
    "Type": "Eigenmode",
    "Verbose": 2,
    "Output": "postpro/egads_box"
  },
  "Model": {
    "Mesh": "$mesh_path",
    "L0": 1.0e-2
  },
  "Domains": {
    "Materials": [
      {
        "Attributes": [1],
        "Permeability": 1.0,
        "Permittivity": 1.0
      }
    ]
  },
  "Boundaries": {
    "PEC": {
      "Attributes": $(face_attrs)
    }
  },
  "Solver": {
    "Order": 1,
    "Device": "CPU",
    "Eigenmode": {
      "N": 2,
      "Tol": 1.0e-8,
      "Target": 2.0,
      "Sensitivity": true,
      "MeshSensitivity": true
    },
    "Linear": {
      "Type": "Default",
      "Tol": 1.0e-8,
      "MaxIts": 200
    }
  }
}"""
    write(filename, config)
    println("  Wrote $filename")
end

# ─── Main ─────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)

    Lx = 2.0  # cm
    Ly = 3.0  # cm
    Lz = 4.0  # cm

    println("═══ EGADS → TetGen → Palace Pipeline Test ═══")
    println("  Box: $Lx × $Ly × $Lz cm")
    println()

    # Step 1: Build EGADS body
    println("Step 1: Building EGADS box...")
    ctx = EGADS.EGADSContext()
    body, box_data = build_box(ctx, Lx, Ly, Lz)
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    println("  Body created: $(length(faces)) faces")

    # Step 2+3: Tessellate + TetGen
    println("\nStep 2+3: Tessellation + TetGen volume meshing...")
    tetgen_out, tess, n_surface_pts = make_volume_mesh(body; max_edge=0.5)

    # Step 4: Write mesh
    println("\nStep 4: Writing Palace mesh...")
    mesh_file = "mesh/egads_box.mesh"
    mkpath("mesh")
    write_palace_mesh(mesh_file, tetgen_out, n_surface_pts)

    # Step 5: Write Palace config
    println("\nStep 5: Writing Palace config...")
    config_file = "egads_box.json"
    write_palace_config(config_file, mesh_file, Lx, Ly, Lz, length(faces))

    # Step 6: Run Palace
    palace = "../build/palace-build/palace-arm64.bin"
    if !isfile(palace)
        println("\n  Palace binary not found at $palace — skipping solve")
        EGADS.close(ctx)
        return
    end

    println("\nStep 6: Running Palace eigensolve...")
    flush(stdout)
    mkpath("postpro/egads_box")
    t0 = time()
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
    dt = time() - t0
    println("  Palace completed in $(round(dt; digits=1))s")

    # Read results
    eig_file = "postpro/egads_box/eig.csv"
    if isfile(eig_file)
        lines = readlines(eig_file)
        for line in lines
            stripped = strip(line)
            startswith(stripped, "m,") && continue
            isempty(stripped) && continue
            tokens = split(stripped, ",")
            m = parse(Int, strip(tokens[1]) |> x -> replace(x, ".00e+00" => ""))
            f = parse(Float64, strip(tokens[2]))
            @printf("  Mode %d: f = %.6f GHz\n", m, f)
        end
    end

    # Check for mesh sensitivity output
    solb_file = "postpro/egads_box/mesh_sensitivity.solb"
    mesh_sens = "postpro/egads_box/mesh_for_sensitivity.mesh"
    if isfile(solb_file) && isfile(mesh_sens)
        println("\n  Mesh sensitivity files found!")
        println("  → $solb_file")
        println("  → $mesh_sens")
    else
        println("\n  WARNING: No mesh sensitivity output files")
        println("  Checking what was produced...")
        for f in readdir("postpro/egads_box")
            println("    $f")
        end
    end

    # ─── Step 7: Chain contraction ────────────────────────────────────────

    solb_file = "postpro/egads_box/mesh_sensitivity.solb"
    mesh_sens_file = "postpro/egads_box/mesh_for_sensitivity.mesh"

    if isfile(solb_file) && isfile(mesh_sens_file)
        println("\nStep 7: Chain contraction (df/dLy via EGADS sensitivity)...")

        # Read Palace mesh sensitivity df/dx
        palace_mesh = read_mesh(mesh_sens_file)
        sol = read_solb(solb_file)
        dfdx = sol.data  # N×3 matrix
        nv_palace = size(palace_mesh.vertices, 1)
        println("  Palace mesh: $nv_palace vertices, sol: $(sol.npts) pts")

        # Set EGADS body derivatives for d/dLy
        # BOX data: [x0, y0, z0, dx, dy, dz], d/dLy: [0, 0, 0, 0, 1, 0]
        # Mode 1 (TE₀₁₁) depends on Ly and Lz, not Lx.
        data_dot_Ly = Float64[0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        EGADS.make_solid_body_dot!(body, EGADS.SolidType.BOX, box_data, data_dot_Ly)

        faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
        println("  Faces: $(length(faces))")

        # For each Palace boundary vertex, project to EGADS face and get evaluate_dot
        is_boundary = classify_boundary_vertices(palace_mesh)
        n_boundary = count(is_boundary)
        println("  Boundary vertices: $n_boundary / $nv_palace")

        dxdp = zeros(Float64, nv_palace, 3)
        n_projected = 0

        for vi in findall(is_boundary)
            pt = palace_mesh.vertices[vi, :]  # Both EGADS body and mesh are in cm
            for fi in eachindex(faces)
                uv, result = EGADS.inv_evaluate(faces[fi], pt)
                dist = norm(result[1:3] - pt)
                if dist < 1e-6  # projection tolerance
                    res = EGADS.evaluate_dot(faces[fi], uv, zeros(length(uv)))
                    dxdp[vi, :] .= res.result_dot.xyz
                    n_projected += 1
                    break
                end
            end
        end

        println("  Projected $n_projected / $n_boundary boundary vertices")

        # Contract: df/dLy = Σ (df/dx_i) · (dx_i/dLy)
        dfdLy = dot(vec(dfdx), vec(dxdp))
        freq = parse(Float64, strip(split(readlines("postpro/egads_box/eig.csv")[2], ",")[2]))

        # Analytical: TE₀₁₁ has f = c/2 * sqrt((1/Ly)² + (1/Lz)²)
        # df/dLy * Ly / f = -(1/Ly²) / ((1/Ly²) + (1/Lz²))
        expected_nondim = -(1/Ly^2) / (1/Ly^2 + 1/Lz^2)
        nondim = dfdLy * Ly / freq
        @printf("\n  df/dLy = %.6e GHz/cm\n", dfdLy)
        @printf("  df/dLy * Ly / f = %.6f  (expect %.4f for TE₀₁₁)\n", nondim, expected_nondim)
        @printf("  Error: %.6f\n", abs(nondim - expected_nondim))

        pass = abs(nondim - expected_nondim) < 0.05
        println("  Status: $(pass ? "PASS ✓" : "FAIL ✗")")
    end

    EGADS.close(ctx)
    println("\n  Done.")
end

main()
