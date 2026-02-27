# Extruded waveguide sensitivity: EGADS extrude → TetGen → Palace → chain contraction.
#
# Builds a rectangular waveguide via EGADS extrude (not a BOX primitive).
# Verifies df/da, df/db, df/dL against analytical waveguide eigenfrequencies.
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev waveguide_extrude_test.jl

using LinearAlgebra, Printf
using EGADS
using TetGen

include(joinpath(@__DIR__, "shared", "palace_io.jl"))

# ─── Step 1: Build extruded waveguide via EGADS ──────────────────────────

function build_waveguide(ctx, a, b, L)
    # Rectangle in XY plane: (0,0,0) to (a,b,0)
    n1 = EGADS.make_node(ctx, 0.0, 0.0, 0.0)
    n2 = EGADS.make_node(ctx, a,   0.0, 0.0)
    n3 = EGADS.make_node(ctx, a,   b,   0.0)
    n4 = EGADS.make_node(ctx, 0.0, b,   0.0)

    e1 = EGADS.make_edge_from_nodes(ctx, n1, n2)
    e2 = EGADS.make_edge_from_nodes(ctx, n2, n3)
    e3 = EGADS.make_edge_from_nodes(ctx, n3, n4)
    e4 = EGADS.make_edge_from_nodes(ctx, n4, n1)

    loop = EGADS.make_wire_loop(ctx, [e1, e2, e3, e4])
    face = EGADS.make_plane_face(ctx, loop)
    src_body = EGADS.make_face_body(ctx, face)

    # Extrude along Z
    dir = [0.0, 0.0, 1.0]
    body = EGADS.extrude(src_body, L, dir)

    return body, src_body, dir
end

# ─── Step 2: Tessellate + TetGen + write .mesh ──────────────────────────

function mesh_and_write(body, mesh_file; max_edge=0.3, sag=0.01, angle=15.0,
                        tetgen_opts="pq1.2a0.01")
    tess = EGADS.tessellate(body; max_edge=max_edge, sag=sag, angle=angle)
    status = EGADS.tess_status(tess)
    npts = status.npts
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)

    # Collect surface vertices
    pts = zeros(Float64, 3, npts)
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end

    # Collect surface triangles with face markers
    tris = Int32[]
    tri_markers = Int32[]
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for t in 1:size(fm.tris, 1)
            for k in 1:3
                ig = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])
                push!(tris, Int32(ig))
            end
            push!(tri_markers, Int32(fi))
        end
    end

    # TetGen
    input = RawTetGenIO{Float64}()
    input.pointlist = pts
    facetlist!(input, reshape(tris, 3, length(tri_markers)))
    input.facetmarkerlist = tri_markers
    output = tetrahedralize(input, tetgen_opts)

    nvol = numberofpoints(output)
    ntet = numberoftetrahedra(output)
    ntri = numberoftrifaces(output)

    # Write .mesh
    mkpath(dirname(mesh_file))
    open(mesh_file, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        pts_out = output.pointlist
        println(io, "Vertices\n$nvol")
        for i in 1:nvol
            @printf(io, "%.15e %.15e %.15e 0\n", pts_out[1,i], pts_out[2,i], pts_out[3,i])
        end
        tets = output.tetrahedronlist
        println(io, "\nTetrahedra\n$ntet")
        for i in 1:ntet
            @printf(io, "%d %d %d %d 1\n", tets[1,i], tets[2,i], tets[3,i], tets[4,i])
        end
        tfaces = output.trifacelist
        tmarkers = output.trifacemarkerlist
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", tfaces[1,i], tfaces[2,i], tfaces[3,i], tmarkers[i])
        end
        println(io, "\nEnd")
    end

    return (nvol=nvol, ntet=ntet, ntri=ntri, npts_surface=npts, tess=tess)
end

# ─── Step 3: Palace config ──────────────────────────────────────────────

function write_config(filename, mesh_path, nfaces; L0=1.0)
    face_attrs = collect(1:nfaces)
    config = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "postpro/waveguide" },
  "Model": { "Mesh": "$mesh_path", "L0": $L0 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $face_attrs } },
  "Solver": {
    "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 3, "Tol": 1.0e-8, "Target": 2.0, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 }
  }
}"""
    write(filename, config)
end

# ─── Step 4: Analytical extrusion sensitivity at tessellation points ─────

function compute_extrude_sensitivity(body, tess, dir, dist, dist_dot;
                                      source_face_dot_fn=nothing)
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    status = EGADS.tess_status(tess)
    npts = status.npts

    # Find the extrusion base plane
    min_proj = Inf
    source_origin = [0.0, 0.0, 0.0]
    for f in faces
        ri = EGADS.get_range(f)
        um = (ri.range[1] + ri.range[2]) / 2
        vm = (ri.range[3] + ri.range[4]) / 2
        result = EGADS.evaluate(f, [um, vm])
        proj = dot(result.xyz, dir)
        if proj < min_proj
            min_proj = proj
            source_origin = copy(result.xyz)
        end
    end

    dxdp = zeros(Float64, npts, 3)

    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pos = fm.xyz[j, :]
            α = clamp(dot(pos - source_origin, dir) / dist, 0.0, 1.0)

            # d/d(dist): linear ramp along extrusion direction
            dxdp[ig, :] .+= (α * dist_dot) .* dir

            # Source face parameter sensitivity (same at all heights)
            if !isnothing(source_face_dot_fn)
                dxdp[ig, :] .+= source_face_dot_fn(pos)
            end
        end
    end

    return dxdp
end

# ─── Main ─────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)

    a = 3.0   # width (cm)
    b = 2.0   # height (cm)
    L = 5.0   # length (cm)
    L0 = 1e-2 # cm → meters

    c_light = 29.9792458  # speed of light in GHz·cm

    println("═══ Extruded Waveguide Sensitivity Test ═══")
    @printf("  a = %.1f cm, b = %.1f cm, L = %.1f cm\n", a, b, L)

    # Analytical eigenfrequencies
    # f_{mnp} = c/(2π) * sqrt((mπ/a)² + (nπ/b)² + (pπ/L)²)
    # but for a cavity: f = c/2 * sqrt((m/a)² + (n/b)² + (p/L)²)
    f_101 = c_light / 2 * sqrt((1/a)^2 + (0/b)^2 + (1/L)^2)
    f_011 = c_light / 2 * sqrt((0/a)^2 + (1/b)^2 + (1/L)^2)
    f_110 = c_light / 2 * sqrt((1/a)^2 + (1/b)^2 + (0/L)^2)
    @printf("  Analytical: f₁₀₁ = %.6f GHz, f₀₁₁ = %.6f GHz, f₁₁₀ = %.6f GHz\n",
            f_101, f_011, f_110)
    println()

    # Step 1: Build geometry
    println("Step 1: Building extruded waveguide...")
    ctx = EGADS.EGADSContext()
    body, src_body, dir = build_waveguide(ctx, a, b, L)
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    println("  Body: $(length(faces)) faces")

    # Step 2: Mesh
    println("\nStep 2: Tessellate + TetGen...")
    mesh_file = "mesh/waveguide.mesh"
    minfo = mesh_and_write(body, mesh_file)
    @printf("  Surface: %d pts, Volume: %d pts, %d tets\n",
            minfo.npts_surface, minfo.nvol, minfo.ntet)

    # Step 3: Palace
    println("\nStep 3: Palace eigensolve...")
    config_file = "waveguide.json"
    write_config(config_file, mesh_file, length(faces); L0=L0)

    palace = "../build/palace-build/palace-arm64.bin"
    @assert isfile(palace) "Palace binary not found"
    mkpath("postpro/waveguide")
    t0 = time()
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
    dt = time() - t0
    @printf("  Palace: %.1f sec\n", dt)

    # Read eigenfrequencies
    eig_file = "postpro/waveguide/eig.csv"
    lines = readlines(eig_file)
    freqs = Float64[]
    for line in lines
        s = strip(line)
        (startswith(s, "m,") || isempty(s)) && continue
        tokens = split(s, ",")
        push!(freqs, parse(Float64, strip(tokens[2])))
    end
    for (i, f) in enumerate(freqs)
        @printf("  Mode %d: f = %.6f GHz\n", i, f)
    end

    # Verify mesh sensitivity output
    solb_file = "postpro/waveguide/mesh_sensitivity.solb"
    mesh_file_palace = "postpro/waveguide/mesh_for_sensitivity.mesh"
    @assert isfile(solb_file) "mesh_sensitivity.solb not found"
    @assert isfile(mesh_file_palace) "mesh_for_sensitivity.mesh not found"

    # Step 4: Chain contraction
    println("\nStep 4: Sensitivity chain contraction...")
    palace_mesh = read_mesh(mesh_file_palace)
    sol = read_solb(solb_file)
    dfdx = sol.data
    nv = size(palace_mesh.vertices, 1)
    is_boundary = classify_boundary_vertices(palace_mesh)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv, count(is_boundary))

    # ─── d/dL (extrusion distance) ────────────────────────────────────────
    println("\n─── d/dL (length) ───")
    dxdp_L = compute_extrude_sensitivity(body, minfo.tess, dir, L, 1.0)

    # Map surface tess vertices to Palace mesh vertices (nearest-neighbor)
    # TetGen preserves surface vertex order, but Palace mesh may reorder
    # Use inv_evaluate on EGADS faces for the mapping
    dxdp_palace_L = zeros(Float64, nv, 3)
    n_mapped = 0
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(faces)
            uv, result = EGADS.inv_evaluate(faces[fi], pt)
            if norm(result[1:3] - pt) < 1e-6
                ig_pos = result[1:3]
                # Find the closest tess vertex
                α = clamp(dot(ig_pos - [0.0,0.0,0.0], dir) / L, 0.0, 1.0)
                # For d/dL: analytical
                dxdp_palace_L[vi, :] .= (α * 1.0) .* dir
                n_mapped += 1
                break
            end
        end
    end
    @printf("  Mapped: %d / %d boundary vertices\n", n_mapped, count(is_boundary))

    dfdL = dot(vec(dfdx), vec(dxdp_palace_L))
    nondim_L = dfdL * L / freqs[1]
    # Analytical for TE₁₀₁: df/dL * L/f = -(pπ/L)² / k² where k² = (mπ/a)²+(pπ/L)²
    k2 = (π/a)^2 + (π/L)^2
    analytical_nondim_L = -(π/L)^2 / k2
    @printf("  df/dL * L / f = %.6f  (analytical: %.6f)\n", nondim_L, analytical_nondim_L)
    @printf("  Error: %.4f\n", abs(nondim_L - analytical_nondim_L))

    # ─── d/da (width) ─────────────────────────────────────────────────────
    println("\n─── d/da (width) ───")
    # For d/da: nodes at x=a move by (1,0,0), nodes at x=0 don't move
    # Source face sensitivity: dx/da = (x/a, 0, 0) for uniform scaling in x
    # This applies at all z levels (extrusion doesn't change xy)
    dxdp_palace_a = zeros(Float64, nv, 3)
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(faces)
            uv, result = EGADS.inv_evaluate(faces[fi], pt)
            if norm(result[1:3] - pt) < 1e-6
                x = result[1]
                dxdp_palace_a[vi, :] .= [x / a, 0.0, 0.0]
                n_mapped += 1
                break
            end
        end
    end

    dfdA = dot(vec(dfdx), vec(dxdp_palace_a))
    nondim_a = dfdA * a / freqs[1]
    analytical_nondim_a = -(π/a)^2 / k2
    @printf("  df/da * a / f = %.6f  (analytical: %.6f)\n", nondim_a, analytical_nondim_a)
    @printf("  Error: %.4f\n", abs(nondim_a - analytical_nondim_a))

    # ─── d/db (height) ────────────────────────────────────────────────────
    println("\n─── d/db (height) ───")
    # For TE₁₀₁, df/db = 0 (mode doesn't depend on b)
    dxdp_palace_b = zeros(Float64, nv, 3)
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(faces)
            uv, result = EGADS.inv_evaluate(faces[fi], pt)
            if norm(result[1:3] - pt) < 1e-6
                y = result[2]
                dxdp_palace_b[vi, :] .= [0.0, y / b, 0.0]
                break
            end
        end
    end

    dfdB = dot(vec(dfdx), vec(dxdp_palace_b))
    nondim_b = dfdB * b / freqs[1]
    @printf("  df/db * b / f = %.6f  (analytical: 0.0 for TE₁₀₁)\n", nondim_b)

    # ─── Summary ──────────────────────────────────────────────────────────
    println("\n═══ Summary ═══")
    pass_L = abs(nondim_L - analytical_nondim_L) < 0.05
    pass_a = abs(nondim_a - analytical_nondim_a) < 0.05
    pass_b = abs(nondim_b) < 0.05
    @printf("  d/dL: %.4f vs %.4f — %s\n", nondim_L, analytical_nondim_L, pass_L ? "PASS ✓" : "FAIL ✗")
    @printf("  d/da: %.4f vs %.4f — %s\n", nondim_a, analytical_nondim_a, pass_a ? "PASS ✓" : "FAIL ✗")
    @printf("  d/db: %.4f vs 0.0     — %s\n", nondim_b, pass_b ? "PASS ✓" : "FAIL ✗")

    # Cleanup
    EGADS.close(ctx)

    all_pass = pass_L && pass_a && pass_b
    all_pass || @warn "Some sensitivity checks failed"
    println("\n  Done.")
end

main()
