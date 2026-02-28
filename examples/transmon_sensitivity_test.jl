# Transmon-like cavity sensitivity: re-entrant cavity with capacitor island.
#
# Builds a PEC box cavity with a central rectangular post ("capacitor island")
# via EGADS extrude + boolean subtraction. The island length (cap_length) is the
# design parameter — it IS the extrusion distance, so sensitivity is analytical.
#
# Pipeline:
#   EGADS extrude + boolean → TetGen → Palace eigensolve (MeshSensitivity)
#   → provenance-based dx/dp → Laplacian propagation → chain contraction
#   → FD Taylor test validation
#
# Geometry (cm, L0=1e-2):
#   Outer PEC box: 3 × 3 × 1.5 cm
#   Island: 0.8 × cap_length × 1.3 cm, centered in XY, on the floor
#   Boolean: cavity = box - island → re-entrant cavity
#
# Usage: julia --project=path/to/DeviceLayout.jl-egads-dev transmon_sensitivity_test.jl

using LinearAlgebra, Printf
using EGADS
using DeviceLayout
using TetGen

const DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

include(joinpath(@__DIR__, "shared", "palace_io.jl"))

# ─── EGADS geometry helpers ─────────────────────────────────────────────────

"""Extrude a rectangle in XY at z=z0 by `dist` in +Z."""
function make_box_z(ctx, x0, y0, z0, dx, dy, dist)
    n1 = EGADS.make_node(ctx, x0,      y0,      z0)
    n2 = EGADS.make_node(ctx, x0 + dx, y0,      z0)
    n3 = EGADS.make_node(ctx, x0 + dx, y0 + dy, z0)
    n4 = EGADS.make_node(ctx, x0,      y0 + dy, z0)
    e1 = EGADS.make_edge_from_nodes(ctx, n1, n2)
    e2 = EGADS.make_edge_from_nodes(ctx, n2, n3)
    e3 = EGADS.make_edge_from_nodes(ctx, n3, n4)
    e4 = EGADS.make_edge_from_nodes(ctx, n4, n1)
    loop = EGADS.make_wire_loop(ctx, [e1, e2, e3, e4])
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    body = EGADS.extrude(src, dist, [0.0, 0.0, 1.0])
    return body, src
end

"""Extrude a rectangle in XZ at y=y0 by `dist` in +Y."""
function make_box_y(ctx, x0, y0, z0, dx, dz, dist)
    n1 = EGADS.make_node(ctx, x0,      y0, z0)
    n2 = EGADS.make_node(ctx, x0 + dx, y0, z0)
    n3 = EGADS.make_node(ctx, x0 + dx, y0, z0 + dz)
    n4 = EGADS.make_node(ctx, x0,      y0, z0 + dz)
    e1 = EGADS.make_edge_from_nodes(ctx, n1, n2)
    e2 = EGADS.make_edge_from_nodes(ctx, n2, n3)
    e3 = EGADS.make_edge_from_nodes(ctx, n3, n4)
    e4 = EGADS.make_edge_from_nodes(ctx, n4, n1)
    loop = EGADS.make_wire_loop(ctx, [e1, e2, e3, e4])
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    body = EGADS.extrude(src, dist, [0.0, 1.0, 0.0])
    return body, src
end

# ─── Build geometry ─────────────────────────────────────────────────────────

const GEOM_DEFAULTS = (
    Lx = 3.0,          # cavity width (cm)
    Ly = 3.0,          # cavity depth (cm)
    Lz = 1.5,          # cavity height (cm)
    w_island = 0.8,    # island width in X (cm)
    island_h = 1.3,    # island height in Z (cm)
    island_z0 = 0.05,  # island bottom offset from floor (cm)
)

function build_geometry(ctx; cap_length=1.0, kw...)
    p = merge(GEOM_DEFAULTS, kw)
    x_c = p.Lx / 2
    y_c = p.Ly / 2
    y_start = y_c - cap_length / 2

    # Outer PEC box
    box, box_src = make_box_z(ctx, 0.0, 0.0, 0.0, p.Lx, p.Ly, p.Lz)

    # Capacitor island: extrude XZ cross-section by cap_length in +Y
    island_x = x_c - p.w_island / 2
    island, island_src = make_box_y(ctx, island_x, y_start, p.island_z0,
                                     p.w_island, p.island_h, cap_length)

    # Boolean: cavity = box - island
    result_vec = EGADS.solid_boolean(box, island, EGADS.BoolOp.SUBTRACTION)
    @assert length(result_vec) >= 1 "Boolean produced no result"
    result = result_vec[1]

    return (box=box, island=island, result=result,
            box_src=box_src, island_src=island_src)
end

# ─── Tessellate + TetGen → .mesh ───────────────────────────────────────────

function mesh_geometry(body, mesh_file; max_edge=0.25, sag=0.01, angle=15.0,
                       tetgen_opts="pq1.2a0.005")
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

    # Collect surface triangles with per-face markers
    tris = Int32[]
    tri_markers = Int32[]
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for t in 1:size(fm.tris, 1)
            for k in 1:3
                push!(tris, Int32(EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])))
            end
            push!(tri_markers, Int32(fi))
        end
    end

    # TetGen volume fill
    input = RawTetGenIO{Float64}()
    input.pointlist = pts
    facetlist!(input, reshape(tris, 3, length(tri_markers)))
    input.facetmarkerlist = tri_markers
    output = tetrahedralize(input, tetgen_opts)

    nvol = numberofpoints(output)
    ntet = numberoftetrahedra(output)
    ntri = numberoftrifaces(output)

    # Write .mesh file
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

    return (nvol=nvol, ntet=ntet, ntri=ntri, npts_surface=npts, tess=tess,
            nfaces=length(faces))
end

# ─── Palace config ──────────────────────────────────────────────────────────

function write_config(filename, mesh_path, nfaces; L0=1e-2, output_dir="postpro/transmon")
    face_attrs = collect(1:nfaces)
    config = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "$output_dir" },
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

# ─── Sensitivity chain contraction ──────────────────────────────────────────

"""
Compute df/d(cap_length) via chain contraction:
  1. Provenance-match result faces to parent bodies
  2. Analytical velocity from island extrude formula at Palace boundary vertices
  3. Laplacian propagation to interior
  4. Contract: df/dp = dot(dfdx, dxdp)
"""
function compute_chain(geo, result_body, palace_mesh, dfdx, is_boundary;
                        cap_length=1.0, kw...)
    p = merge(GEOM_DEFAULTS, kw)

    # Pre-compute provenance: result faces → parent body faces
    parent_bodies = [geo.box, geo.island]
    provenance = DLE._match_provenance(result_body, parent_bodies)

    # Pre-compute island extrusion base for the analytical velocity formula
    island_dir = [0.0, 1.0, 0.0]
    island_faces = EGADS.get_body_topos(geo.island, EGADS.OClass.FACE)
    source_origin = DLE._find_extrusion_base(island_faces, island_dir)

    # Map Palace boundary vertices to result faces and compute dx/dp
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    nv = size(palace_mesh.vertices, 1)
    dxdp_bdy = zeros(Float64, nv, 3)
    n_mapped = 0

    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            uv, eval_result = EGADS.inv_evaluate(result_faces[fi], pt)
            if norm(eval_result[1:3] - pt) < 1e-6
                # Vertex is on face fi — look up provenance
                if haskey(provenance, fi) && provenance[fi].body_idx == 2
                    # Island face: analytical extrude velocity
                    α = clamp(dot(pt - source_origin, island_dir) / cap_length, 0.0, 1.0)
                    dxdp_bdy[vi, :] .= (α * 1.0) .* island_dir  # d/d(cap_length)
                end
                # Box faces (body_idx==1) or unmatched: zero velocity (default)
                n_mapped += 1
                break
            end
        end
    end

    # Propagate boundary sensitivity to interior via Laplacian smoothing
    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)

    # Chain contraction
    dfdp = dot(vec(dfdx), vec(dxdp_full))

    return (dfdp=dfdp, n_mapped=n_mapped, dxdp_bdy=dxdp_bdy, dxdp_full=dxdp_full)
end

# ─── Run a single Palace evaluation ────────────────────────────────────────

function run_palace_eval(ctx, cap_length; label="baseline", palace_bin=nothing)
    if isnothing(palace_bin)
        palace_bin = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    end
    @assert isfile(palace_bin) "Palace binary not found: $palace_bin"

    geo = build_geometry(ctx; cap_length=cap_length)

    mesh_file = "mesh/transmon_$(label).mesh"
    minfo = mesh_geometry(geo.result, mesh_file)

    config_file = "transmon_$(label).json"
    outdir = "postpro/transmon_$(label)"
    write_config(config_file, mesh_file, minfo.nfaces; output_dir=outdir)

    mkpath(outdir)
    run(pipeline(`$palace_bin $config_file`; stdout=devnull, stderr=devnull))

    freq = read_eigenfrequency(joinpath(outdir, "eig.csv"); mode=1)
    return (geo=geo, minfo=minfo, freq=freq, outdir=outdir, config_file=config_file,
            mesh_file=mesh_file)
end

# ─── Main ───────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    cap_length = 1.0  # cm (design parameter)

    println("═══ Transmon Cavity Sensitivity Test ═══")
    println("  Geometry: 3×3×1.5 cm PEC box, 0.8×$(cap_length)×1.3 cm island")
    println()

    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found"

    # ─── Step 1: Build, mesh, solve at baseline ────────────────────────────

    println("Step 1: Baseline geometry + Palace...")
    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)
    baseline = run_palace_eval(ctx, cap_length; label="baseline", palace_bin=palace)
    EGADS.close(ctx)

    @printf("  Mesh: %d verts, %d tets, %d faces\n",
            baseline.minfo.nvol, baseline.minfo.ntet, baseline.minfo.nfaces)
    @printf("  Mode 1: f = %.6f GHz\n", baseline.freq)

    # Read mesh sensitivity output
    solb_file = joinpath(baseline.outdir, "mesh_sensitivity.solb")
    mesh_file_palace = joinpath(baseline.outdir, "mesh_for_sensitivity.mesh")
    @assert isfile(solb_file) "mesh_sensitivity.solb not found"
    @assert isfile(mesh_file_palace) "mesh_for_sensitivity.mesh not found"

    palace_mesh = read_mesh(mesh_file_palace)
    sol = read_solb(solb_file)
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv, count(is_boundary))

    # ─── Step 2: Chain contraction ─────────────────────────────────────────

    println("\nStep 2: Sensitivity chain contraction...")
    ctx2 = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx2, 0)
    geo2 = build_geometry(ctx2; cap_length=cap_length)

    chain = compute_chain(geo2, geo2.result, palace_mesh, dfdx, is_boundary;
                           cap_length=cap_length)
    EGADS.close(ctx2)

    @printf("  Mapped: %d / %d boundary vertices\n", chain.n_mapped, count(is_boundary))
    @printf("  df/d(cap_length) = %.6e GHz/cm\n", chain.dfdp)

    # ─── Step 3: Mesh-perturbation Taylor test ───────────────────────────
    # Perturb the SAME mesh vertices by ε * dx/dp, re-run Palace.
    # This eliminates TetGen re-meshing noise that contaminates geometry-rebuild FD.

    println("\nStep 3: Mesh-perturbation Taylor test...")
    h_values = [1e-1, 5e-2, 2e-2, 1e-2]
    f0 = baseline.freq
    dfdp_chain = chain.dfdp

    println("  Taylor test: |f(p+ε) - f(p) - ε·df/dp| / ε²")
    taylor_errs = Float64[]

    for (i, h) in enumerate(h_values)
        # Write perturbed mesh: x_new = x_old + ε * dxdp
        perturbed_mesh_file = "mesh/transmon_taylor_$(i).mesh"
        mkpath("mesh")
        open(perturbed_mesh_file, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv")
            for vi in 1:nv
                x = palace_mesh.vertices[vi, 1] + h * chain.dxdp_full[vi, 1]
                y = palace_mesh.vertices[vi, 2] + h * chain.dxdp_full[vi, 2]
                z = palace_mesh.vertices[vi, 3] + h * chain.dxdp_full[vi, 3]
                @printf(io, "%.15e %.15e %.15e 0\n", x, y, z)
            end
            ntet_p = size(palace_mesh.tetrahedra, 1)
            if ntet_p > 0
                println(io, "\nTetrahedra\n$ntet_p")
                for e in 1:ntet_p
                    ref = palace_mesh.tet_refs[e]
                    @printf(io, "%d %d %d %d %d\n",
                            palace_mesh.tetrahedra[e,1], palace_mesh.tetrahedra[e,2],
                            palace_mesh.tetrahedra[e,3], palace_mesh.tetrahedra[e,4], ref)
                end
            end
            nhex = size(palace_mesh.hexahedra, 1)
            if nhex > 0
                println(io, "\nHexahedra\n$nhex")
                for e in 1:nhex
                    ref = palace_mesh.hex_refs[e]
                    @printf(io, "%d %d %d %d %d %d %d %d %d\n",
                            palace_mesh.hexahedra[e,1], palace_mesh.hexahedra[e,2],
                            palace_mesh.hexahedra[e,3], palace_mesh.hexahedra[e,4],
                            palace_mesh.hexahedra[e,5], palace_mesh.hexahedra[e,6],
                            palace_mesh.hexahedra[e,7], palace_mesh.hexahedra[e,8], ref)
                end
            end
            ntri_p = size(palace_mesh.triangles, 1)
            if ntri_p > 0
                println(io, "\nTriangles\n$ntri_p")
                for e in 1:ntri_p
                    ref = palace_mesh.tri_refs[e]
                    @printf(io, "%d %d %d %d\n",
                            palace_mesh.triangles[e,1], palace_mesh.triangles[e,2],
                            palace_mesh.triangles[e,3], ref)
                end
            end
            nquad = size(palace_mesh.quadrilaterals, 1)
            if nquad > 0
                println(io, "\nQuadrilaterals\n$nquad")
                for e in 1:nquad
                    ref = palace_mesh.quad_refs[e]
                    @printf(io, "%d %d %d %d %d\n",
                            palace_mesh.quadrilaterals[e,1], palace_mesh.quadrilaterals[e,2],
                            palace_mesh.quadrilaterals[e,3], palace_mesh.quadrilaterals[e,4], ref)
                end
            end
            println(io, "\nEnd")
        end

        # Write config pointing to perturbed mesh
        perturbed_config = "transmon_taylor_$(i).json"
        perturbed_outdir = "postpro/transmon_taylor_$(i)"
        write_config(perturbed_config, perturbed_mesh_file, baseline.minfo.nfaces;
                     output_dir=perturbed_outdir)
        mkpath(perturbed_outdir)

        # Run Palace on perturbed mesh
        run(pipeline(`$palace $perturbed_config`; stdout=devnull, stderr=devnull))
        f_h = read_eigenfrequency(joinpath(perturbed_outdir, "eig.csv"); mode=1)

        linear_pred = f0 + h * dfdp_chain
        taylor_err = abs(f_h - linear_pred)
        push!(taylor_errs, taylor_err)
        @printf("  ε=%.0e: f(p+ε)=%.8f, pred=%.8f, err=%.2e, err/ε²=%.2e\n",
                h, f_h, linear_pred, taylor_err, taylor_err / h^2)

        if i > 1
            rate = log(taylor_errs[i-1] / taylor_errs[i]) / log(h_values[i-1] / h_values[i])
            @printf("  → rate = %.2f\n", rate)
        end

        # Clean up
        rm(perturbed_outdir; recursive=true, force=true)
        rm(perturbed_config; force=true)
        rm(perturbed_mesh_file; force=true)
    end

    # Compute convergence rates
    rates = Float64[]
    for i in 2:length(taylor_errs)
        rate = log(taylor_errs[i-1] / taylor_errs[i]) / log(h_values[i-1] / h_values[i])
        push!(rates, rate)
    end

    # ─── Summary ───────────────────────────────────────────────────────────

    println("\n═══ Summary ═══")
    @printf("  f₀ = %.6f GHz at cap_length = %.2f cm\n", f0, cap_length)
    @printf("  df/d(cap_length) = %.6e GHz/cm\n", dfdp_chain)
    if !isempty(rates)
        @printf("  Taylor rates = %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
        max_rate = maximum(rates)
        pass = max_rate > 1.5
        @printf("  Status: %s (max rate %.2f, need > 1.5 for O(ε²))\n",
                pass ? "PASS ✓" : "FAIL ✗", max_rate)
    end

    # Clean up baseline
    rm(baseline.outdir; recursive=true, force=true)
    rm(baseline.config_file; force=true)
    rm(baseline.mesh_file; force=true)

    println("\nDone.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
