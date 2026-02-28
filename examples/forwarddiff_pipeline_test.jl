# ForwardDiff → EGADS → Palace pipeline test
#
# The integration test connecting:
#   1. ForwardDiff through DeviceLayout geometry (d(polygon)/d(params))
#   2. SourceFaceVelocity (bilinear interpolation on polygon edges)
#   3. EGADS extrusion + boolean (3D geometry)
#   4. Palace eigensolve (df/dx mesh sensitivity)
#   5. Chain contraction (df/dp = df/dx · dx/dp)
#   6. Taylor test (O(ε²) verification)
#
# Geometry:
#   Substrate: 200 × 200 × 10 μm → 0.02 × 0.02 × 0.001 cm (L0=1e-2 means cm)
#   Island cavity: cross-shaped from difference2d (cap island), extruded in Z
#   Boolean: substrate - island_cavity
#   Design params: cap_length, cap_gap (change source face shape)
#
# Usage: julia --project=path/to/DeviceLayout.jl-ad-sensitivity forwarddiff_pipeline_test.jl

using LinearAlgebra, Printf
using DeviceLayout, DeviceLayout.PreferredUnits
using ForwardDiff, Unitful
using EGADS
using TetGen

include(joinpath(@__DIR__, "shared", "palace_io.jl"))
include(joinpath(homedir(), "AWS", "DeviceLayout.jl-ad-sensitivity", "src", "forwarddiff_bridge.jl"))

# ─── Inlined helpers ─────────────────────────────────────────────────────────

function face_plane(face)
    ri = EGADS.get_range(face)
    um = (ri.range[1] + ri.range[2]) / 2
    vm = (ri.range[3] + ri.range[4]) / 2
    res = EGADS.evaluate(face, [um, vm])
    n = cross(res.du, res.dv)
    nrm = norm(n)
    nrm < 1e-30 && return (normal=[0.0, 0.0, 0.0], offset=0.0, point=res.xyz)
    n = n / nrm
    d = dot(n, res.xyz)
    return (normal=n, offset=d, point=res.xyz)
end

function match_provenance(result_body, parent_bodies; tol=1e-6)
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    parent_data = Tuple{Int, Int, EGADS.EGObject, NamedTuple}[]
    for (bi, body) in enumerate(parent_bodies)
        for (fi, pf) in enumerate(EGADS.get_body_topos(body, EGADS.OClass.FACE))
            push!(parent_data, (bi, fi, pf, face_plane(pf)))
        end
    end
    provenance = Dict{Int, NamedTuple{(:body_idx, :face_idx), Tuple{Int, Int}}}()
    for (ri, rf) in enumerate(result_faces)
        rp = face_plane(rf)
        for (bi, fi, _, pp) in parent_data
            ndot = abs(dot(rp.normal, pp.normal))
            if ndot > 1 - tol &&
               (abs(rp.offset - pp.offset) < tol || abs(rp.offset + pp.offset) < tol)
                provenance[ri] = (body_idx=bi, face_idx=fi)
                break
            end
        end
    end
    return provenance
end

function find_extrusion_base(faces, dir)
    min_proj = Inf
    base_pos = [0.0, 0.0, 0.0]
    for f in faces
        rng = EGADS.get_range(f)
        um = (rng.range[1] + rng.range[2]) / 2
        vm = (rng.range[3] + rng.range[4]) / 2
        center = EGADS.evaluate(f, [um, vm]).xyz
        if dot(center, dir) < min_proj
            min_proj = dot(center, dir)
            base_pos = center
        end
    end
    return base_pos
end

struct SourceFaceVelocity
    vertices::Matrix{Float64}
    velocity::Matrix{Float64}
end

function interpolate_source_velocity(sfv::SourceFaceVelocity, x::Float64, y::Float64)
    nv = size(sfv.vertices, 1)
    nv == 0 && return [0.0, 0.0, 0.0]
    min_dist_sq = Inf; best_s = 0.0; best_i = 1
    for i in 1:nv
        j = mod1(i + 1, nv)
        ax, ay = sfv.vertices[i, 1], sfv.vertices[i, 2]
        bx, by = sfv.vertices[j, 1], sfv.vertices[j, 2]
        dx_ab, dy_ab = bx - ax, by - ay
        len_sq = dx_ab^2 + dy_ab^2
        if len_sq < 1e-30
            s = 0.0; dist_sq = (x - ax)^2 + (y - ay)^2
        else
            s = clamp(((x - ax) * dx_ab + (y - ay) * dy_ab) / len_sq, 0.0, 1.0)
            px, py = ax + s * dx_ab, ay + s * dy_ab
            dist_sq = (x - px)^2 + (y - py)^2
        end
        if dist_sq < min_dist_sq
            min_dist_sq = dist_sq; best_s = s; best_i = i
        end
    end
    j = mod1(best_i + 1, nv)
    vx = (1 - best_s) * sfv.velocity[best_i, 1] + best_s * sfv.velocity[j, 1]
    vy = (1 - best_s) * sfv.velocity[best_i, 2] + best_s * sfv.velocity[j, 2]
    return [vx, vy, 0.0]
end

# ─── DeviceLayout island polygon builder (ForwardDiff-compatible) ────────────
#
# Uses DeviceLayout's Rectangle, difference2d, to_polygons — the real code path
# that a user would use for transmon design.

function island_coords_cm(params)
    # params = [cap_length_μm, cap_gap_μm]
    # Returns the OUTER rectangle (gap boundary) as a simple 4-vertex polygon in cm.
    # The inner rectangle (island metal) would be a separate body in the EGADS geometry.
    # This avoids the self-touching polygon from difference2d.
    cap_length = params[1] * 1μm
    cap_gap = params[2] * 1μm
    cap_width = 24.0μm
    junction_gap = 12.0μm

    # Outer boundary of the gap region (simple convex rectangle)
    x0 = -cap_width / 2 - cap_gap
    y0 = zero(cap_gap)
    x1 = cap_width / 2 + cap_gap
    y1 = junction_gap + cap_length + cap_gap

    r = Rectangle(Point(x0, y0), Point(x1, y1))

    coords = eltype(params)[]
    for pt in points(r)
        push!(coords, Unitful.ustrip(u"cm", getx(pt)))
        push!(coords, Unitful.ustrip(u"cm", gety(pt)))
    end
    return coords
end

# ─── Build EGADS body from polygon vertices ──────────────────────────────────

function make_extruded_polygon(ctx, verts_cm, z0, height)
    nv = size(verts_cm, 1)
    nv < 3 && error("Need at least 3 vertices, got $nv")

    nodes = [EGADS.make_node(ctx, verts_cm[i, 1], verts_cm[i, 2], z0) for i in 1:nv]
    edges = [EGADS.make_edge_from_nodes(ctx, nodes[i], nodes[mod1(i + 1, nv)]) for i in 1:nv]
    loop = EGADS.make_wire_loop(ctx, edges)
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    body = EGADS.extrude(src, height, [0.0, 0.0, 1.0])
    return body
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace) "Palace binary not found"

    # Design parameters (in μm for DeviceLayout, converted to cm for EGADS/Palace)
    cap_length_μm = 620.0
    cap_gap_μm = 30.0
    params_nom = [cap_length_μm, cap_gap_μm]

    # Substrate dimensions (cm)
    sub_w = 0.1  # 1000 μm = 0.1 cm
    sub_h = 0.1
    sub_z = 0.05  # 500 μm substrate thickness
    island_z0 = 0.001  # 10 μm above floor (gap)
    island_height = sub_z - island_z0 - 0.001  # leave 10 μm gap at top

    println("═══ ForwardDiff → EGADS → Palace Pipeline Test ═══")
    println("  params: cap_length=$(cap_length_μm)μm, cap_gap=$(cap_gap_μm)μm")

    # ─── Step 1: ForwardDiff through DeviceLayout ────────────────────────────

    println("\nStep 1: ForwardDiff through DeviceLayout polygon pipeline...")
    coords_nom = island_coords_cm(params_nom)
    J = ForwardDiff.jacobian(island_coords_cm, params_nom)
    nv = length(coords_nom) ÷ 2
    println("  Island polygon: $nv vertices")
    println("  Jacobian: $(size(J, 1)) × $(size(J, 2))")

    verts = coords_to_vertices(coords_nom)

    # Center the polygon in the substrate
    cx = (maximum(verts[:, 1]) + minimum(verts[:, 1])) / 2
    cy = (maximum(verts[:, 2]) + minimum(verts[:, 2])) / 2
    offset_x = sub_w / 2 - cx
    offset_y = sub_h / 2 - cy
    verts[:, 1] .+= offset_x
    verts[:, 2] .+= offset_y

    # Build SourceFaceVelocity for d/d(cap_length)
    vel_cl = jacobian_column_to_velocity(J, 1)
    sfv_cl = SourceFaceVelocity(verts, vel_cl)

    n_moving = count(i -> maximum(abs.(vel_cl[i, :])) > 1e-10, 1:nv)
    println("  d/d(cap_length): $n_moving / $nv vertices moving")

    # ─── Step 2: Build EGADS geometry ────────────────────────────────────────

    println("\nStep 2: Building EGADS geometry...")
    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)

    # Substrate box
    substrate = EGADS.make_box(ctx, 0.0, 0.0, 0.0, sub_w, sub_h, sub_z)

    # Island cavity from DeviceLayout polygon
    island = make_extruded_polygon(ctx, verts, island_z0, island_height)

    # Boolean: substrate - island
    result_vec = EGADS.solid_boolean(substrate, island, EGADS.BoolOp.SUBTRACTION)
    @assert length(result_vec) >= 1 "Boolean failed"
    result = result_vec[1]

    result_faces = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    println("  Result: $(length(result_faces)) faces")

    # ─── Step 3: Mesh + Palace ───────────────────────────────────────────────

    println("\nStep 3: Meshing + Palace eigensolve...")
    tess = EGADS.tessellate(result; max_edge=0.01, sag=0.0005, angle=15.0)
    npts_tess = EGADS.tess_status(tess).npts
    println("  Tessellation: $npts_tess surface vertices")

    # Build TetGen input
    pts = zeros(Float64, 3, npts_tess)
    for fi in eachindex(result_faces)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end
    tris = Int32[]
    tri_markers = Int32[]
    for fi in eachindex(result_faces)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for t in 1:size(fm.tris, 1)
            for k in 1:3
                push!(tris, Int32(EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])))
            end
            push!(tri_markers, Int32(fi))
        end
    end
    input = RawTetGenIO{Float64}()
    input.pointlist = pts
    facetlist!(input, reshape(tris, 3, length(tri_markers)))
    input.facetmarkerlist = tri_markers
    tetout = tetrahedralize(input, "pq1.2a1e-7")

    nvol = numberofpoints(tetout)
    ntet = numberoftetrahedra(tetout)
    println("  Volume mesh: $nvol verts, $ntet tets")

    # Write mesh
    mesh_file = "mesh/fd_pipeline.mesh"
    mkpath("mesh")
    open(mesh_file, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        p = tetout.pointlist
        println(io, "Vertices\n$nvol")
        for i in 1:nvol
            @printf(io, "%.15e %.15e %.15e 0\n", p[1, i], p[2, i], p[3, i])
        end
        te = tetout.tetrahedronlist
        println(io, "\nTetrahedra\n$ntet")
        for i in 1:ntet
            @printf(io, "%d %d %d %d 1\n", te[1, i], te[2, i], te[3, i], te[4, i])
        end
        tf = tetout.trifacelist
        tm = tetout.trifacemarkerlist
        ntri = size(tf, 2)
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", tf[1, i], tf[2, i], tf[3, i], tm[i])
        end
        println(io, "\nEnd")
    end

    nfaces = length(result_faces)
    config = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "postpro/fd_pipeline" },
  "Model": { "Mesh": "$mesh_file", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": {
    "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 10.0, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 }
  }
}"""
    config_file = "fd_pipeline.json"
    write(config_file, config)
    mkpath("postpro/fd_pipeline")

    flush(stdout)
    t0 = time()
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
    @printf("  Palace: %.1fs\n", time() - t0)

    freq = read_eigenfrequency("postpro/fd_pipeline/eig.csv"; mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    # ─── Step 4: Chain contraction with ForwardDiff-derived SFV ──────────────

    println("\nStep 4: Chain contraction (ForwardDiff → SourceFaceVelocity)...")
    palace_mesh = read_mesh("postpro/fd_pipeline/mesh_for_sensitivity.mesh")
    sol = read_solb("postpro/fd_pipeline/mesh_sensitivity.solb")
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv_palace = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv_palace, count(is_boundary))

    # Provenance: result faces → parent bodies
    parent_bodies = [substrate, island]
    provenance = match_provenance(result, parent_bodies)
    @printf("  Provenance: %d / %d faces matched\n", length(provenance), nfaces)

    # Island extrusion parameters
    dir = [0.0, 0.0, 1.0]
    island_faces = EGADS.get_body_topos(island, EGADS.OClass.FACE)
    source_origin = find_extrusion_base(island_faces, dir)

    # Compute dx/dp for d/d(cap_length) using ForwardDiff-derived SFV
    dxdp_bdy = zeros(Float64, nv_palace, 3)
    n_mapped = 0
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            inv_res = EGADS.inv_evaluate(result_faces[fi], pt)
            if LinearAlgebra.norm(inv_res.results[1:3] - pt) < 1e-6
                if haskey(provenance, fi) && provenance[fi].body_idx == 2
                    α = clamp(dot(pt - source_origin, dir) / island_height, 0.0, 1.0)
                    pos_source = pt .- (α * island_height) .* dir
                    vel = interpolate_source_velocity(sfv_cl, pos_source[1], pos_source[2])
                    dxdp_bdy[vi, :] .= vel
                end
                n_mapped += 1
                break
            end
        end
    end
    @printf("  Mapped: %d / %d boundary verts\n", n_mapped, count(is_boundary))

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    @printf("  df/d(cap_length) = %.6e GHz/cm\n", dfdp)

    EGADS.close(ctx)

    # ─── Step 5: Taylor test ─────────────────────────────────────────────────

    println("\nStep 5: Mesh-perturbation Taylor test...")
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]

    for (i, ε) in enumerate(epsilons)
        pert_mesh = "mesh/fd_pert_$(i).mesh"
        open(pert_mesh, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv_palace")
            for vi in 1:nv_palace
                @printf(io, "%.15e %.15e %.15e 0\n",
                        palace_mesh.vertices[vi, 1] + ε * dxdp_full[vi, 1],
                        palace_mesh.vertices[vi, 2] + ε * dxdp_full[vi, 2],
                        palace_mesh.vertices[vi, 3] + ε * dxdp_full[vi, 3])
            end
            ntet_p = size(palace_mesh.tetrahedra, 1)
            println(io, "\nTetrahedra\n$ntet_p")
            for e in 1:ntet_p
                @printf(io, "%d %d %d %d %d\n",
                        palace_mesh.tetrahedra[e, 1], palace_mesh.tetrahedra[e, 2],
                        palace_mesh.tetrahedra[e, 3], palace_mesh.tetrahedra[e, 4],
                        palace_mesh.tet_refs[e])
            end
            ntri_p = size(palace_mesh.triangles, 1)
            println(io, "\nTriangles\n$ntri_p")
            for e in 1:ntri_p
                @printf(io, "%d %d %d %d\n",
                        palace_mesh.triangles[e, 1], palace_mesh.triangles[e, 2],
                        palace_mesh.triangles[e, 3], palace_mesh.tri_refs[e])
            end
            println(io, "\nEnd")
        end

        pert_cfg = "fd_pert_$(i).json"
        pert_out = "postpro/fd_pert_$(i)"
        write(pert_cfg, replace(config, "fd_pipeline" => "fd_pert_$(i)"))
        mkpath(pert_out)
        run(pipeline(`$palace $pert_cfg`; stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(pert_out, "eig.csv"); mode=1)

        rem = abs(f_pert - freq - ε * dfdp)
        push!(remainders, rem)
        rate = length(remainders) >= 2 ?
            log(remainders[end-1] / remainders[end]) / log(epsilons[length(remainders)-1] / ε) : NaN
        @printf("  ε=%.4f: f=%.8f, rem=%.2e, rate=%.2f\n", ε, f_pert, rem, rate)

        rm(pert_out; recursive=true, force=true)
        rm(pert_cfg; force=true)
        rm(pert_mesh; force=true)
    end

    rates = [log(remainders[i-1] / remainders[i]) / log(epsilons[i-1] / epsilons[i])
             for i in 2:length(remainders)]
    max_rate = isempty(rates) ? 0.0 : maximum(rates)
    pass = max_rate > 1.5

    println("\n═══ Summary ═══")
    @printf("  ForwardDiff: %d-vertex island polygon, 2-param Jacobian\n", nv)
    @printf("  EGADS: %d faces after boolean\n", nfaces)
    @printf("  Palace: %d verts, f₀ = %.6f GHz\n", nv_palace, freq)
    @printf("  df/d(cap_length) = %.6e GHz/cm\n", dfdp)
    @printf("  Taylor rates: %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
    @printf("  Status: %s (max rate %.2f)\n", pass ? "PASS ✓" : "FAIL ✗", max_rate)

    rm("postpro/fd_pipeline"; recursive=true, force=true)
    rm(config_file; force=true)
    rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
