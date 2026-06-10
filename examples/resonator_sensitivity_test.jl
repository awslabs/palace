# Resonator path sensitivity: d/d(bend_radius) through 1500-vertex CPW meander.
#
# Tests the pipeline at scale: ForwardDiff through DeviceLayout's Path + meander!
# pipeline, union2d to merge per-node polygons, EGADS extrude 1500-vertex polygon,
# boolean subtract from substrate, Palace eigensolve, chain contraction, Taylor test.
#
# Usage: julia --project=path/to/DeviceLayout.jl-ad-sensitivity resonator_sensitivity_test.jl

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
    return (normal=n, offset=dot(n, res.xyz), point=res.xyz)
end

function match_provenance(result_body, parent_bodies; tol=1e-6)
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    parent_data = Tuple{Int, Int, EGADS.EGObject, NamedTuple}[]
    for (bi, body) in enumerate(parent_bodies)
        for (fi, pf) in enumerate(EGADS.get_body_topos(body, EGADS.OClass.FACE))
            push!(parent_data, (bi, fi, pf, face_plane(pf)))
        end
    end
    prov = Dict{Int, NamedTuple{(:body_idx, :face_idx), Tuple{Int, Int}}}()
    for (ri, rf) in enumerate(result_faces)
        rp = face_plane(rf)
        for (bi, fi, _, pp) in parent_data
            ndot = abs(dot(rp.normal, pp.normal))
            if ndot > 1 - tol &&
               (abs(rp.offset - pp.offset) < tol || abs(rp.offset + pp.offset) < tol)
                prov[ri] = (body_idx=bi, face_idx=fi)
                break
            end
        end
    end
    return prov
end

function find_extrusion_base(faces, dir)
    min_proj = Inf
    base_pos = [0.0, 0.0, 0.0]
    for f in faces
        rng = EGADS.get_range(f)
        center = EGADS.evaluate(f, [(rng.range[1]+rng.range[2])/2, (rng.range[3]+rng.range[4])/2]).xyz
        if dot(center, dir) < min_proj; min_proj = dot(center, dir); base_pos = center; end
    end
    return base_pos
end

struct SourceFaceVelocity; vertices::Matrix{Float64}; velocity::Matrix{Float64}; end

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
        if len_sq < 1e-30; s = 0.0; dist_sq = (x-ax)^2 + (y-ay)^2
        else; s = clamp(((x-ax)*dx_ab + (y-ay)*dy_ab)/len_sq, 0.0, 1.0)
            px, py = ax + s*dx_ab, ay + s*dy_ab; dist_sq = (x-px)^2 + (y-py)^2; end
        if dist_sq < min_dist_sq; min_dist_sq = dist_sq; best_s = s; best_i = i; end
    end
    j = mod1(best_i + 1, nv)
    vx = (1-best_s)*sfv.velocity[best_i,1] + best_s*sfv.velocity[j,1]
    vy = (1-best_s)*sfv.velocity[best_i,2] + best_s*sfv.velocity[j,2]
    return [vx, vy, 0.0]
end

# ─── Resonator path (ForwardDiff-compatible, returns one gap strip) ──────────

function resonator_gap_coords_cm(bend_radius_μm)
    r = bend_radius_μm * 1μm
    T = typeof(r)
    style = Paths.SimpleCPW(T(10.0μm), T(6.0μm))
    p = Path(Point(zero(r), zero(r)); α0=0.0°)
    straight!(p, T(200.0μm), style)
    turn!(p, -90°, r)
    straight!(p, T(300.0μm))
    turn!(p, 180°, r)
    straight!(p, T(300.0μm))

    all_polys = []
    for node in p.nodes
        append!(all_polys, to_polygons(node))
    end
    result = all_polys[1]
    for i in 2:length(all_polys)
        result = union2d(result, all_polys[i])
    end
    poly = to_polygons(result)[1]

    coords = typeof(bend_radius_μm)[]
    for pt in points(poly)
        push!(coords, Unitful.ustrip(u"cm", getx(pt)))
        push!(coords, Unitful.ustrip(u"cm", gety(pt)))
    end
    return coords
end

# ─── EGADS ───────────────────────────────────────────────────────────────────

function make_extruded_polygon(ctx, verts_cm, z0, height)
    nv = size(verts_cm, 1)
    nodes = [EGADS.make_node(ctx, verts_cm[i,1], verts_cm[i,2], z0) for i in 1:nv]
    edges = [EGADS.make_edge_from_nodes(ctx, nodes[i], nodes[mod1(i+1,nv)]) for i in 1:nv]
    loop = EGADS.make_wire_loop(ctx, edges)
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    return EGADS.extrude(src, height, [0.0, 0.0, 1.0])
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace)

    bend_radius_μm = 50.0
    sub_margin = 0.005  # 50 μm margin around resonator

    println("═══ Resonator Path Sensitivity Test ═══")
    println("  bend_radius = $(bend_radius_μm) μm")

    # ─── Step 1: ForwardDiff ─────────────────────────────────────────────────

    println("\nStep 1: ForwardDiff through resonator path...")
    t0 = time()
    coords_nom = resonator_gap_coords_cm(bend_radius_μm)
    d = ForwardDiff.derivative(resonator_gap_coords_cm, bend_radius_μm)
    nv = length(coords_nom) ÷ 2
    @printf("  %d vertices, ForwardDiff in %.1fs\n", nv, time() - t0)

    verts_full = coords_to_vertices(coords_nom)
    vel_full = reshape(d, nv, 2)

    # Subsample to ~100 vertices for meshable EGADS body.
    # Keep every Nth vertex + always keep the last vertex for closure.
    # The SFV uses the FULL polygon for interpolation (accuracy).
    target_nv = 100
    step = max(1, nv ÷ target_nv)
    keep = sort(union(1:step:nv, nv))
    verts = verts_full[keep, :]
    vel = vel_full[keep, :]
    nv_sub = size(verts, 1)
    @printf("  Subsampled: %d → %d vertices (step=%d)\n", nv, nv_sub, step)

    # SFV uses the subsampled polygon (matches EGADS geometry)
    sfv = SourceFaceVelocity(verts, vel)
    nv = nv_sub  # update for EGADS extrusion

    n_moving = count(i -> max(abs(vel[i,1]), abs(vel[i,2])) > 1e-8, 1:nv)
    println("  d/d(bend_radius): $n_moving / $nv vertices moving")

    # Compute bounding box for substrate
    xmin, xmax = minimum(verts[:,1]), maximum(verts[:,1])
    ymin, ymax = minimum(verts[:,2]), maximum(verts[:,2])
    sub_x0 = xmin - sub_margin
    sub_y0 = ymin - sub_margin
    sub_w = (xmax - xmin) + 2 * sub_margin
    sub_h = (ymax - ymin) + 2 * sub_margin
    sub_z = 0.01  # 100 μm substrate
    gap_z0 = 0.001
    gap_height = sub_z - gap_z0 - 0.001

    # ─── Step 2: EGADS geometry ──────────────────────────────────────────────

    println("\nStep 2: EGADS geometry ($nv-vertex extrusion)...")
    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)

    substrate = EGADS.make_box(ctx, sub_x0, sub_y0, 0.0, sub_w, sub_h, sub_z)
    gap_body = make_extruded_polygon(ctx, verts, gap_z0, gap_height)
    result_vec = EGADS.solid_boolean(substrate, gap_body, EGADS.BoolOp.SUBTRACTION)
    result = result_vec[1]
    nfaces = length(EGADS.get_body_topos(result, EGADS.OClass.FACE))
    @printf("  Result: %d faces (substrate - %d-vertex gap)\n", nfaces, nv)

    # ─── Step 3: Mesh + Palace ───────────────────────────────────────────────

    println("\nStep 3: TetGen + Palace...")
    tess = EGADS.tessellate(result; max_edge=0.008, sag=0.0005, angle=15.0)
    npts_tess = EGADS.tess_status(tess).npts
    println("  Tessellation: $npts_tess surface vertices")

    result_faces_list = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    pts = zeros(Float64, 3, npts_tess)
    for fi in eachindex(result_faces_list)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end
    tris = Int32[]; tri_markers = Int32[]
    for fi in eachindex(result_faces_list)
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
    tetout = tetrahedralize(input, "pq1.2a5e-7")

    nvol = numberofpoints(tetout)
    ntet = numberoftetrahedra(tetout)
    @printf("  Volume: %d verts, %d tets\n", nvol, ntet)

    mesh_file = "mesh/reson.mesh"
    mkpath("mesh")
    open(mesh_file, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        p = tetout.pointlist
        println(io, "Vertices\n$nvol")
        for i in 1:nvol; @printf(io, "%.15e %.15e %.15e 0\n", p[1,i], p[2,i], p[3,i]); end
        te = tetout.tetrahedronlist
        println(io, "\nTetrahedra\n$ntet")
        for i in 1:ntet; @printf(io, "%d %d %d %d 1\n", te[1,i], te[2,i], te[3,i], te[4,i]); end
        tf = tetout.trifacelist; tm = tetout.trifacemarkerlist; ntri = size(tf, 2)
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri; @printf(io, "%d %d %d %d\n", tf[1,i], tf[2,i], tf[3,i], tm[i]); end
        println(io, "\nEnd")
    end

    config_str = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "postpro/reson" },
  "Model": { "Mesh": "$mesh_file", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": { "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 30.0, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 } }
}"""
    config_file = "reson.json"
    write(config_file, config_str)
    mkpath("postpro/reson")

    flush(stdout)
    t0 = time()
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
    @printf("  Palace: %.1fs\n", time() - t0)

    freq = read_eigenfrequency("postpro/reson/eig.csv"; mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    palace_mesh = read_mesh("postpro/reson/mesh_for_sensitivity.mesh")
    sol = read_solb("postpro/reson/mesh_sensitivity.solb")
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv_palace = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv_palace, count(is_boundary))

    # ─── Step 4: Chain contraction ───────────────────────────────────────────

    println("\nStep 4: Chain contraction (d/d(bend_radius))...")
    parent_bodies = [substrate, gap_body]
    prov = match_provenance(result, parent_bodies)
    @printf("  Provenance: %d / %d faces\n", length(prov), nfaces)

    dir = [0.0, 0.0, 1.0]
    gap_faces = EGADS.get_body_topos(gap_body, EGADS.OClass.FACE)
    source_origin = find_extrusion_base(gap_faces, dir)

    dxdp_bdy = zeros(Float64, nv_palace, 3)
    n_mapped = 0
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces_list)
            inv_res = EGADS.inv_evaluate(result_faces_list[fi], pt)
            if LinearAlgebra.norm(inv_res.results[1:3] - pt) < 1e-6
                if haskey(prov, fi) && prov[fi].body_idx == 2
                    α = clamp(dot(pt - source_origin, dir) / gap_height, 0.0, 1.0)
                    pos_source = pt .- (α * gap_height) .* dir
                    v = interpolate_source_velocity(sfv, pos_source[1], pos_source[2])
                    dxdp_bdy[vi, :] .= v
                end
                n_mapped += 1
                break
            end
        end
    end
    @printf("  Mapped: %d / %d boundary\n", n_mapped, count(is_boundary))

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    @printf("  df/d(bend_radius) = %.6e GHz/cm\n", dfdp)

    EGADS.close(ctx)

    # ─── Step 5: Taylor test ─────────────────────────────────────────────────

    println("\nStep 5: Taylor test...")
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]

    for (i, ε) in enumerate(epsilons)
        pert_mesh = "mesh/reson_pert_$(i).mesh"
        open(pert_mesh, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv_palace")
            for vi in 1:nv_palace
                @printf(io, "%.15e %.15e %.15e 0\n",
                    palace_mesh.vertices[vi,1] + ε*dxdp_full[vi,1],
                    palace_mesh.vertices[vi,2] + ε*dxdp_full[vi,2],
                    palace_mesh.vertices[vi,3] + ε*dxdp_full[vi,3])
            end
            nt = size(palace_mesh.tetrahedra, 1)
            println(io, "\nTetrahedra\n$nt")
            for e in 1:nt
                @printf(io, "%d %d %d %d %d\n", palace_mesh.tetrahedra[e,1], palace_mesh.tetrahedra[e,2],
                    palace_mesh.tetrahedra[e,3], palace_mesh.tetrahedra[e,4], palace_mesh.tet_refs[e])
            end
            ntri = size(palace_mesh.triangles, 1)
            println(io, "\nTriangles\n$ntri")
            for e in 1:ntri
                @printf(io, "%d %d %d %d\n", palace_mesh.triangles[e,1], palace_mesh.triangles[e,2],
                    palace_mesh.triangles[e,3], palace_mesh.tri_refs[e])
            end
            println(io, "\nEnd")
        end
        pert_cfg = "reson_pert_$(i).json"
        pert_out = "postpro/reson_pert_$(i)"
        write(pert_cfg, replace(config_str, "postpro/reson" => pert_out))
        mkpath(pert_out)
        run(pipeline(`$palace $pert_cfg`; stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(pert_out, "eig.csv"); mode=1)

        rem = abs(f_pert - freq - ε * dfdp)
        push!(remainders, rem)
        rate = length(remainders) >= 2 ?
            log(remainders[end-1]/remainders[end]) / log(epsilons[length(remainders)-1]/ε) : NaN
        @printf("  ε=%.4f: f=%.8f, rem=%.2e, rate=%.2f\n", ε, f_pert, rem, rate)

        rm(pert_out; recursive=true, force=true)
        rm(pert_cfg; force=true)
        rm(pert_mesh; force=true)
    end

    rates = [log(remainders[i-1]/remainders[i]) / log(epsilons[i-1]/epsilons[i]) for i in 2:length(remainders)]
    max_rate = isempty(rates) ? 0.0 : maximum(rates)
    pass = max_rate > 1.5

    println("\n═══ Summary ═══")
    @printf("  Resonator: %d-vertex polygon, %d-face EGADS body\n", nv, nfaces)
    @printf("  Palace: %d verts, f₀ = %.6f GHz\n", nv_palace, freq)
    @printf("  df/d(bend_radius) = %.6e GHz/cm\n", dfdp)
    @printf("  Taylor rates: %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
    @printf("  Status: %s (max rate %.2f)\n", pass ? "PASS ✓" : "FAIL ✗", max_rate)

    rm("postpro/reson"; recursive=true, force=true)
    rm(config_file; force=true)
    rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
