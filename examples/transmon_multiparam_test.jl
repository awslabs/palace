# Multi-parameter transmon sensitivity: cap_length AND cap_gap
#
# Geometry (cm):
#   Substrate: 0.1 × 0.1 × 0.05 cm PEC box
#   Outer gap: (cap_width + 2*cap_gap) × (junction_gap + cap_length + cap_gap) rect, Z-extruded
#   Inner island: cap_width × cap_length rect, Z-extruded (smaller, inside outer gap)
#   Boolean: substrate - outer_gap + inner_island (= metal everywhere except gap)
#   Actually simpler: substrate - outer_gap, then later add island back.
#   For PEC cavity: just substrate - outer_gap is a re-entrant cavity.
#
# Design params:
#   cap_length (μm) → changes both outer gap and inner island Y-extent
#   cap_gap (μm) → changes only the outer gap X and Y margins
#
# For this test we use a single gap rectangle (no inner island) parameterized
# by BOTH cap_length and cap_gap. This tests multi-parameter sensitivity
# where both parameters affect the same polygon.

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
        if dist_sq < min_dist_sq; min_dist_sq = dist_sq; best_s = s; best_i = i; end
    end
    j = mod1(best_i + 1, nv)
    vx = (1 - best_s) * sfv.velocity[best_i, 1] + best_s * sfv.velocity[j, 1]
    vy = (1 - best_s) * sfv.velocity[best_i, 2] + best_s * sfv.velocity[j, 2]
    return [vx, vy, 0.0]
end

# ─── DeviceLayout gap rectangle parameterized by cap_length AND cap_gap ──────

const CAP_WIDTH_UM = 24.0
const JUNCTION_GAP_UM = 12.0

function gap_coords_cm(params)
    cap_length = params[1] * 1μm
    cap_gap = params[2] * 1μm
    cap_width = CAP_WIDTH_UM * 1μm
    junction_gap = JUNCTION_GAP_UM * 1μm

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

# ─── EGADS helpers ────────────────────────────────────────────────────────────

function make_extruded_polygon(ctx, verts_cm, z0, height)
    nv = size(verts_cm, 1)
    nodes = [EGADS.make_node(ctx, verts_cm[i, 1], verts_cm[i, 2], z0) for i in 1:nv]
    edges = [EGADS.make_edge_from_nodes(ctx, nodes[i], nodes[mod1(i+1, nv)]) for i in 1:nv]
    loop = EGADS.make_wire_loop(ctx, edges)
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    return EGADS.extrude(src, height, [0.0, 0.0, 1.0])
end

function mesh_body(body, mesh_file; max_edge=0.01, tetgen_opts="pq1.2a1e-7")
    tess = EGADS.tessellate(body; max_edge=max_edge, sag=0.0005, angle=15.0)
    npts = EGADS.tess_status(tess).npts
    faces = EGADS.get_body_topos(body, EGADS.OClass.FACE)
    pts = zeros(Float64, 3, npts)
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end
    tris = Int32[]; tri_markers = Int32[]
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, body, fi)
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
    output = tetrahedralize(input, tetgen_opts)
    nvol = numberofpoints(output)
    ntet = numberoftetrahedra(output)
    mkpath(dirname(mesh_file))
    open(mesh_file, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        p = output.pointlist
        println(io, "Vertices\n$nvol")
        for i in 1:nvol; @printf(io, "%.15e %.15e %.15e 0\n", p[1,i], p[2,i], p[3,i]); end
        te = output.tetrahedronlist
        println(io, "\nTetrahedra\n$ntet")
        for i in 1:ntet; @printf(io, "%d %d %d %d 1\n", te[1,i], te[2,i], te[3,i], te[4,i]); end
        tf = output.trifacelist; tm = output.trifacemarkerlist; ntri = size(tf, 2)
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri; @printf(io, "%d %d %d %d\n", tf[1,i], tf[2,i], tf[3,i], tm[i]); end
        println(io, "\nEnd")
    end
    return (nvol=nvol, ntet=ntet, nfaces=length(faces))
end

# ─── Chain for a given parameter ─────────────────────────────────────────────

function compute_chain_for_param(geo, palace_mesh, dfdx, is_boundary, sfv, gap_height, source_origin)
    parent_bodies = [geo.substrate, geo.gap]
    prov = match_provenance(geo.result, parent_bodies)
    dir = [0.0, 0.0, 1.0]

    result_faces = EGADS.get_body_topos(geo.result, EGADS.OClass.FACE)
    nv = size(palace_mesh.vertices, 1)
    dxdp_bdy = zeros(Float64, nv, 3)
    n_mapped = 0

    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            inv_res = EGADS.inv_evaluate(result_faces[fi], pt)
            if LinearAlgebra.norm(inv_res.results[1:3] - pt) < 1e-6
                if haskey(prov, fi) && prov[fi].body_idx == 2
                    α = clamp(dot(pt - source_origin, dir) / gap_height, 0.0, 1.0)
                    pos_source = pt .- (α * gap_height) .* dir
                    vel = interpolate_source_velocity(sfv, pos_source[1], pos_source[2])
                    dxdp_bdy[vi, :] .= vel
                end
                n_mapped += 1
                break
            end
        end
    end

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    return (dfdp=dfdp, n_mapped=n_mapped, dxdp_full=dxdp_full)
end

# ─── Taylor test helper ──────────────────────────────────────────────────────

function run_taylor_test(palace, palace_mesh, freq, dfdp, dxdp_full, nfaces, label)
    nv = size(palace_mesh.vertices, 1)
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]

    for (i, ε) in enumerate(epsilons)
        pert_mesh = "mesh/mp_pert_$(label)_$(i).mesh"
        open(pert_mesh, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv")
            for vi in 1:nv
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

        pert_cfg = "mp_pert_$(label)_$(i).json"
        pert_out = "postpro/mp_pert_$(label)_$(i)"
        config_str = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "$pert_out" },
  "Model": { "Mesh": "$pert_mesh", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": { "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 10.0 },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 } }
}"""
        write(pert_cfg, config_str)
        mkpath(pert_out)
        run(pipeline(`$palace $pert_cfg`; stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(pert_out, "eig.csv"); mode=1)

        rem = abs(f_pert - freq - ε * dfdp)
        push!(remainders, rem)
        rate = length(remainders) >= 2 ?
            log(remainders[end-1]/remainders[end]) / log(epsilons[length(remainders)-1]/ε) : NaN
        @printf("    ε=%.4f: f=%.8f, rem=%.2e, rate=%.2f\n", ε, f_pert, rem, rate)

        rm(pert_out; recursive=true, force=true)
        rm(pert_cfg; force=true)
        rm(pert_mesh; force=true)
    end

    rates = [log(remainders[i-1]/remainders[i]) / log(epsilons[i-1]/epsilons[i]) for i in 2:length(remainders)]
    return rates
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace)

    cap_length_μm = 620.0
    cap_gap_μm = 30.0
    params_nom = [cap_length_μm, cap_gap_μm]

    sub_w, sub_h, sub_z = 0.1, 0.1, 0.05
    gap_z0 = 0.001
    gap_height = sub_z - gap_z0 - 0.001

    println("═══ Multi-Parameter Transmon Sensitivity Test ═══")
    println("  params: cap_length=$(cap_length_μm)μm, cap_gap=$(cap_gap_μm)μm")

    # ─── ForwardDiff ─────────────────────────────────────────────────────────

    println("\nStep 1: ForwardDiff Jacobian...")
    coords_nom = gap_coords_cm(params_nom)
    J = ForwardDiff.jacobian(gap_coords_cm, params_nom)
    nv_poly = length(coords_nom) ÷ 2
    verts_nom = coords_to_vertices(coords_nom)

    # Center in substrate
    cx = (maximum(verts_nom[:, 1]) + minimum(verts_nom[:, 1])) / 2
    cy = (maximum(verts_nom[:, 2]) + minimum(verts_nom[:, 2])) / 2
    verts_nom[:, 1] .+= (sub_w / 2 - cx)
    verts_nom[:, 2] .+= (sub_h / 2 - cy)

    # Build SFV for each parameter
    vel_cl = jacobian_column_to_velocity(J, 1)  # d/d(cap_length)
    vel_cg = jacobian_column_to_velocity(J, 2)  # d/d(cap_gap)
    sfv_cl = SourceFaceVelocity(verts_nom, vel_cl)
    sfv_cg = SourceFaceVelocity(verts_nom, vel_cg)

    for (name, vel) in [("cap_length", vel_cl), ("cap_gap", vel_cg)]
        n_mov = count(i -> maximum(abs.(vel[i, :])) > 1e-10, 1:nv_poly)
        println("  d/d($name): $n_mov / $nv_poly vertices moving")
    end

    # ─── EGADS geometry ──────────────────────────────────────────────────────

    println("\nStep 2: EGADS geometry...")
    ctx = EGADS.EGADSContext()
    EGADS.set_out_level!(ctx, 0)

    substrate = EGADS.make_box(ctx, 0.0, 0.0, 0.0, sub_w, sub_h, sub_z)
    gap = make_extruded_polygon(ctx, verts_nom, gap_z0, gap_height)
    result_vec = EGADS.solid_boolean(substrate, gap, EGADS.BoolOp.SUBTRACTION)
    result = result_vec[1]
    nfaces = length(EGADS.get_body_topos(result, EGADS.OClass.FACE))
    println("  Result: $nfaces faces")

    geo = (substrate=substrate, gap=gap, result=result)

    # ─── Mesh + Palace ───────────────────────────────────────────────────────

    println("\nStep 3: Mesh + Palace...")
    mesh_file = "mesh/mp_baseline.mesh"
    minfo = mesh_body(result, mesh_file)
    @printf("  Mesh: %d verts, %d tets\n", minfo.nvol, minfo.ntet)

    config_file = "mp_baseline.json"
    outdir = "postpro/mp_baseline"
    config_str = """{
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "$outdir" },
  "Model": { "Mesh": "$mesh_file", "L0": 1.0e-2 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": { "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 10.0, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 } }
}"""
    write(config_file, config_str)
    mkpath(outdir)
    flush(stdout)
    t0 = time()
    run(pipeline(`$palace $config_file`; stdout=devnull, stderr=devnull))
    @printf("  Palace: %.1fs\n", time() - t0)

    freq = read_eigenfrequency(joinpath(outdir, "eig.csv"); mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    palace_mesh = read_mesh(joinpath(outdir, "mesh_for_sensitivity.mesh"))
    sol = read_solb(joinpath(outdir, "mesh_sensitivity.solb"))
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    @printf("  Palace mesh: %d verts, %d boundary\n", size(palace_mesh.vertices, 1), count(is_boundary))

    # Gap extrusion base
    gap_faces = EGADS.get_body_topos(gap, EGADS.OClass.FACE)
    source_origin = find_extrusion_base(gap_faces, [0.0, 0.0, 1.0])

    # ─── Parameter 1: cap_length ─────────────────────────────────────────────

    println("\n═══ Parameter 1: d/d(cap_length) ═══")
    chain_cl = compute_chain_for_param(geo, palace_mesh, dfdx, is_boundary, sfv_cl, gap_height, source_origin)
    @printf("  Mapped: %d / %d\n", chain_cl.n_mapped, count(is_boundary))
    @printf("  df/d(cap_length) = %.6e GHz/cm\n", chain_cl.dfdp)

    rates_cl = run_taylor_test(palace, palace_mesh, freq, chain_cl.dfdp, chain_cl.dxdp_full, nfaces, "cl")

    # ─── Parameter 2: cap_gap ────────────────────────────────────────────────

    println("\n═══ Parameter 2: d/d(cap_gap) ═══")
    chain_cg = compute_chain_for_param(geo, palace_mesh, dfdx, is_boundary, sfv_cg, gap_height, source_origin)
    @printf("  Mapped: %d / %d\n", chain_cg.n_mapped, count(is_boundary))
    @printf("  df/d(cap_gap) = %.6e GHz/cm\n", chain_cg.dfdp)

    rates_cg = run_taylor_test(palace, palace_mesh, freq, chain_cg.dfdp, chain_cg.dxdp_full, nfaces, "cg")

    EGADS.close(ctx)

    # ─── Summary ─────────────────────────────────────────────────────────────

    println("\n═══ Summary ═══")
    @printf("  f₀ = %.6f GHz\n", freq)
    @printf("  df/d(cap_length) = %.6e GHz/cm  rates: %s\n",
            chain_cl.dfdp, join([@sprintf("%.2f", r) for r in rates_cl], ", "))
    @printf("  df/d(cap_gap)    = %.6e GHz/cm  rates: %s\n",
            chain_cg.dfdp, join([@sprintf("%.2f", r) for r in rates_cg], ", "))

    pass_cl = !isempty(rates_cl) && maximum(rates_cl) > 1.5
    pass_cg = !isempty(rates_cg) && maximum(rates_cg) > 1.5
    println("  cap_length: $(pass_cl ? "PASS ✓" : "FAIL ✗")")
    println("  cap_gap:    $(pass_cg ? "PASS ✓" : "FAIL ✗")")
    println("  Overall:    $(pass_cl && pass_cg ? "PASS ✓" : "FAIL ✗")")

    rm(outdir; recursive=true, force=true)
    rm(config_file; force=true)
    rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
