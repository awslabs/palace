# Resonator sensitivity via native EGADS curves + ForwardDiff SFV.
#
# Key insight: use native curves (LINE + ARC) for the 3D geometry (few EGADS faces,
# meshable by TetGen) but the FULL ForwardDiff polygon (1502 vertices) for the source
# face velocity (accurate interpolation). The polygon vertices sit on the same curves.
#
# Approach:
#   1. Build native-curve resonator via EGADSKernel render! (6 solid bodies)
#   2. Union solids, boolean subtract from substrate
#   3. ForwardDiff through DeviceLayout path → SourceFaceVelocity from full polygon
#   4. Chain contraction + Taylor test

using LinearAlgebra, Printf
using DeviceLayout, DeviceLayout.PreferredUnits
using ForwardDiff, Unitful
using EGADS
using TetGen

DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

include(joinpath(@__DIR__, "shared", "palace_io.jl"))
include(joinpath(homedir(), "AWS", "DeviceLayout.jl-ad-sensitivity", "src", "forwarddiff_bridge.jl"))

# ─── Inlined helpers ─────────────────────────────────────────────────────────

function face_plane(face)
    ri = EGADS.get_range(face)
    um = (ri.range[1] + ri.range[2]) / 2
    vm = (ri.range[3] + ri.range[4]) / 2
    res = EGADS.evaluate(face, [um, vm])
    n = cross(res.du, res.dv)
    nrm = norm(n); nrm < 1e-30 && return (normal=[0.0,0.0,0.0], offset=0.0, point=res.xyz)
    n = n / nrm
    return (normal=n, offset=dot(n, res.xyz), point=res.xyz)
end

function match_provenance(result_body, parent_bodies; tol=1e-6)
    result_faces = EGADS.get_body_topos(result_body, EGADS.OClass.FACE)
    pd = Tuple{Int,Int,EGADS.EGObject,NamedTuple}[]
    for (bi, body) in enumerate(parent_bodies)
        for (fi, pf) in enumerate(EGADS.get_body_topos(body, EGADS.OClass.FACE))
            push!(pd, (bi, fi, pf, face_plane(pf)))
        end
    end
    prov = Dict{Int, NamedTuple{(:body_idx,:face_idx), Tuple{Int,Int}}}()
    for (ri, rf) in enumerate(result_faces)
        rp = face_plane(rf)
        for (bi, fi, _, pp) in pd
            ndot = abs(dot(rp.normal, pp.normal))
            if ndot > 1-tol && (abs(rp.offset-pp.offset)<tol || abs(rp.offset+pp.offset)<tol)
                prov[ri] = (body_idx=bi, face_idx=fi); break
            end
        end
    end
    return prov
end

function find_extrusion_base(faces, dir)
    min_proj = Inf; base = [0.0,0.0,0.0]
    for f in faces
        rng = EGADS.get_range(f)
        c = EGADS.evaluate(f, [(rng.range[1]+rng.range[2])/2, (rng.range[3]+rng.range[4])/2]).xyz
        if dot(c,dir) < min_proj; min_proj = dot(c,dir); base = c; end
    end; return base
end

struct SourceFaceVelocity; vertices::Matrix{Float64}; velocity::Matrix{Float64}; end

function interpolate_source_velocity(sfv::SourceFaceVelocity, x::Float64, y::Float64)
    nv = size(sfv.vertices, 1); nv == 0 && return [0.0,0.0,0.0]
    min_d2 = Inf; bs = 0.0; bi = 1
    for i in 1:nv
        j = mod1(i+1,nv); ax,ay = sfv.vertices[i,1],sfv.vertices[i,2]
        bx,by = sfv.vertices[j,1],sfv.vertices[j,2]; dab_x,dab_y = bx-ax,by-ay
        lsq = dab_x^2+dab_y^2
        if lsq < 1e-30; s=0.0; d2=(x-ax)^2+(y-ay)^2
        else; s=clamp(((x-ax)*dab_x+(y-ay)*dab_y)/lsq,0.0,1.0)
            px,py=ax+s*dab_x,ay+s*dab_y; d2=(x-px)^2+(y-py)^2; end
        if d2 < min_d2; min_d2=d2; bs=s; bi=i; end
    end
    j = mod1(bi+1,nv)
    vx = (1-bs)*sfv.velocity[bi,1]+bs*sfv.velocity[j,1]
    vy = (1-bs)*sfv.velocity[bi,2]+bs*sfv.velocity[j,2]
    return [vx, vy, 0.0]
end

# ─── ForwardDiff through DeviceLayout (full-resolution polygon) ──────────────

function resonator_gap_coords_cm(bend_radius_μm)
    r = bend_radius_μm * 1μm; T = typeof(r)
    style = Paths.SimpleCPW(T(10.0μm), T(6.0μm))
    p = Path(Point(zero(r), zero(r)); α0=0.0°)
    straight!(p, T(200.0μm), style)
    turn!(p, -90°, r)
    straight!(p, T(300.0μm))
    turn!(p, 180°, r)
    straight!(p, T(300.0μm))

    all_polys = []
    for node in p.nodes; append!(all_polys, to_polygons(node)); end
    result = all_polys[1]
    for i in 2:length(all_polys); result = union2d(result, all_polys[i]); end
    poly = to_polygons(result)[1]

    coords = typeof(bend_radius_μm)[]
    for pt in points(poly)
        push!(coords, Unitful.ustrip(u"cm", getx(pt)))
        push!(coords, Unitful.ustrip(u"cm", gety(pt)))
    end
    return coords
end

# ─── Build native-curve geometry via EGADSKernel ─────────────────────────────

function build_native_resonator(bend_radius_μm)
    r = Float64(bend_radius_μm) * 1μm
    style = Paths.SimpleCPW(10.0μm, 6.0μm)

    p = Path(Point(0.0μm, 0.0μm); α0=0.0°)
    straight!(p, 200.0μm, style)
    turn!(p, -90°, r)
    straight!(p, 300.0μm)
    turn!(p, 180°, r)
    straight!(p, 300.0μm)

    cs = CoordinateSystem("reson", nm)
    place!(cs, p)

    k = DLE.EGADSKernel(enable_tape=true)
    sm = SolidModel("reson3d", k)
    render!(sm, cs; map_meta=m -> string(DeviceLayout.layer(m)))

    # Extrude all 2D faces
    DLE.start_recording!(k.geometry_tape)
    gap_height_μm = 10.0
    for (_, pg) in DeviceLayout.SolidModels.dimgroupdict(sm, 2)
        DeviceLayout.SolidModels.extrude_z!(pg, gap_height_μm * 1μm)
    end
    DLE.stop_recording!(k.geometry_tape)

    # Collect solid bodies
    solids = EGADS.EGObject[]
    for b in DLE._bodies(k)
        info = EGADS.get_info(b)
        if info.oclass == EGADS.OClass.BODY && info.mtype == EGADS.MType.SOLIDBODY
            push!(solids, b)
        end
    end

    return (k=k, sm=sm, gap_solids=solids, gap_height_μm=gap_height_μm)
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace)

    bend_radius_μm = 50.0
    sub_margin = 0.005

    println("═══ Resonator Native-Curve Sensitivity Test ═══")
    println("  bend_radius = $(bend_radius_μm) μm")

    # ─── Step 1: ForwardDiff → SFV (full polygon) ────────────────────────────

    println("\nStep 1: ForwardDiff through path pipeline...")
    t0 = time()
    coords_nom = resonator_gap_coords_cm(bend_radius_μm)
    d = ForwardDiff.derivative(resonator_gap_coords_cm, bend_radius_μm)
    nv_poly = length(coords_nom) ÷ 2
    @printf("  %d-vertex polygon, ForwardDiff in %.1fs\n", nv_poly, time()-t0)

    verts_poly_cm = coords_to_vertices(coords_nom)
    vel_poly_cm = reshape(d, nv_poly, 2)  # d(cm)/d(μm) — dimensionless in cm per μm
    # Convert to μm for matching with EGADSKernel geometry (STP_UNIT = μm)
    verts_poly_μm = verts_poly_cm .* 1e4   # cm → μm
    vel_poly_μm = vel_poly_cm .* 1e4        # d(cm)/d(μm) * (μm/cm) = d(μm)/d(μm) = dimensionless...
    # Actually: ForwardDiff gives d(coord_cm)/d(param_μm). The mesh is in μm.
    # For the chain: df/dp = dfdx_μm · dxdp_μm where dxdp_μm = d(x_μm)/d(p_μm)
    # vel_poly_cm = d(x_cm)/d(p_μm), so vel_poly_μm = vel_poly_cm * 1e4
    sfv = SourceFaceVelocity(verts_poly_μm, vel_poly_μm)

    # ─── Step 2: Native-curve EGADS geometry ─────────────────────────────────

    println("\nStep 2: Building native-curve geometry via EGADSKernel...")
    reson = build_native_resonator(bend_radius_μm)
    gap_solids = reson.gap_solids
    gap_height_μm = reson.gap_height_μm

    total_gap_faces = sum(length(EGADS.get_body_topos(s, EGADS.OClass.FACE)) for s in gap_solids)
    println("  Gap solids: $(length(gap_solids)), total faces: $total_gap_faces")

    # Bounding box for substrate — in μm (STP_UNIT) to match EGADSKernel
    verts_μm = verts_poly_cm .* 1e4  # cm → μm
    xmin, xmax = extrema(verts_μm[:, 1])
    ymin, ymax = extrema(verts_μm[:, 2])
    margin_μm = 50.0
    sub_x0 = xmin - margin_μm
    sub_y0 = ymin - margin_μm
    sub_w = (xmax - xmin) + 2margin_μm
    sub_h = (ymax - ymin) + 2margin_μm
    sub_z = 100.0  # μm

    ctx = DLE._ctx(reson.k)
    substrate = EGADS.make_box(ctx, sub_x0, sub_y0, 0.0, sub_w, sub_h, sub_z)

    # Sequential subtractions: substrate - gap1 - gap2 - ...
    result = substrate
    for (i, gap) in enumerate(gap_solids)
        result_vec = EGADS.solid_boolean(result, gap, EGADS.BoolOp.SUBTRACTION)
        if isempty(result_vec)
            @warn "Boolean subtraction $i produced empty result, skipping"
            continue
        end
        result = result_vec[1]
    end

    result_faces = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    nfaces = length(result_faces)
    println("  Result: $nfaces faces (substrate - $(length(gap_solids)) gaps)")

    # ─── Step 3: Mesh + Palace ───────────────────────────────────────────────

    println("\nStep 3: TetGen + Palace...")
    tess = EGADS.tessellate(result; max_edge=50.0, sag=5.0, angle=15.0)  # μm scale
    npts_tess = EGADS.tess_status(tess).npts
    println("  Tessellation: $npts_tess surface vertices")

    pts = zeros(Float64, 3, npts_tess)
    for fi in eachindex(result_faces)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end
    tris = Int32[]; tri_markers = Int32[]
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
    tetout = tetrahedralize(input, "pq1.2a50000")  # max tet vol ~50000 μm³
    nvol = numberofpoints(tetout)
    ntet = numberoftetrahedra(tetout)
    @printf("  Volume: %d verts, %d tets\n", nvol, ntet)

    mesh_file = "mesh/reson_native.mesh"
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
  "Problem": { "Type": "Eigenmode", "Verbose": 2, "Output": "postpro/reson_native" },
  "Model": { "Mesh": "$mesh_file", "L0": 1.0e-4 },
  "Domains": { "Materials": [{ "Attributes": [1], "Permeability": 1.0, "Permittivity": 1.0 }] },
  "Boundaries": { "PEC": { "Attributes": $(collect(1:nfaces)) } },
  "Solver": { "Order": 1, "Device": "CPU",
    "Eigenmode": { "N": 2, "Tol": 1.0e-8, "Target": 30.0, "Sensitivity": true, "MeshSensitivity": true },
    "Linear": { "Type": "Default", "Tol": 1.0e-8, "MaxIts": 200 } }
}"""
    write("reson_native.json", config_str)
    mkpath("postpro/reson_native")

    flush(stdout); t0 = time()
    run(pipeline(`$palace reson_native.json`; stdout=devnull, stderr=devnull))
    @printf("  Palace: %.1fs\n", time()-t0)

    freq = read_eigenfrequency("postpro/reson_native/eig.csv"; mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    palace_mesh = read_mesh("postpro/reson_native/mesh_for_sensitivity.mesh")
    sol = read_solb("postpro/reson_native/mesh_sensitivity.solb")
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv_palace = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv_palace, count(is_boundary))

    # ─── Step 4: Chain contraction ───────────────────────────────────────────

    println("\nStep 4: Chain contraction...")
    parent_bodies = vcat([substrate], gap_solids)
    prov = match_provenance(result, parent_bodies)
    @printf("  Provenance: %d / %d faces\n", length(prov), nfaces)

    dir = [0.0, 0.0, 1.0]
    all_gap_faces = vcat([EGADS.get_body_topos(s, EGADS.OClass.FACE) for s in gap_solids]...)
    source_origin = find_extrusion_base(all_gap_faces, dir)

    dxdp_bdy = zeros(Float64, nv_palace, 3)
    n_mapped = 0; n_gap = 0
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            inv_res = EGADS.inv_evaluate(result_faces[fi], pt)
            if LinearAlgebra.norm(inv_res.results[1:3] - pt) < 1e-6
                if haskey(prov, fi) && prov[fi].body_idx > 1  # any gap solid
                    α = clamp(dot(pt - source_origin, dir) / gap_height_μm, 0.0, 1.0)
                    pos_src = pt .- (α * gap_height_μm) .* dir
                    vel = interpolate_source_velocity(sfv, pos_src[1], pos_src[2])
                    dxdp_bdy[vi, :] .= vel
                    if norm(vel) > 1e-10; n_gap += 1; end
                end
                n_mapped += 1; break
            end
        end
    end
    @printf("  Mapped: %d / %d boundary, %d with nonzero velocity\n",
            n_mapped, count(is_boundary), n_gap)

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    @printf("  df/d(bend_radius) = %.6e GHz/cm\n", dfdp)

    # ─── Step 5: Taylor test ─────────────────────────────────────────────────

    println("\nStep 5: Taylor test...")
    epsilons = [0.1, 0.05, 0.025, 0.0125]
    remainders = Float64[]
    for (i, ε) in enumerate(epsilons)
        pert_mesh = "mesh/rn_pert_$(i).mesh"
        open(pert_mesh, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv_palace")
            for vi in 1:nv_palace
                @printf(io, "%.15e %.15e %.15e 0\n",
                    palace_mesh.vertices[vi,1]+ε*dxdp_full[vi,1],
                    palace_mesh.vertices[vi,2]+ε*dxdp_full[vi,2],
                    palace_mesh.vertices[vi,3]+ε*dxdp_full[vi,3])
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
        pert_cfg = "rn_pert_$(i).json"
        pert_out = "postpro/rn_pert_$(i)"
        write(pert_cfg, replace(config_str, "reson_native" => "rn_pert_$(i)"))
        mkpath(pert_out)
        run(pipeline(`$palace $pert_cfg`; stdout=devnull, stderr=devnull))
        f_pert = read_eigenfrequency(joinpath(pert_out, "eig.csv"); mode=1)
        rem = abs(f_pert - freq - ε * dfdp)
        push!(remainders, rem)
        rate = length(remainders) >= 2 ?
            log(remainders[end-1]/remainders[end]) / log(epsilons[length(remainders)-1]/ε) : NaN
        @printf("  ε=%.4f: f=%.8f, rem=%.2e, rate=%.2f\n", ε, f_pert, rem, rate)
        rm(pert_out; recursive=true, force=true); rm(pert_cfg; force=true); rm(pert_mesh; force=true)
    end

    rates = [log(remainders[i-1]/remainders[i]) / log(epsilons[i-1]/epsilons[i]) for i in 2:length(remainders)]
    max_rate = isempty(rates) ? 0.0 : maximum(rates)
    pass = max_rate > 1.5

    println("\n═══ Summary ═══")
    @printf("  Native-curve gap: %d EGADS faces (vs ~1510 from polygon)\n", length(gap_faces))
    @printf("  ForwardDiff SFV: %d-vertex polygon\n", nv_poly)
    @printf("  Palace: %d verts, f₀ = %.6f GHz\n", nv_palace, freq)
    @printf("  df/d(bend_radius) = %.6e GHz/cm\n", dfdp)
    @printf("  Taylor rates: %s\n", join([@sprintf("%.2f", r) for r in rates], ", "))
    @printf("  Status: %s (max rate %.2f)\n", pass ? "PASS ✓" : "FAIL ✗", max_rate)

    rm("postpro/reson_native"; recursive=true, force=true)
    rm("reson_native.json"; force=true); rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
