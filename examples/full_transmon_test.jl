# Full transmon sensitivity: island (cap_length, cap_gap) + resonator (bend_radius)
#
# Combines both geometry types in one cavity:
#   - Island gap: simple rectangle, parameterized by cap_length + cap_gap (ForwardDiff)
#   - Resonator gaps: CPW meander via native EGADS curves (ForwardDiff for SFV)
#   - Substrate: PEC box enclosing everything
#   - Boolean: substrate - island_gap - resonator_gap1 - ... - resonator_gapN
#
# Three design parameters, each differentiated end-to-end through Palace.

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
        if lsq<1e-30; s=0.0; d2=(x-ax)^2+(y-ay)^2
        else; s=clamp(((x-ax)*dab_x+(y-ay)*dab_y)/lsq,0.0,1.0)
            px,py=ax+s*dab_x,ay+s*dab_y; d2=(x-px)^2+(y-py)^2; end
        if d2<min_d2; min_d2=d2; bs=s; bi=i; end
    end
    j = mod1(bi+1,nv)
    vx = (1-bs)*sfv.velocity[bi,1]+bs*sfv.velocity[j,1]
    vy = (1-bs)*sfv.velocity[bi,2]+bs*sfv.velocity[j,2]
    return [vx, vy, 0.0]
end

# ─── Island gap (ForwardDiff-compatible, μm coordinates) ─────────────────────

const CAP_WIDTH_UM = 24.0
const JUNCTION_GAP_UM = 12.0

function island_gap_coords_μm(params)
    # params = [cap_length_μm, cap_gap_μm]
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
        push!(coords, Unitful.ustrip(μm, getx(pt)))
        push!(coords, Unitful.ustrip(μm, gety(pt)))
    end
    return coords
end

# ─── Resonator (ForwardDiff for SFV, EGADSKernel for geometry) ───────────────

function resonator_gap_coords_μm(bend_radius_μm)
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
        push!(coords, Unitful.ustrip(μm, getx(pt)))
        push!(coords, Unitful.ustrip(μm, gety(pt)))
    end
    return coords
end

function build_resonator_solids(bend_radius_μm)
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
    sm = SolidModel("reson", k)
    render!(sm, cs; map_meta=m -> string(DeviceLayout.layer(m)))

    DLE.start_recording!(k.geometry_tape)
    for (_, pg) in DeviceLayout.SolidModels.dimgroupdict(sm, 2)
        DeviceLayout.SolidModels.extrude_z!(pg, 80.0μm)
    end
    DLE.stop_recording!(k.geometry_tape)

    solids = EGADS.EGObject[]
    for b in DLE._bodies(k)
        info = EGADS.get_info(b)
        if info.oclass == EGADS.OClass.BODY && info.mtype == EGADS.MType.SOLIDBODY
            push!(solids, b)
        end
    end
    return (k=k, solids=solids)
end

# ─── EGADS helpers ────────────────────────────────────────────────────────────

function make_extruded_polygon_μm(ctx, verts_μm, z0, height)
    nv = size(verts_μm, 1)
    nodes = [EGADS.make_node(ctx, verts_μm[i,1], verts_μm[i,2], z0) for i in 1:nv]
    edges = [EGADS.make_edge_from_nodes(ctx, nodes[i], nodes[mod1(i+1,nv)]) for i in 1:nv]
    loop = EGADS.make_wire_loop(ctx, edges)
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    return EGADS.extrude(src, height, [0.0, 0.0, 1.0])
end

# ─── Chain contraction for a given parameter ─────────────────────────────────

function compute_chain(result, result_faces, parent_bodies, palace_mesh, dfdx, is_boundary,
                        sfv, gap_height, source_origin; gap_body_indices)
    dir = [0.0, 0.0, 1.0]
    prov = match_provenance(result, parent_bodies)
    nv = size(palace_mesh.vertices, 1)
    dxdp_bdy = zeros(Float64, nv, 3)
    n_mapped = 0; n_gap = 0

    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            inv_res = EGADS.inv_evaluate(result_faces[fi], pt)
            if LinearAlgebra.norm(inv_res.results[1:3] - pt) < 1e-4
                if haskey(prov, fi) && prov[fi].body_idx in gap_body_indices
                    α = clamp(dot(pt - source_origin, dir) / gap_height, 0.0, 1.0)
                    pos_src = pt .- (α * gap_height) .* dir
                    vel = interpolate_source_velocity(sfv, pos_src[1], pos_src[2])
                    dxdp_bdy[vi, :] .= vel
                    if norm(vel) > 1e-10; n_gap += 1; end
                end
                n_mapped += 1; break
            end
        end
    end

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    return (dfdp=dfdp, n_mapped=n_mapped, n_gap=n_gap, dxdp_full=dxdp_full)
end

# ─── Taylor test ─────────────────────────────────────────────────────────────

function run_taylor(palace, palace_mesh, freq, dfdp, dxdp_full, nfaces, label)
    nv = size(palace_mesh.vertices, 1)
    epsilons = [0.1, 0.05, 0.025]
    remainders = Float64[]
    for (i, ε) in enumerate(epsilons)
        pm = "mesh/ft_$(label)_$(i).mesh"
        open(pm, "w") do io
            println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
            println(io, "Vertices\n$nv")
            for vi in 1:nv
                @printf(io, "%.15e %.15e %.15e 0\n",
                    palace_mesh.vertices[vi,1]+ε*dxdp_full[vi,1],
                    palace_mesh.vertices[vi,2]+ε*dxdp_full[vi,2],
                    palace_mesh.vertices[vi,3]+ε*dxdp_full[vi,3])
            end
            nt = size(palace_mesh.tetrahedra, 1)
            println(io, "\nTetrahedra\n$nt")
            for e in 1:nt
                @printf(io, "%d %d %d %d %d\n", palace_mesh.tetrahedra[e,1:4]..., palace_mesh.tet_refs[e])
            end
            ntri = size(palace_mesh.triangles, 1)
            println(io, "\nTriangles\n$ntri")
            for e in 1:ntri
                @printf(io, "%d %d %d %d\n", palace_mesh.triangles[e,1:3]..., palace_mesh.tri_refs[e])
            end
            println(io, "\nEnd")
        end
        pc = "ft_$(label)_$(i).json"; po = "postpro/ft_$(label)_$(i)"
        write(pc, """{
  "Problem":{"Type":"Eigenmode","Verbose":2,"Output":"$po"},
  "Model":{"Mesh":"$pm","L0":1.0e-4},
  "Domains":{"Materials":[{"Attributes":[1],"Permeability":1.0,"Permittivity":1.0}]},
  "Boundaries":{"PEC":{"Attributes":$(collect(1:nfaces))}},
  "Solver":{"Order":1,"Device":"CPU",
    "Eigenmode":{"N":2,"Tol":1.0e-8,"Target":10.0},
    "Linear":{"Type":"Default","Tol":1.0e-8,"MaxIts":200}}}""")
        mkpath(po)
        run(pipeline(`$palace $pc`; stdout=devnull, stderr=devnull))
        fp = read_eigenfrequency(joinpath(po, "eig.csv"); mode=1)
        rem = abs(fp - freq - ε * dfdp)
        push!(remainders, rem)
        rate = length(remainders) >= 2 ?
            log(remainders[end-1]/remainders[end]) / log(epsilons[length(remainders)-1]/ε) : NaN
        @printf("    ε=%.3f: f=%.6f, rem=%.2e, rate=%.2f\n", ε, fp, rem, rate)
        rm(po; recursive=true, force=true); rm(pc; force=true); rm(pm; force=true)
    end
    rates = [log(remainders[i-1]/remainders[i])/log(epsilons[i-1]/epsilons[i]) for i in 2:length(remainders)]
    return rates
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace)

    cap_length_μm = 620.0
    cap_gap_μm = 30.0
    bend_radius_μm = 50.0
    gap_height = 80.0  # μm extrusion height

    println("═══ Full Transmon Sensitivity: 3 Parameters ═══")
    @printf("  cap_length=%.0f μm, cap_gap=%.0f μm, bend_radius=%.0f μm\n",
            cap_length_μm, cap_gap_μm, bend_radius_μm)

    # ─── ForwardDiff: island and resonator SFV ───────────────────────────────

    println("\nStep 1: ForwardDiff Jacobians...")

    # Island: 2-param Jacobian
    island_coords = island_gap_coords_μm([cap_length_μm, cap_gap_μm])
    J_island = ForwardDiff.jacobian(island_gap_coords_μm, [cap_length_μm, cap_gap_μm])
    nv_island = length(island_coords) ÷ 2
    verts_island = coords_to_vertices(island_coords)

    sfv_cl = SourceFaceVelocity(verts_island, jacobian_column_to_velocity(J_island, 1))
    sfv_cg = SourceFaceVelocity(verts_island, jacobian_column_to_velocity(J_island, 2))
    println("  Island: $nv_island vertices, 2-param Jacobian")

    # Resonator: 1-param derivative
    reson_coords = resonator_gap_coords_μm(bend_radius_μm)
    d_reson = ForwardDiff.derivative(resonator_gap_coords_μm, bend_radius_μm)
    nv_reson = length(reson_coords) ÷ 2
    verts_reson = coords_to_vertices(reson_coords)
    vel_reson = reshape(d_reson, nv_reson, 2)
    sfv_br = SourceFaceVelocity(verts_reson, vel_reson)
    println("  Resonator: $nv_reson vertices, 1-param derivative")

    # ─── Build EGADS geometry ────────────────────────────────────────────────

    println("\nStep 2: Building geometry...")

    # Resonator solids (native curves via EGADSKernel)
    reson = build_resonator_solids(bend_radius_μm)
    reson_solids = reson.solids
    ctx = DLE._ctx(reson.k)
    println("  Resonator: $(length(reson_solids)) native-curve solids")

    # Island solid (extruded rectangle)
    # Position island above the resonator (offset in Y)
    # Center island at the resonator's coupling point
    island_offset_y = 50.0  # μm above resonator start
    verts_island_placed = copy(verts_island)
    verts_island_placed[:, 2] .+= island_offset_y

    island_solid = make_extruded_polygon_μm(ctx, verts_island_placed, 0.0, gap_height)
    println("  Island: 1 solid, $(nv_island) vertices")

    # Update island SFV positions to match placed coordinates
    sfv_cl = SourceFaceVelocity(verts_island_placed, jacobian_column_to_velocity(J_island, 1))
    sfv_cg = SourceFaceVelocity(verts_island_placed, jacobian_column_to_velocity(J_island, 2))

    # Substrate bounding box (enclose both island and resonator)
    all_bodies = vcat([island_solid], reson_solids)
    bb_min = [Inf, Inf, Inf]; bb_max = [-Inf, -Inf, -Inf]
    for s in all_bodies
        bb = EGADS.get_bounding_box(s)
        bb_min .= min.(bb_min, [bb.xmin, bb.ymin, bb.zmin])
        bb_max .= max.(bb_max, [bb.xmax, bb.ymax, bb.zmax])
    end
    margin = 50.0
    sub = EGADS.make_box(ctx, bb_min[1]-margin, bb_min[2]-margin, 0.0,
                         bb_max[1]-bb_min[1]+2margin, bb_max[2]-bb_min[2]+2margin,
                         gap_height + 20.0)  # substrate taller than gaps

    # Sequential booleans: substrate - island - reson1 - reson2 - ...
    result = sub
    for body in all_bodies
        rv = EGADS.solid_boolean(result, body, EGADS.BoolOp.SUBTRACTION)
        isempty(rv) && continue
        result = rv[1]
    end

    result_faces = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    nfaces = length(result_faces)
    println("  Result: $nfaces faces")

    # ─── Mesh + Palace ───────────────────────────────────────────────────────

    println("\nStep 3: TetGen + Palace...")
    tess = EGADS.tessellate(result; max_edge=30.0, sag=2.0, angle=12.0)
    npts_tess = EGADS.tess_status(tess).npts

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
                push!(tris, Int32(EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,k])))
            end
            push!(tri_markers, Int32(fi))
        end
    end

    input = RawTetGenIO{Float64}()
    input.pointlist = pts
    facetlist!(input, reshape(tris, 3, length(tri_markers)))
    input.facetmarkerlist = tri_markers
    tetout = tetrahedralize(input, "pq1.2a50000")
    nvol = numberofpoints(tetout)
    ntet = numberoftetrahedra(tetout)
    @printf("  Mesh: %d verts, %d tets\n", nvol, ntet)

    mesh_file = "mesh/ft_baseline.mesh"
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
  "Problem":{"Type":"Eigenmode","Verbose":2,"Output":"postpro/ft_baseline"},
  "Model":{"Mesh":"$mesh_file","L0":1.0e-4},
  "Domains":{"Materials":[{"Attributes":[1],"Permeability":1.0,"Permittivity":1.0}]},
  "Boundaries":{"PEC":{"Attributes":$(collect(1:nfaces))}},
  "Solver":{"Order":1,"Device":"CPU",
    "Eigenmode":{"N":2,"Tol":1.0e-8,"Target":10.0,"Sensitivity":true,"MeshSensitivity":true},
    "Linear":{"Type":"Default","Tol":1.0e-8,"MaxIts":200}}}"""
    write("ft_baseline.json", config_str)
    mkpath("postpro/ft_baseline")
    flush(stdout)
    t0 = time()
    run(pipeline(`$palace ft_baseline.json`; stdout=devnull, stderr=devnull))
    @printf("  Palace: %.1fs\n", time()-t0)

    freq = read_eigenfrequency("postpro/ft_baseline/eig.csv"; mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    palace_mesh = read_mesh("postpro/ft_baseline/mesh_for_sensitivity.mesh")
    sol = read_solb("postpro/ft_baseline/mesh_sensitivity.solb")
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv_palace = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv_palace, count(is_boundary))

    # ─── Chain contraction for each parameter ────────────────────────────────

    # Parent bodies: [substrate, island, reson1, reson2, ...]
    parent_bodies = vcat([sub, island_solid], reson_solids)

    # Find extrusion bases
    island_faces = EGADS.get_body_topos(island_solid, EGADS.OClass.FACE)
    island_base = find_extrusion_base(island_faces, [0.0, 0.0, 1.0])

    all_reson_faces = vcat([EGADS.get_body_topos(s, EGADS.OClass.FACE) for s in reson_solids]...)
    reson_base = find_extrusion_base(all_reson_faces, [0.0, 0.0, 1.0])

    # Island body index = 2, resonator body indices = 3:(2+length(reson_solids))
    island_idx = Set([2])
    reson_indices = Set(3:(2 + length(reson_solids)))

    # --- cap_length ---
    println("\n═══ d/d(cap_length) ═══")
    ch_cl = compute_chain(result, result_faces, parent_bodies, palace_mesh, dfdx, is_boundary,
                           sfv_cl, gap_height, island_base; gap_body_indices=island_idx)
    @printf("  mapped=%d, gap_vel=%d, df/dp=%.6e\n", ch_cl.n_mapped, ch_cl.n_gap, ch_cl.dfdp)
    rates_cl = run_taylor(palace, palace_mesh, freq, ch_cl.dfdp, ch_cl.dxdp_full, nfaces, "cl")

    # --- cap_gap ---
    println("\n═══ d/d(cap_gap) ═══")
    ch_cg = compute_chain(result, result_faces, parent_bodies, palace_mesh, dfdx, is_boundary,
                           sfv_cg, gap_height, island_base; gap_body_indices=island_idx)
    @printf("  mapped=%d, gap_vel=%d, df/dp=%.6e\n", ch_cg.n_mapped, ch_cg.n_gap, ch_cg.dfdp)
    rates_cg = run_taylor(palace, palace_mesh, freq, ch_cg.dfdp, ch_cg.dxdp_full, nfaces, "cg")

    # --- bend_radius ---
    println("\n═══ d/d(bend_radius) ═══")
    ch_br = compute_chain(result, result_faces, parent_bodies, palace_mesh, dfdx, is_boundary,
                           sfv_br, gap_height, reson_base; gap_body_indices=reson_indices)
    @printf("  mapped=%d, gap_vel=%d, df/dp=%.6e\n", ch_br.n_mapped, ch_br.n_gap, ch_br.dfdp)
    rates_br = run_taylor(palace, palace_mesh, freq, ch_br.dfdp, ch_br.dxdp_full, nfaces, "br")

    # ─── Summary ─────────────────────────────────────────────────────────────

    println("\n═══ Full Transmon — 3 Parameter Summary ═══")
    @printf("  f₀ = %.6f GHz\n", freq)
    for (name, ch, rates) in [("cap_length", ch_cl, rates_cl),
                               ("cap_gap", ch_cg, rates_cg),
                               ("bend_radius", ch_br, rates_br)]
        max_r = isempty(rates) ? 0.0 : maximum(rates)
        pass = max_r > 1.5
        @printf("  %-12s: df/dp=%+.4e  rates=[%s]  %s\n",
                name, ch.dfdp, join([@sprintf("%.2f", r) for r in rates], ","), pass ? "PASS ✓" : "FAIL ✗")
    end

    all_pass = all(r -> !isempty(r) && maximum(r) > 1.5, [rates_cl, rates_cg, rates_br])
    println("  Overall: $(all_pass ? "PASS ✓" : "FAIL ✗")")

    rm("postpro/ft_baseline"; recursive=true, force=true)
    rm("ft_baseline.json"; force=true); rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
