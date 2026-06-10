# Full transmon sensitivity: 10-parameter Taylor tests
#
# Extends the 3-param test to verify all 10 design parameters:
#   Island (4):    cap_length, cap_gap, cap_width, junction_gap
#   Resonator (6): bend_radius, trace_width, gap_width, straight_1, straight_2, straight_3
#
# For each parameter, runs a Taylor test:
#   |f(p + ε·δp) - f(p) - ε·df/dp| = O(ε²)
#
# Pass criterion: convergence rate ≥ 1.8 for at least one epsilon pair.

using LinearAlgebra, Printf
using DeviceLayout, DeviceLayout.PreferredUnits
using ForwardDiff, Unitful
using EGADS
using TetGen

DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

include(joinpath(@__DIR__, "shared", "palace_io.jl"))
include(joinpath(homedir(), "AWS", "DeviceLayout.jl-ad-sensitivity", "src", "forwarddiff_bridge.jl"))

# ─── Helpers ──────────────────────────────────────────────────────────────────

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

# ─── Constants ────────────────────────────────────────────────────────────────

const GH = 80.0   # gap height (μm)
const IOY = 50.0   # island Y-offset (μm)

const PARAM_NAMES = ["cap_length", "cap_gap", "cap_width", "junction_gap",
                     "bend_radius", "trace_width", "gap_width",
                     "straight_1", "straight_2", "straight_3"]

# ─── Island gap (ForwardDiff-compatible, 4 params, μm) ───────────────────────

function island_gap_coords_μm(params)
    # params = [cap_length, cap_gap, cap_width, junction_gap]
    cap_length = params[1] * 1μm
    cap_gap = params[2] * 1μm
    cap_width = params[3] * 1μm
    junction_gap = params[4] * 1μm

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

# ─── Resonator gap (ForwardDiff-compatible, 6 params, μm) ────────────────────

function resonator_gap_coords_μm(params)
    # params = [bend_radius, trace_width, gap_width, straight_1, straight_2, straight_3]
    br = params[1] * 1μm
    tw = params[2] * 1μm
    gw = params[3] * 1μm
    s1 = params[4] * 1μm
    s2 = params[5] * 1μm
    s3 = params[6] * 1μm
    style = Paths.SimpleCPW(tw, gw)
    p = Path(Point(zero(br), zero(br)); α0=0.0°)
    straight!(p, s1, style)
    turn!(p, -90°, br)
    straight!(p, s2)
    turn!(p, 180°, br)
    straight!(p, s3)

    all_polys = []
    for node in p.nodes; append!(all_polys, to_polygons(node)); end
    result = all_polys[1]
    for i in 2:length(all_polys); result = union2d(result, all_polys[i]); end
    poly = to_polygons(result)[1]

    coords = eltype(params)[]
    for pt in points(poly)
        push!(coords, Unitful.ustrip(μm, getx(pt)))
        push!(coords, Unitful.ustrip(μm, gety(pt)))
    end
    return coords
end

# ─── EGADS geometry ───────────────────────────────────────────────────────────

function make_extruded_polygon_μm(ctx, verts_μm, z0, height)
    nv = size(verts_μm, 1)
    nodes = [EGADS.make_node(ctx, verts_μm[i,1], verts_μm[i,2], z0) for i in 1:nv]
    edges = [EGADS.make_edge_from_nodes(ctx, nodes[i], nodes[mod1(i+1,nv)]) for i in 1:nv]
    loop = EGADS.make_wire_loop(ctx, edges)
    face = EGADS.make_plane_face(ctx, loop)
    src = EGADS.make_face_body(ctx, face)
    return EGADS.extrude(src, height, [0.0, 0.0, 1.0])
end

function build_resonator_solids(rp)
    # rp = [bend_radius, trace_width, gap_width, straight_1, straight_2, straight_3]
    r = Float64(rp[1]) * 1μm
    tw = Float64(rp[2]) * 1μm
    gw = Float64(rp[3]) * 1μm
    s1 = Float64(rp[4]) * 1μm
    s2 = Float64(rp[5]) * 1μm
    s3 = Float64(rp[6]) * 1μm
    style = Paths.SimpleCPW(tw, gw)
    p = Path(Point(0.0μm, 0.0μm); α0=0.0°)
    straight!(p, s1, style)
    turn!(p, -90°, r)
    straight!(p, s2)
    turn!(p, 180°, r)
    straight!(p, s3)

    cs = CoordinateSystem("reson", nm)
    place!(cs, p)

    k = DLE.EGADSKernel(enable_tape=true)
    sm = SolidModel("reson", k)
    render!(sm, cs; map_meta=m -> string(DeviceLayout.layer(m)))

    DLE.start_recording!(k.geometry_tape)
    for (_, pg) in DeviceLayout.SolidModels.dimgroupdict(sm, 2)
        DeviceLayout.SolidModels.extrude_z!(pg, GH * 1μm)
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

# ─── Chain contraction ────────────────────────────────────────────────────────

function compute_chain(result_faces, prov, palace_mesh, dfdx, is_boundary,
                        bdy_map, sfv, gap_height, source_origin; gap_body_indices)
    dir = [0.0, 0.0, 1.0]
    nv = size(palace_mesh.vertices, 1)
    dxdp_bdy = zeros(Float64, nv, 3)
    n_gap = 0

    for (vii, fi, pt) in bdy_map
        if haskey(prov, fi) && prov[fi].body_idx in gap_body_indices
            α = clamp(dot(pt - source_origin, dir) / gap_height, 0.0, 1.0)
            pos_src = pt .- (α * gap_height) .* dir
            vel = interpolate_source_velocity(sfv, pos_src[1], pos_src[2])
            dxdp_bdy[vii, :] .= vel
            if norm(vel) > 1e-10; n_gap += 1; end
        end
    end

    dxdp_full = propagate_boundary_sensitivity(palace_mesh, dxdp_bdy, is_boundary)
    dfdp = dot(vec(dfdx), vec(dxdp_full))
    return (dfdp=dfdp, n_gap=n_gap, dxdp_full=dxdp_full)
end

# ─── Taylor test ─────────────────────────────────────────────────────────────

function run_taylor(palace, palace_mesh, freq, dfdp, dxdp_full, nfaces, label)
    nv = size(palace_mesh.vertices, 1)
    epsilons = [0.1, 0.05, 0.025]
    remainders = Float64[]
    for (i, ε) in enumerate(epsilons)
        pm = "mesh/ft10_$(label)_$(i).mesh"
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
        pc = "ft10_$(label)_$(i).json"; po = "postpro/ft10_$(label)_$(i)"
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

    # Default parameter values
    cap_length_μm = 620.0
    cap_gap_μm = 30.0
    cap_width_μm = 24.0
    junction_gap_μm = 12.0
    bend_radius_μm = 50.0
    trace_width_μm = 10.0
    gap_width_μm = 6.0
    straight_1_μm = 200.0
    straight_2_μm = 300.0
    straight_3_μm = 300.0

    island_params = [cap_length_μm, cap_gap_μm, cap_width_μm, junction_gap_μm]
    reson_params = [bend_radius_μm, trace_width_μm, gap_width_μm,
                    straight_1_μm, straight_2_μm, straight_3_μm]
    gap_height = GH

    println("═══ Full Transmon Sensitivity: 10 Parameters ═══")
    for (i, (name, val)) in enumerate(zip(PARAM_NAMES,
            vcat(island_params, reson_params)))
        @printf("  [%2d] %-12s = %7.1f μm\n", i, name, val)
    end

    # ─── ForwardDiff: island and resonator Jacobians ──────────────────────────

    println("\nStep 1: ForwardDiff Jacobians...")

    # Island: 4-param Jacobian
    island_coords = island_gap_coords_μm(island_params)
    J_island = ForwardDiff.jacobian(island_gap_coords_μm, island_params)
    nv_island = length(island_coords) ÷ 2
    verts_island = coords_to_vertices(island_coords)
    println("  Island: $nv_island vertices, 4-param Jacobian ($(size(J_island)))")

    # Resonator: 6-param Jacobian
    reson_coords = resonator_gap_coords_μm(reson_params)
    J_reson = ForwardDiff.jacobian(resonator_gap_coords_μm, reson_params)
    nv_reson = length(reson_coords) ÷ 2
    verts_reson = coords_to_vertices(reson_coords)
    println("  Resonator: $nv_reson vertices, 6-param Jacobian ($(size(J_reson)))")

    # Quick sanity: check Jacobian columns vs central FD
    # δ=1e-4 avoids cancellation error for polygon union operations
    println("\n  Jacobian sanity check (vs central FD, δ=1e-4):")
    for (name, func, params, J) in [("island", island_gap_coords_μm, island_params, J_island),
                                      ("reson", resonator_gap_coords_μm, reson_params, J_reson)]
        for j in 1:size(J, 2)
            δ = 1e-4
            pp = copy(params); pp[j] += δ
            pm = copy(params); pm[j] -= δ
            fd_col = (func(pp) - func(pm)) / (2δ)
            err = maximum(abs.(J[:, j] - fd_col))
            @printf("    %s param %d: max|J-FD| = %.2e %s\n", name, j, err, err < 1e-3 ? "✓" : "✗")
        end
    end

    # ─── Build EGADS geometry ─────────────────────────────────────────────────

    println("\nStep 2: Building geometry...")

    # Resonator solids (native curves via EGADSKernel)
    reson = build_resonator_solids(reson_params)
    reson_solids = reson.solids
    ctx = DLE._ctx(reson.k)
    println("  Resonator: $(length(reson_solids)) native-curve solids")

    # Island solid (extruded rectangle)
    verts_island_placed = copy(verts_island)
    verts_island_placed[:, 2] .+= IOY
    island_solid = make_extruded_polygon_μm(ctx, verts_island_placed, 0.0, gap_height)
    println("  Island: 1 solid, $nv_island vertices")

    # Substrate bounding box
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
                         gap_height + 20.0)

    # Sequential booleans
    result = sub
    for body in all_bodies
        rv = EGADS.solid_boolean(result, body, EGADS.BoolOp.SUBTRACTION)
        isempty(rv) && continue
        result = rv[1]
    end

    result_faces = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    nfaces = length(result_faces)
    println("  Result: $nfaces faces")

    # ─── Mesh + Palace ────────────────────────────────────────────────────────

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

    mesh_file = "mesh/ft10_baseline.mesh"
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
  "Problem":{"Type":"Eigenmode","Verbose":2,"Output":"postpro/ft10_baseline"},
  "Model":{"Mesh":"$mesh_file","L0":1.0e-4},
  "Domains":{"Materials":[{"Attributes":[1],"Permeability":1.0,"Permittivity":1.0}]},
  "Boundaries":{"PEC":{"Attributes":$(collect(1:nfaces))}},
  "Solver":{"Order":1,"Device":"CPU",
    "Eigenmode":{"N":2,"Tol":1.0e-8,"Target":10.0,"Sensitivity":true,"MeshSensitivity":true},
    "Linear":{"Type":"Default","Tol":1.0e-8,"MaxIts":200}}}"""
    write("ft10_baseline.json", config_str)
    mkpath("postpro/ft10_baseline")
    flush(stdout)
    t0 = time()
    run(pipeline(`$palace ft10_baseline.json`; stdout=devnull, stderr=devnull))
    @printf("  Palace: %.1fs\n", time()-t0)

    freq = read_eigenfrequency("postpro/ft10_baseline/eig.csv"; mode=1)
    @printf("  f₀ = %.6f GHz\n", freq)

    palace_mesh = read_mesh("postpro/ft10_baseline/mesh_for_sensitivity.mesh")
    sol = read_solb("postpro/ft10_baseline/mesh_sensitivity.solb")
    dfdx = sol.data
    is_boundary = classify_boundary_vertices(palace_mesh)
    nv_palace = size(palace_mesh.vertices, 1)
    @printf("  Palace mesh: %d verts, %d boundary\n", nv_palace, count(is_boundary))

    # ─── Pre-compute shared data ──────────────────────────────────────────────

    # Parent bodies: [substrate, island, reson1, reson2, ...]
    parent_bodies = vcat([sub, island_solid], reson_solids)
    prov = match_provenance(result, parent_bodies)

    # Pre-compute boundary vertex → face mapping (inv_evaluate is expensive)
    println("\nStep 4: Pre-computing boundary mapping...")
    bdy_map = Vector{Tuple{Int, Int, Vector{Float64}}}()
    for vi in findall(is_boundary)
        pt = palace_mesh.vertices[vi, :]
        for fi in eachindex(result_faces)
            inv_res = EGADS.inv_evaluate(result_faces[fi], pt)
            if norm(inv_res.results[1:3] - pt) < 1e-4
                push!(bdy_map, (vi, fi, pt))
                break
            end
        end
    end
    @printf("  Mapped %d / %d boundary vertices\n", length(bdy_map), count(is_boundary))

    # Extrusion bases
    dir = [0.0, 0.0, 1.0]
    island_faces = EGADS.get_body_topos(island_solid, EGADS.OClass.FACE)
    island_base = find_extrusion_base(island_faces, dir)

    all_reson_faces = vcat([EGADS.get_body_topos(s, EGADS.OClass.FACE) for s in reson_solids]...)
    reson_base = find_extrusion_base(all_reson_faces, dir)

    # Body index sets: island=2, resonators=3:(2+N)
    island_idx = Set([2])
    reson_idx = Set(3:(2 + length(reson_solids)))

    # Update island SFV positions to match placed coordinates
    vi_placed = verts_island_placed

    # ─── Chain contraction + Taylor tests for all 10 params ───────────────────

    println("\nStep 5: Chain contraction + Taylor tests (10 parameters)...\n")

    all_rates = Vector{Vector{Float64}}()
    all_dfdp = Float64[]

    for pi in 1:10
        name = PARAM_NAMES[pi]
        if pi <= 4
            # Island param: Ji column pi
            sfv = SourceFaceVelocity(vi_placed, jacobian_column_to_velocity(J_island, pi))
            base = island_base
            gidx = island_idx
        else
            # Resonator param: Jr column (pi - 4)
            sfv = SourceFaceVelocity(verts_reson, jacobian_column_to_velocity(J_reson, pi - 4))
            base = reson_base
            gidx = reson_idx
        end

        println("═══ [$pi] d/d($name) ═══")
        ch = compute_chain(result_faces, prov, palace_mesh, dfdx, is_boundary,
                            bdy_map, sfv, gap_height, base; gap_body_indices=gidx)
        @printf("  gap_vel=%d, df/dp=%.6e\n", ch.n_gap, ch.dfdp)
        push!(all_dfdp, ch.dfdp)

        rates = run_taylor(palace, palace_mesh, freq, ch.dfdp, ch.dxdp_full, nfaces, name)
        push!(all_rates, rates)
        println()
    end

    # ─── Summary ──────────────────────────────────────────────────────────────

    println("═══ Full Transmon — 10 Parameter Summary ═══")
    @printf("  f₀ = %.6f GHz\n\n", freq)
    n_pass = 0
    for pi in 1:10
        rates = all_rates[pi]
        max_r = isempty(rates) ? 0.0 : maximum(rates)
        pass = max_r > 1.5
        pass && (n_pass += 1)
        @printf("  [%2d] %-12s: df/dp=%+.4e  rates=[%s]  %s\n",
                pi, PARAM_NAMES[pi], all_dfdp[pi],
                join([@sprintf("%.2f", r) for r in rates], ","),
                pass ? "PASS ✓" : "FAIL ✗")
    end

    println("\n  Passed: $n_pass / 10")
    all_pass = n_pass == 10
    println("  Overall: $(all_pass ? "PASS ✓" : "FAIL ✗")")

    rm("postpro/ft10_baseline"; recursive=true, force=true)
    rm("ft10_baseline.json"; force=true); rm(mesh_file; force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
