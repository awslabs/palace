# Gradient-based transmon optimization: L-BFGS vs Nelder-Mead (10 parameters)
#
# 10 design parameters:
#   Island (4):    cap_length, cap_gap, cap_width, junction_gap
#   Resonator (6): bend_radius, trace_width, gap_width, straight_1, straight_2, straight_3
#
# Uses ForwardDiff Jacobians for exact df/dp, propagated through the full
# geometry → mesh → eigensolve → sensitivity chain.
#
# Each function evaluation = rebuild geometry + mesh + Palace eigensolve (~15-20s).
# The gradient comes from the sensitivity chain (same Palace solve, no extra cost).

using LinearAlgebra, Printf
using DeviceLayout, DeviceLayout.PreferredUnits
using ForwardDiff, Unitful
using EGADS
using TetGen
using Optim

DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

include(joinpath(@__DIR__, "shared", "palace_io.jl"))
include(joinpath(homedir(), "AWS", "DeviceLayout.jl-ad-sensitivity", "src", "forwarddiff_bridge.jl"))

# ─── Helpers ──────────────────────────────────────────────────────────────────

function face_plane(face)
    ri = EGADS.get_range(face); um=(ri.range[1]+ri.range[2])/2; vm=(ri.range[3]+ri.range[4])/2
    res = EGADS.evaluate(face, [um, vm]); n = cross(res.du, res.dv); nrm = norm(n)
    nrm < 1e-30 && return (normal=[0.0,0.0,0.0], offset=0.0, point=res.xyz)
    n = n / nrm; return (normal=n, offset=dot(n, res.xyz), point=res.xyz)
end

function match_provenance(rb, pbs; tol=1e-6)
    rfs = EGADS.get_body_topos(rb, EGADS.OClass.FACE)
    pd = Tuple{Int,Int,EGADS.EGObject,NamedTuple}[]
    for (bi, b) in enumerate(pbs)
        for (fi, pf) in enumerate(EGADS.get_body_topos(b, EGADS.OClass.FACE))
            push!(pd, (bi, fi, pf, face_plane(pf)))
        end
    end
    prov = Dict{Int, NamedTuple{(:body_idx,:face_idx), Tuple{Int,Int}}}()
    for (ri, rf) in enumerate(rfs)
        rp = face_plane(rf)
        for (bi, fi, _, pp) in pd
            ndot = abs(dot(rp.normal, pp.normal))
            if ndot > 1-tol && (abs(rp.offset-pp.offset)<tol || abs(rp.offset+pp.offset)<tol)
                prov[ri] = (body_idx=bi, face_idx=fi); break; end
        end
    end; return prov
end

function find_extrusion_base(faces, dir)
    mp = Inf; base = [0.0,0.0,0.0]
    for f in faces; rng = EGADS.get_range(f)
        c = EGADS.evaluate(f, [(rng.range[1]+rng.range[2])/2,(rng.range[3]+rng.range[4])/2]).xyz
        if dot(c,dir)<mp; mp=dot(c,dir); base=c; end
    end; return base
end

struct SourceFaceVelocity; vertices::Matrix{Float64}; velocity::Matrix{Float64}; end
function interpolate_source_velocity(sfv::SourceFaceVelocity, x::Float64, y::Float64)
    nv = size(sfv.vertices, 1); nv == 0 && return [0.0,0.0,0.0]
    md2 = Inf; bs = 0.0; bi = 1
    for i in 1:nv; j = mod1(i+1,nv)
        ax,ay = sfv.vertices[i,1],sfv.vertices[i,2]; bx,by = sfv.vertices[j,1],sfv.vertices[j,2]
        dx,dy = bx-ax,by-ay; lsq = dx^2+dy^2
        if lsq<1e-30; s=0.0; d2=(x-ax)^2+(y-ay)^2
        else; s=clamp(((x-ax)*dx+(y-ay)*dy)/lsq,0.0,1.0); px=ax+s*dx; py=ay+s*dy; d2=(x-px)^2+(y-py)^2; end
        if d2<md2; md2=d2; bs=s; bi=i; end
    end; j = mod1(bi+1,nv)
    return [(1-bs)*sfv.velocity[bi,1]+bs*sfv.velocity[j,1], (1-bs)*sfv.velocity[bi,2]+bs*sfv.velocity[j,2], 0.0]
end

# ─── Constants and parameter layout ──────────────────────────────────────────

const GH = 80.0   # gap extrusion height (μm)
const IOY = 50.0   # island Y-offset above resonator (μm)

const PARAM_NAMES = ["cap_length", "cap_gap", "cap_width", "junction_gap",
                     "bend_radius", "trace_width", "gap_width",
                     "straight_1", "straight_2", "straight_3"]
const P0 = [620.0, 30.0, 24.0, 12.0, 50.0, 10.0, 6.0, 200.0, 300.0, 300.0]

# Box constraints (μm)
const LB = [100.0, 5.0, 5.0, 2.0, 20.0, 2.0, 1.0, 50.0, 50.0, 50.0]
const UB = [1000.0, 100.0, 100.0, 50.0, 200.0, 30.0, 20.0, 500.0, 500.0, 500.0]

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

# ─── EGADS geometry builders ──────────────────────────────────────────────────

function mkpoly(ctx, v, z0, h)
    nv=size(v,1); ns=[EGADS.make_node(ctx,v[i,1],v[i,2],z0) for i in 1:nv]
    es=[EGADS.make_edge_from_nodes(ctx,ns[i],ns[mod1(i+1,nv)]) for i in 1:nv]
    l=EGADS.make_wire_loop(ctx,es); f=EGADS.make_plane_face(ctx,l); s=EGADS.make_face_body(ctx,f)
    EGADS.extrude(s, h, [0.0,0.0,1.0])
end

function build_reson_solids(rp)
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

    cs = CoordinateSystem("r", nm)
    place!(cs, p)
    k = DLE.EGADSKernel(enable_tape=true)
    sm = SolidModel("r", k)
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

# ─── Single evaluation (10-param) ────────────────────────────────────────────

function eval_transmon(params; palace_bin, with_grad=false)
    # params = [cap_length, cap_gap, cap_width, junction_gap,
    #           bend_radius, trace_width, gap_width, straight_1, straight_2, straight_3]
    ip = params[1:4]    # island params
    rp = params[5:10]   # resonator params

    # ForwardDiff Jacobians
    ic = island_gap_coords_μm(ip); nvi = length(ic) ÷ 2
    vi = coords_to_vertices(ic); vi[:, 2] .+= IOY
    Ji = with_grad ? ForwardDiff.jacobian(island_gap_coords_μm, ip) : nothing

    rc = resonator_gap_coords_μm(rp); nvr = length(rc) ÷ 2
    vr = coords_to_vertices(rc)
    Jr = with_grad ? ForwardDiff.jacobian(resonator_gap_coords_μm, rp) : nothing

    # Build EGADS geometry
    rs = build_reson_solids(rp); ctx = DLE._ctx(rs.k)
    is_body = mkpoly(ctx, vi, 0.0, GH)
    all_b = vcat([is_body], rs.solids)

    # Substrate bounding box
    bbm = [Inf, Inf, Inf]; bbx = [-Inf, -Inf, -Inf]
    for s in all_b
        bb = EGADS.get_bounding_box(s)
        bbm .= min.(bbm, [bb.xmin, bb.ymin, bb.zmin])
        bbx .= max.(bbx, [bb.xmax, bb.ymax, bb.zmax])
    end
    m = 50.0
    sub = EGADS.make_box(ctx, bbm[1]-m, bbm[2]-m, 0.0,
                         bbx[1]-bbm[1]+2m, bbx[2]-bbm[2]+2m, GH+20)

    # Sequential boolean: substrate - island - resonators
    result = sub
    for b in all_b
        rv = EGADS.solid_boolean(result, b, EGADS.BoolOp.SUBTRACTION)
        isempty(rv) && continue
        result = rv[1]
    end

    # Tessellate + TetGen
    rfs = EGADS.get_body_topos(result, EGADS.OClass.FACE); nf = length(rfs)
    tess = EGADS.tessellate(result; max_edge=30.0, sag=2.0, angle=12.0)
    np = EGADS.tess_status(tess).npts
    pts = zeros(Float64, 3, np)
    for fi in eachindex(rfs)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end
    tris = Int32[]; tmarks = Int32[]
    for fi in eachindex(rfs)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for t in 1:size(fm.tris, 1)
            for k in 1:3
                push!(tris, Int32(EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t, k])))
            end
            push!(tmarks, Int32(fi))
        end
    end
    inp = RawTetGenIO{Float64}()
    inp.pointlist = pts
    facetlist!(inp, reshape(tris, 3, length(tmarks)))
    inp.facetmarkerlist = tmarks
    tetout = tetrahedralize(inp, "pq1.2a50000")
    nvol = numberofpoints(tetout)

    # Write mesh
    mf = "mesh/opt_eval.mesh"; mkpath("mesh")
    open(mf, "w") do io
        println(io, "MeshVersionFormatted 2\n\nDimension 3\n")
        p = tetout.pointlist; println(io, "Vertices\n$nvol")
        for i in 1:nvol
            @printf(io, "%.15e %.15e %.15e 0\n", p[1,i], p[2,i], p[3,i])
        end
        te = tetout.tetrahedronlist; nt = size(te, 2)
        println(io, "\nTetrahedra\n$nt")
        for i in 1:nt
            @printf(io, "%d %d %d %d 1\n", te[1,i], te[2,i], te[3,i], te[4,i])
        end
        tf = tetout.trifacelist; tm = tetout.trifacemarkerlist; ntri = size(tf, 2)
        println(io, "\nTriangles\n$ntri")
        for i in 1:ntri
            @printf(io, "%d %d %d %d\n", tf[1,i], tf[2,i], tf[3,i], tm[i])
        end
        println(io, "\nEnd")
    end

    # Palace eigensolve
    sf = with_grad ? "true" : "false"
    write("opt_eval.json", """{"Problem":{"Type":"Eigenmode","Verbose":2,"Output":"postpro/opt_eval"},
"Model":{"Mesh":"$mf","L0":1.0e-4},
"Domains":{"Materials":[{"Attributes":[1],"Permeability":1.0,"Permittivity":1.0}]},
"Boundaries":{"PEC":{"Attributes":$(collect(1:nf))}},
"Solver":{"Order":1,"Device":"CPU",
"Eigenmode":{"N":2,"Tol":1.0e-8,"Target":10.0,"Sensitivity":$sf,"MeshSensitivity":$sf},
"Linear":{"Type":"Default","Tol":1.0e-8,"MaxIts":200}}}""")
    mkpath("postpro/opt_eval")
    run(pipeline(`$palace_bin opt_eval.json`; stdout=devnull, stderr=devnull))
    freq = read_eigenfrequency("postpro/opt_eval/eig.csv"; mode=1)

    # Sensitivity chain (10 parameters)
    grad = nothing
    if with_grad
        pm = read_mesh("postpro/opt_eval/mesh_for_sensitivity.mesh")
        sl = read_solb("postpro/opt_eval/mesh_sensitivity.solb")
        dfdx = sl.data
        isb = classify_boundary_vertices(pm)
        nvp = size(pm.vertices, 1)
        dir = [0.0, 0.0, 1.0]

        # Provenance (computed once for all params)
        pbs = vcat([sub, is_body], rs.solids)
        prov = match_provenance(result, pbs)

        # Extrusion bases
        ib = find_extrusion_base(EGADS.get_body_topos(is_body, EGADS.OClass.FACE), dir)
        arfs = vcat([EGADS.get_body_topos(s, EGADS.OClass.FACE) for s in rs.solids]...)
        rb = find_extrusion_base(arfs, dir)

        # Body index sets: island=2, resonators=3:(2+N)
        iidx = Set([2])
        ridx = Set(3:(2 + length(rs.solids)))

        # Pre-compute boundary vertex → face mapping (expensive inv_evaluate, do once)
        bdy_map = Vector{Tuple{Int, Int, Vector{Float64}}}()
        for vii in findall(isb)
            pt = pm.vertices[vii, :]
            for fi in eachindex(rfs)
                inv_res = EGADS.inv_evaluate(rfs[fi], pt)
                if norm(inv_res.results[1:3] - pt) < 1e-4
                    push!(bdy_map, (vii, fi, pt))
                    break
                end
            end
        end

        # Build SFVs and chain-contract for all 10 params
        grad = zeros(10)

        # Island params (1-4): Jacobian columns from Ji
        for pi in 1:4
            sfv = SourceFaceVelocity(vi, jacobian_column_to_velocity(Ji, pi))
            dxdp = zeros(Float64, nvp, 3)
            for (vii, fi, pt) in bdy_map
                if haskey(prov, fi) && prov[fi].body_idx in iidx
                    α = clamp(dot(pt - ib, dir) / GH, 0.0, 1.0)
                    ps = pt .- (α * GH) .* dir
                    v = interpolate_source_velocity(sfv, ps[1], ps[2])
                    dxdp[vii, :] .= v
                end
            end
            dxfull = propagate_boundary_sensitivity(pm, dxdp, isb)
            grad[pi] = dot(vec(dfdx), vec(dxfull))
        end

        # Resonator params (5-10): Jacobian columns from Jr
        for pi in 1:6
            sfv = SourceFaceVelocity(vr, jacobian_column_to_velocity(Jr, pi))
            dxdp = zeros(Float64, nvp, 3)
            for (vii, fi, pt) in bdy_map
                if haskey(prov, fi) && prov[fi].body_idx in ridx
                    α = clamp(dot(pt - rb, dir) / GH, 0.0, 1.0)
                    ps = pt .- (α * GH) .* dir
                    v = interpolate_source_velocity(sfv, ps[1], ps[2])
                    dxdp[vii, :] .= v
                end
            end
            dxfull = propagate_boundary_sensitivity(pm, dxdp, isb)
            grad[4 + pi] = dot(vec(dfdx), vec(dxfull))
        end
    end

    rm("postpro/opt_eval"; recursive=true, force=true)
    rm("opt_eval.json"; force=true); rm(mf; force=true)
    return (freq=freq, grad=grad)
end

# ─── Main ────────────────────────────────────────────────────────────────────

function main()
    cd(@__DIR__)
    palace = joinpath(@__DIR__, "..", "build", "palace-build", "palace-arm64.bin")
    @assert isfile(palace)

    p0 = copy(P0)
    p_scale = copy(P0)  # normalize so all params ≈ 1.0 at initial point
    x0 = p0 ./ p_scale  # = ones(10)

    # Normalized bounds
    lb = LB ./ p_scale
    ub = UB ./ p_scale

    println("═══ Transmon Optimization: 10 Parameters ═══")
    println("  Parameters:")
    for i in 1:10
        @printf("    [%2d] %-12s = %7.1f μm  (bounds: [%.0f, %.0f])\n",
                i, PARAM_NAMES[i], p0[i], LB[i], UB[i])
    end

    # Initial evaluation with gradient
    r0 = eval_transmon(p0; palace_bin=palace, with_grad=true)
    f_target = 0.97 * r0.freq  # 3% lower
    @printf("\n  f₀ = %.6f GHz, f_target = %.6f GHz\n", r0.freq, f_target)
    println("  ∇f (raw):")
    for i in 1:10
        @printf("    %-12s: %+.4e  (scaled: %+.4e)\n",
                PARAM_NAMES[i], r0.grad[i], r0.grad[i] * p_scale[i])
    end

    # L-BFGS with box constraints (gradient-based)
    println("\n─── L-BFGS (Fminbox) ───")
    n_lb = Ref(0)

    function f_lbfgs(x_norm)
        n_lb[] += 1
        x_raw = x_norm .* p_scale
        r = eval_transmon(x_raw; palace_bin=palace, with_grad=false)
        obj = (r.freq - f_target)^2
        @printf("  [%d] f=%.6f |Δf|=%.4e\n", n_lb[], r.freq, abs(r.freq - f_target))
        return obj
    end

    function g_lbfgs!(G, x_norm)
        x_raw = x_norm .* p_scale
        r = eval_transmon(x_raw; palace_bin=palace, with_grad=true)
        # Chain rule: ∂obj/∂x_norm = 2(f-f_target) · ∂f/∂x_raw · p_scale
        G .= 2(r.freq - f_target) .* r.grad .* p_scale
        return nothing
    end

    lb_res = optimize(f_lbfgs, g_lbfgs!, lb, ub, x0, Fminbox(LBFGS()),
                       Optim.Options(iterations=10, f_tol=1e-14, g_tol=1e-10))
    p_lb = Optim.minimizer(lb_res) .* p_scale
    f_lb_final = eval_transmon(p_lb; palace_bin=palace).freq

    # Nelder-Mead (derivative-free, baseline)
    println("\n─── Nelder-Mead ───")
    n_nm = Ref(0)

    function f_nm(x_norm)
        n_nm[] += 1
        x_raw = x_norm .* p_scale
        r = eval_transmon(x_raw; palace_bin=palace, with_grad=false)
        obj = (r.freq - f_target)^2
        if n_nm[] <= 20
            @printf("  [%d] f=%.6f |Δf|=%.4e\n", n_nm[], r.freq, abs(r.freq - f_target))
        end
        return obj
    end

    nm_res = optimize(f_nm, x0, NelderMead(),
                       Optim.Options(iterations=30, f_tol=1e-14))
    p_nm = Optim.minimizer(nm_res) .* p_scale
    f_nm_final = eval_transmon(p_nm; palace_bin=palace).freq

    # Results
    println("\n═══ Results ═══")
    @printf("  Target frequency: %.6f GHz\n\n", f_target)

    @printf("  L-BFGS:      f = %.6f GHz  |Δf| = %.4e  evals = %d\n",
            f_lb_final, abs(f_lb_final - f_target), n_lb[])
    @printf("  Nelder-Mead: f = %.6f GHz  |Δf| = %.4e  evals = %d\n",
            f_nm_final, abs(f_nm_final - f_target), n_nm[])
    println()

    # Per-parameter comparison
    println("  Parameter changes:")
    @printf("    %-12s  %10s  %10s  %10s\n", "Parameter", "Initial", "L-BFGS", "NM")
    for i in 1:10
        @printf("    %-12s  %10.1f  %10.1f  %10.1f\n",
                PARAM_NAMES[i], p0[i], p_lb[i], p_nm[i])
    end

    # Gradient magnitude ranking (which params matter most)
    println("\n  Gradient magnitudes (at initial point, normalized):")
    grad_norm = abs.(r0.grad .* p_scale)
    order = sortperm(grad_norm; rev=true)
    for i in order
        bar = repeat("█", round(Int, 30 * grad_norm[i] / maximum(grad_norm)))
        @printf("    %-12s: %+.4e  %s\n", PARAM_NAMES[i], r0.grad[i], bar)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
