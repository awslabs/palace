# ─── P2 curving with continuation ────────────────────────────────────────────
#
# Takes a linear TetGen volume mesh + its CAD classification (from
# `classify_tetgen_vertices`) and produces a valid CURVED order-2 mesh:
#
#   1. Chord-vertex snap: Steiner vertices on curved entities sit on the
#      linear chords; their classification records the true CAD point — the
#      snap displacement.
#   2. P2 mid-edge nodes: every mesh edge gets a midpoint node; edges whose
#      endpoints share a CAD entity (in the topology containment lattice
#      NODE ⊂ EDGE ⊂ FACE) get their mid-node projected onto that entity.
#   3. Continuation: displacements are applied along a homotopy
#      x(t) = x_straight + t·Δx with a per-element validity gate (sampled
#      Jacobian determinant of the isoparametric P2 map) and local Laplacian
#      relaxation of free interior nodes when the gate trips. The mesh is
#      valid at every accepted t — if full curving is unreachable the result
#      reports the highest valid t instead of returning a broken mesh.
#
# The same continuation core serves mesh MORPHING (boundary displacement from
# parameter velocities dx/dp instead of snap+curve); see `Δx` in the result.

"""
    curve_tetgen_mesh(out, cls, body;
                            max_bisect = 10,
                            relax_iters = 30,
                            detj_rel_tol = 0.05,
                            verbose = false)

Build a curved order-2 (P2) mesh from a linear TetGen output `out`, its
vertex classification `cls` (from [`classify_tetgen_vertices`](@ref)), and
the source EGADS `body`.

Returns a NamedTuple:

  * `coords::Matrix{Float64}`  — 3×(nv+nmid) node coordinates (vertices
    first, then mid-edge nodes), at the highest valid homotopy parameter.
  * `tets10::Matrix{Int}`      — 10×ntet connectivity, gmsh ordering
    (corners 1-4; mids (1,2),(2,3),(3,1),(1,4),(3,4),(2,4)).
  * `tris6::Matrix{Int}`       — 6×ntri boundary triangles, gmsh ordering.
  * `tri_markers::Vector`      — facet markers passed through from `out`.
  * `Δx::Matrix{Float64}`      — the full displacement field (3×nnodes) from
    the straight configuration; `coords = straight + t_final·Δx` plus any
    relaxation of free nodes.
  * `t_final::Float64`         — 1.0 when fully curved; < 1 if continuation
    stalled (mesh still valid at `t_final`).
  * `n_curved_mid`, `n_snapped`, `min_detj_ratio` — diagnostics.

Validity gate: at each trial t, every tet's isoparametric P2 Jacobian
determinant is sampled at 15 barycentric points (4 corners, 6 mid-edges,
4 face centroids, centroid) and must exceed `detj_rel_tol` × its straight
P1 value. This is the standard sampled-detJ heuristic (not a Bezier bound);
tighten `detj_rel_tol` or add sample points for paranoid cases.
"""
function curve_tetgen_mesh(out, cls, body::EGObject;
                                 max_bisect::Integer = 10,
                                 relax_iters::Integer = 30,
                                 detj_rel_tol::Real = 0.05,
                                 verbose::Bool = false)
    pts = out.pointlist
    tets = out.tetrahedronlist
    tris = out.trifacelist
    nv = size(pts, 2)
    ntet = size(tets, 2)

    faces = get_body_topos(body, OClass.FACE)
    bedges = get_body_topos(body, OClass.EDGE)
    bnodes = get_body_topos(body, OClass.NODE)

    # ── Topology containment: node → edges, edge → faces ────────────────────
    edge_index = Dict{Ptr{Cvoid}, Int}(e.ptr => i for (i, e) in enumerate(bedges))
    node_index = Dict{Ptr{Cvoid}, Int}(n.ptr => i for (i, n) in enumerate(bnodes))
    edge_faces = [Set{Int}() for _ in bedges]
    node_edges = [Set{Int}() for _ in bnodes]
    for (fi, f) in enumerate(faces)
        for e in get_body_topos(body, f, OClass.EDGE)
            EGADS.get_info(e).mtype == MType.DEGENERATE && continue
            ei = get(edge_index, e.ptr, 0)
            ei == 0 && continue
            push!(edge_faces[ei], fi)
        end
    end
    for (ei, e) in enumerate(bedges)
        for n in EGADS.get_topology(e).children
            ni = get(node_index, n.ptr, 0)
            ni == 0 && continue
            push!(node_edges[ni], ei)
        end
    end

    # ── Per-vertex classification lookup: (dim, idx) ─────────────────────────
    # dim: 0 = NODE, 1 = EDGE, 2 = FACE, -1 = interior/unclassified.
    vdim = fill(Int8(-1), nv)
    vidx = zeros(Int, nv)
    for r in eachrow(cls.on_nodes)
        vdim[r[1]] = 0; vidx[r[1]] = r[2]
    end
    for r in eachrow(cls.on_edges)
        vi = Int(r[1]); vdim[vi] = 1; vidx[vi] = Int(r[2])
    end
    for r in eachrow(cls.on_faces)
        vi = Int(r[1]); vdim[vi] = 2; vidx[vi] = Int(r[2])
    end

    "Edge ids containing vertex v's entity (∅ for FACE/interior)."
    v_edges(vi) = vdim[vi] == 1 ? Set{Int}((vidx[vi],)) :
                  vdim[vi] == 0 ? node_edges[vidx[vi]] : Set{Int}()
    "Face ids containing vertex v's entity (∅ for interior)."
    function v_faces(vi)
        vdim[vi] == 2 && return Set{Int}((vidx[vi],))
        vdim[vi] == 1 && return edge_faces[vidx[vi]]
        if vdim[vi] == 0
            s = Set{Int}()
            for ei in node_edges[vidx[vi]]
                union!(s, edge_faces[ei])
            end
            return s
        end
        return Set{Int}()
    end

    # ── Straight node array + vertex snap displacements ──────────────────────
    coords = zeros(3, nv)  # extended below
    Δ = zeros(3, nv)
    @inbounds for i in 1:nv, d in 1:3
        coords[d, i] = pts[d, i]
    end
    n_snapped = 0
    for r in eachrow(cls.on_edges)
        vi = Int(r[1])
        ev = try EGADS.evaluate(bedges[Int(r[2])], r[3]) catch; nothing end
        ev === nothing && continue
        Δ[:, vi] .= ev.xyz .- coords[:, vi]
        r[4] > 1e-12 && (n_snapped += 1)
    end
    for r in eachrow(cls.on_faces)
        vi = Int(r[1])
        ev = try EGADS.evaluate(faces[Int(r[2])], [r[3], r[4]]) catch; nothing end
        ev === nothing && continue
        Δ[:, vi] .= ev.xyz .- coords[:, vi]
        r[5] > 1e-12 && (n_snapped += 1)
    end

    # ── Mid-edge nodes: enumerate unique tet edges ────────────────────────────
    local_edges = ((1, 2), (2, 3), (3, 1), (1, 4), (3, 4), (2, 4))  # gmsh tet10 order
    edge_mid = Dict{Tuple{Int, Int}, Int}()
    mid_coords = Vector{NTuple{3, Float64}}()
    mid_delta = Vector{NTuple{3, Float64}}()
    mid_constrained = Vector{Bool}()   # CAD-projected (true) vs blend (false)
    n_curved_mid = 0

    # Only edges of constrained boundary facets may be projected onto CAD
    # entities. A through-volume chord whose endpoints both touch a curved
    # face must NOT have its mid pulled onto the surface — that flattens
    # boundary slivers onto the wall (zero/negative detJ) and freezes the
    # very nodes whose freedom the validity cure depends on.
    boundary_edge = Set{Tuple{Int, Int}}()
    for j in 1:size(tris, 2)
        a, b, c = Int(tris[1, j]), Int(tris[2, j]), Int(tris[3, j])
        push!(boundary_edge, minmax(a, b), minmax(b, c), minmax(c, a))
    end

    function _mid_node(a::Int, b::Int)
        key = a < b ? (a, b) : (b, a)
        mi = get(edge_mid, key, 0)
        mi != 0 && return mi
        # Straight midpoint of the SNAPPED endpoints — chord vertices will move
        # to their CAD positions at t=1, so the straight P2 reference uses the
        # original linear coordinates and the mid displacement is measured from
        # the original midpoint toward the curved target.
        mx = ((pts[1, a] + pts[1, b]) / 2, (pts[2, a] + pts[2, b]) / 2,
              (pts[3, a] + pts[3, b]) / 2)
        target = mx
        # Candidate shared CAD entities: edges AND faces. Selection is by
        # minimum projection displacement — this self-resolves the ambiguous
        # case of two vertices on the same CAD edge joined by a face-interior
        # chord (the face projection moves the mid ~0, the curve projection
        # would wrongly pull it to the boundary). A cap at 30% of the mesh
        # edge length rejects pathological projections outright (a real
        # curved segment's sagitta is far below it under any sane
        # tessellation), falling back to the endpoint-displacement blend.
        se = key in boundary_edge ? intersect(v_edges(a), v_edges(b)) : Set{Int}()
        sf = key in boundary_edge ? intersect(v_faces(a), v_faces(b)) : Set{Int}()
        elen = hypot(pts[1, a] - pts[1, b], pts[2, a] - pts[2, b],
                     pts[3, a] - pts[3, b])
        best_d = Inf
        if !isempty(se) || !isempty(sf)
            for ei in se
                r = try EGADS.inv_evaluate(bedges[ei], Float64[mx...]) catch; nothing end
                r === nothing && continue
                d = hypot(r.results[1] - mx[1], r.results[2] - mx[2], r.results[3] - mx[3])
                d < best_d && (best_d = d; target = (r.results[1], r.results[2], r.results[3]))
            end
            for fi in sf
                r = try EGADS.inv_evaluate(faces[fi], Float64[mx...]) catch; nothing end
                r === nothing && continue
                d = hypot(r.results[1] - mx[1], r.results[2] - mx[2], r.results[3] - mx[3])
                d < best_d && (best_d = d; target = (r.results[1], r.results[2], r.results[3]))
            end
        end
        if isfinite(best_d) && best_d <= 0.3 * elen
            best_d > 1e-12 && (n_curved_mid += 1)
            push!(mid_constrained, true)
        else
            # Interior, cross-entity, or capped edge: mid follows the average
            # of the endpoint displacements (keeps elements near snapped
            # boundaries from degenerating — the standard isoparametric blend).
            target = (mx[1] + (Δ[1, a] + Δ[1, b]) / 2,
                      mx[2] + (Δ[2, a] + Δ[2, b]) / 2,
                      mx[3] + (Δ[3, a] + Δ[3, b]) / 2)
            push!(mid_constrained, false)
        end
        push!(mid_coords, mx)
        push!(mid_delta, (target[1] - mx[1], target[2] - mx[2], target[3] - mx[3]))
        mi = nv + length(mid_coords)
        edge_mid[key] = mi
        return mi
    end

    tets10 = zeros(Int, 10, ntet)
    @inbounds for j in 1:ntet
        for c in 1:4
            tets10[c, j] = Int(tets[c, j])
        end
        for (k, (la, lb)) in enumerate(local_edges)
            tets10[4 + k, j] = _mid_node(Int(tets[la, j]), Int(tets[lb, j]))
        end
    end

    ntri = size(tris, 2)
    tris6 = zeros(Int, 6, ntri)
    @inbounds for j in 1:ntri
        a, b, c = Int(tris[1, j]), Int(tris[2, j]), Int(tris[3, j])
        tris6[1, j] = a; tris6[2, j] = b; tris6[3, j] = c
        tris6[4, j] = _mid_node(a, b)
        tris6[5, j] = _mid_node(b, c)
        tris6[6, j] = _mid_node(c, a)
    end

    nmid = length(mid_coords)
    nnodes = nv + nmid
    x0 = zeros(3, nnodes)
    Δx = zeros(3, nnodes)
    x0[:, 1:nv] .= coords
    Δx[:, 1:nv] .= Δ
    @inbounds for m in 1:nmid
        x0[:, nv + m] .= mid_coords[m]
        Δx[:, nv + m] .= mid_delta[m]
    end

    # Fixed nodes are the CAD-CONSTRAINED ones: classified vertices (their
    # position is the CAD point) and CAD-projected mid-nodes. Everything else
    # — interior vertices AND blend mids — may relax. In particular the
    # through-volume mid-edges of boundary slivers (tets with all corners on
    # a curved surface) must stay free: relaxing them is the only way such an
    # element survives full curving.
    is_fixed = falses(nnodes)
    for i in 1:nv
        is_fixed[i] = vdim[i] != Int8(-1)
    end
    for m in 1:nmid
        is_fixed[nv + m] = mid_constrained[m]
    end

    # Harmonic extension: propagate the boundary displacement smoothly into
    # the interior (Jacobi sweeps of Δx over free nodes). Without it, interior
    # nodes adjacent to snapped/curved boundary nodes stay frozen and the
    # first element layer gets pinched — the continuation then stalls short of
    # t = 1. This is the same propagation a parameter-morph uses for dx/dp.
    _harmonic_extend!(Δx, tets10, is_fixed)

    # ── Continuation ──────────────────────────────────────────────────────────
    x = copy(x0)
    t_done = 0.0
    step = 1.0
    straightened = Set{Int}()
    min_detj_ratio = Inf
    while t_done < 1.0 - 1e-12
        t_try = min(1.0, t_done + step)
        trial = x .+ (t_try - t_done) .* Δx   # x already holds x0 + t_done·Δx (+relaxation)
        bad = _invalid_p2_tets(trial, tets10, detj_rel_tol)
        if isempty(bad)
            x = trial
            t_done = t_try
            step = min(1.0, 2 * step)
            continue
        end
        if step > 1.0 / (1 << max_bisect)
            step /= 2
            continue
        end
        # Local relaxation: hill-climb free nodes around the offending
        # elements at the trial configuration, then re-gate.
        relaxed = _relax_around!(trial, tets10, bad, is_fixed, relax_iters)
        bad2 = _invalid_p2_tets(trial, tets10, detj_rel_tol)
        if relaxed > 0 && isempty(bad2)
            x = trial
            t_done = t_try
            step = 1.0 / (1 << (max_bisect - 2))
            continue
        end
        # Straightening fallback (the standard last resort, cf. gmsh): zero
        # the displacement of the offending elements' nodes — EXCEPT snapped
        # chord vertices, whose motion onto the CAD is not optional (shared
        # with neighboring elements and required for geometric fidelity).
        # Only mid-node curvature is sacrificed, and only locally: the
        # offending elements keep straight(er) edges while the rest of the
        # mesh curves fully. Geometry error stays bounded by the chord error
        # the element already had in the linear mesh.
        n_straightened_nodes = 0
        for j in bad2, k in 5:10
            n = tets10[k, j]
            if Δx[1, n] != 0.0 || Δx[2, n] != 0.0 || Δx[3, n] != 0.0
                trial[:, n] .= x0[:, n]   # back to the straight position
                Δx[:, n] .= 0.0           # and never displace it again
                is_fixed[n] = false        # …but relaxation may move it
                n_straightened_nodes += 1
                push!(straightened, n)
            end
        end
        if n_straightened_nodes > 0
            _relax_around!(trial, tets10, bad2, is_fixed, relax_iters)
            bad3 = _invalid_p2_tets(trial, tets10, detj_rel_tol)
            if isempty(bad3)
                x = trial
                t_done = t_try
                step = 1.0 / (1 << (max_bisect - 2))
                continue
            end
        end
        verbose && @warn "curve_tetgen_mesh: continuation stalled" t_done step length(bad)
        break
    end
    # Final quality figure.
    min_detj_ratio = _min_detj_ratio(x, tets10)

    verbose && println("curve_tetgen_mesh: t_final=", round(t_done; digits = 6),
                       "  nodes=", nnodes, " (", nv, "+", nmid, " mid)",
                       "  curved mids=", n_curved_mid, "  snapped verts=", n_snapped,
                       "  min detJ ratio=", round(min_detj_ratio; digits = 4))
    return (coords = x, tets10 = tets10, tris6 = tris6,
            tri_markers = collect(out.trifacemarkerlist),
            Δx = Δx, t_final = t_done,
            n_curved_mid = n_curved_mid, n_snapped = n_snapped,
            n_straightened = length(straightened),
            min_detj_ratio = min_detj_ratio)
end

# ─── P2 tet Jacobian sampling ────────────────────────────────────────────────

# Reference barycentric sample points: corners, edge mids, face centroids, centroid.
const _P2_SAMPLES = let
    pts = Vector{NTuple{4, Float64}}()
    push!(pts, (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1))
    for (i, j) in ((1, 2), (2, 3), (3, 1), (1, 4), (3, 4), (2, 4))
        L = zeros(4); L[i] = 0.5; L[j] = 0.5
        push!(pts, Tuple(L))
    end
    for f in ((1, 2, 3), (1, 2, 4), (2, 3, 4), (1, 3, 4))
        L = zeros(4); for i in f; L[i] = 1 / 3; end
        push!(pts, Tuple(L))
    end
    push!(pts, (0.25, 0.25, 0.25, 0.25))
    pts
end

# ∇N for the 10-node tet at barycentric L (gmsh node order), reference coords
# ξ = (L2, L3, L4): ∇L1 = (-1,-1,-1), ∇L2 = e1, ∇L3 = e2, ∇L4 = e3.
const _GRAD_L = ((-1.0, -1.0, -1.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

function _p2_grads(L::NTuple{4, Float64})
    g = Matrix{Float64}(undef, 3, 10)
    for c in 1:4   # corner nodes: N = L(2L-1), ∇N = (4L-1)∇L
        s = 4 * L[c] - 1
        for d in 1:3
            g[d, c] = s * _GRAD_L[c][d]
        end
    end
    for (k, (i, j)) in enumerate(((1, 2), (2, 3), (3, 1), (1, 4), (3, 4), (2, 4)))
        for d in 1:3   # edge nodes: N = 4 Li Lj
            g[d, 4 + k] = 4 * (L[i] * _GRAD_L[j][d] + L[j] * _GRAD_L[i][d])
        end
    end
    return g
end

const _P2_SAMPLE_GRADS = [_p2_grads(L) for L in _P2_SAMPLES]

# 3×3 determinant (avoids a LinearAlgebra dependency in the extension).
_det3(m::AbstractMatrix) =
    m[1, 1] * (m[2, 2] * m[3, 3] - m[2, 3] * m[3, 2]) -
    m[1, 2] * (m[2, 1] * m[3, 3] - m[2, 3] * m[3, 1]) +
    m[1, 3] * (m[2, 1] * m[3, 2] - m[2, 2] * m[3, 1])

"Sampled min detJ of one P2 tet, normalized by its straight P1 detJ."
function _tet10_min_detj_ratio(x::Matrix{Float64}, tets10::Matrix{Int}, j::Int)
    # Straight P1 reference determinant from the corners.
    a, b, c, d = tets10[1, j], tets10[2, j], tets10[3, j], tets10[4, j]
    e1 = (x[1, b] - x[1, a], x[2, b] - x[2, a], x[3, b] - x[3, a])
    e2 = (x[1, c] - x[1, a], x[2, c] - x[2, a], x[3, c] - x[3, a])
    e3 = (x[1, d] - x[1, a], x[2, d] - x[2, a], x[3, d] - x[3, a])
    J1 = e1[1] * (e2[2] * e3[3] - e2[3] * e3[2]) -
         e1[2] * (e2[1] * e3[3] - e2[3] * e3[1]) +
         e1[3] * (e2[1] * e3[2] - e2[2] * e3[1])
    J1 == 0 && return -Inf
    mn = Inf
    Jm = Matrix{Float64}(undef, 3, 3)
    for g in _P2_SAMPLE_GRADS
        fill!(Jm, 0.0)
        @inbounds for nloc in 1:10
            ncol = tets10[nloc, j]
            for d2 in 1:3, d1 in 1:3
                Jm[d1, d2] += x[d1, ncol] * g[d2, nloc]
            end
        end
        r = _det3(Jm) / J1
        r < mn && (mn = r)
    end
    return mn
end

function _invalid_p2_tets(x, tets10, rel_tol)
    bad = Int[]
    for j in 1:size(tets10, 2)
        _tet10_min_detj_ratio(x, tets10, j) < rel_tol && push!(bad, j)
    end
    return bad
end

function _min_detj_ratio(x, tets10)
    mn = Inf
    for j in 1:size(tets10, 2)
        r = _tet10_min_detj_ratio(x, tets10, j)
        r < mn && (mn = r)
    end
    return mn
end

"""
Cure invalid elements by hill-climbing their FREE nodes: each free node of a
bad tet is nudged along ±axis directions, keeping moves that raise the worst
sampled detJ ratio over ALL elements incident to that node. Direction-free —
works for curved slivers (all corners on a curved surface) where Laplacian
averaging moves nodes the wrong way. Returns the number of nodes moved.
"""
function _relax_around!(x, tets10, bad, is_fixed, iters)
    affected = Set{Int}()
    for j in bad, k in 1:10
        is_fixed[tets10[k, j]] || push!(affected, tets10[k, j])
    end
    isempty(affected) && return 0
    incident = Dict{Int, Vector{Int}}()
    for j in 1:size(tets10, 2), k in 1:10
        n = tets10[k, j]
        n in affected || continue
        push!(get!(incident, n, Int[]), j)
    end
    # Local length scale per node: shortest incident corner edge.
    function _hscale(n)
        h = Inf
        for j in incident[n], (a, b) in ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
            na, nb = tets10[a, j], tets10[b, j]
            h = min(h, hypot(x[1, na] - x[1, nb], x[2, na] - x[2, nb],
                             x[3, na] - x[3, nb]))
        end
        return h
    end
    worst(n) = minimum(_tet10_min_detj_ratio(x, tets10, j) for j in incident[n])
    moved = Set{Int}()
    for _ in 1:iters
        improved = false
        for n in affected
            base = worst(n)
            base > 0.2 && continue   # healthy enough
            h = 0.1 * _hscale(n)
            best_gain, best_d, best_s = 0.0, 0, 0.0
            for d in 1:3, s in (h, -h)
                x[d, n] += s
                gain = worst(n) - base
                x[d, n] -= s
                gain > best_gain && (best_gain = gain; best_d = d; best_s = s)
            end
            if best_d != 0
                x[best_d, n] += best_s
                push!(moved, n)
                improved = true
            end
        end
        improved || break
    end
    return length(moved)
end

# Harmonic (Jacobi) extension of the displacement field into free nodes over
# the P2 node graph. Fixed nodes keep their Δ; free nodes converge to the
# average of their neighbors — a smooth interior continuation of the boundary
# motion (discrete Laplace with Dirichlet data on the constrained set).
function _harmonic_extend!(Δx::Matrix{Float64}, tets10::Matrix{Int},
                           is_fixed::AbstractVector{Bool};
                           max_sweeps::Integer = 400, rtol::Real = 1e-3)
    nnodes = size(Δx, 2)
    nbrs = [Set{Int}() for _ in 1:nnodes]
    for j in 1:size(tets10, 2), k1 in 1:10
        n1 = tets10[k1, j]
        is_fixed[n1] && continue
        for k2 in 1:10
            k2 == k1 && continue
            push!(nbrs[n1], tets10[k2, j])
        end
    end
    free = [i for i in 1:nnodes if !is_fixed[i] && !isempty(nbrs[i])]
    isempty(free) && return Δx
    scale = maximum(abs, Δx)
    scale == 0 && return Δx
    buf = zeros(3, length(free))
    for sweep in 1:max_sweeps
        worst = 0.0
        for (fi, n1) in enumerate(free)
            mx = my = mz = 0.0
            for n2 in nbrs[n1]
                mx += Δx[1, n2]; my += Δx[2, n2]; mz += Δx[3, n2]
            end
            c = length(nbrs[n1])
            buf[1, fi] = mx / c; buf[2, fi] = my / c; buf[3, fi] = mz / c
        end
        for (fi, n1) in enumerate(free)
            worst = max(worst, abs(buf[1, fi] - Δx[1, n1]),
                        abs(buf[2, fi] - Δx[2, n1]), abs(buf[3, fi] - Δx[3, n1]))
            Δx[1, n1] = buf[1, fi]; Δx[2, n1] = buf[2, fi]; Δx[3, n1] = buf[3, fi]
        end
        worst < rtol * scale && break
    end
    return Δx
end

# ─── Order-2 Gmsh MSH v2.2 writer ────────────────────────────────────────────

"""
    write_msh2_p2(path, curved;
                        tet_attr = 1,
                        tri_attr = marker -> marker,
                        physical_names = nothing)

Write the result of [`curve_tetgen_mesh`](@ref) as an ASCII Gmsh MSH v2.2
order-2 mesh (10-node tets, type 11; 6-node boundary triangles, type 9) —
directly readable by Palace/MFEM.

  * `tet_attr`: integer (all tets one region) or a `Vector` of per-tet
    attributes (e.g. from `tetrahedronattributelist`).
  * `tri_attr`: function mapping a facet marker to the physical attribute
    written for that triangle; return `nothing` to drop the triangle.
  * `physical_names`: optional `Vector{Tuple{Int, Int, String}}` of
    `(dim, attr, name)` written to `\$PhysicalNames`.
"""
function write_msh2_p2(path::AbstractString, curved;
                             tet_attr = 1,
                             tri_attr = identity,
                             physical_names = nothing)
    x = curved.coords
    tets10 = curved.tets10
    tris6 = curved.tris6
    nnodes = size(x, 2)
    ntet = size(tets10, 2)
    open(path, "w") do io
        println(io, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat")
        if physical_names !== nothing
            println(io, "\$PhysicalNames\n", length(physical_names))
            for (dim, attr, name) in physical_names
                println(io, dim, " ", attr, " \"", name, "\"")
            end
            println(io, "\$EndPhysicalNames")
        end
        println(io, "\$Nodes\n", nnodes)
        for i in 1:nnodes
            println(io, i, " ", x[1, i], " ", x[2, i], " ", x[3, i])
        end
        println(io, "\$EndNodes")
        kept = Tuple{Int, Int}[]   # (triangle index, attribute)
        for j in 1:size(tris6, 2)
            a = tri_attr(Int(curved.tri_markers[j]))
            a === nothing && continue
            push!(kept, (j, a))
        end
        println(io, "\$Elements\n", ntet + length(kept))
        eid = 0
        for j in 1:ntet
            eid += 1
            attr = tet_attr isa Integer ? Int(tet_attr) : Int(tet_attr[j])
            print(io, eid, " 11 2 ", attr, " ", attr)
            for k in 1:10
                print(io, " ", tets10[k, j])
            end
            println(io)
        end
        for (j, a) in kept
            eid += 1
            print(io, eid, " 9 2 ", a, " ", a)
            for k in 1:6
                print(io, " ", tris6[k, j])
            end
            println(io)
        end
        println(io, "\$EndElements")
    end
    return path
end
