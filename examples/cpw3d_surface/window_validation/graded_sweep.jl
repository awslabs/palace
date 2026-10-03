# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# The graded cross-section sweep (`sweep = :graded`, USER decision 263; design: Q2 option (v),
# reference-quality-20261002/mesher-scoping/DESIGN.md section 4). Included by
# PolygonWindowMesh.jl after structured_band.jl; the tensor sweep stays the default.
#
# The tensor sweep extrudes EVERY plan triangle through EVERY z level, so the 10-nm band rows
# are swept through the 425-um far intervals (needles) and the r-spaced metal levels across the
# 30-um plan interior. Here the structured band's columns carry a GRADED (n, z) cross-section:
#
#   * stacks: Z_0 = the production z levels; row k = 1..K carries Z_k = thin(Z_{k-1}, alpha w_k)
#     (w_k = h_k - h_{k-1} the row's width; a level is kept iff it is a material interface —
#     the trench bottom, substrate surface and metal top of every plane, the box ends, the
#     backsides and, with two planes, the gap midpoint — or at least alpha w_k from the last
#     kept level and from the next interface), nested by construction; the region (Gmsh
#     interior) nodes carry the ladder Z_{K+j} = thin(Z_{K+j-1}, w_K 2^j) (no alpha: alpha
#     grades the band rows only) chosen by their plan size (the shortest incident plan edge,
#     at most the region mesh size), except the ring of region nodes adjacent to the band,
#     which carries Z_K (`region_ring`, decision 275);
#   * row ranges: row line k is active for z in [lo_k, hi_k] = the first level of Z_{k+1} beyond
#     the faces of the fabricated step (trench bottom, substrate surface, metal top) -/+
#     beta h_k, strictly nested (inner rows end first); row K spans the box. With two planes
#     the step is that of the planes the chain's edges belong to (a plan edge is a metal edge
#     of plane k when the conductor changes across it on that plane; a bump footprint edge
#     belongs to both planes), the union over the chain: the band of one plane's edges ends at
#     the gap midpoint at the latest (a gap too small for that is refused), a band whose
#     chain belongs to both planes — a coincident cross-plane run, a chain turning from one
#     plane's edge onto the other's, a bump — keeps every row active across the gap (the
#     design's bump rule). One cross-section per such variant; the stacks are shared;
#   * the cross-section per column: between two neighbouring active lines the cells are the
#     consecutive levels of the COARSER stack; the finer line's extra levels hang on one vertical
#     edge and the end node of the row that ended at the cell's bottom (started at its top)
#     hangs on that horizontal edge; a fan triangulates every cell without a degenerate
#     triangle — from that row-end node where there is one, else from the coarse line's node at
#     the horizontal edge farther from the substrate surface;
#   * the sweep between two consecutive columns of a chain (plan nodes P_0..P_n and Q_0..Q_m,
#     n and m rows) has three kinds of cells:
#       - the (0, k) cross-section cells (from the base line to row k beyond row k-1's range)
#         are swept as triangles: node (k, z) sits on P_min(k, n) and Q_min(k, m); coincident
#         nodes (a capped column: rows above m collapse onto the shorter column's top line; a
#         fan: the seven columns share the base) turn the prism into a pyramid or a tetrahedron,
#         the 3D analogue of the plan band's collapse / fan triangles;
#       - the (k-1, k) strips within row k-1's range are right prisms over the plan band cells
#         (the quad P_{k-1} P_k Q_k Q_{k-1} split by the diagonal through its lowest plan index,
#         or the collapse triangle P_{k-1} P_k Q_m) between consecutive levels of Z_k with the
#         rows' stacks hanging on the vertical edges (the region rule below); clamping the
#         hanging fan triangles of such a strip onto one column top would give three collinear
#         nodes, so the strips are not swept;
#       - the region: plan triangle x interval of the levels common to its corners' lists; the
#         finer nodes' extra levels hang on their vertical edges. The top node of a column with
#         n rows carries the list the band puts on it: Z_n within row n's range and Z_k where
#         row k is the innermost active row beyond (for a full column, Z_K). With two planes
#         the tops of chains of different variants are fine near their own plane and Z_K near
#         the other, so their lists are nested interval by interval, not globally: the common
#         list is the intersection and one corner is hanging-free in every interval;
#   * tetrahedra: every element is a fan from one apex over its boundary triangulation, and
#     every face is triangulated by a rule of the face alone, so neighbours always agree: a quad
#     takes the diagonal through its lowest-ORDERED node; a vertical face with hanging nodes is
#     a merge ladder from the substrate-surface side (a quad where both sides step together, a
#     triangle advancing the nearer side otherwise); the apex is the lowest-ordered among the
#     element's corners on hanging-free vertical edges (the fan over a quad through that
#     corner's diagonal is exact, the ladder from a hanging-free near corner is a fan from it).
#     The node ORDER used by these rules is (z level by distance to the nearest substrate
#     surface, then the plan node by band row DESCENDING, taller columns first, then the
#     welding order): the
#     apex of every cell is then its node nearest the metal, the construction is mirror-
#     symmetric about the surface, and at a row end next to a capped column the apex is the
#     row-end node — the swept corner cell there is a pyramid over the vertical face of the
#     band-top edge whose side quad must be split through that node (with the welding order
#     the rule picks the other diagonal, and below the surface the three diagonals twist into
#     a Schoenhardt prism: not tetrahedralisable). The written mesh uses the tensor sweep's
#     numbering (level - 1) x plan nodes + plan index (remapped at the end), so the surface
#     tagging is shared.
#   * backstop: every tetrahedron is checked positive where it is built (naming the cell), and
#     the swept volume per material must equal the plan's analytic volume to 1e-9 (an overlap
#     or a gap of a single cell is refused).
#
# Scope (decision 263): one or two planes with bumps (M3). alpha and beta are dimensionless
# (defaults 1 and 3): alpha = 1 keeps the r x r edge cell of the recorded family; the
# z-coarsened outer rows are a deliberate deviation to be shown harmless by the acceptance
# solves.

const DEFAULT_GRADED_ALPHA = 1.0
const DEFAULT_GRADED_BETA = 3.0
const LEVEL_TOLERANCE_UM = 1.0e-9
# A row range end is nudged past the previous row's end by this before the level search (strict
# nesting), compared with a tolerance three decades below it.
const RANGE_NUDGE_UM = 1.0e-9
const RANGE_COMPARE_UM = 1.0e-12
# A region node whose plan size is within this (in log2) of a power of two of the last row's
# width takes the coarser ladder step.
const LADDER_LOG2_TOLERANCE = 1.0e-9
const VOLUME_CLOSURE_TOLERANCE = 1.0e-9

struct GradedStacks
    levels::Vector{Vector{Int}} # stack s = 0..K+J as sorted indices into the z levels
    spacings::Vector{Float64} # thinning spacing of stack s = 1..K+J
    rows::Int # K
    ranges::Vector{NTuple{2, Int}} # level indices [lo_k, hi_k] of row k = 1..K
end

is_interface(z::Float64, interfaces::Vector{Float64}) =
    any(abs(z - w) <= LEVEL_TOLERANCE_UM for w in interfaces)

# Nested thinning of a sorted level list: keep an interface always, another level only if it
# is at least `spacing` above the last kept level and below the next interface.
function thin_levels(
    zs::Vector{Float64},
    levels::Vector{Int},
    spacing::Float64,
    interfaces::Vector{Float64}
)
    kept = Int[levels[1]]
    for (position, index) in enumerate(levels)
        position == 1 && continue
        z = zs[index]
        if is_interface(z, interfaces)
            push!(kept, index)
            continue
        end
        next = findnext(
            j -> is_interface(zs[levels[j]], interfaces),
            eachindex(levels),
            position
        )
        room = next === nothing ? Inf : zs[levels[next]] - z
        if z - zs[kept[end]] >= spacing - LEVEL_TOLERANCE_UM &&
           room >= spacing - LEVEL_TOLERANCE_UM
            push!(kept, index)
        end
    end
    return kept
end

"""
    graded_stacks(zs, interfaces, heights, alpha, beta, step_faces, region_size_max;
                  range_cap=(-Inf, Inf)) -> GradedStacks

The nested stacks Z_0..Z_K of the band rows (heights h_1..h_K), the region ladder beyond
Z_K with spacings w_K 2^j up to the region size (no alpha), and the active range of every
row: the first level of Z_{k+1} at or beyond the fabricated step's faces (`step_faces` = the
lowest and highest of trench bottom, substrate surface and metal top over the planes the
band's edges belong to) -/+ beta h_k, strictly nested; row K spans the box. `range_cap`
bounds the limits of the rows k < K (two planes: the band of one plane's edges ends at the
gap midpoint at the latest, DESIGN.md 4.2); a row whose nested range would pass a cap is
refused (two rows would have to end on the same level). A row k < K whose range would reach
a box face is refused (the box must extend more than beta h_{K-1} beyond the step).
"""
function graded_stacks(
    zs::Vector{Float64},
    interfaces::Vector{Float64},
    heights::Vector{Float64},
    alpha::Float64,
    beta::Float64,
    step_faces::NTuple{2, Float64},
    region_size_max::Float64;
    range_cap::NTuple{2, Float64}=(-Inf, Inf)
)
    alpha > 0.0 || error("alpha must be positive")
    beta > 0.0 || error("beta must be positive")
    rows = length(heights)
    widths = [heights[k] - (k == 1 ? 0.0 : heights[k - 1]) for k = 1:rows]
    levels = Vector{Vector{Int}}()
    spacings = Float64[]
    push!(levels, collect(eachindex(zs)))
    for k = 1:rows
        push!(spacings, alpha * widths[k])
        push!(levels, thin_levels(zs, levels[end], spacings[end], interfaces))
    end
    # The region ladder beyond Z_K: spacing w_K 2^j up to the region size — the design's
    # "Z_K thinned to min(plan size, 30 um)" (DESIGN.md 4.1), without alpha, which scales the
    # BAND rows' z spacing only (M1 review MINOR-2).
    ladder = 0
    while widths[rows] * 2.0^(ladder + 1) <= region_size_max
        ladder += 1
        push!(spacings, widths[rows] * 2.0^ladder)
        push!(levels, thin_levels(zs, levels[end], spacings[end], interfaces))
    end
    ranges = Vector{NTuple{2, Int}}(undef, rows)
    for k = 1:rows
        if k == rows
            ranges[k] = (1, length(zs))
            continue
        end
        coarser = levels[k + 2]
        limit_lo = max(step_faces[1] - beta * heights[k], range_cap[1])
        limit_hi = min(step_faces[2] + beta * heights[k], range_cap[2])
        if k > 1
            limit_lo = min(limit_lo, zs[ranges[k - 1][1]] - RANGE_NUDGE_UM)
            limit_hi = max(limit_hi, zs[ranges[k - 1][2]] + RANGE_NUDGE_UM)
        end
        lo = findlast(i -> zs[i] <= limit_lo + RANGE_COMPARE_UM, coarser)
        hi = findfirst(i -> zs[i] >= limit_hi - RANGE_COMPARE_UM, coarser)
        (lo !== nothing && lo > 1 && hi !== nothing && hi < length(coarser)) || error(
            "Graded sweep: row $k of $rows would end on a box face (its range reaches " *
            "beyond [$limit_lo, $limit_hi] um); the box must extend more than beta h_(K-1) = " *
            "$(beta * heights[rows - 1]) um beyond the fabricated step on both sides"
        )
        (
            zs[coarser[lo]] >= range_cap[1] - RANGE_COMPARE_UM &&
            zs[coarser[hi]] <= range_cap[2] + RANGE_COMPARE_UM
        ) || error(
            "Graded sweep: row $k of $rows would end beyond the range cap $range_cap um " *
            "(the flip-chip gap is too small for beta h_k = $(beta * heights[k]) um: rows " *
            "$(k - 1) and $k would both end on the gap midpoint); lower beta or the row count"
        )
        ranges[k] = (coarser[lo], coarser[hi])
    end
    return GradedStacks(levels, spacings, rows, ranges)
end

# The innermost active row at level index i (row K beyond every other range).
function innermost_row(stacks::GradedStacks, i::Int)
    for k = 1:(stacks.rows)
        lo, hi = stacks.ranges[k]
        lo <= i <= hi && return k
    end
    return stacks.rows
end

"""
    column_top_levels(stacks) -> Vector{Vector{Int}}

The level list the band puts on the top node of a column with n = 1..K rows: Z_n within row
n's range, Z_k where row k > n is the innermost active row beyond it (the (0, k) cells collapse
onto the top line there). For n = K this is Z_K.
"""
function column_top_levels(stacks::GradedStacks)
    return [
        [
            i for i in stacks.levels[1] if
            i in stacks.levels[max(n, innermost_row(stacks, i)) + 1]
        ] for n = 1:(stacks.rows)
    ]
end

# ---------------------------------------------------------------------------------------------
# The (n, z) cross-section of a column: nodes (row, level index), fan triangles.

struct CrossSection
    nodes::Vector{NTuple{2, Int}} # (row 0..K, level index)
    triangles::Vector{NTuple{3, Int}}
    triangle_pair::Vector{NTuple{2, Int}} # (line, line) of every triangle's cell
    pair_triangles::Vector{Tuple{NTuple{2, Int}, Int}} # ((line, line), triangles) per pair
end

# `surface_distance[level]`: the distance of a level to the nearest substrate surface (the
# "near" side of every face rule).
function cross_section(
    stacks::GradedStacks,
    zs::Vector{Float64},
    heights::Vector{Float64},
    surface_distance::Vector{Float64}
)
    rows = stacks.rows
    node_index = Dict{NTuple{2, Int}, Int}()
    nodes = NTuple{2, Int}[]
    node(row, level) = get!(node_index, (row, level)) do
        push!(nodes, (row, level))
        return length(nodes)
    end
    position(row) = row == 0 ? 0.0 : heights[row]
    line_range(row) = row == 0 ? (1, length(zs)) : stacks.ranges[row]
    levels_in(row, la, lb) = [i for i in stacks.levels[row + 1] if la <= i <= lb]
    triangles = NTuple{3, Int}[]
    triangle_pair = NTuple{2, Int}[]
    pair_triangles = Tuple{NTuple{2, Int}, Int}[]
    orient2(a, b, c) = begin
        (ra, la), (rb, lb), (rc, lc) = nodes[a], nodes[b], nodes[c]
        (position(rb) - position(ra)) * (zs[lc] - zs[la]) -
        (zs[lb] - zs[la]) * (position(rc) - position(ra))
    end
    function fan!(apex, chain, pair)
        for i = 1:(length(chain) - 1)
            a, b, c = apex, chain[i], chain[i + 1]
            orient2(a, b, c) > 0.0 || ((b, c) = (c, b))
            orient2(a, b, c) > 0.0 || error(
                "Degenerate cross-section triangle (row, z) $(nodes[a]) $(nodes[b]) $(nodes[c])"
            )
            push!(triangles, (a, b, c))
            push!(triangle_pair, pair)
        end
    end
    # Cells between lines L and R (rows) over the level range [la, lb]; `bottom_hanging` /
    # `top_hanging`: the row whose end node hangs on the first cell's bottom / last cell's
    # top. A cell with a hanging row end is fanned FROM that node (every swept element of the
    # cell then contains it, and the node order makes it their apex, which is what keeps the
    # row-end cells next to a capped column star-shaped); any other cell is fanned from the
    # coarse line's node at the horizontal edge farther from the substrate surface (the swept
    # elements' apex is then their base node nearer the surface, on the diagonal of every
    # quad through it).
    function cells_between!(L, R, la, lb, bottom_hanging, top_hanging)
        lb > la || return
        before = length(triangles)
        levels_L, levels_R = levels_in(L, la, lb), levels_in(R, la, lb)
        # The coarser line (ties to the higher row, so the apex of a (0, k) cell is on line k).
        coarse, fine, cl, fl =
            length(levels_L) < length(levels_R) ? (L, R, levels_L, levels_R) :
            (R, L, levels_R, levels_L)
        issubset(cl, fl) ||
            error("Stacks of rows $L and $R are not nested in levels $la..$lb")
        (cl[1] == la && cl[end] == lb) ||
            error("Range ends $la..$lb are not levels of row $coarse")
        for ci = 1:(length(cl) - 1)
            l0, l1 = cl[ci], cl[ci + 1]
            fine_levels = [i for i in fl if l0 <= i <= l1]
            bottom = ci == 1 ? bottom_hanging : Int[]
            top = ci == length(cl) - 1 ? top_hanging : Int[]
            (isempty(bottom) || isempty(top)) ||
                error("Cross-section cell with hanging nodes on both horizontal edges")
            length(bottom) <= 1 && length(top) <= 1 ||
                error("Cross-section cell with two rows ending on one horizontal edge")
            chain = Int[]
            if !isempty(bottom)
                # Around the cell from the coarse bottom corner, leaving out the bottom edge.
                apex = node(bottom[1], l0)
                push!(chain, node(coarse, l0), node(coarse, l1))
                for i in reverse(fine_levels)
                    push!(chain, node(fine, i))
                end
            elseif !isempty(top)
                apex = node(top[1], l1)
                push!(chain, node(coarse, l1), node(coarse, l0))
                for i in fine_levels
                    push!(chain, node(fine, i))
                end
            elseif surface_distance[l1] >= surface_distance[l0]
                apex = node(coarse, l1)
                push!(chain, node(coarse, l0))
                for i in fine_levels
                    push!(chain, node(fine, i))
                end
            else
                apex = node(coarse, l0)
                push!(chain, node(coarse, l1))
                for i in reverse(fine_levels)
                    push!(chain, node(fine, i))
                end
            end
            fan!(apex, chain, (L, R))
        end
        push!(pair_triangles, ((L, R), length(triangles) - before))
        return nothing
    end
    # Every node of every active line exists even where no cell needs it (checked below).
    for row = 0:rows
        la, lb = line_range(row)
        for i in levels_in(row, la, lb)
            node(row, i)
        end
    end
    r1 = line_range(1)
    cells_between!(0, 1, r1[1], r1[2], Int[], Int[])
    for k = 2:rows
        rk, rp = line_range(k), line_range(k - 1)
        cells_between!(0, k, rp[2], rk[2], [k - 1], Int[])
        cells_between!(0, k, rk[1], rp[1], Int[], [k - 1])
    end
    for k = 1:(rows - 1)
        rk = line_range(k)
        cells_between!(k, k + 1, rk[1], rk[2], Int[], Int[])
    end
    area = sum(0.5 * orient2(t...) for t in triangles)
    expected = heights[rows] * (zs[end] - zs[1])
    abs(area - expected) <= 1.0e-9 * expected ||
        error("Cross-section area $area differs from $expected")
    used = falses(length(nodes))
    for t in triangles, n in t
        used[n] = true
    end
    all(used) || error("$(count(!, used)) unused cross-section nodes")
    return CrossSection(nodes, triangles, triangle_pair, pair_triangles)
end

# ---------------------------------------------------------------------------------------------
# Element tetrahedralisation: fans over face triangulations that depend on the face alone.

# A quad's two triangles: the diagonal through its lowest global index.
function quad_triangles(q1::Int32, q2::Int32, q3::Int32, q4::Int32)
    m = min(q1, q2, q3, q4)
    if m == q1 || m == q3
        return (q1, q2, q3), (q1, q3, q4)
    end
    return (q1, q2, q4), (q2, q3, q4)
end

# A prism (a1 a2 a3 with its normal toward b1 b2 b3): the fan from its lowest-ordered node
# over the outward-oriented faces not containing it — 3 tetrahedra, the production conforming
# split; every tetrahedron (apex, u, v, w) is positive when the fan is valid.
function push_prism_fan!(
    tetrahedra::Vector{NTuple{4, Int32}},
    a1::Int32,
    a2::Int32,
    a3::Int32,
    b1::Int32,
    b2::Int32,
    b3::Int32
)
    apex = min(a1, a2, a3, b1, b2, b3)
    faces = (
        (a1, a3, a2),
        (b1, b2, b3),
        quad_triangles(a1, a2, b2, b1)...,
        quad_triangles(a2, a3, b3, b2)...,
        quad_triangles(a3, a1, b1, b3)...
    )
    for (u, v, w) in faces
        (u == apex || v == apex || w == apex) && continue
        push!(tetrahedra, (apex, u, v, w))
    end
    return nothing
end

# A boundary triangle of a fan: skipped when degenerate or when it contains the apex.
function push_face_fan!(
    tetrahedra::Vector{NTuple{4, Int32}},
    apex::Int32,
    u::Int32,
    v::Int32,
    w::Int32
)
    (u == v || v == w || w == u) && return nothing
    (u == apex || v == apex || w == apex) && return nothing
    push!(tetrahedra, (apex, u, v, w))
    return nothing
end

# A swept quad whose cyclically adjacent corners may coincide (a collapsed or a fan edge): a
# quad takes the lowest-index diagonal, a triangle stays, less is no face.
function push_quad_fan!(
    tetrahedra::Vector{NTuple{4, Int32}},
    apex::Int32,
    q1::Int32,
    q2::Int32,
    q3::Int32,
    q4::Int32
)
    corners = (q1, q2, q3, q4)
    distinct = Int32[]
    for i = 1:4
        corners[i] == corners[mod1(i - 1, 4)] || push!(distinct, corners[i])
    end
    if length(distinct) == 4
        t1, t2 = quad_triangles(distinct[1], distinct[2], distinct[3], distinct[4])
        push_face_fan!(tetrahedra, apex, t1...)
        push_face_fan!(tetrahedra, apex, t2...)
    elseif length(distinct) == 3
        push_face_fan!(tetrahedra, apex, distinct[1], distinct[2], distinct[3])
    end
    return nothing
end

# A cross-section triangle swept between two columns (a1 a2 a3 on the first, with its normal
# toward b1 b2 b3 on the second) with coincidences allowed: the fan from the lowest-ordered
# node over the outward-oriented faces not containing it (a prism, a pyramid or a tetrahedron).
function push_swept_fan!(
    tetrahedra::Vector{NTuple{4, Int32}},
    a1::Int32,
    a2::Int32,
    a3::Int32,
    b1::Int32,
    b2::Int32,
    b3::Int32
)
    apex = min(a1, a2, a3, b1, b2, b3)
    push_face_fan!(tetrahedra, apex, a1, a3, a2)
    push_face_fan!(tetrahedra, apex, b1, b2, b3)
    push_quad_fan!(tetrahedra, apex, a1, a2, b2, b1)
    push_quad_fan!(tetrahedra, apex, a2, a3, b3, b2)
    push_quad_fan!(tetrahedra, apex, a3, a1, b1, b3)
    return nothing
end

# Vertical face between plan nodes P and Q with level lists `lp` / `lq` (the same first and
# last level, both listed from the substrate-surface side: ascending above it, descending
# below): the merge ladder from that side.
function ladder_triangles!(
    faces::Vector{NTuple{3, Int32}},
    global_index,
    p::Int,
    q::Int,
    lp::AbstractVector{Int},
    lq::AbstractVector{Int}
)
    ascending = lp[end] > lp[1]
    nearer(a, b) = ascending ? a < b : a > b
    i, j = 1, 1
    while i < length(lp) || j < length(lq)
        if i < length(lp) && j < length(lq) && lp[i + 1] == lq[j + 1]
            t1, t2 = quad_triangles(
                global_index(p, lp[i]),
                global_index(q, lq[j]),
                global_index(q, lq[j + 1]),
                global_index(p, lp[i + 1])
            )
            push!(faces, t1, t2)
            i += 1
            j += 1
        elseif j == length(lq) || (i < length(lp) && nearer(lp[i + 1], lq[j + 1]))
            push!(
                faces,
                (
                    global_index(p, lp[i]),
                    global_index(q, lq[j]),
                    global_index(p, lp[i + 1])
                )
            )
            i += 1
        else
            push!(
                faces,
                (
                    global_index(p, lp[i]),
                    global_index(q, lq[j]),
                    global_index(q, lq[j + 1])
                )
            )
            j += 1
        end
    end
    return nothing
end

# Levels of stack `s` between level indices la and lb inclusive (both are levels of `s`).
function stack_slice(levels::Vector{Int}, la::Int, lb::Int)
    first = searchsortedfirst(levels, la)
    last = searchsortedlast(levels, lb)
    (levels[first] == la && levels[last] == lb) ||
        error("Levels $la..$lb are not both levels of the stack")
    return view(levels, first:last)
end

# A right prism over the counter-clockwise plan triangle `t` between levels la and lb whose
# vertical edges carry the level lists `lists` (each containing la and lb): plain when no list
# has a level inside,
# else the fan from the lowest-ordered hanging-free corner on the substrate-surface side
# (`near_is_bottom`: la is the nearer level) over the merge-ladder faces. Returns (hanging,
# incidences).
function push_hanging_prism!(
    tetrahedra::Vector{NTuple{4, Int32}},
    faces::Vector{NTuple{3, Int32}},
    global_index,
    t::NTuple{3, Int},
    la::Int,
    lb::Int,
    lists::NTuple{3, Vector{Int}},
    near_is_bottom::Bool
)
    slices = (
        stack_slice(lists[1], la, lb),
        stack_slice(lists[2], la, lb),
        stack_slice(lists[3], la, lb)
    )
    # A node is hanging-free in THIS interval when its list has no level inside it; the apex
    # is the lowest-ordered hanging-free corner at the near level, so every adjacent face is a
    # fan from it.
    free = map(slice -> length(slice) == 2, slices)
    if all(free)
        push_prism_fan!(
            tetrahedra,
            global_index(t[1], la),
            global_index(t[2], la),
            global_index(t[3], la),
            global_index(t[1], lb),
            global_index(t[2], lb),
            global_index(t[3], lb)
        )
        return false, 0
    end
    any(free) || error(
        "Prism over plan nodes $t between levels $la and $lb has hanging nodes on all three " *
        "vertical edges (no corner carries the common levels only)"
    )
    # Faces oriented outward (`t` counter-clockwise in plan): the bottom reversed, the top as
    # is, the ladders walked from the near side (edge u -> v below the surface, v -> u above).
    empty!(faces)
    push!(
        faces,
        (global_index(t[1], la), global_index(t[3], la), global_index(t[2], la)),
        (global_index(t[1], lb), global_index(t[2], lb), global_index(t[3], lb))
    )
    near = near_is_bottom ? la : lb
    from_near(slice) = near_is_bottom ? slice : view(slice, length(slice):-1:1)
    for (u, v) in ((1, 2), (2, 3), (3, 1))
        p, q = near_is_bottom ? (u, v) : (v, u)
        ladder_triangles!(
            faces,
            global_index,
            t[p],
            t[q],
            from_near(slices[p]),
            from_near(slices[q])
        )
    end
    apex = minimum(global_index(t[u], near) for u = 1:3 if free[u])
    for (u, v, w) in faces
        (u == apex || v == apex || w == apex) && continue
        push!(tetrahedra, (apex, u, v, w))
    end
    return true, sum(length(slice) - 2 for slice in slices)
end

mutable struct GradedCounts
    band_swept_elements::Int
    band_strip_prisms::Int
    band_strip_prisms_hanging::Int
    band_collapse_prisms::Int
    band_tetrahedra::Int
    region_prisms_plain::Int
    region_prisms_hanging::Int
    region_tetrahedra::Int
    region_hanging_node_incidences::Int
end

# The cross-section variant of a set of planes (DESIGN.md 4.2 / 4.5): the row ranges around
# the fabricated steps of exactly those planes (one plane's band ends at the gap midpoint at
# the latest; a band whose edges belong to both planes — a coincident cross-plane run, a
# chain turning from one plane's edge onto the other's, a bump footprint — keeps every row
# active across the gap), its cross-section, the capped column tops' level lists (ids
# `list_offset + n` for n = 1..K-1) and the material of every cross-section triangle per
# partition class. The stacks are the same in every variant (only the ranges differ).
struct SweepVariant
    planes::Vector{Int}
    stacks::GradedStacks
    section::CrossSection
    top_levels::Vector{Vector{Int}}
    list_offset::Int
    section_material::Vector{Vector{Int8}}
    swept::Vector{Bool}
end

# The planes whose fabricated step a chain's edges belong to: a base edge is a metal edge of
# plane k when the conductor changes across it on that plane; a bump footprint edge belongs
# to both planes (the bump sidewall spans the gap). The union over the chain's base edges.
function chain_planes(chain::GradedChain, spec::PolygonSet, topology::PlanTopology)
    planes = Set{Int}()
    m = length(chain.plan_rows)
    for c = 1:(chain.closed ? m : m - 1)
        a, b = chain.plan_rows[c][1], chain.plan_rows[mod1(c + 1, m)][1]
        a == b && continue # a fan's shared base
        edge = a < b ? (a, b) : (b, a)
        for k in eachindex(spec.planes)
            haskey(topology.plane_edge_conductor[k], edge) && push!(planes, k)
        end
        haskey(topology.bump_edge_conductor, edge) && union!(planes, eachindex(spec.planes))
    end
    isempty(planes) &&
        error("Graded sweep: a band chain whose base edges are no metal edge")
    return sort!(collect(planes))
end

"""
    graded_sweep_elements(spec, plan, topology, stack, radial_um, radial_growth, radial_layers,
                          alpha, beta; region_ring, verbose) -> (tetrahedra, attributes, classes, record)

Build the graded cross-section sweep of a plan mesh (one or two planes, bumps) with its own
structured band: the band cells along the chains' columns (swept (0, k) cells, strip prisms,
collapse prisms) with the cross-section variant of the planes the chain's edges belong to,
and the region prisms with hanging nodes (`region_ring`: the region nodes adjacent to the
band carry Z_K). Node indices are (level - 1) x plan nodes + plan index (the tensor sweep's).
Every tetrahedron is checked positive and the volume per material against the plan's analytic
volume (`record["expected_volume_um3"]` is checked again by the caller on the written mesh).
"""
function graded_sweep_elements(
    spec::PolygonSet,
    plan::PlanMesh,
    topology::PlanTopology,
    stack::ZStack,
    radial_um::Float64,
    radial_growth::Float64,
    radial_layers::Int,
    alpha::Float64,
    beta::Float64;
    region_ring::Bool=true,
    verbose::Bool=true
)
    isempty(plan.chains) &&
        error("Graded sweep needs the own structured band (band_mode own)")
    zs = stack.levels
    n_plan = length(plan.xy)
    n_planes = length(spec.planes)
    # The fabricated step of every plane (trench bottom, substrate surface, metal top); the
    # interfaces every stack keeps: the steps' faces, the box ends, the backsides and, with
    # two planes, the gap midpoint (the gap's own far field, where one plane's band ends).
    step_faces = [
        begin
            s, f = plane.surface_z, plane.facing
            faces = (s - f * spec.overetch, s, s + f * spec.metal_thickness)
            (minimum(faces), maximum(faces))
        end for plane in spec.planes
    ]
    interfaces = Float64[]
    for plane in spec.planes
        s, f = plane.surface_z, plane.facing
        push!(interfaces, s - f * spec.overetch, s, s + f * spec.metal_thickness)
    end
    push!(interfaces, stack.z_bottom, stack.z_top, stack.backsides...)
    midpoint = n_planes == 2 ? stack.gap_midpoint : NaN
    n_planes == 2 && push!(interfaces, midpoint)
    surface_distance =
        [minimum(abs(z - plane.surface_z) for plane in spec.planes) for z in zs]
    heights = band_heights(radial_um, radial_growth, radial_layers)
    rows = radial_layers
    # The shared stacks (no range cap, the first plane's step: the levels do not depend on
    # the ranges) and the variants by plane set, built as the chains need them.
    shared = graded_stacks(
        zs,
        interfaces,
        heights,
        alpha,
        beta,
        step_faces[1],
        REGION_MESH_SIZE_MAX_UM
    )
    n_stacks = length(shared.levels)
    lists = copy(shared.levels)
    variants = Dict{Vector{Int}, SweepVariant}()
    variant_order = Vector{Int}[]
    function variant_of(planes::Vector{Int})
        return get!(variants, planes) do
            step = (
                minimum(step_faces[k][1] for k in planes),
                maximum(step_faces[k][2] for k in planes)
            )
            # One plane's band of a two-plane set ends at the gap midpoint on the gap side.
            cap = (-Inf, Inf)
            if n_planes == 2 && length(planes) == 1
                toward_gap = spec.planes[planes[1]].facing
                cap = toward_gap == 1 ? (-Inf, midpoint) : (midpoint, Inf)
            end
            stacks = graded_stacks(
                zs,
                interfaces,
                heights,
                alpha,
                beta,
                step,
                REGION_MESH_SIZE_MAX_UM;
                range_cap=cap
            )
            stacks.levels == shared.levels ||
                error("Graded sweep: the stacks differ between cross-section variants")
            section = cross_section(stacks, zs, heights, surface_distance)
            top_levels = column_top_levels(stacks)
            list_offset = length(lists)
            append!(lists, top_levels[1:(rows - 1)])
            section_material = [
                begin
                    zc = sum(zs[section.nodes[n][2]] for n in t) / 3
                    [material(spec, class, zc) for class in plan.classes]
                end for t in section.triangles
            ]
            swept = [pair[1] == 0 for pair in section.triangle_pair]
            push!(variant_order, planes)
            return SweepVariant(
                planes,
                stacks,
                section,
                top_levels,
                list_offset,
                section_material,
                swept
            )
        end
    end
    top_list_id(variant::SweepVariant, n::Int) =
        n == rows ? rows + 1 : variant.list_offset + n
    chain_variant =
        [variant_of(chain_planes(chain, spec, topology)) for chain in plan.chains]

    # Plan nodes of the band: row (-1 for a region node), the most rows of a column through
    # the node, and, for column tops, the list id; the capped columns and the fans.
    node_row = fill(-1, n_plan)
    node_column_rows = zeros(Int, n_plan)
    node_list = zeros(Int, n_plan)
    columns, capped_columns, fan_sectors = 0, 0, 0
    for (chain_index, chain) in enumerate(plan.chains)
        variant = chain_variant[chain_index]
        m = length(chain.plan_rows)
        for (c, plan_rows) in enumerate(chain.plan_rows)
            n = length(plan_rows) - 1
            1 <= n <= rows ||
                error("Graded sweep: a column with $n rows (the band has $rows)")
            columns += 1
            n < rows && (capped_columns += 1)
            for (k, p) in enumerate(plan_rows)
                row = k - 1
                node_row[p] in (-1, row) || error(
                    "Plan node $p at $(plan.xy[p]) is row $(node_row[p]) of one column and " *
                    "row $row of another"
                )
                node_row[p] = row
                node_column_rows[p] = max(node_column_rows[p], n)
            end
            top = plan_rows[end]
            node_list[top] in (0, top_list_id(variant, n)) || error(
                "Plan node $top at $(plan.xy[top]) is the top of columns with different " *
                "rows or cross-section variants"
            )
            node_list[top] = top_list_id(variant, n)
            d = mod1(c + 1, m)
            (chain.closed || c < m) &&
                chain.plan_rows[c][1] == chain.plan_rows[d][1] &&
                (fan_sectors += 1)
        end
    end
    node_size = fill(Inf, n_plan)
    for ((a, b), _) in topology.edge_incidence
        d = point_distance(plan.xy[a], plan.xy[b])
        node_size[a] = min(node_size[a], d)
        node_size[b] = min(node_size[b], d)
    end
    width_last = heights[rows] - (rows == 1 ? 0.0 : heights[rows - 1])
    ladder = n_stacks - 1 - rows
    # The region ring (decision 275; "ladder0" in the M2 recovery diagnostic): a region node
    # adjacent to the band (sharing a plan edge with a band node, i.e. a column top) carries
    # the band's outermost stack Z_K, so the ring of region triangles touching the band keeps
    # the band's z resolution; the plan-size ladder applies beyond.
    ring_node = falses(n_plan)
    if region_ring
        for (a, b) in keys(topology.edge_incidence)
            (node_row[a] == -1) == (node_row[b] == -1) && continue
            ring_node[node_row[a] == -1 ? a : b] = true
        end
    end
    region_nodes, ring_nodes = 0, 0
    for p = 1:n_plan
        node_row[p] == -1 || continue
        region_nodes += 1
        size = min(node_size[p], REGION_MESH_SIZE_MAX_UM)
        j =
            size > width_last ?
            floor(Int, log2(size / width_last) + LADDER_LOG2_TOLERANCE) : 0
        if ring_node[p]
            j = 0
            ring_nodes += 1
        end
        node_list[p] = rows + clamp(j, 0, ladder) + 1
    end
    nodes_by_list = [count(==(id), node_list) for id = 1:length(lists)]

    # The node order of the face rules (see the header): levels by distance to the substrate
    # surface, plan nodes by band row descending, taller columns first, then the welding
    # order; region nodes last. The tetrahedra are remapped to the tensor sweep's numbering
    # (level - 1) x plan nodes + plan index before they are returned.
    level_by_rank = sortperm(eachindex(zs); by=i -> (surface_distance[i], -zs[i]))
    rank_of_level = invperm(level_by_rank)
    plan_by_order = sortperm(
        1:n_plan;
        by=p -> (node_row[p] == -1, -node_row[p], -node_column_rows[p], p)
    )
    order_of_plan = invperm(plan_by_order)
    global_index(p::Int, level::Int) =
        Int32((rank_of_level[level] - 1) * n_plan + order_of_plan[p])
    plan_of(index::Int32) = plan_by_order[mod(Int(index) - 1, n_plan) + 1]
    level_of(index::Int32) = level_by_rank[(Int(index) - 1) ÷ n_plan + 1]
    tail_index(index::Int32) = Int32((level_of(index) - 1) * n_plan + plan_of(index))
    point(index::Int32) = begin
        p, level = plan_of(index), level_of(index)
        (plan.xy[p][1], plan.xy[p][2], zs[level])
    end
    near_is_bottom(la::Int, lb::Int) = surface_distance[la] < surface_distance[lb]
    centroid(indices) = begin
        points = map(point, indices)
        ntuple(i -> sum(p[i] for p in points) / length(points), 3)
    end
    ccw(t::NTuple{3, Int}) =
        orient(plan.xy[t[1]], plan.xy[t[2]], plan.xy[t[3]]) > 0.0 ? t : (t[1], t[3], t[2])
    tetrahedra = NTuple{4, Int32}[]
    tetrahedron_attribute = Int8[]
    tetrahedron_class = Int16[]
    counts = GradedCounts(0, 0, 0, 0, 0, 0, 0, 0, 0)
    # Every tetrahedron of a cell is checked positive where it is built (the 3D backstop: the
    # faces are oriented outward, so a non-positive fan tetrahedron means the apex does not
    # see that face from inside — the cell is not star-shaped from it).
    function finish_cell!(before, attribute, class, describe)
        for index = (before + 1):length(tetrahedra)
            t = tetrahedra[index]
            volume = signed_volume(point(t[1]), point(t[2]), point(t[3]), point(t[4]))
            volume > 0.0 || error(
                "Graded sweep: non-positive tetrahedron (volume $volume) in " *
                describe() *
                " (nodes $(point(t[1])) $(point(t[2])) $(point(t[3])) $(point(t[4])))"
            )
            push!(tetrahedron_attribute, attribute)
            push!(tetrahedron_class, class)
        end
        return length(tetrahedra) - before
    end

    # Band cells per column pair, with the chain's cross-section variant.
    faces = NTuple{3, Int32}[]
    for (chain_index, chain) in enumerate(plan.chains)
        variant = chain_variant[chain_index]
        stacks, section, top_levels = variant.stacks, variant.section, variant.top_levels
        section_material, swept = variant.section_material, variant.swept
        m = length(chain.plan_rows)
        class = chain.class
        for c = 1:(chain.closed ? m : m - 1)
            d = mod1(c + 1, m)
            P, Q = chain.plan_rows[c], chain.plan_rows[d]
            n_P, n_Q = length(P) - 1, length(Q) - 1
            describe_pair() = "chain $chain_index columns $c / $d of class $class"
            # (0, k) cells swept with the rows clamped to each column's row count (at a row
            # end next to a capped column the corner triangle (k, l1) (k-1, l0) (k, l0) sweeps
            # a vertical segment onto a triangle: a pyramid from the row-end node, see the
            # header).
            for (ti, t) in enumerate(section.triangles)
                swept[ti] || continue
                attribute = section_material[ti][class]
                attribute == 0 && continue
                before = length(tetrahedra)
                (r1, l1), (r2, l2), (r3, l3) =
                    section.nodes[t[1]], section.nodes[t[2]], section.nodes[t[3]]
                a = (
                    global_index(P[min(r1, n_P) + 1], l1),
                    global_index(P[min(r2, n_P) + 1], l2),
                    global_index(P[min(r3, n_P) + 1], l3)
                )
                b = (
                    global_index(Q[min(r1, n_Q) + 1], l1),
                    global_index(Q[min(r2, n_Q) + 1], l2),
                    global_index(Q[min(r3, n_Q) + 1], l3)
                )
                # The first triangle's normal must point toward the second (outward faces);
                # a collapsed first triangle is judged from the second.
                toward = signed_volume(point(a[1]), point(a[2]), point(a[3]), centroid(b))
                toward == 0.0 && (
                    toward =
                        -signed_volume(point(b[1]), point(b[2]), point(b[3]), centroid(a))
                )
                toward != 0.0 || continue
                if toward < 0.0
                    a, b = (a[1], a[3], a[2]), (b[1], b[3], b[2])
                end
                push_swept_fan!(tetrahedra, a..., b...)
                length(tetrahedra) > before || continue
                counts.band_swept_elements += 1
                counts.band_tetrahedra += finish_cell!(
                    before,
                    attribute,
                    Int16(class),
                    () ->
                        "swept cell (rows, z) $(section.nodes[t[1]]) " *
                        "$(section.nodes[t[2]]) $(section.nodes[t[3]]) of " *
                        describe_pair()
                )
            end
            # (k-1, k) strips within row k-1's range as prisms over the plan band cells.
            function strip_prisms!(triangle, lists_of, k, collapse)
                if ccw(triangle) != triangle
                    triangle = (triangle[1], triangle[3], triangle[2])
                    lists_of = (lists_of[1], lists_of[3], lists_of[2])
                end
                lo, hi = stacks.ranges[k - 1]
                levels = stack_slice(stacks.levels[k + 1], lo, hi)
                for interval = 1:(length(levels) - 1)
                    la, lb = levels[interval], levels[interval + 1]
                    attribute = material(spec, plan.classes[class], 0.5 * (zs[la] + zs[lb]))
                    attribute == 0 && continue
                    before = length(tetrahedra)
                    hanging, _ = push_hanging_prism!(
                        tetrahedra,
                        faces,
                        global_index,
                        triangle,
                        la,
                        lb,
                        lists_of,
                        near_is_bottom(la, lb)
                    )
                    if collapse
                        counts.band_collapse_prisms += 1
                    else
                        counts.band_strip_prisms += 1
                        hanging && (counts.band_strip_prisms_hanging += 1)
                    end
                    counts.band_tetrahedra += finish_cell!(
                        before,
                        attribute,
                        Int16(class),
                        () ->
                            "$(collapse ? "collapse" : "strip") prism over plan nodes " *
                            "$triangle between z $(zs[la]) and $(zs[lb]) (rows $(k - 1), " *
                            "$k) of " *
                            describe_pair()
                    )
                end
                return nothing
            end
            n_min, n_max = minmax(n_P, n_Q)
            for k = 2:n_min
                # The quad P_{k-1} P_k Q_k Q_{k-1} split by the diagonal through its lowest-
                # ordered node (the rule of its horizontal faces).
                corners = (P[k], P[k + 1], Q[k + 1], Q[k])
                t1, t2 = quad_triangles(
                    Int32(order_of_plan[corners[1]]),
                    Int32(order_of_plan[corners[2]]),
                    Int32(order_of_plan[corners[3]]),
                    Int32(order_of_plan[corners[4]])
                )
                list_of = Dict(
                    P[k] => stacks.levels[k],
                    P[k + 1] => stacks.levels[k + 1],
                    Q[k] => stacks.levels[k],
                    Q[k + 1] => stacks.levels[k + 1]
                )
                for t in (t1, t2)
                    triangle =
                        (plan_by_order[t[1]], plan_by_order[t[2]], plan_by_order[t[3]])
                    strip_prisms!(
                        triangle,
                        (list_of[triangle[1]], list_of[triangle[2]], list_of[triangle[3]]),
                        k,
                        false
                    )
                end
            end
            # Collapse triangles of the taller column onto the shorter column's top.
            tall, short = n_P >= n_Q ? (P, Q) : (Q, P)
            for k = (n_min + 1):n_max
                strip_prisms!(
                    (tall[k], tall[k + 1], short[end]),
                    (stacks.levels[k], stacks.levels[k + 1], top_levels[n_min]),
                    k,
                    true
                )
            end
        end
    end

    # Region prisms: plan triangle x interval of the common levels of its corners' lists (one
    # plane: the lists are nested and the coarsest is one of them; two planes: the tops of
    # chains of different variants are fine near their own plane and Z_K near the other, so
    # the common list is their intersection and, interval by interval, one corner is still
    # hanging-free), hanging nodes of the finer lists on their vertical edges.
    for (index, t) in enumerate(plan.triangles)
        plan.band_triangle[index] && continue
        class = plan.triangle_class[index]
        for p in t
            node_list[p] > 0 || error(
                "Region triangle corner $p at $(plan.xy[p]) is a band node that is not a " *
                "column top"
            )
        end
        t = ccw(t)
        node_lists =
            (lists[node_list[t[1]]], lists[node_list[t[2]]], lists[node_list[t[3]]])
        levels = sort!(intersect(node_lists...))
        (levels[1] == 1 && levels[end] == length(zs)) ||
            error("Region triangle $t: the corners' common levels miss a box face")
        for interval = 1:(length(levels) - 1)
            la, lb = levels[interval], levels[interval + 1]
            attribute = material(spec, plan.classes[class], 0.5 * (zs[la] + zs[lb]))
            attribute == 0 && continue
            before = length(tetrahedra)
            hanging, incidences = push_hanging_prism!(
                tetrahedra,
                faces,
                global_index,
                t,
                la,
                lb,
                node_lists,
                near_is_bottom(la, lb)
            )
            if hanging
                counts.region_prisms_hanging += 1
                counts.region_hanging_node_incidences += incidences
            else
                counts.region_prisms_plain += 1
            end
            counts.region_tetrahedra += finish_cell!(
                before,
                attribute,
                Int16(class),
                () ->
                    "region prism over plan nodes $t between z $(zs[la]) and $(zs[lb])"
            )
        end
    end

    # Analytic volume per material of the plan (every plan triangle through every z interval,
    # the tensor sweep's), checked against the swept volume here and on the written mesh.
    material_table = [
        material(spec, class, 0.5 * (zs[i] + zs[i + 1])) for
        class in plan.classes, i = 1:(length(zs) - 1)
    ]
    expected_volume = Dict{String, Float64}()
    for (index, t) in enumerate(plan.triangles)
        area = triangle_area(plan.xy[t[1]], plan.xy[t[2]], plan.xy[t[3]])
        class = plan.triangle_class[index]
        for i = 1:(length(zs) - 1)
            attribute = material_table[class, i]
            attribute == 0 && continue
            key = string(attribute)
            expected_volume[key] =
                get(expected_volume, key, 0.0) + area * (zs[i + 1] - zs[i])
        end
    end
    swept_volume = Dict{String, Float64}()
    for (index, t) in enumerate(tetrahedra)
        key = string(tetrahedron_attribute[index])
        swept_volume[key] =
            get(swept_volume, key, 0.0) +
            signed_volume(point(t[1]), point(t[2]), point(t[3]), point(t[4]))
    end
    for (key, expected) in expected_volume
        actual = get(swept_volume, key, 0.0)
        abs(actual - expected) <= VOLUME_CLOSURE_TOLERANCE * expected || error(
            "Graded sweep: the volume of material $key, $actual um^3, differs from the " *
            "plan's $expected um^3 (an overlap or a gap)"
        )
    end
    keys(swept_volume) == keys(expected_volume) ||
        error("Graded sweep: materials $(keys(swept_volume)) vs $(keys(expected_volume))")
    for index in eachindex(tetrahedra)
        t = tetrahedra[index]
        tetrahedra[index] =
            (tail_index(t[1]), tail_index(t[2]), tail_index(t[3]), tail_index(t[4]))
    end

    # The chains, columns, capped columns and fan sectors per variant.
    variant_counts = Dict(planes => zeros(Int, 4) for planes in variant_order)
    for (chain_index, chain) in enumerate(plan.chains)
        v = variant_counts[chain_variant[chain_index].planes]
        v[1] += 1
        m = length(chain.plan_rows)
        for (c, plan_rows) in enumerate(chain.plan_rows)
            v[2] += 1
            length(plan_rows) - 1 < rows && (v[3] += 1)
            d = mod1(c + 1, m)
            (chain.closed || c < m) &&
                chain.plan_rows[c][1] == chain.plan_rows[d][1] &&
                (v[3 + 1] += 1)
        end
    end
    cross_sections = [
        begin
            variant = variants[planes]
            v = variant_counts[planes]
            Dict{String, Any}(
                "planes" => [spec.planes[k].name for k in planes],
                "step_faces_z_um" => [
                    minimum(step_faces[k][1] for k in planes),
                    maximum(step_faces[k][2] for k in planes)
                ],
                "row_ranges_z_um" =>
                    [[zs[lo], zs[hi]] for (lo, hi) in variant.stacks.ranges],
                "column_top_list_sizes" => [length(l) for l in variant.top_levels],
                "plan_nodes_by_capped_top_list" =>
                    nodes_by_list[(variant.list_offset + 1):(variant.list_offset + rows - 1)],
                "chains" => v[1],
                "columns" => v[2],
                "capped_columns" => v[3],
                "fan_sectors" => v[4],
                "cross_section_nodes" => length(variant.section.nodes),
                "cross_section_triangles" => length(variant.section.triangles),
                "cross_section_triangles_by_line_pair" => [
                    Dict("lines" => collect(pair), "triangles" => n) for
                    (pair, n) in variant.section.pair_triangles
                ]
            )
        end for planes in variant_order
    ]
    record = Dict{String, Any}(
        "alpha" => alpha,
        "beta" => beta,
        "rows" => rows,
        "heights_um" => heights,
        "stack_spacings_um" => shared.spacings,
        "stack_sizes" => [length(l) for l in shared.levels],
        "stack_levels_z_um" => [[zs[i] for i in l] for l in shared.levels],
        "gap_midpoint_z_um" => n_planes == 2 ? midpoint : nothing,
        "cross_sections" => cross_sections,
        "region_ladder_stacks" => ladder,
        "region_size_max_um" => REGION_MESH_SIZE_MAX_UM,
        "region_ring" => region_ring,
        "region_ring_plan_nodes" => ring_nodes,
        "plan_nodes_by_stack" => nodes_by_list[1:n_stacks],
        "region_plan_nodes" => region_nodes,
        "chains" => length(plan.chains),
        "columns" => columns,
        "capped_columns" => capped_columns,
        "fan_sectors" => fan_sectors,
        "band_swept_elements" => counts.band_swept_elements,
        "band_strip_prisms" => counts.band_strip_prisms,
        "band_strip_prisms_hanging" => counts.band_strip_prisms_hanging,
        "band_collapse_prisms" => counts.band_collapse_prisms,
        "band_tetrahedra" => counts.band_tetrahedra,
        "region_prisms_plain" => counts.region_prisms_plain,
        "region_prisms_hanging" => counts.region_prisms_hanging,
        "region_hanging_node_incidences" => counts.region_hanging_node_incidences,
        "region_tetrahedra" => counts.region_tetrahedra,
        "expected_volume_um3" => expected_volume
    )
    verbose && println(
        "Graded sweep: alpha ",
        alpha,
        ", beta ",
        beta,
        ", stacks ",
        record["stack_sizes"],
        ", cross-sections ",
        [
            (section["planes"], section["row_ranges_z_um"], section["chains"]) for
            section in cross_sections
        ],
        ", columns ",
        columns,
        " (capped ",
        capped_columns,
        ", fan sectors ",
        fan_sectors,
        "), swept cells ",
        counts.band_swept_elements,
        ", strip prisms ",
        counts.band_strip_prisms,
        " (hanging ",
        counts.band_strip_prisms_hanging,
        "), collapse prisms ",
        counts.band_collapse_prisms,
        ", region prisms ",
        counts.region_prisms_plain,
        " plain + ",
        counts.region_prisms_hanging,
        " hanging (ring nodes ",
        ring_nodes,
        "); volume closes to ",
        maximum(abs(swept_volume[k] - v) / v for (k, v) in expected_volume; init=0.0)
    )
    return tetrahedra, tetrahedron_attribute, tetrahedron_class, record
end
