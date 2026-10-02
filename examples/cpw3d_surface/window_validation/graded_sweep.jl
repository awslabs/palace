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
#     (w_k = h_k - h_{k-1} the row's width; a level is kept iff it is a material interface or at
#     least alpha w_k from the last kept level and from the next interface), nested by
#     construction; the region (Gmsh interior) nodes carry the ladder Z_{K+j} = thin(Z_{K+j-1},
#     alpha w_K 2^j) chosen by their plan size (the shortest incident plan edge, at most the
#     region mesh size);
#   * row ranges: row line k is active for z in [lo_k, hi_k] = the first level of Z_{k+1} beyond
#     the metal faces -/+ beta h_k, strictly nested (inner rows end first); row K spans the box;
#   * the cross-section per column: between two neighbouring active lines the cells are the
#     consecutive levels of the COARSER stack; the finer line's extra levels hang on one vertical
#     edge and the end node of the row that ended at the cell's bottom (started at its top)
#     hangs on that horizontal edge; a fan from the coarse line's node at the hanging-free
#     horizontal edge triangulates every cell without a degenerate triangle;
#   * the sweep: cross-section node (k, z) is attached to the row-k plan node of every column
#     (exactly the plan band's geometry), consecutive columns are joined into prisms;
#   * the region: plan triangle x interval of the coarsest stack of its three nodes; the finer
#     nodes' extra levels hang on their vertical edges;
#   * tetrahedra: every element is a fan from one apex over its boundary triangulation, and
#     every face is triangulated by a rule of the face alone, so neighbours always agree: a quad
#     takes the diagonal through its lowest global node index; a vertical face with hanging
#     nodes is a merge ladder from the bottom (a quad where both sides step together, a triangle
#     advancing the lower side otherwise); the apex is the lowest global index among the
#     element's corners on hanging-free vertical edges (the fan over a quad through that
#     corner's diagonal is exact, the ladder from a hanging-free bottom corner is a fan from it).
#     Global node index = (level - 1) x plan nodes + plan index, as in the tensor sweep, so a
#     plain prism gets the production three-tetrahedron split and the surface tagging is shared.
#
# M1 scope (decision 263): one plane, no bumps, every column of a chain with the same row count
# and no fans; capped columns / collapse, fans and two planes are refused with a message
# (milestones M2 / M3). alpha and beta are dimensionless (defaults 1 and 3): alpha = 1 keeps the
# r x r edge cell of the recorded family; the z-coarsened outer rows are a deliberate deviation
# to be shown harmless by the acceptance solves.

const DEFAULT_GRADED_ALPHA = 1.0
const DEFAULT_GRADED_BETA = 3.0
const LEVEL_TOLERANCE_UM = 1.0e-9
# A row range end is nudged past the previous row's end by this before the level search (strict
# nesting), compared with a tolerance three decades below it.
const RANGE_NUDGE_UM = 1.0e-9
const RANGE_COMPARE_UM = 1.0e-12

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
        next = findnext(j -> is_interface(zs[levels[j]], interfaces), eachindex(levels), position)
        room = next === nothing ? Inf : zs[levels[next]] - z
        if z - zs[kept[end]] >= spacing - LEVEL_TOLERANCE_UM &&
           room >= spacing - LEVEL_TOLERANCE_UM
            push!(kept, index)
        end
    end
    return kept
end

"""
    graded_stacks(zs, interfaces, heights, alpha, beta, metal_faces, region_size_max) -> GradedStacks

The nested stacks Z_0..Z_K of the band rows (heights h_1..h_K), the region ladder beyond
Z_K up to the spacing alpha w_K 2^J <= alpha region_size_max, and the active range of every
row (the first level of Z_{k+1} at or beyond the metal faces -/+ beta h_k; row K spans the box).
"""
function graded_stacks(
    zs::Vector{Float64},
    interfaces::Vector{Float64},
    heights::Vector{Float64},
    alpha::Float64,
    beta::Float64,
    metal_faces::NTuple{2, Float64},
    region_size_max::Float64
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
    ladder = 0
    while widths[rows] * 2.0^(ladder + 1) <= region_size_max
        ladder += 1
        push!(spacings, alpha * widths[rows] * 2.0^ladder)
        push!(levels, thin_levels(zs, levels[end], spacings[end], interfaces))
    end
    ranges = Vector{NTuple{2, Int}}(undef, rows)
    for k = 1:rows
        if k == rows
            ranges[k] = (1, length(zs))
            continue
        end
        coarser = levels[k + 2]
        limit_lo = metal_faces[1] - beta * heights[k]
        limit_hi = metal_faces[2] + beta * heights[k]
        if k > 1
            limit_lo = min(limit_lo, zs[ranges[k - 1][1]] - RANGE_NUDGE_UM)
            limit_hi = max(limit_hi, zs[ranges[k - 1][2]] + RANGE_NUDGE_UM)
        end
        lo = findlast(i -> zs[i] <= limit_lo + RANGE_COMPARE_UM, coarser)
        hi = findfirst(i -> zs[i] >= limit_hi - RANGE_COMPARE_UM, coarser)
        lo = coarser[lo === nothing ? 1 : lo]
        hi = coarser[hi === nothing ? length(coarser) : hi]
        ranges[k] = (lo, hi)
    end
    return GradedStacks(levels, spacings, rows, ranges)
end

# ---------------------------------------------------------------------------------------------
# The (n, z) cross-section of a column: nodes (row, level index), fan triangles.

struct CrossSection
    nodes::Vector{NTuple{2, Int}} # (row 0..K, level index)
    triangles::Vector{NTuple{3, Int}}
    pair_triangles::Vector{Tuple{NTuple{2, Int}, Int}} # ((line, line), triangles) per pair
end

function cross_section(stacks::GradedStacks, zs::Vector{Float64}, heights::Vector{Float64})
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
    pair_triangles = Tuple{NTuple{2, Int}, Int}[]
    orient2(a, b, c) = begin
        (ra, la), (rb, lb), (rc, lc) = nodes[a], nodes[b], nodes[c]
        (position(rb) - position(ra)) * (zs[lc] - zs[la]) -
        (zs[lb] - zs[la]) * (position(rc) - position(ra))
    end
    function fan!(apex, chain)
        for i = 1:(length(chain) - 1)
            a, b, c = apex, chain[i], chain[i + 1]
            orient2(a, b, c) > 0.0 || ((b, c) = (c, b))
            orient2(a, b, c) > 0.0 || error(
                "Degenerate cross-section triangle (row, z) $(nodes[a]) $(nodes[b]) $(nodes[c])"
            )
            push!(triangles, (a, b, c))
        end
    end
    # Cells between lines L and R (rows) over the level range [la, lb]; `bottom_hanging` /
    # `top_hanging`: rows whose end node hangs on the first cell's bottom / last cell's top.
    function cells_between!(L, R, la, lb, bottom_hanging, top_hanging)
        lb > la || return
        before = length(triangles)
        levels_L, levels_R = levels_in(L, la, lb), levels_in(R, la, lb)
        coarse, fine, cl, fl =
            length(levels_L) <= length(levels_R) ? (L, R, levels_L, levels_R) :
            (R, L, levels_R, levels_L)
        issubset(cl, fl) || error("Stacks of rows $L and $R are not nested in levels $la..$lb")
        (cl[1] == la && cl[end] == lb) ||
            error("Range ends $la..$lb are not levels of row $coarse")
        left = min(L, R)
        for ci = 1:(length(cl) - 1)
            l0, l1 = cl[ci], cl[ci + 1]
            fine_levels = [i for i in fl if l0 <= i <= l1]
            bottom = ci == 1 ? bottom_hanging : Int[]
            top = ci == length(cl) - 1 ? top_hanging : Int[]
            (isempty(bottom) || isempty(top)) ||
                error("Cross-section cell with hanging nodes on both horizontal edges")
            chain = Int[]
            if isempty(top)
                apex = node(coarse, l1)
                if coarse == left
                    push!(chain, node(coarse, l0))
                    for row in sort(bottom)
                        push!(chain, node(row, l0))
                    end
                    for i in fine_levels
                        push!(chain, node(fine, i))
                    end
                else
                    for i in reverse(fine_levels)
                        push!(chain, node(fine, i))
                    end
                    for row in sort(bottom; rev=true)
                        push!(chain, node(row, l0))
                    end
                    push!(chain, node(coarse, l0))
                end
            else
                apex = node(coarse, l0)
                if coarse == left
                    for i in fine_levels
                        push!(chain, node(fine, i))
                    end
                    for row in sort(top; rev=true)
                        push!(chain, node(row, l1))
                    end
                    push!(chain, node(coarse, l1))
                else
                    push!(chain, node(coarse, l1))
                    for row in sort(top)
                        push!(chain, node(row, l1))
                    end
                    for i in reverse(fine_levels)
                        push!(chain, node(fine, i))
                    end
                end
            end
            fan!(apex, chain)
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
    return CrossSection(nodes, triangles, pair_triangles)
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

# A prism (a1 a2 a3 below b1 b2 b3): the fan from its lowest global index over the faces not
# containing it — 3 tetrahedra, the production conforming split.
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
        (a1, a2, a3),
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

# Vertical face between plan nodes P and Q with level lists `lp` / `lq` (the same first and
# last level): the merge ladder from the bottom.
function ladder_triangles!(
    faces::Vector{NTuple{3, Int32}},
    global_index,
    p::Int,
    q::Int,
    lp::AbstractVector{Int},
    lq::AbstractVector{Int}
)
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
        elseif j == length(lq) || (i < length(lp) && lp[i + 1] < lq[j + 1])
            push!(
                faces,
                (global_index(p, lp[i]), global_index(q, lq[j]), global_index(p, lp[i + 1]))
            )
            i += 1
        else
            push!(
                faces,
                (global_index(p, lp[i]), global_index(q, lq[j]), global_index(q, lq[j + 1]))
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

mutable struct GradedCounts
    band_prisms::Int
    band_tetrahedra::Int
    region_prisms_plain::Int
    region_prisms_hanging::Int
    region_tetrahedra::Int
    region_hanging_nodes::Int
end

"""
    graded_sweep_elements(spec, plan, topology, stack, radial_um, radial_growth, radial_layers,
                          alpha, beta; verbose) -> (tetrahedra, attributes, classes, record)

Build the graded cross-section sweep of a one-plane plan mesh with its own structured band:
the band prisms along the chains' columns and the region prisms with hanging nodes. Node
indices are (level - 1) x plan nodes + plan index (the tensor sweep's).
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
    verbose::Bool=true
)
    length(spec.planes) == 1 ||
        error("Graded sweep: two planes are not supported yet (milestone M3)")
    isempty(spec.bumps) || error("Graded sweep: bumps are not supported yet (milestone M3)")
    isempty(plan.chains) &&
        error("Graded sweep needs the own structured band (band_mode own)")
    zs = stack.levels
    n_plan = length(plan.xy)
    plane = spec.planes[1]
    s, f = plane.surface_z, plane.facing
    metal_top = s + f * spec.metal_thickness
    trench = s - f * spec.overetch
    interfaces = [trench, s, metal_top, stack.z_bottom, stack.z_top, stack.backsides...]
    heights = band_heights(radial_um, radial_growth, radial_layers)
    stacks = graded_stacks(
        zs,
        interfaces,
        heights,
        alpha,
        beta,
        (min(s, metal_top), max(s, metal_top)),
        REGION_MESH_SIZE_MAX_UM
    )
    section = cross_section(stacks, zs, heights)
    rows = stacks.rows
    n_stacks = length(stacks.levels)

    # Stack of every plan node: band nodes by row (consistent across chains), region nodes
    # by plan size (the shortest incident edge) on the ladder.
    node_stack = fill(-1, n_plan)
    columns = 0
    for chain in plan.chains
        for plan_rows in chain.plan_rows
            length(plan_rows) == rows + 1 || error(
                "Graded sweep: a column with $(length(plan_rows) - 1) of $rows rows " *
                "(a capped column) is not supported yet (milestone M2)"
            )
            columns += 1
            for (k, p) in enumerate(plan_rows)
                row = k - 1
                node_stack[p] in (-1, row) || error(
                    "Plan node $p at $(plan.xy[p]) is row $(node_stack[p]) of one chain and " *
                    "row $row of another"
                )
                node_stack[p] = row
            end
        end
        for c = 1:(length(chain.plan_rows) - (chain.closed ? 0 : 1))
            d = mod1(c + 1, length(chain.plan_rows))
            chain.plan_rows[c][1] != chain.plan_rows[d][1] ||
                error("Graded sweep: fan columns are not supported yet (milestone M2)")
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
    region_nodes = 0
    for p = 1:n_plan
        node_stack[p] == -1 || continue
        region_nodes += 1
        size = min(node_size[p], REGION_MESH_SIZE_MAX_UM)
        j = size > width_last ? floor(Int, log2(size / width_last) + LEVEL_TOLERANCE_UM) : 0
        node_stack[p] = rows + clamp(j, 0, ladder)
    end
    nodes_by_stack = [count(==(s - 1), node_stack) for s = 1:n_stacks]

    global_index(p::Int, level::Int) = Int32((level - 1) * n_plan + p)
    tetrahedra = NTuple{4, Int32}[]
    tetrahedron_attribute = Int8[]
    tetrahedron_class = Int32[]
    counts = GradedCounts(0, 0, 0, 0, 0, 0)
    function push_element_tetrahedra!(before, attribute, class)
        for _ = (before + 1):length(tetrahedra)
            push!(tetrahedron_attribute, attribute)
            push!(tetrahedron_class, class)
        end
        return length(tetrahedra) - before
    end

    # Band prisms: cross-section triangle x consecutive columns.
    section_material = [
        begin
            zc = sum(zs[section.nodes[n][2]] for n in t) / 3
            [material(spec, class, zc) for class in plan.classes]
        end for t in section.triangles
    ]
    for chain in plan.chains
        m = length(chain.plan_rows)
        for c = 1:(chain.closed ? m : m - 1)
            lower, upper = chain.plan_rows[c], chain.plan_rows[mod1(c + 1, m)]
            for (ti, t) in enumerate(section.triangles)
                attribute = section_material[ti][chain.class]
                attribute == 0 && continue
                counts.band_prisms += 1
                before = length(tetrahedra)
                (r1, l1), (r2, l2), (r3, l3) =
                    section.nodes[t[1]], section.nodes[t[2]], section.nodes[t[3]]
                push_prism_fan!(
                    tetrahedra,
                    global_index(lower[r1 + 1], l1),
                    global_index(lower[r2 + 1], l2),
                    global_index(lower[r3 + 1], l3),
                    global_index(upper[r1 + 1], l1),
                    global_index(upper[r2 + 1], l2),
                    global_index(upper[r3 + 1], l3)
                )
                counts.band_tetrahedra +=
                    push_element_tetrahedra!(before, attribute, chain.class)
            end
        end
    end

    # Region prisms: plan triangle x interval of the coarsest of its nodes' stacks, hanging
    # nodes of the finer stacks on their vertical edges.
    faces = NTuple{3, Int32}[]
    for (index, t) in enumerate(plan.triangles)
        plan.band_triangle[index] && continue
        class = plan.triangle_class[index]
        stack_ids = (node_stack[t[1]], node_stack[t[2]], node_stack[t[3]])
        element_stack = maximum(stack_ids)
        levels = stacks.levels[element_stack + 1]
        for interval = 1:(length(levels) - 1)
            la, lb = levels[interval], levels[interval + 1]
            attribute = material(spec, plan.classes[class], 0.5 * (zs[la] + zs[lb]))
            attribute == 0 && continue
            before = length(tetrahedra)
            slices = (
                stack_slice(stacks.levels[stack_ids[1] + 1], la, lb),
                stack_slice(stacks.levels[stack_ids[2] + 1], la, lb),
                stack_slice(stacks.levels[stack_ids[3] + 1], la, lb)
            )
            # A node is hanging-free in THIS interval when its stack has no level inside it
            # (a finer stack may still have none here); the apex is the lowest global index
            # among the hanging-free bottom corners, so every adjacent face is a fan from it.
            free = map(slice -> length(slice) == 2, slices)
            if all(free)
                counts.region_prisms_plain += 1
                push_prism_fan!(
                    tetrahedra,
                    global_index(t[1], la),
                    global_index(t[2], la),
                    global_index(t[3], la),
                    global_index(t[1], lb),
                    global_index(t[2], lb),
                    global_index(t[3], lb)
                )
            else
                counts.region_prisms_hanging += 1
                counts.region_hanging_nodes += sum(length(sl) - 2 for sl in slices)
                empty!(faces)
                push!(
                    faces,
                    (global_index(t[1], la), global_index(t[2], la), global_index(t[3], la)),
                    (global_index(t[1], lb), global_index(t[2], lb), global_index(t[3], lb))
                )
                for (u, v) in ((1, 2), (2, 3), (3, 1))
                    ladder_triangles!(faces, global_index, t[u], t[v], slices[u], slices[v])
                end
                apex = minimum(global_index(t[u], la) for u = 1:3 if free[u])
                for (u, v, w) in faces
                    (u == apex || v == apex || w == apex) && continue
                    push!(tetrahedra, (apex, u, v, w))
                end
            end
            counts.region_tetrahedra += push_element_tetrahedra!(before, attribute, class)
        end
    end
    record = Dict{String, Any}(
        "alpha" => alpha,
        "beta" => beta,
        "rows" => rows,
        "heights_um" => heights,
        "stack_spacings_um" => stacks.spacings,
        "stack_sizes" => [length(l) for l in stacks.levels],
        "stack_levels_z_um" => [[zs[i] for i in l] for l in stacks.levels],
        "row_ranges_z_um" => [[zs[lo], zs[hi]] for (lo, hi) in stacks.ranges],
        "region_ladder_stacks" => ladder,
        "region_size_max_um" => REGION_MESH_SIZE_MAX_UM,
        "plan_nodes_by_stack" => nodes_by_stack,
        "region_plan_nodes" => region_nodes,
        "chains" => length(plan.chains),
        "columns" => columns,
        "cross_section_nodes" => length(section.nodes),
        "cross_section_triangles" => length(section.triangles),
        "cross_section_triangles_by_line_pair" => [
            Dict("lines" => collect(pair), "triangles" => n) for
            (pair, n) in section.pair_triangles
        ],
        "band_prisms" => counts.band_prisms,
        "band_tetrahedra" => counts.band_tetrahedra,
        "region_prisms_plain" => counts.region_prisms_plain,
        "region_prisms_hanging" => counts.region_prisms_hanging,
        "region_hanging_nodes" => counts.region_hanging_nodes,
        "region_tetrahedra" => counts.region_tetrahedra
    )
    verbose && println(
        "Graded sweep: alpha ",
        alpha,
        ", beta ",
        beta,
        ", stacks ",
        record["stack_sizes"],
        ", row ranges ",
        record["row_ranges_z_um"],
        ", cross-section ",
        length(section.nodes),
        " nodes / ",
        length(section.triangles),
        " triangles, band prisms ",
        counts.band_prisms,
        ", region prisms ",
        counts.region_prisms_plain,
        " plain + ",
        counts.region_prisms_hanging,
        " hanging"
    )
    return tetrahedra, tetrahedron_attribute, tetrahedron_class, record
end
