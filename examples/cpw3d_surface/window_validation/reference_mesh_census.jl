# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Element-quality census of a tensor-sweep reference mesh (fabricated window or transmon),
# attributed to the construction: every tetrahedron is bucketed by its plan distance to the
# nearest metal EDGE (from the vertical metal sidewall faces of the mesh itself, so no polygon
# set is needed) and by its z distance to the nearest metal PLANE (the metal slabs are
# recovered from the sidewall z extents), and by the z spacing of the sweep interval it came
# from. Reports minSICN / gamma / edge-length ratio percentiles per bucket, the share of
# tetrahedra (= the share of DOFs at fixed order, to within the near-constant DOF per tet of
# a conforming p-refinement) of the far-field "needles" (band rows swept far from the plane)
# and of the fine z levels swept across the coarse plan interior. Writes MESH.census.json.
#
#   julia --project=. reference_mesh_census.jl MESH.msh2 [--band-um B] [--slab-um S]
#
# B (default 1.3): plan half-width of the boundary-layer band (r (2^7 - 1) = 1.27 um at r10).
# S (default 1.3): z half-thickness of the "near plane" slab outside the metal thickness.

import Gmsh: gmsh
using JSON
using Statistics
using Printf

const MESH = abspath(ARGS[1])
band_um = 1.3
slab_um = 1.3
let i = 2
    while i <= length(ARGS)
        if ARGS[i] == "--band-um"
            global band_um = parse(Float64, ARGS[i + 1])
            i += 2
        elseif ARGS[i] == "--slab-um"
            global slab_um = parse(Float64, ARGS[i + 1])
            i += 2
        else
            error("Unknown argument $(ARGS[i])")
        end
    end
end

pct(v, s) = isempty(v) ? NaN : sort(v)[clamp(round(Int, s * length(v)), 1, length(v))]
r3(x) = isnan(x) ? nothing : round(x; sigdigits=3)

function summary_row(label, q, g, ar, total)
    n = length(q)
    n == 0 && return nothing
    return Dict(
        "label" => label,
        "tets" => n,
        "share" => round(n / total; digits=4),
        "sicn_min" => r3(minimum(q)),
        "sicn_p01" => r3(pct(q, 0.01)),
        "sicn_median" => r3(pct(q, 0.5)),
        "sicn_p90" => r3(pct(q, 0.9)),
        "gamma_median" => r3(pct(g, 0.5)),
        "edge_ratio_median" => r3(pct(ar, 0.5)),
        "edge_ratio_p90" => r3(pct(ar, 0.9)),
        "edge_ratio_max" => r3(maximum(ar)),
        "share_sicn_below_0p01" => round(count(<(0.01), q) / n; digits=4),
        "share_edge_ratio_above_100" => round(count(>(100), ar) / n; digits=4)
    )
end

function print_row(row)
    row === nothing && return
    @printf(
        "%-46s %9d %6.1f%%  SICN min %-8.3g p1 %-8.3g med %-8.3g | ratio med %-7.3g p90 %-7.3g max %-8.3g | <0.01 %5.1f%%  ratio>100 %5.1f%%\n",
        row["label"],
        row["tets"],
        100 * row["share"],
        row["sicn_min"],
        row["sicn_p01"],
        row["sicn_median"],
        row["edge_ratio_median"],
        row["edge_ratio_p90"],
        row["edge_ratio_max"],
        100 * row["share_sicn_below_0p01"],
        100 * row["share_edge_ratio_above_100"]
    )
end

gmsh.initialize()
gmsh.option.set_number("General.Terminal", 0)
gmsh.open(MESH)

# Nodes.
ntags, coords, _ = gmsh.model.mesh.get_nodes()
X = reshape(coords, 3, :)
node_index = Dict{UInt64, Int32}()
sizehint!(node_index, length(ntags))
for (i, t) in enumerate(ntags)
    node_index[t] = Int32(i)
end

# Tetrahedra.
etypes, etags, enodes = gmsh.model.mesh.get_elements(3)
all(etypes .== 4) || error("Expected first-order tetrahedra only, got types $etypes")
tags = reduce(vcat, etags)
conn = Int32[node_index[k] for k in reduce(vcat, enodes)]
nt = length(tags)
println("mesh ", basename(MESH), "  nodes ", length(ntags), "  tets ", nt)
q = gmsh.model.mesh.get_element_qualities(tags, "minSICN")
g = gmsh.model.mesh.get_element_qualities(tags, "gamma")

# Metal sidewalls -> metal edge segments in plan and the metal z slabs. A sidewall face is a
# tagged 2D face whose three nodes have only two distinct plan positions (metal and trench
# sidewalls, bump columns). The slabs are built from the faces of the r-resolved levels only
# (z extent <= SLAB_FACE_DZ_UM), so bump columns spanning the gap do not merge the planes.
const SLAB_FACE_DZ_UM = 0.02
segments = Dict{NTuple{4, Float64}, Nothing}()
metal_z = Dict{Int, Vector{Float64}}() # physical attribute -> z values of sidewall nodes
for (dim, attribute) in gmsh.model.get_physical_groups(2)
    name = gmsh.model.get_physical_name(dim, attribute)
    (endswith(name, "_air") || endswith(name, "_substrate")) || continue
    for entity in gmsh.model.get_entities_for_physical_group(dim, attribute)
        ftypes, _, fnodes = gmsh.model.mesh.get_elements(dim, entity)
        for (ft, fn) in zip(ftypes, fnodes)
            ft == 2 || continue
            for f = 1:3:length(fn)
                a, b, c = node_index[fn[f]], node_index[fn[f + 1]], node_index[fn[f + 2]]
                pa, pb, pc = (X[1, a], X[2, a]), (X[1, b], X[2, b]), (X[1, c], X[2, c])
                same(u, v) = abs(u[1] - v[1]) < 1.0e-9 && abs(u[2] - v[2]) < 1.0e-9
                pair = same(pa, pb) ? (pa, pc) : same(pa, pc) ? (pa, pb) :
                       same(pb, pc) ? (pa, pb) : nothing
                pair === nothing && continue
                u, v = pair
                key = u < v ? (u[1], u[2], v[1], v[2]) : (v[1], v[2], u[1], u[2])
                segments[key] = nothing
                zlo, zhi = min(X[3, a], X[3, b], X[3, c]), max(X[3, a], X[3, b], X[3, c])
                zhi - zlo <= SLAB_FACE_DZ_UM || continue
                zs = get!(metal_z, attribute, Float64[])
                push!(zs, zlo, zhi)
            end
        end
    end
end
segs = collect(keys(segments))
println("metal edge segments (from sidewalls): ", length(segs))
# Metal slabs in z: the distinct sidewall z ranges, merged.
# Slabs: clusters of the sidewall z values separated by more than 0.5 um.
zvals = sort(unique(round.(reduce(vcat, values(metal_z)); digits=6)))
slabs = Tuple{Float64, Float64}[]
let start = zvals[1]
    for i = 2:length(zvals)
        if zvals[i] - zvals[i - 1] > 0.5
            push!(slabs, (start, zvals[i - 1]))
            start = zvals[i]
        end
    end
    push!(slabs, (start, zvals[end]))
end
println("metal slabs (z): ", slabs)

# Plan distance of every distinct plan node position to the nearest metal edge segment
# (the sweep reuses the plan positions at every level: compute once per plan position).
function segment_distance(px, py, s)
    ax, ay, bx, by = s
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = L2 > 0 ? clamp(((px - ax) * dx + (py - ay) * dy) / L2, 0.0, 1.0) : 0.0
    return hypot(px - ax - t * dx, py - ay - t * dy)
end
plan_positions = Dict{NTuple{2, Float64}, Int}()
plan_of_node = Vector{Int}(undef, length(ntags))
for i = 1:length(ntags)
    key = (round(X[1, i]; digits=9), round(X[2, i]; digits=9))
    plan_of_node[i] = get!(plan_positions, key, length(plan_positions) + 1)
end
println("distinct plan positions: ", length(plan_positions))
plan_distance = fill(Inf, length(plan_positions))
# Uniform grid over the segments for the nearest-segment query.
xs = [s[1] for s in segs]
ys = [s[2] for s in segs]
xmin, xmax = minimum(min.(xs, [s[3] for s in segs])), maximum(max.(xs, [s[3] for s in segs]))
ymin, ymax = minimum(min.(ys, [s[4] for s in segs])), maximum(max.(ys, [s[4] for s in segs]))
cell = max(5.0, (xmax - xmin) / 200)
nx = max(1, ceil(Int, (xmax - xmin) / cell))
ny = max(1, ceil(Int, (ymax - ymin) / cell))
grid = [Int[] for _ = 1:nx, _ = 1:ny]
for (k, s) in enumerate(segs)
    i0 = clamp(floor(Int, (min(s[1], s[3]) - xmin) / cell) + 1, 1, nx)
    i1 = clamp(floor(Int, (max(s[1], s[3]) - xmin) / cell) + 1, 1, nx)
    j0 = clamp(floor(Int, (min(s[2], s[4]) - ymin) / cell) + 1, 1, ny)
    j1 = clamp(floor(Int, (max(s[2], s[4]) - ymin) / cell) + 1, 1, ny)
    for i = i0:i1, j = j0:j1
        push!(grid[i, j], k)
    end
end
for (key, index) in plan_positions
    px, py = key
    ci = clamp(floor(Int, (px - xmin) / cell) + 1, 1, nx)
    cj = clamp(floor(Int, (py - ymin) / cell) + 1, 1, ny)
    best = Inf
    ring = 0
    while true
        for i = max(1, ci - ring):min(nx, ci + ring), j = max(1, cj - ring):min(ny, cj + ring)
            (abs(i - ci) == ring || abs(j - cj) == ring) || continue
            for k in grid[i, j]
                d = segment_distance(px, py, segs[k])
                d < best && (best = d)
            end
        end
        # Every segment closer than (ring) cells has been seen once ring * cell >= best.
        (best <= ring * cell || (ci - ring < 1 && ci + ring > nx && cj - ring < 1 && cj + ring > ny)) && break
        ring += 1
    end
    plan_distance[index] = best
end

# Per tetrahedron: z centroid, z extent (the sweep interval), plan distance (mean of the
# vertex distances), edge-length ratio.
zc = zeros(nt)
dz = zeros(nt)
dedge = zeros(nt)
ar = zeros(nt)
for e = 1:nt
    n1, n2, n3, n4 = conn[4e - 3], conn[4e - 2], conn[4e - 1], conn[4e]
    z1, z2, z3, z4 = X[3, n1], X[3, n2], X[3, n3], X[3, n4]
    zc[e] = 0.25 * (z1 + z2 + z3 + z4)
    dz[e] = max(z1, z2, z3, z4) - min(z1, z2, z3, z4)
    dedge[e] =
        0.25 * (
            plan_distance[plan_of_node[n1]] +
            plan_distance[plan_of_node[n2]] +
            plan_distance[plan_of_node[n3]] +
            plan_distance[plan_of_node[n4]]
        )
    lmin, lmax = Inf, 0.0
    for (a, b) in ((n1, n2), (n1, n3), (n1, n4), (n2, n3), (n2, n4), (n3, n4))
        L = hypot(X[1, a] - X[1, b], X[2, a] - X[2, b], X[3, a] - X[3, b])
        lmin = min(lmin, L)
        lmax = max(lmax, L)
    end
    ar[e] = lmax / lmin
end
function plane_distance(z)
    best = Inf
    for (a, b) in slabs
        best = min(best, z < a ? a - z : z > b ? z - b : 0.0)
    end
    return best
end
dplane = plane_distance.(zc)

total = nt
rows = Any[]
push!(rows, summary_row("ALL", q, g, ar, total))
println()
print_row(rows[end])

println("\n-- by plan distance to the nearest metal edge (um) --")
edge_bins = [(0.0, 0.1, "edge < 0.1 (rows 1-4)"), (0.1, band_um, "0.1 <= edge < $(band_um) (outer band)"),
    (band_um, 10.0, "$(band_um) <= edge < 10"), (10.0, Inf, "edge >= 10 (plan interior)")]
for (lo, hi, label) in edge_bins
    m = (dedge .>= lo) .& (dedge .< hi)
    row = summary_row(label, q[m], g[m], ar[m], total)
    push!(rows, row)
    print_row(row)
end

println("\n-- by z distance to the nearest metal slab (um) --")
plane_bins = [(0.0, 1.0e-12, "inside the metal thickness range"), (1.0e-12, 0.1, "plane < 0.1"),
    (0.1, slab_um, "0.1 <= plane < $(slab_um)"), (slab_um, 10.0, "$(slab_um) <= plane < 10"),
    (10.0, Inf, "plane >= 10 (far field)")]
for (lo, hi, label) in plane_bins
    m = (dplane .>= lo) .& (dplane .< hi)
    row = summary_row(label, q[m], g[m], ar[m], total)
    push!(rows, row)
    print_row(row)
end

println("\n-- cross table: band (edge < $(band_um)) x near plane (plane < $(slab_um)) --")
band = dedge .< band_um
near = dplane .< slab_um
for (m, label) in (
    (band .& near, "band x near plane (the resolved region)"),
    (band .& .!near, "band x far from plane (FAR-FIELD NEEDLES)"),
    (.!band .& near, "interior x near plane (fine z x coarse plan)"),
    (.!band .& .!near, "interior x far from plane")
)
    row = summary_row(label, q[m], g[m], ar[m], total)
    push!(rows, row)
    print_row(row)
end

println("\n-- by sweep interval height dz (um) --")
dz_bins = [(0.0, 0.015, "dz <= 0.01 (metal / trench levels)"), (0.015, 0.11, "0.01 < dz <= 0.1"),
    (0.11, 1.1, "0.1 < dz <= 1"), (1.1, 11.0, "1 < dz <= 10"), (11.0, 110.0, "10 < dz <= 100"),
    (110.0, Inf, "dz > 100")]
for (lo, hi, label) in dz_bins
    m = (dz .> lo) .& (dz .<= hi)
    row = summary_row(label, q[m], g[m], ar[m], total)
    push!(rows, row)
    print_row(row)
end

# Needle attribution: the share of tets whose edge ratio exceeds 100 by construction class.
needles = ar .> 100
println("\n-- tets with edge ratio > 100: ", count(needles), " (", round(100 * count(needles) / nt; digits=1), " %) --")
attribution = Dict{String, Any}()
for (m, label) in (
    (band .& near, "band x near plane"),
    (band .& .!near, "band x far from plane"),
    (.!band .& near, "interior x near plane"),
    (.!band .& .!near, "interior x far from plane")
)
    c = count(needles .& m)
    attribution[label] = Dict("needles" => c, "share_of_needles" => round(c / max(1, count(needles)); digits=4))
    println(rpad(label, 40), lpad(c, 9), "  ", round(100 * c / max(1, count(needles)); digits=1), " % of the needles")
end

levels = sort(unique(round.(X[3, :]; digits=6)))
result = Dict(
    "mesh" => MESH,
    "nodes" => length(ntags),
    "tets" => nt,
    "z_levels" => length(levels),
    "metal_slabs_z" => slabs,
    "metal_edge_segments" => length(segs),
    "plan_positions" => length(plan_positions),
    "band_um" => band_um,
    "slab_um" => slab_um,
    "rows" => filter(!isnothing, rows),
    "needle_attribution" => attribution,
    "share_sicn_below_0p01" => round(count(<(0.01), q) / nt; digits=4),
    "share_sicn_below_0p001" => round(count(<(0.001), q) / nt; digits=4),
    "share_edge_ratio_above_100" => round(count(>(100), ar) / nt; digits=4),
    "share_edge_ratio_above_1000" => round(count(>(1000), ar) / nt; digits=4)
)
open(replace(MESH, r"\.msh2?$" => "") * ".census.json", "w") do io
    JSON.print(io, result, 2)
end
println("\nz levels ", length(levels), "  wrote ", replace(MESH, r"\.msh2?$" => "") * ".census.json")
gmsh.finalize()
