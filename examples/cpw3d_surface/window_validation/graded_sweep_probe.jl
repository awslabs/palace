# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Lane Q2 feasibility PROBE (decision 257; design option (v) "graded cross-section sweep"):
# a straight metal strip across a one-plane window, meshed with the production band columns
# (t-spaced stations along the edge, rows h_k = r (2^k - 1) on both sides of each edge) but
# with a GRADED (n, z) cross-section instead of the tensor product of the plan with every z
# level: band row k carries only the z levels of a nested stack thinned to a spacing >=
# alpha * w_k (w_k = r 2^(k-1), the row's width; the material interfaces are always kept), and
# the row line ends at |z| > beta * h_k (its nodes are not swept further), so far from the plane
# the band collapses onto its base line and its outermost line. The cross-section is a 2D
# triangulation of the resulting (y, z) point set (fans from a hanging-free corner of every
# cell), swept along x into prisms, split into tetrahedra with the production diagonal rule.
# Surfaces are tagged by position (box wall 3, ground_air 4, ground_substrate 5, substrate_air
# 6) exactly as the production attribute table so validate_window_mesh.jl applies unchanged.
# PROBE ONLY: one plane, straight edges, no corners / fans / caps / bumps, no Gmsh.
#
#   julia --project=. graded_sweep_probe.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2
#       [--alpha A] [--beta B] [--interior-max-um M] [--far-growth G]
#
# alpha: z spacing of row k >= alpha * w_k; beta: row k ends at |z| > beta * h_k beyond the metal;
# interior lines grow from 2 w_K by far_growth per line up to interior_max_um (the production
# interior is Gmsh-triangulated at <= 30 um; the probe's structured interior stands in for it).

using JSON
using Printf
using SHA

positional = String[]
alpha = 1.0
beta = 3.0
interior_max_um = 30.0
far_growth = 2.0
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--alpha"
            global alpha = parse(Float64, ARGS[i + 1])
            i += 2
        elseif a == "--beta"
            global beta = parse(Float64, ARGS[i + 1])
            i += 2
        elseif a == "--interior-max-um"
            global interior_max_um = parse(Float64, ARGS[i + 1])
            i += 2
        elseif a == "--far-growth"
            global far_growth = parse(Float64, ARGS[i + 1])
            i += 2
        else
            push!(positional, a)
            i += 1
        end
    end
end
length(positional) == 4 || error("Usage: graded_sweep_probe.jl WINDOW.json R T OUTPUT.msh2")
spec = JSON.parsefile(positional[1])
radial_um = parse(Float64, positional[2])
tangential_um = parse(Float64, positional[3])
output = abspath(positional[4])

box_x = Float64.(spec["Box"]["X"])
box_y = Float64.(spec["Box"]["Y"])
metal_thickness = Float64(get(get(spec, "Process", Dict()), "MetalThickness", 0.1))
overetch = Float64(get(get(spec, "Process", Dict()), "Overetch", 0.05))
length(spec["Planes"]) == 1 || error("Probe: one plane only")
plane = spec["Planes"][1]
plane["Facing"] == "up" || error("Probe: plane facing up only")
surface_z = Float64(plane["SurfaceZ"])
substrate_thickness = Float64(plane["SubstrateThickness"])
vacuum_above = Float64(spec["Vacuum"]["Above"])
vacuum_below = Float64(get(spec["Vacuum"], "Below", 0.0))
vacuum_below == 0.0 || error("Probe: Vacuum.Below 0 only (backside = box wall)")
length(plane["Polygons"]) == 1 || error("Probe: one strip polygon only")
poly = plane["Polygons"][1]
poly["Conductor"] == "ground" || error("Probe: the strip is the ground")
outer = [Float64.(p) for p in poly["Outer"]]
xs_p = [p[1] for p in outer]
ys_p = [p[2] for p in outer]
(minimum(xs_p) == box_x[1] && maximum(xs_p) == box_x[2]) ||
    error("Probe: the strip must span the box in x")
y_metal = (minimum(ys_p), maximum(ys_p))

# ---------------------------------------------------------------------------------------------
# Production z levels (PolygonWindowMesh.z_levels for one plane facing up) and band heights.
const SUBSTRATE_SIDE_OFFSETS_UM = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
const VACUUM_SIDE_OFFSETS_UM =
    [0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]
rows_max = max(1, round(Int, log2(1.0 + 1.27 / 0.01))) # 7 rows as the production band at r10
rows_max = 7
heights = [radial_um * (2.0^k - 1.0) for k = 1:rows_max]
widths = [radial_um * 2.0^(k - 1) for k = 1:rows_max]
metal_layers = max(2, ceil(Int, metal_thickness / radial_um - 1.0e-9))
trench_layers = max(1, ceil(Int, overetch / radial_um - 1.0e-9))
z_all = Float64[]
append!(z_all, range(surface_z - overetch, surface_z; length=trench_layers + 1))
append!(z_all, range(surface_z, surface_z + metal_thickness; length=metal_layers + 1))
for offset in SUBSTRATE_SIDE_OFFSETS_UM
    offset < substrate_thickness || break
    push!(z_all, surface_z - offset)
end
push!(z_all, surface_z - substrate_thickness)
for offset in VACUUM_SIDE_OFFSETS_UM
    push!(z_all, surface_z + offset)
end
z_bottom = surface_z - substrate_thickness - vacuum_below
z_top = surface_z + vacuum_above
push!(z_all, z_bottom, z_top)
filter!(z -> z_bottom - 1.0e-9 <= z <= z_top + 1.0e-9, z_all)
sort!(z_all)
z_all = [z_all[1]; [z_all[i] for i = 2:length(z_all) if z_all[i] - z_all[i - 1] > 1.0e-9]]
interfaces = [surface_z - overetch, surface_z, surface_z + metal_thickness, z_bottom, z_top]
is_interface(z) = any(abs(z - w) <= 1.0e-9 for w in interfaces)

# Nested thinning: keep interfaces; otherwise keep a level only if it is at least `spacing`
# from the last kept level AND from the next interface above (so no sub-spacing sliver).
function thin(levels::Vector{Float64}, spacing::Float64)
    kept = Float64[levels[1]]
    for i = 2:length(levels)
        z = levels[i]
        if is_interface(z)
            push!(kept, z)
        else
            next_interface =
                findfirst(j -> j > i && is_interface(levels[j]), eachindex(levels))
            room = next_interface === nothing ? Inf : levels[next_interface] - z
            if z - kept[end] >= spacing - 1.0e-12 && room >= spacing - 1.0e-12
                push!(kept, z)
            end
        end
    end
    return kept
end
# Row k's stack: thinned from row k-1's (nested); row 0 (the edge base line) = every level.
stacks = Vector{Vector{Float64}}(undef, rows_max + 1)
stacks[1] = copy(z_all)
for k = 1:rows_max
    stacks[k + 1] = thin(stacks[k], alpha * widths[k])
end
# Interior lines (plan distance > h_K from every edge): the outermost row's stack thinned to
# their own width (capped by the far spacing rule below).
interior_stack(width) = thin(stacks[rows_max + 1], alpha * min(width, interior_max_um))

# Row line k is active for z in [lo_k, hi_k]: the first level of the NEXT row's stack beyond
# surface_z -/+ beta * h_k (so the line's end node lies on a level of the coarser neighbour
# it hangs from), strictly nested (inner rows end first); the outermost row spans everything.
ranges = Vector{Tuple{Float64, Float64}}(undef, rows_max)
for k = 1:rows_max
    if k == rows_max
        ranges[k] = (z_bottom, z_top)
        continue
    end
    s = stacks[k + 2]
    limit_lo = surface_z - beta * heights[k]
    limit_hi = surface_z + metal_thickness + beta * heights[k]
    if k > 1
        limit_lo = min(limit_lo, ranges[k - 1][1] - 1.0e-9)
        limit_hi = max(limit_hi, ranges[k - 1][2] + 1.0e-9)
    end
    lo = maximum(z for z in s if z <= limit_lo + 1.0e-12; init=s[1])
    hi = minimum(z for z in s if z >= limit_hi - 1.0e-12; init=s[end])
    ranges[k] = (lo, hi)
end

# ---------------------------------------------------------------------------------------------
# Plan y lines. Each line: y, stack of z levels, active range, and a material side tag.
struct YLine
    y::Float64
    stack::Vector{Float64}
    lo::Float64
    hi::Float64
end
lines = YLine[]
push_line!(y, stack, lo, hi) = push!(lines, YLine(y, stack, lo, hi))
# Interior graded lines from a band top outward to a wall (or to the strip middle), widths
# growing geometrically from the outermost row's width to interior_max_um.
function interior_positions(from, to)
    direction = sign(to - from)
    positions = Float64[]
    w = widths[rows_max] * 2
    y = from
    while true
        remaining = abs(to - y)
        if remaining <= 1.5 * min(w, interior_max_um)
            break
        end
        y += direction * min(w, interior_max_um)
        push!(positions, y)
        w *= far_growth
    end
    return positions
end
# Gap below the strip: wall .. band of edge y_metal[1]; the strip; gap above.
e1, e2 = y_metal
# Edge lines (base), band lines on both sides of both edges.
for (edge, side) in ((e1, -1), (e1, +1), (e2, -1), (e2, +1))
    for k = 1:rows_max
        push_line!(edge + side * heights[k], stacks[k + 1], ranges[k]...)
    end
end
push_line!(e1, stacks[1], z_bottom, z_top)
push_line!(e2, stacks[1], z_bottom, z_top)
# Interior lines: gap 1 (wall y0 .. e1 - h_K), metal (e1 + h_K .. e2 - h_K), gap 2.
for (from, to) in (
    (e1 - heights[end], box_y[1]),
    (e2 + heights[end], box_y[2]),
    (e1 + heights[end], 0.5 * (e1 + e2)),
    (e2 - heights[end], 0.5 * (e1 + e2))
)
    for (i, y) in enumerate(interior_positions(from, to))
        push_line!(y, interior_stack(widths[rows_max] * far_growth^i), z_bottom, z_top)
    end
end
push_line!(box_y[1], interior_stack(interior_max_um), z_bottom, z_top)
push_line!(box_y[2], interior_stack(interior_max_um), z_bottom, z_top)
if 0.5 * (e1 + e2) > e1 + heights[end]
    push_line!(0.5 * (e1 + e2), interior_stack(interior_max_um), z_bottom, z_top)
end
sort!(lines; by=l -> l.y)
ys = [l.y for l in lines]
all(diff(ys) .> 1.0e-9) || error("Duplicate y lines")
println("y lines: ", length(lines), "  z levels: ", length(z_all))
for k = 0:rows_max
    @printf(
        "  row %d  width %.3g  stack %d levels%s\n",
        k,
        k == 0 ? 0.0 : widths[k],
        length(stacks[k + 1]),
        k == 0 ? "" : @sprintf("  active z in [%.3g, %.3g]", ranges[k]...)
    )
end

# ---------------------------------------------------------------------------------------------
# 2D cross-section nodes: (line index, z) with z in the line's stack and active range.
node_id = Dict{Tuple{Int, Float64}, Int}()
node_yz = Tuple{Float64, Float64}[]
function cs_node(li, z)
    return get!(node_id, (li, round(z; digits=9))) do
        push!(node_yz, (lines[li].y, z))
        return length(node_yz)
    end
end
active_levels(li) =
    [z for z in lines[li].stack if lines[li].lo - 1.0e-9 <= z <= lines[li].hi + 1.0e-9]
for li in eachindex(lines), z in active_levels(li)
    cs_node(li, z)
end
println("cross-section nodes: ", length(node_yz))

# Cells. The active lines change only where a band row ends, so the cross-section is
# enumerated per PAIR of lines that are neighbours over a z range: (base, row k) over the range
# where rows < k have ended, (row k, row k+1) over row k's range, (row K, nearest interior
# line) and interior pairs over everything. Within a pair's range the cells are the consecutive
# levels of the COARSER stack (nested); the finer line's extra levels are hanging nodes on one
# vertical edge, the end node of the row that ended at the cell's bottom (started at its top)
# is a hanging node on that horizontal edge; every cell is fanned from the coarse line's node
# at the hanging-free horizontal edge, so no fan triangle is degenerate.
triangles = NTuple{3, Int}[]
orient2(a, b, c) =
    (node_yz[b][1] - node_yz[a][1]) * (node_yz[c][2] - node_yz[a][2]) -
    (node_yz[b][2] - node_yz[a][2]) * (node_yz[c][1] - node_yz[a][1])
function fan!(apex, chain)
    for i = 1:(length(chain) - 1)
        a, b, c = apex, chain[i], chain[i + 1]
        orient2(a, b, c) > 1.0e-18 || ((a, b, c) = (a, c, b))
        orient2(a, b, c) > 1.0e-18 || error(
            "Degenerate cross-section triangle $(node_yz[a]) $(node_yz[b]) $(node_yz[c])"
        )
        push!(triangles, (a, b, c))
    end
end
levels_in(li, za, zb) = [z for z in lines[li].stack if za - 1.0e-9 <= z <= zb + 1.0e-9]
function cells_between!(L, R, za, zb, bottom_hanging, top_hanging)
    zb - za > 1.0e-12 || return
    levels_L, levels_R = levels_in(L, za, zb), levels_in(R, za, zb)
    coarse, fine, cl, fl =
        length(levels_L) <= length(levels_R) ? (L, R, levels_L, levels_R) :
        (R, L, levels_R, levels_L)
    all(any(abs(z - w) <= 1.0e-9 for w in fl) for z in cl) ||
        error("Stacks not nested between y=$(lines[L].y) and y=$(lines[R].y) in [$za, $zb]")
    (abs(cl[1] - za) <= 1.0e-9 && abs(cl[end] - zb) <= 1.0e-9) || error(
        "Range ends [$za, $zb] are not levels of the coarse line y=$(lines[coarse].y)"
    )
    left = min(L, R)
    for ci = 1:(length(cl) - 1)
        z0, z1 = cl[ci], cl[ci + 1]
        fine_levels = [z for z in fl if z0 - 1.0e-9 <= z <= z1 + 1.0e-9]
        bottom = ci == 1 ? bottom_hanging : Int[]
        top = ci == length(cl) - 1 ? top_hanging : Int[]
        (isempty(bottom) || isempty(top)) ||
            error("Cell with hanging nodes on both horizontal edges at z=[$z0, $z1]")
        chain = Int[]
        if isempty(top)
            apex = cs_node(coarse, z1)
            if coarse == left
                push!(chain, cs_node(coarse, z0))
                for li in sort(bottom; by=li -> lines[li].y)
                    push!(chain, cs_node(li, z0))
                end
                for z in fine_levels
                    push!(chain, cs_node(fine, z))
                end
            else
                for z in reverse(fine_levels)
                    push!(chain, cs_node(fine, z))
                end
                for li in sort(bottom; by=li -> -lines[li].y)
                    push!(chain, cs_node(li, z0))
                end
                push!(chain, cs_node(coarse, z0))
            end
        else
            apex = cs_node(coarse, z0)
            if coarse == left
                for z in fine_levels
                    push!(chain, cs_node(fine, z))
                end
                for li in sort(top; by=li -> -lines[li].y)
                    push!(chain, cs_node(li, z1))
                end
                push!(chain, cs_node(coarse, z1))
            else
                push!(chain, cs_node(coarse, z1))
                for li in sort(top; by=li -> lines[li].y)
                    push!(chain, cs_node(li, z1))
                end
                for z in reverse(fine_levels)
                    push!(chain, cs_node(fine, z))
                end
            end
        end
        fan!(apex, chain)
    end
end
line_at(y) = findfirst(l -> abs(l.y - y) <= 1.0e-9, lines)
for (edge, side) in ((e1, -1), (e1, +1), (e2, -1), (e2, +1))
    base = line_at(edge)
    row = [line_at(edge + side * heights[k]) for k = 1:rows_max]
    # (base, row 1) over row 1's range; (base, row k) over the parts of row k's range beyond
    # row k-1's, with row k-1's end node hanging on the shared horizontal edge.
    cells_between!(base, row[1], ranges[1][1], ranges[1][2], Int[], Int[])
    for k = 2:rows_max
        cells_between!(base, row[k], ranges[k - 1][2], ranges[k][2], [row[k - 1]], Int[])
        cells_between!(base, row[k], ranges[k][1], ranges[k - 1][1], Int[], [row[k - 1]])
    end
    for k = 1:(rows_max - 1)
        cells_between!(row[k], row[k + 1], ranges[k][1], ranges[k][2], Int[], Int[])
    end
end
# Full-range pairs: the outermost rows with their interior neighbours and the interior lines.
full = [
    li for li in eachindex(lines) if
    lines[li].lo <= z_bottom + 1.0e-9 && lines[li].hi >= z_top - 1.0e-9
]
for i = 1:(length(full) - 1)
    L, R = full[i], full[i + 1]
    # Skip (base, row K) pairs: they are covered above (the base lines are full-range too).
    any(abs(lines[L].y - e) <= 1.0e-9 || abs(lines[R].y - e) <= 1.0e-9 for e in (e1, e2)) &&
        continue
    cells_between!(L, R, z_bottom, z_top, Int[], Int[])
end
println("cross-section triangles: ", length(triangles))
# Area check of the cross-section.
cs_area = sum(0.5 * orient2(t...) for t in triangles)
expected = (box_y[2] - box_y[1]) * (z_top - z_bottom)
abs(cs_area - expected) <= 1.0e-9 * expected ||
    error("Cross-section area $cs_area != $expected (gap $(cs_area - expected))")
# Node use check: every cross-section node is used.
used_cs = falses(length(node_yz))
for t in triangles, n in t
    used_cs[n] = true
end
all(used_cs) || error("$(count(!, used_cs)) unused cross-section nodes")

# Cross-section census: triangles by region (band = within h_K of an edge in y; near = within
# h_K of the metal slab in z), for the comparison with the tensor sweep's prisms per column.
let band_near = 0, band_far = 0, int_near = 0, int_far = 0
    for t in triangles
        yc = sum(node_yz[n][1] for n in t) / 3
        zc = sum(node_yz[n][2] for n in t) / 3
        band = min(abs(yc - e1), abs(yc - e2)) <= heights[end]
        near = surface_z - heights[end] <= zc <= surface_z + metal_thickness + heights[end]
        band && near && (band_near += 1)
        band && !near && (band_far += 1)
        !band && near && (int_near += 1)
        !band && !near && (int_far += 1)
    end
    println(
        "cross-section triangles: band x near ",
        band_near,
        "  band x far ",
        band_far,
        "  interior x near ",
        int_near,
        "  interior x far ",
        int_far,
        " (tensor sweep per column pair: band rows x intervals = ",
        rows_max,
        " x ",
        length(z_all) - 1,
        " x 2 = ",
        2 * rows_max * (length(z_all) - 1),
        " per side)"
    )
end

# Material of a cross-section triangle by its centroid: 0 excluded (metal), 1 substrate, 2 vacuum.
function material(y, z)
    in_metal_plan = e1 < y < e2
    d = z - surface_z
    if d < -overetch
        return Int8(1)
    elseif d < 0.0
        return in_metal_plan ? Int8(1) : Int8(2)
    elseif d < metal_thickness
        return in_metal_plan ? Int8(0) : Int8(2)
    end
    return Int8(2)
end
cs_material = [
    material(sum(node_yz[n][1] for n in t) / 3, sum(node_yz[n][2] for n in t) / 3) for
    t in triangles
]

# ---------------------------------------------------------------------------------------------
# Sweep along x: stations every tangential_um; node (station s, cs node n).
nx = ceil(Int, (box_x[2] - box_x[1]) / tangential_um - 1.0e-9)
stations = range(box_x[1], box_x[2]; length=nx + 1)
n_cs = length(node_yz)
nodes = Vector{NTuple{3, Float64}}(undef, (nx + 1) * n_cs)
for s = 0:nx, n = 1:n_cs
    nodes[s * n_cs + n] = (stations[s + 1], node_yz[n][1], node_yz[n][2])
end
tetrahedra = NTuple{4, Int32}[]
tetrahedron_attribute = Int8[]
function push_prism!(a1, a2, a3, b1, b2, b3, attribute)
    if a2 < a1 && a2 <= a3
        a1, a2, a3 = a2, a3, a1
        b1, b2, b3 = b2, b3, b1
    elseif a3 < a1 && a3 <= a2
        a1, a2, a3 = a3, a1, a2
        b1, b2, b3 = b3, b1, b2
    end
    if a2 < a3
        push!(tetrahedra, (a1, a2, a3, b3), (a1, a2, b3, b2), (a1, b2, b3, b1))
    else
        push!(tetrahedra, (a1, a2, a3, b2), (a1, a3, b3, b2), (a1, b2, b3, b1))
    end
    push!(tetrahedron_attribute, attribute, attribute, attribute)
    return nothing
end
for (ti, t) in enumerate(triangles)
    attribute = cs_material[ti]
    attribute == 0 && continue
    for s = 0:(nx - 1)
        lower = Int32(s * n_cs)
        upper = Int32((s + 1) * n_cs)
        push_prism!(
            lower + t[1],
            lower + t[2],
            lower + t[3],
            upper + t[1],
            upper + t[2],
            upper + t[3],
            attribute
        )
    end
end
function signed_volume(p1, p2, p3, p4)
    ax, ay, az = p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]
    bx, by, bz = p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3]
    cx, cy, cz = p4[1] - p1[1], p4[2] - p1[2], p4[3] - p1[3]
    return (
        ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)
    ) / 6
end
for index in eachindex(tetrahedra)
    t = tetrahedra[index]
    v = signed_volume(nodes[t[1]], nodes[t[2]], nodes[t[3]], nodes[t[4]])
    if v < 0
        tetrahedra[index] = (t[1], t[3], t[2], t[4])
        v = -v
    end
    v > 0 || error("Degenerate tetrahedron $index")
end
println("stations: ", nx + 1, "  tetrahedra: ", length(tetrahedra))

# Faces and position-based surface tags.
face_data = Dict{NTuple{3, Int32}, Tuple{Int8, Int8}}()
sizehint!(face_data, 2 * length(tetrahedra))
for (index, t) in enumerate(tetrahedra)
    attribute = tetrahedron_attribute[index]
    for face in
        ((t[1], t[2], t[3]), (t[1], t[2], t[4]), (t[1], t[3], t[4]), (t[2], t[3], t[4]))
        a, b, c = face
        a > b && ((a, b) = (b, a))
        b > c && ((b, c) = (c, b))
        a > b && ((a, b) = (b, a))
        previous = get(face_data, (a, b, c), (Int8(0), Int8(0)))
        face_data[(a, b, c)] = (previous[1] + Int8(1), previous[2] + attribute)
    end
end
surface_elements = Tuple{Int, NTuple{3, Int32}}[]
near(a, b) = abs(a - b) <= 1.0e-6
for (face, (count, attribute_sum)) in face_data
    p = (nodes[face[1]], nodes[face[2]], nodes[face[3]])
    if count == 1
        attribute = 0
        if all(near(q[1], box_x[1]) for q in p) ||
           all(near(q[1], box_x[2]) for q in p) ||
           all(near(q[2], box_y[1]) for q in p) ||
           all(near(q[2], box_y[2]) for q in p) ||
           all(near(q[3], z_bottom) for q in p) ||
           all(near(q[3], z_top) for q in p)
            attribute = 3
        elseif all(near(q[3], surface_z + metal_thickness) for q in p)
            attribute = 4 # metal top
        elseif all(near(q[3], surface_z) for q in p)
            attribute = 5 # metal bottom
        elseif all(near(q[2], e1) for q in p) || all(near(q[2], e2) for q in p)
            zmid = sum(q[3] for q in p) / 3
            surface_z < zmid < surface_z + metal_thickness ||
                error("Sidewall face off the metal at z=$zmid")
            attribute = 4 # sidewall
        else
            error("Unclassified boundary face at $(p[1])")
        end
        push!(surface_elements, (attribute, face))
    elseif count == 2 && attribute_sum == 3
        push!(surface_elements, (6, face))
    end
end
# Compact unused nodes (inside the metal).
used = falses(length(nodes))
for t in tetrahedra, n in t
    used[n] = true
end
remap = zeros(Int32, length(nodes))
compacted = NTuple{3, Float64}[]
for i in eachindex(nodes)
    if used[i]
        push!(compacted, nodes[i])
        remap[i] = Int32(length(compacted))
    end
end
tetrahedra = [(remap[t[1]], remap[t[2]], remap[t[3]], remap[t[4]]) for t in tetrahedra]
surface_elements =
    [(a, (remap[f[1]], remap[f[2]], remap[f[3]])) for (a, f) in surface_elements]
nodes = compacted
face_area(p1, p2, p3) = begin
    ux, uy, uz = p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]
    vx, vy, vz = p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3]
    0.5 * hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
end
surface_counts = Dict{String, Int}()
surface_areas = Dict{String, Float64}()
for (a, f) in surface_elements
    surface_counts[string(a)] = get(surface_counts, string(a), 0) + 1
    surface_areas[string(a)] =
        get(surface_areas, string(a), 0.0) +
        face_area(nodes[f[1]], nodes[f[2]], nodes[f[3]])
end
volumes = Dict{String, Float64}()
volume_counts = Dict{String, Int}()
for (i, t) in enumerate(tetrahedra)
    key = string(tetrahedron_attribute[i])
    volume_counts[key] = get(volume_counts, key, 0) + 1
    volumes[key] =
        get(volumes, key, 0.0) +
        signed_volume(nodes[t[1]], nodes[t[2]], nodes[t[3]], nodes[t[4]])
end
names = [
    (3, 1, "substrate"),
    (3, 2, "vacuum"),
    (2, 3, "exterior_boundary"),
    (2, 4, "ground_air"),
    (2, 5, "ground_substrate"),
    (2, 6, "substrate_air")
]
mkpath(dirname(output))
open(output, "w") do stream
    print(stream, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n")
    print(stream, "\$PhysicalNames\n$(length(names))\n")
    for (d, a, name) in names
        print(stream, "$d $a \"$name\"\n")
    end
    print(stream, "\$EndPhysicalNames\n\$Nodes\n$(length(nodes))\n")
    for (i, p) in enumerate(nodes)
        @printf(stream, "%d %.16g %.16g %.16g\n", i, p[1], p[2], p[3])
    end
    print(
        stream,
        "\$EndNodes\n\$Elements\n$(length(surface_elements) + length(tetrahedra))\n"
    )
    e = 0
    for (a, f) in surface_elements
        e += 1
        print(stream, "$e 2 2 $a $a $(f[1]) $(f[2]) $(f[3])\n")
    end
    for (i, t) in enumerate(tetrahedra)
        e += 1
        a = Int(tetrahedron_attribute[i])
        print(stream, "$e 4 2 $a $a $(t[1]) $(t[2]) $(t[3]) $(t[4])\n")
    end
    return print(stream, "\$EndElements\n")
end
manifest = Dict(
    "probe" => "graded cross-section sweep",
    "alpha" => alpha,
    "beta" => beta,
    "interior_max_um" => interior_max_um,
    "far_growth" => far_growth,
    "radial_target_um" => radial_um,
    "tangential_target_um" => tangential_um,
    "rows" => rows_max,
    "heights_um" => heights,
    "z_levels_all" => z_all,
    "stack_sizes_by_row" => [length(s) for s in stacks],
    "active_ranges_by_row" => ranges,
    "y_lines" => length(lines),
    "cross_section_nodes" => n_cs,
    "cross_section_triangles" => length(triangles),
    "stations" => nx + 1,
    "nodes" => length(nodes),
    "tetrahedra" => length(tetrahedra),
    "surface_triangles" => length(surface_elements),
    "surface_attribute_counts" => surface_counts,
    "surface_area_um2" => surface_areas,
    "volume_attribute_counts" => volume_counts,
    "volume_um3" => volumes,
    "sha256" => bytes2hex(open(sha256, output)),
    "bytes" => filesize(output)
)
open(replace(output, r"\.msh2$" => ".json"), "w") do stream
    JSON.print(stream, manifest, 2)
    return println(stream)
end
println(
    "Saved ",
    output,
    ": nodes ",
    length(nodes),
    ", tetrahedra ",
    length(tetrahedra),
    ", surface triangles ",
    length(surface_elements)
)
println(
    "Surface counts: ",
    surface_counts,
    "\nSurface areas (um^2): ",
    surface_areas,
    "\nVolumes (um^3): ",
    volumes
)
