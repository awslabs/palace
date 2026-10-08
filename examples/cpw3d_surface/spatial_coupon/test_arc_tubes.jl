# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
using SHA
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# Arc metal sides in the prism-tube recipe - T2 (block (b) DESIGN sections 1-2, AMENDMENT 1
# A1 (4), A2, A3, A7 MINOR-6 / MINOR-7; supervisor decision 303): the tagged arc runs of the
# plan-view boundary, the revolved ArcTube split into <= 90-degree parts of equal angular
# fraction, the smooth joints sharing one cross-section owned by the arc, the exact revolved
# centroids, the chord-count independence of the CAD and the census.

guard_message(f) =
    try
        f()
        ""
    catch e
        e isa ErrorException ? e.msg : sprint(showerror, e)
    end

# A synthetic L-shaped metal strip of width 1 (R 0.5): the horizontal leg enters through the
# x0 face, bends by `sweep_degrees` about the centre C = (0, 1) along its OUTER edge (a convex
# arc of radius 1, chorded every `chord_degrees`) and leaves through the y1 face (sweep 90) or
# continues straight after the bend; the inner corner at C is a sharp concave 90-degree corner
# (a semantic corner). Both arc joints are exactly tangent (smooth). The signature rows are
# the straight edges and one chord row per chord; the boundary carries the arc tags of
# generate_spatial_response.ARC_BOUNDARY_COLUMNS; the box is pinned by a process library.
# Round 3 class (11) (fix 11; part M 1.5 S1 / S3): `split_chords` = the chord indices after which
# a NEW ArcId starts (the arc serialised as several same-circle entries meeting at smooth joints,
# turn 0, JointSmooth 1), `ids` the ArcId of each member in loop order (default 1, 2, ...: a
# permutation exercises the owner rule "the earlier tube in install order", not the smaller id).
function write_strip_inputs(
    directory;
    chord_degrees=5.0,
    plane=0.0,
    kink_degrees=0.0,
    split_chords=Int[],
    ids=nothing
)
    centre = (0.0, 1.0)
    rho = 1.0
    sweep = 90.0
    n = ceil(Int, sweep / chord_degrees - 1.0e-9)
    kink = (sind(kink_degrees), cosd(kink_degrees))
    # The signature rows first (the straight edges as claims of length 2 reaching 2 R past
    # their free ends under the box rule, one row per chord); the box follows from them.
    rows = NamedTuple[]
    function straight!(a, b, gap_sign)
        d = (b[1] - a[1], b[2] - a[2])
        L = hypot(d...)
        t = (d[1] / L, d[2] / L)
        return push!(
            rows,
            (
                point=(0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), plane),
                tangent=(t[1], t[2], 0.0),
                gap=(gap_sign * t[2], -gap_sign * t[1], 0.0),
                interval=(-0.5 * L, 0.5 * L),
                normal_sign=1.0,
                vertex_arm=false,
                slot=0,
                conductor=1
            )
        )
    end
    chord_points = [
        (
            centre[1] + rho * cosd(-90.0 + sweep * k / n),
            centre[2] + rho * sind(-90.0 + sweep * k / n)
        ) for k = 0:n
    ]
    # The gap points away from the metal: the metal lies to the LEFT of the counterclockwise
    # loop, so gap = (t_y, -t_x) (to the right) for every side of this loop.
    straight!((-1.0, 0.0), (0.0, 0.0), 1.0)
    for k = 1:n
        straight!(chord_points[k], chord_points[k + 1], 1.0)
    end
    straight!((1.0, 1.0), (1.0 + kink[1] / kink[2], 2.0), 1.0)
    straight!((0.0, 2.0), (0.0, 1.0), 1.0)
    straight!((0.0, 1.0), (-1.0, 1.0), 1.0)
    lower, upper = row_coupon_bounds(rows, 0.5, 0.1, 0.05)
    chord_points = [
        (
            centre[1] + rho * cosd(-90.0 + sweep * k / n),
            centre[2] + rho * sind(-90.0 + sweep * k / n)
        ) for k = 0:n
    ]
    # The arc runs from (0, 0) (angle -90) to (1, 1) (angle 0) counterclockwise; a kinked
    # variant rotates the following vertical edge by kink_degrees about (1, 1).
    top_exit = (1.0 + kink[1] / kink[2] * (upper[2] - 1.0), upper[2])
    polygon = Tuple{Float64, Float64}[]
    push!(polygon, (lower[1], 0.0))
    append!(polygon, chord_points)          # (0, 0) ... (1, 1)
    push!(polygon, top_exit)
    push!(polygon, (0.0, upper[2]))
    push!(polygon, (0.0, 1.0))             # the concave corner at the centre
    push!(polygon, (lower[1], 1.0))
    m = length(polygon)
    on_face(p, q) = any(
        (abs(p[d] - lower[d]) <= 1.0e-9 && abs(q[d] - lower[d]) <= 1.0e-9) ||
            (abs(p[d] - upper[d]) <= 1.0e-9 && abs(q[d] - upper[d]) <= 1.0e-9) for
        d = 1:2
    )
    classes =
        [on_face(polygon[i], polygon[i % m + 1]) ? "Continuation" : "Physical" for i = 1:m]
    # Arc tags on the chord rows (vertex 2 .. n + 1 of the polygon); joint tags at both ends.
    all(1 <= k < n for k in split_chords) || error("split_chords must lie inside the arc")
    members = length(split_chords) + 1
    ids = ids === nothing ? collect(1:members) : collect(ids)
    length(ids) == members || error("one id per arc member")
    arcs = Vector{Any}(nothing, m)
    for i = 2:(n + 1)
        member = 1 + count(k -> i - 1 > k, split_chords)
        arcs[i] = (ids[member], centre[1], centre[2], rho, 1)
    end
    joints = Vector{Any}(nothing, m)
    joints[2] = (0.0, 1)
    for k in split_chords
        joints[k + 2] = (0.0, 1)                # the same-circle joint of two members
    end
    joints[n + 2] = (deg2rad(kink_degrees), kink_degrees <= 1.0e-4 * 180 / pi ? 1 : 0)
    # Corners (decision 320 convention): the concave corner; the perpendicular box exits of
    # Physical class ((lower, 0) and (0, upper)) keep the legacy corner; a kinked arc end.
    corners = [(0.0, 1.0, plane), (lower[1], 0.0, plane), (0.0, upper[2], plane)]
    kink_degrees > 1.0e-4 * 180 / pi && push!(corners, (1.0, 1.0, plane))
    open(joinpath(directory, "signature.csv"), "w") do io
        println(io, "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm")
        for (i, row) in enumerate(rows)
            println(
                io,
                join(
                    [
                        i,
                        0,
                        1,
                        row.point...,
                        row.gap...,
                        row.tangent...,
                        1,
                        row.interval...,
                        0
                    ],
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "boundary.csv"), "w") do io
        println(
            io,
            "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y,ArcId,ArcCx,ArcCy,ArcR,ArcSign,JointTurn,JointSmooth"
        )
        for (i, point) in enumerate(polygon)
            arc = arcs[i] === nothing ? ["", "", "", "", ""] : collect(arcs[i])
            joint = joints[i] === nothing ? ["", ""] : collect(joints[i])
            println(
                io,
                join(
                    vcat([1, i, 1, plane, 0, classes[i], point[1], point[2]], arc, joint),
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for point in polygon
            println(io, join([1, 1, plane, point[1], point[2]], ","))
        end
    end
    # Design round 2 F5-A: a kinked arc end is an INVARIANT corner (the arc's end tangent
    # against the kinked side: a non-zero dot product); the contract records it as
    # derive_semantic_contract does.
    contract = Dict{String, Any}(
        "Version" => 1,
        "SemanticCorners" => [[c[1], c[2], c[3]] for c in corners]
    )
    kink_degrees > 1.0e-4 * 180 / pi && (
        contract["Derivation"] = Dict{String, Any}(
            "InvariantCorners" => Dict{String, Any}(
                "Rule" => "test fixture: the kinked arc / line joint",
                "Points" => [[1.0, 1.0, plane]]
            )
        )
    )
    open(joinpath(directory, "semantic.json"), "w") do io
        return write_json(io, contract)
    end
    return (
        signature=joinpath(directory, "signature.csv"),
        boundary=joinpath(directory, "boundary.csv"),
        mask=joinpath(directory, "mask.csv"),
        semantic=joinpath(directory, "semantic.json"),
        polygon=polygon,
        classes=classes,
        corners=corners,
        chords=n,
        lower=lower,
        upper=upper,
        centre=centre,
        rho=rho
    )
end

function build_strip_coupon(
    directory;
    fabricated=true,
    stem="strip",
    chord_degrees=5.0,
    kink_degrees=0.0,
    labels_only=false,
    split_chords=Int[],
    ids=nothing
)
    inputs = write_strip_inputs(
        directory;
        chord_degrees=chord_degrees,
        kink_degrees=kink_degrees,
        split_chords=split_chords,
        ids=ids
    )
    mesh = joinpath(directory, "coupon-$stem.msh")
    census = joinpath(directory, "census-$stem.json")
    generate_spatial_coupon(;
        signature=inputs.signature,
        mask=inputs.mask,
        boundary=inputs.boundary,
        fabricated=fabricated,
        filename=mesh,
        radius=0.5,
        metal_thickness=0.1,
        overetch=0.05,
        sidewall_angle=90.0,
        top_rounding=0.0,
        trench_rounding=0.0,
        lc_fine=0.05,
        lc_tangent=0.1,
        lc_far=0.3,
        max_nodes=3_000_000,
        max_elements=3_000_000,
        semantic_contract=inputs.semantic,
        corner_isotropy_radius=0.1,
        corner_census=census,
        edge_size=0.01,
        edge_growth_ratio=2.0,
        corner_size=0.01,
        prism_tubes=true,
        far_growth=0.5,
        maximum_corner_aspect=4.0,
        minimum_scaled_jacobian=0.01,
        maximum_jacobian_condition=1000.0,
        quality_displacement_over_normal=0.75,
        # The invariant-corner verdict bound (design round 2 F5-A): a kinked arc end is invariant.
        corner_shape_gate=5.0,
        labels_only=labels_only ? census : nothing
    )
    return parse_json(read(census, String)), mesh, inputs
end

# Round 2b (decision 437 (3); CC DESIGN A2 (7) / AMENDMENT 2 B3): a metal block whose bottom edge
# runs from the convex corner (-1, 0) to (0, 0) and continues as a convex arc (centre (0, rho),
# radius rho = x1 / sin theta) that is CUT by the x1 box face at the tilt theta (the angle between
# the arc tangent at the cut and the face normal); the left edge x = -1 leaves the top face
# perpendicularly (a legacy box corner). The box follows from the rows (the bottom edge's claim
# ending at (0, 0) sets x1 = 0 + 2 R + R = 1.5; the left edge's claim sets y1 = 3.5, y0 = -1.2; the chord
# rows stop at x <= x1 - R so they never move the box): the arc's end vertex lies exactly on the
# face. The arc end is a box-face cut end (no corner, class Continuation), the joint at (0, 0)
# exactly tangent (turn 0, smooth).
# Round 3 (DESIGN-part-M 3.3 Fact 1): the default tie rho = x1 / sin theta makes rho (1 - sin theta)
# shrink below the tube envelope at steep tilts (a fixture artefact, not a mesher property), so
# `rho` is a parameter. Given, the arc still starts at (0, 0) tangent to the bottom edge and the
# face moves to x1 = rho sin theta (> 1.5 required); since the bottom edge's extended claim can
# only set x1 = 1.5, a NOTCH in the block's top-right sets the face instead: a horizontal metal
# edge at y_t = y_face + 8 R from (x1 - 3, y_t) to the face (a legacy perpendicular box end) whose
# claim ends at x1 - 1.5 (2 R extension + R padding = x1), closed by a vertical edge up to the top
# face y1 (the left edge's claim is lengthened so that y1 = y_t + 3). The arc face end, 8 R below
# the notch on the same face, keeps its round-2b configuration. The default (rho === nothing)
# writes the round-2b fixture byte for byte.
function write_arc_face_end_inputs(
    directory;
    theta_degrees=45.0,
    chord_degrees=5.0,
    plane=0.0,
    rho=nothing
)
    radius = 0.5
    rows = NamedTuple[]
    function straight!(a, b, gap_sign)
        d = (b[1] - a[1], b[2] - a[2])
        L = hypot(d...)
        t = (d[1] / L, d[2] / L)
        return push!(
            rows,
            (
                point=(0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), plane),
                tangent=(t[1], t[2], 0.0),
                gap=(gap_sign * t[2], -gap_sign * t[1], 0.0),
                interval=(-0.5 * L, 0.5 * L),
                normal_sign=1.0,
                vertex_arm=false,
                slot=0,
                conductor=1
            )
        )
    end
    straight!((-1.0, 0.0), (0.0, 0.0), 1.0)
    notch = rho !== nothing
    if notch
        x1 = rho * sind(theta_degrees)
        x1 > 3.0 * radius + 1.0e-9 || error(
            "rho $rho at $theta_degrees degrees puts the face at x1 = $x1 <= 3 R: the bottom edge sets the box"
        )
        y_face = rho * (1.0 - cosd(theta_degrees))
        y_notch = y_face + 8.0 * radius
        x_notch = x1 - 3.0
    else
        x1 = 0.0 + 3.0 * radius
        rho = x1 / sind(theta_degrees)
        y_notch = 0.0
        x_notch = 0.0
    end
    # The left edge's claim ends 0.3 above the corner so that the box face y0 (-1.2) does not
    # coincide with the 3 R collar of the bottom edge and of the arc (y = -1.5: a tangent collar /
    # face contact leaves sliver tetrahedra between them).
    straight!((-1.0, notch ? y_notch + 1.5 : 2.0), (-1.0, 0.3), 1.0)
    if notch
        # The notch's horizontal edge (claim ending at x1 - 1.5: it sets the face) and its vertical
        # edge (claim 1 long: its extension stays below y1 = y_notch + 3).
        straight!((x_notch, y_notch), (x1 - 1.5, y_notch), 1.0)
        straight!((x_notch, y_notch + 1.0), (x_notch, y_notch), 1.0)
    end
    centre = (0.0, rho)
    sweep = theta_degrees
    n = ceil(Int, sweep / chord_degrees - 1.0e-9)
    chord_points = [
        (
            centre[1] + rho * cosd(-90.0 + sweep * k / n),
            centre[2] + rho * sind(-90.0 + sweep * k / n)
        ) for k = 0:n
    ]
    chord_points[end] = (x1, rho * (1.0 - cosd(theta_degrees)))   # exactly on the face
    # Chord rows while the row's transverse pad (+- R along its gap) and the box pad stay left
    # of x1: max(x) + R |gap_x| + R <= x1. With `rho` given (round 3 9H: the steep representative
    # fixture, rho (1 - sin theta) >> the tube envelope) the LAST chord row carries the chain's
    # outer end and is extended by 2 R along its tangent (extended_interval), so its longitudinal
    # reach 2 R |t_x| joins the test; the default fixture keeps its round-2b rows byte for byte.
    for k = 1:n
        p, q = chord_points[k], chord_points[k + 1]
        t = (q[1] - p[1], q[2] - p[2]) ./ hypot(q[1] - p[1], q[2] - p[2])
        reach = radius * abs(t[2]) + radius + (notch ? 2.0 * radius * abs(t[1]) : 0.0)
        max(p[1], q[1]) + reach <= x1 + 1.0e-12 || break
        straight!(p, q, 1.0)
    end
    lower, upper = row_coupon_bounds(rows, radius, 0.1, 0.05)
    @assert upper[1] == x1 "the box x1 $(upper[1]) is not the designed $x1"
    polygon = Tuple{Float64, Float64}[(-1.0, 0.0)]
    append!(polygon, chord_points)                  # (0, 0) ... the face point
    if notch
        push!(polygon, (x1, y_notch))
        push!(polygon, (x_notch, y_notch))
        push!(polygon, (x_notch, upper[2]))
    else
        push!(polygon, (x1, upper[2]))
    end
    push!(polygon, (-1.0, upper[2]))
    m = length(polygon)
    on_face(p, q) = any(
        (abs(p[d] - lower[d]) <= 1.0e-9 && abs(q[d] - lower[d]) <= 1.0e-9) ||
            (abs(p[d] - upper[d]) <= 1.0e-9 && abs(q[d] - upper[d]) <= 1.0e-9) for
        d = 1:2
    )
    classes =
        [on_face(polygon[i], polygon[i % m + 1]) ? "Continuation" : "Physical" for i = 1:m]
    arcs = Vector{Any}(nothing, m)
    for i = 2:(n + 1)
        arcs[i] = (1, centre[1], centre[2], rho, 1)
    end
    joints = Vector{Any}(nothing, m)
    joints[2] = (0.0, 1)
    # Corners: the convex 90-degree corner (-1, 0) and the legacy perpendicular box corner (-1, y1);
    # the arc end on the face is a cut end (decision 320 for arcs: no corner). With the notch: its
    # legacy perpendicular box end (x1, y_notch), its convex corner (x_notch, y_notch) and its
    # legacy box corner (x_notch, y1).
    corners = [(-1.0, 0.0, plane), (-1.0, upper[2], plane)]
    notch && append!(
        corners,
        [(x1, y_notch, plane), (x_notch, y_notch, plane), (x_notch, upper[2], plane)]
    )
    open(joinpath(directory, "signature.csv"), "w") do io
        println(io, "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm")
        for (i, row) in enumerate(rows)
            println(
                io,
                join(
                    [
                        i,
                        0,
                        1,
                        row.point...,
                        row.gap...,
                        row.tangent...,
                        1,
                        row.interval...,
                        0
                    ],
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "boundary.csv"), "w") do io
        println(
            io,
            "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y,ArcId,ArcCx,ArcCy,ArcR,ArcSign,JointTurn,JointSmooth"
        )
        for (i, point) in enumerate(polygon)
            arc = arcs[i] === nothing ? ["", "", "", "", ""] : collect(arcs[i])
            joint = joints[i] === nothing ? ["", ""] : collect(joints[i])
            println(
                io,
                join(
                    vcat([1, i, 1, plane, 0, classes[i], point[1], point[2]], arc, joint),
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for point in polygon
            println(io, join([1, 1, plane, point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "semantic.json"), "w") do io
        return write_json(
            io,
            Dict{String, Any}(
                "Version" => 1,
                "SemanticCorners" => [[c[1], c[2], c[3]] for c in corners]
            )
        )
    end
    return (
        signature=joinpath(directory, "signature.csv"),
        boundary=joinpath(directory, "boundary.csv"),
        mask=joinpath(directory, "mask.csv"),
        semantic=joinpath(directory, "semantic.json"),
        polygon=polygon,
        classes=classes,
        corners=corners,
        chords=n,
        lower=lower,
        upper=upper,
        centre=centre,
        rho=rho,
        face_point=chord_points[end]
    )
end

# Mesh nodes within `reach` of `point`.
function nodes_near(path, point, reach)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(path)
    _, coordinates, _ = gmsh.model.mesh.getNodes()
    gmsh.finalize()
    xyz = reshape(coordinates, 3, :)
    return [xyz[:, i] for i = 1:size(xyz, 2) if norm(xyz[:, i] .- point) <= reach]
end

# The full build of a fixture written by `writer` (write_strip_inputs / write_arc_face_end_inputs)
# at the test sizes (R 0.5, EdgeSize 0.01: 2-ring tubes) with the production gates; `edge_size`
# 0.00025 / 0.002 gives the production fabricated / thin tubes on the same coupon.
# Round 3 class (7) (G.7.6 S-G7-a / S-G7-b; the C3 trio's and 0c94ec951e10's 3-chord claim arcs):
# a vertical metal lead of width 1 (R 0.5) crossing the box from the y0 to the y1 face whose RIGHT
# side carries a SHORT convex arc of radius `rho` = 52 (104 R) between (0.5, -0.5) and its end,
# `chords` chords of equal angle over `sweep` (default 3 x 0.25 R chords: 7.2e-3 rad), both joints
# exactly tangent (smooth): the lower straight side is vertical, the upper one continues at the
# arc's end tangent (tilted inward by `sweep`) to the y1 face - a straight box-face CUT end. The
# left side is straight; its ends and the right side's bottom end are legacy perpendicular box
# corners (decision 320). Rows: the straight sides as claims of length 3 (the box follows:
# [-1.5, 1.5 + ...] x [-3, 3]) and one row per chord; the boundary carries the arc tags.
# `convex` false: the CONCAVE twin (sigma -1; the metal outside the circle, centre to the right,
# the side turning right and the upper side tilted outward) - S-G7-b's shape, whose two straight
# neighbours never read each other as facing across the metal (F6) around the tiny arc.
function write_short_arc_lead_inputs(
    directory;
    chords=3,
    sweep=3 * 0.125 / 52.0,
    rho=52.0,
    plane=0.0,
    convex=true
)
    sign = convex ? 1 : -1
    centre = (0.5 - sign * rho, -0.5)
    angle(k) = convex ? sweep * k / chords : pi - sweep * k / chords
    chord_points = [
        (centre[1] + rho * cos(angle(k)), centre[2] + rho * sin(angle(k))) for k = 0:chords
    ]
    chord_points[1] = (0.5, -0.5)
    tilt = (-sign * sin(sweep), cos(sweep))             # the arc's end tangent (travel +y)
    rows = NamedTuple[]
    function straight!(a, b)
        d = (b[1] - a[1], b[2] - a[2])
        L = hypot(d...)
        t = (d[1] / L, d[2] / L)
        return push!(
            rows,
            (
                point=(0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), plane),
                tangent=(t[1], t[2], 0.0),
                gap=(t[2], -t[1], 0.0),
                interval=(-0.5 * L, 0.5 * L),
                normal_sign=1.0,
                vertex_arm=false,
                slot=0,
                conductor=1
            )
        )
    end
    straight!((0.5, -1.5), (0.5, -0.5))                 # the right side below the arc
    for k = 1:chords
        straight!(chord_points[k], chord_points[k + 1])
    end
    p_end = chord_points[end]
    upper_row_end = (p_end[1] + tilt[1] * (1.5 - p_end[2]) / tilt[2], 1.5)
    straight!(p_end, upper_row_end)                     # the right side above the arc, tilted
    straight!((-0.5, 1.5), (-0.5, -1.5))                # the left side (travelling down: metal to its left)
    lower, upper = row_coupon_bounds(rows, 0.5, 0.1, 0.05)
    top_right = (p_end[1] + tilt[1] * (upper[2] - p_end[2]) / tilt[2], upper[2])
    polygon = vcat(
        [(-0.5, lower[2]), (0.5, lower[2])],
        chord_points,
        [top_right, (-0.5, upper[2])]
    )
    m = length(polygon)
    on_face(p, q) = any(
        (abs(p[d] - lower[d]) <= 1.0e-9 && abs(q[d] - lower[d]) <= 1.0e-9) ||
            (abs(p[d] - upper[d]) <= 1.0e-9 && abs(q[d] - upper[d]) <= 1.0e-9) for
        d = 1:2
    )
    classes =
        [on_face(polygon[i], polygon[i % m + 1]) ? "Continuation" : "Physical" for i = 1:m]
    arcs = Vector{Any}(nothing, m)
    for i = 3:(2 + chords)                              # the chord sides: vertices 3 .. 2 + chords
        arcs[i] = (1, centre[1], centre[2], rho, sign)
    end
    joints = Vector{Any}(nothing, m)
    joints[3] = (0.0, 1)
    joints[3 + chords] = (0.0, 1)
    corners = [(-0.5, lower[2], plane), (0.5, lower[2], plane), (-0.5, upper[2], plane)]
    open(joinpath(directory, "signature.csv"), "w") do io
        println(io, "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm")
        for (i, row) in enumerate(rows)
            println(
                io,
                join(
                    [
                        i,
                        0,
                        1,
                        row.point...,
                        row.gap...,
                        row.tangent...,
                        1,
                        row.interval...,
                        0
                    ],
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "boundary.csv"), "w") do io
        println(
            io,
            "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y,ArcId,ArcCx,ArcCy,ArcR,ArcSign,JointTurn,JointSmooth"
        )
        for (i, point) in enumerate(polygon)
            arc = arcs[i] === nothing ? ["", "", "", "", ""] : collect(arcs[i])
            joint = joints[i] === nothing ? ["", ""] : collect(joints[i])
            println(
                io,
                join(
                    vcat([1, i, 1, plane, 0, classes[i], point[1], point[2]], arc, joint),
                    ","
                )
            )
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for point in polygon
            println(io, join([1, 1, plane, point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "semantic.json"), "w") do io
        return write_json(
            io,
            Dict{String, Any}(
                "Version" => 1,
                "SemanticCorners" => [[c[1], c[2], c[3]] for c in corners]
            )
        )
    end
    return (
        signature=joinpath(directory, "signature.csv"),
        boundary=joinpath(directory, "boundary.csv"),
        mask=joinpath(directory, "mask.csv"),
        semantic=joinpath(directory, "semantic.json"),
        polygon=polygon,
        classes=classes,
        corners=corners,
        chords=chords,
        lower=lower,
        upper=upper,
        centre=centre,
        rho=rho
    )
end

function build_arc_coupon(
    directory,
    writer;
    fabricated=true,
    stem="arc",
    labels_only=false,
    edge_size=0.01,
    lc_fine=0.05,
    lc_tangent=0.1,
    lc_far=0.3,
    max_elements=3_000_000,
    writer_options...
)
    inputs = writer(directory; writer_options...)
    mesh = joinpath(directory, "coupon-$stem.msh")
    census = joinpath(directory, "census-$stem.json")
    generate_spatial_coupon(;
        signature=inputs.signature,
        mask=inputs.mask,
        boundary=inputs.boundary,
        fabricated=fabricated,
        filename=mesh,
        radius=0.5,
        metal_thickness=0.1,
        overetch=0.05,
        sidewall_angle=90.0,
        top_rounding=0.0,
        trench_rounding=0.0,
        lc_fine=lc_fine,
        lc_tangent=lc_tangent,
        lc_far=lc_far,
        max_nodes=max_elements,
        max_elements=max_elements,
        semantic_contract=inputs.semantic,
        corner_isotropy_radius=0.1,
        corner_census=census,
        edge_size=edge_size,
        edge_growth_ratio=2.0,
        corner_size=edge_size,
        prism_tubes=true,
        far_growth=0.5,
        maximum_corner_aspect=4.0,
        minimum_scaled_jacobian=0.01,
        maximum_jacobian_condition=1000.0,
        quality_displacement_over_normal=0.75,
        corner_shape_gate=5.0,
        labels_only=labels_only ? census : nothing
    )
    return parse_json(read(census, String)), mesh, inputs
end

# Round 3 class (1) (decisions 491 / 510; DESIGN-part-G G.1.6): the CORNER-FREE synthetics.
# S-G1-a: a metal band of width 0.5 (R 0.5) entering through the x0 face and leaving through
# the y1 face with ONE exactly tangent 90-degree bend: its outer edge a convex arc of radius
# 0.75 (1.5 R) and its inner edge the concentric CONCAVE arc of radius 0.25 (0.5 R, below the
# 3 R collar: the fabricated collar offset collapses it) about C = (0, 0.75); every joint is
# tangent (turn 0.0, smooth), every straight end is a box-face cut end: 0 semantic corners.
# S-G1-b (510 MINOR-2; 24be6e90570f's twin): the same band plus a ground region beyond a 0.25
# gap whose edge follows the band with a second CONCAVE arc of radius 1.0 (2 R, below the
# collar too): two collapsed concave arcs in one collar union around the convex one.
function write_bent_band_inputs(
    directory;
    ground=false,
    chord_degrees=5.0,
    plane=0.0,
    tilt_degrees=2.5
)
    radius = 0.5
    centre = (0.0, 0.75)
    inner, outer, ground_radius = 0.25, 0.75, 1.0
    sweep = 90.0
    n = ceil(Int, sweep / chord_degrees - 1.0e-9)
    # The whole outline is rotated by tilt_degrees about the arc centre, so that every
    # straight leg crosses its box face obliquely: a box-face CUT END (decision 320), not the
    # exactly perpendicular legacy corner (24be6e90570f's legs cross at ~2.5 degrees).
    rotate(p) = (
        centre[1] + cosd(tilt_degrees) * (p[1] - centre[1]) -
        sind(tilt_degrees) * (p[2] - centre[2]),
        centre[2] +
        sind(tilt_degrees) * (p[1] - centre[1]) +
        cosd(tilt_degrees) * (p[2] - centre[2])
    )
    west = (-cosd(tilt_degrees), -sind(tilt_degrees))      # the rotated (-1, 0): the x0-face legs
    north = (-sind(tilt_degrees), cosd(tilt_degrees))      # the rotated (0, 1): the y1-face legs
    rows = NamedTuple[]
    function straight!(a, b)
        d = (b[1] - a[1], b[2] - a[2])
        L = hypot(d...)
        t = (d[1] / L, d[2] / L)
        # The metal lies to the LEFT of every side (counterclockwise loops), the gap to the right.
        return push!(
            rows,
            (
                point=(0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), plane),
                tangent=(t[1], t[2], 0.0),
                gap=(t[2], -t[1], 0.0),
                interval=(-0.5 * L, 0.5 * L),
                normal_sign=1.0,
                vertex_arm=false,
                slot=0,
                conductor=1
            )
        )
    end
    # Chord points of a 90-degree arc about the centre from angle -90 (below C) to 0 (right of
    # C), rotated by the tilt.
    arc_points(rho) = [
        rotate((
            centre[1] + rho * cosd(-90.0 + sweep * k / n),
            centre[2] + rho * sind(-90.0 + sweep * k / n)
        )) for k = 0:n
    ]
    outer_points, inner_points, ground_points =
        arc_points(outer), reverse(arc_points(inner)), reverse(arc_points(ground_radius))
    along(p, direction, s) = (p[1] + s * direction[1], p[2] + s * direction[2])
    # The band (loop 1): outer edge east along y = 0, the convex arc counterclockwise, north
    # along x = 0.75; inner edge south along x = 0.25, the concave arc clockwise, west along
    # y = 0.5 (before the tilt); every straight leg's claim is 1.0 = 2 R long.
    straight!(along(outer_points[1], west, 1.0), outer_points[1])
    for k = 1:n
        straight!(outer_points[k], outer_points[k + 1])
    end
    straight!(outer_points[end], along(outer_points[end], north, 1.0))
    straight!(along(inner_points[1], north, 1.0), inner_points[1])
    for k = 1:n
        straight!(inner_points[k], inner_points[k + 1])
    end
    straight!(inner_points[end], along(inner_points[end], west, 1.0))
    if ground
        # The ground (loop 2): south along x = 1.0, the concave arc clockwise, west along y = -0.25.
        straight!(along(ground_points[1], north, 1.0), ground_points[1])
        for k = 1:n
            straight!(ground_points[k], ground_points[k + 1])
        end
        straight!(ground_points[end], along(ground_points[end], west, 1.0))
    end
    lower, upper = row_coupon_bounds(rows, radius, 0.1, 0.05)
    # The legs meet their box faces where their lines cross x = x0 (west legs) / y = y1 (north legs).
    x0_face(p) = (lower[1], p[2] + (lower[1] - p[1]) * west[2] / west[1])
    y1_face(p) = (p[1] + (upper[2] - p[2]) * north[1] / north[2], upper[2])
    loops = Vector{Tuple{Float64, Float64}}[]
    band = Tuple{Float64, Float64}[x0_face(outer_points[1])]
    append!(band, outer_points)                          # (0, 0) ... (0.75, 0.75) before the tilt
    push!(band, y1_face(outer_points[end]))
    push!(band, y1_face(inner_points[1]))
    append!(band, inner_points)                          # (0.25, 0.75) ... (0, 0.5) before the tilt
    push!(band, x0_face(inner_points[end]))
    push!(loops, band)
    if ground
        region = Tuple{Float64, Float64}[y1_face(ground_points[1])]
        append!(region, ground_points)                   # (1, 0.75) ... (0, -0.25) before the tilt
        push!(region, x0_face(ground_points[end]))
        push!(region, (lower[1], lower[2]))
        push!(region, (upper[1], lower[2]))
        push!(region, (upper[1], upper[2]))
        push!(loops, region)
    end
    on_face(p, q) = any(
        (abs(p[d] - lower[d]) <= 1.0e-9 && abs(q[d] - lower[d]) <= 1.0e-9) ||
            (abs(p[d] - upper[d]) <= 1.0e-9 && abs(q[d] - upper[d]) <= 1.0e-9) for
        d = 1:2
    )
    # Arc tags (ArcSign +1: the dielectric outside the circle) on the chord rows, joint tags
    # (turn 0.0, smooth) at both ends of every arc; the band's arcs are 1 (outer) / 2 (inner),
    # the ground's arc 3.
    arcs_of(points, first, rho, sign, id) =
        Dict(first + k - 1 => (id, centre[1], centre[2], rho, sign) for k = 1:n)
    tags = [arcs_of(outer_points, 2, outer, 1, 1)]
    merge!(tags[1], arcs_of(inner_points, n + 5, inner, -1, 2))
    joints = [Set([2, n + 2, n + 5, 2n + 5])]
    if ground
        push!(tags, arcs_of(ground_points, 2, ground_radius, -1, 3))
        push!(joints, Set([2, n + 2]))
    end
    open(joinpath(directory, "signature.csv"), "w") do io
        println(io, "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm")
        for (i, row) in enumerate(rows)
            println(
                io,
                join(
                    [
                        i,
                        row.slot,
                        row.conductor,
                        row.point...,
                        row.gap...,
                        row.tangent...,
                        1,
                        row.interval...,
                        0
                    ],
                    ","
                )
            )
        end
    end
    classes = Vector{String}[]
    open(joinpath(directory, "boundary.csv"), "w") do io
        println(
            io,
            "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y,ArcId,ArcCx,ArcCy,ArcR,ArcSign,JointTurn,JointSmooth"
        )
        for (l, polygon) in enumerate(loops)
            m = length(polygon)
            loop_classes = [
                on_face(polygon[i], polygon[i % m + 1]) ? "Continuation" : "Physical"
                for i = 1:m
            ]
            push!(classes, loop_classes)
            for (i, point) in enumerate(polygon)
                arc = haskey(tags[l], i) ? collect(tags[l][i]) : ["", "", "", "", ""]
                joint = i in joints[l] ? [0.0, 1] : ["", ""]
                println(
                    io,
                    join(
                        vcat(
                            [l, i, 1, plane, 0, loop_classes[i], point[1], point[2]],
                            arc,
                            joint
                        ),
                        ","
                    )
                )
            end
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for (l, polygon) in enumerate(loops), point in polygon
            println(io, join([l, 1, plane, point[1], point[2]], ","))
        end
    end
    arc_count = ground ? 3 : 2
    contract = Dict{String, Any}(
        "Version" => 1,
        "SemanticCorners" => [],
        "Derivation" => Dict{String, Any}(
            "CornerFree" => Dict{String, Any}(
                "Rule" => CORNER_FREE_RULE,
                "ArcInteriorVertices" => arc_count * (n - 1),
                "SmoothJoints" => 2 * arc_count,
                "BoxFaceCutEnds" => ground ? 6 : 4
            )
        )
    )
    open(joinpath(directory, "semantic.json"), "w") do io
        return write_json(io, contract)
    end
    return (
        signature=joinpath(directory, "signature.csv"),
        boundary=joinpath(directory, "boundary.csv"),
        mask=joinpath(directory, "mask.csv"),
        semantic=joinpath(directory, "semantic.json"),
        loops=loops,
        classes=classes,
        chords=n,
        lower=lower,
        upper=upper,
        centre=centre,
        radii=(inner, outer, ground_radius),
        tilt_degrees=tilt_degrees
    )
end

function build_bent_band_coupon(
    directory;
    ground=false,
    fabricated=true,
    stem="band",
    labels_only=false,
    tilt_degrees=2.5
)
    inputs = write_bent_band_inputs(directory; ground=ground, tilt_degrees=tilt_degrees)
    mesh = joinpath(directory, "coupon-$stem.msh")
    census = joinpath(directory, "census-$stem.json")
    generate_spatial_coupon(;
        signature=inputs.signature,
        mask=inputs.mask,
        boundary=inputs.boundary,
        fabricated=fabricated,
        filename=mesh,
        radius=0.5,
        metal_thickness=0.1,
        overetch=0.05,
        sidewall_angle=90.0,
        top_rounding=0.0,
        trench_rounding=0.0,
        lc_fine=0.05,
        lc_tangent=0.1,
        lc_far=0.3,
        max_nodes=3_000_000,
        max_elements=3_000_000,
        semantic_contract=inputs.semantic,
        corner_isotropy_radius=0.1,
        corner_census=census,
        edge_size=0.01,
        edge_growth_ratio=2.0,
        corner_size=0.01,
        prism_tubes=true,
        far_growth=0.5,
        maximum_corner_aspect=4.0,
        minimum_scaled_jacobian=0.01,
        maximum_jacobian_condition=1000.0,
        quality_displacement_over_normal=0.75,
        corner_shape_gate=5.0,
        labels_only=labels_only ? census : nothing
    )
    return parse_json(read(census, String)), mesh, inputs
end

@testset "tagged arc runs: seeded by the tags, cross-checked against the fit, fail closed" begin
    mktempdir() do directory
        inputs = write_strip_inputs(directory; chord_degrees=5.0)
        loops = read_boundary(inputs.boundary)
        @test length(loops) == 1
        loop = loops[1]
        @test loop_has_arcs(loop)
        runs = tagged_arc_runs(loop, 1.0e-7 * 0.5)
        @test length(runs) == 1
        run = runs[1]
        @test run.id == 1 && run.sign == 1
        @test length(run.edge_indices) == inputs.chords
        @test isapprox(run.sweep, pi / 2; atol=1.0e-12)
        @test arc_part_count(run.sweep) == 1
        @test arc_part_count(deg2rad(90.000001)) == 2 &&
              arc_part_count(deg2rad(180.0)) == 2 &&
              arc_part_count(deg2rad(180.000001)) == 3
        # Equal angular fractions: a 90.000001-degree arc splits at 45.0000005, never a sliver.
        angles = arc_split_angles(
            [(1.0, 0.0), (cosd(90.000001), sind(90.000001))],
            [1, 2],
            (0.0, 0.0),
            deg2rad(90.000001)
        )
        @test length(angles) == 3 && isapprox(rad2deg(angles[2]), 45.0000005; atol=1.0e-9)
        # A legacy boundary (no arc columns) carries no tags and takes the untagged fit.
        legacy = joinpath(directory, "legacy.csv")
        open(legacy, "w") do io
            lines = readlines(inputs.boundary)
            for line in lines
                println(io, join(split(line, ",")[1:8], ","))
            end
        end
        legacy_loop = read_boundary(legacy)[1]
        @test !loop_has_arcs(legacy_loop)
        @test isempty(tagged_arc_runs(legacy_loop, 1.0e-7 * 0.5))
        # Round 3 R7 (decision 510 MINOR-6): the optional ArcChain column. The fixture (every stored
        # arc coupon) has none: chain 0, no "Chain" key in the census record; a boundary carrying it
        # reads the chain on every tag and run and records it.
        @test all(arc === nothing || arc.chain == 0 for arc in loop.arcs) && run.chain == 0
        records = metal_loop_records([loop], inputs.lower, inputs.upper, 1.0e-7 * 0.5)
        @test length(records) == 1 && !haskey(records[1]["Arcs"][1], "Chain")
        chained = joinpath(directory, "chained.csv")
        open(chained, "w") do io
            for (i, line) in enumerate(readlines(inputs.boundary))
                cells = split(line, ",")
                println(
                    io,
                    join(
                        vcat(cells, [i == 1 ? "ArcChain" : (cells[9] == "" ? "" : "2")]),
                        ","
                    )
                )
            end
        end
        chained_loop = read_boundary(chained)[1]
        @test all(arc === nothing || arc.chain == 2 for arc in chained_loop.arcs) &&
              count(arc !== nothing for arc in chained_loop.arcs) == inputs.chords
        chained_runs = tagged_arc_runs(chained_loop, 1.0e-7 * 0.5)
        @test length(chained_runs) == 1 &&
              chained_runs[1].chain == 2 &&
              chained_runs[1].point_indices == run.point_indices &&
              chained_runs[1].center == run.center
        records =
            metal_loop_records([chained_loop], inputs.lower, inputs.upper, 1.0e-7 * 0.5)
        @test records[1]["Arcs"][1]["Chain"] == 2
        # A corrupted tag (another centre) fails closed at the residual test (round 3 class (4),
        # G.4.3 Option A: the tagged circle IS the circle; the former three-point-fit comparison
        # is a printed diagnostic).
        broken = joinpath(directory, "broken.csv")
        open(broken, "w") do io
            for line in readlines(inputs.boundary)
                println(io, replace(line, ",1,0.0,1.0,1.0,1," => ",1,0.0,1.001,1.0,1,"))
            end
        end
        message =
            guard_message(() -> tagged_arc_runs(read_boundary(broken)[1], 1.0e-7 * 0.5))
        @test occursin("do not lie on the tagged circle", message)
        @test !occursin("disagrees with the tagged circle", message)
    end
    # Round 3 class (4) (DESIGN-part-G G.4.2 / G.4.5, decision 510): a 3-chord arc of 134 R on
    # the 1e-9 R grid (2bc3d927fda6's context arc 2: R_arc 255.02 um, chords 0.224 R, sagitta
    # 7.98e-4 um, conditioning R_arc / sagitta 3.2e5) is ACCEPTED - every vertex lies within one
    # grid quantum of the tagged circle - although the three-point fit through its quantised
    # vertices is off the tag by ~1e-4 um (> the 5.1e-5 um tolerance): the former check refused it.
    R = 1.9
    rho = 255.02107
    centre = (-246.39, -4.77)
    quantum = 1.0e-9 * R
    quantise(v) = round(v / quantum) * quantum
    theta0 = 0.3
    step = 0.224 * R / rho
    arc_points = [
        (
            quantise(centre[1] + rho * cos(theta0 + step * k)),
            quantise(centre[2] + rho * sin(theta0 + step * k))
        ) for k = 0:3
    ]
    short_points = vcat(
        arc_points,
        [
            (arc_points[end][1] - 3.0, arc_points[end][2] + 1.0),
            (arc_points[1][1] - 3.0, arc_points[1][2] - 1.0)
        ]
    )
    short_arcs = Vector{Union{Nothing, NamedTuple}}(nothing, 6)
    for i = 1:3
        short_arcs[i] = (id=2, centre=centre, radius=rho, sign=1)
    end
    short_loop = (
        conductor=1,
        plane=0.0,
        hole=false,
        points=short_points,
        classes=fill("Physical", 6),
        arcs=short_arcs,
        joints=Vector{Union{Nothing, NamedTuple}}(nothing, 6)
    )
    tolerance = 1.0e-9 * R
    fit_tolerance = arc_fit_tolerance(rho, tolerance)
    residual =
        maximum(abs(hypot(p[1] - centre[1], p[2] - centre[2]) - rho) for p in arc_points)
    @test residual <= quantum * (1.0 + 1.0e-6) && residual <= fit_tolerance
    fit = circle_through(arc_points[1], arc_points[2], arc_points[4], tolerance)
    @test fit !== nothing &&
          hypot(fit.center[1] - centre[1], fit.center[2] - centre[2]) > fit_tolerance
    short_runs = tagged_arc_runs(short_loop, tolerance)
    @test length(short_runs) == 1 &&
          short_runs[1].id == 2 &&
          length(short_runs[1].edge_indices) == 3 &&
          short_runs[1].center == centre &&
          short_runs[1].radius == rho
    diagnostic = arc_fit_diagnostic(
        short_points,
        [1, 2, 3, 4],
        short_arcs[1],
        fit_tolerance,
        tolerance
    )
    @test diagnostic !== nothing &&
          occursin("conditioning factor R_arc / sagitta", diagnostic) &&
          occursin("diagnostic only", diagnostic)
    # One vertex displaced off the tagged circle by three residual tolerances: refused by the
    # residual test (the guard that remains).
    displaced = copy(short_points)
    displaced[2] = (displaced[2][1] + 3.0 * fit_tolerance, displaced[2][2])
    displaced_loop = merge(short_loop, (points=displaced,))
    @test occursin(
        "do not lie on the tagged circle",
        guard_message(() -> tagged_arc_runs(displaced_loop, tolerance))
    )
    # A two-vertex run (one chord) on the circle is admitted; off it (3 tolerances) refused.
    one_chord_points = vcat(
        arc_points[1:2],
        [
            (arc_points[2][1] - 3.0, arc_points[2][2] + 1.0),
            (arc_points[1][1] - 3.0, arc_points[1][2] - 1.0)
        ]
    )
    one_chord_arcs = Vector{Union{Nothing, NamedTuple}}(nothing, 4)
    one_chord_arcs[1] = (id=2, centre=centre, radius=rho, sign=1)
    one_chord_loop = (
        conductor=1,
        plane=0.0,
        hole=false,
        points=one_chord_points,
        classes=fill("Physical", 4),
        arcs=one_chord_arcs,
        joints=Vector{Union{Nothing, NamedTuple}}(nothing, 4)
    )
    one_chord = tagged_arc_runs(one_chord_loop, tolerance)
    @test length(one_chord) == 1 && length(one_chord[1].edge_indices) == 1
    @test arc_fit_diagnostic(
        one_chord_points,
        [1, 2],
        one_chord_arcs[1],
        fit_tolerance,
        tolerance
    ) === nothing
    off_points = copy(one_chord_points)
    off_points[2] = (off_points[2][1] + 3.0 * fit_tolerance, off_points[2][2])
    @test occursin(
        "do not lie on the tagged circle",
        guard_message(
            () -> tagged_arc_runs(merge(one_chord_loop, (points=off_points,)), tolerance)
        )
    )
end

@testset "ArcTube geometry: frame, revolved centroids against OCC, face crossing" begin
    section = TubeSection(0.01, 2.0, 2, [-90.0 + 30.0 * j for j = 0:9], fill(2, 9))
    tube = ArcTube([0.0, 1.0, 0.1], 1.0, 1.0, [0.0, 0.0, 1.0], -pi / 2, 0.0, 0.5 * pi, 0.1)
    @test tube.orientation == -1.0         # sigma +1, b_z +1: e = n x b runs clockwise
    @test isapprox(tube_point(tube, 0.0, 0.0, 0.0), [0.0, 0.0, 0.1]; atol=1.0e-12)
    frame = arc_frame(tube, 0.0)
    @test isapprox(frame.n, [0.0, -1.0, 0.0]; atol=1.0e-12) &&
          isapprox(frame.e, [-1.0, 0.0, 0.0]; atol=1.0e-12)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("arc-tube")
    occ = gmsh.model.occ
    volumes = add_tube_volumes!(occ, tube, section)
    occ.synchronize()
    @test length(volumes) == 1
    volume = volumes[1][1][2]
    matched = match_tube_entities(volume, tube, section, volumes[1][2], 1.0e-9)
    @test length(matched) == length(tube_entities(tube, section, volumes[1][2]))
    gmsh.finalize()
    # A face end on the plane x = 0.3: the node circle of radius rho + u crosses it at the
    # arc length where its own circle meets the plane.
    face = FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        0.3,
        0.0,
        0.0,
        0.04,
        0.01,
        0.1;
        spacing_cap=Inf,
        face_axis=1,
        face_value=0.3
    )
    faced = ArcTube(
        [0.0, 1.0, 0.1],
        1.0,
        1.0,
        [0.0, 0.0, 1.0],
        -pi / 2,
        0.0,
        0.3,
        0.1;
        face_ends=[face]
    )
    for u in (0.0, 0.03, -0.02), w in (0.0, 0.02)
        s = arc_face_station(faced, face, u, w)
        point = tube_point(faced, u, w, s)
        @test isapprox(point[1], 0.3; atol=1.0e-12)
    end
    # Round 3 9H (part M 3.3 Fact 2): the arc's crossing slope max_u |s'(u)| = rho h / (r sqrt(r^2 -
    # h^2)) at the inner node circle r = rho - envelope. At the axis (envelope 0) it is tan theta
    # exactly (h = rho sin theta: the face plane at distance h from the centre is crossed at the
    # tilt theta); it grows with the envelope (the concave side shears more) and with h / rho, and
    # tends to tan theta as rho -> inf at fixed theta and envelope. It is the per-node derivative
    # of arc_face_station: a finite difference of the station across the inner node circle.
    for theta in (deg2rad(2.1), deg2rad(45.0), deg2rad(74.3)), rho in (1.56, 58.5)
        h = rho * sin(theta)
        @test isapprox(arc_crossing_slope(rho, 0.0, h), tan(theta); rtol=1.0e-12)
        slope = arc_crossing_slope(rho, 0.04, h)
        @test slope > tan(theta) && slope > arc_crossing_slope(rho, 0.02, h)
        @test isapprox(
            arc_crossing_slope(1.0e9 * rho, 0.04, 1.0e9 * h),
            tan(theta);
            rtol=1.0e-8
        )
    end
    @test isapprox(
        arc_crossing_slope(58.5, 0.04, 56.324169260983396),
        3.5997093489397574;
        rtol=1.0e-12
    )
    @test_throws ErrorException arc_crossing_slope(1.0, 0.04, 0.97)          # the node circle misses the face
    @test_throws ErrorException arc_crossing_slope(1.0, 1.0, 0.5)
    # the finite difference of the face station across the inner node circle (u = -envelope .. -envelope + du)
    inner = -0.04
    du = 1.0e-6
    s0 = arc_face_station(faced, face, inner, 0.0)
    s1 = arc_face_station(faced, face, inner + du, 0.0)
    @test isapprox(abs(s1 - s0) / du, arc_crossing_slope(1.0, 0.04, 0.3); rtol=1.0e-4)
    # the FaceEnd carries the slope it was given (default |tan theta|) and refuses one below it
    @test face.crossing_slope == abs(tan(0.3))
    sloped = FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        0.3,
        0.0,
        0.0,
        0.04,
        0.01,
        0.1;
        spacing_cap=Inf,
        face_axis=1,
        face_value=0.3,
        crossing_slope=arc_crossing_slope(1.0, 0.04, 0.3)
    )
    @test sloped.crossing_slope > face.crossing_slope &&
          sloped.apex_thickness == 2.0 * 0.01 * sloped.crossing_slope &&
          sloped.envelope_shear == 0.04 * sloped.crossing_slope &&
          sloped.layers >= face.layers
    @test face_end_record(sloped)["CrossingSlope"] == sloped.crossing_slope
    @test_throws ErrorException FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        0.3,
        0.0,
        0.0,
        0.04,
        0.01,
        0.1;
        spacing_cap=Inf,
        crossing_slope=0.5 * tan(0.3)
    )
end

@testset "synthetic rounded strip: fabricated and thin builds with one arc, two smooth joints and a concave corner" begin
    mktempdir() do directory
        census, mesh, inputs = build_strip_coupon(directory; fabricated=true, stem="fab")
        tubes = census["PrismTubes"]
        @test tubes["ArcTubes"]["Count"] == 2            # a top and a bottom arc tube
        @test tubes["ArcTubes"]["SharedSections"] == 4   # two joints x two placements
        @test tubes["ArcTubes"]["JointEnds"] == 4 && tubes["ArcTubes"]["PartSplits"] == 0
        @test census["Scope"]["ExhibitedClasses"] ==
              ["ArcSides", "ContinuationVertices", "ExteriorLoops"]
        loop = census["Scope"]["MetalLoops"][1]
        @test loop["Sides"] == 5 && loop["StraightSides"] == 4 && loop["ArcParts"] == 1
        @test tubes["TubeCount"] == 10
        arc_rows = [row for row in tubes["Tubes"] if haskey(row, "Arc")]
        @test length(arc_rows) == 2
        @test all(
            isapprox(row["Arc"]["SweepDegrees"], 90.0; atol=1.0e-9) for row in arc_rows
        )
        @test all(isapprox(row["Length"], 0.5 * pi; atol=1.0e-9) for row in arc_rows)
        joint_rows = [row for row in tubes["Tubes"] if haskey(row, "Joints")]
        @test length(joint_rows) == 4 && all(
            all(j["TiltRadians"] == 0.0 && !j["PlaneCut"] for j in row["Joints"]) for
            row in joint_rows
        )
        @test tubes["Prisms"] > 0 &&
              tubes["Pyramids"] > 0 &&
              haskey(tubes["Quality"], "Prism")
        @test isfile(mesh)
        thin_census, thin_mesh, _ =
            build_strip_coupon(directory; fabricated=false, stem="thin")
        @test thin_census["PrismTubes"]["ArcTubes"]["Count"] == 1
        @test thin_census["PrismTubes"]["ArcTubes"]["SharedSections"] == 2
        @test thin_census["PrismTubes"]["TubeCount"] == 5
        @test isfile(thin_mesh)
    end
end

@testset "round 3 class (11), fix 11 (decisions 497 / 500 / 510; part M 1.3 / 1.5 S1 / S3): one circle serialised as two or three tagged runs builds fab + thin with the earlier tube owning every arc-arc shared section" begin
    mktempdir() do directory
        # S1: the strip's 90-degree arc (18 chords) split at the non-cardinal angle -55 degrees (after
        # chord 7) into ids 1 and 2 meeting at a smooth joint: the former ScopeGuard[ArcArcJoint]
        # stop (B2) is lifted; the shared section is the owner's (install order), SharedSections
        # gains TubesPerSide per arc-arc joint, every tube volume closes (no KeyError after the mesh).
        census, mesh, _ = build_strip_coupon(
            mkpath(joinpath(directory, "s1-fab"));
            fabricated=true,
            stem="s1-fab",
            split_chords=[7]
        )
        tubes = census["PrismTubes"]
        @test tubes["ArcTubes"]["Count"] == 4                 # two runs x top / bottom
        @test tubes["ArcTubes"]["SharedSections"] == 6        # 2 arc / straight joints + 1 arc-arc, x 2
        @test tubes["ArcTubes"]["JointEnds"] == 4 && tubes["ArcTubes"]["PartSplits"] == 0
        @test tubes["TubeCount"] == 12
        loop = census["Scope"]["MetalLoops"][1]
        @test loop["Sides"] == 6 && loop["ArcParts"] == 2 && length(loop["Arcs"]) == 2
        @test [arc["ArcId"] for arc in loop["Arcs"]] == [1, 2] && [arc["Chords"] for arc in loop["Arcs"]] == [7, 11]
        @test isapprox(sum(arc["SweepDegrees"] for arc in loop["Arcs"]), 90.0; atol=1.0e-9)
        arc_rows = [row for row in tubes["Tubes"] if haskey(row, "Arc")]
        @test length(arc_rows) == 4 &&
              isapprox(sum(row["Length"] for row in arc_rows), 2 * 0.5 * pi; atol=1.0e-9)
        @test tubes["Prisms"] > 0 && tubes["Pyramids"] > 0 && isfile(mesh)
        thin, thin_mesh, _ = build_strip_coupon(
            mkpath(joinpath(directory, "s1-thin"));
            fabricated=false,
            stem="s1-thin",
            split_chords=[7]
        )
        @test thin["PrismTubes"]["ArcTubes"]["Count"] == 2 &&
              thin["PrismTubes"]["ArcTubes"]["SharedSections"] == 3 &&
              thin["PrismTubes"]["TubeCount"] == 6 &&
              isfile(thin_mesh)
        # The control: the same chords under one id (the round-2b strip) - the split build carries
        # the same straight tubes and the same total arc length; its element count is within 5 %.
        control, _, _ = build_strip_coupon(
            mkpath(joinpath(directory, "control"));
            fabricated=true,
            stem="control"
        )
        @test control["PrismTubes"]["ArcTubes"]["SharedSections"] == 4
        elements(c) = c["PrismTubes"]["FarFieldBudgetPolicy"]["Elements"]
        @test abs(elements(census) - elements(control)) <= 0.05 * elements(control)
        # S3: three members (splits after chords 6 and 12) with the ids PERMUTED in loop order
        # (3, 1, 2): the owner of each arc-arc section is the earlier tube in install order (= loop
        # order), whatever the ids; the build equals the (1, 2, 3) build byte for byte - the ids
        # name the tubes, they do not order them.
        ordered, ordered_mesh, _ = build_strip_coupon(
            mkpath(joinpath(directory, "s3"));
            fabricated=true,
            stem="s3",
            split_chords=[6, 12]
        )
        permuted, permuted_mesh, _ = build_strip_coupon(
            mkpath(joinpath(directory, "s3p"));
            fabricated=true,
            stem="s3p",
            split_chords=[6, 12],
            ids=[3, 1, 2]
        )
        for c in (ordered, permuted)
            @test c["PrismTubes"]["ArcTubes"]["Count"] == 6 &&
                  c["PrismTubes"]["ArcTubes"]["SharedSections"] == 8
        end
        @test [arc["ArcId"] for arc in ordered["Scope"]["MetalLoops"][1]["Arcs"]] == [1, 2, 3]
        @test [arc["ArcId"] for arc in permuted["Scope"]["MetalLoops"][1]["Arcs"]] == [3, 1, 2]
        @test elements(ordered) == elements(permuted)
        @test bytes2hex(open(sha256, ordered_mesh)) ==
              bytes2hex(open(sha256, permuted_mesh))
    end
end

@testset "chord-count independence (V2, design A7 MINOR-7): 5-degree and 2.5-degree chords give the same tubes" begin
    mktempdir() do directory
        censuses = Dict{Float64, Any}()
        for step in (5.0, 2.5)
            sub = joinpath(directory, "chords-$step")
            mkpath(sub)
            censuses[step], _, inputs = build_strip_coupon(
                sub;
                fabricated=true,
                stem="fab",
                chord_degrees=step,
                labels_only=true
            )
            @test inputs.chords == round(Int, 90.0 / step)
        end
        a, b = censuses[5.0], censuses[2.5]
        @test a["Scope"]["MetalLoops"][1]["Sides"] ==
              b["Scope"]["MetalLoops"][1]["Sides"] ==
              5
        @test a["Scope"]["MetalLoops"][1]["Arcs"][1]["Chords"] == 18 &&
              b["Scope"]["MetalLoops"][1]["Arcs"][1]["Chords"] == 36
        for key in ("Centre", "Radius", "SweepDegrees", "Parts")
            @test a["Scope"]["MetalLoops"][1]["Arcs"][1][key] ==
                  b["Scope"]["MetalLoops"][1]["Arcs"][1][key]
        end
        @test a["InterfaceAreas"] == b["InterfaceAreas"]
        @test a["SemanticCorners"] == b["SemanticCorners"]
        # Full builds: the tube rows (stations, lengths, layers), the element counts and the
        # mesh file itself are the same under both chordings (the CAD carries no chord vertex).
        digests = Dict{Float64, String}()
        rows = Dict{Float64, Any}()
        counts = Dict{Float64, Any}()
        for step in (5.0, 2.5)
            sub = joinpath(directory, "full-$step")
            mkpath(sub)
            census, mesh, _ =
                build_strip_coupon(sub; fabricated=true, stem="fab", chord_degrees=step)
            digests[step] = bytes2hex(open(sha256, mesh))
            rows[step] = [
                (
                    row["Start"],
                    row["End"],
                    row["Length"],
                    row["Layers"],
                    row["Spacing"],
                    haskey(row, "Arc") ? row["Arc"]["SweepDegrees"] : nothing
                ) for row in census["PrismTubes"]["Tubes"]
            ]
            counts[step] = (
                census["PrismTubes"]["Prisms"],
                census["PrismTubes"]["Pyramids"],
                census["PrismTubes"]["FarFieldBudgetPolicy"]["Elements"]
            )
        end
        @test rows[5.0] == rows[2.5]
        @test counts[5.0] == counts[2.5]
        @test digests[5.0] == digests[2.5]
    end
end

@testset "round 3 class (7) Option B (decision 510 O1; G.7.6 S-G7-a / S-G7-b = G.4.5 S-G4-a): a 3-chord and a 1-chord arc of a 104 R circle build fab + thin; the tag-seeded collar circle equals the refit's to 1e-10" begin
    mktempdir() do directory
        # S-G7-a: the lead's 3-chord convex arc (0.25 R chords, sweep 7.2e-3 rad on 104 R) - the
        # former >= 4-chord guard refused it - builds fab + thin: one arc tube per placement, two
        # smooth joints, every tube volume closed. The SAME arc chorded 6 times takes the untagged-
        # fit collar path (>= 4 chords: today's spelling, asserted against the propagated circle):
        # the thin meshes (no collar: one spelling) are byte-identical, the fabricated tubes equal
        # and the fabricated element counts within 1 % (the two collar circles agree to ~1e-13 R_arc
        # - G.7.3 - not to the last ulp, so the fabricated meshes are not byte-identical: Option B).
        digests = Dict{Tuple{Int, Bool}, String}()
        counts = Dict{Tuple{Int, Bool}, Any}()
        elements = Dict{Tuple{Int, Bool}, Int}()
        for fabricated in (true, false), chords in (3, 6)
            sub = mkpath(joinpath(directory, "lead-$chords-$(fabricated ? "fab" : "thin")"))
            census, mesh, inputs = build_arc_coupon(
                sub,
                write_short_arc_lead_inputs;
                fabricated=fabricated,
                stem="lead",
                chords=chords
            )
            @test inputs.chords == chords
            loop = census["Scope"]["MetalLoops"][1]
            @test loop["Arcs"][1]["Chords"] == chords &&
                  loop["Arcs"][1]["Parts"] == 1 &&
                  isapprox(loop["Arcs"][1]["Radius"], inputs.rho; atol=1.0e-12)
            tubes = census["PrismTubes"]
            @test tubes["ArcTubes"]["Count"] == (fabricated ? 2 : 1)
            @test tubes["ArcTubes"]["SharedSections"] == (fabricated ? 4 : 2)
            @test tubes["Prisms"] > 0 && tubes["Pyramids"] > 0 && isfile(mesh)
            digests[(chords, fabricated)] = bytes2hex(open(sha256, mesh))
            counts[(chords, fabricated)] =
                (tubes["Prisms"], tubes["Pyramids"], tubes["TubeCount"])
            elements[(chords, fabricated)] = tubes["FarFieldBudgetPolicy"]["Elements"]
            if chords == 3 && fabricated
                # The collar record: the propagated circle under the 3-chord arc vs the refit of the
                # 6-chord offset polygon of the same loop - one circle to 1e-10 R_arc (P-G7.2).
                loops = read_boundary(inputs.boundary)
                collar = -3 * 0.5
                record = offset_loop(loops[1], collar, 1.0e-7 * 0.5)
                @test record.short &&
                      length(record.runs) == 1 &&
                      record.runs[1].id == 1 &&
                      length(record.runs[1].edge_indices) == 3
                seeded = fit_carried_run(record.points, record.runs[1], 1.0e-7 * 0.5)
                @test isapprox(seeded.radius, inputs.rho + 3 * 0.5; atol=1.0e-12) &&
                      seeded.center == inputs.centre
                six = write_short_arc_lead_inputs(
                    mkpath(joinpath(directory, "six-offset"));
                    chords=6
                )
                six_record =
                    offset_loop(read_boundary(six.boundary)[1], collar, 1.0e-7 * 0.5)
                @test !six_record.short
                refit = only(circular_arc_runs(six_record.points, 1.0e-7 * 0.5))
                @test hypot(
                    refit.center[1] - seeded.center[1],
                    refit.center[2] - seeded.center[2]
                ) <= 1.0e-10 * inputs.rho
                @test abs(refit.radius - seeded.radius) <= 1.0e-10 * inputs.rho
            end
        end
        for fabricated in (true, false)
            @test counts[(3, fabricated)] == counts[(6, fabricated)]
            @test abs(elements[(3, fabricated)] - elements[(6, fabricated)]) <=
                  0.01 * elements[(6, fabricated)]
        end
        @test digests[(3, false)] == digests[(6, false)]
        # S-G7-b: a 1-chord CONCAVE arc (sweep 0.1 degrees on 104 R: 0.18 R of chord; the collar
        # shrinks its circle to 101 R) builds fab + thin.
        for fabricated in (true, false)
            sub = mkpath(joinpath(directory, "one-$(fabricated ? "fab" : "thin")"))
            census, mesh, inputs = build_arc_coupon(
                sub,
                write_short_arc_lead_inputs;
                fabricated=fabricated,
                stem="one",
                chords=1,
                sweep=deg2rad(0.1),
                convex=false
            )
            @test inputs.chords == 1 &&
                  census["Scope"]["MetalLoops"][1]["Arcs"][1]["Chords"] == 1 &&
                  census["Scope"]["MetalLoops"][1]["Arcs"][1]["Sign"] == -1
            @test census["PrismTubes"]["ArcTubes"]["Count"] == (fabricated ? 2 : 1) &&
                  isfile(mesh)
        end
    end
end

@testset "decision 391 MAJOR-2 (ii) lifted by round 2b (decision 437 (3)): arc face ends and arc joints build inside the TESTED ranges, fail closed beyond" begin
    # Both guards stay in the recipe scope list (mesh_stage_contract.py spells the same list and
    # the same bounds); the loop end's tested turn bound lies inside the smooth range.
    ids = [guard[1] for guard in RECIPE_SCOPE_GUARDS]
    @test "ArcFaceEnds" in ids && "ArcJointTilt" in ids
    # Round 3 B2 (decision 510): the interim guards ArcArcJoint (class 11) and CollarFaceEnd
    # (class 5) are in the list; the arc face-end range is ONE named constant, (lowest built,
    # largest built) = (0.1, 70) degrees (decisions 497 / 510 MAJOR-1 / O6), the thin face-end
    # bound of 8B the largest built thin tilt.
    @test "ArcArcJoint" in ids && "CollarFaceEnd" in ids
    @test ARC_JOINT_TURN_BOUND == 1.6e-6 < ARC_SMOOTH_JOINT_TURN_BOUND == 5.0e-5
    # Round 3 class (6) (decision 556): the corner range is the BUILT range PER KIND - thin 2e-4 rad
    # .. 30 degrees (round 2b), fabricated 2e-4 rad .. 15 degrees (round 3 B4, production sizes).
    @test ARC_CORNER_JOINT_TURN_RANGE ==
          (fabricated=(2.0e-4, deg2rad(15.0)), thin=(2.0e-4, deg2rad(30.0))) &&
          arc_corner_joint_turn_range(true) == ARC_CORNER_JOINT_TURN_RANGE.fabricated &&
          arc_corner_joint_turn_range(false) == ARC_CORNER_JOINT_TURN_RANGE.thin &&
          ARC_FACE_END_TILT_RANGE ==
          (fabricated=(deg2rad(0.1), deg2rad(75.5)), thin=(deg2rad(0.1), deg2rad(70.0))) &&
          arc_face_end_tilt_range(true) == ARC_FACE_END_TILT_RANGE.fabricated &&
          arc_face_end_tilt_range(false) == ARC_FACE_END_TILT_RANGE.thin &&
          THIN_FACE_END_TILT_BOUND == deg2rad(70.0) == ARC_FACE_END_TILT_RANGE.thin[2]
    @test occursin(
        "0.0002 .. 0.2617993877991494 rad on a FABRICATED",
        scope_guard_statement("ArcJointTilt")
    ) && occursin(
        "0.0002 .. 0.5235987755982988 rad on a THIN",
        scope_guard_statement("ArcJointTilt")
    )
    @test occursin(
        string(ARC_SMOOTH_JOINT_TURN_BOUND),
        scope_guard_statement("ArcJointTilt")
    )
    @test occursin(
        "0.1 <= theta <= 70.0 degrees THIN / <= 75.5 degrees FABRICATED",
        scope_guard_statement("ArcFaceEnds")
    ) && occursin("inner node circle", scope_guard_statement("ArcFaceEnds"))
    # Round 3 B3 (fix 11): the ArcArcJoint guard keeps its id for the residual untested class
    # (distinct circles / opposite sigma); two runs of ONE circle build.
    @test occursin("NOT one circle", scope_guard_statement("ArcArcJoint")) &&
          occursin("32b0083dad90", scope_guard_statement("ArcArcJoint")) &&
          occursin("fix 11", scope_guard_statement("ArcArcJoint"))
    @test occursin("rho - h_face < 3 Radius", scope_guard_statement("CollarFaceEnd")) &&
          occursin("CORNER kink", scope_guard_statement("CollarFaceEnd"))
    @test occursin("THIN 70.0 degrees", scope_guard_statement("SteepFaceCrossing")) &&
          occursin("FABRICATED 75.5 degrees", scope_guard_statement("SteepFaceCrossing"))
    clearance(angle) = 0.03 / tan(0.5 * angle) + 0.02
    segments_of(inputs; lower=inputs.lower, upper=inputs.upper, fabricated=false) =
        metal_edge_segments(
            read_boundary(inputs.boundary),
            inputs.corners,
            clearance,
            lower,
            upper,
            1.0e-7 * 0.5;
            edge_size=0.01,
            corner_radius=0.1,
            fabricated=fabricated
        )
    mktempdir() do directory
        # SMOOTH joints: the exactly tangent strip, the loop end's 1e-6 rad and the tested 5e-5 rad
        # pass; 6e-5 rad (smooth-tagged, above the tested bound) fails closed naming the turn.
        for turn in (0.0, 1.0e-6, 1.01 * ARC_JOINT_TURN_BOUND, 5.0e-5)
            inside = write_strip_inputs(
                mkpath(joinpath(directory, "s$turn"));
                kink_degrees=rad2deg(turn)
            )
            @test count(s.kind == :arc for s in segments_of(inside)) == 1
        end
        message = guard_message(
            () -> segments_of(
                write_strip_inputs(
                    mkpath(joinpath(directory, "s6e-5"));
                    kink_degrees=rad2deg(6.0e-5)
                )
            )
        )
        @test occursin("ScopeGuard[ArcJointTilt]", message) &&
              occursin("smooth-joint turn", message) &&
              occursin("arc 1 part 1", message) &&
              occursin("straight side", message)
        @test isapprox(
            parse(Float64, match(r"turn of ([0-9.e+-]+) rad", message)[1]),
            6.0e-5;
            rtol=1.0e-6
        )
        # CORNER joints: 2e-4 rad and 30 degrees (the built ends of the THIN range) pass, with the
        # kinked arc end a corner retracted along the arc (A3 (3)); 1.5e-4 rad (a corner below the
        # range) and 31 degrees fail closed. Round 3 class (6) (decision 556): the FABRICATED range
        # is its own built range, 2e-4 rad .. 15 degrees - 15 passes, 16 fails closed naming the kind.
        for turn in (2.0e-4, deg2rad(30.0))
            kinked = write_strip_inputs(
                mkpath(joinpath(directory, "c$turn"));
                kink_degrees=rad2deg(turn)
            )
            arc = only(s for s in segments_of(kinked) if s.kind == :arc)
            @test arc.corners == (false, true) && arc.s_end < arc.span
            @test isapprox(arc.corner_angles[2], pi - turn; atol=1.0e-9)
        end
        for (turn, builds) in (
            (2.0e-4, true),
            (deg2rad(15.0), true),
            (deg2rad(16.0), false),
            (deg2rad(30.0), false)
        )
            kinked = write_strip_inputs(
                mkpath(joinpath(directory, "cfab$turn"));
                kink_degrees=rad2deg(turn)
            )
            if builds
                arc =
                    only(s for s in segments_of(kinked; fabricated=true) if s.kind == :arc)
                @test arc.corners == (false, true) && arc.s_end < arc.span
            else
                message = guard_message(() -> segments_of(kinked; fabricated=true))
                @test occursin("ScopeGuard[ArcJointTilt]", message) &&
                      occursin("corner-joint turn", message) &&
                      occursin("built FABRICATED range", message)
                @test isapprox(
                    parse(Float64, match(r"turn of ([0-9.e+-]+) rad", message)[1]),
                    turn;
                    rtol=1.0e-6
                )
            end
        end
        # Smooth joints and face ends build on both kinds.
        @test count(
            s.kind == :arc for s in segments_of(
                write_strip_inputs(
                    mkpath(joinpath(directory, "sfab"));
                    kink_degrees=rad2deg(5.0e-5)
                );
                fabricated=true
            )
        ) == 1
        for (turn, label) in ((1.5e-4, "below"), (deg2rad(31.0), "above"))
            message = guard_message(
                () -> segments_of(
                    write_strip_inputs(
                        mkpath(joinpath(directory, "c-$label"));
                        kink_degrees=rad2deg(turn)
                    )
                )
            )
            @test occursin("ScopeGuard[ArcJointTilt]", message) &&
                  occursin("corner-joint turn", message) &&
                  occursin("built THIN range", message)
            @test isapprox(
                parse(Float64, match(r"turn of ([0-9.e+-]+) rad", message)[1]),
                turn;
                rtol=1.0e-6
            )
        end
        # The census of a labels-only build spells both guards among the GuardedClasses.
        census, _, _ = build_strip_coupon(
            mkpath(joinpath(directory, "labels"));
            fabricated=true,
            stem="labels",
            labels_only=true
        )
        @test census["Scope"]["GuardedClasses"] == ids
        @test any(
            guard["Id"] == "ArcJointTilt" && guard["DetectedFrom"] == "build" for
            guard in census["Scope"]["Guards"]
        )
        # FACE ENDS: the arc cut by the x1 face at 0.1 / 0.5 / 2.1 / 8 (round 3 B2) and 15 / 45 / 70
        # degrees (round 2b) is a face end of the arc side (theta exact, the face named); at 71 degrees
        # (above the largest built) and at 0.05 degrees (below the lowest built) it fails closed;
        # an arc whose end tangent runs along the face (the perpendicular / tangential end) and an
        # arc at a box vertex fail closed.
        for (theta, chord) in (
                (0.1, 0.025),
                (0.5, 0.125),
                (2.1, 0.5),
                (8.0, 2.0),
                (15.0, 2.5),
                (45.0, 5.0),
                (70.0, 5.0)
            ),
            fabricated in (true, false)

            inputs = write_arc_face_end_inputs(
                mkpath(joinpath(directory, "fe$theta-$fabricated"));
                theta_degrees=theta,
                chord_degrees=chord
            )
            arc = only(
                s for s in segments_of(inputs; fabricated=fabricated) if s.kind == :arc
            )
            @test arc.face_ends[1] === nothing && arc.face_ends[2] !== nothing
            @test arc.face_ends[2].face == "x1" &&
                  isapprox(rad2deg(arc.face_ends[2].theta), theta; atol=1.0e-9)
            @test arc.joints[1] !== nothing && arc.corners == (false, false)
        end
        # the tops per kind (decisions 563 B / 566; the B4 record run): THIN 70, FABRICATED 75.5 - the
        # 71-degree thin end and the 76-degree fabricated end fail closed, the 74.3-degree fabricated
        # end (32dc558f4810's tilt) is a face end of the arc side
        message = guard_message(
            () -> segments_of(
                write_arc_face_end_inputs(
                    mkpath(joinpath(directory, "fe71"));
                    theta_degrees=71.0
                )
            )
        )
        @test occursin("ScopeGuard[ArcFaceEnds]", message) &&
              occursin("above the largest built THIN 70.0", message)
        message = guard_message(
            () -> segments_of(
                write_arc_face_end_inputs(
                    mkpath(joinpath(directory, "fe76"));
                    theta_degrees=76.0,
                    rho=13.3
                );
                fabricated=true
            )
        )
        @test occursin("ScopeGuard[ArcFaceEnds]", message) &&
              occursin("above the largest built FABRICATED 75.5", message)
        steep_fab = only(
            s for s in segments_of(
                write_arc_face_end_inputs(
                    mkpath(joinpath(directory, "fe74p3-fab"));
                    theta_degrees=74.3,
                    rho=13.3
                );
                fabricated=true
            ) if s.kind == :arc
        )
        @test steep_fab.face_ends[2] !== nothing &&
              isapprox(rad2deg(steep_fab.face_ends[2].theta), 74.3; atol=1.0e-9)
        message = guard_message(
            () -> segments_of(
                write_arc_face_end_inputs(
                    mkpath(joinpath(directory, "fe0p05"));
                    theta_degrees=0.05,
                    chord_degrees=0.0125
                )
            )
        )
        @test occursin("ScopeGuard[ArcFaceEnds]", message) &&
              occursin("below the lowest built 0.1", message)
        @test isapprox(
            parse(Float64, match(r"tilt of ([0-9.e+-]+) degrees", message)[1]),
            0.05;
            rtol=1.0e-6
        )
        # The rho-parametrised fixture (round 3, part M 3.3): rho 3 at 45 degrees moves the face to
        # x1 = 2.12 (the notch sets the box); the arc's face end reads the same tilt, rho (1 - sin
        # theta) = 0.88 against the fixture tie's 0.44.
        rho_inputs = write_arc_face_end_inputs(
            mkpath(joinpath(directory, "fe45-rho3"));
            theta_degrees=45.0,
            rho=3.0
        )
        @test rho_inputs.rho == 3.0 &&
              isapprox(rho_inputs.upper[1], 3.0 * sind(45.0); atol=1.0e-12) &&
              length(rho_inputs.corners) == 5
        rho_arc = only(s for s in segments_of(rho_inputs) if s.kind == :arc)
        @test rho_arc.face_ends[2] !== nothing &&
              rho_arc.face_ends[2].face == "x1" &&
              isapprox(rad2deg(rho_arc.face_ends[2].theta), 45.0; atol=1.0e-9) &&
              rho_arc.arc.rho == 3.0
        chords = 9
        arc_points = [
            (cosd(-90.0 + 90.0 * k / chords), 1.0 + sind(-90.0 + 90.0 * k / chords)) for
            k = 0:chords
        ]
        x_face = arc_points[end][1]
        points = vcat([(-3.0, 0.0)], arc_points, [(x_face, 3.0), (-3.0, 3.0)])
        m = length(points)
        classes =
            [i == m - 2 || i == m - 1 || i == m ? "Continuation" : "Physical" for i = 1:m]
        arcs = Vector{Union{Nothing, NamedTuple}}(nothing, m)
        for i = 2:(chords + 1)
            arcs[i] = (id=1, centre=(0.0, 1.0), radius=1.0, sign=1)
        end
        joints = Vector{Union{Nothing, NamedTuple}}(nothing, m)
        joints[2] = (turn=0.0, smooth=true)
        loop = (
            conductor=1,
            plane=0.0,
            hole=false,
            points=points,
            classes=classes,
            arcs=arcs,
            joints=joints
        )
        message = guard_message(
            () -> metal_edge_segments(
                [loop],
                Tuple{Float64, Float64, Float64}[],
                clearance,
                [-3.0, -3.0],
                [x_face, 3.0],
                1.0e-9;
                edge_size=0.01,
                corner_radius=0.1
            )
        )
        @test occursin("ScopeGuard[ArcFaceEnds]", message) &&
              occursin("exactly perpendicular", message)
        corner_points = vcat([(-3.0, 0.0)], arc_points, [(-3.0, 3.0)])     # the arc end (1, 1) is the box vertex (x1, y1)
        mc = length(corner_points)
        corner_classes = [i == mc - 1 || i == mc ? "Continuation" : "Physical" for i = 1:mc]
        corner_arcs = Vector{Union{Nothing, NamedTuple}}(nothing, mc)
        for i = 2:(chords + 1)
            corner_arcs[i] = (id=1, centre=(0.0, 1.0), radius=1.0, sign=1)
        end
        corner_joints = Vector{Union{Nothing, NamedTuple}}(nothing, mc)
        corner_joints[2] = (turn=0.0, smooth=true)
        corner_loop = (
            conductor=1,
            plane=0.0,
            hole=false,
            points=corner_points,
            classes=corner_classes,
            arcs=corner_arcs,
            joints=corner_joints
        )
        message = guard_message(
            () -> metal_edge_segments(
                [corner_loop],
                Tuple{Float64, Float64, Float64}[],
                clearance,
                [-3.0, -3.0],
                [1.0, 1.0],
                1.0e-9;
                edge_size=0.01,
                corner_radius=0.1
            )
        )
        @test occursin("ScopeGuard[ArcFaceEnds]", message)
        # Round 3 class (11), fix 11 (decisions 497 / 500 / 510; part M 1.3): two DISTINCT tagged arc
        # runs of ONE circle (ids 1 and 2, 4 chords each) meeting at a smooth joint (turn 0,
        # JointSmooth 1) BUILD their sides (two arc sides, each with the joint vertex as a smooth
        # joint: no FreeEdgeEnds, no guard) - the B2 interim ScopeGuard[ArcArcJoint] is lifted for
        # one circle; the same chords under ONE id (the exact part split) stay the control. The
        # residual guard: the second run on a DISTINCT circle (radius 1.001 through the joint and
        # the face-end vertex) fails closed by name; an opposite ArcSign disagrees with the plan-view
        # metal side before the joint is reached.
        split_points =
            [(cosd(-90.0 + 45.0 * k / 8), 1.0 + sind(-90.0 + 45.0 * k / 8)) for k = 0:8]
        x_split = split_points[end][1]            # the arc leaves the x1 face at a 45-degree tilt
        split_loop_points = vcat([(-3.0, 0.0)], split_points, [(x_split, 3.0), (-3.0, 3.0)])
        ms = length(split_loop_points)
        split_classes = [
            i == ms - 2 || i == ms - 1 || i == ms ? "Continuation" : "Physical" for i = 1:ms
        ]
        function split_loop(
            ids;
            second=(centre=(0.0, 1.0), radius=1.0, sign=1),
            points=split_loop_points
        )
            split_arcs = Vector{Union{Nothing, NamedTuple}}(nothing, ms)
            for i = 2:9
                tag = ids[i - 1] == 2 ? second : (centre=(0.0, 1.0), radius=1.0, sign=1)
                split_arcs[i] = (id=ids[i - 1], tag...)
            end
            split_joints = Vector{Union{Nothing, NamedTuple}}(nothing, ms)
            split_joints[2] = (turn=0.0, smooth=true)
            split_joints[6] = (turn=0.0, smooth=true)       # the joint vertex of the two runs
            return (
                conductor=1,
                plane=0.0,
                hole=false,
                points=points,
                classes=split_classes,
                arcs=split_arcs,
                joints=split_joints
            )
        end
        split_segments(loop) = metal_edge_segments(
            [loop],
            Tuple{Float64, Float64, Float64}[],
            clearance,
            [-3.0, -3.0],
            [x_split, 3.0],
            1.0e-9;
            edge_size=0.01,
            corner_radius=0.1
        )
        two_runs = split_segments(split_loop(vcat(fill(1, 4), fill(2, 4))))
        arc_sides = [s for s in two_runs if s.kind == :arc]
        @test length(arc_sides) == 2 && [s.arc.id for s in arc_sides] == [1, 2]
        joint_vertex = split_loop_points[6]
        @test norm(arc_sides[1].stop .- [joint_vertex...]) <= 1.0e-12 &&
              norm(arc_sides[2].start .- [joint_vertex...]) <= 1.0e-12
        @test arc_sides[1].joints[2] !== nothing &&
              arc_sides[2].joints[1] !== nothing &&
              arc_sides[1].face_ends[2] === nothing &&
              arc_sides[2].face_ends[2] !== nothing
        control = split_segments(split_loop(fill(1, 8)))
        @test count(s.kind == :arc for s in control) == 1 &&
              only(s for s in control if s.kind == :arc).face_ends[2] !== nothing
        # Distinct circles: the second run's vertices on the circle of radius 1.001 through the
        # joint vertex and the face-end vertex (its centre on their bisector, towards (0, 1)).
        end_vertex = split_loop_points[10]
        mid = 0.5 .* (joint_vertex .+ end_vertex)
        half = 0.5 * hypot((end_vertex .- joint_vertex)...)
        towards = (0.0, 1.0) .- mid
        towards = towards ./ hypot(towards...)
        c2 = mid .+ sqrt(1.001^2 - half^2) .* towards
        a_j = atan(joint_vertex[2] - c2[2], joint_vertex[1] - c2[1])
        a_e = atan(end_vertex[2] - c2[2], end_vertex[1] - c2[1])
        distinct_points = copy(split_loop_points)
        for k = 1:3
            angle = a_j + (a_e - a_j) * k / 4
            distinct_points[6 + k] =
                (c2[1] + 1.001 * cos(angle), c2[2] + 1.001 * sin(angle))
        end
        message = guard_message(
            () -> split_segments(
                split_loop(
                    vcat(fill(1, 4), fill(2, 4));
                    second=(centre=c2, radius=1.001, sign=1),
                    points=distinct_points
                )
            )
        )
        @test occursin("ScopeGuard[ArcArcJoint]", message) &&
              occursin("DISTINCT circles", message) &&
              occursin("arc 1 part 1", message) &&
              occursin("meets arc 2 part 1", message)
        # Opposite metal sides: the same circle tagged sign -1 on the second run disagrees with the
        # plan-view metal side before the joint is reached (the ArcSign check); the sigma reading of
        # the guard is exercised on the sides directly.
        message = guard_message(
            () -> split_segments(
                split_loop(
                    vcat(fill(1, 4), fill(2, 4));
                    second=(centre=(0.0, 1.0), radius=1.0, sign=-1)
                )
            )
        )
        @test occursin("disagrees with the tagged ArcSign", message)
    end
end

@testset "round 2b (decision 437 (3)): the four synthetic arc FULL builds at the test sizes" begin
    mktempdir() do directory
        # (1) / (2) arc face ends at 15 / 45 / 70 degrees, fabricated and thin: the arc tube ends ON
        # the face (design A2 for arcs) with the face-end record on the ArcTube row, the cap entities
        # on the face plane matched by their conic centroids, every gate passed.
        # (decision 475 MINOR-1: the 15-degree face end closes the shallow end of the lifted range;
        # 2.5-degree chords so that the 15-degree arc keeps the >= 4-chord guard's four chords)
        # Round 3 B2, 9L-b (decisions 497 / 510 O6; DESIGN R9, part M 3.2): the LOW end is BUILT at
        # 0.1 / 0.5 / 2.1 / 8 degrees (4 chords each; 0.1 = the admitted floor, 2.1 = the 32b0083dad90
        # angle), fab + thin; regime I with one sheared layer, continuous towards theta -> 0+ (the
        # production-size record is the round-3 B2 cluster run).
        # Round 3 9H (part M 3.3 Fact 2; decisions 491 / 510 O9): the arc face end's formulas read
        # the arc's CROSSING SLOPE (FaceEnds[].CrossingSlope = max |s'(u)| over the section, >= tan
        # theta, -> tan theta as rho -> inf). The default fixture ties rho = 3 R / sin theta, which
        # at 70 degrees leaves the inner node circle only 2.4 envelopes from the face (part M 3.3
        # Fact 1: a fixture artefact) and, with the slope-aware block, fails the tetrahedral gate
        # (0.0094 < 0.01): the 70-degree case takes the rho-parametrised (notch) fixture. Round 3 B4
        # review (decision 579 MAJOR-3 (b)): the fixture must be DOMINATED by a built-and-passed case
        # (ARC_FACE_END_BUILT_CASES: tilt >= 70 and margin <= the fixture's) - rho 3.5 (5.3 envelopes)
        # is not (the floor is fe75p5r13p3's 10.66 fabricated / fe70r13p3's 10.28 thin), rho 8 (12.1
        # envelopes at the test-size envelope 0.04; slope 1.05 tan theta) is.
        for (theta, chord, rho) in (
                (0.1, 0.025, nothing),
                (0.5, 0.125, nothing),
                (2.1, 0.5, nothing),
                (8.0, 2.0, nothing),
                (15.0, 2.5, nothing),
                (45.0, 5.0, nothing),
                (70.0, 5.0, 8.0)
            ),
            fabricated in (true, false)

            census, mesh, inputs = build_arc_coupon(
                mkpath(joinpath(directory, "fe$theta-$fabricated")),
                write_arc_face_end_inputs;
                fabricated=fabricated,
                stem="fe",
                theta_degrees=theta,
                chord_degrees=chord,
                rho=rho
            )
            tubes = census["PrismTubes"]
            arc_rows = [row for row in tubes["Tubes"] if haskey(row, "Arc")]
            @test length(arc_rows) == (fabricated ? 2 : 1) &&
                  all(haskey(row, "FaceEnds") for row in arc_rows)
            for row in arc_rows
                record = only(row["FaceEnds"])
                @test record["Face"] == "x1" &&
                      isapprox(record["ThetaDegrees"], theta; atol=1.0e-9)
                @test record["Regime"] == "I"
                # 9H: the recorded crossing slope is the arc's (rho, envelope, face distance), above
                # tan theta by the concave side's excess (5.9 % on the default fixture at 45, 5.3 % at
                # 70 / rho 8, within 1.5 % below 15 degrees), and the block's apex / shear read it.
                slope = arc_crossing_slope(
                    row["Arc"]["Radius"],
                    tubes["Section"]["Radius"] + tubes["Section"]["PyramidHeight"],
                    abs(inputs.upper[1] - row["Arc"]["Centre"][1])
                )
                @test record["CrossingSlope"] == slope &&
                      tand(theta) * (1.0 - 1.0e-12) <= slope <= 1.13 * tand(theta)
                theta <= 15.0 && @test slope <= 1.015 * tand(theta)
                # (decision 579: the node-circle margin is recorded and dominated by a built case)
                margin = (row["Arc"]["Radius"] - abs(inputs.upper[1] - row["Arc"]["Centre"][1])) /
                         (tubes["Section"]["Radius"] + tubes["Section"]["PyramidHeight"])
                @test record["NodeCircleMargin"] == margin &&
                      arc_face_end_dominating_case(fabricated, deg2rad(theta), margin) !== nothing
                # (8A, decision 563: at or below 70 degrees the block's pyramid height is the regular one)
                @test record["PyramidHeight"] == tubes["Section"]["PyramidHeight"] &&
                      record["ApexThickness"] == 2.0 * record["PyramidHeight"] * slope &&
                      record["EnvelopeShear"] ==
                      (tubes["Section"]["Radius"] + tubes["Section"]["PyramidHeight"]) *
                      slope
                # The arc tube travels against the loop here (orientation -sigma Nz): the face
                # end is the tube's START station; its axis point lies on the face.
                face_point = record["End"] == "end" ? row["EndPoint"] : row["StartPoint"]
                @test isapprox(face_point[1], inputs.upper[1]; atol=1.0e-9)
            end
            @test tubes["FaceEnds"]["Count"] == length(arc_rows) &&
                  tubes["ArcTubes"]["SharedSections"] == length(arc_rows)
            @test tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
            @test tubes["Quality"]["Prism"]["PositiveOrientation"] &&
                  tubes["Quality"]["Pyramid"]["PositiveOrientation"]
            # Every node of the arc tube's last station lies on the face, none beyond it.
            # (the end polygon of a section: the edge point + Rings x (Sectors + 1) nodes on the face)
            # (the face cut stretches the section by 1 / cos theta along the face)
            near = nodes_near(
                mesh,
                [inputs.face_point[1], inputs.face_point[2], fabricated ? 0.1 : 0.0],
                2.0 * tubes["Section"]["Radius"] / cosd(theta)
            )
            @test count(p -> abs(p[1] - inputs.upper[1]) <= 1.0e-9, near) >=
                  1 + tubes["Section"]["Rings"] * (tubes["Section"]["Sectors"] + 1)
            @test all(p[1] <= inputs.upper[1] + 1.0e-9 for p in near)
        end
        # (3) the 5e-5-rad smooth joint (fabricated): one shared section owned by the arc, the
        # straight tubes' end sheared by the tilt (PlaneCut), no corner at the joint; (4) the
        # 2e-4-rad corner joint and the 30-degree kinked arc / line joint (THIN: the tested kind):
        # the kinked arc end is an INVARIANT corner (kappa_reg under the gate), the arc tube
        # retracted along its arc by the corner clearance, no shared section at that end. (The
        # thin twin of (3) and the 2e-4 thin corner at these TEST sizes meet the decision-353
        # 90-degree seed lottery at another corner on some sizes; the fabricated corner joint fails
        # in Gmsh's surface mesher at the production sizes and stays guarded: round2b-impl REPORT
        # section 3.)
        census, _, _ = build_strip_coupon(
            mkpath(joinpath(directory, "smooth"));
            fabricated=true,
            stem="smooth",
            kink_degrees=rad2deg(5.0e-5)
        )
        tubes = census["PrismTubes"]
        tilted = [
            j for row in tubes["Tubes"] for
            j in get(row, "Joints", []) if j["TiltRadians"] > 0.0
        ]
        @test length(tilted) == 2 && all(
            isapprox(j["TiltRadians"], 5.0e-5; rtol=1.0e-6) && j["PlaneCut"] for j in tilted
        )
        @test tubes["ArcTubes"]["SharedSections"] == 4 &&
              tubes["ArcTubes"]["JointEnds"] == 4
        @test length(census["SemanticCorners"]) == 3 &&
              tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        # (the 2e-4-rad thin corner at these test sizes meets the lottery at a legacy corner, 4.34 > 4.0:
        # its full build is the production-size record)
        # Round 3 class (6) (decisions 491 / 556): the FABRICATED kinked arc / line joint builds too,
        # once the E4 root cause is fixed (install_tube_curves!: a cap ray on the periodic trench wall
        # gets its nodes in increasing curve parameter) - here the top of the built fabricated range,
        # 15 degrees (the production-size builds at 2e-4 rad / 0.1 / 5.33 / 8.53 / 15 degrees are the
        # round-3 B4 record).
        for (fabricated, turn, stem) in
            ((false, deg2rad(30.0), "kink"), (true, deg2rad(15.0), "kink-fab"))
            census, _, _ = build_strip_coupon(
                mkpath(joinpath(directory, stem));
                fabricated=fabricated,
                stem=stem,
                kink_degrees=rad2deg(turn)
            )
            tubes = census["PrismTubes"]
            measures = census["SeedQualityOptimization"]["CornerMeasures"]
            kinked = only(r for r in measures if r["Point"] == [1.0, 1.0, 0.0])
            @test kinked["Kind"] == "Invariant" &&
                  kinked["Passed"] &&
                  kinked["After"] <= 5.0
            @test tubes["ArcTubes"]["SharedSections"] == (fabricated ? 2 : 1)
            for row in tubes["Tubes"]
                haskey(row, "Arc") || continue
                # The arc of radius 1 sweeps 90 degrees (span pi / 2): its kinked end is retracted
                # by the corner clearance along the arc (whichever tube end the travel puts it at).
                @test row["Length"] < 0.5 * pi - 1.0e-6 &&
                      (row["Start"] > 0.0) != (row["End"] < 0.5 * pi - 1.0e-6)
                @test isapprox(row["CornerAngles"][2], pi - turn; atol=1.0e-9)
            end
            @test tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        end
    end
end

@testset "round 3 B4 review (decision 579 MAJOR-3 (b)): an arc face end is admitted only where a built-and-passed case of its kind dominates it (tilt and node-circle margin)" begin
    # The table is PBS 59701's built cases with the margins measured on the identical fixture geometry
    # (every case dominates itself: the margins are rounded DOWN); the tilt range is its projection.
    for fabricated in (true, false)
        cases = arc_face_end_built_cases(fabricated)
        range = arc_face_end_tilt_range(fabricated)
        @test minimum(c.tilt for c in cases) == rad2deg(range[1]) &&
              isapprox(maximum(c.tilt for c in cases), rad2deg(range[2]); atol=1.0e-12)
        @test all(c.margin > 1.0 for c in cases)
        for c in cases
            @test arc_face_end_dominating_case(fabricated, deg2rad(c.tilt), c.margin) === c
        end
    end
    @test length(ARC_FACE_END_BUILT_CASES.fabricated) == 9 &&
          length(ARC_FACE_END_BUILT_CASES.thin) == 7
    # The two measured failures (the default fe70 fixture, 2.42 fabricated / 1.23 thin envelopes) are
    # dominated by nothing; the record's rho-13.3 cases are; a dominated synthetic (45 degrees, 16
    # envelopes: fe45) is admitted, an undominated one (59 degrees, 6.3 envelopes: the C1 0.5 R concave
    # fixture, part M 2.4) and a 74.32-degree end at 10 envelopes (below fe75p5r13p3's 10.66) refused.
    for failed in ARC_FACE_END_FAILED_CASES
        @test arc_face_end_dominating_case(failed.kind == "fabricated", deg2rad(failed.tilt), failed.margin) ===
              nothing
    end
    @test arc_face_end_dominating_case(true, deg2rad(45.0), 16.0).case == "fe45"
    @test arc_face_end_dominating_case(true, deg2rad(74.32), 54.0).case == "fe75p5r13p3"   # 32dc558f4810 fab
    @test arc_face_end_dominating_case(true, deg2rad(74.32), 10.0) === nothing
    @test arc_face_end_dominating_case(true, deg2rad(59.0), 6.3) === nothing
    @test arc_face_end_dominating_case(false, deg2rad(59.0), 6.3) === nothing
    @test arc_face_end_dominating_case(false, deg2rad(45.0), 8.0).case == "fe45"
    @test arc_face_end_dominating_case(false, deg2rad(46.0), 8.0) === nothing
    @test arc_face_end_dominating_case(false, deg2rad(70.0), 10.3).case == "fe70r13p3"
    @test arc_face_end_dominating_case(false, deg2rad(70.0), 10.2) === nothing
    # FaceEnd carries the margin (NaN on a straight tube's end: not recorded, bitwise record); an arc
    # margin at or below 1 is rejected (the Fact-1 condition).
    straight = FaceEnd(1, "x1", [1.0, 0.0, 0.0], deg2rad(45.0), 0.0, 0.0, 0.04, 0.01, 0.1; spacing_cap=Inf)
    @test isnan(straight.node_circle_margin) && !haskey(face_end_record(straight), "NodeCircleMargin")
    arc_end = FaceEnd(1, "x1", [1.0, 0.0, 0.0], deg2rad(45.0), 0.0, 0.0, 0.04, 0.01, 0.1; spacing_cap=Inf,
                      crossing_slope=1.06, node_circle_margin=15.5)
    @test face_end_record(arc_end)["NodeCircleMargin"] == 15.5
    @test_throws ErrorException FaceEnd(1, "x1", [1.0, 0.0, 0.0], deg2rad(45.0), 0.0, 0.0, 0.04, 0.01, 0.1;
                                        spacing_cap=Inf, crossing_slope=1.06, node_circle_margin=0.9)
    mktempdir() do directory
        # The default fe70 fixture (rho 1.596 um) is refused BY NAME on both kinds before any CAD - the
        # record run's two gate failures fail closed at the guard; the record's rho-13.3 fe70 builds
        # pass the guard on both kinds (labels-only census with the margin recorded).
        for fabricated in (true, false)
            message = guard_message(
                () -> build_arc_coupon(
                    mkpath(joinpath(directory, "fe70-default-$fabricated")),
                    write_arc_face_end_inputs;
                    fabricated=fabricated,
                    stem="fe70",
                    labels_only=true,
                    theta_degrees=70.0,
                    chord_degrees=5.0
                )
            )
            @test occursin("ScopeGuard[ArcFaceEnds]", message) &&
                  occursin("dominated by no built-and-passed $(fabricated ? "FABRICATED" : "THIN") case", message) &&
                  occursin("decision 579", message) &&
                  occursin("FAILED fe70 fabricated 70.0 / 2.4218", message)
            margin_text = match(r"margin rho \(1 - sin theta\) / envelope of ([0-9.e+-]+) envelopes", message)
            @test margin_text !== nothing && 2.3 < parse(Float64, margin_text[1]) < 2.5
            census, _, _ = build_arc_coupon(
                mkpath(joinpath(directory, "fe70-rho13p3-$fabricated")),
                write_arc_face_end_inputs;
                fabricated=fabricated,
                stem="fe70r13p3",
                labels_only=true,
                theta_degrees=70.0,
                chord_degrees=5.0,
                rho=13.3
            )
            records = [f for t in census["PrismTubeFaceEnds"]["Tubes"] for f in t["FaceEnds"]]
            @test length(records) == (fabricated ? 2 : 1) &&
                  all(19.0 < f["NodeCircleMargin"] < 21.0 && f["ThetaDegrees"] == 70.0 for f in records)
        end
    end
end

@testset "round 3 class (6) (decision 556): the E4 regression fixture - a fabricated arc corner whose cap ray on the periodic trench wall carries two or more interior nodes builds" begin
    # The E4 root cause (round3 B4 REPORT; Gmsh meshGFace.cpp buildConsecutiveListOfVertices): the
    # periodic surface mesher reads a bounding curve's nodes in storage order and assumes it increases
    # with the curve parameter; the ring-ordered cap-ray nodes ran against it on the fabricated arc
    # corner's trench wall, a zigzag that fails edge recovery as soon as the ray has TWO interior nodes
    # (three rings). The suite's EdgeSize 0.01 gives two rings (one interior node: no zigzag), so this
    # fixture takes EdgeSize 0.004 (three rings under the 0.05 transverse bound) with a 5-degree kink;
    # it fails on the unsorted order ("Impossible to mesh periodic surface") and builds with the rule.
    mktempdir() do directory
        census, _, _ = build_arc_coupon(
            directory,
            write_strip_inputs;
            fabricated=true,
            stem="e4",
            edge_size=0.004,
            kink_degrees=5.0
        )
        tubes = census["PrismTubes"]
        @test tubes["Section"]["Rings"] == 3
        @test tubes["ArcTubes"]["SharedSections"] == 2
        kinked = only(
            r for r in census["SeedQualityOptimization"]["CornerMeasures"] if
            r["Point"] == [1.0, 1.0, 0.0]
        )
        @test kinked["Kind"] == "Invariant" && kinked["Passed"]
        @test tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01 &&
              tubes["Quality"]["Prism"]["PositiveOrientation"] &&
              tubes["Quality"]["Pyramid"]["PositiveOrientation"]
        for row in tubes["Tubes"]
            haskey(row, "Arc") || continue
            # the kinked end is retracted along the arc by the corner clearance
            @test row["Length"] < 0.5 * pi - 1.0e-6
        end
    end
end

@testset "F6 / decision 391 MINOR-6: facing widths of ARC tube intervals (exact circle geometry)" begin
    # Synthetic side records of one plane (the fields metal_facing_widths reads): an arc side's
    # interval runs from the angle of s_start to that of s_end along its travel (theta(s) =
    # theta_start + sign(sweep) s / rho), its outward normal is sigma x the radial direction;
    # a straight side's interval is (start + s_start d, start + s_end d) with its normal.
    function arc_side(
        centre,
        rho,
        sigma,
        theta_a,
        theta_b;
        s_start=0.0,
        s_end=nothing,
        plane=0.0
    )
        sweep = theta_b - theta_a
        span = rho * abs(sweep)
        a = centre .+ rho .* [cos(theta_a), sin(theta_a)]
        b = centre .+ rho .* [cos(theta_b), sin(theta_b)]
        return (
            kind=:arc,
            start=a,
            stop=b,
            direction=sign(sweep) .* [-sin(theta_a), cos(theta_a)],
            normal=sigma .* [cos(theta_a), sin(theta_a)],
            span=span,
            plane=plane,
            untubed=false,
            s_start=s_start,
            s_end=s_end === nothing ? span : s_end,
            arc=(
                id=1,
                centre=centre,
                rho=rho,
                sigma=sigma,
                theta_start=theta_a,
                theta_end=theta_b,
                sweep=sweep,
                part=1,
                parts=1,
                run_sweep=sweep,
                chords=[]
            )
        )
    end
    function straight_side(a, b, normal; s_start=0.0, s_end=nothing, plane=0.0)
        d = b .- a
        span = norm(d)
        return (
            kind=:straight,
            start=a,
            stop=b,
            direction=d ./ span,
            normal=normal,
            span=span,
            plane=plane,
            untubed=false,
            s_start=s_start,
            s_end=s_end === nothing ? span : s_end,
            arc=nothing
        )
    end
    centre = [0.0, 0.0]
    # An annular metal strip of width w between two concentric quarter arcs (the outer arc
    # convex: dielectric outside, sigma +1; the inner one concave about the metal: sigma -1):
    # both sides read exactly w, at every chording (no sampling: the closest points lie at
    # a common angle).
    w = 0.121
    outer = arc_side(centre, 1.0 + w, 1.0, -pi / 2, 0.0)
    inner = arc_side(centre, 1.0, -1.0, 0.0, -pi / 2)
    widths, pairs = metal_facing_widths([outer, inner], 1.0e-9)
    @test widths ≈ [w, w] && length(pairs) == 1 && pairs[1].width ≈ w
    @test metal_facing_width([outer, inner], 1.0e-9) ≈ w
    # The intervals retracted at their ends (a corner clearance): still w where they overlap
    # in angle, Inf once the angular intervals no longer overlap.
    retracted = arc_side(centre, 1.0, -1.0, 0.0, -pi / 2; s_start=0.3, s_end=1.2)
    @test metal_facing_widths([outer, retracted], 1.0e-9)[1] ≈ [w, w]
    far = arc_side(centre, 1.0 + w, 1.0, pi / 4, pi / 2)
    @test metal_facing_widths([far, inner], 1.0e-9)[1] == [Inf, Inf]
    # The same two arcs with the metal OUTSIDE both (a dielectric annulus: the normals point at
    # each other) face across dielectric, not metal: no facing pair.
    dielectric_outer = arc_side(centre, 1.0 + w, -1.0, -pi / 2, 0.0)
    dielectric_inner = arc_side(centre, 1.0, 1.0, 0.0, -pi / 2)
    @test metal_facing_widths([dielectric_outer, dielectric_inner], 1.0e-9)[1] == [Inf, Inf]
    # An arc against a straight side: a concave arc of radius 1 about the origin (the metal
    # outside the circle, sigma -1) and a straight side along y = 1 + w (its metal below it,
    # normal +y) enclose a metal ribbon of width w: the closest pair is the foot of the centre
    # on the straight side against the arc's top point, exactly w; a straight interval starting
    # at x = 1 has its END as the closest point, against the arc point at the end's own angle
    # (48.3 degrees, on the arc): the exact distance from the end to the circle.
    top_arc = arc_side(centre, 1.0, -1.0, pi / 4, 3 * pi / 4)
    line = straight_side([-2.0, 1.0 + w], [2.0, 1.0 + w], [0.0, 1.0])
    widths, pairs = metal_facing_widths([top_arc, line], 1.0e-9)
    @test widths ≈ [w, w] &&
          pairs[1].points[2] ≈ [0.0, 1.0 + w] &&
          pairs[1].points[1] ≈ [0.0, 1.0]
    offset_line = straight_side([1.0, 1.0 + w], [2.0, 1.0 + w], [0.0, 1.0])
    widths, pairs = metal_facing_widths([top_arc, offset_line], 1.0e-9)
    expected = norm([1.0, 1.0 + w]) - 1.0
    @test widths ≈ [expected, expected] &&
          pairs[1].points[2] ≈ [1.0, 1.0 + w] &&
          pairs[1].points[1] ≈ [1.0, 1.0 + w] ./ norm([1.0, 1.0 + w])
    # A straight interval entirely past the arc's angular range (below and to the right, its
    # metal above it so that the metal lies between the two): the arc's END is the closest
    # point.
    beyond_line = straight_side([1.5, 0.5], [2.0, 0.5], [0.0, -1.0])
    widths, pairs = metal_facing_widths([top_arc, beyond_line], 1.0e-9)
    @test length(pairs) == 1 &&
          pairs[1].points[1] ≈ [cos(pi / 4), sin(pi / 4)] &&
          widths ≈ fill(norm([1.5, 0.5] .- [cos(pi / 4), sin(pi / 4)]), 2)
    # Two arcs on different circles: the closest points lie on the line of centres when both
    # arcs contain its direction (two concave arcs - the metal outside both circles, i.e.
    # between them - facing across the metal between them).
    left = arc_side([-1.0, 0.0], 0.5, -1.0, -pi / 4, pi / 4)       # bulging towards +x, metal to its right
    right = arc_side([1.0, 0.0], 0.5, -1.0, 3 * pi / 4, 5 * pi / 4) # bulging towards -x, metal to its left
    widths, pairs = metal_facing_widths([left, right], 1.0e-9)
    @test widths ≈ [1.0, 1.0] &&
          pairs[1].points[1] ≈ [-0.5, 0.0] &&
          pairs[1].points[2] ≈ [0.5, 0.0]
    # Adjacent sides (sharing an end) and untubed sides are not facing pairs.
    leg = straight_side([1.0 + w, 0.0], [1.0 + w, -1.0], [1.0, 0.0])
    @test metal_facing_widths([outer, leg], 1.0e-9)[1] == [Inf, Inf]
    untubed_inner = (inner..., untubed=true)
    @test metal_facing_widths([outer, untubed_inner], 1.0e-9)[1] == [Inf, Inf]
end

@testset "round 3 class (1) (decisions 491 / 510): corner-free coupons build with 0 semantic corners, the corner gates NotApplicable" begin
    # The contract gate (S-G1-c on the Julia side): [] is admitted iff Derivation.CornerFree is
    # recorded; a legacy contract keeps the >= 1 rule; a record beside corners fails closed.
    mktempdir() do directory
        contract = joinpath(directory, "legacy-empty.json")
        write(contract, "{\"Version\": 1, \"SemanticCorners\": []}")
        @test occursin(
            "at least one semantic corner",
            guard_message(() -> read_semantic_corners(contract, IDENTITY_RIGID_TRANSFORM))
        )
        # The record as derive_semantic_contract writes it: the exact rule text and integer counts.
        rule = escape_string(CORNER_FREE_RULE)
        record =
            "\"Derivation\": {\"CornerFree\": {\"Rule\": \"$rule\", \"ArcInteriorVertices\": 34, " *
            "\"SmoothJoints\": 4, \"BoxFaceCutEnds\": 4}}"
        mixed = joinpath(directory, "mixed.json")
        write(mixed, "{\"Version\": 1, \"SemanticCorners\": [[0.0, 0.0, 0.0]], $record}")
        @test occursin(
            "Derivation.CornerFree with 1 semantic corners",
            guard_message(() -> read_semantic_corners(mixed, IDENTITY_RIGID_TRANSFORM))
        )
        free = joinpath(directory, "free.json")
        write(free, "{\"Version\": 1, \"SemanticCorners\": [], $record}")
        @test read_semantic_corners(free, IDENTITY_RIGID_TRANSFORM) == NTuple{3, Float64}[]
        @test contract_is_corner_free(parse_json(read(free, String)))
        @test !contract_is_corner_free(parse_json(read(contract, String)))
        bad = joinpath(directory, "bad.json")
        write(
            bad,
            "{\"Version\": 1, \"SemanticCorners\": [], \"Derivation\": {\"CornerFree\": " *
            "{\"Rule\": \"$rule\", \"ArcInteriorVertices\": 0, \"SmoothJoints\": 0, \"BoxFaceCutEnds\": 4}}}"
        )
        @test occursin(
            "at least one arc vertex",
            guard_message(() -> read_semantic_corners(bad, IDENTITY_RIGID_TRANSFORM))
        )
        # Decision 524 MINOR-2: the Rule text is pinned (one spelling with the Python validator)
        # and a JSON boolean is not a count (Bool <: Integer in Julia), as on the Python side.
        for (name, broken) in (
            ("wrong-rule.json", replace(record, rule => "another rule")),
            (
                "boolean-count.json",
                replace(record, "\"SmoothJoints\": 4" => "\"SmoothJoints\": true")
            ),
            (
                "negative-count.json",
                replace(record, "\"BoxFaceCutEnds\": 4" => "\"BoxFaceCutEnds\": -1")
            ),
            (
                "float-count.json",
                replace(
                    record,
                    "\"ArcInteriorVertices\": 34" => "\"ArcInteriorVertices\": 34.0"
                )
            )
        )
            path = joinpath(directory, name)
            write(path, "{\"Version\": 1, \"SemanticCorners\": [], $broken}")
            @test occursin(
                "CornerFree must record the rule and the non-negative integer counts",
                guard_message(() -> read_semantic_corners(path, IDENTITY_RIGID_TRANSFORM))
            )
        end
    end
    # The empty-safe size laws (G.1.2 sites 5 and 8): the corner law is the identity without
    # corners (lc_tangent, the value it takes far from every corner) and no ball boundary exists.
    grading = CornerGrading(0.01, 2.0, 0.05, 0.1)
    @test corner_curve_size((0.3, 0.2, 0.0), NTuple{3, Float64}[], grading, 0.1, 0.675) ==
          0.1
    @test corner_curve_size((0.3, 0.2, 0.0), [(0.3, 0.2, 0.0)], grading, 0.1, 0.675) == 0.01
    # S-G1-a (the bent band, fab + thin) and S-G1-b (+ the ground's concave arc, fab + thin):
    # every build reaches its census with 0 corners, 0 CornerMeasures rows, 0 cap regions,
    # every joint smooth (one shared section per joint and placement), every straight end a
    # 2.5-degree box-face cut end, the concave arcs collapsed in the fabricated collar.
    mktempdir() do directory
        for (variant, ground) in (("a", false), ("b", true)), fabricated in (true, false)
            kind = fabricated ? "fab" : "thin"
            sub = joinpath(directory, "$variant-$kind")
            mkpath(sub)
            census, mesh, inputs = build_bent_band_coupon(
                sub;
                ground=ground,
                fabricated=fabricated,
                stem="$variant-$kind"
            )
            arcs = ground ? 3 : 2
            placements = fabricated ? 2 : 1
            @test census["SemanticCorners"] == [] &&
                  census["SemanticCornerKinds"] == [] &&
                  census["Corners"] == [] &&
                  census["InvariantCorners"]["Points"] == []
            @test census["CornerGates"] == CORNER_GATES_NOT_APPLICABLE
            optimization = census["SeedQualityOptimization"]
            @test optimization["CornerMeasures"] == [] &&
                  optimization["CornerAspectsAfter"] == [] &&
                  optimization["InvariantCorners"] == 0 &&
                  optimization["RequiredTetrahedra"] == 0 &&
                  optimization["RequiredMinimumScaledJacobianAfter"] === nothing &&
                  optimization["RequiredCellsBelowGateAfter"] == 0
            tubes = census["PrismTubes"]
            @test tubes["ArcTubes"]["Count"] == placements * arcs
            @test tubes["ArcTubes"]["SharedSections"] == 2 * placements * arcs
            @test tubes["TubeCount"] == 3 * placements * arcs
            @test tubes["CapRegions"]["Caps"] == 0 &&
                  tubes["CapRegions"]["Regions"] == [] &&
                  tubes["CapRegions"]["MinimumScaledJacobian"] === nothing
            face_ends = [f for row in tubes["Tubes"] for f in get(row, "FaceEnds", [])]
            @test length(face_ends) == placements * 2 * arcs
            @test all(
                isapprox(f["ThetaDegrees"], inputs.tilt_degrees; atol=1.0e-9) for
                f in face_ends
            )
            @test all(
                j["TiltRadians"] == 0.0 && !j["PlaneCut"] for
                row in tubes["Tubes"] if haskey(row, "Joints") for j in row["Joints"]
            )
            exterior = tubes["SizeLaws"]["Achieved"]["CornerExterior"]
            @test exterior["Corners"] == 0 &&
                  all(shell["Cells"] == 0 for shell in exterior["Shells"])
            @test tubes["Section"]["ReducedSides"] == 0
            @test census["Scope"]["ExhibitedClasses"] == vcat(
                ["ArcSides", "ContinuationVertices", "ExteriorLoops"],
                fabricated ? [] : ["ThinMetal"]
            )
            if fabricated
                # The concave arcs (0.5 R and 2 R, below the 3 R collar) collapse onto the miter of
                # their neighbours' offsets: a simple MiterOffset polygon per loop, no bridge.
                @test [
                    polygon["Construction"] for polygon in census["FootprintPolygons"]
                ] == fill("MiterOffset", ground ? 2 : 1)
            else
                @test census["ThinSheetSeams"]["Count"] == 0
            end
            @test isfile(mesh)
        end
    end
end
