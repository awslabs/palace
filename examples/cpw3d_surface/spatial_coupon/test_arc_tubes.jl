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
function write_strip_inputs(directory; chord_degrees=5.0, plane=0.0, kink_degrees=0.0)
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
    arcs = Vector{Any}(nothing, m)
    for i = 2:(n + 1)
        arcs[i] = (1, centre[1], centre[2], rho, 1)
    end
    joints = Vector{Any}(nothing, m)
    joints[2] = (0.0, 1)
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
    labels_only=false
)
    inputs = write_strip_inputs(
        directory;
        chord_degrees=chord_degrees,
        kink_degrees=kink_degrees
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
function write_arc_face_end_inputs(
    directory;
    theta_degrees=45.0,
    chord_degrees=5.0,
    plane=0.0
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
    # The left edge's claim ends 0.3 above the corner so that the box face y0 (-1.2) does not
    # coincide with the 3 R collar of the bottom edge and of the arc (y = -1.5: a tangent collar /
    # face contact leaves sliver tetrahedra between them).
    straight!((-1.0, 2.0), (-1.0, 0.3), 1.0)
    x1 = 0.0 + 3.0 * radius
    rho = x1 / sind(theta_degrees)
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
    # of x1: max(x) + R |gap_x| + R <= x1.
    for k = 1:n
        p, q = chord_points[k], chord_points[k + 1]
        t = (q[1] - p[1], q[2] - p[2]) ./ hypot(q[1] - p[1], q[2] - p[2])
        max(p[1], q[1]) + radius * abs(t[2]) + radius <= x1 + 1.0e-12 || break
        straight!(p, q, 1.0)
    end
    lower, upper = row_coupon_bounds(rows, radius, 0.1, 0.05)
    @assert upper[1] == x1 "the box x1 $(upper[1]) is not the designed $x1"
    polygon = Tuple{Float64, Float64}[(-1.0, 0.0)]
    append!(polygon, chord_points)                  # (0, 0) ... the face point
    push!(polygon, (x1, upper[2]))
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
    # the arc end on the face is a cut end (decision 320 for arcs: no corner).
    corners = [(-1.0, 0.0, plane), (-1.0, upper[2], plane)]
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
        # A corrupted tag (another centre) fails closed.
        broken = joinpath(directory, "broken.csv")
        open(broken, "w") do io
            for line in readlines(inputs.boundary)
                println(io, replace(line, ",1,0.0,1.0,1.0,1," => ",1,0.0,1.001,1.0,1,"))
            end
        end
        message =
            guard_message(() -> tagged_arc_runs(read_boundary(broken)[1], 1.0e-7 * 0.5))
        @test occursin("disagrees with the tagged circle", message) ||
              occursin("do not lie on the tagged circle", message)
    end
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

@testset "decision 391 MAJOR-2 (ii) lifted by round 2b (decision 437 (3)): arc face ends and arc joints build inside the TESTED ranges, fail closed beyond" begin
    # Both guards stay in the recipe scope list (mesh_stage_contract.py spells the same list and
    # the same bounds); the loop end's tested turn bound lies inside the smooth range.
    ids = [guard[1] for guard in RECIPE_SCOPE_GUARDS]
    @test "ArcFaceEnds" in ids && "ArcJointTilt" in ids
    @test ARC_JOINT_TURN_BOUND == 1.6e-6 < ARC_SMOOTH_JOINT_TURN_BOUND == 5.0e-5
    @test ARC_CORNER_JOINT_TURN_RANGE == (2.0e-4, deg2rad(30.0)) &&
          ARC_FACE_END_TILT_BOUND == deg2rad(70.0)
    @test occursin(
        string(ARC_SMOOTH_JOINT_TURN_BOUND),
        scope_guard_statement("ArcJointTilt")
    )
    @test occursin("70.0", scope_guard_statement("ArcFaceEnds"))
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
        # CORNER joints: 2e-4 rad and 30 degrees (the tested ends of the range) pass, with the
        # kinked arc end a corner retracted along the arc (A3 (3)); 1.5e-4 rad (a corner below the
        # tested range) and 31 degrees fail closed.
        for turn in (2.0e-4, deg2rad(30.0))
            kinked = write_strip_inputs(
                mkpath(joinpath(directory, "c$turn"));
                kink_degrees=rad2deg(turn)
            )
            arc = only(s for s in segments_of(kinked) if s.kind == :arc)
            @test arc.corners == (false, true) && arc.s_end < arc.span
            @test isapprox(arc.corner_angles[2], pi - turn; atol=1.0e-9)
            # The corner joint is tested on THIN coupons only: fabricated fails closed.
            message = guard_message(() -> segments_of(kinked; fabricated=true))
            @test occursin("ScopeGuard[ArcJointTilt]", message) &&
                  occursin("FABRICATED", message)
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
                  occursin("corner-joint turn", message)
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
        # FACE ENDS: the arc cut by the x1 face at 45 and 70 degrees is a face end of the arc side
        # (theta exact, the face named); at 71 degrees it fails closed; an arc whose end tangent
        # runs along the face (the perpendicular / tangential end) and an arc at a box vertex fail closed.
        for (theta, chord) in ((15.0, 2.5), (45.0, 5.0), (70.0, 5.0)),
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
        message = guard_message(
            () -> segments_of(
                write_arc_face_end_inputs(
                    mkpath(joinpath(directory, "fe71"));
                    theta_degrees=71.0
                )
            )
        )
        @test occursin("ScopeGuard[ArcFaceEnds]", message) &&
              occursin("above the tested 70.0", message)
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
    end
end

@testset "round 2b (decision 437 (3)): the four synthetic arc FULL builds at the test sizes" begin
    mktempdir() do directory
        # (1) / (2) arc face ends at 15 / 45 / 70 degrees, fabricated and thin: the arc tube ends ON
        # the face (design A2 for arcs) with the face-end record on the ArcTube row, the cap entities
        # on the face plane matched by their conic centroids, every gate passed.
        # (decision 475 MINOR-1: the 15-degree face end closes the shallow end of the lifted range;
        # 2.5-degree chords so that the 15-degree arc keeps the >= 4-chord guard's four chords)
        for (theta, chord) in ((15.0, 2.5), (45.0, 5.0), (70.0, 5.0)),
            fabricated in (true, false)

            census, mesh, inputs = build_arc_coupon(
                mkpath(joinpath(directory, "fe$theta-$fabricated")),
                write_arc_face_end_inputs;
                fabricated=fabricated,
                stem="fe",
                theta_degrees=theta,
                chord_degrees=chord
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
        for (fabricated, turn, stem) in ((false, deg2rad(30.0), "kink"),)
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
