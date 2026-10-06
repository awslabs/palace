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

@testset "decision 391 MAJOR-2 (ii): arc face ends and arc joints beyond the tested turn fail closed" begin
    # Both guards are in the recipe scope list (mesh_stage_contract.py spells the same list)
    # and the turn bound is the loop end's tested range.
    ids = [guard[1] for guard in RECIPE_SCOPE_GUARDS]
    @test "ArcFaceEnds" in ids && "ArcJointTilt" in ids
    @test ARC_JOINT_TURN_BOUND == 1.6e-6
    @test occursin(string(ARC_JOINT_TURN_BOUND), scope_guard_statement("ArcJointTilt"))
    clearance(angle) = 0.03 / tan(0.5 * angle) + 0.02
    segments_of(inputs; lower=inputs.lower, upper=inputs.upper) = metal_edge_segments(
        read_boundary(inputs.boundary),
        inputs.corners,
        clearance,
        lower,
        upper,
        1.0e-7 * 0.5;
        edge_size=0.01,
        corner_radius=0.1
    )
    mktempdir() do directory
        # POSITIVE: the exactly tangent strip (turn 0 at both joints) and a joint turning by
        # 1e-6 rad (inside the tested range, tagged smooth) pass; the census of a labels-only
        # build spells both guards among the GuardedClasses.
        tangent = write_strip_inputs(mkpath(joinpath(directory, "t0")))
        segments = segments_of(tangent)
        @test count(s.kind == :arc for s in segments) == 1
        @test all(s.face_ends == (nothing, nothing) for s in segments if s.kind == :arc)
        inside = write_strip_inputs(
            mkpath(joinpath(directory, "t1e-6"));
            kink_degrees=rad2deg(1.0e-6)
        )
        @test count(s.kind == :arc for s in segments_of(inside)) == 1
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
        # NEGATIVE (ArcJointTilt): a smooth-tagged joint turning by 1e-5 rad (A3 (1) smooth,
        # above the tested 1.6e-6), a 2e-4-rad corner joint and a 30-degree kinked arc / line
        # joint all fail closed at the same guard, naming the arc and the turn.
        for (turn, label) in
            ((1.0e-5, "smooth"), (2.0e-4, "corner"), (deg2rad(30.0), "kink"))
            kinked = write_strip_inputs(
                mkpath(joinpath(directory, label));
                kink_degrees=rad2deg(turn)
            )
            message = guard_message(() -> segments_of(kinked))
            @test occursin("ScopeGuard[ArcJointTilt]", message) &&
                  occursin("arc 1 part 1", message) &&
                  occursin("straight side", message)
            @test occursin(r"turn of ([0-9.e+-]+) rad", message)
            read_turn = parse(Float64, match(r"turn of ([0-9.e+-]+) rad", message)[1])
            @test isapprox(read_turn, turn; rtol=1.0e-6)
        end
        # The bound is exclusive-above: a turn 1 % above it fails, 1 % below it passes.
        for (factor, fails) in ((1.01, true), (0.99, false))
            edge = write_strip_inputs(
                mkpath(joinpath(directory, "bound-$factor"));
                kink_degrees=rad2deg(factor * ARC_JOINT_TURN_BOUND)
            )
            message = guard_message(() -> segments_of(edge))
            @test occursin("ScopeGuard[ArcJointTilt]", message) == fails
        end
        # NEGATIVE (ArcFaceEnds): a metal strip entering through the x0 face whose outer edge
        # bends about (0, 1) along a convex arc of radius 1 from (0, 0) and leaves the box along
        # the x1 face: with a 45-degree sweep the arc ends ON the x1 face at a 45-degree tilt
        # (a box-face cut end), with a 90-degree sweep its end tangent (0, 1) runs along the
        # face (the box face through the arc end at any angle); both fail closed at the arc
        # (a straight side ending on the box keeps its decision-320 treatment).
        for sweep in (45.0, 90.0)
            chords = 9
            arc_points = [
                (cosd(-90.0 + sweep * k / chords), 1.0 + sind(-90.0 + sweep * k / chords)) for k = 0:chords
            ]
            x_face = arc_points[end][1]
            points = vcat([(-3.0, 0.0)], arc_points, [(x_face, 3.0), (-3.0, 3.0)])
            m = length(points)
            classes = [
                i == m - 2 || i == m - 1 || i == m ? "Continuation" : "Physical" for i = 1:m
            ]
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
                  occursin("arc 1 part 1", message) &&
                  occursin("on the outer box", message)
        end
        # A straight side ending on the box keeps its face end (decision 320): the guard is
        # the arc's alone.
        @test any(
            s.face_ends != (nothing, nothing) || s.legacy_box_corners != (false, false) for
            s in segments if s.kind == :straight
        )
    end
end
