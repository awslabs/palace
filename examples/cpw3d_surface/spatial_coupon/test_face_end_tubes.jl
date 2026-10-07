# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# Oblique / non-axis straight metal edges in the prism-tube recipe (block (b) design
# AMENDMENT 1 A2 / A3 (3) / A6 / A9 family 2; supervisor decisions 302, 304, 320): a
# tube ending on a box face it does not cross perpendicularly ends ON the face (the
# over-long CAD solid intersected with the coupon box, the m-layer sheared end block),
# the box-face vertex class is decided by the vertex (angle-gated: theta == 0 keeps
# the legacy convention bitwise), and the corner clearance carries the pyramid
# envelope margin.

guard_message(f) =
    try
        f()
        ""
    catch e
        e isa ErrorException ? e.msg : sprint(showerror, e)
    end

@testset "FaceEnd parameters: lc_end and m at the production tube sizes (design A2 (4))" begin
    # Fabricated production tube: EdgeSize 0.25 nm, ratio 2, 7 rings -> R 31.75 nm,
    # h_K 16 nm, h_pyr 8 nm, r_env 39.75 nm; TangentialSize 50 nm.
    r_env, h_pyr, lc = 0.03975, 0.008, 0.05
    f45 = FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        deg2rad(45.0),
        1.0,
        0.0,
        r_env,
        h_pyr,
        lc;
        spacing_cap=Inf
    )
    @test f45.spacing == lc && f45.layers == 2 && f45.regime == 1
    @test f45.over_length ≈ r_env + lc && f45.envelope_shear ≈ r_env
    f12 = FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        deg2rad(12.6),
        tan(deg2rad(12.6)),
        0.0,
        r_env,
        h_pyr,
        lc;
        spacing_cap=Inf
    )
    @test f12.spacing == lc && f12.layers == 1
    f74 = FaceEnd(
        0,
        "y0",
        [0.0, -1.0, 0.0],
        deg2rad(74.3),
        tan(deg2rad(74.3)),
        0.0,
        r_env,
        h_pyr,
        lc;
        spacing_cap=Inf
    )
    # Without a cap (an unbounded condition ceiling) the A2 (4) formulas hold at any tilt.
    @test f74.spacing ≈ 4.0 * h_pyr * tan(deg2rad(74.3)) &&
          f74.layers == 3 &&
          f74.regime == 1
    # Thin production tube: EdgeSize 2 nm, 5 rings -> R 62 nm, h_K 32, h_pyr 16, r_env 78 nm.
    t45 = FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        deg2rad(45.0),
        1.0,
        0.0,
        0.078,
        0.016,
        lc;
        spacing_cap=Inf
    )
    @test t45.spacing ≈ 0.064 && t45.layers == 3
    # The shear per layer never exceeds lc_end / 2: thickness in [lc_end / 2, 3 lc_end / 2].
    for face_end in (f45, f12, f74, t45)
        @test face_end.envelope_shear / face_end.layers <=
              0.5 * face_end.spacing * (1.0 + 1.0e-9)
        range = face_end_record(face_end)["LayerThicknessRange"]
        @test range[1] >= 0.5 * face_end.spacing * (1.0 - 1.0e-9) &&
              range[2] <= 1.5 * face_end.spacing * (1.0 + 1.0e-9)
    end
    @test_throws ErrorException FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        0.0,
        0.0,
        0.0,
        r_env,
        h_pyr,
        lc;
        spacing_cap=Inf
    )
end

@testset "F2b: the end-spacing cap from the section's prism frames and the two regimes (design round 2, 3.2)" begin
    # The production sections: fabricated EdgeSize 0.25 nm, ratio 2, 7 rings; thin 2 nm, 5
    # rings; 30-degree sectors (the top section's rays -90 .. 180, the sheet's 180 .. 540).
    fab = TubeSection(0.00025, 2.0, 7, [-90.0 + 30.0 * j for j = 0:9], fill(2, 9))
    thin = TubeSection(
        0.002,
        2.0,
        5,
        [180.0 + 30.0 * j for j = 0:12],
        vcat(fill(1, 6), fill(2, 6))
    )
    # The analytic regular-layer prism condition reproduces the stored S2p 3-edge build
    # censuses (stage2-20261004/coupons-a/registration/sct002-S2p/build/spatial-3-edge-
    # 2c54db92028f[-thin]: Prism.MaximumJacobianCondition at SpacingMaximum) to 1e-9: the
    # measured 589 / 73.5 of the design text are these frames.
    @test section_prism_condition(fab, 0.04998265841897798) ≈ 589.3294590097332 rtol =
        1.0e-9
    @test section_prism_condition(thin, 0.049882932143200165) ≈ 73.51920231844923 rtol =
        1.0e-9
    # The frame with the smallest planar singular value is the second ring's (a, c, d)
    # triangle at its ring-2 vertex: sigma_2 = 0.33925 r_1 at ratio 2 / 30 degrees (the inner
    # triangle's edge-point frame reads sqrt(1 - cos 30) r_1 = 0.36603 r_1).
    frames = section_frame_singular_values(fab)
    @test minimum(last.(frames)) ≈ 0.33925104306172854 * 0.00025 rtol = 1.0e-12
    @test minimum(last.(frames)) < sqrt(1.0 - cosd(30.0)) * 0.00025
    @test minimum(last.(section_frame_singular_values(thin))) ≈ 0.33925104306172854 * 0.002 rtol =
        1.0e-12
    # The cap at the production ceiling 1000: 0.95 x 1000 x sigma_2 - 80.572 nm fabricated,
    # 644.58 nm thin (the design's ~80.5 / ~624 nm were read from measured conditions); its
    # exact value is fixed here. The condition at the cap reads exactly 950, and it is
    # monotone in the spacing above it.
    cap_fab = face_end_spacing_cap(fab, 1000.0)
    cap_thin = face_end_spacing_cap(thin, 1000.0)
    @test cap_fab ≈ 0.08057212272716054 rtol = 1.0e-12
    @test cap_thin ≈ 0.6445769818172843 rtol = 1.0e-12
    @test section_prism_condition(fab, cap_fab) ≈ 950.0 rtol = 1.0e-12
    @test section_prism_condition(thin, cap_thin) ≈ 950.0 rtol = 1.0e-12
    @test section_prism_condition(fab, 1.1 * cap_fab) >
          section_prism_condition(fab, cap_fab) >
          section_prism_condition(fab, 0.9 * cap_fab)
    @test face_end_spacing_cap(fab, 500.0) ≈ 0.5 * cap_fab rtol = 1.0e-12
    @test_throws ErrorException face_end_spacing_cap(fab, 0.0)
    @test_throws ErrorException face_end_spacing_cap(fab, Inf)
    # A ceiling so low that the cap falls under a frame's sigma_1 (not spacing-dominated).
    @test_throws ErrorException face_end_spacing_cap(fab, 1.5)
    # Regimes at the fabricated production sizes (h_pyr 8 nm, r_env 39.75 nm, lc 50 nm; the
    # design's worked values): the boundary 4 h_pyr |tan theta| = lc_cap at 68.3 degrees.
    r_env, h_pyr, lc = 0.03975, 0.008, 0.05
    face(theta; cap=cap_fab, env=r_env, pyr=h_pyr) = FaceEnd(
        1,
        "x1",
        [1.0, 0.0, 0.0],
        deg2rad(theta),
        tand(theta),
        0.0,
        env,
        pyr,
        lc;
        spacing_cap=cap
    )
    @test rad2deg(atan(cap_fab / (4.0 * h_pyr))) ≈ 68.3 atol = 0.05
    for theta in (12.6, 22.5, 45.0, 68.0, 68.3)
        f = face(theta)
        legacy = face(theta; cap=Inf)
        @test f.regime == 1 && f.spacing == legacy.spacing && f.layers == legacy.layers
        @test f.spacing <= cap_fab && f.spacing_cap == cap_fab
        @test f.apex_thickness ≈ 2.0 * h_pyr * tand(theta)
        @test face_end_record(f)["Regime"] == "I"
    end
    @test face(68.3).spacing ≈ 4.0 * h_pyr * tand(68.3) && face(68.3).layers == 3
    # Regime II: 70 deg lc_end = lc_cap, m 3, t_min 44.17 >= 43.96 nm; 74.3 deg (the census
    # crossing 32dc558f4810) m 6 (the uncapped rule: 113.8 nm / m 3); 75.5 deg m 9.
    for (theta, m) in ((68.4, 3), (70.0, 3), (74.3, 6), (75.5, 9))
        f = face(theta)
        slope = tand(theta)
        @test f.regime == 2 && f.spacing == cap_fab && f.layers == m
        @test f.layers >= max(1, ceil(Int, 2.0 * r_env * slope / cap_fab * (1.0 - 1.0e-9)))
        t_min = f.spacing - f.envelope_shear / f.layers
        @test t_min >= f.apex_thickness && f.apex_thickness ≈ 2.0 * h_pyr * slope
        @test f.envelope_shear / f.layers <= 0.5 * f.spacing
        record = face_end_record(f)
        @test record["Regime"] == "II" &&
              record["EndSpacingCap"] == cap_fab &&
              record["ApexThickness"] ≈ 2.0 * h_pyr * slope &&
              record["LayerThicknessRange"][1] ≈ t_min
    end
    @test face(70.0).spacing - face(70.0).envelope_shear / 3 ≈ 0.04417 atol = 2.0e-5
    @test isapprox(face(74.3; cap=Inf).spacing, 0.1138; atol=1.0e-4) &&
          face(74.3; cap=Inf).layers == 3
    # The validity ceiling 2 h_pyr |tan theta| < lc_cap: 78.77 degrees fabricated; 80 degrees
    # has no capped block (the build fails closed at ScopeGuard[SteepFaceCrossing]).
    @test rad2deg(atan(cap_fab / (2.0 * h_pyr))) ≈ 78.77 atol = 0.05
    @test face(78.7).regime == 2
    @test_throws ErrorException face(78.8)
    @test_throws ErrorException face(80.0)
    # Thin: every admissible crossing (<= 75.5 degrees) is regime I under the 644.6-nm cap
    # (the boundary lies at 84.3 degrees), bitwise with the uncapped rule.
    for theta in (45.0, 70.0, 74.3, 75.5)
        f = face(theta; cap=cap_thin, env=0.078, pyr=0.016)
        legacy = face(theta; cap=Inf, env=0.078, pyr=0.016)
        @test f.regime == 1 && f.spacing == legacy.spacing && f.layers == legacy.layers
    end
    @test rad2deg(atan(cap_thin / (4.0 * 0.016))) ≈ 84.3 atol = 0.05
    @test_throws ErrorException face(45.0; cap=0.0)
end

@testset "face-ended tube stations: the sheared end block ends on the face plane" begin
    # The tube axis e = n x b runs towards -x at 40 degrees off the x axis; it reaches the
    # face x = -1 (normal N = (-1, 0, 0)) at s_face = 1 / cos 40; in tube coordinates the
    # face plane is s = s_face - u kappa_u with kappa_u = (N . n) / (N . e) = -tan 40.
    b = [0.0, 0.0, 1.0]
    n = [sind(40.0), -cosd(40.0), 0.0]
    e = cross(n, b)
    @test e ≈ [-cosd(40.0), -sind(40.0), 0.0]
    N = [-1.0, 0.0, 0.0]
    kappa = dot(N, n) / dot(N, e)
    @test kappa ≈ -tand(40.0)
    s_face = 1.0 / cosd(40.0)
    # r_env 0.04, h_pyr 0.01, lc_tangent 0.1 at 40 degrees: lc_end 0.1, m = ceil(0.0671 / 0.1) = 1.
    face_end =
        FaceEnd(1, "x0", N, deg2rad(40.0), kappa, 0.0, 0.04, 0.01, 0.1; spacing_cap=Inf)
    @test face_end.layers == 1 && face_end.spacing == 0.1
    tube = EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, s_face, 0.1; face_ends=[face_end])
    @test tube_cad_interval(tube) == (0.0, s_face + face_end.over_length)
    @test tube_end_station(tube, 1, 0.03, 0.0) ≈ s_face + 0.03 * tand(40.0) &&
          tube_end_station(tube, 0, 0.03, 0.0) == 0.0
    stations, shear_u, shear_w, _, _ = face_ended_tube_stations(
        tube,
        s -> 0.1,
        0.1,
        2.0;
        surface_start=true,
        surface_end=true
    )
    graded = EdgeTube(tube, stations; shear_u=shear_u, shear_w=shear_w)
    @test graded.stations[1] == 0.0 &&
          graded.stations[end] == s_face &&
          all(diff(graded.stations) .> 0.0)
    @test shear_u[end] == kappa && shear_w[end] == 0.0 && all(shear_u[1:(end - 1)] .== 0.0)
    L = graded.layers
    # Every node of the last station lies on the face x = -1; the block's layer thickness at
    # the outer ring (u = +-0.03) is lc_end + u tan 40, within [lc_end / 2, 3 lc_end / 2].
    for (u, w) in ((0.03, 0.0), (-0.03, 0.0), (0.0, 0.03), (0.02, 0.02))
        point = tube_point(graded, u, w, tube_station(graded, L, u, w))
        @test point[1] ≈ -1.0 atol = 1.0e-12
        thickness = tube_station(graded, L, u, w) - tube_station(graded, L - 1, u, w)
        @test thickness ≈ 0.1 + u * tand(40.0)
        @test 0.5 * 0.1 <= thickness <= 1.5 * 0.1
    end
    # The apex station is the layer midpoint at the apex's own (u, w).
    @test tube_station(graded, L - 0.5, 0.03, 0.0) ≈
          0.5 *
          (tube_station(graded, L - 1, 0.03, 0.0) + tube_station(graded, L, 0.03, 0.0))
    # A start face end on the plane x = 0 through the origin (N = (1, 0, 0): the same
    # kappa) mirrors the construction, here with a 70-degree face end (m = 2).
    start_face = FaceEnd(
        0,
        "x1",
        [1.0, 0.0, 0.0],
        deg2rad(70.0),
        -tand(70.0),
        0.0,
        0.04,
        0.01,
        0.1;
        spacing_cap=Inf
    )
    @test start_face.layers == 2 && start_face.spacing ≈ 0.04 * tand(70.0)
    both =
        EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, s_face, 0.1; face_ends=[face_end, start_face])
    stations, shear_u, shear_w, _, _ = face_ended_tube_stations(both, s -> 0.1, 0.1, 2.0)
    @test shear_u[1] == start_face.kappa_u &&
          shear_u[2] ≈ 0.5 * start_face.kappa_u &&
          shear_u[3] == 0.0
    @test stations[1] == 0.0 && stations[3] ≈ 2 * start_face.spacing
    graded = EdgeTube(both, stations; shear_u=shear_u, shear_w=shear_w)
    # With the geometric kappa of the plane x = 0 (-tan 40) the start nodes lie on it; the
    # prescribed 70-degree kappa only exercises the block arithmetic here.
    plane_face = FaceEnd(
        0,
        "x1",
        [1.0, 0.0, 0.0],
        deg2rad(40.0),
        kappa,
        0.0,
        0.04,
        0.01,
        0.1;
        spacing_cap=Inf
    )
    planar =
        EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, s_face, 0.1; face_ends=[face_end, plane_face])
    stations, shear_u, shear_w, _, _ = face_ended_tube_stations(planar, s -> 0.1, 0.1, 2.0)
    planar = EdgeTube(planar, stations; shear_u=shear_u, shear_w=shear_w)
    for (u, w) in ((0.03, 0.0), (-0.03, 0.01))
        @test tube_point(planar, u, w, tube_station(planar, 0, u, w))[1] ≈ 0.0 atol =
            1.0e-12
        @test tube_point(planar, u, w, tube_station(planar, planar.layers, u, w))[1] ≈ -1.0 atol =
            1.0e-12
    end
    # A tube shorter than its end blocks fails closed.
    short =
        EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, 0.25, 0.1; face_ends=[face_end, start_face])
    @test occursin(
        "shorter than its face-end blocks",
        guard_message(() -> face_ended_tube_stations(short, s -> 0.1, 0.1, 2.0))
    )
    # A plain tube: no shear, the (u, w) station is the axis station bitwise, the CAD
    # interval is the tube interval.
    plain = EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, 1.0, 0.1)
    @test isempty(plain.face_ends) && isempty(plain.shear_u)
    @test tube_station(plain, 3, 0.03, 0.02) === tube_station(plain, 3) &&
          tube_station(plain, 2.5, -0.03, 0.0) === tube_station(plain, 2.5)
    @test tube_cad_interval(plain) == (0.0, 1.0)
    @test_throws ErrorException EdgeTube(
        plain,
        plain.stations;
        shear_u=zeros(length(plain.stations)),
        shear_w=zeros(length(plain.stations))
    )
    @test_throws ErrorException EdgeTube(tube, tube.stations)
    @test_throws ErrorException EdgeTube(
        [0.0, 0.0, 0.0],
        n,
        b,
        0.0,
        1.0,
        0.1;
        face_ends=[face_end, face_end]
    )
end

@testset "face-ended tube entities: centroids of the trimmed solid match the OCC ones" begin
    section = TubeSection(0.01, 2.0, 2, [-90.0 + 30.0 * j for j = 0:9], fill(2, 9))
    b = [0.0, 0.0, 1.0]
    # A frame whose extrusion e = n x b runs obliquely towards -x.
    n = [sin(deg2rad(40.0)), -cos(deg2rad(40.0)), 0.0]
    e = cross(n, b)
    @test e[1] < 0.0
    # The tube runs towards the face x = -1 (normal (-1, 0, 0)); kappa_u = (N . n) / (N . e).
    N = [-1.0, 0.0, 0.0]
    theta = acos(abs(dot(N, e)))
    face_end = FaceEnd(
        1,
        "x0",
        N,
        theta,
        dot(N, n) / dot(N, e),
        dot(N, b) / dot(N, e),
        0.04,
        0.01,
        0.1;
        spacing_cap=Inf
    )
    s_end = (1.0 - 0.0) / abs(e[1])     # the axis from (0, 0, 0) reaches x = -1 at this s
    tube = EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, s_end, 0.1; face_ends=[face_end])
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("face-ended-entities")
    occ = gmsh.model.occ
    box = ([-1.0, -2.0, -1.0], [2.0, 2.0, 1.0])
    volumes = add_tube_volumes!(occ, tube, section; box=box)
    @test length(volumes) == 1
    occ.synchronize()
    (volume, group) = volumes[1]
    matched = match_tube_entities(volume[2], tube, section, group, 1.0e-9)
    entities = tube_entities(tube, section, group)
    @test length(matched) == length(entities)
    # Every end-cap entity of the face end lies on the face x = -1; the start cap at s = 0.
    for entity in entities
        entity.id[2] == 1 &&
        entity.kind in (:edge_point, :outer_point, :cap_polygon, :cap_ray, :cap) || continue
        @test entity.centroid[1] ≈ -1.0 atol = 1.0e-12
    end
    bbox = gmsh.model.getBoundingBox(3, volume[2])
    @test bbox[1] >= -1.0 - 1.0e-6       # nothing protrudes past the face (OCC pads its boxes ~1e-7)
    # A plain tube needs no box; a face-ended one fails closed without it.
    @test_throws ErrorException add_tube_volumes!(occ, tube, section)
    plain = EdgeTube([0.0, 0.0, 0.0], n, b, 0.0, 0.5, 0.1)
    @test length(add_tube_volumes!(occ, plain, section)) == 1
    gmsh.finalize()
    @test polygon_centroid_3d([
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [2.0, 1.0, 0.0],
        [0.0, 1.0, 0.0]
    ]) ≈ [1.0, 0.5, 0.0]
end

@testset "metal edge segments: box-face ends by the VERTEX, angle-gated (decision 320)" begin
    # A metal wedge whose two sides leave the box [-4, 4]^2: side A along +x (theta 0,
    # the legacy perpendicular end), side B at 30 degrees off the +y normal (theta 30),
    # meeting at the tip (0, 0); the box sides close the loop (Continuation).
    tip = (0.0, 0.0)
    a_exit = (4.0, 0.0)
    b_dir = (-sin(deg2rad(30.0)), cos(deg2rad(30.0)))
    b_exit = (4.0 * b_dir[1] / b_dir[2], 4.0)
    points = [tip, a_exit, (4.0, 4.0), b_exit]
    classes = ["Physical", "Continuation", "Continuation", "Physical"]
    loop = (conductor=1, plane=0.0, hole=false, points=points, classes=classes)
    lower = [-4.0, -4.0]
    upper = [4.0, 4.0]
    clearance(angle) = 0.03 / tan(0.5 * angle) + 0.02
    # The contract lists the tip only (b_exit is a box-face cut end: theta 30 > 0).
    segments =
        metal_edge_segments([loop], [(0.0, 0.0, 0.0)], clearance, lower, upper, 1.0e-9)
    @test length(segments) == 2
    side_a, side_b = segments
    @test side_a.face_ends[1] === nothing && side_a.face_ends[2] === nothing   # theta == 0 exactly
    @test side_a.s_end == 4.0 && side_a.legacy_box_corners == (false, false)
    @test side_b.face_ends[1] !== nothing && side_b.face_ends[2] === nothing
    @test side_b.face_ends[1].face == "y1" && side_b.face_ends[1].normal == [0.0, 1.0]
    @test side_b.face_ends[1].theta ≈ deg2rad(30.0)
    @test side_b.s_start == 0.0 && side_b.corner_angles[1] == Float64(pi)
    @test side_b.corner_angles[2] ≈ deg2rad(120.0)     # the tip angle between +x and the 120-degree side
    @test side_b.s_end ≈ side_b.span - clearance(deg2rad(120.0))
    # A contract that still lists the oblique box vertex as a corner fails closed.
    message = guard_message(
        () -> metal_edge_segments(
            [loop],
            [(0.0, 0.0, 0.0), (b_exit[1], b_exit[2], 0.0)],
            clearance,
            lower,
            upper,
            1.0e-9
        )
    )
    @test occursin("regenerate the contract", message) && occursin("decision 320", message)
    # The legacy convention at theta == 0: a Physical-class box vertex listed by the contract
    # keeps h_K and is recorded as a LegacyBoxVertexCorner (the mesh is unchanged).
    legacy_classes = ["Physical", "Physical", "Continuation", "Physical"]
    legacy_loop = (loop..., classes=legacy_classes)
    legacy = metal_edge_segments(
        [legacy_loop],
        [(0.0, 0.0, 0.0), (4.0, 0.0, 0.0)],
        clearance,
        lower,
        upper,
        1.0e-9
    )
    @test legacy[1].s_end ≈ 4.0 - clearance(Float64(pi)) &&
          legacy[1].legacy_box_corners == (false, true)
    @test legacy[1].face_ends == (nothing, nothing)
    # Two metal sides meeting at a box vertex stay a corner at any tilt: an island whose
    # two oblique sides meet at (4, 0) on the face x = 4.
    island = [(4.0, 0.0), (2.0, 1.0), (2.0, -1.0)]
    island_loop =
        (conductor=1, plane=0.0, hole=false, points=island, classes=fill("Physical", 3))
    island_segments = metal_edge_segments(
        [island_loop],
        [(p[1], p[2], 0.0) for p in island],
        clearance,
        lower,
        upper,
        1.0e-9
    )
    @test length(island_segments) == 3
    @test all(segment.face_ends == (nothing, nothing) for segment in island_segments)
    @test island_segments[1].start_corner &&
          island_segments[1].s_start ≈ clearance(2 * atan(0.5))
    @test island_segments[1].legacy_box_corners == (false, false)
    # An oblique side through a box corner has no single face plane: fail closed.
    @test occursin(
        "box corner obliquely",
        guard_message(() -> crossed_box_face([4.0, 4.0], [0.6, 0.8], lower, upper, 1.0e-9))
    )
    @test crossed_box_face([4.0, 4.0], [0.0, 1.0], lower, upper, 1.0e-9) == (2, 1)
    @test crossed_box_face([4.0, 1.0], [0.6, 0.8], lower, upper, 1.0e-9) == (1, 1)
    @test crossed_box_face([1.0, 1.0], [0.6, 0.8], lower, upper, 1.0e-9) === nothing
end

# A convex metal tip of angle phi (degrees) at the origin whose side A leaves the coupon
# box through its +x face at the tilt theta (degrees, 0 = perpendicular) and whose side B
# is phi degrees counterclockwise from side A: the V7-c synthetic family. The box is the
# mesher's own (the rows padded by Radius), the metal is the wedge clipped to it, the
# boundary classes are those of the plan-view builder (a vertex carries its outgoing
# side's class) and the contract corners follow decision 320 (the tip, plus the B exit
# when side B is exactly perpendicular to its face).
function write_tip_inputs(
    directory,
    phi,
    theta;
    radius=0.5,
    side_length=2.0,
    side_length_b=side_length,
    plane=0.0
)
    d_a = (cosd(theta), sind(theta))
    d_b = (cosd(theta + phi), sind(theta + phi))
    # Row B may be longer than row A (side_length_b): the box grows along side B, which
    # moves side A's exit onto another face or steepens its crossing.
    rows = [
        (
            point=(0.5 * L * d[1], 0.5 * L * d[2], plane),
            tangent=(d[1], d[2], 0.0),
            gap=(d[2], -d[1], 0.0) .* sign,
            interval=(-0.5 * L, 0.5 * L),
            normal_sign=1.0,
            vertex_arm=false,
            slot=0,
            conductor=1
        ) for (d, sign, L) in ((d_a, 1.0, side_length), (d_b, -1.0, side_length_b))
    ]
    lower, upper = row_coupon_bounds(rows, radius, 0.1, 0.05)
    # The wedge clipped to the box: start from the box rectangle (counterclockwise) and
    # clip by the half-planes left of side A and right of side B.
    polygon = [
        (lower[1], lower[2]),
        (upper[1], lower[2]),
        (upper[1], upper[2]),
        (lower[1], upper[2])
    ]
    function clip(polygon, inside)
        result = Tuple{Float64, Float64}[]
        m = length(polygon)
        for i = 1:m
            p = polygon[i]
            q = polygon[i % m + 1]
            fp, fq = inside(p), inside(q)
            fp >= 0.0 && push!(result, p)
            if (fp > 0.0) != (fq > 0.0) && fp != 0.0 && fq != 0.0
                t = fp / (fp - fq)
                push!(result, (p[1] + t * (q[1] - p[1]), p[2] + t * (q[2] - p[2])))
            end
        end
        return result
    end
    polygon = clip(polygon, p -> d_a[1] * p[2] - d_a[2] * p[1])     # left of side A
    polygon = clip(polygon, p -> d_b[2] * p[1] - d_b[1] * p[2])     # right of side B
    # Snap the exits onto the box faces exactly, quantise to the builder's 1e-9 R grid (so an
    # axis-aligned side reads exactly axis-aligned) and rotate the loop to start at the tip.
    snap(v) = (
        any(abs(v[1] - f) <= 1.0e-9 for f in (lower[1], upper[1])) ?
        (abs(v[1] - lower[1]) <= 1.0e-9 ? lower[1] : upper[1]) :
        round(v[1] / 1.0e-9) * 1.0e-9,
        any(abs(v[2] - f) <= 1.0e-9 for f in (lower[2], upper[2])) ?
        (abs(v[2] - lower[2]) <= 1.0e-9 ? lower[2] : upper[2]) :
        round(v[2] / 1.0e-9) * 1.0e-9
    )
    polygon = [snap(v) for v in polygon]
    k = argmin([hypot(v...) for v in polygon])
    polygon = vcat(polygon[k:end], polygon[1:(k - 1)])
    on_face(p, q) = any(
        (abs(p[d] - lower[d]) <= 1.0e-9 && abs(q[d] - lower[d]) <= 1.0e-9) ||
            (abs(p[d] - upper[d]) <= 1.0e-9 && abs(q[d] - upper[d]) <= 1.0e-9) for
        d = 1:2
    )
    m = length(polygon)
    classes =
        [on_face(polygon[i], polygon[i % m + 1]) ? "Continuation" : "Physical" for i = 1:m]
    @assert count(==("Physical"), classes) == 2
    # Contract corners (decision 320): a Physical vertex after a Continuation side is a
    # corner only when its side is exactly perpendicular to that face.
    corners = Tuple{Float64, Float64, Float64}[]
    for i = 1:m
        classes[i] == "Physical" || continue
        previous = polygon[mod1(i - 1, m)]
        if classes[mod1(i - 1, m)] == "Continuation"
            constant = findfirst(d -> previous[d] == polygon[i][d], 1:2)
            polygon[i][3 - constant] == polygon[i % m + 1][3 - constant] || continue
        end
        push!(corners, (polygon[i][1], polygon[i][2], plane))
    end
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
        println(io, "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y")
        for (i, point) in enumerate(polygon)
            println(io, join([1, i, 1, plane, 0, classes[i], point[1], point[2]], ","))
        end
    end
    open(joinpath(directory, "mask.csv"), "w") do io
        println(io, "Facet,Conductor,Plane,X,Y")
        for point in polygon
            println(io, join([1, 1, plane, point[1], point[2]], ","))
        end
    end
    # Design round 2 F5-A (decision 363): the contract records the INVARIANT corners (the
    # two sides' quantised dot product non-zero: the tip unless phi is exactly 90 with the
    # same arm coordinates; never the theta-0 box vertex) under Derivation.InvariantCorners,
    # as derive_semantic_contract does; the mesher fails closed on a disagreement.
    invariant = Vector{Float64}[]
    for c in corners
        i = findfirst(v -> v[1] == c[1] && v[2] == c[2], polygon)
        a = polygon[mod1(i - 1, m)] .- polygon[i]
        b = polygon[mod1(i + 1, m)] .- polygon[i]
        a[1] * b[1] + a[2] * b[2] == 0.0 || push!(invariant, [c[1], c[2], c[3]])
    end
    contract = Dict{String, Any}(
        "Version" => 1,
        "SemanticCorners" => [[c[1], c[2], c[3]] for c in corners]
    )
    isempty(invariant) || (
        contract["Derivation"] = Dict{String, Any}(
            "InvariantCorners" => Dict{String, Any}(
                "Rule" => "test fixture: the exact dot-product predicate",
                "Points" => invariant
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
        invariant=invariant,
        lower=lower,
        upper=upper
    )
end

function build_tip_coupon(
    directory,
    phi,
    theta;
    fabricated=true,
    stem="tip",
    lc_tangent=0.1
)
    inputs = write_tip_inputs(directory, phi, theta)
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
        lc_tangent=lc_tangent,
        lc_far=0.3,
        max_nodes=2_000_000,
        max_elements=2_000_000,
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
        # The invariant-corner verdict bound (design round 2 F5-A): the tip is invariant.
        corner_shape_gate=5.0
    )
    return parse_json(read(census, String)), mesh, inputs
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

# The tip family at the PRODUCTION tube sizes (EdgeSize 0.25 nm fabricated / 2 nm thin, ratio 2,
# TangentialSize 50 nm, NormalSize 25 nm, CornerIsotropyRadius 0.1 um; the step-4.1 V7-c and the
# round-2 V10 families) on a small coupon (Radius 0.5, FarSize 0.152 as in production so that the
# tangential size stays 50 nm), with a long row B (`factor` x) so that side A meets the +x face
# at the steep crossing theta.
function build_production_tip_coupon(
    directory,
    phi,
    theta,
    factor;
    fabricated=true,
    stem="prod",
    labels_only=false
)
    inputs = write_tip_inputs(
        directory,
        phi,
        theta;
        radius=0.5,
        side_length=0.5,
        side_length_b=factor * 0.5
    )
    mesh = joinpath(directory, "coupon-$stem.msh")
    census = joinpath(directory, "census-$stem.json")
    edge = fabricated ? 0.00025 : 0.002
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
        lc_fine=0.025,
        lc_tangent=0.05,
        lc_far=0.152,
        max_nodes=2_000_000,
        max_elements=2_000_000,
        semantic_contract=inputs.semantic,
        corner_isotropy_radius=0.1,
        corner_census=census,
        edge_size=edge,
        edge_growth_ratio=2.0,
        corner_size=edge,
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

@testset "F2b full builds at the production tube sizes: the 70-degree crossing in regime II (fabricated), regime I (thin), the 80-degree guard" begin
    mktempdir() do directory
        # phi 30, side A at 70 degrees, row B 6x: side A meets the +x face at exactly 70 degrees
        # (the V7-c steep set whose fabricated twin read a prism condition 1045 > 1000 under the
        # uncapped rule, phase2 REPORT 4.1), side B leaves through +y at 10 degrees.
        census, mesh, inputs =
            build_production_tip_coupon(directory, 30.0, 70.0, 6.0; stem="fab70")
        tubes = census["PrismTubes"]
        section = tubes["Section"]
        cap = face_end_spacing_cap(
            TubeSection(
                0.00025,
                2.0,
                section["Rings"],
                section["Top"]["Angles"],
                section["Top"]["Materials"]
            ),
            1000.0
        )
        @test section["FaceEndSpacingCap"] == cap
        @test cap ≈ 0.08057212272716054 rtol = 1.0e-12
        @test occursin("regime II", tubes["FaceEnds"]["Rule"]) &&
              occursin("SteepFaceCrossing", tubes["FaceEnds"]["Rule"])
        records = [f for row in tubes["Tubes"] for f in get(row, "FaceEnds", [])]
        steep = [f for f in records if f["ThetaDegrees"] ≈ 70.0]
        shallow = [f for f in records if f["ThetaDegrees"] ≈ 10.0]
        @test length(steep) == 2 && length(shallow) == 2 && length(records) == 4
        for f in steep
            @test f["Regime"] == "II" &&
                  f["EndSpacing"] == cap &&
                  f["Layers"] == 3 &&
                  f["EndSpacingCap"] == cap
            @test f["ApexThickness"] ≈ 2.0 * section["PyramidHeight"] * tand(70.0)
            @test f["LayerThicknessRange"][1] >= f["ApexThickness"]
            @test f["LayerThicknessRange"][1] ≈ 0.04417 atol = 2.0e-5
        end
        for f in shallow
            @test f["Regime"] == "I" && f["EndSpacing"] == 0.05 && f["Layers"] == 1
        end
        # B14's measured clauses: the maximum prism condition in [900, 1000] (the capped block
        # at 0.95 x 1000 on the analytic frame, the sheared block ~1 % above it), min SJ >= 0.01,
        # every gate passed (the build would have failed otherwise).
        quality = tubes["Quality"]
        @test 900.0 <= quality["Prism"]["MaximumJacobianCondition"] <= 1000.0
        @test quality["Prism"]["MaximumJacobianCondition"] > section_prism_condition(
            TubeSection(
                0.00025,
                2.0,
                section["Rings"],
                section["Top"]["Angles"],
                section["Top"]["Materials"]
            ),
            tubes["TangentialSize"]
        )
        @test quality["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        @test quality["Prism"]["PositiveOrientation"] &&
              quality["Pyramid"]["PositiveOrientation"]
        @test tubes["FaceEndSpacingMaximum"] == cap
        # The thin twin: every thin crossing is regime I (its cap 644.6 nm), lc_end = 4 h_pyr
        # tan 70 = 175.84 nm, m 3; its prisms read ~268 (the design's thin figure).
        thin, _, _ = build_production_tip_coupon(
            directory,
            30.0,
            70.0,
            6.0;
            fabricated=false,
            stem="thin70"
        )
        thin_tubes = thin["PrismTubes"]
        @test thin_tubes["Section"]["FaceEndSpacingCap"] ≈ 0.6445769818172843 rtol = 1.0e-12
        thin_steep = [
            f for row in thin_tubes["Tubes"] for
            f in get(row, "FaceEnds", []) if f["ThetaDegrees"] ≈ 70.0
        ]
        @test length(thin_steep) == 1 &&
              thin_steep[1]["Regime"] == "I" &&
              thin_steep[1]["Layers"] == 3
        @test thin_steep[1]["EndSpacing"] ≈
              4.0 * thin_tubes["Section"]["PyramidHeight"] * tand(70.0)
        @test 250.0 <= thin_tubes["Quality"]["Prism"]["MaximumJacobianCondition"] <= 300.0
        @test thin_tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        # 80 degrees (row B 20x so that side A reaches the +x face): the fabricated tube end lies
        # beyond the validity ceiling 78.77 degrees and fails closed at the guard before any
        # CAD tube exists; the thin twin's 80-degree end is regime I (its ceiling 87.2 degrees).
        message = guard_message(
            () -> build_production_tip_coupon(
                directory,
                30.0,
                80.0,
                20.0;
                stem="fab80",
                labels_only=true
            )
        )
        @test occursin("ScopeGuard[SteepFaceCrossing]", message) &&
              occursin("80.0", message) &&
              occursin("validity ceiling", message)
        # Round 3 class (8) interim 8B (decision 510; DESIGN-part-M 5.2): the thin 80-degree end
        # (regime I under its 87.2-degree ceiling, but above the largest BUILT thin tilt 70 degrees:
        # the thin 74.3 / 75.5-degree crossings fail the tetrahedral gate after the build, E3) now
        # fails closed at the same guard before any CAD tube, naming the tested-range rule.
        @test THIN_FACE_END_TILT_BOUND == deg2rad(70.0)
        message = guard_message(
            () -> build_production_tip_coupon(
                directory,
                30.0,
                80.0,
                20.0;
                fabricated=false,
                stem="thin80",
                labels_only=true
            )
        )
        @test occursin("ScopeGuard[SteepFaceCrossing]", message) &&
              occursin("THIN tube end", message) &&
              occursin("80.0", message) &&
              occursin("largest built thin tilt 70.0", message)
    end
end

@testset "synthetic 45-degree tip, side A at theta 20: fabricated and thin builds end on both faces" begin
    mktempdir() do directory
        for (fabricated, stem) in ((true, "fab"), (false, "thin"))
            census, mesh, inputs =
                build_tip_coupon(directory, 45.0, 20.0; fabricated=fabricated, stem=stem)
            tubes = census["PrismTubes"]
            per_side = tubes["Section"]["TubesPerSide"]
            # Side A exits the +x face at 20 degrees, side B (at 65 degrees) the +y face at a
            # 25-degree tilt: both exits are box-face cut ends, the tip is the only corner.
            @test length(inputs.corners) == 1 && inputs.corners[1] == (0.0, 0.0, 0.0)
            face_ended = [row for row in tubes["Tubes"] if haskey(row, "FaceEnds")]
            @test length(face_ended) == 2 * per_side == length(tubes["Tubes"])
            @test tubes["FaceEnds"]["Count"] == length(face_ended)
            @test tubes["FaceEnds"]["LegacyBoxVertexCorners"] == 0
            @test tubes["FaceEnds"]["EndBlockLayers"] == length(face_ended)
            for row in face_ended
                @test length(row["FaceEnds"]) == 1
                record = row["FaceEnds"][1]
                d = record["Face"] == "x1" ? 1 : 2
                @test record["Face"] in ("x1", "y1")
                @test record["ThetaDegrees"] ≈ (d == 1 ? 20.0 : 25.0)
                @test record["Layers"] == 1 && record["EndSpacing"] == 0.1
                # The axis end lies on the face and the end polygon's nodes with it; nothing of
                # the tube protrudes past the face.
                end_point = record["End"] == "end" ? row["EndPoint"] : row["StartPoint"]
                @test end_point[d] ≈ inputs.upper[d] atol = 1.0e-9
                near = nodes_near(mesh, end_point, 2.0 * tubes["Section"]["Radius"])
                on_face = count(p -> abs(p[d] - inputs.upper[d]) <= 1.0e-9, near)
                @test on_face >= 1 + tubes["Section"]["Rings"] * (fabricated ? 10 : 12)
                @test all(p[d] <= inputs.upper[d] + 1.0e-9 for p in near)
            end
            # Every tube volume is a single fragment descendant (the build would have failed
            # otherwise); every element positively oriented; the quality gates passed.
            quality = tubes["Quality"]
            @test quality["Prism"]["PositiveOrientation"] &&
                  quality["Pyramid"]["PositiveOrientation"]
            @test quality["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
            @test tubes["Prisms"] == sum(row["Prisms"] for row in tubes["Volumes"])
            @test occursin("decision 320", tubes["Section"]["BoxVertexRule"])
            @test occursin("sin(phi / 2)", tubes["Section"]["CornerClearanceRule"])
        end
    end
end

@testset "synthetic 45-degree tip at theta 45 under TangentialSize 0.05 (m = 2 blocks) beside a theta-0 legacy box corner" begin
    mktempdir() do directory
        census, mesh, inputs =
            build_tip_coupon(directory, 45.0, 45.0; stem="sharp", lc_tangent=0.05)
        tubes = census["PrismTubes"]
        face_ended = [row for row in tubes["Tubes"] if haskey(row, "FaceEnds")]
        # Side A leaves through +x at 45 degrees (m = ceil(2 x 0.04 / 0.05) = 2 at lc_end 0.05);
        # side B along +y leaves through +y exactly perpendicularly: a Physical-class box vertex
        # kept as a legacy corner (h_K + ball) on both of its tubes.
        @test length(face_ended) == 2
        for row in face_ended
            record = row["FaceEnds"][1]
            @test record["Face"] == "x1" && record["ThetaDegrees"] ≈ 45.0
            @test record["Layers"] == 2 && record["EndSpacing"] == 0.05
            @test 0.5 * record["EndSpacing"] <= record["LayerThicknessRange"][1] &&
                  record["LayerThicknessRange"][2] <=
                  1.5 * record["EndSpacing"] * (1.0 + 1.0e-9)
            end_point = record["End"] == "end" ? row["EndPoint"] : row["StartPoint"]
            @test end_point[1] ≈ inputs.upper[1] atol = 1.0e-9
            near = nodes_near(mesh, end_point, 2.0 * tubes["Section"]["Radius"])
            @test all(p[1] <= inputs.upper[1] + 1.0e-9 for p in near)
            @test count(p -> abs(p[1] - inputs.upper[1]) <= 1.0e-9, near) >=
                  1 + tubes["Section"]["Rings"] * 10
        end
        @test tubes["FaceEnds"]["EndBlockLayers"] == 4 &&
              tubes["FaceEnds"]["LegacyBoxVertexCorners"] == 2
        @test length(inputs.corners) == 2 && (0.0, inputs.upper[2], 0.0) in inputs.corners
        legacy_rows =
            [row for row in tubes["Tubes"] if haskey(row, "LegacyBoxVertexCorner")]
        @test length(legacy_rows) == 2 && all(
            row["LegacyBoxVertexCorner"] == [[0.0, inputs.upper[2]]] for row in legacy_rows
        )
        # The 45-degree tip clearance carries the pyramid-envelope margin (design A3 (3)):
        # R / tan(22.5) + max(h_K, 1.25 h_pyr / sin(22.5)) - the envelope term wins below 77.4
        # degrees; the tubes ending at the tip start there.
        radius, h_pyr = tubes["Section"]["Radius"], tubes["Section"]["PyramidHeight"]
        h_k = tubes["Section"]["RingSizes"][end]
        expected = radius / tand(22.5) + max(h_k, 1.25 * h_pyr / sind(22.5))
        @test expected > radius / tand(22.5) + h_k
        @test any(abs(row["Start"] - expected) <= 1.0e-9 for row in tubes["Tubes"])
        @test all(
            any(abs(angle - deg2rad(45.0)) <= 1.0e-9 for angle in row["CornerAngles"]) for
            row in tubes["Tubes"]
        )
        @test tubes["Quality"]["Tetrahedron"]["MinimumScaledJacobian"] >= 0.01
        @test tubes["Quality"]["Prism"]["PositiveOrientation"] &&
              tubes["Quality"]["Pyramid"]["PositiveOrientation"]
        # theta 0: side A along +x (the legacy perpendicular end, no FaceEnd), side B at 90
        # degrees along +y: no face end anywhere, the census carries Count 0; side B's exit is a
        # Physical-class box vertex whose side leaves perpendicularly: the legacy convention
        # keeps it a corner with h_K (recorded on the side's two tubes).
        legacy, _, legacy_inputs = build_tip_coupon(directory, 90.0, 0.0; stem="legacy")
        @test legacy["PrismTubes"]["FaceEnds"]["Count"] == 0
        @test all(!haskey(row, "FaceEnds") for row in legacy["PrismTubes"]["Tubes"])
        @test length(legacy_inputs.corners) == 2
        @test legacy["PrismTubes"]["FaceEnds"]["LegacyBoxVertexCorners"] == 2
    end
end
