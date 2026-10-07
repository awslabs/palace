# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Localized prism edge tubes of the Gmsh-only spatial coupon mesher (supervisor
# decisions 37 and 38; first built as the prism-tube feasibility spike).
#
# A straight metal edge is surrounded, on its dielectric side, by a tube whose 2D
# cross-section is meshed once with geometric rings (ring k has radial size
# inner_size x ratio^(k - 1)) and extruded along the edge in layers whose
# thickness follows the composed size field along the edge (at most lc_tangent;
# supervisor decision 40), giving prisms. The tube volumes are OCC polygon-sector prisms fragmented with
# the coupon CAD; their mesh (points, curves, faces, volumes) is installed
# explicitly so that Mesh.MeshOnlyEmpty leaves it alone while Gmsh meshes the
# remaining entities. The tube's outer lateral quadrangles carry explicit
# pyramids, so Gmsh meshes pure tetrahedra against triangles.

import Gmsh: gmsh
using LinearAlgebra

# Angular rays are given in degrees measured from the outward horizontal normal
# (n, the signature "gap" direction) towards +z (b). The extrusion direction is
# e = n x b so that the cross-section triangles (counterclockwise in (n, b)) make
# positively oriented prisms.
struct TubeSection
    ring_radii::Vector{Float64}      # cumulative radii r_1 < ... < r_K (um), r_0 = 0 is the edge
    angles::Vector{Float64}          # rays theta_0 < ... < theta_J (degrees), J sectors
    materials::Vector{Int}           # material attribute of each sector (1 substrate, 2 vacuum)
    closed::Bool                     # theta_J = theta_0 + 360: ray J is ray 0 (a full-turn sheet section)
end

# A section is closed when its rays span a full turn (the thin sheet edge: the metal
# sheet ray is both the first and the last ray); its last ray shares the first ray's
# nodes and CAD entities.
function TubeSection(inner_size, ratio, rings, angles, materials)
    inner_size > 0.0 || error("tube inner size must be positive")
    ratio > 1.0 || error("tube ring ratio must exceed 1")
    rings >= 1 || error("at least one ring")
    length(materials) == length(angles) - 1 || error("one material per angular sector")
    all(diff(angles) .> 0.0) || error("rays must increase")
    span = angles[end] - angles[1]
    span <= 360.0 + 1.0e-9 || error("tube rays span more than a full turn")
    radii = [inner_size * (ratio^k - 1.0) / (ratio - 1.0) for k = 1:rings]
    return TubeSection(
        radii,
        collect(Float64, angles),
        collect(Int, materials),
        abs(span - 360.0) <= 1.0e-9
    )
end

ring_sizes(section::TubeSection) = diff(vcat(0.0, section.ring_radii))
tube_radius(section::TubeSection) = section.ring_radii[end]
ring_count(section::TubeSection) = length(section.ring_radii)
ray_count(section::TubeSection) = length(section.angles)
# Rays with their own nodes: every ray of an open section, all but the last of a closed one.
distinct_ray_count(section::TubeSection) = ray_count(section) - (section.closed ? 1 : 0)

# Cross-section node (k, j): k = 0 the edge point (one node), k >= 1 ring k, ray j
# (0-based; the last ray of a closed section folds onto ray 0). Returns the 1-based
# local index within one cross-section.
function section_node(section::TubeSection, k, j)
    k == 0 && return 1
    j = section.closed ? mod(j, distinct_ray_count(section)) : j
    return 1 + (k - 1) * distinct_ray_count(section) + j + 1
end
section_node_count(section::TubeSection) =
    1 + ring_count(section) * distinct_ray_count(section)

# Local (u, w) coordinates of every cross-section node, u along n, w along b.
function section_coordinates(section::TubeSection)
    uw = zeros(2, section_node_count(section))
    for k = 1:ring_count(section), j = 0:(distinct_ray_count(section) - 1)
        theta = deg2rad(section.angles[j + 1])
        uw[:, section_node(section, k, j)] .=
            section.ring_radii[k] .* (cos(theta), sin(theta))
    end
    return uw
end

# Counterclockwise triangles of the cross-section per sector j (0-based):
# the fan of ring 1 and two triangles per annulus quad.
function sector_triangles(section::TubeSection, j)
    triangles = NTuple{3, Int}[]
    push!(
        triangles,
        (
            section_node(section, 0, 0),
            section_node(section, 1, j),
            section_node(section, 1, j + 1)
        )
    )
    for k = 2:ring_count(section)
        a = section_node(section, k - 1, j)
        b = section_node(section, k, j)
        c = section_node(section, k, j + 1)
        d = section_node(section, k - 1, j + 1)
        push!(triangles, (a, b, c))
        push!(triangles, (a, c, d))
    end
    return triangles
end

# Contiguous sectors of one material: (first sector, last sector, material), 0-based.
function material_groups(section::TubeSection)
    groups = Tuple{Int, Int, Int}[]
    start = 0
    for j = 1:(length(section.materials) - 1)
        if section.materials[j + 1] != section.materials[j]
            push!(groups, (start, j - 1, section.materials[j]))
            start = j
        end
    end
    push!(groups, (start, length(section.materials) - 1, section.materials[end]))
    return groups
end

# Rays that are CAD curves of the cross-section: the two metal-face rays and every
# material interface ray.
function cad_rays(section::TubeSection)
    rays = [0, ray_count(section) - 1]
    for (first, _, _) in material_groups(section)[2:end]
        push!(rays, first)
    end
    return sort!(unique!(rays))
end

# The planar corner frames of a cross-section's prisms (mesher design round 2 F2b,
# section 3.2): every sector triangle of the section (sector_triangles) seen from each
# of its three vertices, the two triangle edges leaving the vertex as the columns of a
# 2 x 2 matrix; a prism layer of axial spacing lc adds the orthogonal third column lc e,
# so its corner frame (VOLUME_CORNER_FRAMES, mixed_volume_quality) has the singular
# values {lc, sigma_1, sigma_2} of this planar block. Returns (sigma_1, sigma_2) per frame
# (closed form for a 2 x 2 matrix). The smallest sigma_2 of the production section is
# the second ring's (a, c, d) triangle at its ring-2 vertex c (0.33925 r_1 at ratio 2 and
# 30-degree sectors), not the inner triangle's edge-point frame (0.36603 r_1).
function section_frame_singular_values(section::TubeSection)
    radii = vcat(0.0, section.ring_radii)
    node(k, j) = radii[k + 1] .* (cosd(section.angles[j + 1]), sind(section.angles[j + 1]))
    frames = Tuple{Float64, Float64}[]
    for j = 0:(length(section.angles) - 2)
        for (a, b, c) in sector_triangles(section, j)
            # sector_triangles indexes the section nodes; their polar coordinates follow
            # from the ring and ray of each local index.
            triangle = [
                node(section_ring_ray(section, local_index, j)...) for
                local_index in (a, b, c)
            ]
            for v = 1:3
                p = triangle[v]
                e1 = triangle[mod1(v + 1, 3)] .- p
                e2 = triangle[mod1(v + 2, 3)] .- p
                n1 = e1[1]^2 + e1[2]^2
                n2 = e2[1]^2 + e2[2]^2
                cross = e1[1] * e2[1] + e1[2] * e2[2]
                half = 0.5 * (n1 + n2)
                root = 0.5 * sqrt((n1 - n2)^2 + 4.0 * cross^2)
                push!(frames, (sqrt(half + root), sqrt(max(half - root, 0.0))))
            end
        end
    end
    return frames
end

# The (ring, ray) of a local section node index on sector j (the inverse of section_node
# restricted to the two rays of one sector; the edge point is ring 0).
function section_ring_ray(section::TubeSection, local_index, j)
    local_index == 1 && return (0, 0)
    for k = 1:ring_count(section), ray in (j, j + 1)
        section_node(section, k, ray) == local_index && return (k, ray)
    end
    return error("section node $local_index does not lie on sector $j")
end

# The largest Jacobian condition over the prism corner frames of one regular layer of
# axial spacing `spacing` of this section: max(spacing, sigma_1) / min(spacing, sigma_2)
# over the planar frames - a closed-form function of the section, non-decreasing in the
# spacing once it exceeds every sigma_1 (the production tubes: 589.33 at 49.98 nm on the
# fabricated 0.25-nm section, 73.52 at 49.88 nm on the thin 2-nm section, equal to the
# stored S2p 3-edge censuses' Prism.MaximumJacobianCondition to 1e-12 relative).
function section_prism_condition(section::TubeSection, spacing)
    spacing > 0.0 || error("a prism layer needs a positive spacing")
    return maximum(
        max(spacing, s1) / min(spacing, s2) for
        (s1, s2) in section_frame_singular_values(section)
    )
end

# The end-spacing cap of a face end (design round 2 F2b 3.2): the largest axial spacing
# at which the section's own prism frames read <= FACE_END_CONDITION_MARGIN x the
# Jacobian-condition ceiling, lc_cap = margin x ceiling x min sigma_2 over the frames
# (80.572 nm on the fabricated production section at the ceiling 1000, 644.58 nm thin).
# Fails closed when the cap does not lie in the spacing-dominated regime (some frame's
# sigma_1 above it), where the condition would not be monotone in the spacing.
const FACE_END_CONDITION_MARGIN = 0.95

function face_end_spacing_cap(section::TubeSection, maximum_jacobian_condition)
    isfinite(maximum_jacobian_condition) && maximum_jacobian_condition > 1.0 || error(
        "a face end's end-spacing cap needs a finite Jacobian-condition ceiling > 1 " *
        "(--maximum-jacobian-condition), got $maximum_jacobian_condition"
    )
    ceiling = FACE_END_CONDITION_MARGIN * maximum_jacobian_condition
    frames = section_frame_singular_values(section)
    cap = ceiling * minimum(s2 for (_, s2) in frames)
    all(s1 <= cap for (s1, _) in frames) || error(
        "the end-spacing cap $cap of the tube section lies below a prism frame's " *
        "largest planar singular value: the prism condition is not spacing-dominated"
    )
    section_prism_condition(section, cap) <= ceiling * (1.0 + 1.0e-12) ||
        error("the end-spacing cap $cap does not meet the condition ceiling $ceiling")
    return cap
end

# A tube end on a face of the coupon box that the edge does not cross
# perpendicularly (block (b) design AMENDMENT 1 A2, supervisor decisions 302 / 320):
# the tube ENDS ON THE FACE. The face plane through the axis point at the end
# station reads, in tube coordinates, s_face(u, w) = s_end - u kappa_u - w kappa_w
# with kappa_u = (N . n) / (N . e), kappa_w = (N . b) / (N . e) for the outward face
# normal N (|kappa| = tan theta for a vertical face and a horizontal tube); theta is
# the angle between the tube axis and N. The CAD solid is extruded over-long by
# over_length = (radius + pyramid height) |tan theta| + the tangential spacing (its
# perpendicular end wholly outside the box) and intersected with the coupon box
# before the fragment; the mesh ends with a block of `layers` = m sheared layers of
# axial spacing `spacing` = lc_end whose last station is the face plane (every end
# node on the face), so every lateral quadrangle stays a planar trapezoid and every
# pyramid apex of the block stays inside the box. Two regimes (mesher design round 2
# F2b, section 3.2; supervisor decisions 358 / 363 / 437), split by the section's
# end-spacing cap lc_cap (face_end_spacing_cap):
#   regime I  (4 h_pyr |tan theta| <= lc_cap): lc_end = max(lc_tangent, 4 h_pyr |tan theta|),
#             m = ceil(2 r_env |tan theta| / lc_end) - the A2 (4) formulas bitwise; every
#             layer's axial thickness lies in [lc_end / 2, 3 lc_end / 2];
#   regime II (4 h_pyr |tan theta| > lc_cap): lc_end = lc_cap and m = max(ceil(r_env |tan
#             theta| / (lc_cap - 2 h_pyr |tan theta|)), the regime-I count at lc_cap), so the
#             thinnest layer t_min = lc_cap - r_env |tan theta| / m keeps the exact apex rule
#             t_min >= 2 h_pyr |tan theta| while the shear per layer stays <= lc_cap / 2;
#   beyond the validity ceiling 2 h_pyr |tan theta| >= lc_cap no block exists: the build
#             fails closed (ScopeGuard[SteepFaceCrossing], checked by the caller).
# A tube end with theta == 0 exactly (every rectilinear coupon) has no FaceEnd and
# takes the unchanged path.
struct FaceEnd
    end_index::Int            # 0: the tube start (s_start) lies on the face, 1: the end
    face::String              # "x0" / "x1" / "y0" / "y1"
    normal::Vector{Float64}   # outward face normal
    theta::Float64            # angle between the tube axis and the face normal (rad), > 0
    kappa_u::Float64
    kappa_w::Float64
    layers::Int               # m, the sheared end-block layers
    spacing::Float64          # lc_end, the axial spacing of the end block
    envelope_shear::Float64   # r_env |tan theta|: the axial shear across the tube envelope
    over_length::Float64      # the CAD over-length beyond the face before the box intersection
    face_axis::Int            # the plan-view axis (1 x, 2 y) the face is normal to (0: unknown)
    face_value::Float64       # the face coordinate on that axis (an ArcTube crosses it per node)
    regime::Int               # 1 or 2 (design round 2 F2b)
    spacing_cap::Float64      # lc_cap of the tube's section
    apex_thickness::Float64   # 2 h_pyr |tan theta|: the thinnest layer the apex rule admits
end

function FaceEnd(
    end_index,
    face,
    normal,
    theta,
    kappa_u,
    kappa_w,
    envelope_radius,
    pyramid_height,
    lc_tangent;
    spacing_cap,
    face_axis=0,
    face_value=NaN
)
    theta > 0.0 || error("a face end needs a positive tilt")
    spacing_cap > 0.0 || error("a face end needs a positive end-spacing cap")
    slope = abs(tan(theta))
    apex_thickness = 2.0 * pyramid_height * slope
    regime_one_spacing = max(lc_tangent, 4.0 * pyramid_height * slope)
    regime_one_layers(spacing) =
        max(1, ceil(Int, 2.0 * envelope_radius * slope / spacing * (1.0 - 1.0e-9)))
    if 4.0 * pyramid_height * slope <= spacing_cap
        regime = 1
        spacing = regime_one_spacing
        layers = regime_one_layers(spacing)
    else
        apex_thickness < spacing_cap || error(
            "face end at $(rad2deg(theta)) degrees lies beyond the validity ceiling of the " *
            "capped end block: 2 h_pyr |tan theta| = $apex_thickness >= lc_cap $spacing_cap"
        )
        regime = 2
        spacing = spacing_cap
        layers = max(
            ceil(Int, envelope_radius * slope / (spacing_cap - apex_thickness)),
            regime_one_layers(spacing_cap)
        )
    end
    # The exact apex rule holds in both regimes (regime I: t_min >= lc_end / 2 >= 2 h_pyr
    # |tan theta| up to the layer count's 1e-9 rounding slack; regime II by the layer
    # count), fail closed otherwise.
    spacing - envelope_radius * slope / layers >= apex_thickness * (1.0 - 1.0e-9) ||
        error("face end block violates the apex rule")
    return FaceEnd(
        end_index,
        face,
        collect(Float64, normal),
        theta,
        kappa_u,
        kappa_w,
        layers,
        spacing,
        envelope_radius * slope,
        envelope_radius * slope + lc_tangent,
        face_axis,
        face_value,
        regime,
        spacing_cap,
        apex_thickness
    )
end

# The census record of a face end: the block's axial layer thickness over the tube
# envelope lies in EndSpacing -+ EnvelopeShear / Layers (within [lc_end / 2, 3 lc_end / 2]
# and at or above ApexThickness, the exact apex rule); Regime "I" / "II" and the section's
# EndSpacingCap bind lc_end <= lc_cap (design round 2 F2b).
face_end_record(face_end::FaceEnd) = Dict{String, Any}(
    "End" => face_end.end_index == 0 ? "start" : "end",
    "Face" => face_end.face,
    "ThetaDegrees" => rad2deg(face_end.theta),
    "Layers" => face_end.layers,
    "EndSpacing" => face_end.spacing,
    "EnvelopeShear" => face_end.envelope_shear,
    "LayerThicknessRange" => [
        face_end.spacing - face_end.envelope_shear / face_end.layers,
        face_end.spacing + face_end.envelope_shear / face_end.layers
    ],
    "OverLength" => face_end.over_length,
    "Kappa" => [face_end.kappa_u, face_end.kappa_w],
    "Regime" => face_end.regime == 1 ? "I" : "II",
    "EndSpacingCap" => face_end.spacing_cap,
    "ApexThickness" => face_end.apex_thickness
)

# A SMOOTH joint of a straight tube with an arc tube (block (b) design A3 (2), decision
# 303): the shared cross-section is OWNED BY THE ARC - its radial plane through the joint
# vertex, with the arc's frame (origin = the joint vertex on the tube's edge, n = the
# arc's outward normal there, b). The straight tube's end nodes ARE the arc tube's nodes
# (TubeMesh adopts their tags), so its last layer is sheared onto the radial plane by the
# post-snap tilt `tilt` (<= the snap over rho: 2e-6 rad on the loop end). The CAD solid of
# the straight tube is extruded over-long by `over_length` and cut by the half-space
# behind the radial plane (`normal` points out of the straight tube) when tilt > 0, so the
# two solids share one identical planar face (no coincident-face tolerance game); at tilt
# == 0 exactly the perpendicular end already is that plane. No cap, no ball, clearance 0.
struct JointEnd
    end_index::Int
    origin::Vector{Float64}
    n::Vector{Float64}
    b::Vector{Float64}
    normal::Vector{Float64}
    tilt::Float64
    over_length::Float64
end

function JointEnd(end_index, origin, n, b, normal, tilt, envelope_radius, lc_tangent)
    tilt >= 0.0 || error("a joint tilt is an angle")
    tilt < 0.25 * pi || error("a joint tilt of $(rad2deg(tilt)) degrees is not smooth")
    return JointEnd(
        end_index,
        collect(Float64, origin),
        collect(Float64, n),
        collect(Float64, b),
        collect(Float64, normal),
        tilt,
        tilt > 0.0 ? envelope_radius * abs(tan(tilt)) + lc_tangent : 0.0
    )
end

joint_end_record(joint::JointEnd) = Dict{String, Any}(
    "End" => joint.end_index == 0 ? "start" : "end",
    "Origin" => copy(joint.origin),
    "Normal" => copy(joint.n),
    "TiltRadians" => joint.tilt,
    "OverLength" => joint.over_length,
    "PlaneCut" => joint.tilt > 0.0
)

# The prism tubes: a straight EdgeTube or a revolved ArcTube (design 1.2 (3)); both carry
# origin-free section geometry through tube_point(tube, u, w, s) (u along the outward
# normal, w along b, s the axis coordinate - the arc length on an ArcTube), the layer
# stations, the face ends and the smooth joints.
abstract type AbstractTube end

# A tube: the edge line origin (a point on the edge), the frame (n, b, e), the
# extrusion interval [s_start, s_end] along e from the origin, the layer
# boundaries (stations) s_start = stations[1] < ... < stations[layers + 1] = s_end
# on the axis, the face ends (at most one per end) and, per station, the shear of
# the face-ended blocks: the station of the cross-section node (u, w) at index i is
# stations[i + 1] - u shear_u[i + 1] - w shear_w[i + 1] (every shear 0 on a tube
# without face ends: the shear vectors are then empty and the station is the axis one).
# `joints`: the smooth joints by end index (design A3 (2)); empty on every tube built so far.
struct EdgeTube <: AbstractTube
    origin::Vector{Float64}
    n::Vector{Float64}
    b::Vector{Float64}
    e::Vector{Float64}
    s_start::Float64
    s_end::Float64
    layers::Int
    stations::Vector{Float64}
    face_ends::Vector{FaceEnd}
    shear_u::Vector{Float64}
    shear_w::Vector{Float64}
    joints::Vector{JointEnd}
end

# Uniform layers: the smallest number of equal layers whose spacing does not
# exceed `spacing` (a tube whose length is a multiple of the spacing keeps it
# exactly; tube_spacing records the largest layer actually used). The face ends
# are recorded; their sheared blocks are installed by face_ended_tube_stations.
function EdgeTube(
    origin,
    n,
    b,
    s_start,
    s_end,
    spacing;
    face_ends=FaceEnd[],
    joints=JointEnd[]
)
    n = collect(Float64, n) ./ norm(n)
    b = collect(Float64, b) ./ norm(b)
    abs(dot(n, b)) < 1.0e-12 || error("tube frame must be orthogonal")
    e = cross(n, b)
    s_end > s_start || error("tube interval must be increasing")
    spacing > 0.0 || error("tube spacing must be positive")
    extent = s_end - s_start
    layers = max(1, ceil(Int, extent / spacing * (1.0 - 1.0e-9)))
    stations = [s_start + extent * i / layers for i = 0:layers]
    face_ends = collect(FaceEnd, face_ends)
    joints = collect(JointEnd, joints)
    length(unique(face_end.end_index for face_end in face_ends)) == length(face_ends) ||
        error("a tube end has two face ends")
    length(unique(joint.end_index for joint in joints)) == length(joints) ||
        error("a tube end has two joints")
    isempty(intersect([f.end_index for f in face_ends], [j.end_index for j in joints])) || error("a tube end is both a face end and a joint")
    return EdgeTube(
        collect(Float64, origin),
        n,
        b,
        e,
        s_start,
        s_end,
        layers,
        stations,
        face_ends,
        Float64[],
        Float64[],
        joints
    )
end

# The same tube with the given layer boundaries (and the end-block shears of a
# face-ended tube, one pair per station).
function EdgeTube(
    tube::EdgeTube,
    stations::AbstractVector;
    shear_u=Float64[],
    shear_w=Float64[]
)
    stations = collect(Float64, stations)
    length(stations) >= 2 && stations[1] == tube.s_start && stations[end] == tube.s_end ||
        error("tube stations must run from s_start to s_end")
    all(diff(stations) .> 0.0) || error("tube stations must increase")
    shear_u = collect(Float64, shear_u)
    shear_w = collect(Float64, shear_w)
    isempty(shear_u) == isempty(shear_w) || error("tube shears come in pairs")
    isempty(shear_u) ||
        length(shear_u) == length(stations) ||
        error("tube shears need one value per station")
    isempty(shear_u) == isempty(tube.face_ends) ||
        error("a face-ended tube needs its end-block shears and a plain tube none")
    return EdgeTube(
        tube.origin,
        tube.n,
        tube.b,
        tube.e,
        tube.s_start,
        tube.s_end,
        length(stations) - 1,
        stations,
        tube.face_ends,
        shear_u,
        shear_w,
        tube.joints
    )
end

tube_layer_thicknesses(tube::AbstractTube) = diff(tube.stations)
# The extrusion spacing actually used: the largest layer (equal to every layer of
# a uniform tube).
tube_spacing(tube::AbstractTube) = maximum(tube_layer_thicknesses(tube))

function tube_point(tube::EdgeTube, u, w, s)
    return tube.origin .+ u .* tube.n .+ w .* tube.b .+ s .* tube.e
end

# Axis coordinate of station i (an integer 0:layers) or of a point between two
# stations (i + 1/2: the pyramid apex station of layer i).
function tube_station(tube::AbstractTube, i)
    k = clamp(floor(Int, i), 0, tube.layers - 1)
    return tube.stations[k + 1] + (i - k) * (tube.stations[k + 2] - tube.stations[k + 1])
end

# ---------------------------------------------------------------------------
# ArcTube (block (b) design 1.2 (3), decision 303): the tube of an ARC metal side, revolved
# about the vertical axis through the arc centre. The axis coordinate s is the ARC LENGTH
# on the edge circle of radius rho: theta(s) = theta0 + orientation s / rho, and the
# cross-section point (u, w) lies at tube_point = centre + (rho + sigma u) radial(theta) +
# w b, u along the outward normal n = sigma radial (sigma = +1 when the dielectric lies
# outside the circle: convex metal; -1 when the metal lies outside), w along b. The travel
# direction e = n x b fixes orientation = -sigma b_z, as the straight tube's e does. The
# OCC volumes are the revolve of each material-group sector face drawn in the radial plane
# at theta(cad_start); their sidewall rays are the cylinder of radius rho about the centre
# = the exact wall of the metal loft and of the trench, so the fragment splits the sectors
# at the wall by construction (the root cause of "2 volume descendants" under chords).
# Face ends (A2 for arcs): the node (u, w) of the end block meets the box face where ITS
# circle (radius rho + sigma u) crosses the face plane (closed form, per node), the block
# fraction grows linearly to 1 at the face; `block_fraction` / `block_end` per station
# (0 / -1 on an interior station). Joints (A3 (2)): the arc owns the shared section.
struct ArcTube <: AbstractTube
    centre::Vector{Float64}    # (cx, cy, z): the axis point at the tube's edge height
    rho::Float64
    sigma::Float64
    b::Vector{Float64}
    theta0::Float64
    orientation::Float64
    s_start::Float64
    s_end::Float64
    layers::Int
    stations::Vector{Float64}
    face_ends::Vector{FaceEnd}
    block_fraction::Vector{Float64}
    block_end::Vector{Int}
    joints::Vector{JointEnd}
end

function ArcTube(
    centre,
    rho,
    sigma,
    b,
    theta0,
    s_start,
    s_end,
    spacing;
    face_ends=FaceEnd[],
    joints=JointEnd[]
)
    rho > 0.0 || error("arc tube radius must be positive")
    abs(abs(sigma) - 1.0) <= 0.0 || error("arc tube sign must be +1 or -1")
    b = collect(Float64, b) ./ norm(b)
    abs(abs(b[3]) - 1.0) <= 1.0e-12 && abs(b[1]) <= 1.0e-12 && abs(b[2]) <= 1.0e-12 ||
        error("an arc tube revolves about the process normal")
    s_end > s_start || error("tube interval must be increasing")
    spacing > 0.0 || error("tube spacing must be positive")
    extent = s_end - s_start
    layers = max(1, ceil(Int, extent / spacing * (1.0 - 1.0e-9)))
    stations = [s_start + extent * i / layers for i = 0:layers]
    face_ends = collect(FaceEnd, face_ends)
    joints = collect(JointEnd, joints)
    length(unique(face_end.end_index for face_end in face_ends)) == length(face_ends) ||
        error("a tube end has two face ends")
    length(unique(joint.end_index for joint in joints)) == length(joints) ||
        error("a tube end has two joints")
    isempty(intersect([f.end_index for f in face_ends], [j.end_index for j in joints])) || error("a tube end is both a face end and a joint")
    all(f.face_axis in (1, 2) && isfinite(f.face_value) for f in face_ends) ||
        error("an arc tube face end needs its face plane (axis, value)")
    return ArcTube(
        collect(Float64, centre),
        Float64(rho),
        Float64(sigma),
        b,
        Float64(theta0),
        -sigma * b[3],
        Float64(s_start),
        Float64(s_end),
        layers,
        stations,
        face_ends,
        Float64[],
        Int[],
        joints
    )
end

# The same arc tube with the given layer boundaries (and the end-block fractions of a
# face-ended tube, one pair per station).
function ArcTube(
    tube::ArcTube,
    stations::AbstractVector;
    block_fraction=Float64[],
    block_end=Int[]
)
    stations = collect(Float64, stations)
    length(stations) >= 2 && stations[1] == tube.s_start && stations[end] == tube.s_end ||
        error("tube stations must run from s_start to s_end")
    all(diff(stations) .> 0.0) || error("tube stations must increase")
    block_fraction = collect(Float64, block_fraction)
    block_end = collect(Int, block_end)
    length(block_fraction) == length(block_end) ||
        error("tube block fractions come with their ends")
    isempty(block_fraction) ||
        length(block_fraction) == length(stations) ||
        error("tube block fractions need one value per station")
    isempty(block_fraction) == isempty(tube.face_ends) ||
        error("a face-ended tube needs its end-block fractions and a plain tube none")
    return ArcTube(
        tube.centre,
        tube.rho,
        tube.sigma,
        tube.b,
        tube.theta0,
        tube.orientation,
        tube.s_start,
        tube.s_end,
        length(stations) - 1,
        stations,
        tube.face_ends,
        block_fraction,
        block_end,
        tube.joints
    )
end

arc_angle(tube::ArcTube, s) = tube.theta0 + tube.orientation * s / tube.rho
arc_radial(theta) = [cos(theta), sin(theta), 0.0]

function tube_point(tube::ArcTube, u, w, s)
    theta = arc_angle(tube, s)
    return tube.centre .+ (tube.rho + tube.sigma * u) .* arc_radial(theta) .+ w .* tube.b
end

# The local frame of the arc tube at axis coordinate s: (origin on the edge, n, b, e).
function arc_frame(tube::ArcTube, s)
    theta = arc_angle(tube, s)
    radial = arc_radial(theta)
    n = tube.sigma .* radial
    e = cross(n, tube.b)
    return (origin=tube.centre .+ tube.rho .* radial, n=n, b=copy(tube.b), e=e)
end

# The arc length on the AXIS circle at which the node circle of radius rho + sigma u crosses
# the face plane of `face_end` nearest to the tube's end angle there (closed form).
function arc_face_station(tube::ArcTube, face_end::FaceEnd, u, w)
    radius = tube.rho + tube.sigma * u
    radius > 0.0 || error("arc tube node radius is not positive")
    s_axis = face_end.end_index == 0 ? tube.s_start : tube.s_end
    theta_axis = arc_angle(tube, s_axis)
    offset = (face_end.face_value - tube.centre[face_end.face_axis]) / radius
    # Every node circle must reach the face plane: |x_face - c| <= rho - (Radius + PyramidHeight),
    # i.e. rho (1 - sin theta) > the tube envelope. The caller fails closed BY NAME before any CAD
    # (ScopeGuard[ArcFaceEnds], mesh_spatial_coupon.jl build_edge_tubes!; mesher design round 3
    # part M 3.3 Fact 1); this is the invariant's assertion.
    abs(offset) <= 1.0 || error(
        "the node circle of radius $radius of an arc tube does not reach its box face " *
        "(|x_face - c| = $(abs(face_end.face_value - tube.centre[face_end.face_axis])) > " *
        "the node radius; the caller guards this at ScopeGuard[ArcFaceEnds])"
    )
    candidates =
        face_end.face_axis == 1 ? (acos(offset), -acos(offset)) :
        (asin(offset), pi - asin(offset))
    theta = argmin(c -> abs(rem2pi(c - theta_axis, RoundNearest)), candidates)
    theta = theta_axis + rem2pi(theta - theta_axis, RoundNearest)
    return tube.orientation * (theta - tube.theta0) * tube.rho
end

# Station of the cross-section point (u, w) at index i on an arc tube: the axis station
# plus the end-block fraction times the per-node offset of the face crossing from the axis
# end (the axis station itself on a tube without face ends).
function tube_station(tube::ArcTube, i, u, w)
    s = tube_station(tube, i)
    isempty(tube.block_fraction) && return s
    k = clamp(floor(Int, i), 0, tube.layers - 1)
    f = i - k
    fraction =
        tube.block_fraction[k + 1] +
        f * (tube.block_fraction[k + 2] - tube.block_fraction[k + 1])
    fraction > 0.0 || return s
    end_index = tube.block_end[k + 1] >= 0 ? tube.block_end[k + 1] : tube.block_end[k + 2]
    face_end = face_end_at(tube, end_index)
    face_end === nothing && return s
    s_axis = end_index == 0 ? tube.s_start : tube.s_end
    return s + fraction * (arc_face_station(tube, face_end, u, w) - s_axis)
end

# Station of the cross-section point (u, w) at index i: the axis station minus the
# (linearly interpolated) end-block shear; the axis station itself on a tube
# without face ends.
function tube_station(tube::EdgeTube, i, u, w)
    s = tube_station(tube, i)
    isempty(tube.shear_u) && return s
    k = clamp(floor(Int, i), 0, tube.layers - 1)
    f = i - k
    su = tube.shear_u[k + 1] + f * (tube.shear_u[k + 2] - tube.shear_u[k + 1])
    sw = tube.shear_w[k + 1] + f * (tube.shear_w[k + 2] - tube.shear_w[k + 1])
    return s - u * su - w * sw
end

face_end_at(tube::AbstractTube, end_index) = (
    i=findfirst(face_end -> face_end.end_index == end_index, tube.face_ends);
    i === nothing ? nothing : tube.face_ends[i]
)

joint_at(tube::AbstractTube, end_index) = (
    i=findfirst(joint -> joint.end_index == end_index, tube.joints);
    i === nothing ? nothing : tube.joints[i]
)

# Whether the tube's end sections are the plain perpendicular ones (no face end, no joint):
# the legacy entity formulas apply bitwise.
plain_ends(tube::AbstractTube) = isempty(tube.face_ends) && isempty(tube.joints)

# Axis coordinate at which the cross-section point (u, w) meets the tube's end
# (end_index 0 the start, 1 the end): the face plane at a face end, the
# perpendicular end otherwise.
function tube_end_station(tube::EdgeTube, end_index, u, w)
    s = end_index == 0 ? tube.s_start : tube.s_end
    face_end = face_end_at(tube, end_index)
    face_end === nothing && return s
    return s - u * face_end.kappa_u - w * face_end.kappa_w
end

function tube_end_station(tube::ArcTube, end_index, u, w)
    s = end_index == 0 ? tube.s_start : tube.s_end
    face_end = face_end_at(tube, end_index)
    face_end === nothing && return s
    return arc_face_station(tube, face_end, u, w)
end

# The position of the cross-section point (u, w) at the tube's end: on a smooth joint the
# arc's radial section (its frame: the joint owns the nodes, design A3 (2)), else the point
# at the end station (a face end's face plane, a perpendicular end).
function tube_end_point(tube::AbstractTube, end_index, u, w)
    joint = joint_at(tube, end_index)
    joint === nothing || return joint.origin .+ u .* joint.n .+ w .* joint.b
    return tube_point(tube, u, w, tube_end_station(tube, end_index, u, w))
end

# The cross-section point (u, w) at station index i: the end section's point at the two
# ends (tube_end_point: the joint frame or the face plane), the station point inside.
function tube_section_point(tube::AbstractTube, i, u, w)
    i == 0 && !plain_ends(tube) && return tube_end_point(tube, 0, u, w)
    i == tube.layers && !plain_ends(tube) && return tube_end_point(tube, 1, u, w)
    return tube_point(tube, u, w, tube_station(tube, i, u, w))
end

# The CAD extrusion interval of the tube before the box / plane intersections: over-long
# at every face end and at every tilted joint.
function tube_cad_interval(tube::AbstractTube)
    start_face = face_end_at(tube, 0)
    end_face = face_end_at(tube, 1)
    start_joint = joint_at(tube, 0)
    end_joint = joint_at(tube, 1)
    return tube.s_start - (start_face === nothing ? 0.0 : start_face.over_length) -
           (start_joint === nothing ? 0.0 : start_joint.over_length),
    tube.s_end +
    (end_face === nothing ? 0.0 : end_face.over_length) +
    (end_joint === nothing ? 0.0 : end_joint.over_length)
end

# Samples per prescribed size of the arclength quadrature that places the layer
# boundaries (a resolution constant, not a mesh target).
const TUBE_LAYER_SAMPLES_PER_SIZE = 8

# The minimum of the sampled, piecewise-linear size over the axis interval [a, b].
function sampled_size_minimum(positions, sizes, a, b)
    function size_at(s)
        i = clamp(searchsortedlast(positions, s), 1, length(positions) - 1)
        fraction = (s - positions[i]) / (positions[i + 1] - positions[i])
        return sizes[i] + fraction * (sizes[i + 1] - sizes[i])
    end
    lowest = min(size_at(a), size_at(b))
    for i = searchsortedfirst(positions, a):searchsortedlast(positions, b)
        lowest = min(lowest, sizes[i])
    end
    return lowest
end

# Layer boundaries following a size field along the tube axis (supervisor decision
# 40): `size_at(s)` is the composed size prescribed at axis coordinate s, capped at
# `spacing`; the field is sampled adaptively (TUBE_LAYER_SAMPLES_PER_SIZE samples
# per local size), gradient-limited along the axis to the slope
# (growth - 1) / growth, and the layers equidistribute the arclength integral of
# 1 / size with ceil(integral) layers, so every layer is at most the largest size
# it spans (strictly below spacing: the cap is spacing x (1 - 1e-9), a tube whose
# length is an exact multiple of a uniform spacing takes one more layer than the
# uniform constructor). At an end that lies on a surface (`surface_start` /
# `surface_end`: the tube ends on the outer box) the end layer is instead the
# size evaluated at the surface, at most the limited field over the layer span
# (iterated t <- min over [s, s + t] from the surface value until stationary), so
# that a field growing away from the surface is resolved from the surface, not
# from the layer midpoint; the remaining interval is equidistributed as above
# (decision 41). The neighbour ratio of consecutive layers is checked against
# `growth` (fail closed; the recorded MaximumNeighbourRatio): with the limited
# slope k a layer against a field growing at k spans at most 2 (e^k - 1) / k
# times the size at its start, so an equidistributed layer is at most that
# factor (1.30 at growth 2) above its predecessor of the same law, and at most
# (1 + k) x that factor (1.95 at growth 2) above a surface layer.
# Returns the stations and the sampled (s, limited size) pairs for the census.
function graded_tube_stations(
    s_start,
    s_end,
    size_at,
    spacing,
    growth;
    surface_start::Bool=false,
    surface_end::Bool=false
)
    s_end > s_start || error("tube interval must be increasing")
    spacing > 0.0 || error("tube spacing must be positive")
    growth > 1.0 || error("tube layer growth must exceed 1")
    cap = spacing * (1.0 - 1.0e-9)
    positions = [Float64(s_start)]
    sizes = Float64[]
    function prescribed(s)
        h = min(cap, size_at(s))
        isfinite(h) && h > 0.0 || error("tube axis size at $s is not positive: $h")
        return h
    end
    push!(sizes, prescribed(s_start))
    while positions[end] < s_end
        step = sizes[end] / TUBE_LAYER_SAMPLES_PER_SIZE
        next = min(s_end, positions[end] + step)
        push!(positions, next)
        push!(sizes, prescribed(next))
    end
    slope = (growth - 1.0) / growth
    function limit!()
        for i = 2:length(sizes)
            sizes[i] =
                min(sizes[i], sizes[i - 1] + slope * (positions[i] - positions[i - 1]))
        end
        for i = (length(sizes) - 1):-1:1
            sizes[i] =
                min(sizes[i], sizes[i + 1] + slope * (positions[i + 1] - positions[i]))
        end
    end
    limit!()
    # The limiter lowers sizes below the sampling density they were sampled at;
    # refine the coarse intervals (midpoints at the prescribed size, then the
    # limiter again) until every interval is resolved.
    while true
        coarse = [
            i for i = 1:(length(positions) - 1) if positions[i + 1] - positions[i] >
            min(sizes[i], sizes[i + 1]) / TUBE_LAYER_SAMPLES_PER_SIZE * (1.0 + 1.0e-9)
        ]
        isempty(coarse) && break
        refined_positions = Float64[]
        refined_sizes = Float64[]
        coarse_set = Set(coarse)
        for i = 1:(length(positions) - 1)
            push!(refined_positions, positions[i]);
            push!(refined_sizes, sizes[i])
            i in coarse_set || continue
            mid = 0.5 * (positions[i] + positions[i + 1])
            push!(refined_positions, mid);
            push!(refined_sizes, prescribed(mid))
        end
        push!(refined_positions, positions[end]);
        push!(refined_sizes, sizes[end])
        positions = refined_positions
        sizes = refined_sizes
        limit!()
    end
    # A surface layer: t <- min of the limited field over [s, s + t] from the
    # surface value until stationary (non-increasing, bounded below by the field).
    function surface_layer(s, direction)
        thickness = sampled_size_minimum(positions, sizes, s, s)
        for _ = 1:1000
            span = sort([s, s + direction * thickness])
            lowest = sampled_size_minimum(positions, sizes, span[1], span[2])
            lowest >= thickness * (1.0 - 1.0e-9) && return thickness
            thickness = lowest
        end
        return error("tube surface layer at $s did not converge")
    end
    start_layer = surface_start ? surface_layer(Float64(s_start), 1.0) : 0.0
    end_layer = surface_end ? surface_layer(Float64(s_end), -1.0) : 0.0
    interior_start = s_start + start_layer
    interior_end = s_end - end_layer
    interior_end > interior_start ||
        error("tube of length $(s_end - s_start) is shorter than its surface layers")
    cumulative = zeros(length(positions))
    for i = 2:length(positions)
        cumulative[i] =
            cumulative[i - 1] +
            (positions[i] - positions[i - 1]) * 0.5 * (1.0 / sizes[i - 1] + 1.0 / sizes[i])
    end
    function cumulative_at(s)
        i = clamp(searchsortedlast(positions, s), 1, length(positions) - 1)
        fraction = (s - positions[i]) / (positions[i + 1] - positions[i])
        return cumulative[i] + fraction * (cumulative[i + 1] - cumulative[i])
    end
    interior_integral = cumulative_at(interior_end) - cumulative_at(interior_start)
    layers = max(1, ceil(Int, interior_integral))
    stations = [Float64(s_start)]
    surface_start && push!(stations, interior_start)
    for k = 1:(layers - 1)
        target = cumulative_at(interior_start) + k * interior_integral / layers
        i = clamp(searchsortedlast(cumulative, target), 1, length(positions) - 1)
        fraction = (target - cumulative[i]) / (cumulative[i + 1] - cumulative[i])
        push!(stations, positions[i] + fraction * (positions[i + 1] - positions[i]))
    end
    surface_end && push!(stations, interior_end)
    push!(stations, Float64(s_end))
    thicknesses = diff(stations)
    all(thicknesses .> 0.0) || error("tube layer placement produced an empty layer")
    ratios = thicknesses[2:end] ./ thicknesses[1:(end - 1)]
    neighbour_ratio = isempty(ratios) ? 1.0 : max(maximum(ratios), 1.0 / minimum(ratios))
    neighbour_ratio <= growth * (1.0 + 1.0e-9) ||
        error("tube layers grow by $neighbour_ratio between neighbours, above $growth")
    return stations, positions, sizes
end

# Stations and end-block shears of a tube with face ends (design A2): the interior
# interval - the tube minus its end blocks of m x lc_end - takes graded_tube_stations
# (a face end replaces the surface layer of its end); each block then adds m layers
# of axial spacing lc_end whose shear grows linearly from 0 at the block's inner
# station to kappa at the face, so station m of the end block (station 0 of a start
# block) IS the face plane. A tube shorter than its end blocks fails closed.
function face_ended_tube_stations(
    tube::EdgeTube,
    size_at,
    spacing,
    growth;
    surface_start::Bool=false,
    surface_end::Bool=false
)
    start_face = face_end_at(tube, 0)
    end_face = face_end_at(tube, 1)
    start_block = start_face === nothing ? 0.0 : start_face.layers * start_face.spacing
    end_block = end_face === nothing ? 0.0 : end_face.layers * end_face.spacing
    a = tube.s_start + start_block
    b = tube.s_end - end_block
    b > a || error(
        "tube of length $(tube.s_end - tube.s_start) is shorter than its face-end blocks " *
        "($start_block at the start, $end_block at the end)"
    )
    interior, positions, sizes = graded_tube_stations(
        a,
        b,
        size_at,
        spacing,
        growth;
        surface_start=start_face === nothing && surface_start,
        surface_end=end_face === nothing && surface_end
    )
    stations = Float64[]
    shear_u = Float64[]
    shear_w = Float64[]
    if start_face !== nothing
        m = start_face.layers
        for k = 0:(m - 1)
            push!(stations, k == 0 ? tube.s_start : tube.s_start + k * start_face.spacing)
            push!(shear_u, (m - k) / m * start_face.kappa_u)
            push!(shear_w, (m - k) / m * start_face.kappa_w)
        end
    end
    append!(stations, interior)
    append!(shear_u, zeros(length(interior)))
    append!(shear_w, zeros(length(interior)))
    if end_face !== nothing
        m = end_face.layers
        for k = 1:m
            push!(stations, k == m ? tube.s_end : b + k * end_face.spacing)
            push!(shear_u, k / m * end_face.kappa_u)
            push!(shear_w, k / m * end_face.kappa_w)
        end
    end
    return stations, shear_u, shear_w, positions, sizes
end

# The arc tube's version (design A2 for arcs): the same interior / block stations, the
# block recorded as a fraction per station (1 at the face, 0 at the block's inner station)
# with the end it belongs to; the per-node face crossing is tube_station's.
function face_ended_tube_stations(
    tube::ArcTube,
    size_at,
    spacing,
    growth;
    surface_start::Bool=false,
    surface_end::Bool=false
)
    start_face = face_end_at(tube, 0)
    end_face = face_end_at(tube, 1)
    start_block = start_face === nothing ? 0.0 : start_face.layers * start_face.spacing
    end_block = end_face === nothing ? 0.0 : end_face.layers * end_face.spacing
    a = tube.s_start + start_block
    b = tube.s_end - end_block
    b > a || error(
        "tube of length $(tube.s_end - tube.s_start) is shorter than its face-end blocks " *
        "($start_block at the start, $end_block at the end)"
    )
    interior, positions, sizes = graded_tube_stations(
        a,
        b,
        size_at,
        spacing,
        growth;
        surface_start=start_face === nothing && surface_start,
        surface_end=end_face === nothing && surface_end
    )
    stations = Float64[]
    fraction = Float64[]
    block_end = Int[]
    if start_face !== nothing
        m = start_face.layers
        for k = 0:(m - 1)
            push!(stations, k == 0 ? tube.s_start : tube.s_start + k * start_face.spacing)
            push!(fraction, (m - k) / m)
            push!(block_end, 0)
        end
    end
    append!(stations, interior)
    append!(fraction, zeros(length(interior)))
    append!(block_end, fill(-1, length(interior)))
    if end_face !== nothing
        m = end_face.layers
        for k = 1:m
            push!(stations, k == m ? tube.s_end : b + k * end_face.spacing)
            push!(fraction, k / m)
            push!(block_end, 1)
        end
    end
    return stations, fraction, block_end, positions, sizes
end

# Per-tube layer statistics against the prescribed (gradient-limited) sizes
# sampled by graded_tube_stations: thickness minimum / P50 / maximum, the
# thickness and the prescribed size at both ends, the achieved-over-prescribed
# ratio (each layer over the size at its midpoint) and the neighbour ratio.
function tube_layer_statistics(tube::AbstractTube, positions, sizes)
    thicknesses = tube_layer_thicknesses(tube)
    # The sheared end blocks of a face-ended tube (design A2) have their own spacing
    # lc_end >= the tangential spacing: the neighbour-ratio law of decision 40 is
    # judged over the interior layers, and each block is recorded with its ratio to
    # the adjacent interior layer (FaceEndBlocks; absent on a plain tube).
    start_face = face_end_at(tube, 0)
    end_face = face_end_at(tube, 1)
    start_block = start_face === nothing ? 0 : start_face.layers
    end_block = end_face === nothing ? 0 : end_face.layers
    interior = thicknesses[(start_block + 1):(length(thicknesses) - end_block)]
    function size_at(s)
        i = clamp(searchsortedlast(positions, s), 1, length(positions) - 1)
        fraction = (s - positions[i]) / (positions[i + 1] - positions[i])
        return sizes[i] + fraction * (sizes[i + 1] - sizes[i])
    end
    midpoints = 0.5 .* (tube.stations[1:(end - 1)] .+ tube.stations[2:end])
    # The end blocks of a face-ended tube lie outside the sampled interior interval:
    # their layers are judged against the size at the nearest sampled position.
    achieved = [
        thicknesses[i] / size_at(clamp(midpoints[i], positions[1], positions[end])) for
        i in eachindex(thicknesses)
    ]
    ratios = interior[2:end] ./ interior[1:(end - 1)]
    median(values) = sort(values)[cld(length(values), 2)]
    statistics = Dict{String, Any}(
        "Minimum" => minimum(thicknesses),
        "P50" => median(thicknesses),
        "Maximum" => maximum(thicknesses),
        "AtStart" => thicknesses[1],
        "AtEnd" => thicknesses[end],
        "PrescribedAtStart" => sizes[1],
        "PrescribedAtEnd" => sizes[end],
        "AchievedOverPrescribed" => Dict{String, Any}(
            "Minimum" => minimum(achieved),
            "P50" => median(achieved),
            "Maximum" => maximum(achieved)
        ),
        "MaximumNeighbourRatio" =>
            isempty(ratios) ? 1.0 : max(maximum(ratios), 1.0 / minimum(ratios))
    )
    blocks = Dict{String, Any}()
    if start_block > 0
        blocks["Start"] = Dict{String, Any}(
            "Layers" => start_block,
            "Thicknesses" => thicknesses[1:start_block],
            "NeighbourRatio" => thicknesses[start_block] / thicknesses[start_block + 1]
        )
    end
    if end_block > 0
        blocks["End"] = Dict{String, Any}(
            "Layers" => end_block,
            "Thicknesses" => thicknesses[(end - end_block + 1):end],
            "NeighbourRatio" =>
                thicknesses[end - end_block + 1] / thicknesses[end - end_block]
        )
    end
    isempty(blocks) || (statistics["FaceEndBlocks"] = blocks)
    return statistics
end

# OCC tube volumes (one polygon-sector prism per material group), returned as
# (dim, tag) pairs with their material group, before synchronization. A tube with
# face ends is extruded over its CAD interval (over-long beyond every face end) and
# each volume is intersected with the coupon box `box` = (lower, upper), so that the
# tube ends exactly on the face (design A2; the single-descendant check of the
# fragment is unchanged).
# The half-space {p : (p - origin) . normal <= 0} as an OCC solid large enough to contain
# the whole coupon (`reach` = a length beyond every coupon dimension): an axis-aligned box
# rotated so that its +z face normal becomes `normal`, then moved to `origin`.
function add_half_space!(occ, origin, normal, reach)
    normal = collect(Float64, normal) ./ norm(normal)
    solid = occ.addBox(-reach, -reach, -2.0 * reach, 2.0 * reach, 2.0 * reach, 2.0 * reach)
    axis = cross([0.0, 0.0, 1.0], normal)
    angle = acos(clamp(normal[3], -1.0, 1.0))
    if norm(axis) > 1.0e-12
        occ.rotate([(3, solid)], 0.0, 0.0, 0.0, axis[1], axis[2], axis[3], angle)
    elseif normal[3] < 0.0
        occ.rotate([(3, solid)], 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, Float64(pi))
    end
    occ.translate([(3, solid)], origin[1], origin[2], origin[3])
    return solid
end

# Intersect the single tube volume with the coupon box (face ends) and with the half-space
# behind every tilted joint's radial plane (design A3 (2)); exactly one volume must remain.
function trim_tube_volume!(occ, tube::AbstractTube, volume, box)
    if !isempty(tube.face_ends)
        lower, upper = box
        coupon_box = occ.addBox(
            lower[1],
            lower[2],
            lower[3],
            upper[1] - lower[1],
            upper[2] - lower[2],
            upper[3] - lower[3]
        )
        trimmed, _ = occ.intersect(volume, [(3, coupon_box)], -1, true, true)
        trimmed = [(dim, tag) for (dim, tag) in trimmed if dim == 3]
        length(trimmed) == 1 || error(
            "the box intersection of a face-ended tube left $(length(trimmed)) volumes"
        )
        volume = trimmed
    end
    for joint in tube.joints
        joint.tilt > 0.0 || continue
        lower, upper = box
        reach = 4.0 * norm(collect(upper) .- collect(lower))
        half_space = add_half_space!(occ, joint.origin, joint.normal, reach)
        trimmed, _ = occ.intersect(volume, [(3, half_space)], -1, true, true)
        trimmed = [(dim, tag) for (dim, tag) in trimmed if dim == 3]
        length(trimmed) == 1 ||
            error("the radial-plane cut of a jointed tube left $(length(trimmed)) volumes")
        volume = trimmed
    end
    return volume
end

function add_tube_volumes!(occ, tube::EdgeTube, section::TubeSection; box=nothing)
    uw = section_coordinates(section)
    K = ring_count(section)
    volumes = Tuple{Tuple{Int32, Int32}, Tuple{Int, Int, Int}}[]
    cad_start, cad_end = tube_cad_interval(tube)
    span = cad_end - cad_start
    plain_ends(tube) ||
        box !== nothing ||
        error("a face-ended or jointed tube needs the coupon box for its intersection")
    for group in material_groups(section)
        first, last, _ = group
        corners = [tube_point(tube, 0.0, 0.0, cad_start)]
        for j = first:(last + 1)
            uwj = uw[:, section_node(section, K, j)]
            push!(corners, tube_point(tube, uwj[1], uwj[2], cad_start))
        end
        points = [occ.addPoint(c[1], c[2], c[3]) for c in corners]
        lines = [
            occ.addLine(points[i], points[i % length(points) + 1]) for
            i in eachindex(points)
        ]
        loop = occ.addCurveLoop(lines)
        face = occ.addPlaneSurface([loop])
        extruded =
            occ.extrude([(2, face)], span * tube.e[1], span * tube.e[2], span * tube.e[3])
        volume = [(dim, tag) for (dim, tag) in extruded if dim == 3]
        length(volume) == 1 || error("tube extrusion produced $(length(volume)) volumes")
        volume = trim_tube_volume!(occ, tube, volume, box)
        push!(volumes, (volume[1], group))
    end
    return volumes
end

# The revolved volumes of an arc tube (design 1.2 (3)): each material-group sector face is
# drawn in the radial plane at the CAD start angle and revolved about the vertical axis
# through the centre by the CAD sweep (over-long at a face end, then box-intersected).
function add_tube_volumes!(occ, tube::ArcTube, section::TubeSection; box=nothing)
    uw = section_coordinates(section)
    K = ring_count(section)
    volumes = Tuple{Tuple{Int32, Int32}, Tuple{Int, Int, Int}}[]
    cad_start, cad_end = tube_cad_interval(tube)
    angle = tube.orientation * (cad_end - cad_start) / tube.rho
    abs(angle) < pi || error(
        "an arc tube revolves by $(rad2deg(abs(angle))) degrees: split the arc into parts below 180"
    )
    plain_ends(tube) ||
        box !== nothing ||
        error("a face-ended or jointed tube needs the coupon box for its intersection")
    for group in material_groups(section)
        first, last, _ = group
        corners = [tube_point(tube, 0.0, 0.0, cad_start)]
        for j = first:(last + 1)
            uwj = uw[:, section_node(section, K, j)]
            push!(corners, tube_point(tube, uwj[1], uwj[2], cad_start))
        end
        points = [occ.addPoint(c[1], c[2], c[3]) for c in corners]
        lines = [
            occ.addLine(points[i], points[i % length(points) + 1]) for
            i in eachindex(points)
        ]
        loop = occ.addCurveLoop(lines)
        face = occ.addPlaneSurface([loop])
        revolved = occ.revolve(
            [(2, face)],
            tube.centre[1],
            tube.centre[2],
            tube.centre[3],
            0.0,
            0.0,
            1.0,
            angle
        )
        volume = [(dim, tag) for (dim, tag) in revolved if dim == 3]
        length(volume) == 1 || error("tube revolve produced $(length(volume)) volumes")
        volume = trim_tube_volume!(occ, tube, volume, box)
        push!(volumes, (volume[1], group))
    end
    return volumes
end

# Structural entities of one tube volume (material group) with their centroids,
# used to match the fragmented OCC entities: kind, identifiers, centroid.
struct TubeEntity
    dim::Int
    kind::Symbol
    id::Tuple{Int, Int}
    centroid::Vector{Float64}
end

function polygon_centroid(points)
    area = 0.0
    cx = 0.0
    cy = 0.0
    m = length(points)
    for i = 1:m
        x1, y1 = points[i]
        x2, y2 = points[i % m + 1]
        w = x1 * y2 - x2 * y1
        area += w
        cx += (x1 + x2) * w
        cy += (y1 + y2) * w
    end
    area *= 0.5
    return (cx / (6area), cy / (6area))
end

# Area centroid of a planar polygon given by its 3D vertices (fan triangulation).
function polygon_centroid_3d(points)
    area = 0.0
    centroid = zeros(3)
    for i = 2:(length(points) - 1)
        a = cross(points[i] .- points[1], points[i + 1] .- points[1])
        weight = 0.5 * norm(a)
        area += weight
        centroid .+= weight .* (points[1] .+ points[i] .+ points[i + 1]) ./ 3.0
    end
    area > 0.0 || error("degenerate tube polygon")
    return centroid ./ area
end

# Structural entities of a FACE-ENDED tube volume (design A2 (5)): the end cap at a
# face end is the planar polygon on the face plane, the trimmed lateral and radial
# faces are planar trapezoids and the longitudinal curves end on the face, so every
# centroid is computed from the actual vertex positions (tube_end_station) instead
# of the perpendicular-end formulas of tube_entities; the fail-closed matching is the
# same.
function face_ended_tube_entities(tube::EdgeTube, section::TubeSection, group)
    first, last, _ = group
    uw = section_coordinates(section)
    K = ring_count(section)
    rays = collect(first:(last + 1))
    entities = TubeEntity[]
    outer(j) = uw[:, section_node(section, K, j)]
    at(end_index, u, w) = tube_end_point(tube, end_index, u, w)
    for end_index in (0, 1)
        push!(entities, TubeEntity(0, :edge_point, (0, end_index), at(end_index, 0.0, 0.0)))
        for j in rays
            push!(
                entities,
                TubeEntity(
                    0,
                    :outer_point,
                    (j, end_index),
                    at(end_index, outer(j)[1], outer(j)[2])
                )
            )
        end
        for j = first:last
            a = at(end_index, outer(j)[1], outer(j)[2])
            b = at(end_index, outer(j + 1)[1], outer(j + 1)[2])
            push!(entities, TubeEntity(1, :cap_polygon, (j, end_index), 0.5 .* (a .+ b)))
        end
        for j in (first, last + 1)
            a = at(end_index, 0.0, 0.0)
            b = at(end_index, outer(j)[1], outer(j)[2])
            push!(entities, TubeEntity(1, :cap_ray, (j, end_index), 0.5 .* (a .+ b)))
        end
        polygon = vcat(
            [at(end_index, 0.0, 0.0)],
            [at(end_index, outer(j)[1], outer(j)[2]) for j in rays]
        )
        push!(entities, TubeEntity(2, :cap, (0, end_index), polygon_centroid_3d(polygon)))
    end
    push!(
        entities,
        TubeEntity(1, :edge_line, (0, 0), 0.5 .* (at(0, 0.0, 0.0) .+ at(1, 0.0, 0.0)))
    )
    for j in rays
        push!(
            entities,
            TubeEntity(
                1,
                :outer_line,
                (j, 0),
                0.5 .* (at(0, outer(j)[1], outer(j)[2]) .+ at(1, outer(j)[1], outer(j)[2]))
            )
        )
    end
    for j = first:last
        a = outer(j)
        b = outer(j + 1)
        push!(
            entities,
            TubeEntity(
                2,
                :lateral,
                (j, 0),
                polygon_centroid_3d([
                    at(0, a[1], a[2]),
                    at(0, b[1], b[2]),
                    at(1, b[1], b[2]),
                    at(1, a[1], a[2])
                ])
            )
        )
    end
    for j in (first, last + 1)
        a = outer(j)
        push!(
            entities,
            TubeEntity(
                2,
                :radial,
                (j, 0),
                polygon_centroid_3d([
                    at(0, 0.0, 0.0),
                    at(0, a[1], a[2]),
                    at(1, a[1], a[2]),
                    at(1, 0.0, 0.0)
                ])
            )
        )
    end
    return entities
end

# Gauss-Legendre nodes and weights on [0, 1] (8 points) for the revolved-face centroids.
const GAUSS_LEGENDRE_8 = let
    x = [
        -0.9602898564975363,
        -0.7966664774136267,
        -0.5255324099163290,
        -0.1834346424956498,
        0.1834346424956498,
        0.5255324099163290,
        0.7966664774136267,
        0.9602898564975363
    ]
    w = [
        0.1012285362903763,
        0.2223810344533745,
        0.3137066458778873,
        0.3626837833783620,
        0.3626837833783620,
        0.3137066458778873,
        0.2223810344533745,
        0.1012285362903763
    ]
    (0.5 .* (x .+ 1.0), 0.5 .* w)
end

# Centroid of the circle arc of the arc tube's node (u, w) between the axis coordinates
# s_a and s_b: centre + radius sinc(dtheta / 2) radial(theta_mid) + w b (exact).
function revolved_curve_centroid(tube::ArcTube, u, w, s_a, s_b)
    radius = tube.rho + tube.sigma * u
    theta_a = arc_angle(tube, s_a)
    theta_b = arc_angle(tube, s_b)
    half = 0.5 * (theta_b - theta_a)
    sinc = abs(half) > 0.0 ? sin(half) / half : 1.0
    return tube.centre .+ radius * sinc .* arc_radial(0.5 * (theta_a + theta_b)) .+
           w .* tube.b
end

# Area centroid of the surface swept by the section segment (u_1, w_1) -> (u_2, w_2) of the
# arc tube between its two end sections (per node: tube_end_station, so a face-ended end's
# trimmed patch is integrated exactly up to the quadrature): dA = radius(t) L |dtheta| dt,
# the angular integral closed-form, the segment integral by Gauss-Legendre (design A7
# MINOR-6: exact moments of the section, no coincidence of centre of mass needed).
function revolved_face_centroid(tube::ArcTube, uw_1, uw_2)
    nodes, weights = GAUSS_LEGENDRE_8
    area = 0.0
    moment = zeros(3)
    for (t, weight) in zip(nodes, weights)
        u = uw_1[1] + t * (uw_2[1] - uw_1[1])
        w = uw_1[2] + t * (uw_2[2] - uw_1[2])
        radius = tube.rho + tube.sigma * u
        theta_a = arc_angle(tube, tube_end_station(tube, 0, u, w))
        theta_b = arc_angle(tube, tube_end_station(tube, 1, u, w))
        sweep = abs(theta_b - theta_a)
        radial_integral =
            sign(theta_b - theta_a) .*
            [sin(theta_b) - sin(theta_a), -(cos(theta_b) - cos(theta_a)), 0.0]
        area += weight * radius * sweep
        moment .+=
            weight .* radius .*
            (tube.centre .* sweep .+ radius .* radial_integral .+ (w * sweep) .* tube.b)
    end
    area > 0.0 || error("degenerate revolved tube face")
    return moment ./ area
end

# The cap of an ARC tube on a box face (round 2b, decision 437 (3); design A2 for arcs): the
# face plane cuts the revolved lateral faces in conics, so the cap curves and the cap face are
# the images of the section's straight segments and triangles under the map (u, w) ->
# P(u, w) = tube_end_point (the node circle of radius rho + sigma u meets the face at its own
# angle). On the face plane the in-plane coordinates are (O(u), w) with O the face-parallel
# horizontal coordinate, (O - c_o)^2 = r(u)^2 - d^2, so dO / du = sigma r(u) / (O - c_o): the
# arc-length centroid of a cap curve and the area centroid of the cap face follow by
# Gauss-Legendre quadrature (8-point; composite over the segment / the sector triangles).
function face_cap_stretch(tube::ArcTube, face_end::FaceEnd, u, w)
    point = tube_end_point(tube, face_end.end_index, u, w)
    other = 3 - face_end.face_axis
    offset = point[other] - tube.centre[other]
    abs(offset) > 0.0 || error("arc tube face cap is tangent to the face at a node")
    return point, tube.sigma * (tube.rho + tube.sigma * u) / offset
end

function face_cap_curve_centroid(tube::ArcTube, face_end::FaceEnd, a, b)
    nodes, weights = GAUSS_LEGENDRE_8
    moment = zeros(3)
    length = 0.0
    du, dw = b[1] - a[1], b[2] - a[2]
    for (t, weight) in zip(nodes, weights)
        point, stretch = face_cap_stretch(tube, face_end, a[1] + t * du, a[2] + t * dw)
        speed = hypot(stretch * du, dw)
        moment .+= weight .* speed .* point
        length += weight * speed
    end
    length > 0.0 || error("degenerate arc tube face cap curve")
    return moment ./ length
end

function face_cap_face_centroid(tube::ArcTube, face_end::FaceEnd, triangles)
    nodes, weights = GAUSS_LEGENDRE_8
    moment = zeros(3)
    area = 0.0
    for (p, q, r) in triangles
        # Duffy map of the unit square onto the (u, w) triangle p q r: (x, y) -> p + x (q - p) +
        # x y (r - q), Jacobian x |(q - p) x (r - p)|.
        twice = abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1]))
        for (x, wx) in zip(nodes, weights), (y, wy) in zip(nodes, weights)
            u = p[1] + x * (q[1] - p[1]) + x * y * (r[1] - q[1])
            w = p[2] + x * (q[2] - p[2]) + x * y * (r[2] - q[2])
            point, stretch = face_cap_stretch(tube, face_end, u, w)
            weight = wx * wy * x * twice * abs(stretch)
            moment .+= weight .* point
            area += weight
        end
    end
    area > 0.0 || error("degenerate arc tube face cap")
    return moment ./ area
end

# Structural entities of an ARC tube volume (design 1.2 (3)): points and cap curves / faces
# on the two end sections (radial planes; at a face end the box face plane, where the cap
# curves are conics and the cap face a conic-bounded region: face_cap_curve_centroid /
# face_cap_face_centroid), the longitudinal curves as circle arcs (exact centroids) and the
# lateral / radial faces as surfaces of revolution (revolved_face_centroid).
function tube_entities(tube::ArcTube, section::TubeSection, group)
    first, last, _ = group
    uw = section_coordinates(section)
    K = ring_count(section)
    rays = collect(first:(last + 1))
    entities = TubeEntity[]
    outer(j) = uw[:, section_node(section, K, j)]
    at(end_index, u, w) = tube_end_point(tube, end_index, u, w)
    for end_index in (0, 1)
        face_end = face_end_at(tube, end_index)
        curve_centroid(a_uw, b_uw) =
            face_end === nothing ?
            0.5 .* (at(end_index, a_uw...) .+ at(end_index, b_uw...)) :
            face_cap_curve_centroid(tube, face_end, a_uw, b_uw)
        push!(entities, TubeEntity(0, :edge_point, (0, end_index), at(end_index, 0.0, 0.0)))
        for j in rays
            push!(
                entities,
                TubeEntity(0, :outer_point, (j, end_index), at(end_index, outer(j)...))
            )
        end
        for j = first:last
            push!(
                entities,
                TubeEntity(
                    1,
                    :cap_polygon,
                    (j, end_index),
                    curve_centroid(outer(j), outer(j + 1))
                )
            )
        end
        for j in (first, last + 1)
            push!(
                entities,
                TubeEntity(
                    1,
                    :cap_ray,
                    (j, end_index),
                    curve_centroid((0.0, 0.0), outer(j))
                )
            )
        end
        cap_centroid = if face_end === nothing
            polygon_centroid_3d(
                vcat([at(end_index, 0.0, 0.0)], [at(end_index, outer(j)...) for j in rays])
            )
        else
            face_cap_face_centroid(
                tube,
                face_end,
                [((0.0, 0.0), outer(j), outer(j + 1)) for j = first:last]
            )
        end
        push!(entities, TubeEntity(2, :cap, (0, end_index), cap_centroid))
    end
    s_a(u, w) = tube_end_station(tube, 0, u, w)
    s_b(u, w) = tube_end_station(tube, 1, u, w)
    push!(
        entities,
        TubeEntity(
            1,
            :edge_line,
            (0, 0),
            revolved_curve_centroid(tube, 0.0, 0.0, s_a(0.0, 0.0), s_b(0.0, 0.0))
        )
    )
    for j in rays
        u, w = outer(j)
        push!(
            entities,
            TubeEntity(
                1,
                :outer_line,
                (j, 0),
                revolved_curve_centroid(tube, u, w, s_a(u, w), s_b(u, w))
            )
        )
    end
    for j = first:last
        push!(
            entities,
            TubeEntity(
                2,
                :lateral,
                (j, 0),
                revolved_face_centroid(tube, outer(j), outer(j + 1))
            )
        )
    end
    for j in (first, last + 1)
        push!(
            entities,
            TubeEntity(
                2,
                :radial,
                (j, 0),
                revolved_face_centroid(tube, [0.0, 0.0], outer(j))
            )
        )
    end
    return entities
end

function tube_entities(tube::EdgeTube, section::TubeSection, group)
    plain_ends(tube) || return face_ended_tube_entities(tube, section, group)
    first, last, _ = group
    uw = section_coordinates(section)
    K = ring_count(section)
    rays = collect(first:(last + 1))
    mid = 0.5 * (tube.s_start + tube.s_end)
    entities = TubeEntity[]
    outer(j) = uw[:, section_node(section, K, j)]
    # Points: edge ends and outer polygon vertices at both caps (end index 0 the
    # start cap, 1 the end cap).
    for (end_index, s) in ((0, tube.s_start), (1, tube.s_end))
        push!(
            entities,
            TubeEntity(0, :edge_point, (0, end_index), tube_point(tube, 0.0, 0.0, s))
        )
        for j in rays
            push!(
                entities,
                TubeEntity(
                    0,
                    :outer_point,
                    (j, end_index),
                    tube_point(tube, outer(j)[1], outer(j)[2], s)
                )
            )
        end
        # Cap curves: polygon edges and the two bounding rays.
        for j = first:last
            a = outer(j)
            b = outer(j + 1)
            push!(
                entities,
                TubeEntity(
                    1,
                    :cap_polygon,
                    (j, end_index),
                    tube_point(tube, 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), s)
                )
            )
        end
        for j in (first, last + 1)
            push!(
                entities,
                TubeEntity(
                    1,
                    :cap_ray,
                    (j, end_index),
                    tube_point(tube, 0.5 * outer(j)[1], 0.5 * outer(j)[2], s)
                )
            )
        end
        # Cap face.
        polygon = vcat([(0.0, 0.0)], [(outer(j)[1], outer(j)[2]) for j in rays])
        c = polygon_centroid(polygon)
        push!(
            entities,
            TubeEntity(2, :cap, (0, end_index), tube_point(tube, c[1], c[2], s))
        )
    end
    # Longitudinal curves: the edge and the outer polygon vertices.
    push!(entities, TubeEntity(1, :edge_line, (0, 0), tube_point(tube, 0.0, 0.0, mid)))
    for j in rays
        push!(
            entities,
            TubeEntity(
                1,
                :outer_line,
                (j, 0),
                tube_point(tube, outer(j)[1], outer(j)[2], mid)
            )
        )
    end
    # Lateral faces and the two radial faces.
    for j = first:last
        a = outer(j)
        b = outer(j + 1)
        push!(
            entities,
            TubeEntity(
                2,
                :lateral,
                (j, 0),
                tube_point(tube, 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), mid)
            )
        )
    end
    for j in (first, last + 1)
        push!(
            entities,
            TubeEntity(
                2,
                :radial,
                (j, 0),
                tube_point(tube, 0.5 * outer(j)[1], 0.5 * outer(j)[2], mid)
            )
        )
    end
    return entities
end

# Match the fragmented OCC entities bounding a tube volume to the structural ones
# by centroid. Every CAD entity of the volume must be matched (otherwise the
# fragment split a tube entity) and every structural entity must be found.
function match_tube_entities(
    volume,
    tube::AbstractTube,
    section::TubeSection,
    group,
    tolerance
)
    entities = tube_entities(tube, section, group)
    cad = Dict{Int, Vector{Int32}}(0 => Int32[], 1 => Int32[], 2 => Int32[])
    faces =
        [tag for (dim, tag) in gmsh.model.getBoundary([(3, volume)], false, false, false)]
    cad[2] = unique(abs.(faces))
    for face in cad[2]
        for (dim, curve) in gmsh.model.getBoundary([(2, face)], false, false, false)
            push!(cad[1], abs(curve))
        end
    end
    unique!(cad[1])
    for curve in cad[1]
        _, points = gmsh.model.getAdjacencies(1, curve)
        append!(cad[0], points)
    end
    unique!(cad[0])
    matched = Dict{Tuple{Symbol, Tuple{Int, Int}}, Int32}()
    used = Set{Tuple{Int, Int32}}()
    for entity in entities
        candidates = cad[entity.dim]
        distances = [
            norm(
                collect(
                    entity.dim == 0 ? gmsh.model.getValue(0, tag, Float64[]) :
                    gmsh.model.occ.getCenterOfMass(entity.dim, tag)
                ) .- entity.centroid
            ) for tag in candidates
        ]
        isempty(distances) &&
            error("tube volume $volume has no dimension-$(entity.dim) entities")
        best = argmin(distances)
        distances[best] <= tolerance || error(
            "tube $(entity.kind) $(entity.id) not found among the CAD entities of " *
            "volume $volume (nearest $(distances[best]); $(tube_description(tube)), group $group, " *
            "expected centroid $(entity.centroid))"
        )
        matched[(entity.kind, entity.id)] = candidates[best]
        push!(used, (entity.dim, candidates[best]))
    end
    for dim = 0:2, tag in cad[dim]
        (dim, tag) in used || error(
            "CAD entity ($dim, $tag) of tube volume $volume is not a tube entity: " *
            "the fragment split the tube"
        )
    end
    length(used) == length(entities) || error("tube entity matching is not one-to-one")
    return matched
end

# One line naming a tube in an error message.
tube_description(tube::EdgeTube) =
    "straight tube at $(tube.origin) along $(tube.e), s $(tube.s_start)..$(tube.s_end)" *
    (isempty(tube.joints) ? "" : ", joints $([j.end_index for j in tube.joints])") *
    (isempty(tube.face_ends) ? "" : ", face ends $([f.end_index for f in tube.face_ends])")
tube_description(tube::ArcTube) =
    "arc tube about $(tube.centre) rho $(tube.rho) sigma $(tube.sigma) theta0 $(rad2deg(tube.theta0)) deg " *
    "orientation $(tube.orientation), s $(tube.s_start)..$(tube.s_end)" *
    (isempty(tube.joints) ? "" : ", joints $([j.end_index for j in tube.joints])") *
    (isempty(tube.face_ends) ? "" : ", face ends $([f.end_index for f in tube.face_ends])")

# Explicit tube mesh state. Gmsh's GenerateMesh deletes every face mesh before
# its 1D pass and every volume mesh before its 2D and 3D passes, and
# Mesh.MeshOnlyEmpty protects only entities of the dimension being meshed, so the
# tube mesh is installed in three phases around the generator:
#   1. install_tube_curves!  (points, curves)      -> generate(2)
#   2. install_tube_faces!   (clear + faces)       -> remove_tube_volumes!, generate(3)
#   3. finalize_tube_volumes! (discrete volumes, interior nodes, prisms)
# The OCC tube volumes are removed (faces stay) before generate(3) because the
# Delaunay pyramid closure of the lateral quadrangles refuses quadrangles whose
# far side belongs to another model volume; discrete volumes replace them.
mutable struct TubeMesh
    tube::AbstractTube
    section::TubeSection
    volumes::Vector{Tuple{Int32, Tuple{Int, Int, Int}, Dict}}   # (OCC volume, group, matched)
    pyramid_height::Float64                                     # apex distance above each lateral quad
    tags::Dict{Tuple{Int, Int}, Int}                            # (local section node, station) -> node tag
    coordinates::Dict{Int, Vector{Float64}}
    apex::Dict{Tuple{Int, Int}, Int}                            # (sector j, station i) -> apex node tag
    faces::Dict{Int32, Vector{Int32}}                           # OCC volume -> boundary faces
    discrete::Dict{Int32, Int32}                                # OCC volume -> discrete volume
    joints::Vector{Tuple{Int, TubeMesh, Int}}                   # (own end, the owning state, its end)
end

# A smooth joint (design A3 (2)): the end section `end_index` of `state` IS the end section
# `owner_end` of `owner` (the arc tube, or the earlier part of one arc). `owner` is installed
# before `state` in every phase, so adopt_joint_tags! copies the owner's section nodes that
# exist at that point (the CAD points after the curves phase, the ray nodes, the cap interior
# nodes after the faces phase) into the state's own tags; the shared CAD entities are meshed
# once (the shared `meshed` sets of the install phases).
function register_joint!(state::TubeMesh, end_index, owner::TubeMesh, owner_end)
    end_index in (0, 1) && owner_end in (0, 1) || error("a joint end is 0 or 1")
    owner !== state || error("a tube cannot joint itself")
    state.section.ring_radii == owner.section.ring_radii &&
    state.section.angles == owner.section.angles ||
        error("a smooth joint needs the same tube section on both tubes")
    push!(state.joints, (end_index, owner, owner_end))
    return state
end

function adopt_joint_tags!(state::TubeMesh)
    for (end_index, owner, owner_end) in state.joints
        own_station = end_index == 0 ? 0 : state.tube.layers
        owner_station = owner_end == 0 ? 0 : owner.tube.layers
        for ((local_index, station), tag) in owner.tags
            station == owner_station || continue
            haskey(state.tags, (local_index, own_station)) && continue
            state.tags[(local_index, own_station)] = tag
            state.coordinates[tag] = owner.coordinates[tag]
        end
    end
    return state
end

# The tube's outer surface presented to the tetrahedral mesher is made of
# triangles: every lateral quadrangle (outer prism face) carries an explicit
# pyramid whose apex lies pyramid_height outside the quad on the sector bisector,
# so Gmsh meshes pure tetrahedra against a closed triangulated boundary (its
# tetrahedral optimizer was observed to create overlapping tetrahedra when it
# had to close quadrangles with its own pyramids).
function TubeMesh(tube::AbstractTube, section::TubeSection, volumes; pyramid_height)
    pyramid_height > 0.0 || error("pyramid height must be positive")
    faces = Dict{Int32, Vector{Int32}}()
    for (volume, _, _) in volumes
        faces[volume] = unique(
            abs(tag) for
            (dim, tag) in gmsh.model.getBoundary([(3, volume)], false, false, false)
        )
    end
    return TubeMesh(
        tube,
        section,
        volumes,
        pyramid_height,
        Dict{Tuple{Int, Int}, Int}(),
        Dict{Int, Vector{Float64}}(),
        Dict{Tuple{Int, Int}, Int}(),
        faces,
        Dict{Int32, Int32}(),
        Tuple{Int, TubeMesh, Int}[]
    )
end

function apex_node!(state::TubeMesh, next_node, j, i)
    get!(state.apex, (j, i)) do
        theta = deg2rad(0.5 * (state.section.angles[j + 1] + state.section.angles[j + 2]))
        radius =
            tube_radius(state.section) *
            cos(0.5 * deg2rad(state.section.angles[j + 2] - state.section.angles[j + 1])) +
            state.pyramid_height
        next_node[] += 1
        u = radius * cos(theta)
        w = radius * sin(theta)
        state.coordinates[next_node[]] =
            tube_point(state.tube, u, w, tube_station(state.tube, i + 0.5, u, w))
        return next_node[]
    end
end

function tube_node!(state::TubeMesh, next_node, local_index, i)
    get!(state.tags, (local_index, i)) do
        uw = section_coordinates(state.section)
        next_node[] += 1
        state.coordinates[next_node[]] =
            tube_section_point(state.tube, i, uw[1, local_index], uw[2, local_index])
        return next_node[]
    end
end

# Phase 1: nodes and line elements of the tube's CAD points and curves. Node tags
# are shared with the caller's explicit curve meshes through next_node and
# point_nodes (CAD point tag -> node tag). `meshed` (shared across the tubes of one
# build) keeps a CAD curve two jointed tubes share meshed once.
function install_tube_curves!(state::TubeMesh, next_node, point_nodes; meshed=Set{Int32}())
    adopt_joint_tags!(state)
    section = state.section
    K = ring_count(section)
    L = state.tube.layers
    lines = Dict{Int32, Vector{Int}}()
    curve_nodes = Dict{Int32, Vector{Int}}()
    function curve_node!(curve, local_index, i)
        haskey(state.tags, (local_index, i)) && return state.tags[(local_index, i)]
        t = tube_node!(state, next_node, local_index, i)
        push!(get!(curve_nodes, curve, Int[]), t)
        return t
    end
    for (_, group, matched) in state.volumes
        first, last, _ = group
        rays = first:(last + 1)
        for (end_index, i) in ((0, 0), (1, L))
            for (kind, id, local_index) in vcat(
                [(:edge_point, (0, end_index), section_node(section, 0, 0))],
                [(:outer_point, (j, end_index), section_node(section, K, j)) for j in rays]
            )
                point = matched[(kind, id)]
                tag = get!(point_nodes, point) do
                    t = tube_node!(state, next_node, local_index, i)
                    gmsh.model.mesh.addNodes(0, point, [t], state.coordinates[t])
                    return t
                end
                state.tags[(local_index, i)] = tag
            end
        end
        for (kind, id, local_index) in vcat(
            [(:edge_line, (0, 0), section_node(section, 0, 0))],
            [(:outer_line, (j, 0), section_node(section, K, j)) for j in rays]
        )
            curve = matched[(kind, id)]
            curve in meshed && continue
            push!(meshed, curve)
            for i = 0:(L - 1)
                append!(
                    get!(lines, curve, Int[]),
                    (
                        curve_node!(curve, local_index, i),
                        curve_node!(curve, local_index, i + 1)
                    )
                )
            end
        end
        for (end_index, i) in ((0, 0), (1, L))
            for j = first:last
                curve = matched[(:cap_polygon, (j, end_index))]
                curve in meshed && continue
                push!(meshed, curve)
                append!(
                    get!(lines, curve, Int[]),
                    (
                        state.tags[(section_node(section, K, j), i)],
                        state.tags[(section_node(section, K, j + 1), i)]
                    )
                )
            end
            for j in (first, last + 1)
                curve = matched[(:cap_ray, (j, end_index))]
                curve in meshed && continue
                push!(meshed, curve)
                for k = 1:K
                    append!(
                        get!(lines, curve, Int[]),
                        (
                            curve_node!(curve, section_node(section, k - 1, j), i),
                            curve_node!(curve, section_node(section, k, j), i)
                        )
                    )
                end
            end
        end
    end
    for (curve, nodes) in curve_nodes
        parameters = [
            gmsh.model.getParametrization(1, curve, state.coordinates[t])[1] for t in nodes
        ]
        gmsh.model.mesh.addNodes(
            1,
            curve,
            nodes,
            reduce(vcat, state.coordinates[t] for t in nodes),
            parameters
        )
    end
    for (curve, connectivity) in lines
        gmsh.model.mesh.addElementsByType(curve, 1, Int[], connectivity)
    end
    return Dict{String, Any}(
        "Curves" => length(lines),
        "Lines" => sum(length(v) ÷ 2 for v in values(lines); init=0)
    )
end

# Phase 2 (after generate(2)): replace Gmsh's triangulation of the tube faces by
# the explicit cap triangles and lateral/radial quadrangles.
function install_tube_faces!(state::TubeMesh, next_node; meshed=Set{Int32}())
    adopt_joint_tags!(state)
    section = state.section
    K = ring_count(section)
    L = state.tube.layers
    all_faces = unique(reduce(vcat, values(state.faces); init=Int32[]))
    # A face shared with a jointed tube installed earlier keeps that tube's mesh.
    gmsh.model.mesh.clear([(2, face) for face in all_faces if !(face in meshed)])
    # generate(2) numbered its nodes after the explicit curve nodes.
    next_node[] = max(next_node[], Int(gmsh.model.mesh.getMaxNodeTag()))
    face_nodes = Dict{Int32, Vector{Int}}()
    triangles = Dict{Int32, Vector{Int}}()
    quads = Dict{Int32, Vector{Int}}()
    function face_node!(face, local_index, i)
        haskey(state.tags, (local_index, i)) && return state.tags[(local_index, i)]
        t = tube_node!(state, next_node, local_index, i)
        push!(get!(face_nodes, face, Int[]), t)
        return t
    end
    for (_, group, matched) in state.volumes
        first, last, _ = group
        for (end_index, i) in ((0, 0), (1, L))
            face = matched[(:cap, (0, end_index))]
            face in meshed && continue
            push!(meshed, face)
            for j = first:last, (a, b, c) in sector_triangles(section, j)
                append!(
                    get!(triangles, face, Int[]),
                    (
                        face_node!(face, a, i),
                        face_node!(face, b, i),
                        face_node!(face, c, i)
                    )
                )
            end
        end
        for j = first:last
            face = matched[(:lateral, (j, 0))]
            push!(meshed, face)
            a = section_node(section, K, j)
            b = section_node(section, K, j + 1)
            for i = 0:(L - 1)
                apex = apex_node!(state, next_node, j, i)
                push!(get!(face_nodes, face, Int[]), apex)
                a0, b0, b1, a1 = state.tags[(a, i)],
                state.tags[(b, i)],
                state.tags[(b, i + 1)],
                state.tags[(a, i + 1)]
                append!(
                    get!(triangles, face, Int[]),
                    (a0, b0, apex, b0, b1, apex, b1, a1, apex, a1, a0, apex)
                )
            end
        end
        for j in (first, last + 1)
            face = matched[(:radial, (j, 0))]
            face in meshed && continue
            push!(meshed, face)
            for k = 1:K, i = 0:(L - 1)
                a = section_node(section, k - 1, j)
                b = section_node(section, k, j)
                append!(
                    get!(quads, face, Int[]),
                    (
                        face_node!(face, a, i),
                        face_node!(face, b, i),
                        face_node!(face, b, i + 1),
                        face_node!(face, a, i + 1)
                    )
                )
            end
        end
    end
    all(face in meshed for face in all_faces) ||
        error("tube faces and matched faces differ")
    for (face, nodes) in face_nodes
        gmsh.model.mesh.addNodes(
            2,
            face,
            nodes,
            reduce(vcat, state.coordinates[t] for t in nodes)
        )
    end
    for (face, connectivity) in triangles
        gmsh.model.mesh.addElementsByType(face, 2, Int[], connectivity)
    end
    for (face, connectivity) in quads
        gmsh.model.mesh.addElementsByType(face, 3, Int[], connectivity)
    end
    return Dict{String, Any}(
        "Faces" => length(all_faces),
        "PyramidApexes" => length(state.apex),
        "Triangles" => sum(length(v) ÷ 3 for v in values(triangles); init=0),
        "Quads" => sum(length(v) ÷ 4 for v in values(quads); init=0)
    )
end

# Remove the OCC tube volumes (non-recursively: faces, curves, points and their
# meshes stay) so that the remaining volumes can be meshed with pyramids.
function remove_tube_volumes!(states::Vector{TubeMesh})
    return gmsh.model.removeEntities(
        [(3, volume) for state in states for (volume, _, _) in state.volumes],
        false
    )
end

# Phase 3 (after generate(3)): discrete volumes with the interior nodes and the
# prisms. Returns the discrete volume tags per material and the census.
function finalize_tube_volumes!(states::Vector{TubeMesh})
    volumes = Dict{Int, Vector{Int32}}()
    census = Dict{String, Any}[]
    next_node = Ref(Int(gmsh.model.mesh.getMaxNodeTag()))
    for state in states
        adopt_joint_tags!(state)
        section = state.section
        K = ring_count(section)
        L = state.tube.layers
        for (volume, group, _) in state.volumes
            first, last, material = group
            rays = first:(last + 1)
            discrete = gmsh.model.addDiscreteEntity(3, -1, state.faces[volume])
            state.discrete[volume] = discrete
            interior = Int[]
            for i = 1:(L - 1),
                local_index in vcat(
                    [section_node(section, 0, 0)],
                    [section_node(section, k, j) for k = 1:K for j in rays]
                )

                haskey(state.tags, (local_index, i)) && continue
                push!(interior, tube_node!(state, next_node, local_index, i))
            end
            isempty(interior) || gmsh.model.mesh.addNodes(
                3,
                discrete,
                interior,
                reduce(vcat, state.coordinates[t] for t in interior)
            )
            prisms = Int[]
            for i = 0:(L - 1), j = first:last, (a, b, c) in sector_triangles(section, j)
                append!(
                    prisms,
                    (
                        state.tags[(a, i)],
                        state.tags[(b, i)],
                        state.tags[(c, i)],
                        state.tags[(a, i + 1)],
                        state.tags[(b, i + 1)],
                        state.tags[(c, i + 1)]
                    )
                )
            end
            gmsh.model.mesh.addElementsByType(discrete, 6, Int[], prisms)
            pyramids = Int[]
            for i = 0:(L - 1), j = first:last
                a = section_node(section, K, j)
                b = section_node(section, K, j + 1)
                append!(
                    pyramids,
                    (
                        state.tags[(a, i)],
                        state.tags[(b, i)],
                        state.tags[(b, i + 1)],
                        state.tags[(a, i + 1)],
                        state.apex[(j, i)]
                    )
                )
            end
            gmsh.model.mesh.addElementsByType(discrete, 7, Int[], pyramids)
            push!(get!(volumes, material, Int32[]), discrete)
            push!(
                census,
                Dict{String, Any}(
                    "Volume" => Int(volume),
                    "DiscreteVolume" => Int(discrete),
                    "Material" => material,
                    "Sectors" => last - first + 1,
                    "Layers" => L,
                    "Prisms" => length(prisms) ÷ 6,
                    "Pyramids" => length(pyramids) ÷ 5,
                    "InteriorNodes" => length(interior)
                )
            )
        end
    end
    return volumes, census
end

# Build every tube volume in the OCC model (before the fragment) and return the
# records needed afterwards.
struct TubeRecord
    tube::AbstractTube
    section::TubeSection
    group::Tuple{Int, Int, Int}
    tool::Tuple{Int32, Int32}
end

# Fragment map lookup: the single volume descendant of a tube tool.
function tube_volume_after_fragment(
    record::TubeRecord,
    tool_index,
    fragment_map,
    object_count
)
    descendants =
        [tag for (dim, tag) in fragment_map[object_count + tool_index] if dim == 3]
    length(descendants) == 1 ||
        error("tube tool $(record.tool) has $(length(descendants)) volume descendants")
    return descendants[1]
end

# Anisotropy of the tube prisms: the largest longitudinal spacing over the
# smallest cross-section edge (the innermost arc), and the ring sizes.
function tube_anisotropy(tube::AbstractTube, section::TubeSection, pyramid_height)
    spacing = tube_spacing(tube)
    inner_arc = section.ring_radii[1] * deg2rad(minimum(diff(section.angles)))
    outer_arc =
        2.0 * tube_radius(section) * sin(0.5 * deg2rad(maximum(diff(section.angles))))
    return Dict{String, Any}(
        "LongitudinalSpacing" => spacing,
        "RingSizes" => ring_sizes(section),
        "OuterArc" => outer_arc,
        "PyramidHeight" => pyramid_height,
        "SectorDegrees" => diff(section.angles),
        "RingRadii" => copy(section.ring_radii),
        "InnermostArc" => inner_arc,
        "MaximumEdgeAspect" => spacing / min(inner_arc, ring_sizes(section)[1])
    )
end
