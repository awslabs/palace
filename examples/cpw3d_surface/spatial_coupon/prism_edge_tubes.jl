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
end

function TubeSection(inner_size, ratio, rings, angles, materials)
    inner_size > 0.0 || error("tube inner size must be positive")
    ratio > 1.0 || error("tube ring ratio must exceed 1")
    rings >= 1 || error("at least one ring")
    length(materials) == length(angles) - 1 || error("one material per angular sector")
    all(diff(angles) .> 0.0) || error("rays must increase")
    radii = [inner_size * (ratio^k - 1.0) / (ratio - 1.0) for k in 1:rings]
    return TubeSection(radii, collect(Float64, angles), collect(Int, materials))
end

ring_sizes(section::TubeSection) = diff(vcat(0.0, section.ring_radii))
tube_radius(section::TubeSection) = section.ring_radii[end]
ring_count(section::TubeSection) = length(section.ring_radii)
ray_count(section::TubeSection) = length(section.angles)

# Cross-section node (k, j): k = 0 the edge point (one node), k >= 1 ring k, ray j
# (0-based). Returns the 1-based local index within one cross-section.
function section_node(section::TubeSection, k, j)
    k == 0 && return 1
    return 1 + (k - 1) * ray_count(section) + j + 1
end
section_node_count(section::TubeSection) = 1 + ring_count(section) * ray_count(section)

# Local (u, w) coordinates of every cross-section node, u along n, w along b.
function section_coordinates(section::TubeSection)
    uw = zeros(2, section_node_count(section))
    for k in 1:ring_count(section), j in 0:(ray_count(section) - 1)
        theta = deg2rad(section.angles[j + 1])
        uw[:, section_node(section, k, j)] .= section.ring_radii[k] .* (cos(theta), sin(theta))
    end
    return uw
end

# Counterclockwise triangles of the cross-section per sector j (0-based):
# the fan of ring 1 and two triangles per annulus quad.
function sector_triangles(section::TubeSection, j)
    triangles = NTuple{3, Int}[]
    push!(triangles, (section_node(section, 0, 0), section_node(section, 1, j),
                      section_node(section, 1, j + 1)))
    for k in 2:ring_count(section)
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
    for j in 1:(length(section.materials) - 1)
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

# A tube: the edge line origin (a point on the edge), the frame (n, b, e), the
# extrusion interval [s_start, s_end] along e from the origin, and the layer
# boundaries (stations) s_start = stations[1] < ... < stations[layers + 1] = s_end.
struct EdgeTube
    origin::Vector{Float64}
    n::Vector{Float64}
    b::Vector{Float64}
    e::Vector{Float64}
    s_start::Float64
    s_end::Float64
    layers::Int
    stations::Vector{Float64}
end

# Uniform layers: the smallest number of equal layers whose spacing does not
# exceed `spacing` (a tube whose length is a multiple of the spacing keeps it
# exactly; tube_spacing records the largest layer actually used).
function EdgeTube(origin, n, b, s_start, s_end, spacing)
    n = collect(Float64, n) ./ norm(n)
    b = collect(Float64, b) ./ norm(b)
    abs(dot(n, b)) < 1.0e-12 || error("tube frame must be orthogonal")
    e = cross(n, b)
    s_end > s_start || error("tube interval must be increasing")
    spacing > 0.0 || error("tube spacing must be positive")
    extent = s_end - s_start
    layers = max(1, ceil(Int, extent / spacing * (1.0 - 1.0e-9)))
    stations = [s_start + extent * i / layers for i in 0:layers]
    return EdgeTube(collect(Float64, origin), n, b, e, s_start, s_end, layers, stations)
end

# The same tube with the given layer boundaries.
function EdgeTube(tube::EdgeTube, stations::AbstractVector)
    stations = collect(Float64, stations)
    length(stations) >= 2 && stations[1] == tube.s_start && stations[end] == tube.s_end ||
        error("tube stations must run from s_start to s_end")
    all(diff(stations) .> 0.0) || error("tube stations must increase")
    return EdgeTube(tube.origin, tube.n, tube.b, tube.e, tube.s_start, tube.s_end,
                    length(stations) - 1, stations)
end

tube_layer_thicknesses(tube::EdgeTube) = diff(tube.stations)
# The extrusion spacing actually used: the largest layer (equal to every layer of
# a uniform tube).
tube_spacing(tube::EdgeTube) = maximum(tube_layer_thicknesses(tube))

function tube_point(tube::EdgeTube, u, w, s)
    return tube.origin .+ u .* tube.n .+ w .* tube.b .+ s .* tube.e
end

# Axis coordinate of station i (an integer 0:layers) or of a point between two
# stations (i + 1/2: the pyramid apex station of layer i).
function tube_station(tube::EdgeTube, i)
    k = clamp(floor(Int, i), 0, tube.layers - 1)
    return tube.stations[k + 1] + (i - k) * (tube.stations[k + 2] - tube.stations[k + 1])
end

# Samples per prescribed size of the arclength quadrature that places the layer
# boundaries (a resolution constant, not a mesh target).
const TUBE_LAYER_SAMPLES_PER_SIZE = 8

# Layer boundaries following a size field along the tube axis (supervisor decision
# 40): `size_at(s)` is the composed size prescribed at axis coordinate s, capped at
# `spacing`; the field is sampled adaptively (TUBE_LAYER_SAMPLES_PER_SIZE samples
# per local size), gradient-limited along the axis to the slope
# (growth - 1) / growth, and the layers equidistribute the arclength integral of
# 1 / size with ceil(integral) layers, so every layer is at most the largest size
# it spans (strictly below spacing: the cap is spacing x (1 - 1e-9), a tube whose
# length is an exact multiple of a uniform spacing takes one more layer than the
# uniform constructor) and consecutive layers differ by a factor of at most
# exp((growth - 1) / growth) < growth (the recorded MaximumNeighbourRatio is
# checked against `growth`, fail closed). Returns the stations and the sampled
# (s, limited size) pairs for the census.
function graded_tube_stations(s_start, s_end, size_at, spacing, growth)
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
        for i in 2:length(sizes)
            sizes[i] = min(sizes[i], sizes[i - 1] + slope * (positions[i] - positions[i - 1]))
        end
        for i in (length(sizes) - 1):-1:1
            sizes[i] = min(sizes[i], sizes[i + 1] + slope * (positions[i + 1] - positions[i]))
        end
    end
    limit!()
    # The limiter lowers sizes below the sampling density they were sampled at;
    # refine the coarse intervals (midpoints at the prescribed size, then the
    # limiter again) until every interval is resolved.
    while true
        coarse = [i for i in 1:(length(positions) - 1)
                  if positions[i + 1] - positions[i] >
                     min(sizes[i], sizes[i + 1]) / TUBE_LAYER_SAMPLES_PER_SIZE * (1.0 + 1.0e-9)]
        isempty(coarse) && break
        refined_positions = Float64[]
        refined_sizes = Float64[]
        coarse_set = Set(coarse)
        for i in 1:(length(positions) - 1)
            push!(refined_positions, positions[i]); push!(refined_sizes, sizes[i])
            i in coarse_set || continue
            mid = 0.5 * (positions[i] + positions[i + 1])
            push!(refined_positions, mid); push!(refined_sizes, prescribed(mid))
        end
        push!(refined_positions, positions[end]); push!(refined_sizes, sizes[end])
        positions = refined_positions
        sizes = refined_sizes
        limit!()
    end
    cumulative = zeros(length(positions))
    for i in 2:length(positions)
        cumulative[i] = cumulative[i - 1] + (positions[i] - positions[i - 1]) *
                                            0.5 * (1.0 / sizes[i - 1] + 1.0 / sizes[i])
    end
    layers = max(1, ceil(Int, cumulative[end]))
    stations = [Float64(s_start)]
    for k in 1:(layers - 1)
        target = k * cumulative[end] / layers
        i = clamp(searchsortedlast(cumulative, target), 1, length(positions) - 1)
        fraction = (target - cumulative[i]) / (cumulative[i + 1] - cumulative[i])
        push!(stations, positions[i] + fraction * (positions[i + 1] - positions[i]))
    end
    push!(stations, Float64(s_end))
    thicknesses = diff(stations)
    all(thicknesses .> 0.0) || error("tube layer placement produced an empty layer")
    ratios = thicknesses[2:end] ./ thicknesses[1:(end - 1)]
    neighbour_ratio = isempty(ratios) ? 1.0 : max(maximum(ratios), 1.0 / minimum(ratios))
    neighbour_ratio <= growth * (1.0 + 1.0e-9) ||
        error("tube layers grow by $neighbour_ratio between neighbours, above $growth")
    return stations, positions, sizes
end

# Per-tube layer statistics against the prescribed (gradient-limited) sizes
# sampled by graded_tube_stations: thickness minimum / P50 / maximum, the
# thickness and the prescribed size at both ends, the achieved-over-prescribed
# ratio (each layer over the size at its midpoint) and the neighbour ratio.
function tube_layer_statistics(tube::EdgeTube, positions, sizes)
    thicknesses = tube_layer_thicknesses(tube)
    function size_at(s)
        i = clamp(searchsortedlast(positions, s), 1, length(positions) - 1)
        fraction = (s - positions[i]) / (positions[i + 1] - positions[i])
        return sizes[i] + fraction * (sizes[i + 1] - sizes[i])
    end
    midpoints = 0.5 .* (tube.stations[1:(end - 1)] .+ tube.stations[2:end])
    achieved = [thicknesses[i] / size_at(midpoints[i]) for i in eachindex(thicknesses)]
    ratios = thicknesses[2:end] ./ thicknesses[1:(end - 1)]
    median(values) = sort(values)[cld(length(values), 2)]
    return Dict{String, Any}(
        "Minimum" => minimum(thicknesses), "P50" => median(thicknesses),
        "Maximum" => maximum(thicknesses),
        "AtStart" => thicknesses[1], "AtEnd" => thicknesses[end],
        "PrescribedAtStart" => sizes[1], "PrescribedAtEnd" => sizes[end],
        "AchievedOverPrescribed" => Dict{String, Any}(
            "Minimum" => minimum(achieved), "P50" => median(achieved),
            "Maximum" => maximum(achieved)),
        "MaximumNeighbourRatio" => isempty(ratios) ? 1.0 :
                                   max(maximum(ratios), 1.0 / minimum(ratios)))
end

# OCC tube volumes (one polygon-sector prism per material group), returned as
# (dim, tag) pairs with their material group, before synchronization.
function add_tube_volumes!(occ, tube::EdgeTube, section::TubeSection)
    uw = section_coordinates(section)
    K = ring_count(section)
    volumes = Tuple{Tuple{Int32, Int32}, Tuple{Int, Int, Int}}[]
    span = tube.s_end - tube.s_start
    for group in material_groups(section)
        first, last, _ = group
        corners = [tube_point(tube, 0.0, 0.0, tube.s_start)]
        for j in first:(last + 1)
            uwj = uw[:, section_node(section, K, j)]
            push!(corners, tube_point(tube, uwj[1], uwj[2], tube.s_start))
        end
        points = [occ.addPoint(c[1], c[2], c[3]) for c in corners]
        lines = [occ.addLine(points[i], points[i % length(points) + 1])
                 for i in eachindex(points)]
        loop = occ.addCurveLoop(lines)
        face = occ.addPlaneSurface([loop])
        extruded = occ.extrude([(2, face)], span * tube.e[1], span * tube.e[2],
                               span * tube.e[3])
        volume = [(dim, tag) for (dim, tag) in extruded if dim == 3]
        length(volume) == 1 || error("tube extrusion produced $(length(volume)) volumes")
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
    for i in 1:m
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

function tube_entities(tube::EdgeTube, section::TubeSection, group)
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
        push!(entities, TubeEntity(0, :edge_point, (0, end_index), tube_point(tube, 0.0, 0.0, s)))
        for j in rays
            push!(entities, TubeEntity(0, :outer_point, (j, end_index),
                                       tube_point(tube, outer(j)[1], outer(j)[2], s)))
        end
        # Cap curves: polygon edges and the two bounding rays.
        for j in first:last
            a = outer(j)
            b = outer(j + 1)
            push!(entities, TubeEntity(1, :cap_polygon, (j, end_index),
                                       tube_point(tube, 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), s)))
        end
        for j in (first, last + 1)
            push!(entities, TubeEntity(1, :cap_ray, (j, end_index),
                                       tube_point(tube, 0.5 * outer(j)[1], 0.5 * outer(j)[2], s)))
        end
        # Cap face.
        polygon = vcat([(0.0, 0.0)], [(outer(j)[1], outer(j)[2]) for j in rays])
        c = polygon_centroid(polygon)
        push!(entities, TubeEntity(2, :cap, (0, end_index), tube_point(tube, c[1], c[2], s)))
    end
    # Longitudinal curves: the edge and the outer polygon vertices.
    push!(entities, TubeEntity(1, :edge_line, (0, 0), tube_point(tube, 0.0, 0.0, mid)))
    for j in rays
        push!(entities, TubeEntity(1, :outer_line, (j, 0),
                                   tube_point(tube, outer(j)[1], outer(j)[2], mid)))
    end
    # Lateral faces and the two radial faces.
    for j in first:last
        a = outer(j)
        b = outer(j + 1)
        push!(entities, TubeEntity(2, :lateral, (j, 0),
                                   tube_point(tube, 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]), mid)))
    end
    for j in (first, last + 1)
        push!(entities, TubeEntity(2, :radial, (j, 0),
                                   tube_point(tube, 0.5 * outer(j)[1], 0.5 * outer(j)[2], mid)))
    end
    return entities
end

# Match the fragmented OCC entities bounding a tube volume to the structural ones
# by centroid. Every CAD entity of the volume must be matched (otherwise the
# fragment split a tube entity) and every structural entity must be found.
function match_tube_entities(volume, tube::EdgeTube, section::TubeSection, group, tolerance)
    entities = tube_entities(tube, section, group)
    cad = Dict{Int, Vector{Int32}}(0 => Int32[], 1 => Int32[], 2 => Int32[])
    faces = [tag for (dim, tag) in gmsh.model.getBoundary([(3, volume)], false, false, false)]
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
        distances = [norm(collect(entity.dim == 0 ? gmsh.model.getValue(0, tag, Float64[]) :
                                  gmsh.model.occ.getCenterOfMass(entity.dim, tag)) .-
                          entity.centroid) for tag in candidates]
        isempty(distances) && error("tube volume $volume has no dimension-$(entity.dim) entities")
        best = argmin(distances)
        distances[best] <= tolerance ||
            error("tube $(entity.kind) $(entity.id) not found among the CAD entities of " *
                  "volume $volume (nearest $(distances[best]))")
        matched[(entity.kind, entity.id)] = candidates[best]
        push!(used, (entity.dim, candidates[best]))
    end
    for dim in 0:2, tag in cad[dim]
        (dim, tag) in used ||
            error("CAD entity ($dim, $tag) of tube volume $volume is not a tube entity: " *
                  "the fragment split the tube")
    end
    length(used) == length(entities) || error("tube entity matching is not one-to-one")
    return matched
end

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
    tube::EdgeTube
    section::TubeSection
    volumes::Vector{Tuple{Int32, Tuple{Int, Int, Int}, Dict}}   # (OCC volume, group, matched)
    pyramid_height::Float64                                     # apex distance above each lateral quad
    tags::Dict{Tuple{Int, Int}, Int}                            # (local section node, station) -> node tag
    coordinates::Dict{Int, Vector{Float64}}
    apex::Dict{Tuple{Int, Int}, Int}                            # (sector j, station i) -> apex node tag
    faces::Dict{Int32, Vector{Int32}}                           # OCC volume -> boundary faces
    discrete::Dict{Int32, Int32}                                # OCC volume -> discrete volume
end

# The tube's outer surface presented to the tetrahedral mesher is made of
# triangles: every lateral quadrangle (outer prism face) carries an explicit
# pyramid whose apex lies pyramid_height outside the quad on the sector bisector,
# so Gmsh meshes pure tetrahedra against a closed triangulated boundary (its
# tetrahedral optimizer was observed to create overlapping tetrahedra when it
# had to close quadrangles with its own pyramids).
function TubeMesh(tube::EdgeTube, section::TubeSection, volumes; pyramid_height)
    pyramid_height > 0.0 || error("pyramid height must be positive")
    faces = Dict{Int32, Vector{Int32}}()
    for (volume, _, _) in volumes
        faces[volume] = unique(abs(tag) for (dim, tag) in
                               gmsh.model.getBoundary([(3, volume)], false, false, false))
    end
    return TubeMesh(tube, section, volumes, pyramid_height, Dict{Tuple{Int, Int}, Int}(),
                    Dict{Int, Vector{Float64}}(), Dict{Tuple{Int, Int}, Int}(), faces,
                    Dict{Int32, Int32}())
end

function apex_node!(state::TubeMesh, next_node, j, i)
    get!(state.apex, (j, i)) do
        theta = deg2rad(0.5 * (state.section.angles[j + 1] + state.section.angles[j + 2]))
        radius = tube_radius(state.section) * cos(0.5 * deg2rad(state.section.angles[j + 2] -
                                                                state.section.angles[j + 1])) +
                 state.pyramid_height
        next_node[] += 1
        state.coordinates[next_node[]] = tube_point(state.tube, radius * cos(theta),
                                                    radius * sin(theta),
                                                    tube_station(state.tube, i + 0.5))
        next_node[]
    end
end

function tube_node!(state::TubeMesh, next_node, local_index, i)
    get!(state.tags, (local_index, i)) do
        uw = section_coordinates(state.section)
        next_node[] += 1
        state.coordinates[next_node[]] = tube_point(state.tube, uw[1, local_index],
                                                    uw[2, local_index],
                                                    tube_station(state.tube, i))
        next_node[]
    end
end

# Phase 1: nodes and line elements of the tube's CAD points and curves. Node tags
# are shared with the caller's explicit curve meshes through next_node and
# point_nodes (CAD point tag -> node tag).
function install_tube_curves!(state::TubeMesh, next_node, point_nodes)
    section = state.section
    K = ring_count(section)
    L = state.tube.layers
    meshed = Set{Int32}()
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
                    [(:outer_point, (j, end_index), section_node(section, K, j)) for j in rays])
                point = matched[(kind, id)]
                tag = get!(point_nodes, point) do
                    t = tube_node!(state, next_node, local_index, i)
                    gmsh.model.mesh.addNodes(0, point, [t], state.coordinates[t])
                    t
                end
                state.tags[(local_index, i)] = tag
            end
        end
        for (kind, id, local_index) in vcat(
                [(:edge_line, (0, 0), section_node(section, 0, 0))],
                [(:outer_line, (j, 0), section_node(section, K, j)) for j in rays])
            curve = matched[(kind, id)]
            curve in meshed && continue
            push!(meshed, curve)
            for i in 0:(L - 1)
                append!(get!(lines, curve, Int[]),
                        (curve_node!(curve, local_index, i), curve_node!(curve, local_index, i + 1)))
            end
        end
        for (end_index, i) in ((0, 0), (1, L))
            for j in first:last
                curve = matched[(:cap_polygon, (j, end_index))]
                curve in meshed && continue
                push!(meshed, curve)
                append!(get!(lines, curve, Int[]), (state.tags[(section_node(section, K, j), i)],
                                                    state.tags[(section_node(section, K, j + 1), i)]))
            end
            for j in (first, last + 1)
                curve = matched[(:cap_ray, (j, end_index))]
                curve in meshed && continue
                push!(meshed, curve)
                for k in 1:K
                    append!(get!(lines, curve, Int[]),
                            (curve_node!(curve, section_node(section, k - 1, j), i),
                             curve_node!(curve, section_node(section, k, j), i)))
                end
            end
        end
    end
    for (curve, nodes) in curve_nodes
        parameters = [gmsh.model.getParametrization(1, curve, state.coordinates[t])[1] for t in nodes]
        gmsh.model.mesh.addNodes(1, curve, nodes, reduce(vcat, state.coordinates[t] for t in nodes),
                                 parameters)
    end
    for (curve, connectivity) in lines
        gmsh.model.mesh.addElementsByType(curve, 1, Int[], connectivity)
    end
    return Dict{String, Any}("Curves" => length(lines),
                             "Lines" => sum(length(v) ÷ 2 for v in values(lines); init=0))
end

# Phase 2 (after generate(2)): replace Gmsh's triangulation of the tube faces by
# the explicit cap triangles and lateral/radial quadrangles.
function install_tube_faces!(state::TubeMesh, next_node)
    section = state.section
    K = ring_count(section)
    L = state.tube.layers
    all_faces = unique(reduce(vcat, values(state.faces); init=Int32[]))
    gmsh.model.mesh.clear([(2, face) for face in all_faces])
    # generate(2) numbered its nodes after the explicit curve nodes.
    next_node[] = max(next_node[], Int(gmsh.model.mesh.getMaxNodeTag()))
    meshed = Set{Int32}()
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
            push!(meshed, face)
            for j in first:last, (a, b, c) in sector_triangles(section, j)
                append!(get!(triangles, face, Int[]),
                        (face_node!(face, a, i), face_node!(face, b, i), face_node!(face, c, i)))
            end
        end
        for j in first:last
            face = matched[(:lateral, (j, 0))]
            push!(meshed, face)
            a = section_node(section, K, j)
            b = section_node(section, K, j + 1)
            for i in 0:(L - 1)
                apex = apex_node!(state, next_node, j, i)
                push!(get!(face_nodes, face, Int[]), apex)
                a0, b0, b1, a1 = state.tags[(a, i)], state.tags[(b, i)], state.tags[(b, i + 1)],
                                 state.tags[(a, i + 1)]
                append!(get!(triangles, face, Int[]),
                        (a0, b0, apex, b0, b1, apex, b1, a1, apex, a1, a0, apex))
            end
        end
        for j in (first, last + 1)
            face = matched[(:radial, (j, 0))]
            face in meshed && continue
            push!(meshed, face)
            for k in 1:K, i in 0:(L - 1)
                a = section_node(section, k - 1, j)
                b = section_node(section, k, j)
                append!(get!(quads, face, Int[]),
                        (face_node!(face, a, i), face_node!(face, b, i), face_node!(face, b, i + 1),
                         face_node!(face, a, i + 1)))
            end
        end
    end
    Set(all_faces) == meshed || error("tube faces and matched faces differ")
    for (face, nodes) in face_nodes
        gmsh.model.mesh.addNodes(2, face, nodes, reduce(vcat, state.coordinates[t] for t in nodes))
    end
    for (face, connectivity) in triangles
        gmsh.model.mesh.addElementsByType(face, 2, Int[], connectivity)
    end
    for (face, connectivity) in quads
        gmsh.model.mesh.addElementsByType(face, 3, Int[], connectivity)
    end
    return Dict{String, Any}(
        "Faces" => length(all_faces), "PyramidApexes" => length(state.apex),
        "Triangles" => sum(length(v) ÷ 3 for v in values(triangles); init=0),
        "Quads" => sum(length(v) ÷ 4 for v in values(quads); init=0))
end

# Remove the OCC tube volumes (non-recursively: faces, curves, points and their
# meshes stay) so that the remaining volumes can be meshed with pyramids.
function remove_tube_volumes!(states::Vector{TubeMesh})
    gmsh.model.removeEntities([(3, volume) for state in states for (volume, _, _) in state.volumes],
                              false)
end

# Phase 3 (after generate(3)): discrete volumes with the interior nodes and the
# prisms. Returns the discrete volume tags per material and the census.
function finalize_tube_volumes!(states::Vector{TubeMesh})
    volumes = Dict{Int, Vector{Int32}}()
    census = Dict{String, Any}[]
    next_node = Ref(Int(gmsh.model.mesh.getMaxNodeTag()))
    for state in states
        section = state.section
        K = ring_count(section)
        L = state.tube.layers
        for (volume, group, _) in state.volumes
            first, last, material = group
            rays = first:(last + 1)
            discrete = gmsh.model.addDiscreteEntity(3, -1, state.faces[volume])
            state.discrete[volume] = discrete
            interior = Int[]
            for i in 1:(L - 1), local_index in vcat([section_node(section, 0, 0)],
                                                     [section_node(section, k, j) for k in 1:K for j in rays])
                haskey(state.tags, (local_index, i)) && continue
                push!(interior, tube_node!(state, next_node, local_index, i))
            end
            isempty(interior) ||
                gmsh.model.mesh.addNodes(3, discrete, interior,
                                         reduce(vcat, state.coordinates[t] for t in interior))
            prisms = Int[]
            for i in 0:(L - 1), j in first:last, (a, b, c) in sector_triangles(section, j)
                append!(prisms, (state.tags[(a, i)], state.tags[(b, i)], state.tags[(c, i)],
                                 state.tags[(a, i + 1)], state.tags[(b, i + 1)], state.tags[(c, i + 1)]))
            end
            gmsh.model.mesh.addElementsByType(discrete, 6, Int[], prisms)
            pyramids = Int[]
            for i in 0:(L - 1), j in first:last
                a = section_node(section, K, j)
                b = section_node(section, K, j + 1)
                append!(pyramids, (state.tags[(a, i)], state.tags[(b, i)], state.tags[(b, i + 1)],
                                   state.tags[(a, i + 1)], state.apex[(j, i)]))
            end
            gmsh.model.mesh.addElementsByType(discrete, 7, Int[], pyramids)
            push!(get!(volumes, material, Int32[]), discrete)
            push!(census, Dict{String, Any}(
                "Volume" => Int(volume), "DiscreteVolume" => Int(discrete), "Material" => material,
                "Sectors" => last - first + 1, "Layers" => L, "Prisms" => length(prisms) ÷ 6,
                "Pyramids" => length(pyramids) ÷ 5, "InteriorNodes" => length(interior)))
        end
    end
    return volumes, census
end

# Build every tube volume in the OCC model (before the fragment) and return the
# records needed afterwards.
struct TubeRecord
    tube::EdgeTube
    section::TubeSection
    group::Tuple{Int, Int, Int}
    tool::Tuple{Int32, Int32}
end

# Fragment map lookup: the single volume descendant of a tube tool.
function tube_volume_after_fragment(record::TubeRecord, tool_index, fragment_map, object_count)
    descendants = [tag for (dim, tag) in fragment_map[object_count + tool_index] if dim == 3]
    length(descendants) == 1 ||
        error("tube tool $(record.tool) has $(length(descendants)) volume descendants")
    return descendants[1]
end

# Anisotropy of the tube prisms: the largest longitudinal spacing over the
# smallest cross-section edge (the innermost arc), and the ring sizes.
function tube_anisotropy(tube::EdgeTube, section::TubeSection, pyramid_height)
    spacing = tube_spacing(tube)
    inner_arc = section.ring_radii[1] * deg2rad(minimum(diff(section.angles)))
    outer_arc = 2.0 * tube_radius(section) * sin(0.5 * deg2rad(maximum(diff(section.angles))))
    return Dict{String, Any}(
        "LongitudinalSpacing" => spacing, "RingSizes" => ring_sizes(section),
        "OuterArc" => outer_arc, "PyramidHeight" => pyramid_height,
        "SectorDegrees" => diff(section.angles),
        "RingRadii" => copy(section.ring_radii), "InnermostArc" => inner_arc,
        "MaximumEdgeAspect" => spacing / min(inner_arc, ring_sizes(section)[1]))
end
