# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabricated (finite-metal, overetched) reference mesher for validation WINDOWS given as a
# POLYGON SET (validation plan section (c), USER decision 184): the single transmon's
# shared-plan tensor-sweep generator generalised from the DeviceLayout geometry to
#
#   * an arbitrary number of conductors per plane, each labelled; the label `ground` is shared
#     by every plane and every bump, every other label is a terminal with its own metal-air /
#     metal-substrate attribute pair;
#   * one or two metal levels (planes): a plane facing `up` carries its metal on top of a
#     substrate below its surface z (the transmon chip), a plane facing `down` carries its
#     metal below a substrate above its surface (the flip-chip L2). Both planes' metal edges
#     constrain ONE plan mesh, swept through explicit z levels with each plane's metal band,
#     overetch band and the graded gap between the levels;
#   * bumps: plan footprints extruded as excluded metal columns from the L1 metal top to the
#     L2 metal bottom, joining the two grounds (their footprint must lie on metal of BOTH
#     planes);
#   * a truncated box: the plan rectangle, each plane's substrate thickness and optional
#     vacuum beyond the outermost substrate backsides.
#
# Process geometry as the transmon reference: metal thickness 0.1 um with vertical sidewalls,
# 50-nm overetch recess of the exposed substrate, metal volume excluded from the domain, the
# complete conforming substrate-air interface incl. the overetch step, separate metal-air /
# metal-substrate shells per conductor; SA / MS / MA participations are surface integrals on
# those shells (no meshed 2-nm layers). Resolution: edge-normal boundary layers of first
# height r (geometric growth 2, band ~1.55 um) on every metal edge of every plane and bump,
# tangent spacing t along the edges, r-resolved metal and overetch bands in z.
#
# Method (unchanged from generate_shared_plan_mesh.jl): the plan rectangle is fragmented by
# every polygon (all planes) and bump footprint; each resulting partition is classified by
# the fragment map (conductor per plane, bump); each partition is copied and meshed alone
# with a one-sided Gmsh BoundaryLayer field on its metal-edge curves (transfinite tangent
# spacing), the copies' coincident nodes are welded by coordinate, the plan is swept through
# the z levels, prisms are split conformally into tetrahedra, the metal is omitted, unused
# nodes compacted and the physical attributes written directly to an ASCII MSH2 file with a
# JSON manifest (counts, per-attribute areas / volumes, first-layer heights, z levels).
#
# Deliberate deviation from the recorded transmon generator (supervisor decision 188): the
# BoundaryLayer Thickness is the geometric sum r (2^n - 1) times (1 + 1e-6) by default. The
# recorded generator passes the exact sum, so floating-point rounding decides whether a
# column gets n or n - 1 rows (~8 % of the columns short at r10, ~20 % at r50, also in the
# recorded meshes), and on Linux that mix fails Gmsh's edge recovery at r10. With the margin
# every column has exactly n rows on every platform. `exact_band_thickness=true`
# (`--exact-band-thickness`) reproduces the recorded meshes.
#
# Polygon-set JSON (micrometres; `Version` 1):
#   {"Name": ..., "Box": {"X": [x0, x1], "Y": [y0, y1]},
#    "Process": {"MetalThickness": 0.1, "Overetch": 0.05},
#    "Planes": [{"Name": "L1", "SurfaceZ": 0.0, "Facing": "up", "SubstrateThickness": 525.0,
#                "Polygons": [{"Conductor": "ground", "Outer": [[x, y], ...],
#                              "Holes": [[[x, y], ...], ...]}, ...]}, ...],
#    "Bumps": [{"Conductor": "ground", "Footprint": [[x, y], ...]}],
#    "Vacuum": {"Below": 0.0, "Above": 0.0},
#    "Terminals": ["island", ...] (optional attribute order; default sorted labels)}
# Polygons of one plane must not overlap; vertices may lie on the box wall (a ground plane
# reaching the wall); wall edges are not metal edges (no boundary layer).
#
# Attributes: 3D 1 substrate, 2 vacuum; 2D 3 exterior_boundary, 4 ground_air,
# 5 ground_substrate, 6 substrate_air, 7 <terminal 1>_air, 8 <terminal 1>_substrate,
# 9 substrate_backside (a substrate face against vacuum beyond the box: SA diagnostic),
# 10 / 11 <terminal 2>_air / _substrate, ... (the transmon reference's table for one
# terminal). The manifest records the table.

module PolygonWindowMesh

import Gmsh: gmsh
using JSON
using Printf
using SHA
using Statistics

export read_polygon_set, mesh_polygon_window, plan_partitions

const Point2 = NTuple{2, Float64}

struct Polygon
    conductor::String
    outer::Vector{Point2}
    holes::Vector{Vector{Point2}}
end

struct Plane
    name::String
    surface_z::Float64
    facing::Int # +1: metal above the surface (substrate below); -1: metal below (flip-chip)
    substrate_thickness::Float64
    polygons::Vector{Polygon}
end

struct Bump
    conductor::String
    footprint::Vector{Point2}
end

struct PolygonSet
    name::String
    box::NTuple{4, Float64} # xmin, xmax, ymin, ymax
    metal_thickness::Float64
    overetch::Float64
    planes::Vector{Plane}
    bumps::Vector{Bump}
    vacuum_below::Float64
    vacuum_above::Float64
    terminals::Vector{String}
end

const GROUND = "ground"

point_list(values) = [(Float64(p[1]), Float64(p[2])) for p in values]

function polygon_signed_area(points::Vector{Point2})
    n = length(points)
    return 0.5 * sum(
        points[i][1] * points[mod1(i + 1, n)][2] - points[mod1(i + 1, n)][1] * points[i][2]
        for i = 1:n
    )
end

function read_polygon_set(data::AbstractDict)
    get(data, "Version", 1) == 1 || error("Unsupported polygon-set version")
    box = data["Box"]
    x, y = box["X"], box["Y"]
    x[1] < x[2] && y[1] < y[2] || error("Degenerate box")
    process = get(data, "Process", Dict())
    planes = Plane[]
    labels = Set{String}()
    for plane in data["Planes"]
        facing = plane["Facing"]
        facing in ("up", "down") || error("Plane facing must be up or down")
        polygons = Polygon[]
        for polygon in plane["Polygons"]
            outer = point_list(polygon["Outer"])
            length(outer) >= 3 || error("Polygon outer loop needs 3 vertices")
            holes = [point_list(h) for h in get(polygon, "Holes", [])]
            all(h -> length(h) >= 3, holes) || error("Polygon hole needs 3 vertices")
            label = String(polygon["Conductor"])
            label == GROUND || push!(labels, label)
            push!(polygons, Polygon(label, outer, holes))
        end
        isempty(polygons) && error("Plane $(plane["Name"]) has no polygons")
        push!(
            planes,
            Plane(
                String(plane["Name"]),
                Float64(plane["SurfaceZ"]),
                facing == "up" ? 1 : -1,
                Float64(plane["SubstrateThickness"]),
                polygons
            )
        )
    end
    1 <= length(planes) <= 2 || error("One or two planes are supported")
    sort!(planes; by=p -> p.surface_z)
    if length(planes) == 2
        planes[1].facing == 1 && planes[2].facing == -1 || error(
            "Two planes must be a lower plane facing up and an upper plane facing down"
        )
        planes[2].surface_z - planes[1].surface_z >
        2 * Float64(get(process, "MetalThickness", 0.1)) ||
            error("Planes too close for two metal bands")
    end
    bumps = Bump[]
    for bump in get(data, "Bumps", [])
        footprint = point_list(bump["Footprint"])
        length(footprint) >= 3 || error("Bump footprint needs 3 vertices")
        label = String(get(bump, "Conductor", GROUND))
        label == GROUND || push!(labels, label)
        push!(bumps, Bump(label, footprint))
    end
    isempty(bumps) || length(planes) == 2 || error("Bumps need two planes")
    vacuum = get(data, "Vacuum", Dict())
    terminals = String.(get(data, "Terminals", sort!(collect(labels))))
    Set(terminals) == labels || error("Terminals must list every non-ground conductor once")
    length(unique(terminals)) == length(terminals) || error("Duplicate terminal label")
    return PolygonSet(
        String(get(data, "Name", "window")),
        (Float64(x[1]), Float64(x[2]), Float64(y[1]), Float64(y[2])),
        Float64(get(process, "MetalThickness", 0.1)),
        Float64(get(process, "Overetch", 0.05)),
        planes,
        bumps,
        Float64(get(vacuum, "Below", 0.0)),
        Float64(get(vacuum, "Above", 0.0)),
        terminals
    )
end

read_polygon_set(path::AbstractString) = read_polygon_set(JSON.parsefile(path))

# Attribute table: 4 / 5 ground, 7 / 8 the first terminal, 10 / 11 the second, ...
function attribute_table(spec::PolygonSet)
    table = Dict{String, Tuple{Int, Int}}(GROUND => (4, 5))
    for (index, label) in enumerate(spec.terminals)
        table[label] = index == 1 ? (7, 8) : (10 + 2 * (index - 2), 11 + 2 * (index - 2))
    end
    return table
end

function physical_names(spec::PolygonSet)
    table = attribute_table(spec)
    names = [
        (3, 1, "substrate"),
        (3, 2, "vacuum"),
        (2, 3, "exterior_boundary"),
        (2, 4, "ground_air"),
        (2, 5, "ground_substrate"),
        (2, 6, "substrate_air"),
        (2, 9, "substrate_backside")
    ]
    for label in spec.terminals
        air, substrate = table[label]
        push!(names, (2, air, label * "_air"), (2, substrate, label * "_substrate"))
    end
    return sort!(names; by=entry -> (entry[1], entry[2]))
end

# ---------------------------------------------------------------------------------------------
# Plan geometry: the rectangle fragmented by every polygon of every plane and every bump.

# Gmsh's OCC plane surface expects every hole wire with the SAME orientation as the outer
# wire (it reverses the holes itself), so every loop is passed counterclockwise.
function add_polygon_surface(occ, outer::Vector{Point2}, holes::Vector{Vector{Point2}})
    function loop(points)
        polygon_signed_area(points) < 0.0 && (points = reverse(points))
        tags = [occ.add_point(p[1], p[2], 0.0) for p in points]
        lines = [
            occ.add_line(tags[i], tags[mod1(i + 1, length(tags))]) for i in eachindex(tags)
        ]
        return occ.add_curve_loop(lines)
    end
    return occ.add_plane_surface(vcat(loop(outer), [loop(h) for h in holes]))
end

# Partition classification: conductor label per plane ("" = none) and bump index (0 = none).
struct PartitionClass
    conductors::Vector{String}
    bump::Int
end
Base.:(==)(a::PartitionClass, b::PartitionClass) =
    a.conductors == b.conductors && a.bump == b.bump
Base.hash(c::PartitionClass, h::UInt) = hash((c.conductors, c.bump), h)

function on_box_wall(box, bbox; tolerance=1.0e-5)
    return (abs(bbox[1] - box[1]) <= tolerance && abs(bbox[4] - box[1]) <= tolerance) ||
           (abs(bbox[1] - box[2]) <= tolerance && abs(bbox[4] - box[2]) <= tolerance) ||
           (abs(bbox[2] - box[3]) <= tolerance && abs(bbox[5] - box[3]) <= tolerance) ||
           (abs(bbox[2] - box[4]) <= tolerance && abs(bbox[5] - box[4]) <= tolerance)
end

"""
    plan_partitions(spec) -> (surfaces, classes, metal_edge_curves)

Build the fragmented plan in the current Gmsh model. Returns the partition surfaces
(dimtags), their classes, and the set of curve tags that are metal edges (a curve whose two
adjacent partitions differ in conductor on some plane or in bump membership).
"""
function plan_partitions(spec::PolygonSet)
    occ = gmsh.model.occ
    box = spec.box
    rectangle = occ.add_rectangle(box[1], box[3], 0.0, box[2] - box[1], box[4] - box[3])
    tools = Tuple{Int32, Int32}[]
    tool_meaning = Tuple{Int, Int, String}[] # (plane index or 0, bump index or 0, label)
    for (plane_index, plane) in enumerate(spec.planes), polygon in plane.polygons
        push!(tools, (2, add_polygon_surface(occ, polygon.outer, polygon.holes)))
        push!(tool_meaning, (plane_index, 0, polygon.conductor))
    end
    for (bump_index, bump) in enumerate(spec.bumps)
        push!(tools, (2, add_polygon_surface(occ, bump.footprint, Vector{Point2}[])))
        push!(tool_meaning, (0, bump_index, bump.conductor))
    end
    fragments, fragment_map = occ.fragment([(2, rectangle)], tools)
    occ.synchronize()
    all(entity -> entity[1] == 2, fragments) ||
        error("Plan fragmentation produced non-surface entities")
    rectangle_fragments = Set(fragment_map[1])
    surfaces = sort!(collect(rectangle_fragments))
    n_planes = length(spec.planes)
    classes = Dict(surface => PartitionClass(fill("", n_planes), 0) for surface in surfaces)
    for (tool_index, (plane_index, bump_index, label)) in enumerate(tool_meaning)
        for fragment in fragment_map[tool_index + 1]
            haskey(classes, fragment) || error(
                "A polygon of $(label) leaves the plan rectangle (fragment $fragment)"
            )
            class = classes[fragment]
            if plane_index > 0
                isempty(class.conductors[plane_index]) || error(
                    "Polygons overlap on plane $(spec.planes[plane_index].name) " *
                    "($(class.conductors[plane_index]) and $label)"
                )
                class.conductors[plane_index] = label
            else
                class.bump == 0 || error("Bump footprints overlap")
                classes[fragment] = PartitionClass(class.conductors, bump_index)
            end
        end
    end
    for surface in surfaces
        class = classes[surface]
        if class.bump > 0
            all(!isempty, class.conductors) ||
                error("Bump $(class.bump) footprint is not on metal of both planes")
        end
    end
    # Metal edges: curves between partitions of different class, strictly inside the box.
    metal_edge_curves = Set{Int32}()
    for (_, tag) in unique(gmsh.model.get_boundary(surfaces, false, false, false))
        tag = abs(tag)
        upward, _ = gmsh.model.get_adjacencies(1, tag)
        owners = [(2, Int32(s)) for s in upward if haskey(classes, (2, Int32(s)))]
        if length(owners) == 1
            on_box_wall(box, gmsh.model.get_bounding_box(1, tag)) ||
                error("Plan curve $tag bounds one partition but is not on the box wall")
        elseif length(owners) == 2
            classes[owners[1]] == classes[owners[2]] || push!(metal_edge_curves, tag)
        else
            error("Plan curve $tag bounds $(length(owners)) partitions")
        end
    end
    isempty(metal_edge_curves) && error("No metal edges found")
    return surfaces, classes, metal_edge_curves
end

# Curves are matched between the original and its copies by geometry.
function curve_signature(dimension, tag)
    entity = (dimension, abs(tag))
    bbox = gmsh.model.get_bounding_box(entity...)
    length_um = gmsh.model.occ.get_mass(entity...)
    bounds = gmsh.model.get_parametrization_bounds(entity...)
    mid = gmsh.model.get_value(entity..., [0.5 * (bounds[1][1] + bounds[2][1])])
    return Tuple(round.(collect((bbox..., length_um, mid[1], mid[2])); digits=6))
end

function triangle_area(a::Point2, b::Point2, c::Point2)
    return 0.5 * abs((b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]))
end

function metric_summary(values)
    return Dict(
        "minimum" => minimum(values),
        "p01" => quantile(values, 0.01),
        "median" => median(values),
        "p99" => quantile(values, 0.99),
        "maximum" => maximum(values)
    )
end

# ---------------------------------------------------------------------------------------------
# z levels

const SUBSTRATE_SIDE_OFFSETS_UM = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
const VACUUM_SIDE_OFFSETS_UM =
    [0.15, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0]

struct ZStack
    levels::Vector{Float64}
    z_bottom::Float64
    z_top::Float64
    backsides::Vector{Float64} # substrate faces against vacuum beyond the box
end

function z_levels(spec::PolygonSet, metal_layers::Int, trench_layers::Int)
    m, o = spec.metal_thickness, spec.overetch
    planes = spec.planes
    levels = Float64[]
    backsides = Float64[]
    metal_tops = Float64[]
    for plane in planes
        s, f = plane.surface_z, plane.facing
        # The r-resolved trench and metal bands, then the graded substrate side.
        append!(levels, range(s - f * o, s; length=trench_layers + 1))
        append!(levels, range(s, s + f * m; length=metal_layers + 1))
        push!(metal_tops, s + f * m)
        for offset in SUBSTRATE_SIDE_OFFSETS_UM
            offset < plane.substrate_thickness || break
            push!(levels, s - f * offset)
        end
        push!(levels, s - f * plane.substrate_thickness)
    end
    if length(planes) == 2
        # Gap between the metal tops: each plane grades toward the midpoint.
        midpoint = 0.5 * (metal_tops[1] + metal_tops[2])
        for plane in planes
            s, f = plane.surface_z, plane.facing
            for offset in VACUUM_SIDE_OFFSETS_UM
                z = s + f * offset
                f * (midpoint - z) > 0.0 || break
                push!(levels, z)
            end
        end
        push!(levels, midpoint)
        lower, upper = planes
        z_bottom = lower.surface_z - lower.substrate_thickness - spec.vacuum_below
        z_top = upper.surface_z + upper.substrate_thickness + spec.vacuum_above
        spec.vacuum_below > 0.0 &&
            push!(backsides, lower.surface_z - lower.substrate_thickness)
        spec.vacuum_above > 0.0 &&
            push!(backsides, upper.surface_z + upper.substrate_thickness)
    else
        plane = planes[1]
        s, f = plane.surface_z, plane.facing
        for offset in VACUUM_SIDE_OFFSETS_UM
            push!(levels, s + f * offset)
        end
        backside = s - f * plane.substrate_thickness
        if f == 1
            spec.vacuum_above > m || error("Vacuum.Above must exceed the metal thickness")
            z_bottom = backside - spec.vacuum_below
            z_top = s + spec.vacuum_above
            spec.vacuum_below > 0.0 && push!(backsides, backside)
        else
            spec.vacuum_below > m || error("Vacuum.Below must exceed the metal thickness")
            z_bottom = s - spec.vacuum_below
            z_top = backside + spec.vacuum_above
            spec.vacuum_above > 0.0 && push!(backsides, backside)
        end
    end
    push!(levels, z_bottom, z_top)
    filter!(z -> z_bottom - 1.0e-9 <= z <= z_top + 1.0e-9, levels)
    sort!(levels)
    merged = Float64[levels[1]]
    for z in levels[2:end]
        z - merged[end] > 1.0e-9 && push!(merged, z)
    end
    return ZStack(merged, z_bottom, z_top, backsides)
end

# Material of a prism at height z for a partition class: 1 substrate, 2 vacuum, 0 excluded.
function material(spec::PolygonSet, class::PartitionClass, z::Float64)
    m, o = spec.metal_thickness, spec.overetch
    for (k, plane) in enumerate(spec.planes)
        s, f = plane.surface_z, plane.facing
        d = f * (z - s) # signed distance from the surface toward the metal side
        metal = !isempty(class.conductors[k])
        if -plane.substrate_thickness < d < -o
            return Int8(1)
        elseif -o < d < 0.0
            return metal ? Int8(1) : Int8(2)
        elseif 0.0 < d < m
            metal && return Int8(0)
        end
    end
    if class.bump > 0
        lower, upper = spec.planes
        lower.surface_z + m < z < upper.surface_z - m && return Int8(0)
    end
    return Int8(2)
end

# ---------------------------------------------------------------------------------------------

# Relative margin on the BoundaryLayer Thickness (decision 188): every column gets exactly
# `radial_layers` rows instead of n or n - 1 by rounding of the exact geometric sum.
const BAND_THICKNESS_MARGIN = 1.0e-6

struct PlanMesh
    xy::Vector{Point2}
    triangles::Vector{NTuple{3, Int}}
    triangle_class::Vector{Int} # index into classes
    classes::Vector{PartitionClass}
end

function mesh_plan(
    spec::PolygonSet,
    radial_um::Float64,
    tangential_um::Float64;
    radial_growth::Float64=2.0,
    radial_band_um::Float64=1.55,
    exact_band_thickness::Bool=false,
    verbose::Bool=true
)
    surfaces, class_by_surface, metal_edge_curves = plan_partitions(spec)
    metal_signatures = Set(curve_signature(1, tag) for tag in metal_edge_curves)
    class_list = unique(class_by_surface[s] for s in surfaces)
    class_index = Dict(class => i for (i, class) in enumerate(class_list))

    # One-sided boundary layers need every source curve to bound a single surface: copy each
    # partition and mesh the copies one at a time (coincident copies never share a Delaunay
    # problem); the copies' interface nodes are welded by coordinate afterwards.
    copies = Tuple{Int32, Int32}[]
    copy_class = Int[]
    for surface in surfaces
        append!(copies, gmsh.model.occ.copy([surface]))
        push!(copy_class, class_index[class_by_surface[surface]])
    end
    gmsh.model.occ.synchronize()
    gmsh.model.set_visibility(gmsh.model.get_entities(), 0, true)
    gmsh.option.set_number("Mesh.MeshOnlyVisible", 1)
    copy_curves = unique(gmsh.model.get_boundary(copies, false, false, false))
    metal_copy_curves = Int32[]
    for (_, tag) in copy_curves
        if curve_signature(1, tag) in metal_signatures
            push!(metal_copy_curves, abs(tag))
            length_um = gmsh.model.occ.get_mass(1, abs(tag))
            segments = max(1, ceil(Int, length_um / tangential_um - 1.0e-6))
            gmsh.model.mesh.set_transfinite_curve(abs(tag), segments + 1)
        end
    end
    metal_copy_set = Set(metal_copy_curves)
    length(metal_copy_curves) == 2 * length(metal_edge_curves) || error(
        "Metal-edge curve copies: $(length(metal_copy_curves)) for " *
        "$(length(metal_edge_curves)) curves"
    )

    radial_layers = max(1, round(Int, log2(1.0 + radial_band_um / radial_um)))
    radial_thickness =
        radial_um * (radial_growth^radial_layers - 1.0) / (radial_growth - 1.0) *
        (exact_band_thickness ? 1.0 : 1.0 + BAND_THICKNESS_MARGIN)
    for option in (
        "General.NumThreads",
        "Mesh.MaxNumThreads1D",
        "Mesh.MaxNumThreads2D",
        "Mesh.MaxNumThreads3D"
    )
        gmsh.option.set_number(option, 1)
    end
    gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.set_number("Mesh.MeshSizeMin", radial_um)
    gmsh.option.set_number("Mesh.MeshSizeMax", 30.0)
    gmsh.option.set_number("Mesh.BoundaryLayerFanElements", 7)
    gmsh.option.set_number("Mesh.Algorithm", 6)
    verbose && println(
        "Plan partitions: ",
        length(copies),
        ", metal-edge curves: ",
        length(metal_edge_curves),
        ", boundary layer first=",
        radial_um,
        " um, layers=",
        radial_layers,
        ", thickness=",
        radial_thickness,
        " um (",
        exact_band_thickness ? "exact geometric sum" : "geometric sum x (1 + 1e-6)",
        ")"
    )

    raw_triangles = Tuple{NTuple{3, Point2}, Int}[]
    for (copy_index, surface) in enumerate(copies)
        gmsh.model.set_visibility(copies, 0, true)
        gmsh.model.set_visibility([surface], 1, true)
        sources = [
            Float64(abs(tag)) for
            (_, tag) in gmsh.model.get_boundary([surface], false, false, false) if
            abs(tag) in metal_copy_set
        ]
        field = 0
        if !isempty(sources)
            field = gmsh.model.mesh.field.add("BoundaryLayer")
            gmsh.model.mesh.field.set_numbers(field, "CurvesList", sources)
            gmsh.model.mesh.field.set_number(field, "Size", radial_um)
            gmsh.model.mesh.field.set_number(field, "Ratio", radial_growth)
            gmsh.model.mesh.field.set_number(field, "Thickness", radial_thickness)
            gmsh.model.mesh.field.set_number(field, "Quads", 1)
            gmsh.model.mesh.field.set_number(field, "IntersectMetrics", 1)
            gmsh.model.mesh.field.set_as_boundary_layer(field)
        end
        gmsh.model.mesh.generate(2)
        node_tags, node_coordinates, _ = gmsh.model.mesh.get_nodes()
        coordinate_by_tag = Dict{UInt64, Point2}()
        for (index, tag) in enumerate(node_tags)
            abs(node_coordinates[3index]) <= 1.0e-6 || error("Plan node is not on z=0")
            coordinate_by_tag[tag] =
                (node_coordinates[3index - 2], node_coordinates[3index - 1])
        end
        area = 0.0
        element_types, _, nodes_by_type = gmsh.model.mesh.get_elements(2, surface[2])
        for (element_type, element_nodes) in zip(element_types, nodes_by_type)
            _, dimension, _, node_count, _, primary =
                gmsh.model.mesh.get_element_properties(element_type)
            dimension == 2 || continue
            primary in (3, 4) || error("Expected triangular or quadrilateral plan elements")
            for offset = 0:node_count:(length(element_nodes) - node_count)
                corners = [coordinate_by_tag[element_nodes[offset + i]] for i = 1:primary]
                split =
                    primary == 3 ? ((corners[1], corners[2], corners[3]),) :
                    (
                        (corners[1], corners[2], corners[3]),
                        (corners[1], corners[3], corners[4])
                    )
                for triangle in split
                    area += triangle_area(triangle...)
                    push!(raw_triangles, (triangle, copy_class[copy_index]))
                end
            end
        end
        occ_area = gmsh.model.occ.get_mass(surface...)
        abs(area - occ_area) <= 1.0e-6 * max(1.0, occ_area) ||
            error("Partition $(surface[2]) mesh area $area differs from OCC area $occ_area")
        verbose && println(
            "  partition ",
            copy_index,
            "/",
            length(copies),
            " surface ",
            surface[2],
            " class ",
            class_list[copy_class[copy_index]].conductors,
            " bump ",
            class_list[copy_class[copy_index]].bump,
            " area ",
            area,
            " sources ",
            length(sources)
        )
        field == 0 || gmsh.model.mesh.field.remove(field)
        gmsh.model.mesh.clear(copies)
    end
    isempty(raw_triangles) && error("No plan triangles extracted")

    merge_tolerance = 1.0e-6
    coordinate_index = Dict{NTuple{2, Int64}, Int}()
    xy = Point2[]
    function node_index(p::Point2)
        key = (round(Int64, p[1] / merge_tolerance), round(Int64, p[2] / merge_tolerance))
        index = get(coordinate_index, key, 0)
        if index == 0
            push!(xy, p)
            index = length(xy)
            coordinate_index[key] = index
        else
            q = xy[index]
            hypot(p[1] - q[1], p[2] - q[2]) <= 2merge_tolerance ||
                error("Plan-node merge tolerance collision")
        end
        return index
    end
    triangles = NTuple{3, Int}[]
    triangle_class = Int[]
    for (triangle, class) in raw_triangles
        push!(
            triangles,
            (node_index(triangle[1]), node_index(triangle[2]), node_index(triangle[3]))
        )
        push!(triangle_class, class)
    end
    return PlanMesh(xy, triangles, triangle_class, class_list),
    radial_layers,
    radial_thickness
end

# Plan edge -> owning triangles; metal perimeter edges per plane / bump with their conductor.
struct PlanTopology
    edge_incidence::Dict{Tuple{Int, Int}, Vector{Int}}
    open_edges::Int
    plane_edge_conductor::Vector{Dict{Tuple{Int, Int}, String}} # per plane
    bump_edge_conductor::Dict{Tuple{Int, Int}, String}
end

function plan_topology(spec::PolygonSet, plan::PlanMesh)
    edge_incidence = Dict{Tuple{Int, Int}, Vector{Int}}()
    for (index, t) in enumerate(plan.triangles)
        for (a, b) in ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
            push!(get!(edge_incidence, a < b ? (a, b) : (b, a), Int[]), index)
        end
    end
    any(owners -> length(owners) > 2, values(edge_incidence)) &&
        error("Plan mesh has nonmanifold edges")
    box = spec.box
    function on_wall(edge)
        a, b = plan.xy[edge[1]], plan.xy[edge[2]]
        tolerance = 1.0e-5
        return (abs(a[1] - box[1]) <= tolerance && abs(b[1] - box[1]) <= tolerance) ||
               (abs(a[1] - box[2]) <= tolerance && abs(b[1] - box[2]) <= tolerance) ||
               (abs(a[2] - box[3]) <= tolerance && abs(b[2] - box[3]) <= tolerance) ||
               (abs(a[2] - box[4]) <= tolerance && abs(b[2] - box[4]) <= tolerance)
    end
    open_edges = [edge for (edge, owners) in edge_incidence if length(owners) == 1]
    unexpected = count(!on_wall, open_edges)
    unexpected == 0 || error("Plan mesh has $unexpected open edges off the box wall")
    plane_edge_conductor = [Dict{Tuple{Int, Int}, String}() for _ in spec.planes]
    bump_edge_conductor = Dict{Tuple{Int, Int}, String}()
    for (edge, owners) in edge_incidence
        length(owners) == 2 || continue
        c1 = plan.classes[plan.triangle_class[owners[1]]]
        c2 = plan.classes[plan.triangle_class[owners[2]]]
        for k in eachindex(spec.planes)
            a, b = c1.conductors[k], c2.conductors[k]
            a == b && continue
            (isempty(a) || isempty(b)) ||
                error("Conductors $a and $b touch on plane $(spec.planes[k].name)")
            plane_edge_conductor[k][edge] = isempty(a) ? b : a
        end
        if c1.bump != c2.bump
            (c1.bump == 0 || c2.bump == 0) || error("Bumps touch")
            bump_edge_conductor[edge] = spec.bumps[max(c1.bump, c2.bump)].conductor
        end
    end
    return PlanTopology(
        edge_incidence,
        length(open_edges),
        plane_edge_conductor,
        bump_edge_conductor
    )
end

function edge_metrics(plan::PlanMesh, topology::PlanTopology, radial_um, tangential_um)
    metal_edges = Set{Tuple{Int, Int}}()
    for edges in topology.plane_edge_conductor
        union!(metal_edges, keys(edges))
    end
    union!(metal_edges, keys(topology.bump_edge_conductor))
    tangents = Float64[]
    heights = Float64[]
    outliers = 0
    for edge in metal_edges
        a, b = plan.xy[edge[1]], plan.xy[edge[2]]
        tangent = hypot(b[1] - a[1], b[2] - a[2])
        push!(tangents, tangent)
        for index in topology.edge_incidence[edge]
            t = plan.triangles[index]
            third = only(setdiff((t[1], t[2], t[3]), edge))
            height = 2 * triangle_area(a, b, plan.xy[third]) / tangent
            push!(heights, height)
            height > 1.5radial_um && (outliers += 1)
        end
    end
    tangent_summary = metric_summary(tangents)
    height_summary = metric_summary(heights)
    tangent_summary["maximum"] <= 1.01tangential_um ||
        error("Metal-edge tangent target was not enforced: $tangent_summary")
    height_summary["maximum"] <= 1.1radial_um ||
        error("Metal-edge first-layer normal target was not enforced: $height_summary")
    return length(metal_edges), tangent_summary, height_summary, outliers
end

# Connected components of every conductor label per plane (a diagnostic: an island split in
# two pieces or a ground split by the window is reported, not refused).
function conductor_components(spec::PolygonSet, plan::PlanMesh, topology::PlanTopology)
    result = Dict{String, Any}()
    for (k, plane) in enumerate(spec.planes)
        parent = collect(1:length(plan.triangles))
        function root(i)
            while parent[i] != i
                parent[i] = parent[parent[i]]
                i = parent[i]
            end
            return i
        end
        for owners in values(topology.edge_incidence)
            length(owners) == 2 || continue
            a = plan.classes[plan.triangle_class[owners[1]]].conductors[k]
            b = plan.classes[plan.triangle_class[owners[2]]].conductors[k]
            (a == b && !isempty(a)) && (parent[root(owners[1])] = root(owners[2]))
        end
        areas = Dict{String, Dict{Int, Float64}}()
        for (index, t) in enumerate(plan.triangles)
            label = plan.classes[plan.triangle_class[index]].conductors[k]
            isempty(label) && continue
            component = get!(areas, label, Dict{Int, Float64}())
            component[root(index)] =
                get(component, root(index), 0.0) +
                triangle_area(plan.xy[t[1]], plan.xy[t[2]], plan.xy[t[3]])
        end
        result[plane.name] = Dict(
            label => Dict(
                "components" => length(components),
                "area_um2" => sum(values(components)),
                "component_areas_um2" => sort!(collect(values(components)); rev=true)
            ) for (label, components) in areas
        )
    end
    return result
end

# ---------------------------------------------------------------------------------------------

signed_volume(p1, p2, p3, p4) =
    (
        (p2[1] - p1[1]) *
        ((p3[2] - p1[2]) * (p4[3] - p1[3]) - (p3[3] - p1[3]) * (p4[2] - p1[2])) -
        (p2[2] - p1[2]) *
        ((p3[1] - p1[1]) * (p4[3] - p1[3]) - (p3[3] - p1[3]) * (p4[1] - p1[1])) +
        (p2[3] - p1[3]) *
        ((p3[1] - p1[1]) * (p4[2] - p1[2]) - (p3[2] - p1[2]) * (p4[1] - p1[1]))
    ) / 6

function face_area(p1, p2, p3)
    ux, uy, uz = p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]
    vx, vy, vz = p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3]
    return 0.5 * hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
end

"""
    mesh_polygon_window(spec, radial_um, tangential_um, output; verbose=true,
                        plan_only=false, exact_band_thickness=false) -> manifest

Generate the fabricated reference mesh of a polygon set and write `output` (ASCII MSH2) with
its JSON manifest next to it. `exact_band_thickness=true` passes the exact geometric sum as
the boundary-layer Thickness (the recorded transmon generator's formula; see the header).
"""
function mesh_polygon_window(
    spec::PolygonSet,
    radial_um::Float64,
    tangential_um::Float64,
    output::AbstractString;
    verbose::Bool=true,
    plan_only::Bool=false,
    exact_band_thickness::Bool=false
)
    output = abspath(output)
    metal_layers = max(2, ceil(Int, spec.metal_thickness / radial_um - 1.0e-9))
    trench_layers = max(1, ceil(Int, spec.overetch / radial_um - 1.0e-9))
    gmsh.initialize()
    plan, radial_layers, radial_thickness = try
        gmsh.option.set_number("General.Verbosity", 2)
        gmsh.model.add(spec.name)
        mesh_plan(
            spec,
            radial_um,
            tangential_um;
            exact_band_thickness=exact_band_thickness,
            verbose=verbose
        )
    finally
        gmsh.finalize()
    end
    topology = plan_topology(spec, plan)
    box = spec.box
    plan_area = sum(
        triangle_area(plan.xy[t[1]], plan.xy[t[2]], plan.xy[t[3]]) for t in plan.triangles
    )
    expected_area = (box[2] - box[1]) * (box[4] - box[3])
    abs(plan_area - expected_area) / expected_area < 1.0e-8 ||
        error("Plan area mismatch: $plan_area vs $expected_area")
    metal_edge_count, tangent_summary, height_summary, height_outliers =
        edge_metrics(plan, topology, radial_um, tangential_um)
    components = conductor_components(spec, plan, topology)
    verbose && println(
        "Plan nodes: ",
        length(plan.xy),
        ", triangles: ",
        length(plan.triangles),
        ", metal perimeter edges: ",
        metal_edge_count,
        "\nTangent (um): ",
        tangent_summary,
        "\nFirst-layer normal height (um): ",
        height_summary,
        ", outliers > 1.5 r: ",
        height_outliers,
        "\nConductor components: ",
        components
    )
    stack = z_levels(spec, metal_layers, trench_layers)
    zs = stack.levels
    manifest = Dict{String, Any}(
        "mesh" => output,
        "polygon_set" => spec.name,
        "radial_target_um" => radial_um,
        "tangential_target_um" => tangential_um,
        "radial_growth" => 2.0,
        "radial_layers" => radial_layers,
        "radial_band_thickness_um" => radial_thickness,
        "radial_band_thickness_mode" =>
            exact_band_thickness ? "exact_geometric_sum" : "geometric_sum_x_1p000001",
        "metal_thickness_um" => spec.metal_thickness,
        "overetch_um" => spec.overetch,
        "metal_layers" => metal_layers,
        "trench_layers" => trench_layers,
        "planes" => [
            Dict(
                "name" => p.name,
                "surface_z_um" => p.surface_z,
                "facing" => p.facing == 1 ? "up" : "down",
                "substrate_thickness_um" => p.substrate_thickness
            ) for p in spec.planes
        ],
        "bumps" => length(spec.bumps),
        "plan_nodes" => length(plan.xy),
        "plan_triangles" => length(plan.triangles),
        "plan_partition_classes" => length(plan.classes),
        "plan_open_outer_edges" => topology.open_edges,
        "plan_perimeter_edges" => metal_edge_count,
        "perimeter_tangent_length_um" => tangent_summary,
        "first_layer_normal_height_um" => height_summary,
        "first_layer_normal_outliers_above_1p5x_target" => height_outliers,
        "conductor_components" => components,
        "z_bounds_um" => Dict("bottom" => stack.z_bottom, "top" => stack.z_top),
        "substrate_backsides_z_um" => stack.backsides,
        "z_levels" => zs,
        "attributes" => Dict(
            name => Dict("dimension" => d, "attribute" => a) for
            (d, a, name) in physical_names(spec)
        )
    )
    plan_only && return manifest

    # Sweep: every plan triangle x every z interval is a prism of one material (or excluded).
    n_plan = length(plan.xy)
    n_levels = length(zs)
    material_table = [
        material(spec, class, 0.5 * (zs[i] + zs[i + 1])) for
        class in plan.classes, i = 1:(n_levels - 1)
    ]
    nodes = Vector{NTuple{3, Float64}}(undef, n_plan * n_levels)
    for level = 1:n_levels, index = 1:n_plan
        nodes[(level - 1) * n_plan + index] =
            (plan.xy[index][1], plan.xy[index][2], zs[level])
    end
    tetrahedra = NTuple{4, Int32}[]
    tetrahedron_attribute = Int8[]
    sizehint!(tetrahedra, length(plan.triangles) * (n_levels - 1) * 3)
    function push_prism!(a1, a2, a3, b1, b2, b3, attribute)
        # Conforming split: the diagonal choice on every quad side follows global node order.
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
    for (index, t) in enumerate(plan.triangles)
        class = plan.triangle_class[index]
        for interval = 1:(n_levels - 1)
            attribute = material_table[class, interval]
            attribute == 0 && continue
            lower = Int32((interval - 1) * n_plan)
            upper = Int32(interval * n_plan)
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
    for index in eachindex(tetrahedra)
        t = tetrahedra[index]
        volume = signed_volume(nodes[t[1]], nodes[t[2]], nodes[t[3]], nodes[t[4]])
        if volume < 0
            tetrahedra[index] = (t[1], t[3], t[2], t[4])
            volume = -volume
        end
        volume > 0 || error("Degenerate tetrahedron $index")
    end
    verbose && println("z levels: ", n_levels, ", tetrahedra: ", length(tetrahedra))

    # Faces: count of owning tetrahedra and the sum of their attributes.
    face_data = Dict{NTuple{3, Int32}, Tuple{Int8, Int8}}()
    sizehint!(face_data, 2 * length(tetrahedra) + length(plan.triangles) * 4)
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

    table = attribute_table(spec)
    plane_by_z = Dict{Float64, Tuple{Int, Symbol}}() # metal top / metal bottom z per plane
    for (k, plane) in enumerate(spec.planes)
        plane_by_z[plane.surface_z] = (k, :substrate)
        plane_by_z[plane.surface_z + plane.facing * spec.metal_thickness] = (k, :air)
    end
    level_of_z = Dict(z => i for (i, z) in enumerate(zs))
    plan_face_class = Dict{NTuple{3, Int}, Int}()
    for (index, t) in enumerate(plan.triangles)
        plan_face_class[Tuple(sort([t[1], t[2], t[3]]))] = plan.triangle_class[index]
    end
    base_index(node) = mod(Int(node) - 1, n_plan) + 1
    level_index(node) = (Int(node) - 1) ÷ n_plan + 1
    lower_metal_top =
        spec.planes[1].surface_z + spec.planes[1].facing * spec.metal_thickness
    upper_metal_bottom =
        length(spec.planes) == 2 ?
        spec.planes[2].surface_z + spec.planes[2].facing * spec.metal_thickness : NaN

    surface_elements = Tuple{Int, NTuple{3, Int32}}[]
    for (face, (count, attribute_sum)) in face_data
        if count == 1
            levels = (level_index(face[1]), level_index(face[2]), level_index(face[3]))
            coordinates = (nodes[face[1]], nodes[face[2]], nodes[face[3]])
            attribute = 0
            if levels[1] == levels[2] == levels[3]
                z = zs[levels[1]]
                if levels[1] == 1 || levels[1] == n_levels
                    attribute = 3
                else
                    haskey(plane_by_z, z) || error("Unclassified horizontal face at z=$z")
                    k, side = plane_by_z[z]
                    class_index = get(
                        plan_face_class,
                        Tuple(sort(unique(base_index.(collect(face))))),
                        0
                    )
                    class_index > 0 ||
                        error("Horizontal face is not a plan triangle at z=$z")
                    class = plan.classes[class_index]
                    label = class.conductors[k]
                    isempty(label) &&
                        error("Horizontal metal face without conductor at z=$z")
                    attribute = side == :air ? table[label][1] : table[label][2]
                end
            else
                plan_nodes = unique(base_index.(collect(face)))
                length(plan_nodes) == 2 || error("Non-vertical boundary face")
                edge = Tuple(sort(plan_nodes))
                if all(abs(p[1] - box[1]) <= 1.0e-6 for p in coordinates) ||
                   all(abs(p[1] - box[2]) <= 1.0e-6 for p in coordinates) ||
                   all(abs(p[2] - box[3]) <= 1.0e-6 for p in coordinates) ||
                   all(abs(p[2] - box[4]) <= 1.0e-6 for p in coordinates)
                    attribute = 3
                else
                    zmid = sum(p[3] for p in coordinates) / 3
                    label = ""
                    for (k, plane) in enumerate(spec.planes)
                        d = plane.facing * (zmid - plane.surface_z)
                        if 0.0 < d < spec.metal_thickness
                            label = get(topology.plane_edge_conductor[k], edge, "")
                            isempty(label) &&
                                error("Sidewall face off the metal perimeter at z=$zmid")
                        end
                    end
                    if isempty(label) && lower_metal_top < zmid < upper_metal_bottom
                        label = get(topology.bump_edge_conductor, edge, "")
                        isempty(label) &&
                            error("Bump sidewall face off a bump perimeter at z=$zmid")
                    end
                    isempty(label) && error("Unclassified vertical face at z=$zmid")
                    attribute = table[label][1]
                end
            end
            push!(surface_elements, (attribute, face))
        elseif count == 2 && attribute_sum == 3
            levels = (level_index(face[1]), level_index(face[2]), level_index(face[3]))
            z = levels[1] == levels[2] == levels[3] ? zs[levels[1]] : NaN
            attribute = any(b -> abs(z - b) <= 1.0e-9, stack.backsides) ? 9 : 6
            push!(surface_elements, (attribute, face))
        end
    end
    face_data = nothing
    GC.gc()

    # Compact the nodes inside the excluded metal.
    uncompacted = length(nodes)
    used = falses(length(nodes))
    for t in tetrahedra, node in t
        used[node] = true
    end
    remap = zeros(Int32, length(nodes))
    compacted = NTuple{3, Float64}[]
    sizehint!(compacted, count(used))
    for index in eachindex(nodes)
        if used[index]
            push!(compacted, nodes[index])
            remap[index] = Int32(length(compacted))
        end
    end
    tetrahedra = [(remap[t[1]], remap[t[2]], remap[t[3]], remap[t[4]]) for t in tetrahedra]
    surface_elements = [
        (attribute, (remap[f[1]], remap[f[2]], remap[f[3]])) for
        (attribute, f) in surface_elements
    ]
    nodes = compacted

    names = physical_names(spec)
    valid_attributes = Set(a for (d, a, _) in names if d == 2)
    all(element -> element[1] in valid_attributes, surface_elements) ||
        error("Surface element with an attribute outside the table")
    surface_counts = Dict{String, Int}()
    surface_areas = Dict{String, Float64}()
    for (attribute, f) in surface_elements
        key = string(attribute)
        surface_counts[key] = get(surface_counts, key, 0) + 1
        surface_areas[key] =
            get(surface_areas, key, 0.0) + face_area(nodes[f[1]], nodes[f[2]], nodes[f[3]])
    end
    volume_counts = Dict{String, Int}()
    volumes = Dict{String, Float64}()
    for (index, t) in enumerate(tetrahedra)
        key = string(tetrahedron_attribute[index])
        volume_counts[key] = get(volume_counts, key, 0) + 1
        volumes[key] =
            get(volumes, key, 0.0) +
            signed_volume(nodes[t[1]], nodes[t[2]], nodes[t[3]], nodes[t[4]])
    end
    for (d, a, name) in names
        d == 2 &&
            a != 9 &&
            !haskey(surface_counts, string(a)) &&
            error("Physical surface $name (attribute $a) is empty")
    end

    mkpath(dirname(output))
    open(output, "w") do stream
        print(stream, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n")
        written_names = filter(
            entry -> entry[1] == 3 || haskey(surface_counts, string(entry[2])),
            names
        )
        print(stream, "\$PhysicalNames\n$(length(written_names))\n")
        for (dimension, tag, name) in written_names
            print(stream, "$dimension $tag \"$name\"\n")
        end
        print(stream, "\$EndPhysicalNames\n")
        print(stream, "\$Nodes\n$(length(nodes))\n")
        for (index, point) in enumerate(nodes)
            @printf(stream, "%d %.16g %.16g %.16g\n", index, point[1], point[2], point[3])
        end
        print(stream, "\$EndNodes\n")
        print(stream, "\$Elements\n$(length(surface_elements) + length(tetrahedra))\n")
        element = 0
        for (attribute, f) in surface_elements
            element += 1
            print(stream, "$element 2 2 $attribute $attribute $(f[1]) $(f[2]) $(f[3])\n")
        end
        for (index, t) in enumerate(tetrahedra)
            element += 1
            attribute = Int(tetrahedron_attribute[index])
            print(
                stream,
                "$element 4 2 $attribute $attribute $(t[1]) $(t[2]) $(t[3]) $(t[4])\n"
            )
        end
        return print(stream, "\$EndElements\n")
    end
    merge!(
        manifest,
        Dict(
            "bytes" => filesize(output),
            "sha256" => bytes2hex(open(sha256, output)),
            "uncompacted_nodes" => uncompacted,
            "removed_unused_nodes" => uncompacted - length(nodes),
            "nodes" => length(nodes),
            "tetrahedra" => length(tetrahedra),
            "nonpositive" => 0,
            "volume_attribute_counts" => volume_counts,
            "volume_um3" => volumes,
            "surface_triangles" => length(surface_elements),
            "surface_attribute_counts" => surface_counts,
            "surface_area_um2" => surface_areas
        )
    )
    open(replace(output, r"\.msh2$" => ".json"), "w") do stream
        JSON.print(stream, manifest, 2)
        return println(stream)
    end
    verbose && println(
        "Saved ",
        output,
        ": nodes ",
        length(nodes),
        ", tetrahedra ",
        length(tetrahedra),
        ", surface triangles ",
        length(surface_elements),
        "\nSurface counts: ",
        surface_counts,
        "\nSurface areas (um^2): ",
        surface_areas,
        "\nVolumes (um^3): ",
        volumes
    )
    return manifest
end

end # module
