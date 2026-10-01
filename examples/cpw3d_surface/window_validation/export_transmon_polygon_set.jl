# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Export the single transmon's fabricated-reference footprint as a polygon set for
# mesh_polygon_window.jl (the validation of the generalised mesher against the recorded
# transmon meshes, USER decision 184): the same DeviceLayout geometry, the `metal`, `port_1`
# and `port_2` plan groups fused into the two physical conductors (`ground` = the large
# component incl. the readout port patches, `island` = the small one) exactly as
# generate_shared_plan_mesh.jl, every fused boundary curve exported as a chain of vertices:
# a straight curve by its two endpoints (the mesher re-splits it at the tangent target) and a
# curved one by the nodes of its transfinite discretisation at the tangent target (the nodes
# the recorded generator placed on it: ceil(length / t) equal arc-length segments), so the
# polygon set carries the recorded plan constraints. The box is the DeviceLayout substrate
# footprint with its 525-um substrate, the 475-um vacuum below and the 1000-um vacuum above.
#
# Run in the SingleTransmon example environment (branch simlapointe/transmon-thick):
#   julia --project=../../../../DeviceLayout.jl/examples/SingleTransmon \
#       export_transmon_polygon_set.jl TANGENTIAL_UM OUTPUT.json

using DeviceLayout
using DeviceLayout.PreferredUnits
using JSON

const DEVICE_LAYOUT_ROOT = get(
    ENV,
    "DEVICE_LAYOUT_ROOT",
    normpath(joinpath(@__DIR__, "..", "..", "..", "..", "DeviceLayout.jl"))
)
include(joinpath(DEVICE_LAYOUT_ROOT, "examples", "SingleTransmon", "SingleTransmon.jl"))
using .SingleTransmon

length(ARGS) == 2 ||
    error("Usage: export_transmon_polygon_set.jl TANGENTIAL_UM OUTPUT.json")
const TANGENTIAL_UM = parse(Float64, ARGS[1])
const OUTPUT = abspath(ARGS[2])

_, solid_model = SingleTransmon.single_transmon(;
    include_bridges=false,
    surface_participation=true,
    mesh_order=1
)
const gmsh = SolidModels.gmsh
occ = gmsh.model.occ

group_surfaces(name) = SolidModels.dimtags(solid_model[name, 2])
inside_surfaces =
    unique(reduce(vcat, (group_surfaces(name) for name in ("metal", "port_1", "port_2"))))
bounds(name) = reduce(
    (a, b) -> (min.(a[1:3], b[1:3])..., max.(a[4:6], b[4:6])...),
    (
        gmsh.model.get_bounding_box(entity...) for
        entity in SolidModels.dimtags(solid_model[name, 3])
    )
)
substrate_bounds = bounds("substrate")
vacuum_bounds = bounds("vacuum")
box = (substrate_bounds[1], substrate_bounds[4], substrate_bounds[2], substrate_bounds[5])
substrate_thickness = 0.0 - substrate_bounds[3]
vacuum_below = substrate_bounds[3] - vacuum_bounds[3]
vacuum_above = vacuum_bounds[6]

tools = occ.copy(inside_surfaces)
fused, _ = occ.fuse([first(tools)], tools[2:end])
occ.synchronize()
all(entity -> entity[1] == 2, fused) || error("Fuse produced non-surface entities")
length(fused) == 2 || error("Expected two fused conductor faces, found $(length(fused))")

# Transfinite discretisation of every curved boundary curve at the tangent target (the
# recorded generator's rule), meshed in 1D so Gmsh places the nodes.
gmsh.model.set_visibility(gmsh.model.get_entities(), 0, true)
gmsh.model.set_visibility(fused, 1, true)
gmsh.option.set_number("Mesh.MeshOnlyVisible", 1)
gmsh.option.set_number("General.Verbosity", 2)
curves = unique(gmsh.model.get_boundary(fused, false, false, false))
curve_types = Dict(abs(tag) => gmsh.model.get_type(1, abs(tag)) for (_, tag) in curves)
for (_, tag) in curves
    tag = abs(tag)
    length_um = occ.get_mass(1, tag)
    segments = max(1, ceil(Int, length_um / TANGENTIAL_UM))
    gmsh.model.mesh.set_transfinite_curve(tag, segments + 1)
end
gmsh.model.mesh.generate(1)

function curve_endpoints(tag)
    points = gmsh.model.get_boundary([(1, tag)], false, false, false)
    length(points) == 2 || error("Curve $tag is closed or degenerate")
    return abs(points[1][2]), abs(points[2][2])
end
point_xy(tag) = Tuple(round.(gmsh.model.get_value(0, tag, Float64[])[1:2]; digits=9))

# Vertices of a curve from its start point to its end point (oriented as the OCC curve).
function curve_vertices(tag)
    start_point, end_point = curve_endpoints(tag)
    start_xy, end_xy = point_xy(start_point), point_xy(end_point)
    curve_types[tag] == "Line" && return [start_xy, end_xy]
    _, coordinates, parameters = gmsh.model.mesh.get_nodes(1, tag, false, true)
    isempty(parameters) && error("Curved curve $tag has no interior transfinite nodes")
    order = sortperm(parameters)
    interior = [
        (round(coordinates[3i - 2]; digits=9), round(coordinates[3i - 1]; digits=9)) for
        i in order
    ]
    # The curve parameter increases from its start to its end for OCC curves.
    lower, _ = gmsh.model.get_parametrization_bounds(1, tag)
    start_first = gmsh.model.get_value(1, tag, [lower[1]])
    hypot(start_first[1] - start_xy[1], start_first[2] - start_xy[2]) <= 1.0e-6 ||
        error("Curve $tag parametrisation does not start at its start point")
    return vcat([start_xy], interior, [end_xy])
end

signed_area(points) =
    0.5 * sum(
        points[i][1] * points[mod1(i + 1, length(points))][2] -
        points[mod1(i + 1, length(points))][1] * points[i][2] for i in eachindex(points)
    )

function face_loops(face)
    face_curves =
        [abs(tag) for (_, tag) in gmsh.model.get_boundary([face], false, false, false)]
    endpoints = Dict(tag => curve_endpoints(tag) for tag in face_curves)
    incident = Dict{Int32, Vector{Int32}}()
    for (tag, (a, b)) in endpoints
        push!(get!(incident, a, Int32[]), tag)
        push!(get!(incident, b, Int32[]), tag)
    end
    all(list -> length(list) == 2, values(incident)) ||
        error("Face $(face[2]) boundary is not a set of simple loops")
    remaining = Set(face_curves)
    loops = Vector{Vector{NTuple{2, Float64}}}()
    while !isempty(remaining)
        tag = first(remaining)
        loop = NTuple{2, Float64}[]
        point = endpoints[tag][1]
        while true
            delete!(remaining, tag)
            vertices = curve_vertices(tag)
            a, b = endpoints[tag]
            if point == b
                vertices = reverse(vertices)
                point = a
            else
                point = b
            end
            append!(loop, vertices[1:(end - 1)])
            candidates = filter(in(remaining), incident[point])
            isempty(candidates) && break
            tag = only(candidates)
        end
        length(loop) >= 3 || error("Degenerate loop on face $(face[2])")
        push!(loops, loop)
    end
    return loops
end

polygons = []
areas = [occ.get_mass(face...) for face in fused]
island_index = argmin(areas)
for (index, face) in enumerate(fused)
    loops = face_loops(face)
    order = sortperm([abs(signed_area(loop)) for loop in loops]; rev=true)
    outer = loops[order[1]]
    holes = loops[order[2:end]]
    loop_area = abs(signed_area(outer)) - sum(abs(signed_area(h)) for h in holes; init=0.0)
    println(
        "Face ",
        face[2],
        ": OCC area ",
        areas[index],
        ", polygon area ",
        loop_area,
        ", outer vertices ",
        length(outer),
        ", holes ",
        length(holes)
    )
    push!(
        polygons,
        Dict(
            "Conductor" => index == island_index ? "island" : "ground",
            "Outer" => [[p[1], p[2]] for p in outer],
            "Holes" => [[[p[1], p[2]] for p in h] for h in holes]
        )
    )
end
curve_type_counts = Dict{String, Int}()
for (_, kind) in curve_types
    curve_type_counts[kind] = get(curve_type_counts, kind, 0) + 1
end
println("Boundary curve types: ", curve_type_counts)

polygon_set = Dict(
    "Version" => 1,
    "Name" => "single-transmon-fabricated",
    "Source" => Dict(
        "DeviceLayoutRoot" => DEVICE_LAYOUT_ROOT,
        "TangentialUm" => TANGENTIAL_UM,
        "CurveTypes" => curve_type_counts,
        "SubstrateBounds" => collect(substrate_bounds),
        "VacuumBounds" => collect(vacuum_bounds)
    ),
    "Box" => Dict("X" => [box[1], box[2]], "Y" => [box[3], box[4]]),
    "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
    "Planes" => [
        Dict(
            "Name" => "L1",
            "SurfaceZ" => 0.0,
            "Facing" => "up",
            "SubstrateThickness" => substrate_thickness,
            "Polygons" => polygons
        )
    ],
    "Bumps" => [],
    "Vacuum" => Dict("Below" => vacuum_below, "Above" => vacuum_above),
    "Terminals" => ["island"]
)
mkpath(dirname(OUTPUT))
open(OUTPUT, "w") do stream
    JSON.print(stream, polygon_set, 2)
    return println(stream)
end
gmsh.finalize()
println("Saved ", OUTPUT)
