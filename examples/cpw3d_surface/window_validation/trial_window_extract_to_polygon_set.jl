# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# TRIAL converter (sizing / robustness trial only, supervisor direction on USER decision 184):
# a lane-W WindowExtract (window_extract.py: chip-wide perimeter Loops per plane with Parent /
# Depth / MetalInside / Body, Bodies with their chip conductor, paired bump Footprints) to the
# mesher's polygon set (SCHEMA.md). The chip's loops are clipped at the window box HERE (Gmsh
# OCC intersection of every metal body with the box rectangle). Lane W owns the real polygon
# set (termination of cut traces, bump-aware bodies): meshes from this converter are TRIAL
# meshes, never a stage-1 reference.
#
#   julia --project=. trial_window_extract_to_polygon_set.jl EXTRACT.json OUTPUT.json
#       [--l1-substrate 525.0] [--l2-substrate 525.0] [--metal-planes 0.0,4.8]
#
# Metal planes: the two planes (by z) whose loops carry metal bodies; the lower one faces up,
# the upper one faces down. A plane's unbounded ground (GroundBodies) is the box minus its
# depth-0 loops with MetalInside false; a closed loop with MetalInside true is a body whose
# metal is the loop minus its direct child loops. Ground bodies (Kind Ground) are `ground`;
# every other body is the terminal `body_<chip conductor>`. Bumps: footprints with Inside
# true, one per partner pair, the lower plane's polygon; a footprint straddling the wall is
# refused.

import Gmsh: gmsh
using JSON

const Point2 = NTuple{2, Float64}

function parse_arguments(args)
    options = Dict{String, String}()
    positional = String[]
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            options[args[i][3:end]] = args[i + 1]
            i += 2
        else
            push!(positional, args[i])
            i += 1
        end
    end
    length(positional) == 2 || error(
        "Usage: trial_window_extract_to_polygon_set.jl EXTRACT.json OUTPUT.json " *
        "[--l1-substrate T] [--l2-substrate T] [--metal-planes z1,z2]"
    )
    return positional, options
end

points(values) = Point2[(Float64(p[1]), Float64(p[2])) for p in values]

function signed_area(polygon::Vector{Point2})
    n = length(polygon)
    return 0.5 * sum(
        polygon[i][1] * polygon[mod1(i + 1, n)][2] -
        polygon[mod1(i + 1, n)][1] * polygon[i][2] for i = 1:n
    )
end

function add_face(polygon::Vector{Point2}, holes::Vector{Vector{Point2}})
    occ = gmsh.model.occ
    function loop(vertices)
        signed_area(vertices) < 0.0 && (vertices = reverse(vertices))
        tags = [occ.add_point(p[1], p[2], 0.0) for p in vertices]
        lines = [
            occ.add_line(tags[i], tags[mod1(i + 1, length(tags))]) for i in eachindex(tags)
        ]
        return occ.add_curve_loop(lines)
    end
    return occ.add_plane_surface(vcat(loop(polygon), [loop(h) for h in holes]))
end

# Vertices of a face's loops in order (straight edges only), snapped to the box wall.
function face_polygons(surface::Int32, box; wall_tolerance=1.0e-6)
    snap(v, lo, hi) =
        abs(v - lo) <= wall_tolerance ? lo : (abs(v - hi) <= wall_tolerance ? hi : v)
    loops = Vector{Point2}[]
    _, curve_tags = gmsh.model.occ.get_curve_loops(surface)
    for curves in curve_tags
        edges = Tuple{Point2, Point2}[]
        for tag in curves
            bounds = gmsh.model.get_parametrization_bounds(1, abs(tag))
            a = gmsh.model.get_value(1, abs(tag), [bounds[1][1]])
            b = gmsh.model.get_value(1, abs(tag), [bounds[2][1]])
            push!(edges, ((a[1], a[2]), (b[1], b[2])))
        end
        # Chain the edges by endpoint coincidence.
        used = falses(length(edges))
        chain = Point2[edges[1][1], edges[1][2]]
        used[1] = true
        while count(used) < length(edges)
            tail = chain[end]
            found = false
            for (k, (p, q)) in enumerate(edges)
                used[k] && continue
                if hypot(p[1] - tail[1], p[2] - tail[2]) <= 1.0e-7
                    push!(chain, q)
                    used[k] = true
                    found = true
                    break
                elseif hypot(q[1] - tail[1], q[2] - tail[2]) <= 1.0e-7
                    push!(chain, p)
                    used[k] = true
                    found = true
                    break
                end
            end
            found || error("Curve loop of surface $surface does not chain")
        end
        hypot(chain[end][1] - chain[1][1], chain[end][2] - chain[1][2]) <= 1.0e-7 ||
            error("Curve loop of surface $surface is not closed")
        pop!(chain)
        polygon =
            Point2[(snap(p[1], box[1], box[2]), snap(p[2], box[3], box[4])) for p in chain]
        # Drop repeated vertices produced by snapping.
        cleaned = Point2[]
        for p in polygon
            isempty(cleaned) ||
                hypot(p[1] - cleaned[end][1], p[2] - cleaned[end][2]) > 1.0e-9 ||
                continue
            push!(cleaned, p)
        end
        while length(cleaned) > 1 &&
            hypot(cleaned[end][1] - cleaned[1][1], cleaned[end][2] - cleaned[1][2]) <=
            1.0e-9
            pop!(cleaned)
        end
        length(cleaned) >= 3 || error("Degenerate loop on surface $surface")
        push!(loops, cleaned)
    end
    order = sortperm([abs(signed_area(l)) for l in loops]; rev=true)
    return loops[order[1]], loops[order[2:end]]
end

# Clip a body (outer loop minus hole loops) at the box; returns (outer, holes) per component.
function clip_body(outer::Vector{Point2}, holes::Vector{Vector{Point2}}, box)
    gmsh.initialize()
    try
        gmsh.option.set_number("General.Verbosity", 1)
        gmsh.model.add("clip")
        occ = gmsh.model.occ
        rectangle = occ.add_rectangle(box[1], box[3], 0.0, box[2] - box[1], box[4] - box[3])
        face = add_face(outer, holes)
        pieces, _ = occ.intersect([(2, face)], [(2, rectangle)], -1, true, true)
        occ.synchronize()
        return [face_polygons(tag, box) for (dimension, tag) in pieces if dimension == 2]
    finally
        gmsh.finalize()
    end
end

function bbox_meets(bbox, box)
    return bbox[1] < box[2] && bbox[2] > box[1] && bbox[3] < box[4] && bbox[4] > box[3]
end

function main()
    positional, options = parse_arguments(ARGS)
    extract = JSON.parsefile(positional[1])
    window = extract["WindowExtract"]
    box = (
        Float64(window["Box"][1]),
        Float64(window["Box"][2]),
        Float64(window["Box"][3]),
        Float64(window["Box"][4])
    )
    loops = window["Loops"]
    bodies = window["Bodies"]
    ground_bodies = Dict(parse(Float64, k) => v for (k, v) in window["GroundBodies"])
    metal_planes = if haskey(options, "metal-planes")
        sort!(parse.(Float64, split(options["metal-planes"], ",")))
    else
        zs = sort!(unique(Float64(b["Plane"]) for b in values(bodies) if b["Closed"] === true))
        length(zs) == 2 ||
            error("Give --metal-planes z1,z2 (closed bodies on planes $zs)")
        zs
    end
    length(metal_planes) == 2 || error("Two metal planes are needed for the trial")
    substrates = [
        parse(Float64, get(options, "l1-substrate", "525.0")),
        parse(Float64, get(options, "l2-substrate", "525.0"))
    ]
    by_index = Dict(l["Index"] => l for l in loops)
    children = Dict{Any, Vector{Any}}()
    for l in loops
        push!(get!(children, l["Parent"], Any[]), l)
    end
    terminals = String[]
    planes = Any[]
    for (k, z) in enumerate(metal_planes)
        plane_loops = [l for l in loops if isapprox(Float64(l["Plane"]), z; atol=1.0e-6)]
        polygons = Any[]
        function emit(label, outer, holes)
            pieces = clip_body(outer, holes, box)
            for (clipped_outer, clipped_holes) in pieces
                push!(
                    polygons,
                    Dict(
                        "Conductor" => label,
                        "Outer" => [[p[1], p[2]] for p in clipped_outer],
                        "Holes" => [[[p[1], p[2]] for p in h] for h in clipped_holes]
                    )
                )
            end
            return length(pieces)
        end
        # The unbounded ground: the box minus the depth-0 non-metal loops.
        if haskey(ground_bodies, z)
            outer = Point2[
                (box[1], box[3]),
                (box[2], box[3]),
                (box[2], box[4]),
                (box[1], box[4])
            ]
            holes = Vector{Point2}[]
            for l in plane_loops
                l["Parent"] === nothing && l["MetalInside"] == false || continue
                bbox_meets(l["BBox"], box) || continue
                # A depth-0 hole larger than the box: clip it separately by subtracting.
                push!(holes, points(l["Polygon"]))
            end
            # Holes may cross the box wall: clip the ground as (box minus holes) via OCC cut.
            gmsh.initialize()
            try
                gmsh.option.set_number("General.Verbosity", 1)
                gmsh.model.add("ground")
                occ = gmsh.model.occ
                rectangle =
                    occ.add_rectangle(box[1], box[3], 0.0, box[2] - box[1], box[4] - box[3])
                tools = [(2, add_face(h, Vector{Point2}[])) for h in holes]
                pieces =
                    isempty(tools) ? [(Int32(2), rectangle)] :
                    occ.cut([(2, rectangle)], tools, -1, true, true)[1]
                occ.synchronize()
                for (dimension, tag) in pieces
                    dimension == 2 || continue
                    clipped_outer, clipped_holes = face_polygons(tag, box)
                    push!(
                        polygons,
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[p[1], p[2]] for p in clipped_outer],
                            "Holes" => [[[p[1], p[2]] for p in h] for h in clipped_holes]
                        )
                    )
                end
            finally
                gmsh.finalize()
            end
        end
        for l in plane_loops
            l["Closed"] == true && l["MetalInside"] == true || continue
            bbox_meets(l["BBox"], box) || continue
            body = bodies[string(l["Body"])]
            label = body["Kind"] == "Ground" ? "ground" : "body_$(body["Conductor"])"
            holes = Vector{Point2}[
                points(c["Polygon"]) for c in get(children, l["Index"], Any[]) if
                c["MetalInside"] == false && bbox_meets(c["BBox"], box)
            ]
            pieces = emit(label, points(l["Polygon"]), holes)
            pieces > 0 &&
                label != "ground" &&
                !(label in terminals) &&
                push!(terminals, label)
        end
        isempty(polygons) && error("No metal on plane $z inside the box")
        push!(
            planes,
            Dict(
                "Name" => k == 1 ? "L1" : "L2",
                "SurfaceZ" => z,
                "Facing" => k == 1 ? "up" : "down",
                "SubstrateThickness" => substrates[k],
                "Polygons" => polygons
            )
        )
    end
    bumps = Any[]
    for f in window["Footprints"]
        isapprox(Float64(f["Plane"]), metal_planes[1]; atol=1.0e-6) || continue
        f["Partner"] === nothing && continue
        f["Straddles"] == true &&
            error("Bump footprint $(f["Index"]) straddles the box wall")
        f["Inside"] == true || continue
        body = bodies[string(f["Body"])]
        label = body["Kind"] == "Ground" ? "ground" : "body_$(body["Conductor"])"
        push!(bumps, Dict("Conductor" => label, "Footprint" => f["Polygon"]))
    end
    sort!(terminals)
    set = Dict(
        "Version" => 1,
        "Name" => "TRIAL-$(window["Chip"])-$(window["Window"])",
        "Box" => Dict("X" => [box[1], box[2]], "Y" => [box[3], box[4]]),
        "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
        "Planes" => planes,
        "Bumps" => bumps,
        "Vacuum" => Dict("Below" => 0.0, "Above" => 0.0),
        "Terminals" => terminals,
        "Source" => Dict(
            "Extract" => abspath(positional[1]),
            "Trial" => "lane-F sizing / robustness trial; lane W owns the real polygon set"
        )
    )
    open(positional[2], "w") do io
        return JSON.print(io, set, 1)
    end
    for plane in planes
        println(
            plane["Name"],
            " z=",
            plane["SurfaceZ"],
            ": ",
            [
                (p["Conductor"], length(p["Outer"]), length(p["Holes"])) for
                p in plane["Polygons"]
            ]
        )
    end
    return println(
        "Bumps: ",
        length(bumps),
        ", terminals: ",
        terminals,
        " -> ",
        positional[2]
    )
end

main()
