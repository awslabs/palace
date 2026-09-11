# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Export the actual planar etched footprint of a retained coupon. This preserves
# its discretized artificial collar instead of silently replacing it with a
# different analytic offset. No new process definition or ownership is inferred.
using Gmsh: gmsh
using SHA
using TOML

function main()
    length(ARGS) == 4 || error("retained.msh floor_z plane_z output-boundary.csv")
    input = abspath(ARGS[1])
    floor_z, plane_z = parse.(Float64, ARGS[2:3])
    all(isfinite, (floor_z, plane_z)) && floor_z != plane_z || error("Invalid process planes")
    output = abspath(ARGS[4])
    isfile(output) && error("Refuse to overwrite an etch footprint")
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity", 1)
        gmsh.open(input)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        coordinates = Dict(tag => ntuple(d -> xyz[3i - 3 + d], 3) for (i, tag) in enumerate(tags))
        scale = max(1.0, maximum(abs, xyz))
        tolerance = 1e-11 * scale
        edges = Dict{Tuple{UInt64,UInt64},Vector{Tuple{UInt64,UInt64}}}()
        area = 0.0
        for (_, attribute) in gmsh.model.getPhysicalGroups(2)
            3100 <= attribute < 4000 || continue
            for entity in gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
                types, _, connectivity = gmsh.model.mesh.getElements(2, entity)
                for (type, nodes) in zip(types, connectivity)
                    _, _, _, nnode, local_nodes, primary = gmsh.model.mesh.getElementProperties(type)
                    primary in (3, 4) || error("Unsupported floor face geometry")
                    for i in 1:nnode:length(nodes)
                        ids = collect(nodes[i:i + primary - 1])
                        points = [coordinates[id] for id in ids]
                        all(abs(p[3] - floor_z) <= tolerance for p in points) || continue
                        for j in 1:nnode
                            u, v = local_nodes[2j - 1], local_nodes[2j]
                            weights = primary == 3 ? (1 - u - v, u, v) :
                                      ((1 - u) * (1 - v) / 4, (1 + u) * (1 - v) / 4,
                                       (1 + u) * (1 + v) / 4, (1 - u) * (1 + v) / 4)
                            actual = coordinates[nodes[i + j - 1]]
                            all(abs(actual[d] - sum(weights[k] * points[k][d] for k in 1:primary)) <= tolerance
                                for d in 1:3) || error("Curved floor faces require a curved-footprint exporter")
                        end
                        signed = sum(points[j][1] * points[mod1(j + 1, primary)][2] -
                                     points[mod1(j + 1, primary)][1] * points[j][2]
                                     for j in 1:primary) / 2
                        abs(signed) > 0 || error("Degenerate floor face")
                        signed < 0 && reverse!(ids)
                        area += abs(signed)
                        for j in 1:primary
                            u, v = ids[j], ids[mod1(j + 1, primary)]
                            push!(get!(edges, minmax(u, v), Tuple{UInt64,UInt64}[]), (u, v))
                        end
                    end
                end
            end
        end
        isempty(edges) && error("No etched floor faces at requested elevation")
        next_node = Dict{UInt64,UInt64}()
        for adjacent in values(edges)
            if length(adjacent) == 2
                adjacent[1] == reverse(adjacent[2]) || error("Inconsistent floor orientation")
            elseif length(adjacent) == 1
                u, v = only(adjacent)
                haskey(next_node, u) && error("Nonmanifold or touching footprint loops")
                next_node[u] = v
            else
                error("Nonmanifold etched floor edge")
            end
        end
        Set(keys(next_node)) == Set(values(next_node)) || error("Open footprint boundary")
        loops = Vector{NTuple{2,Float64}}[]
        while !isempty(next_node)
            first = minimum(keys(next_node))
            node = first
            loop = NTuple{2,Float64}[]
            while true
                p = coordinates[node]
                push!(loop, (p[1], p[2]))
                node = pop!(next_node, node)
                node == first && break
            end
            # Drop only geometrically collinear vertices; retain every actual collar turn.
            changed = true
            while changed && length(loop) > 3
                changed = false
                for i in eachindex(loop)
                    a, b, c = loop[mod1(i - 1, length(loop))], loop[i], loop[mod1(i + 1, length(loop))]
                    ab, ac = (b[1] - a[1], b[2] - a[2]), (c[1] - a[1], c[2] - a[2])
                    cross = abs(ab[1] * ac[2] - ab[2] * ac[1])
                    if cross <= tolerance * hypot(ac...) &&
                       ab[1] * ac[1] + ab[2] * ac[2] >= 0 &&
                       ab[1]^2 + ab[2]^2 <= ac[1]^2 + ac[2]^2
                        deleteat!(loop, i)
                        changed = true
                        break
                    end
                end
            end
            push!(loops, loop)
        end
        areas = [sum(p[i][1] * p[mod1(i + 1, length(p))][2] -
                     p[mod1(i + 1, length(p))][1] * p[i][2] for i in eachindex(p)) / 2
                 for p in loops]
        abs(sum(areas) - area) <= 1e-9 * area || error("Footprint area changed during extraction")
        open(output, "w") do stream
            println(stream, "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y")
            for (index, loop) in enumerate(loops), (vertex, p) in enumerate(loop)
                println(stream, join((index, vertex, 1, plane_z, Int(areas[index] < 0), "Physical", p...), ','))
            end
        end
        open(output * ".provenance.toml", "w") do stream
            TOML.print(stream, Dict("Input" => input, "InputSHA256" => bytes2hex(sha256(read(input))),
                "OutputSHA256" => bytes2hex(sha256(read(output))), "FloorZ" => floor_z,
                "PlaneZ" => plane_z, "FloorArea" => area, "LoopArea" => sum(areas),
                "Loops" => length(loops), "Vertices" => sum(length, loops),
                "Scope" => "Exact retained mesh footprint, with collinear vertices removed"))
        end
    finally
        gmsh.finalize()
    end
    return nothing
end

main()
