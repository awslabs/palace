# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Independent validation of a fabricated window mesh (or of the transmon reference's mesh):
# every tagged surface face is checked against the volumes on its two sides by its physical
# NAME — `exterior_boundary`: exactly one volume; `<conductor>_air`: one vacuum, no substrate;
# `<conductor>_substrate`: one substrate, no vacuum; `substrate_air` / `substrate_backside`:
# one of each; `lumped_element`: two vacuum — plus first-order simplices only, no face in two
# physical groups, no empty group, no nonpositive tetrahedron (minSICN), and the node /
# element / per-attribute counts of the generation manifest (MESH.json) when present. The
# per-attribute areas and volumes are recomputed from the node coordinates so two meshes of
# one geometry can be compared attribute by attribute. Writes MESH.validation.json.
#
#   julia --project=. validate_window_mesh.jl MESH.msh2

import Gmsh: gmsh
using JSON
using SHA
using Statistics

length(ARGS) == 1 || error("Usage: validate_window_mesh.jl MESH.msh2")
const MESH = abspath(ARGS[1])

function sorted_face(a::UInt64, b::UInt64, c::UInt64)
    a > b && ((a, b) = (b, a))
    b > c && ((b, c) = (c, b))
    a > b && ((a, b) = (b, a))
    return (a, b, c)
end

function element_blocks(dimension, physical_attribute)
    blocks = Tuple{Vector{UInt64}, Int}[]
    for entity in gmsh.model.get_entities_for_physical_group(dimension, physical_attribute)
        element_types, _, nodes_by_type = gmsh.model.mesh.get_elements(dimension, entity)
        for (element_type, element_nodes) in zip(element_types, nodes_by_type)
            _, element_dimension, _, node_count, _, primary =
                gmsh.model.mesh.get_element_properties(element_type)
            element_dimension == dimension || continue
            primary == dimension + 1 ||
                error("Expected first-order simplices for attribute $physical_attribute")
            node_count == primary ||
                error("Expected first-order elements for attribute $physical_attribute")
            push!(blocks, (element_nodes, node_count))
        end
    end
    return blocks
end

# Adjacency rule of a physical surface by its name: (substrate count, vacuum count) or the
# exterior rule (one volume of either kind).
function adjacency_rule(name)
    name == "exterior_boundary" && return :exterior
    name in ("substrate_air", "substrate_backside") && return (1, 1)
    name == "lumped_element" && return (0, 2)
    endswith(name, "_air") && return (0, 1)
    endswith(name, "_substrate") && return (1, 0)
    return error("No adjacency rule for physical surface $name")
end

function face_area(p1, p2, p3)
    ux, uy, uz = p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]
    vx, vy, vz = p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3]
    return 0.5 * hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx)
end

function tetrahedron_volume(p1, p2, p3, p4)
    ax, ay, az = p2[1] - p1[1], p2[2] - p1[2], p2[3] - p1[3]
    bx, by, bz = p3[1] - p1[1], p3[2] - p1[2], p3[3] - p1[3]
    cx, cy, cz = p4[1] - p1[1], p4[2] - p1[2], p4[3] - p1[3]
    return (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx)) /
           6
end

gmsh.initialize()
try
    gmsh.option.set_number("General.Verbosity", 2)
    gmsh.open(MESH)

    groups = Dict(
        gmsh.model.get_physical_name(dimension, tag) => (dimension, tag) for
        (dimension, tag) in gmsh.model.get_physical_groups()
    )
    get(groups, "substrate", nothing) == (3, 1) || error("Expected substrate = 3D attribute 1")
    get(groups, "vacuum", nothing) == (3, 2) || error("Expected vacuum = 3D attribute 2")
    surface_groups = sort!([(tag, name) for (name, (dimension, tag)) in groups if dimension == 2])
    isempty(surface_groups) && error("No physical surfaces")
    rules = Dict(tag => adjacency_rule(name) for (tag, name) in surface_groups)
    name_of = Dict(tag => name for (tag, name) in surface_groups)

    node_tags, coordinates, _ = gmsh.model.mesh.get_nodes()
    node_xyz = Dict{UInt64, NTuple{3, Float64}}()
    sizehint!(node_xyz, length(node_tags))
    for (index, tag) in enumerate(node_tags)
        node_xyz[tag] =
            (coordinates[3index - 2], coordinates[3index - 1], coordinates[3index])
    end

    surface_index = Dict{NTuple{3, UInt64}, Int32}()
    surface_attribute = Int32[]
    substrate_adjacent = UInt8[]
    vacuum_adjacent = UInt8[]
    surface_counts = Dict{String, Int}()
    surface_areas = Dict{String, Float64}()
    for (tag, name) in surface_groups
        count = 0
        area = 0.0
        for (element_nodes, node_count) in element_blocks(2, tag)
            for offset = 0:node_count:(length(element_nodes) - node_count)
                a, b, c = element_nodes[offset + 1], element_nodes[offset + 2], element_nodes[offset + 3]
                face = sorted_face(a, b, c)
                haskey(surface_index, face) &&
                    error("Surface face in more than one physical group: $face ($name)")
                push!(surface_attribute, tag)
                push!(substrate_adjacent, 0)
                push!(vacuum_adjacent, 0)
                surface_index[face] = Int32(length(surface_attribute))
                area += face_area(node_xyz[a], node_xyz[b], node_xyz[c])
                count += 1
            end
        end
        count > 0 || error("Physical surface $name (attribute $tag) is empty")
        surface_counts[string(tag)] = count
        surface_areas[string(tag)] = area
    end

    volume_counts = Dict{String, Int}()
    volumes = Dict{String, Float64}()
    for attribute = 1:2
        count = 0
        volume = 0.0
        for (element_nodes, node_count) in element_blocks(3, attribute)
            for offset = 0:node_count:(length(element_nodes) - node_count)
                a, b, c, d = element_nodes[offset + 1], element_nodes[offset + 2],
                element_nodes[offset + 3], element_nodes[offset + 4]
                for face in (
                    sorted_face(a, b, c),
                    sorted_face(a, b, d),
                    sorted_face(a, c, d),
                    sorted_face(b, c, d)
                )
                    index = get(surface_index, face, Int32(0))
                    index == 0 && continue
                    attribute == 1 ? (substrate_adjacent[index] += 1) :
                    (vacuum_adjacent[index] += 1)
                end
                volume += abs(tetrahedron_volume(node_xyz[a], node_xyz[b], node_xyz[c], node_xyz[d]))
                count += 1
            end
        end
        volume_counts[string(attribute)] = count
        volumes[string(attribute)] = volume
    end

    adjacency_errors = Dict{String, Int}(name => 0 for (_, name) in surface_groups)
    for index in eachindex(surface_attribute)
        rule = rules[surface_attribute[index]]
        s, v = substrate_adjacent[index], vacuum_adjacent[index]
        ok = rule == :exterior ? (s + v == 1) : ((s, v) == rule)
        ok || (adjacency_errors[name_of[surface_attribute[index]]] += 1)
    end
    all(iszero, values(adjacency_errors)) ||
        error("Physical surface adjacency errors: $adjacency_errors")

    _, volume_tags_by_type, _ = gmsh.model.mesh.get_elements(3)
    volume_tags = reduce(vcat, volume_tags_by_type; init=UInt64[])
    qualities = gmsh.model.mesh.get_element_qualities(volume_tags, "minSICN")

    manifest_path = replace(MESH, r"\.msh2$" => ".json")
    manifest = isfile(manifest_path) ? JSON.parsefile(manifest_path) : Dict()
    length(node_tags) == get(manifest, "nodes", length(node_tags)) ||
        error("Node count does not match the generation manifest")
    length(volume_tags) == get(manifest, "tetrahedra", length(volume_tags)) ||
        error("Tetrahedron count does not match the generation manifest")
    manifest_surface_counts = get(manifest, "surface_attribute_counts", surface_counts)
    all(
        get(surface_counts, attribute, -1) == count for
        (attribute, count) in manifest_surface_counts
    ) || error("Surface counts do not match the generation manifest")
    volume_counts == get(manifest, "volume_attribute_counts", volume_counts) ||
        error("Volume counts do not match the generation manifest")

    validation = Dict(
        "mesh" => MESH,
        "bytes" => filesize(MESH),
        "sha256" => bytes2hex(open(sha256, MESH)),
        "nodes" => length(node_tags),
        "tetrahedra" => length(volume_tags),
        "surface_triangles" => length(surface_attribute),
        "physical_attributes" => Dict(
            name => Dict("dimension" => value[1], "attribute" => value[2]) for
            (name, value) in groups
        ),
        "volume_attribute_counts" => volume_counts,
        "volume_um3" => volumes,
        "surface_attribute_counts" => surface_counts,
        "surface_area_um2" => surface_areas,
        "surface_adjacency_errors" => adjacency_errors,
        "mesh_quality" => Dict(
            "metric" => "minSICN",
            "minimum" => minimum(qualities),
            "p01" => quantile(qualities, 0.01),
            "median" => median(qualities),
            "below_0.001" => count(<(0.001), qualities),
            "below_0.01" => count(<(0.01), qualities),
            "nonpositive" => count(<=(0.0), qualities)
        )
    )
    validation["mesh_quality"]["nonpositive"] == 0 ||
        error("Mesh contains nonpositive tetrahedra")
    output = replace(MESH, r"\.msh2$" => ".validation.json")
    open(output, "w") do stream
        JSON.print(stream, validation, 2)
        println(stream)
    end
    println("Validated: ", MESH)
    println("Nodes: ", length(node_tags), ", tetrahedra: ", length(volume_tags))
    println("Surface adjacencies: ", adjacency_errors)
    println("Surface areas (um^2): ", surface_areas)
    println("Volumes (um^3): ", volumes)
    println("Mesh quality: ", validation["mesh_quality"])
    println("Saved: ", output)
finally
    gmsh.finalize()
end
