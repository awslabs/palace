# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Census of a hybrid Gmsh mesh held in the current model: element counts by type,
# quality per type (Gmsh minSJ / minSICN), orientation (positive Jacobian of every
# linear element by its own volume formula), physical-group areas. Report only, no
# gate. SPIKE code (decision 37), not production.

import Gmsh: gmsh
using LinearAlgebra

const GMSH_TYPE_NAMES = Dict(1 => "Line", 2 => "Triangle", 3 => "Quadrangle", 4 => "Tetrahedron",
                             5 => "Hexahedron", 6 => "Prism", 7 => "Pyramid")

function tetrahedron_volume(p)
    return dot(p[:, 2] .- p[:, 1], cross(p[:, 3] .- p[:, 1], p[:, 4] .- p[:, 1])) / 6.0
end

# Signed volume of a linear element by decomposition into tetrahedra, each of
# which must itself be positive for the element to count as positively oriented.
function element_orientation_volumes(type, p)
    if type == 4
        return [tetrahedron_volume(p)]
    elseif type == 6
        # Prism (0,1,2 | 3,4,5): the three tetrahedra of the standard split.
        return [tetrahedron_volume(p[:, [1, 2, 3, 4]]), tetrahedron_volume(p[:, [2, 3, 4, 5]]),
                tetrahedron_volume(p[:, [3, 4, 5, 6]])]
    elseif type == 7
        # Pyramid (0,1,2,3 | 4): the two tetrahedra of each base diagonal (the
        # four together cover the pyramid twice).
        return [tetrahedron_volume(p[:, [1, 2, 3, 5]]), tetrahedron_volume(p[:, [1, 3, 4, 5]]),
                tetrahedron_volume(p[:, [1, 2, 4, 5]]), tetrahedron_volume(p[:, [2, 3, 4, 5]])]
    end
    error("unsupported volume element type $type")
end

function hybrid_volume_census()
    node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
    index = Dict{Int, Int}(Int(tag) => i for (i, tag) in enumerate(node_tags))
    xyz = reshape(coordinates, 3, :)
    types, tags, nodes = gmsh.model.mesh.getElements(3)
    census = Dict{String, Any}()
    total = 0
    for (type, element_tags, connectivity) in zip(types, tags, nodes)
        _, _, _, width, _, _ = gmsh.model.mesh.getElementProperties(type)
        conn = reshape(connectivity, Int(width), :)
        element_count = size(conn, 2)
        total += element_count
        negative = 0
        minimum_volume = Inf
        for c in axes(conn, 2)
            p = xyz[:, [index[Int(conn[r, c])] for r in axes(conn, 1)]]
            volumes = element_orientation_volumes(Int(type), p)
            minimum_volume = min(minimum_volume, sum(volumes))
            all(>(0.0), volumes) || (negative += 1)
        end
        sj = gmsh.model.mesh.getElementQualities(element_tags, "minSJ")
        sicn = gmsh.model.mesh.getElementQualities(element_tags, "minSICN")
        gamma = gmsh.model.mesh.getElementQualities(element_tags, "gamma")
        census[GMSH_TYPE_NAMES[Int(type)]] = Dict{String, Any}(
            "Count" => element_count, "NonPositive" => negative, "MinimumVolume" => minimum_volume,
            "MinScaledJacobian" => minimum(sj), "MedianScaledJacobian" => sort(sj)[(end + 1) ÷ 2],
            "MinSICN" => minimum(sicn), "MinGamma" => minimum(gamma),
            "ScaledJacobianBelow0.01" => count(<(0.01), sj))
    end
    census["Total"] = total
    census["Nodes"] = length(node_tags)
    return census
end

function physical_surface_areas()
    areas = Dict{String, Any}()
    for (dim, tag) in gmsh.model.getPhysicalGroups(2)
        area = 0.0
        counts = Dict{String, Int}()
        for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            types, element_tags, _ = gmsh.model.mesh.getElements(2, entity)
            for (type, tags) in zip(types, element_tags)
                counts[GMSH_TYPE_NAMES[Int(type)]] = get(counts, GMSH_TYPE_NAMES[Int(type)], 0) +
                                                     length(tags)
                # The "volume" quality of a surface element is its area.
                area += sum(gmsh.model.mesh.getElementQualities(tags, "volume"))
            end
        end
        areas[string(tag)] = Dict{String, Any}("Area" => area, "Elements" => counts)
    end
    return areas
end

# Gmsh's tetrahedral optimizer (edge swaps) has been observed to leave exact
# duplicate tetrahedra next to pyramids; remove them and report how many. The
# conformity of the result is checked independently (check_hybrid_mesh.py).
function remove_duplicate_volume_elements!()
    before = sum(length(tags) for tags in gmsh.model.mesh.getElements(3)[2]; init=0)
    gmsh.model.mesh.removeDuplicateElements([(3, tag) for (dim, tag) in gmsh.model.getEntities(3)])
    after = sum(length(tags) for tags in gmsh.model.mesh.getElements(3)[2]; init=0)
    return before - after
end

# Minimal JSON writer (Dict / Vector / numbers / strings / Bool / nothing).
function write_census_json(io, value, indent=0)
    pad = " "^indent
    if value isa Dict
        keys_sorted = sort(collect(keys(value)); by=string)
        println(io, "{")
        for (i, key) in enumerate(keys_sorted)
            print(io, pad, "  \"", key, "\": ")
            write_census_json(io, value[key], indent + 2)
            println(io, i < length(keys_sorted) ? "," : "")
        end
        print(io, pad, "}")
    elseif value isa AbstractVector
        print(io, "[")
        for (i, item) in enumerate(value)
            write_census_json(io, item, indent + 2)
            i < length(value) && print(io, ", ")
        end
        print(io, "]")
    elseif value isa AbstractString || value isa Symbol
        print(io, "\"", value, "\"")
    elseif value isa Bool
        print(io, value ? "true" : "false")
    elseif value === nothing
        print(io, "null")
    elseif value isa Integer
        print(io, value)
    elseif value isa Real
        print(io, isfinite(value) ? repr(Float64(value)) : "null")
    else
        print(io, "\"", string(value), "\"")
    end
end

