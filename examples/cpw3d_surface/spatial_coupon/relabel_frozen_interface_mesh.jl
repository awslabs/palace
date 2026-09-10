# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Apply the exact element-wise interface ownership partition to an already generated
# all-tetrahedral mesh. This permits fixed-surface volume generation to remain
# independent of bookkeeping slots while retaining the original response attributes.
using Gmsh: gmsh
using DelimitedFiles
using SHA
using TOML
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))
include(joinpath(@__DIR__, "label_interface_patches.jl"))

function file_sha256(path)
    return bytes2hex(sha256(read(path)))
end

function node_coordinates()
    tags, coordinates, _ = gmsh.model.mesh.getNodes()
    return Dict(tag => ntuple(d -> coordinates[3i - 3 + d], 3) for (i, tag) in enumerate(tags))
end

function dimension_elements(dimension)
    result = Dict{UInt64,Tuple{Int,Vector{UInt64}}}()
    types, tags, connectivity = gmsh.model.mesh.getElements(dimension)
    for (type, elements, nodes) in zip(types, tags, connectivity)
        _, _, _, nnode, _, _ = gmsh.model.mesh.getElementProperties(type)
        for (i, element) in enumerate(elements)
            haskey(result, element) && error("Duplicate element tag")
            result[element] = (Int(type), collect(nodes[(i - 1) * nnode + 1:i * nnode]))
        end
    end
    return result
end

function canonical_nodes()
    return sort!(collect(values(node_coordinates())))
end

function canonical_elements(dimension)
    coordinates = node_coordinates()
    records = Tuple{Int,Tuple}[]
    types, _, connectivity = gmsh.model.mesh.getElements(dimension)
    for (type, nodes) in zip(types, connectivity)
        _, _, order, nnode, _, primary = gmsh.model.mesh.getElementProperties(type)
        order == 1 && nnode == primary || error("Frozen relabeling requires linear elements")
        for i in 1:nnode:length(nodes)
            points = Tuple(sort!([coordinates[node] for node in nodes[i:i + nnode - 1]]))
            push!(records, (Int(type), points))
        end
    end
    return sort!(records)
end

function surface_elements()
    result = Dict{UInt64,Tuple{Int,Vector{UInt64}}}()
    for (_, attribute) in gmsh.model.getPhysicalGroups(2)
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
            types, tags, connectivity = gmsh.model.mesh.getElements(2, entity)
            for (type, elements, nodes) in zip(types, tags, connectivity)
                _, _, _, nnode, _, primary = gmsh.model.mesh.getElementProperties(type)
                primary == 3 || error("Physical surface contains a nontriangular element")
                for (i, element) in enumerate(elements)
                    haskey(result, element) && error("A surface element has multiple physical assignments")
                    result[element] = (Int(type), collect(nodes[(i - 1) * nnode + 1:i * nnode]))
                end
            end
        end
    end
    return result
end

function physical_family(attribute)
    family = div(attribute, 1000)
    family == 3 && return attribute >= 3100 ? 3100 : 3000
    family in (4, 5, 6) && return 1000family + mod(attribute, 100)
    return attribute
end

function grouped_surface_measures()
    result = Dict{Int,Float64}()
    for (_, attribute) in gmsh.model.getPhysicalGroups(2)
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2, attribute)
            types, tags, _ = gmsh.model.mesh.getElements(2, entity)
            for (type, elements) in zip(types, tags)
                points, weights = gmsh.model.mesh.getIntegrationPoints(type, "Gauss8")
                _, determinants, _ = gmsh.model.mesh.getJacobians(type, points, entity)
                nq = length(weights)
                length(determinants) == nq * length(elements) || error("Unexpected surface Jacobian data")
                for i in eachindex(elements)
                    measure = sum(weights[q] * determinants[(i - 1) * nq + q] for q in 1:nq)
                    key = physical_family(Int(attribute))
                    result[key] = get(result, key, 0.0) + measure
                end
            end
        end
    end
    return result
end

function check_measures(before, after)
    keys(before) == keys(after) || error("Grouped physical-family measures have different keys")
    worst = maximum(abs(after[key] - value) / max(abs(value), 1e-300) for (key, value) in before)
    worst <= 1e-11 || error("Grouped physical-family area changed by $worst")
    return worst
end

function main()
    length(ARGS) in (4, 6) ||
        error("signature_directory thin|fabricated input.msh output.msh [--expected-measures measures.csv]")
    root, input, output = abspath.((ARGS[1], ARGS[3], ARGS[4]))
    kind = ARGS[2]
    kind in ("thin", "fabricated") || error("Kind must be thin or fabricated")
    isfile(input) || error("Input mesh does not exist")
    isfile(output) && error("Refuse to overwrite output mesh")
    expected_path = nothing
    if length(ARGS) == 6
        ARGS[5] == "--expected-measures" || error("Unknown option $(ARGS[5])")
        expected_path = abspath(ARGS[6])
        isfile(expected_path) || error("Expected-measures file does not exist")
    end
    process_path = joinpath(root, "process.toml")
    process = TOML.parsefile(process_path)
    process["Units"] == "um" || error("Process units must be um")
    radius = Float64(process["Radius"])
    thickness = Float64(process["MetalThickness"])
    overetch = Float64(process["Overetch"])
    edges = read_edges(joinpath(root, "mesh-signature.csv"))
    loops = read_boundary(joinpath(root, "plan-view-boundary.csv"))
    expected = nothing
    if expected_path !== nothing
        data, _ = readdlm(expected_path, ',', header = true)
        expected = Set(Int(data[i, 2]) for i in axes(data, 1)
                       if Int(data[i, 1]) == 2 && Int(data[i, 2]) != 1)
        isempty(expected) && error("Expected interface-attribute set is empty")
    end

    gmsh.initialize()
    try
        gmsh.open(input)
        before_nodes = node_coordinates()
        before_volume = dimension_elements(3)
        before_surface = surface_elements()
        canonical_node_geometry = canonical_nodes()
        canonical_volume_geometry = canonical_elements(3)
        canonical_surface_geometry = canonical_elements(2)
        before_measures = grouped_surface_measures()
        volume_tags = collect(keys(before_volume))
        before_quality = minimum(gmsh.model.mesh.getElementQualities(volume_tags, "minSICN"))
        before_quality > 1e-10 || error("Input has a nonpositive or near-singular volume element")
        report = output * ".interface-partition.csv"
        label_interface_patches(
            edges,
            loops,
            radius,
            report;
            minimum_size = 0.0,
            fabricated = kind == "fabricated",
            metal_thickness = thickness,
            overetch = overetch
        )
        before_nodes == node_coordinates() || error("Node coordinates or IDs changed")
        before_volume == dimension_elements(3) || error("Volume connectivity changed")
        before_surface == surface_elements() || error("Surface element IDs or connectivity changed")
        measure_difference = check_measures(before_measures, grouped_surface_measures())
        actual = Set(Int(attribute) for (_, attribute) in gmsh.model.getPhysicalGroups(2)
                     if attribute != 1)
        expected === nothing || actual == expected ||
            error("Interface coverage differs: missing=$(setdiff(expected, actual)) extra=$(setdiff(actual, expected))")
        gmsh.option.setNumber("Mesh.Renumber", 0)
        # Palace/MFEM consumes MSH 2.2. Gmsh can renumber discrete entities while
        # writing that format, so the round-trip check below binds geometry and
        # connectivity by coordinates in addition to the stronger in-memory ID check.
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(output)
        output_sha256 = file_sha256(output)
        certificate = output * ".interface-partition.csv.elements.csv"
        open(output * ".partition-certificate.toml", "w") do stream
            TOML.print(stream, Dict(
                "Version" => 1,
                "Method" => "Whole-element Lipschitz ownership on a frozen surface",
                "MeshSHA256" => output_sha256,
                "ElementCertificateSHA256" => file_sha256(certificate),
                "ParentMeshSHA256" => file_sha256(input),
                "SignatureSHA256" => file_sha256(joinpath(root, "mesh-signature.csv")),
                "BoundarySHA256" => file_sha256(joinpath(root, "plan-view-boundary.csv")),
                "Radius" => radius,
                "Fabricated" => kind == "fabricated"
            ); sorted = true)
        end
        gmsh.clear()
        gmsh.open(output)
        canonical_node_geometry == canonical_nodes() || error("Serialized node geometry changed")
        canonical_volume_geometry == canonical_elements(3) || error("Serialized volume geometry/connectivity changed")
        canonical_surface_geometry == canonical_elements(2) || error("Serialized surface geometry/connectivity changed")
        check_measures(before_measures, grouped_surface_measures())
        serialized_volume_tags = collect(keys(dimension_elements(3)))
        after_quality = minimum(gmsh.model.mesh.getElementQualities(serialized_volume_tags, "minSICN"))
        abs(after_quality - before_quality) <= 1e-14 || error("Serialized volume quality changed")
        metadata = Dict(
            "Version" => 2,
            "Method" => "Exact element-wise interface ownership after frozen-surface volume generation",
            "Kind" => kind,
            "Input" => input,
            "InputSHA256" => file_sha256(input),
            "OutputSHA256" => output_sha256,
            "SignatureSHA256" => file_sha256(joinpath(root, "mesh-signature.csv")),
            "BoundarySHA256" => file_sha256(joinpath(root, "plan-view-boundary.csv")),
            "ProcessSHA256" => file_sha256(process_path),
            "InMemoryNodeAndElementIDsPreserved" => true,
            "SerializedGeometryAndConnectivityPreserved" => true,
            "SerializedIDsMayBeRenumberedByMSH22" => true,
            "SerializedRoundTripVerified" => true,
            "MinimumSignedInverseCondition" => after_quality,
            "MaximumGroupedPhysicalFamilyAreaDifference" => measure_difference,
            "ExpectedInterfaceAttributes" => expected === nothing ? nothing : sort!(collect(expected)),
            "ActualInterfaceAttributes" => sort!(collect(actual))
        )
        open(output * ".relabel.toml", "w") do stream
            TOML.print(stream, metadata; sorted = true)
        end
    finally
        gmsh.isInitialized() != 0 && gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
