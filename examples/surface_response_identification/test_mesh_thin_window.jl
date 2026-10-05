# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests of the thin window mesher's physical groups (decision 321):
#   julia --project=../cpw3d_surface/window_validation test_mesh_thin_window.jl
# A tiny synthetic two-plane window (a terminal pad on both planes joined by a bump, a ground
# bump) is meshed at lc 4 / 40 um: the terminal bump's shells form `bump_<label>` (footprint
# discs + sidewalls inside the footprint column, none of them grounded), the ground bump's
# shells keep `bump_surface` with the chip's bump attribute; the ground-only variant has
# `bump_surface` alone; the manifest `Attributes` table lists exactly the groups the mesh
# carries (a plane that is one ground sheet has no `gap_<plane>`); a bump whose conductor is
# neither ground nor a terminal is refused.

using Test
using JSON
import Gmsh: gmsh

include(joinpath(@__DIR__, "mesh_thin_window.jl"))

const SQUARE = [[0.0, 0.0], [40.0, 0.0], [40.0, 40.0], [0.0, 40.0]]
const HOLE = [[10.0, 10.0], [30.0, 10.0], [30.0, 30.0], [10.0, 30.0]]
const PAD = [[12.0, 12.0], [28.0, 12.0], [28.0, 28.0], [12.0, 28.0]]
const GROUND_BUMP = [[2.0, 2.0], [5.0, 2.0], [5.0, 5.0], [2.0, 5.0]]
const PAD_BUMP = [[18.0, 18.0], [22.0, 18.0], [22.0, 22.0], [18.0, 22.0]]

function synthetic_window(; pad_bump=true, l2_full_ground=false)
    l1 = [
        Dict("Conductor" => "ground", "Outer" => SQUARE, "Holes" => [HOLE]),
        Dict("Conductor" => "pad", "Outer" => PAD)
    ]
    l2 =
        l2_full_ground ? [Dict("Conductor" => "ground", "Outer" => SQUARE)] :
        [
            Dict("Conductor" => "ground", "Outer" => SQUARE, "Holes" => [HOLE]),
            Dict("Conductor" => "pad", "Outer" => PAD)
        ]
    bumps = [Dict("Conductor" => "ground", "Footprint" => GROUND_BUMP)]
    pad_bump && push!(bumps, Dict("Conductor" => "pad", "Footprint" => PAD_BUMP))
    return Dict(
        "Version" => 1,
        "Name" => "synthetic_bumps",
        "Box" => Dict("X" => [0.0, 40.0], "Y" => [0.0, 40.0]),
        "Planes" => [
            Dict(
                "Name" => "L1",
                "SurfaceZ" => 0.0,
                "Facing" => "up",
                "SubstrateThickness" => 20.0,
                "Polygons" => l1
            ),
            Dict(
                "Name" => "L2",
                "SurfaceZ" => 4.8,
                "Facing" => "down",
                "SubstrateThickness" => 20.0,
                "Polygons" => l2
            )
        ],
        "Bumps" => bumps,
        "Vacuum" => Dict("Below" => 0.0, "Above" => 0.0),
        "Terminals" => ["pad"]
    )
end

const CHIP = Dict(
    "Name" => "synthetic",
    "Planes" => [
        Dict("Name" => "L1", "Attributes" => [95], "Gap" => [124]),
        Dict("Name" => "L2", "Attributes" => [106], "Gap" => [26])
    ],
    "Bump" => [131],
    "Exterior" => [132],
    "Volumes" => Dict("Substrate" => [2, 1], "Vacuum" => 3),
    "TerminalAttributeBase" => 1000
)

function mesh_synthetic(directory, name, window)
    window_path = joinpath(directory, "$(name).json")
    chip_path = joinpath(directory, "chip.json")
    output = joinpath(directory, "$(name).msh2")
    open(window_path, "w") do stream
        return JSON.print(stream, window)
    end
    open(chip_path, "w") do stream
        return JSON.print(stream, CHIP)
    end
    return mesh_thin_window(window_path, chip_path, output, 4.0, 40.0), output
end

# Bounding box of every element of a physical group of a written mesh (the mesh re-read).
function group_bounds(output, dimension, attribute)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(output)
    lo = [Inf, Inf, Inf]
    hi = [-Inf, -Inf, -Inf]
    for entity in gmsh.model.getEntitiesForPhysicalGroup(dimension, attribute)
        _, _, node_tags = gmsh.model.mesh.getElements(dimension, entity)
        for tag in unique!(vcat(node_tags...))
            coordinates, _, _, _ = gmsh.model.mesh.getNode(tag)
            lo .= min.(lo, coordinates)
            hi .= max.(hi, coordinates)
        end
    end
    gmsh.finalize()
    return lo, hi
end

@testset "mesh_thin_window physical groups" begin
    @testset "terminal bump shells form bump_<label>" begin
        mktempdir() do directory
            manifest, output = mesh_synthetic(directory, "pad_bump", synthetic_window())
            attributes = manifest["Attributes"]
            groups = manifest["Groups"]
            @test Set(keys(attributes)) == Set(keys(groups))
            @test attributes["bump_surface"] == 131
            @test attributes["bump_pad"] == 1003
            @test attributes["pad_L1"] == 1001 && attributes["pad_L2"] == 1002
            @test groups["bump_pad"]["Elements"] > 0
            @test groups["bump_surface"]["Elements"] > 0
            # The terminal bump's shells lie in its footprint column, the ground bump's in its own.
            lo, hi = group_bounds(output, 2, 1003)
            @test all(lo .>= [18.0, 18.0, 0.0] .- 1e-9) &&
                  all(hi .<= [22.0, 22.0, 4.8] .+ 1e-9)
            @test hi[3] - lo[3] > 4.0  # both footprint discs and the sidewalls
            lo, hi = group_bounds(output, 2, 131)
            @test all(lo .>= [2.0, 2.0, 0.0] .- 1e-9) && all(hi .<= [5.0, 5.0, 4.8] .+ 1e-9)
            @test hi[3] - lo[3] > 4.0
            # The pad sheets no longer carry the footprint disc faces, but keep the rest of the pad.
            lo, hi = group_bounds(output, 2, 1001)
            @test all(lo .>= [12.0, 12.0, 0.0] .- 1e-9) &&
                  all(hi .<= [28.0, 28.0, 0.0] .+ 1e-9)
        end
    end

    @testset "ground-only bumps keep the single bump_surface group" begin
        mktempdir() do directory
            manifest, _ = mesh_synthetic(
                directory,
                "ground_bumps",
                synthetic_window(; pad_bump=false)
            )
            attributes = manifest["Attributes"]
            @test Set(keys(attributes)) == Set(keys(manifest["Groups"]))
            @test attributes["bump_surface"] == 131
            @test filter(startswith("bump_"), collect(keys(attributes))) == ["bump_surface"]
            @test haskey(attributes, "gap_L1") && haskey(attributes, "gap_L2")
        end
    end

    @testset "a plane without a gap sheet has no gap attribute" begin
        mktempdir() do directory
            manifest, _ = mesh_synthetic(
                directory,
                "full_l2",
                synthetic_window(; pad_bump=false, l2_full_ground=true)
            )
            attributes = manifest["Attributes"]
            @test Set(keys(attributes)) == Set(keys(manifest["Groups"]))
            @test haskey(attributes, "gap_L1")
            @test !haskey(attributes, "gap_L2")
            @test !haskey(manifest["Groups"], "gap_L2")
        end
    end

    @testset "a bump of an unknown conductor is refused" begin
        mktempdir() do directory
            window = synthetic_window(; pad_bump=false)
            push!(window["Bumps"], Dict("Conductor" => "stray", "Footprint" => PAD_BUMP))
            @test_throws ErrorException mesh_synthetic(directory, "stray", window)
            return gmsh.isInitialized() == 1 && gmsh.finalize()
        end
    end
end
