# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests of the thin window mesher's physical groups (decision 321) and of its
# ring-orientation contract + gates (decision 1022, census finding F-1):
#   julia --project=../cpw3d_surface/window_validation test_mesh_thin_window.jl
# A tiny synthetic two-plane window (a terminal pad on both planes joined by a bump, a ground
# bump) is meshed at lc 4 / 40 um: the terminal bump's shells form `bump_<label>` (footprint
# discs + sidewalls inside the footprint column, none of them grounded), the ground bump's
# shells keep `bump_surface` with the chip's bump attribute; the ground-only variant has
# `bump_surface` alone; the manifest `Attributes` table lists exactly the groups the mesh
# carries (a plane that is one ground sheet has no `gap_<plane>`); a bump whose conductor is
# neither ground nor a terminal is refused.
# The F-1 reproducer (a single-plane ground square with ONE clockwise hole holding an
# island; the mesher of record left the ground as a free sheet: 2,012 orphan + 330 duplicate
# triangles, rc 0): the hole is reoriented (manifest `ReorientedLoops`), every triangle of
# the written mesh is the face of a tetrahedron, and with the normalisation disabled the
# pre-fragment face gate refuses the invalid face (area outer PLUS hole); a hole outside its
# outer (a malformed face) and two coincident metal sheets are refused too; a metal sheet
# reaching beyond the box wall (its image outside the box is adjacent to no volume) is
# refused by the post-fragment conformity gate.

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

    @testset "co-oriented record inputs reorient nothing" begin
        mktempdir() do directory
            manifest, _ =
                mesh_synthetic(directory, "cooriented", synthetic_window(; pad_bump=false))
            @test manifest["ReorientedLoops"] == Dict("Outers" => 0, "Holes" => 0)
            @test manifest["ToolImagesNotShared"] == 0
        end
    end
end

# The F-1 reproducer (decision 1022): a 1400-um ground square whose single 1300-um hole is
# given CLOCKWISE against the counter-clockwise outer and holds a 600-um island (single
# plane, LC 40 / 400).
const RING_OUTER = [[300.0, 300.0], [1700.0, 300.0], [1700.0, 1700.0], [300.0, 1700.0]]
const RING_HOLE_CW = [[350.0, 1650.0], [1650.0, 1650.0], [1650.0, 350.0], [350.0, 350.0]]
const RING_ISLAND = [[700.0, 700.0], [1300.0, 700.0], [1300.0, 1300.0], [700.0, 1300.0]]

function ring_island_window(; hole=RING_HOLE_CW)
    return Dict(
        "Version" => 1,
        "Name" => "ring_island",
        "Box" => Dict("X" => [0.0, 2000.0], "Y" => [0.0, 2000.0]),
        "Planes" => [
            Dict(
                "Name" => "L1",
                "SurfaceZ" => 0.0,
                "Facing" => "up",
                "SubstrateThickness" => 525.0,
                "Polygons" => [
                    Dict("Conductor" => "ground", "Outer" => RING_OUTER, "Holes" => [hole]),
                    Dict("Conductor" => "island", "Outer" => RING_ISLAND)
                ]
            )
        ],
        "Bumps" => [],
        "Vacuum" => Dict("Below" => 0.0, "Above" => 200.0),
        "Terminals" => ["island"]
    )
end

const SINGLE_PLANE_CHIP = Dict(
    "Name" => "single_plane",
    "Planes" => [Dict("Name" => "L1", "Attributes" => [4], "Gap" => [8])],
    "Bump" => [131],
    "Exterior" => [3],
    "Volumes" => Dict("Substrate" => [1], "Vacuum" => 2),
    "TerminalAttributeBase" => 1000
)

function mesh_ring_island(directory, name, window; normalise_loops=true)
    window_path = joinpath(directory, "$(name).json")
    chip_path = joinpath(directory, "chip.json")
    output = joinpath(directory, "$(name).msh2")
    open(stream -> JSON.print(stream, window), window_path, "w")
    open(stream -> JSON.print(stream, SINGLE_PLANE_CHIP), chip_path, "w")
    manifest =
        mesh_thin_window(window_path, chip_path, output, 40.0, 400.0; normalise_loops)
    return manifest, output
end

# The triangles of a written mesh that are the face of no tetrahedron, and the duplicated
# ones (the census conformity gate: both counts must be 0 for MFEM to load the mesh).
function nonconforming_triangles(output)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(output)
    tet_faces = Set{NTuple{3, Int}}()
    for (dim, tag) in gmsh.model.getEntities(3)
        types, _, node_tags = gmsh.model.mesh.getElements(dim, tag)
        for (t, nodes) in zip(types, node_tags)
            t == 4 || error("element type $(t) is not a linear tetrahedron")
            for e = 0:(length(nodes) ÷ 4 - 1)
                a, b, c, d = Int.(nodes[(4e + 1):(4e + 4)])
                for face in ((a, b, c), (a, b, d), (a, c, d), (b, c, d))
                    push!(tet_faces, Tuple(sort!(collect(face))))
                end
            end
        end
    end
    orphans = 0
    seen = Set{NTuple{3, Int}}()
    duplicates = 0
    for (dim, tag) in gmsh.model.getEntities(2)
        types, _, node_tags = gmsh.model.mesh.getElements(dim, tag)
        for (t, nodes) in zip(types, node_tags)
            t == 2 || error("element type $(t) is not a linear triangle")
            for e = 0:(length(nodes) ÷ 3 - 1)
                face = Tuple(sort!(Int.(nodes[(3e + 1):(3e + 3)])))
                face in tet_faces || (orphans += 1)
                face in seen ? (duplicates += 1) : push!(seen, face)
            end
        end
    end
    gmsh.finalize()
    return orphans, duplicates
end

@testset "mesh_thin_window ring orientation and conformity gate (decision 1022)" begin
    @testset "a clockwise hole holding an island is reoriented and meshes conforming" begin
        mktempdir() do directory
            manifest, output = mesh_ring_island(directory, "ring_cw", ring_island_window())
            @test manifest["ReorientedLoops"] == Dict("Outers" => 0, "Holes" => 1)
            @test manifest["ToolImagesNotShared"] == 0
            @test manifest["ToolImagesShared"] == 0
            @test manifest["Groups"]["ground_L1"]["Elements"] > 0
            @test manifest["Groups"]["island_L1"]["Elements"] > 0
            @test manifest["Groups"]["gap_L1"]["Elements"] > 0
            @test nonconforming_triangles(output) == (0, 0)
        end
    end

    @testset "the same hole given counter-clockwise reorients nothing" begin
        mktempdir() do directory
            manifest, output = mesh_ring_island(
                directory,
                "ring_ccw",
                ring_island_window(; hole=reverse(RING_HOLE_CW))
            )
            @test manifest["ReorientedLoops"] == Dict("Outers" => 0, "Holes" => 0)
            @test nonconforming_triangles(output) == (0, 0)
        end
    end

    @testset "a clockwise outer is reoriented too" begin
        mktempdir() do directory
            window = ring_island_window()
            window["Planes"][1]["Polygons"][2]["Outer"] = reverse(RING_ISLAND)
            manifest, output = mesh_ring_island(directory, "ring_cw_outer", window)
            @test manifest["ReorientedLoops"] == Dict("Outers" => 1, "Holes" => 1)
            @test nonconforming_triangles(output) == (0, 0)
        end
    end

    @testset "a hole outside its outer is a malformed face: refused" begin
        mktempdir() do directory
            # the island's single "hole" lies in the ground metal, as in the census CHR-012
            # set
            window = ring_island_window()
            window["Planes"][1]["Polygons"][2]["Holes"] =
                [[[400.0, 400.0], [500.0, 400.0], [500.0, 500.0], [400.0, 500.0]]]
            err = try
                mesh_ring_island(directory, "hole_outside", window)
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test occursin("malformed metal face for island", err.msg)
            @test !isfile(joinpath(directory, "hole_outside.msh2"))
            return gmsh.isInitialized() == 1 && gmsh.finalize()
        end
    end

    @testset "coincident metal sheets are refused after the fragment" begin
        mktempdir() do directory
            window = ring_island_window()
            push!(
                window["Planes"][1]["Polygons"],
                Dict("Conductor" => "twin", "Outer" => RING_ISLAND)
            )
            push!(window["Terminals"], "twin")
            err = try
                mesh_ring_island(directory, "twin", window)
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException && occursin("more than one metal sheet", err.msg)
            @test !isfile(joinpath(directory, "twin.msh2"))
            return gmsh.isInitialized() == 1 && gmsh.finalize()
        end
    end

    @testset "without the normalisation the face gate refuses the invalid face" begin
        mktempdir() do directory
            # gmsh reverses the hole wire itself: the clockwise hole as given makes a face
            # whose area is the outer PLUS the hole, caught before the fragment by
            # check_face_area
            err = try
                mesh_ring_island(
                    directory,
                    "ring_unnormalised",
                    ring_island_window();
                    normalise_loops=false
                )
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test occursin("malformed metal face for ground", err.msg)
            @test !isfile(joinpath(directory, "ring_unnormalised.msh2"))
            return gmsh.isInitialized() == 1 && gmsh.finalize()
        end
    end

    @testset "a metal sheet beyond the box wall: the conformity gate refuses" begin
        mktempdir() do directory
            # the island straddles the wall x = 2000 (a polygon not clipped to the box): its
            # face has the polygon's net area (the face gate passes), the fragment leaves
            # the part outside the box as an image adjacent to no volume
            window = ring_island_window(; hole=reverse(RING_HOLE_CW))
            window["Planes"][1]["Polygons"][2]["Outer"] =
                [[1750.0, 700.0], [2250.0, 700.0], [2250.0, 1300.0], [1750.0, 1300.0]]
            err = try
                mesh_ring_island(directory, "beyond_wall", window)
                nothing
            catch caught
                caught
            end
            @test err isa ErrorException
            @test occursin("conformity gate", err.msg)
            @test occursin("island (plane 1)", err.msg)
            @test occursin("0 adjacent volume(s)", err.msg)
            @test !isfile(joinpath(directory, "beyond_wall.msh2"))
            return gmsh.isInitialized() == 1 && gmsh.finalize()
        end
    end
end
