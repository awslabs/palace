# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests of the graded cross-section sweep (graded_sweep.jl; Block G milestone 1):
#   julia --project=. test_graded_sweep.jl
# The nested stacks and row ranges (the Q2 probe's values on the strip window), the 2D
# cross-section (area closure, every node used), the graded sweep of synthetic one-plane
# windows (a straight strip, a CPW with a terminal, a slot with a backside) against the tensor
# sweep: areas / volumes equal to 1e-9, boundary closure and tag consistency in process, every
# node used, the maximum edge ratio <= t / r (+0.1 %), deterministic bytes, alpha / beta
# economics, the refusals (Gmsh band, two planes), and the default tensor output unchanged
# (the recorded probe sha of the strip window). WINDOW_MESH_LONG_TESTS=1 adds the external
# validator and census runs (Julia subprocesses, ~10 s each).

using Test
using JSON
using SHA

include(joinpath(@__DIR__, "PolygonWindowMesh.jl"))
using .PolygonWindowMesh
include(joinpath(@__DIR__, "synthetic_two_level_window.jl"))

const PWM = PolygonWindowMesh
const LONG_TESTS = get(ENV, "WINDOW_MESH_LONG_TESTS", "0") == "1"

function one_plane_window(name, box_x, box_y, polygons; substrate, below=0.0, above)
    return read_polygon_set(
        Dict(
            "Version" => 1,
            "Name" => name,
            "Box" => Dict("X" => box_x, "Y" => box_y),
            "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
            "Planes" => [
                Dict(
                    "Name" => "L1",
                    "SurfaceZ" => 0.0,
                    "Facing" => "up",
                    "SubstrateThickness" => substrate,
                    "Polygons" => polygons
                )
            ],
            "Vacuum" => Dict("Below" => below, "Above" => above)
        )
    )
end

# The Q2 probe's strip window (reference-quality-20261002/mesher-scoping/probe/strip_window.json).
strip_window() = one_plane_window(
    "strip_probe",
    [0.0, 50.0],
    [0.0, 60.0],
    [Dict("Conductor" => "ground", "Outer" => [[0, 20], [50, 20], [50, 40], [0, 40]])];
    substrate=525.0,
    above=525.0
)
# A CPW cut by the window on both sides: two grounds and a terminal trace (open chains with
# wall ends, two conductors).
cpw_window() = one_plane_window(
    "cpw",
    [0.0, 30.0],
    [0.0, 60.0],
    [
        Dict("Conductor" => "ground", "Outer" => [[0, 0], [30, 0], [30, 24], [0, 24]]),
        Dict("Conductor" => "trace", "Outer" => [[0, 28], [30, 28], [30, 32], [0, 32]]),
        Dict("Conductor" => "ground", "Outer" => [[0, 36], [30, 36], [30, 60], [0, 60]])
    ];
    substrate=100.0,
    above=100.0
)
# A rectangular slot in a ground (a closed hole chain and a closed outer chain) with vacuum
# below the substrate (a `substrate_backside` group).
slot_window() = one_plane_window(
    "slot",
    [0.0, 20.0],
    [0.0, 40.0],
    [
        Dict(
            "Conductor" => "ground",
            "Outer" => [[0, 0], [20, 0], [20, 40], [0, 40]],
            "Holes" => [[[5, 10], [15, 10], [15, 30], [5, 30]]]
        )
    ];
    substrate=50.0,
    below=20.0,
    above=50.0
)

# ASCII MSH2 reader for the in-process checks.
function read_msh2(path)
    lines = readlines(path)
    i = findfirst(==("\$Nodes"), lines)
    n_nodes = parse(Int, lines[i + 1])
    nodes = Vector{NTuple{3, Float64}}(undef, n_nodes)
    for k = 1:n_nodes
        fields = split(lines[i + 1 + k])
        nodes[k] = (
            parse(Float64, fields[2]),
            parse(Float64, fields[3]),
            parse(Float64, fields[4])
        )
    end
    i = findfirst(==("\$Elements"), lines)
    n_elements = parse(Int, lines[i + 1])
    triangles = Tuple{Int, NTuple{3, Int}}[]
    tetrahedra = Tuple{Int, NTuple{4, Int}}[]
    for k = 1:n_elements
        fields = parse.(Int, split(lines[i + 1 + k]))
        if fields[2] == 2
            push!(triangles, (fields[4], (fields[6], fields[7], fields[8])))
        elseif fields[2] == 4
            push!(tetrahedra, (fields[4], (fields[6], fields[7], fields[8], fields[9])))
        else
            error("Unexpected element type $(fields[2])")
        end
    end
    return nodes, triangles, tetrahedra
end

sorted3(a, b, c) = Tuple(sort([a, b, c]))

# Boundary closure and tag consistency, every node used, the maximum edge-length ratio
# overall and off the metal / trench slab (tetrahedra not inside z in [-0.05, 0.1]: the
# interfaces are kept in every stack, so a 30-um region triangle over the 0.05-um trench is
# an inherent 600 in both sweeps; everything else must stay within the designed t / r).
function mesh_checks(path; slab=(-0.05, 0.1))
    nodes, triangles, tetrahedra = read_msh2(path)
    owners = Dict{NTuple{3, Int}, Vector{Int}}()
    used = falses(length(nodes))
    max_ratio, max_ratio_off_slab = 0.0, 0.0
    for (attribute, t) in tetrahedra
        for n in t
            used[n] = true
        end
        for (a, b, c) in ((1, 2, 3), (1, 2, 4), (1, 3, 4), (2, 3, 4))
            push!(get!(owners, sorted3(t[a], t[b], t[c]), Int[]), attribute)
        end
        edges = [
            hypot((nodes[t[a]] .- nodes[t[b]])...) for
            (a, b) in ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
        ]
        ratio = maximum(edges) / minimum(edges)
        max_ratio = max(max_ratio, ratio)
        z = [nodes[n][3] for n in t]
        in_slab = minimum(z) >= slab[1] - 1.0e-9 && maximum(z) <= slab[2] + 1.0e-9
        in_slab || (max_ratio_off_slab = max(max_ratio_off_slab, ratio))
    end
    tagged = Dict(sorted3(t...) => attribute for (attribute, t) in triangles)
    length(tagged) == length(triangles) || error("A surface face is tagged twice")
    untagged = count(length(o) == 1 && !haskey(tagged, f) for (f, o) in owners)
    inconsistent = 0
    for (face, attribute) in tagged
        o = get(owners, face, Int[])
        ok =
            attribute in (6, 9) ? sort(o) == [1, 2] :
            attribute == 3 ? length(o) == 1 : attribute in (5, 8, 11) ? o == [1] : o == [2] # metal-substrate / metal-air
        ok || (inconsistent += 1)
    end
    return (
        nodes=length(nodes),
        tetrahedra=length(tetrahedra),
        untagged=untagged,
        inconsistent=inconsistent,
        unused_nodes=count(!, used),
        max_edge_ratio=max_ratio,
        max_edge_ratio_off_slab=max_ratio_off_slab
    )
end

# JIT warm-up outside the test sets (the first mesher call compiles for ~10 s): a tiny strip in
# both sweeps and the in-process reader, so every test set below measures its own work.
let directory = mktempdir()
    tiny = one_plane_window(
        "warmup",
        [0.0, 10.0],
        [0.0, 20.0],
        [Dict("Conductor" => "ground", "Outer" => [[0, 6], [10, 6], [10, 14], [0, 14]])];
        substrate=20.0,
        above=20.0
    )
    mesh_polygon_window(tiny, 0.05, 5.0, joinpath(directory, "tensor.msh2"); verbose=false)
    mesh_polygon_window(
        tiny,
        0.05,
        5.0,
        joinpath(directory, "graded.msh2");
        verbose=false,
        sweep=:graded
    )
    mesh_checks(joinpath(directory, "graded.msh2"))
end

@testset "graded stacks and row ranges (the Q2 probe's strip values)" begin
    spec = strip_window()
    stack = PWM.z_levels(spec, 10, 5)
    zs = stack.levels
    @test length(zs) == 40
    interfaces = [-0.05, 0.0, 0.1, stack.z_bottom, stack.z_top]
    heights = PWM.band_heights(0.01, 2.0, 7)
    stacks = PWM.graded_stacks(zs, interfaces, heights, 1.0, 3.0, (0.0, 0.1), 30.0)
    @test stacks.rows == 7
    @test [length(l) for l in stacks.levels] == [40, 40, 32, 28, 25, 23, 22, 20, 18, 16, 14, 12, 10]
    @test stacks.spacings ≈
          [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64, 1.28, 2.56, 5.12, 10.24, 20.48]
    # Nested, interfaces kept everywhere, the spacing rule between kept non-interface levels.
    for s = 2:length(stacks.levels)
        levels = stacks.levels[s]
        @test issubset(levels, stacks.levels[s - 1])
        @test all(any(i -> abs(zs[i] - w) <= 1.0e-9, levels) for w in interfaces)
        for (a, b) in zip(levels[1:(end - 1)], levels[2:end])
            PWM.is_interface(zs[b], interfaces) && continue
            @test zs[b] - zs[a] >= stacks.spacings[s - 1] - 1.0e-9
        end
    end
    ranges = [(zs[lo], zs[hi]) for (lo, hi) in stacks.ranges]
    probe_ranges =
        [(-0.03, 0.15), (-0.1, 0.2), (-0.5, 0.5), (-1.0, 1.0), (-2.0, 2.0), (-5.0, 5.0)]
    @test all(all(isapprox.(ranges[k], probe_ranges[k]; atol=1.0e-9)) for k = 1:6)
    @test ranges[7] == (-525.0, 525.0)
    # Strictly nested; every range end is a level of the next row's stack.
    for k = 2:7
        @test ranges[k][1] < ranges[k - 1][1] && ranges[k][2] > ranges[k - 1][2]
    end
    for k = 1:6
        @test stacks.ranges[k][1] in stacks.levels[k + 2] &&
              stacks.ranges[k][2] in stacks.levels[k + 2]
    end
    # alpha 2: coarser stacks (nested thinning at twice the spacing); beta 2 at the same
    # alpha: rows end at the same or a nearer level.
    coarser = PWM.graded_stacks(zs, interfaces, heights, 2.0, 3.0, (0.0, 0.1), 30.0)
    @test all(length.(coarser.levels) .<= length.(stacks.levels))
    @test sum(length.(coarser.levels)) < sum(length.(stacks.levels))
    shorter = PWM.graded_stacks(zs, interfaces, heights, 1.0, 2.0, (0.0, 0.1), 30.0)
    @test shorter.levels == stacks.levels
    @test all(
        zs[shorter.ranges[k][1]] >= ranges[k][1] &&
        zs[shorter.ranges[k][2]] <= ranges[k][2] for k = 1:6
    )
    @test any(shorter.ranges[k] != stacks.ranges[k] for k = 1:6)
    @test_throws ErrorException PWM.graded_stacks(
        zs,
        interfaces,
        heights,
        0.0,
        3.0,
        (0.0, 0.1),
        30.0
    )
    section = PWM.cross_section(stacks, zs, heights)
    @test length(section.nodes) == 128 && length(section.triangles) == 194
    @test sum(n for (_, n) in section.pair_triangles) == 194
    @test length(section.pair_triangles) == 1 + 2 * 6 + 6
end

function compare_sweeps(spec, radial_um, tangential_um, directory; kwargs...)
    tensor = mesh_polygon_window(
        spec,
        radial_um,
        tangential_um,
        joinpath(directory, "tensor.msh2");
        verbose=false
    )
    graded = mesh_polygon_window(
        spec,
        radial_um,
        tangential_um,
        joinpath(directory, "graded.msh2");
        verbose=false,
        sweep=:graded,
        kwargs...
    )
    @test tensor["sweep"] == "tensor" && graded["sweep"] == "graded"
    @test graded["nonpositive"] == 0 && tensor["nonpositive"] == 0
    @test graded["tetrahedra"] < 0.5 * tensor["tetrahedra"]
    @test keys(graded["surface_area_um2"]) == keys(tensor["surface_area_um2"])
    for (attribute, area) in tensor["surface_area_um2"]
        @test graded["surface_area_um2"][attribute] ≈ area rtol = 1.0e-9
    end
    for (attribute, volume) in tensor["volume_um3"]
        @test graded["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
    end
    # Metal faces are the same plan quads / triangles in both sweeps.
    for attribute in ("4", "5", "7", "8")
        haskey(tensor["surface_attribute_counts"], attribute) || continue
        @test graded["surface_attribute_counts"][attribute] ==
              tensor["surface_attribute_counts"][attribute]
    end
    checks = mesh_checks(joinpath(directory, "graded.msh2"))
    @test checks.nodes == graded["nodes"] && checks.tetrahedra == graded["tetrahedra"]
    @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
    # The designed anisotropy t / r off the slab: + 2 for the mitre column of a convex corner
    # (its row-k node lies h_k < 2 w_k farther along the edge, the same plan edge as in the
    # tensor sweep) and + 0.1 % for the prism diagonals; the slab region prisms are bounded by
    # the region size over the overetch (30 / 0.05 = 600).
    @test checks.max_edge_ratio_off_slab <=
          (tangential_um / radial_um + 2.0) * (1.0 + 1.0e-3)
    @test checks.max_edge_ratio <=
          max(tangential_um / radial_um, PWM.REGION_MESH_SIZE_MAX_UM / spec.overetch) *
          (1.0 + 1.0e-3)
    record = graded["graded_sweep"]
    @test record["alpha"] == get(kwargs, :alpha, 1.0) &&
          record["beta"] == get(kwargs, :beta, 3.0)
    @test record["band_tetrahedra"] + record["region_tetrahedra"] == graded["tetrahedra"]
    @test record["band_tetrahedra"] == 3 * record["band_prisms"]
    return tensor, graded, checks
end

@testset "strip window: graded vs tensor at the probe's r10 / t5" begin
    directory = mktempdir()
    # Recorded: tensor 63,768 tets (the probe's baseline); the probe's graded alpha 1 / beta 3
    # had 32,460 with a structured interior finer than Gmsh's region (kept here).
    tensor, graded, checks = compare_sweeps(strip_window(), 0.01, 5.0, directory)
    @test tensor["tetrahedra"] == 63768 && graded["tetrahedra"] == 23216
    @test graded["graded_sweep"]["chains"] == 4 && graded["graded_sweep"]["columns"] == 44
    @test graded["graded_sweep"]["region_prisms_hanging"] > 0
    @test checks.max_edge_ratio ≈ 500.0 rtol = 1.0e-3
end

@testset "strip window: deterministic bytes and alpha / beta" begin
    directory = mktempdir()
    first = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "graded.msh2");
        verbose=false,
        sweep=:graded
    )
    again = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "graded_again.msh2");
        verbose=false,
        sweep=:graded
    )
    @test again["sha256"] == first["sha256"] && again["bytes"] == first["bytes"]
    @test first["sha256"] == bytes2hex(open(sha256, joinpath(directory, "graded.msh2")))
    # alpha 2 / beta 2: coarser stacks, shorter rows, fewer tetrahedra, the same geometry.
    _, cheaper, _ = compare_sweeps(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "a2b2");
        alpha=2.0,
        beta=2.0
    )
    @test cheaper["tetrahedra"] < first["tetrahedra"]
    @test all(
        cheaper["graded_sweep"]["stack_sizes"] .<= first["graded_sweep"]["stack_sizes"]
    )
end

@testset "CPW window with a terminal (open chains, wall ends)" begin
    _, cpw, _ = compare_sweeps(cpw_window(), 0.02, 5.0, mktempdir())
    @test cpw["graded_sweep"]["chains"] == 8
    @test haskey(cpw["surface_area_um2"], "7") && haskey(cpw["surface_area_um2"], "8")
    @test cpw["surface_area_um2"]["8"] ≈ 30.0 * 4.0 rtol = 1.0e-9
end

@testset "slot window (closed chains, substrate backside)" begin
    _, slot, _ = compare_sweeps(slot_window(), 0.02, 5.0, mktempdir())
    @test slot["graded_sweep"]["chains"] == 2
    @test haskey(slot["surface_area_um2"], "9")
    @test slot["surface_area_um2"]["9"] ≈ 800.0 rtol = 1.0e-9
end

@testset "graded sweep refusals and the default tensor output" begin
    directory = mktempdir()
    @test_throws ErrorException mesh_polygon_window(
        strip_window(),
        0.05,
        5.0,
        joinpath(directory, "gmsh.msh2");
        verbose=false,
        sweep=:graded,
        band_mode=:gmsh
    )
    @test_throws ErrorException mesh_polygon_window(
        strip_window(),
        0.05,
        5.0,
        joinpath(directory, "bad.msh2");
        verbose=false,
        sweep=:slanted
    )
    message = try
        mesh_polygon_window(
            read_polygon_set(synthetic_two_level_window()),
            0.2,
            5.0,
            joinpath(directory, "two.msh2");
            verbose=false,
            sweep=:graded
        )
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("two planes", message) && occursin("M3", message)
    # The default sweep reproduces the recorded tensor mesh of the probe window bit for bit
    # (reference-quality-20261002/mesher-scoping/probe/strip_tensor_r10_t5.json, generated on
    # macOS; Gmsh's region triangulation may differ on another platform, so the sha is pinned
    # there only and the counts everywhere).
    tensor = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "strip_tensor.msh2");
        verbose=false
    )
    @test tensor["sweep"] == "tensor" && !haskey(tensor, "graded_sweep")
    @test tensor["nodes"] == 12418 && tensor["tetrahedra"] == 63768
    @test tensor["bytes"] == 2635089
    if Sys.isapple()
        @test tensor["sha256"] ==
              "8f234e5b0933f90aaef640399f991e73f950b5fb102766f840565d9b934980f6"
    end
end

if LONG_TESTS
    @testset "validator and census on the strip (long)" begin
        directory = mktempdir()
        meshes = Dict(
            "tensor" => mesh_polygon_window(
                strip_window(),
                0.01,
                5.0,
                joinpath(directory, "tensor.msh2");
                verbose=false
            ),
            "graded" => mesh_polygon_window(
                strip_window(),
                0.01,
                5.0,
                joinpath(directory, "graded.msh2");
                verbose=false,
                sweep=:graded
            )
        )
        julia = `$(Base.julia_cmd()) --project=$(@__DIR__)`
        census = Dict{String, Any}()
        for (name, manifest) in meshes
            mesh = joinpath(directory, "$name.msh2")
            run(pipeline(`$julia validate_window_mesh.jl $mesh`; stdout=devnull))
            validation = JSON.parsefile(joinpath(directory, "$name.validation.json"))
            @test all(iszero, values(validation["surface_adjacency_errors"]))
            @test validation["untagged_boundary_faces"] == 0
            @test validation["mesh_quality"]["nonpositive"] == 0
            @test validation["nodes"] == manifest["nodes"]
            run(pipeline(`$julia reference_mesh_census.jl $mesh`; stdout=devnull))
            census[name] = JSON.parsefile(joinpath(directory, "$name.census.json"))
        end
        overall(name) = only(filter(row -> row["label"] == "ALL", census[name]["rows"]))
        all_tensor, all_graded = overall("tensor"), overall("graded")
        @test all_graded["edge_ratio_max"] <= 500.0 * (1.0 + 1.0e-3)
        @test all_tensor["edge_ratio_max"] > 1.0e4
        @test all_graded["sicn_min"] > 10 * all_tensor["sicn_min"]
    end
end
