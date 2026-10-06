# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tests of the graded cross-section sweep (graded_sweep.jl; Block G milestones 1-3):
#   julia --project=. test_graded_sweep.jl
# The nested stacks and row ranges (the Q2 probe's strip values with the fabricated step's
# faces; two planes: the gap midpoint kept, one plane's rows capped there), the 2D
# cross-section (area closure, every node used), the graded sweep of synthetic windows against
# the tensor sweep — one plane: a straight strip, a CPW with a terminal, a slot with a backside,
# a facing-down strip, a metal spike (a fan, clearance-capped columns, an inward mitre scaled
# 4x, collapse) and a 1-um channel (facing-capped columns, collapse in both directions, a
# region strip between capped tops); two planes: the synthetic two-level window (crossing
# edges, a bump), a flip-chip window with pure-L1, pure-L2 and mixed chains, facing tops of
# both planes 2 um apart and a bump, and a coincident cross-plane run: areas / volumes equal
# to 1e-9, metal face counts equal (same plan), boundary closure and tag consistency in
# process, every node used, the maximum edge ratio <= t / r (+ the scaled columns' reach,
# + 0.1 %) off the step slabs, deterministic bytes, alpha / beta economics, the M3 defaults
# (region size grading of the graded plan, the region ring at Z_K; the tensor plan unchanged),
# the refusals (Gmsh band, a gap too small for beta, a box too short for the rows), and the
# default tensor output unchanged (the recorded probe sha of the strip window), and the CLI
# (mesh_polygon_window.jl --plan-only on the strip in both sweeps, a Julia subprocess each)
# against the in-process plan manifest.
# WINDOW_MESH_LONG_TESTS=1 adds the external validator and census runs (Julia subprocesses,
# ~10 s each).

using Test
using JSON
using SHA

include(joinpath(@__DIR__, "PolygonWindowMesh.jl"))
using .PolygonWindowMesh
include(joinpath(@__DIR__, "synthetic_two_level_window.jl"))

const PWM = PolygonWindowMesh
const LONG_TESTS = get(ENV, "WINDOW_MESH_LONG_TESTS", "0") == "1"

# The polygon-set JSON document (SCHEMA.md) of a one-plane window; `one_plane_window` reads
# it in process, the CLI smoke test writes it to a file.
function one_plane_window_data(
    name,
    box_x,
    box_y,
    polygons;
    substrate,
    below=0.0,
    above,
    facing="up"
)
    return Dict(
        "Version" => 1,
        "Name" => name,
        "Box" => Dict("X" => box_x, "Y" => box_y),
        "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
        "Planes" => [
            Dict(
                "Name" => "L1",
                "SurfaceZ" => 0.0,
                "Facing" => facing,
                "SubstrateThickness" => substrate,
                "Polygons" => polygons
            )
        ],
        "Vacuum" => Dict("Below" => below, "Above" => above)
    )
end
one_plane_window(name, box_x, box_y, polygons; kwargs...) =
    read_polygon_set(one_plane_window_data(name, box_x, box_y, polygons; kwargs...))

# The Q2 probe's strip window (reference-quality-20261002/mesher-scoping/probe/strip_window.json).
strip_window_data() = one_plane_window_data(
    "strip_probe",
    [0.0, 50.0],
    [0.0, 60.0],
    [Dict("Conductor" => "ground", "Outer" => [[0, 20], [50, 20], [50, 40], [0, 40]])];
    substrate=525.0,
    above=525.0
)
strip_window() = read_polygon_set(strip_window_data())
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
# The strip facing down (the metal below the surface, the substrate above, a backside against
# the vacuum above).
strip_down_window() = one_plane_window(
    "strip_down",
    [0.0, 50.0],
    [0.0, 60.0],
    [Dict("Conductor" => "ground", "Outer" => [[0, 20], [50, 20], [50, 40], [0, 40]])];
    substrate=100.0,
    below=100.0,
    above=20.0,
    facing="down"
)
# A ground with a 28-degree spike: for the vacuum a 152-degree outward turn (a fan of 7
# columns) and two inward mitres; for the metal an inward mitre scaled 1 / cos(76 deg) = 4.1
# with the clearance caps along both edges (4 rows at the tip, 5 at the base corners, 6
# elsewhere: collapse in both directions); four wall ends.
spike_window() = one_plane_window(
    "spike",
    [0.0, 30.0],
    [0.0, 30.0],
    [
        Dict(
            "Conductor" => "ground",
            "Outer" =>
                [[0, 0], [30, 0], [30, 10], [17, 10], [15, 18], [13, 10], [0, 10]]
        )
    ];
    substrate=50.0,
    above=50.0
)
# A hole with a 1-um channel: the facing rule caps the channel's columns at 4 rows (0.4 um)
# against 6 on the square, so the hole chain collapses at the channel mouth and the region in
# the channel has both corners on capped tops; eight inward / outward 90-degree mitres.
channel_window() = one_plane_window(
    "channel",
    [0.0, 30.0],
    [0.0, 30.0],
    [
        Dict(
            "Conductor" => "ground",
            "Outer" => [[0, 0], [30, 0], [30, 30], [0, 30]],
            "Holes" => [[
                [10, 10],
                [18, 10],
                [18, 18],
                [14.5, 18],
                [14.5, 26],
                [13.5, 26],
                [13.5, 18],
                [10, 18]
            ]]
        )
    ];
    substrate=50.0,
    above=50.0
)

# A flip-chip window with chains of every cross-section variant: L1 (z 0, up) and L2 (z 4.8,
# down) grounds with nested rectangular holes whose edges run 2 um apart in plan (facing tops
# of the two planes with a region strip between them: capped columns on both, region triangles
# with three top corners of two variants), an L1 trace inside both holes (a pure L1 chain,
# capped at the gap midpoint), an L2 trace inside both holes (pure L2), and a ground bump on the
# common ground frame (a bump chain: every row active across the gap).
function two_plane_window(; gap_z=4.8, beta_probe=false)
    return read_polygon_set(
        Dict(
            "Version" => 1,
            "Name" => "two_plane",
            "MatchingRadius" => 1.9,
            "Box" => Dict("X" => [0.0, 80.0], "Y" => [0.0, 60.0]),
            "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
            "Planes" => [
                Dict(
                    "Name" => "L1",
                    "SurfaceZ" => 0.0,
                    "Facing" => "up",
                    "SubstrateThickness" => 20.0,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0, 0], [80, 0], [80, 60], [0, 60]],
                            "Holes" => [[[14, 14], [66, 14], [66, 46], [14, 46]]]
                        ),
                        Dict(
                            "Conductor" => "trace_l1",
                            "Outer" => [[20, 28], [40, 28], [40, 32], [20, 32]]
                        )
                    ]
                ),
                Dict(
                    "Name" => "L2",
                    "SurfaceZ" => gap_z,
                    "Facing" => "down",
                    "SubstrateThickness" => 20.0,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0, 0], [80, 0], [80, 60], [0, 60]],
                            "Holes" => [[[12, 12], [68, 12], [68, 48], [12, 48]]]
                        ),
                        Dict(
                            "Conductor" => "trace_l2",
                            "Outer" => [[44, 28], [60, 28], [60, 32], [44, 32]]
                        )
                    ]
                )
            ],
            "Bumps" => [
                Dict(
                    "Conductor" => "ground",
                    "Footprint" => regular_polygon(6.0, 6.0, 4.0, 16)
                )
            ],
            "Vacuum" => Dict("Below" => 0.0, "Above" => 0.0),
            "Terminals" => ["trace_l1", "trace_l2"]
        )
    )
end
# Two grounds whose plan edges coincide exactly (x = 40 on both planes): one plan curve that
# is a metal edge of both planes (a "both" chain without a bump), and the empty half.
function coincident_edge_window()
    return read_polygon_set(
        Dict(
            "Version" => 1,
            "Name" => "coincident",
            "MatchingRadius" => 1.9,
            "Box" => Dict("X" => [0.0, 80.0], "Y" => [0.0, 30.0]),
            "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
            "Planes" => [
                Dict(
                    "Name" => "L1",
                    "SurfaceZ" => 0.0,
                    "Facing" => "up",
                    "SubstrateThickness" => 20.0,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0, 0], [40, 0], [40, 30], [0, 30]]
                        )
                    ]
                ),
                Dict(
                    "Name" => "L2",
                    "SurfaceZ" => 4.8,
                    "Facing" => "down",
                    "SubstrateThickness" => 20.0,
                    "Polygons" => [
                        Dict(
                            "Conductor" => "ground",
                            "Outer" => [[0, 0], [40, 0], [40, 30], [0, 30]]
                        )
                    ]
                )
            ],
            "Vacuum" => Dict("Below" => 0.0, "Above" => 0.0)
        )
    )
end

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

# The z ranges of the fabricated steps (trench bottom to metal top) of a set's planes.
function step_slabs(spec)
    return [
        minmax(
            plane.surface_z - plane.facing * spec.overetch,
            plane.surface_z + plane.facing * spec.metal_thickness
        ) for plane in spec.planes
    ]
end

# Boundary closure and tag consistency, every node used, the maximum edge-length ratio
# overall and off the metal / trench slabs (tetrahedra not inside a step's z range: the
# interfaces are kept in every stack, so a 30-um region triangle over the 0.05-um trench is
# an inherent 600 in both sweeps; everything else must stay within the designed t / r).
# `metal_substrate`: the attributes of the metal-substrate faces (from the manifest's
# attribute table: `*_substrate` surfaces).
function mesh_checks(path, slabs, metal_substrate)
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
        in_slab = any(
            minimum(z) >= slab[1] - 1.0e-9 && maximum(z) <= slab[2] + 1.0e-9 for
            slab in slabs
        )
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
            attribute == 3 ? length(o) == 1 :
            attribute in metal_substrate ? o == [1] : o == [2]
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
    mesh_checks(joinpath(directory, "graded.msh2"), step_slabs(tiny), Set([5]))
end

metal_substrate_attributes(manifest) =
    Set(
        entry["attribute"] for
        (name, entry) in manifest["attributes"] if entry["dimension"] == 2 &&
            endswith(name, "_substrate") &&
            name != "ground_substrate"
    ) ∪ Set([5])

@testset "graded stacks and row ranges (the Q2 probe's strip values)" begin
    spec = strip_window()
    stack = PWM.z_levels(spec, 10, 5)
    zs = stack.levels
    @test length(zs) == 40
    interfaces = [-0.05, 0.0, 0.1, stack.z_bottom, stack.z_top]
    heights = PWM.band_heights(0.01, 2.0, 7)
    step = (-0.05, 0.1)
    stacks = PWM.graded_stacks(zs, interfaces, heights, 1.0, 3.0, step, 30.0)
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
    # Row ranges from the fabricated step's faces (M1 review MAJOR-1): the trench foot at
    # -0.05 lies strictly inside row 1's range; with the metal faces alone (the probe's rule)
    # row 1 ended at -0.03 and the foot fell to row 2.
    ranges = [(zs[lo], zs[hi]) for (lo, hi) in stacks.ranges]
    step_ranges =
        [(-0.1, 0.15), (-0.2, 0.2), (-0.5, 0.5), (-1.0, 1.0), (-2.0, 2.0), (-5.0, 5.0)]
    @test all(all(isapprox.(ranges[k], step_ranges[k]; atol=1.0e-9)) for k = 1:6)
    @test ranges[7] == (-525.0, 525.0)
    @test ranges[1][1] < -0.05 - 3 * heights[1]
    metal_only = PWM.graded_stacks(zs, interfaces, heights, 1.0, 3.0, (0.0, 0.1), 30.0)
    @test zs[metal_only.ranges[1][1]] ≈ -0.03 atol = 1.0e-9
    @test zs[metal_only.ranges[2][1]] ≈ -0.1 atol = 1.0e-9
    @test metal_only.ranges[3:7] == stacks.ranges[3:7]
    # The top list of a column with n rows: Z_n within row n's range, the innermost active
    # row's stack beyond; a full column's is Z_K.
    tops = PWM.column_top_levels(stacks)
    @test tops[7] == stacks.levels[8]
    for n = 1:6
        @test issubset(stacks.levels[8], tops[n]) && issubset(tops[n], stacks.levels[n + 1])
        lo, hi = stacks.ranges[n]
        @test filter(i -> lo <= i <= hi, tops[n]) ==
              filter(i -> lo <= i <= hi, stacks.levels[n + 1])
        @test filter(i -> i > stacks.ranges[6][2], tops[n]) ==
              filter(i -> i > stacks.ranges[6][2], stacks.levels[8])
    end
    @test all(issubset(tops[n + 1], tops[n]) for n = 1:6)
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
    coarser = PWM.graded_stacks(zs, interfaces, heights, 2.0, 3.0, step, 30.0)
    @test all(length.(coarser.levels) .<= length.(stacks.levels))
    @test sum(length.(coarser.levels)) < sum(length.(stacks.levels))
    shorter = PWM.graded_stacks(zs, interfaces, heights, 1.0, 2.0, step, 30.0)
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
        step,
        30.0
    )
    # A box shorter than beta h_(K-1) beyond the step is refused (M1 review MINOR-6).
    short = filter(z -> -3.0 <= z <= 3.0, zs)
    message = try
        PWM.graded_stacks(short, [-0.05, 0.0, 0.1, -3.0, 3.0], heights, 1.0, 3.0, step, 30.0)
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("box face", message)
    section = PWM.cross_section(stacks, zs, heights, abs.(zs))
    @test length(section.nodes) == 132 && length(section.triangles) == 202
    @test sum(n for (_, n) in section.pair_triangles) == 202
    @test length(section.pair_triangles) == 1 + 2 * 6 + 6
    @test length(section.triangle_pair) == 202 &&
          count(pair -> pair[1] == 0, section.triangle_pair) ==
          sum(n for (pair, n) in section.pair_triangles if pair[1] == 0)
end

# The graded sweep against the tensor sweep of the same polygon set. The graded call gets
# `kwargs` on top of the plan-preserving options (region_grading / region_ring off), so the
# plan, the band and the z levels are those of the tensor mesh unless a test asks otherwise;
# with the plan unchanged the metal faces must be the same plan quads / triangles. The graded
# mesh must have fewer than `fewer_than` x the tensor's tetrahedra (0.5 on the one-plane windows
# with their 50-525-um boxes; the small two-plane test boxes save less).
function compare_sweeps(
    spec,
    radial_um,
    tangential_um,
    directory;
    fewer_than=0.5,
    kwargs...
)
    tensor = mesh_polygon_window(
        spec,
        radial_um,
        tangential_um,
        joinpath(directory, "tensor.msh2");
        verbose=false
    )
    options = merge((; region_grading=false, region_ring=false), (; kwargs...))
    graded = mesh_polygon_window(
        spec,
        radial_um,
        tangential_um,
        joinpath(directory, "graded.msh2");
        verbose=false,
        sweep=:graded,
        options...
    )
    @test tensor["sweep"] == "tensor" && graded["sweep"] == "graded"
    @test graded["nonpositive"] == 0 && tensor["nonpositive"] == 0
    @test graded["tetrahedra"] < fewer_than * tensor["tetrahedra"]
    @test keys(graded["surface_area_um2"]) == keys(tensor["surface_area_um2"])
    for (attribute, area) in tensor["surface_area_um2"]
        @test graded["surface_area_um2"][attribute] ≈ area rtol = 1.0e-9
    end
    for (attribute, volume) in tensor["volume_um3"]
        @test graded["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
    end
    @test (tensor["band"]["region_grading"] == false) &&
          (graded["band"]["region_grading"] == false) == !options.region_grading
    if !options.region_grading
        @test graded["plan_triangles"] == tensor["plan_triangles"]
        # One plane: the metal faces (the attributes from 4 on except 6 / 9) are the same plan
        # elements. Two planes: the metal faces of one plane under the OTHER plane's band are
        # the (0, K) cells' faces there (one quad per station: that plane's rows have ended at
        # the gap midpoint), not the tensor's row-split quads, so only the areas agree.
        if length(spec.planes) == 1
            for (attribute, count) in tensor["surface_attribute_counts"]
                parse(Int, attribute) >= 4 && attribute != "6" && attribute != "9" ||
                    continue
                @test graded["surface_attribute_counts"][attribute] == count
            end
        end
    end
    checks = mesh_checks(
        joinpath(directory, "graded.msh2"),
        step_slabs(spec),
        metal_substrate_attributes(graded)
    )
    @test checks.nodes == graded["nodes"] && checks.tetrahedra == graded["tetrahedra"]
    @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
    # The designed anisotropy t / r off the slabs: the longest band plan edge is a station t
    # plus the spread of the scaled corner columns at both ends (a column scaled by c puts its
    # row-K node up to c h_K along the band from its base: the top edge of a short segment
    # between two outward corners is t + 2 c h_K, the same plan edge as in the tensor sweep),
    # as a cell diagonal with h_K (the strip prisms split the band quads), plus 0.1 %; the slab
    # region prisms are bounded by the region size over the overetch (30 / 0.05 = 600). The
    # shortest edge off the slabs is a z spacing >= r (two planes: the other plane's 0.05-um
    # offsets hang on the metal edge line, >= r at the production r).
    band = graded["band"]
    thickness = maximum(band["heights_um"])
    scale =
        max(band["max_mitre_scale"], band["max_inward_scale"], band["max_wall_end_scale"])
    longest = hypot(tangential_um + 2.0 * scale * thickness, thickness)
    @test checks.max_edge_ratio_off_slab <= longest / radial_um * (1.0 + 1.0e-3)
    @test checks.max_edge_ratio <=
          max(longest / radial_um, PWM.REGION_MESH_SIZE_MAX_UM / spec.overetch) *
          (1.0 + 1.0e-3)
    record = graded["graded_sweep"]
    @test record["alpha"] == get(kwargs, :alpha, 1.0) &&
          record["beta"] == get(kwargs, :beta, 3.0)
    @test record["region_ring"] == options.region_ring
    @test record["band_tetrahedra"] + record["region_tetrahedra"] == graded["tetrahedra"]
    @test sum(c["chains"] for c in record["cross_sections"]) == record["chains"]
    @test sum(c["columns"] for c in record["cross_sections"]) == record["columns"]
    @test sum(c["capped_columns"] for c in record["cross_sections"]) ==
          record["capped_columns"]
    @test sum(c["fan_sectors"] for c in record["cross_sections"]) == record["fan_sectors"]
    length(spec.planes) == 1 && @test length(record["cross_sections"]) == 1
    @test record["band_swept_elements"] > 0 && record["band_strip_prisms"] > 0
    @test record["region_prisms_plain"] + record["region_prisms_hanging"] > 0
    # The swept volume closes on the plan's analytic volume (checked in process too).
    for (attribute, volume) in record["expected_volume_um3"]
        @test graded["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
    end
    return tensor, graded, checks
end

@testset "strip window: graded vs tensor at the probe's r10 / t5" begin
    directory = mktempdir()
    # Recorded: tensor 63,768 tets (the probe's baseline); the probe's graded alpha 1 / beta 3
    # had 32,460 with a structured interior finer than Gmsh's region (kept here).
    tensor, graded, checks = compare_sweeps(strip_window(), 0.01, 5.0, directory)
    # Gmsh's region triangulation (and so the counts) is pinned where the probe ran (macOS);
    # M1 had 23,216 graded tetrahedra with the metal-only range rule.
    Sys.isapple() && @test tensor["tetrahedra"] == 63768 && graded["tetrahedra"] == 24176
    record = graded["graded_sweep"]
    @test record["chains"] == 4 && record["columns"] == 44
    @test record["capped_columns"] == 0 && record["fan_sectors"] == 0
    @test record["band_collapse_prisms"] == 0 && record["band_strip_prisms_hanging"] > 0
    @test record["region_prisms_hanging"] > 0 &&
          record["region_hanging_node_incidences"] > 0
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

@testset "strip facing down (the mirrored construction, backside above)" begin
    tensor, down, _ = compare_sweeps(strip_down_window(), 0.02, 5.0, mktempdir())
    @test down["planes"][1]["facing"] == "down"
    @test haskey(down["surface_area_um2"], "9") && down["surface_area_um2"]["9"] ≈ 3000.0
    @test down["graded_sweep"]["capped_columns"] == 0
    # The row ranges are mirrored about the surface: the step is [-0.1, 0.05].
    ranges = only(down["graded_sweep"]["cross_sections"])["row_ranges_z_um"]
    @test ranges[1][2] > 0.05 + 3 * 0.02 && ranges[1][1] < -0.1 - 3 * 0.02
end

@testset "spike window (a fan, clearance caps, a 4x inward mitre, collapse, wall ends)" begin
    tensor, spike, _ = compare_sweeps(spike_window(), 0.02, 5.0, mktempdir())
    band = tensor["band"]
    @test band["fans"] == 1 && band["clearance_capped_columns"] > 0
    @test band["collapse_triangles"] > 0 && band["max_inward_scale"] > 4.0
    record = spike["graded_sweep"]
    @test record["fan_sectors"] == PWM.FAN_COLUMNS - 1
    @test record["capped_columns"] > 0 && record["band_collapse_prisms"] > 0
    @test sum(only(record["cross_sections"])["plan_nodes_by_capped_top_list"]) ==
          record["capped_columns"]
    @test record["chains"] == 2
end

@testset "channel window (facing caps, collapse both ways, capped tops in the region)" begin
    tensor, channel, _ = compare_sweeps(channel_window(), 0.02, 5.0, mktempdir())
    band = tensor["band"]
    @test band["collapse_triangles"] > 0 && band["fans"] == 0
    @test band["inward_corners"] == 8 && band["mitre_outward_corners"] == 8
    record = channel["graded_sweep"]
    @test record["capped_columns"] >= 8 && record["band_collapse_prisms"] > 0
    @test record["fan_sectors"] == 0
    # Region prisms with a capped top as a corner carry its hanging levels.
    @test record["region_prisms_hanging"] > 0
end

@testset "M3 defaults on the strip: region size grading and the region ring (graded only)" begin
    directory = mktempdir()
    tensor = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "tensor.msh2");
        verbose=false
    )
    plain = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "plain.msh2");
        verbose=false,
        sweep=:graded,
        region_grading=false,
        region_ring=false
    )
    ring = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "ring.msh2");
        verbose=false,
        sweep=:graded,
        region_grading=false
    )
    graded = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "graded.msh2");
        verbose=false,
        sweep=:graded
    )
    # The ring: the region nodes adjacent to the band move to Z_K (the same plan); the
    # grading: a finer region plan next to the band (more plan triangles and region nodes),
    # the band and the z levels unchanged; areas / volumes equal to the tensor's throughout.
    @test ring["plan_triangles"] == plain["plan_triangles"] == tensor["plan_triangles"]
    @test ring["graded_sweep"]["region_ring"] && !plain["graded_sweep"]["region_ring"]
    @test ring["graded_sweep"]["region_ring_plan_nodes"] > 0 &&
          plain["graded_sweep"]["region_ring_plan_nodes"] == 0
    @test ring["graded_sweep"]["plan_nodes_by_stack"][8] >
          plain["graded_sweep"]["plan_nodes_by_stack"][8]
    @test ring["tetrahedra"] > plain["tetrahedra"]
    @test graded["plan_triangles"] > tensor["plan_triangles"]
    @test graded["graded_sweep"]["region_plan_nodes"] >
          plain["graded_sweep"]["region_plan_nodes"]
    @test graded["band"]["region_grading"]["size_min_um"] == 5.0 &&
          graded["band"]["region_grading"]["size_max_um"] == PWM.REGION_MESH_SIZE_MAX_UM &&
          graded["band"]["region_grading"]["slope"] == PWM.REGION_GRADING_SLOPE
    @test graded["band"]["heights_um"] == tensor["band"]["heights_um"] &&
          graded["band"]["quads"] == tensor["band"]["quads"]
    @test graded["z_levels"] == tensor["z_levels"]
    for m in (ring, graded), (attribute, area) in tensor["surface_area_um2"]
        @test m["surface_area_um2"][attribute] ≈ area rtol = 1.0e-9
    end
    for m in (ring, graded), (attribute, volume) in tensor["volume_um3"]
        @test m["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
    end
    checks = mesh_checks(
        joinpath(directory, "graded.msh2"),
        step_slabs(strip_window()),
        metal_substrate_attributes(graded)
    )
    @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
    @test checks.max_edge_ratio <= 600.0 * (1.0 + 1.0e-3)
    @test graded["tetrahedra"] < 0.5 * tensor["tetrahedra"]
    # The tensor sweep's plan cannot be graded (the recorded family's plan is fixed).
    @test_throws ErrorException mesh_polygon_window(
        strip_window(),
        0.05,
        5.0,
        joinpath(directory, "bad.msh2");
        verbose=false,
        region_grading=true
    )
end

@testset "region z-grading near the metal planes (decision 276 option (ii), off by default)" begin
    directory = mktempdir()
    tensor = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "tensor.msh2");
        verbose=false
    )
    ring = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "ring.msh2");
        verbose=false,
        sweep=:graded,
        region_grading=false
    )
    graded = mesh_polygon_window(
        strip_window(),
        0.01,
        5.0,
        joinpath(directory, "zg.msh2");
        verbose=false,
        sweep=:graded,
        region_grading=false,
        region_z_grading=true
    )
    @test ring["graded_sweep"]["region_z_grading"] == false
    record = graded["graded_sweep"]["region_z_grading"]
    @test record["first_cell_widths"] == PWM.REGION_Z_FIRST_CELL_WIDTHS &&
          record["growth"] == PWM.REGION_Z_GROWTH
    levels = record["levels_z_um"]
    zk = graded["graded_sweep"]["stack_levels_z_um"][8]
    # Contains Z_K; the first region cell beyond the step is at most w_K = 0.64 um tall (the
    # farthest production level within it: 0.5 um above the metal top, 0.5 below the trench);
    # no level inside the step is added; geometric growth beyond (1, 2, 5, 10, ...).
    @test all(any(z -> abs(z - w) <= 1.0e-9, levels) for w in zk)
    @test any(z -> abs(z - 0.5) <= 1.0e-9, levels) &&
          any(z -> abs(z + 0.5) <= 1.0e-9, levels)
    @test !any(z -> -0.05 + 1.0e-9 < z < 0.1 - 1.0e-9 && abs(z) > 1.0e-9, levels)
    above = sort(filter(z -> z > 0.1 + 1.0e-9, levels))
    @test above[1:5] ≈ [0.5, 1.0, 2.0, 5.0, 10.0]
    # Every region node carries that stack; the plan, areas and volumes are the tensor's.
    @test graded["plan_triangles"] == tensor["plan_triangles"]
    @test graded["tetrahedra"] > ring["tetrahedra"]
    for (attribute, volume) in tensor["volume_um3"]
        @test graded["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
    end
    checks = mesh_checks(
        joinpath(directory, "zg.msh2"),
        step_slabs(strip_window()),
        metal_substrate_attributes(graded)
    )
    @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
end

@testset "two planes: stacks keep the gap midpoint, one plane's rows end there at the latest" begin
    spec = two_plane_window()
    stack = PWM.z_levels(spec, 2, 1)
    zs = stack.levels
    @test stack.gap_midpoint ≈ 2.4
    @test any(z -> abs(z - 2.4) <= 1.0e-9, zs)
    interfaces = [-0.05, 0.0, 0.1, 4.85, 4.8, 4.7, 2.4, stack.z_bottom, stack.z_top]
    heights = PWM.band_heights(0.05, 2.0, 5)
    # L1's band (facing up): capped at the midpoint from above; L2's from below; a band of both
    # planes' edges spans both steps.
    l1 = PWM.graded_stacks(
        zs,
        interfaces,
        heights,
        1.0,
        3.0,
        (-0.05, 0.1),
        30.0;
        range_cap=(-Inf, 2.4)
    )
    l2 = PWM.graded_stacks(
        zs,
        interfaces,
        heights,
        1.0,
        3.0,
        (4.7, 4.85),
        30.0;
        range_cap=(2.4, Inf)
    )
    both = PWM.graded_stacks(zs, interfaces, heights, 1.0, 3.0, (-0.05, 4.85), 30.0)
    @test l1.levels == l2.levels == both.levels
    for s = 2:length(l1.levels)
        @test all(any(i -> abs(zs[i] - w) <= 1.0e-9, l1.levels[s]) for w in interfaces)
    end
    @test all(zs[l1.ranges[k][2]] <= 2.4 + 1.0e-9 for k = 1:4)
    @test all(zs[l2.ranges[k][1]] >= 2.4 - 1.0e-9 for k = 1:4)
    @test zs[l1.ranges[4][2]] ≈ 2.4 && zs[l2.ranges[4][1]] ≈ 2.4
    @test all(zs[both.ranges[k][1]] < -0.05 && zs[both.ranges[k][2]] > 4.85 for k = 1:4)
    @test all(zs[l1.ranges[k][1]] == zs[both.ranges[k][1]] for k = 1:4)
    # beta 4: row 4's limit 3.1 um passes the midpoint and is capped there; beta 8 would need
    # rows 3 and 4 both on the midpoint: refused.
    capped = PWM.graded_stacks(
        zs,
        interfaces,
        heights,
        1.0,
        4.0,
        (-0.05, 0.1),
        30.0;
        range_cap=(-Inf, 2.4)
    )
    @test zs[capped.ranges[4][2]] ≈ 2.4 && zs[capped.ranges[3][2]] < 2.4
    message = try
        PWM.graded_stacks(
            zs,
            interfaces,
            heights,
            1.0,
            8.0,
            (-0.05, 0.1),
            30.0;
            range_cap=(-Inf, 2.4)
        )
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("gap midpoint", message)
    distance = [min(abs(z), abs(z - 4.8)) for z in zs]
    for stacks in (l1, l2, both)
        section = PWM.cross_section(stacks, zs, heights, distance)
        @test length(section.nodes) > 0 &&
              sum(n for (_, n) in section.pair_triangles) == length(section.triangles)
    end
end

@testset "two planes: the synthetic two-level window (crossing edges, a bump; every chain of both planes)" begin
    spec = read_polygon_set(synthetic_two_level_window())
    # Every chain keeps every row across the 4.8-um gap and the box is 20 um of substrate on
    # either side: the saving is in the far field only (1.6x here).
    tensor, graded, checks = compare_sweeps(spec, 0.05, 5.0, mktempdir(); fewer_than=0.7)
    record = graded["graded_sweep"]
    @test record["gap_midpoint_z_um"] ≈ 2.4
    @test length(record["cross_sections"]) == 1
    section = only(record["cross_sections"])
    @test section["planes"] == ["L1", "L2"]
    @test section["step_faces_z_um"] ≈ [-0.05, 4.85]
    @test all(r[1] < -0.05 && r[2] > 4.85 for r in section["row_ranges_z_um"][1:(end - 1)])
    @test record["capped_columns"] > 0 && record["band_collapse_prisms"] > 0
    @test haskey(graded["surface_area_um2"], "10") &&
          haskey(graded["surface_area_um2"], "11")
    @test graded["bumps"] == 1
end

@testset "two planes: pure L1 / pure L2 / both chains, facing tops of both planes, a bump" begin
    spec = two_plane_window()
    tensor, graded, checks = compare_sweeps(spec, 0.05, 5.0, mktempdir(); fewer_than=0.7)
    record = graded["graded_sweep"]
    sections = Dict(join(c["planes"], "+") => c for c in record["cross_sections"])
    @test Set(keys(sections)) == Set(["L1", "L2", "L1+L2"])
    l1, l2, both = sections["L1"], sections["L2"], sections["L1+L2"]
    @test l1["chains"] >= 1 && l2["chains"] >= 1 && both["chains"] >= 1
    @test all(r[2] <= 2.4 + 1.0e-9 for r in l1["row_ranges_z_um"][1:(end - 1)])
    @test all(r[1] >= 2.4 - 1.0e-9 for r in l2["row_ranges_z_um"][1:(end - 1)])
    @test all(r[1] < -0.05 && r[2] > 4.85 for r in both["row_ranges_z_um"][1:(end - 1)])
    @test l1["step_faces_z_um"] ≈ [-0.05, 0.1] && l2["step_faces_z_um"] ≈ [4.7, 4.85]
    # The 2-um-apart hole edges cap both planes' columns (0.4 x 2 um = 0.8 um -> 4 of 5 rows)
    # and the region strip between their tops has triangles whose three corners are tops of
    # different variants: nested interval by interval only (the intersection rule).
    @test l1["capped_columns"] > 0 && l2["capped_columns"] > 0
    @test record["region_prisms_hanging"] > 0
    @test graded["bumps"] == 1 && graded["surface_area_um2"]["4"] > 0
end

@testset "two planes: a coincident cross-plane run (one plan curve, a metal edge of both planes)" begin
    tensor, graded, checks =
        compare_sweeps(coincident_edge_window(), 0.05, 5.0, mktempdir(); fewer_than=0.7)
    record = graded["graded_sweep"]
    @test length(record["cross_sections"]) == 1
    @test only(record["cross_sections"])["planes"] == ["L1", "L2"]
    @test tensor["cross_plane_reconciliation"]["coincident_segments"] > 0
    @test graded["surface_attribute_counts"]["4"] == tensor["surface_attribute_counts"]["4"]
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
    # A flip-chip gap too small for beta h_(K-2): two rows would end on the midpoint.
    message = try
        mesh_polygon_window(
            two_plane_window(; gap_z=1.4),
            0.05,
            5.0,
            joinpath(directory, "narrow.msh2");
            verbose=false,
            sweep=:graded
        )
        ""
    catch err
        sprint(showerror, err)
    end
    @test occursin("gap midpoint", message)
    # The default sweep reproduces the recorded tensor mesh of the probe window bit for bit
    # (reference-quality-20261002/mesher-scoping/probe/strip_tensor_r10_t5.json, generated on
    # macOS; Gmsh's region triangulation may differ on another platform, so the sha is pinned
    # there only and the counts everywhere). The M3 options (region grading / ring) do not
    # touch the tensor sweep: the same sha with them named explicitly off.
    for (name, kwargs) in
        (("strip_tensor", (;)), ("strip_tensor_off", (; region_ring=false)))
        tensor = mesh_polygon_window(
            strip_window(),
            0.01,
            5.0,
            joinpath(directory, "$name.msh2");
            verbose=false,
            kwargs...
        )
        @test tensor["sweep"] == "tensor" && !haskey(tensor, "graded_sweep")
        @test tensor["band"]["region_grading"] == false
        @test tensor["plan_perimeter_edges"] == 20 && tensor["radial_layers"] == 7
        if Sys.isapple()
            @test tensor["nodes"] == 12418 && tensor["tetrahedra"] == 63768
            @test tensor["bytes"] == 2635089
            @test tensor["sha256"] ==
                  "8f234e5b0933f90aaef640399f991e73f950b5fb102766f840565d9b934980f6"
        end
    end
end

# A square hole in a ground: four plan vertices, each a 90-degree outward corner of the gap
# partition and an inward corner of the ground partition (every incident curve graded from
# both sides); the box corners lie on the window wall (no vertices).
hole_window() = one_plane_window(
    "hole",
    [0.0, 30.0],
    [0.0, 30.0],
    [
        Dict(
            "Conductor" => "ground",
            "Outer" => [[0, 0], [30, 0], [30, 30], [0, 30]],
            "Holes" => [[[10, 10], [20, 10], [20, 20], [10, 20]]]
        )
    ];
    substrate=50.0,
    above=50.0
)

@testset "D1: step_face_z_levels and the mirrored z stack (decision 406; legacy default)" begin
    spec = strip_window()
    # r10: metal top 0.1 + 0.01 / 0.03 (0.07 reaches the fixed 0.15 level), trench bottom
    # -0.05 - 0.01 / 0.03; r5: the probe variant's levels 0.105 / 0.115 / 0.135 and -0.055 /
    # -0.065 / -0.085 (reference-quality REPORT section 4); r50: the first cell is already r.
    @test PWM.step_face_z_levels(spec, 0.01) ≈ [-0.08, -0.06, 0.11, 0.13]
    @test PWM.step_face_z_levels(spec, 0.005) ≈
          [-0.085, -0.065, -0.055, 0.105, 0.115, 0.135]
    @test isempty(PWM.step_face_z_levels(spec, 0.05))
    # Facing down: mirrored about the surface (metal top -0.1, trench bottom +0.05).
    @test PWM.step_face_z_levels(strip_down_window(), 0.01) ≈ [-0.13, -0.11, 0.06, 0.08]
    # Two planes: both steps, the L2 (z 4.8, facing down) levels mirrored; none beyond the
    # gap midpoint (2.4) or the fixed offsets.
    two = PWM.step_face_z_levels(two_plane_window(), 0.01)
    @test two ≈ [-0.08, -0.06, 0.11, 0.13, 4.67, 4.69, 4.86, 4.88]
    # The legacy stack is untouched by the keyword's default; the mirrored stack merges the
    # levels (sorted, unique) and keeps every legacy level.
    legacy = PWM.z_levels(spec, 10, 5)
    mirrored =
        PWM.z_levels(spec, 10, 5; step_face_levels=PWM.step_face_z_levels(spec, 0.01))
    @test legacy.levels == PWM.z_levels(spec, 10, 5; step_face_levels=Float64[]).levels
    @test length(mirrored.levels) == length(legacy.levels) + 4
    @test all(any(z -> abs(z - w) <= 1.0e-9, mirrored.levels) for w in legacy.levels)
    @test issorted(mirrored.levels) && all(diff(mirrored.levels) .> 1.0e-9)
    @test (mirrored.z_bottom, mirrored.z_top, mirrored.backsides) ==
          (legacy.z_bottom, legacy.z_top, legacy.backsides)
    @test isequal(mirrored.gap_midpoint, legacy.gap_midpoint)
    @test_throws ErrorException mesh_polygon_window(
        spec,
        0.01,
        5.0,
        joinpath(mktempdir(), "bad.msh2");
        verbose=false,
        step_face_z_grading=:graded
    )
end

@testset "D1: the mirrored stack in both sweeps; legacy named explicitly = the default bytes" begin
    directory = mktempdir()
    for (sweep, kwargs) in ((:tensor, (;)), (:graded, (;)))
        default = mesh_polygon_window(
            strip_window(),
            0.01,
            5.0,
            joinpath(directory, "$(sweep)_default.msh2");
            verbose=false,
            sweep=sweep
        )
        legacy = mesh_polygon_window(
            strip_window(),
            0.01,
            5.0,
            joinpath(directory, "$(sweep)_legacy.msh2");
            verbose=false,
            sweep=sweep,
            step_face_z_grading=:legacy,
            vertex_column_grading=false
        )
        mirrored = mesh_polygon_window(
            strip_window(),
            0.01,
            5.0,
            joinpath(directory, "$(sweep)_mirrored.msh2");
            verbose=false,
            sweep=sweep,
            step_face_z_grading=:mirrored
        )
        # The options named at their defaults reproduce the default bytes (the S2.1 family).
        @test legacy["sha256"] == default["sha256"] && legacy["bytes"] == default["bytes"]
        @test default["step_face_z_grading"] == "legacy" &&
              isempty(default["step_face_z_levels_um"])
        @test default["band"]["vertex_column_grading"] == false
        @test mirrored["step_face_z_grading"] == "mirrored"
        @test mirrored["step_face_z_levels_um"] ≈ [-0.08, -0.06, 0.11, 0.13]
        @test mirrored["sha256"] != default["sha256"]
        # The z stack gains exactly the mirrored levels; the plan is unchanged.
        @test length(mirrored["z_levels"]) == length(default["z_levels"]) + 4
        @test all(
            any(z -> abs(z - w) <= 1.0e-9, mirrored["z_levels"]) for w in [0.11, 0.13]
        )
        @test mirrored["plan_triangles"] == default["plan_triangles"]
        # The first cell beyond the metal top is now r tall (it was 50 nm).
        above = sort(filter(z -> z > 0.1 + 1.0e-9, mirrored["z_levels"]))
        @test above[1] ≈ 0.11
        @test sort(filter(z -> z > 0.1 + 1.0e-9, default["z_levels"]))[1] ≈ 0.15
        # Areas and volumes are the default's; the mesh is closed, tagged and non-degenerate.
        for (attribute, area) in default["surface_area_um2"]
            @test mirrored["surface_area_um2"][attribute] ≈ area rtol = 1.0e-9
        end
        for (attribute, volume) in default["volume_um3"]
            @test mirrored["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
        end
        @test mirrored["nonpositive"] == 0
        checks = mesh_checks(
            joinpath(directory, "$(sweep)_mirrored.msh2"),
            step_slabs(strip_window()),
            metal_substrate_attributes(mirrored)
        )
        @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
        # Cost: the tensor sweep pays four full plan layers; the graded sweep thins the levels
        # away from the edge (the probe variant: +3.7 % at r20, +1.5 % at r5, -3.7 % at r1.25).
        if sweep == :tensor
            # Four more intervals of 3 tets per plan triangle (less the excluded metal).
            @test default["tetrahedra"] <
                  mirrored["tetrahedra"] <=
                  default["tetrahedra"] + 4 * 3 * default["plan_triangles"]
        else
            @test 0.9 * default["tetrahedra"] <
                  mirrored["tetrahedra"] <
                  1.15 * default["tetrahedra"]
            @test mirrored["graded_sweep"]["stack_sizes"][1] ==
                  default["graded_sweep"]["stack_sizes"][1] + 4
        end
    end
end

@testset "D3: vertex ladder stations and the graded curve nodes" begin
    r, t = 0.02, 5.0
    # A right angle: the first station is 3r (2 r / tan(45 deg) = 2 r, the row-1 node of the
    # inward mitre lies r along the arm), then 7r, 15r, ... while the spacing stays within t
    # and half of it fits before the curve's middle (L = 10: 2.54 + 1.28 <= 5, 5.1 + 2.56 > 5).
    right = PWM.vertex_ladder_stations(10.0, t, r, pi / 2)
    @test right ≈ r .* [3, 7, 15, 31, 63, 127]
    @test PWM.vertex_ladder_stations(100.0, t, r, pi / 2) ≈
          r .* [3, 7, 15, 31, 63, 127, 255]
    # An acute wedge (52.8 deg, the S5 tip): 2 / tan(26.4) = 4.03 -> the first station 7r; an
    # obtuse one (135 deg): 2 / tan(67.5) = 0.83 -> r.
    @test PWM.vertex_ladder_stations(100.0, t, r, deg2rad(52.8))[1] ≈ 7r
    @test PWM.vertex_ladder_stations(100.0, t, r, deg2rad(135.0))[1] ≈ r
    # A curve too short for any station keeps none.
    @test isempty(PWM.vertex_ladder_stations(0.1, t, r, pi / 2))
    # Every graded column at a right angle carries j - 1 rows under the clearance rule.
    heights = PWM.band_heights(r, 2.0, 7)
    for (j, s) in enumerate(right)
        rows, clamped = PWM.rows_within(heights, PWM.CORNER_CLEARANCE_FRACTION * s)
        @test rows == j && !clamped
    end
    # The curve nodes: both ends graded, the middle uniform within t, the legacy nodes when
    # neither end is a vertex (bit for bit).
    curve = PWM.PlanCurve(Int32(1), (0.0, 0.0), (12.0, 0.0), (Int32(1), Int32(2)), true)
    legacy = PWM.curve_transfinite_nodes(curve, t)
    nodes, added = PWM.vertex_graded_curve_nodes(curve, t, r, NaN, NaN)
    @test nodes == legacy && added == 0
    nodes, added = PWM.vertex_graded_curve_nodes(curve, t, r, pi / 2, pi / 2)
    xs = [n[1] for n in nodes]
    @test added == 12 && length(nodes) == 2 + 12 + 1
    @test xs[1] == 0.0 && xs[end] == 12.0 && issorted(xs) && all(diff(xs) .> 0.0)
    @test xs[2:7] ≈ r .* [3, 7, 15, 31, 63, 127]
    @test xs[(end - 6):(end - 1)] ≈ 12.0 .- r .* [127, 63, 31, 15, 7, 3]
    @test maximum(diff(xs)) <= t * (1.0 + 1.0e-6)
    # Only one vertex end: the other end keeps the regular spacing.
    nodes, added = PWM.vertex_graded_curve_nodes(curve, t, r, NaN, pi / 2)
    @test added == 6 && nodes[1] == (0.0, 0.0) && length(nodes) == 2 + 1 + 6
    # The middle (0 .. 12 - 2.54 = 9.46) in two segments of 4.73 <= t.
    @test nodes[2][1] ≈ 4.73 && nodes[3][1] ≈ 9.46
end

@testset "D3: vertex column grading on the hole, the spike and the two-plane window" begin
    directory = mktempdir()
    # The hole: 4 vertices, 8 graded curve ends, 6 stations each (r 0.02, t 5, L 10).
    legacy = mesh_polygon_window(
        hole_window(),
        0.02,
        5.0,
        joinpath(directory, "hole_legacy.msh2");
        verbose=false,
        sweep=:graded
    )
    graded = mesh_polygon_window(
        hole_window(),
        0.02,
        5.0,
        joinpath(directory, "hole_d3.msh2");
        verbose=false,
        sweep=:graded,
        vertex_column_grading=true
    )
    @test legacy["band"]["vertex_column_grading"] == false
    record = graded["band"]["vertex_column_grading"]
    @test record["vertices"] == 4 && record["graded_curve_ends"] == 8
    @test record["ladder_stations"] == 48 && record["min_wedge_angle_deg"] ≈ 90.0
    @test record["min_turn_deg"] == PWM.DEFAULT_VERTEX_GRADING_MIN_TURN_DEG
    # Columns on both sides of every curve: 2 + 12 ladder nodes and no middle node (the
    # remaining 4.92 um is one segment), 13 segments per curve instead of 2.
    @test graded["graded_sweep"]["columns"] ==
          2 * 4 * 13 ==
          legacy["graded_sweep"]["columns"] + 2 * 48 - 2 * 4
    # The inward side's clearance rule caps the graded columns (j - 1 rows), the outward
    # side's mitres are unchanged; no clamped column, no fan.
    @test graded["band"]["clearance_capped_columns"] >
          legacy["band"]["clearance_capped_columns"]
    @test graded["band"]["clearance_clamped_columns"] == 0 && graded["band"]["fans"] == 0
    @test graded["band"]["inward_corners"] == 4 &&
          graded["band"]["mitre_outward_corners"] == 4
    @test graded["graded_sweep"]["capped_columns"] >
          legacy["graded_sweep"]["capped_columns"]
    @test graded["band"]["quad_min_abs_sin"] > 0.0 && graded["nonpositive"] == 0
    # Same geometry: areas and volumes equal the legacy mesh's; the z stack is untouched.
    for (attribute, area) in legacy["surface_area_um2"]
        @test graded["surface_area_um2"][attribute] ≈ area rtol = 1.0e-9
    end
    for (attribute, volume) in legacy["volume_um3"]
        @test graded["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
    end
    @test graded["z_levels"] == legacy["z_levels"]
    checks = mesh_checks(
        joinpath(directory, "hole_d3.msh2"),
        step_slabs(hole_window()),
        metal_substrate_attributes(graded)
    )
    @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
    # Below the turn threshold nothing is graded (the hole's corners turn by 90 deg).
    coarse = mesh_polygon_window(
        hole_window(),
        0.02,
        5.0,
        joinpath(directory, "hole_t95.msh2");
        verbose=false,
        sweep=:graded,
        vertex_column_grading=true,
        vertex_grading_min_turn_deg=95.0
    )
    @test coarse["band"]["vertex_column_grading"]["vertices"] == 0
    @test coarse["sha256"] == legacy["sha256"]
    # The spike (an acute 28-degree tip: 15r first station; a fan on the vacuum side) and the
    # two-plane window (16 vertices of 90 deg; the bump's 16-gon joints turn by 22.5 deg and are
    # not vertices) with D1 and D3 together: built, closed, volumes equal to the tensor's.
    for (name, spec, r) in
        (("spike", spike_window(), 0.02), ("two_plane", two_plane_window(), 0.02))
        tensor = mesh_polygon_window(
            spec,
            r,
            5.0,
            joinpath(directory, "$(name)_tensor.msh2");
            verbose=false
        )
        both = mesh_polygon_window(
            spec,
            r,
            5.0,
            joinpath(directory, "$(name)_d1d3.msh2");
            verbose=false,
            sweep=:graded,
            step_face_z_grading=:mirrored,
            vertex_column_grading=true
        )
        record = both["band"]["vertex_column_grading"]
        if name == "spike"
            @test record["vertices"] == 3
            @test record["min_wedge_angle_deg"] ≈ 28.07 atol = 0.01
            @test both["band"]["fans"] == 1
        else
            @test record["vertices"] == 16 && record["min_wedge_angle_deg"] ≈ 90.0
            @test length(both["step_face_z_levels_um"]) == 4 # r20: 0.02 / 0.06 beyond each face
        end
        @test both["nonpositive"] == 0
        for (attribute, volume) in tensor["volume_um3"]
            @test both["volume_um3"][attribute] ≈ volume rtol = 1.0e-9
        end
        for (attribute, area) in tensor["surface_area_um2"]
            @test both["surface_area_um2"][attribute] ≈ area rtol = 1.0e-9
        end
        checks = mesh_checks(
            joinpath(directory, "$(name)_d1d3.msh2"),
            step_slabs(spec),
            metal_substrate_attributes(both)
        )
        @test checks.untagged == 0 && checks.inconsistent == 0 && checks.unused_nodes == 0
    end
    # Two planes whose edges CROSS in plan (the synthetic two-level window): a crossing is four
    # incident curves, two collinear per plane — no plane has a corner there, so it is not a
    # vertex; the vertices are exactly the off-wall polygon corners of both planes (turn >= 30).
    crossing = read_polygon_set(synthetic_two_level_window())
    corners = 0
    for plane in crossing.planes,
        polygon in plane.polygons,
        ring in (polygon.outer, polygon.holes...)

        n = length(ring)
        for i = 1:n
            a, b, c = ring[mod1(i - 1, n)], ring[i], ring[mod1(i + 1, n)]
            any(
                abs(b[j] - crossing.box[k]) <= 1.0e-5 for
                (j, k) in ((1, 1), (1, 2), (2, 3), (2, 4))
            ) && continue
            u = PWM.unit((b[1] - a[1], b[2] - a[2]))
            v = PWM.unit((c[1] - b[1], c[2] - b[2]))
            acos(clamp(u[1] * v[1] + u[2] * v[2], -1.0, 1.0)) >= deg2rad(30.0) - 1.0e-9 &&
                (corners += 1)
        end
    end
    plan = mesh_polygon_window(
        crossing,
        0.05,
        5.0,
        joinpath(directory, "crossing.msh2");
        verbose=false,
        sweep=:graded,
        plan_only=true,
        vertex_column_grading=true
    )
    # 16 polygon corners; the L1 and L2 traces share the corners (40, 28) and (40, 32) -> one
    # plan point each: 14 vertices.
    @test corners == 16 && plan["band"]["vertex_column_grading"]["vertices"] == 14
    @test plan["band"]["inward_corners"] > 2 * corners # the crossings are band corners, not vertices
    # The option needs the own band.
    @test_throws ErrorException mesh_polygon_window(
        hole_window(),
        0.05,
        5.0,
        joinpath(directory, "gmsh.msh2");
        verbose=false,
        band_mode=:gmsh,
        vertex_column_grading=true
    )
end

# The production entry point (the stage-1 lanes and the PBS scripts call it, no suite loads
# it otherwise): the CLI must load and build the plan in both sweeps, and its printed manifest
# must match the in-process one. --plan-only keeps each subprocess to the Julia start-up plus
# the plan (no volume mesh).
@testset "CLI smoke test: mesh_polygon_window.jl --plan-only in both sweeps" begin
    directory = mktempdir()
    window = joinpath(directory, "strip.json")
    open(io -> JSON.print(io, strip_window_data()), window, "w")
    julia = `$(Base.julia_cmd()) --project=$(@__DIR__)`
    script = joinpath(@__DIR__, "mesh_polygon_window.jl")
    for (sweep, options) in (
        ("tensor", String[]),
        ("graded", ["--sweep", "graded"]),
        (
            "graded_d1d3",
            [
                "--sweep",
                "graded",
                "--step-face-z-grading",
                "mirrored",
                "--vertex-column-grading",
                "on",
                "--vertex-grading-min-turn-deg",
                "45"
            ]
        )
    )
        d1d3 = sweep == "graded_d1d3"
        sweep = d1d3 ? "graded" : sweep
        output = joinpath(directory, "$sweep.msh2")
        log = joinpath(directory, "$sweep.log")
        command = `$julia $script $window 0.05 5.0 $output --plan-only $options`
        process = run(pipeline(command; stdout=log, stderr=log); wait=false)
        wait(process)
        @test success(process)
        # The manifest is the last JSON object printed (the mesher's progress lines precede it).
        lines = readlines(log)
        start = findlast(==("{"), lines)
        @test start !== nothing
        success(process) && start !== nothing || continue
        printed = JSON.parse(join(lines[start:end], '\n'))
        expected = mesh_polygon_window(
            strip_window(),
            0.05,
            5.0,
            output;
            verbose=false,
            plan_only=true,
            sweep=Symbol(sweep),
            step_face_z_grading=d1d3 ? :mirrored : :legacy,
            vertex_column_grading=d1d3,
            vertex_grading_min_turn_deg=d1d3 ? 45.0 : 30.0
        )
        @test printed["sweep"] == sweep == expected["sweep"]
        @test printed["step_face_z_grading"] ==
              expected["step_face_z_grading"] ==
              (d1d3 ? "mirrored" : "legacy")
        @test printed["step_face_z_levels_um"] == expected["step_face_z_levels_um"]
        @test printed["band"]["vertex_column_grading"] ==
              expected["band"]["vertex_column_grading"]
        d1d3 && @test printed["band"]["vertex_column_grading"]["min_turn_deg"] == 45.0
        # The region size grading is on for the graded sweep only (a record of its parameters)
        # and off (false) for the tensor sweep.
        if sweep == "graded"
            @test printed["band"]["region_grading"]["size_min_um"] == 5.0
        else
            @test printed["band"]["region_grading"] == false
        end
        @test printed["band"]["region_grading"] == expected["band"]["region_grading"]
        for key in ("plan_nodes", "plan_triangles", "plan_perimeter_edges", "radial_layers")
            @test printed[key] == expected[key]
        end
        @test printed["z_levels"] == expected["z_levels"]
        @test printed["plan_perimeter_edges"] == 20 && printed["radial_layers"] == 5
    end
end

if LONG_TESTS
    @testset "validator and census on the strip and the two-plane window (long)" begin
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
                sweep=:graded,
                region_grading=false,
                region_ring=false
            ),
            "defaults" => mesh_polygon_window(
                strip_window(),
                0.01,
                5.0,
                joinpath(directory, "defaults.msh2");
                verbose=false,
                sweep=:graded
            ),
            # r 0.02: the census takes the metal slabs from the sidewall faces of the
            # r-resolved levels (<= 0.02 um tall).
            "two_plane" => mesh_polygon_window(
                two_plane_window(),
                0.02,
                5.0,
                joinpath(directory, "two_plane.msh2");
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
        # The region grading shortens the slab region prisms' plan edges (30 -> ~t um).
        @test overall("defaults")["edge_ratio_max"] <=
              all_graded["edge_ratio_max"] * (1.0 + 1.0e-3)
        # Two planes: the census finds both metal slabs (the bump columns do not merge them).
        @test length(census["two_plane"]["metal_slabs_z"]) == 2
    end
end
