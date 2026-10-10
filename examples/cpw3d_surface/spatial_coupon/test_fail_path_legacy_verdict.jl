# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# F5T (supervisor decisions 674 (1) / (2) and 675; the diagnosis of record
# round4/thincornerF5/DIAGNOSIS.md 3 / 3.1): the descent record (R-1) and the fail-path
# legacy verdict (R-7) on the STORED star of the 39e4806d547d thin corner
# [-2.7004415, 15.3200762] (PBS 61338, descent call 11: 14 corner cells, 87 guard cells, 36
# referenced vertices, 12 movable). The fixture is the dump's record byte for byte; the
# replay set is the diagnosis replay's (its (A) reproduced the production descent
# bitwise). T-F5T-2 = the 3.1 O-E numbers to the digit.

const FIXTURE = joinpath(@__DIR__, "testdata", "fail-path-legacy",
                         "descent-call11-39e4806d547d-thin.json")
const D = parse_json(read(FIXTURE, String))
asint(x) = Int(round(x))

function replay_set()
    referenced = [asint(i) for i in D["ReferencedVertices"]]
    loc = Dict(g => n for (n, g) in enumerate(referenced))
    nv = length(referenced)
    start = hcat([Vector{Float64}(p) for p in D["ReferencedStart"]]...)
    original = hcat([Vector{Float64}(p) for p in D["ReferencedOriginal"]]...)
    endpts = hcat([Vector{Float64}(p) for p in D["ReferencedEnd"]]...)
    guards = D["GuardCells"]
    cells = [ntuple(i -> loc[asint(g["Vertices"][i])], 4) for g in guards]
    floors = [Float64(g["Floor"]) for g in guards]
    ceilings = [g["Ceiling"] === nothing ? Inf : Float64(g["Ceiling"]) for g in guards]
    condition_ceilings = [g["ConditionCeiling"] === nothing ? Inf : Float64(g["ConditionCeiling"])
                          for g in guards]
    slot = Dict(Tuple(sort!(collect(c))) => k for (k, c) in enumerate(cells))
    targets = [slot[Tuple(sort!([loc[asint(v)] for v in c]))] for c in D["TargetCells"]]
    incident = [Int[] for _ in 1:nv]
    for (k, c) in enumerate(cells), i in c
        push!(incident[i], k)
    end
    bases = Dict{Int, Matrix{Float64}}(i => zeros(3, 0) for i in 1:nv)
    bounds = zeros(nv); edge_floors = zeros(nv)
    for m in D["Movable"]
        i = loc[asint(m["Vertex"])]
        bases[i] = hcat([Vector{Float64}(c) for c in m["Basis"]]...)
        bounds[i] = Float64(m["Bound"]); edge_floors[i] = Float64(m["EdgeFloor"])
    end
    return (; start, original, endpts, cells, floors, ceilings, condition_ceilings, targets,
            incident, bases, bounds, edge_floors)
end

const S = replay_set()
kappa_v0_max(points) = maximum(tetrahedron_aspect([points[:, i] for i in S.cells[k]]) for k in S.targets)
kappa_reg_max(points) = maximum(tetrahedron_regular_condition([points[:, i] for i in S.cells[k]])
                                for k in S.targets)

@testset "the replay set is the record's" begin
    @test length(S.targets) == 14
    @test length(S.cells) == 87
    @test count(i -> size(S.bases[i], 2) > 0, axes(S.start, 2)) == 12
    @test kappa_v0_max(S.start) == Float64(D["StartAchieved"])
    @test kappa_v0_max(S.endpts) == Float64(D["EndAchieved"])
    @test Float64(D["EndAchieved"]) == 4.547610478844198
end

@testset "R-1: the descent record leaves the arithmetic bitwise and names the stop" begin
    plain = copy(S.start)
    achieved_plain, moves_plain = optimize_seed_cells!(plain, S.original, S.cells, S.incident,
                                                       S.bases, S.floors, S.bounds, S.targets,
                                                       :aspect, Float64(D["Goal"]);
                                                       ceilings=S.ceilings,
                                                       condition_ceilings=S.condition_ceilings,
                                                       edge_floors=S.edge_floors)
    recorded = copy(S.start)
    record = Dict{String, Any}()
    achieved, moves = optimize_seed_cells!(recorded, S.original, S.cells, S.incident, S.bases,
                                           S.floors, S.bounds, S.targets, :aspect,
                                           Float64(D["Goal"]); ceilings=S.ceilings,
                                           condition_ceilings=S.condition_ceilings,
                                           edge_floors=S.edge_floors, record=record)
    # (A) of the diagnosis: the production descent reproduces the build of record bitwise.
    @test achieved_plain == 4.547610478844198
    @test moves_plain == 259
    @test plain == S.endpts
    # The record changes nothing: the same end state, value and moves.
    @test achieved == achieved_plain
    @test moves == moves_plain
    @test recorded == plain
    # The stop reason and the best true value seen = the instrumented dump of record.
    @test record["StopReason"] == "no improving move"
    @test record["Sweeps"] == asint(D["Sweeps"]) == 12
    @test record["BestAchieved"] == Float64(D["BestAchieved"]) == 4.492409188971175
    @test record["BestSweep"] == asint(D["BestSweep"]) == 6
    @test record["BestMove"] == asint(D["BestMove"]) == 206
    @test record["StartAchieved"] == Float64(D["StartAchieved"])
    @test record["EndAchieved"] == achieved
    @test record["Moves"] == moves
    @test record["Objective"] == "aspect"
    @test record["Goal"] == 3.8
    @test record["SweepLimit"] == 60
    # The best state is above the gate: a best-state return would not pass 4.0 (669 (1) / 674).
    @test record["BestAchieved"] > 4.0
    # A descent that reaches its goal stops on "goal" (the kappa_reg descent of T-F5T-2 below).
end

@testset "T-F5T-2: the fail-path legacy verdict on the stored star (3.1 O-E to the digit)" begin
    points = copy(S.endpts)
    measure = fail_path_legacy_verdict!(points, S.original, S.cells, S.incident, S.bases, S.floors,
                                        S.bounds, S.targets, 5.0; ceilings=S.ceilings,
                                        condition_ceilings=S.condition_ceilings,
                                        edge_floors=S.edge_floors)
    @test measure["Measure"] == "RegularCondition"
    @test measure["Rule"] == FAIL_PATH_LEGACY_VERDICT_RULE
    @test measure["Gate"] == 5.0
    @test measure["Target"] == INVARIANT_CORNER_TARGET == 3.8
    # The end state of record read by the invariant measure (DIAGNOSIS 1.2 (c)).
    @test measure["Before"] == 4.374913381572155
    @test measure["VertexFrameConditionBefore"] == 4.547610478844198
    # The kappa_reg descent from the end state (DIAGNOSIS 3.1 O-E): 3.775 in 48 moves.
    @test measure["After"] == 3.7752797439279706
    @test measure["Moves"] == 48
    @test measure["Passed"] === true
    @test measure["After"] <= INVARIANT_CORNER_TARGET
    @test measure["VertexFrameConditionAfter"] == 5.576930814895972
    @test measure["Descent"]["StopReason"] == "goal"
    @test measure["Descent"]["Objective"] == "regular_condition"
    @test measure["Descent"]["Moves"] == 48
    @test measure["Descent"]["EndAchieved"] == measure["After"]
    @test measure["Descent"]["BestAchieved"] <= measure["After"]
    # `points` holds the re-descended state: the readings are its readings.
    @test kappa_reg_max(points) == measure["After"]
    @test kappa_v0_max(points) == measure["VertexFrameConditionAfter"]
    # The guards held through the descent: every target cell above its floor, no vertex beyond
    # its bound, no edge of a moved vertex below its floor (the production acceptance rules).
    for k in S.targets
        @test tetrahedron_scaled_jacobian([points[:, i] for i in S.cells[k]]) >= S.floors[k] - 1.0e-9
    end
    for i in axes(points, 2)
        @test norm(points[:, i] .- S.original[:, i]) <= S.bounds[i] * (1.0 + 1.0e-12) ||
              size(S.bases[i], 2) == 0
    end
    # The anti-correlation of the two measures at this corner (DIAGNOSIS 3.1): the
    # stored-order frame rises above the legacy gate while kappa_reg reaches the goal.
    @test measure["VertexFrameConditionAfter"] > 4.0
    # The verdict is the gate's: below the achieved kappa_reg it fails, with the same descent.
    tight = copy(S.endpts)
    failed = fail_path_legacy_verdict!(tight, S.original, S.cells, S.incident, S.bases, S.floors,
                                       S.bounds, S.targets, 3.5; ceilings=S.ceilings,
                                       condition_ceilings=S.condition_ceilings,
                                       edge_floors=S.edge_floors)
    @test failed["After"] == measure["After"]
    @test failed["Passed"] === false
    @test tight == points
    # No gate (a build without --corner-shape-gate) is refused: today's fail-closed path.
    @test_throws ErrorException fail_path_legacy_verdict!(copy(S.endpts), S.original, S.cells,
                                                           S.incident, S.bases, S.floors, S.bounds,
                                                           S.targets, 0.0; ceilings=S.ceilings,
                                                           condition_ceilings=S.condition_ceilings,
                                                           edge_floors=S.edge_floors)
end

@testset "the rule text is a named constant (decision 670 discipline)" begin
    @test occursin("FAIL PATH only", FAIL_PATH_LEGACY_VERDICT_RULE)
    @test occursin("CornerShapeGate", FAIL_PATH_LEGACY_VERDICT_RULE)
    @test occursin("INVARIANT_CORNER_TARGET", FAIL_PATH_LEGACY_VERDICT_RULE)
    @test !occursin("fail path", CORNER_RECONNECTION_RULE)
    @test !occursin("fail path", INVARIANT_CORNER_RULE)
end
