# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

# WedgeSplit (WEDGE_SPLIT_RULE; supervisor decisions 649 (2) / 651 S-1 / 669 (4) / 670; the design
# of record round4/wedgesplit/DESIGN.md, pre-registered test T-W-1) on the two stored star dumps
# of the design lane (testdata/wedge-split, byte-identical copies of the diagnosis / design
# records): the 7fd4 fab corner 8 (the REAL failing cell 4276470, F5-B's only move refused by
# the improvement guard) and the synthetic 151.19-degree theta-0 fab with F5-B off (the class
# with a legal edge removal), plus two no-trigger controls (the passed star of the synthetic
# draw; 7fd4's corner 1). A star carries the post-descent / post-F5-B positions, the original
# Gmsh positions, the cells with their volume entities, the surface triangles / lines and the
# corner's parameters: the production per-cell arrays of optimize_required_region! at the
# reconnection rounds' exit are rebuilt from it exactly as the design tool did (the vertices
# with incomplete incidence and the CAD points frozen), and wedge_split_corner! runs on them.

const STAR_DIR = joinpath(@__DIR__, "testdata", "wedge-split")

xyz(points, cell) = [points[:, i] for i in cell]
kappa(points, cell) = tetrahedron_regular_condition(xyz(points, cell))

# The state of optimize_required_region! at the exit of the F5-B rounds, from a star dump. The
# point matrix already carries the spare column `s_index` (as the caller grows it before the
# pass); `corner_incident` reads that very matrix.
function star_state(name)
    star = parse_json(read(joinpath(STAR_DIR, name), String))
    parameters = star["Parameters"]
    n = length(star["Points"])
    points = hcat(reduce(hcat, [Vector{Float64}(p) for p in star["Points"]]), zeros(3))
    original = hcat(reduce(hcat, [Vector{Float64}(p) for p in star["OriginalPoints"]]), zeros(3))
    tetrahedra = [NTuple{4, Int}(Tuple(Int.(round.(c)))) for c in star["Cells"]]
    entity_of = Int.(round.(star["CellEntities"]))
    triangles = [NTuple{3, Int}(Tuple(Int.(round.(t)))) for t in star["Triangles"]]
    lines = [NTuple{2, Int}(Tuple(Int.(round.(l)))) for l in star["Lines"]]
    fixed = Vector{Bool}(star["Fixed"])
    core = Vector{Bool}(star["Core"])
    center = Vector{Float64}(star["Corner"])
    sides = (walls=[Vector{Float64}(w) for w in star["Sides"]["Walls"]],
             metal=Vector{Float64}(star["Sides"]["Metal"]))
    grading = CornerGrading(Float64(parameters["CornerGrading"]["CornerSize"]),
                            Float64(parameters["CornerGrading"]["Ratio"]),
                            Float64(parameters["CornerGrading"]["LcFine"]),
                            Float64(parameters["CornerGrading"]["Radius"]))
    gate = Float64(parameters["CornerShapeGate"])
    goal = Float64(parameters["InvariantCornerTarget"])
    scaled_target = 2.0 * Float64(parameters["MinimumScaledJacobian"])
    condition_target = 0.95 * Float64(parameters["MaximumJacobianCondition"])
    collapse_threshold = Float64(parameters["CollapseSizeFraction"]) * grading.corner_size
    tolerance = Float64(parameters["Tolerance"])
    displacement_ratio = Float64(parameters["DisplacementRatio"])
    incident = [Int[] for _ in 1:(n + 1)]
    for (k, cell) in enumerate(tetrahedra), i in cell
        push!(incident[i], k)
    end
    floors = [min(tetrahedron_scaled_jacobian(xyz(original, cell)), scaled_target) for cell in tetrahedra]
    ceilings = fill(Inf, length(tetrahedra))
    condition_ceilings = [max(tetrahedron_aspect(xyz(original, cell)), condition_target) for cell in tetrahedra]
    # The design tool's movable set: the vertices of the core cells (complete incidence in the
    # star) that are not CAD points; every other star vertex is frozen.
    complete = falses(n)
    for (k, cell) in enumerate(tetrahedra)
        core[k] && foreach(i -> complete[i] = true, cell)
    end
    bases = surface_movement_bases(points, triangles, lines)
    for i in 1:n
        (complete[i] && !fixed[i]) || (bases[i] = zeros(3, 0))
    end
    corner_distance = [norm(original[:, i] .- center) for i in 1:n]
    bound_of(point) = (distance = norm(point .- center);
                       displacement_ratio * min(grading.lc_fine,
                                                distance <= grading.radius ? corner_ball_size(grading, distance) :
                                                                             grading.lc_fine))
    bounds = vcat([bound_of(original[:, i]) for i in 1:n], 0.0)
    edge_floors = vcat([corner_distance[i] <= grading.radius ? collapse_threshold : 0.0 for i in 1:n], 0.0)
    surface_faces = Set(Tuple(sort!(collect(t))) for t in triangles)
    masks = sidewall_vertex_masks(points, triangles, center, sides, grading.radius)
    dead = Set{Int}()
    corner_incident(c) = [k for (k, cell) in enumerate(tetrahedra)
                          if !(k in dead) && any(norm(points[:, i] .- c) <= tolerance for i in cell)]
    mesh = (tetrahedra=tetrahedra, incident=incident, floors=floors, ceilings=ceilings,
            condition_ceilings=condition_ceilings, origin=collect(1:length(tetrahedra)),
            entity_of=entity_of, reconnection_replaced=Int[],
            reconnection_added=Dict{Int, Tuple{Int, NTuple{4, Int}}}(), dead=dead)
    targets = corner_incident(center)
    bridging = bridging_sliver_cells(points, tetrahedra, targets, masks, center, sides)
    trigger = [k for k in bridging if kappa(points, tetrahedra[k]) > gate]
    return (star=star, points=points, original=original, s_index=n + 1, mesh=mesh, bases=bases, bounds=bounds,
            edge_floors=edge_floors, bound_of=bound_of, surface_faces=surface_faces, triangles=triangles,
            center=center, sides=sides, masks=masks, trigger=trigger, targets=targets,
            corner_incident=corner_incident, gate=gate, goal=goal, scaled_target=scaled_target,
            condition_target=condition_target, grading=grading, collapse_threshold=collapse_threshold,
            tolerance=tolerance)
end

function run_wedge_split!(state; remapped=Set{Int}())
    return wedge_split_corner!(state.points, state.original, state.s_index, state.mesh, state.bases,
                               state.bounds, state.edge_floors, state.bound_of, remapped, state.surface_faces,
                               state.center, state.sides, state.masks, state.trigger, state.corner_incident,
                               state.gate, state.goal, state.scaled_target, state.condition_target,
                               state.grading.corner_size, state.collapse_threshold, state.tolerance)
end

# The number of cells owning each surface triangle (a conformity witness: the insertion may
# change no surface element, so every triangle keeps its owner count).
function surface_owner_counts(state)
    cells = [state.mesh.tetrahedra[k] for k in eachindex(state.mesh.tetrahedra) if !(k in state.mesh.dead)]
    return [count(cell -> all(v in cell for v in t), cells) for t in state.triangles]
end

# The wedge-spanning cells of the entity owning the edge (l1, l2) after the pass.
function spanning_after(state, l1, l2)
    masks = sidewall_vertex_masks(state.points, state.triangles, state.center, state.sides, state.grading.radius)
    return count(k -> !(k in state.mesh.dead) && l2 in state.mesh.tetrahedra[k] &&
                      wedge_spanning_cell(state.points, state.mesh.tetrahedra[k], masks) !== nothing,
                 state.mesh.incident[l1])
end

@testset "WEDGE_SPLIT_RULE is the design's 3.7 sentence, recorded only where it fired (decision 670)" begin
    @test startswith(WEDGE_SPLIT_RULE, "When the reconnection rounds leave a wedge-SPANNING bridging sliver")
    @test occursin("sub-size collapse remapped (erratum S-1)", WEDGE_SPLIT_RULE)
    @test occursin("up to WedgeSplitRounds (2) per corner", WEDGE_SPLIT_RULE)
    @test WEDGE_SPLIT_ROUNDS == 2 && WEDGE_SPLIT_RADIUS_STEPS == 8 && WEDGE_SPLIT_Z_STEPS == 8
    # The F5-B rule's census leaf is byte-unchanged by the design (decision 670 (A)).
    @test !occursin("wedge-SPANNING", CORNER_RECONNECTION_RULE) && !occursin("WedgeSplit", CORNER_RECONNECTION_RULE)
    @test endswith(CORNER_RECONNECTION_RULE, "decisions 363 / 365 / 392)")
end

@testset "7fd4 fab corner 8: the real cell 4276470 (DESIGN 6.1)" begin
    state = star_state("7fd4-fab-corner8-failed-star.json")
    tetrahedra = state.mesh.tetrahedra
    @test length(tetrahedra) == 399 && size(state.points, 2) == 122
    @test all(tetrahedron_signed_volume(state.points, cell) > 0.0 for cell in tetrahedra)
    # The trigger: star cell 377 (the design's 0-based 376), kappa_reg 5.505949976910859 = the
    # build log to the digit, the dihedral at O-c = the 151.1887-degree wall angle.
    @test state.trigger == [377]
    k_t = 377
    @test kappa(state.points, tetrahedra[k_t]) ≈ 5.505949976910859 atol=1e-12
    span = wedge_spanning_cell(state.points, tetrahedra[k_t], state.masks)
    @test span !== nothing
    @test span.line == (1, 13) && Set(span.walls) == Set([72, 59])
    dihedral = rad2deg(tetrahedron_edge_dihedral(state.points, span.line..., span.walls...))
    wall_angle = acosd(clamp(dot(state.sides.walls[1], state.sides.walls[2]), -1.0, 1.0))
    @test dihedral ≈ wall_angle atol=1e-6
    @test dihedral ≈ 151.1887 atol=1e-3
    @test kappa_reg_dihedral_floor(deg2rad(dihedral)) ≈ 5.5058 atol=2e-4
    # F5-B's replay: no 2-3 flip; the edge removal of a-b (a 4-cell ring) is REFUSED by the
    # improvement guard (its best re-fan 5.7815 > 5.50595), every other guard passed.
    entity_of(k) = state.mesh.entity_of[k]
    surface_edges = Set(minmax(t[i], t[mod1(i + 1, 3)]) for t in state.triangles for i in 1:3)
    flip = best_two_three_flip(state.points, tetrahedra, state.mesh.incident, state.surface_faces, entity_of, k_t,
                               state.mesh.floors[k_t], state.mesh.condition_ceilings[k_t])
    removal = best_edge_removal(state.points, tetrahedra, state.mesh.incident, surface_edges, entity_of, k_t,
                                state.mesh.floors[k_t], state.mesh.condition_ceilings[k_t])
    @test flip === nothing
    @test removal !== nothing && length(removal.ring) == 4 && Set(removal.edge) == Set([59, 72])
    @test removal.quality ≈ 5.7815415177843645 atol=1e-9
    @test !(removal.quality < kappa(state.points, tetrahedra[k_t]))
    owners_before = surface_owner_counts(state)
    cells_before = length(tetrahedra)
    record, applied, after, moves = run_wedge_split!(state)
    # The candidate census and the selection of record (DESIGN 6.1): 112 candidates, 58 valid;
    # SELECTED minimal, r = 5 / 8 x CornerSize, z = 3 / 4: cavity 4 -> 8, shortest new edge
    # 0.171 nm, at the insertion corner 2.5054 / new cells 20.42, after the re-descent 3.0832.
    @test applied
    @test record["Candidates"] == 112 && record["Valid"] == 58
    @test length(record["TriggerCells"]) == 1
    @test record["TriggerCells"][1]["KappaReg"] ≈ 5.505949976910859 atol=1e-12
    @test sort(record["TriggerCells"][1]["Masks"]) == [1, 2, 3, 3]
    selected = record["Selected"]
    @test selected["Form"] == "minimal" && selected["RadiusOverCornerSize"] == 0.625 && selected["ZFraction"] == 0.75
    @test selected["CavityCells"] == 4 && selected["NewCells"] == 8
    @test selected["ShortestNewEdge"] ≈ 0.171e-3 atol=1e-6
    @test selected["KappaRegAtInsertion"]["Corner"] ≈ 2.5054 atol=1e-4
    @test selected["KappaRegAtInsertion"]["NewCells"] ≈ 20.4205 atol=1e-3
    @test selected["Score"] ≈ 3.08322992566899 atol=1e-9
    @test after ≈ 3.08322992566899 atol=1e-9
    @test selected["NewCellsKappaRegMaxAfter"] ≈ 3.08322992566899 atol=1e-9
    @test selected["MinScaledJacobianNew"] ≈ 0.238 atol=2e-3
    @test selected["MinScaledJacobianNewAfter"] ≈ 0.2885 atol=1e-3
    @test selected["VolumeConservation"] <= 1.0e-9
    @test moves == 94
    @test selected["Displacement"] ≈ 0.132e-3 atol=2e-6
    # The rejections of record: 44 below the edge floor, 7 no improvement, 2 the trigger's own
    # face not visible, 1 the scaled-Jacobian floor.
    rejections = record["Rejections"]
    @test rejections["a new edge below the sub-size collapse threshold"] == 44
    @test rejections["no improvement on the trigger"] == 7
    @test rejections["the trigger's own boundary face is not visible from s"] == 2
    @test rejections["the scaled-Jacobian floor"] == 1
    @test sum(values(rejections)) + record["Valid"] == 112
    # The mesh after: the cavity's 4 slots reused + 4 appended (no dead slot), every cell
    # positive, every surface triangle keeps its owner count (no surface element changed),
    # the corner-line edge O-c owned by 2 cells of the entity and 0 wedge-spanning cells.
    @test length(tetrahedra) == cells_before + 4 && isempty(state.mesh.dead)
    @test length(record["NewCellSlots"]) == 8 && record["CavityCells"] == sort(record["CavityCells"])
    @test k_t in record["CavityCells"]
    @test all(tetrahedron_signed_volume(state.points, cell) > 0.0 for cell in tetrahedra)
    @test surface_owner_counts(state) == owners_before
    @test all(state.s_index in tetrahedra[p] for p in record["NewCellSlots"])
    @test all(entity_of(p) == 61 for p in record["NewCellSlots"])
    l1, l2 = span.line
    @test count(k -> l2 in tetrahedra[k] && entity_of(k) == 61, state.mesh.incident[l1]) == 2
    @test spanning_after(state, l1, l2) == 0
    @test maximum(kappa(state.points, tetrahedra[k]) for k in state.corner_incident(state.center)) ≈ 3.08322992566899 atol=1e-9
    # Bookkeeping as F5-B's: the 4 replaced cells by their original index, 8 added slots.
    @test sort(state.mesh.reconnection_replaced) == record["CavityCells"]
    @test sort(collect(keys(state.mesh.reconnection_added))) == sort(record["NewCellSlots"])
    @test all(state.mesh.origin[p] == 0 for p in record["NewCellSlots"])
    # Decision 34: the inserted vertex carries the edge floor and the corner law's bound.
    @test state.edge_floors[state.s_index] == state.collapse_threshold
    @test state.bounds[state.s_index] > 0.0
    @test norm(state.points[:, state.s_index] .- state.original[:, state.s_index]) <= state.bounds[state.s_index] * (1 + 1e-12)
end

@testset "synthetic 151.19-degree theta-0 fab, F5-B off: cell 106 (DESIGN 5.2)" begin
    state = star_state("tip-151.19-th0.0-fab-noF5B-corner1-failed-star.json")
    tetrahedra = state.mesh.tetrahedra
    @test length(tetrahedra) == 427
    @test state.trigger == [107]
    k_t = 107
    @test kappa(state.points, tetrahedra[k_t]) ≈ 5.506049923945327 atol=1e-12
    span = wedge_spanning_cell(state.points, tetrahedra[k_t], state.masks)
    @test span !== nothing && span.line == (1, 2) && Set(span.walls) == Set([34, 62])
    # On this star F5-B was ABSENT, not refused: its edge removal of a-b is LEGAL (ring 3, 2.0806).
    entity_of(k) = state.mesh.entity_of[k]
    surface_edges = Set(minmax(t[i], t[mod1(i + 1, 3)]) for t in state.triangles for i in 1:3)
    removal = best_edge_removal(state.points, tetrahedra, state.mesh.incident, surface_edges, entity_of, k_t,
                                state.mesh.floors[k_t], state.mesh.condition_ceilings[k_t])
    @test removal !== nothing && length(removal.ring) == 3
    @test removal.quality ≈ 2.0806137348227107 atol=1e-9
    @test best_two_three_flip(state.points, tetrahedra, state.mesh.incident, state.surface_faces, entity_of, k_t,
                              state.mesh.floors[k_t], state.mesh.condition_ceilings[k_t]) === nothing
    owners_before = surface_owner_counts(state)
    record, applied, after, moves = run_wedge_split!(state)
    # 112 candidates, 30 valid; SELECTED delaunay r = CornerSize / 2, z = 5 / 8: cavity 18 -> 24,
    # shortest new edge 0.155 nm, at the insertion 2.5147 / 7.92, after the re-descent 3.2927
    # (corner) and 3.2923 (new cells), 255 moves.
    @test applied
    @test record["Candidates"] == 112 && record["Valid"] == 30
    selected = record["Selected"]
    @test selected["Form"] == "delaunay" && selected["RadiusOverCornerSize"] == 0.5 && selected["ZFraction"] == 0.625
    @test selected["CavityCells"] == 18 && selected["NewCells"] == 24
    @test selected["ShortestNewEdge"] ≈ 0.155e-3 atol=1e-6
    @test selected["KappaRegAtInsertion"]["Corner"] ≈ 2.5147 atol=1e-4
    @test selected["KappaRegAtInsertion"]["NewCells"] ≈ 7.9228 atol=1e-3
    @test selected["Score"] ≈ 3.2926843126450995 atol=1e-9
    @test selected["NewCellsKappaRegMaxAfter"] ≈ 3.2923174534878346 atol=1e-9
    @test selected["MinScaledJacobianNewAfter"] ≈ 0.2654 atol=1e-3
    @test moves == 255
    rejections = record["Rejections"]
    @test rejections["a new edge below the sub-size collapse threshold"] == 60
    @test rejections["the trigger cell is not in the Delaunay cavity of s"] == 16
    @test rejections["no improvement on the trigger"] == 4
    @test sum(values(rejections)) + record["Valid"] == 112
    @test all(tetrahedron_signed_volume(state.points, cell) > 0.0 for k in eachindex(tetrahedra) if !(k in state.mesh.dead) for cell in (tetrahedra[k],))
    @test surface_owner_counts(state) == owners_before
    l1, l2 = span.line
    @test count(k -> !(k in state.mesh.dead) && l2 in tetrahedra[k] && entity_of(k) == 9, state.mesh.incident[l1]) == 2
    @test spanning_after(state, l1, l2) == 0
    @test after ≈ 3.2926843126450995 atol=1e-9
end

@testset "no trigger: the passed synthetic star and 7fd4's corner 1 are no-ops" begin
    for name in ("tip-151.19-th0.0-fab-corner1-passed-star.json", "7fd4-fab-corner1-control-star.json")
        state = star_state(name)
        @test isempty(state.trigger)
        @test maximum(kappa(state.points, state.mesh.tetrahedra[k]) for k in state.targets) <= state.gate
        points_before = copy(state.points); cells_before = copy(state.mesh.tetrahedra)
        record, applied, after, moves = run_wedge_split!(state)
        @test !applied && record["Candidates"] == 0 && record["Valid"] == 0 && isempty(record["TriggerCells"])
        @test record["Selected"] === nothing && moves == 0
        @test state.points == points_before && state.mesh.tetrahedra == cells_before
        @test isempty(state.mesh.reconnection_replaced) && isempty(state.mesh.reconnection_added) && isempty(state.mesh.dead)
    end
end

@testset "erratum S-1 (vii): a collapse-remapped cell in the cavity rejects the candidate" begin
    # The selected cavity of record on 7fd4 holds the trigger and its three face neighbours.
    reference = star_state("7fd4-fab-corner8-failed-star.json")
    record0, _, _, _ = run_wedge_split!(reference)
    cavity = record0["CavityCells"]
    neighbour = first(k for k in cavity if k != 377)
    # (a) that neighbour remapped: every candidate whose cavity holds it is rejected with the
    # S-1 reason; the selection changes or vanishes, never silently keeps the cavity.
    state = star_state("7fd4-fab-corner8-failed-star.json")
    record, applied, after, moves = run_wedge_split!(state; remapped=Set([neighbour]))
    reason = "a cavity cell the sub-size collapse remapped (erratum S-1 (vii))"
    @test haskey(record["Rejections"], reason) && record["Rejections"][reason] >= 1
    @test record["Valid"] < record0["Valid"]
    @test record["Candidates"] == 112
    if applied
        @test !(neighbour in record["CavityCells"])
        @test record["Selected"]["Score"] >= record0["Selected"]["Score"] - 1e-12
    end
    # (b) the trigger itself remapped: treated as F5-B treats it (skipped), nothing inserted.
    state = star_state("7fd4-fab-corner8-failed-star.json")
    record, applied, _, _ = run_wedge_split!(state; remapped=Set([377]))
    @test !applied && isempty(record["TriggerCells"]) && record["Candidates"] == 0
    @test length(state.mesh.tetrahedra) == 399
end

@testset "cavity forms on the stored star: conformity of every valid candidate" begin
    # Every candidate of the grid whose cavity exists: the children carry their parents' sign,
    # the volumes sum to the cavity's, the cavity is face-connected and single-entity.
    state = star_state("7fd4-fab-corner8-failed-star.json")
    tetrahedra = state.mesh.tetrahedra
    entity_of(k) = state.mesh.entity_of[k]
    span = wedge_spanning_cell(state.points, tetrahedra[377], state.masks)
    l1, l2 = span.line
    bisector = state.sides.walls[1] .+ state.sides.walls[2]
    bisector = vcat(bisector ./ norm(bisector), 0.0)
    cavities = 0
    for form in WEDGE_SPLIT_FORMS, j in 1:WEDGE_SPLIT_RADIUS_STEPS, i in 1:(WEDGE_SPLIT_Z_STEPS - 1)
        s = state.points[:, l1] .+ (i / 8) .* (state.points[:, l2] .- state.points[:, l1]) .+
            (j / 8 * state.grading.corner_size) .* bisector
        state.points[:, state.s_index] = s
        new_cells, parents, cavity, reason = wedge_split_cavity(state.points, state.s_index, tetrahedra,
                                                                state.mesh.incident, state.surface_faces,
                                                                entity_of, 377, form)
        new_cells === nothing && continue
        cavities += 1
        @test 377 in cavity && length(new_cells) == length(parents)
        @test all(entity_of(k) == 61 for k in cavity)
        @test all(sign(tetrahedron_signed_volume(state.points, c)) == sign(tetrahedron_signed_volume(state.points, tetrahedra[p]))
                  for (c, p) in zip(new_cells, parents))
        @test all(tetrahedron_signed_volume(state.points, c) != 0.0 for c in new_cells)
        volume(c) = abs(tetrahedron_signed_volume(state.points, c))
        @test abs(sum(volume(c) for c in new_cells) - sum(volume(tetrahedra[k]) for k in cavity)) <=
              1e-9 * sum(volume(tetrahedra[k]) for k in cavity)
        # every boundary face of the cavity is a face of exactly one child; every child holds s
        @test all(state.s_index in c for c in new_cells)
        @test length(unique(new_cells)) == length(new_cells)
    end
    @test cavities >= 58
end
