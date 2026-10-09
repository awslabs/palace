# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "prism_edge_tubes.jl"))

# The BOUND RING of the fabricated arc-side tube sections (round 4 F-3 redesign, decisions 629 /
# 630; DESIGN.md sections 4 / 8, the lane's DESIGN AMENDMENT): covered rays / radii / sector
# kinds / pivots at the production process, the two cell forms (A: prism + pyramid + tet and the
# shoulder pyramid + tet; B: the centroid split of an unpaired sector) with their conformity,
# the side-pair node map across z = T / 2, and the unchanged sections without a bound ring.

# The production section pair (inner 0.25 nm, ratio 2, 7 rings, 30-degree sectors; Overetch
# 0.05, MetalThickness 0.1).
production_bottom(; planes=[-0.05, 0.05]) = TubeSection(
    0.00025, 2.0, 7, [180.0 + 30.0 * j for j = 0:9], vcat(fill(1, 3), fill(2, 6)); bound_planes=planes)
production_top() = TubeSection(0.00025, 2.0, 7, [-90.0 + 30.0 * j for j = 0:9], fill(2, 9);
                               bound_planes=[-0.05])

# Every face of a cell as a sorted node tuple with its multiplicity sign: a face is interior
# when two cells share it (as a whole), boundary when one does. Triangles are compared as
# sets of 3 nodes; quads as sets of 4.
function cell_faces(kind, nodes)
    if kind == :tet
        a, b, c, d = nodes
        return [Set([a, b, c]), Set([a, b, d]), Set([a, c, d]), Set([b, c, d])]
    elseif kind == :pyramid
        a, b, c, d, e = nodes
        return [Set([a, b, c, d]), Set([a, b, e]), Set([b, c, e]), Set([c, d, e]), Set([d, a, e])]
    else
        a, b, c, d, e, f = nodes
        return [Set([a, b, c]), Set([d, e, f]), Set([a, b, e, d]), Set([b, c, f, e]), Set([c, a, d, f])]
    end
end

# The 3D coordinates of a bound-cell node reference over one straight layer of spacing lc
# (u, w, s): the centroid node is the mean of the hex corners.
function reference_point(section, uw, j, reference, lc)
    local_index, offset = reference
    if local_index == BOUND_CENTROID_NODE
        a, d = bound_ray_nodes(section, j)
        b, c = bound_ray_nodes(section, j + 1)
        corners = [[uw[1, n], uw[2, n], o * lc] for n in (a, b, c, d) for o in (0.0, 1.0)]
        return sum(corners) ./ 8
    end
    return [uw[1, local_index], uw[2, local_index], offset * lc]
end

@testset "bound ring: covered rays, radii, kinds and pivots at the production process (DESIGN 8.2)" begin
    bottom = production_bottom()
    top = production_top()
    @test has_bound_ring(bottom) && has_bound_ring(top)
    @test !has_bound_ring(TubeSection(0.00025, 2.0, 7, bottom.angles, bottom.materials))
    record = bound_ring_record(bottom)
    @test record["Rays"] == [240.0, 270.0, 300.0, 420.0, 450.0]
    @test isapprox(record["Radii"], [0.05 / sind(60.0), 0.05, 0.05 / sind(60.0), 0.05 / sind(60.0), 0.05]; rtol=1.0e-12)
    @test record["Planes"] == [-0.05, -0.05, -0.05, 0.05, 0.05]
    @test record["CoveredSectors"] == [[240.0, 270.0], [270.0, 300.0], [420.0, 450.0]]
    @test record["ShoulderSectors"] == [[210.0, 240.0], [300.0, 330.0], [390.0, 420.0]]
    @test record["Pivots"] == [240.0, 240.0, 300.0, 300.0, 420.0, 420.0]
    @test isempty(record["UnpairedSectors"])
    @test record["Reach"] == 0.03175 + 2.0 * 0.016 && record["MaximumRadius"] == 0.05 / sind(60.0)
    @test section_reach(bottom) == 0.05 / sind(60.0)
    top_record = bound_ring_record(top)
    @test top_record["Rays"] == [-90.0, -60.0] && top_record["CoveredSectors"] == [[-90.0, -60.0]] &&
          top_record["ShoulderSectors"] == [[-60.0, -30.0]] && top_record["Pivots"] == [-60.0, -60.0]
    # The bound nodes sit ON their planes (DESIGN 8.2 geometry: bound240 (-0.02887, -0.05)).
    uw = section_coordinates(bottom)
    @test size(uw, 2) == section_node_count(bottom) == 1 + 7 * 10 + 5
    @test isapprox(uw[:, bound_node(bottom, 2)], [-0.05 / tand(60.0), -0.05]; atol=1.0e-12)
    @test uw[:, bound_node(bottom, 3)] == [0.0, -0.05] || isapprox(uw[:, bound_node(bottom, 3)], [0.0, -0.05]; atol=1.0e-15)
    @test isapprox(uw[:, bound_node(bottom, 8)], [0.05 / tand(60.0), 0.05]; atol=1.0e-12)
    @test outer_node(bottom, 0) == section_node(bottom, 7, 0) && outer_node(bottom, 2) == bound_node(bottom, 2)
    @test section_node(bottom, 8, 2) == bound_node(bottom, 2)
    @test_throws ErrorException bound_node(bottom, 0)
    # The presentations: the shoulder rays split, the shared ray 270 and the sidewall ray 450 quads.
    @test [bound_ray_presentation(bottom, j) for j in (2, 3, 4, 8, 9)] == [:tri, :quad, :tri, :tri, :quad]
    @test [bound_ray_presentation(top, j) for j in (0, 1)] == [:quad, :tri]
    # The free sectors keep their pyramids; the kinds per sector.
    @test [sector_kind(bottom, j) for j = 0:8] ==
          [:free, :shoulder, :covered, :covered, :shoulder, :free, :free, :shoulder, :covered]
    @test [sector_kind(top, j) for j = 0:8] == [:covered, :shoulder, :free, :free, :free, :free, :free, :free, :free]
    # The thin sheet section and a plane inside ring K fail closed.
    sheet = TubeSection(0.002, 2.0, 3, [180.0 + 30.0 * j for j = 0:12], vcat(fill(1, 6), fill(2, 6)))
    @test sheet.closed && !has_bound_ring(sheet)
    @test_throws ErrorException bound_ring(sheet, [-0.05], 2.0)
    @test_throws ErrorException TubeSection(0.00025, 2.0, 7, bottom.angles, bottom.materials; bound_planes=[-0.02])
    # No ray within reach: an empty bound ring (every sector free).
    far = TubeSection(0.00025, 2.0, 7, bottom.angles, bottom.materials; bound_planes=[-0.5])
    @test has_bound_ring(far) && isempty(bound_ring_record(far)["Rays"]) &&
          all(sector_kind(far, j) == :free for j = 0:8) && section_reach(far) == 0.03175
end

@testset "bound ring: form (A) cells, the shoulder cells and their conformity (decision 630 (1))" begin
    bottom = production_bottom()
    uw = section_coordinates(bottom)
    lc = 0.04
    counts = Dict(:prism => 0, :pyramid => 0, :tet => 0)
    for j = 0:8
        cells = bound_sector_cells(bottom, j)
        kind = sector_kind(bottom, j)
        kind == :free && (@test isempty(cells); continue)
        kinds = [c[1] for c in cells]
        @test kinds == (kind == :covered ? [:prism, :pyramid, :tet] : [:pyramid, :tet])
        for (cell_kind, references) in cells
            counts[cell_kind] += 1
            nodes = [(j, r) for r in references]
            points = Dict(n => reference_point(bottom, uw, j, n[2], lc) for n in nodes)
            oriented = orient_cell(cell_kind, nodes, n -> points[n])
            # A positive volume after the orientation fix.
            x = [points[n] for n in oriented]
            if cell_kind == :tet
                @test dot(cross(x[2] - x[1], x[3] - x[1]), x[4] - x[1]) > 0.0
            elseif cell_kind == :pyramid
                @test dot(cross(x[2] - x[1], x[4] - x[1]), x[5] - x[1]) > 0.0
            else
                @test dot(cross(x[2] - x[1], x[3] - x[1]), x[4] - x[1]) > 0.0
            end
        end
    end
    # Per station per bottom section: +4 prisms, 8 pyramids (the free sectors lose theirs), 8 tets
    # across the whole arc side = bottom 3 covered + 3 shoulders, top 1 + 1 (DESIGN AMENDMENT (3)).
    @test counts == Dict(:prism => 3, :pyramid => 6, :tet => 6)
    top_counts = Dict(k => count(c -> c[1] == k, reduce(vcat, [bound_sector_cells(production_top(), j) for j = 0:8]))
                      for k in (:prism, :pyramid, :tet))
    @test top_counts == Dict(:prism => 1, :pyramid => 2, :tet => 2)
    # CONFORMITY over one layer of the whole section beyond ring K: every face of the bound cells
    # is either shared by exactly two cells (interior) or lies on the outer boundary - the plane /
    # slant triangles (bound_outer_triangles), the ring-K quads (the regular prisms' faces), the
    # station cross-section triangles (bound_cap_triangles), the metal-face bound segments
    # (bound_ray_faces of the first / last rays).
    faces = Dict{Set{Tuple{Int, Int}}, Int}()
    for j = 0:8, (kind, references) in bound_sector_cells(bottom, j)
        for face in cell_faces(kind, references)
            faces[face] = get(faces, face, 0) + 1
        end
    end
    boundary = Set{Set{Tuple{Int, Int}}}()
    K = ring_count(bottom)
    for j = 0:8
        sector_kind(bottom, j) == :free && continue
        for t in bound_outer_triangles(bottom, j)
            push!(boundary, Set(t))
        end
        for t in bound_cap_triangles(bottom, j), offset in (0, 1)
            push!(boundary, Set((n, offset) for n in t))
        end
        a, b = section_node(bottom, K, j), section_node(bottom, K, j + 1)
        push!(boundary, Set([(a, 0), (b, 0), (b, 1), (a, 1)]))
    end
    for r in (0, 9)
        covered_ray(bottom, r) || continue
        _, elements = bound_ray_faces(bottom, r)
        for e in elements
            push!(boundary, Set(e))
        end
    end
    @test all(m in (1, 2) for m in values(faces))
    @test all(m == 1 for (f, m) in faces if f in boundary)
    @test all(f in boundary for (f, m) in faces if m == 1)
    # The shoulder / covered interface quads are split on the same diagonal (bound^0, ring_K^1).
    for r in (2, 4, 8)
        _, elements = bound_ray_faces(bottom, r)
        ring, bnd = bound_ray_nodes(bottom, r)
        @test length(elements) == 2 && all((bnd, 0) in e && (ring, 1) in e for e in elements)
    end
    # The plane diagonal of 240-270 starts at bound240 (the pivot) at station 0 and ends at
    # bound270 at station 1; 270-300 mirrors it (DESIGN AMENDMENT (1)).
    t = bound_outer_triangles(bottom, 2)
    @test all((bound_node(bottom, 2), 0) in s && (bound_node(bottom, 3), 1) in s for s in (Set(t[1]), Set(t[2])))
    t = bound_outer_triangles(bottom, 3)
    @test all((bound_node(bottom, 4), 0) in s && (bound_node(bottom, 3), 1) in s for s in (Set(t[1]), Set(t[2])))
    # The quality bounds of the amendment (4): the shoulder pyramid's height (bound240 over the
    # ring-K quad 210-240) is the triangle's altitude 0.0251; the covered pyramid's 0.024.
    a, b, c = section_node(bottom, K, 1), section_node(bottom, K, 2), bound_node(bottom, 2)
    p, q, r = uw[:, a], uw[:, b], uw[:, c]
    altitude = abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / norm(q - p)
    @test 0.024 <= altitude <= 0.026
    # (the covered pyramid of 240-270: bound240 over the diagonal quad ring240 - bound270)
    p, q, r = uw[:, section_node(bottom, K, 2)], uw[:, bound_node(bottom, 3)], uw[:, bound_node(bottom, 2)]
    altitude = abs((q[1] - p[1]) * (r[2] - p[2]) - (q[2] - p[2]) * (r[1] - p[1])) / norm(q - p)
    @test 0.023 <= altitude <= 0.025
end

@testset "bound ring: the side pair across z = T / 2 and the shared mid-plane triangles (DESIGN 8.1)" begin
    bottom = production_bottom()
    top = production_top()
    node_map = bound_side_pair_map(top, bottom, 0.1)
    @test node_map == Dict(bound_node(top, 0) => bound_node(bottom, 9), bound_node(top, 1) => bound_node(bottom, 8))
    uw_top, uw_bottom = section_coordinates(top), section_coordinates(bottom)
    for (t, b) in node_map
        # (cos(-60 deg) and cos(420 deg) agree to roundoff; the mesh shares the node tags themselves)
        @test isapprox(uw_top[1, t], uw_bottom[1, b]; atol=1.0e-15) && uw_top[2, t] == -0.05 &&
              uw_bottom[2, b] == 0.05
    end
    # The two tubes' mid-plane triangles coincide: the bottom's diagonal (bound420^0, bound450^1)
    # is the top's (bound-60^0, bound-90^1) under the node map.
    mapped(t) = Set((get(node_map, n, n), o) for (n, o) in t)
    bottom_triangles = Set(Set(t) for t in bound_outer_triangles(bottom, 8))
    top_triangles = Set(mapped(t) for t in bound_outer_triangles(top, 0))
    @test bottom_triangles == top_triangles
    # A top without the mid plane or a bottom without a partner fails closed.
    @test_throws ErrorException bound_side_pair_map(top, production_bottom(; planes=[-0.05]), 0.1)
    @test_throws ErrorException bound_side_pair_map(TubeSection(0.00025, 2.0, 7, top.angles, top.materials), bottom, 0.1)
end

@testset "bound ring: form (B) for an unpaired sector (decision 630 (2); an odd run of three)" begin
    # 20-degree sectors: the floor rays 240 / 260 / 280 / 300 are within reach (sin >= 0.784), a
    # run of three covered sectors between the shoulders 220-240 and 300-320.
    odd = TubeSection(0.00025, 2.0, 7, [180.0 + 20.0 * j for j = 0:13], vcat(fill(1, 4), fill(2, 9));
                      bound_planes=[-0.05, 0.05])
    record = bound_ring_record(odd)
    @test record["CoveredSectors"] == [[240.0, 260.0], [260.0, 280.0], [280.0, 300.0], [420.0, 440.0]]
    @test record["UnpairedSectors"] == [[260.0, 280.0]]
    @test record["Pivots"] == [240.0, 240.0, 300.0, 300.0, 420.0, 420.0]
    middle = 4   # the sector 260-280 (0-based)
    @test sector_pivot(odd, middle) == -1 && sector_pivot(odd, 3) == 3 && sector_pivot(odd, 5) == 6
    # Its neighbours pivot at the shoulder rays 240 / 300, so they present QUADS at the shared
    # rays 260 / 280: the middle's side faces are pyramids (3 pyramids: the ring-K quad and the two
    # sides; 6 tets: the plane and the two station faces split).
    @test bound_ray_presentation(odd, 4) == :quad && bound_ray_presentation(odd, 5) == :quad
    cells = bound_sector_cells(odd, middle)
    @test count(c -> c[1] == :pyramid, cells) == 3 && count(c -> c[1] == :tet, cells) == 6 &&
          !any(c -> c[1] == :prism, cells)
    @test any(any(r == (BOUND_CENTROID_NODE, 0) for r in c[2]) for c in cells)
    # Conformity of the whole run (the three covered sectors and their shoulders) over one layer.
    uw = section_coordinates(odd)
    faces = Dict{Set, Int}()
    sectors = [j for j = 0:12 if sector_kind(odd, j) != :free]
    for j in sectors, (kind, references) in bound_sector_cells(odd, j)
        for face in cell_faces(kind, [(j, r) for r in references])
            faces[face] = get(faces, face, 0) + 1
        end
    end
    # (the centroid references are per sector; every other reference is shared across sectors)
    normalise(face) = Set(r[2][1] == BOUND_CENTROID_NODE ? r : (0, r[2]) for r in face)
    merged = Dict{Set, Int}()
    for (f, m) in faces
        merged[normalise(f)] = get(merged, normalise(f), 0) + m
    end
    boundary = Set{Set}()
    K = ring_count(odd)
    for j in sectors
        for t in bound_outer_triangles(odd, j)
            push!(boundary, Set((0, n) for n in t))
        end
        for t in bound_cap_triangles(odd, j), offset in (0, 1)
            push!(boundary, Set((0, (n, offset)) for n in t))
        end
        a, b = section_node(odd, K, j), section_node(odd, K, j + 1)
        push!(boundary, Set((0, n) for n in [(a, 0), (b, 0), (b, 1), (a, 1)]))
    end
    for r in (0, 13)
        covered_ray(odd, r) || continue
        for e in bound_ray_faces(odd, r)[2]
            push!(boundary, Set((0, n) for n in e))
        end
    end
    @test all(m in (1, 2) for m in values(merged))
    @test all(m == 1 for (f, m) in merged if f in boundary)
    @test all(f in boundary for (f, m) in merged if m == 1)
    # Every form-(B) cell has a positive volume about the centroid over a straight layer.
    for (kind, references) in cells
        points = Dict(r => reference_point(odd, uw, middle, r, 0.04) for r in references)
        oriented = orient_cell(kind, references, r -> points[r])
        x = [points[r] for r in oriented]
        @test (kind == :tet ? dot(cross(x[2] - x[1], x[3] - x[1]), x[4] - x[1]) :
               dot(cross(x[2] - x[1], x[4] - x[1]), x[5] - x[1])) > 0.0
    end
    # A 1-sector run between two shoulders is unpaired too (45-degree sectors reach only ray 270
    # with the floor at 0.05 x 1.1 reach ... a reach that admits 225 / 315 too gives a 2-run):
    # reach 0.0638 covers 270 only -> shoulders 225-270 and 270-315, no covered sector; a plane at
    # -0.035 covers 225 / 270 / 315 (0.0495 / 0.035 / 0.0495) -> an even run of two.
    coarse = TubeSection(0.00025, 2.0, 7, [180.0 + 45.0 * j for j = 0:6], vcat(fill(1, 2), fill(2, 4));
                         bound_planes=[-0.05, 0.05])
    @test bound_ring_record(coarse)["CoveredSectors"] == Vector{Float64}[] &&
          bound_ring_record(coarse)["ShoulderSectors"] == [[225.0, 270.0], [270.0, 315.0], [405.0, 450.0]]
    @test bound_ray_presentation(coarse, 2) == :tri
end
