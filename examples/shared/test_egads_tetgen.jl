#=
    test_egads_tetgen.jl

Tests of egads_tetgen.jl (moved from EGADS.jl's test_tetgen.jl): PLC construction from BODY,
direct tetrahedralisation, face-marker preservation, multi-body region
seeding, and shared-interface stitching.
=#

using Test
using EGADS
using TetGen
include(joinpath(@__DIR__, "egads_tetgen.jl"))
using .EGADSTetGen

@testset "EGADS → TetGen integration" begin
    @testset "single-body PLC: box (triangulated)" begin
        ctx = EGADS.EGADSContext()
        try
            box = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)

            io = EGADSTetGen.tetgen_io(box; angle = 20.0)
            @test size(io.pointlist, 1) == 3
            @test size(io.pointlist, 2) == 8
            @test length(io.facetlist) == 12                  # 6 faces × 2 tris
            @test length(io.facetmarkerlist) == 12
            @test sort(unique(io.facetmarkerlist)) == Int32.(1:6)
            @test count(==(Int32(1)), io.facetmarkerlist) == 2
        finally
            close(ctx)
        end
    end

    @testset "direct tetrahedralize(body, flags)" begin
        ctx = EGADS.EGADSContext()
        try
            box = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            out = EGADSTetGen.tetrahedralize(box, "pq1.414Qa0.01")
            @test size(out.pointlist, 1) == 3
            @test size(out.tetrahedronlist, 1) == 4
            @test size(out.tetrahedronlist, 2) > 0
            @test sort(unique(out.trifacemarkerlist)) == Int32.(1:6)
        finally
            close(ctx)
        end
    end

    @testset "face_marker_offset (single body)" begin
        ctx = EGADS.EGADSContext()
        try
            box = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            io = EGADSTetGen.tetgen_io(box; angle = 30.0, face_marker_offset = 500)
            @test sort(unique(io.facetmarkerlist)) == Int32.(501:506)
        finally
            close(ctx)
        end
    end

    @testset "multi-body disjoint: two boxes (polygon path)" begin
        ctx = EGADS.EGADSContext()
        try
            b1 = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            b2 = EGADS.make_box(ctx, 2.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            bodies = [b1, b2]

            io = EGADSTetGen.tetgen_io(bodies; angle = 30.0)
            # Two disjoint boxes, all 6 faces per body emitted as polygon
            # facets — 12 total; planar → each facet has exactly one 4-vertex
            # polygon.
            @test length(io.facetlist) == 12
            @test all(length(f.polygonlist) == 1 for f in io.facetlist)
            @test all(length(f.polygonlist[1]) == 4 for f in io.facetlist)

            u = sort(unique(io.facetmarkerlist))
            @test u == vcat(Int32.(1:6), Int32.(1001:1006))

            # One region seed per body.
            @test size(io.regionlist, 2) == 2
            @test sort(io.regionlist[4, :]) == [1.0, 2.0]
            @test EGADS.in_topology(b1, io.regionlist[1:3, 1])
            @test EGADS.in_topology(b2, io.regionlist[1:3, 2])

            out = EGADSTetGen.tetrahedralize(bodies, "pq1.414AQa")
            @test size(out.tetrahedronlist, 2) > 0
            @test size(out.tetrahedronattributelist, 1) == 1
            @test sort(unique(out.tetrahedronattributelist)) == [1.0, 2.0]
        finally
            close(ctx)
        end
    end

    @testset "multi-body stitched: two abutting boxes (shared face)" begin
        # Two boxes sharing their x=1 face. With polygon-facet emission +
        # point / polygon dedup, the shared interface collapses to one PLC
        # facet and TetGen produces a single conformal multi-region mesh.
        ctx = EGADS.EGADSContext()
        try
            b1 = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            b2 = EGADS.make_box(ctx, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0)
            bodies = [b1, b2]

            io = EGADSTetGen.tetgen_io(bodies; angle = 30.0)
            # 12 unique points (8 per box minus 4 shared corners)
            @test size(io.pointlist, 2) == 12
            # 11 unique facets (12 faces − 1 collapsed interface)
            @test length(io.facetlist) == 11

            @test EGADS.in_topology(b1, io.regionlist[1:3, 1])
            @test EGADS.in_topology(b2, io.regionlist[1:3, 2])

            out = TetGen.tetrahedralize(io, "pq1.414AQa0.05")
            @test size(out.tetrahedronlist, 2) > 0
            attrs = sort(unique(out.tetrahedronattributelist))
            @test attrs == [1.0, 2.0]
            n1 = count(==(1.0), out.tetrahedronattributelist)
            n2 = count(==(2.0), out.tetrahedronattributelist)
            @test n1 > 0 && n2 > 0
            @test 0.5 < n1 / n2 < 2.0   # identical box volumes → similar tet counts
        finally
            close(ctx)
        end
    end

    @testset "hole_seeds: box with interior cavity" begin
        ctx = EGADS.EGADSContext()
        try
            outer = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 4.0, 4.0, 4.0)
            inner = EGADS.make_box(ctx, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0)
            shell_bodies = EGADS.solid_boolean(outer, inner, EGADS.BoolOp.SUBTRACTION)
            @test length(shell_bodies) == 1
            shell_body = first(shell_bodies)

            io = EGADSTetGen.tetgen_io(shell_body; angle = 30.0)
            # Triangulated single-body path: 12 faces × 2 tris = 24 facets
            @test length(sort(unique(io.facetmarkerlist))) == 12

            io.holelist = reshape(Cdouble[2.0, 2.0, 2.0], 3, 1)
            out = TetGen.tetrahedralize(io, "pq1.414Qa0.5")
            @test size(out.tetrahedronlist, 2) > 0
            @test length(unique(out.trifacemarkerlist)) == 12
        finally
            close(ctx)
        end
    end

    @testset "tetrahedralize_size_field: graded field (:callback)" begin
        ctx = EGADS.EGADSContext()
        try
            box = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0)
            # Fine target near the (0,0,0) corner, coarse far away — a
            # gmsh-style graded edge-length field.
            h(x, y, z) = 0.5 + 0.4 * sqrt(x^2 + y^2 + z^2)

            out = EGADSTetGen.tetrahedralize_size_field([box], "pAYQ";
                meshsize_callback = h, method = :callback,
                max_iterations = 5, verbose = false)
            ntets = size(out.tetrahedronlist, 2)
            @test ntets > 0

            # Every tet respects the field: longest edge ≤ h(centroid)
            # (allow numerical slack on the comparison itself).
            unsuit = EGADSTetGen._make_unsuitable(h, 1.0)
            nbad = EGADSTetGen._count_unsuitable(out, unsuit)
            @test nbad == 0

            # Grading is real: mean edge length near the fine corner is
            # smaller than near the coarse corner.
            pts = out.pointlist
            tets = out.tetrahedronlist
            near_sum = 0.0; near_n = 0
            far_sum = 0.0; far_n = 0
            for j in 1:ntets
                a, b = tets[1, j], tets[2, j]
                e = sqrt((pts[1, a] - pts[1, b])^2 + (pts[2, a] - pts[2, b])^2 +
                         (pts[3, a] - pts[3, b])^2)
                cx = sum(pts[1, tets[i, j]] for i in 1:4) / 4
                cy = sum(pts[2, tets[i, j]] for i in 1:4) / 4
                cz = sum(pts[3, tets[i, j]] for i in 1:4) / 4
                r = sqrt(cx^2 + cy^2 + cz^2)
                if r < 3.0
                    near_sum += e; near_n += 1
                elseif r > 12.0
                    far_sum += e; far_n += 1
                end
            end
            @test near_n > 0 && far_n > 0
            @test (near_sum / near_n) < 0.6 * (far_sum / far_n)
        finally
            close(ctx)
        end
    end

    @testset "tetrahedralize_size_field: legacy :metric path still runs" begin
        ctx = EGADS.EGADSContext()
        try
            box = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0)
            h(x, y, z) = 1.0 + 0.5 * sqrt(x^2 + y^2 + z^2)
            out = EGADSTetGen.tetrahedralize_size_field([box], "pAYQ";
                meshsize_callback = h, method = :metric,
                max_iterations = 3, tolerance = 5.0, verbose = false)
            @test size(out.tetrahedronlist, 2) > 0
        finally
            close(ctx)
        end
    end

    @testset "classify_tetgen_vertices: box (planar exactness)" begin
        ctx = EGADS.EGADSContext()
        try
            box = EGADS.make_box(ctx, 0.0, 0.0, 0.0, 10.0, 10.0, 10.0)
            # Force facet/segment Steiner points so all three classes appear.
            out = EGADSTetGen.tetrahedralize_size_field([box], "pAYQ";
                meshsize_callback = (x, y, z) -> 3.0, method = :callback,
                max_iterations = 3, verbose = false)
            cls = EGADSTetGen.classify_tetgen_vertices(out, box)

            # Full coverage: every constrained vertex classified.
            @test cls.n_unclassified == 0
            @test cls.n_constrained ==
                  size(cls.on_nodes, 1) + size(cls.on_edges, 1) + size(cls.on_faces, 1)
            # All 8 corners found as geometric NODEs.
            @test size(cls.on_nodes, 1) == 8
            @test sort(cls.on_nodes[:, 2]) == collect(1:8)
            # Box edges and faces both populated (Steiner points forced by sizing).
            @test size(cls.on_edges, 1) > 0
            @test size(cls.on_faces, 1) > 0
            @test all(1 .<= cls.on_edges[:, 2] .<= 12)
            @test all(1 .<= cls.on_faces[:, 2] .<= 6)
            # Planar geometry: every classification distance is exact.
            @test cls.max_dist < 1e-9
            # FACE (u,v) round-trip: evaluating the face at (u,v) recovers the vertex.
            faces = EGADS.get_body_topos(box, EGADS.OClass.FACE)
            for r in eachrow(cls.on_faces[1:min(end, 20), :])
                vi, fid = Int(r[1]), Int(r[2])
                ev = EGADS.evaluate(faces[fid], [r[3], r[4]])
                p = out.pointlist[:, vi]
                @test hypot((ev.xyz .- p)...) < 1e-9
            end
        finally
            close(ctx)
        end
    end

    @testset "classify_tetgen_vertices: cylinder (curved chord error)" begin
        ctx = EGADS.EGADSContext()
        try
            cyl = EGADS.make_cylinder(ctx, [0.0, 0.0, 0.0], [0.0, 0.0, 5.0], 2.0)
            out = EGADSTetGen.tetrahedralize_size_field([cyl], "pAYQ";
                meshsize_callback = (x, y, z) -> 1.0, method = :callback,
                max_iterations = 3, angle = 30.0, verbose = false)
            cls = EGADSTetGen.classify_tetgen_vertices(out, cyl)

            @test cls.n_unclassified == 0
            @test size(cls.on_faces, 1) > 0
            # Chord error bound: tessellation at angle=30° gives sagitta
            # ≤ r·(1−cos(15°)) ≈ 0.0341·r; TetGen Steiner points on the linear
            # facets cannot exceed it (with slack for refinement subdivision).
            sag = 2.0 * (1 - cosd(15.0))
            @test cls.max_dist < 2 * sag
            # Curved-side vertices: projecting along (u,v) must recover a point
            # on the true cylinder (radius 2 exactly). Side faces are TRIMMED
            # surfaces wrapping the cylinder, so select by "not planar".
            faces = EGADS.get_body_topos(cyl, EGADS.OClass.FACE)
            n_curved = 0
            for r in eachrow(cls.on_faces)
                fid = Int(r[2])
                gi = EGADS.get_info(EGADS.get_topology(faces[fid]).geom)
                gi.mtype == EGADS.MType.PLANE && continue
                ev = EGADS.evaluate(faces[fid], [r[3], r[4]])
                @test hypot(ev.xyz[1], ev.xyz[2]) ≈ 2.0 atol = 1e-9
                n_curved += 1
                n_curved >= 25 && break
            end
            @test n_curved > 0
        finally
            close(ctx)
        end
    end

    @testset "curve_tetgen_mesh: cylinder → valid curved P2" begin
        ctx = EGADS.EGADSContext()
        try
            cyl = EGADS.make_cylinder(ctx, [0.0, 0.0, 0.0], [0.0, 0.0, 5.0], 2.0)
            out = EGADSTetGen.tetrahedralize_size_field([cyl], "pAYQ";
                meshsize_callback = (x, y, z) -> 1.0, method = :callback,
                max_iterations = 3, angle = 30.0, verbose = false)
            cls = EGADSTetGen.classify_tetgen_vertices(out, cyl)
            curved = EGADSTetGen.curve_tetgen_mesh(out, cls, cyl)

            # Continuation reaches the fully curved configuration, and the
            # curved mesh is valid by the sampled-detJ gate.
            @test curved.t_final == 1.0
            @test curved.min_detj_ratio > 0.05
            @test curved.n_curved_mid > 0
            @test curved.n_snapped > 0

            nv = size(out.pointlist, 2)
            # Snapped chord vertices now lie on the true cylinder: every
            # side-face-classified vertex has radius 2.
            for r in eachrow(cls.on_faces)
                vi = Int(r[1])
                fid = Int(r[2])
                gi = EGADS.get_info(EGADS.get_topology(
                    EGADS.get_body_topos(cyl, EGADS.OClass.FACE)[fid]).geom)
                gi.mtype == EGADS.MType.PLANE && continue
                @test hypot(curved.coords[1, vi], curved.coords[2, vi]) ≈ 2.0 atol = 1e-7
            end
            # Boundary tri mid-nodes on the curved side: radius 2 as well.
            # (Find them via tris6 rows whose corners are all side-classified.)
            side_verts = Set{Int}(Int(r[1]) for r in eachrow(cls.on_faces)
                                  if begin
                                      f = EGADS.get_body_topos(cyl, EGADS.OClass.FACE)[Int(r[2])]
                                      EGADS.get_info(EGADS.get_topology(f).geom).mtype !=
                                      EGADS.MType.PLANE
                                  end)
            n_checked = 0
            for j in 1:size(curved.tris6, 2)
                all(curved.tris6[k, j] in side_verts for k in 1:3) || continue
                for k in 4:6
                    m = curved.tris6[k, j]
                    @test hypot(curved.coords[1, m], curved.coords[2, m]) ≈ 2.0 atol = 1e-7
                    n_checked += 1
                end
                n_checked >= 12 && break
            end
            @test n_checked > 0

            # Writer: order-2 gmsh file with the expected counts.
            path = tempname() * ".msh2"
            EGADSTetGen.write_msh2_p2(path, curved;
                tet_attr = 1, tri_attr = m -> m,
                physical_names = [(3, 1, "volume")])
            txt = read(path, String)
            @test occursin("\$MeshFormat\n2.2 0 8", txt)
            @test occursin("\$Nodes\n$(size(curved.coords, 2))\n", txt)
            @test count(" 11 2 ", txt) == size(curved.tets10, 2)
            @test count(" 9 2 ", txt) == size(curved.tris6, 2)
            rm(path; force = true)
        finally
            close(ctx)
        end
    end
end
