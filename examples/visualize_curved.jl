# Generate a curved (P2) Gmsh MSH file from the EGADS resonator geometry.
# Mid-edge nodes are evaluated on the EGADS surface → true curvature preserved.

using DeviceLayout, DeviceLayout.PreferredUnits, Unitful, EGADS, Printf

DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

function main()
    # Build the resonator geometry (same as resonator_native_test.jl)
    r = 50.0μm
    style = Paths.SimpleCPW(10.0μm, 6.0μm)
    p = Path(Point(0.0μm, 0.0μm); α0=0.0°)
    straight!(p, 200.0μm, style)
    turn!(p, -90°, r)
    straight!(p, 300.0μm)
    turn!(p, 180°, r)
    straight!(p, 300.0μm)

    cs = CoordinateSystem("reson", nm)
    place!(cs, p)

    k = DLE.EGADSKernel(enable_tape=true)
    sm = SolidModel("reson3d", k)
    render!(sm, cs; map_meta=m -> string(DeviceLayout.layer(m)))

    DLE.start_recording!(k.geometry_tape)
    for (_, pg) in DeviceLayout.SolidModels.dimgroupdict(sm, 2)
        DeviceLayout.SolidModels.extrude_z!(pg, 10.0μm)
    end
    DLE.stop_recording!(k.geometry_tape)

    solids = EGADS.EGObject[]
    for b in DLE._bodies(k)
        info = EGADS.get_info(b)
        if info.oclass == EGADS.OClass.BODY && info.mtype == EGADS.MType.SOLIDBODY
            push!(solids, b)
        end
    end

    # Build substrate + sequential subtraction
    ctx = DLE._ctx(k)
    bb_min = [Inf, Inf, Inf]; bb_max = [-Inf, -Inf, -Inf]
    for s in solids
        bb = EGADS.get_bounding_box(s)
        bb_min .= min.(bb_min, [bb.xmin, bb.ymin, bb.zmin])
        bb_max .= max.(bb_max, [bb.xmax, bb.ymax, bb.zmax])
    end
    margin = 50.0
    sub = EGADS.make_box(ctx, bb_min[1]-margin, bb_min[2]-margin, 0.0,
                         bb_max[1]-bb_min[1]+2margin, bb_max[2]-bb_min[2]+2margin, 100.0)
    result = sub
    for gap in solids
        rv = EGADS.solid_boolean(result, gap, EGADS.BoolOp.SUBTRACTION)
        isempty(rv) && continue
        result = rv[1]
    end

    faces = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    edges = EGADS.get_body_topos(result, EGADS.OClass.EDGE)
    println("Result: $(length(faces)) faces, $(length(edges)) edges")

    # Count edge types
    n_circle = count(edges) do edge
        topo = EGADS.get_topology(edge)
        geom = topo.geom
        !isnothing(geom) && Int(EGADS.get_info(geom).mtype) == 2
    end
    println("  CIRCLE edges: $n_circle / $(length(edges))")

    # Tessellate with fine resolution to show curves
    tess = EGADS.tessellate(result; max_edge=10.0, sag=0.5, angle=8.0)
    npts = EGADS.tess_status(tess).npts
    println("Tessellation: $npts surface vertices")

    # Collect P1 mesh
    pts = zeros(Float64, 3, npts)
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end

    all_tris = Vector{NTuple{3,Int}}()
    tri_face_ids = Int[]
    for fi in eachindex(faces)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for t in 1:size(fm.tris, 1)
            g1 = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,1])
            g2 = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,2])
            g3 = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,3])
            push!(all_tris, (g1, g2, g3))
            push!(tri_face_ids, fi)
        end
    end
    println("Triangles: $(length(all_tris))")

    # Build P2 mid-edge nodes by evaluating EGADS faces at edge midpoints
    # For each triangle edge, compute the midpoint in physical space, then
    # project onto the EGADS face to get the exact curved position.
    edge_to_midnode = Dict{Tuple{Int,Int}, Int}()
    midnodes = Vector{Vector{Float64}}()
    next_id = npts + 1

    for (ti, (g1, g2, g3)) in enumerate(all_tris)
        fi = tri_face_ids[ti]
        face = faces[fi]
        for (a, b) in [(g1,g2), (g2,g3), (g3,g1)]
            ekey = a < b ? (a, b) : (b, a)
            haskey(edge_to_midnode, ekey) && continue

            # Physical midpoint
            mid_phys = (pts[:, a] .+ pts[:, b]) ./ 2

            # Project onto EGADS face → get the curved position
            inv_res = EGADS.inv_evaluate(face, mid_phys)
            mid_curved = inv_res.results[1:3]

            push!(midnodes, mid_curved)
            edge_to_midnode[ekey] = next_id
            next_id += 1
        end
    end
    println("P2 mid-edge nodes: $(length(midnodes))")
    println("Total nodes (P1+P2): $(npts + length(midnodes))")

    # Show curvature deviation: max distance between chord midpoint and curved midpoint
    max_dev = 0.0
    for (ti, (g1, g2, g3)) in enumerate(all_tris)
        for (a, b) in [(g1,g2), (g2,g3), (g3,g1)]
            ekey = a < b ? (a, b) : (b, a)
            mid_id = edge_to_midnode[ekey]
            mid_curved = midnodes[mid_id - npts]
            mid_chord = (pts[:, a] .+ pts[:, b]) ./ 2
            dev = sqrt(sum((mid_curved .- mid_chord).^2))
            max_dev = max(max_dev, dev)
        end
    end
    @printf("Max curvature deviation (P1→P2): %.4f μm\n", max_dev)

    # Write Gmsh MSH v4.1 with P2 triangles (6-node, type 9)
    msh_file = joinpath(@__DIR__, "resonator_curved.msh")
    total_nodes = npts + length(midnodes)
    open(msh_file, "w") do io
        println(io, "\$MeshFormat")
        println(io, "4.1 0 8")
        println(io, "\$EndMeshFormat")

        # Nodes
        println(io, "\$Nodes")
        println(io, "1 $total_nodes 1 $total_nodes")  # 1 entity block
        println(io, "2 1 0 $total_nodes")  # dim=2, tag=1, parametric=0, numNodes
        for i in 1:total_nodes
            println(io, i)
        end
        for i in 1:npts
            @printf(io, "%.10f %.10f %.10f\n", pts[1,i], pts[2,i], pts[3,i])
        end
        for mn in midnodes
            @printf(io, "%.10f %.10f %.10f\n", mn[1], mn[2], mn[3])
        end
        println(io, "\$EndNodes")

        # Elements: P2 triangles (type 9 = 6-node triangle)
        println(io, "\$Elements")
        println(io, "1 $(length(all_tris)) 1 $(length(all_tris))")
        println(io, "2 1 9 $(length(all_tris))")  # dim=2, tag=1, type=9
        for (ti, (g1, g2, g3)) in enumerate(all_tris)
            e12 = edge_to_midnode[g1 < g2 ? (g1,g2) : (g2,g1)]
            e23 = edge_to_midnode[g2 < g3 ? (g2,g3) : (g3,g2)]
            e31 = edge_to_midnode[g3 < g1 ? (g3,g1) : (g1,g3)]
            println(io, "$ti $g1 $g2 $g3 $e12 $e23 $e31")
        end
        println(io, "\$EndElements")
    end
    println("\nWrote: $msh_file")
    println("Open with: gmsh $msh_file")
end

main()
