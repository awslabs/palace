# Visualize the resonator mesh geometry.
# Outputs an OBJ surface mesh + STEP file for inspection.

using DeviceLayout, DeviceLayout.PreferredUnits, Unitful, EGADS, Printf

DLE = Base.get_extension(DeviceLayout, :DeviceLayoutEGADSExt)

function main()
    # Build the native-curve resonator
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
    println("$(length(solids)) solid bodies")

    # Inspect edge types
    for (si, solid) in enumerate(solids)
        edges = EGADS.get_body_topos(solid, EGADS.OClass.EDGE)
        counts = Dict{Int, Int}()
        for edge in edges
            topo = EGADS.get_topology(edge)
            geom = topo.geom
            if !isnothing(geom)
                gi = EGADS.get_info(geom)
                mt = Int(gi.mtype)
                counts[mt] = get(counts, mt, 0) + 1
            end
        end
        # EGADS MType: 1=LINE, 2=CIRCLE, 3=ELLIPSE, 4=PARABOLA, 5=HYPERBOLA, 6=TRIMMED, 7=BEZIER, 8=BSPLINE
        type_names = Dict(1=>"LINE", 2=>"CIRCLE", 3=>"ELLIPSE", 8=>"BSPLINE")
        desc = join(["$(get(type_names, k, "TYPE_$k")): $v" for (k, v) in sort(collect(counts))], ", ")
        println("  Solid $si: $(length(edges)) edges — $desc")
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

    result_faces = EGADS.get_body_topos(result, EGADS.OClass.FACE)
    result_edges = EGADS.get_body_topos(result, EGADS.OClass.EDGE)
    println("\nResult: $(length(result_faces)) faces, $(length(result_edges)) edges")

    # Edge types on final body
    counts = Dict{Int, Int}()
    for edge in result_edges
        topo = EGADS.get_topology(edge)
        geom = topo.geom
        if !isnothing(geom)
            gi = EGADS.get_info(geom)
            counts[Int(gi.mtype)] = get(counts, Int(gi.mtype), 0) + 1
        end
    end
    type_names = Dict(1=>"LINE", 2=>"CIRCLE", 3=>"ELLIPSE", 8=>"BSPLINE")
    println("Edge types: ", join(["$(get(type_names, k, "TYPE_$k")): $v" for (k, v) in sort(collect(counts))], ", "))

    # Face surface types
    face_types = Dict{Int, Int}()
    for face in result_faces
        topo = EGADS.get_topology(face)
        geom = topo.geom
        if !isnothing(geom)
            gi = EGADS.get_info(geom)
            face_types[Int(gi.mtype)] = get(face_types, Int(gi.mtype), 0) + 1
        end
    end
    # Surface MType: 1=PLANE, 2=SPHERICAL, 3=CYLINDRICAL, 4=REVOLUTION, 5=TOROIDAL, 6=TRIMMED, 7=BEZIER, 8=BSPLINE, 9=OFFSET, 10=CONICAL, 11=EXTRUSION
    surf_names = Dict(1=>"PLANE", 2=>"SPHERICAL", 3=>"CYLINDRICAL", 5=>"TOROIDAL", 11=>"EXTRUSION")
    println("Surface types: ", join(["$(get(surf_names, k, "SURF_$k")): $v" for (k, v) in sort(collect(face_types))], ", "))

    # Tessellate (fine enough to see curves)
    tess = EGADS.tessellate(result; max_edge=15.0, sag=1.0, angle=10.0)
    npts = EGADS.tess_status(tess).npts
    println("\nTessellation: $npts surface vertices")

    pts = zeros(Float64, 3, npts)
    for fi in eachindex(result_faces)
        fm = EGADS.get_face_mesh(tess, result, fi)
        for j in 1:size(fm.xyz, 1)
            ig = EGADS.eg_local_to_global(tess.ptr, fi, j)
            pts[:, ig] .= fm.xyz[j, :]
        end
    end

    # Write OBJ
    obj_file = joinpath(@__DIR__, "resonator_mesh.obj")
    open(obj_file, "w") do io
        println(io, "# Resonator cavity mesh (substrate - CPW gaps)")
        println(io, "# Units: μm (STP_UNIT)")
        for i in 1:npts
            @printf(io, "v %.6f %.6f %.6f\n", pts[1,i], pts[2,i], pts[3,i])
        end
        for fi in eachindex(result_faces)
            fm = EGADS.get_face_mesh(tess, result, fi)
            for t in 1:size(fm.tris, 1)
                g1 = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,1])
                g2 = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,2])
                g3 = EGADS.eg_local_to_global(tess.ptr, fi, fm.tris[t,3])
                println(io, "f $g1 $g2 $g3")
            end
        end
    end
    println("Wrote: $obj_file")

    # Save as STEP
    step_file = joinpath(@__DIR__, "resonator_result.step")
    EGADS.save_model(result, step_file)
    println("Wrote: $step_file")

    # Save as EGADS
    egads_file = joinpath(@__DIR__, "resonator_result.egads")
    EGADS.save_model(result, egads_file)
    println("Wrote: $egads_file")
end

main()
