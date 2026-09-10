# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Isolated experiment harness, never writes the reference mesh.
# Args: original_mesher.jl signature_dir exported_plan_prefix output_dir kind variant
include(abspath(ARGS[1]))
include(joinpath(@__DIR__, "repair_plan_triangles.jl"))
include(joinpath(@__DIR__, "relocate_plan_rows.jl"))
include(joinpath(@__DIR__, "align_corner_rows.jl"))

function main()
length(ARGS) in (6, 7) || error("Expected mesher, signature dir, plan prefix, output dir, kind, variant")
plan_only = length(ARGS) == 7 && ARGS[7] == "--plan-only"
length(ARGS) == 6 || plan_only || error("Unknown extra argument")
signature_dir, prefix, output = abspath.(ARGS[2:4])
kind, variant = ARGS[5:6]
kind in ("thin","fabricated") || error("invalid kind")
variant in ("baseline","repair","repair-coarsen-z","repair-relocate","repair-relocate-coarsen-z","aligned","aligned-coarsen-z","aligned-smooth","aligned-smooth-coarsen-z","aligned-nc-seed") || error("invalid variant")
fabricated = kind == "fabricated"
radius, thickness, overetch, normal, tangent, far, core, growth = 2.0,0.1,0.05,0.002,0.1,0.5,0.8,1.4
edges = read_edges(joinpath(signature_dir,"mesh-signature.csv"))
facets = read_mask(joinpath(signature_dir,"plan-view-mask.csv"))
loops = read_boundary(joinpath(signature_dir,"plan-view-boundary.csv"))
lower,upper = coupon_bounds(edges,radius,thickness,overetch)
nd,_ = readdlm(prefix*"-nodes.csv",',',header=true)
td,_ = readdlm(prefix*"-triangles.csv",',',header=true)
coords = Dict{UInt64,NTuple{2,Float64}}(UInt64(nd[i,1]) => (nd[i,2],nd[i,3]) for i in axes(nd,1))
triangles = [Tuple(UInt64.(td[i,2:4])) for i in axes(td,1)]
plan = (node_tags=sort!(collect(keys(coords))),coordinates=coords,triangles=triangles,
        active_segments=active_plan_segments(edges,radius),
        fabrication_primitives=classified_fabrication_primitives(loops,1e-9radius),
        normal_distances=geometric_distances(core,normal,growth),
        maximum_tangent_spacing=tangent,maximum_far_spacing=far,
        embedded_curve_count=0,plan_surface_count=0,unmatched_internal_edge_count=0)
mkpath(output)
metrics = nothing
if variant != "baseline"
    plan,metrics = repair_plan_triangles(plan,facets,edges,radius)
end
if startswith(variant,"aligned")
    plan,alignment = align_corner_rows(plan,facets,edges,radius)
    metrics=(initial=metrics,alignment=alignment)
end
if startswith(variant,"aligned-smooth") || variant=="aligned-nc-seed"
    plan,smoothed = relocate_plan_rows(plan,facets,edges,radius;
                                      edge_growth=1.5,infer_interior_axis=true)
    plan,final_flips = repair_plan_triangles(plan,facets,edges,radius;max_passes=40)
    metrics=(alignment=metrics,smoothed=smoothed,final_flips=final_flips)
end
if startswith(variant,"repair-relocate")
    plan,relocation = relocate_plan_rows(plan,facets,edges,radius)
    metrics=(flips=metrics,relocation=relocation)
end
z = edge_cluster_z_coordinates(edges[1].point[3],lower[3],upper[3],fabricated,
                               thickness,overetch,normal,growth)
original_z = copy(z)
if variant=="aligned-nc-seed"
    planes = fabricated ? [edges[1].point[3]-overetch,edges[1].point[3],edges[1].point[3]+thickness] : [edges[1].point[3]]
    z = [value for value in z if minimum(abs(value-plane) for plane in planes)>0.2 ||
                                minimum(abs(value-plane) for plane in planes)<1.0e-10]
    println("NC seed vertical intervals: $(length(original_z)-1) -> $(length(z)-1)")
end
if endswith(variant,"coarsen-z")
    # Conservative first coarsening experiment: merge only far vertical intervals.
    # This does NOT implement different plan triangulations at different heights.
    planes = fabricated ? [edges[1].point[3]-overetch,edges[1].point[3],edges[1].point[3]+thickness] : [edges[1].point[3]]
    keep = trues(length(z))
    i = 2
    while i < length(z)
        if all(minimum(abs(value-plane) for plane in planes) > 0.2 for value in z[i-1:i+1]) &&
           z[i+1]-z[i-1] <= 0.65
            keep[i] = false
            i += 2
        else
            i += 1
        end
    end
    z = z[keep]
    println("Far-z interval merging: $(length(original_z)-1) -> $(length(z)-1) intervals")
end
open(joinpath(output,"experiment.txt"),"w") do f
    println(f,"variant=$variant kind=$kind")
    println(f,"repair=$metrics")
    println(f,"original_z=$original_z")
    println(f,"z=$z")
end
open(joinpath(output,"plan-nodes.csv"),"w") do f
    println(f,"node,x,y")
    for id in sort!(collect(keys(plan.coordinates)))
        x,y=plan.coordinates[id];println(f,"$id,$x,$y")
    end
end
open(joinpath(output,"plan-triangles.csv"),"w") do f
    println(f,"triangle,a,b,c")
    for (i,t) in enumerate(plan.triangles)
        println(f,"$i,$(t[1]),$(t[2]),$(t[3])")
    end
end
plan_only && return
filename = joinpath(output,"$kind.msh")
isfile(filename) && error("Refuse to overwrite candidate mesh")
gmsh.initialize()
try
    gmsh.option.setNumber("General.Verbosity",2)
    gmsh.model.add("coupon_mesh_experiment")
    evidence = add_discrete_edge_cluster_mesh!(plan,edges,facets,fabricated,z,lower,upper,
                radius,thickness,overetch,30_000_000,20_000_000,2)
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    gmsh.option.setNumber("Mesh.Binary",1)
    gmsh.write(filename)
    gmsh.clear();gmsh.open(filename)
    all_tags = reduce(vcat,gmsh.model.mesh.getElements(3)[2];init=UInt64[])
    minimum(gmsh.model.mesh.getElementQualities(all_tags,"minSJ")) > 0 || error("Invalid serialized Jacobian")
    prism_face_counts() == evidence.face_counts || error("Serialized topology differs")
    write_edge_cluster_metadata(filename*".metadata.json",evidence,plan,z,normal,tangent,far,core,growth,2,fabricated,length(edges))
finally
    gmsh.finalize()
end
println(filename)
end

main()
