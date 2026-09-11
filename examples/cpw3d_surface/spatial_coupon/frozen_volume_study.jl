# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Paired volume-meshing experiment: one fixed surface mesh, several independent
# 3D grading profiles. No prism extrusion, remeshing of the interfaces, or slot
# integration changes are permitted between candidates.
function frozen_surface_fingerprint()
    tags,xyz,_=gmsh.model.mesh.getNodes()
    coordinates=Dict(tag=>ntuple(d->xyz[3i-3+d],3) for (i,tag) in enumerate(tags))
    records=Tuple{Int,NTuple{3,NTuple{3,Float64}}}[]
    used=Set{UInt64}()
    attributes=Int[]
    for (_,attribute) in gmsh.model.getPhysicalGroups(2)
        push!(attributes,Int(attribute))
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2,attribute)
            types,_,connectivity=gmsh.model.mesh.getElements(2,entity)
            for (type,conn) in zip(types,connectivity)
                _,_,order,nnode,_,primary=gmsh.model.mesh.getElementProperties(type)
                order==1 && nnode==3 && primary==3 || error("Frozen surface study requires linear triangles")
                for i in 1:3:length(conn)
                    nodes=conn[i:i+2]
                    union!(used,nodes)
                    points=Tuple(sort([coordinates[node] for node in nodes]))
                    push!(records,(Int(attribute),points))
                end
            end
        end
    end
    length(unique(coordinates[node] for node in used))==length(used) ||
        error("Coincident distinct surface nodes are not supported by this study")
    sort!(records)
    length(unique(records))==length(records) || error("Duplicate physical surface triangles")
    buffer=IOBuffer()
    write(buffer,"frozen-surface-v1")
    for (attribute,points) in records
        write(buffer,Int64(attribute))
        for p in points,coordinate in p;write(buffer,coordinate);end
    end
    return Dict("SHA256"=>bytes2hex(sha256(take!(buffer))),"Triangles"=>length(records),
                "Nodes"=>length(used),"Attributes"=>sort!(attributes))
end

# Reconstruct the same material shells from a saved surface without rerunning
# expensive surface meshing. The surface file remains a Gmsh 2.2 checkpoint.
function restore_frozen_surface(resume)
    Set(keys(resume))==Set(["Mesh","SHA256"]) || error("Invalid ResumeSurface settings")
    path=String(resume["Mesh"])
    bytes2hex(sha256(read(path)))==resume["SHA256"] || error("Saved surface hash mismatch")
    shells=[(Int(volume),[Int(tag) for (_,tag) in gmsh.model.getBoundary([(3,volume)],false,true,false)])
            for (_,volume) in gmsh.model.getEntities(3)]
    materials=[(Int(attribute),Int.(gmsh.model.getEntitiesForPhysicalGroup(3,attribute)))
               for (_,attribute) in gmsh.model.getPhysicalGroups(3)]
    gmsh.clear()
    gmsh.open(path)
    isempty(gmsh.model.getEntities(3)) || error("ResumeSurface must contain only the saved boundary")
    surfaces=Set(Int(tag) for (_,tag) in gmsh.model.getEntities(2))
    Set(abs(tag) for (_,shell) in shells for tag in shell)==surfaces ||
        error("Saved surface entities differ from the regenerated CAD shells")
    for (volume,shell) in shells
        loop=gmsh.model.geo.addSurfaceLoop(shell)
        gmsh.model.geo.addVolume([loop],volume)
    end
    gmsh.model.geo.synchronize()
    for (attribute,volumes) in materials
        gmsh.model.addPhysicalGroup(3,volumes,attribute)
    end
    callback=@cfunction(exact_grading_callback,Cdouble,(Cint,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cvoid}))
    push!(SIZE_CALLBACK_ROOTS,callback)
    ierr=Ref{Cint}()
    ccall((:gmshModelMeshSetSizeCallback,gmsh.lib),Cvoid,
          (Ptr{Cvoid},Ptr{Cvoid},Ptr{Cint}),callback,C_NULL,ierr)
    ierr[]==0 || error(gmsh.logger.getLastError())
    return nothing
end

function run_frozen_volume_study(settings,output,lower,upper,fine)
    Set(keys(settings))<=Set(["Variants","MaxElements","MaxNodes","OptimizeVolume","ResumeSurface"]) ||
        error("Unknown volume-study setting")
    variants=get(settings,"Variants",nothing)
    variants isa Vector && !isempty(variants) || error("Volume study requires Variants")
    maximum_elements=get(settings,"MaxElements",5_000_000)
    maximum_nodes=get(settings,"MaxNodes",1_500_000)
    maximum_elements>0 && maximum_nodes>0 || error("Invalid volume-study mesh budget")
    optimize=get(settings,"OptimizeVolume",true)
    names=String[]
    profiles=NamedTuple[]
    for variant in variants
        required=Set(["Name","NearGrowth","FarGrowth","TransitionDistance","MaximumSize"])
        required<=Set(keys(variant)) && Set(keys(variant))<=union(required,Set(["MinimumSize"])) ||
            error("Invalid volume variant keys")
        name=String(variant["Name"])
        occursin(r"^[A-Za-z0-9_-]+$",name) || error("Unsafe volume variant name")
        name in names && error("Duplicate volume variant")
        profile=(near_growth=Float64(variant["NearGrowth"]),
                 far_growth=Float64(variant["FarGrowth"]),
                 transition_distance=Float64(variant["TransitionDistance"]),
                 maximum_size=Float64(variant["MaximumSize"]),
                 minimum_size=Float64(get(variant,"MinimumSize",fine)))
        all(isfinite,values(profile)) && profile.near_growth>0 && profile.far_growth>0 &&
            profile.transition_distance>=0 && 0<profile.minimum_size<=profile.maximum_size || error("Invalid volume grading profile")
        filename=splitext(output)[1]*"-"*name*".msh"
        isfile(filename) && error("Refuse to overwrite volume candidate")
        push!(names,name);push!(profiles,profile)
    end
    gmsh.option.setNumber("Mesh.Renumber",0)
    started=time()
    if haskey(settings,"ResumeSurface")
        restore_frozen_surface(settings["ResumeSurface"])
    else
        gmsh.model.mesh.generate(2)
    end
    surface_seconds=time()-started
    surface=frozen_surface_fingerprint()
    gmsh.write(splitext(output)[1]*"-surface.msh")
    report=Dict{String,Any}("Scope"=>"Fixed-surface volume mesh comparison; no accuracy qualification",
        "Boundary"=>surface,"Lower"=>collect(lower),"Upper"=>collect(upper),
        "SurfaceSeconds"=>surface_seconds,"GeometryOrder"=>1,
        "ResumedSurface"=>haskey(settings,"ResumeSurface"),
        "SurfaceSizing"=>Dict("MinimumSize"=>fine,"MaximumSize"=>GRADING_FAR[],
                              "Growth"=>GRADING_GROWTH[],"TangentialEdgeSize"=>GRADING_TANGENT[]),
        "GmshVersion"=>gmsh.option.getString("General.Version"),"Variants"=>Dict{String,Any}[])
    report_path=output*".volume-study.toml"
    function save_report()
        open(report_path,"w") do stream;TOML.print(stream,report;sorted=true);end
    end
    save_report()
    volumes=gmsh.model.getEntities(3)
    try
        for (name,profile) in zip(names,profiles)
            gmsh.model.mesh.clear(volumes)
            frozen_surface_fingerprint()==surface || error("Clearing the volume altered the boundary mesh")
            GRADING_VOLUME_PROFILE[]=profile
            fill!(GRADING_QUERY_COUNTS,0)
            gmsh.option.setNumber("Mesh.MeshSizeMin",min(fine,profile.minimum_size,
                GRADING_TRACE_SIZE[]>0 ? GRADING_TRACE_SIZE[] : fine))
            gmsh.option.setNumber("Mesh.MeshSizeMax",max(GRADING_FAR[],profile.maximum_size))
            started=time()
            println("Volume variant $name: $profile");flush(stdout)
            gmsh.model.mesh.generate(3)
            generated=time()
            types,element_tags,_=gmsh.model.mesh.getElements(3)
            types==[4] || error("Volume study produced non-tetrahedral elements")
            elements=sum(length,element_tags)
            nodes=length(gmsh.model.mesh.getNodes()[1])
            println("Volume counts before optimization: elements=$elements nodes=$nodes");flush(stdout)
            elements<=maximum_elements && nodes<=maximum_nodes ||
                error("Volume study exceeded its mesh budget: elements=$elements/$maximum_elements nodes=$nodes/$maximum_nodes")
            frozen_surface_fingerprint()==surface || error("Volume mesher changed an input boundary triangle/node")
            optimize && gmsh.model.mesh.optimize("Netgen")
            frozen_surface_fingerprint()==surface || error("Volume optimization changed the fixed boundary")
            elements=sum(length,gmsh.model.mesh.getElements(3)[2])
            nodes=length(gmsh.model.mesh.getNodes()[1])
            elements<=maximum_elements && nodes<=maximum_nodes || error("Optimized volume exceeds mesh budget")
            quality=gmsh.model.mesh.getElementQualities(reduce(vcat,gmsh.model.mesh.getElements(3)[2]),"minSICN")
            minimum(quality)>1e-10 || error("Nonpositive/near-singular volume element")
            GRADING_QUERY_COUNTS[5]>0 || error("The volume mesher did not query the 3D grading profile")
            filename=splitext(output)[1]*"-"*name*".msh"
            gmsh.write(filename)
            entry=Dict{String,Any}("Name"=>name,"Mesh"=>filename,"MeshSHA256"=>bytes2hex(sha256(read(filename))),
                "Elements"=>elements,"Nodes"=>nodes,"VolumeGenerationSeconds"=>generated-started,
                "VolumeTotalSeconds"=>time()-started,"MinimumSICN"=>minimum(quality),
                "BoundarySHA256"=>surface["SHA256"],"SizeQueryCounts"=>copy(GRADING_QUERY_COUNTS),
                "NearGrowth"=>profile.near_growth,"FarGrowth"=>profile.far_growth,
                "TransitionDistance"=>profile.transition_distance,"MaximumSize"=>profile.maximum_size,
                "MinimumVolumeSize"=>profile.minimum_size)
            push!(report["Variants"],entry);save_report()
            println("Completed $name: $elements tets, $nodes nodes; boundary unchanged");flush(stdout)
        end
    finally
        GRADING_VOLUME_PROFILE[]=nothing
    end
    return report
end
