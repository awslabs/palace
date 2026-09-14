# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Read-only extraction. Preserve Gmsh element/node IDs, quadratic geometry and physical
# matching attribute. Ground segments are proved by original physical face incidence,
# not guessed from the open boundary of the extracted surface.
import Gmsh: gmsh
using SHA, TOML, Printf

function main(source, output)
    ispath(output * ".msh") && error("Refusing to overwrite extraction")
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.NumThreads", 1)
        gmsh.open(source)
        tags, xyz, _ = gmsh.model.mesh.getNodes()
        coordinates = Dict(Int(t) => xyz[3i-2:3i] for (i,t) in enumerate(tags))
        # Supported complete linear/quadratic triangle and quad types only.
        properties = Dict(2=>(3,3), 3=>(4,4), 9=>(3,6), 10=>(4,9))
        faces = Tuple{Int,Int,Vector{Int}}[]
        matching_edges = Dict{Tuple{Int,Int},Vector{Vector{Int}}}()
        ground_edges = Dict{Tuple{Int,Int},Vector{Tuple{Int,Vector{Int}}}}()
        matching_vertices = Set{Int}(); ground_vertices = Set{Int}()
        counts = Dict{String,Int}()
        for attr in (1,5001,6001)
            entities = gmsh.model.getEntitiesForPhysicalGroup(2,attr)
            for entity in entities
                types, ids, connections = gmsh.model.mesh.getElements(2,entity)
                for (type, element_ids, conn) in zip(types,ids,connections)
                    haskey(properties,Int(type)) || error("Unsupported surface type $type")
                    corners, nodes = properties[Int(type)]
                    counts["attribute_$(attr)_type_$(type)"] = get(counts,"attribute_$(attr)_type_$(type)",0)+length(element_ids)
                    for (i,id) in enumerate(element_ids)
                        c = Int.(conn[nodes*(i-1)+1:nodes*i])
                        if attr == 1
                            push!(faces,(Int(id),Int(type),c)); union!(matching_vertices,c[1:corners])
                        else
                            union!(ground_vertices,c[1:corners])
                        end
                        for j in 1:corners
                            a,b=c[j],c[mod1(j+1,corners)]
                            edge = nodes == corners ? [a,b] : [a,b,c[corners+j]]
                            key = minmax(a,b)
                            if attr == 1
                                push!(get!(matching_edges,key,Vector{Int}[]),edge)
                            else
                                push!(get!(ground_edges,key,Tuple{Int,Vector{Int}}[]),(attr,edge))
                            end
                        end
                    end
                end
            end
        end
        segments = Tuple{Int,Vector{Int}}[]
        open_edges = 0
        covered = Set{Int}()
        for (key, edges) in sort!(collect(matching_edges);by=first)
            length(edges) in (1,2) || error("Nonmanifold matching edge")
            isopen = length(edges)==1
            open_edges += isopen
            isground = haskey(ground_edges,key)
            isopen == isground || error("Open-edge / physical-ground mismatch at $key")
            if isground
                ground = ground_edges[key]
                all(sort(e)==sort(edges[1]) for (_,e) in ground) || error("Ground geometry edge mismatch")
                attr = minimum(first.(ground))
                push!(segments,(attr,edges[1])); union!(covered,key)
            end
        end
        shared = intersect(matching_vertices,ground_vertices)
        shared == covered || error("Isolated ground vertex contacts require explicit handling")
        used = sort!(collect(union((Set(c) for (_,_,c) in faces)...)))
        open(output * ".msh","w") do io
            println(io,"\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n",length(used))
            for t in used
                @printf(io,"%d %.17g %.17g %.17g\n",t,coordinates[t]...)
            end
            println(io,"\$EndNodes\n\$Elements\n",length(faces)+length(segments))
            for (id,type,c) in faces
                println(io,join([id,type,2,1,1,c...]," "))
            end
            first_id = maximum(first.(faces))
            for (i,(attr,c)) in enumerate(segments)
                type = length(c)==2 ? 1 : 8
                println(io,join([first_id+i,type,2,attr,attr,c...]," "))
            end
            println(io,"\$EndElements")
        end
        report = Dict("SourceMesh"=>abspath(source),"SourceSHA256"=>bytes2hex(open(sha256,source)),
            "ExtractSHA256"=>bytes2hex(open(sha256,output*".msh")),"MatchingFaces"=>length(faces),
            "MatchingGeometryNodes"=>length(used),"GroundContactEdges"=>length(segments),
            "OpenEdges"=>open_edges,"SharedGroundCornerVertices"=>length(shared),
            "OpenEdgesExactlyPhysicalGroundContacts"=>true,"OriginalSurfaceCounts"=>counts)
        open(output * ".toml","w") do io; TOML.print(io,report); end
        println(report)
    finally
        gmsh.finalize()
    end
end
length(ARGS)==2 || error("usage: source.msh output-prefix")
main(ARGS...)
