# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Experimental tetrahedral bulk replacement with a fixed triangular boundary.
# Uses a built-in volume bounded by an already-meshed discrete surface. There is
# no OCC remeshing, geometric simplification, or modification of the input shell.
using Gmsh: gmsh
using LinearAlgebra
using DelimitedFiles

triangle_key(t) = Tuple(sort(collect(t)))

function check_closed_shell(points, triangles)
    isempty(triangles) && error("Empty shell")
    edges = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    seen = Set{NTuple{3,Int}}()
    volume = 0.0
    for (i,t) in enumerate(triangles)
        length(unique(t)) == 3 || error("Degenerate shell connectivity")
        all(1 <= v <= length(points) for v in t) || error("Invalid shell node")
        key = triangle_key(t)
        key in seen && error("Duplicate shell triangle")
        push!(seen,key)
        a,b,c = (collect(points[v]) for v in t)
        norm(cross(b-a,c-a)) > 0 || error("Zero-area shell triangle")
        volume += dot(a,cross(b,c))/6
        for (u,v) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
            push!(get!(edges,minmax(u,v),Tuple{Int,Int}[]),(i,u<v ? 1 : -1))
        end
    end
    all(length(v)==2 && v[1][2]+v[2][2]==0 for v in values(edges)) ||
        error("Shell must be closed and consistently oriented")
    neighbors = [Int[] for _ in triangles]
    for v in values(edges)
        a,b=v[1][1],v[2][1];push!(neighbors[a],b);push!(neighbors[b],a)
    end
    visited=Set([1]);pending=[1]
    while !isempty(pending)
        for j in neighbors[pop!(pending)]
            if !(j in visited);push!(visited,j);push!(pending,j);end
        end
    end
    length(visited)==length(triangles) || error("Expected one connected shell")
    volume > 0 || error("Shell must have positive oriented volume")
    return volume
end

# Preserve the complete input XY triangulation on both end caps. Side quads are
# split consistently; this changes their FE trace space and is NOT a promise of
# equivalence to a reference tensor-product boundary discretization.
function extruded_plan_shell(xy, triangles, heights)
    length(heights)>=2 && issorted(heights) && all(diff(heights).>0) ||
        error("Heights must be strictly increasing")
    n=length(xy)
    points=[(p[1],p[2],z) for z in heights for p in xy]
    face=NTuple{3,Int}[]
    edges=Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    for t in triangles
        a,b,c=(xy[v] for v in t)
        area=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
        area>0 || error("Plan triangles must be counterclockwise")
        push!(face,(t[1],t[3],t[2]))
        offset=n*(length(heights)-1)
        push!(face,Tuple(v+offset for v in t))
        for (u,v) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
            push!(get!(edges,minmax(u,v),Tuple{Int,Int}[]),(u,v))
        end
    end
    for adj in values(edges)
        length(adj)<=2 || error("Nonmanifold plan edge")
        if length(adj)==1
            u,v=adj[1]
            for k in 0:length(heights)-2
                a,b,c,d=u+k*n,v+k*n,v+(k+1)*n,u+(k+1)*n
                push!(face,(a,b,c),(a,c,d))
            end
        else
            adj[1]==reverse(adj[2]) || error("Inconsistent plan orientation")
        end
    end
    used=sort!(unique!(collect(Iterators.flatten(face))))
    index=Dict(old=>new for (new,old) in enumerate(used))
    points=points[used]
    face=[Tuple(index[v] for v in t) for t in face]
    check_closed_shell(points,face)
    return points,face
end

function cap_layer_seeds(points,triangles)
    zmin,zmax=extrema(p[3] for p in points)
    seeds=NTuple{3,Float64}[]
    for t in triangles
        a,b,c=(collect(points[v]) for v in t)
        all(p[3]==zmin for p in (a,b,c)) || all(p[3]==zmax for p in (a,b,c)) || continue
        center=(a+b+c)/3
        longest=max(norm(b-a),norm(c-a),norm(c-b))
        height=min(0.5longest,0.2(zmax-zmin))
        push!(seeds,(center[1],center[2],center[3]+(center[3]<(zmin+zmax)/2 ? height : -height)))
    end
    return seeds
end

function mesh_fixed_shell(points,triangles,output; maximum_size=0.5,algorithm=10,
                          optimize_threshold=0.3, seeds=NTuple{3,Float64}[])
    isfile(output) && error("Refuse to overwrite existing mesh")
    isfinite(maximum_size) && maximum_size>0 || error("Invalid size")
    isfinite(optimize_threshold) && 0<=optimize_threshold<1 || error("Invalid quality target")
    expected_volume=check_closed_shell(points,triangles)
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.option.setNumber("Mesh.MaxNumThreads3D",1)
        gmsh.option.setNumber("Mesh.Algorithm3D",algorithm)
        gmsh.option.setNumber("Mesh.OptimizeThreshold",optimize_threshold)
        gmsh.option.setNumber("Mesh.Renumber",0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints",0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",0)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary",0)
        gmsh.option.setNumber("Mesh.MeshSizeMax",maximum_size)
        gmsh.model.add("fixed_shell_bulk")
        surface=gmsh.model.addDiscreteEntity(2)
        tags=collect(1:length(points))
        gmsh.model.mesh.addNodes(2,surface,tags,collect(Iterators.flatten(points)))
        gmsh.model.mesh.addElementsByType(surface,2,[],collect(Iterators.flatten(triangles)))
        shell=gmsh.model.geo.addSurfaceLoop([surface])
        volume=gmsh.model.geo.addVolume([shell])
        seed_tags=[gmsh.model.geo.addPoint(p...) for p in seeds]
        gmsh.model.geo.synchronize()
        isempty(seed_tags) || gmsh.model.mesh.embed(0,seed_tags,3,volume)
        gmsh.model.addPhysicalGroup(2,[surface],1)
        gmsh.model.addPhysicalGroup(3,[volume],1)
        gmsh.model.mesh.generate(3)
        types,elements,connectivity=gmsh.model.mesh.getElements(3)
        types==[4] || error("Expected only linear tetrahedra")
        # Fixed node identities and triangle connectivity are essential for merging.
        actual=gmsh.model.mesh.getElementsByType(2,surface)[2]
        actual_faces=Set(triangle_key(Tuple(actual[i:i+2])) for i in 1:3:length(actual))
        expected_faces=Set(triangle_key(t) for t in triangles)
        actual_faces==expected_faces || error("Boundary was remeshed: expected=$(length(expected_faces)) actual=$(length(actual_faces)) missing=$(collect(setdiff(expected_faces,actual_faces))[1:min(end,3)]) added=$(collect(setdiff(actual_faces,expected_faces))[1:min(end,3)])")
        for (id,p) in enumerate(points)
            actual_p=gmsh.model.mesh.getNode(id)[1]
            norm(actual_p-collect(p))<=1e-12*max(1.,norm(collect(p))) || error("Boundary node moved")
        end
        # Verify every shell facet bounds exactly one tet and every other tet face two.
        counts=Dict{NTuple{3,Int},Int}()
        conn=connectivity[1]
        for k in 1:4:length(conn)
            a,b,c,d=Int.(conn[k:k+3])
            for t in ((a,b,c),(a,b,d),(a,c,d),(b,c,d))
                key=triangle_key(t);counts[key]=get(counts,key,0)+1
            end
        end
        all(v==(k in actual_faces ? 1 : 2) for (k,v) in counts) || error("Nonconforming volume")
        all(get(counts,k,0)==1 for k in actual_faces) || error("Missing shell facet")
        minimum(gmsh.model.mesh.getElementQualities(elements[1],"minSJ"))>0 || error("Invalid tet Jacobian")
        # minSJ alone is 1 for any positive affine tetrahedron, even a sliver.
        inverse_condition=minimum(gmsh.model.mesh.getElementQualities(elements[1],"minSICN"))
        inverse_condition>=1e-6 || error("Near-singular tetrahedra: minimum signed inverse condition=$inverse_condition")
        actual_volume=sum(gmsh.model.mesh.getElementQualities(elements[1],"volume"))
        isapprox(actual_volume,expected_volume;rtol=1e-9) || error("Volume changed")
        gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
        gmsh.option.setNumber("Mesh.Binary",1)
        gmsh.write(output)
        return (boundary_nodes=length(points),nodes=length(gmsh.model.mesh.getNodes()[1]),
                tetrahedra=length(elements[1]),volume=actual_volume,
                minimum_inverse_condition=inverse_condition)
    finally
        gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS) in (5,6,7) || error("plan_prefix lower_z upper_z maximum_size output.msh [quality_threshold] [cap-seeds]")
    length(ARGS)<7 || ARGS[7]=="cap-seeds" || error("Unknown option")
    nd,_=readdlm(ARGS[1]*"-nodes.csv",',',header=true)
    td,_=readdlm(ARGS[1]*"-triangles.csv",',',header=true)
    index=Dict(Int(nd[i,1])=>i for i in axes(nd,1))
    xy=[(Float64(nd[i,2]),Float64(nd[i,3])) for i in axes(nd,1)]
    tri=[Tuple(index[Int(td[i,j])] for j in 2:4) for i in axes(td,1)]
    p,t=extruded_plan_shell(xy,tri,parse.(Float64,ARGS[2:3]))
    println(mesh_fixed_shell(p,t,ARGS[5];maximum_size=parse(Float64,ARGS[4]),
                            optimize_threshold=length(ARGS)>=6 ? parse(Float64,ARGS[6]) : 0.3,
                            seeds=length(ARGS)==7 ? cap_layer_seeds(p,t) : NTuple{3,Float64}[]))
end
