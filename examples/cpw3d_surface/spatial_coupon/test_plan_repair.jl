# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
include(joinpath(@__DIR__, "repair_plan_triangles.jl"))
include(joinpath(@__DIR__, "align_corner_rows.jl"))

triangle_phase(t,c,facets,p,plane,r,tol) = isempty(facets) ? (:vacuum,0) : facets[t]
triangle_centroid(t,c) = (sum(c[v][1] for v in t)/3,sum(c[v][2] for v in t)/3)
edge_cluster_owner(args...) = (slot=1,)
function boundary_edges(triangles)
    counts=Dict{Tuple{UInt64,UInt64},Int}()
    for t in triangles, (a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
        k=a<b ? (a,b) : (b,a);counts[k]=get(counts,k,0)+1
    end
    @test maximum(values(counts))<=2
    return Set(k for (k,n) in counts if n==1)
end
@testset "Plan repair invariants" begin
    coords=Dict{UInt64,NTuple{2,Float64}}(1=>(0.,0.),2=>(1.,0.),3=>(1.,1.),4=>(.001,.01))
    tris=[(UInt64(1),UInt64(2),UInt64(3)),(UInt64(1),UInt64(3),UInt64(4))]
    plan=(coordinates=coords,triangles=tris,fabrication_primitives=[])
    edges=[(point=(0.,0.,0.),)]
    p,m=repair_plan_triangles(plan,[],edges,1.0)
    @test m.flips==1
    @test m.minimum_angle_after>m.minimum_angle_before
    @test p.coordinates==coords
    @test boundary_edges(p.triangles)==boundary_edges(tris)
    labels=Dict(tris[1]=>(:vacuum,0),tris[2]=>(:metal,1))
    p2,m2=repair_plan_triangles(plan,labels,edges,1.0)
    @test m2.flips==0
    @test p2.triangles==tris

    # A refined physical boundary opposite a sparsely sampled internal row creates
    # skinny fans. Enrichment must preserve every old point and outer boundary edge.
    xs=[0.,.002,.004,.006,.008,.01,.1]
    nodes=Dict{UInt64,NTuple{2,Float64}}(UInt64(i)=>(x,0.) for (i,x) in enumerate(xs))
    nodes[8]=(0.,.002);nodes[9]=(.1,.002);nodes[10]=(0.,.004);nodes[11]=(.1,.004)
    triangles=NTuple{3,UInt64}[(UInt64(i),UInt64(i+1),UInt64(9)) for i in 1:6]
    append!(triangles,[(UInt64(1),UInt64(9),UInt64(8)),
                      (UInt64(8),UInt64(9),UInt64(11)),
                      (UInt64(8),UInt64(11),UInt64(10))])
    plan=(coordinates=nodes,triangles=triangles,node_tags=sort!(collect(keys(nodes))),
          fabrication_primitives=[],active_segments=[])
    p,m=align_corner_rows(plan,[],edges,1.0)
    @test m.added_nodes>0
    @test all(p.coordinates[k]==v for (k,v) in nodes)
    @test boundary_edges(p.triangles)==boundary_edges(triangles)
    @test m.flips.minimum_angle_after>=m.minimum_angle_before-1e-10
end
