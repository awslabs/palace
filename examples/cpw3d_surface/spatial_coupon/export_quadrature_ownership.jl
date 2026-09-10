# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Export the SAME polygonal ownership definition used by element labeling, not
# a nearest infinite-line or chord approximation. One group has one process plane.
include(joinpath(@__DIR__,"mesh_spatial_coupon.jl"))
using TOML
using Random
length(ARGS) in (2,3) || error("signature_directory output.csv [validation_queries.csv]")
root,output=abspath.(ARGS[1:2])
isfile(output) && error("Refuse to overwrite ownership data")
process=TOML.parsefile(joinpath(root,"process.toml"))
process["Units"]=="um" || error("This exporter expects micrometre input geometry")
radius=Float64(process["Radius"])
edges=read_edges(joinpath(root,"mesh-signature.csv"))
loops=read_boundary(joinpath(root,"plan-view-boundary.csv"))
length(unique(e.point[3] for e in edges))==1 ||
    error("Export layers separately with distinct interface attributes; shared-layer groups are unsupported")
primitives=NamedTuple[]
for loop in loops,p in physical_segments([loop],0.,1e-9radius)
    p.kind==:line || error("Quadrature ownership currently requires polygonal plan-view boundaries")
    push!(primitives,(;p...,conductor=loop.conductor,plane=loop.plane))
end
open(output,"w") do f
    println(f,"Group,Conductor,Slot,Role,Radius,X0,Y0,Z0,X1,Y1,Z1")
    for p in primitives
        println(f,join((0,p.conductor,0,"Boundary",radius,p.first...,p.plane,p.last...,p.plane),","))
    end
    for edge in edges
        first,last=extended_interval(edge,radius)
        a=add(edge.point,scale(first,edge.tangent));b=add(edge.point,scale(last,edge.tangent))
        for group in (0,edge.conductor)
            println(f,join((group,edge.conductor,edge.slot,"Signature",radius,a...,b...),","))
        end
    end
end
if length(ARGS)==3
    include(joinpath(@__DIR__,"interface_ownership.jl"))
    ownership=build_interface_ownership(edges,loops,radius;fabricated=true,
        metal_thickness=Float64(process["MetalThickness"]),overetch=Float64(process["Overetch"]))
    lower,upper=coupon_bounds(edges,radius,Float64(process["MetalThickness"]),Float64(process["Overetch"]))
    rng=MersenneTwister(20260910)
    groups=vcat([0],sort!(unique(e.conductor for e in edges)))
    isfile(ARGS[3]) && error("Refuse to overwrite query file")
    open(ARGS[3],"w") do f
        println(f,"Group,Slot,X,Y,Z")
        for _ in 1:1000
            x=lower[1]+rand(rng)*(upper[1]-lower[1])
            y=lower[2]+rand(rng)*(upper[2]-lower[2])
            z=edges[1].point[3]+edges[1].normal_sign*0.5Float64(process["MetalThickness"])
            for group in groups
                attribute=group==0 ? 3100 : 5000+group
                result=ownership.classify(attribute,(x,y,z))
                slot=group==0 ? result-3100 : (result-5000-group)÷100
                println(f,join((group,slot,x,y,z),","))
            end
        end
    end
end
