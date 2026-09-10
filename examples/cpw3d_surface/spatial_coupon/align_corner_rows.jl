# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Add matching tangential stations on adjacent artificial rows near skinny triangle
# fans, then flip interior diagonals. This addresses row-node staggering without
# moving any original vertex, coarsening the 2 nm boundary layer, or moving a boundary.
function align_corner_rows(plan, facets, edges, radius; angle_limit=1.0, spacing=0.002)
    coords=copy(plan.coordinates);triangles=copy(plan.triangles)
    tolerance=1.0e-9radius
    edgekey(a,b)=a<b ? (a,b) : (b,a)
    cross2(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
    orient(t)=cross2(coords[t[1]],coords[t[2]],coords[t[3]])>0 ? t : (t[1],t[3],t[2])
    function angle(t)
        p=[coords[v] for v in t]
        ds=[hypot(p[2][1]-p[3][1],p[2][2]-p[3][2]),
            hypot(p[1][1]-p[3][1],p[1][2]-p[3][2]),
            hypot(p[1][1]-p[2][1],p[1][2]-p[2][2])]
        return minimum(acosd(clamp((ds[mod1(i+1,3)]^2+ds[mod1(i+2,3)]^2-ds[i]^2)/
                           (2ds[mod1(i+1,3)]*ds[mod1(i+2,3)]),-1.0,1.0)) for i in 1:3)
    end
    function label(t)
        phase=triangle_phase(t,coords,facets,plan.fabrication_primitives,
                             edges[1].point[3],radius,tolerance)
        owner=edge_cluster_owner(edges,plan.active_segments,triangle_centroid(t,coords),phase[2])
        return (phase[1],phase[1]==:metal ? phase[2] : 0),owner.slot
    end
    owners=Dict{Tuple{UInt64,UInt64},Set{Int}}()
    function add_triangle(i,t)
        for (a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
            push!(get!(owners,edgekey(a,b),Set{Int}()),i)
        end
    end
    function remove_triangle(i,t)
        for (a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
            delete!(owners[edgekey(a,b)],i)
        end
    end
    for (i,t) in enumerate(triangles)
        add_triangle(i,t)
    end
    rows=Dict{Tuple{Int,Int},Vector{Tuple{UInt64,UInt64}}}()
    for (a,b) in keys(owners)
        for axis in 1:2
            normal=3-axis
            if abs(coords[a][normal]-coords[b][normal])<tolerance
                push!(get!(rows,(axis,round(Int,coords[a][normal]/tolerance)),Tuple{UInt64,UInt64}[]),(a,b))
            end
        end
    end
    candidates=Dict{Tuple{UInt64,UInt64},Vector{NTuple{2,Float64}}}()
    original_min=minimum(angle(t) for t in triangles)
    for t in triangles
        angle(t)<angle_limit || continue
        for (a,b,c) in ((t[1],t[2],t[3]),(t[2],t[3],t[1]),(t[3],t[1],t[2]))
            hypot(coords[a][1]-coords[b][1],coords[a][2]-coords[b][2])<10spacing || continue
            for axis in 1:2
                normal=3-axis
                abs(coords[a][normal]-coords[b][normal])<tolerance || continue
                for vertex in (a,b)
                    p=axis==1 ? (coords[vertex][1],coords[c][2]) : (coords[c][1],coords[vertex][2])
                    for edge in get(rows,(axis,round(Int,p[normal]/tolerance)),Tuple{UInt64,UInt64}[])
                        x,y=coords[edge[1]][axis],coords[edge[2]][axis]
                        min(x,y)+spacing-tolerance<=p[axis]<=max(x,y)-spacing+tolerance || continue
                        length(owners[edge])==2 || continue
                        adjacent=collect(owners[edge])
                        label(triangles[adjacent[1]])==label(triangles[adjacent[2]]) || continue
                        push!(get!(candidates,edge,NTuple{2,Float64}[]),p)
                        break
                    end
                end
            end
        end
    end
    added=0;next_id=maximum(keys(coords))+UInt64(1)
    for (a,b) in sort!(collect(keys(candidates)))
        axis=abs(coords[a][1]-coords[b][1])>abs(coords[a][2]-coords[b][2]) ? 1 : 2
        direction=coords[b][axis]>coords[a][axis] ? 1 : -1
        points=sort!(unique(candidates[(a,b)]);by=p->direction*p[axis])
        previous=a
        for p in points
            hypot(p[1]-coords[previous][1],p[2]-coords[previous][2])>=spacing-tolerance || continue
            hypot(p[1]-coords[b][1],p[2]-coords[b][2])>=spacing-tolerance || continue
            adjacent=collect(get(owners,edgekey(previous,b),Set{Int}()))
            length(adjacent)==2 || continue
            node=next_id;coords[node]=p
            replacements=[]
            valid=true
            for i in adjacent
                t=triangles[i];c=only([v for v in t if v!=previous && v!=b])
                one,two=orient((previous,node,c)),orient((node,b,c))
                if label(one)!=label(t) || label(two)!=label(t)
                    valid=false;break
                end
                push!(replacements,(i,t,one,two))
            end
            if !valid
                delete!(coords,node);continue
            end
            for (i,t,one,two) in replacements
                remove_triangle(i,t);triangles[i]=one;add_triangle(i,one)
                push!(triangles,two);add_triangle(length(triangles),two)
            end
            previous=node;next_id+=1;added+=1
        end
    end
    enriched=merge(plan,(coordinates=coords,triangles=triangles,node_tags=sort!(collect(keys(coords)))))
    repaired,flip_stats=repair_plan_triangles(enriched,facets,edges,radius;max_passes=40,angle_limit=angle_limit)
    println("Corner row alignment: added_nodes=$added triangles=$(length(plan.triangles)) -> $(length(triangles)) original_min_angle=$original_min final_min_angle=$(flip_stats.minimum_angle_after)")
    return repaired,(added_nodes=added,triangles_before=length(plan.triangles),triangles_after=length(triangles),minimum_angle_before=original_min,flips=flip_stats)
end
