# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Experimental tangential relocation of artificial row nodes. Physical boundaries,
# material/slot boundaries, curve intersections, and row-normal coordinates stay fixed.
function relocate_plan_rows(plan, facets, edges, radius; normal_spacing=0.002,
                            angle_limit=1.0, max_passes=8, edge_growth=1.01,
                            infer_interior_axis=false)
    coords=copy(plan.coordinates)
    result=merge(plan,(coordinates=coords,))
    triangles=result.triangles
    tolerance=1.0e-9radius
    cross2(a,b,c)=(b[1]-a[1])*(c[2]-a[2])-(b[2]-a[2])*(c[1]-a[1])
    function min_angle(t)
        p=[coords[v] for v in t]
        lengths=[hypot(p[2][1]-p[3][1],p[2][2]-p[3][2]),
                 hypot(p[1][1]-p[3][1],p[1][2]-p[3][2]),
                 hypot(p[1][1]-p[2][1],p[1][2]-p[2][2])]
        return minimum(acosd(clamp((lengths[mod1(i+1,3)]^2+lengths[mod1(i+2,3)]^2-lengths[i]^2)/
                         (2lengths[mod1(i+1,3)]*lengths[mod1(i+2,3)]),-1.0,1.0)) for i in 1:3)
    end
    function label(t)
        phase=triangle_phase(t,coords,facets,result.fabrication_primitives,
                             edges[1].point[3],radius,tolerance)
        owner=edge_cluster_owner(edges,result.active_segments,triangle_centroid(t,coords),phase[2])
        physical_phase=(phase[1],phase[1]==:metal ? phase[2] : 0)
        return physical_phase,owner.slot
    end
    incident=Dict{UInt64,Vector{Int}}()
    edge_counts=Dict{Tuple{UInt64,UInt64},Int}()
    for (i,t) in enumerate(triangles)
        for a in t
            push!(get!(incident,a,Int[]),i)
        end
        for (a,b) in ((t[1],t[2]),(t[2],t[3]),(t[3],t[1]))
            key=a<b ? (a,b) : (b,a)
            edge_counts[key]=get(edge_counts,key,0)+1
        end
    end
    boundary=Set{UInt64}()
    for ((a,b),count) in edge_counts
        if count==1
            push!(boundary,a);push!(boundary,b)
        end
    end
    labels=[label(t) for t in triangles]
    distance(p)=minimum(point_primitive_distance(p,primitive,tolerance) for primitive in result.fabrication_primitives)
    before=minimum(min_angle(t) for t in triangles)
    moves=0;maximum_move=0.0
    for pass in 1:max_passes
        changed=0
        candidates=Set{UInt64}()
        for t in triangles
            min_angle(t)<angle_limit && union!(candidates,t)
        end
        for vertex in sort!(collect(candidates))
            vertex in boundary && continue
            star=incident[vertex]
            all(labels[i]==labels[star[1]] for i in star) || continue
            original=coords[vertex]
            distance(original)>tolerance || continue
            neighbors=sort!(collect(setdiff(Set(v for i in star for v in triangles[i]),Set([vertex]))))
            same_x=count(v->abs(coords[v][1]-original[1])<tolerance,neighbors)
            same_y=count(v->abs(coords[v][2]-original[2])<tolerance,neighbors)
            axis=same_y>=2 && same_x<2 ? 1 : same_x>=2 && same_y<2 ? 2 : 0
            if axis==0 && infer_interior_axis && same_x<2 && same_y<2
                # A face-interior node has no prescribed row. Choose the tangential
                # direction of an aligned opposite edge in its worst incident triangle.
                worst=star[argmin([min_angle(triangles[i]) for i in star])]
                other=[v for v in triangles[worst] if v!=vertex]
                if abs(coords[other[1]][2]-coords[other[2]][2])<tolerance
                    axis=1
                elseif abs(coords[other[1]][1]-coords[other[2]][1])<tolerance
                    axis=2
                end
            end
            axis==0 && continue
            old_min=minimum(min_angle(triangles[i]) for i in star)
            old_min<angle_limit || continue
            old_edge_max=maximum(hypot(coords[v][1]-original[1],coords[v][2]-original[2]) for v in neighbors)
            choices=Float64[]
            for v in neighbors
                push!(choices,coords[v][axis])
                push!(choices,(coords[v][axis]+original[axis])/2)
            end
            best=original;best_angle=old_min
            for value in unique(choices)
                abs(value-original[axis])<=0.5old_edge_max || continue
                proposed=axis==1 ? (value,original[2]) : (original[1],value)
                distance(proposed)>=min(normal_spacing,distance(original))-tolerance || continue
                minimum(hypot(coords[v][1]-proposed[1],coords[v][2]-proposed[2]) for v in neighbors)>tolerance || continue
                maximum(hypot(coords[v][1]-proposed[1],coords[v][2]-proposed[2]) for v in neighbors)<=edge_growth*old_edge_max || continue
                coords[vertex]=proposed
                if all(cross2(coords[triangles[i][1]],coords[triangles[i][2]],coords[triangles[i][3]])>0 for i in star) &&
                   all(label(triangles[i])==labels[i] for i in star)
                    score=minimum(min_angle(triangles[i]) for i in star)
                    if score>max(best_angle*1.01,best_angle+1.0e-5)
                        best=proposed;best_angle=score
                    end
                end
                coords[vertex]=original
            end
            if best!=original
                coords[vertex]=best;moves+=1;changed+=1
                maximum_move=max(maximum_move,hypot(best[1]-original[1],best[2]-original[2]))
            end
        end
        changed==0 && break
    end
    after=minimum(min_angle(t) for t in triangles)
    after>=before-1.0e-10 || error("Relocation worsened minimum angle")
    println("Row relocation: moves=$moves max_move=$maximum_move min_angle=$before -> $after")
    return result,(moves=moves,maximum_move=maximum_move,minimum_angle_before=before,minimum_angle_after=after)
end
