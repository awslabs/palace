# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Experimental connectivity-only repair. No vertex moves and no physical/material
# boundary edge flips: only convex unions of two triangles with identical phase labels.
function repair_plan_triangles(plan, facets, edges, radius; max_passes=12, angle_limit=1.0)
    coordinates = plan.coordinates
    triangles = copy(plan.triangles)
    function label(t)
        phase = triangle_phase(t, coordinates, facets, plan.fabrication_primitives,
                               edges[1].point[3], radius, 1.0e-9radius)
        if hasproperty(plan, :active_segments)
            owner = edge_cluster_owner(edges, plan.active_segments,
                                       triangle_centroid(t, coordinates), phase[2])
            # Different nearest conductors in one trench slot do not define a
            # physical boundary: etched material and SA attributes are identical.
            physical_phase = (phase[1], phase[1] == :metal ? phase[2] : 0)
            return (physical_phase, owner.slot)
        end
        return phase
    end
    labels = [label(t) for t in triangles]
    cross2(a, b, c) = (b[1]-a[1])*(c[2]-a[2]) - (b[2]-a[2])*(c[1]-a[1])
    function min_angle(t)
        p = [coordinates[v] for v in t]
        lengths = [hypot(p[2][1]-p[3][1], p[2][2]-p[3][2]),
                   hypot(p[1][1]-p[3][1], p[1][2]-p[3][2]),
                   hypot(p[1][1]-p[2][1], p[1][2]-p[2][2])]
        angle = 180.0
        for i in 1:3
            a, b, c = lengths[i], lengths[mod1(i+1,3)], lengths[mod1(i+2,3)]
            angle = min(angle, acosd(clamp((b*b+c*c-a*a)/(2b*c), -1.0, 1.0)))
        end
        return angle
    end
    function orient(t)
        return cross2(coordinates[t[1]], coordinates[t[2]], coordinates[t[3]]) > 0 ?
               t : (t[1], t[3], t[2])
    end
    area(t) = abs(cross2(coordinates[t[1]], coordinates[t[2]], coordinates[t[3]]))/2
    function max_length(t)
        return maximum(hypot(coordinates[t[i]][1]-coordinates[t[mod1(i+1,3)]][1],
                             coordinates[t[i]][2]-coordinates[t[mod1(i+1,3)]][2]) for i in 1:3)
    end
    before = minimum(min_angle(t) for t in triangles)
    before_area = sum(area(t) for t in triangles)
    flips = 0
    passes = 0
    for pass in 1:max_passes
        owners = Dict{Tuple{UInt64,UInt64},Vector{Tuple{Int,UInt64}}}()
        for (i,t) in enumerate(triangles), (a,b,c) in ((t[1],t[2],t[3]),
                                                     (t[2],t[3],t[1]),
                                                     (t[3],t[1],t[2]))
            key = a < b ? (a,b) : (b,a)
            push!(get!(owners,key,Tuple{Int,UInt64}[]),(i,c))
        end
        changed = 0
        touched = Set{Int}()
        for (a,b) in sort!(collect(keys(owners)))
            adjacent = owners[(a,b)]
            length(adjacent)==2 || continue
            (i,c),(j,d) = adjacent
            (i in touched || j in touched || labels[i]!=labels[j] || c==d) && continue
            old_min = min(min_angle(triangles[i]), min_angle(triangles[j]))
            old_min < angle_limit || continue
            pc,pd,pa,pb = coordinates[c],coordinates[d],coordinates[a],coordinates[b]
            # The proposed diagonal must split a convex quadrilateral, not cross outside it.
            cross2(pc,pd,pa)*cross2(pc,pd,pb) < 0 || continue
            cross2(pa,pb,pc)*cross2(pa,pb,pd) < 0 || continue
            new_edge = c < d ? (c,d) : (d,c)
            haskey(owners,new_edge) && continue
            first,second = orient((c,d,a)),orient((d,c,b))
            label(first) == labels[i] && label(second) == labels[i] || continue
            old_area = area(triangles[i])+area(triangles[j])
            abs(area(first)+area(second)-old_area) <= 1.0e-10old_area || continue
            # Do not enlarge the local mesh scale while improving the shape tail.
            max(max_length(first),max_length(second)) <=
                1.01max(max_length(triangles[i]),max_length(triangles[j])) || continue
            min(min_angle(first),min_angle(second)) > max(1.01old_min,old_min+1.0e-5) || continue
            triangles[i],triangles[j] = first,second
            push!(touched,i);push!(touched,j);changed+=1
        end
        flips += changed
        passes = pass
        changed==0 && break
    end
    after = minimum(min_angle(t) for t in triangles)
    after_area = sum(area(t) for t in triangles)
    abs(after_area/before_area-1) < 1.0e-11 || error("Plan repair changed total area")
    after >= before-1.0e-10 || error("Plan repair worsened the minimum angle")
    # Boundary edges are unchanged because every accepted flip affects only an interior
    # edge with equal phase labels; outer and material-interface edges cannot be removed.
    println("Plan repair: flips=$flips passes=$passes min_angle=$before -> $after")
    return merge(plan,(triangles=triangles,)),
           (flips=flips, passes=passes, minimum_angle_before=before,
            minimum_angle_after=after, area_before=before_area, area_after=after_area)
end
