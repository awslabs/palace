# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Exclusive geometric ownership and a conservative per-element certificate.
# Distance to a segment/arc, and the minimum of such distances, are 1-Lipschitz.
# A competing-label distance margin > 2*enclosing_radius therefore certifies
# the selected label on the entire element, not just at sampled points.
include(joinpath(@__DIR__,"ownership_bernstein.jl"))

function build_interface_ownership(edges,loops,radius;
                                   fabricated=false,metal_thickness=0.1,overetch=0.05)
    tolerance=1e-9radius
    grouped=Dict{Tuple{Int,Float64},Vector{Float64}}()
    for edge in edges
        sign=hasproperty(edge,:normal_sign) ? edge.normal_sign : 1.
        push!(get!(grouped,(round(Int,edge.point[3]/tolerance),sign),Float64[]),edge.point[3])
    end
    layers=sort!([(sum(z)/length(z),key[2]) for (key,z) in grouped])
    bands=[fabricated ? minmax(z-sign*overetch,z+sign*metal_thickness) : (z,z)
           for (z,sign) in layers]
    for i in eachindex(bands),j in (i+1):length(bands)
        max(bands[i][1],bands[j][1])<=min(bands[i][2],bands[j][2])+tolerance &&
            error("Interface layer bands overlap; ownership would be ambiguous")
    end
    layer_edges=[[e for e in edges if abs(e.point[3]-plane)<=tolerance] for (plane,_) in layers]
    primitives=NamedTuple[]
    for loop in loops, primitive in physical_segments([loop],0.,tolerance)
        push!(primitives,(;primitive...,conductor=loop.conductor,plane=loop.plane))
    end
    edge_segments=Dict(e=>ownership_edge_segment(e,radius) for e in edges)
    layer_primitives=[[p for p in primitives if abs(p.plane-plane)<=tolerance] for (plane,_) in layers]
    function layer_at(point)
        distances=[max(lo-point[3],0.,point[3]-hi) for (lo,hi) in bands]
        best=minimum(distances)
        best<=10tolerance || error("Interface point $point lies outside process layer bands $bands")
        candidates=findall(d->d<=best+tolerance,distances)
        length(candidates)==1 || error("Ambiguous interface layer")
        return only(candidates)
    end
    function sa_selection(point,layer)
        ps=layer_primitives[layer]
        isempty(ps) && return layer_edges[layer],0,Float64[]
        distances=[point_primitive_distance((point[1],point[2]),p,tolerance) for p in ps]
        best=minimum(distances)
        best>3radius+tolerance && return layer_edges[layer],0,distances
        conductor=minimum(ps[i].conductor for i in eachindex(ps) if distances[i]<=best+tolerance)
        return [e for e in layer_edges[layer] if e.conductor==conductor],conductor,distances
    end
    function candidates_at(attribute,point,layer)
        family=div(attribute,1000)
        family==3 && return first(sa_selection(point,layer))
        family in (4,5,6) || error("Unsupported physical surface family $attribute")
        conductor=mod(attribute,100)
        candidates=[e for e in layer_edges[layer] if e.conductor==conductor]
        isempty(candidates) && error("No edge signature for CAD conductor $conductor")
        return candidates
    end
    function classify(attribute,point)
        layer=layer_at(point)
        owner=nearest_edge(candidates_at(attribute,point,layer),point,radius)
        family=div(attribute,1000)
        return family==3 ? (attribute>=3100 ? 3100 : 3000)+owner.slot :
               metal_surface_attribute(1000family,owner.slot,mod(attribute,100))
    end
    function stable_slots(candidates,point,r,slot)
        own=Inf;other=Inf
        for edge in candidates
            distance=segment_distance(edge,point,radius)
            if edge.slot==slot;own=min(own,distance);else;other=min(other,distance);end
        end
        return isfinite(own) && other-own>2r+4tolerance
    end
    function certify(attribute,center,hull)
        layer=layer_at(center)
        lo,hi=bands[layer]
        all(lo-10tolerance<=p[3]<=hi+10tolerance for p in hull) || return false
        r=maximum(sqrt(sum((p[d]-center[d])^2 for d in 1:3)) for p in hull)
        candidates=candidates_at(attribute,center,layer)
        slot=nearest_edge(candidates,center,radius).slot
        if div(attribute,1000)==3
            # Conductor selection is a first-level nearest-primitive partition.
            # Certify it too, unless every candidate on this layer has one slot.
            length(unique(e.slot for e in layer_edges[layer]))==1 && return true
            _,conductor,distances=sa_selection(center,layer)
            if !isempty(distances)
                ps=layer_primitives[layer]
                best=minimum(distances)
                if conductor==0
                    best-r>3radius+tolerance || return false
                else
                    best+r<3radius-tolerance || return false
                    own=minimum(distances[i] for i in eachindex(ps) if ps[i].conductor==conductor)
                    other=minimum((distances[i] for i in eachindex(ps) if ps[i].conductor!=conductor);init=Inf)
                    if !(other-own>2r+4tolerance)
                        all(p.kind==:line for p in ps) || return false
                        segments=[ownership_segment((p.first[1],p.first[2],0.),(p.last[1],p.last[2],0.)) for p in ps]
                        xy_center=(center[1],center[2],0.)
                        xy_hull=[(p[1],p[2],0.) for p in hull]
                        certify_segment_label(segments,[p.conductor for p in ps],conductor,xy_center,xy_hull,tolerance) || return false
                    end
                end
            end
        end
        return stable_slots(candidates,center,r,slot) ||
               certify_segment_label([edge_segments[e] for e in candidates],
                                     [e.slot for e in candidates],slot,center,hull,tolerance)
    end
    return (classify=classify,certify=certify)
end

# A quadratic triangular map lies in the convex hull of its Bernstein control
# points. Lagrange mid-edge nodes alone do not provide that bound.
function interface_triangle_hull(points)
    length(points)==3 && return points
    if length(points)==6
        hull=copy(points[1:3])
        for (mid,a,b) in ((4,1,2),(5,2,3),(6,3,1))
            push!(hull,ntuple(d->2points[mid][d]-(points[a][d]+points[b][d])/2,3))
        end
        return hull
    end
    return nothing # No certificate for unsupported geometric orders.
end
