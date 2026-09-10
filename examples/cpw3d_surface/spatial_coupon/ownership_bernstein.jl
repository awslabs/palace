# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Sharper sufficient ownership test for segment families. Composing squared
# distance to an affine subspace with a linear/quadratic triangle map gives a
# degree-2/4 polynomial. Its Bernstein coefficients bound it over the triangle.
# Unknown segment-projection regimes use a point upper bound for the selected
# segment, and an infinite-line lower bound for competitors: never an unsafe
# continuation of a clamped distance formula.
function ownership_product_weights(alphas)
    degree=sum(first(alphas))
    multi(a)=factorial(sum(a))/prod(factorial(x) for x in a)
    sums=sort!(unique(ntuple(d->a[d]+b[d],3) for a in alphas for b in alphas))
    pairs=Tuple{Int,Int,Int,Float64}[]
    for (i,a) in enumerate(alphas),(j,b) in enumerate(alphas)
        c=ntuple(d->a[d]+b[d],3)
        push!(pairs,(i,j,findfirst(==(c),sums),multi(a)*multi(b)/multi(c)))
    end
    return (pairs=pairs,count=length(sums),degree=2degree)
end
const OWNERSHIP_WEIGHTS_LINEAR=ownership_product_weights([(1,0,0),(0,1,0),(0,0,1)])
const OWNERSHIP_WEIGHTS_QUADRATIC=ownership_product_weights([(2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1)])

function ownership_segment(a,b)
    direction=ntuple(d->b[d]-a[d],3)
    norm2=sum(x*x for x in direction)
    return (first=a,last=b,direction=direction,norm2=norm2)
end
function ownership_edge_segment(edge,radius)
    first,last=extended_interval(edge,radius)
    return ownership_segment(add(edge.point,scale(first,edge.tangent)),
                             add(edge.point,scale(last,edge.tangent)))
end
function ownership_residuals(segment,center,hull,upper)
    segment.norm2>0 || return [ntuple(d->p[d]-segment.first[d],3) for p in hull]
    projection(p)=sum((p[d]-segment.first[d])*segment.direction[d] for d in 1:3)/segment.norm2
    values=projection.(hull)
    fixed=nothing
    if maximum(values)<=0
        fixed=segment.first
    elseif minimum(values)>=1
        fixed=segment.last
    elseif upper && !(minimum(values)>=0 && maximum(values)<=1)
        t=clamp(projection(center),0.,1.)
        fixed=ntuple(d->segment.first[d]+t*segment.direction[d],3)
    end
    if fixed!==nothing
        return [ntuple(d->p[d]-fixed[d],3) for p in hull]
    end
    return [ntuple(d->p[d]-segment.first[d]-projection(p)*segment.direction[d],3) for p in hull]
end
function certify_segment_label(segments,labels,target,center,hull,tolerance)
    weights=length(hull)==3 ? OWNERSHIP_WEIGHTS_LINEAR :
            length(hull)==6 ? OWNERSHIP_WEIGHTS_QUADRATIC : nothing
    weights===nothing && return false
    selected=findall(==(target),labels)
    isempty(selected) && return false
    competitors=findall(!=(target),labels)
    isempty(competitors) && return true
    lower=[ownership_residuals(segments[j],center,hull,false) for j in competitors]
    # The distance-tie allowance is included in the squared-distance margin.
    scale=max(1.,maximum(sqrt(sum(x*x for x in p)) for residuals in lower for p in residuals))
    margin=8tolerance*scale
    for i in selected
        upper=ownership_residuals(segments[i],center,hull,true)
        proven=true
        for other in lower
            coefficients=zeros(weights.count)
            for (a,b,c,w) in weights.pairs
                coefficients[c]+=w*(dot(upper[a],upper[b])-dot(other[a],other[b]))
            end
            if maximum(coefficients)>=-margin
                proven=false;break
            end
        end
        proven && return true
    end
    return false
end
