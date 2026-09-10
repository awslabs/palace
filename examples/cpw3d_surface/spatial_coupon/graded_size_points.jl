# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Exact min_i(size_i + distance(point, point_i)), accelerated by a bounding tree.
# Bounds include the minimum source size; unweighted nearest-neighbor lookup
# would be incorrect when the requested sizes differ.
struct SizePointNode
    lower::NTuple{3,Float64}
    upper::NTuple{3,Float64}
    minimum_size::Float64
    left::Int
    right::Int
    first::Int
    last::Int
end
struct SizePointTree
    points::Vector{NTuple{4,Float64}}
    order::Vector{Int}
    nodes::Vector{SizePointNode}
end
function SizePointTree(points)
    data=NTuple{4,Float64}[Tuple(Float64.(p)) for p in points]
    all(all(isfinite,p) && p[4]>0 for p in data) || error("Invalid refinement point")
    tree=SizePointTree(data,collect(eachindex(data)),SizePointNode[])
    function build(first,last)
        lower=ntuple(d->minimum(data[tree.order[i]][d] for i in first:last),3)
        upper=ntuple(d->maximum(data[tree.order[i]][d] for i in first:last),3)
        size=minimum(data[tree.order[i]][4] for i in first:last)
        index=length(tree.nodes)+1
        push!(tree.nodes,SizePointNode(lower,upper,size,0,0,first,last))
        if last-first>=8
            axis=argmax(ntuple(d->upper[d]-lower[d],3))
            sort!(view(tree.order,first:last);by=i->data[i][axis])
            middle=(first+last)÷2
            left=build(first,middle);right=build(middle+1,last)
            tree.nodes[index]=SizePointNode(lower,upper,size,left,right,first,last)
        end
        return index
    end
    isempty(data) || build(1,length(data))
    return tree
end
function size_point_bound(node,point)
    return node.minimum_size+sqrt(sum(max(node.lower[d]-point[d],0.,point[d]-node.upper[d])^2 for d in 1:3))
end
function query_size_point(tree::SizePointTree,point,best,node_index=1)
    isempty(tree.nodes) && return best
    node=tree.nodes[node_index]
    size_point_bound(node,point)>=best && return best
    if node.left==0
        @inbounds for i in node.first:node.last
            p=tree.points[tree.order[i]]
            best=min(best,p[4]+sqrt(sum((p[d]-point[d])^2 for d in 1:3)))
        end
    else
        first,second=node.left,node.right
        if size_point_bound(tree.nodes[first],point)>size_point_bound(tree.nodes[second],point)
            first,second=second,first
        end
        best=query_size_point(tree,point,best,first)
        best=query_size_point(tree,point,best,second)
    end
    return best
end
