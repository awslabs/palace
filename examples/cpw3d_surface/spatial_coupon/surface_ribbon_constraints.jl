# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Explicit in-surface rows: geometry constrains transverse spacing, while a
# separate scalar target may keep the smooth edge-tangent direction coarser.
# Rows are clipped to existing CAD faces before fragmentation. They are not
# extruded into the far volume and do not alter any material boundary.
function ribbon_distances(first_size, count; ratio=1.5)
    first_size > 0 && count >= 0 && ratio > 1 || error("Invalid ribbon spacing")
    result = Float64[]
    distance = 0.0
    step = first_size
    for _ in 1:count
        distance += step
        push!(result, distance)
        step *= ratio
    end
    return result
end

function surface_ribbon_lines(occ, loops, planes, fabricated, metal_thickness,
                              overetch, first_size, count, tolerance)
    isempty(loops) && error("Surface ribbons require explicit polygonal mask boundaries")
    any(!isempty(circular_arc_runs(loop.points, tolerance)) for loop in loops) &&
        error("Curved ribbons require a curved offset implementation")
    occ.synchronize()
    plane_surfaces = Dict{Float64,Vector{Tuple{Int32,Int32}}}()
    vertical_surfaces = Tuple{Int32,Int32}[]
    for (_, surface) in gmsh.model.getEntities(2)
        center=collect(occ.getCenterOfMass(2,surface))
        uv=gmsh.model.getParametrization(2,surface,center)
        normal=gmsh.model.getNormal(surface,uv)
        lo,hi=gmsh.model.getParametrizationBounds(2,surface)
        for a in (.2,.5,.8),b in (.2,.5,.8)
            point=gmsh.model.getValue(2,surface,[lo[1]+a*(hi[1]-lo[1]),lo[2]+b*(hi[2]-lo[2])])
            abs(dot(normal,point-center))<=10tolerance || error("Ribbon scout supports planar geometry only")
        end
        if abs(normal[3])>1-1e-10
            z=Float64(center[3])
            push!(get!(plane_surfaces,z,Tuple{Int32,Int32}[]),(2,surface))
        else
            push!(vertical_surfaces,(2,surface))
        end
    end
    distances=ribbon_distances(first_size,count)
    tools=Tuple{Int32,Int32}[]
    function clipped_lines(points,z,surfaces; classes=nothing)
        isempty(points) && return
        lines=Tuple{Int32,Int32}[]
        vertices=[occ.addPoint(p[1],p[2],z) for p in points]
        for i in eachindex(points)
            classes !== nothing && classes[i] != "Physical" && continue
            a,b=points[i],points[mod1(i+1,length(points))]
            hypot(a[1]-b[1],a[2]-b[2])>tolerance || continue
            push!(lines,(1,occ.addLine(vertices[i],vertices[mod1(i+1,length(points))])))
        end
        if isempty(surfaces) || isempty(lines)
            isempty(lines) || occ.remove(lines,true)
            return
        end
        clipped,_=occ.intersect(lines,occ.copy(surfaces),-1,true,true)
        append!(tools,[(dim,tag) for (dim,tag) in clipped if dim==1])
        return nothing
    end
    for (plane,sign) in planes
        local_loops=[loop for loop in loops if abs(loop.plane-plane)<=tolerance]
        horizontal = fabricated ? [(plane,1.),(plane+sign*metal_thickness,1.),
                                   (plane-sign*overetch,-1.)] : [(plane,1.),(plane,-1.)]
        for (z,side) in horizontal
            surfaces=reduce(vcat,[s for (height,s) in plane_surfaces if abs(height-z)<=tolerance];init=Tuple{Int32,Int32}[])
            for loop in local_loops,distance in distances
                shifted=loop.hole ? offset_hole_points(loop,side*distance,tolerance) :
                                    offset_loop_points(loop,side*distance,tolerance)
                clipped_lines(shifted,z,surfaces)
            end
        end
        if fabricated
            for loop in local_loops, distance in distances
                if distance < metal_thickness/2
                    clipped_lines(loop.points,plane+sign*distance,vertical_surfaces;classes=loop.classes)
                    clipped_lines(loop.points,plane+sign*(metal_thickness-distance),vertical_surfaces;classes=loop.classes)
                end
                if 0<distance<overetch/2
                    clipped_lines(loop.points,plane-sign*distance,vertical_surfaces;classes=loop.classes)
                    clipped_lines(loop.points,plane-sign*(overetch-distance),vertical_surfaces;classes=loop.classes)
                end
            end
        end
    end
    isempty(tools) && error("No valid surface ribbon constraints were produced")
    println("Surface ribbon constraints: $(length(tools)) clipped row segments; distances=$distances")
    flush(stdout)
    return tools
end
