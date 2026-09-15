# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

include(joinpath(@__DIR__,"interface_ownership.jl"))

# Physical material/interface families are retained from CAD. Coarse per-element
# slot labels are non-authoritative visualization/routing tags when a triangle
# crosses the response-ownership partition. Science uses the positive-weight
# quadrature-point partition recorded below; the whole-element Lipschitz result is
# retained separately as an ambiguity diagnostic.
function label_interface_patches(edges,loops,radius,report_path;minimum_size=0.0,
                                 fabricated=false,metal_thickness=0.1,overetch=0.05,
                                 ownership_coordinates=identity)
    ownership=build_interface_ownership(edges,loops,radius;fabricated=fabricated,
                                        metal_thickness=metal_thickness,overetch=overetch)
    classify=ownership.classify
    tags,xyz,_=gmsh.model.mesh.getNodes()
    coords=Dict(tag=>(xyz[3i-2],xyz[3i-1],xyz[3i]) for (i,tag) in enumerate(tags))
    groups=Dict{Int,Vector{Tuple{Int,UInt64,Vector{UInt64}}}}()
    areas=Dict{Int,Float64}();ambiguous=Dict{Int,Float64}();counts=Dict{Int,Int}()
    unresolved_area=Dict{Int,Float64}();unresolved_count=Dict{Int,Int}()
    quadrature_area=Dict{Int,Float64}();quadrature_whole=0.;quadrature_points=0
    refinement_points=NTuple{4,Float64}[]
    certificates=Tuple{UInt64,Int,Bool,NTuple{3,UInt64}}[]
    old_groups=Tuple{Int32,Int32}[];old_entities=Set{Int32}()
    for (_,attribute) in gmsh.model.getPhysicalGroups(2)
        attribute==1 && continue
        push!(old_groups,(Int32(2),Int32(attribute)))
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2,attribute)
            entity in old_entities && error("A CAD face has multiple physical interface assignments")
            push!(old_entities,entity)
            types,element_tags,connectivity=gmsh.model.mesh.getElements(2,entity)
            for (type,etags,enodes) in zip(types,element_tags,connectivity)
                name,_,_,nnode,_,primary=gmsh.model.mesh.getElementProperties(type)
                startswith(name,"Triangle") && primary==3 || error("Interface labeling requires triangular faces")
                quadrature_rule="Gauss4"
                integration_points,integration_weights=gmsh.model.mesh.getIntegrationPoints(type,quadrature_rule)
                all(weight>0 for weight in integration_weights) ||
                    error("Response ownership requires positive quadrature weights")
                _,jacobian_measures,quadrature_coordinates=
                    gmsh.model.mesh.getJacobians(type,integration_points,entity)
                nq=length(integration_weights)
                length(jacobian_measures)==nq*length(etags) || error("Unexpected interface Jacobian data")
                for (i,etag) in enumerate(etags)
                    nodes=collect(enodes[(i-1)*nnode+1:i*nnode])
                    points=[coords[node] for node in nodes]
                    center=ntuple(d->sum(points[j][d] for j in 1:3)/3,3)
                    ownership_center=ownership_coordinates(center)
                    target=classify(attribute,ownership_center)
                    ab=[points[2][d]-points[1][d] for d in 1:3]
                    ac=[points[3][d]-points[1][d] for d in 1:3]
                    corner_area=norm(cross(ab,ac))/2
                    corner_area>0 || error("Degenerate interface triangle")
                    area=sum(integration_weights[q]*jacobian_measures[(i-1)*nq+q] for q in 1:nq)
                    area>0 || error("Nonpositive integrated interface area")
                    quadrature_whole+=area
                    for q in 1:nq
                        coordinate_index=3*((i-1)*nq+q-1)
                        point=(quadrature_coordinates[coordinate_index+1],
                               quadrature_coordinates[coordinate_index+2],
                               quadrature_coordinates[coordinate_index+3])
                        owner=classify(attribute,ownership_coordinates(point))
                        measure=integration_weights[q]*jacobian_measures[(i-1)*nq+q]
                        measure>0 || error("Nonpositive response-ownership measure")
                        quadrature_area[owner]=get(quadrature_area,owner,0.)+measure
                        quadrature_points+=1
                    end
                    areas[target]=get(areas,target,0.)+area
                    counts[target]=get(counts,target,0)+1
                    samples=[ntuple(d->0.8p[d]+0.2center[d],3) for p in points[1:3]]
                    if any(classify(attribute,ownership_coordinates(p))!=target for p in samples)
                        ambiguous[target]=get(ambiguous,target,0.)+area
                    end
                    hull=interface_triangle_hull(points)
                    ownership_hull = hull === nothing ? nothing : ownership_coordinates.(hull)
                    certified=ownership_hull!==nothing &&
                              ownership.certify(attribute,ownership_center,ownership_hull)
                    push!(certificates,(etag,target,certified,Tuple(nodes[1:3])))
                    if !certified
                        unresolved_area[target]=get(unresolved_area,target,0.)+area
                        unresolved_count[target]=get(unresolved_count,target,0)+1
                        target_size=max(minimum_size,sqrt(area)/4)
                        push!(refinement_points,(center[1],center[2],center[3],target_size))
                    end
                    push!(get!(groups,target,Tuple{Int,UInt64,Vector{UInt64}}[]),(Int(type),etag,nodes))
                end
            end
        end
    end
    # No geometry, volume connectivity, boundary connectivity, or element ID changes.
    for entity in old_entities
        gmsh.model.mesh.removeElements(2,entity)
    end
    gmsh.model.removePhysicalGroups(old_groups)
    for attribute in sort!(collect(keys(groups)))
        entity=gmsh.model.addDiscreteEntity(2)
        for type in unique(r[1] for r in groups[attribute])
            records=[r for r in groups[attribute] if r[1]==type]
            gmsh.model.mesh.addElementsByType(entity,type,[r[2] for r in records],reduce(vcat,[r[3] for r in records]))
        end
        gmsh.model.addPhysicalGroup(2,[entity],attribute,"surface_$attribute")
    end
    owned_measure=sum(values(quadrature_area))
    closure_tolerance=1e-12
    relative_closure=abs(owned_measure-quadrature_whole)/quadrature_whole
    relative_closure<=closure_tolerance || error("Response-ownership quadrature does not close")
    open(report_path,"w") do f
        println(f,"attribute,elements,area,ambiguous_area,ambiguous_fraction,unresolved_elements,unresolved_area,unresolved_fraction,quadrature_rule,quadrature_order,quadrature_points,quadrature_whole_measure,quadrature_owned_measure,quadrature_relative_closure,quadrature_closure_tolerance,quadrature_unmatched,quadrature_overlaps,quadrature_positive_weights")
        for attribute in sort!(collect(keys(areas)))
            a=get(ambiguous,attribute,0.);u=get(unresolved_area,attribute,0.)
            println(f,"$attribute,$(counts[attribute]),$(areas[attribute]),$a,$(a/areas[attribute]),$(get(unresolved_count,attribute,0)),$u,$(u/areas[attribute]),Gauss4,4,$quadrature_points,$quadrature_whole,$owned_measure,$relative_closure,$closure_tolerance,0,0,1")
        end
    end
    open(report_path*".quadrature.csv","w") do f
        println(f,"attribute,measure")
        for attribute in sort!(collect(keys(quadrature_area)))
            println(f,"$attribute,$(quadrature_area[attribute])")
        end
    end
    open(report_path*".refine.csv","w") do f
        println(f,"x,y,z,size")
        for p in refinement_points;println(f,join(p,","));end
    end
    open(report_path*".elements.csv","w") do f
        println(f,"element,attribute,certified,node_a,node_b,node_c")
        for (element,attribute,certified,nodes) in sort!(certificates;by=first)
            println(f,"$element,$attribute,$(Int(certified)),$(nodes[1]),$(nodes[2]),$(nodes[3])")
        end
    end
    println("Element-wise interface slots: ",join(sort!(collect(keys(areas))),","),
            "; certified elements=",count(x->x[3],certificates),"/",length(certificates))
    return nothing
end
