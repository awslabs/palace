# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Separate bounded tensor-metric scout. Unlike a scalar size callback, MathEvalAniso
# supplies the complete metric to BAMG/MMG3D. No production defaults are changed.
# Gmsh field references: https://gmsh.info/doc/texinfo/gmsh.html#index-MathEvalAniso
# and https://gmsh.info/doc/texinfo/gmsh.html#index-MinAniso
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))
include(joinpath(@__DIR__, "label_interface_patches.jl"))
using TOML
using SHA

function edge_metric(point, segment, minimum_size, maximum_size, growth, aspect, end_growth)
    a, b = collect(segment[1:3]), collect(segment[4:6])
    direction = b - a
    length_squared = dot(direction, direction)
    length_squared > 0 || error("Zero-length metric edge")
    tangent = direction / sqrt(length_squared)
    p = collect(point)
    distance = norm(p - a - clamp(dot(p - a, direction) / length_squared, 0., 1.) * direction)
    normal_size = min(maximum_size, minimum_size + growth * distance)
    end_size = minimum_size + end_growth * min(norm(p - a), norm(p - b))
    tangent_size = max(normal_size, min(maximum_size, aspect * normal_size, end_size))
    return Matrix{Float64}(I, 3, 3) / normal_size^2 +
           (1 / tangent_size^2 - 1 / normal_size^2) * tangent * tangent'
end

function tensor_field(segments, minimum_size, maximum_size, growth, aspect, end_growth)
    fields = Int[]
    number = 2000
    for segment in segments
        a, b = segment[1:3], segment[4:6]
        v = b .- a
        length_squared = sum(x^2 for x in v)
        tangent = v ./ sqrt(length_squared)
        delta = ["($(axis)-($(a[d])))" for (d,axis) in enumerate(("x","y","z"))]
        projection = "Min(1,Max(0,(" * join(["$(delta[d])*($(v[d]))" for d in 1:3], "+") * ")/($length_squared)))"
        distance = "Sqrt(" * join(["($(delta[d])-($projection)*($(v[d])))^2" for d in 1:3], "+") * ")"
        enda = join(["$(delta[d])^2" for d in 1:3], "+")
        endb = join(["($(axis)-($(b[d])))^2" for (d,axis) in enumerate(("x","y","z"))], "+")
        normal_size = "Min($maximum_size,$minimum_size+$growth*($distance))"
        tangent_size = "Max(($normal_size),Min($maximum_size,Min($aspect*($normal_size),$minimum_size+$end_growth*Sqrt(Min($enda,$endb)))))"
        # Inline expressions: nested MathEvalAniso -> MathEval references can
        # deadlock in Gmsh's expression-evaluation critical section.
        metric = number
        gmsh.model.mesh.field.add("MathEvalAniso", metric)
        for (name,i,j) in (("M11",1,1),("M12",1,2),("M13",1,3),("M22",2,2),("M23",2,3),("M33",3,3))
            if tangent[i]*tangent[j] == 0
                expression = i==j ? "1/(($normal_size)^2)" : "0"
            elseif i==j && abs(tangent[i])==1
                expression = "1/(($tangent_size)^2)"
            else
                diagonal = i==j ? "1/(($normal_size)^2)" : "0"
                expression = "$diagonal+(1/(($tangent_size)^2)-1/(($normal_size)^2))*($(tangent[i]*tangent[j]))"
            end
            gmsh.model.mesh.field.setString(metric,name,expression)
        end
        push!(fields,metric)
        number += 1
    end
    gmsh.model.mesh.field.add("MinAniso",number)
    gmsh.model.mesh.field.setNumbers(number,"FieldsList",Float64.(fields))
    gmsh.model.mesh.field.setAsBackgroundMesh(number)
    return fields
end

function main()
    length(ARGS)==7 || error("signature_dir thin|fabricated output.msh minimum_size maximum_size aspect reference-measures.csv")
    root=abspath(ARGS[1]);kind=ARGS[2];output=abspath(ARGS[3])
    minimum_size,maximum_size,aspect=parse.(Float64,ARGS[4:6])
    kind in ("thin","fabricated") || error("Invalid kind")
    all(isfinite,(minimum_size,maximum_size,aspect)) &&
        0<minimum_size<maximum_size && 1<=aspect<=8 || error("Invalid metric scout settings")
    isfile(output) && error("Refuse to overwrite a metric scout")
    process=TOML.parsefile(joinpath(root,"process.toml"))
    process["Units"]=="um" && process["SidewallAngle"]==90 &&
        process["TopRounding"]==process["TrenchRounding"]==0 ||
        error("This tensor scout requires sharp vertical planar geometry")
    growth=1.0;end_growth=0.25
    reference,_=readdlm(ARGS[7],',',header=true)
    expected=Dict((Int(reference[i,1]),Int(reference[i,2]))=>Float64(reference[i,3]) for i in axes(reference,1))
    # Same five-argument CAD control contract as the scalar experiment; this
    # existing tensor scout does not apply any scalar trace sizing.
    function controls(curves,groups,matching_surfaces,lower,upper)
        isempty(matching_surfaces) && error("No matching surfaces")
        occursin("Mmg",gmsh.option.getString("General.BuildOptions")) || error("Gmsh was built without MMG")
        actual=Dict{Tuple{Int,Int},Float64}()
        for dim in (2,3),(_,attribute) in gmsh.model.getPhysicalGroups(dim)
            actual[(dim,attribute)]=sum(gmsh.model.occ.getMass(dim,entity) for entity in gmsh.model.getEntitiesForPhysicalGroup(dim,attribute))
        end
        Set(keys(actual))==Set(keys(expected)) || error("Physical families changed")
        worst=maximum(abs(actual[k]/v-1) for (k,v) in expected)
        worst<=1e-8 || error("CAD geometry changed by $worst")
        segments=NTuple{6,Float64}[]
        for curve in curves
            lo,hi=gmsh.model.getParametrizationBounds(1,curve)
            a=gmsh.model.getValue(1,curve,[lo[1]]);b=gmsh.model.getValue(1,curve,[hi[1]])
            norm(b-a)>0 || error("Nonlinear closed metric edge")
            for t in (.2,.4,.6,.8)
                p=gmsh.model.getValue(1,curve,[lo[1]+t*(hi[1]-lo[1])])
                norm(cross(p-a,b-a))/norm(b-a)<1e-10 || error("Nonlinear curve needs an explicitly curved metric")
            end
            push!(segments,(a...,b...))
        end
        tensor_field(segments,minimum_size,maximum_size,growth,aspect,end_growth)
        for (name,value) in (("General.NumThreads",1),("General.Verbosity",4),("Mesh.MaxNumThreads1D",1),
                             ("Mesh.MaxNumThreads2D",1),("Mesh.MaxNumThreads3D",1),
                             ("Mesh.Algorithm",7),("Mesh.Algorithm3D",7),
                             ("Mesh.MeshSizeMin",minimum_size),("Mesh.MeshSizeMax",maximum_size))
            gmsh.option.setNumber(name,value)
        end
        println("Tensor fields ready: aspect=$aspect segments=$(length(segments))");flush(stdout)
        open(output*".metric.toml","w") do stream
            TOML.print(stream,Dict("ExperimentalOnly"=>true,"Metric"=>"Exact line-distance tensors with endpoint isotropy and MinAniso intersection",
                "AspectLimit"=>aspect,"NormalMinimum"=>minimum_size,"MaximumSize"=>maximum_size,
                "NormalGrowth"=>growth,"EndpointGrowth"=>end_growth,"Segments"=>length(segments),
                "CADRelativeMeasureDifference"=>worst,"Mesher3D"=>"MMG3D","SurfaceMesher"=>"BAMG",
                "GmshVersion"=>gmsh.option.getString("General.Version"),"SolverAccuracyQualified"=>false))
        end
        if get(ENV,"METRIC_REGIONWISE","0")=="1"
            # Diagnostic workaround for Gmsh builds whose MMG wrapper processes
            # an empty second region after adapting a multi-volume component.
            gmsh.model.mesh.generate(2)
            volumes=gmsh.model.getEntities(3)
            gmsh.option.setNumber("Mesh.MeshOnlyVisible",1)
            for volume in volumes
                gmsh.model.setVisibility(volumes,0,false)
                gmsh.model.setVisibility([volume],1,false)
                gmsh.model.mesh.generate(3)
            end
            gmsh.model.setVisibility(volumes,1,false)
            gmsh.option.setNumber("Mesh.MeshOnlyVisible",0)
        end
    end
    function label(edges,loops,radius)
        label_interface_patches(edges,loops,radius,output*".interface-partition.csv";
            fabricated=kind=="fabricated",metal_thickness=Float64(process["MetalThickness"]),
            overetch=Float64(process["Overetch"]))
    end
    try
        generate_spatial_coupon(signature=joinpath(root,"mesh-signature.csv"),
            mask=joinpath(root,"plan-view-mask.csv"),boundary=joinpath(root,"plan-view-boundary.csv"),
            fabricated=kind=="fabricated",radius=Float64(process["Radius"]),
            metal_thickness=Float64(process["MetalThickness"]),overetch=Float64(process["Overetch"]),
            sidewall_angle=90.,top_rounding=0.,trench_rounding=0.,lc_fine=minimum_size,
            lc_far=maximum_size,mesh_order=1,mesh_control=controls,mesh_postprocess=label,
            optimize_volume=false,max_elements=2_000_000,max_nodes=500_000,filename=output)
    finally
        gmsh.isInitialized()!=0 && gmsh.finalize()
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
