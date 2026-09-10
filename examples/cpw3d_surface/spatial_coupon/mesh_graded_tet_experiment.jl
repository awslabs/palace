# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fresh all-tetrahedral geometry/meshing experiment. No prism splitting/extrusion.
# Args: signature_dir kind output.msh growth [minimum_size=0.002] [maximum_size=0.5]
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))
include(joinpath(@__DIR__, "label_interface_patches.jl"))
include(joinpath(@__DIR__, "graded_curve_distance.jl"))
include(joinpath(@__DIR__, "graded_size_points.jl"))
using TOML
using SHA
const SIZE_CALLBACK_ROOTS = Any[]
const GRADING_SEGMENTS = NTuple{6,Float64}[]
const GRADING_ARCS = GradingArc[]
const GRADING_CONICS = GradingConic[]
const GRADING_FINE = Ref(0.002)
const GRADING_FAR = Ref(0.5)
const GRADING_GROWTH = Ref(1.0)
const GRADING_VOLUME_PROFILE = Ref{Union{Nothing,NamedTuple}}(nothing)
const GRADING_QUERY_COUNTS = zeros(Int,5) # dimensions -1 through 3

function volume_grading_size(distance,fine,profile)
    minimum_size=hasproperty(profile,:minimum_size) ? profile.minimum_size : fine
    return min(profile.maximum_size,
               minimum_size+profile.near_growth*min(distance,profile.transition_distance)+
               profile.far_growth*max(0.,distance-profile.transition_distance))
end
const GRADING_TANGENT = Ref(0.0)
const GRADING_HORIZONTAL_CURVES = Dict{Cint,NTuple{6,Float64}}()
const SLOT_REFINEMENT_POINTS = NTuple{4,Float64}[]
const SLOT_SIZE_TREE = Ref{Union{Nothing,SizePointTree}}(nothing)
const GRADING_CAP_SIZE = Ref(0.0)
const GRADING_CAP_Z = Ref((0.0,0.0))

# A named callback avoids closure trampolines, unsupported by Julia on aarch64.
function exact_grading_callback(dim,tag,x,y,z,lc,data)::Cdouble
    -1<=dim<=3 && (GRADING_QUERY_COUNTS[Int(dim)+2]+=1)
    tangent_override=nothing
    if dim == 1 && GRADING_TANGENT[] > 0 && haskey(GRADING_HORIZONTAL_CURVES, tag)
        s = GRADING_HORIZONTAL_CURVES[tag]
        endpoint_distance = min(sqrt((x-s[1])^2+(y-s[2])^2+(z-s[3])^2),
                                sqrt((x-s[4])^2+(y-s[5])^2+(z-s[6])^2))
        tangent_override=min(GRADING_TANGENT[], GRADING_FINE[] + 0.5endpoint_distance)
    end
    best=Inf
    @inbounds for s in GRADING_SEGMENTS
        dx,dy,dz=s[4]-s[1],s[5]-s[2],s[6]-s[3]
        t=clamp(((x-s[1])*dx+(y-s[2])*dy+(z-s[3])*dz)/(dx*dx+dy*dy+dz*dz),0.,1.)
        d2=(x-s[1]-t*dx)^2+(y-s[2]-t*dy)^2+(z-s[3]-t*dz)^2
        best=min(best,d2)
    end
    @inbounds for arc in GRADING_ARCS
        best=min(best,grading_arc_distance2((x,y,z),arc))
    end
    @inbounds for conic in GRADING_CONICS
        best=min(best,grading_conic_distance2((x,y,z),conic))
    end
    distance=sqrt(best)
    size=dim==3 && GRADING_VOLUME_PROFILE[]!==nothing ?
         volume_grading_size(distance,GRADING_FINE[],GRADING_VOLUME_PROFILE[]) :
         min(GRADING_FAR[],GRADING_FINE[]+GRADING_GROWTH[]*distance)
    tangent_override === nothing || (size=tangent_override)
    if GRADING_CAP_SIZE[] > 0
        distance=min(abs(z-GRADING_CAP_Z[][1]),abs(z-GRADING_CAP_Z[][2]))
        size=min(size,GRADING_CAP_SIZE[]+0.5distance)
    end
    SLOT_SIZE_TREE[]===nothing || (size=query_size_point(SLOT_SIZE_TREE[],(x,y,z),size))
    return size
end

function main()
    get(ENV,"TET_ANISOTROPIC_SURFACE","0") in ("0","1") ||
        error("TET_ANISOTROPIC_SURFACE must be 0 or 1")
    get(ENV,"TET_TRACE_CONSTRAINT_MODE","all") in ("all","levels","sides") ||
        error("Unsupported trace constraint mode")
    geometry_only = "--geometry-only" in ARGS
    element_slots = "--element-interface-slots" in ARGS
    args = filter(arg -> arg != "--geometry-only" && arg != "--element-interface-slots", ARGS)
    volume_study=nothing
    study_flag=findfirst(==("--volume-study"),args)
    if study_flag!==nothing
        study_flag<length(args) || error("Missing volume-study TOML file")
        volume_study=TOML.parsefile(args[study_flag+1])
        deleteat!(args,study_flag:study_flag+1)
        geometry_only && error("Volume study and geometry-only are mutually exclusive")
        element_slots && error("Frozen-boundary volume studies retain CAD family labels; do not relabel between variants")
    end
    GRADING_VOLUME_PROFILE[]=nothing
    process=Dict("Units"=>"um", "Radius"=>2.0, "MetalThickness"=>0.1,
                 "Overetch"=>0.05, "SidewallAngle"=>90.0, "TopRounding"=>0.0,
                 "TrenchRounding"=>0.0)
    process_flag=findfirst(==("--process"),args)
    if process_flag!==nothing
        process_flag<length(args) || error("Missing process TOML path")
        supplied=TOML.parsefile(args[process_flag+1])
        haskey(supplied,"Units") || error("Process file must explicitly declare Units = \"um\"")
        isempty(setdiff(keys(supplied),keys(process))) || error("Unknown process settings")
        merge!(process,supplied)
        deleteat!(args,process_flag:process_flag+1)
    end
    process["Units"]=="um" || error("Process inputs currently require explicit um units")
    all(v isa Real && !(v isa Bool) && isfinite(v) for (k,v) in process if k!="Units") ||
        error("Process dimensions must be finite numbers")
    radius=Float64(process["Radius"])
    thickness=Float64(process["MetalThickness"])
    etch=Float64(process["Overetch"])
    expected_interfaces=nothing
    interfaces_flag=findfirst(==("--expected-interfaces"),args)
    if interfaces_flag!==nothing
        interfaces_flag<length(args) || error("Missing expected-interface CSV path")
        data,_=readdlm(args[interfaces_flag+1],',',header=true)
        expected_interfaces=Set(Int.(vec(data)))
        isempty(expected_interfaces) && error("Empty expected interface set")
        deleteat!(args,interfaces_flag:interfaces_flag+1)
    end
    matching_trace = nothing
    trace_flag = findfirst(==("--matching-trace"), args)
    if trace_flag !== nothing
        trace_flag < length(args) || error("Missing matching-trace path")
        matching_trace = abspath(args[trace_flag+1])
        deleteat!(args,trace_flag:(trace_flag+1))
    end
    slot_refine_path = get(ENV,"TET_SLOT_REFINEMENT_POINTS","")
    empty!(SLOT_REFINEMENT_POINTS)
    if !isempty(slot_refine_path)
        data,_=readdlm(slot_refine_path,',',header=true)
        for i in axes(data,1)
            p=ntuple(d->Float64(data[i,d]),4)
            all(isfinite,p) && p[4]>0 || error("Invalid slot refinement point")
            push!(SLOT_REFINEMENT_POINTS,p)
        end
    end
    SLOT_SIZE_TREE[]=SizePointTree(SLOT_REFINEMENT_POINTS)
    reference_measures = nothing
    reference_flag = findfirst(==("--reference-measures"), args)
    if reference_flag !== nothing
        reference_flag < length(args) || error("Missing reference-measures CSV path")
        reference_measures = abspath(args[reference_flag + 1])
        deleteat!(args, reference_flag:(reference_flag + 1))
    end
    length(args) in (4,5,6) || error("signature_dir thin|fabricated output.msh growth [hmin] [hmax] [--geometry-only] [--reference-measures CSV]")
    geometry_only || reference_measures !== nothing ||
        error("A full meshing experiment requires --reference-measures; audit geometry first")
    geometry_order=parse(Int,get(ENV,"TET_GEOMETRY_ORDER","2"))
    geometry_order in (1,2) || error("Unsupported geometry order")
    volume_study===nothing || geometry_order==1 || error("Frozen-boundary study currently requires exact planar linear geometry")
    root=abspath(args[1]);kind=args[2];output=abspath(args[3]);growth=parse(Float64,args[4])
    kind in ("thin","fabricated") || error("invalid kind")
    fine=length(args)>=5 ? parse(Float64,args[5]) : 0.002
    slot_minimum=parse(Float64,get(ENV,"TET_SLOT_MINIMUM_SIZE",string(fine)))
    isfinite(slot_minimum) && 0<slot_minimum<=fine || error("Invalid slot-refinement minimum size")
    far=length(args)>=6 ? parse(Float64,args[6]) : 0.5
    all(isfinite, (growth,fine,far)) && growth>0 && fine>0 && far>fine || error("invalid size parameters")
    isfile(output) && error("Refuse to overwrite an existing mesh")
    mkpath(dirname(output))
    function controls(curves,boundary_groups,lower,upper)
        # Preserve callback roots for the complete native meshing call; the generated
        # Julia Gmsh binding does not retain the closure CFunction after returning.
        segments=NTuple{6,Float64}[]
        empty!(GRADING_ARCS)
        empty!(GRADING_CONICS)
        GRADING_CAP_SIZE[]=parse(Float64,get(ENV,"TET_CAP_SIZE","0"))
        isfinite(GRADING_CAP_SIZE[]) && GRADING_CAP_SIZE[]>=0 || error("Invalid cap size")
        GRADING_CAP_Z[]=(lower[3],upper[3])
        empty!(GRADING_HORIZONTAL_CURVES)
        GRADING_TANGENT[] = parse(Float64, get(ENV, "TET_EDGE_TANGENT_SIZE", "0"))
        isfinite(GRADING_TANGENT[]) || error("Non-finite tangential size")
        GRADING_TANGENT[] == 0 || GRADING_TANGENT[] >= fine || error("Tangential size is smaller than minimum size")
        for curve in curves
            bounds=gmsh.model.getParametrizationBounds(1,curve)
            lo,hi=bounds[1][1],bounds[2][1]
            points=[gmsh.model.getValue(1,curve,[lo+t*(hi-lo)]) for t in (0.,0.25,0.5,0.75,1.)]
            a,b=points[1],points[end]
            curve_kind=gmsh.model.getType(1,curve)
            scale=max(1.,maximum(norm(p-a) for p in points))
            chord=norm(b-a)
            # OCC lofts often encode exact straight edges as Bezier/BSpline curves.
            # Check the spatial curve, not only the CAD representation's type name.
            straight=chord>0 && abs(gmsh.model.occ.getMass(1,curve)-chord)<=1e-11scale &&
                all(grading_segment_distance2(Tuple(p),Tuple(vcat(a,b)))<=1e-20scale^2 for p in points)
            if curve_kind=="Line" || straight
                chord>0 || error("zero-length feature curve")
                push!(segments,Tuple(vcat(a,b)))
            elseif curve_kind=="Circle"
                arc=grading_arc(a,points[2],points[3],b;tolerance=1e-8scale)
                all(grading_arc_distance2(Tuple(p),arc)<=1e-16scale^2 for p in points) ||
                    error("Inconsistent circular feature curve")
                push!(GRADING_ARCS,arc)
                geometry_order>=2 || error("Circular geometry requires geometric order >= 2")
            elseif curve_kind in ("Ellipse","TrimmedCurve")
                # OCC fillets can wrap a conic in a generic TrimmedCurve. Validate
                # its analytic position and derivative model before using it.
                parameters=collect(range(lo,hi;length=17))
                sampled=[Tuple(gmsh.model.getValue(1,curve,[t])) for t in parameters]
                conic,center,u,v=grading_conic(parameters,sampled,fine)
                for t in (lo+(hi-lo)*0.1234567,lo+(hi-lo)*0.6180339)
                    p=gmsh.model.getValue(1,curve,[t])
                    derivative=gmsh.model.getDerivative(1,curve,[t])
                    norm(center+u*cos(t)+v*sin(t)-p)<=1e-9scale &&
                        norm(-u*sin(t)+v*cos(t)-derivative)<=1e-8scale ||
                        error("Unsupported nonlinear trimmed curve")
                end
                push!(GRADING_CONICS,conic)
                geometry_order>=2 || error("Conic geometry requires geometric order >= 2")
            else
                error("Unsupported nonlinear grading curve $curve_kind: supported geometry is straight segments and validated conic arcs")
            end
            if maximum(p[3] for p in points)-minimum(p[3] for p in points)<1e-10
                GRADING_HORIZONTAL_CURVES[Cint(curve)] = Tuple(vcat(a,b))
            end
        end
        isempty(segments) && isempty(GRADING_ARCS) && isempty(GRADING_CONICS) && error("No feature curves")
        empty!(GRADING_SEGMENTS);append!(GRADING_SEGMENTS,segments)
        GRADING_FINE[]=fine;GRADING_FAR[]=far;GRADING_GROWTH[]=growth
        # Unit checks include segment interiors: avoid sampled-distance aliasing on
        # long edges with a 2 nm target size.
        for s in segments
            size=exact_grading_callback(2,0,(s[1]+s[4])/2,(s[2]+s[5])/2,(s[3]+s[6])/2,far,C_NULL)
            @assert 0<size<=fine+1e-10 # Optional partition hints may request a smaller size.
        end
        callback=@cfunction(exact_grading_callback,Cdouble,(Cint,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cvoid}))
        push!(SIZE_CALLBACK_ROOTS,callback)
        ierr=Ref{Cint}()
        ccall((:gmshModelMeshSetSizeCallback,gmsh.lib),Cvoid,
              (Ptr{Cvoid},Ptr{Cvoid},Ptr{Cint}),callback,C_NULL,ierr)
        ierr[]==0 || error(gmsh.logger.getLastError())
        gmsh.model.mesh.field.add("MathEval",999)
        gmsh.model.mesh.field.setString(999,"F",string(far))
        gmsh.model.mesh.field.setAsBackgroundMesh(999)
        gmsh.option.setNumber("Mesh.MeshSizeMin",min(fine,slot_minimum))
        gmsh.option.setNumber("General.Verbosity", parse(Int, get(ENV,"TET_VERBOSITY","4")))
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.option.setNumber("Mesh.MaxNumThreads1D",1)
        gmsh.option.setNumber("Mesh.MaxNumThreads2D",1)
        gmsh.option.setNumber("Mesh.MaxNumThreads3D",1)
        surface_algorithm=parse(Int,get(ENV,"TET_SURFACE_ALGORITHM","6"))
        surface_algorithm in (5,6,7) || error("Unsupported surface algorithm")
        gmsh.option.setNumber("Mesh.Algorithm",surface_algorithm)
        algorithm3d=parse(Int,get(ENV,"TET_ALGORITHM3D","10"))
        algorithm3d in (1,10) || error("Unsupported experimental 3D algorithm")
        gmsh.option.setNumber("Mesh.Algorithm3D",algorithm3d)
        hxt_quality=parse(Float64,get(ENV,"TET_HXT_QUALITY","0.3"))
        isfinite(hxt_quality) && 0<hxt_quality<1 || error("Invalid HXT quality target")
        gmsh.option.setNumber("Mesh.OptimizeThreshold",hxt_quality)
        gmsh.option.setNumber("Mesh.RecombineAll",0)
        layer_width = parse(Float64, get(ENV, "TET_SURFACE_LAYER_WIDTH", "0"))
        isfinite(layer_width) && layer_width >= 0 || error("Invalid surface-layer width")
        if layer_width > 0
            gmsh.model.mesh.field.add("BoundaryLayer", 998)
            gmsh.model.mesh.field.setNumbers(998,"CurvesList",Float64.(curves))
            gmsh.model.mesh.field.setNumber(998,"Size",fine)
            gmsh.model.mesh.field.setNumber(998,"SizeFar",far)
            gmsh.model.mesh.field.setNumber(998,"Ratio",1.4)
            gmsh.model.mesh.field.setNumber(998,"Thickness",layer_width)
            gmsh.model.mesh.field.setNumber(998,"Quads",0)
            gmsh.model.mesh.field.setAsBoundaryLayer(998)
        end
        measures = Dict{Tuple{Int,Int},Float64}()
        open(output*".cad-measures.csv","w") do f
            println(f,"dimension,attribute,measure")
            for dim in (2,3),(_,attribute) in gmsh.model.getPhysicalGroups(dim)
                measure=sum(gmsh.model.occ.getMass(dim,entity) for entity in gmsh.model.getEntitiesForPhysicalGroup(dim,attribute))
                measures[(dim,attribute)] = measure
                println(f,"$dim,$attribute,$measure")
            end
        end
        if reference_measures !== nothing
            data, _ = readdlm(reference_measures, ',', header=true)
            expected = Dict((Int(data[i,1]),Int(data[i,2]))=>Float64(data[i,3]) for i in axes(data,1))
            if element_slots
                # CAD must already have the right physical geometry. Bookkeeping slot
                # areas are checked/reported separately after element-wise labeling.
                function physical_key(key)
                    dim,attribute=key
                    dim==3 && return key
                    family=div(attribute,1000)
                    family in (4,5,6) && return (dim,1000family+mod(attribute,100))
                    family==3 && return (dim,attribute>=3100 ? 3100 : 3000)
                    return key
                end
                function grouped(input)
                    grouped_measures=Dict{Tuple{Int,Int},Float64}()
                    for (key,value) in input
                        k=physical_key(key);grouped_measures[k]=get(grouped_measures,k,0.)+value
                    end
                    return grouped_measures
                end
                measures=grouped(measures);expected=grouped(expected)
            end
            Set(keys(expected)) == Set(keys(measures)) ||
                error("Geometry gate: physical attributes differ; missing=$(setdiff(Set(keys(expected)),Set(keys(measures)))) extra=$(setdiff(Set(keys(measures)),Set(keys(expected))))")
            errors = Dict(key=>abs(measures[key]-value)/max(abs(value),1e-300) for (key,value) in expected)
            worst = maximum(Base.values(errors))
            open(output*".geometry-gate.txt","w") do f
                println(f,"scope=$(element_slots ? "physical-family measures only; slot partition not qualified" : "per-attribute measures")")
                println(f,"maximum_relative_measure_difference=$worst tolerance=1e-5")
                for key in sort!(collect(keys(expected)))
                    println(f,"$key reference=$(expected[key]) actual=$(measures[key]) relative=$(errors[key])")
                end
            end
            worst <= 1e-5 || error("Geometry gate failed: maximum relative area/volume difference $worst")
            println("Geometry gate passed: maximum relative area/volume difference $worst")
        end
        open(output*".process.toml","w") do f
            TOML.print(f,process;sorted=true)
        end
        open(output*".sizing.txt","w") do f
            println(f,"h=min($far,$fine+$growth*exact_distance_to_supported_CAD_curves)")
            println(f,"physical_segments=$(length(segments)) circular_arcs=$(length(GRADING_ARCS)) conic_arcs=$(length(GRADING_CONICS)) algorithm3d=$algorithm3d surface_algorithm=$surface_algorithm threads=1")
            println(f,"matching_trace=$matching_trace mode=$(get(ENV,"TET_TRACE_CONSTRAINT_MODE","all"))")
            println(f,"cap_size=$(GRADING_CAP_SIZE[]) slot_refinement_points=$(length(SLOT_REFINEMENT_POINTS)) slot_minimum_size=$slot_minimum")
            println(f,"surface_boundary_layer_width=$layer_width hxt_quality_target=$hxt_quality")
            println(f,"edge_tangent_size=$(GRADING_TANGENT[]) (0 means fully isotropic sizing)")
            println(f,"conic_distance_model_bound=$(maximum((c.distance_error for c in GRADING_CONICS);init=0.))")
            println(f,"gmsh_version=$(gmsh.option.getString("General.Version"))")
            println(f,"lower=$lower upper=$upper")
        end
        println("Exact 3D grading: $(length(segments)) segments, $(length(GRADING_ARCS)) arcs, hmin=$fine hmax=$far growth=$growth")
        flush(stdout)
        if volume_study!==nothing
            isempty(GRADING_ARCS) && isempty(GRADING_CONICS) || error("Frozen-boundary study requires planar geometry")
            isempty(SLOT_REFINEMENT_POINTS) || error("Do not vary partition refinement in a volume study")
            get(ENV,"TET_ANISOTROPIC_SURFACE","0")=="0" || error("Do not vary surface meshing in a volume study")
            run_frozen_volume_study(volume_study,output,lower,upper,fine)
        end
        if get(ENV,"TET_ANISOTROPIC_SURFACE","0")=="1" && !geometry_only
            GRADING_TANGENT[]>0 || error("Anisotropic surface meshing needs a tangential size")
            # Mesh curves with the original tangential controls, then let BAMG
            # resolve the surface-normal direction independently. Restore the exact
            # scalar 3D field for HXT; no anisotropic 3D algorithm is assumed here.
            gmsh.model.mesh.generate(1)
            gmsh.model.mesh.removeSizeCallback()
            gmsh.model.mesh.field.add("AttractorAnisoCurve",997)
            gmsh.model.mesh.field.setNumbers(997,"CurvesList",Float64.(curves))
            gmsh.model.mesh.field.setNumber(997,"DistMin",0.)
            gmsh.model.mesh.field.setNumber(997,"DistMax",min(radius,far/growth))
            gmsh.model.mesh.field.setNumber(997,"SizeMinNormal",fine)
            gmsh.model.mesh.field.setNumber(997,"SizeMaxNormal",far)
            gmsh.model.mesh.field.setNumber(997,"SizeMinTangent",GRADING_TANGENT[])
            gmsh.model.mesh.field.setNumber(997,"SizeMaxTangent",far)
            gmsh.model.mesh.field.setNumber(997,"Sampling",25000)
            gmsh.model.mesh.field.setAsBackgroundMesh(997)
            gmsh.option.setNumber("Mesh.Algorithm",7)
            gmsh.model.mesh.generate(2)
            gmsh.write(output*".surface.msh")
            gmsh.model.mesh.field.setAsBackgroundMesh(999)
            ccall((:gmshModelMeshSetSizeCallback,gmsh.lib),Cvoid,
                  (Ptr{Cvoid},Ptr{Cvoid},Ptr{Cint}),callback,C_NULL,ierr)
            ierr[]==0 || error(gmsh.logger.getLastError())
            println("Anisotropic surface mesh complete; restoring exact 3D grading")
            flush(stdout)
        end
    end
    function relabel_and_check(edges,loops,radius)
        report=output*".interface-partition.csv"
        label_interface_patches(edges,loops,radius,report;minimum_size=slot_minimum,
                                fabricated=kind=="fabricated",metal_thickness=thickness,overetch=etch)
        if reference_measures !== nothing
            expected_data,_=readdlm(reference_measures,',',header=true)
            expected=expected_interfaces===nothing ?
                Set(Int(expected_data[i,2]) for i in axes(expected_data,1)
                    if Int(expected_data[i,1])==2 && Int(expected_data[i,2])!=1) : expected_interfaces
            actual_data,_=readdlm(report,',',header=true)
            actual=Set(Int(actual_data[i,1]) for i in axes(actual_data,1))
            actual==expected || error("Element-wise interface coverage differs from reference: missing=$(setdiff(expected,actual)), extra=$(setdiff(actual,expected))")
        end
    end
    postprocess = element_slots ? relabel_and_check : nothing
    try
        generate_spatial_coupon(signature=joinpath(root,"mesh-signature.csv"),
        mask=joinpath(root,"plan-view-mask.csv"),boundary=joinpath(root,"plan-view-boundary.csv"),
        fabricated=kind=="fabricated",radius=radius,metal_thickness=thickness,overetch=etch,
        sidewall_angle=Float64(process["SidewallAngle"]),
        top_rounding=Float64(process["TopRounding"]),
        trench_rounding=Float64(process["TrenchRounding"]),lc_fine=fine,
        lc_tangent=0.,lc_far=far,process_core_width=0.4radius,process_fine_width=0.,
        process_grading_power=1.,max_nodes=10_000_000,max_elements=5_000_000,
        mesh_order=geometry_order,mesh_control=controls,mesh_postprocess=postprocess,
        optimize_volume=get(ENV,"TET_ANISOTROPIC_SURFACE","0")!="1",
        geometry_only=geometry_only || volume_study!==nothing,matching_trace=matching_trace,
        matching_trace_mode=get(ENV,"TET_TRACE_CONSTRAINT_MODE","all"),filename=output)
    finally
        gmsh.isInitialized() != 0 && gmsh.finalize()
    end
    if !geometry_only && volume_study===nothing
        metadata_path = output * ".metadata.json"
        metadata = read(metadata_path, String)
        trace_mode=get(ENV,"TET_TRACE_CONSTRAINT_MODE","all")
        metadata = replace(metadata, "\"Version\": 1," =>
            "\"Version\": 1,\n  \"Experimental3DGrading\": true,\n  \"AnisotropicSurfaceMeshing\": $(get(ENV,"TET_ANISOTROPIC_SURFACE","0")=="1"),\n  \"SurfaceMeshingAlgorithm\": $(get(ENV,"TET_ANISOTROPIC_SURFACE","0")=="1" ? 7 : parse(Int,get(ENV,"TET_SURFACE_ALGORITHM","6"))),\n  \"MatchingTraceConstraints\": $(matching_trace !== nothing),\n  \"TraceConstraintMode\": \"$trace_mode\",\n  \"CapSize\": $(GRADING_CAP_SIZE[]),\n  \"SlotRefinementPoints\": $(length(SLOT_REFINEMENT_POINTS)),\n  \"SlotMinimumSize\": $slot_minimum,\n  \"ElementInterfaceSlots\": $element_slots,\n  \"InterfacePartitionQualified\": false,\n  \"ExpectedInterfaceAttributes\": $(expected_interfaces===nothing ? "null" : "["*join(sort!(collect(expected_interfaces)),",")*"]"),\n  \"Growth\": $growth,\n  \"HXTQualityTarget\": $(parse(Float64,get(ENV,"TET_HXT_QUALITY","0.3"))),\n  \"VolumeMeshingAlgorithm\": $(parse(Int,get(ENV,"TET_ALGORITHM3D","10"))),\n  \"EdgeCurveTangentialSize\": $(GRADING_TANGENT[]),\n  \"NormalLayerEnforced\": false,\n  \"NetgenOptimization\": $(get(ENV,"TET_ANISOTROPIC_SURFACE","0")!="1"),")
        write(metadata_path, metadata)
        if element_slots
            certificate=output*".interface-partition.csv.elements.csv"
            open(output*".partition-certificate.toml","w") do f
                TOML.print(f,Dict("Version"=>1,"Method"=>"Whole-element Lipschitz ownership",
                    "MeshSHA256"=>bytes2hex(sha256(read(output))),
                    "ElementCertificateSHA256"=>bytes2hex(sha256(read(certificate))),
                    "SignatureSHA256"=>bytes2hex(sha256(read(joinpath(root,"mesh-signature.csv")))),
                    "BoundarySHA256"=>bytes2hex(sha256(read(joinpath(root,"plan-view-boundary.csv")))),
                    "Radius"=>radius,"Fabricated"=>kind=="fabricated",
                    "MetalThickness"=>thickness,"Overetch"=>etch))
            end
        end
    end
end

include(joinpath(@__DIR__,"frozen_volume_study.jl"))

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
