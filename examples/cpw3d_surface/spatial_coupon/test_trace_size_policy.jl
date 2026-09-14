# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
include(joinpath(@__DIR__, "mesh_graded_tet_experiment.jl"))

@testset "Trace supplied versus CAD constraints applied" begin
    @test !trace_constraints_applied("complete-trace.csv", "none")
    for mode in ("all", "sides", "levels")
        @test trace_constraints_applied("complete-trace.csv", mode)
        @test !trace_constraints_applied(nothing, mode)
    end
    @test !trace_constraints_applied(nothing, "none")
end

@testset "Explicit trace sizing policy and fail-closed defaults" begin
    @test trace_size_policy(Dict()).scope == :off
    @test trace_size_policy(Dict()).maximum_size == 0
    @test trace_size_policy(Dict("TET_TRACE_CONSTRAINT_MODE"=>"none")).mode == "none"
    @test_throws ErrorException trace_size_policy(Dict("TET_TRACE_SIZE"=>".01"))
    for scope in ("matching","matching-and-volume","legacy-global")
        env=Dict("TET_TRACE_SIZE"=>".01","TET_TRACE_SIZE_SCOPE"=>scope)
        @test_throws ErrorException trace_size_policy(env) # Missing explicit mode.
        env["TET_TRACE_CONSTRAINT_MODE"]="levels"
        @test trace_size_policy(env).scope == Symbol(replace(scope,"-"=>"_"))
        @test trace_size_policy(env).surface_growth == (scope=="legacy-global" ? .5 : 4.)
        env["TET_TRACE_CONSTRAINT_MODE"]="none"
        @test_throws ErrorException trace_size_policy(env)
        env["TET_TRACE_SIZE"]="0"
        @test_throws ErrorException trace_size_policy(env)
    end
    for (key,value) in (("TET_TRACE_SIZE","-1"),("TET_TRACE_SIZE","NaN"),
                        ("TET_TRACE_RELATIVE_SIZE",".2"),("TET_TRACE_RELATIVE_SIZE","Inf"),
                        ("TET_TRACE_SIZE_SCOPE","global"),("TET_TRACE_CONSTRAINT_MODE","bad"),
                        ("TET_TRACE_SURFACE_GROWTH","0"),("TET_TRACE_SURFACE_GROWTH","Inf"))
        @test_throws ErrorException trace_size_policy(Dict(key=>value))
    end
    @test_throws ErrorException trace_size_policy(Dict("TET_TRACE_SIZE"=>".01",
        "TET_TRACE_SIZE_SCOPE"=>"legacy-global","TET_TRACE_CONSTRAINT_MODE"=>"all",
        "TET_TRACE_SURFACE_GROWTH"=>"4"))
end

@testset "Trace callback scopes, entity restrictions, dim1 and volume opt-in" begin
    empty!(GRADING_SEGMENTS);push!(GRADING_SEGMENTS,(10.,0.,0.,11.,0.,0.))
    empty!(GRADING_ARCS);empty!(GRADING_CONICS)
    empty!(GRADING_MATCHING_SURFACES);push!(GRADING_MATCHING_SURFACES,Cint(71))
    empty!(GRADING_TRACE_SEGMENTS);push!(GRADING_TRACE_SEGMENTS,(0.,0.,0.,1.,0.,0.))
    empty!(GRADING_TRACE_SEGMENT_SIZES)
    SLOT_SIZE_TREE[]=nothing;GRADING_CAP_SIZE[]=0.;GRADING_TANGENT[]=0.
    GRADING_RIBBON_ASPECT[]=1.
    GRADING_FINE[]=.002;GRADING_FAR[]=.5;GRADING_GROWTH[]=1.
    GRADING_TRACE_SIZE[]=.01;GRADING_TRACE_SURFACE_GROWTH[]=4.
    GRADING_VOLUME_PROFILE[]=nothing
    query(dim,tag,z=0.)=exact_grading_callback(Cint(dim),Cint(tag),.5,0.,z,1.,C_NULL)
    for weighted in (false,true)
        empty!(GRADING_TRACE_SEGMENT_SIZES)
        weighted && push!(GRADING_TRACE_SEGMENT_SIZES,.005)
        h=weighted ? .005 : .01
        GRADING_TRACE_SIZE_SCOPE[]=:off # Even stale positive size/segments cannot leak.
        for dim in -1:3, tag in (71,4001)
            @test query(dim,tag)==.5
        end
        GRADING_TRACE_SIZE_SCOPE[]=:matching
        @test query(2,71)≈h
        @test query(2,71,.04)≈h+.16
        for dim in (-1,0,1,3), tag in (71,4001)
            @test query(dim,tag)==.5
        end
        @test query(2,4001)==.5 # Physical interface at identical coordinates.
        @test query(2,72)==.5 # Any other surface, not an entity whitelist wildcard.
        profile=(near_growth=.5,far_growth=2.,transition_distance=.03,
                 minimum_size=.002,maximum_size=.75,trace_near_growth=.5,
                 trace_far_growth=2.,trace_transition_distance=.03)
        GRADING_VOLUME_PROFILE[]=profile
        @test query(3,71,.1)==.75 # Matching-only does not change volume profile.
        @test query(2,71,.04)≈h+.16
        GRADING_TRACE_SIZE_SCOPE[]=:matching_and_volume
        @test query(3,71,.1)≈h+.015+.14
        @test query(3,4001,.02)≈h+.01
        @test query(2,4001)==.5
        @test query(1,71)==.5
        GRADING_TRACE_SIZE_SCOPE[]=:legacy_global
        for dim in -1:3, tag in (71,4001)
            @test query(dim,tag,.1)≈h+.05 # Replay ignores revised growth/profile.
        end
        GRADING_VOLUME_PROFILE[]=nothing
    end
    GRADING_TRACE_SIZE_SCOPE[]=:off
    GRADING_TRACE_SIZE[]=0.;empty!(GRADING_TRACE_SEGMENTS);empty!(GRADING_TRACE_SEGMENT_SIZES)
    empty!(GRADING_MATCHING_SURFACES)
end

@testset "Skinny triangles and explicit legacy minimum-altitude replay" begin
    # Short edge altitude is large; the two long edges still get the tiny scale.
    triangles=Dict(1=>[(0.,0.,0.),(10.,0.,0.),(0.,.00001,0.)])
    before=deepcopy(triangles)
    segments,sizes=trace_segment_sizes(triangles,.01,.2)
    short=findfirst(s->s[1:3]==(0.,0.,0.) && s[4:6]==(0.,.00001,0.),segments)
    @test sizes[short]==.01
    @test count(<(3e-6),sizes)==2
    legacy_segments,legacy_sizes=trace_segment_sizes(triangles,.01,.2;legacy=true,mode="levels")
    @test legacy_segments==segments
    @test all(h->h≈.2*.0001/sqrt(100+1e-10),legacy_sizes)
    @test triangles==before
end

# Complete synthetic box, with two side strips and unchanged nodal values/IDs.
function trace_policy_box()
    nodes=[(x,y,z) for z in (0.,.01,1.) for (x,y) in ((0.,0.),(1.,0.),(1.,1.),(0.,1.))]
    connectivity=[(1,3,2),(1,4,3),(9,10,11),(9,11,12)]
    for offset in (0,4), i in 1:4
        a=offset+i;b=offset+mod1(i+1,4)
        append!(connectivity,[(a,b,b+4),(a,b+4,a+4)])
    end
    return nodes,connectivity
end

@testset "All modes preserve full trace nodes, triangles, DOFs, values and bytes" begin
    nodes,connectivity=trace_policy_box()
    triangles=Dict(i=>[nodes[n] for n in conn] for (i,conn) in enumerate(connectivity))
    original=deepcopy(triangles)
    @test length(triangles)==20
    all_segments,_=trace_segment_sizes(triangles,.01,0.;mode="all")
    side_segments,_=trace_segment_sizes(triangles,.01,0.;mode="sides")
    level_segments,level_sizes=trace_segment_sizes(triangles,.01,.2;mode="levels")
    @test length(all_segments)==30
    @test length(side_segments)==28
    @test length(level_segments)==4
    @test all(s->s[3]==s[6]==.01,level_segments)
    @test all(==(.01),level_sizes) # Levels do not invent a function scale.
    @test triangles==original
    mktempdir() do root
        gmsh.initialize()
        try
            gmsh.option.setNumber("General.Verbosity",0)
            for basis in eachindex(nodes)
                path=joinpath(root,"basis-$basis.csv")
                open(path,"w") do stream
                    println(stream,"x,y,z,V,triangle,node,dof")
                    for (i,conn) in enumerate(connectivity), n in conn
                        println(stream,join((nodes[n]...,Float64(n==basis),i,n,n),','))
                    end
                end
                bytes=read(path);digest=sha256(bytes)
                data,header=readdlm(path,',',header=true)
                @test length(unique(data[:,6]))==length(nodes)==12
                for mode in ("none","levels","sides","all")
                    gmsh.model.add("$basis-$mode")
                    before=gmsh.model.occ.getEntities()
                    lines=matching_trace_lines(gmsh.model.occ,path,(0.,0.,0.),(1.,1.,1.),1e-9;mode=mode)
                    @test length(lines)==Dict("none"=>0,"levels"=>4,"sides"=>28,"all"=>30)[mode]
                    mode=="none" && @test gmsh.model.occ.getEntities()==before
                    readback=read_matching_trace_triangles(path)
                    @test readback==original
                    mode=="none" || trace_segment_sizes(readback,.01,.2;mode=mode)
                    @test sha256(read(path))==digest
                    @test read(path)==bytes
                    after,after_header=readdlm(path,',',header=true)
                    @test after==data && after_header==header
                    gmsh.model.remove()
                end
            end
        finally
            gmsh.finalize()
        end
    end
end

@testset "Malformed traces fail before every mode, even levels/none" begin
    valid=[(0.,0.,0.),(1.,0.,0.),(0.,1.,.01)]
    for malformed in (Dict{Int,Vector{NTuple{3,Float64}}}(),Dict(0=>valid),
                      Dict(1=>valid[1:2]),Dict(1=>vcat(valid,valid[1:1])),
                      Dict(1=>[valid[1],valid[1],valid[3]]),
                      Dict(1=>[(0.,0.,0.),(1.,0.,0.),(2.,0.,0.)]),
                      Dict(1=>[(NaN,0.,0.),valid[2],valid[3]]),
                      Dict(1=>[(0.,0.),(1.,0.),(0.,1.)]))
        for mode in ("all","sides","levels"), legacy in (false,true)
            @test_throws ErrorException trace_segment_sizes(malformed,.01,.2;mode=mode,legacy=legacy)
        end
    end
    mktempdir() do root
        path=joinpath(root,"trace.csv")
        for text in ("x,y,z,V,triangle\n", "x,y,z,V\n0,0,0,1\n",
                     "x,x,z,triangle\n0,0,0,1\n",
                     "x,y,z,triangle\n0,0,0,1\n1,0,0,1\n",
                     "x,y,z,triangle\n0,0,0,1\n1,0,0,1\n2,0,0,1\n",
                     "x,y,z,triangle\n0,0,0,1.5\n1,0,0,1.5\n0,1,1,1.5\n",
                     "x,y,z,triangle\n0,0,0,-1\n1,0,0,-1\n0,1,1,-1\n",
                     "x,y,z,triangle\nNaN,0,0,1\n1,0,0,1\n0,1,1,1\n")
            write(path,text)
            @test_throws Exception read_matching_trace_triangles(path)
            for mode in ("all","sides","levels","none")
                # `nothing` proves no CAD API is reached for malformed input.
                @test_throws Exception matching_trace_lines(nothing,path,(0.,0.,0.),(1.,1.,1.),1e-9;mode=mode)
            end
            @test read(path,String)==text
        end
    end
end

@testset "Physical-only CAD diagnostic exercises the five-argument scalar control" begin
    mktempdir() do root
        for (source,destination) in (("one-edge-masked-signature.csv","mesh-signature.csv"),
                                     ("one-edge-masked-mask.csv","plan-view-mask.csv"),
                                     ("one-edge-masked-boundary.csv","plan-view-boundary.csv"))
            cp(joinpath(@__DIR__,"testdata",source),joinpath(root,destination))
        end
        trace=joinpath(root,"unchanged-trace.csv")
        write(trace,"x,y,z,V,triangle\n0,0,0,1,1\n1,0,0,0,1\n0,1,1,0,1\n")
        source=read(trace)
        original_args=copy(ARGS)
        empty!(ARGS)
        append!(ARGS,[root,"thin",joinpath(root,"diagnostic.msh"),"1",".02",".5",
                      "--geometry-only","--matching-trace",trace])
        try
            withenv("TET_TRACE_SIZE"=>"0","TET_TRACE_SIZE_SCOPE"=>"off",
                    "TET_TRACE_RELATIVE_SIZE"=>"0","TET_TRACE_CONSTRAINT_MODE"=>"none",
                    "TET_GEOMETRY_ORDER"=>"1","TET_VERBOSITY"=>"0") do
                main() # Constructs only tiny fixture CAD; never generates a mesh.
            end
            @test !isfile(joinpath(root,"diagnostic.msh"))
            @test !isempty(GRADING_MATCHING_SURFACES)
            @test isempty(GRADING_TRACE_SEGMENTS)
            @test read(trace)==source
            report=read(joinpath(root,"diagnostic.msh.sizing.txt"),String)
            @test occursin("mode=none",report) && occursin("scope=off",report)
        finally
            empty!(ARGS);append!(ARGS,original_args)
            gmsh.isInitialized()!=0 && gmsh.finalize()
        end
    end
end
