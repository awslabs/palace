# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using Random
include(joinpath(@__DIR__, "mesh_graded_tet_experiment.jl"))

@testset "Exact mask membership" begin
    edge=(point=(0.,0.,0.),conductor=1,tangent=(1.,0.,0.),gap=(0.,-1.,0.),
          interval=(0.,1.),vertex_arm=false)
    facets=[(conductor=1,plane=0.,points=[(0.,0.),(10.,0.),(10.,10.),(0.,10.)])]
    @test point_in_metal(edge,(8.,8.,0.),0.1,1e-9,facets)
    @test !point_in_metal(edge,(11.,8.,0.),0.1,1e-9,facets)
    @test !point_in_metal(edge,(8.,8.,0.1),0.1,1e-9,facets)
end

@testset "Etch respects the union of conductors" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.model.add("etch_mask_union")
        loops=[(conductor=c,plane=0.,hole=false,
                points=[(x,0.),(x+1.,0.),(x+1.,1.),(x,1.)],
                classes=fill("Physical",4)) for (c,x) in [(1,0.),(2,1.1)]]
        trench=boundary_strips(gmsh.model.occ,loops,0.2,0.,-0.05,0.,1e-9)
        @test sum(gmsh.model.occ.getMass(3,t) for (_,t) in trench) ≈ 0.263 atol=1e-10
        retained=loft_mask(gmsh.model.occ,loops,0.,-0.05,0.,1e-9)
        intersection,_=gmsh.model.occ.intersect(trench,retained)
        @test sum((gmsh.model.occ.getMass(3,t) for (d,t) in intersection if d==3);init=0.)<1e-10
    finally
        gmsh.finalize()
    end
end

@testset "Fabrication holes respect material side and collapse" begin
    outer=(conductor=1,plane=0.,hole=false,points=[(-.4,-.4),(.4,-.4),(.4,.4),(-.4,.4)],classes=fill("Physical",4))
    hole=(conductor=1,plane=0.,hole=true,points=[(-.1,-.1),(.1,-.1),(.1,.1),(-.1,.1)],classes=fill("Physical",4))
    for points in (hole.points,reverse(hole.points))
        h=merge(hole,(points=points,))
        grown=offset_hole_points(h,0.02,1e-9)
        @test maximum(p[1] for p in grown) ≈ 0.12
        @test isempty(offset_hole_points(h,-0.3,1e-9))
        @test isempty(offset_hole_points(h,-0.1,1e-9))
    end
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.model.add("etched_hole")
        trenches=boundary_strips(gmsh.model.occ,[outer,hole],0.1,0.,-0.05,0.,1e-9)
        @test sum(gmsh.model.occ.getMass(3,t) for (_,t) in trenches) ≈ (1.4^2-(0.8^2-0.2^2))*0.05 atol=1e-10
        @test_throws ErrorException loft_mask(gmsh.model.occ,[outer,merge(hole,(conductor=2,))],0.,.1,0.,1e-9)
    finally
        gmsh.finalize()
    end
end

@testset "Planar mask holes and rigid transforms" begin
    for angle in (0.,pi/6)
        gmsh.initialize()
        try
            gmsh.option.setNumber("General.Verbosity",0)
            gmsh.model.add("mask_holes")
            transform(p)=(cos(angle)*p[1]-sin(angle)*p[2]+3.,
                          sin(angle)*p[1]+cos(angle)*p[2]-2.)
            outer=map(transform,[(0.,0.),(4.,0.),(4.,3.),(0.,3.)])
            hole=map(transform,[(1.,1.),(2.,1.),(2.,2.),(1.,2.)])
            loops=[(conductor=1,plane=0.,hole=false,points=outer,classes=fill("Physical",4)),
                   (conductor=1,plane=0.,hole=true,points=hole,classes=fill("Physical",4))]
            faces=planar_mask_surfaces(gmsh.model.occ,loops,1e-9)
            @test sum(gmsh.model.occ.getMass(2,t) for (_,t) in faces) ≈ 11.0 atol=1e-10
            orphan=(conductor=2,plane=0.,hole=true,points=hole,classes=fill("Physical",4))
            @test_throws ErrorException planar_mask_surfaces(gmsh.model.occ,[loops[1],orphan],1e-9)
        finally
            gmsh.finalize()
        end
    end
end

@testset "Element slots preserve a coplanar face" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.option.setNumber("General.NumThreads",1)
        gmsh.option.setNumber("Mesh.MeshSizeMin",0.1)
        gmsh.option.setNumber("Mesh.MeshSizeMax",0.1)
        gmsh.model.add("slot_partition")
        face=gmsh.model.occ.addRectangle(0.,0.,0.,1.,1.)
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(2,[face],4001)
        gmsh.model.mesh.generate(2)
        types,tags,nodes=gmsh.model.mesh.getElements(2)
        before=Dict(tag=>Tuple(conn[(i-1)*3+1:i*3]) for (et,conn) in zip(tags,nodes) for (i,tag) in enumerate(et))
        edges=[(point=(x,0.,0.),tangent=(0.,1.,0.),gap=(1.,0.,0.),
                interval=(0.,1.),vertex_arm=false,conductor=1,slot=s)
               for (x,s) in [(0.,0),(1.,1)]]
        mktempdir() do root
            report=joinpath(root,"partition.csv")
            label_interface_patches(edges,NamedTuple[],1.0,report)
            groups=sort([Int(a) for (_,a) in gmsh.model.getPhysicalGroups(2)])
            @test groups==[4001,4101]
            data,_=readdlm(report,',',header=true)
            @test sum(data[:,3]) ≈ 1.0 atol=1e-12
            @test all(0 .<= data[:,5] .<= 1)
            after=Dict{UInt64,Tuple}()
            for (_,a) in gmsh.model.getPhysicalGroups(2), entity in gmsh.model.getEntitiesForPhysicalGroup(2,a)
                _,ts,cs=gmsh.model.mesh.getElements(2,entity)
                for (et,conn) in zip(ts,cs), (i,tag) in enumerate(et)
                    after[tag]=Tuple(conn[(i-1)*3+1:i*3])
                end
            end
            @test before==after
        end
    finally
        gmsh.finalize()
    end
end

@testset "Trace level constraints and malformed input" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.model.add("trace_levels")
        mktempdir() do root
            path=joinpath(root,"trace.csv")
            write(path,"x,y,z,V,triangle\n0,0,0.4,0,1\n1,0,0.4,0,1\n0,1,0.6,0,1\n")
            lines=matching_trace_lines(gmsh.model.occ,path,(0.,0.,0.),(1.,1.,1.),1e-9;mode="levels")
            @test length(lines)==8
            write(path,"x,y,z,V,triangle\n0,0,0.4,0,1\n1,0,0.4,0,1\n")
            @test_throws ErrorException matching_trace_lines(gmsh.model.occ,path,(0.,0.,0.),(1.,1.,1.),1e-9;mode="levels")
        end
    finally
        gmsh.finalize()
    end
end

@testset "Signature validation and deterministic ownership" begin
    mktempdir() do root
        original=read(joinpath(@__DIR__,"testdata","one-edge-masked-signature.csv"),String)
        path=joinpath(root,"signature.csv")
        write(path,original)
        @test length(read_edges(path))==1
        write(path,replace(original,"1,0,1,0,0,0"=>"1,0.2,1,0,0,0"))
        @test_throws ErrorException read_edges(path)
        write(path,replace(original,"0,0,0,1,0,0,0,1,0,1"=>"0,0,0,2,0,0,0,1,0,1"))
        @test_throws ErrorException read_edges(path)
    end
    edges=[(point=(x,0.,0.),tangent=(0.,1.,0.),gap=(1.,0.,0.),
            interval=(0.,1.),vertex_arm=false,conductor=1,slot=slot)
           for (x,slot) in ((-1.,1),(1.,0))]
    @test nearest_edge(edges,(0.,0.5,0.),1.).slot==0
    @test nearest_edge(reverse(edges),(0.,0.5,0.),1.).slot==0
end

@testset "Requested submicron fillet is not silently ignored" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.model.add("small_fillet")
        box=(Int32(3),gmsh.model.occ.addBox(0.,0.,0.,0.5,0.2,0.08))
        before=gmsh.model.occ.getMass(3,box[2])
        rounded=fillet_plane_edges(gmsh.model.occ,[box],0.005,0.08,5e-8)
        @test sum(gmsh.model.occ.getMass(3,t) for (_,t) in rounded)<before*(1-1e-6)
        gmsh.model.occ.synchronize()
        @test any(gmsh.model.getType(2,t)!="Plane" for (_,t) in gmsh.model.getEntities(2))
        @test_throws ErrorException fillet_plane_edges(gmsh.model.occ,rounded,0.005,1.,5e-8)
    finally
        gmsh.finalize()
    end
end

@testset "Exact circular grading, including rigid transforms" begin
    point(t)=(cos(t),sin(t),0.)
    arc=grading_arc(point(0.),point(pi/8),point(pi/4),point(pi/2))
    @test grading_arc_distance2(point(pi/6),arc)<1e-25
    @test grading_arc_distance2((-1.,0.,0.),arc) ≈ 2.
    @test grading_arc_distance2((sqrt(.5),sqrt(.5),2.),arc) ≈ 4.
    circle=grading_arc(point(0.),point(pi/2),point(pi),point(2pi))
    @test grading_arc_distance2((0.,-2.,0.),circle) ≈ 1.
    for angle in (0.3,1.2)
        transform(p)=(cos(angle)*p[1]+sin(angle)*p[3]+4.,p[2]-3.,
                      -sin(angle)*p[1]+cos(angle)*p[3]+2.)
        rotated=grading_arc(transform(point(0.)),transform(point(pi/8)),
                            transform(point(pi/4)),transform(point(pi/2)))
        @test grading_arc_distance2(transform((-1.,0.,0.)),rotated) ≈ 2.
        @test grading_arc_distance2(transform((sqrt(.5),sqrt(.5),2.)),rotated) ≈ 4.
    end
    empty!(GRADING_SEGMENTS);empty!(GRADING_ARCS);push!(GRADING_ARCS,arc)
    GRADING_FINE[]=0.002;GRADING_GROWTH[]=1.;GRADING_FAR[]=0.5
    @test exact_grading_callback(2,0,1.,0.,0.01,1.,C_NULL) ≈ 0.012
    empty!(GRADING_ARCS)
end

@testset "Conic sizing is conservative without changing geometry" begin
    parameters=collect(range(0.,pi/2;length=17))
    point(t)=(2cos(t),sin(t),0.)
    conic,_,_,_=grading_conic(parameters,point.(parameters),0.02)
    @test conic.distance_error<=0.0002
    for t in (.1,.4,.9,1.4)
        normal=[cos(t)/2,sin(t),0.]
        normal/=norm(normal)
        p=Tuple(collect(point(t))+0.1normal)
        distance=sqrt(grading_conic_distance2(p,conic))
        @test 0.1-2conic.distance_error-1e-12<=distance<=0.1+1e-12
    end
    @test_throws ErrorException grading_conic(parameters,[(t,t^2,t^3) for t in parameters],0.02)
end

@testset "Explicit layer ownership, including repeated conductor labels" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.option.setNumber("Mesh.MeshSizeMin",0.2)
        gmsh.option.setNumber("Mesh.MeshSizeMax",0.2)
        gmsh.model.add("layer_labels")
        faces=[gmsh.model.occ.addRectangle(0.,0.,z,1.,1.) for z in (0.,0.01)]
        gmsh.model.occ.synchronize()
        gmsh.model.addPhysicalGroup(2,faces,4001)
        gmsh.model.mesh.generate(2)
        edges=[(point=(x,0.,z),tangent=(0.,1.,0.),gap=(1.,0.,0.),
                interval=(0.,1.),vertex_arm=false,conductor=1,slot=s,normal_sign=sign)
               for (x,z,s,sign) in ((0.,0.,0,1.),(1.,0.01,1,-1.))]
        mktempdir() do root
            label_interface_patches(edges,NamedTuple[],1.,joinpath(root,"slots.csv"))
            @test sort([Int(a) for (_,a) in gmsh.model.getPhysicalGroups(2)])==[4001,4101]
            for (attr,z) in ((4001,0.),(4101,0.01))
                for entity in gmsh.model.getEntitiesForPhysicalGroup(2,attr)
                    _,_,connect=gmsh.model.mesh.getElements(2,entity)
                    @test all(abs(gmsh.model.mesh.getNode(node)[1][3]-z)<1e-12 for node in unique(vcat(connect...)))
                end
            end
            @test_throws ErrorException label_interface_patches(edges,NamedTuple[],1.,
                joinpath(root,"overlap.csv");fabricated=true)
        end
    finally
        gmsh.finalize()
    end
end

@testset "Whole-element ownership certificate" begin
    edges=[(point=(x,0.,0.),tangent=(0.,1.,0.),gap=(1.,0.,0.),
            interval=(0.,1.),vertex_arm=false,conductor=1,slot=slot)
           for (x,slot) in ((0.,0),(1.,1))]
    ownership=build_interface_ownership(edges,NamedTuple[],1.)
    points=[(.4,.2,0.),(.51,.2,0.),(.4,.3,0.)]
    center=ntuple(d->sum(p[d] for p in points)/3,3)
    label=ownership.classify(4001,center)
    samples=[ntuple(d->.8p[d]+.2center[d],3) for p in points]
    @test all(ownership.classify(4001,p)==label for p in samples)
    @test ownership.classify(4001,points[2])!=label
    @test !ownership.certify(4001,center,points)
    stable=[(.1,.2,0.),(.2,.2,0.),(.1,.3,0.)]
    center=ntuple(d->sum(p[d] for p in stable)/3,3)
    @test ownership.certify(4001,center,stable)
    curved=vcat(stable,[(.49,.2,0.),(.15,.25,0.),(.1,.25,0.)])
    hull=interface_triangle_hull(curved)
    @test maximum(p[1] for p in hull)>.5
    @test !ownership.certify(4001,center,hull)
    @test interface_triangle_hull(vcat(curved,[(0.,0.,0.)]))===nothing
    # Long along the bisector, narrow across it: a radius-only test is needlessly
    # pessimistic, while the squared-distance polynomial certifies the label.
    narrow=[(.40,.1,0.),(.49,.1,0.),(.40,.9,0.)]
    center=ntuple(d->sum(p[d] for p in narrow)/3,3)
    @test ownership.certify(4001,center,narrow)
    rng=MersenneTwister(131)
    for _ in 1:60
        triangle=[(rand(rng),rand(rng),0.) for _ in 1:3]
        center=ntuple(d->sum(p[d] for p in triangle)/3,3)
        if ownership.certify(4001,center,triangle)
            label=ownership.classify(4001,center)
            @test all(ownership.classify(4001,ntuple(d->(i*triangle[1][d]+j*triangle[2][d]+(8-i-j)*triangle[3][d])/8,3))==label
                      for i in 0:8 for j in 0:8-i)
        end
    end
end

@testset "Bernstein ownership bounds and clamped endpoints" begin
    for weights in (OWNERSHIP_WEIGHTS_LINEAR,OWNERSHIP_WEIGHTS_QUADRATIC)
        totals=zeros(weights.count)
        for (_,_,k,w) in weights.pairs;totals[k]+=w;end
        @test all(abs.(totals.-1).<1e-14)
    end
    selected=ownership_segment((0.,0.,0.),(1.,0.,0.))
    other=ownership_segment((1.,-1.,1.),(1.,1.,1.))
    hull=[(.98,.05,.1),(1.02,.05,.1),(.98,.1,.1),
          (1.,.06,.1),(1.,.09,.11),(.98,.08,.1)]
    center=Tuple((collect(hull[1])+collect(hull[2])+collect(hull[3])+2collect(hull[4])+2collect(hull[5])+2collect(hull[6]))/9)
    @test certify_segment_label([selected,other],[0,1],0,center,hull,1e-9)
    rng=MersenneTwister(21)
    for _ in 1:60
        bary=rand(rng,3);bary/=sum(bary)
        a,b,c=bary
        weights=(a*a,b*b,c*c,2a*b,2b*c,2c*a)
        p=ntuple(d->sum(weights[i]*hull[i][d] for i in 1:6),3)
        d0=grading_segment_distance2(p,(selected.first...,selected.last...))
        d1=grading_segment_distance2(p,(other.first...,other.last...))
        @test d0<d1
    end
end

@testset "Exact weighted refinement-point query" begin
    rng=MersenneTwister(1729)
    points=[(rand(rng),rand(rng),rand(rng),.001+.1rand(rng)) for _ in 1:200]
    tree=SizePointTree(points)
    for _ in 1:100
        p=(2rand(rng)-.5,2rand(rng)-.5,2rand(rng)-.5)
        reference=min(.5,minimum(q[4]+sqrt(sum((p[d]-q[d])^2 for d in 1:3)) for q in points))
        @test query_size_point(tree,p,.5) ≈ reference atol=1e-14
    end
    @test query_size_point(SizePointTree(NTuple{4,Float64}[]),(0.,0.,0.),.5)==.5
    @test_throws ErrorException SizePointTree([(0.,0.,0.,-1.)])
end

@testset "Volume-only growth preserves the surface size field" begin
    empty!(GRADING_SEGMENTS);push!(GRADING_SEGMENTS,(0.,0.,0.,1.,0.,0.))
    empty!(GRADING_ARCS);empty!(GRADING_CONICS)
    SLOT_SIZE_TREE[]=nothing;GRADING_CAP_SIZE[]=0.;GRADING_TANGENT[]=0.
    GRADING_FINE[]=.002;GRADING_FAR[]=.5;GRADING_GROWTH[]=1.
    profile=(near_growth=.5,far_growth=2.,transition_distance=.03,maximum_size=.75)
    GRADING_VOLUME_PROFILE[]=profile
    @test volume_grading_size(.03,.002,profile) ≈ .017
    @test volume_grading_size(.10,.002,profile) ≈ .157
    @test volume_grading_size(1.,.002,profile)==.75
    independent=merge(profile,(minimum_size=.002,))
    @test volume_grading_size(0.,.001,independent)==.002
    @test volume_grading_size(.03,.001,independent) ≈ .017
    @test exact_grading_callback(2,0,.5,.10,0.,1.,C_NULL) ≈ .102
    @test exact_grading_callback(3,0,.5,.10,0.,1.,C_NULL) ≈ .157
    GRADING_VOLUME_PROFILE[]=nothing
    @test exact_grading_callback(3,0,.5,.10,0.,1.,C_NULL) ≈ .102
end

@testset "Named ARM-compatible size callback" begin
    empty!(GRADING_SEGMENTS);push!(GRADING_SEGMENTS,(0.,0.,0.,10.,0.,0.))
    empty!(GRADING_HORIZONTAL_CURVES);GRADING_TANGENT[]=0.
    GRADING_FINE[]=0.002;GRADING_FAR[]=0.5;GRADING_GROWTH[]=1.
    callback=@cfunction(exact_grading_callback,Cdouble,(Cint,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cvoid}))
    @test ccall(callback,Cdouble,(Cint,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cvoid}),1,1,4.123456,0.,0.,1.,C_NULL) ≈ 0.002 atol=1e-12
    @test ccall(callback,Cdouble,(Cint,Cint,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cvoid}),2,1,4.123456,.01,0.,1.,C_NULL) ≈ 0.012 atol=1e-12
    @test exact_grading_callback(3,1,4.,1.,1.,1.,C_NULL)==0.5
    GRADING_TANGENT[]=0.1;GRADING_HORIZONTAL_CURVES[Cint(1)]=GRADING_SEGMENTS[1]
    @test exact_grading_callback(1,Cint(1),5.,0.,0.,1.,C_NULL)==0.1
    @test exact_grading_callback(2,Cint(1),5.,0.,0.,1.,C_NULL)==0.002
    @test exact_grading_callback(1,Cint(1),0.,0.,0.,1.,C_NULL)==0.002
end
