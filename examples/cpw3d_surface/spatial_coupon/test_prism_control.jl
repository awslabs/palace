# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
include(joinpath(@__DIR__,"mesh_spatial_coupon_swept.jl"))

@testset "Explicit full-gap process option" begin
    coordinates=Dict{UInt64,NTuple{2,Float64}}(1=>(10.,0.),2=>(11.,0.),3=>(10.,1.))
    triangle=(UInt64(1),UInt64(2),UInt64(3))
    primitives=[(kind=:line,first=(0.,0.),last=(0.,1.),conductor=1)]
    @test triangle_phase(triangle,coordinates,NamedTuple[],primitives,0.,1.,1e-9)==(:ordinary,0)
    @test triangle_phase(triangle,coordinates,NamedTuple[],primitives,0.,1.,1e-9;etch_full_gap=true)==(:trench,1)
    @test interval_material(true,:ordinary,-.025,0.,.1,.05)==:substrate
    @test interval_material(true,:trench,-.025,0.,.1,.05)==:vacuum
end

@testset "Unused metal-interior nodes are not exported" begin
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.model.add("prism_control")
        coordinates=Dict{UInt64,NTuple{2,Float64}}(1=>(0.,0.),2=>(1.,0.),3=>(1.,1.),4=>(0.,1.))
        triangles=[(UInt64(1),UInt64(2),UInt64(3)),(UInt64(1),UInt64(3),UInt64(4))]
        edges=[(point=(0.,0.,0.),conductor=1,slot=0,tangent=(0.,1.,0.),gap=(1.,0.,0.),
                interval=(0.,1.),vertex_arm=false,normal_sign=1.)]
        facets=[(conductor=1,plane=0.,points=[(0.,0.),(1.,0.),(1.,1.),(0.,1.)])]
        plan=(node_tags=UInt64[1,2,3,4],coordinates=coordinates,triangles=triangles,
              fabrication_primitives=[(kind=:line,first=(0.,0.),last=(0.,1.),conductor=1)],
              active_segments=active_plan_segments(edges,1.))
        result=add_discrete_edge_cluster_mesh!(plan,edges,facets,true,[-.1,0.,.05,.1,.2],
            (0.,0.,-.1),(1.,1.,.2),1.,.1,.05,1000,1000,1)
        @test result.node_count==16
        @test result.volume_count==4
        mktempdir() do directory
            path=joinpath(directory,"control.msh")
            gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
            gmsh.write(path);gmsh.clear();gmsh.open(path)
            @test length(gmsh.model.mesh.getNodes()[1])==16
        end
    finally
        gmsh.finalize()
    end
end

@testset "Ten-edge immutable source has no ordinary fabricated support" begin
    source = joinpath(@__DIR__, "testdata", "ten-edge-6791f1c84123")
    process = read(joinpath(source, "process-library.json"), String)
    source_number(name) = parse(Float64, only(match(
        Regex("\\\"$name\\\"\\s*:\\s*([0-9.eE+-]+)"), process).captures))
    phases = audit_source_plan_phases(
        joinpath(source, "mesh-signature.csv"),
        joinpath(source, "plan-view-mask.csv"),
        joinpath(source, "plan-view-boundary.csv");
        radius=2.0,
        metal_thickness=source_number("MetalThickness"),
        overetch=source_number("OveretchDepth"),
        lc_normal=0.025,
        lc_tangent=0.1,
        lc_far=0.16,
        process_core_width=0.2,
        normal_growth_ratio=1.4
    )
    @test phases[:metal] > 0
    @test phases[:trench] > 0
    @test phases[:ordinary] == 0 # Therefore attributes 3000/3001 cannot be emitted.
end
