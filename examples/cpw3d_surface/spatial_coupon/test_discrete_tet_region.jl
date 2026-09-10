# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
include(joinpath(@__DIR__,"mesh_discrete_tet_region.jl"))

@testset "Fixed shell bulk meshing" begin
    xy=[(0.,0.),(1.,0.),(1.,1.),(0.,1.)]
    tri=[(1,2,3),(1,3,4)]
    p,t=extruded_plan_shell(xy,tri,[0.,0.3,1.])
    @test check_closed_shell(p,t) ≈ 1.
    @test_throws ErrorException check_closed_shell(p,t[2:end])
    @test_throws ErrorException check_closed_shell(p,vcat(t,[t[1]]))
    @test_throws ErrorException extruded_plan_shell(xy,tri,[0.,0.])
    @test_throws ErrorException extruded_plan_shell(xy,[(1,3,2),(1,3,4)],[0.,1.])
    mktempdir() do root
        result=mesh_fixed_shell(p,t,joinpath(root,"bulk.msh");maximum_size=0.25)
        @test result.volume ≈ 1.
        @test result.tetrahedra>0
        @test result.boundary_nodes==12
        @test_throws ErrorException mesh_fixed_shell(p,t,joinpath(root,"bulk.msh"))
    end
end
