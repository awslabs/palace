# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Test
using LinearAlgebra
include(joinpath(@__DIR__, "mesh_metric_tet_experiment.jl"))

@testset "Exact edge tensor geometry" begin
    segment=(0.,0.,0.,1.,0.,0.)
    metric=edge_metric((.5,0.,0.),segment,.005,.4,1.,4.,.25)
    @test eigvals(Symmetric(metric)) ≈ [2500.,40000.,40000.]
    @test edge_metric((0.,0.,0.),segment,.005,.4,1.,4.,.25) ≈ Matrix{Float64}(I,3,3)*40000
    angle=.37
    rotation=[cos(angle) -sin(angle) 0.;sin(angle) cos(angle) 0.;0. 0. 1.]
    shift=[2.,-1.,.4]
    point=rotation*[.5,0.,0.]+shift
    first=rotation*collect(segment[1:3])+shift
    last=rotation*collect(segment[4:6])+shift
    transformed=edge_metric(Tuple(point),(first...,last...),.005,.4,1.,4.,.25)
    @test transformed ≈ rotation*metric*rotation'
    @test minimum(eigvals(Symmetric(transformed)))>0
end
