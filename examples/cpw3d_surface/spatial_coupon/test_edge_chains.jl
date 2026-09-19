# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
# Decision 47: the coupon box of a straight metal edge does not depend on how the
# source CAD subdivides it into signature rows; rows that merely share a line but do
# not touch (two edges separated by a slot, as on the ten-edge coupon) never chain.
using Test
include(joinpath(@__DIR__, "mesh_spatial_coupon.jl"))

const HEADER = "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm"

function signature(rows)
    path = tempname() * ".csv"
    open(path, "w") do io
        println(io, HEADER)
        for (i, row) in enumerate(rows)
            println(io, string(i, ",", row))
        end
    end
    return path
end

@testset "CAD subdivision of a straight edge keeps the coupon box" begin
    straight = read_edges(signature(["0,1,0,0,0,1,0,0,0,1,0,1,-2,2,0"]))
    halves = read_edges(signature(["0,1,0,-1,0,1,0,0,0,1,0,1,-1,1,0", "0,1,0,1,0,1,0,0,0,1,0,1,-1,1,0"]))
    @test isempty(EDGE_CHAIN_RECORDS) == false
    @test length(EDGE_CHAIN_RECORDS) == 1 && EDGE_CHAIN_RECORDS[1]["Rows"] == [1, 2]
    @test EDGE_CHAIN_RECORDS[1]["UnionLength"] ≈ 4.0
    @test collect(extended_interval(halves[1], 2.0)) ≈ [-5.0, 1.0]     # outer end extended by 2R
    @test collect(extended_interval(halves[2], 2.0)) ≈ [-1.0, 5.0]
    subdivided_box = coupon_bounds(halves, 2.0, 0.1, 0.05)
    read_edges(signature(["0,1,0,0,0,1,0,0,0,1,0,1,-2,2,0"]))            # re-register the straight row
    @test isempty(EDGE_CHAIN_RECORDS)
    @test collect(extended_interval(straight[1], 2.0)) ≈ [-6.0, 6.0]
    straight_box = coupon_bounds(straight, 2.0, 0.1, 0.05)
    @test all(collect(subdivided_box[1]) .≈ collect(straight_box[1]))
    @test all(collect(subdivided_box[2]) .≈ collect(straight_box[2]))
    @test collect(straight_box[1]) ≈ [-4.0, -8.0, -2.05] && collect(straight_box[2]) ≈ [4.0, 8.0, 2.1]
    # A chain shorter than 2 R is not extended (the single-row rule on the union).
    short = read_edges(signature(["0,1,0,-0.5,0,1,0,0,0,1,0,1,-0.5,0.5,0", "0,1,0,0.5,0,1,0,0,0,1,0,1,-0.5,0.5,0"]))
    @test length(EDGE_CHAIN_RECORDS) == 1
    @test collect(extended_interval(short[1], 2.0)) ≈ [-0.5, 0.5]
end

@testset "collinear rows separated by a gap or differing in conductor do not chain" begin
    gapped = read_edges(signature(["0,1,-6.5,-0.6,0,1,0,0,0,-1,0,1,0,2,0", "0,1,-6.5,0.4,0,1,0,0,0,-1,0,1,-2,0,0"]))
    @test isempty(EDGE_CHAIN_RECORDS)
    @test collect(extended_interval(gapped[1], 2.0)) ≈ [0.0, 6.0]
    other_conductor = read_edges(signature(["0,1,0,-1,0,1,0,0,0,1,0,1,-1,1,0", "0,2,0,1,0,1,0,0,0,1,0,1,-1,1,0"]))
    @test isempty(EDGE_CHAIN_RECORDS)
    vertex_arm = read_edges(signature(["0,1,0,-1,0,1,0,0,0,1,0,1,-1,1,0", "0,1,0,0,0,1,0,0,0,1,0,1,0,2,1"]))
    @test isempty(EDGE_CHAIN_RECORDS)
    @test collect(extended_interval(vertex_arm[2], 2.0)) ≈ [0.0, 6.0]
end
