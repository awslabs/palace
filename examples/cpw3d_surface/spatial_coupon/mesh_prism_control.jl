# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fresh prism control for matched mesh benchmarks. Never modifies a retained mesh.
include(joinpath(@__DIR__,"mesh_spatial_coupon_swept.jl"))
using TOML
length(ARGS) in (6,7) || error("input_directory thin|fabricated output.msh normal_size tangent_size far_size [--full-gap-etch]")
full_gap=length(ARGS)==7 && ARGS[7]=="--full-gap-etch"
length(ARGS)==6 || full_gap || error("Unknown prism control option")
root=abspath(ARGS[1]);kind=ARGS[2];output=abspath(ARGS[3])
kind in ("thin","fabricated") || error("Invalid coupon kind")
isfile(output) && error("Refuse to overwrite prism control")
process=TOML.parsefile(joinpath(root,"process.toml"))
process["Units"]=="um" || error("Prism control expects mesh dimensions in um")
fine,tangent,far=parse.(Float64,ARGS[4:6])
generate_masked_edge_cluster_coupon(
    signature=joinpath(root,"mesh-signature.csv"),mask=joinpath(root,"plan-view-mask.csv"),
    boundary=joinpath(root,"plan-view-boundary.csv"),fabricated=kind=="fabricated",
    radius=Float64(process["Radius"]),metal_thickness=Float64(process["MetalThickness"]),
    overetch=Float64(process["Overetch"]),sidewall_angle=Float64(process["SidewallAngle"]),
    top_rounding=Float64(process["TopRounding"]),trench_rounding=Float64(process["TrenchRounding"]),
    lc_normal=fine,lc_tangent=tangent,lc_far=far,process_core_width=0.4Float64(process["Radius"]),
    process_grading_power=1.,normal_growth_ratio=1.4,max_nodes=2_000_000,
    max_elements=3_000_000,mesh_order=1,filename=output,etch_full_gap=full_gap)
