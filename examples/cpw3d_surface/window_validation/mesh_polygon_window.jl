# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabricated reference mesh of a validation window given as a polygon set (see
# PolygonWindowMesh.jl for the JSON schema, the process geometry and the attribute table).
#
#   julia --project=. mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2
#       [--plan-only]
#
# Writes OUTPUT.msh2 and OUTPUT.json (the manifest). `--plan-only` builds and checks the plan
# mesh only and prints the manifest (no volume mesh).

include(joinpath(@__DIR__, "PolygonWindowMesh.jl"))
using .PolygonWindowMesh
using JSON

positional = filter(argument -> !startswith(argument, "--"), ARGS)
plan_only = "--plan-only" in ARGS
length(positional) == 4 || error(
    "Usage: mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2 [--plan-only]"
)
spec = read_polygon_set(positional[1])
manifest = mesh_polygon_window(
    spec,
    parse(Float64, positional[2]),
    parse(Float64, positional[3]),
    positional[4];
    plan_only=plan_only
)
if plan_only
    JSON.print(stdout, manifest, 2)
    println()
end
