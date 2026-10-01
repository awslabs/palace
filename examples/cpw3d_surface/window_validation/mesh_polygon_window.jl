# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabricated reference mesh of a validation window given as a polygon set (see
# PolygonWindowMesh.jl for the JSON schema, the process geometry and the attribute table).
#
#   julia --project=. mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2
#       [--plan-only] [--exact-band-thickness] [--cross-plane-snap-um DELTA]
#
# Writes OUTPUT.msh2 and OUTPUT.json (the manifest). `--plan-only` builds and checks the plan
# mesh only and prints the manifest (no volume mesh). `--exact-band-thickness` passes the
# exact geometric sum as the boundary-layer Thickness (the recorded transmon generator's
# formula, platform-fragile; the default adds a 1e-6 margin, decision 188).
# `--cross-plane-snap-um` overrides the cross-plane snap distance (default 0.025 x the set's
# MatchingRadius; recorded in the manifest).

include(joinpath(@__DIR__, "PolygonWindowMesh.jl"))
using .PolygonWindowMesh
using JSON

positional = String[]
plan_only = false
exact_band_thickness = false
cross_plane_snap_um = NaN
let i = 1
    while i <= length(ARGS)
        argument = ARGS[i]
        if argument == "--plan-only"
            global plan_only = true
        elseif argument == "--exact-band-thickness"
            global exact_band_thickness = true
        elseif argument == "--cross-plane-snap-um"
            global cross_plane_snap_um = parse(Float64, ARGS[i + 1])
            i += 1
        elseif startswith(argument, "--")
            error(
                "Unknown option $argument; known: --plan-only --exact-band-thickness " *
                "--cross-plane-snap-um DELTA"
            )
        else
            push!(positional, argument)
        end
        i += 1
    end
end
length(positional) == 4 || error(
    "Usage: mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2 " *
    "[--plan-only] [--exact-band-thickness] [--cross-plane-snap-um DELTA]"
)
spec = read_polygon_set(positional[1])
manifest = mesh_polygon_window(
    spec,
    parse(Float64, positional[2]),
    parse(Float64, positional[3]),
    positional[4];
    plan_only=plan_only,
    exact_band_thickness=exact_band_thickness,
    cross_plane_snap_um=cross_plane_snap_um
)
if plan_only
    JSON.print(stdout, manifest, 2)
    println()
end
