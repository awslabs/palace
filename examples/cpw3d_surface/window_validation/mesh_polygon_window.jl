# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Fabricated reference mesh of a validation window given as a polygon set (see
# PolygonWindowMesh.jl for the JSON schema, the process geometry and the attribute table).
#
#   julia --project=. mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2
#       [--plan-only] [--band-mode own|gmsh] [--fan-turn-angle-deg A]
#       [--cross-plane-snap-um DELTA] [--exact-band-thickness] [--band-cap none|partition|curve]
#       [--sweep tensor|graded] [--alpha A] [--beta B] [--region-grading on|off]
#       [--region-ring on|off]
#
# Writes OUTPUT.msh2 and OUTPUT.json (the manifest). `--sweep tensor` (default) sweeps every
# plan triangle through every z level (the recorded reference family); `--sweep graded` (own
# band only) sweeps the graded (n, z) cross-section along the band columns (graded_sweep.jl):
# `--alpha A` (default 1) sets row k's z spacing to at least A x its width, `--beta B` (default 3)
# ends row k at B h_k beyond the fabricated step; `--region-grading` (graded only, default on)
# grades the Gmsh region's plan size from the station size t at the band tops up to the region
# size; `--region-ring` (graded only, default on) puts the region nodes adjacent to the band on
# the band's outermost z stack. `--plan-only` builds and checks the plan
# mesh only and prints the manifest (no volume mesh). `--band-mode own` (default) builds the
# structured boundary-layer band itself with the per-segment, per-side band cap
# (structured_band.jl); `gmsh` uses Gmsh's BoundaryLayer field as the recorded transmon
# generator did. `--fan-turn-angle-deg A` (default 90, at most 120): outward corners turning
# more than A degrees get a fan of columns instead of the scaled bisector column.
# `--cross-plane-snap-um` overrides the cross-plane snap distance (default 0.05 x the set's
# MatchingRadius; recorded in the manifest). Gmsh mode only: `--exact-band-thickness` passes
# the exact geometric sum as the boundary-layer Thickness (the recorded transmon generator's
# formula, platform-fragile; the default adds a 1e-6 margin, decision 188); `--band-cap`
# applies the local band cap rule (default none: recorded only; see PolygonWindowMesh.jl).

include(joinpath(@__DIR__, "PolygonWindowMesh.jl"))
using .PolygonWindowMesh
using JSON

positional = String[]
plan_only = false
exact_band_thickness = false
cross_plane_snap_um = NaN
band_cap = :none
band_mode = :own
fan_turn_angle_deg = 90.0
sweep = :tensor
alpha = PolygonWindowMesh.DEFAULT_GRADED_ALPHA
beta = PolygonWindowMesh.DEFAULT_GRADED_BETA
region_grading = nothing # nothing: the sweep's default (graded on, tensor off)
region_ring = true
function parse_switch(option, value)
    value in ("on", "off") || error("$option takes on or off, not $value")
    return value == "on"
end
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
        elseif argument == "--band-cap"
            global band_cap = Symbol(ARGS[i + 1])
            i += 1
        elseif argument == "--band-mode"
            global band_mode = Symbol(ARGS[i + 1])
            i += 1
        elseif argument == "--fan-turn-angle-deg"
            global fan_turn_angle_deg = parse(Float64, ARGS[i + 1])
            i += 1
        elseif argument == "--sweep"
            global sweep = Symbol(ARGS[i + 1])
            i += 1
        elseif argument == "--alpha"
            global alpha = parse(Float64, ARGS[i + 1])
            i += 1
        elseif argument == "--beta"
            global beta = parse(Float64, ARGS[i + 1])
            i += 1
        elseif argument == "--region-grading"
            global region_grading = parse_switch(argument, ARGS[i + 1])
            i += 1
        elseif argument == "--region-ring"
            global region_ring = parse_switch(argument, ARGS[i + 1])
            i += 1
        elseif startswith(argument, "--")
            error(
                "Unknown option $argument; known: --plan-only --band-mode own|gmsh " *
                "--fan-turn-angle-deg A --cross-plane-snap-um DELTA --exact-band-thickness " *
                "--band-cap none|partition|curve --sweep tensor|graded --alpha A --beta B " *
                "--region-grading on|off --region-ring on|off"
            )
        else
            push!(positional, argument)
        end
        i += 1
    end
end
length(positional) == 4 || error(
    "Usage: mesh_polygon_window.jl WINDOW.json RADIAL_UM TANGENTIAL_UM OUTPUT.msh2 " *
    "[--plan-only] [--band-mode own|gmsh] [--fan-turn-angle-deg A] " *
    "[--cross-plane-snap-um DELTA] [--exact-band-thickness] [--band-cap none|partition|curve] " *
    "[--sweep tensor|graded] [--alpha A] [--beta B] [--region-grading on|off] " *
    "[--region-ring on|off]"
)
spec = read_polygon_set(positional[1])
manifest = mesh_polygon_window(
    spec,
    parse(Float64, positional[2]),
    parse(Float64, positional[3]),
    positional[4];
    plan_only=plan_only,
    exact_band_thickness=exact_band_thickness,
    cross_plane_snap_um=cross_plane_snap_um,
    band_cap=band_cap,
    band_mode=band_mode,
    fan_turn_angle_deg=fan_turn_angle_deg,
    sweep=sweep,
    alpha=alpha,
    beta=beta,
    region_grading=region_grading === nothing ? sweep == :graded : region_grading,
    region_ring=region_ring
)
if plan_only
    JSON.print(stdout, manifest, 2)
    println()
end
