# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Tiny synthetic two-level window for the mesher's checks: two CPW stubs facing each other
# across the flip-chip gap (an L1 stub open-ended under the L2 stub's open end, so the two
# planes' metal edges cross in plan) and one bump joining the two grounds. Every conductor is
# labelled: `ground` (both planes + the bump), `trace_l1` and `trace_l2` (terminals).
#
#   julia --project=. synthetic_two_level_window.jl OUTPUT.json [--scale 1.0]
#
# Box 80 x 60 um; L1 (facing up, z = 0, substrate 20 um) trace [10, 40] x [28, 32] in the gap
# hole [4, 46] x [22, 38]; L2 (facing down, z = 4.8, substrate 20 um) trace [40, 70] x
# [28, 32] in the hole [34, 76] x [22, 38]; bump: 16-gon of radius 4 at (15, 8), on ground
# of both planes. No vacuum beyond the substrates (truncated box, Neumann backsides).

using JSON

function rectangle(x0, x1, y0, y1)
    return [[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
end

function regular_polygon(cx, cy, radius, n)
    return [
        [cx + radius * cos(2pi * (k - 1) / n), cy + radius * sin(2pi * (k - 1) / n)] for
        k = 1:n
    ]
end

function synthetic_two_level_window(; scale=1.0, gap_z=4.8, substrate=20.0)
    s = scale
    return Dict(
        "Version" => 1,
        "Name" => "synthetic-two-level",
        "Box" => Dict("X" => [0.0, 80.0s], "Y" => [0.0, 60.0s]),
        "Process" => Dict("MetalThickness" => 0.1, "Overetch" => 0.05),
        "Planes" => [
            Dict(
                "Name" => "L1",
                "SurfaceZ" => 0.0,
                "Facing" => "up",
                "SubstrateThickness" => substrate,
                "Polygons" => [
                    Dict(
                        "Conductor" => "ground",
                        "Outer" => rectangle(0.0, 80.0s, 0.0, 60.0s),
                        "Holes" => [rectangle(4.0s, 46.0s, 22.0s, 38.0s)]
                    ),
                    Dict(
                        "Conductor" => "trace_l1",
                        "Outer" => rectangle(10.0s, 40.0s, 28.0s, 32.0s),
                        "Holes" => []
                    )
                ]
            ),
            Dict(
                "Name" => "L2",
                "SurfaceZ" => gap_z,
                "Facing" => "down",
                "SubstrateThickness" => substrate,
                "Polygons" => [
                    Dict(
                        "Conductor" => "ground",
                        "Outer" => rectangle(0.0, 80.0s, 0.0, 60.0s),
                        "Holes" => [rectangle(34.0s, 76.0s, 22.0s, 38.0s)]
                    ),
                    Dict(
                        "Conductor" => "trace_l2",
                        "Outer" => rectangle(40.0s, 70.0s, 28.0s, 32.0s),
                        "Holes" => []
                    )
                ]
            )
        ],
        "Bumps" => [
            Dict("Conductor" => "ground", "Footprint" => regular_polygon(15.0s, 8.0s, 4.0s, 16))
        ],
        "Vacuum" => Dict("Below" => 0.0, "Above" => 0.0),
        "Terminals" => ["trace_l1", "trace_l2"]
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    positional = filter(argument -> !startswith(argument, "--"), ARGS)
    length(positional) == 1 ||
        error("Usage: synthetic_two_level_window.jl OUTPUT.json [--scale S]")
    index = findfirst(==("--scale"), ARGS)
    scale = isnothing(index) ? 1.0 : parse(Float64, ARGS[index + 1])
    open(positional[1], "w") do stream
        JSON.print(stream, synthetic_two_level_window(; scale=scale), 2)
        println(stream)
    end
end
