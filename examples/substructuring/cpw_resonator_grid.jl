# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses DeviceLayout.jl to create the mesh of a chip with a grid of coplanar
waveguide (CPW) resonators, for the driven substructuring example (see
cpw_resonator_grid_substructuring.py and cpw_resonator_grid*.json). The resonators of each row
are capacitively coupled to a shared feedline with a lumped port at each end.

Each resonator is an open-ended (half-wave) CPW: a coupling section parallel to the feedline
followed by a meander. The resonators have different lengths, as for frequency-multiplexed
readout, and are shorter than typical readout resonators to keep the mesh small.

## Prerequisites

This script requires DeviceLayout.jl and its dependencies. If you don't already have it
installed, you can install it with (from this directory)

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From this directory, run:

```bash
julia --project -e 'include("cpw_resonator_grid.jl"); generate_cpw_resonator_grid()'
```

which writes mesh/cpw_resonator_grid.msh2.
=#

using FileIO
using DeviceLayout, DeviceLayout.SchematicDrivenLayout, DeviceLayout.PreferredUnits
import .SchematicDrivenLayout.ExamplePDK
import .SchematicDrivenLayout.ExamplePDK: LayerVocabulary

"""
    generate_cpw_resonator_grid(;
        rows = 3,
        cols = 4,
        lengths = collect(3000:100:4100)μm,
        coupling_gaps = fill(6μm, rows * cols),
        coupling_length = 300μm,
        bend_radius = 20μm,
        cell_x = 900μm,
        cell_y = 800μm,
        mesh_scale = 1.0,
        mesh_filename = "cpw_resonator_grid.msh2",
        mesh_order = 2
    )

Generate the mesh of `rows` x `cols` open-ended CPW resonators. The cell of resonator
k = (r - 1) * `cols` + c (row r, column c) is centered at x = (c - (cols + 1) / 2) * `cell_x`,
above the feedline of row r, along x at y = (r - 1) * `cell_y`. Resonator k has total length
`lengths[k]`, and its coupling section (of length `coupling_length`) is separated from the
feedline by a ground strip of width `coupling_gaps[k]`. The meander runs have the length of
the coupling section and turn with radius `bend_radius`. Each feedline ends on the ground
plane with a square lumped port near each end, as in DeviceLayout's SingleTransmon example: the
ports of row r are the physical groups "port_(2r - 1)" (left) and "port_(2r)" (right).
`mesh_scale` > 1 coarsens the mesh. The mesh is written to `mesh/mesh_filename` with physical
groups "metal", "substrate", "vacuum", "exterior_boundary" and the ports.
"""
function generate_cpw_resonator_grid(;
    rows=3,
    cols=4,
    lengths=collect(3000:100:4100)μm,
    coupling_gaps=fill(6μm, rows * cols),
    coupling_length=300μm,
    bend_radius=20μm,
    cell_x=900μm,
    cell_y=800μm,
    mesh_scale=1.0,
    mesh_filename="cpw_resonator_grid.msh2",
    mesh_order=2
)
    reset_uniquename!()
    @assert length(lengths) == rows * cols == length(coupling_gaps)
    cpw_width, cpw_gap = 10μm, 6μm
    style = Paths.SimpleCPW(cpw_width, cpw_gap)
    paths = Path[]
    for r = 1:rows
        origin, len = Point(0μm, (r - 1) * cell_y), cols * cell_x - 200μm
        push!(paths, feedline_path("feedline_$r", origin, len, style, cpw_width))
    end
    for r = 1:rows, c = 1:cols
        k = (r - 1) * cols + c
        origin = Point((c - (cols + 1) / 2) * cell_x, (r - 1) * cell_y)
        push!(
            paths,
            resonator_path(
                "resonator_$k",
                origin,
                lengths[k],
                coupling_gaps[k],
                coupling_length,
                bend_radius,
                style
            )
        )
    end
    # 300 μm of ground below the first feedline, 700 μm above the last one (the meanders
    # reach about 500 μm above their feedline).
    chip_x, chip_y = cols * cell_x + 200μm, (rows - 1) * cell_y + 1000μm
    area = centered(
        Rectangle(chip_x, chip_y),
        on_pt=Point(0μm, (rows - 1) * cell_y / 2 + 200μm)
    )
    render_mesh(
        paths,
        area,
        ["port_$i" for i = 1:(2 * rows)];
        mesh_filename,
        mesh_order,
        mesh_scale
    )
    return nothing
end

# A feedline along x centered at `origin`, ending on the ground plane with a square lumped
# port near each end (centers `cpw_width` from the ends, to avoid corner effects).
function feedline_path(name, origin, len, style, cpw_width)
    p = Path(
        origin - Point(len / 2, 0μm);
        α0=0,
        name=name,
        metadata=LayerVocabulary.METAL_NEGATIVE
    )
    straight!(p, len, style)
    csport = CoordinateSystem(uniquename("port"), nm)
    render!(
        csport,
        only_simulated(centered(Rectangle(cpw_width, cpw_width))),
        LayerVocabulary.PORT
    )
    attach!(p, sref(csport), cpw_width, i=1)
    attach!(p, sref(csport), len - cpw_width, i=1)
    return p
end

# An open-ended CPW resonator coupled to a feedline along x through `origin`: a coupling
# section of length `coupling_length` centered above `origin` and separated from the
# feedline by a ground strip of width `gap`, then a meander away from it.
function resonator_path(name, origin, len, gap, coupling_length, bend_radius, style)
    cpw_width, cpw_gap = Paths.trace(style), Paths.gap(style)
    y_c = cpw_width / 2 + 2 * cpw_gap + cpw_width / 2 + gap
    p = Path(
        origin + Point(-coupling_length / 2, y_c);
        α0=0,
        name=name,
        metadata=LayerVocabulary.METAL_NEGATIVE
    )
    straight!(p, coupling_length, style)
    turn!(p, π, bend_radius)
    meander!(p, len - coupling_length - π * bend_radius, coupling_length, bend_radius, -π)
    terminate!(p)
    terminate!(p; initial=true)
    return p
end

# Render the paths on a chip `area` with ExamplePDK's single-chip target (without the
# airbridge steps), keeping the given port physical groups, and save the mesh.
function render_mesh(paths, area, port_names; mesh_filename, mesh_order, mesh_scale)
    g = SchematicGraph("cpw-resonator-grid")
    for p in paths
        add_node!(g, p)
    end
    floorplan = plan(g)
    render!(floorplan.coordinate_system, area, LayerVocabulary.SIMULATED_AREA)
    render!(floorplan.coordinate_system, area, LayerVocabulary.WRITEABLE_AREA)
    render!(floorplan.coordinate_system, area, LayerVocabulary.CHIP_AREA)
    check!(floorplan)
    target = ExamplePDK.singlechip_solidmodel_target(port_names...)
    is_bridge(x) = x isa AbstractString && (startswith(x, "_") || occursin("bridge", x))
    filter!(op -> !any(is_bridge, (op[1], op[3]...)), target.postrenderer)

    sm = SolidModel("cpw-resonator-grid", overwrite=true)
    SolidModels.set_gmsh_option("General.Verbosity", 1)
    SolidModels.mesh_order(mesh_order)
    SolidModels.mesh_scale(mesh_scale)
    render!(sm, floorplan, target)
    SolidModels.gmsh.model.mesh.generate(3)
    mesh_path = joinpath(@__DIR__, "mesh", mesh_filename)
    mkpath(dirname(mesh_path))
    save(mesh_path, sm)
    println("Generated mesh: $(mesh_path)")
    return nothing
end
