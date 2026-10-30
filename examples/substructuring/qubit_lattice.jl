# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#=
# README

This Julia script uses DeviceLayout.jl to create the mesh of a chip with an n x n lattice of
simplified qubits, for the substructuring example (see qubit_lattice_substructuring.py and
qubit_lattice_*.json). Each qubit is a pair of rectangular metal pads in a rectangular
cutout of the ground plane.

## Prerequisites

This script requires DeviceLayout.jl and its dependencies. If you don't already have it
installed, you can install it with (from this directory)

```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## How to run

From this directory, run:

```bash
julia --project -e 'include("qubit_lattice.jl"); generate_qubit_lattice()'
```

which writes mesh/qubit_lattice.msh2. To change the design of the middle qubit, for example
its pad gap (see the substructuring documentation for re-meshing the region of the saved
model):

```bash
julia --project -e 'include("qubit_lattice.jl"); using DeviceLayout: μm;
    generate_qubit_lattice(center_pad_gap=60μm,
                           mesh_filename="qubit_lattice_redesign.msh2")'
```
=#

using FileIO
using DeviceLayout, DeviceLayout.SchematicDrivenLayout, DeviceLayout.PreferredUnits
import .SchematicDrivenLayout.ExamplePDK
import .SchematicDrivenLayout.ExamplePDK: LayerVocabulary

"""
    generate_qubit_lattice(;
        n = 5,
        pitch = 1000μm,
        cutout = (600μm, 400μm),
        pad = (240μm, 150μm),
        pad_gap = 40μm,
        center_pad_gap = pad_gap,
        mesh_size = 20μm,
        mesh_filename = "qubit_lattice.msh2",
        mesh_order = 2
    )

Generate the mesh of an `n` x `n` lattice of qubits with spacing `pitch`, centered at the
origin. Each qubit is a pair of pads of size `pad`, separated along x by `pad_gap`, in a
ground-plane cutout of size `cutout`; the middle qubit uses `center_pad_gap`. The metal
edges have mesh size `mesh_size`. The mesh is written to `mesh/mesh_filename` with physical
groups "metal", "substrate", "vacuum", and "exterior_boundary".
"""
function generate_qubit_lattice(;
    n=5,
    pitch=1000μm,
    cutout=(600μm, 400μm),
    pad=(240μm, 150μm),
    pad_gap=40μm,
    center_pad_gap=pad_gap,
    mesh_size=20μm,
    mesh_filename="qubit_lattice.msh2",
    mesh_order=2
)
    reset_uniquename!()
    chip = n * pitch
    g = SchematicGraph("qubit-lattice")
    floorplan = plan(g)
    cs = floorplan.coordinate_system
    sized(r) = MeshSized(mesh_size)(r)
    for i = 1:n, j = 1:n
        c = Point((i - (n + 1) / 2) * pitch, (j - (n + 1) / 2) * pitch)
        gap = (i == j == (n + 1) ÷ 2) ? center_pad_gap : pad_gap
        render!(
            cs,
            sized(centered(Rectangle(cutout...), on_pt=c)),
            LayerVocabulary.METAL_NEGATIVE
        )
        for s in (-1, 1)
            p = c + Point(s * (gap + pad[1]) / 2, 0μm)
            render!(
                cs,
                sized(centered(Rectangle(pad...), on_pt=p)),
                LayerVocabulary.METAL_POSITIVE
            )
        end
    end
    area = centered(Rectangle(chip, chip))
    render!(cs, area, LayerVocabulary.SIMULATED_AREA)
    render!(cs, area, LayerVocabulary.WRITEABLE_AREA)
    render!(cs, area, LayerVocabulary.CHIP_AREA)
    check!(floorplan)

    # Example single-chip target, without the airbridge post-rendering steps (no bridges).
    target = ExamplePDK.singlechip_solidmodel_target()
    is_bridge(x) = x isa AbstractString && (startswith(x, "_") || occursin("bridge", x))
    filter!(op -> !any(is_bridge, (op[1], op[3]...)), target.postrenderer)

    sm = SolidModel("qubit-lattice", overwrite=true)
    SolidModels.set_gmsh_option("General.Verbosity", 1)
    SolidModels.mesh_order(mesh_order)
    render!(sm, floorplan, target)
    SolidModels.gmsh.model.mesh.generate(3)
    mesh_path = joinpath(@__DIR__, "mesh", mesh_filename)
    mkpath(dirname(mesh_path))
    save(mesh_path, sm)
    println("Generated mesh: $(mesh_path)")
    return nothing
end
