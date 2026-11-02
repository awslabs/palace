<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# PalaceSensitivity.jl

Read a *Palace* sensitivity export and contract its gradients with design velocities, giving
derivatives with respect to design parameters: `dJ/dp = Σₖ gₖ ⋅ V(xₖ)`, where `gₖ` is the
gradient of objective `J` at node `xₖ` and `V` is how the geometry moves with the parameter
`p` (for example from `DeviceLayout.SolidModels.design_velocities`).

```julia
using PalaceSensitivity
s = read_sensitivity("postpro/sensitivity")
objective(s, 1).value          # the value of objective 1
J = jacobian(s, velocities)    # objectives × parameters
```

`velocities` is a collection, such as a `NamedTuple`, of callables mapping a `3 × N` matrix
of points to a `3 × N` matrix of velocities, with lengths in μm by default.
`derivative(s, index, V)` gives one derivative, `check_velocity(s, V; fixed, rigid)` checks
that a velocity vanishes or is rigid on the given boundary attributes, and
`write_sensitivity` writes the format.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/awslabs/palace", subdir="scripts/julia/PalaceSensitivity")
```

## Format (version 1)

| File                      | Content                                                                                                    |
|:------------------------- |:---------------------------------------------------------------------------------------------------------- |
| `sensitivity.json`        | metadata, below                                                                                            |
| `points.bin`              | N × 3 Float64, row-major: the node positions, in mesh units                                                |
| `gradient-<index>.bin`    | N × 3 Float64 per objective: its gradient, in its unit per mesh unit                                       |
| `boundary.bin` (optional) | N + 1 Int64 offsets, then the Int32 boundary attributes of each point; −1 marks a substructuring interface |

`sensitivity.json` holds `"Format": "PalaceSensitivity"`, `"Version": 1`, `"ProblemType"`,
`"LengthUnit"` (the mesh unit, `Model.L0`, in m), `"NumPoints"`, `"Points"`,
`"ByteOrder": "LittleEndian"`, the optional `"Boundary"` and `"Region"` (the configured
export region), and `"Objectives"`. Each objective has `"Index"`, `"Type"`, `"Value"`,
`"Unit"`, `"Gradient"` (the file name) and the rest of its configuration, for example
`"Terminals"`, and the optional `"MaterialProperties"`: for each material, its `"Index"` in
`Domains.Materials`, its `"Attributes"`, and the derivatives with respect to a relative
change of its properties, such as `"Permittivity"` (ε ∂J/∂ε, in the objective's unit).
Without `"Shape"` in the configuration, there are no points and no gradient files. Other
entries are read as metadata.

The points are the nodes of the mesh *Palace* solves on, each once, in an order that depends
on the partition: a gradient is only contracted with velocities evaluated at its own points.
These are the nodes of the input mesh, except where *Palace* refines it or, for high-order
Gmsh meshes, re-interpolates the nodes onto its own basis. There, morphing the input mesh by
`V` moves the nodes by the interpolant of `V`, and the contraction agrees with finite
differences of morphed input meshes to the interpolation error of `V`.

## Tests

```bash
julia --project=scripts/julia/PalaceSensitivity -e 'using Pkg; Pkg.test()'
```
