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
`derivative(s, index, V)` gives one derivative, and `write_sensitivity` writes the format.

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
`"ByteOrder": "LittleEndian"`, the optional `"Boundary"`, and `"Objectives"`. Each objective
has `"Index"`, `"Type"`, `"Value"`, `"Unit"`, `"Gradient"` (the file name) and the rest of
its configuration, for example `"Terminals"`. Other entries, such as `"PalaceVersion"` or
material sensitivities, are read as metadata.

Each node of the export region appears once. Their order depends on *Palace*'s partition, so
a gradient is only contracted with velocities evaluated at its own points.

## Tests

```bash
julia --project=scripts/julia/PalaceSensitivity -e 'using Pkg; Pkg.test()'
```
