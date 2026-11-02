<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# PalaceSensitivity.jl

Read the sensitivity export of a *Palace* run, and contract it with design velocities.

A run with a `"Sensitivity"` configuration writes, for each objective (for example an
entry of the capacitance matrix, an eigenfrequency or an S-parameter), its value and its
gradient with respect to the positions of the mesh nodes in an export region. Contracting
that gradient with a design velocity `V`, which says how the geometry moves with a design
parameter `p`, gives the derivative `dJ/dp = Σₖ gₖ ⋅ V(xₖ)`. Design velocities come from the
layout tool, for example `DeviceLayout.SolidModels.design_velocities`; this package needs no
knowledge of the geometry.

| Function                                             | Purpose                                                      |
|:---------------------------------------------------- |:------------------------------------------------------------ |
| `read_sensitivity(dir)`                              | read the export in `<Problem.Output>/sensitivity`            |
| `objective(s, index)`                                | an objective's type, configuration, value, unit and gradient |
| `jacobian(s, velocities)`, `derivative(s, index, V)` | `dJ/dp` for design velocities `V`                            |
| `write_sensitivity(dir; ...)`                        | write the format, for tests and other tools                  |

## Installation

The package lives in the *Palace* repository and is not registered. Add it to an
environment with

```julia
using Pkg
Pkg.add(url="https://github.com/awslabs/palace", subdir="scripts/julia/PalaceSensitivity")
```

or, from a clone of the repository, `Pkg.develop(path="scripts/julia/PalaceSensitivity")`.

## Usage

```julia
using PalaceSensitivity
s = read_sensitivity("postpro/sensitivity")
objective(s, 1).value                # the value of objective 1
J = jacobian(s, velocities)          # objectives × parameters
```

`velocities` is a collection of design velocities: callables mapping a `3 × N` matrix of
points to a `3 × N` matrix of velocities, with lengths in μm by default
(`velocity_length_unit`). For a `NamedTuple` of velocities, the columns follow its keys.

## Format (version 1)

| File                      | Content                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
|:------------------------- |:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sensitivity.json`        | `"Format": "PalaceSensitivity"`, `"Version": 1`, `"ProblemType"`, `"LengthUnit"` (the mesh unit `Model.L0`, in m), `"NumPoints"`, `"Points"`, `"ByteOrder": "LittleEndian"`, and `"Objectives"`: for each, `"Index"`, `"Type"`, its configuration (for example `"Terminals"`), `"Value"`, `"Unit"` and `"Gradient"` (the file name). Optional: `"Boundary"`, `"Region"`, `"PalaceVersion"`, `"MeshFingerprint"`, and scalar sensitivities in `"Materials"` and `"LumpedElements"`. |
| `points.bin`              | N × 3 Float64, row-major: the positions of the exported nodes, in mesh units                                                                                                                                                                                                                                                                                                                                                                                                       |
| `gradient-<index>.bin`    | N × 3 Float64 per objective: the gradient with respect to the node positions, in the objective's unit per mesh unit                                                                                                                                                                                                                                                                                                                                                                |
| `boundary.bin` (optional) | N + 1 Int64 offsets (0-based, compressed sparse rows), then the Int32 boundary attributes of each point; −1 marks the interface of a substructuring region                                                                                                                                                                                                                                                                                                                         |

Each node of the export region appears once. The set of points does not depend on *Palace*'s
partition, but their order does: a gradient is a covector at its points, and is only ever
contracted with velocities evaluated at the same points, never interpolated.

## Tests

```bash
julia --project=scripts/julia/PalaceSensitivity -e 'using Pkg; Pkg.test()'
```
