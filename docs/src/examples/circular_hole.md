```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

```@setup include_example
function include_example_file(example_path, filename)
    print(read(joinpath(@__DIR__, "..", "..", "..", "test", "examples", "ref", example_path, filename), String))
end
```

# Flux Trapping Analysis

## Problem description

This example demonstrates the use of flux boundary conditions in magnetostatics for
analyzing magnetic flux trapping in superconducting structures. Specifically, we consider
cases where fixed amounts of magnetic flux are prescribed to holes in metallic planes and
compute the magnetic response of the system.

The flux conditions are imposed as integral boundary conditions on the hole perimeters using:

```math
\oint_h \mathbf{A} \cdot d\boldsymbol{\ell} = \Phi,
```

where ``\mathbf{A}`` is the magnetic vector potential and ``\Phi`` is the prescribed flux
through hole ``h``. The solver first solves a 2D surface curl problem on the metallic plane
to find the tangential component of ``\mathbf{A}`` satisfying the integral constraint, then
uses this as a boundary condition for the full 3D magnetostatic problem.

The inductance matrix is extracted from the stored magnetic energy:

```math
M_{ij} = \frac{\mathbf{A}_j^T K \mathbf{A}_i}{\Phi_i \Phi_j},
```

where ``K`` is the curl-curl stiffness matrix. For multiple flux loops, the mutual
inductance between loops characterizes the magnetic coupling.

Four configurations are provided in the
[`examples/circular_hole/`](https://github.com/awslabs/palace/blob/main/examples/circular_hole)
directory:

  - **Single hole** (`circular_hole.json`): One circular hole in a circular metal plate.
  - **Two holes, single flux loop** (`double_circular_hole.json`): Two holes on one metal
    plate with opposite flux prescribed through a single excitation.
  - **Two holes, separate flux loops** (`double_circular_hole_multi_flux.json`): Two holes
    on one metal plate with independent flux loop excitations.
  - **Two holes on separate planes** (`double_circular_hole_multi_planes.json`): Two holes
    on spatially separated metal plates with independent excitations.

The mesh length unit is ``\mu\text{m}`` in all configurations.

## Configuration

All configurations use
[`"Problem": {"Type": "Magnetostatic"}`](../config/problem.md#config-problem) and specify
flux loop boundaries via the
[`"FluxLoop"`](../config/boundaries.md#boundaries-fluxloop) keyword. The key solver
settings shared across all examples are:

```json
"Solver":
{
    "Order": 2,
    "Device": "CPU",
    "Magnetostatic":
    {
        "Save": 2
    },
    "Linear":
    {
        "Type": "AMS",
        "KSPType": "CG",
        "Tol": 1.0e-8,
        "MaxIts": 200
    }
}
```

The AMS (Auxiliary-space Maxwell Solver) preconditioner is used with CG iteration, which is
well-suited for the symmetric positive-definite curl-curl system arising in magnetostatics.

### Single hole

The first configuration (`circular_hole.json`) considers a circular metal plate with radius
``R = 3\,\mu\text{m}`` containing a concentric circular hole with radius
``r = 1\,\mu\text{m}``. A single flux quantum ``\Phi_0 = 2.068 \times 10^{-15}\,\text{Wb}``
is prescribed through the hole.

The `"FluxLoop"` boundary specification is:

```json
"FluxLoop":
[
    {
        "Index": 1,
        "FluxLoopPEC": [8],
        "HoleAttributes": [9],
        "FluxAmounts": [1.0],
        "Direction": "+Z"
    }
]
```

Here `"FluxLoopPEC"` specifies the boundary attributes of the metal surface on which the
2D surface curl problem is solved. `"HoleAttributes"` lists the boundary attributes of the
hole perimeters, and `"FluxAmounts"` gives the prescribed flux (in units of ``\Phi_0``)
through each hole. The `"Direction"` sets the normal direction for flux orientation.

### Two holes with single flux loop

The second configuration (`double_circular_hole.json`) uses a rectangular metal plate with
two circular holes separated by ``5\,\mu\text{m}``. Both holes are assigned to a single
flux loop excitation with opposite flux:

```json
"FluxLoop":
[
    {
        "Index": 1,
        "FluxLoopPEC": [8],
        "HoleAttributes": [9, 10],
        "FluxAmounts": [1.0, -1.0],
        "Direction": "+Z"
    }
]
```

This models the case where one flux quantum enters through one hole and exits through the
other, as might occur when a vortex-antivortex pair is trapped in a superconducting film.
Since both holes are part of a single excitation, the solver computes a single self-inductance
value.

### Two holes with separate flux loops

The third configuration (`double_circular_hole_multi_flux.json`) uses the same geometry but
treats each hole as an independent flux loop excitation:

```json
"FluxLoop":
[
    {
        "Index": 1,
        "FluxLoopPEC": [8],
        "HoleAttributes": [9],
        "FluxAmounts": [1.0],
        "Direction": [0.0, 0.0, 1.0]
    },
    {
        "Index": 2,
        "FluxLoopPEC": [8],
        "HoleAttributes": [10],
        "FluxAmounts": [1.0],
        "Direction": [0.0, 0.0, 1.0]
    }
]
```

With two independent excitations, the solver computes a full ``2 \times 2`` inductance
matrix including self-inductances ``M_{11}``, ``M_{22}`` and mutual inductance ``M_{12}``.
Note that `"Direction"` can be specified as either a string (`"+Z"`) or a numeric array
(`[0.0, 0.0, 1.0]`).

### Two holes on separate planes

The fourth configuration (`double_circular_hole_multi_planes.json`) places the two holes on
spatially separated metal plates, each with its own `"FluxLoopPEC"` attribute:

```json
"FluxLoop":
[
    {
        "Index": 1,
        "FluxLoopPEC": [8],
        "HoleAttributes": [9],
        "FluxAmounts": [1.0],
        "Direction": [0.0, 0.0, 1.0]
    },
    {
        "Index": 2,
        "FluxLoopPEC": [10],
        "HoleAttributes": [11],
        "FluxAmounts": [1.0],
        "Direction": [0.0, 0.0, 1.0]
    }
]
```

This tests the solver's ability to handle multiple independent metal surfaces. Since the
planes are spatially separated, the mutual inductance between the two loops is expected to
be smaller than in the shared-plane configuration.

## Mesh

The meshes are generated using Julia scripts with the Gmsh package, located in the `mesh/`
subdirectory. For example, the two-hole geometry can be generated with:

```bash
julia -e 'include("mesh/sheet_w_two_holes.jl"); generate_sheet_with_two_holes_mesh()'
```

The figure below shows the mesh for the rectangular plate with two circular holes and the
outer domain boundary:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/sheet_w_two_holes.png" width="60%" />
</p><br/>
```

## Results

### Single hole

For the single-hole configuration, the solver extracts a self-inductance of
``M = 2.808 \times 10^{-12}\,\text{H}`` for one flux quantum trapped in the
hole. The magnetic energy stored in the system is
``E_\text{mag} = 1.679 \times 10^{-13}\,\text{J}``.

The figures below show the magnetic vector potential amplitude ``|\mathbf{A}|``, the 
in-plane components of the vector potential ``A_x``, and ``A_y``, the magnetic field 
magnitude ``|\mathbf{B}|``, and surface currents ``J_x`` and ``J_y`` on the metal 
surface for the single-hole case:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Amagnitude_singlehole.png" width="27%" />
  <img src="../../assets/examples/Ax_inplane_singlehole.png" width="27%" />
  <img src="../../assets/examples/Ay_inplane_singlehole.png" width="27%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Bz_surface_singlehole.png" width="27%" />
  <img src="../../assets/examples/Js_x_singlehole.png" width="27%" />
  <img src="../../assets/examples/Js_y_singlehole.png" width="27%" />
</p><br/>
```

### Two holes with single flux loop


The figures below show the magnetic vector potential amplitude ``|\mathbf{A}|``, the 
in-plane components of the vector potential ``A_x``, and ``A_y``, the magnetic field 
``B_z``, and surface currents ``J_x`` and ``J_y`` on the metal 
surface for the two-hole geometry:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Amagnitude_doublehole.png" width="27%" />
  <img src="../../assets/examples/Ax_inplane_doublehole.png" width="27%" />
  <img src="../../assets/examples/Ay_inplane_doublehole.png" width="27%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Bz_surface_doublehole.png" width="27%" />
  <img src="../../assets/examples/Js_x_doublehole.png" width="27%" />
  <img src="../../assets/examples/Js_y_doublehole.png" width="27%" />
</p><br/>
```

### Two holes with separate flux loops

For the two-hole configuration with independent excitations, the computed inductance matrix
is:

```math
M = \begin{pmatrix}
1.372 & -1.733 \times 10^{-6} \\
-1.733 \times 10^{-6} & 1.389
\end{pmatrix} \text{pH}
```

The self-inductances are approximately equal (as expected by symmetry), and the mutual
inductance is negligible compared to the self-inductance. This indicates weak magnetic
coupling between the two holes on the same plane at this separation distance.

