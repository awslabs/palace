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

This example demonstrates Palace's flux boundary conditions for magnetostatic analysis of
magnetic flux trapping in superconducting structures. The problem considers metallic planes
containing holes through which fixed amounts of magnetic flux are prescribed, and computes
the resulting magnetic field distribution and inductance matrix.

Flux conditions are imposed as integral constraints on the hole perimeters:

```math
\oint_h \mathbf{A} \cdot d\boldsymbol{\ell} = \Phi,
```

where ``\mathbf{A}`` is the magnetic vector potential and ``\Phi`` is the prescribed flux
through hole ``h``. The solution proceeds in two stages: first, a 2D surface curl problem
is solved on the metallic plane to determine the tangential component of ``\mathbf{A}``
satisfying the integral constraint; then, this surface field serves as a Dirichlet boundary
condition for the full 3D magnetostatic problem.

Four configurations of increasing complexity are provided in the
[`examples/circular_hole/`](https://github.com/awslabs/palace/blob/main/examples/circular_hole)
directory:

  - **Single hole** (`circular_hole.json`): A circular hole in a circular metal plate.
  - **Two holes, single flux loop** (`double_circular_hole.json`): Two holes on one plate
    with equal and opposite flux prescribed through a single excitation.
  - **Two holes, separate flux loops** (`double_circular_hole_multi_flux.json`): Two holes
    on one plate with independent flux loop excitations.
  - **Two holes on separate planes** (`double_circular_hole_multi_planes.json`): Two holes
    on spatially separated plates with independent excitations.

All configurations use a mesh length unit of ``\mu\text{m}``.

## Configuration

Each configuration uses
[`"Problem": {"Type": "Magnetostatic"}`](../config/reference.md#config-problem) and specifies
flux loop boundaries via the
[`"FluxLoop"`](../config/reference.md#config-boundaries) keyword. The shared solver
settings are:

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

The AMS (Auxiliary-space Maxwell Solver) preconditioner with CG iteration is well-suited for
the symmetric positive-definite curl-curl system arising in magnetostatics.

### Single hole

The first configuration (`circular_hole.json`) models a circular metal plate of radius
``R = 3\,\mu\text{m}`` with a concentric hole of radius ``r = 1\,\mu\text{m}``. One
flux quantum ``\Phi_0 = 2.068 \times 10^{-15}\,\text{Wb}`` is prescribed through the hole.

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

The fields have the following meaning:

  - `"FluxLoopPEC"`: boundary attributes of the metal surface on which the 2D surface curl
    problem is solved.
  - `"HoleAttributes"`: boundary attributes of the hole perimeters where integral
    constraints are applied.
  - `"FluxAmounts"`: prescribed flux through each hole, in units of ``\Phi_0``.
  - `"Direction"`: surface normal direction for flux orientation.

### Two holes with single flux loop

The second configuration (`double_circular_hole.json`) uses a rectangular plate with two
circular holes separated by ``5\,\mu\text{m}``. Both holes belong to a single flux loop
excitation with opposite flux values:

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

This configuration models the scenario where one flux quantum enters through one hole and
exits through the other, as occurs when a vortex-antivortex pair is trapped in a
superconducting film. Since both holes share a single excitation index, the solver computes
one self-inductance value for the combined configuration.

### Two holes with separate flux loops

The third configuration (`double_circular_hole_multi_flux.json`) uses the same two-hole
geometry but treats each hole as an independent excitation:

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

With two independent excitations, the solver computes the full ``2 \times 2`` inductance
matrix, including self-inductances ``M_{11}``, ``M_{22}`` and mutual inductance ``M_{12}``.
Note that `"Direction"` accepts either a string shorthand (`"+Z"`) or an explicit numeric
array (`[0.0, 0.0, 1.0]`).

### Two holes on separate planes

The fourth configuration (`double_circular_hole_multi_planes.json`) places each hole on its
own spatially separated metal plate, each with a distinct `"FluxLoopPEC"` attribute:

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

Because the plates are physically separated, this configuration tests the solver's handling
of multiple independent metal surfaces and provides a comparison case where inter-hole
coupling is expected to be reduced.

## Mesh

The meshes are generated using Julia scripts with the Gmsh package, located in the `mesh/`
subdirectory. For example, the two-hole rectangular plate mesh can be generated with:

```bash
julia -e 'include("mesh/sheet_w_two_holes.jl"); generate_sheet_with_two_holes_mesh()'
```

The figures below show the three mesh geometries. From left to right: the single circular
hole on a circular plate, the rectangular plate with two holes, and the two spatially
separated square plates each containing a hole:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/circular_hole_mesh.png" width="30%" />
  <img src="../../assets/examples/sheet_w_two_holes.png" width="30%" />
  <img src="../../assets/examples/two_square_sheets.png" width="30%" />
</p><br/>
```

## Results

### Single hole

For the single-hole configuration, the solver extracts a self-inductance of
``M = 2.808\,\text{pH}`` for one flux quantum trapped in the hole, with a stored magnetic
energy of ``E_\text{mag} = 1.680 \times 10^{-13}\,\text{J}``.

The figures below show the magnetic vector potential amplitude ``|\mathbf{A}|``, its
in-plane components ``A_x`` and ``A_y``, the out-of-plane magnetic field ``B_z``, and the
surface current components ``J_x`` and ``J_y`` on the metal surface:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Amagnitude_singlehole.png" width="30%" />
  <img src="../../assets/examples/Ax_inplane_singlehole.png" width="30%" />
  <img src="../../assets/examples/Ay_inplane_singlehole.png" width="30%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Bz_surface_singlehole.png" width="30%" />
  <img src="../../assets/examples/Js_x_singlehole.png" width="30%" />
  <img src="../../assets/examples/Js_y_singlehole.png" width="30%" />
</p><br/>
```

As a verification step, we compute the flux threading through the hole by evaluating the
surface integral ``\int_h \mathbf{B} \cdot d\mathbf{S}`` over the hole area. The computed
flux agrees with the prescribed value to high accuracy:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/singlehole_flux_postpro.png" width="95%" />
</p><br/>
```

### Two holes with single flux loop

The figures below show the same field quantities for the two-hole geometry with opposing
flux (``+\Phi_0`` through one hole, ``-\Phi_0`` through the other):

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Amagnitude_doublehole.png" width="30%" />
  <img src="../../assets/examples/Ax_inplane_doublehole.png" width="30%" />
  <img src="../../assets/examples/Ay_inplane_doublehole.png" width="30%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/Bz_surface_doublehole.png" width="30%" />
  <img src="../../assets/examples/Js_x_doublehole.png" width="30%" />
  <img src="../../assets/examples/Js_y_doublehole.png" width="30%" />
</p><br/>
```

Again, the computed flux through each hole matches the prescribed values, confirming the
accuracy of the solution:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/doublehole_flux_postpro.png" width="95%" />
</p><br/>
```

### Two holes with separate flux loops

When independent flux excitations are configured, the solver extracts the full inductance
matrix from the stored magnetic energy:

```math
M_{ij} = \frac{\mathbf{A}_j^T K \mathbf{A}_i}{\Phi_i \Phi_j},
```

where ``K`` is the curl-curl stiffness matrix. The off-diagonal entries ``M_{ij}``
(``i \neq j``) represent the mutual inductance between loops, quantifying their magnetic
coupling. For the two-hole configuration on a shared plate, the computed inductance matrix
is:

```math
M = \begin{pmatrix}
2.8076 & 0.4654 \\
0.4654 & 2.8080
\end{pmatrix} \text{pH}
```

The nearly equal self-inductances reflect the geometric symmetry, while the mutual
inductance (about 17% of the self-inductance) indicates moderate magnetic coupling at this
hole separation.

### Two holes on separate planes

For the configuration with holes on spatially separated plates, the inductance matrix is:

```math
M = \begin{pmatrix}
2.9989 & -0.0510 \\
-0.0510 & 2.9988
\end{pmatrix} \text{pH}
```

The mutual inductance is an order of magnitude smaller than in the shared-plate case, despite
the same center-to-center separation between holes. This reduction arises because the
surface currents generated by each trapped flux are confined to their respective plates and
cannot overlap with the other hole's current distribution, greatly suppressing the inductive
coupling.
