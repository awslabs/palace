```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Shape and Material Sensitivities

*Palace* can compute the derivatives of these quantities, for these
[problem types](../guide/problem.md), with respect to design parameters of the geometry,
material properties and lumped port properties:

| Problem type    | Quantities                                    |
|:--------------- |:--------------------------------------------- |
| `Electrostatic` | the Maxwell capacitance matrix                |
| `Magnetostatic` | the inductance matrix                         |
| `Eigenmode`     | the frequency and quality factor of each mode |

A design parameter ``p`` moves the nodes of the input mesh with velocities
``\mathbf{V}_a``, given in a table, for example as computed by the layout tool which makes
the mesh. These define the velocity field ``\mathbf{V}`` as the nodal interpolant on the
input mesh, which *Palace* evaluates at the nodes ``\mathbf{x}_k`` of the mesh it solves
on: these move with ``\mathbf{V}(\mathbf{x}_k)`` when the input mesh is morphed, through
refinement and the curving of the elements. A quantity ``J`` then changes by

```math
\frac{dJ}{dp} = \sum_k \mathbf{g}_k \cdot \mathbf{V}(\mathbf{x}_k)
```

where ``\mathbf{g}_k`` is its gradient with respect to the position of node ``k``. The
derivatives are those of the discrete quantities *Palace* computes on morphed input meshes,
so they agree with finite differences of simulations on them, and none needs an adjoint
solve. Their cost is the evaluation of the velocities at the mesh nodes, a point search in
the input mesh on one process, and a few operator applications per quantity and
parameter, plus one correction of the eigenvector of each mode of an eigenmode simulation
with losses.

  - The capacitances are ``C_{ij} = V_j^T K V_i`` for the terminal potentials ``V_i``. As
    each ``V_i`` minimizes the field energy for its terminal voltages, the derivative is
    that of the stiffness form at fixed potentials.
  - The inductances are ``M_{ij} = A_j^T K A_i / (I_i I_j)`` for the vector potentials
    ``A_i`` of the surface current sources, which solve ``K A_i = b_i``, so
    ``dM_{ij} = (b_i'(A_j) + b_j'(A_i) - K'(A_i, A_j)) / (I_i I_j)``.
  - An eigenmode solves ``T(\lambda) x = 0`` with ``\lambda = i\omega``, where ``T`` sums
    the stiffness, damping, mass and frequency-dependent terms,
    ``T(\lambda) = K + \lambda C + \lambda^2 M + A_2(\lambda)``. As ``T`` is complex
    symmetric, a simple eigenvalue changes by
    ``d\lambda = -x^T dT\, x / x^T T'(\lambda) x``, which gives the derivatives of
    ``f = \mathrm{Re}\{\omega\} / 2\pi`` and ``Q = |\omega| / 2|\mathrm{Im}\{\omega\}|``.

The resistance, inductance and capacitance of the lumped ports, and the currents of the
surface current sources, are fixed. The impedances per square of the elements of a port
scale with their width-to-length ratio, and the current density of an element of a source
with the inverse of its width, both of which *Palace* differentiates. For example, a
dilation by ``1 + t`` scales the inductances by ``1 + t``.

## Configuration

The [`config["Sensitivity"]`](../config/reference.md#config-sensitivity) section requests
the derivatives with respect to the material and lumped port properties, and to the design
parameters of an optional table of velocities:

```json
"Sensitivity":
{
  "Velocities": "velocities.csv"
}
```

The table has the columns `x`, `y` and `z` of the coordinates of nodes of the input mesh,
in its length units, then the columns `<name>.x`, `<name>.y` and `<name>.z` of the velocity
of each parameter, in mesh length units per unit of the parameter:

```
x, y, z, gap.x, gap.y, gap.z, radius.x, radius.y, radius.z
0.5, 0.25, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0
...
```

Its rows must be nodes of the mesh file, with the full precision of their coordinates,
including the high-order nodes of a curved mesh. The velocity is zero at the nodes it does
not list, and it may be compressed with gzip. With
[`config["Model"]["RemoveCurvature"]`](../config/reference.md#config-model),
the mesh moves with its vertices, so the velocities of the other nodes are ignored. The
table is read when the simulation starts, and its evaluation requires *Palace* built with
GSLIB.

The derivatives with respect to the properties of each material are those with respect to
relative changes, such as ``\varepsilon \, \partial J / \partial \varepsilon``, of its
`"Permittivity"` (electrostatic and eigenmode simulations), its `"Permeability"` and
`"LondonDepth"` (magnetostatic and eigenmode simulations), and its `"LossTan"` and
`"Conductivity"` (eigenmode simulations), where nonzero. Eigenmode simulations also give
those with respect to relative changes of the `"R"`, `"L"` and `"C"` of each lumped port.

## Output

The values and derivatives are written to `sensitivity.json` in
[`config["Problem"]["Output"]`](../config/reference.md#config-problem-output), with the
derivatives of each quantity keyed by the name of the parameter, or the material or lumped
port property:

```json
{
  "Capacitance":
  {
    "Unit": "F",
    "Terminals": [1, 2],
    "Value": [[1.2e-12, -3.4e-13], [-3.4e-13, 1.1e-12]],
    "Derivatives":
    {
      "gap": [[-2.1e-14, 5.6e-15], [5.6e-15, -1.9e-14]],
      "Materials[1].Permittivity": [[1.1e-12, -3.2e-13], [-3.2e-13, 1.0e-12]]
    }
  }
}
```

  - Electrostatic simulations write the `"Capacitance"` matrix and magnetostatic ones the
    `"Inductance"` matrix, in F and H, with the indices of their `"Terminals"` or
    `"SurfaceCurrent"` boundaries.
  - Eigenmode simulations write the `"Frequency"`, in GHz, and the `"QualityFactor"` of each
    of the `"Modes"` of `eig.csv`. The quality factor of a mode without losses, and the
    derivatives of a degenerate mode, are `null`.

The derivatives with respect to a parameter are in the unit of the quantity per unit of
the parameter, and those with respect to a property per relative change of it.

The configurations
[`spheres_sensitivity.json`](https://github.com/awslabs/palace/blob/main/examples/spheres/spheres_sensitivity.json),
[`rings_sensitivity.json`](https://github.com/awslabs/palace/blob/main/examples/rings/rings_sensitivity.json)
and
[`cavity_sensitivity.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_sensitivity.json)
of the examples read a table `velocities.csv`. For example, for the
[spheres](../examples/spheres.md), the derivatives with respect to the radius of sphere A,
of radius 1 cm and centered at ``\mathbf{c}``, follow from a radial velocity which vanishes
1 cm away from it:

```julia
using Gmsh: gmsh
using LinearAlgebra
gmsh.initialize()
gmsh.open("mesh/spheres.msh")
_, x, _ = gmsh.model.mesh.getNodes(3, -1, true)  # The nodes of the volume elements
c = [-2.5, 0.0, 0.0]
open("velocities.csv", "w") do io
    println(io, "x, y, z, radius.x, radius.y, radius.z")
    for xk in eachcol(reshape(x, 3, :))
        r = norm(xk - c)
        v = clamp(2.0 - r, 0.0, 1.0) * (xk - c) / r  # cm per cm
        println(io, join(repr.([xk; v]), ", "))
    end
end
gmsh.finalize()
```

## Limitations

  - Velocities are not available with periodic boundaries or coaxial surface current
    sources, eigenmode sensitivities not with wave ports, Floquet ports, a Floquet wave
    vector, PML regions or frequency-dependent permittivities, and magnetostatic ones not
    with flux loops, superconductors or shorted inactive ports.
  - With domain conductivity, eigenmode sensitivities require the divergence-free
    projection to be off (`config["Solver"]["Linear"]["DivFreeMaxIts"]` of 0), as it does
    not account for it.
  - The derivatives of an eigenvalue are defined only if it is simple: for a degenerate
    mode, *Palace* warns and writes `null`.
  - A velocity must keep each lumped port and surface current source planar, and a source
    with its edges off conductors along its current, for its excitation to remain
    compatible.
