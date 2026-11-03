```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Shape and Material Sensitivities

*Palace* can compute the derivatives of objectives with respect to the positions of the mesh
nodes in a region, and with respect to material properties, for these
[problem types](../guide/problem.md):

| Problem type    | Objectives (`"Type"`)                                           |
|:--------------- |:--------------------------------------------------------------- |
| `Electrostatic` | `"Capacitance"`: an entry of the Maxwell capacitance matrix     |
| `Magnetostatic` | `"Inductance"`: an entry of the inductance matrix               |
| `Eigenmode`     | `"Frequency"`, `"QualityFactor"`: the frequency and Q of a mode |

The node gradients give the derivative with respect to any design parameter of the
geometry: if a parameter ``p`` moves the geometry with velocity ``\mathbf{V}(\mathbf{x})``,
for example as computed by a layout tool, an objective ``J`` changes by

```math
\frac{dJ}{dp} = \sum_k \mathbf{g}_k \cdot \mathbf{V}(\mathbf{x}_k)
```

where ``\mathbf{g}_k`` is the gradient of ``J`` with respect to the position
``\mathbf{x}_k`` of node ``k``. *Palace* needs no knowledge of the parametrization, and one
simulation gives the derivatives with respect to every parameter.

The derivatives are those of the discrete objectives *Palace* computes, so they agree with
finite differences of simulations on moved meshes, and none needs an additional solve: the
cost is a few operator applications per objective.

  - The capacitances are ``C_{ij} = V_j^T K V_i`` for the terminal potentials ``V_i``. As
    each ``V_i`` minimizes the field energy for its terminal voltages, the derivative is
    that of the stiffness form at fixed potentials.
  - The inductances are ``M_{ij} = A_j^T K A_i / (I_i I_j)`` for the vector potentials
    ``A_i`` of the surface current sources, which solve ``K A_i = b_i``, so
    ``dM_{ij} = (b_i'(A_j) + b_j'(A_i) - K'(A_i, A_j)) / (I_i I_j)``, with ``b'`` the
    derivative of an excitation at a fixed current density.
  - An eigenmode solves ``T(\lambda) x = 0`` with ``\lambda = i\omega``, where ``T`` sums
    the stiffness, damping, mass and frequency-dependent terms,
    ``T(\lambda) = K + \lambda C + \lambda^2 M + A_2(\lambda)``. As ``T`` is complex
    symmetric, a simple eigenvalue changes by
    ``d\lambda = -x^T dT\, x / x^T T'(\lambda) x``, which gives the derivatives of
    ``f = \mathrm{Re}\{\omega\} / 2\pi`` and ``Q = |\omega| / 2|\mathrm{Im}\{\omega\}|``.

## Configuration

The [`config["Sensitivity"]`](../config/reference.md#config-sensitivity) section lists the
objectives: capacitances and inductances by the `"Terminal"` or `"SurfaceCurrent"` indices
of their entries, and eigenmode objectives by the index of their mode in `eig.csv`:

```json
"Sensitivity":
{
  "Objectives":
  [
    { "Index": 1, "Type": "Frequency", "Mode": 1 },
    { "Index": 2, "Type": "QualityFactor", "Mode": 1 }
  ],
  "Shape":
  {
    "BoundingBoxMin": [-5.0, -5.0, -5.0],
    "BoundingBoxMax": [5.0, 5.0, 5.0],
    "BoundaryTags": true
  },
  "MaterialProperties": true,
  "LumpedPorts": true
}
```

`"Shape"` requests the gradients with respect to the node positions, for the nodes in a box
(in mesh length units) and of the elements with the given domain `"Attributes"`, or of the
whole mesh. The region must contain every node a velocity moves.

`"MaterialProperties"` requests the derivatives with respect to relative changes of the
properties of each material, such as ``\varepsilon \, \partial J / \partial \varepsilon``:
its `"Permittivity"` (electrostatic and eigenmode simulations), its `"Permeability"`
(magnetostatic and eigenmode simulations), and its `"LossTan"`, `"Conductivity"` and
`"LondonDepth"` (eigenmode simulations), where nonzero. In eigenmode simulations,
`"LumpedPorts"` requests the derivatives with respect to relative changes of the `"R"`,
`"L"` and `"C"` of each lumped port. As the impedances per square of a port scale with its
width-to-length ratio ``q``, they also give the derivative with respect to it:

```math
\frac{\partial J}{\partial \log q} = R \frac{\partial J}{\partial R} + L \frac{\partial J}{\partial L} - C \frac{\partial J}{\partial C}
```

## Output

The derivatives are written to the `sensitivity` folder of
[`config["Problem"]["Output"]`](../config/reference.md#config-problem-output):

  - `sensitivity.json`, with each objective's value, unit, and material and lumped port
    derivatives;
  - `points.bin`, the positions of the exported nodes in mesh length units;
  - `gradient-<index>.bin`, each objective's gradient, in its unit (F, H, GHz, or none for
    Q) per mesh length unit;
  - `boundary.bin`, with `"BoundaryTags"`, the boundary attributes of each node.

The Julia package
[`scripts/julia/PalaceSensitivity`](https://github.com/awslabs/palace/tree/main/scripts/julia/PalaceSensitivity)
reads this export, whose format its README specifies, and contracts it with velocities. For
example, for the [spheres](../examples/spheres.md) with
[`spheres_sensitivity.json`](https://github.com/awslabs/palace/blob/main/examples/spheres/spheres_sensitivity.json),
the derivatives with respect to the radius of sphere A, of radius 1 cm and centered at
``\mathbf{c}``, follow from a radial velocity which vanishes 1 cm away from it:

```julia
using LinearAlgebra, PalaceSensitivity
s = read_sensitivity("postpro/spheres_sensitivity/sensitivity")
c = [-2.5, 0.0, 0.0]
radial(x) = (r=norm(x - c); max(0.0, min(1.0, 2.0 - r)) * (x - c) / r)
V(X) = reduce(hcat, radial.(eachcol(X)))
dC = jacobian(s, (radius=V,); velocity_length_unit=s.length_unit)  # F/cm
```

For the cylinder cavity with
[`cavity_sensitivity.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_sensitivity.json),
of radius ``a = 2.74`` cm and height ``d = 5.48`` cm along ``z``, the derivatives of the
frequencies and the quality factor with respect to ``a`` and ``d`` follow from velocities
which scale them:

```julia
s = read_sensitivity("postpro/cavity_sensitivity/sensitivity")
radius(X) = X .* [1.0, 1.0, 0.0] ./ 2.74
height(X) = X .* [0.0, 0.0, 1.0] ./ 5.48
J = jacobian(s, (; radius, height); velocity_length_unit=s.length_unit)  # per cm
```

`check_velocity` checks that a velocity vanishes, or is rigid, on given boundary
attributes, from `boundary.bin`.

## Limitations

  - Node gradients are not available with periodic boundaries, eigenmode sensitivities not
    with wave ports, Floquet ports or a Floquet wave vector, and magnetostatic ones not with
    flux loops, superconductors or shorted inactive ports.
  - The derivatives of an eigenvalue are defined only if it is simple. For a degenerate
    mode, *Palace* warns, and the result is for the eigenvector it computed. The quality
    factor of a lossless mode is not defined.
  - The gradients are at fixed impedances per square of the lumped ports, and at fixed
    current densities of the surface current sources. A velocity which changes the width or
    length of a port must add the derivative with respect to its width-to-length ratio
    above; others should be rigid on the ports and the sources, which `check_velocity` can
    check. For example, a dilation by ``1 + t`` scales the inductances by ``(1 + t)^3`` at
    fixed current densities, rather than ``1 + t``. A velocity must also keep each source
    planar, with its edges off conductors along its current, for its excitation to remain
    compatible.
  - The nodes are those of the mesh *Palace* solves on. Where *Palace* refines the input
    mesh or, for high-order meshes read from Gmsh, re-interpolates its nodes, morphing the
    input mesh by ``\mathbf{V}`` moves these nodes by the interpolant of ``\mathbf{V}``:
    the derivatives then agree with finite differences of simulations on morphed input
    meshes to the interpolation error of ``\mathbf{V}``.
