```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Shape and Material Sensitivities

*Palace* can compute the derivatives of objectives with respect to the positions of the mesh
nodes in a region, and with respect to material properties:

| [Problem type](../guide/problem.md) | Objectives (`"Type"`)                                                    |
|:----------------------------------- |:------------------------------------------------------------------------ |
| `Electrostatic`                     | `"Capacitance"`: an entry of the Maxwell capacitance matrix              |
| `Eigenmode`                         | `"Frequency"` and `"QualityFactor"`: the frequency and Q of an eigenmode |

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
  - An eigenmode solves ``T(\lambda) x = 0`` with ``T(\lambda) = K + \lambda C + \lambda^2 M + A_2(\lambda)`` and ``\lambda = i\omega``. As ``T`` is complex
    symmetric, a simple eigenvalue changes by ``d\lambda = -x^T dT\, x / x^T T'(\lambda) x``, which gives the derivatives of ``f = \mathrm{Re}\{\omega\} / 2\pi``
    and ``Q = |\omega| / 2|\mathrm{Im}\{\omega\}|``.

## Configuration

The [`config["Sensitivity"]`](../config/reference.md#config-sensitivity) section lists the
objectives: capacitances by the terminal indices of their entries, and eigenmode objectives
by the index of their mode in `eig.csv`, for example
`{ "Index": 1, "Type": "Frequency", "Mode": 1 }`:

```json
"Sensitivity":
{
  "Objectives":
  [
    { "Index": 1, "Type": "Capacitance", "Terminals": [1, 1] },
    { "Index": 2, "Type": "Capacitance", "Terminals": [1, 2] }
  ],
  "Shape":
  {
    "BoundingBoxMin": [-5.0, -5.0, -5.0],
    "BoundingBoxMax": [5.0, 5.0, 5.0],
    "BoundaryTags": true
  },
  "MaterialProperties": true
}
```

`"Shape"` requests the gradients with respect to the node positions, for the nodes in a box
(in mesh length units) and of the elements with the given domain `"Attributes"`, or of the
whole mesh. The region must contain every node a velocity moves. `"MaterialProperties"`
requests the derivative with respect to a relative change of the permittivity of each
material, ``\varepsilon \, \partial C_{ij} / \partial \varepsilon``.

## Output

The derivatives are written to the `sensitivity` folder of
[`config["Problem"]["Output"]`](../config/reference.md#config-problem-output):

  - `sensitivity.json`, with each objective's value, unit and material derivatives;
  - `points.bin`, the positions of the exported nodes in mesh length units;
  - `gradient-<index>.bin`, each objective's gradient, in its unit (F, GHz, or none for Q)
    per mesh length unit;
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

  - Node gradients are not available with periodic boundaries, and eigenmode sensitivities
    not with wave ports, Floquet ports or a Floquet wave vector.
  - The derivatives of an eigenvalue are defined only if it is simple. For a degenerate
    mode, *Palace* warns, and the result is for the eigenvector it computed. The quality
    factor of a lossless mode is not defined.
  - The gradients are at fixed impedances per square of the lumped ports, so velocities
    should be rigid on lumped ports, which `check_velocity` can check.
  - The nodes are those of the mesh *Palace* solves on. Where *Palace* refines the input
    mesh or, for high-order meshes read from Gmsh, re-interpolates its nodes, morphing the
    input mesh by ``\mathbf{V}`` moves these nodes by the interpolant of ``\mathbf{V}``:
    the derivatives then agree with finite differences of simulations on morphed input
    meshes to the interpolation error of ``\mathbf{V}``.
