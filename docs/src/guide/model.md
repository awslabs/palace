```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Simulation Models

## Supported mesh formats

The [`config["Model"]`](../config/reference.md#config-model) object is used to
specify the mesh for the discretized computational domain. In general, inputs are expected
to be dimensional nondimensionalized internally. A length scale, specified under
[`config["Model"]["L0"]`](../config/reference.md#config-model), describes the length
units of the mesh relative to 1 meter (i.e. `config["Model"]["L0"]: 1.0e-6` if the mesh
coordinates are in ``\mu``m, this is the default value). All other entries in the
configuration file which have units of length should be specified in units of
`config["Model"]["L0"]` m.

MFEM supports a [wide variety](https://mfem.org/mesh-formats/) of mesh formats. In
addition, *Palace* has built-in support for
[Nastran (`.nas`, `.bdf`)](https://docs.plm.automation.siemens.com/tdoc/scnastran/2020_1/help/#uid:index_element)
and
[COMSOL (`.mphtxt`, `.mphbin`)](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ProgrammingReferenceManual.pdf)
format mesh files, for both linear and high-order curved elements.

Geometric attributes for domains and boundaries in the mesh are used to define material
properties and boundary conditions on the desired model regions and surfaces (see
[`config["Domains"]`](../config/reference.md#config-domains) and
[`config["Boundaries"]`](../config/reference.md#config-boundaries)). These attribute integers correspond to
tags for the domain and boundary elements in the mesh, and should be non-negative and start
at 1. They do not need to be contiguous in the mesh file. Throughout the configuration
file, the `"Attributes"` keyword is used to indicate which domain or boundary attributes
are relevant to the material properties or boundary conditions being specified.

## Mesh quality

Degenerate tetrahedra slow down the iterative linear solvers: slivers, whose four vertices
are nearly coplanar (minimum dihedral angles of a few degrees or less), and needles, which
are much longer than they are wide. The basis functions of such an element are close to
linearly dependent, the smoothers of geometric multigrid cannot correct its local modes, and
every linear solve needs more iterations. A small fraction of such elements in a large
mesh can cause a large part of the iterations of a linear solve, and even a few slivers
can matter.

After the mesh information of the initial mesh and of every refined mesh, *Palace* reports
the tetrahedra with a minimum dihedral angle below 5 degrees, if there are any, with the
numbers below 2 and 1 degrees and the smallest angle in degrees:

```
 Tetrahedra with a minimum dihedral angle below 5 degrees: 1520 (214 below 2, 37 below 1, smallest 0.412)
```

This line is informational: the simulation proceeds, and for driven and eigenmode
simulations the patch solves described below remove much of the cost of the slivers.
Removing slivers and needles in the mesh generator helps every simulation type and solver.

### Optimizing meshes generated with Gmsh

[Gmsh](https://gmsh.info) can improve the quality of tetrahedral meshes after generating
them. In the Gmsh script of a model (the options and functions are the same in the Python
and Julia APIs):

  - Keep the default optimization of tetrahedral meshes (`Mesh.Optimize`, on by default),
    and add Netgen's optimizer with `gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)`
    before `gmsh.model.mesh.generate(3)`, or with `gmsh.model.mesh.optimize("Netgen")`
    after it. `Mesh.OptimizeThreshold` (default `0.3`) is the quality below which
    tetrahedra are optimized: a larger value optimizes more of them.
  - For high-order (curved) meshes, also optimize after `gmsh.model.mesh.setOrder(p)`, for
    example with `gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)` (elastic smoothing
    and optimization) or `gmsh.model.mesh.optimize("HighOrder")`, which untangle invalid
    curved elements and improve distorted ones.
  - Optimization cannot remove the slivers which the geometry forces, such as in thin
    layers, between nearly coincident surfaces, or at small features next to large
    elements. Remove these in the geometry or with mesh size fields, and check the line
    above again for the new mesh.

### Sliver patch solves

For the complex-valued operators of frequency domain driven and eigenmode simulations,
geometric multigrid also solves exactly, with dense factorizations, on small patches of the
degrees of freedom around the tetrahedra with a minimum dihedral angle below
[`config["Solver"]["Linear"]["MGSliverAngle"]`](../config/reference.md#config-solver-linear-mgsliverangle)
(2 degrees by default), on every level above the coarsest. This removes most of the extra
iterations of the slivers at a small cost. The setup reports the patches of each level, for
example:

```
  Level 1 (p = 2): 152 patches with 38410 true dofs (largest 1124, 2961 on other processes) for 214 slivers
   Patch layers: 1 (patch work ratio 0.62)
```

A level gets no patch solves when they would be expensive compared with its smoothing
(see
[`config["Solver"]["Linear"]["MGSliverLayers"]`](../config/reference.md#config-solver-linear-mgsliverlayers)),
and the patch solves are off for the real-valued operators of electrostatic, magnetostatic
and transient simulations, whose iterations they do not reduce. They also do not treat
needles.

## Mesh refinement

Refinement of the input mesh file can be performed using levels of global uniform refinement
or region-based refinement, specified using the
[`config["Model"]["Refinement"]`](../config/reference.md#config-model-refinement) object.
The user can specify any combination of uniform refinement levels as well as local
refinement regions which refines the elements inside of a certain box or sphere-shaped
region. For simplex meshes, the refinement maintains a conforming mesh but meshes
containing hexahedra, prism, or pyramid elements will be nonconforming after local
refinement (this is not supported at this time).

[Adaptive mesh refinement (AMR)](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement)
according to error estimates calculated from the computed solution can also be specified
using the [`config["Model"]["Refinement"]`](../config/reference.md#config-model-refinement)
object. Nonconformal refinement is supported for all mesh types, and additionally conformal
refinement is supported for simplex meshes. AMR is available for all problem types apart
from driven problems in the time domain. The error estimates account for the
discontinuity of the fields across [interior boundaries](boundaries.md#Interior-boundaries)
such as thin metal sheets, without modifying the mesh.

## Material models

Material properties are handled by the
[`config["Domains"]["Materials"]`](../config/reference.md#config-domains-materials)
object. *Palace* supports linear, frequency independent constitutive laws for material
modeling.

Materials with scalar or symmetric matrix-valued material properties are supported. For most
simulation types, each material in the model requires a specified relative permittivity and
relative permeability (for electrostatic simulations, only the permittivity is required,
while for magnetostatics, only the permeability is required). For dielectric domains, a
loss tangent may be specified. Alternatively, for normal conducting domains, an electrical
conductivity may be specified which is used to relate the current density and electric
field via Ohm's law.

Modeling of superconducting domains is performed using the current-field constitutive
relations given by the London equations. The user can specify a London penetration depth to
activate this model. It can also be used in conjunction with a material conductivity when
wishing to model both superconducting and normal currents.
