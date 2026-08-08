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
from driven problems in the time domain.

## Material models

Material properties are handled by the
[`config["Domains"]["Materials"]`](../config/reference.md#config-domains-materials)
object. *Palace* supports linear material models; scalar pole-residue permittivity is also
available for supported frequency-domain simulations.

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

### Scalar pole-residue permittivity

For driven frequency-domain and three-dimensional eigenmode simulations, a scalar material
can additionally specify `PermittivityPoles`. This uses the causal pole-residue model

```math
\epsilon_r(s) = \epsilon_\infty + \sum_j \frac{r_j}{s-p_j}, \qquad s=i\omega,
```

where `Permittivity` is ``\epsilon_\infty``, and `Pole` and `Residue` are SI angular rates
in rad/s. Each value is either a real number or `[real, imag]`. Real poles have real
residues; a complex pole must be supplied in the upper half-plane and *Palace* adds its
conjugate term automatically. Poles require scalar/isotropic `Permittivity` and cannot be
combined with `LossTan`; ordinary `Conductivity` remains additive.

The high-frequency term uses the ordinary permittivity mass matrix. Exact zero poles are
folded into the ordinary conductivity coefficient, while each material support with a
nonzero pole reuses one unit vector-mass operator whose complete coefficient
``\sum_j r_j s^2/(s-p_j)`` is evaluated in the frequency-dependent system term. No
pole-derived term enters the stiffness matrix.

Transient, electrostatic, magnetostatic, and boundary-mode problems, Absorbing boundaries,
Floquet or numeric wave ports, a nonzero `Periodic.FloquetWaveVector`, and adaptive circuit
synthesis are not supported with this material model. The eigenmode divergence-free
projection is disabled because longitudinal plasma modes can be physical.

Frequency-aware postprocessing is not yet available. Domain electric energy, EPR, and
quality factor; permittivity-based error indicators; electric `SurfaceFlux`; and automatic
electric boundary charge ``Q_s`` all use ``\epsilon_\infty``, not ``\epsilon(s)``, for pole
materials. Do not use these outputs as physical validation of a dispersive model.
