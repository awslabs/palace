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
object. *Palace* supports linear constitutive laws for material modeling, which are
frequency independent except for the complex conductivity models of superconducting domains
described below.

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

Alternatively, a superconducting domain can be described by a complex conductivity
``\sigma(\omega) = \sigma_1(\omega) - i\sigma_2(\omega)``, specified with
[`config["Domains"]["Materials"][]["ComplexConductivity"]`](../config/reference.md#config-domains-materials):

  - `"TwoFluid"`: frequency-independent ``\sigma_1`` and London penetration depth, computed
    from the normal-state conductivity, the zero-temperature penetration depth, and the
    temperature with the classical (Gorter-Casimir) or enhanced (Vendik) temperature law.
    This is equivalent to `"Conductivity"` and `"LondonDepth"`.
  - `"MattisBardeen"`: the local (dirty-limit) Mattis-Bardeen conductivity, from the
    normal-state conductivity, the temperature, and the gap or the critical temperature
    (weak-coupling BCS gap). It captures the frequency dependence up to and above the gap
    frequency ``2\Delta/h``, as needed for picosecond pulses on superconducting lines.
  - `"Table"`: tabulated ``\sigma_1`` and ``\sigma_2``, interpolated in frequency.

The dispersive models (`"MattisBardeen"` and `"Table"`) are supported for driven
simulations, with uniform or adaptive frequency sweeps and lumped or wave ports, and for
boundary mode simulations. They are not yet supported for eigenmode or transient
simulations, or adaptive circuit synthesis. The kinetic energy and conduction loss of
superconducting domains are postprocessed in driven and eigenmode simulations (see
[Domain postprocessing](postprocessing.md#Domain-postprocessing)).
