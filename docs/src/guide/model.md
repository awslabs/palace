```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Simulation Models

## Supported mesh formats

The [`config["Model"]`](../config/model.md#config%5B%22Model%22%5D) object is used to
specify the mesh for the discretized computational domain. In general, inputs are expected
to be dimensional nondimensionalized internally. A length scale, specified under[`config ["Model"]["L0"]`](../config/model.md#config%5B%22Model%22%5D), describes the length units
of the mesh relative to 1 meter (i.e. `config["Model"]["L0"]: 1.0e-6` if the mesh
coordinates are in ``\mu``m, this is the default value). All other entries in the
configuration file which have units of length should be specified in units of `config ["Model"]["L0"]` m.

MFEM supports a [wide variety](https://mfem.org/mesh-formats/) of mesh formats. In
addition, *Palace* has built-in support for [Nastran (`.nas`, `.bdf`)]
(https://docs.plm.automation.siemens.com/tdoc/scnastran/2020_1/help/#uid:index_element) and
[COMSOL (`.mphtxt`, `.mphbin`)]
(https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ProgrammingReferenceManual.pdf)
format mesh files, for both linear and high-order curved elements.

Geometric attributes for domains and boundaries in the mesh are used to define material
properties and boundary conditions on the desired model regions and surfaces (see
[`config["Domains"]`](../config/domains.md) and [`config["Boundaries"]`]
(../config/boundaries.md)). These attribute integers correspond to tags for the domain and
boundary elements in the mesh, and should be non-negative and start at 1. They do not need
to be contiguous in the mesh file. Throughout the configuration file, the `"Attributes"`
keyword is used to indicate which domain or boundary attributes are relevant to the
material properties or boundary conditions being specified.

## Mesh refinement

Refinement of the input mesh file can be performed using levels of global uniform refinement
or region-based refinement, specified using the [`config["Model"]["Refinement"]`]
(../config/model.md#model["Refinement"]) object. The user can specify any combination of
uniform refinement levels as well as local refinement regions which refines the elements
inside of a certain box or sphere-shaped region. For simplex meshes, the refinement
maintains a conforming mesh but meshes containing hexahedra, prism, or pyramid elements
will be non-conforming after local refinement (this is not supported at this time).

[Adaptive mesh refinement (AMR)](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement)
according to error estimates in the computed solution is a work in progress for all
simulation types.

## Material models

Material properties are handled by the [`config["Domains"]["Materials"]`]
(../config/domains.md#domains["Materials"]) object. *Palace* supports linear, frequency
independent constitutive laws for material modeling.

Materials with scalar or general matrix-valued properties are supported. For most simulation
types, each material in the model requires a specified relative permittivity and relative
permeability (for electrostatic simulations, only the permittivity is required, while for
magnetostatics, only the permeability is required). For dielectric domains, a loss tangent
may be specified. Alternatively, for normal conducting domains, an electrical conductivity
may be specified which is used to relate the current density and electric field via Ohm's
law.

Modeling of superconducting domains is performed using the current-field constitutive
relations given by the London equations. The user can specify a London penetration depth to
activate this model. It can also be used in conjunction with a materal conductivity when
wishing to model both superconducting and normal currents.

## Boundary conditions

### Lumped and wave ports

  - [`config["Boundaries"]["LumpedPort"]`]
    (../config/boundaries.md#boundaries["LumpedPort"]) :  Applies a similar boundary
    condition to a [surface impedance]
    (../config/boundaries.md#boundaries%5B%22Impedance%22%5D) boundary, but takes on a
    special meaning for each simulation type. For frequency domain driven simulations,
    ports are used to provide a lumped port excitation and postprocess voltages, currents,
    and scattering parameters. For transient simulations, they perform a similar purpose
    but for time domain computed quantities. For eigenmode simulations, ports are used to
    postprocess energy-participation ratios (EPRs) corresponding to linearized circuit
    elements.

  - [`config["Boundaries"]["WavePort"]`]
    (../config/boundaries.md#boundaries["WavePort"]) :  Numeric wave ports are available for
    frequency domain driven simulations, and apply a port boundary condition and optional
    excitation using a modal field shape which is computed by solving a 2D boundary mode
    eigenproblem on each wave port boundary. This allows for more accurate scattering
    parameter calculations when modeling waveguides or transmission lines with arbitrary
    cross sections. The homogenous Dirichlet boundary conditions for the wave port boundary
    mode analysis are taken from the `"PEC"` boundaries of the full 3D model, as well as
    any optional additional boundary attributes given under `"WavePortPEC"`. Any boundary
    of the wave port not labeled with with an essential PEC condition has the natural
    boundary condition for zero tangential magnetic field prescribe for the purpose of port
    mode calculation.

The incident field excitation at a port is controlled by setting
[`config["Boundaries"]["LumpedPort"][]["Excitation"]: true`]
(../config/boundaries.md#boundaries["LumpedPort"]) or
[`config["WavePort"][]["Excitation"]: true`]
(../config/boundaries.md#boundaries%5B%22WavePort%22%5D) for that port. The excitation for
each port is defined to have unit incident power over the port boundary surface. Finally, a
single lumped port (given by a single integer `"Index"`) can be made up of multiple boundary
attributes in the mesh to model, for example, a multielement lumped port.

### Natural boundary conditions

The perfect magnetic conductor (PMC) boundary condition (zero tangential magnetic field) is
the natural boundary condition for the frequency or time domain finite element formulation,
as well as the magnetostatic formulation. Thus, it has the same effect as not specifying any
additional boundary condition on external boundary surfaces. They can also be explicitly
specified using the `"PMC"` boundary keyword under
[`config["Boundaries"]`](../config/boundaries.md#boundaries%5B%22PMC%22%5D).

Likewise, for electrostatic simulations, the natrual boundary condition implies a
zero-charge boundary, and thus zero gradient of the voltage in the direction normal to the
boundary. This is specified using the `"ZeroCharge"` boundary keyword under
[`config["Boundaries"]`](../config/boundaries.md#boundaries%5B%22ZeroCharge%22%5D).
