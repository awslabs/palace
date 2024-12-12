```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Boundary Conditions

## Perfect electric conductor (PEC) boundary

The perfect electric conductor (PEC) boundary condition (zero tangential electric field) is
specified using the `"PEC"` boundary keyword under
[`config["Boundaries"]`](../config/boundaries.md#boundaries%5B%22PEC%22%5D). It is a
homogeneous Dirichlet boundary condition for the frequency or time domain finite element
formulation, as well as the magnetostatic formulation.

For electrostatic simulations, the homogeneous Dirichlet boundary condition is prescribed
using the [`"Ground"`](../config/boundaries.md#boundaries%5B%22Ground%22%5D) boundary
keyword which prescribes zero voltage at the boundary.

## Perfect magnetic conductor (PMC) boundary

The perfect magnetic conductor (PMC) boundary condition (zero tangential magnetic field) is
a homogeneous Neumann boundary condition for the frequency or time domain finite element
formulation, as well as the magnetostatic formulation. It is the natural boundary condition
and thus it has the same effect as not specifying any additional boundary condition on
external boundary surfaces. It can also be explicitly specified using the `"PMC"` boundary
keyword under [`config["Boundaries"]`](../config/boundaries.md#boundaries%5B%22PMC%22%5D).

Likewise, for electrostatic simulations, the homogeneous Neumann boundary condition implies
a zero-charge boundary, and thus zero gradient of the voltage in the direction normal to the
boundary. This is specified using the `"ZeroCharge"` boundary keyword under
[`config["Boundaries"]`](../config/boundaries.md#boundaries%5B%22ZeroCharge%22%5D).

## Impedance boundary

The impedance boundary condition is a mixed (Robin) boundary condition and is available for
the frequency or time domain finite element formulations and thus for eigenmode or frequency
or time domain driven simulation types. It is specified using the
[`"Impedance"`](../config/boundaries.md#boundaries%5B%22Impedance%22%5D) boundary keyword.
The surface impedance relating the tangential electric and magnetic fields on the boundary
is computed from the parallel impedances due to the specified resistance, inductance, and
capacitance per square.

## Absorbing (scattering) boundary

Absorbing boundary conditions at farfield boundaries, also referred to as scattering
boundary conditions, can be applied using the `"Absorbing"` boundary keyword under
[`config["Boundaries"]`](../config/boundaries.md#boundaries%5B%22Absorbing%22%5D). The
first-order absorbing boundary condition is a special case of the above impedance boundary
and is available for eigenmode or frequency or time domain driven simulation types. The
second-order absorbing boundary condition is only available for frequency domain driven
simulations.

[Perfectly matched layer (PML)](https://en.wikipedia.org/wiki/Perfectly_matched_layer)
boundaries for frequency and time domain electromagnetic formulations are not yet
implemented, but are
[common](https://www.sciencedirect.com/science/article/abs/pii/S0021999112000344) in solvers
for computational electromagnetics and will be a useful addition.

## Finite conductivity boundary

A finite conductivity boundary condition can be specified using the
[`"Conductivity"`](../config/boundaries.md#boundaries%5B%22Conductivity%22%5D) boundary
keyword. This boundary condition models the effect of a boundary with non-infinite
conductivity (an imperfect conductor) for conductors with thickness much larger than the
skin depth. It is available only for frequency domain driven simulations. For more
information see the
[Other boundary conditions](../reference.md#Other-boundary-conditions) section of the
reference.

## Periodic boundary

Periodic boundary conditions on an existing mesh can be specified using the
["Periodic"](../config/boundaries.md#boundaries%5B%22Periodic%22%5D) boundary keyword. This
boundary condition enforces that the solution on the specified boundaries be exactly equal,
and requires that the surface meshes on the donor and receiver boundaries be identical up to
translation or rotation. Periodicity in *Palace* is also supported through meshes generated
incorporating periodicity as part of the meshing process.

## Lumped and wave port excitation

  - [`config["Boundaries"]["LumpedPort"]`](../config/boundaries.md#boundaries%5B%22LumpedPort%22%5D) :
    A lumped port applies a similar boundary condition to a
    [surface impedance](#Impedance-boundary) boundary, but takes on a special meaning for
    each simulation type.
    
    For frequency domain driven simulations, ports are used to provide a lumped port
    excitation and postprocess voltages, currents, and scattering parameters. Likewise, for
    transient simulations, they perform a similar purpose but for time domain computed
    quantities.
    
    For eigenmode simulations where there is no excitation, lumped ports are used to specify
    properties and postprocess energy-participation ratios (EPRs) corresponding to
    linearized circuit elements.
    
    Note that a single lumped port (given by a single integer `"Index"`) can be made up of
    multiple boundary attributes in the mesh in order to model, for example, a multielement
    lumped port. To use this functionality, use the `"Elements"` object under
    [`"LumpedPort"`](../config/boundaries.md#boundaries%5B%22LumpedPort%22%5D).

  - [`config["Boundaries"]["WavePort"]`](../config/boundaries.md#boundaries%5B%22WavePort%22%5D) :
    Numeric wave ports are available for frequency domain driven simulations. In this case,
    a port boundary condition is applied with an optional excitation using a modal field
    shape which is computed by solving a 2D boundary mode eigenproblem on each wave port
    boundary. This allows for more accurate scattering parameter calculations when modeling
    waveguides or transmission lines with arbitrary cross sections.
    
    The homogeneous Dirichlet boundary conditions for the wave port boundary mode analysis
    are taken from the `"PEC"` boundaries of the full 3D model, as well as any optional
    additional boundary attributes given under `"WavePortPEC"`. Any boundary of the wave
    port not labeled with with a PEC condition has the natural boundary condition for zero
    tangential magnetic field prescribed for the purpose of port mode calculation.
    
    Unlike lumped ports, wave port boundaries cannot be defined internal to the
    computational domain and instead must exist only on the outer boundary of the domain
    (they are to be "one-sided" in the sense that mesh elements only exist on one side of
    the boundary).
    
    Wave ports are not currently compatible with nonconformal mesh refinement.

The incident field excitation at a lumped or wave port is controlled by setting
[`config["Boundaries"]["LumpedPort"][]["Excitation"]: true`](../config/boundaries.md#boundaries%5B%22LumpedPort%22%5D)
or
[`config["WavePort"][]["Excitation"]: true`](../config/boundaries.md#boundaries%5B%22WavePort%22%5D)
for that port. The excitation for each port is defined to have unit incident power over the
port boundary surface.

## Surface current excitation

An alternative source excitation to lumped or wave ports for frequency and time domain
driven simulations is a surface current excitation, specified under
[`config["Boundaries"]["SurfaceCurrent"]`](../config/boundaries.md#boundaries%5B%22SurfaceCurrent%22%5D).
This is the excitation used for magnetostatic simulation types as well. This option
prescribes a unit source surface current excitation on the given boundary in order to
excite the model. It does does not prescribe any boundary condition to the model and only
affects the source term on the right hand side.
