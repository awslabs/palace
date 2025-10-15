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
and eigenmode simulations.

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
skin depth. It is available only for frequency domain driven and eigenmode simulations. For more
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

*Palace* also supports Floquet periodic boundary conditions, where a phase shift is imposed
between the fields on the donor and receiver boundaries. The phase shift is
``e^{-i \bm{k}_p \cdot (\bm{x}_{\textrm{receiver}}-\bm{x}_{\textrm{donor}})}``, where
``\bm{k}_p`` is the Floquet wave vector and ``\bm{x}`` is the position vector. See
[Floquet periodic boundary conditions](../reference.md#Floquet-periodic-boundary-conditions)
for implementation details.

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
    Numeric wave ports are available for frequency domain driven and eigenmode simulations. In this case,
    a port boundary condition is applied with an optional excitation using a modal field
    shape which is computed by solving a 2D boundary mode eigenproblem on each wave port
    boundary. This allows for more accurate scattering parameter calculations when modeling
    waveguides or transmission lines with arbitrary cross sections.
    
    The 2D wave port eigenproblem only supports PEC and PMC boundary conditions. Boundaries
    that are specified as `"PEC"` or `"Conductivity"` in the full 3D model and intersect the wave port
    boundary will be considered as PEC in the 2D boundary mode analysis, as well as any additional
    boundary attributes given under `"WavePortPEC"`. [`config["Boundaries"]["WavePortPEC"`](../config/boundaries.md#boundaries%5B%22WavePortPEC%22%5D)
    allows to assign non-PEC attributes from the 3D model (e.g. impedance or absorbing boundary conditions)
    as a PEC boundary condition for the 2D wave port solve. In addition, boundaries of wave ports other
    than the wave port currently being considered, in the case wave ports are touching and share one or
    more edges, are also considered as PEC for the wave port boundary mode analysis. Boundaries of the
    wave port not labeled with a `"PEC"`, `"Conductivity"`, `"WavePortPEC"`, or `"WavePort"` condition
    have the natural boundary condition of zero tangential magnetic field (PMC) prescribed for the purpose
    of port mode calculation.
    
    Unlike lumped ports, wave port boundaries cannot be defined internal to the
    computational domain and instead must exist only on the outer boundary of the domain
    (they are to be "one-sided" in the sense that mesh elements only exist on one side of
    the boundary).

For each port, the excitation is normalized to have unit incident power over the port boundary
surface.

The presence of an incident excitation at a port is controlled by the settings
[`config["Boundaries"]["LumpedPort"][]["Excitation"]`](../config/boundaries.md#boundaries%5B%22LumpedPort%22%5D)
and [`config["WavePort"][]["Excitation"]`](../config/boundaries.md#boundaries%5B%22WavePort%22%5D).
The `Excitation` settings can either be specified as non-negative integers or booleans.

  - *Boolean setting*: `true`/`false` indicates the presence / absence of an incident excitation.
    Usually, only a single port will be marked as excited. In that case, the `"Excitation"` will promoted to the port `"Index"`. If there are multiple excited ports, the `"Excitation"` is `1`.

  - *Integer setting*: Here the user manually assigns excitation indices to ports. The value `0`
    corresponds to no excitation. A positive integer `i` means that port is excited during
    excitation `i`. If multiple ports share an excitation index `i`, they will be excited at the
    same time. In the special, but common, case that each excitation consists of only a single port,
    the port index and excitation index must be equal. This avoids ambiguity in the scattering
    matrix.

For frequency domain driven simulations only, it is possible to specify multiple excitations in the
same simulation using different positive integers ("multi-excitation"). These excitations are
simulated consecutively during the Palace run. The results are printed to shared csv files. When
there are multiple excitations, the columns of the csv files are post-indexed by the excitation
index (e.g. `Φ_elec[1][5] (C)` denoting the flux through surface 1 of excitation 5). Note that a
port can only be part of one excitation.

!!! warning "Indexing"
    
    Any `"Index"` of [`"LumpedPort"`](../config/boundaries.md#boundaries%5B%22LumpedPort%22%5D),
    [`"WavePort"`](../config/boundaries.md#boundaries%5B%22WavePort%22%5D),
    [`"SurfaceCurrent"`](../config/boundaries.md#boundaries%5B%22SurfaceCurrent%22%5D), or
    [`"Terminal"`](../config/boundaries.md#boundaries%5B%22Terminal%22%5D) must be unique, including between
    different boundary conditions types (e.g. you can not have an lumped port and wave port both with
    `Index: 5`).

## Surface current excitation

An alternative source excitation to lumped or wave ports for frequency and time domain
driven simulations is a surface current excitation, specified under
[`config["Boundaries"]["SurfaceCurrent"]`](../config/boundaries.md#boundaries%5B%22SurfaceCurrent%22%5D).
This is the excitation used for magnetostatic simulation types as well. This option
prescribes a unit source surface current excitation on the given boundary in order to
excite the model. It does does not prescribe any boundary condition to the model and only
affects the source term on the right hand side.
