```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Postprocessing and Visualization

As described in the section [Problem Types](problem.md), each simulation type writes
relevant postprocessed scalar quantities to disk in the directory specified by
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D), including
but not limited to computed values like eigenfrequencies, scattering parameters, or lumped
element parameters. In addition, each simulation type will write a file called
`domain-E.csv`, which includes information about the electric and magnetic field energies,
as well as lumped element energies, for each step of the simulation (eigenmode, frequency,
or time step, for examples).

Models containing lumped or wave port boundaries or surface current excitations will
automatically postprocess quantities related to those boundaries. This is described in
[Ports and surface currents](#Ports-and-surface-currents).

The participation ratios for bulk dielectrics and interface dielectric layers can be
computed for simulations involving the electric field. For model boundaries, the integrated
surface charge or magnetic flux can also be postprocessed. These features are described
in [Domain postprocessing](#Domain-postprocessing) and in
[Boundary postprocessing](#Boundary-postprocessing).

Additionally, the computed fields can be automatically probed for their vector values at one
or more points in space. This probe functionality is also described in
[Domain postprocessing](#Domain-postprocessing).

Finally, as described further in [Visualization](#Visualization), various field quantities
on the 3D computational domain as well as 2D domain boundaries and material interfaces are
written to disk when requested using the relevant parameters under
[`config["Solver"]`](../config/solver.md). These fields are meant to be visualized with
[ParaView](https://www.paraview.org/) or [GLVis](https://glvis.org/).

## Ports and surface currents

When lumped ports are present in a model, the lumped port voltages and currents computed for
each step of the simulation (eigenmode, frequency, or time step) are written to ASCII files
named `port-V.csv` and `port-I.csv`, respectively. These files also include the excitation
voltage and current corresponding to the incident wave on excited port boundaries.

Additionally, when surface current excitations are present, the excitations are written to
`surface-I.csv`.

For frequency domain problems, the values output are the complex-valued peak voltages and
currents, computed from the field phasors.

## Domain postprocessing

Domain postprocessing capabilities are enabled by including objects under
[`config["Domains"]["Postprocessing"]`](../config/domains.md) in the configuration file.
These include:

  - [`config["Domains"]["Postprocessing"]["Energy"]`](../config/domains.md#domains%5B%22Postprocessing%22%5D%5B%22Energy%22%5D) :
    Postprocessess the electric and magnetic field energy inside of a given domain
    (associated with the specified domain attributes and indexed by the specified integer
    `"Index"`). These are from the electric and magnetic field solutions and written to the
    same `domain-E.csv` file in the specified postprocessing output directory used for the
    global energies (described above).
  - [`config["Domains"]["Postprocessing"]["Probe"]`](../config/domains.md#domains%5B%22Postprocessing%22%5D%5B%22Probe%22%5D) :
    Probe the values of the computed electric field and magnetic flux density solutions at
    specified locations in the computational domain. The availability of the ``\bm{E}`` and
    ``\bm{B}`` fields depends on the problem type (for example, for magnetostatic problems,
    only ``\bm{B}`` is output and ``\bm{E}`` is not computed, whereas the inverse is true
    for electrostatics). For each computed field, the postprocessed values are written to
    `probe-E.csv` and `probe-B.csv` in the specified output directory.

## Boundary postprocessing

Boundary postprocessing capabilities are enabled by including objects under
[`config["Boundaries"]["Postprocessing"]`](../config/boundaries.md) in the configuration
file. These include:

  - [`config["Boundaries"]["Postprocessing"]["SurfaceFlux"]`](../config/boundaries.md#boundaries%5B%22Postprocessing%22%5D%5B%22SurfaceFlux%22%5D) :
    Postprocess the integrated flux through a surface defined by a list of boundary
    attributes. Electric, magnetic, and power flux are all supported. Surface capacitance
    can be computed by dividing the computed electric flux by the excitation voltage, while
    inductance can be computed by dividing the computed magnetic flux by the excitation
    current. The resulting fluxes are written to `surface-F.csv` in the specified output
    directory.
  - [`config["Boundaries"]["Postprocessing"]["Dielectric"]`](../config/boundaries.md#boundaries%5B%22Postprocessing%22%5D%5B%22Dielectric%22%5D) :
    Postprocesses interface dielectric loss at surfaces of the model by specifying the
    interface thickness, permittivity, and loss tangent. See the
    [Bulk and interface dielectric loss](../reference.md#Bulk-and-interface-dielectric-loss)
    section of the reference, or
    [https://arxiv.org/pdf/1509.01854.pdf](https://arxiv.org/pdf/1509.01854.pdf) or
    [https://aip.scitation.org/doi/10.1063/1.3637047](https://aip.scitation.org/doi/10.1063/1.3637047)
    for more information. The participation ratios and associated quality factors are
    written to the file `surface-Q.csv` in the specified output directory.

## Visualization

When specified in the configuration file, the electric field and magnetic flux density
solutions are written to disk for 3D visualization with [ParaView](https://www.paraview.org/)
or [GLVis](https://glvis.org/). Various other postprocessed fields are also written to the ParaView
or grid function (GLVis) database as available, including electric and magnetic energy density,
surface currents, and charge density. These files are found in the `paraview/` or `gridfunction/`
directories located in the output directory specified under
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D). The output
formats are specified in [`config["Problem"]["OutputFormats"]`](../config/problem.md#config%5B%22Problem%22%5D).

ParaView is recommended to visualize large simulations in parallel. The grid function (GLVis)
format can be useful to embed visualizations in webpages with its
[Javascript version](https://github.com/GLVis/glvis-js/).

All fields are written out in SI units and the post-processing mesh has the same units of `config["Model"]["L0"]` m
as the input mesh. The specific quantities available vary by [simulation type](problem.md#Problem-Types),
but the variable names and corresponding units for various possible postprocessed scalar and vector are:

  - Electric field: `E`, `E_real`, and `E_imag` (V/m)
  - Magnetic flux density: `B`, `B_real`, and `B_imag` (Wb/m²)
  - Electric potential: `V` (V)
  - Magnetic vector potential : `A`, `A_real`, and `A_imag` (A)
  - Electric energy density : `U_e` (J/m³)
  - Magnetic energy density : `U_m` (J/m³)
  - Poynting vector: `S` (W/m²)

Also, at the final step of the simulation the following element-wise quantities are written
for visualization:

  - Mesh partitioning (1-based): `Rank`
  - Error indicator: `Indicator`

When saving fields in the grid function (GLVis) format, the file names have the format
`Field_xxxxxx.gf.yyyyyy` where `Field` is the variable name of the postprocessed scalar
or vector field, `xxxxxx` is the six-digit index of the terminal index (electrostatic
or magnetostatic), time step index (transient), or frequency index (driven or eigenmode),
and `yyyyyy` is the six-digit index of the rank of the corresponding MPI process.

In addition to the full 3D fields, a ParaView data collection for the boundary mesh and
fields is also written to disk. The boundary mesh includes all surfaces with prescribed
boundary conditions as well as any material interfaces in the computational domain. It is
located in the same `paraview/` directory, with suffix `_boundary`. The boundary data
collection is only available for the ParaView output format.

The boundary data collection includes the 3D field values sampled on the boundary mesh as
well as:

  - Surface charge density: `Q_s`, `Q_s_real`, `Q_s_imag` (Wb/m²)
  - Surface current density: `J_s`, `J_s_real`, `J_s_imag` (A/m)
  - Wave port boundary mode electric field: `E0_real`, `E0_imag` (V/m)

## Adaptive mesh refinement

At the start of an adaptive mesh refinement (AMR) iteration, if
[`config["Model"]["Refinement"]["SaveAdaptIterations"]`](../config/model.md#model%5B%22Refinement%22%5D)
is enabled, the postprocessing results from the solve on the previous mesh will be saved off
within a subdirectory denoted `iterationX`, where `X` is the (1-based) iteration number.
The results in the top level directory will always be those from the most recent successful
solve.
