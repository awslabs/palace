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

The participation ratios for bulk dielectrics and interface dielectric layers can be
computed for simulations involving the electric field. For model boundaries, the integrated
surface charge or magnetic flux can also be postprocessed. These features are described
in [Domain postprocessing](#Domain-postprocessing) and in [Boundary postprocessing]
(#Boundary-postprocessing).

Additionally, the computed fields can be automatically probed for their vector values at one
or more points in space. This probe functionality is also described in
[Domain postprocessing](#Domain-postprocessing).

Finally, as described further in [Visualization](#Visualization), various field quantities
on the 3D computational domain as well as 2D domain boundaries and material interfaces are
written to disk when requested using the relevant parameters under [`config["Solver"]`]
(../config/solver.md). These fields are meant to be visualized with [ParaView]
(https://www.paraview.org/).

## Domain postprocessing

Domain postprocessing capabilities are enabled by including objects under
[`config["Domains"]["Postprocessing"]`](../config/domains.md) in the configuration file.
These include:

  - [`config["Domains"]["Postprocessing"]["Dielectric"]`]
    (../config/domains.md#domains["Postprocessing"]["Dielectric"]) :  Postprocessess bulk
    dielectric loss based on the participation ratio of the electric field in a lossy
    region. The respective participation ratios and quality factors for each domain
    (associated with the specified domain attributes and indexed by the specified integer
    `"Index"`) are computed using the material properties provided and are written to
    `domain-Q.csv` in the specified postprocessing output directory.
  - [`config["Domains"]["Postprocessing"]["Probe"]`]
    (../config/domains.md#domains["Postprocessing"]["Probe"]) :  Probe the values of the
    computed electric field and magnetic flux density solutions at specified locations in
    the computational domain. The availability of the ``\bm{E}`` and ``\bm{B}`` fields
    depends on the problem type (for example, for magnetostatic problems, only ``\bm{B}``
    is output and ``\bm{E}`` is not computed, whereas the inverse is true for
    electrostatics). For each computed field, the postprocessed values are written to
    `probe-E.csv` and `probe-B.csv` in the specified output directory.

## Boundary postprocessing

Boundary postprocessing capabilities are enabled by including objects under
`config["Boundaries"]["Postprocessing"]` in the configuration file. These include:

  - [`config["Boundaries"]["Postprocessing"]["Capacitance"]`]
    (../config/boundaries.md#boundaries["Postprocessing"]["Capacitance"]) :  Postprocess the
    integral of the surface charge on a surface defined by a list of boundary attributes,
    and divide by the excitation voltage to get the capacitive coupling. The resulting
    capcitances are written to `surface-C.csv` in the specified output directory.
  - [`config["Boundaries"]["Postprocessing"]["Inductance"]`]
    (../config/boundaries.md#boundaries["Postprocessing"]["Inductance"]) :  Postprocess the
    magnetic flux through a surface defined by a list of boundary attributes, and divide by
    the excitation current to the inductive coupling. The resulting inductances are written
    to `surface-M.csv` in the specified output directory.
  - [`config["Boundaries"]["Postprocessing"]["Dielectric"]`]
    (../config/boundaries.md#boundaries["Postprocessing"]["Dielectric"]) :  Postprocesses
    interface dielectric loss at surfaces of the model by specifying the interface
    thickness, permittivity, and loss tangent. See [https://arxiv.org/pdf/1509.01854.pdf]
    (https://arxiv.org/pdf/1509.01854.pdf) or
    [https://aip.scitation.org/doi/10.1063/1.3637047]
    (https://aip.scitation.org/doi/10.1063/1.3637047) for more information. The
    participation ratios and associated quality factors are written to the file
    `surface-Q.csv` in the specified output directory.

## Visualization

When specified in the configuration file, the electric field and magnetic flux density
solutions are written to disk for 3D visualization with [ParaView]
(https://www.paraview.org/). Various other postprocessed fields are also written to the
ParaView database as available, including electric and magnetic energy density, surface
currents, and charge density. These files are found in the `paraview/` directory located in
the output directory specified under [`config["Problem"]["Output"]`]
(../config/problem.md#config["Problem"]).

In addition to the full 3D fields, a ParaView data collection for the boundary mesh is also
written to disk. The boundary mesh includes all surfaces with prescribed boundary
conditions as well as any material interfaces in the computational domain.
