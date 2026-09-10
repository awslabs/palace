```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Postprocessing and Visualization

As described in the section [Problem Types](problem.md), each simulation type writes
relevant postprocessed scalar quantities to disk in the directory specified by
[`config["Problem"]["Output"]`](../config/reference.md#config-problem), including
but not limited to computed values like eigenfrequencies, scattering parameters, or lumped
element parameters. In addition, each simulation type will write a file called
`domain-E.csv`, which includes information about the electric and magnetic field energies,
as well as lumped element energies, for each step of the simulation (eigenmode, frequency,
or time step, for examples).

Models containing lumped or wave port boundaries or surface current excitations will
automatically postprocess quantities related to those boundaries. This is described in
[Ports and surface currents](#Ports-and-surface-currents).

## Resolved configuration record

At startup *Palace* writes a fully-resolved copy of the run configuration to the output
directory specified by
[`config["Problem"]["Output"]`](../config/reference.md#config-problem-output), named after
the input file with a `_resolved.json` suffix (for example, `cavity.json` produces
`cavity_resolved.json`). Every implicit default and `"Default"` sentinel is replaced with the
concrete value *Palace* actually used, so the file is a complete, self-contained record of
the run. It passes schema validation and can be supplied directly to *Palace* to reproduce
the same simulation deterministically. A small number of options whose default is delegated
to an external library (for example the sparse direct solver `"ColumnOrdering"` and the GMRES
`"PCSide"`) remain `"Default"`, since *Palace* does not select a concrete value for them. The
same resolved configuration can be printed without running a simulation by passing
`--dry-run`.

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
[`config["Solver"]`](../config/reference.md#config-solver). These fields are meant to be visualized with
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
[`config["Domains"]["Postprocessing"]`](../config/reference.md#config-domains) in the configuration file.
These include:

  - [`config["Domains"]["Postprocessing"]["Energy"]`](../config/reference.md#config-domains-postprocessing-energy) :
    Postprocessess the electric and magnetic field energy inside of a given domain
    (associated with the specified domain attributes and indexed by the specified integer
    `"Index"`). These are from the electric and magnetic field solutions and written to the
    same `domain-E.csv` file in the specified postprocessing output directory used for the
    global energies (described above).
  - [`config["Domains"]["Postprocessing"]["Probe"]`](../config/reference.md#config-domains-postprocessing-probe) :
    Probe the values of the computed electric field and magnetic flux density solutions at
    specified locations in the computational domain. The availability of the ``\bm{E}`` and
    ``\bm{B}`` fields depends on the problem type (for example, for magnetostatic problems,
    only ``\bm{B}`` is output and ``\bm{E}`` is not computed, whereas the inverse is true
    for electrostatics). For each computed field, the postprocessed values are written to
    `probe-E.csv` and `probe-B.csv` in the specified output directory.

## Boundary postprocessing

Boundary postprocessing capabilities are enabled by including objects under
[`config["Boundaries"]["Postprocessing"]`](../config/reference.md#config-boundaries) in the configuration
file. These include:

  - [`config["Boundaries"]["Postprocessing"]["SurfaceFlux"]`](../config/reference.md#config-boundaries-postprocessing-surfaceflux) :
    Postprocess the integrated flux through a surface defined by a list of boundary
    attributes. Electric, magnetic, and power flux are all supported. Surface capacitance
    can be computed by dividing the computed electric flux by the excitation voltage, while
    inductance can be computed by dividing the computed magnetic flux by the excitation
    current. The resulting fluxes are written to `surface-F.csv` in the specified output
    directory.
  - [`config["Boundaries"]["Postprocessing"]["Dielectric"]`](../config/reference.md#config-boundaries-postprocessing-dielectric) :
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
[`config["Problem"]["Output"]`](../config/reference.md#config-problem). The output
formats are specified in [`config["Problem"]["OutputFormats"]`](../config/reference.md#config-problem).

*Palace* writes ParaView data using raw-appended binary VTK XML. Use
[ParaView 5.12.1](https://www.kitware.com/paraview-5-12-1-release-notes/) or newer to
read these files. Older versions may work with a compatible VTK/Expat combination,
but are not supported because some distribution builds fail to parse raw-appended data.

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

### Visualization field spaces and interface values

Domain fields available in both the ParaView and grid function formats are represented in
different ways. ParaView stores values sampled at visualization points, whereas an MFEM grid
function stores the finite element space and its degrees of freedom. Primary solution fields
retain their native finite element spaces in grid function output. Derived fields are
evaluated directly into discontinuous, interpolatory L2 spaces of the solution order:

| Fields                                           | Grid function space | Continuity represented by the space                                                                        |
|:------------------------------------------------ |:------------------- |:---------------------------------------------------------------------------------------------------------- |
| `V`, `En`, and other scalar H1 solution fields   | H1                  | The scalar value is continuous across conforming element interfaces.                                       |
| `E`, `E_real`, `E_imag`, `A`, `A_real`, `A_imag` | ND (H(curl))        | The tangential trace is continuous; the normal component may differ between the two sides of an interface. |
| `B`, `B_real`, `B_imag` in 3D                    | RT (H(div))         | The normal trace is continuous; the tangential component may differ between the two sides of an interface. |
| `B`, `B_real`, `B_imag` in 2D                    | Scalar L2           | Element-local scalar values are retained.                                                                  |
| `U_e`, `U_m`, `S`, `Sn`                          | Scalar or vector L2 | All element-local values are retained, including both traces at an interface.                              |

The L2 representation of the derived fields is important even when a primary solution field
belongs to a conforming space. H1, H(curl), and H(div) impose different notions of
continuity; neither H(curl) nor H(div) makes the complete vector field single-valued at an
interface. In addition, the constitutive parameters can be discontinuous between materials.
For example,

```math
U_e = \frac{1}{2} \bm{E}^{*} \bm{\epsilon} \bm{E}, \qquad
U_m = \frac{1}{2} \bm{B}^{*} \bm{\mu}^{-1} \bm{B}, \qquad
\bm{S} = \operatorname{Re}\!\left\{\bm{E} \times
\left(\bm{\mu}^{-1}\bm{B}\right)^{*}\right\}.
```

A jump in ``\bm{\epsilon}`` or ``\bm{\mu}^{-1}`` can therefore produce a jump in an energy
density or the Poynting vector, even where the conforming trace of a primary field is
continuous. Similarly, at a material interface:

  - H(curl) continuity constrains only the tangential electric field; its normal component
    can jump.
  - H(div) continuity constrains only the normal magnetic flux density; its tangential
    component can jump.
  - ``\bm{D} = \bm{\epsilon}\bm{E}`` and
    ``\bm{H} = \bm{\mu}^{-1}\bm{B}`` inherit jumps from both the fields and the material
    tensors.
  - Surface charge ``Q_s = \bm{D} \cdot \bm{n}`` and surface current
    ``\bm{J}_s = \bm{n} \times \bm{H}`` are side-dependent traces when the corresponding
    interface charge or current is nonzero.

Forcing these quantities into a continuous H1 space would select, average, or smooth values
across interfaces and would discard the distinction between the two element-side traces.
The L2 domain fields instead preserve the element-local result. The separate boundary
visualization has one tuple per geometric interface point: primary fields are the average of
the two traces, energy density and Poynting output average the corresponding per-side
quantities, and surface charge/current use their oriented jump definitions. Exterior
boundaries use the single adjacent element trace.

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
[`config["Model"]["Refinement"]["SaveAdaptIterations"]`](../config/reference.md#config-model-refinement)
is enabled, the postprocessing results from the solve on the previous mesh will be saved off
within a subdirectory denoted `iterationX`, where `X` is the (1-based) iteration number.
The results in the top level directory will always be those from the most recent successful
solve.

When [`config["Model"]["Refinement"]["SaveAdaptMesh"]`](../config/reference.md#config-model-refinement)
is enabled and adaptation was performed, *Palace* saves the adapted mesh to disk and also
records a `SavedAdaptedMesh` block in `palace.json` summarizing its topology. The block always
describes the mesh saved next to it: the top-level `palace.json` describes the final adapted
mesh, and (when `SaveAdaptIterations` is also enabled) each `iterationX/palace.json` describes
that iteration's saved mesh, so a saved mesh and its block can always be paired. The block
contains:

  - `Dimension` : Spatial dimension of the mesh, as an integer (`2` or `3`).
  - `TrueVertices`, `TrueEdges` : True (order-independent, conforming) counts of the mesh
    vertices and edges.
  - `TrueFaces` : True (conforming) counts of the codimension-one faces, broken down by face
    geometry (`triangle`, `quadrilateral`); only the geometries with a nonzero count are
    listed. Reported in 3D only. In 2D the codimension-one faces are edges and are already
    given by `TrueEdges`, so this field is omitted.
  - `Cells` : Counts of the mesh elements, broken down by element geometry (`tetrahedron`,
    `hexahedron`, `prism`, `pyramid` in 3D; `triangle`, `quadrilateral` in 2D); only the
    geometries with a nonzero count are listed. The sum matches `Problem`'s `MeshElements`.
  - `DomainAttributes`, `BoundaryAttributes` : Sorted lists of the unique domain (material)
    and boundary attributes present in the mesh.

The face and cell counts are reported per geometry because the number of degrees of freedom
contributed by each entity depends on its geometry at a given finite-element order. Splitting
the counts this way makes the exact global degree-of-freedom size recoverable for mixed or
non-simplicial meshes, where triangle and quadrilateral faces (and tetrahedral,
hexahedral, prism, and pyramid cells) each contribute a different per-order count.
