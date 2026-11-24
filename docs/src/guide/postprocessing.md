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
or time step, for examples). For a frequency-dependent dielectric, electric energy includes
polarization storage according to the [material model](model.md#Energy-and-postprocessing),
evaluated at the real part of the frequency for eigenmodes.

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
voltage and current corresponding to the incident wave on excited port boundaries. When the
weight of an excited lumped port is not real, that is, its
[`"ExcitationPhase"`](../config/reference.md#config-boundaries-lumpedport) is not a multiple
of 180 degrees, the incident voltage and current are complex and written as real and
imaginary parts, columns `Re{V_inc[...]}` and `Im{V_inc[...]}` in `port-V.csv`, and
likewise for `port-I.csv`; see
[Lumped and wave port excitation](boundaries.md#Lumped-and-wave-port-excitation).

Additionally, when surface current excitations are present, the excitations are written to
`surface-I.csv`.

For frequency domain problems, the values output are the complex-valued peak voltages and
currents, computed from the field phasors. The current of a lumped port is ``V/Z``, the current
through its own `"R"`, `"L"`, and `"C"`, and not the current delivered to the model. Quadratic
quantities (powers, energies, Poynting vectors) are time averages over one period, and each
excited port delivers 1 W of time-averaged incident power.

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

### Limits of stored-energy accounting

Boundary-mode energies contain only the transverse electric and longitudinal magnetic
components, not the full modal energy. The electric energy includes the polarization
energy of a frequency-dependent material that has a stored energy; for a material without
one, it is not reported (see [Energy and postprocessing](model.md#Energy-and-postprocessing)).
Lumped L/C energies use the circuit voltage and branch current, which can differ from the
energy stored in the boundary sheet when the port field is not uniform. The energy stored
in `"Impedance"`, `"RationalImpedance"` and `"Conductivity"` boundaries and, in the
frequency domain, the kinetic energy of London superconductors are not included. These
columns therefore do not give a complete energy balance.

## Boundary postprocessing

Boundary postprocessing capabilities are enabled by including objects under
[`config["Boundaries"]["Postprocessing"]`](../config/reference.md#config-boundaries) in the configuration
file. These include:

  - [`config["Boundaries"]["Postprocessing"]["SurfaceFlux"]`](../config/reference.md#config-boundaries-postprocessing-surfaceflux) :
    Postprocess the integrated flux through a surface defined by a list of boundary
    attributes. Electric, magnetic, and power flux are all supported. Surface capacitance
    can be computed by dividing the computed electric flux by the excitation voltage (for a
    lossy dielectric the flux uses the complex permittivity, so take the real part), while
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
    written to the file `surface-Q.csv` in the specified output directory. When
    `EdgeDistances` and either `EdgeAttributes` or `AutomaticEdges` are specified, Palace
    also writes edge-distance diagnostics to `surface-Q-edge.csv`. For each matching
    radius ``R``, these contain the interface energy and participation outside ``R`` and
    in the annulus ``R \leq d < 2R``, where ``d`` is the distance to the selected metal
    perimeter. In 3D, `AutomaticEdges` discovers PEC, terminal, conductivity, and
    impedance metal surfaces, removes perimeter segments caused by simulation-domain
    truncation, and associates the remaining physical segments with each explicit SA,
    MS, or MA interface. The manual `EdgeAttributes` path remains available;
    `EdgeExcludeAttributes` can remove its perimeter segments on artificial boundaries
    such as extrusion end faces, and `EdgeExcludeSegments` removes the perimeter segments
    coincident with explicitly listed 3D segments (a coupon's foreign edges). For
    localized automatic diagnostics, Palace infers each
    segment's process normal from its supporting metal faces and orients it from the
    lower-wave-speed substrate side toward the higher-wave-speed air side.
    `EdgeFrameNormal` remains an optional fallback for materially ambiguous surfaces.
    `EdgeDistanceSmoothing` can
    replace the hard distance indicators with complementary cubic transitions to reduce
    mesh-to-mesh cutoff noise while preserving the inside/outside partition. With
    `LocalizeEdgeEnergy`, Palace additionally writes `surface-Q-edge-local.csv`, which
    assigns the total interface energy and resolves the inside, annular, and surrounding
    volume-tube energies by nearest physical perimeter segment for use in geometry-aware
    edge corrections. For automatically extracted edges, it also integrates the portions
    of the surface core and volume annulus within an along-edge distance ``R`` of a
    corner, endpoint, or junction. These exact vertex-neighborhood energies support
    confidence estimates for geometry not represented by a straight-edge calibration.
    `EdgeRefinement` can enforce a geometry-based element-size target
    throughout the correction tube before the first solve and downweight field-error
    indicators wholly inside the replaced singular core during later AMR. For MA, MS,
    and SA interfaces,
    `FluxRecovery` evaluates the normal displacement field using an H(div)-conforming
    recovered flux; the SA tangential term continues to use the native electric field.
    This option is rejected on interior boundary sheets, whether or not the mesh is
    cracked along them, and on cracked internal boundaries, because a volume L2
    projection does not provide a controlled normal trace on the singular PEC surface.
    Electrostatic fabrication-coupon simulations can combine `LocalizeEdgeEnergy` with
    prescribed-potential basis excitations and
    [`ResponseMatrix`](../config/reference.md#config-solver-electrostatic-responsematrix).
    Palace then writes `surface-response-matrix.csv`, the upper triangles of the localized
    quadratic core- and full-patch surface-energy operators for every interface, edge, and
    matching radius. Coupon generation enables `AggregateResponseMatrix` to sum physical
    edges by interface and matching radius before writing. It also writes
    `domain-response-matrix.csv`, the domain-energy response over the same prescribed
    potential basis. Driven, eigenmode, and
    two-dimensional boundary-mode simulations can apply a compatible
    fabrication-process library with
    [`Solver.SurfaceResponseCorrection`](../config/reference.md#config-solver-surfaceresponsecorrection).
    This postprocessing-only correction writes `surface-Q-corrected.csv` and
    `surface-response-confidence.csv`; the solved field, AMR sequence, and standard
    `surface-Q.csv` remain the uncorrected thin-metal result. Post-processing is the
    default and the designer's path; the self-consistent corrected solve and the other
    options set aside from it are developer options, refused by name unless the block
    sets `"Developer": true` (see the developer switch in
    [the thin-metal response correction](../reference.md#Self-consistent-thin-metal-response-correction)).
    PEC, conductivity, impedance, and rational-impedance boundaries require coupons with matching
    boundary-law metadata; the default path serves PEC and impedance (Ls) sheets, while
    conductivity and rational-impedance coupons and the boundary-mode correction are developer
    options. See
    [Maxwell response correction](../reference.md#Maxwell-response-correction)
    for the supported geometry and confidence limits.
    The energies (J) of these outputs are time averages, as the other energies of a
    frequency-domain solution: in eigenmode and boundary-mode simulations they are one
    half of the values of *Palace* before version 0.20.0, in port-driven simulations they
    are unchanged (the ports deliver 1 W of time-averaged incident power); the
    participations, shares and statuses do not change. A frequency-dependent material
    `"Permittivity"` is refused with the correction, and with the interface keys it builds
    on (`"EdgeDistances"`, `"LocalizeEdgeEnergy"`, `"FluxRecovery"`, the ownership keys and
    `"Region"`), in this version.
    Every corrected run also reports an edge-resolution census of the trace contour
    (`surface-response-resolution.csv`, `SurfaceResponse.Resolution` in `palace.json`, a
    PREFLIGHT line and a banner; see
    [Resolution census of the response correction](../reference.md#resolution-census-of-the-response-correction)):
    RESOLVED is necessary, not sufficient; read it only from a block with `"Final": true`
    and the intended `Level`.
    Every corrected run (both solve families) also reports its coverage: ``W`` is the
    raw surface energy within the matching radius of an interface and ``U`` the part of
    it where no correction is applied. `surface-response-coverage.csv` gives per cell
    (mode, frequency and excitation, or source) and interface `Type` the uncovered share
    ``U/W`` split into the missing share (features whose library model is `Missing`) and
    the method gap, the pending-model and no-status-model shares, the covered correction
    factor, an interim coverage status and a regime status;
    `surface-response-uncovered-energy-by-key.csv` lists the uncovered energy per missing
    key, `surface-response-near-edge-energy-budget.csv` gives per cell where every part of
    ``W`` went (each bin of the coverage partition once: replaced by a model, kept raw, or
    dropped by the corrected total: uncovered is not the same as kept raw), and
    `terminal-C-corrected-coverage.csv` the status of each corrected capacitance
    entry. The same columns are appended to `surface-Q-corrected.csv` and
    `surface-response-confidence.csv`, and every cell that is not `OK`, or whose regime
    status is not `IN` (computed, NOT VALIDATED in this regime), is named by the end-of-run
    banner (at most 40 cell lines, then a line counting the rest) and in `palace.json`
    (`SurfaceResponse.Coverage.Flagged`: every such cell). A coverage status is about coverage only; the
    regime status says whether the covered part's correction factor lies in the range
    validated for this library and excitation class; NEITHER says that the mesh of the
    solve resolves the field at the correction radius. See
    [Coverage of the response correction](../reference.md#Coverage-of-the-response-correction).

## Visualization

When specified in the configuration file, the electric field and magnetic flux density
solutions are written to disk for 3D visualization with [ParaView](https://www.paraview.org/)
or [GLVis](https://glvis.org/). Various other postprocessed fields are also written to the ParaView
or grid function (GLVis) database as available, including electric and magnetic energy density,
surface currents, and charge density. These files are found in the `paraview/` or `gridfunction/`
directories located in the output directory specified under
[`config["Problem"]["Output"]`](../config/reference.md#config-problem). The output
formats are specified in [`config["Problem"]["OutputFormats"]`](../config/reference.md#config-problem).
Driven simulations with several excitations write the fields of each excitation to a separate
subfolder, see [Lumped and wave port excitation](boundaries.md#Lumped-and-wave-port-excitation).

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
U_e = \frac{1}{4} \bm{E}^{*} \bm{W} \bm{E}, \qquad
U_m = \frac{1}{4} \bm{B}^{*} \bm{\mu}^{-1} \bm{B}, \qquad
\bm{S} = \frac{1}{2} \operatorname{Re}\!\left\{\bm{E} \times
\left(\bm{\mu}^{-1}\bm{B}\right)^{*}\right\}
```

for the complex peak phasors of a frequency domain solution (time-averaged quantities), or
``U_e = \frac{1}{2} \bm{E}\cdot\bm{\epsilon}\bm{E}``, ``U_m = \frac{1}{2} \bm{B}\cdot\bm{\mu}^{-1}\bm{B}``,
and ``\bm{S} = \bm{E}\times\bm{\mu}^{-1}\bm{B}`` for the real-valued instantaneous fields of a
time domain or static solution.

Here ``\bm{W}=\operatorname{Re}\bm{\epsilon}`` for constant materials, while a
frequency-dependent model supplies its stored-energy weight, which includes polarization.
When *Palace* cannot report the stored electric energy of a frequency-dependent material
on the mesh, ``U_e`` is not written (see [the model guide](model.md)).
Electric flux and charge still use the complete complex ``\bm{\epsilon}(i\omega)``,
not ``\bm{W}``.

A jump in ``\bm{W}`` or ``\bm{\mu}^{-1}`` can therefore produce a jump in an energy
density or the Poynting vector, even where the conforming trace of a primary field is
continuous. Similarly, at a material interface:

  - H(curl) continuity constrains only the tangential electric field; its normal component
    can jump.
  - H(div) continuity constrains only the normal magnetic flux density; its tangential
    component can jump.
  - ``\bm{D} = \bm{\epsilon}\bm{E}`` and
    ``\bm{H} = \bm{\mu}^{-1}\bm{B}`` inherit jumps from both the fields and the material
    tensors. For a lossy dielectric in the frequency domain, ``\bm{\epsilon}`` is the complex
    permittivity ``\bm{\epsilon}_r(1 - i\tan\delta)``, so the electric flux and surface
    charge of the complex field include the lagging polarization.
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

When Maxwell surface-response correction supplies a self-consistent corrected field (a
developer option: `CorrectionMode` `SelfConsistent` or `Both` with `"Developer": true`),
Palace preserves the ordinary ``E`` output as the raw thin-metal solution and additionally
writes ``E_corrected_real`` and ``E_corrected_imag``. If recovered electric flux is needed
by interface postprocessing, ``D_corrected_real`` and ``D_corrected_imag`` are also written.
These names are used by both ParaView and MFEM GridFunction output.

The boundary data collection includes the 3D field values sampled on the boundary mesh as
well as:

  - Surface charge density: `Q_s`, `Q_s_real`, `Q_s_imag` (C/m²)
  - Surface current density: `J_s`, `J_s_real`, `J_s_imag` (A/m)
  - Wave port boundary mode electric field: `E0_<i>_real`, `E0_<i>_imag` for wave port `i`
    (V/m), the port mode normalized to 1 W of time-averaged power, which is the incident
    field of an excitation of the port

## Adaptive mesh refinement

At the start of an adaptive mesh refinement (AMR) iteration, if
[`config["Model"]["Refinement"]["SaveAdaptIterations"]`](../config/reference.md#config-model-refinement)
is enabled, the postprocessing results from the solve on the previous mesh will be saved off
within a subdirectory denoted `iterationX`, where `X` is the (1-based) iteration number.
The results in the top level directory will always be those from the most recent successful
solve. When AMR stops because the predicted size of a refined mesh exceeds
[`config["Model"]["Refinement"]["MaxSize"]`](../config/reference.md#config-model-refinement),
no `iterationX` subdirectory is created for the skipped solve: the top level directory
holds the results of the last solve and `palace.json` records the stop under
`"AdaptiveMeshRefinement"`/`"MaxSizeStop"` (the AMR iteration, the last solved and the
predicted numbers of degrees of freedom, the element counts before and after the refinement
step, and `"MaxSize"`).

When [`config["Model"]["Refinement"]["SaveAdaptMesh"]`](../config/reference.md#config-model-refinement)
is enabled and adaptation was performed, *Palace* saves the adapted mesh to disk and also
records a `SavedAdaptedMesh` block in `palace.json` summarizing its topology. When *Palace*
is built with `zlib`, the saved mesh is gzip-compressed and written with a `.meshgz`
extension; without `zlib` it is written uncompressed as `.mesh`. MFEM detects the compression
from the file contents when reading, so either file can be used directly as a restart mesh. A
nonconforming mesh is written by all processes in parallel: the file keeps the refinement
tree of the input mesh, with its roots listed process by process, so its numbering depends on
the partition of the mesh, but not the mesh it describes. When such a mesh is used as the mesh
of a restart, *Palace* reads it in parallel, which takes much less time and memory for large
meshes, unless the restart modifies the mesh (for example with `"MakeSimplex"`, refinement
regions, or periodic boundaries). The root process decompresses a compressed file and sends
each process its part; for very large meshes, decompressing the file beforehand (for example
with `gunzip`) is faster, since every process then reads its part of the file directly.
Meshes which cannot be written in parallel (for example NURBS meshes, or three-dimensional
meshes with anisotropic refinement) are gathered on the root process instead, with a warning.
The block always describes the mesh saved next to it: the top-level `palace.json` describes the final adapted
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

### Following the corrected values over adaptive rounds

A run with a surface-response correction (electrostatic, eigenmode or driven) writes one
table, `surface-participation.csv`, at the top level of the output directory and appends to
it at the end of every adaptive level, so you can start from a coarse mesh and follow the
values over the rounds as you would follow raw values. Each row is one interface of one
point (terminal, mode or output frequency) at one level: the raw and the corrected
participation side by side, then in the same row the resolution verdict (with the element
size near the metal edges over ``R``), the coverage and regime verdicts, the confidence
flag, the change since the previous level and the convergence indicator. A box target (an
entry with a `Region`) gets verdicts of its own patches, so a box on coarse elements is
flagged even when the rest of the device is resolved. A box's resolution words read
`-UNTESTED` for now (`MARGINAL-UNTESTED`, for example): the bars were measured with the
census of the box's whole cell, which reads coarser than a junction box's own patches, and
are re-measured on the box's own census before its words lose the suffix. Nothing is
withheld: a value whose verdicts are not met is printed with them listed in `flags`.

Read the rows whose `scope` is `Type` (the corrected participation of a `Type`, summed over
its interfaces) or `region` (a box value): their verdicts are of that value; an `interface`
row shows one interface's share and carries its `Type`'s verdicts. The values are the
fixed-trace corrected participations (`p corrected (fixed-trace)`).

The convergence indicator is the change per doubling of the unknowns, read against the
latest earlier level with at least 20 % fewer unknowns and only when the mesh at that
value's own edges was refined since then (otherwise `MOVING: CELL-NOT-REFINED`: on a large
device the unknowns can grow elsewhere). A value is SETTLED below 1 % per doubling, and
RESOLVED only when it is SETTLED and the mesh near its edges is within the bar of its order.
A SINGLE SOLVE NEVER READS RESOLVED, so run adaptive rounds (order 3 with the default
adaptive rounds is the cheapest route measured). Stop on RESOLVED rows of the values you
need, not on a small plain change. What the words mean, measured on electrostatic values of
three designs against thick-metal references: RESOLVED means the value stopped moving inside
the tested domain, not that it matches thick: the worst value printed RESOLVED was 4.0 % off
thick and 26 of the 39 RESOLVED values were beyond 2 %, so a 2 % target is not supported.
SETTLED cannot see a bias of the library or of the reference: one design settled at -2.8 %
at all three orders. MARGINAL values were within 5 % in 61 % of the cases and outside in 39
%; RESOLVED came 1-2 rounds after a value was within 5 %, at 3-13 times its unknowns. The
bars are the edge of the tested domain, not sizes that give 5 %. No eigenmode value has yet
settled or read RESOLVED in the tests: eigenmode and driven rows read `SETTLED-UNTESTED` /
`RESOLVED-UNTESTED` until the rule is tested there.

For an eigenmode the history must follow the same physical mode across rounds, not the same
index. The table joins a mode to the previous level's mode nearest in frequency only when
that match is clear: the two are each other's nearest, the frequency moved by less than half
the spacing to the next mode, their participations (with the junction ports' energy
participations, the sharpest mark of a qubit mode) agree, no second mode of the same kind is
nearly as close, and the modes of that kind keep their order. A mode whose match is not
unique reads `UNMATCHED` for that level and its indicator `MOVING: UNMATCHED`: the safe
side, it never borrows a neighbour's history. A mode of the same kind entering or leaving
the solved window (for example a mode below it moving out of the `Target` range) unmatches
its neighbours for one level. The log lists every join and every `UNMATCHED` mode with both
frequencies. Data marked interim make the words read `-INTERIM`. The thresholds are data
(see [the reference](../reference.md#The-per-iteration-table-of-the-response-correction)).
The log prints the values in short form after each level. `PreflightOnly` does not evaluate
the resolution census; every corrected run prints it before its first solve and in the first
rows of this table. A run restarted from a saved adapted mesh starts a new history (its
first level reads `MOVING: FIRST-ITERATION`) and replaces the table in its output directory:
give it another `Output` directory to keep the earlier table.
