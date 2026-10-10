```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Boundary Conditions

## Perfect electric conductor (PEC) boundary

The perfect electric conductor (PEC) boundary condition (zero tangential electric field) is
specified using the `"PEC"` boundary keyword under
[`config["Boundaries"]`](../config/reference.md#config-boundaries-pec). It is a
homogeneous Dirichlet boundary condition for the frequency or time domain finite element
formulation, as well as the magnetostatic formulation.

For electrostatic simulations, the homogeneous Dirichlet boundary condition is prescribed
using the [`"Ground"`](../config/reference.md#config-boundaries-ground) boundary
keyword which prescribes zero voltage at the boundary.

## Perfect magnetic conductor (PMC) boundary

The perfect magnetic conductor (PMC) boundary condition (zero tangential magnetic field) is
a homogeneous Neumann boundary condition for the frequency or time domain finite element
formulation, as well as the magnetostatic formulation. It is the natural boundary condition
and thus it has the same effect as not specifying any additional boundary condition on
external boundary surfaces. It can also be explicitly specified using the `"PMC"` boundary
keyword under [`config["Boundaries"]`](../config/reference.md#config-boundaries-pmc).

Likewise, for electrostatic simulations, the homogeneous Neumann boundary condition implies
a zero-charge boundary, and thus zero gradient of the voltage in the direction normal to the
boundary. This is specified using the `"ZeroCharge"` boundary keyword under
[`config["Boundaries"]`](../config/reference.md#config-boundaries-zerocharge).

## Impedance boundary

The impedance boundary condition is a mixed (Robin) boundary condition and is available for
the frequency or time domain finite element formulations and thus for eigenmode or frequency
or time domain driven simulation types. It is specified using the
[`"Impedance"`](../config/reference.md#config-boundaries-impedance) boundary keyword.
The surface impedance relating the tangential electric and magnetic fields on the boundary
is computed from the parallel impedances due to the specified resistance, inductance, and
capacitance per square.

## Rational impedance boundary

A surface impedance boundary whose per-square impedance is an arbitrary rational function
of frequency can be specified using the
[`"RationalImpedance"`](../config/reference.md#config-boundaries-rationalimpedance)
boundary keyword. The surface impedance ``Z_s(s) = N(s)/D(s)``, with ``s = i\omega``, is
given by the real coefficients of the numerator and denominator polynomials, ordered from
highest to lowest degree (the roots of ``N`` and ``D`` are the zeros and poles of ``Z_s``).
This generalizes the parallel-RLC [impedance boundary](#Impedance-boundary) to any passive
lumped-network response, for example a fitted rational approximation of a measured or
synthesized surface impedance. Warnings are emitted for inputs that cannot correspond to a
passive (positive-real) impedance. It is available only for frequency domain driven,
eigenmode, and boundary mode simulations. In eigenmode simulations, when the resulting
contribution to the system matrix is not polynomial in frequency, it is handled by the
nonlinear eigenvalue solver.

## Absorbing (scattering) boundary

Absorbing boundary conditions at farfield boundaries, also referred to as scattering
boundary conditions, can be applied using the `"Absorbing"` boundary keyword under
[`config["Boundaries"]`](../config/reference.md#config-boundaries-absorbing). The
first-order absorbing boundary condition is a special case of the above impedance boundary
and is available for eigenmode or frequency or time domain driven simulation types. The
second-order absorbing boundary condition is only available for frequency domain
simulations.

## Perfectly matched layer

A [perfectly matched layer (PML)](https://en.wikipedia.org/wiki/Perfectly_matched_layer) is
an absorbing region surrounding the physical domain, which typically reflects much less
than an absorbing boundary condition, at the cost of additional mesh elements. In *Palace*,
the domains of the PML regions are specified with the `"Attributes"` of
[`config["Domains"]["PML"]`](../config/reference.md#config-domains-pml), and a PML is
available for 3D frequency domain driven and eigenmode simulations. The material properties
of these domains in
[`config["Domains"]["Materials"]`](../config/reference.md#config-domains-materials)
(relative permittivity and permeability, possibly anisotropic, and loss tangent) define the
background medium matched by the layer. The waves are attenuated before they reach the outer
boundary of the layer, which can be left as a PEC or natural (PMC) boundary.

In most cases, the default parameters are appropriate, and it is sufficient to mesh the layer
as a shell of a few elements around a box-shaped physical domain and to list its domains:

```json
"Domains":
{
  "Materials":
  [
    { "Attributes": [1, 2], "Permittivity": 1.0 }
  ],
  "PML": { "Attributes": [2] }
}
```

The coordinate stretch is the same in all PML regions: the PML is only reflectionless if the
stretch is the same function of position in the whole layer. The domains extending into the
layer, for example a substrate and the vacuum above it, must be PML regions in the layer
(which is checked), each with the background material properties of its material.

The PML terms replace the bulk material terms of the PML regions in the system matrices.
Otherwise, the PML regions have the material properties of their background material: for
the postprocessed fields, the error estimators, and the boundary conditions on their
boundaries, which are not transformed by the stretch (a warning is printed for boundary
conditions other than PEC on the boundaries of the PML regions). The domain energies exclude
the PML regions.

The uniaxial PML terminates an axis-aligned, box-shaped physical domain by stretching the
coordinate normal to each face of the box with the complex factor
``s = \kappa + \sigma / (\varepsilon_0 (\alpha + i\omega))``, where ``\kappa``,
``\sigma``, and ``\alpha`` are graded polynomially from their values at the interface with
the physical domain (1, 0, and 0) to `"KappaMax"`, `"SigmaMax"`, and ``2\pi`` `"AlphaMax"`
at the outer edge of the layer. Layers on several faces, including the edge and corner
regions of the box, are supported. By default, the faces and thicknesses of the layer are
detected by comparing the bounding boxes of the physical (non-PML) region and of the whole
mesh, and `"SigmaMax"` is computed from the target reflection coefficient at normal
incidence `"ReflectionTarget"`, for the smallest refractive index of the materials of the
PML regions. Alternatively, they can be specified
using `"Direction"`, `"Thickness"`, and `"SigmaMax"`. A real stretch `"KappaMax"` > 1
improves the absorption of evanescent waves, when near fields reach the layer, but reduces
the accuracy for propagating waves on a given mesh.

By default, the stretch factors are evaluated at a fixed `"ReferenceFrequency"`: the lowest
frequency for driven simulations, or the target frequency for eigenmode simulations. This
static PML keeps the system matrices independent of frequency, and so eigenmode problems
remain linear, while the absorption in the layer grows with frequency, so that the
reflection target is met at all frequencies. With `"FrequencyDependent": true`, the stretch
factors are instead evaluated at the solve frequency, for an absorption independent of
frequency, and a complex frequency shift `"AlphaMax"` can be added (CFS-PML), which reduces
the absorption below this frequency. For eigenmode simulations, the solve frequency is the
complex eigenfrequency, and the resulting nonlinear eigenvalue problem is solved with the
nonlinear eigenvalue solver as for other frequency-dependent boundary conditions. For
adaptive frequency sweeps, the frequency-dependent PML terms of the reduced-order model are
evaluated from a rational fit on the frequency band, and for circuit synthesis they are
approximated by a quadratic polynomial in frequency and, if required by the tolerance, poles
on the imaginary frequency axis, which add decaying auxiliary states to the synthesized
circuit.

The PML equations are difficult for the iterative solvers of *Palace*. In a strongly absorbing
layer, where the conductivity exceeds about ``\omega \varepsilon_0``, the anisotropic PML
tensors have components with complex phases of opposite signs, for which the polynomial
smoothers of the geometric multigrid preconditioner do not converge. By default
([`config["Solver"]["Linear"]["PMLSubdomainSolver"]`](../config/reference.md#config-solver-linear-pmlsubdomainsolver)),
the multigrid preconditioner is therefore complemented by a sparse direct solve of the
equations of the unknowns of the PML elements (except those shared with the physical
region), which restores the convergence of the multigrid preconditioner without PML. This
is affordable when the PML regions hold a small fraction of the unknowns, as is typical when
the mesh is coarse in the layer compared to the physical region. Otherwise, a sparse direct solve of the whole system (`"MGMaxLevels": 1`
and `"ComplexCoarseSolve": true`) costs about as much. Simulations with PML regions
therefore require a sparse direct solver (SuperLU_DIST, STRUMPACK, or MUMPS): the
real-valued approximation of the system matrix used by the AMS solver is indefinite in
strongly absorbing layers, and AMS is not supported with PML regions.

Sample configurations are provided in the
[`examples/pml_waveguide`](https://github.com/awslabs/palace/blob/main/examples/pml_waveguide),
[`examples/pml_layered`](https://github.com/awslabs/palace/blob/main/examples/pml_layered),
and [`examples/pml_oblique`](https://github.com/awslabs/palace/blob/main/examples/pml_oblique)
directories: a rectangular waveguide terminated by a PML (driven, adaptive sweep with circuit
synthesis, and leaky-cavity eigenmode simulations), a layered substrate and vacuum
parallel-plate waveguide whose material interface crosses the PML, and a periodic cell with a
Floquet port for plane-wave absorption at oblique incidence.

## Finite conductivity boundary

A finite conductivity boundary condition can be specified using the
[`"Conductivity"`](../config/reference.md#config-boundaries-conductivity) boundary
keyword. This boundary condition models the effect of a boundary with non-infinite
conductivity (an imperfect conductor) for conductors with thickness much larger than the
skin depth. It is available only for frequency domain driven and eigenmode simulations. For more
information see the
[Other boundary conditions](../reference.md#Other-boundary-conditions) section of the
reference.

## Periodic boundary

Periodic boundary conditions on an existing mesh can be specified using the
["Periodic"](../config/reference.md#config-boundaries-periodic) boundary keyword. This
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

  - [`config["Boundaries"]["LumpedPort"]`](../config/reference.md#config-boundaries-lumpedport) :
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
    [`"LumpedPort"`](../config/reference.md#config-boundaries-lumpedport).

  - [`config["Boundaries"]["WavePort"]`](../config/reference.md#config-boundaries-waveport) :
    Numeric wave ports are available for frequency domain driven and eigenmode simulations. In this case,
    a port boundary condition is applied with an optional excitation using a modal field
    shape which is computed by solving a 2D boundary mode eigenproblem on each wave port
    boundary. This allows for more accurate scattering parameter calculations when modeling
    waveguides or transmission lines with arbitrary cross sections.

    The 2D wave port eigenproblem supports PEC, PMC, impedance, absorbing, and conductivity
    boundary conditions. Boundaries that are specified as `"PEC"` in the full 3D model and
    intersect the wave port boundary will be considered as PEC in the 2D boundary mode
    analysis. Impedance (`"Impedance"`), absorbing (`"Absorbing"`), and conductivity
    (`"Conductivity"`) boundaries are treated as Robin boundary conditions in the wave port
    eigenvalue problem, matching the standalone boundary mode solver.
    [`config["Boundaries"]["WavePortPEC"`](../config/reference.md#config-boundaries-waveportpec)
    allows forcing specific boundary attributes to act as PEC in the wave port solve,
    overriding any other boundary condition (e.g. impedance or absorbing) that may be
    assigned to the same attributes. In addition, boundaries of wave ports other than the
    wave port currently being considered, in the case wave ports are touching and share one
    or more edges, are also considered as PEC for the wave port boundary mode analysis.

    Unlike lumped ports, wave port boundaries cannot be defined internal to the
    computational domain and instead must exist only on the outer boundary of the domain
    (they are to be "one-sided" in the sense that mesh elements only exist on one side of
    the boundary). A wave port boundary must also be planar; both the 2D boundary mode
    formulation and coordinate-based `"VoltagePath"` postprocessing use a single tangent
    frame and surface normal for the complete port face.

    The overall sign of the wave port mode E-field is internally fixed by an arbitrary
    convention that does not necessarily match the polarity convention of lumped ports
    in the same simulation. As a result, when mixing lumped and wave ports in a driven
    simulation, the cross-type S-parameters (e.g. `S_{ij}` where one port is lumped and
    the other is a wave port) may appear 180° out of phase relative to what would be
    obtained with all-lumped or all-wave ports. To pin the wave-port polarity, specify
    one of the following on each wave port — both list the **signal** (high-potential)
    terminal first and the **ground** (low-potential) terminal second:

      + [`"VoltagePath"`](../config/reference.md#config-boundaries-waveport-voltagepath): an
        ordered list of coordinate points across the port face, directed signal → ground.
        The mode is flipped so that `\int E_{\text{mode}} \cdot dl > 0` along this path.
        Also enables ``Z_{PV}`` postprocessing. Requires GSLIB.
      + [`"PolarityAttributes"`](../config/reference.md#config-boundaries-waveport-polarityattributes):
        a pair of parent-mesh boundary attributes `[signal, ground]` (e.g. distinct PEC
        attributes for the two terminals). The mode is flipped so that the modal E-field
        points from the signal attribute toward the ground attribute. Lightweight
        polarity-only alternative to `"VoltagePath"` (no GSLIB).

  - [`config["Boundaries"]["FloquetPort"]`](../config/reference.md#config-boundaries-floquetport) :
    Floquet ports are available for frequency domain driven simulations on periodic
    structures (gratings, metasurfaces, photonic crystals). They provide an absorbing
    boundary condition that decomposes the scattered field into Floquet diffraction orders
    and extracts power-normalized S-parameters for each propagating order.

    Floquet ports require periodic boundary conditions to be configured under
    [`config["Boundaries"]["Periodic"]`](../config/reference.md#config-boundaries-periodic).
    The `"FloquetWaveVector"` in the periodic configuration specifies the tangential
    component of the incident wave vector, which determines the angle of incidence. For
    normal incidence, set the wave vector to zero. For frequency sweeps at a fixed angle
    of incidence, set `"FloquetReferenceFrequency"` to the frequency (in GHz) at which the
    wave vector is defined. The wave vector then scales linearly with frequency according
    to ``\bm{k}_F(f) = \bm{k}_{F,\mathrm{ref}} f / f_\mathrm{ref}``, where
    ``\bm{k}_{F,\mathrm{ref}}`` is `"FloquetWaveVector"` and ``f_\mathrm{ref}`` is
    `"FloquetReferenceFrequency"`, maintaining a constant incidence angle across the sweep.

    The incident field is a plane wave in the specular (0,0) diffraction order with
    user-specified polarization (TE, TM, or circular RHC/LHC). S-parameters are extracted
    for all propagating diffraction orders within `"MaxOrder"` and reported in the
    `port-floquet-S.csv` output file. Each mode is labeled as
    `S[P<port>(<m>,<n>)<pol>][<exc>]` where `<port>` is the port index, `(<m>,<n>)` is
    the diffraction order, `<pol>` is the polarization (TE/TM or RHC/LHC for circular
    excitation), and `<exc>` is the excitation index. Values of `nan` are given to non-propagating
    modes.

    The port boundary must be planar and on the true boundary of the computational domain.
    The medium adjacent to the port must be homogeneous and isotropic.

For each port, the excitation is normalized to have unit incident power over the port boundary
surface. Lumped ports with reactive elements (`"L"` and/or `"C"`, or the surface-parameter
equivalents `"Ls"` and/or `"Cs"`) may also be excited, including purely reactive ports with
`"R": 0` (`"Rs": 0`). For a purely reactive port (`R = 0`) the power normalization
references an internal unit impedance; see the [mathematical reference](../reference.md#Lumped-ports-and-wave-ports)
for the S-parameter interpretation in this case.

The presence of an incident excitation at a port is controlled by the settings
[`config["Boundaries"]["LumpedPort"][]["Excitation"]`](../config/reference.md#config-boundaries-lumpedport)
and [`config["WavePort"][]["Excitation"]`](../config/reference.md#config-boundaries-waveport).
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
index (e.g. `Φ_elec[1][5] (C)` denoting the flux through surface 1 of excitation 5). The far-field
file (`farfield-rE.csv`) is the exception: it has one row per (frequency, angle) pair so the
excitation index is encoded as a row column (`exc`) rather than a column suffix. Note that a port
can only be part of one excitation.

!!! warning "Indexing"

    Any `"Index"` of [`"LumpedPort"`](../config/reference.md#config-boundaries-lumpedport),
    [`"WavePort"`](../config/reference.md#config-boundaries-waveport),
    [`"FloquetPort"`](../config/reference.md#config-boundaries-floquetport),
    [`"SurfaceCurrent"`](../config/reference.md#config-boundaries-surfacecurrent), or
    [`"Terminal"`](../config/reference.md#config-boundaries-terminal) must be unique, including between
    different boundary conditions types (e.g. you can not have a lumped port and wave port both with
    `Index: 5`).

## Surface current excitation

An alternative source excitation to lumped or wave ports for frequency and time domain
driven simulations is a surface current excitation, specified under
[`config["Boundaries"]["SurfaceCurrent"]`](../config/reference.md#config-boundaries-surfacecurrent).
This is the excitation used for magnetostatic simulation types as well. This option
prescribes a unit source surface current excitation on the given boundary in order to
excite the model. It does does not prescribe any boundary condition to the model and only
affects the source term on the right hand side.

## Flux boundary

Flux loop boundary conditions are available for magnetostatic simulations and are specified
using the [`"FluxLoop"`](../config/reference.md#config-boundaries-fluxloop) boundary
keyword. This boundary condition prescribes magnetic flux through specified holes in
conducting surfaces, enabling inductance matrix extraction for flux-based excitations. The
flux loop boundary condition works by:

 1. **Identifying holes**: Mesh boundary attributes specify holes through which flux is
    prescribed
 2. **Constraining the fluxoid**: The fluxoid through each hole (the magnetic flux plus the
    London kinetic term, reducing to the magnetic flux in the perfect-conductor limit) is set
    to the specified nondimensional excitation amplitude
 3. **Building the drive**: A curl-free cohomology generator on the film carries that fluxoid
    around the hole and drives the shifted London penalty
 4. **Computing inductance**: The resulting 3D field solutions enable inductance matrix
    extraction

Flux-loop excitations can be combined with surface-current excitations in the same
magnetostatic simulation, in which case Palace reports a single inductance matrix covering
both port types. Every surface-current element must provide an oriented `"Aperture"` with
`"Attributes"` and a Cartesian `"Direction"`, and all current ports must be `"Open"` when
inactive; see [Magnetostatic problems](problem.md#Magnetostatic-problems) for details. The
`"FluxAmounts"` entries are nondimensional internal excitation amplitudes, analogous to
unit-current excitations for current-driven magnetostatic solves. Palace reports
`terminal-Phi.csv` in physical webers.

!!! note "Flux loop requirements"

    Flux loop boundaries require:

      - Metal surface attributes defining the conducting surface containing the holes
      - Hole attributes specifying the boundaries through which flux is prescribed
      - Flux amounts defining the magnetic flux through each hole
      - Loop normal vector defining the flux orientation

The mesh must be topologically compatible with the flux loop geometry, with holes properly
defined as boundary surfaces within the conducting region. Currently, only planar holes are
supported; the hole axis is taken from the flux loop `"Direction"` and may point along any
Cartesian direction. Both conformal and nonconformal adaptive mesh refinement are supported.
