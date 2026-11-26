```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Substructuring for Region Redesign

When a design iterates on one part of a large device (one qubit on a chip, say), every
iteration normally solves the whole device again, although most of it never changes.
*Palace*'s substructuring splits the domain into a **region** (the part being redesigned) and
an **environment** (everything else). The environment is condensed onto its interface with the
region, once; each redesign then only solves the region against the condensed environment.

Substructuring is available for [electrostatic](../guide/problem.md#Electrostatic-problems) capacitance
extraction, [magnetostatic](../guide/problem.md#Magnetostatic-problems) inductance extraction with
flux-loop and surface-current excitations, and
[driven](../guide/problem.md#Driven-problems-in-the-frequency-domain) frequency sweeps (see
[Driven problems](#Driven-problems)).

## How it works

Write the discrete problem on the degrees of freedom of the region interior ``R``, the shared
interface ``\Gamma``, and the environment interior ``E``. Eliminating ``E`` leaves the region
problem with one extra term on the interface, the environment's Schur complement

```math
\bm{S}_E = \bm{A}_{\Gamma\Gamma} - \bm{A}_{\Gamma E}\,\bm{A}_{EE}^{-1}\,\bm{A}_{E\Gamma},
```

a dense ``|\Gamma| \times |\Gamma|`` matrix (the discrete Dirichlet-to-Neumann map of the
environment). The condensation is an exact algebraic rearrangement: the region solution, and the
capacitance or inductance matrix, are those of the full solve (for magnetostatics up to a small
regularization, see [Magnetostatics](#Magnetostatics)). Electrostatic terminals and `"PEC"`
boundaries (grounded, as in a regular simulation) may lie in the region, in the environment, or
span both; the environment's coupling to its terminals is condensed together
with ``\bm{S}_E``, so terminal energies need no environment solve.

The environment model (``\bm{S}_E`` and the terminal couplings) can be saved to a file. A later
run loads it and solves only the region, possibly after re-meshing the region: its solves then
follow the size of the region rather than of the device. The run checks that its environment is
the model's: its elements and boundaries, the configuration on its attributes (materials, boundary
conditions, order and quadrature) exactly, and on the full mesh its operator, applied matrix-free
without assembling it; on the mesh of the region written with the model, only the region is read.
The formulation is written out in the [reference](../reference.md#Substructuring).

## Configuration

Substructuring is enabled with a
[`config["Solver"]["Substructuring"]`](../config/reference.md#config-solver-substructuring) block. The
region and the environment are given by domain attributes, which together must cover the mesh
(`"Environment"` defaults to the domain attributes not in `"Region"`):

```json
"Substructuring":
{
  "Region": { "Attributes": [1, 2] },
  "Environment": { "Attributes": [10, 11] }
}
```

  - A run without `"ModelFile"` is offline: it condenses the environment and saves its model as
    `environment.model` in
    [`config["Problem"]["Output"]`](../config/reference.md#config-problem-output) (none with
    `"SaveModel": false`).
  - A run with `"ModelFile"`, the path of a saved model, is online: it loads the model and solves
    only the region. The region may be re-meshed between runs, as long as the environment and the
    interface are unchanged.

An offline run that saves a model on a conforming mesh (not one from adaptive refinement) also
writes a mesh of the region next to it, `environment_region.mesh`, with the few environment elements
along the boundary conditions of the environment that reach the interface and, in magnetostatics, at
the surface currents of the environment and at the flux loops whose film or hole reaches the region
or the interface (with these films whole and, if these elements reach a `"TwoSided"` film where it
crosses the interface, the environment along every such crossing) and next to the films on the
interface. The flux loops entirely in the environment are left out: their excitations and the
fluxoids they impose come from the model. An online run on it (as `"Model"/"Mesh"`, with the
configuration otherwise unchanged) reads and partitions the region only, and tells it from the full
mesh by the environment elements it holds. The environment is then that of the model: a
configuration that changes it is rejected (an excitation of the environment that changed needs the
full mesh, and in driven problems a new offline run), and the mesh must hold the environment
elements written with it; the interface is found from the saved interface DOFs, and the
nondimensionalization is that of the model. Every excitation with sources in the environment must be
in the model, the fields are written in the region and on the interface only, and a priori and
adaptive refinement need the full mesh (the mesh of the region is written refined).
A redesign of the region is re-meshed inside the mesh of the region with
[`examples/substructuring/remesh_region.py`](https://github.com/awslabs/palace/blob/main/examples/substructuring/remesh_region.py)
and `--region-mesh` (see [Example: qubit lattice](@ref substructuring-lattice-example)).

In electrostatics and magnetostatics (driven runs are described in
[Driven problems](#Driven-problems)), the remaining settings (terminals, materials,
`"Solver"/"Order"`, ...) are those of the regular simulation, except `"KSPType"`, `"Tol"` and
`"MaxIts"` of `config["Solver"]["Linear"]`, which are not used; its `"Type"` only selects the sparse
direct solver (see [Solvers](#Solvers)), and another type is not used (*Palace* warns). The
capacitance and inductance matrices are written to the files of a regular simulation
(`terminal-C.csv`, `terminal-M.csv`, and the others). The domain and surface postprocessing of a
regular simulation is that of the region, from its field and the total energy of the condensation,
so that the participation ratios are those of the whole model: the energies of the domains of the
region and the probes in it or on the interface (`config["Domains"]["Postprocessing"]`, a probe on
the interface with the field of the region), and the `"SurfaceFlux"` and `"Dielectric"` entries of
`config["Boundaries"]["Postprocessing"]` on surfaces in it. The others are left out with a warning,
and `"FarField"` is rejected. In magnetostatics, the fields keep the regularization of the solves
(see the [reference](../reference.md#Substructuring)): without superconducting sheets or domains,
the quantities from them, such as the energies of the domains, the probes and the fluxes, are within
about ``10^{-4}`` of a regular simulation, while the inductance matrix, from energies stationary in
the regularization, agrees to its square. Field output follows `"Solver"/"Electrostatic"/"Save"` and
`"Solver"/"Magnetostatic"/"Save"`, in the formats of `config["Problem"]["OutputFormats"]`; the
environment part of a saved field is recovered only when fields are written.

### Choosing the region

Most meshes do not come with a region/environment split. Any split into two sets of domain
attributes works: the interface follows the element faces between them and need not be planar. The
cost of the condensation grows with the number of interface degrees of freedom, so a region boundary
in a coarsely meshed part of the domain (substrate or vacuum away from the metal) is cheaper than
one that cuts through a finely meshed area. Its memory does too: ``\bm{S}_E`` is dense
(``|\Gamma|^2`` values, 7 GB at ``|\Gamma| = 30{,}000``), distributed over a grid of processes
(STRUMPACK and SuperLU_DIST form it whole on one process first), and the region solve holds as many
values again. A saved model holds its lower triangle, and a driven model one per frequency or
sample, complex (see [Driven problems](#Driven-problems)). One process holds the whole mesh to write
the mesh of the region. The [transmon example](@ref substructuring-transmon-example) retags the
elements of an existing mesh by the position of their centroids.

### Adaptive mesh refinement

In electrostatics and magnetostatics (driven substructuring refines no mesh), the environment is
condensed once, so its mesh resolution is fixed in the saved model. Meshes from layout tools are
often coarse and meant to be refined by
[adaptive mesh refinement](../reference.md#Error-estimation-and-adaptive-mesh-refinement-(AMR))
(AMR), so the refinement is done in two steps:

 1. Refine the full model before condensing it: run the simulation without
    `"Substructuring"`, with AMR and
    [`config["Model"]["Refinement"]["SaveAdaptMesh"]`](../config/reference.md#config-model-refinement)
    set to `true`, then condense the environment offline with the saved mesh as
    `config["Model"]["Mesh"]`. The refined elements keep the region and environment
    attributes.
 2. Refine the region only in online runs, for example after a redesign: with
    `config["Model"]["Refinement"]["MaxIts"]` greater than zero, the error indicators are
    computed in the region and the environment is left as condensed. This requires
    nonconforming refinement without a level constraint (`"Nonconformal": true`,
    `"MaxNCLevels": 0`), so that the refinement does not spread into the environment, and the
    full mesh.

Adaptive refinement in an offline run is rejected (use step 1). The indicators are those of
the regular simulation's estimator (the gradient flux in electrostatics, the curl flux in
magnetostatics), with the flux recovery restricted to the region, so that the unknown
environment field does not affect them near the interface, and normalized as there by the energy
of the field in the volume (without the kinetic energy of the London sheets and domains).

### Solvers

Substructuring also runs with a GPU device (`config["Solver"]["Device"]` `"GPU"`), with the dense
interface operator and the sparse direct factorizations on the CPU. In electrostatics and
magnetostatics, the environment and the region are each factored for their Schur complement onto the
interface by a sparse direct solver (the one
[`config["Solver"]["Linear"]["Type"]`](../config/reference.md#config-solver-linear-type) names, or
else MUMPS, STRUMPACK or SuperLU_DIST, preferred in this order among those Palace is built with), so
*Palace* must be built with one of them. MUMPS and STRUMPACK (with the Schur complement interface
that the superbuild and the Spack package patch in) give it by a partial factorization, the others
by one solve per interface unknown; MUMPS compresses its fronts at round-off (block low-rank), which
is faster than full rank at no loss of accuracy. By default, on 32 processes or more, the Schur
complement of a region with at most 10,000 interface unknowns comes from STRUMPACK first, with its
Schur complement interface: STRUMPACK gathers it on one process, and MUMPS distributes it but with
a workspace of a few tenths of a GB per process, so that from about 30 processes MUMPS takes more
memory (on 192 processes, twice as much), though it stays faster by up to a third. On fewer
processes, MUMPS is both faster and smaller. The region solve adds ``\bm{S}_E`` to the Schur
complement of the region and factors this dense interface system by Cholesky, on a grid of processes
(by ScaLAPACK): its cost grows as the cube of the number of interface unknowns, and each solve
condenses onto the interface, solves there and expands. A saved magnetostatic model also holds the
energy operator of the environment, which takes an environment solve and a condensation (the forward
elimination of the partial factorization, half a solve) per interface unknown, and with London
sheets in the environment their kinetic energy operator, for another condensation per interface
unknown and as many values again in the model. The factorizations follow the
[`"ColumnOrdering"`](../config/reference.md#config-solver-linear-columnordering) of
`config["Solver"]["Linear"]` (the partial factorizations are analyzed with the sequential
counterpart of `"ParMETIS"` or `"PTScotch"`), and take no other option of it: they are exact, so the
compression of STRUMPACK and the 3D algorithm of SuperLU_DIST are not used (*Palace* warns). When
running MUMPS with MPI, setting `OPENBLAS_NUM_THREADS=1` (or the equivalent for the BLAS in use)
avoids oversubscribing cores. Driven runs follow the same `"ColumnOrdering"`, with the sparse direct
solvers of [Driven problems](#Driven-problems).

## Magnetostatics

Magnetostatic substructuring extracts the inductance matrix from
[`config["Boundaries"]["FluxLoop"]`](../config/reference.md#config-boundaries-fluxloop) and
[`config["Boundaries"]["SurfaceCurrent"]`](../config/reference.md#config-boundaries-surfacecurrent)
excitations, alone or together (with flux loops grouped in an excitation or held at zero fluxoid, as
in a regular simulation), including the kinetic inductance of superconducting films. The `"PEC"`
boundaries are the Dirichlet boundaries of the condensation. The flux-loop films and the
[`config["Boundaries"]["Superconductor"]`](../config/reference.md#config-boundaries-superconductor)
films are London sheets, as in a regular simulation: a film without a `"Superconductor"` entry has
the penetration depth `"PecPenetrationDepth"`, and a `"Superconductor"` film its own. The sheet term
of each film face goes to one side, that of its neighboring element with the lower domain attribute
(as for every boundary term of *Palace*'s operators), so films may lie in the region, in the
environment, or cross the interface, and films that are not `"TwoSided"` may also lie on it; the
surface-current ports and their apertures may also cross it. A
[`"TwoSided"`](../config/reference.md#config-boundaries-superconductor-twosided) film, across which
the solution space is broken as in a regular simulation, may also cross the interface (the interface
DOFs on it then have a version on each side of the film, and a saved model needs a conforming mesh),
but not lie on it, and the region of an adaptive refinement cannot hold one. A saved model includes
the condensed excitations of the environment, so an online run needs no environment solve as long as
the excitations in the environment are unchanged; excitations inside the region, such as a port
added or moved in a redesign, need none at all. A small mass regularization keeps the curl-curl
operator definite. It is weighted by the inverse permeability and scaled with the extent of the
mesh, so that it scales as the curl-curl operator and the results for the part of a current that
closes depend on neither the permeability nor
[`config["Model"]["Lc"]`](../config/reference.md#config-model); the extracted
energies use the unregularized operator, in a form whose error is quadratic in the regularization.
Each surface current must close: it has to start and end on `"PEC"`
boundaries or superconducting films, as in a physical circuit, with no net current into any
connected conductor. A current that does not close has no magnetostatic solution and is rejected.

!!! note

    With substructuring, inactive surface-current ports must be `"Open"` (see
    [`config["Solver"]["Magnetostatic"]["InactivePorts"]`](../config/reference.md#config-solver-magnetostatic-inactiveports)):
    every excitation is solved against the same condensed environment.

## Driven problems

For a [driven](../guide/problem.md#Driven-problems-in-the-frequency-domain) simulation the
environment operator depends on the frequency, and substructuring condenses it exactly at each
frequency of a sweep (`"Solver"/"Driven"/"Samples"`): one partial factorization of the environment
gives ``\bm{S}_E(\omega)``, the condensation of the environment's sources onto the interface, and
the condensed voltages of the environment's lumped ports and wave port overlaps (see the
[reference](../reference.md#Driven-substructuring)). The region is then solved against them,
factored once per frequency for all excitations. Lumped ports, impedance, absorbing and conductivity
boundaries, and lossy materials may lie in the region, in the environment, or on both sides, and
wave ports inside either side, away from the interface (the modal terms of the environment's border
its factorization).

An offline sweep saves these per frequency (unless `"SaveModel": false`). An online sweep at saved
frequencies (any subset of them, matched by value) then solves only the region, which may be
redesigned (materials, re-meshing) and run on a different number of processes, as long as the
environment and the interface are unchanged and the excitations with sources in the environment
are among the model's. Offline and online, the field is known in the region and on the interface
only (a regular simulation gives the rest):

  - S-parameters of all ports, voltages and currents of the lumped ports, and the energies of
    lumped elements, are those of a regular simulation;
  - domain energies are written for the `"Domains"/"Postprocessing"/"Energy"` entries inside
    the region, without the total energies and participation ratios (and the power of the
    environment's ports, and the voltage of its wave ports, are not available);
  - surface fluxes (`"Boundaries"/"Postprocessing"/"SurfaceFlux"`) on boundaries that border only
    region elements (on the mesh of the region, those away from the environment's), and probes in
    the region or on the interface (there, with the field of the region), are those of a regular
    simulation; the others are skipped;
  - interface dielectrics (`"Boundaries"/"Postprocessing"/"Dielectric"`) on such boundaries give
    their energies in `surface-Q.csv` (`E_surf`), without the participation ratios and quality
    factors, which need the total energy;
  - fields (`"Save"`) are written in the region and on the interface, zero elsewhere;
  - `"FarField"` is rejected.

The mesh of the region written by the offline sweep (see [Configuration](#Configuration)) also
holds the environment elements at the environment's ports, whose excitations and outputs need
them. Impedance and conductivity boundaries on the interface itself need the full mesh online,
since the mesh of the region would hold them in the region alone. On large designs, reading and
partitioning the full mesh take most of an online sweep's time.

A restarted sweep ([`"Restart": x`](../config/reference.md#config-solver-driven-restart)) continues
in the order of its solves, all excitations at a frequency before the next one, after its first
``x - 1`` results. With more than one excitation, a regular sweep counts its solves the other way
(all frequencies of an excitation before the next one), so the same ``x`` resumes elsewhere. An
offline sweep continues the model it was saving, so it runs in the same configuration on as many
processes.

An [adaptive sweep](../guide/problem.md#Driven-problems-in-the-frequency-domain) (nonzero
`"AdaptiveTol"`) instead condenses the environment at frequencies chosen greedily, as the
adaptive fast frequency sweep chooses its samples, for a rational model of the condensed data: a
barycentric interpolant with one denominator for all of it. Sampling stops when the model
predicts `"AdaptiveConvergenceMemory"` consecutive new samples to `"AdaptiveTol"` (relative to
the condensed data), or at `"AdaptiveMaxSamples"` samples. The model's passivity is then checked
between the samples (the imaginary part of the condensed environment operator, the power the
environment absorbs, is positive semidefinite), and a violation beyond `"AdaptiveTol"` adds a
sample there. The adaptive fast frequency sweep then reduces the region against the model, as it
reduces a whole device, by Galerkin projection onto its fields at greedily chosen frequencies,
which gives every frequency of the sweep with the outputs of an online sweep. Saved, such a model
serves online sweeps at any frequency between its lowest and highest samples, reduced the same
way with a nonzero `"AdaptiveTol"`, or solved frequency by frequency. An online sweep warns when the
sampling of the model did not converge (it stopped at `"AdaptiveMaxSamples"` samples, or the
offline run stopped before it ended), and when the model is not passive between its samples.

A saved model, and an adaptive sweep, need the lumped ports and the other sources (surface
currents, current dipoles) away from the interface: the excitation of a source on it, and the
voltage of a lumped port, would depend on both sides.

!!! note

    Offline, driven substructuring condenses the environment with MUMPS or, without it, STRUMPACK,
    so *Palace* must be built with one of them: their partial factorizations give its Schur
    complement onto the interface directly (STRUMPACK with the Schur complement interface that the
    superbuild and the Spack package patch in, else with one solve per interface unknown).
    SuperLU_DIST, whose condensation goes wrong through a resonance, is not used offline. FGMRES
    solves the region, preconditioned by the region's own Schur complement onto the interface with
    the environment's added, a dense matrix factored across the processes by ScaLAPACK (which MUMPS
    and STRUMPACK build with; without it, on one process). On 32 processes or more, the region's
    Schur complement is STRUMPACK's, gathered on one process, up to 10,000 interface unknowns (with
    its Schur complement interface), and otherwise MUMPS's, distributed, or, without it, STRUMPACK's
    (or, without them, from SuperLU_DIST's solves). Near a resonance of the region with the
    interface held fixed, where this preconditioner is nearly singular, FGMRES takes more
    iterations, a few with MUMPS's Schur complement of the region, about ten with STRUMPACK's and
    tens with SuperLU_DIST's, within a
    limit of its own (500, or `"MaxIts"` if larger; it is used whatever
    [`config["Solver"]["Linear"]["KSPType"]`](../config/reference.md#config-solver-linear-ksptype)
    names, as GMRES loses digits there). A sparse direct solver that
    [`config["Solver"]["Linear"]["Type"]`](../config/reference.md#config-solver-linear-type) names
    is used for all of them instead. At a resonance of the environment with the interface held
    fixed, where the environment's Schur complement has a pole, its exact condensation loses digits
    as the frequency nears it, and an adaptive sweep places samples there. Floquet ports (which make
    the system non-symmetric) and periodic boundaries (on whose seam the interface degrees of
    freedom, matched by their coordinates, depend on the partition), frequency-dependent
    permittivity, PML regions, `"AdaptiveCircuitSynthesis"` and the adaptive refinement of the
    region in online sweeps are not supported yet. The model holds one dense
    ``|\Gamma| \times |\Gamma|`` matrix per frequency or sample (16 bytes per entry of its lower
    triangle). The processes hold it distributed, each a share of its lower triangle, which they
    read in parallel online and which one process writes in chunks offline; the samples of an
    adaptive sweep are distributed alike, and its passivity checks factor ``\mathrm{Im}\,\bm{S}_E``
    on a grid of the processes.

## [Example: transmon capacitance](@id substructuring-transmon-example)

The [transmon example](../examples/transmon.md) mesh contains a transmon island, a feedline, and
a ground plane with a readout resonator, all as a single metal boundary. The script
[`examples/transmon/transmon_substructuring.py`](https://github.com/awslabs/palace/blob/main/examples/transmon/transmon_substructuring.py)
splits the metal into its three conductors, the ground plane (a `"PEC"` boundary), the feedline and
the island (the terminals), and assigns the elements in a box around the qubit to the region. It
also writes the configurations of the later runs from that of the offline run,
`transmon_substructuring_offline.json`: they reuse its model (`"ModelFile"`), each with an output
folder of its own, and the redesign below takes the mesh of the redesigned region:

```bash
cd examples/transmon
python3 transmon_substructuring.py            # the mesh, and transmon_substructuring_<run>.json
palace transmon_substructuring_offline.json   # condense the environment, save the model
palace transmon_substructuring_online.json    # reuse the model: only the region is solved
```

Both runs write the ``2 \times 2`` Maxwell capacitance matrix to `terminal-C.csv`; the online
run solves only the region. The box can be changed with `--box`.

### Redesigning the region

A new design of the region needs a mesh with the same environment and interface as the
offline run. The script
[`examples/substructuring/remesh_region.py`](https://github.com/awslabs/palace/blob/main/examples/substructuring/remesh_region.py)
builds one from the split mesh of the offline run and a full mesh of the new design, as produced
by the layout tool: it keeps the environment and the interface, takes the layout inside the
region from the new design, and re-meshes the region with [Gmsh](https://gmsh.info), following the
element sizes of the new design's mesh. It checks that the new design matches the saved one
outside the region and at the interface. For the transmon, with a wider island:

```bash
julia --project -e 'include("transmon.jl"); using DeviceLayout: μm;
    generate_transmon(cap_width=30μm, mesh_filename="transmon_redesign.msh2",
                      config_filename="transmon_redesign.json")'
python3 transmon_substructuring.py --input mesh/transmon_redesign.msh2 \
    --output mesh/transmon_redesign_labeled.msh2
python3 ../substructuring/remesh_region.py --model mesh/transmon_substructuring.msh2 \
    --design mesh/transmon_redesign_labeled.msh2 \
    --output mesh/transmon_substructuring_redesign.msh2
palace transmon_substructuring_redesign.json  # reuses the saved model
```

The script assumes a region made of two materials separated by a layout plane (substrate and vacuum)
with planar features (metal and other boundaries) on that plane, and needs the `gmsh` and `numpy`
Python packages. It reads Gmsh 2.2 meshes only, so it does not apply to a model condensed from a
mesh saved by adaptive refinement (`"SaveAdaptMesh"` writes an MFEM mesh). The online run on the
redesigned mesh matches a full simulation on the same mesh. It differs from a simulation on the
layout tool's own mesh of the new design only by the change of mesh, as any two meshes of the same
geometry do.

## [Example: qubit lattice](@id substructuring-lattice-example)

The script
[`examples/substructuring/qubit_lattice.jl`](https://github.com/awslabs/palace/blob/main/examples/substructuring/qubit_lattice.jl)
generates the mesh of a ``5 \times 5`` lattice of simplified qubits with DeviceLayout.jl: each
qubit is a pair of rectangular pads in a cutout of the ground plane. The script
[`examples/substructuring/qubit_lattice_substructuring.py`](https://github.com/awslabs/palace/blob/main/examples/substructuring/qubit_lattice_substructuring.py)
splits the metal into the ground plane (a `"PEC"` boundary) and the 50 pads (the terminals),
and assigns the elements in a box around the middle qubit, about 4% of the mesh, to the region.
As for the transmon, it writes the configurations of the later runs from that of the offline
run, `qubit_lattice_offline.json`: they reuse its model, and the redesign and the runs on the
mesh of the region below take their meshes:

```bash
cd examples/substructuring
julia --project -e 'include("qubit_lattice.jl"); generate_qubit_lattice()'
python3 qubit_lattice_substructuring.py  # the mesh, and qubit_lattice_<run>.json
palace qubit_lattice_offline.json   # condense the environment, save the model
palace qubit_lattice_online.json    # reuse the model: only the region is solved
```

Both runs write the ``50 \times 50`` Maxwell capacitance matrix to `terminal-C.csv`. With many
terminals and a small region, the online run is much faster than a regular simulation of the
device: each terminal costs a solve in the region instead of one in the whole device. The middle
qubit can be redesigned as in the transmon example, for example with a wider pad gap:

```bash
julia --project -e 'include("qubit_lattice.jl"); using DeviceLayout: μm;
    generate_qubit_lattice(center_pad_gap=60μm, mesh_filename="qubit_lattice_redesign.msh2")'
python3 qubit_lattice_substructuring.py --input mesh/qubit_lattice_redesign.msh2 \
    --output mesh/qubit_lattice_redesign_labeled.msh2
python3 remesh_region.py --model mesh/qubit_lattice_substructuring.msh2 \
    --design mesh/qubit_lattice_redesign_labeled.msh2 \
    --output mesh/qubit_lattice_substructuring_redesign.msh2
palace qubit_lattice_redesign.json  # reuses the saved model
```

The offline run also writes the mesh of the region,
`postpro/qubit_lattice_offline/environment_region.mesh`, on which the online runs read only the
region (10,205 of the mesh's 249,307 elements). With `--region-mesh`, the re-meshing script writes
the redesign in the same way:

```bash
python3 remesh_region.py --model mesh/qubit_lattice_substructuring.msh2 \
    --design mesh/qubit_lattice_redesign_labeled.msh2 \
    --region-mesh postpro/qubit_lattice_offline/environment_region.mesh \
    --output mesh/qubit_lattice_redesign_region.msh2
palace qubit_lattice_online_region.json    # the original region, on the mesh of the region
palace qubit_lattice_redesign_region.json  # the redesign, on its mesh of the region
```

The region box stops above the bottom of the substrate, so that the region contains a single
substrate-vacuum interface, the layout plane, as the re-meshing script requires.

## [Example: flux-loop lattice](@id substructuring-flux-lattice-example)

The script
[`examples/substructuring/flux_lattice.jl`](https://github.com/awslabs/palace/blob/main/examples/substructuring/flux_lattice.jl)
generates the mesh of a ``5 \times 5`` lattice of superconducting rings with Gmsh: each ring is
an annular London film with a flux hole. The elements in a box around the middle ring are the
region. Each ring is a [`"FluxLoop"`](../config/reference.md#config-boundaries-fluxloop)
excitation, and all films are `"Superconductor"` boundaries:

```bash
cd examples/substructuring
julia --project -e 'include("flux_lattice.jl"); generate_flux_lattice()'
palace flux_lattice_offline.json   # condense the environment, save the model
palace flux_lattice_online.json    # reuse the model: only the region is solved
```

The script also writes the configuration of the online run, `flux_lattice_online.json`, from that
of the offline run, which it differs from by the model it reuses and its output folder. Both runs
write the ``25 \times 25`` inductance matrix to `terminal-M.csv`. The saved model
includes the condensed flux-loop excitations of the environment, so the online run needs no
environment solve.

The script
[`examples/substructuring/flux_sheet.jl`](https://github.com/awslabs/palace/blob/main/examples/substructuring/flux_sheet.jl)
generates a variant with one film: a square plate with a ``5 \times 5`` lattice of flux holes,
each hole a `"FluxLoop"` of the shared film. The film crosses the interface between the region
and the environment, and its screening currents couple the holes:

```bash
julia --project -e 'include("flux_sheet.jl"); generate_flux_sheet()'  # and flux_sheet_online.json
palace flux_sheet_offline.json
palace flux_sheet_online.json
```

## [Example: CPW resonator grid](@id substructuring-cpw-example)

The script
[`examples/substructuring/cpw_resonator_grid.jl`](https://github.com/awslabs/palace/blob/main/examples/substructuring/cpw_resonator_grid.jl)
generates the mesh of a ``3 \times 4`` grid of open-ended CPW resonators with DeviceLayout.jl:
the resonators of each row are coupled to a shared feedline with a lumped port at each end. The
script
[`examples/substructuring/cpw_resonator_grid_substructuring.py`](https://github.com/awslabs/palace/blob/main/examples/substructuring/cpw_resonator_grid_substructuring.py)
assigns the elements in a box around resonator 6, in the middle row, about 8% of the mesh, to the
region, and writes the configurations of the substructured runs from that of the regular one,
`cpw_resonator_grid.json`. The middle feedline crosses the interface, and all ports are in the
environment:

```bash
cd examples/substructuring
julia --project -e 'include("cpw_resonator_grid.jl"); generate_cpw_resonator_grid()'
python3 cpw_resonator_grid_substructuring.py  # the mesh, and cpw_resonator_grid_<run>.json
palace cpw_resonator_grid.json                # regular sweep, for comparison
palace cpw_resonator_grid_offline.json        # condense the environment, save the model
palace cpw_resonator_grid_online.json         # reuse the model: only the region is solved
palace cpw_resonator_grid_online_region.json  # the same on the mesh of the region
palace cpw_resonator_grid_redesign.json       # region substrate permittivity 0.4% higher
```

The sweep resolves the resonance of resonator 6 near 18.28 GHz in the transmission between the ports
of the middle feedline (`port-S.csv`), and the runs agree to the tolerance of the regular
simulation's linear solver. The redesign moves the resonance to 18.25 GHz, as a regular simulation
of the redesigned device does. Each online frequency factors only the region: on a 192-core node, an
online sweep is about 5 times faster than a regular one at the example's linear solver tolerance and
4 times at Palace's default, and about 6 times on the mesh of the region. The mesh has about 1.3
million unknowns, so the runs are best on a cluster node. The mesh, and so these frequencies, depend
on the versions of DeviceLayout.jl and Gmsh, which `examples/Manifest.toml` pins.

The script also writes `cpw_resonator_grid_offline_adaptive.json` and
`cpw_resonator_grid_online_adaptive.json`, which sample a rational model of the environment over
15–21.5 GHz, and sweep 1301 frequencies from it online, with the region reduced as in the adaptive
fast frequency sweep.
