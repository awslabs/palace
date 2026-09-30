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
extraction and [magnetostatic](../guide/problem.md#Magnetostatic-problems) inductance extraction with
flux-loop excitations.

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
regularization, see [Magnetostatics](#Magnetostatics)). Terminals and `"PEC"`
boundaries (grounded, as in a regular electrostatic simulation) may lie in the region, in the
environment, or span both; the environment's coupling to its terminals is condensed together
with ``\bm{S}_E``, so terminal energies need no environment solve.

The environment model (``\bm{S}_E`` and the terminal couplings) can be saved to a file. A later
run loads it and solves only the region, possibly after re-meshing the region: its solves then
follow the size of the region rather than of the device (the whole mesh is still read and the
environment operator assembled, to check that the model applies).

## Configuration

Substructuring is enabled with a
[`config["Solver"]["Substructuring"]`](../config/reference.md#config-solver-substructuring) block. The
region and the environment are given by domain attributes, which together must cover the mesh:

```json
"Substructuring":
{
  "Region": { "Attributes": [1, 2] },
  "Environment": { "Attributes": [10, 11] },
  "Mode": "Offline",
  "SaveModel": "postpro/environment.model"
}
```

  - `"Mode": "Offline"` condenses the environment and, if `"SaveModel"` is set, writes the
    environment model to that file.
  - `"Mode": "Online"` loads the model from `"SaveModel"` and solves only the region. The region
    may be re-meshed between runs, as long as the environment and the interface are unchanged.

The remaining settings (terminals, materials, `"Solver"/"Order"`, ...) are those of the regular
electrostatic or magnetostatic simulation. Field output follows
`"Solver"/"Electrostatic"/"Save"` and `"Solver"/"Magnetostatic"/"Save"`; the environment part
of a saved field is recovered only when fields are written.

### Choosing the region

Most meshes do not come with a region/environment split. Any split into two sets of domain
attributes works: the interface follows the element faces between them and need not be planar.
The cost of the condensation grows with the number of interface degrees of freedom, so a region
boundary in a coarsely meshed part of the domain (substrate or vacuum away from the metal) is
cheaper than one that cuts through a finely meshed area. The
[transmon example](@ref substructuring-transmon-example) retags the elements of an existing mesh by
the position of their centroids.

### Optional approximations

Two optional settings trade a controlled loss of accuracy for resources. Both default to `0`,
which keeps the condensation exact.

  - `"FactorizationTol"`: relative tolerance of a block low-rank (BLR) approximation of the
    environment factorization, which reduces the time to condense the environment. It is used
    when *Palace* is built with MUMPS. The error in the extracted matrices follows the tolerance
    but is not bounded by it; `1e-12` is a conservative choice.
  - `"InterfaceOffdiagTol"`: relative tolerance of a hierarchical (HODLR) compression of
    ``\bm{S}_E``, which reduces its memory for large interfaces. It does not reduce the time to
    condense the environment.

### Adaptive mesh refinement

The environment is condensed once, so its mesh resolution is fixed in the saved model. Meshes
from layout tools are often coarse and meant to be refined by
[adaptive mesh refinement](../reference.md#Error-estimation-and-adaptive-mesh-refinement-(AMR))
(AMR), so the refinement is done in two steps:

 1. Refine the full model before condensing it: run the electrostatic simulation without
    `"Substructuring"`, with AMR and
    [`config["Model"]["Refinement"]["SaveAdaptMesh"]`](../config/reference.md#config-model-refinement)
    set to `true`, then run the `"Offline"` condensation with the saved mesh as
    `config["Model"]["Mesh"]`. The refined elements keep the region and environment
    attributes.
 2. Refine the region only in `"Online"` runs, for example after a redesign: with
    `config["Model"]["Refinement"]["MaxIts"]` greater than zero, the error indicators are
    computed in the region and the environment is left as condensed. This requires
    nonconforming refinement without a level constraint (`"Nonconformal": true`,
    `"MaxNCLevels": 0`), so that the refinement does not spread into the environment.

Adaptive refinement in an `"Offline"` run is rejected (use step 1), and adaptive refinement with
substructuring is only available for electrostatics.

### Solvers

The environment and the region are factored with a sparse direct solver (SuperLU_DIST,
STRUMPACK or MUMPS, whichever Palace is built with) when they fit, and solved iteratively
otherwise. With MUMPS, the Schur complement of an electrostatic environment comes from a single
partial factorization. When running MUMPS with MPI, setting `OPENBLAS_NUM_THREADS=1` (or the
equivalent for the BLAS in use) avoids oversubscribing cores.

## Magnetostatics

Magnetostatic substructuring extracts the inductance matrix from
[`config["Boundaries"]["FluxLoop"]`](../config/reference.md#config-boundaries-fluxloop)
excitations, including the kinetic inductance of superconducting films. The `"PEC"` boundaries
are the Dirichlet boundaries of the condensation. The flux-loop films and the
[`config["Boundaries"]["Superconductor"]`](../config/reference.md#config-boundaries-superconductor)
films are London sheets, as in a regular simulation: a film without a `"Superconductor"` entry
has the penetration depth `"PecPenetrationDepth"`, and a `"Superconductor"` film its own. The
sheet term of each film face is condensed with the part of the domain the face bounds, so films
may lie in the region, in the environment, or cross the interface. A saved model includes the
condensed flux-loop excitations of the environment, so an online run needs no environment
solve as long as the excitations in the environment are unchanged. A small mass regularization
keeps the curl-curl operator definite; the extracted energies use the unregularized operator.

!!! note

    `"SurfaceCurrent"` excitations are not supported with substructuring yet: a configuration
    that uses them is rejected.

## [Example: transmon capacitance](@id substructuring-transmon-example)

The [transmon example](../examples/transmon.md) mesh contains a transmon island, a feedline, and
a ground plane with a readout resonator, all as a single metal boundary. The script
[`examples/transmon/transmon_substructuring.py`](https://github.com/awslabs/palace/blob/main/examples/transmon/transmon_substructuring.py)
splits the metal into its three conductors (the terminals) and assigns the elements in a box
around the qubit to the region:

```bash
cd examples/transmon
python3 transmon_substructuring.py
palace transmon_substructuring_offline.json   # condense the environment, save the model
palace transmon_substructuring_online.json    # reuse the model: only the region is solved
```

Both runs write the ``3 \times 3`` Maxwell capacitance matrix to `terminal-C.csv`; the online
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

The script assumes a region made of two materials separated by a layout plane (substrate and
vacuum) with planar features (metal and other boundaries) on that plane, and needs the `gmsh`
and `numpy` Python packages. The online run on the redesigned mesh matches a full simulation
on the same mesh. It differs from a simulation on the layout tool's own mesh of the new design
only by the change of mesh, as any two meshes of the same geometry do.
