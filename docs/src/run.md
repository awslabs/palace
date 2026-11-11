```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Running *Palace*

Once installed into a directory `<INSTALL_DIR>`, a parallel simulation using *Palace* can
be started with the following command:

```bash
<INSTALL_DIR>/bin/palace -np <NUM_PROCS> config.json
```

where

  - The installed [`palace`](https://github.com/awslabs/palace/blob/main/scripts/palace)
    script wraps a call to the desired MPI launcher (`mpirun` by default).
  - `<NUM_PROCS>` is the number of MPI processes to use for the simulation.
  - `config.json` is the JSON format configuration file used to specify the simulation
    parameters. The structure of this configuration file is outlined in detail in the
    section [Configuration File](config/config.md).

A full list of available script options is available using the `-h` or `--help` flag.

During the course of a simulation, the solver will write a number of useful statistics and
logging information to standard output. It is often helpful to save this information to a
file, for example with:

```bash
<INSTALL_DIR>/bin/palace ... | tee log.out
```

Of course, the interested user can explicitly run the *Palace* binary in parallel,
supplying options directly to their MPI launcher of choice, as:

```bash
<MPI_RUN> [OPTIONS] <INSTALL_DIR>/bin/palace-<ARCH>.bin config.json
```

where `<MPI_RUN>` is the MPI launcher command, `[OPTIONS]` is a list of command line options
passed to the MPI launcher, and `<ARCH>` is the machine architecture (`x86_64` or
`arm64`).

## Environment variables

The following environment variables, set for all MPI processes, adjust how *Palace* and its
sparse direct solver STRUMPACK run. The defaults suit most simulations, and the variables
are meant for troubleshooting and performance tuning: apart from the STRUMPACK settings
marked below, which change the factorization or its accuracy, they do not change the
results.

  - `PALACE_SETUP_THREADS`: the number of threads of the setup of the serial mesh on the
    processes which hold it: its preprocessing (the marking of tetrahedra for refinement and
    the cracking of interior boundaries) and, on the root process, the extraction of the
    parts of the other processes when it distributes the initial mesh. By default, they use
    one thread per MPI process of their node, on the cores of the other processes, which
    wait without using their CPUs meanwhile. Values below `1` count as `1`; the extraction
    uses at most one thread per MPI process, and a single one for meshes with high-order
    (curved) nodes. The threads which write meshes are not affected. Set it to the number of
    cores of the node when running fewer MPI processes per node than cores (for example one
    process per node with OpenMP threads), or to `1` to rule out the threads when diagnosing
    a problem.
  - `PALACE_GMSH_FAST`: `0` reads Gmsh meshes with MFEM's reader. By default, *Palace* reads
    meshes in Gmsh format 2.2 (ASCII or binary) with first-order elements and without
    periodicity with its own, faster reader, which builds the same mesh (the same order of
    the vertices, elements and boundary elements, and the same attributes) as MFEM's reader.
    Other mesh files, including other Gmsh format versions, compressed files, high-order
    elements and periodic meshes, are always read by MFEM. Set it to `0` if a Gmsh mesh
    fails to load, or to compare with MFEM's reader.

The STRUMPACK variables below exist only when STRUMPACK is built with the patches of
*Palace*, which the CMake superbuild and the Spack package of `palace@develop` apply.

  - `STRUMPACK_PTSCOTCH_NPROCS` (changes the factorization): the number of processes on
    which PT-Scotch computes the nested dissection separators of the parallel column
    ordering, after which METIS orders the subgraph of each process. Default: the number of
    MPI processes, at most 256; values are clamped to between 1 and the number of processes.
    It is used when STRUMPACK (`"Type": "STRUMPACK"`; `"STRUMPACK-MP"` ignores the ordering
    options) orders the matrix with PT-Scotch, by default whenever it is built with
    PT-Scotch (see
    [`config["Solver"]["Linear"]["ColumnOrdering"]`](config/reference.md#config-solver-linear-columnordering)).
    The time of PT-Scotch grows with the number of processes, while the subgraphs left to
    METIS grow as it decreases: lower it if PT-Scotch dominates the time of the ordering at
    large process counts, and raise it if the orderings of the subgraphs do.
  - `STRUMPACK_SHARED_PERM`: `0` stores the fill-reducing permutation and its inverse (two
    integers per matrix row) on every process. By default, without compression, the
    processes of a node share one copy in an MPI shared-memory window, after an ordering on
    a single process as after the parallel PT-Scotch ordering, which saves memory with many
    processes per node. Set it to `0` if the shared window cannot be allocated, for example
    in containers with a small `/dev/shm`.
  - `STRUMPACK_MEMORY_VERIFY`: a nonzero value checks, whenever the distributed elimination
    tree is built, the compact table which maps the matrix rows to the fronts of the
    distributed factorization against the full table, and prints the number of mismatches
    (zero expected), and checks the shared permutation (no invalid inverse entries, the same
    copy on every node). Off by default: it allocates the arrays of the global matrix size
    on every process that the compact table avoids.
  - `STRUMPACK_SI_MAX_COND` (changes the accuracy): the selective inversion is chosen by
    [`config["Solver"]["Linear"]["STRUMPACKSelectiveInversion"]`](config/reference.md#config-solver-linear-strumpackselectiveinversion),
    on by default for eigenmode, boundary mode and transient simulations and off for driven,
    electrostatic and magnetostatic ones, and has no effect with `"STRUMPACK-MP"`. With it,
    the solve phase replaces the triangular solves of a distributed separator block by a
    product with its explicit inverse only if the 1-norm condition number κ of the block,
    computed with the inverse, is at most this value, since a solve with the inverse is
    accurate to about u κ, with u the unit round-off; the processes of each front take the
    same decision. Default `9.0e10` (u κ ≤ `1e-5` in double precision). Lower values keep
    more triangular solves, which are slower; higher values invert more blocks, which can
    make the sparse direct solves of ill-conditioned matrices less accurate and increase the
    iterations of the Krylov solver they precondition. Any value below `1` skips the
    inversion entirely, so that `STRUMPACK_SI_MAX_COND=0` costs the same as
    `"STRUMPACKSelectiveInversion": false`.
  - `STRUMPACK_SI_VERBOSE`: a positive value prints, after each factorization, how many
    distributed fronts were inverted, the distribution of their condition numbers by decade,
    and up to 12 of the blocks kept with triangular solves. Off by default.
  - `STRUMPACK_PM_FRONT_COST`, `STRUMPACK_PM_DP` and `STRUMPACK_PM_MINMAX` (change the
    factorization only by round-off): the mapping of the elimination tree to the processes,
    which *Palace* bases on the memory of the factors to balance the solve phase.
    `STRUMPACK_PM_FRONT_COST` (default `1000`, at least `0`) is the cost added to each
    front, in factor entries, for the work of the solve phase per front, so that subtrees of
    many small fronts are mapped to enough processes. For the subtrees mapped to at most
    `STRUMPACK_PM_DP` processes (default `16`; values below `2` disable it), the split of
    the processes between the two children of each front minimizes the largest work per
    process of the whole subtree. For larger subtrees, the split minimizes the larger work
    per process of the two children; `STRUMPACK_PM_MINMAX=0` instead rounds their
    proportional shares (default `1`). They are meant for performance studies of the solve
    phase.

## Singularity/Apptainer

Assuming *Palace* was built using Singularity/Apptainer to `palace.sif`, running:

```bash
singularity run palace.sif <ARGS...>
```

corresponds to running a *Palace* simulation with command line arguments `<ARGS...>` using:

```bash
<INSTALL_DIR>/bin/palace <ARGS...>
```

as described above.
