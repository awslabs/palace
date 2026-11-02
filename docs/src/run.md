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

On Linux, *Palace* tunes two process settings at startup, each of which can be disabled by
setting an environment variable for all MPI processes:

  - When the kernel enables transparent huge pages only on request (`madvise` mode) and the
    C library is glibc 2.35 or newer, *Palace* restarts itself once with the glibc tunable
    `glibc.malloc.hugetlb=1`, so that its heap memory uses transparent huge pages. Set
    `PALACE_DISABLE_THP` to disable this.
  - When the MPI launcher binds the processes of a node to its first cores, one core each,
    and this leaves the NUMA domains of the node unevenly loaded, the processes are rebound
    to spread them evenly across the NUMA domains (as with `--map-by numa` for Open MPI).
    Set `PALACE_DISABLE_REBIND` to disable this.

The AMS preconditioner gathers the coarse levels of its auxiliary space multigrid
hierarchies onto fewer processes during the solve, and computes its auxiliary space
corrections concurrently. Set the environment variable `HYPRE_AMG_REDIST=0` for all MPI
processes to disable this.

With Open MPI on AWS Elastic Fabric Adapter (EFA) networks, exclude the OFI BTL (Open MPI
option `--mca btl ^ofi`, or `OMPI_MCA_btl=^ofi` in the environment): it opens additional
network domains, which can make the MPI initialization fail intermittently when many
processes per node are launched. The *Palace* containers set this by default.

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
