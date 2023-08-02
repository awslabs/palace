```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Running *Palace*

Once installed into a directory `<INSTALL_DIR>`, a parallel simulation using *Palace* can
be started with the following command:

```bash
<INSTALL_DIR>/bin/palace -np <NUM_PROCS> config.json
```

where

  - The installed [`palace`]
    (https://github.com/awslabs/palace/blob/main/scripts/palace) script wraps
    a call to the desired MPI launcher (`mpirun` by default).
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
