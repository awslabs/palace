<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
# Apptainer/Singularity Container

This directory contains a [definition file](singularity.def) for building *Palace*
containers using [Apptainer/Singularity](https://apptainer.org/).

## Quick start

Assuming you have installed Apptainer/Singularity and the `singularity` executable is on
your path, the container can be built with:

```bash
singularity build palace.sif singularity.def
```

and run with:

```bash
singularity run palace.sif <ARGS...>
```

where `<ARGS...>` is a list of command line arguments provided to the `palace` executable.

For detailed instructions, see the documentation specific to
[building](https://awslabs.github.io/palace/dev/install/#Build-using-Singularity/Apptainer)
and [running](https://awslabs.github.io/palace/dev/run/#Singularity/Apptainer) *Palace*
with Apptainer/Singularity.
