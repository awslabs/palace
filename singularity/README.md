# _Palace_ Singularity / Apptainer container

This folder contains a [definition file](singularity.def) for building
containers using Singularity or Apptainer, which are widely used for HPC.


## Quick start

The container may be built with:

```bash
singularity build palace.sif singularity.def
```

and run with:

```bash
singularity run palace.sif [ARGS ...]
```

For detailed instructions, see [Installation – Build using Singularity / Apptainer](https://awslabs.github.io/palace/stable/install/#Build-using-Singularity-/-Apptainer) and
[Running _Palace_ – Singularity / Apptainer](https://awslabs.github.io/palace/stable/run/index.html#Singularity-/-Apptainer).