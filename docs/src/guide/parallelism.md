```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Parallelism and GPU Support

*Palace* employs multiple types of parallelism in an attempt to maximize performance across
a wide range of deployment possibilities. The first is MPI-based distributed-memory
parallelism. This is controlled using the `-np` command line flag as outlined in
[Running *Palace*](../run.md).

Shared-memory parallelism using OpenMP is also available. To enable this, the
`-DPALACE_WITH_OPENMP=ON` option should be specified at configure time. At runtime, the
number of threads is configured with the `-nt` argument to the `palace` executable, or by
setting the [`OMP_NUM_THREADS`](https://www.openmp.org/spec-html/5.0/openmpse50.html)
environment variable.

Lastly, *Palace* supports GPU-acceleration using NVIDIA and AMD GPUs, activated with the
build options `-DPALACE_WITH_CUDA=ON` and `-DPALACE_WITH_HIP=ON`, respectively. At runtime,
the [`config["Solver"]["Device"]`](../config/solver.md#config%5B%22Solver%22%5D) parameter
in the configuration file can be set to `"CPU"` (the default) or `"GPU"` in order to
configure *Palace* and MFEM to use the available GPU device or devices. The
[`config["Solver"]["Backend"]`](../config/solver.md#config%5B%22Solver%22%5D) parameter, on
the other hand, controls the
[libCEED backend](https://libceed.org/en/latest/gettingstarted/#backends). Users typically
do not need to provide a value for this option and can instead rely on *Palace*'s default,
which selects the most appropriate backend for the given value of
[`config["Solver"]["Device"]`](../config/solver.md#config%5B%22Solver%22%5D).

In order to take full advantage of the performance benefits made available by GPU-
acceleration, it is recommended to make use of
[operator partial assembly](https://mfem.org/performance/), activated when the value of
[`config["Solver"]["PartialAssemblyOrder"]`](../config/solver.md#config%5B%22Solver%22%5D)
is less than [`config["Solver"]["Order"]`](../config/solver.md#config%5B%22Solver%22%5D).
This feature avoids assembling a global sparse matrix and instead makes use of data
structures for operators which lend themselves to more efficient asymptotic storage and
application costs. See also
[https://libceed.org/en/latest/intro/](https://libceed.org/en/latest/intro/) for more
details. Partial assembly in *Palace* supports mixed meshes including both tensor product
elements (hexahedra and quadrilaterals) as well as non-tensor product elements
(tetrahedra, prisms, pyramids, and triangles).
