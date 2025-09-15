```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Single transmon with read-out resonator from [`DeviceLayout.jl`](https://aws-cqc.github.io/DeviceLayout.jl/)

!!! note
    
    The files for this example can be found in the
    [`examples/transmon/`](https://github.com/awslabs/palace/blob/main/examples/transmon)
    directory of the *Palace* source code. The example is run with `Order = 1` for computational efficiency. For best results, run with `Order = 2`.

In this example, we simulate a
[transmon](https://en.wikipedia.org/wiki/Transmon) coupled to a read-out
resonator with mesh generated with
[`DeviceLayout.jl`](https://aws-cqc.github.io/DeviceLayout.jl). This example
closely follows `DeviceLayout.jl`'s
[example](https://aws-cqc.github.io/DeviceLayout.jl/stable/examples/singletransmon/).

This example demonstrates the use of adaptive mesh refinement and anisotropic
materials in *Palace*.

### Overview

The system is depicted below. Physical details and properties can be customized
by passing keywords to `generate_transmon`.

We can see the transmon and the coplanar resonator. There are two ports at the
end of the line.

The substrate is sapphire, a material with low dielectric loss and anisotropic
properties. The properties were taken from reference [1]. The entire system is
in a box in vacuum. At the boundary of the box, we place `"Absorbing"` boundary
conditions to capture the radiative processes.

The metal lines are assumed to be perfect electric conductors, and we have three
lumped ports, two at the edges of the feedline, and one modeling a linearized
Jospheson junction.

We only focus on the first two eignemodes, which give us transmon and readout
resonator mode frequencies, decay rates, and corresponding electric and magnetic
field modes.

We are interested in running a `"Eigenmode"` simulation here to identify
operating frequencies and measure quality factors. The `DeviceLayout.jl`'s
[example](https://aws-cqc.github.io/DeviceLayout.jl/stable/examples/singletransmon/)
takes this one step further, and introduces an optimization loop to optimize the
transmon properties to match target design goals.

### Results

The simulation takes one minute or less and requires 3 GB of memory.

Note that this simulation was run with `Order = 1`. This greatly increases speed
(in this case, by a factor of 30), but severely affects accuracy. For best
results, run with `Order = 2`.

When we look at the log, we see the following. The problem is solved twice.
First, with the given mesh, second with an adaptively refined mesh.

We can plot the modes with ParaView. Below, we show the magnetic energy density.

<!-- Plot generated with transmon/plot_slice.py -->
```@raw html
<br/><p align="center">
  <img src="../../assets/examples/transmon-2.png" width="100%" />
</p><br/>

We see that the two modes indeed represent the transmon and the resonator modes.

Frequencies are in the `postpro/eig.csv` file.

The first column corresponds to the mode number, the second and third to the
real and imaginary parts of the eigenfrequency, and the forth to the computed
quality factor.

`postpro/port-EPR.csv` contains the energy participation ratios for the
inductive lumped elements (`LumpedPort`s where we specified the `L` parameter).
The only one in the JSON file was marked with `Index` 3, that's why we see
`p[3]` in the header.

The bulk dielectric loss is in `domain-E`.

Comparing `error-indicators.csv` in iteration 1 we can see that norm of the
error improved by more than 10 percent in one iteration (note that in this case
using Order 2 would lead to a 2x improvement in error with the same mesh).

## References

[1] Where the sapphire properties come from?
