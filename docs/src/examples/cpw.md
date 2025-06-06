```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Crosstalk Between Coplanar Waveguides

!!! note
    
    The files for this example can be found in the
    [`examples/cpw/`](https://github.com/awslabs/palace/blob/main/examples/cpw)
    directory of the *Palace* source code.

In this example, we construct a frequency domain model to analyze the wave transmission,
reflection, near-end crosstalk, and far-end crosstalk for a four-port system comprised of
two side-by-side coplanar waveguides (CPW). Each CPW is characterized by a trace width
``w = 30\text{ μm}`` and gap width ``s = 18\text{ μm}``. The metal is modeled as an
infinitely thin, perfectly conducting boundary surface on top of a sapphire dielectric
substrate (parallel to C-axis: ``\varepsilon_r = 11.5``,
``\tan\delta = 8.6\times 10^{-5}``, perpendicular to C-axis: ``\varepsilon_r = 9.3``,
``\tan\delta = 3.0\times 10^{-5}``) of ``500\text{ μm}`` thickness with the
C-axis in the z-direction. This yields a characteristic impedance
``Z_0 = 56.02\text{ }\Omega`` for each of the lines [[1]](#References). The center-to-center
separating distance between the transmission lines on the substrate is ``266\text{ μm}``,
which means there is exactly ``200\text{ μm}`` of ground plane between them.

A visualization of the computational domain is shown below.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cpw-1.png" width="60%" />
</p><br/>
```

There are two different options for modeling the termination at the ends of the CPW:

  - Lumped port: A multielement uniform lumped port can be used to terminate the CPW by
    connecting the center conductor to the ground plane on each side with impedance
    ``Z = 2Z_0``.
  - Wave port: We can solve a 2D boundary eigenvalue problem for the mode shape and
    propagation constants for the characteristic CPW mode, and use this to terminate the
    transmission line.

Views of the mesh boundaries for these two configurations are shown below. In both cases the
computational domain is discretized using an unstructured tetrahedral mesh. The mesh files
are
[`mesh/cpw_wave_0.msh`](https://github.com/awslabs/palace/blob/main/examples/cpw/mesh/cpw_wave_0.msh)
and
[`mesh/cpw_lumped_0.msh`](https://github.com/awslabs/palace/blob/main/examples/cpw/mesh/cpw_lumped_0.msh),
respectively. In addition, this example includes two mesh files which include the thickness
of the metal trace:
[`mesh/cpw_wave.msh`](https://github.com/awslabs/palace/blob/main/examples/cpw/mesh/cpw_wave.msh)
and
[`mesh/cpw_lumped.msh`](https://github.com/awslabs/palace/blob/main/examples/cpw/mesh/cpw_lumped.msh).

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cpw-2.png" width="45%" />
  <img src="../../assets/examples/cpw-3.png" width="45%" />
</p><br/>
```

Likewise, there are two different options for how the system response is calculated over the
desired frequency band:

  - Uniform: Sample the frequency band with the full-fidelity model at equally spaced
    frequencies over the desired range.
  - Adaptive: Use the full-fidelity model to sample the solution at a few adaptively
    selected frequency points in the desired band, and then construct a low-cost surrogate
    model which is used to compute the response over the entire band.

This leads to four possible configurations, for which there are four configuration files in
the example directory:
[`cpw_lumped_uniform.json`](https://github.com/awslabs/palace/blob/main/examples/cpw/cpw_lumped_uniform.json),
[`cpw_lumped_adaptive.json`](https://github.com/awslabs/palace/blob/main/examples/cpw/cpw_lumped_adaptive.json),
[`cpw_wave_uniform.json`](https://github.com/awslabs/palace/blob/main/examples/cpw/cpw_wave_uniform.json),
and
[`cpw_wave_adaptive.json`](https://github.com/awslabs/palace/blob/main/examples/cpw/cpw_wave_adaptive.json).

The frequency response is computed for the band ``f\in[2.0,32.0]\text{ GHz}``. For the
uniform sweep, a step size of ``\Delta f=6.0\text{ GHz}`` is used, while the adaptive sweep
employs a much finer step size ``\Delta f=0.1\text{ GHz}``. Additionally both sweeps have an
explicit sample placed at ``17.0\text{ GHz}``. The adaptive fast frequency
sweep algorithm is given a tolerance of ``1\times10^{-3}`` for choosing the sampling
points; the simulation with uniform ports uses ``9`` frequency samples and that with wave
ports uses ``10``. Despite the much finer frequency resolution, the adaptive frequency
sweep simulations take roughly the same amount of time as the uniform ones where the
resulting resolution is worse by a factor of ``20``. Lastly, for all simulations, a
second-order finite element approximation for the solution is used.

The results from the four different simulations are presented in the plots below. Note that
here, ``\text{dB}`` means ``20\log_{10}(|S_{ij}|)``:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cpw-p2-11.png" width="70%" />
  <img src="../../assets/examples/cpw-p2-21.png" width="70%" />
  <img src="../../assets/examples/cpw-p2-31.png" width="70%" />
  <img src="../../assets/examples/cpw-p2-41.png" width="70%" />
</p><br/>
```

The first remark is that in both the lumped port and wave port cases, the adaptive fast
frequency sweep results are very close to the true solutions sampled by the uniform
sweeps.

Second, there is a discrepancy between the results using lumped ports and those with wave
ports, namely the lumped port excitation exhibits much higher reflection than for wave
ports. This is expected when using a lumped port to approximate the termination of a CPW,
and refining the mesh or increasing the order of the solution approximation leads to less
reflection. See below for the results with again ``p = 4`` for the order of the solution
space, effectively doubling the spatial resolution from ``p = 2``. For the adaptive solver
in these plots, we have also reduced the adaptive tolerance to ``1\times10^{-5}`` due to the
small value of ``|S_{41}|``.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cpw-p4-11.png" width="70%" />
  <img src="../../assets/examples/cpw-p4-21.png" width="70%" />
  <img src="../../assets/examples/cpw-p4-31.png" width="70%" />
  <img src="../../assets/examples/cpw-p4-41.png" width="70%" />
</p><br/>
```

!!! note
    
    The examples files for uniform sampling in `examples/cpw` actually specify excitations
    on two ports ("multi-excitation"). The two excitation are run in sequence during a
    single palace simulation.

## References

[1] H. J. Visser, _Antenna Theory and Applications_, Wiley, Hoboken, NJ, 2012.
