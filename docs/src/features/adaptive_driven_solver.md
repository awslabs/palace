```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

```@raw html
<!---
MAINTAINER NOTE

The SVGs referenced from this page (docs/src/assets/examples/driven_ua_*.svg)
are committed snapshots, not auto-generated artifacts. The underlying Palace
runs take several hours and are not feasible to re-run in CI or at docs build
time. The plots will drift relative to the current solver behaviour over time,
which is acceptable for a feature guide; if you make a change that meaningfully
alters the qualitative shape of these curves, please refresh
them by hand.

To regenerate, from the repository root:

  julia --project=examples -e 'include("examples/cpw/cpw_tutorial_lumped_driven.jl"); \
      generate_cpw_lumped_driven_data(num_processors=4)'
  julia --project=examples -e 'include("examples/transmon/transmon_tutorial_driven.jl"); \
      generate_transmon_driven_data(num_processors=4)'
  julia --project=examples examples/cpw/cpw_tutorial_lumped_driven_plots.jl
  julia --project=examples examples/transmon/transmon_tutorial_driven_plots.jl

then commit the updated docs/src/assets/examples/driven_ua_*.svg.
--->
```

# Adaptive Frequency Sweeps for Driven Simulations

*Palace* has two modes for frequency-domain driven simulations: uniform and adaptive. The
uniform solver runs one full finite-element solve per output frequency. The adaptive solver
runs a much smaller number of full solves, constructs a projection-based reduced-order model
(PROM), and evaluates that model over the requested output grid.

The adaptive solver can be much faster than the uniform solver on fine frequency grids, but
it introduces an approximation error that should be controlled and validated for the
quantities of interest. The implementation details are described in the
[adaptive driven solver reference](../reference.md#Adaptive-driven-solver-and-reduced-order-modeling).
This page focuses on practical setup and validation.

!!! warning "Advanced feature"

    Adaptive driven sweeps rely on algorithmic choices inside *Palace*. They are useful and
    efficient, but users should validate them against uniform sweeps when studying a new
    model or a new output quantity.

## Quick start

The only configuration difference between uniform and adaptive driven sweeps is the value of
[`config["Solver"]["Driven"]["AdaptiveTol"]`](../config/reference.md#config-solver-driven).
The default value is `0.0`, which uses the uniform solver. Any positive value enables the
adaptive solver.

```json
"Solver": {
  "Driven": {
    "Samples": [ {"Type": "Linear", "MinFreq": 3.5, "MaxFreq": 7.0, "FreqStep": 0.1} ],
    "AdaptiveTol": 1e-3
  }
}
```

The `"Samples"` specification defines the output frequency grid for both solvers. For the
adaptive solver this grid can usually be fine, since adding output frequencies is much
cheaper than adding full finite-element solves. The output files, such as `domain-E.csv` and
`port-S.csv`, have the same format in both modes.

When studying a new model, a good workflow is:

 1. Run a uniform sweep on a coarse output grid or at key frequencies.
 2. Run the adaptive solver on the same grid and compare the quantities you care about.
 3. Tighten tolerances, add validation frequencies, or fall back to the uniform solver if the
    adaptive result is not accurate enough for the application.
 4. Once validated, run the adaptive solver on the desired fine output grid.

## Tuning parameters

The most important configuration options are:

  - [`config["Solver"]["Driven"]["AdaptiveTol"]`](../config/reference.md#config-solver-driven):
    adaptive convergence tolerance. A value around `1e-3` is often a reasonable starting
    point for S-parameter sweeps, but the right value is model- and quantity-dependent.
  - [`config["Solver"]["Linear"]["Tol"]`](../config/reference.md#config-solver-linear):
    linear solver tolerance for the underlying full solves. This should usually be
    substantially smaller than `"AdaptiveTol"`. If *Palace* logs warnings about a
    rank-deficient minimal rational interpolation, try tightening this tolerance or using a
    looser `"AdaptiveTol"`.
  - [`config["Solver"]["Driven"]["AdaptiveMaxSamples"]`](../config/reference.md#config-solver-driven):
    maximum number of full samples per excitation. If this cap is reached before
    convergence, *Palace* warns that the outputs may not be accurate to `"AdaptiveTol"`:
    increase it or revisit the tolerances. If *Palace* warns that the PROM did not reproduce
    some of the samples even after their iterative refinement, `"AdaptiveTol"` is likely
    below the accuracy of the full solves: loosen it.
  - [`config["Solver"]["Driven"]["AdaptiveConvergenceMemory"]`](../config/reference.md#config-solver-driven):
    number of consecutive successful sample checks required before convergence. Increasing
    this can make early termination less likely, at the cost of more full solves.

`"AdaptiveTol"` controls an error criterion on the finite-element electric-field solution,
not a strict relative-error bound on every derived output. For example, a scattering
parameter that is nearly zero can have a large pointwise relative error even when the
absolute error is small. Validate the derived quantities that matter for your analysis.

## Fields at every output frequency

Writing fields with `"Save"` costs one full field file per saved frequency. The reduced-order
model already contains the field at every frequency of the sweep:
``\bm{E}(\omega) = \sum_k y_k(\omega)\,\bm{E}_k``, with ``n`` real basis fields
``\bm{E}_k`` and complex coefficients ``y_k(\omega)``. With
[`config["Solver"]["Driven"]["AdaptiveSaveBasis"]`](../config/reference.md#config-solver-driven)
set to `true`, *Palace* writes:

  - the basis fields ``\bm{E}_k`` and their curls ``\nabla\times\bm{E}_k`` as MFEM grid
    functions in `gridfunction/driven/basis/` (`E_<k>.gf.<rank>` and `curlE_<k>.gf.<rank>`,
    in V/m and V/m², at full double precision), with the mesh files `mesh.<rank>` of the
    same run, in the input length unit. The grid functions of each MPI rank go with the
    `mesh.<rank>` of that rank in `basis/`;
  - the coefficients at every output frequency in `rom-basis-coefficients.csv`, one row per
    frequency, with columns `Re{y[k][e]}` and `Im{y[k][e]}` for basis vector `k` and
    excitation `e`.

The magnetic flux density follows from the same coefficients: with the ``e^{i\omega t}``
time convention of *Palace*, ``\nabla\times\bm{E} = -i\omega\bm{B}``, so
``\bm{B}(\omega) = \frac{i}{\omega}\sum_k y_k(\omega)\,\nabla\times\bm{E}_k``, in T for
``\omega`` in rad/s. In 2D, the curl is the scalar ``z`` component. The option is not
available with a Floquet wave vector, for which ``\bm{B}`` has an additional
``\bm{k}\times\bm{E}`` term.

The storage grows with ``n``, not with the number of output frequencies, and the online
sweep is cheap, so a fine frequency grid costs little. The fields at one frequency are
linear combinations of degrees of freedom, so they need only NumPy (the files of one MPI
rank shown; repeat for each rank):

```python
import numpy as np

def read_gf(path):  # MFEM grid function: header, blank line, degrees of freedom
    header, values = open(path).read().split("\n\n", 1)
    return header, np.array(values.split(), dtype=float)

csv = "postpro/rom-basis-coefficients.csv"
names = [s.strip() for s in open(csv).readline().split(",")]
c = np.loadtxt(csv, delimiter=",", skiprows=1, ndmin=2)
i = np.argmin(np.abs(c[:, 0] - 5.0))  # row at 5 GHz
e = 1  # excitation
n = sum(s.startswith("Re{y[") and s.endswith(f"][{e}]}}") for s in names)
y = [c[i, names.index(f"Re{{y[{k}][{e}]}}")] + 1j * c[i, names.index(f"Im{{y[{k}][{e}]}}")]
     for k in range(1, n + 1)]
basis = "postpro/gridfunction/driven/basis/{}_{:06d}.gf.000000"
E = sum(yk * read_gf(basis.format("E", k))[1] for k, yk in enumerate(y, start=1))
omega = 2 * np.pi * 1e9 * c[i, 0]  # rad/s
B = 1j / omega * sum(yk * read_gf(basis.format("curlE", k))[1]
                     for k, yk in enumerate(y, start=1))
```

The real and imaginary parts can be written back with the same headers for visualization
with `mesh.000000`, or evaluated at points with PyMFEM. The accuracy is that of the
reduced-order model, set by `"AdaptiveTol"`.

## Coplanar waveguide example

The files for the CPW part of this guide are in
[`examples/cpw/`](https://github.com/awslabs/palace/blob/main/examples/cpw):

  - `cpw_tutorial_lumped_uniform.json`: uniform reference sweep.
  - `cpw_tutorial_lumped_adaptive.json`: adaptive sweep template.
  - `cpw_tutorial_lumped_driven.jl`: runs the uniform and adaptive simulations.
  - `cpw_tutorial_lumped_driven_plots.jl`: regenerates the plots below from the simulation
    outputs.

This example revisits the lumped-port version of the
[CPW crosstalk example](../examples/cpw.md). The uniform reference performs a full solve at
each output frequency, while the adaptive sweeps use the same output grid but construct a
PROM from adaptively selected full solves.

The uniform solver gives the following total electric energy and S-parameter magnitude:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_energy_uniform.svg" width="80%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_sparam_uniform.svg" width="80%" />
</p><br/>
```

For the adaptive solver, the diamond markers show the internal sample frequencies used to
build the PROM. These samples need not coincide with the output grid, except for required
boundary or user-forced samples.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_energy_adaptive_single.svg" width="80%" />
</p><br/>
```

As `"AdaptiveTol"` is tightened, *Palace* adds more internal samples and the error against
the uniform reference generally decreases.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_energy_adaptive_sweep.svg" width="80%" />
</p><br/>
```

The pointwise relative error for S-parameters can look large when the reference value is
small. In the CPW example, some scattering parameters are tiny over much of the band, so
relative error alone can exaggerate visually small absolute differences.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_sparam_adaptive_pointwise.svg" width="80%" />
</p><br/>
```

A complementary view is the absolute error normalized by an RMS scale for the reference
S-parameter over the frequency range:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_sparam_adaptive_rms.svg" width="80%" />
</p><br/>
```

The adaptive error indicator does not have to decrease monotonically. The convergence memory
parameter exists because an estimator can occasionally look converged before the next samples
would reveal a larger error.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_adaptive_convergence_curve.svg" width="80%" />
</p><br/>
```

## Transmon example

The transmon part of this guide uses the same geometry as the
[transmon eigenmode example](../examples/transmon.md), but runs driven simulations by
exciting the two resistive feedline ports. The files are in
[`examples/transmon/`](https://github.com/awslabs/palace/blob/main/examples/transmon):

  - `transmon_tutorial_driven.json`: driven simulation template.
  - `transmon_tutorial_driven.jl`: runs the uniform and adaptive simulations.
  - `transmon_tutorial_driven_plots.jl`: regenerates the plots below.

The model contains resonant modes near the driven frequency band, so the response is more
structured than in the CPW example. This is a useful case for adaptive sampling: the PROM can
place full solves near the frequencies that most influence the response.

The uniform driven sweep resolves the electric energy and the feedline S-parameters:

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_transmon_energy_uniform.svg" width="80%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_transmon_sparam_uniform.svg" width="80%" />
</p><br/>
```

The adaptive error is largest near resonant features and decreases as the tolerance is
tightened. The sample locations become more structured than in the CPW example because the
response is influenced by nearby poles.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_transmon_energy_adaptive_sweep.svg" width="80%" />
</p><br/>
```

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_transmon_sparam_adaptive_rms.svg" width="80%" />
</p><br/>
```

The uniform reference can also lose accuracy near high-``Q`` poles because the full system is
more poorly conditioned there. Tightening only the adaptive tolerance is not useful if the
underlying full solves are themselves inaccurate.

## Regenerating the plots

The committed plots are snapshots. To regenerate them, first run the simulations, then run
the plot scripts:

```bash
julia --project=examples -e 'include("examples/cpw/cpw_tutorial_lumped_driven.jl"); generate_cpw_lumped_driven_data(num_processors=4)'
julia --project=examples -e 'include("examples/transmon/transmon_tutorial_driven.jl"); generate_transmon_driven_data(num_processors=4)'
julia --project=examples examples/cpw/cpw_tutorial_lumped_driven_plots.jl
julia --project=examples examples/transmon/transmon_tutorial_driven_plots.jl
```

The scripts write SVGs to `docs/src/assets/examples/` by default.
