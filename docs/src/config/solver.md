```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# `config["Solver"]`

```json
"Solver":
{
    "Order": <int>,
    "PartialAssemblyOrder": <int>,
    "Device": <string>,
    "Backend": <string>,
    "Eigenmode":
    {
        ...
    },
    "Driven":
    {
        ...
    },
    "Transient":
    {
        ...
    },
    "Electrostatic":
    {
        ...
    },
    "Magnetostatic":
    {
        ...
    },
    "Linear":
    {
        ...
    }
}
```

with

`"Order" [1]` :  Finite element order (degree). Arbitrary high-order spaces are supported.

`"PartialAssemblyOrder" [1]` :  Order at which to switch from full assembly of finite
element operators to [partial assembly](https://mfem.org/howto/assembly_levels/). Setting
this parameter equal to 1 will fully activate operator partial assembly on all levels, while
setting it to some large number (greater than the finite element order) will result in
fully assembled operators as sparse matrices.

`"Device" ["CPU"]` :  The runtime device configuration passed to
[MFEM](https://mfem.org/howto/assembly_levels/) in order to activate different options
specified during configuration. The available options are:

  - `"CPU"`
  - `"GPU"`
  - `"Debug"`

The `"GPU"` option will automatically activate the `cuda` or `hip` device based on whether
MFEM is built with CUDA (`MFEM_USE_CUDA=ON`) or HIP (`MFEM_USE_HIP=ON`) support. When
*Palace* is built with OpenMP support (`PALACE_WITH_OPENMP=ON`), `omp` is automatically
added to the list of activated MFEM devices. The `"Debug"` option for MFEM's `debug` device
is useful for debugging issues associated with GPU-based runs of *Palace*.

`"Backend" [""]` :  Specifies the
[libCEED backend](https://libceed.org/en/latest/gettingstarted/#backends) to use for the
simulation. If no backend is specified, a suitable default backend is selected based on the
given `config["Solver"]["Device"]`.

`"Eigenmode"` :  Top-level object for configuring the eigenvalue solver for the eigenmode
simulation type. Thus, this object is only relevant for
[`config["Problem"]["Type"]: "Eigenmode"`](problem.md#config%5B%22Problem%22%5D).

`"Driven"` :  Top-level object for configuring the frequency domain driven simulation type.
Thus, this object is only relevant for
[`config["Problem"]["Type"]: "Driven"`](problem.md#config%5B%22Problem%22%5D).

`"Transient"` :  Top-level object for configuring the time domain driven simulation type.
Thus, this object is only relevant for
[`config["Problem"]["Type"]: "Transient"`](problem.md#config%5B%22Problem%22%5D).

`"Electrostatic"` :  Top-level object for configuring the electrostatic simulation type.
Thus, this object is only relevant for
[`config["Problem"]["Type"]: "Electrostatic"`](problem.md#config%5B%22Problem%22%5D).

`"Magnetostatic"` :  Top-level object for configuring the magnetostatic simulation type.
Thus, this object is only relevant for
[`config["Problem"]["Type"]: "Magnetostatic"`](problem.md#config%5B%22Problem%22%5D).

`"Linear"` :  Top-level object for configuring the linear solver employed by all simulation
types.

### Advanced solver options

  - `"QuadratureOrderJacobian" [false]`
  - `"ExtraQuadratureOrder" [0]`

## `solver["Eigenmode"]`

```json
"Eigenmode":
{
    "Target": <float>,
    "Tol": <float>,
    "MaxIts": <int>,
    "MaxSize": <int>,
    "N": <int>,
    "Save": <int>,
    "Type": <int>,
    "ContourTargetUpper": <float>,
    "ContourAspectRatio": <float>,
    "ContourNPoints": <int>
}
```

with

`"Target" [None]` :  (Nonzero) frequency target above which to search for eigenvalues, GHz.

`"Tol" [1.0e-6]` :  Relative convergence tolerance for the eigenvalue solver.

`"MaxIts" [0]` :  Maximum number of iterations for the iterative eigenvalue solver. A value
less than 1 uses the solver default.

`"MaxSize" [0]` :  Maximum subspace dimension for eigenvalue solver. A value less than 1
uses the solver default.

`"N" [1]` :  Number of eigenvalues to compute.

`"Save" [0]` :  Number of computed field modes to save to disk for
[visualization with ParaView](../guide/postprocessing.md#Visualization). Files are saved in
the `paraview/` directory under the directory specified by
[`config["Problem"]["Output"]`](problem.md#config%5B%22Problem%22%5D).

`"Type" ["Default"]` :  Specifies the eigenvalue solver to be used in computing the given
number of eigenmodes of the problem. The available options are:

  - `"SLEPc"`
  - `"ARPACK"`
  - `"Default"` :  Use the default eigensolver. Currently, this is the Krylov-Schur
    eigenvalue solver from `"SLEPc"`.

### Advanced eigenmode solver options

  - `"PEPLinear" [true]`
  - `"Scaling" [true]`
  - `"StartVector" [true]`
  - `"StartVectorConstant" [false]`
  - `"MassOrthogonal" [false]`

## `solver["Driven"]`

```json
"Driven":
{
    "MinFreq": <float>,
    "MaxFreq": <float>,
    "FreqStep": <float>,
    "SaveStep": <int>,
    "Restart": <int>,
    "AdaptiveTol": <float>,
    "AdaptiveMaxSamples": <int>,
    "AdaptiveConvergenceMemory": <int>
}
```

with

`"MinFreq" [None]` :  Lower bound of frequency sweep interval, GHz.

`"MaxFreq" [None]` :  Upper bound of frequency sweep interval, GHz.

`"FreqStep" [None]` :  Frequency step size for frequency sweep, GHz.

`"SaveStep" [0]` :  Controls how often, in number of frequency steps, to save computed
fields to disk for [visualization with ParaView](../guide/postprocessing.md#Visualization).
Files are saved in the `paraview/` directory under the directory specified by
[`config["Problem"]["Output"]`](problem.md#config%5B%22Problem%22%5D).

`"Restart" [1]` :  Iteration (1-based) from which to restart for a partial frequency sweep
simulation. That is, the initial frequency will be computed as
`"MinFreq" + ("Restart" - 1) * "FreqStep"`.

`"AdaptiveTol" [0.0]` :  Relative error convergence tolerance for adaptive frequency sweep.
If zero, adaptive frequency sweep is disabled and the full-order model is solved at each
frequency step in the specified interval. If positive, this tolerance is used to ensure the
reliability of the reduced-order model relative to the full-order one in the frequency band
of interest.

`"AdaptiveMaxSamples" [20]` :  Maximum number of frequency samples used to construct the
reduced-order model for adaptive fast frequency sweep, if the specified tolerance
(`"AdaptiveTol"`) is not met first.

`"AdaptiveConvergenceMemory" [2]` :  Memory used for assessing convergence of the adaptive
sampling algorithm for constructing the reduced-order model for adaptive fast frequency
sweep. For example, a memory of "2" requires two consecutive samples which satisfy the
error tolerance.

## `solver["Transient"]`

```json
"Transient":
{
    "Type": <string>,
    "Excitation": <string>,
    "ExcitationFreq": <float>,
    "ExcitationWidth": <float>,
    "MaxTime": <float>,
    "TimeStep": <float>,
    "SaveStep": <int>,
    "Order": <int>,
    "RelTol": <float>,
    "AbsTol": <float>
}
```

with

`"Type" ["Default"]` :  Specifies the time integration scheme used for the discretization of
the second-order system of differential equations. The available options are:

  - `"GeneralizedAlpha"` :  The second-order implicit generalized-``\alpha`` method with
    ``\rho_{\inf} = 1.0``. This scheme is unconditionally stable.
  - `"ARKODE"` :  SUNDIALS ARKode implicit Runge-Kutta scheme applied to the first-order
    ODE system for the electric field with adaptive time-stepping. This option is only available when *Palace* has been [built with SUNDIALS support](../install.md#Configuration-options).
  - `"CVODE"` :  SUNDIALS CVODE implicit multistep method scheme applied to the first-order
    ODE system for the electric field with adaptive time-stepping. This option is only available when *Palace* has been [built with SUNDIALS support](../install.md#Configuration-options).
  - `"RungeKutta"` : Two stage, singly diagonal implicit Runge-Kutta (SDIRK) method. Second order and L-stable.
  - `"Default"` :  Use the default `"GeneralizedAlpha"` time integration scheme.

`"Excitation" [None]` :  Controls the time dependence of the source excitation. The
available options are:

  - `"Sinusoidal"` :  A sinusoidal excitation at a user specified frequency.
  - `"Gaussian"` :  A Gaussian pulse with a user specified width which defines the
    bandwidth.
  - `"DifferentiatedGaussian"` :  A differentiated Gaussian pulse with a user specified
    width which defines the bandwidth.
  - `"ModulatedGaussian"` :  A modulated Gaussian pulse at a user specified center frequency
    and width used to excite the system without any DC component.
  - `"Ramp"` :  A differentiable unit step function to model the ramp up to a DC signal.
  - `"SmoothStep"` :  A smoother many-times differentiable unit step function to model the
    ramp up to a DC signal over a specified width of time.

`"ExcitationFreq" [None]` :  Center frequency used for harmonic source excitations, GHz.
Only relevant when `"Excitation"` is one of `"Sinusoidal"`, `"Gaussian"`,
`"DifferentiatedGaussian"`, or `"ModulatedGaussian"`.

`"ExcitationWidth" [None]` :  Pulse width for Gaussian-type source excitations, ns. Only
relevant when `"Excitation"` is one of `"Gaussian"`, `"DifferentiatedGaussian"`,
`"ModulatedGaussian"`, or `"SmoothStep"`.

`"MaxTime" [None]` :  End of simulation time interval, ns. Transient simulations always
start from rest at ``t = 0.0``.

`"TimeStep" [None]` :  Uniform time step size for time integration, ns.

`"SaveStep" [0]` :  Controls how often, in number of time steps, to save computed fields to
disk for [visualization with ParaView](../guide/postprocessing.md#Visualization). Files are
saved in the `paraview/` directory under the directory specified by
[`config["Problem"]["Output"]`](problem.md#config%5B%22Problem%22%5D).

`"Order" [2]` :  Order of the adaptive Runge-Kutta integrators or maximum order of the
multistep method, must be within `[2,5]`. Should only be specified if `"Type"` is `"ARKODE"`
or `"CVODE"`.

`"RelTol" [1e-4]` :  Relative tolerance used in adaptive time-stepping schemes. Should only
be specified if `"Type"` is `"ARKODE"` or `"CVODE"`.

`"AbsTol" [1e-9]` :  Absolute tolerance used in adaptive time-stepping schemes. Should only
be specified if `"Type"` is `"ARKODE"` or `"CVODE"`.

## `solver["Electrostatic"]`

```json
"Electrostatic":
{
    "Save": <int>
}
```

with

`"Save" [0]` :  Number of computed electric field solutions to save to disk for
[visualization with ParaView](../guide/postprocessing.md#Visualization), ordered by the
entries in the computed capacitance matrix. Files are saved in the `paraview/` directory
under the directory specified by
[`config["Problem"]["Output"]`](problem.md#config%5B%22Problem%22%5D).

## `solver["Magnetostatic"]`

```json
"Magnetostatic":
{
    "Save": <int>
}
```

with

`"Save" [0]` :  Number of computed magnetic field solutions to save to disk for
[visualization with ParaView](../guide/postprocessing.md#Visualization)), ordered by the
entries in the computed inductance matrix. Files are saved in the `paraview/` directory
under the directory specified by
[`config["Problem"]["Output"]`](problem.md#config%5B%22Problem%22%5D).

## `solver["Linear"]`

```json
"Linear":
{
    "Type": <string>,
    "KSPType": <string>,
    "Tol": <float>,
    "MaxIts": <int>,
    "MaxSize": <int>,
    "MGMaxLevels": <int>,
    "MGCoarsenType": <string>,
    "MGCycleIts": <int>,
    "MGSmoothIts": <int>,
    "MGSmoothOrder": <int>,
    "PCMatReal": <bool>,
    "PCMatShifted": <bool>,
    "ComplexCoarseSolve": <bool>,
    "PCSide": <string>,
    "DivFreeTol": <float>,
    "DivFreeMaxIts": <float>,
    "EstimatorTol": <float>,
    "EstimatorMaxIts": <float>,
    "EstimatorMG": <bool>,
    "GSOrthogonalization": <string>
}
```

with

`"Type" ["Default"]` :  Specifies the solver used for
[preconditioning](https://en.wikipedia.org/wiki/Preconditioner) the linear system of
equations to be solved for each simulation type. The available options are:

  - `"SuperLU"` :  The [SuperLU_DIST](https://github.com/xiaoyeli/superlu_dist) sparse
    direct solver in real double precision is used to factorize the system matrix. For
    frequency domain problems this uses a real approximation to the true complex linear
    system matrix. This option is only available when *Palace* has been
    [built with SuperLU_DIST support](../install.md#Configuration-options).
  - `"STRUMPACK"` :  The [STRUMPACK](https://portal.nersc.gov/project/sparse/strumpack)
    sparse direct solver in real double precision is used to factorize the system matrix.
    For frequency domain problems this uses a real approximation to the true complex linear
    system matrix. This option is only available when *Palace* has been
    [built with STRUMPACK support](../install.md#Configuration-options).
  - `"MUMPS"` :  The [MUMPS](http://mumps.enseeiht.fr/) sparse direct solver in real double
    precision is used to factorize the system matrix. For frequency domain problems this
    uses a real approximation to the true complex linear system matrix. This option is only
    available when *Palace* has been
    [built with MUMPS support](../install.md#Configuration-options).
  - `"AMS"` :  Hypre's
    [Auxiliary-space Maxwell Solver (AMS)](https://hypre.readthedocs.io/en/latest/solvers-ams.html),
    an algebraic multigrid (AMG)-based preconditioner.
  - `"BoomerAMG"` :  The
    [BoomerAMG](https://hypre.readthedocs.io/en/latest/solvers-boomeramg.html) AMG solver
    from Hypre.
  - `"Jacobi"` :  Diagonal scaling with a simple Jacobi preconditioner (not recommended in
    general).
  - `"Default"` :  Use the default `"AMS"` solver for simulation types involving definite or
    semi-definite curl-curl operators (time domain problems as well as magnetostatics). For
    frequency domain problems, use a sparse direct solver if available, otherwise uses
    `"AMS"`. For electrostatic problems, uses `"BoomerAMG"`.

`"KSPType" ["Default"]` :  Specifies the iterative
[Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace) solver type for solving
linear systems of equations arising for each simulation type. The available options are:

  - `"CG"`
  - `"GMRES"`
  - `"FGMRES"`
  - `"Default"` :  Use the default `"GMRES"` Krylov subspace solver for frequency domain
    problems, that is when
    [`config["Problem"]["Type"]`](problem.md#config%5B%22Problem%22%5D) is `"Eigenmode"` or
    `"Driven"`. For the other simulation types, the linear system matrix is always real and
    symmetric positive definite (SPD) and the preconditioned conjugate gradient method
    (`"CG"`) is used as the Krylov solver.

`"Tol" [1.0e-6]` :  Relative residual convergence tolerance for the iterative linear solver.

`"MaxIts" [100]` :  Maximum number of iterations for the iterative linear solver.

`"MaxSize" [0]` :  Maximum Krylov space size for the GMRES and FGMRES solvers. A value less
than 1 defaults to the value specified by `"MaxIts"`.

`"MGMaxLevels" [100]` :  Chose whether to enable
[geometric multigrid preconditioning](https://en.wikipedia.org/wiki/Multigrid_method) which
uses p- and h-multigrid coarsening as available to construct the multigrid hierarchy. The
solver specified by `"Type"` is used on the coarsest level. Relaxation on the fine levels
is performed with Chebyshev smoothing.

`"MGCoarsenType" ["Logarithmic"]` :  Coarsening to create p-multigrid levels.

  - `"Logarithmic"`
  - `"Linear"`

`"MGCycleIts" [1]` :  Number of V-cycle iterations per preconditioner application for
multigrid preconditioners (when `"UseMultigrid"` is `true` or `"Type"` is `"AMS"` or
`"BoomerAMG"`).

`"MGSmoothIts" [1]` :  Number of pre- and post-smooth iterations used for multigrid
preconditioners (when `"UseMultigrid"` is `true` or `"Type"` is `"AMS"` or `"BoomerAMG"`).

`"MGSmoothOrder" [0]` :  Order of polynomial smoothing for geometric multigrid
preconditioning (when `"UseMultigrid"` is `true`). A value less than 1 defaults to twice
the solution order given in
[`config["Solver"]["Order"]`](problem.md#config%5B%22Solver%22%5D) or 4, whichever is
larger.

`"PCMatReal" [false]` :  When set to `true`, constructs the preconditioner for frequency
domain problems using a real-valued approximation of the system matrix. This is always
performed for the coarsest multigrid level regardless of the setting of `"PCMatReal"`.

`"PCMatShifted" [false]` :  When set to `true`, constructs the preconditioner for frequency
domain problems using a positive definite approximation of the system matrix by flipping
the sign for the mass matrix contribution, which can help performance at high frequencies
(relative to the lowest nonzero eigenfrequencies of the model).

`"ComplexCoarseSolve" [false]` : When set to `true`, the coarse-level solver uses the true
complex-valued system matrix. When set to `false`, the real-valued approximation is used.

`"PCSide" ["Default"]` :  Side for preconditioning. Not all options are available for all
iterative solver choices, and the default choice depends on the iterative solver used.

  - `"Left"`
  - `"Right"`
  - `"Default"`

`"DivFreeTol" [1.0e-12]` :  Relative tolerance for divergence-free cleaning used in the
eigenmode simulation type. Ignored if non-zero Floquet wave vector is specified in
[`config["Boundaries"]["Periodic"]["FloquetWaveVector"]`](boundaries.md##boundaries%5B%%22Periodic%22%5D%22FloquetWaveVector%22%5D)
or
[`config["Boundaries"]["FloquetWaveVector"]`](boundaries.md##boundaries%5B%%22FloquetWaveVector%22%5D),
or non-zero
[`config["Domains"]["Materials"]["LondonDepth"]`](domains.md##domains%5B%22Materials%22%5D%5B%22LondonDepth%22%5D)
is specified.

`"DivFreeMaxIts" [1000]` :  Maximum number of iterations for divergence-free cleaning use in
the eigenmode simulation type. Ignored if non-zero Floquet wave vector is specified in
[`config["Boundaries"]["Periodic"]["FloquetWaveVector"]`](boundaries.md##boundaries%5B%%22Periodic%22%5D%22FloquetWaveVector%22%5D)
or
[`config["Boundaries"]["FloquetWaveVector"]`](boundaries.md##boundaries%5B%%22FloquetWaveVector%22%5D),
or non-zero
[`config["Domains"]["Materials"]["LondonDepth"]`](domains.md##domains%5B%22Materials%22%5D%5B%22LondonDepth%22%5D)
is specified.

`"EstimatorTol" [1.0e-6]` :  Relative tolerance for flux projection used in the
error estimate calculation.

`"EstimatorMaxIts" [10000]` :  Maximum number of iterations for flux projection use in the
error estimate calculation.

`"EstimatorMG" [false]` :  Set to true in order to enable multigrid preconditioner with AMG
coarse solve for the error estimate linear solver, instead of just Jacobi.

`"GSOrthogonalization" ["MGS"]` :  Gram-Schmidt variant used to explicitly orthogonalize
vectors in Krylov subspace methods or other parts of the code.

  - `"MGS"` :  Modified Gram-Schmidt
  - `"CGS"` :  Classical Gram-Schmidt
  - `"CGS2"` :  Two-step classical Gram-Schmidt with reorthogonalization

### Advanced linear solver options

  - `"InitialGuess" [true]`
  - `"MGUseMesh" [true]`
  - `"MGAuxiliarySmoother" [true]`
  - `"MGSmoothEigScaleMax" [1.0]`
  - `"MGSmoothEigScaleMin" [0.0]`
  - `"MGSmoothChebyshev4th" [true]`
  - `"ColumnOrdering" ["Default"]` :  `"METIS"`, `"ParMETIS"`,`"Scotch"`, `"PTScotch"`,
    `"PORD"`, `"AMD"`, `"RCM"`, `"Default"`
  - `"STRUMPACKCompressionType" ["None"]` :  `"None"`, `"BLR"`, `"HSS"`, `"HODLR"`, `"ZFP"`,
    `"BLR-HODLR"`, `"ZFP-BLR-HODLR"`
  - `"STRUMPACKCompressionTol" [1.0e-3]`
  - `"STRUMPACKLossyPrecision" [16]`
  - `"STRUMPACKButterflyLevels" [1]`
  - `"SuperLU3DCommunicator" [false]`
  - `"AMSVectorInterpolation" [false]`
  - `"AMSSingularOperator" [false]`
  - `"AMGAggressiveCoarsening" [false]`
