```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Driven Solver: Uniform vs Adaptive

In this tutorial, we will discuss the [driven solver](../config/solver.md#solver%5B%22Driven%22%5D),
which computes the frequency-domain response (steady-state) of a system driven by external
excitations. We will especially focus on the “adaptive” driven solver, which uses reduced-order
modelling (ROM) techniques. The ROM approach can dramatically increase the speed of a
frequency-domain simulation, but requires more careful usage. The uniform and adaptive driven
solvers are also discussed in the example on [cross-talk between coplanar waveguides](cpw.md) and we
will assume familiarity with that example.

In this tutorial, we will go into more of the algorithmic detail, to help users effectively use the
driven solvers. Additionally, the adaptive solver is heavily used in *Palace*'s [circuit extraction
feature](transmon_circuit.md). We will consider the [CPW example](cpw.md) in more detail as well as
apply the driven solvers to the [single transmon model](transmon.md). The adaptive solver works
particularly well for the transmon case.

!!! warning

    This tutorial is advanced. The algorithmic choices described here are considered “internal” to
    *Palace* and may change at any time if the developers decide that alternative approaches are better.
    Please proceed with caution.

!!! note

    The simulations below were performed version ABC ... on date XYZ.

    TODO: Where to get files and folders.

## Driven Solver Quick-Start

<!-- 
To use the adpative

The following configurations allow you to tune the adaptive driven solver.

  - **`"AdaptiveTol"`**: target relative error tolerance $\varepsilon$. Setting this to a
    positive value activates the adaptive sweep; zero (the default) gives the plain uniform
    solver. This tolerance should be set *larger* than `"Solver.Linear.Tol"`, since the
    ROM introduces an additional approximation on top of the finite-element discretisation.
    Values of `1e-2` to `1e-3` are a good starting point.

  - **`"AdaptiveMaxSamples"`**: hard ceiling on HDM solves per excitation (default: 20).
    If the algorithm reaches this limit before the tolerance is met, it exits using whatever
    ROM has been built and prints a warning. Increase this limit if the default is hit.

  - **`"AdaptiveConvergenceMemory"`**: number of consecutive greedy iterations that must
    all individually satisfy the tolerance before convergence is declared (default: 2).
    This guards against premature termination when the error estimate fluctuates near
    $\varepsilon$. Increasing to 3 or 4 makes the criterion stricter.

  - **`"AddToPROM"`**: number of consecutive greedy iterations that must
    all individually satisfy the tolerance before convergence is declared (default: 2).
    This guards against premature termination when the error estimate fluctuates near
    $\varepsilon$. Increasing to 3 or 4 makes the criterion stricter. -->

## Revisiting the Two Co-planar Waveguide Example (Part I)

We start by re-visiting the driven response of the [two CPW example](cpw.md), using the mesh
[`mesh/cpw_lumped_0.msh`](https://github.com/awslabs/palace/blob/main/examples/cpw/mesh/cpw_lumped_0.msh)
that places lumped ports at the end of the CPWs.

### Uniform Solver

First, we set up a uniform driven simulation using the configuration file
[`cpw_lumped_uniform_convergence.json`](https://github.com/awslabs/palace/blob/main/examples/cpw/cpw_lumped_uniform_convergence.json).
In this example, we will only be driving the single lumped port with `"Index": 1`. The `Solver`
section of the configuration is

```json
"Solver" : {
  "Order" : 2,
  "Device": "CPU",
  "Driven": { "Samples": [ {"Type": "Linear", "MinFreq": 2.0, "MaxFreq": 32.0, "FreqStep": 0.1} ] },
  "Linear": {"Type": "Default", "Tol": 1.0e-12, "MaxIts": 1000, "ComplexCoarseSolve": true}
}
```

*Palace* calls the uniform driven solver when
[`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D) in `Solver/Driven` is `0.0` or not
specified. The uniform solver will iterate over all sample points specified in
`Solver/Driven/Samples` and independently solve the linear equation for the electric field at the
sample frequency point $f$. We will also refer to these solves as "full" or "High-Dimensional
Solves" (HDM) since they solve the full linear problem on the finite element space. Here we are
sampling a linear grid between $2~\mathrm{GHz}$ and $32~\mathrm{GHz}$ in steps of $\Delta f =
0.1~\mathrm{GHz}$, but [users can provide more complicated sample
specification](../config/solver.md#solver%5B%22Driven%22%5D%5B%22Samples%22%5D).

The linear solver tolerance `"Tol": 1.0e-12` is chosen to be very small. Such a small tolerance is
near the limit of what *Palace*'s solvers can reach using double precision. The error of the uniform
solver also depends on the condition number of the system matrix $A(\omega) = \bm{K} + i\omega
\bm{C} - \omega^2 \bm{M} + \bm{A}_{2}(\omega)$ at each frequency $\omega = 2 \pi f$ that it solves.
Importantly, if the system has resonances (poles) close to the real axis, the linear solve may
become very poorly conditioned. This leads to a loss of numerical accuracy.

The total electric energy $E_{\mathrm{elec}}$ is printed out as one of the columns in
`domain-E.csv`.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_energy_uniform.svg" width="80%" />
</p><br/>
```

The scattering matrix $S_{ij}$ is printed out in `port-S.csv`. Since we drive port 1 only in this
simulation, *Palace* will only compute the column $S_{i1}$. We plot the magnitude below.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_sparam_uniform.svg" width="80%" />
</p><br/>
```

This data is equivalent to the data previously shown in the [CPW example](cpw.md), since we consider
the same physical model although we have used slightly different solver parameters as shown above.

### Adaptive Solver

The primary downside of the uniform driven solver is its slow performance. Every output sample has
to be individually calculated, which becomes prohibitive on a fine output grid. This is the problem
the adaptive solver addresses.

To perform a driven simulation with the adaptive driven solver, we set
[`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D) to `1e-1` in the configuration. We
will define this quantity more precisely later. Simulation output, like `domain-E.csv`, has the same
structure as the uniform solver — they are evaluated on the same frequency grid specified by
`"Samples"` in the configuration file. However, the adaptive solver constructs a different
"internal" frequency grid on which it does HDM solves. This internal grid is automatically chosen by
the solver and generally contains far fewer points. From this small internal grid, *Palace*
constructs a reduced order model (ROM) that interpolates the solution at the output grid. We
emphasize that the number of HDM solves done is controlled by the
[`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D), not the output frequency grid.
Adding extra points to the output grid incurs only a small incremental cost compared to performing a
full HDM solve.

Let us plot the point-wise relative error $\vert E_\mathrm{elec,adaptive} -
E_\mathrm{elec,uniform}\vert / \vert E_\mathrm{elec,uniform}\vert$ between the electric
energy $E_{\mathrm{elec}}$ of the uniform driven solver above and the adaptive solver with tol
`1e-1`.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_energy_adaptive_single.svg" width="80%" />
</p><br/>
```

In the plot, the dotted black line at `1e-12` corresponds to the linear residual tolerance `Tol` of
both uniform and adaptive solvers. This is the best-case accuracy floor, below which the difference
is dominated by solver noise rather than systematic adaptive solver tolerance. In fact, the error on
the actual solution can be higher since it depends on the condition number of the linear system
$\kappa(A(\omega))$.

The dashed line at `1e-1` is the adaptive tolerance set by
[`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D). We will define this quantity below;
however, it is important to note that the adaptive tolerance controls the error in the finite
element electric field solution. This does not strictly bound the relative error of a derived
quantity, although it is a decent heuristic guide.

The ROM sample points $f_\mathrm{sample} = \{2.0, 11.08, 30.81, 32.0\}~\mathrm{GHz}$ are marked with
the diamond symbols in the plot. Note that these internal samples are not on the same grid as the
output measurements, although they match at the frequency boundaries $f_{\mathrm{min}}$ and
$f_{\mathrm{max}}$. We see that near the frequencies of the internal samples, the relative error has
a local minimum. Since *Palace* is only performing 4 HDM solves, rather than the 301 of the uniform
solver, the adaptive computation is substantially faster since the interpolation adds only a small
overhead.

Let us now look at a sweep of adaptive driven solvers with different values of
[`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D).

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_energy_adaptive_sweep.svg" width="80%" />
</p><br/>
```

As the tolerance gets smaller, the relative error compared to the uniform solver gets smaller. As
expected, the error is generally above the solver tolerance of `1e-12`. We also see that *Palace*
adds more internal sample frequencies $f_\mathrm{sample}$ (coloured diamonds) as the tolerance
becomes smaller. For this CPW example, the location of the sample frequencies is spread out
relatively evenly in the interval $[f_\mathrm{min}, f_\mathrm{max}]$. Finally, we see that the
solution has systematic minima near the sample points. This will turn out to be true by construction
of the adaptive solver. Local minima in accuracy, not associated with sample points, disappear as
the tolerance decreases and the error landscape changes.

To understand the origin of these results, it is important to understand the internals of the
adaptive solver in more detail.

## Details of the Reduced Order Modeling

### Projective Construction

Let us now consider the ROM construction in *Palace*. As discussed in the [theory
reference](../reference.md), the matrix equation of the driven solver has the form

$$\left[\bm{K} + i\omega \bm{C} - \omega^2 \bm{M} + \bm{A}_{2}(\omega)\right] \bm{x} = i \omega
\bm{b} + \bm{b}_2(\omega).$$

The matrix $\bm{K}$ represents the discretized curl-curl operator, $\bm{M}$ the discretized
displacement term, and $\bm{C}$ the dissipative contributions from ports or finite conductivity. The
vector $\bm{x}$ is the solution of the discretized electric field we want to solve for and $\bm{b}$
is the drive on ports that is exciting the system. The dimension of these vectors $N$ is the
dimension of the finite element space. For the discussion below we will ignore $\bm{A}_{2}(\omega)$
and $\bm{b}_2(\omega)$, which arise from boundary conditions that are non-quadratic in frequency
such as wave ports.

The idea of the adaptive solver is to transform this problem with a large dimension $N$ into a
linear problem with much smaller dimension $n$ that approximates the high-dimensional model well.

```math
\left[\bm{K}_r + i\omega \bm{C}_r - \omega^2 \bm{M}_r \right] \bm{x}_r = i \omega \bm{b}_r
```

This also requires a method to efficiently recover the high-dimensional field solution $\bm{x}$ from
the reduced $\bm{x}_r$. If we have found such a reduced model, we can quickly solve the small linear
system at a specific $\omega$ and recover the full solution and all output measurements. We refer to
the procedure of finding the model as the “offline” or “training” phase, and to using it to
calculate many output samples as the “online” phase.

There is a large literature on constructing dimensional reductions of linear systems and we refer to
the literature for an introduction [1,2]. *Palace* finds a set of $n$ orthogonal vectors $\bm{V}$
that capture the full response of the system in the frequency range of interest $[f_\mathrm{min},
f_\mathrm{max}]$. The projection $\bm{K}_r = \bm{V}^T \bm{K}\bm{V}$, $\bm{b}_r = \bm{V}^T \bm{b}$,
etc., and the prolongation $\bm{x} = \bm{V} \bm{x}_r$ are particularly simple.

The set of basis vectors $\bm{V}$ that *Palace* uses are the orthogonalized components of the HDM
solutions $\bm{x}^*$ at the "internal" sampling frequencies $f_\mathrm{sample} \in [f_\mathrm{min},
f_\mathrm{max}]$. In fact, *Palace* uses $\mathrm{Re}\bm{x}^*$, $\mathrm{Im}\bm{x}^*$ as separate
vectors so that the basis $\bm{V}$ is real. As we add more sampling frequencies, the basis of the
reduced order model grows and increases the accuracy of the projection, until we satisfy a
convergence criterion. This construction is essentially a type of rational interpolation on the
whole domain. It is also related to the mathematics of rational Krylov spaces [2].

The Gram-Schmidt orthogonalization algorithm used to construct $\bm{V}$ from the HDM solutions is
determined by the user option
[`"AdaptiveGSOrthogonalization"`](../config/solver.md#solver%5B%22Driven%22%5D), although a user
should generally not need to change this. By default, it uses classical Gram-Schmidt with
reorthogonalization which is a good high-precision option.

### Choosing Internal Sample Frequencies and Convergence

The above construction does not tell us how to efficiently choose the location of internal sample
frequencies $f_\mathrm{sample}$. For this we use the greedy sampling algorithm developed by
Pradovera in [3]. This is based on a type of rational (barycentric) interpolation for the solution
$\bm{u}$ of linear systems in frequency

$$\bm{u}(\omega) = \frac{\sum_i w_i \bm{u}(\omega_i) / (\omega - \omega_i)}{\sum_i w_i / (\omega -
\omega_i)},$$

where the $\omega_i$ are sample frequencies, $\bm{u}(\omega_i)$ solutions at those frequencies, and
$w_i$ some weight coefficients fitted. To use this interpolation for our driven electromagnetic
solver, we linearize the quadratic $KCM$ operator as a linear operator and identify $\bm{u}^T =
(\bm{x}^T, i \omega \bm{x}^T)$. This type of linearization is familiar to physicists as going from
Newton's equations (second order in time) to Hamilton's equations (first order in time) and doubling
the state dimension from $N$ to $2N$.

We construct a separate barycentric interpolation for each excitation pattern. We initialize the
rational interpolation with HDM solutions at the frequency boundary points $f_\mathrm{min}$ and
$f_\mathrm{max}$ as well as any frequencies explicitly set by the user using `"AddToPROM": true` in
the [sample specification](../config/solver.md#solver%5B%22Driven%22%5D%5B%22Samples%22%5D).
Afterwards, the Pradovera algorithm suggests new sample locations, based on the maximum expected
error of the interpolation. The insight of this approach is that this maximum error can be quickly
calculated from the denominator of the barycentric interpolation itself, without any additional HDM
solves. *Palace* prints the location of the sampling points found during the adaptive solver simulations
to the log-file.

Convergence of the procedure is based on the relative error of the HDM electric field solution
$\bm{x}_\mathrm{HDM}$ compared to the solution calculated by the ROM:

$$\varepsilon = \frac{\vert\vert \bm{x}_\mathrm{HDM}(f^*) -
\bm{x}_\mathrm{ROM}(f^*)\vert\vert}{\vert\vert \bm{x}_\mathrm{HDM}(f^*) \vert\vert}$$

Specifically, we declare convergence when this error is below the user-defined tolerance
$\varepsilon < \varepsilon_\mathrm{tol}$, set by
[`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D), for a number of consecutive
frequency samples, set by
[`"AdaptiveConvergenceMemory"`](../config/solver.md#solver%5B%22Driven%22%5D). Note that the norm
$\vert\vert \bm{x} \vert\vert$ is the L2 norm in the finite element space. So this is similar to,
but not the same as, the physical norm computed from the overlap in space $\int \bm{E}^*(\bm{r})
\cdot \bm{E}(\bm{r}) d^3r$, which corresponds to a mass-weighting in finite-elements
$\bm{x}^T\bm{M}\bm{x}$. This difference between L2 and mass-weighted overlap can be significant on
highly anisotropic meshes.

As mentioned above, the relative error is of the electric fields over the global domain. However,
the relative error on a derived measurement need not satisfy that error bound. This occurs when the
measurement is difficult, such as integrating the electric field in a region where the field is
tiny.

### Adaptive Solver Problems and User Guidance

The rational interpolation approach above is extremely efficient at building a basis using only a
very small number of HDM solves. However, there are two challenges that we encounter in practical
use.

**Early termination** due to accidental convergence. As discussed above, the convergence criterion
is that the error from the HDM solve at the next suggested sample point $f^*$ is below the threshold
$\varepsilon_\mathrm{tol}$. Although the algorithm picks $f^*$ to be the position where the error
indicator is largest, the error evaluated there can sometimes be smaller than the error expected
from the interpolation. This means that the solution looks more converged than it is and the
algorithm terminates early. This is discussed in the paper by Pradovera [3].

To avoid this, we introduce a memory. The algorithm only terminates once a certain number of
consecutive samples in a row are below tolerance. This number is defined by the
[`"AdaptiveConvergenceMemory"`](../config/solver.md#solver%5B%22Driven%22%5D) configuration. The
default is `2`, which is small. Increasing this number to say `3` or `4`, especially in models where
the user does not understand the convergence properties, can be helpful. However, each HDM sample
adds a substantial computational cost to the training phase.

**Numerical problems** due to finite precision. The mathematical formulation of the rational
interpolation for the error indicator assumes that the solution at the sample points is exact.
However, in reality there is the linear solver tolerance
[`"Linear"/"Tol"`](../config/solver.md#solver%5B%22Linear%22%5D) that sets the error on this and can
propagate that error to the prediction of the samples. A detailed discussion of this error
propagation is beyond the scope of this tutorial. In practice this is not a problem, provided that
your value of [`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D) is much larger than
the linear solver [`"Tol"`](../config/solver.md#solver%5B%22Linear%22%5D). What “much larger” means
is a model-dependent quantity, although we have found that a factor $10^4$ or larger works well.

What happens if the [`"AdaptiveTol"`](../config/solver.md#solver%5B%22Driven%22%5D) is too close to
the linear solver [`"Tol"`](../config/solver.md#solver%5B%22Linear%22%5D)? In practice, we have
observed that the rational interpolation then cannot be accurately fit, which *Palace* will log as
warnings. This breaks the next best frequency prediction, which causes the adaptive solver to
frequently oversample certain frequency regions and sometimes undersample other regions. While this
does not break the ROM structure as such — it will accommodate a larger number of basis states — the
quality of that basis can be poor and the desired tolerance not reached.

The number of HDM samples for the ROM is capped by the configuration
[`"AdaptiveMaxSamples"`](../config/solver.md#solver%5B%22Driven%22%5D) (default: 20), after which
the solver just uses the capped ROM to do the computation. That cap is helpful to bound the
oversampling cases due to numerical issues. However, under normal conditions the user should
increase [`"AdaptiveMaxSamples"`](../config/solver.md#solver%5B%22Driven%22%5D) to guarantee
accuracy.

As discussed above, a user can also force certain frequency points to be added to the PROM before
the adaptive sampling starts by using [`"AddToPROM": true`](../config/solver.md#solver%5B%22Driven%22%5D%5B%22Samples%22%5D). This is primarily a
debugging tool. Routine usage is not recommended — adding sample by hand negates the major benefit
of the adaptive solver to “choose the best” samples for a given accuracy. However, this option might
help investigate convergence issues. If adding points, the user should be careful not to accidentally
set [`"AddToPROM": true`](../config/solver.md#solver%5B%22Driven%22%5D%5B%22Samples%22%5D) and a dense grid, since this will perform a very large number of HDM
solves and store them in memory.

### Multi-excitations

In the case of [multi-excitations](../guide/boundaries.md#Lumped-and-wave-port-excitation), the
adaptive solver will construct a barycentric rational interpolation of each excitation separately
and perform frequency prediction/sampling separately for each excitation. The boundary frequencies
$f_{\mathrm{min}}$, $f_{\mathrm{max}}$ and the [`"AddToPROM": true`](../config/solver.md#solver%5B%22Driven%22%5D%5B%22Samples%22%5D) user-samples will be added
separately to the interpolation for each excitation.

However, the ROM itself and basis $\bm{V}$ will be the combination of samples from all excitations
and that will enter the stopping criterion $\bm{x}_\mathrm{ROM}(f^*)$. It can sometimes happen that
later excitations already have a good enough basis from the samples of previous excitations and
return earlier than if it was simulated by itself.

Note that `"AdaptiveMaxSamples"` refers to the number of samples *per excitation* in the
multi-excitation case.

### Non-Quadratic Boundary Conditions in the Adaptive Solver

In the above discussion, we have neglected the possibility of an $\bm{A}_2(\omega)$ term. This is a
non-quadratic frequency term, which can arise from `WavePort` / `WavePortPEC`,  `Conductivity`, or
`Farfield` boundary conditions. If such terms are present, the above algorithms will mostly proceed
as described above — there will be frequency samples with HDM solves of the full system now
containing $\bm{A}_2(\omega)$, which will be added to the basis $\bm{V}$. During the online phase,
*Palace* uses the projected $\bm{K}_r$, $\bm{C}_r$, $\bm{M}_r$, $\bm{A}_{2,r}(\omega)$ to solve the
reduced problem efficiently. Having a $\bm{A}_2(\omega)$ term will be slower, since *Palace* will
have to recompute $\bm{A}_{2,r}(\omega)$ for every output $\omega$.

However, the interplay between these non-quadratic boundary conditions and the frequency prediction
algorithm is more delicate. The frequency prediction algorithm does not naturally take these
contributions into account, since they cannot be linearized into a system of size $2N$. The error of
linearizing $\bm{A}_2(\omega)$ should be small if its non-$KCM$ contributions are small as compared
to the $KCM$ contribution of the rest of the system. However, if non-$KCM$ contributions of
$\bm{A}_2(\omega)$ are large, it is possible that the frequency sampling will choose a sub-optimal
basis. The tuning described above in [Adaptive Solver Problems and User
Guidance](#adaptive-solver-problems-and-user-guidance) could be helpful for this case too.

## Revisiting the Two Coplanar Waveguide Example (Part II)

With the knowledge of how the adaptive solver works internally, let us revisit the CPW example
discussed above. We can now interpret several features of the convergence plot better:

 1. The origin of the local minima in the relative error near the sampling frequencies. Because the
    adaptive solver performs a full HDM solution at the sample frequency and uses that as a basis
    vector for the ROM, the response at that frequency would be exactly reproduced. This is seen
    exactly at $f_\mathrm{min}$ and $f_\mathrm{max}$. Output frequency points near the sample
    frequency are still much more accurate than expected.
 2. The fact that the sampling points used are spaced roughly evenly in frequency. Because the
    algorithm for finding the next best frequency is based on a rational interpolation, the next
    error is largest near any unaccounted pole in the response of the system. The CPW model does not
    have high-Q poles near the real frequency domain of interest $[f_\mathrm{min}, f_\mathrm{max}]$,
    so the error estimator is not particularly structured in frequency. We will discuss this idea in
    more detail in the next section on a transmon model.
 3. The difference between the `"AdaptiveTol"` and the relative error of the output quantity. The
    adaptive tolerance is a relative convergence error for the finite element coefficients of the
    electric field. Any derived quantity will have an absolute error bounded proportionally by this
    accuracy and the sensitivity to the electric field. `"AdaptiveTol"` is at best a heuristic guide
    for the relative error of a derived quantity. The true relative error could be far better or far
    worse than `"AdaptiveTol"`.

Point 3 above is really about the nature of solver convergence and errors generally, rather than
being specific to the adaptive solver. However, it is helpful to illustrate this point in more striking detail by looking at the convergence of the scattering parameters with different `"AdaptiveTol"`.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_sparam_adaptive_pointwise.svg" width="80%" />
</p><br/>
```

The plot above shows the relative error $\vert S_{\mathrm{adaptive}} - S_{\mathrm{uniform}}\vert /
\vert S_{\mathrm{uniform}}\vert$ evaluated at each output frequency (point-wise error). We see that
the relative error can be very bad! However, this should not be surprising. The only scattering
matrix that has a large value (of order 1) is $S_{21}$ and we see that the relative error there is
well behaved. The bad relative error in the other plots is because the magnitude of $S_{i1}$ is tiny
to begin with, which magnifies the relative error. For example, the magnitude $\vert S_{11} \vert$
drops down to around $-60~\mathrm{dB} = 10^{-3}$ near $15~\mathrm{GHz}$. That is exactly the
location where the relative error of $S_{11}$ spikes in the error plot above. The relative error in
$S_{31}$, $S_{41}$ is always large since their values are small.

A detailed discussion of error estimates and error propagation from the electric field to derived
quantities is beyond the scope of this tutorial. From a practical perspective, we encourage users to
validate some of the adaptive driven simulations against uniform driven simulations in order to get
a practical sense of the error for different quantities, to help identify an appropriate adaptive
tolerance and to diagnose convergence error.

A useful and cheap visualization may also be to plot the absolute error normalized by an
”appropriate” scale factor. For example, below we plot the $\vert S_{\mathrm{adaptive}} -
S_{\mathrm{uniform}}\vert /  \vert\vert S_{\mathrm{uniform}}\vert\vert_{\mathrm{RMS}}$. Here
$\vert\vert S_{\mathrm{uniform}}\vert\vert_{\mathrm{RMS}}$ is the root-mean-square average of $\vert
S \vert$ over the frequency range. While this does not substitute for a proper error analysis, it
gives a sense of convergence while removing the extreme spikes in the pointwise-relative error.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_domain_sparam_adaptive_rms.svg" width="80%" />
</p><br/>
```

Finally, let us look at how the error indicator $\varepsilon$ calculated by the adaptive solver
evolves as new sample points are added.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/driven_ua_cpw_adaptive_convergence_curve.svg" width="80%" />
</p><br/>
```

The light grey star points $1,2$ correspond to initial samples at $f_\mathrm{min}$ and
$f_\mathrm{max}$. We actually exclude these points from convergence (i.e. they have error "`inf`").
They are drawn here as points at $1.0$ only for illustration. As the adaptive algorithm continues,
more sample points are found that decrease the error indicator. In our simulation, we used the
default `"AdaptiveConvergenceMemory"` value of 2. So the simulations converge after two consecutive
points are below `"AdaptiveTol"`. As we can see from the plot above, the error indicator does not
have to decrease monotonically. This means that the error indicator can go above the target
tolerance before dipping below it again. See the discussion above in [Adaptive Solver Problems and
User Guidance](#adaptive-solver-problems-and-user-guidance).


## Uniform Driven Solver

We use exactly the same mesh and basic set up as in the tutorial. But we modify our configuration files as appropriate for a driven solver.

Specifically 


When we run this configuration, *Palace* will iterate over each requested frequency point and solve the full linear equation $$[]$$, see the reference section. 


Driven Solver Config file: Uniform solve (what does it do)
- Single-Excitation & Multi-excitation

We can look at — for example — the scattering parameters:



We see resonance dips corresponding to the eigenmodes as well as the background signal.

These scattering parameters are currently only calculated on purely dissipative ports with the assumption that the reference impedance $Z_R$ and surface impedance values $R$ specified the configuration file.

Of course, we can analyze this data

Fitting the uniform solve using Vector Fit and AAA with determinant surrugate, plot. 



Can get response and eigenmodes, but requires more detailed checks and tuning.
Validate against eigenmodes. We will return to this point in our tutorial on 
synthisizing circuits (LINK).

## Adaptive Driven Solver and Reduced Order Model

- The uniform solver is great baseline, but quite slow since it solves indepedantly for every frequenccy.
- If we make a finer mesh, this becomes really bad
- As we can see from our vector fitting examples, we can extract much of the information from a much smaller sample of data. This leads us to the idea of reduced-order modeling.
- The idea is to solve the full system only at a few key frequency points and the recustruct the rest of the response at other frequncies bases on that. This is implemented in the "adaptive" driven solver of *Palace*.
- We will first show and example of running with this solver and then return to discuss the
  algorithmic details.


### A first run of the adaptive tol


### Background reduced order model based on solutions

- Adaptive Solve based on Reduced order modeling
	- Give good references at end
	- How does it work? Mode shapes
- Only manages to fit what it can see — think of it as an rational matrix interpolation on the real axis

### More details reduced order modeling solver

- Tolerances and what they mean
- How frequency points are chose
- How does it pick frequencies: manual "AddToPROM" and "AdaptiveTol"
- AdaptiveConvergenceMemory

### Tuning tolerances and convergence

What is the error tolerance.


## Summary and take-aways

We have now reached a natural break if this tutorial 

- Works really well since there are few poles in or close to the region of interest.
- Takeaway: rom_tol >> solver_tol but there are limits.
- It is a greedy algorithm — default AdaptiveConvergenceMemory is low at two — if it is a difficult case might need to increase for challenging problem.
- When not to use this — (small number of samples compared to number of poles needed). High precision in difficult region.

A more difficult example is here.

---

## Obtaining a circuit out of the Adaptive Driven Solver

- How to interpret the circuits
- Is this a normal circuit?
- How to post-process these circuits

## CPW Line with and LC Port

!! 

- Convergence: "AdaptiveTol"


- Getting the circuit parameters
- Changing circuit parameters — warning

- CPW Line with port — a more difficult example
	- Why?


## Literature & References

[1] P. Benner, D. C. Sorensen, and V. Mehrmann, Eds., Dimension Reduction of Large-Scale Systems: Proceedings of a Workshop held in Oberwolfach, Germany, October 19–25, 2003. in Lecture Notes in Computational Science and Engineering, no. 45. Berlin, Heidelberg: Springer Berlin Heidelberg, 2005. doi: 10.1007/3-540-27909-1.

[2] A. C. Antoulas, Approximation of Large-Scale Dynamical Systems. Society for Industrial and Applied Mathematics, 2005. doi: 10.1137/1.9780898718713.

[3] D. Pradovera, “Toward a certified greedy Loewner framework with minimal sampling,” Adv Comput Math, vol. 49, no. 6, p. 92, Dec. 2023, doi: 10.1007/s10444-023-10091-7.

[4] D. Pradovera, “Interpolatory rational model order reduction of parametric problems lacking uniform inf-sup stability,” SIAM J. Numer. Anal., vol. 58, no. 4, pp. 2265–2293, Jan. 2020, doi: 10.1137/19M1269695.
