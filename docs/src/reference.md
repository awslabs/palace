```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Reference

## Mathematical background

The solver computes a finite element approximation to the three-dimensional, time-harmonic
Maxwell's equations in second-order form. The nondimensionalized, source-free, boundary
value problem for ``\bm{E}(\bm{x})\in\mathbb{C}^3``, ``\bm{x}\in\Omega``,
``\partial\Omega = \Gamma``, where
``\bm{\mathscr{E}}(\bm{x},t) = \text{Re}\{\bm{E}(\bm{x})e^{i\omega t}\}`` denotes the
electric field, is written as

```math
\begin{aligned}
\nabla\times\mu_r^{-1}\nabla\times\bm{E} + i\omega\sigma\bm{E}
    - \omega^2\varepsilon_r\bm{E} &= 0 \,,\; \bm{x}\in\Omega \\
\bm{n}\times\bm{E} &= 0 \,,\; \bm{x}\in\Gamma_{PEC} \\
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{E}) &= 0 \,,\; \bm{x}\in\Gamma_{PMC} \\
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{E})
    + \gamma\bm{n}\times(\bm{n}\times\bm{E}) &= \bm{U}^{inc} \,,\; \bm{x}\in\Gamma_{Z}
\end{aligned}
```

where the nondimensionalization has been performed with respect to a characteristic length
``L_0``, time ``L_0/c_0``, magnetic field strength ``H_0``, and electric field strength
``Z_0 H_0``. Here, ``c_0`` and ``Z_0`` are the speed of light and impedance of free space,
respectively. This nondimensionalization will be used throughout this entire reference. For
more details, see [[1]](#References) and [[2]](#References).

Given the electric field solution, the time-harmonic magnetic flux density can be calculated
as

```math
\bm{B} = -\frac{1}{i\omega}\nabla\times\bm{E} \,.
```

The flux density is related to the magnetic field, ``\bm{H}``, by the standard linear
constitutive relationship ``\bm{H} = \mu_r^{-1}\bm{B}``.

In general, the material property coefficients may be scalar- or matrix-valued. In the
matrix-valued (anisotropic) case, the material property coefficients should still always be
symmetric.

For a general isotropic lossy dielectric, the relative permittivity ``\varepsilon_r`` is a
complex-valued quantity:

```math
\varepsilon_r = \varepsilon_r' (1-i\tan{\delta})
```

where ``\varepsilon_r'`` is the real relative permittivity and ``\tan{\delta}`` is the loss
tangent. Alternatively, conductor loss is modeled by Ohm's law ``\bm{J} = \sigma\bm{E}``
with electrical conductivity ``\sigma > 0.0``. For a superconducting domain, the constitive
current-field relationship given by Ohm's law is replaced by that given by the London
equations:

```math
\frac{\partial \bm{J}}{\partial t} = \frac{1}{\mu_r\lambda_L^2}\bm{E}
```

where ``\lambda_L = \sqrt{m/\mu n_s e^2}/L_0`` is the nondimensionalized London penetration
depth. In this case, the term ``+i\omega\sigma \bm{E}`` arising for a normal conductor in
the time-harmonic Maxwell's equations becomes ``+(\mu_r \lambda_L^2)^{-1}\bm{E}``.

The domain boundary ``\Gamma = \Gamma_{PEC}\cup\Gamma_{PMC}\cup\Gamma_{Z}``, is separated
into perfect electric conductor (PEC), perfect magnetic conductor (PMC), and impedance
boundaries, respectively. The PEC boundary condition is a homogeneous Dirichlet condition,
while the PMC boundary condition is the natural boundary condition for the problem and is
satisfied at all exterior boundaries by the finite element formulation. Impedance
boundaries are modeled using a Robin boundary condition with ``\gamma = i\omega/Z_s``, in
which ``Z_s`` the surface impedance of the boundary, with units of impedance per square.

## Floquet periodic boundary conditions

When applying Floquet periodic boundary conditions, the phase delay is incorporated into
the time-harmonic Maxwell equations and exact periodic boundary conditions are applied.
The modified Maxwell equations are obtained by substituting
``\bm{E}(\bm{x}) = \bm{E}_p(\bm{x})e^{-i \bm{k}_p \cdot \bm{x}}``, where ``\bm{E}_p`` is
the periodic electric field and  ``\bm{k}_p`` is the user-specified Bloch wavevector.
The resulting equation is

```math
\begin{aligned}
\nabla\times\mu_r^{-1}\nabla\times\bm{E}_p
- i\bm{k}_p\times\mu_r^{-1}\nabla\times\bm{E}_p
- i\nabla\times(\mu_r^{-1}\bm{k}_p\times\bm{E}_p) & \\
- \bm{k}_p\times\mu_r^{-1}\bm{k}_p\times\bm{E}_p
+ i\omega\sigma\bm{E}_p
- \omega^2\varepsilon_r\bm{E}_p &= 0 \,,\; \bm{x}\in\Omega
\end{aligned}
```

and given the electric field solution, the time-harmonic magnetic flux density can be calculated
as

```math
\bm{B}_p = -\frac{1}{i\omega}\nabla\times\bm{E}_p + \frac{1}{\omega} \bm{k}_p \times \bm{E}_p \,.
```

## Time domain formulation

A time-dependent formulation is also available to compute the electric field response
``\bm{E}(\bm{x},t)`` for a given time-dependent source excitation
``\bm{U}^{inc}(\bm{x},t)``. The governing equations in this case are

```math
\nabla\times\mu_r^{-1}\nabla\times\bm{E} + \sigma\frac{\partial\bm{E}}{\partial t}
    + \varepsilon_r\frac{\partial^2\bm{E}}{\partial t^2} = 0 \,,\; \bm{x}\in\Omega
```

subject to the same boundary conditions as the frequency-dependent case except for the Robin
boundary condition which is written for a lumped resistive port boundary, for example, as

```math
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{E})
    + Z_s^{-1}\bm{n}\times\left(\bm{n}\times\frac{\partial\bm{E}}{\partial t}\right)
    = \bm{U}^{inc} \,,\; \bm{x}\in\Gamma_{Z} \,.
```

The second-order electric field differential equation is transformed into a first-order
ODE system which is solved along with the equation for the magnetic flux density

```math
\left(\begin{matrix} \varepsilon_r & 0 & 0 \\ 0 & I & 0 \\ 0 & 0 & I\end{matrix}\right)
  \left(\begin{matrix} \ddot{\bm{E}} \\ \dot{\bm{E}} \\ \dot{\bm{B}}\end{matrix} \right)
  = \left(\begin{matrix} -\sigma & -\nabla\times\mu_r^{-1}\nabla\times & 0\\ I & 0 & 0 \\ 0 & -\nabla\times & 0\end{matrix}\right)
    \left(\begin{matrix}\dot{\bm{E}}\\ \bm{E} \\ \bm{B} \end{matrix}\right) \,.
```

The first-order ODE system formulation is chosen to take advantage of implicit adaptive
time-stepping integration schemes. The ``3 \times 3`` system can be block-eliminated to
avoid an expensive coupled block system solve. It offers the additional benefit
of sharing many similarities in the spatial discretization as the frequency domain
formulation outlined above.

## Eigenmode calculations

For eigenmode problems, the source term is zero and a quadratic eigenvalue problem for the
eigenvalues ``\omega`` is solved:

```math
(\bm{K}+i\omega\bm{C}-\omega^2\bm{M})\bm{x} = 0
```

where the matrix ``\bm{K}`` represents the discretized curl-curl operator, ``\bm{M}`` the
mass term, and ``\bm{C}`` the port impedance boundary conditions. The damped frequency
``\omega_d`` and quality factor ``Q`` are postprocessed from each of the resulting
eigenvalues as

```math
\omega_d = \text{Re}\{\omega\} \,, \qquad Q = \frac{|\omega|}{2|\text{Im}\{\omega\}|} \,.
```

When wave port, surface conductivity, or second-order absorbing boundary conditions are used,
a nonlinear eigenvalue problem is solved:

```math
(\bm{K}+i\omega\bm{C}-\omega^2\bm{M}+\bm{A}_2(\omega))\bm{x} = 0
```

where the matrix ``\bm{A}_2`` represents the nonlinear frequency-dependent boundary conditions.

The eigenmodes are normalized such that they have unit norm and their mean phase is a positive real number.

## Lumped ports and wave ports

For lumped port boundaries, the surface impedance can be related to an equivalent circuit
impedance, ``Z``. There are two common cases:

 1. *Rectangular ports*: ``Z = Z_s l / w``, where ``l`` and ``w`` are the length and width
    of the port, respectively (length here is defined as the distance between the two
    conductors).

 2. *Coaxial ports*: ``Z = Z_s \ln(b/a) / 2\pi``, where ``a`` and ``b`` denote the inner and
    outer radii of the port, respectively.

A lumped parallel RLC circuit boundary has a circuit impedance

```math
\frac{1}{Z} = \frac{1}{R}+\frac{1}{i\omega L}+i\omega C \,.
```

Thus, the relationships between the circuit and surface element parameters for the user to
specify are given by ``R_s = \alpha R``, ``L_s = \alpha L``, and ``C_s = C/\alpha``, where
``\alpha = w/l`` for a rectangular port or ``\alpha = 2\pi / \ln(b/a)`` for a coaxial
port.

For multielement lumped ports, the effective circuit impedance is given by

```math
\frac{1}{Z} = \sum_k \frac{1}{Z_k} \,.
```

That is, the circuit impedances of each port contributing to the multielement port add in
parallel. For the specific case of a two element multielement port with two identical
lumped elements, we have ``Z = (1/Z_1 + 1/Z_2)^{-1} = Z_k / 2``, where ``Z_k`` is the
circuit impedance of a single port element.

The source term ``\bm{U}^{inc}`` in a driven frequency-response problem is related to the
incident field at an excited port boundary by

```math
\bm{U}^{inc} = -2\gamma_{exc}(\bm{n}\times\bm{E}^{inc})\times\bm{n} \,, \qquad
\gamma_{exc} = \frac{i\omega}{R_{ref}}
```

where ``(\bm{n}\times\bm{E}^{inc})\times\bm{n}`` is just the projection of the excitation
field onto the port surface, and ``R_{ref}`` is the real reference resistance of the
excited port: the port resistance ``R`` when ``R > 0``, or the impedance of free space
``Z_0`` for a purely reactive port with ``R = 0`` (see the discussion of reactive port
excitation below). The excitation coefficient ``\gamma_{exc}`` is in general distinct from
the termination coefficient ``\gamma = i\omega/Z_s`` appearing in the Robin boundary
condition, and the two coincide for a purely resistive port. The incident fields for lumped
ports depend on the port shape:

 1. *Rectangular ports*: ``\bm{E}^{inc} = E_0 \, \hat{\bm{l}}``, where ``E_0`` is a uniform
    constant field strength and ``\hat{\bm{l}}`` a unit vector defining the direction of
    polarization on the port (typically should be the direction between the two conductors).

 2. *Coaxial ports*: ``\bm{E}^{inc} = \frac{E_0 r_0}{r} \, \hat{\bm{r}}``, where ``E_0`` is
    again a uniform constant field strength, ``r_0`` is a characteristic length for the
    port, ``r`` is the distance from the port center, and ``\hat{\bm{r}}`` a unit vector
    specifying the port radial direction.

In the time domain formulation, the source term ``\bm{U}^{inc}`` appears as

```math
\bm{U}^{inc} = -2 R_{ref}^{-1}\left(\bm{n}\times\frac{\partial\bm{E}^{inc}}{\partial t}\right)
    \times\bm{n} \,.
```

The incident field ``\bm{E}^{inc}(\bm{x},t)`` is

```math
\bm{E}^{inc}(\bm{x},t) = p(t)\bm{E}^{inc}(\bm{x})
```

where ``\bm{E}^{inc}(\bm{x})`` is identical to the spatial excitation in the frequency
domain formulation, and ``p(t)`` describes the temporal shape of the excitation. Possible
options include a sinusoidal, Gaussian, modulated Gaussian, or step excitation.

In the frequency domain, the scattering parameters can be postprocessed from the computed
electric field for each lumped port with boundary ``\Gamma_i`` as

```math
S_{ij} = \frac{\displaystyle\int_{\Gamma_i}\bm{E}\cdot\bm{E}^{inc}_i\,dS}
    {\displaystyle\int_{\Gamma_i}\bm{E}^{inc}_i\cdot\bm{E}^{inc}_i\,dS} - \delta_{ij} \,.
```

For a resistive excited port, the incident field ``\bm{E}^{inc}`` is normalized so that the
power integrated over the port boundary is unity, referenced to the port resistance
``R = R_{ref}``, and the excitation and termination coefficients coincide. A lumped port
with nonzero inductance and/or capacitance may also be excited. In that case the port's
full ``R``, ``L``, and ``C`` still enter the system matrix as a physical parallel
termination through ``\gamma = i\omega/Z_s`` with the complex ``Z_s(\omega)``, while the
source term is assembled with the real ``\gamma_{exc} = i\omega/R_{ref}`` defined above
(``R_{ref} = Z_0 \approx 376.7~\Omega`` for a purely reactive port with ``R = 0``); the
port reactance therefore acts on the solution through the termination, not through the
drive normalization. The
scattering parameters at resistive ports (observation or drive) retain their standard power-wave
interpretation. For a purely reactive drive port (``R = 0``), the reported ``S_{ii}`` is related
to the loaded port impedance by

```math
S_{ii} + 1 = \frac{2\,Z_{\text{loaded}}}{Z_0}
```

where ``Z_{\text{loaded}} = Z_p \parallel Z_{\text{struct}}`` is the parallel combination of
the port impedance ``Z_p = (1/R + 1/(i\omega L) + i\omega C)^{-1}`` (zero-valued branches
are omitted from the sum: the ``1/R`` term is dropped when ``R = 0``, and likewise for
``L`` and ``C``) and the structure's input impedance at the port plane. This is a
reciprocal, frequency-dependent quantity — but it is not bounded by unity and should not
be interpreted as a power reflection coefficient. Transmission parameters ``S_{ji}``
(``j \ne i``) between a reactive drive and a resistive observation port are
real-reference-normalized transmission amplitudes: each side is normalized to its own real
reference (``Z_0`` at the reactive drive, ``R_j`` at the resistive observation port), so a
fixed ``\sqrt{R_j/Z_0}`` conversion factor separates them from a voltage-transfer ratio.
They satisfy reciprocity (``S_{ji} = S_{ij}``). The same unit-reference projection applies
to *passive* purely reactive lumped ports (``R = 0`` with ``L`` and/or ``C``) appearing as
observation ports: their S-parameter rows are finite and referenced to ``Z_0``, rather
than identically zero. Likewise, the incident voltage and current columns of the port
postprocessing output (``V^{inc}``, ``I^{inc}``) at a purely reactive excited port are
finite and referenced to ``Z_0`` — like the S-parameter row, they are bookkeeping
quantities consistent with the assembled drive normalization, not physical incident
waves (an ``R = 0`` port supports no incident traveling wave).

In the time domain, the time histories of the port voltages can be Fourier-transformed to
get their frequency domain representation for scattering parameter calculation.

Numeric wave ports assume a field with known normal-direction dependence
``\bm{E}(\bm{x}) = \bm{e}(\bm{x}_t)e^{ik_n x_n}`` where ``k_n`` is the propagation constant.
For each operating frequency ``\omega``, a two-dimensional eigenvalue problem is solved on
the port yielding the mode shapes ``\bm{e}_m`` and associated propagation constants
``k_{n,m}``. These are used in the full 3D model where the Robin port boundary condition has
coefficient ``\gamma = i\text{Re}\{k_{n,m}\}/\mu_r`` and the computed mode is used to
compute the incident field in the source term ``\bm{U}^{inc}`` at excited ports. Scattering
parameter postprocessing takes the same form as the lumped port counterpart using the
computed modal solutions. Since the propagation constants are known for each wave port,
scattering parameter de-embedding can be performed by specifying an offset distance ``d``
for each port:

```math
\tilde{S}_{ij} = S_{ij}e^{ik_{n,i}d_i}e^{ik_{n,j}d_j} \,.
```

For more information on the implementation of numeric wave ports, see [[3]](#References).

## Adaptive driven solver and reduced-order modeling

The driven solver assembles a frequency-dependent linear system of the form

```math
\bm{A}(\omega)\bm{x} =
\left[\bm{K} + i\omega\bm{C} - \omega^2\bm{M} + \bm{A}_2(\omega)\right]\bm{x}
= i\omega\bm{b} + \bm{b}_2(\omega),
```

where ``\bm{x}`` contains the finite-element electric-field degrees of freedom.
The matrices ``\bm{K}``, ``\bm{M}``, and ``\bm{C}`` represent the discretized curl-curl,
displacement, and dissipative terms, respectively. The additional terms
``\bm{A}_2(\omega)`` and ``\bm{b}_2(\omega)`` account for boundary conditions that are not
quadratic in frequency, such as numeric wave ports or far-field boundaries.

The adaptive driven solver constructs a projection-based reduced-order model (PROM) for this
system. It seeks a real basis ``\bm{Q}\in\mathbb{R}^{N\times n}``, with ``n \ll N``, and
projects the high-dimensional system onto that basis:

```math
\left[\bm{K}_r + i\omega\bm{C}_r - \omega^2\bm{M}_r\right]\bm{x}_r
= i\omega\bm{b}_r,
```

with, for example, ``\bm{K}_r = \bm{Q}^T\bm{K}\bm{Q}`` and
``\bm{b}_r = \bm{Q}^T\bm{b}``. Once the reduced system is solved, the high-dimensional field
is recovered as ``\bm{x}\approx\bm{Q}\bm{x}_r``. The expensive part is the offline
construction of ``\bm{Q}``; the online evaluation of the reduced system is cheap for many
output frequencies. This follows the standard projection-based model-order-reduction
framework [[7]](#References), [[8]](#References).

The PROM basis is built from full high-dimensional model (HDM) solutions at selected
internal sample frequencies. For each complex HDM solution ``\bm{x}^*``, *Palace* adds the
orthogonalized real and imaginary components, ``\mathrm{Re}(\bm{x}^*)`` and
``\mathrm{Im}(\bm{x}^*)``, as separate real basis vectors. The Gram-Schmidt
orthogonalization variant is controlled by
[`config["Solver"]["Driven"]["AdaptiveGSOrthogonalization"]`](config/reference.md#config-solver-driven),
though most users should not need to change it.

### Internal sample selection

The remaining question is how to choose the internal sample frequencies. *Palace* uses a
greedy sampling strategy based on minimal rational interpolation, following Pradovera
[[9]](#References), [[10]](#References). For a linearized state ``\bm{u}(\omega)``, the
interpolation has the barycentric form

```math
\bm{u}(\omega) =
\frac{\sum_i w_i \bm{u}(\omega_i) / (\omega - \omega_i)}
     {\sum_i w_i / (\omega - \omega_i)},
```

where ``\omega_i`` are already-sampled frequencies and ``w_i`` are fitted weights. For the
quadratic ``KCM`` part of the driven system, *Palace* uses the linearized state
``\bm{u}^T = (\bm{x}^T, i\omega\bm{x}^T)``. This is analogous to rewriting a second-order
time system as a first-order system with twice as many state variables.

Each excitation pattern gets its own rational interpolation for selecting future samples.
The initial samples are the frequency-domain endpoints and any user-requested samples with
[`"AddToPROM": true`](config/reference.md#config-solver-driven-samples). After
that, the interpolation suggests the next sample from an error indicator that can be
evaluated without another HDM solve. The selected HDM solution is then added to the shared
PROM basis.

### Convergence criterion

At a proposed sample frequency ``f^*``, *Palace* compares the HDM solution against the PROM
solution using

```math
\varepsilon =
\frac{\|\bm{x}_\mathrm{HDM}(f^*) - \bm{x}_\mathrm{ROM}(f^*)\|}
     {\|\bm{x}_\mathrm{HDM}(f^*)\|}.
```

Convergence is declared when this quantity is below
[`config["Solver"]["Driven"]["AdaptiveTol"]`](config/reference.md#config-solver-driven) for
[`config["Solver"]["Driven"]["AdaptiveConvergenceMemory"]`](config/reference.md#config-solver-driven)
consecutive samples. The norm here is the finite-element coefficient-space L2 norm. It is
related to, but not the same as, a physical energy norm such as
``\int \bm{E}^*(\bm{r})\cdot\bm{E}(\bm{r})\,d^3r``. Consequently, the adaptive tolerance is
not a strict relative-error bound for every derived quantity, such as S-parameters or domain
energies.

The indicator is also not a mathematical certificate of the maximum error over the entire
frequency interval. The convergence memory mitigates accidental early termination when one
suggested sample happens to fall below tolerance before the broader approximation has fully
settled.

### Multi-excitation simulations

For multi-excitation driven simulations, *Palace* constructs a separate rational
interpolation for each excitation, but the PROM basis ``\bm{Q}`` is shared across all
excitations. Samples selected for earlier excitations can therefore improve the basis used
for later excitations. `"AdaptiveMaxSamples"` is interpreted per excitation.

### Non-quadratic boundary terms

When non-quadratic frequency terms ``\bm{A}_2(\omega)`` are present, the HDM samples still
include the full operator and the projected reduced problem uses the projected
``\bm{A}_{2,r}(\omega)`` during online evaluation. This is more expensive than the purely
quadratic case because the projected non-quadratic contribution has to be updated at output
frequencies.

The sample-selection interpolation is based on the quadratic linearization and does not
fully represent arbitrary non-quadratic frequency dependence. If the non-quadratic
contribution is large compared with the ``KCM`` part of the operator, the selected samples
can be less effective and the adaptive solve may require tighter tolerances, more samples,
or additional validation against a uniform sweep.

### Finite-precision effects

The rational interpolation assumes that the sampled HDM solutions are accurate. In practice,
they inherit the error of the linear solver. If
[`config["Solver"]["Linear"]["Tol"]`](config/reference.md#config-solver-linear) is too close
to `"AdaptiveTol"`, the interpolation can become ill-conditioned and *Palace* may warn
about rank-deficient minimal rational interpolation matrices. In that case, tighten the
linear solver tolerance, loosen the adaptive tolerance, or validate with a uniform solve.

## Other boundary conditions

The first-order absorbing boundary condition, also referred to as a scattering boundary
condition, is a special case of the general impedance boundary condition described above:

```math
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{E})
    + i\omega\sqrt{\mu_r^{-1}\varepsilon_r}\bm{n}\times(\bm{n}\times\bm{E}) = 0 \,.
```

This is also known as the Sommerfeld radiation condition, and one can recognize the
dependence on the impedance of free space ``Z_0^{-1} = \sqrt{\mu_r^{-1}\varepsilon_r}``. The
second-order absorbing boundary condition is

```math
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{E})
    + i\omega\sqrt{\mu_r^{-1}\varepsilon_r}\bm{n}\times(\bm{n}\times\bm{E})
    - \beta\nabla\times[(\nabla\times\bm{E})_n\bm{n}] = 0
```

where assuming an infinite radius of curvature ``\beta = \mu_r^{-1}c_0/(2i\omega)``, and the
contribution depending on ``(\nabla\cdot\bm{E}_t)`` has been neglected.

Additionally, while metals with finite conductivity can be modeled using an impedance
boundary condition with constant impedance ``Z_s``, a more accurate model taking into
account the frequency dependence of the skin depth is

```math
Z_s = \frac{1+i}{\delta\sigma}
```

where ``\delta = \sqrt{2/\mu_r\sigma\omega}`` is the skin depth and ``\sigma`` is the
conductivity of the metal. Another model, which takes into account finite thickness effects,
is given by

```math
Z_s = \frac{1}{\delta\sigma}\left(\frac{\sinh{\nu}+\sin{\nu}}{\cosh{\nu}+\cos{\nu}}
    + i\frac{\sinh{\nu}-\sin{\nu}}{\cosh{\nu}+\cos{\nu}}\right)
```

where ``\nu = h/\delta`` and ``h`` is the layer thickness. This model correctly produces the
DC limit when ``h\ll\delta``.

## Energy-participation ratios

The energy-participation ratio (EPR) for lumped inductive elements is computed from the
electric and magnetic fields corresponding to eigenmode ``m``, ``\bm{E}_m`` and
``\bm{H}_m``, using the formula

```math
p_{mj} = \frac{1}{\mathcal{E}^{elec}_m} \, \frac{1}{2} \, L_j I_{mj}^2
```

where ``p_{mj}\in[-1,1]`` denotes the signed participation ratio for junction ``j`` in mode
``m``, ``L_j`` is the provided junction circuit inductance, ``I_ {mj}`` is the peak
junction current for mode ``m``, and ``\mathcal{E}^{elec}_m`` is the electric energy in
mode ``m``. The junction current is computed using the mean voltage across the port,
``\overline{V}_{mj}``, as ``I_{mj} = \overline{V}_{mj}/Z_{mj}``, where
``Z_{mj} = 1/(i\omega_m L_j)`` is the impedance of the inductive branch of the lumped
circuit. The mean port voltage depends on the computed electric field mode and the shape of
the port:

 1. *Rectangular ports*:
    ``\overline{V}_{mj} = \frac{1}{w_j}\int_{\Gamma_j}\bm{E}_m\cdot\hat{\bm{l}}_j\,dS``.

 2. *Coaxial ports*:
    ``\overline{V}_{mj} = \frac{1}{2\pi}\int_{\Gamma_j}\frac{\bm{E}_m}{r}\cdot\hat{\bm{r}}_j\,dS``.

Finally, the total electric energy in mode ``m`` is

```math
\mathcal{E}^{elec}_m
    = \frac{1}{2} \, \text{Re}\left\{\int_\Omega\bm{D}_m^*\cdot\bm{E}_m\,dV\right\}
    + \sum_j \frac{1}{2} \, C_jV_{mj}^2
```

where ``\bm{D}_m = \varepsilon_r\bm{E}_m`` is the electric flux density for mode ``m`` and
the second term on the right-hand side accounts for any lumped capacitive boundaries with
nonzero circuit capacitance ``C_j``.

The EPR can also be used to estimate mode quality factors due to input-output (I-O) line
coupling. The mode coupling quality factor due to the ``j``-th I-O port is given by

```math
Q_{mj} = \frac{\omega_m}{\kappa_{mj}}
```

where the port coupling rate ``\kappa_{mj}`` is calculated as

```math
\kappa_{mj} = \frac{1}{\mathcal{E}^{elec}_m} \, \frac{1}{2}\,R_j I_{mj}^2 \,.
```

## Bulk and interface dielectric loss

The quality factor due to bulk dielectric loss resulting from an electric field ``\bm{E}``
present in domain ``j`` with associated loss tangent ``\tan{\delta}_j`` is given by

```math
\frac{1}{Q_j} = p_j \tan{\delta}_j =
    \frac{1}{\mathcal{E}^{elec}} \, \frac{1}{2} \, \tan{\delta}_j \,
    \text{Re}\left\{\int_{\Omega_j}\bm{D}^*\cdot\bm{E}\,dV\right\}
```

where, as above, ``\mathcal{E}^{elec}`` is the total electric field energy in the domain,
including the contributions due to capacitive lumped elements.

Likewise, the quality factor due to surface interface dielectric loss for interface ``j`` is
given by

```math
\frac{1}{Q_j} = p_j \tan{\delta}_j =
    \frac{1}{\mathcal{E}^{elec}} \, \frac{1}{2} \, t_j\tan{\delta}_j \,
    \text{Re}\left\{\int_{\Gamma_j}\bm{D}^*\cdot\bm{E}\,dS\right\}
```

where ``t_j`` is the thickness of the layer and ``\bm{D} = \varepsilon_{r,j}\bm{E}`` is the
electric displacement field in the layer evaluated using the relative permittivity of the
interface ``\varepsilon_{r,j}``. For an internal boundary, this integral is evaluated on a
single side to resolve ambiguity due to the discontinuity of ``\bm{E}`` across the boundary.

The above formula for interface dielectric loss can be specialized for the case of a
metal-air, metal-substrate, or substrate-air interface [[4]](#References). In each case, the
quality factor for interface ``j`` is given by

  - *Metal-air*:

```math
\frac{1}{Q^{MA}_j} =
    \frac{1}{\mathcal{E}^{elec}} \, \frac{1}{2} \,
    \frac{t_j\tan{\delta}_j}{\varepsilon_{r,j}^{MA}} \,
    \text{Re}\left\{\int_{\Gamma_j}\bm{E}_n^*\cdot\bm{E}_n\,dS\right\}
```

  - *Metal-substrate*:

```math
\frac{1}{Q^{MS}_j} =
    \frac{1}{\mathcal{E}^{elec}} \, \frac{1}{2} \,
    \frac{t_j\tan{\delta}_j(\varepsilon_{r,j}^{S})^2}{\varepsilon_{r,j}^{MS}} \,
    \text{Re}\left\{\int_{\Gamma_j}\bm{E}_n^*\cdot\bm{E}_n\,dS\right\}
```

  - *Substrate-air*:

```math
\frac{1}{Q^{SA}_j} =
    \frac{1}{\mathcal{E}^{elec}} \, \frac{1}{2} \,
    t_j\tan{\delta}_j\left(\varepsilon_{r,j}^{SA} \,
    \text{Re}\left\{\int_{\Gamma_j}\bm{E}_t^*\cdot\bm{E}_t\,dS\right\}
    + \frac{1}{\varepsilon_{r,j}^{SA}} \,
    \text{Re}\left\{\int_{\Gamma_j}\bm{E}_n^*\cdot\bm{E}_n\,dS\right\}\right)
```

where ``\bm{E}_n`` denotes the normal field to the interface and
``\bm{E}_t = \bm{E}-\bm{E}_n`` denotes the tangential field.

For an infinitely thin metal model, the fields in these integrals are singular at the
metal perimeter. Interface dielectric entries may specify one or more `EdgeDistances`
together with exactly one edge source. `EdgeAttributes` manually identifies boundary
surfaces whose geometric perimeter defines the edges. Alternatively, `AutomaticEdges`
extracts the physical perimeter associated with that interface from the configured metal
boundary conditions. Palace then evaluates, for each matching radius ``R``,

```math
\mathcal{E}^{surf}_j(d \geq R), \qquad
\mathcal{E}^{surf}_j(R \leq d < 2R),
```

and writes both energies and their participation ratios to `surface-Q-edge.csv`. These
quantities expose a matching region that is separated from the unresolved singular edge.
They are diagnostics for constructing an edge subgrid correction; they do not by
themselves replace the raw participation in `surface-Q.csv`. The matching radii should be
large compared with the local surface mesh size and small compared with nearby geometric
features.

`AutomaticEdges` is currently supported on 3D meshes for interfaces with explicit `SA`,
`MS`, or `MA` type. It recognizes PEC, terminal or prescribed-potential, conductivity,
impedance, and rational-impedance metal surfaces. Palace forms their global perimeter
graph, removes segments supported by exterior simulation cut surfaces, and selects the
remaining physical segments that geometrically coincide with the requested interface.
This does not require additional boundary splitting beyond the existing interface
attributes. Manual `EdgeAttributes` and automatic extraction are mutually exclusive.

An optional `EdgeDistanceSmoothing` value ``f`` replaces each hard indicator at ``R`` by
a cubic transition from zero to one over
``R(1-f) \leq d \leq R(1+f)``. The reported outside energy uses this smooth indicator,
the inside energy uses its complement, and the annular energy uses the difference
between the indicators at ``R`` and ``2R``. This preserves the inside/outside partition
while reducing changes caused by quadrature points crossing a hard cutoff during mesh
refinement. The default ``f=0`` retains the hard partition.

Setting `LocalizeEdgeEnergy` writes `surface-Q-edge-local.csv`. For every physical
perimeter segment and matching radius, this table reports the segment endpoints and
length, the total surface energy assigned to its nearest-segment region, surface energy
inside ``R``, surface energy in the matching annulus, and electric energy in the
surrounding volume annulus. The latter scales approximately as ``R`` near an ideal
thin-sheet edge and provides an interface-independent local singular-amplitude proxy.
Localized volume integration adds postprocessing cost and can produce large tables on
detailed 3D perimeters.

For automatically extracted edges, the table also reports `p_vertex_in` and
`p_bulk_vertex_ann`. The first is the surface-core participation at points with
``d_{edge}<R`` whose graph distance along the physical edge chain to a non-regular
vertex is less than ``R``. The second applies the same along-edge condition to the
volume annulus ``R \leq d_{edge}<2R``. Non-regular vertices are corners, open endpoints,
and junctions. The graph distance continues across regular mesh vertices, so the
neighborhood does not depend on how a physical chain is subdivided into mesh segments.
These quantities measure how much of the correction input lies near topology which is
not represented by a locally straight, translation-invariant edge coupon. Automatic
closed loops with no non-regular vertices and manual edge selections report zero for
both quantities. When `EdgeDistanceSmoothing` is nonzero, its cubic transition is
applied to both the radial and along-edge cutoffs.

The local table also records the edge geometry used by the automatic path. `component`
labels a connected physical perimeter and `chain` labels a maximal path which continues
through locally straight vertices but stops at corners, endpoints, and junctions.
`v0_type` and `v1_type` classify the segment endpoints as regular (0), corner (1),
endpoint (2), or junction (3). `process_nx`, `process_ny`, and `process_nz` contain the
unit process normal used for polarization, and `automatic` is one when this metadata was
inferred. Manual edge selections use `-1` for unavailable topology labels but still
report their configured, edge-projected process normal. These columns allow correction
postprocessing to pool short mesh segments on one physical edge chain and to treat
corner or junction neighborhoods separately.

Corner classification is necessarily a mesh-geometry heuristic because the finite-element
mesh does not retain whether a polygonal turn came from a sharp CAD vertex or a tessellated
smooth curve. Palace treats local turns below 30 degrees as smooth curve faceting and larger
turns as corners. This removes dependence on ordinary curve tessellation and preserves common
45- and 90-degree layout corners, but a very coarsely faceted curve or an intentional shallow
corner may require inspection of the local table.

For manual edge selection, `LocalizeEdgeEnergy` also requires `EdgeFrameNormal`, a vector
which points from the substrate or fabrication-process side of the metal toward air or
vacuum. With `AutomaticEdges`, Palace instead infers a separate process normal for every
segment from its supporting metal faces. Adjacent material wave speeds orient the normal
from the process or substrate side toward the air or vacuum side, so opposing flip-chip
layers acquire opposite directions without additional boundary splitting.
`EdgeFrameNormal` is optional in this mode and acts as a fallback where the adjacent
materials do not distinguish the two sides.

Palace projects the selected or inferred normal perpendicular to each physical edge
segment and constructs a local right-handed frame: `n` is the projected process normal,
`l` follows the edge, and `m = n x l` is transverse to the edge in the metal plane. In
2D, `l` is the out-of-plane direction. A configured fallback must therefore not be
parallel to any selected 3D edge segment.

The local table resolves surface energy into normal and tangential terms with columns
such as `p_total_n`, `p_in_n`, and `p_ann_n` (and corresponding `_t` columns). It also
resolves volume-annulus energy into
`p_bulk_top_n_ann`, `p_bulk_top_m_ann`, `p_bulk_top_l_ann` and matching `bottom`
columns. `top` is the side toward the configured or inferred process normal; reversing a
configured normal swaps the top and bottom channels. The polarized terms sum exactly to their
unpolarized surface or volume quantity. These channels allow a process calibration to
distinguish fields on the two sides of a zero-thickness sheet and to separate
cross-sectional from edge-longitudinal polarization.

Setting `EdgeRefinement` resolves a correction tube independently of the field-based
error indicator. Before the first field solve, Palace refines every element whose
bounding sphere intersects ``d \leq \alpha R`` until its diameter is at most ``R/N``.
Here ``R`` is `Radius`, ``N`` is `ElementsPerRadius`, and ``\alpha`` is
`OuterRadiusFactor` (default 2). This refinement also occurs when the configured number
of AMR iterations is zero. Identical requests shared by several interface entries are
deduplicated. `Radius` must match one of the entry's `EdgeDistances`, ensuring the
geometrically resolved and suppressed core is the same core used by postprocessing.

After the tube is resolved, `CoreIndicatorWeight` multiplies the field-based AMR
indicator for elements wholly inside ``d<R``. The default value zero prevents the
thin-sheet singular core, which is replaced by the subgrid correction, from dominating
subsequent marking and the AMR stopping criterion. Elements crossing ``R`` retain full
weight so the inside/outside partition remains resolved. This weighting changes mesh
selection only; the field solve and raw participation still use the complete domain.

For MA, MS, and SA interfaces, setting `FluxRecovery` evaluates the normal
electric displacement using the same H(div)-conforming flux recovery used by the
finite-element error estimator. This supplies a single-valued normal trace and can
improve convergence of the outside and annular interface energies. The tangential
SA contribution remains evaluated from the native H(curl) electric field.
`FluxRecovery` is not supported on cracked internal boundaries, including the
zero-thickness PEC sheets used for thin-metal simulations. The recovery is a volume
L2 projection and its boundary trace is not controlled for the singular thin-sheet
field.

### Local fabrication response matrices

An electrostatic fabrication-resolved coupon can be driven by a set of
`PrescribedPotential` boundary traces. A four-column `x,y,z,V` CSV describes an
ordered closed polyline. For a three-dimensional coupon, a five-column
`x,y,z,V,triangle` CSV describes a piecewise-linear surface: each positive integer
triangle index must occur in exactly three rows, and the potential is interpolated
barycentrically on the nearest triangle. Enabling
`Solver.Electrostatic.ResponseMatrix` treats these traces as basis functions
``\phi_i`` and writes `surface-response-matrix.csv`. For every dielectric interface,
physical edge, and matching radius, the file contains the upper triangle of symmetric
matrices ``Q`` such that the interface energy inside the matching radius is

```math
\mathcal{E}^{surf}_{in}(\bm{a}) =
\sum_i a_i^2 Q_{ii} + 2 \sum_{i<j} a_i a_j Q_{ij}
= \bm{a}^T Q \bm{a}.
```

Here ``\sum_i a_i\phi_i`` is the imposed coupon potential trace. The coefficients
``a_i`` scale the voltages stored in the corresponding trace files, and the reported
matrix entries are dimensional energies in joules. Separate `Q_ij normal` and
`Q_ij tangential` columns provide the same operator for the polarized interface terms.
The `Q_total_ij` columns provide total surface energy assigned to the edge's
nearest-segment patch, without the matching-radius cutoff, and its normal and tangential
components. These full-patch operators can be differenced between matched fabricated and
thin coupons to obtain a process-local surface-energy defect operator.
The basis indices are the `Index` values of the `PrescribedPotential` entries.
For a coupled coupon containing two independent conductors, the final excitation combines
a conforming conductor-state trace on `Attributes` with `TerminalAttributes` held at one
volt. This adds the conductor-voltage basis field while all other essential boundaries
remain at zero:

```json
{
  "Index": 97,
  "Attributes": [1],
  "TerminalAttributes": [5],
  "DataFile": "trace-conductor-b.csv"
}
```

When the matching boundary intersects conductor B, the trace must equal one at that
junction; a zero trace would impose incompatible Dirichlet values and create an
artificial singularity. The thin and fabricated coupons must use the same conductor-cut
trace topology: include both endpoints of every finite-thickness cut, remove all
cut-constrained knots from the free basis, and make the final trace zero on conductor A
and one on conductor B.

Palace also writes `domain-response-matrix.csv`. Its matrix ``Q_domain`` satisfies

```math
\mathcal{E}^{domain}(\bm{a}) = \bm{a}^T Q_{domain} \bm{a}
```

for the complete coupon electric-field energy. Equivalently, ``2 Q_domain`` is the
energy-form Dirichlet-to-Neumann operator on the matching contour. Differencing this
matrix between geometrically matched fabricated and thin coupons measures how the
fabrication details change the local contour response independently of the surrounding
device geometry.

This mode requires at least one dielectric entry with `LocalizeEdgeEnergy`. Palace
performs one electrostatic solve for each of the ``N`` basis traces. It obtains cross
terms by applying the existing surface postprocessor to summed basis fields, so no
pairwise PDE solves are required; the remaining surface-integration work scales as
``N^2``. Only the upper triangle is written because ``Q`` is symmetric.

The basis must be converged in the induced interface-energy norm, not only in an L2 norm
of the potential trace. Small trace errors can produce large relative errors in a weak
SA, MS, or MA channel. Local finite-element or spline bases on the matching contour are
therefore generally safer than assuming that a low-order global Fourier truncation is
adequate. Validate the reconstructed trace against direct coupon solves before reusing a
process matrix.

When any participating interface sets `FluxRecovery`, the recovered normal flux also
enters the response matrix. In that case `Solver.Linear.EstimatorTol` controls a physical
postprocessing solve and must remain tight even if the simulation performs no AMR.

### Self-consistent thin-metal response correction

The experimental `Solver.Electrostatic.ResponseCorrection` option uses matched
fabricated- and thin-coupon domain response matrices to alter the global thin-metal
field solve:

```math
K_{corr} = K_{thin} + P^T (S_{fabricated} - S_{thin}) P,
\qquad S = 2 Q_{domain}.
```

Here ``P`` samples the global potential at the coupon contour knots relative to the
potential at `Reference`. Each `Models` entry defines one reusable coupon response, and
each `Patches` entry places that model in the global mesh. `Origin`, `AxisU`, and
`AxisV` map the coupon coordinates into the global frame. For a one-edge coupon,
`AxisU` points from metal into the gap and `AxisV` points from the process side toward
vacuum. The fabricated and thin response matrices must use the contour ordering in
`BasisPoints`, followed by any explicitly documented conductor-state coefficients.

For repeated use of one fabrication process, `Library` replaces the explicit `Models`
and `Patches` lists:

```json
"ResponseCorrection": {
  "Library": "process-library.json",
  "TargetInterfaces": [1, 2, 3],
  "UnmatchedPolicy": "Warn",
  "SolveTol": 1.0e-6
}
```

`TargetInterfaces` is optional. If omitted, Palace considers every typed MA, MS, and SA
interface with edge-distance postprocessing. In 2D, all targets with the same
`EdgeAttributes` form one physical interface group and must use the same
`EdgeFrameNormal`. The largest `EdgeDistances` value must equal the library matching
radius. Palace intersects the perimeter of each group with the perimeter of the union
of all configured metal boundary conditions. This removes points introduced only by a
change in boundary attribute, such as a ground plane continuing into a bump bond. It
also removes endpoints on exterior simulation-domain truncations.

`SolveTol` controls only the optional self-consistent corrected-field solve. The raw
thin-metal field and AMR estimator retain `Solver.Linear.Tol` and `EstimatorTol`,
respectively. The default `1.0e-6` is normally well below the accuracy of the local
fabrication response model; it can be tightened independently when needed.

The remaining sites are assigned to their electrostatic conductor and clustered when
their separation is less than ``2R``. A one-site cluster uses an `IsolatedEdge` model.
A two-site cluster is classified as `SameConductorGap`, `DifferentConductorGap`, or
`SameConductorStrip` from its conductor ownership and the directions in which the metal
continues. Palace first looks for a model within its declared separation tolerance. If
none matches, it linearly combines the nearest lower- and upper-separation models with
the same topology and compatible conductor-state representation. It does not
extrapolate. A cluster with three or more collinear sites requires an exact
`ParallelEdgeCluster` signature. Noncanonical orientations, opposing process normals,
missing cluster signatures, and unbracketed pair separations are unsupported. With
`UnmatchedPolicy` set to `Warn`, Palace reports the unsupported topology and disables
correction for the complete interface group, avoiding a partial correction. `Error`
instead terminates the run.

In 3D, targets selecting the same automatically extracted physical metal segments form
one interface group. Palace infers the process normal and the in-plane direction from
metal toward gap directly from the supporting metal faces. It partitions straight,
parallel edge segments wherever the active set of edges within ``2R`` changes.
Longitudinally overlapping pairs use the same three coupled topologies as 2D, and active
sets with three or more edges require an exact `ParallelEdgeCluster` signature; the
remaining one-edge intervals use `IsolatedEdge`. The 2D coupon response is integrated
along each interval using a Gauss rule selected from the global FEM order. Separation
interpolation evaluates both pair-coupon traces at every longitudinal quadrature point
and combines their domain and surface responses with convex weights. The normalized
bracket width contributes to `max library distance`. A matching `ConvexCorner`,
`ConcaveCorner`, `Endpoint`, or `Junction` response replaces a physical vertex and trims
a distance ``R`` from each
incident straight arm. A rounded corner first uses an angle- and radius-tolerance match.
Otherwise, Palace linearly combines the nearest lower- and upper-radius models with the
same topology and compatible angle. Both bracket radii must be positive: a sharp-corner
model is not an interpolation endpoint, and Palace does not extrapolate outside the
available rounded-radius range. A rounded corner with neither an exact match nor a
positive-radius bracket remains uncorrected and contributes to the unmatched-corner
diagnostics. The response operators and conductor reference are combined with convex
weights; the normalized bracket width contributes to `max library distance`. Locally
connected arms are not treated as an unsupported nearby-edge pair. An exact
`SpatialEdgeCluster` model can replace a localized nonparallel or endpoint-adjacent
interaction among physical edge chains, including an interaction across target-interface
groups when the model declares explicit interface slots. Nonparallel nearby edges
without that model, cross-layer interactions without an exact multi-slot cluster,
multi-edge active sets without an exact cluster signature, a nearby edge outside the
selected target interfaces, and partially overlapping interface groups fail closed
locally. With `UnmatchedPolicy` set to `Warn`, Palace omits the affected physical edge
segments, reports the reduced per-interface matched fraction, and continues to correct
supported segments in the same interface group. `Error` instead terminates at the first
unsupported interaction. This prevents Palace from silently combining incompatible
local responses without discarding unrelated edge neighborhoods.

A version-1 process library has the following form:

```json
{
  "Version": 1,
  "Name": "100nm-Al-50nm-overetch",
  "MatchingRadius": 2.0,
  "CouponDepth": 1055.0,
  "Models": [
    {
      "Name": "isolated-edge",
      "Topology": "IsolatedEdge",
      "FabricatedMatrix": "isolated/fabricated/domain-response-matrix.csv",
      "ThinMatrix": "isolated/thin/domain-response-matrix.csv",
      "FabricatedSurfaceMatrix":
        "isolated/fabricated/surface-response-matrix.csv",
      "ThinSurfaceMatrix": "isolated/thin/surface-response-matrix.csv",
      "BasisPoints": "isolated/basis-points.csv",
      "Interfaces": [
        {"Type": "SA", "Coupon": 1},
        {"Type": "MS", "Coupon": 2},
        {"Type": "MA", "Coupon": 3}
      ]
    },
    {
      "Name": "same-conductor-gap-2um",
      "Topology": "SameConductorGap",
      "Separation": 2.0,
      "SeparationTolerance": 0.01,
      "Reference": [-1.0, 0.0, 0.0],
      "FabricatedMatrix": "gap-2um/fabricated/domain-response-matrix.csv",
      "ThinMatrix": "gap-2um/thin/domain-response-matrix.csv",
      "BasisPoints": "gap-2um/basis-points.csv"
    }
  ]
}
```

All lengths use the mesh coordinate unit of the application config. Matrix and basis
paths are resolved relative to the library file. `CouponDepth` is the implicit extrusion
depth represented by each 2D response matrix and should equal `Model.Lc` from the coupon
simulation. Palace divides longitudinal 3D quadrature weights by this depth. It also
uses the value to transfer a 2D coupon between applications with different
characteristic lengths. `CouponDepth` is required for 3D; a model-level value overrides
the library default. Omitting it preserves the legacy 2D assumption that coupon and
application use the same `Lc`.

Paired-edge models require a positive `Separation`. Palace prefers the closest model
within its nonnegative `SeparationTolerance`. Otherwise it uses the narrowest bracket
formed by lower- and upper-separation entries of the same topology. It combines the two
responses and conductor references with linear convex weights. A
`DifferentConductorGap` bracket cannot mix a legacy one-reference entry with a
two-conductor-state entry. Palace never extrapolates beyond the library range. The
bracket span divided by `MatchingRadius` contributes to `max library distance`, so a
sparse but technically complete bracket can still fail confidence.

`Reference` is the local coupon point whose sampled potential is subtracted from every
contour coefficient. Surface matrices and `Interfaces` are optional as a set, but are
required for corrected participation output.

A version-2 library can give a `DifferentConductorGap` coupon an independent voltage
state for its second conductor:

```json
{
  "Version": 2,
  "MatchingRadius": 2.0,
  "Models": [
    {
      "Name": "different-conductor-gap-2um",
      "Topology": "DifferentConductorGap",
      "Separation": 2.0,
      "SeparationTolerance": 0.01,
      "ConductorReferences": [
        [-3.0, 0.0, 0.0],
        [ 3.0, 0.0, 0.0]
      ],
      "OpenContourPaths": [
        {
          "Indices": [7, 8, 1, 2],
          "StartConductor": 1,
          "EndConductor": 2
        },
        {
          "Indices": [6, 5, 4, 3],
          "StartConductor": 1,
          "EndConductor": 2
        }
      ],
      "FabricatedMatrix": "gap/fabricated/domain-response-matrix.csv",
      "ThinMatrix": "gap/thin/domain-response-matrix.csv",
      "BasisPoints": "gap/basis-points.csv"
    }
  ]
}
```

The two references must be distinct points on conductors A and B in coupon coordinates.
For ``N`` free contour knots, the matrices are ``(N+1)`` by ``(N+1)`` and their
coefficient ordering is

```math
[V(\bm{x}_1)-V_A,\ldots,V(\bm{x}_N)-V_A,V_B-V_A].
```

Generate the last row and column by appending the conforming `TerminalAttributes` coupon
excitation described above after the free contour basis excitations. A version-2 library
may also contain legacy `DifferentConductorGap` entries without `ConductorReferences`,
which retain the original single-reference, contour-only approximation. This permits
version-1 and version-2 model libraries to be combined without rewriting their matrices.

Version 3 records the dielectric layers used when generating the surface-response
matrices:

```json
{
  "Version": 3,
  "Fabrication": {
    "InterfaceLayers": {
      "SA": {"Thickness": 0.002, "Permittivity": 4.0},
      "MS": {"Thickness": 0.002, "Permittivity": 11.47},
      "MA": {"Thickness": 0.002, "Permittivity": 10.0}
    }
  }
}
```

Thickness uses the application mesh coordinate unit, like `MatchingRadius`. Palace
requires an entry for every interface type represented by a version-3 model. Before
automatic 2D or 3D matching, it also requires every selected runtime interface to have
the same type, thickness, and relative permittivity. A mismatch is fatal because the
surface-response matrices already include both layer properties. Version-1 and
version-2 libraries remain usable, but Palace warns that their interface properties
cannot be verified when they omit this metadata.

`OpenContourPaths` is required when conductor cuts interrupt the free matching contour.
Each one-based `Indices` list orders a disjoint subset of `BasisPoints` from
`StartConductor` to `EndConductor`; the paths must partition all free contour points.
Their conductor graph must also connect every `ConductorReferences` entry. The coupon
generator emits this metadata automatically. Electrostatics samples the same point
coefficients directly. Maxwell uses a spanning tree of open paths to reconstruct every
conductor voltage relative to conductor 1, then reconstructs the contour knots from
their starting conductor. Additional paths close independent loops and contribute to
the loop-residual diagnostic. No integration segment crosses an omitted metal cut.

Three or more nearby collinear edges use a `ParallelEdgeCluster` entry:

```json
{
  "Name": "cpw-four-edge-cluster",
  "Topology": "ParallelEdgeCluster",
  "EdgeOffsetTolerance": 0.01,
  "Edges": [
    {"Offset": 0.0,  "GapDirection":  1, "Conductor": 1},
    {"Offset": 12.0, "GapDirection": -1, "Conductor": 2},
    {"Offset": 32.0, "GapDirection":  1, "Conductor": 2},
    {"Offset": 44.0, "GapDirection": -1, "Conductor": 3}
  ],
  "ConductorReferences": [
    [0.0, 0.0, 0.0],
    [12.0, 0.0, 0.0],
    [44.0, 0.0, 0.0]
  ],
  "OpenContourPaths": [
    {"Indices": [1, 2], "StartConductor": 1, "EndConductor": 2},
    {"Indices": [3, 4], "StartConductor": 2, "EndConductor": 3}
  ]
}
```

An exact cluster entry is self-contained over its matched longitudinal span. The
library does not also need redundant pair entries for every pair of edges inside
that span; pair coupons are required only on residual spans where exactly two
edges interact.

`Offset` starts at zero and increases along the canonical in-plane cluster axis.
`GapDirection` is the sign of the local metal-to-dielectric direction along that axis.
`Conductor` labels are contiguous and assigned in order of first occurrence. Palace
matches edge count, directions, conductor ownership, and every offset within
`EdgeOffsetTolerance`. In 3D it applies the cluster over each longitudinal span where
the active close-edge set is constant. With ``N`` contour knots and ``C`` conductor
references, response matrices are ``(N+C-1)`` by ``(N+C-1)``.

A localized three-dimensional interaction among nonparallel or endpoint-adjacent edges
uses a `SpatialEdgeCluster` entry:

```json
{
  "Name": "offset-corner-pair",
  "Topology": "SpatialEdgeCluster",
  "EdgePositionTolerance": 0.01,
  "EdgeAngleTolerance": 1.0,
  "Edges": [
    {
      "Point": [0.0, 0.0, 0.0],
      "GapDirection": [0.0, -1.0, 0.0],
      "ProcessNormal": [0.0, 0.0, 1.0],
      "Interval": [0.0, 2.0],
      "Conductor": 1,
      "BoundaryCondition": "PEC"
    },
    {
      "Point": [0.0, 0.0, 0.0],
      "GapDirection": [1.0, 0.0, 0.0],
      "ProcessNormal": [0.0, 0.0, 1.0],
      "Interval": [-2.0, 0.0],
      "Conductor": 1,
      "BoundaryCondition": "PEC"
    },
    {
      "Point": [1.25, -1.25, 0.0],
      "GapDirection": [-1.0, 0.0, 0.0],
      "ProcessNormal": [0.0, 0.0, 1.0],
      "Interval": [-2.0, 0.0],
      "Conductor": 2,
      "BoundaryCondition": "PEC"
    },
    {
      "Point": [1.25, -1.25, 0.0],
      "GapDirection": [0.0, 1.0, 0.0],
      "ProcessNormal": [0.0, 0.0, 1.0],
      "Interval": [0.0, 2.0],
      "Conductor": 2,
      "BoundaryCondition": "PEC"
    }
  ],
  "ConductorReferences": [
    [0.0, 0.0, 0.0],
    [1.25, -1.25, 0.0]
  ]
}
```

`Point`, `GapDirection`, and `ProcessNormal` define each edge in one rigid coupon
frame. The vectors must be orthonormal unit vectors. `Interval` gives the signed portion
of the physical edge, along
``\bm{t}=\text{GapDirection}\mathbin{\times}\text{ProcessNormal}``, replaced by the
spatial coupon. `BoundaryCondition` uses the same string or object form described below.
`EdgeAngleTolerance` is in degrees. `InterfaceSlot` is an optional nonnegative integer
and defaults to zero. Edges belonging to one physical fabrication layer use the same
slot.

Palace forms interaction events from the minimum-distance approach of each pair of
physical edge chains and groups nearby events into local neighborhoods. It matches the
complete edge count using one rigid transformation, a one-to-one edge assignment, and a
bijection between canonical and physical conductor labels. Every point, orientation,
process normal, and boundary law must match within tolerance. It does not
interpolate or extrapolate spatial clusters. A successful model creates one spatial
response patch, trims each declared interval from the straight-edge response, and
suppresses corner or endpoint models inside the replaced neighborhood. A mismatch
omits only that local interaction neighborhood under `Warn`.

For a cross-layer coupon, each `Interfaces` entry also declares the corresponding slot:

```json
"Interfaces": [
  {"Slot": 1, "Type": "SA", "Coupon": 1},
  {"Slot": 2, "Type": "SA", "Coupon": 2}
]
```

Palace globally matches the exact spatial cluster, binds each model slot to the complete
target-interface signature on its physical edges, and maps the slot's coupon output to
those target indices. This permits repeated interface types such as separate L1 and L2
SA outputs without additional correction-specific boundary splitting. Slot-to-signature
matching is one-to-one, and the coupon interface mappings in every slot must cover all
corresponding physical target-interface types. Extra compatible type aliases remain
permitted so one coupon can support different selected subsets. Absent, incomplete,
ambiguous, or geometrically mismatched cross-layer coupons retain the normal
`UnmatchedPolicy` behavior. `Slot` defaults to zero in `Interfaces`, preserving existing
single-group libraries.

A three-dimensional spatial-vertex coupon is added to the same library. For example, a
corner model has the form:

```json
{
  "Name": "convex-corner-90deg",
  "Topology": "ConvexCorner",
  "Angle": 90.0,
  "AngleTolerance": 2.0,
  "Reference": [0.0, 0.0, 0.0],
  "FabricatedMatrix": "corner-90/fabricated/domain-response-matrix.csv",
  "ThinMatrix": "corner-90/thin/domain-response-matrix.csv",
  "FabricatedSurfaceMatrix":
    "corner-90/fabricated/surface-response-matrix.csv",
  "ThinSurfaceMatrix": "corner-90/thin/surface-response-matrix.csv",
  "BasisPoints": "corner-90/basis-points.csv",
  "ContourGroups": [16, 16],
  "ZeroTraceIndices": [17, 18, 19],
  "Interfaces": [
    {"Type": "SA", "Coupon": 1},
    {"Type": "MS", "Coupon": 2},
    {"Type": "MA", "Coupon": 3}
  ]
}
```

`Angle` and `AngleTolerance` are in degrees; the angle is between the two incident arm
directions and must lie strictly between zero and 180 degrees. The topology distinguishes
an exterior metal corner from a reentrant corner. `BasisPoints` contains spatial
``(u,v,w)`` coordinates in an orthonormal frame whose first axis follows one oriented
incident arm, whose second axis lies in the fabrication plane, and whose third axis is
the process normal. The other arm lies at the model's `Angle` in this frame.
`ContourGroups` partitions those points into independent closed Maxwell voltage
contours. It is optional for electrostatics and defaults to one contour containing all
basis points. `ZeroTraceIndices` optionally lists one-based contour knots constrained to
the reference-conductor potential by PEC. Palace fixes every listed Maxwell trace to
zero and reconstructs each free contour arc independently between consecutive PEC knots.
It creates no line functional on an arc whose endpoints and interior knots are all PEC
constrained. This is required when a rounded-corner `Reference` lies inside the thin PEC
footprint and either the direct anchor path or part of the matching contour lies on the
conductor. Corner matrices are generated by paired thin and fabrication-resolved 3D
coupons and already represent one complete vertex, so Palace applies them once with no
`CouponDepth` scaling.

Endpoint and junction coupons use the same spatial basis:

```json
{
  "Name": "four-arm-junction",
  "Topology": "Junction",
  "ArmAngles": [0.0, 90.0, 180.0, 270.0],
  "ArmAngleTolerance": 2.0,
  "Reference": [0.0, 0.0, 0.0],
  "FabricatedMatrix": "junction/fabricated/domain-response-matrix.csv",
  "ThinMatrix": "junction/thin/domain-response-matrix.csv",
  "BasisPoints": "junction/basis-points.csv",
  "ContourGroups": [16, 16]
}
```

Every non-`SpatialEdgeCluster` model may set `BoundaryCondition`; the default is `PEC`.
The recommended object forms are:

```json
{"Type": "PEC"}
{"Type": "Conductivity", "Conductivity": 5.8e7, "Permeability": 1.0,
 "Thickness": 1.0e-7, "External": false}
{"Type": "Impedance", "Rs": 0.0, "Ls": 1.0e-13, "Cs": 0.0}
{"Type": "RationalImpedance",
 "Numerator": [5.0e-8, 0.0],
 "Denominator": [5.0e-20, 1.0e-9, 50.0]}
```

These fields have the same SI units and defaults as the corresponding entries in
`Boundaries`: conductivity in S/m, thickness in m, relative permeability dimensionless,
and surface `Rs`, `Ls`, and `Cs` in ohm/sq, H/sq, and F/sq. Rational numerator and
denominator coefficients are in highest-degree-first order for
``Z_s(s)=N(s)/D(s)``, with ``s`` in rad/s and ``Z_s`` in ohm/sq. Palace removes leading
zeros and a common polynomial scale before comparison. It compares conductivity,
permeability, and effective thickness; all three RLC parameters; or the complete
canonical rational function with a relative tolerance of ``10^{-10}``. Unknown object
keys are rejected.

The legacy strings `PEC`, `Conductivity`, `Impedance`, and `RationalImpedance` remain
accepted. A non-PEC string verifies only the class, not its numerical parameters, and
therefore fails the `boundary-law parameters verified` confidence gate. `PEC` has no
numerical parameters and is fully verified in either form.

Palace selects only a model whose boundary law matches every incident metal arm. A
vertex or parallel-edge neighborhood whose arms have different laws is left unmatched.
A `SpatialEdgeCluster` declares the law separately on every `Edges` entry, permitting an
exact mixed-boundary or mixed-parameter coupon.

An `Endpoint` model omits `ArmAngles`; its first local axis points away from the vertex
along the incident edge, its second local axis points from metal into the local gap, and
its third axis is the process normal. A `Junction` model requires at least three strictly
increasing `ArmAngles` in degrees, beginning with zero. Palace compares the complete
detected arm signature up to cyclic rotation while preserving the handedness fixed by
the process normal. Its first local axis follows the matched zero-degree arm, its second
axis lies in the fabrication plane, and its third axis is the process normal.
`ArmAngleTolerance` is applied independently to every arm. A reflected junction
signature therefore requires its own library entry when it is not rotationally
equivalent to the original one.

Like corner matrices, endpoint and junction matrices represent one complete vertex and
are not scaled by `CouponDepth`. Palace trims every incident arm by ``R`` before
integrating the remaining straight-edge response. Spatial endpoint and junction models
support PEC and non-PEC Maxwell metal. A non-PEC spatial model cannot use
`ZeroTraceIndices`; its voltage reference is the physical vertex because the sheet is
not equipotential. A non-PEC rounded-corner model must instead place
`Reference` on the resolved curved metal surface, since the reconstructed virtual sharp
corner and the fillet center are not physical surface points.

A model may contain more than one interacting edge. This is required when two one-edge
matching regions would overlap, for example when the width between two edges is less
than ``2R`` for matching radius ``R``. Palace rejects overlapping patch contours instead
of adding two independent corrections over the same region. In that case, generate one
coupled coupon whose contour encloses all interacting edges and place it as a single
patch. Independent one-edge coupons remain reusable wherever their patch contours do not
overlap.

The correction is applied matrix-free, while the assembled thin-metal Laplace operator
remains the preconditioner. Prescribed terminal values are included through consistent
right-hand-side elimination. Explicit `Models` and `Patches` placement supports 2D
electrostatic meshes. Automatic `Library` matching supports both 2D and 3D electrostatic
meshes. In 3D, every quadrature patch contributes

```math
\frac{ds}{L_{coupon}}\,
P(s)^T (S_{fabricated} - S_{thin}) P(s),
```

where ``L_{coupon}`` is `CouponDepth`.

Palace reports every 3D corner, endpoint, or junction for which no vertex model matched.
It integrates the straight-edge response through an unmatched neighborhood so that the
complete thin core is still replaced consistently. This fallback is not a corner model.
The reported unmodeled vertex-neighborhood fraction excludes successfully matched
spatial-vertex neighborhoods.

Palace solves both the original thin-metal system and the response-corrected system.
The standard saved fields, terminal matrices, `domain-E.csv`, `surface-Q.csv`, and AMR
indicator retain the original thin-metal solution. Enabling response correction
therefore does not change the historical raw result or mesh-refinement sequence.

When the models provide `FabricatedSurfaceMatrix`, `ThinSurfaceMatrix`, and `Interfaces`,
Palace additionally writes `surface-Q-corrected.csv`. This file contains the raw and
corrected total electric energies and, for every configured interface, both raw and
corrected surface energies, participation ratios, and quality factors. The raw columns
are the same quantities reported by the standard output. The corrected interface energy
replaces the complete thin coupon patch:

```math
\mathcal{E}^{surf}_{corr}
= \mathcal{E}^{surf}_{thin,outside\ coupon}
+ \bm{a}^T Q^{surf}_{fabricated,total} \bm{a},
```

where Palace uses the largest configured `EdgeDistances` value as ``R``. Each target
interface therefore requires an edge source and `EdgeDistances` containing the coupon
matching radius. The corrected total electric energy adds the matched
fabricated-minus-thin domain-response defect.

`EdgeExcludeAttributes` may be used to remove perimeter segments that lie on artificial
boundaries, such as the front and back faces of an extruded geometry. In 3D, the
perimeter is represented by straight segments joining the endpoints of the selected
surface mesh edges; curved high-order edges are therefore approximated by their
piecewise-linear mesh geometry.

### Maxwell response correction

Driven, eigenmode, and two-dimensional boundary-mode simulations can apply the same
fabrication-process library to a complex Maxwell field:

```json
"Solver": {
  "SurfaceResponseCorrection": {
    "Library": "process-library.json",
    "TargetInterfaces": [1, 2, 3],
    "UnmatchedPolicy": "Error"
  }
}
```

The standard `surface-Q.csv`, saved fields, error indicator, and AMR sequence remain the
raw thin-metal result. `surface-Q-corrected.csv` is an additional long-form table
containing raw, fixed-trace, fixed-flux, and, when available, self-consistent
normalization energies, interface energies, participations, and quality factors for each
frequency and excitation, eigenmode, or boundary mode.

In 3D, Palace maps the selected coupon contour into the plane transverse to the physical
edge at every longitudinal edge quadrature point. A corner coupon instead maps its full
spatial basis at one vertex. In a 2D boundary-mode solve, each detected cross-section
edge receives one planar coupon patch. Palace reconstructs each complex
quasi-electrostatic contour voltage using finite-element Nedelec line integrals. The
first knot of every closed contour is referenced to an automatically placed point on the
metal:

```math
a_i = -\int_{\bm{x}_{metal}}^{\bm{x}_i} \bm{E}\cdot d\bm{l}.
```

For PEC, an isolated-edge or gap anchor lies one matching radius inside the metal, and a
same-conductor strip anchor lies inside the strip. For a conductivity, impedance, or
rational-impedance boundary, the anchor is the local physical metal-surface point because
the tangential electric field need not vanish inside a finite-impedance sheet model.
Palace obtains the remaining coefficients by integrating around the coupon contour.
Fixing this reference is essential: contour differences alone leave an arbitrary
constant which changes the coupon response. The complex coefficients are evaluated with
Hermitian domain- and surface-response quadratic forms.

The corrected interface energy keeps the resolved raw energy outside the matching radius
and replaces the complete thin core with the fabricated coupon response. The
normalization energy receives the corresponding fabricated-minus-thin domain-response
defect. Every 3D target interface must therefore enable `AutomaticEdges`; every 2D
boundary-mode target must specify `EdgeAttributes` and `EdgeFrameNormal`. In both cases,
`EdgeDistances` must include the library `MatchingRadius`.

The Maxwell implementation supports:

  - a two-dimensional boundary-mode or three-dimensional driven/eigenmode simulation,
  - PEC, conductivity, impedance, or rational-impedance metal on isolated or parallel
    target edges, with a compatible boundary law in the coupon library,
  - any number of independent planar fabrication layers,
  - planar isolated-, paired-edge, or `ParallelEdgeCluster` coupon matches in 2D, or
  - isolated edges, paired edges, or longitudinal `ParallelEdgeCluster` spans in 3D,
    with optional spatial corner, endpoint, junction, and exact `SpatialEdgeCluster`
    coupons.

The paired-edge topology and conductor ownership are inferred without additional user
boundary splitting. In 2D boundary-mode meshes, Palace identifies PEC and
finite-impedance conductors from connected metal-boundary components, so disconnected
CPW grounds may share one boundary attribute without producing an anchor path through
the omitted metal.
`SameConductorGap`, `SameConductorStrip`, and
`DifferentConductorGap` library entries use the same separation matching as the
electrostatic path. Locally connected corner arms are handled by the matching corner
model or the unmodeled-corner fallback. Nonparallel close edges can use an exact
`SpatialEdgeCluster`; clusters without an exact library signature, interacting
cross-layer neighborhoods without an exact multi-slot library signature, and unbracketed
library separations follow `UnmatchedPolicy`. Non-PEC spatial-vertex models must declare
a compatible `BoundaryCondition` and cannot specify `ZeroTraceIndices`, because their
contour traces are not constrained by PEC. Object-form metadata verifies the numerical
boundary law against the simulation. Palace does not interpolate coupon matrices over
frequency, so a frequency-dependent non-PEC coupon must still be calibrated and
validated over its intended frequency range. Explicit `Models` or `Patches` are not
supported by the Maxwell path.

Before a field solve, the same production geometry classifier can audit process-library
coverage:

```text
palace --surface-response-preflight config.json
```

This loads, preprocesses, partitions, and applies a priori refinement to the solve mesh,
but does not assemble a Maxwell operator or solve a field problem. It writes
`surface-response-requirements.json` under `Problem.Output`. The deterministic manifest
aggregates equivalent neighborhoods and reports their topology, interface slots, metal
boundary law, count, total edge length, geometric signature, and exact, interpolated, or
missing library selection. Its `Complete` field is true only when no requirement is
missing. During preflight only, `UnmatchedPolicy: "Error"` is treated as `"Warn"` so all
missing requirements can be reported in one pass.

A multi-conductor coupon appends one ``V_i-V_1`` coefficient for every conductor after
the first. Electrostatics uses H1 point values. Maxwell integrates the declared open
paths and accumulates their signed voltage differences along a conductor spanning tree.
Any redundant path contributes an independent loop-closure check. Version-1
different-conductor coupons retain one PEC reference and no independent second
conductor voltage, so their use should be restricted to validated field families.
Same-conductor pairs and isolated PEC edges do not need an additional state.

The fixed-trace result applies the fabricated coupon matrices directly to the
reconstructed thin-field trace. The fixed-flux result transforms that trace by
``S_{fabricated}^{-1}S_{thin}`` before evaluating the fabricated response. Their relative
surface-energy difference, after summing each target interface, is the interface
`trace closure spread`; a large spread means a postprocessing-only correction is not
determined by the thin field. Palace also computes the spread independently for every
coupon patch and target interface. It reports the RMS of those local spreads weighted by
the larger-magnitude fixed-trace or fixed-flux response, and the fraction of that response
weight whose local spread exceeds 0.05. These patch-local diagnostics prevent opposite
closure errors in different neighborhoods from canceling in an interface total.
For a PEC coupon with `ZeroTraceIndices`, the fixed-flux equation is solved only on the
free trace subspace and the constrained coefficients remain exactly zero. Matrix rows
associated with those constrained calibration traces do not enter the physical closure.

For a uniform driven sweep, Palace also solves

```math
A_{corr}(\omega)
= A_{thin}(\omega)
- \omega^2 P^T(S_{fabricated}-S_{thin})P
```

with the raw thin operator as the preconditioner and the raw field as the initial guess.
The trace action ``P`` and its exact transpose use the same Nedelec line-integral
representation as postprocessing. Adaptive PROM sweeps do not currently provide a
self-consistent corrected field.

For a linear eigenmode problem without damping or frequency-dependent operators, Palace
also solves

```math
K\bm{E}_{corr}
= \omega_{corr}^2
\left[M + P^T(S_{fabricated}-S_{thin})P\right]\bm{E}_{corr}.
```

The raw eigenproblem remains unchanged and supplies the ordinary output and AMR error
indicator. The corrected eigenproblem uses the same target and the raw shifted operator
as its preconditioner. Palace pairs raw and corrected modes using a maximum-total-weight
one-to-one assignment of their normalized raw-``M`` overlaps. The corrected frequency,
paired-mode index, and overlap are reported separately; a self-consistent corrected
eigenmode requires overlap at least 0.8 to pass confidence. Damped and
frequency-dependent eigenproblems retain raw, fixed-trace, and fixed-flux output and
warn that no self-consistent corrected mode is available.

Boundary modes also report raw, fixed-trace, and fixed-flux results only. Their correction
uses the transverse field ``\bm{E}_t`` and is intended for quasi-TEM modes whose
longitudinal electric-field contribution is negligible. Palace does not solve a second
corrected boundary-mode eigenproblem.

Palace writes one row per frequency and excitation, eigenmode, or boundary mode to
`surface-response-confidence.csv`. `confidence pass` is one only when all of the
gating limits in the following table hold. `surface-Q-corrected.csv` also reports the
matched and unmodeled-corner fractions for each target interface; nontarget interfaces
receive `NaN` in those columns.

When a uniform driven or supported linear eigenmode solve provides a self-consistent
corrected field,
`self-consistent confidence pass` applies the same limits except for trace closure.
The global corrected solve supplies the exterior coupling which fixed-trace and
fixed-flux postprocessing leave undetermined. Its corrected-field contour diagnostics
are reported separately. A failed trace-closure gate therefore invalidates the two
postprocessing-only estimates, but does not by itself invalidate a self-consistent
result. Electrostatic correction has the same distinction: the fixed-trace/fixed-flux
spread is diagnostic of the postprocessed thin field, while the separately reported
corrected solution is globally self-consistent.

| Column | Meaning | Current limit |
|:---|:---|:---|
| `kR` | Electrical size of the coupon, ``|\omega|R/v_{min}`` | ``\leq 0.1`` |
| `loop residual` | Maximum ``|\oint \bm{E}\cdot d\bm{l}|`` divided by the sum of absolute contour-segment integrals over all contour paths | diagnostic only |
| `response-weighted loop residual` | RMS loop residual weighted by the fabricated surface response carried by each contour path | ``\leq 0.05`` |
| `loop-response fraction above limit` | Fraction of fabricated surface response assigned to paths whose loop residual exceeds 0.05 | ``\leq 0.01`` |
| `matched edge fraction` | Minimum, over target interfaces, of the selected physical edge length assigned a library model | effectively 1 |
| `unmodeled corner neighborhood fraction` | Maximum, over target interfaces, of the matched edge length within ``R`` of an unmatched corner, endpoint, or junction | ``\leq 0.1`` |
| `max R/rho` | Largest matching-radius to local-curvature-radius ratio | ``\leq 0.25`` |
| `max library distance` | Largest normalized paired-edge mismatch, separation-interpolation span, corner mismatch, or rounded-radius interpolation span | ``\leq 0.8`` |
| `boundary-law parameters verified` | One when every selected coupon declares PEC or object-form numerical boundary parameters matching the simulation | 1 |
| `max trace closure spread` | Largest relative difference between interface-aggregated fixed-trace and fixed-flux surface response | ``\leq 0.05`` |
| `response-weighted local trace closure spread` | RMS patch/interface closure spread weighted by the larger-magnitude fixed-trace or fixed-flux response | ``\leq 0.05`` |
| `trace-closure response fraction above limit` | Fraction of patch/interface response weight whose local closure spread exceeds 0.05 | ``\leq 0.01`` |
| `self-consistent M-overlap` | Raw-``M`` normalized overlap between paired raw and corrected eigenmodes | ``\geq 0.8`` (eigenmode only) |

The maximum loop residual identifies the worst local violation of the quasi-electrostatic
approximation or contour resolution. It does not gate confidence because an arbitrarily
small trace can have a large relative residual without affecting the corrected energy.
For response weighting, each coupon patch is assigned the larger of its fixed-trace and
fixed-flux fabricated surface energies, summed over target interfaces. That weight is
split among a patch's contour paths in proportion to the squared sum of absolute segment
integrals. The weighted RMS and failed-response fraction gate whether nonconservative
paths carry a material fraction of the modeled response. The remaining columns test
whether the geometry is
represented by the isolated-, paired-, or corner-response library. The interface-wise
extrema prevent a large-perimeter interface from hiding poor coverage of a smaller one.
The closure weighting uses each patch/interface response directly, before summation over
patches, so a low-response neighborhood does not dominate confidence while a material
local closure failure cannot be hidden by cancellation. Electrostatic
`surface-Q-corrected.csv` reports the same three closure diagnostics and confidence pass
for each source. Corrected values are still written when a limit fails, but Palace warns
that they are not validated. The unmodeled corner fraction is a geometric length
fraction, not a statistical confidence interval.

## Lumped parameter extraction

For electrostatic simulations, the Maxwell capacitance matrix is computed in the following
manner. First, the Laplace equation subject to Dirichlet boundary conditions is solved for
each terminal with boundary ``\Gamma_i`` in the model, yielding an associated voltage field
``V_i(\bm{x})``:

```math
\begin{aligned}
\nabla\cdot(\varepsilon_r\nabla V_i) &= 0 \,,\; \bm{x}\in\Omega \\
V_i &= 1 \,,\; \bm{x}\in\Gamma_i \\
V_i &= 0 \,,\; \bm{x}\in\Gamma_j \,,\; j\neq i \,.
\end{aligned}
```

The energy of the electric field associated with any solution is

```math
\mathcal{E}(V_i) = \frac{1}{2}\int_\Omega\varepsilon_r\bm{E}_i\cdot\bm{E}_i\,dV
```

where ``\bm{E}_i=-\nabla V_i`` is the electric field. Then, the entries of the Maxwell
capacitance matrix, ``\bm{C}``, are given by

```math
\bm{C}_{ij} = \mathcal{E}(V_i+V_j)-\frac{1}{2}(\bm{C}_{ii}+\bm{C}_{jj}) \,.
```

Magnetostatic problems for inductance matrix extraction are based on the magnetic vector
potential formulation:

```math
\begin{aligned}
\nabla\times(\mu_r^{-1}\nabla\times\bm{A}_i) &= 0 \,,\; \bm{x}\in\Omega \\
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{A}_i) =
    \bm{n}\times\bm{H}_i &= \bm{J}_s^{inc} \,,\; \bm{x}\in\Gamma_i \\
\bm{n}\times(\mu_r^{-1}\nabla\times\bm{A}_i) &= 0 \,,\; \bm{x}\in\Gamma_j \,,\; j\neq i \,.
\end{aligned}
```

For each port with boundary ``\Gamma_i``, a unit source surface current ``\bm{J}_s^{inc}``
is applied, yielding an associated vector potential solution ``\bm{A}_i(\bm{x})``.
Homogeneous Dirichlet boundary conditions ``\bm{n}\times\bm{A}_i=0`` are also imposed on
specified surfaces of the model. The magnetic field energy associated with any solution is

```math
\mathcal{E}(\bm{A}_i) = \frac{1}{2}\int_\Omega\mu_r^{-1}\bm{B}_i\cdot\bm{B}_i\,dV
```

where ``\bm{B}_i = \nabla\times\bm{A}_i`` is the magnetic flux density. Then, the entries of
the inductance matrix, ``\bm{M}``, are given by

```math
\bm{M}_{ij} = \frac{1}{I_i I_j}\mathcal{E}(\bm{A}_i+\bm{A}_j)
    - \frac{1}{2}\left(\frac{I_i}{I_j}\bm{M}_{ii}+\frac{I_j}{I_i}\bm{M}_{jj}\right)
```

where ``I_i`` is the excitation current for port ``i``, computed by integrating the source
surface current ``\bm{J}_s^{inc}`` over the surface of the port.

## Error estimation and adaptive mesh refinement (AMR)

Error estimation is used to provide element-wise error estimates for AMR, as well as a
global error indicator used to terminate AMR iterations or provide an estimate for solution
accuracy. A Zienkiewicz–Zhu (ZZ) error estimator based on [[5]](#References) is
implemented, which measures the error in the recovered magnetic field and electric flux
density. On element ``K``, we have

```math
\eta^2_K = \eta_{m,2}^2+\eta_{e,K}^2 =
    \|\mu_r^{1/2}\bm{R}_{ND}(\mu^{-1}\bm{B})
    - (\mu_r^{-1/2}\bm{B})\|_{L^2(\Omega_K)}^2
    + \|\varepsilon_r^{-1/2}\bm{R}_{RT}(\varepsilon_r\bm{E})
    - (\varepsilon_r^{1/2}\bm{E})\|_{L^2(\Omega_K)}^2
```

where ``\bm{R}_{ND}`` and ``\bm{R}_{RT}`` are the smooth-space recovery operators which
orthogonally project their argument onto ``H(\text{curl})`` and ``H(\text{div})``,
discretized by Nédélec and Raviart-Thomas elements, respectively.

## Far-field extraction

This feature is based upon Stratton-Chu's transformations [6] in the limit of ``kr \gg 1``
(with ``k`` wave number and ``r`` observation distance). One can show (see below) that, in
this limit,

```math
r \mathbf{E}_p(\mathbf{r}_0) = \frac{ik}{4\pi} \mathbf{r}_0 \times \int_S [\mathbf{n} \times \mathbf{E} - Z \mathbf{r}_0 \times (\mathbf{n} \times \mathbf{H})] \exp(ik\mathbf{r} \cdot \mathbf{r}_0) dS
```

where:

  - ``E_p`` is the electric field at the observation point
  - ``k`` is the wave number
  - ``r₀`` is the unit vector from source to observation point, parameterized by ``(\theta, \phi)``
  - ``n`` is the surface normal (to ``S``)
  - ``E, H`` are the tangential fields on the surface
  - ``Z`` is the impedance

The integral is over the exterior surface ``S``.

Note, we obtain ``r \mathbf{E}_p`` because the electric field decays with
``exp(ikr)/r``, so multiplying it by ``r`` ensures that the quantity is finite.
Note also that the solution is defined up to a global phase factor.

This equation relies on an analytic form for Green's function and is only valid
in 3D and if ``S`` only crosses isotropic materials.

From ``r \mathbf{E}_p``, one can obtain the magnetic field assuming that the
waves are propagating in free space,

```math
r \mathbf{H}_p = \frac{r_0 \times r \mathbf{E}_p}{Z_0}\,,
```

with ``Z_0`` impedance of free space.

With this, one can immediately compute the far-field relative radiation pattern
as ``|r \mathbf{E}_p|``.

#### How to obtain the equation above from Stratton-Chu's original equations

Let us start from Stratton-Chu's transformation for the electric field (we will get the magnetic field from ``E``):

```math
\mathbf{E}(\mathbf{r}_0) = \int_S \left[ i \omega \mu (\mathbf{n} \times \mathbf{H}) g(\mathbf{r}, \mathbf{r}_0) +
(\mathbf{n} \times \mathbf{E}) \times \nabla g(\mathbf{r}, \mathbf{r}_0) + (\mathbf{n} \cdot \mathbf{E}) \nabla g(\mathbf{r}, \mathbf{r}_0) \right] dS
```

with Green's function (here is where the assumption of isotropicity comes in):

```math
g(\mathbf{r}, \mathbf{r}_0) = \frac{e^{-i k |\mathbf{r} - \mathbf{r}_0|}}{4 \pi |\mathbf{r} - \mathbf{r}_0|}.
```

Let us take the limit for ``r \to \infty`` and define ``R = |\mathbf{r} - \mathbf{r}_0|`` (``R \to \infty`` when ``r \to \infty``).
For ``r \gg r_0`` (far-field approximation):

```math
R \approx r - \mathbf{r}\cdot\mathbf{r}_0
```

where ``\mathbf{r}_0 = \mathbf{r}/r`` is the unit vector in the direction of ``\mathbf{r}``.

The far-field approximation for Green's function becomes:

```math
g(\mathbf{r}, \mathbf{r}_0) \approx \frac{e^{-i k r}}{4 \pi r} e^{i k \mathbf{r}_0\cdot\mathbf{r}}.
```

For the gradient of ``g``, we start with the exact expression and expand phase and magnitude to reach:

```math
\nabla g(\mathbf{r}, \mathbf{r}_0) = -\frac{e^{-i k R}}{4 \pi R}\left(\frac{1}{R} + i k\right)\hat{R}
```

where ``\hat{R} = (\mathbf{r} - \mathbf{r}_0)/R`` is the unit vector pointing from ``\mathbf{r}_0`` to ``\mathbf{r}``.

In the far-field limit, ``R \approx r`` and ``\hat{R} \approx \mathbf{r}_0``, so:

```math
\nabla g(\mathbf{r}, \mathbf{r}_0) \approx -i k \mathbf{r}_0 g(\mathbf{r}, \mathbf{r}_0)
```

where we've neglected the ``1/R`` term since ``k R \gg 1`` in the far-field.

With these ingredients, one then uses the triple vector product rule and drops
the radial terms (i.e., those proportional to ``\mathbf{r}_0``, in the wave zone
there are only transverse fields) to arrive at the equation presented in the
previous section and implemented in *Palace*.

## References

[1] J.-M. Jin, _The Finite Element Method in Electromagnetics_, Wiley-IEEE Press, Hoboken,
NJ, Third edition, 2014.\
[2] P. Monk, _Finite Element Methods for Maxwell's Equations_, Oxford University Press,
Oxford, 2003.\
[3] L. Vardapetyan and L. Demkowicz, Full-wave analysis of dielectric waveguides at a given
frequency, _Mathematics of Computation_ 72 (2003) 105-129.\
[4] J. Wenner, R. Barends, R. C. Bialczak, et al., Surface loss of superconducting coplanar
waveguide resonators, _Applied Physics Letters_ 99, 113513 (2011).\
[5] S. Nicaise, On Zienkiewicz-Zhu error estimators for Maxwell’s equations, _Comptes Rendus
Mathematique_ 340 (2005) 697-702.\
[6] J. A, Stratton and L. J. Chu, Diffraction theory of Electromagnetic
Waves, _Physical Review_, 56, 1, (1939), 99-107.\
[7] P. Benner, D. C. Sorensen, and V. Mehrmann, Eds., _Dimension Reduction of
Large-Scale Systems_, Lecture Notes in Computational Science and Engineering,
vol. 45, Springer, 2005. doi: [10.1007/3-540-27909-1](https://doi.org/10.1007/3-540-27909-1).\
[8] A. C. Antoulas, _Approximation of Large-Scale Dynamical Systems_, SIAM,
2005. doi: [10.1137/1.9780898718713](https://doi.org/10.1137/1.9780898718713).\
[9] D. Pradovera, Toward a certified greedy Loewner framework with minimal
sampling, _Advances in Computational Mathematics_ 49, 92 (2023). doi:
[10.1007/s10444-023-10091-7](https://doi.org/10.1007/s10444-023-10091-7).\
[10] D. Pradovera, Interpolatory rational model order reduction of parametric
problems lacking uniform inf-sup stability, _SIAM Journal on Numerical
Analysis_ 58 (2020) 2265-2293. doi:
[10.1137/19M1269695](https://doi.org/10.1137/19M1269695).
