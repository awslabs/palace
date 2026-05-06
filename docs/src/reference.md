```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Reference

*Palace* calculate solutions of Maxwell's equation for five different "Problem Types":

  - Electrostatics: time-independent voltage sources,
  - Magnetostatics: time-independent current sources,
  - Driven: monochromatic excitations ``\bm{U}^{inc}(f)``,
  - Eigenmode: resonances of the system,
  - Transient: time-dependant excitations ``\bm{U}^{inc}(t)``.

Here, we will summarize the mathematics of each, provide a description of boundary conditions, and
specify convention choices. This should be read in conjunction with the "User Guide" and
"Configuration File" documentation.

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

Ports have two effects. First, they impose a boundary condition that relates the electric and
magnetic fields on the surface. Second, for the driven and transient solver, they can excite the
system with an incident field.

For lumped ports, the boundary condition is an impedance condition as described above — tangential
fields satisfy ``\bm{E}_t = Z_s \bm{H}_t \times \bm{n}``. Here ``Z_s`` is the physical surface
impedance (per unit area) and ``\bm{n}`` the port normal. The excitation modes ``\bm{E}^{inc}`` are
analytically specified.

For wave ports, a solver performs eigenmode simulations on the boundary to find supported modes of a
transmission line with the wave-port cross section. The boundary modes will be used to relate
``\bm{E}`` and ``\bm{H}``, which will only be impedance-like for TEM, TE or TM modes. Additionally,
users can choose which of the found modes to excite.

We normalize port fields so the total power flow is ``\vert P^{inc} \vert = 1~\mathrm{W}``, where

[TODO: Right now this is peak power, not the average. This is the Palace convention, but needs
adjusting].

```math
P^{inc} = V^{inc} [I^{inc}]^*  = \sum_e \int_{\Gamma_e} dS_e \,  \bm{n}_e \cdot (\bm{E}_e^{inc} \times [\bm{H}_e^{inc}]^*).
```

Here we have allowed for the possibility that a port can be made of disjoint surface elements ``e``,
which separately contribute to the Poynting vector.

There is an additional normalization related to voltages and currents. For AC simulations of
microwave networks, circuit ``V`` and ``I`` are not generally well defined and have to be fixed by
convention [7,8]. This is tantamount to relating the circuit characteristic impedance

```math
Z = \frac{V^{inc}}{I^{inc}} = \frac{\vert V^{inc} \vert^2}{[P^{inc}]^*},
```

to the physical surface (wave) impedance ``Z_s``. We discuss this for each port type below.

### Lumped Ports

Lumped ports in *Palace* can be composed of two types of elements:

 1. *Rectangular*: The incident field is ``\bm{E}^{inc} = E_0 \, \hat{\bm{l}}``. Here ``E_0`` is a
    normalization constant and ``\hat{\bm{l}}`` is the user-defined polarization direction, which
    should point between two conductors. We denote the length of the port along the polarization as
    ``l`` and the perpendicular width as ``w``.

 2. *Coaxial*: ``\bm{E}^{inc} = {E_0 r_0} \hat{\bm{r}} / {r}``, where ``E_0 r_0`` is a
    normalization constant, ``r`` is the distance from the port center, and ``\hat{\bm{r}}`` is the
    unit radial vector. We denote the inner and outer radii as ``a`` and ``b``.

We define the voltage across each element as averages over the area, rather than line integrals. For
an electric field ``\bm{E}``, the voltage across element ``e`` is:

 1. *Rectangular*: ``V_{e} = \int_{\Gamma_e}dS\, \bm{E} \cdot\hat{\bm{l}}_e / w_e``,

 2. *Coaxial*: ``V_{e} = \int_{\Gamma_e}dS\, \bm{E} \cdot\hat{\bm{r}}_e / (2 \pi r)``.

This is a convention choice. For electrostatics or for TEM and TM transmission lines, the voltage in
the port plane is path-independent between different conductors; in general, it is not.

The power normalization and voltage convention now fix the coefficients ``E_0, r_0`` in the port
definition and relate the physical surface impedance to the circuit impedance. For each element:

 1. *Rectangular*: ``Z = Z_s l / w``,

 2. *Coaxial*: ``Z = Z_s \ln(b/a) / (2\pi)``,

which we also write as ``Z = Z_s / \alpha``  where ``\alpha = w/l`` or ``\alpha = 2\pi / \ln(b/a)``
respectively. If ``Z`` is a combination of LRC responses, these add in parallel

```math
\frac{1}{Z} = \frac{1}{R}+\frac{1}{i\omega L}+i\omega C.
```

Specifically, ``R_s = \alpha R``, ``L_s = \alpha L``, and ``C_s = C/\alpha``.

In the configuration file, a user specifies either ``L``, ``R``, ``C``, corresponding to ``Z``, or
``L_s`` ,``R_s``, ``C_s`` corresponding to ``Z_s``. For a single-element port, these will just be
converted using the single scale factor ``\alpha``. For a multi-element port, we require impedances
to add in parallel ``{1}/{Z} = \sum_e {1}/{Z_e}``, so that each element sees the same voltage.
Specifying a lumped ``Z`` in the config means that each element ``e`` can have a different physical
surface impedance depending on their shape ``\alpha_e``:

```math
Z_{s, e} = n_\mathrm{elem} \alpha_{e} Z.
```

Conversely, specifying ``Z_s`` means elements of different shape will have different ``Z_e``.

The source term ``\bm{U}^{inc}`` in a driven frequency-response problem is related to the
tangential component of the incident field at an excited port boundary by

```math
\bm{U}^{inc} = - 2 \gamma \bm{E}^{inc}_t = -2\gamma(\bm{n}\times\bm{E}^{inc})\times\bm{n}.
```

The excitation amplitude is fixed at the unit-power mode normalization described above
``1~\mathrm{W}``.

In the time-domain formulation, the source term ``\bm{U}^{inc}`` appears as

```math
\bm{U}^{inc} = -2 Z_s^{-1}\left(\bm{n}\times\frac{\partial\bm{E}^{inc}}{\partial t}\right)
    \times\bm{n} \,.
```

The incident field ``\bm{E}^{inc}(\bm{x},t)`` is

```math
\bm{E}^{inc}(\bm{x},t) = p(t)\bm{E}^{inc}(\bm{x})
```

where ``\bm{E}^{inc}(\bm{x})`` is identical to the spatial excitation in the frequency domain
formulation, and ``p(t)`` describes the temporal shape of the excitation. Possible options include a
sinusoidal, Gaussian, modulated Gaussian, or step excitation.

### Wave ports

Numeric wave ports assume a field with known normal-direction dependence
``\bm{E}(\bm{x}) = \bm{e}(\bm{x}_t)e^{ik_n x_n}`` where ``k_n`` is the propagation constant. For each operating
frequency ``\omega``, a two-dimensional eigenvalue problem is solved on the port yielding the mode
shapes ``\bm{e}_m`` and associated propagation constants ``k_{n,m}``. These are used in the full 3D
model where the Robin port boundary condition has coefficient
``\gamma = i\text{Re}\{k_{n,m}\}/\mu_r`` and the computed mode is used to compute the incident field in the
source term ``\bm{U}^{inc}`` at excited ports.

For more information on the implementation of numeric wave ports, see [[3]](#References).

To Do:

  - Discuss V / I / Z convention choice here.
  - Fix dembedding sign and normal sign convention (outward vs inward normal)

### Scattering Parameters

In the frequency domain, the scattering parameters can be postprocessed from the computed
electric field for each lumped port with boundary ``\Gamma_i`` as

```math
S_{ij} = \frac{\displaystyle\int_{\Gamma_i}\bm{E}\cdot\bm{E}^{inc}_i\,dS}
    {\displaystyle\int_{\Gamma_i}\bm{E}^{inc}_i\cdot\bm{E}^{inc}_i\,dS} - \delta_{ij} \,.
```

For wave ports, since the propagation constants are known, we allow additional de-embedding by
specifying an offset distance ``d`` for each wave port:

```math
\tilde{S}_{ij} = S_{ij}e^{ik_{n,i}d_i}e^{ik_{n,j}d_j} \,.
```

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

## Electrostatic Simulations

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

## Magnetostatic Simulations

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

[1] J.-M. Jin, *The Finite Element Method in Electromagnetics*, Wiley-IEEE Press, Hoboken,
NJ, Third edition, 2014.\
[2] P. Monk, *Finite Element Methods for Maxwell's Equations*, Oxford University Press,
Oxford, 2003.\
[3] L. Vardapetyan and L. Demkowicz, Full-wave analysis of dielectric waveguides at a given
frequency, *Mathematics of Computation* 72 (2003) 105-129.\
[4] J. Wenner, R. Barends, R. C. Bialczak, et al., Surface loss of superconducting coplanar
waveguide resonators, *Applied Physics Letters* 99, 113513 (2011).\
[5] S. Nicaise, On Zienkiewicz-Zhu error estimators for Maxwell’s equations, *Comptes Rendus Mathematique* 340 (2005) 697-702.\
[6] J. A, Stratton and L. J. Chu, Diffraction theory of Electromagnetic
Waves, *Physical Review*, 56, 1, (1939), 99-107.
[7] R. B. Marks and D. F. Williams, A general waveguide circuit theory, *J. Res. Natl. Inst. Stan.* 97, 533 (1992).\
[8] D. M. Pozar, *Microwave engineering*, Fourth edition. John Wiley & Sons, Hoboken, NJ, 2012.
