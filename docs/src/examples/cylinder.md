```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

```@setup include_example
function include_example_file(example_path, filename)
    print(read(joinpath(@__DIR__, "..", "..", "..", "test", "data", "regression", "ref", example_path, filename), String))
end
```

# Eigenmodes of a Cylinder

!!! note

    The files for this example can be found in the
    [`examples/cylinder/`](https://github.com/awslabs/palace/blob/main/examples/cylinder)
    directory of the *Palace* source code.

## Cavity

This example demonstrates *Palace*'s eigenmode simulation type to solve for the lowest
frequency modes of a cylindrical cavity resonator. In particular, we consider a cylindrical
cavity filled with Teflon (``\varepsilon_r = 2.08``,
``\tan\delta = 4\times 10^{-4}``), with radius ``a = 2.74\text{ cm}`` and height
``d = 2a``. From [[1]](#References), the frequencies of the ``\text{TE}_{nml}`` and
``\text{TM}_{nml}`` modes are given by

```math
\begin{aligned}
f_{\text{TE},nml} &= \frac{1}{2\pi\sqrt{\mu\varepsilon}}
    \sqrt{\left(\frac{p'_{nm}}{a}\right)^2 +
    \left(\frac{l\pi}{d}\right)^2} \\
f_{\text{TM},nml} &= \frac{1}{2\pi\sqrt{\mu\varepsilon}}
    \sqrt{\left(\frac{p_{nm}}{a}\right)^2 +
    \left(\frac{l\pi}{d}\right)^2} \\
\end{aligned}
```

where  ``p_{nm}`` and ``p'_{nm}`` denote the ``m``-th root (``m\geq 1``) of the ``n``-th
order Bessel function (``n\geq 0``) of the first kind, ``J_n``, and its derivative,
``J'_n``, respectively.

In addition, we have analytic expressions for the unloaded quality factors due to dielectric
loss, ``Q_d``, and imperfectly conducting walls, ``Q_c``. In particular,

```math
Q_d = \frac{1}{\tan\delta}
```

and, for a surface resistance ``R_s``,

```math
Q_c = \frac{(ka)^3\eta ad}{4(p'_{nm})^2 R_s}
    \left[1-\left(\frac{n}{p'_{nm}}\right)^2\right]
    \left\{\frac{ad}{2}
        \left[1+\left(\frac{\beta an}{(p'_{nm})^2}\right)^2\right] +
        \left(\frac{\beta a^2}{p'_{nm}}\right)^2
        \left(1-\frac{n^2}{(p'_{nm})^2}\right)\right\}^{-1}
```

where ``k=\omega\sqrt{\mu\varepsilon}``, ``\eta=\sqrt{\mu/\varepsilon}``, and
``\beta=l\pi/d``.

The initial Gmsh mesh for this problem, from
[`mesh/cavity_prism.msh`](https://github.com/awslabs/palace/blob/main/examples/cylinder/mesh/cavity_prism.msh),
is shown below. We use quadratic triangular prism elements. There are also two other
included mesh files,
[`mesh/cavity_tet.msh`](https://github.com/awslabs/palace/blob/main/examples/cylinder/mesh/cavity_tet.msh)
and
[`mesh/cavity_hex.msh`](https://github.com/awslabs/palace/blob/main/examples/cylinder/mesh/cavity_hex.msh),
which use curved tetrahedral and hexahedral elements, respectively.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cavity-1.png" width="60%" />
</p><br/>
```

There are two configuration files for this problem,
[`cavity_pec.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_pec.json)
and
[`cavity_impedance.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_impedance.json).

In both, the [`config["Problem"]["Type"]`](../config/reference.md#config-problem)
field is set to `"Eigenmode"`, and we use the mesh shown above. The material properties for
Teflon are entered under
[`config["Domains"]["Materials"]`](../config/reference.md#config-domains-materials). The
[`config["Domains"]["Postprocessing"]["Energy"]`](../config/reference.md#config-domains-postprocessing-energy)
object is used to extract the quality factor due to bulk dielectric loss; in this problem
since there is only one domain this is trivial, but in problems with multiple material
domains this feature can be used to isolate the energy-participation ratio (EPR) and
associated quality factor due to different domains in the model.

The only difference between the two configuration files is in the `"Boundaries"` object:
`cavity_pec.json` prescribes a perfect electric conductor (`"PEC"`) boundary condition to
the cavity boundary surfaces, while `cavity_impedance.json` prescribes a surface impedance
condition with the surface resistance ``R_s = 0.0184\text{ }\Omega\text{/sq}``, for copper
at ``5\text{ GHz}``.

In both cases, we configure the eigenvalue solver to solve for the ``15`` lowest frequency
modes above ``2.0\text{ GHz}`` (the dominant mode frequencies for both the
``\text{TE}`` and ``\text{TM}`` cases fall around ``2.9\text{ GHz}`` frequency for this
problem). A sparse direct solver is used for the solutions of the linear system resulting
from the spatial discretization of the governing equations, using in this case a
fourth-order finite element space.

The frequencies for the lowest-order ``\text{TE}`` and ``\text{TM}`` modes computed using
the above formula for this problem are listed in the table below.

| ``(n,m,l)`` | ``f_{\text{TE}}``       | ``f_{\text{TM}}``       |
|:----------- | -----------------------:| -----------------------:|
| ``(0,1,0)`` | ----                    | ``2.903605\text{ GHz}`` |
| ``(1,1,0)`` | ----                    | ``4.626474\text{ GHz}`` |
| ``(2,1,0)`` | ----                    | ``6.200829\text{ GHz}`` |
| ``(0,1,1)`` | ``5.000140\text{ GHz}`` | ``3.468149\text{ GHz}`` |
| ``(1,1,1)`` | ``2.922212\text{ GHz}`` | ``5.000140\text{ GHz}`` |
| ``(2,1,1)`` | ``4.146842\text{ GHz}`` | ``6.484398\text{ GHz}`` |
| ``(0,1,2)`` | ``5.982709\text{ GHz}`` | ``4.776973\text{ GHz}`` |
| ``(1,1,2)`` | ``4.396673\text{ GHz}`` | ``5.982709\text{ GHz}`` |
| ``(2,1,2)`` | ``5.290341\text{ GHz}`` | ``7.269033\text{ GHz}`` |

First, we examine the output of the `cavity_pec.json` simulation. The file
`postpro/cavity_pec/eig.csv` contains information about the computed eigenfrequencies and
associated quality factors:

```@example include_example
include_example_file("cylinder/cavity_pec", "eig.csv") # hide
```

Indeed we can find a correspondence between the analytic modes predicted and the solutions
obtained by *Palace*. Since the only source of loss in the simulation is the nonzero
dielectric loss tangent, we have ``Q = Q_d = 1/0.0004 = 2.50\times 10^3`` in all cases.

Next, we run `cavity_impedance.json`, which  adds the surface impedance boundary condition.
Examining `postpro/cavity_impedance/eig.csv` we see that the mode frequencies are roughly
unchanged but the quality factors have fallen due to the addition of imperfectly conducting
walls to the model:

```@example include_example
include_example_file("cylinder/cavity_impedance", "eig.csv") # hide
```

However, the bulk dielectric loss postprocessing results, computed from the energies written
to `postpro/cavity_impedance/domain-E.csv`, still give ``Q_d = 1/0.004 = 2.50\times 10^3`` for
every mode as expected.

Focusing on the ``\text{TE}_{011}`` mode with ``f_{\text{TE},010} = 5.00\text{ GHz}``, we
can read the mode quality factor ``Q = 2.30\times 10^3``. Subtracting out the contribution
of dielectric losses, we have

```math
Q_c = \left(\frac{1}{Q}-\frac{1}{Q_d}\right)^{-1} = 2.94\times 10^4
```

which is the same as the analytical result given in Example 6.4 from [[1]](#References) for
this geometry.

Finally, a clipped view of the electric field (left) and magnetic flux density magnitudes
for the ``\text{TE}_{011}`` mode is shown below.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cavity-2a.png" width="45%" />
  <img src="../../assets/examples/cavity-2b.png" width="45%" />
</p>
```

### Dispersive cavity material

The compact
[`cavity_drude.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_drude.json)
variant uses the hexahedral mesh with a damped Drude permittivity,

```math
\varepsilon_r(s) = \varepsilon_\infty +
    \frac{\omega_p^2}{s(s+\gamma)}, \qquad s=i\omega.
```

Here ``\varepsilon_\infty=2.08``, ``\omega_p=2\pi\times 1\text{ GHz}``, and
``\gamma=2\pi\times 0.1\text{ GHz}``. The
[`Permittivity`](../config/reference.md#config-domains-materials-permittivity) object uses a
`Drude` term whose `PlasmaFrequency` and `CollisionFrequency` are entered directly in GHz.

For a cavity mode with geometric angular frequency ``\Omega_0``, the expected root
satisfies ``(\Omega_0^2+\varepsilon_\infty s^2)(s+\gamma)+\omega_p^2s=0``. For the
lowest mode this predicts ``f=2.9852+0.0026947i\text{ GHz}``; the computed result is
``2.9866+0.0026921i\text{ GHz}``. The reported electric energy includes the kinetic
energy of the Drude carriers, evaluated at the real part of the mode frequency.

### Polaritons in a Lorentz-filled cavity

The
[`cavity_lorentz.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_lorentz.json)
variant fills the cavity with a damped Lorentz oscillator, the classical model of a bound
polarization resonance [[3]](#References),

```math
\varepsilon_r(s) = \varepsilon_\infty +
    \frac{\Delta\varepsilon\,\omega_0^2}{s^2+\gamma s+\omega_0^2}, \qquad s=i\omega,
```

with ``\varepsilon_\infty=2.08``, ``\Delta\varepsilon=0.01``, and
``\gamma=2\pi\times 0.02\text{ GHz}``. The resonance ``\omega_0`` is tuned to the
``\text{TM}_{010}`` frequency ``\omega_c = 2\pi\times 2.903636\text{ GHz}`` of the cavity filled
with ``\varepsilon_\infty`` alone, so that light and matter are degenerate.

The medium is homogeneous, so every cavity mode keeps its field profile and only its
frequency changes: a mode at ``\omega_c`` for ``\varepsilon_\infty`` requires
``\omega^2\varepsilon_r(\omega) = \varepsilon_\infty\omega_c^2``, or

```math
(\varepsilon_\infty\omega_c^2+\varepsilon_\infty s^2)(s^2+\gamma s+\omega_0^2)
    + \Delta\varepsilon\,\omega_0^2 s^2 = 0.
```

This quartic in ``s`` has two roots with positive frequency, so each cavity mode splits into
a pair of mixed light-matter modes, the polaritons of Hopfield [[4]](#References). Panel (b)
below shows the lossless case graphically, with ``\Delta\varepsilon`` exaggerated. ``\omega^2\varepsilon_r`` has a pole at
``\omega_0`` and a zero at the longitudinal frequency
``\omega_L = \omega_0\sqrt{1+\Delta\varepsilon/\varepsilon_\infty}``; in between,
``\varepsilon_r<0`` and no transverse mode exists, so one root lies below the gap and one
above. On resonance the roots are
``\omega_\pm = \tfrac{1}{2}\omega_0(\sqrt{\kappa+4}\pm\sqrt{\kappa})`` with
``\kappa=\Delta\varepsilon/\varepsilon_\infty``, so the splitting is exactly
``\omega_0\sqrt{\Delta\varepsilon/\varepsilon_\infty}``. Since ``\Delta\varepsilon`` is
proportional to the density ``N`` of oscillators, this is the ``\sqrt{N}`` scaling of the
vacuum Rabi splitting, which linear dispersion theory reproduces for any number of
oscillators [[5]](#References). Each polariton is half light and half matter on resonance,
so each inherits half the damping of the oscillator, whose own decay rate is
``\operatorname{Im}f = f_\gamma/2`` for ``\gamma = 2\pi f_\gamma``: each branch has
``\operatorname{Im}f\approx f_\gamma/4 = 5\text{ MHz}``.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cavity-4.svg" width="100%" />
</p><br/>
```

The ``\text{TE}_{111}`` pair at ``2.922197\text{ GHz}`` is only 0.6% above
``\text{TM}_{010}``, so it splits too. With `Target` at ``2.8\text{ GHz}``, the lowest four
modes are the lower branches of ``\text{TM}_{010}``, the ``\text{TE}_{111}`` pair, and
``\text{TM}_{011}``, all below the gap. The lower branches of all higher modes accumulate
below ``\omega_0``, so `N` counts them too when `Target` is close below the resonance. The
analytic values below use the uncoupled frequency ``f_c = \omega_c/2\pi`` both from the
closed form [[1]](#References) and as computed on the same mesh with
``\varepsilon_r=\varepsilon_\infty``. The first comparison includes the discretization error
of the cavity; the second removes it, leaving only the error in the dispersion model:

| Mode                | ``f_c``, exact (GHz) | ``f_c``, mesh (GHz) | *Palace* (GHz)         | Analytic, exact ``f_c`` | Analytic, mesh ``f_c`` |
|:------------------- | --------------------:| -------------------:| ----------------------:| -----------------------:| ----------------------:|
| ``\text{TM}_{010}`` | ``2.903636``         | ``2.905077``        | ``2.805548+0.004863i`` | ``2.804831+0.004827i``  | ``2.805548+0.004863i`` |
| ``\text{TE}_{111}`` | ``2.922197``         | ``2.926602``        | ``2.815652+0.005396i`` | ``2.813677+0.005288i``  | ``2.815652+0.005396i`` |
| ``\text{TM}_{011}`` | ``3.468175``         | ``3.469647``        | ``2.887844+0.009656i`` | ``2.887802+0.009654i``  | ``2.887844+0.009656i`` |

The script
[`cavity_polariton.jl`](https://github.com/awslabs/palace/blob/main/examples/cylinder/cavity_polariton.jl)
reproduces this table, then sweeps ``f_0`` from ``2.55`` to ``3.3\text{ GHz}``. At each
of 31 steps it solves once with `Target` below the lower ``\text{TM}_{010}`` branch and
once above the gap, for three modes each. All 186 computed modes are
``\text{TM}_{010}`` or ``\text{TE}_{111}`` branches and lie on the dispersion relation for
the mesh's uncoupled frequencies to a relative ``2\times 10^{-9}``. The branches
anticross around the uncoupled modes and approach the gap edges away from them, while the
damping passes from the matter-like branch to the light-like one, crossing at
``f_\gamma/4`` on resonance. The reported electric energy includes the polarization
oscillator's potential and kinetic energy, evaluated at the real part of the mode
frequency; it is not simply the energy at ``\operatorname{Re}\varepsilon_r``.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cavity-5.png" width="80%" />
</p><br/>
```

### Mesh convergence

The effect of mesh size can be investigated for the cylindrical cavity resonator using
[`convergence_study.jl`](https://github.com/awslabs/palace/blob/main/examples/cylinder/convergence_study.jl).
For a polynomial order of solution and refinement level, a mesh is generated using Gmsh
using polynomials of the same order to resolve the boundary geometry. The eigenvalue
problem is then solved for ``f_{\text{TM},010}`` and ``f_{\text{TE},111}``, and the
relative error, ``\frac{f-f_{\text{true}}}{f_{\text{true}}}``, of each mode plotted against
``\text{DOF}^{-\frac{1}{3}}``, a notional mesh size. Three different element types are
considered: tetrahedra, prisms and hexahedra, and the results are plotted below. The
``x``-axis is a notional measure of the overall cost of the solve, accounting for
polynomial order.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/cavity-3a.png" width="70%" />
  <img src="../../assets/examples/cavity-3b.png" width="70%" />
  <img src="../../assets/examples/cavity-3c.png" width="70%" />
</p><br/>
```

The observed rate of convergence for the eigenvalues are ``p+1`` for odd polynomials and
``p+2`` for even polynomials. Given the eigenmodes are analytic functions, the theoretical
maximum convergence rate is ``2p`` [[2]](#References). The figures demonstrate that
increasing the polynomial order of the solution will give reduced error, however the effect
may only become significant on sufficiently refined meshes.

## Waveguide

This example demonstrates the eigenmode simulation type in  *Palace* to solve for the
cutoff-frequencies of a circular waveguide. As with the cavity the interior material to be
Silicon (``\varepsilon_r = 2.08``, ``\tan\delta = 4\times 10^{-4}``), with cylindrical
domain radius ``a = 2.74\text{ cm}``, and length ``d=2a = 5.48\text{ cm}``, however now
periodic boundary conditions are applied in the $z$-direction. According to
[[1]](#References), the cutoff frequencies for the transverse electric and magnetic modes
are given by the formulae:

```math
\begin{aligned}
f_{\text{TE},nm} &= \frac{1}{2\pi\sqrt{\mu\varepsilon}} \frac{p'_{nm}}{a}\\
f_{\text{TM},nm} &= \frac{1}{2\pi\sqrt{\mu\varepsilon}} \frac{p_{nm}}{a}
\end{aligned}
```

which are identical to those for the cavity modes, in the special case of ``l=0``.

In addition to these pure waveguide modes, there are aliasing cavity
modes corresponding to a full wavelength in the computational domain (``l=2``). In a
practical problem these can be suppressed by choosing a smaller value of ``d`` which shifts
such modes to higher frequencies. The relevant modes are tabulated as

| ``(n,m,l)`` | ``f_{\text{TE}}``       | ``f_{\text{TM}}``       |
|:----------- | -----------------------:| -----------------------:|
| ``(0,1,0)`` | ``4.626481\text{ GHz}`` | ``2.903636\text{ GHz}`` |
| ``(1,1,0)`` | ``2.223083\text{ GHz}`` | ``4.626481\text{ GHz}`` |
| ``(2,1,0)`` | ``3.687749\text{ GHz}`` | ``6.200856\text{ GHz}`` |
| ``(3,1,0)`` | ``5.072602\text{ GHz}`` | ``7.703539\text{ GHz}`` |
| ``(0,1,2)`` | ``5.982715\text{ GHz}`` | ``4.776992\text{ GHz}`` |
| ``(1,1,2)`` | ``4.396663\text{ GHz}`` | ``5.982715\text{ GHz}`` |
| ``(2,1,2)`` | ``5.290372\text{ GHz}`` | ``7.269056\text{ GHz}`` |
| ``(3,1,2)`` | ``6.334023\text{ GHz}`` | ``8.586796\text{ GHz}`` |

For this problem, we use curved tetrahedral elements from the mesh file
[`mesh/cavity_tet.msh`](https://github.com/awslabs/palace/blob/main/examples/cylinder/mesh/cavity_tet.msh),
and the configuration files
[`waveguide.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/waveguide.json) and
[`floquet.json`](https://github.com/awslabs/palace/blob/main/examples/cylinder/floquet.json).

The main difference between these configuration files and those used in the cavity example is
in the `"Boundaries"` object: `waveguide.json` specifies a perfect electric conductor
(`"PEC"`) boundary condition for the exterior surface and a periodic boundary condition
(`"Periodic"`) on the cross-sections of the cylinder (in the $z-$ direction). The periodic
attribute pairs are defined by `"DonorAttributes"` and `"ReceiverAttributes"`, and the
distance between them is given by the `"Translation"` vector in mesh units. In `floquet.json`,
an additional `"FloquetWaveVector"` specifies the phase delay between the donor and receiver
boundaries in the X/Y/Z directions.

The file `postpro/waveguide/eig.csv` contains information about the computed eigenfrequencies and
associated quality factors:

```@example include_example
include_example_file("cylinder/waveguide", "eig.csv") # hide
```

In common with the PEC cavity ``Q = Q_d = 1/0.0004 = 2.50\times 10^3`` in all cases, and all
the anticipated waveguide modes are recovered with ``\text{TE}_{1,1}`` having the lowest
cutoff frequency followed by ``\text{TM}_{0,1}`` and ``\text{TE}_{2,1}``, while the aliasing
mode ``\text{TE}_{1,1,2}`` has marginally lower frequency than the waveguide modes
``\text{TE}_{0,1}`` and ``\text{TM}_{1,1}`` (``4.397\text{ GHz}`` compared to ``4.627\text{ GHz}``) and is thus found first.

## References

[1] D. M. Pozar, _Microwave Engineering_, Wiley, Hoboken, NJ, 2012.\
[2] A. Buffa, P. Houston, I. Perugia, Discontinuous Galerkin computation of the Maxwell
eigenvalues on simplicial meshes, _Journal of Computational and Applied Mathematics_ 204
(2007) 317-333.\
[3] J. D. Jackson, _Classical Electrodynamics_, 3rd ed., Wiley, New York, NY, 1999,
Section 7.5.\
[4] J. J. Hopfield, Theory of the contribution of excitons to the complex dielectric
constant of crystals, _Physical Review_ 112 (1958) 1555-1567.\
[5] Y. Zhu, D. J. Gauthier, S. E. Morin, Q. Wu, H. J. Carmichael, T. W. Mossberg, Vacuum
Rabi splitting as a feature of linear-dispersion theory: Analysis and experimental
observations, _Physical Review Letters_ 64 (1990) 2499-2502.
