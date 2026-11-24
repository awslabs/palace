```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Simulation Models

## Supported mesh formats

The [`config["Model"]`](../config/reference.md#config-model) object is used to
specify the mesh for the discretized computational domain. In general, inputs are expected
to be dimensional nondimensionalized internally. A length scale, specified under
[`config["Model"]["L0"]`](../config/reference.md#config-model), describes the length
units of the mesh relative to 1 meter (i.e. `config["Model"]["L0"]: 1.0e-6` if the mesh
coordinates are in ``\mu``m, this is the default value). All other entries in the
configuration file which have units of length should be specified in units of
`config["Model"]["L0"]` m.

MFEM supports a [wide variety](https://mfem.org/mesh-formats/) of mesh formats. In
addition, *Palace* has built-in support for
[Nastran (`.nas`, `.bdf`)](https://docs.plm.automation.siemens.com/tdoc/scnastran/2020_1/help/#uid:index_element)
and
[COMSOL (`.mphtxt`, `.mphbin`)](https://doc.comsol.com/6.0/doc/com.comsol.help.comsol/COMSOL_ProgrammingReferenceManual.pdf)
format mesh files, for both linear and high-order curved elements.

Geometric attributes for domains and boundaries in the mesh are used to define material
properties and boundary conditions on the desired model regions and surfaces (see
[`config["Domains"]`](../config/reference.md#config-domains) and
[`config["Boundaries"]`](../config/reference.md#config-boundaries)). These attribute integers correspond to
tags for the domain and boundary elements in the mesh, and should be non-negative and start
at 1. They do not need to be contiguous in the mesh file. Throughout the configuration
file, the `"Attributes"` keyword is used to indicate which domain or boundary attributes
are relevant to the material properties or boundary conditions being specified.

## Mesh refinement

Refinement of the input mesh file can be performed using levels of global uniform refinement
or region-based refinement, specified using the
[`config["Model"]["Refinement"]`](../config/reference.md#config-model-refinement) object.
The user can specify any combination of uniform refinement levels as well as local
refinement regions which refines the elements inside of a certain box or sphere-shaped
region. For simplex meshes, the refinement maintains a conforming mesh but meshes
containing hexahedra, prism, or pyramid elements will be nonconforming after local
refinement (this is not supported at this time).

[Adaptive mesh refinement (AMR)](https://en.wikipedia.org/wiki/Adaptive_mesh_refinement)
according to error estimates calculated from the computed solution can also be specified
using the [`config["Model"]["Refinement"]`](../config/reference.md#config-model-refinement)
object. Nonconformal refinement is supported for all mesh types, and additionally conformal
refinement is supported for simplex meshes. AMR is available for all problem types apart
from driven problems in the time domain. The error estimates account for the
discontinuity of the fields across [interior boundaries](boundaries.md#Interior-boundaries)
such as thin metal sheets, without modifying the mesh.

AMR stops when the error indicator falls below `"Tol"`, after `"MaxIts"` iterations, or when
the last solved size exceeds `"MaxSize"` degrees of freedom, so the final solve can exceed
`"MaxSize"` by the growth of one refinement step. With `"MaxSizePredicted": true`, the size
is also checked after each refinement step and before its solve: the number of degrees of freedom of the refined mesh
is predicted as the last solved size scaled by the growth of the element count in that step
(the count of a fixed-order space scales with the element count; on graded nonconforming
meshes the prediction is within a few percent on fine meshes, while coarse early cycles
over-predict by up to ~15 % and stop conservatively), and when the prediction exceeds
`"MaxSize"` the solve is skipped and the last solved mesh is the final result (the log
reports the prediction and `palace.json` records it under
`"AdaptiveMeshRefinement"`/`"MaxSizeStop"`). This guarantees that no solve runs on a mesh
predicted larger than `"MaxSize"`; one Dörfler step can otherwise grow the size by a factor
2–3. When the prediction under-reads by a few percent the solve can still land slightly
above `"MaxSize"`, in which case the check against the last solved size ends AMR as before.
The predictive check is off by default (the rule of upstream Palace); the configuration
writers of the fabrication-response workflow under `examples/` set it explicitly.

An element-size floor, `"MinElementSize"` (a length in mesh units, scaled by `"L0"`; the
default `0` disables it), keeps AMR from chasing a geometric singularity to element sizes
that the linear solver can no longer resolve well. With a floor, an element is eligible for
marking only when its size `h` is at least twice the floor, so that one refinement of a
marked element stays at or above the floor; `h` is the volume-based element size whose
minimum and maximum the mesh statistics in the log report as `h`. The Dörfler marking then
selects the `"UpdateFraction"` of the error among the eligible elements only, while the
error estimate and the `"Tol"` test keep every element. Each AMR iteration logs the floor,
the number of ineligible elements with their share of the error, and, after the refinement
step, the number of elements below the floor (the conformity closure of a nonconforming or
conforming refinement can create elements smaller than a marked element's children; they
are counted, not prevented). When no element is eligible, or every eligible element has a
zero error indicator, AMR stops before the refinement step and the last solved mesh is the
final result.

## Material models

Material properties are handled by the
[`config["Domains"]["Materials"]`](../config/reference.md#config-domains-materials)
object. *Palace* supports linear material models and scalar frequency-dependent
permittivity for supported frequency-domain simulations.

Materials with scalar or symmetric matrix-valued material properties are supported. For most
simulation types, each material in the model requires a specified relative permittivity and
relative permeability (for electrostatic simulations, only the permittivity is required,
while for magnetostatics, only the permeability is required). For dielectric domains, a
loss tangent may be specified. Alternatively, for normal conducting domains, an electrical
conductivity may be specified which is used to relate the current density and electric
field via Ohm's law.

Modeling of superconducting domains is performed using the current-field constitutive
relations given by the London equations. The user can specify a London penetration depth to
activate this model. It can also be used in conjunction with a material conductivity when
wishing to model both superconducting and normal currents. In magnetostatic simulations, the
supercurrent stores kinetic energy, which is included in the extracted inductance (see
[Magnetostatic problems](problem.md#Magnetostatic-problems)).

### Frequency-dependent permittivity

Scalar additive permittivity models use the object form

```json
"Permittivity": {
  "HighFrequencyLimit": 2.08,
  "Terms": [
    {"Type": "Drude", "PlasmaFrequency": 1.0, "CollisionFrequency": 0.1}
  ]
}
```

`HighFrequencyLimit` is the scalar ``\epsilon_\infty`` and `Terms` is a nonempty additive array,
so that ``\epsilon(s) = \epsilon_\infty + \sum_k \chi_k(s)`` with ``s = i\omega`` (*Palace*
uses an ``e^{i\omega t}`` time dependence). Frequency inputs are ordinary frequencies ``f`` in
GHz, converted to angular rates as ``\omega = 2\pi f``. In particular, `CollisionFrequency`
and `DampingFrequency` give ``\gamma = 2\pi f``, so a damping specified as an inverse time
``1/\tau`` corresponds to ``f = 1/(2\pi\tau)``. The available terms are:

  - `Drude`: positive `PlasmaFrequency` and `CollisionFrequency`, both in GHz, with
    ``\chi(s) = \omega_p^2/(s(s+\gamma))``.
  - `Debye`: positive `DeltaPermittivity` and `RelaxationTime` ``\tau`` in ns, with
    ``\chi(s) = \Delta\epsilon/(1+s\tau)``.
  - `Lorentz`: positive `DeltaPermittivity` and `ResonanceFrequency`, and nonnegative
    `DampingFrequency`, both frequencies in GHz, with
    ``\chi(s) = \Delta\epsilon\,\omega_0^2/(s^2+\gamma s+\omega_0^2)``. The critically
    damped case ``\gamma = 2\omega_0``, to within a relative ``10^{-6}``, is not supported.
  - `PoleResidue`: a general term ``r/(s-p)`` with `Pole` and `Residue` in SI angular rates
    (rad/s). Values are real or `[real, imag]`, with ``\operatorname{Re}p \le 0``,
    ``\operatorname{Im}p \ge 0``, and a nonzero residue; a real pole needs a real residue.
    A pole in the upper half-plane includes its conjugate term automatically. Terms with
    the same pole are one term with the sum of their residues, which must not cancel
    exactly.
  - `DjordjevicSarkar`: positive `Strength`, `LowerFrequency`, and
    `UpperFrequency` bounds in GHz, with the upper bound larger than the lower. This native
    term is ``\mathrm{Strength}\,\ln((s+\omega_\mathrm{hi})/(s+\omega_\mathrm{lo}))``, with
    the natural logarithm, so a static increment ``\Delta\epsilon`` corresponds to
    ``\mathrm{Strength} = \Delta\epsilon/\ln(\omega_\mathrm{hi}/\omega_\mathrm{lo})``.

The `PoleResidue` rates are in rad/s, unlike the GHz and ns of the typed terms. For example,
the Debye term ``\Delta\epsilon/(1+s\tau)`` is ``r/(s-p)`` with ``p = -1/\tau`` and
``r = \Delta\epsilon/\tau``: with ``\Delta\epsilon = 0.5`` and ``\tau = 0.01`` ns
``= 10^{-11}`` s,

```json
{"Type": "Debye", "DeltaPermittivity": 0.5, "RelaxationTime": 0.01}
```

is the same term as

```json
{"Type": "PoleResidue", "Pole": -1.0e11, "Residue": 5.0e10}
```

The two forms give the same ``\epsilon(\omega)`` and the same field solve, but the
`PoleResidue` form specifies no stored energy (see below): `E_elec` and the outputs
normalized by it are omitted, and the error indicators are normalized differently, so that
adaptive mesh refinement can differ.

A pole and residue fit in GHz must be multiplied by ``2\pi \times 10^9`` first: the check for
gain does not depend on the scale of the rates, and does not detect a fit in other units.

Each term is the polarization of a model driven by the field: Debye is a relaxation,
`DjordjevicSarkar` a continuum of them, and Lorentz and Drude are oscillators. With the
parameter bounds above, these typed terms store nonnegative energy and dissipate
nonnegative power at every frequency. A `PoleResidue` term may describe gain, for instance
a fit that is only accurate over the frequencies of interest. It is allowed, with a
warning identifying the detected frequency intervals where the complete model, all its
terms together, has gain (a negative dissipated power), and whether they overlap the
requested band. The field solve uses the model as given. This interval search is sampled,
not an all-frequency passivity certificate.

A `PoleResidue` term specifies a susceptibility, not a model of the polarization, and so
not its stored energy: the same susceptibility has realizations with different stored
energies, and a termwise energy would depend on how the response is split into poles. For
example, ``1/(s+1) + 1/(s+2)`` stores ``0.9`` at ``\omega = 1`` as two typed Debye terms,
and other energies in other realizations. A material with a `PoleResidue` term therefore has
no stored electric energy, even for a single Debye, Drude, or Lorentz section written as
`PoleResidue` terms, which has the response of the typed term. An exact zero pole ``r/s``
with real ``r`` has zero additional stored energy, by convention, and does not remove the
stored energy of the material: its harmonic current is non-reactive, and its signed power
is dissipated (``r > 0``) or supplied (``r < 0``). It is assembled in conductivity form,
but it remains configured polarization in ``\epsilon`` and ``\mathbf{D}``; a user
`Conductivity` is separate and is added once. Use typed terms where the stored energy is
needed. *Palace* warns about each material with no stored energy when it reads the
configuration; this is not a statement that the material has gain, and it does not change
the field solve.

The object form requires scalar permittivity and cannot be combined with `LossTan`;
ordinary `Conductivity` remains additive. The high-frequency value uses the ordinary
permittivity mass matrix. Exact zero poles use the ordinary conductivity coefficient, and
all nonlinear terms on a material share one unit vector-mass operator.

Transient problems are not supported with this material model. Electrostatic problems use
the static permittivity ``\epsilon(0)``, so they reject a pole at zero frequency (a `Drude`
term or a `PoleResidue` term with a zero pole), and magnetostatic
problems do not use the permittivity. Numeric wave ports, boundary-mode problems, Floquet
ports, and Absorbing boundaries next to the material use ``\epsilon(i\omega)`` at each
frequency. The eigenmode divergence-free projection is disabled: its fixed
high-frequency weighting does not impose the frequency-dependent interface condition, and
longitudinal modes at permittivity zeros can be physical. With the default solver
settings, eigenmode problems whose only frequency-dependent terms are the pole terms of the
materials are solved as linear eigenvalue problems, including the modes that accumulate near the resonance of a partially
filled region. With other frequency-dependent terms the nonlinear eigensolver is used,
which does not guarantee the modes nearest the target or all the requested modes; see
[Eigenmode problems](problem.md#Eigenmode-problems).
In adaptive circuit synthesis, the pole terms of the volume permittivity are realized
exactly by auxiliary states, in the reduced-basis directions that the synthesis rank
cutoff retains. `DjordjevicSarkar` terms are fit over the frequency sweep, and so are
wave ports and `Absorbing` boundaries next to the material: the complete complex
wave-port propagation constant, and the absorbing-boundary scalars of the wave
permittivity (``\epsilon(i\omega)`` with the material `Conductivity` included). These
responses are in general not rational, and near a material resonance the fits may not
meet `AdaptiveTol`. *Palace* prints each fit residual and warns when it exceeds the
effective synthesis fit tolerance. Only the rational volume terms are realized exactly;
the material-dependent port and absorbing-boundary responses and native
`DjordjevicSarkar` response are fits. Their real-band fit control does not certify
complex-plane or resonance accuracy.

### Energy and postprocessing

The domain electric energy and the ``U_e`` visualization field include the stored energy of
the field and polarization, including an oscillator's kinetic energy. Their material
weight is the model's stored-energy weight ``W(\omega)``, not generally
``\operatorname{Re}\epsilon(i\omega)``. Eigenmode postprocessing evaluates material
properties at the real part of the mode frequency.

The participation ratios of lumped inductive ports and of interface dielectrics use the
electric-side normalization
``U = (E_\mathrm{stored} + E_{\operatorname{Re}\epsilon})/2 + E_\mathrm{cap}``, the
time-averaged electric energy of the mode, and the mode coupling rates its total stored
energy ``2U``. The domain participations `p_elec[k]` in `domain-E.csv` divide the stored
electric energy of each domain by the total stored electric energy instead.
Here ``E_{\operatorname{Re}\epsilon}`` is the field energy weighted by the real part of
the permittivity; it equals the magnetic field energy for a lossless mode. This avoids the
slower convergence of magnetic-field reconstruction and reduces to
``E_\mathrm{elec} + E_\mathrm{cap}`` for constant materials. For lossy eigenmodes it is
an approximation evaluated at the real frequency, not an exact complex-frequency energy
identity.

For a material on the mesh with no stored electric energy (see above), *Palace* does not
report a value in its place. It omits, for the whole simulation and
with a warning that lists the actual omissions: the total
electric energy ``E_\mathrm{elec}`` and the electric participations in `domain-E.csv`, the
electric energy of each postprocessing domain that contains the material, the ``U_e``
field, and the outputs normalized by ``U``: `surface-Q.csv`, `port-Q.csv`, and
`port-EPR.csv`. The field solution, S-parameters, port voltages, currents, and powers,
fluxes, boundary charges, eigenvalues and their quality factors, magnetic energies, and the
electric energy of postprocessing domains without the material are reported as usual.

Electric `SurfaceFlux` and boundary charge ``Q_s`` use the complete complex
``\epsilon(i\omega)``, including zero-pole polarization. Ordinary `Conductivity` remains
free current and is not included in displacement. Metal-substrate interface participation
uses ``|\epsilon(i\omega)|^2`` for a dispersive substrate's normal displacement; constant
substrates retain the existing real-permittivity approximation. Point probes still sample
``\mathbf{E}`` and ``\mathbf{B}``; there is no separate displacement or energy-density
probe output.

The Stratton–Chu far field retains its homogeneous, isotropic, lossless-exterior
approximation. For a dispersive exterior its wave speed uses
``\mu\operatorname{Re}\epsilon(i\operatorname{Re}\omega)`` and requires positive real
permittivity next to the observation surfaces. This does not restrict unrelated interior
materials. The wave number retains the mode's complex frequency.

### Error indicators

With frequency-dependent materials, the electric flux estimator recovers the complex
physical displacement ``\mathbf{D}=\epsilon(i\omega)\mathbf{E}`` and weights its residual by
the inverse stored-energy tensor ``W^{-1}(\omega)``. The recovery is the projection weighted
by ``\epsilon_\infty^{-1}``, the same at every frequency. It is therefore not the exact
minimizer of the residual at the interfaces of a frequency-dependent material, where
``W(\omega)`` and ``\epsilon_\infty`` differ, most of all near a resonance. For a material
whose stored electric energy *Palace* cannot report, the indicator uses ``\epsilon_\infty``
in place of ``W`` and is normalized by the energy norm of the recovered fluxes in place of
the total field energy, with a warning. These are numerical norms, which keep adaptive
refinement usable; they are not physical energies and do not validate the physics of a fit
outside its useful band.
