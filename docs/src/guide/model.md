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
from driven problems in the time domain.

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
wishing to model both superconducting and normal currents.

### Frequency-dependent permittivity

Scalar additive permittivity models use the object form

```json
"Permittivity": {
  "HighFrequency": 2.08,
  "Terms": [
    {"Type": "Drude", "PlasmaFrequency": 1.0, "CollisionFrequency": 0.1}
  ]
}
```

`HighFrequency` is the scalar ``\epsilon_\infty`` and `Terms` is a nonempty additive array,
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
    (rad/s). Values are real or `[real, imag]`; upper-half-plane poles include their
    conjugates automatically.
  - `DjordjevicSarkar`: positive `Strength`, `LowerFrequency`, and
    `UpperFrequency` bounds in GHz, with the upper bound larger than the lower. This native
    term is ``\mathrm{Strength}\,\log((s+\omega_\mathrm{hi})/(s+\omega_\mathrm{lo}))``, so a
    static increment ``\Delta\epsilon`` corresponds to
    ``\mathrm{Strength} = \Delta\epsilon/\log(\omega_\mathrm{hi}/\omega_\mathrm{lo})``.

Each term is the polarization of a model driven by the field: Debye is a relaxation,
`DjordjevicSarkar` a continuum of them, and Lorentz and Drude are oscillators. With the
parameter bounds above, these typed terms store nonnegative energy and dissipate
nonnegative power at every frequency. A `PoleResidue` term may describe gain or negative
stored energy, for instance a fit that is only accurate over the frequencies of interest.
It is allowed, with warnings identifying detected frequency intervals where the complete
model has gain or negative stored energy, and whether they overlap the requested band.
This interval search is sampled, not an all-frequency passivity certificate.

The stored energy of a `PoleResidue` term is that of a relaxation for a real pole, and of
an oscillator ``\ddot{x} + a_1\dot{x} + a_0 x = E`` with polarization
``b_0 x + b_1 \dot{x}`` for a conjugate pair. When ``b_1 \neq 0``, the constitutive
response does not uniquely determine microscopic stored energy; *Palace* uses its canonical
auxiliary-state midpoint convention.

The object form requires scalar permittivity and cannot be combined with `LossTan`;
ordinary `Conductivity` remains additive. The high-frequency value uses the ordinary
permittivity mass matrix. Exact zero poles use the ordinary conductivity coefficient, and
all nonlinear terms on a material share one unit vector-mass operator.

Transient problems are not supported with this material model. Electrostatic problems use
the static permittivity ``\epsilon(0)``, so they reject a `Drude` term, and magnetostatic
problems do not use the permittivity. Numeric wave ports, boundary-mode problems, Floquet
ports, and Absorbing boundaries next to the material use ``\epsilon(i\omega)`` at each
frequency. The eigenmode divergence-free projection is disabled because
longitudinal plasma modes can be physical. Eigenmodes of a region partially filled with a
pole material accumulate at the material resonance, so choose the eigenmode `Target` away
from material resonances; the eigensolver warns when it cannot resolve the requested modes.
In adaptive circuit synthesis, pole terms are realized exactly, while `DjordjevicSarkar`
terms are fit over the frequency sweep.

### Energy and postprocessing

The domain electric energy and the ``U_e`` visualization field include the stored energy of
the field and polarization, including an oscillator's kinetic energy. Their material
weight is the model's stored-energy weight ``W(\omega)``, not generally
``\operatorname{Re}\epsilon(i\omega)``. Eigenmode postprocessing evaluates material
properties at the real part of the mode frequency.

Mode participation ratios and coupling rates use the electric-side normalization
``U = (E_\mathrm{stored} + E_{\operatorname{Re}\epsilon})/2 + E_\mathrm{cap}``.
Here ``E_{\operatorname{Re}\epsilon}`` is the field energy weighted by the real part of
the permittivity; it equals the magnetic field energy for a lossless mode. This avoids the
slower convergence of magnetic-field reconstruction and reduces to
``E_\mathrm{elec} + E_\mathrm{cap}`` for constant materials. For lossy eigenmodes it is
an approximation evaluated at the real frequency, not an exact complex-frequency energy
identity.

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
physical displacement ``\mathbf{D}=\epsilon(i\omega)\mathbf{E}`` and weights its
residual by the inverse stored-energy tensor ``W^{-1}(\omega)``. If a general model has
nonpositive stored-energy weight at the requested frequency, the indicator alone falls
back to ``\epsilon_\infty`` with a warning; signed reported energies are unchanged. A
nonpositive total-energy normalizer similarly uses a positive recovered-flux norm with a
warning. These fallbacks keep unusual fitted models usable, but do not validate their
physics outside the fit's useful band.
