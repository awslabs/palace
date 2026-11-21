```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Plane Wave Excitation and Radar Cross Section

Scattering problems, such as the radar cross section of an object, are driven by an
incident plane wave that comes from infinity. *Palace* excites such a wave inside a
bounded domain with the total-field/scattered-field (TF/SF) formulation
[[1]](#References): a closed surface ``S`` encloses the scatterer, the unknown is the
total field ``\mathbf{E}`` inside ``S`` and the scattered field
``\mathbf{E}_s = \mathbf{E} - \mathbf{E}_\text{inc}`` outside, and the incident field
enters only on ``S``. The truncation of the domain (an absorbing boundary or a perfectly
matched layer) therefore only sees the outgoing scattered field. From the solution,
*Palace* reconstructs the total and the scattered fields everywhere, for the outputs and
the far field.

## Formulation

The surface ``S`` is the boundary of a total-field region ``T``, given by volume
attributes. The unknowns are unchanged: on the degrees of freedom of the closure of ``T``
(those on ``S`` included) they are the total field, and on all the other degrees of
freedom ``O`` they are the scattered field. With the system matrix ``A`` of the total
field, which is not modified, and the interpolant ``\mathbf{g}`` of the incident field,
the right-hand side is

```math
\mathbf{b}_T = -(A\,\mathbf{g}_O)_T, \qquad \mathbf{b}_O = (A\,\mathbf{g}_T)_O,
```

where ``\mathbf{g}_T`` and ``\mathbf{g}_O`` are the restrictions of ``\mathbf{g}`` to
the degrees of freedom of ``T`` and ``O``. Only the elements outside of ``T`` which touch
``S`` couple the two sets of degrees of freedom, so the incident field is only evaluated
in that layer of elements: this is the discrete counterpart of the equivalent (Love)
electric and magnetic currents ``\mathbf{J} = \mathbf{n} \times \mathbf{H}_\text{inc}``
and ``\mathbf{M} = -\mathbf{n} \times \mathbf{E}_\text{inc}`` on ``S``, which radiate the
incident field inside ``S`` and nothing outside. The incident field is known in closed
form only in a homogeneous medium, so these elements must have one isotropic, lossless
material without a frequency-dependent permittivity, the background, which *Palace*
checks at setup; inside ``S``, the materials and boundary conditions are arbitrary
(conductors, lossy, dispersive, magnetic or anisotropic materials).

The incident field is

```math
\mathbf{E}_\text{inc}(\mathbf{x}) = E_0\, \hat{\mathbf{p}}\,
e^{-i k \hat{\mathbf{k}} \cdot (\mathbf{x} - \mathbf{x}_0)}
```

with the propagation direction ``\hat{\mathbf{k}}``, the polarization ``\hat{\mathbf{p}}``
(its component along ``\hat{\mathbf{k}}`` is removed), the amplitude ``E_0``, the
wavenumber ``k = \omega \sqrt{\varepsilon_r \mu_r} / c_0`` of the background and the point
``\mathbf{x}_0`` of zero phase.

## Setup

A plane wave is an excitation of
[`config["Domains"]["PlaneWave"]`](../config/reference.md#config-domains):

```json
"PlaneWave":
[
    {
        "Index": 1,
        "Attributes": [1, 2],    // Volume attributes of the total-field region T
        "Direction": "+Z",       // Propagation direction
        "Polarization": "+X",    // Polarization of the electric field
        "Amplitude": 1.0,        // V/m
        "Origin": [0.0, 0.0, 0.0]  // Point of zero phase
    }
]
```

`"Index"` is the excitation index: several plane waves with different indices (for example
several directions of incidence) are solved one after the other, as port excitations are.
The mesh needs:

  - a total-field region ``T`` containing the scatterer, separated from the truncation of
    the domain by at least one layer of elements of the background;
  - for the far field and the radar cross section, a closed surface
    ([`config["Boundaries"]["Postprocessing"]["FarField"]`](../config/reference.md#config-boundaries))
    around the scatterer, in the background material. It is usually an interior surface
    between ``S`` and the truncation, but, since the far field is computed from the
    scattered field everywhere, it may also lie inside ``T`` or cross ``S``.

Plane wave excitations are supported for driven simulations, with a uniform or an
[adaptive](adaptive_driven_solver.md) frequency sweep. The excitation
depends on the frequency through the phase of the incident field, which is not linear in
the frequency: the reduced-order model of the adaptive sweep takes its exact projection at
each frequency, from the incident field interpolated on the layer of elements around
``S`` only. This requires that the frequency-dependent terms of the system matrix (wave
and Floquet ports, `"Conductivity"` and rational `"Impedance"` boundaries, second-order
`"Absorbing"` boundaries) do not couple ``T`` and the rest of the domain, that is, that
these boundaries do not cross ``S``, which *Palace* checks at setup. The outputs of plane
waves are computed from the full solution at every frequency of the online sweep, and
circuit synthesis (`"AdaptiveCircuitSynthesis"`) does not support them.

On a periodic cell, the incident field must be quasi-periodic: with a
[Floquet wave vector](../config/reference.md#config-boundaries-periodic-floquetwavevector)
``\bm{k}_F``, the unknown is the periodic envelope ``\bm{E}\,e^{i\bm{k}_F\cdot\bm{x}}``
and *Palace* interpolates the envelope of the incident field, and ``\bm{k}_F`` (zero if
not given) must be the part of the wave vector ``k\hat{\bm{k}}`` of the plane wave along
the periodic directions, up to a reciprocal lattice vector. Oblique incidence therefore
always needs a Floquet wave vector. *Palace* checks this at setup for all the frequencies
of the sweep, against the translation of each periodic boundary pair: it stops with an
error for a phase mismatch above ``10^{-2}`` rad over a lattice translation, and warns
above ``10^{-6}`` rad, since the source then leaks into the scattered field at the level of
the mismatch. Periodic boundary pairs which are not translations are rejected, and for
pairs without an explicit `"Translation"` (detected from the mesh) the check is done on
the interpolated incident field at each frequency instead. A fixed ``\bm{k}_F`` matches an
oblique plane wave only at isolated frequencies; with
[`config["Boundaries"]["Periodic"]["FloquetReferenceFrequency"]`](../config/reference.md#config-boundaries-periodic)
it scales with the frequency, as the wave vector of the plane wave does, and must then
match it exactly, so a frequency sweep at a fixed angle of incidence stays consistent.
Floquet port S-parameters (`port-floquet-S.csv`) are normalized to an incident Floquet
mode, and are not written for excitations with a plane wave.

## Accuracy

The discrete incident field is not exactly a solution of the discrete problem, so part of
it leaks out of ``S``: without a scatterer, the scattered field outside of ``S``, which
should vanish, measures this leakage. It comes from the dispersion error of the discrete
plane wave across the total-field region, and decreases quickly with the element size and
the order. A simulation of the empty total-field region (same mesh without the scatterer)
is the simplest check of the setup: the leakage should be well below the scattered field
of interest. This matters most for weak scattered fields, such as the backscattered field
of a scatterer whose forward scattering is much stronger.

## Output

The discrete solution ``\mathbf{x}`` is the total field on the degrees of freedom of the
closure of ``T``, where ``M = 1``, and the scattered field elsewhere. *Palace* reconstructs
the total field ``\mathbf{E} = \mathbf{x} + (1 - M)\,\mathbf{g}`` and the scattered
field ``\mathbf{E}_s = \mathbf{x} - M\,\mathbf{g}`` everywhere, and the magnetic flux
densities from them. All the outputs (fields, energies, ports, probes, surface fluxes) use
the total field, as for any other excitation; in the ParaView and grid function outputs,
the scattered fields are also written as `Es_real`, `Es_imag`, `Bs_real` and `Bs_imag`.
The error estimates take the total field in ``T``, where the mesh must resolve it, and the
scattered field outside, where only the scattered field propagates, both normalized by the
energy of the same fields; adaptive mesh refinement therefore does not refine the
scattered-field region to resolve the incident wave.

The far field in `farfield-rE.csv` is the far field of the scattered field. When all the
excitations are plane waves, *Palace* also writes the bistatic radar cross section in
`farfield-rcs.csv` and the monostatic radar cross section in `rcs-monostatic.csv`; see
[Radar cross section](farfield.md#Radar-cross-section) for their columns. The
[coated sphere example](../examples/rcs_coated_sphere.md) shows a complete setup and its
validation against the Mie series.

## References

[1] A. Taflove and S. C. Hagness, *Computational Electrodynamics: The Finite-Difference
Time-Domain Method*, 3rd ed., Artech House, 2005, Chapter 5.
