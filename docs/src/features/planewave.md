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
[[1]](#References): a closed surface ``S`` encloses the scatterer, the solution is the
total field ``\mathbf{E}`` inside ``S`` and the scattered field
``\mathbf{E}_s = \mathbf{E} - \mathbf{E}_\text{inc}`` outside, and the incident field
enters only on ``S``. The truncation of the domain (an absorbing boundary or a perfectly
matched layer) therefore only sees the outgoing scattered field, and the far field
computed outside of ``S`` is the far field of the scattered field.

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
form only in a homogeneous medium, so these elements must be vacuum; inside ``S``, the
materials and boundary conditions are arbitrary (conductors, lossy, magnetic or
anisotropic materials).

The incident field is

```math
\mathbf{E}_\text{inc}(\mathbf{x}) = E_0\, \hat{\mathbf{p}}\,
e^{-i k \hat{\mathbf{k}} \cdot (\mathbf{x} - \mathbf{x}_0)}
```

with the propagation direction ``\hat{\mathbf{k}}``, the polarization ``\hat{\mathbf{p}}``
(its component along ``\hat{\mathbf{k}}`` is removed), the amplitude ``E_0``, the
wavenumber ``k`` of vacuum and the point ``\mathbf{x}_0`` of zero phase.

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
    the domain by at least one layer of vacuum elements;
  - for the far field and the radar cross section, an interior closed surface
    ([`config["Boundaries"]["Postprocessing"]["FarField"]`](../config/reference.md#config-boundaries)
    with a `"Center"`) in the scattered-field region, outside of the elements which touch
    ``S``. *Palace* stops with an error if the far-field surface touches the total-field
    region, where the far field would include the incident wave.

Plane wave excitations are supported for driven simulations with a uniform frequency sweep.

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

The far field in `farfield-rE.csv` is the far field of the scattered field. When all the
excitations are plane waves, *Palace* also writes the bistatic radar cross section in
`farfield-rcs.csv` and the monostatic radar cross section in `rcs-monostatic.csv`; see
[Radar cross section](farfield.md#Radar-cross-section) for their columns. The
[coated sphere example](../examples/rcs_coated_sphere.md) shows a complete setup and its
validation against the Mie series.

## References

[1] A. Taflove and S. C. Hagness, *Computational Electrodynamics: The Finite-Difference
Time-Domain Method*, 3rd ed., Artech House, 2005, Chapter 5.
