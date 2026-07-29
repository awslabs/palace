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

# 2D Coplanar Waveguide Mode Analysis

## Problem description

This example demonstrates the [`"BoundaryMode"`](../config/reference.md#config-solver-boundarymode)
simulation type by computing propagation
constants, effective indices, and characteristic impedance for a coplanar waveguide (CPW)
cross-section. It uses a similar CPW geometry as the
[3D crosstalk example](cpw.md) but solves the 2D eigenvalue problem directly on the
cross-section mesh, making it much cheaper than a full 3D driven simulation.

Two baseline configurations are provided in the
[`examples/cpw2d/`](https://github.com/awslabs/palace/blob/main/examples/cpw2d) directory:

  - **Thin metal** (`cpw2d_thin.json`): Zero-thickness PEC traces on a dielectric substrate.
  - **Thick metal with impedance BC** (`cpw2d_thick_impedance.json`): Finite-thickness
    traces with a surface impedance boundary condition modeling kinetic inductance.

The CPW has a center trace on a silicon substrate (``\varepsilon_r = 11.47``,
``\tan\delta = 1.2 \times 10^{-7}``). The mesh length unit is ``\mu\text{m}``.

## Configuration

Both configurations use
[`"Problem": {"Type": "BoundaryMode"}`](../config/reference.md#config-problem-type) and request 2
modes at 5 GHz. The key
[solver settings](../config/reference.md#config-solver-boundarymode) are:

```json
"BoundaryMode":
{
    "Freq": 5.0,
    "N": 2,
    "Save": 2,
    "Target": 2.497,
    "Tol": 1.0e-8
}
```

The `"Target"` parameter sets the effective index ``n_\text{eff} = k_n / k_0`` near which
the eigenvalue solver searches via shift-and-invert, where ``k_n`` is the propagation
constant and ``k_0 = \omega / c_0`` is the free-space wavenumber. The solver finds modes
with ``n_\text{eff}`` closest to (but not necessarily above) the target. A value near the
expected ``n_\text{eff}`` (between 1 for air and ``\sqrt{\varepsilon_r} \approx 3.39`` for
the substrate) helps the solver converge to the desired propagating modes.

### Thin metal configuration

The thin metal case (`cpw2d_thin.json`) uses PEC boundary conditions on the trace and ground
boundaries. It also specifies impedance and voltage postprocessing with a coordinate path
across the CPW gap:

```json
"Postprocessing":
{
    "Impedance":
    [
        {
            "Index": 1,
            "VoltagePath": [[518.5, 0.0], [522.0, 0.0]],
            "NSamples": 200
        }
    ],
    "Voltage":
    [
        {
            "Index": 1,
            "VoltagePath": [[518.5, 0.0], [522.0, 0.0]],
            "NSamples": 200
        }
    ]
}
```

Additionally, domain energy postprocessing and interface dielectric loss (SA, MS, MA types)
are configured for energy participation ratio analysis.

### Experimental singular-element configurations

Five configurations exercise the additive singular-element implementation on the
quadratic Gmsh mesh used by the baseline thin-metal example:

  - `cpw2d_thin_electrostatic_standard.json`
  - `cpw2d_thin_electrostatic.json`
  - `cpw2d_thin_boundarymode_standard.json`
  - `cpw2d_thin_boundarymode_singular.json`
  - `cpw2d_thin_driven_singular.json`

The finite-thickness PEC comparison uses
`cpw2d_thick_boundarymode_standard.json` and
`cpw2d_thick_boundarymode_singular.json` on the same quadratic thick-metal mesh
with standard polynomial order two. The singular configuration selects both PEC
faces at each metal corner. Palace walks the complete incident triangle fan,
groups its ordered isotropic material sectors, and computes the smallest positive
Dirichlet transmission-wedge exponent. For the substrate-vacuum corners in this
mesh, the exponent is approximately `nu = 0.5255535`; a homogeneous 270-degree
dielectric corner instead gives `nu = 2/3`.

The two standard configurations and their singular counterparts use the same mesh and
materials. The singular configurations add first-order triangular basis functions at the
four endpoints of the zero-thickness trace and ground lines:

```json
"SingularElements":
{
    "Attributes": [1, 2],
    "Order": 1,
    "QuadratureOrder": 8,
    "AbsTol": 2.0e-6,
    "RelTol": 2.0e-6,
    "MaxSubdivisions": 6
}
```

Generate the mesh with
`julia --project=examples examples/cpw2d/mesh/mesh_thin.jl`. A high-order nodal triangle
whose map is actually affine uses the same pretabulated reference-integral path as a
linear mesh. A genuinely curved enriched triangle instead uses singularity-aligned
reference quadrature with its physical Jacobian evaluated pointwise. Curved elements
away from singular features are unrestricted. Selected PEC line segments themselves
must remain geometrically straight and regular, but may use a nonuniform high-order
parametrization along that straight image. Conforming and nonconforming AMR may refine
through enriched elements, after which Palace reconstructs the physical singular
features and enrichment degrees of freedom; mesh rebalancing is supported. Parallel
uniform refinement, GPU execution, and triangular singular order greater than one are
not supported.

Selected singular attributes are automatically excluded from Palace's internal-boundary
cracking pass because their conforming two-sided topology defines the singular feature.
The example files also set `Model.CrackInternalBoundaryElements` to `false` for clarity,
but that global setting is not required; nonselected boundary attributes continue to
follow it normally.

The electrostatic solve reports capacitance, combined-field domain energy, canonical
singular coefficients, fitted tip slopes, and combined-field dielectric surface
participation. The singular `BoundaryMode` solve supports isotropic materials, including
bulk dielectric loss, with PEC or natural PMC boundaries. It reports propagation
constants, combined-field energies, canonical ND/H1 coefficients, fitted transverse field
slopes, and dielectric surface participation. Voltage, impedance, and saved grid-function
output remain disabled on the enriched `BoundaryMode` path. Its AMR indicator is computed
from the standard-space smooth remainder of the combined solution.

The thin-sheet examples specify an explicit cutoff on each dielectric interface:

```json
"EdgeCutoff": 2.0e-3
```

The value is in mesh length units, so `2.0e-3` is 2 nm for these micrometer meshes.
The reported value in `surface-Q.csv` is the outer thin-layer integral with the
neighborhood closer than `EdgeCutoff` to each singular line tip excluded. It includes the
complete standard-plus-singular field. Decreasing the cutoff increases the result
logarithmically. This is a cutoff-dependent diagnostic, not a finite-thickness response
correction, and it does not include the energy in the excluded neighborhood. Omitting the
cutoff for an active `nu <= 1/2` trace is rejected rather than returning a quadrature- or
mesh-dependent finite number.

The finite-thickness thick-metal corner has `nu > 1/2`, so its surface trace is integrable.
The paired thick `BoundaryMode` configurations therefore leave `EdgeCutoff` at zero and
the singular case computes the untruncated combined-field MS and MA participation.

The driven configuration is a restricted full-wave validation case with a surface-current
source. Singular driven and eigenmode simulations currently require CPU execution and
isotropic materials without bulk conductivity or London penetration depth. Conforming and
nonconforming AMR, including refinement through enriched elements and subsequent mesh
rebalancing, is supported. Bulk dielectric loss, lumped-port `R/L/C` terms, and
combined-field dielectric surface postprocessing are supported.

Reported electric and magnetic field energies exclude reactive lumped-boundary matrix
terms, matching standard Palace domain postprocessing. Eigenmode balance still uses the
complete augmented matrices, while surface participation and lumped EPR use bulk electric
plus measured capacitor energy for normalization.

For polynomial order greater than one, the outer Krylov method remains the configured
solver and may use polynomial multigrid. The transfer at every level is
`diag(P_standard, I_enrichment)`, the paired H1-to-ND gradient is retained at every level,
and the coarse operator contains the complete standard-plus-singular coupling. Sparse
direct solver types such as SuperLU_DIST or STRUMPACK therefore factor only the combined
`p = 1` coarse matrix, not the finest system. AMS uses AMS on the standard ND principal
block and an enrichment correction; this path is validated for driven solves but is not
yet certified for eigenmode shift-and-invert.

Parallel uniform refinement, wave and Floquet ports, periodic boundaries, impedance or
absorbing boundaries, saved standard-grid fields, general field or probe postprocessing,
adaptive frequency sweeps, and PROM are not yet supported on the singular full-wave path.

The zero-thickness surface integral of ``|E|^2`` still diverges logarithmically at a sheet
edge. A fabricated-device participation requires resolved finite-thickness interface
geometry or a separately validated response correction which combines the converged outer
field with an inner edge model.

### Thick metal with impedance BC

The thick metal case (`cpw2d_thick_impedance.json`) replaces PEC boundaries with a surface
impedance boundary condition that models the kinetic inductance of a superconducting film:

```json
"Impedance":
[
    {
        "Attributes": [1],
        "Ls": 1.332e-13
    }
]
```

The surface inductance ``L_s`` adds an inductive contribution to the boundary condition,
which increases the effective index of the guided modes compared to the ideal PEC case.

The figure below shows the ``E_x`` component of the first propagating mode for the thick
metal configuration, zoomed into the trace and gap region. The field is concentrated in the
gaps between the center conductor and ground planes, with singularities at the metal corners.

```@raw html
<!-- Plot generated with examples/cpw2d/plot_cpw2d_docs.py -->
<br/><p align="center">
  <img src="../../assets/examples/cpw2d-1.png" width="70%" />
</p><br/>
```

## Results

The mode analysis solver writes propagation constants and effective indices to `mode-kn.csv`
and characteristic impedance to `mode-Z.csv`.

### Effective index

For the thin metal (PEC) case:

```@example include_example
include_example_file("cpw2d/thin", "mode-kn.csv") # hide
```

For the thick metal with impedance BC:

```@example include_example
include_example_file("cpw2d/thick_impedance", "mode-kn.csv") # hide
```

The impedance BC shifts the effective index upward for both modes. This is expected: the
surface inductance ``L_s`` increases the effective path length seen by the wave, raising
``n_\text{eff}``. The shift is larger for mode 1, which has more field energy concentrated
near the conductor surfaces.

### Characteristic impedance

The power-voltage characteristic impedance ``Z_\text{PV} = |V|^2 / (2P)`` is computed from
the voltage line integral across the CPW gap and the mode power. For the thin metal case:

```@example include_example
include_example_file("cpw2d/thin", "mode-Z.csv") # hide
```

For the thick metal with impedance BC:

```@example include_example
include_example_file("cpw2d/thick_impedance", "mode-Z.csv") # hide
```

The two modes correspond to the even and odd CPW modes, which have different impedance
values. Note that the mode ordering (by propagation constant) can differ between the thin
and thick configurations.
