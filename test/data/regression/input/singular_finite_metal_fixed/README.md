<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Fixed-wedge finite-metal regression

This small cavity contains a one-millimeter metal cube inside a four-millimeter PEC box.
The complete six-face metal boundary is selected for fixed-wedge extraction. Its edge
skeleton contains 12 physical mesh segments and eight degree-three junctions.

The case is intentionally small enough to validate the paper-compatible approximation locally:

| Configuration | Standard order `p` | Singular order `s` | Setup time | First two modes (GHz) |
|---|---:|---:|---:|---|
| `singular_finite_metal_fixed_regression.json` | 2 | 1 | 0.204 s | 47.84454, 47.86461 |
| `singular_finite_metal_fixed_s2.json` | 2 | 2 | 0.767 s | 47.90210, 47.92214 |
| `singular_finite_metal_fixed_s3.json` | 2 | 3 | 2.876 s | 47.94771, 47.95786 |

The singular cases use `AbsTol = 5e-3`, `RelTol = 1e-3`, and ARPACK. Tight entrywise
quadrature tolerances such as `2e-6` are intentionally not used in this smoke/convergence
fixture: adaptive integration of every mixed high-order singular matrix entry becomes much
more expensive without changing feature topology. Separate quadrature-tolerance studies
should use a targeted element-level test rather than a large device solve.

The fixed-wedge solve is lossless. Interface and other dissipative physics must be evaluated
in postprocessing until low-loss complex-material support is implemented for this model.
