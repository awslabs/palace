<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# General graded-mesh regression results — 2026-09-10

## Acceptance policy

0.1% interface agreement is retained as an optional strict diagnostic, not a
universal requirement. The provisional engineering profile uses 0.05% domain
energy and 0.5% fabricated interface-energy agreement. Neither profile grants
library qualification. Reference uncertainty, basis/operator checks, defect
cancellation, and final-device observables remain separate requirements.

## Geometry/mesh suite

All **19** bounded scouts passed in
`/tmp/general-coupon-mesh/regression-final/summary.json`. The directory contains
hashed snapshots of the tools and audit executables used for the run.

| Geometry | Thin tets | Fabricated tets |
|---|---:|---:|
| Closed strip | 5,998 | 8,285 |
| Concave multislot mask | 9,646 | 15,367 |
| Rotated/translated concave mask | 10,673 | 24,313 |
| Two conductors | 8,642 | 13,460 |
| Mask with a hole | 13,388 | 19,904 |
| Opposed layers, repeated conductor ID | 16,185 | 39,986 |
| T junction | 8,606 | 15,307 |
| Three conductors | 10,374 | 17,219 |
| Rounded strip | — | 9,421 |
| 80-degree sidewalls | 6,639 | 8,687 |

These are coarse 20 nm scouts, **not production response meshes**. They exercise
arbitrary in-plane orientation, concavity, holes, different edge counts, conductor
counts, slots, layer signs, thicknesses and etch depths. No case-specific mesher
branch is used.

Checks passed:

- independent analytic metal-interface area gates (including the sloped strip);
- expected physical label sets and no duplicate physical face assignments;
- coverage of every exterior/material-interface face;
- positive quadrature-point Jacobians;
- quadrature-integrated mesh/CAD volume and physical-family area agreement.

Maximum kappa across these scouts was about 17.03. Planar area differences were
at floating-point scale. The genuinely rounded, quadratic mesh differed from CAD
in surface area by about **0.018%**, below the explicit 0.1% curved scout target.
The rounded-case gate now checks that rounding actually changed the geometry:
an earlier exploratory run exposed a silently skipped fillet, which was fixed.

The hole test also led to a material-side offset fix and handling of collapsed
convex holes. An independent trench-volume unit test verifies etching inside the
hole. A final legacy three-edge geometry-only check still matches the retained
CAD measures to about **3.0e-11 relative**.

### Remaining partition qualification

Coarse multislot scouts have maximum sampled ambiguous-area fractions of roughly
10–24%. Every triangle is exclusively assigned and all expected regions exist,
but this does not establish an accurate slot boundary or an interface-energy
error bound. Energy-weighted partition convergence remains required before using
these scouts to regenerate a production library.

## Numerical covariance check

A two-conductor fabricated scout was solved at p5, then rigidly rotated by 0.63
radians and translated without changing connectivity or physical labels. Both
independent conductor excitations were solved on two local MPI ranks, each under
60 s / 6 GiB guards.

- Each solve took about 18.2 s, peak process-tree RSS about 2.65 GB.
- Domain energies were identical at CSV precision.
- The complete 2-by-2 terminal capacitance matrix was byte-identical.
- Worst interface-energy relative difference: **4.50e-10**.

Artifacts: `/tmp/general-coupon-mesh/rotation-probe2`.
This verifies coordinate covariance of a multi-conductor solve on the same FE
mesh, not mesh convergence or the full spatial trace response operator.

## Automated checks

- Geometry/grading/ownership/offset tests: 56 assertions.
- Probe comparison tests: 5 tests, including engineering versus strict profiles,
  thin-SPR exclusion, and ambiguous matrix-schema rejection.
- Integrated mesh-measure tests: 2 tests, including missing-boundary rejection.
- Portable audit-tool makefile compiled successfully against the local MFEM build.

All existing reference libraries and production meshes remain unchanged.
See `graded-mesh-qualification.md` for supported cases and explicit limitations.
