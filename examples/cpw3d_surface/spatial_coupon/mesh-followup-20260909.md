<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Mesh follow-up: finer fabricated surface resolution

## Candidate retained for p5 validation

`three-fabricated-t6-del2d.msh` uses the existing exact 3D distance grading,
6 nm horizontal edge spacing, 2 nm size-field minimum, growth 1, 500 nm maximum,
and matching-trace height constraints. Surface algorithm: Delaunay (5); volume
algorithm: HXT (10), quality target 0.1, followed by Netgen optimization.

The geometric mesh is linear because all CAD faces/edges in this experiment are
planar/straight (90-degree walls, zero rounding). This is exact geometric
representation, **not a reduction of the p5 solution order**.

- Local meshing: 110.99 s, 3,250,864,128 bytes peak process-tree RSS.
- 1,110,638 tetrahedra, 270,911 geometric vertices.
- Maximum kappa 17.17; median 3.47; p99 7.03.
- Median first surface-element altitudes: SA 5.68 nm, MS 6.13 nm, MA 5.73 nm.
  These are actual corner-geometry altitudes, not h/p. A 2 nm normal layer is
  still NOT enforced; most incident length does not meet a 2 nm altitude gate.
- Physical CAD area/volume gate passed against the reference.
- Candidate and results are separate from all reference libraries.

This is a mesh-resolution/quality improvement, not an accuracy certificate.
The dedicated-host p5 worker/reducer comparison uses original material data,
excitation, and solve tolerance. See
`/data/home/simlap/coupon_mesh_followup_20260909`.

P5 held-out result: 24,860,845 H1 DOFs (versus about 147.9 million in the reference),
7 PCG iterations. On 32 dedicated-host ranks: worker 117.95 s, reducer 80.59 s,
peak process-tree RSS 69,210,783,744 bytes. Differences against retained reference:

| Quantity | Signed difference |
|---|---:|
| Domain energy | -0.011702% |
| MA energy | +3.633553% |
| MS energy | +0.117823% |
| SA energy | +2.625582% |

The configured comparison gate failed. This mesh is **not installed as a replacement**.

## Reference convergence evidence

The retained `7f03270dca8e/p4-p5.csv` describes a **six-probe** comparison, not an
exhaustive proof of full-matrix convergence. Both configurations use the same
`mesh-n2/fabricated.msh`; all six prescribed-potential CSV files are byte-identical
between p4 and p5. Reported fabricated changes are MA matrix norm 2.618957%
(worst energy 5.530516%), MS norm 0.040506% (worst 2.020517%), and SA norm 1.016442%
(worst 1.452175%). These changes do not establish 0.1% reference accuracy.

Thus a growing difference from the reference under tetrahedral refinement cannot
by itself establish which discretization is more accurate. A p6 check on the
retained 10 nm-capped tet mesh completed as an additional convergence diagnostic,
not a change to the p5 production target. It used 30,179,419 H1 DOFs and 11 PCG
iterations; worker/reducer times were 198.75/124.39 s with peak RSS about 82.6 GB.

For this one held-out excitation, the 6 nm-capped p5 mesh and 10 nm-capped p6 mesh
agree in domain energy to 0.000690%, MA to 0.276747%, MS to 0.016956%, and SA to
0.032450%. The MA comparison still fails the 0.1% diagnostic gate. The p6 result
also lies above the old reference in MA (+3.35%) and SA (+2.59%). This is supporting
convergence evidence across two discretization choices, not proof of an exact
continuum limit or full-operator qualification. Thin raw SPR remains excluded.


## Fixed-boundary hybrid-bulk investigation

Added `mesh_discrete_tet_region.jl` and nine smoke-test assertions. The helper checks
closed connected oriented shells, boundary triangle identity, node coordinates,
face multiplicities, positive Jacobians, minimum signed inverse condition, and
volume conservation. It does not claim full-coupon or interface-energy validity.

A far slab of the ten-edge coupon, using the previously aligned/smoothed plan,
was meshed with 198,694 tetrahedra in 7.7 s and roughly 1 GiB. Every fixed boundary
triangle was preserved. However, max kappa was about 409,000, median 42, p99 1,288.
It is unsuitable as a replacement for the reference bulk without further work.

Negative outcomes retained:

- Original plan: HXT inserted a boundary point and replaced two triangles; the
  fixed-boundary check rejected it instead of silently allowing a nonconforming join.
- Aggressive isotropic HXT quality targets repeatedly spent the bounded run in
  optimization of a strongly anisotropic fixed boundary.
- No optimization: near-degenerate slivers with kappa about 1.46 billion.
- Centroid-based cap point layers increased cost and introduced near-singular
  elements; no resulting mesh is approved for a solve.
- Netgen optimization on this discrete-shell construction crashed natively.

Side quads are triangulated by the helper. Their p5 trace space is not identical
to the reference tensor-product trace space. This is documented rather than
claimed to be an exact full-response replacement.

## Anisotropic surface meshing investigation

An opt-in two-stage experiment meshes curves first, uses a BAMG anisotropic
surface metric, then restores scalar 3D distance sizing for the volume mesher.
Surface snapshots and settings are retained. It has NOT produced a qualified
volume mesh: native HXT crashes, Delaunay boundary-recovery failures, and Netgen
crashes occurred. HXT failure was reproduced with Gmsh 4.13.1 on the dedicated
host, so this is not merely a local 4.15-git issue.

The Gmsh BoundaryLayer field was also retested after removing thin-mask vertical
CAD partitions. It still rejects an edge shared by two surfaces. No normal-layer
claim is made from either path. The experimental paths remain opt-in; the retained
6 nm candidate does not use them.

## Validation infrastructure

The surface-altitude auditor now uses a cross product rather than subtracting
nearly equal squared lengths. Three analytic tests check known normal/tangential
lengths, matching-boundary exclusion, and an extreme aspect ratio that previously
rounded the altitude to zero.

The existing geometry/callback tests pass 21 assertions. The fixed-shell tests
pass nine assertions. No physics is disabled to make a mesh or solve pass.

The separately frozen `palace-lazy` binary avoids the unused RT estimator setup
in archive-only/reduce-only modes. An identical-mesh, 32-rank control is checked
before interpreting the new candidate; ordinary raw-only AMR is unchanged.
The 17.7M-DOF control completed in 103.26 s + 77.48 s, versus 200.06 s + 174.87 s
before the change, at the same 32-rank count and on the same mesh. Peak RSS dropped
from approximately 100.2 GB to 49.7 GB. All 32 potential archive files and the domain
and surface response matrix CSVs were byte-identical. This is a measured setup-cost
improvement, not a claim of a universal full-library speedup.
