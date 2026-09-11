<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Graded-library campaign: preparation and source-contract gate

Campaign root on SOCA:
`/data/home/simlap/transmon_library_graded_20260910`.

The original-mesh campaign remains untouched at
`/data/home/simlap/transmon_library_speed_20260909`.
Its 9 h 19 min 31 s numerical-generation result remains valid as a reproduction
of the retained discrete library, not as a new continuum-accuracy qualification.

## Intended comparison

The six spatial models retain their original source counts and mixed orders:
115/130/80 sources at p4, and 180/161/68 at p5, per thin/fabricated case.
The initial `mesh-plan.json` binds original config, trace and geometry hashes.
The four inexpensive non-spatial model pairs must also be retained or rerun when
a complete ten-model library is assembled; they have not been silently dropped.

All Palace-consumed meshes are **Gmsh MSH 2.2 binary**. A serialized multislotted
mesh passed both MFEM geometry/coverage auditing and an actual two-rank Palace
solve. MSH 4.1 was tested and rejected by this MFEM build; it is not a candidate
output format. Gmsh's 2.2 writer can renumber nodes/elements after physical
relabeling, so certificates are remapped to the serialized IDs and checked by
the downstream diagnostic tagger.

The primary timing comparison uses the same frozen executable as the prior
campaign (`e1d14c903c6b0b41ce65221a57c569ab18c9c804e9b0ac0a6e1dfcffd8cf6cbd`)
and original r8g.48xlarge/rank allocations. Dedicated-host diagnostic timings
are not compared with that campaign. Meshing, relabeling and auditing must be
reported separately: the original numerical campaign reused existing meshes.

## First full-scale geometry checks

Three-edge model `419576fdab24`:

- Thin, physical-feature grading alone: 722,093 tetrahedra, maximum kappa 15.6.
  Its material volumes and physical-family areas match the retained prism mesh
  within approximately 7e-13 relative.
- Thin with explicit side-trace constraints and 50 nm cap sizing: 1,655,754
  tetrahedra, 20,396,628 p4 H1 DOFs, maximum kappa 193.5.
- Fabricated, 0.5 nm surface / 2 nm volume minimum: 2,330,406 surface triangles.
  Surface generation on Gmsh 4.13.1 took 3428.9 seconds. The first volume attempt
  exceeded the deliberately conservative 4-million-element/1.2-million-node cap.
- Resuming the saved surface, without changing resolution or geometry, produced
  5,017,186 tetrahedra and 1,270,847 nodes. Additional Netgen optimization exceeded
  a 600-second guard. A separate HXT-only candidate completed the volume stage
  in 126 seconds, had maximum kappa 82.3, and passed full physical-boundary
  coverage and material/area checks (about 1e-12 relative to the retained mesh).
  This is explicitly a different optimization variant, not a hidden replacement
  of the timed attempt. All frozen surfaces and failed-attempt evidence remain.

No full 20-case response build has been launched on these candidates.

## Why the original-source pilot blocked acceptance

The original six-source pilot used indices 1, 32, 66, 96, 128 and 130. Source 66
was underrepresented by about 99% on the unconstrained boundary and about 90%
with side constraints. **66 and 96 are already in the retained model's
`ZeroTraceIndices`; they are not runtime accuracy gates.** No source was removed
from the planned generation workload.

The active sources still showed a genuine issue: with side constraints and
50 nm caps, source 128 differed by +31.65% and source 130 by +8.83%. That pilot
(job 41600) completed on one r8g node/192 ranks with the frozen reference binary.
It converged in 17–22 iterations per source but does not qualify the mesh.

An independent trace audit found **16 discontinuous runtime-active sources** in
the original three-edge traces: 19–26 and 123–130. For example, source 128 is
one volt at `(-8, -0.5, 2.1)` on a side triangle, but zero at the same point when
interpolated on cap triangle 244. Source 130 likewise has a 0.638 V mismatch at
another collinear cap-edge point. This is a declared trace-function discrepancy,
not solver nonconvergence or a geometry-volume mismatch.

The old `cap_ring` ear clipping could remove the last non-collinear corner and
then silently abandon remaining collinear boundary nodes. Their side hats had
no corresponding cap support. The fix disallows ears whose diagonal skips any
remaining boundary vertex and fails rather than dropping a chain.

A separate corrected trace set:

- retains all 130 original nodes, nodal values and source indices;
- retains all original side triangles;
- increases total matching triangles from 244 to 256;
- passes continuity auditing for every source;
- changes the cap functions and explicitly records `SourceDefinitionChanged`.

Original trace files and the old library are not edited. This is **not** a
byte-identical source migration and must not be described as a pure mesh-only
comparison to the historical nine-hour campaign.

## Corrected-cap controls and remaining gates

Job 41610 completed a same-allocation, same-binary, same-p4/tolerance comparison
of the retained prism and graded tet meshes using the **same corrected cap
traces**, for active indices 1, 19, 26, 32, 123, 128 and 130:

| Case | H1 DOFs | Iterations | Seven-source pilot time |
|---|---:|---|---:|
| Retained prism | 12,808,312 | 47/83/87/135/78/108/90 | 210.21 s |
| Graded tet, side constraints + 50 nm cap | 20,396,628 | 18/15/16/20/13/17/13 | 139.77 s |

That is 1.50x for this pilot only, despite the larger tet DOF count and memory
use. It is **not** a full-library speedup or an accuracy pass. The domain
submatrix differs by 1.541% Frobenius and 15.10% worst energy. Source 130's
individual discrepancy fell from 8.83% to 0.0344% after correcting the cap
functions, but the sharp source 128 remains sensitive to discretization.
Job 41651 repeats these seven active sources at p5 on both meshes to distinguish
reference p4 uncertainty from candidate error. Production orders are not changed
by that diagnostic.

All 12 spatial CAD candidates were constructed. Against the retained mesh,
material/family geometry checks passed for cases 05/06/08 and thin 07. Fabricated
07 exposes a separate mismatch: the prism's finite Euclidean etch collar differs
from the CAD polygon-offset collar. The substrate volume differs by 0.077% and
the ordinary-SA area by 66.46%. It must be reconciled before a mesh-only comparison;
it must not be hidden by aggregating ordinary and recessed interfaces. The large
09 reference's serial quadrature audit hit its 180-second guard; 09/10 reference
checks remain incomplete, not passed by inference.

Acceptance still requires full original-dimensional matrices, fabricated
MA/MS/SA accuracy, domain defect checks, strict transmon preflight, and a device
smoke comparison. These pilots are not a replacement for full-basis generation.
**The candidate library is not yet generated or qualified.**
