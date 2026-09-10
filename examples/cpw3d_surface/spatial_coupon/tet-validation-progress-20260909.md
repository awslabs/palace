<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Graded tetrahedral coupon validation — work in progress

The meshing experiments are not production-qualified. Mesh construction/debugging runs
locally with time/RSS guards. P5 validations run directly on the owner's dedicated
`soca-green-job` c7g.16xlarge host, not new PBS allocations. Existing libraries are unchanged.

## Correct acceptance quantities

Thin zero-thickness-metal edge-inclusive raw SPR is mesh-cutoff dependent and is not an
accuracy gate. Thin coupon validation concerns the domain response operator and boundary
excitation representation. Fabricated coupons additionally require converged interface
energies. `compare_coupon_probe.py` implements this distinction; a probe passing is not
full-operator qualification.

## Findings so far (three-edge coupon, p5)

All meshes use the same physical geometry, material data, original source definitions,
and solve tolerance. The selected source sample is 1, 21, 52, and 180; the independent
held-out excitation is separate.

| Candidate | Thin held-out domain error | Fabricated interface-energy errors, MA/MS/SA |
|---|---:|---:|
| growth 1, horizontal edge cap 100 nm | -0.0312% | -4.570% / -0.266% / -0.0315% |
| growth 0.5, horizontal edge cap 100 nm | -0.0520% | not tested |
| growth 1, horizontal edge cap 20 nm | not tested | +0.831% / -0.0294% / +1.509% |
| full matching-trace edge constraints | +0.00148% | timed out; no valid result |
| matching-trace level constraints | +0.00298% | +0.768% / -0.0166% / +1.515% (20 nm edge cap) |
| matching-trace levels, 10 nm edge cap | not tested | +2.403% / +0.0715% / +2.220% |

These are energy differences against the retained prism reference, not proofs of error
against an exact continuum solution. The thin raw SPRs are deliberately excluded.

The requested 2 nm size-field minimum was NOT an enforced normal layer. Measured median
surface-element altitudes adjacent to physical edges were 50–83 nm in relevant families
with the 100 nm edge cap, and about 17 nm with the 20 nm cap. Good bulk kappa or domain
energy does not establish interface-energy accuracy.

## Boundary representation matters

- Thin metal is now represented as planar surface tools, not full-height mask columns.
  This avoids unnecessary CAD partitions throughout the bulk while keeping the exact mask.
- Matching-trace constraints are derived from the original trace triangulation, not model
  names or edge counts.
- Embedding all matching-trace edges improved the thin held-out domain energy but introduced
  kappa around 10,490 on the far caps and 170 PCG iterations. Seven off-box trace segments
  were skipped/reported; the original trace representation is not everywhere coincident
  with the box surface.
- Constraining only trace height levels retained kappa around 13.6, 8 PCG iterations, and
  0.00298% held-out domain agreement. However, selected thin-domain diagonal matrix entries
  still differed by -0.047%, +0.565%, -4.263%, and +1.187%. The selected submatrix relative
  L2 difference was 1.12%. This candidate is therefore NOT operator-qualified.
- A new candidate retains side-face trace edges but avoids cap diagonals and instead limits
  cap mesh size to 100 nm. It has 453,637 tets and max kappa about 242. It completed
  five sources in 272 s and reduction in 112 s on 32 ranks of the dedicated host.
  Held-out domain error was +0.00255%, but the four original diagonal entries differed
  by +0.163%, -1.196%, -2.369%, and +16.514%; selected-submatrix relative L2 error was
  7.11%. It is NOT operator-qualified.
- Inspection of the prescribed trace mesh itself found cap triangles only about 1 nm
  in altitude, at caps roughly 2 um from the metal plane. Forcing every trace diagonal
  into CAD can create additional poorly shaped elements. Boundary-data discretization
  needs explicit validation; neither a smooth probe nor material-volume agreement is
  sufficient to establish equivalence of the response operators.

## General multi-slot labeling

`label_interface_patches.jl` labels individual boundary triangles, retaining their physical
family/conductor, instead of assigning a whole CAD face from one sample. Boundary element
connectivity/tags and volume meshes remain unchanged. Missing expected interface attributes
fail closed. Multi-plane labeling currently requires explicit layer ownership and is rejected.

The ten-edge scout now contains all six expected SA/thin-metal attribute groups, including
previously missing slot-0 regions. It reports sample disagreement within each triangle and
emits local refinement hints. These are diagnostics, not rigorous bounds on all possible
unresolved regions. Cumulative hint refinement reduced ambiguous fractions in the small
metal regions from about 26% to 6–7% with a modest increase from 84,408 to 92,766 tets.
The partition remains unqualified. Reference slot areas themselves arise from centroid
assignment on another mesh and are not automatically exact continuum boundaries.

## Next gates

- Investigate matching-boundary representation and convergence across both mesh families;
  finer tetrahedral resolution has not yet established agreement with the retained operator.
- Preserve all reference libraries and mesh variants. No tested candidate is qualified yet.
- Continue interface partition refinement or explicit patch-boundary support before claiming
  general multi-slot accuracy.
- Extend the independent excitation set and perform full-matrix/device validation only after
  these pilots pass; no candidate is installed automatically.

## Avoidable archive-worker/reducer memory

The archive-only path still constructed the RT error-estimation hierarchy even with no
flux recovery, and the reducer constructed it despite never recovering a field. This
unused allocation accounted for about 47 GiB in the 10 nm-capped pilot. Construction is
now conditional: ordinary solves retain their estimator, streaming workers retain it
when flux recovery is required, and reducers do not build it. Small two-rank regression
experiments produced byte-identical potential/flux archives (including a fabricated
flux-recovery case) and byte-identical domain/surface matrices. The modified build is
local only at this point; the reported large-pilot timings used the earlier binary.

## Tests

Regression coverage includes exact-mask membership beyond finite edge strips, exclusion of
all conductor footprints from etch collars, holes and rigidly rotated masks, preservation of
boundary connectivity during slot labeling, malformed trace input, native ARM-compatible
callbacks, and the thin-vs-fabricated probe acceptance policy.
