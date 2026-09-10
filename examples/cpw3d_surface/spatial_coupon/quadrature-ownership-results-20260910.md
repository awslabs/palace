<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Quadrature ownership pilot — 2026-09-10

## Why change the partition integration?

Whole-triangle ownership refinement remained expensive: the two-conductor study
reached 123,856 tets while its conservative unresolved-energy bound was still
4.64%. This bound is not the actual error, but it cannot certify the requested
0.5% target. Its final per-slot inter-mesh change was still about 4.37%.

The new opt-in electrostatic postprocessor evaluates the same polygonal ownership
rule at integration points instead of assigning an entire triangle from its
centroid. It preserves all slots and the complete fields/basis supplied to it.
It does not modify geometry, boundary conditions, the PDE solution, or raw AMR.

## Fixed-field results

The quadrature sweeps reused archived p5 conductor-state fields. They did not
repeat PDE solves or change mesh connectivity. Geometry:

- Two conductors: 13,460 tets, two independent conductor states, full 2-by-2 matrices.
- T junction: 15,307 tets, one conductor state.

All ownership slots of each physical family/conductor use the same positive
quadrature rule. The table reports the **largest relative slot-matrix Frobenius
change**, not merely a change in a smooth held-out total.

| Case | Order 12 to 20 | Order 20 to 32 | Order 32 to 48 |
|---|---:|---:|---:|
| Two conductors | 0.05493% | 0.04670% | 0.01295% |
| T junction | 0.09694% | 0.03349% | 0.02659% |

These are empirical convergence observations, not rigorous error bounds for the
continuum field or a complete response library.

Controls:

- Domain response matrices were unchanged at CSV precision.
- Sum of slot matrices agrees with the unpartitioned physical-family matrix to
  at worst about 4.2e-13 relative.
- No negative eigenvalue was found in the slot matrices; the implementation uses
  nonnegative quadrature weights in the Gram construction.
- Independent low-order Dunavant order-20 versus positive Duffy order-48 comparison
  differed by at most 0.06663% for the two-conductor slot matrices.
- Ordinary scalar postprocessing agrees with the response-matrix diagonals to
  about 5.5e-13 relative.
- Changing the characteristic length `Lc` to 7.3 mesh units left every tested
  slot matrix unchanged at CSV precision, checking CSV/mesh nondimensionalization.
- Julia and C++ selectors agree on 12,000 generated queries across the generic
  two-conductor and real ten-edge signatures at two coordinate scales.
- A p5 legacy case with ownership disabled reproduced domain energy, surface
  participation, and capacitance CSVs byte-for-byte.

The first implementation rejected MFEM's order-32 simplex rule because it had
negative weights. It now uses a positive Gauss-product rule with a symmetrized
Duffy map for triangles; polynomial-moment and positive-weight tests cover orders
through 48. Negative weights are not hidden with absolute values.

## What remains a mesh error?

On the coarse two-conductor field, whole-triangle labels differed from the
order-48 quadrature partition by as much as 2.66% in an SA slot matrix, 1.57%
in MS, and 0.43% in MA. This is a measurable bookkeeping/integration discrepancy
on **identical fields**, not a PDE convergence estimate.

Comparing quadrature-partitioned fields on 13,460 and 28,560 tet meshes gives:

- SA/MS slot-matrix changes no larger than about 0.142%.
- MA changes up to about 3.58%.
- On the 28,560 tet mesh, order 20 to 32 changes were at most 0.0367%.

Thus MA still needs physical field/near-metal mesh convergence. The improved
partition integration is not a reason to accept the current coarse mesh.

## Cost and retained evidence

Two-rank local reductions for the coarse two-conductor mesh took approximately
12.0 / 19.6 / 39.3 / 79.1 seconds at orders 12 / 20 / 32 / 48, with peak process-tree
RSS about 1.1–1.3 GiB. The T-junction reductions took 9.3–45.3 seconds. The refined
order-32 reduction exceeded its initial 90-second cap; an independent retry from
the retained archives completed in 114.6 seconds under a 150-second cap, without
re-solving sources.

Evidence is retained under `/tmp/quadrature-ownership-study/`, including per-run
resource JSON, configurations, ownership CSVs, matrix outputs and summary JSON.
The old simplex-rule failure and the timed-out reduction remain as separate
artifacts. No reference library was overwritten.

## Automated verification

- Build and test target completed with no more than six build jobs.
- Config tests: 699 assertions; schema tests: 362 assertions (one pre-existing skip).
- Ownership selector, positive-rule moments, scalar/Gram conservation, unit scaling,
  missing/duplicate slots, and malformed input: 154 assertions, including four-rank MPI.
- Existing surface-response MPI suite: 1,466 assertions passed on two ranks.
- Existing nondimensionalization suite: 96 assertions passed.

## Scope and gates

This is an experimental surface-postprocessing path, currently electrostatic.
The native selector supports straight polygonal ownership segments and optional
conductor-first SA selection. Groups are process-plane-specific. The exporter
fails on unsupported curved plan-view ownership or unsplit multi-plane inputs.
Missing slots, duplicate slots, unequal quadrature orders, absent attributes, and
invalid ownership rows fail closed.

Positive local surface matrices do not prove stability/coercivity of a global
`F-T` correction. Full trace-basis qualification, near-metal field convergence,
curved/multilayer ownership integration, and corrected-device observables remain
required before any production library replacement.
