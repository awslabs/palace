<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Interface-slot ownership: certificates, energy bounds, and refinement

For the newer fixed-field quadrature partition pilot, see
[`quadrature-ownership-results-20260910.md`](quadrature-ownership-results-20260910.md).
It avoids requiring a 3D mesh to align with every bookkeeping boundary; its
quadrature convergence and PDE-mesh convergence are deliberately separate checks.

A mesh can have correct physical areas and all expected labels while assigning
significant energy to the wrong slot. Three interior samples per triangle were
not sufficient to detect every crossing. The new workflow separates ownership
uncertainty from field/PDE discretization error.

## Whole-element certificate

`interface_ownership.jl` keeps the same exclusive nearest-feature assignment,
including conductor and process-layer selection. It adds sufficient tests for
constant ownership over the **complete** boundary element:

1. Distance to a fixed segment/arc and minimum distance to a label's feature set
   are 1-Lipschitz. A competing-label margin greater than twice an enclosing
   radius certifies constant ownership.
2. For segment sets, `ownership_bernstein.jl` uses a sharper squared-distance
   test. Linear/quadratic geometric maps yield degree-2/4 polynomials. Bernstein
   coefficients bound those polynomials over the entire triangle. Clamped
   segment regimes are handled by safe point upper bounds for the selected
   segment and infinite-line lower bounds for competitors when necessary.
3. Curved quadratic faces use their Bernstein control hull, not merely their
   corner or mid-edge nodes. Unsupported geometric orders remain unresolved.
4. Layer and, for SA, conductor selection must also be certified. A small slot
   distance alone is insufficient when another conductor/layer can win.

These are sufficient tests: an unresolved triangle may still have one label.
Conversely, interior samples can all agree while a corner region has another
owner; a regression test covers precisely that case. Geometric certification
is not electromagnetic accuracy qualification.

## Disjoint diagnostic mesh copy

`tag_partition_uncertainty.jl` creates a separate mesh with each original interface
split into certified and unresolved subsets. Unresolved tags use an offset of
10,000. Every boundary triangle belongs to exactly one subset. Coordinates,
volume connectivity, and boundary connectivity are unchanged.

The diagnostic solver configuration includes **both subsets** in every original
conductor boundary condition. Interface energies are queried on disjoint subsets
and summed back into canonical slots and physical families. Nothing is disabled,
omitted, or counted twice.

The element-certificate CSV is bound to the serialized mesh by SHA-256, with
signature/process provenance recorded. Missing, stale, tampered, or mismatched
certificates fail before producing the diagnostic copy. A tag-only control solve
checks that fields and reconstructed slot energies are unchanged.

## What the energy bound means

Let `U_g` be the positive interface energy in all unresolved elements of physical
family/conductor group `g`, evaluated with the current FE field. Then each slot's
absolute ownership error is at most `U_g`. The sum of absolute errors across the
group's slots is at most `2 U_g`. Family totals themselves do not change under
redistribution.

The reported group-normalized bound is `U_g / E_g`. Slot-relative bounds use the
computed slot energy as denominator. They are **not** a bound on the PDE field,
the domain defect `F-T`, a different interface model, or final device quantities.
Quadrature accuracy remains a separate obligation; the planar p5 studies use the
solver's whole-interface energy integrals. Small observed inter-mesh changes are
reported separately and are not substituted for a bound.

## Bounded refinement study

`run_partition_refinement.py` runs multi-conductor or single-conductor p5 probes
on generic meshes. It retains physical sizing, process geometry, and solver
tolerance while adding cumulative refinement hints on unresolved regions.

- `TET_SLOT_MINIMUM_SIZE` is independent of physical-edge sizing.
- `graded_size_points.jl` evaluates the exact weighted point-size envelope using
  bounding-tree pruning; it does not replace variable-size sources by an
  unweighted nearest point.
- Every mesher, tagger, worker, and reducer has a time/process-tree RSS limit.
- Tools, inputs, artifacts, and the Palace binary hash are retained.
- `--archive-states` uses prescribed one-volt conductor states with zero matching
  trace, exact streaming archives, and reduction. This avoids unused estimator
  construction without lowering FEM order or tolerance.
- The response-matrix API currently requires localization metadata. The study
  explicitly reads `Q_total_ij`, not the optional localized `Q_ij` term.
- A completed level-zero baseline can be reused only with matching geometry and
  solve mode. Failed or incomplete stages are not promoted to successful data.

Example (each case is a directory prepared by the general-mesh suite):

```sh
python3 examples/cpw3d_surface/spatial_coupon/run_partition_refinement.py \
  --inputs /tmp/general-study --root /tmp/new-partition-study \
  --palace /path/to/palace --julia-project /path/to/julia/environment \
  --case two-conductors --case t-junction --iterations 3 --ranks 2 \
  --archive-states
```

The default study target is 0.5% group-normalized ownership uncertainty, with
separate 0.05% domain and 0.5% group/slot inter-mesh checks. These are provisional
engineering diagnostics. Reports always retain `LibraryQualified: false`.
