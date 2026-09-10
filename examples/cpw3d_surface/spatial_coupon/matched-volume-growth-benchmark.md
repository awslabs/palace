<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Matched volume-growth mesh benchmark

## Purpose and scope

The primary optimization target is lower library-generation cost through fewer
FEM unknowns and better iterative-solver convergence. Partition integration is
held fixed here: the pilot uses standard whole-family MA/MS/SA matrices, not a
new ownership-quadrature experiment.

The complete **pilot** basis is one state per conductor plus two affine matching
traces (x and y). These boundary data are exactly representable in both the prism
and tetrahedral p5 trace spaces. Conductors must lie strictly inside the matching
box; the preparation tool rejects contact with the outer boundary rather than
silently using an incompatible zero trace. This pilot is not the complete nodal
trace basis of a production process library, and no modal truncation is installed.

## Mesh experiment

`mesh_graded_tet_experiment.jl --volume-study volume-growth-example.toml` first
creates one surface triangulation, then regenerates only the volume mesh for each
profile. The 3D size field is

```
h(d) = min(hmax, hmin + g_near*min(d,d_switch)
                     + g_far*max(0,d-d_switch))
```

where `d` is the true 3D distance to the physical feature curves. `MinimumSize`
can optionally set `hmin` for the volume independently of the surface size field;
otherwise it uses the surface minimum. This permits a finer interface mesh without
extending that minimum spacing through the volume. Surface sizing is unchanged
between volume profiles. This is not a fixed XY plan swept through
the bulk. The experimental study currently uses exact linear geometry for planar,
zero-rounding coupons; curved geometry requires separate qualification.

Each candidate must preserve an exact physical surface fingerprint. The fingerprint
includes all interface and matching-boundary triangles and their coordinates.
`verify_frozen_volume_meshes.jl` checks it again after serialization. Surface source
basis and interface sampling are therefore fixed within a volume-profile family.
Two independent baseline generations were byte-identical in the tested strip case.

A separate stronger-near-metal candidate changes surface growth/tangential spacing
as well as volume grading. That change is explicit, not confused with the fixed-
surface growth experiment. Its matching excitations remain the same affine/constant
functions.

## Matched prism control and geometry checks

`mesh_prism_control.jl` produces an independent control; reference libraries are
never overwritten. Linear geometry is used for both controls and candidates, with
p5 solution order in Palace.

Two issues were handled before comparing solves:

- Unused prism vertices inside excluded metal cells are no longer registered.
  They were not FEM unknowns and were dropped during serialization, invalidating
  the old linear-geometry node-count check.
- The legacy prism mesh's finite-distance etch collar left small unetched corner
  regions in this closed-strip geometry, unlike the CAD candidate. The control
  explicitly uses `--full-gap-etch`, a new **opt-in** process choice. The legacy
  default and retained library meshes are unchanged. Integrated material volumes
  and interface-family areas then agreed to about 1e-13 relative.

The tested strip uses 100 nm metal, 50 nm overetch, 90-degree walls, zero rounding,
substrate permittivity 11.47, and a 0.5 um coupon radius.

## Execution and measurement

`prepare_matched_mesh_benchmark.py` creates a portable bundle with fixed traces,
configs, mesh/input hashes and audit results. `run_matched_mesh_benchmark.py` runs
it on a dedicated host/allocation, one solve at a time. It fixes:

- one executable hash and host;
- one MPI rank count and binding policy;
- p5, PCG, p-multigrid/BoomerAMG, relative tolerance 1e-8;
- identical source definitions and interface parameters;
- exact streaming workers and reducers, with initial-guess recycling disabled.

`PALACE_RESPONSE_SOURCE_TIMING=1` adds maximum-over-ranks per-source solve and
source-total timings to streaming workers. Source-total includes RHS preparation,
solve, field construction and archiving. The normal elapsed-time report retains
setup/preconditioner breakdowns. Timing on/off produced byte-identical archives in
a three-source, two-rank check.

The benchmark records p5 H1 DOFs, input element/Jacobian-quality statistics,
process-tree RSS, source iterations/times, worker time, reducer time, and their
sum. Mesh-generation stage timings are separate; a production-sized full-library
wall time is not inferred from this small pilot without qualification.

Thin raw edge-inclusive SPR is retained but is not an accuracy gate. Comparisons
include thin/fabricated domain matrices, individual source diagonals, fabricated
interface matrices, and `F-T` (both relative to itself and to the thin domain).

## Bounded failures and recovery

Every command has a time and aggregate process-tree RSS limit. A completed source
is certified by its source-timing line after collective archive completion. If a
worker hits its time cap, `recover_matched_mesh_case.py` can preserve completed
sources and solve only the missing ones, provided the essential-boundary contract
is unchanged. Partial source files are retained separately.

Recovered matrix data are usable for accuracy checks. Interrupted timing is **not**
reported as an uninterrupted total: the original time cap gives a lower bound,
and restart cost is reported separately. A speedup involving that baseline must
be labeled as a lower bound, not an exact ratio.

The 2026-09-10 evidence is under `/tmp/coupon-volume-growth-20260910` locally and
`/data/home/simlap/coupon_volume_growth_matched64_20260910` plus
`/data/home/simlap/coupon_volume_growth_near_20260910` on the dedicated host. A separate
32-rank interrupted run is retained and is not mixed into 64-rank speed ratios.
