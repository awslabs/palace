<!-- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Original-mesh streaming regeneration audit — 2026-09-10

This report concerns the SOCA campaign at
`/data/home/simlap/transmon_library_speed_20260909`, **not** the later graded-mesh
experiments. Meshes, source traces, coupon orders, solver settings and precision
were preserved from the retained library provenance.

## Result

All 20 thin/fabricated numerical cases are complete: **10 real models and 2,138
source solves**. Every worker and reducer returned zero; all expected source
indices were observed and all 2,138 PCG solves converged. The frozen executable
hash and all 1,089 recorded mesh/trace input hashes matched. Original solver-linear
settings and excitation definitions also matched.

The four corner PBS wrappers initially exited with status 1 because their matrix
comparison expected the legacy compact corner format. Their numerical worker and
reducer stages succeeded. The documented core-energy conversion recovered them
without PDE/reducer reruns. The other 16 PBS jobs exited zero.

Independent full symmetric matrix checks verified dimensions, complete basis-index
sets, finite values, metadata identity and energy-matrix sanity. Worst relative
matrix difference from the reference was about **6.38e-7** (0.0000638%). Domain
comparisons were about 3.6e-13 or smaller in the campaign's stored comparison metric.

This establishes faithful regeneration, **not** a new continuum accuracy claim.
The library retains its original mixed orders: five models at p4 and five at p5.
It is not an all-p5 coupon library merely because the historical directory name
contains `p5`.

## Time

Submitted: **2026-09-09 15:18:00 UTC**.
Final manifest originally published: **2026-09-10 00:37:31 UTC**.
Submission-to-publication duration: **9 h 19 min 31 s**.

The table lists worker plus reducer wall time for each job; jobs ran concurrently.
Meshing and prior qualification were reused and are not included as fresh work.

| Model | Coupon order | Sources per thin/fab case | Thin | Fabricated |
|---|---:|---:|---:|---:|
| Isolated edge | 5 | 95 | 00:01:00 | 00:01:02 |
| Same-conductor strip, 2 um | 5 | 96 | 00:01:28 | 00:01:30 |
| Convex 90-degree corner | 4 | 72 | 00:02:25 | 00:02:39 |
| Concave 90-degree corner | 4 | 72 | 00:02:24 | 00:02:43 |
| Two-edge cluster `3f8992613e95` | 4 | 115 | 00:53:19 | 01:28:32 |
| Three-edge cluster `419576fdab24` | 4 | 130 | 01:11:33 | 02:00:24 |
| Four-edge cluster `9d2cb9bbb3fe` | 4 | 80 | 01:17:45 | 02:11:02 |
| Large three-edge cluster `7f03270dca8e` | 5 | 180 | 09:10:35 | 08:37:21 |
| Ten-edge cluster `6791f1c84123` | 5 | 161 | 08:35:43 | 08:00:52 |
| Two-edge cluster `8dd4bc70f183` | 5 | 68 | 00:23:20 | 00:30:59 |

All jobs used r8g.48xlarge nodes. Most used one node; the large three-edge thin/fab
cases used 2/4 nodes, and the ten-edge cases used 4/8 nodes. Peak concurrent node
count was 26. Timed response stages used about **161.5 node-hours / 31,016 allocated
core-hours**. These are allocation-weighted stage times, not measured CPU utilization.

The critical case was large three-edge thin: 8:44:45 worker + 0:25:49 reducer,
with total case overhead bringing it to 9:13:23. The campaign's extra time includes
pilot/control execution, provisioning, comparisons, copying and cleanup. PBS reports
zero `resources_used.walltime` for several large jobs; the application stage clocks
and scheduler start/end timestamps establish the elapsed times instead.

The only direct ordinary-versus-streaming control was the small isolated-edge
thin case: **46.76 s ordinary versus 48.44 + 11.56 = 60.01 s streaming/reduction**.
Streaming was slower for that case. There is no matched whole-campaign ordinary
baseline from which to claim a universal speedup factor. The unchanged large meshes
and iterative solves still dominated this roughly nine-hour campaign.

## Packaging and readiness checks

The generated library is now usable at:

```
/data/home/simlap/transmon_library_speed_20260909/new-library/process-library.json
```

The device smoke test caught a CSV compatibility problem not detected by Python's
numeric comparison: four compact corner matrices used CRLF line endings, so the
native reader saw the wrong final column name. Those four files were normalized
to LF. CSV values were verified unchanged, and originals/hashes are retained in
`audit-20260910/csv-format-backup` and `csv-format-repair.json`. No coupon solves
were repeated. Today's format repair/readiness audit is separate from the numerical
campaign duration above.

The manifest also contains three **geometry-discovery closure descriptors** with
intentionally unmaterialized dummy matrix paths. The original reference has the same
entries and no dummy files. Removing them changes closure and produces missing
requirements, despite a zero preflight process exit code. They must remain as
metadata: the runtime loads matrices only for selected models. Both preflight and
the smoke run confirmed that no discovery descriptor was selected for evaluation.
No dummy numerical response files were invented or added.

Validation with the campaign's frozen executable:

- Fresh, uncached strict preflight: **Complete=true; 3,120 exact matches, zero
  interpolated and zero missing**. Requirements and statistics match the old library.
- Coarse p5 device smoke, 3,071,387 H1 DOFs and 32 ranks: **13 raw iterations and
  four corrected iterations**, matching the reference run.
- Identical nine runtime model instances, identical patch data, no active placeholders.
- Raw domain energy, raw surface participation and capacitance CSVs byte-identical.
- Fixed-trace, fixed-flux and self-consistent participation differences below roughly
  **9.6e-8 relative** (0.0000096%). All corrected observables were finite.
- All 73 files of the original reference-library snapshot still match the original
  retained library. It was not modified.

The coarse runs both retain the known fixed-trace/fixed-flux closure-confidence
warning; the self-consistent solves converged. This validates replacement-library
compatibility and numerical equivalence, not final AMR convergence or a new
high-resolution device result.

## Ready device configuration

```
/data/home/simlap/transmon_library_speed_20260909/single-transmon-generated.json
```

This is the prior p5, AMR-11, full-surface-mortar configuration pointed at the new
library, with `UnmatchedPolicy: Error`, `MortarOversampling: 2`, and fixed-trace
translational domain correction. Raw-only AMR is retained. No full AMR device job
was launched as part of this audit. Any node-local staging/cache must be prepared
for the new library path rather than silently reusing old-library matrix paths.

Detailed machine-readable evidence lives in `audit-20260910/`: `campaign.json`,
`matrix-audit.json`, `pbs-jobs.json`, `device-smoke-comparison.json`,
`original-library-integrity.json`, and `readiness.json`. Failed exploratory removal
of discovery descriptors is retained as an audit artifact and is not a deployment
candidate. The original staging `candidate-library.json` is not the runnable manifest.
